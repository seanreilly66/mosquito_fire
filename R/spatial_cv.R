suppressPackageStartupMessages({
  library(sf)
  library(bit)
  library(rsample)
  library(tibble)
  library(terra)
  library(predicts)
  library(concaveman)
})

# ------------------------------------------------------------------------------
# Spatial block membership
# ------------------------------------------------------------------------------
# Divide ROI into n blocks using predicts::divider (seeded),
# buffer each block by buffer_dist,
# assign pts to blocks and to buffered blocks,
# return incl_bits and excl_bits (bitsets per block).

spatial_block_membership <- function(roi = NULL,
                                     pts,
                                     n_blocks,
                                     buffer_dist = 2000,
                                     seed = 1) {
  # --- Block prep ---
  
  if (is.null(roi)) {
    roi <- derive_roi_from_pts(pts)
  }
  
  if (!inherits(roi, "SpatVector")) {
    roi = terra::vect(roi)
  }
  
  set.seed(seed)
  blocks = predicts::divider(roi, n = n_blocks) |>
    st_as_sf()
  
  blocks_buf = blocks |>
    st_buffer(dist = buffer_dist)
  
  n_pts <- nrow(pts)
  
  # --- Inclusion: block contains point ---
  
  incl_hits <- sf::st_contains(blocks, pts, sparse = TRUE)
  
  incl_bits <- lapply(incl_hits, function(i) {
    b <- bit::bit(n_pts)
    if (length(i)) b[i] <- TRUE
    b
  })
  
  # --- Exclusion: buffered block contains point ---
  
  excl_hits <- sf::st_contains(blocks_buf, pts, sparse = TRUE)
  
  excl_bits <- lapply(excl_hits, function(i) {
    b <- bit::bit(n_pts)
    if (length(i)) b[i] <- TRUE
    b
  })
  
  list(
    incl_bits = incl_bits, 
    excl_bits = excl_bits,
    blocks = blocks,
    n_blocks = length(incl_bits),
    n_pts = n_pts)
  
}


# ------------------------------------------------------------------------------
# Spatial block folds
# ------------------------------------------------------------------------------

# Generate repeated CV fold assignments at block-level
# Return bit mapping of pts to blocks + buffer
# Return integer vector of block to fold assignments
# n_blocks = v * blocks_per_fold.

spatial_block_folds <- function(roi = NULL,
                                pts,
                                v,
                                blocks_per_fold,
                                repeats,
                                buffer_dist,
                                seed) {
  
  if (is.null(roi)) {
    roi <- derive_roi_from_pts(pts)
  }
  
  theo_n_blocks <- v * blocks_per_fold
  
  block_mem <- spatial_block_membership(
    roi = roi,
    pts = pts,
    n_blocks = theo_n_blocks,
    buffer_dist = buffer_dist,
    seed = seed
  )
  
  n_blocks <- block_mem$n_blocks
  
  set.seed(seed)
  fold_maps <- vector('list', repeats)
  
  for (r in seq_len(repeats)) {
    perm <- sample(seq_len(n_blocks))
    folds <- rep(seq_len(v), each = blocks_per_fold)
    fold_map <- integer(n_blocks)
    fold_map[perm] <- folds
    fold_maps[[r]] <- fold_map
  }
  
  list(
    incl_bits = block_mem$incl_bits,
    excl_bits = block_mem$excl_bits,
    fold_maps = fold_maps,
    buffer = buffer_dist,
    blocks_per_fold = blocks_per_fold
  )
  
}

# ------------------------------------------------------------------------------
# spatial_block_rset()
# ------------------------------------------------------------------------------

# Convert incl_bits / excl_bits / fold_maps into an rset compatible with tidymodels
# and spatialsample::autoplot() (spatial_rset class).
#
# Inputs:
#   fold_plan: output of spatial_block_folds() with fields:
#              incl_bits, excl_bits, fold_maps
#   pts:    the sf point data used to build bit indices (row order must match)
#   v:         number of folds
#   subsample_repeats: optional integer; if provided, sample that many repeats
#   subsample_seed: seed for repeat subsampling only
#
# Notes:
# - Fold labels and repeat labels are stable and zero-padded for readability.
# - Repeats are sampled *without replacement* and returned in ascending order.

spatial_block_rset <- function(fold_plan, pts, v, subsample_repeats = NULL, subsample_seed = 1) {
  
  incl_bits  <- fold_plan$incl_bits
  excl_bits  <- fold_plan$excl_bits
  fold_maps  <- fold_plan$fold_maps
  
  repeats <- length(fold_maps)
  n_blocks <- length(incl_bits)
  
  # Select repeats
  reps_to_use <- seq_len(repeats)
  if (!is.null(subsample_repeats) && subsample_repeats < repeats) {
    set.seed(subsample_seed)
    reps_to_use <- sort(sample(reps_to_use, subsample_repeats))
  }
  
  n_reps_used <- length(reps_to_use)
  
  n_splits <- n_reps_used * v
  splits <- vector("list", n_splits)
  id  <- character(n_splits)
  id2 <- character(n_splits)
  
  idx <- 0L
  for (i in seq_along(reps_to_use)) {
    r <- reps_to_use[i]
    fold_map <- fold_maps[[r]]
    rep_lab <- sprintf("Repeat%02d", r)
    
    for (k in seq_len(v)) {
      idx <- idx + 1L
      
      assess_blocks <- which(fold_map == k)

      # Compute masks
      incl_mask <- Reduce(`|`, incl_bits[assess_blocks])
      excl_mask <- Reduce(`|`, excl_bits[assess_blocks])
      
      # Convert to indices
      assess_idx <- as.integer(bit::as.which(incl_mask))
      train_idx  <- as.integer(bit::as.which(!excl_mask))
      
      # Create split
      splits[[idx]] <- rsample::make_splits(
        x = list(analysis = train_idx, assessment = assess_idx),
        data = pts,
        class = 'spatial_rsplit'
      )
      
      id[[idx]]  <- sprintf("Fold%02d", k)
      id2[[idx]] <- rep_lab
    }
  }
  
  cv_att <- list(
    v = v,
    repeats = n_reps_used,
    n_blocks = n_blocks,
    blocks_per_fold = fold_plan$blocks_per_fold,
    buffer = fold_plan$buffer,
    method = 'random'
  )
  
  rsample::new_rset(
    splits = splits,
    ids = tibble::tibble(id2 = id2, id = id),
    attrib = cv_att,
    subclass = c("spatial_block_cv", "spatial_rset", "rset")
  )
}

# ------------------------------------------------------------------------------
# spatial_init_split
# ------------------------------------------------------------------------------
# Generate spatial initial split
# Generates n + 1 blocks needed to capture proportion (1 / 1 - prop)
# returns an rsplit like rsample::initial_split.
# If anchor_edge = TRUE, selects from subset of blocks with edge proportion 
# greater than edge_prop


spatial_init_split <- function(roi = NULL,
                               pts,
                               buffer_dist,
                               test_prop,
                               seed = 1,
                               anchor_edge = TRUE,
                               edge_prop = 0.3) {
  
  stopifnot(test_prop > 0, test_prop < 1)
  
  if (is.null(roi)) {
    roi <- derive_roi_from_pts(pts)
  }
  
  theo_n_blocks <- ceiling(1 / test_prop)
  
  bits <- spatial_block_membership(
    roi = roi,
    pts = pts,
    n_blocks = theo_n_blocks,
    buffer_dist = buffer_dist,
    seed = seed
  )
  
  blocks <- bits$blocks
  n_blocks <- bits$n_blocks
  
  if (anchor_edge) {
    
    roi_bnd <- roi |>
      st_as_sf() |>
      st_geometry() |>
      st_boundary()
    
    blk_bnd <- blocks |>
      st_geometry() |>
      st_boundary()
    
    blk_len <- blk_bnd |>
      st_length() |>
      as.numeric()
      
    touch_len <- st_intersection(blk_bnd, roi_bnd) |>
      st_length() |>
      as.numeric()
    
    edge_ids <- which(touch_len/blk_len >= edge_prop)

    stopifnot(length(edge_ids) > 0)
    
    cand_ids <- edge_ids
  
  } else {
    
    cand_ids <- seq_len(n_blocks)
    
  }
  
  set.seed(seed)
  rand_block <- sample(cand_ids, 1)
  
  assess_mask <- bits$incl_bits[[rand_block]]
  excl_mask <- bits$excl_bits[[rand_block]]
  
  assess_idx <- as.integer(which(as.logical(assess_mask)))
  train_idx <- as.integer(which(!as.logical(excl_mask)))
  
  init_split <- rsample::make_splits(
    x = list(analysis = train_idx, assessment = assess_idx),
    data = pts,
    class = 'spatial_rsplit'
  )

  class(init_split) <- c('initial_split', class(init_split))
  init_split
  
}

# ==============================================================================

# ==============================================================================
# Helper functions
# ==============================================================================

derive_roi_from_pts <- function(pts) {
  roi_sf <- concaveman::concaveman(pts)
  roi_sf <- st_make_valid(roi_sf)
  roi_sf <- st_union(roi_sf)
  terra::vect(roi_sf)
}