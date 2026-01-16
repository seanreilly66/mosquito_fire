suppressPackageStartupMessages({
  library(sf)
  library(bit)
  library(rsample)
  library(tibble)
})

# =============================================================================
# Minimal exact spatial resampling with point-level buffering (predicts + sf + bit + rsample)
# Four functions: fun_1 .. fun_4
# Assumptions: ROI/points are valid sf objects in a projected CRS (meters),
#              points all lie within ROI, and divider partitions ROI cleanly.
# =============================================================================

# ------------------------------- fun_1 ---------------------------------------
# Divide ROI into n blocks using predicts::divider (seeded),
# buffer each block by buffer_dist,
# assign points to blocks and to buffered blocks,
# return incl_bits and excl_bits (bitsets per block).

fun_1 <- function(roi, points, n_blocks, buffer_dist = 2000, seed = 1) {
  
  # Block and point prep
  
  if (!inherits(roi, "SpatVector")) {
    roi = terra::vect(roi)
  }
  
  set.seed(seed)
  blocks = predicts::divider(roi, n = n_blocks) |>
    st_as_sf() 
  
  blocks$block_id <- seq_len(nrow(blocks))
  
  blocks_buf = blocks %>%
    st_buffer(dist = buffer)
  
  points$.pt_id <- seq_len(nrow(points))
  n_points <- nrow(points)
  
  # Inclusion: point to block assignment
  
  tictoc::tic()
  incl_list <- st_contains(blocks, points, sparse = TRUE) |>
    lapply(function(i) points$.pt_id[i])
  
  names(incl_list) <- as.character(blocks$block_id)
  tictoc::toc()
  
  # Exclusion: point is within buffer(block)

  tictoc::tic()
  excl_list <- st_contains(blocks_buf, points, sparse = TRUE) |>
    lapply(function(i) points$.pt_id[i])
  
  names(excl_list) <- as.character(blocks_buf$block_id)
  tictoc::toc()

  
  
  # Convert lists to packed bitsets (1 bit per point per block)
  incl_bits <- vector("list", n_blocks); names(incl_bits) <- as.character(blocks$block_id)
  excl_bits <- vector("list", n_blocks); names(excl_bits) <- as.character(blocks$block_id)
  
  for (b in names(incl_bits)) {
    bv <- bit(n_points)
    idx <- incl_list[[b]]
    if (length(idx)) bv[idx] <- TRUE
    incl_bits[[b]] <- bv
  }
  
  for (b in names(excl_bits)) {
    bv <- bit(n_points)
    idx <- excl_list[[b]]
    if (length(idx)) bv[idx] <- TRUE
    excl_bits[[b]] <- bv
  }
  
  points$.pt_id <- NULL
  
  list(
    points = points,
    blocks = blocks,
    blocks_buffer = blocks_buf,
    n_points = n_points,
    n_blocks = n_blocks,
    buffer_dist = buffer_dist,
    seed = seed,
    incl_bits = incl_bits,
    excl_bits = excl_bits
  )
}

# ------------------------------- fun_2 ---------------------------------------
# initial_split-like function:
# calls fun_1, selects ~ (1-prop) points worth of blocks as assessment,
# returns an rsplit like rsample::initial_split.
fun_2 <- function(roi, points, n_blocks, buffer_dist = 2000, prop = 0.75, seed = 1) {
  bits <- fun_1(roi, points, n_blocks, buffer_dist, seed)
  
  set.seed(seed)
  target <- round((1 - prop) * bits$n_points)
  
  perm <- sample(seq_len(bits$n_blocks))
  test_blocks <- integer()
  test_mask <- bit(bits$n_points)
  
  for (b in perm) {
    cand <- test_mask | bits$incl_bits[[as.character(b)]]
    if (sum(cand) <= target || length(test_blocks) == 0) {
      test_mask <- cand
      test_blocks <- c(test_blocks, b)
    }
    if (sum(test_mask) >= target) break
  }
  
  assess_idx <- which(as.logical(test_mask))
  
  excl_mask <- Reduce(`|`, bits$excl_bits[as.character(test_blocks)])
  train_idx <- which(!as.logical(excl_mask))
  
  sp <- rsample::make_splits(
    x = list(analysis = train_idx, assessment = assess_idx),
    data = bits$points,
    class = c("spatial_rsplit")
  )
  class(sp) <- c("initial_split", class(sp))
  sp
}

# ------------------------------- fun_3 ---------------------------------------
# Generate repeated CV fold assignments at block-level using fun_1.
# n_blocks = v * blocks_per_fold.
# Returns: bits + fold_maps[[rep]] = integer vector (block_id -> fold 1..v).
fun_3 <- function(roi, points, v = 10, blocks_per_fold, repeats = 1,
                  buffer_dist = 2000, seed = 1) {
  n_blocks <- v * blocks_per_fold
  bits <- fun_1(roi, points, n_blocks, buffer_dist, seed)
  
  set.seed(seed)
  fold_maps <- vector("list", repeats)
  
  for (r in seq_len(repeats)) {
    perm <- sample(seq_len(n_blocks))
    folds <- rep(seq_len(v), each = blocks_per_fold)
    fold_map <- integer(n_blocks)
    fold_map[perm] <- folds
    fold_maps[[r]] <- fold_map
  }
  
  list(
    bits = bits,
    v = v,
    blocks_per_fold = blocks_per_fold,
    repeats = repeats,
    seed = seed,
    fold_maps = fold_maps
  )
}

# ------------------------------- fun_4 ---------------------------------------
# Convert fold assignments from fun_3 into an rset compatible with tidymodels,
# and with spatialsample::autoplot() (spatial_rset class).
# Optionally subsample repeats (for racing).
fun_4 <- function(fold_plan, subsample_repeats = NULL, subsample_seed = 1) {
  bits <- fold_plan$bits
  v <- fold_plan$v
  repeats <- fold_plan$repeats
  
  reps_to_use <- seq_len(repeats)
  if (!is.null(subsample_repeats) && subsample_repeats < repeats) {
    set.seed(subsample_seed)
    reps_to_use <- sort(sample(reps_to_use, subsample_repeats))
  }
  
  splits <- list()
  id <- character()
  id2 <- character()
  
  for (rr in seq_along(reps_to_use)) {
    r <- reps_to_use[[rr]]
    fold_map <- fold_plan$fold_maps[[r]]
    
    for (k in seq_len(v)) {
      assess_blocks <- which(fold_map == k)
      
      assess_mask <- Reduce(`|`, bits$incl_bits[as.character(assess_blocks)])
      excl_mask   <- Reduce(`|`, bits$excl_bits[as.character(assess_blocks)])
      
      assess_idx <- which(as.logical(assess_mask))
      train_idx  <- which(!as.logical(excl_mask))
      
      sp <- rsample::make_splits(
        x = list(analysis = train_idx, assessment = assess_idx),
        data = bits$points,
        class = c("spatial_rsplit")
      )
      
      splits[[length(splits) + 1L]] <- sp
      id[[length(id) + 1L]] <- paste0("Fold", k)
      id2[[length(id2) + 1L]] <- paste0("Repeat", rr)
    }
  }
  
  rsample::new_rset(
    splits = splits,
    ids = tibble(id2 = id2, id = id),
    attrib = list(v = v, repeats = length(reps_to_use), buffer = bits$buffer_dist),
    subclass = c("spatial_rset", "rset")
  )
}

# =============================================================================
# Example
# plan <- fun_3(roi, pts, v=10, blocks_per_fold=20, repeats=100, buffer_dist=2000, seed=1)
# rset_small <- fun_4(plan, subsample_repeats=15, subsample_seed=99)
# rset_full  <- fun_4(plan)
# finetune::tune_race_anova(wf, resamples=rset_small, grid=grid, ...)
# spatialsample::autoplot(rset_small)
# =============================================================================
