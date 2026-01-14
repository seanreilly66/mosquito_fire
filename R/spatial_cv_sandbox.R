suppressPackageStartupMessages({
  library(tidyverse)
  library(terra)
  # library(tidymodels)
  # library(spatialsample)
  library(sf)

})

# ==============================================================================
# Testing data
# ==============================================================================

burn_file  <- 'data/burn_severity/mosquito_rdnbr_mtbs.tif'     # RdNBR response
metric_dir <- "data/fuelrasters"            # predictor rasters

roi_gpkg = 'data/mosquito_fire.gpkg'
roi_layer = 'mosquito_study_area'

initial_prop = 0.2
v_folds = 10
blocks_per_fold = 5
tune_reps = 30
fit_reps = 100
buffer_dist = 2000

# ------------------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------------------

roi <- vect(roi_gpkg, roi_layer)

burn_rast <- rast(burn_file) %>%
  crop(roi) %>%
  mask(roi)
names(burn_rast) <- 'rdnbr'

pred_rast <- rast(list.files(metric_dir, pattern="\\.tif$", full.names=TRUE)) %>%
  crop(roi) %>%
  mask(roi)

input_df <- c(burn_rast, pred_rast) %>%
  crop(roi) %>%
  mask(roi) %>%
  as.points() %>%
  st_as_sf() %>%
  mutate(st_coordinates(geometry) |> as_tibble()) %>%
  slice_sample(n = 10000)

v = v_folds
bpf = blocks_per_fold
buffer = 2000
seeds = c(7127)

# ==============================================================================

cellsize = 2*(sqrt(2*st_area(roi)/(3*sqrt(3)*v*bpf)))


spatial_cv <- function(sf_df, roi, v, bpf, buffer, seeds) {
  
  if (!inherits(roi, "SpatVector")) {
    roi = vect(roi)
  }
  
  set.seed(seeds[1])
  blocks = predicts::divider(roi, n = v * bpf) %>%
    st_as_sf() %>%
    mutate(block_id = row_number())
  
  buffer_blocks = blocks %>%
    st_buffer(dist = buffer)
  
  pts_id
  
  
  
  
  
  
}




# =============================================================================
# Exact spatial CV with point-level 2 km exclusion buffer (sf + terra)
# Precompute once; generate ~1300 splits with set algebra (no per-split geometry)
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(data.table)
})

# -----------------------------------------------------------------------------
# 0) Helpers
# -----------------------------------------------------------------------------

assert_projected_meters <- function(x) {
  cr <- sf::st_crs(x)
  if (is.na(cr$wkt)) stop("CRS missing.")
  if (sf::st_is_longlat(x)) {
    stop("Geometry is lon/lat. Reproject to a meter-based projected CRS (UTM/Albers).")
  }
  invisible(TRUE)
}

as_sorted_unique_int <- function(x) sort(unique(as.integer(x)))

union_int_lists <- function(lst) {
  # Fast-enough baseline: concatenate -> sort -> unique
  # For bigger workloads, see bitset option further down.
  if (length(lst) == 0) return(integer())
  as_sorted_unique_int(unlist(lst, use.names = FALSE))
}

int_setdiff <- function(universe, remove) {
  # universe must be sorted unique
  # remove should be sorted unique
  if (length(remove) == 0) return(universe)
  universe[!fmatch::fmatch(universe, remove, nomatch = 0L) > 0L]
}

# -----------------------------------------------------------------------------
# 1) Build blocks (grid clipped to ROI)
# -----------------------------------------------------------------------------
build_blocks <- function(roi_sf, block_size_m, block_id_col = "block_id") {
  assert_projected_meters(roi_sf)
  
  # Create a square grid over ROI bbox; clip to ROI
  g <- sf::st_make_grid(roi_sf, cellsize = block_size_m, square = TRUE)
  blocks <- sf::st_sf(geometry = g)
  blocks <- sf::st_intersection(blocks, sf::st_union(roi_sf))
  blocks <- blocks[!sf::st_is_empty(blocks), , drop = FALSE]
  blocks[[block_id_col]] <- seq_len(nrow(blocks))
  
  # Ensure valid geometries
  blocks <- sf::st_make_valid(blocks)
  
  blocks
}

# -----------------------------------------------------------------------------
# 2) Generate point centers from a terra SpatRaster (30 m) restricted to ROI
#    (If you already have points, skip this and provide an sf POINT object)
# -----------------------------------------------------------------------------
raster_to_points_sf <- function(r, roi_sf = NULL, keep_cell_id = TRUE) {
  # r: terra SpatRaster (any layers; we only need geometry)
  # Returns sf POINT with pt_id and optional cell index
  if (!inherits(r, "SpatRaster")) stop("r must be a terra SpatRaster.")
  
  # Convert raster cells to points (centers); keep only non-NA in first layer
  # For a pure index raster, set all cells non-NA first.
  pts <- terra::as.points(r, values = FALSE, na.rm = TRUE)
  
  if (!is.null(roi_sf)) {
    assert_projected_meters(roi_sf)
    roi_v <- terra::vect(roi_sf)
    pts <- terra::crop(pts, roi_v)
    pts <- terra::mask(pts, roi_v)
  }
  
  # Build sf points
  sf_pts <- sf::st_as_sf(pts)
  
  if (keep_cell_id) {
    # terra points include cell index in "cell" only if values kept;
    # If absent, create pt_id only.
    sf_pts$pt_id <- seq_len(nrow(sf_pts))
  } else {
    sf_pts$pt_id <- seq_len(nrow(sf_pts))
  }
  
  sf_pts
}

# -----------------------------------------------------------------------------
# 3) One-time: point -> block membership (assessment membership)
# -----------------------------------------------------------------------------
precompute_pts_in_block <- function(points_sf, blocks_sf,
                                    pt_id_col = "pt_id",
                                    block_id_col = "block_id") {
  assert_projected_meters(points_sf)
  assert_projected_meters(blocks_sf)
  
  # Spatial join: each point within exactly one block (after clipping)
  joined <- sf::st_join(
    points_sf[, pt_id_col, drop = FALSE],
    blocks_sf[, block_id_col, drop = FALSE],
    join = sf::st_within,
    left = FALSE
  )
  
  if (nrow(joined) == 0) stop("No point-block matches. Check CRS / extents / ROI clipping.")
  
  dt <- as.data.table(sf::st_drop_geometry(joined))
  setnames(dt, c(pt_id_col, block_id_col), c("pt_id", "block_id"))
  
  # pts_in_block: list(block_id -> integer pt_ids)
  pts_in_block <- dt[, .(pt_id = as_sorted_unique_int(pt_id)), by = block_id][
    order(block_id)
  ]
  
  # Return list + a point->block vector if needed
  pt_block <- dt[order(pt_id), block_id]
  names(pt_block) <- dt[order(pt_id), pt_id]
  
  list(
    pts_in_block = setNames(pts_in_block$pt_id, pts_in_block$block_id),
    pt_block = pt_block
  )
}

# -----------------------------------------------------------------------------
# 4) One-time: point -> buffered-block membership (exclusion engine)
#     pts_in_buf[[b]] = all points within buffer(block_b, 2000 m)
# -----------------------------------------------------------------------------
precompute_pts_in_buf <- function(points_sf, blocks_sf,
                                  buffer_m = 2000,
                                  pt_id_col = "pt_id",
                                  block_id_col = "block_id") {
  assert_projected_meters(points_sf)
  assert_projected_meters(blocks_sf)
  
  blocks_buf <- sf::st_buffer(sf::st_make_valid(blocks_sf), dist = buffer_m)
  
  # st_intersects returns sparse list: for each point, which buffered blocks contain it
  hits <- sf::st_intersects(points_sf, blocks_buf, sparse = TRUE)
  
  pt_ids <- points_sf[[pt_id_col]]
  if (is.null(pt_ids)) stop("pt_id_col not found on points_sf.")
  
  # Build a long table: (pt_id, block_id_hit)
  # This inverts the point->blocks list into block->points list.
  lens <- lengths(hits)
  if (sum(lens) == 0) stop("No point-buffer hits. Check buffer/CRS/geometry.")
  
  dt <- data.table(
    pt_id = rep.int(pt_ids, lens),
    block_id = as.integer(unlist(hits, use.names = FALSE))
  )
  
  # Map numeric index back to actual block_id values
  # Because blocks_buf is in the same row order as blocks_sf
  block_ids <- blocks_sf[[block_id_col]]
  dt[, block_id := block_ids[block_id]]
  
  pts_in_buf <- dt[, .(pt_id = as_sorted_unique_int(pt_id)), by = block_id][
    order(block_id)
  ]
  
  setNames(pts_in_buf$pt_id, pts_in_buf$block_id)
}

# -----------------------------------------------------------------------------
# 5) Repeated fold assignment at block level
# -----------------------------------------------------------------------------
make_block_folds <- function(block_ids, v = 10L, reps = 1L, seed = 1L) {
  set.seed(seed)
  block_ids <- as.integer(block_ids)
  n <- length(block_ids)
  
  # For each rep, permute blocks and assign folds as evenly as possible
  folds_by_rep <- vector("list", reps)
  
  for (r in seq_len(reps)) {
    perm <- sample(block_ids, size = n, replace = FALSE)
    fold <- rep.int(seq_len(v), length.out = n)
    # map: block_id -> fold
    m <- integer(max(block_ids))
    m[perm] <- fold
    folds_by_rep[[r]] <- m
  }
  
  folds_by_rep
}

# -----------------------------------------------------------------------------
# 6) Split generation (exact) using precomputed lists
# -----------------------------------------------------------------------------
make_split <- function(folds_map, k,
                       pts_in_block, pts_in_buf,
                       all_pt_ids = NULL) {
  # folds_map: integer vector indexed by block_id giving fold id
  # k: fold number
  # pts_in_block: list(block_id -> integer pt_ids)
  # pts_in_buf:   list(block_id -> integer pt_ids within 2 km buffer)
  # all_pt_ids: optional sorted unique integer vector of all point ids
  
  # Assessment blocks for fold k
  block_ids <- which(folds_map == k)
  if (length(block_ids) == 0) stop("No blocks assigned to fold ", k)
  
  # Assessment points: union of points in those blocks
  assess_pts <- union_int_lists(pts_in_block[as.character(block_ids)])
  
  # Excluded points: union of points in buffers of those blocks
  excluded_pts <- union_int_lists(pts_in_buf[as.character(block_ids)])
  
  if (is.null(all_pt_ids)) {
    # derive from pts_in_block if not supplied
    all_pt_ids <- as_sorted_unique_int(unlist(pts_in_block, use.names = FALSE))
  }
  
  # Training points are all points minus excluded points
  train_pts <- setdiff(all_pt_ids, excluded_pts)
  
  list(
    assess_blocks = block_ids,
    assess_pts = assess_pts,
    excluded_pts = excluded_pts,
    train_pts = train_pts
  )
}

# -----------------------------------------------------------------------------
# 7) Optional: build a compact split manifest for downstream frameworks
# -----------------------------------------------------------------------------
split_manifest <- function(rep_id, fold_id, split_obj) {
  data.table(
    rep = rep_id,
    fold = fold_id,
    pt_id = c(split_obj$train_pts, split_obj$assess_pts),
    role = c(rep.int("train", length(split_obj$train_pts)),
             rep.int("assess", length(split_obj$assess_pts)))
  )
}

# =============================================================================
# Usage pattern (end-to-end)
# =============================================================================
# Inputs you provide:
#   roi_sf: sf POLYGON/MULTIPOLYGON (projected, meters)
#   r:      terra SpatRaster (30 m) or points_sf directly
#
# Example skeleton:
#
# roi_sf <- st_read("roi.gpkg") |> st_transform(26910)   # example
# r <- rast("your_30m_raster.tif") |> project("EPSG:26910")
#
# blocks_sf <- build_blocks(roi_sf, block_size_m = 2000)  # block size is your choice
# points_sf <- raster_to_points_sf(r, roi_sf)
#
# pre_block <- precompute_pts_in_block(points_sf, blocks_sf)
# pts_in_block <- pre_block$pts_in_block
#
# pts_in_buf <- precompute_pts_in_buf(points_sf, blocks_sf, buffer_m = 2000)
#
# all_pt_ids <- as_sorted_unique_int(points_sf$pt_id)
# block_ids  <- blocks_sf$block_id
#
# folds_by_rep <- make_block_folds(block_ids, v = 10, reps = 130, seed = 123)
#
# # Generate all splits without further geometry:
# manifests <- vector("list", length(folds_by_rep) * 10L)
# idx <- 1L
# for (r in seq_along(folds_by_rep)) {
#   fm <- folds_by_rep[[r]]
#   for (k in 1:10) {
#     sp <- make_split(fm, k, pts_in_block, pts_in_buf, all_pt_ids)
#     manifests[[idx]] <- split_manifest(r, k, sp)
#     idx <- idx + 1L
#   }
# }
# manifest_dt <- rbindlist(manifests, use.names = TRUE)
#
# manifest_dt now gives exact point-level train/assess membership with a strict 2 km rule.

# =============================================================================
# Scaling note: bitset acceleration (optional)
# =============================================================================
# If unions become your bottleneck at ~1300 splits, switch pts_in_buf (and/or
# pts_in_block) to bitsets:
#
#   library(bit)
#   n <- length(all_pt_ids)
#   buf_bitsets <- lapply(pts_in_buf, function(ids) { b <- bit(n); b[ids] <- TRUE; b })
#   # for fold:
#   excluded_b <- Reduce(`|`, buf_bitsets[as.character(block_ids_assess)])
#   train_ids  <- which(!excluded_b)
#
# This is typically very fast; memory cost is ~n/8 bytes per block bitset.




