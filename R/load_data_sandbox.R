suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
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
  mutate(st_coordinates(geometry) |> as_tibble())

# %>%
#   slice_sample(n = 10000)

v = v_folds
bpf = blocks_per_fold
buffer = 2000
seeds = c(7127)

points = input_df
seed = 1
buffer_dist = 2000
n_blocks = 30
# ==============================================================================