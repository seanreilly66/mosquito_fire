# ==============================================================================
#
# Raster grid metric extraction
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 12 Oct 2024
# Last commit: 12 Dec 2024
#
# Status: Needs documentation
#
# ==============================================================================
#
# Description:
#
# ==============================================================================
#
# User inputs:
#
# ==============================================================================
#
# Package dependencies:
#
# ==============================================================================
#
# Known problems:
#
# ==============================================================================

# library(sf)
library(tidyverse)
library(glue)
library(lidR)
library(terra)
library(yaml)
library(here)
library(fs)
library(future)

# ==============================================================================
# ================================= User inputs ================================
# ==============================================================================

config <- read_yaml(here('config.yaml'))

data_root <- config$data_root_ms

las_folder <- path(data_root, config$norm_lidar_folder)

burn_file <- path(data_root, config$burn_file)

# ==============================================================================
# ====================== Extract metrics per merged group ======================
# ==============================================================================

source(config$las_metric_function)

fire_rast <- rast(burn_file)

las_files <- list.files(las_folder, pattern = '.las$', full.names = T)

ctg <- readLAScatalog(folder = las_folder,filter = "-drop_withheld")

# lidR:::catalog_laxindex(ctg)

opt_filter(ctg) <- '-drop_withheld -drop_class 7 18'
opt_select(ctg) <- 'xyzr'

future::plan(multisession)

metrics <- pixel_metrics(ctg,
                         func = ~ las_cld_metrics(z = Z, r = ReturnNumber),
                         res = fire_rast)
# 
# 
# for (las_i in las_files) {
#   
#   las <- readLAS(las_i) %>%
#     filter_duplicates() %>%
#     filter_poi(Classification != 7, Classification != 18) %>% # Remove errors
#     normalize_height(tin())
#   
#   metrics <- pixel_metrics(las,
#                   # func = ~ las_cld_metrics(z = Z, r = ReturnNumber),
#                   func = ~ las_ladder_evenness(z = Z, r = ReturnNumber),
#                   res = fire_rast)
#   
#   for (i in names(metrics)) {
#     
#     layer <- metrics[[i]]
#     output <- str_replace(las_i, '\\.las', glue('_{i}.tif')) %>%
#       str_replace('las/als_2022', 'rast_metrics')
#     
#     writeRaster(layer, output, overwrite = T)
#     
#     message('Exporting: ', output)
#     
#   }
# 
# }
