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

library(sf)
library(tidyverse)
library(glue)
library(lidR)
library(doParallel)
library(terra)

# ==============================================================================
# ================================= User inputs ================================
# ==============================================================================

# File parameters

las_folder <- 'data/mosquito/las/als_2022'

raster_file <- 'data/mosquito/mosquito_reproj.tif'

# ==============================================================================
# ====================== Extract metrics per merged group ======================
# ==============================================================================

# source('R/las_metric_function.R')
source('R/las_metric_ladder_evenness.R')

fire_rast <- rast(raster_file)

las_files <- list.files(las_folder, full.names = T)

for (las_i in las_files) {
  
  las <- readLAS(las_i) %>%
    filter_duplicates() %>%
    filter_poi(Classification != 7, Classification != 18) %>% # Remove errors
    normalize_height(tin())
  
  metrics <- pixel_metrics(las,
                  # func = ~ las_cld_metrics(z = Z, r = ReturnNumber),
                  func = ~ las_ladder_evenness(z = Z, r = ReturnNumber),
                  res = fire_rast)
  
  for (i in names(metrics)) {
    
    layer <- metrics[[i]]
    output <- str_replace(las_i, '\\.las', glue('_{i}.tif')) %>%
      str_replace('las/als_2022', 'rast_metrics')
    
    writeRaster(layer, output, overwrite = T)
    
    message('Exporting: ', output)
    
  }

}
