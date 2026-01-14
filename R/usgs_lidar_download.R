# ==============================================================================
#
# ALS data download
#
# ==============================================================================
#
# Author: Sean Reilly, sean.reilly66@gmail.com
#
# Created: 20 May 2025
#
# Status: Needs documentation
#
# ==============================================================================
#
# Description:
#
# Downloads ALS files from rockyweb for fire event 
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
library(USGSlidar)

# ================================= User inputs ================================

parent_url <- 'https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Elevation/LPC/Projects/CA_SierraNevada_B22/'

fire_perim_gdb <- 'data/frap_perim/fire24_1.gdb'
als_tiles_folder <- 'data/als/lidar_tiles'
als_bndry_file <- 'data/als/USGS_CA_Sierra_Nevada_22_bndry.shp'

fire_name <- 'MOSQUITO'
fire_year <- 2022 

# ================================= Data prep ================================== 

als_bndry <- st_read(als_bndry_file)

fire_perim <- fire_perim_gdb %>%
  st_read('firep24_1') %>%
  st_transform(st_crs(als_bndry)) %>%
  filter(FIRE_NAME == fire_name) %>%
  filter(YEAR_ == fire_year)

read_shp <- function(shp_file, fire_perim) {
  shp = shp_file %>%
    st_read() %>%
    st_transform(st_crs(fire_perim)) %>%
    st_filter(fire_perim) %>%
    rename_with(tolower)
}

als_tiles <- list.files(als_tiles_folder, 
                        pattern = '\\.shp$',
                        full.names = T) %>%
  lapply(FUN = read_shp, fire_perim = fire_perim) %>%
  bind_rows()

# ================================= Rockeyweb setup ================================== 

options(timeout = max(5000, getOption("timeout")))
fetchUSGSProjectIndex('data/als/usgs_lidar_index/wesm_proj_index.gpkg')
fetchUSGSTileIndex('data/als/usgs_lidar_index/wesm_tile_index.gpkg')
