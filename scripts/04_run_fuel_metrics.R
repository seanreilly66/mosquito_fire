#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(lidR); library(terra); library(sf); library(data.table)
})
# Read config via env (set by you when running) or default paths
DATA_ROOT <- Sys.getenv("DATA_ROOT", "/mnt/data/mosquito_fire")
CRS  <- Sys.getenv("CRS", "EPSG:6339")
las_dir <- file.path(DATA_ROOT, "02_normalized")
dtm_mos <- file.path(DATA_ROOT, "03_derivatives", "dtm_mosaic.tif")
burn_aln <- file.path(DATA_ROOT, "10_burn", "burn_severity_utm.tif")
out_dir <- file.path(DATA_ROOT, "04_metrics"); dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

cat("lidR fuel metrics — placeholder wiring\n")
cat(paste("LAS:",las_dir,"\nDTM:",dtm_mos,"\nBURN:",burn_aln,"\nOUT:",out_dir,"\n"))

# TODO: Paste your lidR pipeline here.
# Example outline:
# ctg <- readLAScatalog(las_dir); opt_cores(ctg) <- parallel::detectCores()-1
# dtm <- rast(dtm_mos); burn <- rast(burn_aln)
# ctg_norm <- normalize_height(ctg, dtm)
# metrics <- pixel_metrics(ctg_norm, ~list(p10=quantile(Z,0.1,na.rm=TRUE),
#                                          p50=quantile(Z,0.5,na.rm=TRUE),
#                                          p95=quantile(Z,0.95,na.rm=TRUE),
#                                          cover2=sum(Z>2,na.rm=TRUE)/length(Z),
#                                          mean=mean(Z,na.rm=TRUE), sd=sd(Z,na.rm=TRUE), n=length(Z)),
#                           res=res(burn)[1])
# metrics <- resample(metrics, burn, method='bilinear')
# stack <- c(burn, metrics)
# writeRaster(stack, file.path(out_dir, "fuel_metrics_stack.tif"), overwrite=TRUE,
#             gdal=c("COMPRESS=LZW","TILED=YES"))
