# mosquito_fire — Python‑First LiDAR Fuel‑Metrics Pipeline (Conda + stepped scripts)

A clean, current workflow to: (1) download RockyWeb LiDAR, (2) build DTMs and normalized point clouds, (3) acquire/align Mosquito Fire burn‑severity from FRAP, (4) run **your** lidR fuel‑metrics script, and (5) package results as a Cloud‑Optimized GeoTIFF (COG) + minimal STAC.

**Assumptions**
- Data root (external drive): `/mnt/data/mosquito_fire`
- CRS: **EPSG:6339** — NAD83(2011) / UTM 10N
- You’re using **Conda/Mamba** on Ubuntu

## Quickstart
```bash
mamba env create -f environment.yml
conda activate mosquito_fire

# Fill 00_manifests/wesm.csv under your data root with RockyWeb URLs (or identifiers with a URL column)
# Place FRAP Mosquito Fire raster in 10_burn/ or set `frap_url` in config.yaml

python scripts/01_download_lidar.py
python scripts/02_build_dtm_and_normalize.py
python scripts/03_get_burn_frap.py
Rscript scripts/04_run_fuel_metrics.R   # paste your lidR code into the script
python scripts/05_package_outputs.py
```

**Data directories expected (on your external drive):**
```
/mnt/data/mosquito_fire/
  00_manifests/   # wesm.csv + URL list
  01_raw/         # downloaded *.laz/*.las
  02_normalized/  # normalized LAZ with HAG
  03_derivatives/ # DTM tiles + dtm_mosaic.tif
  04_metrics/     # final stack outputs
  06_stac/        # STAC catalog
  10_burn/        # FRAP burn‑severity (input + aligned)
```

## Notes
- Keep CRS consistent (EPSG:6339). 
- Use an external SSD for best performance.
- Set `frap_url` in `config.yaml` to auto-download FRAP, or drop a local file named in `burn_filename`.
