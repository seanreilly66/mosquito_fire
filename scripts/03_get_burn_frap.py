#!/usr/bin/env python3
import sys, subprocess
from pathlib import Path
import yaml
import rioxarray as rxr

cfg = yaml.safe_load(open('config.yaml'))
DATA = Path(cfg['data_root'])
CRS  = cfg['crs']
RES  = float(cfg['resolution'])
RESAMP = cfg.get('burn_resampling','nearest').lower()
BURN_DIR = DATA/'10_burn'
BURN_DIR.mkdir(parents=True, exist_ok=True)

src_name = cfg.get('burn_filename','mosquito_burn_severity.tif')
BURN_IN  = BURN_DIR/src_name
BURN_OUT = BURN_DIR/'burn_severity_utm.tif'

url = cfg.get('frap_url','').strip()
if url and not BURN_IN.exists():
    print(f"Downloading FRAP raster from {url} → {BURN_IN}")
    subprocess.run(['aria2c','-c','-x8','-s8','-o',str(BURN_IN.name),url], cwd=str(BURN_DIR), check=True)

if not BURN_IN.exists():
    sys.exit(f"Missing burn raster: {BURN_IN}")

print("Reprojecting/snapping burn raster…")
arr = rxr.open_rasterio(BURN_IN, masked=True)
method = {'nearest':'nearest','bilinear':'bilinear'}.get(RESAMP,'nearest')
arr = arr.rio.reproject(CRS, resolution=RES, resampling=method)
arr.rio.to_raster(BURN_OUT, compress='LZW', tiled=True)
print(f"Aligned burn raster → {BURN_OUT}")
