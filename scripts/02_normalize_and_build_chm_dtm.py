#!/usr/bin/env python3
"""
Step 2 — Normalize points, build DTM and CHM (separate PDAL runs), parallel per tile.

Inputs:  <data_root>/01_raw/*.las|*.laz
Outputs:
  - <data_root>/02_normalized/*_norm.las            (LAS 1.4 / pointfmt 6, uncompressed, keeps HeightAboveGround)
  - <data_root>/03_derivatives/*_dtm.tif            (per-tile DTM GeoTIFF)
  - <data_root>/03_derivatives/*_chm.tif            (per-tile CHM GeoTIFF, smoothed if enabled)
  - <data_root>/03_derivatives/{dtm,chm}.vrt        (virtual mosaics)
  - <data_root>/03_derivatives/{dtm,chm}_mosaic.tif (COG mosaics)

Config (config.yaml):
  data_root: /mnt/d/mosquito_fire
  crs: EPSG:6339
  dtm_resolution: 1.0
  make_chm: true
  chm_resolution: 1.0
  chm_smooth: true
"""
import argparse
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import rasterio
import yaml
from numpy.lib.stride_tricks import sliding_window_view
from scipy.ndimage import median_filter  # used only if you switch to a simpler median
from tqdm import tqdm


# ---------------- utilities ----------------

def run(cmd, **kw):
    kw.setdefault("check", True)
    return subprocess.run(cmd, **kw)


def smooth_chm_fill_pits(infile: Path):
    """
    3×3 median smoothing that IGNORES NoData and FILLS pits:
      - If center is NoData but neighbors exist, replace with median(neighbors)
      - Only keep NoData if all 3×3 values are NoData
    """
    with rasterio.open(infile) as src:
        band = src.read(1)
        meta = src.meta.copy()
        nodata = src.nodata if src.nodata is not None else -9999

    arr = band.astype(np.float32)
    # turn NoData into NaN for nan-aware median
    arr[arr == nodata] = np.nan

    # Sliding 3x3 windows (shape -> (H-2, W-2, 3, 3))
    win = sliding_window_view(arr, (3, 3))
    med = np.nanmedian(win, axis=(-2, -1))  # median ignoring NaN

    # Pad back to original shape
    smoothed = np.pad(med, ((1, 1), (1, 1)), mode="edge")

    # Where even neighbors were all NaN, keep nodata
    smoothed = np.where(np.isnan(smoothed), nodata, smoothed)

    meta.update(dtype="float32", compress="lzw", tiled=True)
    with rasterio.open(infile, "w", **meta) as dst:
        dst.write(smoothed.astype(np.float32), 1)


# ------------- PDAL pipeline builders (one writer each) -------------

def build_pdal_dtm_pipeline(in_laz, out_tif, crs, dtm_res):
    return {
        "pipeline": [
            {"type": "readers.las", "filename": str(in_laz)},
            {"type": "filters.reprojection", "out_srs": crs},
            {"type": "filters.smrf",
             "scalar": 1.2, "slope": 0.2, "threshold": 0.45, "window": 16},
            {"type": "filters.range", "limits": "Classification[2:2]"},
            {"type": "writers.gdal",
             "filename": str(out_tif),
             "gdaldriver": "GTiff",
             "output_type": "idw",
             "resolution": float(dtm_res),
             "nodata": -9999,
             "override_srs": crs,
             "gdalopts": ["COMPRESS=LZW", "TILED=YES"]}
        ]
    }


def build_pdal_norm_pipeline(in_laz, out_las, crs):
    return {
        "pipeline": [
            {"type": "readers.las", "filename": str(in_laz)},
            {"type": "filters.reprojection", "out_srs": crs},
            {"type": "filters.hag_delaunay"},
            {"type": "writers.las",
             "filename": str(out_las),     # .las (uncompressed for fast re-reads)
             "minor_version": 4,
             "dataformat_id": 6,
             # persist HAG as Extra Byte so downstream steps can use it
             "extra_dims": "HeightAboveGround=float32"}
        ]
    }


def build_pdal_chm_from_norm_pipeline(in_norm_las, out_tif, crs, chm_res):
    return {
        "pipeline": [
            {"type": "readers.las", "filename": str(in_norm_las)},
            {"type": "filters.reprojection", "out_srs": crs},
            {"type": "writers.gdal",
             "filename": str(out_tif),
             "gdaldriver": "GTiff",
             "resolution": float(chm_res),
             "output_type": "max",
             "dimension": "HeightAboveGround",
             "nodata": -9999,              # keep NoData so smoothing can fill intelligently
             "override_srs": crs,
             "gdalopts": ["COMPRESS=LZW", "TILED=YES"]}
        ]
    }


# ---------------- per-tile worker ----------------

def process_tile(tile_path: Path, norm_dir: Path, der_dir: Path, crs: str,
                 dtm_res: float, make_chm: bool, chm_res: float, chm_smooth: bool, force: bool):
    stem = tile_path.stem
    out_dtm  = der_dir / f"{stem}_dtm.tif"
    out_norm = norm_dir / f"{stem}_norm.las"   # uncompressed LAS
    out_chm  = der_dir / f"{stem}_chm.tif"

    # Skip if everything we need already exists
    if (not force
        and out_dtm.exists()
        and out_norm.exists()
        and (not make_chm or out_chm.exists())):
        return {"tile": tile_path.name, "status": "skipped"}

    # 1) DTM
    if force or not out_dtm.exists():
        dtm_pipe = build_pdal_dtm_pipeline(tile_path, out_dtm, crs, dtm_res)
        run(["pdal", "pipeline", "--stdin"], input=json.dumps(dtm_pipe).encode())

    # 2) Normalize to LAS with HAG persisted
    if force or not out_norm.exists():
        norm_pipe = build_pdal_norm_pipeline(tile_path, out_norm, crs)
        run(["pdal", "pipeline", "--stdin"], input=json.dumps(norm_pipe).encode())

    # 3) CHM from normalized LAS
    if make_chm and (force or not out_chm.exists()):
        chm_pipe = build_pdal_chm_from_norm_pipeline(out_norm, out_chm, crs, chm_res)
        run(["pdal", "pipeline", "--stdin"], input=json.dumps(chm_pipe).encode())
        if chm_smooth:
            smooth_chm_fill_pits(out_chm)

    return {"tile": tile_path.name, "status": "done"}


# ---------------- mosaics ----------------

def build_vrt_and_cog(pattern: str, der_dir: Path, out_base: str):
    tiles = sorted(der_dir.glob(pattern))
    if not tiles:
        return None, None
    vrt = der_dir / f"{out_base}.vrt"
    cog = der_dir / f"{out_base}_mosaic.tif"
    run(["gdalbuildvrt", str(vrt)] + [str(t) for t in tiles])
    run([
        "gdal_translate", "-of", "COG",
        "-co", "COMPRESS=LZW",
        "-co", "RESAMPLING=AVERAGE",
        str(vrt), str(cog)
    ])
    return vrt, cog


# ---------------- main ----------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None, help="Parallel workers (default ~= half CPUs, max 8)")
    ap.add_argument("--force", action="store_true", help="Recompute even if outputs exist")
    args = ap.parse_args()

    cfg = yaml.safe_load(open("config.yaml"))
    data_root = Path(cfg["data_root"])
    crs = cfg["crs"]
    dtm_res = float(cfg.get("dtm_resolution", 1.0))
    make_chm = bool(cfg.get("make_chm", True))
    chm_res = float(cfg.get("chm_resolution", dtm_res))
    chm_smooth = bool(cfg.get("chm_smooth", False))

    raw_dir = data_root / "01_raw"
    norm_dir = data_root / "02_normalized"
    der_dir = data_root / "03_derivatives"
    for d in (norm_dir, der_dir):
        d.mkdir(parents=True, exist_ok=True)

    laz_files = sorted(list(raw_dir.glob("*.la?")))
    if not laz_files:
        raise SystemExit(f"No LAS/LAZ found in {raw_dir}")

    cpu = os.cpu_count() or 4
    workers = args.workers or min(8, max(2, cpu // 2))

    print(f"Processing {len(laz_files)} tiles with {workers} workers…")
    results = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {
            ex.submit(process_tile, p, norm_dir, der_dir, crs,
                      dtm_res, make_chm, chm_res, chm_smooth, args.force): p
            for p in laz_files
        }
        for fut in tqdm(as_completed(futs), total=len(futs), desc="Tiles"):
            try:
                results.append(fut.result())
            except subprocess.CalledProcessError:
                results.append({"tile": futs[fut].name, "status": "error"})

    done = sum(1 for r in results if r["status"] == "done")
    skipped = sum(1 for r in results if r["status"] == "skipped")
    errors = [r for r in results if r["status"] == "error"]
    print(f"Completed: {done}, Skipped: {skipped}, Errors: {len(errors)}")

    # Build mosaics
    print("Building DTM mosaic…")
    dtm_vrt, dtm_cog = build_vrt_and_cog("*_dtm.tif", der_dir, "dtm")
    if dtm_vrt:
        print(f"DTM VRT  → {dtm_vrt}\nDTM COG  → {dtm_cog}")

    if make_chm:
        print("Building CHM mosaic…")
        chm_vrt, chm_cog = build_vrt_and_cog("*_chm.tif", der_dir, "chm")
        if chm_vrt:
            print(f"CHM VRT  → {chm_vrt}\nCHM COG  → {chm_cog}")

    if errors:
        print("Some tiles failed — rerun with --force after addressing issues.")


if __name__ == "__main__":
    main()
