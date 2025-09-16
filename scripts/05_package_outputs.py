#!/usr/bin/env python3
import hashlib, subprocess, datetime as dt
from pathlib import Path
import yaml, pystac

cfg = yaml.safe_load(open('config.yaml'))
DATA = Path(cfg['data_root'])
OUT  = DATA/'04_metrics'
stack_in  = OUT/'fuel_metrics_stack.tif'
stack_cog = OUT/'fuel_metrics_stack_cog.tif'
if not stack_in.exists():
    raise SystemExit(f"Expected output not found: {stack_in}")

subprocess.run(['gdal_translate','-of','COG','-co','COMPRESS=LZW','-co','RESAMPLING=AVERAGE',
                str(stack_in), str(stack_cog)], check=True)
print(f"COG → {stack_cog}")

sha256 = hashlib.sha256(stack_cog.read_bytes()).hexdigest()
cat_dir = DATA/'06_stac'; cat_dir.mkdir(parents=True, exist_ok=True)
collection = pystac.Collection(
    id=f"{cfg['project_name']}_fuelmetrics_v1",
    description=f"Fuel metrics + burn severity for {cfg['project_name']}",
    extent=pystac.Extent.from_dict({"spatial":{"bbox":[[-180,-90,180,90]]},"temporal":{"interval":[[None,None]]}}),
    license="CC-BY-4.0")
item = pystac.Item(id='fuel_metrics_stack', geometry=None, bbox=None, datetime=dt.datetime.utcnow(), properties={
    "project": cfg['project_name'], "crs": cfg['crs'], "sha256": sha256
})
item.add_asset('stack', pystac.Asset(href=str(stack_cog.relative_to(cat_dir.parent)), roles=['data'], media_type='image/tiff; application=geotiff; profile=cloud-optimized'))
collection.add_item(item)
collection.normalize_hrefs(str(cat_dir))
collection.save(catalog_type=pystac.CatalogType.SELF_CONTAINED)
print(f"STAC → {cat_dir}")
