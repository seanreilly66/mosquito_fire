#!/usr/bin/env python3
import sys, subprocess, re, io
from pathlib import Path
import pandas as pd
import yaml
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

# ---------------- helpers ----------------

def get_cfg_col(df, preferred_names, fallback=None):
    # Return first existing column name matching (case-insensitive) any in preferred_names
    lower = {c.lower(): c for c in df.columns}
    for name in preferred_names:
        if name.lower() in lower:
            return lower[name.lower()]
    return fallback

def http_get_text(url: str, timeout=60):
    # fetch text file with a real UA to avoid some 403s
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req, timeout=timeout) as r:
        data = r.read()
    # USGS link lists are simple ASCII/UTF-8
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError:
        return data.decode("latin-1")

def normalize_tile_id(s: str) -> str:
    if s is None:
        return ""
    s = str(s).strip()
    # drop extension if someone pasted a filename
    s = re.sub(r"\.(las|laz|zip)$", "", s, flags=re.IGNORECASE)
    # keep alnum only; tile IDs like 10SFJ7740 remain
    s = re.sub(r"[^0-9A-Za-z]+", "", s)
    return s.upper()

def url_matches_tile(url: str, tile: str) -> bool:
    """
    Match strategy:
      - extract basename (e.g., USGS_LPC_CA_SierraNevada_B22_10SFJ7740.laz)
      - match if tile appears before extension as its own token or suffix
    """
    name = re.split(r"[\\/]", url.split("?")[0])[-1]
    base = re.sub(r"\.(las|laz|zip)$", "", name, flags=re.IGNORECASE)
    return tile in base.upper()

# ---------------- load config & IO ----------------

cfg = yaml.safe_load(open("config.yaml"))
DATA = Path(cfg["data_root"])
MAN  = DATA / "00_manifests"
RAW  = DATA / "01_raw"
REP  = MAN / "reports"
for d in (MAN, RAW, REP):
    d.mkdir(parents=True, exist_ok=True)

wesm_path = Path(cfg.get("wesm_file", MAN / "WESM.csv"))
sel_path  = Path(cfg.get("selection_file", MAN / "mosquito_usgs_ca_sierra_nevada_tiles.csv"))
if not wesm_path.exists():
    sys.exit(f"Missing WESM file: {wesm_path}")
if not sel_path.exists():
    sys.exit(f"Missing selection file: {sel_path}")

wesm_wu_col   = cfg.get("wesm_workunit_col", "WU_NAME")
wesm_link_col = cfg.get("wesm_lpc_link_col", "lpc_link")
sel_wu_col    = cfg.get("sel_workunit_col", "WU_NAME")
sel_id_col    = cfg.get("sel_tileid_col", "Tile_ID")

# ---------------- read tables ----------------

print(f"Reading WESM: {wesm_path}")
wesm = pd.read_csv(wesm_path)
# resolve flexible column names if needed
wesm_wu_col   = get_cfg_col(wesm, [wesm_wu_col, "workunit", "wu_name"])
wesm_link_col = get_cfg_col(wesm, [wesm_link_col, "lpc_link", "link", "url"])
if wesm_wu_col is None or wesm_link_col is None:
    sys.exit(f"Could not find required columns in WESM (got WU={wesm_wu_col}, LINK={wesm_link_col})")

print(f"Reading selection: {sel_path}")
sel = pd.read_csv(sel_path)
sel_wu_col = get_cfg_col(sel, [sel_wu_col, "workunit", "wu_name"])
sel_id_col = get_cfg_col(sel, [sel_id_col, "tile_id", "tile", "name", "filename"])
if sel_wu_col is None or sel_id_col is None:
    sys.exit(f"Could not find required columns in selection (got WU={sel_wu_col}, ID={sel_id_col})")

# normalize tiles per work unit
sel["_tile_norm"] = sel[sel_id_col].apply(normalize_tile_id)
sel["_wu_norm"] = sel[sel_wu_col].astype(str).str.strip()

# build mapping: workunit -> set(tile_ids)
want_by_wu = (
    sel[["_wu_norm", "_tile_norm"]]
    .dropna()
    .query("_tile_norm != ''")
    .groupby("_wu_norm")["_tile_norm"]
    .apply(lambda s: set(s.unique()))
    .to_dict()
)

if not want_by_wu:
    sys.exit("No tile IDs found in selection file after normalization.")

# ---------------- per-work-unit link fetching ----------------

all_urls = []
per_wu_counts = []
missing_tiles_report = []

wesm["_wu_norm"] = wesm[wesm_wu_col].astype(str).str.strip()
wesm["_lpc_link"] = wesm[wesm_link_col].astype(str).str.strip()

# Only keep WESM rows for work units we actually want
wesm_sub = wesm[wesm["_wu_norm"].isin(want_by_wu.keys())].copy()

if wesm_sub.empty:
    sys.exit("No matching work units between WESM and selection file.")

for wu, group in wesm_sub.groupby("_wu_norm"):
    tiles_wanted = want_by_wu.get(wu, set())
    if not tiles_wanted:
        continue

    # There might be multiple rows/links per WU; loop them
    wu_urls = set()
    this_wu_missing = set(tiles_wanted)

    for link in group["_lpc_link"].dropna().unique():
        base = link.rstrip("/") + "/"
        txt_url = base + "0_file_download_links.txt"
        print(f"[{wu}] Fetching link list: {txt_url}")
        try:
            content = http_get_text(txt_url)
        except (HTTPError, URLError) as e:
            print(f"  !! Failed to fetch {txt_url}: {e}")
            continue

        # parse lines to URLs
        lines = [ln.strip() for ln in content.splitlines() if ln.strip()]
        # Some USGS files have comments or relative paths; keep only http(s)
        lines = [ln for ln in lines if ln.lower().startswith("http")]
        if not lines:
            print(f"  !! No http URLs found in {txt_url}")
            continue

        # Filter by tile IDs
        for u in lines:
            for t in list(this_wu_missing):
                if url_matches_tile(u, t):
                    wu_urls.add(u)
                    this_wu_missing.discard(t)

        print(f"  + matched {len(wu_urls)} so far (remaining tiles for this WU: {len(this_wu_missing)})")

        # Optional micro-optimization: if we’ve matched all tiles for this WU, stop scanning more rows
        if not this_wu_missing:
            break

    # record results
    per_wu_counts.append({"work_unit": wu, "tiles_requested": len(tiles_wanted), "urls_matched": len(wu_urls)})
    if this_wu_missing:
        for t in sorted(this_wu_missing):
            missing_tiles_report.append({"work_unit": wu, "missing_tile_id": t})

    all_urls.extend(sorted(wu_urls))

# ---------------- write outputs & download ----------------

urls_txt = MAN / "rockyweb_laz_urls.txt"
pd.Series(all_urls).to_csv(urls_txt, index=False, header=False)
print(f"\nWrote consolidated URL list: {urls_txt}  (total {len(all_urls)} URLs)")

pd.DataFrame(per_wu_counts).to_csv(REP / "per_workunit_counts.csv", index=False)
if missing_tiles_report:
    pd.DataFrame(missing_tiles_report).to_csv(REP / "missing_tiles_by_workunit.csv", index=False)
    print(f"Some requested tiles were not found. See: {REP / 'missing_tiles_by_workunit.csv'}")
else:
    print("All requested tiles matched from link lists.")

# Download with aria2c
if all_urls:
    conn = int(cfg.get("aria2_connections", 16))
    cmd = [
        "aria2c", "-c", f"-x{conn}", f"-s{conn}",
        "--auto-file-renaming=false",
        "-i", str(urls_txt),
        "--dir", str(RAW),
    ]
    print("Downloading with:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"Done. Files in {RAW}")
else:
    print("No URLs to download (check inputs and reports).")
