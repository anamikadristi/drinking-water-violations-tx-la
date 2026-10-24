#!/usr/bin/env bash
# Runs the full pipeline. Raw downloads (about 6 GB after unzipping) go to ./raw_data unless RAW_DIR is set.
set -euo pipefail
RAW_DIR="${RAW_DIR:-raw_data}"
python scripts/01_download_data.py --raw-dir "$RAW_DIR"
python scripts/02_build_tables.py --raw-dir "$RAW_DIR"
python scripts/03_analysis.py
python scripts/04_county_maps.py --raw-dir "$RAW_DIR"
