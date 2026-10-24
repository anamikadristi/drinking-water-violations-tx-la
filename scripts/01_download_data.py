"""Step 1: download public source data and subset SDWIS records to Texas and Louisiana.

Downloads (no account or API key needed):
  - EPA ECHO SDWA bulk download (SDWIS extract, about 560 MB zip, 5 GB unzipped)
  - CDC/ATSDR Social Vulnerability Index 2022, US county CSV
  - USDA ERS Rural-Urban Continuum Codes 2023
  - US Census Bureau 2022 cartographic boundary file, counties (1:500k)

Large raw files go to a raw-data directory OUTSIDE version control
(default: ./raw_data, which is listed in .gitignore). Only the Texas and
Louisiana subsets needed by later steps are written there as well.

Usage:
    python scripts/01_download_data.py [--raw-dir PATH]
"""
import argparse
import csv
import datetime as dt
import os
import urllib.request
import zipfile

SOURCES = {
    "SDWA_latest_downloads.zip": "https://echo.epa.gov/files/echodownloads/SDWA_latest_downloads.zip",
    "SVI_2022_US_county.csv": "https://svi.cdc.gov/Documents/Data/2022/csv/states_counties/SVI_2022_US_county.csv",
    "Ruralurbancontinuumcodes2023.csv": "https://www.ers.usda.gov/media/5768/2023-rural-urban-continuum-codes.csv",
    "cb_2022_us_county_500k.zip": "https://www2.census.gov/geo/tiger/GENZ2022/shp/cb_2022_us_county_500k.zip",
}
SDWA_MEMBERS = [
    "SDWA_PUB_WATER_SYSTEMS.csv",
    "SDWA_GEOGRAPHIC_AREAS.csv",
    "SDWA_VIOLATIONS_ENFORCEMENT.csv",
    "SDWA_REF_CODE_VALUES.csv",
]
STATES = ("TX", "LA")


def download(url, dest):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (research script)"})
    with urllib.request.urlopen(req) as r:
        remote_len = int(r.headers.get("Content-Length", 0))
        last_mod = r.headers.get("Last-Modified", "")
        if os.path.exists(dest) and os.path.getsize(dest) == remote_len:
            print(f"  exists, same size, skipping: {os.path.basename(dest)}")
            return remote_len, last_mod
        tmp = dest + ".part"
        with open(tmp, "wb") as f:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
        os.replace(tmp, dest)
    return os.path.getsize(dest), last_mod


def subset_by_pwsid(src, dest, pwsid_col="PWSID"):
    """Stream a large SDWA CSV and keep rows whose PWSID starts with TX or LA."""
    n_in = n_out = 0
    with open(src, newline="", encoding="utf-8", errors="replace") as fi, \
            open(dest, "w", newline="", encoding="utf-8") as fo:
        rd = csv.reader(fi)
        wr = csv.writer(fo)
        header = next(rd)
        wr.writerow(header)
        k = header.index(pwsid_col)
        for row in rd:
            n_in += 1
            if row[k][:2] in STATES:
                wr.writerow(row)
                n_out += 1
    return n_in, n_out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", default="raw_data")
    ap.add_argument("--log", default=os.path.join("data", "download_log.csv"))
    a = ap.parse_args()
    os.makedirs(a.raw_dir, exist_ok=True)
    today = dt.date.today().isoformat()
    log_rows = []
    for name, url in SOURCES.items():
        print("Downloading", url)
        size, last_mod = download(url, os.path.join(a.raw_dir, name))
        log_rows.append([name, url, size, last_mod, today])

    print("Extracting SDWA tables")
    with zipfile.ZipFile(os.path.join(a.raw_dir, "SDWA_latest_downloads.zip")) as z:
        for m in SDWA_MEMBERS:
            out = os.path.join(a.raw_dir, m)
            if not os.path.exists(out) or os.path.getsize(out) != z.getinfo(m).file_size:
                z.extract(m, a.raw_dir)

    print("Subsetting to Texas and Louisiana (this reads about 4 GB)")
    for m in ["SDWA_PUB_WATER_SYSTEMS.csv", "SDWA_GEOGRAPHIC_AREAS.csv", "SDWA_VIOLATIONS_ENFORCEMENT.csv"]:
        n_in, n_out = subset_by_pwsid(os.path.join(a.raw_dir, m),
                                      os.path.join(a.raw_dir, m.replace(".csv", "_TXLA.csv")))
        print(f"  {m}: {n_in:,} rows read, {n_out:,} TX/LA rows kept")
        log_rows.append([m.replace(".csv", "_TXLA.csv"), "subset of SDWA zip", n_out, "rows kept", today])

    os.makedirs(os.path.dirname(a.log), exist_ok=True)
    with open(a.log, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["file", "source_url", "bytes_or_rows", "last_modified_or_note", "access_date"])
        w.writerows(log_rows)
    print("Wrote", a.log)


if __name__ == "__main__":
    main()
