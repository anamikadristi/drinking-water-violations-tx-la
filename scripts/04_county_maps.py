"""Step 4: county summaries, choropleth maps, and a GeoPackage for ArcGIS Pro / QGIS.

Each community water system is assigned to the first county it lists as served
in SDWIS (its primary county). Counties with fewer than 3 systems are shown in
gray because their percentages are unstable.

Usage:
    python scripts/04_county_maps.py [--raw-dir PATH]
"""
import argparse
import os
import zipfile

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

INK, INK2, SURF, NODATA = "#0b0b0b", "#52514e", "#fcfcfb", "#d9d8d4"
BLUES = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
ORANGES = ["#fde3d6", "#f9bfa3", "#f39a72", "#eb6834", "#c4521f", "#963d14", "#6b2a0c"]
plt.rcParams.update({"figure.facecolor": SURF, "savefig.facecolor": SURF, "text.color": INK, "font.size": 9})


def choropleth(ax, gdf, col, bins, colors, title, label_fmt, small_mask):
    cmap = ListedColormap(colors[:len(bins) - 1])
    norm = BoundaryNorm(bins, cmap.N)
    show = gdf[~small_mask & gdf[col].notna()]
    show.plot(column=col, cmap=cmap, norm=norm, ax=ax, edgecolor=SURF, linewidth=0.25)
    rest = gdf[small_mask | gdf[col].isna()]
    if len(rest):
        rest.plot(ax=ax, color=NODATA, edgecolor=SURF, linewidth=0.25)
    ax.set_axis_off()
    ax.set_title(title, loc="left", fontsize=10)
    handles = [Patch(facecolor=cmap(i), label=label_fmt(bins[i], bins[i + 1])) for i in range(cmap.N)]
    if small_mask.any() or gdf[col].isna().any():
        handles.append(Patch(facecolor=NODATA, label="Fewer than 3 systems / no data"))
    ax.legend(handles=handles, loc="lower left", frameon=False, fontsize=7.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", default="raw_data")
    a = ap.parse_args()
    shp_dir = os.path.join(a.raw_dir, "cb_2022_us_county_500k")
    if not os.path.exists(shp_dir):
        with zipfile.ZipFile(os.path.join(a.raw_dir, "cb_2022_us_county_500k.zip")) as z:
            z.extractall(shp_dir)
    co = gpd.read_file(os.path.join(shp_dir, "cb_2022_us_county_500k.shp"))
    co = co[co.STATEFP.isin(["48", "22"])][["GEOID", "NAMELSAD", "STUSPS", "geometry"]].to_crs(5070)

    c = pd.read_csv("data/cws_systems_txla_2015_2025.csv", dtype={"county_fips": str})
    c["county_fips"] = c.county_fips.str.replace(r"\.0$", "", regex=True).str.zfill(5)
    g = c.dropna(subset=["county_fips"]).groupby("county_fips")
    s = pd.DataFrame({
        "n_cws": g.size(),
        "pop_served": g.population.sum(),
        "sys_years": g.exposure_years.sum(),
        "hb_viol": g.n_hb_total.sum(),
        "cws_any_hb": g.n_hb_total.apply(lambda x: (x > 0).sum()),
        "as_viol": g.n_arsenic.sum(),
        "as_cws": g.n_arsenic.apply(lambda x: (x > 0).sum()),
        "no3_viol": g.n_nitrate.sum(),
        "no3_cws": g.n_nitrate.apply(lambda x: (x > 0).sum()),
        "dbp_viol": g.n_dbp.sum(),
    })
    s["pct_any_hb"] = 100 * s.cws_any_hb / s.n_cws
    s["rate_100sy"] = 100 * s.hb_viol / s.sys_years
    cov = pd.read_csv("data/county_covariates_txla.csv", dtype={"county_fips": str})
    gdf = co.merge(s, left_on="GEOID", right_index=True, how="left").merge(
        cov[["county_fips", "RPL_THEMES", "EP_POV150", "rucc_2023", "metro"]],
        left_on="GEOID", right_on="county_fips", how="left").drop(columns="county_fips")
    for k in ["n_cws", "hb_viol", "cws_any_hb", "as_viol", "as_cws", "no3_viol", "no3_cws", "dbp_viol"]:
        gdf[k] = gdf[k].fillna(0).astype(int)
    gdf = gdf.rename(columns={"RPL_THEMES": "svi_2022"})
    num = gdf.select_dtypes("number").columns
    gdf[num] = gdf[num].round(3)

    os.makedirs("outputs/gis", exist_ok=True)
    gdf.drop(columns="geometry").to_csv("outputs/tables/t11_county_summary.csv", index=False)
    out = gdf.copy()
    out["geometry"] = out.geometry.simplify(100)  # 100 m tolerance keeps the file small
    out.to_file("outputs/gis/county_violation_summary_txla.gpkg", layer="county_summary", driver="GPKG")
    # KML (WGS84) for Google Earth
    kml = out.to_crs(4326)[["GEOID", "NAMELSAD", "STUSPS", "n_cws", "pct_any_hb", "as_viol", "no3_viol", "svi_2022", "geometry"]]
    kml.to_file("outputs/gis/county_violation_summary_txla.kml", driver="KML")

    small = gdf.n_cws < 3
    pct_fmt = lambda lo, hi: f"{lo:.0f} to {hi:.0f}%"
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    choropleth(axes[0], gdf, "pct_any_hb", [0, 10, 20, 30, 50, 70, 100.01], BLUES,
               "% of community water systems with >=1\nhealth-based violation, 2015-2025", pct_fmt, small)
    choropleth(axes[1], gdf, "svi_2022", [0, .2, .4, .6, .8, 1.0001], ORANGES,
               "CDC/ATSDR SVI 2022, overall percentile\n(national county ranking)",
               lambda lo, hi: f"{lo:.1f} to {min(hi, 1):.1f}", pd.Series(False, index=gdf.index))
    fig.text(0.01, 0.01, "Systems assigned to their first listed county served. Sources: EPA SDWIS via ECHO (Oct 2026 extract), CDC/ATSDR SVI 2022, Census 2022 county boundaries.",
             fontsize=7, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig("outputs/figures/fig06_map_any_violation_and_svi.png", dpi=170)
    plt.close(fig)

    cnt_fmt = lambda lo, hi: f"{lo:.0f}" if hi - lo <= 1 else f"{lo:.0f} to {hi - 1:.0f}"
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    none = pd.Series(False, index=gdf.index)
    choropleth(axes[0], gdf, "as_viol", [0, 1, 10, 50, 150, 400], BLUES,
               "Arsenic MCL violations, 2015-2025", cnt_fmt, none)
    choropleth(axes[1], gdf, "no3_viol", [0, 1, 10, 50, 150, 300], BLUES,
               "Nitrate/nitrite MCL violations, 2015-2025", cnt_fmt, none)
    fig.text(0.01, 0.01, "Counts of violation records (one system can have many quarterly violations). Source: EPA SDWIS via ECHO.",
             fontsize=7, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig("outputs/figures/fig07_map_arsenic_nitrate.png", dpi=170)
    plt.close(fig)

    # county-level correlation (descriptive only)
    ok = gdf[gdf.n_cws >= 3]
    rho = ok[["pct_any_hb", "svi_2022"]].corr(method="spearman").iloc[0, 1]
    print(f"Counties with >=3 CWS: {len(ok)}; Spearman rho(pct_any_hb, SVI) = {rho:.3f}")
    for st in ["TX", "LA"]:
        o = ok[ok.STUSPS == st]
        print(f"  {st}: n={len(o)}, rho={o[['pct_any_hb', 'svi_2022']].corr(method='spearman').iloc[0, 1]:.3f}")
    print("Counties with any arsenic violation:", (gdf.as_viol > 0).sum(), "; any nitrate:", (gdf.no3_viol > 0).sum())


if __name__ == "__main__":
    main()
