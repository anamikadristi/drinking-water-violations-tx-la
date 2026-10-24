"""Step 2: build the analysis tables from the TX/LA SDWIS subsets and county covariates.

Outputs (small, kept in the repo):
  data/cws_systems_txla_2015_2025.csv      one row per community water system (CWS)
  data/hb_violations_cws_txla_2015_2025.csv one row per health-based violation
  data/county_covariates_txla.csv          SVI 2022 and RUCC 2023 for TX and LA counties

Usage:
    python scripts/02_build_tables.py [--raw-dir PATH]
"""
import argparse
import os

import numpy as np
import pandas as pd

START, END = pd.Timestamp("2015-01-01"), pd.Timestamp("2025-12-31")
STATE_FIPS = {"TX": "48", "LA": "22"}

# Contaminant groups follow SDWIS RULE_CODE (EPA reference table SDWA_REF_CODE_VALUES)
RULE_GROUP = {
    "210": "Disinfection byproducts", "220": "Disinfection byproducts", "230": "Disinfection byproducts",
    "110": "Coliform (TCR/RTCR)", "111": "Coliform (TCR/RTCR)",
    "121": "Surface/ground water treatment", "122": "Surface/ground water treatment",
    "123": "Surface/ground water treatment", "140": "Surface/ground water treatment",
    "331": "Nitrate/nitrite",
    "332": "Arsenic",
    "333": "Other inorganics",
    "340": "Radionuclides",
    "350": "Lead and copper",
    "310": "Organic chemicals", "320": "Organic chemicals",
    "351": "LCRR service line inventory",
}
GROUP_SLUG = {
    "Disinfection byproducts": "dbp", "Coliform (TCR/RTCR)": "coliform",
    "Surface/ground water treatment": "treatment", "Nitrate/nitrite": "nitrate",
    "Arsenic": "arsenic", "Other inorganics": "other_ioc", "Radionuclides": "rad",
    "Lead and copper": "leadcopper", "Organic chemicals": "organics",
    "LCRR service line inventory": "lcrr_inventory",
}


def size_class(p):
    # EPA SDWIS size categories by population served
    bins = [0, 500, 3300, 10000, 100000, np.inf]
    labels = ["Very small (<=500)", "Small (501-3,300)", "Medium (3,301-10,000)",
              "Large (10,001-100,000)", "Very large (>100,000)"]
    return pd.cut(p, bins=bins, labels=labels)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", default="raw_data")
    a = ap.parse_args()
    R = a.raw_dir
    d = lambda s: pd.to_datetime(s, format="%m/%d/%Y", errors="coerce")

    # ---- water systems
    pws = pd.read_csv(os.path.join(R, "SDWA_PUB_WATER_SYSTEMS_TXLA.csv"), dtype=str)
    cws = pws[pws.PWS_TYPE_CODE == "CWS"].copy()
    cws["state"] = cws.PWSID.str[:2]
    cws["deact"] = d(cws.PWS_DEACTIVATION_DATE)
    cws["first_rep"] = d(cws.FIRST_REPORTED_DATE)
    cws["start"] = cws.first_rep.clip(lower=START).fillna(START)
    cws["end"] = cws.deact.clip(upper=END).fillna(END)
    # systems marked inactive or non-public without a deactivation date cannot be placed in time
    no_date = cws.PWS_ACTIVITY_CODE.ne("A") & cws.deact.isna()
    cws["exposure_years"] = ((cws.end - cws.start).dt.days + 1) / 365.25
    cws.loc[no_date, "exposure_years"] = np.nan
    cws["population"] = pd.to_numeric(cws.POPULATION_SERVED_COUNT, errors="coerce")
    n_all = len(cws)
    cws = cws[cws.exposure_years > 0.25].copy()
    print(f"CWS records in TX/LA: {n_all:,}; operating >= 3 months in 2015-2025: {len(cws):,}")

    src = cws.PRIMARY_SOURCE_CODE.fillna(cws.GW_SW_CODE)
    base = np.where(src.str.startswith("SW") | src.str.startswith("GU"), "Surface water", "Groundwater")
    purch = np.where(src.str.endswith("P"), " (purchased)", "")
    cws["source"] = pd.Series(base, index=cws.index) + purch
    cws["size_class"] = size_class(cws.population)

    # ---- primary county (first listed county served)
    geo = pd.read_csv(os.path.join(R, "SDWA_GEOGRAPHIC_AREAS_TXLA.csv"), dtype=str)
    geo = geo[(geo.AREA_TYPE_CODE == "CN") & geo.ANSI_ENTITY_CODE.notna()]
    geo = geo.sort_values(["PWSID", "ANSI_ENTITY_CODE"]).drop_duplicates("PWSID")
    geo["county_fips"] = geo.PWSID.str[:2].map(STATE_FIPS) + geo.ANSI_ENTITY_CODE.str.zfill(3)
    cws = cws.merge(geo[["PWSID", "county_fips", "COUNTY_SERVED"]], on="PWSID", how="left")

    # ---- county covariates
    svi = pd.read_csv(os.path.join(R, "SVI_2022_US_county.csv"), dtype={"FIPS": str, "STCNTY": str})
    svi = svi[svi.ST_ABBR.isin(["TX", "LA"])]
    keep = ["FIPS", "ST_ABBR", "COUNTY", "E_TOTPOP", "RPL_THEMES", "RPL_THEME1", "RPL_THEME2",
            "RPL_THEME3", "RPL_THEME4", "EP_POV150", "EP_MINRTY"]
    svi = svi[keep].rename(columns={"FIPS": "county_fips"})
    for c in keep[4:]:
        svi.loc[svi[c] < 0, c] = np.nan  # -999 = missing
    ru = pd.read_csv(os.path.join(R, "Ruralurbancontinuumcodes2023.csv"), dtype={"FIPS": str}, encoding="latin-1")
    ru = ru[ru.Attribute == "RUCC_2023"][["FIPS", "Value"]].rename(columns={"FIPS": "county_fips", "Value": "rucc_2023"})
    ru["rucc_2023"] = pd.to_numeric(ru.rucc_2023)
    cov = svi.merge(ru, on="county_fips", how="left")
    cov["metro"] = np.where(cov.rucc_2023 <= 3, "Metro", "Nonmetro")
    cov.to_csv("data/county_covariates_txla.csv", index=False)
    print(f"County covariates: {len(cov)} counties/parishes; SVI missing: {cov.RPL_THEMES.isna().sum()}")

    # ---- health-based violations
    v = pd.read_csv(os.path.join(R, "SDWA_VIOLATIONS_ENFORCEMENT_TXLA.csv"), dtype=str,
                    usecols=["PWSID", "VIOLATION_ID", "COMPL_PER_BEGIN_DATE", "VIOLATION_CODE",
                             "VIOLATION_CATEGORY_CODE", "IS_HEALTH_BASED_IND", "CONTAMINANT_CODE",
                             "VIOL_MEASURE", "UNIT_OF_MEASURE", "FEDERAL_MCL", "RULE_CODE",
                             "VIOLATION_STATUS"])
    v = v[v.IS_HEALTH_BASED_IND == "Y"].drop_duplicates(["PWSID", "VIOLATION_ID"])
    v["begin"] = d(v.COMPL_PER_BEGIN_DATE)
    v = v[v.begin.between(START, END) & v.PWSID.isin(cws.PWSID)].copy()
    v["year"] = v.begin.dt.year
    v["group"] = v.RULE_CODE.map(RULE_GROUP).fillna("Other")
    ref = pd.read_csv(os.path.join(R, "SDWA_REF_CODE_VALUES.csv"), dtype=str)
    cname = ref[ref.VALUE_TYPE == "CONTAMINANT_CODE"].set_index("VALUE_CODE").VALUE_DESCRIPTION
    v["contaminant"] = v.CONTAMINANT_CODE.map(cname)
    vname = ref[ref.VALUE_TYPE == "VIOLATION_CODE"].set_index("VALUE_CODE").VALUE_DESCRIPTION
    v["violation_desc"] = v.VIOLATION_CODE.map(vname)
    v = v.drop(columns=["IS_HEALTH_BASED_IND", "COMPL_PER_BEGIN_DATE"])
    v["state"] = v.PWSID.str[:2]
    v.to_csv("data/hb_violations_cws_txla_2015_2025.csv", index=False)
    print("Health-based violations at cohort CWS, 2015-2025:", len(v))
    print(v.group.value_counts().to_string())

    # ---- counts per system
    v["slug"] = v.group.map(GROUP_SLUG).fillna("other")
    cnt = v.pivot_table(index="PWSID", columns="slug", values="VIOLATION_ID", aggfunc="count", fill_value=0)
    cnt.columns = ["n_" + c for c in cnt.columns]
    cws = cws.merge(cnt, left_on="PWSID", right_index=True, how="left")
    ncols = [c for c in cws.columns if c.startswith("n_")]
    cws[ncols] = cws[ncols].fillna(0).astype(int)
    cws["n_hb_total"] = cws[[c for c in ncols if c != "n_lcrr_inventory"]].sum(axis=1)
    out = cws[["PWSID", "PWS_NAME", "state", "PWS_ACTIVITY_CODE", "OWNER_TYPE_CODE", "population",
               "size_class", "source", "county_fips", "COUNTY_SERVED", "start", "end",
               "exposure_years"] + sorted(ncols) + ["n_hb_total"]].copy()
    out["start"] = out.start.dt.date
    out["end"] = out.end.dt.date
    out["exposure_years"] = out.exposure_years.round(3)
    out.to_csv("data/cws_systems_txla_2015_2025.csv", index=False)
    print(f"Wrote {len(out):,} systems; missing county: {out.county_fips.isna().sum()}; "
          f"missing population: {out.population.isna().sum()}")


if __name__ == "__main__":
    main()
