"""Step 3: descriptive rates, regression models, and the nitrate/arsenic section.

Reads the small tables in data/ (built by step 2) and writes
outputs/tables/*.csv and outputs/figures/*.png.

Usage:
    python scripts/03_analysis.py
"""
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore", category=FutureWarning)
T, F = "outputs/tables", "outputs/figures"
os.makedirs(T, exist_ok=True)
os.makedirs(F, exist_ok=True)

# palette (validated categorical slots 1-2; sequential blue; neutral ink)
TX_C, LA_C = "#2a78d6", "#eb6834"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
})
SIZE_ORDER = ["Very small (<=500)", "Small (501-3,300)", "Medium (3,301-10,000)",
              "Large (10,001-100,000)", "Very large (>100,000)", "Population not reported (0)"]
SRC_ORDER = ["Groundwater", "Groundwater (purchased)", "Surface water", "Surface water (purchased)"]
GROUPS = [("dbp", "Disinfection byproducts"), ("treatment", "Surface/ground water treatment"),
          ("rad", "Radionuclides"), ("arsenic", "Arsenic"), ("other_ioc", "Other inorganics"),
          ("nitrate", "Nitrate/nitrite"), ("coliform", "Coliform (TCR/RTCR)"),
          ("leadcopper", "Lead and copper"), ("organics", "Organic chemicals")]

c = pd.read_csv("data/cws_systems_txla_2015_2025.csv", dtype={"county_fips": str})
c["county_fips"] = c.county_fips.str.replace(r"\.0$", "", regex=True).str.zfill(5)
c.loc[c.county_fips.str.contains("nan", na=True), "county_fips"] = np.nan
c["size_class"] = c.size_class.fillna("Population not reported (0)")
v = pd.read_csv("data/hb_violations_cws_txla_2015_2025.csv", dtype=str)
v["year"] = v.year.astype(int)
cov = pd.read_csv("data/county_covariates_txla.csv", dtype={"county_fips": str})
c = c.merge(cov[["county_fips", "RPL_THEMES", "EP_POV150", "rucc_2023", "metro"]], on="county_fips", how="left")
for k, _ in GROUPS + [("lcrr_inventory", "")]:
    if "n_" + k not in c:
        c["n_" + k] = 0


def rate_table(df, by, count="n_hb_total"):
    g = df.groupby(by, observed=True)
    out = pd.DataFrame({
        "systems": g.size(),
        "population_served": g.population.sum(),
        "system_years": g.exposure_years.sum().round(1),
        "hb_violations": g[count].sum(),
        "systems_with_any": g[count].apply(lambda s: (s > 0).sum()),
    })
    out["rate_per_100_system_years"] = (100 * out.hb_violations / out.system_years).round(2)
    out["pct_systems_with_any"] = (100 * out.systems_with_any / out.systems).round(1)
    out["population_at_systems_with_any"] = g.apply(lambda d: d.loc[d[count] > 0, "population"].sum())
    out["pct_population_at_systems_with_any"] = (100 * out.population_at_systems_with_any / out.population_served).round(1)
    return out.reset_index()


# ---------- 1. descriptive rates
c["All"] = "All"
t_state = rate_table(c, ["state"])
t_all = rate_table(c, ["All"]).rename(columns={"All": "state"})
pd.concat([t_state, t_all]).to_csv(f"{T}/t01_rates_by_state.csv", index=False)
t_size = rate_table(c, ["state", "size_class"])
t_size["size_class"] = pd.Categorical(t_size.size_class, SIZE_ORDER)
t_size.sort_values(["state", "size_class"]).to_csv(f"{T}/t02_rates_by_state_size.csv", index=False)
t_src = rate_table(c, ["state", "source"])
t_src.to_csv(f"{T}/t03_rates_by_state_source.csv", index=False)

rows = []
for k, lab in GROUPS + [("lcrr_inventory", "LCRR service line inventory (reported separately)")]:
    for st, d in list(c.groupby("state")) + [("TX+LA", c)]:
        r = rate_table(d.assign(x=1), ["x"], count="n_" + k).iloc[0]
        rows.append({"group": lab, "state": st, "hb_violations": r.hb_violations,
                     "systems_with_any": r.systems_with_any, "pct_systems_with_any": r.pct_systems_with_any,
                     "rate_per_100_system_years": r.rate_per_100_system_years,
                     "population_at_systems_with_any": r.population_at_systems_with_any,
                     "pct_population_at_systems_with_any": r.pct_population_at_systems_with_any})
t_grp = pd.DataFrame(rows)
t_grp.to_csv(f"{T}/t04_rates_by_contaminant_group_state.csv", index=False)

# group by source (TX+LA pooled)
rows = []
for k, lab in GROUPS:
    for s, d in c.groupby("source"):
        rows.append({"group": lab, "source": s, "systems": len(d), "hb_violations": d["n_" + k].sum(),
                     "pct_systems_with_any": round(100 * (d["n_" + k] > 0).mean(), 1),
                     "rate_per_100_system_years": round(100 * d["n_" + k].sum() / d.exposure_years.sum(), 2)})
pd.DataFrame(rows).to_csv(f"{T}/t05_rates_by_contaminant_group_source.csv", index=False)

ann = v[v.group != "LCRR service line inventory"].pivot_table(index="year", columns=["state", "group"],
                                                               values="VIOLATION_ID", aggfunc="count", fill_value=0)
ann.to_csv(f"{T}/t06_annual_hb_violations_by_state_group.csv")

# ---------- figures: size and group
fig, ax = plt.subplots(figsize=(8, 4.2))
sizes = SIZE_ORDER[:5]
x = np.arange(len(sizes))
for i, (st, col) in enumerate([("TX", TX_C), ("LA", LA_C)]):
    d = t_size[t_size.state == st].set_index("size_class").reindex(sizes)
    b = ax.bar(x + (i - 0.5) * 0.38, d.pct_systems_with_any, 0.36, color=col, label=st,
               edgecolor=SURF, linewidth=1)
    for xi, val, n in zip(x, d.pct_systems_with_any, d.systems):
        ax.text(xi + (i - 0.5) * 0.38, val + 1, f"{val:.0f}%\nn={n}", ha="center", va="bottom", fontsize=7.5, color=INK2)
ax.set_xticks(x, [s.replace(" (", "\n(") for s in sizes])
ax.set_ylabel("% of community water systems with\n>=1 health-based violation, 2015-2025")
ax.set_ylim(0, 100)
ax.legend(frameon=False, title="State")
ax.set_title("Share of systems with a health-based violation, by system size", loc="left", fontsize=11)
fig.tight_layout()
fig.savefig(f"{F}/fig01_any_violation_by_size_state.png", dpi=170)
plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 4.6))
labs = [lab for _, lab in GROUPS]
y = np.arange(len(labs))
for i, (st, col) in enumerate([("TX", TX_C), ("LA", LA_C)]):
    d = t_grp[t_grp.state == st].set_index("group").reindex(labs)
    ax.barh(y + (i - 0.5) * 0.38, d.pct_systems_with_any, 0.36, color=col, label=st, edgecolor=SURF, linewidth=1)
    for yi, val in zip(y, d.pct_systems_with_any):
        ax.text(val + 0.4, yi + (i - 0.5) * 0.38, f"{val:.1f}", va="center", fontsize=7.5, color=INK2)
ax.set_yticks(y, labs)
ax.invert_yaxis()
ax.set_xlabel("% of systems with >=1 health-based violation in the group, 2015-2025")
ax.legend(frameon=False, title="State", loc="lower right")
ax.set_title("Health-based violations by contaminant or rule group", loc="left", fontsize=11)
fig.tight_layout()
fig.savefig(f"{F}/fig02_contaminant_groups_by_state.png", dpi=170)
plt.close(fig)

fig, axes = plt.subplots(2, 3, figsize=(9, 5), sharex=True)
for ax, (_, lab) in zip(axes.flat, GROUPS[:6]):
    for st, col in [("TX", TX_C), ("LA", LA_C)]:
        s = ann.get((st, lab), pd.Series(0, index=ann.index))
        ax.plot(ann.index, s, color=col, lw=2, marker="o", ms=3.5, label=st)
    ax.set_title(lab, fontsize=9.5, loc="left")
    ax.set_ylim(bottom=0)
    ax.set_xticks([2015, 2018, 2021, 2024])
axes[0, 0].legend(frameon=False, fontsize=8)
fig.suptitle("Health-based violations per year (compliance period start), 2015-2025", x=0.01, ha="left", fontsize=11)
fig.text(0.01, 0.005, "2025 may be incomplete because of reporting lag.", fontsize=8, color=INK2)
fig.tight_layout(rect=(0, 0.02, 1, 1))
fig.savefig(f"{F}/fig03_annual_trends_by_group.png", dpi=170)
plt.close(fig)

# most common violation types within each group and state
t12 = (v.groupby(["state", "group", "contaminant", "VIOLATION_CODE", "violation_desc"]).size()
       .reset_index(name="violations").sort_values(["state", "violations"], ascending=[True, False]))
t12.groupby("state").head(12).to_csv(f"{T}/t12_top_violation_types_by_state.csv", index=False)

# ---------- 2. models
m = c[(c.population > 0) & c.RPL_THEMES.notna()].copy()
m["size3"] = pd.Categorical(m.size_class, SIZE_ORDER[:5])
m["source"] = pd.Categorical(m.source, SRC_ORDER)
m["state"] = pd.Categorical(m.state, ["TX", "LA"])
m["metro"] = pd.Categorical(m.metro, ["Metro", "Nonmetro"])
m["svi10"] = m.RPL_THEMES * 10  # one decile of national county SVI rank
m["log_exp"] = np.log(m.exposure_years)
rhs = "C(size3) + C(source) + C(state) + svi10 + C(metro)"
print(f"Model sample: {len(m):,} systems in {m.county_fips.nunique()} counties")


def tidy(res, name, expo=True):
    ci = res.conf_int()
    out = pd.DataFrame({"term": res.params.index, "estimate": res.params.values,
                        "ratio": np.exp(res.params.values), "ci_low": np.exp(ci[0].values),
                        "ci_high": np.exp(ci[1].values), "p_value": res.pvalues.values})
    out["model"] = name
    return out[~out.term.isin(["Intercept", "alpha"])]


nb = smf.negativebinomial(f"n_hb_total ~ {rhs}", data=m, offset=m.log_exp).fit(
    cov_type="cluster", cov_kwds={"groups": pd.factorize(m.county_fips)[0]}, maxiter=200, disp=0)
alpha = nb.params["alpha"]
nb_t = tidy(nb, "NB2 count, all health-based (IRR)")
full = m[m.exposure_years > 10.9].copy()
full["any_hb"] = (full.n_hb_total > 0).astype(int)
lg = smf.logit(f"any_hb ~ {rhs}", data=full).fit(
    cov_type="cluster", cov_kwds={"groups": pd.factorize(full.county_fips)[0]}, disp=0)
lg_t = tidy(lg, "Logistic, any health-based (OR)")
# sensitivity: SVI quartiles instead of linear
m["svi_q"] = pd.cut(m.RPL_THEMES, [0, .25, .5, .75, 1], include_lowest=True,
                    labels=["Q1 (least)", "Q2", "Q3", "Q4 (most)"])
nbq = smf.negativebinomial("n_hb_total ~ C(size3) + C(source) + C(state) + C(svi_q) + C(metro)",
                           data=m, offset=m.log_exp).fit(
    cov_type="cluster", cov_kwds={"groups": pd.factorize(m.county_fips)[0]}, maxiter=200, disp=0)
nbq_t = tidy(nbq, "NB2 count, SVI quartiles (IRR)")
# sensitivity: drop surface/ground water treatment-rule violations (Louisiana GWR deficiency records dominate them)
m["n_hb_no_tt"] = m.n_hb_total - m.n_treatment
nbx = smf.negativebinomial(f"n_hb_no_tt ~ {rhs}", data=m, offset=m.log_exp).fit(
    cov_type="cluster", cov_kwds={"groups": pd.factorize(m.county_fips)[0]}, maxiter=200, disp=0)
nbx_t = tidy(nbx, "NB2 count, excluding SWTR/GWR treatment-rule violations (IRR)")
# within-state models for the SVI term
st_t = []
for st in ["TX", "LA"]:
    ms = m[m.state == st]
    r = smf.negativebinomial("n_hb_total ~ C(size3) + C(source) + svi10 + C(metro)", data=ms,
                             offset=ms.log_exp).fit(cov_type="cluster",
                                                    cov_kwds={"groups": pd.factorize(ms.county_fips)[0]},
                                                    maxiter=200, disp=0)
    st_t.append(tidy(r, f"NB2 count, {st} only (IRR)"))
# group-specific logistic models for arsenic and nitrate (systems in operation whole window)
grp_t = []
full["log10pop"] = np.log10(full.population)
conv = {"NB2 main": nb.mle_retvals["converged"], "NB2 no TT": nbx.mle_retvals["converged"], "NB2 SVI quartiles": nbq.mle_retvals["converged"],
        "Logistic any": lg.mle_retvals["converged"]}
# Louisiana had no nitrate MCL violations at CWS in 2015-2025, so the nitrate model is Texas only
for k, lab, data, f in [
        ("arsenic", "Arsenic", full, "C(state) + "),
        ("nitrate", "Nitrate/nitrite (Texas only)", full[full.state == "TX"].copy(), ""),
        ("dbp", "Disinfection byproducts", full, "C(state) + ")]:
    data["y"] = (data["n_" + k] > 0).astype(int)
    r = smf.logit(f"y ~ log10pop + C(source) + {f}svi10 + C(metro)", data=data).fit(
        cov_type="cluster", cov_kwds={"groups": pd.factorize(data.county_fips)[0]}, disp=0, maxiter=200)
    conv[f"Logistic {lab}"] = r.mle_retvals["converged"]
    t = tidy(r, f"Logistic, any {lab} (OR)")
    t["n"], t["events"] = int(r.nobs), int(data.y.sum())
    grp_t.append(t)
print("Converged:", conv)
models = pd.concat([nb_t, lg_t, nbq_t, nbx_t] + st_t + grp_t)
models[["estimate", "ratio", "ci_low", "ci_high"]] = models[["estimate", "ratio", "ci_low", "ci_high"]].round(3)
models["p_value"] = models.p_value.map(lambda p: f"{p:.3g}")
models.to_csv(f"{T}/t07_model_estimates.csv", index=False)
with open(f"{T}/t07_model_notes.txt", "w") as fh:
    fh.write(f"NB2 model: n={int(nb.nobs)} systems, {m.county_fips.nunique()} counties, alpha={alpha:.3f}, "
             f"county-clustered SE, offset log(system-years).\n")
    fh.write(f"Logistic model: n={int(lg.nobs)} systems in operation for the whole 2015-2025 window, "
             f"events={int(full.any_hb.sum())}, pseudo R2={lg.prsquared:.3f}.\n")
    fh.write(f"Convergence: {conv}\n")
    fh.write("Reference levels: very small size, groundwater, Texas, metro county. svi10 = per 0.1 of national county SVI percentile rank.\n")

# forest plot of the main NB model
d = nb_t.copy()
pretty = (d.term.str.replace(r"C\(size3\)\[T\.", "Size: ", regex=True)
          .str.replace(r"C\(source\)\[T\.", "Source: ", regex=True)
          .str.replace(r"C\(state\)\[T\.", "State: ", regex=True)
          .str.replace(r"C\(metro\)\[T\.", "County: ", regex=True).str.replace("]", "", regex=False)
          .str.replace("svi10", "County SVI (+0.1 rank)", regex=False))
fig, ax = plt.subplots(figsize=(7.5, 4.4))
yy = np.arange(len(d))[::-1]
ax.hlines(yy, d.ci_low, d.ci_high, color=TX_C, lw=2)
ax.plot(d.ratio, yy, "o", color=TX_C, ms=6, mec=SURF, mew=1.5)
ax.axvline(1, color=INK2, lw=1, ls="--")
ax.set_xscale("log")
ax.set_yticks(yy, pretty)
ax.set_xlabel("Incidence rate ratio (95% CI, county-clustered), log scale")
ax.set_title("Negative binomial model of health-based violation counts", loc="left", fontsize=11)
fig.text(0.01, 0.01, "Reference levels: very small system, groundwater, Texas, metro county. Offset: log(system-years).",
         fontsize=7.5, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(f"{F}/fig04_nb_model_irr.png", dpi=170)
plt.close(fig)

# ---------- 3. nitrate and arsenic
na = v[v.group.isin(["Arsenic", "Nitrate/nitrite"])].copy()
na["measure"] = pd.to_numeric(na.VIOL_MEASURE, errors="coerce")
na.loc[na.UNIT_OF_MEASURE.str.upper() == "UG/L", "measure"] /= 1000
na["mcl"] = na.FEDERAL_MCL.str.extract(r"([\d.]+)")[0].astype(float)
na["ratio_to_mcl"] = na.measure / na.mcl
na = na.merge(c[["PWSID", "source", "size_class", "county_fips", "COUNTY_SERVED", "population"]], on="PWSID")
rows = []
for (g, st), d in na.groupby(["contaminant", "state"]):
    q = d.measure.quantile([.5, .9])
    rows.append({"contaminant": g, "state": st, "violations": len(d), "systems": d.PWSID.nunique(),
                 "counties": d.county_fips.nunique(), "median_reported_mg_L": round(q[.5], 4),
                 "p90_reported_mg_L": round(q[.9], 4), "mcl_mg_L": d.mcl.iloc[0],
                 "median_ratio_to_mcl": round(d.ratio_to_mcl.median(), 2),
                 "pct_groundwater_systems": round(100 * d.drop_duplicates("PWSID").source.str.startswith("Groundwater").mean(), 1)})
pd.DataFrame(rows).to_csv(f"{T}/t08_nitrate_arsenic_summary.csv", index=False)
top = (na.groupby(["contaminant", "state", "county_fips", "COUNTY_SERVED"])
       .agg(violations=("VIOLATION_ID", "size"), systems=("PWSID", "nunique"))
       .reset_index().sort_values(["contaminant", "violations"], ascending=[True, False]))
top.groupby("contaminant").head(10).to_csv(f"{T}/t09_nitrate_arsenic_top_counties.csv", index=False)
# all groundwater systems: share with any arsenic / nitrate violation by size
gw = c[c.source.str.startswith("Groundwater") & (c.population > 0)]
t10 = gw.groupby(["state", "size_class"]).agg(
    gw_systems=("PWSID", "size"),
    pct_any_arsenic=("n_arsenic", lambda s: round(100 * (s > 0).mean(), 1)),
    pct_any_nitrate=("n_nitrate", lambda s: round(100 * (s > 0).mean(), 1))).reset_index()
t10.to_csv(f"{T}/t10_groundwater_arsenic_nitrate_by_size.csv", index=False)

fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
for ax, g, mcl in [(axes[0], "Arsenic", 0.010), (axes[1], "Nitrate", 10)]:
    for st, col in [("TX", TX_C), ("LA", LA_C)]:
        d = na[(na.contaminant == g) & (na.state == st)].ratio_to_mcl.dropna().clip(upper=5)
        if len(d):
            ax.hist(d, bins=np.arange(1, 5.01, 0.1), color=col, alpha=0.75, label=f"{st} (n={len(d)})",
                    edgecolor=SURF, linewidth=0.5)
    ax.set_title(f"{g}: reported value / MCL ({mcl:g} mg/L)", loc="left", fontsize=10)
    ax.set_xlabel("Ratio to federal MCL (values above 5 shown at 5)")
    ax.legend(frameon=False, fontsize=8)
axes[0].set_ylabel("Violations, 2015-2025")
fig.tight_layout()
fig.savefig(f"{F}/fig05_arsenic_nitrate_ratio_to_mcl.png", dpi=170)
plt.close(fig)

# console summary
print(pd.concat([t_state, t_all]).to_string())
print(t_grp[t_grp.state == "TX+LA"].to_string())
print(nb_t.round(3).to_string())
print(lg_t.round(3).to_string())
print(open(f"{T}/t07_model_notes.txt").read())
