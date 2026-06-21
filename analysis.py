"""
Reproducibility script for "The Geography of Blockchain Infrastructure" (course edition).

Reads the assembled country-level dataset (master.gpkg) and reproduces every
statistic reported in the paper: exploratory spatial data analysis (Moran's I,
LISA with FDR correction, Getis-Ord Gi*), the concentration-vs-clustering
comparison, the negative binomial count model, and GWR/MGWR.

Run:  python analysis.py
Needs the packages in requirements.txt (geopandas, libpysal, esda, statsmodels, mgwr, scipy).

All permutation tests use seed 42. Deterministic statistics (Moran's I,
coefficients, bandwidths, concentration indices) reproduce exactly; permutation
p-values are robustly significant and may differ negligibly from the manuscript,
which chained the random state across several scripts.
"""
import os
import json
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import geopandas as gpd
import libpysal
from esda.moran import Moran, Moran_Local
from esda.getisord import G_Local
from esda import fdr
import statsmodels.api as sm
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
RES = {}
QUAD = np.array(["-", "HH", "LH", "LL", "HL"])  # esda local-Moran quadrant codes


def moran_global(y, w):
    """Global Moran's I with a fixed seed (reseed so the result is order-independent)."""
    np.random.seed(42)
    return Moran(y, w, permutations=9999)


# ---------------- data + spatial weights ----------------
m = gpd.read_file(os.path.join(HERE, "master.gpkg"))
s = m[m["complete_cov"] == 1].copy().reset_index(drop=True)   # 140-country analysis sample
cent = s.geometry.to_crs(6933).centroid                        # equal-area centroids
coords = np.array(list(zip(cent.x, cent.y)))
w = libpysal.weights.KNN.from_array(coords, k=6)               # k-nearest-neighbour weights
w.transform = "R"                                              # row-standardize
print(f"analysis sample: n={len(s)} countries\n")


# ---------------- 1. ESDA: Moran's I, LISA (FDR), Getis-Ord Gi* ----------------
def esda_block(var, label):
    y = np.log1p(s[var].values)                                # ln(1 + nodes per million)
    mi = moran_global(y, w)
    li = Moran_Local(y, w, permutations=9999, seed=42)
    cut = fdr(li.p_sim, 0.05)                                  # Benjamini-Hochberg cutoff
    sig = li.p_sim < cut
    cls = np.where(sig, QUAD[li.q], "ns")
    counts = {q: int((cls == q).sum()) for q in ["HH", "LL", "HL", "LH"]}
    gi = G_Local(y, w, star=True, permutations=9999, seed=42)
    hot = int(((gi.Zs > 0) & (gi.p_sim < 0.05)).sum())
    cold = int(((gi.Zs < 0) & (gi.p_sim < 0.05)).sum())
    RES[f"esda_{label}"] = {
        "moran_I": round(float(mi.I), 4), "moran_z": round(float(mi.z_sim), 2),
        "moran_p": float(mi.p_sim), "fdr_cutoff": round(float(cut), 4),
        "lisa_fdr_counts": counts, "gistar_hot": hot, "gistar_cold": cold,
        "HH_countries": s.loc[cls == "HH", "ne_name"].tolist(),
    }
    print(f"[ESDA {label}] Moran's I={mi.I:.3f} (z={mi.z_sim:.1f}, p={mi.p_sim}) | "
          f"FDR LISA {counts} | Gi* hot={hot} cold={cold}")

esda_block("btc_rate", "btc")
esda_block("eth_rate", "eth")
RES["btc_eth_rate_pearson"] = round(float(np.corrcoef(
    np.log1p(s["btc_rate"]), np.log1p(s["eth_rate"]))[0, 1]), 3)
print(f"  Bitcoin-Ethereum log-rate correlation r={RES['btc_eth_rate_pearson']}\n")


# ---------------- 2. concentration vs clustering ----------------
def concentration(counts):
    c = np.sort(np.asarray(counts, dtype=float))[::-1]
    shares = c / c.sum()
    hhi = float((shares ** 2).sum())
    nakamoto = int((np.cumsum(shares) < 0.5).sum() + 1)
    cs = np.sort(c)
    n = len(c)
    gini = float((2 * np.arange(1, n + 1) - n - 1).dot(cs) / (n * cs.sum()))
    return {"nakamoto": nakamoto, "gini": round(gini, 3),
            "hhi": round(hhi, 3), "top5_share": round(float(shares[:5].sum()), 3)}

RES["concentration_btc"] = concentration(m.loc[m.btc_nodes > 0, "btc_nodes"])
RES["concentration_eth"] = concentration(m.loc[m.eth_nodes > 0, "eth_nodes"])
RES["concentration_mining"] = concentration(m.loc[m.mining_share > 0, "mining_share"])

# Moran's I on log SHARES (matched transform: lets nodes be compared to mining like-for-like)
for var, lab in [("btc_nodes", "btc"), ("eth_nodes", "eth")]:
    sh = s[var].values / s[var].values.sum() * 100
    mi = moran_global(np.log1p(sh), w)
    RES[f"shares_moran_{lab}"] = {"moran_I": round(float(mi.I), 4), "p": float(mi.p_sim)}
mi_mine = moran_global(np.log1p(s["mining_share"].values * 100), w)
RES["shares_moran_mining"] = {"moran_I": round(float(mi_mine.I), 4), "p": float(mi_mine.p_sim)}
both = m[(m.mining_share > 0) | (m.btc_nodes > 0)]
rho, _ = spearmanr(both["mining_share"], both["btc_nodes"])
RES["mining_vs_nodes_spearman"] = round(float(rho), 3)
print("[Concentration] btc:", RES["concentration_btc"], "| eth:", RES["concentration_eth"],
      "| mining:", RES["concentration_mining"])
print(f"  Moran log-shares: btc={RES['shares_moran_btc']['moran_I']} "
      f"eth={RES['shares_moran_eth']['moran_I']} mining={RES['shares_moran_mining']['moran_I']} "
      f"| mining-vs-nodes Spearman rho={RES['mining_vs_nodes_spearman']}\n")


# ---------------- 3. negative binomial GLM (population offset) ----------------
d = s.copy()
d["log_gdp_pc"] = np.log(d["gdp_pc"])
d["partial_ban"] = (d["reg_status"] == "partial_ban").astype(int)
d["general_ban"] = (d["reg_status"] == "general_ban").astype(int)
X = sm.add_constant(d[["log_gdp_pc", "internet_pct", "usd_per_kwh", "partial_ban", "general_ban"]])
off = np.log(d["population"].values)

def profile_alpha_mle(y, lo=1e-3, hi=50.0):
    """Profile-likelihood MLE of the NB2 dispersion alpha (more honest here than the
    Cameron-Trivedi auxiliary regression, which is dominated by the few large-mean countries)."""
    def nll(a):
        try:
            return -sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=a), offset=off).fit().llf
        except Exception:
            return 1e12
    return float(minimize_scalar(nll, bounds=(lo, hi), method="bounded", options={"xatol": 1e-4}).x)

def nb_glm(count_var, label):
    y = d[count_var].values
    pois = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit()
    pois_disp = float(pois.pearson_chi2 / pois.df_resid)        # Poisson overdispersion check
    alpha = profile_alpha_mle(y)
    nb = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=alpha), offset=off).fit()
    np.random.seed(42)
    rm = Moran(nb.resid_deviance, w, permutations=9999)          # residual spatial autocorrelation
    RES[f"glm_{label}"] = {
        "poisson_dispersion": round(pois_disp, 1),
        "alpha_mle": round(alpha, 3),
        "pearson_chi2_df": round(float(nb.pearson_chi2 / nb.df_resid), 2),
        "deviance_explained": round(1 - nb.deviance / nb.null_deviance, 3),
        "irr": {k: round(float(np.exp(v)), 3) for k, v in nb.params.items()},
        "pvalues": {k: round(float(v), 4) for k, v in nb.pvalues.items()},
        "irr_ci95": {k: [round(float(np.exp(nb.conf_int().loc[k, 0])), 3),
                          round(float(np.exp(nb.conf_int().loc[k, 1])), 3)] for k in nb.params.index},
        "gdp_elasticity": round(float(nb.params["log_gdp_pc"]), 3),
        "resid_moran_I": round(float(rm.I), 3), "resid_moran_p": float(rm.p_sim),
    }
    print(f"[NB-GLM {label}] Poisson disp={pois_disp:.0f} -> NB alpha={alpha:.3f} | "
          f"chi2/df={RES[f'glm_{label}']['pearson_chi2_df']} devExpl={RES[f'glm_{label}']['deviance_explained']}")
    print(f"  IRR: {RES[f'glm_{label}']['irr']}")
    print(f"  p:   {RES[f'glm_{label}']['pvalues']}")
    print(f"  residual Moran's I={rm.I:.3f} (p={rm.p_sim})\n")

nb_glm("btc_nodes", "btc")
nb_glm("eth_nodes", "eth")


# ---------------- 4. GWR / MGWR ----------------
from mgwr.gwr import GWR, MGWR
from mgwr.sel_bw import Sel_BW

def gwr_block(yvar, label):
    d2 = s.copy()
    d2["log_gdp_pc"] = np.log(d2["gdp_pc"])
    yv = np.log1p(d2[yvar].values).reshape(-1, 1)
    Xv = d2[["log_gdp_pc", "internet_pct", "usd_per_kwh"]].values
    yz = (yv - yv.mean()) / yv.std()                            # standardize (z-scores)
    Xz = (Xv - Xv.mean(axis=0)) / Xv.std(axis=0)
    co = coords
    names = ["intercept", "log_gdp_pc", "internet_pct", "usd_per_kwh"]
    # global OLS
    ols = sm.OLS(yz, sm.add_constant(Xz)).fit()
    p = Xz.shape[1] + 2
    ols_aicc = 2 * p - 2 * ols.llf + (2 * p * (p + 1)) / (len(d2) - p - 1)
    # GWR (single bandwidth)
    bw = Sel_BW(co, yz, Xz, kernel="bisquare", fixed=False).search(criterion="AICc")
    gwr = GWR(co, yz, Xz, bw, kernel="bisquare", fixed=False).fit()
    # MGWR (per-covariate bandwidth)
    msel = Sel_BW(co, yz, Xz, multi=True, kernel="bisquare", fixed=False)
    msel.search(criterion="AICc", multi_bw_min=[10], multi_bw_max=[len(d2)])
    mg = MGWR(co, yz, Xz, selector=msel, kernel="bisquare", fixed=False).fit()
    tv = mg.filter_tvals()                                      # corrected-t local significance
    RES[f"gwr_{label}"] = {
        "ols_aicc": round(float(ols_aicc), 1), "ols_adjR2": round(float(ols.rsquared_adj), 3),
        "gwr_bw": int(bw), "gwr_aicc": round(float(gwr.aicc), 1),
        "gwr_adjR2": round(float(gwr.adj_R2), 3), "gwr_ENP": round(float(gwr.tr_S), 1),
        "mgwr_aicc": round(float(mg.aicc), 1), "mgwr_adjR2": round(float(mg.adj_R2), 3),
        "mgwr_ENP": round(float(mg.tr_S), 1),
        "mgwr_bandwidths": {nm: int(b) for nm, b in zip(names, msel.bw[0])},
        "mgwr_ENP_j": {nm: round(float(e), 1) for nm, e in zip(names, np.asarray(mg.ENP_j).ravel())},
        "mgwr_sig_countries": {nm: int((tv[:, i] != 0).sum()) for i, nm in enumerate(names)},
    }
    r = RES[f"gwr_{label}"]
    print(f"[(M)GWR {label}] AICc OLS={r['ols_aicc']} GWR={r['gwr_aicc']}(bw {r['gwr_bw']}) MGWR={r['mgwr_aicc']} "
          f"| MGWR adjR2={r['mgwr_adjR2']} ENP={r['mgwr_ENP']}")
    print(f"  MGWR bandwidths: {r['mgwr_bandwidths']}")
    print(f"  MGWR significant countries / 140: {r['mgwr_sig_countries']}\n")

gwr_block("btc_rate", "btc")
gwr_block("eth_rate", "eth")


# ---------------- write results ----------------
with open(os.path.join(HERE, "results_course.json"), "w") as f:
    json.dump(RES, f, indent=1)
print("wrote results_course.json")
