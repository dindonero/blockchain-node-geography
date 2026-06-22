# Reproducibility bundle: The Geography of Blockchain Infrastructure

Minimal data and code to reproduce every statistic in the paper *The geography of blockchain infrastructure: spatial clustering, wealth, and the limits of aspatial decentralization metrics*. The paper itself (`paper.pdf`) is included in this repository.

## Contents

| File | What it is |
|---|---|
| `paper.pdf` | The paper. |
| `master.csv` | The assembled country-level dataset, human-readable (one row per country). |
| `master.gpkg` | The same dataset with country boundary geometries, used to build the spatial weights. |
| `analysis.py` | One script that reproduces all reported statistics from `master.gpkg`. |
| `requirements.txt` | Python package versions. |
| `figures/` | The maps used in the paper (per-capita node rates, FDR-corrected LISA clusters, MGWR local GDP surfaces). |

## How to run

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
python analysis.py
```

The script prints every reported number and writes `results.json`. Runtime is a few minutes (the MGWR bandwidth search is the slow part).

## Dataset columns (`master.csv`)

`ne_name`, `iso2`, `iso3`, `CONTINENT`; node counts `btc_nodes`, `eth_el_nodes`, `eth_cl_nodes`, `eth_nodes` (EL+CL); covariates `gdp_pc`, `population`, `internet_pct`, `usd_per_kwh`, `reg_status`, `mining_share`; derived per-capita rates `btc_rate`, `eth_rate`; and `complete_cov` (1 for the 140 countries with complete covariates that form the analysis sample).

## Data sources and provenance

All sources are public and free. Counts are of *reachable* nodes (accepting inbound connections, therefore visible to crawlers), captured June 2026.

| Variable | Source | Retrieved | Notes |
|---|---|---|---|
| Bitcoin nodes | Luke Dashjr DNS seeder dump (`seeds.txt`), geolocated locally with the openly licensed (CDLA) ip-location-db IP-to-country databases | 2026-06-11 | 4,255 good addresses (incl. onion/unknown); 4,227 IP-geolocatable to 87 countries; 4,226 join to 86 countries in the full table (4,192 in the 140-country analysis sample). Full country coverage (long single-node tail), so absent countries are true zeros, not top-N truncation. |
| Ethereum execution-layer nodes | ethernodes.org, live browser scrape | 2026-06-11 | 5,668 nodes, 85 countries (the source's duplicate name variants for the Netherlands and Türkiye are merged). |
| Ethereum consensus-layer nodes | ChainSafe Nodewatch GraphQL API (`nodewatch.chainsafe.io/query`) | 2026-06-11 | 8,160 beacon nodes, 72 identified countries. Summing EL+CL gives 13,828 raw records; four single-node territories without a Natural Earth polygon are dropped, so `eth_nodes` totals 13,824 (13,779 in the analysis sample). |
| GDP per capita, population, internet users | World Bank Indicators API | 2020-2024 (most recent) | GDP per capita in current USD; population 2024 (the count-model offset). |
| Household electricity price | globalpetrolprices.com | 2023-2026 average | USD per kWh. |
| Cryptocurrency regulation status | Library of Congress global survey | 2021 | Coded `legal` (no documented ban), `partial_ban`, or `general_ban`. Countries absent from the survey are coded `legal`; all such cases are minor hosts. |
| Bitcoin mining shares | Cambridge CBECI mining map, latest published month | 2022-01 | Used only for the concentration-vs-clustering comparison. |
| Country boundaries | Natural Earth 50m Admin-0 | static | For centroids and the maps. |

`master.gpkg`/`master.csv` are produced by joining these layers on ISO codes (for Natural Earth, `ISO_A3_EH`/`ADM0_A3` are preferred over `ISO_A3` to avoid the -99 sentinels; Namibia's ISO2 "NA" is read as a literal string, not a missing value). Taiwan is excluded because the World Bank publishes no data for it.

## What `analysis.py` reproduces

- **ESDA.** Global Moran's I on log per-capita rates (Bitcoin 0.484, Ethereum 0.506; 9,999 permutations); FDR-corrected LISA clusters (17 / 22 High-High European countries, 25 / 39 Low-Low); Getis-Ord Gi* hot/cold spots; the Bitcoin-Ethereum rate correlation (r = 0.91); weight-matrix sensitivity (Moran's I across k = 4..10); and the Empirical Bayes-adjusted Moran's I (Bitcoin 0.24, Ethereum 0.09), the conservative bound for unstable small-population rates.
- **Concentration vs clustering.** Nakamoto coefficient, Gini, HHI, and top-5 share for Bitcoin nodes, Ethereum nodes, and Bitcoin mining; Moran's I on log shares (the matched transform) showing nodes cluster more than mining despite mining being more concentrated.
- **Negative binomial GLM** with a population offset: VIFs (max 5.2, all < 7.5), Poisson overdispersion check, the likelihood-ratio test rejecting Poisson for the negative binomial, zero-count tally (60 / 58 of 140), profile-likelihood dispersion, incidence-rate ratios with confidence intervals, deviance explained, and residual Moran's I.
- **GWR and MGWR** on standardized log rates: AICc model comparison (OLS, GWR, MGWR), per-covariate bandwidths, effective parameters, corrected-t local significance counts, local condition numbers (max < 15, none above 30), GWR Cook's distance (max 0.33), and residual Moran's I after MGWR (Bitcoin -0.06, Ethereum -0.09).

All permutation tests use seed 42. Deterministic statistics (Moran's I, coefficients, bandwidths, concentration indices) reproduce exactly; permutation p-values are robustly significant.
