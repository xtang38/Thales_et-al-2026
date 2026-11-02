# Thales et al. — Measurement variance and objective-function formulation in synthetic gene circuit model parameterization

MATLAB code, experimental data, and fitted results for the manuscript.
The analysis covers three RNA-based gene circuits:

- **AND gate** and **CRISPRi** — objective-function comparison with profile
  likelihood for practical identifiability.
- **STAR** — objective-function comparison with GenSSI structural
  identifiability (original 8-parameter vs reduced models).

## Repository layout

```
ANDgate Files/     scripts, data, and fitted results for the AND gate
CRISPRi Files/     scripts, data, and fitted results for CRISPRi
STAR Files/        scripts, data, GenSSI, and results for STAR
Shared Files/      utilities used by multiple systems
```

### ANDgate Files/
- `*.m` — model, objective, fit, and profile-likelihood scripts
- `ANDgate_Data_Short2.xlsx` — experimental data (Time + replicates, 10–480 min)
- `Fitted_Results/` — GA fits for the nine scenarios (Raw/20%/40% noise ×
  No-Variance/Variance/Relative-Variance), ~30 restarts each

### CRISPRi Files/
- `*.m` — model, objective, fit, and profile-likelihood scripts
- `CRISPRdata_short.xlsx` — experimental data (Time + replicates, 0–240 min)
- `Fitted_Results/` — GA fits for the nine scenarios

### STAR Files/
- `scripts/` — GA run/plot scripts; `scripts/star_rebuild/` holds the reduced
  and minimal model builds and their GA + GenSSI runs
- `data/` — `STAR_Data.csv` (experimental), `Simulated_Data_Jan22_2025.mat`
  (100 synthetic datasets), `variance_Data.xlsx` (replicate variance)
- `GenSSI/` — GenSSI engine + `STAR_models/` (structural-identifiability
  outputs for the Single8nM, Reduced6, and Minimal3 models)
- `Results/` — summary CSVs, final plots, and fitted `*_results.mat`

### Shared Files/
- `utils_generate_noise.m`, `utils_get_mean_std.m`

## Notes
- CRISPRi: `CRISPRdata.xlsx` (5-min step) is the same measurements at finer
  resolution; the figures use the 10-min `CRISPRdata_short.xlsx` included here.
- `ANDgate_Fit_Raw.m` loads `ANDgate_data_Short.xlsx` (lowercase "d") — update
  it to `ANDgate_Data_Short2.xlsx` to match the included data.
- The GenSSI toolbox's unrelated demo models (Arabidopsis, HIV, etc.) are not
  included; only the engine and the STAR-specific outputs are here.
