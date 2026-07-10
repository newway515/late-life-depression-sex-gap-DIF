# Late-life depression sex gap — DIF & robustness analysis

Analysis code for the study:

> **Is the sex gap in late-life depressive symptoms robust to measurement non-invariance? A cross-cohort differential item functioning and decomposition study in CHARLS, HRS, and ELSA.**

This repository contains the R and Python code that reproduces all results, tables, and
figures. It tests whether the higher burden of depressive symptoms (CES-D) among older
women is robust to differential item functioning (DIF) and structural adjustment across
three harmonized ageing cohorts (CHARLS, HRS, ELSA).

**No participant-level data are included in this repository** (see *Data access* below).

---

## Repository layout

```
R/        Item-response (GRM/2PL) DIF, alignment, longitudinal & sensitivity analyses
python/   Descriptive stats, effect sizes, decomposition, cluster-robust SEs, figures
```

### R scripts (run in order)

| Script | Purpose |
|--------|---------|
| `00_export_from_db.R` | Build the analytic CHARLS wave-4 sample from the SQLite item database |
| `01_charls_dif_gap.R` | Sex DIF (GRM), ESSD effect-size gate, three-estimand gap (Δ_raw / Δ_latent / Δ_adj) |
| `01b_dimensionality.R` | Single-factor vs bifactor; ECV, ωH (justifies unidimensional scoring) |
| `02_charls_longitudinal.R` | Temporal stability of the gap across CHARLS waves 1–4 |
| `03_crosscohort_alignment.R` | Cross-cohort factor-analysis alignment (common six items) |
| `04_external_validity.R` | External anchors: handgrip strength and CIDI-SF probable depression |
| `05_sensitivity_matrix.R` | Nine-cell sensitivity matrix (anchors × scale version × wave) |
| `06_dd_bootstrap.R` | Joint bootstrap CIs for the decomposition components (ΔΔ_meas, ΔΔ_struct) |
| `07_essd_threshold_sensitivity.R` | Sensitivity of Δ_latent to the ESSD threshold (0.15 / 0.20 / 0.25) |
| `08_mh_essd_concordance.R` | Concordance between the Mantel–Haenszel screen and the IRT ESSD gate |
| `09_missingness.R`, `11_missingness.R` | Item missingness, predictors, complete-case vs multiple imputation |
| `10_grip_multicov.R` | Grip-strength models with additional covariates |
| `12_export_extvalid.R` | Export corrected latent scores for the pooled external-validity analyses |
| `13_sensitivity_extra.R` | Contemporary-HRS-wave and bifactor-general-factor sensitivities |
| `export_for_python.R` | Fit the IRT models and export CSVs consumed by the Python scripts |
| `*.inp` | Mplus alignment / MNLFA input files (optional cross-check) |

### Python scripts

| Script | Purpose |
|--------|---------|
| `threecohort_descriptive_dif.py` | Table 1 descriptives, Δ_raw with bootstrap CIs, nonparametric DIF, cluster-robust external-validity SEs |
| `charls_analysis_phaseA.py` | CHARLS descriptives, sex×education collinearity, nonparametric DIF screen |
| `charls_phaseA_longitudinal.py` | Longitudinal nonparametric DIF and gap across waves |
| `gen_fig1.py`, `gen_sensitivity.py` | Figure generation (study-flow diagram; sensitivity matrix) |

---

## Requirements

**R (≥ 4.2):** `mirt`, `sirt`, `mice`, `DBI`, `RSQLite`, `dplyr`
**Python (≥ 3.9):** `numpy`, `pandas`, `matplotlib`

```r
install.packages(c("mirt","sirt","mice","DBI","RSQLite","dplyr"))
```
```bash
pip install numpy pandas matplotlib
```

## How to reproduce

1. Obtain the cohort data (see *Data access*) and build the item-level SQLite database
   expected by `00_export_from_db.R` (path set via the `CESD_DB` environment variable).
2. Run the R scripts in numeric order (`00_` → `13_`); `export_for_python.R` writes the
   CSVs used by the Python step.
3. Run the Python scripts to produce the effect sizes, decomposition, cluster-robust
   external-validity estimates, and figures.

## Data access

The analysis uses third-party data that we are not permitted to redistribute. They are
freely available to registered/authorized researchers from the original providers:

- **CHARLS** — China Health and Retirement Longitudinal Study: http://charls.pku.edu.cn
- **HRS** — Health and Retirement Study (University of Michigan): https://hrs.isr.umich.edu
- **ELSA** — English Longitudinal Study of Ageing (UK Data Service): https://www.elsa-project.ac.uk
- Harmonized versions: Gateway to Global Aging Data — https://g2aging.org

## Citation

If you use this code, please cite the associated article (see `CITATION.cff`).

## License

Released under the MIT License (see `LICENSE`).
