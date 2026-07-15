# Late-life depression sex gap: DIF and robustness analysis

Versioned analysis code accompanying the manuscript:

> **Is the sex gap in late-life depressive symptoms robust to measurement non-invariance? A cross-cohort differential item functioning and decomposition study in CHARLS, HRS, and ELSA**

The repository implements the item-response, differential item functioning (DIF),
decomposition, longitudinal, cross-cohort, external-validity, and sensitivity analyses
reported in the manuscript.

## Scope and data boundary

This repository contains **analysis software only**. It does not contain participant-level
CHARLS, HRS, or ELSA data because those third-party data are distributed under the
respective cohort providers' access terms and may not be redistributed by the authors.
The Zenodo DOI for this repository therefore identifies the versioned software archive,
not a public copy of the cohort microdata.

Researchers with authorized access to the source cohorts can reconstruct the inputs
described in [`INPUT_SCHEMA.md`](INPUT_SCHEMA.md). The repository starts from three
harmonized person-wave tables in a local SQLite database. It does not provide a turnkey
converter from every provider-specific raw release because those files and variable
layouts are access-controlled and version-specific. See [`DATA_ACCESS.md`](DATA_ACCESS.md)
for provider links and the exact reproducibility boundary.

## Repository layout

```text
R/          IRT, DIF, alignment, longitudinal, decomposition, and sensitivity analyses
python/     Descriptives, nonparametric DIF, clustered models, and figure generation
```

### Main R scripts

| Script | Purpose |
|---|---|
| `00_export_from_db.R` | Validate and export analysis-ready files from the three SQLite input tables |
| `01_charls_dif_gap.R` | CHARLS sex DIF, ESSD gate, and three-estimand gap |
| `01b_dimensionality.R` | Single-factor and bifactor dimensionality checks |
| `02_charls_longitudinal.R` | Stability across CHARLS waves 1–4 |
| `03_crosscohort_alignment.R` | Alignment of the six common items across cohort-by-sex groups |
| `04_external_validity.R` | Initial CIDI-SF and handgrip external-validity models |
| `05_sensitivity_matrix.R` | Nine-setting robustness matrix |
| `06_dd_bootstrap.R` | Joint bootstrap for decomposition components |
| `07_essd_threshold_sensitivity.R` | ESSD threshold sensitivity |
| `08_mh_essd_concordance.R` | Mantel–Haenszel and IRT-ESSD concordance |
| `09_missingness.R`, `11_missingness.R` | Missingness diagnostics and multiple-imputation sensitivity |
| `10_grip_multicov.R` | Handgrip models with available additional covariates |
| `12_export_extvalid.R` | Export row-level scores for person-clustered models |
| `13_sensitivity_extra.R` | Contemporary-HRS-wave and bifactor sensitivities |
| `export_for_python.R` | Export intermediate results consumed by Python |

### Python scripts

| Script | Purpose |
|---|---|
| `threecohort_descriptive_dif.py` | Table 1 descriptives, bootstrap CIs, and nonparametric DIF |
| `charls_analysis_phaseA.py` | CHARLS descriptives and nonparametric DIF screen |
| `charls_phaseA_longitudinal.py` | Longitudinal nonparametric DIF and purified-score sensitivity |
| `clustered_external_validity.py` | Person-clustered HRS logistic and ELSA OLS models used in Table 6 |
| `gen_fig1.py`, `gen_sensitivity.py` | Figure generation |

## Requirements

- R 4.2 or newer: `DBI`, `RSQLite`, `dplyr`, `tidyr`, `mirt`, `sirt`,
  `boot`, `mice`, `pROC`, and `ggplot2`.
- Python 3.9 or newer. Install the Python packages with:

```bash
python -m pip install -r requirements.txt
```

## Reproduction workflow

Run commands from the repository root.

1. Obtain authorized copies of CHARLS, HRS, and ELSA/Harmonized data from the
   providers listed in [`DATA_ACCESS.md`](DATA_ACCESS.md).
2. Construct a local SQLite database containing the three tables and coding described
   in [`INPUT_SCHEMA.md`](INPUT_SCHEMA.md).
3. Set `CESD_DB` to that database. No local path is hard-coded in the release.
4. Run the analysis scripts in the documented order.

PowerShell example:

```powershell
$env:CESD_DB = "C:\path\to\cesd_analysis.db"
Rscript R/00_export_from_db.R
python python/charls_analysis_phaseA.py --out output
Rscript R/01_charls_dif_gap.R
Rscript R/01b_dimensionality.R
Rscript R/02_charls_longitudinal.R
Rscript R/03_crosscohort_alignment.R
Rscript R/04_external_validity.R
Rscript R/05_sensitivity_matrix.R
Rscript R/06_dd_bootstrap.R
Rscript R/07_essd_threshold_sensitivity.R
Rscript R/08_mh_essd_concordance.R
Rscript R/11_missingness.R
Rscript R/12_export_extvalid.R
Rscript R/13_sensitivity_extra.R
python python/clustered_external_validity.py --input-dir export --output output/clustered_external_validity.csv
```

The scripts use the fixed random seed `20260709` where stochastic computation is
performed. Some R steps are computationally intensive and require complete-case or
imputed participant-level inputs that cannot be included in this public archive.

## Data access

- [CHARLS](https://charls.pku.edu.cn/)
- [Health and Retirement Study](https://hrs.isr.umich.edu/)
- [English Longitudinal Study of Ageing](https://www.elsa-project.ac.uk/)
- [Gateway to Global Aging Data](https://g2aging.org/)

## Citation and license

Use the repository's `CITATION.cff` or Zenodo citation for the archived software
version. The code is released under the MIT License. Third-party cohort data remain
subject to their providers' terms and are not covered by this software license.
