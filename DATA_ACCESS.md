# Data access and redistribution boundary

## Why participant-level data are not included

This study is a secondary analysis of third-party longitudinal cohort data. The authors
did not create or own the underlying participant-level records and are not authorized to
redistribute them through GitHub or Zenodo. No participant-level database, CSV, RDS,
exported scores, or other row-level intermediate file is included in this repository.

The repository and its Zenodo DOI archive the analysis **software and documentation**.
They do not constitute a public deposit of the CHARLS, HRS, or ELSA microdata.

## Obtaining the same source data

Qualified researchers can request or register for data access directly from the cohort
providers, on the same terms available to the authors:

- CHARLS: <https://charls.pku.edu.cn/>
- Health and Retirement Study: <https://hrs.isr.umich.edu/>
- English Longitudinal Study of Ageing: <https://www.elsa-project.ac.uk/>
- Harmonized variables: <https://g2aging.org/>

Access requirements, data-use agreements, and available releases are controlled by the
providers and may change. Researchers should cite the specific cohort releases they use.

## What is needed to run the code

After authorized data access, construct the three local SQLite input tables described in
[`INPUT_SCHEMA.md`](INPUT_SCHEMA.md). The public scripts then export analytic files,
fit the models, and generate the reported numerical outputs. Provider raw files and the
local SQLite database must remain outside the repository.

## Reproducibility statement

The public archive documents the statistical implementation used for the manuscript.
Exact end-to-end reproduction additionally requires authorized copies of the same cohort
releases and the harmonized variable construction described by the cohort documentation.
The authors had no special access privileges.
