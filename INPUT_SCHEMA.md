# Input schema for the local analysis database

## Overview

Set the environment variable `CESD_DB` to a local SQLite database containing these
tables:

| Cohort | Table | Person identifier | Waves represented in the analysis database |
|---|---|---|---|
| CHARLS | `charls_cesd_items_long` | `ID` | 1–4 |
| ELSA | `elsa_cesd_items_long` | `idauniq` | 1–9 |
| HRS | `hrs_cesd_items_long` | `hhidpn` | 2–15 |

Each row represents a cohort observation with harmonized depressive-symptom items and
covariates. The database itself is not distributed.

## Required common fields

| Field | Meaning and coding used by the scripts |
|---|---|
| person ID | Stable cohort-specific person identifier listed above |
| `wave` | Integer harmonized wave number |
| `cohort` | Cohort label (`CHARLS`, `ELSA`, or `HRS`) |
| `ragender` | `1` = male, `2` = female; other values treated as missing |
| `agey` | Age in years |
| `cesd` | Official/harmonized CES-D total used for raw-score analyses |
| `*_d` | Item recoded in the depressive direction; higher values indicate more symptoms |

CHARLS uses ten four-level items (`0`–`3`): `depres_d`, `effort_d`, `sleepr_d`,
`whappy_d`, `flone_d`, `going_d`, `bother_d`, `mindts_d`, `fhope_d`, and `fear_d`.

ELSA and HRS use eight binary items (`0`/`1`): `depres_d`, `effort_d`, `sleepr_d`,
`whappy_d`, `flone_d`, `fsad_d`, `going_d`, and `enlife_d`.

The six common items used for cross-cohort alignment are `depres`, `effort`, `sleepr`,
`whappy`, `flone`, and `going` (using their depressive-direction versions before the
binary common-item export).

## Cohort-specific fields used by the release

- CHARLS: `raeducl`, `raeduc_c`, `mstat`, `grip_max`, and optional health/function
  variables used in sensitivity analyses.
- ELSA: `raeducl`, `raeduc_e`, `mstat`, `grip_max`, `crp`, and `fibrinogen`.
- HRS: `raeduc`, `raedyrs`, `mstat`, `cidi_mde3`, `cidi_mde5`, `cidi_symp`,
  `cidi_dep`, and `cidi_anh`.

The export script uses `any_of()` for optional external-validity fields but requires all
item, age, sex, and cohort-specific education fields used in the main analyses.

## Integrity checks performed for the archived release

- CHARLS and HRS had no duplicate person-wave keys in the local analysis database.
- ELSA contained two duplicated person-wave keys at wave 8 because the source merge
  retained multiple nurse/biomarker observations. One of these persons remained in the
  age-eligible pooled handgrip export. The published pooled handgrip model clusters
  standard errors by person, so all within-person observations—including this repeated
  within-wave measurement—are handled within the same cluster. The main ELSA wave-6
  analysis is unaffected.
- The primary main-wave samples are complete on the required CES-D items and use adults
  aged 60 years or older.

Researchers reconstructing the database should document their source releases, missing
codes, value recoding, and any rule used to collapse multiple measurements within a wave.
