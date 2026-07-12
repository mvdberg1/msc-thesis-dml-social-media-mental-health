# Replication Guide

This repository contains the core R replication package for the thesis
analysis. All commands below are run from
the project root.

## Data

The analysis uses *Understanding Society: Waves 1-15, 2009-2024 and
Harmonised BHPS: Waves 1-18, 1991-2009*, 20th Edition, UK Data Service
Study Number 6614. Access is subject to the UK Data Service licence, so the
raw files cannot be redistributed with the replication code.

After obtaining the Stata release, place the youth files at:

```text
UKDA-6614-stata/stata/stata14_se/ukhls/a_youth.dta
...
UKDA-6614-stata/stata/stata14_se/ukhls/o_youth.dta
```

Alternatively, set `UKHLS_RAW_DIR` to the directory containing these files.
Set `THESIS_PROJECT_DIR` if the scripts are not launched from the project root.
The cleaning script audits variable availability across all fifteen waves and
records the resulting wave and variable documentation under `data/analysis/`.

## Software

The analysis was run with R 4.4 and the following principal packages used by
the core scripts:
`DoubleML`, `mlr3`, `mlr3learners`, `data.table`, `dplyr`, `haven`, `plm`,
`sandwich`, `lmtest`, `glmnet`, `ranger`, `nnet`, and `lgr`. The exact versions
used for the final analysis are recorded in `R/package_versions.csv`. 

All stochastic DML folds, learners, simulations, and bootstrap procedures use
fixed seeds in the corresponding scripts. Respondent identifiers, rather than
person-wave records, define folds and bootstrap clusters.

## Workflows

The core workflow rebuilds the cleaned analysis data and runs the empirical
methods:

```sh
Rscript R/run_replication.R
```

The workflow re-estimates the main linear, panel, PLR-DML, and IRM-DML models
and saves numeric CSV outputs under `tables/core_numeric/`. Random seeds and
computational settings are fixed in the corresponding scripts.

## Script order

1. `R/core files/01_data_core.R` constructs the cleaned analysis data from
   licensed UKHLS youth files.
2. `R/core files/02_methodology_core.R` records the DML setup: treatments,
   outcomes, controls, respondent clusters, cross-fitting, and nuisance
   learners.
3. `R/core files/03_results_core.R` estimates the main pooled OLS, panel,
   PLR-DML, IRM-DML, and weekend robustness specifications and writes numeric
   CSV files.

Generated analysis data are stored locally under `data/analysis/`; numeric
estimates under `tables/core_numeric/`. These generated folders are not
redistributed in the GitHub repository. The raw UKHLS files are never modified.
Licensed users can regenerate the local outputs from the scripts above.
