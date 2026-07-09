# Replication Guide

This repository contains the compact R replication package for the thesis
analysis. It is not a full thesis archive: the submitted thesis PDF and LaTeX
source are kept outside the GitHub repository. All commands below are run from
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

The analysis was run with R 4.4 and the following principal packages:
`DoubleML`, `mlr3`, `mlr3learners`, `data.table`, `dplyr`, `haven`, `plm`,
`fixest`, `sandwich`, `lmtest`, `glmnet`, `ranger`, `nnet`, `ggplot2`, and
`cowplot`. The exact versions used for the final analysis are recorded in
`R/package_versions.csv`. A TeX installation with `latexmk` is needed to
rebuild the PDF.

All stochastic DML folds, learners, simulations, and bootstrap procedures use
fixed seeds in the corresponding scripts. Respondent identifiers, rather than
person-wave records, define folds and bootstrap clusters.

## Workflows

The compact workflow rebuilds the cleaned analysis data and runs the
chapter-oriented core scripts:

```sh
Rscript R/run_replication.R core
```

The production workflow reruns the original data-construction and main
empirical-estimation scripts:

```sh
Rscript R/run_replication.R production
```

The production workflow is slower because it re-estimates the Double Machine
Learning learner comparisons and linear panel benchmarks. Random seeds and
computational settings are fixed in the corresponding scripts.

## Script order

1. `R/thesis files/01_clean_ukhls_youth.R` constructs the cleaned analysis data
   from licensed UKHLS youth files.
2. `R/core files/01_data_core.R` reproduces the Chapter 3 descriptive data
   checks and figures in compact form.
3. `R/core files/02_methodology_core.R` reproduces the Chapter 4 methodology
   illustrations in compact form.
4. `R/core files/03_results_core.R` reproduces the Chapter 5 main empirical
   tables in compact form.
5. In the optional production workflow, `R/thesis files/02_preliminary_ols.R`
   and `R/thesis files/03_linear_panel_benchmarks.R` rerun the original DML and
   linear panel estimation scripts used for the final thesis outputs.

Generated analysis data are stored locally under `data/analysis/`; numeric
estimates under `tables/`; and figures under `figures/`. These generated
folders are not redistributed in the GitHub repository. The raw UKHLS files are
never modified. Licensed users can regenerate the local outputs from the
scripts above.
