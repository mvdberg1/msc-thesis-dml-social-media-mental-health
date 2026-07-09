# Replication Guide

This repository contains the R and LaTeX files used for the thesis analysis.
All commands below are run from the project root.

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

## One-command workflow

The core workflow rebuilds the cleaned data, descriptive outputs, principal
linear and DML tables, compact numeric inputs, summary figures, and PDF:

```sh
Rscript R/run_replication.R core
```

The full workflow additionally re-estimates all outcome-specific paired
FE-DML comparisons and robustness analyses:

```sh
Rscript R/run_replication.R full
```

Full replication uses 500 respondent-cluster bootstrap replications for the
paired comparisons and 100 repeated cross-fitting partitions for the split
stability analysis. It is therefore substantially slower than the core run.
The replication counts and computational settings can be changed through the
environment variables documented near the beginning of the relevant scripts.

## Script order

1. `R/thesis files/01_clean_ukhls_youth.R` constructs the analysis data,
   variable audit, descriptive tables, and descriptive figures.
2. `R/thesis files/02_preliminary_ols.R` estimates the PLR- and IRM-DML learner
   comparisons.
3. `R/thesis files/03_linear_panel_benchmarks.R` estimates the linear panel
   models and specification tests.
4. `R/core files/03_prepare_core_inputs.R` creates the compact numeric files
   consumed by later comparisons and summary figures.
5. `R/core files/04_fe_dml_comparison.R` and
   `R/core files/05_fe_irm_comparison.R` run the outcome-specific paired tests.
6. `R/core files/06_robustness_design_checks.R` through
   `R/core files/09_robustness_paired_comparisons.R` produce the robustness
   analyses.
7. `R/core files/10_results_summary_figure.R` creates the chapter summary
   figures.
8. `latexmk -pdf thesis.tex` compiles the thesis.

Generated analysis data are stored under `data/analysis/`; numeric estimates
and LaTeX tables under `tables/`; figures under `figures/`; and the final
document as `thesis.pdf`. The raw UKHLS files are never modified. In a public
GitHub or Canvas replication package, `data/analysis/` should be regenerated
locally by licensed users rather than redistributed, because it can contain
cleaned individual-level extracts from the restricted survey data.
