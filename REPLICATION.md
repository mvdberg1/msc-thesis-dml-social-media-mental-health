# Replication Guide

This repository contains the core R replication package for the thesis
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
`fixest`, `sandwich`, `lmtest`, `glmnet`, `ranger`, `nnet`, `future`, `lgr`,
`ggplot2`, `gridExtra`, `cowplot`, `nlme`, `rpart`, and `scales`. The exact
versions used for the final analysis are recorded in `R/package_versions.csv`.
A TeX installation is only needed if the thesis PDF itself is rebuilt outside
this replication repository.

All stochastic DML folds, learners, simulations, and bootstrap procedures use
fixed seeds in the corresponding scripts. Respondent identifiers, rather than
person-wave records, define folds and bootstrap clusters.

## Workflows

The complete core workflow rebuilds the cleaned analysis data and regenerates
the thesis tables and figures produced by the R scripts:

```sh
Rscript R/run_replication.R core
```

For a faster main-text check that skips the appendix robustness scripts:

```sh
Rscript R/run_replication.R main
```

The complete workflow is slower because it re-estimates Double Machine
Learning learner comparisons, paired bootstrap tests, and robustness checks.
Random seeds and computational settings are fixed in the corresponding scripts.

## Script order

1. `R/core files/01_data_core.R` constructs the cleaned analysis data from
   licensed UKHLS youth files and writes the Chapter 3 descriptive outputs.
2. `R/core files/02_methodology_core.R` writes the Chapter 4 methodology
   figures, including the PLR causal diagram and bias-comparison figure.
3. `R/core files/03_dml_results_core.R` estimates the pooled OLS, PLR-DML, and
   IRM-DML result tables.
4. `R/core files/04_linear_panel_core.R` estimates the linear panel benchmarks,
   panel-model tests, and panel schematic.
5. `R/core files/05_prepare_core_inputs.R` prepares numeric inputs used by the
   paired comparisons and summary figures.
6. `R/core files/06_fe_dml_comparison_core.R` and
   `R/core files/07_fe_irm_comparison_core.R` run the fixed-effects versus DML
   paired comparisons. `R/run_replication.R core` repeats these scripts for all
   four outcomes.
7. `R/core files/08_robustness_design_core.R` to
   `R/core files/11_robustness_paired_comparisons_core.R` run the robustness
   checks used in the main text and appendix.
8. `R/core files/12_summary_figures_core.R` writes the Chapter 5 summary
   figures.

Generated analysis data are stored locally under `data/analysis/`; numeric
estimates under `tables/`; and figures under `figures/`. These generated
folders are not redistributed in the GitHub repository. The raw UKHLS files are
never modified. Licensed users can regenerate the local outputs from the
scripts above.
