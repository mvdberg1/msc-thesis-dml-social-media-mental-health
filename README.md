# Double Machine Learning for Social Media and Adolescent Mental Well-Being

This repository contains the replication package for:

**Double Machine Learning for Estimating the Effects of Social Media Use on Mental Health**

Marloes van den Berg, MSc Econometrics, University of Amsterdam.

## Repository contents

- `R/run_replication.R`: entry point for the ordered replication workflow.
- `R/core files/`: the core scripts used to regenerate the tables and figures
  referenced from the thesis text.
  - `01_data_core.R`: data construction and Chapter 3 descriptives.
  - `02_methodology_core.R`: Chapter 4 methodology figures.
  - `03_dml_results_core.R`: pooled OLS, PLR-DML, and IRM-DML results.
  - `04_linear_panel_core.R`: linear panel benchmarks and model tests.
  - `05_prepare_core_inputs.R`: numeric inputs for comparisons and figures.
  - `06_fe_dml_comparison_core.R` to `12_summary_figures_core.R`: paired
    comparisons, robustness checks, and Chapter 5 summary figures.
- `R/package_versions.csv`: R package versions used for the final analysis.
- `REPLICATION.md`: step-by-step data and software instructions.

The repository intentionally does not include the thesis PDF, LaTeX source,
raw UKHLS data, cleaned individual-level analysis data, or generated output
folders. Those files are either submitted separately, restricted by data
licence, or regenerated locally by the scripts.

## Data access

The raw survey data and cleaned individual-level analysis datasets are not
included in this repository.

The analysis uses *Understanding Society: Waves 1-15, 2009-2024 and
Harmonised BHPS: Waves 1-18, 1991-2009*, 20th Edition, UK Data Service Study
Number 6614. These data are restricted and subject to the UK Data Service
licence. Users who want to reproduce the full pipeline must obtain their own
licensed copy and place the files as described in `REPLICATION.md`. The
cleaning scripts then regenerate the local files under `data/analysis/`.

## Reproducibility

The complete replication entry point is:

```sh
Rscript R/run_replication.R core
```

For a faster main-text check that skips the appendix robustness scripts:

```sh
Rscript R/run_replication.R main
```

See `REPLICATION.md` for required software, directory structure, script order,
random seeds, and notes on computation time.

## Notes for submission

This repository is intended as a transparent replication package, not as a full
thesis archive. It should not contain raw UKHLS data, cleaned individual-level
UKHLS extracts, licensed documentation, temporary build files, local
machine-specific files, or thesis submission files. The `.gitignore` file is
configured to exclude those materials.
