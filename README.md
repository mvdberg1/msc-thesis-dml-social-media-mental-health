# Double Machine Learning for Social Media and Adolescent Mental Well-Being

This repository contains a compact replication package for:

**Double Machine Learning for Estimating the Effects of Social Media Use on Mental Health**

Marloes van den Berg, MSc Econometrics, University of Amsterdam.

## Repository contents

- `R/run_replication.R`: entry point for the compact replication workflow.
- `R/thesis files/`: production scripts used for data construction and the main empirical estimates.
- `R/core files/01_data_core.R`: compact Chapter 3 data and descriptive analysis code.
- `R/core files/02_methodology_core.R`: compact Chapter 4 methodology-figure code.
- `R/core files/03_results_core.R`: compact Chapter 5 main-results code.
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

The compact replication entry point is:

```sh
Rscript R/run_replication.R core
```

To rerun the production scripts that generated the main empirical tables:

```sh
Rscript R/run_replication.R production
```

See `REPLICATION.md` for required software, directory structure, script order,
random seeds, and notes on computation time.

## Notes for submission

This repository is intended as a transparent replication package, not as a full
thesis archive. It should not contain raw UKHLS data, cleaned individual-level
UKHLS extracts, licensed documentation, temporary build files, local
machine-specific files, or thesis submission files. The `.gitignore` file is
configured to exclude those materials.
