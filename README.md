# Double Machine Learning for Estimating the Effects of Social Media Use on Mental Health

This repository contains the thesis source files, analysis code, generated
tables and figures, and replication documentation for:

**Double Machine Learning for Estimating the Effects of Social Media Use on Mental Health**

Marloes van den Berg, MSc Econometrics, University of Amsterdam.

## Repository contents

- `thesis.pdf`: final compiled thesis.
- `thesis.tex`: main LaTeX source file.
- `BibTeXTemplate.bib`: bibliography database used by the thesis.
- `R/`: R scripts for data cleaning, estimation, robustness checks, and figures.
- `tables/`: generated numeric tables and LaTeX table inputs.
- `figures/`: generated thesis figures and included image files.
- `REPLICATION.md`: step-by-step instructions for reproducing the analysis.
- `R/package_versions.csv`: package versions used for the final analysis.

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

The main replication entry point is:

```sh
Rscript R/run_replication.R core
```

The full workflow, including all paired comparisons and robustness analyses, is:

```sh
Rscript R/run_replication.R full
```

See `REPLICATION.md` for required software, directory structure, script order,
random seeds, and notes on computation time.

## Notes for submission

This repository is intended as a transparent replication package. It should not
contain raw UKHLS data, cleaned individual-level UKHLS extracts, licensed
documentation, temporary build files, or local machine-specific files. The
`.gitignore` file is configured to exclude those materials.
