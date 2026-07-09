# Reproduce the thesis outputs from the project root.
#
# Usage:
#   Rscript R/run_replication.R core
#   Rscript R/run_replication.R full

mode <- if (length(commandArgs(trailingOnly = TRUE))) {
  commandArgs(trailingOnly = TRUE)[1]
} else {
  "core"
}
if (!(mode %in% c("core", "full"))) {
  stop("Mode must be 'core' or 'full'.")
}

project_dir <- normalizePath(
  Sys.getenv("THESIS_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
Sys.setenv(THESIS_PROJECT_DIR = project_dir)

run_script <- function(path, env = character()) {
  message("\nRunning ", path)
  status <- system2(
    command = file.path(R.home("bin"), "Rscript"),
    args = shQuote(file.path(project_dir, path)),
    env = env
  )
  if (!identical(status, 0L)) {
    stop("Replication stopped in ", path, call. = FALSE)
  }
}

required_packages <- c(
  "DoubleML", "cowplot", "data.table", "dplyr", "fixest", "future",
  "ggplot2", "glmnet", "gridExtra", "haven", "lmtest", "mlr3",
  "mlr3learners", "nlme", "nnet", "plm", "ranger", "rpart", "sandwich",
  "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the missing R packages before replication: ",
    paste(missing_packages, collapse = ", ")
  )
}

run_script("R/thesis files/01_clean_ukhls_youth.R")
run_script("R/thesis files/02_preliminary_ols.R")
run_script("R/thesis files/03_linear_panel_benchmarks.R")
run_script("R/core files/03_prepare_core_inputs.R")

if (mode == "full") {
  outcomes <- c(
    "life_dissatisfaction",
    "loneliness",
    "schoolwork_dissatisfaction",
    "school_dissatisfaction"
  )
  for (outcome in outcomes) {
    run_script(
      "R/core files/04_fe_dml_comparison.R",
      paste0("FE_DML_OUTCOME=", outcome)
    )
    run_script(
      "R/core files/05_fe_irm_comparison.R",
      paste0("FE_IRM_OUTCOME=", outcome)
    )
  }
  run_script("R/core files/06_robustness_design_checks.R")
  run_script("R/core files/07_robustness_dml_diagnostics.R")
  run_script("R/core files/08_weekend_fe_dml_comparison.R")
  run_script("R/core files/09_robustness_paired_comparisons.R")
}

run_script("R/core files/10_results_summary_figure.R")

latexmk <- Sys.which("latexmk")
if (nzchar(latexmk)) {
  message("\nCompiling thesis.tex")
  status <- system2(
    latexmk,
    c("-pdf", "-interaction=nonstopmode", "thesis.tex"),
    stdout = "",
    stderr = ""
  )
  if (!identical(status, 0L)) {
    stop("LaTeX compilation failed.", call. = FALSE)
  }
} else {
  message("latexmk was not found; tables and figures were generated without PDF compilation.")
}

message("\nReplication completed in '", mode, "' mode.")
