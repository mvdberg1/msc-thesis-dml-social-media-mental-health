# Reproduce the thesis replication code from the project root.
#
# Usage:
#   Rscript R/run_replication.R core
#   Rscript R/run_replication.R main

mode <- if (length(commandArgs(trailingOnly = TRUE))) {
  commandArgs(trailingOnly = TRUE)[1]
} else {
  "core"
}
if (!(mode %in% c("core", "main"))) {
  stop("Mode must be 'core' or 'main'.")
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
  "ggplot2", "glmnet", "gridExtra", "haven", "lgr", "lmtest", "mlr3",
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

run_script("R/core files/01_data_core.R")
run_script("R/core files/02_methodology_core.R")
run_script("R/core files/03_dml_results_core.R")
run_script("R/core files/04_linear_panel_core.R")
run_script("R/core files/05_prepare_core_inputs.R")

outcomes <- c(
  "life_dissatisfaction",
  "loneliness",
  "schoolwork_dissatisfaction",
  "school_dissatisfaction"
)
if (mode == "main") {
  outcomes <- "life_dissatisfaction"
}

for (outcome in outcomes) {
  run_script(
    "R/core files/06_fe_dml_comparison_core.R",
    env = paste0("FE_DML_OUTCOME=", outcome)
  )
  run_script(
    "R/core files/07_fe_irm_comparison_core.R",
    env = paste0("FE_IRM_OUTCOME=", outcome)
  )
}

if (mode == "core") {
  run_script("R/core files/08_robustness_design_core.R")
  run_script("R/core files/09_robustness_dml_diagnostics_core.R")
  run_script("R/core files/10_weekend_comparison_core.R")
  run_script("R/core files/11_robustness_paired_comparisons_core.R")
  run_script("R/core files/12_summary_figures_core.R")
}

message("\nReplication completed in '", mode, "' mode.")
