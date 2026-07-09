# Reproduce the compact thesis replication code from the project root.
#
# Usage:
#   Rscript R/run_replication.R core
#   Rscript R/run_replication.R production

mode <- if (length(commandArgs(trailingOnly = TRUE))) {
  commandArgs(trailingOnly = TRUE)[1]
} else {
  "core"
}
if (!(mode %in% c("core", "production"))) {
  stop("Mode must be 'core' or 'production'.")
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
run_script("R/core files/01_data_core.R")
run_script("R/core files/02_methodology_core.R")
run_script("R/core files/03_results_core.R")

if (mode == "production") {
  run_script("R/thesis files/02_preliminary_ols.R")
  run_script("R/thesis files/03_linear_panel_benchmarks.R")
}

message("\nReplication completed in '", mode, "' mode.")
