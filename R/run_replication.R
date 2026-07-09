# Run the three core replication scripts from the project root.
#
# Usage:
#   Rscript R/run_replication.R

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
  "DoubleML", "data.table", "dplyr", "glmnet", "haven", "lgr", "lmtest",
  "mlr3", "mlr3learners", "nnet", "plm", "ranger", "sandwich"
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
run_script("R/core files/03_results_core.R")

message("\nCore replication completed.")
