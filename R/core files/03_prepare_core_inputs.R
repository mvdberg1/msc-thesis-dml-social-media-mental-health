# Build the compact numeric inputs used by the paired comparisons and figures.

suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- normalizePath(
  Sys.getenv("THESIS_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
tables_dir <- file.path(project_dir, "tables")
output_dir <- file.path(tables_dir, "chapter5_core")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

outcome_labels <- c(
  loneliness = "Loneliness",
  life_dissatisfaction = "Life dissatisfaction",
  schoolwork_dissatisfaction = "School-work dissatisfaction",
  school_dissatisfaction = "School dissatisfaction"
)
outcome_stubs <- c(
  loneliness = "loneliness",
  life_dissatisfaction = "life",
  schoolwork_dissatisfaction = "schoolwork",
  school_dissatisfaction = "school"
)

panel_tests <- fread(
  file.path(tables_dir, "table_panel_model_tests_numeric.csv")
)
panel_models <- rbindlist(list(
  panel_tests[, .(
    outcome,
    outcome_label,
    model = "Pooled OLS",
    estimate = pooled_estimate,
    std_error = pooled_se,
    p_value = pooled_p_value,
    n
  )],
  panel_tests[, .(
    outcome,
    outcome_label,
    model = "Random effects",
    estimate = random_estimate,
    std_error = random_se,
    p_value = random_p_value,
    n
  )],
  panel_tests[, .(
    outcome,
    outcome_label,
    model = "Fixed effects + wave FE",
    estimate = fixed_estimate,
    std_error = fixed_se,
    p_value = fixed_p_value,
    n
  )]
))
fwrite(
  panel_models,
  file.path(output_dir, "table_panel_benchmarks_numeric.csv")
)

read_framework <- function(prefix, framework) {
  rbindlist(lapply(names(outcome_stubs), function(outcome_name) {
    result <- fread(file.path(
      tables_dir,
      paste0(prefix, outcome_stubs[[outcome_name]], "_methods_numeric.csv")
    ))
    result[, outcome_label := unname(outcome_labels[[outcome_name]])]
    result[, .(
      framework,
      outcome,
      outcome_label,
      learner,
      estimate,
      std_error,
      p_value,
      n
    )]
  }))
}

fwrite(
  read_framework("table_dml_", "PLR"),
  file.path(output_dir, "table_dml_plr_learner_comparison_numeric.csv")
)
fwrite(
  read_framework("table_irm_", "IRM"),
  file.path(output_dir, "table_dml_irm_learner_comparison_numeric.csv")
)

file.copy(
  file.path(tables_dir, "table_panel_model_tests_numeric.csv"),
  file.path(output_dir, "table_panel_model_tests_numeric.csv"),
  overwrite = TRUE
)

cat("Wrote compact Chapter 4 numeric inputs to", output_dir, "\n")
