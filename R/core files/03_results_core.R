# Core file 3: empirical methods and numeric results.
#
# Purpose: run the main estimators used in the thesis without any LaTeX table
# formatting. Outputs are numeric CSV files only.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(DoubleML)
  library(lmtest)
  library(mlr3)
  library(mlr3learners)
  library(plm)
  library(sandwich)
})

project_dir <- normalizePath(
  Sys.getenv("THESIS_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
data_path <- file.path(project_dir, "data", "analysis", "ukhls_youth_l_to_o_clean.rds")
output_dir <- file.path(project_dir, "tables", "core_numeric")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(RUN_METHODOLOGY_SUMMARY = "0")
source(file.path(project_dir, "R", "core files", "02_methodology_core.R"))

dt <- readRDS(data_path)
data.table::setDTthreads(1L)
lgr::get_logger("mlr3")$set_threshold("warn")

estimate_row <- function(section, outcome, treatment, model, estimate, std_error, p_value, n) {
  data.frame(
    section = section,
    outcome = outcome,
    treatment = treatment,
    model = model,
    estimate = estimate,
    std_error = std_error,
    p_value = p_value,
    n = n,
    stringsAsFactors = FALSE
  )
}

extract_term <- function(model, vcov_matrix, term, n, section, outcome, treatment, label) {
  test <- lmtest::coeftest(model, vcov. = vcov_matrix)
  if (!(term %in% rownames(test))) {
    return(estimate_row(section, outcome, treatment, label, NA_real_, NA_real_, NA_real_, n))
  }
  p_col <- grep("^Pr\\(", colnames(test), value = TRUE)[1]
  estimate_row(
    section,
    outcome,
    treatment,
    label,
    test[term, "Estimate"],
    test[term, "Std. Error"],
    test[term, p_col],
    n
  )
}

fit_pooled_ols <- function(data, outcome, treatment, controls, label) {
  vars <- unique(c("pidp", outcome, treatment, controls))
  model_data <- droplevels(data[complete.cases(data[, vars]), vars, drop = FALSE])
  formula <- as.formula(paste(outcome, "~", paste(c(treatment, controls), collapse = " + ")))
  model <- lm(formula, data = model_data)
  vcov_matrix <- sandwich::vcovCL(model, cluster = model_data$pidp, type = "HC1")
  extract_term(model, vcov_matrix, treatment, nrow(model_data), "linear", outcome, treatment, label)
}

fit_panel_models <- function(data, outcome, treatment) {
  controls <- unique(c(rich_controls, "wave_f"))
  vars <- unique(c("pidp", "wave_number", outcome, treatment, controls))
  model_data <- droplevels(data[complete.cases(data[, vars]), vars, drop = FALSE])
  pdata <- pdata.frame(model_data, index = c("pidp", "wave_number"), drop.index = FALSE)
  formula <- as.formula(paste(outcome, "~", paste(c(treatment, controls), collapse = " + ")))

  fixed <- plm(formula, data = pdata, model = "within", effect = "individual")
  fixed_vcov <- vcovHC(fixed, method = "arellano", type = "HC1", cluster = "group")
  random <- plm(formula, data = pdata, model = "random", random.method = "swar")
  random_vcov <- vcovHC(random, method = "arellano", type = "HC1", cluster = "group")

  bind_rows(
    extract_term(fixed, fixed_vcov, treatment, nrow(model_data), "panel", outcome, treatment, "Individual fixed effects + wave FE"),
    extract_term(random, random_vcov, treatment, nrow(model_data), "panel", outcome, treatment, "Random effects + wave FE")
  )
}

fit_plr <- function(data, outcome, treatment = continuous_treatment) {
  design <- prepare_dml_design(data, outcome, treatment, rich_controls, include_wave = TRUE)
  dml_data <- make_clustered_doubleml_data(design)
  learners <- make_plr_learners(length(design$x_cols))

  bind_rows(lapply(seq_along(learners), function(i) {
    learner_name <- names(learners)[i]
    set.seed(1000L + i)
    fit <- DoubleMLPLR$new(
      dml_data,
      learners[[i]]$clone(),
      learners[[i]]$clone(),
      n_folds = 5
    )
    fit$fit()
    summary <- as.data.frame(fit$summary())
    estimate_row(
      "PLR-DML",
      outcome,
      treatment,
      learner_name,
      unname(summary[1, 1]),
      unname(summary[1, 2]),
      unname(summary[1, 4]),
      design$n
    )
  }))
}

fit_irm <- function(data, outcome, treatment = binary_treatment) {
  design <- prepare_dml_design(data, outcome, treatment, rich_controls, include_wave = TRUE)
  dml_data <- make_clustered_doubleml_data(design)
  learners <- make_irm_learners(length(design$x_cols))

  bind_rows(lapply(seq_along(learners), function(i) {
    learner_name <- names(learners)[i]
    set.seed(2000L + i)
    fit <- DoubleMLIRM$new(
      dml_data,
      learners[[i]]$ml_g$clone(),
      learners[[i]]$ml_m$clone(),
      score = "ATE",
      n_folds = 5
    )
    fit$fit()
    summary <- as.data.frame(fit$summary())
    estimate_row(
      "IRM-DML",
      outcome,
      treatment,
      learner_name,
      unname(summary[1, 1]),
      unname(summary[1, 2]),
      unname(summary[1, 4]),
      design$n
    )
  }))
}

linear_results <- bind_rows(lapply(names(outcomes), function(outcome) {
  bind_rows(
    fit_pooled_ols(dt, outcome, continuous_treatment, character(0), "Pooled OLS"),
    fit_pooled_ols(dt, outcome, continuous_treatment, basic_controls, "Pooled OLS + basic controls"),
    fit_pooled_ols(dt, outcome, continuous_treatment, c(rich_controls, "wave_f"), "Pooled OLS + rich controls + wave FE"),
    fit_panel_models(dt, outcome, continuous_treatment)
  )
}))

dml_results <- bind_rows(lapply(names(outcomes), function(outcome) {
  bind_rows(
    fit_plr(dt, outcome, continuous_treatment),
    fit_irm(dt, outcome, binary_treatment)
  )
}))

# Main robustness check used in the thesis: repeat the life-dissatisfaction DML
# specifications with weekend social-media use.
weekend_results <- bind_rows(
  fit_plr(dt, "life_dissatisfaction", weekend_continuous_treatment),
  fit_irm(dt, "life_dissatisfaction", weekend_binary_treatment)
)
weekend_results$section <- paste("weekend robustness", weekend_results$section)

all_results <- bind_rows(linear_results, dml_results, weekend_results)
write.csv(linear_results, file.path(output_dir, "linear_panel_results.csv"), row.names = FALSE)
write.csv(dml_results, file.path(output_dir, "dml_results.csv"), row.names = FALSE)
write.csv(weekend_results, file.path(output_dir, "weekend_robustness_results.csv"), row.names = FALSE)
write.csv(all_results, file.path(output_dir, "all_core_results.csv"), row.names = FALSE)

message("Saved numeric core results to: ", output_dir)
print(all_results)
