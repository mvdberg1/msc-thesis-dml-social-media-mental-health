# Core file 2: methodology choices.
#
# Purpose: document the analysis design used in the thesis code: outcomes,
# treatments, controls, clustered folds, and the nuisance learners used for PLR
# and IRM Double Machine Learning.

suppressPackageStartupMessages({
  library(data.table)
  library(DoubleML)
  library(mlr3)
  library(mlr3learners)
})

project_dir <- normalizePath(
  Sys.getenv("THESIS_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
data_path <- file.path(project_dir, "data", "analysis", "ukhls_youth_l_to_o_clean.rds")

outcomes <- c(
  loneliness = "Loneliness",
  life_dissatisfaction = "Life dissatisfaction",
  schoolwork_dissatisfaction = "School-work dissatisfaction",
  school_dissatisfaction = "School dissatisfaction"
)

continuous_treatment <- "social_media_weekday"
binary_treatment <- "high_social_media_weekday"
weekend_continuous_treatment <- "social_media_weekend"
weekend_binary_treatment <- "high_social_media_weekend"

basic_controls <- c("age_dv", "sex_f")
rich_controls <- c(
  "age_dv", "sex_f", "ethnicity_broad_f", "urban_f", "region_f",
  "npns_dv", "ngrp_dv", "nnssib_dv", "close_friends_log",
  "ypeatlivu", "yptvvidhrs", paste0("ypdevice", 1:6)
)

prepare_dml_design <- function(
  data,
  outcome,
  treatment,
  controls = rich_controls,
  include_wave = TRUE
) {
  x_vars <- unique(c(controls, if (include_wave) "wave_f"))
  vars <- unique(c("pidp", outcome, treatment, x_vars))
  model_data <- droplevels(data[complete.cases(data[, vars]), vars, drop = FALSE])
  x_matrix <- model.matrix(~ . - 1, data = model_data[, x_vars, drop = FALSE])

  dml_data <- data.table(
    pidp = model_data$pidp,
    y = model_data[[outcome]],
    d = model_data[[treatment]]
  )
  dml_data <- cbind(dml_data, as.data.table(x_matrix))

  list(
    data = dml_data,
    y_col = "y",
    d_col = "d",
    x_cols = setdiff(names(dml_data), c("pidp", "y", "d")),
    n = nrow(dml_data),
    clusters = length(unique(dml_data$pidp))
  )
}

make_plr_learners <- function(n_features) {
  list(
    "Elastic net" = mlr3::lrn("regr.cv_glmnet", alpha = 0.5, s = "lambda.min"),
    "Lasso" = mlr3::lrn("regr.cv_glmnet", alpha = 1, s = "lambda.min"),
    "Random forest" = mlr3::lrn(
      "regr.ranger",
      num.trees = 300L,
      mtry = max(1L, floor(sqrt(n_features))),
      min.node.size = 5L,
      num.threads = 1L
    ),
    "Neural net" = mlr3::lrn(
      "regr.nnet",
      size = 5L,
      decay = 0.01,
      maxit = 300L,
      trace = FALSE
    )
  )
}

make_irm_learners <- function(n_features) {
  list(
    "Elastic net" = list(
      ml_g = mlr3::lrn("regr.cv_glmnet", alpha = 0.5, s = "lambda.min"),
      ml_m = mlr3::lrn(
        "classif.cv_glmnet",
        alpha = 0.5,
        s = "lambda.min",
        predict_type = "prob"
      )
    ),
    "Lasso" = list(
      ml_g = mlr3::lrn("regr.cv_glmnet", alpha = 1, s = "lambda.min"),
      ml_m = mlr3::lrn(
        "classif.cv_glmnet",
        alpha = 1,
        s = "lambda.min",
        predict_type = "prob"
      )
    ),
    "Random forest" = list(
      ml_g = mlr3::lrn(
        "regr.ranger",
        num.trees = 300L,
        mtry = max(1L, floor(sqrt(n_features))),
        min.node.size = 5L,
        num.threads = 1L
      ),
      ml_m = mlr3::lrn(
        "classif.ranger",
        predict_type = "prob",
        num.trees = 300L,
        mtry = max(1L, floor(sqrt(n_features))),
        min.node.size = 5L,
        num.threads = 1L
      )
    )
  )
}

make_clustered_doubleml_data <- function(design) {
  DoubleMLClusterData$new(
    design$data,
    y_col = design$y_col,
    d_cols = design$d_col,
    x_cols = design$x_cols,
    cluster_cols = "pidp"
  )
}

if (Sys.getenv("RUN_METHODOLOGY_SUMMARY", unset = "1") == "1") {
  cat("\nCore methodology design\n")
  print(data.frame(
    item = c(
      "Continuous treatment",
      "Binary treatment",
      "Cluster variable",
      "Cross-fitting folds",
      "Main controls"
    ),
    choice = c(
      continuous_treatment,
      binary_treatment,
      "pidp",
      "5",
      paste(rich_controls, collapse = ", ")
    )
  ), row.names = FALSE)

  if (file.exists(data_path)) {
    data <- readRDS(data_path)
    design <- prepare_dml_design(data, "life_dissatisfaction", continuous_treatment)
    cat("\nLife-dissatisfaction PLR design:\n")
    print(data.frame(
      person_wave_rows = design$n,
      respondent_clusters = design$clusters,
      nuisance_features = length(design$x_cols)
    ), row.names = FALSE)
    cat("\nPLR learners:", paste(names(make_plr_learners(length(design$x_cols))), collapse = ", "), "\n")
    cat("IRM learners:", paste(names(make_irm_learners(length(design$x_cols))), collapse = ", "), "\n")
  }
}
