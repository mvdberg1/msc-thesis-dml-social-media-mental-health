#  run 01_clean_ukhls_youth.R first
#  this file builds the preliminary chapter 5 regression tables

#  load the packages used for linear, fixed-effects, and dml estimation
library(dplyr)
library(sandwich)
library(lmtest)
library(fixest)
library(DoubleML)
library(mlr3)
library(mlr3learners)
library(data.table)

#  keep mlr3 output quiet in the terminal
lgr::get_logger("mlr3")$set_threshold("warn")

#  define the cleaned input data and the folder for regression tables
data_path <- "data/analysis/ukhls_youth_l_to_o_clean.rds"
tables_dir <- "tables"
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

#  load the cleaned chapter 3 dataset
dt <- readRDS(data_path)

#  map the outcome variable names to readable table labels
outcomes <- c(
  loneliness = "Loneliness",
  life_dissatisfaction = "Life dissatisfaction",
  schoolwork_dissatisfaction = "School-work dissatisfaction",
  school_dissatisfaction = "School dissatisfaction"
)

continuous_treatment <- "social_media_weekday"
binary_treatment <- "high_social_media_weekday"

#  define the basic control set used in the staged regressions
some_controls <- c("age_dv", "sex_f")

#  define the richer baseline control set used in the main specifications
many_controls <- c(
  "age_dv",
  "sex_f",
  "ethnicity_broad_f",
  "urban_f",
  "region_f",
  "npns_dv",
  "ngrp_dv",
  "nnssib_dv",
  "close_friends_log",
  "ypeatlivu",
  "yptvvidhrs",
  paste0("ypdevice", 1:6)
)

#  define the wave-specific, pooled, and within-person samples
sample_defs <- list(
  "Wave l" = list(data = dt[dt$wave == "l", ], wave_fe = FALSE, respondent_fe = FALSE),
  "Wave m" = list(data = dt[dt$wave == "m", ], wave_fe = FALSE, respondent_fe = FALSE),
  "Wave n" = list(data = dt[dt$wave == "n", ], wave_fe = FALSE, respondent_fe = FALSE),
  "Wave o" = list(data = dt[dt$wave == "o", ], wave_fe = FALSE, respondent_fe = FALSE),
  "Pooled" = list(data = dt, wave_fe = FALSE, respondent_fe = FALSE),
  "Pooled + wave FE" = list(data = dt, wave_fe = TRUE, respondent_fe = FALSE),
  "Within FE" = list(data = dt, wave_fe = FALSE, respondent_fe = TRUE),
  "Within FE + wave FE" = list(data = dt, wave_fe = TRUE, respondent_fe = TRUE)
)

#  build a standard linear formula with optional wave fixed effects
make_formula <- function(outcome, treatment_var, controls, wave_fe = FALSE) {
  rhs <- c(treatment_var, controls, if (wave_fe) "wave_f")
  base <- paste(outcome, "~", paste(rhs, collapse = " + "))
  as.formula(base)
}

#  build a fixest formula with respondent fixed effects
make_fe_formula <- function(outcome, treatment_var, controls, wave_fe = FALSE) {
  rhs <- c(treatment_var, controls)
  rhs_text <- if (length(rhs) == 0) "1" else paste(rhs, collapse = " + ")
  fe_text <- if (wave_fe) "pidp + wave_f" else "pidp"
  as.formula(paste0(outcome, " ~ ", rhs_text, " | ", fe_text))
}

#  keep only the columns and complete cases needed for one model
prepare_model_data <- function(data, outcome, treatment_var, controls, wave_fe = FALSE, keep_pidp = TRUE) {
  vars_needed <- unique(c(outcome, treatment_var, controls, if (wave_fe) "wave_f", if (keep_pidp) "pidp"))
  vars_needed <- intersect(vars_needed, names(data))

  data_model <- data[complete.cases(data[vars_needed]), vars_needed, drop = FALSE]
  droplevels(data_model)
}

#  estimate one ols model with the right standard-error rule
fit_ols <- function(data, outcome, treatment_var, controls, wave_fe = FALSE) {
  data_model <- prepare_model_data(
    data,
    outcome,
    treatment_var,
    controls,
    wave_fe = wave_fe,
    keep_pidp = TRUE
  )

  if (nrow(data_model) < 30 || length(unique(data_model[[treatment_var]])) < 2) {
    return(NULL)
  }

#  cluster by respondent when the sample contains repeat observations
  vcov_choice <- if (any(duplicated(data_model$pidp))) {
    ~pidp
  } else {
    "hetero"
  }

  formula <- make_formula(outcome, treatment_var, controls, wave_fe)
  model <- lm(formula, data = data_model)

  list(model = model, data = data_model, vcov = vcov_choice)
}

#  estimate one within-person fixed-effects model when variation allows it
fit_fe <- function(data, outcome, treatment_var, controls, wave_fe = FALSE) {
  data_model <- prepare_model_data(
    data,
    outcome,
    treatment_var,
    controls,
    wave_fe = wave_fe,
    keep_pidp = TRUE
  )

  if (nrow(data_model) < 50 || length(unique(data_model[[treatment_var]])) < 2) {
    return(NULL)
  }

  if (!any(duplicated(data_model$pidp))) {
    return(NULL)
  }

#  skip fixed-effects models when treatment has no within-person variation
  treatment_within_sd <- tapply(
    data_model[[treatment_var]],
    data_model$pidp,
    function(x) stats::sd(x, na.rm = TRUE)
  )

  if (all(is.na(treatment_within_sd)) || all(treatment_within_sd == 0, na.rm = TRUE)) {
    return(NULL)
  }

  formula <- make_fe_formula(outcome, treatment_var, controls, wave_fe = wave_fe)
  model <- tryCatch(
    suppressWarnings(
      fixest::feols(
        formula,
        data = data_model,
        cluster = ~pidp,
        notes = FALSE
      )
    ),
    error = function(e) NULL
  )

  if (is.null(model)) {
    return(NULL)
  }

  list(model = model, data = data_model)
}

#  define the random-forest learner used in plr nuisance estimation
make_dml_learner <- function(n_features) {
  mlr3::lrn(
    "regr.ranger",
    num.trees = 100L,
    mtry = max(1L, min(20L, as.integer(n_features))),
    min.node.size = 2L,
    max.depth = 5L,
    num.threads = 1L
  )
}

#  estimate one dml plr model for the continuous treatment
fit_dml_plr <- function(data, outcome, treatment_var, controls, wave_fe = FALSE) {
  x_cols <- unique(c(controls, if (wave_fe) "wave_f"))
  data_model <- prepare_model_data(
    data,
    outcome,
    treatment_var,
    x_cols,
    wave_fe = FALSE,
    keep_pidp = TRUE
  )

  if (nrow(data_model) < 100 || length(unique(data_model[[treatment_var]])) < 2 || length(x_cols) == 0) {
    return(NULL)
  }

#  switch to clustered dml objects when respondents appear multiple times
  use_cluster <- any(duplicated(data_model$pidp))
  learner <- make_dml_learner(length(x_cols))

  dml_df <- as.data.table(data_model)
  dml_data <- if (use_cluster) {
    DoubleMLClusterData$new(
      dml_df,
      y_col = outcome,
      d_cols = treatment_var,
      x_cols = x_cols,
      cluster_cols = "pidp"
    )
  } else {
    DoubleMLData$new(
      dml_df[, c(outcome, treatment_var, x_cols), with = FALSE],
      y_col = outcome,
      d_cols = treatment_var,
      x_cols = x_cols
    )
  }

  set.seed(1111)
  fit <- tryCatch(
    {
      obj <- DoubleMLPLR$new(
        dml_data,
        learner$clone(),
        learner$clone(),
        n_folds = 5
      )
      obj$fit()
      obj
    },
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(NULL)
  }

  list(model = fit, data = data_model, clustered = use_cluster)
}

#  define the outcome and propensity learners used in irm
make_irm_learners <- function(n_features) {
  list(
    ml_g = mlr3::lrn(
      "regr.ranger",
      num.trees = 100L,
      mtry = max(1L, min(10L, as.integer(n_features))),
      min.node.size = 2L,
      max.depth = 5L,
      num.threads = 1L
    ),
    ml_m = mlr3::lrn(
      "classif.ranger",
      predict_type = "prob",
      num.trees = 100L,
      mtry = max(1L, min(10L, as.integer(n_features))),
      min.node.size = 2L,
      max.depth = 5L,
      num.threads = 1L
    )
  )
}

#  estimate one dml irm model for the binary high-use treatment
fit_dml_irm <- function(data, outcome, treatment_var, controls, wave_fe = FALSE) {
  x_cols <- unique(c(controls, if (wave_fe) "wave_f"))
  data_model <- prepare_model_data(
    data,
    outcome,
    treatment_var,
    x_cols,
    wave_fe = FALSE,
    keep_pidp = TRUE
  )

  treatment_values <- sort(unique(data_model[[treatment_var]]))
  if (
    nrow(data_model) < 100 ||
    length(x_cols) == 0 ||
    length(treatment_values) != 2 ||
    !all(treatment_values %in% c(0, 1))
  ) {
    return(NULL)
  }

  treatment_counts <- table(data_model[[treatment_var]])
  if (any(treatment_counts < 30)) {
    return(NULL)
  }

#  switch to clustered dml objects when respondents appear multiple times
  use_cluster <- any(duplicated(data_model$pidp))
  learners <- make_irm_learners(length(x_cols))
  dml_df <- as.data.table(data_model)

  dml_data <- if (use_cluster) {
    DoubleMLClusterData$new(
      dml_df,
      y_col = outcome,
      d_cols = treatment_var,
      x_cols = x_cols,
      cluster_cols = "pidp"
    )
  } else {
    DoubleMLData$new(
      dml_df[, c(outcome, treatment_var, x_cols), with = FALSE],
      y_col = outcome,
      d_cols = treatment_var,
      x_cols = x_cols
    )
  }

  set.seed(3333)
  fit <- tryCatch(
    {
      obj <- DoubleMLIRM$new(
        dml_data,
        learners$ml_g$clone(),
        learners$ml_m$clone(),
        score = "ATE",
        n_folds = 5
      )
      obj$fit()
      obj
    },
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(NULL)
  }

  list(model = fit, data = data_model, clustered = use_cluster)
}

#  extract one treatment coefficient and its standard error from ols
extract_ols_treatment <- function(model, treatment_var) {
  if (is.null(model)) {
    return(list(cell = "", estimate = NA_real_, se = NA_real_, p = NA_real_, n = NA_integer_))
  }

  vcov_matrix <- suppressWarnings(
    if (identical(model$vcov, "hetero")) {
      vcovHC(model$model, type = "HC1")
    } else {
      vcovCL(model$model, cluster = model$data$pidp, type = "HC1")
    }
  )

  ct <- coeftest(model$model, vcov. = vcov_matrix)
  if (!(treatment_var %in% rownames(ct))) {
    return(list(cell = "", estimate = NA_real_, se = NA_real_, p = NA_real_, n = nobs(model$model)))
  }

  estimate <- unname(ct[treatment_var, "Estimate"])
  se <- unname(ct[treatment_var, "Std. Error"])
  p_value <- unname(ct[treatment_var, "Pr(>|t|)"])
  stars <- ifelse(
    p_value < 0.01,
    "***",
    ifelse(p_value < 0.05, "**", ifelse(p_value < 0.1, "*", ""))
  )

#  add significance stars in the same style used in the thesis tables
  list(
    cell = sprintf("%.3f%s (%.3f)", estimate, stars, se),
    estimate = estimate,
    se = se,
    p = p_value,
    n = nobs(model$model)
  )
}

#  extract one treatment coefficient and its standard error from fixest
extract_fe_treatment <- function(model, treatment_var) {
  if (is.null(model)) {
    return(list(cell = "", estimate = NA_real_, se = NA_real_, p = NA_real_, n = NA_integer_))
  }

  ct <- as.data.frame(fixest::coeftable(model$model))
  if (!(treatment_var %in% rownames(ct))) {
    return(list(cell = "", estimate = NA_real_, se = NA_real_, p = NA_real_, n = nobs(model$model)))
  }

  estimate <- unname(ct[treatment_var, "Estimate"])
  se <- unname(ct[treatment_var, "Std. Error"])
  p_value <- unname(ct[treatment_var, "Pr(>|t|)"])
  stars <- ifelse(
    p_value < 0.01,
    "***",
    ifelse(p_value < 0.05, "**", ifelse(p_value < 0.1, "*", ""))
  )

#  add significance stars in the same style used in the thesis tables
  list(
    cell = sprintf("%.3f%s (%.3f)", estimate, stars, se),
    estimate = estimate,
    se = se,
    p = p_value,
    n = nobs(model$model)
  )
}

#  extract the treatment estimate and standard error from a dml object
extract_dml_treatment <- function(model) {
  if (is.null(model)) {
    return(list(cell = "", estimate = NA_real_, se = NA_real_, p = NA_real_, n = NA_integer_))
  }

  summary_df <- as.data.frame(model$model$summary())
  estimate <- unname(summary_df[1, 1])
  se <- unname(summary_df[1, 2])
  p_value <- unname(summary_df[1, 4])
  stars <- ifelse(
    p_value < 0.01,
    "***",
    ifelse(p_value < 0.05, "**", ifelse(p_value < 0.1, "*", ""))
  )

#  add significance stars in the same style used in the thesis tables
  list(
    cell = sprintf("%.3f%s (%.3f)", estimate, stars, se),
    estimate = estimate,
    se = se,
    p = p_value,
    n = nrow(model$data)
  )
}

#  escape latex-sensitive characters before writing tables
escape_latex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x
}

#  write one single-outcome latex table
write_latex_table <- function(tab, file, caption, label, note_text) {
  columns <- names(tab)
  align <- paste0("l", paste(rep("c", length(columns) - 1), collapse = ""))
  n_cols <- ncol(tab)
  grouped_header <- identical(
    columns,
    c(
      "Specification",
      "Wave l",
      "Wave m",
      "Wave n",
      "Wave o",
      "Pooled",
      "Pooled + wave FE",
      "Within FE",
      "Within FE + wave FE"
    )
  )

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    paste0("\\caption{", escape_latex(caption), "}"),
    paste0("\\label{", label, "}"),
    "\\scriptsize",
    "\\resizebox{\\textwidth}{!}{%",
    paste0("\\begin{tabular}{", align, "}"),
    "\\hline"
  )

  if (grouped_header) {
    lines <- c(
      lines,
      " & \\multicolumn{4}{c}{Wave-specific OLS} & \\multicolumn{2}{c}{Pooled OLS} & \\multicolumn{2}{c}{Within-respondent FE} \\\\",
      "\\cline{2-5} \\cline{6-7} \\cline{8-9}",
      "Specification & l & m & n & o & No wave FE & + wave FE & No wave FE & + wave FE \\\\",
      "\\hline"
    )
  } else {
    lines <- c(
      lines,
      paste(escape_latex(columns), collapse = " & "),
      "\\\\",
      "\\hline"
    )
  }

#  write the body rows one by one into the latex table
  for (i in seq_len(nrow(tab))) {
    row <- as.character(tab[i, ])
    lines <- c(lines, paste(escape_latex(row), collapse = " & "), "\\\\")
  }

  lines <- c(
    lines,
    "\\hline",
    "\\end{tabular}",
    "}",
    "\\begin{minipage}{0.95\\textwidth}",
    paste0("\\footnotesize Notes: ", note_text),
    "\\end{minipage}",
    "\\end{table}"
  )

  writeLines(lines, file)
}

#  write the combined multi-panel latex table across outcomes
write_latex_panel_table <- function(panel_tables, file, caption, label, note_text) {
  sample_columns <- names(panel_tables[[1]])
  n_cols <- length(sample_columns)
  align <- paste0("l", paste(rep("c", n_cols - 1), collapse = ""))
  grouped_header <- identical(
    sample_columns,
    c(
      "Specification",
      "Wave l",
      "Wave m",
      "Wave n",
      "Wave o",
      "Pooled",
      "Pooled + wave FE",
      "Within FE",
      "Within FE + wave FE"
    )
  )

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    paste0("\\caption{", escape_latex(caption), "}"),
    paste0("\\label{", label, "}"),
    "\\scriptsize",
    "\\resizebox{\\textwidth}{!}{%",
    paste0("\\begin{tabular}{", align, "}"),
    "\\hline"
  )

  if (grouped_header) {
    lines <- c(
      lines,
      " & \\multicolumn{4}{c}{Wave-specific OLS} & \\multicolumn{2}{c}{Pooled OLS} & \\multicolumn{2}{c}{Within-respondent FE} \\\\",
      "\\cline{2-5} \\cline{6-7} \\cline{8-9}",
      "Specification & l & m & n & o & No wave FE & + wave FE & No wave FE & + wave FE \\\\",
      "\\hline"
    )
  } else {
    lines <- c(
      lines,
      paste(escape_latex(sample_columns), collapse = " & "),
      "\\\\",
      "\\hline"
    )
  }

#  add one panel per outcome to the combined latex table
  panel_names <- names(panel_tables)
  panel_labels <- LETTERS[seq_along(panel_tables)]

  for (i in seq_along(panel_tables)) {
    panel_table <- panel_tables[[i]]
    panel_title <- paste0("Panel ", panel_labels[i], ": ", panel_names[i])
    lines <- c(
      lines,
      paste0("\\multicolumn{", n_cols, "}{l}{\\textit{", escape_latex(panel_title), "}} \\\\"),
      "\\hline"
    )

    for (j in seq_len(nrow(panel_table))) {
      row <- as.character(panel_table[j, ])
      lines <- c(lines, paste(escape_latex(row), collapse = " & "), "\\\\")
    }

    lines <- c(lines, "\\hline")
  }

  lines <- c(
    lines,
    "\\end{tabular}",
    "}",
    "\\begin{minipage}{0.95\\textwidth}",
    paste0("\\footnotesize Notes: ", note_text),
    "\\end{minipage}",
    "\\end{table}"
  )

  writeLines(lines, file)
}

#  prepare a numeric dml design matrix so all learner classes use the same sample
prepare_dml_numeric_sample <- function(data, outcome, treatment_var, controls, include_wave = TRUE) {
  x_cols <- unique(c(controls, if (include_wave) "wave_f"))
  vars_needed <- unique(c(outcome, treatment_var, x_cols, "pidp"))
  vars_needed <- intersect(vars_needed, names(data))

  data_model <- droplevels(data[complete.cases(data[vars_needed]), vars_needed, drop = FALSE])
  x_matrix <- model.matrix(~ . - 1, data = data_model[, x_cols, drop = FALSE])

  dml_df <- data.table(
    outcome_value = data_model[[outcome]],
    treatment_value = data_model[[treatment_var]],
    pidp = data_model$pidp
  )
  names(dml_df)[1:2] <- c(outcome, treatment_var)
  dml_df <- cbind(dml_df, as.data.table(x_matrix))

  list(
    data = dml_df,
    x_cols = setdiff(names(dml_df), c(outcome, treatment_var, "pidp"))
  )
}

#  define the learner set used in the compact plr-dml comparison table
make_plr_comparison_learners <- function(n_features) {
  list(
    "Elastic net" = mlr3::lrn(
      "regr.cv_glmnet",
      s = "lambda.min",
      alpha = 0.5
    ),
    "Lasso" = mlr3::lrn(
      "regr.cv_glmnet",
      s = "lambda.min",
      alpha = 1
    ),
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

#  define the learner set used in the compact irm-dml comparison table
make_irm_comparison_learners <- function(n_features) {
  list(
    "Elastic net" = list(
      ml_g = mlr3::lrn(
        "regr.cv_glmnet",
        s = "lambda.min",
        alpha = 0.5
      ),
      ml_m = mlr3::lrn(
        "classif.cv_glmnet",
        s = "lambda.min",
        alpha = 0.5,
        predict_type = "prob"
      )
    ),
    "Lasso" = list(
      ml_g = mlr3::lrn(
        "regr.cv_glmnet",
        s = "lambda.min",
        alpha = 1
      ),
      ml_m = mlr3::lrn(
        "classif.cv_glmnet",
        s = "lambda.min",
        alpha = 1,
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

#  compact significance helper for the dml learner-comparison table
star_string <- function(p_value) {
  if (is.na(p_value)) {
    return("")
  }
  if (p_value < 0.01) {
    return("***")
  }
  if (p_value < 0.05) {
    return("**")
  }
  if (p_value < 0.10) {
    return("*")
  }
  ""
}

#  estimate one pooled plr-dml coefficient under several nuisance learners
build_plr_method_comparison <- function(
  outcome = "life_dissatisfaction",
  treatment_var = continuous_treatment,
  controls = many_controls,
  include_wave = TRUE
) {
  prepared <- prepare_dml_numeric_sample(
    data = dt,
    outcome = outcome,
    treatment_var = treatment_var,
    controls = controls,
    include_wave = include_wave
  )

  dml_data <- DoubleMLClusterData$new(
    prepared$data,
    y_col = outcome,
    d_cols = treatment_var,
    x_cols = prepared$x_cols,
    cluster_cols = "pidp"
  )

  learners <- make_plr_comparison_learners(length(prepared$x_cols))
  numeric_rows <- list()

  for (i in seq_along(learners)) {
    learner_name <- names(learners)[i]
    set.seed(100 + i)

    fit <- tryCatch(
      {
        obj <- DoubleMLPLR$new(
          dml_data,
          learners[[i]]$clone(),
          learners[[i]]$clone(),
          n_folds = 5
        )
        obj$fit()
        obj
      },
      error = function(e) NULL
    )

    if (is.null(fit)) {
      numeric_rows[[length(numeric_rows) + 1]] <- data.frame(
        outcome = outcome,
        learner = learner_name,
        estimate = NA_real_,
        std_error = NA_real_,
        p_value = NA_real_,
        n = nrow(prepared$data),
        stringsAsFactors = FALSE
      )
      next
    }

    summary_df <- as.data.frame(fit$summary())
    numeric_rows[[length(numeric_rows) + 1]] <- data.frame(
      framework = "PLR",
      outcome = outcome,
      treatment = treatment_var,
      learner = learner_name,
      estimate = unname(summary_df[1, 1]),
      std_error = unname(summary_df[1, 2]),
      p_value = unname(summary_df[1, 4]),
      n = nrow(prepared$data),
      stringsAsFactors = FALSE
    )
  }

  bind_rows(numeric_rows)
}

#  estimate one pooled irm-dml coefficient under several nuisance learners
build_irm_method_comparison <- function(
  outcome = "life_dissatisfaction",
  treatment_var = binary_treatment,
  controls = many_controls,
  include_wave = TRUE
) {
  prepared <- prepare_dml_numeric_sample(
    data = dt,
    outcome = outcome,
    treatment_var = treatment_var,
    controls = controls,
    include_wave = include_wave
  )

  treatment_values <- sort(unique(prepared$data[[treatment_var]]))
  if (
    length(prepared$x_cols) == 0 ||
    length(treatment_values) != 2 ||
    !all(treatment_values %in% c(0, 1))
  ) {
    return(data.frame(
      framework = "IRM",
      outcome = outcome,
      treatment = treatment_var,
      learner = names(make_irm_comparison_learners(1)),
      estimate = NA_real_,
      std_error = NA_real_,
      p_value = NA_real_,
      n = nrow(prepared$data),
      stringsAsFactors = FALSE
    ))
  }

  treatment_counts <- table(prepared$data[[treatment_var]])
  if (any(treatment_counts < 30)) {
    return(data.frame(
      framework = "IRM",
      outcome = outcome,
      treatment = treatment_var,
      learner = names(make_irm_comparison_learners(length(prepared$x_cols))),
      estimate = NA_real_,
      std_error = NA_real_,
      p_value = NA_real_,
      n = nrow(prepared$data),
      stringsAsFactors = FALSE
    ))
  }

  dml_data <- DoubleMLClusterData$new(
    prepared$data,
    y_col = outcome,
    d_cols = treatment_var,
    x_cols = prepared$x_cols,
    cluster_cols = "pidp"
  )

  learners <- make_irm_comparison_learners(length(prepared$x_cols))
  numeric_rows <- list()

  for (i in seq_along(learners)) {
    learner_name <- names(learners)[i]
    set.seed(500 + i)

    fit <- tryCatch(
      {
        obj <- DoubleMLIRM$new(
          dml_data,
          learners[[i]]$ml_g$clone(),
          learners[[i]]$ml_m$clone(),
          score = "ATE",
          n_folds = 5
        )
        obj$fit()
        obj
      },
      error = function(e) NULL
    )

    if (is.null(fit)) {
      numeric_rows[[length(numeric_rows) + 1]] <- data.frame(
        framework = "IRM",
        outcome = outcome,
        treatment = treatment_var,
        learner = learner_name,
        estimate = NA_real_,
        std_error = NA_real_,
        p_value = NA_real_,
        n = nrow(prepared$data),
        stringsAsFactors = FALSE
      )
      next
    }

    summary_df <- as.data.frame(fit$summary())
    numeric_rows[[length(numeric_rows) + 1]] <- data.frame(
      framework = "IRM",
      outcome = outcome,
      treatment = treatment_var,
      learner = learner_name,
      estimate = unname(summary_df[1, 1]),
      std_error = unname(summary_df[1, 2]),
      p_value = unname(summary_df[1, 4]),
      n = nrow(prepared$data),
      stringsAsFactors = FALSE
    )
  }

  bind_rows(numeric_rows)
}

#  write the compact learner-comparison table in the style of dmlvsols
write_dml_method_table <- function(results, file_name, caption, label, variable_label, note_text) {
  learners <- results$learner
  estimate_cells <- vapply(seq_len(nrow(results)), function(i) {
    if (any(is.na(results[i, c("estimate", "p_value")]))) {
      return("")
    }
    sprintf("%.3f%s", results$estimate[i], star_string(results$p_value[i]))
  }, character(1))

  se_cells <- vapply(seq_len(nrow(results)), function(i) {
    if (is.na(results$std_error[i])) {
      return("")
    }
    sprintf("(%.3f)", results$std_error[i])
  }, character(1))

  n_cells <- format(results$n, big.mark = ",")

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    paste0("\\begin{tabular}{l", paste(rep("c", length(learners)), collapse = ""), "}"),
    "\\hline",
    paste(c("Variable", learners), collapse = " & "),
    "\\\\",
    "\\hline",
    paste(c(variable_label, estimate_cells), collapse = " & "),
    "\\\\",
    paste(c("", se_cells), collapse = " & "),
    "\\\\",
    paste(c("N", n_cells), collapse = " & "),
    "\\\\",
    "\\hline",
    "\\end{tabular}",
    "\\begin{minipage}{0.95\\textwidth}",
    paste0("\\footnotesize Notes: ", note_text),
    "\\end{minipage}",
    "\\end{table}"
  )

  writeLines(lines, file.path(tables_dir, file_name))
  write.csv(
    results,
    file.path(tables_dir, sub("\\.tex$", "_numeric.csv", file_name)),
    row.names = FALSE
  )
}

write_plr_method_tables <- function(outcome_names) {
  plr_note <- paste(
    "All columns report pooled double machine learning estimates from the partially linear regression specification for the coefficient on weekday social-media use.",
    "The same common outcome-specific sample with rich controls is used throughout, and wave fixed effects are included among the controls.",
    "Elastic net uses alpha = 0.5, while Lasso uses alpha = 1.",
    "Standard errors come from respondent-clustered DoubleML inference with 5-fold cross-fitting.",
    "Significance: *** \\(p < 0.01\\), ** \\(p < 0.05\\), * \\(p < 0.10\\)."
  )

  for (outcome_name in outcome_names) {
    file_stub <- switch(
      outcome_name,
      loneliness = "loneliness",
      life_dissatisfaction = "life",
      schoolwork_dissatisfaction = "schoolwork",
      school_dissatisfaction = "school",
      gsub("[^a-z0-9_]+", "_", outcome_name)
    )
    results <- build_plr_method_comparison(
      outcome = outcome_name,
      treatment_var = continuous_treatment,
      controls = many_controls,
      include_wave = TRUE
    )

    write_dml_method_table(
      results = results,
      file_name = paste0("table_dml_", file_stub, "_methods.tex"),
      caption = paste0("Partially linear double machine learning estimates for ", tolower(outcomes[[outcome_name]]), " across nuisance learners"),
      label = paste0("tab:dml_", file_stub, "_methods"),
      variable_label = "Weekday social-media use",
      note_text = plr_note
    )
  }
}

write_irm_method_tables <- function(outcome_names) {
  irm_note <- paste(
    "All columns report pooled double machine learning average treatment effect estimates from the interactive regression model for the binary indicator of at least four hours of weekday social-media use.",
    "The same common outcome-specific sample with rich controls is used throughout, and wave fixed effects are included among the controls.",
    "Elastic net uses alpha = 0.5, while Lasso uses alpha = 1 for both the outcome and propensity nuisances.",
    "The interactive-regression comparison is restricted to Elastic net, Lasso, and Random forest because these were the learners that remained numerically stable in the binary-treatment setting.",
    "Standard errors come from respondent-clustered DoubleML inference with 5-fold cross-fitting.",
    "Significance: *** \\(p < 0.01\\), ** \\(p < 0.05\\), * \\(p < 0.10\\)."
  )

  for (outcome_name in outcome_names) {
    file_stub <- switch(
      outcome_name,
      loneliness = "loneliness",
      life_dissatisfaction = "life",
      schoolwork_dissatisfaction = "schoolwork",
      school_dissatisfaction = "school",
      gsub("[^a-z0-9_]+", "_", outcome_name)
    )
    results <- build_irm_method_comparison(
      outcome = outcome_name,
      treatment_var = binary_treatment,
      controls = many_controls,
      include_wave = TRUE
    )

    write_dml_method_table(
      results = results,
      file_name = paste0("table_irm_", file_stub, "_methods.tex"),
      caption = paste0("Interactive regression model double machine learning estimates for ", tolower(outcomes[[outcome_name]]), " across nuisance learners"),
      label = paste0("tab:irm_", file_stub, "_methods"),
      variable_label = "High social-media use",
      note_text = irm_note
    )
  }
}

#  run the full set of models and export tables for one treatment definition
build_results_suite <- function(
  treatment_var,
  specs,
  table_prefix,
  combined_caption,
  per_outcome_caption_prefix,
  per_outcome_label_prefix,
  combined_label,
  note_text
) {
  all_results <- list()
  all_display_tables <- list()

#  loop over outcomes and estimate each specification in each sample
  for (outcome in names(outcomes)) {
    display <- data.frame(
      Specification = vapply(specs, function(x) x$label, character(1)),
      stringsAsFactors = FALSE
    )
    numeric_rows <- list()

#  fill one display column for each sample definition
    for (sample_name in names(sample_defs)) {
      cells <- character(length(specs))

      for (i in seq_along(specs)) {
        spec <- specs[[i]]
        sample_def <- sample_defs[[sample_name]]
        model <- NULL
        extracted <- list(cell = "", estimate = NA_real_, se = NA_real_, p = NA_real_, n = NA_integer_)
        actual_estimator <- spec$estimator

        if (identical(spec$estimator, "linear")) {
          if (isTRUE(sample_def$respondent_fe)) {
            model <- fit_fe(
              data = sample_def$data,
              outcome = outcome,
              treatment_var = treatment_var,
              controls = spec$controls,
              wave_fe = sample_def$wave_fe
            )
            extracted <- extract_fe_treatment(model, treatment_var)
            actual_estimator <- "fe"
          } else {
            model <- fit_ols(
              data = sample_def$data,
              outcome = outcome,
              treatment_var = treatment_var,
              controls = spec$controls,
              wave_fe = sample_def$wave_fe
            )
            extracted <- extract_ols_treatment(model, treatment_var)
            actual_estimator <- "ols"
          }
        } else if (identical(spec$estimator, "dml_plr")) {
          if (!isTRUE(sample_def$respondent_fe)) {
            model <- fit_dml_plr(
              data = sample_def$data,
              outcome = outcome,
              treatment_var = treatment_var,
              controls = spec$controls,
              wave_fe = sample_def$wave_fe
            )
            extracted <- extract_dml_treatment(model)
          }
        } else if (identical(spec$estimator, "dml_irm")) {
          if (!isTRUE(sample_def$respondent_fe)) {
            model <- fit_dml_irm(
              data = sample_def$data,
              outcome = outcome,
              treatment_var = treatment_var,
              controls = spec$controls,
              wave_fe = sample_def$wave_fe
            )
            extracted <- extract_dml_treatment(model)
          }
        }

#  store both the printable cell and the numeric result underneath
        cells[i] <- extracted$cell

        numeric_rows[[length(numeric_rows) + 1]] <- data.frame(
          outcome = outcome,
          outcome_label = outcomes[[outcome]],
          treatment = treatment_var,
          specification = spec$label,
          estimator = actual_estimator,
          sample = sample_name,
          wave_fe = sample_def$wave_fe,
          respondent_fe = sample_def$respondent_fe,
          estimate = extracted$estimate,
          std_error = extracted$se,
          p_value = extracted$p,
          n = extracted$n,
          stringsAsFactors = FALSE
        )
      }

      display[[sample_name]] <- cells
    }

#  save the per-outcome csv, numeric csv, and latex table
    safe_outcome <- gsub("[^a-z0-9_]+", "_", outcome)
    csv_path <- file.path(tables_dir, paste0(table_prefix, "_", safe_outcome, ".csv"))
    tex_path <- file.path(tables_dir, paste0(table_prefix, "_", safe_outcome, ".tex"))
    numeric_path <- file.path(tables_dir, paste0(table_prefix, "_", safe_outcome, "_numeric.csv"))

    numeric_out <- bind_rows(numeric_rows)
    write.csv(display, csv_path, row.names = FALSE)
    write.csv(numeric_out, numeric_path, row.names = FALSE)
    write_latex_table(
      display,
      tex_path,
      caption = paste(per_outcome_caption_prefix, outcomes[[outcome]]),
      label = paste0(per_outcome_label_prefix, safe_outcome),
      note_text = note_text
    )

    all_results[[outcome]] <- numeric_out
    all_display_tables[[outcomes[[outcome]]]] <- display
    print(display)
    message("Saved: ", csv_path)
    message("Saved: ", tex_path)
  }

#  save the stacked numeric results and the combined panel table
  all_results_df <- bind_rows(all_results)
  write.csv(
    all_results_df,
    file.path(tables_dir, paste0(table_prefix, "_all_outcomes_numeric.csv")),
    row.names = FALSE
  )

  combined_csv <- bind_rows(lapply(names(all_display_tables), function(outcome_label) {
    tab <- all_display_tables[[outcome_label]]
    cbind(Outcome = outcome_label, tab, stringsAsFactors = FALSE)
  }))

  write.csv(
    combined_csv,
    file.path(tables_dir, paste0(table_prefix, "_combined.csv")),
    row.names = FALSE
  )

  write_latex_panel_table(
    all_display_tables,
    file.path(tables_dir, paste0(table_prefix, "_combined.tex")),
    caption = combined_caption,
    label = combined_label,
    note_text = note_text
  )
}

#  define the staged specifications for the continuous-treatment benchmark
plr_specs <- list(
  list(label = "Bivariate", estimator = "linear", controls = character(0)),
  list(label = "+ Basic controls", estimator = "linear", controls = some_controls),
  list(label = "+ Rich controls", estimator = "linear", controls = many_controls),
  list(label = "DML PLR (+ rich controls)", estimator = "dml_plr", controls = many_controls)
)

#  define the staged specifications for the binary-treatment benchmark
irm_specs <- list(
  list(label = "Bivariate", estimator = "linear", controls = character(0)),
  list(label = "+ Basic controls", estimator = "linear", controls = some_controls),
  list(label = "+ Rich controls", estimator = "linear", controls = many_controls),
  list(label = "DML IRM ATE (+ rich controls)", estimator = "dml_irm", controls = many_controls)
)

#  write the table notes shown below the plr results
plr_note <- paste(
  "cells report the coefficient on weekday social-media use with standard errors in parentheses.",
  "The first three rows are linear specifications; wave-specific and pooled columns use OLS, while the last two columns use within-respondent fixed effects.",
  "Wave-specific OLS columns use HC1 robust standard errors; pooled OLS and within-FE columns cluster standard errors at the respondent level.",
  "The DML row uses DoubleMLPLR with partialling out, 5-fold cross-fitting, and \\\\texttt{regr.ranger} nuisance learners.",
  "Pooled DML columns use respondent clustering via \\\\texttt{DoubleMLClusterData}; DML with respondent fixed effects is not yet reported and those cells are left blank.",
  "Significance: *** \\(p < 0.01\\), ** \\(p < 0.05\\), * \\(p < 0.10\\)."
)

#  write the table notes shown below the irm results
irm_note <- paste(
  "cells report the coefficient on the binary indicator for at least four hours of weekday social-media use with standard errors in parentheses.",
  "The first three rows are linear specifications; wave-specific and pooled columns use OLS, while the last two columns use within-respondent fixed effects.",
  "Wave-specific OLS columns use HC1 robust standard errors; pooled OLS and within-FE columns cluster standard errors at the respondent level.",
  "The DML row uses DoubleMLIRM with the ATE score, 5-fold cross-fitting, \\\\texttt{regr.ranger} for the outcome nuisance, and \\\\texttt{classif.ranger} for the propensity score.",
  "Pooled DML columns use respondent clustering via \\\\texttt{DoubleMLClusterData}; DML with respondent fixed effects is not yet reported and those cells are left blank.",
  "Significance: *** \\(p < 0.01\\), ** \\(p < 0.05\\), * \\(p < 0.10\\)."
)

#  run the continuous-treatment results suite and export its tables
build_results_suite(
  treatment_var = continuous_treatment,
  specs = plr_specs,
  table_prefix = "table1_preliminary_ols",
  combined_caption = "Preliminary baseline estimates across wellbeing outcomes",
  per_outcome_caption_prefix = "Preliminary partially linear estimates for",
  per_outcome_label_prefix = "tab:prelim_ols_",
  combined_label = "tab:prelim_ols_combined",
  note_text = plr_note
)

#  build the compact plr- and irm-dml learner-comparison tables
write_plr_method_tables(names(outcomes))
write_irm_method_tables(names(outcomes))
