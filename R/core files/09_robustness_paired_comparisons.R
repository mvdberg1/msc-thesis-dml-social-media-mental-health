# Paired robustness comparisons for weekend exposure and extended waves.
#
# Complete respondent histories are resampled. The fixed-effects and DML
# estimators are re-estimated in each bootstrap sample so their covariance is
# retained when testing equality of the estimates.

suppressPackageStartupMessages({
  library(data.table)
  library(DoubleML)
  library(lmtest)
  library(mlr3)
  library(mlr3learners)
  library(plm)
  library(sandwich)
})

lgr::get_logger("mlr3")$set_threshold("warn")
data.table::setDTthreads(1L)
future::plan(future::sequential)

project_dir <- normalizePath(
  Sys.getenv("THESIS_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
output_dir <- file.path(project_dir, "tables", "chapter5_robustness")
tables_dir <- file.path(project_dir, "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

outcome <- "life_dissatisfaction"
plr_learners <- c("Elastic net", "Lasso", "Random forest", "Neural net")
irm_learners <- c("Elastic net", "Lasso", "Random forest")
table_learners <- plr_learners

rich_controls <- c(
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
common_controls <- rich_controls[!grepl("^ypdevice", rich_controls)]

bootstrap_reps <- as.integer(
  Sys.getenv("ROBUSTNESS_BOOTSTRAP_REPS", unset = "500")
)
detected_cores <- suppressWarnings(parallel::detectCores())
if (is.na(detected_cores) || detected_cores < 2L) detected_cores <- 2L
bootstrap_cores <- as.integer(
  Sys.getenv(
    "ROBUSTNESS_BOOTSTRAP_CORES",
    unset = as.character(min(4L, max(1L, detected_cores - 1L)))
  )
)
bootstrap_batch_size <- as.integer(
  Sys.getenv("ROBUSTNESS_BOOTSTRAP_BATCH_SIZE", unset = "10")
)
bootstrap_seed <- as.integer(
  Sys.getenv("ROBUSTNESS_BOOTSTRAP_SEED", unset = "26717869")
)

star_string <- function(p_value) {
  if (is.na(p_value)) return("")
  if (p_value < 0.01) return("***")
  if (p_value < 0.05) return("**")
  if (p_value < 0.10) return("*")
  ""
}

format_number <- function(x) {
  ifelse(is.na(x), "--", sprintf("%.3f", x))
}

format_p_value <- function(x) {
  ifelse(
    is.na(x),
    "--",
    ifelse(x < 0.001, "$<0.001$", sprintf("%.3f", x))
  )
}

format_estimate <- function(estimate, p_value) {
  if (!is.finite(estimate)) return("--")
  sprintf("%.3f%s", estimate, star_string(p_value))
}

format_se <- function(std_error) {
  ifelse(is.finite(std_error), sprintf("(%.3f)", std_error), "--")
}

prepare_sample <- function(data, ordered_treatment, binary_treatment, controls) {
  data <- as.data.table(data)
  data[, (binary_treatment) := as.integer(get(ordered_treatment) >= 4)]
  vars <- unique(c(
    outcome,
    ordered_treatment,
    binary_treatment,
    controls,
    "wave_f",
    "pidp",
    "wave"
  ))
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0L) {
    stop("Missing variables: ", paste(missing_vars, collapse = ", "))
  }
  sample <- data[complete.cases(data[, ..vars]), ..vars]
  droplevels(as.data.frame(sample))
}

draw_cluster_bootstrap <- function(data, seed) {
  set.seed(seed)
  data <- as.data.table(data)
  respondent_ids <- unique(data$pidp)
  sampled_ids <- sample(
    respondent_ids,
    size = length(respondent_ids),
    replace = TRUE
  )
  mapping <- data.table(
    pidp = sampled_ids,
    bootstrap_pidp = seq_along(sampled_ids)
  )
  bootstrap_data <- data[mapping, on = "pidp", allow.cartesian = TRUE]
  bootstrap_data[, pidp := as.integer(bootstrap_pidp)]
  bootstrap_data[, bootstrap_pidp := NULL]
  droplevels(as.data.frame(bootstrap_data))
}

fit_fixed_effects <- function(
  data,
  treatment,
  controls,
  include_inference = FALSE
) {
  formula <- as.formula(
    paste(
      outcome,
      "~",
      paste(c(treatment, controls, "wave_f"), collapse = " + ")
    )
  )
  pdata <- pdata.frame(
    data,
    index = c("pidp", "wave"),
    drop.index = TRUE
  )
  fit <- plm(
    formula,
    data = pdata,
    model = "within",
    effect = "individual"
  )
  estimate <- unname(coef(fit)[[treatment]])
  if (!include_inference) return(estimate)

  vcov_mat <- vcovHC(
    fit,
    method = "arellano",
    type = "HC1",
    cluster = "group"
  )
  test <- coeftest(fit, vcov. = vcov_mat)
  p_col <- grep("Pr\\(", colnames(test), value = TRUE)[1]
  list(
    estimate = estimate,
    std_error = unname(test[treatment, "Std. Error"]),
    p_value = unname(test[treatment, p_col])
  )
}

prepare_dml_data <- function(data, treatment, controls) {
  x_matrix <- model.matrix(
    ~ . - 1,
    data = data[, c(controls, "wave_f"), drop = FALSE]
  )
  varying <- apply(x_matrix, 2L, function(x) length(unique(x)) > 1L)
  x_matrix <- x_matrix[, varying, drop = FALSE]

  dml_df <- data.table(
    pidp = data$pidp,
    outcome_value = data[[outcome]],
    treatment_value = data[[treatment]]
  )
  setnames(
    dml_df,
    c("outcome_value", "treatment_value"),
    c(outcome, treatment)
  )
  dml_df <- cbind(dml_df, as.data.table(x_matrix))
  list(
    data = dml_df,
    x_cols = setdiff(names(dml_df), c("pidp", outcome, treatment))
  )
}

make_plr_learners <- function(n_features) {
  list(
    "Elastic net" = lrn(
      "regr.cv_glmnet",
      s = "lambda.min",
      alpha = 0.5
    ),
    "Lasso" = lrn(
      "regr.cv_glmnet",
      s = "lambda.min",
      alpha = 1
    ),
    "Random forest" = lrn(
      "regr.ranger",
      num.trees = 300L,
      mtry = max(1L, floor(sqrt(n_features))),
      min.node.size = 5L,
      num.threads = 1L
    ),
    "Neural net" = lrn(
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
      ml_g = lrn("regr.cv_glmnet", s = "lambda.min", alpha = 0.5),
      ml_m = lrn(
        "classif.cv_glmnet",
        s = "lambda.min",
        alpha = 0.5,
        predict_type = "prob"
      )
    ),
    "Lasso" = list(
      ml_g = lrn("regr.cv_glmnet", s = "lambda.min", alpha = 1),
      ml_m = lrn(
        "classif.cv_glmnet",
        s = "lambda.min",
        alpha = 1,
        predict_type = "prob"
      )
    ),
    "Random forest" = list(
      ml_g = lrn(
        "regr.ranger",
        num.trees = 300L,
        mtry = max(1L, floor(sqrt(n_features))),
        min.node.size = 5L,
        num.threads = 1L
      ),
      ml_m = lrn(
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

fit_plr_learners <- function(
  data,
  treatment,
  controls,
  seed,
  include_inference = FALSE
) {
  prepared <- prepare_dml_data(data, treatment, controls)
  dml_data <- DoubleMLClusterData$new(
    prepared$data,
    y_col = outcome,
    d_cols = treatment,
    x_cols = prepared$x_cols,
    cluster_cols = "pidp"
  )
  learners <- make_plr_learners(length(prepared$x_cols))

  rbindlist(lapply(seq_along(learners), function(j) {
    fit <- NULL
    for (attempt in 0:2) {
      set.seed(seed + 1000L * j + 10000000L * attempt)
      fit <- tryCatch(
        {
          object <- DoubleMLPLR$new(
            dml_data,
            learners[[j]]$clone(),
            learners[[j]]$clone(),
            n_folds = 5
          )
          object$fit()
          object
        },
        error = function(e) NULL
      )
      if (!is.null(fit)) break
    }
    data.frame(
      learner = names(learners)[[j]],
      estimate = if (is.null(fit)) NA_real_ else unname(fit$coef[[treatment]]),
      std_error = if (is.null(fit) || !include_inference) {
        NA_real_
      } else {
        unname(fit$se[[treatment]])
      },
      p_value = if (is.null(fit) || !include_inference) {
        NA_real_
      } else {
        unname(fit$pval[[treatment]])
      }
    )
  }))
}

fit_irm_learners <- function(
  data,
  treatment,
  controls,
  seed,
  include_inference = FALSE
) {
  prepared <- prepare_dml_data(data, treatment, controls)
  treatment_counts <- table(prepared$data[[treatment]])
  if (length(treatment_counts) != 2L || any(treatment_counts < 30L)) {
    return(data.table(
      learner = irm_learners,
      estimate = NA_real_,
      std_error = NA_real_,
      p_value = NA_real_
    ))
  }
  dml_data <- DoubleMLClusterData$new(
    prepared$data,
    y_col = outcome,
    d_cols = treatment,
    x_cols = prepared$x_cols,
    cluster_cols = "pidp"
  )
  learners <- make_irm_learners(length(prepared$x_cols))

  rbindlist(lapply(seq_along(learners), function(j) {
    fit <- NULL
    for (attempt in 0:2) {
      set.seed(seed + 1000L * j + 10000000L * attempt)
      fit <- tryCatch(
        {
          object <- DoubleMLIRM$new(
            dml_data,
            learners[[j]]$ml_g$clone(),
            learners[[j]]$ml_m$clone(),
            score = "ATE",
            n_folds = 5
          )
          object$fit()
          object
        },
        error = function(e) NULL
      )
      if (!is.null(fit)) break
    }
    data.frame(
      learner = names(learners)[[j]],
      estimate = if (is.null(fit)) NA_real_ else unname(fit$coef[[treatment]]),
      std_error = if (is.null(fit) || !include_inference) {
        NA_real_
      } else {
        unname(fit$se[[treatment]])
      },
      p_value = if (is.null(fit) || !include_inference) {
        NA_real_
      } else {
        unname(fit$pval[[treatment]])
      }
    )
  }))
}

fit_point_estimates <- function(design) {
  ordered_fe <- fit_fixed_effects(
    design$data,
    design$ordered_treatment,
    design$controls,
    include_inference = TRUE
  )
  binary_fe <- fit_fixed_effects(
    design$data,
    design$binary_treatment,
    design$controls,
    include_inference = TRUE
  )
  plr <- fit_plr_learners(
    design$data,
    design$ordered_treatment,
    design$controls,
    seed = design$seed + 1000L,
    include_inference = TRUE
  )
  irm <- fit_irm_learners(
    design$data,
    design$binary_treatment,
    design$controls,
    seed = design$seed + 2000L,
    include_inference = TRUE
  )
  list(
    ordered_fe = ordered_fe,
    binary_fe = binary_fe,
    plr = plr,
    irm = irm
  )
}

scale_cluster_contributions <- function(
  observation_influence,
  cluster_ids,
  target_std_error
) {
  n <- length(observation_influence)
  contributions <- data.table(
    pidp = cluster_ids,
    influence = observation_influence
  )[
    ,
    .(contribution = sum(influence) / n),
    by = pidp
  ]
  centered <- contributions$contribution - mean(contributions$contribution)
  n_clusters <- nrow(contributions)
  raw_std_error <- sqrt(
    n_clusters / (n_clusters - 1) * sum(centered^2)
  )
  contributions[, contribution := contribution * target_std_error / raw_std_error]
  contributions
}

fit_fixed_effects_influence <- function(data, treatment, controls) {
  formula <- as.formula(
    paste(
      outcome,
      "~",
      paste(c(treatment, controls, "wave_f"), collapse = " + ")
    )
  )
  pdata <- pdata.frame(
    data,
    index = c("pidp", "wave"),
    drop.index = TRUE
  )
  fit <- plm(
    formula,
    data = pdata,
    model = "within",
    effect = "individual"
  )
  vcov_mat <- vcovHC(
    fit,
    method = "arellano",
    type = "HC1",
    cluster = "group"
  )
  test <- coeftest(fit, vcov. = vcov_mat)
  p_col <- grep("Pr\\(", colnames(test), value = TRUE)[1]
  x_matrix <- model.matrix(fit)
  retained_columns <- names(coef(fit))
  x_matrix <- x_matrix[, retained_columns, drop = FALSE]
  target_column <- match(treatment, colnames(x_matrix))
  bread <- solve(crossprod(x_matrix) / nrow(x_matrix))
  observation_influence <- as.numeric(
    (x_matrix * as.numeric(residuals(fit))) %*% bread[, target_column]
  )
  cluster_ids <- as.character(index(fit)[[1]])
  std_error <- unname(test[treatment, "Std. Error"])
  list(
    estimate = unname(test[treatment, "Estimate"]),
    std_error = std_error,
    p_value = unname(test[treatment, p_col]),
    contributions = scale_cluster_contributions(
      observation_influence,
      cluster_ids,
      std_error
    )
  )
}

extract_dml_contributions <- function(fit, treatment, cluster_ids) {
  psi <- as.numeric(fit$psi[, 1L, 1L])
  psi_deriv <- as.numeric(fit$psi_deriv[, 1L, 1L])
  if (length(psi_deriv) == 0L && !is.null(fit$psi_a)) {
    psi_deriv <- as.numeric(fit$psi_a[, 1L, 1L])
  }
  if (length(psi_deriv) == 0L && !is.null(fit$psi_elements$psi_a)) {
    psi_deriv <- as.numeric(fit$psi_elements$psi_a[, 1L, 1L])
  }
  if (length(psi) != length(cluster_ids)) {
    stop(
      "DML score length (",
      length(psi),
      ") differs from observation count (",
      length(cluster_ids),
      ")."
    )
  }
  observation_influence <- -psi / mean(psi_deriv)
  if (any(!is.finite(observation_influence))) {
    stop(
      "Non-finite DML influence scores: psi non-finite=",
      sum(!is.finite(psi)),
      ", derivative non-finite=",
      sum(!is.finite(psi_deriv)),
      ", mean derivative=",
      mean(psi_deriv),
      ", fit members=",
      paste(names(fit), collapse = ",")
    )
  }
  scale_cluster_contributions(
    observation_influence,
    cluster_ids,
    unname(fit$se[[treatment]])
  )
}

fit_plr_influence <- function(data, treatment, controls, seed) {
  prepared <- prepare_dml_data(data, treatment, controls)
  original_cluster_ids <- as.character(prepared$data$pidp)
  dml_data <- DoubleMLClusterData$new(
    prepared$data,
    y_col = outcome,
    d_cols = treatment,
    x_cols = prepared$x_cols,
    cluster_cols = "pidp"
  )
  learners <- make_plr_learners(length(prepared$x_cols))

  lapply(seq_along(learners), function(j) {
    set.seed(seed + 1000L * j)
    fit <- DoubleMLPLR$new(
      dml_data,
      learners[[j]]$clone(),
      learners[[j]]$clone(),
      n_folds = 5
    )
    fit$fit()
    list(
      learner = names(learners)[[j]],
      estimate = unname(fit$coef[[treatment]]),
      std_error = unname(fit$se[[treatment]]),
      p_value = unname(fit$pval[[treatment]]),
      contributions = extract_dml_contributions(
        fit,
        treatment,
        original_cluster_ids
      )
    )
  })
}

fit_irm_influence <- function(data, treatment, controls, seed) {
  prepared <- prepare_dml_data(data, treatment, controls)
  original_cluster_ids <- as.character(prepared$data$pidp)
  dml_data <- DoubleMLClusterData$new(
    prepared$data,
    y_col = outcome,
    d_cols = treatment,
    x_cols = prepared$x_cols,
    cluster_cols = "pidp"
  )
  learners <- make_irm_learners(length(prepared$x_cols))

  lapply(seq_along(learners), function(j) {
    set.seed(seed + 1000L * j)
    fit <- DoubleMLIRM$new(
      dml_data,
      learners[[j]]$ml_g$clone(),
      learners[[j]]$ml_m$clone(),
      score = "ATE",
      n_folds = 5
    )
    fit$fit()
    list(
      learner = names(learners)[[j]],
      estimate = unname(fit$coef[[treatment]]),
      std_error = unname(fit$se[[treatment]]),
      p_value = unname(fit$pval[[treatment]]),
      contributions = extract_dml_contributions(
        fit,
        treatment,
        original_cluster_ids
      )
    )
  })
}

paired_influence_summary <- function(
  dml_fit,
  fe_fit,
  cluster_ids,
  multiplier_counts
) {
  fe_contributions <- setNames(
    fe_fit$contributions$contribution,
    as.character(fe_fit$contributions$pidp)
  )[cluster_ids]
  dml_contributions <- setNames(
    dml_fit$contributions$contribution,
    as.character(dml_fit$contributions$pidp)
  )[cluster_ids]
  # Singleton respondents have zero within-FE contribution and can be absent
  # from the transformed plm model matrix.
  fe_contributions[is.na(fe_contributions)] <- 0
  if (anyNA(dml_contributions)) {
    missing_ids <- cluster_ids[is.na(dml_contributions)]
    stop(
      "DML influence contributions could not be aligned for ",
      length(missing_ids),
      " respondents; examples: ",
      paste(head(missing_ids), collapse = ", "),
      ". Available contribution IDs: ",
      paste(head(as.character(dml_fit$contributions$pidp)), collapse = ", "),
      ". ID overlap: ",
      length(intersect(
        cluster_ids,
        as.character(dml_fit$contributions$pidp)
      ))
    )
  }
  centered_counts <- multiplier_counts - 1
  fe_perturbations <- as.numeric(centered_counts %*% fe_contributions)
  dml_perturbations <- as.numeric(centered_counts %*% dml_contributions)
  difference_draws <- (
    fe_fit$estimate + fe_perturbations
  ) - (
    dml_fit$estimate + dml_perturbations
  )
  difference <- fe_fit$estimate - dml_fit$estimate
  bootstrap_se <- sd(difference_draws)
  t_statistic <- difference / bootstrap_se
  data.frame(
    learner = dml_fit$learner,
    fe_estimate = fe_fit$estimate,
    fe_std_error = fe_fit$std_error,
    fe_p_value = fe_fit$p_value,
    dml_estimate = dml_fit$estimate,
    dml_std_error = dml_fit$std_error,
    dml_p_value = dml_fit$p_value,
    difference = difference,
    bootstrap_se = bootstrap_se,
    t_statistic = t_statistic,
    p_value = 2 * pnorm(-abs(t_statistic)),
    successful_reps = bootstrap_reps,
    requested_reps = bootstrap_reps
  )
}

run_influence_design <- function(design) {
  ordered_fe <- fit_fixed_effects_influence(
    design$data,
    design$ordered_treatment,
    design$controls
  )
  binary_fe <- fit_fixed_effects_influence(
    design$data,
    design$binary_treatment,
    design$controls
  )
  plr_fits <- fit_plr_influence(
    design$data,
    design$ordered_treatment,
    design$controls,
    design$seed + 1000L
  )
  irm_fits <- fit_irm_influence(
    design$data,
    design$binary_treatment,
    design$controls,
    design$seed + 2000L
  )

  cluster_ids <- sort(unique(as.character(design$data$pidp)))
  n_clusters <- length(cluster_ids)
  set.seed(bootstrap_seed + design$seed)
  multiplier_counts <- t(replicate(
    bootstrap_reps,
    tabulate(
      sample.int(n_clusters, n_clusters, replace = TRUE),
      nbins = n_clusters
    )
  ))

  plr_results <- rbindlist(lapply(
    plr_fits,
    paired_influence_summary,
    fe_fit = ordered_fe,
    cluster_ids = cluster_ids,
    multiplier_counts = multiplier_counts
  ))
  plr_results[, panel := "PLR"]
  irm_results <- rbindlist(lapply(
    irm_fits,
    paired_influence_summary,
    fe_fit = binary_fe,
    cluster_ids = cluster_ids,
    multiplier_counts = multiplier_counts
  ))
  irm_results[, panel := "IRM"]
  results <- rbindlist(list(plr_results, irm_results), use.names = TRUE)
  results[, `:=`(
    design = design$name,
    n = nrow(design$data),
    n_clusters = n_clusters,
    bootstrap_type = "paired cluster influence-function bootstrap"
  )]
  fwrite(
    results,
    file.path(
      output_dir,
      paste0(
        "robustness_",
        design$name,
        "_paired_comparison_numeric.csv"
      )
    )
  )
  list(
    results = results,
    ordered_fe = ordered_fe,
    binary_fe = binary_fe,
    plr_fits = plr_fits,
    irm_fits = irm_fits
  )
}

run_replication <- function(replication, design, run_plr) {
  seed <- bootstrap_seed + design$seed + 100000L * replication
  bootstrap_data <- draw_cluster_bootstrap(design$data, seed)

  ordered_fe <- if (run_plr) {
    tryCatch(
      fit_fixed_effects(
        bootstrap_data,
        design$ordered_treatment,
        design$controls
      ),
      error = function(e) NA_real_
    )
  } else {
    NA_real_
  }
  plr <- if (run_plr) {
    fit_plr_learners(
      bootstrap_data,
      design$ordered_treatment,
      design$controls,
      seed + 10000L
    )
  } else {
    data.table()
  }

  binary_fe <- tryCatch(
    fit_fixed_effects(
      bootstrap_data,
      design$binary_treatment,
      design$controls
    ),
    error = function(e) NA_real_
  )
  irm <- fit_irm_learners(
    bootstrap_data,
    design$binary_treatment,
    design$controls,
    seed + 20000L
  )

  ordered_rows <- if (run_plr) {
    plr[, .(
      replication,
      panel = "PLR",
      learner,
      fe_estimate = ordered_fe,
      dml_estimate = estimate,
      difference = ordered_fe - estimate
    )]
  } else {
    data.table()
  }
  binary_rows <- irm[, .(
    replication,
    panel = "IRM",
    learner,
    fe_estimate = binary_fe,
    dml_estimate = estimate,
    difference = binary_fe - estimate
  )]
  rbindlist(list(ordered_rows, binary_rows), use.names = TRUE)
}

summarize_panel <- function(draws, point_estimates, panel_name) {
  learner_names <- if (panel_name == "PLR") plr_learners else irm_learners
  fe_point <- if (panel_name == "PLR") {
    point_estimates$ordered_fe
  } else {
    point_estimates$binary_fe
  }
  dml_points <- if (panel_name == "PLR") {
    point_estimates$plr
  } else {
    point_estimates$irm
  }

  rbindlist(lapply(learner_names, function(learner_name) {
    learner_draws <- draws[
      panel == panel_name &
        learner == learner_name &
        is.finite(fe_estimate) &
        is.finite(dml_estimate)
    ]
    point_row <- dml_points[learner == learner_name]
    difference <- fe_point$estimate - point_row$estimate
    bootstrap_se <- sd(learner_draws$difference)
    t_statistic <- difference / bootstrap_se
    data.frame(
      panel = panel_name,
      learner = learner_name,
      fe_estimate = fe_point$estimate,
      fe_std_error = fe_point$std_error,
      fe_p_value = fe_point$p_value,
      dml_estimate = point_row$estimate,
      dml_std_error = point_row$std_error,
      dml_p_value = point_row$p_value,
      difference = difference,
      bootstrap_se = bootstrap_se,
      t_statistic = t_statistic,
      p_value = 2 * pnorm(-abs(t_statistic)),
      successful_reps = nrow(learner_draws),
      requested_reps = bootstrap_reps
    )
  }))
}

get_table_values <- function(results, variable, learners = table_learners) {
  values <- setNames(results[[variable]], results$learner)
  unname(values[learners])
}

write_weekend_table <- function(results, n) {
  plr <- results[panel == "PLR"]
  irm <- results[panel == "IRM"]
  ordered_estimate <- c(
    format_estimate(plr$fe_estimate[[1]], plr$fe_p_value[[1]]),
    mapply(format_estimate, plr$dml_estimate, plr$dml_p_value)
  )
  ordered_se <- c(
    format_se(plr$fe_std_error[[1]]),
    format_se(plr$dml_std_error)
  )
  binary_estimate <- c(
    format_estimate(irm$fe_estimate[[1]], irm$fe_p_value[[1]]),
    mapply(format_estimate, irm$dml_estimate, irm$dml_p_value),
    "--"
  )
  binary_se <- c(
    format_se(irm$fe_std_error[[1]]),
    format_se(irm$dml_std_error),
    "--"
  )

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Weekend social-media use and paired fixed-effects--DML comparisons}",
    "\\label{tab:robustness_weekend}",
    "\\small",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lccccc}",
    "\\hline",
    paste(c("Statistic", "Fixed effects", table_learners), collapse = " & "),
    "\\\\",
    "\\hline",
    "\\multicolumn{6}{l}{\\textit{Panel A: Ordered-treatment FE versus PLR-DML}} \\\\",
    paste(c("Weekend-use estimate", ordered_estimate), collapse = " & "),
    "\\\\",
    paste(c("", ordered_se), collapse = " & "),
    "\\\\",
    paste(
      c("$t$-statistic", "--", format_number(plr$t_statistic)),
      collapse = " & "
    ),
    "\\\\",
    paste(
      c("$p$-value", "--", format_p_value(plr$p_value)),
      collapse = " & "
    ),
    "\\\\[0.5ex]",
    "\\multicolumn{6}{l}{\\textit{Panel B: Binary-treatment FE versus IRM-DML}} \\\\",
    paste(c("High weekend use estimate", binary_estimate), collapse = " & "),
    "\\\\",
    paste(c("", binary_se), collapse = " & "),
    "\\\\",
    paste(
      c("$t$-statistic", "--", format_number(irm$t_statistic), "--"),
      collapse = " & "
    ),
    "\\\\",
    paste(
      c("$p$-value", "--", format_p_value(irm$p_value), "--"),
      collapse = " & "
    ),
    "\\\\",
    "\\hline",
    "\\end{tabular*}",
    "\\begin{minipage}{0.97\\textwidth}",
    paste0(
      "\\footnotesize Notes: Standard errors are reported in parentheses ",
      "below the coefficients. Panel A uses the ordered five-category ",
      "weekend-use measure. Panel B defines high weekend use as at least ",
      "four hours and reports the IRM average treatment effect. The fixed-",
      "effects models include respondent and wave effects; DML uses five-fold ",
      "respondent-level cross-fitting. The neural-net IRM is not reported. ",
      "Panel A uses ",
      bootstrap_reps,
      " full respondent-cluster bootstrap re-estimations. Panel B uses ",
      bootstrap_reps,
      " paired respondent-cluster influence-function bootstrap draws. The ",
      "complete-case sample ",
      "contains $N=",
      format(n, big.mark = ","),
      "$ person-wave observations. $^{***}p<0.01$, $^{**}p<0.05$, ",
      "$^{*}p<0.10$."
    ),
    "\\end{minipage}",
    "\\end{table}"
  )
  writeLines(lines, file.path(tables_dir, "table_robustness_weekend.tex"))
}

write_extended_table <- function(
  results,
  extended_estimates,
  n_all_wave
) {
  plr <- results[panel == "PLR"]
  irm <- results[panel == "IRM"]
  extended_estimates <- as.data.table(extended_estimates)

  specification_order <- c("lo_rich", "lo_common", "ao_common")
  specification_labels <- c(
    lo_rich = "Waves l--o, rich controls",
    lo_common = "Waves l--o, common-core controls",
    ao_common = "Waves a--o, common-core controls"
  )
  estimate_rows <- unlist(lapply(specification_order, function(specification_name) {
    rows <- extended_estimates[specification == specification_name]
    values <- setNames(rows$estimate, rows$model)
    p_values <- setNames(rows$p_value, rows$model)
    std_errors <- setNames(rows$std_error, rows$model)
    estimates <- c(
      format_estimate(values[["Fixed effects"]], p_values[["Fixed effects"]]),
      vapply(
        table_learners,
        function(model) format_estimate(values[[model]], p_values[[model]]),
        character(1)
      )
    )
    ses <- c(
      format_se(std_errors[["Fixed effects"]]),
      vapply(
        table_learners,
        function(model) format_se(std_errors[[model]]),
        character(1)
      )
    )
    c(
      paste(
        c(specification_labels[[specification_name]], estimates),
        collapse = " & "
      ),
      "\\\\",
      paste(c("", ses), collapse = " & "),
      "\\\\"
    )
  }))

  binary_estimate <- c(
    format_estimate(irm$fe_estimate[[1]], irm$fe_p_value[[1]]),
    mapply(format_estimate, irm$dml_estimate, irm$dml_p_value),
    "--"
  )
  binary_se <- c(
    format_se(irm$fe_std_error[[1]]),
    format_se(irm$dml_std_error),
    "--"
  )

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Wave coverage, controls, and paired all-wave DML comparisons}",
    "\\label{tab:robustness_extended_waves}",
    "\\small",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lccccc}",
    "\\hline",
    paste(c("Specification", "Fixed effects", table_learners), collapse = " & "),
    "\\\\",
    "\\hline",
    "\\multicolumn{6}{l}{\\textit{Panel A: Ordered-treatment FE versus PLR-DML}} \\\\",
    estimate_rows,
    paste(
      c(
        "$t$-statistic, Waves a--o",
        "--",
        format_number(plr$t_statistic)
      ),
      collapse = " & "
    ),
    "\\\\",
    paste(
      c(
        "$p$-value, Waves a--o",
        "--",
        format_p_value(plr$p_value)
      ),
      collapse = " & "
    ),
    "\\\\[0.5ex]",
    "\\multicolumn{6}{l}{\\textit{Panel B: Binary-treatment FE versus IRM-DML}} \\\\",
    paste(c("Waves a--o, high weekday use", binary_estimate), collapse = " & "),
    "\\\\",
    paste(c("", binary_se), collapse = " & "),
    "\\\\",
    paste(
      c("$t$-statistic", "--", format_number(irm$t_statistic), "--"),
      collapse = " & "
    ),
    "\\\\",
    paste(
      c("$p$-value", "--", format_p_value(irm$p_value), "--"),
      collapse = " & "
    ),
    "\\\\",
    "\\hline",
    "\\end{tabular*}",
    "\\begin{minipage}{0.97\\textwidth}",
    paste0(
      "\\footnotesize Notes: Standard errors are reported in parentheses ",
      "below the coefficients. Panel A preserves the control-set and wave-",
      "coverage comparison and reports paired tests for the Waves a--o ",
      "common-core specification. Panel B uses the same all-wave sample and ",
      "defines high weekday use as at least four hours. The neural-net IRM is ",
      "not reported. Standard errors of the paired differences use ",
      bootstrap_reps,
      " paired respondent-cluster influence-function bootstrap draws. The ",
      "Waves l--o rich-control sample contains $N=4,316$, the Waves l--o ",
      "common-core sample contains $N=5,001$, and the Waves a--o common-core ",
      "sample contains $N=",
      format(n_all_wave, big.mark = ","),
      "$ person-wave observations. $^{***}p<0.01$, $^{**}p<0.05$, ",
      "$^{*}p<0.10$."
    ),
    "\\end{minipage}",
    "\\end{table}"
  )
  writeLines(
    lines,
    file.path(tables_dir, "table_robustness_extended_waves.tex")
  )
}

run_design <- function(design, run_plr = TRUE) {
  draws_path <- file.path(
    output_dir,
    paste0("robustness_", design$name, "_paired_bootstrap_draws.csv")
  )
  numeric_path <- file.path(
    output_dir,
    paste0("robustness_", design$name, "_paired_comparison_numeric.csv")
  )

  point_estimates <- fit_point_estimates(design)
  existing_draws <- if (file.exists(draws_path)) {
    fread(draws_path)
  } else {
    data.table()
  }
  required_panels <- if (run_plr) c("PLR", "IRM") else "IRM"
  required_rows <- if (run_plr) {
    length(plr_learners) + length(irm_learners)
  } else {
    length(irm_learners)
  }

  if (nrow(existing_draws) > 0L) {
    existing_draws <- existing_draws[
      replication <= bootstrap_reps & panel %in% required_panels
    ]
    complete_reps <- existing_draws[
      ,
      .(
        n_rows = .N,
        n_panels = uniqueN(panel),
        all_finite = all(is.finite(fe_estimate) & is.finite(dml_estimate))
      ),
      by = replication
    ][
      n_rows == required_rows &
        n_panels == length(required_panels) &
        all_finite,
      replication
    ]
    existing_draws <- existing_draws[replication %in% complete_reps]
  } else {
    complete_reps <- integer(0)
  }

  replications_to_run <- setdiff(seq_len(bootstrap_reps), complete_reps)
  draws <- existing_draws
  batches <- split(
    replications_to_run,
    ceiling(seq_along(replications_to_run) / bootstrap_batch_size)
  )
  cat(sprintf(
    "%s: target %d replications; %d already complete.\n",
    design$name,
    bootstrap_reps,
    length(complete_reps)
  ))

  for (batch in batches) {
    batch_results <- parallel::mclapply(
      batch,
      run_replication,
      design = design,
      run_plr = run_plr,
      mc.cores = bootstrap_cores,
      mc.preschedule = FALSE,
      mc.set.seed = FALSE
    )
    draws <- rbindlist(
      list(draws, rbindlist(batch_results)),
      use.names = TRUE
    )
    setorder(draws, replication, panel, learner)
    fwrite(draws, draws_path)
    cat(sprintf(
      "%s: completed %d of %d replications.\n",
      design$name,
      uniqueN(draws$replication),
      bootstrap_reps
    ))
  }

  summaries <- list()
  if (run_plr) {
    summaries$plr <- summarize_panel(draws, point_estimates, "PLR")
  }
  summaries$irm <- summarize_panel(draws, point_estimates, "IRM")
  results <- rbindlist(summaries, use.names = TRUE)
  results[, `:=`(
    design = design$name,
    n = nrow(design$data),
    n_clusters = uniqueN(design$data$pidp)
  )]
  fwrite(results, numeric_path)
  list(
    results = results,
    point_estimates = point_estimates
  )
}

baseline <- readRDS(file.path(
  project_dir,
  "data",
  "analysis",
  "ukhls_youth_l_to_o_clean.rds"
))
all_wave <- readRDS(file.path(
  project_dir,
  "data",
  "analysis",
  "ukhls_youth_a_to_o_common_core.rds"
))

weekend_design <- list(
  name = "weekend",
  data = prepare_sample(
    baseline,
    "social_media_weekend",
    "high_social_media_weekend",
    rich_controls
  ),
  ordered_treatment = "social_media_weekend",
  binary_treatment = "high_social_media_weekend",
  controls = rich_controls,
  seed = 10000L
)
extended_design <- list(
  name = "extended_waves",
  data = prepare_sample(
    all_wave,
    "social_media_weekday",
    "high_social_media_weekday",
    common_controls
  ),
  ordered_treatment = "social_media_weekday",
  binary_treatment = "high_social_media_weekday",
  controls = common_controls,
  seed = 20000L
)

# The existing weekend PLR bootstrap already contains 500 full re-estimation
# draws. Reuse its summary and estimate the missing IRM panel with the paired
# cluster influence-function bootstrap.
weekend_run <- run_influence_design(weekend_design)
weekend_plr <- fread(file.path(
  output_dir,
  "robustness_weekend_fe_dml_comparison_numeric.csv"
))
weekend_plr[, `:=`(
  panel = "PLR",
  fe_p_value = 2 * pnorm(-abs(fe_estimate / fe_std_error)),
  dml_p_value = 2 * pnorm(-abs(dml_estimate / dml_std_error)),
  design = "weekend",
  n_clusters = uniqueN(weekend_design$data$pidp)
)]
weekend_results <- rbindlist(
  list(weekend_plr, weekend_run$results[panel == "IRM"]),
  use.names = TRUE,
  fill = TRUE
)
fwrite(
  weekend_results,
  file.path(
    output_dir,
    "robustness_weekend_combined_comparison_numeric.csv"
  )
)
write_weekend_table(weekend_results, nrow(weekend_design$data))

extended_run <- run_influence_design(extended_design)
extended_estimates <- fread(file.path(
  output_dir,
  "robustness_extended_waves_numeric.csv"
))
all_wave_plr <- extended_run$results[panel == "PLR"]
all_wave_fe <- all_wave_plr$fe_estimate[[1]]
all_wave_fe_se <- all_wave_plr$fe_std_error[[1]]
all_wave_fe_p <- all_wave_plr$fe_p_value[[1]]
extended_estimates <- extended_estimates[specification != "ao_common"]
extended_estimates <- rbindlist(
  list(
    extended_estimates,
    data.table(
      model = c("Fixed effects", all_wave_plr$learner),
      estimate = c(all_wave_fe, all_wave_plr$dml_estimate),
      std_error = c(all_wave_fe_se, all_wave_plr$dml_std_error),
      p_value = c(all_wave_fe_p, all_wave_plr$dml_p_value),
      n = nrow(extended_design$data),
      n_clusters = uniqueN(extended_design$data$pidp),
      specification = "ao_common"
    )
  ),
  use.names = TRUE,
  fill = TRUE
)
write_extended_table(
  extended_run$results,
  extended_estimates,
  nrow(extended_design$data)
)

print(weekend_results)
print(extended_run$results)
cat("Wrote paired robustness comparisons and LaTeX tables.\n")
