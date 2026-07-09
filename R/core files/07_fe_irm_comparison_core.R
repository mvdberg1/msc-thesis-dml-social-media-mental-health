# Paired binary-treatment fixed-effects versus IRM-DML comparison.
#
# Respondents are resampled as complete clusters. The binary-treatment FE model
# and all IRM learners are re-estimated on every bootstrap sample.

suppressPackageStartupMessages({
  library(data.table)
  library(DoubleML)
  library(mlr3)
  library(mlr3learners)
  library(plm)
})

lgr::get_logger("mlr3")$set_threshold("warn")
data.table::setDTthreads(1L)
future::plan(future::sequential)

project_dir <- normalizePath(
  Sys.getenv("THESIS_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
data_path <- file.path(
  project_dir,
  "data",
  "analysis",
  "ukhls_youth_l_to_o_clean.rds"
)
supported_outcomes <- c(
  "loneliness",
  "life_dissatisfaction",
  "schoolwork_dissatisfaction",
  "school_dissatisfaction"
)
outcome <- Sys.getenv(
  "FE_IRM_OUTCOME",
  unset = "life_dissatisfaction"
)
if (!(outcome %in% supported_outcomes)) {
  stop(
    "FE_IRM_OUTCOME must be one of: ",
    paste(supported_outcomes, collapse = ", ")
  )
}
outcome_labels <- c(
  loneliness = "Loneliness",
  life_dissatisfaction = "Life dissatisfaction",
  schoolwork_dissatisfaction = "School-work dissatisfaction",
  school_dissatisfaction = "School dissatisfaction"
)
outcome_label <- unname(outcome_labels[[outcome]])
output_suffix <- if (outcome == "life_dissatisfaction") {
  ""
} else {
  paste0("_", outcome)
}
irm_results_path <- file.path(
  project_dir,
  "tables",
  "chapter5_core",
  "table_dml_irm_learner_comparison_numeric.csv"
)
plr_comparison_path <- file.path(
  project_dir,
  "tables",
  "chapter5_core",
  paste0(
    "table_fe_dml_comparison",
    output_suffix,
    "_numeric.csv"
  )
)
output_dir <- file.path(project_dir, "tables", "chapter5_core")
draws_path <- file.path(
  output_dir,
  paste0(
    "table_fe_irm_comparison",
    output_suffix,
    "_bootstrap_draws.csv"
  )
)
numeric_path <- file.path(
  output_dir,
  paste0(
    "table_fe_irm_comparison",
    output_suffix,
    "_numeric.csv"
  )
)
table_path <- file.path(
  project_dir,
  "tables",
  if (outcome == "life_dissatisfaction") {
    "table_fe_dml_comparison.tex"
  } else {
    paste0("table_fe_dml_comparison_", outcome, ".tex")
  }
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

treatment <- "high_social_media_weekday"
controls <- c(
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
learner_order <- c("Elastic net", "Lasso", "Random forest")
table_learner_order <- c(
  "Elastic net",
  "Lasso",
  "Random forest",
  "Neural net"
)

bootstrap_reps <- as.integer(
  Sys.getenv("FE_IRM_BOOTSTRAP_REPS", unset = "500")
)
bootstrap_cores <- as.integer(
  Sys.getenv(
    "FE_IRM_BOOTSTRAP_CORES",
    unset = as.character(min(4L, max(1L, parallel::detectCores() - 1L)))
  )
)
bootstrap_seed <- as.integer(
  Sys.getenv("FE_IRM_BOOTSTRAP_SEED", unset = "14717868")
)
bootstrap_batch_size <- as.integer(
  Sys.getenv("FE_IRM_BOOTSTRAP_BATCH_SIZE", unset = "20")
)

prepare_analysis_sample <- function() {
  data <- as.data.table(readRDS(data_path))
  if (!("wave_f" %in% names(data))) {
    wave_levels <- sort(unique(as.character(na.omit(data$wave))))
    data[, wave_f := factor(
      wave,
      levels = wave_levels,
      labels = paste("Wave", toupper(wave_levels))
    )]
  }

  vars_needed <- unique(c(
    outcome,
    treatment,
    controls,
    "wave_f",
    "pidp",
    "wave"
  ))
  sample <- data[
    complete.cases(data[, ..vars_needed]),
    ..vars_needed
  ]
  sample <- droplevels(as.data.frame(sample))
  setDT(sample)
  sample
}

draw_cluster_bootstrap <- function(data, seed) {
  set.seed(seed)
  respondent_ids <- unique(data$pidp)
  sampled_ids <- sample(
    respondent_ids,
    length(respondent_ids),
    replace = TRUE
  )
  pieces <- lapply(seq_along(sampled_ids), function(draw_id) {
    rows <- data[pidp == sampled_ids[[draw_id]]]
    rows[, bootstrap_pidp := draw_id]
    rows
  })
  bootstrap_data <- rbindlist(pieces, use.names = TRUE)
  bootstrap_data[, pidp := NULL]
  setnames(bootstrap_data, "bootstrap_pidp", "pidp")
  bootstrap_data[, pidp := as.integer(pidp)]
  droplevels(as.data.frame(bootstrap_data))
}

fit_fixed_effects <- function(data) {
  rhs <- c(treatment, controls, "wave_f")
  formula <- as.formula(
    paste(outcome, "~", paste(rhs, collapse = " + "))
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
  unname(coef(fit)[[treatment]])
}

prepare_dml_data <- function(data) {
  x_cols <- c(controls, "wave_f")
  x_matrix <- model.matrix(
    ~ . - 1,
    data = data[, x_cols, drop = FALSE]
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

make_learners <- function(n_features) {
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

fit_irm_learners <- function(data, seed) {
  prepared <- prepare_dml_data(data)
  treatment_counts <- table(prepared$data[[treatment]])
  if (length(treatment_counts) != 2L || any(treatment_counts < 30L)) {
    return(setNames(rep(NA_real_, length(learner_order)), learner_order))
  }

  dml_data <- DoubleMLClusterData$new(
    prepared$data,
    y_col = outcome,
    d_cols = treatment,
    x_cols = prepared$x_cols,
    cluster_cols = "pidp"
  )
  learners <- make_learners(length(prepared$x_cols))

  estimates <- vapply(seq_along(learners), function(j) {
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
      if (!is.null(fit)) {
        break
      }
    }
    if (is.null(fit)) {
      return(NA_real_)
    }
    unname(fit$coef[[treatment]])
  }, numeric(1))
  names(estimates) <- names(learners)
  estimates
}

run_replication <- function(replication, data) {
  seed <- bootstrap_seed + 100000L * replication
  bootstrap_data <- draw_cluster_bootstrap(data, seed)
  fe_estimate <- tryCatch(
    fit_fixed_effects(bootstrap_data),
    error = function(e) NA_real_
  )
  irm_estimates <- fit_irm_learners(
    bootstrap_data,
    seed + 50000L
  )
  data.frame(
    replication = replication,
    learner = names(irm_estimates),
    fe_estimate = fe_estimate,
    dml_estimate = unname(irm_estimates),
    difference = fe_estimate - unname(irm_estimates),
    stringsAsFactors = FALSE
  )
}

read_point_estimates <- function(data) {
  irm_results <- as.data.frame(fread(irm_results_path))
  irm_rows <- irm_results[
    irm_results$outcome == outcome &
      irm_results$learner %in% learner_order,
    ,
    drop = FALSE
  ]
  irm_rows$learner <- factor(irm_rows$learner, levels = learner_order)
  irm_rows <- irm_rows[order(irm_rows$learner), , drop = FALSE]
  fe_estimate <- fit_fixed_effects(as.data.frame(data))

  data.frame(
    learner = as.character(irm_rows$learner),
    fe_estimate = fe_estimate,
    dml_estimate = irm_rows$estimate,
    n = irm_rows$n,
    stringsAsFactors = FALSE
  )
}

summarize_bootstrap <- function(draws, point_estimates) {
  rbindlist(lapply(learner_order, function(learner_name) {
    learner_draws <- draws[
      draws$learner == learner_name &
        is.finite(draws$fe_estimate) &
        is.finite(draws$dml_estimate),
    ]
    point_row <- point_estimates[
      point_estimates$learner == learner_name,
    ]
    bootstrap_se <- sd(learner_draws$difference)
    difference <- point_row$fe_estimate - point_row$dml_estimate
    t_statistic <- difference / bootstrap_se

    data.frame(
      learner = learner_name,
      fe_estimate = point_row$fe_estimate,
      dml_estimate = point_row$dml_estimate,
      difference = difference,
      bootstrap_se = bootstrap_se,
      t_statistic = t_statistic,
      p_value = 2 * pnorm(-abs(t_statistic)),
      fe_variance = var(learner_draws$fe_estimate),
      dml_variance = var(learner_draws$dml_estimate),
      covariance = cov(
        learner_draws$fe_estimate,
        learner_draws$dml_estimate
      ),
      successful_reps = nrow(learner_draws),
      requested_reps = bootstrap_reps,
      n = point_row$n
    )
  }))
}

format_number <- function(x) {
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = 3L))
}

write_combined_table <- function(irm_results) {
  plr_results <- fread(plr_comparison_path)
  dml_table_label <- paste0("tab:dml_", sub("_dissatisfaction$", "", outcome), "_methods")
  irm_table_label <- paste0("tab:irm_", sub("_dissatisfaction$", "", outcome), "_methods")
  if (outcome == "life_dissatisfaction") {
    dml_table_label <- "tab:dml_life_methods"
    irm_table_label <- "tab:irm_life_methods"
  }
  comparison_table_label <- if (outcome == "life_dissatisfaction") {
    "tab:fe_dml_comparison"
  } else {
    paste0("tab:fe_dml_comparison_", outcome)
  }

  get_values <- function(results, variable) {
    setDT(results)
    results[, learner := factor(learner, levels = table_learner_order)]
    setorder(results, learner)
    values <- setNames(results[[variable]], as.character(results$learner))
    format_number(unname(values[table_learner_order]))
  }

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    paste0(
      "\\caption{Paired comparisons of fixed-effects and DML estimates for ",
      tolower(outcome_label),
      "}"
    ),
    paste0("\\label{", comparison_table_label, "}"),
    "\\small",
    "\\setlength{\\tabcolsep}{6pt}",
    "\\begin{tabular}{lcccc}",
    "\\hline",
    paste(c("Statistic", table_learner_order), collapse = " & "),
    "\\\\",
    "\\hline",
    "\\multicolumn{5}{l}{\\textit{Panel A: Ordered-treatment FE versus PLR-DML}} \\\\",
    paste(
      c("$t$-statistic", get_values(plr_results, "t_statistic")),
      collapse = " & "
    ),
    "\\\\",
    paste(
      c("$p$-value", get_values(plr_results, "p_value")),
      collapse = " & "
    ),
    "\\\\[0.4ex]",
    "\\multicolumn{5}{l}{\\textit{Panel B: Binary-treatment FE versus IRM-DML}} \\\\",
    paste(
      c("$t$-statistic", get_values(irm_results, "t_statistic")),
      collapse = " & "
    ),
    "\\\\",
    paste(
      c("$p$-value", get_values(irm_results, "p_value")),
      collapse = " & "
    ),
    "\\\\",
    "\\hline",
    "\\end{tabular}",
    "\\begin{minipage}{0.96\\textwidth}",
    paste0(
      "\\footnotesize Notes: Panel A compares the fixed-effects coefficient on the ordered five-category weekday-use measure with the PLR-DML coefficients in Table~\\ref{",
      dml_table_label,
      "}. Panel B compares a separately estimated fixed-effects coefficient on the binary high-use indicator with the IRM-DML average treatment effects in Table~\\ref{",
      irm_table_label,
      "}. ",
      "The neural-net IRM cell is not reported because that learner was not retained in the binary-treatment analysis. ",
      "For both panels, standard errors of the paired differences use ",
      bootstrap_reps,
      " respondent-cluster bootstrap replications, and the two-sided $p$-values use the standard-normal approximation. ",
      "All comparisons use the same outcome-specific sample with $N=",
      format(irm_results$n[[1]], big.mark = ","),
      "$ person-wave observations."
    ),
    "\\end{minipage}",
    "\\end{table}"
  )
  writeLines(lines, table_path)
}

analysis_sample <- prepare_analysis_sample()
point_estimates <- read_point_estimates(analysis_sample)

existing_draws <- if (file.exists(draws_path)) {
  fread(draws_path)
} else {
  data.table()
}
if (nrow(existing_draws) > 0L) {
  complete_reps <- existing_draws[
    ,
    .(
      n_learners = uniqueN(learner),
      all_finite = all(is.finite(fe_estimate) & is.finite(dml_estimate))
    ),
    by = replication
  ][
    n_learners == length(learner_order) & all_finite,
    replication
  ]
  existing_draws <- existing_draws[replication %in% complete_reps]
} else {
  complete_reps <- integer(0)
}

replications_to_run <- setdiff(seq_len(bootstrap_reps), complete_reps)
cat(
  sprintf(
    "Target: %d paired binary FE-IRM replications on %d cores; %d complete.\n",
    bootstrap_reps,
    bootstrap_cores,
    length(complete_reps)
  )
)

bootstrap_draws <- existing_draws
batches <- split(
  replications_to_run,
  ceiling(seq_along(replications_to_run) / bootstrap_batch_size)
)
for (batch in batches) {
  batch_results <- parallel::mclapply(
    batch,
    run_replication,
    data = analysis_sample,
    mc.cores = bootstrap_cores,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
  bootstrap_draws <- rbindlist(
    list(bootstrap_draws, rbindlist(batch_results)),
    use.names = TRUE
  )
  setorder(bootstrap_draws, replication, learner)
  fwrite(bootstrap_draws, draws_path)
  cat(
    sprintf(
      "Completed %d of %d replications.\n",
      uniqueN(bootstrap_draws$replication),
      bootstrap_reps
    )
  )
}

comparison_results <- summarize_bootstrap(
  as.data.frame(bootstrap_draws),
  point_estimates
)
fwrite(comparison_results, numeric_path)
write_combined_table(comparison_results)
print(comparison_results)
