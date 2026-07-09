# Paired fixed-effects versus PLR-DML comparison for Chapter 5.
#
# The bootstrap resamples respondents rather than person-wave records, so each
# adolescent's complete observed history remains together. Both estimators are
# re-estimated on every bootstrap sample, which preserves their covariance.

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
panel_results_path <- file.path(
  project_dir,
  "tables",
  "chapter5_core",
  "table_panel_benchmarks_numeric.csv"
)
dml_results_path <- file.path(
  project_dir,
  "tables",
  "chapter5_core",
  "table_dml_plr_learner_comparison_numeric.csv"
)
output_dir <- file.path(project_dir, "tables", "chapter5_core")
supported_outcomes <- c(
  "loneliness",
  "life_dissatisfaction",
  "schoolwork_dissatisfaction",
  "school_dissatisfaction"
)
outcome <- Sys.getenv(
  "FE_DML_OUTCOME",
  unset = "life_dissatisfaction"
)
if (!(outcome %in% supported_outcomes)) {
  stop(
    "FE_DML_OUTCOME must be one of: ",
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
table_path <- file.path(
  output_dir,
  paste0("table_fe_plr_comparison_", outcome, "_standalone.tex")
)
draws_path <- file.path(
  output_dir,
  paste0(
    "table_fe_dml_comparison",
    output_suffix,
    "_bootstrap_draws.csv"
  )
)
numeric_path <- file.path(
  output_dir,
  paste0(
    "table_fe_dml_comparison",
    output_suffix,
    "_numeric.csv"
  )
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

treatment <- "social_media_weekday"
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
learner_order <- c("Elastic net", "Lasso", "Random forest", "Neural net")

bootstrap_reps <- as.integer(
  Sys.getenv("FE_DML_BOOTSTRAP_REPS", unset = "500")
)
bootstrap_cores <- as.integer(
  Sys.getenv(
    "FE_DML_BOOTSTRAP_CORES",
    unset = as.character(min(4L, max(1L, parallel::detectCores() - 1L)))
  )
)
bootstrap_seed <- as.integer(
  Sys.getenv("FE_DML_BOOTSTRAP_SEED", unset = "14717867")
)
bootstrap_batch_size <- as.integer(
  Sys.getenv("FE_DML_BOOTSTRAP_BATCH_SIZE", unset = "25")
)

if (is.na(bootstrap_reps) || bootstrap_reps < 2L) {
  stop("FE_DML_BOOTSTRAP_REPS must be an integer of at least 2.")
}
if (is.na(bootstrap_cores) || bootstrap_cores < 1L) {
  stop("FE_DML_BOOTSTRAP_CORES must be a positive integer.")
}
if (is.na(bootstrap_batch_size) || bootstrap_batch_size < 1L) {
  stop("FE_DML_BOOTSTRAP_BATCH_SIZE must be a positive integer.")
}

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
  missing_vars <- setdiff(vars_needed, names(data))
  if (length(missing_vars) > 0L) {
    stop("Missing analysis variables: ", paste(missing_vars, collapse = ", "))
  }

  sample <- data[
    complete.cases(data[, ..vars_needed]),
    ..vars_needed
  ]
  sample <- droplevels(as.data.frame(sample))
  setDT(sample)
  sample
}

make_learners <- function(n_features) {
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

draw_cluster_bootstrap <- function(data, seed) {
  set.seed(seed)
  respondent_ids <- unique(data$pidp)
  sampled_ids <- sample(
    respondent_ids,
    size = length(respondent_ids),
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
  x_matrix <- x_matrix[, apply(x_matrix, 2L, function(x) {
    length(unique(x)) > 1L
  }), drop = FALSE]

  dml_df <- data.table(
    pidp = data$pidp,
    outcome_value = data[[outcome]],
    treatment_value = data[[treatment]]
  )
  setnames(dml_df, c("outcome_value", "treatment_value"), c(outcome, treatment))
  dml_df <- cbind(dml_df, as.data.table(x_matrix))

  list(
    data = dml_df,
    x_cols = setdiff(names(dml_df), c("pidp", outcome, treatment))
  )
}

fit_dml_learners <- function(data, seed) {
  prepared <- prepare_dml_data(data)
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

run_bootstrap_replication <- function(replication, data) {
  replication_seed <- bootstrap_seed + 100000L * replication
  bootstrap_data <- draw_cluster_bootstrap(data, replication_seed)

  fe_estimate <- tryCatch(
    fit_fixed_effects(bootstrap_data),
    error = function(e) NA_real_
  )
  dml_estimates <- fit_dml_learners(
    bootstrap_data,
    replication_seed + 50000L
  )

  data.frame(
    replication = replication,
    learner = names(dml_estimates),
    fe_estimate = fe_estimate,
    dml_estimate = unname(dml_estimates),
    difference = fe_estimate - unname(dml_estimates),
    stringsAsFactors = FALSE
  )
}

read_point_estimates <- function() {
  panel_results <- as.data.frame(fread(panel_results_path))
  dml_results <- as.data.frame(fread(dml_results_path))

  fe_row <- panel_results[
    panel_results[["outcome"]] == outcome &
      panel_results[["model"]] == "Fixed effects + wave FE",
    ,
    drop = FALSE
  ]
  dml_rows <- dml_results[
    dml_results[["outcome"]] == outcome &
      dml_results[["learner"]] %in% learner_order,
    ,
    drop = FALSE
  ]
  dml_rows$learner <- factor(dml_rows$learner, levels = learner_order)
  dml_rows <- dml_rows[order(dml_rows$learner), , drop = FALSE]

  if (nrow(fe_row) != 1L || nrow(dml_rows) != length(learner_order)) {
    stop("Could not recover the baseline FE and DML point estimates.")
  }
  if (any(fe_row$n != dml_rows$n)) {
    stop("The FE and DML point estimates are not based on the same sample size.")
  }

  data.frame(
    learner = as.character(dml_rows$learner),
    fe_estimate = rep(fe_row$estimate[[1]], nrow(dml_rows)),
    dml_estimate = dml_rows$estimate,
    n = dml_rows$n,
    stringsAsFactors = FALSE
  )
}

summarize_bootstrap <- function(draws, point_estimates) {
  summaries <- lapply(learner_order, function(learner_name) {
    learner_draws <- draws[
      draws$learner == learner_name &
        complete.cases(draws[, c("fe_estimate", "dml_estimate")]),
    ]
    point_row <- point_estimates[
      point_estimates$learner == learner_name,
    ]

    if (nrow(learner_draws) < 2L) {
      return(data.frame(
        learner = learner_name,
        fe_estimate = point_row$fe_estimate,
        dml_estimate = point_row$dml_estimate,
        difference = point_row$fe_estimate - point_row$dml_estimate,
        bootstrap_se = NA_real_,
        t_statistic = NA_real_,
        p_value = NA_real_,
        fe_variance = NA_real_,
        dml_variance = NA_real_,
        covariance = NA_real_,
        successful_reps = nrow(learner_draws),
        requested_reps = bootstrap_reps,
        n = point_row$n
      ))
    }

    fe_variance <- var(learner_draws$fe_estimate)
    dml_variance <- var(learner_draws$dml_estimate)
    covariance <- cov(
      learner_draws$fe_estimate,
      learner_draws$dml_estimate
    )
    bootstrap_se <- sqrt(
      fe_variance + dml_variance - 2 * covariance
    )
    point_difference <- point_row$fe_estimate - point_row$dml_estimate
    t_statistic <- point_difference / bootstrap_se
    p_value <- 2 * pnorm(-abs(t_statistic))

    data.frame(
      learner = learner_name,
      fe_estimate = point_row$fe_estimate,
      dml_estimate = point_row$dml_estimate,
      difference = point_difference,
      bootstrap_se = bootstrap_se,
      t_statistic = t_statistic,
      p_value = p_value,
      fe_variance = fe_variance,
      dml_variance = dml_variance,
      covariance = covariance,
      successful_reps = nrow(learner_draws),
      requested_reps = bootstrap_reps,
      n = point_row$n
    )
  })

  rbindlist(summaries)
}

format_number <- function(x, digits = 3L) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

format_p_value <- function(x) {
  ifelse(
    is.na(x),
    "",
    ifelse(x < 0.001, "$< 0.001$", formatC(x, format = "f", digits = 3L))
  )
}

write_latex_table <- function(results) {
  results <- as.data.table(results)
  results[, learner := factor(learner, levels = learner_order)]
  setorder(results, learner)

  row_line <- function(label, values) {
    paste(c(label, values), collapse = " & ")
  }

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    paste0(
      "\\caption{Paired comparison of fixed-effects and PLR-DML estimates for ",
      tolower(outcome_label),
      "}"
    ),
    paste0("\\label{tab:fe_plr_comparison_", outcome, "}"),
    "\\small",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\begin{tabular}{lcccc}",
    "\\hline",
    paste(c("Statistic", learner_order), collapse = " & "),
    "\\\\",
    "\\hline",
    row_line(
      "Fixed-effects estimate",
      format_number(results$fe_estimate)
    ),
    "\\\\",
    row_line(
      "PLR-DML estimate",
      format_number(results$dml_estimate)
    ),
    "\\\\",
    row_line(
      "Difference (FE $-$ DML)",
      format_number(results$difference)
    ),
    "\\\\",
    row_line(
      "Bootstrap standard error",
      format_number(results$bootstrap_se)
    ),
    "\\\\",
    row_line(
      "$t$-statistic",
      format_number(results$t_statistic)
    ),
    "\\\\",
    row_line(
      "$p$-value",
      format_p_value(results$p_value)
    ),
    "\\\\",
    "\\hline",
    "\\end{tabular}",
    "\\begin{minipage}{0.96\\textwidth}",
    paste0(
      "\\footnotesize Notes: The null hypothesis is equality of the fixed-effects and PLR-DML coefficients. ",
      "The fixed-effects estimate is column (5) of Table~\\ref{tab:pooled_ols_life}; the DML estimates are from Table~\\ref{tab:dml_life_methods}. ",
      "The standard error of each difference is estimated from ",
      bootstrap_reps,
      " paired respondent-cluster bootstrap replications. Complete respondent histories are resampled, and both estimators are re-estimated in every replication. ",
      "The two-sided $p$-values use the standard-normal approximation. All estimates use the same outcome-specific sample with $N=",
      format(results$n[[1]], big.mark = ","),
      "$ person-wave observations."
    ),
    "\\end{minipage}",
    "\\end{table}"
  )

  writeLines(lines, table_path)
}

analysis_sample <- prepare_analysis_sample()
point_estimates <- read_point_estimates()

existing_draws <- if (file.exists(draws_path)) {
  fread(draws_path)
} else {
  data.table()
}
if (nrow(existing_draws) > 0L) {
  existing_draws <- existing_draws[
    replication <= bootstrap_reps & learner %in% learner_order
  ]
  complete_replications <- existing_draws[
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
  existing_draws <- existing_draws[
    replication %in% complete_replications
  ]
} else {
  complete_replications <- integer(0)
}
replications_to_run <- setdiff(
  seq_len(bootstrap_reps),
  complete_replications
)

cat(
  sprintf(
    "Target: %d paired respondent-cluster bootstrap replications on %d cores; %d already complete.\n",
    bootstrap_reps,
    bootstrap_cores,
    length(complete_replications)
  )
)

bootstrap_draws <- existing_draws
replication_batches <- split(
  replications_to_run,
  ceiling(seq_along(replications_to_run) / bootstrap_batch_size)
)

for (batch in replication_batches) {
  batch_results <- parallel::mclapply(
    X = batch,
    FUN = run_bootstrap_replication,
    data = analysis_sample,
    mc.cores = bootstrap_cores,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
  batch_draws <- rbindlist(batch_results)
  bootstrap_draws <- rbindlist(
    list(bootstrap_draws, batch_draws),
    use.names = TRUE
  )
  setorder(bootstrap_draws, replication, learner)
  fwrite(bootstrap_draws, draws_path)
  cat(
    sprintf(
      "Completed %d of %d bootstrap replications.\n",
      uniqueN(bootstrap_draws$replication),
      bootstrap_reps
    )
  )
}
comparison_results <- summarize_bootstrap(
  as.data.frame(bootstrap_draws),
  point_estimates
)

fwrite(
  comparison_results,
  numeric_path
)
write_latex_table(comparison_results)

print(comparison_results)
cat("Wrote paired comparison outputs and LaTeX table.\n")
