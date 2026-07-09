# Weekend-use robustness and paired fixed-effects versus PLR-DML tests.
#
# Respondents, rather than person-wave records, are resampled so that complete
# panel histories remain together. Both estimators are re-estimated in every
# bootstrap replication, preserving their sampling covariance.

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
data_path <- file.path(
  project_dir,
  "data",
  "analysis",
  "ukhls_youth_l_to_o_clean.rds"
)
output_dir <- file.path(project_dir, "tables", "chapter5_robustness")
table_path <- file.path(project_dir, "tables", "table_robustness_weekend.tex")
draws_path <- file.path(
  output_dir,
  "robustness_weekend_fe_dml_bootstrap_draws.csv"
)
numeric_path <- file.path(
  output_dir,
  "robustness_weekend_fe_dml_comparison_numeric.csv"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

outcome <- "life_dissatisfaction"
treatment <- "social_media_weekend"
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
  Sys.getenv("WEEKEND_BOOTSTRAP_REPS", unset = "500")
)
detected_cores <- suppressWarnings(parallel::detectCores())
if (is.na(detected_cores) || detected_cores < 2L) detected_cores <- 2L
bootstrap_cores <- as.integer(
  Sys.getenv(
    "WEEKEND_BOOTSTRAP_CORES",
    unset = as.character(min(4L, max(1L, detected_cores - 1L)))
  )
)
bootstrap_seed <- as.integer(
  Sys.getenv("WEEKEND_BOOTSTRAP_SEED", unset = "14717867")
)
bootstrap_batch_size <- as.integer(
  Sys.getenv("WEEKEND_BOOTSTRAP_BATCH_SIZE", unset = "25")
)

prepare_analysis_sample <- function() {
  data <- as.data.table(readRDS(data_path))
  vars <- unique(c(
    outcome,
    treatment,
    controls,
    "wave_f",
    "pidp",
    "wave"
  ))
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0L) {
    stop("Missing analysis variables: ", paste(missing_vars, collapse = ", "))
  }
  sample <- data[complete.cases(data[, ..vars]), ..vars]
  droplevels(as.data.frame(sample))
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

fit_fixed_effects <- function(data, include_inference = FALSE) {
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

prepare_dml_data <- function(data) {
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

fit_dml_learners <- function(data, seed, include_inference = FALSE) {
  prepared <- prepare_dml_data(data)
  dml_data <- DoubleMLClusterData$new(
    prepared$data,
    y_col = outcome,
    d_cols = treatment,
    x_cols = prepared$x_cols,
    cluster_cols = "pidp"
  )
  learners <- make_learners(length(prepared$x_cols))

  rows <- lapply(seq_along(learners), function(j) {
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

    if (is.null(fit)) {
      return(data.frame(
        learner = names(learners)[[j]],
        estimate = NA_real_,
        std_error = NA_real_,
        p_value = NA_real_
      ))
    }

    data.frame(
      learner = names(learners)[[j]],
      estimate = unname(fit$coef[[treatment]]),
      std_error = if (include_inference) {
        unname(fit$se[[treatment]])
      } else {
        NA_real_
      },
      p_value = if (include_inference) {
        unname(fit$pval[[treatment]])
      } else {
        NA_real_
      }
    )
  })

  rbindlist(rows)
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
  pieces <- lapply(seq_along(sampled_ids), function(draw_id) {
    rows <- copy(data[pidp == sampled_ids[[draw_id]]])
    rows[, pidp := draw_id]
    rows
  })
  droplevels(as.data.frame(rbindlist(pieces, use.names = TRUE)))
}

run_bootstrap_replication <- function(replication, data) {
  replication_seed <- bootstrap_seed + 100000L * replication
  bootstrap_data <- draw_cluster_bootstrap(data, replication_seed)
  fe_estimate <- tryCatch(
    fit_fixed_effects(bootstrap_data),
    error = function(e) NA_real_
  )
  dml_rows <- fit_dml_learners(
    bootstrap_data,
    replication_seed + 50000L
  )
  dml_rows[, `:=`(
    replication = replication,
    fe_estimate = fe_estimate,
    dml_estimate = estimate,
    difference = fe_estimate - estimate
  )]
  dml_rows[, .(
    replication,
    learner,
    fe_estimate,
    dml_estimate,
    difference
  )]
}

summarize_bootstrap <- function(draws, fe_point, dml_points, n) {
  rbindlist(lapply(learner_order, function(learner_name) {
    learner_draws <- draws[
      learner == learner_name &
        is.finite(fe_estimate) &
        is.finite(dml_estimate)
    ]
    dml_point <- dml_points[learner == learner_name]
    difference <- fe_point$estimate - dml_point$estimate
    bootstrap_se <- sd(learner_draws$difference)
    t_statistic <- difference / bootstrap_se
    p_value <- 2 * pnorm(-abs(t_statistic))
    data.frame(
      learner = learner_name,
      fe_estimate = fe_point$estimate,
      fe_std_error = fe_point$std_error,
      dml_estimate = dml_point$estimate,
      dml_std_error = dml_point$std_error,
      difference = difference,
      bootstrap_se = bootstrap_se,
      t_statistic = t_statistic,
      p_value = p_value,
      successful_reps = nrow(learner_draws),
      requested_reps = bootstrap_reps,
      n = n
    )
  }))
}

star_string <- function(p_value) {
  if (is.na(p_value)) return("")
  if (p_value < 0.01) return("***")
  if (p_value < 0.05) return("**")
  if (p_value < 0.10) return("*")
  ""
}

format_estimate <- function(estimate, std_error, p_value) {
  sprintf(
    "%.3f%s (%.3f)",
    estimate,
    star_string(p_value),
    std_error
  )
}

format_number <- function(x) {
  ifelse(is.na(x), "--", sprintf("%.3f", x))
}

write_latex_table <- function(results, fe_point, dml_points, n) {
  results <- as.data.table(results)
  dml_points <- as.data.table(dml_points)
  results[, learner := factor(learner, levels = learner_order)]
  dml_points[, learner := factor(learner, levels = learner_order)]
  setorder(results, learner)
  setorder(dml_points, learner)

  estimate_cells <- c(
    format_estimate(
      fe_point$estimate,
      fe_point$std_error,
      fe_point$p_value
    ),
    mapply(
      format_estimate,
      dml_points$estimate,
      dml_points$std_error,
      dml_points$p_value
    )
  )
  t_cells <- c("--", format_number(results$t_statistic))
  p_cells <- c("--", format_number(results$p_value))

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Weekend social-media use and paired fixed-effects--DML comparisons}",
    "\\label{tab:robustness_weekend}",
    "\\small",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{lccccc}",
    "\\hline",
    paste(
      c(
        "Statistic",
        "Fixed effects",
        "Elastic net",
        "Lasso",
        "Random forest",
        "Neural net"
      ),
      collapse = " & "
    ),
    "\\\\",
    "\\hline",
    paste(c("Weekend-use estimate", estimate_cells), collapse = " & "),
    "\\\\",
    paste(c("$t$-statistic", t_cells), collapse = " & "),
    "\\\\",
    paste(c("$p$-value", p_cells), collapse = " & "),
    "\\\\",
    "\\hline",
    "\\end{tabular}",
    "}",
    "\\begin{minipage}{0.97\\textwidth}",
    paste0(
      "\\footnotesize Notes: The first row reports coefficients with ",
      "respondent-clustered standard errors in parentheses. The fixed-effects ",
      "model includes respondent and wave effects; the DML columns use ",
      "five-fold respondent-level cross-fitting. The $t$-statistics test ",
      "equality of the fixed-effects coefficient and each PLR-DML coefficient. ",
      "Standard errors of the paired differences use ",
      bootstrap_reps,
      " respondent-cluster bootstrap replications in which both estimators ",
      "are re-estimated. Two-sided $p$-values use the standard-normal ",
      "approximation. The weekend complete-case sample contains $N=",
      format(n, big.mark = ","),
      "$ person-wave observations. ",
      "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$."
    ),
    "\\end{minipage}",
    "\\end{table}"
  )
  writeLines(lines, table_path)
}

analysis_sample <- prepare_analysis_sample()
fe_point <- fit_fixed_effects(analysis_sample, include_inference = TRUE)
dml_points <- fit_dml_learners(
  analysis_sample,
  seed = 2100L,
  include_inference = TRUE
)

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
  existing_draws <- existing_draws[replication %in% complete_replications]
} else {
  complete_replications <- integer(0)
}

replications_to_run <- setdiff(seq_len(bootstrap_reps), complete_replications)
bootstrap_draws <- existing_draws
replication_batches <- split(
  replications_to_run,
  ceiling(seq_along(replications_to_run) / bootstrap_batch_size)
)

cat(sprintf(
  "Target: %d weekend bootstrap replications; %d already complete.\n",
  bootstrap_reps,
  length(complete_replications)
))

for (batch in replication_batches) {
  batch_results <- parallel::mclapply(
    batch,
    run_bootstrap_replication,
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
  cat(sprintf(
    "Completed %d of %d weekend bootstrap replications.\n",
    uniqueN(bootstrap_draws$replication),
    bootstrap_reps
  ))
}

comparison_results <- summarize_bootstrap(
  bootstrap_draws,
  fe_point,
  dml_points,
  nrow(analysis_sample)
)
fwrite(comparison_results, numeric_path)
write_latex_table(
  comparison_results,
  fe_point,
  dml_points,
  nrow(analysis_sample)
)

print(comparison_results)
cat("Wrote weekend robustness comparison outputs.\n")
