# DML-specific robustness checks for life dissatisfaction.
#
# The script evaluates:
# 1. 100 repeated five-fold cross-fitting partitions for PLR;
# 2. propensity-score overlap and IRM propensity truncation;
# 3. PLR omitted-variable robustness values from cross-fitted residuals.

suppressPackageStartupMessages({
  library(data.table)
  library(DoubleML)
  library(ggplot2)
  library(mlr3)
  library(mlr3learners)
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
tables_dir <- file.path(project_dir, "tables")
figures_dir <- file.path(project_dir, "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

outcome <- "life_dissatisfaction"
ordered_treatment <- "social_media_weekday"
binary_treatment <- "high_social_media_weekday"
repetitions <- as.integer(
  Sys.getenv("DML_ROBUSTNESS_REPETITIONS", unset = "100")
)
reuse_repeated <- identical(
  Sys.getenv("DML_REUSE_REPEATED", unset = "0"),
  "1"
)

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

plr_learner_order <- c(
  "Elastic net",
  "Lasso",
  "Random forest",
  "Neural net"
)
irm_learner_order <- c(
  "Elastic net",
  "Lasso",
  "Random forest"
)

data <- readRDS(data_path)

prepare_sample <- function(treatment) {
  vars <- unique(c(
    outcome,
    treatment,
    controls,
    "wave_f",
    "pidp"
  ))
  droplevels(
    data[
      complete.cases(data[, vars, drop = FALSE]),
      vars,
      drop = FALSE
    ]
  )
}

prepare_dml_data <- function(sample, treatment) {
  x_matrix <- model.matrix(
    ~ . - 1,
    data = sample[, c(controls, "wave_f"), drop = FALSE]
  )
  varying <- apply(x_matrix, 2L, function(x) length(unique(x)) > 1L)
  x_matrix <- x_matrix[, varying, drop = FALSE]

  dml_data <- data.table(
    pidp = sample$pidp,
    outcome_value = sample[[outcome]],
    treatment_value = sample[[treatment]]
  )
  setnames(
    dml_data,
    c("outcome_value", "treatment_value"),
    c(outcome, treatment)
  )
  dml_data <- cbind(dml_data, as.data.table(x_matrix))
  list(
    data = dml_data,
    x_cols = setdiff(
      names(dml_data),
      c("pidp", outcome, treatment)
    )
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
      ml_g = lrn(
        "regr.cv_glmnet",
        s = "lambda.min",
        alpha = 0.5
      ),
      ml_m = lrn(
        "classif.cv_glmnet",
        s = "lambda.min",
        alpha = 0.5,
        predict_type = "prob"
      )
    ),
    "Lasso" = list(
      ml_g = lrn(
        "regr.cv_glmnet",
        s = "lambda.min",
        alpha = 1
      ),
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

make_cluster_data <- function(prepared, treatment) {
  DoubleMLClusterData$new(
    prepared$data,
    y_col = outcome,
    d_cols = treatment,
    x_cols = prepared$x_cols,
    cluster_cols = "pidp"
  )
}

star_string <- function(p_value) {
  if (is.na(p_value)) return("")
  if (p_value < 0.01) return("***")
  if (p_value < 0.05) return("**")
  if (p_value < 0.10) return("*")
  ""
}

format_cell <- function(estimate, std_error, p_value) {
  sprintf(
    "%.3f%s (%.3f)",
    estimate,
    star_string(p_value),
    std_error
  )
}

# 1. Repeated cross-fitting -------------------------------------------------
plr_sample <- prepare_sample(ordered_treatment)
plr_prepared <- prepare_dml_data(plr_sample, ordered_treatment)
plr_learners <- make_plr_learners(length(plr_prepared$x_cols))

fit_repeated_plr <- function(i) {
  data.table::setDTthreads(1L)
  future::plan(future::sequential)
  learner_name <- names(plr_learners)[[i]]
  dml_data <- make_cluster_data(plr_prepared, ordered_treatment)
  set.seed(7100L + i)
  fit <- DoubleMLPLR$new(
    dml_data,
    plr_learners[[i]]$clone(),
    plr_learners[[i]]$clone(),
    n_folds = 5,
    n_rep = repetitions
  )
  fit$fit()
  split_estimates <- as.numeric(fit$all_coef[1, ])
  data.frame(
    learner = learner_name,
    repeated_estimate = unname(fit$coef[[ordered_treatment]]),
    repeated_se = unname(fit$se[[ordered_treatment]]),
    p_value = unname(fit$pval[[ordered_treatment]]),
    split_sd = sd(split_estimates),
    split_p05 = unname(quantile(split_estimates, 0.05)),
    split_p95 = unname(quantile(split_estimates, 0.95)),
    split_min = min(split_estimates),
    split_max = max(split_estimates),
    n_rep = repetitions,
    n = nrow(plr_sample),
    stringsAsFactors = FALSE
  )
}

repeated_output_path <- file.path(
  output_dir,
  "robustness_repeated_cross_fitting_numeric.csv"
)

if (reuse_repeated && file.exists(repeated_output_path)) {
  repeated_results <- fread(repeated_output_path)
} else {
  repeated_results <- rbindlist(
    parallel::mclapply(
      seq_along(plr_learners),
      fit_repeated_plr,
      mc.cores = min(4L, length(plr_learners)),
      mc.preschedule = FALSE,
      mc.set.seed = FALSE
    )
  )

  single_results <- fread(
    file.path(
      project_dir,
      "tables",
      "chapter5_core",
      "table_dml_plr_learner_comparison_numeric.csv"
    )
  )[
    outcome == "life_dissatisfaction",
    .(
      learner,
      single_estimate = estimate,
      single_se = std_error
    )
  ]

  repeated_results <- merge(
    repeated_results,
    single_results,
    by = "learner",
    all.x = TRUE
  )
  repeated_results[, learner := factor(
    learner,
    levels = plr_learner_order
  )]
  setorder(repeated_results, learner)
  repeated_results[, learner := as.character(learner)]
  fwrite(repeated_results, repeated_output_path)
}

repeat_lines <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Repeated cross-fitting robustness for life dissatisfaction}",
  "\\label{tab:robustness_repeated_cross_fitting}",
  "\\small",
  "\\setlength{\\tabcolsep}{5pt}",
  "\\begin{tabular}{lcccc}",
  "\\hline",
  paste0(
    "Learner & Single split & Repeated estimate & Repeated SE & ",
    "5th--95th percentile"
  ),
  "\\\\",
  "\\hline",
  unlist(lapply(seq_len(nrow(repeated_results)), function(i) {
    row <- repeated_results[i]
    c(
      paste(
        row$learner,
        sprintf("%.3f", row$single_estimate),
        sprintf("%.3f", row$repeated_estimate),
        sprintf("%.3f", row$repeated_se),
        sprintf("[%.3f, %.3f]", row$split_p05, row$split_p95),
        sep = " & "
      ),
      "\\\\"
    )
  })),
  "\\hline",
  "\\end{tabular}",
  "\\begin{minipage}{0.96\\textwidth}",
  paste0(
    "\\footnotesize Notes: Repeated estimates aggregate ",
    repetitions,
    " independent five-fold respondent-level partitions by their median. ",
    "The final column describes variation in the treatment estimate across ",
    "the individual sample splits; it is not a confidence interval."
  ),
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(
  repeat_lines,
  file.path(tables_dir, "table_robustness_repeated_cross_fitting.tex")
)

# 2. IRM overlap and propensity truncation ---------------------------------
irm_sample <- prepare_sample(binary_treatment)
irm_prepared <- prepare_dml_data(irm_sample, binary_treatment)
irm_learners <- make_irm_learners(length(irm_prepared$x_cols))
trimming_thresholds <- c(1e-12, 0.01, 0.025, 0.05, 0.10)

fit_irm_diagnostics <- function(i) {
  learner_name <- names(irm_learners)[[i]]
  threshold_rows <- list()
  propensity_data <- NULL

  for (j in seq_along(trimming_thresholds)) {
    threshold <- trimming_thresholds[[j]]
    dml_data <- make_cluster_data(irm_prepared, binary_treatment)
    # Reuse the same fold assignment and stochastic learner seed so that
    # differences across rows isolate the propensity truncation threshold.
    set.seed(500L + i)
    fit <- DoubleMLIRM$new(
      dml_data,
      irm_learners[[i]]$ml_g$clone(),
      irm_learners[[i]]$ml_m$clone(),
      n_folds = 5,
      score = "ATE",
      trimming_rule = "truncate",
      trimming_threshold = threshold
    )
    fit$fit(store_predictions = j == 1L)

    threshold_rows[[j]] <- data.frame(
      learner = learner_name,
      threshold = threshold,
      estimate = unname(fit$coef[[binary_treatment]]),
      std_error = unname(fit$se[[binary_treatment]]),
      p_value = unname(fit$pval[[binary_treatment]]),
      n = nrow(irm_sample),
      stringsAsFactors = FALSE
    )

    if (j == 1L) {
      propensity_data <- data.frame(
        learner = learner_name,
        treatment = irm_sample[[binary_treatment]],
        propensity = as.numeric(fit$predictions$ml_m),
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    trimming = rbindlist(threshold_rows),
    propensity = propensity_data
  )
}

irm_diagnostics <- parallel::mclapply(
  seq_along(irm_learners),
  fit_irm_diagnostics,
  mc.cores = min(3L, length(irm_learners)),
  mc.preschedule = FALSE,
  mc.set.seed = FALSE
)
trimming_results <- rbindlist(lapply(
  irm_diagnostics,
  `[[`,
  "trimming"
))
propensity_results <- rbindlist(lapply(
  irm_diagnostics,
  `[[`,
  "propensity"
))

overlap_summary <- propensity_results[
  ,
  .(
    minimum = min(propensity),
    p01 = quantile(propensity, 0.01),
    p05 = quantile(propensity, 0.05),
    median = median(propensity),
    p95 = quantile(propensity, 0.95),
    p99 = quantile(propensity, 0.99),
    maximum = max(propensity),
    share_below_001 = mean(propensity < 0.01),
    share_below_005 = mean(propensity < 0.05),
    share_above_095 = mean(propensity > 0.95)
  ),
  by = learner
]

fwrite(
  trimming_results,
  file.path(output_dir, "robustness_irm_trimming_numeric.csv")
)
fwrite(
  overlap_summary,
  file.path(output_dir, "robustness_irm_overlap_numeric.csv")
)
fwrite(
  propensity_results,
  file.path(output_dir, "robustness_irm_propensities.csv")
)

trimming_results[, threshold_label := fifelse(
  threshold < 1e-6,
  "Default",
  paste0(formatC(100 * threshold, format = "fg"), "\\%")
)]
trimming_results[, cell := mapply(
  format_cell,
  estimate,
  std_error,
  p_value
)]

trimming_lines <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{IRM propensity-truncation robustness for life dissatisfaction}",
  "\\label{tab:robustness_irm_trimming}",
  "\\small",
  "\\setlength{\\tabcolsep}{7pt}",
  "\\begin{tabular}{lccc}",
  "\\hline",
  "Propensity threshold & Elastic net & Lasso & Random forest",
  "\\\\",
  "\\hline",
  unlist(lapply(
    c("Default", "1\\%", "2.5\\%", "5\\%", "10\\%"),
    function(threshold_value) {
      row <- trimming_results[threshold_label == threshold_value]
      values <- vapply(irm_learner_order, function(learner_name) {
        value <- row[learner == learner_name, cell]
        if (length(value) == 0L) "--" else value[[1]]
      }, character(1))
      c(
        paste(c(threshold_value, values), collapse = " & "),
        "\\\\"
      )
    }
  )),
  "\\hline",
  "\\end{tabular}",
  "\\begin{minipage}{0.96\\textwidth}",
  paste0(
    "\\footnotesize Notes: Estimated propensity scores are truncated to ",
    "the interval defined by each threshold before evaluating the ATE score. ",
    "The default threshold is $10^{-12}$ and reproduces the main IRM folds; ",
    "all rows within a learner use the same fold partition. ",
    "Cells report IRM-DML estimates with respondent-clustered standard errors ",
    "in parentheses. $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$."
  ),
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(
  trimming_lines,
  file.path(tables_dir, "table_robustness_irm_trimming.tex")
)

overlap_plot <- ggplot(
  propensity_results,
  aes(
    x = propensity,
    color = factor(
      treatment,
      levels = c(0, 1),
      labels = c("Lower use", "High use")
    ),
    fill = factor(
      treatment,
      levels = c(0, 1),
      labels = c("Lower use", "High use")
    )
  )
) +
  geom_density(alpha = 0.24, linewidth = 0.8) +
  geom_vline(
    xintercept = c(0.05, 0.95),
    linetype = "dashed",
    color = "#7F6B91",
    linewidth = 0.45
  ) +
  facet_wrap(~ learner, ncol = 1) +
  scale_color_manual(values = c(
    "Lower use" = "#8F82D8",
    "High use" = "#DE8FB8"
  )) +
  scale_fill_manual(values = c(
    "Lower use" = "#DCD5F8",
    "High use" = "#F8D7E3"
  )) +
  labs(
    x = "Cross-fitted propensity score",
    y = "Density",
    color = "Observed treatment",
    fill = "Observed treatment"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#EEE9F1", linewidth = 0.35),
    strip.text = element_text(face = "bold", color = "#403A47"),
    axis.text = element_text(color = "#5A5165"),
    axis.title = element_text(color = "#403A47"),
    legend.title = element_text(color = "#403A47"),
    legend.text = element_text(color = "#5A5165"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  file.path(figures_dir, "robustness_irm_overlap_life.png"),
  overlap_plot,
  width = 7.2,
  height = 8.0,
  dpi = 320
)

# 3. PLR omitted-variable robustness values -------------------------------
fit_plr_sensitivity <- function(i) {
  learner_name <- names(plr_learners)[[i]]
  dml_data <- make_cluster_data(plr_prepared, ordered_treatment)
  set.seed(9100L + i)
  fit <- DoubleMLPLR$new(
    dml_data,
    plr_learners[[i]]$clone(),
    plr_learners[[i]]$clone(),
    n_folds = 5
  )
  fit$fit(store_predictions = TRUE)

  beta <- unname(fit$coef[[ordered_treatment]])
  d <- plr_sample[[ordered_treatment]]
  y <- plr_sample[[outcome]]
  l_hat <- as.numeric(fit$predictions$ml_l)
  m_hat <- as.numeric(fit$predictions$ml_m)
  v_hat <- d - m_hat
  u_hat <- y - l_hat - beta * v_hat
  sigma2 <- mean(u_hat^2)
  nu2 <- 1 / mean(v_hat^2)
  c_value <- beta^2 / (sigma2 * nu2)
  robustness_value <- (
    sqrt(c_value^2 + 4 * c_value) - c_value
  ) / 2

  scenario_strength <- 0.05
  bias_bound <- sqrt(
    sigma2 *
      nu2 *
      scenario_strength^2 /
      (1 - scenario_strength)
  )

  data.frame(
    learner = learner_name,
    estimate = beta,
    std_error = unname(fit$se[[ordered_treatment]]),
    sigma2 = sigma2,
    nu2 = nu2,
    robustness_value = robustness_value,
    lower_bound_005 = beta - bias_bound,
    upper_bound_005 = beta + bias_bound,
    n = nrow(plr_sample),
    stringsAsFactors = FALSE
  )
}

sensitivity_results <- rbindlist(
  parallel::mclapply(
    seq_along(plr_learners),
    fit_plr_sensitivity,
    mc.cores = min(4L, length(plr_learners)),
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
)
sensitivity_results[, learner := factor(
  learner,
  levels = plr_learner_order
)]
setorder(sensitivity_results, learner)
sensitivity_results[, learner := as.character(learner)]
fwrite(
  sensitivity_results,
  file.path(output_dir, "robustness_plr_sensitivity_numeric.csv")
)

sensitivity_lines <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Omitted-variable sensitivity of the PLR-DML estimates}",
  "\\label{tab:robustness_plr_sensitivity}",
  "\\small",
  "\\setlength{\\tabcolsep}{7pt}",
  "\\begin{tabular}{lccc}",
  "\\hline",
  "Learner & Estimate & Robustness value & 5\\% confounding bound",
  "\\\\",
  "\\hline",
  unlist(lapply(seq_len(nrow(sensitivity_results)), function(i) {
    row <- sensitivity_results[i]
    c(
      paste(
        row$learner,
        sprintf("%.3f", row$estimate),
        sprintf("%.2f\\%%", 100 * row$robustness_value),
        sprintf(
          "[%.3f, %.3f]",
          row$lower_bound_005,
          row$upper_bound_005
        ),
        sep = " & "
      ),
      "\\\\"
    )
  })),
  "\\hline",
  "\\end{tabular}",
  "\\begin{minipage}{0.96\\textwidth}",
  paste0(
    "\\footnotesize Notes: The robustness value is the common nonparametric ",
    "partial $R^2$ strength with treatment and outcome required for an ",
    "adversarial omitted confounder to move the point estimate to zero. ",
    "The final column sets both sensitivity parameters to 5\\% and uses ",
    "$\\rho=1$. These are point-identification bounds and do not remove the ",
    "need for the unconfoundedness assumption."
  ),
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(
  sensitivity_lines,
  file.path(tables_dir, "table_robustness_plr_sensitivity.tex")
)

cat("Wrote DML robustness diagnostics.\n")
print(repeated_results)
print(overlap_summary)
print(trimming_results)
print(sensitivity_results)
