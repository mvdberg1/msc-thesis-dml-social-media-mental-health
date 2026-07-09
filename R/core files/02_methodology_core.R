# Chapter 4
#
# Summary:
# 1. This script rebuilds the main Chapter 4 methodology figures from a compact
#    baseline sample based on UKHLS youth waves l--o.
# 2. It shows how the DML setup is constructed: the causal diagram, learner
#    workflow, nuisance tuning, orthogonalization idea, and bias simulation.
# 3. All outputs are shown interactively in the console and plot pane, so the
#    methodology can be explained step by step without exporting thesis files.

library(haven)
library(ggplot2)
library(glmnet)
library(rpart)

project_dir <- normalizePath(
  Sys.getenv("THESIS_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
root <- Sys.getenv(
  "UKHLS_RAW_DIR",
  unset = file.path(
    project_dir,
    "UKDA-6614-stata",
    "stata",
    "stata14_se",
    "ukhls"
  )
)
miss_codes <- c(-1, -2, -3, -4, -5, -6, -7, -8, -9)

show_table <- function(title, x) {
  cat("\n", title, "\n", sep = "")
  if (inherits(x, "data.frame")) {
    print(as.data.frame(x), row.names = FALSE)
  } else {
    print(x)
  }
}

show_plot <- function(title, plot_obj) {
  cat("\n", title, "\n", sep = "")
  print(plot_obj)
}

read_wave <- function(w) {
  ff <- file.path(root, paste0(w, "_youth.dta"))
  vars <- c(
    "pidp",
    paste0(
      w,
      c(
        "_ypnetcht", "_yphlf", "_age_dv", "_sex_dv", "_ethn_dv", "_country",
        "_gor_dv", "_urban_dv", "_ypnpal", "_ypeatlivu", "_yptvvidhrs",
        "_ypdevice1", "_ypdevice3", "_ypdevice5", "_ypdevice6"
      )
    )
  )
  x <- read_dta(ff, col_select = any_of(vars))
  names(x) <- sub(paste0("^", w, "_"), "", names(x))
  x$wave <- match(w, c("l", "m", "n", "o"))

  for (nm in setdiff(names(x), c("pidp", "wave"))) {
    x[[nm]] <- as.numeric(x[[nm]])
    x[[nm]][x[[nm]] %in% miss_codes] <- NA_real_
  }

  # Higher values indicate more life dissatisfaction in the methodology plots.
  x$life_distress <- ifelse(is.na(x$yphlf), NA_real_, 8 - x$yphlf)
  x
}

build_baseline_sample <- function() {
  df <- do.call(rbind, lapply(c("l", "m", "n", "o"), read_wave))
  keep <- c(
    "pidp", "wave", "ypnetcht", "life_distress", "age_dv", "sex_dv", "ethn_dv",
    "country", "gor_dv", "urban_dv", "ypnpal", "ypeatlivu", "yptvvidhrs",
    "ypdevice1", "ypdevice3", "ypdevice5", "ypdevice6"
  )
  df[complete.cases(df[, keep]), keep]
}

baseline_formula <- ~ wave + age_dv + sex_dv + ethn_dv + country + gor_dv +
  urban_dv + ypnpal + ypeatlivu + yptvvidhrs + ypdevice1 + ypdevice3 +
  ypdevice5 + ypdevice6

fit_tree <- function(formula, data, cp = 0.005, minsplit = 35, maxdepth = 4) {
  rpart(
    formula,
    data = data,
    method = "anova",
    control = rpart.control(cp = cp, minsplit = minsplit, maxdepth = maxdepth)
  )
}

make_causal_diagram <- function() {
  nodes <- data.frame(
    x = c(0.9, 3.4, 3.4, 5.9),
    y = c(1.0, 1.0, 0.0, 1.0),
    label = c("Y", "D", "X", "V"),
    fill = c("#d9c8f0", "#f6c8d8", "#cbb7e8", "#f2d9ea")
  )

  edges <- data.frame(
    x = c(3.08, 5.56, 3.04, 3.40),
    y = c(1.00, 1.00, 0.14, 0.24),
    xend = c(1.20, 3.72, 1.20, 3.40),
    yend = c(1.00, 1.00, 0.84, 0.76)
  )

  ggplot() +
    geom_segment(
      data = edges,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.8,
      color = "#61546f",
      arrow = grid::arrow(length = grid::unit(0.16, "inches"), type = "closed")
    ) +
    geom_point(
      data = nodes,
      aes(x = x, y = y, fill = fill),
      shape = 21,
      size = 20,
      stroke = 1.1,
      color = "#7f6b91",
      show.legend = FALSE
    ) +
    geom_text(
      data = nodes,
      aes(x = x, y = y, label = label),
      size = 7.2,
      color = "#33263d"
    ) +
    annotate(
      "text",
      x = 3.4,
      y = -0.48,
      label = "Causal diagram",
      size = 7,
      color = "#5a5165"
    ) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0.2, 6.6), ylim = c(-0.7, 1.45), expand = FALSE) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

make_learner_workflow_plot <- function() {
  boxes <- data.frame(
    xmin = c(0.3, 3.0, 5.7, 8.4),
    xmax = c(2.5, 5.2, 7.9, 10.6),
    ymin = 0.45,
    ymax = 2.25,
    fill = c("#d9edf0", "#e5dbf5", "#fae7cf", "#e5f2d8")
  )

  labels <- data.frame(
    x = c(1.4, 4.1, 6.8, 9.5),
    y = 1.40,
    label = c(
      "1. Build sample\nPooled Waves l-o\nChoose Y, D and X",
      "2. Split by pidp\n5 folds, one partition\nAll waves stay together",
      "3. Fit nuisance models\nUse four training folds\nLearner-specific settings",
      "4. Predict hold-out fold\nStack out-of-sample fits\nEstimate orthogonal score"
    )
  )

  arrows <- data.frame(
    x = c(2.5, 5.2, 7.9),
    y = 1.35,
    xend = c(3.0, 5.7, 8.4),
    yend = 1.35
  )

  fold_blocks <- data.frame(
    xmin = c(3.30, 3.62, 3.94, 4.26, 4.58),
    xmax = c(3.56, 3.88, 4.20, 4.52, 4.84),
    ymin = 0.68,
    ymax = 0.98,
    fill = c("#b8d9df", "#cdb7e8", "#f1bfd2", "#cfe7b8", "#f6d5ad")
  )

  ggplot() +
    geom_rect(
      data = boxes,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      color = "#7f6b91",
      linewidth = 0.9,
      alpha = 0.95,
      show.legend = FALSE
    ) +
    geom_segment(
      data = arrows,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.9,
      color = "#8d85a5",
      lineend = "round",
      arrow = grid::arrow(length = grid::unit(0.15, "inches"), type = "closed")
    ) +
    geom_rect(
      data = fold_blocks,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      color = "#ffffff",
      linewidth = 0.3,
      show.legend = FALSE
    ) +
    geom_text(
      data = labels,
      aes(x = x, y = y, label = label),
      size = 3.7,
      lineheight = 1.08,
      color = "#2f2940"
    ) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0, 10.9), ylim = c(0.2, 2.55), expand = FALSE) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

make_glmnet_tuning_plot <- function(dat) {
  tuning_df <- dat
  factor_vars <- c(
    "wave", "sex_dv", "ethn_dv", "country", "gor_dv", "urban_dv",
    "ypnpal", "ypeatlivu", "ypdevice1", "ypdevice3", "ypdevice5", "ypdevice6"
  )

  for (nm in factor_vars) {
    tuning_df[[nm]] <- factor(tuning_df[[nm]])
  }

  x_mat <- model.matrix(
    ~ wave + age_dv + sex_dv + ethn_dv + country + gor_dv + urban_dv +
      ypnpal + ypeatlivu + yptvvidhrs + ypdevice1 + ypdevice3 +
      ypdevice5 + ypdevice6,
    data = tuning_df
  )[, -1]

  cv_fit <- cv.glmnet(
    x = x_mat,
    y = tuning_df$ypnetcht,
    family = "gaussian",
    nfolds = 5
  )

  plot_df <- data.frame(
    log_lambda = log(cv_fit$lambda),
    cvm = cv_fit$cvm,
    cvup = cv_fit$cvup,
    cvlo = cv_fit$cvlo
  )

  ref_df <- data.frame(
    log_lambda = log(c(cv_fit$lambda.min, cv_fit$lambda.1se)),
    label = c("lambda.min", "lambda.1se")
  )

  ggplot(plot_df, aes(x = log_lambda, y = cvm)) +
    geom_ribbon(aes(ymin = cvlo, ymax = cvup), fill = "#ead7ef", alpha = 0.5) +
    geom_line(color = "#7e6aa8", linewidth = 1.0) +
    geom_vline(
      data = ref_df,
      aes(xintercept = log_lambda, linetype = label),
      color = "#d27b95",
      linewidth = 0.8,
      show.legend = TRUE
    ) +
    scale_linetype_manual(
      values = c("lambda.min" = "solid", "lambda.1se" = "dashed"),
      name = NULL
    ) +
    labs(
      x = expression(log(lambda)),
      y = "Cross-validated mean squared error"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

make_orthogonalization_plot <- function(dat) {
  dat$y_obs <- dat$life_distress

  fit_outcome <- fit_tree(update(baseline_formula, y_obs ~ .), dat, cp = 0.003, minsplit = 25, maxdepth = 5)
  fit_treat <- fit_tree(update(baseline_formula, ypnetcht ~ .), dat, cp = 0.003, minsplit = 25, maxdepth = 5)

  dat$y_resid <- dat$y_obs - predict(fit_outcome, newdata = dat)
  dat$d_resid <- dat$ypnetcht - predict(fit_treat, newdata = dat)

  q <- quantile(dat$d_resid, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
  q[duplicated(q)] <- q[duplicated(q)] + 1e-8
  dat$bin <- cut(dat$d_resid, breaks = q, include.lowest = TRUE)

  raw_df <- aggregate(life_distress ~ ypnetcht, data = dat, FUN = mean)
  names(raw_df) <- c("x", "y")
  raw_df$panel <- "Observed association"

  resid_df <- aggregate(cbind(y_resid, d_resid) ~ bin, data = dat, FUN = mean)
  resid_df <- data.frame(x = resid_df$d_resid, y = resid_df$y_resid)
  resid_df$panel <- "After orthogonalization"

  plot_df <- rbind(raw_df, resid_df)
  plot_df$panel <- factor(
    plot_df$panel,
    levels = c("Observed association", "After orthogonalization")
  )

  ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(size = 3.0, color = "#8578bc") +
    geom_smooth(method = "lm", se = FALSE, color = "#de98ab", linewidth = 1.1) +
    facet_wrap(~panel, scales = "free_x") +
    labs(x = NULL, y = "Mean life dissatisfaction") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 13),
      panel.grid.minor = element_blank()
    )
}

calibrate_empirical_dgp <- function(dat, beta0 = 0.08) {
  fit_m0 <- fit_tree(update(baseline_formula, ypnetcht ~ .), dat, cp = 0.005, minsplit = 35, maxdepth = 4)
  m0_hat <- predict(fit_m0, newdata = dat)

  age_z <- as.numeric(scale(dat$age_dv))
  tv_z <- as.numeric(scale(dat$yptvvidhrs))
  group_indicator <- as.numeric(dat$sex_dv == 2 & dat$country == 2)

  g0 <- 0.45 * m0_hat + 0.30 * age_z^2 + 0.20 * sin(tv_z) + 0.20 * group_indicator

  list(beta0 = beta0, g0 = g0, sigma_u = 0.8)
}

simulate_outcome <- function(dat, dgp, seed) {
  set.seed(seed)
  dgp$beta0 * dat$ypnetcht + dgp$g0 + rnorm(nrow(dat), sd = dgp$sigma_u)
}

estimate_prediction_ml <- function(dat, y_sim, seed) {
  set.seed(seed)
  ids <- unique(dat$pidp)
  train_ids <- sample(ids, ceiling(length(ids) / 2))
  train_rows <- dat$pidp %in% train_ids

  train <- dat[train_rows, ]
  test <- dat[!train_rows, ]
  train$y_sim <- y_sim[train_rows]

  fit_l <- fit_tree(update(baseline_formula, y_sim ~ .), train, cp = 0.01, minsplit = 50, maxdepth = 3)
  l_hat <- predict(fit_l, newdata = test)

  y_test <- y_sim[!train_rows]
  d_test <- test$ypnetcht

  beta_hat <- sum(d_test * (y_test - l_hat)) / sum(d_test^2)
  psi <- d_test * (y_test - l_hat - beta_hat * d_test)
  j_hat <- -mean(d_test^2)
  se_hat <- sqrt(mean(psi^2) / (length(d_test) * j_hat^2))

  c(beta = beta_hat, se = se_hat)
}

estimate_dml_partialling <- function(dat, y_sim, seed, k_folds = 5) {
  set.seed(seed)
  ids <- unique(dat$pidp)
  fold_id <- sample(rep(seq_len(k_folds), length.out = length(ids)))
  names(fold_id) <- ids

  l_hat <- rep(NA_real_, nrow(dat))
  m_hat <- rep(NA_real_, nrow(dat))

  for (k in seq_len(k_folds)) {
    test_rows <- fold_id[as.character(dat$pidp)] == k
    train <- dat[!test_rows, ]
    test <- dat[test_rows, ]

    train$y_sim <- y_sim[!test_rows]

    fit_l <- fit_tree(update(baseline_formula, y_sim ~ .), train, cp = 0.005, minsplit = 35, maxdepth = 4)
    fit_m <- fit_tree(update(baseline_formula, ypnetcht ~ .), train, cp = 0.005, minsplit = 35, maxdepth = 4)

    l_hat[test_rows] <- predict(fit_l, newdata = test)
    m_hat[test_rows] <- predict(fit_m, newdata = test)
  }

  v_hat <- dat$ypnetcht - m_hat
  w_hat <- y_sim - l_hat

  beta_hat <- sum(v_hat * w_hat) / sum(v_hat^2)
  psi <- (w_hat - beta_hat * v_hat) * v_hat
  j_hat <- -mean(v_hat^2)
  se_hat <- sqrt(mean(psi^2) / (length(v_hat) * j_hat^2))

  c(beta = beta_hat, se = se_hat)
}

make_bias_replications <- function(dat, n_rep = 150, beta0 = 0.08) {
  dgp <- calibrate_empirical_dgp(dat, beta0 = beta0)
  res <- data.frame(
    rep = seq_len(n_rep),
    beta_pred = NA_real_,
    se_pred = NA_real_,
    beta_dml = NA_real_,
    se_dml = NA_real_
  )

  for (b in seq_len(n_rep)) {
    y_sim <- simulate_outcome(dat, dgp, seed = 1000 + b)
    pred_est <- estimate_prediction_ml(dat, y_sim, seed = 2000 + b)
    dml_est <- estimate_dml_partialling(dat, y_sim, seed = 3000 + b, k_folds = 5)

    res$beta_pred[b] <- pred_est["beta"]
    res$se_pred[b] <- pred_est["se"]
    res$beta_dml[b] <- dml_est["beta"]
    res$se_dml[b] <- dml_est["se"]
  }

  res$beta0 <- dgp$beta0
  res
}

summarise_bias_methods <- function(rep_df) {
  data.frame(
    method = c("Prediction-focused ML", "Double ML"),
    mean_beta = c(mean(rep_df$beta_pred), mean(rep_df$beta_dml)),
    empirical_sd = c(stats::sd(rep_df$beta_pred), stats::sd(rep_df$beta_dml)),
    mean_reported_se = c(mean(rep_df$se_pred), mean(rep_df$se_dml)),
    bias = c(mean(rep_df$beta_pred - rep_df$beta0), mean(rep_df$beta_dml - rep_df$beta0)),
    rmse = c(
      sqrt(mean((rep_df$beta_pred - rep_df$beta0)^2)),
      sqrt(mean((rep_df$beta_dml - rep_df$beta0)^2))
    ),
    stringsAsFactors = FALSE
  )
}

make_bias_comparison_plot <- function(rep_df) {
  plot_df <- rbind(
    data.frame(
      t_stat = (rep_df$beta_pred - rep_df$beta0) / rep_df$se_pred,
      method = "Prediction-focused ML\n(non-orthogonal)"
    ),
    data.frame(
      t_stat = (rep_df$beta_dml - rep_df$beta0) / rep_df$se_dml,
      method = "Double ML\n(orthogonal, cross-fitted)"
    )
  )
  plot_df$method <- factor(
    plot_df$method,
    levels = c(
      "Prediction-focused ML\n(non-orthogonal)",
      "Double ML\n(orthogonal, cross-fitted)"
    )
  )

  ggplot(plot_df, aes(x = t_stat)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 26,
      fill = "#f3c4d3",
      color = "#8e7dbe",
      alpha = 0.75
    ) +
    geom_vline(xintercept = 0, color = "#6d597a", linewidth = 0.7) +
    geom_function(fun = dnorm, color = "#c67aa1", linewidth = 1.0) +
    facet_wrap(~method, nrow = 1) +
    coord_cartesian(xlim = c(-15, 5)) +
    labs(
      x = "Studentized statistic",
      y = "Density"
    ) +
    theme_minimal(base_size = 12)
}

baseline_df <- build_baseline_sample()
show_table(
  "Methodology baseline sample overview",
  data.frame(
    person_wave_rows = nrow(baseline_df),
    unique_adolescents = length(unique(baseline_df$pidp)),
    waves_included = length(unique(baseline_df$wave)),
    stringsAsFactors = FALSE
  )
)
show_table("Methodology baseline sample preview (first 10 rows)", utils::head(baseline_df, 10))

show_plot("Figure: causal diagram", make_causal_diagram())
show_plot("Figure: learner workflow", make_learner_workflow_plot())
show_plot("Figure: orthogonalization", make_orthogonalization_plot(baseline_df))
show_plot(
  "Figure: elastic-net tuning for weekday social-media use",
  make_glmnet_tuning_plot(baseline_df)
)

bias_df <- make_bias_replications(baseline_df, n_rep = 150, beta0 = 0.08)
show_table("Bias simulation preview (first 10 replications)", utils::head(bias_df, 10))
show_table("Bias simulation summary", summarise_bias_methods(bias_df))
show_plot("Figure: bias comparison", make_bias_comparison_plot(bias_df))

cat("\nMethodology core outputs were printed to the console and plot pane.\n")
