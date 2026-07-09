# Design-based robustness checks for life dissatisfaction.
#
# The script evaluates:
# 1. weekday versus weekend social-media use on one common l--o sample;
# 2. rich l--o controls versus common-core controls in l--o and a--o;
# 3. leave-one-wave-out estimates for the four baseline waves.

suppressPackageStartupMessages({
  library(data.table)
  library(DoubleML)
  library(haven)
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
raw_root <- file.path(
  project_dir,
  "UKDA-6614-stata",
  "stata",
  "stata14_se",
  "ukhls"
)
baseline_path <- file.path(
  project_dir,
  "data",
  "analysis",
  "ukhls_youth_l_to_o_clean.rds"
)
all_wave_path <- file.path(
  project_dir,
  "data",
  "analysis",
  "ukhls_youth_a_to_o_common_core.rds"
)
output_dir <- file.path(project_dir, "tables", "chapter5_robustness")
tables_dir <- file.path(project_dir, "tables")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

outcome <- "life_dissatisfaction"
weekday_treatment <- "social_media_weekday"
weekend_treatment <- "social_media_weekend"

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

common_controls <- c(
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
  "yptvvidhrs"
)

learner_order <- c(
  "Fixed effects",
  "Elastic net",
  "Lasso",
  "Random forest",
  "Neural net"
)

star_string <- function(p_value) {
  if (is.na(p_value)) return("")
  if (p_value < 0.01) return("***")
  if (p_value < 0.05) return("**")
  if (p_value < 0.10) return("*")
  ""
}

format_cell <- function(estimate, std_error, p_value) {
  if (any(!is.finite(c(estimate, std_error, p_value)))) return("--")
  sprintf("%.3f%s (%.3f)", estimate, star_string(p_value), std_error)
}

escape_latex <- function(x) {
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x
}

clean_numeric <- function(x) {
  x <- as.numeric(haven::zap_labels(x))
  x[x < 0] <- NA_real_
  x
}

factor_from_codes <- function(x, labels) {
  factor(labels[as.character(x)], levels = unname(labels))
}

build_all_wave_common_core <- function() {
  all_waves <- letters[1:15]
  wave_numbers <- setNames(seq_along(all_waves), all_waves)
  raw_vars <- c(
    "ypnetcht",
    "yphlf",
    "age_dv",
    "sex_dv",
    "ethn_dv",
    "urban_dv",
    "gor_dv",
    "npns_dv",
    "ngrp_dv",
    "nnssib_dv",
    "ypnpal",
    "ypeatlivu",
    "yptvvidhrs"
  )

  pieces <- lapply(all_waves, function(wave) {
    file_path <- file.path(raw_root, paste0(wave, "_youth.dta"))
    keep <- c("pidp", paste0(wave, "_", raw_vars))
    data <- as.data.frame(
      read_dta(file_path, col_select = all_of(keep))
    )
    names(data) <- sub(paste0("^", wave, "_"), "", names(data))
    numeric_cols <- names(data)
    data[numeric_cols] <- lapply(data[numeric_cols], clean_numeric)
    data$wave <- wave
    data$wave_number <- unname(wave_numbers[[wave]])
    data
  })

  data <- rbindlist(pieces, use.names = TRUE, fill = TRUE)
  data[!(ypnetcht %in% 1:5), ypnetcht := NA_real_]
  data[!(yphlf %in% 1:7), yphlf := NA_real_]
  data[ypnpal > 100, ypnpal := NA_real_]

  data[, social_media_weekday := ypnetcht]
  data[, life_dissatisfaction := yphlf]
  data[, close_friends_log := log1p(ypnpal)]
  data[, wave_f := factor(
    wave,
    levels = all_waves,
    labels = paste("Wave", toupper(all_waves))
  )]
  data[, sex_f := factor_from_codes(
    sex_dv,
    c(`1` = "Male", `2` = "Female")
  )]
  data[, urban_f := factor_from_codes(
    urban_dv,
    c(`1` = "Urban", `2` = "Rural")
  )]
  data[, region_f := factor(gor_dv)]
  data[, ethnicity_broad_f := factor(
    fifelse(
      ethn_dv %in% 1:4,
      "White",
      fifelse(
        ethn_dv %in% 5:8,
        "Mixed",
        fifelse(
          ethn_dv %in% 9:13,
          "Asian",
          fifelse(
            ethn_dv %in% 14:16,
            "Black",
            fifelse(
              ethn_dv == 17,
              "Arab",
              fifelse(ethn_dv == 97, "Other", NA_character_)
            )
          )
        )
      )
    ),
    levels = c("White", "Mixed", "Asian", "Black", "Arab", "Other")
  )]

  saveRDS(as.data.frame(data), all_wave_path)
  as.data.frame(data)
}

load_all_wave_common_core <- function() {
  if (file.exists(all_wave_path)) {
    readRDS(all_wave_path)
  } else {
    build_all_wave_common_core()
  }
}

make_plr_learners <- function(n_features, selected = NULL) {
  learners <- list(
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
  if (!is.null(selected)) learners <- learners[selected]
  learners
}

prepare_sample <- function(
  data,
  treatment,
  controls,
  require_extra = character(0)
) {
  vars <- unique(c(
    outcome,
    treatment,
    require_extra,
    controls,
    "wave_f",
    "pidp",
    "wave"
  ))
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0L) {
    stop("Missing variables: ", paste(missing_vars, collapse = ", "))
  }
  sample <- data[
    complete.cases(data[, vars, drop = FALSE]),
    vars,
    drop = FALSE
  ]
  droplevels(sample)
}

fit_fixed_effects <- function(data, treatment, controls) {
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
  vcov_mat <- vcovHC(
    fit,
    method = "arellano",
    type = "HC1",
    cluster = "group"
  )
  test <- coeftest(fit, vcov. = vcov_mat)
  p_col <- grep("Pr\\(", colnames(test), value = TRUE)[1]
  data.frame(
    model = "Fixed effects",
    estimate = unname(test[treatment, "Estimate"]),
    std_error = unname(test[treatment, "Std. Error"]),
    p_value = unname(test[treatment, p_col]),
    n = nrow(data),
    n_clusters = length(unique(data$pidp)),
    stringsAsFactors = FALSE
  )
}

prepare_dml_data <- function(data, treatment, controls) {
  x_cols <- c(controls, "wave_f")
  x_matrix <- model.matrix(
    ~ . - 1,
    data = data[, x_cols, drop = FALSE]
  )
  varying <- apply(x_matrix, 2L, function(x) length(unique(x)) > 1L)
  x_matrix <- x_matrix[, varying, drop = FALSE]

  dml_data <- data.table(
    pidp = data$pidp,
    outcome_value = data[[outcome]],
    treatment_value = data[[treatment]]
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

fit_plr <- function(
  data,
  treatment,
  controls,
  selected_learners = NULL,
  seed_offset = 0L
) {
  prepared <- prepare_dml_data(data, treatment, controls)
  dml_data <- DoubleMLClusterData$new(
    prepared$data,
    y_col = outcome,
    d_cols = treatment,
    x_cols = prepared$x_cols,
    cluster_cols = "pidp"
  )
  learners <- make_plr_learners(
    length(prepared$x_cols),
    selected = selected_learners
  )

  rbindlist(lapply(seq_along(learners), function(i) {
    learner_name <- names(learners)[[i]]
    set.seed(100L + seed_offset + i)
    fit <- DoubleMLPLR$new(
      dml_data,
      learners[[i]]$clone(),
      learners[[i]]$clone(),
      n_folds = 5
    )
    fit$fit()
    data.frame(
      model = learner_name,
      estimate = unname(fit$coef[[treatment]]),
      std_error = unname(fit$se[[treatment]]),
      p_value = unname(fit$pval[[treatment]]),
      n = nrow(data),
      n_clusters = length(unique(data$pidp)),
      stringsAsFactors = FALSE
    )
  }))
}

fit_all_models <- function(
  data,
  treatment,
  controls,
  selected_learners = NULL,
  seed_offset = 0L
) {
  rbindlist(list(
    fit_fixed_effects(data, treatment, controls),
    fit_plr(
      data,
      treatment,
      controls,
      selected_learners = selected_learners,
      seed_offset = seed_offset
    )
  ))
}

write_estimate_table <- function(
  results,
  row_var,
  row_order,
  row_labels,
  caption,
  label,
  notes,
  path,
  models = learner_order
) {
  results <- as.data.table(results)
  results[, cell := mapply(
    format_cell,
    estimate,
    std_error,
    p_value
  )]

  rows <- lapply(row_order, function(row_value) {
    row_data <- results[get(row_var) == row_value]
    cells <- vapply(models, function(model_name) {
      value <- row_data[model == model_name, cell]
      if (length(value) == 0L) "--" else value[[1]]
    }, character(1))
    n_value <- unique(row_data$n)
    c(
      row_labels[[row_value]],
      cells,
      format(n_value[[1]], big.mark = ",")
    )
  })

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\small",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\resizebox{\\textwidth}{!}{%",
    paste0(
      "\\begin{tabular}{l",
      paste(rep("c", length(models)), collapse = ""),
      "r}"
    ),
    "\\hline",
    paste(
      c("Specification", models, "$N$"),
      collapse = " & "
    ),
    "\\\\",
    "\\hline",
    unlist(lapply(rows, function(row) {
      c(paste(row, collapse = " & "), "\\\\")
    })),
    "\\hline",
    "\\end{tabular}",
    "}",
    "\\begin{minipage}{0.97\\textwidth}",
    paste0("\\footnotesize Notes: ", notes),
    "\\end{minipage}",
    "\\end{table}"
  )
  writeLines(lines, path)
}

baseline <- readRDS(baseline_path)
all_wave <- load_all_wave_common_core()

# 1. Weekday versus weekend use on the same complete-case sample -----------
weekend_sample <- prepare_sample(
  baseline,
  weekday_treatment,
  rich_controls,
  require_extra = weekend_treatment
)

weekday_same_sample <- fit_all_models(
  weekend_sample,
  weekday_treatment,
  rich_controls,
  seed_offset = 1000L
)
weekday_same_sample[, specification := "weekday"]

weekend_same_sample <- fit_all_models(
  weekend_sample,
  weekend_treatment,
  rich_controls,
  # Match the weekday fold partition learner by learner. The samples have
  # identical row order, so differences isolate the treatment definition.
  seed_offset = 1000L
)
weekend_same_sample[, specification := "weekend"]

weekend_results <- rbindlist(list(
  weekday_same_sample,
  weekend_same_sample
))
fwrite(
  weekend_results,
  file.path(output_dir, "robustness_weekday_weekend_numeric.csv")
)
write_estimate_table(
  weekend_results,
  row_var = "specification",
  row_order = c("weekday", "weekend"),
  row_labels = c(
    weekday = "Weekday use (common sample)",
    weekend = "Weekend use (common sample)"
  ),
  caption = paste0(
    "Weekday and weekend social-media use estimates for life ",
    "dissatisfaction"
  ),
  label = "tab:robustness_weekday_weekend_common_sample",
  notes = paste0(
    "Both rows use the same complete-case sample requiring non-missing ",
    "weekday and weekend use and the rich baseline control set. This ",
    "common sample contains 25 fewer observations than the weekday-only ",
    "sample in Table~\\ref{tab:dml_life_methods}. Within each learner, ",
    "both rows use the same cross-fitting partition. Cells report ",
    "coefficients with respondent-clustered ",
    "standard errors in parentheses. The fixed-effects model includes ",
    "respondent and wave effects; the DML columns use five-fold ",
    "respondent-level cross-fitting. ",
    "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$."
  ),
  path = file.path(tables_dir, "table_robustness_weekday_weekend_common_sample.tex")
)

# 2. Extended wave coverage with common-core controls ----------------------
panel_results <- fread(
  file.path(
    project_dir,
    "tables",
    "chapter5_core",
    "table_panel_benchmarks_numeric.csv"
  )
)
plr_results <- fread(
  file.path(
    project_dir,
    "tables",
    "chapter5_core",
    "table_dml_plr_learner_comparison_numeric.csv"
  )
)

baseline_reported <- rbindlist(list(
  panel_results[
    outcome == "life_dissatisfaction" &
      model == "Fixed effects + wave FE",
    .(
      model = "Fixed effects",
      estimate,
      std_error,
      p_value,
      n,
      n_clusters = NA_integer_
    )
  ],
  plr_results[
    outcome == "life_dissatisfaction",
    .(
      model = learner,
      estimate,
      std_error,
      p_value,
      n,
      n_clusters = NA_integer_
    )
  ]
))
baseline_reported[, specification := "lo_rich"]

lo_common_sample <- prepare_sample(
  baseline,
  weekday_treatment,
  common_controls
)
lo_common <- fit_all_models(
  lo_common_sample,
  weekday_treatment,
  common_controls,
  seed_offset = 3000L
)
lo_common[, specification := "lo_common"]

ao_common_sample <- prepare_sample(
  all_wave,
  weekday_treatment,
  common_controls
)
ao_common <- fit_all_models(
  ao_common_sample,
  weekday_treatment,
  common_controls,
  seed_offset = 4000L
)
ao_common[, specification := "ao_common"]

extended_results <- rbindlist(list(
  baseline_reported,
  lo_common,
  ao_common
), fill = TRUE)
fwrite(
  extended_results,
  file.path(output_dir, "robustness_extended_waves_numeric.csv")
)
write_estimate_table(
  extended_results,
  row_var = "specification",
  row_order = c("lo_rich", "lo_common", "ao_common"),
  row_labels = c(
    lo_rich = "Waves l--o, rich controls",
    lo_common = "Waves l--o, common-core controls",
    ao_common = "Waves a--o, common-core controls"
  ),
  caption = paste0(
    "Wave coverage and control-set robustness for life dissatisfaction"
  ),
  label = "tab:robustness_extended_waves",
  notes = paste0(
    "The first row reproduces the reported baseline estimates. The second ",
    "row uses the same four waves but removes device controls, isolating ",
    "the effect of the common-core specification. The third row extends ",
    "the common-core specification to Waves a--o. Cells report ",
    "coefficients with respondent-clustered standard errors in parentheses. ",
    "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$."
  ),
  path = file.path(
    tables_dir,
    "table_robustness_extended_waves.tex"
  )
)

# 3. Leave-one-wave-out estimates ------------------------------------------
leave_out_results <- rbindlist(lapply(
  c("l", "m", "n", "o"),
  function(omitted_wave) {
    reduced <- baseline[baseline$wave != omitted_wave, , drop = FALSE]
    sample <- prepare_sample(
      reduced,
      weekday_treatment,
      rich_controls
    )
    estimates <- fit_all_models(
      sample,
      weekday_treatment,
      rich_controls,
      selected_learners = c("Elastic net", "Random forest"),
      seed_offset = 5000L + match(omitted_wave, c("l", "m", "n", "o")) * 100L
    )
    estimates[, specification := omitted_wave]
    estimates
  }
))

fwrite(
  leave_out_results,
  file.path(output_dir, "robustness_leave_one_wave_out_numeric.csv")
)
write_estimate_table(
  leave_out_results,
  row_var = "specification",
  row_order = c("l", "m", "n", "o"),
  row_labels = c(
    l = "Exclude Wave l",
    m = "Exclude Wave m",
    n = "Exclude Wave n",
    o = "Exclude Wave o"
  ),
  caption = "Leave-one-wave-out estimates for life dissatisfaction",
  label = "tab:robustness_leave_one_wave_out",
  notes = paste0(
    "Each row removes one baseline wave and rebuilds the complete-case ",
    "sample. Elastic net and random forest represent a penalized linear ",
    "and a flexible tree-based nuisance specification. Cells report ",
    "coefficients with respondent-clustered standard errors in parentheses. ",
    "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$."
  ),
  path = file.path(
    tables_dir,
    "table_robustness_leave_one_wave_out.tex"
  ),
  models = c("Fixed effects", "Elastic net", "Random forest")
)

cat("Wrote design robustness outputs.\n")
print(weekend_results)
print(extended_results)
print(leave_out_results)
