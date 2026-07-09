library(dplyr)
library(lmtest)
library(sandwich)
library(plm)
library(nlme)
library(ggplot2)
library(gridExtra)

data_path <- "data/analysis/ukhls_youth_l_to_o_clean.rds"
tables_dir <- "tables"
figures_dir <- "figures"

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

dt <- readRDS(data_path)
dt <- dt %>%
  mutate(
    female_ind = if_else(as.character(sex_f) == "Female", 1, 0, missing = NA_real_),
    urban_ind = if_else(as.character(urban_f) == "Urban", 1, 0, missing = NA_real_)
  )

outcomes <- c(
  loneliness = "Loneliness",
  life_dissatisfaction = "Life dissatisfaction",
  schoolwork_dissatisfaction = "School-work dissatisfaction",
  school_dissatisfaction = "School dissatisfaction"
)

treatment_var <- "social_media_weekday"
main_outcome <- "life_dissatisfaction"
some_controls <- c("age_dv", "sex_f")
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
life_basic_controls <- c("age_dv", "female_ind", "urban_ind")
life_rich_controls <- c(
  "age_dv",
  "female_ind",
  "urban_ind",
  "ethnicity_broad_f",
  "region_f",
  "npns_dv",
  "ngrp_dv",
  "nnssib_dv",
  "close_friends_log",
  "ypeatlivu",
  "yptvvidhrs",
  paste0("ypdevice", 1:6)
)

escape_latex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x
}

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

format_coef_cell <- function(estimate, std_error, p_value) {
  if (any(is.na(c(estimate, std_error, p_value)))) {
    return("")
  }
  sprintf("%.3f%s (%.3f)", estimate, star_string(p_value), std_error)
}

format_coef_cell_stack <- function(estimate, std_error, p_value) {
  if (any(is.na(c(estimate, std_error, p_value)))) {
    return("")
  }
  sprintf("\\shortstack[c]{%.3f%s \\\\[-0.35ex] (%.3f)}", estimate, star_string(p_value), std_error)
}

format_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return("")
  }
  if (p_value < 0.001) {
    return("$< 0.001$")
  }
  sprintf("%.3f", p_value)
}

format_yes_no <- function(flag) {
  if (isTRUE(flag)) "Yes" else "No"
}

prepare_common_sample <- function(data, outcome, include_wave = TRUE) {
  vars_needed <- unique(c(
    outcome,
    treatment_var,
    many_controls,
    if (include_wave) "wave_f",
    "pidp",
    "wave"
  ))
  vars_needed <- intersect(vars_needed, names(data))
  out <- data[complete.cases(data[vars_needed]), vars_needed, drop = FALSE]
  droplevels(out)
}

fit_pooled_lm <- function(data, outcome, controls, wave_fe = FALSE) {
  rhs <- c(treatment_var, controls, if (wave_fe) "wave_f")
  formula <- as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
  lm(formula, data = data)
}

extract_lm_term <- function(model, data, term) {
  vcov_mat <- vcovCL(model, cluster = data$pidp, type = "HC1")
  ct <- coeftest(model, vcov. = vcov_mat)
  if (!(term %in% rownames(ct))) {
    return(c(estimate = NA_real_, se = NA_real_, p_value = NA_real_))
  }
  c(
    estimate = unname(ct[term, "Estimate"]),
    se = unname(ct[term, "Std. Error"]),
    p_value = unname(ct[term, "Pr(>|t|)"])
  )
}

write_pooled_ols_table <- function(outcome_subset, file_name, caption, label) {
  specs <- list(
    list(label = "(1)", controls = character(0), wave_fe = FALSE),
    list(label = "(2)", controls = some_controls, wave_fe = FALSE),
    list(label = "(3)", controls = many_controls, wave_fe = FALSE),
    list(label = "(4)", controls = many_controls, wave_fe = TRUE)
  )

  term_rows <- list(
    list(label = "Weekday social-media use", term = treatment_var),
    list(label = "Age", term = "age_dv"),
    list(label = "Female", term = "sex_fFemale"),
    list(label = "Log(1 + close friends)", term = "close_friends_log"),
    list(label = "Family meal frequency", term = "ypeatlivu"),
    list(label = "TV/video time", term = "yptvvidhrs")
  )

  selected_outcomes <- outcomes[outcome_subset]
  panel_letters <- LETTERS[seq_along(selected_outcomes)]
  numeric_rows <- list()

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\renewcommand{\\arraystretch}{0.88}",
    "\\begin{tabular}{lcccc}",
    "\\hline",
    "Variable & (1) & (2) & (3) & (4) \\\\",
    "\\hline"
  )

  for (i in seq_along(selected_outcomes)) {
    outcome_name <- names(selected_outcomes)[i]
    outcome_label <- selected_outcomes[[i]]
    panel_title <- paste0("Panel ", panel_letters[i], ": ", outcome_label)
    panel_data <- prepare_common_sample(dt, outcome_name, include_wave = TRUE)
    models <- lapply(specs, function(spec) {
      fit_pooled_lm(
        data = panel_data,
        outcome = outcome_name,
        controls = spec$controls,
        wave_fe = spec$wave_fe
      )
    })

    if (length(selected_outcomes) > 1) {
      lines <- c(
        lines,
        paste0("\\multicolumn{5}{l}{\\textit{", escape_latex(panel_title), "}} \\\\"),
        "\\hline"
      )
    }

    for (term_row in term_rows) {
      cells <- vapply(seq_along(models), function(j) {
        extracted <- extract_lm_term(models[[j]], panel_data, term_row$term)
        format_coef_cell_stack(extracted["estimate"], extracted["se"], extracted["p_value"])
      }, character(1))

      lines <- c(
        lines,
        paste(c(term_row$label, cells), collapse = " & "),
        "\\\\"
      )

      for (j in seq_along(models)) {
        extracted <- extract_lm_term(models[[j]], panel_data, term_row$term)
        numeric_rows[[length(numeric_rows) + 1]] <- data.frame(
          outcome = outcome_name,
          outcome_label = outcome_label,
          row_label = term_row$label,
          specification = specs[[j]]$label,
          estimate = extracted["estimate"],
          std_error = extracted["se"],
          p_value = extracted["p_value"],
          n = nobs(models[[j]]),
          stringsAsFactors = FALSE
        )
      }
    }

    adj_r2 <- vapply(models, function(model) sprintf("%.3f", summary(model)$adj.r.squared), character(1))
    n_values <- vapply(models, function(model) format(nobs(model), big.mark = ","), character(1))

    lines <- c(
      lines,
      "\\hline",
      paste(c("Adjusted $R^2$", adj_r2), collapse = " & "),
      "\\\\",
      paste(c("N", n_values), collapse = " & "),
      "\\\\",
      "\\hline"
    )
  }

  lines <- c(
    lines,
    "\\end{tabular}",
    "\\renewcommand{\\arraystretch}{1}",
    "\\begin{minipage}{0.95\\textwidth}",
    "\\footnotesize Notes: Each column reports a pooled person-wave OLS regression estimated on a common outcome-specific sample, so the sample size is held fixed within each panel. Standard errors are clustered at the respondent level and reported in parentheses below the coefficient. Column (1) is bivariate, column (2) adds age and sex, column (3) adds the rich baseline control set, and column (4) adds wave fixed effects. The rich control set adds broad ethnicity, urban residence, region, household-composition counts, close friends, family-meal frequency, television/video time, and device ownership. Significance: *** \\(p < 0.01\\), ** \\(p < 0.05\\), * \\(p < 0.10\\).",
    "\\end{minipage}",
    "\\end{table}"
  )

  writeLines(lines, file.path(tables_dir, file_name))
  write.csv(
    bind_rows(numeric_rows),
    file.path(tables_dir, sub("\\.tex$", "_numeric.csv", file_name)),
    row.names = FALSE
  )
}

write_linear_model_table <- function(outcome_name, file_name, caption, label) {
  outcome_label <- outcomes[[outcome_name]]
  vars_needed <- unique(c(
    outcome_name,
    treatment_var,
    life_rich_controls,
    "wave_f",
    "pidp",
    "wave"
  ))
  panel_data <- droplevels(dt[complete.cases(dt[vars_needed]), vars_needed, drop = FALSE])
  pdata <- pdata.frame(panel_data, index = c("pidp", "wave"), drop.index = TRUE)

  specs <- list(
    list(label = "(1)", kind = "lm", controls = character(0), wave_fe = FALSE),
    list(label = "(2)", kind = "lm", controls = life_basic_controls, wave_fe = FALSE),
    list(label = "(3)", kind = "lm", controls = life_rich_controls, wave_fe = FALSE),
    list(label = "(4)", kind = "fe", controls = life_rich_controls, wave_fe = FALSE),
    list(label = "(5)", kind = "fe", controls = life_rich_controls, wave_fe = TRUE),
    list(label = "(6)", kind = "re", controls = life_rich_controls, wave_fe = FALSE)
  )

  term_rows <- list(
    list(label = "Weekday social-media use", term = treatment_var),
    list(label = "Age", term = "age_dv"),
    list(label = "Female", term = "female_ind"),
    list(label = "Urban area", term = "urban_ind"),
    list(label = "Close friends", term = "close_friends_log"),
    list(label = "Parents/step-parents in household", term = "npns_dv"),
    list(label = "Grandparents in household", term = "ngrp_dv"),
    list(label = "Siblings/step-siblings in household", term = "nnssib_dv"),
    list(label = "Family meals", term = "ypeatlivu"),
    list(label = "Television/video hours", term = "yptvvidhrs"),
    list(label = "Smartphone", term = "ypdevice1"),
    list(label = "Tablet", term = "ypdevice3"),
    list(label = "Gaming console", term = "ypdevice5"),
    list(label = "Laptop/desktop computer", term = "ypdevice6")
  )

  fit_one <- function(spec) {
    rhs <- c(treatment_var, spec$controls, if (spec$wave_fe) "wave_f")
    formula <- as.formula(paste(outcome_name, "~", paste(rhs, collapse = " + ")))

    if (spec$kind == "lm") {
      return(lm(formula, data = panel_data))
    }
    if (spec$kind == "fe") {
      return(plm(formula, data = pdata, model = "within", effect = "individual"))
    }
    plm(
      formula,
      data = pdata,
      model = "random",
      effect = "individual",
      random.method = "swar"
    )
  }

  extract_term <- function(model, spec, term) {
    if (spec$kind == "lm") {
      return(extract_lm_term(model, panel_data, term))
    }
    extract_plm_term(model, term)
  }

  model_fit_value <- function(model, spec) {
    if (spec$kind == "lm") {
      return(sprintf("%.3f", summary(model)$r.squared))
    }
    rsq <- summary(model)$r.squared["rsq"]
    if (is.na(rsq)) {
      return("")
    }
    sprintf("%.3f", rsq)
  }

  models <- lapply(specs, fit_one)
  numeric_rows <- list()

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3.2pt}",
    "\\renewcommand{\\arraystretch}{0.9}",
    "\\begin{tabular}{p{0.34\\textwidth}cccccc}",
    "\\hline",
    "Variable & (1) & (2) & (3) & (4) & (5) & (6) \\\\",
    "\\hline"
  )

  for (term_row in term_rows) {
    cells <- vapply(seq_along(models), function(j) {
      extracted <- extract_term(models[[j]], specs[[j]], term_row$term)
      format_coef_cell_stack(extracted["estimate"], extracted["se"], extracted["p_value"])
    }, character(1))

    lines <- c(
      lines,
      paste(c(term_row$label, cells), collapse = " & "),
      "\\\\"
    )

    for (j in seq_along(models)) {
      extracted <- extract_term(models[[j]], specs[[j]], term_row$term)
      numeric_rows[[length(numeric_rows) + 1]] <- data.frame(
        outcome = outcome_name,
        outcome_label = outcome_label,
        row_label = term_row$label,
        specification = specs[[j]]$label,
        model_kind = specs[[j]]$kind,
        estimate = extracted["estimate"],
        std_error = extracted["se"],
        p_value = extracted["p_value"],
        n = nobs(models[[j]]),
        stringsAsFactors = FALSE
      )
    }
  }

  r2_values <- vapply(seq_along(models), function(j) model_fit_value(models[[j]], specs[[j]]), character(1))
  n_values <- vapply(models, function(model) format(nobs(model), big.mark = ","), character(1))
  model_labels <- c("Bivariate", "Basic", "Rich pooled", "FE", "FE + wave", "RE")

  lines <- c(
    lines,
    "\\hline",
    paste(c("Model", model_labels), collapse = " & "),
    "\\\\",
    paste(c("$R^2$", r2_values), collapse = " & "),
    "\\\\",
    paste(c("N", n_values), collapse = " & "),
    "\\\\",
    "\\hline",
    "\\end{tabular}",
    "\\renewcommand{\\arraystretch}{1}",
    "\\begin{minipage}{0.97\\textwidth}",
    "\\footnotesize Notes: All columns use the same common outcome-specific sample and report respondent-clustered standard errors in parentheses below the coefficient. Column (1) is a bivariate pooled OLS model. Column (2) adds age, female, and urban-area controls. Column (3) adds the rich baseline control set in a standard pooled linear regression. Column (4) is a fixed-effects model with rich controls and individual intercepts. Column (5) adds wave fixed effects to the fixed-effects model. Column (6) is a random-effects model with rich controls. Blank cells indicate covariates that are absorbed by individual fixed effects or otherwise not identified in a given specification. The rich-control models also include broad ethnicity and region dummies, which are omitted from display for compactness. The row labeled Close friends refers to the logged control \\texttt{log(1 + close\\_friends)} used in the regressions. Significance: *** \\(p < 0.01\\), ** \\(p < 0.05\\), * \\(p < 0.10\\).",
    "\\end{minipage}",
    "\\end{table}"
  )

  writeLines(lines, file.path(tables_dir, file_name))
  write.csv(
    bind_rows(numeric_rows),
    file.path(tables_dir, sub("\\.tex$", "_numeric.csv", file_name)),
    row.names = FALSE
  )
}

fit_panel_models <- function(outcome_name) {
  panel_data <- prepare_common_sample(dt, outcome_name, include_wave = TRUE)
  pdata <- pdata.frame(panel_data, index = c("pidp", "wave"), drop.index = TRUE)

  rhs_nowave <- c(treatment_var, many_controls)
  rhs_wave <- c(treatment_var, many_controls, "wave_f")
  formula_nowave <- as.formula(paste(outcome_name, "~", paste(rhs_nowave, collapse = " + ")))
  formula_wave <- as.formula(paste(outcome_name, "~", paste(rhs_wave, collapse = " + ")))

  pooled_display <- lm(formula_nowave, data = panel_data)
  fixed_display <- plm(formula_wave, data = pdata, model = "within", effect = "individual")
  random_display <- plm(
    formula_nowave,
    data = pdata,
    model = "random",
    effect = "individual",
    random.method = "swar"
  )

  pooled_nowave <- plm(formula_nowave, data = pdata, model = "pooling")
  fixed_nowave <- plm(formula_nowave, data = pdata, model = "within", effect = "individual")
  pooled_wave <- plm(formula_wave, data = pdata, model = "pooling")
  random_nowave <- plm(
    formula_nowave,
    data = pdata,
    model = "random",
    effect = "individual",
    random.method = "swar"
  )

  list(
    data = panel_data,
    pooled_display = pooled_display,
    fixed_display = fixed_display,
    random_display = random_display,
    pooled_nowave = pooled_nowave,
    fixed_nowave = fixed_nowave,
    pooled_wave = pooled_wave,
    random_nowave = random_nowave
  )
}

extract_plm_term <- function(model, term) {
  vcov_mat <- vcovHC(model, method = "arellano", type = "HC1", cluster = "group")
  ct <- coeftest(model, vcov. = vcov_mat)
  if (!(term %in% rownames(ct))) {
    return(c(estimate = NA_real_, se = NA_real_, p_value = NA_real_))
  }
  p_col <- grep("Pr\\(", colnames(ct), value = TRUE)[1]
  c(
    estimate = unname(ct[term, "Estimate"]),
    se = unname(ct[term, "Std. Error"]),
    p_value = unname(ct[term, p_col])
  )
}

preferred_model_label <- function(f_p, lm_p, hausman_p) {
  if (!is.na(hausman_p) && hausman_p < 0.05) {
    return("Fixed effects")
  }
  if (!is.na(lm_p) && lm_p < 0.05) {
    return("Random effects")
  }
  if (!is.na(f_p) && f_p < 0.05) {
    return("Fixed effects")
  }
  "Pooled OLS"
}

write_panel_tests_table <- function() {
  ordered_outcomes <- c(
    life_dissatisfaction = "Life dissatisfaction",
    loneliness = "Loneliness",
    schoolwork_dissatisfaction = "School-work dissatisfaction",
    school_dissatisfaction = "School dissatisfaction"
  )

  row_labels <- c(
    "Pooled classical linear regression model",
    "Fixed effects model",
    "Random effects model",
    "F test p-value",
    "LM test p-value",
    "Hausman p-value",
    "N"
  )

  display <- data.frame(
    row_label = row_labels,
    stringsAsFactors = FALSE
  )
  numeric_rows <- list()

  for (outcome_name in names(ordered_outcomes)) {
    fitted <- fit_panel_models(outcome_name)

    pooled_term <- extract_lm_term(fitted$pooled_display, fitted$data, treatment_var)
    random_term <- extract_plm_term(fitted$random_display, treatment_var)
    fixed_term <- extract_plm_term(fitted$fixed_display, treatment_var)

    f_test <- tryCatch(pFtest(fitted$fixed_nowave, fitted$pooled_nowave), error = function(e) NULL)
    lm_test <- tryCatch(plmtest(fitted$pooled_nowave, type = "bp"), error = function(e) NULL)
    hausman_test <- tryCatch(phtest(fitted$fixed_nowave, fitted$random_nowave), error = function(e) NULL)

    f_p <- if (is.null(f_test)) NA_real_ else f_test$p.value
    lm_p <- if (is.null(lm_test)) NA_real_ else lm_test$p.value
    hausman_p <- if (is.null(hausman_test)) NA_real_ else hausman_test$p.value

    display[[ordered_outcomes[[outcome_name]]]] <- c(
      format_coef_cell(pooled_term["estimate"], pooled_term["se"], pooled_term["p_value"]),
      format_coef_cell(fixed_term["estimate"], fixed_term["se"], fixed_term["p_value"]),
      format_coef_cell(random_term["estimate"], random_term["se"], random_term["p_value"]),
      format_p_value(f_p),
      format_p_value(lm_p),
      format_p_value(hausman_p),
      format(nrow(fitted$data), big.mark = ",")
    )

    numeric_rows[[length(numeric_rows) + 1]] <- data.frame(
      outcome = outcome_name,
      outcome_label = outcomes[[outcome_name]],
      pooled_estimate = pooled_term["estimate"],
      pooled_se = pooled_term["se"],
      pooled_p_value = pooled_term["p_value"],
      fixed_estimate = fixed_term["estimate"],
      fixed_se = fixed_term["se"],
      fixed_p_value = fixed_term["p_value"],
      random_estimate = random_term["estimate"],
      random_se = random_term["se"],
      random_p_value = random_term["p_value"],
      f_test_p_value = f_p,
      lm_test_p_value = lm_p,
      hausman_p_value = hausman_p,
      preferred_model = preferred_model_label(f_p, lm_p, hausman_p),
      n = nrow(fitted$data),
      stringsAsFactors = FALSE
    )
  }

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Panel-model comparison and specification tests}",
    "\\label{tab:panel_model_tests}",
    "\\scriptsize",
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{lcccc}",
    "\\hline",
    paste(c("Model or test", escape_latex(unname(ordered_outcomes))), collapse = " & "),
    "\\\\",
    "\\hline"
  )

  for (i in seq_len(nrow(display))) {
    row <- as.character(display[i, ])
    lines <- c(lines, paste(escape_latex(row), collapse = " & "), "\\\\")
  }

  lines <- c(
    lines,
    "\\hline",
    "\\end{tabular}",
    "}",
    "\\begin{minipage}{0.95\\textwidth}",
    "\\footnotesize Notes: The first three rows reproduce the weekday social-media coefficient from the outcome-specific linear-results tables, specifically columns (3), (5), and (6). All rows use the same outcome-specific estimation samples as those tables. Robust standard errors are clustered at the respondent level. The F test evaluates the pooled-versus-fixed-effects comparison on the matched rich-control specifications without wave dummies. The Breusch-Pagan LM test evaluates the pooled-versus-random-effects comparison for the matched rich-control specifications, using the unbalanced-panel implementation in \\texttt{plm}. The Hausman test evaluates the fixed-versus-random-effects comparison on the matched rich-control specifications without wave dummies. Small p-values therefore reject the simpler or stronger-assumption model in favour of the more flexible alternative. Significance: *** \\(p < 0.01\\), ** \\(p < 0.05\\), * \\(p < 0.10\\).",
    "\\end{minipage}",
    "\\end{table}"
  )

  writeLines(lines, file.path(tables_dir, "table_panel_model_tests.tex"))
  write.csv(bind_rows(numeric_rows), file.path(tables_dir, "table_panel_model_tests_numeric.csv"), row.names = FALSE)
}

make_panel_schematic <- function() {
  set.seed(20260518)

  toy_data <- data.frame(
    adolescent = c(
      rep("Adolescent 1", 4),
      rep("Adolescent 2", 3),
      rep("Adolescent 3", 2),
      rep("Adolescent 4", 3)
    ),
    wave_index = c(
      0, 1, 2, 3,
      0, 1, 3,
      1, 2,
      0, 2, 3
    ),
    stringsAsFactors = FALSE
  ) %>%
    mutate(adolescent = factor(adolescent, levels = paste("Adolescent", 1:4))) %>%
    arrange(adolescent, wave_index)

  wave_labels <- c("l", "m", "n", "o")
  intercepts <- c(2.0, 3.7, 2.8, 4.3)
  toy_data$mental_health <- intercepts[as.integer(toy_data$adolescent)] +
    0.26 * toy_data$wave_index +
    rnorm(nrow(toy_data), sd = 0.22)
  toy_data$mental_health <- pmin(7, pmax(1, toy_data$mental_health))
  toy_data$wave_label <- factor(
    wave_labels[toy_data$wave_index + 1],
    levels = wave_labels
  )

  pooled_fit <- lm(mental_health ~ wave_index, data = toy_data)
  fixed_fit <- lm(mental_health ~ wave_index + adolescent, data = toy_data)
  random_fit <- lme(
    mental_health ~ wave_index,
    random = ~1 | adolescent,
    data = toy_data,
    method = "REML"
  )

  pooled_grid <- data.frame(
    wave_index = seq(0, 3, length.out = 100)
  )
  pooled_grid$fit <- predict(pooled_fit, newdata = pooled_grid)

  grid <- expand.grid(
    wave_index = seq(0, 3, length.out = 100),
    adolescent = levels(toy_data$adolescent),
    KEEP.OUT.ATTRS = FALSE
  )
  grid$adolescent <- factor(grid$adolescent, levels = levels(toy_data$adolescent))

  fixed_grid <- grid
  fixed_grid$fit <- predict(fixed_fit, newdata = fixed_grid)

  random_grid <- grid
  random_grid$fit_subject <- predict(random_fit, newdata = random_grid, level = 1)
  random_population <- data.frame(
    wave_index = seq(0, 3, length.out = 100)
  )
  random_population$fit <- predict(random_fit, newdata = random_population, level = 0)

  base_theme <- theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey88"),
      plot.title = element_text(face = "bold", size = 11),
      axis.title = element_text(size = 10)
    )

  palette_values <- c("#de8fb8", "#efc7e8", "#c9c3f4", "#8f82d8")
  common_line_color <- "#4b5563"
  x_breaks <- 0:3
  x_labels <- c("l", "m", "n", "o")

  p_a <- ggplot(toy_data, aes(x = wave_index, y = mental_health, color = adolescent)) +
    geom_point(size = 2.1) +
    geom_line(aes(group = adolescent), alpha = 0.35) +
    geom_line(
      data = pooled_grid,
      aes(x = wave_index, y = fit),
      inherit.aes = FALSE,
      color = common_line_color,
      linewidth = 0.95
    ) +
    scale_color_manual(values = palette_values) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = c(1, 7), breaks = 1:7) +
    labs(title = "A. Pooled OLS", x = "Survey wave", y = "Life dissatisfaction") +
    base_theme

  p_b <- ggplot(toy_data, aes(x = wave_index, y = mental_health, color = adolescent)) +
    geom_point(size = 2.1) +
    geom_line(
      data = fixed_grid,
      aes(x = wave_index, y = fit, group = adolescent, color = adolescent),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    scale_color_manual(values = palette_values) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = c(1, 7), breaks = 1:7) +
    labs(title = "B. Fixed effects", x = "Survey wave", y = NULL) +
    base_theme

  p_c <- ggplot(toy_data, aes(x = wave_index, y = mental_health, color = adolescent)) +
    geom_point(size = 2.1) +
    geom_line(
      data = random_grid,
      aes(x = wave_index, y = fit_subject, group = adolescent, color = adolescent),
      inherit.aes = FALSE,
      linewidth = 0.8,
      alpha = 0.7
    ) +
    geom_line(
      data = random_population,
      aes(x = wave_index, y = fit),
      inherit.aes = FALSE,
      color = common_line_color,
      linewidth = 0.95
    ) +
    scale_color_manual(values = palette_values) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = c(1, 7), breaks = 1:7) +
    labs(title = "C. Random effects", x = "Survey wave", y = NULL) +
    base_theme

  png(
    filename = file.path(figures_dir, "panel_model_schematic.png"),
    width = 3300,
    height = 1100,
    res = 300
  )
  grid.arrange(p_a, p_b, p_c, ncol = 3)
  dev.off()
}

write_linear_model_table(
  outcome_name = main_outcome,
  file_name = "table_linear_benchmark_ols_life.tex",
  caption = "Linear model estimates for life dissatisfaction",
  label = "tab:pooled_ols_life"
)
write_linear_model_table(
  outcome_name = "loneliness",
  file_name = "table_linear_benchmark_ols_loneliness.tex",
  caption = "Linear model estimates for loneliness",
  label = "tab:linear_models_loneliness"
)
write_linear_model_table(
  outcome_name = "schoolwork_dissatisfaction",
  file_name = "table_linear_benchmark_ols_schoolwork.tex",
  caption = "Linear model estimates for school-work dissatisfaction",
  label = "tab:linear_models_schoolwork"
)
write_linear_model_table(
  outcome_name = "school_dissatisfaction",
  file_name = "table_linear_benchmark_ols_school.tex",
  caption = "Linear model estimates for school dissatisfaction",
  label = "tab:linear_models_school"
)
write_panel_tests_table()
make_panel_schematic()
