# Chapter 5

# Extended summary for explaining Chapter 5:
# 1. This script rebuilds only the Chapter 5 result tables that are
#    actually discussed in thesis.tex, and prints them directly in the console.
# 2. The cleaned input data are the pooled UKHLS youth panel for Waves l--o.
#    The panel identifier is pidp, so the same adolescent can appear in
#    multiple waves. That is why we cluster standard errors at respondent
#    level whenever repeated observations per pidp matter.
# 3. Table 5.1 is the linear benchmark table. For each outcome we first build
#    one common complete-case sample, and then estimate six specifications on
#    exactly that same sample:
#    (1) bivariate pooled OLS
#    (2) pooled OLS + age, female, urban
#    (3) pooled OLS + rich controls
#    (4) individual fixed effects + rich controls
#    (5) individual fixed effects + rich controls + wave fixed effects
#    (6) random effects + rich controls
#    Using one common sample per outcome means that differences across columns
#    come from the specification, not from a changing estimation sample.
# 4. Tables 5.2 and A5--A7 use partially linear double machine learning (PLR)
#    with the ordered weekday social-media treatment. In the PLR theory there
#    are two nuisance parts: the treatment model m(X) and the outcome-related
#    nuisance. In DoubleMLPLR this is implemented through two learners, usually
#    interpreted as l(X) = E[Y | X] and m(X) = E[D | X]. In this script we use
#    the same learner class twice in each comparison column, so for example the
#    elastic-net column uses elastic net both for the outcome nuisance and for
#    the treatment nuisance.
# 5. Tables 5.3 and A8--A10 use interactive regression model DML (IRM) with
#    the binary high-use treatment. Here we again need two nuisance functions,
#    but now they are conceptually different:
#    ml_g = outcome regression nuisance
#    ml_m = treatment propensity nuisance
#    Because the treatment is binary in IRM, ml_m must be a classifier rather
#    than a regressor. That is why the code uses pairs like regr.cv_glmnet for
#    ml_g and classif.cv_glmnet for ml_m.
# 6. IRM is numerically more fragile than PLR because it relies on a smaller
#    treated group and a binary treatment split. For that reason the thesis
#    comparison keeps only the learners that were stable in this setting:
#    Elastic net, Lasso, and Random forest. The neural net is used in the PLR
#    comparison, but not in the IRM comparison.
# 7. All DML fits use 5-fold cross-fitting. This means the sample is split into
#    five parts. In each round, the nuisance models are trained on four folds
#    and evaluated on the held-out fold. This is rotated over all five folds.
#    The goal is to avoid overfitting the nuisance functions on the same data
#    used to construct the treatment effect estimate, which is one of the key
#    ideas behind DML.
# 8. The learner comparison itself is not meant to prove that one learner is
#    globally optimal. The point is to see whether the estimated treatment
#    effect is reasonably robust across different plausible nuisance-model
#    classes: penalized linear models, tree-based models, and in PLR also a
#    small neural net benchmark.

library(dplyr)
library(data.table)
library(lmtest)
library(sandwich)
library(DoubleML)
library(mlr3)
library(mlr3learners)
library(plm)

lgr::get_logger("mlr3")$set_threshold("warn")

# 0. Setup -----------------------------------------------------------------
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

show_table <- function(title, x) {
  cat("\n", title, "\n", sep = "")
  if (inherits(x, "data.frame")) {
    print(as.data.frame(x), row.names = FALSE)
  } else {
    print(x)
  }
}

load_results_core_data <- function() {
  readRDS(data_path) %>%
    mutate(
      female_ind = if_else(as.character(sex_f) == "Female", 1, 0, missing = NA_real_),
      urban_ind = if_else(as.character(urban_f) == "Urban", 1, 0, missing = NA_real_)
    )
}

dt <- load_results_core_data()

# 1. Variables used in the results chapter ---------------------------------
outcomes <- c(
  loneliness = "Loneliness",
  life_dissatisfaction = "Life dissatisfaction",
  schoolwork_dissatisfaction = "School-work dissatisfaction",
  school_dissatisfaction = "School dissatisfaction"
)

continuous_treatment <- "social_media_weekday"
binary_treatment <- "high_social_media_weekday"

# These are the rich baseline controls used in the pooled DML tables.
# The factor variables are expanded later with model.matrix().
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

# In Table 5.1, the six columns are:
# (1) bivariate pooled OLS
# (2) pooled OLS + age, female, urban area
# (3) pooled OLS + rich baseline controls
# (4) fixed effects + rich controls
# (5) fixed effects + rich controls + wave fixed effects
# (6) random effects + rich controls

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

# 2. Helpers for the tables -------------------------
format_coef_cell_stack <- function(estimate, std_error, p_value) {
  if (any(is.na(c(estimate, std_error, p_value)))) {
    return("")
  }
  sprintf("%.3f%s (%.3f)", estimate, star_string(p_value), std_error)
}

# Make the core helpers safe to call interactively.
# If you run one function by itself in the console, this block recreates the
# small derived variables that the results code expects.
prepare_results_core_data <- function(data) {
  data <- as.data.frame(data)

  if (!("female_ind" %in% names(data)) && ("sex_f" %in% names(data))) {
    data$female_ind <- dplyr::if_else(
      as.character(data$sex_f) == "Female",
      1,
      0,
      missing = NA_real_
    )
  }

  if (!("urban_ind" %in% names(data)) && ("urban_f" %in% names(data))) {
    data$urban_ind <- dplyr::if_else(
      as.character(data$urban_f) == "Urban",
      1,
      0,
      missing = NA_real_
    )
  }

  if (!("wave_f" %in% names(data)) && ("wave" %in% names(data))) {
    wave_levels <- sort(unique(as.character(stats::na.omit(data$wave))))
    data$wave_f <- factor(
      data$wave,
      levels = wave_levels,
      labels = paste("Wave", toupper(wave_levels))
    )
  }

  data
}

# LM model Explained NL: vcovCL(model, cluster = data$pidp, type = "HC1")
# maakt een cluster-robuste variantie-covariantie matrix.
# Hier cluster op pidp, dus op respondentniveau.
# Dat betekent: standaardfouten worden gecorrigeerd voor het feit 
# dat dezelfde jongere meerdere observaties kan hebben over waves heen.
# coeftest(model, vcov. = vcov_mat) maakt een nette regressietabel met:
# coefficient estimates, standard errors, t-values, p-values
# maar dan met die cluster-robuste standaardfouten.
# if (!(term %in% rownames(ct)))
  #checkt of die variabele echt in het model zit.
# Soms zit een variabele er niet in, bijvoorbeeld door collineariteit of omdat 
# die in een bepaalde specificatie wegvalt.
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

extract_plm_term <- function(model, term) {
  # For panel models, use HC1/Arellano robust SE clustered at the group level.
  vcov_mat <- vcovHC(model, method = "arellano", type = "HC1", cluster = "group")
  ct <- coeftest(model, vcov. = vcov_mat)

  # Some rows disappear because of collinearity or model structure.
  if (!(term %in% rownames(ct))) {
    return(c(estimate = NA_real_, se = NA_real_, p_value = NA_real_))
  }

  # The exact p-value column name can differ slightly across plm outputs,
  # so search for the column that starts with "Pr(" rather than hard-coding it.
  p_col <- grep("Pr\\(", colnames(ct), value = TRUE)[1]
  c(
    estimate = unname(ct[term, "Estimate"]),
    se = unname(ct[term, "Std. Error"]),
    p_value = unname(ct[term, p_col])
  )
}

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

# This rebuilds Table 5.1 for life dissatisfaction and Appendix A2--A4 for the
# other outcomes.
# Explained:
# - first choose one common complete-case sample for this outcome
# - then fit the six specifications on exactly that same sample
# - then pull out the rows we want to show in the table
# - finally stack those numbers into a display data frame for printing
build_linear_model_display <- function(outcome_name) {
  # Rebuild the small helper columns if this function is called on its own.
  data_use <- prepare_results_core_data(dt)

  # Keep one common sample for this outcome so all six columns are directly
  # comparable and differences come from the model, not from changing N!
  vars_needed <- unique(c(
    outcome_name,
    continuous_treatment,
    life_rich_controls,
    "wave_f",
    "pidp",
    "wave"
  ))

  # Give a clean error message if required columns are missing,
  # instead of the vague "undefined columns selected".
  missing_vars <- setdiff(vars_needed, names(data_use))

  # If dt in the current R session is somehow stale or incomplete,
  # reload the cleaned analysis file directly and try once more.
  if (length(missing_vars) > 0) {
    data_use <- prepare_results_core_data(load_results_core_data())
    missing_vars <- setdiff(vars_needed, names(data_use))
  }

  if (length(missing_vars) > 0) {
    stop(
      paste(
        "build_linear_model_display() mist deze kolommen in dt:",
        paste(missing_vars, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  panel_data <- droplevels(
    data_use[
      complete.cases(data_use[, vars_needed, drop = FALSE]),
      vars_needed,
      drop = FALSE
    ]
  )

  # plm needs an explicit panel structure: person id + wave.
  pdata <- pdata.frame(panel_data, index = c("pidp", "wave"), drop.index = TRUE)

  # These six specs correspond to Table 5.1:
  # (1) bivariate pooled OLS
  # (2) pooled OLS + age/female/urban
  # (3) pooled OLS + rich controls
  # (4) individual fixed effects + rich controls
  # (5) individual fixed effects + rich controls + wave FE
  # (6) random effects + rich controls
  specs <- list(
    list(label = "(1)", kind = "lm", controls = character(0), wave_fe = FALSE),
    list(label = "(2)", kind = "lm", controls = life_basic_controls, wave_fe = FALSE),
    list(label = "(3)", kind = "lm", controls = life_rich_controls, wave_fe = FALSE),
    list(label = "(4)", kind = "fe", controls = life_rich_controls, wave_fe = FALSE),
    list(label = "(5)", kind = "fe", controls = life_rich_controls, wave_fe = TRUE),
    list(label = "(6)", kind = "re", controls = life_rich_controls, wave_fe = FALSE)
  )

  # These are the rows that are shown in the printed table.
  # Ethnicity and region are still in the regressions when relevant,
  # but are omitted from display to keep the table compact.
  term_rows <- list(
    list(label = "Weekday social-media use", term = continuous_treatment),
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
    # Build the right-hand side of the regression:
    # treatment + chosen controls + optional wave fixed effects.
    rhs <- c(continuous_treatment, spec$controls, if (spec$wave_fe) "wave_f")
    formula <- as.formula(paste(outcome_name, "~", paste(rhs, collapse = " + ")))

    # Pooled columns use lm().
    if (spec$kind == "lm") {
      return(lm(formula, data = panel_data))
    }

    # FE columns use the within estimator in plm().
    if (spec$kind == "fe") {
      return(plm(formula, data = pdata, model = "within", effect = "individual"))
    }

    # The last column is a random-effects model.
    plm(
      formula,
      data = pdata,
      model = "random",
      effect = "individual",
      # Swamy-Arora methode om de variantiecomponenten (indiv effect+idiosyncratisch term)
      # in het random-effects model te schatten
      random.method = "swar"
    )
  }

  extract_term <- function(model, spec, term) {
    # Use the matching extractor for the model class:
    # lm columns use extract_lm_term(), panel columns use extract_plm_term().
    if (spec$kind == "lm") {
      return(extract_lm_term(model, panel_data, term))
    }
    extract_plm_term(model, term)
  }

  model_fit_value <- function(model, spec) {
    # For pooled OLS, use the regular R-squared from summary(lm).
    if (spec$kind == "lm") {
      # Mooi afronden
      return(sprintf("%.3f", summary(model)$r.squared))
    }

    # For FE/RE panel models, use the reported panel R-squared when available.
    rsq <- summary(model)$r.squared["rsq"]
    if (is.na(rsq)) {
      return("")
    }
    sprintf("%.3f", rsq)
  }

  # Fit all six models first.
  models <- lapply(specs, fit_one)

  # Start the printed table with one first column listing the row labels.
  display <- data.frame(
    Variable = c(vapply(term_rows, `[[`, character(1), "label"), "$R^2$", "N"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  for (j in seq_along(models)) {
    # For each model/column:
    # 1. extract the coefficient row for every displayed variable
    # 2. format it as "estimate (se)" with stars
    # 3. add R-squared and N at the bottom
    column_cells <- c(
      vapply(term_rows, function(term_row) {
        extracted <- extract_term(models[[j]], specs[[j]], term_row$term)
        format_coef_cell_stack(extracted["estimate"], extracted["se"], extracted["p_value"])
      }, character(1)),
      model_fit_value(models[[j]], specs[[j]]),
      format(nobs(models[[j]]), big.mark = ",")
    )
    display[[specs[[j]]$label]] <- column_cells
  }

  display
}

# Print main linear table and appendix linear tables
table_5_1_linear <- build_linear_model_display("life_dissatisfaction")
show_table(
  "Table 5.1: Linear model estimates for life dissatisfaction",
  table_5_1_linear
)

appendix_a2_linear <- build_linear_model_display("loneliness")
show_table(
  "Appendix A2: Linear model estimates for loneliness",
  appendix_a2_linear
)

appendix_a3_linear <- build_linear_model_display("schoolwork_dissatisfaction")
show_table(
  "Appendix A3: Linear model estimates for school-work dissatisfaction",
  appendix_a3_linear
)

appendix_a4_linear <- build_linear_model_display("school_dissatisfaction")
show_table(
  "Appendix A4: Linear model estimates for school dissatisfaction",
  appendix_a4_linear
)


####################################### 

# This helper turns the numeric learner results into the simple three-row layout
# used in Tables 5.2, 5.3, and Appendix A5--A10.
build_dml_display_table <- function(results, variable_label) {
  # Each learner becomes one column; rows are estimate, SE, and N.
  learners <- results$learner
  display <- data.frame(
    Row = c(variable_label, "Standard error", "N"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(results))) {
    # Print the coefficient with stars on the first row,
    # then the standard error in parentheses, then the sample size.
    display[[learners[i]]] <- c(
      if (is.na(results$estimate[i])) "" else sprintf("%.3f%s", results$estimate[i], star_string(results$p_value[i])),
      if (is.na(results$std_error[i])) "" else sprintf("(%.3f)", results$std_error[i]),
      format(results$n[i], big.mark = ",")
    )
  }

  display
}

# Suppress the verbose DoubleML summary printout and keep only the numbers.
extract_dml_summary <- function(fit) {
  # fit$summary() prints to the console by default.
  # capture.output() hides that text, while summary_df keeps the actual values.
  invisible(capture.output(summary_df <- as.data.frame(fit$summary())))
  summary_df
}

# 3. Learner comparisons used in Tables 5.2, 5.3, and Appendix A5--A10 ----
prepare_dml_numeric_sample <- function(data, outcome, treatment_var, controls, include_wave = TRUE) {
  data <- prepare_results_core_data(data)

  # DML needs a clean outcome, treatment, and X-matrix.
  # Optionally add wave fixed effects to the X-set.
  x_cols <- unique(c(controls, if (include_wave) "wave_f"))
  vars_needed <- unique(c(outcome, treatment_var, x_cols, "pidp"))

  missing_vars <- setdiff(vars_needed, names(data))
  if (length(missing_vars) > 0) {
    stop(
      paste(
        "prepare_dml_numeric_sample() mist deze kolommen in dt:",
        paste(missing_vars, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Use one complete-case sample for the chosen outcome/treatment setup.
  data_model <- droplevels(
    data[
      complete.cases(data[, vars_needed, drop = FALSE]),
      vars_needed,
      drop = FALSE
    ]
  )

  # model.matrix() turns factor controls like region or wave into dummy columns,
  # which is what the learners need.
  x_matrix <- model.matrix(~ . - 1, data = data_model[, x_cols, drop = FALSE])

  # Build a numeric data.table with pid, outcome, treatment, and all X columns.
  dml_df <- data.table(
    pidp = data_model$pidp,
    outcome_value = data_model[[outcome]],
    treatment_value = data_model[[treatment_var]]
  )
  names(dml_df)[2:3] <- c(outcome, treatment_var)
  dml_df <- cbind(dml_df, as.data.table(x_matrix))

  list(
    data = dml_df,
    x_cols = setdiff(names(dml_df), c(outcome, treatment_var, "pidp"))
  )
}

make_plr_comparison_learners <- function(n_features) {
  # These are exactly the nuisance learners compared in Table 5.2:
  # Elastic net, Lasso, Random forest, and Neural net.
  list(
    # Explained NL: alpha = 0.5 betekent: half tussen ridge (alpha = 0) en lasso (alpha = 1)
    # dus: zowel shrinkage als wat variable selection
    # handig als covariaten onderling gecorreleerd zijn
    # s = "lambda.min" betekent: neem de lambda die via cross-validation de laagste predictiefout gaf
    "Elastic net" = mlr3::lrn("regr.cv_glmnet", s = "lambda.min", alpha = 0.5),
    # Explained NL:alpha = 1 is pure lasso (strengere benchmark)
    # die doet sterkere variable selection dan elastic net
    # handig als je wilt zien of de resultaten overeind blijven onder een “sparser” model
    # weer lambda.min omdat de strafparameter dan door CV gekozen wordt
    "Lasso" = mlr3::lrn("regr.cv_glmnet", s = "lambda.min", alpha = 1),
    "Random forest" = mlr3::lrn(
      "regr.ranger",
      # genoeg bomen voor redelijke stabiliteit
      # maar niet onnodig zwaar qua rekentijd
      num.trees = 300L,
      # heel gebruikelijke standaardheuristiek in random forests
      # per split kijkt het model maar naar een subset van features
      # dat helpt tegen overfit en maakt bomen diverser
      mtry = max(1L, floor(sqrt(n_features))),
      # laat bladeren niet té klein worden
      # dus iets meer regularisatie, minder ruis-fit
      min.node.size = 5L,
      # puur praktisch: stabieler/reproduceerbaarder
      num.threads = 1L
    ),
    "Neural net" = mlr3::lrn(
      "regr.nnet",
      # klein hidden layer, expres niet te complex, 
      # want sample is niet enorm en nnet kan snel instabiel worden
      size = 5L,
      # weight decay = regularisatie, voorkomt dat het netwerk te wild gaat fitten
      decay = 0.01,
      # geef het model genoeg iteraties om te convergeren
      maxit = 300L,
      trace = FALSE
    )
  )
}

make_irm_comparison_learners <- function(n_features) {
  # IRM needs two learners:
  # ml_g for the outcome regression and ml_m for the treatment propensity.
  # ml_g is een regressiemodel, want de outcome is numeriek
  # ml_m is een classificatiemodel, want de treatment in IRM is binair: 0/1
  list(
    "Elastic net" = list(
      ml_g = mlr3::lrn("regr.cv_glmnet", s = "lambda.min", alpha = 0.5),
      ml_m = mlr3::lrn(
        "classif.cv_glmnet",
        s = "lambda.min",
        alpha = 0.5,
        predict_type = "prob"
      )
    ),
    "Lasso" = list(
      ml_g = mlr3::lrn("regr.cv_glmnet", s = "lambda.min", alpha = 1),
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
        # geen Neural Net: dat binary-treatment probleem is vaak instabieler
        # zeker bij deze data is de high-use groep kleiner
        # neural nets zijn dan sneller gevoelig voor instabiliteit, c
        # onvergence issues, of rare propensity predictions
      )
    )
  )
}

build_plr_method_comparison <- function(outcome, controls = rich_controls, include_wave = TRUE) {
  # Step 1: prepare the pooled numeric DML design for this outcome.
  prepared <- prepare_dml_numeric_sample(
    data = dt,
    outcome = outcome,
    treatment_var = continuous_treatment,
    controls = controls,
    include_wave = include_wave
  )

  # Step 2: create the clustered DoubleML data object.
  # Clustering is at respondent level because pidp appears in multiple waves.
  # pidp is de unieke respondent-ID
  # dezelfde persoon komt meerdere keren terug in de paneldata
  # dus observaties van dezelfde persoon over verschillende waves zijn 
  # niet onafhankelijk van elkaar
  # zonder clustering behandel je die herhaalde observaties te veel als volledig los van elkaar
  # dan worden standaardfouten vaak te klein
  # en lijken resultaten sneller significant dan terecht is
  dml_data <- DoubleMLClusterData$new(
    prepared$data,
    y_col = outcome,
    d_cols = continuous_treatment,
    x_cols = prepared$x_cols,
    cluster_cols = "pidp"
  )

  # Step 3: define the learner menu and an empty results collector.
  learners <- make_plr_comparison_learners(length(prepared$x_cols))
  numeric_rows <- list()

  for (i in seq_along(learners)) {
    # Run the same PLR setup with a different nuisance learner each time.
    # zet een random seed zodat de schatting reproduceerbaar is
    # sommige learners en cross-fitting gebruiken random splits / randomisatie
    # 100 + i zorgt dat elke learner een andere, maar vaste seed krijgt
    # Dus: Elastic net krijgt bijvoorbeeld seed 101, Lasso 102
    learner_name <- names(learners)[i]
    set.seed(100 + i)

    # dml_data: data in DoubleML-formaat
    # learners[[i]]$clone(): learner voor de outcome nuisance
    # learners[[i]]$clone(): learner voor de treatment nuisance
    # n_folds = 5: gebruik 5-fold cross-fitting
    #  Waarom twee keer dezelfde learner? Omdat in deze learner-comparison het idee is:
    # gebruik één bepaald type learner en laat datzelfde type zowel de outcome-functie 
    # als de treatment-functie schatten. Dus bijvoorbeeld:
    # Elastic net voor beide nuisance functies, Random forest voor beide nuisance functies
    # Waarom $clone()? Omdat je twee aparte learner-objecten nodig hebt:
    # één voor de ene nuisance stap, één voor de andere nuisance stap
    # Wil niet exact hetzelfde object hergebruiken in memory.
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

    # If a learner fails, keep the row but mark the estimate as missing.
    if (is.null(fit)) {
      numeric_rows[[length(numeric_rows) + 1]] <- data.frame(
        framework = "PLR",
        outcome = outcome,
        outcome_label = outcomes[[outcome]],
        learner = learner_name,
        estimate = NA_real_,
        std_error = NA_real_,
        p_value = NA_real_,
        n = nrow(prepared$data),
        stringsAsFactors = FALSE
      )
      next
    }

    # Keep only the estimate, SE, p-value, and N that the printed table needs.
    summary_df <- extract_dml_summary(fit)
    numeric_rows[[length(numeric_rows) + 1]] <- data.frame(
      framework = "PLR",
      outcome = outcome,
      outcome_label = outcomes[[outcome]],
      learner = learner_name,
      estimate = unname(summary_df[1, 1]),
      std_error = unname(summary_df[1, 2]),
      p_value = unname(summary_df[1, 4]),
      n = nrow(prepared$data),
      stringsAsFactors = FALSE
    )
  }

  # Bind all learner rows into one table for this outcome.
  bind_rows(numeric_rows)
}

build_irm_method_comparison <- function(outcome, controls = rich_controls, include_wave = TRUE) {
  # Step 1: prepare the pooled numeric DML design, but now for the binary
  # high-use treatment instead of the ordered weekday-use treatment.
  prepared <- prepare_dml_numeric_sample(
    data = dt,
    outcome = outcome,
    treatment_var = binary_treatment,
    controls = controls,
    include_wave = include_wave
  )

  # Step 2: define the learner menu and check that treatment is truly binary.
  learners <- make_irm_comparison_learners(length(prepared$x_cols))
  treatment_values <- sort(unique(prepared$data[[binary_treatment]]))

  if (
    length(prepared$x_cols) == 0 ||
    length(treatment_values) != 2 ||
    !all(treatment_values %in% c(0, 1))
  ) {
    # If treatment is not 0/1 in the prepared sample, return missing rows
    # instead of crashing.
    return(bind_rows(lapply(names(learners), function(learner_name) {
      data.frame(
        framework = "IRM",
        outcome = outcome,
        outcome_label = outcomes[[outcome]],
        learner = learner_name,
        estimate = NA_real_,
        std_error = NA_real_,
        p_value = NA_real_,
        n = nrow(prepared$data),
        stringsAsFactors = FALSE
      )
    })))
  }

  treatment_counts <- table(prepared$data[[binary_treatment]])
  if (any(treatment_counts < 30)) {
    # IRM becomes unstable if one treatment arm is too small,
    # so again return missing rows rather than forcing estimation.
    return(bind_rows(lapply(names(learners), function(learner_name) {
      data.frame(
        framework = "IRM",
        outcome = outcome,
        outcome_label = outcomes[[outcome]],
        learner = learner_name,
        estimate = NA_real_,
        std_error = NA_real_,
        p_value = NA_real_,
        n = nrow(prepared$data),
        stringsAsFactors = FALSE
      )
    })))
  }

  # Step 3: create the clustered DoubleML object for the binary-treatment setup.
  dml_data <- DoubleMLClusterData$new(
    prepared$data,
    y_col = outcome,
    d_cols = binary_treatment,
    x_cols = prepared$x_cols,
    cluster_cols = "pidp"
  )

  numeric_rows <- list()

  for (i in seq_along(learners)) {
    # Run the IRM ATE estimator for each learner pair.
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

    # Keep a missing row if that learner pair fails.
    if (is.null(fit)) {
      numeric_rows[[length(numeric_rows) + 1]] <- data.frame(
        framework = "IRM",
        outcome = outcome,
        outcome_label = outcomes[[outcome]],
        learner = learner_name,
        estimate = NA_real_,
        std_error = NA_real_,
        p_value = NA_real_,
        n = nrow(prepared$data),
        stringsAsFactors = FALSE
      )
      next
    }

    # Keep only the summary numbers that feed the printed appendix/main table.
    summary_df <- extract_dml_summary(fit)
    numeric_rows[[length(numeric_rows) + 1]] <- data.frame(
      framework = "IRM",
      outcome = outcome,
      outcome_label = outcomes[[outcome]],
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

# Precompute all learner-comparison results once,
# then subset them below when printing each main-table or appendix-table output.
plr_method_results <- bind_rows(lapply(names(outcomes), build_plr_method_comparison))
irm_method_results <- bind_rows(lapply(names(outcomes), build_irm_method_comparison))

# Main text: Table 5.2 in thesis.tex.
show_table(
  "Table 5.2: PLR-DML learner comparison for life dissatisfaction",
  build_dml_display_table(
    plr_method_results[plr_method_results$outcome == "life_dissatisfaction", ],
    "Weekday social-media use"
  )
)

# Main text: Table 5.3 in thesis.tex.
show_table(
  "Table 5.3: IRM-DML learner comparison for life dissatisfaction",
  build_dml_display_table(
    irm_method_results[irm_method_results$outcome == "life_dissatisfaction", ],
    "High social-media use"
  )
)

# Appendix A5--A7: PLR-DML learner tables for the other outcomes.
show_table(
  "Appendix A5: PLR-DML learner comparison for loneliness",
  build_dml_display_table(
    plr_method_results[plr_method_results$outcome == "loneliness", ],
    "Weekday social-media use"
  )
)
show_table(
  "Appendix A6: PLR-DML learner comparison for school-work dissatisfaction",
  build_dml_display_table(
    plr_method_results[plr_method_results$outcome == "schoolwork_dissatisfaction", ],
    "Weekday social-media use"
  )
)
show_table(
  "Appendix A7: PLR-DML learner comparison for school dissatisfaction",
  build_dml_display_table(
    plr_method_results[plr_method_results$outcome == "school_dissatisfaction", ],
    "Weekday social-media use"
  )
)

# Appendix A8--A10: IRM-DML learner tables for the other outcomes.
show_table(
  "Appendix A8: IRM-DML learner comparison for loneliness",
  build_dml_display_table(
    irm_method_results[irm_method_results$outcome == "loneliness", ],
    "High social-media use"
  )
)
show_table(
  "Appendix A9: IRM-DML learner comparison for school-work dissatisfaction",
  build_dml_display_table(
    irm_method_results[irm_method_results$outcome == "schoolwork_dissatisfaction", ],
    "High social-media use"
  )
)
show_table(
  "Appendix A10: IRM-DML learner comparison for school dissatisfaction",
  build_dml_display_table(
    irm_method_results[irm_method_results$outcome == "school_dissatisfaction", ],
    "High social-media use"
  )
)

cat("\nResults core tables were printed to the console.\n")
