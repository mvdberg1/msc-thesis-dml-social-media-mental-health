# Chapter 3 data construction and descriptive outputs.
#
# This script reads the licensed UKHLS youth files, builds the cleaned Waves
# l--o analysis panel, and regenerates the descriptive tables and figures used
# in the thesis.

library(haven)
library(dplyr)
library(ggplot2)
library(scales)

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
out_dir <- file.path(project_dir, "data", "analysis")
tables_dir <- file.path(project_dir, "tables")
figures_dir <- file.path(project_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

all_waves <- letters[1:15]
wave_numbers <- setNames(seq_along(all_waves), all_waves)

#  read one raw variable directly from a chosen wave file
get_raw_var <- function(wave, var) {
  raw_name <- paste0(wave, "_", var)
  read_dta(file.path(root, paste0(wave, "_youth.dta")), col_select = all_of(raw_name))[[1]]
}

#  pull the questionnaire label attached to a raw variable
get_var_label <- function(wave, var) {
  label <- attr(get_raw_var(wave, var), "label")
  if (is.null(label)) "" else label
}

#  escape special characters before writing latex tables
escape_latex <- function(x) {
  x <- gsub("#", "\\#", x, fixed = TRUE)
  x <- gsub("_", "\\_", x, fixed = TRUE)
  x <- gsub("%", "\\%", x, fixed = TRUE)
  x <- gsub("&", "\\&", x, fixed = TRUE)
  x
}

#  format numeric output for printed tables
format_number <- function(x, digits = 2) {
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits))
}

#  use wave l as the rich questionnaire file for candidate discovery
metadata_search_wave <- "l"

#  load the metadata header used to search for candidate variables
baseline_header <- read_dta(file.path(root, paste0(metadata_search_wave, "_youth.dta")), n_max = 0)
baseline_metadata <- data.frame(
  variable = sub(paste0("^", metadata_search_wave, "_"), "", names(baseline_header)),
  survey_label = vapply(baseline_header, function(x) {
    label <- attr(x, "label")
    if (is.null(label)) "" else label
  }, character(1)),
  stringsAsFactors = FALSE
)

#  define the keyword searches used to find treatments, outcomes, and controls
search_plan <- data.frame(
  search_group = c(
    "Treatment / exposure candidates",
    "Outcome candidates",
    "Baseline covariate candidates"
  ),
  keywords_used = c(
    "social media; profile",
    "lonely; happiness",
    "friend; family; tv; device; age; sex; ethnic; country; urban; region; parent; grandparent; sibling; household"
  ),
  regex = c(
    "social media|profile",
    "lonely|happiness",
    "friend|family|tv|device|age|sex|ethnic|country|urban|region|parent|grandparent|sibling|household|\\bhh\\b"
  ),
  stringsAsFactors = FALSE
)

#  search the baseline metadata with those keyword groups
search_results <- bind_rows(lapply(seq_len(nrow(search_plan)), function(i) {
  hits <- baseline_metadata %>%
    filter(
      grepl(search_plan$regex[i], variable, ignore.case = TRUE) |
        grepl(search_plan$regex[i], survey_label, ignore.case = TRUE)
    ) %>%
    mutate(
      search_group = search_plan$search_group[i],
      keywords_used = search_plan$keywords_used[i]
    ) %>%
    select(search_group, keywords_used, variable, survey_label)

  hits
}))

cat("\nStep 1/8: search the youth questionnaire metadata to find candidate variables\n")
cat("Search wave used for candidate discovery:", metadata_search_wave, "\n\n")
print(search_plan[, c("search_group", "keywords_used")], row.names = FALSE)
cat("\nSearch hits from the variable names and survey labels:\n")
print(search_results, row.names = FALSE)

#  record which variables we keep after inspecting the candidate hits
chosen_variables <- data.frame(
  analysis_role = c(
    rep("treatment / exposure", 3),
    rep("outcome", 4),
    rep("baseline covariate", 18),
    rep("additional social variable", 3)
  ),
  source_variable = c(
    "ypnetcht",
    "ypnetchtw",
    "ypsocweb",
    "yplonely",
    "yphlf",
    "yphsw",
    "yphsc",
    "age_dv",
    "sex_dv",
    "ethn_dv",
    "country",
    "urban_dv",
    "gor_dv",
    "npns_dv",
    "ngrp_dv",
    "nnssib_dv",
    "ypnpal",
    "ypeatlivu",
    "yptvvidhrs",
    paste0("ypdevice", 1:6),
    "ypvirfnd",
    "ypfndmeet",
    "ypfndonl"
  ),
  why_kept = c(
    "main weekday treatment",
    "weekend treatment for robustness checks",
    "social-media profile indicator",
    "baseline loneliness outcome",
    "life wellbeing outcome",
    "school-work wellbeing outcome",
    "school wellbeing outcome",
    "age control",
    "sex control",
    "ethnicity control",
    "country control",
    "urban-rural control",
    "region control",
    "all-parent-figures-in-household control",
    "grandparents-in-household control",
    "all-siblings-in-household control",
    "close-friends control",
    "family-meals control",
    "screen-time control",
    "smartphone control",
    "mobile-phone control",
    "tablet control",
    "television control",
    "gaming-console control",
    "laptop-desktop control",
    "online-only-friends measure",
    "in-person-friend-contact measure",
    "online-friend-contact measure"
  ),
  stringsAsFactors = FALSE
) %>%
  left_join(baseline_metadata, by = c("source_variable" = "variable")) %>%
  select(analysis_role, source_variable, survey_label, why_kept)

cat("\nVariables kept after inspecting the search hits:\n")
print(chosen_variables, row.names = FALSE)

#  derive the selected variable groups from the explicit post-search choice table
selected_treatment_vars <- chosen_variables$source_variable[chosen_variables$analysis_role == "treatment / exposure"]
selected_outcome_vars <- chosen_variables$source_variable[chosen_variables$analysis_role == "outcome"]
selected_baseline_covariate_vars <- chosen_variables$source_variable[chosen_variables$analysis_role == "baseline covariate"]
selected_additional_social_vars <- chosen_variables$source_variable[chosen_variables$analysis_role == "additional social variable"]
required_vars <- c(
  selected_treatment_vars,
  selected_outcome_vars,
  selected_baseline_covariate_vars,
  selected_additional_social_vars
)

#  keep the survey-design fields together with the chosen analysis variables
# psu = primary sampling unit
# strata = strata
# ythscui_xw = youth survey weight
stems <- c("psu", "strata", "ythscui_xw", required_vars)

#  check for each wave which chosen variables are present
wave_availability <- bind_rows(lapply(all_waves, function(wave) {
  wave_names <- names(read_dta(file.path(root, paste0(wave, "_youth.dta")), n_max = 0))
  wave_names <- sub(paste0("^", wave, "_"), "", wave_names)

  data.frame(
    wave = wave,
    variable = required_vars,
    available = required_vars %in% wave_names,
    stringsAsFactors = FALSE
  )
}))

#  record the first wave where each required variable appears
first_available_wave <- setNames(
  vapply(required_vars, function(var) {
    min(wave_availability$wave[wave_availability$variable == var & wave_availability$available])
  }, character(1)),
  required_vars
)

#  keep only chosen variables that have labelled raw response codes
coded_survey_vars <- required_vars[
  vapply(required_vars, function(var) {
    !is.null(attr(get_raw_var(first_available_wave[[var]], var), "labels"))
  }, logical(1))
]

#  summarise the first-appearance pattern used in the wave choice
wave_choice_summary <- data.frame(
  variable = required_vars,
  survey_label = vapply(required_vars, function(var) {
    get_var_label(first_available_wave[[var]], var)
  }, character(1)),
  first_available_wave = unname(first_available_wave),
  availability_group = ifelse(
    match(unname(first_available_wave), all_waves) == 1,
    "Available from wave a onward",
    paste0("First appears in wave ", unname(first_available_wave))
  ),
  stringsAsFactors = FALSE
)

#  compare all consecutive wave blocks that end in wave o
candidate_wave_blocks <- bind_rows(lapply(seq_along(all_waves), function(i) {
  candidate_waves <- all_waves[i:length(all_waves)]
  vars_available_in_block <- vapply(required_vars, function(var) {
    all(wave_availability$available[wave_availability$variable == var & wave_availability$wave %in% candidate_waves])
  }, logical(1))

  data.frame(
    start_wave = all_waves[i],
    end_wave = tail(candidate_waves, 1),
    n_waves = length(candidate_waves),
    n_required_vars_available = sum(vars_available_in_block),
    all_required_vars_available = all(vars_available_in_block),
    stringsAsFactors = FALSE
  )
}))

#  pick the first feasible block that contains every required variable
selected_start_wave <- candidate_wave_blocks$start_wave[candidate_wave_blocks$all_required_vars_available][1]
if (is.na(selected_start_wave)) {
  stop("No feasible consecutive wave block found for the required variables.", call. = FALSE)
}

waves <- all_waves[match(selected_start_wave, all_waves):length(all_waves)]

cat("\nStep 2/8: inspect all waves before selecting a baseline block\n")
cat("Chosen variables used for wave selection:", length(required_vars), "\n\n")
cat("First wave in which each required variable appears:\n")
print(wave_choice_summary, row.names = FALSE)
cat("\nCandidate consecutive blocks ending in wave o:\n")
print(candidate_wave_blocks, row.names = FALSE)
cat(
  "\nChosen baseline block:",
  paste(waves, collapse = ", "),
  "(Wave", wave_numbers[waves[1]], "to Wave", wave_numbers[tail(waves, 1)], ")\n"
)

#  build one codebook row from the first wave where a variable appears
get_codebook_entry <- function(var) {
  first_wave <- first_available_wave[[var]]
  x <- get_raw_var(first_wave, var)
  labels <- attr(x, "labels")
  labels <- labels[order(unname(labels))]

  special_label_names <- c("missing", "inapplicable", "refusal", "don't know", "inconsistent")
  special_codes <- unname(labels[tolower(names(labels)) %in% special_label_names])

  observed_codes <- sort(unique(as.numeric(zap_labels(x))))
  observed_codes <- observed_codes[!is.na(observed_codes)]
  valid_codes <- observed_codes[!(observed_codes %in% special_codes)]
  valid_labels <- vapply(valid_codes, function(code) {
    idx <- match(code, unname(labels))
    if (is.na(idx)) {
      as.character(code)
    } else {
      names(labels)[idx]
    }
  }, character(1))

  bounds_rule <- if (any(valid_codes == 0)) {
    "Keep all observed non-special response codes, including 0 when it is a substantive category."
  } else {
    "Keep all observed non-special response codes after dropping labelled missing/refusal/inapplicable categories."
  }

  data.frame(
    source_variable = var,
    survey_label = get_var_label(first_wave, var),
    first_wave_used = first_wave,
    all_labeled_codes = paste0(unname(labels), "=", names(labels), collapse = "; "),
    dropped_special_codes = paste0(special_codes, "=", names(labels)[match(special_codes, unname(labels))], collapse = "; "),
    valid_codes_used = paste(valid_codes, collapse = ", "),
    valid_value_labels = paste0(valid_codes, "=", valid_labels, collapse = "; "),
    bounds_rule = bounds_rule,
    stringsAsFactors = FALSE
  )
}

#  collect the raw labels and valid codes that survive data cleaning
raw_scale_reference <- bind_rows(lapply(coded_survey_vars, get_codebook_entry))
bounded_vars <- setNames(
  lapply(coded_survey_vars, function(var) {
    first_wave <- first_available_wave[[var]]
    x <- get_raw_var(first_wave, var)
    labels <- attr(x, "labels")
    labels <- labels[order(unname(labels))]
    special_label_names <- c("missing", "inapplicable", "refusal", "don't know", "inconsistent")
    special_codes <- unname(labels[tolower(names(labels)) %in% special_label_names])

    observed_codes <- sort(unique(as.numeric(zap_labels(x))))
    observed_codes <- observed_codes[!is.na(observed_codes)]
    as.integer(observed_codes[!(observed_codes %in% special_codes)])
  }),
  coded_survey_vars
)

cat("\nStep 3/8: inspect all labelled answer options and derive the valid response codes\n")
print(raw_scale_reference, row.names = FALSE)
cat(
  "\nFor yphlf, yphsw, and yphsc, the raw labels run from 1 = completely happy to 7 = not at all happy.\n",
  "So higher raw values indicate worse wellbeing, and the cleaned outcomes keep the raw scale",
  sep = ""
)

#  convert labelled survey fields to plain numeric values
clean_numeric <- function(x) {
  x <- as.numeric(zap_labels(x))
  x[x < 0] <- NA_real_
  x
}

#  map numeric survey codes into readable factor labels
factor_from_codes <- function(x, labels) {
  factor(labels[as.character(x)], levels = unname(labels))
}

#  read one wave, keep the needed columns, and strip the wave prefix
read_one_wave <- function(wave) {
  df <- as.data.frame(read_dta(file.path(root, paste0(wave, "_youth.dta"))))
  keep <- c("pidp", paste0(wave, "_", stems))
  df <- df[intersect(keep, names(df))]
  names(df) <- sub(paste0("^", wave, "_"), "", names(df))
  df$wave <- wave
  df$wave_number <- unname(wave_numbers[wave])
  df
}

#  stack the selected baseline waves into one person-wave dataset
youth <- bind_rows(lapply(waves, read_one_wave))

cat("\nStep 4/8: stack the selected baseline waves and clean the raw responses\n")
print(data.frame(wave = waves, wave_number = unname(wave_numbers[waves])), row.names = FALSE)
cat("Rows before cleaning:", nrow(youth), "\n")

#  turn all non-wave columns into numeric values before cleaning
numeric_cols <- setdiff(names(youth), "wave")
youth[numeric_cols] <- lapply(youth[numeric_cols], clean_numeric)

#  drop any remaining response codes outside the valid raw ranges
for (var in names(bounded_vars)) {
  allowed <- bounded_vars[[var]]
  bad <- !(youth[[var]] %in% allowed) & !is.na(youth[[var]])
  youth[[var]][bad] <- NA_real_
}

cat("Negative special codes are set to NA, and only the valid codes shown above are retained.\n")

#  construct the treatment, outcome, and close-friends variables
youth$social_media_weekday <- youth$ypnetcht
youth$social_media_weekend <- youth$ypnetchtw
youth$high_social_media_weekday <- ifelse(
  is.na(youth$social_media_weekday),
  NA_integer_,
  as.integer(youth$social_media_weekday >= 4)
)
youth$has_social_profile <- ifelse(
  is.na(youth$ypsocweb),
  NA_integer_,
  as.integer(youth$ypsocweb == 1)
)
youth$loneliness <- youth$yplonely
youth$life_dissatisfaction <- youth$yphlf
youth$schoolwork_dissatisfaction <- youth$yphsw
youth$school_dissatisfaction <- youth$yphsc
#  keep zero defined and reduce the leverage of the right tail in richer models
youth$close_friends <- ifelse(youth$ypnpal > 100, NA_real_, youth$ypnpal)
youth$close_friends_log <- log1p(youth$close_friends)

#  rename the survey weight so it is easier to reuse later
names(youth)[names(youth) == "ythscui_xw"] <- "youth_weight"

#  create readable factor versions of the key demographic variables
youth$wave_f <- factor(youth$wave, levels = waves, labels = paste("Wave", toupper(waves)))
youth$sex_f <- factor_from_codes(youth$sex_dv, c(`1` = "Male", `2` = "Female"))
#  keep country for documentation and possible robustness checks, even though the
#  richer main regression later uses region_f instead of including both jointly
youth$country_f <- factor_from_codes(
  youth$country,
  c(`1` = "England", `2` = "Wales", `3` = "Scotland", `4` = "Northern Ireland")
)
youth$urban_f <- factor_from_codes(youth$urban_dv, c(`1` = "Urban", `2` = "Rural"))
youth$region_f <- factor(youth$gor_dv)
youth$ethnicity_f <- factor(youth$ethn_dv)
youth$ethnicity_broad_f <- factor(
  case_when(
    youth$ethn_dv %in% 1:4 ~ "White",
    youth$ethn_dv %in% 5:8 ~ "Mixed",
    youth$ethn_dv %in% 9:13 ~ "Asian",
    youth$ethn_dv %in% 14:16 ~ "Black",
    youth$ethn_dv == 17 ~ "Arab",
    youth$ethn_dv == 97 ~ "Other",
    TRUE ~ NA_character_
  ),
  levels = c("White", "Mixed", "Asian", "Black", "Arab", "Other")
)

#  store the main variable-construction choices for terminal output
derived_variable_spec <- data.frame(
  source_variable = c(
    "ypnetcht",
    "ypnetchtw",
    "ypnetcht",
    "ypsocweb",
    "yplonely",
    "yphlf",
    "yphsw",
    "yphsc",
    "ypnpal",
    "close_friends"
  ),
  source_label = c(
    get_var_label(first_available_wave[["ypnetcht"]], "ypnetcht"),
    get_var_label(first_available_wave[["ypnetchtw"]], "ypnetchtw"),
    get_var_label(first_available_wave[["ypnetcht"]], "ypnetcht"),
    get_var_label(first_available_wave[["ypsocweb"]], "ypsocweb"),
    get_var_label(first_available_wave[["yplonely"]], "yplonely"),
    get_var_label(first_available_wave[["yphlf"]], "yphlf"),
    get_var_label(first_available_wave[["yphsw"]], "yphsw"),
    get_var_label(first_available_wave[["yphsc"]], "yphsc"),
    get_var_label(first_available_wave[["ypnpal"]], "ypnpal"),
    "Cleaned number of close friends"
  ),
  analysis_variable = c(
    "social_media_weekday",
    "social_media_weekend",
    "high_social_media_weekday",
    "has_social_profile",
    "loneliness",
    "life_dissatisfaction",
    "schoolwork_dissatisfaction",
    "school_dissatisfaction",
    "close_friends",
    "close_friends_log"
  ),
  transformation = c(
    "copy ypnetcht",
    "copy ypnetchtw",
    "1(ypnetcht >= 4)",
    "1(ypsocweb == 1)",
    "copy yplonely",
    "copy yphlf",
    "copy yphsw",
    "copy yphsc",
    "copy ypnpal, then set values > 100 to NA",
    "log1p(close_friends)"
  ),
  why_used = c(
    "Main weekday treatment",
    "Weekend treatment for robustness checks",
    "Binary high-use treatment",
    "Indicator for having a social-media account/profile",
    "Baseline loneliness outcome",
    "Raw scale already has higher values = worse wellbeing",
    "Raw scale already has higher values = worse wellbeing",
    "Raw scale already has higher values = worse wellbeing",
    "Control for number of close friends",
    "Log-transformed control used in the regression tables"
  ),
  stringsAsFactors = FALSE
)

#  store the retained baseline covariates and how they enter the cleaned dataset
baseline_covariate_reference <- data.frame(
  source_variable = c(
    "age_dv",
    "sex_dv",
    "ethn_dv",
    "country",
    "urban_dv",
    "gor_dv",
    "npns_dv",
    "ngrp_dv",
    "nnssib_dv",
    "ypnpal",
    "ypeatlivu",
    "yptvvidhrs",
    paste0("ypdevice", 1:6)
  ),
  source_label = c(
    get_var_label(first_available_wave[["age_dv"]], "age_dv"),
    get_var_label(first_available_wave[["sex_dv"]], "sex_dv"),
    get_var_label(first_available_wave[["ethn_dv"]], "ethn_dv"),
    get_var_label(first_available_wave[["country"]], "country"),
    get_var_label(first_available_wave[["urban_dv"]], "urban_dv"),
    get_var_label(first_available_wave[["gor_dv"]], "gor_dv"),
    get_var_label(first_available_wave[["npns_dv"]], "npns_dv"),
    get_var_label(first_available_wave[["ngrp_dv"]], "ngrp_dv"),
    get_var_label(first_available_wave[["nnssib_dv"]], "nnssib_dv"),
    get_var_label(first_available_wave[["ypnpal"]], "ypnpal"),
    get_var_label(first_available_wave[["ypeatlivu"]], "ypeatlivu"),
    get_var_label(first_available_wave[["yptvvidhrs"]], "yptvvidhrs"),
    vapply(paste0("ypdevice", 1:6), function(var) {
      get_var_label(first_available_wave[[var]], var)
    }, character(1))
  ),
  analysis_use = c(
    "age_dv",
    "sex_f",
    "ethnicity_broad_f",
    "country_f (retained separately from the rich main regression)",
    "urban_f",
    "region_f",
    "npns_dv",
    "ngrp_dv",
    "nnssib_dv",
    "close_friends and close_friends_log",
    "ypeatlivu",
    "yptvvidhrs",
    paste0("ypdevice", 1:6)
  ),
  transformation = c(
    "copy age_dv",
    "factor: 1 = Male, 2 = Female",
    "collapse detailed ethnicity into broad groups",
    "factor: 1 = England, 2 = Wales, 3 = Scotland, 4 = Northern Ireland",
    "factor: 1 = Urban, 2 = Rural",
    "factor from government office region code",
    "copy npns_dv",
    "copy ngrp_dv",
    "copy nnssib_dv",
    "copy ypnpal, then set values > 100 to NA and take log1p for the regression control",
    "copy ypeatlivu",
    "copy yptvvidhrs",
    rep("keep 0/1 mentioned indicator", 6)
  ),
  why_used = c(
    "Baseline demographic control",
    "Baseline demographic control",
    "Baseline demographic control",
    "Retained national geography reference",
    "Baseline demographic control",
    "Baseline geographic control",
    "Baseline household-composition control",
    "Baseline household-composition control",
    "Baseline household-composition control",
    "Baseline social control",
    "Baseline family-environment control",
    "Baseline screen-time control",
    rep("Baseline device control", 6)
  ),
  stringsAsFactors = FALSE
)

cat("\nStep 5/8: record the analyst-defined transformations after inspecting the raw coding\n")
cat("These transformation choices do not come directly from the data; they are the explicit construction decisions used in the thesis.\n")
cat("\nMain constructed treatment and outcome variables:\n")
print(derived_variable_spec, row.names = FALSE)
cat("\nBaseline covariates retained in the Chapter 3 dataset:\n")
cat("The rich main regression later uses region_f instead of entering country_f and region_f jointly.\n")
print(baseline_covariate_reference, row.names = FALSE)

#  move the main analysis variables to the front of the cleaned dataset
front_cols <- c(
  "pidp",
  "wave",
  "wave_number",
  "wave_f",
  "social_media_weekday",
  "social_media_weekend",
  "high_social_media_weekday",
  "has_social_profile",
  "loneliness",
  "life_dissatisfaction",
  "schoolwork_dissatisfaction",
  "school_dissatisfaction",
  "close_friends",
  "close_friends_log"
)
youth <- youth[c(front_cols, setdiff(names(youth), front_cols))]

#  save the cleaned person-wave dataset used by the later scripts
analysis_path_rds <- file.path(out_dir, paste0("ukhls_youth_", waves[1], "_to_", tail(waves, 1), "_clean.rds"))
saveRDS(youth, analysis_path_rds)

#  list the baseline outcomes used in the descriptive summaries
outcomes <- c(
  "loneliness",
  "life_dissatisfaction",
  "schoolwork_dissatisfaction",
  "school_dissatisfaction"
)

#  count usable observations and treatment coverage by wave
summary_by_wave <- youth %>%
  group_by(wave) %>%
  summarise(
    n_rows = n(),
    n_persons = n_distinct(pidp),
    n_treatment = sum(!is.na(social_media_weekday)),
    mean_treatment = mean(social_media_weekday, na.rm = TRUE),
    across(all_of(outcomes), ~ sum(!is.na(.x)), .names = "n_{.col}"),
    .groups = "drop"
  )

#  count how many waves each adolescent appears in
panel_counts <- youth %>%
  count(pidp, name = "n_waves") %>%
  count(n_waves, name = "n_adolescents") %>%
  mutate(share = n_adolescents / sum(n_adolescents))

#  record the person-wave and unique-person sample sizes by wave
wave_sample_sizes <- youth %>%
  group_by(wave) %>%
  summarise(
    n_person_wave_rows = n(),
    n_unique_adolescents = n_distinct(pidp),
    .groups = "drop"
  )

cat("\nStep 6/8: save the cleaned baseline dataset and summarise the selected sample\n")
cat(
  nrow(youth),
  "person-wave observations for",
  n_distinct(youth$pidp),
  "unique adolescents.\n\n"
)
print(wave_sample_sizes)
print(panel_counts)
print(summary_by_wave)

#  create a working copy for the descriptive tables and figures
dt <- youth

#  build simple binary indicators used in the descriptives
dt$female <- ifelse(is.na(dt$sex_dv), NA_real_, as.numeric(dt$sex_dv == 2))
dt$urban <- ifelse(is.na(dt$urban_dv), NA_real_, as.numeric(dt$urban_dv == 1))

#  define the rows shown in the descriptive-statistics table
desc_vars <- data.frame(
  section = c(
    rep("Sample characteristics", 3),
    rep("Social-media exposure", 4),
    rep("Mental wellbeing", 4),
    rep("Baseline controls", 10)
  ),
  variable = c(
    "Age",
    "Female",
    "Urban area",
    "Has social-media account/profile",
    "Weekday social-media use",
    "Weekend social-media interaction",
    "High weekday social-media use",
    "Loneliness",
    "Life dissatisfaction",
    "School-work dissatisfaction",
    "School dissatisfaction",
    "Number of close friends",
    "Parents/step-parents in household",
    "Grandparents in household",
    "Siblings/step-siblings in household",
    "Family meals",
    "Television/video hours",
    "Smartphone",
    "Tablet",
    "Gaming console",
    "Laptop/desktop computer"
  ),
  range = c(
    "9--16",
    "0--1",
    "0--1",
    "0--1",
    "1--5",
    "1--5",
    "0--1",
    "1--3",
    "1--7",
    "1--7",
    "1--7",
    "0--100",
    "0--2",
    "0--4",
    "0--20",
    "1--4",
    "1--5",
    "0--1",
    "0--1",
    "0--1",
    "0--1"
  ),
  source = c(
    "age_dv",
    "female",
    "urban",
    "has_social_profile",
    "social_media_weekday",
    "social_media_weekend",
    "high_social_media_weekday",
    "loneliness",
    "life_dissatisfaction",
    "schoolwork_dissatisfaction",
    "school_dissatisfaction",
    "close_friends",
    "npns_dv",
    "ngrp_dv",
    "nnssib_dv",
    "ypeatlivu",
    "yptvvidhrs",
    "ypdevice1",
    "ypdevice3",
    "ypdevice5",
    "ypdevice6"
  ),
  stringsAsFactors = FALSE
)

#  summarise one variable for the pooled sample
summarise_one <- function(source) {
  x <- dt[[source]]
  x <- x[!is.na(x)]

  data.frame(
    N = length(x),
    Mean = ifelse(length(x) == 0, NA_real_, mean(x)),
    SD = ifelse(length(x) <= 1, NA_real_, stats::sd(x)),
    stringsAsFactors = FALSE
  )
}

desc_simple <- bind_cols(
  desc_vars,
  bind_rows(lapply(desc_vars$source, summarise_one))
) %>%
  mutate(
    variable_with_range = paste0(variable, " (", range, ")"),
    N = format_number(N, digits = 0),
    Mean = format_number(Mean),
    SD = format_number(SD)
  ) %>%
  select(section, variable_with_range, N, Mean, SD)

wave_rows <- data.frame(
  statistic = c(
    "Wave l",
    "Wave m",
    "Wave n",
    "Wave o",
    "Pooled person-wave observations"
  ),
  count = c(
    wave_sample_sizes$n_person_wave_rows[match(c("l", "m", "n", "o"), wave_sample_sizes$wave)],
    nrow(dt)
  ),
  stringsAsFactors = FALSE
)

panel_rows <- panel_counts %>%
  mutate(
    statistic = paste0("Observed in ", n_waves, ifelse(n_waves == 1, " wave", " waves")),
    count = n_adolescents
  ) %>%
  select(statistic, count) %>%
  bind_rows(
    data.frame(
      statistic = "Unique adolescents",
      count = n_distinct(dt$pidp),
      stringsAsFactors = FALSE
    )
  )

write_sample_composition_latex <- function(wave_tab, panel_tab, file) {
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Sample composition of the baseline panel}",
    "\\label{tab:sample_composition}",
    "\\small",
    "\\renewcommand{\\arraystretch}{1.12}",
    "\\begin{tabular}{lr}",
    "\\toprule",
    "Statistic & N \\\\",
    "\\midrule",
    "\\multicolumn{2}{l}{\\textit{Panel A: Wave sample sizes}} \\\\"
  )

  for (i in seq_len(nrow(wave_tab))) {
    if (wave_tab$statistic[i] == "Pooled person-wave observations") {
      lines <- c(lines, "\\cmidrule(lr){1-2}")
    }
    lines <- c(
      lines,
      paste0(
        paste(
          escape_latex(wave_tab$statistic[i]),
          format_number(wave_tab$count[i], digits = 0),
          sep = " & "
        ),
        " \\\\"
      )
    )
  }

  lines <- c(lines, "\\addlinespace", "\\multicolumn{2}{l}{\\textit{Panel B: Number of observed waves per adolescent}} \\\\")

  for (i in seq_len(nrow(panel_tab))) {
    if (panel_tab$statistic[i] == "Unique adolescents") {
      lines <- c(lines, "\\cmidrule(lr){1-2}")
    }
    lines <- c(
      lines,
      paste0(
        paste(
          escape_latex(panel_tab$statistic[i]),
          format_number(panel_tab$count[i], digits = 0),
          sep = " & "
        ),
        " \\\\"
      )
    )
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{minipage}{0.9\\textwidth}",
    "\\footnotesize Notes: The sample pools youth observations from Waves \\texttt{l}--\\texttt{o}. The pooled person-wave count equals the sum of the four wave-specific counts.",
    "\\end{minipage}",
    "\\end{table}"
  )

  writeLines(lines, file)
}

write_descriptive_latex <- function(desc, file) {
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Descriptive statistics for the main analysis variables}",
    "\\label{tab:descriptive_statistics}",
    "\\small",
    "\\renewcommand{\\arraystretch}{1.12}",
    "\\begin{tabular}{p{0.62\\textwidth}rrr}",
    "\\toprule",
    "Variable (possible range) & N & Mean & SD \\\\",
    "\\midrule"
  )

  for (section in unique(desc$section)) {
    lines <- c(lines, paste0("\\multicolumn{4}{l}{\\textit{", escape_latex(section), "}} \\\\"))
    section_rows <- desc[desc$section == section, ]

    for (i in seq_len(nrow(section_rows))) {
      row <- section_rows[i, ]
      lines <- c(
        lines,
        paste0(
          paste(
            escape_latex(row$variable_with_range),
            row$N,
            row$Mean,
            row$SD,
            sep = " & "
          ),
          " \\\\"
        )
      )
    }

    if (section != tail(unique(desc$section), 1)) {
      lines <- c(lines, "\\addlinespace")
    }
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{minipage}{0.9\\textwidth}",
    "\\footnotesize Notes: The sample stacks youth observations from Waves \\texttt{l}--\\texttt{o}. Higher values of the recoded wellbeing outcomes indicate worse wellbeing. Binary variables are coded 0/1, so their means are proportions. Reported means and standard deviations use the raw variables; \\texttt{close\\_friends\\_log} is used only in the regression control set.",
    "\\end{minipage}",
    "\\end{table}"
  )

  writeLines(lines, file)
}

cat("\nStep 7/8: generate the Chapter 3 descriptive tables\n")
print(desc_simple)
write_sample_composition_latex(
  wave_rows,
  panel_rows,
  file.path(tables_dir, "table_sample_composition.tex")
)
write_descriptive_latex(
  desc_simple,
  file.path(tables_dir, "table1_descriptive_statistics.tex")
)

#  define readable labels for the weekday-use categories in the figures
use_labels <- c(
  `1` = "None",
  `2` = "<1 hour",
  `3` = "1-3 hours",
  `4` = "4-6 hours",
  `5` = "7+ hours"
)

#  build the shared plotting
dt$use_cat <- factor(use_labels[as.character(dt$social_media_weekday)], levels = unname(use_labels))
dt$wave_display <- factor(dt$wave, levels = c("l", "m", "n", "o"), labels = c("Wave l", "Wave m", "Wave n", "Wave o"))

base_theme <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

use_palette <- c(
  "None" = "#f8d7e3",
  "<1 hour" = "#efc7e8",
  "1-3 hours" = "#dcc8f2",
  "4-6 hours" = "#c9c3f4",
  "7+ hours" = "#b8c7f7"
)

presentation_use_labels <- c(
  `1` = "None",
  `2` = "<1h",
  `3` = "1-3h",
  `4` = "4-6h",
  `5` = "7+h"
)

dt$use_cat_presentation <- factor(
  presentation_use_labels[as.character(dt$social_media_weekday)],
  levels = unname(presentation_use_labels)
)

presentation_use_palette <- c(
  "None" = "#FFE3F1",
  "<1h" = "#F1C9FF",
  "1-3h" = "#D7B6FF",
  "4-6h" = "#A783FF",
  "7+h" = "#FF4FB3"
)

presentation_theme <- theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", color = "#352044", size = 16),
    legend.text = element_text(color = "#352044", size = 15),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E7E0EE", linewidth = 0.55),
    strip.text = element_text(face = "bold", size = 18, color = "#352044"),
    axis.text = element_text(color = "#352044", size = 15),
    axis.title = element_text(face = "bold", color = "#352044", size = 17),
    panel.spacing = grid::unit(1.1, "lines"),
    plot.margin = margin(6, 10, 4, 6)
  )

save_presentation_figure <- function(file_name, plot, width, height) {
  ggsave(
    filename = file.path(figures_dir, file_name),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 320,
    bg = "white"
  )
}

#  prepare and save the distribution-of-use figure by wave
distribution_df <- dt %>%
  filter(!is.na(use_cat), !is.na(wave_display)) %>%
  count(wave_display, use_cat) %>%
  group_by(wave_display) %>%
  mutate(share = n / sum(n))

fig_distribution <- ggplot(distribution_df, aes(x = wave_display, y = share, fill = use_cat)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = use_palette, name = "Weekday social-media use") +
  labs(
    x = NULL,
    y = "Share within wave"
  ) +
  base_theme

ggsave(
  filename = file.path(figures_dir, "descriptive_treatment_distribution_by_wave.png"),
  plot = fig_distribution,
  width = 8.4,
  height = 5.2,
  dpi = 300
)

distribution_presentation_df <- dt %>%
  filter(!is.na(use_cat_presentation), !is.na(wave_display)) %>%
  count(wave_display, use_cat_presentation) %>%
  group_by(wave_display) %>%
  mutate(share = n / sum(n))

fig_distribution_presentation <- ggplot(
  distribution_presentation_df,
  aes(x = wave_display, y = share, fill = use_cat_presentation)
) +
  geom_col(width = 0.74, color = "white", linewidth = 0.35) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.04))
  ) +
  scale_fill_manual(values = presentation_use_palette, name = "Weekday use") +
  labs(x = NULL, y = "Share within wave") +
  presentation_theme

save_presentation_figure(
  "descriptive_treatment_distribution_by_wave_presentation.png",
  fig_distribution_presentation,
  width = 10.4,
  height = 6.0
)

#  prepare and save the high-use-by-age-and-sex figure
age_sex_df <- dt %>%
  filter(
    !is.na(age_dv),
    age_dv >= 10,
    age_dv <= 15,
    !is.na(sex_f),
    !is.na(high_social_media_weekday)
  ) %>%
  group_by(age_dv, sex_f) %>%
  summarise(
    N = n(),
    share_high_use = mean(high_social_media_weekday),
    se = sqrt(share_high_use * (1 - share_high_use) / N),
    lower = pmax(0, share_high_use - 1.96 * se),
    upper = pmin(1, share_high_use + 1.96 * se),
    .groups = "drop"
  )

sex_line_palette <- c(
  "Male" = "#8f82d8",
  "Female" = "#de8fb8"
)

sex_fill_palette <- c(
  "Male" = "#dcd5f8",
  "Female" = "#f8d7e3"
)

fig_high_use_age_sex <- ggplot(
  age_sex_df,
  aes(x = age_dv, y = share_high_use, color = sex_f, fill = sex_f, group = sex_f)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.28, linewidth = 0) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  scale_x_continuous(breaks = sort(unique(age_sex_df$age_dv))) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_color_manual(values = sex_line_palette, name = NULL) +
  scale_fill_manual(values = sex_fill_palette, name = NULL) +
  labs(
    x = "Age",
    y = "Share using social media 4+ hours"
  ) +
  base_theme

ggsave(
  filename = file.path(figures_dir, "descriptive_high_use_by_age_and_sex.png"),
  plot = fig_high_use_age_sex,
  width = 8.2,
  height = 5.2,
  dpi = 300
)

presentation_sex_line_palette <- c(
  "Male" = "#8C1AFF",
  "Female" = "#FF4FB3"
)

presentation_sex_fill_palette <- c(
  "Male" = "#D7B6FF",
  "Female" = "#FFD0E8"
)

fig_high_use_age_sex_presentation <- ggplot(
  age_sex_df,
  aes(x = age_dv, y = share_high_use, color = sex_f, fill = sex_f, group = sex_f)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.24, linewidth = 0) +
  geom_line(linewidth = 1.35) +
  geom_point(size = 3.2) +
  scale_x_continuous(breaks = sort(unique(age_sex_df$age_dv))) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.03, 0.10))
  ) +
  scale_color_manual(values = presentation_sex_line_palette, name = NULL) +
  scale_fill_manual(values = presentation_sex_fill_palette, name = NULL) +
  labs(x = "Age", y = "Share using social media 4+ hours") +
  presentation_theme +
  theme(legend.position = "top")

save_presentation_figure(
  "descriptive_high_use_by_age_and_sex_presentation.png",
  fig_high_use_age_sex_presentation,
  width = 10.0,
  height = 6.2
)

#  prepare and save the high-use-by-controls figure
family_meals_labels <- c(
  `1` = "Hardly ever",
  `2` = "Less than weekly",
  `3` = "Several times weekly",
  `4` = "Every day"
)

tv_hours_labels <- c(
  `1` = "None",
  `2` = "<1 hour",
  `3` = "1-3 hours",
  `4` = "4-6 hours",
  `5` = "7+ hours"
)

high_use_controls_df <- bind_rows(
  dt %>%
    filter(!is.na(ypeatlivu), !is.na(high_social_media_weekday)) %>%
    mutate(
      control = "Family meals",
      control_group = factor(
        family_meals_labels[as.character(ypeatlivu)],
        levels = unname(family_meals_labels)
      )
    ) %>%
    group_by(control, control_group) %>%
    summarise(
      N = n(),
      share_high_use = mean(high_social_media_weekday),
      se = sqrt(share_high_use * (1 - share_high_use) / N),
      lower = pmax(0, share_high_use - 1.96 * se),
      upper = pmin(1, share_high_use + 1.96 * se),
      .groups = "drop"
    ),
  dt %>%
    filter(!is.na(yptvvidhrs), !is.na(high_social_media_weekday)) %>%
    mutate(
      control = "Television/video time",
      control_group = factor(
        tv_hours_labels[as.character(yptvvidhrs)],
        levels = unname(tv_hours_labels)
      )
    ) %>%
    group_by(control, control_group) %>%
    summarise(
      N = n(),
      share_high_use = mean(high_social_media_weekday),
      se = sqrt(share_high_use * (1 - share_high_use) / N),
      lower = pmax(0, share_high_use - 1.96 * se),
      upper = pmin(1, share_high_use + 1.96 * se),
      .groups = "drop"
    )
)

fig_high_use_controls <- ggplot(
  high_use_controls_df,
  aes(x = control_group, y = share_high_use, group = 1)
) +
  geom_col(width = 0.72, fill = "#dcd5f8", color = "#8f82d8", linewidth = 0.25) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12, color = "#5f4fa3", linewidth = 0.55) +
  facet_wrap(~ control, scales = "free_x", ncol = 2) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    x = NULL,
    y = "Share using social media 4+ hours"
  ) +
  base_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

ggsave(
  filename = file.path(figures_dir, "descriptive_high_use_by_controls.png"),
  plot = fig_high_use_controls,
  width = 9.2,
  height = 5.6,
  dpi = 300
)

high_use_controls_presentation_df <- high_use_controls_df %>%
  mutate(
    control_raw = as.character(control),
    control = factor(
      ifelse(control_raw == "Television/video time", "TV/video time", control_raw),
      levels = c("Family meals", "TV/video time")
    ),
    control_group = factor(
      case_when(
        control_raw == "Family meals" & as.character(control_group) == "Hardly ever" ~ "Hardly ever",
        control_raw == "Family meals" & as.character(control_group) == "Less than weekly" ~ "< weekly",
        control_raw == "Family meals" & as.character(control_group) == "Several times weekly" ~ "Several\n/week",
        control_raw == "Family meals" & as.character(control_group) == "Every day" ~ "Daily",
        control_raw == "Television/video time" & as.character(control_group) == "<1 hour" ~ "<1h",
        control_raw == "Television/video time" & as.character(control_group) == "1-3 hours" ~ "1-3h",
        control_raw == "Television/video time" & as.character(control_group) == "4-6 hours" ~ "4-6h",
        control_raw == "Television/video time" & as.character(control_group) == "7+ hours" ~ "7+h",
        TRUE ~ as.character(control_group)
      ),
      levels = c(
        "Hardly ever", "< weekly", "Several\n/week", "Daily",
        "None", "<1h", "1-3h", "4-6h", "7+h"
      )
    )
  )

fig_high_use_controls_presentation <- ggplot(
  high_use_controls_presentation_df,
  aes(x = control_group, y = share_high_use, group = 1)
) +
  geom_col(width = 0.72, fill = "#D7B6FF", color = "#8C1AFF", linewidth = 0.45) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.14, color = "#39224F", linewidth = 0.75) +
  facet_wrap(~ control, scales = "free_x", ncol = 2) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  labs(x = NULL, y = "Share using social media 4+ hours") +
  presentation_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 13.5)
  )

save_presentation_figure(
  "descriptive_high_use_by_controls_presentation.png",
  fig_high_use_controls_presentation,
  width = 10.4,
  height = 5.8
)

#  prepare and save the high-use-by-household-composition figure
high_use_household_df <- bind_rows(
  dt %>%
    filter(!is.na(npns_dv), !is.na(high_social_media_weekday)) %>%
    mutate(
      control = "Parent figures in household",
      control_group = factor(
        as.character(npns_dv),
        levels = c("0", "1", "2")
      )
    ) %>%
    group_by(control, control_group) %>%
    summarise(
      N = n(),
      share_high_use = mean(high_social_media_weekday),
      se = sqrt(share_high_use * (1 - share_high_use) / N),
      lower = pmax(0, share_high_use - 1.96 * se),
      upper = pmin(1, share_high_use + 1.96 * se),
      .groups = "drop"
    ),
  dt %>%
    filter(!is.na(ngrp_dv), !is.na(high_social_media_weekday)) %>%
    mutate(
      control = "Grandparent in household",
      control_group = factor(
        ifelse(ngrp_dv > 0, "Yes", "No"),
        levels = c("No", "Yes")
      )
    ) %>%
    group_by(control, control_group) %>%
    summarise(
      N = n(),
      share_high_use = mean(high_social_media_weekday),
      se = sqrt(share_high_use * (1 - share_high_use) / N),
      lower = pmax(0, share_high_use - 1.96 * se),
      upper = pmin(1, share_high_use + 1.96 * se),
      .groups = "drop"
    ),
  dt %>%
    filter(!is.na(nnssib_dv), !is.na(high_social_media_weekday)) %>%
    mutate(
      control = "Siblings in household",
      control_group = factor(
        case_when(
          nnssib_dv == 0 ~ "0",
          nnssib_dv == 1 ~ "1",
          nnssib_dv == 2 ~ "2",
          nnssib_dv >= 3 ~ "3+"
        ),
        levels = c("0", "1", "2", "3+")
      )
    ) %>%
    group_by(control, control_group) %>%
    summarise(
      N = n(),
      share_high_use = mean(high_social_media_weekday),
      se = sqrt(share_high_use * (1 - share_high_use) / N),
      lower = pmax(0, share_high_use - 1.96 * se),
      upper = pmin(1, share_high_use + 1.96 * se),
      .groups = "drop"
    )
) %>%
  mutate(
    control = factor(
      control,
      levels = c(
        "Parent figures in household",
        "Grandparent in household",
        "Siblings in household"
      )
    )
  )

fig_high_use_household <- ggplot(
  high_use_household_df,
  aes(x = control_group, y = share_high_use, group = 1)
) +
  geom_col(width = 0.72, fill = "#dcd5f8", color = "#8f82d8", linewidth = 0.25) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12, color = "#5f4fa3", linewidth = 0.55) +
  facet_wrap(~ control, scales = "free_x", ncol = 3) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    x = NULL,
    y = "Share using social media 4+ hours"
  ) +
  base_theme +
  theme(legend.position = "none")

ggsave(
  filename = file.path(figures_dir, "descriptive_high_use_by_household.png"),
  plot = fig_high_use_household,
  width = 10.2,
  height = 4.8,
  dpi = 300
)

high_use_household_presentation_df <- high_use_household_df %>%
  mutate(
    control = factor(
      as.character(control),
      levels = c(
        "Parent figures in household",
        "Grandparent in household",
        "Siblings in household"
      ),
      labels = c("Parent figures", "Grandparent", "Siblings")
    )
  )

fig_high_use_household_presentation <- ggplot(
  high_use_household_presentation_df,
  aes(x = control_group, y = share_high_use, group = 1)
) +
  geom_col(width = 0.72, fill = "#D7B6FF", color = "#8C1AFF", linewidth = 0.45) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.14, color = "#39224F", linewidth = 0.75) +
  facet_wrap(~ control, scales = "free_x", ncol = 3) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  labs(x = NULL, y = "Share using social media 4+ hours") +
  presentation_theme +
  theme(legend.position = "none")

save_presentation_figure(
  "descriptive_high_use_by_household_presentation.png",
  fig_high_use_household_presentation,
  width = 10.4,
  height = 5.0
)

#  define the outcome labels used in the final descriptive figure
outcome_map <- c(
  loneliness = "Loneliness",
  life_dissatisfaction = "Life dissatisfaction",
  schoolwork_dissatisfaction = "School-work dissatisfaction",
  school_dissatisfaction = "School dissatisfaction"
)

#  prepare and save the outcome-means-by-treatment figure
outcome_df <- bind_rows(lapply(names(outcome_map), function(var) {
  data.frame(
    use_cat = dt$use_cat,
    outcome = outcome_map[[var]],
    value = dt[[var]],
    stringsAsFactors = FALSE
  )
})) %>%
  filter(!is.na(use_cat), !is.na(value)) %>%
  group_by(outcome, use_cat) %>%
  summarise(
    N = n(),
    mean_value = mean(value),
    se = sd(value) / sqrt(N),
    lower = mean_value - 1.96 * se,
    upper = mean_value + 1.96 * se,
    .groups = "drop"
  )

fig_outcomes <- ggplot(outcome_df, aes(x = use_cat, y = mean_value, group = 1)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#dcc8f2", alpha = 0.35) +
  geom_line(color = "#8f82d8", linewidth = 0.9) +
  geom_point(color = "#de8fb8", size = 2.1) +
  facet_wrap(~ outcome, scales = "free_y", ncol = 2) +
  labs(
    x = NULL,
    y = "Mean outcome"
  ) +
  base_theme +
  theme(legend.position = "none")

ggsave(
  filename = file.path(figures_dir, "descriptive_outcomes_by_treatment.png"),
  plot = fig_outcomes,
  width = 9.4,
  height = 6.8,
  dpi = 300
)

outcome_presentation_df <- outcome_df %>%
  mutate(
    outcome = factor(
      outcome,
      levels = c(
        "Life dissatisfaction",
        "Loneliness",
        "School dissatisfaction",
        "School-work dissatisfaction"
      )
    ),
    use_cat_presentation = factor(
      case_when(
        as.character(use_cat) == "<1 hour" ~ "<1h",
        as.character(use_cat) == "1-3 hours" ~ "1-3h",
        as.character(use_cat) == "4-6 hours" ~ "4-6h",
        as.character(use_cat) == "7+ hours" ~ "7+h",
        TRUE ~ as.character(use_cat)
      ),
      levels = unname(presentation_use_labels)
    )
  )

fig_outcomes_presentation <- ggplot(
  outcome_presentation_df,
  aes(x = use_cat_presentation, y = mean_value, group = 1)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "#D7B6FF",
    alpha = 0.28
  ) +
  geom_line(color = "#8C1AFF", linewidth = 1.35) +
  geom_point(color = "#FF4FB3", size = 3.2) +
  facet_wrap(~ outcome, scales = "free_y", ncol = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.22))) +
  labs(x = NULL, y = "Mean outcome") +
  presentation_theme +
  theme(legend.position = "none")

save_presentation_figure(
  "descriptive_outcomes_by_treatment_presentation.png",
  fig_outcomes_presentation,
  width = 10.4,
  height = 6.9
)

cat("\nStep 8/8: generate the Chapter 3 descriptive figures\n")
cat("Saved:\n")
cat("  -", file.path(figures_dir, "descriptive_treatment_distribution_by_wave.png"), "\n")
cat("  -", file.path(figures_dir, "descriptive_treatment_distribution_by_wave_presentation.png"), "\n")
cat("  -", file.path(figures_dir, "descriptive_high_use_by_age_and_sex.png"), "\n")
cat("  -", file.path(figures_dir, "descriptive_high_use_by_age_and_sex_presentation.png"), "\n")
cat("  -", file.path(figures_dir, "descriptive_high_use_by_controls.png"), "\n")
cat("  -", file.path(figures_dir, "descriptive_high_use_by_controls_presentation.png"), "\n")
cat("  -", file.path(figures_dir, "descriptive_high_use_by_household.png"), "\n")
cat("  -", file.path(figures_dir, "descriptive_high_use_by_household_presentation.png"), "\n")
cat("  -", file.path(figures_dir, "descriptive_outcomes_by_treatment.png"), "\n")
cat("  -", file.path(figures_dir, "descriptive_outcomes_by_treatment_presentation.png"), "\n")
cat("  -", file.path(tables_dir, "table_sample_composition.tex"), "\n")
cat("  -", file.path(tables_dir, "table1_descriptive_statistics.tex"), "\n")

message("Saved clean data to: ", analysis_path_rds)
