# Chapter 3 

# Summary for Chapter 3:
# 1. This script identifies which UKHLS youth waves contain all treatment,
#    outcome, and baseline-control variables needed for the analysis.
# 2. It reads the selected youth files, cleans special missing codes, and keeps
#    only valid survey response categories based on the UKHLS value labels.
# 3. It constructs the main analysis variables used in Chapter 3:
#    social-media treatment measures, wellbeing outcomes, and core controls.
# 4. It then shows the main descriptive tables and figures in the
#    console and plot pane.

library(haven)
library(dplyr)
library(ggplot2)
library(scales)

# 0. Setup -----------------------------------------------------------------
# Raw UKHLS youth files.
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

all_waves <- letters[1:15]
wave_numbers <- setNames(seq_along(all_waves), all_waves)

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

# 1. Variables used in Chapter 3 -------------------------------------------
# These are the treatment, outcome, and control variables
treatment_vars <- c("ypnetcht", "ypnetchtw", "ypsocweb")
outcome_vars <- c("yplonely", "yphlf", "yphsw", "yphsc")
baseline_control_vars <- c(
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
)

# These variables are not in the baseline control set, but they matter for the
# wave-availability discussion in Section 3.2 and for possible extensions.?
additional_social_vars <- c("ypvirfnd", "ypfndmeet", "ypfndonl")

required_vars <- c(treatment_vars, outcome_vars, baseline_control_vars, additional_social_vars)
raw_keep_vars <- c("psu", "strata", "ythscui_xw", required_vars)

# 2. Wave availability check behind Section 3.2 ----------------------------
# Helper: read one raw variable directly from one wave.
get_raw_var <- function(wave, var) {
  raw_name <- paste0(wave, "_", var)
  read_dta(file.path(root, paste0(wave, "_youth.dta")), col_select = all_of(raw_name))[[1]]
}

# Helper: pull the questionnaire label shown in the UKHLS files.
get_var_label <- function(wave, var) {
  label <- attr(get_raw_var(wave, var), "label")
  if (is.null(label)) "" else label
}

# Check in which waves each required variable is available.
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

first_available_wave <- setNames(
  vapply(required_vars, function(var) {
    min(wave_availability$wave[wave_availability$variable == var & wave_availability$available])
  }, character(1)),
  required_vars
)

wave_choice_summary <- data.frame(
  variable = required_vars,
  questionnaire_item = vapply(required_vars, function(var) {
    get_var_label(first_available_wave[[var]], var)
  }, character(1)),
  first_available_wave = unname(first_available_wave),
  stringsAsFactors = FALSE
)

# Compare all consecutive wave blocks ending in wave o.
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

selected_start_wave <- candidate_wave_blocks$start_wave[candidate_wave_blocks$all_required_vars_available][1]
if (is.na(selected_start_wave)) {
  stop("No feasible consecutive wave block found for the selected Chapter 3 variables.", call. = FALSE)
}

waves <- all_waves[match(selected_start_wave, all_waves):length(all_waves)]

show_table("Wave choice summary", wave_choice_summary)
show_table("Candidate wave blocks", candidate_wave_blocks)

cat("Selected baseline block:", paste(waves, collapse = ", "), "\n")

# 3. Read and clean the raw youth files ------------------------------------
# Helper: convert labelled fields to plain numeric values and drop negative special codes.
clean_numeric <- function(x) {
  x <- as.numeric(zap_labels(x))
  x[x < 0] <- NA_real_
  x
}

# Helper: keep only valid non-special codes for each raw survey item.
get_valid_codes <- function(var) {
  first_wave <- first_available_wave[[var]]
  x <- get_raw_var(first_wave, var)
  labels <- attr(x, "labels")

  if (is.null(labels)) {
    return(NULL)
  }

  labels <- labels[order(unname(labels))]
  special_names <- c("missing", "inapplicable", "refusal", "don't know", "inconsistent")
  special_codes <- unname(labels[tolower(names(labels)) %in% special_names])

  observed_codes <- sort(unique(as.numeric(zap_labels(x))))
  observed_codes <- observed_codes[!is.na(observed_codes)]
  as.integer(observed_codes[!(observed_codes %in% special_codes)])
}

valid_codes <- setNames(lapply(required_vars, get_valid_codes), required_vars)

# Helper: read one wave, keep the selected columns, and remove the wave prefix.
# Explained NL: leest precies één wave in, bijvoorbeeld wave l. 
# read_dta(...) opent dan l_youth.dta, 
# keep <- c("pidp", paste0(wave, "_", raw_keep_vars)) maakt de lijst van kolommen die je wilt bewaren, 
# intersect(keep, names(df)) zorgt dat alleen bestaande kolommen worden gekozen, 
# sub(paste0("^", wave, "_"), "", names(df)) haalt de prefix weg zodat 
# l_age_dv gewoon age_dv wordt, en daarna voegt de functie nog wave en 
# wave_number toe. Het resultaat is dus: één nette data frame voor één wave, 
# klaar om later met bind_rows(...) onder elkaar te plakken
read_one_wave <- function(wave) {
  df <- as.data.frame(read_dta(file.path(root, paste0(wave, "_youth.dta"))))
  keep <- c("pidp", paste0(wave, "_", raw_keep_vars))
  df <- df[intersect(keep, names(df))]
  names(df) <- sub(paste0("^", wave, "_"), "", names(df))
  df$wave <- wave
  df$wave_number <- unname(wave_numbers[wave])
  df
}

# Stack waves l--o into one person-wave file.
youth <- bind_rows(lapply(waves, read_one_wave))

# Convert labelled survey fields to numeric values.
# Explained NL: names(youth) pakt alle kolomnamen.
# setdiff(..., "wave") haalt alleen de kolom wave eruit, 
# omdat die letters bevat zoals "l", "m", "n", "o" en dus niet 
# numeriek moet worden gemaakt. lapply(..., clean_numeric) past de functie 
# clean_numeric toe op alle andere kolommen. En clean_numeric deed eerder dit:
  # zap_labels(x): haalt de Stata/UKHLS labels eraf.
  # as.numeric(...): maakt de variabele een gewone numerieke vector.
  # x[x < 0] <- NA_real_: zet alle negatieve special codes om naar missing.
# Dus bijvoorbeeld:
#  -8 = “don’t know” wordt NA
#  -1 = “missing” wordt NA
# een echte antwoordcode zoals 3 blijft gewoon 3
numeric_cols <- setdiff(names(youth), "wave")
youth[numeric_cols] <- lapply(youth[numeric_cols], clean_numeric)

# Enforce the valid code ranges derived from the UKHLS value labels. (Extra check)
for (var in names(valid_codes)) {
  allowed <- valid_codes[[var]]
  if (is.null(allowed)) next

  bad <- !(youth[[var]] %in% allowed) & !is.na(youth[[var]])
  youth[[var]][bad] <- NA_real_
}

# 4. Construct treatment, outcomes, and controls for Section 3.3 -----------
# Treatment variables used in the baseline and robustness checks.
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

# Outcome variables: higher values always mean worse wellbeing.
youth$loneliness <- youth$yplonely
youth$life_dissatisfaction <- youth$yphlf
youth$schoolwork_dissatisfaction <- youth$yphsw
youth$school_dissatisfaction <- youth$yphsc

# Social controls: close friends is cleaned, then log-transformed for richer models
# so zero remains defined and the right tail has less leverage.
youth$close_friends <- ifelse(youth$ypnpal > 100, NA_real_, youth$ypnpal)
youth$close_friends_log <- log1p(youth$close_friends)

# Make factor versions of the main categorical controls.
# Explained NL: zet numerieke surveycodes om naar leesbare categorieën. 
# Als x bijvoorbeeld 1 2 2 1 is en labels is c(`1` = "Male", `2` = "Female"), 
# dan maakt labels[as.character(x)] daar "Male" "Female" "Female" "Male" van. 
# factor(..., levels = unname(labels)) maakt er daarna een factor van met precies die volgorde van categorieën.
factor_from_codes <- function(x, labels) {
  factor(labels[as.character(x)], levels = unname(labels))
}

# Explained NL: youth$wave bevat nu letters zoals "l", "m", "n", "o"
# factor(...) zegt: behandel dit niet als gewone tekst, maar als categorieën
# levels = waves bepaalt de volgorde van die categorieën
# labels = paste("Wave", toupper(waves)) geeft mooiere namen: "l" wordt "Wave L"
# Dus als waves = c("l","m","n","o"), dan krijg je: ruwe waarde: l, factorlabel: Wave L
youth$wave_f <- factor(youth$wave, levels = waves, labels = paste("Wave", toupper(waves)))
youth$sex_f <- factor_from_codes(youth$sex_dv, c(`1` = "Male", `2` = "Female"))
youth$country_f <- factor_from_codes(
  youth$country,
  c(`1` = "England", `2` = "Wales", `3` = "Scotland", `4` = "Northern Ireland")
)
youth$urban_f <- factor_from_codes(youth$urban_dv, c(`1` = "Urban", `2` = "Rural"))
youth$region_f <- factor(youth$gor_dv)
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

# 5. Organise the cleaned Chapter 3 panel ----------------------------------
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

# Explained: front_cols is een lijst met de belangrijkste variabelen die vooraan wilt hebben, 
# zoals pidp, wave, treatment en outcomes,
# setdiff(names(youth), front_cols) pakt alle andere kolommen die nog over zijn
# c(front_cols, ...) plakt die twee lijsten achter elkaar:
# eerst de belangrijke kolommen, daarna alle rest
# youth[...] herschikt dan de dataset in die nieuwe kolomvolgorde
youth <- youth[c(front_cols, setdiff(names(youth), front_cols))]

cat("\nCleaned Chapter 3 panel dimensions:", nrow(youth), "rows x", ncol(youth), "columns\n")
show_table(
  "Cleaned Chapter 3 panel preview (first 10 rows)",
  utils::head(youth[, front_cols], 10)
)

# 6. Plain outputs for Chapter 3 tables ------------------------------------
# Table tab:sample_composition in Section 3.4.
wave_sample_sizes <- youth %>%
  group_by(wave) %>%
  summarise(
    n_person_wave_rows = n(),
    .groups = "drop"
  )

panel_counts <- youth %>%
  count(pidp, name = "n_waves") %>%
  count(n_waves, name = "n_adolescents") %>%
  mutate(share = n_adolescents / sum(n_adolescents))

wave_sample_sizes_export <- bind_rows(
  data.frame(
    statistic = paste("Wave", wave_sample_sizes$wave),
    count = wave_sample_sizes$n_person_wave_rows,
    stringsAsFactors = FALSE
  ),
  data.frame(
    statistic = "Pooled person-wave observations",
    count = nrow(youth),
    stringsAsFactors = FALSE
  )
)

panel_counts_export <- panel_counts %>%
  mutate(
    statistic = paste0("Observed in ", n_waves, ifelse(n_waves == 1, " wave", " waves")),
    count = n_adolescents
  ) %>%
  select(statistic, count) %>%
  bind_rows(
    data.frame(
      statistic = "Unique adolescents",
      count = n_distinct(youth$pidp),
      stringsAsFactors = FALSE
    )
  )

# Explained NL: maakt gewoon een frequentietabel van leeftijd. 
# count(age_dv, name = "n_observations") telt hoeveel rijen er per leeftijd zijn, 
# en arrange(age_dv) sorteert die leeftijden oplopend. 
# Belangrijk detail: dit zijn observations, dus person-wave rijen, niet per se unieke jongeren!!!
age_counts <- youth %>%
  count(age_dv, name = "n_observations") %>%
  arrange(age_dv)

show_table("Sample composition: wave sizes", wave_sample_sizes_export)
show_table("Sample composition: panel counts", panel_counts_export)
show_table("Sample composition: age counts", age_counts)

# Table tab:descriptive_statistics in Section 3.4.
youth$female <- ifelse(is.na(youth$sex_dv), NA_real_, as.integer(youth$sex_dv == 2))
youth$urban <- ifelse(is.na(youth$urban_dv), NA_real_, as.integer(youth$urban_dv == 1))

desc_spec <- data.frame(
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

# Explained NL: maakt voor één variabele een mini-samenvatting. 
# Eerst worden missings weggehaald, daarna berekent de functie N, mean, sd, min en max. 
# De ifelse(...)-stukjes zijn er om nette NAs terug te geven als er te weinig data is, 
# bijvoorbeeld geen standaarddeviatie bij maar 1 observatie. Later wordt die functie 
# op veel variabelen achter elkaar toegepast om de descriptieve tabel te bouwen.
summarise_one <- function(x) {
  x <- x[!is.na(x)]
  data.frame(
    N = length(x),
    mean = ifelse(length(x) == 0, NA_real_, mean(x)),
    sd = ifelse(length(x) <= 1, NA_real_, stats::sd(x)),
    min = ifelse(length(x) == 0, NA_real_, min(x)),
    max = ifelse(length(x) == 0, NA_real_, max(x)),
    stringsAsFactors = FALSE
  )
}

descriptive_statistics <- bind_cols(
  desc_spec,
  bind_rows(lapply(desc_spec$source, function(v) summarise_one(youth[[v]])))
)

show_table("Descriptive statistics", descriptive_statistics)

# 7. Figures for Section 3.4 -----------------------------------------------
# Shared labels for the treatment categories.
use_labels <- c(
  `1` = "None",
  `2` = "<1h",
  `3` = "1-3h",
  `4` = "4-6h",
  `5` = "7+h"
)

youth$use_cat <- factor(use_labels[as.character(youth$social_media_weekday)], levels = unname(use_labels))
youth$wave_display <- factor(youth$wave, levels = c("l", "m", "n", "o"), labels = c("Wave l", "Wave m", "Wave n", "Wave o"))

base_theme <- theme_minimal(base_size = 18) +
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

use_palette <- c(
  "None" = "#FFE3F1",
  "<1h" = "#F1C9FF",
  "1-3h" = "#D7B6FF",
  "4-6h" = "#A783FF",
  "7+h" = "#FF4FB3"
)

# Figure fig:descriptive_treatment_distribution.
distribution_df <- youth %>%
  filter(!is.na(use_cat), !is.na(wave_display)) %>%
  count(wave_display, use_cat) %>%
  group_by(wave_display) %>%
  mutate(share = n / sum(n))

fig_distribution <- ggplot(distribution_df, aes(x = wave_display, y = share, fill = use_cat)) +
  geom_col(width = 0.74, color = "white", linewidth = 0.35) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.04))) +
  scale_fill_manual(values = use_palette, name = "Weekday use") +
  labs(x = NULL, y = "Share within wave") +
  base_theme

show_plot("Figure: treatment distribution by wave", fig_distribution)

# Figure fig:descriptive_high_use_age_sex.
age_sex_df <- youth %>%
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


fig_high_use_age_sex <- ggplot(
  age_sex_df,
  aes(x = age_dv, y = share_high_use, color = sex_f, fill = sex_f, group = sex_f)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.24, linewidth = 0) +
  geom_line(linewidth = 1.35) +
  geom_point(size = 3.2) +
  scale_x_continuous(breaks = sort(unique(age_sex_df$age_dv))) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.03, 0.10))) +
  scale_color_manual(values = c("Male" = "#8C1AFF", "Female" = "#FF4FB3"), name = NULL) +
  scale_fill_manual(values = c("Male" = "#D7B6FF", "Female" = "#FFD0E8"), name = NULL) +
  labs(x = "Age", y = "Share using social media 4+ hours") +
  base_theme +
  theme(legend.position = "top")

show_plot("Figure: high social-media use by age and sex", fig_high_use_age_sex)

# Figure fig:descriptive_high_use_controls.
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
  youth %>%
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
  youth %>%
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

fig_high_use_controls <- ggplot(high_use_controls_df, aes(x = control_group, y = share_high_use, group = 1)) +
  geom_col(width = 0.72, fill = "#D7B6FF", color = "#8C1AFF", linewidth = 0.45) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.14, color = "#39224F", linewidth = 0.75) +
  facet_wrap(~ control, scales = "free_x", ncol = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.04))) +
  labs(x = NULL, y = "Share using social media 4+ hours") +
  base_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1, size = 13.5)
  )

show_plot("Figure: high social-media use by family meals and TV time", fig_high_use_controls)

# Figure fig:descriptive_high_use_household.
high_use_household_df <- bind_rows(
  youth %>%
    filter(!is.na(npns_dv), !is.na(high_social_media_weekday)) %>%
    mutate(
      control = "Parent figures in household",
      control_group = factor(as.character(npns_dv), levels = c("0", "1", "2"))
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
  youth %>%
    filter(!is.na(ngrp_dv), !is.na(high_social_media_weekday)) %>%
    mutate(
      control = "Grandparent in household",
      control_group = factor(ifelse(ngrp_dv > 0, "Yes", "No"), levels = c("No", "Yes"))
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
  youth %>%
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

fig_high_use_household <- ggplot(high_use_household_df, aes(x = control_group, y = share_high_use, group = 1)) +
  geom_col(width = 0.72, fill = "#D7B6FF", color = "#8C1AFF", linewidth = 0.45) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.14, color = "#39224F", linewidth = 0.75) +
  facet_wrap(~ control, scales = "free_x", ncol = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.04))) +
  labs(x = NULL, y = "Share using social media 4+ hours") +
  base_theme +
  theme(legend.position = "none")

show_plot("Figure: high social-media use by household composition", fig_high_use_household)

# Figure fig:descriptive_outcomes_by_treatment.
outcome_map <- c(
  loneliness = "Loneliness",
  life_dissatisfaction = "Life dissatisfaction",
  schoolwork_dissatisfaction = "School-work dissatisfaction",
  school_dissatisfaction = "School dissatisfaction"
)

outcome_df <- bind_rows(lapply(names(outcome_map), function(var) {
  data.frame(
    use_cat = youth$use_cat,
    outcome = outcome_map[[var]],
    value = youth[[var]],
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
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#D7B6FF", alpha = 0.28) +
  geom_line(color = "#8C1AFF", linewidth = 1.35) +
  geom_point(color = "#FF4FB3", size = 3.2) +
  facet_wrap(~ outcome, scales = "free_y", ncol = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.22))) +
  labs(x = NULL, y = "Mean outcome") +
  base_theme +
  theme(legend.position = "none")

show_plot("Figure: outcomes by treatment category", fig_outcomes)
