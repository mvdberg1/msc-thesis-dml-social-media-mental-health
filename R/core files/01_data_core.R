# Core file 1: data selection and cleaning.
#
# Purpose: read the licensed UKHLS youth files, keep the variables used in the
# thesis analysis, construct the treatment/outcome/control variables, and save a
# cleaned person-wave dataset for the methodology and results scripts.

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
})

project_dir <- normalizePath(
  Sys.getenv("THESIS_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
raw_dir <- Sys.getenv(
  "UKHLS_RAW_DIR",
  unset = file.path(project_dir, "UKDA-6614-stata", "stata", "stata14_se", "ukhls")
)
analysis_dir <- file.path(project_dir, "data", "analysis")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

waves <- c("l", "m", "n", "o")

treatment_vars <- c("ypnetcht", "ypnetchtw", "ypsocweb")
outcome_vars <- c("yplonely", "yphlf", "yphsw", "yphsc")
control_vars <- c(
  "age_dv", "sex_dv", "ethn_dv", "country", "urban_dv", "gor_dv",
  "npns_dv", "ngrp_dv", "nnssib_dv", "ypnpal", "ypeatlivu",
  "yptvvidhrs", paste0("ypdevice", 1:6)
)
design_vars <- c("psu", "strata", "ythscui_xw")
raw_stems <- c(design_vars, treatment_vars, outcome_vars, control_vars)

clean_numeric <- function(x) {
  x <- as.numeric(zap_labels(x))
  x[x < 0] <- NA_real_
  x
}

factor_from_codes <- function(x, labels) {
  factor(labels[as.character(x)], levels = unname(labels))
}

read_youth_wave <- function(wave) {
  path <- file.path(raw_dir, paste0(wave, "_youth.dta"))
  keep <- c("pidp", paste0(wave, "_", raw_stems))
  df <- as.data.frame(read_dta(path, col_select = all_of(keep)))
  names(df) <- sub(paste0("^", wave, "_"), "", names(df))
  df$wave <- wave
  df$wave_number <- match(wave, letters)
  df
}

youth <- bind_rows(lapply(waves, read_youth_wave))
numeric_cols <- setdiff(names(youth), "wave")
youth[numeric_cols] <- lapply(youth[numeric_cols], clean_numeric)

# Treatments.
youth$social_media_weekday <- youth$ypnetcht
youth$social_media_weekend <- youth$ypnetchtw
youth$high_social_media_weekday <- ifelse(
  is.na(youth$social_media_weekday),
  NA_integer_,
  as.integer(youth$social_media_weekday >= 4)
)
youth$high_social_media_weekend <- ifelse(
  is.na(youth$social_media_weekend),
  NA_integer_,
  as.integer(youth$social_media_weekend >= 4)
)
youth$has_social_profile <- ifelse(
  is.na(youth$ypsocweb),
  NA_integer_,
  as.integer(youth$ypsocweb == 1)
)

# Outcomes. In these UKHLS items, higher retained values indicate worse
# wellbeing, so no reverse-coding is applied.
youth$loneliness <- youth$yplonely
youth$life_dissatisfaction <- youth$yphlf
youth$schoolwork_dissatisfaction <- youth$yphsw
youth$school_dissatisfaction <- youth$yphsc

# Controls.
youth$close_friends <- ifelse(youth$ypnpal > 100, NA_real_, youth$ypnpal)
youth$close_friends_log <- log1p(youth$close_friends)
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
names(youth)[names(youth) == "ythscui_xw"] <- "youth_weight"

front_cols <- c(
  "pidp", "wave", "wave_number", "wave_f", "youth_weight", "psu", "strata",
  "social_media_weekday", "social_media_weekend",
  "high_social_media_weekday", "high_social_media_weekend",
  "has_social_profile", "loneliness", "life_dissatisfaction",
  "schoolwork_dissatisfaction", "school_dissatisfaction",
  "age_dv", "sex_f", "ethnicity_broad_f", "country_f", "urban_f", "region_f",
  "npns_dv", "ngrp_dv", "nnssib_dv", "close_friends", "close_friends_log",
  "ypeatlivu", "yptvvidhrs", paste0("ypdevice", 1:6)
)
youth <- youth[c(front_cols, setdiff(names(youth), front_cols))]

clean_path <- file.path(analysis_dir, "ukhls_youth_l_to_o_clean.rds")
saveRDS(youth, clean_path)

sample_summary <- youth %>%
  summarise(
    person_wave_rows = n(),
    adolescents = n_distinct(pidp),
    waves = n_distinct(wave),
    weekday_treatment_nonmissing = sum(!is.na(social_media_weekday)),
    high_weekday_nonmissing = sum(!is.na(high_social_media_weekday)),
    life_outcome_nonmissing = sum(!is.na(life_dissatisfaction))
  )
write.csv(
  sample_summary,
  file.path(analysis_dir, "core_data_summary.csv"),
  row.names = FALSE
)

message("Saved cleaned analysis data to: ", clean_path)
print(sample_summary)
