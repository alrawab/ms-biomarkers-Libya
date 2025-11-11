#**********************************************************************************************#
#                                            MS Biomarkers
#  Copyright (C) 2024-2025 Osamah Alrouwab <rawab@uoz.edu.ly> Faculty of Medicine Zintan-university-Libya
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
#************************************************************************************************#
# ============================================================
# 10_descriptives.R
# Baseline descriptives for MS cohort
# - Table 1: Clinical & Demographic by EDSS (includes Age at Onset)
# - Table 2: Biochemical by EDSS
# - Exports .docx via flextable/officer (no Pandoc needed) + CSVs
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
  library(yaml)
  library(gtsummary)
  library(flextable)
  library(officer)
})

# ---------------- Project setup ----------------
# Requires:
#  - analysis/00_setup/config.yml  (with data_file and column mappings)
#  - analysis/00_setup/functions.R (derive_vars, clean_levels)
source("analysis/00_setup/functions.R")
cfg <- yaml::read_yaml("analysis/00_setup/config.yml")
cn  <- cfg$columns

# ---------------- Load & prep data ----------------
df <- readr::read_csv(cfg$data_file) %>% janitor::clean_names()
names(df) <- gsub("\\.", "_", names(df))          # align with helpers
df <- derive_vars(df) %>% clean_levels()

# Ensure/derive EDSS categories
if (!"EDSS_cat" %in% names(df)) {
  edss_v <- suppressWarnings(as.numeric(df[[cn$edss]]))
  df$EDSS_cat <- dplyr::case_when(
    is.na(edss_v) ~ NA_character_,
    edss_v <= cfg$edss_categories$mild_max ~ "Mild (0–2.5)",
    edss_v <= cfg$edss_categories$moderate_max ~ "Moderate (3–5.5)",
    is.finite(edss_v) ~ "Severe (≥6.0)"
  )
}

# ---------------- Helpers ----------------
add_headers_safely <- function(tbl, overall_label = "**Overall**", edss_label = "**EDSS category**") {
  # add_overall() creates stat_0 for "Overall" in gtsummary (guard across versions)
  if (!is.null(tbl$table_styling) && "header" %in% names(tbl$table_styling)) {
    hdr_cols <- tryCatch(colnames(tbl$table_styling$header), error = function(e) character())
    if ("stat_0" %in% hdr_cols) {
      tbl <- tbl %>% modify_header(stat_0 ~ overall_label)
    }
  }
  tbl %>% modify_header(all_stat_cols() ~ edss_label)
}

# Ensure output dir
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# ============================================================
#                       TABLE 1
# Clinical & Demographic by EDSS (with Age at Onset)
# ============================================================

keep1 <- c(
  cn$age,               # Age at assessment
  cn$age_onset,         # Age at onset (your header -> age_of_onset_of_disease_years)
  cn$sex,               # Sex
  cn$phenotype,         # MS phenotype
  cn$duration_years,    # Duration of illness (years)
  "arr",                # Annualized relapse rate (present or derived)
  cn$bmi                # BMI (kg/m^2)
)
keep1 <- keep1[keep1 %in% names(df)]

t1 <- df %>%
  select(EDSS_cat, all_of(keep1)) %>%
  gtsummary::tbl_summary(
    by = EDSS_cat,
    missing = "no",
    type = list(where(is.numeric) ~ "continuous2"),
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}\n{median} ({p25}–{p75})",
      all_categorical() ~ "{n} ({p}%)"
    )
  ) %>%
  add_overall(last = TRUE) %>%
  add_headers_safely(overall_label = "**Overall**", edss_label = "**EDSS category**") %>%
  modify_header(label ~ "**Characteristic**") %>%
  bold_labels()

# DOCX export
doc1 <- read_docx()
doc1 <- body_add_par(doc1, "Table 1. Clinical and demographic characteristics by EDSS category", style = "heading 2")
doc1 <- body_add_flextable(doc1, as_flex_table(t1))
print(doc1, target = "outputs/tables/Table1_Clinical_by_EDSS.docx")

# QA CSV: quick means by group
t1_qc <- df %>%
  group_by(EDSS_cat) %>%
  summarize(across(all_of(keep1), \(x) suppressWarnings(mean(as.numeric(x), na.rm = TRUE))), .groups = "drop")
readr::write_csv(t1_qc, "outputs/tables/Table1_Clinical_by_EDSS.csv")

message("✓ Wrote Table1_Clinical_by_EDSS.docx and CSV")

# ============================================================
#                       TABLE 2
# Biochemical by EDSS (Mean ± SD + Median (IQR))
# ============================================================

bio <- c(
  cn$vitd, cn$crp, cn$hba1c, cn$tc, cn$ldl, cn$hdl, cn$tg,
  "chol_hdl_ratio","ldl_hdl_ratio","aip"
)
bio <- bio[bio %in% names(df)]

t2 <- df %>%
  select(EDSS_cat, all_of(bio)) %>%
  gtsummary::tbl_summary(
    by = EDSS_cat,
    missing = "no",
    type = list(where(is.numeric) ~ "continuous2"),
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}\n{median} ({p25}–{p75})"
    )
  ) %>%
  add_overall(last = TRUE) %>%
  add_headers_safely(overall_label = "**Overall**", edss_label = "**EDSS category**") %>%
  modify_header(label ~ "**Biomarker**") %>%
  bold_labels()

doc2 <- read_docx()
doc2 <- body_add_par(doc2, "Table 2. Biochemical markers by EDSS category", style = "heading 2")
doc2 <- body_add_flextable(doc2, as_flex_table(t2))
print(doc2, target = "outputs/tables/Table2_Biochemical_by_EDSS.docx")

# Raw dump (for checks or downstream enhanced tables)
readr::write_csv(
  df %>% select(EDSS_cat, all_of(bio)),
  "outputs/tables/Table2_Biochemical_by_EDSS_raw.csv"
)

message(" Wrote Table2_Biochemical_by_EDSS.docx and raw CSV")
