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
# Table 2 (EDSS-stratified): Biochemical markers by disability
# - Columns: Mild / Moderate / Severe / Overall
# - Rows: Sex distribution; Vit D (with categories); TC/TG/HDL/LDL (+ thresholds);
#         Chol/HDL, LDL/HDL, AIP
# - Adds Female/Male n(%) for every category/threshold row
# - Outputs: DOCX + CSV
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
  library(yaml)
  library(flextable)
  library(officer)
})

# --- Project setup (uses your repo files) ---
source("analysis/00_setup/functions.R")
cfg <- yaml::read_yaml("analysis/00_setup/config.yml")
cn  <- cfg$columns

# --- Load & prep data ---
df <- readr::read_csv(cfg$data_file) %>% clean_names()
names(df) <- gsub("\\.", "_", names(df))
df <- derive_vars(df) %>% clean_levels()

# --- Ensure EDSS numeric & EDSS group ---
stopifnot(cn$edss %in% names(df))
edss_num <- suppressWarnings(as.numeric(df[[cn$edss]]))
df <- df %>%
  mutate(
    edss_num = edss_num,
    EDSS_group = case_when(
      !is.na(edss_num) & edss_num <= 2.5 ~ "Mild (0–2.5)",
      !is.na(edss_num) & edss_num <= 5.5 ~ "Moderate (3–5.5)",
      !is.na(edss_num) & edss_num > 5.5  ~ "Severe (≥6)",
      TRUE ~ NA_character_
    )
  )

# --- Normalize sex for counts ---
normalize_sex <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  case_when(
    x_chr %in% c("f","female","امرأة","انثى","أنثى") ~ "Female",
    x_chr %in% c("m","male","رجل","ذكر") ~ "Male",
    TRUE ~ NA_character_
  )
}
sex_vec <- normalize_sex(df[[cn$sex]])

# --- EDSS-group sizes (n) ---
n_mild    <- sum(df$EDSS_group == "Mild (0–2.5)",    na.rm = TRUE)
n_moder   <- sum(df$EDSS_group == "Moderate (3–5.5)",na.rm = TRUE)
n_severe  <- sum(df$EDSS_group == "Severe (≥6)",     na.rm = TRUE)
n_overall <- nrow(df)

# --- Helper: format continuous cell ---
cell_cont <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return("NA")
  m  <- mean(x); s <- stats::sd(x)
  med <- stats::median(x); q1 <- as.numeric(stats::quantile(x, 0.25)); q3 <- as.numeric(stats::quantile(x, 0.75))
  sprintf("%.2f ± %.2f [%.2f (%.2f–%.2f)]", m, s, med, q1, q3)
}

# --- Helper: build a row of continuous stats across 4 columns (by groups + overall) ---
row_cont <- function(label, vec) {
  tibble(
    Characteristic = label,
    `Mild (0–2.5)\nN = {n_mild}`    := cell_cont(vec[df$EDSS_group == "Mild (0–2.5)"]),
    `Moderate (3–5.5)\nN = {n_moder}`:= cell_cont(vec[df$EDSS_group == "Moderate (3–5.5)"]),
    `Severe (≥6)\nN = {n_severe}`    := cell_cont(vec[df$EDSS_group == "Severe (≥6)"]),
    `Overall\nN = {n_overall}`       := cell_cont(vec)
  )
}

# --- Helper: frequency cell with sex split for a logical flag (TRUE = category present) ---
cell_freq_sex <- function(flag, universe_mask) {
  # flag: logical vector (TRUE = counted). universe_mask: logical rows belonging to stratum.
  idx <- which(universe_mask & !is.na(flag))
  if (!length(idx)) return(c("n(%)"="0 (0.0%)","F"="","M"=""))
  n_eval <- length(idx)
  n_yes  <- sum(flag[idx], na.rm = TRUE)
  pct    <- if (n_eval > 0) 100 * n_yes / n_eval else NA_real_
  
  s <- sex_vec[idx]
  f_n <- sum(s == "Female" & flag[idx], na.rm = TRUE)
  m_n <- sum(s == "Male"   & flag[idx], na.rm = TRUE)
  # percents within the abnormal subset:
  f_p <- if (n_yes > 0) 100 * f_n / n_yes else NA_real_
  m_p <- if (n_yes > 0) 100 * m_n / n_yes else NA_real_
  
  out1 <- sprintf("%d (%.1f%%)", n_yes, pct)
  outF <- if (!is.na(f_p)) sprintf("%d (%.1f%%)", f_n, f_p) else ""
  outM <- if (!is.na(m_p)) sprintf("%d (%.1f%%)", m_n, m_p) else ""
  c(out1, outF, outM)
}

# --- Helper: build a 3-column block for a categorical/threshold line across 4 strata ---
row_cat <- function(label, flag) {
  mild_m  <- df$EDSS_group == "Mild (0–2.5)"
  moder_m <- df$EDSS_group == "Moderate (3–5.5)"
  severe_m<- df$EDSS_group == "Severe (≥6)"
  all_m   <- rep(TRUE, nrow(df))
  
  c_mild   <- cell_freq_sex(flag, mild_m)
  c_moder  <- cell_freq_sex(flag, moder_m)
  c_severe <- cell_freq_sex(flag, severe_m)
  c_over   <- cell_freq_sex(flag, all_m)
  
  tibble(
    Characteristic = paste0("  ", label),
    `Mild (0–2.5)\nN = {n_mild}`     = sprintf("%s | F: %s | M: %s", c_mild[1],  c_mild[2],  c_mild[3]),
    `Moderate (3–5.5)\nN = {n_moder}`= sprintf("%s | F: %s | M: %s", c_moder[1], c_moder[2], c_moder[3]),
    `Severe (≥6)\nN = {n_severe}`     = sprintf("%s | F: %s | M: %s", c_severe[1],c_severe[2],c_severe[3]),
    `Overall\nN = {n_overall}`        = sprintf("%s | F: %s | M: %s", c_over[1],  c_over[2],  c_over[3])
  )
}

# --- Convenience getters (numeric) ---
getv <- function(name) if (name %in% names(df)) suppressWarnings(as.numeric(df[[name]])) else rep(NA_real_, nrow(df))
vitd <- getv(cn$vitd)
tc   <- getv(cn$tc)
tg   <- getv(cn$tg)
hdl  <- getv(cn$hdl)
ldl  <- getv(cn$ldl)
chr  <- if ("chol_hdl_ratio" %in% names(df)) suppressWarnings(as.numeric(df[["chol_hdl_ratio"]])) else rep(NA_real_, nrow(df))
lhr  <- if ("ldl_hdl_ratio"  %in% names(df)) suppressWarnings(as.numeric(df[["ldl_hdl_ratio"]]))  else rep(NA_real_, nrow(df))
aip  <- if ("aip"            %in% names(df)) suppressWarnings(as.numeric(df[["aip"]]))            else rep(NA_real_, nrow(df))

# --- Sex distribution row (counts & %) ---
sex_row_cell <- function(mask) {
  idx <- which(mask)
  if (!length(idx)) return("—")
  n <- length(idx)
  f_n <- sum(sex_vec[idx] == "Female", na.rm = TRUE)
  m_n <- sum(sex_vec[idx] == "Male",   na.rm = TRUE)
  sprintf("%d F (%.1f%%), %d M (%.1f%%)", f_n, 100*f_n/n, m_n, 100*m_n/n)
}
row_sex <- tibble(
  Characteristic = "Sex distribution",
  `Mild (0–2.5)\nN = {n_mild}`      = sex_row_cell(df$EDSS_group == "Mild (0–2.5)"),
  `Moderate (3–5.5)\nN = {n_moder}` = sex_row_cell(df$EDSS_group == "Moderate (3–5.5)"),
  `Severe (≥6)\nN = {n_severe}`      = sex_row_cell(df$EDSS_group == "Severe (≥6)"),
  `Overall\nN = {n_overall}`         = sex_row_cell(rep(TRUE, nrow(df)))
)

# --- Build table blocks ---

tbl <- list()

# Sex distribution
tbl[["sex"]] <- row_sex

# Vitamin D continuous + categories
if (any(is.finite(vitd))) {
  tbl[["vitd_main"]] <- row_cont("Vitamin D (ng/mL)", vitd)
  
  tbl[["vitd_def"]] <- row_cat("Deficient (<10 ng/mL)",       vitd < 10)
  tbl[["vitd_ins"]] <- row_cat("Insufficient (10–30 ng/mL)",  vitd >= 10 & vitd <= 30)
  tbl[["vitd_suf"]] <- row_cat("Sufficient (>30 ng/mL)",      vitd > 30)
}

# Lipids continuous + thresholds
if (any(is.finite(tc))) {
  tbl[["tc_main"]]  <- row_cont("Total Cholesterol (mg/dL)", tc)
  tbl[["tc_high"]]  <- row_cat("High TC (≥240 mg/dL)", tc >= 240)
}
if (any(is.finite(tg))) {
  tbl[["tg_main"]]  <- row_cont("Triglycerides (mg/dL)", tg)
  tbl[["tg_high"]]  <- row_cat("High TG (≥200 mg/dL)", tg >= 200)
}
if (any(is.finite(hdl))) {
  tbl[["hdl_main"]] <- row_cont("HDL (mg/dL)", hdl)
  tbl[["hdl_low"]]  <- row_cat("Low HDL (<40 mg/dL)", hdl < 40)
}
if (any(is.finite(ldl))) {
  tbl[["ldl_main"]] <- row_cont("LDL (mg/dL)", ldl)
  tbl[["ldl_high"]] <- row_cat("High LDL (≥160 mg/dL)", ldl >= 160)
}
if (any(is.finite(chr))) tbl[["chr_main"]] <- row_cont("Chol/HDL ratio", chr)
if (any(is.finite(lhr))) tbl[["lhr_main"]] <- row_cont("LDL/HDL ratio",  lhr)
if (any(is.finite(aip))) tbl[["aip_main"]] <- row_cont("Atherogenic Index of Plasma (AIP)", aip)

# Bind all blocks
out_tbl <- bind_rows(tbl)

# Make pretty column headers with n
colnames(out_tbl) <- c(
  "Characteristic",
  sprintf("Mild (0–2.5)\nN = %d", n_mild),
  sprintf("Moderate (3–5.5)\nN = %d", n_moder),
  sprintf("Severe (≥6)\nN = %d", n_severe),
  sprintf("Overall\nN = %d", n_overall)
)

# Write CSV
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(out_tbl, "outputs/tables/Table2_Biochemical_by_EDSS.csv")

# --- Word export (flextable) ---
ft <- flextable(out_tbl)
ft <- autofit(ft)
# style: bold main lines (no leading spaces); indent sub-lines
main_rows <- which(!grepl("^\\s", out_tbl$Characteristic))
sub_rows  <- which(grepl("^\\s", out_tbl$Characteristic))
if (length(main_rows)) ft <- bold(ft, i = main_rows, bold = TRUE)
if (length(sub_rows)) {
  ft <- padding(ft, i = sub_rows, j = 1, padding.left = 18)
  ft <- bg(ft, i = sub_rows, bg = "#F5F5F5")
}
ft <- align(ft, j = 2:5, align = "center")
ft <- theme_booktabs(ft)

doc <- read_docx()
doc <- body_add_par(doc, "Table 2. Biochemical markers stratified by EDSS severity", style = "heading 2")
doc <- body_add_par(doc, "Cell format for continuous variables: Mean ± SD [Median (Q1–Q3)]. Category rows show n (%) with female/male split ‘n (%)’.", style = "Normal")
doc <- flextable::body_add_flextable(doc, value = ft)

print(doc, target = "outputs/tables/Table2_Biochemical_by_EDSS.docx")

message("✅ Wrote: outputs/tables/Table2_Biochemical_by_EDSS.docx and outputs/tables/Table2_Biochemical_by_EDSS.csv")
