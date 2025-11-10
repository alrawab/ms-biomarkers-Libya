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
# ===============================================
# Table 2 (Enhanced): Vitamin D & Lipid Profile
# - Overall cohort summary with clinical thresholds
# - Adds Female/Male n (%) for each row (main + threshold)
# - Outputs: DOCX + CSV
# ===============================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
  library(yaml)
  library(flextable)
  library(officer)
})

# ---- Project helpers ----
source("analysis/00_setup/functions.R")
cfg <- yaml::read_yaml("analysis/00_setup/config.yml")
cn  <- cfg$columns

# ---- Load & prepare data ----
df  <- readr::read_csv(cfg$data_file) %>% janitor::clean_names()
names(df) <- gsub("\\.", "_", names(df))
df  <- derive_vars(df) %>% clean_levels()

# ---- Gender normalization ----
# Map any variants to "Female" / "Male"; everything else -> NA
normalize_sex <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    x_chr %in% c("f", "female", "امرأة", "انثى", "أنثى") ~ "Female",
    x_chr %in% c("m", "male", "رجل", "ذكر") ~ "Male",
    TRUE ~ NA_character_
  )
}
sex_vec <- if (cn$sex %in% names(df)) normalize_sex(df[[cn$sex]]) else rep(NA_character_, nrow(df))

# ---- Helpers ----
fmt_mean_sd_range <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(mean_sd = "NA", range = "NA"))
  c(
    mean_sd = sprintf("%.2f ± %.2f", mean(x), stats::sd(x)),
    range   = sprintf("%s–%s",
                      format(round(min(x), 2), nsmall = 2),
                      format(round(max(x), 2), nsmall = 2))
  )
}

count_pct <- function(flag) {
  n_all <- sum(!is.na(flag))
  n_yes <- sum(flag %in% TRUE, na.rm = TRUE)
  pct   <- if (n_all > 0) 100 * n_yes / n_all else NA_real_
  list(n = n_all, n_yes = n_yes, pct = pct)
}

sex_counts <- function(mask) {
  # mask = logical vector selecting the subset for this row
  idx <- which(mask %in% TRUE)
  s   <- sex_vec[idx]
  n   <- length(s)
  f_n <- sum(s == "Female", na.rm = TRUE)
  m_n <- sum(s == "Male", na.rm = TRUE)
  f_p <- if (n > 0) 100 * f_n / n else NA_real_
  m_p <- if (n > 0) 100 * m_n / n else NA_real_
  list(n = n, f_n = f_n, f_p = f_p, m_n = m_n, m_p = m_p)
}

# Build a two-row block:
# - Main row: stats over non-missing biomarker values; Female/Male within that non-missing set
# - Threshold row: Abnormal n (%); Female/Male within abnormal subset
row_block <- function(label_main, vec, abnormal_label = NULL, abn_flag = NULL) {
  # Main line (non-missing values)
  nonmiss_mask <- !is.na(vec)
  msr <- fmt_mean_sd_range(vec)
  sx  <- sex_counts(nonmiss_mask)
  
  top <- tibble::tibble(
    Variable             = label_main,
    n                    = sum(nonmiss_mask),
    `Mean ± SD (Range)`  = sprintf("%s (%s)", msr["mean_sd"], msr["range"]),
    `Abnormal, n (%)`    = "",
    `Female, n (%)`      = if (!is.na(sx$f_p)) sprintf("%d (%.1f%%)", sx$f_n, sx$f_p) else "",
    `Male, n (%)`        = if (!is.na(sx$m_p)) sprintf("%d (%.1f%%)", sx$m_n, sx$m_p) else ""
  )
  
  # Threshold line
  if (!is.null(abn_flag) && !is.null(abnormal_label)) {
    # Restrict to rows where threshold evaluable (non-missing flag)
    abn_evaluable <- !is.na(abn_flag)
    # Within evaluable, define subgroup = abnormal TRUE
    abn_mask <- abn_flag %in% TRUE
    abn_cts  <- count_pct(abn_flag)
    sx_abn   <- sex_counts(abn_mask & abn_evaluable)
    
    sub <- tibble::tibble(
      Variable             = paste0("  ", abnormal_label),
      n                    = abn_cts$n,                          # evaluable N
      `Mean ± SD (Range)`  = "",
      `Abnormal, n (%)`    = sprintf("%d (%.1f%%)", abn_cts$n_yes, abn_cts$pct),
      `Female, n (%)`      = if (!is.na(sx_abn$f_p)) sprintf("%d (%.1f%%)", sx_abn$f_n, sx_abn$f_p) else "",
      `Male, n (%)`        = if (!is.na(sx_abn$m_p)) sprintf("%d (%.1f%%)", sx_abn$m_n, sx_abn$m_p) else ""
    )
    dplyr::bind_rows(top, sub)
  } else {
    top
  }
}

# ---- Build table ----
tab_list <- list()

# Convenient getter; returns numeric vector or NA
getv <- function(name) if (name %in% names(df)) suppressWarnings(as.numeric(df[[name]])) else rep(NA_real_, nrow(df))

vitd <- getv(cn$vitd)
tc   <- getv(cn$tc)
tg   <- getv(cn$tg)
hdl  <- getv(cn$hdl)
ldl  <- getv(cn$ldl)

# Vitamin D (ng/mL): insufficient < 30
# ---- Vitamin D categories ----
if (any(!is.na(vitd))) {
  # Main line (mean ± SD, range)
  tab_vitd <- row_block("Vitamin D, ng/mL", vitd)
  
  # Category masks
  def_flag  <- vitd < 10
  ins_flag  <- vitd >= 10 & vitd <= 30
  suf_flag  <- vitd > 30
  
  # Helper for category lines
  vitd_cat_line <- function(label, flag) {
    cts <- count_pct(flag)
    sx  <- sex_counts(flag)
    tibble::tibble(
      Variable = paste0("  ", label),
      n        = cts$n,
      `Mean ± SD (Range)` = "",
      `Abnormal, n (%)`   = sprintf("%d (%.1f%%)", cts$n_yes, cts$pct),
      `Female, n (%)`     = if (!is.na(sx$f_p)) sprintf("%d (%.1f%%)", sx$f_n, sx$f_p) else "",
      `Male, n (%)`       = if (!is.na(sx$m_p)) sprintf("%d (%.1f%%)", sx$m_n, sx$m_p) else ""
    )
  }
  
  # Bind all rows together
  tab_vitd <- bind_rows(
    tab_vitd,
    vitd_cat_line("Deficient (<10 ng/mL)", def_flag),
    vitd_cat_line("Insufficient (10–30 ng/mL)", ins_flag),
    vitd_cat_line("Sufficient (>30 ng/mL)", suf_flag)
  )
  
  tab_list[["vitd"]] <- tab_vitd
}


# Total Cholesterol (mg/dL): high ≥ 240
if (any(!is.na(tc))) {
  tab_list[["tc"]] <- row_block(
    "Total Cholesterol, mg/dL",
    tc,
    abnormal_label = "High TC (≥240 mg/dL)",
    abn_flag       = tc >= 240
  )
}

# TG (mg/dL): high ≥ 200
if (any(!is.na(tg))) {
  tab_list[["tg"]] <- row_block(
    "TG, mg/dL",
    tg,
    abnormal_label = "High TG (≥200 mg/dL)",
    abn_flag       = tg >= 200
  )
}

# HDL (mg/dL): low < 40
if (any(!is.na(hdl))) {
  tab_list[["hdl"]] <- row_block(
    "HDL, mg/dL",
    hdl,
    abnormal_label = "Low HDL (<40 mg/dL)",
    abn_flag       = hdl < 40
  )
}

# LDL (mg/dL): high ≥ 160
if (any(!is.na(ldl))) {
  tab_list[["ldl"]] <- row_block(
    "LDL, mg/dL",
    ldl,
    abnormal_label = "High LDL (≥160 mg/dL)",
    abn_flag       = ldl >= 160
  )
}

# Ratios / AIP (main lines only)
if ("chol_hdl_ratio" %in% names(df)) {
  chr <- suppressWarnings(as.numeric(df[["chol_hdl_ratio"]]))
  tab_list[["chr"]] <- row_block("Chol/HDL ratio", chr)
}
if ("ldl_hdl_ratio" %in% names(df)) {
  lhr <- suppressWarnings(as.numeric(df[["ldl_hdl_ratio"]]))
  tab_list[["lhr"]] <- row_block("LDL/HDL ratio", lhr)
}
if ("aip" %in% names(df)) {
  aip <- suppressWarnings(as.numeric(df[["aip"]]))
  tab_list[["aip"]] <- row_block("Atherogenic Index of Plasma (AIP)", aip)
}

tbl <- bind_rows(tab_list)

# ---- Export CSV ----
readr::write_csv(tbl, "outputs/tables/Table2_Biochemical_Enhanced.csv")

# ---- Flextable formatting ----
ft <- flextable(tbl)
ft <- autofit(ft)
ft <- align(ft, j = c("n", "Abnormal, n (%)", "Female, n (%)", "Male, n (%)"),
            align = "center", part = "body")

main_rows <- which(!grepl("^\\s", tbl$Variable))
sub_rows  <- which(grepl("^\\s", tbl$Variable))

if (length(main_rows)) ft <- bold(ft, i = main_rows, bold = TRUE)
if (length(sub_rows)) {
  ft <- padding(ft, i = sub_rows, j = "Variable", padding.left = 18)
  ft <- bg(ft, i = sub_rows, bg = "#FBE5E5")
}

ft <- set_header_labels(
  ft,
  Variable = "Variable",
  n = "n",
  `Mean ± SD (Range)` = "Mean ± SD (Range)",
  `Abnormal, n (%)` = "Abnormal, n (%)",
  `Female, n (%)` = "Female, n (%)",
  `Male, n (%)`   = "Male, n (%)"
)
ft <- theme_booktabs(ft)

# ---- Create Word doc ----
N <- sum(
  !is.na(vitd) | !is.na(tc) | !is.na(tg) | !is.na(hdl) | !is.na(ldl)
)

doc <- read_docx()
doc <- body_add_par(
  doc,
  sprintf("Table 2. Vitamin D and lipid profile levels (overall; n=%d)", N),
  style = "heading 2"
)
doc <- body_add_par(
  doc,
  "HDL: high-density lipoprotein; TC: total cholesterol; TG: triglycerides; SD: standard deviation.",
  style = "Normal"
)
doc <- flextable::body_add_flextable(doc, value = ft)

print(doc, target = "outputs/tables/Table2_Biochemical_Enhanced.docx")

message("✅ Table2_Biochemical_Enhanced.docx and CSV created successfully.")
