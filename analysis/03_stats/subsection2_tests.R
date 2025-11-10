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

# =====================================================================
# Subsection 2 — EDSS Associations (v2 robust)
# Table 3: Kruskal–Wallis + Jonckheere–Terpstra (trend) by EDSS category
# Table 4: Spearman correlations (crude) + optional partial (age±sex)
# Exports: CSV + DOCX
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
  library(rstatix)    # Kruskal–Wallis helper
  library(flextable)  # DOCX export
})

# ---------- Optional trend test (clinfun) ----------
jt_available <- requireNamespace("clinfun", quietly = TRUE)
jt_p <- function(y, g) {
  if (!jt_available) return(NA_real_)
  g <- factor(g, levels = c("Mild (0–2.5)","Moderate (3–5.5)","Severe (≥6)"), ordered = TRUE)
  out <- suppressWarnings(clinfun::jonckheere.test(y, g, alternative = "two.sided"))
  unname(out$p.value)
}

# ---------- Paths ----------
data_file <- "data/Libya_MS.csv"               # adjust if needed
out_dir   <- "outputs/tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- Load & clean ----------
df <- read_csv(data_file, show_col_types = FALSE) %>% clean_names()
names(df) <- gsub("\\.", "_", names(df))

# ---------- Helpers ----------
col_or <- function(...) {
  cands <- c(...)
  cands <- cands[!is.na(cands)]
  cands[cands %in% names(df)][1]
}
to_num <- function(x) suppressWarnings(as.numeric(x))

# ---------- EDSS numeric + category ----------
edss_col <- col_or("edss")
if (is.na(edss_col)) stop("No 'edss' column found in data.")
df <- df %>%
  mutate(
    edss_num = to_num(.data[[edss_col]]),
    edss_cat = case_when(
      is.na(edss_num) ~ NA_character_,
      edss_num <= 2.5 ~ "Mild (0–2.5)",
      edss_num <= 5.5 ~ "Moderate (3–5.5)",
      edss_num >  5.5 ~ "Severe (≥6)"
    ),
    edss_cat = factor(edss_cat, levels = c("Mild (0–2.5)","Moderate (3–5.5)","Severe (≥6)"))
  )

# ---------- Candidate variables (robust to your names) ----------
vars <- list(
  age            = col_or("age"),
  age_onset      = col_or("age_of_onset_of_disease_years","age_at_onset","age_at_onset_years"),
  duration_years = col_or("duration_of_illness_years","duration_years"),
  arr            = col_or("arr","annualized_relapse_rate"),
  bmi            = col_or("bmi_kg_m2","bmi"),
  vit_d          = col_or("vitamin_d","vit_d","vitamin_d_ng_ml","vitd"),
  tc             = col_or("total_cholesterol","tc"),
  ldl            = col_or("ldl"),
  hdl            = col_or("hdl"),
  tg             = col_or("tg","triglycerides"),
  chr            = col_or("chol_hdl_ratio","total_cholesterol_to_hdl_ratio"),
  lhr            = col_or("ldl_hdl_ratio"),
  aip            = col_or("aip")
)
vars <- vars[!is.na(unlist(vars))]

# Numeric copy
df_num <- df %>%
  mutate(across(all_of(unlist(vars)), to_num))

# ---------- Pretty labels via safe join (no recode masking) ----------
var_labels <- c(
  age="Age (years)", age_onset="Age at onset (years)", duration_years="Disease duration (years)",
  arr="Annualized relapse rate (ARR)", bmi="BMI (kg/m²)",
  vit_d="25(OH) Vitamin D (ng/mL)", tc="Total cholesterol (mg/dL)",
  ldl="LDL (mg/dL)", hdl="HDL (mg/dL)", tg="Triglycerides (mg/dL)",
  chr="Total cholesterol/HDL", lhr="LDL/HDL", aip="Atherogenic Index of Plasma (AIP)"
)
label_df <- tibble(variable = names(var_labels),
                   Variable = unname(var_labels))

# ===================== TABLE 3: EDSS group tests ======================
kw_results <- map_dfr(names(vars), function(vname) {
  vcol <- vars[[vname]]
  x <- df_num[[vcol]]
  g <- df_num$edss_cat
  ok <- is.finite(x) & !is.na(g)
  if (!any(ok)) {
    return(tibble(variable = vname, n = 0,
                  kruskal_p = NA_real_, jt_trend_p = NA_real_, direction = NA_character_))
  }
  
  # Kruskal–Wallis
  kw <- tryCatch(rstatix::kruskal_test(data.frame(x = x[ok], g = g[ok]), x ~ g), error = function(e) NULL)
  kw_p <- if (!is.null(kw)) kw$p else NA_real_
  
  # Trend p (JT)
  jt_pv <- tryCatch(jt_p(x[ok], g[ok]), error = function(e) NA_real_)
  
  # Direction (Spearman vs ordered EDSS levels)
  rho <- suppressWarnings(cor(x[ok], as.numeric(g[ok]), method = "spearman"))
  dir <- case_when(
    is.na(rho)    ~ NA_character_,
    rho >  0.05   ~ "↑ with EDSS",
    rho < -0.05   ~ "↓ with EDSS",
    TRUE          ~ "flat"
  )
  
  tibble(variable = vname,
         n        = sum(ok),
         kruskal_p  = kw_p,
         jt_trend_p = jt_pv,
         direction  = dir)
})

kw_results <- kw_results %>%
  mutate(variable = as.character(variable)) %>%
  left_join(label_df, by = "variable") %>%
  mutate(Variable = coalesce(Variable, variable),
         `Kruskal–Wallis p` = ifelse(is.na(kruskal_p), NA, format.pval(kruskal_p, digits = 3, eps = 1e-4)),
         `Trend p (JT)`     = ifelse(is.na(jt_trend_p), NA, format.pval(jt_trend_p, digits = 3, eps = 1e-4))) %>%
  select(Variable, n, `Kruskal–Wallis p`, `Trend p (JT)`, direction)

# Save Table 3
write_csv(kw_results, file.path(out_dir, "Table3_EDSS_group_tests.csv"))

ft3 <- flextable(kw_results) |>
  autofit() |>
  align(align = "left", j = 1) |>
  bold(j = 1) |>
  set_caption("Table 3. Comparisons across EDSS categories: Kruskal–Wallis and Jonckheere–Terpstra trend tests.")
save_as_docx(`Table 3` = ft3, path = file.path(out_dir, "Table3_EDSS_group_tests.docx"))

# ===================== TABLE 4: Spearman correlations =================
sex_col <- col_or("gender","sex")
has_sex <- !is.na(sex_col)

spearman_results <- map_dfr(names(vars), function(vname) {
  vcol <- vars[[vname]]
  x <- df_num[[vcol]]         # variable
  y <- df_num$edss_num        # EDSS numeric
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) {
    return(tibble(variable=vname, n=0, rho=NA_real_, p_value=NA_real_,
                  rho_partial=NA_real_, p_partial=NA_real_))
  }
  
  # crude Spearman
  cr <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  rho <- unname(cr$estimate); pv <- cr$p.value; nn <- sum(ok)
  
  # partial (age±sex) via residual-on-residual if age present
  age_col <- vars[["age"]] %||% NA_character_
  if (!is.na(age_col)) {
    df_adj <- tibble(
      x = x, y = y,
      age = to_num(df_num[[age_col]]),
      sex = if (has_sex) factor(df_num[[sex_col]]) else NULL
    )
    
    covars <- df_adj %>% select(age, any_of("sex"))
    get_resid <- function(resp, covars) {
      mod <- tryCatch(lm(resp ~ ., data = covars), error = function(e) NULL)
      if (is.null(mod)) return(rep(NA_real_, nrow(covars)))
      resid(mod)
    }
    rx <- get_resid(df_adj$x, covars)
    ry <- get_resid(df_adj$y, covars)
    okp <- is.finite(rx) & is.finite(ry)
    pr <- suppressWarnings(cor.test(rx[okp], ry[okp], method = "spearman", exact = FALSE))
    rho_p <- unname(pr$estimate); pv_p <- pr$p.value
  } else {
    rho_p <- NA_real_; pv_p <- NA_real_
  }
  
  tibble(variable=vname, n=nn, rho=rho, p_value=pv,
         rho_partial=rho_p, p_partial=pv_p)
})

spearman_results <- spearman_results %>%
  mutate(variable = as.character(variable)) %>%
  left_join(label_df, by = "variable") %>%
  mutate(Variable = coalesce(Variable, variable),
         `Spearman `          = ifelse(is.na(rho), NA, sprintf("%.2f", rho)),
         `p-value`                  = ifelse(is.na(p_value), NA, format.pval(p_value, digits = 3, eps = 1e-4)),
         `Partial  (age)` = ifelse(is.na(rho_partial), NA, sprintf("%.2f", rho_partial)),
         `Partial p`                = ifelse(is.na(p_partial), NA, format.pval(p_partial, digits = 3, eps = 1e-4))) %>%
  select(Variable, n, `Spearman `, `p-value`, `Partial  (age)`, `Partial p`)

# Save Table 4
write_csv(spearman_results, file.path(out_dir, "Table4_Spearman_EDSS.csv"))

ft4 <- flextable(spearman_results) |>
  autofit() |>
  align(align = "left", j = 1) |>
  bold(j = 1) |>
  set_caption("Table 4. Spearman correlations between EDSS (numeric) and clinical/biochemical variables; optional partial correlations adjust for age and sex.")
save_as_docx(`Table 4` = ft4, path = file.path(out_dir, "Table4_Spearman_EDSS.docx"))

message("✅ Exported Table 3 and Table 4 (CSV + DOCX) to: ", out_dir)
if (!jt_available) message("ℹ Install trend test package: install.packages('clinfun')")
