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
# Subsection 2 — EDSS Associations (Q1-ready)
# Table 3: Kruskal–Wallis + Jonckheere–Terpstra (trend) by EDSS category
# Table 4: Spearman correlations (crude) + partial (age±sex)
# - p-values formatted to 4 decimals; "< 0.0001" for very small p
# - Direction arrows ↑ / ↓
# - Exports CSV + Q1-styled DOCX
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
  library(rstatix)
  library(clinfun)     # jonckheere.test
  library(flextable)
  library(officer)
  library(stringr)
})

# --------------------------- Paths ------------------------------------
# Edit base_dir if needed; by default use current working directory.
base_dir <- getwd()
data_file <- file.path(base_dir, "Libya_MS.csv")
out_dir   <- file.path(base_dir, "outputs", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_file)) {
  stop("Dataset not found at: ", data_file,
       "\nPlace 'Libya_MS.csv' in the working folder or update 'data_file'.")
}

# ------------------------- Load & clean --------------------------------
df <- read_csv(data_file, show_col_types = FALSE) %>% clean_names()
names(df) <- gsub("\\.", "_", names(df))

# Helpers
col_or <- function(...) {
  cands <- c(...)
  cands <- cands[!is.na(cands)]
  hits  <- cands[cands %in% names(df)]
  if (length(hits)) hits[1] else NA_character_
}
to_num <- function(x) suppressWarnings(as.numeric(x))
p_fmt4 <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 1e-4) "< 0.0001" else sprintf("%.4f", p)
}

# ---------------------- EDSS numeric + groups --------------------------
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
    edss_cat = factor(edss_cat,
                      levels = c("Mild (0–2.5)", "Moderate (3–5.5)", "Severe (≥6)"),
                      ordered = TRUE)
  )

# -------------------- Variables to analyze -----------------------------
vars <- list(
  age            = col_or("age"),
  age_onset      = col_or("age_of_onset_of_disease_years","age_at_onset","age_at_onset_years"),
  duration_years = col_or("duration_of_illness_years","duration_years"),
  arr            = col_or("arr","annualized_relapse_rate"),
  bmi            = col_or("bmi_kg_m2","bmi"),
  vit_d          = col_or("vitamin_d","vit_d"),
  tc             = col_or("total_cholesterol","tc"),
  ldl            = col_or("ldl"),
  hdl            = col_or("hdl"),
  tg             = col_or("tg","triglycerides"),
  chr            = col_or("chol_hdl_ratio","total_cholesterol_to_hdl_ratio"),
  lhr            = col_or("ldl_hdl_ratio"),
  aip            = col_or("aip")
)
vars <- vars[!is.na(unlist(vars))]

# labels for display
var_labels <- c(
  age="Age (years)", age_onset="Age at onset (years)", duration_years="Disease duration (years)",
  arr="Annualized relapse rate (ARR)", bmi="BMI (kg/m²)",
  vit_d="25(OH) Vitamin D (ng/mL)", tc="Total cholesterol (mg/dL)",
  ldl="LDL (mg/dL)", hdl="HDL (mg/dL)", tg="Triglycerides (mg/dL)",
  chr="Total cholesterol/HDL ratio", lhr="LDL/HDL ratio", aip="Atherogenic Index of Plasma (AIP)"
)
label_df <- tibble(variable = names(var_labels), Variable = unname(var_labels))

# numeric copy
df_num <- df %>% mutate(across(all_of(unlist(vars)), to_num))

# ======================================================================
#                              TABLE 3
# ======================================================================
kw_jt_results <- map_dfr(names(vars), function(vname) {
  vcol <- vars[[vname]]
  x <- df_num[[vcol]]
  g <- df_num$edss_cat
  ok <- is.finite(x) & !is.na(g)
  if (!any(ok)) {
    return(tibble(variable=vname, n=0, kw_p=NA_real_, jt_p=NA_real_, direction=NA_character_))
  }
  # Kruskal–Wallis (omnibus)
  kw <- tryCatch(rstatix::kruskal_test(data.frame(x=x[ok], g=g[ok]), x ~ g), error=function(e) NULL)
  kw_p <- if (is.null(kw)) NA_real_ else kw$p
  
  # Jonckheere–Terpstra (monotonic trend)
  jt <- tryCatch(clinfun::jonckheere.test(x[ok], as.numeric(g[ok])), error=function(e) NULL)
  jt_p <- if (is.null(jt)) NA_real_ else unname(jt$p.value)
  
  # Direction via Spearman vs ordered levels
  rho <- suppressWarnings(cor(x[ok], as.numeric(g[ok]), method="spearman"))
  dir <- case_when(
    is.na(rho)     ~ NA_character_,
    rho >  0.05    ~ "\u2191 with EDSS",  # ↑
    rho < -0.05    ~ "\u2193 with EDSS",  # ↓
    TRUE           ~ "flat"
  )
  
  tibble(variable=vname, n=sum(ok), kw_p=kw_p, jt_p=jt_p, direction=dir)
})

table3 <- kw_jt_results %>%
  left_join(label_df, by="variable") %>%
  mutate(
    Variable = coalesce(Variable, variable),
    `Kruskal–Wallis p value` = map_chr(kw_p, p_fmt4),
    `Trend p (Jonckheere–Terpstra)` = map_chr(jt_p, p_fmt4),
    Direction = direction
  ) %>%
  select(Variable, n, `Kruskal–Wallis p value`, `Trend p (Jonckheere–Terpstra)`, Direction)

# save CSV
write_csv(table3, file.path(out_dir, "Table3_EDSS_group_tests_Q1.csv"))

# DOCX (Q1 style)
is_sig <- function(pstr) {
  # treat "< 0.0001" as significant
  if (is.na(pstr)) return(FALSE)
  if (grepl("<", pstr)) return(TRUE)
  suppressWarnings(as.numeric(pstr)) < 0.05
}

ft3 <- flextable(table3) |>
  autofit() |>
  align(j = 2:ncol(table3), align = "center") |>
  bold(i = which(sapply(table3$`Kruskal–Wallis p value`, is_sig)),
       j = "Kruskal–Wallis p value") |>
  bold(i = which(sapply(table3$`Trend p (Jonckheere–Terpstra)`, is_sig)),
       j = "Trend p (Jonckheere–Terpstra)") |>
  set_caption("Table 3. Group-wise EDSS comparisons with Kruskal–Wallis (omnibus) and Jonckheere–Terpstra (monotonic trend) tests; p-values shown to four decimals.") |>
  theme_booktabs()

save_as_docx(`Table 3` = ft3, path = file.path(out_dir, "Table3_EDSS_group_tests_Q1.docx"))

# ======================================================================
#                              TABLE 4
# ======================================================================
sex_col <- col_or("gender","sex")
has_sex <- !is.na(sex_col)

spearman_results <- map_dfr(names(vars), function(vname) {
  vcol <- vars[[vname]]
  x <- df_num[[vcol]]
  y <- df_num$edss_num
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(tibble(variable=vname, n=0, rho=NA_real_, pv=NA_real_, rho_p=NA_real_, pv_p=NA_real_))
  
  # crude spearman
  cr <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman", exact=FALSE))
  rho <- unname(cr$estimate); pv <- cr$p.value; nn <- sum(ok)
  
  # partial (age±sex) via residuals if age present
  age_col <- vars[["age"]] %||% NA_character_
  if (!is.na(age_col)) {
    df_adj <- tibble(
      x = x, y = y,
      age = to_num(df_num[[age_col]]),
      sex = if (has_sex) factor(df_num[[sex_col]]) else NULL
    )
    covars <- df_adj %>% select(age, any_of("sex"))
    get_res <- function(resp, covars) {
      mod <- tryCatch(lm(resp ~ ., data=covars), error=function(e) NULL)
      if (is.null(mod)) rep(NA_real_, nrow(covars)) else resid(mod)
    }
    rx <- get_res(df_adj$x, covars); ry <- get_res(df_adj$y, covars)
    okp <- is.finite(rx) & is.finite(ry)
    pr <- suppressWarnings(cor.test(rx[okp], ry[okp], method="spearman", exact=FALSE))
    rho_p <- unname(pr$estimate); pv_p <- pr$p.value
  } else {
    rho_p <- NA_real_; pv_p <- NA_real_
  }
  
  tibble(variable=vname, n=nn, rho=rho, pv=pv, rho_p=rho_p, pv_p=pv_p)
})

table4 <- spearman_results %>%
  left_join(label_df, by="variable") %>%
  mutate(
    Variable = coalesce(Variable, variable),
    `Spearman rho` = ifelse(is.na(rho), NA_character_, sprintf("%.2f", rho)),
    `p value` = map_chr(pv, p_fmt4),
    `Partial rho (age ± sex)` = ifelse(is.na(rho_p), NA_character_, sprintf("%.2f", rho_p)),
    `Partial p value` = map_chr(pv_p, p_fmt4)
  ) %>%
  select(Variable, n, `Spearman rho`, `p value`, `Partial rho (age ± sex)`, `Partial p value`)

write_csv(table4, file.path(out_dir, "Table4_Spearman_EDSS_Q1.csv"))

ft4 <- flextable(table4) |>
  autofit() |>
  align(j = 2:ncol(table4), align = "center") |>
  bold(i = which(sapply(table4$`p value`, is_sig)),           j = "p value") |>
  bold(i = which(sapply(table4$`Partial p value`, is_sig)),   j = "Partial p value") |>
  set_caption("Table 4. Spearman correlations between EDSS and clinical/biochemical variables; partial rho adjusts for age and sex. P-values shown to four decimals.") |>
  theme_booktabs()

save_as_docx(`Table 4` = ft4, path = file.path(out_dir, "Table4_Spearman_EDSS_Q1.docx"))

message("Done. Files saved in: ", out_dir,
        "\n- Table3_EDSS_group_tests_Q1.csv / .docx",
        "\n- Table4_Spearman_EDSS_Q1.csv / .docx")
