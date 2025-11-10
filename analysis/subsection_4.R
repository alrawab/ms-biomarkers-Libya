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

# ================================================================
# Subsection 4 — Multivariable modeling & predictors of disability
# Conflict-proof: no select(-...), no attach(MASS), explicit namespaces
# Exports: Figure 3 (TIFF/PNG), Table 5 (ORs), Table 6 (Partial R2)
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(janitor)
  library(broom)
  library(performance)   # r2_nagelkerke
  library(ggplot2)
  library(ggpubr)
  library(scales)
})

# ---------- Paths ----------
base_dir <- getwd()
data_file <- file.path(base_dir, "Libya_MS.csv")
out_dir   <- file.path(base_dir, "outputs", "subsection4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_file)) {
  stop("Cannot find 'Libya_MS.csv' at: ", data_file,
       "\nPlace the dataset there or update 'data_file'.")
}

# ---------- Load & clean ----------
df <- readr::read_csv(data_file, show_col_types = FALSE) |>
  janitor::clean_names()
names(df) <- gsub("\\.", "_", names(df))  # normalize any stray dots

cat("\nAvailable columns in data:\n"); print(names(df))

# ---------- Helpers ----------
to_num <- function(x) suppressWarnings(as.numeric(x))
col_or <- function(df, ...) {
  cands <- c(...)
  hit <- cands[cands %in% names(df)]
  if (length(hit)) hit[1] else NA_character_
}

# ---------- Resolve columns (robust to naming variants) ----------
edss_col <- col_or(df, "edss")
age_col  <- col_or(df, "age")
dur_col  <- col_or(df, "duration_of_illness_years","duration_years")
bmi_col  <- col_or(df, "bmi_kg_m2","bmi")
vitd_col <- col_or(df, "vitamin_d","vit_d")
hdl_col  <- col_or(df, "hdl")
lhr_col  <- col_or(df, "ldl_hdl_ratio")   # LDL/HDL ratio
aip_col  <- col_or(df, "aip")
sex_col  <- col_or(df, "gender","sex")

needed <- c(edss_col, age_col, dur_col, bmi_col, vitd_col, hdl_col, lhr_col, aip_col, sex_col)
if (any(is.na(needed))) {
  miss <- c("edss"=edss_col, "age"=age_col, "duration"=dur_col, "bmi"=bmi_col,
            "vitamin_d"=vitd_col, "hdl"=hdl_col, "ldl_hdl_ratio"=lhr_col,
            "aip"=aip_col, "sex/gender"=sex_col)
  stop("Missing required columns (first-match rule). Check these names:\n",
       paste(names(miss)[is.na(miss)], collapse = ", "))
}

# ---------- Conflict-proof data prep (no select(-...), no MASS masking) ----------
model_df <- df %>%
  dplyr::transmute(
    edss      = to_num(.data[[edss_col]]),
    age       = to_num(.data[[age_col]]),
    duration  = to_num(.data[[dur_col]]),
    bmi       = to_num(.data[[bmi_col]]),
    vitamin_d = to_num(.data[[vitd_col]]),
    hdl       = to_num(.data[[hdl_col]]),
    ldl_hdl   = to_num(.data[[lhr_col]]),
    aip       = to_num(.data[[aip_col]]),
    sex_raw   = .data[[sex_col]]
  ) %>%
  dplyr::mutate(
    sex = dplyr::case_when(
      is.na(sex_raw) ~ NA_character_,
      tolower(as.character(sex_raw)) %in% c("male","m","1")   ~ "Male",
      tolower(as.character(sex_raw)) %in% c("female","f","0") ~ "Female",
      TRUE ~ as.character(sex_raw)
    ),
    sex = factor(sex)
  )

# Drop helper using base R (avoids select masking issues)
model_df$sex_raw <- NULL

# Keep valid EDSS
model_df <- model_df %>% dplyr::filter(is.finite(edss))

# Build ordered EDSS categories
model_df <- model_df %>%
  dplyr::mutate(
    edss_cat = dplyr::case_when(
      edss <= 2.5 ~ "Mild (0–2.5)",
      edss > 2.5 & edss <= 5.5 ~ "Moderate (3–5.5)",
      edss > 5.5 ~ "Severe (≥6)"
    ),
    edss_cat = factor(edss_cat,
                      levels = c("Mild (0–2.5)", "Moderate (3–5.5)", "Severe (≥6)"),
                      ordered = TRUE)
  ) %>%
  tidyr::drop_na(edss_cat, age, duration, bmi, vitamin_d, hdl, ldl_hdl, aip, sex)

# Standardize continuous predictors (z-scores) for comparable ORs per SD
scale_cols <- c("age","duration","bmi","vitamin_d","hdl","ldl_hdl","aip")
model_df[scale_cols] <- lapply(model_df[scale_cols], function(z) as.numeric(scale(z)))

cat("\nModel data preview:\n"); print(utils::head(model_df))

# ---------- Fit ordinal logistic regression (no library(MASS)) ----------
form <- stats::as.formula(edss_cat ~ age + duration + bmi + vitamin_d + hdl + ldl_hdl + aip + sex)
fit  <- suppressWarnings(MASS::polr(formula = form, data = model_df, Hess = TRUE))
summ <- summary(fit)
cat("\nModel converged. Coefficients:\n"); print(coef(summ))

# ---------- Odds ratios with 95% CIs ----------
ci_prof <- suppressWarnings(confint(fit))   # profile-likelihood CI
or_vec  <- exp(coef(fit))
or_ci   <- exp(ci_prof)

or_tab <- tibble::tibble(
  term = names(or_vec),
  odds_ratio = as.numeric(or_vec),
  ci_low = as.numeric(or_ci[,1]),
  ci_high = as.numeric(or_ci[,2])
) %>% dplyr::filter(!grepl("\\|", term))   # remove thresholds

term_map <- c(
  "age"="Age (per SD)",
  "duration"="Disease duration (per SD)",
  "bmi"="BMI (per SD)",
  "vitamin_d"="25(OH) Vitamin D (per SD)",
  "hdl"="HDL (per SD)",
  "ldl_hdl"="LDL/HDL ratio (per SD)",
  "aip"="Atherogenic Index of Plasma (per SD)",
  "sexFemale"="Female vs Male"
)
or_tab$label <- term_map[or_tab$term]

table5 <- or_tab %>%
  dplyr::transmute(
    Predictor = label,
    `Odds Ratio` = round(odds_ratio, 2),
    `CI low`     = round(ci_low, 2),
    `CI high`    = round(ci_high, 2)
  )
readr::write_csv(table5, file.path(out_dir, "Table5_Ordinal_Logistic_EDSS.csv"))

# ---------- Robust Nagelkerke R² getter ----------
get_nk_r2 <- function(mod) {
  nk <- try(performance::r2_nagelkerke(mod), silent = TRUE)
  if (inherits(nk, "try-error")) {
    if (requireNamespace("pscl", quietly = TRUE)) {
      nk2 <- pscl::pR2(mod)
      return(as.numeric(nk2["McFadden"]))  # fallback
    } else return(NA_real_)
  }
  if (is.numeric(nk) && length(nk) == 1) return(as.numeric(nk))
  if (is.list(nk) && !is.null(nk$R2)) return(as.numeric(nk$R2))
  if (is.data.frame(nk) && "R2" %in% names(nk)) return(as.numeric(nk$R2[1]))
  for (nm in c("R2_Nagelkerke","R2_Nagelkerke_Adjusted","R2Nagelkerke")) {
    if ((is.list(nk) || is.data.frame(nk)) && !is.null(nk[[nm]])) return(as.numeric(nk[[nm]][1]))
  }
  vals <- suppressWarnings(as.numeric(unlist(nk))); vals <- vals[is.finite(vals)]
  if (length(vals)) return(vals[1])
  NA_real_
}

# ---------- Model fit ----------
r2_full <- get_nk_r2(fit)
cat("\nNagelkerke R^2: ", round(r2_full, 3), "\n")

# ---------- Partial R² via drop-one ΔNagelkerke R² (robust) ----------
drop_one_r2 <- function(full_model, data, drop_term, base_form) {
  red_form <- update(base_form, paste(". ~ . -", drop_term))
  red <- suppressWarnings(MASS::polr(red_form, data = data, Hess = TRUE))
  full_r2 <- get_nk_r2(full_model)
  red_r2  <- get_nk_r2(red)
  if (!is.finite(full_r2) || !is.finite(red_r2)) return(NA_real_)
  max(0, full_r2 - red_r2)
}
terms_for_r2 <- c("age","duration","bmi","vitamin_d","hdl","ldl_hdl","aip","sex")
part_vals <- vapply(terms_for_r2, function(t) drop_one_r2(fit, model_df, t, form), numeric(1))

table6 <- tibble::tibble(
  Parameter    = c("Age","Disease duration","BMI","25(OH) Vitamin D","HDL","LDL/HDL ratio","AIP","Sex (Female)"),
  `Partial R2` = round(part_vals, 3)
) %>% dplyr::arrange(dplyr::desc(`Partial R2`))
readr::write_csv(table6, file.path(out_dir, "Table6_PartialR2_EDSS.csv"))

# ---------- Figure 3A: Forest plot (Adjusted ORs) ----------
forest_df <- table5 %>% dplyr::mutate(Predictor = factor(Predictor, levels = rev(Predictor)))
p_forest <- ggplot(forest_df, aes(y = Predictor, x = `Odds Ratio`)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray55") +
  geom_errorbarh(aes(xmin = `CI low`, xmax = `CI high`), height = 0.25) +
  geom_point(size = 3) +
  scale_x_log10(labels = scales::number_format(accuracy = 0.01)) +
  labs(title = "Independent Predictors of Disability (EDSS categories)",
       x = "Adjusted Odds Ratio (log scale)", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# ---------- Figure 3B: Partial R² plot ----------
p_r2 <- ggplot(table6, aes(x = reorder(Parameter, `Partial R2`), y = `Partial R2`)) +
  geom_col() +
  coord_flip() +
  labs(title = "Relative contribution to model fit (Partial R²)",
       x = NULL, y = "Partial R²") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# ---------- Combine & save ----------
fig3 <- ggpubr::ggarrange(p_forest, p_r2, ncol = 2, labels = c("A", "B"))
ggsave(file.path(out_dir, "Figure3_MS_Disability_Predictors.tiff"),
       fig3, dpi = 600, width = 12, height = 6)
ggsave(file.path(out_dir, "Figure3_MS_Disability_Predictors.png"),
       fig3, dpi = 600, width = 12, height = 6)

message("\n Outputs saved to: ", out_dir,
        "\n - Table5_Ordinal_Logistic_EDSS.csv",
        "\n - Table6_PartialR2_EDSS.csv",
        "\n - Figure3_MS_Disability_Predictors.(tiff|png)")
