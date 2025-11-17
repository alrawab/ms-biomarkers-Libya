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
# Subsection 3.4 — Multivariable modeling & predictors of disability
# Outputs:
#   - Table5_Ordinal_Logistic_EDSS.(csv|docx)
#   - Table6_PartialR2_EDSS.(csv|docx)
#   - Figure3_MS_Disability_Predictors.(png|tiff)
#   - model_diagnostics.txt
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(janitor)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(scales)
  library(broom)
  library(performance)     # r2_nagelkerke
  library(MASS)            # polr
  library(flextable)       # tables to docx
})

# ---------- Paths ----------
base_dir <- getwd()
data_file <- file.path(base_dir, "Libya_MS.csv")
out_dir   <- file.path(base_dir, "outputs", "subsection3_4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_file)) {
  stop("Cannot find 'Libya_MS.csv' at: ", data_file,
       "\nUpdate 'data_file' or set your working directory correctly.")
}

# ---------- Load & clean ----------
df <- readr::read_csv(data_file, show_col_types = FALSE) |>
  janitor::clean_names()
names(df) <- gsub("\\.", "_", names(df))

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
col_or <- function(df, ...) {
  cands <- c(...)
  hit <- cands[cands %in% names(df)]
  if (length(hit)) hit[1] else NA_character_
}

# Resolve column names robustly
edss_col <- col_or(df, "edss")
age_col  <- col_or(df, "age")
dur_col  <- col_or(df, "duration_of_illness_years","duration_years")
bmi_col  <- col_or(df, "bmi_kg_m2","bmi")
vitd_col <- col_or(df, "vitamin_d","vit_d")
hdl_col  <- col_or(df, "hdl")
lhr_col  <- col_or(df, "ldl_hdl_ratio")
aip_col  <- col_or(df, "aip")
sex_col  <- col_or(df, "gender","sex")

needed <- c(edss_col, age_col, dur_col, bmi_col, vitd_col, hdl_col, lhr_col, aip_col, sex_col)
if (any(is.na(needed))) {
  miss <- c("edss"=edss_col, "age"=age_col, "duration"=dur_col, "bmi"=bmi_col,
            "vitamin_d"=vitd_col, "hdl"=hdl_col, "ldl_hdl_ratio"=lhr_col,
            "aip"=aip_col, "sex/gender"=sex_col)
  stop("Missing required columns: ", paste(names(miss)[is.na(miss)], collapse = ", "))
}

# Prepare modeling dataframe
model_df <- df %>%
  transmute(
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
  mutate(
    sex = case_when(
      is.na(sex_raw) ~ NA_character_,
      tolower(as.character(sex_raw)) %in% c("male","m","1")   ~ "Male",
      tolower(as.character(sex_raw)) %in% c("female","f","0") ~ "Female",
      TRUE ~ as.character(sex_raw)
    ),
    sex = factor(sex)
  )

model_df$sex_raw <- NULL
model_df <- model_df |> filter(is.finite(edss))

# Build ordered EDSS categories
model_df <- model_df %>%
  mutate(
    edss_cat = case_when(
      edss <= 2.5 ~ "Mild (0–2.5)",
      edss > 2.5 & edss <= 5.5 ~ "Moderate (3–5.5)",
      edss > 5.5 ~ "Severe (≥6)"
    ),
    edss_cat = factor(edss_cat,
                      levels = c("Mild (0–2.5)", "Moderate (3–5.5)", "Severe (≥6)"),
                      ordered = TRUE)
  ) %>%
  tidyr::drop_na(edss_cat, age, duration, bmi, vitamin_d, hdl, ldl_hdl, aip, sex)

# Standardize continuous predictors (per SD) for interpretable ORs
scale_cols <- c("age","duration","bmi","vitamin_d","hdl","ldl_hdl","aip")
model_df[scale_cols] <- lapply(model_df[scale_cols], function(z) as.numeric(scale(z)))

# ---------- Fit proportional-odds model ----------
form <- as.formula(edss_cat ~ age + duration + bmi + vitamin_d + hdl + ldl_hdl + aip + sex)
fit <- suppressWarnings(MASS::polr(formula = form, data = model_df, Hess = TRUE))

# ---------- Odds ratios with 95% CIs (Table 5) ----------
ci_prof <- suppressWarnings(confint(fit))
or_vec  <- exp(coef(fit))
or_ci   <- exp(ci_prof)

or_tab <- tibble(
  term = names(or_vec),
  OR   = as.numeric(or_vec),
  CI_lo = as.numeric(or_ci[,1]),
  CI_hi = as.numeric(or_ci[,2])
) |> filter(!grepl("\\|", term))

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

table5 <- or_tab %>%
  mutate(Predictor = term_map[term],
         `OR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", OR, CI_lo, CI_hi)) %>%
  select(Predictor, `OR (95% CI)`) %>%
  arrange(Predictor)

write_csv(table5, file.path(out_dir, "Table5_Ordinal_Logistic_EDSS.csv"))

ft5 <- flextable(table5)
ft5 <- autofit(ft5)
flextable::save_as_docx("Table 5. Adjusted odds ratios for EDSS categories" = ft5,
                        path = file.path(out_dir, "Table5_Ordinal_Logistic_EDSS.docx"))

# ---------- Model fit & partial R2 (Table 6) ----------
safe_nagelkerke <- function(m) {
  rr <- try(performance::r2_nagelkerke(m), silent = TRUE)
  if (inherits(rr, "try-error")) return(NA_real_)
  # handle possible structures
  v <- try(suppressWarnings(unlist(rr)), silent = TRUE)
  if (!inherits(v,"try-error")) {
    # prefer element with 'Nagelkerke' in name
    idx <- which(grepl("Nagelkerke", names(v), ignore.case = TRUE))
    if (length(idx)) return(as.numeric(v[idx[1]]))
    return(as.numeric(v[1]))
  }
  NA_real_
}

full_r2 <- safe_nagelkerke(fit)

drop_one <- function(term) {
  red <- suppressWarnings(MASS::polr(update(form, paste(". ~ . -", term)),
                                     data = model_df, Hess = TRUE))
  r2 <- safe_nagelkerke(red)
  max(0, full_r2 - r2)
}

terms_for_r2 <- c("age","duration","bmi","vitamin_d","hdl","ldl_hdl","aip","sex")
part_vals <- sapply(terms_for_r2, drop_one)

table6 <- tibble(
  Parameter = c("Age","Disease duration","BMI","25(OH) Vitamin D",
                "HDL","LDL/HDL ratio","AIP","Sex (Female)"),
  `Partial R2` = round(as.numeric(part_vals), 3)
) |> arrange(desc(`Partial R2`))

write_csv(table6, file.path(out_dir, "Table6_PartialR2_EDSS.csv"))

ft6 <- flextable(table6)
ft6 <- autofit(ft6)
flextable::save_as_docx("Table 6. Relative contribution to model fit (partial R²)" = ft6,
                        path = file.path(out_dir, "Table6_PartialR2_EDSS.docx"))

# ---------- Figure 3A (forest of ORs) ----------
forest_df <- or_tab %>%
  transmute(
    Predictor = factor(term_map[term],
                       levels = rev(term_map[or_tab$term])),
    OR = OR, lo = CI_lo, hi = CI_hi
  )

pA <- ggplot(forest_df, aes(y = Predictor, x = OR)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, size = 0.6, color = "#3B5BA5") +
  geom_point(size = 3, color = "#3B5BA5") +
  scale_x_log10(labels = label_number(accuracy = 0.01)) +
  labs(title = "A) Independent predictors of EDSS category",
       x = "Adjusted Odds Ratio (log scale)", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# ---------- Figure 3B (partial R2 bars) ----------
pB <- ggplot(table6, aes(x = reorder(Parameter, `Partial R2`), y = `Partial R2`)) +
  geom_col(fill = "#E15759") +
  coord_flip() +
  labs(title = "B) Relative contribution to model fit (Partial R²)",
       x = NULL, y = "Partial R²") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

fig3 <- ggpubr::ggarrange(pA, pB, ncol = 2, widths = c(1.2, 1))

ggsave(file.path(out_dir, "Figure3_MS_Disability_Predictors.png"),
       fig3, dpi = 600, width = 12, height = 6)
ggsave(file.path(out_dir, "Figure3_MS_Disability_Predictors.tiff"),
       fig3, dpi = 600, width = 12, height = 6, compression = "lzw")

# ---------- Diagnostics ----------
has_brant <- requireNamespace("brant", quietly = TRUE)
brant_text <- if (has_brant) {
  capture.output(try(brant::brant(fit), silent = TRUE))
} else {
  "Package 'brant' not installed; proportional-odds assumption not formally tested here."
}

diag_lines <- c(
  "Model: polr(EDSS category ~ age + duration + bmi + vitamin_d + hdl + ldl_hdl + aip + sex)",
  paste0("Nagelkerke R^2 (full model): ", round(full_r2, 3)),
  "",
  "Brant test (if available):",
  brant_text
)
writeLines(diag_lines, con = file.path(out_dir, "model_diagnostics.txt"))

message("\n Subsection 3.4 outputs saved to: ", out_dir,
        "\n - Table5_Ordinal_Logistic_EDSS.(csv|docx)",
        "\n - Table6_PartialR2_EDSS.(csv|docx)",
        "\n - Figure3_MS_Disability_Predictors.(png|tiff)",
        "\n - model_diagnostics.txt")
