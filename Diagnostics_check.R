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
# Subsection 3.6 — Model Diagnostics & Robustness Checks
# - Proportional-odds assumption (Brant; fallback to ordinal::nominal_test)
# - Multicollinearity (custom VIF on design matrix)
# - Fit quality (Nagelkerke R²) + partial R² via drop-one ΔR²
# - Calibration plot & residual diagnostics
# - Exports a compact Supplementary Table S1 (DOCX) + figures
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(janitor)
  library(broom)
  library(performance)    # r2_nagelkerke()
  library(officer)
  library(flextable)
})

# ---------- Paths ----------
base_dir <- getwd()
data_file <- file.path(base_dir, "Libya_MS.csv")
out_dir   <- file.path(base_dir, "outputs", "subsection_diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(data_file))

# ---------- Load & normalize ----------
df <- read_csv(data_file, show_col_types = FALSE) %>% clean_names()
names(df) <- gsub("\\.", "_", names(df))

to_num <- function(x) suppressWarnings(as.numeric(x))
col_or <- function(df, ...) {cands <- c(...); hit <- cands[cands %in% names(df)]; if (length(hit)) hit[1] else NA_character_}

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
if (any(is.na(needed))) stop("Missing required variable(s). Check column names.")

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
  ) %>%
  select(-sex_raw) %>%
  filter(is.finite(edss))

# EDSS categories
model_df <- model_df %>%
  mutate(
    edss_cat = case_when(
      edss <= 2.5 ~ "Mild (0–2.5)",
      edss > 2.5 & edss <= 5.5 ~ "Moderate (3–5.5)",
      edss > 5.5 ~ "Severe (≥6)"
    ),
    edss_cat = factor(edss_cat, levels = c("Mild (0–2.5)","Moderate (3–5.5)","Severe (≥6)"), ordered = TRUE)
  ) %>%
  drop_na(edss_cat, age, duration, bmi, vitamin_d, hdl, ldl_hdl, aip, sex)

# Standardize continuous predictors
scale_cols <- c("age","duration","bmi","vitamin_d","hdl","ldl_hdl","aip")
model_df[scale_cols] <- lapply(model_df[scale_cols], function(z) as.numeric(scale(z)))

# ---------- Fit ordinal model ----------
form <- as.formula(edss_cat ~ age + duration + bmi + vitamin_d + hdl + ldl_hdl + aip + sex)
fit  <- suppressWarnings(MASS::polr(formula = form, data = model_df, Hess = TRUE))

# ---------- Proportional-odds assumption (Brant / ordinal fallback) ----------
po_test <- NULL
po_method <- NULL

if (requireNamespace("brant", quietly = TRUE)) {
  tmp <- try(brant::brant(fit), silent = TRUE)
  if (!inherits(tmp, "try-error")) {
    po_test <- list(statistic = unname(tmp$"Test statistic"[1]), p = unname(tmp$"p value"[1]))
    po_method <- "Brant (polr)"
  }
}
if (is.null(po_test) && requireNamespace("ordinal", quietly = TRUE)) {
  # Refit with ordinal::clm and run nominal_test as a PO check
  fit_clm <- ordinal::clm(form, data = model_df, link = "logit")
  nt <- ordinal::nominal_test(fit_clm)  # tests non-parallelism
  # Combine p-values (max p as global acceptance proxy)
  p_glob <- tryCatch(max(nt[["Pr(>Chi)"]], na.rm = TRUE), error = function(e) NA_real_)
  po_test <- list(statistic = NA_real_, p = p_glob)
  po_method <- "ordinal::nominal_test (CLM fallback)"
}

# ---------- VIF (custom on design matrix) ----------
# Build design matrix used by the model (no intercept)
X <- model.matrix(fit)  # includes intercept col "(Intercept)"
X <- X[, colnames(X) != "(Intercept)", drop = FALSE]

compute_vif <- function(Xmat) {
  vifs <- numeric(ncol(Xmat))
  for (j in seq_len(ncol(Xmat))) {
    yj <- Xmat[, j]
    Xj <- Xmat[, -j, drop = FALSE]
    r2 <- tryCatch(summary(lm(yj ~ Xj))$r.squared, error = function(e) NA_real_)
    vifs[j] <- if (is.na(r2) || r2 >= 0.9999) Inf else 1 / (1 - r2)
  }
  tibble::tibble(term = colnames(Xmat), VIF = vifs)
}
vif_tab <- compute_vif(X)
readr::write_csv(vif_tab, file.path(out_dir, "diag_vif.csv"))

# ---------- Nagelkerke R² ----------
nk_obj <- performance::r2_nagelkerke(fit)
nk_val <- if (is.list(nk_obj) && !is.null(nk_obj$R2)) nk_obj$R2 else as.numeric(nk_obj)

# ---------- Partial R² via drop-one ΔNagelkerke R² ----------
drop_one_r2 <- function(full_model, data, drop_term, base_form) {
  red_form <- update(base_form, paste(". ~ . -", drop_term))
  red <- suppressWarnings(MASS::polr(red_form, data = data, Hess = TRUE))
  full_r2 <- as.numeric(performance::r2_nagelkerke(full_model))
  red_r2  <- as.numeric(performance::r2_nagelkerke(red))
  max(0, full_r2 - red_r2)
}
terms_for_r2 <- c("age","duration","bmi","vitamin_d","hdl","ldl_hdl","aip","sex")
part_vals <- vapply(terms_for_r2, function(t) drop_one_r2(fit, model_df, t, form), numeric(1))
part_df <- tibble::tibble(
  Parameter = c("Age","Disease duration","BMI","25(OH) Vitamin D","HDL","LDL/HDL ratio","AIP","Sex (Female)"),
  Partial_R2 = round(part_vals, 3)
) %>% arrange(desc(Partial_R2))
readr::write_csv(part_df, file.path(out_dir, "diag_partialR2.csv"))

# ---------- Calibration plot ----------
# Predicted class probabilities
pred_probs <- as_tibble(predict(fit, type = "probs"))
colnames(pred_probs) <- levels(model_df$edss_cat)
lp <- as.numeric(predict(fit, type = "link"))  # linear predictor
cal_df <- bind_cols(model_df %>% select(edss_cat), pred_probs) %>%
  mutate(bin = ntile(lp, 10)) %>%
  group_by(bin) %>%
  summarise(
    n = n(),
    obs_mild     = mean(edss_cat == levels(edss_cat)[1]),
    obs_moderate = mean(edss_cat == levels(edss_cat)[2]),
    obs_severe   = mean(edss_cat == levels(edss_cat)[3]),
    pred_mild     = mean(`Mild (0–2.5)`),
    pred_moderate = mean(`Moderate (3–5.5)`),
    pred_severe   = mean(`Severe (≥6)`)
  ) %>%
  ungroup() %>%
  pivot_longer(-c(bin,n), names_to = c("type","class"), names_sep = "_", values_to = "prop") %>%
  mutate(type = ifelse(type == "obs", "Observed", "Predicted"),
         class = factor(class, levels = c("mild","moderate","severe"),
                        labels = c("Mild (0–2.5)","Moderate (3–5.5)","Severe (≥6)")))

p_cal <- ggplot(cal_df, aes(x = bin, y = prop, color = type, group = type)) +
  geom_line(size = 1) + geom_point() +
  facet_wrap(~ class, nrow = 1) +
  scale_color_manual(values = c("Observed" = "#1f78b4", "Predicted" = "#e31a1c")) +
  labs(x = "Risk decile (by linear predictor)", y = "Proportion",
       title = "Calibration of ordinal model across EDSS categories") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "Figure_S2_Calibration.png"), p_cal, width = 12, height = 4.2, dpi = 400)

# ---------- Residual diagnostics ----------
res_dev <- residuals(fit, type = "deviance")
res_std <- as.numeric(scale(res_dev))
p_hist <- ggplot(data.frame(res_std), aes(x = res_std)) +
  geom_histogram(bins = 30, fill = "#6baed6") +
  labs(title = "Standardized deviance residuals", x = "Residual", y = "Count") +
  theme_minimal(base_size = 13)

qq <- qqnorm(res_std, plot.it = FALSE)
p_qq <- ggplot(data.frame(x = qq$x, y = qq$y), aes(x, y)) +
  geom_point(color = "#636363") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "QQ plot of standardized residuals", x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal(base_size = 13)

p_res <- ggpubr::ggarrange(p_hist, p_qq, ncol = 2, labels = c("A", "B"))
ggsave(file.path(out_dir, "Figure_S3_Residuals.png"), p_res, width = 10, height = 4.2, dpi = 400)

# ---------- Summarize for Supplementary Table S1 ----------
po_p <- if (is.null(po_test)) NA_real_ else po_test$p
po_stat <- if (is.null(po_test)) NA_real_ else po_test$statistic
po_row <- tibble::tibble(
  Item = "Proportional-odds assumption",
  Method = ifelse(is.null(po_method), "Not available", po_method),
  Statistic = sprintf("%.2f", ifelse(is.na(po_stat), NA, po_stat)),
  `p value` = ifelse(is.na(po_p), "NA", sprintf("%.4f", po_p))
)

vif_summary <- vif_tab %>%
  summarise(`VIF min` = min(VIF, na.rm = TRUE),
            `VIF max` = max(VIF, na.rm = TRUE),
            `VIF median` = median(VIF, na.rm = TRUE))

fit_row <- tibble::tibble(
  Item = "Model fit (Nagelkerke R²)",
  Method = "performance::r2_nagelkerke",
  Statistic = "",
  `p value` = sprintf("%.3f", nk_val)
)

# Partial R² top contributors (sorted already)
part_top <- part_df

# Build DOCX
doc <- read_docx() %>%
  body_add_par(value = "Supplementary Table S1. Model diagnostics and robustness checks", style = "heading 1") %>%
  body_add_flextable(
    flextable(bind_rows(po_row, fit_row)) %>%
      autofit()
  ) %>%
  body_add_par(value = "") %>%
  body_add_par(value = "Multicollinearity summary (design-matrix VIFs):", style = "heading 2") %>%
  body_add_flextable(flextable(vif_summary) %>% autofit()) %>%
  body_add_par(value = "") %>%
  body_add_par(value = "Predictor importance (partial R² from drop-one ΔR²):", style = "heading 2") %>%
  body_add_flextable(
    flextable(part_top %>% rename(`Partial R²` = Partial_R2)) %>% autofit()
  ) %>%
  body_add_par(value = "") %>%
  body_add_par(value = "Calibration and residual diagnostics are provided in Figure S2 and Figure S3.", style = "Normal")

print(doc, target = file.path(out_dir, "Diag_S1_ModelDiagnostics.docx"))

message("\n✅ Diagnostics written to: ", out_dir,
        "\n - Diag_S1_ModelDiagnostics.docx",
        "\n - Figure_S2_Calibration.png",
        "\n - Figure_S3_Residuals.png",
        "\n - diag_vif.csv, diag_partialR2.csv")
