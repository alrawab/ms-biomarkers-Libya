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
# Subsection 3.3 — Correlation analysis with continuous EDSS
# Safe Spearman + Partial (age, sex) correlations
# Exports: Table4_Spearman_EDSS_Q1.docx, Figure2E-F_Correlations.png
# ================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(janitor)
  library(ppcor)
  library(flextable); library(officer)
  library(ggplot2); library(ggpubr); library(scales); library(tidyr)
})

# ---------- Paths ----------
base_dir <- getwd()
data_file <- file.path(base_dir, "Libya_MS.csv")
out_dir   <- file.path(base_dir, "outputs", "subsection3_3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- Load & clean ----------
df <- read_csv(data_file, show_col_types = FALSE) |> janitor::clean_names()
names(df) <- gsub("\\.", "_", names(df))

# Helper
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# Variables
vars <- c(
  "age",
  "age_of_onset_of_disease_years",
  "duration_of_illness_years",
  "arr",
  "bmi_kg_m2",
  "vitamin_d",
  "total_cholesterol",
  "ldl",
  "hdl",
  "tg",
  "chol_hdl_ratio",   # <- correct column name in your dataset
  "ldl_hdl_ratio",
  "aip"
)
edss_var <- "edss"
sex_var  <- "gender"

# ---------- Prepare analysis frame ----------
if (!edss_var %in% names(df)) stop("EDSS column not found.")

df2 <- df %>%
  mutate(
    across(all_of(intersect(vars, names(df))), to_num),
    edss = to_num(.data[[edss_var]]),
    sex = case_when(
      tolower(as.character(.data[[sex_var]])) %in% c("male","m","1") ~ 1,
      tolower(as.character(.data[[sex_var]])) %in% c("female","f","0") ~ 0,
      TRUE ~ NA_real_
    ),
    age = if ("age" %in% names(.)) to_num(age) else NA_real_
  ) %>%
  filter(is.finite(edss))

vars <- vars[vars %in% names(df2)]
if (length(vars) == 0) stop("None of the candidate variables were found after cleaning.")

# ---------- Safe correlation helpers ----------
safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(unique(x)) < 3 || length(unique(y)) < 3) return(list(rho=NA_real_, p=NA_real_))
  ct <- suppressWarnings(cor.test(x, y, method="spearman", exact=FALSE))
  list(rho = unname(ct$estimate), p = unname(ct$p.value))
}

safe_partial <- function(x, y, age, sex) {
  ok <- is.finite(x) & is.finite(y) & is.finite(age) & is.finite(sex)
  x <- x[ok]; y <- y[ok]; age <- age[ok]; sex <- sex[ok]
  if (length(x) < 10 ||
      length(unique(x)) < 3 || length(unique(y)) < 3 ||
      length(unique(age)) < 3 || length(unique(sex)) < 2) {
    return(list(rho=NA_real_, p=NA_real_))
  }
  # singularity guard
  if (any(apply(cbind(x,y,age,sex), 2, function(v) sd(v)==0))) return(list(rho=NA_real_, p=NA_real_))
  out <- try(ppcor::pcor.test(x, y, cbind(age, sex), method="spearman"), silent=TRUE)
  if (inherits(out, "try-error")) return(list(rho=NA_real_, p=NA_real_))
  list(rho = unname(out$estimate), p = unname(out$p.value))
}

# ---------- Compute correlations ----------
results <- lapply(vars, function(v) {
  x <- df2[[v]]; y <- df2$edss
  sp <- safe_spearman(x, y)
  pr <- safe_partial(x, y, df2$age, df2$sex)
  data.frame(
    Variable = v,
    Spearman_rho = sp$rho,
    Spearman_p   = sp$p,
    Partial_rho  = pr$rho,
    Partial_p    = pr$p,
    stringsAsFactors = FALSE
  )
}) |> bind_rows()

results_fmt <- results %>%
  mutate(
    Spearman_rho = round(Spearman_rho, 2),
    Partial_rho  = round(Partial_rho, 2),
    Spearman_p   = ifelse(is.na(Spearman_p), NA, formatC(Spearman_p, digits=4, format="f")),
    Partial_p    = ifelse(is.na(Partial_p), NA, formatC(Partial_p, digits=4, format="f"))
  ) %>%
  arrange(desc(abs(Spearman_rho)))

# ---------- Export Table (DOCX) ----------
ft <- flextable(results_fmt) |>
  set_header_labels(
    Variable     = "Variable",
    Spearman_rho = "Spearman ρ",
    Spearman_p   = "p value",
    Partial_rho  = "Partial ρ (age, sex)",
    Partial_p    = "Partial p value"
  ) |>
  autofit()

doc <- read_docx() %>%
  body_add_par("Table 4. Correlation between biochemical markers and EDSS (Spearman and partial)", style="heading 1") %>%
  body_add_flextable(ft)

print(doc, target=file.path(out_dir, "Table4_Spearman_EDSS_Q1.docx"))

# ---------- Visualization ----------
res_long <- results_fmt %>%
  mutate(Variable = gsub("_", " ", Variable)) %>%
  mutate(Variable = factor(Variable, levels = rev(Variable))) %>%
  pivot_longer(cols = c(Spearman_rho, Partial_rho),
               names_to = "Type", values_to = "rho") %>%
  mutate(Type = recode(Type,
                       Spearman_rho = "Spearman ρ",
                       Partial_rho  = "Partial ρ (age, sex)"))

p_corr <- ggplot(res_long, aes(x = Variable, y = rho, fill = Type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  coord_flip() +
  scale_fill_manual(values = c("Spearman ρ"="#2C7BB6", "Partial ρ (age, sex)"="#D7191C")) +
  geom_hline(yintercept = 0, color = "gray40", linetype = "dashed") +
  labs(y = "Correlation coefficient (ρ)", x = NULL,
       title = "Correlations between biochemical markers and disability (EDSS)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(out_dir, "Figure2E-F_Correlations.png"),
       p_corr, dpi = 600, width = 10, height = 7)

message("\n✅ Subsection 3.3 outputs saved to: ", out_dir,
        "\n - Table4_Spearman_EDSS_Q1.docx",
        "\n - Figure2E-F_Correlations.png")
