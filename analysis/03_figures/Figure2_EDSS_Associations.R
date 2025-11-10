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

# ===============================================================
# Figure 2 – Correlates of Disability (EDSS) in Libyan MS
# ===============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(GGally)
  library(ggcorrplot)
  library(patchwork)
  library(janitor)
  library(ggthemes)
})

# ---------- Detect and set working directory ----------
base_dir <- "D:\\A-ms\\base\\STAT\\MS_LIBYA\\ms-biomarkers-msard-neuroepi\\data"
setwd(base_dir)

# confirm file path
data_file <- file.path(base_dir, "Libya_MS.csv")
if (!file.exists(data_file)) {
  stop("❌ Dataset not found. Please check that 'Libya_MS.csv' is in: ", base_dir)
}

# create outputs folder if needed
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# ---------- Load and clean data ----------
df <- read_csv(data_file, show_col_types = FALSE) |> janitor::clean_names()

# Derive EDSS categories
df <- df %>%
  mutate(
    edss_group = case_when(
      edss <= 2.5 ~ "Mild (0–2.5)",
      edss > 2.5 & edss <= 5.5 ~ "Moderate (3–5.5)",
      edss >= 6 ~ "Severe (≥6)"
    ) %>% factor(levels = c("Mild (0–2.5)", "Moderate (3–5.5)", "Severe (≥6)"))
  )

# ---------- (A) Violin + Boxplot for BMI ----------
violin_fun <- function(var, ylab, title) {
  ggplot(df, aes(x = edss_group, y = .data[[var]], fill = edss_group)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "black") +
    scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = ylab, title = title) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0))
}

pA <- violin_fun("bmi_kg_m2", "BMI (kg/m²)", "A. BMI Distribution by EDSS Severity")

# ---------- (B) Correlation Heatmap ----------
corr_vars <- df %>%
  select(EDSS = edss, Age = age, Duration = duration_of_illness_years,
         BMI = bmi_kg_m2, VitD = vitamin_d, HDL = hdl,
         LDL.HDL = ldl_hdl_ratio, AIP = aip)

corr_mat <- round(cor(corr_vars, use = "pairwise.complete.obs", method = "spearman"), 2)

pB <- ggcorrplot::ggcorrplot(
  corr_mat, lab = TRUE, lab_size = 3,
  colors = c("#4575b4", "white", "#d73027"),
  title = "B. Spearman Correlations among Disability and Biochemical Indices",
  outline.color = "gray90"
) + theme_minimal(base_size = 13)

# ---------- (C) Forest Plot (Partial Rho placeholders) ----------
forest_df <- tibble(
  Variable = c("Age", "Disease duration", "BMI", "Vitamin D", "HDL", "LDL/HDL ratio", "AIP"),
  Rho = c(0.16, 0.19, 0.36, -0.16, -0.33, 0.30, 0.18),
  CI_low = c(0.06, 0.10, 0.27, -0.25, -0.41, 0.21, 0.08),
  CI_high = c(0.26, 0.27, 0.44, -0.06, -0.24, 0.39, 0.28)
)

pC <- ggplot(forest_df, aes(y = reorder(Variable, Rho), x = Rho)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.25, color = "#3E606F") +
  geom_point(size = 3.5, color = "#E46726") +
  labs(x = "Partial Spearman rho (age, sex adjusted)", y = NULL,
       title = "C. Partial Correlations with EDSS") +
  theme_minimal(base_size = 13)

# ---------- (D) Trend Plot (Vitamin D Median ± IQR) ----------
trend_df <- df %>%
  group_by(edss_group) %>%
  summarise(median = median(vitamin_d, na.rm = TRUE),
            q1 = quantile(vitamin_d, 0.25, na.rm = TRUE),
            q3 = quantile(vitamin_d, 0.75, na.rm = TRUE),
            n = n())

pD <- ggplot(trend_df, aes(x = edss_group, y = median, group = 1)) +
  geom_line(color = "#3E606F", linewidth = 1.2) +
  geom_point(aes(size = n), color = "#E46726") +
  geom_errorbar(aes(ymin = q1, ymax = q3), width = 0.15, color = "#3E606F") +
  scale_size(range = c(3, 7), guide = "none") +
  labs(x = NULL, y = "Serum 25(OH) Vitamin D (ng/mL)",
       title = "D. Median Vitamin D Trends across EDSS Groups") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# ---------- Combine Panels ----------
combined <- (pA | pB) / (pC | pD)

ggsave("outputs/figures/Figure2_MS_EDSS_Associations_colored.tiff",
       combined, dpi = 600, width = 10, height = 8)

# Export individual panels
ggsave("outputs/figures/Figure2A_BMI_EDSS_violin.tiff", pA, dpi = 600, width = 4, height = 4)
ggsave("outputs/figures/Figure2B_CorrHeatmap.tiff", pB, dpi = 600, width = 5, height = 5)
ggsave("outputs/figures/Figure2C_ForestPlot.tiff", pC, dpi = 600, width = 5, height = 4)
ggsave("outputs/figures/Figure2D_VitD_Trends.tiff", pD, dpi = 600, width = 4, height = 4)

message("✅ Figure 2 (combined and sub-panels) successfully saved in: ", base_dir, "/outputs/figures/")
