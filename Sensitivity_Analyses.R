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
# Section 3.5 – Sensitivity and Interaction Analyses
# ============================================================

# --- Libraries ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)
  library(rstatix)
  library(performance)
  library(ggplot2)
  library(ggpubr)
})

# --- Load data ---
data <- read.csv("Libya_MS.csv")

# --- Recreate baseline model (from Section 3.4) ---
fit_main <- MASS::polr(
  factor(EDSS_cat, ordered = TRUE) ~ age + sex + duration_of_illness_years +
    bmi_kg_m2 + vitamin_d + hdl + ldl_hdl_ratio + aip,
  data = data,
  Hess = TRUE
)

# ============================================================
# 3.5.1 Sensitivity Analyses
# ============================================================

# --- A. Model excluding outliers in biochemical variables ---
data_no_outliers <- data %>%
  filter(between(vitamin_d, quantile(vitamin_d, 0.01), quantile(vitamin_d, 0.99))) %>%
  filter(between(ldl_hdl_ratio, quantile(ldl_hdl_ratio, 0.01), quantile(ldl_hdl_ratio, 0.99)))

fit_sens1 <- MASS::polr(
  factor(EDSS_cat, ordered = TRUE) ~ age + sex + duration_of_illness_years +
    bmi_kg_m2 + vitamin_d + hdl + ldl_hdl_ratio + aip,
  data = data_no_outliers,
  Hess = TRUE
)

# --- B. Model restricted to RRMS phenotype only ---
fit_sens2 <- MASS::polr(
  factor(EDSS_cat, ordered = TRUE) ~ age + sex + duration_of_illness_years +
    bmi_kg_m2 + vitamin_d + hdl + ldl_hdl_ratio + aip,
  data = filter(data, ms_phenotype == "RRMS"),
  Hess = TRUE
)

# --- Compare ORs across sensitivity models ---
extract_or <- function(model) {
  broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    select(term, estimate, conf.low, conf.high) %>%
    mutate(OR = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high)) %>%
    select(term, OR)
}

sensitivity_summary <- bind_rows(
  extract_or(fit_main) %>% mutate(model = "Main model"),
  extract_or(fit_sens1) %>% mutate(model = "Excl. outliers"),
  extract_or(fit_sens2) %>% mutate(model = "RRMS only")
)

# Save summary table
write.csv(sensitivity_summary, "outputs/subsection3_5/Table7_Sensitivity_ORs.csv", row.names = FALSE)

# ============================================================
# 3.5.2 Interaction Analyses
# ============================================================

# --- Define interactions of interest ---
# Biological rationale: vitamin D may modify lipid–disability associations
# and BMI may modify vitamin D effect.

fit_int1 <- MASS::polr(
  factor(EDSS_cat, ordered = TRUE) ~ age + sex + duration_of_illness_years +
    bmi_kg_m2 * vitamin_d + hdl + ldl_hdl_ratio + aip,
  data = data,
  Hess = TRUE
)

fit_int2 <- MASS::polr(
  factor(EDSS_cat, ordered = TRUE) ~ age + sex + duration_of_illness_years +
    vitamin_d * hdl + bmi_kg_m2 + ldl_hdl_ratio + aip,
  data = data,
  Hess = TRUE
)

# --- Extract interaction coefficients ---
int_summary <- bind_rows(
  broom::tidy(fit_int1, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "BMI × Vit D"),
  broom::tidy(fit_int2, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "Vit D × HDL")
) %>%
  filter(str_detect(term, ":")) %>%
  select(model, term, estimate, conf.low, conf.high, p.value) %>%
  mutate(OR = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high))

write.csv(int_summary, "outputs/subsection3_5/Table8_Interactions.csv", row.names = FALSE)

# --- Visualization of interactions ---
ggplot(data, aes(x = vitamin_d, y = as.numeric(EDSS_cat), color = cut(bmi_kg_m2, 3))) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  labs(
    title = "Interaction between BMI and Vitamin D on disability (EDSS)",
    x = "Serum 25(OH) Vitamin D (ng/mL)",
    y = "EDSS category (numeric)",
    color = "BMI tertile"
  ) +
  theme_pubr(base_size = 14) +
  ggsave("outputs/subsection3_5/Figure4_Interaction_BMI_VitD.png", dpi = 600, width = 8, height = 6)
