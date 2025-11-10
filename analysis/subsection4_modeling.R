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

library(clinfun)   # for jonckheere.test
library(dplyr)

# Ensure EDSS group is ordered
df <- df %>%
  mutate(edss_group = factor(edss_group,
                             levels = c("Mild (0–2.5)", "Moderate (3–5.5)", "Severe (≥6)"),
                             ordered = TRUE
  ))

vars <- c("age", "age_of_onset_of_disease_years", "duration_of_illness_years",
          "arr", "bmi_kg_m2", "vitamin_d", "total_cholesterol", "ldl", "hdl",
          "tg", "chol_hdl_ratio", "ldl_hdl_ratio", "aip")

results <- lapply(vars, function(v) {
  sub <- df[, c(v, "edss_group")] %>% drop_na()
  kw_p <- kruskal.test(as.formula(paste(v, "~ edss_group")), data = sub)$p.value
  jt_p <- jonckheere.test(sub[[v]], as.numeric(sub$edss_group))$p.value
  data.frame(
    Variable = v,
    n = nrow(sub),
    `Kruskal–Wallis p` = sprintf("%.4f", kw_p),
    `Trend p (JT)` = sprintf("%.4f", jt_p),
    Direction = ifelse(mean(sub[[v]][sub$edss_group == "Severe (≥6)"]) >
                         mean(sub[[v]][sub$edss_group == "Mild (0–2.5)"]),
                       "↑ with EDSS", "↓ with EDSS")
  )
}) %>% bind_rows()
