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

# analysis/02_models/20_models.R
library(tidyverse); library(readr); library(broom); library(broom.helpers)
library(effectsize); library(MASS); library(ordinal); library(rsq); library(sandwich); library(lmtest)
source("analysis/00_setup/functions.R")
cfg <- yaml::read_yaml("analysis/00_setup/config.yml")

df <- readr::read_csv(cfg$data_file) %>% janitor::clean_names()
names(df) <- gsub("\\.", "_", names(df))
df <- derive_vars(df) %>% clean_levels()

# Outcome as ordered factor
df$EDSS_cat <- factor(df$EDSS_cat, levels=c("Mild (0–2.5)","Moderate (3–5.5)","Severe (≥6.0)"), ordered = TRUE)

# Model formula
form <- as.formula(paste0("EDSS_cat ~ ", paste(na.omit(c(cfg$columns$vitd, cfg$columns$crp, cfg$columns$hba1c,
                                            cfg$columns$tc, cfg$columns$hdl, cfg$columns$tg,
                                            "Chol_HDL_Ratio","LDL_HDL_Ratio","AIP",
                                            cfg$columns$bmi, cfg$columns$age, cfg$columns$sex,
                                            cfg$columns$phenotype, cfg$columns$duration_years)), collapse = " + ")))

# Ordinal regression
m_polr <- MASS::polr(form, data=df, Hess=TRUE, na.action = na.omit)
co <- broom::tidy(m_polr, conf.int = TRUE, exponentiate = TRUE)
std <- effectsize::standardize_parameters(m_polr)
out <- dplyr::left_join(co, std, by="term")
out <- out %>% rename(OR = estimate, CI_low = conf.low, CI_high = conf.high, beta_std = Std_Coefficient)

# Partial R2 (McFadden-like)
r2 <- rsq::rsq(m_polr, adj=TRUE, type="v")

# Multiplicity correction
out$p_adj <- p.adjust(out$p.value, method = "BH")

readr::write_csv(out, "outputs/tables/Table3_Ordinal_EDSS_predictors.csv")
readr::write_lines(paste0("Adjusted pseudo-R2: ", round(r2, 3)), "outputs/tables/Table3_Ordinal_EDSS_R2.txt")

# Linear EDSS (sensitivity)
f_lin <- as.formula(paste0(cfg$columns$edss, " ~ ", paste(na.omit(c(cfg$columns$vitd, cfg$columns$crp, cfg$columns$hba1c,
                                            cfg$columns$tc, cfg$columns$hdl, cfg$columns$tg,
                                            "Chol_HDL_Ratio","LDL_HDL_Ratio","AIP",
                                            cfg$columns$bmi, cfg$columns$age, cfg$columns$sex,
                                            cfg$columns$phenotype, cfg$columns$duration_years)), collapse = " + ")))
m_lin <- lm(f_lin, data=df)
lin_co <- broom::tidy(m_lin, conf.int = TRUE)
lin_std <- effectsize::standardize_parameters(m_lin)
lin_out <- dplyr::left_join(lin_co, lin_std, by="term") %>%
  rename(beta = estimate, CI_low = conf.low, CI_high = conf.high, beta_std = Std_Coefficient)
lin_out$p_adj <- p.adjust(lin_out$p.value, method = "BH")
readr::write_csv(lin_out, "outputs/tables/TableS_Linear_EDSS.csv")
