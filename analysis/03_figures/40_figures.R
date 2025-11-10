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

# analysis/03_figures/40_figures.R
library(tidyverse); library(readr); library(ggplot2); library(patchwork); library(broom)
source("analysis/00_setup/functions.R")
cfg <- yaml::read_yaml("analysis/00_setup/config.yml")
df <- readr::read_csv(cfg$data_file) %>% janitor::clean_names()
names(df) <- gsub("\\.", "_", names(df))
df <- derive_vars(df) %>% clean_levels()

# Correlation heatmap for biochemical markers
bio <- c(cfg$columns$vitd, cfg$columns$crp, cfg$columns$hba1c, cfg$columns$tc, cfg$columns$ldl, cfg$columns$hdl, cfg$columns$tg, "Chol_HDL_Ratio","LDL_HDL_Ratio","AIP")
bio <- bio[bio %in% names(df)]
M <- cor(df[, bio], use="pairwise.complete.obs", method="spearman")
Md <- reshape2::melt(M, varnames=c("x","y"), value.name="rho")
p1 <- ggplot(Md, aes(x, y, fill=rho)) +
  geom_tile() + scale_fill_gradient2(limits=c(-1,1)) +
  geom_text(aes(label=sprintf("%.2f", rho)), size=3) +
  coord_equal() + theme_minimal() + labs(title="Biochemical markers — Spearman correlation")
ggsave("outputs/figures/Fig2_CorrelationHeatmap.pdf", p1, width=7, height=6, dpi=600)
ggsave("outputs/figures/Fig2_CorrelationHeatmap.png", p1, width=7, height=6, dpi=600)

# Forest plot from ordinal model output
ord_file <- "outputs/tables/Table3_Ordinal_EDSS_predictors.csv"
if (file.exists(ord_file)) {
  ord <- readr::read_csv(ord_file) %>% filter(!is.na(OR))
  ord <- ord %>% mutate(term = forcats::fct_reorder(term, OR))
  p2 <- ggplot(ord, aes(x=term, y=OR)) +
    geom_point() +
    geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0.1) +
    coord_flip() + theme_minimal() +
    labs(title="Adjusted predictors of EDSS category (OR, 95% CI)", x="", y="Odds Ratio")
  ggsave("outputs/figures/Fig1_ForestPlot.pdf", p2, width=7, height=6, dpi=600)
  ggsave("outputs/figures/Fig1_ForestPlot.png", p2, width=7, height=6, dpi=600)
}
