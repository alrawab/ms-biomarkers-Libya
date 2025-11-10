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

# ================= Figure 2E: Unadjusted Spearman correlations =================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(janitor)
  library(ggplot2); library(scales)
})

# ---- Paths ----
base_dir <- getwd()
data_file <- file.path(base_dir, "Libya_MS.csv")
out_dir   <- file.path(base_dir, "outputs", "subsection3_3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load & clean ----
df <- readr::read_csv(data_file, show_col_types = FALSE) |> janitor::clean_names()
names(df) <- gsub("\\.", "_", names(df))

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

vars <- c("age","age_of_onset_of_disease_years","duration_of_illness_years","arr",
          "bmi_kg_m2","vitamin_d","total_cholesterol","ldl","hdl","tg",
          "chol_hdl_ratio","ldl_hdl_ratio","aip")
edss_var <- "edss"

stopifnot(edss_var %in% names(df))
df2 <- df |> mutate(across(all_of(intersect(vars, names(df))), to_num),
                    edss = to_num(.data[[edss_var]])) |>
  filter(is.finite(edss))

vars <- vars[vars %in% names(df2)]

safe_spearman <- function(x, y){
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(unique(x)) < 3 || length(unique(y)) < 3) return(NA_real_)
  unname(suppressWarnings(cor(x, y, method="spearman")))
}

rho <- sapply(vars, function(v) safe_spearman(df2[[v]], df2$edss))

lab <- c(
  age="Age", age_of_onset_of_disease_years="Age at onset",
  duration_of_illness_years="Disease duration", arr="ARR",
  bmi_kg_m2="BMI", vitamin_d="Vitamin D",
  total_cholesterol="Total cholesterol", ldl="LDL", hdl="HDL", tg="Triglycerides",
  chol_hdl_ratio="Total/HDL ratio", ldl_hdl_ratio="LDL/HDL ratio",
  aip="Atherogenic Index (AIP)"
)

plot_df <- data.frame(Variable = lab[names(rho)], rho = as.numeric(rho), stringsAsFactors = FALSE)
# Order by absolute effect (and keep same order for panel F later)
plot_df <- plot_df |> arrange(desc(abs(rho))) |>
  mutate(Variable = factor(Variable, levels = rev(Variable)))

# Fix common x-limits to make panel matching easy
xmax <- max(0.35, ceiling(max(abs(plot_df$rho))*10)/10)  # at least ±0.35
lims <- c(-xmax, xmax)

pE <- ggplot(plot_df, aes(x = Variable, y = rho)) +
  geom_col(fill = "#2F6CBA", width = 0.65) +
  coord_flip(ylim = lims) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray55") +
  labs(
    title = "E) Unadjusted Spearman correlations",
    subtitle = "Biochemical and clinical markers vs EDSS",
    x = NULL, y = "Correlation coefficient (ρ)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "Figure2E_Correlations_Unadjusted.png"),
       pE, dpi = 600, width = 8.5, height = 6.5)
ggsave(file.path(out_dir, "Figure2E_Correlations_Unadjusted.tiff"),
       pE, dpi = 600, width = 8.5, height = 6.5, compression = "lzw")
