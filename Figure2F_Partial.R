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

# ================= Figure 2F: Partial correlations (Age & Sex adjusted) =================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(janitor)
  library(ppcor)
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
edss_var <- "edss"; sex_var <- "gender"

stopifnot(edss_var %in% names(df))
df2 <- df %>%
  mutate(
    across(all_of(intersect(vars, names(df))), to_num),
    edss = to_num(.data[[edss_var]]),
    Age  = if ("age" %in% names(.)) to_num(age) else NA_real_,
    Sex  = case_when(
      tolower(as.character(.data[[sex_var]])) %in% c("male","m","1") ~ 1,
      tolower(as.character(.data[[sex_var]])) %in% c("female","f","0") ~ 0,
      TRUE ~ NA_real_
    )
  ) %>% filter(is.finite(edss))

vars <- vars[vars %in% names(df2)]

safe_partial <- function(x, y, age, sex){
  ok <- is.finite(x) & is.finite(y) & is.finite(age) & is.finite(sex)
  x <- x[ok]; y <- y[ok]; age <- age[ok]; sex <- sex[ok]
  if (length(x) < 10 || length(unique(x))<3 || length(unique(y))<3 ||
      length(unique(age))<3 || length(unique(sex))<2) return(NA_real_)
  out <- try(ppcor::pcor.test(x, y, cbind(age, sex), method = "spearman"), silent = TRUE)
  if (inherits(out, "try-error")) return(NA_real_)
  unname(out$estimate)
}

rho_p <- sapply(vars, function(v) safe_partial(df2[[v]], df2$edss, df2$Age, df2$Sex))

lab <- c(
  age="Age", age_of_onset_of_disease_years="Age at onset",
  duration_of_illness_years="Disease duration", arr="ARR",
  bmi_kg_m2="BMI", vitamin_d="Vitamin D",
  total_cholesterol="Total cholesterol", ldl="LDL", hdl="HDL", tg="Triglycerides",
  chol_hdl_ratio="Total/HDL ratio", ldl_hdl_ratio="LDL/HDL ratio",
  aip="Atherogenic Index (AIP)"
)

plot_df <- data.frame(Variable = lab[names(rho_p)], rho = as.numeric(rho_p), stringsAsFactors = FALSE)

# To keep identical order to panel E, sort by absolute *unadjusted* effect if you saved it:
# If you don't have panel E data handy, we sort by abs(partial rho), but you can pass the
# order manually by reading Figure2E table to match exactly.
plot_df <- plot_df |> arrange(desc(abs(rho))) |> # replace 'rho' with your unadjusted vector to lock ordering
  mutate(Variable = factor(Variable, levels = rev(Variable)))

# Match the same x-limits you used for panel E (set manually if needed)
xmax <- max(0.35, ceiling(max(abs(plot_df$rho))*10)/10)
lims <- c(-xmax, xmax)

pF <- ggplot(plot_df, aes(x = Variable, y = rho)) +
  geom_col(fill = "#E15759", width = 0.65) +
  coord_flip(ylim = lims) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray55") +
  labs(
    title = "F) Partial correlations (adjusted for Age and Sex)",
    subtitle = "Biochemical and clinical markers vs EDSS",
    x = NULL, y = "Correlation coefficient (ρ)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "Figure2F_Correlations_Partial.png"),
       pF, dpi = 600, width = 8.5, height = 6.5)
ggsave(file.path(out_dir, "Figure2F_Correlations_Partial.tiff"),
       pF, dpi = 600, width = 8.5, height = 6.5, compression = "lzw")
