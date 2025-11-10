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
# Figure 1 (A–D): Clinical & biochemical gradients across EDSS
# - A: Age & Age at Onset by EDSS
# - B: Disease Duration & ARR by EDSS
# - C: Vitamin D by EDSS (with 10 and 30 ng/mL cut lines)
# - D: Lipid ratios (TC/HDL, LDL/HDL, AIP) by EDSS
# Exports: PNG + TIFF at 600 DPI
# ================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
  library(ggplot2)
  library(patchwork)   # for multi-panel layout
  library(yaml)
})

# ---------- Config (edit path if needed) ----------
# If you already have config.yml, this will read column names from there.
cfg_path <- "analysis/00_setup/config.yml"
if (file.exists(cfg_path)) {
  cfg <- yaml::read_yaml(cfg_path)
  data_file <- cfg$data_file
  cn <- cfg$columns
} else {
  # Fallback if no config file present
  data_file <- "data/Libya_MS.csv"
  cn <- list(
    edss = "edss",
    age  = "age",
    age_onset = "age_of_onset_of_disease_years",
    duration_years = "duration_of_illness_years",
    arr = "arr",
    vitamin_d = "vitamin_d",
    chol_hdl_ratio = "chol_hdl_ratio",
    ldl_hdl_ratio  = "ldl_hdl_ratio",
    aip = "aip"
  )
}

# ---------- Load & clean data ----------
df <- readr::read_csv(data_file, show_col_types = FALSE) %>% clean_names()
names(df) <- gsub("\\.", "_", names(df))

# Helper to safely get a column by config alias
col <- function(name) {
  alias <- cn[[name]]
  if (!is.null(alias) && alias %in% names(df)) df[[alias]] else NA
}

# Ensure EDSS category exists
if (!"edss_cat" %in% names(df)) {
  edss_num <- suppressWarnings(as.numeric(col("edss")))
  df <- df %>%
    mutate(
      edss_num = edss_num,
      edss_cat = case_when(
        is.na(edss_num) ~ NA_character_,
        edss_num <= 2.5 ~ "Mild (0–2.5)",
        edss_num <= 5.5 ~ "Moderate (3–5.5)",
        edss_num >  5.5 ~ "Severe (≥6)"
      )
    )
} else {
  names(df)[names(df) == "edss_cat"] <- "edss_cat"
}

# Make factors ordered for consistent plotting
df <- df %>%
  mutate(
    edss_cat = factor(
      edss_cat,
      levels = c("Mild (0–2.5)", "Moderate (3–5.5)", "Severe (≥6)")
    )
  )

# Numeric getters
num <- function(x) suppressWarnings(as.numeric(x))

age         <- num(col("age"))
age_onset   <- num(col("age_onset"))
duration    <- num(col("duration_years"))
arr         <- num(col("arr"))
vitd        <- num(col("vitamin_d"))
chr         <- num(col("chol_hdl_ratio"))
lhr         <- num(col("ldl_hdl_ratio"))
aip         <- num(col("aip"))

df <- df %>%
  mutate(
    age = age, age_onset = age_onset, duration_years = duration, arr = arr,
    vitamin_d = vitd, chol_hdl_ratio = chr, ldl_hdl_ratio = lhr, aip = aip
  )

# ---------- Plotting helpers ----------
p_theme <- theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13, margin = margin(b = 6)),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

# Small jitter for points
point_jitter <- position_jitter(width = 0.12, height = 0)

# ---------- Panel A: Age & Age at Onset ----------
dfA <- df %>%
  select(edss_cat, age, age_onset) %>%
  pivot_longer(cols = c(age, age_onset),
               names_to = "measure", values_to = "value") %>%
  mutate(
    measure = factor(measure,
                     levels = c("age", "age_onset"),
                     labels = c("Age (years)", "Age at onset (years)"))
  )

pA <- ggplot(dfA, aes(x = edss_cat, y = value)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(position = point_jitter, alpha = 0.25, size = 1) +
  stat_summary(fun = median, geom = "crossbar", width = 0.45) +
  facet_wrap(~ measure, nrow = 1, scales = "free_y") +
  labs(title = "A) Age and age at onset by EDSS category",
       x = "EDSS category", y = "Years") +
  p_theme

# ---------- Panel B: Disease Duration & ARR ----------
dfB <- df %>%
  select(edss_cat, duration_years, arr) %>%
  pivot_longer(cols = c(duration_years, arr),
               names_to = "measure", values_to = "value") %>%
  mutate(
    measure = factor(measure,
                     levels = c("duration_years", "arr"),
                     labels = c("Disease duration (years)", "Annualized relapse rate (ARR)"))
  )

pB <- ggplot(dfB, aes(x = edss_cat, y = value)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(position = point_jitter, alpha = 0.25, size = 1) +
  stat_summary(fun = median, geom = "crossbar", width = 0.45) +
  facet_wrap(~ measure, nrow = 1, scales = "free_y") +
  labs(title = "B) Disease duration increases while ARR declines",
       x = "EDSS category", y = NULL) +
  p_theme

# ---------- Panel C: Vitamin D with cutoffs ----------
pC <- ggplot(df, aes(x = edss_cat, y = vitamin_d)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(position = point_jitter, alpha = 0.25, size = 1) +
  stat_summary(fun = median, geom = "crossbar", width = 0.45) +
  geom_hline(yintercept = 10, linetype = "dashed") +
  geom_hline(yintercept = 30, linetype = "dashed") +
  annotate("text", x = 3.05, y = 10, label = "10 ng/mL", hjust = 1, vjust = -0.4, size = 3) +
  annotate("text", x = 3.05, y = 30, label = "30 ng/mL", hjust = 1, vjust = -0.4, size = 3) +
  labs(title = "C) Serum 25(OH) vitamin D by EDSS category",
       x = "EDSS category", y = "Vitamin D (ng/mL)") +
  p_theme

# ---------- Panel D: Lipid ratios ----------
dfD <- df %>%
  select(edss_cat, chol_hdl_ratio, ldl_hdl_ratio, aip) %>%
  pivot_longer(cols = c(chol_hdl_ratio, ldl_hdl_ratio, aip),
               names_to = "ratio", values_to = "value") %>%
  mutate(
    ratio = factor(ratio,
                   levels = c("chol_hdl_ratio", "ldl_hdl_ratio", "aip"),
                   labels = c("TC / HDL", "LDL / HDL", "AIP"))
  )

pD <- ggplot(dfD, aes(x = edss_cat, y = value)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(position = point_jitter, alpha = 0.25, size = 1) +
  stat_summary(fun = median, geom = "crossbar", width = 0.45) +
  facet_wrap(~ ratio, nrow = 1, scales = "free_y") +
  labs(title = "D) Lipid ratios rise with increasing disability",
       x = "EDSS category", y = NULL) +
  p_theme

# ---------- Assemble 2x2 composite ----------
# Top row: A + B; Bottom row: C + D
composite <- (pA | pB) / (pC | pD)
composite <- composite + plot_annotation(
  theme = theme(plot.margin = margin(t = 6, r = 6, b = 6, l = 6))
)

# ---------- Export ----------
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("outputs/figures/Figure1_MS_Baseline.png", composite,
       width = 10, height = 8, units = "in", dpi = 600)
ggsave("outputs/figures/Figure1_MS_Baseline.tiff", composite,
       width = 10, height = 8, units = "in", dpi = 600, compression = "lzw")

message("✅ Saved Figure1_MS_Baseline.{png,tiff} at 600 DPI in outputs/figures")
