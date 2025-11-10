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
# =====================================================================
# Figure 1 (A–D) — Colored version + individual panels (A,B,C,D)
# Clinical & biochemical gradients across EDSS categories
# Exports: PNG + TIFF at 600 DPI (composite + separate panels)
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)  # colorblind-safe palette (via scale_*_viridis_d)
  library(yaml)
})

# ------------------------ Config / Columns ----------------------------
cfg_path <- "analysis/00_setup/config.yml"
if (file.exists(cfg_path)) {
  cfg <- yaml::read_yaml(cfg_path)
  data_file <- cfg$data_file
  cn <- cfg$columns
} else {
  # Fallback if config absent; adjust if needed
  data_file <- "data/Libya_MS.csv"
  cn <- list(
    edss = "edss",
    age  = "age",
    age_onset = "age_of_onset_of_disease_years",
    duration_years = "duration_of_illness_years",
    arr = "arr",
    vitd = "vitamin_d",
    chol_hdl_ratio = "chol_hdl_ratio",
    ldl_hdl_ratio  = "ldl_hdl_ratio",
    aip = "aip"
  )
}

# ------------------------ Load & clean data ---------------------------
df <- readr::read_csv(data_file, show_col_types = FALSE) %>% clean_names()
names(df) <- gsub("\\.", "_", names(df))

# Helper to fetch column by alias name from config
col_ <- function(key) {
  alias <- cn[[key]]
  if (!is.null(alias) && alias %in% names(df)) df[[alias]] else NULL
}

# Ensure EDSS category
edss_num <- suppressWarnings(as.numeric(col_("edss") %||% df$edss))
df <- df %>%
  mutate(
    edss_num = edss_num,
    edss_cat = case_when(
      is.na(edss_num) ~ NA_character_,
      edss_num <= 2.5 ~ "Mild (0–2.5)",
      edss_num <= 5.5 ~ "Moderate (3–5.5)",
      edss_num >  5.5 ~ "Severe (≥6)"
    ),
    edss_cat = factor(edss_cat, levels = c("Mild (0–2.5)", "Moderate (3–5.5)", "Severe (≥6)"))
  )

# Numeric getters
num <- function(x) suppressWarnings(as.numeric(x))
df <- df %>%
  mutate(
    age            = num(col_("age")        %||% df$age),
    age_onset      = num(col_("age_onset")  %||% df$age_of_onset_of_disease_years),
    duration_years = num(col_("duration_years") %||% df$duration_of_illness_years),
    arr            = num(col_("arr")        %||% df$arr),
    chol_hdl_ratio = num(col_("chol_hdl_ratio") %||% df$chol_hdl_ratio),
    ldl_hdl_ratio  = num(col_("ldl_hdl_ratio")  %||% df$ldl_hdl_ratio),
    aip            = num(col_("aip")            %||% df$aip)
  )

# Robust vitamin D detection (handles several name variants)
vitd_candidates <- c(
  cn$vitd, cn$vitamin_d,
  "vitamin_d", "vit_d", "vitamin_d_ng_ml", "vitamin_d_ng_ml_", "vitd"
)
vitd_name <- vitd_candidates[vitd_candidates %in% names(df)][1]
if (!is.na(vitd_name) && length(vitd_name)) {
  df$vitamin_d <- num(df[[vitd_name]])
} else {
  df$vitamin_d <- NA_real_
  message("⚠️ No vitamin D column found among: ", paste(vitd_candidates, collapse = ", "))
}

# ------------------------ Plotting helpers ----------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

edss_labels <- c("Mild (0–2.5)"="Mild", "Moderate (3–5.5)"="Moderate", "Severe (≥6)"="Severe")

# Color palette per EDSS group
edss_colors <- setNames(viridisLite::viridis(3, option = "D"), levels(df$edss_cat))

p_theme <- theme_bw(base_size = 12) +
  theme(
    legend.position   = "none",
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey90", colour = "grey60"),
    strip.text        = element_text(face = "bold", size = 12, margin = margin(b = 4)),
    plot.title        = element_text(face = "bold", size = 14, margin = margin(b = 8)),
    axis.title        = element_text(size = 12),
    axis.text.x       = element_text(size = 11, margin = margin(t = 3))
  )

point_jitter <- position_jitter(width = 0.12, height = 0)

annotate_n <- function(d, xvar, y_pos) {
  d %>% filter(!is.na(.data[[xvar]])) %>%
    count(.data[[xvar]]) %>%
    mutate(label = paste0("", NULL), y = 0)
}

save_both <- function(p, path_base, w = 10, h = 8, dpi = 600) {
  ggsave(paste0(path_base, ".png"),  p, width = w, height = h, units = "in", dpi = dpi)
  ggsave(paste0(path_base, ".tiff"), p, width = w, height = h, units = "in", dpi = dpi, compression = "lzw")
}

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------ Panel A: Age & Age at Onset ------------------
dfA <- df %>%
  select(edss_cat, age, age_onset) %>%
  pivot_longer(c(age, age_onset), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = c("age","age_onset"),
                          labels = c("Age (years)", "Age at onset (years)")))

nA <- annotate_n(dfA %>% filter(is.finite(value)), "edss_cat", y_pos = max(dfA$value, na.rm = TRUE)*0.02)

pA <- ggplot(dfA, aes(x = edss_cat, y = value, fill = edss_cat)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85, colour = "grey20") +
  geom_jitter(position = point_jitter, alpha = 0.30, size = 1, colour = "grey20") +
  stat_summary(fun = median, geom = "crossbar", width = 0.45, colour = "black") +
  geom_text(data = nA, aes(x = edss_cat, y = y, label = label), inherit.aes = FALSE, size = 3.2, vjust = 1.1) +
  facet_wrap(~ measure, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = edss_colors) +
  scale_x_discrete(labels = edss_labels) +
  labs(title = "A) Age and age at onset by EDSS category", x = "EDSS category", y = "Years") +
  p_theme

# ------------------------ Panel B: Duration & ARR ----------------------
dfB <- df %>%
  select(edss_cat, duration_years, arr) %>%
  pivot_longer(c(duration_years, arr), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure,
                          levels = c("duration_years","arr"),
                          labels = c("Disease duration (years)", "ARR")))

nB <- annotate_n(dfB %>% filter(is.finite(value)), "edss_cat", y_pos = max(dfB$value, na.rm = TRUE)*0.02)

pB <- ggplot(dfB, aes(x = edss_cat, y = value, fill = edss_cat)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85, colour = "grey20") +
  geom_jitter(position = point_jitter, alpha = 0.30, size = 1, colour = "grey20") +
  stat_summary(fun = median, geom = "crossbar", width = 0.45, colour = "black") +
  geom_text(data = nB, aes(x = edss_cat, y = y, label = label), inherit.aes = FALSE, size = 3.2, vjust = 1.1) +
  facet_wrap(~ measure, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = edss_colors) +
  scale_x_discrete(labels = edss_labels) +
  labs(title = "B) Disease duration increases while ARR declines", x = "EDSS category", y = NULL) +
  p_theme

# ------------------------ Panel C: Vitamin D ---------------------------
has_vtd <- any(is.finite(df$vitamin_d))
vtd_top <- if (has_vtd) max(df$vitamin_d, na.rm = TRUE) else 35
vtd_nlab <- annotate_n(df %>% filter(is.finite(vitamin_d)), "edss_cat", y_pos = vtd_top*0.04)

pC <- ggplot(df, aes(x = edss_cat, y = vitamin_d, fill = edss_cat)) +
  { if (has_vtd) geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85, colour = "grey20") else geom_blank() } +
  { if (has_vtd) geom_jitter(position = point_jitter, alpha = 0.30, size = 1, colour = "grey20") else geom_blank() } +
  { if (has_vtd) stat_summary(fun = median, geom = "crossbar", width = 0.45, colour = "black") else geom_blank() } +
  geom_hline(yintercept = 10, linetype = "dashed") +
  geom_hline(yintercept = 30, linetype = "dashed") +
  annotate("text", x = 3.05, y = 10, label = "", hjust = 1, vjust = -0.4, size = 3) +
  annotate("text", x = 3.05, y = 30, label = "", hjust = 1, vjust = -0.4, size = 3) +
  { if (has_vtd && nrow(vtd_nlab)) geom_text(data = vtd_nlab, aes(x = edss_cat, y = y, label = label),
                                             inherit.aes = FALSE, size = 3.2, vjust = 1.1) else geom_blank() } +
  scale_fill_manual(values = edss_colors) +
  scale_x_discrete(labels = edss_labels) +
  coord_cartesian(ylim = c(0, max(35, vtd_top + 3))) +
  labs(title = "C) Serum 25(OH) vitamin D by EDSS category", x = "EDSS category", y = "Vitamin D (ng/mL)") +
  p_theme +
  { if (!has_vtd) annotate("text", x = 2, y = 18, label = "No vitamin D data detected", size = 4, fontface = "bold") else NULL }

# ------------------------ Panel D: Lipid ratios ------------------------
dfD <- df %>%
  select(edss_cat, chol_hdl_ratio, ldl_hdl_ratio, aip) %>%
  pivot_longer(c(chol_hdl_ratio, ldl_hdl_ratio, aip),
               names_to = "ratio", values_to = "value") %>%
  mutate(ratio = factor(ratio,
                        levels = c("chol_hdl_ratio","ldl_hdl_ratio","aip"),
                        labels = c("TC / HDL", "LDL / HDL", "AIP")))

nD <- dfD %>%
  filter(is.finite(value)) %>%
  group_by(ratio) %>%
  group_map(~annotate_n(.x, "edss_cat", y_pos = max(.x$value, na.rm = TRUE)*0.04)) %>%
  list_rbind()

pD <- ggplot(dfD, aes(x = edss_cat, y = value, fill = edss_cat)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85, colour = "grey20") +
  geom_jitter(position = point_jitter, alpha = 0.30, size = 1, colour = "grey20") +
  stat_summary(fun = median, geom = "crossbar", width = 0.45, colour = "black") +
  geom_text(data = nD, aes(x = edss_cat, y = y, label = label), inherit.aes = FALSE, size = 3.0, vjust = 1.1) +
  facet_wrap(~ ratio, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = edss_colors) +
  scale_x_discrete(labels = edss_labels) +
  labs(title = "D) Lipid ratios rise with increasing disability", x = "EDSS category", y = NULL) +
  p_theme

# ------------------------ Composite (2x2) ------------------------------
composite <- (pA | pB) / (pC | pD)
composite <- composite + plot_annotation(
  theme = theme(plot.margin = margin(t = 6, r = 6, b = 6, l = 6))
)

# ------------------------ Save: composite + panels ---------------------
save_both(composite, "outputs/figures/Figure1_MS_Baseline_colored", w = 10, h = 8, dpi = 600)
save_both(pA,        "outputs/figures/Figure1A_Age_Onset_colored",  w = 8,  h = 4.5, dpi = 600)
save_both(pB,        "outputs/figures/Figure1B_Duration_ARR_colored", w = 8, h = 4.5, dpi = 600)
save_both(pC,        "outputs/figures/Figure1C_VitD_colored",        w = 8,  h = 4.5, dpi = 600)
save_both(pD,        "outputs/figures/Figure1D_LipidRatios_colored", w = 8,  h = 4.5, dpi = 600)

message("✅ Saved: composite + individual panels (PNG & TIFF, 600 DPI) in outputs/figures/")
