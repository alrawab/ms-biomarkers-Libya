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
# Figure 1 (A–D) — Colored + individual panels, polished (v4)
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
  library(yaml)
  library(rlang)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# --------------------- Styling knobs (easy to tune) -------------------
BASE_SIZE   <- 12
TITLE_SIZE  <- 16
STRIP_SIZE  <- 13
POINT_ALPHA <- 0.25
POINT_SIZE  <- 0.9
MEDIAN_W    <- 0.9     # median crossbar linewidth
N_SIZE      <- 3.1
N_OFFSET    <- 0.075   # fraction below min per facet for n-labels
BOTTOM_PAD  <- 18      # bottom margin (px) to avoid clipping

# ------------------------ Config / Columns ----------------------------
cfg_path <- "analysis/00_setup/config.yml"
if (file.exists(cfg_path)) {
  cfg <- yaml::read_yaml(cfg_path)
  data_file <- cfg$data_file
  cn <- cfg$columns
} else {
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

col_ <- function(key) { alias <- cn[[key]]; if (!is.null(alias) && alias %in% names(df)) df[[alias]] else NULL }
num  <- function(x) suppressWarnings(as.numeric(x))

# EDSS category
edss_num <- num(col_("edss") %||% df$edss)
df <- df %>%
  mutate(
    edss_num = edss_num,
    edss_cat = case_when(
      is.na(edss_num) ~ NA_character_,
      edss_num <= 2.5 ~ "Mild (0–2.5)",
      edss_num <= 5.5 ~ "Moderate (3–5.5)",
      TRUE            ~ "Severe (≥6)"
    ),
    edss_cat = factor(edss_cat, levels = c("Mild (0–2.5)", "Moderate (3–5.5)", "Severe (≥6)")),
    age            = num(col_("age")            %||% df$age),
    age_onset      = num(col_("age_onset")      %||% df$age_of_onset_of_disease_years),
    duration_years = num(col_("duration_years") %||% df$duration_of_illness_years),
    arr            = num(col_("arr")            %||% df$arr),
    chol_hdl_ratio = num(col_("chol_hdl_ratio") %||% df$chol_hdl_ratio),
    ldl_hdl_ratio  = num(col_("ldl_hdl_ratio")  %||% df$ldl_hdl_ratio),
    aip            = num(col_("aip")            %||% df$aip)
  )

# Vitamin D mapping (robust to common variants)
vitd_candidates <- c(cn$vitd, cn$vitamin_d, "vitamin_d","vit_d","vitamin_d_ng_ml","vitd")
vitd_name <- vitd_candidates[vitd_candidates %in% names(df)][1]
df$vitamin_d <- if (!is.na(vitd_name) && length(vitd_name)) num(df[[vitd_name]]) else NA_real_

# -------------------------- Visual helpers ----------------------------
edss_labels <- c("Mild (0–2.5)"="Mild", "Moderate (3–5.5)"="Moderate", "Severe (≥6)"="Severe")
edss_colors <- setNames(viridisLite::viridis(3, option = "D"), levels(df$edss_cat))

p_theme <- theme_bw(base_size = BASE_SIZE) +
  theme(
    legend.position   = "none",
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey90", colour = "grey60"),
    strip.text        = element_text(face = "bold", size = STRIP_SIZE, margin = margin(b = 5)),
    plot.title        = element_text(face = "bold", size = TITLE_SIZE, margin = margin(b = 8)),
    axis.title        = element_text(size = BASE_SIZE + 1),
    axis.text.x       = element_text(size = BASE_SIZE + 1, margin = margin(t = 3))
  )

point_jitter <- position_jitter(width = 0.12, height = 0)
box_layer    <- list(geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.88, colour = "grey20"))
dots_layer   <- list(geom_jitter(position = point_jitter, alpha = POINT_ALPHA, size = POINT_SIZE, colour = "grey20"))
median_bar   <- list(stat_summary(fun = median, geom = "crossbar", width = 0.45, linewidth = MEDIAN_W, colour = "black"))

# One clean "n=" per x-level (and per facet when provided), placed just below data
n_labels_per_facet <- function(data, xvar, yvar, facet_var = NULL, offset_frac = N_OFFSET) {
  x <- rlang::sym(xvar); y <- rlang::sym(yvar)
  d <- data %>% dplyr::filter(is.finite(!!y), !is.na(!!x))
  if (!nrow(d)) return(tibble())
  if (is.null(facet_var)) {
    rng <- range(d %>% dplyr::pull(!!y), na.rm = TRUE)
    y_pos <- rng[1] - offset_frac * diff(rng)
    out <- d %>% dplyr::count(!!x, name = "n") %>% dplyr::mutate(y = y_pos, lbl = paste0("n=", n))
    names(out)[names(out) == rlang::as_name(x)] <- rlang::as_name(x)
    out %>% dplyr::select(!!x, y, lbl)
  } else {
    f <- rlang::sym(facet_var)
    yb <- d %>% dplyr::group_by(!!f) %>%
      dplyr::summarise(ymin = min(!!y, na.rm = TRUE), ymax = max(!!y, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(y = ymin - offset_frac * (ymax - ymin)) %>% dplyr::select(!!f, y)
    out <- d %>% dplyr::count(!!f, !!x, name="n") %>%
      dplyr::left_join(yb, by = setNames(rlang::as_name(f), rlang::as_name(f))) %>%
      dplyr::mutate(lbl = paste0("n=", n))
    names(out)[names(out) == rlang::as_name(x)] <- rlang::as_name(x)
    names(out)[names(out) == rlang::as_name(f)] <- rlang::as_name(f)
    out %>% dplyr::select(!!f, !!x, y, lbl)
  }
}

save_both <- function(p, path_base, w = 10, h = 8, dpi = 600) {
  ggsave(paste0(path_base, ".png"),  p, width = w, height = h, units = "in", dpi = dpi)
  ggsave(paste0(path_base, ".tiff"), p, width = w, height = h, units = "in", dpi = dpi, compression = "lzw")
}

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------ Panel A: Age & Onset -------------------------
dfA <- df %>%
  select(edss_cat, age, age_onset) %>%
  pivot_longer(c(age, age_onset), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = c("age","age_onset"),
                          labels = c("Age (years)", "Age at onset (years)")))
nA <- n_labels_per_facet(dfA, "edss_cat", "value", facet_var = "measure", offset_frac = N_OFFSET)

pA <- ggplot(dfA, aes(x = edss_cat, y = value, fill = edss_cat)) +
  box_layer + dots_layer + median_bar +
  geom_text(data = nA, aes(x = edss_cat, y = y, label = lbl),
            inherit.aes = FALSE, size = N_SIZE, vjust = 1, alpha = 0.9) +
  facet_wrap(~ measure, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = edss_colors) +
  scale_x_discrete(labels = edss_labels) +
  labs(title = "A) Age and age at onset by EDSS category", x = "EDSS category", y = "Years") +
  p_theme +
  coord_cartesian(clip = "off", expand = FALSE) +
  theme(plot.margin = margin(t = 6, r = 6, b = BOTTOM_PAD, l = 6))

# ------------------------ Panel B: Duration & ARR ----------------------
dfB <- df %>%
  select(edss_cat, duration_years, arr) %>%
  pivot_longer(c(duration_years, arr), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure,
                          levels = c("duration_years","arr"),
                          labels = c("Disease duration (years)", "Annualized relapse rate (ARR)")))
nB <- n_labels_per_facet(dfB, "edss_cat", "value", facet_var = "measure", offset_frac = N_OFFSET)

pB <- ggplot(dfB, aes(x = edss_cat, y = value, fill = edss_cat)) +
  box_layer + dots_layer + median_bar +
  geom_text(data = nB, aes(x = edss_cat, y = y, label = lbl),
            inherit.aes = FALSE, size = N_SIZE, vjust = 1, alpha = 0.9) +
  facet_wrap(~ measure, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = edss_colors) +
  scale_x_discrete(labels = edss_labels) +
  labs(title = "B) Disease duration increases while ARR declines", x = "EDSS category", y = NULL) +
  p_theme +
  coord_cartesian(clip = "off", expand = FALSE) +
  theme(plot.margin = margin(t = 6, r = 6, b = BOTTOM_PAD, l = 6))

# ------------------------ Panel C: Vitamin D ---------------------------
has_vtd <- any(is.finite(df$vitamin_d))
nC <- if (has_vtd) n_labels_per_facet(df, "edss_cat", "vitamin_d", facet_var = NULL, offset_frac = N_OFFSET) else tibble()

pC <- ggplot(df, aes(x = edss_cat, y = vitamin_d, fill = edss_cat)) +
  { if (has_vtd) box_layer else list() } +
  { if (has_vtd) dots_layer else list() } +
  { if (has_vtd) median_bar else list() } +
  geom_hline(yintercept = 10, linetype = "dashed") +
  geom_hline(yintercept = 30, linetype = "dashed") +
  annotate("text", x = 3.05, y = 10, label = "10 ng/mL", hjust = 1, vjust = -0.4, size = 3) +
  annotate("text", x = 3.05, y = 30, label = "30 ng/mL", hjust = 1, vjust = -0.4, size = 3) +
  { if (has_vtd && nrow(nC)) geom_text(data = nC, aes(x = edss_cat, y = y, label = lbl),
                                       inherit.aes = FALSE, size = N_SIZE, vjust = 1, alpha = 0.9) else list() } +
  scale_fill_manual(values = edss_colors) +
  scale_x_discrete(labels = edss_labels) +
  labs(title = "C) Serum 25(OH) vitamin D by EDSS category", x = "EDSS category", y = "Vitamin D (ng/mL)") +
  p_theme +
  coord_cartesian(clip = "off", expand = FALSE) +
  theme(plot.margin = margin(t = 6, r = 6, b = BOTTOM_PAD, l = 6)) +
  { if (!has_vtd) annotate("text", x = 2, y = 18, label = "No vitamin D data detected",
                           size = 4, fontface = "bold") else list() }

# ------------------------ Panel D: Lipid ratios ------------------------
dfD <- df %>%
  select(edss_cat, chol_hdl_ratio, ldl_hdl_ratio, aip) %>%
  pivot_longer(c(chol_hdl_ratio, ldl_hdl_ratio, aip), names_to = "ratio", values_to = "value") %>%
  mutate(ratio = factor(ratio,
                        levels = c("chol_hdl_ratio","ldl_hdl_ratio","aip"),
                        labels = c("Total cholesterol / HDL", "LDL / HDL", "Atherogenic Index of Plasma (AIP)")))
nD <- n_labels_per_facet(dfD, "edss_cat", "value", facet_var = "ratio", offset_frac = N_OFFSET)

pD <- ggplot(dfD, aes(x = edss_cat, y = value, fill = edss_cat)) +
  box_layer + dots_layer + median_bar +
  geom_text(data = nD, aes(x = edss_cat, y = y, label = lbl),
            inherit.aes = FALSE, size = N_SIZE, vjust = 1, alpha = 0.9) +
  facet_wrap(~ ratio, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = edss_colors) +
  scale_x_discrete(labels = edss_labels) +
  labs(title = "D) Lipid ratios rise with increasing disability", x = "EDSS category", y = NULL) +
  p_theme +
  coord_cartesian(clip = "off", expand = FALSE) +
  theme(plot.margin = margin(t = 6, r = 6, b = BOTTOM_PAD, l = 6))

# ------------------------ Optional: p-for-trend in titles --------------
trend_p <- function(y, g) {
  if (!requireNamespace("clinfun", quietly = TRUE)) return(NA_real_)
  g <- factor(g, levels = c("Mild (0–2.5)","Moderate (3–5.5)","Severe (≥6)"), ordered = TRUE)
  o <- suppressWarnings(clinfun::jonckheere.test(y, g, alternative = "two.sided"))
  unname(o$p.value)
}
fmtp <- function(p) ifelse(is.na(p), "NA", format.pval(p, eps = 1e-4, digits = 2))

# Compute once; safe if vars missing
p_dur <- trend_p(df$duration_years, df$edss_cat)
p_arr <- trend_p(df$arr,            df$edss_cat)
p_vit <- trend_p(df$vitamin_d,      df$edss_cat)
p_chr <- trend_p(df$chol_hdl_ratio, df$edss_cat)
p_lhr <- trend_p(df$ldl_hdl_ratio,  df$edss_cat)
p_aip <- trend_p(df$aip,            df$edss_cat)

pB <- pB + labs(title = paste0("B) Disease duration increases while ARR declines",
                               "  (p-trend: duration=", fmtp(p_dur),
                               ", ARR=", fmtp(p_arr), ")"))
pC <- pC + labs(title = paste0("C) Serum 25(OH) vitamin D by EDSS category  (p-trend=",
                               fmtp(p_vit), ")"))
pD <- pD + labs(title = paste0("D) Lipid ratios rise with increasing disability  (p-trend: CHR=",
                               fmtp(p_chr), ", LDL/HDL=", fmtp(p_lhr), ", AIP=", fmtp(p_aip), ")"))

# ------------------------ Composite (2x2) ------------------------------
composite <- (pA | pB) / (pC | pD)
composite <- composite + plot_annotation(
  theme = theme(plot.margin = margin(t = 6, r = 6, b = 6, l = 6))
)

# ------------------------ Save: composite + panels ---------------------
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
save_both(composite, "outputs/figures/Figure1_MS_Baseline_colored", w = 10, h = 8, dpi = 600)
save_both(pA, "outputs/figures/Figure1A_Age_Onset_colored",         w = 8, h = 4.8, dpi = 600)
save_both(pB, "outputs/figures/Figure1B_Duration_ARR_colored",      w = 8, h = 4.8, dpi = 600)
save_both(pC, "outputs/figures/Figure1C_VitD_colored",              w = 8, h = 4.8, dpi = 600)
save_both(pD, "outputs/figures/Figure1D_LipidRatios_colored",       w = 8, h = 4.8, dpi = 600)

message("✅ Saved composite + individual panels (PNG & TIFF, 600 DPI) in outputs/figures/")
