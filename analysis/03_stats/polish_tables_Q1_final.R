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
# Q1 polish for Table 3 & 4 — precompute indices (no tidy-eval)
# ===============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(flextable)
  library(officer)
  library(stringr)
  library(readr)
})

out_dir <- "outputs/tables"

# ---------- load ----------
t3 <- read_csv(file.path(out_dir, "Table3_EDSS_group_tests.csv"), show_col_types = FALSE)
t4 <- read_csv(file.path(out_dir, "Table4_Spearman_EDSS.csv"),     show_col_types = FALSE)

# ---------- normalize names ----------
norm_names <- function(nm) {
  nm %>%
    str_replace_all("[–—]", "-") %>%     # en/em dashes to hyphen
    str_replace_all("\\s+", " ") %>%     # collapse spaces
    str_trim()
}
names(t3) <- norm_names(names(t3))
names(t4) <- norm_names(names(t4))

# ---------- helpers ----------
find_col <- function(df, candidates) {
  nm <- names(df)
  for (pat in candidates) {
    hit <- nm[ str_detect(nm, regex(pat, ignore_case = TRUE)) ]
    if (length(hit)) return(hit[1])
  }
  NA_character_
}
p_to_num <- function(x) {
  # turn strings like "<1e-04" into numeric
  x <- str_replace_all(x, "<\\s*1e-0*4", "0.0001")
  suppressWarnings(as.numeric(x))
}
is_sig_vec <- function(x) {
  y <- p_to_num(x)
  !is.na(y) & y < 0.05
}
num_fmt <- function(x, d = 2) ifelse(is.na(x), NA, formatC(as.numeric(x), digits = d, format = "f"))

# =========================================================
#                          TABLE 3
# =========================================================
# locate columns
kw_col   <- find_col(t3, c("^Kruskal-?Wallis\\s*p$","^Kruskal","Wallis"))
jt_col   <- find_col(t3, c("^Trend\\s*p","\\(JT\\)"))
dir_col  <- find_col(t3, c("^Direction$","^direction$"))
var_col  <- find_col(t3, c("^Variable$"))
n_col    <- find_col(t3, c("^n$","^N$"))

# safe defaults if not found
stopifnot(!is.na(var_col), !is.na(n_col), !is.na(kw_col), !is.na(jt_col))
if (is.na(dir_col)) t3$Direction <- "flat" else names(t3)[names(t3)==dir_col] <- "Direction"

# make arrowed direction
t3$Direction <- case_when(
  str_detect(t3$Direction, "↑|up|\\+|increase") ~ "\u2191 with EDSS",
  str_detect(t3$Direction, "↓|down|\\-|decrease") ~ "\u2193 with EDSS",
  TRUE ~ "flat"
)

# precompute indices to bold
idx_kw_sig <- which(is_sig_vec(t3[[kw_col]]))
idx_jt_sig <- which(is_sig_vec(t3[[jt_col]]))

# assemble display table (keep stable column order)
t3_disp <- t3 %>%
  transmute(
    !!var_col := .data[[var_col]],
    !!n_col   := .data[[n_col]],
    !!kw_col  := .data[[kw_col]],
    !!jt_col  := .data[[jt_col]],
    Direction = Direction
  )

# build flextable
ft3 <- flextable(t3_disp) |>
  autofit() |>
  align(j = 2:ncol(t3_disp), align = "center") |>
  bold(i = idx_kw_sig, j = kw_col) |>
  bold(i = idx_jt_sig, j = jt_col) |>
  set_caption("Table 3. Group-wise EDSS comparisons: Kruskal–Wallis and Jonckheere–Terpstra trend tests.") |>
  theme_booktabs()

save_as_docx(`Table 3` = ft3, path = file.path(out_dir, "Table3_EDSS_group_tests_Q1.docx"))

# =========================================================
#                          TABLE 4
# =========================================================
# locate columns (robust to naming variants)
rho_col       <- find_col(t4, c("^Spearman\\s*rho$","^Spearman","rho"))
p_col         <- find_col(t4, c("^p-value$","^p value$","^p$"))
rho_part_col  <- find_col(t4, c("^Partial\\s*rho","Partial.*rho"))
p_part_col    <- find_col(t4, c("^Partial\\s*p","Partial.*p"))
var_col4      <- find_col(t4, c("^Variable$"))
n_col4        <- find_col(t4, c("^n$","^N$"))

stopifnot(!is.na(var_col4), !is.na(n_col4), !is.na(rho_col), !is.na(p_col))
# if optional partial columns miss, create placeholders
if (is.na(rho_part_col)) { t4$`Partial rho (age±sex)` <- NA_real_; rho_part_col <- "Partial rho (age±sex)"} 
if (is.na(p_part_col))   { t4$`Partial p` <- NA_character_;           p_part_col  <- "Partial p" }

# format rho columns to 2dp
t4[[rho_col]]      <- num_fmt(t4[[rho_col]], 2)
t4[[rho_part_col]] <- num_fmt(t4[[rho_part_col]], 2)

# precompute indices to bold
idx_p_sig     <- which(str_detect(t4[[p_col]], "<|^0\\.|^0,"))
idx_ppart_sig <- which(str_detect(t4[[p_part_col]], "<|^0\\.|^0,"))

# assemble display table
t4_disp <- t4 %>%
  transmute(
    !!var_col4     := .data[[var_col4]],
    !!n_col4       := .data[[n_col4]],
    !!rho_col      := .data[[rho_col]],
    !!p_col        := .data[[p_col]],
    !!rho_part_col := .data[[rho_part_col]],
    !!p_part_col   := .data[[p_part_col]]
  )

ft4 <- flextable(t4_disp) |>
  autofit() |>
  align(j = 2:ncol(t4_disp), align = "center") |>
  bold(i = idx_p_sig,     j = p_col) |>
  bold(i = idx_ppart_sig, j = p_part_col) |>
  set_caption("Table 4. Spearman correlations of EDSS with clinical/biochemical variables; partial rho adjusts for age and sex.") |>
  theme_booktabs()

save_as_docx(`Table 4` = ft4, path = file.path(out_dir, "Table4_Spearman_EDSS_Q1.docx"))

message("Q1-styled DOCX exported to: ", out_dir)
