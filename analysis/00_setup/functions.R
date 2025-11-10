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

# analysis/00_setup/functions.R
library(tidyverse)
library(yaml)

cfg <- yaml::read_yaml("analysis/00_setup/config.yml")
cn <- cfg$columns

col <- function(df, key) {
  nm <- cn[[key]]
  if (!nm %in% names(df)) warning(sprintf("Column '%s' (%s) not found in data", key, nm))
  nm
}

derive_vars <- function(df) {
  # ARR
  if (!"ARR" %in% names(df) && all(c(col(df,"relapses_total"), col(df,"duration_years")) %in% names(df))) {
    dur <- df[[ col(df,"duration_years") ]]
    df$ARR <- with(df, df[[ col(df,"relapses_total") ]] / ifelse(dur > 0, dur, NA_real_))
  }
  # EDSS categories
  edss_v <- df[[ col(df,"edss") ]]
  df$EDSS_cat <- dplyr::case_when(
    edss_v <= cfg$edss_categories$mild_max ~ "Mild (0–2.5)",
    edss_v <= cfg$edss_categories$moderate_max ~ "Moderate (3–5.5)",
    is.finite(edss_v) ~ "Severe (≥6.0)",
    TRUE ~ NA_character_
  )
  # Ratios/AIP if absent
  if (!"Chol_HDL_Ratio" %in% names(df) && all(c(col(df,"tc"), col(df,"hdl")) %in% names(df)))
    df$Chol_HDL_Ratio <- df[[ col(df,"tc") ]] / df[[ col(df,"hdl") ]]
  if (!"LDL_HDL_Ratio" %in% names(df) && all(c(col(df,"ldl"), col(df,"hdl")) %in% names(df)))
    df$LDL_HDL_Ratio <- df[[ col(df,"ldl") ]] / df[[ col(df,"hdl") ]]
  if (!"AIP" %in% names(df) && all(c(col(df,"tg"), col(df,"hdl")) %in% names(df)))
    df$AIP <- log10(df[[ col(df,"tg") ]]/ df[[ col(df,"hdl") ]])
  df
}

clean_levels <- function(df) {
  if (col(df,"sex") %in% names(df)) {
    df[[ col(df,"sex") ]] <- forcats::fct_relevel(as.factor(df[[ col(df,"sex") ]]), "Female")
  }
  if (col(df,"phenotype") %in% names(df)) {
    df[[ col(df,"phenotype") ]] <- as.factor(df[[ col(df,"phenotype") ]])
  }
  df
}

num_summ <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)==0) return(tibble::tibble(n=0, mean=NA_real_, sd=NA_real_, median=NA_real_, IQR=NA_real_, min=NA_real_, max=NA_real_))
  q <- stats::quantile(x, probs=c(.25,.75), na.rm=TRUE)
  tibble::tibble(
    n=length(x),
    mean=mean(x),
    sd=stats::sd(x),
    median=stats::median(x),
    IQR=diff(q),
    min=min(x),
    max=max(x)
  )
}
