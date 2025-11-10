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

# analysis/00_setup/01_readiness_check.R
library(tidyverse); library(readr); library(janitor); library(yaml)
source("analysis/00_setup/functions.R")

cfg <- yaml::read_yaml("analysis/00_setup/config.yml")
df <- readr::read_csv(cfg$data_file) %>% janitor::clean_names()
names(df) <- gsub("\\.", "_", names(df))

# Check critical columns by key
required <- c("age","sex","phenotype","edss","duration_years","vitd","crp","bmi")
missing <- required[!unlist(lapply(required, function(k) cn[[k]] %in% names(df)))]
if (length(missing)) {
  message("Missing required columns (keys): ", paste(missing, collapse=", "))
}

df <- derive_vars(df) %>% clean_levels()

# Completeness table
vars <- c("edss","vitd","crp","hba1c","tc","ldl","hdl","tg","AIP","Chol_HDL_Ratio","LDL_HDL_Ratio","bmi","ARR","duration_years","age","sex","phenotype")
var_map <- tibble::tibble(
  key = vars,
  column = sapply(vars, function(k) ifelse(k %in% names(cn), cn[[k]], k))
)
qc <- var_map %>%
  mutate(
    present = column %in% names(df),
    non_missing_n = ifelse(present, sapply(column, \(x) sum(!is.na(df[[x]]))), NA),
    min = ifelse(
      present,
      sapply(column, \(x) {
        v <- df[[x]]
        if (is.numeric(v)) suppressWarnings(min(v, na.rm = TRUE)) else NA
      }),
      NA
    ),
    max = ifelse(
      present,
      sapply(column, \(x) {
        v <- df[[x]]
        if (is.numeric(v)) suppressWarnings(max(v, na.rm = TRUE)) else NA
      }),
      NA
    )
  )

readr::write_csv(qc, "outputs/supplement/data_quality_completeness.csv")
message("Wrote: outputs/supplement/data_quality_completeness.csv")
