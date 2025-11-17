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
# Subsection 3.2 – Group differences across EDSS
# Tests: Kruskal–Wallis + Jonckheere–Terpstra (trend)
# Exports: Table3_EDSS_group_tests_Q1.docx + Figure2A–D
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(janitor)
  library(ggplot2)
  library(ggpubr)
  library(DescTools)    # for JonckheereTerpstraTest
  library(flextable)
  library(officer)
})

# ---------- Paths ----------
base_dir <- getwd()
data_file <- file.path(base_dir, "Libya_MS.csv")
out_dir   <- file.path(base_dir, "outputs", "subsection3_2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- Load & clean ----------
df <- read_csv(data_file, show_col_types = FALSE) |> janitor::clean_names()
names(df) <- gsub("\\.", "_", names(df))

# ---------- Define variables ----------
kw_vars <- c("age", "age_of_onset_of_disease_years", "duration_of_illness_years",
             "arr", "bmi_kg_m2", "vitamin_d", "total_cholesterol",
             "ldl", "hdl", "tg", "total_cholesterol_hdl_ratio",
             "ldl_hdl_ratio", "aip")

# Check EDSS and group variable
if(!"edss" %in% names(df)) stop("No EDSS column found.")
df <- df %>%
  mutate(edss_cat = case_when(
    edss <= 2.5 ~ "Mild (0–2.5)",
    edss > 2.5 & edss <= 5.5 ~ "Moderate (3–5.5)",
    edss > 5.5 ~ "Severe (≥6)"
  )) %>%
  filter(!is.na(edss_cat))

# ---------- Run Kruskal–Wallis and Jonckheere–Terpstra ----------
res_list <- lapply(kw_vars, function(v) {
  x <- df[[v]]
  g <- df$edss_cat
  kw_p  <- tryCatch(kruskal.test(x ~ g)$p.value, error = function(e) NA)
  jt_p  <- tryCatch(JonckheereTerpstraTest(x, g, alternative="two.sided")$p.value, error = function(e) NA)
  trend <- ifelse(is.na(jt_p), "", ifelse(kw_p < 0.05 & jt_p < 0.05,
                                          ifelse(mean(x[g=="Mild (0–2.5)"], na.rm=TRUE) < mean(x[g=="Severe (≥6)"], na.rm=TRUE),
                                                 "↑ with EDSS", "↓ with EDSS"), "flat"))
  data.frame(Variable=v,
             Kruskal_Wallis_p=kw_p,
             Trend_p_JT=jt_p,
             Direction=trend)
}) |> bind_rows()

# Round p-values to 4 digits
res_list <- res_list %>%
  mutate(across(contains("p"), ~formatC(., digits=4, format="f")))

# ---------- Export Table to DOCX ----------
ft <- flextable(res_list)
ft <- autofit(ft)
doc <- read_docx() %>%
  body_add_par("Table 3. Kruskal–Wallis and Jonckheere–Terpstra tests across EDSS categories", style="heading 1") %>%
  body_add_flextable(ft)
print(doc, target=file.path(out_dir, "Table3_EDSS_group_tests_Q1.docx"))

# ---------- Visualization ----------
num_vars <- c("age", "duration_of_illness_years", "bmi_kg_m2", "vitamin_d",
              "hdl", "ldl_hdl_ratio", "aip")
plots <- lapply(num_vars, function(v){
  ggboxplot(df, x="edss_cat", y=v, color="edss_cat", palette="jco",
            add="jitter", ylab=v, xlab="EDSS category") +
    theme_minimal(base_size=13) +
    theme(legend.position="none")
})
fig <- ggarrange(plotlist=plots[1:4], ncol=2, nrow=2, labels=c("A","B","C","D"))
ggsave(file.path(out_dir,"Figure2A-D_GroupTests.png"), fig, dpi=600, width=10, height=8)

message("\nSubsection 3.2 outputs saved to: ", out_dir)
