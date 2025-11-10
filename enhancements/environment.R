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

# Environment ----
# Use the same R you cite in Methods (choose ONE and keep it everywhere)
R.version.string
# Target: "R version 4.4.1 (2024-06-14)"  # <- if different, update Methods OR re-run on 4.4.1

set.seed(1234)

# Packages
library(dplyr); library(readr); library(janitor)
library(MASS)       # polr
library(ordinal)    # clm (optional; for proportional-odds test)
library(performance)# r2_nagelkerke
library(broom)
library(rstatix)    # kruskal_test, cor_test helpers
library(clinfun)    # jonckheere.test
library(flextable); library(officer) # for docx export via flextable::save_as_docx
# Data ----
df <- read_csv("Libya_MS.csv", show_col_types = FALSE) |> clean_names()

# Columns (robust)
pick <- function(...) { x <- c(...); x[x %in% names(df)][1] }
edss_col <- pick("edss")
age_col  <- pick("age")
dur_col  <- pick("duration_of_illness_years","duration_years")
bmi_col  <- pick("bmi_kg_m2","bmi")
vitd_col <- pick("vitamin_d","vit_d","vitamin_d_ng_ml")
hdl_col  <- pick("hdl","hdl_mg_dl")
ldl_col  <- pick("ldl","ldl_mg_dl")
tg_col   <- pick("tg","triglycerides","triglycerides_mg_dl")
sex_col  <- pick("gender","sex")

stopifnot(!is.na(edss_col), !is.na(age_col), !is.na(dur_col), !is.na(bmi_col),
          !is.na(vitd_col), !is.na(hdl_col), !is.na(ldl_col), !is.na(tg_col), !is.na(sex_col))

# Derived: LDL/HDL, AIP
df <- df |>
  mutate(
    edss = as.numeric(.data[[edss_col]]),
    age  = as.numeric(.data[[age_col]]),
    duration = as.numeric(.data[[dur_col]]),
    bmi = as.numeric(.data[[bmi_col]]),
    vitamin_d = as.numeric(.data[[vitd_col]]),
    hdl = as.numeric(.data[[hdl_col]]),
    ldl = as.numeric(.data[[ldl_col]]),
    tg  = as.numeric(.data[[tg_col]]),
    ldl_hdl = ldl/hdl,
    aip = log10(tg/hdl),
    sex_raw = .data[[sex_col]],
    sex = case_when(
      tolower(as.character(sex_raw)) %in% c("male","m","1") ~ "Male",
      tolower(as.character(sex_raw)) %in% c("female","f","0") ~ "Female",
      TRUE ~ as.character(sex_raw)
    )
  ) |>
  mutate(
    edss_cat = case_when(
      edss <= 2.5 ~ "Mild (0–2.5)",
      edss > 2.5 & edss <= 5.5 ~ "Moderate (3–5.5)",
      edss > 5.5 ~ "Severe (≥6)"
    ),
    edss_cat = factor(edss_cat, levels = c("Mild (0–2.5)","Moderate (3–5.5)","Severe (≥6)"), ordered = TRUE),
    sex = factor(sex, levels = c("Male","Female"))
  ) |>
  filter(is.finite(edss), is.finite(age), is.finite(duration), is.finite(bmi),
         is.finite(vitamin_d), is.finite(hdl), is.finite(ldl_hdl), is.finite(aip), !is.na(sex))


cont_vars <- c("age","duration","bmi","vitamin_d","hdl","ldl","tg","ldl_hdl","aip")
kw <- lapply(cont_vars, \(v) {
  k <- kruskal_test(as.formula(paste(v, "~ edss_cat")), data = df)
  tibble(Variable = v, `Kruskal–Wallis p` = signif(k$p, 4))
}) |> bind_rows()

jt <- lapply(cont_vars, \(v) {
  # Convert EDSS order to scores 1,2,3 for JT
  edss_score <- as.numeric(df$edss_cat)
  j <- jonckheere.test(df[[v]], edss_score, alternative = "two.sided")
  tibble(Variable = v, `Trend p (Jonckheere–Terpstra)` = signif(j$p.value, 4),
         Direction = dplyr::case_when(
           cor(df[[v]], edss_score, use="complete.obs", method="spearman") >  0.05 ~ "↑ with EDSS",
           cor(df[[v]], edss_score, use="complete.obs", method="spearman") < -0.05 ~ "↓ with EDSS",
           TRUE ~ "flat"
         ))
}) |> bind_rows()

table3 <- left_join(kw, jt, by = "Variable")


ft3 <- flextable(table3 |> mutate(
  Variable = recode(Variable,
                    age="Age (years)", duration="Disease duration (years)", bmi="BMI (kg/m²)",
                    vitamin_d="25(OH) Vitamin D (ng/mL)", hdl="HDL (mg/dL)", ldl="LDL (mg/dL)",
                    tg="Triglycerides (mg/dL)", ldl_hdl="LDL/HDL ratio", aip="AIP"
  )
))
save_as_docx("Table 3. EDSS group tests and trends" = ft3, path = "outputs/tables/Table3_EDSS_group_tests_Q1.docx")

vars <- c("age","duration","bmi","vitamin_d","hdl","ldl","tg","ldl_hdl","aip")
spearman <- lapply(vars, \(v) {
  s <- cor.test(df[[v]], df$edss, method = "spearman", exact = FALSE)
  tibble(Variable = v, `Spearman rho` = round(as.numeric(s$estimate), 2),
         `p value` = format.pval(s$p.value, digits = 4))
}) |> bind_rows()

# Partial Spearman via residualization (simple & robust):
residize <- function(y, X) lm(y ~ ., data = as.data.frame(X))$residuals
edss_res <- residize(df$edss, model.matrix(~ age + sex, df)[, -1])
par_list <- lapply(vars, \(v) {
  x_res <- residize(df[[v]], model.matrix(~ age + sex, df)[, -1])
  s <- cor.test(x_res, edss_res, method = "spearman", exact = FALSE)
  tibble(Variable = v, `Partial rho (age+sex)` = round(as.numeric(s$estimate), 2),
         `Partial p value` = format.pval(s$p.value, digits = 4))
}) |> bind_rows()

table4 <- left_join(spearman, par_list, by = "Variable") |>
  mutate(Variable = recode(Variable,
    age="Age (years)", duration="Disease duration (years)", bmi="BMI (kg/m²)",
    vitamin_d="25(OH) Vitamin D (ng/mL)", hdl="HDL (mg/dL)", ldl="LDL (mg/dL)",
    tg="Triglycerides (mg/dL)", ldl_hdl="LDL/HDL ratio", aip="AIP"
  ))
############################
ft4 <- flextable(table4)
save_as_docx("Table 4. Correlations with EDSS" = ft4, path = "outputs/tables/Table4_Spearman_EDSS_Q1.docx")


##########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
vars <- c("age","duration","bmi","vitamin_d","hdl","ldl","tg","ldl_hdl","aip")
spearman <- lapply(vars, \(v) {
  s <- cor.test(df[[v]], df$edss, method = "spearman", exact = FALSE)
  tibble(Variable = v, `Spearman rho` = round(as.numeric(s$estimate), 2),
         `p value` = format.pval(s$p.value, digits = 4))
}) |> bind_rows()

# Partial Spearman via residualization (simple & robust):
residize <- function(y, X) lm(y ~ ., data = as.data.frame(X))$residuals
edss_res <- residize(df$edss, model.matrix(~ age + sex, df)[, -1])
par_list <- lapply(vars, \(v) {
  x_res <- residize(df[[v]], model.matrix(~ age + sex, df)[, -1])
  s <- cor.test(x_res, edss_res, method = "spearman", exact = FALSE)
  tibble(Variable = v, `Partial rho (age+sex)` = round(as.numeric(s$estimate), 2),
         `Partial p value` = format.pval(s$p.value, digits = 4))
}) |> bind_rows()

table4 <- left_join(spearman, par_list, by = "Variable") |>
  mutate(Variable = recode(Variable,
                           age="Age (years)", duration="Disease duration (years)", bmi="BMI (kg/m²)",
                           vitamin_d="25(OH) Vitamin D (ng/mL)", hdl="HDL (mg/dL)", ldl="LDL (mg/dL)",
                           tg="Triglycerides (mg/dL)", ldl_hdl="LDL/HDL ratio", aip="AIP"
  ))

ft4 <- flextable(table4)
save_as_docx("Table 4. Correlations with EDSS" = ft4, path = "outputs/tables/Table4_Spearman_EDSS_Q1.docx")



#%%%%%%%%%%%%%%%%%%%%%%%
# Standardize continuous predictors
std <- function(x) as.numeric(scale(x))
mod_df <- df |>
  transmute(
    edss_cat, sex,
    age = std(age), duration = std(duration), bmi = std(bmi),
    vitamin_d = std(vitamin_d), hdl = std(hdl),
    ldl_hdl = std(ldl_hdl), aip = std(aip)
  )

# Primary model: choose ONE lipid representation set
# Option A (recommended): HDL + LDL/HDL (drop AIP to reduce collinearity)
form <- edss_cat ~ age + duration + bmi + vitamin_d + hdl + ldl_hdl + sex

fit <- suppressWarnings(MASS::polr(form, data = mod_df, Hess = TRUE))
nk <- performance::r2_nagelkerke(fit)$R2

# ORs and profile-likelihood 95% CI
ci <- confint(fit)
or <- exp(coef(fit)); ci_or <- exp(ci[!grepl("\\|", rownames(ci)), , drop=FALSE])

table5 <- tibble(
  term = names(or),
  OR = as.numeric(or),
  CI_low = ci_or[,1],
  CI_high = ci_or[,2]
) |>
  mutate(Predictor = recode(term,
                            "age"="Age (per SD)","duration"="Disease duration (per SD)","bmi"="BMI (per SD)",
                            "vitamin_d"="25(OH) Vitamin D (per SD)","hdl"="HDL (per SD)",
                            "ldl_hdl"="LDL/HDL ratio (per SD)","sexFemale"="Female vs Male"
  )) |>
  transmute(Predictor,
            `OR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", OR, CI_low, CI_high))

ft5 <- flextable(table5)
save_as_docx("Table 5. Ordinal logistic regression (EDSS categories)" = ft5,
             path = "outputs/tables/Table5_Ordinal_Logistic_EDSS.docx")

# Partial R² by drop-one ΔNagelkerke
drop_r2 <- function(term){
  rf <- update(form, paste(". ~ . -", term))
  rfit <- suppressWarnings(MASS::polr(rf, data = mod_df, Hess = TRUE))
  full <- performance::r2_nagelkerke(fit)$R2
  red  <- performance::r2_nagelkerke(rfit)$R2
  max(0, full - red)
}
terms <- c("age","duration","bmi","vitamin_d","hdl","ldl_hdl","sex")
part <- sapply(terms, drop_r2)
table6 <- tibble(
  Parameter = c("Age","Disease duration","BMI","25(OH) Vitamin D","HDL","LDL/HDL ratio","Sex (Female)"),
  `Partial R²` = round(part, 3)
) |> arrange(desc(`Partial R²`))

ft6 <- flextable(table6)
save_as_docx("Table 6. Variable importance (Partial R²)" = ft6,
             path = "outputs/tables/Table6_PartialR2_EDSS.docx")

# Forest + Partial R² figure (colored)
library(ggplot2); library(ggpubr); library(scales)
fdf <- table5 |> mutate(OR = sub(" .*","", `OR (95% CI)`)) # not used for plotting; rebuild from raw:

plot_df <- tibble(
  Predictor = c("Age (per SD)","Disease duration (per SD)","BMI (per SD)",
                "25(OH) Vitamin D (per SD)","HDL (per SD)","LDL/HDL ratio (per SD)","Female vs Male"),
  OR = as.numeric(or),
  L = ci_or[,1], H = ci_or[,2]
)

pA <- ggplot(plot_df, aes(y=fct_reorder(Predictor, OR), x=OR))+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_errorbarh(aes(xmin=L, xmax=H), height=.2)+
  geom_point(size=3)+
  scale_x_log10(labels=label_number(accuracy=0.01))+
  labs(title="Independent Predictors of Disability (EDSS categories)",
       x="Adjusted Odds Ratio (log scale)", y=NULL)+
  theme_minimal(base_size=12)

pB <- ggplot(table6, aes(x=reorder(Parameter, `Partial R²`), y=`Partial R²`))+
  geom_col() + coord_flip()+
  labs(title="Relative Contribution to Model Fit (Partial R²)", x=NULL, y="Partial R²")+
  theme_minimal(base_size=12)

fig3 <- ggarrange(pA, pB, ncol=2, labels=c("A","B"))
ggsave("outputs/figures/Figure3_MS_Disability_Predictors.png", fig3, dpi=600, width=12, height=6)
