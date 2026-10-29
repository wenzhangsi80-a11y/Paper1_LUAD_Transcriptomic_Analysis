##GitHub2_Paper1_tumor-only WGCNA_5

# Resume Workflow-04-2D WGCNA
# ============================================================

work_dir <- "input_dir/GitHub2_Paper1_tumor-only WGCNA"

setwd(work_dir)

getwd()

load(
  "input_dir/Paper1_ME_brown_Gender_checkpoint_Final.RData"
)

ls()

# ============================================================
# Step 12A. Check Stage_trait structure, sample size,
#            and event distribution
# ============================================================

# ------------------------------------------------------------
# 12A-1. Check variable class and factor levels
# ------------------------------------------------------------

class(cox_design_complete$Stage_trait)

levels(cox_design_complete$Stage_trait)

is.factor(cox_design_complete$Stage_trait)


# ------------------------------------------------------------
# 12A-2. Overall distribution of Stage_trait
# ------------------------------------------------------------

table(cox_design_complete$Stage_trait)

prop.table(
  table(cox_design_complete$Stage_trait)
)


# ------------------------------------------------------------
# 12A-3. Stage × OS event table
# ------------------------------------------------------------

stage_event_table <- table(
  Stage = cox_design_complete$Stage_trait,
  OS_event = cox_design_complete$OS_event
)

stage_event_table


# ------------------------------------------------------------
# 12A-4. Stage-specific event rate
# ------------------------------------------------------------

stage_event_summary <- data.frame(
  Stage = levels(cox_design_complete$Stage_trait),
  N = as.numeric(
    table(cox_design_complete$Stage_trait)
  ),
  Events = as.numeric(
    table(
      cox_design_complete$Stage_trait,
      cox_design_complete$OS_event
    )[, "1"]
  ),
  Non_events = as.numeric(
    table(
      cox_design_complete$Stage_trait,
      cox_design_complete$OS_event
    )[, "0"]
  )
)

stage_event_summary$Event_rate <-
  stage_event_summary$Events /
  stage_event_summary$N

stage_event_summary


# ------------------------------------------------------------
# 12A-5. Check the reference level
# ------------------------------------------------------------

contrasts(cox_design_complete$Stage_trait)


# ------------------------------------------------------------
# 12A-6. Check missing values
# ------------------------------------------------------------

sum(is.na(cox_design_complete$Stage_trait))

table(
  is.na(cox_design_complete$Stage_trait)
)


# ------------------------------------------------------------
# 12A-7. Check consistency with the current Cox analysis set
# ------------------------------------------------------------

c(
  Total_n = nrow(cox_design_complete),
  Total_events = sum(cox_design_complete$OS_event == 1),
  Total_non_events = sum(cox_design_complete$OS_event == 0)
)


# ============================================================
# Step 12A-2. Confirm Stage distribution in the actual
#             multivariable Cox analysis set (n = 485)
# ============================================================

# ------------------------------------------------------------
# 12A-2.1 Identify the exact observations used by the
#         finalized Gender model
# ------------------------------------------------------------

gender_model_rows <- as.integer(
  rownames(model.frame(fit_gender_sqrt))
)

length(gender_model_rows)

# ------------------------------------------------------------
# 12A-2.2 Create the exact 485-patient analysis dataset
# ------------------------------------------------------------

stage_analysis_data <- cox_design_complete[
  gender_model_rows,
]

# ------------------------------------------------------------
# 12A-2.3 Check dimensions and events
# ------------------------------------------------------------

c(
  N = nrow(stage_analysis_data),
  Events = sum(stage_analysis_data$OS_event == 1),
  Non_events = sum(stage_analysis_data$OS_event == 0)
)

# ------------------------------------------------------------
# 12A-2.4 Check Stage missingness in the exact analysis set
# ------------------------------------------------------------

sum(is.na(stage_analysis_data$Stage_trait))

table(
  is.na(stage_analysis_data$Stage_trait)
)

# ------------------------------------------------------------
# 12A-2.5 Stage distribution in the exact 485-patient set
# ------------------------------------------------------------

table(
  stage_analysis_data$Stage_trait,
  useNA = "ifany"
)

# ------------------------------------------------------------
# 12A-2.6 Stage × event table in the exact 485-patient set
# ------------------------------------------------------------

stage_analysis_event_table <- table(
  Stage = stage_analysis_data$Stage_trait,
  OS_event = stage_analysis_data$OS_event,
  useNA = "ifany"
)

stage_analysis_event_table

# ------------------------------------------------------------
# 12A-2.7 Stage-specific summary
# ------------------------------------------------------------

stage_analysis_summary <- data.frame(
  Stage = levels(stage_analysis_data$Stage_trait),
  N = as.numeric(
    table(stage_analysis_data$Stage_trait)
  ),
  Events = as.numeric(
    table(
      stage_analysis_data$Stage_trait,
      stage_analysis_data$OS_event
    )[, "1"]
  ),
  Non_events = as.numeric(
    table(
      stage_analysis_data$Stage_trait,
      stage_analysis_data$OS_event
    )[, "0"]
  )
)

stage_analysis_summary$Event_rate <-
  stage_analysis_summary$Events /
  stage_analysis_summary$N

stage_analysis_summary


# ============================================================
# Step 12A-3. Define the exact complete-case dataset
#             for Stage multivariable Cox analysis
# ============================================================

# ------------------------------------------------------------
# 12A-3.1 Check missingness of all variables required for
#         the Stage Cox model
# ------------------------------------------------------------

stage_required_vars <- c(
  "OS_time",
  "OS_event",
  "ME_brown",
  "Age_trait",
  "Gender_trait",
  "Stage_trait"
)

sapply(
  cox_design_complete[, stage_required_vars],
  function(x) sum(is.na(x))
)

# ------------------------------------------------------------
# 12A-3.2 Construct true complete-case dataset
# ------------------------------------------------------------

stage_complete_idx <- complete.cases(
  cox_design_complete[, stage_required_vars]
)

stage_complete_data <- cox_design_complete[
  stage_complete_idx,
]

# ------------------------------------------------------------
# 12A-3.3 Check final Stage analysis sample
# ------------------------------------------------------------

c(
  N = nrow(stage_complete_data),
  Events = sum(stage_complete_data$OS_event),
  Non_events = sum(stage_complete_data$OS_event == 0)
)


# ------------------------------------------------------------
# 12A-3.4 Stage distribution
# ------------------------------------------------------------

table(
  stage_complete_data$Stage_trait
)

# ------------------------------------------------------------
# 12A-3.5 Stage × event
# ------------------------------------------------------------

stage_complete_event_table <- table(
  Stage = stage_complete_data$Stage_trait,
  OS_event = stage_complete_data$OS_event
)

stage_complete_event_table

# ------------------------------------------------------------
# 12A-3.6 Stage-specific summary
# ------------------------------------------------------------

stage_complete_summary <- data.frame(
  Stage = levels(stage_complete_data$Stage_trait),
  N = as.numeric(
    table(stage_complete_data$Stage_trait)
  ),
  Events = as.numeric(
    table(
      stage_complete_data$Stage_trait,
      stage_complete_data$OS_event
    )[, "1"]
  ),
  Non_events = as.numeric(
    table(
      stage_complete_data$Stage_trait,
      stage_complete_data$OS_event
    )[, "0"]
  )
)

stage_complete_summary$Event_rate <-
  stage_complete_summary$Events /
  stage_complete_summary$N

stage_complete_summary


# ============================================================
# Step 12B. Stage-specific PH assessment
# ============================================================
# ============================================================
# Step 12B-1. Prepare exact Stage analysis dataset
# ============================================================

stage_complete_data$ME_brown_SD <-
  as.numeric(scale(
    stage_complete_data$ME_brown
  ))

stage_complete_data$Gender_Male <-
  ifelse(
    stage_complete_data$Gender_trait == "Male",
    1,
    0
  )

# Check
c(
  N = nrow(stage_complete_data),
  Events = sum(stage_complete_data$OS_event),
  Non_events = sum(stage_complete_data$OS_event == 0)
)

table(stage_complete_data$Stage_trait)

# ============================================================
# Step 12B-2. Conventional Cox model for Stage
# ============================================================

stage_cox_conventional <- coxph(
  Surv(OS_time, OS_event) ~
    ME_brown_SD +
    tt(ME_brown_SD) +
    Age_trait +
    Gender_Male +
    tt(Gender_Male) +
    Stage_trait,
  data = stage_complete_data,
  ties = "efron",
  tt = function(x, t, ...) {
    x * (sqrt(t) / sqrt(1000))
  }
)

summary(stage_cox_conventional)

# ============================================================
# Step 12B-3. Conventional Cox model for Stage PH diagnosis
# ============================================================

stage_ph_diag_model <- coxph(
  Surv(OS_time, OS_event) ~
    ME_brown_SD +
    Age_trait +
    Gender_Male +
    Stage_trait,
  data = stage_complete_data,
  ties = "efron"
)

summary(stage_ph_diag_model)

# ============================================================
# Step 12B-4. Schoenfeld residual PH test
# ============================================================

stage_zph <- cox.zph(
  stage_ph_diag_model,
  transform = "km"
)

stage_zph

# ============================================================
# Step 12B-5. PH test P values
# ============================================================

stage_zph$table

#Step 12B-6：Pay special attention to these 4 lines when running:

stage_zph

# ============================================================
# Step 12B-7. Schoenfeld residual plots
# ============================================================

par(mfrow = c(2, 2))

plot(
  stage_zph,
  var = "Stage_traitStage II"
)

plot(
  stage_zph,
  var = "Stage_traitStage III"
)

plot(
  stage_zph,
  var = "Stage_traitStage IV"
)

plot(
  stage_zph,
  var = "GLOBAL"
)

par(mfrow = c(1, 1))


# ============================================================
# Step 12B-8. Plot the overall Stage PH pattern
# ============================================================
library(ggplot2)
plot(
  stage_zph,
  var = "Stage_trait"
)

# View the variable names in the `cox.zph` object.
rownames(stage_zph$table)

# Examine the structure of the cox.zph object
names(stage_zph)

# ============================================================
# Step 12B-8. Inspect Stage PH diagnostic object
# ============================================================

rownames(stage_zph$table)

names(stage_zph)

# Overall Stage Schoenfeld residual pattern
plot(
  stage_zph,
  var = "Stage_trait"
)

# ============================================================
# Step 12B-9. Inspect the underlying Schoenfeld residuals
# ============================================================

dim(stage_zph$y)

colnames(stage_zph$y)

head(stage_zph$y)

# ============================================================
# Step 12B-10. Inspect transformed event-time information
# ============================================================

head(stage_zph$x)

summary(stage_zph$x)

### The most robust approach is to construct three separate PH diagnostic models, each based on a binary stage contrast.
# ============================================================
# Step 12B-11A revised. True Stage II vs Stage I
# ============================================================

stage_II_I_data <- subset(
  stage_complete_data,
  Stage_trait %in% c("Stage I", "Stage II")
)

stage_II_I_data$Stage_II <- ifelse(
  stage_II_I_data$Stage_trait == "Stage II",
  1,
  0
)

table(stage_II_I_data$Stage_trait)
table(stage_II_I_data$Stage_II)

stage_II_I_diag <- coxph(
  Surv(OS_time, OS_event) ~
    ME_brown_SD +
    Age_trait +
    Gender_Male +
    Stage_II,
  data = stage_II_I_data,
  ties = "efron"
)

stage_II_I_zph <- cox.zph(
  stage_II_I_diag,
  transform = "km"
)

stage_II_I_zph

# ============================================================
# Step 12B-11B revised. True Stage III vs Stage I
# ============================================================

stage_III_I_data <- subset(
  stage_complete_data,
  Stage_trait %in% c("Stage I", "Stage III")
)

stage_III_I_data$Stage_III <- ifelse(
  stage_III_I_data$Stage_trait == "Stage III",
  1,
  0
)

table(stage_III_I_data$Stage_trait)
table(stage_III_I_data$Stage_III)

stage_III_I_diag <- coxph(
  Surv(OS_time, OS_event) ~
    ME_brown_SD +
    Age_trait +
    Gender_Male +
    Stage_III,
  data = stage_III_I_data,
  ties = "efron"
)

stage_III_I_zph <- cox.zph(
  stage_III_I_diag,
  transform = "km"
)

stage_III_I_zph

# ============================================================
# Step 12B-11C revised. True Stage IV vs Stage I
# ============================================================

stage_IV_I_data <- subset(
  stage_complete_data,
  Stage_trait %in% c("Stage I", "Stage IV")
)

stage_IV_I_data$Stage_IV <- ifelse(
  stage_IV_I_data$Stage_trait == "Stage IV",
  1,
  0
)

table(stage_IV_I_data$Stage_trait)
table(stage_IV_I_data$Stage_IV)

stage_IV_I_diag <- coxph(
  Surv(OS_time, OS_event) ~
    ME_brown_SD +
    Age_trait +
    Gender_Male +
    Stage_IV,
  data = stage_IV_I_data,
  ties = "efron"
)

stage_IV_I_zph <- cox.zph(
  stage_IV_I_diag,
  transform = "km"
)

stage_IV_I_zph


###Step 12C：Stage III time-dependent candidate models
# ============================================================
# Step 12C-1. Prepare Stage dummy variables
# ============================================================

stage_complete_data$Stage_II <- ifelse(
  stage_complete_data$Stage_trait == "Stage II",
  1,
  0
)

stage_complete_data$Stage_III <- ifelse(
  stage_complete_data$Stage_trait == "Stage III",
  1,
  0
)

stage_complete_data$Stage_IV <- ifelse(
  stage_complete_data$Stage_trait == "Stage IV",
  1,
  0
)

# Check
table(
  stage_complete_data$Stage_trait,
  stage_complete_data$Stage_II
)

table(
  stage_complete_data$Stage_trait,
  stage_complete_data$Stage_III
)

table(
  stage_complete_data$Stage_trait,
  stage_complete_data$Stage_IV
)

# ============================================================
# Step 12C-2A. Stage III linear time-dependent model
# ============================================================

fit_stageIII_linear <- coxph(
  Surv(OS_time, OS_event) ~
    ME_brown_SD +
    tt(ME_brown_SD) +
    Age_trait +
    Gender_Male +
    tt(Gender_Male) +
    Stage_II +
    Stage_III +
    tt(Stage_III) +
    Stage_IV,
  data = stage_complete_data,
  ties = "efron",
  tt = function(x, t, ...) {
    x * (t / 1000)
  }
)

summary(fit_stageIII_linear)

# ============================================================
# Step 12C-2B. Stage III log1p time-dependent model
# ============================================================

fit_stageIII_log <- coxph(
  Surv(OS_time, OS_event) ~
    ME_brown_SD +
    tt(ME_brown_SD) +
    Age_trait +
    Gender_Male +
    tt(Gender_Male) +
    Stage_II +
    Stage_III +
    tt(Stage_III) +
    Stage_IV,
  data = stage_complete_data,
  ties = "efron",
  tt = function(x, t, ...) {
    x * log1p(t)
  }
)

summary(fit_stageIII_log)

# ============================================================
# Step 12C-2C. Stage III sqrt time-dependent model
# ============================================================

fit_stageIII_sqrt <- coxph(
  Surv(OS_time, OS_event) ~
    ME_brown_SD +
    tt(ME_brown_SD) +
    Age_trait +
    Gender_Male +
    tt(Gender_Male) +
    Stage_II +
    Stage_III +
    tt(Stage_III) +
    Stage_IV,
  data = stage_complete_data,
  ties = "efron",
  tt = function(x, t, ...) {
    x * (sqrt(t) / sqrt(1000))
  }
)

summary(fit_stageIII_sqrt)

# ============================================================
# Step 12C-3. Compare Stage III time-dependent models
# ============================================================

stageIII_model_comparison <- data.frame(
  Model = c(
    "Linear",
    "Log1p",
    "Sqrt"
  ),
  AIC = c(
    AIC(fit_stageIII_linear),
    AIC(fit_stageIII_log),
    AIC(fit_stageIII_sqrt)
  ),
  LogLik = c(
    as.numeric(logLik(fit_stageIII_linear)),
    as.numeric(logLik(fit_stageIII_log)),
    as.numeric(logLik(fit_stageIII_sqrt))
  )
)

stageIII_model_comparison$Delta_AIC <-
  stageIII_model_comparison$AIC -
  min(stageIII_model_comparison$AIC)

stageIII_model_comparison


# ============================================================
# Step 12C-4. Stage III vs Stage I time-specific HR
# ============================================================

beta_stageIII <- coef(fit_stageIII_sqrt)["Stage_III"]
gamma_stageIII <- coef(fit_stageIII_sqrt)["tt(Stage_III)"]

vcov_stageIII <- vcov(fit_stageIII_sqrt)

var_beta <- vcov_stageIII["Stage_III", "Stage_III"]
var_gamma <- vcov_stageIII["tt(Stage_III)", "tt(Stage_III)"]
cov_bg <- vcov_stageIII["Stage_III", "tt(Stage_III)"]

time_points <- c(
  365,
  730,
  1095,
  1460,
  1825
)

stageIII_HR_time <- data.frame(
  Time_days = time_points
)

stageIII_HR_time$Time_years <-
  stageIII_HR_time$Time_days / 365

stageIII_HR_time$f_time <-
  sqrt(stageIII_HR_time$Time_days / 1000)

stageIII_HR_time$Beta_t <-
  beta_stageIII +
  gamma_stageIII * stageIII_HR_time$f_time

stageIII_HR_time$HR <-
  exp(stageIII_HR_time$Beta_t)

stageIII_HR_time$SE_Beta_t <-
  sqrt(
    var_beta +
      stageIII_HR_time$f_time^2 * var_gamma +
      2 * stageIII_HR_time$f_time * cov_bg
  )

stageIII_HR_time$Lower95 <-
  exp(
    stageIII_HR_time$Beta_t -
      1.96 * stageIII_HR_time$SE_Beta_t
  )

stageIII_HR_time$Upper95 <-
  exp(
    stageIII_HR_time$Beta_t +
      1.96 * stageIII_HR_time$SE_Beta_t
  )

stageIII_HR_time

# ============================================================
# Step 12C-5. Point-estimate HR=1 crossing
# ============================================================

t_cross_stageIII <-
  1000 *
  (beta_stageIII / (-gamma_stageIII))^2

t_cross_stageIII

t_cross_stageIII / 365

# ============================================================
# Step 12C-6. Stage III sqrt model stability
# ============================================================

stageIII_vcov <- vcov(fit_stageIII_sqrt)

stageIII_eigen <- eigen(
  stageIII_vcov,
  symmetric = TRUE
)$values

stageIII_eigen

kappa(stageIII_vcov)

stageIII_cor <- cov2cor(stageIII_vcov)

stageIII_cor[
  c("Stage_III", "tt(Stage_III)"),
  c("Stage_III", "tt(Stage_III)")
]

# ============================================================
# Step 12D. Incremental contribution of Stage III time effect
# ============================================================

anova(
  fit_gender_sqrt,
  fit_stageIII_sqrt,
  test = "LRT"
)

AIC(
  fit_gender_sqrt,
  fit_stageIII_sqrt
)

logLik(fit_gender_sqrt)
logLik(fit_stageIII_sqrt)


###Step 13：Final integrated time-dependent Cox model validation/reporting。
# ============================================================
# Step 13A. Final integrated time-dependent Cox model
# Corrected version
# ============================================================

final_td_model <- fit_stageIII_sqrt

# ------------------------------------------------------------
# 13A-1. Display final model
# ------------------------------------------------------------

final_td_model

# ------------------------------------------------------------
# 13A-2. Sample size
# ------------------------------------------------------------

n_final <- nrow(
  model.frame(final_td_model)
)

# ------------------------------------------------------------
# 13A-3. Number of events
# Use the original analysis dataset
# ------------------------------------------------------------

events_final <- sum(
  stage_complete_data$OS_event == 1,
  na.rm = TRUE
)

non_events_final <- sum(
  stage_complete_data$OS_event == 0,
  na.rm = TRUE
)

cat(
  "\nFinal model N =", n_final,
  "\nFinal model events =", events_final,
  "\nFinal model non-events =", non_events_final,
  "\n"
)

# ------------------------------------------------------------
# 13A-4. Final model formula
# ------------------------------------------------------------

formula(final_td_model)

# ============================================================
# 13A-5. Verify model data consistency
# ============================================================

model_data <- model.frame(final_td_model)

cat(
  "Rows in model.frame:",
  nrow(model_data),
  "\n"
)

cat(
  "Rows in stage_complete_data:",
  nrow(stage_complete_data),
  "\n"
)

cat(
  "Rows match:",
  nrow(model_data) == nrow(stage_complete_data),
  "\n"
)

# ============================================================
# 13A-6. Verify survival outcome
# ============================================================

table(stage_complete_data$OS_event)

summary(stage_complete_data$OS_time)

cat(
  "Events:",
  sum(stage_complete_data$OS_event == 1),
  "\n"
)

cat(
  "Censored:",
  sum(stage_complete_data$OS_event == 0),
  "\n"
)

# ============================================================
# Step 13B. Final model coefficient summary
# ============================================================

final_summary <- summary(final_td_model)

final_summary

# Coefficient table
final_coef_table <- as.data.frame(
  final_summary$coefficients
)

final_coef_table

# 95% CI
final_ci <- as.data.frame(
  confint(final_td_model)
)

final_ci

# ============================================================
# Step 13C. Global model fit statistics
# ============================================================

final_logLik <- as.numeric(
  logLik(final_td_model)
)

final_AIC <- AIC(final_td_model)

final_LRT <- final_summary$logtest["test"]
final_LRT_df <- final_summary$logtest["df"]
final_LRT_P <- final_summary$logtest["pvalue"]

final_Wald <- final_summary$waldtest["test"]
final_Wald_df <- final_summary$waldtest["df"]
final_Wald_P <- final_summary$waldtest["pvalue"]

final_Score <- final_summary$sctest["test"]
final_Score_df <- final_summary$sctest["df"]
final_Score_P <- final_summary$sctest["pvalue"]

final_concordance <- final_summary$concordance[1]
final_concordance_se <- final_summary$concordance[2]

cat("\n================ FINAL MODEL FIT ================\n")

cat("N =", n_final, "\n")
cat("Events =", events_final, "\n")
cat("LogLik =", final_logLik, "\n")
cat("AIC =", final_AIC, "\n")

cat(
  "LRT =",
  final_LRT,
  "df =",
  final_LRT_df,
  "P =",
  final_LRT_P,
  "\n"
)

cat(
  "Wald =",
  final_Wald,
  "df =",
  final_Wald_df,
  "P =",
  final_Wald_P,
  "\n"
)

cat(
  "Score =",
  final_Score,
  "df =",
  final_Score_df,
  "P =",
  final_Score_P,
  "\n"
)

cat(
  "Concordance =",
  final_concordance,
  "SE =",
  final_concordance_se,
  "\n"
)


# ============================================================
# Step 13D. ME_brown time-specific HR
# Final model: 1-SD increase
# ============================================================

time_points <- c(365, 730, 1095, 1460, 1825)

time_df <- data.frame(
  Time_days = time_points,
  Time_years = time_points / 365,
  f_time = sqrt(time_points / 1000)
)

# ------------------------------------------------------------
# ME_brown coefficients
# ------------------------------------------------------------

beta_ME <- coef(final_td_model)["ME_brown_SD"]

gamma_ME <- coef(
  final_td_model
)["tt(ME_brown_SD)"]

vc <- vcov(final_td_model)

var_beta_ME <- vc[
  "ME_brown_SD",
  "ME_brown_SD"
]

var_gamma_ME <- vc[
  "tt(ME_brown_SD)",
  "tt(ME_brown_SD)"
]

cov_ME <- vc[
  "ME_brown_SD",
  "tt(ME_brown_SD)"
]

# ------------------------------------------------------------
# Time-specific HR
# ------------------------------------------------------------

ME_HR_final <- time_df

ME_HR_final$Beta_t <-
  beta_ME +
  gamma_ME * ME_HR_final$f_time

ME_HR_final$HR <-
  exp(ME_HR_final$Beta_t)

ME_HR_final$SE_Beta_t <-
  sqrt(
    var_beta_ME +
      ME_HR_final$f_time^2 * var_gamma_ME +
      2 * ME_HR_final$f_time * cov_ME
  )

ME_HR_final$Lower95 <-
  exp(
    ME_HR_final$Beta_t -
      1.96 * ME_HR_final$SE_Beta_t
  )

ME_HR_final$Upper95 <-
  exp(
    ME_HR_final$Beta_t +
      1.96 * ME_HR_final$SE_Beta_t
  )

ME_HR_final

# ============================================================
# Step 13E
# Final integrated model:
# Time-specific HR for Gender_Male
# ============================================================

# ------------------------------------------------------------
# 13E.1 Define time points
# ------------------------------------------------------------

time_points <- c(365, 730, 1095, 1460, 1825)

time_df_gender <- data.frame(
  Time_days = time_points,
  Time_years = time_points / 365,
  f_time = sqrt(time_points / 1000)
)


# ------------------------------------------------------------
# 13E.2 Extract Gender coefficients and covariance matrix
# ------------------------------------------------------------

beta_gender <- coef(final_td_model)["Gender_Male"]

gamma_gender <- coef(final_td_model)["tt(Gender_Male)"]

vc_gender <- vcov(final_td_model)


# Variance of baseline Gender coefficient
var_beta_gender <-
  vc_gender["Gender_Male", "Gender_Male"]


# Variance of time-dependent Gender coefficient
var_gamma_gender <-
  vc_gender["tt(Gender_Male)", "tt(Gender_Male)"]


# Covariance between baseline and time-dependent coefficients
cov_gender <-
  vc_gender["Gender_Male", "tt(Gender_Male)"]


# ------------------------------------------------------------
# 13E.3 Calculate time-specific coefficient
# ------------------------------------------------------------

Gender_HR_final <- time_df_gender

Gender_HR_final$Beta_t <-
  beta_gender +
  gamma_gender * Gender_HR_final$f_time


# ------------------------------------------------------------
# 13E.4 Calculate time-specific HR
# ------------------------------------------------------------

Gender_HR_final$HR <-
  exp(Gender_HR_final$Beta_t)


# ------------------------------------------------------------
# 13E.5 Calculate SE of Beta(t)
# ------------------------------------------------------------

Gender_HR_final$SE_Beta_t <-
  sqrt(
    var_beta_gender +
      Gender_HR_final$f_time^2 * var_gamma_gender +
      2 * Gender_HR_final$f_time * cov_gender
  )


# ------------------------------------------------------------
# 13E.6 Calculate 95% CI
# ------------------------------------------------------------

Gender_HR_final$Lower95 <-
  exp(
    Gender_HR_final$Beta_t -
      1.96 * Gender_HR_final$SE_Beta_t
  )

Gender_HR_final$Upper95 <-
  exp(
    Gender_HR_final$Beta_t +
      1.96 * Gender_HR_final$SE_Beta_t
  )


# ------------------------------------------------------------
# 13E.7 Display final results
# ------------------------------------------------------------

Gender_HR_final


Gender_HR_final_report <-
  Gender_HR_final[, c(
    "Time_years",
    "HR",
    "Lower95",
    "Upper95"
  )]

Gender_HR_final_report

#If you wish to round the values directly for use in a formal paper or report:
Gender_HR_final_report_rounded <-
  transform(
    Gender_HR_final_report,
    HR = round(HR, 3),
    Lower95 = round(Lower95, 3),
    Upper95 = round(Upper95, 3)
  )

Gender_HR_final_report_rounded

# ============================================================
# Step 13F
# Final integrated model:
# Time-specific HR for Stage III vs Stage I
# ============================================================


# ------------------------------------------------------------
# 13F.1 Define time points
# ------------------------------------------------------------

time_points <- c(365, 730, 1095, 1460, 1825)

time_df_stageIII <- data.frame(
  Time_days = time_points,
  Time_years = time_points / 365,
  f_time = sqrt(time_points / 1000)
)


# ------------------------------------------------------------
# 13F.2 Extract Stage III coefficients
# ------------------------------------------------------------

beta_stageIII <-
  coef(final_td_model)["Stage_III"]

gamma_stageIII <-
  coef(final_td_model)["tt(Stage_III)"]


# ------------------------------------------------------------
# 13F.3 Extract covariance matrix
# ------------------------------------------------------------

vc_stageIII <- vcov(final_td_model)


# Variance of baseline Stage III coefficient
var_beta_stageIII <-
  vc_stageIII[
    "Stage_III",
    "Stage_III"
  ]


# Variance of time-dependent Stage III coefficient
var_gamma_stageIII <-
  vc_stageIII[
    "tt(Stage_III)",
    "tt(Stage_III)"
  ]


# Covariance between baseline and
# time-dependent Stage III coefficients
cov_stageIII <-
  vc_stageIII[
    "Stage_III",
    "tt(Stage_III)"
  ]


# ------------------------------------------------------------
# 13F.4 Calculate time-specific coefficient
# ------------------------------------------------------------

StageIII_HR_final <-
  time_df_stageIII

StageIII_HR_final$Beta_t <-
  beta_stageIII +
  gamma_stageIII *
  StageIII_HR_final$f_time


# ------------------------------------------------------------
# 13F.5 Calculate time-specific HR
# ------------------------------------------------------------

StageIII_HR_final$HR <-
  exp(
    StageIII_HR_final$Beta_t
  )


# ------------------------------------------------------------
# 13F.6 Calculate SE of Beta(t)
# ------------------------------------------------------------

StageIII_HR_final$SE_Beta_t <-
  sqrt(
    var_beta_stageIII +
      StageIII_HR_final$f_time^2 *
      var_gamma_stageIII +
      2 *
      StageIII_HR_final$f_time *
      cov_stageIII
  )


# ------------------------------------------------------------
# 13F.7 Calculate 95% CI
# ------------------------------------------------------------

StageIII_HR_final$Lower95 <-
  exp(
    StageIII_HR_final$Beta_t -
      1.96 *
      StageIII_HR_final$SE_Beta_t
  )

StageIII_HR_final$Upper95 <-
  exp(
    StageIII_HR_final$Beta_t +
      1.96 *
      StageIII_HR_final$SE_Beta_t
  )


# ------------------------------------------------------------
# 13F.8 Display complete results
# ------------------------------------------------------------

StageIII_HR_final

# ============================================================
# Step 13F.9
# Compact reporting table
# ============================================================

StageIII_HR_final_report <-
  StageIII_HR_final[, c(
    "Time_years",
    "HR",
    "Lower95",
    "Upper95"
  )]

StageIII_HR_final_report


# ============================================================
# Step 13F.10
# Rounded reporting table
# ============================================================

StageIII_HR_final_report_rounded <-
  transform(
    StageIII_HR_final_report,
    HR = round(HR, 3),
    Lower95 = round(Lower95, 3),
    Upper95 = round(Upper95, 3)
  )

StageIII_HR_final_report_rounded

###Proceed to Step 13G to systematically summarize all the results of the final integrated model.
Surv(OS_time, OS_event) ~
  ME_brown_SD + tt(ME_brown_SD) +
  Age_trait +
  Gender_Male + tt(Gender_Male) +
  Stage_II + Stage_III + tt(Stage_III) + Stage_IV



###Step 13A–13G （frozen checkpoint）。
save.image(
  file = "result_dir/Step13G_FINAL_checkpoint.RData"
)

#If you wish to place it in a dedicated directory：

dir.create(
  "Step13G_FINAL_checkpoint",
  showWarnings = FALSE
)


#1. Save the final Cox model
saveRDS(
  final_td_model,
  file = "result_dir/Step13G_FINAL_checkpoint/final_td_model.rds"
)

#Can be restored independently later:

final_td_model <- readRDS(
  "result_dir/Step13G_FINAL_checkpoint/final_td_model.rds"
)
#2. Save the final, complete data
saveRDS(
  stage_complete_data,
  file = "result_dir/Step13G_FINAL_checkpoint/stage_complete_data.rds"
)
#3. Save core WGCNA objects
saveRDS(
  MEs,
  file = "result_dir/Step13G_FINAL_checkpoint/MEs_tumor_only.rds"
)

saveRDS(
  datExpr,
  file = "result_dir/Step13G_FINAL_checkpoint/datExpr_tumor_only.rds"
)

saveRDS(
  moduleColors,
  file = "result_dir/Step13G_FINAL_checkpoint/moduleColors_tumor_only.rds"
)



