# ============================================================
# Script: GitHub2_Paper1_tumor_only_WGCNA_3.R
# Analysis module: ME_brown and Gender time-dependent Cox analysis
#
# Study:
# Integrated transcriptomic and single-cell characterization
# of AOC3, MS4A7, and IL7R in LUAD
#
# 1. Research objective
# To further investigate the proportional hazards (PH)
# assumption of the ME_brown module and Gender in the
# overall survival analysis, and to evaluate candidate
# time-dependent Cox regression models.
#
# 2. Data sources
# TCGA-LUAD primary tumor cohort.
# This script resumes the analysis from the previously
# generated ME_brown time-dependent Cox checkpoint.
#
# 3. Input files
# data/Paper1_ME_brown_checkpoint_Step10G.RData
#
# 4. Sample inclusion criteria
# Uses the complete-case Cox analysis data stored in the
# upstream checkpoint.
# The final sample size and event count are determined
# from the fitted model objects and should be verified
# after execution.
#
# 5. Analysis methods
# - Schoenfeld residual-based PH diagnostics
# - Gender variable coding and verification
# - Linear time-dependent Cox regression
# - log1p(time)-dependent Cox regression
# - Square-root time-dependent Cox regression
# - AIC and log-likelihood model comparisons
# - Time-specific HR estimation with 95% confidence intervals
# - Coefficient and variance-covariance stability diagnostics
#
# 6. Software environment
# R version: To be recorded from the actual runtime
# Key packages: here, WGCNA, survival, survminer
# Package versions: To be recorded from the actual runtime
#
# 7. Random seed
# No random seed is explicitly set in this script.
# The upstream checkpoint is treated as a fixed input.
#
# 8. Output files
# - Paper1_ME_brown_Gender_checkpoint_Final.RData
# - Paper1_Gender_HR_Male_vs_Female_1to5y_Final.csv
#
# 9. Key analytical results
# Numerical results are not hard-coded in this script.
# They should be reported from the fitted model summaries
# and generated output tables after execution.
#
# 10. Reproducibility notes
# - Requires the upstream ME_brown checkpoint.
# - Uses project-relative paths through the here package.
# - Cox models use Efron handling of tied event times.
# - Candidate time transformations are specified in the script.
# - The final Gender HR estimates are based on the
#   square-root time-dependent Cox model.
# - This script does not independently reconstruct the
#   upstream expression matrix, WGCNA modules, or clinical cohort.
# ============================================================


# ============================================================
# 00. Project setup
# ============================================================

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Please install the 'here' package.")
}

project_root <- here::here()

input_dir  <- file.path(project_root, "data")
result_dir <- file.path(project_root, "results")
figure_dir <- file.path(project_root, "figures")

dir.create(
  result_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  figure_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# Step 1. Resume the ME_brown time-dependent Cox analysis
# ============================================================

checkpoint_file <- file.path(
  input_dir,
  "Paper1_ME_brown_checkpoint_Step10G.RData"
)

if (!file.exists(checkpoint_file)) {
  stop(
    "Required checkpoint file was not found: ",
    checkpoint_file
  )
}

load(checkpoint_file)

# Inspect loaded objects
ls()


# ============================================================
# Required packages
# ============================================================

required_packages <- c(
  "WGCNA",
  "survival",
  "survminer"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(WGCNA)
library(survival)
library(survminer)


# ============================================================
# Step 11A. Gender_trait PH assumption
# ============================================================

gender_zph <- cox.zph(
  brown_cox_model,
  transform = "km"
)

print(
  gender_zph,
  digits = 4
)

gender_zph$table


# ============================================================
# Step 11B. Schoenfeld residual plot for Gender_trait
# ============================================================

plot(
  gender_zph,
  var = "Gender_trait",
  resid = TRUE,
  se = TRUE,
  df = 4,
  main = "Schoenfeld residuals for Gender_trait",
  xlab = "Time",
  ylab = "Scaled Schoenfeld residuals"
)

abline(
  h = 0,
  lty = 2
)


# ============================================================
# Step 11C. Check Gender_trait coding
# ============================================================

class(cox_design_complete$Gender_trait)

levels(cox_design_complete$Gender_trait)

table(
  cox_design_complete$Gender_trait,
  useNA = "ifany"
)

summary(
  cox_design_complete$Gender_trait
)


# ============================================================
# Step 11D-1. Create binary indicator for Gender
# ============================================================

brown_sd_model_data$Gender_Male <- ifelse(
  brown_sd_model_data$Gender_trait == "Male",
  1,
  0
)

table(
  brown_sd_model_data$Gender_trait,
  brown_sd_model_data$Gender_Male,
  useNA = "ifany"
)


# ============================================================
# Step 11D-2. Linear time-dependent Gender model
# ============================================================

fit_gender_linear <- coxph(
  Surv(OS_time, OS_event) ~
    ME_brown_SD +
    tt(ME_brown_SD) +
    Age_trait +
    Gender_Male +
    tt(Gender_Male) +
    Stage_trait,
  data = brown_sd_model_data,
  ties = "efron",
  tt = function(x, t, ...) {
    x * (t / 1000)
  }
)

summary(
  fit_gender_linear
)


# ============================================================
# Step 11E. Log-time-dependent Gender model
# ============================================================

fit_gender_log <- coxph(
  Surv(OS_time, OS_event) ~
    ME_brown_SD +
    tt(ME_brown_SD) +
    Age_trait +
    Gender_Male +
    tt(Gender_Male) +
    Stage_trait,
  data = brown_sd_model_data,
  ties = "efron",
  tt = function(x, t, ...) {
    x * log1p(t)
  }
)

summary(
  fit_gender_log
)


# ============================================================
# Step 11F. Square-root time-dependent Gender model
# ============================================================

fit_gender_sqrt <- coxph(
  Surv(OS_time, OS_event) ~
    ME_brown_SD +
    tt(ME_brown_SD) +
    Age_trait +
    Gender_Male +
    tt(Gender_Male) +
    Stage_trait,
  data = brown_sd_model_data,
  ties = "efron",
  tt = function(x, t, ...) {
    x * (sqrt(t) / sqrt(1000))
  }
)

summary(
  fit_gender_sqrt
)


# ============================================================
# Recheck Gender binary coding
# ============================================================

brown_sd_model_data$Gender_Male <- ifelse(
  brown_sd_model_data$Gender_trait == "Male",
  1,
  0
)

table(
  brown_sd_model_data$Gender_trait,
  brown_sd_model_data$Gender_Male,
  useNA = "ifany"
)


# ============================================================
# Step 11G. Compare Gender time-dependent specifications
# ============================================================

gender_model_comparison <- data.frame(
  Model = c(
    "Linear time-dependent",
    "log1p time-dependent",
    "sqrt time-dependent"
  ),
  
  AIC = c(
    AIC(fit_gender_linear),
    AIC(fit_gender_log),
    AIC(fit_gender_sqrt)
  ),
  
  LogLik = c(
    as.numeric(logLik(fit_gender_linear)),
    as.numeric(logLik(fit_gender_log)),
    as.numeric(logLik(fit_gender_sqrt))
  )
)

gender_model_comparison$Delta_AIC <-
  gender_model_comparison$AIC -
  min(gender_model_comparison$AIC)

gender_model_comparison[
  order(gender_model_comparison$AIC),
]


# ============================================================
# Step 11H. Sample-size consistency
# ============================================================

c(
  Linear_n = fit_gender_linear$n,
  Log_n = fit_gender_log$n,
  Sqrt_n = fit_gender_sqrt$n,
  
  Linear_events = fit_gender_linear$nevent,
  Log_events = fit_gender_log$nevent,
  Sqrt_events = fit_gender_sqrt$nevent
)


# ============================================================
# Step 11I. Time-specific HR for Gender (Male vs Female)
# Final candidate: sqrt-time-dependent model
# ============================================================

gender_coef <- coef(fit_gender_sqrt)
gender_vcov <- vcov(fit_gender_sqrt)

# Extract coefficients
b1 <- gender_coef["Gender_Male"]
b2 <- gender_coef["tt(Gender_Male)"]

# Extract variance-covariance components
v11 <- gender_vcov[
  "Gender_Male",
  "Gender_Male"
]

v22 <- gender_vcov[
  "tt(Gender_Male)",
  "tt(Gender_Male)"
]

v12 <- gender_vcov[
  "Gender_Male",
  "tt(Gender_Male)"
]

# Time points: 1–5 years
time_years <- 1:5
time_days <- time_years * 365

# sqrt-time transformation
f_time <- sqrt(time_days) / sqrt(1000)

# Time-dependent log-HR
beta_t <- b1 + b2 * f_time

# Standard error of beta(t)
se_beta_t <- sqrt(
  v11 +
    f_time^2 * v22 +
    2 * f_time * v12
)

# HR and 95% CI
HR_t <- exp(beta_t)

Lower95_t <- exp(
  beta_t - 1.96 * se_beta_t
)

Upper95_t <- exp(
  beta_t + 1.96 * se_beta_t
)

gender_HR_time <- data.frame(
  Time_years = time_years,
  Time_days = time_days,
  Beta = beta_t,
  SE = se_beta_t,
  HR_Male_vs_Female = HR_t,
  Lower95 = Lower95_t,
  Upper95 = Upper95_t
)

gender_HR_time


# ============================================================
# Step 11J. Stability diagnostics for final Gender model
# ============================================================

# 1. Coefficients
coef(fit_gender_sqrt)

# 2. Variance-covariance matrix
gender_vcov

# 3. Eigenvalues
eigen(
  gender_vcov,
  symmetric = TRUE
)$values

# 4. Condition number
kappa(gender_vcov)

# 5. Correlation between Gender baseline and time-dependent term
gender_cor <- cov2cor(gender_vcov)

gender_cor[
  "Gender_Male",
  "tt(Gender_Male)"
]


# ============================================================
# Save finalized Gender time-dependent Cox checkpoint
# ============================================================

gender_output_dir <- file.path(
  result_dir,
  "02_tumor_only_WGCNA",
  "Gender_time_dependent_Cox"
)

dir.create(
  gender_output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

save(
  cox_design_complete,
  brown_cox_model,
  
  # ME_brown finalized models
  fit_brown_linear,
  fit_brown_log,
  fit_brown_sqrt,
  brown_final_td,
  brown_final_td_SD,
  brown_HR_1SD,
  
  # Gender models
  fit_gender_linear,
  fit_gender_log,
  fit_gender_sqrt,
  gender_HR_time,
  gender_vcov,
  gender_cor,
  gender_model_comparison,
  
  # Gender binary variable
  brown_sd_model_data,
  
  file = file.path(
    gender_output_dir,
    "Paper1_ME_brown_Gender_checkpoint_Final.RData"
  )
)


# ============================================================
# Save the Gender HR table
# ============================================================

write.csv(
  gender_HR_time,
  file = file.path(
    gender_output_dir,
    "Paper1_Gender_HR_Male_vs_Female_1to5y_Final.csv"
  ),
  row.names = FALSE
)

