
# 00_setup.R

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Please install the 'here' package.")
}

project_root <- here::here()

input_dir <- file.path(project_root, "data")
result_dir <- file.path(project_root, "results")
figure_dir <- file.path(project_root, "figures")

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# Step 1：# Resume Workflow-04-2D WGCNA

# ============================================================

work_dir <- "input_dir/GitHub2_Paper1_tumor-only WGCNA"

setwd(work_dir)

getwd()

load(
  "input_dir/Workflow-04-2N_WGCNA_PH_multivariable_checkpoint_20260911.RData"
)

ls()

dim(datExpr_MAD)

dim(MEs_biological)

dim(clinical_WGCNA_traits)

dim(clinical_cox)

zph_grey60
zph_brown

# ============================================================
# Step 8C-2K
# ME_brown: assessment of time-dependent functional form
# ============================================================

# ------------------------------------------------------------
# 1. Check the existing PH test object
# ------------------------------------------------------------

zph_brown

# ------------------------------------------------------------
# 2. Extract Schoenfeld residuals for ME_brown
# ------------------------------------------------------------

brown_schoenfeld <- zph_brown$y[, "ME_brown"]

brown_time <- as.numeric(rownames(zph_brown$y))

# ------------------------------------------------------------
# 3. Basic QC
# ------------------------------------------------------------

length(brown_schoenfeld)
length(brown_time)

summary(brown_schoenfeld)
summary(brown_time)

sum(is.na(brown_schoenfeld))
sum(is.na(brown_time))


# ------------------------------------------------------------
# 4. Plot Schoenfeld residuals for ME_brown
# ------------------------------------------------------------

plot(
  zph_brown,
  var = "ME_brown",
  main = "Schoenfeld residuals: ME_brown"
)

#If you would like the graphics to be clearer, you can use:

plot(
  zph_brown,
  var = "ME_brown",
  main = "Schoenfeld residuals: ME_brown",
  xlab = "Time",
  ylab = "Scaled Schoenfeld residuals"
)

# ------------------------------------------------------------
# 5. Prepare diagnostic data
# ------------------------------------------------------------

brown_diag <- data.frame(
  time = brown_time,
  residual = brown_schoenfeld
)

brown_diag <- brown_diag[
  is.finite(brown_diag$time) &
    is.finite(brown_diag$residual) &
    brown_diag$time > 0,
]

# ------------------------------------------------------------
# 6. Candidate time transformations
# ------------------------------------------------------------

brown_diag$log_time <- log(brown_diag$time)
brown_diag$sqrt_time <- sqrt(brown_diag$time)
brown_diag$linear_time <- brown_diag$time
brown_diag$log1p_time <- log1p(brown_diag$time)

# ------------------------------------------------------------
# 7. Exploratory Spearman correlations
# ------------------------------------------------------------

brown_time_form_results <- data.frame(
  Time_function = character(),
  Spearman_rho = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

time_functions <- c(
  "log_time",
  "sqrt_time",
  "linear_time",
  "log1p_time"
)

for (tf in time_functions) {
  
  test <- cor.test(
    brown_diag$residual,
    brown_diag[[tf]],
    method = "spearman",
    exact = FALSE
  )
  
  brown_time_form_results <- rbind(
    brown_time_form_results,
    data.frame(
      Time_function = tf,
      Spearman_rho = unname(test$estimate),
      P_value = test$p.value
    )
  )
}

brown_time_form_results

#BH correction is merely supplementary information.

brown_time_form_results$FDR <- p.adjust(
  brown_time_form_results$P_value,
  method = "BH"
)

brown_time_form_results

# ------------------------------------------------------------
# 8. Visual comparison of candidate time functions
# ------------------------------------------------------------

par(mfrow = c(2, 2))

plot(
  brown_diag$log_time,
  brown_diag$residual,
  xlab = "log(Time)",
  ylab = "Scaled Schoenfeld residual",
  main = "ME_brown: log(Time)"
)

abline(
  lm(residual ~ log_time, data = brown_diag)
)

plot(
  brown_diag$sqrt_time,
  brown_diag$residual,
  xlab = "sqrt(Time)",
  ylab = "Scaled Schoenfeld residual",
  main = "ME_brown: sqrt(Time)"
)

abline(
  lm(residual ~ sqrt_time, data = brown_diag)
)

plot(
  brown_diag$linear_time,
  brown_diag$residual,
  xlab = "Time",
  ylab = "Scaled Schoenfeld residual",
  main = "ME_brown: Time"
)

abline(
  lm(residual ~ linear_time, data = brown_diag)
)

plot(
  brown_diag$log1p_time,
  brown_diag$residual,
  xlab = "log1p(Time)",
  ylab = "Scaled Schoenfeld residual",
  main = "ME_brown: log1p(Time)"
)

abline(
  lm(residual ~ log1p_time, data = brown_diag)
)

par(mfrow = c(1, 1))


# ============================================================
# Step 8C-2L
# ME_brown: candidate time-dependent Cox models
# ============================================================
library(survival)
packageVersion("survival")

exists("coxph")
exists("Surv")
# ------------------------------------------------------------
# 1. Check the primary model data
# ------------------------------------------------------------

dim(cox_design_complete)

table(cox_design_complete$OS_event)

summary(cox_design_complete$OS_time)

# ------------------------------------------------------------
# 2. Check required variables
# ------------------------------------------------------------

required_vars <- c(
  "OS_time",
  "OS_event",
  "ME_brown",
  "Age_trait",
  "Gender_trait",
  "Stage_trait"
)

all(required_vars %in% colnames(cox_design_complete))

# ------------------------------------------------------------
# 3. Make sure factor reference levels are correct
# ------------------------------------------------------------

cox_design_complete$Gender_trait <- factor(
  cox_design_complete$Gender_trait,
  levels = c("Female", "Male")
)

cox_design_complete$Stage_trait <- factor(
  cox_design_complete$Stage_trait,
  levels = c(
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IV"
  )
)

# ------------------------------------------------------------
# 4. Candidate 1: linear time-varying coefficient
# ------------------------------------------------------------

fit_brown_linear <- coxph(
  Surv(OS_time, OS_event) ~
    ME_brown +
    tt(ME_brown) +
    Age_trait +
    Gender_trait +
    Stage_trait,
  data = cox_design_complete,
  ties = "efron",
  tt = function(x, t, ...) {
    x * (t / 1000)
  }
)

# ------------------------------------------------------------
# 5. Candidate 2: log(time)-varying coefficient
# ------------------------------------------------------------

fit_brown_log <- coxph(
  Surv(OS_time, OS_event) ~
    ME_brown +
    tt(ME_brown) +
    Age_trait +
    Gender_trait +
    Stage_trait,
  data = cox_design_complete,
  ties = "efron",
  tt = function(x, t, ...) {
    x * log1p(t)
  }
)

# ------------------------------------------------------------
# 6. Candidate 3: sqrt(time)-varying coefficient
# ------------------------------------------------------------

fit_brown_sqrt <- coxph(
  Surv(OS_time, OS_event) ~
    ME_brown +
    tt(ME_brown) +
    Age_trait +
    Gender_trait +
    Stage_trait,
  data = cox_design_complete,
  ties = "efron",
  tt = function(x, t, ...) {
    x * (sqrt(t) / sqrt(1000))
  }
)

# ------------------------------------------------------------
# 7. Model summaries
# ------------------------------------------------------------

summary(fit_brown_linear)

summary(fit_brown_log)

summary(fit_brown_sqrt)


brown_time_model_comparison <- data.frame(
  Model = c(
    "Linear time",
    "log1p(time)",
    "sqrt(time)"
  ),
  
  AIC = c(
    AIC(fit_brown_linear),
    AIC(fit_brown_log),
    AIC(fit_brown_sqrt)
  ),
  
  LogLik = c(
    as.numeric(logLik(fit_brown_linear)),
    as.numeric(logLik(fit_brown_log)),
    as.numeric(logLik(fit_brown_sqrt))
  )
)

brown_time_model_comparison <- brown_time_model_comparison[
  order(brown_time_model_comparison$AIC),
]

brown_time_model_comparison$Delta_AIC <-
  brown_time_model_comparison$AIC -
  min(brown_time_model_comparison$AIC)

print(brown_time_model_comparison)

### It is also recommended to perform another very important comparison: standard Cox vs. time-dependent Cox.
brown_model_comparison_all <- data.frame(
  Model = c(
    "Conventional Cox",
    "Linear time-dependent",
    "log1p time-dependent",
    "sqrt time-dependent"
  ),
  
  AIC = c(
    AIC(primary_models$ME_brown),
    AIC(fit_brown_linear),
    AIC(fit_brown_log),
    AIC(fit_brown_sqrt)
  ),
  
  LogLik = c(
    as.numeric(logLik(primary_models$ME_brown)),
    as.numeric(logLik(fit_brown_linear)),
    as.numeric(logLik(fit_brown_log)),
    as.numeric(logLik(fit_brown_sqrt))
  )
)

brown_model_comparison_all$Delta_AIC <-
  brown_model_comparison_all$AIC -
  min(brown_model_comparison_all$AIC)

brown_model_comparison_all <-
  brown_model_comparison_all[
    order(brown_model_comparison_all$AIC),
  ]

print(brown_model_comparison_all)


# ============================================================
# Confirm conventional Cox model for ME_brown
# ============================================================

brown_cox_model <- primary_models$ME_brown

print(brown_cox_model)

# ============================================================
# Schoenfeld PH test for ME_brown
# ============================================================

brown_zph <- cox.zph(
  brown_cox_model,
  transform = "km"
)

print(
  brown_zph,
  digits = 4
)

# ============================================================
# Schoenfeld residual plot for ME_brown
# ============================================================

plot(
  brown_zph,
  var = "ME_brown",
  resid = TRUE,
  se = TRUE,
  df = 4,
  main = "Schoenfeld residuals for ME_brown",
  xlab = "Time",
  ylab = "Scaled Schoenfeld residuals"
)

abline(
  h = 0,
  lty = 2
)
rownames(brown_zph$table)


brown_model_comparison_all <- data.frame(
  
  Model = c(
    "Conventional Cox",
    "Linear time-dependent",
    "log1p time-dependent",
    "sqrt time-dependent"
  ),
  
  AIC = c(
    AIC(primary_models$ME_brown),
    AIC(fit_brown_linear),
    AIC(fit_brown_log),
    AIC(fit_brown_sqrt)
  ),
  
  LogLik = c(
    as.numeric(logLik(primary_models$ME_brown)),
    as.numeric(logLik(fit_brown_linear)),
    as.numeric(logLik(fit_brown_log)),
    as.numeric(logLik(fit_brown_sqrt))
  )
)

#Conventional Cox model:
primary_models$ME_brown

#Three time-dependent models:
fit_brown_linear
fit_brown_log
fit_brown_sqrt

# ============================================================
# Step 1. Compare the analysis samples
# Conventional Cox vs ME_brown sqrt(time)-dependent Cox
# ============================================================

# 1.1 Conventional Cox model
brown_cox_model <- primary_models$ME_brown

# 1.2 Final candidate: sqrt(time)-dependent model
brown_sqrt_model <- fit_brown_sqrt

# 1.3 Number of observations used by each model
nobs(brown_cox_model)
nobs(brown_sqrt_model)

# 1.4 Number of events
brown_cox_model$nevent
brown_sqrt_model$nevent

# ============================================================
# Step 1A. Confirm the actual number of observations
# ============================================================

# Conventional Cox
brown_cox_model$n
brown_cox_model$nevent

# Final sqrt(time)-dependent Cox
brown_final_td$n
brown_final_td$nevent

# ============================================================
# Step 2. Define the final ME_brown time-dependent model
# ============================================================

brown_final_td <- fit_brown_sqrt

summary(brown_final_td)


# ============================================================
# Step 3. Extract coefficients from the final model
# ============================================================

brown_td_coef <- summary(brown_final_td)$coefficients

brown_td_coef[
  c("ME_brown", "tt(ME_brown)"),
  ,
  drop = FALSE
]

# ============================================================
# Step 4. Time-specific HR for ME_brown
# ============================================================

b1 <- coef(brown_final_td)["ME_brown"]
b2 <- coef(brown_final_td)["tt(ME_brown)"]

time_points <- c(
  100,
  365,
  500,
  1000,
  1500,
  2000,
  3000
)

f_t <- sqrt(time_points / 1000)

beta_t <- b1 + b2 * f_t

brown_hr_time <- data.frame(
  Time_days = time_points,
  Beta_t = beta_t,
  HR = exp(beta_t)
)

brown_hr_time

# ============================================================
# Step 5. Time-specific HR and 95% CI
# ============================================================

V <- vcov(brown_final_td)

var_b1 <- V["ME_brown", "ME_brown"]
var_b2 <- V["tt(ME_brown)", "tt(ME_brown)"]
cov_b1b2 <- V["ME_brown", "tt(ME_brown)"]

se_beta_t <- sqrt(
  var_b1 +
    (f_t^2 * var_b2) +
    (2 * f_t * cov_b1b2)
)

brown_hr_time <- data.frame(
  Time_days = time_points,
  Beta_t = beta_t,
  SE = se_beta_t,
  HR = exp(beta_t),
  Lower95 = exp(beta_t - 1.96 * se_beta_t),
  Upper95 = exp(beta_t + 1.96 * se_beta_t)
)

brown_hr_time

# ============================================================
# Step 6. Estimate the time point at which HR(t) = 1
# ============================================================

t_hr1 <- 1000 * (
  b1 / (-b2)
)^2

t_hr1

# ============================================================
# Step 7. PH diagnostics after time-dependent correction
# ============================================================

brown_final_zph <- cox.zph(
  brown_final_td,
  transform = "km"
)

print(brown_final_zph, digits = 4)

###The error message at Step 7 actually indicates that a different verification method is required.

# ============================================================
# Step 8. Examine the observed follow-up time distribution
# ============================================================

summary(cox_design_complete$OS_time)

quantile(
  cox_design_complete$OS_time,
  probs = c(
    0,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.975,
    1
  ),
  na.rm = TRUE
)

# ============================================================
# Step 8B. Event-time distribution
# ============================================================

event_times <- cox_design_complete$OS_time[
  cox_design_complete$OS_event == 1
]

summary(event_times)

quantile(
  event_times,
  probs = c(
    0,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    1
  ),
  na.rm = TRUE
)

# ============================================================
# Step 8C. Number of events across follow-up intervals
# ============================================================

cut(
  event_times,
  breaks = c(
    0,
    365,
    730,
    1095,
    1460,
    1825,
    2190,
    2555,
    2920,
    Inf
  ),
  right = FALSE
) |> table()


# ============================================================
# Step 9. Time-specific HR at 1–5 years
# ============================================================

time_points_1to5y <- c(
  365,
  730,
  1095,
  1460,
  1825
)

f_t_1to5y <- sqrt(time_points_1to5y / 1000)

beta_t_1to5y <- b1 + b2 * f_t_1to5y

se_beta_t_1to5y <- sqrt(
  var_b1 +
    (f_t_1to5y^2 * var_b2) +
    (2 * f_t_1to5y * cov_b1b2)
)

brown_hr_1to5y <- data.frame(
  Time_years = 1:5,
  Time_days = time_points_1to5y,
  Beta_t = beta_t_1to5y,
  SE = se_beta_t_1to5y,
  HR = exp(beta_t_1to5y),
  Lower95 = exp(beta_t_1to5y - 1.96 * se_beta_t_1to5y),
  Upper95 = exp(beta_t_1to5y + 1.96 * se_beta_t_1to5y)
)

brown_hr_1to5y

# ============================================================
# Step 9B. Cumulative number of events by time point
# ============================================================

brown_hr_1to5y$Cumulative_events <- sapply(
  time_points_1to5y,
  function(x) sum(event_times <= x, na.rm = TRUE)
)

brown_hr_1to5y

# ============================================================
# Step 10A. Coefficient stability
# ============================================================

# Final model
brown_final_td <- fit_brown_sqrt

# 10A.1 Coefficients
brown_final_td$coefficients


# 10A.2 Full coefficient table
brown_coef_table <- coef(summary(brown_final_td))

brown_coef_table


# 10A.3 Create compact stability table
brown_coef_check <- data.frame(
  Term = names(coef(brown_final_td)),
  Estimate = coef(brown_final_td),
  SE = sqrt(diag(vcov(brown_final_td)))
)

brown_coef_check$Z <- with(
  brown_coef_check,
  Estimate / SE
)

brown_coef_check$Abs_Z <- abs(brown_coef_check$Z)

brown_coef_check

# ============================================================
# Step 10B. Variance-covariance matrix
# ============================================================

brown_vcov <- vcov(brown_final_td)

brown_vcov

# 10B.1 Symmetry
isSymmetric(brown_vcov)

# 10B.2 Eigenvalues
brown_vcov_eigen <- eigen(brown_vcov)

brown_vcov_eigen$values

# 10B.3 Numerical validity
anyNA(brown_vcov)
any(!is.finite(brown_vcov))

# 10B.4 Condition number
brown_vcov_condition <- kappa(brown_vcov)

brown_vcov_condition

# ============================================================
# Step 10B-5. Parameter correlation matrix
# ============================================================

brown_cor <- cov2cor(brown_vcov)

brown_cor

#Then specifically extract two core parameters:
brown_cor[
  c("ME_brown", "tt(ME_brown)"),
  c("ME_brown", "tt(ME_brown)")
]

# ============================================================
# Step 10C. Check time-varying beta(t)
# ============================================================

# Time points for numerical stability check
beta_check_times <- c(
  0,
  90,
  180,
  365,
  730,
  1095,
  1460,
  1825,
  2190,
  2555,
  2920,
  3650
)

# Extract coefficients
b1 <- coef(brown_final_td)["ME_brown"]
b2 <- coef(brown_final_td)["tt(ME_brown)"]

# Time transformation
f_beta <- sqrt(beta_check_times / 1000)

# Calculate beta(t)
beta_check <- data.frame(
  Time_days = beta_check_times,
  Time_years = beta_check_times / 365,
  Beta_t = as.numeric(b1 + b2 * f_beta)
)

# Calculate time-specific HR
beta_check$HR_t <- exp(beta_check$Beta_t)

beta_check

# ============================================================
# Step 10C-2. Monotonicity of beta(t)
# ============================================================

diff(beta_check$Beta_t)

all(diff(beta_check$Beta_t) <= 0)

# ============================================================
# Step 10C-3. 95% CI for beta(t)
# ============================================================

var_b1 <- brown_vcov["ME_brown", "ME_brown"]

var_b2 <- brown_vcov[
  "tt(ME_brown)",
  "tt(ME_brown)"
]

cov_b1b2 <- brown_vcov[
  "ME_brown",
  "tt(ME_brown)"
]

beta_check$SE_beta <- sqrt(
  var_b1 +
    (f_beta^2 * var_b2) +
    (2 * f_beta * cov_b1b2)
)

beta_check$Lower95_beta <- (
  beta_check$Beta_t -
    1.96 * beta_check$SE_beta
)

beta_check$Upper95_beta <- (
  beta_check$Beta_t +
    1.96 * beta_check$SE_beta
)

beta_check$Lower95_HR <- exp(
  beta_check$Lower95_beta
)

beta_check$Upper95_HR <- exp(
  beta_check$Upper95_beta
)

beta_check

# ============================================================
# Step 10D. ME_brown distribution and scale
# ============================================================

brown_mean <- mean(
  cox_design_complete$ME_brown,
  na.rm = TRUE
)

brown_sd <- sd(
  cox_design_complete$ME_brown,
  na.rm = TRUE
)

brown_mean
brown_sd

summary(
  cox_design_complete$ME_brown
)

quantile(
  cox_design_complete$ME_brown,
  probs = c(
    0,
    0.01,
    0.05,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99,
    1
  ),
  na.rm = TRUE
)

# Boxplot
boxplot(
  cox_design_complete$ME_brown,
  horizontal = TRUE,
  main = "Distribution of ME_brown",
  xlab = "ME_brown"
)


# ============================================================
# Step 10E. Reparameterize ME_brown per 1-SD increase
# ============================================================

brown_sd_model_data <- cox_design_complete

brown_sd_model_data$ME_brown_SD <- as.numeric(
  scale(
    brown_sd_model_data$ME_brown
  )
)

brown_final_td_SD <- coxph(
  Surv(OS_time, OS_event) ~
    ME_brown_SD +
    tt(ME_brown_SD) +
    Age_trait +
    Gender_trait +
    Stage_trait,
  data = brown_sd_model_data,
  ties = "efron",
  tt = function(x, t, ...) {
    x * (sqrt(t) / sqrt(1000))
  }
)

summary(brown_final_td_SD)

# ============================================================
# Step 10F. Time-specific HR per 1-SD increase
# ============================================================

# Extract coefficients
b1_SD <- coef(brown_final_td_SD)["ME_brown_SD"]
b2_SD <- coef(brown_final_td_SD)["tt(ME_brown_SD)"]

# Variance-covariance matrix
vcov_SD <- vcov(brown_final_td_SD)

var_b1_SD <- vcov_SD[
  "ME_brown_SD",
  "ME_brown_SD"
]

var_b2_SD <- vcov_SD[
  "tt(ME_brown_SD)",
  "tt(ME_brown_SD)"
]

cov_b1b2_SD <- vcov_SD[
  "ME_brown_SD",
  "tt(ME_brown_SD)"
]

# Time points: 1–5 years
time_days_SD <- c(
  365,
  730,
  1095,
  1460,
  1825
)

time_years_SD <- time_days_SD / 365

# sqrt(time / 1000)
f_SD <- sqrt(time_days_SD / 1000)

# Time-specific beta
beta_t_SD <- b1_SD + b2_SD * f_SD

# Standard error of beta(t)
SE_beta_SD <- sqrt(
  var_b1_SD +
    (f_SD^2 * var_b2_SD) +
    (2 * f_SD * cov_b1b2_SD)
)

# 95% CI on beta scale
Lower95_beta_SD <- beta_t_SD - 1.96 * SE_beta_SD
Upper95_beta_SD <- beta_t_SD + 1.96 * SE_beta_SD

# Convert to HR scale
HR_t_SD <- exp(beta_t_SD)

Lower95_HR_SD <- exp(Lower95_beta_SD)
Upper95_HR_SD <- exp(Upper95_beta_SD)

# Cumulative number of events
event_times <- cox_design_complete$OS_time[
  cox_design_complete$OS_event == 1
]

Cumulative_events_SD <- sapply(
  time_days_SD,
  function(x) sum(event_times <= x)
)

# Final table
brown_HR_1SD <- data.frame(
  Time_years = time_years_SD,
  Time_days = time_days_SD,
  Beta_t = beta_t_SD,
  SE = SE_beta_SD,
  HR_1SD = HR_t_SD,
  Lower95 = Lower95_HR_SD,
  Upper95 = Upper95_HR_SD,
  Cumulative_events = Cumulative_events_SD
)

brown_HR_1SD

# ============================================================
# Step 10G. Final equivalence check
# ============================================================

c(
  Original_AIC = AIC(brown_final_td),
  SD_AIC = AIC(brown_final_td_SD),
  
  Original_LogLik = as.numeric(logLik(brown_final_td)),
  SD_LogLik = as.numeric(logLik(brown_final_td_SD)),
  
  Original_n = brown_final_td$n,
  SD_n = brown_final_td_SD$n,
  
  Original_events = brown_final_td$nevent,
  SD_events = brown_final_td_SD$nevent
)

# Coefficient transformation check
c(
  Original_ME_brown =
    coef(brown_final_td)["ME_brown"],
  
  Expected_SD_ME_brown =
    coef(brown_final_td)["ME_brown"] * brown_sd,
  
  Observed_SD_ME_brown =
    coef(brown_final_td_SD)["ME_brown"],
  
  Original_tt =
    coef(brown_final_td)["tt(ME_brown)"],
  
  Expected_SD_tt =
    coef(brown_final_td)["tt(ME_brown)"] * brown_sd,
  
  Observed_SD_tt =
    coef(brown_final_td_SD)["tt(ME_brown_SD)"]
)


##I. What is recommended to be saved at this time?

# Raw data/design matrix
cox_design_complete

# ME_brown Standard Cox
brown_cox_model

# ME_brown time-dependent candidate models
fit_brown_linear
fit_brown_log
fit_brown_sqrt

# ME_brown Final Original-Scale Model
brown_final_td

# ME_brown Final 1-SD Reporting Model
brown_final_td_SD

# ME_brown 1-SD time-specific HR
brown_HR_1SD

# Statistical results related to ME_brown
brown_zph
brown_vcov
brown_cor
beta_check

#The most important of these is:

brown_final_td
brown_final_td_SD
brown_HR_1SD
cox_design_complete

##① R workspace checkpoint
save(
  cox_design_complete,
  brown_cox_model,
  fit_brown_linear,
  fit_brown_log,
  fit_brown_sqrt,
  brown_final_td,
  brown_final_td_SD,
  brown_HR_1SD,
  brown_zph,
  brown_vcov,
  brown_cor,
  beta_check,
  file = "Paper1_ME_brown_checkpoint_Step10G.RData"
)

#② Simultaneously save a "readable results table."
write.csv(
  brown_HR_1SD,
  "Paper1_ME_brown_HR_1SD_Final.csv",
  row.names = FALSE
)

