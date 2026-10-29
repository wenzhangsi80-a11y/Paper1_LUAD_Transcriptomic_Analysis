#Paper 1 – Survival visualization package
##GitHub3_Paper1_survival_analysis_3

## 00_setup.R

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Please install the 'here' package.")
}

project_root <- here::here()

input_dir <- file.path(project_root, "data")
result_dir <- file.path(project_root, "results")
figure_dir <- file.path(project_root, "figures")

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

work_dir <- input_dir 
setwd(work_dir)

getwd()

load(
  "input_dir/Paper1_Step6C5_adjusted_MS4A7_time_varying_Cox_checkpoint.RData"
)

ls()

# ============================================================
# Paper 1
# Survival analysis visualization package
#
# Figure 4
# Supplementary Figure/Table 
#
# Based on LOCKED Step 6B / 6C results
# Do NOT overwrite any locked checkpoint
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)


# ============================================================
# 1. Paths
# ============================================================

survival_dir <- "input_dir/Paper 1_Step 6-survival analysis"

checkpoint_6B5 <- file.path(
  survival_dir,
  "Paper1_Step6B5_continuous_Cox_PH_checkpoint.RData"
)

checkpoint_6C5 <- file.path(
  survival_dir,
  "Paper1_Step6C5_adjusted_MS4A7_time_varying_Cox_checkpoint.RData"
)

output_dir <- file.path(
  survival_dir,
  "Figure4D_TableSy_output"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ============================================================
# 2. Load LOCKED checkpoints
# ============================================================

load(checkpoint_6B5)
load(checkpoint_6C5)


# ============================================================
# 3. Basic object QC
# ============================================================

required_objects_6B5 <- c(
  "survival_gene_primary_z",
  "cox_univariable_models",
  "cox_univariable_ph"
)

missing_6B5 <- required_objects_6B5[
  !vapply(required_objects_6B5, exists, logical(1))
]

if (length(missing_6B5) > 0) {
  stop(
    "Missing required Step 6B5 objects: ",
    paste(missing_6B5, collapse = ", ")
  )
}


required_objects_6C <- c(
  "cox_modelA",
  "cox_modelB"
)

missing_6C <- required_objects_6C[
  !vapply(required_objects_6C, exists, logical(1))
]

if (length(missing_6C) > 0) {
  stop(
    "Missing required Step 6C objects: ",
    paste(missing_6C, collapse = ", ")
  )
}


# ============================================================
# 4. Locked cohort QC
# ============================================================

stopifnot(
  nrow(survival_gene_primary_z) == 485
)

stopifnot(
  sum(survival_gene_primary_z$OS_event) == 178
)

stopifnot(
  cox_modelA$n == 485
)

stopifnot(
  cox_modelA$nevent == 178
)

stopifnot(
  cox_modelB$n == 336
)

stopifnot(
  cox_modelB$nevent == 137
)


cat("\n============================================\n")
cat("LOCKED SURVIVAL COHORT QC PASSED\n")
cat("============================================\n")
cat("Primary cohort: N =", cox_modelA$n,
    " Events =", cox_modelA$nevent, "\n")
cat("TNM sensitivity: N =", cox_modelB$n,
    " Events =", cox_modelB$nevent, "\n")


# ============================================================
# 5. Extract univariable Cox results
# ============================================================

extract_cox_result <- function(model, gene, model_name, ph_p = NA_real_) {
  
  sm <- summary(model)
  
  coef_row <- sm$coefficients[1, ]
  ci_row   <- sm$conf.int[1, ]
  
  data.frame(
    Gene = gene,
    Model = model_name,
    N = model$n,
    Events = model$nevent,
    Beta = unname(coef_row["coef"]),
    HR = unname(ci_row["exp(coef)"]),
    CI_low = unname(ci_row["lower .95"]),
    CI_high = unname(ci_row["upper .95"]),
    P = unname(coef_row["Pr(>|z|)"]),
    PH_P = ph_p,
    stringsAsFactors = FALSE
  )
}


univariable_results <- bind_rows(
  extract_cox_result(
    cox_univariable_models$AOC3,
    "AOC3",
    "Univariable",
    cox_univariable_ph$AOC3$table["GLOBAL", "p"]
  ),
  extract_cox_result(
    cox_univariable_models$MS4A7,
    "MS4A7",
    "Univariable",
    cox_univariable_ph$MS4A7$table["GLOBAL", "p"]
  ),
  extract_cox_result(
    cox_univariable_models$IL7R,
    "IL7R",
    "Univariable",
    cox_univariable_ph$IL7R$table["GLOBAL", "p"]
  )
)


# ============================================================
# 6. Function to extract adjusted model results
# ============================================================

extract_gene_from_model <- function(model,
                                    gene,
                                    model_name,
                                    variable_pattern = paste0("^", gene, "_z$")) {
  
  sm <- summary(model)
  
  idx <- grep(
    variable_pattern,
    rownames(sm$coefficients)
  )
  
  if (length(idx) != 1) {
    
    # fallback: exact gene name
    idx <- which(
      rownames(sm$coefficients) == paste0(gene, "_z")
    )
  }
  
  if (length(idx) != 1) {
    stop(
      "Could not uniquely identify ",
      gene,
      " in model ",
      model_name
    )
  }
  
  coef_row <- sm$coefficients[idx, ]
  ci_row   <- sm$conf.int[idx, ]
  
  data.frame(
    Gene = gene,
    Model = model_name,
    N = model$n,
    Events = model$nevent,
    Beta = unname(coef_row["coef"]),
    HR = unname(ci_row["exp(coef)"]),
    CI_low = unname(ci_row["lower .95"]),
    CI_high = unname(ci_row["upper .95"]),
    P = unname(coef_row["Pr(>|z|)"]),
    PH_P = NA_real_,
    stringsAsFactors = FALSE
  )
}


# ============================================================
# 7. IMPORTANT:
#    Refit adjusted models for AOC3 / MS4A7 / IL7R
#
#    Cohorts are NOT reconstructed.
#    Clinical datasets are extracted from the locked
#    Model A / Model B complete-case objects.
# ============================================================


# ------------------------------------------------------------
# Model A complete-case design
# ------------------------------------------------------------

modelA_data <- model.frame(cox_modelA)

names(modelA_data)


# ------------------------------------------------------------
# Identify the response and clinical covariates
# ------------------------------------------------------------

# The locked Model A already contains the correct:
#   OS.time
#   OS_event
#   Age
#   Gender
#   Stage
#
# We only replace the gene predictor.
# ------------------------------------------------------------


# ------------------------------------------------------------
# Model B complete-case design
# ------------------------------------------------------------

modelB_data <- model.frame(cox_modelB)

names(modelB_data)


# ============================================================
# 8. Build gene-specific Model A
# ============================================================

fit_modelA_gene <- function(gene) {
  
  dat <- modelA_data
  
  if (!paste0(gene, "_z") %in% names(survival_gene_primary_z)) {
    
    stop(
      "Missing standardized expression variable: ",
      gene,
      "_z"
    )
  }
  
  # Match samples using Sample when available
  if ("Sample" %in% names(dat) &&
      "Sample" %in% names(survival_gene_primary_z)) {
    
    idx <- match(
      dat$Sample,
      survival_gene_primary_z$Sample
    )
    
    if (anyNA(idx)) {
      stop(
        "Sample matching failed for Model A / ",
        gene
      )
    }
    
    dat[[paste0(gene, "_z")]] <-
      survival_gene_primary_z[[paste0(gene, "_z")]][idx]
    
  } else {
    
    # Model frame preserves the locked complete-case order.
    if (nrow(dat) != 485) {
      stop(
        "Model A row number is not 485 and Sample matching is unavailable."
      )
    }
    
    dat[[paste0(gene, "_z")]] <-
      survival_gene_primary_z[[paste0(gene, "_z")]]
  }
  
  
  formula_gene <- as.formula(
    paste0(
      "Surv(OS.time, OS_event) ~ ",
      gene, "_z + Age + Gender + Stage"
    )
  )
  
  coxph(
    formula_gene,
    data = dat,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
}


# ============================================================
# 9. Build gene-specific Model B
# ============================================================

fit_modelB_gene <- function(gene) {
  
  dat <- modelB_data
  
  if ("Sample" %in% names(dat) &&
      "Sample" %in% names(survival_gene_primary_z)) {
    
    idx <- match(
      dat$Sample,
      survival_gene_primary_z$Sample
    )
    
    if (anyNA(idx)) {
      stop(
        "Sample matching failed for Model B / ",
        gene
      )
    }
    
    dat[[paste0(gene, "_z")]] <-
      survival_gene_primary_z[[paste0(gene, "_z")]][idx]
    
  } else {
    
    # Locked Model B complete-case order
    # is used only if Sample is unavailable.
    if (nrow(dat) != 336) {
      stop(
        "Model B row number is not 336 and Sample matching is unavailable."
      )
    }
    
    # We must identify the 336 Model-B observations
    # from the model frame itself.
    #
    # Expression is recovered using the model-frame row names
    # when available.
    
    if (!is.null(rownames(dat)) &&
        !is.null(rownames(survival_gene_primary_z))) {
      
      idx <- match(
        rownames(dat),
        rownames(survival_gene_primary_z)
      )
      
      if (!anyNA(idx)) {
        
        dat[[paste0(gene, "_z")]] <-
          survival_gene_primary_z[
            idx,
            paste0(gene, "_z")
          ]
        
      } else {
        
        stop(
          "Unable to safely recover Model B expression for ",
          gene
        )
      }
      
    } else {
      
      stop(
        "No safe sample-level mapping available for Model B / ",
        gene
      )
    }
  }
  
  
  formula_gene <- as.formula(
    paste0(
      "Surv(OS.time, OS_event) ~ ",
      gene, "_z + Age + Gender + T_stage + N_stage + M_stage"
    )
  )
  
  coxph(
    formula_gene,
    data = dat,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
}

# ============================================================
# Diagnostic check
# ============================================================

cat("\n============================================\n")
cat("CHECK survival_gene_primary_z\n")
cat("============================================\n\n")

cat("Dimensions:\n")
print(dim(survival_gene_primary_z))

cat("\nColumn names:\n")
print(colnames(survival_gene_primary_z))

cat("\nFirst rows:\n")
print(head(survival_gene_primary_z))

cat("\nStandardization parameters:\n")
print(z_parameters)

# ============================================================
# FIX: use the already-standardized expression columns
# in survival_gene_primary_z
#
# DO NOT recalculate Z-scores
# DO NOT modify the locked object
# ============================================================

genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

survival_gene_primary_z_downstream <-
  survival_gene_primary_z

# Create explicit _z aliases for downstream reporting only
for (g in genes) {
  
  survival_gene_primary_z_downstream[[paste0(g, "_z")]] <-
    survival_gene_primary_z_downstream[[g]]
}

# ============================================================
# QC: verify that aliases reproduce the locked Z-scores
# ============================================================

cat("\n============================================\n")
cat("DOWNSTREAM Z-SCORE QC\n")
cat("============================================\n")

for (g in genes) {
  
  original_z <-
    survival_gene_primary_z[[g]]
  
  alias_z <-
    survival_gene_primary_z_downstream[[paste0(g, "_z")]]
  
  cat("\n", g, "\n", sep = "")
  
  cat(
    "  N =",
    sum(!is.na(original_z)),
    "\n"
  )
  
  cat(
    "  Mean =",
    mean(original_z, na.rm = TRUE),
    "\n"
  )
  
  cat(
    "  SD =",
    sd(original_z, na.rm = TRUE),
    "\n"
  )
  
  cat(
    "  Alias identical =",
    isTRUE(
      all.equal(
        original_z,
        alias_z,
        check.attributes = FALSE
      )
    ),
    "\n"
  )
}


# ============================================================
# STEP 6D PRE-CHECK
# Verify locked multivariable Cox model structures
#
# NO MODEL REFITTING
# NO OBJECT MODIFICATION
# ============================================================

cat("\n============================================================\n")
cat("STEP 6D PRE-CHECK: LOCKED COX MODEL STRUCTURE\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. Locked standardized survival object
# ------------------------------------------------------------

cat("\n[1] survival_gene_primary_z\n")
cat("Dimensions: ",
    paste(dim(survival_gene_primary_z), collapse = " x "),
    "\n")

cat("Columns:\n")
print(colnames(survival_gene_primary_z))


# ------------------------------------------------------------
# 2. Model A
# ------------------------------------------------------------

cat("\n[2] cox_modelA\n")

cat("N = ", cox_modelA$n, "\n", sep = "")
cat("Events = ", cox_modelA$nevent, "\n", sep = "")

cat("\nModel formula:\n")
print(formula(cox_modelA))

cat("\nModel-frame columns:\n")
print(colnames(model.frame(cox_modelA)))

cat("\nFirst 3 rows of model frame:\n")
print(head(model.frame(cox_modelA), 3))


# ------------------------------------------------------------
# 3. Model B
# ------------------------------------------------------------

cat("\n[3] cox_modelB\n")

cat("N = ", cox_modelB$n, "\n", sep = "")
cat("Events = ", cox_modelB$nevent, "\n", sep = "")

cat("\nModel formula:\n")
print(formula(cox_modelB))

cat("\nModel-frame columns:\n")
print(colnames(model.frame(cox_modelB)))

cat("\nFirst 3 rows of model frame:\n")
print(head(model.frame(cox_modelB), 3))


# ------------------------------------------------------------
# 4. Important locked objects from Step 6C-5
# ------------------------------------------------------------

cat("\n[4] Step 6C-5 object availability\n")

check_objects <- c(
  "modelA_cc",
  "cox_modelA",
  "cox_modelB",
  "ms4a7_tv_adjusted_models",
  "ms4a7_adjusted_hr_time"
)

print(
  setNames(
    check_objects %in% ls(),
    check_objects
  )
)


# ------------------------------------------------------------
# 5. Final QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("EXPECTED LOCKED COHORTS\n")
cat("============================================================\n")

cat("Primary Model A: N = 485, Events = 178\n")
cat("TNM Model B:    N = 336, Events = 137\n")

cat("\nActual:\n")
cat("Model A: N = ", cox_modelA$n,
    ", Events = ", cox_modelA$nevent, "\n", sep = "")

cat("Model B: N = ", cox_modelB$n,
    ", Events = ", cox_modelB$nevent, "\n", sep = "")

##STEP 6D-1：Survival Reporting Table
# ============================================================
# PAPER 1
# STEP 6D-1: SURVIVAL REPORTING TABLE
#
# Purpose:
#   Generate manuscript-ready survival reporting tables
#   from the locked continuous-expression survival framework.
#
# IMPORTANT:
#   - Do NOT reconstruct the primary cohort.
#   - Do NOT recalculate Z-scores.
#   - Do NOT modify locked Step 6B5 / 6B6 / 6C5 objects.
#   - Primary cohort remains N = 485, Events = 178.
#   - TNM sensitivity cohort remains N = 336, Events = 137.
#   - BH-FDR is calculated separately within each model
#     across AOC3, MS4A7, and IL7R.
#
# Output:
#   1. Three-gene continuous Cox reporting table
#   2. MS4A7 time-varying Cox reporting table
#   3. MS4A7 time-specific HR table
#   4. Figure 4D forest-plot dataset
#   5. Excel workbook
#   6. CSV files
#   7. New Step 6D-1 checkpoint
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

required_packages <- c(
  "survival",
  "dplyr",
  "tibble",
  "openxlsx"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nPlease install them before running STEP 6D-1."
  )
}

library(survival)
library(dplyr)
library(tibble)
library(openxlsx)


# ============================================================
# 1. OUTPUT DIRECTORY
# ============================================================

output_dir <- file.path(
  getwd(),
  "STEP_6D1_survival_reporting"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ============================================================
# 2. REQUIRED OBJECT CHECK
# ============================================================

required_objects <- c(
  "survival_gene_primary_z",
  "survival_gene_primary_z_downstream",
  "z_parameters",
  "OS_event_locked",
  "cox_univariable_models",
  "cox_univariable_ph",
  "survival_cox_univariable_continuous",
  "ms4a7_tv_models",
  "ms4a7_tv_comparison",
  "best_time_function",
  "best_tv_model",
  "beta0",
  "beta1",
  "beta0_se",
  "beta1_se",
  "beta0_p",
  "beta1_p",
  "ms4a7_constant_model",
  "ms4a7_LRT",
  "ms4a7_AIC_comparison",
  "ms4a7_AIC_delta",
  "ms4a7_hr_time",
  "cox_modelA",
  "cox_modelB",
  "modelA_cc",
  "ms4a7_tv_adjusted_models",
  "ms4a7_tv_adjusted_comparison",
  "best_time_function_adjusted",
  "best_tv_adjusted_model",
  "ms4a7_adjusted_LRT",
  "ms4a7_adjusted_hr_time",
  "coef_adjusted",
  "vcov_adjusted",
  "summary_adjusted"
)

missing_objects <- required_objects[
  !sapply(required_objects, exists, inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    "The following required objects are missing:\n",
    paste(missing_objects, collapse = "\n")
  )
}


# ============================================================
# 3. LOCKED DATASET QC
# ============================================================

cat("\n============================================================\n")
cat("STEP 6D-1: LOCKED DATASET QC\n")
cat("============================================================\n")

genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

# Use the already-standardized downstream aliases.
reporting_data <- survival_gene_primary_z_downstream

# Confirm dimensions
stopifnot(
  nrow(reporting_data) == 485
)

# Confirm required columns
required_reporting_columns <- c(
  "Sample",
  "OS.time",
  "OS.status",
  "OS_event",
  "AOC3_z",
  "MS4A7_z",
  "IL7R_z"
)

missing_reporting_columns <- setdiff(
  required_reporting_columns,
  colnames(reporting_data)
)

if (length(missing_reporting_columns) > 0) {
  stop(
    "Missing reporting columns: ",
    paste(missing_reporting_columns, collapse = ", ")
  )
}

# Confirm complete primary cohort
for (g in paste0(genes, "_z")) {
  if (sum(!is.na(reporting_data[[g]])) != 485) {
    stop("Unexpected missing values in ", g)
  }
}

if (sum(reporting_data$OS_event == 1) != 178) {
  stop("Unexpected number of events in reporting_data.")
}

cat("Primary reporting cohort: N = ",
    nrow(reporting_data),
    ", Events = ",
    sum(reporting_data$OS_event == 1),
    "\n",
    sep = "")

cat("Z-score recalculation: NO\n")
cat("Locked cohort reconstruction: NO\n")


# ============================================================
# 4. VERIFY Z-SCORE ALIASES
# ============================================================

cat("\n============================================================\n")
cat("Z-SCORE ALIAS QC\n")
cat("============================================================\n")

for (g in genes) {
  
  original_z <- survival_gene_primary_z[[g]]
  alias_z <- reporting_data[[paste0(g, "_z")]]
  
  alias_ok <- isTRUE(
    all.equal(
      original_z,
      alias_z,
      check.attributes = FALSE
    )
  )
  
  cat(
    g,
    ": N = ",
    sum(!is.na(alias_z)),
    ", mean = ",
    signif(mean(alias_z), 6),
    ", SD = ",
    signif(sd(alias_z), 6),
    ", identical = ",
    alias_ok,
    "\n",
    sep = ""
  )
  
  if (!alias_ok) {
    stop("Z-score alias QC failed for ", g)
  }
}


# ============================================================
# 5. HELPER FUNCTION
#    EXTRACT CONSTANT COX RESULTS
# ============================================================

extract_cox_result <- function(
    fit,
    gene,
    model_name,
    ph_p = NA_real_
) {
  
  s <- summary(fit)
  
  coef_table <- as.data.frame(s$coefficients)
  conf_table <- as.data.frame(s$conf.int)
  
  # Find the gene coefficient robustly
  gene_row <- grep(
    paste0("^", gene, "_?z?$"),
    rownames(coef_table)
  )
  
  if (length(gene_row) != 1) {
    
    # Fallback: exact expected naming
    expected_name <- paste0(gene, "_z")
    
    if (expected_name %in% rownames(coef_table)) {
      gene_row <- which(
        rownames(coef_table) == expected_name
      )
    } else {
      stop(
        "Could not uniquely identify coefficient for ",
        gene,
        " in model ",
        model_name
      )
    }
  }
  
  tibble(
    Gene = gene,
    Model = model_name,
    N = fit$n,
    Events = fit$nevent,
    Beta = unname(coef_table[gene_row, "coef"]),
    HR = unname(conf_table[gene_row, "exp(coef)"]),
    CI_low = unname(conf_table[gene_row, "lower .95"]),
    CI_high = unname(conf_table[gene_row, "upper .95"]),
    P = unname(coef_table[gene_row, "Pr(>|z|)"]),
    PH_P = ph_p
  )
}


# ============================================================
# 6. UNIVARIABLE CONTINUOUS COX
#
# Use the LOCKED Step 6B5 models directly.
# No refitting.
# ============================================================

cat("\n============================================================\n")
cat("UNIVARIABLE CONTINUOUS COX\n")
cat("============================================================\n")

univariable_rows <- list()

for (g in genes) {
  
  fit <- cox_univariable_models[[g]]
  
  ph_value <- NA_real_
  
  if (!is.null(cox_univariable_ph[[g]])) {
    
    zph_obj <- cox_univariable_ph[[g]]
    
    zph_table <- as.data.frame(zph_obj$table)
    
    gene_row <- grep(
      paste0("^", g),
      rownames(zph_table)
    )
    
    if (length(gene_row) == 1) {
      ph_value <- zph_table[gene_row, "p"]
    }
  }
  
  univariable_rows[[g]] <-
    extract_cox_result(
      fit = fit,
      gene = g,
      model_name = "Univariable",
      ph_p = ph_value
    )
}

univariable_table <- bind_rows(
  univariable_rows
)


# ============================================================
# 7. TARGETED BH-FDR FOR UNIVARIABLE ANALYSIS
# ============================================================

univariable_table <- univariable_table %>%
  mutate(
    FDR = p.adjust(P, method = "BH")
  )


# ============================================================
# 8. VERIFY LOCKED UNIVARIABLE RESULTS
# ============================================================

cat("\n============================================================\n")
cat("LOCKED UNIVARIABLE RESULTS\n")
cat("============================================================\n")

print(
  univariable_table %>%
    mutate(
      HR = round(HR, 3),
      CI_low = round(CI_low, 3),
      CI_high = round(CI_high, 3),
      P = signif(P, 4),
      FDR = signif(FDR, 4),
      PH_P = signif(PH_P, 4)
    )
)


# ============================================================
# 9. MODEL A
#
# Primary adjusted model:
#   MS4A7_z + Age + Gender + Stage
#
# For complete three-gene reporting, the same locked
# Model-A cohort/coding is applied to all three genes.
#
# NOTE:
#   This does NOT alter the locked MS4A7 Model-A result.
# ============================================================

cat("\n============================================================\n")
cat("MODEL A: AGE + SEX + STAGE\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 9.1 Verify locked Model-A cohort
# ------------------------------------------------------------

if (cox_modelA$n != 485 ||
    cox_modelA$nevent != 178) {
  
  stop(
    "Locked Model A cohort mismatch. Expected N=485, Events=178."
  )
}

# ============================================================
# DIAGNOSTIC ONLY
# Inspect the exact values of locked pathologic_stage_2
# ============================================================

cat("\n============================================================\n")
cat("EXACT pathologic_stage_2 VALUES\n")
cat("============================================================\n")

print(
  sort(
    table(
      clinical_aligned$pathologic_stage_2,
      useNA = "ifany"
    )
  )
)

cat("\nUnique non-missing values:\n")

print(
  sort(
    unique(
      clinical_aligned$pathologic_stage_2[
        !is.na(clinical_aligned$pathologic_stage_2)
      ]
    )
  )
)

# ============================================================
# DIAGNOSTIC ONLY
# Compare Stage distribution in the LOCKED 485-sample cohort
# ============================================================

cat("\n============================================================\n")
cat("LOCKED 485-SAMPLE STAGE DISTRIBUTION\n")
cat("============================================================\n")

# Locked Model-A sample IDs
locked_modelA_samples <- rownames(
  model.frame(cox_modelA)
)

cat(
  "\nNumber of rownames in model.frame(cox_modelA): ",
  length(locked_modelA_samples),
  "\n",
  sep = ""
)

# IMPORTANT:
# model.frame() rownames should correspond to the model observations,
# but we will also compare against the formal locked sample vector
# if available.

if (exists("locked_primary_samples", inherits = TRUE)) {
  
  cat(
    "locked_primary_samples available: ",
    length(locked_primary_samples),
    "\n",
    sep = ""
  )
  
  stage_check_samples <- locked_primary_samples
  
} else {
  
  stage_check_samples <- locked_modelA_samples
  
}

stage_check <- clinical_aligned %>%
  filter(
    Sample %in% stage_check_samples
  ) %>%
  mutate(
    Stage_locked_raw = pathologic_stage_2
  )

cat(
  "\nMatched clinical rows: ",
  nrow(stage_check),
  "\n",
  sep = ""
)

cat("\nRaw pathologic_stage_2 distribution:\n")

print(
  table(
    stage_check$Stage_locked_raw,
    useNA = "ifany"
  )
)

cat("\nSample-count check:\n")

cat(
  "Expected locked cohort N = ",
  length(stage_check_samples),
  "\n",
  sep = ""
)

cat(
  "Matched clinical N = ",
  nrow(stage_check),
  "\n",
  sep = ""
)

# ============================================================
# 9.2 FINAL REPORTING-ONLY MODEL-A DATASET
#     Restrict to locked 485-sample cohort FIRST
# ============================================================

cat("\n============================================================\n")
cat("9.2 FINAL MODEL-A REPORTING DATASET\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 9.2.1 Verify locked sample vector
# ------------------------------------------------------------

if (!exists("locked_primary_samples", inherits = TRUE)) {
  stop(
    "locked_primary_samples is not available."
  )
}

if (length(locked_primary_samples) != 485) {
  stop(
    "Expected 485 locked primary samples, found ",
    length(locked_primary_samples),
    "."
  )
}

# ------------------------------------------------------------
# 9.2.2 Required raw clinical variables
# ------------------------------------------------------------

required_modelA_raw <- c(
  "Sample",
  "Age",
  "Gender",
  "pathologic_stage_2"
)

missing_modelA_raw <- setdiff(
  required_modelA_raw,
  colnames(clinical_aligned)
)

if (length(missing_modelA_raw) > 0) {
  stop(
    "Missing required clinical variables: ",
    paste(missing_modelA_raw, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 9.2.3 Restrict clinical data to the LOCKED 485 cohort
# ------------------------------------------------------------

clinical_modelA_locked <- clinical_aligned %>%
  filter(
    Sample %in% locked_primary_samples
  )

if (nrow(clinical_modelA_locked) != 485) {
  stop(
    "Clinical data matched ",
    nrow(clinical_modelA_locked),
    " samples instead of the locked 485."
  )
}

if (anyDuplicated(clinical_modelA_locked$Sample) > 0) {
  stop(
    "Duplicate Sample IDs detected in locked clinical cohort."
  )
}

# ------------------------------------------------------------
# 9.2.4 Preserve locked sample order
# ------------------------------------------------------------

clinical_modelA_locked <- clinical_modelA_locked[
  match(
    locked_primary_samples,
    clinical_modelA_locked$Sample
  ),
]

if (!identical(
  clinical_modelA_locked$Sample,
  locked_primary_samples
)) {
  stop(
    "Locked clinical sample order does not match locked_primary_samples."
  )
}

# ------------------------------------------------------------
# 9.2.5 Reconstruct Stage using the LOCKED Step 6C-2 rules
# ------------------------------------------------------------

clinical_modelA_vars <- clinical_modelA_locked %>%
  transmute(
    Sample = Sample,
    Age = as.numeric(Age),
    Gender = Gender,
    
    Stage = case_when(
      
      pathologic_stage_2 %in% c(
        "Stage I",
        "Stage IA",
        "Stage IB"
      ) ~ "Stage I",
      
      pathologic_stage_2 %in% c(
        "Stage II",
        "Stage IIA",
        "Stage IIB"
      ) ~ "Stage II",
      
      pathologic_stage_2 %in% c(
        "Stage IIIA",
        "Stage IIIB"
      ) ~ "Stage III",
      
      pathologic_stage_2 %in% c(
        "Stage IV"
      ) ~ "Stage IV",
      
      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------
# 9.2.6 Apply locked factor reference levels
# ------------------------------------------------------------

clinical_modelA_vars <- clinical_modelA_vars %>%
  mutate(
    
    Gender = factor(
      Gender,
      levels = c(
        "female",
        "male"
      )
    ),
    
    Stage = factor(
      Stage,
      levels = c(
        "Stage I",
        "Stage II",
        "Stage III",
        "Stage IV"
      )
    )
  )

# ------------------------------------------------------------
# 9.2.7 STRICT QC
# ------------------------------------------------------------

cat("\nAge missing:\n")
print(sum(is.na(clinical_modelA_vars$Age)))

cat("\nGender distribution:\n")
print(
  table(
    clinical_modelA_vars$Gender,
    useNA = "ifany"
  )
)

cat("\nStage distribution in LOCKED 485 cohort:\n")
print(
  table(
    clinical_modelA_vars$Stage,
    useNA = "ifany"
  )
)

expected_stage_counts <- c(
  "Stage I" = 263,
  "Stage II" = 117,
  "Stage III" = 79,
  "Stage IV" = 26
)

observed_stage_counts <- table(
  clinical_modelA_vars$Stage
)

if (!identical(
  as.integer(
    observed_stage_counts[
      names(expected_stage_counts)
    ]
  ),
  as.integer(
    expected_stage_counts
  )
)) {
  
  stop(
    "LOCKED Stage distribution mismatch."
  )
}

if (sum(is.na(clinical_modelA_vars$Stage)) != 0) {
  stop(
    "Unexpected missing Stage values in locked 485 cohort."
  )
}

if (sum(is.na(clinical_modelA_vars$Age)) != 0) {
  stop(
    "Unexpected missing Age values in locked 485 cohort."
  )
}

if (sum(is.na(clinical_modelA_vars$Gender)) != 0) {
  stop(
    "Unexpected missing Gender values in locked 485 cohort."
  )
}

cat(
  "\nLOCKED MODEL-A CLINICAL CODING QC: PASSED\n",
  "N = ", nrow(clinical_modelA_vars),
  "\n",
  sep = ""
)

# ============================================================
# 9.3 BUILD MODEL-A REPORTING DATASET
# ============================================================

modelA_reporting <- reporting_data %>%
  left_join(
    clinical_modelA_vars,
    by = "Sample"
  )

if (nrow(modelA_reporting) != 485) {
  stop(
    "Model-A reporting dataset does not contain 485 observations."
  )
}

if (anyDuplicated(modelA_reporting$Sample) > 0) {
  stop(
    "Duplicate Sample IDs detected."
  )
}

if (sum(modelA_reporting$OS_event) != 178) {
  stop(
    "Expected 178 events."
  )
}

if (any(is.na(modelA_reporting$Age))) {
  stop("Missing Age detected.")
}

if (any(is.na(modelA_reporting$Gender))) {
  stop("Missing Gender detected.")
}

if (any(is.na(modelA_reporting$Stage))) {
  stop("Missing Stage detected.")
}

cat(
  "\nMODEL-A REPORTING DATASET QC: PASSED\n",
  "N = ", nrow(modelA_reporting),
  "\nEvents = ", sum(modelA_reporting$OS_event),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9.4 Fit Model A for all three genes
#
# These fits reproduce the locked cohort/coding and are
# used only to construct a unified three-gene reporting table.
# ------------------------------------------------------------


modelA_fits <- list()

for (g in genes) {
  
  formula_A <- as.formula(
    paste(
      "Surv(OS.time, OS_event) ~",
      paste0(g, "_z"),
      "+ Age + Gender + Stage"
    )
  )
  
  modelA_fits[[g]] <- coxph(
    formula_A,
    data = modelA_reporting,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
}

# ------------------------------------------------------------
# 9.5 Extract Model-A results
# ------------------------------------------------------------

modelA_table <- bind_rows(
  lapply(
    genes,
    function(g) {
      
      extract_cox_result(
        fit = modelA_fits[[g]],
        gene = g,
        model_name = "Adjusted: age + sex + stage",
        ph_p = NA_real_
      )
      
    }
  )
)


# ============================================================
# 10. MODEL-A TARGETED BH-FDR
# ============================================================

modelA_table <- modelA_table %>%
  mutate(
    FDR = p.adjust(P, method = "BH")
  )

modelA_table 



# ============================================================
# 11. MODEL B
#
# Sensitivity model:
#   gene + Age + Gender + T + N + M
#
# Complete-case cohort must remain:
#   N = 336
#   Events = 137
# ============================================================

cat("\n============================================================\n")
cat("MODEL B: AGE + SEX + T/N/M\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 11.1 Verify locked Model-B cohort
# ------------------------------------------------------------

if (cox_modelB$n != 336 ||
    cox_modelB$nevent != 137) {
  
  stop(
    "Locked Model B cohort mismatch. Expected N=336, Events=137."
  )
}


# ------------------------------------------------------------
# 11.2 Prepare Model-B reporting dataset
# ------------------------------------------------------------

cat("\n--- Reconstructing TNM variables for reporting-only Model B ---\n")

# ------------------------------------------------------------
# 11.2.1 Check required raw clinical variables
# ------------------------------------------------------------

required_raw_TNM <- c(
  "Sample",
  "Age",
  "Gender",
  "pathologic_t",
  "pathologic_n",
  "pathologic_m"
)

missing_raw_TNM <- setdiff(
  required_raw_TNM,
  colnames(clinical_aligned)
)

if (length(missing_raw_TNM) > 0) {
  stop(
    paste0(
      "Missing required raw clinical variables: ",
      paste(missing_raw_TNM, collapse = ", ")
    )
  )
}

cat("Required raw TNM variables: PASSED\n")


# ------------------------------------------------------------
# 11.2.2 Restrict to the locked primary 485-sample cohort
# ------------------------------------------------------------

clinical_modelB_raw <- clinical_aligned %>%
  filter(Sample %in% locked_primary_samples) %>%
  mutate(
    .sample_order = match(Sample, locked_primary_samples)
  ) %>%
  arrange(.sample_order) %>%
  select(
    Sample,
    Age,
    Gender,
    pathologic_t,
    pathologic_n,
    pathologic_m
  )

if (nrow(clinical_modelB_raw) != 485) {
  stop(
    paste0(
      "Model-B raw clinical cohort mismatch. Expected 485, got ",
      nrow(clinical_modelB_raw)
    )
  )
}

if (anyDuplicated(clinical_modelB_raw$Sample) > 0) {
  stop("Duplicate Sample IDs detected in Model-B raw clinical cohort.")
}

if (!identical(
  clinical_modelB_raw$Sample,
  locked_primary_samples
)) {
  stop("Model-B raw clinical sample order does not match locked primary cohort.")
}

cat("Locked primary cohort alignment: PASSED\n")
cat("N =", nrow(clinical_modelB_raw), "\n")


# ------------------------------------------------------------
# 11.2.3 Reconstruct T stage
#
# Locked Step 6C-2 rules:
#   T1/T1a/T1b -> T1
#   T2/T2a/T2b -> T2
#   T3         -> T3
#   T4         -> T4
#   TX         -> NA
# ------------------------------------------------------------

clinical_modelB_raw <- clinical_modelB_raw %>%
  mutate(
    T_stage = case_when(
      pathologic_t %in% c("T1", "T1a", "T1b") ~ "T1",
      pathologic_t %in% c("T2", "T2a", "T2b") ~ "T2",
      pathologic_t == "T3" ~ "T3",
      pathologic_t == "T4" ~ "T4",
      pathologic_t == "TX" ~ NA_character_,
      TRUE ~ NA_character_
    )
  )


# ------------------------------------------------------------
# 11.2.4 Reconstruct N stage
#
# Locked Step 6C-2 rules:
#   N0 -> N0
#   N1 -> N1
#   N2 -> N2
#   N3 -> N3
#   NX -> NA
# ------------------------------------------------------------

clinical_modelB_raw <- clinical_modelB_raw %>%
  mutate(
    N_stage = case_when(
      pathologic_n == "N0" ~ "N0",
      pathologic_n == "N1" ~ "N1",
      pathologic_n == "N2" ~ "N2",
      pathologic_n == "N3" ~ "N3",
      pathologic_n == "NX" ~ NA_character_,
      TRUE ~ NA_character_
    )
  )


# ------------------------------------------------------------
# 11.2.5 Reconstruct M stage
#
# Locked Step 6C-2 rules:
#   M0/M1/M1a/M1b retained as M0/M1
#   MX/NA -> NA
# ------------------------------------------------------------

clinical_modelB_raw <- clinical_modelB_raw %>%
  mutate(
    M_stage = case_when(
      pathologic_m == "M0" ~ "M0",
      pathologic_m %in% c("M1", "M1a", "M1b") ~ "M1",
      pathologic_m %in% c("MX", NA_character_) ~ NA_character_,
      TRUE ~ NA_character_
    )
  )


# ------------------------------------------------------------
# 11.2.6 Convert to frozen factor levels
# ------------------------------------------------------------

clinical_modelB_raw <- clinical_modelB_raw %>%
  mutate(
    Gender = factor(
      Gender,
      levels = c("female", "male")
    ),
    
    T_stage = factor(
      T_stage,
      levels = c("T1", "T2", "T3", "T4")
    ),
    
    N_stage = factor(
      N_stage,
      levels = c("N0", "N1", "N2", "N3")
    ),
    
    M_stage = factor(
      M_stage,
      levels = c("M0", "M1")
    )
  )


# ------------------------------------------------------------
# 11.2.7 QC reconstructed TNM variables
# ------------------------------------------------------------

cat("\n--- Reconstructed TNM distribution in locked 485 cohort ---\n")

cat("\nT_stage:\n")
print(table(clinical_modelB_raw$T_stage, useNA = "ifany"))

cat("\nN_stage:\n")
print(table(clinical_modelB_raw$N_stage, useNA = "ifany"))

cat("\nM_stage:\n")
print(table(clinical_modelB_raw$M_stage, useNA = "ifany"))


# ------------------------------------------------------------
# 11.2.8 Construct complete-case Model-B reporting cohort
# ------------------------------------------------------------

modelB_reporting <- clinical_modelB_raw %>%
  select(
    Sample,
    Age,
    Gender,
    T_stage,
    N_stage,
    M_stage
  ) %>%
  left_join(
    reporting_data %>%
      select(
        Sample,
        OS.time,
        OS.status,
        OS_event,
        AOC3_z,
        MS4A7_z,
        IL7R_z
      ),
    by = "Sample"
  )


# ------------------------------------------------------------
# 11.2.9 Verify expression/survival alignment
# ------------------------------------------------------------

if (nrow(modelB_reporting) != 485) {
  stop(
    paste0(
      "Model-B reporting dataset mismatch. Expected 485 rows, got ",
      nrow(modelB_reporting)
    )
  )
}

if (anyDuplicated(modelB_reporting$Sample) > 0) {
  stop("Duplicate Sample IDs detected after Model-B join.")
}

if (anyNA(modelB_reporting$OS_event)) {
  stop("Missing OS_event detected in Model-B reporting dataset.")
}


# ------------------------------------------------------------
# 11.2.10 Define complete-case Model-B cohort
# ------------------------------------------------------------

modelB_complete <- modelB_reporting %>%
  filter(
    complete.cases(
      Age,
      Gender,
      T_stage,
      N_stage,
      M_stage,
      OS.time,
      OS_event
    )
  )


# ------------------------------------------------------------
# 11.2.11 Locked cohort QC
# ------------------------------------------------------------

modelB_N <- nrow(modelB_complete)
modelB_events <- sum(modelB_complete$OS_event == 1)

cat("\n============================================================\n")
cat("MODEL B COMPLETE-CASE QC\n")
cat("============================================================\n")

cat("N =", modelB_N, "\n")
cat("Events =", modelB_events, "\n")
cat("Censored =", modelB_N - modelB_events, "\n")
cat("Excluded =", 485 - modelB_N, "\n")


if (modelB_N != 336) {
  stop(
    paste0(
      "Model-B complete-case N mismatch. Expected 336, got ",
      modelB_N
    )
  )
}

if (modelB_events != 137) {
  stop(
    paste0(
      "Model-B event count mismatch. Expected 137, got ",
      modelB_events
    )
  )
}

cat("\nMODEL B REPORTING COHORT QC: PASSED\n")
cat("Locked complete-case cohort = N=336, Events=137\n")


# ------------------------------------------------------------
# 11.3 Apply the SAME complete-case rule
# ------------------------------------------------------------

modelB_required <- c(
  "OS.time",
  "OS_event",
  "Age",
  "Gender",
  "T_stage",
  "N_stage",
  "M_stage"
)

modelB_complete <- complete.cases(
  modelB_reporting[, modelB_required]
)

modelB_reporting_cc <- modelB_reporting[
  modelB_complete,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 11.4 Model-B cohort QC
# ------------------------------------------------------------

modelB_N <- nrow(modelB_reporting_cc)

modelB_events <- sum(
  modelB_reporting_cc$OS_event == 1
)

cat(
  "Reproduced Model-B complete-case cohort: N = ",
  modelB_N,
  ", Events = ",
  modelB_events,
  "\n",
  sep = ""
)

if (modelB_N != 336 ||
    modelB_events != 137) {
  
  stop(
    "Model-B complete-case reconstruction does not match ",
    "the locked cohort (N=336, Events=137)."
  )
}


# ------------------------------------------------------------
# 11.5 Fit Model B for all three genes
# ------------------------------------------------------------

modelB_fits <- list()

for (g in genes) {
  
  formula_B <- as.formula(
    paste(
      "Surv(OS.time, OS_event) ~",
      paste0(g, "_z"),
      "+ Age + Gender + T_stage + N_stage + M_stage"
    )
  )
  
  modelB_fits[[g]] <- coxph(
    formula_B,
    data = modelB_reporting_cc,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
}


# ------------------------------------------------------------
# 11.6 Extract Model-B results
# ------------------------------------------------------------

modelB_table <- bind_rows(
  lapply(
    genes,
    function(g) {
      
      extract_cox_result(
        fit = modelB_fits[[g]],
        gene = g,
        model_name = "Sensitivity: age + sex + T/N/M",
        ph_p = NA_real_
      )
      
    }
  )
)


# ============================================================
# 12. MODEL-B TARGETED BH-FDR
# ============================================================

modelB_table <- modelB_table %>%
  mutate(
    FDR = p.adjust(P, method = "BH")
  )


# ============================================================
# 13. COMBINE THREE-GENE CONSTANT-COEFFICIENT TABLE
# ============================================================

continuous_cox_table <- bind_rows(
  univariable_table,
  modelA_table,
  modelB_table
) %>%
  select(
    Gene,
    Model,
    N,
    Events,
    Beta,
    HR,
    CI_low,
    CI_high,
    P,
    FDR,
    PH_P
  )


# ============================================================
# 14. ADD DISPLAY COLUMNS
# ============================================================

continuous_cox_table_display <-
  continuous_cox_table %>%
  mutate(
    HR_CI = sprintf(
      "%.3f (%.3f–%.3f)",
      HR,
      CI_low,
      CI_high
    ),
    P_display = format.pval(
      P,
      digits = 3,
      eps = 0.001
    ),
    FDR_display = format.pval(
      FDR,
      digits = 3,
      eps = 0.001
    ),
    PH_P_display = ifelse(
      is.na(PH_P),
      NA_character_,
      format.pval(
        PH_P,
        digits = 3,
        eps = 0.001
      )
    )
  )


# ============================================================
# 15. MS4A7 TIME-VARYING COX
#
# Read locked results.
# Do NOT refit.
# ============================================================

cat("\n============================================================\n")
cat("MS4A7 TIME-VARYING COX\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 15.1 Unadjusted time-varying model
# ------------------------------------------------------------

ms4a7_tv_summary <- tibble(
  Analysis = "Unadjusted",
  Time_function = best_time_function,
  N = 485,
  Events = 178,
  Beta0 = beta0,
  SE_Beta0 = beta0_se,
  P_Beta0 = beta0_p,
  HR0 = exp(beta0),
  Beta1 = beta1,
  SE_Beta1 = beta1_se,
  P_Beta1 = beta1_p,
  LRT = as.numeric(ms4a7_LRT),
  LRT_df = 1,
  LRT_P = pchisq(
    as.numeric(ms4a7_LRT),
    df = 1,
    lower.tail = FALSE
  ),
  Delta_AIC = as.numeric(ms4a7_AIC_delta)
)


# ------------------------------------------------------------
# 15.2 Adjusted time-varying model
# ------------------------------------------------------------

adjusted_lrt_value <- as.numeric(
  ms4a7_adjusted_LRT
)

adjusted_lrt_p <- pchisq(
  adjusted_lrt_value,
  df = 1,
  lower.tail = FALSE
)

ms4a7_tv_adjusted_summary <- tibble(
  Analysis = "Adjusted: age + sex + stage",
  Time_function = best_time_function_adjusted,
  N = 485,
  Events = 178,
  Beta0 = beta0_adjusted_locked,
  SE_Beta0 = se0_adjusted_locked,
  P_Beta0 = p0_adjusted_locked,
  HR0 = HR0_adjusted_locked,
  Beta1 = beta1_adjusted_locked,
  SE_Beta1 = se1_adjusted_locked,
  P_Beta1 = p1_adjusted_locked,
  LRT = adjusted_lrt_value,
  LRT_df = 1,
  LRT_P = adjusted_lrt_p,
  Delta_AIC = 2.222337
)



# ============================================================
# 16. TIME-SPECIFIC HR TABLE
# ============================================================


extract_time_hr <- function(
    x,
    analysis_name
) {
  
  x <- as.data.frame(x)
  
  # Standardize column names
  names(x) <- make.names(names(x))
  
  # ----------------------------------------------------------
  # Locked time-specific HR object structure:
  #
  # Time_years
  # Time_days
  # HR_per_1SD
  # CI_low
  # CI_high
  # ----------------------------------------------------------
  
  required_columns <- c(
    "Time_years",
    "Time_days",
    "HR_per_1SD",
    "CI_low",
    "CI_high"
  )
  
  if (!all(required_columns %in% names(x))) {
    
    stop(
      "Unexpected structure of time-specific HR object for ",
      analysis_name,
      ".\nAvailable columns: ",
      paste(names(x), collapse = ", "),
      "\nExpected columns: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  x %>%
    as_tibble() %>%
    mutate(
      Analysis = analysis_name
    ) %>%
    select(
      Analysis,
      Time_years,
      Time_days,
      HR_per_1SD,
      CI_low,
      CI_high
    )
}

# ------------------------------------------------------------
# 16.1 Unadjusted HR
# ------------------------------------------------------------

ms4a7_time_hr_unadjusted <-
  extract_time_hr(
    ms4a7_hr_time,
    "Unadjusted"
  )

print(ms4a7_time_hr_unadjusted)

# ------------------------------------------------------------
# 16.2 Adjusted HR
# ------------------------------------------------------------

ms4a7_time_hr_adjusted <-
  extract_time_hr(
    ms4a7_adjusted_hr_time,
    "Adjusted: age + sex + stage"
  )


ms4a7_time_specific_hr <- bind_rows(
  ms4a7_time_hr_unadjusted,
  ms4a7_time_hr_adjusted
)


# ============================================================
# 17. VERIFY EXPECTED MS4A7 TIME-SPECIFIC RESULTS
# ============================================================

cat("\n============================================================\n")
cat("MS4A7 TIME-SPECIFIC HR\n")
cat("============================================================\n")

print(ms4a7_time_specific_hr)


# ============================================================
# 18. FIGURE 4D DATASET
#
# Constant-coefficient models only.
#
# MS4A7 time-varying analysis is reported separately.
# ============================================================

figure4D_data <- continuous_cox_table %>%
  filter(
    Model %in% c(
      "Univariable",
      "Adjusted: age + sex + stage",
      "Sensitivity: age + sex + T/N/M"
    )
  ) %>%
  mutate(
    Model = factor(
      Model,
      levels = c(
        "Univariable",
        "Adjusted: age + sex + stage",
        "Sensitivity: age + sex + T/N/M"
      )
    ),
    Gene = factor(
      Gene,
      levels = c(
        "AOC3",
        "MS4A7",
        "IL7R"
      )
    )
  ) %>%
  arrange(
    Gene,
    Model
  )


# ============================================================
# 19. REPORTING QC
# ============================================================

cat("\n============================================================\n")
cat("FINAL REPORTING QC\n")
cat("============================================================\n")


# Primary univariable cohort
stopifnot(
  all(univariable_table$N == 485),
  all(univariable_table$Events == 178)
)

# Model A cohort
stopifnot(
  all(modelA_table$N == 485),
  all(modelA_table$Events == 178)
)

# Model B cohort
stopifnot(
  all(modelB_table$N == 336),
  all(modelB_table$Events == 137)
)

# Gene completeness
stopifnot(
  identical(
    sort(unique(as.character(
      continuous_cox_table$Gene
    ))),
    sort(genes)
  )
)

cat("PASS: Univariable N = 485 / Events = 178\n")
cat("PASS: Model A N = 485 / Events = 178\n")
cat("PASS: Model B N = 336 / Events = 137\n")
cat("PASS: Three genes present in all constant-coefficient models\n")
cat("PASS: Targeted BH-FDR calculated separately within each model\n")
cat("PASS: No Z-score recalculation\n")
cat("PASS: No locked checkpoint overwritten\n")


# ============================================================
# 20. PRINT FINAL CONSTANT-COEFFICIENT TABLE
# ============================================================

cat("\n============================================================\n")
cat("FINAL THREE-GENE CONTINUOUS COX TABLE\n")
cat("============================================================\n")

print(
  continuous_cox_table_display %>%
    select(
      Gene,
      Model,
      N,
      Events,
      HR_CI,
      P_display,
      FDR_display,
      PH_P_display
    )
)

### Step 6D-1 reporting checkpoint
step6d1_checkpoint_metadata <- list(
  status = "LOCKED_REPORTING_CHECKPOINT",
  cohort_primary_N = 485,
  cohort_primary_events = 178,
  cohort_modelB_N = 336,
  cohort_modelB_events = 137,
  genes = genes,
  primary_model = "MS4A7 + Age + Gender + Stage",
  sensitivity_model = "MS4A7 + Age + Gender + T/N/M",
  standardization = "within locked 485-sample cohort; mean=0, SD=1",
  fdr_method = "BH within each three-gene model",
  time_varying_gene = "MS4A7",
  time_functions = c("linear", "sqrt", "log1p"),
  selected_time_function_unadjusted = best_time_function,
  selected_time_function_adjusted = best_time_function_adjusted,
  locked_objects_overwritten = FALSE
)

step6d1_checkpoint_objects <- c(
  "continuous_cox_table",
  "continuous_cox_table_display",
  "univariable_table",
  "modelA_table",
  "modelB_table",
  "figure4D_data",
  "ms4a7_tv_summary",
  "ms4a7_tv_adjusted_summary",
  "ms4a7_time_specific_hr",
  "modelB_reporting_cc"
)

save(
  continuous_cox_table,
  continuous_cox_table_display,
  univariable_table,
  modelA_table,
  modelB_table,
  figure4D_data,
  ms4a7_tv_summary,
  ms4a7_tv_adjusted_summary,
  ms4a7_time_specific_hr,
  modelB_reporting_cc,
  step6d1_checkpoint_metadata,
  step6d1_checkpoint_objects,
  file = file.path(
    output_dir,
    "Paper1_Step6D1_continuous_Cox_reporting_checkpoint.RData"
  )
)

cat(
  "\nSTEP 6D-1 REPORTING CHECKPOINT SAVED:\n",
  file.path(
    output_dir,
    "Paper1_Step6D1_continuous_Cox_reporting_checkpoint.RData"
  ),
  "\n"
)


# ============================================================
# 21. WRITE CSV FILES
# ============================================================

write.csv(
  continuous_cox_table,
  file = file.path(
    output_dir,
    "Step6D1_three_gene_continuous_Cox_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  continuous_cox_table_display,
  file = file.path(
    output_dir,
    "Step6D1_three_gene_continuous_Cox_results_display.csv"
  ),
  row.names = FALSE
)

write.csv(
  ms4a7_tv_summary,
  file = file.path(
    output_dir,
    "Step6D1_MS4A7_time_varying_unadjusted.csv"
  ),
  row.names = FALSE
)

write.csv(
  ms4a7_tv_adjusted_summary,
  file = file.path(
    output_dir,
    "Step6D1_MS4A7_time_varying_adjusted.csv"
  ),
  row.names = FALSE
)

write.csv(
  ms4a7_time_specific_hr,
  file = file.path(
    output_dir,
    "Step6D1_MS4A7_time_specific_HR.csv"
  ),
  row.names = FALSE
)

write.csv(
  figure4D_data,
  file = file.path(
    output_dir,
    "Step6D1_Figure4D_forest_plot_data.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. WRITE EXCEL WORKBOOK
# ============================================================

xlsx_file <- file.path(
  output_dir,
  "Step6D1_Survival_Reporting_Tables.xlsx"
)

wb <- createWorkbook()

addWorksheet(
  wb,
  "Three_gene_Cox"
)

writeData(
  wb,
  "Three_gene_Cox",
  continuous_cox_table_display
)

addWorksheet(
  wb,
  "MS4A7_TV_summary"
)

writeData(
  wb,
  "MS4A7_TV_summary",
  bind_rows(
    ms4a7_tv_summary,
    ms4a7_tv_adjusted_summary
  )
)

addWorksheet(
  wb,
  "MS4A7_time_specific_HR"
)

writeData(
  wb,
  "MS4A7_time_specific_HR",
  ms4a7_time_specific_hr
)

addWorksheet(
  wb,
  "Figure4D_data"
)

writeData(
  wb,
  "Figure4D_data",
  figure4D_data
)

addWorksheet(
  wb,
  "QC"
)

qc_table <- tibble(
  Item = c(
    "Primary cohort N",
    "Primary cohort events",
    "Model A N",
    "Model A events",
    "Model B N",
    "Model B events",
    "Genes",
    "Z-score recalculation",
    "Locked objects modified"
  ),
  Result = c(
    485,
    178,
    485,
    178,
    336,
    137,
    paste(genes, collapse = ", "),
    "No",
    "No"
  )
)

writeData(
  wb,
  "QC",
  qc_table
)

saveWorkbook(
  wb,
  xlsx_file,
  overwrite = TRUE
)


# ============================================================
# 23. SAVE NEW STEP 6D-1 CHECKPOINT
#
# This checkpoint contains reporting objects only.
# Existing locked checkpoints are NOT overwritten.
# ============================================================

step6d1_checkpoint <- list(
  
  checkpoint_name =
    "Paper1_Step6D1_Survival_Reporting",
  
  checkpoint_status =
    "LOCKED_REPORTING",
  
  analysis_scope =
    "Three-gene continuous-expression survival reporting",
  
  primary_cohort_N =
    485,
  
  primary_cohort_events =
    178,
  
  TNM_sensitivity_N =
    336,
  
  TNM_sensitivity_events =
    137,
  
  genes =
    genes,
  
  z_score_source =
    "survival_gene_primary_z",
  
  z_score_alias_object =
    "survival_gene_primary_z_downstream",
  
  no_z_score_recalculation =
    TRUE,
  
  univariable_results =
    univariable_table,
  
  modelA_results =
    modelA_table,
  
  modelB_results =
    modelB_table,
  
  continuous_cox_table =
    continuous_cox_table,
  
  continuous_cox_table_display =
    continuous_cox_table_display,
  
  ms4a7_tv_unadjusted =
    ms4a7_tv_summary,
  
  ms4a7_tv_adjusted =
    ms4a7_tv_adjusted_summary,
  
  ms4a7_time_specific_hr =
    ms4a7_time_specific_hr,
  
  figure4D_data =
    figure4D_data,
  
  output_directory =
    output_dir,
  
  excel_file =
    xlsx_file,
  
  generated_at =
    Sys.time()
)


save(
  step6d1_checkpoint,
  file = file.path(
    output_dir,
    "Paper1_Step6D1_Survival_Reporting_checkpoint.RData"
  )
)


# ============================================================
# 24. FINAL MESSAGE
# ============================================================

cat("\n============================================================\n")
cat("STEP 6D-1 COMPLETED\n")
cat("============================================================\n")

cat(
  "Output directory:\n",
  normalizePath(output_dir, winslash = "/"),
  "\n\n"
)

cat(
  "Excel workbook:\n",
  normalizePath(xlsx_file, winslash = "/"),
  "\n\n"
)

cat(
  "Checkpoint:\n",
  normalizePath(
    file.path(
      output_dir,
      "Paper1_Step6D1_Survival_Reporting_checkpoint.RData"
    ),
    winslash = "/"
  ),
  "\n"
)



# ============================================================
# 17. Prepare Figure 4
# ============================================================

cat("\n============================================================\n")
cat("PREPARE FIGURE 4 DATA\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 17.1 Use the validated continuous Cox reporting table
# ------------------------------------------------------------

figure4_data <- continuous_cox_table %>%
  filter(
    Model %in% c(
      "Univariable",
      "Adjusted: age + sex + stage",
      "Sensitivity: age + sex + T/N/M"
    )
  ) %>%
  mutate(
    
    # Reverse gene order for forest-plot display
    Gene = factor(
      Gene,
      levels = rev(c(
        "AOC3",
        "MS4A7",
        "IL7R"
      ))
    ),
    
    # Keep the three model types in the intended order
    Model = factor(
      Model,
      levels = c(
        "Univariable",
        "Adjusted: age + sex + stage",
        "Sensitivity: age + sex + T/N/M"
      )
    )
  ) %>%
  arrange(
    Gene,
    Model
  )


# ------------------------------------------------------------
# 17.2 QC Figure 4D dataset
# ------------------------------------------------------------

cat("\nFigure 4 dataset:\n")
print(figure4_data)


# Expected:
# 3 genes × 3 models = 9 rows

if (nrow(figure4_data) != 9) {
  stop(
    paste0(
      "Figure 4 dataset should contain 9 rows, got ",
      nrow(figure4_data)
    )
  )
}

if (anyNA(figure4_data$Gene)) {
  stop("Unexpected NA values in Figure 4 Gene factor.")
}

if (anyNA(figure4_data$Model)) {
  stop("Unexpected NA values in Figure 4 Model factor.")
}

if (any(!is.finite(figure4_data$HR))) {
  stop("Non-finite HR detected in Figure 4 dataset.")
}

if (any(!is.finite(figure4_data$CI_low)) ||
    any(!is.finite(figure4_data$CI_high))) {
  stop("Non-finite confidence interval detected in Figure 4 dataset.")
}


# ------------------------------------------------------------
# 17.3 Verify cohort sizes
# ------------------------------------------------------------

figure4_expected_N <- c(
  "Univariable" = 485,
  "Adjusted: age + sex + stage" = 485,
  "Sensitivity: age + sex + T/N/M" = 336
)

figure4_expected_events <- c(
  "Univariable" = 178,
  "Adjusted: age + sex + stage" = 178,
  "Sensitivity: age + sex + T/N/M" = 137
)

for (m in names(figure4_expected_N)) {
  
  tmp <- figure4_data[
    as.character(figure4_data$Model) == m,
    ,
    drop = FALSE
  ]
  
  if (any(tmp$N != figure4_expected_N[[m]])) {
    stop(
      paste0(
        "Unexpected N for Figure 4 model: ",
        m
      )
    )
  }
  
  if (any(tmp$Events != figure4_expected_events[[m]])) {
    stop(
      paste0(
        "Unexpected Events for Figure 4D model: ",
        m
      )
    )
  }
}


cat("\nFIGURE 4 DATASET QC: PASSED\n")
cat("Genes: AOC3, MS4A7, IL7R\n")
cat("Models: 3\n")
cat("Rows: 9\n")
cat("Univariable: N=485 / Events=178\n")
cat("Adjusted stage: N=485 / Events=178\n")
cat("TNM sensitivity: N=336 / Events=137\n")

# ============================================================
# 18. Figure 4 visual parameters
# ============================================================
# Paper 1
# Figure 4
# FINAL — Continuous Cox Forest Plot
#
# Style:
#   - No top legend
#   - Direct model labels
#   - Direct HR / 95% CI / P / FDR annotation
#   - Arial / Cairo PDF
# ============================================================

cat("\n============================================================\n")
cat("FIGURE 4 FINAL — CONTINUOUS COX FOREST PLOT\n")
cat("============================================================\n\n")

# ============================================================
# FONT DIAGNOSTIC
# ============================================================

cat("\n============================================================\n")
cat("FONT DIAGNOSTIC\n")
cat("============================================================\n")

cat("\n1. Windows Arial registration:\n")
print(windowsFonts())

cat("\n2. Arial mapping:\n")
print(windowsFonts("Arial"))

cat("\n3. Current graphics device:\n")
print(dev.cur())

cat("\n4. R version:\n")
print(R.version.string)

cat("\n5. ggplot2 version:\n")
print(as.character(packageVersion("ggplot2")))

cat("\n============================================================\n")

# ============================================================
# REGISTER ARIAL FOR R / WINDOWS
# ============================================================

windowsFonts(
  Arial = windowsFont("TT Arial")
)

cat("\nArial mapping after registration:\n")
print(windowsFonts("Arial"))

library(grid)

grid.newpage()

grid.text(
  "Arial font test: AOC3 MS4A7 IL7R",
  gp = gpar(
    family = "Arial",
    fontsize = 14
  )
)


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(grid)


# ------------------------------------------------------------
# 2. Validate input
# ------------------------------------------------------------

required_cols <- c(
  "Gene",
  "Model",
  "N",
  "Events",
  "HR",
  "CI_low",
  "CI_high",
  "P",
  "FDR"
)

missing_cols <- setdiff(
  required_cols,
  names(figure4_data)
)

if (length(missing_cols) > 0) {
  
  stop(
    paste0(
      "Missing required columns in figure4_data: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 3. Create plotting-only dataset
# ------------------------------------------------------------

fig4_plot <- figure4_data %>%
  
  mutate(
    
    Gene = as.character(Gene),
    Model = as.character(Model)
    
  ) %>%
  
  arrange(
    factor(
      Gene,
      levels = c(
        "AOC3",
        "MS4A7",
        "IL7R"
      )
    ),
    factor(
      Model,
      levels = c(
        "Univariable",
        "Adjusted: age + sex + stage",
        "Sensitivity: age + sex + T/N/M"
      )
    )
  )


# ------------------------------------------------------------
# 4. Locked-value QC
# ------------------------------------------------------------

if (nrow(fig4_plot) != 9) {
  stop("Figure 4 must contain exactly 9 rows.")
}

if (
  !setequal(
    unique(fig4_plot$Gene),
    c("AOC3", "MS4A7", "IL7R")
  )
) {
  stop("Unexpected genes in Figure 4.")
}

if (
  !setequal(
    unique(fig4_plot$Model),
    c(
      "Univariable",
      "Adjusted: age + sex + stage",
      "Sensitivity: age + sex + T/N/M"
    )
  )
) {
  stop("Unexpected models in Figure 4D.")
}

if (
  any(!is.finite(fig4_plot$HR)) ||
  any(!is.finite(fig4_plot$CI_low)) ||
  any(!is.finite(fig4_plot$CI_high))
) {
  stop("Non-finite HR or CI detected.")
}


# ------------------------------------------------------------
# 5. Create vertical offsets
#
# Three models are slightly separated within each gene.
# ------------------------------------------------------------

gene_order <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

model_order <- c(
  "Univariable",
  "Adjusted: age + sex + stage",
  "Sensitivity: age + sex + T/N/M"
)

fig4_plot <- fig4_plot %>%
  
  mutate(
    
    gene_index = match(
      Gene,
      gene_order
    ),
    
    model_index = match(
      Model,
      model_order
    ),
    
    y = gene_index +
      c(
        -0.18,
        0,
        0.18
      )[model_index]
    
  )


# ------------------------------------------------------------
# 6. Direct display labels
# ------------------------------------------------------------

fig4_plot <- fig4_plot %>%
  
  mutate(
    
    HR_CI_label = sprintf(
      "%.3f (%.3f–%.3f)",
      HR,
      CI_low,
      CI_high
    ),
    
    P_label = format.pval(
      P,
      digits = 3,
      eps = 0.001
    ),
    
    FDR_label = format.pval(
      FDR,
      digits = 3,
      eps = 0.001
    )
    
  )


# ------------------------------------------------------------
# 7. Model display labels
# ------------------------------------------------------------

fig4_plot <- fig4_plot %>%
  
  mutate(
    
    Model_short = case_when(
      
      Model == "Univariable" ~
        "Univariable",
      
      Model == "Adjusted: age + sex + stage" ~
        "Adjusted stage",
      
      Model == "Sensitivity: age + sex + T/N/M" ~
        "TNM sensitivity",
      
      TRUE ~ Model
      
    )
    
  )


# ------------------------------------------------------------
# 8. Colors
# ------------------------------------------------------------

fig4_model_colors <- c(
  
  "Univariable" = "#4DBBD5",
  
  "Adjusted: age + sex + stage" = "#E64B35",
  
  "Sensitivity: age + sex + T/N/M" = "grey45"
  
)

# ============================================================
# Figure 4— DATA STRUCTURE CHECK
# ============================================================

cat("\n============================================================\n")
cat("FIGURE 4 DATA STRUCTURE\n")
cat("============================================================\n")

print(names(figure4_data))

cat("\nModel values:\n")
print(unique(figure4_data$Model))

cat("\nGene values:\n")
print(unique(figure4_data$Gene))

cat("\nDimensions:\n")
print(dim(figure4_data))


# ============================================================
# Figure 4 — PREPARE FINAL PLOT DATA
# IMPORTANT:
# Do NOT modify the locked figure4_data object
# ============================================================

library(dplyr)

fig4_plot <- figure4_data

# ------------------------------------------------------------
# 1. Required-column QC
# ------------------------------------------------------------

required_cols <- c(
  "Gene",
  "Model",
  "N",
  "Events",
  "Beta",
  "HR",
  "CI_low",
  "CI_high",
  "P",
  "FDR",
  "PH_P"
)

missing_cols <- setdiff(
  required_cols,
  names(fig4_plot)
)

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# 2. Fix gene order
# ------------------------------------------------------------

fig4_plot$Gene <- factor(
  fig4_plot$Gene,
  levels = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )
)

# ------------------------------------------------------------
# 3. Compact model labels
# ------------------------------------------------------------

fig4_plot$Model_short <- case_when(
  
  as.character(fig4_plot$Model) ==
    "Univariable" ~
    "Univariable",
  
  as.character(fig4_plot$Model) ==
    "Adjusted: age + sex + stage" ~
    "Adjusted stage",
  
  as.character(fig4_plot$Model) ==
    "Sensitivity: age + sex + T/N/M" ~
    "TNM sensitivity",
  
  TRUE ~ as.character(fig4_plot$Model)
)

# ------------------------------------------------------------
# 4. Model order
# ------------------------------------------------------------

fig4_plot$Model_short <- factor(
  fig4_plot$Model_short,
  levels = c(
    "Univariable",
    "Adjusted stage",
    "TNM sensitivity"
  )
)

# ------------------------------------------------------------
# 5. Create y positions
#
# Each gene occupies one block:
#
# AOC3
#   Univariable
#   Adjusted stage
#   TNM sensitivity
#
# MS4A7
#   Univariable
#   Adjusted stage
#   TNM sensitivity
#
# IL7R
#   Univariable
#   Adjusted stage
#   TNM sensitivity
# ------------------------------------------------------------

gene_base <- c(
  "AOC3"  = 3.0,
  "MS4A7" = 2.0,
  "IL7R"  = 1.0
)

model_offset <- c(
  "Univariable"    =  0.18,
  "Adjusted stage" =  0.00,
  "TNM sensitivity" = -0.18
)

fig4_plot$y_position <-
  gene_base[
    as.character(fig4_plot$Gene)
  ] +
  model_offset[
    as.character(fig4_plot$Model_short)
  ]

# ------------------------------------------------------------
# 6. Create display labels
# ------------------------------------------------------------

fig4_plot$HR_CI_label <- sprintf(
  "%.3f (%.3f–%.3f)",
  fig4_plot$HR,
  fig4_plot$CI_low,
  fig4_plot$CI_high
)

fig4_plot$P_label <- sprintf(
  "%.4f",
  fig4_plot$P
)

fig4_plot$FDR_label <- sprintf(
  "%.4f",
  fig4_plot$FDR
)

# ------------------------------------------------------------
# 7. QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FIGURE 4 PLOT DATA QC\n")
cat("============================================================\n")

cat("\nDimensions:\n")
print(dim(fig4_plot))

cat("\nGenes:\n")
print(levels(fig4_plot$Gene))

cat("\nModels:\n")
print(levels(fig4_plot$Model_short))

cat("\nY positions:\n")
print(
  fig4_plot[
    c(
      "Gene",
      "Model_short",
      "y_position"
    )
  ]
)

cat("\nDisplay labels:\n")
print(
  fig4_plot[
    c(
      "Gene",
      "Model_short",
      "HR_CI_label",
      "P_label",
      "FDR_label"
    )
  ]
)

cat("\nOriginal figure4D_data unchanged:\n")
print(
  identical(
    names(figure4_data),
    c(
      "Gene",
      "Model",
      "N",
      "Events",
      "Beta",
      "HR",
      "CI_low",
      "CI_high",
      "P",
      "FDR",
      "PH_P"
    )
  )
)

# ============================================================
# Figure 4 — COMPACT FOREST PLOT VERSION
# ============================================================

library(ggplot2)

p_Figure4_OPT <- ggplot(
  fig4_plot,
  aes(
    x = HR,
    y = y_position,
    xmin = CI_low,
    xmax = CI_high
  )
) +
  
  # ----------------------------------------------------------
# HR = 1 reference line
# ----------------------------------------------------------
geom_vline(
  xintercept = 1,
  linetype = "dashed",
  linewidth = 0.50,
  colour = "grey45"
) +
  
  # ----------------------------------------------------------
# Light vertical reference guides
# ----------------------------------------------------------
geom_vline(
  xintercept = c(0.5, 0.75, 1.5, 2),
  linetype = "dotted",
  linewidth = 0.25,
  colour = "grey85"
) +
  
  # ----------------------------------------------------------
# Confidence intervals
# ----------------------------------------------------------
geom_errorbar(
  orientation = "y",
  linewidth = 0.65,
  width = 0.045
) +
  
  # ----------------------------------------------------------
# Point estimates
# Solid color circles
# ----------------------------------------------------------
geom_point(
  aes(colour = Model_short),
  shape = 16,
  size = 3.8
) +
  
  # ----------------------------------------------------------
# Gene labels
# ----------------------------------------------------------
geom_text(
  data = subset(
    fig4_plot,
    Model_short == "Univariable"
  ),
  aes(
    x = 0.435,
    y = y_position,
    label = Gene
  ),
  inherit.aes = FALSE,
  hjust = 0,
  size = 4.0,
  fontface = "bold",
  family = "Arial"
) +
  
  # ----------------------------------------------------------
# Model labels
# ----------------------------------------------------------
geom_text(
  aes(
    x = 0.435,
    y = y_position,
    label = Model_short
  ),
  inherit.aes = FALSE,
  hjust = 0,
  nudge_y = -0.105,
  size = 2.75,
  family = "Arial"
) +
  
  # ----------------------------------------------------------
# HR and 95% CI
# ----------------------------------------------------------
geom_text(
  aes(
    x = 1.95,
    y = y_position,
    label = HR_CI_label
  ),
  inherit.aes = FALSE,
  hjust = 0,
  size = 3.05,
  family = "Arial"
) +
  
  # ----------------------------------------------------------
# P
# ----------------------------------------------------------
geom_text(
  aes(
    x = 2.95,
    y = y_position,
    label = P_label
  ),
  inherit.aes = FALSE,
  hjust = 0,
  size = 3.05,
  family = "Arial"
) +
  
  # ----------------------------------------------------------
# FDR
# ----------------------------------------------------------
geom_text(
  aes(
    x = 3.45,
    y = y_position,
    label = FDR_label
  ),
  inherit.aes = FALSE,
  hjust = 0,
  size = 3.05,
  family = "Arial"
) +
  
  # ----------------------------------------------------------
# Headers
# ----------------------------------------------------------
annotate(
  "text",
  x = 1.95,
  y = 3.55,
  label = "HR (95% CI)",
  hjust = 0,
  size = 3.15,
  fontface = "bold",
  family = "Arial"
) +
  
  annotate(
    "text",
    x = 2.95,
    y = 3.55,
    label = "P",
    hjust = 0,
    size = 3.15,
    fontface = "bold",
    family = "Arial"
  ) +
  
  annotate(
    "text",
    x = 3.45,
    y = 3.55,
    label = "FDR",
    hjust = 0,
    size = 3.15,
    fontface = "bold",
    family = "Arial"
  ) +
  
  # ----------------------------------------------------------
# Log HR axis
# ----------------------------------------------------------
scale_x_log10(
  limits = c(0.42, 3.85),
  breaks = c(
    0.5,
    0.75,
    1,
    1.5,
    2,
    3
  ),
  labels = c(
    "0.5",
    "0.75",
    "1",
    "1.5",
    "2",
    "3"
  ),
  expand = c(0, 0)
) +
  
  # ----------------------------------------------------------
# Solid point colors
# ----------------------------------------------------------
scale_colour_manual(
  values = c(
    "Univariable" = "#3DB7D3",
    "Adjusted stage" = "#E64B35",
    "TNM sensitivity" = "#7F7F7F"
  ),
  guide = "none"
) +
  
  # ----------------------------------------------------------
# Y axis
# ----------------------------------------------------------
scale_y_continuous(
  breaks = NULL,
  limits = c(0.60, 3.70),
  expand = c(0, 0)
) +
  
  # ----------------------------------------------------------
# Theme
# ----------------------------------------------------------
theme_classic(
  base_family = "Arial",
  base_size = 10.5
) +
  
  theme(
    axis.title = element_blank(),
    
    axis.text.x = element_text(
      size = 9.5,
      colour = "black"
    ),
    
    axis.line.y = element_blank(),
    
    axis.ticks.y = element_blank(),
    
    axis.ticks.x = element_line(
      linewidth = 0.40
    ),
    
    panel.grid = element_blank(),
    
    plot.margin = margin(
      t = 10,
      r = 6,
      b = 8,
      l = 5
    )
  )


# ------------------------------------------------------------
# 10. Display
# ------------------------------------------------------------
print(p_Figure4_OPT)



# ============================================================
# Figure 4 — FINAL HIGH-RESOLUTION EXPORT
# ============================================================

ggsave(
  filename = figure4_opt_pdf,
  plot = p_Figure4_OPT,
  width = 8.8,
  height = 4.6,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

ggsave(
  filename = figure4_opt_tiff,
  plot = p_Figure4_OPT,
  width = 8.8,
  height = 4.6,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = figure4_opt_png,
  plot = p_Figure4_OPT,
  width = 8.8,
  height = 4.6,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = figure4_opt_jpeg,
  plot = p_Figure4_OPT,
  width = 8.8,
  height = 4.6,
  units = "in",
  dpi = 600,
  quality = 100,
  bg = "white"
)

# ============================================================
# Figure 4 — EXPORT QC
# ============================================================

cat("\n============================================================\n")
cat("FIGURE 4 FINAL PLOT QC\n")
cat("============================================================\n")

cat("\nPlot object exists:\n")
print(exists("p_Figure4_OPT"))

cat("\nPlot data dimensions:\n")
print(dim(fig4_plot))

cat("\nNumber of rows:\n")
print(nrow(fig4_plot))

cat("\nGene labels:\n")
print(unique(as.character(fig4_plot$Gene)))

cat("\nModel labels:\n")
print(levels(fig4_plot$Model_short))

cat("\nHR range:\n")
print(range(fig4_plot$HR))

cat("\nCI range:\n")
print(range(fig4_plot$CI_low))
print(range(fig4_plot$CI_high))

cat("\nAll display values complete:\n")
print(
  all(
    !is.na(fig4_plot$HR_CI_label),
    !is.na(fig4_plot$P_label),
    !is.na(fig4_plot$FDR_label),
    !is.na(fig4_plot$y_position)
  )
)


# ============================================================
# Paper 1
# STEP 6X
# Genome-wide survival association analysis
#
# Independent analysis
# Does NOT modify Step 6B / 6C / 6D1 locked objects
# ============================================================

library(survival)
library(dplyr)
library(tibble)
library(ggplot2)
library(openxlsx)

cat("\n============================================================\n")
cat("STEP 6X: GENOME-WIDE SURVIVAL ASSOCIATION ANALYSIS\n")
cat("============================================================\n")

step6X_root <- file.path(
  "input_dir/Paper 1_Step 6-survival analysis",
  "STEP_6X_genome_wide_survival"
)

step6X_checkpoint <- file.path(
  step6X_root,
  "checkpoint"
)

step6X_results <- file.path(
  step6X_root,
  "results"
)

step6X_figure <- file.path(
  step6X_root,
  "figure"
)

step6X_audit <- file.path(
  step6X_root,
  "audit"
)

for (x in c(
  step6X_checkpoint,
  step6X_results,
  step6X_figure,
  step6X_audit
)) {
  if (!dir.exists(x)) {
    dir.create(x, recursive = TRUE)
  }
}

# ============================================================
# 6X-1
# Recover locked primary survival cohort
# ============================================================

if (!exists("survival_gene_primary_z")) {
  stop(
    "Locked object 'survival_gene_primary_z' is not available."
  )
}

if (!exists("locked_primary_samples")) {
  stop(
    "Locked object 'locked_primary_samples' is not available."
  )
}

locked_samples_6X <- as.character(
  locked_primary_samples
)

cat("\nLocked survival cohort:\n")
cat("N =", length(locked_samples_6X), "\n")

if (length(locked_samples_6X) != 485) {
  stop(
    "Locked survival cohort is not N=485. ",
    "Do not proceed without resolving the cohort."
  )
}

if (anyDuplicated(locked_samples_6X)) {
  stop("Duplicate samples detected in locked survival cohort.")
}

cat("Duplicate samples: NONE\n")


# ============================================================
# 6X-2
# Load genome-wide expression matrix
# ============================================================

wgcna_expr_file <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt"
)

if (!file.exists(wgcna_expr_file)) {
  stop(
    "Genome-wide expression file not found:\n",
    wgcna_expr_file
  )
}

expr_6X <- read.table(
  wgcna_expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\nRaw expression dimensions:\n")
print(dim(expr_6X))


# ============================================================
# 6X-3
# Identify gene and sample columns
# ============================================================

cat("\nFirst 10 columns:\n")
print(names(expr_6X)[1:min(10, ncol(expr_6X))])

cat("\nNumber of columns:\n")
print(ncol(expr_6X))

gene_col_6X <- names(expr_6X)[1]

sample_cols_6X <- setdiff(
  names(expr_6X),
  gene_col_6X
)

cat("\nGene column:\n")
print(gene_col_6X)

cat("\nNumber of sample columns:\n")
print(length(sample_cols_6X))


# ============================================================
# 6X-4
# Construct gene × sample expression matrix
# ============================================================

gene_ids_6X <- as.character(
  expr_6X[[gene_col_6X]]
)

if (anyDuplicated(gene_ids_6X)) {
  stop(
    "Duplicated gene identifiers detected in genome-wide expression matrix."
  )
}

expr_mat_6X <- as.matrix(
  expr_6X[, sample_cols_6X, drop = FALSE]
)

rownames(expr_mat_6X) <- gene_ids_6X

storage.mode(expr_mat_6X) <- "numeric"

cat("\nExpression matrix:\n")
print(dim(expr_mat_6X))


# ============================================================
# 6X-5
# Align expression matrix to locked survival cohort
# ============================================================

common_samples_6X <- intersect(
  locked_samples_6X,
  colnames(expr_mat_6X)
)

cat("\nCommon samples:\n")
print(length(common_samples_6X))

if (length(common_samples_6X) != 485) {
  stop(
    "Expression/survival overlap is not exactly N=485."
  )
}

expr_mat_6X <- expr_mat_6X[
  ,
  locked_samples_6X,
  drop = FALSE
]

if (!identical(
  colnames(expr_mat_6X),
  locked_samples_6X
)) {
  stop(
    "Sample order mismatch after alignment."
  )
}

cat("Sample alignment: PASSED\n")


# ============================================================
# 6X-6
# Construct locked survival phenotype
# ============================================================

surv_6X <- survival_gene_primary_z[
  ,
  c(
    "Sample",
    "OS.time",
    "OS.status",
    "OS_event"
  )
]

surv_6X$Sample <- as.character(
  surv_6X$Sample
)

surv_6X <- surv_6X[
  match(
    locked_samples_6X,
    surv_6X$Sample
  ),
  ,
  drop = FALSE
]

if (!identical(
  surv_6X$Sample,
  locked_samples_6X
)) {
  stop(
    "Survival data order does not match locked cohort."
  )
}

cat("\nSurvival cohort QC:\n")
cat("N =", nrow(surv_6X), "\n")
cat("Events =", sum(surv_6X$OS_event), "\n")
cat("Censored =", sum(surv_6X$OS_event == 0), "\n")


# ============================================================
# STEP 6X-8A
# Genome-wide expression matrix preparation
# Do NOT modify locked survival objects
# ============================================================

cat("\n============================================================\n")
cat("STEP 6X-8A: Genome-wide expression preparation\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. Check current objects
# ------------------------------------------------------------

cat("\nObjects related to expression:\n")

print(
  grep(
    "expr|expression|wgcna|candidate",
    ls(),
    ignore.case = TRUE,
    value = TRUE
  )
)

# ------------------------------------------------------------
# 2. Check locked survival cohort
# ------------------------------------------------------------

stopifnot(exists("surv_6X"))

cat("\nSurvival cohort:\n")
cat("N =", nrow(surv_6X), "\n")
cat("Events =", sum(surv_6X$OS_event), "\n")
cat("Censored =", sum(surv_6X$OS_event == 0), "\n")

stopifnot(nrow(surv_6X) == 485)
stopifnot(sum(surv_6X$OS_event) == 178)
stopifnot(!anyDuplicated(surv_6X$Sample))

# ------------------------------------------------------------
# 3. Locate the genome-wide expression file
# ------------------------------------------------------------

expr_file <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt"
)

stopifnot(file.exists(expr_file))

cat("\nExpression file:\n")
cat(expr_file, "\n")

# ------------------------------------------------------------
# 4. Read expression matrix
# ------------------------------------------------------------

expr_raw_6X <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

cat("\nRaw expression dimensions:\n")
print(dim(expr_raw_6X))

cat("\nFirst columns:\n")
print(head(colnames(expr_raw_6X), 10))

cat("\nFirst rows:\n")
print(expr_raw_6X[1:min(3, nrow(expr_raw_6X)), 1:min(5, ncol(expr_raw_6X))])

# ------------------------------------------------------------
# 5. Identify gene ID column
# ------------------------------------------------------------

gene_id_col <- colnames(expr_raw_6X)[1]

cat("\nGene ID column =", gene_id_col, "\n")

gene_ids_6X <- as.character(expr_raw_6X[[gene_id_col]])

cat(
  "Number of genes =",
  length(gene_ids_6X),
  "\n"
)

cat(
  "Duplicated gene IDs =",
  sum(duplicated(gene_ids_6X)),
  "\n"
)

stopifnot(!anyDuplicated(gene_ids_6X))

# ------------------------------------------------------------
# 6. Build gene x sample matrix
# ------------------------------------------------------------

expr_mat_6X <- as.matrix(
  expr_raw_6X[, -1, drop = FALSE]
)

rownames(expr_mat_6X) <- gene_ids_6X

storage.mode(expr_mat_6X) <- "numeric"

cat("\nExpression matrix dimensions:\n")
print(dim(expr_mat_6X))

# ------------------------------------------------------------
# 7. Check whether samples are columns
# ------------------------------------------------------------

cat("\nNumber of expression samples =", ncol(expr_mat_6X), "\n")

# ------------------------------------------------------------
# 8. Align expression matrix to locked 485 samples
# ------------------------------------------------------------

locked_samples_6X <- as.character(surv_6X$Sample)

common_samples_6X <- intersect(
  locked_samples_6X,
  colnames(expr_mat_6X)
)

cat("\nCommon samples =", length(common_samples_6X), "\n")

if (length(common_samples_6X) != 485) {
  cat("\nMissing from expression matrix:\n")
  print(setdiff(locked_samples_6X, colnames(expr_mat_6X)))
  
  cat("\nExtra expression samples not in survival cohort:\n")
  print(setdiff(colnames(expr_mat_6X), locked_samples_6X))
}

stopifnot(length(common_samples_6X) == 485)

# Exact ordering according to survival cohort
expr_mat_6X <- expr_mat_6X[
  ,
  locked_samples_6X,
  drop = FALSE
]

# ------------------------------------------------------------
# 9. Final alignment QC
# ------------------------------------------------------------

stopifnot(
  identical(
    colnames(expr_mat_6X),
    locked_samples_6X
  )
)

cat("\nFinal expression matrix:\n")
cat(
  "Genes =",
  nrow(expr_mat_6X),
  "\n"
)

cat(
  "Samples =",
  ncol(expr_mat_6X),
  "\n"
)

cat(
  "Sample order identical =",
  identical(
    colnames(expr_mat_6X),
    locked_samples_6X
  ),
  "\n"
)

cat(
  "Any NA =",
  anyNA(expr_mat_6X),
  "\n"
)

cat(
  "Any Inf =",
  any(!is.finite(expr_mat_6X)),
  "\n"
)

# ------------------------------------------------------------
# 10. Check expression scale
# ------------------------------------------------------------

cat("\nExpression summary:\n")

print(
  summary(
    as.numeric(expr_mat_6X[
      sample(
        seq_len(nrow(expr_mat_6X)),
        min(1000, nrow(expr_mat_6X))
      ),
      ,
      drop = FALSE
    ])
  )
)

cat("\n============================================================\n")
cat("STEP 6X-8A QC completed\n")
cat("============================================================\n")


# ============================================================
# STEP 6X-8B
# Define genome-wide survival gene universe
#
# Purpose:
#   Restrict the survival analysis to the same
#   MAD-positive gene universe used for the revised
#   tumor-only WGCNA analysis.
#
# IMPORTANT:
#   This step does NOT modify:
#   - Step 6B objects
#   - Step 6C objects
#   - locked survival cohort
#   - WGCNA objects
# ============================================================

cat("\n============================================================\n")
cat("STEP 6X-8B: Define genome-wide survival gene universe\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Check candidate WGCNA objects
# ------------------------------------------------------------

cat("\nChecking WGCNA objects...\n")

cat(
  "datExpr_MAD dimensions: ",
  paste(dim(datExpr_MAD), collapse = " x "),
  "\n"
)

cat(
  "datExpr_MAD_recovery dimensions: ",
  paste(dim(datExpr_MAD_recovery), collapse = " x "),
  "\n"
)


# ------------------------------------------------------------
# 2. Identify WGCNA gene universe
# ------------------------------------------------------------

# datExpr_MAD is expected to be:
# samples × genes

wgcna_gene_universe_6X <- colnames(datExpr_MAD)

if (is.null(wgcna_gene_universe_6X)) {
  stop("datExpr_MAD has no column names; cannot define gene universe.")
}

cat(
  "\nWGCNA gene universe size =",
  length(wgcna_gene_universe_6X),
  "\n"
)

cat(
  "Duplicated WGCNA gene IDs =",
  sum(duplicated(wgcna_gene_universe_6X)),
  "\n"
)

stopifnot(!anyDuplicated(wgcna_gene_universe_6X))


# ------------------------------------------------------------
# 3. Check expected size
# ------------------------------------------------------------

stopifnot(
  length(wgcna_gene_universe_6X) == 31609
)


# ------------------------------------------------------------
# 4. Compare with genome-wide expression matrix
# ------------------------------------------------------------

expr_gene_ids_6X <- rownames(expr_mat_6X)

cat(
  "\nExpression matrix gene count =",
  length(expr_gene_ids_6X),
  "\n"
)

cat(
  "WGCNA gene count =",
  length(wgcna_gene_universe_6X),
  "\n"
)


genes_missing_from_expression_6X <- setdiff(
  wgcna_gene_universe_6X,
  expr_gene_ids_6X
)

genes_extra_in_expression_6X <- setdiff(
  expr_gene_ids_6X,
  wgcna_gene_universe_6X
)

cat(
  "\nWGCNA genes missing from expression matrix =",
  length(genes_missing_from_expression_6X),
  "\n"
)

cat(
  "Expression genes outside WGCNA universe =",
  length(genes_extra_in_expression_6X),
  "\n"
)


# ------------------------------------------------------------
# 5. Require complete overlap
# ------------------------------------------------------------

stopifnot(
  length(genes_missing_from_expression_6X) == 0
)


# ------------------------------------------------------------
# 6. Restrict expression matrix to WGCNA universe
# ------------------------------------------------------------

expr_mat_genome_6X <- expr_mat_6X[
  wgcna_gene_universe_6X,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 7. Final dimensions
# ------------------------------------------------------------

cat("\nFinal genome-wide survival matrix:\n")

cat(
  "Genes =",
  nrow(expr_mat_genome_6X),
  "\n"
)

cat(
  "Samples =",
  ncol(expr_mat_genome_6X),
  "\n"
)


# ------------------------------------------------------------
# 8. Sample-order QC
# ------------------------------------------------------------

stopifnot(
  identical(
    colnames(expr_mat_genome_6X),
    locked_samples_6X
  )
)

cat(
  "Sample order identical =",
  identical(
    colnames(expr_mat_genome_6X),
    locked_samples_6X
  ),
  "\n"
)


# ------------------------------------------------------------
# 9. Gene-order QC
# ------------------------------------------------------------

stopifnot(
  identical(
    rownames(expr_mat_genome_6X),
    wgcna_gene_universe_6X
  )
)

cat(
  "Gene order identical to WGCNA universe =",
  identical(
    rownames(expr_mat_genome_6X),
    wgcna_gene_universe_6X
  ),
  "\n"
)


# ------------------------------------------------------------
# 10. Numeric / missing-value QC
# ------------------------------------------------------------

cat(
  "Any NA =",
  anyNA(expr_mat_genome_6X),
  "\n"
)

cat(
  "Any Inf =",
  any(!is.finite(expr_mat_genome_6X)),
  "\n"
)

stopifnot(!anyNA(expr_mat_genome_6X))
stopifnot(all(is.finite(expr_mat_genome_6X)))


# ------------------------------------------------------------
# 11. Verify the three core genes
# ------------------------------------------------------------

core_genes_6X <- c(
  "ENSG00000131471",  # AOC3
  "ENSG00000166927",  # MS4A7
  "ENSG00000168685"   # IL7R
)

cat("\nCore gene presence:\n")

print(
  data.frame(
    Ensembl = core_genes_6X,
    Present = core_genes_6X %in%
      rownames(expr_mat_genome_6X),
    stringsAsFactors = FALSE
  )
)

stopifnot(
  all(
    core_genes_6X %in%
      rownames(expr_mat_genome_6X)
  )
)


# ------------------------------------------------------------
# 12. Save a lightweight audit object
# ------------------------------------------------------------

genome_survival_gene_universe_6X <- list(
  gene_universe = wgcna_gene_universe_6X,
  n_genes = length(wgcna_gene_universe_6X),
  n_samples = ncol(expr_mat_genome_6X),
  survival_N = nrow(surv_6X),
  survival_events = sum(surv_6X$OS_event),
  survival_censored = sum(surv_6X$OS_event == 0),
  core_genes = core_genes_6X,
  expression_source = expr_file,
  wgcna_gene_source = "datExpr_MAD",
  sample_order_locked = TRUE,
  gene_order_locked = TRUE
)


# ------------------------------------------------------------
# 13. Final checkpoint message
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 6X-8B QC SUMMARY\n")
cat("============================================================\n")

cat(
  "Genome-wide genes =",
  nrow(expr_mat_genome_6X),
  "\n"
)

cat(
  "Survival samples =",
  ncol(expr_mat_genome_6X),
  "\n"
)

cat(
  "Events =",
  sum(surv_6X$OS_event),
  "\n"
)

cat(
  "Censored =",
  sum(surv_6X$OS_event == 0),
  "\n"
)

cat(
  "Sample alignment =",
  identical(
    colnames(expr_mat_genome_6X),
    locked_samples_6X
  ),
  "\n"
)

cat(
  "Gene alignment =",
  identical(
    rownames(expr_mat_genome_6X),
    wgcna_gene_universe_6X
  ),
  "\n"
)

cat(
  "NA-free =",
  !anyNA(expr_mat_genome_6X),
  "\n"
)

cat(
  "Inf-free =",
  all(is.finite(expr_mat_genome_6X)),
  "\n"
)

cat("\nSTEP 6X-8B: PASSED\n")
cat("============================================================\n")


##STEP 6X-8C：Genome-wide Cox
# ============================================================
# PAPER 1
# STEP 6X-8C
# Genome-wide univariable Cox survival analysis
#
# Purpose:
#   Evaluate genome-wide survival associations across the
#   31,609-gene MAD-positive universe used for revised
#   tumor-only WGCNA.
#
# Locked survival cohort:
#   N = 485
#   Events = 178
#   Censored = 307
#
# IMPORTANT:
#   - Independent from Step 6B / 6C
#   - Does NOT modify any locked objects
#   - HR is per 1-SD increase in gene expression
#   - Expression is already log2-transformed
#   - No additional log transformation
#   - No genome-wide PH testing at this stage
#   - Warnings and failures are explicitly audited
# ============================================================

cat("\n============================================================\n")
cat("STEP 6X-8C: Genome-wide univariable Cox analysis\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages_6X <- c(
  "survival",
  "dplyr"
)

missing_packages_6X <- required_packages_6X[
  !vapply(
    required_packages_6X,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages_6X) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages_6X, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 2. Verify required objects
# ------------------------------------------------------------

required_objects_6X <- c(
  "surv_6X",
  "expr_mat_genome_6X",
  "wgcna_gene_universe_6X"
)

missing_objects_6X <- required_objects_6X[
  !vapply(
    required_objects_6X,
    exists,
    logical(1),
    inherits = FALSE
  )
]

if (length(missing_objects_6X) > 0) {
  stop(
    "Missing required objects: ",
    paste(missing_objects_6X, collapse = ", ")
  )
}
##First, check the current state of the R session.
cat("\nCurrent STEP 6X objects:\n")

print(
  intersect(
    c(
      "surv_6X",
      "expr_mat_6X",
      "expr_mat_genome_6X",
      "wgcna_gene_universe_6X",
      "locked_primary_samples_6X",
      "genome_survival_gene_universe_6X"
    ),
    ls()
  )
)

# ============================================================
# STEP 6X-8B-R
# Recover missing helper object and re-verify locked objects
# ============================================================

cat("\n============================================================\n")
cat("STEP 6X-8B-R: Recover STEP 6X helper object\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Recover locked sample vector
# ------------------------------------------------------------

locked_primary_samples_6X <- as.character(
  surv_6X$Sample
)


# ------------------------------------------------------------
# 2. Survival cohort QC
# ------------------------------------------------------------

stopifnot(
  nrow(surv_6X) == 485
)

stopifnot(
  sum(surv_6X$OS_event) == 178
)

stopifnot(
  sum(surv_6X$OS_event == 0) == 307
)

stopifnot(
  !anyDuplicated(surv_6X$Sample)
)


# ------------------------------------------------------------
# 3. Expression matrix QC
# ------------------------------------------------------------

stopifnot(
  nrow(expr_mat_genome_6X) == 31609
)

stopifnot(
  ncol(expr_mat_genome_6X) == 485
)

stopifnot(
  identical(
    colnames(expr_mat_genome_6X),
    locked_primary_samples_6X
  )
)

stopifnot(
  identical(
    rownames(expr_mat_genome_6X),
    wgcna_gene_universe_6X
  )
)

stopifnot(
  !anyNA(expr_mat_genome_6X)
)

stopifnot(
  all(is.finite(expr_mat_genome_6X))
)


# ------------------------------------------------------------
# 4. Verify core genes
# ------------------------------------------------------------

core_genes_6X <- c(
  "ENSG00000131471",  # AOC3
  "ENSG00000166927",  # MS4A7
  "ENSG00000168685"   # IL7R
)

stopifnot(
  all(
    core_genes_6X %in%
      rownames(expr_mat_genome_6X)
  )
)


# ------------------------------------------------------------
# 5. Final recovery report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 6X-8B-R QC SUMMARY\n")
cat("============================================================\n")

cat(
  "Survival N =",
  nrow(surv_6X),
  "\n"
)

cat(
  "Events =",
  sum(surv_6X$OS_event),
  "\n"
)

cat(
  "Censored =",
  sum(surv_6X$OS_event == 0),
  "\n"
)

cat(
  "Genome-wide genes =",
  nrow(expr_mat_genome_6X),
  "\n"
)

cat(
  "Genome-wide samples =",
  ncol(expr_mat_genome_6X),
  "\n"
)

cat(
  "Sample alignment =",
  identical(
    colnames(expr_mat_genome_6X),
    locked_primary_samples_6X
  ),
  "\n"
)

cat(
  "Gene alignment =",
  identical(
    rownames(expr_mat_genome_6X),
    wgcna_gene_universe_6X
  ),
  "\n"
)

cat(
  "NA-free =",
  !anyNA(expr_mat_genome_6X),
  "\n"
)

cat(
  "Inf-free =",
  all(is.finite(expr_mat_genome_6X)),
  "\n"
)

cat(
  "Core genes present =",
  all(
    core_genes_6X %in%
      rownames(expr_mat_genome_6X)
  ),
  "\n"
)

cat("\nSTEP 6X-8B-R: PASSED\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 3. Locked cohort QC
# ------------------------------------------------------------

stopifnot(nrow(surv_6X) == 485)
stopifnot(sum(surv_6X$OS_event) == 178)
stopifnot(sum(surv_6X$OS_event == 0) == 307)
stopifnot(!anyDuplicated(surv_6X$Sample))


# ------------------------------------------------------------
# 4. Expression matrix QC
# ------------------------------------------------------------

stopifnot(
  nrow(expr_mat_genome_6X) == 31609
)

stopifnot(
  ncol(expr_mat_genome_6X) == 485
)

stopifnot(
  identical(
    colnames(expr_mat_genome_6X),
    as.character(surv_6X$Sample)
  )
)

stopifnot(
  identical(
    rownames(expr_mat_genome_6X),
    wgcna_gene_universe_6X
  )
)

stopifnot(!anyNA(expr_mat_genome_6X))
stopifnot(all(is.finite(expr_mat_genome_6X)))


# ------------------------------------------------------------
# 5. Prepare survival data
# ------------------------------------------------------------

time_6X <- as.numeric(surv_6X$OS.time)
event_6X <- as.numeric(surv_6X$OS_event)

stopifnot(length(time_6X) == 485)
stopifnot(length(event_6X) == 485)
stopifnot(all(is.finite(time_6X)))
stopifnot(all(event_6X %in% c(0, 1)))


# ------------------------------------------------------------
# 6. Gene-wise Cox function
# ------------------------------------------------------------

run_one_gene_cox_6X <- function(
    expr_gene,
    gene_id,
    time,
    event) {
  
  result <- data.frame(
    Gene = gene_id,
    N = NA_integer_,
    Events = NA_integer_,
    Mean = NA_real_,
    SD = NA_real_,
    Beta = NA_real_,
    HR = NA_real_,
    CI_low = NA_real_,
    CI_high = NA_real_,
    P = NA_real_,
    Status = "FAILED",
    Warning = NA_character_,
    Error = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # ----------------------------------------------------------
  # Expression QC
  # ----------------------------------------------------------
  
  if (length(expr_gene) != length(time)) {
    result$Error <- "Expression length does not match survival cohort"
    return(result)
  }
  
  if (any(!is.finite(expr_gene))) {
    result$Error <- "Non-finite expression values"
    return(result)
  }
  
  gene_sd <- stats::sd(expr_gene)
  
  result$Mean <- mean(expr_gene)
  result$SD <- gene_sd
  
  if (!is.finite(gene_sd) || gene_sd <= 0) {
    result$Error <- "Zero or non-finite expression SD"
    return(result)
  }
  
  # ----------------------------------------------------------
  # Standardize within the locked 485-observation cohort
  # ----------------------------------------------------------
  
  gene_z <- as.numeric(
    (expr_gene - mean(expr_gene)) / gene_sd
  )
  
  # ----------------------------------------------------------
  # Cox model
  # ----------------------------------------------------------
  
  dat_gene <- data.frame(
    OS_time = time,
    OS_event = event,
    gene_z = gene_z
  )
  
  warning_messages <- character(0)
  
  fit <- tryCatch(
    
    withCallingHandlers(
      
      survival::coxph(
        survival::Surv(OS_time, OS_event) ~ gene_z,
        data = dat_gene,
        ties = "efron",
        x = FALSE,
        y = FALSE,
        model = FALSE
      ),
      
      warning = function(w) {
        
        warning_messages <<- c(
          warning_messages,
          conditionMessage(w)
        )
        
        invokeRestart("muffleWarning")
      }
    ),
    
    error = function(e) {
      
      result$Error <<- conditionMessage(e)
      
      return(NULL)
    }
  )
  
  # ----------------------------------------------------------
  # Model failure
  # ----------------------------------------------------------
  
  if (is.null(fit)) {
    
    if (length(warning_messages) > 0) {
      result$Warning <- paste(
        unique(warning_messages),
        collapse = " | "
      )
    }
    
    return(result)
  }
  
  # ----------------------------------------------------------
  # Extract model statistics
  # ----------------------------------------------------------
  
  coef_table <- summary(fit)$coefficients
  conf_table <- summary(fit)$conf.int
  
  beta <- unname(coef(fit)[1])
  
  HR <- unname(
    conf_table[1, "exp(coef)"]
  )
  
  CI_low <- unname(
    conf_table[1, "lower .95"]
  )
  
  CI_high <- unname(
    conf_table[1, "upper .95"]
  )
  
  P <- unname(
    coef_table[1, "Pr(>|z|)"]
  )
  
  # ----------------------------------------------------------
  # Final numerical QC
  # ----------------------------------------------------------
  
  result$N <- nrow(dat_gene)
  result$Events <- sum(dat_gene$OS_event)
  
  result$Beta <- beta
  result$HR <- HR
  result$CI_low <- CI_low
  result$CI_high <- CI_high
  result$P <- P
  
  if (
    !all(
      is.finite(
        c(beta, HR, CI_low, CI_high, P)
      )
    )
  ) {
    
    result$Status <- "FAILED"
    
    result$Error <- paste(
      c(
        result$Error,
        "Non-finite model statistic"
      )[!is.na(
        c(result$Error, "Non-finite model statistic")
      )],
      collapse = " | "
    )
    
  } else {
    
    result$Status <- "OK"
  }
  
  if (length(warning_messages) > 0) {
    
    result$Warning <- paste(
      unique(warning_messages),
      collapse = " | "
    )
  }
  
  return(result)
}


# ------------------------------------------------------------
# 7. Run genome-wide Cox models
# ------------------------------------------------------------

gene_ids_run_6X <- rownames(expr_mat_genome_6X)

n_genes_6X <- length(gene_ids_run_6X)

cat(
  "\nNumber of genes to analyze =",
  n_genes_6X,
  "\n"
)

cat(
  "Expected = 31609\n"
)


genome_cox_results_list_6X <- vector(
  "list",
  n_genes_6X
)


start_time_6X <- Sys.time()


for (i in seq_len(n_genes_6X)) {
  
  gene_id <- gene_ids_run_6X[i]
  
  expr_gene <- as.numeric(
    expr_mat_genome_6X[i, ]
  )
  
  genome_cox_results_list_6X[[i]] <-
    run_one_gene_cox_6X(
      expr_gene = expr_gene,
      gene_id = gene_id,
      time = time_6X,
      event = event_6X
    )
  
  if (
    i %% 1000 == 0 ||
    i == 1 ||
    i == n_genes_6X
  ) {
    
    elapsed <- difftime(
      Sys.time(),
      start_time_6X,
      units = "mins"
    )
    
    cat(
      "Progress:",
      i,
      "/",
      n_genes_6X,
      "| elapsed:",
      round(
        as.numeric(elapsed),
        2
      ),
      "min\n"
    )
  }
}


# ------------------------------------------------------------
# 8. Combine results
# ------------------------------------------------------------

genome_cox_results_6X <- dplyr::bind_rows(
  genome_cox_results_list_6X
)


# ------------------------------------------------------------
# 9. Genome-wide BH-FDR
# ------------------------------------------------------------

genome_cox_results_6X$FDR <- NA_real_

ok_idx_6X <- which(
  genome_cox_results_6X$Status == "OK" &
    is.finite(genome_cox_results_6X$P)
)

genome_cox_results_6X$FDR[
  ok_idx_6X
] <- p.adjust(
  genome_cox_results_6X$P[
    ok_idx_6X
  ],
  method = "BH"
)

# ------------------------------------------------------------
# 10. Derived quantities for Supplementary Figure X
# ------------------------------------------------------------

genome_cox_results_6X$log2HR <- NA_real_

genome_cox_results_6X$neglog10FDR <- NA_real_

genome_cox_results_6X$log2HR[
  ok_idx_6X
] <- log2(
  genome_cox_results_6X$HR[
    ok_idx_6X
  ]
)

genome_cox_results_6X$neglog10FDR[
  ok_idx_6X
] <- -log10(
  pmax(
    genome_cox_results_6X$FDR[
      ok_idx_6X
    ],
    .Machine$double.xmin
  )
)


# ------------------------------------------------------------
# 11. Map core genes
# ------------------------------------------------------------

core_gene_map_6X <- c(
  "ENSG00000131471" = "AOC3",
  "ENSG00000166927" = "MS4A7",
  "ENSG00000168685" = "IL7R"
)

genome_cox_results_6X$Gene_symbol <- unname(
  core_gene_map_6X[
    genome_cox_results_6X$Gene
  ]
)

genome_cox_results_6X$Core_gene <- ifelse(
  genome_cox_results_6X$Gene %in%
    names(core_gene_map_6X),
  "Core gene",
  "Other gene"
)


# ------------------------------------------------------------
# 12. Multiple-testing significance
# ------------------------------------------------------------

genome_cox_results_6X$FDR_significant <- FALSE

genome_cox_results_6X$FDR_significant[
  ok_idx_6X
] <- genome_cox_results_6X$FDR[
  ok_idx_6X
] < 0.05


# ------------------------------------------------------------
# 13. Basic result QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 6X-8C RESULT QC\n")
cat("============================================================\n")

cat(
  "Total genes =",
  nrow(genome_cox_results_6X),
  "\n"
)

cat(
  "Successful Cox models =",
  sum(genome_cox_results_6X$Status == "OK"),
  "\n"
)

cat(
  "Failed Cox models =",
  sum(genome_cox_results_6X$Status == "FAILED"),
  "\n"
)

cat(
  "Models with warnings =",
  sum(
    !is.na(
      genome_cox_results_6X$Warning
    )
  ),
  "\n"
)

cat(
  "Genome-wide FDR < 0.05 =",
  sum(
    genome_cox_results_6X$FDR_significant,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Finite HR =",
  sum(
    is.finite(
      genome_cox_results_6X$HR
    )
  ),
  "\n"
)

cat(
  "Finite FDR =",
  sum(
    is.finite(
      genome_cox_results_6X$FDR
    )
  ),
  "\n"
)


# ------------------------------------------------------------
# 14. Core-gene results
# ------------------------------------------------------------

cat("\nCore gene genome-wide survival results:\n")

core_result_6X <- genome_cox_results_6X[
  genome_cox_results_6X$Gene %in%
    names(core_gene_map_6X),
  c(
    "Gene",
    "Gene_symbol",
    "N",
    "Events",
    "SD",
    "HR",
    "CI_low",
    "CI_high",
    "P",
    "FDR",
    "log2HR",
    "neglog10FDR",
    "Status",
    "Warning"
  ),
  drop = FALSE
]

core_result_6X <- core_result_6X[
  match(
    names(core_gene_map_6X),
    core_result_6X$Gene
  ),
  ,
  drop = FALSE
]

print(core_result_6X)


# ------------------------------------------------------------
# 15. Compare genome-wide core results with locked
#     targeted Step 6B univariable results
# ------------------------------------------------------------

cat(
  "\nNote: core-gene estimates should be numerically close to\n",
  "the locked Step 6B univariable Cox results because both\n",
  "analyses use the same 485-observation cohort and continuous\n",
  "1-SD expression framework.\n"
)


# ------------------------------------------------------------
# 16. Save intermediate checkpoint
# ------------------------------------------------------------

step6X8C_checkpoint_dir <- file.path(
  "result_dir/Paper 1_Step 6-survival analysis",
  "STEP_6X_genome_wide_survival",
  "checkpoint"
)

dir.create(
  step6X8C_checkpoint_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


save(
  genome_cox_results_6X,
  core_result_6X,
  genome_survival_gene_universe_6X,
  file = file.path(
    step6X8C_checkpoint_dir,
    "Paper1_STEP6X8C_genome_wide_Cox_checkpoint.RData"
  )
)


# ------------------------------------------------------------
# 17. Final message
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 6X-8C COMPLETED\n")
cat("============================================================\n")

cat(
  "Genome-wide genes analyzed =",
  nrow(genome_cox_results_6X),
  "\n"
)

cat(
  "Successful models =",
  sum(genome_cox_results_6X$Status == "OK"),
  "\n"
)

cat(
  "Failed models =",
  sum(genome_cox_results_6X$Status == "FAILED"),
  "\n"
)

cat(
  "Warnings =",
  sum(!is.na(genome_cox_results_6X$Warning)),
  "\n"
)

cat(
  "Genome-wide FDR < 0.05 =",
  sum(
    genome_cox_results_6X$FDR_significant,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Checkpoint saved:\n",
  file.path(
    step6X8C_checkpoint_dir,
    "Paper1_STEP6X8C_genome_wide_Cox_checkpoint.RData"
  ),
  "\n"
)

cat("\n============================================================\n")

# ============================================================
# Paper 1
# STEP 6X-8D
# Genome-wide survival Cox result audit
# ============================================================

cat("\n============================================================\n")
cat("STEP 6X-8D: Genome-wide survival Cox result audit\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Basic result structure
# ------------------------------------------------------------

cat("\n[1] Result dimensions\n")

print(dim(genome_cox_results_6X))

cat("\nColumn names:\n")
print(colnames(genome_cox_results_6X))


# ------------------------------------------------------------
# 2. Status audit
# ------------------------------------------------------------

cat("\n[2] Model status\n")

print(table(
  genome_cox_results_6X$Status,
  useNA = "ifany"
))


# ------------------------------------------------------------
# 3. Warning audit
# ------------------------------------------------------------

cat("\n[3] Warning audit\n")

warning_idx_6X <- which(
  !is.na(genome_cox_results_6X$Warning) &
    genome_cox_results_6X$Warning != ""
)

cat(
  "Number of genes with warnings =",
  length(warning_idx_6X),
  "\n"
)

if (length(warning_idx_6X) > 0) {
  
  warning_table_6X <- genome_cox_results_6X[
    warning_idx_6X,
    c(
      "Gene",
      "Beta",
      "HR",
      "CI_low",
      "CI_high",
      "P",
      "FDR",
      "Warning"
    ),
    drop = FALSE
  ]
  
  print(
    warning_table_6X,
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 4. Numerical validity audit
# ------------------------------------------------------------

cat("\n[4] Numerical validity audit\n")

cat(
  "NA in HR =",
  sum(is.na(genome_cox_results_6X$HR)),
  "\n"
)

cat(
  "NA in P =",
  sum(is.na(genome_cox_results_6X$P)),
  "\n"
)

cat(
  "NA in FDR =",
  sum(is.na(genome_cox_results_6X$FDR)),
  "\n"
)

cat(
  "Non-finite HR =",
  sum(!is.finite(genome_cox_results_6X$HR)),
  "\n"
)

cat(
  "Non-finite P =",
  sum(!is.finite(genome_cox_results_6X$P)),
  "\n"
)

cat(
  "Non-finite FDR =",
  sum(!is.finite(genome_cox_results_6X$FDR)),
  "\n"
)


# ------------------------------------------------------------
# 5. FDR distribution
# ------------------------------------------------------------

cat("\n[5] Genome-wide FDR distribution\n")

print(
  summary(genome_cox_results_6X$FDR)
)

cat(
  "FDR < 0.05 =",
  sum(
    genome_cox_results_6X$FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "FDR < 0.01 =",
  sum(
    genome_cox_results_6X$FDR < 0.01,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "FDR < 0.001 =",
  sum(
    genome_cox_results_6X$FDR < 0.001,
    na.rm = TRUE
  ),
  "\n"
)


# ------------------------------------------------------------
# 6. Direction of association
# ------------------------------------------------------------

cat("\n[6] Direction of genome-wide associations\n")

direction_table_6X <- with(
  genome_cox_results_6X,
  table(
    ifelse(
      HR < 1,
      "HR < 1",
      ifelse(HR > 1, "HR > 1", "HR = 1")
    )
  )
)

print(direction_table_6X)


# ------------------------------------------------------------
# 7. Significant genes by direction
# ------------------------------------------------------------

sig_6X <- subset(
  genome_cox_results_6X,
  FDR < 0.05
)

cat("\n[7] Significant genes by direction\n")

cat(
  "FDR < 0.05 and HR < 1 =",
  sum(sig_6X$HR < 1, na.rm = TRUE),
  "\n"
)

cat(
  "FDR < 0.05 and HR > 1 =",
  sum(sig_6X$HR > 1, na.rm = TRUE),
  "\n"
)


# ------------------------------------------------------------
# 8. Top genome-wide associations
# ------------------------------------------------------------

cat("\n[8] Top 20 genes by FDR\n")

top20_FDR_6X <- genome_cox_results_6X[
  order(
    genome_cox_results_6X$FDR,
    genome_cox_results_6X$P
  ),
  c(
    "Gene",
    "Beta",
    "HR",
    "CI_low",
    "CI_high",
    "P",
    "FDR"
  )
]

print(
  head(top20_FDR_6X, 20),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Top genes with strongest inverse HR
# ------------------------------------------------------------

cat("\n[9] Strongest inverse associations\n")

top_inverse_6X <- genome_cox_results_6X[
  order(genome_cox_results_6X$HR),
  c(
    "Gene",
    "Beta",
    "HR",
    "CI_low",
    "CI_high",
    "P",
    "FDR"
  )
]

print(
  head(top_inverse_6X, 20),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Top genes with strongest positive HR
# ------------------------------------------------------------

cat("\n[10] Strongest positive associations\n")

top_positive_6X <- genome_cox_results_6X[
  order(
    genome_cox_results_6X$HR,
    decreasing = TRUE
  ),
  c(
    "Gene",
    "Beta",
    "HR",
    "CI_low",
    "CI_high",
    "P",
    "FDR"
  )
]

print(
  head(top_positive_6X, 20),
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Core gene mapping
# ------------------------------------------------------------

cat("\n[11] Core gene audit\n")

core_gene_map_6X <- c(
  "ENSG00000131471" = "AOC3",
  "ENSG00000166927" = "MS4A7",
  "ENSG00000168685" = "IL7R"
)

core_idx_6X <- match(
  names(core_gene_map_6X),
  genome_cox_results_6X$Gene
)

core_results_6X <- genome_cox_results_6X[
  core_idx_6X,
  c(
    "Gene",
    "N",
    "Events",
    "Mean",
    "SD",
    "Beta",
    "HR",
    "CI_low",
    "CI_high",
    "P",
    "FDR",
    "log2HR",
    "neglog10FDR"
  ),
  drop = FALSE
]

core_results_6X$Gene_symbol <- unname(
  core_gene_map_6X[
    as.character(core_results_6X$Gene)
  ]
)

core_results_6X <- core_results_6X[
  c(
    "Gene_symbol",
    "Gene",
    "N",
    "Events",
    "Mean",
    "SD",
    "Beta",
    "HR",
    "CI_low",
    "CI_high",
    "P",
    "FDR",
    "log2HR",
    "neglog10FDR"
  )
]

print(
  core_results_6X,
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Compare core genes with locked Step 6B results
# ------------------------------------------------------------

cat("\n[12] Core-gene consistency check\n")

core_locked_expected_6X <- data.frame(
  Gene_symbol = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  HR_locked = c(
    0.8694600,
    0.8530846,
    0.8469172
  ),
  P_locked = c(
    0.04829482,
    0.02127990,
    0.01986458
  )
)

core_compare_6X <- merge(
  core_results_6X[
    c("Gene_symbol", "HR", "P", "FDR")
  ],
  core_locked_expected_6X,
  by = "Gene_symbol",
  sort = FALSE
)

core_compare_6X$HR_difference <-
  core_compare_6X$HR -
  core_compare_6X$HR_locked

core_compare_6X$P_difference <-
  core_compare_6X$P -
  core_compare_6X$P_locked

print(
  core_compare_6X,
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Save audit object
# ------------------------------------------------------------

step6X8D_audit <- list(
  warning_table = warning_table_6X,
  top20_FDR = head(top20_FDR_6X, 20),
  top_inverse = head(top_inverse_6X, 20),
  top_positive = head(top_positive_6X, 20),
  core_results = core_results_6X,
  core_comparison = core_compare_6X
)

save(
  step6X8D_audit,
  file = file.path(
    step6X8D_audit_dir,
    "Paper1_STEP6X8D_genome_wide_Cox_audit.RData"
  )
)

cat("\n============================================================\n")
cat("STEP 6X-8D AUDIT COMPLETED\n")
cat("============================================================\n")


### ============================================================
# STEP 6X-8D
# Define audit directory and save completed audit object
# ============================================================

step6X8D_base_dir <- file.path(
  "result_dir/Paper 1_Step 6-survival analysis",
  "STEP_6X_genome_wide_survival"
)

step6X8D_audit_dir <- file.path(
  step6X8D_base_dir,
  "audit"
)

if (!dir.exists(step6X8D_audit_dir)) {
  dir.create(
    step6X8D_audit_dir,
    recursive = TRUE
  )
}

save(
  step6X8D_audit,
  file = file.path(
    step6X8D_audit_dir,
    "Paper1_STEP6X8D_genome_wide_Cox_audit.RData"
  )
)

cat("\n============================================================\n")
cat("STEP 6X-8D AUDIT SAVED\n")
cat("============================================================\n")

cat(
  "Audit file:\n",
  file.path(
    step6X8D_audit_dir,
    "Paper1_STEP6X8D_genome_wide_Cox_audit.RData"
  ),
  "\n"
)

cat(
  "File exists =",
  file.exists(
    file.path(
      step6X8D_audit_dir,
      "Paper1_STEP6X8D_genome_wide_Cox_audit.RData"
    )
  ),
  "\n"
)

# ============================================================
# Paper 1
# STEP 6X-8E-R2
# Recover and audit Figure X source objects
# WITHOUT rerunning genome-wide Cox analysis
# ============================================================

cat("\n============================================================\n")
cat("STEP 6X-8E-R2: Recover and audit Figure X source objects\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Define directories
# ------------------------------------------------------------

step6X_base_dir <- file.path(
  "input_dir/Paper 1_Step 6-survival analysis",
  "STEP_6X_genome_wide_survival"
)

step6X_checkpoint_dir <- file.path(
  step6X_base_dir,
  "checkpoint"
)

step6X_audit_dir <- file.path(
  step6X_base_dir,
  "audit"
)


step6X8C_checkpoint_file <- file.path(
  step6X_checkpoint_dir,
  "Paper1_STEP6X8C_genome_wide_Cox_checkpoint.RData"
)

step6X8D_audit_file <- file.path(
  step6X_audit_dir,
  "Paper1_STEP6X8D_genome_wide_Cox_audit.RData"
)


# ------------------------------------------------------------
# 2. Load 8C checkpoint into a temporary environment
# ------------------------------------------------------------

cat("\n[1] Loading 6X-8C checkpoint into temporary environment\n")

env_6X8C <- new.env(parent = emptyenv())

checkpoint_objects_6X8E_R2 <- load(
  step6X8C_checkpoint_file,
  envir = env_6X8C
)

print(checkpoint_objects_6X8E_R2)


# ------------------------------------------------------------
# 3. Inspect exact object classes and dimensions
# ------------------------------------------------------------

cat("\n[2] Objects recovered from 6X-8C\n")

for (obj_name in checkpoint_objects_6X8E_R2) {
  
  obj <- get(obj_name, envir = env_6X8C)
  
  cat("\nObject:", obj_name, "\n")
  cat("Class :", paste(class(obj), collapse = ", "), "\n")
  
  if (is.data.frame(obj) || is.matrix(obj)) {
    cat(
      "Dim   :",
      paste(dim(obj), collapse = " x "),
      "\n"
    )
  }
  
  if (is.vector(obj) && !is.list(obj)) {
    cat(
      "Length:",
      length(obj),
      "\n"
    )
  }
}


# ------------------------------------------------------------
# 4. Recover genome-wide Cox result
# ------------------------------------------------------------

if (!exists(
  "genome_cox_results_6X",
  envir = env_6X8C,
  inherits = FALSE
)) {
  stop(
    "genome_cox_results_6X is not present in the 6X-8C checkpoint."
  )
}

genome_cox_results_6X <- get(
  "genome_cox_results_6X",
  envir = env_6X8C
)

cat("\n[3] Genome-wide Cox result recovered\n")

cat(
  "Dimensions:",
  paste(dim(genome_cox_results_6X), collapse = " x "),
  "\n"
)

cat(
  "Columns:",
  paste(
    colnames(genome_cox_results_6X),
    collapse = ", "
  ),
  "\n"
)


# ------------------------------------------------------------
# 5. Recover genome survival gene universe
# ------------------------------------------------------------

if (!exists(
  "genome_survival_gene_universe_6X",
  envir = env_6X8C,
  inherits = FALSE
)) {
  stop(
    "genome_survival_gene_universe_6X is not present."
  )
}

genome_survival_gene_universe_6X <- get(
  "genome_survival_gene_universe_6X",
  envir = env_6X8C
)

cat("\n[4] Genome survival gene universe\n")

cat(
  "Length:",
  length(genome_survival_gene_universe_6X),
  "\n"
)

cat(
  "Unique:",
  length(unique(genome_survival_gene_universe_6X)),
  "\n"
)


# ------------------------------------------------------------
# 6. Verify genome-wide result integrity
# ------------------------------------------------------------

cat("\n[5] Genome-wide result integrity\n")

cat(
  "Genes analyzed =",
  nrow(genome_cox_results_6X),
  "\n"
)

cat(
  "All Status == OK =",
  all(genome_cox_results_6X$Status == "OK"),
  "\n"
)

cat(
  "NA in HR =",
  sum(is.na(genome_cox_results_6X$HR)),
  "\n"
)

cat(
  "NA in P =",
  sum(is.na(genome_cox_results_6X$P)),
  "\n"
)

cat(
  "NA in FDR =",
  sum(is.na(genome_cox_results_6X$FDR)),
  "\n"
)

cat(
  "FDR < 0.05 =",
  sum(genome_cox_results_6X$FDR < 0.05),
  "\n"
)


# ------------------------------------------------------------
# 7. Verify survival cohort information directly
# ------------------------------------------------------------

cat("\n[6] Survival cohort information encoded in Cox results\n")

cat(
  "Unique N values =",
  paste(
    sort(unique(genome_cox_results_6X$N)),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Unique Events values =",
  paste(
    sort(unique(genome_cox_results_6X$Events)),
    collapse = ", "
  ),
  "\n"
)

if (!all(genome_cox_results_6X$N == 485)) {
  stop("Genome-wide Cox result does not consistently report N = 485.")
}

if (!all(genome_cox_results_6X$Events == 178)) {
  stop("Genome-wide Cox result does not consistently report Events = 178.")
}

cat(
  "Locked survival N = 485: PASSED\n"
)

cat(
  "Locked survival Events = 178: PASSED\n"
)

cat(
  "Locked survival Censored =",
  485 - 178,
  "\n"
)


# ------------------------------------------------------------
# 8. Verify gene universe against Cox result
# ------------------------------------------------------------

cat("\n[7] Gene-universe consistency\n")

cox_genes_6X <- as.character(
  genome_cox_results_6X$Gene
)

universe_genes_6X <- as.character(
  genome_survival_gene_universe_6X
)

cat(
  "Cox genes:",
  length(cox_genes_6X),
  "\n"
)

cat(
  "Universe genes:",
  length(universe_genes_6X),
  "\n"
)

cat(
  "Cox genes unique =",
  length(unique(cox_genes_6X)) == length(cox_genes_6X),
  "\n"
)

cat(
  "Universe unique =",
  length(unique(universe_genes_6X)) == length(universe_genes_6X),
  "\n"
)

cat(
  "Same gene universe =",
  setequal(cox_genes_6X, universe_genes_6X),
  "\n"
)

if (!setequal(cox_genes_6X, universe_genes_6X)) {
  
  cat("\nGenes in Cox result but not universe:\n")
  
  print(
    setdiff(
      cox_genes_6X,
      universe_genes_6X
    )
  )
  
  cat("\nGenes in universe but not Cox result:\n")
  
  print(
    setdiff(
      universe_genes_6X,
      cox_genes_6X
    )
  )
  
  stop(
    "Genome-wide Cox gene universe does not match saved genome survival universe."
  )
}


# ------------------------------------------------------------
# 9. Core gene audit
# ------------------------------------------------------------

core_genes_6X <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_ensembl_6X <- c(
  AOC3  = "ENSG00000131471",
  MS4A7 = "ENSG00000166927",
  IL7R  = "ENSG00000168685"
)

cat("\n[8] Core-gene audit\n")

core_result_6X_R2 <- genome_cox_results_6X[
  genome_cox_results_6X$Gene %in% core_ensembl_6X,
  ,
  drop = FALSE
]

core_result_6X_R2$Gene_symbol <- names(
  core_ensembl_6X[
    match(
      core_result_6X_R2$Gene,
      core_ensembl_6X
    )
  ]
)

core_result_6X_R2 <- core_result_6X_R2[
  match(
    core_genes_6X,
    core_result_6X_R2$Gene_symbol
  ),
  ,
  drop = FALSE
]

print(
  core_result_6X_R2[
    ,
    c(
      "Gene",
      "Gene_symbol",
      "N",
      "Events",
      "HR",
      "CI_low",
      "CI_high",
      "P",
      "FDR",
      "log2HR",
      "neglog10FDR"
    )
  ]
)


# ------------------------------------------------------------
# 10. Confirm exact core-gene values
# ------------------------------------------------------------

expected_core_HR_6X <- c(
  AOC3  = 0.8694600,
  MS4A7 = 0.8530846,
  IL7R  = 0.8469172
)

expected_core_P_6X <- c(
  AOC3  = 0.04829482,
  MS4A7 = 0.02127990,
  IL7R  = 0.01986458
)

for (g in core_genes_6X) {
  
  observed_hr <- core_result_6X_R2$HR[
    core_result_6X_R2$Gene_symbol == g
  ]
  
  observed_p <- core_result_6X_R2$P[
    core_result_6X_R2$Gene_symbol == g
  ]
  
  cat(
    g,
    ": HR difference =",
    observed_hr - expected_core_HR_6X[g],
    "; P difference =",
    observed_p - expected_core_P_6X[g],
    "\n"
  )
}


# ------------------------------------------------------------
# 11. Load 8D audit separately if available
# ------------------------------------------------------------

cat("\n[9] 8D audit availability\n")

if (file.exists(step6X8D_audit_file)) {
  
  env_6X8D <- new.env(parent = emptyenv())
  
  audit_objects_6X8E_R2 <- load(
    step6X8D_audit_file,
    envir = env_6X8D
  )
  
  cat(
    "8D audit objects:\n"
  )
  
  print(audit_objects_6X8E_R2)
  
} else {
  
  cat(
    "8D audit file not found.\n"
  )
}


# ------------------------------------------------------------
# 12. Final recovery conclusion
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 6X-8E-R2: RECOVERY AUDIT COMPLETED\n")
cat("============================================================\n")

cat(
  "Genome-wide Cox results available = TRUE\n"
)

cat(
  "Genome survival gene universe available = TRUE\n"
)

cat(
  "Genome-wide genes =",
  nrow(genome_cox_results_6X),
  "\n"
)

cat(
  "Survival N = 485\n"
)

cat(
  "Events = 178\n"
)

cat(
  "Censored = 307\n"
)

cat(
  "FDR < 0.05 =",
  sum(genome_cox_results_6X$FDR < 0.05),
  "\n"
)

cat(
  "Core genes recovered =",
  all(core_genes_6X %in% core_result_6X_R2$Gene_symbol),
  "\n"
)

cat(
  "No genome-wide Cox rerun performed.\n"
)

cat(
  "No locked 6X-8C/8D results modified.\n"
)

cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 6X-8E-R3
# Parse saved genome-survival universe metadata
# NO Cox rerun
# NO expression reload
# NO modification of locked results
# ============================================================

cat("\n============================================================\n")
cat("STEP 6X-8E-R3: Parse saved genome-survival universe metadata\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Check object
# ------------------------------------------------------------

cat("\n[1] Object structure\n")

cat(
  "Class =",
  paste(
    class(genome_survival_gene_universe_6X),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Length =",
  length(genome_survival_gene_universe_6X),
  "\n"
)

cat("\nNames:\n")
print(names(genome_survival_gene_universe_6X))


# ------------------------------------------------------------
# 2. Inspect every list element safely
# ------------------------------------------------------------

cat("\n[2] Element-level structure\n")

for (i in seq_along(genome_survival_gene_universe_6X)) {
  
  x <- genome_survival_gene_universe_6X[[i]]
  
  cat(
    "\nElement",
    i,
    "\n"
  )
  
  cat(
    "Class:",
    paste(class(x), collapse = ", "),
    "\n"
  )
  
  cat(
    "Length:",
    length(x),
    "\n"
  )
  
  if (length(x) <= 10) {
    
    print(x)
    
  } else {
    
    cat(
      "First 5 values:\n"
    )
    
    print(
      head(x, 5)
    )
    
    cat(
      "Last 5 values:\n"
    )
    
    print(
      tail(x, 5)
    )
  }
}


# ------------------------------------------------------------
# 3. Identify the actual 31,609-gene vector
# ------------------------------------------------------------

gene_vector_candidates_6X <- which(
  vapply(
    genome_survival_gene_universe_6X,
    function(x) {
      
      is.character(x) &&
        length(x) == 31609 &&
        all(
          grepl(
            "^ENSG[0-9]+$",
            x
          )
        )
      
    },
    logical(1)
  )
)

cat(
  "\n[3] Candidate 31,609-gene vector index(es):",
  paste(
    gene_vector_candidates_6X,
    collapse = ", "
  ),
  "\n"
)

if (length(gene_vector_candidates_6X) != 1) {
  
  stop(
    paste(
      "Expected exactly one 31,609-gene Ensembl vector, found",
      length(gene_vector_candidates_6X)
    )
  )
}

# ============================================================
# STEP 6X-8E-R3
# Continue from Section 4
# ============================================================

cat("\n============================================================\n")
cat("STEP 6X-8E-R3: Continue from Section 4\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 4. Extract actual gene universe
# ------------------------------------------------------------

genome_gene_universe_vector_6X <-
  genome_survival_gene_universe_6X[[gene_vector_candidates_6X]]

cat(
  "\n[4] Extracted genome gene universe\n"
)

cat(
  "Class =",
  class(genome_gene_universe_vector_6X),
  "\n"
)

cat(
  "Length =",
  length(genome_gene_universe_vector_6X),
  "\n"
)

cat(
  "Unique =",
  length(unique(genome_gene_universe_vector_6X)),
  "\n"
)

cat(
  "All Ensembl format =",
  all(
    grepl(
      "^ENSG[0-9]+$",
      genome_gene_universe_vector_6X
    )
  ),
  "\n"
)


# ------------------------------------------------------------
# 5. Compare with genome-wide Cox result
# ------------------------------------------------------------

cox_genes_6X <- as.character(
  genome_cox_results_6X$Gene
)

cat("\n[5] Compare gene universe with Cox results\n")

cat(
  "Cox genes =",
  length(cox_genes_6X),
  "\n"
)

cat(
  "Cox genes unique =",
  length(unique(cox_genes_6X)),
  "\n"
)

cat(
  "Universe genes =",
  length(genome_gene_universe_vector_6X),
  "\n"
)

cat(
  "Universe genes unique =",
  length(unique(genome_gene_universe_vector_6X)),
  "\n"
)

cat(
  "identical order =",
  identical(
    cox_genes_6X,
    genome_gene_universe_vector_6X
  ),
  "\n"
)

cat(
  "same gene set =",
  setequal(
    cox_genes_6X,
    genome_gene_universe_vector_6X
  ),
  "\n"
)


# ------------------------------------------------------------
# 6. Explicit set-difference audit
# ------------------------------------------------------------

cox_not_in_universe_6X <- setdiff(
  cox_genes_6X,
  genome_gene_universe_vector_6X
)

universe_not_in_cox_6X <- setdiff(
  genome_gene_universe_vector_6X,
  cox_genes_6X
)

cat(
  "\nGenes in Cox results but not in universe =",
  length(cox_not_in_universe_6X),
  "\n"
)

cat(
  "Genes in universe but not in Cox results =",
  length(universe_not_in_cox_6X),
  "\n"
)


# ------------------------------------------------------------
# 7. Core-gene audit
# ------------------------------------------------------------

core_genes_6X <- c(
  AOC3  = "ENSG00000131471",
  MS4A7 = "ENSG00000166927",
  IL7R  = "ENSG00000168685"
)

cat("\n[7] Core-gene audit\n")

print(
  data.frame(
    Gene = names(core_genes_6X),
    Ensembl = unname(core_genes_6X),
    Present_in_universe =
      unname(core_genes_6X) %in%
      genome_gene_universe_vector_6X,
    Present_in_Cox =
      unname(core_genes_6X) %in%
      cox_genes_6X,
    row.names = NULL
  )
)


# ------------------------------------------------------------
# 8. Check saved metadata against Cox results
# ------------------------------------------------------------

cat("\n[8] Metadata consistency audit\n")

cat(
  "Metadata n_genes =",
  genome_survival_gene_universe_6X$n_genes,
  "\n"
)

cat(
  "Metadata n_samples =",
  genome_survival_gene_universe_6X$n_samples,
  "\n"
)

cat(
  "Metadata survival_N =",
  genome_survival_gene_universe_6X$survival_N,
  "\n"
)

cat(
  "Metadata survival_events =",
  genome_survival_gene_universe_6X$survival_events,
  "\n"
)

cat(
  "Metadata survival_censored =",
  genome_survival_gene_universe_6X$survival_censored,
  "\n"
)

cat(
  "Metadata wgcna_gene_source =",
  genome_survival_gene_universe_6X$wgcna_gene_source,
  "\n"
)

cat(
  "Metadata sample_order_locked =",
  genome_survival_gene_universe_6X$sample_order_locked,
  "\n"
)

cat(
  "Metadata gene_order_locked =",
  genome_survival_gene_universe_6X$gene_order_locked,
  "\n"
)


# ------------------------------------------------------------
# 9. Final PASS conditions
# ------------------------------------------------------------

pass_gene_number_6X <-
  length(genome_gene_universe_vector_6X) == 31609

pass_gene_unique_6X <-
  length(unique(genome_gene_universe_vector_6X)) == 31609

pass_gene_format_6X <-
  all(
    grepl(
      "^ENSG[0-9]+$",
      genome_gene_universe_vector_6X
    )
  )

pass_same_gene_set_6X <-
  setequal(
    cox_genes_6X,
    genome_gene_universe_vector_6X
  )

pass_core_genes_6X <-
  all(
    unname(core_genes_6X) %in%
      genome_gene_universe_vector_6X
  )

pass_core_cox_6X <-
  all(
    unname(core_genes_6X) %in%
      cox_genes_6X
  )

pass_metadata_genes_6X <-
  identical(
    as.integer(
      genome_survival_gene_universe_6X$n_genes
    ),
    31609L
  )

pass_metadata_samples_6X <-
  identical(
    as.integer(
      genome_survival_gene_universe_6X$n_samples
    ),
    485L
  )

pass_metadata_events_6X <-
  identical(
    as.integer(
      genome_survival_gene_universe_6X$survival_events
    ),
    178L
  )

pass_metadata_censored_6X <-
  identical(
    as.integer(
      genome_survival_gene_universe_6X$survival_censored
    ),
    307L
  )


# ------------------------------------------------------------
# 10. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 6X-8E-R3 FINAL QC\n")
cat("============================================================\n")

cat(
  "Gene number =",
  pass_gene_number_6X,
  "\n"
)

cat(
  "Gene unique =",
  pass_gene_unique_6X,
  "\n"
)

cat(
  "Gene format =",
  pass_gene_format_6X,
  "\n"
)

cat(
  "Same gene set as Cox =",
  pass_same_gene_set_6X,
  "\n"
)

cat(
  "Core genes in universe =",
  pass_core_genes_6X,
  "\n"
)

cat(
  "Core genes in Cox =",
  pass_core_cox_6X,
  "\n"
)

cat(
  "Metadata n_genes =",
  pass_metadata_genes_6X,
  "\n"
)

cat(
  "Metadata n_samples =",
  pass_metadata_samples_6X,
  "\n"
)

cat(
  "Metadata events =",
  pass_metadata_events_6X,
  "\n"
)

cat(
  "Metadata censored =",
  pass_metadata_censored_6X,
  "\n"
)


all_pass_6X <- all(
  pass_gene_number_6X,
  pass_gene_unique_6X,
  pass_gene_format_6X,
  pass_same_gene_set_6X,
  pass_core_genes_6X,
  pass_core_cox_6X,
  pass_metadata_genes_6X,
  pass_metadata_samples_6X,
  pass_metadata_events_6X,
  pass_metadata_censored_6X
)

cat("\n")

if (all_pass_6X) {
  
  cat("============================================================\n")
  cat("STEP 6X-8E-R3: PASSED\n")
  cat("============================================================\n")
  
} else {
  
  cat("============================================================\n")
  cat("STEP 6X-8E-R3: FAILED\n")
  cat("============================================================\n")
  
  stop(
    "STEP 6X-8E-R3 failed. Do not proceed to Figure X freeze."
  )
}


# ============================================================
# Paper 1
# STEP 6X-8E-final
# Figure X genome-wide survival data freeze
# ============================================================

cat("\n============================================================\n")
cat("STEP 6X-8E-final: Figure X data freeze\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 0. Output directory
# ------------------------------------------------------------

base_dir_6X <- 
  "input_dir/Paper 1_Step 6-survival analysis/STEP_6X_genome_wide_survival"

checkpoint_dir_6X <- 
  file.path(base_dir_6X, "checkpoint")

if (!dir.exists(checkpoint_dir_6X)) {
  dir.create(
    checkpoint_dir_6X,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ------------------------------------------------------------
# 1. Statistical source: DO NOT MODIFY
# ------------------------------------------------------------

stopifnot(
  exists("genome_cox_results_6X")
)

stopifnot(
  is.data.frame(genome_cox_results_6X)
)

figureX_source_6X <- genome_cox_results_6X


cat("\n[1] Statistical source\n")

cat(
  "Rows =",
  nrow(figureX_source_6X),
  "\n"
)

cat(
  "Columns =",
  ncol(figureX_source_6X),
  "\n"
)


# ------------------------------------------------------------
# 2. Gene universe: use R3-verified vector
# ------------------------------------------------------------

stopifnot(
  exists("genome_gene_universe_vector_6X")
)

stopifnot(
  length(genome_gene_universe_vector_6X) == 31609
)

stopifnot(
  length(unique(genome_gene_universe_vector_6X)) == 31609
)

stopifnot(
  setequal(
    figureX_source_6X$Gene,
    genome_gene_universe_vector_6X
  )
)


# ------------------------------------------------------------
# 3. Create plotting copy
# ------------------------------------------------------------

figureX_data_6X <- figureX_source_6X


# ------------------------------------------------------------
# 4. Add plotting variables only
# ------------------------------------------------------------

figureX_data_6X$Direction <- ifelse(
  figureX_data_6X$HR < 1,
  "HR < 1",
  ifelse(
    figureX_data_6X$HR > 1,
    "HR > 1",
    "HR = 1"
  )
)


figureX_data_6X$Plot_class <- "Non-significant"


figureX_data_6X$Plot_class[
  figureX_data_6X$FDR < 0.05 &
    figureX_data_6X$HR < 1
] <- "Significant_favorable_direction"


figureX_data_6X$Plot_class[
  figureX_data_6X$FDR < 0.05 &
    figureX_data_6X$HR > 1
] <- "Significant_unfavorable_direction"


figureX_data_6X$Plot_class[
  figureX_data_6X$Core_gene == TRUE
] <- "Core_gene"


# ------------------------------------------------------------
# 5. Fixed Figure X thresholds
# ------------------------------------------------------------

figureX_FDR_threshold_6X <- 0.05

figureX_log2HR_reference_6X <- 0

figureX_neglog10FDR_threshold_6X <-
  -log10(
    figureX_FDR_threshold_6X
  )


# ------------------------------------------------------------
# 6. Metadata
# ------------------------------------------------------------

figureX_metadata_6X <- list(
  
  figure_title =
    "Genome-wide survival associations across 31,609 genes in LUAD",
  
  analysis_type =
    "Genome-wide univariable Cox proportional hazards analysis",
  
  gene_universe_n =
    31609L,
  
  survival_N =
    485L,
  
  survival_events =
    178L,
  
  survival_censored =
    307L,
  
  fdr_method =
    "Benjamini-Hochberg",
  
  fdr_threshold =
    figureX_FDR_threshold_6X,
  
  log2HR_reference =
    figureX_log2HR_reference_6X,
  
  neglog10FDR_threshold =
    figureX_neglog10FDR_threshold_6X,
  
  expression_source =
    genome_survival_gene_universe_6X$expression_source,
  
  wgcna_gene_source =
    genome_survival_gene_universe_6X$wgcna_gene_source,
  
  core_genes =
    c(
      "AOC3",
      "MS4A7",
      "IL7R"
    ),
  
  core_gene_ensembl =
    c(
      AOC3  = "ENSG00000131471",
      MS4A7 = "ENSG00000166927",
      IL7R  = "ENSG00000168685"
    ),
  
  source_object =
    "genome_cox_results_6X",
  
  source_checkpoint =
    "Paper1_STEP6X8C_genome_wide_Cox_checkpoint.RData"
)


# ------------------------------------------------------------
# 7. Statistical integrity checks
# ------------------------------------------------------------

cat("\n[7] Statistical integrity checks\n")


check_HR_6X <- identical(
  figureX_source_6X$HR,
  figureX_data_6X$HR
)

check_P_6X <- identical(
  figureX_source_6X$P,
  figureX_data_6X$P
)

check_FDR_6X <- identical(
  figureX_source_6X$FDR,
  figureX_data_6X$FDR
)

check_log2HR_6X <- identical(
  figureX_source_6X$log2HR,
  figureX_data_6X$log2HR
)

check_neglog10FDR_6X <- identical(
  figureX_source_6X$neglog10FDR,
  figureX_data_6X$neglog10FDR
)

check_gene_6X <- identical(
  figureX_source_6X$Gene,
  figureX_data_6X$Gene
)


cat(
  "HR unchanged =",
  check_HR_6X,
  "\n"
)

cat(
  "P unchanged =",
  check_P_6X,
  "\n"
)

cat(
  "FDR unchanged =",
  check_FDR_6X,
  "\n"
)

cat(
  "log2HR unchanged =",
  check_log2HR_6X,
  "\n"
)

cat(
  "neglog10FDR unchanged =",
  check_neglog10FDR_6X,
  "\n"
)

cat(
  "Gene order unchanged =",
  check_gene_6X,
  "\n"
)


# ------------------------------------------------------------
# 8. Figure-level counts
# ------------------------------------------------------------

n_genes_6X <- nrow(
  figureX_data_6X
)

n_sig_6X <- sum(
  figureX_data_6X$FDR < 0.05,
  na.rm = TRUE
)

n_sig_hr_lt1_6X <- sum(
  figureX_data_6X$FDR < 0.05 &
    figureX_data_6X$HR < 1,
  na.rm = TRUE
)

n_sig_hr_gt1_6X <- sum(
  figureX_data_6X$FDR < 0.05 &
    figureX_data_6X$HR > 1,
  na.rm = TRUE
)

n_core_6X <- sum(
  figureX_data_6X$Core_gene == TRUE,
  na.rm = TRUE
)


cat("\n[8] Figure-level counts\n")

cat(
  "Total genes =",
  n_genes_6X,
  "\n"
)

cat(
  "Genome-wide FDR < 0.05 =",
  n_sig_6X,
  "\n"
)

cat(
  "FDR < 0.05 & HR < 1 =",
  n_sig_hr_lt1_6X,
  "\n"
)

cat(
  "FDR < 0.05 & HR > 1 =",
  n_sig_hr_gt1_6X,
  "\n"
)

cat(
  "Core genes =",
  n_core_6X,
  "\n"
)


# ------------------------------------------------------------
# 9. Core-gene statistical status
# ------------------------------------------------------------

core_check_6X <- figureX_data_6X[
  figureX_data_6X$Core_gene == TRUE,
  c(
    "Gene",
    "Gene_symbol",
    "HR",
    "CI_low",
    "CI_high",
    "P",
    "FDR",
    "log2HR",
    "neglog10FDR",
    "FDR_significant",
    "Plot_class"
  )
]

cat("\n[9] Core-gene Figure X status\n")

print(
  core_check_6X,
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Final freeze checks
# ------------------------------------------------------------

all_freeze_checks_6X <- all(
  n_genes_6X == 31609,
  n_sig_6X == 1940,
  n_sig_hr_lt1_6X == 1318,
  n_sig_hr_gt1_6X == 622,
  n_core_6X == 3,
  check_HR_6X,
  check_P_6X,
  check_FDR_6X,
  check_log2HR_6X,
  check_neglog10FDR_6X,
  check_gene_6X,
  setequal(
    figureX_data_6X$Gene,
    genome_gene_universe_vector_6X
  )
)


# ------------------------------------------------------------
# 11. Save final Figure X checkpoint
# ------------------------------------------------------------

if (all_freeze_checks_6X) {
  
  save(
    figureX_source_6X,
    figureX_data_6X,
    figureX_metadata_6X,
    genome_gene_universe_vector_6X,
    figureX_FDR_threshold_6X,
    figureX_log2HR_reference_6X,
    figureX_neglog10FDR_threshold_6X,
    core_check_6X,
    n_genes_6X,
    n_sig_6X,
    n_sig_hr_lt1_6X,
    n_sig_hr_gt1_6X,
    n_core_6X,
    file = file.path(
      checkpoint_dir_6X,
      "Paper1_STEP6X8E_FigureX_data_freeze_checkpoint.RData"
    )
  )
  
  cat("\n============================================================\n")
  cat("STEP 6X-8E-final: PASSED\n")
  cat("Figure X data freeze completed successfully.\n")
  cat("============================================================\n")
  
} else {
  
  cat("\n============================================================\n")
  cat("STEP 6X-8E-final: FAILED\n")
  cat("Figure X data freeze was NOT saved.\n")
  cat("============================================================\n")
  
  stop(
    "Figure X freeze QC failed. Do not proceed to plotting."
  )
}


cat("\n============================================================\n")
cat("Diagnose Figure X Core_gene annotation\n")
cat("============================================================\n")

cat("\n[1] Core_gene class:\n")
print(class(figureX_source_6X$Core_gene))

cat("\n[2] Core_gene table:\n")
print(table(
  figureX_source_6X$Core_gene,
  useNA = "ifany"
))

cat("\n[3] Core_gene head:\n")
print(
  head(
    figureX_source_6X[, c(
      "Gene",
      "Gene_symbol",
      "Core_gene"
    )],
    10
  ),
  row.names = FALSE
)

cat("\n[4] Search core Ensembl IDs directly:\n")

core_ensembl_6X <- c(
  "ENSG00000131471",
  "ENSG00000166927",
  "ENSG00000168685"
)

print(
  figureX_source_6X[
    figureX_source_6X$Gene %in% core_ensembl_6X,
    c(
      "Gene",
      "Gene_symbol",
      "HR",
      "CI_low",
      "CI_high",
      "P",
      "FDR",
      "log2HR",
      "neglog10FDR",
      "Core_gene"
    )
  ],
  row.names = FALSE
)

### ============================================================
# STEP 6X-8E-final
# Correct Figure X core-gene annotation
# ============================================================

cat("\n============================================================\n")
cat("STEP 6X-8E-final: Correct core-gene annotation\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Preserve original statistical source
# ------------------------------------------------------------

figureX_source_6X <- genome_cox_results_6X


# ------------------------------------------------------------
# 2. Create plotting copy
# ------------------------------------------------------------

figureX_data_6X <- figureX_source_6X


# ------------------------------------------------------------
# 3. Define locked core-gene Ensembl IDs
# ------------------------------------------------------------

core_ensembl_6X <- c(
  "ENSG00000131471",  # AOC3
  "ENSG00000166927",  # MS4A7
  "ENSG00000168685"   # IL7R
)


# ------------------------------------------------------------
# 4. Rebuild Figure X-specific logical Core_gene flag
# ------------------------------------------------------------

figureX_data_6X$Core_gene <- 
  figureX_data_6X$Gene %in% core_ensembl_6X


# ------------------------------------------------------------
# 5. Rebuild Plot_class
# ------------------------------------------------------------

figureX_data_6X$Plot_class <- "Non-significant"


figureX_data_6X$Plot_class[
  figureX_data_6X$FDR < 0.05 &
    figureX_data_6X$HR < 1
] <- "Significant_favorable_direction"


figureX_data_6X$Plot_class[
  figureX_data_6X$FDR < 0.05 &
    figureX_data_6X$HR > 1
] <- "Significant_unfavorable_direction"


figureX_data_6X$Plot_class[
  figureX_data_6X$Core_gene
] <- "Core_gene"


# ------------------------------------------------------------
# 6. Core-gene count
# ------------------------------------------------------------

n_core_6X <- sum(
  figureX_data_6X$Core_gene,
  na.rm = TRUE
)


cat("\nCore genes identified =",
    n_core_6X,
    "\n")


# ------------------------------------------------------------
# 7. Core-gene audit
# ------------------------------------------------------------

core_check_6X <- figureX_data_6X[
  figureX_data_6X$Core_gene,
  c(
    "Gene",
    "Gene_symbol",
    "HR",
    "CI_low",
    "CI_high",
    "P",
    "FDR",
    "log2HR",
    "neglog10FDR",
    "FDR_significant",
    "Plot_class"
  )
]


cat("\nCore-gene Figure X status:\n")

print(
  core_check_6X,
  row.names = FALSE
)


# ============================================================
# 8. Figure-level counts
# ============================================================

n_genes_6X <- nrow(
  figureX_data_6X
)

n_sig_6X <- sum(
  figureX_data_6X$FDR < 0.05,
  na.rm = TRUE
)

n_sig_hr_lt1_6X <- sum(
  figureX_data_6X$FDR < 0.05 &
    figureX_data_6X$HR < 1,
  na.rm = TRUE
)

n_sig_hr_gt1_6X <- sum(
  figureX_data_6X$FDR < 0.05 &
    figureX_data_6X$HR > 1,
  na.rm = TRUE
)


cat("\n============================================================\n")
cat("Figure X counts\n")
cat("============================================================\n")

cat("Total genes =", n_genes_6X, "\n")
cat("Genome-wide FDR < 0.05 =", n_sig_6X, "\n")
cat("FDR < 0.05 & HR < 1 =", n_sig_hr_lt1_6X, "\n")
cat("FDR < 0.05 & HR > 1 =", n_sig_hr_gt1_6X, "\n")
cat("Core genes =", n_core_6X, "\n")


# ============================================================
# 9. Final Figure X freeze QC
# ============================================================

check_HR_6X <- identical(
  figureX_source_6X$HR,
  figureX_data_6X$HR
)

check_P_6X <- identical(
  figureX_source_6X$P,
  figureX_data_6X$P
)

check_FDR_6X <- identical(
  figureX_source_6X$FDR,
  figureX_data_6X$FDR
)

check_log2HR_6X <- identical(
  figureX_source_6X$log2HR,
  figureX_data_6X$log2HR
)

check_neglog10FDR_6X <- identical(
  figureX_source_6X$neglog10FDR,
  figureX_data_6X$neglog10FDR
)

check_gene_6X <- identical(
  figureX_source_6X$Gene,
  figureX_data_6X$Gene
)


all_freeze_checks_6X <- all(
  
  n_genes_6X == 31609,
  
  n_sig_6X == 1940,
  
  n_sig_hr_lt1_6X == 1318,
  
  n_sig_hr_gt1_6X == 622,
  
  n_core_6X == 3,
  
  check_HR_6X,
  
  check_P_6X,
  
  check_FDR_6X,
  
  check_log2HR_6X,
  
  check_neglog10FDR_6X,
  
  check_gene_6X,
  
  setequal(
    figureX_data_6X$Gene,
    genome_gene_universe_vector_6X
  ),
  
  all(
    core_ensembl_6X %in%
      figureX_data_6X$Gene
  )
  
)


cat("\n============================================================\n")
cat("STEP 6X-8E-final FINAL QC\n")
cat("============================================================\n")

cat(
  "31,609 genes =",
  n_genes_6X == 31609,
  "\n"
)

cat(
  "FDR < 0.05 =",
  n_sig_6X == 1940,
  "\n"
)

cat(
  "HR < 1 significant =",
  n_sig_hr_lt1_6X == 1318,
  "\n"
)

cat(
  "HR > 1 significant =",
  n_sig_hr_gt1_6X == 622,
  "\n"
)

cat(
  "Core genes = 3 =",
  n_core_6X == 3,
  "\n"
)

cat(
  "HR unchanged =",
  check_HR_6X,
  "\n"
)

cat(
  "P unchanged =",
  check_P_6X,
  "\n"
)

cat(
  "FDR unchanged =",
  check_FDR_6X,
  "\n"
)

cat(
  "log2HR unchanged =",
  check_log2HR_6X,
  "\n"
)

cat(
  "neglog10FDR unchanged =",
  check_neglog10FDR_6X,
  "\n"
)

cat(
  "Gene order unchanged =",
  check_gene_6X,
  "\n"
)

cat(
  "Gene universe unchanged =",
  setequal(
    figureX_data_6X$Gene,
    genome_gene_universe_vector_6X
  ),
  "\n"
)


# ============================================================
# Paper 1
# STEP 6X-9A
# Figure X plotting distribution audit
# ============================================================
#
# Purpose:
#   Audit the frozen Figure X data before final plotting.
#
# IMPORTANT:
#   - No Cox model recalculation
#   - No P-value recalculation
#   - No FDR recalculation
#   - No modification of frozen analysis objects
#   - Read-only audit of figureX_data_6X
#
# Input:
#   figureX_data_6X
#
# Output:
#   Distribution / range / core-gene / reference-line audit
#
# ============================================================

cat("\n============================================================\n")
cat("STEP 6X-9A: Figure X plotting distribution audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Basic object check
# ------------------------------------------------------------

stopifnot(exists("figureX_data_6X"))
stopifnot(is.data.frame(figureX_data_6X))

required_cols_9A <- c(
  "Gene",
  "Gene_symbol",
  "HR",
  "FDR",
  "log2HR",
  "neglog10FDR",
  "Core_gene",
  "Plot_class"
)

missing_cols_9A <- setdiff(
  required_cols_9A,
  colnames(figureX_data_6X)
)

if (length(missing_cols_9A) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols_9A, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 1. Basic dimensions and numerical QC
# ------------------------------------------------------------

cat("1. BASIC DIMENSIONS\n")
cat("------------------------------------------------------------\n")

cat("Rows:", nrow(figureX_data_6X), "\n")
cat("Columns:", ncol(figureX_data_6X), "\n")

cat(
  "NA in log2HR:",
  sum(is.na(figureX_data_6X$log2HR)),
  "\n"
)

cat(
  "NA in neglog10FDR:",
  sum(is.na(figureX_data_6X$neglog10FDR)),
  "\n"
)

cat(
  "Non-finite log2HR:",
  sum(!is.finite(figureX_data_6X$log2HR)),
  "\n"
)

cat(
  "Non-finite neglog10FDR:",
  sum(!is.finite(figureX_data_6X$neglog10FDR)),
  "\n"
)

cat("\n")


# ------------------------------------------------------------
# 2. log2HR distribution
# ------------------------------------------------------------

cat("2. log2HR DISTRIBUTION\n")
cat("------------------------------------------------------------\n")

log2HR_summary_9A <- summary(
  figureX_data_6X$log2HR
)

print(log2HR_summary_9A)

log2HR_quantiles_9A <- quantile(
  figureX_data_6X$log2HR,
  probs = c(
    0,
    0.001,
    0.005,
    0.01,
    0.025,
    0.05,
    0.25,
    0.50,
    0.75,
    0.95,
    0.975,
    0.99,
    0.995,
    0.999,
    1
  ),
  na.rm = TRUE
)

cat("\nSelected quantiles:\n")
print(log2HR_quantiles_9A)

cat("\n")


# ------------------------------------------------------------
# 3. neglog10FDR distribution
# ------------------------------------------------------------

cat("3. -log10(FDR) DISTRIBUTION\n")
cat("------------------------------------------------------------\n")

neglog10FDR_summary_9A <- summary(
  figureX_data_6X$neglog10FDR
)

print(neglog10FDR_summary_9A)

neglog10FDR_quantiles_9A <- quantile(
  figureX_data_6X$neglog10FDR,
  probs = c(
    0,
    0.001,
    0.005,
    0.01,
    0.025,
    0.05,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.975,
    0.99,
    0.995,
    0.999,
    1
  ),
  na.rm = TRUE
)

cat("\nSelected quantiles:\n")
print(neglog10FDR_quantiles_9A)

cat("\n")


# ------------------------------------------------------------
# 4. Significant genes
# ------------------------------------------------------------

cat("4. GENOME-WIDE SIGNIFICANT GENE DISTRIBUTION\n")
cat("------------------------------------------------------------\n")

sig_fav_9A <- subset(
  figureX_data_6X,
  FDR < 0.05 & HR < 1
)

sig_unfav_9A <- subset(
  figureX_data_6X,
  FDR < 0.05 & HR > 1
)

cat(
  "FDR < 0.05:",
  sum(figureX_data_6X$FDR < 0.05),
  "\n"
)

cat(
  "Significant HR < 1:",
  nrow(sig_fav_9A),
  "\n"
)

cat(
  "Significant HR > 1:",
  nrow(sig_unfav_9A),
  "\n"
)

cat(
  "Significant HR < 1 and > 1 total:",
  nrow(sig_fav_9A) + nrow(sig_unfav_9A),
  "\n"
)

cat("\n")


# ------------------------------------------------------------
# 5. Significant gene x-axis distribution
# ------------------------------------------------------------

cat("5. SIGNIFICANT GENE x-AXIS DISTRIBUTION\n")
cat("------------------------------------------------------------\n")

cat("HR < 1 significant genes:\n")
print(
  summary(sig_fav_9A$log2HR)
)

cat("\nHR > 1 significant genes:\n")
print(
  summary(sig_unfav_9A$log2HR)
)

cat("\n")


# ------------------------------------------------------------
# 6. Significant gene y-axis distribution
# ------------------------------------------------------------

cat("6. SIGNIFICANT GENE y-AXIS DISTRIBUTION\n")
cat("------------------------------------------------------------\n")

cat("HR < 1 significant genes:\n")
print(
  summary(sig_fav_9A$neglog10FDR)
)

cat("\nHR > 1 significant genes:\n")
print(
  summary(sig_unfav_9A$neglog10FDR)
)

cat("\n")


# ------------------------------------------------------------
# 7. Extreme y-values
# ------------------------------------------------------------

cat("7. EXTREME -log10(FDR) VALUES\n")
cat("------------------------------------------------------------\n")

top_y_9A <- figureX_data_6X[
  order(
    figureX_data_6X$neglog10FDR,
    decreasing = TRUE
  ),
]

top_y_9A <- head(top_y_9A, 20)

print(
  top_y_9A[
    ,
    c(
      "Gene",
      "Gene_symbol",
      "HR",
      "FDR",
      "log2HR",
      "neglog10FDR"
    )
  ],
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 8. Core gene positions
# ------------------------------------------------------------

cat("8. CORE GENE POSITIONS\n")
cat("------------------------------------------------------------\n")

core_genes_9A <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_position_9A <- figureX_data_6X[
  figureX_data_6X$Gene_symbol %in% core_genes_9A,
  c(
    "Gene",
    "Gene_symbol",
    "HR",
    "FDR",
    "log2HR",
    "neglog10FDR",
    "Core_gene",
    "Plot_class"
  )
]

core_position_9A <- core_position_9A[
  match(
    core_genes_9A,
    core_position_9A$Gene_symbol
  ),
]

print(
  core_position_9A,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 9. Core gene distance from significance threshold
# ------------------------------------------------------------

cat("9. CORE GENE DISTANCE FROM FDR THRESHOLD\n")
cat("------------------------------------------------------------\n")

fdr_threshold_y_9A <- -log10(0.05)

core_position_9A$Y_distance_from_threshold <-
  core_position_9A$neglog10FDR -
  fdr_threshold_y_9A

core_position_9A$X_distance_from_HR1 <-
  core_position_9A$log2HR

print(
  core_position_9A[
    ,
    c(
      "Gene_symbol",
      "log2HR",
      "neglog10FDR",
      "Y_distance_from_threshold",
      "X_distance_from_HR1"
    )
  ],
  row.names = FALSE
)

cat(
  "\nFDR threshold y =",
  round(fdr_threshold_y_9A, 6),
  "\n\n"
)


# ------------------------------------------------------------
# 10. Reference-line checks
# ------------------------------------------------------------

cat("10. REFERENCE LINE CHECK\n")
cat("------------------------------------------------------------\n")

cat(
  "Vertical reference: log2(HR) = 0\n"
)

cat(
  "Equivalent HR = 1\n"
)

cat(
  "Horizontal reference: -log10(FDR) =",
  round(fdr_threshold_y_9A, 6),
  "\n"
)

cat(
  "Equivalent FDR = 0.05\n"
)

cat("\n")


# ------------------------------------------------------------
# 11. Candidate plotting limits
# ------------------------------------------------------------

cat("11. CANDIDATE PLOTTING LIMITS\n")
cat("------------------------------------------------------------\n")

x_min_9A <- min(
  figureX_data_6X$log2HR,
  na.rm = TRUE
)

x_max_9A <- max(
  figureX_data_6X$log2HR,
  na.rm = TRUE
)

y_min_9A <- min(
  figureX_data_6X$neglog10FDR,
  na.rm = TRUE
)

y_max_9A <- max(
  figureX_data_6X$neglog10FDR,
  na.rm = TRUE
)

cat(
  "Raw x range:",
  x_min_9A,
  "to",
  x_max_9A,
  "\n"
)

cat(
  "Raw y range:",
  y_min_9A,
  "to",
  y_max_9A,
  "\n"
)

cat("\n")


# ------------------------------------------------------------
# 12. Quantify extreme-y compression
# ------------------------------------------------------------

cat("12. EXTREME-y COMPRESSION AUDIT\n")
cat("------------------------------------------------------------\n")

compression_cutoffs_9A <- c(
  5,
  7.5,
  10,
  12.5,
  15
)

compression_table_9A <- data.frame(
  Y_cutoff = compression_cutoffs_9A,
  N_above = sapply(
    compression_cutoffs_9A,
    function(x)
      sum(
        figureX_data_6X$neglog10FDR > x,
        na.rm = TRUE
      )
  )
)

compression_table_9A$Percent_above <-
  100 *
  compression_table_9A$N_above /
  nrow(figureX_data_6X)

print(
  compression_table_9A,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 13. Core genes vs significant cloud
# ------------------------------------------------------------

cat("13. CORE GENE vs SIGNIFICANT CLOUD\n")
cat("------------------------------------------------------------\n")

core_position_9A$Within_genomewide_FDR05 <-
  core_position_9A$FDR < 0.05

core_position_9A$HR_direction <-
  ifelse(
    core_position_9A$HR < 1,
    "HR < 1",
    ifelse(
      core_position_9A$HR > 1,
      "HR > 1",
      "HR = 1"
    )
  )

print(
  core_position_9A[
    ,
    c(
      "Gene_symbol",
      "HR",
      "FDR",
      "Within_genomewide_FDR05",
      "HR_direction"
    )
  ],
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 14. Direct-label overlap audit
# ------------------------------------------------------------

cat("14. DIRECT LABEL AUDIT\n")
cat("------------------------------------------------------------\n")

core_x_range_9A <- range(
  core_position_9A$log2HR,
  na.rm = TRUE
)

core_y_range_9A <- range(
  core_position_9A$neglog10FDR,
  na.rm = TRUE
)

cat(
  "Core gene log2HR range:",
  paste(
    round(core_x_range_9A, 4),
    collapse = " to "
  ),
  "\n"
)

cat(
  "Core gene -log10(FDR) range:",
  paste(
    round(core_y_range_9A, 4),
    collapse = " to "
  ),
  "\n"
)

cat("\n")


# ------------------------------------------------------------
# 15. Plot-class audit
# ------------------------------------------------------------

cat("15. PLOT CLASS AUDIT\n")
cat("------------------------------------------------------------\n")

print(
  table(
    figureX_data_6X$Plot_class,
    useNA = "ifany"
  )
)

cat("\n")


# ------------------------------------------------------------
# 16. Core-gene plotting class audit
# ------------------------------------------------------------

cat("16. CORE GENE PLOT CLASS\n")
cat("------------------------------------------------------------\n")

print(
  figureX_data_6X[
    figureX_data_6X$Core_gene,
    c(
      "Gene_symbol",
      "HR",
      "FDR",
      "log2HR",
      "neglog10FDR",
      "Plot_class"
    )
  ],
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 17. Preliminary interpretation for plotting
# ------------------------------------------------------------

cat("17. PRELIMINARY PLOTTING INTERPRETATION\n")
cat("------------------------------------------------------------\n")

if (y_max_9A > 10) {
  cat(
    "NOTE: The y-axis contains high -log10(FDR) values.\n",
    "A direct full-range y-axis may visually compress the main point cloud.\n",
    sep = ""
  )
}

if (nrow(core_position_9A) == 3) {
  cat(
    "Core genes AOC3, MS4A7 and IL7R are all available for direct labeling.\n"
  )
}

if (
  all(
    !core_position_9A$Within_genomewide_FDR05
  )
) {
  cat(
    "All three core genes are below the genome-wide FDR < 0.05 threshold.\n"
  )
}

if (
  all(
    core_position_9A$HR < 1
  )
) {
  cat(
    "All three core genes have HR < 1 in the genome-wide univariable analysis.\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 18. Save audit object
# ------------------------------------------------------------

audit_6X9A <- list(
  input_object = "figureX_data_6X",
  n_genes = nrow(figureX_data_6X),
  log2HR_summary = log2HR_summary_9A,
  log2HR_quantiles = log2HR_quantiles_9A,
  neglog10FDR_summary = neglog10FDR_summary_9A,
  neglog10FDR_quantiles = neglog10FDR_quantiles_9A,
  significant_favorable = sig_fav_9A,
  significant_unfavorable = sig_unfav_9A,
  top_y = top_y_9A,
  core_position = core_position_9A,
  compression_table = compression_table_9A,
  fdr_threshold_y = fdr_threshold_y_9A
)

cat("Audit object created: audit_6X9A\n")

cat("\n============================================================\n")
cat("STEP 6X-9A: AUDIT COMPLETED\n")
cat("============================================================\n")

##STEP 6X-9B — Figure X final plotting-data preparation

# ============================================================
# Paper 1
# STEP 6X-9B
# Figure X final plotting-data preparation
# ============================================================
#
# Purpose:
#   Prepare a frozen plotting copy for Supplementary Figure X.
#
# IMPORTANT:
#   - No Cox recalculation
#   - No P-value recalculation
#   - No FDR recalculation
#   - No modification of figureX_data_6X
#   - No modification of genome-wide analysis objects
#
# Input:
#   figureX_data_6X
#
# Output:
#   figureX_plot_6X
#   figureX_core_6X
#   figureX_reference_6X
#
# Figure:
#   "Genome-wide survival associations across 31,609 genes in LUAD"
#
# ============================================================

cat("\n============================================================\n")
cat("STEP 6X-9B: Figure X final plotting-data preparation\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Input object check
# ------------------------------------------------------------

stopifnot(exists("figureX_data_6X"))
stopifnot(is.data.frame(figureX_data_6X))

required_cols_9B <- c(
  "Gene",
  "Gene_symbol",
  "HR",
  "FDR",
  "log2HR",
  "neglog10FDR",
  "Core_gene",
  "Plot_class"
)

missing_cols_9B <- setdiff(
  required_cols_9B,
  colnames(figureX_data_6X)
)

if (length(missing_cols_9B) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols_9B, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 1. Preserve original object integrity
# ------------------------------------------------------------

figureX_data_6X_nrow_before <- nrow(
  figureX_data_6X
)

figureX_data_6X_gene_order_before <- figureX_data_6X$Gene

figureX_data_6X_log2HR_before <- figureX_data_6X$log2HR

figureX_data_6X_FDR_before <- figureX_data_6X$FDR

figureX_data_6X_HR_before <- figureX_data_6X$HR


# ------------------------------------------------------------
# 2. Create plotting copy
# ------------------------------------------------------------

figureX_plot_6X <- figureX_data_6X


# ------------------------------------------------------------
# 3. Explicitly reconstruct plotting class
# ------------------------------------------------------------
#
# This is a plotting annotation only.
#
# Priority:
#   Core_gene
#       >
#   Significant_favorable_direction
#       >
#   Significant_unfavorable_direction
#       >
#   Non-significant
#
# Core genes remain a separate class even though
# all three have FDR > 0.05 in the genome-wide analysis.
#
# ------------------------------------------------------------

core_ensembl_6X <- c(
  "ENSG00000131471",  # AOC3
  "ENSG00000166927",  # MS4A7
  "ENSG00000168685"   # IL7R
)

figureX_plot_6X$Core_gene_plot <-
  figureX_plot_6X$Gene %in% core_ensembl_6X

figureX_plot_6X$Plot_class_final <-
  "Non-significant"

figureX_plot_6X$Plot_class_final[
  figureX_plot_6X$FDR < 0.05 &
    figureX_plot_6X$HR < 1 &
    !figureX_plot_6X$Core_gene_plot
] <- "Significant_favorable_direction"

figureX_plot_6X$Plot_class_final[
  figureX_plot_6X$FDR < 0.05 &
    figureX_plot_6X$HR > 1 &
    !figureX_plot_6X$Core_gene_plot
] <- "Significant_unfavorable_direction"

figureX_plot_6X$Plot_class_final[
  figureX_plot_6X$Core_gene_plot
] <- "Core_gene"


# ------------------------------------------------------------
# 4. Factor order
# ------------------------------------------------------------

figureX_plot_6X$Plot_class_final <- factor(
  figureX_plot_6X$Plot_class_final,
  levels = c(
    "Non-significant",
    "Significant_favorable_direction",
    "Significant_unfavorable_direction",
    "Core_gene"
  )
)


# ------------------------------------------------------------
# 5. Fixed plotting coordinates
# ------------------------------------------------------------

figureX_plot_6X$x <- figureX_plot_6X$log2HR

figureX_plot_6X$y <- figureX_plot_6X$neglog10FDR


# ------------------------------------------------------------
# 6. Fixed reference-line values
# ------------------------------------------------------------

figureX_reference_6X <- data.frame(
  reference = c(
    "HR = 1",
    "FDR = 0.05"
  ),
  xintercept = c(
    0,
    NA_real_
  ),
  yintercept = c(
    NA_real_,
    -log10(0.05)
  )
)

fdr_threshold_6X <- -log10(0.05)


# ------------------------------------------------------------
# 7. Fixed plotting limits
# ------------------------------------------------------------

figureX_limits_6X <- list(
  x = c(-0.75, 0.75),
  y = c(0, 4.60)
)


# ------------------------------------------------------------
# 8. Core-gene plotting data
# ------------------------------------------------------------

figureX_core_6X <- figureX_plot_6X[
  figureX_plot_6X$Core_gene_plot,
  c(
    "Gene",
    "Gene_symbol",
    "HR",
    "FDR",
    "log2HR",
    "neglog10FDR",
    "x",
    "y",
    "Core_gene_plot",
    "Plot_class_final"
  )
]

figureX_core_6X <- figureX_core_6X[
  match(
    c("AOC3", "MS4A7", "IL7R"),
    figureX_core_6X$Gene_symbol
  ),
]


# ------------------------------------------------------------
# 9. Fixed label order
# ------------------------------------------------------------

figureX_core_6X$Gene_symbol <- factor(
  figureX_core_6X$Gene_symbol,
  levels = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )
)


# ------------------------------------------------------------
# 10. Direct-label helper columns
# ------------------------------------------------------------
#
# Initial label nudges are deliberately conservative.
# They can be adjusted in STEP 6X-9C without modifying
# any analytical values.
#
# ------------------------------------------------------------

figureX_core_6X$label_x <- c(
  -0.245,   # AOC3
  -0.285,   # MS4A7
  -0.300    # IL7R
)

figureX_core_6X$label_y <- c(
  0.57,     # AOC3
  0.78,     # MS4A7
  1.00      # IL7R
)


# ------------------------------------------------------------
# 11. Plotting annotation
# ------------------------------------------------------------

figureX_plot_6X$Point_size <- NA_real_

figureX_plot_6X$Point_size[
  figureX_plot_6X$Plot_class_final ==
    "Non-significant"
] <- 1.15

figureX_plot_6X$Point_size[
  figureX_plot_6X$Plot_class_final ==
    "Significant_favorable_direction"
] <- 1.55

figureX_plot_6X$Point_size[
  figureX_plot_6X$Plot_class_final ==
    "Significant_unfavorable_direction"
] <- 1.55

figureX_plot_6X$Point_size[
  figureX_plot_6X$Plot_class_final ==
    "Core_gene"
] <- 3.6


# ------------------------------------------------------------
# 12. Plotting alpha annotation
# ------------------------------------------------------------

figureX_plot_6X$Point_alpha <- NA_real_

figureX_plot_6X$Point_alpha[
  figureX_plot_6X$Plot_class_final ==
    "Non-significant"
] <- 0.42

figureX_plot_6X$Point_alpha[
  figureX_plot_6X$Plot_class_final ==
    "Significant_favorable_direction"
] <- 0.78

figureX_plot_6X$Point_alpha[
  figureX_plot_6X$Plot_class_final ==
    "Significant_unfavorable_direction"
] <- 0.78

figureX_plot_6X$Point_alpha[
  figureX_plot_6X$Plot_class_final ==
    "Core_gene"
] <- 1


# ------------------------------------------------------------
# 13. Verify original analytical values unchanged
# ------------------------------------------------------------

qc_rows_9B <- nrow(
  figureX_plot_6X
) == figureX_data_6X_nrow_before

qc_gene_order_9B <- identical(
  figureX_plot_6X$Gene,
  figureX_data_6X_gene_order_before
)

qc_HR_9B <- isTRUE(
  all.equal(
    figureX_plot_6X$HR,
    figureX_data_6X_HR_before
  )
)

qc_FDR_9B <- isTRUE(
  all.equal(
    figureX_plot_6X$FDR,
    figureX_data_6X_FDR_before
  )
)

qc_log2HR_9B <- isTRUE(
  all.equal(
    figureX_plot_6X$log2HR,
    figureX_data_6X_log2HR_before
  )
)


# ------------------------------------------------------------
# 14. Plotting-class QC
# ------------------------------------------------------------

plot_class_table_9B <- table(
  figureX_plot_6X$Plot_class_final,
  useNA = "ifany"
)

print(
  plot_class_table_9B
)

qc_core_n_9B <- sum(
  figureX_plot_6X$Core_gene_plot
) == 3

qc_sig_fav_n_9B <- sum(
  figureX_plot_6X$Plot_class_final ==
    "Significant_favorable_direction"
) == 1318

qc_sig_unfav_n_9B <- sum(
  figureX_plot_6X$Plot_class_final ==
    "Significant_unfavorable_direction"
) == 622

qc_non_sig_n_9B <- sum(
  figureX_plot_6X$Plot_class_final ==
    "Non-significant"
) == 29666


# ------------------------------------------------------------
# 15. Core-gene QC
# ------------------------------------------------------------

cat("\nCORE GENE PLOTTING DATA\n")
cat("------------------------------------------------------------\n")

print(
  figureX_core_6X,
  row.names = FALSE
)

qc_core_order_9B <- identical(
  as.character(figureX_core_6X$Gene_symbol),
  c("AOC3", "MS4A7", "IL7R")
)

qc_core_FDR_9B <- isTRUE(
  all(
    figureX_core_6X$FDR > 0.05
  )
)

qc_core_HR_9B <- isTRUE(
  all(
    figureX_core_6X$HR < 1
  )
)


# ------------------------------------------------------------
# 16. Reference-line QC
# ------------------------------------------------------------

qc_xref_9B <- identical(
  figureX_reference_6X$xintercept[1],
  0
)

qc_yref_9B <- isTRUE(
  all.equal(
    figureX_reference_6X$yintercept[2],
    -log10(0.05)
  )
)

qc_limits_9B <- identical(
  figureX_limits_6X$x,
  c(-0.75, 0.75)
) &&
  identical(
    figureX_limits_6X$y,
    c(0, 4.60)
  )


# ------------------------------------------------------------
# 17. Print complete QC summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 6X-9B QC SUMMARY\n")
cat("============================================================\n")

cat(
  "Rows unchanged:",
  qc_rows_9B,
  "\n"
)

cat(
  "Gene order unchanged:",
  qc_gene_order_9B,
  "\n"
)

cat(
  "HR unchanged:",
  qc_HR_9B,
  "\n"
)

cat(
  "FDR unchanged:",
  qc_FDR_9B,
  "\n"
)

cat(
  "log2HR unchanged:",
  qc_log2HR_9B,
  "\n"
)

cat(
  "Core genes = 3:",
  qc_core_n_9B,
  "\n"
)

cat(
  "Significant favorable = 1318:",
  qc_sig_fav_n_9B,
  "\n"
)

cat(
  "Significant unfavorable = 622:",
  qc_sig_unfav_n_9B,
  "\n"
)

cat(
  "Non-significant = 29666:",
  qc_non_sig_n_9B,
  "\n"
)

cat(
  "Core-gene order AOC3/MS4A7/IL7R:",
  qc_core_order_9B,
  "\n"
)

cat(
  "All core genes FDR > 0.05:",
  qc_core_FDR_9B,
  "\n"
)

cat(
  "All core genes HR < 1:",
  qc_core_HR_9B,
  "\n"
)

cat(
  "Vertical reference x = 0:",
  qc_xref_9B,
  "\n"
)

cat(
  "Horizontal reference FDR = 0.05:",
  qc_yref_9B,
  "\n"
)

cat(
  "Fixed plotting limits:",
  qc_limits_9B,
  "\n"
)


# ------------------------------------------------------------
# 18. Final pass/fail
# ------------------------------------------------------------

step6X9B_pass <- all(
  qc_rows_9B,
  qc_gene_order_9B,
  qc_HR_9B,
  qc_FDR_9B,
  qc_log2HR_9B,
  qc_core_n_9B,
  qc_sig_fav_n_9B,
  qc_sig_unfav_n_9B,
  qc_non_sig_n_9B,
  qc_core_order_9B,
  qc_core_FDR_9B,
  qc_core_HR_9B,
  qc_xref_9B,
  qc_yref_9B,
  qc_limits_9B
)

cat("\n============================================================\n")

if (step6X9B_pass) {
  cat("STEP 6X-9B: PASSED\n")
} else {
  cat("STEP 6X-9B: FAILED\n")
}

cat("============================================================\n")


# ------------------------------------------------------------
# 19. Save plotting-preparation checkpoint
# ------------------------------------------------------------

step6X9B_checkpoint_dir <-
  "result_dir/Paper 1_Step 6-survival analysis/STEP_6X_genome_wide_survival/checkpoint"

if (!dir.exists(step6X9B_checkpoint_dir)) {
  dir.create(
    step6X9B_checkpoint_dir,
    recursive = TRUE
  )
}

save(
  figureX_plot_6X,
  figureX_core_6X,
  figureX_reference_6X,
  figureX_limits_6X,
  fdr_threshold_6X,
  step6X9B_pass,
  file = file.path(
    step6X9B_checkpoint_dir,
    "Paper1_STEP6X9B_FigureX_plotting_preparation_checkpoint.RData"
  )
)

cat(
  "\nCheckpoint saved:\n",
  file.path(
    step6X9B_checkpoint_dir,
    "Paper1_STEP6X9B_FigureX_plotting_preparation_checkpoint.RData"
  ),
  "\n"
)

cat("\n============================================================\n")
cat("STEP 6X-9B COMPLETE\n")
cat("============================================================\n")


# ============================================================
# STEP 6X-9C-R1
# Inspect the exact contents of the locked STEP 6X-9B checkpoint
# ============================================================

checkpoint_9B <- file.path(
  "result_dir/Paper 1_Step 6-survival analysis/STEP_6X_genome_wide_survival",
  "checkpoint",
  "Paper1_STEP6X9B_FigureX_plotting_preparation_checkpoint.RData"
)

cat("\n============================================================\n")
cat("STEP 6X-9C-R1: Inspect STEP 6X-9B checkpoint\n")
cat("============================================================\n\n")

# Load into an isolated environment
env_9B <- new.env(parent = emptyenv())

loaded_names_9B <- load(
  checkpoint_9B,
  envir = env_9B
)

cat("Objects stored in STEP 6X-9B checkpoint:\n\n")
print(loaded_names_9B)

cat("\n------------------------------------------------------------\n")
cat("Detailed object structure\n")
cat("------------------------------------------------------------\n\n")

for (nm in loaded_names_9B) {
  
  obj <- get(
    nm,
    envir = env_9B
  )
  
  cat("\nObject:", nm, "\n")
  cat("Class :", paste(class(obj), collapse = ", "), "\n")
  
  if (is.data.frame(obj)) {
    cat(
      "Dim   :",
      paste(dim(obj), collapse = " x "),
      "\n"
    )
    cat(
      "Cols  :",
      paste(
        names(obj),
        collapse = ", "
      ),
      "\n"
    )
  } else if (is.list(obj)) {
    cat(
      "Length:",
      length(obj),
      "\n"
    )
    cat(
      "Names :",
      paste(
        names(obj),
        collapse = ", "
      ),
      "\n"
    )
  } else {
    cat(
      "Length:",
      length(obj),
      "\n"
    )
  }
}

cat("\n============================================================\n")
cat("STEP 6X-9C-R1 COMPLETE\n")
cat("============================================================\n")

# ============================================================
# STEP 6X-9C-R2
# Restore locked STEP 6X-9B objects into the current workspace
# ============================================================

cat("\n============================================================\n")
cat("STEP 6X-9C-R2: Restore STEP 6X-9B objects\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Confirm isolated checkpoint environment
# ------------------------------------------------------------

required_objects_9B <- c(
  "figureX_plot_6X",
  "figureX_core_6X",
  "figureX_reference_6X",
  "figureX_limits_6X",
  "fdr_threshold_6X",
  "step6X9B_pass"
)

missing_in_env_9B <- required_objects_9B[
  !vapply(
    required_objects_9B,
    exists,
    logical(1),
    envir = env_9B,
    inherits = FALSE
  )
]

if (length(missing_in_env_9B) > 0) {
  
  stop(
    "STEP 6X-9B checkpoint is missing object(s): ",
    paste(
      missing_in_env_9B,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 2. Restore objects without recalculation
# ------------------------------------------------------------

for (nm in required_objects_9B) {
  
  assign(
    nm,
    get(
      nm,
      envir = env_9B
    ),
    envir = .GlobalEnv
  )
}


# ------------------------------------------------------------
# 3. Verify restoration
# ------------------------------------------------------------

missing_after_restore <- required_objects_9B[
  !vapply(
    required_objects_9B,
    exists,
    logical(1),
    envir = .GlobalEnv,
    inherits = FALSE
  )
]

if (length(missing_after_restore) > 0) {
  
  stop(
    "Restoration failed for object(s): ",
    paste(
      missing_after_restore,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 4. Structural QC
# ------------------------------------------------------------

cat("Restored objects:\n\n")

for (nm in required_objects_9B) {
  
  obj <- get(
    nm,
    envir = .GlobalEnv
  )
  
  cat(
    nm,
    " | class = ",
    paste(class(obj), collapse = ", "),
    " | length = ",
    length(obj),
    "\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 5. Critical dimensional checks
# ------------------------------------------------------------

qc_9C_R2 <- list(
  
  figureX_plot_rows =
    nrow(figureX_plot_6X) == 31609,
  
  figureX_core_rows =
    nrow(figureX_core_6X) == 3,
  
  figureX_reference_rows =
    nrow(figureX_reference_6X) == 2,
  
  figureX_limits_is_list =
    is.list(figureX_limits_6X),
  
  fdr_threshold_is_numeric =
    is.numeric(fdr_threshold_6X),
  
  step6X9B_pass =
    isTRUE(step6X9B_pass)
  
)

cat("\n------------------------------------------------------------\n")
cat("STEP 6X-9C-R2 QC\n")
cat("------------------------------------------------------------\n\n")

print(qc_9C_R2)


# ------------------------------------------------------------
# 6. Verify exact core genes
# ------------------------------------------------------------

cat("\nCore genes restored from locked 9B object:\n")

print(
  figureX_core_6X[
    ,
    c(
      "Gene",
      "Gene_symbol",
      "HR",
      "FDR",
      "log2HR",
      "neglog10FDR"
    )
  ]
)


# ------------------------------------------------------------
# 7. Verify fixed plotting limits
# ------------------------------------------------------------

cat("\nLocked plotting limits:\n")

print(figureX_limits_6X)

cat("\nLocked FDR threshold:\n")

print(fdr_threshold_6X)


# ------------------------------------------------------------
# 8. Overall R2 result
# ------------------------------------------------------------

step6X9C_R2_pass <- all(
  unlist(qc_9C_R2)
)

cat("\n============================================================\n")

if (isTRUE(step6X9C_R2_pass)) {
  
  cat(
    "STEP 6X-9C-R2: PASSED\n"
  )
  
  cat(
    "All locked STEP 6X-9B objects restored successfully.\n"
  )
  
  cat(
    "No statistical recalculation was performed.\n"
  )
  
} else {
  
  cat(
    "STEP 6X-9C-R2: FAILED\n"
  )
  
}

cat("============================================================\n")


# ============================================================
# Paper 1
# STEP 6X-9C-R4
# Figure X final visual refinement:
# Triangular core-gene annotation + FDR threshold label
#
# IMPORTANT:
#   - No statistical recalculation
#   - No modification of gene-level statistical results
#   - No modification of x/y coordinates of gene points
#   - Only visual annotation helper coordinates are updated
# ============================================================


cat("\n============================================================\n")
cat("STEP 6X-9C-R4: Figure X visual refinement\n")
cat("Triangular core-gene labels + FDR threshold annotation\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Required STEP 6X-9B objects
# ------------------------------------------------------------

required_objects_9C_R4 <- c(
  "figureX_plot_6X",
  "figureX_core_6X",
  "figureX_reference_6X",
  "figureX_limits_6X",
  "fdr_threshold_6X",
  "step6X9B_pass"
)

missing_objects_9C_R4 <- required_objects_9C_R4[
  !vapply(
    required_objects_9C_R4,
    exists,
    logical(1),
    envir = .GlobalEnv,
    inherits = FALSE
  )
]

if (length(missing_objects_9C_R4) > 0) {
  stop(
    paste0(
      "Missing required STEP 6X-9B object(s): ",
      paste(missing_objects_9C_R4, collapse = ", ")
    )
  )
}

if (!isTRUE(step6X9B_pass)) {
  stop("STEP 6X-9B QC is not TRUE. Do not proceed.")
}


# ------------------------------------------------------------
# 1. Load plotting packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(grid)
})


# ------------------------------------------------------------
# 2. Create plotting copies
#
# IMPORTANT:
# Original locked objects are NOT modified.
# ------------------------------------------------------------

figureX_plot_R4 <- figureX_plot_6X
figureX_core_R4 <- figureX_core_6X


# ------------------------------------------------------------
# 3. Confirm statistical coordinates are unchanged
# ------------------------------------------------------------

required_plot_columns_R4 <- c(
  "Gene",
  "Gene_symbol",
  "HR",
  "FDR",
  "log2HR",
  "neglog10FDR",
  "x",
  "y",
  "Plot_class_final"
)

missing_plot_columns_R4 <- setdiff(
  required_plot_columns_R4,
  colnames(figureX_plot_R4)
)

if (length(missing_plot_columns_R4) > 0) {
  stop(
    paste0(
      "Missing required plotting column(s): ",
      paste(missing_plot_columns_R4, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 4. Verify the 31,609-gene plotting dataset
# ------------------------------------------------------------

if (nrow(figureX_plot_R4) != 31609) {
  stop("Unexpected number of genes in Figure X plotting dataset.")
}

if (anyDuplicated(figureX_plot_R4$Gene) > 0) {
  stop("Duplicate genes detected in Figure X plotting dataset.")
}

if (any(!is.finite(figureX_plot_R4$x))) {
  stop("Non-finite x values detected.")
}

if (any(!is.finite(figureX_plot_R4$y))) {
  stop("Non-finite y values detected.")
}


# ------------------------------------------------------------
# 5. Verify core genes
# ------------------------------------------------------------

core_symbols_R4 <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_ensembl_R4 <- c(
  "ENSG00000131471",
  "ENSG00000166927",
  "ENSG00000168685"
)

core_check_R4 <- figureX_plot_R4[
  figureX_plot_R4$Gene %in% core_ensembl_R4,
  ,
  drop = FALSE
]

if (nrow(core_check_R4) != 3) {
  stop("Exactly three core genes are required.")
}

if (!setequal(core_check_R4$Gene_symbol, core_symbols_R4)) {
  stop("Core-gene identity check failed.")
}


# ------------------------------------------------------------
# 6. Freeze actual gene-point coordinates
#
# These coordinates MUST NOT be changed.
# ------------------------------------------------------------

core_point_coordinates_R4 <- core_check_R4 %>%
  dplyr::select(
    Gene,
    Gene_symbol,
    x,
    y,
    HR,
    FDR,
    log2HR,
    neglog10FDR
  ) %>%
  arrange(match(Gene_symbol, core_symbols_R4))


cat("\nLocked core-gene point coordinates:\n")
print(core_point_coordinates_R4)


# ------------------------------------------------------------
# 7. Triangular label layout
#
# IMPORTANT:
# Only label_x / label_y are changed.
# The actual core-gene point coordinates remain untouched.
#
# Layout:
#
#             MS4A7
#                /
#               /
#      AOC3    /
#          \  /
#           \/       IL7R
#
# More precisely:
#   AOC3  -> left/down
#   MS4A7 -> left/up
#   IL7R  -> right/up
# ------------------------------------------------------------
triangle_labels_R4 <- data.frame(
  Gene_symbol = c("AOC3", "MS4A7", "IL7R"),
  
  label_x = c(
    -0.40,   # AOC3：left/down
    -0.50,   # MS4A7：left/up
    0.02    # IL7R：right/up
  ),
  
  label_y = c(
    0.42,    # AOC3
    0.98,    # MS4A7
    1.13     # IL7R
  ),
  
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 8. Verify label order
# ------------------------------------------------------------

if (!identical(
  triangle_labels_R4$Gene_symbol,
  core_symbols_R4
)) {
  stop("Triangle label order is incorrect.")
}


# ------------------------------------------------------------
# 9. Merge visual-only label coordinates
# ------------------------------------------------------------

figureX_core_R4 <- figureX_core_R4 %>%
  dplyr::select(
    -dplyr::any_of(c("label_x", "label_y"))
  ) %>%
  left_join(
    triangle_labels_R4,
    by = "Gene_symbol"
  )


# ------------------------------------------------------------
# 10. Check that all labels received coordinates
# ------------------------------------------------------------

if (any(is.na(figureX_core_R4$label_x)) ||
    any(is.na(figureX_core_R4$label_y))) {
  stop("One or more core genes are missing label coordinates.")
}


# ------------------------------------------------------------
# 11. FDR threshold
# ------------------------------------------------------------

fdr_threshold_R4 <- -log10(0.05)


if (!isTRUE(all.equal(
  as.numeric(fdr_threshold_R4),
  as.numeric(fdr_threshold_6X)
))) {
  stop("FDR threshold does not match locked STEP 6X-9B value.")
}


# ------------------------------------------------------------
# 12. Reference lines
# ------------------------------------------------------------

figureX_reference_R4 <- data.frame(
  reference = c(
    "HR = 1",
    "FDR = 0.05"
  ),
  xintercept = c(
    0,
    NA_real_
  ),
  yintercept = c(
    NA_real_,
    fdr_threshold_R4
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 13. Output directory
# ------------------------------------------------------------

figureX_output_dir_R4 <-
  "figure_dir/Paper 1_Step 6-survival analysis/STEP_6X_genome_wide_survival/figure/FigureX_FINAL_EXPORT"

if (!dir.exists(figureX_output_dir_R4)) {
  dir.create(
    figureX_output_dir_R4,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ------------------------------------------------------------
# 14. Colors
# ------------------------------------------------------------

col_nonsig_R4 <- "#BDBDBD"
col_sig_fav_R4 <- "#3DB7D3"
col_sig_unfav_R4 <- "#E64B35"
col_core_R4 <- "#7A3E9D"


# ------------------------------------------------------------
# 15. Plot class checks
# ------------------------------------------------------------

if (!"Plot_class_final" %in% colnames(figureX_plot_R4)) {
  stop("Plot_class_final is missing.")
}

figureX_plot_R4$Plot_class_final <- factor(
  figureX_plot_R4$Plot_class_final,
  levels = c(
    "Non-significant",
    "Significant_favorable_direction",
    "Significant_unfavorable_direction",
    "Core_gene"
  )
)

# ------------------------------------------------------------
# Core-gene leader-line geometry
#
# Start slightly outside the hollow circle so that the
# leader line does not pass through the marker.
# ------------------------------------------------------------

figureX_core_R4 <- figureX_core_R4 %>%
  mutate(
    dx = label_x - x,
    dy = label_y - y,
    distance = sqrt(dx^2 + dy^2),
    
    start_gap = 0.030,
    
    leader_x =
      x + dx / distance * start_gap,
    
    leader_y =
      y + dy / distance * start_gap,
    
    end_gap = case_when(
      Gene_symbol == "AOC3"  ~ 0.040,
      Gene_symbol == "MS4A7" ~ 0.050,
      Gene_symbol == "IL7R"  ~ 0.050,
      TRUE ~ 0.045
    ),
    
    leader_xend =
      label_x -
      dx / distance * end_gap,
    
    leader_yend =
      label_y -
      dy / distance * end_gap
  )

# ------------------------------------------------------------
# 16. Build Figure X
# ------------------------------------------------------------

p_FigureX_R4 <- ggplot(
  figureX_plot_R4,
  aes(
    x = x,
    y = y
  )
) +
  
  # ----------------------------------------------------------
# Background / non-significant genes
# ----------------------------------------------------------

geom_point(
  data = subset(
    figureX_plot_R4,
    Plot_class_final == "Non-significant"
  ),
  aes(
    size = Point_size,
    alpha = Point_alpha
  ),
  color = col_nonsig_R4,
  shape = 16
) +
  
  # ----------------------------------------------------------
# Significant genes: HR < 1
# ----------------------------------------------------------

geom_point(
  data = subset(
    figureX_plot_R4,
    Plot_class_final ==
      "Significant_favorable_direction"
  ),
  aes(
    size = Point_size,
    alpha = Point_alpha
  ),
  color = col_sig_fav_R4,
  shape = 16
) +
  
  # ----------------------------------------------------------
# Significant genes: HR > 1
# ----------------------------------------------------------

geom_point(
  data = subset(
    figureX_plot_R4,
    Plot_class_final ==
      "Significant_unfavorable_direction"
  ),
  aes(
    size = Point_size,
    alpha = Point_alpha
  ),
  color = col_sig_unfav_R4,
  shape = 16
) +
  
  # ----------------------------------------------------------
# Reference line: HR = 1
#
# x = 0 because x = log2(HR)
# ----------------------------------------------------------

geom_vline(
  xintercept = 0,
  linetype = "dashed",
  linewidth = 0.45,
  color = "grey40"
) +
  
  # ----------------------------------------------------------
# Reference line: FDR = 0.05
# ----------------------------------------------------------

geom_hline(
  yintercept = fdr_threshold_R4,
  linetype = "dashed",
  linewidth = 0.45,
  color = "grey40"
) +
  
  # ----------------------------------------------------------
# FDR threshold annotation
#
# Positioned just above the horizontal dashed line.
# ----------------------------------------------------------

annotate(
  "text",
  x = 0.70,
  y = fdr_threshold_R4 + 0.085,
  label = "FDR = 0.05",
  hjust = 1,
  vjust = 0,
  family = "Arial",
  fontface = "plain",
  size = 3.25,
  color = "grey25"
) +
  
  # ----------------------------------------------------------
# Core genes:
# purple hollow circles
# ----------------------------------------------------------

geom_point(
  data = figureX_core_R4,
  aes(
    x = x,
    y = y
  ),
  shape = 21,
  size = 4.1,
  stroke = 1.15,
  fill = "white",
  color = col_core_R4
) +
  
  # ----------------------------------------------------------
# Core-gene leader lines
#
# segment starts at actual gene point and ends at
# the triangular label position.
# ----------------------------------------------------------
geom_segment(
  data = figureX_core_R4,
  aes(
    x = leader_x,
    y = leader_y,
    xend = leader_xend,
    yend = leader_yend
  ),
  linewidth = 0.45,
  color = col_core_R4,
  lineend = "round"
) +
  
  # ----------------------------------------------------------
# Core-gene labels
# ----------------------------------------------------------

geom_text(
  data = figureX_core_R4,
  aes(
    x = label_x,
    y = label_y,
    label = Gene_symbol
  ),
  family = "Arial",
  fontface = "bold",
  size = 3.55,
  color = col_core_R4,
  hjust = 0.5,
  vjust = 0.5
) +
  
  # ----------------------------------------------------------
# Scales
# ----------------------------------------------------------

scale_size_identity() +
  
  scale_alpha_identity() +
  
  # ----------------------------------------------------------
# Coordinate limits
# ----------------------------------------------------------

coord_cartesian(
  xlim = figureX_limits_6X$x,
  ylim = figureX_limits_6X$y,
  expand = FALSE,
  clip = "off"
) +
  
  # ----------------------------------------------------------
# Labels
# ----------------------------------------------------------

labs(
  title =
    "Genome-wide survival associations across 31,609 genes in LUAD",
  
  x = expression(log[2](HR)),
  
  y = expression(-log[10](FDR))
) +
  
  # ----------------------------------------------------------
# Theme
# ----------------------------------------------------------

theme_classic(
  base_family = "Arial"
) +
  
  theme(
    
    # title
    plot.title = element_text(
      family = "Arial",
      face = "bold",
      size = 12,
      hjust = 0.5,
      margin = margin(
        b = 8
      )
    ),
    
    # axis titles
    axis.title.x = element_text(
      family = "Arial",
      size = 10.5,
      margin = margin(
        t = 7
      )
    ),
    
    axis.title.y = element_text(
      family = "Arial",
      size = 10.5,
      margin = margin(
        r = 7
      )
    ),
    
    # axis text
    axis.text = element_text(
      family = "Arial",
      size = 9
    ),
    
    # panel border
    panel.border = element_rect(
      fill = NA,
      linewidth = 0.55,
      color = "black"
    ),
    
    # ticks
    axis.ticks = element_line(
      linewidth = 0.45,
      color = "black"
    ),
    
    axis.ticks.length = unit(
      2.2,
      "pt"
    ),
    
    # plot margins
    plot.margin = margin(
      t = 10,
      r = 15,
      b = 10,
      l = 12
    ),
    
    # background
    panel.background = element_blank(),
    plot.background = element_blank()
  )


# ------------------------------------------------------------
# 17. Display plot
# ------------------------------------------------------------

print(p_FigureX_R4)


# ------------------------------------------------------------
# 18. Save visual-only annotation object
#
# This object contains only updated label positions.
# Statistical values remain inherited from locked STEP 6X-9B.
# ------------------------------------------------------------

figureX_visual_annotation_R4 <- list(
  core_labels = figureX_core_R4 %>%
    dplyr::select(
      Gene,
      Gene_symbol,
      x,
      y,
      label_x,
      label_y
    ),
  
  fdr_threshold = fdr_threshold_R4,
  
  fdr_threshold_label =
    "FDR = 0.05",
  
  hr_reference =
    "HR = 1",
  
  x_limits =
    figureX_limits_6X$x,
  
  y_limits =
    figureX_limits_6X$y
)


# ------------------------------------------------------------
# 19. Export paths
# ------------------------------------------------------------

pdf_file_R4 <- file.path(
  figureX_output_dir_R4,
  "FigureX_genome_wide_survival_associations_FINAL_R4.pdf"
)

tiff_file_R4 <- file.path(
  figureX_output_dir_R4,
  "FigureX_genome_wide_survival_associations_FINAL_R4_600dpi.tiff"
)

png_file_R4 <- file.path(
  figureX_output_dir_R4,
  "FigureX_genome_wide_survival_associations_FINAL_R4_600dpi.png"
)

jpeg_file_R4 <- file.path(
  figureX_output_dir_R4,
  "FigureX_genome_wide_survival_associations_FINAL_R4_600dpi.jpg"
)


# ------------------------------------------------------------
# 20. PDF export
# ------------------------------------------------------------

ggsave(
  filename = pdf_file_R4,
  plot = p_FigureX_R4,
  device = cairo_pdf,
  width = 8.8,
  height = 5.8,
  units = "in",
  bg = "white"
)


# ------------------------------------------------------------
# 21. TIFF export
# ------------------------------------------------------------

ggsave(
  filename = tiff_file_R4,
  plot = p_FigureX_R4,
  device = "tiff",
  width = 8.8,
  height = 5.8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 22. PNG export
# ------------------------------------------------------------

ggsave(
  filename = png_file_R4,
  plot = p_FigureX_R4,
  device = "png",
  width = 8.8,
  height = 5.8,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ------------------------------------------------------------
# 23. JPEG export
# ------------------------------------------------------------

ggsave(
  filename = jpeg_file_R4,
  plot = p_FigureX_R4,
  device = "jpeg",
  width = 8.8,
  height = 5.8,
  units = "in",
  dpi = 600,
  quality = 100,
  bg = "white"
)


# ------------------------------------------------------------
# 24. Final visual QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 6X-9C-R4 FINAL VISUAL QC\n")
cat("============================================================\n")

cat(
  "Gene rows:",
  nrow(figureX_plot_R4),
  " (expected 31609)\n"
)

cat(
  "Core genes:",
  nrow(figureX_core_R4),
  " (expected 3)\n"
)

cat(
  "Core order:",
  paste(
    figureX_core_R4$Gene_symbol,
    collapse = " / "
  ),
  "\n"
)

cat(
  "FDR threshold:",
  format(
    fdr_threshold_R4,
    digits = 8
  ),
  "\n"
)

cat(
  "X limits:",
  paste(
    figureX_limits_6X$x,
    collapse = " to "
  ),
  "\n"
)

cat(
  "Y limits:",
  paste(
    figureX_limits_6X$y,
    collapse = " to "
  ),
  "\n"
)

cat(
  "Core point coordinates unchanged:",
  identical(
    core_point_coordinates_R4 %>%
      dplyr::select(
        Gene,
        x,
        y
      ),
    core_check_R4 %>%
      dplyr::select(
        Gene,
        x,
        y
      )
  ),
  "\n"
)

cat(
  "FDR threshold label:",
  "FDR = 0.05",
  "\n"
)

cat(
  "PDF exists:",
  file.exists(pdf_file_R4),
  "\n"
)

cat(
  "TIFF exists:",
  file.exists(tiff_file_R4),
  "\n"
)

cat(
  "PNG exists:",
  file.exists(png_file_R4),
  "\n"
)

cat(
  "JPEG exists:",
  file.exists(jpeg_file_R4),
  "\n"
)


# ------------------------------------------------------------
# 25. Save R4 visual checkpoint
# ------------------------------------------------------------

checkpoint_R4_dir <-
  "result_dir/Paper 1_Step 6-survival analysis/STEP_6X_genome_wide_survival/checkpoint"

if (!dir.exists(checkpoint_R4_dir)) {
  dir.create(
    checkpoint_R4_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

checkpoint_R4_file <- file.path(
  checkpoint_R4_dir,
  "Paper1_STEP6X9C_R4_FigureX_visual_refinement_checkpoint.RData"
)

save(
  figureX_core_R4,
  figureX_reference_R4,
  figureX_visual_annotation_R4,
  fdr_threshold_R4,
  p_FigureX_R4,
  pdf_file_R4,
  tiff_file_R4,
  png_file_R4,
  jpeg_file_R4,
  file = checkpoint_R4_file
)


# ------------------------------------------------------------
# 26. Final message
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 6X-9C-R4: COMPLETED\n")
cat("============================================================\n")

cat(
  "Triangle annotation: AOC3 left/down; ",
  "MS4A7 left/up; IL7R right/up.\n"
)

cat(
  "FDR threshold annotation added: FDR = 0.05.\n"
)

cat(
  "HR = 1 reference line retained without text annotation.\n"
)

cat(
  "No statistical recalculation performed.\n"
)

cat(
  "No gene-level statistical values modified.\n"
)

cat(
  "No actual gene-point coordinates modified.\n"
)

cat(
  "All exports written to:\n",
  figureX_output_dir_R4,
  "\n"
)

cat(
  "Checkpoint:\n",
  checkpoint_R4_file,
  "\n"
)

cat("\nSTEP 6X-9C-R4: PASSED\n")
# ------------------------------------------------------------
# 17. Display plot
# ------------------------------------------------------------

print(p_FigureX_R4)


# ------------------------------------------------------------
# 18. Save visual-only annotation object
#
# This object contains only updated label positions.
# Statistical values remain inherited from locked STEP 6X-9B.
# ------------------------------------------------------------

figureX_visual_annotation_R4 <- list(
  core_labels = figureX_core_R4 %>%
    dplyr::select(
      Gene,
      Gene_symbol,
      x,
      y,
      label_x,
      label_y
    ),
  
  fdr_threshold = fdr_threshold_R4,
  
  fdr_threshold_label =
    "FDR = 0.05",
  
  hr_reference =
    "HR = 1",
  
  x_limits =
    figureX_limits_6X$x,
  
  y_limits =
    figureX_limits_6X$y
)


# ------------------------------------------------------------
# 19. Export paths
# ------------------------------------------------------------

pdf_file_R4 <- file.path(
  figureX_output_dir_R4,
  "FigureX_genome_wide_survival_associations_FINAL_R4.pdf"
)

tiff_file_R4 <- file.path(
  figureX_output_dir_R4,
  "FigureX_genome_wide_survival_associations_FINAL_R4_600dpi.tiff"
)

png_file_R4 <- file.path(
  figureX_output_dir_R4,
  "FigureX_genome_wide_survival_associations_FINAL_R4_600dpi.png"
)

jpeg_file_R4 <- file.path(
  figureX_output_dir_R4,
  "FigureX_genome_wide_survival_associations_FINAL_R4_600dpi.jpg"
)




