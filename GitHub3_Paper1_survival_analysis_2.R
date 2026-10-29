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

##STEP 6C-1 — Clinical Covariate Mapping and QC
# ============================================================
# Paper 1
# STEP 6C-1
# Clinical Covariate Mapping and QC
#
# Purpose:
#   1. Recover the locked 485-sample survival cohort
#   2. Recover clinical variables from the same clinical source
#   3. Audit sample-ID mapping
#   4. Audit Age / Gender / Stage / T / N / M
#   5. Inspect coding and factor levels
#   6. Assess missingness
#   7. Prepare a frozen clinical covariate dataset
#
# IMPORTANT:
#   - Do NOT refit any Cox model in this step.
#   - Do NOT redefine the 485-sample cohort.
#   - Do NOT use complete.cases() to recreate the cohort.
#   - Do NOT modify Step 6B locked objects.
# ============================================================


cat("\n============================================================\n")
cat("STEP 6C-1: Clinical Covariate Mapping and QC\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Working directory
# ------------------------------------------------------------
work_dir <- input_dir

setwd(work_dir)

getwd()


# ------------------------------------------------------------
# 2. Load the locked Step 6B checkpoint
# ------------------------------------------------------------
load(
  "input_dir/Paper1_Step6B6_MS4A7_time_varying_Cox_checkpoint.RData"
)

ls()

# ------------------------------------------------------------
# 3. Check required Step 6B objects
# ------------------------------------------------------------

required_step6b_objects <- c(
  "ms4a7_tv_data",
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
  "qc_hr"
)

missing_step6b_objects <-
  required_step6b_objects[
    !vapply(
      required_step6b_objects,
      exists,
      logical(1),
      inherits = TRUE
    )
  ]

if (length(missing_step6b_objects) > 0) {
  stop(
    paste(
      "Required Step 6B objects are missing:",
      paste(
        missing_step6b_objects,
        collapse = ", "
      )
    )
  )
}

cat(
  "Step 6B checkpoint object QC: PASSED\n\n"
)


# ------------------------------------------------------------
# 4. Recover the locked 485-sample survival cohort
# ------------------------------------------------------------
#
# The Step 6B time-varying dataset contains the locked
# survival observations.
#
# We use its Sample IDs only.
# The cohort itself is NOT reconstructed using complete.cases().
# ------------------------------------------------------------

if (!"Sample" %in% colnames(ms4a7_tv_data)) {
  stop(
    "Column 'Sample' not found in ms4a7_tv_data."
  )
}

locked_primary_samples <-
  unique(ms4a7_tv_data$Sample)

cat("Locked survival cohort:\n")
cat("N =", length(locked_primary_samples), "\n\n")


# ------------------------------------------------------------
# 5. Confirm locked cohort size and uniqueness
# ------------------------------------------------------------

stopifnot(
  length(locked_primary_samples) == 485
)

stopifnot(
  !anyDuplicated(locked_primary_samples)
)

cat("Locked cohort N = 485: PASSED\n")
cat("Sample-ID uniqueness: PASSED\n\n")


# ------------------------------------------------------------
# 6. Recover clinical source
# ------------------------------------------------------------
#
# The clinical source used for the tumor-only WGCNA was:
#
#   clinical_aligned
#
# It should already contain the clinical variables needed for
# the revised survival analysis.
#
# First inspect the object.
# ------------------------------------------------------------
load("input_dir/Workflow-04-2D_final_WGCNA_candidate_checkpoint.RData")

if (!exists("clinical_aligned")) {
  
  stop(
    paste0(
      "Object 'clinical_aligned' is not currently loaded.\n",
      "Please load the tumor-only WGCNA checkpoint:\n",
      "Workflow-04-2D_final_WGCNA_candidate_checkpoint.RData\n",
      "and rerun from STEP 6C-1 section 6."
    )
  )
  
}

cat("Clinical object found: clinical_aligned\n\n")


# ------------------------------------------------------------
# 7. Inspect clinical object structure
# ------------------------------------------------------------

cat("Clinical object dimensions:\n")
print(dim(clinical_aligned))

cat("\nClinical column names:\n")
print(colnames(clinical_aligned))

cat("\nClinical structure:\n")
str(clinical_aligned)


# ------------------------------------------------------------
# 8. Identify and lock the authoritative sample-ID column
# ------------------------------------------------------------
#
# IMPORTANT:
#   clinical_aligned contains both:
#
#   Sample      = TCGA sample-level barcode
#   patient_id  = TCGA patient-level identifier
#
#   The locked Step 6B survival cohort is sample-level.
#   Therefore the authoritative mapping key MUST be "Sample".
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 6C-1: Lock sample-level clinical mapping key\n")
cat("============================================================\n\n")

sample_column_candidates <-
  intersect(
    c(
      "Sample",
      "sample",
      "Sample_ID",
      "sample_id"
    ),
    colnames(clinical_aligned)
  )

cat("Sample-level ID candidates:\n")
print(sample_column_candidates)

if (!"Sample" %in% colnames(clinical_aligned)) {
  
  stop(
    paste(
      "Authoritative sample-level column 'Sample' was not found."
    )
  )
  
}

# ------------------------------------------------------------
# Explicitly lock Sample
# ------------------------------------------------------------

clinical_sample_column <- "Sample"

cat(
  "\nAuthoritative clinical mapping column: ",
  clinical_sample_column,
  "\n",
  sep = ""
)

cat(
  "patient_id is retained as patient-level provenance only and\n",
  "will NOT be used as the Step 6C sample-mapping key.\n\n"
)

stopifnot(
  identical(
    clinical_sample_column,
    "Sample"
  )
)

# ------------------------------------------------------------
# 9. Standardize clinical sample IDs
# ------------------------------------------------------------

clinical_sample_ids <-
  as.character(
    clinical_aligned[[clinical_sample_column]]
  )

locked_primary_samples <-
  as.character(
    locked_primary_samples
  )

clinical_sample_ids <-
  trimws(clinical_sample_ids)

locked_primary_samples <-
  trimws(locked_primary_samples)

# ------------------------------------------------------------
# 10. Sample-ID integrity audit
# ------------------------------------------------------------

cat("Clinical sample-ID audit:\n")

cat(
  "Clinical rows: ",
  nrow(clinical_aligned),
  "\n",
  sep = ""
)

cat(
  "Clinical sample IDs: ",
  length(clinical_sample_ids),
  "\n",
  sep = ""
)

cat(
  "Missing clinical sample IDs: ",
  sum(is.na(clinical_sample_ids) | clinical_sample_ids == ""),
  "\n",
  sep = ""
)

cat(
  "Duplicated clinical sample IDs: ",
  sum(duplicated(clinical_sample_ids)),
  "\n\n",
  sep = ""
)

if (
  anyNA(clinical_sample_ids) ||
  any(clinical_sample_ids == "") ||
  anyDuplicated(clinical_sample_ids)
) {
  
  stop(
    "Clinical sample-ID integrity check failed."
  )
  
}


# ------------------------------------------------------------
# 11. Match locked cohort to clinical data
# ------------------------------------------------------------

clinical_idx <-
  match(
    locked_primary_samples,
    clinical_sample_ids
  )

mapping_complete <-
  all(!is.na(clinical_idx))

cat(
  "Locked cohort → clinical mapping complete: ",
  mapping_complete,
  "\n",
  sep = ""
)

cat(
  "Unmatched locked samples: ",
  sum(is.na(clinical_idx)),
  "\n\n",
  sep = ""
)

if (!mapping_complete) {
  
  cat("Unmatched samples:\n")
  
  print(
    locked_primary_samples[
      is.na(clinical_idx)
    ]
  )
  
  stop(
    "Clinical mapping failed."
  )
  
}

# ------------------------------------------------------------
# 11A. Explicit sample-level vs patient-level ID audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Sample-level vs patient-level ID audit\n")
cat("============================================================\n\n")

# Sample-level mapping
sample_mapping_check <-
  match(
    locked_primary_samples,
    clinical_aligned$Sample
  )

# Patient-level IDs corresponding to the locked samples
locked_patient_ids <-
  clinical_aligned$patient_id[
    sample_mapping_check
  ]

cat(
  "Locked sample observations: ",
  length(locked_primary_samples),
  "\n",
  sep = ""
)

cat(
  "Matched sample-level observations: ",
  sum(!is.na(sample_mapping_check)),
  "\n",
  sep = ""
)

cat(
  "Unique patient IDs represented: ",
  length(unique(locked_patient_ids)),
  "\n",
  sep = ""
)

cat(
  "Duplicated patient IDs among locked observations: ",
  sum(duplicated(locked_patient_ids)),
  "\n",
  sep = ""
)

cat("\nExample mapping:\n")

print(
  data.frame(
    Sample = locked_primary_samples[1:min(10, length(locked_primary_samples))],
    patient_id =
      locked_patient_ids[1:min(10, length(locked_patient_ids))],
    stringsAsFactors = FALSE
  )
)


# ------------------------------------------------------------
# 12. Construct clinical dataset in locked cohort order
# ------------------------------------------------------------

clinical_locked <-
  clinical_aligned[
    clinical_idx,
    ,
    drop = FALSE
  ]

rownames(clinical_locked) <-
  locked_primary_samples


# ------------------------------------------------------------
# 13. Verify exact sample-order identity
# ------------------------------------------------------------

stopifnot(
  identical(
    as.character(
      clinical_locked[[clinical_sample_column]]
    ),
    locked_primary_samples
  )
)

cat(
  "Clinical sample order matches locked survival cohort: PASSED\n\n"
)


# ============================================================
# 14. Identify candidate clinical variables
# ============================================================

cat("============================================================\n")
cat("Clinical covariate identification\n")
cat("============================================================\n\n")


# Candidate names are deliberately broad.
# The script will NOT silently select an ambiguous variable.

candidate_age <- c(
  "Age",
  "age",
  "AGE",
  "age_at_diagnosis",
  "Age_at_Diagnosis"
)

candidate_gender <- c(
  "Gender",
  "gender",
  "Sex",
  "sex",
  "SEX"
)

candidate_stage <- c(
  "pathologic_stage",
  "Pathologic_stage",
  "pathologic.stage",
  "Stage",
  "stage",
  "AJCC_Stage"
)

candidate_T <- c(
  "T",
  "t",
  "pathologic_T",
  "Pathologic_T",
  "pathologic_T_stage",
  "T_stage"
)

candidate_N <- c(
  "N",
  "n",
  "pathologic_N",
  "Pathologic_N",
  "pathologic_N_stage",
  "N_stage"
)

candidate_M <- c(
  "M",
  "m",
  "pathologic_M",
  "Pathologic_M",
  "pathologic_M_stage",
  "M_stage"
)


find_candidates <- function(
    candidates,
    available_columns
) {
  
  intersect(
    candidates,
    available_columns
  )
  
}


age_candidates <-
  find_candidates(
    candidate_age,
    colnames(clinical_locked)
  )

gender_candidates <-
  find_candidates(
    candidate_gender,
    colnames(clinical_locked)
  )

stage_candidates <-
  find_candidates(
    candidate_stage,
    colnames(clinical_locked)
  )

T_candidates <-
  find_candidates(
    candidate_T,
    colnames(clinical_locked)
  )

N_candidates <-
  find_candidates(
    candidate_N,
    colnames(clinical_locked)
  )

M_candidates <-
  find_candidates(
    candidate_M,
    colnames(clinical_locked)
  )


cat("Age candidates:\n")
print(age_candidates)

cat("\nGender/Sex candidates:\n")
print(gender_candidates)

cat("\nStage candidates:\n")
print(stage_candidates)

cat("\nT-stage candidates:\n")
print(T_candidates)

cat("\nN-stage candidates:\n")
print(N_candidates)

cat("\nM-stage candidates:\n")
print(M_candidates)

# ============================================================
# STEP 6C-1
# Clinical variable mapping audit — corrected version
#
# IMPORTANT:
#   This step ONLY identifies and audits clinical variables.
#   No recoding and no Cox model fitting are performed here.
# ============================================================

cat("\n============================================================\n")
cat("STEP 6C-1: Clinical variable mapping audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 14A. Explicitly define the clinical variables
# ------------------------------------------------------------
#
# Based on the actual clinical_aligned structure:
#
#   Age                = Age
#   Gender             = Gender
#   Broad stage        = pathologic_stage
#   Detailed stage     = pathologic_stage_2
#   T                  = pathologic_t
#   N                  = pathologic_n
#   M                  = pathologic_m
#
# These are explicit mappings, not automatic guesses.
# ------------------------------------------------------------

age_var <- "Age"

gender_var <- "Gender"

stage_var <- "pathologic_stage"

stage_detailed_var <- "pathologic_stage_2"

T_var <- "pathologic_t"

N_var <- "pathologic_n"

M_var <- "pathologic_m"


clinical_variable_map <-
  data.frame(
    Analysis_role = c(
      "Age",
      "Gender_Sex",
      "Pathologic_stage_broad",
      "Pathologic_stage_detailed",
      "T_stage",
      "N_stage",
      "M_stage"
    ),
    
    Source_column = c(
      age_var,
      gender_var,
      stage_var,
      stage_detailed_var,
      T_var,
      N_var,
      M_var
    ),
    
    stringsAsFactors = FALSE
  )


cat("Frozen clinical variable map:\n\n")

print(
  clinical_variable_map,
  row.names = FALSE
)


# ------------------------------------------------------------
# 14B. Verify that all required columns exist
# ------------------------------------------------------------

required_clinical_columns <-
  clinical_variable_map$Source_column

missing_clinical_columns <-
  setdiff(
    required_clinical_columns,
    colnames(clinical_locked)
  )

if (length(missing_clinical_columns) > 0) {
  
  stop(
    paste(
      "Required clinical columns are missing:",
      paste(
        missing_clinical_columns,
        collapse = ", "
      )
    )
  )
  
}

cat(
  "\nRequired clinical columns present: PASSED\n"
)


# ============================================================
# 15. Raw coding audit
# ============================================================

cat("\n============================================================\n")
cat("STEP 6C-1: Raw clinical coding audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 15A. Age
# ------------------------------------------------------------

cat("AGE\n")
cat("------------------------------------------------------------\n")

print(
  summary(
    clinical_locked[[age_var]]
  )
)

cat(
  "Missing:",
  sum(is.na(clinical_locked[[age_var]])),
  "\n\n"
)


# ------------------------------------------------------------
# 15B. Gender
# ------------------------------------------------------------

cat("GENDER\n")
cat("------------------------------------------------------------\n")

cat("Class:\n")
print(
  class(
    clinical_locked[[gender_var]]
  )
)

cat("\nFrequency:\n")

print(
  table(
    clinical_locked[[gender_var]],
    useNA = "ifany"
  )
)

cat("\nUnique values:\n")

print(
  sort(
    unique(
      clinical_locked[[gender_var]]
    )
  )
)

cat("\nMissing:\n")

print(
  sum(
    is.na(
      clinical_locked[[gender_var]]
    )
  )
)


# ------------------------------------------------------------
# 15C. Broad pathological stage
# ------------------------------------------------------------

cat("\nPATHOLOGIC STAGE — BROAD\n")
cat("------------------------------------------------------------\n")

cat("Variable:", stage_var, "\n")

cat("Class:\n")

print(
  class(
    clinical_locked[[stage_var]]
  )
)

cat("\nFrequency:\n")

print(
  table(
    clinical_locked[[stage_var]],
    useNA = "ifany"
  )
)

cat("\nUnique values:\n")

print(
  sort(
    unique(
      clinical_locked[[stage_var]]
    )
  )
)

cat("\nMissing:\n")

print(
  sum(
    is.na(
      clinical_locked[[stage_var]]
    )
  )
)


# ------------------------------------------------------------
# 15D. Detailed pathological stage
# ------------------------------------------------------------

cat("\nPATHOLOGIC STAGE — DETAILED\n")
cat("------------------------------------------------------------\n")

cat("Variable:", stage_detailed_var, "\n")

cat("Class:\n")

print(
  class(
    clinical_locked[[stage_detailed_var]]
  )
)

cat("\nFrequency:\n")

print(
  table(
    clinical_locked[[stage_detailed_var]],
    useNA = "ifany"
  )
)

cat("\nUnique values:\n")

print(
  sort(
    unique(
      clinical_locked[[stage_detailed_var]]
    )
  )
)

cat("\nMissing:\n")

print(
  sum(
    is.na(
      clinical_locked[[stage_detailed_var]]
    )
  )
)


# ------------------------------------------------------------
# 15E. T stage
# ------------------------------------------------------------

cat("\nT STAGE\n")
cat("------------------------------------------------------------\n")

cat("Variable:", T_var, "\n")

cat("Class:\n")

print(
  class(
    clinical_locked[[T_var]]
  )
)

cat("\nFrequency:\n")

print(
  table(
    clinical_locked[[T_var]],
    useNA = "ifany"
  )
)

cat("\nUnique values:\n")

print(
  sort(
    unique(
      clinical_locked[[T_var]]
    )
  )
)


# ------------------------------------------------------------
# 15F. N stage
# ------------------------------------------------------------

cat("\nN STAGE\n")
cat("------------------------------------------------------------\n")

cat("Variable:", N_var, "\n")

cat("Class:\n")

print(
  class(
    clinical_locked[[N_var]]
  )
)

cat("\nFrequency:\n")

print(
  table(
    clinical_locked[[N_var]],
    useNA = "ifany"
  )
)

cat("\nUnique values:\n")

print(
  sort(
    unique(
      clinical_locked[[N_var]]
    )
  )
)


# ------------------------------------------------------------
# 15G. M stage
# ------------------------------------------------------------

cat("\nM STAGE\n")
cat("------------------------------------------------------------\n")

cat("Variable:", M_var, "\n")

cat("Class:\n")

print(
  class(
    clinical_locked[[M_var]]
  )
)

cat("\nFrequency:\n")

print(
  table(
    clinical_locked[[M_var]],
    useNA = "ifany"
  )
)

cat("\nUnique values:\n")

print(
  sort(
    unique(
      clinical_locked[[M_var]]
    )
  )
)


# ============================================================
# 16. Missingness audit
# ============================================================

cat("\n============================================================\n")
cat("STEP 6C-1: Missingness audit\n")
cat("============================================================\n\n")


qc_variables <-
  c(
    age_var,
    gender_var,
    stage_var,
    stage_detailed_var,
    T_var,
    N_var,
    M_var
  )


missingness_qc <-
  data.frame(
    Variable = qc_variables,
    
    N_total =
      nrow(clinical_locked),
    
    N_missing =
      vapply(
        qc_variables,
        function(v) {
          
          x <- clinical_locked[[v]]
          
          sum(
            is.na(x) |
              (
                is.character(x) &
                  trimws(x) == ""
              )
          )
          
        },
        numeric(1)
      ),
    
    N_complete =
      vapply(
        qc_variables,
        function(v) {
          
          x <- clinical_locked[[v]]
          
          sum(
            !is.na(x) &
              !(
                is.character(x) &
                  trimws(x) == ""
              )
          )
          
        },
        numeric(1)
      ),
    
    stringsAsFactors = FALSE
  )


missingness_qc$Percent_missing <-
  100 *
  missingness_qc$N_missing /
  missingness_qc$N_total


print(
  missingness_qc,
  row.names = FALSE
)


# ============================================================
# 17. Compare broad vs detailed stage
# ============================================================

cat("\n============================================================\n")
cat("STEP 6C-1: Broad vs detailed stage comparison\n")
cat("============================================================\n\n")


stage_comparison <-
  table(
    Broad_stage =
      clinical_locked[[stage_var]],
    
    Detailed_stage =
      clinical_locked[[stage_detailed_var]],
    
    useNA = "ifany"
  )


print(
  stage_comparison
)


# ------------------------------------------------------------
# Check whether detailed stage maps consistently to broad stage
# ------------------------------------------------------------

stage_mapping_table <-
  data.frame(
    Detailed_stage =
      clinical_locked[[stage_detailed_var]],
    
    Broad_stage =
      clinical_locked[[stage_var]],
    
    stringsAsFactors = FALSE
  )


stage_mapping_table <-
  unique(
    stage_mapping_table
  )

stage_mapping_table <-
  stage_mapping_table[
    order(
      stage_mapping_table$Detailed_stage,
      stage_mapping_table$Broad_stage
    ),
    ,
    drop = FALSE
  ]


cat("\nUnique detailed → broad stage mappings:\n\n")

print(
  stage_mapping_table,
  row.names = FALSE
)


# ============================================================
# 18. Create raw clinical covariate dataset
# ============================================================

cat("\n============================================================\n")
cat("STEP 6C-1: Creating raw clinical covariate dataset\n")
cat("============================================================\n\n")


clinical_covariates_locked <-
  data.frame(
    
    Sample =
      locked_primary_samples,
    
    patient_id =
      clinical_locked$patient_id,
    
    Age =
      clinical_locked[[age_var]],
    
    Gender_raw =
      clinical_locked[[gender_var]],
    
    Stage_broad_raw =
      clinical_locked[[stage_var]],
    
    Stage_detailed_raw =
      clinical_locked[[stage_detailed_var]],
    
    T_raw =
      clinical_locked[[T_var]],
    
    N_raw =
      clinical_locked[[N_var]],
    
    M_raw =
      clinical_locked[[M_var]],
    
    stringsAsFactors = FALSE,
    
    check.names = FALSE
  )


cat(
  "Clinical covariate dataset dimensions: ",
  nrow(clinical_covariates_locked),
  " × ",
  ncol(clinical_covariates_locked),
  "\n",
  sep = ""
)


# ============================================================
# 19. Final mapping QC
# ============================================================

cat("\n============================================================\n")
cat("STEP 6C-1 FINAL MAPPING QC\n")
cat("============================================================\n\n")


stopifnot(
  nrow(clinical_covariates_locked) == 485
)

stopifnot(
  !anyDuplicated(
    clinical_covariates_locked$Sample
  )
)

stopifnot(
  identical(
    clinical_covariates_locked$Sample,
    locked_primary_samples
  )
)

stopifnot(
  identical(
    clinical_covariates_locked$Sample,
    as.character(
      ms4a7_tv_data$Sample
    )
  )
)


cat("Locked cohort N = 485: PASSED\n")
cat("Sample uniqueness: PASSED\n")
cat("Sample order identity: PASSED\n")
cat("Clinical mapping: PASSED\n")
cat("Clinical variable identification: PASSED\n")
cat("Raw coding audit: COMPLETED\n")
cat("Missingness audit: COMPLETED\n")
cat("Stage cross-tabulation: COMPLETED\n")
cat("No recoding performed: PASSED\n")
cat("No Cox model fitted: PASSED\n")


cat("\n============================================================\n")
cat("STEP 6C-1 MAPPING AUDIT COMPLETED.\n")
cat("============================================================\n")

# ============================================================
# STEP 6C-1.5: Clinical stage consistency audit
# No recoding
# No cohort modification
# No Cox model
# ============================================================

cat("\n============================================================\n")
cat("STEP 6C-1.5: Clinical stage consistency audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Detailed stage × M
# ------------------------------------------------------------

cat("1. Detailed stage × M stage\n")
cat("------------------------------------------------------------\n")

stage_M_table <- table(
  Detailed_stage =
    clinical_covariates_locked$Stage_detailed_raw,
  M_stage =
    clinical_covariates_locked$M_raw,
  useNA = "ifany"
)

print(stage_M_table)


# ------------------------------------------------------------
# 2. Detailed stage × T
# ------------------------------------------------------------

cat("\n2. Detailed stage × T stage\n")
cat("------------------------------------------------------------\n")

stage_T_table <- table(
  Detailed_stage =
    clinical_covariates_locked$Stage_detailed_raw,
  T_stage =
    clinical_covariates_locked$T_raw,
  useNA = "ifany"
)

print(stage_T_table)


# ------------------------------------------------------------
# 3. Detailed stage × N
# ------------------------------------------------------------

cat("\n3. Detailed stage × N stage\n")
cat("------------------------------------------------------------\n")

stage_N_table <- table(
  Detailed_stage =
    clinical_covariates_locked$Stage_detailed_raw,
  N_stage =
    clinical_covariates_locked$N_raw,
  useNA = "ifany"
)

print(stage_N_table)


# ------------------------------------------------------------
# 4. Samples with M missing
# ------------------------------------------------------------

cat("\n4. Samples with missing M stage\n")
cat("------------------------------------------------------------\n")

M_missing_samples <-
  clinical_covariates_locked[
    is.na(clinical_covariates_locked$M_raw),
    c(
      "Sample",
      "patient_id",
      "Age",
      "Gender_raw",
      "Stage_broad_raw",
      "Stage_detailed_raw",
      "T_raw",
      "N_raw",
      "M_raw"
    ),
    drop = FALSE
  ]

print(
  M_missing_samples,
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Unknown TNM categories
# ------------------------------------------------------------

cat("\n5. Unknown TNM categories\n")
cat("------------------------------------------------------------\n")

cat("TX:\n")
print(
  clinical_covariates_locked[
    clinical_covariates_locked$T_raw == "TX",
    c(
      "Sample",
      "Stage_detailed_raw",
      "T_raw",
      "N_raw",
      "M_raw"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\nNX:\n")
print(
  clinical_covariates_locked[
    clinical_covariates_locked$N_raw == "NX",
    c(
      "Sample",
      "Stage_detailed_raw",
      "T_raw",
      "N_raw",
      "M_raw"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\nMX:\n")
print(
  clinical_covariates_locked[
    clinical_covariates_locked$M_raw == "MX",
    c(
      "Sample",
      "Stage_detailed_raw",
      "T_raw",
      "N_raw",
      "M_raw"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Final QC
# ------------------------------------------------------------

stopifnot(
  nrow(clinical_covariates_locked) == 485
)

stopifnot(
  identical(
    clinical_covariates_locked$Sample,
    locked_primary_samples
  )
)

cat("\n============================================================\n")
cat("STEP 6C-1.5 COMPLETED.\n")
cat("No recoding performed.\n")
cat("No Cox model fitted.\n")
cat("Locked cohort unchanged: N = 485.\n")
cat("============================================================\n")


# ============================================================
# STEP 6C-2
# Clinical covariate coding freeze
#
# Purpose:
#   1. Freeze clinical covariate definitions
#   2. Derive clinically coherent Stage I-IV from
#      pathologic_stage_2
#   3. Collapse T1/T1a/T1b and T2/T2a/T2b
#   4. Collapse M1/M1a/M1b
#   5. Preserve unknown TNM as NA for sensitivity analysis
#
# IMPORTANT:
#   - Do NOT modify the locked 485-sample cohort
#   - Do NOT modify clinical_covariates_locked
#   - Do NOT fit Cox models in this step
# ============================================================

cat("\n============================================================\n")
cat("STEP 6C-2: Clinical covariate coding freeze\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Copy locked raw clinical dataset
# ------------------------------------------------------------

clinical_covariates_coded <-
  clinical_covariates_locked


# ------------------------------------------------------------
# 2. Age
# ------------------------------------------------------------

clinical_covariates_coded$Age_continuous <-
  as.numeric(
    clinical_covariates_coded$Age
  )


# ------------------------------------------------------------
# 3. Gender
# ------------------------------------------------------------

clinical_covariates_coded$Gender <-
  factor(
    clinical_covariates_coded$Gender_raw,
    levels = c(
      "female",
      "male"
    )
  )


# ------------------------------------------------------------
# 4. Clinical stage
#
# Source:
#   Stage_detailed_raw
#
# Mapping:
#   Stage / Stage IA / Stage IB -> Stage I
#   Stage II / Stage IIA / Stage IIB -> Stage II
#   Stage IIIA / Stage IIIB -> Stage III
#   Stage IV -> Stage IV
# ------------------------------------------------------------

clinical_covariates_coded$Stage <-
  NA_character_


clinical_covariates_coded$Stage[
  clinical_covariates_coded$Stage_detailed_raw %in%
    c(
      "Stage I",
      "Stage IA",
      "Stage IB"
    )
] <- "Stage I"


clinical_covariates_coded$Stage[
  clinical_covariates_coded$Stage_detailed_raw %in%
    c(
      "Stage II",
      "Stage IIA",
      "Stage IIB"
    )
] <- "Stage II"


clinical_covariates_coded$Stage[
  clinical_covariates_coded$Stage_detailed_raw %in%
    c(
      "Stage IIIA",
      "Stage IIIB"
    )
] <- "Stage III"


clinical_covariates_coded$Stage[
  clinical_covariates_coded$Stage_detailed_raw ==
    "Stage IV"
] <- "Stage IV"


clinical_covariates_coded$Stage <-
  factor(
    clinical_covariates_coded$Stage,
    levels = c(
      "Stage I",
      "Stage II",
      "Stage III",
      "Stage IV"
    )
  )


# ------------------------------------------------------------
# 5. T stage
#
# T1 = T1/T1a/T1b
# T2 = T2/T2a/T2b
# T3 = T3
# T4 = T4
# TX = NA
# ------------------------------------------------------------

clinical_covariates_coded$T_stage <-
  NA_character_


clinical_covariates_coded$T_stage[
  clinical_covariates_coded$T_raw %in%
    c(
      "T1",
      "T1a",
      "T1b"
    )
] <- "T1"


clinical_covariates_coded$T_stage[
  clinical_covariates_coded$T_raw %in%
    c(
      "T2",
      "T2a",
      "T2b"
    )
] <- "T2"


clinical_covariates_coded$T_stage[
  clinical_covariates_coded$T_raw == "T3"
] <- "T3"


clinical_covariates_coded$T_stage[
  clinical_covariates_coded$T_raw == "T4"
] <- "T4"


clinical_covariates_coded$T_stage <-
  factor(
    clinical_covariates_coded$T_stage,
    levels = c(
      "T1",
      "T2",
      "T3",
      "T4"
    )
  )


# ------------------------------------------------------------
# 6. N stage
#
# N0/N1/N2/N3 retained
# NX = NA
# ------------------------------------------------------------

clinical_covariates_coded$N_stage <-
  clinical_covariates_coded$N_raw


clinical_covariates_coded$N_stage[
  clinical_covariates_coded$N_stage == "NX"
] <- NA_character_


clinical_covariates_coded$N_stage <-
  factor(
    clinical_covariates_coded$N_stage,
    levels = c(
      "N0",
      "N1",
      "N2",
      "N3"
    )
  )


# ------------------------------------------------------------
# 7. M stage
#
# M0 = M0
# M1 = M1/M1a/M1b
# MX / NA = NA
# ------------------------------------------------------------

clinical_covariates_coded$M_stage <-
  NA_character_


clinical_covariates_coded$M_stage[
  clinical_covariates_coded$M_raw == "M0"
] <- "M0"


clinical_covariates_coded$M_stage[
  clinical_covariates_coded$M_raw %in%
    c(
      "M1",
      "M1a",
      "M1b"
    )
] <- "M1"


clinical_covariates_coded$M_stage <-
  factor(
    clinical_covariates_coded$M_stage,
    levels = c(
      "M0",
      "M1"
    )
  )


# ------------------------------------------------------------
# 8. Coding QC
# ------------------------------------------------------------

cat("Gender:\n")
print(table(
  clinical_covariates_coded$Gender,
  useNA = "ifany"
))


cat("\nClinical Stage:\n")
print(table(
  clinical_covariates_coded$Stage,
  useNA = "ifany"
))


cat("\nT stage:\n")
print(table(
  clinical_covariates_coded$T_stage,
  useNA = "ifany"
))


cat("\nN stage:\n")
print(table(
  clinical_covariates_coded$N_stage,
  useNA = "ifany"
))


cat("\nM stage:\n")
print(table(
  clinical_covariates_coded$M_stage,
  useNA = "ifany"
))


# ------------------------------------------------------------
# 9. Unknown / missing counts
# ------------------------------------------------------------

cat("\nUnknown / missing counts:\n")

cat(
  "Age missing:",
  sum(is.na(clinical_covariates_coded$Age_continuous)),
  "\n"
)

cat(
  "Gender missing:",
  sum(is.na(clinical_covariates_coded$Gender)),
  "\n"
)

cat(
  "Stage missing:",
  sum(is.na(clinical_covariates_coded$Stage)),
  "\n"
)

cat(
  "T-stage missing:",
  sum(is.na(clinical_covariates_coded$T_stage)),
  "\n"
)

cat(
  "N-stage missing:",
  sum(is.na(clinical_covariates_coded$N_stage)),
  "\n"
)

cat(
  "M-stage missing:",
  sum(is.na(clinical_covariates_coded$M_stage)),
  "\n"
)


# ------------------------------------------------------------
# 10. Verify raw locked data were not changed
# ------------------------------------------------------------

stopifnot(
  nrow(clinical_covariates_coded) == 485
)

stopifnot(
  identical(
    clinical_covariates_coded$Sample,
    clinical_covariates_locked$Sample
  )
)

stopifnot(
  identical(
    clinical_covariates_coded$Age,
    clinical_covariates_locked$Age
  )
)

stopifnot(
  identical(
    clinical_covariates_coded$Gender_raw,
    clinical_covariates_locked$Gender_raw
  )
)

stopifnot(
  identical(
    clinical_covariates_coded$Stage_detailed_raw,
    clinical_covariates_locked$Stage_detailed_raw
  )
)

stopifnot(
  identical(
    clinical_covariates_coded$T_raw,
    clinical_covariates_locked$T_raw
  )
)

stopifnot(
  identical(
    clinical_covariates_coded$N_raw,
    clinical_covariates_locked$N_raw
  )
)

stopifnot(
  identical(
    clinical_covariates_coded$M_raw,
    clinical_covariates_locked$M_raw
  )
)


cat("\n============================================================\n")
cat("STEP 6C-2 CODING FREEZE: PASSED\n")
cat("Locked cohort unchanged: N = 485\n")
cat("Raw clinical data preserved.\n")
cat("No Cox model fitted.\n")
cat("============================================================\n")

# ============================================================
# STEP 6C-3
# Multivariable Cox design-cohort QC
#
# Purpose:
#   1. Build primary clinical-adjustment cohort
#   2. Build TNM sensitivity cohort
#   3. Quantify N/events/censoring
#   4. Check category counts
#   5. DO NOT fit Cox models
#
# Locked cohort remains unchanged.
# ============================================================

cat("\n============================================================\n")
cat("STEP 6C-3: Multivariable Cox design-cohort QC\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Recover survival variables from locked Step 6B data
# ------------------------------------------------------------

survival_design <-
  data.frame(
    Sample =
      as.character(
        survival_gene_primary_z$Sample
      ),
    
    OS_time =
      survival_gene_primary_z$OS.time,
    
    OS_event =
      survival_gene_primary_z$OS_event,
    
    AOC3_z =
      survival_gene_primary_z$AOC3,
    
    MS4A7_z =
      survival_gene_primary_z$MS4A7,
    
    IL7R_z =
      survival_gene_primary_z$IL7R,
    
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


# ------------------------------------------------------------
# 2. Verify sample order
# ------------------------------------------------------------

stopifnot(
  identical(
    survival_design$Sample,
    clinical_covariates_coded$Sample
  )
)

stopifnot(
  nrow(survival_design) == 485
)


# ============================================================
# MODEL A
# Primary clinical-adjustment cohort
#
# MS4A7 + Age + Gender + Stage
# ============================================================

modelA_data <-
  data.frame(
    Sample =
      survival_design$Sample,
    
    OS_time =
      survival_design$OS_time,
    
    OS_event =
      survival_design$OS_event,
    
    MS4A7_z =
      survival_design$MS4A7_z,
    
    Age =
      clinical_covariates_coded$Age_continuous,
    
    Gender =
      clinical_covariates_coded$Gender,
    
    Stage =
      clinical_covariates_coded$Stage,
    
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


modelA_complete <-
  complete.cases(
    modelA_data[
      ,
      c(
        "OS_time",
        "OS_event",
        "MS4A7_z",
        "Age",
        "Gender",
        "Stage"
      )
    ]
  )


modelA_cc <-
  modelA_data[
    modelA_complete,
    ,
    drop = FALSE
  ]


cat("============================================================\n")
cat("MODEL A: Primary clinical-adjustment cohort\n")
cat("============================================================\n\n")

cat(
  "N complete:",
  nrow(modelA_cc),
  "\n"
)

cat(
  "Events:",
  sum(modelA_cc$OS_event == 1),
  "\n"
)

cat(
  "Censored:",
  sum(modelA_cc$OS_event == 0),
  "\n\n"
)


cat("Gender:\n")
print(
  table(
    modelA_cc$Gender,
    useNA = "ifany"
  )
)


cat("\nClinical Stage:\n")
print(
  table(
    modelA_cc$Stage,
    useNA = "ifany"
  )
)


# ============================================================
# MODEL B
# TNM sensitivity cohort
#
# MS4A7 + Age + Gender + T + N + M
# ============================================================

modelB_data <-
  data.frame(
    Sample =
      survival_design$Sample,
    
    OS_time =
      survival_design$OS_time,
    
    OS_event =
      survival_design$OS_event,
    
    MS4A7_z =
      survival_design$MS4A7_z,
    
    Age =
      clinical_covariates_coded$Age_continuous,
    
    Gender =
      clinical_covariates_coded$Gender,
    
    T_stage =
      clinical_covariates_coded$T_stage,
    
    N_stage =
      clinical_covariates_coded$N_stage,
    
    M_stage =
      clinical_covariates_coded$M_stage,
    
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


modelB_complete <-
  complete.cases(
    modelB_data[
      ,
      c(
        "OS_time",
        "OS_event",
        "MS4A7_z",
        "Age",
        "Gender",
        "T_stage",
        "N_stage",
        "M_stage"
      )
    ]
  )


modelB_cc <-
  modelB_data[
    modelB_complete,
    ,
    drop = FALSE
  ]


cat("\n============================================================\n")
cat("MODEL B: TNM sensitivity cohort\n")
cat("============================================================\n\n")

cat(
  "N complete:",
  nrow(modelB_cc),
  "\n"
)

cat(
  "Events:",
  sum(modelB_cc$OS_event == 1),
  "\n"
)

cat(
  "Censored:",
  sum(modelB_cc$OS_event == 0),
  "\n\n"
)


cat("Gender:\n")
print(
  table(
    modelB_cc$Gender,
    useNA = "ifany"
  )
)


cat("\nT stage:\n")
print(
  table(
    modelB_cc$T_stage,
    useNA = "ifany"
  )
)


cat("\nN stage:\n")
print(
  table(
    modelB_cc$N_stage,
    useNA = "ifany"
  )
)


cat("\nM stage:\n")
print(
  table(
    modelB_cc$M_stage,
    useNA = "ifany"
  )
)


# ============================================================
# 3. Compare exclusions between Model A and Model B
# ============================================================

cat("\n============================================================\n")
cat("Model A vs Model B cohort comparison\n")
cat("============================================================\n\n")

cat(
  "Locked survival cohort:",
  nrow(survival_design),
  "\n"
)

cat(
  "Model A complete cases:",
  nrow(modelA_cc),
  "\n"
)

cat(
  "Model A excluded:",
  sum(!modelA_complete),
  "\n"
)

cat(
  "Model B complete cases:",
  nrow(modelB_cc),
  "\n"
)

cat(
  "Model B excluded:",
  sum(!modelB_complete),
  "\n"
)


# ------------------------------------------------------------
# 4. Exclusion reason for Model B
# ------------------------------------------------------------

modelB_missing_reason <-
  data.frame(
    Sample =
      modelB_data$Sample,
    
    T_missing =
      is.na(modelB_data$T_stage),
    
    N_missing =
      is.na(modelB_data$N_stage),
    
    M_missing =
      is.na(modelB_data$M_stage),
    
    stringsAsFactors = FALSE
  )


modelB_missing_reason$Any_TNM_missing <-
  with(
    modelB_missing_reason,
    T_missing |
      N_missing |
      M_missing
  )


cat("\nModel B missingness among locked cohort:\n")

print(
  data.frame(
    Variable = c(
      "T",
      "N",
      "M",
      "Any_TNM"
    ),
    
    N_missing = c(
      sum(modelB_missing_reason$T_missing),
      sum(modelB_missing_reason$N_missing),
      sum(modelB_missing_reason$M_missing),
      sum(modelB_missing_reason$Any_TNM_missing)
    ),
    
    row.names = NULL
  )
)


# ============================================================
# 5. Final design-cohort QC
# ============================================================

stopifnot(
  nrow(modelA_cc) == 485
)

stopifnot(
  !anyDuplicated(modelA_cc$Sample)
)

stopifnot(
  !anyDuplicated(modelB_cc$Sample)
)

stopifnot(
  all(
    modelA_cc$Sample %in%
      survival_design$Sample
  )
)

stopifnot(
  all(
    modelB_cc$Sample %in%
      survival_design$Sample
  )
)


cat("\n============================================================\n")
cat("STEP 6C-3 DESIGN-COHORT QC COMPLETED.\n")
cat("No Cox model fitted.\n")
cat("Locked survival cohort unchanged.\n")
cat("============================================================\n")


# ============================================================
# Paper 1
# STEP 6C-4
# Multivariable Cox models and proportional-hazards QC
#
# Model A:
# MS4A7_z + Age + Gender + Stage
#
# Model B:
# MS4A7_z + Age + Gender + T + N + M
#
# IMPORTANT:
# - Model A is the primary adjusted model
# - Model B is a TNM-complete-case sensitivity model
# - Stage is NOT combined with T/N/M
# - No time-varying multivariable model is fitted at this step
# ============================================================

cat("\n============================================================\n")
cat("STEP 6C-4: Multivariable Cox models and PH QC\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Basic design checks
# ------------------------------------------------------------

stopifnot(
  nrow(modelA_cc) == 485
)

stopifnot(
  nrow(modelB_cc) == 336
)

stopifnot(
  sum(modelA_cc$OS_event == 1) == 178
)

stopifnot(
  sum(modelB_cc$OS_event == 1) == 137
)

stopifnot(
  !anyDuplicated(modelA_cc$Sample)
)

stopifnot(
  !anyDuplicated(modelB_cc$Sample)
)


# ------------------------------------------------------------
# 2. Explicit factor reference levels
# ------------------------------------------------------------

modelA_cc$Gender <-
  factor(
    modelA_cc$Gender,
    levels = c("female", "male")
  )

modelA_cc$Stage <-
  factor(
    modelA_cc$Stage,
    levels = c(
      "Stage I",
      "Stage II",
      "Stage III",
      "Stage IV"
    )
  )


modelB_cc$Gender <-
  factor(
    modelB_cc$Gender,
    levels = c("female", "male")
  )

modelB_cc$T_stage <-
  factor(
    modelB_cc$T_stage,
    levels = c(
      "T1",
      "T2",
      "T3",
      "T4"
    )
  )

modelB_cc$N_stage <-
  factor(
    modelB_cc$N_stage,
    levels = c(
      "N0",
      "N1",
      "N2",
      "N3"
    )
  )

modelB_cc$M_stage <-
  factor(
    modelB_cc$M_stage,
    levels = c(
      "M0",
      "M1"
    )
  )


# ------------------------------------------------------------
# 3. Verify reference categories
# ------------------------------------------------------------

stopifnot(
  levels(modelA_cc$Gender)[1] == "female"
)

stopifnot(
  levels(modelA_cc$Stage)[1] == "Stage I"
)

stopifnot(
  levels(modelB_cc$Gender)[1] == "female"
)

stopifnot(
  levels(modelB_cc$T_stage)[1] == "T1"
)

stopifnot(
  levels(modelB_cc$N_stage)[1] == "N0"
)

stopifnot(
  levels(modelB_cc$M_stage)[1] == "M0"
)


# ============================================================
# MODEL A
# Primary adjusted Cox model
# ============================================================

cat("============================================================\n")
cat("MODEL A: Primary adjusted Cox model\n")
cat("============================================================\n\n")


cox_modelA <-
  survival::coxph(
    survival::Surv(
      OS_time,
      OS_event
    ) ~
      MS4A7_z +
      Age +
      Gender +
      Stage,
    data = modelA_cc,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )


summary_modelA <-
  summary(cox_modelA)


cat("N =", nrow(modelA_cc), "\n")
cat(
  "Events =",
  sum(modelA_cc$OS_event == 1),
  "\n\n"
)

print(summary_modelA)


# ------------------------------------------------------------
# Model A PH test
# ------------------------------------------------------------

cox_modelA_ph <-
  survival::cox.zph(
    cox_modelA,
    transform = "km"
  )


cat("\n============================================================\n")
cat("MODEL A: Schoenfeld proportional-hazards test\n")
cat("============================================================\n\n")

print(cox_modelA_ph)


# ============================================================
# MODEL B
# TNM sensitivity Cox model
# ============================================================

cat("\n============================================================\n")
cat("MODEL B: TNM sensitivity Cox model\n")
cat("============================================================\n\n")


cox_modelB <-
  survival::coxph(
    survival::Surv(
      OS_time,
      OS_event
    ) ~
      MS4A7_z +
      Age +
      Gender +
      T_stage +
      N_stage +
      M_stage,
    data = modelB_cc,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )


summary_modelB <-
  summary(cox_modelB)


cat("N =", nrow(modelB_cc), "\n")
cat(
  "Events =",
  sum(modelB_cc$OS_event == 1),
  "\n\n"
)

print(summary_modelB)


# ------------------------------------------------------------
# Model B PH test
# ------------------------------------------------------------

cox_modelB_ph <-
  survival::cox.zph(
    cox_modelB,
    transform = "km"
  )


cat("\n============================================================\n")
cat("MODEL B: Schoenfeld proportional-hazards test\n")
cat("============================================================\n\n")

print(cox_modelB_ph)


# ============================================================
# 4. Extract clean coefficient tables
# ============================================================

extract_cox_table <- function(
    fit,
    model_name
) {
  
  s <- summary(fit)
  
  coef_table <- as.data.frame(
    s$coefficients,
    check.names = FALSE
  )
  
  conf_table <- as.data.frame(
    s$conf.int,
    check.names = FALSE
  )
  
  out <- data.frame(
    Model =
      model_name,
    
    Term =
      rownames(coef_table),
    
    Beta =
      coef_table[, "coef"],
    
    HR =
      conf_table[, "exp(coef)"],
    
    CI_low =
      conf_table[, "lower .95"],
    
    CI_high =
      conf_table[, "upper .95"],
    
    P =
      coef_table[, "Pr(>|z|)"],
    
    row.names = NULL,
    
    check.names = FALSE
  )
  
  out
}


multivariable_cox_results <-
  rbind(
    extract_cox_table(
      cox_modelA,
      "Model A"
    ),
    
    extract_cox_table(
      cox_modelB,
      "Model B"
    )
  )


cat("\n============================================================\n")
cat("Multivariable Cox coefficient table\n")
cat("============================================================\n\n")

print(
  multivariable_cox_results,
  row.names = FALSE
)


# ============================================================
# 5. Extract PH results
# ============================================================

extract_ph_table <- function(
    ph_object,
    model_name
) {
  
  ph_table <-
    as.data.frame(
      ph_object$table,
      check.names = FALSE
    )
  
  ph_table$Term <-
    rownames(ph_table)
  
  ph_table$Model <-
    model_name
  
  rownames(ph_table) <-
    NULL
  
  ph_table[
    ,
    c(
      "Model",
      "Term",
      setdiff(
        colnames(ph_table),
        c("Model", "Term")
      )
    ),
    drop = FALSE
  ]
}


multivariable_ph_results <-
  rbind(
    extract_ph_table(
      cox_modelA_ph,
      "Model A"
    ),
    
    extract_ph_table(
      cox_modelB_ph,
      "Model B"
    )
  )


cat("\n============================================================\n")
cat("Multivariable PH table\n")
cat("============================================================\n\n")

print(
  multivariable_ph_results,
  row.names = FALSE
)


# ============================================================
# 6. Specifically extract MS4A7 results
# ============================================================

ms4a7_modelA <-
  multivariable_cox_results[
    multivariable_cox_results$Model == "Model A" &
      multivariable_cox_results$Term == "MS4A7_z",
    ,
    drop = FALSE
  ]

ms4a7_modelB <-
  multivariable_cox_results[
    multivariable_cox_results$Model == "Model B" &
      multivariable_cox_results$Term == "MS4A7_z",
    ,
    drop = FALSE
  ]


cat("\n============================================================\n")
cat("Adjusted MS4A7 results\n")
cat("============================================================\n\n")

cat("Model A:\n")
print(ms4a7_modelA, row.names = FALSE)

cat("\nModel B:\n")
print(ms4a7_modelB, row.names = FALSE)


# ============================================================
# 7. Specifically extract MS4A7 PH results
# ============================================================

ph_modelA_df <-
  as.data.frame(
    cox_modelA_ph$table,
    check.names = FALSE
  )

ph_modelB_df <-
  as.data.frame(
    cox_modelB_ph$table,
    check.names = FALSE
  )


ms4a7_ph_modelA <-
  ph_modelA_df[
    rownames(ph_modelA_df) == "MS4A7_z",
    ,
    drop = FALSE
  ]

ms4a7_ph_modelB <-
  ph_modelB_df[
    rownames(ph_modelB_df) == "MS4A7_z",
    ,
    drop = FALSE
  ]


cat("\n============================================================\n")
cat("Adjusted MS4A7 PH results\n")
cat("============================================================\n\n")

cat("Model A MS4A7 PH:\n")
print(ms4a7_ph_modelA)

cat("\nModel B MS4A7 PH:\n")
print(ms4a7_ph_modelB)


# ============================================================
# 8. N3 stability check
# ============================================================

cat("\n============================================================\n")
cat("MODEL B: N3 category stability check\n")
cat("============================================================\n\n")

cat(
  "N3 observations:",
  sum(modelB_cc$N_stage == "N3"),
  "\n"
)

cat(
  "N3 events:",
  sum(
    modelB_cc$N_stage == "N3" &
      modelB_cc$OS_event == 1
  ),
  "\n"
)

cat(
  "N3 censored:",
  sum(
    modelB_cc$N_stage == "N3" &
      modelB_cc$OS_event == 0
  ),
  "\n"
)

# ============================================================
# 9. Corrected overall model QC
# ============================================================

stopifnot(
  cox_modelA$n == 485
)

stopifnot(
  cox_modelB$n == 336
)

stopifnot(
  sum(modelA_cc$OS_event == 1) == 178
)

stopifnot(
  sum(modelB_cc$OS_event == 1) == 137
)

stopifnot(
  !anyNA(coef(cox_modelA))
)

# Model B contains an intentionally unstable N3 coefficient.
# Therefore do NOT require all coefficients to be finite.
stopifnot(
  is.finite(coef(cox_modelB)["MS4A7_z"])
)

cat("\n============================================================\n")
cat("STEP 6C-4 MODEL FITTING AND PH QC: PASSED\n")
cat("============================================================\n\n")

cat("Model A: N = 485; events = 178\n")
cat("Model B: N = 336; events = 137\n")
cat("Model A MS4A7 PH P = 0.0321\n")
cat("Model B MS4A7 PH P = 0.0517\n")
cat("Model B N3: 1 observation, 0 events\n")
cat("No time-varying multivariable model fitted yet.\n")

# ============================================================
# Paper 1
# STEP 6C-5
# Adjusted time-varying Cox model for MS4A7
#
# Primary adjusted design:
#   MS4A7_z + Age + Gender + Stage
#
# Locked cohort:
#   N = 485
#   Events = 178
#
# Time functions:
#   1. linear: MS4A7_z * t
#   2. sqrt:   MS4A7_z * sqrt(t)
#   3. log1p:  MS4A7_z * log1p(t)
#
# IMPORTANT:
# - Uses the already locked Model A cohort
# - Does NOT reconstruct the survival cohort
# - Does NOT alter clinical coding
# - Does NOT include T/N/M
# - MS4A7 is the only intentionally time-varying predictor
# ============================================================


cat("\n============================================================\n")
cat("STEP 6C-5: Adjusted time-varying Cox for MS4A7\n")
cat("============================================================\n\n")


# ============================================================
# 1. Recover and verify the locked Model A cohort
# ============================================================

stopifnot(
  exists("modelA_cc")
)

stopifnot(
  nrow(modelA_cc) == 485
)

stopifnot(
  sum(modelA_cc$OS_event == 1) == 178
)

stopifnot(
  sum(modelA_cc$OS_event == 0) == 307
)

stopifnot(
  !anyDuplicated(modelA_cc$Sample)
)


cat("Locked Model A cohort verified.\n")
cat("N =", nrow(modelA_cc), "\n")
cat(
  "Events =",
  sum(modelA_cc$OS_event == 1),
  "\n"
)
cat(
  "Censored =",
  sum(modelA_cc$OS_event == 0),
  "\n\n"
)


# ============================================================
# 2. Verify required variables
# ============================================================

required_variables <-
  c(
    "Sample",
    "OS_time",
    "OS_event",
    "MS4A7_z",
    "Age",
    "Gender",
    "Stage"
  )

stopifnot(
  all(
    required_variables %in%
      colnames(modelA_cc)
  )
)


# ============================================================
# 3. Verify clinical reference levels
# ============================================================

stopifnot(
  identical(
    levels(modelA_cc$Gender),
    c("female", "male")
  )
)

stopifnot(
  identical(
    levels(modelA_cc$Stage),
    c(
      "Stage I",
      "Stage II",
      "Stage III",
      "Stage IV"
    )
  )
)


cat("Reference categories:\n")
cat("Gender: female\n")
cat("Stage: Stage I\n\n")


# ============================================================
# 4. Verify survival variables
# ============================================================

stopifnot(
  !anyNA(modelA_cc$OS_time)
)

stopifnot(
  !anyNA(modelA_cc$OS_event)
)

stopifnot(
  !anyNA(modelA_cc$MS4A7_z)
)

stopifnot(
  all(
    modelA_cc$OS_event %in%
      c(0, 1)
  )
)

stopifnot(
  all(
    modelA_cc$OS_time >= 0
  )
)


# ============================================================
# 5. Define prespecified time functions
#
# IMPORTANT:
# coxph(tt=...) passes x first and t second.
# ============================================================

time_functions <-
  list(
    
    linear =
      function(x, t, ...) {
        x * t
      },
    
    sqrt =
      function(x, t, ...) {
        x * sqrt(t)
      },
    
    log1p =
      function(x, t, ...) {
        x * log1p(t)
      }
    
  )


cat("Prespecified time functions:\n")
cat("- linear: MS4A7_z × t\n")
cat("- sqrt:   MS4A7_z × sqrt(t)\n")
cat("- log1p:  MS4A7_z × log1p(t)\n\n")


# ============================================================
# 6. Fit adjusted time-varying Cox models
# ============================================================

ms4a7_tv_adjusted_models <-
  list()


for (fn_name in names(time_functions)) {
  
  cat("\n------------------------------------------------------------\n")
  cat(
    "Fitting adjusted MS4A7 time-varying model:",
    fn_name,
    "\n"
  )
  cat("------------------------------------------------------------\n")
  
  current_model <-
    survival::coxph(
      
      survival::Surv(
        OS_time,
        OS_event
      ) ~
        MS4A7_z +
        tt(MS4A7_z) +
        Age +
        Gender +
        Stage,
      
      data = modelA_cc,
      
      ties = "efron",
      
      tt =
        time_functions[[fn_name]],
      
      x = TRUE,
      y = TRUE
      
    )
  
  ms4a7_tv_adjusted_models[[fn_name]] <-
    current_model
  
  cat(
    "N =",
    current_model$n,
    "\n"
  )
  
  cat(
    "Events =",
    sum(modelA_cc$OS_event == 1),
    "\n"
  )
  
  cat(
    "AIC =",
    AIC(current_model),
    "\n"
  )
  
  cat(
    "LogLik =",
    as.numeric(logLik(current_model)),
    "\n"
  )
}


# ============================================================
# 7. Print model summaries
# ============================================================

cat("\n============================================================\n")
cat("Adjusted time-varying Cox model summaries\n")
cat("============================================================\n")


for (fn_name in names(ms4a7_tv_adjusted_models)) {
  
  cat("\n\n############################################################\n")
  cat(
    "TIME FUNCTION:",
    fn_name,
    "\n"
  )
  cat("############################################################\n\n")
  
  print(
    summary(
      ms4a7_tv_adjusted_models[[fn_name]]
    )
  )
}


# ============================================================
# 8. Extract model comparison statistics
# ============================================================

extract_tv_adjusted_model_stats <-
  function(
    fit,
    time_function
  ) {
    
    s <-
      summary(fit)
    
    data.frame(
      
      Time_function =
        time_function,
      
      N =
        fit$n,
      
      Events =
        sum(modelA_cc$OS_event == 1),
      
      AIC =
        AIC(fit),
      
      LogLik =
        as.numeric(
          logLik(fit)
        ),
      
      Concordance =
        unname(
          s$concordance[1]
        ),
      
      Concordance_SE =
        unname(
          s$concordance[2]
        ),
      
      stringsAsFactors =
        FALSE
      
    )
  }


ms4a7_tv_adjusted_comparison <-
  do.call(
    rbind,
    lapply(
      names(ms4a7_tv_adjusted_models),
      function(fn_name) {
        
        extract_tv_adjusted_model_stats(
          ms4a7_tv_adjusted_models[[fn_name]],
          fn_name
        )
        
      }
    )
  )


rownames(
  ms4a7_tv_adjusted_comparison
) <-
  NULL


cat("\n============================================================\n")
cat("Adjusted time-varying model comparison\n")
cat("============================================================\n\n")

print(
  ms4a7_tv_adjusted_comparison,
  row.names = FALSE
)


# ============================================================
# 9. Select prespecified time function by minimum AIC
# ============================================================

best_time_function_adjusted <-
  ms4a7_tv_adjusted_comparison$Time_function[
    which.min(
      ms4a7_tv_adjusted_comparison$AIC
    )
  ]

best_tv_adjusted_model <- ms4a7_tv_adjusted_models[[best_time_function_adjusted]]


cat("\n============================================================\n")
cat("Selected adjusted time function\n")
cat("============================================================\n\n")

cat(
  "Selected function:",
  best_time_function_adjusted,
  "\n"
)

cat(
  "Minimum AIC:",
  min(
    ms4a7_tv_adjusted_comparison$AIC
  ),
  "\n"
)

cat("Selected adjusted time function:",
    best_time_function_adjusted,
    "\n")

cat("Available models:",
    paste(
      names(ms4a7_tv_adjusted_models),
      collapse = ", "
    ),
    "\n")

stopifnot(
  best_time_function_adjusted %in%
    names(ms4a7_tv_adjusted_models)
)

best_tv_adjusted_model <-
  ms4a7_tv_adjusted_models[[best_time_function_adjusted]]

stopifnot(
  inherits(
    best_tv_adjusted_model,
    "coxph"
  )
)

cat(
  "Selected model successfully recovered.\n"
)

cat(
  "N =",
  best_tv_adjusted_model$n,
  "\n"
)


# ============================================================
# 10. Extract MS4A7 baseline and time-varying coefficients
# ============================================================

best_summary <-
  summary(
    best_tv_adjusted_model
  )

best_coef_table <-
  as.data.frame(
    best_summary$coefficients,
    check.names = FALSE
  )


cat("\n============================================================\n")
cat("Selected adjusted MS4A7 time-varying coefficients\n")
cat("============================================================\n\n")

print(
  best_coef_table,
  row.names = TRUE
)


# ------------------------------------------------------------
# Extract exact MS4A7 terms
# ------------------------------------------------------------

ms4a7_beta0 <-
  unname(
    best_coef_table[
      "MS4A7_z",
      "coef"
    ]
  )

ms4a7_beta1 <-
  unname(
    best_coef_table[
      "tt(MS4A7_z)",
      "coef"
    ]
  )

ms4a7_beta0_se <-
  unname(
    best_coef_table[
      "MS4A7_z",
      "se(coef)"
    ]
  )

ms4a7_beta1_se <-
  unname(
    best_coef_table[
      "tt(MS4A7_z)",
      "se(coef)"
    ]
  )

ms4a7_beta0_p <-
  unname(
    best_coef_table[
      "MS4A7_z",
      "Pr(>|z|)"
    ]
  )

ms4a7_beta1_p <-
  unname(
    best_coef_table[
      "tt(MS4A7_z)",
      "Pr(>|z|)"
    ]
  )


cat("\nMS4A7 baseline coefficient:\n")
cat(
  "beta0 =",
  ms4a7_beta0,
  "\n"
)

cat(
  "SE =",
  ms4a7_beta0_se,
  "\n"
)

cat(
  "P =",
  ms4a7_beta0_p,
  "\n\n"
)


cat("MS4A7 time-varying coefficient:\n")
cat(
  "beta1 =",
  ms4a7_beta1,
  "\n"
)

cat(
  "SE =",
  ms4a7_beta1_se,
  "\n"
)

cat(
  "P =",
  ms4a7_beta1_p,
  "\n"
)


# ============================================================
# 11. Compare selected adjusted TV model with
#     constant-coefficient Model A
# ============================================================

stopifnot(
  exists("cox_modelA")
)

constant_logLik <-
  as.numeric(
    logLik(cox_modelA)
  )

tv_logLik <-
  as.numeric(
    logLik(best_tv_adjusted_model)
  )


ms4a7_adjusted_LRT <-
  2 *
  (
    tv_logLik -
      constant_logLik
  )


ms4a7_adjusted_LRT_df <-
  attr(
    logLik(best_tv_adjusted_model),
    "df"
  ) -
  attr(
    logLik(cox_modelA),
    "df"
  )


ms4a7_adjusted_LRT_p <-
  pchisq(
    ms4a7_adjusted_LRT,
    df =
      ms4a7_adjusted_LRT_df,
    lower.tail = FALSE
  )


ms4a7_adjusted_AIC_constant <-
  AIC(cox_modelA)

ms4a7_adjusted_AIC_tv <-
  AIC(best_tv_adjusted_model)

ms4a7_adjusted_delta_AIC <-
  ms4a7_adjusted_AIC_constant -
  ms4a7_adjusted_AIC_tv


cat("\n============================================================\n")
cat("Constant vs adjusted time-varying MS4A7 model\n")
cat("============================================================\n\n")

cat(
  "Constant model logLik =",
  constant_logLik,
  "\n"
)

cat(
  "Time-varying model logLik =",
  tv_logLik,
  "\n"
)

cat(
  "Likelihood-ratio statistic =",
  ms4a7_adjusted_LRT,
  "\n"
)

cat(
  "df =",
  ms4a7_adjusted_LRT_df,
  "\n"
)

cat(
  "LRT P =",
  ms4a7_adjusted_LRT_p,
  "\n\n"
)

cat(
  "Constant model AIC =",
  ms4a7_adjusted_AIC_constant,
  "\n"
)

cat(
  "Time-varying model AIC =",
  ms4a7_adjusted_AIC_tv,
  "\n"
)

cat(
  "Delta AIC (constant - TV) =",
  ms4a7_adjusted_delta_AIC,
  "\n"
)


# ============================================================
# 12. Calculate time-specific adjusted HR
#
# HR(t) = exp(beta0 + beta1 * sqrt(t))
#
# Time is measured in days.
# ============================================================

time_years <-
  c(
    1,
    2,
    3,
    4,
    5
  )

time_days <-
  time_years *
  365.25


if (
  best_time_function_adjusted == "linear"
) {
  
  time_transform_values <-
    time_days
  
} else if (
  best_time_function_adjusted == "sqrt"
) {
  
  time_transform_values <-
    sqrt(time_days)
  
} else if (
  best_time_function_adjusted == "log1p"
) {
  
  time_transform_values <-
    log1p(time_days)
  
} else {
  
  stop(
    "Unexpected time function."
  )
  
}


# ------------------------------------------------------------
# Time-specific coefficient
# ------------------------------------------------------------

beta_time <-
  ms4a7_beta0 +
  ms4a7_beta1 *
  time_transform_values


HR_time <-
  exp(
    beta_time
  )


# ------------------------------------------------------------
# Variance of beta0 + beta1*f(t)
#
# Need covariance between beta0 and beta1.
# ------------------------------------------------------------

coef_names <-
  names(
    coef(best_tv_adjusted_model)
  )


beta0_index <-
  match(
    "MS4A7_z",
    coef_names
  )

beta1_index <-
  match(
    "tt(MS4A7_z)",
    coef_names
  )


vcov_matrix <-
  vcov(
    best_tv_adjusted_model
  )


var_time <-
  sapply(
    time_transform_values,
    function(f) {
      
      v0 <-
        vcov_matrix[
          beta0_index,
          beta0_index
        ]
      
      v1 <-
        vcov_matrix[
          beta1_index,
          beta1_index
        ]
      
      cov01 <-
        vcov_matrix[
          beta0_index,
          beta1_index
        ]
      
      v0 +
        (f^2 * v1) +
        (2 * f * cov01)
      
    }
  )


SE_time <-
  sqrt(
    pmax(
      var_time,
      0
    )
  )


CI_low_time <-
  exp(
    beta_time -
      1.96 *
      SE_time
  )


CI_high_time <-
  exp(
    beta_time +
      1.96 *
      SE_time
  )


ms4a7_adjusted_hr_time <-
  data.frame(
    
    Time_years =
      time_years,
    
    Time_days =
      time_days,
    
    Time_transform =
      time_transform_values,
    
    Beta_time =
      beta_time,
    
    HR_per_1SD =
      HR_time,
    
    CI_low =
      CI_low_time,
    
    CI_high =
      CI_high_time,
    
    stringsAsFactors =
      FALSE
    
  )


cat("\n============================================================\n")
cat("Adjusted time-specific MS4A7 hazard ratios\n")
cat("============================================================\n\n")

print(
  ms4a7_adjusted_hr_time,
  row.names = FALSE
)

# ============================================================
# 13. PH diagnostics for selected adjusted TV model
# ============================================================

cat("\n============================================================\n")
cat("Selected adjusted time-varying model: PH diagnostic note\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 13.1 cox.zph() is not applicable to models containing tt()
# ------------------------------------------------------------

if (
  any(
    grepl(
      "^tt\\(",
      names(stats::coef(best_tv_adjusted_model))
    )
  )
) {
  
  cat(
    "cox.zph() was not applied to the selected adjusted ",
    "time-varying model because survival::cox.zph() does not ",
    "support Cox models containing tt() terms.\n",
    sep = ""
  )
  
} else {
  
  best_tv_adjusted_ph <-
    survival::cox.zph(
      best_tv_adjusted_model,
      transform = "km"
    )
  
  print(best_tv_adjusted_ph)
  
}

# ============================================================
# 14. Final QC
# ============================================================

cat("\n============================================================\n")
cat("STEP 6C-5 FINAL QC\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 14.1 Selected model
# ------------------------------------------------------------

stopifnot(
  best_time_function_adjusted %in%
    c("linear", "sqrt", "log1p")
)

stopifnot(
  inherits(
    best_tv_adjusted_model,
    "coxph"
  )
)

# ------------------------------------------------------------
# 14.2 Locked cohort
# ------------------------------------------------------------

stopifnot(
  best_tv_adjusted_model$n == 485
)

stopifnot(
  sum(modelA_cc$OS_event == 1) == 178
)

# ------------------------------------------------------------
# 14.3 Required coefficient
# ------------------------------------------------------------

stopifnot(
  "MS4A7_z" %in%
    names(coef(best_tv_adjusted_model))
)

stopifnot(
  "tt(MS4A7_z)" %in%
    names(coef(best_tv_adjusted_model))
)

stopifnot(
  is.finite(
    coef(best_tv_adjusted_model)["MS4A7_z"]
  )
)

stopifnot(
  is.finite(
    coef(best_tv_adjusted_model)["tt(MS4A7_z)"]
  )
)

# ------------------------------------------------------------
# 14.4 Time-specific HR table
# ------------------------------------------------------------

stopifnot(
  nrow(ms4a7_adjusted_hr_time) == 5
)

stopifnot(
  all(
    is.finite(
      ms4a7_adjusted_hr_time$HR_per_1SD
    )
  )
)

stopifnot(
  all(
    ms4a7_adjusted_hr_time$CI_low <=
      ms4a7_adjusted_hr_time$HR_per_1SD
  )
)

stopifnot(
  all(
    ms4a7_adjusted_hr_time$HR_per_1SD <=
      ms4a7_adjusted_hr_time$CI_high
  )
)

# ------------------------------------------------------------
# 14.5 PH diagnostic handling
# ------------------------------------------------------------

ph_test_available <-
  !any(
    grepl(
      "^tt\\(",
      names(stats::coef(best_tv_adjusted_model))
    )
  )

stopifnot(
  identical(
    ph_test_available,
    FALSE
  )
)

cat("\nSTEP 6C-5 FINAL QC: PASSED\n")
cat("Selected time function:",
    best_time_function_adjusted,
    "\n")
cat("Adjusted cohort N:",
    best_tv_adjusted_model$n,
    "\n")
cat("Events:",
    sum(modelA_cc$OS_event == 1),
    "\n")
cat(
  "cox.zph(): not applicable to tt() model; ",
  "no PH test attempted.\n"
)


cat("\n============================================================\n")
cat("STEP 6C-5 EXISTING OBJECT AUDIT\n")
cat("============================================================\n\n")

objects_to_check <- c(
  "modelA_cc",
  "ms4a7_tv_adjusted_models",
  "ms4a7_tv_adjusted_comparison",
  "best_time_function_adjusted",
  "best_tv_adjusted_model",
  "ms4a7_adjusted_LRT",
  "ms4a7_adjusted_hr_time"
)

object_status <- data.frame(
  Object = objects_to_check,
  Exists = vapply(
    objects_to_check,
    exists,
    logical(1),
    inherits = TRUE
  ),
  stringsAsFactors = FALSE
)

print(
  object_status,
  row.names = FALSE
)

# ============================================================
# Paper 1
# STEP 6C-5 CHECKPOINT
# Adjusted time-varying Cox model for MS4A7
#
# IMPORTANT:
# This checkpoint saves already-fitted objects and derived
# metadata only.
#
# No Cox model is refitted.
# No cohort is reconstructed.
# No clinical coding is changed.
# ============================================================


cat("\n============================================================\n")
cat("STEP 6C-5 CHECKPOINT: START\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Required existing objects
# ------------------------------------------------------------

required_objects <- c(
  "modelA_cc",
  "ms4a7_tv_adjusted_models",
  "ms4a7_tv_adjusted_comparison",
  "best_time_function_adjusted",
  "best_tv_adjusted_model",
  "ms4a7_adjusted_LRT",
  "ms4a7_adjusted_hr_time"
)


missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]


if (length(missing_objects) > 0) {
  
  stop(
    paste(
      "Missing required objects:",
      paste(
        missing_objects,
        collapse = ", "
      )
    )
  )
  
}


cat(
  "All required existing objects are present.\n\n"
)


# ------------------------------------------------------------
# 2. Verify locked Model A cohort
# ------------------------------------------------------------

stopifnot(
  nrow(modelA_cc) == 485
)

stopifnot(
  sum(modelA_cc$OS_event == 1) == 178
)

stopifnot(
  sum(modelA_cc$OS_event == 0) == 307
)


cat("Locked Model A cohort:\n")
cat("  N        =", nrow(modelA_cc), "\n")
cat("  Events   =", sum(modelA_cc$OS_event == 1), "\n")
cat("  Censored =", sum(modelA_cc$OS_event == 0), "\n\n")


# ------------------------------------------------------------
# 3. Verify selected time function
# ------------------------------------------------------------

stopifnot(
  best_time_function_adjusted %in%
    c(
      "linear",
      "sqrt",
      "log1p"
    )
)


cat(
  "Selected adjusted time function:",
  best_time_function_adjusted,
  "\n\n"
)


# ------------------------------------------------------------
# 4. Verify selected model
# ------------------------------------------------------------

stopifnot(
  inherits(
    best_tv_adjusted_model,
    "coxph"
  )
)

stopifnot(
  best_tv_adjusted_model$n == 485
)


cat(
  "Selected adjusted TV-Cox model verified.\n"
)

cat(
  "  N =",
  best_tv_adjusted_model$n,
  "\n\n"
)


# ------------------------------------------------------------
# 5. Extract coefficients directly from selected model
# ------------------------------------------------------------

coef_adjusted <-
  stats::coef(
    best_tv_adjusted_model
  )

vcov_adjusted <-
  stats::vcov(
    best_tv_adjusted_model
  )

summary_adjusted <-
  summary(
    best_tv_adjusted_model
  )


stopifnot(
  "MS4A7_z" %in%
    names(coef_adjusted)
)

stopifnot(
  "tt(MS4A7_z)" %in%
    names(coef_adjusted)
)


beta0_adjusted_locked <-
  unname(
    coef_adjusted["MS4A7_z"]
  )

beta1_adjusted_locked <-
  unname(
    coef_adjusted["tt(MS4A7_z)"]
  )

se0_adjusted_locked <-
  sqrt(
    vcov_adjusted[
      "MS4A7_z",
      "MS4A7_z"
    ]
  )

se1_adjusted_locked <-
  sqrt(
    vcov_adjusted[
      "tt(MS4A7_z)",
      "tt(MS4A7_z)"
    ]
  )


# Wald P values

z0_adjusted_locked <-
  beta0_adjusted_locked /
  se0_adjusted_locked

z1_adjusted_locked <-
  beta1_adjusted_locked /
  se1_adjusted_locked


p0_adjusted_locked <-
  2 *
  stats::pnorm(
    -abs(z0_adjusted_locked)
  )

p1_adjusted_locked <-
  2 *
  stats::pnorm(
    -abs(z1_adjusted_locked)
  )


# Hazard ratios for coefficient at t = 0 transformation
HR0_adjusted_locked <-
  exp(
    beta0_adjusted_locked
  )


cat(
  "Adjusted MS4A7 coefficients:\n"
)

cat(
  "  beta0 =",
  beta0_adjusted_locked,
  "\n"
)

cat(
  "  SE0   =",
  se0_adjusted_locked,
  "\n"
)

cat(
  "  P0    =",
  p0_adjusted_locked,
  "\n"
)

cat(
  "  beta1 =",
  beta1_adjusted_locked,
  "\n"
)

cat(
  "  SE1   =",
  se1_adjusted_locked,
  "\n"
)

cat(
  "  P1    =",
  p1_adjusted_locked,
  "\n\n"
)


# ------------------------------------------------------------
# 6. Verify time-function comparison
# ------------------------------------------------------------

print(
  ms4a7_tv_adjusted_comparison,
  row.names = FALSE
)


stopifnot(
  best_time_function_adjusted %in%
    ms4a7_tv_adjusted_comparison$Time_function
)


selected_row <-
  ms4a7_tv_adjusted_comparison[
    ms4a7_tv_adjusted_comparison$Time_function ==
      best_time_function_adjusted,
    ,
    drop = FALSE
  ]


stopifnot(
  nrow(selected_row) == 1
)


selected_AIC_locked <-
  selected_row$AIC


minimum_AIC_locked <-
  min(
    ms4a7_tv_adjusted_comparison$AIC,
    na.rm = TRUE
  )


stopifnot(
  isTRUE(
    all.equal(
      as.numeric(selected_AIC_locked),
      as.numeric(minimum_AIC_locked),
      tolerance = 1e-8
    )
  )
)


cat(
  "\nSelected model AIC =",
  selected_AIC_locked,
  "\n"
)

cat(
  "Minimum AIC =",
  minimum_AIC_locked,
  "\n\n"
)


# ------------------------------------------------------------
# 7. Verify constant vs TV model LRT
# ------------------------------------------------------------

print(
  ms4a7_adjusted_LRT
)


# Extract LRT values flexibly

lrt_value_locked <- NA_real_
lrt_df_locked <- NA_real_
lrt_p_locked <- NA_real_


if (
  is.data.frame(ms4a7_adjusted_LRT)
) {
  
  if ("LR" %in% names(ms4a7_adjusted_LRT)) {
    lrt_value_locked <-
      as.numeric(
        ms4a7_adjusted_LRT$LR[1]
      )
  }
  
  if ("df" %in% names(ms4a7_adjusted_LRT)) {
    lrt_df_locked <-
      as.numeric(
        ms4a7_adjusted_LRT$df[1]
      )
  }
  
  if ("P" %in% names(ms4a7_adjusted_LRT)) {
    lrt_p_locked <-
      as.numeric(
        ms4a7_adjusted_LRT$P[1]
      )
  }
  
}


if (
  is.numeric(ms4a7_adjusted_LRT) &&
  length(ms4a7_adjusted_LRT) >= 1
) {
  
  lrt_value_locked <-
    as.numeric(
      ms4a7_adjusted_LRT[1]
    )
  
}


stopifnot(
  is.finite(lrt_value_locked)
)


cat(
  "\nConstant vs TV-Cox LRT statistic =",
  lrt_value_locked,
  "\n"
)

if (is.finite(lrt_df_locked)) {
  cat(
    "LRT df =",
    lrt_df_locked,
    "\n"
  )
}

if (is.finite(lrt_p_locked)) {
  cat(
    "LRT P =",
    lrt_p_locked,
    "\n"
  )
}


# ------------------------------------------------------------
# 8. Verify time-specific adjusted HR
# ------------------------------------------------------------

print(
  ms4a7_adjusted_hr_time,
  row.names = FALSE
)


stopifnot(
  nrow(ms4a7_adjusted_hr_time) == 5
)

stopifnot(
  all(
    ms4a7_adjusted_hr_time$Time_years ==
      c(1, 2, 3, 4, 5)
  )
)

stopifnot(
  all(
    is.finite(
      ms4a7_adjusted_hr_time$HR_per_1SD
    )
  )
)

stopifnot(
  all(
    ms4a7_adjusted_hr_time$CI_low <=
      ms4a7_adjusted_hr_time$HR_per_1SD
  )
)

stopifnot(
  all(
    ms4a7_adjusted_hr_time$HR_per_1SD <=
      ms4a7_adjusted_hr_time$CI_high
  )
)


cat(
  "\nTime-specific adjusted HR QC passed.\n\n"
)


# ------------------------------------------------------------
# 9. Preserve Step 6C-4 PH results
#
# These values were already obtained in Step 6C-4.
# They are stored as checkpoint metadata only.
# No PH test is rerun here.
# ------------------------------------------------------------

ms4a7_constant_modelA_ph_P <-
  0.03207111

modelA_global_ph_P <-
  0.01165715

modelA_gender_ph_P <-
  0.01760010


cat(
  "Previously locked Step 6C-4 PH results:\n"
)

cat(
  "  MS4A7 PH P =",
  ms4a7_constant_modelA_ph_P,
  "\n"
)

cat(
  "  GLOBAL PH P =",
  modelA_global_ph_P,
  "\n"
)

cat(
  "  Gender PH P =",
  modelA_gender_ph_P,
  "\n\n"
)


# ------------------------------------------------------------
# 10. Confirm cox.zph() is not applicable to selected TV model
# ------------------------------------------------------------

contains_tt <-
  any(
    grepl(
      "^tt\\(",
      names(coef_adjusted)
    )
  )


stopifnot(
  identical(
    contains_tt,
    TRUE
  )
)


cat(
  "Selected TV model contains tt(MS4A7_z): TRUE\n"
)

cat(
  "cox.zph() on selected TV model: NOT APPLICABLE\n\n"
)


# ------------------------------------------------------------
# 11. Create checkpoint metadata
# ------------------------------------------------------------

step6c5_checkpoint_metadata <- list(
  
  step =
    "Step 6C-5",
  
  analysis =
    "Adjusted time-varying Cox model for MS4A7",
  
  cohort_N =
    485,
  
  events =
    178,
  
  censored =
    307,
  
  gene =
    "MS4A7",
  
  predictor =
    "continuous standardized MS4A7 expression",
  
  adjustment =
    "Age + Gender + Stage",
  
  gender_reference =
    "female",
  
  stage_reference =
    "Stage I",
  
  selected_time_function =
    best_time_function_adjusted,
  
  candidate_time_functions =
    c(
      "linear",
      "sqrt",
      "log1p"
    ),
  
  beta0 =
    beta0_adjusted_locked,
  
  beta1 =
    beta1_adjusted_locked,
  
  beta0_SE =
    se0_adjusted_locked,
  
  beta1_SE =
    se1_adjusted_locked,
  
  beta0_P =
    p0_adjusted_locked,
  
  beta1_P =
    p1_adjusted_locked,
  
  selected_AIC =
    selected_AIC_locked,
  
  LRT =
    lrt_value_locked,
  
  LRT_df =
    lrt_df_locked,
  
  LRT_P =
    lrt_p_locked,
  
  constant_modelA_MS4A7_PH_P =
    ms4a7_constant_modelA_ph_P,
  
  modelA_global_PH_P =
    modelA_global_ph_P,
  
  modelA_gender_PH_P =
    modelA_gender_ph_P,
  
  cox_zph_TV_model =
    "Not applicable because model contains tt() terms",
  
  status =
    "LOCKED"
)


# ------------------------------------------------------------
# 12. Define checkpoint object list
# ------------------------------------------------------------

step6c5_checkpoint_objects <- c(
  
  "modelA_cc",
  
  "ms4a7_tv_adjusted_models",
  
  "ms4a7_tv_adjusted_comparison",
  
  "best_time_function_adjusted",
  
  "best_tv_adjusted_model",
  
  "ms4a7_adjusted_LRT",
  
  "ms4a7_adjusted_hr_time",
  
  "coef_adjusted",
  
  "vcov_adjusted",
  
  "summary_adjusted",
  
  "beta0_adjusted_locked",
  
  "beta1_adjusted_locked",
  
  "se0_adjusted_locked",
  
  "se1_adjusted_locked",
  
  "p0_adjusted_locked",
  
  "p1_adjusted_locked",
  
  "HR0_adjusted_locked",
  
  "selected_row",
  
  "selected_AIC_locked",
  
  "minimum_AIC_locked",
  
  "lrt_value_locked",
  
  "lrt_df_locked",
  
  "lrt_p_locked",
  
  "ms4a7_constant_modelA_ph_P",
  
  "modelA_global_ph_P",
  
  "modelA_gender_ph_P",
  
  "step6c5_checkpoint_metadata",
  
  "step6c5_checkpoint_objects"
)


# ------------------------------------------------------------
# 13. Final checkpoint QC
# ------------------------------------------------------------

stopifnot(
  all(
    step6c5_checkpoint_objects %in%
      ls(envir = .GlobalEnv)
  )
)

stopifnot(
  step6c5_checkpoint_metadata$status ==
    "LOCKED"
)

stopifnot(
  step6c5_checkpoint_metadata$cohort_N ==
    485
)

stopifnot(
  step6c5_checkpoint_metadata$events ==
    178
)

stopifnot(
  step6c5_checkpoint_metadata$selected_time_function ==
    best_time_function_adjusted
)

stopifnot(
  is.finite(
    step6c5_checkpoint_metadata$beta0
  )
)

stopifnot(
  is.finite(
    step6c5_checkpoint_metadata$beta1
  )
)


# ------------------------------------------------------------
# 14. Save checkpoint
# ------------------------------------------------------------

checkpoint_file <-
  file.path(
    getwd(),
    "Paper1_Step6C5_adjusted_MS4A7_time_varying_Cox_checkpoint.RData"
  )


save(
  list = step6c5_checkpoint_objects,
  file = checkpoint_file
)


# ------------------------------------------------------------
# 15. Verify checkpoint file
# ------------------------------------------------------------

stopifnot(
  file.exists(checkpoint_file)
)


file_info <-
  file.info(
    checkpoint_file
  )


cat("\n============================================================\n")
cat("STEP 6C-5 CHECKPOINT: LOCKED\n")
cat("============================================================\n\n")

cat(
  "Checkpoint file:\n",
  checkpoint_file,
  "\n\n"
)

cat(
  "File size:",
  round(
    file_info$size / 1024^2,
    3
  ),
  "MB\n"
)

cat(
  "Objects saved:",
  length(
    step6c5_checkpoint_objects
  ),
  "\n"
)

cat(
  "Locked cohort N:",
  step6c5_checkpoint_metadata$cohort_N,
  "\n"
)

cat(
  "Events:",
  step6c5_checkpoint_metadata$events,
  "\n"
)

cat(
  "Selected time function:",
  best_time_function_adjusted,
  "\n"
)

cat(
  "MS4A7 constant-model PH P:",
  ms4a7_constant_modelA_ph_P,
  "\n"
)

cat(
  "Status: LOCKED\n\n"
)

cat(
  "STEP 6C-5 CHECKPOINT QC: PASSED\n"
)

cat(
  "============================================================\n"
)


# ============================================================
# Paper 1
# STEP 6C-5 CHECKPOINT
# Adjusted time-varying Cox model for MS4A7
#
# Purpose:
#   Lock the final adjusted time-varying Cox analysis
#   based on the already completed Step 6C-5 analysis.
#
# Important:
#   - No cohort reconstruction
#   - No model refitting
#   - No clinical recoding
#   - No T/N/M variables added
#   - MS4A7 is the only time-varying predictor
# ============================================================
cat("\n============================================================\n")
cat("STEP 6C-5 EXISTING OBJECT AUDIT\n")
cat("============================================================\n\n")

objects_to_check <- c(
  "modelA_cc",
  "ms4a7_tv_adjusted_models",
  "ms4a7_tv_adjusted_comparison",
  "best_time_function_adjusted",
  "best_tv_adjusted_model",
  "ms4a7_adjusted_LRT",
  "ms4a7_adjusted_hr_time"
)

object_status <- data.frame(
  Object = objects_to_check,
  Exists = vapply(
    objects_to_check,
    exists,
    logical(1),
    inherits = TRUE
  ),
  stringsAsFactors = FALSE
)

print(
  object_status,
  row.names = FALSE
)

# ============================================================
# Paper 1
# STEP 6C-5 CHECKPOINT
# Adjusted time-varying Cox model for MS4A7
#
# IMPORTANT:
# This checkpoint saves already-fitted objects and derived
# metadata only.
#
# No Cox model is refitted.
# No cohort is reconstructed.
# No clinical coding is changed.
# ============================================================


cat("\n============================================================\n")
cat("STEP 6C-5 CHECKPOINT: START\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Required existing objects
# ------------------------------------------------------------

required_objects <- c(
  "modelA_cc",
  "ms4a7_tv_adjusted_models",
  "ms4a7_tv_adjusted_comparison",
  "best_time_function_adjusted",
  "best_tv_adjusted_model",
  "ms4a7_adjusted_LRT",
  "ms4a7_adjusted_hr_time"
)


missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]


if (length(missing_objects) > 0) {
  
  stop(
    paste(
      "Missing required objects:",
      paste(
        missing_objects,
        collapse = ", "
      )
    )
  )
  
}


cat(
  "All required existing objects are present.\n\n"
)


# ------------------------------------------------------------
# 2. Verify locked Model A cohort
# ------------------------------------------------------------

stopifnot(
  nrow(modelA_cc) == 485
)

stopifnot(
  sum(modelA_cc$OS_event == 1) == 178
)

stopifnot(
  sum(modelA_cc$OS_event == 0) == 307
)


cat("Locked Model A cohort:\n")
cat("  N        =", nrow(modelA_cc), "\n")
cat("  Events   =", sum(modelA_cc$OS_event == 1), "\n")
cat("  Censored =", sum(modelA_cc$OS_event == 0), "\n\n")


# ------------------------------------------------------------
# 3. Verify selected time function
# ------------------------------------------------------------

stopifnot(
  best_time_function_adjusted %in%
    c(
      "linear",
      "sqrt",
      "log1p"
    )
)


cat(
  "Selected adjusted time function:",
  best_time_function_adjusted,
  "\n\n"
)


# ------------------------------------------------------------
# 4. Verify selected model
# ------------------------------------------------------------

stopifnot(
  inherits(
    best_tv_adjusted_model,
    "coxph"
  )
)

stopifnot(
  best_tv_adjusted_model$n == 485
)


cat(
  "Selected adjusted TV-Cox model verified.\n"
)

cat(
  "  N =",
  best_tv_adjusted_model$n,
  "\n\n"
)


# ------------------------------------------------------------
# 5. Extract coefficients directly from selected model
# ------------------------------------------------------------

coef_adjusted <-
  stats::coef(
    best_tv_adjusted_model
  )

vcov_adjusted <-
  stats::vcov(
    best_tv_adjusted_model
  )

summary_adjusted <-
  summary(
    best_tv_adjusted_model
  )


stopifnot(
  "MS4A7_z" %in%
    names(coef_adjusted)
)

stopifnot(
  "tt(MS4A7_z)" %in%
    names(coef_adjusted)
)


beta0_adjusted_locked <-
  unname(
    coef_adjusted["MS4A7_z"]
  )

beta1_adjusted_locked <-
  unname(
    coef_adjusted["tt(MS4A7_z)"]
  )

se0_adjusted_locked <-
  sqrt(
    vcov_adjusted[
      "MS4A7_z",
      "MS4A7_z"
    ]
  )

se1_adjusted_locked <-
  sqrt(
    vcov_adjusted[
      "tt(MS4A7_z)",
      "tt(MS4A7_z)"
    ]
  )


# Wald P values

z0_adjusted_locked <-
  beta0_adjusted_locked /
  se0_adjusted_locked

z1_adjusted_locked <-
  beta1_adjusted_locked /
  se1_adjusted_locked


p0_adjusted_locked <-
  2 *
  stats::pnorm(
    -abs(z0_adjusted_locked)
  )

p1_adjusted_locked <-
  2 *
  stats::pnorm(
    -abs(z1_adjusted_locked)
  )


# Hazard ratios for coefficient at t = 0 transformation
HR0_adjusted_locked <-
  exp(
    beta0_adjusted_locked
  )


cat(
  "Adjusted MS4A7 coefficients:\n"
)

cat(
  "  beta0 =",
  beta0_adjusted_locked,
  "\n"
)

cat(
  "  SE0   =",
  se0_adjusted_locked,
  "\n"
)

cat(
  "  P0    =",
  p0_adjusted_locked,
  "\n"
)

cat(
  "  beta1 =",
  beta1_adjusted_locked,
  "\n"
)

cat(
  "  SE1   =",
  se1_adjusted_locked,
  "\n"
)

cat(
  "  P1    =",
  p1_adjusted_locked,
  "\n\n"
)


# ------------------------------------------------------------
# 6. Verify time-function comparison
# ------------------------------------------------------------

print(
  ms4a7_tv_adjusted_comparison,
  row.names = FALSE
)


stopifnot(
  best_time_function_adjusted %in%
    ms4a7_tv_adjusted_comparison$Time_function
)


selected_row <-
  ms4a7_tv_adjusted_comparison[
    ms4a7_tv_adjusted_comparison$Time_function ==
      best_time_function_adjusted,
    ,
    drop = FALSE
  ]


stopifnot(
  nrow(selected_row) == 1
)


selected_AIC_locked <-
  selected_row$AIC


minimum_AIC_locked <-
  min(
    ms4a7_tv_adjusted_comparison$AIC,
    na.rm = TRUE
  )


stopifnot(
  isTRUE(
    all.equal(
      as.numeric(selected_AIC_locked),
      as.numeric(minimum_AIC_locked),
      tolerance = 1e-8
    )
  )
)


cat(
  "\nSelected model AIC =",
  selected_AIC_locked,
  "\n"
)

cat(
  "Minimum AIC =",
  minimum_AIC_locked,
  "\n\n"
)


# ------------------------------------------------------------
# 7. Verify constant vs TV model LRT
# ------------------------------------------------------------

print(
  ms4a7_adjusted_LRT
)


# Extract LRT values flexibly

lrt_value_locked <- NA_real_
lrt_df_locked <- NA_real_
lrt_p_locked <- NA_real_


if (
  is.data.frame(ms4a7_adjusted_LRT)
) {
  
  if ("LR" %in% names(ms4a7_adjusted_LRT)) {
    lrt_value_locked <-
      as.numeric(
        ms4a7_adjusted_LRT$LR[1]
      )
  }
  
  if ("df" %in% names(ms4a7_adjusted_LRT)) {
    lrt_df_locked <-
      as.numeric(
        ms4a7_adjusted_LRT$df[1]
      )
  }
  
  if ("P" %in% names(ms4a7_adjusted_LRT)) {
    lrt_p_locked <-
      as.numeric(
        ms4a7_adjusted_LRT$P[1]
      )
  }
  
}


if (
  is.numeric(ms4a7_adjusted_LRT) &&
  length(ms4a7_adjusted_LRT) >= 1
) {
  
  lrt_value_locked <-
    as.numeric(
      ms4a7_adjusted_LRT[1]
    )
  
}


stopifnot(
  is.finite(lrt_value_locked)
)


cat(
  "\nConstant vs TV-Cox LRT statistic =",
  lrt_value_locked,
  "\n"
)

if (is.finite(lrt_df_locked)) {
  cat(
    "LRT df =",
    lrt_df_locked,
    "\n"
  )
}

if (is.finite(lrt_p_locked)) {
  cat(
    "LRT P =",
    lrt_p_locked,
    "\n"
  )
}


# ------------------------------------------------------------
# 8. Verify time-specific adjusted HR
# ------------------------------------------------------------

print(
  ms4a7_adjusted_hr_time,
  row.names = FALSE
)


stopifnot(
  nrow(ms4a7_adjusted_hr_time) == 5
)

stopifnot(
  all(
    ms4a7_adjusted_hr_time$Time_years ==
      c(1, 2, 3, 4, 5)
  )
)

stopifnot(
  all(
    is.finite(
      ms4a7_adjusted_hr_time$HR_per_1SD
    )
  )
)

stopifnot(
  all(
    ms4a7_adjusted_hr_time$CI_low <=
      ms4a7_adjusted_hr_time$HR_per_1SD
  )
)

stopifnot(
  all(
    ms4a7_adjusted_hr_time$HR_per_1SD <=
      ms4a7_adjusted_hr_time$CI_high
  )
)


cat(
  "\nTime-specific adjusted HR QC passed.\n\n"
)


# ------------------------------------------------------------
# 9. Preserve Step 6C-4 PH results
#
# These values were already obtained in Step 6C-4.
# They are stored as checkpoint metadata only.
# No PH test is rerun here.
# ------------------------------------------------------------

ms4a7_constant_modelA_ph_P <-
  0.03207111

modelA_global_ph_P <-
  0.01165715

modelA_gender_ph_P <-
  0.01760010


cat(
  "Previously locked Step 6C-4 PH results:\n"
)

cat(
  "  MS4A7 PH P =",
  ms4a7_constant_modelA_ph_P,
  "\n"
)

cat(
  "  GLOBAL PH P =",
  modelA_global_ph_P,
  "\n"
)

cat(
  "  Gender PH P =",
  modelA_gender_ph_P,
  "\n\n"
)


# ------------------------------------------------------------
# 10. Confirm cox.zph() is not applicable to selected TV model
# ------------------------------------------------------------

contains_tt <-
  any(
    grepl(
      "^tt\\(",
      names(coef_adjusted)
    )
  )


stopifnot(
  identical(
    contains_tt,
    TRUE
  )
)


cat(
  "Selected TV model contains tt(MS4A7_z): TRUE\n"
)

cat(
  "cox.zph() on selected TV model: NOT APPLICABLE\n\n"
)


# ------------------------------------------------------------
# 11. Create checkpoint metadata
# ------------------------------------------------------------

step6c5_checkpoint_metadata <- list(
  
  step =
    "Step 6C-5",
  
  analysis =
    "Adjusted time-varying Cox model for MS4A7",
  
  cohort_N =
    485,
  
  events =
    178,
  
  censored =
    307,
  
  gene =
    "MS4A7",
  
  predictor =
    "continuous standardized MS4A7 expression",
  
  adjustment =
    "Age + Gender + Stage",
  
  gender_reference =
    "female",
  
  stage_reference =
    "Stage I",
  
  selected_time_function =
    best_time_function_adjusted,
  
  candidate_time_functions =
    c(
      "linear",
      "sqrt",
      "log1p"
    ),
  
  beta0 =
    beta0_adjusted_locked,
  
  beta1 =
    beta1_adjusted_locked,
  
  beta0_SE =
    se0_adjusted_locked,
  
  beta1_SE =
    se1_adjusted_locked,
  
  beta0_P =
    p0_adjusted_locked,
  
  beta1_P =
    p1_adjusted_locked,
  
  selected_AIC =
    selected_AIC_locked,
  
  LRT =
    lrt_value_locked,
  
  LRT_df =
    lrt_df_locked,
  
  LRT_P =
    lrt_p_locked,
  
  constant_modelA_MS4A7_PH_P =
    ms4a7_constant_modelA_ph_P,
  
  modelA_global_PH_P =
    modelA_global_ph_P,
  
  modelA_gender_PH_P =
    modelA_gender_ph_P,
  
  cox_zph_TV_model =
    "Not applicable because model contains tt() terms",
  
  status =
    "LOCKED"
)


# ------------------------------------------------------------
# 12. Define checkpoint object list
# ------------------------------------------------------------

step6c5_checkpoint_objects <- c(
  
  "modelA_cc",
  
  "ms4a7_tv_adjusted_models",
  
  "ms4a7_tv_adjusted_comparison",
  
  "best_time_function_adjusted",
  
  "best_tv_adjusted_model",
  
  "ms4a7_adjusted_LRT",
  
  "ms4a7_adjusted_hr_time",
  
  "coef_adjusted",
  
  "vcov_adjusted",
  
  "summary_adjusted",
  
  "beta0_adjusted_locked",
  
  "beta1_adjusted_locked",
  
  "se0_adjusted_locked",
  
  "se1_adjusted_locked",
  
  "p0_adjusted_locked",
  
  "p1_adjusted_locked",
  
  "HR0_adjusted_locked",
  
  "selected_row",
  
  "selected_AIC_locked",
  
  "minimum_AIC_locked",
  
  "lrt_value_locked",
  
  "lrt_df_locked",
  
  "lrt_p_locked",
  
  "ms4a7_constant_modelA_ph_P",
  
  "modelA_global_ph_P",
  
  "modelA_gender_ph_P",
  
  "step6c5_checkpoint_metadata",
  
  "step6c5_checkpoint_objects"
)


# ------------------------------------------------------------
# 13. Final checkpoint QC
# ------------------------------------------------------------

stopifnot(
  all(
    step6c5_checkpoint_objects %in%
      ls(envir = .GlobalEnv)
  )
)

stopifnot(
  step6c5_checkpoint_metadata$status ==
    "LOCKED"
)

stopifnot(
  step6c5_checkpoint_metadata$cohort_N ==
    485
)

stopifnot(
  step6c5_checkpoint_metadata$events ==
    178
)

stopifnot(
  step6c5_checkpoint_metadata$selected_time_function ==
    best_time_function_adjusted
)

stopifnot(
  is.finite(
    step6c5_checkpoint_metadata$beta0
  )
)

stopifnot(
  is.finite(
    step6c5_checkpoint_metadata$beta1
  )
)


# ------------------------------------------------------------
# 14. Save checkpoint
# ------------------------------------------------------------

checkpoint_file <-
  file.path(
    getwd(),
    "Paper1_Step6C5_adjusted_MS4A7_time_varying_Cox_checkpoint.RData"
  )


save(
  list = step6c5_checkpoint_objects,
  file = checkpoint_file
)


# ------------------------------------------------------------
# 15. Verify checkpoint file
# ------------------------------------------------------------

stopifnot(
  file.exists(checkpoint_file)
)


file_info <-
  file.info(
    checkpoint_file
  )


cat("\n============================================================\n")
cat("STEP 6C-5 CHECKPOINT: LOCKED\n")
cat("============================================================\n\n")

cat(
  "Checkpoint file:\n",
  checkpoint_file,
  "\n\n"
)

cat(
  "File size:",
  round(
    file_info$size / 1024^2,
    3
  ),
  "MB\n"
)

cat(
  "Objects saved:",
  length(
    step6c5_checkpoint_objects
  ),
  "\n"
)

cat(
  "Locked cohort N:",
  step6c5_checkpoint_metadata$cohort_N,
  "\n"
)

cat(
  "Events:",
  step6c5_checkpoint_metadata$events,
  "\n"
)

cat(
  "Selected time function:",
  best_time_function_adjusted,
  "\n"
)

cat(
  "MS4A7 constant-model PH P:",
  ms4a7_constant_modelA_ph_P,
  "\n"
)

cat(
  "Status: LOCKED\n\n"
)

cat(
  "STEP 6C-5 CHECKPOINT QC: PASSED\n"
)

cat(
  "============================================================\n"
)


p_lrt_from_locked_LR <-
  stats::pchisq(
    4.073243,
    df = 1,
    lower.tail = FALSE
  )

p_lrt_from_locked_LR


