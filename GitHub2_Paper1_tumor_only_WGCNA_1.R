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


############################
# 0. Set working directory
############################
setwd("input_dir/GitHub2_Paper1_tumor-only WGCNA")
getwd()


############################################################
# Workflow-04-1
#
# Reconstruction of the tumor-only WGCNA cohort
#
# Purpose:
#   1. Read the original TCGA expression matrix
#   2. Identify Primary Tumor samples (-01)
#   3. Match Primary Tumor samples to TCGA clinical data
#   4. Quantify clinical matching and missingness
#   5. Determine the final reliable tumor cohort
#
# IMPORTANT:
#   - Expression source = Merge_RNA_seq_FPKM_UQ.txt
#   - Do NOT use expr_symbol_counts.csv
#   - No WGCNA is performed in this workflow
#   - Original data are not modified
############################################################


############################################################
# 0. Paths
############################################################

project_dir <- input_dir

expression_file <- file.path(
  project_dir,
  "Merge_RNA_seq_FPKM_UQ.txt"
)

############################################################
# IMPORTANT:
# Replace the following filename with the actual
# TCGA clinical Excel/CSV file used for the original
# clinical matching.
############################################################

clinical_file <- file.path(
  project_dir,
  "YOUR_TCGA_CLINICAL_FILE.xlsx"
)


############################################################
# 1. Read expression matrix
############################################################

expr <- read.table(
  expression_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


############################################################
# 2. Basic expression structure
############################################################

cat("\n========================================\n")
cat("ORIGINAL EXPRESSION MATRIX\n")
cat("========================================\n\n")

cat(
  "Genes: ",
  nrow(expr),
  "\n",
  sep = ""
)

cat(
  "Samples: ",
  ncol(expr) - 1,
  "\n",
  sep = ""
)


############################################################
# 3. Extract sample IDs
############################################################

sample_ids <- colnames(expr)[-1]


############################################################
# 4. Identify TCGA sample types
#
# TCGA barcode:
#   -01 = Primary Tumor
#   -02 = Recurrent/Metastatic
#   -11 = Solid Tissue Normal
############################################################

sample_type_code <- substr(
  sample_ids,
  14,
  15
)


############################################################
# 5. Summarize original sample composition
############################################################

sample_type_table <- table(
  sample_type_code,
  useNA = "ifany"
)

cat("\n========================================\n")
cat("ORIGINAL SAMPLE COMPOSITION\n")
cat("========================================\n\n")

print(sample_type_table)


############################################################
# 6. Extract Primary Tumor samples only
############################################################

tumor_sample_ids <- sample_ids[
  sample_type_code == "01"
]

cat("\n========================================\n")
cat("PRIMARY TUMOR SAMPLES\n")
cat("========================================\n\n")

cat(
  "Primary Tumor samples: ",
  length(tumor_sample_ids),
  "\n",
  sep = ""
)


############################################################
# 7. Convert sample barcode to patient barcode
#
# TCGA-XX-YYYY-01
#       ↓
# TCGA-XX-YYYY
############################################################

tumor_patient_ids <- substr(
  tumor_sample_ids,
  1,
  12
)


############################################################
# 8. Check whether multiple Primary Tumor samples
#    exist for the same patient
############################################################

patient_frequency <- table(
  tumor_patient_ids
)

multiple_tumor_patients <- names(
  patient_frequency[
    patient_frequency > 1
  ]
)

cat("\n========================================\n")
cat("PATIENT-LEVEL DUPLICATION CHECK\n")
cat("========================================\n\n")

cat(
  "Unique Primary Tumor patients: ",
  length(
    unique(tumor_patient_ids)
  ),
  "\n",
  sep = ""
)

cat(
  "Patients with >1 Primary Tumor sample: ",
  length(
    multiple_tumor_patients
  ),
  "\n",
  sep = ""
)

if (
  length(multiple_tumor_patients) > 0
) {
  
  cat("\nPatients with multiple Primary Tumor samples:\n")
  
  print(
    patient_frequency[
      patient_frequency > 1
    ]
  )
}


############################################################
# 9. Read clinical data
############################################################

library(readxl)

clinical <- read_excel(
  clinical_file
)

clinical <- as.data.frame(
  clinical,
  stringsAsFactors = FALSE
)


############################################################
# 10. Inspect clinical structure
############################################################

cat("\n========================================\n")
cat("CLINICAL DATA STRUCTURE\n")
cat("========================================\n\n")

cat(
  "Clinical rows: ",
  nrow(clinical),
  "\n",
  sep = ""
)

cat(
  "Clinical columns: ",
  ncol(clinical),
  "\n",
  sep = ""
)

cat("\nClinical column names:\n")

print(
  colnames(clinical)
)


############################################################
# 11. IMPORTANT:
#     Identify the TCGA sample/patient ID column
############################################################

cat("\n========================================\n")
cat("POSSIBLE ID COLUMNS\n")
cat("========================================\n\n")

id_candidates <- colnames(clinical)[
  grepl(
    "sample|barcode|patient|submitter",
    colnames(clinical),
    ignore.case = TRUE
  )
]

print(
  id_candidates
) 

############################################################
# Workflow-04-1
# Part B: Exact sample-level clinical matching
############################################################


############################################################
# 12. Check whether clinical Sample IDs are unique
############################################################

cat("\n========================================\n")
cat("CLINICAL SAMPLE ID QUALITY CHECK\n")
cat("========================================\n\n")

cat(
  "Clinical rows: ",
  nrow(clinical),
  "\n",
  sep = ""
)

cat(
  "Unique clinical Sample IDs: ",
  length(
    unique(clinical$Sample)
  ),
  "\n",
  sep = ""
)

cat(
  "Duplicated clinical Sample IDs: ",
  sum(
    duplicated(clinical$Sample)
  ),
  "\n",
  sep = ""
)


############################################################
# 13. Check missing Sample IDs
############################################################

cat(
  "Missing clinical Sample IDs: ",
  sum(
    is.na(clinical$Sample) |
      clinical$Sample == ""
  ),
  "\n",
  sep = ""
)


############################################################
# 14. Exact matching of expression samples
#     to clinical Sample IDs
############################################################

expression_in_clinical <- sample_ids %in%
  clinical$Sample

cat("\n========================================\n")
cat("EXPRESSION → CLINICAL MATCHING\n")
cat("========================================\n\n")

cat(
  "Expression samples: ",
  length(sample_ids),
  "\n",
  sep = ""
)

cat(
  "Expression samples with clinical match: ",
  sum(expression_in_clinical),
  "\n",
  sep = ""
)

cat(
  "Expression samples without clinical match: ",
  sum(!expression_in_clinical),
  "\n",
  sep = ""
)


############################################################
# 15. List unmatched expression samples
############################################################

unmatched_expression <- sample_ids[
  !expression_in_clinical
]

if (
  length(unmatched_expression) > 0
) {
  
  cat("\nUnmatched expression samples:\n")
  
  print(
    unmatched_expression
  )
}


############################################################
# 16. Check clinical samples not present
#     in expression matrix
############################################################

clinical_in_expression <- clinical$Sample %in%
  sample_ids

cat("\n========================================\n")
cat("CLINICAL → EXPRESSION MATCHING\n")
cat("========================================\n\n")

cat(
  "Clinical samples: ",
  nrow(clinical),
  "\n",
  sep = ""
)

cat(
  "Clinical samples with expression match: ",
  sum(clinical_in_expression),
  "\n",
  sep = ""
)

cat(
  "Clinical samples without expression match: ",
  sum(!clinical_in_expression),
  "\n",
  sep = ""
)


############################################################
# 17. Build exact matched clinical table
############################################################

clinical_matched_sample <- clinical[
  match(
    sample_ids,
    clinical$Sample
  ),
]

rownames(
  clinical_matched_sample
) <- sample_ids


############################################################
# 18. Verify one-to-one matching
############################################################

cat("\n========================================\n")
cat("MATCHING CONSISTENCY CHECK\n")
cat("========================================\n\n")

cat(
  "Rows in matched clinical table: ",
  nrow(clinical_matched_sample),
  "\n",
  sep = ""
)

cat(
  "Rows expected: ",
  length(sample_ids),
  "\n",
  sep = ""
)

cat(
  "Number of NA Sample IDs after matching: ",
  sum(
    is.na(
      clinical_matched_sample$Sample
    )
  ),
  "\n",
  sep = ""
)


############################################################
# 19. Compare tumor classification
############################################################

cat("\n========================================\n")
cat("TUMOR/NORMAL CLASSIFICATION\n")
cat("========================================\n\n")

print(
  table(
    clinical_matched_sample$tumor_normal,
    useNA = "ifany"
  )
)


############################################################
# 20. Primary Tumor only
############################################################

tumor_clinical <- clinical_matched_sample[
  sample_type_code == "01",
  ,
  drop = FALSE
]

cat("\n========================================\n")
cat("PRIMARY TUMOR CLINICAL COHORT\n")
cat("========================================\n\n")

cat(
  "Primary Tumor samples: ",
  nrow(tumor_clinical),
  "\n",
  sep = ""
)

cat(
  "Primary Tumor classified as tumor: ",
  sum(
    tumor_clinical$tumor_normal == "tumor",
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Primary Tumor classified as normal: ",
  sum(
    tumor_clinical$tumor_normal == "normal",
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Primary Tumor with NA tumor_normal: ",
  sum(
    is.na(
      tumor_clinical$tumor_normal
    )
  ),
  "\n",
  sep = ""
) 

############################################################
# Workflow-04-1 Part C
# Primary Tumor clinical completeness audit
############################################################

cat("\n========================================\n")
cat("PRIMARY TUMOR CLINICAL COMPLETENESS AUDIT\n")
cat("========================================\n\n")

############################################################
# 21. Define Primary Tumor cohort
############################################################

primary_tumor_ids <- sample_ids[
  sample_type_code == "01"
]

primary_tumor_clinical <- clinical_matched_sample[
  primary_tumor_ids,
  ,
  drop = FALSE
]

cat(
  "Primary Tumor samples: ",
  nrow(primary_tumor_clinical),
  "\n",
  sep = ""
)

############################################################
# 22. Check NA / empty values for all clinical variables
############################################################

clinical_vars <- c(
  "tumor_normal",
  "patient_id",
  "Gender",
  "Age",
  "OS.time",
  "OS.status",
  "pathologic_stage_2",
  "pathologic_t",
  "pathologic_n",
  "pathologic_m",
  "pathologic_stage"
)

clinical_missing_summary <- data.frame(
  Variable = clinical_vars,
  N = sapply(
    primary_tumor_clinical[, clinical_vars, drop = FALSE],
    function(x) {
      sum(
        is.na(x) |
          (is.character(x) & trimws(x) == "")
      )
    }
  ),
  Percent = sapply(
    primary_tumor_clinical[, clinical_vars, drop = FALSE],
    function(x) {
      100 * sum(
        is.na(x) |
          (is.character(x) & trimws(x) == "")
      ) / length(x)
    }
  )
)

print(
  clinical_missing_summary,
  row.names = FALSE
)

############################################################
# 23. Detailed distribution of major clinical variables
############################################################

cat("\n========================================\n")
cat("GENDER DISTRIBUTION\n")
cat("========================================\n\n")

print(
  table(
    primary_tumor_clinical$Gender,
    useNA = "ifany"
  )
)

cat("\n========================================\n")
cat("OS.STATUS DISTRIBUTION\n")
cat("========================================\n\n")

print(
  table(
    primary_tumor_clinical$OS.status,
    useNA = "ifany"
  )
)

cat("\n========================================\n")
cat("PATHOLOGIC STAGE DISTRIBUTION\n")
cat("========================================\n\n")

print(
  table(
    primary_tumor_clinical$pathologic_stage,
    useNA = "ifany"
  )
)

cat("\n========================================\n")
cat("PATHOLOGIC STAGE 2 DISTRIBUTION\n")
cat("========================================\n\n")

print(
  table(
    primary_tumor_clinical$pathologic_stage_2,
    useNA = "ifany"
  )
)

############################################################
# 24. T / N / M distributions
############################################################

cat("\n========================================\n")
cat("PATHOLOGIC T DISTRIBUTION\n")
cat("========================================\n\n")

print(
  table(
    primary_tumor_clinical$pathologic_t,
    useNA = "ifany"
  )
)

cat("\n========================================\n")
cat("PATHOLOGIC N DISTRIBUTION\n")
cat("========================================\n\n")

print(
  table(
    primary_tumor_clinical$pathologic_n,
    useNA = "ifany"
  )
)

cat("\n========================================\n")
cat("PATHOLOGIC M DISTRIBUTION\n")
cat("========================================\n\n")

print(
  table(
    primary_tumor_clinical$pathologic_m,
    useNA = "ifany"
  )
)

############################################################
# 25. OS.time numerical audit
############################################################

cat("\n========================================\n")
cat("OS.TIME AUDIT\n")
cat("========================================\n\n")

cat(
  "OS.time class: ",
  class(primary_tumor_clinical$OS.time),
  "\n",
  sep = ""
)

cat(
  "OS.time NA: ",
  sum(is.na(primary_tumor_clinical$OS.time)),
  "\n",
  sep = ""
)

if (is.numeric(primary_tumor_clinical$OS.time)) {
  
  cat(
    "OS.time minimum: ",
    min(
      primary_tumor_clinical$OS.time,
      na.rm = TRUE
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "OS.time median: ",
    median(
      primary_tumor_clinical$OS.time,
      na.rm = TRUE
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "OS.time maximum: ",
    max(
      primary_tumor_clinical$OS.time,
      na.rm = TRUE
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "OS.time <= 0: ",
    sum(
      primary_tumor_clinical$OS.time <= 0,
      na.rm = TRUE
    ),
    "\n",
    sep = ""
  )
}

############################################################
# 26. OS.status coding
############################################################

cat("\n========================================\n")
cat("OS.STATUS CODING\n")
cat("========================================\n\n")

print(
  unique(
    primary_tumor_clinical$OS.status
  )
)

############################################################
# 27. Save clinical completeness summary
############################################################

write.table(
  clinical_missing_summary,
  file =
    "result_dir/Workflow-04-1C_clinical_missing_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nSaved:\n",
  "Workflow-04-1C_clinical_missing_summary.txt\n"
) 


############################################################
# Workflow-04-1 Part D
# Final clinical missing-value encoding audit
############################################################

cat("\n========================================\n")
cat("FINAL CLINICAL MISSING-VALUE AUDIT\n")
cat("========================================\n\n")

############################################################
# 28. Variable classes
############################################################

audit_vars <- c(
  "Age",
  "OS.time",
  "OS.status",
  "pathologic_stage_2",
  "pathologic_t",
  "pathologic_n",
  "pathologic_m",
  "pathologic_stage"
)

cat("VARIABLE CLASSES\n\n")

print(
  sapply(
    primary_tumor_clinical[, audit_vars, drop = FALSE],
    class
  )
)

############################################################
# 29. Detailed missing-value encoding
############################################################

audit_missing <- function(x) {
  
  x_chr <- trimws(
    as.character(x)
  )
  
  c(
    True_NA = sum(is.na(x)),
    
    Empty = sum(
      !is.na(x) &
        x_chr == ""
    ),
    
    Literal_NA = sum(
      !is.na(x) &
        toupper(x_chr) == "NA"
    ),
    
    N_A = sum(
      !is.na(x) &
        toupper(x_chr) == "N/A"
    ),
    
    NULL = sum(
      !is.na(x) &
        toupper(x_chr) == "NULL"
    ),
    
    Dash = sum(
      !is.na(x) &
        x_chr == "-"
    ),
    
    Dot = sum(
      !is.na(x) &
        x_chr == "."
    )
  )
}

missing_audit <- do.call(
  rbind,
  lapply(
    audit_vars,
    function(v) {
      audit_missing(
        primary_tumor_clinical[[v]]
      )
    }
  )
)

missing_audit <- as.data.frame(
  missing_audit
)

missing_audit$Variable <- rownames(
  missing_audit
)

rownames(
  missing_audit
) <- NULL

missing_audit <- missing_audit[
  ,
  c(
    "Variable",
    "True_NA",
    "Empty",
    "Literal_NA",
    "N_A",
    "NULL",
    "Dash",
    "Dot"
  )
]

cat("\nMISSING-VALUE ENCODING SUMMARY\n\n")

print(
  missing_audit,
  row.names = FALSE
)

############################################################
# 30. Inspect actual problematic values
############################################################

problem_values <- function(x) {
  
  x_chr <- trimws(
    as.character(x)
  )
  
  idx <- is.na(x) |
    x_chr == "" |
    toupper(x_chr) %in%
    c(
      "NA",
      "N/A",
      "NULL",
      "-",
      "."
    )
  
  unique(
    x_chr[idx]
  )
}

for (v in c(
  "Age",
  "OS.time",
  "pathologic_stage_2"
)) {
  
  cat("\n========================================\n")
  cat(v, "\n")
  cat("========================================\n\n")
  
  print(
    problem_values(
      primary_tumor_clinical[[v]]
    )
  )
}

############################################################
# 31. Identify samples with missing Age
############################################################

age_chr <- trimws(
  as.character(
    primary_tumor_clinical$Age
  )
)

age_missing <- is.na(age_chr) |
  age_chr == "" |
  toupper(age_chr) %in%
  c(
    "NA",
    "N/A",
    "NULL",
    "-",
    "."
  )

cat("\n========================================\n")
cat("SAMPLES WITH MISSING AGE\n")
cat("========================================\n\n")

cat(
  "Number: ",
  sum(age_missing),
  "\n\n",
  sep = ""
)

print(
  primary_tumor_clinical$Sample[
    age_missing
  ]
)

############################################################
# 32. Identify samples with missing OS.time
############################################################

os_chr <- trimws(
  as.character(
    primary_tumor_clinical$OS.time
  )
)

os_missing <- is.na(os_chr) |
  os_chr == "" |
  toupper(os_chr) %in%
  c(
    "NA",
    "N/A",
    "NULL",
    "-",
    "."
  )

cat("\n========================================\n")
cat("SAMPLES WITH MISSING OS.TIME\n")
cat("========================================\n\n")

cat(
  "Number: ",
  sum(os_missing),
  "\n\n",
  sep = ""
)

print(
  primary_tumor_clinical$Sample[
    os_missing
  ]
)

############################################################
# 33. Identify samples with missing Stage 2
############################################################

stage2_chr <- trimws(
  as.character(
    primary_tumor_clinical$pathologic_stage_2
  )
)

stage2_missing <- is.na(stage2_chr) |
  stage2_chr == "" |
  toupper(stage2_chr) %in%
  c(
    "NA",
    "N/A",
    "NULL",
    "-",
    "."
  )

cat("\n========================================\n")
cat("SAMPLES WITH MISSING PATHOLOGIC_STAGE_2\n")
cat("========================================\n\n")

cat(
  "Number: ",
  sum(stage2_missing),
  "\n\n",
  sep = ""
)

print(
  primary_tumor_clinical$Sample[
    stage2_missing
  ]
) 


############################################################
# Workflow-04-1 Part E
# Clinical data cleaning and consistency check
#
# IMPORTANT:
# Do NOT delete samples here.
# Do NOT overwrite primary_tumor_clinical.
############################################################

cat("\n========================================\n")
cat("WORKFLOW-04-1 PART E\n")
cat("CLINICAL DATA CLEANING AND QC\n")
cat("========================================\n\n")


############################################################
# 34. Create a clean copy
############################################################

clinical_clean <- primary_tumor_clinical


############################################################
# 35. Convert all recognized missing-value encodings to NA
############################################################

clean_missing <- function(x) {
  
  x_chr <- trimws(
    as.character(x)
  )
  
  x_chr[
    is.na(x_chr) |
      x_chr == "" |
      toupper(x_chr) %in%
      c(
        "NA",
        "N/A",
        "NULL",
        "-",
        "."
      )
  ] <- NA_character_
  
  x_chr
}


clinical_clean[] <- lapply(
  clinical_clean,
  clean_missing
)


############################################################
# 36. Check missing values after cleaning
############################################################

cat("\nMISSING VALUES AFTER STANDARDIZATION\n\n")

missing_clean <- sapply(
  clinical_clean,
  function(x) sum(is.na(x))
)

print(
  missing_clean
)


############################################################
# 37. Convert Age to numeric
############################################################

clinical_clean$Age <- suppressWarnings(
  as.numeric(
    clinical_clean$Age
  )
)


############################################################
# 38. Convert OS.time to numeric
############################################################

clinical_clean$OS.time <- suppressWarnings(
  as.numeric(
    clinical_clean$OS.time
  )
)


############################################################
# 39. Standardize OS.status
############################################################

clinical_clean$OS.status <- trimws(
  clinical_clean$OS.status
)

clinical_clean$OS.status[
  clinical_clean$OS.status == ""
] <- NA

clinical_clean$OS.status <- factor(
  clinical_clean$OS.status,
  levels = c(
    "Alive",
    "Dead"
  )
)


############################################################
# 40. Inspect resulting variable classes
############################################################

cat("\n========================================\n")
cat("VARIABLE CLASSES AFTER CLEANING\n")
cat("========================================\n\n")

print(
  sapply(
    clinical_clean,
    class
  )
)


############################################################
# 41. Check whether OS.time contains non-numeric values
############################################################

cat("\n========================================\n")
cat("OS.TIME QC\n")
cat("========================================\n\n")

cat(
  "OS.time class: ",
  class(clinical_clean$OS.time),
  "\n",
  sep = ""
)

cat(
  "Number of NA: ",
  sum(is.na(clinical_clean$OS.time)),
  "\n",
  sep = ""
)

cat(
  "Minimum: ",
  min(
    clinical_clean$OS.time,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Maximum: ",
  max(
    clinical_clean$OS.time,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Median: ",
  median(
    clinical_clean$OS.time,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


############################################################
# 42. Check Age
############################################################

cat("\n========================================\n")
cat("AGE QC\n")
cat("========================================\n\n")

cat(
  "Number of NA: ",
  sum(is.na(clinical_clean$Age)),
  "\n",
  sep = ""
)

cat(
  "Minimum: ",
  min(
    clinical_clean$Age,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Maximum: ",
  max(
    clinical_clean$Age,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Median: ",
  median(
    clinical_clean$Age,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


############################################################
# 43. OS.status distribution
############################################################

cat("\n========================================\n")
cat("OS.STATUS DISTRIBUTION\n")
cat("========================================\n\n")

print(
  table(
    clinical_clean$OS.status,
    useNA = "always"
  )
)


############################################################
# 44. Stage distribution
############################################################

cat("\n========================================\n")
cat("PATHOLOGIC STAGE DISTRIBUTION\n")
cat("========================================\n\n")

print(
  table(
    clinical_clean$pathologic_stage,
    useNA = "always"
  )
)

cat("\nDetailed stage:\n\n")

print(
  table(
    clinical_clean$pathologic_stage_2,
    useNA = "always"
  )
)


############################################################
# 45. T / N / M distributions
############################################################

cat("\n========================================\n")
cat("T / N / M DISTRIBUTION\n")
cat("========================================\n\n")

cat("T:\n")
print(
  table(
    clinical_clean$pathologic_t,
    useNA = "always"
  )
)

cat("\nN:\n")
print(
  table(
    clinical_clean$pathologic_n,
    useNA = "always"
  )
)

cat("\nM:\n")
print(
  table(
    clinical_clean$pathologic_m,
    useNA = "always"
  )
)


############################################################
# 46. Check duplicate sample IDs
############################################################

cat("\n========================================\n")
cat("SAMPLE ID QC\n")
cat("========================================\n\n")

cat(
  "Total rows: ",
  nrow(clinical_clean),
  "\n",
  sep = ""
)

cat(
  "Unique Sample IDs: ",
  length(
    unique(
      clinical_clean$Sample
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Duplicated Sample IDs: ",
  sum(
    duplicated(
      clinical_clean$Sample
    )
  ),
  "\n",
  sep = ""
)


############################################################
# 47. Check primary tumor samples
############################################################

cat("\n========================================\n")
cat("PRIMARY TUMOR COHORT QC\n")
cat("========================================\n\n")

cat(
  "Primary Tumor samples: ",
  sum(
    substr(
      clinical_clean$Sample,
      14,
      15
    ) == "01"
  ),
  "\n",
  sep = ""
)

cat(
  "Tumor classification among Primary Tumor:\n"
)

print(
  table(
    clinical_clean$tumor_normal[
      substr(
        clinical_clean$Sample,
        14,
        15
      ) == "01"
    ],
    useNA = "always"
  )
)


############################################################
# 48. Save cleaned clinical table
############################################################

write.table(
  clinical_clean,
  file = "result_dir/Workflow-04-1E_clinical_clean.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nClean clinical table saved as:\n")
cat("Workflow-04-1E_clinical_clean.txt\n")


############################################################
# Workflow-04-1F-1
# OS consistency QC
############################################################

cat("\n========================================\n")
cat("OS CONSISTENCY QC\n")
cat("========================================\n\n")


############################################################
# 49. Samples with OS.time = 0
############################################################

os_zero <- !is.na(clinical_clean$OS.time) &
  clinical_clean$OS.time == 0

cat("OS.time = 0 samples: ",
    sum(os_zero),
    "\n\n",
    sep = "")

if (sum(os_zero) > 0) {
  
  print(
    clinical_clean[
      os_zero,
      c(
        "Sample",
        "patient_id",
        "OS.time",
        "OS.status",
        "Age",
        "pathologic_stage_2",
        "pathologic_stage"
      )
    ]
  )
  
}


############################################################
# 50. Sample with missing OS.status
############################################################

os_status_missing <- is.na(
  clinical_clean$OS.status
)

cat("\n========================================\n")
cat("MISSING OS.STATUS\n")
cat("========================================\n\n")

cat(
  "Number: ",
  sum(os_status_missing),
  "\n\n",
  sep = ""
)

if (sum(os_status_missing) > 0) {
  
  print(
    clinical_clean[
      os_status_missing,
      c(
        "Sample",
        "patient_id",
        "OS.time",
        "OS.status",
        "Age",
        "pathologic_stage_2",
        "pathologic_stage"
      )
    ]
  )
  
}


############################################################
# 51. Basic survival-data completeness
############################################################

survival_complete <-
  !is.na(clinical_clean$OS.time) &
  !is.na(clinical_clean$OS.status)

cat("\n========================================\n")
cat("SURVIVAL DATA COMPLETENESS\n")
cat("========================================\n\n")

cat(
  "Complete OS.time + OS.status: ",
  sum(survival_complete),
  "\n",
  sep = ""
)

cat(
  "Incomplete survival data: ",
  sum(!survival_complete),
  "\n",
  sep = ""
)


############################################################
# 52. Event distribution among complete survival cases
############################################################

cat("\n========================================\n")
cat("SURVIVAL EVENT DISTRIBUTION\n")
cat("========================================\n\n")

print(
  table(
    clinical_clean$OS.status[
      survival_complete
    ],
    useNA = "always"
  )
)

###Workflow-04-2A
## Extract 513 Primary Tumor samples and verify bidirectional sample matching

############################################################
# Workflow-04-2A
#
# Purpose:
#   1. Read Merge_RNA_seq_FPKM_UQ.txt
#   2. Identify TCGA Primary Tumor samples (-01)
#   3. Extract exactly the 513 Primary Tumor samples
#   4. Match them against clinical_clean$Sample
#   5. Perform bidirectional consistency checks
#   6. Save the tumor-only expression matrix
#
# IMPORTANT:
#   - No gene filtering
#   - No Ensembl -> Symbol conversion
#   - No log2 transformation
#   - No normalization
#   - No WGCNA
############################################################


############################################################
# 1. Define input file
############################################################

expr_file <- "input_dir /Merge_RNA_seq_FPKM_UQ.txt"


############################################################
# 2. Check input file
############################################################

cat("\n========================================\n")
cat("WORKFLOW-04-2A\n")
cat("PRIMARY TUMOR EXPRESSION EXTRACTION\n")
cat("========================================\n\n")

cat(
  "Input file:\n",
  expr_file,
  "\n\n",
  sep = ""
)

if (!file.exists(expr_file)) {
  stop(
    "ERROR: Input file not found: ",
    expr_file
  )
}

cat("Input file exists: YES\n")


############################################################
# 3. Read expression matrix
############################################################

cat("\n========================================\n")
cat("READING EXPRESSION MATRIX\n")
cat("========================================\n\n")

expr_raw <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = ""
)

cat(
  "Rows: ",
  nrow(expr_raw),
  "\n",
  sep = ""
)

cat(
  "Columns: ",
  ncol(expr_raw),
  "\n",
  sep = ""
)

cat(
  "Expected:\n",
  "60,483 rows × 575 columns\n",
  sep = ""
)


############################################################
# 4. Basic structure check
############################################################

cat("\n========================================\n")
cat("BASIC MATRIX STRUCTURE\n")
cat("========================================\n\n")

cat(
  "First column name: ",
  colnames(expr_raw)[1],
  "\n",
  sep = ""
)

cat(
  "Number of gene rows: ",
  nrow(expr_raw),
  "\n",
  sep = ""
)

cat(
  "Number of sample columns: ",
  ncol(expr_raw) - 1,
  "\n",
  sep = ""
)

cat("\nFirst 10 sample IDs:\n")

print(
  head(
    colnames(expr_raw)[-1],
    10
  )
)


############################################################
# 5. Extract sample IDs
############################################################

expr_sample_ids <- colnames(
  expr_raw
)[-1]


############################################################
# 6. Verify sample IDs are TCGA 15-character barcodes
############################################################

tcga_barcode_pattern <- "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}$"

barcode_valid <- grepl(
  tcga_barcode_pattern,
  expr_sample_ids
)

cat("\n========================================\n")
cat("TCGA SAMPLE ID FORMAT CHECK\n")
cat("========================================\n\n")

cat(
  "Total sample IDs: ",
  length(expr_sample_ids),
  "\n",
  sep = ""
)

cat(
  "Valid 15-character TCGA sample IDs: ",
  sum(barcode_valid),
  "\n",
  sep = ""
)

cat(
  "Invalid sample IDs: ",
  sum(!barcode_valid),
  "\n",
  sep = ""
)

if (any(!barcode_valid)) {
  
  cat("\nInvalid sample IDs:\n")
  
  print(
    expr_sample_ids[
      !barcode_valid
    ]
  )
  
  stop(
    "ERROR: Invalid TCGA sample IDs detected."
  )
}


############################################################
# 7. Determine TCGA sample type
#
# Barcode structure:
# TCGA-XX-XXXX-01
#              ^^
#              sample type code
############################################################

sample_type_code <- substr(
  expr_sample_ids,
  14,
  15
)


############################################################
# 8. Check complete sample-type composition
############################################################

cat("\n========================================\n")
cat("TCGA SAMPLE-TYPE COMPOSITION\n")
cat("========================================\n\n")

sample_type_table <- table(
  sample_type_code,
  useNA = "always"
)

print(
  sample_type_table
)


############################################################
# 9. Explicitly identify Primary Tumor samples
############################################################

primary_tumor_ids <- expr_sample_ids[
  sample_type_code == "01"
]


############################################################
# 10. Primary Tumor cohort size
############################################################

cat("\n========================================\n")
cat("PRIMARY TUMOR COHORT\n")
cat("========================================\n\n")

cat(
  "Primary Tumor samples: ",
  length(primary_tumor_ids),
  "\n",
  sep = ""
)

cat(
  "Unique Primary Tumor samples: ",
  length(unique(primary_tumor_ids)),
  "\n",
  sep = ""
)

if (length(primary_tumor_ids) != 513) {
  
  stop(
    "ERROR: Expected 513 Primary Tumor samples, but found ",
    length(primary_tumor_ids),
    "."
  )
}

if (
  length(unique(primary_tumor_ids)) !=
  length(primary_tumor_ids)
) {
  
  stop(
    "ERROR: Duplicate Primary Tumor sample IDs detected."
  )
}


############################################################
# 11. Extract the 513 Primary Tumor expression columns
############################################################

primary_tumor_column_indices <- which(
  colnames(expr_raw) %in%
    primary_tumor_ids
)


############################################################
# 12. Check extraction indices
############################################################

cat("\n========================================\n")
cat("EXPRESSION COLUMN EXTRACTION CHECK\n")
cat("========================================\n\n")

cat(
  "Number of selected expression columns: ",
  length(primary_tumor_column_indices),
  "\n",
  sep = ""
)

if (
  length(primary_tumor_column_indices) !=
  513
) {
  
  stop(
    "ERROR: Expression extraction did not return exactly 513 samples."
  )
}


############################################################
# 13. Create tumor-only expression matrix
############################################################

expr_primary_tumor <- expr_raw[
  ,
  c(
    1,
    primary_tumor_column_indices
  ),
  drop = FALSE
]


############################################################
# 14. Verify extracted matrix dimensions
############################################################

cat("\n========================================\n")
cat("PRIMARY TUMOR EXPRESSION MATRIX\n")
cat("========================================\n\n")

cat(
  "Rows (genes): ",
  nrow(expr_primary_tumor),
  "\n",
  sep = ""
)

cat(
  "Columns (including gene ID): ",
  ncol(expr_primary_tumor),
  "\n",
  sep = ""
)

cat(
  "Expression sample columns: ",
  ncol(expr_primary_tumor) - 1,
  "\n",
  sep = ""
)


############################################################
# 15. Verify gene ID column
############################################################

cat("\n========================================\n")
cat("GENE ID CHECK\n")
cat("========================================\n\n")

cat(
  "Gene ID column name: ",
  colnames(expr_primary_tumor)[1],
  "\n",
  sep = ""
)

cat(
  "Unique gene IDs: ",
  length(
    unique(
      expr_primary_tumor[[1]]
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Missing gene IDs: ",
  sum(
    is.na(
      expr_primary_tumor[[1]]
    )
  ),
  "\n",
  sep = ""
)


############################################################
# 16. Expression -> Clinical matching
############################################################

clinical_sample_ids <- clinical_clean$Sample


############################################################
# 17. Check clinical sample IDs
############################################################

cat("\n========================================\n")
cat("CLINICAL SAMPLE ID CHECK\n")
cat("========================================\n\n")

cat(
  "Clinical rows: ",
  nrow(clinical_clean),
  "\n",
  sep = ""
)

cat(
  "Unique clinical Sample IDs: ",
  length(
    unique(
      clinical_sample_ids
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Duplicated clinical Sample IDs: ",
  sum(
    duplicated(
      clinical_sample_ids
    )
  ),
  "\n",
  sep = ""
)

if (
  length(unique(clinical_sample_ids)) !=
  nrow(clinical_clean)
) {
  
  stop(
    "ERROR: Duplicate clinical Sample IDs detected."
  )
}


############################################################
# 18. Expression -> Clinical matching
############################################################

expr_not_in_clinical <- setdiff(
  primary_tumor_ids,
  clinical_sample_ids
)


############################################################
# 19. Clinical -> Expression matching
############################################################

clinical_not_in_expr <- setdiff(
  clinical_sample_ids,
  primary_tumor_ids
)


############################################################
# 20. Matching summary
############################################################

cat("\n========================================\n")
cat("BIDIRECTIONAL SAMPLE MATCHING\n")
cat("========================================\n\n")

cat(
  "Expression Primary Tumor samples: ",
  length(primary_tumor_ids),
  "\n",
  sep = ""
)

cat(
  "Clinical Sample IDs: ",
  length(clinical_sample_ids),
  "\n",
  sep = ""
)

cat(
  "Expression samples NOT in clinical table: ",
  length(expr_not_in_clinical),
  "\n",
  sep = ""
)

cat(
  "Clinical samples NOT in Primary Tumor expression cohort: ",
  length(clinical_not_in_expr),
  "\n",
  sep = ""
)


############################################################
# 21. Print unmatched samples if any
############################################################

if (length(expr_not_in_clinical) > 0) {
  
  cat(
    "\nExpression -> Clinical unmatched samples:\n"
  )
  
  print(
    expr_not_in_clinical
  )
}


if (length(clinical_not_in_expr) > 0) {
  
  cat(
    "\nClinical -> Expression unmatched samples:\n"
  )
  
  print(
    clinical_not_in_expr
  )
}


############################################################
# 22. Require exact bidirectional matching
############################################################

if (
  length(expr_not_in_clinical) != 0 ||
  length(clinical_not_in_expr) != 0
) {
  
  stop(
    "\nERROR: Bidirectional sample matching is NOT exact."
  )
}

cat(
  "\nRESULT: Exact bidirectional sample matching confirmed.\n"
)


############################################################
# 23. Check identical sample sets
############################################################

sample_sets_identical <- identical(
  sort(primary_tumor_ids),
  sort(clinical_sample_ids)
)

cat(
  "Identical sample sets: ",
  sample_sets_identical,
  "\n",
  sep = ""
)

if (!sample_sets_identical) {
  
  stop(
    "ERROR: Expression and clinical sample sets are not identical."
  )
}


############################################################
# 24. Check expression matrix column order
############################################################

expr_matrix_sample_ids <- colnames(
  expr_primary_tumor
)[-1]


############################################################
# 25. Verify extracted expression columns
############################################################

column_order_identical <- identical(
  expr_matrix_sample_ids,
  primary_tumor_ids
)

cat(
  "Expression extraction preserved original column order: ",
  column_order_identical,
  "\n",
  sep = ""
)


############################################################
# 26. Create clinical table in expression-column order
############################################################

clinical_clean_wgcna <- clinical_clean[
  match(
    expr_matrix_sample_ids,
    clinical_clean$Sample
  ),
  ,
  drop = FALSE
]


############################################################
# 27. Verify expression-clinical row alignment
############################################################

alignment_confirmed <- identical(
  expr_matrix_sample_ids,
  clinical_clean_wgcna$Sample
)

cat("\n========================================\n")
cat("EXPRESSION–CLINICAL ALIGNMENT\n")
cat("========================================\n\n")

cat(
  "Expression samples: ",
  length(expr_matrix_sample_ids),
  "\n",
  sep = ""
)

cat(
  "Clinical rows after matching: ",
  nrow(clinical_clean_wgcna),
  "\n",
  sep = ""
)

cat(
  "Sample order aligned: ",
  alignment_confirmed,
  "\n",
  sep = ""
)

if (!alignment_confirmed) {
  
  stop(
    "ERROR: Expression and clinical sample order is not aligned."
  )
}


############################################################
# 28. Verify Primary Tumor classification again
############################################################

cat("\n========================================\n")
cat("FINAL PRIMARY TUMOR CLASSIFICATION CHECK\n")
cat("========================================\n\n")

print(
  table(
    clinical_clean_wgcna$tumor_normal,
    useNA = "always"
  )
)

if (
  any(
    is.na(
      clinical_clean_wgcna$tumor_normal
    )
  )
) {
  
  stop(
    "ERROR: Missing tumor_normal classification detected."
  )
}

if (
  any(
    clinical_clean_wgcna$tumor_normal != "tumor"
  )
) {
  
  stop(
    "ERROR: Non-tumor samples detected in WGCNA cohort."
  )
}


############################################################
# 29. Final cohort summary
############################################################

cat("\n========================================\n")
cat("WORKFLOW-04-2A FINAL SUMMARY\n")
cat("========================================\n\n")

cat(
  "Input expression matrix: ",
  nrow(expr_raw),
  " genes × ",
  ncol(expr_raw) - 1,
  " samples\n",
  sep = ""
)

cat(
  "Primary Tumor samples extracted: ",
  ncol(expr_primary_tumor) - 1,
  "\n",
  sep = ""
)

cat(
  "Clinical samples matched: ",
  nrow(clinical_clean_wgcna),
  "\n",
  sep = ""
)

cat(
  "Expression -> Clinical unmatched: ",
  length(expr_not_in_clinical),
  "\n",
  sep = ""
)

cat(
  "Clinical -> Expression unmatched: ",
  length(clinical_not_in_expr),
  "\n",
  sep = ""
)

cat(
  "Sample order aligned: ",
  alignment_confirmed,
  "\n",
  sep = ""
)

cat(
  "Tumor classification: ",
  paste(
    unique(
      clinical_clean_wgcna$tumor_normal
    ),
    collapse = ", "
  ),
  "\n",
  sep = ""
)


############################################################
# 30. Save tumor-only expression matrix
############################################################

write.table(
  expr_primary_tumor,
  file = "result_dir/Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 31. Save aligned clinical table
############################################################

write.table(
  clinical_clean_wgcna,
  file = "result_dir/Workflow-04-2A_clinical_aligned.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 32. Save Primary Tumor sample list
############################################################

write.table(
  data.frame(
    Sample = expr_matrix_sample_ids
  ),
  file = "result_dir/Workflow-04-2A_primary_tumor_sample_list.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 33. Final message
############################################################

cat("\n========================================\n")
cat("WORKFLOW-04-2A COMPLETED\n")
cat("========================================\n\n")

cat(
  "Saved files:\n",
  "1. Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt\n",
  "2. Workflow-04-2A_clinical_aligned.txt\n",
  "3. Workflow-04-2A_primary_tumor_sample_list.txt\n\n"
)

cat(
  "No expression transformation or WGCNA preprocessing was performed.\n"
)  

#### Workflow-04-2B：FPKM-UQ expression data structure and distribution audit
############################################################
# Workflow-04-2B
# FPKM-UQ Expression Data Structure and Distribution Audit
#
# Input:
#   Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt
#
#   513 Primary Tumor samples
#   60,483 Ensembl gene IDs
#
# Purpose:
#   1. Verify expression matrix structure
#   2. Check numeric data type
#   3. Check NA / NaN / Inf
#   4. Summarize global expression distribution
#   5. Evaluate sample-level distributions
#   6. Evaluate gene-level distributions
#   7. Quantify zero / low-expression proportions
#
# IMPORTANT:
#   No transformation
#   No normalization
#   No gene filtering
#   No PCA
#   No WGCNA
############################################################


############################
# 0. Environment
############################

rm(list = ls())

options(stringsAsFactors = FALSE)

cat("\n========================================\n")
cat("Workflow-04-2B started\n")
cat("FPKM-UQ expression distribution audit\n")
cat("========================================\n\n")


############################
# 1. Set working directory
############################

setwd("input_dir/GitHub2_Paper1_tumor-only WGCNA")


############################
# 2. Input files
############################

expr_file <- "Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt"
clinical_file <- "Workflow-04-2A_clinical_aligned.txt"


############################
# 3. Read expression matrix
############################

expr_raw <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

cat("Expression matrix loaded.\n")


############################
# 4. Basic structure audit
############################

cat("\n----------------------------------------\n")
cat("PART A: Basic expression matrix structure\n")
cat("----------------------------------------\n\n")

cat("Rows:", nrow(expr_raw), "\n")
cat("Columns:", ncol(expr_raw), "\n")

cat(
  "Number of genes:",
  nrow(expr_raw),
  "\n"
)

cat(
  "Number of samples:",
  ncol(expr_raw) - 1,
  "\n"
)

cat(
  "First column name:",
  colnames(expr_raw)[1],
  "\n"
)

cat("\nFirst 10 gene IDs:\n")
print(head(expr_raw[[1]], 10))

cat("\nFirst 10 sample IDs:\n")
print(head(colnames(expr_raw)[-1], 10))


############################
# 5. Gene ID audit
############################

gene_id <- expr_raw[[1]]

cat("\nGene ID audit:\n")

cat(
  "Missing gene IDs:",
  sum(is.na(gene_id)),
  "\n"
)

cat(
  "Empty gene IDs:",
  sum(trimws(gene_id) == ""),
  "\n"
)

cat(
  "Unique gene IDs:",
  length(unique(gene_id)),
  "\n"
)

cat(
  "Duplicated gene IDs:",
  sum(duplicated(gene_id)),
  "\n"
)


############################
# 6. Sample ID audit
############################

sample_ids <- colnames(expr_raw)[-1]

cat("\nSample ID audit:\n")

cat(
  "Number of sample IDs:",
  length(sample_ids),
  "\n"
)

cat(
  "Unique sample IDs:",
  length(unique(sample_ids)),
  "\n"
)

cat(
  "Duplicated sample IDs:",
  sum(duplicated(sample_ids)),
  "\n"
)


############################
# 7. Convert expression part
#    to numeric matrix
############################

expr_mat <- as.matrix(
  expr_raw[, -1, drop = FALSE]
)

storage.mode(expr_mat) <- "numeric"


############################
# 8. Numeric data type audit
############################

cat("\n----------------------------------------\n")
cat("PART B: Numeric data type and missing values\n")
cat("----------------------------------------\n\n")

cat(
  "Matrix dimensions:",
  nrow(expr_mat),
  "genes x",
  ncol(expr_mat),
  "samples\n"
)

cat(
  "Storage mode:",
  storage.mode(expr_mat),
  "\n"
)

cat(
  "Class:",
  class(expr_mat),
  "\n"
)


############################
# 9. NA / NaN / Inf audit
############################

n_na <- sum(is.na(expr_mat))

n_nan <- sum(
  is.nan(expr_mat)
)

n_inf <- sum(
  is.infinite(expr_mat)
)

cat("\nMissing / invalid value audit:\n")

cat("NA:", n_na, "\n")
cat("NaN:", n_nan, "\n")
cat("Inf:", n_inf, "\n")


############################
# 10. Negative values
############################

n_negative <- sum(
  expr_mat < 0,
  na.rm = TRUE
)

cat(
  "Negative values:",
  n_negative,
  "\n"
)


############################
# 11. Zero values
############################

n_zero <- sum(
  expr_mat == 0,
  na.rm = TRUE
)

n_total <- sum(
  !is.na(expr_mat)
)

zero_fraction <- n_zero / n_total

cat(
  "Zero values:",
  n_zero,
  "\n"
)

cat(
  "Zero fraction:",
  sprintf("%.4f", zero_fraction),
  "\n"
)


############################
# 12. Global expression
#     distribution
############################

cat("\n----------------------------------------\n")
cat("PART C: Global FPKM-UQ expression distribution\n")
cat("----------------------------------------\n\n")


expr_values <- as.numeric(expr_mat)

expr_values <- expr_values[
  is.finite(expr_values)
]


global_summary <- quantile(
  expr_values,
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

print(global_summary)


cat("\nMean:\n")
print(mean(expr_values))

cat("\nSD:\n")
print(sd(expr_values))


############################
# 13. Expression value
#     summary table
############################

global_summary_df <- data.frame(
  Statistic = names(global_summary),
  Value = as.numeric(global_summary)
)

print(global_summary_df)


############################
# 14. Per-sample summary
############################

cat("\n----------------------------------------\n")
cat("PART D: Per-sample expression distribution\n")
cat("----------------------------------------\n\n")


sample_summary <- data.frame(
  Sample = sample_ids,
  
  Min = apply(
    expr_mat,
    2,
    min,
    na.rm = TRUE
  ),
  
  Q1 = apply(
    expr_mat,
    2,
    quantile,
    probs = 0.25,
    na.rm = TRUE
  ),
  
  Median = apply(
    expr_mat,
    2,
    median,
    na.rm = TRUE
  ),
  
  Mean = apply(
    expr_mat,
    2,
    mean,
    na.rm = TRUE
  ),
  
  Q3 = apply(
    expr_mat,
    2,
    quantile,
    probs = 0.75,
    na.rm = TRUE
  ),
  
  Max = apply(
    expr_mat,
    2,
    max,
    na.rm = TRUE
  ),
  
  ZeroFraction = colMeans(
    expr_mat == 0,
    na.rm = TRUE
  )
)


############################
# 15. Print sample summary
############################

cat("\nSummary of sample-level medians:\n")

print(
  summary(sample_summary$Median)
)

cat("\nSummary of sample-level means:\n")

print(
  summary(sample_summary$Mean)
)

cat("\nSummary of sample-level zero fractions:\n")

print(
  summary(sample_summary$ZeroFraction)
)


############################
# 16. Identify samples with
#     extreme median expression
############################

median_q <- quantile(
  sample_summary$Median,
  probs = c(
    0.01,
    0.99
  ),
  na.rm = TRUE
)

sample_median_extreme <- subset(
  sample_summary,
  Median < median_q[1] |
    Median > median_q[2]
)

cat("\nSamples outside 1st-99th percentile of sample median:\n")

print(
  sample_median_extreme
)


############################
# 17. Identify samples with
#     extreme zero fraction
############################

zero_q <- quantile(
  sample_summary$ZeroFraction,
  probs = c(
    0.01,
    0.99
  ),
  na.rm = TRUE
)

sample_zero_extreme <- subset(
  sample_summary,
  ZeroFraction < zero_q[1] |
    ZeroFraction > zero_q[2]
)

cat("\nSamples outside 1st-99th percentile of zero fraction:\n")

print(
  sample_zero_extreme
)


############################
# 18. Per-gene summary
############################

cat("\n----------------------------------------\n")
cat("PART E: Per-gene expression distribution\n")
cat("----------------------------------------\n\n")


gene_summary <- data.frame(
  GeneID = gene_id,
  
  Min = apply(
    expr_mat,
    1,
    min,
    na.rm = TRUE
  ),
  
  Q1 = apply(
    expr_mat,
    1,
    quantile,
    probs = 0.25,
    na.rm = TRUE
  ),
  
  Median = apply(
    expr_mat,
    1,
    median,
    na.rm = TRUE
  ),
  
  Mean = apply(
    expr_mat,
    1,
    mean,
    na.rm = TRUE
  ),
  
  Q3 = apply(
    expr_mat,
    1,
    quantile,
    probs = 0.75,
    na.rm = TRUE
  ),
  
  Max = apply(
    expr_mat,
    1,
    max,
    na.rm = TRUE
  ),
  
  ZeroFraction = rowMeans(
    expr_mat == 0,
    na.rm = TRUE
  )
)


############################
# 19. Gene zero fraction
############################

cat("\nGene-level zero fraction summary:\n")

print(
  summary(gene_summary$ZeroFraction)
)


############################
# 20. Number of genes with
#     different zero thresholds
############################

zero_thresholds <- c(
  0.25,
  0.50,
  0.75,
  0.90,
  0.95,
  0.99
)

gene_zero_counts <- data.frame(
  ZeroFractionThreshold = zero_thresholds,
  NumberOfGenes = sapply(
    zero_thresholds,
    function(x) {
      sum(
        gene_summary$ZeroFraction >= x,
        na.rm = TRUE
      )
    }
  )
)

cat("\nGenes exceeding different zero-fraction thresholds:\n")

print(
  gene_zero_counts
)


############################
# 21. Highly expressed genes
############################

cat("\nTop 20 genes by median FPKM-UQ:\n")

top_genes_median <- gene_summary[
  order(
    gene_summary$Median,
    decreasing = TRUE
  ),
]

print(
  head(
    top_genes_median,
    20
  )
)


############################
# 22. Low-expression genes
############################

cat("\nGenes with median FPKM-UQ = 0:\n")

n_median_zero <- sum(
  gene_summary$Median == 0,
  na.rm = TRUE
)

cat(
  n_median_zero,
  "\n"
)

cat(
  "Percentage:",
  sprintf(
    "%.2f%%",
    100 * n_median_zero / nrow(gene_summary)
  ),
  "\n"
)


############################
# 23. Distribution plots
############################

cat("\n----------------------------------------\n")
cat("PART F: Distribution plots\n")
cat("----------------------------------------\n\n")


############################
# 23.1 Histogram of raw FPKM-UQ
############################

hist(
  expr_values,
  breaks = 100,
  main = "Global FPKM-UQ Expression Distribution",
  xlab = "FPKM-UQ",
  ylab = "Frequency"
)


############################
# 23.2 Histogram focusing on
#      lower expression range
############################

hist(
  expr_values[
    expr_values <= quantile(
      expr_values,
      0.99
    )
  ],
  breaks = 100,
  main = "FPKM-UQ Distribution (<= 99th Percentile)",
  xlab = "FPKM-UQ",
  ylab = "Frequency"
)


############################
# 23.3 Sample median distribution
############################

hist(
  sample_summary$Median,
  breaks = 50,
  main = "Distribution of Sample Median FPKM-UQ",
  xlab = "Sample Median FPKM-UQ",
  ylab = "Number of Samples"
)


############################
# 23.4 Sample zero fraction
############################

hist(
  sample_summary$ZeroFraction,
  breaks = 50,
  main = "Distribution of Sample Zero Fraction",
  xlab = "Fraction of Zero Expression Values",
  ylab = "Number of Samples"
)


############################
# 24. Boxplot of sample distributions
############################

boxplot(
  as.data.frame(expr_mat),
  outline = FALSE,
  las = 2,
  cex.axis = 0.25,
  main = "Raw FPKM-UQ Distribution Across Samples",
  ylab = "FPKM-UQ"
)


############################
# 25. Sample median ranking
############################

sample_order <- order(
  sample_summary$Median
)

plot(
  sample_summary$Median[sample_order],
  type = "p",
  pch = 16,
  main = "Sample Median FPKM-UQ",
  xlab = "Samples ranked by median expression",
  ylab = "Median FPKM-UQ"
)


############################
# 26. Sample zero-fraction ranking
############################

sample_order_zero <- order(
  sample_summary$ZeroFraction
)

plot(
  sample_summary$ZeroFraction[
    sample_order_zero
  ],
  type = "p",
  pch = 16,
  main = "Sample Zero Fraction",
  xlab = "Samples ranked by zero fraction",
  ylab = "Zero fraction"
)


############################
# 27. Save audit tables
############################

write.table(
  global_summary_df,
  "result_dir/Workflow-04-2B_global_expression_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sample_summary,
  "result_dir/Workflow-04-2B_sample_expression_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sample_median_extreme,
  "result_dir/Workflow-04-2B_sample_median_extreme.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sample_zero_extreme,
  "result_dir/Workflow-04-2B_sample_zero_fraction_extreme.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_summary,
  "result_dir/Workflow-04-2B_gene_expression_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_zero_counts,
  "result_dir/Workflow-04-2B_gene_zero_fraction_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 28. Final audit report
############################

cat("\n========================================\n")
cat("Workflow-04-2B FINAL AUDIT\n")
cat("========================================\n\n")

cat(
  "Genes:",
  nrow(expr_mat),
  "\n"
)

cat(
  "Samples:",
  ncol(expr_mat),
  "\n"
)

cat(
  "NA:",
  n_na,
  "\n"
)

cat(
  "NaN:",
  n_nan,
  "\n"
)

cat(
  "Inf:",
  n_inf,
  "\n"
)

cat(
  "Negative values:",
  n_negative,
  "\n"
)

cat(
  "Zero fraction:",
  sprintf(
    "%.4f",
    zero_fraction
  ),
  "\n"
)

cat(
  "Median of all expression values:",
  median(
    expr_values
  ),
  "\n"
)

cat(
  "Mean of all expression values:",
  mean(
    expr_values
  ),
  "\n"
)

cat(
  "Maximum expression value:",
  max(
    expr_values
  ),
  "\n"
)

cat("\nWorkflow-04-2B completed.\n")
cat("NO transformation / filtering / WGCNA performed.\n") 


############################################################
# Workflow-04-2B-1
# Extreme-expression feature and sample provenance audit
#
# Input:
#   Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt
#
# Purpose:
#   1. Identify the gene/sample responsible for extreme values
#   2. Characterize the highest-expression features
#   3. Determine whether extreme values are mitochondrial,
#      ribosomal, protein-coding, or other features
#   4. Investigate TCGA-44 samples with unusually high
#      sample-level median expression
#   5. Compare sample-level distribution characteristics
#   6. Identify whether high sample medians are driven by
#      a small number of extreme features or by broader
#      transcriptome-wide elevation
#
# IMPORTANT:
#   NO transformation
#   NO normalization
#   NO gene filtering
#   NO sample removal
#   NO WGCNA
############################################################


############################
# 0. Environment
############################

options(stringsAsFactors = FALSE)

cat("\n========================================\n")
cat("Workflow-04-2B-1 started\n")
cat("Extreme-expression feature and sample provenance audit\n")
cat("========================================\n\n")


############################
# 1. Working directory
############################

setwd("input_dir/GitHub2_Paper1_tumor-only WGCNA")


############################
# 2. Input
############################

expr_file <- "Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt"


############################
# 3. Read expression matrix
############################

expr_raw <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

gene_id <- expr_raw[[1]]

sample_ids <- colnames(expr_raw)[-1]

expr_mat <- as.matrix(
  expr_raw[, -1, drop = FALSE]
)

storage.mode(expr_mat) <- "numeric"


cat("Expression matrix:\n")
cat("Genes:", nrow(expr_mat), "\n")
cat("Samples:", ncol(expr_mat), "\n\n")


############################################################
# PART A
# Identify absolute maximum expression values
############################################################

cat("----------------------------------------\n")
cat("PART A: Extreme expression values\n")
cat("----------------------------------------\n\n")


############################
# A1. Global maximum
############################

max_value <- max(
  expr_mat,
  na.rm = TRUE
)

max_index <- which(
  expr_mat == max_value,
  arr.ind = TRUE
)

max_info <- data.frame(
  GeneID = gene_id[max_index[, 1]],
  Sample = sample_ids[max_index[, 2]],
  Expression = expr_mat[max_index]
)

cat("Global maximum expression:\n")
print(max_info)


############################
# A2. Top 100 individual
#     gene × sample values
############################

expr_vec <- as.vector(expr_mat)

top100_index <- order(
  expr_vec,
  decreasing = TRUE
)[1:100]

top100_row <- ((top100_index - 1) %% nrow(expr_mat)) + 1

top100_col <- ((top100_index - 1) %/% nrow(expr_mat)) + 1

top100_extreme <- data.frame(
  Rank = 1:100,
  GeneID = gene_id[top100_row],
  Sample = sample_ids[top100_col],
  Expression = expr_vec[top100_index]
)

cat("\nTop 20 individual expression values:\n")
print(
  head(top100_extreme, 20)
)


############################################################
# PART B
# Extreme values by sample
############################################################

cat("\n----------------------------------------\n")
cat("PART B: Sample-level extreme values\n")
cat("----------------------------------------\n\n")


############################
# B1. Maximum expression
#     for every sample
############################

sample_max <- apply(
  expr_mat,
  2,
  max,
  na.rm = TRUE
)

sample_max_df <- data.frame(
  Sample = sample_ids,
  MaxExpression = sample_max
)

sample_max_df <- sample_max_df[
  order(
    sample_max_df$MaxExpression,
    decreasing = TRUE
  ),
]


cat("Top 20 samples by maximum expression:\n")

print(
  head(sample_max_df, 20)
)


############################
# B2. Number of genes above
#     extreme thresholds
############################

thresholds <- c(
  1e6,
  1e7,
  1e8,
  1e9
)

sample_extreme_counts <- data.frame(
  Sample = sample_ids
)

for (thr in thresholds) {
  
  sample_extreme_counts[
    paste0(
      "N_ge_",
      format(
        thr,
        scientific = FALSE
      )
    )
  ] <- colSums(
    expr_mat >= thr,
    na.rm = TRUE
  )
}


cat("\nNumber of genes exceeding extreme thresholds:\n")

print(
  head(
    sample_extreme_counts[
      order(
        sample_extreme_counts$N_ge_1000000000,
        decreasing = TRUE
      ),
    ],
    20
  )
)


############################################################
# PART C
# Investigate TCGA-44 samples
############################################################

cat("\n----------------------------------------\n")
cat("PART C: TCGA-44 sample audit\n")
cat("----------------------------------------\n\n")


tcga44_samples <- grep(
  "^TCGA-44-",
  sample_ids,
  value = TRUE
)

cat(
  "Number of TCGA-44 samples:",
  length(tcga44_samples),
  "\n\n"
)


############################
# C1. TCGA-44 summary
############################

tcga44_summary <- data.frame(
  Sample = tcga44_samples,
  Median = apply(
    expr_mat[, tcga44_samples, drop = FALSE],
    2,
    median,
    na.rm = TRUE
  ),
  Mean = apply(
    expr_mat[, tcga44_samples, drop = FALSE],
    2,
    mean,
    na.rm = TRUE
  ),
  Q1 = apply(
    expr_mat[, tcga44_samples, drop = FALSE],
    2,
    quantile,
    probs = 0.25,
    na.rm = TRUE
  ),
  Q3 = apply(
    expr_mat[, tcga44_samples, drop = FALSE],
    2,
    quantile,
    probs = 0.75,
    na.rm = TRUE
  ),
  Max = apply(
    expr_mat[, tcga44_samples, drop = FALSE],
    2,
    max,
    na.rm = TRUE
  ),
  ZeroFraction = colMeans(
    expr_mat[, tcga44_samples, drop = FALSE] == 0,
    na.rm = TRUE
  )
)

tcga44_summary <- tcga44_summary[
  order(
    tcga44_summary$Median,
    decreasing = TRUE
  ),
]

cat("TCGA-44 samples ranked by median expression:\n")

print(
  tcga44_summary
)


############################################################
# PART D
# Determine whether extreme sample medians
# are driven by many genes or a small number
############################################################

cat("\n----------------------------------------\n")
cat("PART D: Breadth of extreme expression\n")
cat("----------------------------------------\n\n")


############################
# D1. Selected TCGA-44 samples
############################

selected_samples <- c(
  "TCGA-44-2662-01",
  "TCGA-44-2656-01",
  "TCGA-44-6147-01",
  "TCGA-44-4112-01",
  "TCGA-44-6146-01",
  "TCGA-44-3917-01"
)

selected_samples <- intersect(
  selected_samples,
  sample_ids
)

cat(
  "Selected TCGA-44 samples found:",
  length(selected_samples),
  "\n"
)

print(selected_samples)


############################
# D2. For each selected sample,
#     calculate distribution
#     percentiles
############################

selected_distribution <- do.call(
  rbind,
  lapply(
    selected_samples,
    function(s) {
      
      x <- expr_mat[, s]
      
      data.frame(
        Sample = s,
        P50 = quantile(
          x,
          0.50,
          na.rm = TRUE
        ),
        P75 = quantile(
          x,
          0.75,
          na.rm = TRUE
        ),
        P90 = quantile(
          x,
          0.90,
          na.rm = TRUE
        ),
        P95 = quantile(
          x,
          0.95,
          na.rm = TRUE
        ),
        P99 = quantile(
          x,
          0.99,
          na.rm = TRUE
        ),
        P999 = quantile(
          x,
          0.999,
          na.rm = TRUE
        ),
        Max = max(
          x,
          na.rm = TRUE
        ),
        ZeroFraction = mean(
          x == 0,
          na.rm = TRUE
        )
      )
    }
  )
)

cat("\nDistribution of selected TCGA-44 samples:\n")

print(
  selected_distribution
)


############################################################
# PART E
# Identify the genes responsible for the extreme
# expression in selected TCGA-44 samples
############################################################

cat("\n----------------------------------------\n")
cat("PART E: Extreme genes in selected samples\n")
cat("----------------------------------------\n\n")


top_genes_selected <- list()


for (s in selected_samples) {
  
  x <- expr_mat[, s]
  
  idx <- order(
    x,
    decreasing = TRUE
  )[1:30]
  
  top_genes_selected[[s]] <- data.frame(
    Rank = 1:30,
    GeneID = gene_id[idx],
    Expression = x[idx]
  )
  
  cat("\n====================================\n")
  cat("Sample:", s, "\n")
  cat("====================================\n")
  
  print(
    top_genes_selected[[s]]
  )
}


############################################################
# PART F
# Determine whether extreme genes are mitochondrial
# or other genomic categories
############################################################

cat("\n----------------------------------------\n")
cat("PART F: Gene-category provenance\n")
cat("----------------------------------------\n\n")


############################
# F1. Ensembl gene IDs
#     starting with ENSG
############################

ensg_fraction <- mean(
  grepl(
    "^ENSG",
    gene_id
  )
)

cat(
  "Fraction of gene IDs matching ENSG:",
  sprintf(
    "%.4f",
    ensg_fraction
  ),
  "\n"
)


############################
# F2. Mitochondrial genes
#
# Ensembl mitochondrial genes
# commonly begin with ENSG000001988...
# but this is NOT sufficient for exact annotation.
#
# Here we only flag IDs matching the
# mitochondrial Ensembl pattern for initial audit.
############################

mito_flag <- grepl(
  "^ENSG000001988|^ENSG000002100",
  gene_id
)

cat(
  "Genes flagged by mitochondrial Ensembl pattern:",
  sum(mito_flag),
  "\n"
)


############################
# F3. Count extreme values
#     attributable to flagged genes
############################

top100_extreme$MitoPattern <- mito_flag[
  match(
    top100_extreme$GeneID,
    gene_id
  )
]

cat("\nTop 100 extreme values:\n")

print(
  table(
    top100_extreme$MitoPattern,
    useNA = "ifany"
  )
)


############################################################
# PART G
# Compare selected samples against all 513 samples
############################################################

cat("\n----------------------------------------\n")
cat("PART G: Relative position among all samples\n")
cat("----------------------------------------\n\n")


all_sample_summary <- data.frame(
  Sample = sample_ids,
  Median = apply(
    expr_mat,
    2,
    median,
    na.rm = TRUE
  ),
  Mean = apply(
    expr_mat,
    2,
    mean,
    na.rm = TRUE
  ),
  P75 = apply(
    expr_mat,
    2,
    quantile,
    probs = 0.75,
    na.rm = TRUE
  ),
  P90 = apply(
    expr_mat,
    2,
    quantile,
    probs = 0.90,
    na.rm = TRUE
  ),
  P95 = apply(
    expr_mat,
    2,
    quantile,
    probs = 0.95,
    na.rm = TRUE
  ),
  P99 = apply(
    expr_mat,
    2,
    quantile,
    probs = 0.99,
    na.rm = TRUE
  ),
  Max = apply(
    expr_mat,
    2,
    max,
    na.rm = TRUE
  ),
  ZeroFraction = colMeans(
    expr_mat == 0,
    na.rm = TRUE
  )
)


############################
# G1. Percentile rank of
#     selected samples
############################

selected_relative <- all_sample_summary[
  match(
    selected_samples,
    all_sample_summary$Sample
  ),
]

selected_relative$MedianPercentile <- sapply(
  selected_relative$Median,
  function(x) {
    mean(
      all_sample_summary$Median <= x
    )
  }
)

selected_relative$MeanPercentile <- sapply(
  selected_relative$Mean,
  function(x) {
    mean(
      all_sample_summary$Mean <= x
    )
  }
)

cat("\nRelative position of selected TCGA-44 samples:\n")

print(
  selected_relative
)


############################################################
# PART H
# Check whether the extreme pattern is restricted
# to TCGA-44
############################################################

cat("\n----------------------------------------\n")
cat("PART H: Sample provenance comparison\n")
cat("----------------------------------------\n\n")


############################
# H1. Compare maximum values
#     across TCGA patient prefixes
############################

patient_prefix <- sub(
  "^((TCGA-[0-9A-Z]+)-[0-9A-Z]+)-.*$",
  "\\1",
  sample_ids
)

sample_prefix_summary <- data.frame(
  Sample = sample_ids,
  Prefix = patient_prefix,
  Median = all_sample_summary$Median,
  Max = all_sample_summary$Max,
  ZeroFraction = all_sample_summary$ZeroFraction
)


############################
# H2. Patient prefix summary
############################

prefix_summary <- aggregate(
  cbind(
    Median,
    Max,
    ZeroFraction
  ) ~ Prefix,
  data = sample_prefix_summary,
  FUN = median
)

prefix_summary <- prefix_summary[
  order(
    prefix_summary$Median,
    decreasing = TRUE
  ),
]

cat("\nPatient-prefix summary ranked by median sample median:\n")

print(
  head(
    prefix_summary,
    20
  )
)


############################################################
# PART I
# Save audit results
############################################################

cat("\n----------------------------------------\n")
cat("PART I: Save audit tables\n")
cat("----------------------------------------\n\n")


write.table(
  max_info,
  "result_dir/Workflow-04-2B-1_global_maximum.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  top100_extreme,
  "result_dir/Workflow-04-2B-1_top100_extreme_values.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sample_max_df,
  "result_dir/Workflow-04-2B-1_sample_max_expression.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sample_extreme_counts,
  "result_dir/Workflow-04-2B-1_sample_extreme_expression_counts.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  tcga44_summary,
  "result_dir/Workflow-04-2B-1_TCGA44_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  selected_distribution,
  "result_dir/Workflow-04-2B-1_selected_sample_distribution.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  selected_relative,
  "result_dir/Workflow-04-2B-1_selected_sample_relative_position.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  prefix_summary,
  "result_dir/Workflow-04-2B-1_patient_prefix_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# PART J
# Final report
############################################################

cat("\n========================================\n")
cat("Workflow-04-2B-1 FINAL REPORT\n")
cat("========================================\n\n")

cat(
  "Global maximum:",
  max_value,
  "\n"
)

cat(
  "Number of TCGA-44 samples:",
  length(tcga44_samples),
  "\n"
)

cat(
  "Selected TCGA-44 samples analyzed:",
  length(selected_samples),
  "\n"
)

cat(
  "Mitochondrial-pattern genes:",
  sum(mito_flag),
  "\n"
)

cat("\nNo samples or genes were removed.\n")
cat("No transformation was performed.\n")
cat("No normalization was performed.\n")
cat("No WGCNA was performed.\n")

cat("\nWorkflow-04-2B-1 completed.\n") 


############################################################
# Workflow-04-2B-2
# Exact annotation of extreme Ensembl features
#
# Purpose:
# 1. Precisely annotate Ensembl IDs responsible for
#    extreme expression values
# 2. Determine gene symbol and biotype where available
# 3. Quantify mitochondrial / non-mitochondrial features
# 4. Specifically investigate:
#      ENSG00000202198
#      ENSG00000274012
# 5. No transformation
# 6. No normalization
# 7. No gene filtering
# 8. No sample removal
# 9. No WGCNA
############################################################


############################################################
# 1. Set working directory
############################################################

setwd("input_dir/GitHub2_Paper1_tumor-only WGCNA")


############################################################
# 2. Load required packages
############################################################

if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  stop("Package 'AnnotationDbi' is not installed.")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  stop("Package 'org.Hs.eg.db' is not installed.")
}

library(AnnotationDbi)
library(org.Hs.eg.db)


############################################################
# 3. Read expression matrix
############################################################

expr_file <- "result_dir/Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt"

expr_raw <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

gene_id <- expr_raw[[1]]

sample_ids <- colnames(expr_raw)[-1]

expr_mat <- as.matrix(
  expr_raw[, -1, drop = FALSE]
)

storage.mode(expr_mat) <- "numeric"


cat("========================================\n")
cat("Workflow-04-2B-2\n")
cat("Exact annotation of extreme Ensembl features\n")
cat("========================================\n\n")

cat("Genes:", nrow(expr_mat), "\n")
cat("Samples:", ncol(expr_mat), "\n")
cat("Unique Ensembl IDs:", length(unique(gene_id)), "\n\n")


############################################################
# 4. Read Top 100 extreme values from Workflow-04-2B-1
############################################################

top100_file <- "result_dir/Workflow-04-2B-1_top100_extreme_values.txt"

top100_extreme <- read.table(
  top100_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

cat("----------------------------------------\n")
cat("PART A: Top 100 extreme features\n")
cat("----------------------------------------\n\n")

cat(
  "Number of extreme observations:",
  nrow(top100_extreme),
  "\n"
)


############################################################
# 5. Extract unique Ensembl IDs from Top 100
############################################################

extreme_ensembl <- unique(
  top100_extreme$GeneID
)

cat(
  "Unique Ensembl IDs among Top 100:",
  length(extreme_ensembl),
  "\n\n"
)


############################################################
# 6. Add specifically selected high-impact features
############################################################

key_features <- c(
  "ENSG00000202198",
  "ENSG00000274012",
  "ENSG00000251562",
  "ENSG00000202538",
  "ENSG00000210082",
  "ENSG00000200087",
  "ENSG00000200156",
  "ENSG00000199568",
  "ENSG00000245532",
  "ENSG00000198804"
)

annotation_targets <- unique(
  c(
    extreme_ensembl,
    key_features
  )
)

cat(
  "Total annotation targets:",
  length(annotation_targets),
  "\n\n"
)


############################################################
# 7. Remove version suffix if present
############################################################

annotation_targets_clean <- sub(
  "\\..*$",
  "",
  annotation_targets
)

gene_id_clean <- sub(
  "\\..*$",
  "",
  gene_id
)


############################################################
# 8. Exact Ensembl annotation using org.Hs.eg.db
############################################################

cat("----------------------------------------\n")
cat("PART B: Ensembl annotation\n")
cat("----------------------------------------\n\n")


annotation_raw <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = annotation_targets_clean,
  keytype = "ENSEMBL",
  columns = c(
    "ENSEMBL",
    "SYMBOL",
    "GENENAME",
    "ENTREZID"
  )
)


############################################################
# 9. Remove duplicated annotation rows
############################################################

annotation_raw <- annotation_raw[
  !duplicated(
    annotation_raw[, c(
      "ENSEMBL",
      "SYMBOL",
      "GENENAME",
      "ENTREZID"
    )]
  ),
]


############################################################
# 10. Create one-row-per-Ensembl annotation table
############################################################

annotation_summary <- do.call(
  rbind,
  lapply(
    annotation_targets_clean,
    function(id) {
      
      x <- annotation_raw[
        annotation_raw$ENSEMBL == id,
        ,
        drop = FALSE
      ]
      
      if (nrow(x) == 0) {
        
        data.frame(
          ENSEMBL = id,
          SYMBOL = NA_character_,
          GENENAME = NA_character_,
          ENTREZID = NA_character_,
          AnnotationStatus = "Not_annotated",
          stringsAsFactors = FALSE
        )
        
      } else {
        
        symbol_values <- unique(
          na.omit(x$SYMBOL)
        )
        
        gene_name_values <- unique(
          na.omit(x$GENENAME)
        )
        
        entrez_values <- unique(
          na.omit(x$ENTREZID)
        )
        
        data.frame(
          ENSEMBL = id,
          SYMBOL = ifelse(
            length(symbol_values) > 0,
            paste(symbol_values, collapse = ";"),
            NA_character_
          ),
          GENENAME = ifelse(
            length(gene_name_values) > 0,
            paste(gene_name_values, collapse = ";"),
            NA_character_
          ),
          ENTREZID = ifelse(
            length(entrez_values) > 0,
            paste(entrez_values, collapse = ";"),
            NA_character_
          ),
          AnnotationStatus = "Annotated",
          stringsAsFactors = FALSE
        )
      }
    }
  )
)


############################################################
# 11. Add mitochondrial pattern annotation
#
# IMPORTANT:
# This is only a supplementary flag.
# It is NOT used as the formal annotation.
############################################################

annotation_summary$MitoPattern <- grepl(
  "^ENSG000001988|^ENSG000002100",
  annotation_summary$ENSEMBL
)


############################################################
# 12. Display annotation results
############################################################

cat("\nAnnotation summary:\n")

print(
  annotation_summary
)


############################################################
# 13. Specifically inspect key features
############################################################

cat("\n----------------------------------------\n")
cat("PART C: Key extreme features\n")
cat("----------------------------------------\n\n")

key_annotation <- annotation_summary[
  annotation_summary$ENSEMBL %in%
    sub(
      "\\..*$",
      "",
      key_features
    ),
  ,
  drop = FALSE
]

print(
  key_annotation
)


############################################################
# 14. Check whether key features are mitochondrial-pattern
############################################################

cat("\nMitochondrial-pattern status of key features:\n")

print(
  key_annotation[
    ,
    c(
      "ENSEMBL",
      "SYMBOL",
      "GENENAME",
      "MitoPattern"
    )
  ]
)

############################################################
# Workflow-04-2B-2
# PART 15–21: Robust annotation of Top 100 extreme observations
#
# IMPORTANT:
# Use match() instead of merge()
# to avoid column-name conflicts and preserve original Rank.
############################################################


############################################################
# 15. Map Top 100 extreme observations to annotation
############################################################

cat("\n----------------------------------------\n")
cat("PART E: Annotating Top 100 extreme observations\n")
cat("----------------------------------------\n\n")


# Start from the original Top 100 table
top100_annotated <- top100_extreme


# Clean Ensembl IDs
top100_annotated$ENSEMBL_clean <- sub(
  "\\..*$",
  "",
  top100_annotated$GeneID
)


# Match each Top 100 Ensembl ID
# to annotation_summary
match_idx <- match(
  top100_annotated$ENSEMBL_clean,
  annotation_summary$ENSEMBL
)


# Add annotation columns directly
top100_annotated$SYMBOL <- annotation_summary$SYMBOL[match_idx]

top100_annotated$GENENAME <- annotation_summary$GENENAME[match_idx]

top100_annotated$ENTREZID <- annotation_summary$ENTREZID[match_idx]

top100_annotated$MitoPattern <- annotation_summary$MitoPattern[match_idx]

top100_annotated$AnnotationStatus <-
  annotation_summary$AnnotationStatus[match_idx]


############################################################
# 16. Restore / verify original extreme-value rank
############################################################

top100_annotated <- top100_annotated[
  order(top100_annotated$Rank),
  ,
  drop = FALSE
]

rownames(top100_annotated) <- NULL


cat("Number of Top 100 observations:",
    nrow(top100_annotated), "\n")

cat("Unique Ensembl features:",
    length(unique(top100_annotated$ENSEMBL_clean)),
    "\n")

cat("Rank sequence preserved:",
    identical(
      top100_annotated$Rank,
      sort(top100_annotated$Rank)
    ),
    "\n\n")


############################################################
# 17. Count Top 100 by mitochondrial pattern
############################################################

cat("----------------------------------------\n")
cat("PART F: Top 100 category summary\n")
cat("----------------------------------------\n\n")


cat(
  "Mitochondrial-pattern TRUE:",
  sum(
    top100_annotated$MitoPattern %in% TRUE,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Mitochondrial-pattern FALSE:",
  sum(
    top100_annotated$MitoPattern %in% FALSE,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Mitochondrial-pattern NA:",
  sum(
    is.na(top100_annotated$MitoPattern)
  ),
  "\n\n"
)


############################################################
# 18. Show Top 100 annotated extreme observations
############################################################

cat("----------------------------------------\n")
cat("PART G: Top 100 annotated extreme observations\n")
cat("----------------------------------------\n\n")


print(
  top100_annotated[
    ,
    c(
      "Rank",
      "GeneID",
      "SYMBOL",
      "GENENAME",
      "Sample",
      "Expression",
      "MitoPattern",
      "AnnotationStatus"
    )
  ],
  row.names = FALSE
)


############################################################
# 19. Count unique extreme features by category
############################################################

cat("\n----------------------------------------\n")
cat("PART H: Unique extreme feature summary\n")
cat("----------------------------------------\n\n")


unique_extreme_annotation <- unique(
  top100_annotated[
    ,
    c(
      "ENSEMBL_clean",
      "SYMBOL",
      "GENENAME",
      "MitoPattern",
      "AnnotationStatus"
    )
  ]
)


cat(
  "Unique extreme Ensembl features:",
  nrow(unique_extreme_annotation),
  "\n\n"
)


print(
  unique_extreme_annotation,
  row.names = FALSE
)


############################################################
# 20. Identify unannotated / incompletely annotated features
############################################################

cat("\n----------------------------------------\n")
cat("PART I: Features requiring additional annotation\n")
cat("----------------------------------------\n\n")


unannotated_features <- unique_extreme_annotation[
  is.na(unique_extreme_annotation$SYMBOL) |
    unique_extreme_annotation$SYMBOL == "",
  ,
  drop = FALSE
]


cat(
  "Features without SYMBOL:",
  nrow(unannotated_features),
  "\n\n"
)


if (nrow(unannotated_features) > 0) {
  
  print(
    unannotated_features,
    row.names = FALSE
  )
  
} else {
  
  cat("No unannotated SYMBOLs among extreme features.\n")
}


############################################################
# 21. Save corrected Top 100 annotation
############################################################

write.table(
  top100_annotated,
  "result_dir/Workflow-04-2B-2_top100_annotated.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  unique_extreme_annotation,
  "result_dir/Workflow-04-2B-2_unique_extreme_feature_annotation.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat("\n========================================\n")
cat("PART E–I completed successfully.\n")
cat("No samples removed.\n")
cat("No genes removed.\n")
cat("No transformation performed.\n")
cat("No normalization performed.\n")
cat("No WGCNA performed.\n")
cat("========================================\n") 

############################################################
# Workflow-04-2B-2A
# Exact GENCODE / Ensembl biotype annotation
#
# Target features:
#   ENSG00000210082
#   ENSG00000251562
#
# Purpose:
#   Determine official gene name and gene biotype
#   for unresolved extreme expression features.
#
# IMPORTANT:
#   No samples removed
#   No genes removed
#   No expression transformation
#   No normalization
#   No WGCNA
############################################################


############################################################
# 1. Working directory
############################################################

setwd(
  "input_dir/GitHub2_Paper1_tumor-only WGCNA"
)


############################################################
# 2. Target Ensembl IDs
############################################################

target_genes <- c(
  "ENSG00000210082",
  "ENSG00000251562"
)


cat("\n========================================\n")
cat("Workflow-04-2B-2A\n")
cat("Exact GENCODE / Ensembl biotype annotation\n")
cat("========================================\n\n")

cat(
  "Target genes:\n",
  paste(target_genes, collapse = "\n"),
  "\n\n"
)


############################################################
# 3. Check biomaRt availability
############################################################

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  
  stop(
    paste0(
      "Package 'biomaRt' is not installed.\n",
      "Please install it first using:\n\n",
      "BiocManager::install('biomaRt')"
    )
  )
}

BiocManager::install('biomaRt')

library(biomaRt)


############################################################
# 4. Connect to Ensembl
############################################################

cat("----------------------------------------\n")
cat("PART A: Connecting to Ensembl BioMart\n")
cat("----------------------------------------\n\n")


ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)


cat("Ensembl BioMart connection established.\n\n")


############################################################
# 5. Check available attributes
############################################################

cat("----------------------------------------\n")
cat("PART B: Checking annotation attributes\n")
cat("----------------------------------------\n\n")


available_attributes <- listAttributes(
  ensembl
)


required_attributes <- c(
  "ensembl_gene_id",
  "external_gene_name",
  "description",
  "gene_biotype",
  "chromosome_name",
  "start_position",
  "end_position",
  "strand"
)


missing_attributes <- setdiff(
  required_attributes,
  available_attributes$name
)


if (length(missing_attributes) > 0) {
  
  cat(
    "WARNING: The following attributes are unavailable:\n"
  )
  
  print(missing_attributes)
  
  stop(
    "Required BioMart attributes are unavailable."
  )
  
}


cat("All required attributes are available.\n\n")


############################################################
# 6. Query Ensembl / GENCODE annotation
############################################################

cat("----------------------------------------\n")
cat("PART C: Querying Ensembl annotation\n")
cat("----------------------------------------\n\n")


annotation_ensembl <- getBM(
  
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "description",
    "gene_biotype",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand"
  ),
  
  filters = "ensembl_gene_id",
  
  values = target_genes,
  
  mart = ensembl
  
)


############################################################
# 7. Display raw annotation
############################################################

cat("\nRaw Ensembl annotation:\n\n")

print(
  annotation_ensembl,
  row.names = FALSE
)


############################################################
# 8. Check whether all target genes were found
############################################################

found_genes <- unique(
  annotation_ensembl$ensembl_gene_id
)


missing_genes <- setdiff(
  target_genes,
  found_genes
)


cat("\n----------------------------------------\n")
cat("PART D: Annotation completeness\n")
cat("----------------------------------------\n\n")


cat(
  "Target genes:",
  length(target_genes),
  "\n"
)

cat(
  "Successfully annotated:",
  length(found_genes),
  "\n"
)

cat(
  "Not found:",
  length(missing_genes),
  "\n\n"
)


if (length(missing_genes) > 0) {
  
  cat("Genes not found:\n")
  print(missing_genes)
  
} else {
  
  cat(
    "All target Ensembl IDs were successfully annotated.\n"
  )
}


############################################################
# 9. Create clean final annotation table
############################################################

annotation_final <- annotation_ensembl[
  ,
  c(
    "ensembl_gene_id",
    "external_gene_name",
    "description",
    "gene_biotype",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand"
  )
]


colnames(annotation_final) <- c(
  "ENSEMBL",
  "SYMBOL",
  "DESCRIPTION",
  "GENE_BIOTYPE",
  "CHROMOSOME",
  "START",
  "END",
  "STRAND"
)


############################################################
# 10. Add interpretation category
############################################################

annotation_final$Category <- NA_character_


annotation_final$Category[
  annotation_final$GENE_BIOTYPE == "Mt_rRNA"
] <- "Mitochondrial ribosomal RNA"


annotation_final$Category[
  annotation_final$GENE_BIOTYPE == "Mt_tRNA"
] <- "Mitochondrial transfer RNA"


annotation_final$Category[
  annotation_final$GENE_BIOTYPE == "protein_coding"
] <- "Protein-coding gene"


annotation_final$Category[
  annotation_final$GENE_BIOTYPE == "lncRNA"
] <- "Long non-coding RNA"


annotation_final$Category[
  annotation_final$GENE_BIOTYPE == "rRNA"
] <- "Ribosomal RNA"


annotation_final$Category[
  annotation_final$GENE_BIOTYPE == "snRNA"
] <- "Small nuclear RNA"


annotation_final$Category[
  annotation_final$GENE_BIOTYPE == "snoRNA"
] <- "Small nucleolar RNA"


############################################################
# 11. Display final annotation
############################################################

cat("\n----------------------------------------\n")
cat("PART E: Final biotype annotation\n")
cat("----------------------------------------\n\n")


print(
  annotation_final,
  row.names = FALSE
)


############################################################
# 12. Specific interpretation
############################################################

cat("\n----------------------------------------\n")
cat("PART F: Interpretation\n")
cat("----------------------------------------\n\n")


for (i in seq_len(nrow(annotation_final))) {
  
  cat(
    "\n",
    annotation_final$ENSEMBL[i],
    "\n",
    sep = ""
  )
  
  cat(
    "  SYMBOL: ",
    annotation_final$SYMBOL[i],
    "\n",
    sep = ""
  )
  
  cat(
    "  Gene biotype: ",
    annotation_final$GENE_BIOTYPE[i],
    "\n",
    sep = ""
  )
  
  cat(
    "  Category: ",
    annotation_final$Category[i],
    "\n",
    sep = ""
  )
  
  cat(
    "  Description: ",
    annotation_final$DESCRIPTION[i],
    "\n",
    sep = ""
  )
  
}


############################################################
# 13. Save result
############################################################

write.table(
  annotation_final,
  "result_dir/Workflow-04-2B-2A_GENCODE_Ensembl_biotype_annotation.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat("\n========================================\n")
cat("Workflow-04-2B-2A completed.\n")
cat("No samples removed.\n")
cat("No genes removed.\n")
cat("No transformation performed.\n")
cat("No normalization performed.\n")
cat("No WGCNA performed.\n")
cat("========================================\n")  


##Workflow-04-2C
##Expression Transformation + Gene Filtering 
############################################################
# Workflow-04-2C
# Expression Transformation + Gene Filtering
#
# Purpose:
#   1. Read tumor-only FPKM-UQ expression matrix
#   2. Verify matrix structure
#   3. Apply log2(FPKM-UQ + 1) transformation
#   4. Filter low-expression genes by prevalence
#   5. Annotate gene biotypes
#   6. Generate WGCNA candidate expression matrices
#
# Input:
#   Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt
#
# Output:
#   - log2 transformed expression matrix
#   - prevalence-filtered matrix
#   - gene filtering statistics
#   - gene biotype annotation
#   - final WGCNA candidate matrix
#
# IMPORTANT:
#   No sample removal
#   No PCA
#   No clustering
#   No WGCNA
############################################################


############################################################
# 1. Working directory
############################################################

setwd(
  "input_dir/GitHub2_Paper1_tumor-only WGCNA"
)


############################################################
# 2. Workflow information
############################################################

cat("\n========================================\n")
cat("Workflow-04-2C\n")
cat("Expression Transformation + Gene Filtering\n")
cat("========================================\n\n")


############################################################
# 3. Required packages
############################################################

required_packages <- c(
  "data.table",
  "biomaRt"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    stop(
      paste0(
        "Required package '",
        pkg,
        "' is not installed.\n"
      )
    )
    
  }
}


library(data.table)
library(biomaRt)


############################################################
# 4. Input file
############################################################

input_file <-
  "result_dir/Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt"


if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Input file not found:\n",
      input_file
    )
  )
  
}


cat("----------------------------------------\n")
cat("PART A: Reading input expression matrix\n")
cat("----------------------------------------\n\n")


expr_raw <- fread(
  input_file,
  data.table = FALSE,
  check.names = FALSE
)


############################################################
# 5. Basic structure check
############################################################

cat("Input dimensions:\n")
cat(
  "Rows (features): ",
  nrow(expr_raw),
  "\n",
  sep = ""
)

cat(
  "Columns: ",
  ncol(expr_raw),
  "\n",
  sep = ""
)

cat(
  "Samples: ",
  ncol(expr_raw) - 1,
  "\n\n",
  sep = ""
)


############################################################
# 6. Verify gene ID column
############################################################

if (colnames(expr_raw)[1] != "Tag") {
  
  stop(
    paste0(
      "Unexpected first column name: ",
      colnames(expr_raw)[1],
      "\nExpected: Tag"
    )
  )
  
}


gene_id <- expr_raw$Tag


############################################################
# 7. Verify gene IDs
############################################################

cat("----------------------------------------\n")
cat("PART B: Gene ID quality control\n")
cat("----------------------------------------\n\n")


cat(
  "Missing gene IDs: ",
  sum(is.na(gene_id)),
  "\n",
  sep = ""
)

cat(
  "Empty gene IDs: ",
  sum(gene_id == ""),
  "\n",
  sep = ""
)

cat(
  "Unique gene IDs: ",
  length(unique(gene_id)),
  "\n",
  sep = ""
)

cat(
  "Duplicated gene IDs: ",
  sum(duplicated(gene_id)),
  "\n\n",
  sep = ""
)


if (
  any(is.na(gene_id)) ||
  any(gene_id == "") ||
  any(duplicated(gene_id))
) {
  
  stop(
    "Gene ID quality control failed."
  )
  
}


############################################################
# 8. Extract expression matrix
############################################################

expr_mat <- as.matrix(
  expr_raw[, -1, drop = FALSE]
)

rownames(expr_mat) <- gene_id


############################################################
# 9. Verify sample IDs
############################################################

sample_ids <- colnames(expr_mat)


cat("----------------------------------------\n")
cat("PART C: Sample ID quality control\n")
cat("----------------------------------------\n\n")


cat(
  "Number of samples: ",
  length(sample_ids),
  "\n",
  sep = ""
)

cat(
  "Unique sample IDs: ",
  length(unique(sample_ids)),
  "\n",
  sep = ""
)

cat(
  "Duplicated sample IDs: ",
  sum(duplicated(sample_ids)),
  "\n\n",
  sep = ""
)


if (any(duplicated(sample_ids))) {
  
  stop(
    "Duplicated sample IDs detected."
  )
  
}


############################################################
# 10. Verify expression values
############################################################

cat("----------------------------------------\n")
cat("PART D: Expression matrix QC\n")
cat("----------------------------------------\n\n")


cat(
  "Storage mode: ",
  storage.mode(expr_mat),
  "\n",
  sep = ""
)

cat(
  "NA values: ",
  sum(is.na(expr_mat)),
  "\n",
  sep = ""
)

cat(
  "NaN values: ",
  sum(is.nan(expr_mat)),
  "\n",
  sep = ""
)

cat(
  "Inf values: ",
  sum(is.infinite(expr_mat)),
  "\n",
  sep = ""
)

cat(
  "Negative values: ",
  sum(expr_mat < 0),
  "\n\n",
  sep = ""
)


if (
  any(is.na(expr_mat)) ||
  any(is.nan(expr_mat)) ||
  any(is.infinite(expr_mat)) ||
  any(expr_mat < 0)
) {
  
  stop(
    "Expression matrix contains invalid values."
  )
  
}


############################################################
# 11. Record original dimensions
############################################################

n_genes_original <- nrow(expr_mat)
n_samples <- ncol(expr_mat)


############################################################
# 12. Log2 transformation
############################################################

cat("----------------------------------------\n")
cat("PART E: log2(FPKM-UQ + 1) transformation\n")
cat("----------------------------------------\n\n")


expr_log2 <- log2(
  expr_mat + 1
)


cat(
  "Transformation applied:\n",
  "log2(FPKM-UQ + 1)\n\n"
)


############################################################
# 13. Verify transformation
############################################################

cat(
  "Transformed matrix dimensions: ",
  nrow(expr_log2),
  " genes x ",
  ncol(expr_log2),
  " samples\n\n",
  sep = ""
)


cat("Transformed expression summary:\n")

print(
  summary(
    as.vector(expr_log2)
  )
)


############################################################
# 14. Save log2-transformed matrix
############################################################

log2_output <-
  "result_dir/Workflow-04-2C_expr_primary_tumor_log2_FPKM_UQ.txt"


write.table(
  cbind(
    Tag = rownames(expr_log2),
    as.data.frame(
      expr_log2,
      check.names = FALSE
    )
  ),
  file = log2_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat(
  "\nSaved log2-transformed matrix:\n",
  log2_output,
  "\n\n"
)


############################################################
# 15. Define expression prevalence threshold
############################################################

cat("----------------------------------------\n")
cat("PART F: Expression prevalence filtering\n")
cat("----------------------------------------\n\n")


# Expression threshold on original FPKM-UQ scale
expression_cutoff <- 1

# Minimum prevalence
prevalence_cutoff <- 0.20

# Minimum number of samples
min_samples <- ceiling(
  prevalence_cutoff * n_samples
)


cat(
  "Expression threshold: FPKM-UQ > ",
  expression_cutoff,
  "\n",
  sep = ""
)

cat(
  "Required prevalence: ",
  prevalence_cutoff * 100,
  "%\n",
  sep = ""
)

cat(
  "Number of samples required: ",
  min_samples,
  " / ",
  n_samples,
  "\n\n",
  sep = ""
)


############################################################
# 16. Calculate gene prevalence
############################################################

gene_detected_count <- rowSums(
  expr_mat > expression_cutoff
)


gene_prevalence <- gene_detected_count /
  n_samples


############################################################
# 17. Generate gene filtering table
############################################################

gene_filtering_table <- data.frame(
  
  ENSEMBL = rownames(expr_mat),
  
  Detected_Samples =
    gene_detected_count,
  
  Prevalence =
    gene_prevalence,
  
  Median_FPKM_UQ =
    apply(
      expr_mat,
      1,
      median
    ),
  
  Mean_FPKM_UQ =
    rowMeans(expr_mat),
  
  stringsAsFactors = FALSE
  
)


############################################################
# 18. Filtering decision
############################################################

gene_filtering_table$Keep <-
  gene_filtering_table$Detected_Samples >=
  min_samples


############################################################
# 19. Filtering statistics
############################################################

n_genes_kept <- sum(
  gene_filtering_table$Keep
)

n_genes_removed <- n_genes_original -
  n_genes_kept


cat(
  "Original genes: ",
  n_genes_original,
  "\n",
  sep = ""
)

cat(
  "Genes retained: ",
  n_genes_kept,
  "\n",
  sep = ""
)

cat(
  "Genes removed: ",
  n_genes_removed,
  "\n",
  sep = ""
)

cat(
  "Retention rate: ",
  round(
    100 * n_genes_kept / n_genes_original,
    2
  ),
  "%\n\n",
  sep = ""
)


############################################################
# 20. Distribution of gene prevalence
############################################################

cat("Gene prevalence summary:\n")

print(
  summary(
    gene_filtering_table$Prevalence
  )
)


cat("\n")


############################################################
# 21. Save gene filtering statistics
############################################################

filtering_output <-
  "result_dir/Workflow-04-2C_gene_filtering_statistics.txt"


write.table(
  gene_filtering_table,
  file = filtering_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat(
  "Saved gene filtering table:\n",
  filtering_output,
  "\n\n"
)


############################################################
# 22. Extract retained genes
############################################################

genes_keep <- gene_filtering_table$ENSEMBL[
  gene_filtering_table$Keep
]


############################################################
# 23. Create prevalence-filtered log2 matrix
############################################################

expr_log2_filtered <-
  expr_log2[
    genes_keep,
    ,
    drop = FALSE
  ]


cat("----------------------------------------\n")
cat("PART G: Prevalence-filtered expression matrix\n")
cat("----------------------------------------\n\n")


cat(
  "Filtered matrix dimensions: ",
  nrow(expr_log2_filtered),
  " genes x ",
  ncol(expr_log2_filtered),
  " samples\n\n",
  sep = ""
)


############################################################
# 24. Save prevalence-filtered matrix
############################################################

filtered_output <-
  "result_dir/Workflow-04-2C_expr_primary_tumor_log2_FPKM_UQ_prevalence_filtered.txt"


write.table(
  cbind(
    Tag = rownames(expr_log2_filtered),
    as.data.frame(
      expr_log2_filtered,
      check.names = FALSE
    )
  ),
  file = filtered_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat(
  "Saved prevalence-filtered matrix:\n",
  filtered_output,
  "\n\n"
)


############################################################
# 25. Gene biotype annotation
############################################################

cat("----------------------------------------\n")
cat("PART H: Ensembl gene biotype annotation\n")
cat("----------------------------------------\n\n")


genes_for_annotation <-
  rownames(expr_log2_filtered)


cat(
  "Genes requiring annotation: ",
  length(genes_for_annotation),
  "\n\n",
  sep = ""
)


############################################################
# 26. Connect to Ensembl BioMart
############################################################
library(biomaRt)

cat(
  "Connecting to Ensembl BioMart...\n"
)


ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)


############################################################
# 27. Query annotation
############################################################

annotation_biomart <- getBM(
  
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "description",
    "gene_biotype",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand"
  ),
  
  filters = "ensembl_gene_id",
  
  values = genes_for_annotation,
  
  mart = ensembl
  
)


############################################################
# 28. Check annotation completeness
############################################################

annotation_ids <-
  unique(
    annotation_biomart$ensembl_gene_id
  )


n_annotated <-
  length(
    annotation_ids
  )


n_unannotated <-
  length(
    setdiff(
      genes_for_annotation,
      annotation_ids
    )
  )


cat(
  "\nAnnotation results:\n"
)

cat(
  "Genes submitted: ",
  length(genes_for_annotation),
  "\n",
  sep = ""
)

cat(
  "Genes annotated: ",
  n_annotated,
  "\n",
  sep = ""
)

cat(
  "Genes not annotated: ",
  n_unannotated,
  "\n\n",
  sep = ""
)


############################################################
# 29. Deal with possible one-to-many annotation
############################################################

annotation_biomart <-
  annotation_biomart[
    !duplicated(
      annotation_biomart$ensembl_gene_id
    ),
    ,
    drop = FALSE
  ]


############################################################
# 30. Merge annotation with filtering statistics
############################################################

gene_annotation_final <-
  merge(
    
    gene_filtering_table,
    
    annotation_biomart,
    
    by.x = "ENSEMBL",
    
    by.y = "ensembl_gene_id",
    
    all.x = TRUE,
    
    sort = FALSE
    
  )


############################################################
# 31. Add annotation status
############################################################

gene_annotation_final$Annotation_Status <-
  ifelse(
    is.na(
      gene_annotation_final$gene_biotype
    ),
    "Not_annotated",
    "Annotated"
  )


############################################################
# 32. Add broad category
############################################################

gene_annotation_final$Category <-
  "Other"


gene_annotation_final$Category[
  gene_annotation_final$gene_biotype ==
    "protein_coding"
] <- "Protein-coding"


gene_annotation_final$Category[
  gene_annotation_final$gene_biotype ==
    "lncRNA"
] <- "Long non-coding RNA"


gene_annotation_final$Category[
  gene_annotation_final$gene_biotype ==
    "Mt_rRNA"
] <- "Mitochondrial rRNA"


gene_annotation_final$Category[
  gene_annotation_final$gene_biotype ==
    "Mt_tRNA"
] <- "Mitochondrial tRNA"


gene_annotation_final$Category[
  gene_annotation_final$gene_biotype ==
    "rRNA"
] <- "Ribosomal RNA"


gene_annotation_final$Category[
  gene_annotation_final$gene_biotype ==
    "snRNA"
] <- "Small nuclear RNA"


gene_annotation_final$Category[
  gene_annotation_final$gene_biotype ==
    "snoRNA"
] <- "Small nucleolar RNA"


############################################################
# 33. Display biotype distribution
############################################################

cat("----------------------------------------\n")
cat("PART I: Gene biotype distribution\n")
cat("----------------------------------------\n\n")


cat("Gene biotype counts:\n\n")

print(
  sort(
    table(
      gene_annotation_final$gene_biotype,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)


cat("\nBroad category counts:\n\n")

print(
  sort(
    table(
      gene_annotation_final$Category,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)


############################################################
# 34. Specifically inspect MT-RNR2 and MALAT1
############################################################

cat("\n----------------------------------------\n")
cat("PART J: Specific QC for MT-RNR2 and MALAT1\n")
cat("----------------------------------------\n\n")


target_ids <- c(
  "ENSG00000210082",
  "ENSG00000251562"
)


target_annotation <-
  gene_annotation_final[
    gene_annotation_final$ENSEMBL %in%
      target_ids,
    ,
    drop = FALSE
  ]


print(
  target_annotation,
  row.names = FALSE
)


############################################################
# 35. Create final WGCNA candidate set
############################################################

# IMPORTANT:
# At this stage we DO NOT automatically remove
# lncRNA or other non-protein-coding genes.
#
# Only the prevalence filter has been applied.
#
# Biotype information is retained for the next
# methodological decision.


wgcna_candidate_ids <-
  gene_annotation_final$ENSEMBL[
    gene_annotation_final$Keep
  ]


expr_wgcna_candidate <-
  expr_log2_filtered[
    wgcna_candidate_ids,
    ,
    drop = FALSE
  ]


############################################################
# 36. Verify final candidate matrix
############################################################

cat("\n----------------------------------------\n")
cat("PART K: WGCNA candidate matrix\n")
cat("----------------------------------------\n\n")


cat(
  "Genes: ",
  nrow(expr_wgcna_candidate),
  "\n",
  sep = ""
)

cat(
  "Samples: ",
  ncol(expr_wgcna_candidate),
  "\n",
  sep = ""
)

cat(
  "NA values: ",
  sum(is.na(expr_wgcna_candidate)),
  "\n",
  sep = ""
)

cat(
  "NaN values: ",
  sum(is.nan(expr_wgcna_candidate)),
  "\n",
  sep = ""
)

cat(
  "Inf values: ",
  sum(is.infinite(expr_wgcna_candidate)),
  "\n\n",
  sep = ""
)


############################################################
# 37. Save complete annotation table
############################################################

annotation_output <-
  "result_dir/Workflow-04-2C_gene_biotype_annotation.txt"


write.table(
  gene_annotation_final,
  file = annotation_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat(
  "Saved complete gene annotation table:\n",
  annotation_output,
  "\n\n"
)


############################################################
# 38. Save WGCNA candidate matrix
############################################################

candidate_output <-
  "result_dir/Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt"


write.table(
  cbind(
    Tag = rownames(expr_wgcna_candidate),
    as.data.frame(
      expr_wgcna_candidate,
      check.names = FALSE
    )
  ),
  file = candidate_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat(
  "Saved WGCNA candidate matrix:\n",
  candidate_output,
  "\n\n"
)


############################################################
# 39. Save concise workflow summary
############################################################

workflow_summary <- data.frame(
  
  Item = c(
    "Original genes",
    "Primary tumor samples",
    "Transformation",
    "Expression cutoff",
    "Prevalence cutoff",
    "Minimum samples required",
    "Genes retained after prevalence filtering",
    "Genes removed after prevalence filtering",
    "Retention rate (%)",
    "Annotation submitted",
    "Annotation successful",
    "Annotation not found",
    "Final WGCNA candidate genes",
    "Final WGCNA candidate samples"
  ),
  
  Value = c(
    n_genes_original,
    n_samples,
    "log2(FPKM-UQ + 1)",
    paste0(
      "FPKM-UQ > ",
      expression_cutoff
    ),
    paste0(
      prevalence_cutoff * 100,
      "%"
    ),
    min_samples,
    n_genes_kept,
    n_genes_removed,
    round(
      100 * n_genes_kept /
        n_genes_original,
      2
    ),
    length(genes_for_annotation),
    n_annotated,
    n_unannotated,
    nrow(expr_wgcna_candidate),
    ncol(expr_wgcna_candidate)
    
  ),
  
  stringsAsFactors = FALSE
)


summary_output <-
  "result_dir/Workflow-04-2C_summary.txt"


write.table(
  workflow_summary,
  file = summary_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 40. Final message
############################################################

############################################################
# Workflow-04-2C final verification
############################################################

cat("\n========================================\n")
cat("Workflow-04-2C FINAL VERIFICATION\n")
cat("========================================\n\n")

cat(
  "Original genes: ",
  n_genes_original,
  "\n",
  sep = ""
)

cat(
  "Retained genes: ",
  n_genes_kept,
  "\n",
  sep = ""
)

cat(
  "Removed genes: ",
  n_genes_removed,
  "\n",
  sep = ""
)

cat(
  "Retention rate: ",
  round(
    100 * n_genes_kept /
      n_genes_original,
    2
  ),
  "%\n",
  sep = ""
)

cat(
  "Samples: ",
  n_samples,
  "\n",
  sep = ""
)

cat(
  "Final candidate genes: ",
  nrow(expr_wgcna_candidate),
  "\n",
  sep = ""
)

cat(
  "Final candidate samples: ",
  ncol(expr_wgcna_candidate),
  "\n",
  sep = ""
)

cat(
  "NA: ",
  sum(is.na(expr_wgcna_candidate)),
  "\n",
  sep = ""
)

cat(
  "NaN: ",
  sum(is.nan(expr_wgcna_candidate)),
  "\n",
  sep = ""
)

cat(
  "Inf: ",
  sum(is.infinite(expr_wgcna_candidate)),
  "\n",
  sep = ""
)

cat("\n========================================\n")
cat("Verification completed.\n")
cat("========================================\n")

############################################################
# Workflow-04-2C-1
# Complete Ensembl / GENCODE Annotation Provenance Recovery
#
# Purpose:
#   1. Read the 04-2C prevalence-filtered candidate genes
#   2. Recover Ensembl gene annotation using GENCODE GTF
#   3. Resolve the 22,971 genes with NA biotype from BioMart
#   4. Cross-check BioMart and GENCODE annotations
#   5. Identify unresolved Ensembl IDs
#   6. Produce annotation statistics for deciding the final
#      WGCNA gene universe
#
# IMPORTANT:
#   - No samples are removed
#   - No genes are removed
#   - No expression transformation is performed
#   - No PCA / clustering / WGCNA is performed
#   - This workflow is annotation/provenance recovery only
############################################################


############################
# 0. Working directory
############################

setwd("input_dir/GitHub2_Paper1_tumor-only WGCNA")

cat("\n============================================\n")
cat("Workflow-04-2C-1\n")
cat("Ensembl / GENCODE annotation provenance recovery\n")
cat("============================================\n\n")


############################
# 1. Required packages
############################

required_packages <- c(
  "data.table",
  "dplyr",
  "stringr"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  
  library(
    pkg,
    character.only = TRUE
  )
}


############################
# 2. Input files
############################

expr_file <- paste0(
  "Workflow-04-2C_",
  "expr_primary_tumor_log2_FPKM_UQ_prevalence_filtered.txt"
)

biotype_file <- paste0(
  "Workflow-04-2C_",
  "gene_biotype_annotation.txt"
)


############################
# 3. Check input files
############################

if (!file.exists(expr_file)) {
  
  stop(
    paste0(
      "\nERROR: Expression file not found:\n",
      expr_file
    )
  )
}

if (!file.exists(biotype_file)) {
  
  stop(
    paste0(
      "\nERROR: BioMart annotation file not found:\n",
      biotype_file
    )
  )
}


cat("Input expression file:\n")
cat(expr_file, "\n\n")

cat("Input BioMart annotation file:\n")
cat(biotype_file, "\n\n")


############################
# 4. Read expression matrix
############################

cat("Reading prevalence-filtered expression matrix...\n")

expr <- fread(
  expr_file,
  data.table = FALSE,
  check.names = FALSE
)

cat("Expression matrix dimensions:\n")
cat(
  nrow(expr),
  "genes x",
  ncol(expr) - 1,
  "samples\n\n"
)


############################
# 5. Basic expression matrix QC
############################

if (!"Tag" %in% colnames(expr)) {
  
  stop(
    "\nERROR: First column 'Tag' was not found."
  )
}

gene_ids <- as.character(expr$Tag)

cat("Gene ID QC:\n")

cat(
  "Missing gene IDs:",
  sum(is.na(gene_ids)),
  "\n"
)

cat(
  "Empty gene IDs:",
  sum(gene_ids == ""),
  "\n"
)

cat(
  "Duplicated gene IDs:",
  sum(duplicated(gene_ids)),
  "\n"
)

cat(
  "Unique gene IDs:",
  length(unique(gene_ids)),
  "\n\n"
)

if (any(is.na(gene_ids))) {
  stop("ERROR: Missing gene IDs detected.")
}

if (any(gene_ids == "")) {
  stop("ERROR: Empty gene IDs detected.")
}

if (any(duplicated(gene_ids))) {
  stop("ERROR: Duplicated gene IDs detected.")
}


############################
# 6. Read previous BioMart annotation
############################

cat("Reading previous BioMart annotation...\n")

biomart_annot <- fread(
  biotype_file,
  data.table = FALSE,
  check.names = FALSE
)

cat(
  "BioMart annotation rows:",
  nrow(biomart_annot),
  "\n"
)

cat(
  "BioMart annotation columns:",
  ncol(biomart_annot),
  "\n\n"
)

cat("BioMart annotation columns:\n")
print(colnames(biomart_annot))
cat("\n")


############################
# 7. Identify Ensembl ID column
############################

possible_id_columns <- c(
  "ensembl_gene_id",
  "GeneID",
  "gene_id",
  "Ensembl",
  "ENSEMBL",
  "Tag"
)

id_col <- possible_id_columns[
  possible_id_columns %in% colnames(biomart_annot)
][1]

if (is.na(id_col)) {
  
  stop(
    paste0(
      "\nERROR: Could not identify Ensembl gene ID column in BioMart file.\n",
      "Available columns:\n",
      paste(colnames(biomart_annot), collapse = ", ")
    )
  )
}

cat(
  "BioMart Ensembl ID column:",
  id_col,
  "\n\n"
)


############################
# 8. Standardize BioMart ID
############################

biomart_annot$Ensembl_ID <- as.character(
  biomart_annot[[id_col]]
)

# Remove Ensembl version if present
biomart_annot$Ensembl_ID <- sub(
  "\\..*$",
  "",
  biomart_annot$Ensembl_ID
)

biomart_annot <- biomart_annot[
  !is.na(biomart_annot$Ensembl_ID) &
    biomart_annot$Ensembl_ID != "",
  ,
  drop = FALSE
]


############################
# 9. Restrict BioMart annotation
#    to current WGCNA candidate genes
############################

biomart_candidate <- biomart_annot[
  biomart_annot$Ensembl_ID %in% gene_ids,
  ,
  drop = FALSE
]

cat(
  "Candidate genes:",
  length(gene_ids),
  "\n"
)

cat(
  "Candidate genes found in BioMart annotation:",
  length(
    unique(
      biomart_candidate$Ensembl_ID
    )
  ),
  "\n\n"
)


############################
# 10. Inspect BioMart biotype
############################

biotype_col_candidates <- c(
  "gene_biotype",
  "Gene_Biotype",
  "biotype",
  "Biotype"
)

biotype_col <- biotype_col_candidates[
  biotype_col_candidates %in% colnames(biomart_candidate)
][1]

if (is.na(biotype_col)) {
  
  stop(
    paste0(
      "\nERROR: Could not identify BioMart biotype column.\n",
      "Available columns:\n",
      paste(
        colnames(biomart_candidate),
        collapse = ", "
      )
    )
  )
}

cat(
  "BioMart biotype column:",
  biotype_col,
  "\n\n"
)


############################
# 11. Count BioMart NA biotypes
############################

biomart_biotype <- as.character(
  biomart_candidate[[biotype_col]]
)

biomart_na <- is.na(biomart_biotype) |
  biomart_biotype == ""

cat("BioMart annotation status:\n")

cat(
  "Total candidate genes:",
  length(gene_ids),
  "\n"
)

cat(
  "BioMart biotype available:",
  sum(!biomart_na),
  "\n"
)

cat(
  "BioMart biotype NA/empty:",
  sum(biomart_na),
  "\n"
)

cat(
  "BioMart biotype available (%):",
  round(
    100 * sum(!biomart_na) /
      length(gene_ids),
    2
  ),
  "%\n"
)

cat(
  "BioMart biotype NA/empty (%):",
  round(
    100 * sum(biomart_na) /
      length(gene_ids),
    2
  ),
  "%\n\n"
)


############################
# 12. GENCODE GTF
#
# IMPORTANT:
# Download a GENCODE GTF corresponding to
# the annotation release used for the TCGA/GDC
# expression data if possible.
#
# Put the GTF file in:
#
#input_dir/GitHub2_Paper1_tumor-only WGCNA/
#
# Example:
# gencode.vXX.annotation.gtf.gz
#
############################

gencode_files <- list.files(
  path = ".",
  pattern = "^gencode.*\\.gtf(\\.gz)?$",
  full.names = TRUE,
  ignore.case = TRUE
)

cat("Searching for local GENCODE GTF...\n")

if (length(gencode_files) == 0) {
  
  cat("\n")
  cat("====================================================\n")
  cat("NO GENCODE GTF FOUND\n")
  cat("====================================================\n\n")
  
  cat(
    "Please download the appropriate GENCODE GTF and place it in:\n"
  )
  
  cat(
    getwd(),
    "\n\n"
  )
  
  cat(
    "Recommended filename pattern:\n"
  )
  
  cat(
    "gencode.vXX.annotation.gtf.gz\n\n"
  )
  
  cat(
    "The script will STOP here intentionally.\n"
  )
  
  cat(
    "No expression data have been modified.\n\n"
  )
  
  stop(
    "Workflow-04-2C-1 requires a GENCODE GTF file."
  )
}


cat(
  "GENCODE GTF candidate files found:\n"
)

print(gencode_files)

cat("\n")


############################################################
# Workflow-04-2C-1A
#
# GENCODE v22 Annotation Provenance Verification
#
# Purpose:
#   Determine whether the 60,483 Ensembl gene IDs in
#   Merge_RNA_seq_FPKM_UQ.txt are compatible with
#   GENCODE v22.
#
# IMPORTANT:
#   - No samples removed
#   - No genes removed
#   - No expression values modified
#   - No normalization
#   - No PCA
#   - No WGCNA
############################################################


############################
# 0. Working directory
############################

setwd(
  "input_dir/GitHub2_Paper1_tumor-only WGCNA"
)


cat("\n")
cat("====================================================\n")
cat("Workflow-04-2C-1A\n")
cat("GENCODE v22 annotation provenance verification\n")
cat("====================================================\n\n")


############################
# 1. Required packages
############################

required_packages <- c(
  "data.table",
  "dplyr",
  "stringr"
)


for (pkg in required_packages) {
  
  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {
    
    install.packages(pkg)
  }
  
  library(
    pkg,
    character.only = TRUE
  )
}


############################
# 2. Input expression file
############################

expr_file <- "input_dir/Merge_RNA_seq_FPKM_UQ.txt"

if (!file.exists(expr_file)) {
  
  stop(
    paste0(
      "ERROR: Expression file not found:\n",
      expr_file
    )
  )
}


############################
# 3. GENCODE v22 GTF
############################

gencode_file <- "gencode.v22.annotation.gtf.gz"


if (!file.exists(gencode_file)) {
  
  stop(
    paste0(
      "\nERROR: GENCODE v22 GTF not found:\n",
      gencode_file,
      "\n\nPlease place the official GENCODE v22 GTF in:\n",
      getwd()
    )
  )
}


cat(
  "Expression file:\n",
  expr_file,
  "\n\n"
)

cat(
  "GENCODE annotation:\n",
  gencode_file,
  "\n\n"
)


############################
# 4. Read expression matrix
############################

cat(
  "Reading expression matrix...\n"
)


expr <- fread(
  expr_file,
  data.table = FALSE,
  check.names = FALSE
)


cat(
  "Expression dimensions:",
  nrow(expr),
  "genes x",
  ncol(expr) - 1,
  "samples\n\n"
)


############################
# 5. Extract Ensembl IDs
############################

gene_ids <- as.character(
  expr[[1]]
)


gene_ids_clean <- sub(
  "\\..*$",
  "",
  gene_ids
)


cat(
  "Total gene IDs:",
  length(gene_ids_clean),
  "\n"
)

cat(
  "Unique gene IDs:",
  length(unique(gene_ids_clean)),
  "\n"
)

cat(
  "Duplicated gene IDs:",
  sum(duplicated(gene_ids_clean)),
  "\n\n"
)


############################
# 6. Check Ensembl format
############################

ensembl_valid <- grepl(
  "^ENSG[0-9]+$",
  gene_ids_clean
)


cat(
  "Valid Ensembl gene IDs:",
  sum(ensembl_valid),
  "\n"
)

cat(
  "Non-standard IDs:",
  sum(!ensembl_valid),
  "\n\n"
)


############################
# 7. Read GENCODE v22 GTF
############################

cat(
  "Reading GENCODE v22 GTF...\n"
)

cat(
  "This may take some time.\n\n"
)


gtf <- fread(
  gencode_file,
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  data.table = FALSE,
  quote = ""
)


if (ncol(gtf) != 9) {
  
  stop(
    paste0(
      "ERROR: GTF should have 9 columns.\n",
      "Detected: ",
      ncol(gtf)
    )
  )
}


colnames(gtf) <- c(
  "seqname",
  "source",
  "feature",
  "start",
  "end",
  "score",
  "strand",
  "frame",
  "attribute"
)


cat(
  "Total GTF rows:",
  nrow(gtf),
  "\n\n"
)


############################
# 8. Keep gene-level records
############################

gtf_gene <- gtf[
  gtf$feature == "gene",
  ,
  drop = FALSE
]


cat(
  "GENCODE gene-level records:",
  nrow(gtf_gene),
  "\n\n"
)


############################
# 9. Extract gene_id
############################

extract_attribute <- function(
    x,
    key
) {
  
  pattern <- paste0(
    key,
    " \"([^\"]+)\""
  )
  
  stringr::str_match(
    x,
    pattern
  )[, 2]
}


gtf_gene$Ensembl_ID <- extract_attribute(
  gtf_gene$attribute,
  "gene_id"
)


gtf_gene$gene_name <- extract_attribute(
  gtf_gene$attribute,
  "gene_name"
)


gtf_gene$gene_type <- extract_attribute(
  gtf_gene$attribute,
  "gene_type"
)


############################
# 10. Remove Ensembl version
############################

gtf_gene$Ensembl_ID <- sub(
  "\\..*$",
  "",
  gtf_gene$Ensembl_ID
)


############################
# 11. GENCODE annotation QC
############################

cat(
  "GENCODE annotation QC:\n\n"
)


cat(
  "Missing Ensembl IDs:",
  sum(
    is.na(gtf_gene$Ensembl_ID)
  ),
  "\n"
)


cat(
  "Missing gene_name:",
  sum(
    is.na(gtf_gene$gene_name)
  ),
  "\n"
)


cat(
  "Missing gene_type:",
  sum(
    is.na(gtf_gene$gene_type)
  ),
  "\n"
)


cat(
  "Unique Ensembl IDs:",
  length(
    unique(
      gtf_gene$Ensembl_ID
    )
  ),
  "\n\n"
)


############################
# 12. Match expression IDs
#     against GENCODE v22
############################

gencode_ids <- unique(
  gtf_gene$Ensembl_ID
)


matched <- gene_ids_clean %in%
  gencode_ids


n_matched <- sum(
  matched
)


n_unmatched <- sum(
  !matched
)


pct_matched <- round(
  100 *
    n_matched /
    length(gene_ids_clean),
  4
)


pct_unmatched <- round(
  100 *
    n_unmatched /
    length(gene_ids_clean),
  4
)


############################
# 13. Print main result
############################

cat("\n")
cat("====================================================\n")
cat("GENCODE v22 MATCHING RESULT\n")
cat("====================================================\n\n")


cat(
  "Expression genes:",
  length(gene_ids_clean),
  "\n"
)


cat(
  "Matched to GENCODE v22:",
  n_matched,
  "\n"
)


cat(
  "Unmatched:",
  n_unmatched,
  "\n"
)


cat(
  "Match rate:",
  pct_matched,
  "%\n"
)


cat(
  "Unmatched rate:",
  pct_unmatched,
  "%\n\n"
)


############################
# 14. Create annotation table
############################

gencode_annotation <- gtf_gene[
  ,
  c(
    "Ensembl_ID",
    "gene_name",
    "gene_type",
    "seqname",
    "start",
    "end",
    "strand"
  ),
  drop = FALSE
]


# Make sure one row per Ensembl ID
gencode_annotation <- gencode_annotation[
  !duplicated(
    gencode_annotation$Ensembl_ID
  ),
  ,
  drop = FALSE
]


############################
# 15. Merge annotation
############################

annotation_final <- data.frame(
  Ensembl_ID = gene_ids_clean,
  stringsAsFactors = FALSE
)


annotation_final <- annotation_final %>%
  
  left_join(
    gencode_annotation,
    by = "Ensembl_ID"
  )


############################
# 16. Annotation status
############################

annotation_final$Annotation_status <- ifelse(
  
  is.na(
    annotation_final$gene_type
  ),
  
  "Unresolved",
  
  "GENCODE_v22_annotated"
)


############################
# 17. Summary
############################

annotation_summary <- annotation_final %>%
  
  count(
    Annotation_status,
    name = "N_genes"
  )


annotation_summary$Percentage <- round(
  100 *
    annotation_summary$N_genes /
    nrow(annotation_final),
  2
)


cat(
  "Annotation status:\n\n"
)

print(
  annotation_summary,
  row.names = FALSE
)


############################
# 18. gene_type distribution
############################

cat(
  "\nGENCODE v22 gene_type distribution:\n\n"
)


gene_type_summary <- annotation_final %>%
  
  count(
    gene_type,
    sort = TRUE,
    name = "N_genes"
  )


gene_type_summary$Percentage <- round(
  100 *
    gene_type_summary$N_genes /
    nrow(annotation_final),
  2
)


print(
  gene_type_summary,
  row.names = FALSE
)


############################
# 19. Check important extreme IDs
############################

important_ids <- c(
  
  "ENSG00000274012",
  "ENSG00000202198",
  "ENSG00000210082",
  "ENSG00000251562",
  
  "ENSG00000198804",
  "ENSG00000198886",
  "ENSG00000198899",
  "ENSG00000198763",
  "ENSG00000198938",
  "ENSG00000198712",
  "ENSG00000198888",
  "ENSG00000198840",
  
  "ENSG00000161055",
  "ENSG00000087086",
  "ENSG00000124107",
  "ENSG00000211895",
  "ENSG00000163220",
  "ENSG00000185303",
  "ENSG00000211592",
  "ENSG00000202538",
  "ENSG00000200087",
  "ENSG00000200156",
  "ENSG00000199568",
  "ENSG00000245532"
)


important_annotation <- annotation_final %>%
  
  filter(
    Ensembl_ID %in% important_ids
  )


cat("\n")
cat("====================================================\n")
cat("IMPORTANT EXTREME-FEATURE ANNOTATION\n")
cat("====================================================\n\n")


print(
  important_annotation,
  row.names = FALSE
)


############################
# 20. Specifically check
#     MT-RNR2 and MALAT1
############################

special_ids <- c(
  "ENSG00000210082",
  "ENSG00000251562"
)


special_annotation <- annotation_final %>%
  
  filter(
    Ensembl_ID %in% special_ids
  )


cat("\n")
cat(
  "MT-RNR2 / MALAT1 verification:\n\n"
)


print(
  special_annotation,
  row.names = FALSE
)


############################
# 21. Save complete annotation
############################

output_annotation <- paste0(
  "Workflow-04-2C-1A_",
  "GENCODE_v22_complete_annotation.txt"
)


write.table(
  annotation_final,
  output_annotation,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)


############################
# 22. Save gene type summary
############################

output_gene_type <- paste0(
  "Workflow-04-2C-1A_",
  "GENCODE_v22_gene_type_summary.txt"
)


write.table(
  gene_type_summary,
  output_gene_type,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)


############################
# 23. Save important features
############################

output_important <- paste0(
  "Workflow-04-2C-1A_",
  "GENCODE_v22_important_feature_annotation.txt"
)


write.table(
  important_annotation,
  output_important,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)


############################
# 24. Save unmatched IDs
############################

unmatched_ids <- data.frame(
  Ensembl_ID =
    gene_ids_clean[
      !matched
    ]
)


output_unmatched <- paste0(
  "Workflow-04-2C-1A_",
  "GENCODE_v22_unmatched_Ensembl_IDs.txt"
)


write.table(
  unmatched_ids,
  output_unmatched,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 25. Save summary
############################

summary_lines <- c(
  
  "Workflow-04-2C-1A",
  "GENCODE v22 Annotation Provenance Verification",
  "============================================",
  "",
  
  paste(
    "Expression genes:",
    length(gene_ids_clean)
  ),
  
  paste(
    "GENCODE v22 matched:",
    n_matched
  ),
  
  paste(
    "GENCODE v22 unmatched:",
    n_unmatched
  ),
  
  paste(
    "Match rate:",
    pct_matched,
    "%"
  ),
  
  paste(
    "Unmatched rate:",
    pct_unmatched,
    "%"
  ),
  
  "",
  
  "No samples removed.",
  "No genes removed.",
  "No expression values modified.",
  "No normalization performed.",
  "No WGCNA performed."
)


output_summary <- paste0(
  "Workflow-04-2C-1A_summary.txt"
)


writeLines(
  summary_lines,
  output_summary
)


############################
# 26. Final output
############################

cat("\n")
cat("====================================================\n")
cat("Workflow-04-2C-1A COMPLETED\n")
cat("====================================================\n\n")


cat(
  "GENCODE v22 match rate:",
  pct_matched,
  "%\n"
)


cat(
  "Unmatched:",
  n_unmatched,
  "\n\n"
)


cat(
  "Files saved:\n"
)

cat(
  output_annotation,
  "\n"
)

cat(
  output_gene_type,
  "\n"
)

cat(
  output_important,
  "\n"
)

cat(
  output_unmatched,
  "\n"
)

cat(
  output_summary,
  "\n\n"
)


cat(
  "NEXT STEP:\n"
)

cat(
  "Review the GENCODE v22 matching result before\n"
)

cat(
  "constructing the final WGCNA gene universe.\n"
)

cat("\n")

###Workflow-04-2C-2

### GENCODE v22–based gene-type composition and preliminary variance/MAD QC
############################################################
# Workflow-04-2C-2
# GENCODE v22-based gene-type composition
# + mitochondrial/rRNA/non-coding annotation
# + preliminary variance/MAD QC
#
# IMPORTANT:
# This workflow is DESCRIPTIVE ONLY.
# NO GENES ARE REMOVED.
# NO SAMPLE IS REMOVED.
# NO WGCNA IS PERFORMED.
############################################################


############################################################
# 0. Environment
############################################################

rm(list = ls())

options(stringsAsFactors = FALSE)

setwd("input_dir/GitHub2_Paper1_tumor-only WGCNA")

cat("\n")
cat("====================================================\n")
cat("Workflow-04-2C-2\n")
cat("GENCODE v22-based gene composition and QC\n")
cat("====================================================\n\n")


############################################################
# 1. Load packages
############################################################

required_packages <- c(
  "data.table",
  "dplyr"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    install.packages(pkg)
    
  }
  
}

library(data.table)
library(dplyr)


############################################################
# 2. Input files
############################################################

expr_file <- file.path(
  getwd(),
  "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt"
)

annotation_file <- file.path(
  getwd(),
  "Workflow-04-2C-1A_GENCODE_v22_complete_annotation.txt"
)


############################################################
# 3. Check input files
############################################################

cat("Checking input files...\n\n")

if (!file.exists(expr_file)) {
  
  stop(
    "ERROR: Expression input file not found:\n",
    expr_file
  )
  
}

if (!file.exists(annotation_file)) {
  
  stop(
    "ERROR: GENCODE v22 annotation file not found:\n",
    annotation_file
  )
  
}

cat("Expression file:\n", expr_file, "\n\n")
cat("Annotation file:\n", annotation_file, "\n\n")


############################################################
# 4. Read expression matrix
############################################################

cat("Reading expression matrix...\n")

expr_df <- fread(
  expr_file,
  data.table = FALSE,
  check.names = FALSE
)

cat(
  "Expression dimensions:",
  nrow(expr_df),
  "genes x",
  ncol(expr_df) - 1,
  "samples\n\n"
)


############################################################
# 5. Basic expression matrix QC
############################################################

gene_ids <- expr_df[[1]]

expr_mat <- as.matrix(
  expr_df[, -1, drop = FALSE]
)

rownames(expr_mat) <- gene_ids

storage.mode(expr_mat) <- "numeric"


cat("Basic expression QC:\n")

cat(
  "Number of genes:",
  nrow(expr_mat),
  "\n"
)

cat(
  "Number of samples:",
  ncol(expr_mat),
  "\n"
)

cat(
  "Unique gene IDs:",
  length(unique(rownames(expr_mat))),
  "\n"
)

cat(
  "NA values:",
  sum(is.na(expr_mat)),
  "\n"
)

cat(
  "NaN values:",
  sum(is.nan(expr_mat)),
  "\n"
)

cat(
  "Inf values:",
  sum(is.infinite(expr_mat)),
  "\n\n"
)


############################################################
# 6. Read GENCODE v22 annotation
############################################################

cat("Reading GENCODE v22 annotation...\n")

annotation <- fread(
  annotation_file,
  data.table = FALSE,
  check.names = FALSE
)

cat(
  "Annotation rows:",
  nrow(annotation),
  "\n\n"
)

cat("Annotation columns:\n")

print(
  colnames(annotation)
)

cat("\n")


############################################################
# 7. Check required annotation columns
############################################################

required_annotation_columns <- c(
  "Ensembl_ID",
  "gene_name",
  "gene_type"
)

missing_annotation_columns <-
  setdiff(
    required_annotation_columns,
    colnames(annotation)
  )

if (length(missing_annotation_columns) > 0) {
  
  stop(
    "ERROR: Missing annotation columns: ",
    paste(
      missing_annotation_columns,
      collapse = ", "
    )
  )
  
}


############################################################
# 8. Remove Ensembl version suffix if present
############################################################

annotation$Ensembl_ID_clean <-
  sub(
    "\\..*$",
    "",
    annotation$Ensembl_ID
  )

expr_gene_ids_clean <-
  sub(
    "\\..*$",
    "",
    rownames(expr_mat)
  )


############################################################
# 9. Check expression ↔ GENCODE matching
############################################################

annotation_match_index <-
  match(
    expr_gene_ids_clean,
    annotation$Ensembl_ID_clean
  )

n_matched <-
  sum(
    !is.na(annotation_match_index)
  )

n_unmatched <-
  sum(
    is.na(annotation_match_index)
  )

match_rate <-
  100 * n_matched / length(expr_gene_ids_clean)


cat("====================================================\n")
cat("EXPRESSION ↔ GENCODE v22 MATCHING\n")
cat("====================================================\n\n")

cat(
  "Expression genes:",
  length(expr_gene_ids_clean),
  "\n"
)

cat(
  "Matched:",
  n_matched,
  "\n"
)

cat(
  "Unmatched:",
  n_unmatched,
  "\n"
)

cat(
  "Match rate:",
  round(match_rate, 4),
  "%\n\n"
)

if (n_unmatched > 0) {
  
  warning(
    "Some expression genes were not matched to GENCODE v22."
  )
  
}


############################################################
# 10. Build annotation for the 39,280 retained genes
############################################################

gene_annotation <-
  annotation[
    annotation_match_index,
    c(
      "Ensembl_ID",
      "gene_name",
      "gene_type"
    ),
    drop = FALSE
  ]

rownames(gene_annotation) <-
  rownames(expr_mat)


############################################################
# 11. Verify annotation alignment
############################################################

alignment_check <-
  all(
    sub(
      "\\..*$",
      "",
      gene_annotation$Ensembl_ID
    ) ==
      expr_gene_ids_clean
  )

cat(
  "Annotation order aligned to expression:",
  alignment_check,
  "\n\n"
)

if (!alignment_check) {
  
  stop(
    "ERROR: Annotation order does not match expression matrix."
  )
  
}


############################################################
# 12. Calculate gene-level variance and MAD
############################################################

cat("Calculating gene-level variance and MAD...\n")

gene_variance <-
  apply(
    expr_mat,
    1,
    var,
    na.rm = TRUE
  )

gene_sd <-
  apply(
    expr_mat,
    1,
    sd,
    na.rm = TRUE
  )

gene_mad <-
  apply(
    expr_mat,
    1,
    mad,
    na.rm = TRUE
  )

gene_mean <-
  rowMeans(
    expr_mat,
    na.rm = TRUE
  )

gene_median <-
  apply(
    expr_mat,
    1,
    median,
    na.rm = TRUE
  )


############################################################
# 13. Detection / prevalence information
############################################################

# Because this matrix was already filtered using
# FPKM-UQ > 1 in >=20% samples before log2 transformation,
# detection prevalence is calculated here only for QC reporting.

# log2(FPKM-UQ + 1) > log2(2) corresponds to FPKM-UQ > 1

detection_threshold <-
  log2(2)

gene_detected_n <-
  rowSums(
    expr_mat > detection_threshold,
    na.rm = TRUE
  )

gene_detected_fraction <-
  gene_detected_n / ncol(expr_mat)


############################################################
# 14. Construct QC table
############################################################

gene_qc <-
  data.frame(
    
    Ensembl_ID =
      rownames(expr_mat),
    
    gene_name =
      gene_annotation$gene_name,
    
    gene_type =
      gene_annotation$gene_type,
    
    Mean =
      gene_mean,
    
    Median =
      gene_median,
    
    SD =
      gene_sd,
    
    Variance =
      gene_variance,
    
    MAD =
      gene_mad,
    
    Detected_N =
      gene_detected_n,
    
    Detected_fraction =
      gene_detected_fraction,
    
    stringsAsFactors = FALSE
    
  )


############################################################
# 15. Add broad annotation categories
############################################################

gene_qc$Broad_category <- "Other"


############################################################
# Protein-coding
############################################################

gene_qc$Broad_category[
  gene_qc$gene_type ==
    "protein_coding"
] <- "Protein-coding"


############################################################
# Long non-coding RNA
############################################################

lncRNA_types <- c(
  "lincRNA",
  "antisense",
  "sense_intronic",
  "sense_overlapping",
  "3prime_overlapping_ncrna",
  "processed_transcript",
  "macro_lncRNA"
)

gene_qc$Broad_category[
  gene_qc$gene_type %in% lncRNA_types
] <- "Long non-coding RNA"


############################################################
# Small RNA
############################################################

small_RNA_types <- c(
  "miRNA",
  "snRNA",
  "snoRNA",
  "scaRNA",
  "sRNA"
)

gene_qc$Broad_category[
  gene_qc$gene_type %in% small_RNA_types
] <- "Small RNA"


############################################################
# Mitochondrial RNA
############################################################

gene_qc$Broad_category[
  gene_qc$gene_type %in%
    c(
      "Mt_tRNA",
      "Mt_rRNA"
    )
] <- "Mitochondrial RNA"


############################################################
# rRNA
############################################################

gene_qc$Broad_category[
  gene_qc$gene_type %in%
    c(
      "rRNA"
    )
] <- "Ribosomal RNA"


############################################################
# Pseudogene
############################################################

pseudogene_pattern <-
  grepl(
    "pseudogene",
    gene_qc$gene_type,
    ignore.case = TRUE
  )

gene_qc$Broad_category[
  pseudogene_pattern
] <- "Pseudogene"


############################################################
# Immune receptor genes
############################################################

immune_receptor_pattern <-
  grepl(
    "^(IG_|TR_)",
    gene_qc$gene_type
  )

gene_qc$Broad_category[
  immune_receptor_pattern
] <- "Immune receptor"


############################################################
# TEC
############################################################

gene_qc$Broad_category[
  gene_qc$gene_type == "TEC"
] <- "TEC"


############################################################
# 16. Recalculate priority categories
############################################################

# Mitochondrial genes are given priority over other categories.

gene_qc$Mitochondrial <-
  gene_qc$gene_type %in%
  c(
    "Mt_tRNA",
    "Mt_rRNA"
  ) |
  grepl(
    "^MT-",
    gene_qc$gene_name
  )


gene_qc$rRNA_related <-
  gene_qc$gene_type %in%
  c(
    "rRNA",
    "Mt_rRNA"
  )


gene_qc$Non_protein_coding <-
  gene_qc$gene_type !=
  "protein_coding"


############################################################
# 17. Gene-type composition
############################################################

gene_type_summary <-
  gene_qc %>%
  
  count(
    gene_type,
    sort = TRUE,
    name = "N_genes"
  ) %>%
  
  mutate(
    Percentage =
      100 * N_genes / nrow(gene_qc)
  ) %>%
  
  mutate(
    Percentage =
      round(
        Percentage,
        2
      )
  )


cat("\n")
cat("====================================================\n")
cat("GENCODE v22 GENE-TYPE COMPOSITION\n")
cat("====================================================\n\n")

print(
  gene_type_summary,
  row.names = FALSE
)


############################################################
# 18. Broad-category composition
############################################################

broad_category_summary <-
  gene_qc %>%
  
  count(
    Broad_category,
    sort = TRUE,
    name = "N_genes"
  ) %>%
  
  mutate(
    Percentage =
      round(
        100 * N_genes / nrow(gene_qc),
        2
      )
  )


cat("\n")
cat("====================================================\n")
cat("BROAD GENE CATEGORY COMPOSITION\n")
cat("====================================================\n\n")

print(
  broad_category_summary,
  row.names = FALSE
)


############################################################
# 19. Mitochondrial / rRNA / protein-coding summary
############################################################

special_category_summary <-
  data.frame(
    
    Category = c(
      "All retained genes",
      "Protein-coding",
      "Non-protein-coding",
      "Mitochondrial genes",
      "rRNA-related genes",
      "Long non-coding RNA",
      "Small RNA",
      "Pseudogenes",
      "TEC"
    ),
    
    N_genes = c(
      
      nrow(gene_qc),
      
      sum(
        gene_qc$gene_type ==
          "protein_coding"
      ),
      
      sum(
        gene_qc$Non_protein_coding
      ),
      
      sum(
        gene_qc$Mitochondrial
      ),
      
      sum(
        gene_qc$rRNA_related
      ),
      
      sum(
        gene_qc$Broad_category ==
          "Long non-coding RNA"
      ),
      
      sum(
        gene_qc$Broad_category ==
          "Small RNA"
      ),
      
      sum(
        gene_qc$Broad_category ==
          "Pseudogene"
      ),
      
      sum(
        gene_qc$Broad_category ==
          "TEC"
      )
      
    ),
    
    stringsAsFactors = FALSE
    
  )


special_category_summary$Percentage <-
  round(
    100 *
      special_category_summary$N_genes /
      nrow(gene_qc),
    2
  )


cat("\n")
cat("====================================================\n")
cat("SPECIAL CATEGORY SUMMARY\n")
cat("====================================================\n\n")

print(
  special_category_summary,
  row.names = FALSE
)


############################################################
# 20. Variance distribution summary
############################################################

variance_summary <-
  data.frame(
    
    Statistic = c(
      "Minimum",
      "1st percentile",
      "5th percentile",
      "10th percentile",
      "25th percentile",
      "Median",
      "75th percentile",
      "90th percentile",
      "95th percentile",
      "99th percentile",
      "Maximum",
      "Mean",
      "SD"
    ),
    
    Value = c(
      
      min(
        gene_variance,
        na.rm = TRUE
      ),
      
      quantile(
        gene_variance,
        0.01,
        na.rm = TRUE
      ),
      
      quantile(
        gene_variance,
        0.05,
        na.rm = TRUE
      ),
      
      quantile(
        gene_variance,
        0.10,
        na.rm = TRUE
      ),
      
      quantile(
        gene_variance,
        0.25,
        na.rm = TRUE
      ),
      
      median(
        gene_variance,
        na.rm = TRUE
      ),
      
      quantile(
        gene_variance,
        0.75,
        na.rm = TRUE
      ),
      
      quantile(
        gene_variance,
        0.90,
        na.rm = TRUE
      ),
      
      quantile(
        gene_variance,
        0.95,
        na.rm = TRUE
      ),
      
      quantile(
        gene_variance,
        0.99,
        na.rm = TRUE
      ),
      
      max(
        gene_variance,
        na.rm = TRUE
      ),
      
      mean(
        gene_variance,
        na.rm = TRUE
      ),
      
      sd(
        gene_variance,
        na.rm = TRUE
      )
      
    )
    
  )


cat("\n")
cat("====================================================\n")
cat("GENE VARIANCE DISTRIBUTION\n")
cat("====================================================\n\n")

print(
  variance_summary,
  row.names = FALSE
)


############################################################
# 21. MAD distribution summary
############################################################

mad_summary <-
  data.frame(
    
    Statistic = c(
      "Minimum",
      "1st percentile",
      "5th percentile",
      "10th percentile",
      "25th percentile",
      "Median",
      "75th percentile",
      "90th percentile",
      "95th percentile",
      "99th percentile",
      "Maximum",
      "Mean",
      "SD"
    ),
    
    Value = c(
      
      min(
        gene_mad,
        na.rm = TRUE
      ),
      
      quantile(
        gene_mad,
        0.01,
        na.rm = TRUE
      ),
      
      quantile(
        gene_mad,
        0.05,
        na.rm = TRUE
      ),
      
      quantile(
        gene_mad,
        0.10,
        na.rm = TRUE
      ),
      
      quantile(
        gene_mad,
        0.25,
        na.rm = TRUE
      ),
      
      median(
        gene_mad,
        na.rm = TRUE
      ),
      
      quantile(
        gene_mad,
        0.75,
        na.rm = TRUE
      ),
      
      quantile(
        gene_mad,
        0.90,
        na.rm = TRUE
      ),
      
      quantile(
        gene_mad,
        0.95,
        na.rm = TRUE
      ),
      
      quantile(
        gene_mad,
        0.99,
        na.rm = TRUE
      ),
      
      max(
        gene_mad,
        na.rm = TRUE
      ),
      
      mean(
        gene_mad,
        na.rm = TRUE
      ),
      
      sd(
        gene_mad,
        na.rm = TRUE
      )
      
    )
    
  )


cat("\n")
cat("====================================================\n")
cat("GENE MAD DISTRIBUTION\n")
cat("====================================================\n\n")

print(
  mad_summary,
  row.names = FALSE
)


############################################################
# 22. Zero-variance genes
############################################################

zero_variance_genes <-
  gene_qc[
    gene_qc$Variance == 0,
    ,
    drop = FALSE
  ]


cat("\n")
cat("====================================================\n")
cat("ZERO-VARIANCE GENES\n")
cat("====================================================\n\n")

cat(
  "Number of zero-variance genes:",
  nrow(zero_variance_genes),
  "\n\n"
)


############################################################
# 23. Very-low-MAD genes
############################################################

mad_quantiles <-
  quantile(
    gene_qc$MAD,
    probs = c(
      0.01,
      0.05,
      0.10,
      0.20,
      0.25
    ),
    na.rm = TRUE
  )


mad_threshold_summary <-
  data.frame(
    
    Percentile =
      names(mad_quantiles),
    
    MAD =
      as.numeric(mad_quantiles),
    
    N_genes_below_or_equal =
      sapply(
        mad_quantiles,
        function(x)
          sum(
            gene_qc$MAD <= x,
            na.rm = TRUE
          )
      )
    
  )


cat("\n")
cat("====================================================\n")
cat("LOW-MAD DISTRIBUTION\n")
cat("====================================================\n\n")

print(
  mad_threshold_summary,
  row.names = FALSE
)


############################################################
# 24. Top variable genes
############################################################

top_variable_genes <-
  gene_qc %>%
  
  arrange(
    desc(Variance)
  ) %>%
  
  slice_head(
    n = 100
  )


############################################################
# 25. Top MAD genes
############################################################

top_mad_genes <-
  gene_qc %>%
  
  arrange(
    desc(MAD)
  ) %>%
  
  slice_head(
    n = 100
  )


############################################################
# 26. Top variable genes by gene type
############################################################

top_variable_gene_type_summary <-
  top_variable_genes %>%
  
  count(
    gene_type,
    sort = TRUE,
    name = "N_genes"
  )


############################################################
# 27. Top MAD genes by gene type
############################################################

top_mad_gene_type_summary <-
  top_mad_genes %>%
  
  count(
    gene_type,
    sort = TRUE,
    name = "N_genes"
  )


############################################################
# 28. Important genes / features
############################################################

important_genes <- c(
  
  # Paper 1 candidate genes
  "AOC3",
  "MS4V7",
  "IL7R",
  
  # housekeeping
  "GAPDH",
  "ACTB",
  
  # mitochondrial
  "MT-ND1",
  "MT-ND2",
  "MT-ND3",
  "MT-ND4",
  "MT-CO1",
  "MT-CO2",
  "MT-CO3",
  "MT-ATP6",
  
  # high-expression RNA
  "MALAT1",
  "MT-RNR2",
  "RN7SL2",
  "RN7SK",
  
  # immune-related genes previously observed
  "S100A9",
  "IGKC",
  "IGHA1",
  
  # epithelial / lung-associated
  "SFTPA2",
  "SCGB3A1",
  "SLPI"
)


important_gene_annotation <-
  gene_qc[
    gene_qc$gene_name %in%
      important_genes,
    ,
    drop = FALSE
  ]


############################################################
# 29. Also search by Ensembl ID for previously identified
# extreme features
############################################################

important_ensembl <- c(
  
  "ENSG00000274012",
  "ENSG00000202198",
  "ENSG00000210082",
  "ENSG00000251562",
  "ENSG00000198888",
  "ENSG00000198899",
  "ENSG00000198938",
  "ENSG00000198763",
  "ENSG00000198804",
  "ENSG00000198840",
  "ENSG00000198712",
  "ENSG00000198886"
  
)


important_ensembl_annotation <-
  gene_qc[
    gene_qc$Ensembl_ID %in%
      important_ensembl,
    ,
    drop = FALSE
  ]


############################################################
# 30. Combined important-feature annotation
############################################################

important_features_combined <-
  unique(
    rbind(
      important_gene_annotation,
      important_ensembl_annotation
    )
  )


############################################################
# 31. Check Paper 1 candidate genes
############################################################

candidate_gene_check <-
  gene_qc[
    gene_qc$gene_name %in%
      c(
        "AOC3",
        "MS4V7",
        "IL7R"
      ),
    ,
    drop = FALSE
  ]


cat("\n")
cat("====================================================\n")
cat("PAPER 1 CANDIDATE GENE QC\n")
cat("====================================================\n\n")

print(
  candidate_gene_check,
  row.names = FALSE
)


############################################################
# 32. Save full gene QC table
############################################################

write.table(
  gene_qc,
  file =
    "result_dir/Workflow-04-2C-2_GENCODE_v22_gene_QC_full.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 33. Save gene-type summary
############################################################

write.table(
  gene_type_summary,
  file =
    "result_dir/Workflow-04-2C-2_gene_type_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 34. Save broad-category summary
############################################################

write.table(
  broad_category_summary,
  file =
    "result_dir/Workflow-04-2C-2_broad_gene_category_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 35. Save special-category summary
############################################################

write.table(
  special_category_summary,
  file =
    "result_dir/Workflow-04-2C-2_special_category_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 36. Save variance summary
############################################################

write.table(
  variance_summary,
  file =
    "result_dir/Workflow-04-2C-2_variance_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 37. Save MAD summary
############################################################

write.table(
  mad_summary,
  file =
    "result_dir/Workflow-04-2C-2_MAD_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 38. Save low-MAD summary
############################################################

write.table(
  mad_threshold_summary,
  file =
    "result_dir/Workflow-04-2C-2_low_MAD_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 39. Save zero-variance genes
############################################################

write.table(
  zero_variance_genes,
  file =
    "result_dir/Workflow-04-2C-2_zero_variance_genes.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 40. Save top variable genes
############################################################

write.table(
  top_variable_genes,
  file =
    "result_dir/Workflow-04-2C-2_top100_variable_genes.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 41. Save top MAD genes
############################################################

write.table(
  top_mad_genes,
  file =
    "result_dir/Workflow-04-2C-2_top100_MAD_genes.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 42. Save top variable gene-type summary
############################################################

write.table(
  top_variable_gene_type_summary,
  file =
    "result_dir/Workflow-04-2C-2_top100_variable_gene_type_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 43. Save top MAD gene-type summary
############################################################

write.table(
  top_mad_gene_type_summary,
  file =
    "result_dir/Workflow-04-2C-2_top100_MAD_gene_type_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 44. Save important feature annotation
############################################################

write.table(
  important_features_combined,
  file =
    "result_dir/Workflow-04-2C-2_important_feature_annotation.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 45. Save candidate gene QC
############################################################

write.table(
  candidate_gene_check,
  file =
    "result_dir/Workflow-04-2C-2_candidate_gene_QC.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 46. Create compact overall summary
############################################################

overall_summary <-
  data.frame(
    
    Metric = c(
      
      "Expression genes",
      "Expression samples",
      "GENCODE v22 matched genes",
      "GENCODE v22 match rate (%)",
      
      "Protein-coding genes",
      "Non-protein-coding genes",
      "Mitochondrial genes",
      "rRNA-related genes",
      
      "Long non-coding RNA",
      "Small RNA",
      "Pseudogenes",
      "TEC",
      
      "Zero-variance genes"
      
    ),
    
    Value = c(
      
      nrow(expr_mat),
      
      ncol(expr_mat),
      
      n_matched,
      
      round(
        match_rate,
        4
      ),
      
      sum(
        gene_qc$gene_type ==
          "protein_coding"
      ),
      
      sum(
        gene_qc$Non_protein_coding
      ),
      
      sum(
        gene_qc$Mitochondrial
      ),
      
      sum(
        gene_qc$rRNA_related
      ),
      
      sum(
        gene_qc$Broad_category ==
          "Long non-coding RNA"
      ),
      
      sum(
        gene_qc$Broad_category ==
          "Small RNA"
      ),
      
      sum(
        gene_qc$Broad_category ==
          "Pseudogene"
      ),
      
      sum(
        gene_qc$Broad_category ==
          "TEC"
      ),
      
      nrow(
        zero_variance_genes
      )
      
    )
    
  )


############################################################
# 47. Save overall summary
############################################################

write.table(
  overall_summary,
  file =
    "result_dir/Workflow-04-2C-2_overall_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 48. Final console summary
############################################################

cat("\n")
cat("====================================================\n")
cat("Workflow-04-2C-2 COMPLETED\n")
cat("====================================================\n\n")

cat(
  "Genes analyzed:",
  nrow(expr_mat),
  "\n"
)

cat(
  "Samples analyzed:",
  ncol(expr_mat),
  "\n"
)

cat(
  "GENCODE v22 match rate:",
  round(match_rate, 4),
  "%\n"
)

cat(
  "Protein-coding genes:",
  sum(
    gene_qc$gene_type ==
      "protein_coding"
  ),
  "\n"
)

cat(
  "Non-protein-coding genes:",
  sum(
    gene_qc$Non_protein_coding
  ),
  "\n"
)

cat(
  "Mitochondrial genes:",
  sum(
    gene_qc$Mitochondrial
  ),
  "\n"
)

cat(
  "rRNA-related genes:",
  sum(
    gene_qc$rRNA_related
  ),
  "\n"
)

cat(
  "Zero-variance genes:",
  nrow(
    zero_variance_genes
  ),
  "\n\n"
)

cat(
  "IMPORTANT:\n"
)

cat(
  "No genes were removed in Workflow-04-2C-2.\n"
)

cat(
  "No samples were removed in Workflow-04-2C-2.\n"
)

cat(
  "No WGCNA network was constructed.\n\n"
)

cat(
  "All results saved to:\n",
  getwd(),
  "\n"
)

cat(
  "\nNext step should be determined only after QC results are reviewed.\n"
)

############################################################
# END
############################################################ 


############################################################
# Workflow-04-2C-3
# MAD=0 gene composition and candidate-universe comparison
#
# Purpose:
# 1. Characterize the 7,671 MAD=0 genes
# 2. Compare all 39,280 genes vs MAD>0 genes
# 3. Evaluate composition and Paper 1 candidate genes
#
# No genes are removed from the original files.
# No WGCNA is performed in this workflow.
############################################################


############################################################
# 0. Setup
############################################################

rm(list = ls())

options(stringsAsFactors = FALSE)

project_dir <- "input_dir/GitHub2_Paper1_tumor-only WGCNA"

setwd(project_dir)


############################################################
# 1. Input files
############################################################

qc_file <- file.path(
  project_dir,
  "Workflow-04-2C-2_GENCODE_v22_gene_QC_full.txt"
)

candidate_file <- file.path(
  project_dir,
  "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt"
)


############################################################
# 2. Read GENCODE v22 gene QC table
############################################################

gene_qc <- read.table(
  qc_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n==============================\n")
cat("Part A: Input QC table\n")
cat("==============================\n")

cat("Rows:", nrow(gene_qc), "\n")
cat("Columns:", ncol(gene_qc), "\n")

cat("\nColumn names:\n")
print(colnames(gene_qc))


############################################################
# 3. Basic integrity checks
############################################################

required_cols <- c(
  "Ensembl_ID",
  "gene_type",
  "Variance",
  "MAD"
)

missing_cols <- setdiff(
  required_cols,
  colnames(gene_qc)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

if (anyDuplicated(gene_qc$Ensembl_ID) > 0) {
  stop("Duplicated Ensembl IDs detected.")
}

if (anyNA(gene_qc$MAD)) {
  stop("NA values detected in MAD.")
}

if (anyNA(gene_qc$Variance)) {
  stop("NA values detected in Variance.")
}


############################################################
# 4. Define MAD groups
############################################################

gene_qc$MAD_group <- ifelse(
  gene_qc$MAD == 0,
  "MAD = 0",
  "MAD > 0"
)

gene_qc$MAD_group <- factor(
  gene_qc$MAD_group,
  levels = c("MAD = 0", "MAD > 0")
)


############################################################
# 5. Basic MAD composition
############################################################

mad_counts <- table(gene_qc$MAD_group)

mad_summary <- data.frame(
  MAD_group = names(mad_counts),
  Gene_N = as.numeric(mad_counts),
  Percentage = round(
    as.numeric(mad_counts) / nrow(gene_qc) * 100,
    2
  )
)

cat("\n==============================\n")
cat("Part B: MAD group composition\n")
cat("==============================\n")

print(mad_summary)


############################################################
# 6. Extract MAD=0 and MAD>0 gene universes
############################################################

genes_mad0 <- gene_qc[
  gene_qc$MAD == 0,
]

genes_mad_gt0 <- gene_qc[
  gene_qc$MAD > 0,
]

cat("\nMAD = 0 genes:", nrow(genes_mad0), "\n")
cat("MAD > 0 genes:", nrow(genes_mad_gt0), "\n")
cat("Total genes:", nrow(gene_qc), "\n")


############################################################
# 7. Verify expected counts
############################################################

if (
  nrow(genes_mad0) + nrow(genes_mad_gt0) != nrow(gene_qc)
) {
  stop("MAD groups do not sum to total gene number.")
}


############################################################
# 8. Gene-type composition
############################################################

gene_type_all <- as.data.frame(
  table(gene_qc$gene_type),
  stringsAsFactors = FALSE
)

colnames(gene_type_all) <- c(
  "gene_type",
  "Gene_N"
)

gene_type_all$Percentage <- round(
  gene_type_all$Gene_N /
    sum(gene_type_all$Gene_N) * 100,
  2
)


gene_type_mad0 <- as.data.frame(
  table(genes_mad0$gene_type),
  stringsAsFactors = FALSE
)

colnames(gene_type_mad0) <- c(
  "gene_type",
  "MAD0_N"
)


gene_type_madgt0 <- as.data.frame(
  table(genes_mad_gt0$gene_type),
  stringsAsFactors = FALSE
)

colnames(gene_type_madgt0) <- c(
  "gene_type",
  "MADgt0_N"
)


############################################################
# 9. Combined gene-type comparison table
############################################################

gene_type_comparison <- merge(
  gene_type_all,
  gene_type_mad0,
  by = "gene_type",
  all = TRUE
)

gene_type_comparison <- merge(
  gene_type_comparison,
  gene_type_madgt0,
  by = "gene_type",
  all = TRUE
)

gene_type_comparison$MAD0_N[
  is.na(gene_type_comparison$MAD0_N)
] <- 0

gene_type_comparison$MADgt0_N[
  is.na(gene_type_comparison$MADgt0_N)
] <- 0


############################################################
# 10. Calculate within-group percentages
############################################################

total_mad0 <- nrow(genes_mad0)
total_madgt0 <- nrow(genes_mad_gt0)

gene_type_comparison$MAD0_Percentage <- round(
  gene_type_comparison$MAD0_N /
    total_mad0 * 100,
  2
)

gene_type_comparison$MADgt0_Percentage <- round(
  gene_type_comparison$MADgt0_N /
    total_madgt0 * 100,
  2
)


############################################################
# 11. Calculate MAD=0 enrichment ratio
#
# Ratio > 1:
# gene type is relatively enriched in MAD=0 genes
#
# Ratio < 1:
# gene type is relatively depleted in MAD=0 genes
############################################################

gene_type_comparison$MAD0_Enrichment_Ratio <- round(
  gene_type_comparison$MAD0_Percentage /
    gene_type_comparison$Percentage,
  3
)


############################################################
# 12. Calculate fraction of each gene type that is MAD=0
############################################################

gene_type_comparison$Fraction_of_gene_type_MAD0 <- round(
  gene_type_comparison$MAD0_N /
    gene_type_comparison$Gene_N * 100,
  2
)


############################################################
# 13. Sort by MAD=0 gene number
############################################################

gene_type_comparison <- gene_type_comparison[
  order(
    -gene_type_comparison$MAD0_N
  ),
]


############################################################
# 14. Print gene-type comparison
############################################################

cat("\n==============================\n")
cat("Part C: Gene-type composition\n")
cat("==============================\n")

print(
  gene_type_comparison,
  row.names = FALSE
)


############################################################
# 15. Broad category composition
#
# Broad_category should already exist in QC table.
############################################################

if ("Broad_category" %in% colnames(gene_qc)) {
  
  broad_all <- as.data.frame(
    table(gene_qc$Broad_category),
    stringsAsFactors = FALSE
  )
  
  colnames(broad_all) <- c(
    "Broad_category",
    "Gene_N"
  )
  
  broad_mad0 <- as.data.frame(
    table(genes_mad0$Broad_category),
    stringsAsFactors = FALSE
  )
  
  colnames(broad_mad0) <- c(
    "Broad_category",
    "MAD0_N"
  )
  
  broad_madgt0 <- as.data.frame(
    table(genes_mad_gt0$Broad_category),
    stringsAsFactors = FALSE
  )
  
  colnames(broad_madgt0) <- c(
    "Broad_category",
    "MADgt0_N"
  )
  
  broad_comparison <- merge(
    broad_all,
    broad_mad0,
    by = "Broad_category",
    all = TRUE
  )
  
  broad_comparison <- merge(
    broad_comparison,
    broad_madgt0,
    by = "Broad_category",
    all = TRUE
  )
  
  broad_comparison$MAD0_N[
    is.na(broad_comparison$MAD0_N)
  ] <- 0
  
  broad_comparison$MADgt0_N[
    is.na(broad_comparison$MADgt0_N)
  ] <- 0
  
  broad_comparison$All_Percentage <- round(
    broad_comparison$Gene_N /
      sum(broad_comparison$Gene_N) * 100,
    2
  )
  
  broad_comparison$MAD0_Percentage <- round(
    broad_comparison$MAD0_N /
      total_mad0 * 100,
    2
  )
  
  broad_comparison$MADgt0_Percentage <- round(
    broad_comparison$MADgt0_N /
      total_madgt0 * 100,
    2
  )
  
  broad_comparison$MAD0_Enrichment_Ratio <- round(
    broad_comparison$MAD0_Percentage /
      broad_comparison$All_Percentage,
    3
  )
  
  broad_comparison <- broad_comparison[
    order(-broad_comparison$MAD0_N),
  ]
  
} else {
  
  broad_comparison <- NULL
  
  cat(
    "\nWARNING: Broad_category column not found.\n"
  )
}


############################################################
# 16. Special-category comparison
############################################################

special_categories <- data.frame(
  Category = c(
    "Protein-coding",
    "Non-protein-coding",
    "Mitochondrial",
    "rRNA-related",
    "Long non-coding RNA",
    "Small RNA",
    "Pseudogene",
    "TEC"
  ),
  All_N = NA_integer_,
  MAD0_N = NA_integer_,
  MADgt0_N = NA_integer_
)


# Protein coding
special_categories$All_N[1] <-
  sum(
    gene_qc$Broad_category ==
      "Protein-coding"
  )

special_categories$MAD0_N[1] <-
  sum(
    genes_mad0$Broad_category ==
      "Protein-coding"
  )

special_categories$MADgt0_N[1] <-
  sum(
    genes_mad_gt0$Broad_category ==
      "Protein-coding"
  )


# Non-protein-coding
special_categories$All_N[2] <-
  sum(
    gene_qc$Broad_category !=
      "Protein-coding"
  )

special_categories$MAD0_N[2] <-
  sum(
    genes_mad0$Broad_category !=
      "Protein-coding"
  )

special_categories$MADgt0_N[2] <-
  sum(
    genes_mad_gt0$Broad_category !=
      "Protein-coding"
  )


# Mitochondrial
special_categories$All_N[3] <-
  sum(
    gene_qc$Broad_category ==
      "Mitochondrial RNA"
  )

special_categories$MAD0_N[3] <-
  sum(
    genes_mad0$Broad_category ==
      "Mitochondrial RNA"
  )

special_categories$MADgt0_N[3] <-
  sum(
    genes_mad_gt0$Broad_category ==
      "Mitochondrial RNA"
  )


# rRNA
special_categories$All_N[4] <-
  sum(
    gene_qc$Broad_category ==
      "Ribosomal RNA"
  )

special_categories$MAD0_N[4] <-
  sum(
    genes_mad0$Broad_category ==
      "Ribosomal RNA"
  )

special_categories$MADgt0_N[4] <-
  sum(
    genes_mad_gt0$Broad_category ==
      "Ribosomal RNA"
  )


# Long non-coding RNA
special_categories$All_N[5] <-
  sum(
    gene_qc$Broad_category ==
      "Long non-coding RNA"
  )

special_categories$MAD0_N[5] <-
  sum(
    genes_mad0$Broad_category ==
      "Long non-coding RNA"
  )

special_categories$MADgt0_N[5] <-
  sum(
    genes_mad_gt0$Broad_category ==
      "Long non-coding RNA"
  )


# Small RNA
special_categories$All_N[6] <-
  sum(
    gene_qc$Broad_category ==
      "Small RNA"
  )

special_categories$MAD0_N[6] <-
  sum(
    genes_mad0$Broad_category ==
      "Small RNA"
  )

special_categories$MADgt0_N[6] <-
  sum(
    genes_mad_gt0$Broad_category ==
      "Small RNA"
  )


# Pseudogene
special_categories$All_N[7] <-
  sum(
    gene_qc$Broad_category ==
      "Pseudogene"
  )

special_categories$MAD0_N[7] <-
  sum(
    genes_mad0$Broad_category ==
      "Pseudogene"
  )

special_categories$MADgt0_N[7] <-
  sum(
    genes_mad_gt0$Broad_category ==
      "Pseudogene"
  )


# TEC
special_categories$All_N[8] <-
  sum(
    gene_qc$Broad_category ==
      "TEC"
  )

special_categories$MAD0_N[8] <-
  sum(
    genes_mad0$Broad_category ==
      "TEC"
  )

special_categories$MADgt0_N[8] <-
  sum(
    genes_mad_gt0$Broad_category ==
      "TEC"
  )


############################################################
# 17. Percentages
############################################################

special_categories$All_Percentage <- round(
  special_categories$All_N /
    nrow(gene_qc) * 100,
  2
)

special_categories$MAD0_Percentage <- round(
  special_categories$MAD0_N /
    total_mad0 * 100,
  2
)

special_categories$MADgt0_Percentage <- round(
  special_categories$MADgt0_N /
    total_madgt0 * 100,
  2
)


############################################################
# 18. Paper 1 candidate genes
############################################################

paper1_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


############################################################
# 19. Candidate gene check
############################################################

candidate_check <- gene_qc[
  gene_qc$gene_symbol %in% paper1_genes,
]

cat("\n==============================\n")
cat("Part D: Paper 1 candidate genes\n")
cat("==============================\n")

print(
  candidate_check[
    ,
    intersect(
      c(
        "Ensembl_ID",
        "gene_symbol",
        "gene_type",
        "Broad_category",
        "Variance",
        "MAD",
        "MAD_group"
      ),
      colnames(candidate_check)
    ),
    drop = FALSE
  ],
  row.names = FALSE
)


############################################################
# 20. Candidate genes missing from QC table
############################################################

candidate_found <- unique(
  candidate_check$gene_symbol
)

candidate_missing <- setdiff(
  paper1_genes,
  candidate_found
)

cat("\nPaper 1 candidate genes found:\n")
print(candidate_found)

cat("\nPaper 1 candidate genes missing:\n")
print(candidate_missing)


###21. candidate gene check AOC3–MS4A7–IL7R 

colnames(gene_qc) 
gene_qc[
  gene_qc$Ensembl_ID %in% c(
    "ENSG00000131471",
    "ENSG00000168685"
  ),
  ,
  drop = FALSE
]

grep(
  "MS4A7",
  gene_qc$gene_name,
  value = TRUE,
  ignore.case = TRUE
)  

############################################################
# Corrected Paper 1 candidate gene check
############################################################

paper1_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

candidate_check <- gene_qc[
  gene_qc$gene_name %in% paper1_genes,
  ,
  drop = FALSE
]

candidate_check <- candidate_check[
  match(
    paper1_genes,
    candidate_check$gene_name
  ),
  ,
  drop = FALSE
]

print(
  candidate_check[
    ,
    c(
      "Ensembl_ID",
      "gene_name",
      "gene_type",
      "Mean",
      "Median",
      "SD",
      "Variance",
      "MAD",
      "Detected_N",
      "Detected_fraction",
      "Broad_category",
      "MAD_group"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)  

############################################################
# Final MAD status of Paper 1 core genes
############################################################

candidate_status <- data.frame(
  Gene = paper1_genes,
  stringsAsFactors = FALSE
)

candidate_status$Found <- paper1_genes %in% gene_qc$gene_name

candidate_status$MAD_group <- NA_character_

idx <- match(
  candidate_status$Gene,
  gene_qc$gene_name
)

candidate_status$MAD_group <- as.character(
  gene_qc$MAD_group[idx]
)

candidate_status$Variance <- gene_qc$Variance[idx]

candidate_status$MAD <- gene_qc$MAD[idx]

candidate_status$Detected_N <- gene_qc$Detected_N[idx]

candidate_status$Detected_fraction <-
  gene_qc$Detected_fraction[idx]

print(
  candidate_status,
  row.names = FALSE
)


############################################################
# 22. Variance comparison between MAD groups
############################################################

variance_mad_summary <- data.frame(
  MAD_group = c(
    "MAD = 0",
    "MAD > 0"
  ),
  Gene_N = c(
    nrow(genes_mad0),
    nrow(genes_mad_gt0)
  ),
  Variance_median = c(
    median(
      genes_mad0$Variance,
      na.rm = TRUE
    ),
    median(
      genes_mad_gt0$Variance,
      na.rm = TRUE
    )
  ),
  Variance_mean = c(
    mean(
      genes_mad0$Variance,
      na.rm = TRUE
    ),
    mean(
      genes_mad_gt0$Variance,
      na.rm = TRUE
    )
  ),
  Variance_min = c(
    min(
      genes_mad0$Variance,
      na.rm = TRUE
    ),
    min(
      genes_mad_gt0$Variance,
      na.rm = TRUE
    )
  ),
  Variance_max = c(
    max(
      genes_mad0$Variance,
      na.rm = TRUE
    ),
    max(
      genes_mad_gt0$Variance,
      na.rm = TRUE
    )
  )
)


############################################################
# 23. Compare candidate universes
############################################################

universe_comparison <- data.frame(
  Universe = c(
    "All 39,280 genes",
    "MAD > 0 genes"
  ),
  Gene_N = c(
    nrow(gene_qc),
    nrow(genes_mad_gt0)
  ),
  Protein_coding_N = c(
    sum(
      gene_qc$Broad_category ==
        "Protein-coding"
    ),
    sum(
      genes_mad_gt0$Broad_category ==
        "Protein-coding"
    )
  ),
  Protein_coding_percentage = c(
    round(
      sum(
        gene_qc$Broad_category ==
          "Protein-coding"
      ) /
        nrow(gene_qc) * 100,
      2
    ),
    round(
      sum(
        genes_mad_gt0$Broad_category ==
          "Protein-coding"
      ) /
        nrow(genes_mad_gt0) * 100,
      2
    )
  ),
  Non_protein_coding_percentage = c(
    round(
      sum(
        gene_qc$Broad_category !=
          "Protein-coding"
      ) /
        nrow(gene_qc) * 100,
      2
    ),
    round(
      sum(
        genes_mad_gt0$Broad_category !=
          "Protein-coding"
      ) /
        nrow(genes_mad_gt0) * 100,
      2
    )
  ),
  Mitochondrial_percentage = c(
    round(
      sum(
        gene_qc$Broad_category ==
          "Mitochondrial RNA"
      ) /
        nrow(gene_qc) * 100,
      2
    ),
    round(
      sum(
        genes_mad_gt0$Broad_category ==
          "Mitochondrial RNA"
      ) /
        nrow(genes_mad_gt0) * 100,
      2
    )
  ),
  rRNA_percentage = c(
    round(
      sum(
        gene_qc$Broad_category ==
          "Ribosomal RNA"
      ) /
        nrow(gene_qc) * 100,
      2
    ),
    round(
      sum(
        genes_mad_gt0$Broad_category ==
          "Ribosomal RNA"
      ) /
        nrow(genes_mad_gt0) * 100,
      2
    )
  )
)


############################################################
# 24. Save outputs
############################################################

write.table(
  mad_summary,
  file.path(
    project_dir,
    "Workflow-04-2C-3_MAD_group_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  gene_type_comparison,
  file.path(
    project_dir,
    "Workflow-04-2C-3_gene_type_composition_MAD_comparison.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


if (!is.null(broad_comparison)) {
  
  write.table(
    broad_comparison,
    file.path(
      project_dir,
      "Workflow-04-2C-3_broad_category_MAD_comparison.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}


write.table(
  special_categories,
  file.path(
    project_dir,
    "Workflow-04-2C-3_special_category_MAD_comparison.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  candidate_check,
  file.path(
    project_dir,
    "Workflow-04-2C-3_Paper1_candidate_MAD_check.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  candidate_status,
  file.path(
    project_dir,
    "Workflow-04-2C-3_Paper1_candidate_status.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  variance_mad_summary,
  file.path(
    project_dir,
    "Workflow-04-2C-3_variance_by_MAD_group.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  universe_comparison,
  file.path(
    project_dir,
    "Workflow-04-2C-3_candidate_universe_comparison.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 25. Save complete gene lists
############################################################

write.table(
  genes_mad0,
  file.path(
    project_dir,
    "Workflow-04-2C-3_MAD_zero_genes_full.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  genes_mad_gt0,
  file.path(
    project_dir,
    "Workflow-04-2C-3_MAD_positive_genes_full.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 26. Overall summary
############################################################

summary_lines <- c(
  "Workflow-04-2C-3 summary",
  "=========================",
  "",
  paste("Total genes:", nrow(gene_qc)),
  paste("MAD = 0 genes:", nrow(genes_mad0)),
  paste("MAD > 0 genes:", nrow(genes_mad_gt0)),
  paste(
    "MAD = 0 percentage:",
    round(nrow(genes_mad0) /
            nrow(gene_qc) * 100, 2)
  ),
  paste(
    "MAD > 0 percentage:",
    round(nrow(genes_mad_gt0) /
            nrow(gene_qc) * 100, 2)
  ),
  "",
  "Paper 1 candidate genes:",
  paste(paper1_genes, collapse = ", "),
  paste(
    "Missing from QC table:",
    ifelse(
      length(candidate_missing) == 0,
      "None",
      paste(candidate_missing, collapse = ", ")
    )
  ),
  "",
  paste(
    "All genes universe:",
    nrow(gene_qc)
  ),
  paste(
    "MAD > 0 universe:",
    nrow(genes_mad_gt0)
  )
)

writeLines(
  summary_lines,
  file.path(
    project_dir,
    "Workflow-04-2C-3_summary.txt"
  )
)


############################################################
# 27. Final console summary
############################################################

cat("\n\n")
cat("====================================================\n")
cat("Workflow-04-2C-3 COMPLETED\n")
cat("====================================================\n")

cat("\nTotal genes:", nrow(gene_qc), "\n")
cat("MAD = 0 genes:", nrow(genes_mad0), "\n")
cat("MAD > 0 genes:", nrow(genes_mad_gt0), "\n")

cat(
  "\nPaper 1 genes:",
  paste(paper1_genes, collapse = ", "),
  "\n"
)

cat("\nCandidate status:\n")
print(candidate_status)

cat("\nCandidate universe comparison:\n")
print(universe_comparison)

cat("\nOutput files saved to:\n")
cat(project_dir, "\n")

cat("\n====================================================\n")  




############################################################
# Workflow-04-2D
# WGCNA pre-analysis:
# Final input matrix construction and sample-level QC
#
# Input:
# Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt
# Workflow-04-2C-2_GENCODE_v22_gene_QC_full.txt
# Workflow-04-2A_clinical_aligned.txt
#
# No WGCNA network construction is performed in this workflow.
############################################################


############################
# 0. Basic settings
############################

rm(list = ls())

options(stringsAsFactors = FALSE)

set.seed(20260909)


############################
# 1. Working directory
############################

setwd("input_dir/GitHub2_Paper1_tumor-only WGCNA")


############################
# 2. Load packages
############################

required_packages <- c(
  "WGCNA",
  "ggplot2"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "Package '", pkg,
        "' is not installed. Please install it before running this workflow."
      )
    )
  }
}

library(WGCNA)
library(ggplot2)

options(stringsAsFactors = FALSE)

# Allow WGCNA to use multiple cores if available
allowWGCNAThreads()


############################
# 3. Define input files
############################

expr_file <- 
  "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt"

gene_qc_file <- 
  "Workflow-04-2C-2_GENCODE_v22_gene_QC_full.txt"

clinical_file <- 
  "Workflow-04-2A_clinical_aligned.txt"


############################
# 4. Read expression matrix
############################

expr_df <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n========== Expression input ==========\n")

cat(
  "Rows:",
  nrow(expr_df),
  "\n"
)

cat(
  "Columns:",
  ncol(expr_df),
  "\n"
)

cat(
  "Number of genes:",
  nrow(expr_df),
  "\n"
)

cat(
  "Number of samples:",
  ncol(expr_df) - 1,
  "\n"
)


############################
# 5. Basic expression structure QC
############################

if (!"Tag" %in% colnames(expr_df)) {
  
  stop(
    "The first column 'Tag' was not found."
  )
}

gene_ids <- expr_df$Tag

sample_ids <- colnames(expr_df)[-1]


############################
# 6. Check gene IDs
############################

cat("\n========== Gene ID QC ==========\n")

cat(
  "Missing gene IDs:",
  sum(is.na(gene_ids)),
  "\n"
)

cat(
  "Empty gene IDs:",
  sum(gene_ids == ""),
  "\n"
)

cat(
  "Duplicated gene IDs:",
  sum(duplicated(gene_ids)),
  "\n"
)


if (
  any(is.na(gene_ids)) ||
  any(gene_ids == "") ||
  any(duplicated(gene_ids))
) {
  
  stop(
    "Gene ID QC failed. Please resolve duplicated/missing gene IDs."
  )
}


############################
# 7. Check sample IDs
############################

cat("\n========== Sample ID QC ==========\n")

cat(
  "Number of samples:",
  length(sample_ids),
  "\n"
)

cat(
  "Duplicated sample IDs:",
  sum(duplicated(sample_ids)),
  "\n"
)

if (any(duplicated(sample_ids))) {
  
  stop(
    "Duplicated sample IDs detected."
  )
}


############################
# 8. Convert expression data
############################

expr_mat <- as.matrix(
  expr_df[, -1, drop = FALSE]
)

rownames(expr_mat) <- gene_ids


############################
# 9. Matrix numeric QC
############################

cat("\n========== Numeric matrix QC ==========\n")

cat(
  "Storage mode:",
  storage.mode(expr_mat),
  "\n"
)

cat(
  "NA:",
  sum(is.na(expr_mat)),
  "\n"
)

cat(
  "NaN:",
  sum(is.nan(expr_mat)),
  "\n"
)

cat(
  "Inf:",
  sum(is.infinite(expr_mat)),
  "\n"
)

cat(
  "Negative values:",
  sum(expr_mat < 0, na.rm = TRUE),
  "\n"
)


if (
  any(is.na(expr_mat)) ||
  any(is.nan(expr_mat)) ||
  any(is.infinite(expr_mat))
) {
  
  stop(
    "Expression matrix contains NA/NaN/Inf values."
  )
}

if (any(expr_mat < 0)) {
  
  stop(
    "Negative expression values detected."
  )
}


############################
# 10. Read gene QC annotation
############################

gene_qc <- read.table(
  gene_qc_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n========== Gene annotation ==========\n")

cat(
  "Annotated genes:",
  nrow(gene_qc),
  "\n"
)

cat(
  "Gene QC columns:\n"
)

print(
  colnames(gene_qc)
)


############################################################
# Verify numeric MAD filtering
############################################################

cat("\n========== Verify MAD filtering ==========\n")

cat(
  "Total genes:",
  nrow(gene_qc),
  "\n"
)

cat(
  "MAD = 0:",
  sum(gene_qc$MAD == 0, na.rm = TRUE),
  "\n"
)

cat(
  "MAD > 0:",
  sum(gene_qc$MAD > 0, na.rm = TRUE),
  "\n"
)

cat(
  "MAD missing:",
  sum(is.na(gene_qc$MAD)),
  "\n"
) 
############################################################
# Correct MAD > 0 filtering
############################################################

mad_genes <- gene_qc$Ensembl_ID[
  !is.na(gene_qc$MAD) &
    gene_qc$MAD > 0
]

mad_genes <- unique(mad_genes)

cat(
  "MAD > 0 genes:",
  length(mad_genes),
  "\n"
) 

############################################################
# Check MAD genes against expression matrix
############################################################

mad_genes_in_matrix <- intersect(
  mad_genes,
  rownames(expr_mat)
)

cat(
  "MAD > 0 genes found in expression matrix:",
  length(mad_genes_in_matrix),
  "\n"
)

cat(
  "MAD > 0 genes missing from expression matrix:",
  length(
    setdiff(
      mad_genes,
      rownames(expr_mat)
    )
  ),
  "\n"
)

############################################################
# Construct preliminary WGCNA candidate matrix
############################################################

expr_mad <- expr_mat[
  mad_genes_in_matrix,
  ,
  drop = FALSE
]

cat(
  "\n========== Preliminary WGCNA candidate matrix ==========\n"
)

cat(
  "Genes:",
  nrow(expr_mad),
  "\n"
)

cat(
  "Samples:",
  ncol(expr_mad),
  "\n"
)

cat(
  "Dimension:",
  paste(
    dim(expr_mad),
    collapse = " × "
  ),
  "\n"
)  

############################################################
# Confirm Paper 1 candidate genes
############################################################

candidate_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

candidate_check <- gene_qc[
  gene_qc$gene_name %in% candidate_genes,
  c(
    "Ensembl_ID",
    "gene_name",
    "gene_type",
    "Variance",
    "MAD",
    "Detected_N",
    "Detected_fraction"
  )
]

print(candidate_check)

############################################################
# Final matrix integrity check
############################################################

cat("\n========== Final matrix integrity ==========\n")

cat(
  "Genes:",
  nrow(expr_mad),
  "\n"
)

cat(
  "Samples:",
  ncol(expr_mad),
  "\n"
)

cat(
  "NA:",
  sum(is.na(expr_mad)),
  "\n"
)

cat(
  "NaN:",
  sum(is.nan(expr_mad)),
  "\n"
)

cat(
  "Inf:",
  sum(is.infinite(expr_mad)),
  "\n"
)

cat(
  "Duplicated genes:",
  sum(duplicated(rownames(expr_mad))),
  "\n"
)

cat(
  "Duplicated samples:",
  sum(duplicated(colnames(expr_mad))),
  "\n"
) 

############################################################
# 11. Identify MAD > 0 genes
############################################################

cat("\n========== MAD filtering ==========\n")

# Use numeric MAD directly.
# Do NOT rely on a text-based MAD_group column.

if (!"MAD" %in% colnames(gene_qc)) {
  
  stop(
    "Column 'MAD' not found in gene QC table."
  )
}

mad_genes <- gene_qc$Ensembl_ID[
  !is.na(gene_qc$MAD) &
    gene_qc$MAD > 0
]

mad_genes <- unique(mad_genes)

cat(
  "Total genes:",
  nrow(gene_qc),
  "\n"
)

cat(
  "MAD = 0 genes:",
  sum(
    gene_qc$MAD == 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "MAD > 0 genes:",
  length(mad_genes),
  "\n"
)

cat(
  "MAD missing:",
  sum(
    is.na(gene_qc$MAD)
  ),
  "\n"
)


############################################################
# 12. Check MAD genes against expression matrix
############################################################

mad_genes_in_matrix <- intersect(
  mad_genes,
  rownames(expr_mat)
)

cat(
  "MAD > 0 genes found in expression matrix:",
  length(mad_genes_in_matrix),
  "\n"
)

cat(
  "MAD > 0 genes missing from expression matrix:",
  length(
    setdiff(
      mad_genes,
      rownames(expr_mat)
    )
  ),
  "\n"
)


if (
  length(mad_genes_in_matrix) !=
  length(mad_genes)
) {
  
  stop(
    "Some MAD > 0 genes are missing from the expression matrix."
  )
}


############################################################
# 13. Construct preliminary WGCNA candidate matrix
############################################################

expr_mad <- expr_mat[
  mad_genes_in_matrix,
  ,
  drop = FALSE
]

cat(
  "\n========== Preliminary WGCNA candidate matrix ==========\n"
)

cat(
  "Genes:",
  nrow(expr_mad),
  "\n"
)

cat(
  "Samples:",
  ncol(expr_mad),
  "\n"
)

cat(
  "Dimension:",
  paste(
    dim(expr_mad),
    collapse = " × "
  ),
  "\n"
)


if (
  nrow(expr_mad) != length(mad_genes) ||
  ncol(expr_mad) != length(sample_ids)
) {
  
  stop(
    "Unexpected dimensions in expr_mad."
  )
} 


############################
# 14. Save preliminary matrix
############################

write.table(
  cbind(
    Tag = rownames(expr_mad),
    as.data.frame(
      expr_mad,
      check.names = FALSE
    )
  ),
  file =
    "Workflow-04-2D_WGCNA_candidate_MADgt0.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# Workflow-04-2D
# WGCNA package/function verification
############################################################

cat("\n========== WGCNA package verification ==========\n")

cat(
  "WGCNA installed:",
  "WGCNA" %in% rownames(installed.packages()),
  "\n"
)

cat(
  "WGCNA version:"
)

print(
  packageVersion("WGCNA")
)

cat(
  "\nLoading WGCNA...\n"
)

library(WGCNA)

cat(
  "\nWGCNA loaded:",
  "WGCNA" %in% loadedNamespaces(),
  "\n"
)

cat(
  "goodSamplesGenes exists:",
  exists(
    "goodSamplesGenes",
    where = asNamespace("WGCNA"),
    inherits = FALSE
  ),
  "\n"
)

cat(
  "\nFunction location:\n"
)

print(
  getAnywhere(goodSamplesGenes)
) 

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("impute", ask = FALSE, update = FALSE) 

cat("\n========== Verify impute ==========\n")

cat(
  "impute installed:",
  "impute" %in% rownames(installed.packages()),
  "\n"
)

cat("impute version:\n")
print(packageVersion("impute")) 

library(impute)

cat(
  "\nimpute loaded:",
  "impute" %in% loadedNamespaces(),
  "\n"
) 

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(
  "preprocessCore",
  ask = FALSE,
  update = FALSE
)

cat("\n========== Verify preprocessCore ==========\n")

cat(
  "preprocessCore installed:",
  "preprocessCore" %in% rownames(installed.packages()),
  "\n"
)

cat("preprocessCore version:\n")
print(packageVersion("preprocessCore"))

library(preprocessCore)

cat(
  "\npreprocessCore loaded:",
  "preprocessCore" %in% loadedNamespaces(),
  "\n"
) 

cat("\n========== Reload WGCNA ==========\n")

library(WGCNA)

cat(
  "WGCNA loaded:",
  "WGCNA" %in% loadedNamespaces(),
  "\n"
)

cat(
  "goodSamplesGenes exists:",
  exists(
    "goodSamplesGenes",
    where = asNamespace("WGCNA"),
    inherits = FALSE
  ),
  "\n"
) 

cat("\n========== Workflow-04-2D: WGCNA pre-QC ==========\n")

# ============================================================
# 1. Build WGCNA input matrix
# ============================================================

datExpr_pre <- t(expr_mad)

cat("\n[1] WGCNA input matrix\n")
cat("Samples:", nrow(datExpr_pre), "\n")
cat("Genes:", ncol(datExpr_pre), "\n")

# ============================================================
# 2. Matrix integrity check
# ============================================================

cat("\n[2] Matrix integrity\n")
cat("NA:", sum(is.na(datExpr_pre)), "\n")
cat("NaN:", sum(is.nan(datExpr_pre)), "\n")
cat("Inf:", sum(is.infinite(datExpr_pre)), "\n")

cat(
  "Duplicated samples:",
  sum(duplicated(rownames(datExpr_pre))),
  "\n"
)

cat(
  "Duplicated genes:",
  sum(duplicated(colnames(datExpr_pre))),
  "\n"
)

# ============================================================
# 3. Confirm clinical alignment
# ============================================================

cat("\n[3] Clinical alignment\n")

cat(
  "Expression samples:",
  nrow(datExpr_pre),
  "\n"
)

cat(
  "Clinical samples:",
  nrow(clinical_clean),
  "\n"
)

cat(
  "Same sample IDs:",
  identical(
    rownames(datExpr_pre),
    clinical_clean$Sample
  ),
  "\n"
)

cat(
  "Expression sample IDs missing in clinical:",
  sum(
    !rownames(datExpr_pre) %in% clinical_clean$Sample
  ),
  "\n"
)

cat(
  "Clinical sample IDs missing in expression:",
  sum(
    !clinical_clean$Sample %in% rownames(datExpr_pre)
  ),
  "\n"
)

# ============================================================
# 4. Confirm tumor-only status
# ============================================================

cat("\n[4] Tumor-only verification\n")

if ("tumor_normal" %in% colnames(clinical_clean)) {
  
  print(
    table(
      clinical_clean$tumor_normal,
      useNA = "ifany"
    )
  )
  
} else {
  
  cat(
    "WARNING: tumor_normal column not found.\n"
  )
}

# ============================================================
# 5. WGCNA goodSamplesGenes QC
# ============================================================

cat("\n[5] Running goodSamplesGenes()\n")

gsg <- goodSamplesGenes(
  datExpr_pre,
  verbose = 3
)

# ============================================================
# 6. Summarize QC result
# ============================================================

cat("\n========== goodSamplesGenes RESULT ==========\n")

cat(
  "allOK:",
  gsg$allOK,
  "\n"
)

cat(
  "Good samples:",
  sum(gsg$goodSamples),
  "\n"
)

cat(
  "Bad samples:",
  sum(!gsg$goodSamples),
  "\n"
)

cat(
  "Good genes:",
  sum(gsg$goodGenes),
  "\n"
)

cat(
  "Bad genes:",
  sum(!gsg$goodGenes),
  "\n"
)

# ============================================================
# 7. List bad samples / genes if present
# ============================================================

if (any(!gsg$goodSamples)) {
  
  cat("\n---------- Bad samples ----------\n")
  print(
    rownames(datExpr_pre)[!gsg$goodSamples]
  )
}

if (any(!gsg$goodGenes)) {
  
  cat("\n---------- Bad genes ----------\n")
  print(
    colnames(datExpr_pre)[!gsg$goodGenes]
  )
}

cat("\n========== Workflow-04-2D Step 1 completed ==========\n")

cat("\n========== Workflow-04-2D Step 2: PCA + Hierarchical Clustering ==========\n")

# ============================================================
# 1. Prepare expression matrix
# ============================================================

datExpr_qc <- datExpr_pre

cat("\n[1] Input matrix\n")
cat("Samples:", nrow(datExpr_qc), "\n")
cat("Genes:", ncol(datExpr_qc), "\n")

# ============================================================
# 2. Select highly variable genes for sample-level QC
#    IMPORTANT:
#    This does NOT change the final WGCNA gene set.
# ============================================================

gene_mad <- apply(
  datExpr_qc,
  2,
  mad,
  na.rm = TRUE
)

top_n <- min(5000, length(gene_mad))

top_genes <- names(
  sort(
    gene_mad,
    decreasing = TRUE
  )
)[seq_len(top_n)]

datExpr_qc_top <- datExpr_qc[, top_genes, drop = FALSE]

cat("\n[2] PCA / clustering QC gene set\n")
cat("Genes used for QC:", ncol(datExpr_qc_top), "\n")
cat("Samples:", nrow(datExpr_qc_top), "\n")

# ============================================================
# 3. PCA
# ============================================================

cat("\n[3] Running PCA...\n")

pca <- prcomp(
  datExpr_qc_top,
  center = TRUE,
  scale. = FALSE
)

pca_var <- pca$sdev^2
pca_var_percent <- 100 * pca_var / sum(pca_var)

pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],
  PC4 = pca$x[, 4],
  PC5 = pca$x[, 5],
  stringsAsFactors = FALSE
)

cat("\nPCA variance explained:\n")
print(
  round(
    pca_var_percent[1:10],
    3
  )
)

# ============================================================
# 4. Clinical annotation for PCA
# ============================================================

clinical_pca <- clinical_clean[
  match(
    pca_df$Sample,
    clinical_clean$Sample
  ),
  ,
  drop = FALSE
]

cat(
  "\nClinical annotation matched:",
  identical(
    pca_df$Sample,
    clinical_pca$Sample
  ),
  "\n"
)

# ============================================================
# 5. Hierarchical clustering
# ============================================================

cat("\n[4] Running hierarchical clustering...\n")

sample_dist <- dist(
  datExpr_qc_top,
  method = "euclidean"
)

sample_tree <- hclust(
  sample_dist,
  method = "average"
)

cat(
  "Hierarchical clustering completed.\n"
)

# ============================================================
# 6. Save QC results
# ============================================================

write.table(
  pca_df,
  file = "result_dir/Workflow-04-2D-2_PCA_scores.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    PC = paste0(
      "PC",
      seq_along(pca_var_percent)
    ),
    Variance_percent = pca_var_percent
  ),
  file = "result_dir/Workflow-04-2D-2_PCA_variance.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    Sample = sample_tree$labels,
    Cluster_order = seq_along(
      sample_tree$order
    )
  ),
  file = "result_dir/Workflow-04-2D-2_hierarchical_clustering_order.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========== Workflow-04-2D Step 2 completed ==========\n")


cat("\n========== Workflow-04-2D Step 3: Sample-level QC ==========\n")

# ============================================================
# 1. Load required packages
# ============================================================

library(WGCNA)
library(ggplot2)

# ============================================================
# 2. PCA data
# ============================================================

pca_plot_df <- pca_df

# Add broad clinical stage
pca_plot_df$Stage <- clinical_pca$pathologic_stage

# Add tumor/normal status
pca_plot_df$TumorNormal <- clinical_pca$tumor_normal

# ============================================================
# 3. PCA plot: PC1 vs PC2
# ============================================================

p_pca_12 <- ggplot(
  pca_plot_df,
  aes(
    x = PC1,
    y = PC2
  )
) +
  geom_point(
    size = 2,
    alpha = 0.7
  ) +
  theme_classic() +
  labs(
    title = "PCA of Primary LUAD Samples",
    x = paste0(
      "PC1 (",
      round(pca_var_percent[1], 2),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(pca_var_percent[2], 2),
      "%)"
    )
  )

print(p_pca_12)

ggsave(
  "figure_dir/Workflow-04-2D-3_PCA_PC1_PC2.pdf",
  p_pca_12,
  width = 7,
  height = 6
)

ggsave(
  "figure_dir/Workflow-04-2D-3_PCA_PC1_PC2.png",
  p_pca_12,
  width = 7,
  height = 6,
  dpi = 300
)

# ============================================================
# 4. PCA plot colored by broad pathological stage
# ============================================================

p_pca_stage <- ggplot(
  pca_plot_df,
  aes(
    x = PC1,
    y = PC2,
    color = Stage
  )
) +
  geom_point(
    size = 2,
    alpha = 0.7
  ) +
  theme_classic() +
  labs(
    title = "PCA of Primary LUAD Samples",
    x = paste0(
      "PC1 (",
      round(pca_var_percent[1], 2),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(pca_var_percent[2], 2),
      "%)"
    ),
    color = "Pathologic stage"
  )

print(p_pca_stage)

ggsave(
  "figure_dir/Workflow-04-2D-3_PCA_PC1_PC2_stage.pdf",
  p_pca_stage,
  width = 8,
  height = 6
)

ggsave(
  "figure_dir/Workflow-04-2D-3_PCA_PC1_PC2_stage.png",
  p_pca_stage,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# 5. Hierarchical clustering
# ============================================================

pdf(
  "Workflow-04-2D-3_Sample_Hierarchical_Clustering.pdf",
  width = 12,
  height = 8
)

plot(
  sample_tree,
  main = "Hierarchical Clustering of Primary LUAD Samples",
  xlab = "",
  sub = "",
  cex = 0.45
)

dev.off()

# ============================================================
# 6. Sample connectivity
# ============================================================

cat("\n[6] Calculating sample connectivity...\n")

# Pearson correlation among samples
sample_cor <- cor(
  t(datExpr_qc_top),
  use = "pairwise.complete.obs",
  method = "pearson"
)

# Connectivity = sum of correlations to all other samples
sample_connectivity <- rowSums(
  sample_cor
) - 1

connectivity_df <- data.frame(
  Sample = rownames(sample_cor),
  Connectivity = sample_connectivity,
  stringsAsFactors = FALSE
)

# Rank samples
connectivity_df$Connectivity_rank <- rank(
  connectivity_df$Connectivity,
  ties.method = "average"
)

# Percentile
connectivity_df$Connectivity_percentile <- 100 *
  connectivity_df$Connectivity_rank /
  nrow(connectivity_df)

# Add clinical information
connectivity_df$Stage <- clinical_clean$pathologic_stage[
  match(
    connectivity_df$Sample,
    clinical_clean$Sample
  )
]

connectivity_df$OS_status <- clinical_clean$OS.status[
  match(
    connectivity_df$Sample,
    clinical_clean$Sample
  )
]

# ============================================================
# 7. Connectivity summary
# ============================================================

cat("\n========== Sample Connectivity Summary ==========\n")

print(
  summary(
    connectivity_df$Connectivity
  )
)

cat(
  "\nLowest 10 connectivity samples:\n"
)

print(
  head(
    connectivity_df[
      order(
        connectivity_df$Connectivity
      ),
    ],
    10
  )
)

cat(
  "\nHighest 10 connectivity samples:\n"
)

print(
  head(
    connectivity_df[
      order(
        -connectivity_df$Connectivity
      ),
    ],
    10
  )
)

# ============================================================
# 8. Save connectivity results
# ============================================================

write.table(
  connectivity_df,
  file = "input_dir/Workflow-04-2D-3_sample_connectivity.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nConnectivity results saved.\n"
)

cat(
  "\n========== Workflow-04-2D Step 3 completed ==========\n"
) 


cat("\n========== Workflow-04-2D Step 4: TCGA-44 Connectivity Audit ==========\n")

# ============================================================
# 1. Extract TCGA Tissue Source Site (TSS)
# ============================================================

connectivity_df$TSS <- substr(
  connectivity_df$Sample,
  6,
  7
)

cat(
  "\nNumber of TSS groups:",
  length(unique(connectivity_df$TSS)),
  "\n"
)

# ============================================================
# 2. Generate TSS-level connectivity summary
# ============================================================

tss_split <- split(
  connectivity_df$Connectivity,
  connectivity_df$TSS
)

tss_summary_df <- data.frame(
  TSS = names(tss_split),
  N = sapply(
    tss_split,
    length
  ),
  Median = sapply(
    tss_split,
    median
  ),
  Mean = sapply(
    tss_split,
    mean
  ),
  Min = sapply(
    tss_split,
    min
  ),
  Q1 = sapply(
    tss_split,
    function(x) quantile(x, 0.25)
  ),
  Q3 = sapply(
    tss_split,
    function(x) quantile(x, 0.75)
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# 3. Sort by median connectivity
# ============================================================

tss_summary_df <- tss_summary_df[
  order(
    tss_summary_df$Median
  ),
  ,
  drop = FALSE
]

rownames(tss_summary_df) <- NULL

cat("\n========== Connectivity by TSS ==========\n")

print(
  tss_summary_df,
  row.names = FALSE
)

# ============================================================
# 4. Specifically inspect TCGA-44
# ============================================================

tcga44 <- connectivity_df[
  connectivity_df$TSS == "44",
  ,
  drop = FALSE
]

cat("\n========== TCGA-44 Connectivity ==========\n")

cat(
  "Number of TCGA-44 samples:",
  nrow(tcga44),
  "\n"
)

if (nrow(tcga44) > 0) {
  
  cat(
    "Median connectivity:",
    round(
      median(tcga44$Connectivity),
      3
    ),
    "\n"
  )
  
  cat(
    "Mean connectivity:",
    round(
      mean(tcga44$Connectivity),
      3
    ),
    "\n"
  )
  
  cat(
    "Minimum connectivity:",
    round(
      min(tcga44$Connectivity),
      3
    ),
    "\n"
  )
  
  cat(
    "Maximum connectivity:",
    round(
      max(tcga44$Connectivity),
      3
    ),
    "\n"
  )
  
  cat(
    "\nTCGA-44 samples sorted by connectivity:\n"
  )
  
  print(
    tcga44[
      order(
        tcga44$Connectivity
      ),
      c(
        "Sample",
        "Connectivity",
        "Connectivity_percentile",
        "Stage",
        "OS_status"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No TCGA-44 samples found.\n"
  )
}

# ============================================================
# 5. Identify the lowest 10% connectivity samples
#    This is descriptive only.
#    DO NOT remove samples based on this threshold.
# ============================================================

connectivity_cutoff_10 <- quantile(
  connectivity_df$Connectivity,
  0.10
)

low10 <- connectivity_df[
  connectivity_df$Connectivity <= connectivity_cutoff_10,
  ,
  drop = FALSE
]

cat(
  "\n========== Lowest 10% Connectivity ==========\n"
)

cat(
  "10th percentile cutoff:",
  round(
    connectivity_cutoff_10,
    3
  ),
  "\n"
)

cat(
  "Number of samples:",
  nrow(low10),
  "\n"
)

cat(
  "Number of TCGA-44 samples:",
  sum(
    low10$TSS == "44"
  ),
  "\n"
)

cat(
  "Percentage of lowest 10% that are TCGA-44:",
  round(
    100 *
      mean(
        low10$TSS == "44"
      ),
    2
  ),
  "%\n"
)

# ============================================================
# 6. Save results
# ============================================================

write.table(
  tss_summary_df,
  file = "result_dir/Workflow-04-2D-4_TSS_connectivity_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  tcga44[
    order(
      tcga44$Connectivity
    ),
    ,
    drop = FALSE
  ],
  file = "result_dir/Workflow-04-2D-4_TCGA44_connectivity.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  low10[
    order(
      low10$Connectivity
    ),
    ,
    drop = FALSE
  ],
  file = "result_dir/Workflow-04-2D-4_low10_connectivity_samples.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\n========== Workflow-04-2D Step 4 completed ==========\n"
)

save.image(
  file = "result_dir/Workflow-04-2D_current_workspace.RData"
)  

savehistory(
  file = "result_dir/Workflow-04-2D_R_history.Rhistory"
)

load("O:/Rproject/Paper1_tumor-only WGCNA/Workflow-04-2D_current_workspace.RData")
