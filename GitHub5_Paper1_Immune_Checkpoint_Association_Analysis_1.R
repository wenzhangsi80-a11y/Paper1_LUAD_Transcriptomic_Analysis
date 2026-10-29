###GitHub5_Paper1_Immune_Checkpoint_Association_Analysis_1

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

# Resume Workflow-08-Immune Checkpoint Association Analysis
# ============================================================

work_dir <- "input_dir/Paper 1_Step 8_Immune Checkpoint"

setwd(work_dir)

getwd()

rm(list = ls())
gc()
getwd()

# ============================================================
# Paper 1
# WORKFLOW-08
# STEP 8A
# Data source inventory
# ============================================================

cat("\n============================================================\n")
cat("WORKFLOW-08 / STEP 8A\n")
cat("Data source inventory\n")
cat("============================================================\n\n")

wgcna_dir <- "input_dir/Paper1_tumor-only WGCNA"

step7_dir <- "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis/STEP7B-7E-10"

workflow08_dir <- "input_dir/Paper1_WORKFLOW08_Immune_Checkpoint/STEP8A_Data_Lock"

cat("WGCNA directory:\n")
cat(wgcna_dir, "\n\n")

cat("STEP7 source directory:\n")
cat(step7_dir, "\n\n")

cat("WORKFLOW08 directory:\n")
cat(workflow08_dir, "\n\n")

cat("------------------------------------------------------------\n")
cat("WGCNA files\n")
cat("------------------------------------------------------------\n")

wgcna_files <- list.files(
  wgcna_dir,
  recursive = TRUE,
  full.names = TRUE
)

print(wgcna_files)

cat("\n------------------------------------------------------------\n")
cat("STEP7 files\n")
cat("------------------------------------------------------------\n")

step7_files <- list.files(
  step7_dir,
  recursive = TRUE,
  full.names = TRUE
)

print(step7_files)

# ============================================================
# Paper 1
# WORKFLOW-08
# STEP 8A-1
# Confirm tumor-only WGCNA source files
# ============================================================

rm(list = ls())
gc()

wgcna_dir <- "input_dir/Paper1_tumor-only WGCNA"

expr_file <- file.path(
  wgcna_dir,
  "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt"
)

clinical_file <- file.path(
  wgcna_dir,
  "Workflow-04-2A_clinical_aligned.txt"
)

cat("\n============================================================\n")
cat("WORKFLOW-08 / STEP 8A-1\n")
cat("Tumor-only WGCNA source file confirmation\n")
cat("============================================================\n\n")

cat("Expression file exists: ",
    file.exists(expr_file), "\n")

cat("Clinical file exists: ",
    file.exists(clinical_file), "\n\n")

cat("Expression file:\n")
cat(expr_file, "\n\n")

cat("Clinical file:\n")
cat(clinical_file, "\n\n")


## ------------------------------------------------------------
# Expression file: inspect header
# ------------------------------------------------------------

expr_header <- read.delim(
  expr_file,
  header = TRUE,
  sep = "\t",
  nrows = 3,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("Expression file preview\n")
cat("============================================================\n")

print(dim(expr_header))
print(colnames(expr_header)[1:min(10, ncol(expr_header))])
print(expr_header[1:min(3, nrow(expr_header)),
                  1:min(6, ncol(expr_header))])


# ------------------------------------------------------------
# Clinical file: inspect header
# ------------------------------------------------------------

clinical_preview <- read.delim(
  clinical_file,
  header = TRUE,
  sep = "\t",
  nrows = 5,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("Clinical file preview\n")
cat("============================================================\n")

print(dim(clinical_preview))
print(colnames(clinical_preview))
print(clinical_preview)

# ============================================================
# Full dimensions
# ============================================================

expr_all <- read.delim(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

clinical_all <- read.delim(
  clinical_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("FULL DIMENSIONS\n")
cat("============================================================\n")

cat("Expression:\n")
print(dim(expr_all))

cat("\nClinical:\n")
print(dim(clinical_all))


# ============================================================
# Paper 1
# WORKFLOW-08
# STEP 8A-2
# Strict sample ID consistency check
# ============================================================

cat("\n============================================================\n")
cat("WORKFLOW-08 / STEP 8A-2\n")
cat("Strict expression–clinical sample ID consistency check\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define sample IDs
# ------------------------------------------------------------

expr_samples <- colnames(expr_all)[-1]
clinical_samples <- clinical_all$Sample


# ------------------------------------------------------------
# 2. Basic dimensions
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("1. BASIC SAMPLE DIMENSIONS\n")
cat("------------------------------------------------------------\n")

cat("Expression sample count: ", length(expr_samples), "\n")
cat("Clinical sample count:   ", length(clinical_samples), "\n\n")


# ------------------------------------------------------------
# 3. Check missing sample IDs
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("2. MISSING SAMPLE ID CHECK\n")
cat("------------------------------------------------------------\n")

expr_missing <- sum(
  is.na(expr_samples) |
    trimws(expr_samples) == ""
)

clinical_missing <- sum(
  is.na(clinical_samples) |
    trimws(clinical_samples) == ""
)

cat("Expression missing/blank IDs: ", expr_missing, "\n")
cat("Clinical missing/blank IDs:   ", clinical_missing, "\n\n")


# ------------------------------------------------------------
# 4. Check duplicated sample IDs
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("3. DUPLICATE SAMPLE ID CHECK\n")
cat("------------------------------------------------------------\n")

expr_dup <- unique(expr_samples[duplicated(expr_samples)])
clinical_dup <- unique(clinical_samples[duplicated(clinical_samples)])

cat("Expression duplicated IDs: ", length(expr_dup), "\n")

if (length(expr_dup) > 0) {
  cat("Duplicated expression IDs:\n")
  print(expr_dup)
}

cat("\nClinical duplicated IDs: ", length(clinical_dup), "\n")

if (length(clinical_dup) > 0) {
  cat("Duplicated clinical IDs:\n")
  print(clinical_dup)
}

cat("\n")


# ------------------------------------------------------------
# 5. Expression-only samples
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("4. EXPRESSION-ONLY SAMPLE CHECK\n")
cat("------------------------------------------------------------\n")

expr_only <- setdiff(
  expr_samples,
  clinical_samples
)

cat("Expression-only samples: ", length(expr_only), "\n")

if (length(expr_only) > 0) {
  print(expr_only)
}

cat("\n")


# ------------------------------------------------------------
# 6. Clinical-only samples
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("5. CLINICAL-ONLY SAMPLE CHECK\n")
cat("------------------------------------------------------------\n")

clinical_only <- setdiff(
  clinical_samples,
  expr_samples
)

cat("Clinical-only samples: ", length(clinical_only), "\n")

if (length(clinical_only) > 0) {
  print(clinical_only)
}

cat("\n")


# ------------------------------------------------------------
# 7. Exact set equality
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("6. SAMPLE ID SET EQUALITY\n")
cat("------------------------------------------------------------\n")

same_sample_set <- setequal(
  expr_samples,
  clinical_samples
)

cat("Expression and clinical sample ID sets identical: ",
    same_sample_set,
    "\n\n")


# ------------------------------------------------------------
# 8. Check sample order
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("7. SAMPLE ORDER CHECK\n")
cat("------------------------------------------------------------\n")

same_sample_order <- identical(
  expr_samples,
  clinical_samples
)

cat("Expression and clinical sample order identical: ",
    same_sample_order,
    "\n\n")


# ------------------------------------------------------------
# 9. Show first and last sample IDs for visual confirmation
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("8. SAMPLE ID ORDER PREVIEW\n")
cat("------------------------------------------------------------\n")

cat("\nFirst 10 expression sample IDs:\n")
print(expr_samples[1:min(10, length(expr_samples))])

cat("\nFirst 10 clinical sample IDs:\n")
print(clinical_samples[1:min(10, length(clinical_samples))])

cat("\nLast 10 expression sample IDs:\n")
print(
  tail(
    expr_samples,
    min(10, length(expr_samples))
  )
)

cat("\nLast 10 clinical sample IDs:\n")
print(
  tail(
    clinical_samples,
    min(10, length(clinical_samples))
  )
)

cat("\n")


# ------------------------------------------------------------
# 10. Check tumor_normal variable
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("9. TUMOR/NORMAL STATUS CHECK\n")
cat("------------------------------------------------------------\n")

cat("Unique tumor_normal values:\n")
print(unique(clinical_all$tumor_normal))

cat("\nFrequency table:\n")
print(table(
  clinical_all$tumor_normal,
  useNA = "ifany"
))

all_tumor <- all(
  !is.na(clinical_all$tumor_normal) &
    clinical_all$tumor_normal == "tumor"
)

cat("\nAll clinical samples labeled 'tumor': ",
    all_tumor,
    "\n\n")


# ------------------------------------------------------------
# 11. Check duplicate gene IDs
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("10. GENE ID CHECK\n")
cat("------------------------------------------------------------\n")

gene_ids <- expr_all[[1]]

gene_missing <- sum(
  is.na(gene_ids) |
    trimws(gene_ids) == ""
)

gene_dup <- unique(
  gene_ids[duplicated(gene_ids)]
)

cat("Total gene IDs:             ", length(gene_ids), "\n")
cat("Missing/blank gene IDs:     ", gene_missing, "\n")
cat("Duplicated gene IDs:        ", length(gene_dup), "\n")

if (length(gene_dup) > 0) {
  cat("\nDuplicated gene IDs:\n")
  print(gene_dup)
}

cat("\n")


# ------------------------------------------------------------
# 12. Final QC decision
# ------------------------------------------------------------

cat("============================================================\n")
cat("FINAL STEP 8A-2 QC SUMMARY\n")
cat("============================================================\n")

cat(
  "Expression samples = Clinical samples: ",
  length(expr_samples) == length(clinical_samples),
  "\n"
)

cat(
  "No missing expression sample IDs:      ",
  expr_missing == 0,
  "\n"
)

cat(
  "No missing clinical sample IDs:        ",
  clinical_missing == 0,
  "\n"
)

cat(
  "No duplicated expression sample IDs:   ",
  length(expr_dup) == 0,
  "\n"
)

cat(
  "No duplicated clinical sample IDs:     ",
  length(clinical_dup) == 0,
  "\n"
)

cat(
  "No expression-only samples:            ",
  length(expr_only) == 0,
  "\n"
)

cat(
  "No clinical-only samples:              ",
  length(clinical_only) == 0,
  "\n"
)

cat(
  "Sample ID sets identical:              ",
  same_sample_set,
  "\n"
)

cat(
  "Sample order identical:                ",
  same_sample_order,
  "\n"
)

cat(
  "All clinical samples are tumor:        ",
  all_tumor,
  "\n"
)

cat(
  "No missing gene IDs:                   ",
  gene_missing == 0,
  "\n"
)

cat(
  "No duplicated gene IDs:                ",
  length(gene_dup) == 0,
  "\n"
)


# ------------------------------------------------------------
# 13. Overall pass/fail
# ------------------------------------------------------------

step8A2_pass <- (
  length(expr_samples) == length(clinical_samples) &&
    expr_missing == 0 &&
    clinical_missing == 0 &&
    length(expr_dup) == 0 &&
    length(clinical_dup) == 0 &&
    length(expr_only) == 0 &&
    length(clinical_only) == 0 &&
    same_sample_set &&
    same_sample_order &&
    all_tumor &&
    gene_missing == 0 &&
    length(gene_dup) == 0
)

cat("\n============================================================\n")
cat(
  "STEP 8A-2 OVERALL RESULT: ",
  ifelse(step8A2_pass, "PASS", "FAIL"),
  "\n"
)
cat("============================================================\n")


# ============================================================
# Paper 1
# WORKFLOW-08
# STEP 8A-3
# Numerical integrity and clinical variable QC
# ============================================================

cat("\n============================================================\n")
cat("WORKFLOW-08 / STEP 8A-3\n")
cat("Numerical integrity and clinical variable quality control\n")
cat("============================================================\n\n")


# ============================================================
# PART A
# EXPRESSION MATRIX NUMERICAL QC
# ============================================================


# ------------------------------------------------------------
# A1. Basic structure
# ------------------------------------------------------------

cat("============================================================\n")
cat("PART A: EXPRESSION MATRIX NUMERICAL QC\n")
cat("============================================================\n\n")

cat("Expression matrix dimensions:\n")
print(dim(expr_all))

cat("\nExpression column classes:\n")
print(table(sapply(expr_all, class)))

cat("\nGene ID column class:\n")
print(class(expr_all[[1]]))

cat("\n")


# ------------------------------------------------------------
# A2. Check expression columns are numeric
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("A2. EXPRESSION COLUMN TYPE CHECK\n")
cat("------------------------------------------------------------\n")

expr_numeric <- sapply(
  expr_all[, -1, drop = FALSE],
  is.numeric
)

cat("Number of expression columns: ",
    length(expr_numeric), "\n")

cat("Numeric expression columns: ",
    sum(expr_numeric), "\n")

cat("Non-numeric expression columns: ",
    sum(!expr_numeric), "\n")

if (any(!expr_numeric)) {
  cat("\nNon-numeric columns:\n")
  print(names(expr_numeric)[!expr_numeric])
}

cat("\n")


# ------------------------------------------------------------
# A3. Missing values
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("A3. MISSING VALUE CHECK\n")
cat("------------------------------------------------------------\n")

expr_numeric_matrix <- as.matrix(
  expr_all[, -1, drop = FALSE]
)

expr_na_total <- sum(
  is.na(expr_numeric_matrix)
)

cat("Total NA values: ",
    expr_na_total, "\n")

cat("Rows containing NA: ",
    sum(rowSums(is.na(expr_numeric_matrix)) > 0),
    "\n")

cat("Columns containing NA: ",
    sum(colSums(is.na(expr_numeric_matrix)) > 0),
    "\n")

cat("\n")


# ------------------------------------------------------------
# A4. NaN check
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("A4. NaN CHECK\n")
cat("------------------------------------------------------------\n")

expr_nan_total <- sum(
  is.nan(expr_numeric_matrix)
)

cat("Total NaN values: ",
    expr_nan_total, "\n\n")


# ------------------------------------------------------------
# A5. Infinite values
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("A5. INFINITE VALUE CHECK\n")
cat("------------------------------------------------------------\n")

expr_inf_total <- sum(
  is.infinite(expr_numeric_matrix)
)

cat("Total Inf/-Inf values: ",
    expr_inf_total, "\n\n")


# ------------------------------------------------------------
# A6. Finite-value check
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("A6. FINITE VALUE CHECK\n")
cat("------------------------------------------------------------\n")

expr_finite <- all(
  is.finite(expr_numeric_matrix)
)

cat("All expression values are finite: ",
    expr_finite,
    "\n\n")


# ------------------------------------------------------------
# A7. Expression value range
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("A7. EXPRESSION VALUE RANGE\n")
cat("------------------------------------------------------------\n")

expr_min <- min(
  expr_numeric_matrix,
  na.rm = TRUE
)

expr_max <- max(
  expr_numeric_matrix,
  na.rm = TRUE
)

expr_mean <- mean(
  expr_numeric_matrix,
  na.rm = TRUE
)

expr_median <- median(
  expr_numeric_matrix,
  na.rm = TRUE
)

cat("Minimum expression value: ",
    expr_min, "\n")

cat("Maximum expression value: ",
    expr_max, "\n")

cat("Mean expression value:    ",
    expr_mean, "\n")

cat("Median expression value:  ",
    expr_median, "\n\n")


# ------------------------------------------------------------
# A8. Negative expression values
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("A8. NEGATIVE VALUE CHECK\n")
cat("------------------------------------------------------------\n")

expr_negative_total <- sum(
  expr_numeric_matrix < 0,
  na.rm = TRUE
)

cat("Total negative expression values: ",
    expr_negative_total,
    "\n")

cat("Negative values present: ",
    expr_negative_total > 0,
    "\n\n")


# ------------------------------------------------------------
# A9. Zero-expression burden per sample
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("A9. ZERO-EXPRESSION BURDEN PER SAMPLE\n")
cat("------------------------------------------------------------\n")

zero_count_per_sample <- colSums(
  expr_numeric_matrix == 0,
  na.rm = TRUE
)

zero_fraction_per_sample <- zero_count_per_sample /
  nrow(expr_numeric_matrix)

cat("Minimum zero fraction: ",
    min(zero_fraction_per_sample),
    "\n")

cat("Median zero fraction:  ",
    median(zero_fraction_per_sample),
    "\n")

cat("Maximum zero fraction: ",
    max(zero_fraction_per_sample),
    "\n")

cat("\nSamples with >95% zero expression:\n")
high_zero_samples <- names(
  zero_fraction_per_sample[
    zero_fraction_per_sample > 0.95
  ]
)

cat("Number of samples: ",
    length(high_zero_samples),
    "\n")

if (length(high_zero_samples) > 0) {
  print(high_zero_samples)
}

cat("\n")


# ------------------------------------------------------------
# A10. Gene-level zero burden
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("A10. ZERO-EXPRESSION BURDEN PER GENE\n")
cat("------------------------------------------------------------\n")

zero_count_per_gene <- rowSums(
  expr_numeric_matrix == 0,
  na.rm = TRUE
)

zero_fraction_per_gene <- zero_count_per_gene /
  ncol(expr_numeric_matrix)

cat("Minimum gene zero fraction: ",
    min(zero_fraction_per_gene),
    "\n")

cat("Median gene zero fraction:  ",
    median(zero_fraction_per_gene),
    "\n")

cat("Maximum gene zero fraction: ",
    max(zero_fraction_per_gene),
    "\n")

cat("\nGenes with >95% zero expression:\n")

high_zero_genes <- sum(
  zero_fraction_per_gene > 0.95
)

cat("Number of genes: ",
    high_zero_genes,
    "\n\n")


# ============================================================
# PART B
# CLINICAL VARIABLE QC
# ============================================================

cat("============================================================\n")
cat("PART B: CLINICAL VARIABLE QC\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# B1. Clinical dimensions
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("B1. CLINICAL DATA DIMENSIONS\n")
cat("------------------------------------------------------------\n")

print(dim(clinical_all))

cat("\n")


# ------------------------------------------------------------
# B2. Missing values by clinical variable
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("B2. MISSING VALUES BY CLINICAL VARIABLE\n")
cat("------------------------------------------------------------\n")

clinical_missing_summary <- data.frame(
  Variable = names(clinical_all),
  Missing_N = sapply(
    clinical_all,
    function(x) sum(is.na(x) | trimws(as.character(x)) == "")
  ),
  Missing_Percent = sapply(
    clinical_all,
    function(x)
      100 * sum(
        is.na(x) |
          trimws(as.character(x)) == ""
      ) / length(x)
  ),
  row.names = NULL
)

print(clinical_missing_summary)

cat("\n")


# ------------------------------------------------------------
# B3. OS.time QC
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("B3. OS.time QC\n")
cat("------------------------------------------------------------\n")

cat("Class:\n")
print(class(clinical_all$OS.time))

cat("\nSummary:\n")
print(summary(clinical_all$OS.time))

cat("\nMissing OS.time: ",
    sum(is.na(clinical_all$OS.time)),
    "\n")

cat("Non-finite OS.time: ",
    sum(!is.finite(clinical_all$OS.time)),
    "\n")

cat("Negative OS.time: ",
    sum(
      clinical_all$OS.time < 0,
      na.rm = TRUE
    ),
    "\n")

cat("Zero OS.time: ",
    sum(
      clinical_all$OS.time == 0,
      na.rm = TRUE
    ),
    "\n\n")


# ------------------------------------------------------------
# B4. OS.status QC
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("B4. OS.status QC\n")
cat("------------------------------------------------------------\n")

cat("Class:\n")
print(class(clinical_all$OS.status))

cat("\nObserved values:\n")
print(
  table(
    clinical_all$OS.status,
    useNA = "ifany"
  )
)

cat("\n")


# Convert status to normalized character representation
os_status_chr <- trimws(
  tolower(
    as.character(clinical_all$OS.status)
  )
)

cat("Normalized OS.status values:\n")
print(
  table(
    os_status_chr,
    useNA = "ifany"
  )
)

cat("\n")


# ------------------------------------------------------------
# B5. Verify expected survival status coding
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("B5. OS.status CODING CHECK\n")
cat("------------------------------------------------------------\n")

expected_status <- c(
  "alive",
  "dead"
)

unexpected_status <- setdiff(
  unique(
    os_status_chr[
      !is.na(os_status_chr)
    ]
  ),
  expected_status
)

cat("Unexpected OS.status values: ",
    length(unexpected_status),
    "\n")

if (length(unexpected_status) > 0) {
  print(unexpected_status)
}

cat("\n")


# ------------------------------------------------------------
# B6. Gender QC
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("B6. GENDER QC\n")
cat("------------------------------------------------------------\n")

print(
  table(
    clinical_all$Gender,
    useNA = "ifany"
  )
)

cat("\n")


# ------------------------------------------------------------
# B7. Age QC
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("B7. AGE QC\n")
cat("------------------------------------------------------------\n")

cat("Class:\n")
print(class(clinical_all$Age))

cat("\nSummary:\n")
print(summary(clinical_all$Age))

cat("\nMissing Age: ",
    sum(is.na(clinical_all$Age)),
    "\n")

cat("Non-finite Age: ",
    sum(!is.finite(clinical_all$Age)),
    "\n")

cat("Negative Age: ",
    sum(
      clinical_all$Age < 0,
      na.rm = TRUE
    ),
    "\n")

cat("\n")


# ------------------------------------------------------------
# B8. Pathologic stage variables
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("B8. PATHOLOGIC STAGE VARIABLES QC\n")
cat("------------------------------------------------------------\n")

stage_variables <- c(
  "pathologic_stage_2",
  "pathologic_t",
  "pathologic_n",
  "pathologic_m",
  "pathologic_stage"
)

for (v in stage_variables) {
  
  cat("\nVariable: ", v, "\n")
  print(
    table(
      clinical_all[[v]],
      useNA = "ifany"
    )
  )
}


# ------------------------------------------------------------
# B9. Patient-level clinical ID consistency
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("B9. PATIENT ID QC\n")
cat("------------------------------------------------------------\n")

patient_ids <- clinical_all$patient_id

cat("Missing patient IDs: ",
    sum(
      is.na(patient_ids) |
        trimws(patient_ids) == ""
    ),
    "\n")

cat("Duplicated patient IDs: ",
    sum(duplicated(patient_ids)),
    "\n")

cat("\nNumber of unique patient IDs: ",
    length(unique(patient_ids)),
    "\n")


# ============================================================
# PART C
# CROSS-CHECK SURVIVAL VARIABLES
# ============================================================

cat("\n============================================================\n")
cat("PART C: SURVIVAL VARIABLE CROSS-CHECK\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# C1. OS.time / OS.status completeness
# ------------------------------------------------------------

survival_complete <- complete.cases(
  clinical_all[, c("OS.time", "OS.status")]
)

cat("Complete OS.time + OS.status cases: ",
    sum(survival_complete),
    "\n")

cat("Incomplete OS.time + OS.status cases: ",
    sum(!survival_complete),
    "\n\n")


# ------------------------------------------------------------
# C2. Cross-tabulation OS.status vs OS.time = 0
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("C2. OS.status vs OS.time = 0\n")
cat("------------------------------------------------------------\n")

zero_time <- clinical_all$OS.time == 0

print(
  table(
    OS_status = os_status_chr,
    OS_time_zero = zero_time,
    useNA = "ifany"
  )
)

cat("\n")


# ============================================================
# PART D
# FINAL STEP 8A-3 QC SUMMARY
# ============================================================

cat("============================================================\n")
cat("FINAL STEP 8A-3 QC SUMMARY\n")
cat("============================================================\n\n")


# Expression checks
expr_type_pass <- all(expr_numeric)
expr_na_pass <- expr_na_total == 0
expr_nan_pass <- expr_nan_total == 0
expr_inf_pass <- expr_inf_total == 0
expr_finite_pass <- expr_finite
expr_negative_pass <- expr_negative_total == 0


# Clinical checks
clinical_os_time_pass <- (
  sum(is.na(clinical_all$OS.time)) == 0 &&
    sum(!is.finite(clinical_all$OS.time)) == 0 &&
    sum(clinical_all$OS.time < 0, na.rm = TRUE) == 0
)

clinical_os_status_pass <- (
  sum(is.na(os_status_chr)) == 0 &&
    length(unexpected_status) == 0
)

clinical_tumor_pass <- all(
  !is.na(clinical_all$tumor_normal) &
    clinical_all$tumor_normal == "tumor"
)


cat("Expression columns numeric:             ",
    expr_type_pass, "\n")

cat("No expression NA:                       ",
    expr_na_pass, "\n")

cat("No expression NaN:                      ",
    expr_nan_pass, "\n")

cat("No expression Inf/-Inf:                 ",
    expr_inf_pass, "\n")

cat("All expression values finite:           ",
    expr_finite_pass, "\n")

cat("No negative expression values:          ",
    expr_negative_pass, "\n")

cat("OS.time numerically valid:              ",
    clinical_os_time_pass, "\n")

cat("OS.status coding valid:                 ",
    clinical_os_status_pass, "\n")

cat("All clinical samples are tumor:         ",
    clinical_tumor_pass, "\n")


# ------------------------------------------------------------
# Overall PASS/FAIL
# ------------------------------------------------------------

step8A3_pass <- (
  expr_type_pass &&
    expr_na_pass &&
    expr_nan_pass &&
    expr_inf_pass &&
    expr_finite_pass &&
    expr_negative_pass &&
    clinical_os_time_pass &&
    clinical_os_status_pass &&
    clinical_tumor_pass
)

cat("\n============================================================\n")
cat(
  "STEP 8A-3 OVERALL RESULT: ",
  ifelse(step8A3_pass, "PASS", "FAIL"),
  "\n"
)
cat("============================================================\n")

# ============================================================
# Paper 1
# WORKFLOW-08
# STEP 8A-3 Diagnostic
# Investigate clinical QC failures
# ============================================================

cat("\n============================================================\n")
cat("WORKFLOW-08 / STEP 8A-3 DIAGNOSTIC\n")
cat("Investigate OS.time and OS.status QC failures\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. OS.time detailed diagnostic
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("1. OS.time DETAILED DIAGNOSTIC\n")
cat("------------------------------------------------------------\n")

cat("Class:\n")
print(class(clinical_all$OS.time))

cat("\nLength:\n")
print(length(clinical_all$OS.time))

cat("\nNA count:\n")
print(sum(is.na(clinical_all$OS.time)))

cat("\nNaN count:\n")
print(sum(is.nan(clinical_all$OS.time)))

cat("\nInfinite count:\n")
print(sum(!is.finite(clinical_all$OS.time)))

cat("\nNegative count:\n")
print(sum(
  clinical_all$OS.time < 0,
  na.rm = TRUE
))

cat("\nZero count:\n")
print(sum(
  clinical_all$OS.time == 0,
  na.rm = TRUE
))

cat("\nPositive count:\n")
print(sum(
  clinical_all$OS.time > 0,
  na.rm = TRUE
))

cat("\nSummary:\n")
print(summary(clinical_all$OS.time))


# ------------------------------------------------------------
# 2. Identify problematic OS.time rows
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("2. PROBLEMATIC OS.time ROWS\n")
cat("------------------------------------------------------------\n")

os_time_problem <- (
  is.na(clinical_all$OS.time) |
    is.nan(clinical_all$OS.time) |
    !is.finite(clinical_all$OS.time) |
    clinical_all$OS.time < 0
)

cat("Number of problematic OS.time rows: ",
    sum(os_time_problem),
    "\n\n")

if (sum(os_time_problem) > 0) {
  
  print(
    clinical_all[
      os_time_problem,
      c(
        "Sample",
        "patient_id",
        "OS.time",
        "OS.status",
        "pathologic_stage"
      ),
      drop = FALSE
    ]
  )
  
} else {
  
  cat("No missing, non-finite, or negative OS.time values.\n")
  cat("Therefore any previous FALSE must be investigated separately.\n")
  
}


# ------------------------------------------------------------
# 3. OS.status raw values
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("3. OS.status RAW VALUE DIAGNOSTIC\n")
cat("------------------------------------------------------------\n")

cat("Raw class:\n")
print(class(clinical_all$OS.status))

cat("\nRaw unique values:\n")
print(unique(clinical_all$OS.status))

cat("\nRaw frequency table:\n")
print(
  table(
    clinical_all$OS.status,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 4. OS.status normalized values
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("4. OS.status NORMALIZED VALUE DIAGNOSTIC\n")
cat("------------------------------------------------------------\n")

os_status_chr <- trimws(
  tolower(
    as.character(clinical_all$OS.status)
  )
)

cat("Normalized unique values:\n")
print(unique(os_status_chr))

cat("\nNormalized frequency table:\n")
print(
  table(
    os_status_chr,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 5. Unexpected OS.status values
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("5. UNEXPECTED OS.status VALUES\n")
cat("------------------------------------------------------------\n")

expected_status <- c(
  "alive",
  "dead"
)

unexpected_status <- setdiff(
  unique(
    os_status_chr[
      !is.na(os_status_chr)
    ]
  ),
  expected_status
)

cat("Expected values:\n")
print(expected_status)

cat("\nUnexpected values:\n")

if (length(unexpected_status) == 0) {
  cat("NONE\n")
} else {
  print(unexpected_status)
}


# ------------------------------------------------------------
# 6. Identify problematic OS.status rows
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("6. PROBLEMATIC OS.status ROWS\n")
cat("------------------------------------------------------------\n")

os_status_problem <- (
  is.na(os_status_chr) |
    !(os_status_chr %in% expected_status)
)

cat("Number of problematic OS.status rows: ",
    sum(os_status_problem),
    "\n\n")

if (sum(os_status_problem) > 0) {
  
  print(
    clinical_all[
      os_status_problem,
      c(
        "Sample",
        "patient_id",
        "OS.time",
        "OS.status",
        "pathologic_stage"
      ),
      drop = FALSE
    ]
  )
  
} else {
  
  cat("No problematic OS.status rows.\n")
  
}

# ============================================================
# STEP 8A-3
# Corrected OS.time / OS.status cross-check
# ============================================================

cat("\n============================================================\n")
cat("STEP 8A-3 / OS.time - OS.status CROSS-CHECK\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. OS.time category
# ------------------------------------------------------------

os_time_category <- ifelse(
  is.na(clinical_all$OS.time),
  "NA",
  ifelse(
    clinical_all$OS.time == 0,
    "0",
    ifelse(
      clinical_all$OS.time < 0,
      "<0",
      ">0"
    )
  )
)


# ------------------------------------------------------------
# 2. Normalized OS.status
# ------------------------------------------------------------

os_status_chr <- trimws(
  tolower(
    as.character(clinical_all$OS.status)
  )
)


# ------------------------------------------------------------
# 3. Cross-tabulation
# ------------------------------------------------------------

cat("OS.time category × OS.status:\n\n")

print(
  table(
    os_time_category,
    os_status_chr,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 4. OS.time detailed counts
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("OS.time detailed counts\n")
cat("------------------------------------------------------------\n")

cat("NA:       ",
    sum(is.na(clinical_all$OS.time)), "\n")

cat("NaN:      ",
    sum(is.nan(clinical_all$OS.time)), "\n")

cat("Inf/-Inf: ",
    sum(!is.finite(clinical_all$OS.time)), "\n")

cat("< 0:     ",
    sum(clinical_all$OS.time < 0, na.rm = TRUE), "\n")

cat("= 0:     ",
    sum(clinical_all$OS.time == 0, na.rm = TRUE), "\n")

cat("> 0:     ",
    sum(clinical_all$OS.time > 0, na.rm = TRUE), "\n")


# ------------------------------------------------------------
# 5. OS.status detailed counts
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("OS.status detailed counts\n")
cat("------------------------------------------------------------\n")

cat("Raw values:\n")
print(
  table(
    clinical_all$OS.status,
    useNA = "ifany"
  )
)

cat("\nNormalized values:\n")
print(
  table(
    os_status_chr,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 6. Unexpected status values
# ------------------------------------------------------------

expected_status <- c(
  "alive",
  "dead"
)

unexpected_status <- setdiff(
  unique(
    os_status_chr[
      !is.na(os_status_chr)
    ]
  ),
  expected_status
)

cat("\nUnexpected OS.status values:\n")

if (length(unexpected_status) == 0) {
  cat("NONE\n")
} else {
  print(unexpected_status)
}


# ------------------------------------------------------------
# 7. Problematic OS.time rows
# ------------------------------------------------------------

os_time_problem <- (
  is.na(clinical_all$OS.time) |
    is.nan(clinical_all$OS.time) |
    !is.finite(clinical_all$OS.time) |
    clinical_all$OS.time < 0
)

cat("\nNumber of problematic OS.time rows: ",
    sum(os_time_problem),
    "\n")

if (sum(os_time_problem) > 0) {
  
  cat("\nProblematic OS.time records:\n")
  
  print(
    clinical_all[
      os_time_problem,
      c(
        "Sample",
        "patient_id",
        "OS.time",
        "OS.status",
        "pathologic_stage"
      ),
      drop = FALSE
    ]
  )
  
}


# ------------------------------------------------------------
# 8. Problematic OS.status rows
# ------------------------------------------------------------

os_status_problem <- (
  is.na(os_status_chr) |
    !(os_status_chr %in% expected_status)
)

cat("\nNumber of problematic OS.status rows: ",
    sum(os_status_problem),
    "\n")

if (sum(os_status_problem) > 0) {
  
  cat("\nProblematic OS.status records:\n")
  
  print(
    clinical_all[
      os_status_problem,
      c(
        "Sample",
        "patient_id",
        "OS.time",
        "OS.status",
        "pathologic_stage"
      ),
      drop = FALSE
    ]
  )
  
}


# ============================================================
# STEP 8A-3B
# Survival-variable completeness audit
# NO DATA MODIFICATION
# ============================================================

cat("\n============================================================\n")
cat("STEP 8A-3B: SURVIVAL VARIABLE COMPLETENESS AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define valid OS.time
# ------------------------------------------------------------

os_time_valid <- (
  !is.na(clinical_all$OS.time) &
    is.finite(clinical_all$OS.time) &
    clinical_all$OS.time >= 0
)


# ------------------------------------------------------------
# 2. Define valid OS.status
# ------------------------------------------------------------

os_status_valid <- (
  !is.na(os_status_chr) &
    os_status_chr %in% c("alive", "dead")
)


# ------------------------------------------------------------
# 3. Combined survival completeness
# ------------------------------------------------------------

os_complete <- (
  os_time_valid &
    os_status_valid
)


# ------------------------------------------------------------
# 4. Summary
# ------------------------------------------------------------

cat("Total clinical samples: ",
    nrow(clinical_all), "\n")

cat("Valid OS.time:           ",
    sum(os_time_valid), "\n")

cat("Missing/invalid OS.time: ",
    sum(!os_time_valid), "\n")

cat("Valid OS.status:         ",
    sum(os_status_valid), "\n")

cat("Missing/invalid OS.status:",
    sum(!os_status_valid), "\n")

cat("Complete OS.time + OS.status:",
    sum(os_complete), "\n")


# ------------------------------------------------------------
# 5. Cross-tab of survival completeness
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("OS.time validity × OS.status validity\n")
cat("------------------------------------------------------------\n\n")

print(
  table(
    OS_time_valid = os_time_valid,
    OS_status_valid = os_status_valid,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 6. Samples excluded from complete survival analysis
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Samples without complete survival information\n")
cat("------------------------------------------------------------\n\n")

survival_incomplete <- !os_complete

print(
  clinical_all[
    survival_incomplete,
    c(
      "Sample",
      "patient_id",
      "OS.time",
      "OS.status",
      "pathologic_stage"
    ),
    drop = FALSE
  ]
)


# ============================================================
# PAPER 1 — WORKFLOW-08
# STEP 8A-4
# Build analysis-ready survival cohort
#
# Purpose:
#   1. Do NOT modify original clinical input
#   2. Create complete-case survival cohort
#   3. Preserve original sample order
#   4. Create explicit OS event variable
#   5. Verify expression/clinical sample alignment
#   6. Save a new analysis-ready clinical object
#
# Original input:
#   Workflow-04-2A_clinical_aligned.txt
#
# Expected survival cohort:
#   N = 503
# ============================================================

cat("\n============================================================\n")
cat("STEP 8A-4: ANALYSIS-READY SURVIVAL COHORT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Safety check
# ------------------------------------------------------------

stopifnot(exists("clinical_all"))
stopifnot(exists("expr_all"))

cat("clinical_all exists: TRUE\n")
cat("expr_all exists: TRUE\n")


# ------------------------------------------------------------
# 1. Reconstruct / verify survival validity flags
# ------------------------------------------------------------

os_status_chr <- trimws(
  tolower(
    as.character(clinical_all$OS.status)
  )
)

os_time_valid <- (
  !is.na(clinical_all$OS.time) &
    is.finite(clinical_all$OS.time) &
    clinical_all$OS.time >= 0
)

os_status_valid <- (
  !is.na(os_status_chr) &
    os_status_chr %in% c("alive", "dead")
)

os_complete <- (
  os_time_valid &
    os_status_valid
)


# ------------------------------------------------------------
# 2. Create analysis-ready survival cohort
#
# IMPORTANT:
#   Original clinical_all is NOT modified.
# ------------------------------------------------------------

clinical_survival_503 <- clinical_all[
  os_complete,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 3. Add normalized OS.status
# ------------------------------------------------------------

clinical_survival_503$OS.status_normalized <- os_status_chr[
  os_complete
]


# ------------------------------------------------------------
# 4. Create explicit OS event variable
#
# alive = 0
# dead  = 1
# ------------------------------------------------------------

clinical_survival_503$OS.event <- ifelse(
  clinical_survival_503$OS.status_normalized == "dead",
  1L,
  0L
)


# ------------------------------------------------------------
# 5. Preserve original sample order
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Cohort size\n")
cat("------------------------------------------------------------\n\n")

cat("Original clinical cohort: ",
    nrow(clinical_all), "\n")

cat("Analysis-ready survival cohort: ",
    nrow(clinical_survival_503), "\n")

cat("Excluded from survival analyses: ",
    nrow(clinical_all) - nrow(clinical_survival_503), "\n")


# ------------------------------------------------------------
# 6. Verify expected cohort size
# ------------------------------------------------------------

expected_survival_n <- 503

survival_n_pass <- (
  nrow(clinical_survival_503) == expected_survival_n
)

cat("\nExpected survival cohort N = ",
    expected_survival_n, "\n")

cat("Observed survival cohort N = ",
    nrow(clinical_survival_503), "\n")

cat("Survival cohort size check: ",
    survival_n_pass, "\n")


# ------------------------------------------------------------
# 7. Verify OS.time
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("OS.time QC in analysis-ready cohort\n")
cat("------------------------------------------------------------\n\n")

cat("Missing OS.time: ",
    sum(is.na(clinical_survival_503$OS.time)), "\n")

cat("Non-finite OS.time: ",
    sum(!is.finite(clinical_survival_503$OS.time)), "\n")

cat("Negative OS.time: ",
    sum(clinical_survival_503$OS.time < 0, na.rm = TRUE), "\n")

cat("Zero OS.time: ",
    sum(clinical_survival_503$OS.time == 0, na.rm = TRUE), "\n")

cat("Positive OS.time: ",
    sum(clinical_survival_503$OS.time > 0, na.rm = TRUE), "\n")


# ------------------------------------------------------------
# 8. Verify OS.status
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("OS.status QC in analysis-ready cohort\n")
cat("------------------------------------------------------------\n\n")

print(
  table(
    clinical_survival_503$OS.status_normalized,
    useNA = "ifany"
  )
)

cat("\nOS.event distribution:\n")

print(
  table(
    clinical_survival_503$OS.event,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 9. Verify Sample IDs
# ------------------------------------------------------------

survival_sample_ids <- clinical_survival_503$Sample
expression_sample_ids <- colnames(expr_all)[-1]


# ------------------------------------------------------------
# 10. Check survival sample IDs are unique
# ------------------------------------------------------------

survival_id_unique <- (
  !anyDuplicated(survival_sample_ids)
)

cat("\n------------------------------------------------------------\n")
cat("Survival sample ID integrity\n")
cat("------------------------------------------------------------\n\n")

cat("Number of survival samples: ",
    length(survival_sample_ids), "\n")

cat("Missing survival sample IDs: ",
    sum(is.na(survival_sample_ids) |
          trimws(survival_sample_ids) == ""), "\n")

cat("Duplicated survival sample IDs: ",
    sum(duplicated(survival_sample_ids)), "\n")

cat("Unique survival sample IDs: ",
    survival_id_unique, "\n")


# ------------------------------------------------------------
# 11. Check survival cohort against expression cohort
#
# IMPORTANT:
#   Expression matrix contains all 513 tumor samples.
#   Therefore every survival sample must exist in expression.
# ------------------------------------------------------------

survival_not_in_expression <- setdiff(
  survival_sample_ids,
  expression_sample_ids
)

expression_not_in_survival <- setdiff(
  expression_sample_ids,
  survival_sample_ids
)

cat("\n------------------------------------------------------------\n")
cat("Expression ↔ survival sample-ID check\n")
cat("------------------------------------------------------------\n\n")

cat("Survival samples not found in expression: ",
    length(survival_not_in_expression), "\n")

cat("Expression samples not included in survival cohort: ",
    length(expression_not_in_survival), "\n")


# ------------------------------------------------------------
# 12. Check that survival cohort is a subset of expression cohort
# ------------------------------------------------------------

survival_subset_expression <- (
  length(survival_not_in_expression) == 0
)

cat("\nSurvival cohort is subset of expression cohort: ",
    survival_subset_expression, "\n")


# ------------------------------------------------------------
# 13. Check sample order
#
# The 503 survival samples should occur in the same
# relative order as in the original expression matrix.
# ------------------------------------------------------------

expression_positions <- match(
  survival_sample_ids,
  expression_sample_ids
)

survival_order_correct <- (
  all(diff(expression_positions) > 0)
)

cat("Survival sample order follows expression order: ",
    survival_order_correct, "\n")


# ------------------------------------------------------------
# 14. Check tumor status
# ------------------------------------------------------------

survival_tumor_pass <- all(
  !is.na(clinical_survival_503$tumor_normal) &
    clinical_survival_503$tumor_normal == "tumor"
)

cat("All survival samples are tumor: ",
    survival_tumor_pass, "\n")


# ------------------------------------------------------------
# 15. Check patient IDs
# ------------------------------------------------------------

patient_id_missing <- sum(
  is.na(clinical_survival_503$patient_id) |
    trimws(as.character(
      clinical_survival_503$patient_id
    )) == ""
)

patient_id_duplicate <- sum(
  duplicated(clinical_survival_503$patient_id)
)

cat("\n------------------------------------------------------------\n")
cat("Patient ID QC\n")
cat("------------------------------------------------------------\n\n")

cat("Missing patient IDs: ",
    patient_id_missing, "\n")

cat("Duplicated patient IDs: ",
    patient_id_duplicate, "\n")


# ------------------------------------------------------------
# 16. Final survival data integrity check
# ------------------------------------------------------------

final_survival_qc <- (
  survival_n_pass &
    sum(is.na(clinical_survival_503$OS.time)) == 0 &
    sum(!is.finite(clinical_survival_503$OS.time)) == 0 &
    sum(clinical_survival_503$OS.time < 0) == 0 &
    sum(is.na(clinical_survival_503$OS.status_normalized)) == 0 &
    all(
      clinical_survival_503$OS.status_normalized %in%
        c("alive", "dead")
    ) &
    sum(is.na(clinical_survival_503$OS.event)) == 0 &
    all(
      clinical_survival_503$OS.event %in% c(0L, 1L)
    ) &
    survival_id_unique &
    survival_subset_expression &
    survival_order_correct &
    survival_tumor_pass &
    patient_id_missing == 0
)


# ------------------------------------------------------------
# 17. Final report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 8A-4 FINAL QC\n")
cat("============================================================\n\n")

cat("Original clinical cohort:          ",
    nrow(clinical_all), "\n")

cat("Analysis-ready survival cohort:    ",
    nrow(clinical_survival_503), "\n")

cat("Complete OS.time + OS.status:      ",
    sum(os_complete), "\n")

cat("OS.event coding valid:             ",
    all(clinical_survival_503$OS.event %in% c(0L, 1L)),
    "\n")

cat("Survival IDs unique:               ",
    survival_id_unique, "\n")

cat("Survival IDs all in expression:    ",
    survival_subset_expression, "\n")

cat("Survival order aligned:             ",
    survival_order_correct, "\n")

cat("All survival samples are tumor:     ",
    survival_tumor_pass, "\n")

cat("Patient IDs complete:               ",
    patient_id_missing == 0, "\n")

cat("\n------------------------------------------------------------\n")

if (final_survival_qc) {
  
  cat("STEP 8A-4 OVERALL RESULT: PASS\n")
  
} else {
  
  cat("STEP 8A-4 OVERALL RESULT: FAIL\n")
  
}

cat("------------------------------------------------------------\n\n")


# ------------------------------------------------------------
# 18. Inspect first rows of final cohort
# ------------------------------------------------------------

cat("First 6 rows of analysis-ready survival cohort:\n\n")

print(
  clinical_survival_503[
    1:min(6, nrow(clinical_survival_503)),
    c(
      "Sample",
      "patient_id",
      "OS.time",
      "OS.status",
      "OS.status_normalized",
      "OS.event",
      "pathologic_stage"
    ),
    drop = FALSE
  ]
)


# ------------------------------------------------------------
# 19. Inspect OS event distribution
# ------------------------------------------------------------

cat("\nOS.event summary:\n")

print(
  table(
    clinical_survival_503$OS.event
  )
)


# ------------------------------------------------------------
# 20. Save analysis-ready object
#
# This creates a NEW file.
# The original Workflow-04-2A clinical file is untouched.
# ------------------------------------------------------------

saveRDS(
  clinical_survival_503,
  file = "result_dir/Paper1_tumor-only WGCNA/Workflow-08-STEP8A4_clinical_survival_503.rds"
)

cat("\nAnalysis-ready survival cohort saved as:\n")

cat(
  "result_dir/Paper1_tumor-only WGCNA/",
  "Workflow-08-STEP8A4_clinical_survival_503.rds\n",
  sep = ""
)

cat("\n============================================================\n")
cat("STEP 8A-4 COMPLETED\n")
cat("============================================================\n")


# ============================================================
# PAPER 1 — WORKFLOW-08
# STEP 8B
# Final expression–survival cohort alignment
#
# Purpose:
#   1. Load / verify the locked 503-sample survival cohort
#   2. Match the 503 clinical samples to the expression matrix
#   3. Verify exact sample-ID correspondence
#   4. Verify expression sample order after alignment
#   5. Create an analysis-ready 503-sample expression matrix
#   6. Do NOT modify the original 513-sample expression matrix
#
# INPUTS:
#   Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt
#   Workflow-08-STEP8A4_clinical_survival_503.rds
#
# OUTPUT:
#   Workflow-08-STEP8B_expr_survival_503.rds
#
# Expected:
#   39,280 genes × 503 samples
#   + 1 gene-ID column if exported as table
# ============================================================

cat("\n============================================================\n")
cat("STEP 8B: FINAL EXPRESSION–SURVIVAL COHORT ALIGNMENT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Safety checks
# ------------------------------------------------------------

stopifnot(exists("expr_all"))

cat("expr_all exists: TRUE\n")


# ------------------------------------------------------------
# 1. Load locked 503-sample survival cohort
# ------------------------------------------------------------

survival_file <- paste0(
  "input_dir/Paper1_tumor-only WGCNA/",
  "Workflow-08-STEP8A4_clinical_survival_503.rds"
)

if (!file.exists(survival_file)) {
  stop(
    "STEP 8B STOP: Survival cohort file not found:\n",
    survival_file
  )
}

clinical_survival_503 <- readRDS(
  survival_file
)

cat("Loaded analysis-ready survival cohort:\n")
cat(survival_file, "\n")


# ------------------------------------------------------------
# 2. Basic cohort verification
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("2. Basic cohort verification\n")
cat("------------------------------------------------------------\n\n")

cat("Clinical survival cohort dimensions: ",
    nrow(clinical_survival_503),
    " rows × ",
    ncol(clinical_survival_503),
    " columns\n",
    sep = "")

cat("Expression matrix dimensions: ",
    nrow(expr_all),
    " genes × ",
    ncol(expr_all) - 1,
    " samples\n",
    sep = "")


# ------------------------------------------------------------
# 3. Verify expected sample number
# ------------------------------------------------------------

expected_n_survival <- 503

if (nrow(clinical_survival_503) != expected_n_survival) {
  stop(
    "STEP 8B STOP: Expected 503 survival samples, but found ",
    nrow(clinical_survival_503)
  )
}

cat("\nExpected survival cohort N = 503: PASS\n")


# ------------------------------------------------------------
# 4. Extract sample IDs
# ------------------------------------------------------------

clinical_ids <- as.character(
  clinical_survival_503$Sample
)

expression_ids <- colnames(
  expr_all
)[-1]


# ------------------------------------------------------------
# 5. Basic sample-ID QC
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("5. Sample-ID QC\n")
cat("------------------------------------------------------------\n\n")

clinical_missing_ids <- sum(
  is.na(clinical_ids) |
    trimws(clinical_ids) == ""
)

expression_missing_ids <- sum(
  is.na(expression_ids) |
    trimws(expression_ids) == ""
)

clinical_duplicate_ids <- sum(
  duplicated(clinical_ids)
)

expression_duplicate_ids <- sum(
  duplicated(expression_ids)
)

cat("Clinical missing IDs:   ",
    clinical_missing_ids, "\n")

cat("Clinical duplicated IDs:",
    clinical_duplicate_ids, "\n")

cat("Expression missing IDs:",
    expression_missing_ids, "\n")

cat("Expression duplicated IDs:",
    expression_duplicate_ids, "\n")


# ------------------------------------------------------------
# 6. Exact set comparison
# ------------------------------------------------------------

clinical_not_in_expression <- setdiff(
  clinical_ids,
  expression_ids
)

expression_not_in_clinical <- setdiff(
  expression_ids,
  clinical_ids
)

cat("\n------------------------------------------------------------\n")
cat("6. Exact set comparison\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Clinical IDs not in expression: ",
  length(clinical_not_in_expression),
  "\n",
  sep = ""
)

cat(
  "Expression IDs not in clinical survival cohort: ",
  length(expression_not_in_clinical),
  "\n",
  sep = ""
)

if (length(clinical_not_in_expression) > 0) {
  cat("\nClinical IDs not in expression:\n")
  print(clinical_not_in_expression)
}

if (length(expression_not_in_clinical) > 0) {
  cat("\nExpression IDs not in clinical survival cohort:\n")
  print(expression_not_in_clinical)
}

cat("\nExact clinical ↔ expression ID set comparison: ",
    length(clinical_not_in_expression) == 0 &&
      length(expression_not_in_clinical) == 10,
    "\n",
    sep = "")


# ------------------------------------------------------------
# 7. Match clinical survival IDs to expression IDs
# ------------------------------------------------------------

expression_position <- match(
  clinical_ids,
  expression_ids
)

match_complete <- !anyNA(expression_position)

cat("\n------------------------------------------------------------\n")
cat("7. Clinical → expression matching\n")
cat("------------------------------------------------------------\n\n")

cat(
  "All 503 clinical IDs matched to expression: ",
  match_complete,
  "\n",
  sep = ""
)

if (!match_complete) {
  stop(
    "STEP 8B STOP: One or more clinical survival IDs could not be matched."
  )
}

# ------------------------------------------------------------
# 8. Verify matched order
# ------------------------------------------------------------

matched_expression_ids <- expression_ids[
  expression_position
]

exact_order_match <- identical(
  clinical_ids,
  matched_expression_ids
)

cat(
  "Clinical ID order equals matched expression order: ",
  exact_order_match,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 9. Build aligned 503-sample expression matrix
# ------------------------------------------------------------

expr_survival_503 <- expr_all[
  ,
  c(1, expression_position + 1),
  drop = FALSE
]

cat("\nAligned expression matrix dimensions: ",
    nrow(expr_survival_503),
    " genes × ",
    ncol(expr_survival_503) - 1,
    " samples\n",
    sep = "")

# ------------------------------------------------------------
# 10. Final STEP 8B alignment QC
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("10. Final STEP 8B alignment QC\n")
cat("------------------------------------------------------------\n\n")

# ------------------------------------------------------------
# 10.1 Verify aligned sample IDs
# ------------------------------------------------------------

aligned_expression_ids <- colnames(expr_survival_503)[-1]

aligned_ids_identical <- identical(
  aligned_expression_ids,
  clinical_ids
)

cat(
  "Aligned expression IDs identical to clinical IDs: ",
  aligned_ids_identical,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10.2 Verify gene number preservation
# ------------------------------------------------------------

gene_number_pass <- identical(
  nrow(expr_survival_503),
  nrow(expr_all)
)

cat(
  "Gene number preserved: ",
  gene_number_pass,
  " (",
  nrow(expr_survival_503),
  " genes)\n",
  sep = ""
)


# ------------------------------------------------------------
# 10.3 Verify gene IDs are unchanged
# ------------------------------------------------------------

gene_ids_original <- as.character(
  expr_all[[1]]
)

gene_ids_aligned <- as.character(
  expr_survival_503[[1]]
)

gene_ids_identical <- identical(
  gene_ids_original,
  gene_ids_aligned
)

cat(
  "Gene IDs preserved exactly: ",
  gene_ids_identical,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10.4 Verify expression values are preserved
# ------------------------------------------------------------

expression_value_check <- TRUE

for (i in seq_along(expression_position)) {
  
  original_col <- expression_position[i] + 1
  aligned_col <- i + 1
  
  if (!identical(
    expr_all[[original_col]],
    expr_survival_503[[aligned_col]]
  )) {
    
    expression_value_check <- FALSE
    
    cat(
      "Expression value mismatch detected at aligned sample: ",
      clinical_ids[i],
      "\n",
      sep = ""
    )
    
    break
  }
}

cat(
  "Expression values preserved exactly: ",
  expression_value_check,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10.5 Verify row-by-row clinical/expression alignment
# ------------------------------------------------------------

row_by_row_alignment <- identical(
  clinical_survival_503$Sample,
  colnames(expr_survival_503)[-1]
)

cat(
  "Row-by-row clinical ↔ expression alignment: ",
  row_by_row_alignment,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10.6 Verify survival variables
# ------------------------------------------------------------

survival_os_complete <- (
  sum(is.na(clinical_survival_503$OS.time)) == 0 &&
    sum(!is.finite(clinical_survival_503$OS.time)) == 0 &&
    sum(clinical_survival_503$OS.time < 0) == 0 &&
    sum(is.na(clinical_survival_503$OS.event)) == 0 &&
    all(clinical_survival_503$OS.event %in% c(0L, 1L))
)

cat(
  "Survival variables complete and valid: ",
  survival_os_complete,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10.7 Verify tumor-only status
# ------------------------------------------------------------

tumor_status_pass <- all(
  !is.na(clinical_survival_503$tumor_normal) &
    clinical_survival_503$tumor_normal == "tumor"
)

cat(
  "All survival samples are tumor samples: ",
  tumor_status_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10.8 Verify unique sample IDs
# ------------------------------------------------------------

aligned_sample_unique <- (
  length(aligned_expression_ids) == 503 &&
    !anyDuplicated(aligned_expression_ids)
)

clinical_sample_unique <- (
  length(clinical_ids) == 503 &&
    !anyDuplicated(clinical_ids)
)

cat(
  "Aligned expression sample IDs unique: ",
  aligned_sample_unique,
  "\n",
  sep = ""
)

cat(
  "Clinical survival sample IDs unique: ",
  clinical_sample_unique,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Final STEP 8B decision
# ------------------------------------------------------------

step8b_pass <- (
  nrow(clinical_survival_503) == 503 &&
    ncol(expr_survival_503) - 1 == 503 &&
    length(clinical_not_in_expression) == 0 &&
    length(expression_not_in_clinical) == 10 &&
    match_complete &&
    exact_order_match &&
    aligned_ids_identical &&
    gene_number_pass &&
    gene_ids_identical &&
    expression_value_check &&
    row_by_row_alignment &&
    survival_os_complete &&
    tumor_status_pass &&
    aligned_sample_unique &&
    clinical_sample_unique
)

cat("\n============================================================\n")
cat("STEP 8B FINAL RESULT\n")
cat("============================================================\n\n")

cat(
  "Clinical survival cohort: ",
  nrow(clinical_survival_503),
  " samples\n",
  sep = ""
)

cat(
  "Original expression cohort: ",
  ncol(expr_all) - 1,
  " samples\n",
  sep = ""
)

cat(
  "Aligned expression cohort: ",
  ncol(expr_survival_503) - 1,
  " samples\n",
  sep = ""
)

cat(
  "Expression genes preserved: ",
  nrow(expr_survival_503),
  "\n",
  sep = ""
)

cat(
  "Clinical → expression IDs all matched: ",
  match_complete,
  "\n",
  sep = ""
)

cat(
  "Aligned IDs identical: ",
  aligned_ids_identical,
  "\n",
  sep = ""
)

cat(
  "Expression values preserved: ",
  expression_value_check,
  "\n",
  sep = ""
)

cat(
  "Row-by-row alignment: ",
  row_by_row_alignment,
  "\n",
  sep = ""
)

cat(
  "Survival variables complete: ",
  survival_os_complete,
  "\n",
  sep = ""
)

cat(
  "All samples tumor: ",
  tumor_status_pass,
  "\n",
  sep = ""
)

cat("\n------------------------------------------------------------\n")

if (step8b_pass) {
  
  cat("STEP 8B OVERALL RESULT: PASS\n")
  
} else {
  
  cat("STEP 8B OVERALL RESULT: FAIL\n")
  
}

cat("------------------------------------------------------------\n")


# ------------------------------------------------------------
# STEP 8B-12
# Save locked 503-sample expression cohort
# ------------------------------------------------------------

step8b_expression_file <- paste0(
  "result_dir/Paper1_tumor-only WGCNA/",
  "Workflow-08-STEP8B_expr_survival_503.rds"
)

saveRDS(
  expr_survival_503,
  file = step8b_expression_file
)

cat("\n============================================================\n")
cat("STEP 8B SAVE\n")
cat("============================================================\n\n")

cat(
  "Saved file:\n",
  step8b_expression_file,
  "\n\n",
  sep = ""
)

cat(
  "File exists: ",
  file.exists(step8b_expression_file),
  "\n",
  sep = ""
)

if (file.exists(step8b_expression_file)) {
  cat("STEP 8B SAVE RESULT: PASS\n")
} else {
  cat("STEP 8B SAVE RESULT: FAIL\n")
}


# ============================================================
# Paper 1 — WORKFLOW-08
# STEP 8C-0
# Gene ID namespace audit
# ============================================================

cat("\n============================================================\n")
cat("STEP 8C-0: Gene ID namespace audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Safety checks
# ------------------------------------------------------------

stopifnot(exists("expr_survival_503"))

gene_ids_raw <- as.character(
  expr_survival_503[[1]]
)

cat(
  "Total gene IDs: ",
  length(gene_ids_raw),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 1. Show first 30 gene IDs
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("1. First 30 gene IDs\n")
cat("------------------------------------------------------------\n\n")

print(
  head(gene_ids_raw, 30)
)


# ------------------------------------------------------------
# 2. Show last 30 gene IDs
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("2. Last 30 gene IDs\n")
cat("------------------------------------------------------------\n\n")

print(
  tail(gene_ids_raw, 30)
)


# ------------------------------------------------------------
# 3. Basic formatting characteristics
# ------------------------------------------------------------

gene_id_has_na <- any(
  is.na(gene_ids_raw)
)

gene_id_has_blank <- any(
  trimws(gene_ids_raw) == ""
)

gene_id_has_version <- any(
  grepl(
    "^ENSG[0-9]+\\.[0-9]+$",
    gene_ids_raw
  )
)

gene_id_ensembl_like <- sum(
  grepl(
    "^ENSG[0-9]+",
    gene_ids_raw
  )
)

gene_id_entrez_like <- sum(
  grepl(
    "^[0-9]+$",
    gene_ids_raw
  )
)

gene_id_refseq_like <- sum(
  grepl(
    "^(NM_|NR_|XM_|XR_)",
    gene_ids_raw
  )
)

cat("\n------------------------------------------------------------\n")
cat("3. Gene ID pattern audit\n")
cat("------------------------------------------------------------\n\n")

cat(
  "NA gene IDs: ",
  sum(is.na(gene_ids_raw)),
  "\n",
  sep = ""
)

cat(
  "Blank gene IDs: ",
  sum(trimws(gene_ids_raw) == ""),
  "\n",
  sep = ""
)

cat(
  "Ensembl-like IDs: ",
  gene_id_ensembl_like,
  "\n",
  sep = ""
)

cat(
  "Ensembl versioned IDs: ",
  sum(gene_id_has_version),
  "\n",
  sep = ""
)

cat(
  "Entrez-like numeric IDs: ",
  gene_id_entrez_like,
  "\n",
  sep = ""
)

cat(
  "RefSeq-like IDs: ",
  gene_id_refseq_like,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Check whether any standard gene symbols are present
# ------------------------------------------------------------

known_symbols <- c(
  "AOC3",
  "MS4A7",
  "IL7R",
  "PDCD1",
  "CD274",
  "CTLA4",
  "LAG3"
)

symbol_presence <- known_symbols %in% gene_ids_raw

cat("\n------------------------------------------------------------\n")
cat("4. Known gene-symbol presence check\n")
cat("------------------------------------------------------------\n\n")

symbol_check_table <- data.frame(
  Gene = known_symbols,
  Present = symbol_presence,
  stringsAsFactors = FALSE
)

print(
  symbol_check_table,
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Basic conclusion
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 8C-0 PRELIMINARY RESULT\n")
cat("============================================================\n\n")

if (
  gene_id_ensembl_like > 0 &&
  sum(gene_id_has_version) > 0
) {
  
  cat(
    "The expression matrix appears to use versioned Ensembl-like gene IDs.\n"
  )
  
} else if (
  gene_id_ensembl_like > 0
) {
  
  cat(
    "The expression matrix appears to use Ensembl-like gene IDs.\n"
  )
  
} else if (
  gene_id_entrez_like > 0
) {
  
  cat(
    "The expression matrix appears to use numeric/Entrez-like gene IDs.\n"
  )
  
} else if (
  gene_id_refseq_like > 0
) {
  
  cat(
    "The expression matrix appears to use RefSeq-like gene IDs.\n"
  )
  
} else {
  
  cat(
    "Gene ID namespace is not yet identifiable from the tested patterns.\n"
  )
}

# ============================================================
# Paper 1 — WORKFLOW-08
# STEP 8C-0-2
# Check available human gene annotation resources
# ============================================================

cat("\n============================================================\n")
cat("STEP 8C-0-2: Annotation resource check\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Check whether org.Hs.eg.db is installed
# ------------------------------------------------------------

org_hs_available <- requireNamespace(
  "org.Hs.eg.db",
  quietly = TRUE
)

cat(
  "org.Hs.eg.db available: ",
  org_hs_available,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 2. Check whether AnnotationDbi is installed
# ------------------------------------------------------------

annotation_dbi_available <- requireNamespace(
  "AnnotationDbi",
  quietly = TRUE
)

cat(
  "AnnotationDbi available: ",
  annotation_dbi_available,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 3. Check whether biomaRt is installed
# ------------------------------------------------------------

biomart_available <- requireNamespace(
  "biomaRt",
  quietly = TRUE
)

cat(
  "biomaRt available: ",
  biomart_available,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 4. Current Ensembl ID integrity
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Current Ensembl ID integrity\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Total IDs: ",
  length(gene_ids_raw),
  "\n",
  sep = ""
)

cat(
  "Unique IDs: ",
  length(unique(gene_ids_raw)),
  "\n",
  sep = ""
)

cat(
  "Duplicated IDs: ",
  sum(duplicated(gene_ids_raw)),
  "\n",
  sep = ""
)

cat(
  "All IDs match Ensembl pattern: ",
  all(grepl("^ENSG[0-9]+$", gene_ids_raw)),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 5. Safety conclusion
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 8C-0-2 PRELIMINARY RESULT\n")
cat("============================================================\n\n")

if (org_hs_available && annotation_dbi_available) {
  
  cat(
    "Recommended route: use org.Hs.eg.db / AnnotationDbi for\n",
    "Ensembl-to-symbol mapping.\n",
    sep = ""
  )
  
} else if (biomart_available) {
  
  cat(
    "org.Hs.eg.db is not fully available, but biomaRt is available.\n",
    "Annotation mapping can be considered using biomaRt.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "No suitable local annotation resource detected.\n",
    "Do not perform manual ID replacement.\n",
    sep = ""
  )
}


# ============================================================
# Paper 1 — WORKFLOW-08
# STEP 8C-0-3
# Build auditable Ensembl -> Gene Symbol annotation
# ============================================================

cat("\n============================================================\n")
cat("STEP 8C-0-3: Ensembl-to-Symbol annotation audit\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Safety checks
# ------------------------------------------------------------

stopifnot(exists("expr_survival_503"))
stopifnot(exists("gene_ids_raw"))

stopifnot(
  requireNamespace("org.Hs.eg.db", quietly = TRUE)
)

stopifnot(
  requireNamespace("AnnotationDbi", quietly = TRUE)
)

# ------------------------------------------------------------
# 1. Load annotation database
# ------------------------------------------------------------

library(
  AnnotationDbi
)

library(
  org.Hs.eg.db
)

# ------------------------------------------------------------
# 2. Prepare unique Ensembl IDs
# ------------------------------------------------------------

ensembl_ids <- unique(
  as.character(gene_ids_raw)
)

cat(
  "Input Ensembl IDs: ",
  length(ensembl_ids),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 3. Perform Ensembl -> Symbol mapping
# ------------------------------------------------------------

annotation_raw <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_ids,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)

# ------------------------------------------------------------
# 4. Basic annotation output inspection
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("4. Raw annotation result\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Rows returned: ",
  nrow(annotation_raw),
  "\n",
  sep = ""
)

cat(
  "Columns returned: ",
  paste(
    colnames(annotation_raw),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat("\nFirst 20 annotation rows:\n")
print(
  head(annotation_raw, 20),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Mapping coverage
# ------------------------------------------------------------

mapping_has_symbol <- !is.na(
  annotation_raw$SYMBOL
) &
  trimws(
    annotation_raw$SYMBOL
  ) != ""

mapped_ensembl_ids <- unique(
  annotation_raw$ENSEMBL[
    mapping_has_symbol
  ]
)

unmapped_ensembl_ids <- setdiff(
  ensembl_ids,
  mapped_ensembl_ids
)

cat("\n------------------------------------------------------------\n")
cat("5. Mapping coverage\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Total input Ensembl IDs: ",
  length(ensembl_ids),
  "\n",
  sep = ""
)

cat(
  "Ensembl IDs with SYMBOL mapping: ",
  length(mapped_ensembl_ids),
  "\n",
  sep = ""
)

cat(
  "Ensembl IDs without SYMBOL mapping: ",
  length(unmapped_ensembl_ids),
  "\n",
  sep = ""
)

cat(
  "Mapping coverage: ",
  round(
    100 *
      length(mapped_ensembl_ids) /
      length(ensembl_ids),
    2
  ),
  "%\n",
  sep = ""
)

# ------------------------------------------------------------
# 6. Check one-to-many Ensembl -> SYMBOL mappings
# ------------------------------------------------------------

mapping_per_ensembl <- aggregate(
  SYMBOL ~ ENSEMBL,
  data = annotation_raw[
    mapping_has_symbol,
    c("ENSEMBL", "SYMBOL")
  ],
  FUN = function(x) {
    unique(
      x[
        !is.na(x) &
          trimws(x) != ""
      ]
    )
  }
)

one_to_many_ensembl <- mapping_per_ensembl[
  lengths(
    mapping_per_ensembl$SYMBOL
  ) > 1,
  ,
  drop = FALSE
]

cat("\n------------------------------------------------------------\n")
cat("6. One-to-many Ensembl -> SYMBOL mappings\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Ensembl IDs mapping to >1 SYMBOL: ",
  nrow(one_to_many_ensembl),
  "\n",
  sep = ""
)

if (nrow(one_to_many_ensembl) > 0) {
  
  print(
    head(
      one_to_many_ensembl,
      20
    ),
    row.names = FALSE
  )
  
}

# ------------------------------------------------------------
# 7. Check many-to-one SYMBOL mappings
# ------------------------------------------------------------

symbol_to_ensembl <- split(
  annotation_raw$ENSEMBL[
    mapping_has_symbol
  ],
  annotation_raw$SYMBOL[
    mapping_has_symbol
  ]
)

many_to_one_symbol <- symbol_to_ensembl[
  lengths(symbol_to_ensembl) > 1
]

cat("\n------------------------------------------------------------\n")
cat("7. Many-to-one SYMBOL mappings\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Gene symbols associated with >1 Ensembl ID: ",
  length(many_to_one_symbol),
  "\n",
  sep = ""
)

# Show a small audit sample
if (length(many_to_one_symbol) > 0) {
  
  example_symbols <- head(
    many_to_one_symbol,
    20
  )
  
  example_table <- data.frame(
    SYMBOL = names(example_symbols),
    ENSEMBL_IDs = vapply(
      example_symbols,
      function(x) paste(
        unique(x),
        collapse = "; "
      ),
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  
  print(
    example_table,
    row.names = FALSE
  )
  
}

# ------------------------------------------------------------
# 8. Explicit audit of Paper 1 core genes
# ------------------------------------------------------------

core_genes_8c0 <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_annotation <- annotation_raw[
  !is.na(annotation_raw$SYMBOL) &
    annotation_raw$SYMBOL %in% core_genes_8c0,
  c("ENSEMBL", "SYMBOL"),
  drop = FALSE
]

core_annotation <- unique(
  core_annotation
)

cat("\n------------------------------------------------------------\n")
cat("8. Paper 1 core gene annotation audit\n")
cat("------------------------------------------------------------\n\n")

print(
  core_annotation,
  row.names = FALSE
)

cat(
  "\nAOC3 mapped: ",
  "AOC3" %in% core_annotation$SYMBOL,
  "\n",
  sep = ""
)

cat(
  "MS4A7 mapped: ",
  "MS4A7" %in% core_annotation$SYMBOL,
  "\n",
  sep = ""
)

cat(
  "IL7R mapped: ",
  "IL7R" %in% core_annotation$SYMBOL,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 9. Final preliminary conclusion
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 8C-0-3 PRELIMINARY RESULT\n")
cat("============================================================\n\n")

if (
  nrow(core_annotation) > 0 &&
  length(mapped_ensembl_ids) > 0
) {
  
  cat(
    "Ensembl-to-SYMBOL annotation was successfully generated.\n"
  )
  
  cat(
    "The mapping will be retained as an independent annotation\n",
    "layer; the original expression matrix will not be modified.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "Annotation mapping requires further investigation.\n"
  )
}



# ============================================================
# Paper 1 — WORKFLOW-08
# STEP 8C-1
# Immune checkpoint panel — Ensembl mapping audit
# ============================================================

cat("\n============================================================\n")
cat("STEP 8C-1: Immune checkpoint panel mapping audit\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Safety checks
# ------------------------------------------------------------

stopifnot(exists("expr_survival_503"))
stopifnot(exists("annotation_raw"))
stopifnot(exists("gene_ids_raw"))

# ------------------------------------------------------------
# 1. Predefined checkpoint-associated gene panel
# ------------------------------------------------------------

checkpoint_panel <- c(
  "PDCD1",
  "CD274",
  "PDCD1LG2",
  "CTLA4",
  "CD80",
  "CD86",
  "LAG3",
  "HAVCR2",
  "TIGIT",
  "CD96",
  "VSIR",
  "VTCN1",
  "BTLA",
  "TNFRSF4",
  "TNFRSF9",
  "ICOS",
  "LGALS9",
  "CD200"
)

cat(
  "Predefined checkpoint-associated genes: ",
  length(checkpoint_panel),
  "\n",
  sep = ""
)

cat(
  "Panel duplicates: ",
  sum(duplicated(checkpoint_panel)),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 2. Extract panel mapping from annotation
# ------------------------------------------------------------

checkpoint_annotation <- annotation_raw[
  !is.na(annotation_raw$SYMBOL) &
    annotation_raw$SYMBOL %in% checkpoint_panel,
  c("ENSEMBL", "SYMBOL"),
  drop = FALSE
]

checkpoint_annotation <- unique(
  checkpoint_annotation
)

# ------------------------------------------------------------
# 3. Mapping summary by gene symbol
# ------------------------------------------------------------

checkpoint_mapping_summary <- do.call(
  rbind,
  lapply(
    checkpoint_panel,
    function(gene) {
      
      ids <- unique(
        checkpoint_annotation$ENSEMBL[
          checkpoint_annotation$SYMBOL == gene
        ]
      )
      
      data.frame(
        Gene = gene,
        Ensembl_IDs = if (
          length(ids) == 0
        ) {
          NA_character_
        } else {
          paste(
            ids,
            collapse = "; "
          )
        },
        N_Ensembl = length(ids),
        stringsAsFactors = FALSE
      )
    }
  )
)

checkpoint_mapping_summary$Mapping_Status <- ifelse(
  checkpoint_mapping_summary$N_Ensembl == 0,
  "NOT_MAPPED",
  ifelse(
    checkpoint_mapping_summary$N_Ensembl == 1,
    "UNIQUE",
    "MULTIPLE_ENSEMBL"
  )
)

# ------------------------------------------------------------
# 4. Print mapping audit
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("4. Checkpoint panel mapping audit\n")
cat("------------------------------------------------------------\n\n")

print(
  checkpoint_mapping_summary,
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Summary counts
# ------------------------------------------------------------

n_not_mapped <- sum(
  checkpoint_mapping_summary$Mapping_Status ==
    "NOT_MAPPED"
)

n_unique <- sum(
  checkpoint_mapping_summary$Mapping_Status ==
    "UNIQUE"
)

n_multiple <- sum(
  checkpoint_mapping_summary$Mapping_Status ==
    "MULTIPLE_ENSEMBL"
)

cat("\n------------------------------------------------------------\n")
cat("5. Mapping summary\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Panel genes uniquely mapped: ",
  n_unique,
  "\n",
  sep = ""
)

cat(
  "Panel genes with multiple Ensembl IDs: ",
  n_multiple,
  "\n",
  sep = ""
)

cat(
  "Panel genes not mapped: ",
  n_not_mapped,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 6. Check whether mapped Ensembl IDs are actually
#    present in the expression matrix
# ------------------------------------------------------------

expression_ensembl_ids <- gene_ids_raw

checkpoint_mapping_summary$Expression_Present <- vapply(
  seq_len(nrow(checkpoint_mapping_summary)),
  function(i) {
    
    ids <- checkpoint_annotation$ENSEMBL[
      checkpoint_annotation$SYMBOL ==
        checkpoint_mapping_summary$Gene[i]
    ]
    
    any(
      ids %in% expression_ensembl_ids
    )
  },
  logical(1)
)

# ------------------------------------------------------------
# 7. Final gene-level status
# ------------------------------------------------------------

checkpoint_mapping_summary$Final_Status <- vapply(
  seq_len(nrow(checkpoint_mapping_summary)),
  function(i) {
    
    gene <- checkpoint_mapping_summary$Gene[i]
    
    ids <- checkpoint_annotation$ENSEMBL[
      checkpoint_annotation$SYMBOL == gene
    ]
    
    ids_in_expression <- ids[
      ids %in% expression_ensembl_ids
    ]
    
    if (length(ids_in_expression) == 0) {
      
      "NOT_AVAILABLE_IN_EXPRESSION"
      
    } else if (length(ids_in_expression) == 1) {
      
      "AVAILABLE_UNIQUE"
      
    } else {
      
      "AVAILABLE_MULTIPLE"
      
    }
  },
  character(1)
)

cat("\n------------------------------------------------------------\n")
cat("7. Final expression availability audit\n")
cat("------------------------------------------------------------\n\n")

print(
  checkpoint_mapping_summary[
    ,
    c(
      "Gene",
      "Ensembl_IDs",
      "N_Ensembl",
      "Expression_Present",
      "Final_Status"
    )
  ],
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Final preliminary result
# ------------------------------------------------------------

n_available_unique <- sum(
  checkpoint_mapping_summary$Final_Status ==
    "AVAILABLE_UNIQUE"
)

n_available_multiple <- sum(
  checkpoint_mapping_summary$Final_Status ==
    "AVAILABLE_MULTIPLE"
)

n_not_available <- sum(
  checkpoint_mapping_summary$Final_Status ==
    "NOT_AVAILABLE_IN_EXPRESSION"
)

cat("\n============================================================\n")
cat("STEP 8C-1 PRELIMINARY RESULT\n")
cat("============================================================\n\n")

cat(
  "Predefined panel size: ",
  length(checkpoint_panel),
  "\n",
  sep = ""
)

cat(
  "Available with unique Ensembl ID: ",
  n_available_unique,
  "\n",
  sep = ""
)

cat(
  "Available with multiple Ensembl IDs: ",
  n_available_multiple,
  "\n",
  sep = ""
)

cat(
  "Not available in expression matrix: ",
  n_not_available,
  "\n",
  sep = ""
)

cat(
  "\nNo expression matrix was modified in this step.\n"
)



# ============================================================
# Paper 1 — WORKFLOW-08
# STEP 8C-2
# Build checkpoint expression matrix
# ============================================================

cat("\n============================================================\n")
cat("STEP 8C-2: Build checkpoint expression matrix\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Safety checks
# ------------------------------------------------------------

stopifnot(exists("expr_survival_503"))
stopifnot(exists("checkpoint_mapping_summary"))
stopifnot(exists("checkpoint_annotation"))
stopifnot(exists("clinical_survival_503"))

# ------------------------------------------------------------
# 1. Verify all checkpoint genes are uniquely available
# ------------------------------------------------------------

checkpoint_status_ok <- all(
  checkpoint_mapping_summary$Final_Status ==
    "AVAILABLE_UNIQUE"
)

cat(
  "All checkpoint genes uniquely available: ",
  checkpoint_status_ok,
  "\n",
  sep = ""
)

stopifnot(checkpoint_status_ok)

# ------------------------------------------------------------
# 2. Expression matrix sample IDs
# ------------------------------------------------------------

expr_sample_ids <- colnames(
  expr_survival_503
)[-1]

clinical_sample_ids <- clinical_survival_503$Sample

cat(
  "Expression sample number: ",
  length(expr_sample_ids),
  "\n",
  sep = ""
)

cat(
  "Clinical survival sample number: ",
  length(clinical_sample_ids),
  "\n",
  sep = ""
)

stopifnot(
  identical(
    expr_sample_ids,
    clinical_sample_ids
  )
)

# ------------------------------------------------------------
# 3. Retrieve Ensembl IDs in panel order
# ------------------------------------------------------------

panel_ensembl_ids <- vapply(
  checkpoint_panel,
  function(gene) {
    
    ids <- checkpoint_annotation$ENSEMBL[
      checkpoint_annotation$SYMBOL == gene
    ]
    
    ids_in_expression <- ids[
      ids %in% gene_ids_raw
    ]
    
    stopifnot(
      length(ids_in_expression) == 1
    )
    
    ids_in_expression
    
  },
  character(1)
)

# ------------------------------------------------------------
# 4. Locate Ensembl rows in expression matrix
# ------------------------------------------------------------

gene_row_position <- match(
  panel_ensembl_ids,
  gene_ids_raw
)

stopifnot(
  all(!is.na(gene_row_position))
)

# ------------------------------------------------------------
# 5. Extract expression values
# ------------------------------------------------------------

checkpoint_expr_503 <- expr_survival_503[
  gene_row_position,
  -1,
  drop = FALSE
]

# ------------------------------------------------------------
# 6. Assign gene symbols as row names
# ------------------------------------------------------------

rownames(
  checkpoint_expr_503
) <- checkpoint_panel

# ------------------------------------------------------------
# 7. Convert to numeric matrix explicitly
# ------------------------------------------------------------

checkpoint_expr_503 <- as.matrix(
  checkpoint_expr_503
)

storage.mode(
  checkpoint_expr_503
) <- "numeric"

# ------------------------------------------------------------
# 8. Basic dimensions
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("8. Checkpoint expression matrix dimensions\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Rows: ",
  nrow(checkpoint_expr_503),
  "\n",
  sep = ""
)

cat(
  "Columns: ",
  ncol(checkpoint_expr_503),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 9. Gene order audit
# ------------------------------------------------------------

gene_order_ok <- identical(
  rownames(checkpoint_expr_503),
  checkpoint_panel
)

cat(
  "Gene order matches predefined panel: ",
  gene_order_ok,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 10. Sample order audit
# ------------------------------------------------------------

sample_order_ok <- identical(
  colnames(checkpoint_expr_503),
  clinical_sample_ids
)

cat(
  "Sample order matches clinical survival cohort: ",
  sample_order_ok,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 11. Numerical integrity
# ------------------------------------------------------------

checkpoint_na <- sum(
  is.na(checkpoint_expr_503)
)

checkpoint_nan <- sum(
  is.nan(checkpoint_expr_503)
)

checkpoint_inf <- sum(
  is.infinite(checkpoint_expr_503)
)

checkpoint_negative <- sum(
  checkpoint_expr_503 < 0,
  na.rm = TRUE
)

cat("\n------------------------------------------------------------\n")
cat("11. Numerical integrity\n")
cat("------------------------------------------------------------\n\n")

cat(
  "NA values: ",
  checkpoint_na,
  "\n",
  sep = ""
)

cat(
  "NaN values: ",
  checkpoint_nan,
  "\n",
  sep = ""
)

cat(
  "Inf/-Inf values: ",
  checkpoint_inf,
  "\n",
  sep = ""
)

cat(
  "Negative values: ",
  checkpoint_negative,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 12. Gene-level expression summary
# ------------------------------------------------------------

checkpoint_gene_summary <- data.frame(
  Gene = rownames(checkpoint_expr_503),
  Ensembl_ID = panel_ensembl_ids,
  Min = apply(
    checkpoint_expr_503,
    1,
    min,
    na.rm = TRUE
  ),
  Median = apply(
    checkpoint_expr_503,
    1,
    median,
    na.rm = TRUE
  ),
  Mean = apply(
    checkpoint_expr_503,
    1,
    mean,
    na.rm = TRUE
  ),
  Max = apply(
    checkpoint_expr_503,
    1,
    max,
    na.rm = TRUE
  ),
  SD = apply(
    checkpoint_expr_503,
    1,
    sd,
    na.rm = TRUE
  ),
  stringsAsFactors = FALSE
)

cat("\n------------------------------------------------------------\n")
cat("12. Gene-level expression summary\n")
cat("------------------------------------------------------------\n\n")

print(
  checkpoint_gene_summary,
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Final QC
# ------------------------------------------------------------

final_pass <- (
  nrow(checkpoint_expr_503) == 18 &&
    ncol(checkpoint_expr_503) == 503 &&
    gene_order_ok &&
    sample_order_ok &&
    checkpoint_na == 0 &&
    checkpoint_nan == 0 &&
    checkpoint_inf == 0 &&
    checkpoint_negative == 0
)

cat("\n============================================================\n")
cat("STEP 8C-2 FINAL RESULT\n")
cat("============================================================\n\n")

if (final_pass) {
  
  cat(
    "STEP 8C-2 OVERALL RESULT: PASS\n"
  )
  
  cat(
    "18 checkpoint-associated genes × 503 survival-complete LUAD samples.\n"
  )
  
  cat(
    "Original expression matrix was not modified.\n"
  )
  
} else {
  
  cat(
    "STEP 8C-2 OVERALL RESULT: FAIL\n"
  )
  
}

# ------------------------------------------------------------
# 14. Save only after PASS
# ------------------------------------------------------------

if (final_pass) {
  
  checkpoint_expr_path <- file.path(
    "result_dir/Paper1_tumor-only WGCNA",
    "Workflow-08-STEP8C2_checkpoint_expr_18x503.rds"
  )
  
  saveRDS(
    checkpoint_expr_503,
    checkpoint_expr_path
  )
  
  cat(
    "\nSaved checkpoint expression matrix:\n",
    checkpoint_expr_path,
    "\n",
    sep = ""
  )
  
}

# ============================================================
# Paper 1 — WORKFLOW-08
# STEP 8C-3
# Build core-gene expression matrix
# ============================================================

cat("\n============================================================\n")
cat("STEP 8C-3: Build core-gene expression matrix\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Safety checks
# ------------------------------------------------------------

stopifnot(exists("expr_survival_503"))
stopifnot(exists("clinical_survival_503"))
stopifnot(exists("gene_ids_raw"))

# ------------------------------------------------------------
# 1. Predefined core genes and locked Ensembl IDs
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_ensembl_ids <- c(
  "ENSG00000131471",
  "ENSG00000166927",
  "ENSG00000168685"
)

names(core_ensembl_ids) <- core_genes

# ------------------------------------------------------------
# 2. Verify Ensembl IDs are present in expression matrix
# ------------------------------------------------------------

core_id_available <- core_ensembl_ids %in% gene_ids_raw

cat("Core-gene Ensembl availability:\n")
print(
  data.frame(
    Gene = core_genes,
    Ensembl_ID = unname(core_ensembl_ids),
    Available = core_id_available,
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)

stopifnot(all(core_id_available))

# ------------------------------------------------------------
# 3. Verify expression sample order
# ------------------------------------------------------------

expr_sample_ids <- colnames(expr_survival_503)[-1]
clinical_sample_ids <- clinical_survival_503$Sample

stopifnot(
  identical(
    expr_sample_ids,
    clinical_sample_ids
  )
)

cat(
  "\nExpression sample number: ",
  length(expr_sample_ids),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 4. Locate core-gene rows
# ------------------------------------------------------------

core_row_position <- match(
  unname(core_ensembl_ids),
  gene_ids_raw
)

stopifnot(
  all(!is.na(core_row_position))
)

# ------------------------------------------------------------
# 5. Extract expression matrix
# ------------------------------------------------------------

core_expr_503 <- expr_survival_503[
  core_row_position,
  -1,
  drop = FALSE
]

# ------------------------------------------------------------
# 6. Assign gene symbols
# ------------------------------------------------------------

rownames(core_expr_503) <- core_genes

core_expr_503 <- as.matrix(core_expr_503)

storage.mode(core_expr_503) <- "numeric"

# ------------------------------------------------------------
# 7. Dimension and order QC
# ------------------------------------------------------------

gene_order_ok <- identical(
  rownames(core_expr_503),
  core_genes
)

sample_order_ok <- identical(
  colnames(core_expr_503),
  clinical_sample_ids
)

cat("\n------------------------------------------------------------\n")
cat("Dimension and order QC\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Rows: ",
  nrow(core_expr_503),
  "\n",
  sep = ""
)

cat(
  "Columns: ",
  ncol(core_expr_503),
  "\n",
  sep = ""
)

cat(
  "Gene order correct: ",
  gene_order_ok,
  "\n",
  sep = ""
)

cat(
  "Sample order correct: ",
  sample_order_ok,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 8. Numerical integrity
# ------------------------------------------------------------

core_na <- sum(is.na(core_expr_503))
core_nan <- sum(is.nan(core_expr_503))
core_inf <- sum(is.infinite(core_expr_503))
core_negative <- sum(
  core_expr_503 < 0,
  na.rm = TRUE
)

cat("\n------------------------------------------------------------\n")
cat("Numerical integrity\n")
cat("------------------------------------------------------------\n\n")

cat("NA values: ", core_na, "\n", sep = "")
cat("NaN values: ", core_nan, "\n", sep = "")
cat("Inf/-Inf values: ", core_inf, "\n", sep = "")
cat("Negative values: ", core_negative, "\n", sep = "")

# ------------------------------------------------------------
# 9. Expression summary
# ------------------------------------------------------------

core_gene_summary <- data.frame(
  Gene = rownames(core_expr_503),
  Ensembl_ID = unname(core_ensembl_ids),
  Min = apply(core_expr_503, 1, min),
  Median = apply(core_expr_503, 1, median),
  Mean = apply(core_expr_503, 1, mean),
  Max = apply(core_expr_503, 1, max),
  SD = apply(core_expr_503, 1, sd),
  stringsAsFactors = FALSE
)

cat("\n------------------------------------------------------------\n")
cat("Core-gene expression summary\n")
cat("------------------------------------------------------------\n\n")

print(
  core_gene_summary,
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Final QC
# ------------------------------------------------------------

final_pass <- (
  nrow(core_expr_503) == 3 &&
    ncol(core_expr_503) == 503 &&
    gene_order_ok &&
    sample_order_ok &&
    core_na == 0 &&
    core_nan == 0 &&
    core_inf == 0 &&
    core_negative == 0
)

cat("\n============================================================\n")
cat("STEP 8C-3 FINAL RESULT\n")
cat("============================================================\n\n")

if (final_pass) {
  
  cat("STEP 8C-3 OVERALL RESULT: PASS\n")
  cat("3 core genes × 503 survival-complete LUAD samples.\n")
  cat("Original expression matrix was not modified.\n")
  
} else {
  
  cat("STEP 8C-3 OVERALL RESULT: FAIL\n")
  
}

# ------------------------------------------------------------
# 11. Save only after PASS
# ------------------------------------------------------------

if (final_pass) {
  
  core_expr_path <- file.path(
    "result_dir/Paper1_tumor-only WGCNA",
    "Workflow-08-STEP8C3_core_expr_3x503.rds"
  )
  
  saveRDS(
    core_expr_503,
    core_expr_path
  )
  
  cat(
    "\nSaved core-gene expression matrix:\n",
    core_expr_path,
    "\n",
    sep = ""
  )
}

# ============================================================
# Paper 1 — WORKFLOW-08
# STEP 8C-4
# Spearman association analysis:
# 3 core genes × 18 checkpoint-associated genes
# ============================================================

cat("\n============================================================\n")
cat("STEP 8C-4: Core-gene / checkpoint-gene expression associations\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Safety checks
# ------------------------------------------------------------

stopifnot(exists("core_expr_503"))
stopifnot(exists("checkpoint_expr_503"))
stopifnot(exists("clinical_survival_503"))

# ------------------------------------------------------------
# 1. Define locked analysis panels
# ------------------------------------------------------------

core_genes_locked <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

checkpoint_genes_locked <- c(
  "PDCD1",
  "CD274",
  "PDCD1LG2",
  "CTLA4",
  "CD80",
  "CD86",
  "LAG3",
  "HAVCR2",
  "TIGIT",
  "CD96",
  "VSIR",
  "VTCN1",
  "BTLA",
  "TNFRSF4",
  "TNFRSF9",
  "ICOS",
  "LGALS9",
  "CD200"
)

# ------------------------------------------------------------
# 2. Verify matrix dimensions and order
# ------------------------------------------------------------

stopifnot(
  identical(
    rownames(core_expr_503),
    core_genes_locked
  )
)

stopifnot(
  identical(
    rownames(checkpoint_expr_503),
    checkpoint_genes_locked
  )
)

stopifnot(
  ncol(core_expr_503) == 503
)

stopifnot(
  ncol(checkpoint_expr_503) == 503
)

stopifnot(
  identical(
    colnames(core_expr_503),
    colnames(checkpoint_expr_503)
  )
)

stopifnot(
  identical(
    colnames(core_expr_503),
    clinical_survival_503$Sample
  )
)

cat(
  "Core-gene matrix: ",
  nrow(core_expr_503),
  " × ",
  ncol(core_expr_503),
  "\n",
  sep = ""
)

cat(
  "Checkpoint matrix: ",
  nrow(checkpoint_expr_503),
  " × ",
  ncol(checkpoint_expr_503),
  "\n",
  sep = ""
)

cat(
  "Sample order identical across both matrices: TRUE\n"
)

# ------------------------------------------------------------
# 3. Numerical integrity
# ------------------------------------------------------------

stopifnot(
  sum(is.na(core_expr_503)) == 0,
  sum(is.nan(core_expr_503)) == 0,
  sum(is.infinite(core_expr_503)) == 0,
  sum(is.na(checkpoint_expr_503)) == 0,
  sum(is.nan(checkpoint_expr_503)) == 0,
  sum(is.infinite(checkpoint_expr_503)) == 0
)

# ------------------------------------------------------------
# 4. Predefine 54 analysis pairs
# ------------------------------------------------------------

association_grid <- expand.grid(
  Core_Gene = core_genes_locked,
  Checkpoint_Gene = checkpoint_genes_locked,
  stringsAsFactors = FALSE
)

# Preserve core-gene × checkpoint-gene order
association_grid <- association_grid[
  order(
    match(
      association_grid$Core_Gene,
      core_genes_locked
    ),
    match(
      association_grid$Checkpoint_Gene,
      checkpoint_genes_locked
    )
  ),
  ,
  drop = FALSE
]

rownames(association_grid) <- NULL

cat(
  "\nPredefined association tests: ",
  nrow(association_grid),
  "\n",
  sep = ""
)

stopifnot(
  nrow(association_grid) == 54
)

# ------------------------------------------------------------
# 5. Run Spearman correlations
# ------------------------------------------------------------

association_results <- vector(
  "list",
  nrow(association_grid)
)

for (i in seq_len(nrow(association_grid))) {
  
  core_gene <- association_grid$Core_Gene[i]
  checkpoint_gene <- association_grid$Checkpoint_Gene[i]
  
  x <- as.numeric(
    core_expr_503[core_gene, ]
  )
  
  y <- as.numeric(
    checkpoint_expr_503[checkpoint_gene, ]
  )
  
  complete_pair <- complete.cases(x, y)
  
  n_complete <- sum(complete_pair)
  
  test <- suppressWarnings(
    cor.test(
      x[complete_pair],
      y[complete_pair],
      method = "spearman",
      exact = FALSE
    )
  )
  
  association_results[[i]] <- data.frame(
    Core_Gene = core_gene,
    Checkpoint_Gene = checkpoint_gene,
    N = n_complete,
    Spearman_rho = unname(test$estimate),
    P_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

checkpoint_association_results <- do.call(
  rbind,
  association_results
)

rownames(
  checkpoint_association_results
) <- NULL

# ------------------------------------------------------------
# 6. BH-FDR correction across all 54 tests
# ------------------------------------------------------------

checkpoint_association_results$FDR_BH <- p.adjust(
  checkpoint_association_results$P_value,
  method = "BH"
)

# ------------------------------------------------------------
# 7. Association direction
# ------------------------------------------------------------

checkpoint_association_results$Direction <- ifelse(
  checkpoint_association_results$Spearman_rho > 0,
  "Positive",
  ifelse(
    checkpoint_association_results$Spearman_rho < 0,
    "Negative",
    "Zero"
  )
)

# ------------------------------------------------------------
# 8. Final ordering
# ------------------------------------------------------------

checkpoint_association_results$Core_Gene_Order <- match(
  checkpoint_association_results$Core_Gene,
  core_genes_locked
)

checkpoint_association_results$Checkpoint_Gene_Order <- match(
  checkpoint_association_results$Checkpoint_Gene,
  checkpoint_genes_locked
)

checkpoint_association_results <- checkpoint_association_results[
  order(
    checkpoint_association_results$Core_Gene_Order,
    checkpoint_association_results$Checkpoint_Gene_Order
  ),
  ,
  drop = FALSE
]

checkpoint_association_results$Core_Gene_Order <- NULL
checkpoint_association_results$Checkpoint_Gene_Order <- NULL

rownames(
  checkpoint_association_results
) <- NULL

# ------------------------------------------------------------
# 9. QC summary
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Association analysis QC\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Total tests: ",
  nrow(checkpoint_association_results),
  "\n",
  sep = ""
)

cat(
  "Unique core genes: ",
  length(unique(
    checkpoint_association_results$Core_Gene
  )),
  "\n",
  sep = ""
)

cat(
  "Unique checkpoint genes: ",
  length(unique(
    checkpoint_association_results$Checkpoint_Gene
  )),
  "\n",
  sep = ""
)

cat(
  "N range: ",
  paste(
    range(checkpoint_association_results$N),
    collapse = " - "
  ),
  "\n",
  sep = ""
)

cat(
  "Missing P values: ",
  sum(
    is.na(
      checkpoint_association_results$P_value
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Missing FDR values: ",
  sum(
    is.na(
      checkpoint_association_results$FDR_BH
    )
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 10. Print complete result table
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Complete 54-test association results\n")
cat("------------------------------------------------------------\n\n")

print(
  checkpoint_association_results,
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Significant results by nominal P < 0.05
# ------------------------------------------------------------

nominal_sig <- checkpoint_association_results[
  checkpoint_association_results$P_value < 0.05,
  ,
  drop = FALSE
]

cat("\n------------------------------------------------------------\n")
cat("Nominal associations: P < 0.05\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Number of nominal associations: ",
  nrow(nominal_sig),
  "\n",
  sep = ""
)

if (nrow(nominal_sig) > 0) {
  print(
    nominal_sig,
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 12. FDR-significant results
# ------------------------------------------------------------

fdr_sig <- checkpoint_association_results[
  checkpoint_association_results$FDR_BH < 0.05,
  ,
  drop = FALSE
]

cat("\n------------------------------------------------------------\n")
cat("FDR-significant associations: BH-FDR < 0.05\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Number of FDR-significant associations: ",
  nrow(fdr_sig),
  "\n",
  sep = ""
)

if (nrow(fdr_sig) > 0) {
  print(
    fdr_sig,
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 13. Final QC
# ------------------------------------------------------------

final_pass <- (
  nrow(checkpoint_association_results) == 54 &&
    length(unique(
      checkpoint_association_results$Core_Gene
    )) == 3 &&
    length(unique(
      checkpoint_association_results$Checkpoint_Gene
    )) == 18 &&
    all(
      checkpoint_association_results$N == 503
    ) &&
    all(
      !is.na(
        checkpoint_association_results$P_value
      )
    ) &&
    all(
      !is.na(
        checkpoint_association_results$FDR_BH
      )
    ) &&
    all(
      is.finite(
        checkpoint_association_results$Spearman_rho
      )
    ) &&
    all(
      is.finite(
        checkpoint_association_results$P_value
      )
    ) &&
    all(
      is.finite(
        checkpoint_association_results$FDR_BH
      )
    )
)

cat("\n============================================================\n")
cat("STEP 8C-4 FINAL RESULT\n")
cat("============================================================\n\n")

if (final_pass) {
  
  cat(
    "STEP 8C-4 OVERALL RESULT: PASS\n"
  )
  
  cat(
    "54 predefined Spearman associations completed.\n"
  )
  
  cat(
    "BH-FDR correction applied across all 54 tests.\n"
  )
  
  cat(
    "All associations used N = 503 samples.\n"
  )
  
} else {
  
  cat(
    "STEP 8C-4 OVERALL RESULT: FAIL\n"
  )
}

# ------------------------------------------------------------
# 14. Save after PASS
# ------------------------------------------------------------

if (final_pass) {
  
  checkpoint_association_path <- file.path(
    "result_dir/Paper1_tumor-only WGCNA",
    "Workflow-08-STEP8C4_core_checkpoint_Spearman_54tests.rds"
  )
  
  saveRDS(
    checkpoint_association_results,
    checkpoint_association_path
  )
  
  cat(
    "\nSaved association results:\n",
    checkpoint_association_path,
    "\n",
    sep = ""
  )
}

# ============================================================
# Paper 1
# WORKFLOW-08
# STEP 8C-5A
# ImmuneScore alignment and QC before partial Spearman analysis
# ============================================================
# RECOVERY-A
# Locate previously saved ESTIMATE/WGCNA objects
# ============================================================

cat("\n============================================================\n")
cat("RECOVERY-A: Locate locked ESTIMATE/WGCNA files\n")
cat("============================================================\n")

base_dir <- "input_dir/Paper1_tumor-only WGCNA"

cat("\n--- RData files ---\n")
rdata_files <- list.files(
  base_dir,
  pattern = "\\.(RData|rda|rds)$",
  full.names = TRUE,
  ignore.case = TRUE
)

print(basename(rdata_files))


cat("\n--- Files containing 'estimate' ---\n")
estimate_files <- list.files(
  base_dir,
  pattern = "estimate",
  full.names = TRUE,
  ignore.case = TRUE
)

print(basename(estimate_files))


cat("\n--- Files containing 'workflow-03' ---\n")
workflow03_files <- list.files(
  base_dir,
  pattern = "workflow-03",
  full.names = TRUE,
  ignore.case = TRUE
)

print(basename(workflow03_files))


cat("\n--- Files containing '14' ---\n")
step14_files <- list.files(
  base_dir,
  pattern = "14",
  full.names = TRUE,
  ignore.case = TRUE
)

print(basename(step14_files))


cat("\n============================================================\n")
cat("RECOVERY-A COMPLETE\n")
cat("No files were modified.\n")
cat("============================================================\n")


# ============================================================
# Paper 1
# WORKFLOW-08
# RECOVERY-B
# Inspect locked WORKFLOW-03 final RData
# ============================================================

cat("\n============================================================\n")
cat("RECOVERY-B: Inspect WORKFLOW-03 final RData\n")
cat("============================================================\n")

base_dir <- "input_dir/Paper1_tumor-only WGCNA"

estimate_rdata <- file.path(
  base_dir,
  "Paper1_WORKFLOW03_14F_4D_2E_FINAL.RData"
)

if (!file.exists(estimate_rdata)) {
  stop("ERROR: WORKFLOW-03 final RData file does not exist.")
}

# ------------------------------------------------------------
# 1. Load into isolated temporary environment
# ------------------------------------------------------------

recovery_env <- new.env(parent = emptyenv())

loaded_objects <- load(
  estimate_rdata,
  envir = recovery_env
)

cat("\nObjects stored in the RData file:\n")
print(loaded_objects)


# ------------------------------------------------------------
# 2. Inspect object classes and dimensions
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Object structure summary\n")
cat("============================================================\n")

for (obj_name in loaded_objects) {
  
  obj <- get(obj_name, envir = recovery_env)
  
  cat("\n------------------------------------------------------------\n")
  cat("Object:", obj_name, "\n")
  cat("Class :", paste(class(obj), collapse = ", "), "\n")
  
  if (!is.null(dim(obj))) {
    cat("Dim   :", paste(dim(obj), collapse = " x "), "\n")
  } else {
    cat("Length:", length(obj), "\n")
  }
  
  if (!is.null(names(obj))) {
    cat("Names preview:\n")
    print(head(names(obj), 10))
  }
  
  if (!is.null(rownames(obj))) {
    cat("Row names preview:\n")
    print(head(rownames(obj), 10))
  }
  
  if (!is.null(colnames(obj))) {
    cat("Column names preview:\n")
    print(head(colnames(obj), 10))
  }
}


# ------------------------------------------------------------
# 3. Specifically search for ESTIMATE-related objects
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ESTIMATE-related objects\n")
cat("============================================================\n")

estimate_candidates <- loaded_objects[
  grepl(
    "estimate|immune|stromal|score|tme",
    loaded_objects,
    ignore.case = TRUE
  )
]

print(estimate_candidates)


# ------------------------------------------------------------
# 4. Inspect likely ESTIMATE objects in detail
# ------------------------------------------------------------

if (length(estimate_candidates) > 0) {
  
  for (obj_name in estimate_candidates) {
    
    obj <- get(obj_name, envir = recovery_env)
    
    cat("\n------------------------------------------------------------\n")
    cat("Candidate:", obj_name, "\n")
    cat("Class:", paste(class(obj), collapse = ", "), "\n")
    
    if (!is.null(dim(obj))) {
      cat(
        "Dimension:",
        paste(dim(obj), collapse = " x "),
        "\n"
      )
    }
    
    cat("\nStructure:\n")
    print(str(obj, max.level = 2))
    
    if (!is.null(rownames(obj))) {
      cat("\nRow names:\n")
      print(rownames(obj))
    }
    
    if (!is.null(colnames(obj))) {
      cat("\nFirst 10 column names:\n")
      print(head(colnames(obj), 10))
    }
  }
}


cat("\n============================================================\n")
cat("RECOVERY-B COMPLETE\n")
cat("No objects were loaded into the global environment.\n")
cat("No files were modified.\n")
cat("============================================================\n")


# ============================================================
# Paper 1
# WORKFLOW-08
# RECOVERY-C
# Restore locked estimate_scores_wgcna
# and verify 513 -> 503 sample alignment
# ============================================================

cat("\n============================================================\n")
cat("RECOVERY-C: Restore locked ESTIMATE object\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Confirm recovery environment still exists
# ------------------------------------------------------------

if (!exists("recovery_env")) {
  stop(
    "ERROR: recovery_env is not available. Please rerun RECOVERY-B first."
  )
}

if (!exists(
  "estimate_scores_wgcna",
  envir = recovery_env,
  inherits = FALSE
)) {
  stop(
    "ERROR: estimate_scores_wgcna not found in recovery_env."
  )
}

cat("PASS: Locked estimate_scores_wgcna found in recovery_env.\n")


# ------------------------------------------------------------
# 2. Extract only the locked ESTIMATE object
# ------------------------------------------------------------

estimate_scores_wgcna_locked <- get(
  "estimate_scores_wgcna",
  envir = recovery_env
)


# ------------------------------------------------------------
# 3. Inspect restored object
# ------------------------------------------------------------

cat("\nRestored object:\n")
cat(
  "Class:",
  paste(class(estimate_scores_wgcna_locked), collapse = ", "),
  "\n"
)

cat(
  "Dimension:",
  paste(dim(estimate_scores_wgcna_locked), collapse = " x "),
  "\n"
)

cat("Row names:\n")
print(rownames(estimate_scores_wgcna_locked))

cat("\nFirst 10 sample IDs:\n")
print(head(colnames(estimate_scores_wgcna_locked), 10))


# ------------------------------------------------------------
# 4. Structural checkpoint
# ------------------------------------------------------------

if (
  !identical(
    dim(estimate_scores_wgcna_locked),
    c(3L, 513L)
  )
) {
  stop(
    "ERROR: Locked ESTIMATE object is not exactly 3 x 513."
  )
}

if (
  !identical(
    rownames(estimate_scores_wgcna_locked),
    c(
      "StromalScore",
      "ImmuneScore",
      "ESTIMATEScore"
    )
  )
) {
  stop(
    "ERROR: ESTIMATE row names do not match the locked expected order."
  )
}

cat("PASS: Locked ESTIMATE object has expected 3 x 513 structure.\n")


# ------------------------------------------------------------
# 5. Verify 513 sample IDs
# ------------------------------------------------------------

estimate_ids_513 <- colnames(
  estimate_scores_wgcna_locked
)

if (
  length(estimate_ids_513) != 513L
) {
  stop(
    "ERROR: ESTIMATE sample count is not 513."
  )
}

if (
  anyNA(estimate_ids_513) ||
  any(!nzchar(estimate_ids_513))
) {
  stop(
    "ERROR: ESTIMATE sample IDs contain NA or blank values."
  )
}

if (
  anyDuplicated(estimate_ids_513)
) {
  stop(
    "ERROR: ESTIMATE sample IDs contain duplicates."
  )
}

cat("PASS: 513 ESTIMATE sample IDs are complete and unique.\n")


# ------------------------------------------------------------
# 6. Confirm current 503 survival cohort exists
# ------------------------------------------------------------

if (!exists("clinical_survival_503")) {
  stop(
    "ERROR: clinical_survival_503 is missing."
  )
}

clinical_ids_503 <- as.character(
  clinical_survival_503$Sample
)

if (
  length(clinical_ids_503) != 503L
) {
  stop(
    "ERROR: clinical_survival_503 does not contain exactly 503 samples."
  )
}

if (
  anyNA(clinical_ids_503) ||
  anyDuplicated(clinical_ids_503)
) {
  stop(
    "ERROR: clinical_survival_503 IDs contain NA or duplicates."
  )
}

cat("PASS: 503-sample clinical cohort is complete and unique.\n")


# ------------------------------------------------------------
# 7. Check whether all 503 samples exist in locked ESTIMATE
# ------------------------------------------------------------

estimate_match_position <- match(
  clinical_ids_503,
  estimate_ids_513
)

cat("\n503 -> 513 matching QC:\n")
cat(
  "Matched:",
  sum(!is.na(estimate_match_position)),
  "\n"
)

cat(
  "Unmatched:",
  sum(is.na(estimate_match_position)),
  "\n"
)

if (anyNA(estimate_match_position)) {
  cat("\nUnmatched sample IDs:\n")
  print(
    clinical_ids_503[
      is.na(estimate_match_position)
    ]
  )
  
  stop(
    "ERROR: Some 503 clinical samples are absent from locked ESTIMATE object."
  )
}

cat("PASS: All 503 clinical samples are present in locked ESTIMATE object.\n")


# ------------------------------------------------------------
# 8. Extract exactly the 503 samples in clinical order
# ------------------------------------------------------------

estimate_scores_503 <- estimate_scores_wgcna_locked[
  ,
  estimate_match_position,
  drop = FALSE
]


# ------------------------------------------------------------
# 9. Verify exact sample order
# ------------------------------------------------------------

if (
  !identical(
    colnames(estimate_scores_503),
    clinical_ids_503
  )
) {
  stop(
    "ERROR: Extracted 503-sample ESTIMATE object does not match clinical order."
  )
}

cat("PASS: ESTIMATE 503-sample order exactly matches clinical_survival_503.\n")


# ------------------------------------------------------------
# 10. Verify ImmuneScore
# ------------------------------------------------------------

immune_score_503 <- as.numeric(
  estimate_scores_503[
    "ImmuneScore",
  ]
)

names(immune_score_503) <- colnames(
  estimate_scores_503
)

cat("\nImmuneScore QC:\n")

cat(
  "Length:",
  length(immune_score_503),
  "\n"
)

cat(
  "NA:",
  sum(is.na(immune_score_503)),
  "\n"
)

cat(
  "NaN:",
  sum(is.nan(immune_score_503)),
  "\n"
)

cat(
  "Inf:",
  sum(is.infinite(immune_score_503)),
  "\n"
)

cat("\nImmuneScore summary:\n")
print(summary(immune_score_503))


if (
  anyNA(immune_score_503) ||
  any(!is.finite(immune_score_503))
) {
  stop(
    "ERROR: ImmuneScore contains NA/NaN/Inf values."
  )
}

cat("PASS: ImmuneScore is complete and finite for all 503 samples.\n")


# ------------------------------------------------------------
# 11. Confirm sample order against existing expression objects
# ------------------------------------------------------------

if (!exists("checkpoint_expr_503")) {
  stop(
    "ERROR: checkpoint_expr_503 is missing."
  )
}

if (!exists("core_expr_503")) {
  stop(
    "ERROR: core_expr_503 is missing."
  )
}

if (
  !identical(
    colnames(checkpoint_expr_503),
    clinical_ids_503
  )
) {
  stop(
    "ERROR: checkpoint_expr_503 order differs from clinical cohort."
  )
}

if (
  !identical(
    colnames(core_expr_503),
    clinical_ids_503
  )
) {
  stop(
    "ERROR: core_expr_503 order differs from clinical cohort."
  )
}

if (
  !identical(
    names(immune_score_503),
    clinical_ids_503
  )
) {
  stop(
    "ERROR: ImmuneScore order differs from clinical cohort."
  )
}

cat("PASS: Core expression, checkpoint expression, and ImmuneScore share identical 503-sample order.\n")


# ------------------------------------------------------------
# 12. Restore only the required object to global environment
# ------------------------------------------------------------

estimate_scores_wgcna <- estimate_scores_wgcna_locked


# ------------------------------------------------------------
# 13. Final recovery checkpoint
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("RECOVERY-C COMPLETE\n")
cat("============================================================\n")

cat("Locked ESTIMATE object restored:\n")
cat("  Original:", paste(dim(estimate_scores_wgcna), collapse = " x "), "\n")
cat("  Analysis:", paste(dim(estimate_scores_503), collapse = " x "), "\n")
cat("  ImmuneScore samples:", length(immune_score_503), "\n")

cat("\nSample alignment:\n")
cat("  ESTIMATE 513 -> clinical 503: PASS\n")
cat("  checkpoint_expr_503: PASS\n")
cat("  core_expr_503: PASS\n")
cat("  ImmuneScore completeness: PASS\n")

cat("\nThe locked ESTIMATE baseline was NOT recalculated.\n")
cat("No source files were modified.\n")

cat("\n============================================================\n")
cat("READY FOR STEP 8C-5A RE-RUN\n")
cat("============================================================\n")


# ============================================================
# Paper 1
# WORKFLOW-08
# STEP 8C-5A
# ImmuneScore-adjusted partial Spearman analysis
# 3 core genes × 18 immune-checkpoint genes = 54 tests
# ============================================================

cat("\n============================================================\n")
cat("STEP 8C-5A: ImmuneScore-adjusted partial Spearman analysis\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Required objects
# ------------------------------------------------------------

required_objects <- c(
  "clinical_survival_503",
  "checkpoint_expr_503",
  "core_expr_503",
  "immune_score_503"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1)
  )
]

if (length(missing_objects) > 0) {
  stop(
    paste(
      "ERROR: Required object(s) missing:",
      paste(missing_objects, collapse = ", ")
    )
  )
}

cat("PASS: All required objects are available.\n")


# ------------------------------------------------------------
# 2. Basic structural checks
# ------------------------------------------------------------

if (!identical(dim(core_expr_503), c(3L, 503L))) {
  stop("ERROR: core_expr_503 is not exactly 3 x 503.")
}

if (!identical(dim(checkpoint_expr_503), c(18L, 503L))) {
  stop("ERROR: checkpoint_expr_503 is not exactly 18 x 503.")
}

if (length(immune_score_503) != 503L) {
  stop("ERROR: immune_score_503 does not contain exactly 503 samples.")
}

clinical_ids_503 <- as.character(
  clinical_survival_503$Sample
)

if (!identical(
  colnames(core_expr_503),
  clinical_ids_503
)) {
  stop(
    "ERROR: core_expr_503 sample order differs from clinical cohort."
  )
}

if (!identical(
  colnames(checkpoint_expr_503),
  clinical_ids_503
)) {
  stop(
    "ERROR: checkpoint_expr_503 sample order differs from clinical cohort."
  )
}

if (!identical(
  names(immune_score_503),
  clinical_ids_503
)) {
  stop(
    "ERROR: immune_score_503 sample order differs from clinical cohort."
  )
}

cat("PASS: All three data components share identical 503-sample order.\n")


# ------------------------------------------------------------
# 3. Numerical integrity
# ------------------------------------------------------------

if (
  anyNA(core_expr_503) ||
  any(!is.finite(as.matrix(core_expr_503)))
) {
  stop(
    "ERROR: core_expr_503 contains NA/NaN/Inf."
  )
}

if (
  anyNA(checkpoint_expr_503) ||
  any(!is.finite(as.matrix(checkpoint_expr_503)))
) {
  stop(
    "ERROR: checkpoint_expr_503 contains NA/NaN/Inf."
  )
}

if (
  anyNA(immune_score_503) ||
  any(!is.finite(immune_score_503))
) {
  stop(
    "ERROR: immune_score_503 contains NA/NaN/Inf."
  )
}

cat("PASS: All analysis variables are complete and finite.\n")


# ------------------------------------------------------------
# 4. Define partial Spearman function
#
# Method:
#   1. Rank X, Y, and Z
#   2. Regress rank(X) on rank(Z)
#   3. Regress rank(Y) on rank(Z)
#   4. Correlate the two residual vectors
#
# One covariate => df = N - 3
# ------------------------------------------------------------

partial_spearman_test <- function(
    x,
    y,
    z
) {
  
  complete <- complete.cases(x, y, z)
  
  x <- x[complete]
  y <- y[complete]
  z <- z[complete]
  
  n <- length(x)
  
  if (n < 4L) {
    stop("Insufficient complete observations.")
  }
  
  rx <- rank(
    x,
    ties.method = "average"
  )
  
  ry <- rank(
    y,
    ties.method = "average"
  )
  
  rz <- rank(
    z,
    ties.method = "average"
  )
  
  rx_res <- resid(
    lm(
      rx ~ rz
    )
  )
  
  ry_res <- resid(
    lm(
      ry ~ rz
    )
  )
  
  rho <- cor(
    rx_res,
    ry_res,
    method = "pearson"
  )
  
  df <- n - 3L
  
  if (abs(rho) >= 1) {
    
    p_value <- if (abs(rho) == 1) {
      0
    } else {
      NA_real_
    }
    
  } else {
    
    t_value <- rho * sqrt(
      df / (1 - rho^2)
    )
    
    p_value <- 2 * pt(
      -abs(t_value),
      df = df
    )
  }
  
  data.frame(
    N = n,
    rho = as.numeric(rho),
    P_value = as.numeric(p_value),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 5. Define analysis panel
# ------------------------------------------------------------

core_genes <- rownames(
  core_expr_503
)

checkpoint_genes <- rownames(
  checkpoint_expr_503
)

if (!identical(
  core_genes,
  c("AOC3", "MS4A7", "IL7R")
)) {
  stop(
    "ERROR: Core-gene order is not AOC3 -> MS4A7 -> IL7R."
  )
}

if (length(checkpoint_genes) != 18L) {
  stop(
    "ERROR: Checkpoint panel does not contain exactly 18 genes."
  )
}

cat("\nCore genes:\n")
print(core_genes)

cat("\nCheckpoint genes:\n")
print(checkpoint_genes)


# ------------------------------------------------------------
# 6. Run 54 predefined partial Spearman tests
# ------------------------------------------------------------

partial_results_list <- vector(
  "list",
  length = length(core_genes) *
    length(checkpoint_genes)
)

counter <- 0L

for (core_gene in core_genes) {
  
  x <- as.numeric(
    core_expr_503[
      core_gene,
    ]
  )
  
  for (checkpoint_gene in checkpoint_genes) {
    
    y <- as.numeric(
      checkpoint_expr_503[
        checkpoint_gene,
      ]
    )
    
    result <- partial_spearman_test(
      x = x,
      y = y,
      z = immune_score_503
    )
    
    counter <- counter + 1L
    
    partial_results_list[[counter]] <- data.frame(
      Core_Gene = core_gene,
      Checkpoint_Gene = checkpoint_gene,
      N = result$N,
      Partial_Spearman_rho = result$rho,
      P_value = result$P_value,
      stringsAsFactors = FALSE
    )
  }
}


# ------------------------------------------------------------
# 7. Combine results
# ------------------------------------------------------------

core_checkpoint_partial_ImmuneScore <- do.call(
  rbind,
  partial_results_list
)

rownames(
  core_checkpoint_partial_ImmuneScore
) <- NULL


# ------------------------------------------------------------
# 8. Global BH-FDR correction across all 54 tests
# ------------------------------------------------------------

core_checkpoint_partial_ImmuneScore$FDR <- p.adjust(
  core_checkpoint_partial_ImmuneScore$P_value,
  method = "BH"
)


# ------------------------------------------------------------
# 9. Direction
# ------------------------------------------------------------

core_checkpoint_partial_ImmuneScore$Direction <- ifelse(
  core_checkpoint_partial_ImmuneScore$Partial_Spearman_rho > 0,
  "Positive",
  ifelse(
    core_checkpoint_partial_ImmuneScore$Partial_Spearman_rho < 0,
    "Negative",
    "Zero"
  )
)


# ------------------------------------------------------------
# 10. Adjustment label
# ------------------------------------------------------------

core_checkpoint_partial_ImmuneScore$Adjustment <- (
  "ImmuneScore"
)


# ------------------------------------------------------------
# 11. Structural validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 8C-5A QC\n")
cat("============================================================\n")

cat(
  "Total tests:",
  nrow(core_checkpoint_partial_ImmuneScore),
  "\n"
)

cat(
  "Expected tests: 54\n"
)

if (
  nrow(core_checkpoint_partial_ImmuneScore) != 54L
) {
  stop(
    "ERROR: Number of tests is not 54."
  )
}

if (
  any(
    duplicated(
      core_checkpoint_partial_ImmuneScore[
        ,
        c("Core_Gene", "Checkpoint_Gene")
      ]
    )
  )
) {
  stop(
    "ERROR: Duplicate core/checkpoint test pairs detected."
  )
}

cat(
  "N values:\n"
)

print(
  table(
    core_checkpoint_partial_ImmuneScore$N
  )
)

if (
  !all(
    core_checkpoint_partial_ImmuneScore$N == 503L
  )
) {
  stop(
    "ERROR: Not all partial Spearman tests have N=503."
  )
}

if (
  anyNA(
    core_checkpoint_partial_ImmuneScore$Partial_Spearman_rho
  )
) {
  stop(
    "ERROR: Partial Spearman rho contains NA."
  )
}

if (
  anyNA(
    core_checkpoint_partial_ImmuneScore$P_value
  )
) {
  stop(
    "ERROR: P values contain NA."
  )
}

if (
  anyNA(
    core_checkpoint_partial_ImmuneScore$FDR
  )
) {
  stop(
    "ERROR: FDR values contain NA."
  )
}

cat("PASS: 54/54 tests completed successfully.\n")
cat("PASS: All tests have N=503.\n")
cat("PASS: All rho/P/FDR values are finite.\n")
cat("PASS: BH-FDR applied globally across all 54 tests.\n")


# ------------------------------------------------------------
# 12. Summary counts
# ------------------------------------------------------------

n_nominal <- sum(
  core_checkpoint_partial_ImmuneScore$P_value < 0.05
)

n_fdr <- sum(
  core_checkpoint_partial_ImmuneScore$FDR < 0.05
)

n_positive <- sum(
  core_checkpoint_partial_ImmuneScore$Partial_Spearman_rho > 0
)

n_negative <- sum(
  core_checkpoint_partial_ImmuneScore$Partial_Spearman_rho < 0
)

cat("\nSummary:\n")

cat(
  "Nominal P < 0.05:",
  n_nominal,
  "/ 54\n"
)

cat(
  "BH-FDR < 0.05:",
  n_fdr,
  "/ 54\n"
)

cat(
  "Positive partial rho:",
  n_positive,
  "\n"
)

cat(
  "Negative partial rho:",
  n_negative,
  "\n"
)


# ------------------------------------------------------------
# 13. Display complete result table
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ImmuneScore-adjusted partial Spearman results\n")
cat("============================================================\n")

print(
  core_checkpoint_partial_ImmuneScore,
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Save result
# ------------------------------------------------------------

step8c5_file <- file.path(
  "result_dir/Paper1_tumor-only WGCNA",
  "Workflow-08-STEP8C5_core_checkpoint_partialSpearman_ImmuneScore_54tests.rds"
)

saveRDS(
  core_checkpoint_partial_ImmuneScore,
  step8c5_file
)

cat("\nSaved:\n")
cat(
  normalizePath(
    step8c5_file,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)


# ------------------------------------------------------------
# 15. Final checkpoint
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 8C-5A COMPLETE\n")
cat("============================================================\n")

cat("Analysis: Partial Spearman correlation\n")
cat("Adjustment: ImmuneScore\n")
cat("Core genes: 3\n")
cat("Checkpoint genes: 18\n")
cat("Total tests: 54\n")
cat("Complete-case N per test: 503\n")
cat("FDR: BH across all 54 tests\n")
cat("Source ESTIMATE object: locked Step 14 baseline\n")
cat("Source files: NOT modified\n")
cat("STEP 8C-4 unadjusted results: NOT overwritten\n")

# ============================================================
# Paper 1
# WORKFLOW-08
# STEP 8C-5B
# Paired audit of unadjusted vs ImmuneScore-adjusted
# 3 core genes × 18 immune-checkpoint genes = 54 pairs
#
# Purpose:
#   Compare locked STEP 8C-4 unadjusted Spearman results
#   with locked STEP 8C-5A ImmuneScore-adjusted
#   partial Spearman results.
#
# IMPORTANT:
#   No correlations are recalculated.
#   STEP 8C-4 and STEP 8C-5A source files are NOT overwritten.
# ============================================================


cat("\n============================================================\n")
cat("STEP 8C-5B: Paired audit of unadjusted vs adjusted associations\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Define source files
# ------------------------------------------------------------

step8c4_file <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "Workflow-08-STEP8C4_core_checkpoint_Spearman_54tests.rds"
)

step8c5_file <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "Workflow-08-STEP8C5_core_checkpoint_partialSpearman_ImmuneScore_54tests.rds"
)


# ------------------------------------------------------------
# 2. Check source files
# ------------------------------------------------------------

if (!file.exists(step8c4_file)) {
  stop(
    "ERROR: STEP 8C-4 RDS file not found:\n",
    step8c4_file
  )
}

if (!file.exists(step8c5_file)) {
  stop(
    "ERROR: STEP 8C-5A RDS file not found:\n",
    step8c5_file
  )
}

cat("PASS: Both locked source RDS files exist.\n")


# ------------------------------------------------------------
# 3. Read locked results
# ------------------------------------------------------------

step8c4_results <- readRDS(
  step8c4_file
)

step8c5_results <- readRDS(
  step8c5_file
)


cat("\nSource dimensions:\n")

cat(
  "STEP 8C-4:",
  paste(dim(step8c4_results), collapse = " x "),
  "\n"
)

cat(
  "STEP 8C-5A:",
  paste(dim(step8c5_results), collapse = " x "),
  "\n"
)
# ============================================================
# STEP 8C-5B RECOVERY-A
# Inspect actual STEP 8C-4 / STEP 8C-5A column structures
# No modification / no recalculation
# ============================================================

cat("\n============================================================\n")
cat("STEP 8C-5B RECOVERY-A: Source structure audit\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. STEP 8C-4 actual column names
# ------------------------------------------------------------

cat("\nSTEP 8C-4 actual column names:\n")
print(
  colnames(step8c4_results)
)


# ------------------------------------------------------------
# 2. STEP 8C-5A actual column names
# ------------------------------------------------------------

cat("\nSTEP 8C-5A actual column names:\n")
print(
  colnames(step8c5_results)
)


# ------------------------------------------------------------
# 3. STEP 8C-4 structure
# ------------------------------------------------------------

cat("\nSTEP 8C-4 structure:\n")
str(
  step8c4_results
)


# ------------------------------------------------------------
# 4. STEP 8C-5A structure
# ------------------------------------------------------------

cat("\nSTEP 8C-5A structure:\n")
str(
  step8c5_results
)


# ------------------------------------------------------------
# 5. First 6 rows of STEP 8C-4
# ------------------------------------------------------------

cat("\nFirst 6 rows of STEP 8C-4:\n")
print(
  head(
    step8c4_results
  )
)


# ------------------------------------------------------------
# 6. First 6 rows of STEP 8C-5A
# ------------------------------------------------------------

cat("\nFirst 6 rows of STEP 8C-5A:\n")
print(
  head(
    step8c5_results
  )
)


# ------------------------------------------------------------
# 7. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("RECOVERY-A COMPLETE\n")
cat("============================================================\n")

cat("No source files modified.\n")
cat("No correlation recalculated.\n")
cat("Waiting for actual column structure before rebuilding STEP 8C-5B.\n")


# ============================================================
# Paper 1 - Workflow 08
# STEP 8C-5B
# Paired audit:
# Unadjusted Spearman vs ImmuneScore-adjusted partial Spearman
#
# Purpose:
#   Compare the LOCKED results from:
#   STEP 8C-4  : unadjusted Spearman
#   STEP 8C-5A : ImmuneScore-adjusted partial Spearman
#
# IMPORTANT:
#   - No correlation recalculation
#   - No source RDS modification
#   - Only pair, compare, classify, and summarize
# ============================================================

cat("\n============================================================\n")
cat("STEP 8C-5B: Paired audit of 8C-4 vs 8C-5A\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 0. Output path
# ------------------------------------------------------------

output_file <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "Workflow-08-STEP8C5B_core_checkpoint_pairedAudit_unadjusted_vs_ImmuneScore.rds"
)


# ------------------------------------------------------------
# 1. Source files
# ------------------------------------------------------------

file_8c4 <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "Workflow-08-STEP8C4_core_checkpoint_Spearman_54tests.rds"
)

file_8c5a <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "Workflow-08-STEP8C5_core_checkpoint_partialSpearman_ImmuneScore_54tests.rds"
)


# ------------------------------------------------------------
# 2. Read locked source results
# ------------------------------------------------------------

step8c4_results <- readRDS(file_8c4)
step8c5_results <- readRDS(file_8c5a)

cat("\nSource dimensions:\n")
cat("STEP 8C-4 :", nrow(step8c4_results), "x",
    ncol(step8c4_results), "\n")
cat("STEP 8C-5A:", nrow(step8c5_results), "x",
    ncol(step8c5_results), "\n")


# ------------------------------------------------------------
# 3. Verify exact source structure
# ------------------------------------------------------------

required_c4 <- c(
  "Core_Gene",
  "Checkpoint_Gene",
  "N",
  "Spearman_rho",
  "P_value",
  "FDR_BH",
  "Direction"
)

required_c5a <- c(
  "Core_Gene",
  "Checkpoint_Gene",
  "N",
  "Partial_Spearman_rho",
  "P_value",
  "FDR",
  "Direction",
  "Adjustment"
)

if (!all(required_c4 %in% colnames(step8c4_results))) {
  stop(
    "ERROR: STEP 8C-4 structure does not match expected locked structure."
  )
}

if (!all(required_c5a %in% colnames(step8c5_results))) {
  stop(
    "ERROR: STEP 8C-5A structure does not match expected locked structure."
  )
}

cat("\nPASS: Source column structures verified.\n")


# ------------------------------------------------------------
# 4. Verify number of tests
# ------------------------------------------------------------

if (nrow(step8c4_results) != 54) {
  stop("ERROR: STEP 8C-4 does not contain exactly 54 tests.")
}

if (nrow(step8c5_results) != 54) {
  stop("ERROR: STEP 8C-5A does not contain exactly 54 tests.")
}

cat("PASS: Both source tables contain exactly 54 tests.\n")


# ------------------------------------------------------------
# 5. Verify N
# ------------------------------------------------------------

if (!all(step8c4_results$N == 503)) {
  stop("ERROR: STEP 8C-4 contains N != 503.")
}

if (!all(step8c5_results$N == 503)) {
  stop("ERROR: STEP 8C-5A contains N != 503.")
}

cat("PASS: All tests use N = 503.\n")


# ------------------------------------------------------------
# 6. Verify core genes and checkpoint genes
# ------------------------------------------------------------

expected_core <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

expected_checkpoint <- c(
  "PDCD1",
  "CD274",
  "PDCD1LG2",
  "CTLA4",
  "CD80",
  "CD86",
  "LAG3",
  "HAVCR2",
  "TIGIT",
  "CD96",
  "VSIR",
  "VTCN1",
  "BTLA",
  "TNFRSF4",
  "TNFRSF9",
  "ICOS",
  "LGALS9",
  "CD200"
)

if (!setequal(unique(step8c4_results$Core_Gene), expected_core)) {
  stop("ERROR: STEP 8C-4 core-gene set mismatch.")
}

if (!setequal(unique(step8c5_results$Core_Gene), expected_core)) {
  stop("ERROR: STEP 8C-5A core-gene set mismatch.")
}

if (!setequal(unique(step8c4_results$Checkpoint_Gene),
              expected_checkpoint)) {
  stop("ERROR: STEP 8C-4 checkpoint-gene set mismatch.")
}

if (!setequal(unique(step8c5_results$Checkpoint_Gene),
              expected_checkpoint)) {
  stop("ERROR: STEP 8C-5A checkpoint-gene set mismatch.")
}

cat("PASS: Core-gene and checkpoint-gene sets verified.\n")


# ------------------------------------------------------------
# 7. Create unique paired-test key
# ------------------------------------------------------------

step8c4_results$Pair_Key <- paste(
  step8c4_results$Core_Gene,
  step8c4_results$Checkpoint_Gene,
  sep = "__"
)

step8c5_results$Pair_Key <- paste(
  step8c5_results$Core_Gene,
  step8c5_results$Checkpoint_Gene,
  sep = "__"
)


# ------------------------------------------------------------
# 8. Verify uniqueness
# ------------------------------------------------------------

if (anyDuplicated(step8c4_results$Pair_Key)) {
  stop("ERROR: Duplicate Pair_Key detected in STEP 8C-4.")
}

if (anyDuplicated(step8c5_results$Pair_Key)) {
  stop("ERROR: Duplicate Pair_Key detected in STEP 8C-5A.")
}

if (!setequal(
  step8c4_results$Pair_Key,
  step8c5_results$Pair_Key
)) {
  stop("ERROR: Pair_Key sets do not match between 8C-4 and 8C-5A.")
}

cat("PASS: All 54 pair keys match uniquely.\n")


# ------------------------------------------------------------
# 9. Align STEP 8C-5A to STEP 8C-4
# ------------------------------------------------------------

idx <- match(
  step8c4_results$Pair_Key,
  step8c5_results$Pair_Key
)

if (anyNA(idx)) {
  stop("ERROR: Some STEP 8C-4 tests could not be matched to STEP 8C-5A.")
}

step8c5_aligned <- step8c5_results[idx, ]


# ------------------------------------------------------------
# 10. Verify exact pair alignment
# ------------------------------------------------------------

if (!all(
  step8c4_results$Core_Gene ==
  step8c5_aligned$Core_Gene
)) {
  stop("ERROR: Core-gene alignment failed.")
}

if (!all(
  step8c4_results$Checkpoint_Gene ==
  step8c5_aligned$Checkpoint_Gene
)) {
  stop("ERROR: Checkpoint-gene alignment failed.")
}

if (!all(
  step8c4_results$N ==
  step8c5_aligned$N
)) {
  stop("ERROR: N alignment failed.")
}

cat("PASS: Exact pair alignment verified.\n")


# ------------------------------------------------------------
# 11. Construct paired audit table
# ------------------------------------------------------------

paired_audit <- data.frame(
  
  Core_Gene =
    step8c4_results$Core_Gene,
  
  Checkpoint_Gene =
    step8c4_results$Checkpoint_Gene,
  
  N =
    step8c4_results$N,
  
  Unadjusted_Spearman_rho =
    step8c4_results$Spearman_rho,
  
  Unadjusted_P =
    step8c4_results$P_value,
  
  Unadjusted_FDR_BH =
    step8c4_results$FDR_BH,
  
  Unadjusted_Direction =
    step8c4_results$Direction,
  
  Adjusted_Partial_Spearman_rho =
    step8c5_aligned$Partial_Spearman_rho,
  
  Adjusted_P =
    step8c5_aligned$P_value,
  
  Adjusted_FDR =
    step8c5_aligned$FDR,
  
  Adjusted_Direction =
    step8c5_aligned$Direction,
  
  Adjustment =
    step8c5_aligned$Adjustment,
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 12. Calculate descriptive changes
# ------------------------------------------------------------

paired_audit$Delta_rho <-
  paired_audit$Adjusted_Partial_Spearman_rho -
  paired_audit$Unadjusted_Spearman_rho

paired_audit$Absolute_Delta_rho <-
  abs(paired_audit$Delta_rho)

paired_audit$Absolute_rho_change <-
  abs(paired_audit$Adjusted_Partial_Spearman_rho) -
  abs(paired_audit$Unadjusted_Spearman_rho)


# ------------------------------------------------------------
# 13. Direction comparison
# ------------------------------------------------------------

paired_audit$Direction_Transition <- paste(
  paired_audit$Unadjusted_Direction,
  "->",
  paired_audit$Adjusted_Direction
)

paired_audit$Direction_Switch <-
  paired_audit$Unadjusted_Direction !=
  paired_audit$Adjusted_Direction


# ------------------------------------------------------------
# 14. FDR significance comparison
# ------------------------------------------------------------

paired_audit$Unadjusted_FDR_Significant <-
  paired_audit$Unadjusted_FDR_BH < 0.05

paired_audit$Adjusted_FDR_Significant <-
  paired_audit$Adjusted_FDR < 0.05


paired_audit$FDR_Transition <- paste(
  ifelse(
    paired_audit$Unadjusted_FDR_Significant,
    "FDR<0.05",
    "FDR>=0.05"
  ),
  "->",
  ifelse(
    paired_audit$Adjusted_FDR_Significant,
    "FDR<0.05",
    "FDR>=0.05"
  )
)


# ------------------------------------------------------------
# 15. Combined descriptive classification
# ------------------------------------------------------------

paired_audit$Audit_Class <- "Other"

paired_audit$Audit_Class[
  paired_audit$Unadjusted_FDR_Significant &
    paired_audit$Adjusted_FDR_Significant &
    !paired_audit$Direction_Switch
] <- "Significant_both_same_direction"

paired_audit$Audit_Class[
  paired_audit$Unadjusted_FDR_Significant &
    paired_audit$Adjusted_FDR_Significant &
    paired_audit$Direction_Switch
] <- "Significant_both_direction_switch"

paired_audit$Audit_Class[
  paired_audit$Unadjusted_FDR_Significant &
    !paired_audit$Adjusted_FDR_Significant
] <- "Significance_lost_after_adjustment"

paired_audit$Audit_Class[
  !paired_audit$Unadjusted_FDR_Significant &
    paired_audit$Adjusted_FDR_Significant
] <- "Significance_gained_after_adjustment"


# ------------------------------------------------------------
# 16. QC: finite numerical values
# ------------------------------------------------------------

numeric_columns <- c(
  "Unadjusted_Spearman_rho",
  "Unadjusted_P",
  "Unadjusted_FDR_BH",
  "Adjusted_Partial_Spearman_rho",
  "Adjusted_P",
  "Adjusted_FDR",
  "Delta_rho",
  "Absolute_Delta_rho",
  "Absolute_rho_change"
)

for (v in numeric_columns) {
  
  if (any(!is.finite(paired_audit[[v]]))) {
    
    stop(
      paste(
        "ERROR: Non-finite values detected in",
        v
      )
    )
  }
}

cat("PASS: All numerical audit values are finite.\n")


# ------------------------------------------------------------
# 17. Overall summary
# ------------------------------------------------------------

overall_summary <- data.frame(
  
  Metric = c(
    "Total paired tests",
    "Unadjusted FDR < 0.05",
    "Adjusted FDR < 0.05",
    "Unadjusted positive",
    "Unadjusted negative",
    "Adjusted positive",
    "Adjusted negative",
    "Direction switches",
    "Significance lost",
    "Significance gained"
  ),
  
  N = c(
    nrow(paired_audit),
    
    sum(paired_audit$Unadjusted_FDR_Significant),
    
    sum(paired_audit$Adjusted_FDR_Significant),
    
    sum(paired_audit$Unadjusted_Direction == "Positive"),
    
    sum(paired_audit$Unadjusted_Direction == "Negative"),
    
    sum(paired_audit$Adjusted_Direction == "Positive"),
    
    sum(paired_audit$Adjusted_Direction == "Negative"),
    
    sum(paired_audit$Direction_Switch),
    
    sum(
      paired_audit$Unadjusted_FDR_Significant &
        !paired_audit$Adjusted_FDR_Significant
    ),
    
    sum(
      !paired_audit$Unadjusted_FDR_Significant &
        paired_audit$Adjusted_FDR_Significant
    )
    
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 18. Direction-transition summary
# ------------------------------------------------------------

direction_transition_summary <- as.data.frame(
  table(
    paired_audit$Direction_Transition
  ),
  stringsAsFactors = FALSE
)

colnames(direction_transition_summary) <- c(
  "Direction_Transition",
  "N"
)


# ------------------------------------------------------------
# 19. Per-core-gene summary
# ------------------------------------------------------------

core_gene_summary <- do.call(
  rbind,
  lapply(
    expected_core,
    function(g) {
      
      x <- paired_audit[
        paired_audit$Core_Gene == g,
      ]
      
      data.frame(
        
        Core_Gene = g,
        
        Total_tests = nrow(x),
        
        Unadjusted_FDR_sig =
          sum(x$Unadjusted_FDR_Significant),
        
        Adjusted_FDR_sig =
          sum(x$Adjusted_FDR_Significant),
        
        Unadjusted_positive =
          sum(x$Unadjusted_Direction == "Positive"),
        
        Unadjusted_negative =
          sum(x$Unadjusted_Direction == "Negative"),
        
        Adjusted_positive =
          sum(x$Adjusted_Direction == "Positive"),
        
        Adjusted_negative =
          sum(x$Adjusted_Direction == "Negative"),
        
        Direction_switches =
          sum(x$Direction_Switch),
        
        Significance_lost =
          sum(
            x$Unadjusted_FDR_Significant &
              !x$Adjusted_FDR_Significant
          ),
        
        Significance_gained =
          sum(
            !x$Unadjusted_FDR_Significant &
              x$Adjusted_FDR_Significant
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)


# ------------------------------------------------------------
# 20. Direction-switch table
# ------------------------------------------------------------

direction_switches <- paired_audit[
  paired_audit$Direction_Switch,
]

direction_switches <- direction_switches[
  order(
    direction_switches$Core_Gene,
    direction_switches$Checkpoint_Gene
  ),
]


# ------------------------------------------------------------
# 21. Significance-lost table
# ------------------------------------------------------------

significance_lost <- paired_audit[
  paired_audit$Unadjusted_FDR_Significant &
    !paired_audit$Adjusted_FDR_Significant,
]

significance_lost <- significance_lost[
  order(
    significance_lost$Adjusted_FDR,
    significance_lost$Core_Gene
  ),
]


# ------------------------------------------------------------
# 22. Significant in both analyses
# ------------------------------------------------------------

significant_both <- paired_audit[
  paired_audit$Unadjusted_FDR_Significant &
    paired_audit$Adjusted_FDR_Significant,
]

significant_both <- significant_both[
  order(
    significant_both$Adjusted_FDR,
    significant_both$Core_Gene
  ),
]


# ------------------------------------------------------------
# 23. Largest absolute changes in rho
# ------------------------------------------------------------

largest_absolute_changes <- paired_audit[
  order(
    -paired_audit$Absolute_Delta_rho
  ),
]

largest_absolute_changes <- head(
  largest_absolute_changes,
  15
)


# ------------------------------------------------------------
# 24. Key QC counts
# ------------------------------------------------------------

qc_summary <- list(
  
  source_8c4 =
    normalizePath(file_8c4, winslash = "/", mustWork = FALSE),
  
  source_8c5a =
    normalizePath(file_8c5a, winslash = "/", mustWork = FALSE),
  
  output_file =
    normalizePath(output_file, winslash = "/", mustWork = FALSE),
  
  n_8c4 =
    nrow(step8c4_results),
  
  n_8c5a =
    nrow(step8c5_results),
  
  n_paired =
    nrow(paired_audit),
  
  all_N_503 =
    all(paired_audit$N == 503),
  
  all_pairs_unique =
    !anyDuplicated(
      paste(
        paired_audit$Core_Gene,
        paired_audit$Checkpoint_Gene,
        sep = "__"
      )
    ),
  
  all_numeric_finite =
    all(
      vapply(
        paired_audit[numeric_columns],
        function(x) all(is.finite(x)),
        logical(1)
      )
    ),
  
  source_results_unchanged =
    TRUE
)


# ------------------------------------------------------------
# 25. Assemble final audit object
# ------------------------------------------------------------

step8c5b_audit <- list(
  
  metadata = list(
    
    workflow = "Paper 1 - Workflow 08",
    
    step = "STEP 8C-5B",
    
    purpose =
      "Paired audit of unadjusted Spearman versus ImmuneScore-adjusted partial Spearman",
    
    sample_size = 503,
    
    total_tests = 54,
    
    adjustment = "ImmuneScore",
    
    significance_threshold = 0.05,
    
    multiple_testing =
      "BH FDR preserved from locked source analyses",
    
    recalculation_performed = FALSE,
    
    source_files_modified = FALSE
  ),
  
  paired_audit =
    paired_audit,
  
  overall_summary =
    overall_summary,
  
  direction_transition_summary =
    direction_transition_summary,
  
  core_gene_summary =
    core_gene_summary,
  
  direction_switches =
    direction_switches,
  
  significance_lost =
    significance_lost,
  
  significant_both =
    significant_both,
  
  largest_absolute_changes =
    largest_absolute_changes,
  
  qc_summary =
    qc_summary
)


# ------------------------------------------------------------
# 26. Save audit object
# ------------------------------------------------------------

saveRDS(
  step8c5b_audit,
  output_file
)


# ------------------------------------------------------------
# 27. Final verification
# ------------------------------------------------------------

if (!file.exists(output_file)) {
  stop("ERROR: STEP 8C-5B output file was not created.")
}

test_reload <- readRDS(output_file)

if (!is.list(test_reload)) {
  stop("ERROR: Saved STEP 8C-5B object could not be reloaded correctly.")
}

if (nrow(test_reload$paired_audit) != 54) {
  stop("ERROR: Reloaded paired audit does not contain 54 tests.")
}


# ------------------------------------------------------------
# 28. Print key results
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 8C-5B COMPLETE\n")
cat("============================================================\n")

cat("\nOverall paired-audit summary:\n")
print(overall_summary)

cat("\nDirection transitions:\n")
print(direction_transition_summary)

cat("\nPer-core-gene summary:\n")
print(core_gene_summary)

cat("\nDirection switches:\n")
print(direction_switches)

cat("\nSignificance lost after ImmuneScore adjustment:\n")
print(significance_lost)

cat("\nLargest absolute rho changes (Top 15):\n")
print(largest_absolute_changes)

cat("\nOutput file:\n")
cat(output_file, "\n")

cat("\nQC:\n")
print(qc_summary)

cat("\n============================================================\n")
cat("No correlation recalculated.\n")
cat("STEP 8C-4 source unchanged.\n")
cat("STEP 8C-5A source unchanged.\n")
cat("============================================================\n")


# ============================================================
# Paper 1
# WORKFLOW-08
# STEP 8C-6
# Manuscript-level checkpoint association extraction
# and stratification
#
# Purpose:
#   1. Extract the locked 54-pair results
#   2. Summarize the effect of ImmuneScore adjustment
#   3. Stratify associations descriptively
#   4. Prepare manuscript-level tables
#
# IMPORTANT:
#   - No new statistical calculation
#   - No new Spearman correlation
#   - No new partial Spearman correlation
#   - No new checkpoint selection
#   - No source RDS modification
# ============================================================


cat("\n============================================================\n")
cat("STEP 8C-6: Manuscript-level checkpoint extraction\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

workdir <- "input_dir/Paper1_tumor-only WGCNA"

setwd(workdir)

cat("\nWorking directory:\n")
cat(getwd(), "\n")


# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

file_8c4 <- file.path(
  workdir,
  "Workflow-08-STEP8C4_core_checkpoint_Spearman_54tests.rds"
)

file_8c5a <- file.path(
  workdir,
  "Workflow-08-STEP8C5_core_checkpoint_partialSpearman_ImmuneScore_54tests.rds"
)

file_8c5b <- file.path(
  workdir,
  "Workflow-08-STEP8C5B_core_checkpoint_pairedAudit_unadjusted_vs_ImmuneScore.rds"
)


# ------------------------------------------------------------
# 2. Check files
# ------------------------------------------------------------

stopifnot(file.exists(file_8c4))
stopifnot(file.exists(file_8c5a))
stopifnot(file.exists(file_8c5b))

cat("\nInput files found: TRUE\n")


# ------------------------------------------------------------
# 3. Read locked results
# ------------------------------------------------------------

res_8c4 <- readRDS(file_8c4)
res_8c5a <- readRDS(file_8c5a)
paired_audit <- readRDS(file_8c5b)


# ------------------------------------------------------------
# 4. Inspect object structures
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("4. Object structures\n")
cat("------------------------------------------------------------\n")

cat("\n8C-4:\n")
print(dim(res_8c4))
print(names(res_8c4))

cat("\n8C-5A:\n")
print(dim(res_8c5a))
print(names(res_8c5a))

cat("\n8C-5B:\n")
print(dim(paired_audit))
print(names(paired_audit))


# ------------------------------------------------------------
# 5. Basic QC
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("5. Basic QC\n")
cat("------------------------------------------------------------\n")

stopifnot(nrow(res_8c4) == 54)
stopifnot(nrow(res_8c5a) == 54)

stopifnot(
  all(res_8c4$N == 503),
  all(res_8c5a$N == 503)
)

stopifnot(
  length(unique(
    paste(res_8c4$Core_Gene, res_8c4$Checkpoint_Gene, sep = "__")
  )) == 54
)

stopifnot(
  length(unique(
    paste(res_8c5a$Core_Gene, res_8c5a$Checkpoint_Gene, sep = "__")
  )) == 54
)

cat("8C-4 tests: ", nrow(res_8c4), "\n")
cat("8C-5A tests: ", nrow(res_8c5a), "\n")
cat("N per test: 503\n")


# ============================================================
# STEP 8C-6 recovery
# Confirm the actual 8C-5B paired-audit structure
# ============================================================

cat("\n============================================================\n")
cat("STEP 8C-6 RECOVERY: inspect 8C-5B paired audit\n")
cat("============================================================\n")

cat("\nClass of paired_audit:\n")
print(class(paired_audit))

cat("\nNames of paired_audit:\n")
print(names(paired_audit))

cat("\nClass of paired_audit$paired_audit:\n")
print(class(paired_audit$paired_audit))

cat("\nDimensions of paired_audit$paired_audit:\n")
print(dim(paired_audit$paired_audit))

cat("\nColumns of paired_audit$paired_audit:\n")
print(names(paired_audit$paired_audit))

cat("\nFirst 6 rows:\n")
print(
  head(
    paired_audit$paired_audit
  )
)


# ============================================================
# Paper 1
# WORKFLOW-08
# STEP 8C-6
# Manuscript-level checkpoint result extraction and stratification
#
# RECOVERY VERSION
#
# IMPORTANT:
#   - Uses locked STEP 8C-5B paired_audit$paired_audit
#   - No new correlation calculation
#   - No new statistical test
#   - No new checkpoint-gene selection
#   - No modification of locked source objects
# ============================================================


cat("\n============================================================\n")
cat("STEP 8C-6: Manuscript-level checkpoint extraction\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Set working directory
# ------------------------------------------------------------

workdir <- "input_dir/Paper1_tumor-only WGCNA"

setwd(workdir)

cat("\nWorking directory:\n")
cat(getwd(), "\n")


# ============================================================
# 2. Use the actual locked 8C-5B paired-audit table
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("2. Load locked paired-audit table\n")
cat("------------------------------------------------------------\n")


complete_54 <- paired_audit$paired_audit


# ------------------------------------------------------------
# 2.1 Confirm structure
# ------------------------------------------------------------

cat("\nDimensions:\n")
print(dim(complete_54))

cat("\nColumns:\n")
print(names(complete_54))


stopifnot(
  nrow(complete_54) == 54,
  ncol(complete_54) == 21
)


# ============================================================
# 3. Basic structural QC
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("3. Basic structural QC\n")
cat("------------------------------------------------------------\n")


core_order <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


checkpoint_order <- c(
  "PDCD1",
  "CD274",
  "PDCD1LG2",
  "CTLA4",
  "CD80",
  "CD86",
  "LAG3",
  "HAVCR2",
  "TIGIT",
  "CD96",
  "VSIR",
  "VTCN1",
  "BTLA",
  "TNFRSF4",
  "TNFRSF9",
  "ICOS",
  "LGALS9",
  "CD200"
)


stopifnot(
  all(complete_54$N == 503)
)

stopifnot(
  setequal(
    unique(complete_54$Core_Gene),
    core_order
  )
)

stopifnot(
  setequal(
    unique(complete_54$Checkpoint_Gene),
    checkpoint_order
  )
)

stopifnot(
  length(
    unique(
      paste(
        complete_54$Core_Gene,
        complete_54$Checkpoint_Gene,
        sep = "__"
      )
    )
  ) == 54
)


cat("\nQC:\n")
cat("Total paired tests:", nrow(complete_54), "\n")
cat("N per test:", unique(complete_54$N), "\n")
cat("Core genes:", paste(core_order, collapse = ", "), "\n")
cat("Checkpoint genes:", length(checkpoint_order), "\n")
cat("Unique pairs: TRUE\n")


# ============================================================
# 4. Define significance states
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("4. Define significance states\n")
cat("------------------------------------------------------------")


complete_54$Sig_Unadjusted <-
  complete_54$Unadjusted_FDR_Significant

complete_54$Sig_Adjusted <-
  complete_54$Adjusted_FDR_Significant


# Explicit verification against FDR values
stopifnot(
  identical(
    complete_54$Sig_Unadjusted,
    complete_54$Unadjusted_FDR_BH < 0.05
  )
)

stopifnot(
  identical(
    complete_54$Sig_Adjusted,
    complete_54$Adjusted_FDR < 0.05
  )
)

cat("\nSignificance definitions verified.\n")


# ============================================================
# 5. Define direction strata
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("5. Define direction strata\n")
cat("------------------------------------------------------------\n")


complete_54$Association_Stratum <- NA_character_


# Stable positive
complete_54$Association_Stratum[
  complete_54$Unadjusted_Direction == "Positive" &
    complete_54$Adjusted_Direction == "Positive"
] <- "Stable positive"


# Stable negative
complete_54$Association_Stratum[
  complete_54$Unadjusted_Direction == "Negative" &
    complete_54$Adjusted_Direction == "Negative"
] <- "Stable negative"


# Positive -> Negative
complete_54$Association_Stratum[
  complete_54$Unadjusted_Direction == "Positive" &
    complete_54$Adjusted_Direction == "Negative"
] <- "Positive to negative"


# Negative -> Positive
complete_54$Association_Stratum[
  complete_54$Unadjusted_Direction == "Negative" &
    complete_54$Adjusted_Direction == "Positive"
] <- "Negative to positive"


# Verify no missing strata
stopifnot(
  !any(is.na(complete_54$Association_Stratum))
)


# Verify against locked direction switch
stopifnot(
  all(
    (complete_54$Association_Stratum %in%
       c("Positive to negative",
         "Negative to positive"))
    ==
      complete_54$Direction_Switch
  )
)


# ============================================================
# 6. Define significance-transition strata
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("6. Define significance-transition strata\n")
cat("------------------------------------------------------------\n")


complete_54$Significance_Stratum <- NA_character_


# Significant before AND after
complete_54$Significance_Stratum[
  complete_54$Sig_Unadjusted &
    complete_54$Sig_Adjusted
] <- "Significant before and after adjustment"


# Lost significance
complete_54$Significance_Stratum[
  complete_54$Sig_Unadjusted &
    !complete_54$Sig_Adjusted
] <- "Significance lost after adjustment"


# Gained significance
complete_54$Significance_Stratum[
  !complete_54$Sig_Unadjusted &
    complete_54$Sig_Adjusted
] <- "Significance gained after adjustment"


# Neither
complete_54$Significance_Stratum[
  !complete_54$Sig_Unadjusted &
    !complete_54$Sig_Adjusted
] <- "Not significant before or after adjustment"


stopifnot(
  !any(is.na(complete_54$Significance_Stratum))
)


# ============================================================
# 7. Combined manuscript classification
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("7. Combined manuscript classification\n")
cat("------------------------------------------------------------\n")


complete_54$Manuscript_Class <- paste(
  complete_54$Association_Stratum,
  "|",
  complete_54$Significance_Stratum
)


# ============================================================
# 8. Core-gene summary
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("8. Core-gene summary\n")
cat("------------------------------------------------------------\n")


core_summary_3 <- do.call(
  rbind,
  lapply(core_order, function(g) {
    
    d <- complete_54[
      complete_54$Core_Gene == g,
    ]
    
    data.frame(
      
      Core_Gene = g,
      
      N_Checkpoint_Tests =
        nrow(d),
      
      Unadjusted_FDR_Significant =
        sum(d$Sig_Unadjusted),
      
      Adjusted_FDR_Significant =
        sum(d$Sig_Adjusted),
      
      Significant_Before_And_After =
        sum(
          d$Sig_Unadjusted &
            d$Sig_Adjusted
        ),
      
      Significance_Lost =
        sum(
          d$Sig_Unadjusted &
            !d$Sig_Adjusted
        ),
      
      Significance_Gained =
        sum(
          !d$Sig_Unadjusted &
            d$Sig_Adjusted
        ),
      
      Unadjusted_Positive =
        sum(
          d$Unadjusted_Direction == "Positive"
        ),
      
      Unadjusted_Negative =
        sum(
          d$Unadjusted_Direction == "Negative"
        ),
      
      Adjusted_Positive =
        sum(
          d$Adjusted_Direction == "Positive"
        ),
      
      Adjusted_Negative =
        sum(
          d$Adjusted_Direction == "Negative"
        ),
      
      Direction_Switches =
        sum(d$Direction_Switch),
      
      stringsAsFactors = FALSE
    )
  })
)

rownames(core_summary_3) <- NULL


cat("\nCore-gene summary:\n")
print(core_summary_3)


# ============================================================
# 9. Overall direction transition summary
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("9. Overall direction transition summary\n")
cat("------------------------------------------------------------\n")


transition_summary <- as.data.frame(
  table(
    Before =
      complete_54$Unadjusted_Direction,
    
    After =
      complete_54$Adjusted_Direction
  )
)

names(transition_summary) <- c(
  "Direction_Before",
  "Direction_After",
  "N"
)

transition_summary <-
  transition_summary[
    transition_summary$N > 0,
  ]

rownames(transition_summary) <- NULL


print(transition_summary)


# ============================================================
# 10. Overall significance transition summary
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("10. Overall significance transition summary\n")
cat("------------------------------------------------------------\n")


significance_summary <- data.frame(
  
  Transition = c(
    "Significant before and after adjustment",
    "Significance lost after adjustment",
    "Significance gained after adjustment",
    "Not significant before or after adjustment"
  ),
  
  N = c(
    
    sum(
      complete_54$Sig_Unadjusted &
        complete_54$Sig_Adjusted
    ),
    
    sum(
      complete_54$Sig_Unadjusted &
        !complete_54$Sig_Adjusted
    ),
    
    sum(
      !complete_54$Sig_Unadjusted &
        complete_54$Sig_Adjusted
    ),
    
    sum(
      !complete_54$Sig_Unadjusted &
        !complete_54$Sig_Adjusted
    )
  ),
  
  stringsAsFactors = FALSE
)


print(significance_summary)


# ============================================================
# 11. Adjusted-significant associations
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("11. Adjusted-significant associations\n")
cat("------------------------------------------------------------\n")


adjusted_significant_38 <- complete_54[
  complete_54$Sig_Adjusted,
]


adjusted_significant_38 <-
  adjusted_significant_38[
    order(
      match(
        adjusted_significant_38$Core_Gene,
        core_order
      ),
      match(
        adjusted_significant_38$Checkpoint_Gene,
        checkpoint_order
      )
    ),
  ]


cat(
  "\nAdjusted FDR-significant associations: ",
  nrow(adjusted_significant_38),
  "\n"
)

print(adjusted_significant_38)


# ============================================================
# 12. Direction-switching associations
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("12. Direction-switching associations\n")
cat("------------------------------------------------------------\n")


direction_switch_32 <- complete_54[
  complete_54$Direction_Switch,
]


direction_switch_32 <-
  direction_switch_32[
    order(
      match(
        direction_switch_32$Core_Gene,
        core_order
      ),
      match(
        direction_switch_32$Checkpoint_Gene,
        checkpoint_order
      )
    ),
  ]


cat(
  "\nDirection switches: ",
  nrow(direction_switch_32),
  "\n"
)

print(direction_switch_32)


# ============================================================
# 13. Significance lost
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("13. Significance lost after adjustment\n")
cat("------------------------------------------------------------\n")


significance_lost <- complete_54[
  complete_54$Sig_Unadjusted &
    !complete_54$Sig_Adjusted,
]


significance_lost <-
  significance_lost[
    order(
      match(
        significance_lost$Core_Gene,
        core_order
      ),
      match(
        significance_lost$Checkpoint_Gene,
        checkpoint_order
      )
    ),
  ]


cat(
  "\nSignificance lost: ",
  nrow(significance_lost),
  "\n"
)

print(significance_lost)


# ============================================================
# 14. Significance gained
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("14. Significance gained after adjustment\n")
cat("------------------------------------------------------------\n")


significance_gained <- complete_54[
  !complete_54$Sig_Unadjusted &
    complete_54$Sig_Adjusted,
]


significance_gained <-
  significance_gained[
    order(
      match(
        significance_gained$Core_Gene,
        core_order
      ),
      match(
        significance_gained$Checkpoint_Gene,
        checkpoint_order
      )
    ),
  ]


cat(
  "\nSignificance gained: ",
  nrow(significance_gained),
  "\n"
)

print(significance_gained)


# ============================================================
# 15. Significant before and after adjustment
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("15. Significant before and after adjustment\n")
cat("------------------------------------------------------------\n")


significant_both_34 <- complete_54[
  complete_54$Sig_Unadjusted &
    complete_54$Sig_Adjusted,
]


significant_both_34 <-
  significant_both_34[
    order(
      match(
        significant_both_34$Core_Gene,
        core_order
      ),
      match(
        significant_both_34$Checkpoint_Gene,
        checkpoint_order
      )
    ),
  ]


cat(
  "\nSignificant before and after: ",
  nrow(significant_both_34),
  "\n"
)

print(significant_both_34)


# ============================================================
# 16. Stable positive associations
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("16. Stable positive associations\n")
cat("------------------------------------------------------------\n")


stable_positive <- complete_54[
  complete_54$Association_Stratum ==
    "Stable positive",
]


stable_positive <-
  stable_positive[
    order(
      match(
        stable_positive$Core_Gene,
        core_order
      ),
      match(
        stable_positive$Checkpoint_Gene,
        checkpoint_order
      )
    ),
  ]


cat(
  "\nStable positive: ",
  nrow(stable_positive),
  "\n"
)

print(stable_positive)


# ============================================================
# 17. Stable negative associations
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("17. Stable negative associations\n")
cat("------------------------------------------------------------\n")


stable_negative <- complete_54[
  complete_54$Association_Stratum ==
    "Stable negative",
]


stable_negative <-
  stable_negative[
    order(
      match(
        stable_negative$Core_Gene,
        core_order
      ),
      match(
        stable_negative$Checkpoint_Gene,
        checkpoint_order
      )
    ),
  ]


cat(
  "\nStable negative: ",
  nrow(stable_negative),
  "\n"
)

print(stable_negative)


# ============================================================
# STEP 8C-6
# 18. Positive -> Negative switches that remain significant
# RECOVERY / ROBUST VERSION
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("18. Positive-to-negative switches remaining significant\n")
cat("------------------------------------------------------------\n")


switch_and_significant <- complete_54[
  complete_54$Direction_Transition == "Positive -> Negative" &
    complete_54$Adjusted_FDR_Significant == TRUE,
]


switch_and_significant <-
  switch_and_significant[
    order(
      match(
        switch_and_significant$Core_Gene,
        core_order
      ),
      match(
        switch_and_significant$Checkpoint_Gene,
        checkpoint_order
      )
    ),
  ]


cat(
  "\nPositive-to-negative switches remaining significant: ",
  nrow(switch_and_significant),
  "\n"
)


print(
  switch_and_significant
)


# ============================================================
# STEP 8C-6
# Internal consistency QC
# ============================================================

cat("\n============================================================\n")
cat("STEP 8C-6 INTERNAL CONSISTENCY QC\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Direction transition consistency
# ------------------------------------------------------------

expected_direction_transition <- paste(
  complete_54$Unadjusted_Direction,
  complete_54$Adjusted_Direction,
  sep = " -> "
)


cat("\nDirection transition consistency:\n")

print(
  table(
    complete_54$Direction_Transition ==
      expected_direction_transition
  )
)


stopifnot(
  all(
    complete_54$Direction_Transition ==
      expected_direction_transition
  )
)


# ------------------------------------------------------------
# 2. Direction switch consistency
# ------------------------------------------------------------

expected_direction_switch <-
  complete_54$Unadjusted_Direction !=
  complete_54$Adjusted_Direction


cat("\nDirection switch consistency:\n")

print(
  table(
    complete_54$Direction_Switch ==
      expected_direction_switch
  )
)


stopifnot(
  all(
    complete_54$Direction_Switch ==
      expected_direction_switch
  )
)


# ------------------------------------------------------------
# 3. Unadjusted significance consistency
# ------------------------------------------------------------

expected_unadjusted_sig <-
  complete_54$Unadjusted_FDR_BH < 0.05


cat("\nUnadjusted significance consistency:\n")

print(
  table(
    complete_54$Unadjusted_FDR_Significant ==
      expected_unadjusted_sig
  )
)


stopifnot(
  all(
    complete_54$Unadjusted_FDR_Significant ==
      expected_unadjusted_sig
  )
)


# ------------------------------------------------------------
# 4. Adjusted significance consistency
# ------------------------------------------------------------

expected_adjusted_sig <-
  complete_54$Adjusted_FDR < 0.05


cat("\nAdjusted significance consistency:\n")

print(
  table(
    complete_54$Adjusted_FDR_Significant ==
      expected_adjusted_sig
  )
)


stopifnot(
  all(
    complete_54$Adjusted_FDR_Significant ==
      expected_adjusted_sig
  )
)


# ------------------------------------------------------------
# 5. Global locked counts
# ------------------------------------------------------------

cat("\nGlobal counts:\n")

cat(
  "Total pairs = ",
  nrow(complete_54),
  "\n"
)

cat(
  "Unadjusted significant = ",
  sum(
    complete_54$Unadjusted_FDR_Significant
  ),
  "\n"
)

cat(
  "Adjusted significant = ",
  sum(
    complete_54$Adjusted_FDR_Significant
  ),
  "\n"
)

cat(
  "Direction switches = ",
  sum(
    complete_54$Direction_Switch
  ),
  "\n"
)

cat(
  "Positive -> Negative = ",
  sum(
    complete_54$Direction_Transition ==
      "Positive -> Negative"
  ),
  "\n"
)

cat(
  "Negative -> Positive = ",
  sum(
    complete_54$Direction_Transition ==
      "Negative -> Positive"
  ),
  "\n"
)

cat(
  "Significance lost = ",
  sum(
    complete_54$Unadjusted_FDR_Significant &
      !complete_54$Adjusted_FDR_Significant
  ),
  "\n"
)

cat(
  "Significance gained = ",
  sum(
    !complete_54$Unadjusted_FDR_Significant &
      complete_54$Adjusted_FDR_Significant
  ),
  "\n"
)


cat("\nQC STATUS: PASS\n")


# ============================================================
# 19. Core × direction-stratum summary
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("19. Core × direction-stratum summary\n")
cat("------------------------------------------------------------\n")


core_stratum_summary <- as.data.frame(
  table(
    Core_Gene =
      complete_54$Core_Gene,
    
    Association_Stratum =
      complete_54$Association_Stratum
  )
)

names(core_stratum_summary) <- c(
  "Core_Gene",
  "Association_Stratum",
  "N"
)


core_stratum_summary$Core_Gene <-
  factor(
    core_stratum_summary$Core_Gene,
    levels = core_order
  )


core_stratum_summary <-
  core_stratum_summary[
    order(
      core_stratum_summary$Core_Gene,
      core_stratum_summary$Association_Stratum
    ),
  ]


core_stratum_summary$Core_Gene <-
  as.character(
    core_stratum_summary$Core_Gene
  )


print(core_stratum_summary)


# ============================================================
# 20. Core × significance-transition summary
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("20. Core × significance-transition summary\n")
cat("------------------------------------------------------------\n")


core_significance_summary <- do.call(
  rbind,
  lapply(core_order, function(g) {
    
    d <- complete_54[
      complete_54$Core_Gene == g,
    ]
    
    data.frame(
      
      Core_Gene = g,
      
      Significant_Before =
        sum(d$Sig_Unadjusted),
      
      Significant_After =
        sum(d$Sig_Adjusted),
      
      Significant_Both =
        sum(
          d$Sig_Unadjusted &
            d$Sig_Adjusted
        ),
      
      Lost =
        sum(
          d$Sig_Unadjusted &
            !d$Sig_Adjusted
        ),
      
      Gained =
        sum(
          !d$Sig_Unadjusted &
            d$Sig_Adjusted
        ),
      
      stringsAsFactors = FALSE
    )
  })
)


rownames(core_significance_summary) <- NULL


print(core_significance_summary)


# ============================================================
# 21. Manuscript-oriented 54-pair table
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("21. Manuscript-oriented 54-pair table\n")
cat("------------------------------------------------------------\n")


manuscript_table_54 <- complete_54[, c(
  
  "Core_Gene",
  "Checkpoint_Gene",
  "N",
  
  "Unadjusted_Spearman_rho",
  "Unadjusted_P",
  "Unadjusted_FDR_BH",
  "Unadjusted_Direction",
  
  "Adjusted_Partial_Spearman_rho",
  "Adjusted_P",
  "Adjusted_FDR",
  "Adjusted_Direction",
  
  "Delta_rho",
  "Absolute_Delta_rho",
  
  "Direction_Transition",
  "Direction_Switch",
  
  "Unadjusted_FDR_Significant",
  "Adjusted_FDR_Significant",
  
  "FDR_Transition",
  "Audit_Class",
  
  "Association_Stratum",
  "Significance_Stratum"
)]


manuscript_table_54 <-
  manuscript_table_54[
    order(
      match(
        manuscript_table_54$Core_Gene,
        core_order
      ),
      match(
        manuscript_table_54$Checkpoint_Gene,
        checkpoint_order
      )
    ),
  ]


# ============================================================
# 22. Overall manuscript strata
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("22. Overall manuscript strata\n")
cat("------------------------------------------------------------\n")


manuscript_strata <- data.frame(
  
  Stratum = c(
    "Stable positive",
    "Stable negative",
    "Positive to negative",
    "Negative to positive"
  ),
  
  N = c(
    
    sum(
      complete_54$Association_Stratum ==
        "Stable positive"
    ),
    
    sum(
      complete_54$Association_Stratum ==
        "Stable negative"
    ),
    
    sum(
      complete_54$Association_Stratum ==
        "Positive to negative"
    ),
    
    sum(
      complete_54$Association_Stratum ==
        "Negative to positive"
    )
  ),
  
  Adjusted_FDR_Significant = c(
    
    sum(
      complete_54$Association_Stratum ==
        "Stable positive" &
        complete_54$Sig_Adjusted
    ),
    
    sum(
      complete_54$Association_Stratum ==
        "Stable negative" &
        complete_54$Sig_Adjusted
    ),
    
    sum(
      complete_54$Association_Stratum ==
        "Positive to negative" &
        complete_54$Sig_Adjusted
    ),
    
    sum(
      complete_54$Association_Stratum ==
        "Negative to positive" &
        complete_54$Sig_Adjusted
    )
  ),
  
  stringsAsFactors = FALSE
)


print(manuscript_strata)


# ============================================================
# 23. Final QC
# ============================================================

cat("\n============================================================\n")
cat("23. FINAL QC\n")
cat("============================================================\n")


qc <- list(
  
  n_pairs =
    nrow(complete_54),
  
  n_core_genes =
    length(unique(complete_54$Core_Gene)),
  
  n_checkpoint_genes =
    length(unique(complete_54$Checkpoint_Gene)),
  
  all_N_503 =
    all(complete_54$N == 503),
  
  all_pairs_unique =
    length(
      unique(
        paste(
          complete_54$Core_Gene,
          complete_54$Checkpoint_Gene,
          sep = "__"
        )
      )
    ) == 54,
  
  unadjusted_significant =
    sum(complete_54$Sig_Unadjusted),
  
  adjusted_significant =
    sum(complete_54$Sig_Adjusted),
  
  significant_both =
    sum(
      complete_54$Sig_Unadjusted &
        complete_54$Sig_Adjusted
    ),
  
  significance_lost =
    sum(
      complete_54$Sig_Unadjusted &
        !complete_54$Sig_Adjusted
    ),
  
  significance_gained =
    sum(
      !complete_54$Sig_Unadjusted &
        complete_54$Sig_Adjusted
    ),
  
  direction_switches =
    sum(complete_54$Direction_Switch),
  
  stable_positive =
    sum(
      complete_54$Association_Stratum ==
        "Stable positive"
    ),
  
  stable_negative =
    sum(
      complete_54$Association_Stratum ==
        "Stable negative"
    ),
  
  positive_to_negative =
    sum(
      complete_54$Association_Stratum ==
        "Positive to negative"
    ),
  
  negative_to_positive =
    sum(
      complete_54$Association_Stratum ==
        "Negative to positive"
    )
)


print(qc)


# ------------------------------------------------------------
# Expected locked values
# ------------------------------------------------------------

stopifnot(
  qc$n_pairs == 54,
  qc$n_core_genes == 3,
  qc$n_checkpoint_genes == 18,
  qc$all_N_503,
  qc$all_pairs_unique,
  
  qc$unadjusted_significant == 48,
  
  qc$adjusted_significant == 38,
  
  qc$significant_both == 34,
  
  qc$significance_lost == 14,
  
  qc$significance_gained == 4,
  
  qc$direction_switches == 32,
  
  qc$stable_positive == 21,
  
  qc$stable_negative == 1,
  
  qc$positive_to_negative == 32,
  
  qc$negative_to_positive == 0
)


cat("\nQC STATUS: PASS\n")


# ============================================================
# 24. Save STEP 8C-6 master object
# ============================================================

output_rds <- file.path(
  workdir,
  "Workflow-08-STEP8C6_manuscript_checkpoint_extraction.rds"
)


step8c6 <- list(
  
  complete_54 =
    complete_54,
  
  manuscript_table_54 =
    manuscript_table_54,
  
  core_summary_3 =
    core_summary_3,
  
  transition_summary =
    transition_summary,
  
  significance_summary =
    significance_summary,
  
  adjusted_significant_38 =
    adjusted_significant_38,
  
  significant_both_34 =
    significant_both_34,
  
  direction_switch_32 =
    direction_switch_32,
  
  significance_lost =
    significance_lost,
  
  significance_gained =
    significance_gained,
  
  stable_positive =
    stable_positive,
  
  stable_negative =
    stable_negative,
  
  switch_and_significant =
    switch_and_significant,
  
  manuscript_strata =
    manuscript_strata,
  
  core_stratum_summary =
    core_stratum_summary,
  
  core_significance_summary =
    core_significance_summary,
  
  qc =
    qc
)


saveRDS(
  step8c6,
  output_rds
)


cat("\nSTEP 8C-6 master RDS saved:\n")
cat(output_rds, "\n")


# ============================================================
# 25. Save manuscript-level TSV files
# ============================================================


write.table(
  manuscript_table_54,
  file = file.path(
    workdir,
    "Workflow-08-STEP8C6_manuscript_table_54.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  core_summary_3,
  file = file.path(
    workdir,
    "Workflow-08-STEP8C6_core_summary_3.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  significance_summary,
  file = file.path(
    workdir,
    "Workflow-08-STEP8C6_significance_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  adjusted_significant_38,
  file = file.path(
    workdir,
    "Workflow-08-STEP8C6_adjusted_significant_38.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  significant_both_34,
  file = file.path(
    workdir,
    "Workflow-08-STEP8C6_significant_both_34.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  direction_switch_32,
  file = file.path(
    workdir,
    "Workflow-08-STEP8C6_direction_switch_32.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  significance_lost,
  file = file.path(
    workdir,
    "Workflow-08-STEP8C6_significance_lost.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  significance_gained,
  file = file.path(
    workdir,
    "Workflow-08-STEP8C6_significance_gained.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat("\n============================================================\n")
cat("STEP 8C-6 COMPLETED\n")
cat("============================================================\n")

cat("\nNo new statistical analysis was performed.\n")
cat("All manuscript-level outputs were extracted from locked 8C-5B results.\n")


# ============================================================
# Paper 1 — Workflow 08
# STEP 8C-7A
# Manuscript result-definition lock
#
# PURPOSE
# ------------------------------------------------------------
# Lock the manuscript-level numerical definitions using ONLY
# the already locked STEP 8C-6 master object.
#
# NO new Spearman correlation is calculated.
# NO new partial Spearman correlation is calculated.
# NO expression matrix is re-read.
# NO ImmuneScore is recalculated.
#
# Locked manuscript-level values:
#
#   54 = total predefined core × checkpoint pairs
#   48 = unadjusted FDR < 0.05
#   38 = adjusted FDR < 0.05
#   34 = significant both before and after adjustment
#   14 = significance lost after adjustment
#    4 = significance gained after adjustment
#   32 = direction switches
#   21 = positive → negative switches remaining
#        significant after adjustment
#
# Core genes:
#   AOC3
#   MS4A7
#   IL7R
#
# Source:
#   Workflow-08-STEP8C6_manuscript_checkpoint_extraction.rds
#
# ============================================================


cat("\n============================================================\n")
cat("STEP 8C-7A: Manuscript result-definition lock\n")
cat("============================================================\n\n")


# ============================================================
# 0. WORKING DIRECTORY
# ============================================================

workdir <- "input_dir/Paper1_tumor-only WGCNA"

setwd(workdir)

cat("Working directory:\n")
cat(normalizePath(getwd(), winslash = "/", mustWork = TRUE), "\n\n")


# ============================================================
# 1. LOAD LOCKED STEP 8C-6 OBJECT
# ============================================================

source_file <- file.path(
  workdir,
  "Workflow-08-STEP8C6_manuscript_checkpoint_extraction.rds"
)

stopifnot(file.exists(source_file))

step8c6 <- readRDS(source_file)

cat("STEP 8C-6 source file loaded:\n")
cat(basename(source_file), "\n\n")


# ============================================================
# 2. VERIFY ACTUAL STEP 8C-6 OBJECT STRUCTURE
# ============================================================

cat("STEP 8C-6 object structure:\n")
print(names(step8c6))
cat("\n")


required_objects <- c(
  "complete_54",
  "core_summary_3",
  "significance_summary"
)

missing_objects <- setdiff(
  required_objects,
  names(step8c6)
)

if (length(missing_objects) > 0) {
  stop(
    "Required STEP 8C-6 objects are missing: ",
    paste(missing_objects, collapse = ", ")
  )
}

cat("Required STEP 8C-6 objects: PASS\n\n")


# ============================================================
# 3. EXTRACT THE LOCKED 54-PAIR TABLE
# ============================================================

complete_54 <- step8c6$complete_54

stopifnot(
  is.data.frame(complete_54),
  nrow(complete_54) == 54
)

cat("complete_54 dimensions:\n")
print(dim(complete_54))
cat("\n")

cat("54-pair master table: PASS\n\n")


# ============================================================
# 4. VERIFY CORE GENES
# ============================================================

expected_core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

observed_core_genes <- unique(
  complete_54$Core_Gene
)

cat("Observed core genes:\n")
print(observed_core_genes)
cat("\n")

stopifnot(
  setequal(
    observed_core_genes,
    expected_core_genes
  )
)

cat("Core gene identity: PASS\n\n")


# ============================================================
# 5. VERIFY 18 CHECKPOINT GENES
# ============================================================

observed_checkpoint_genes <- unique(
  complete_54$Checkpoint_Gene
)

cat("Number of checkpoint genes:\n")
cat(length(observed_checkpoint_genes), "\n\n")

stopifnot(
  length(observed_checkpoint_genes) == 18
)

cat("18-checkpoint panel: PASS\n\n")


# ============================================================
# 6. VERIFY PAIR UNIQUENESS AND N
# ============================================================

pair_id <- paste(
  complete_54$Core_Gene,
  complete_54$Checkpoint_Gene,
  sep = "___"
)

cat("Number of unique core × checkpoint pairs:\n")
cat(length(unique(pair_id)), "\n\n")

stopifnot(
  length(unique(pair_id)) == 54,
  all(complete_54$N == 503)
)

cat("Pair uniqueness and N=503: PASS\n\n")


# ============================================================
# 7. VERIFY REQUIRED COLUMNS
# ============================================================

required_columns <- c(
  "Core_Gene",
  "Checkpoint_Gene",
  "N",
  "Unadjusted_Spearman_rho",
  "Unadjusted_P",
  "Unadjusted_FDR_BH",
  "Unadjusted_Direction",
  "Adjusted_Partial_Spearman_rho",
  "Adjusted_P",
  "Adjusted_FDR",
  "Adjusted_Direction",
  "Direction_Transition",
  "Direction_Switch",
  "Unadjusted_FDR_Significant",
  "Adjusted_FDR_Significant",
  "FDR_Transition"
)

missing_columns <- setdiff(
  required_columns,
  colnames(complete_54)
)

if (length(missing_columns) > 0) {
  stop(
    "Required columns missing from complete_54: ",
    paste(missing_columns, collapse = ", ")
  )
}

cat("Required source columns: PASS\n\n")


# ============================================================
# 8. LOCK THE EIGHT MANUSCRIPT-LEVEL NUMBERS
# ============================================================

# ------------------------------------------------------------
# 8.1 Total predefined associations
# ------------------------------------------------------------

LOCK_54 <- nrow(complete_54)


# ------------------------------------------------------------
# 8.2 Unadjusted FDR < 0.05
# ------------------------------------------------------------

LOCK_48 <- sum(
  complete_54$Unadjusted_FDR_Significant == TRUE,
  na.rm = TRUE
)


# ------------------------------------------------------------
# 8.3 Adjusted FDR < 0.05
# ------------------------------------------------------------

LOCK_38 <- sum(
  complete_54$Adjusted_FDR_Significant == TRUE,
  na.rm = TRUE
)


# ------------------------------------------------------------
# 8.4 Significant both
# ------------------------------------------------------------

LOCK_34 <- sum(
  complete_54$Unadjusted_FDR_Significant == TRUE &
    complete_54$Adjusted_FDR_Significant == TRUE,
  na.rm = TRUE
)


# ------------------------------------------------------------
# 8.5 Significance lost
# ------------------------------------------------------------

LOCK_14 <- sum(
  complete_54$Unadjusted_FDR_Significant == TRUE &
    complete_54$Adjusted_FDR_Significant == FALSE,
  na.rm = TRUE
)


# ------------------------------------------------------------
# 8.6 Significance gained
# ------------------------------------------------------------

LOCK_4 <- sum(
  complete_54$Unadjusted_FDR_Significant == FALSE &
    complete_54$Adjusted_FDR_Significant == TRUE,
  na.rm = TRUE
)


# ------------------------------------------------------------
# 8.7 Direction switches
# ------------------------------------------------------------

LOCK_32 <- sum(
  complete_54$Direction_Switch == TRUE,
  na.rm = TRUE
)


# ------------------------------------------------------------
# 8.8 Positive → negative switches remaining significant
# ------------------------------------------------------------

LOCK_21 <- sum(
  complete_54$Direction_Transition ==
    "Positive -> Negative" &
    complete_54$Adjusted_FDR_Significant == TRUE,
  na.rm = TRUE
)


# ============================================================
# 9. MANUSCRIPT RESULT-DEFINITION TABLE
# ============================================================

manuscript_lock <- data.frame(
  Metric = c(
    "Total predefined core-checkpoint pairs",
    "Unadjusted FDR < 0.05",
    "Adjusted FDR < 0.05",
    "Significant both before and after adjustment",
    "Significance lost after adjustment",
    "Significance gained after adjustment",
    "Direction switches",
    "Positive-to-negative switches remaining significant"
  ),
  
  Locked_Value = c(
    LOCK_54,
    LOCK_48,
    LOCK_38,
    LOCK_34,
    LOCK_14,
    LOCK_4,
    LOCK_32,
    LOCK_21
  ),
  
  Definition = c(
    "All predefined core gene × checkpoint gene pairs",
    
    "Unadjusted Spearman association with BH-FDR < 0.05",
    
    "ImmuneScore-adjusted partial Spearman association with BH-FDR < 0.05",
    
    "FDR significant in both unadjusted and adjusted analyses",
    
    "FDR significant before but not after ImmuneScore adjustment",
    
    "FDR not significant before but significant after ImmuneScore adjustment",
    
    "Unadjusted and adjusted association directions differ",
    
    "Positive-to-negative direction switch with adjusted BH-FDR < 0.05"
  ),
  
  stringsAsFactors = FALSE
)


cat("\n============================================================\n")
cat("LOCKED MANUSCRIPT-LEVEL RESULTS\n")
cat("============================================================\n\n")

print(
  manuscript_lock,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 10. HARD LOCK QC
# ============================================================

stopifnot(
  LOCK_54 == 54,
  LOCK_48 == 48,
  LOCK_38 == 38,
  LOCK_34 == 34,
  LOCK_14 == 14,
  LOCK_4  == 4,
  LOCK_32 == 32,
  LOCK_21 == 21
)

cat("Expected manuscript values:\n")
cat("54 / 48 / 38 / 34 / 14 / 4 / 32 / 21\n")
cat("QC STATUS: PASS\n\n")


# ============================================================
# 11. ARITHMETIC CONSISTENCY
# ============================================================

n_neither <- sum(
  complete_54$Unadjusted_FDR_Significant == FALSE &
    complete_54$Adjusted_FDR_Significant == FALSE,
  na.rm = TRUE
)

stopifnot(
  LOCK_48 == LOCK_34 + LOCK_14,
  LOCK_38 == LOCK_34 + LOCK_4,
  LOCK_54 ==
    LOCK_34 +
    LOCK_14 +
    LOCK_4 +
    n_neither
)

cat("Arithmetic consistency:\n")
cat("  48 = 34 + 14              : PASS\n")
cat("  38 = 34 + 4               : PASS\n")
cat(
  "  54 = 34 + 14 + 4 + ",
  n_neither,
  " : PASS\n",
  sep = ""
)

cat("\n")


# ============================================================
# 12. USE THE ACTUAL LOCKED DIRECTION LABELS
#
# IMPORTANT:
# STEP 8C-6 actually stores:
#
#   "Negative -> Negative"
#   "Positive -> Negative"
#   "Positive -> Positive"
#
# Therefore we DO NOT use the non-existent labels:
#   "Stable negative"
#   "Stable positive"
#
# ============================================================

cat("Actual Direction_Transition values:\n")
print(
  table(
    complete_54$Direction_Transition,
    useNA = "ifany"
  )
)

cat("\n")


# ------------------------------------------------------------
# Define stable/switch strata using the ACTUAL labels
# ------------------------------------------------------------

negative_to_negative <- sum(
  complete_54$Direction_Transition ==
    "Negative -> Negative",
  na.rm = TRUE
)

positive_to_positive <- sum(
  complete_54$Direction_Transition ==
    "Positive -> Positive",
  na.rm = TRUE
)

positive_to_negative <- sum(
  complete_54$Direction_Transition ==
    "Positive -> Negative",
  na.rm = TRUE
)

negative_to_positive <- sum(
  complete_54$Direction_Transition ==
    "Negative -> Positive",
  na.rm = TRUE
)


# ------------------------------------------------------------
# Direction QC
# ------------------------------------------------------------

stopifnot(
  negative_to_negative == 1,
  positive_to_positive == 21,
  positive_to_negative == 32,
  negative_to_positive == 0,
  
  LOCK_32 == positive_to_negative,
  
  positive_to_negative +
    positive_to_positive +
    negative_to_negative +
    negative_to_positive ==
    LOCK_54
)

cat("Direction transition summary:\n")
cat(
  "  Negative -> Negative = ",
  negative_to_negative,
  "\n",
  sep = ""
)

cat(
  "  Positive -> Positive = ",
  positive_to_positive,
  "\n",
  sep = ""
)

cat(
  "  Positive -> Negative = ",
  positive_to_negative,
  "\n",
  sep = ""
)

cat(
  "  Negative -> Positive = ",
  negative_to_positive,
  "\n\n",
  sep = ""
)

cat("Direction transition QC: PASS\n\n")


# ============================================================
# 13. THREE CORE-GENE SUMMARY
# ============================================================

core_order <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


core_summary_locked <- do.call(
  rbind,
  lapply(
    core_order,
    function(g) {
      
      x <- complete_54[
        complete_54$Core_Gene == g,
        ,
        drop = FALSE
      ]
      
      data.frame(
        
        Core_Gene = g,
        
        N_pairs = nrow(x),
        
        Unadjusted_FDR_significant =
          sum(
            x$Unadjusted_FDR_Significant == TRUE,
            na.rm = TRUE
          ),
        
        Adjusted_FDR_significant =
          sum(
            x$Adjusted_FDR_Significant == TRUE,
            na.rm = TRUE
          ),
        
        Significant_both =
          sum(
            x$Unadjusted_FDR_Significant == TRUE &
              x$Adjusted_FDR_Significant == TRUE,
            na.rm = TRUE
          ),
        
        Significance_lost =
          sum(
            x$Unadjusted_FDR_Significant == TRUE &
              x$Adjusted_FDR_Significant == FALSE,
            na.rm = TRUE
          ),
        
        Significance_gained =
          sum(
            x$Unadjusted_FDR_Significant == FALSE &
              x$Adjusted_FDR_Significant == TRUE,
            na.rm = TRUE
          ),
        
        Direction_switches =
          sum(
            x$Direction_Switch == TRUE,
            na.rm = TRUE
          ),
        
        # ACTUAL STEP 8C-6 LABELS
        Stable_positive =
          sum(
            x$Direction_Transition ==
              "Positive -> Positive",
            na.rm = TRUE
          ),
        
        Stable_negative =
          sum(
            x$Direction_Transition ==
              "Negative -> Negative",
            na.rm = TRUE
          ),
        
        Positive_to_negative =
          sum(
            x$Direction_Transition ==
              "Positive -> Negative",
            na.rm = TRUE
          ),
        
        Negative_to_positive =
          sum(
            x$Direction_Transition ==
              "Negative -> Positive",
            na.rm = TRUE
          ),
        
        Positive_to_negative_adjusted_significant =
          sum(
            x$Direction_Transition ==
              "Positive -> Negative" &
              x$Adjusted_FDR_Significant == TRUE,
            na.rm = TRUE
          ),
        
        Adjusted_positive =
          sum(
            x$Adjusted_Direction == "Positive",
            na.rm = TRUE
          ),
        
        Adjusted_negative =
          sum(
            x$Adjusted_Direction == "Negative",
            na.rm = TRUE
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(core_summary_locked) <- NULL


cat("============================================================\n")
cat("LOCKED CORE-GENE SUMMARY\n")
cat("============================================================\n\n")

print(
  core_summary_locked,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 14. EXPECTED CORE-GENE SUMMARY
# ============================================================

expected_core_summary <- data.frame(
  
  Core_Gene =
    c(
      "AOC3",
      "MS4A7",
      "IL7R"
    ),
  
  N_pairs =
    c(
      18,
      18,
      18
    ),
  
  Unadjusted_FDR_significant =
    c(
      14,
      17,
      17
    ),
  
  Adjusted_FDR_significant =
    c(
      13,
      12,
      13
    ),
  
  Significant_both =
    c(
      9,
      12,
      13
    ),
  
  Significance_lost =
    c(
      5,
      5,
      4
    ),
  
  Significance_gained =
    c(
      4,
      0,
      0
    ),
  
  Direction_switches =
    c(
      13,
      12,
      7
    ),
  
  Stable_positive =
    c(
      4,
      6,
      11
    ),
  
  Stable_negative =
    c(
      1,
      0,
      0
    ),
  
  Positive_to_negative =
    c(
      13,
      12,
      7
    ),
  
  Negative_to_positive =
    c(
      0,
      0,
      0
    ),
  
  Positive_to_negative_adjusted_significant =
    c(
      9,
      8,
      4
    ),
  
  Adjusted_positive =
    c(
      4,
      6,
      11
    ),
  
  Adjusted_negative =
    c(
      14,
      12,
      7
    ),
  
  stringsAsFactors = FALSE
)


# ============================================================
# 15. EXACT CORE SUMMARY QC
# ============================================================

core_summary_check <- isTRUE(
  all.equal(
    core_summary_locked,
    expected_core_summary,
    check.attributes = TRUE
  )
)

cat(
  "Core-gene summary exact match:",
  core_summary_check,
  "\n\n"
)

stopifnot(
  core_summary_check
)

cat("Core-gene summary QC: PASS\n\n")


# ============================================================
# 16. GLOBAL ↔ CORE SUMMARY CONSISTENCY
# ============================================================

stopifnot(
  
  sum(
    core_summary_locked$N_pairs
  ) == LOCK_54,
  
  sum(
    core_summary_locked$Unadjusted_FDR_significant
  ) == LOCK_48,
  
  sum(
    core_summary_locked$Adjusted_FDR_significant
  ) == LOCK_38,
  
  sum(
    core_summary_locked$Significant_both
  ) == LOCK_34,
  
  sum(
    core_summary_locked$Significance_lost
  ) == LOCK_14,
  
  sum(
    core_summary_locked$Significance_gained
  ) == LOCK_4,
  
  sum(
    core_summary_locked$Direction_switches
  ) == LOCK_32,
  
  sum(
    core_summary_locked$
      Positive_to_negative_adjusted_significant
  ) == LOCK_21
)

cat("Global ↔ core-gene consistency: PASS\n\n")


# ============================================================
# 17. LOCK THE 21 POSITIVE → NEGATIVE SIGNIFICANT PAIRS
# ============================================================

switch_and_significant_21 <- complete_54[
  complete_54$Direction_Transition ==
    "Positive -> Negative" &
    complete_54$Adjusted_FDR_Significant == TRUE,
  ,
  drop = FALSE
]

stopifnot(
  nrow(switch_and_significant_21) == 21
)

cat("Positive -> negative switches remaining significant:\n")
cat(
  "N = ",
  nrow(switch_and_significant_21),
  "\n\n",
  sep = ""
)

print(
  switch_and_significant_21[
    ,
    c(
      "Core_Gene",
      "Checkpoint_Gene",
      "N",
      "Unadjusted_Spearman_rho",
      "Unadjusted_FDR_BH",
      "Adjusted_Partial_Spearman_rho",
      "Adjusted_FDR"
    )
  ],
  row.names = FALSE
)

cat("\n")
cat("21-pair manuscript switch set: PASS\n\n")


# ============================================================
# 18. IMPORTANT DISTINCTION BETWEEN THE TWO '21' COUNTS
# ============================================================

adjusted_positive <- sum(
  complete_54$Adjusted_Direction == "Positive",
  na.rm = TRUE
)

adjusted_negative <- sum(
  complete_54$Adjusted_Direction == "Negative",
  na.rm = TRUE
)

stopifnot(
  adjusted_positive == 21,
  adjusted_negative == 33,
  adjusted_positive + adjusted_negative == 54
)

cat("Adjusted direction totals:\n")
cat(
  "  Adjusted positive = ",
  adjusted_positive,
  "\n",
  sep = ""
)

cat(
  "  Adjusted negative = ",
  adjusted_negative,
  "\n\n",
  sep = ""
)

cat(
  "IMPORTANT:\n",
  "  21 adjusted-positive associations\n",
  "  !=\n",
  "  21 positive-to-negative switches remaining significant.\n",
  "\n",
  sep = ""
)

cat("Adjusted-direction QC: PASS\n\n")


# ============================================================
# 19. CREATE FINAL RESULT-DEFINITION LOCK OBJECT
# ============================================================

result_definition_lock <- list(
  
  metadata = list(
    
    workflow =
      "Paper 1 — Workflow 08",
    
    step =
      "STEP 8C-7A",
    
    purpose =
      "Manuscript result-definition lock",
    
    source =
      basename(source_file),
    
    source_status =
      "LOCKED STEP 8C-6",
    
    recalculation_performed =
      FALSE,
    
    note =
      "All values extracted from locked STEP 8C-6 results."
  ),
  
  
  manuscript_lock =
    manuscript_lock,
  
  
  core_gene_summary =
    core_summary_locked,
  
  
  direction_summary =
    data.frame(
      
      Direction_Transition =
        c(
          "Positive -> Positive",
          "Negative -> Negative",
          "Positive -> Negative",
          "Negative -> Positive"
        ),
      
      N =
        c(
          positive_to_positive,
          negative_to_negative,
          positive_to_negative,
          negative_to_positive
        ),
      
      stringsAsFactors = FALSE
    ),
  
  
  significance_summary =
    data.frame(
      
      Transition =
        c(
          "Significant both",
          "Significance lost",
          "Significance gained",
          "Not significant both"
        ),
      
      N =
        c(
          34,
          14,
          4,
          n_neither
        ),
      
      stringsAsFactors = FALSE
    ),
  
  
  switch_and_significant_21 =
    switch_and_significant_21,
  
  
  qc =
    list(
      
      total_pairs =
        LOCK_54,
      
      unadjusted_significant =
        LOCK_48,
      
      adjusted_significant =
        LOCK_38,
      
      significant_both =
        LOCK_34,
      
      significance_lost =
        LOCK_14,
      
      significance_gained =
        LOCK_4,
      
      direction_switches =
        LOCK_32,
      
      positive_to_negative_adjusted_significant =
        LOCK_21,
      
      adjusted_positive =
        adjusted_positive,
      
      adjusted_negative =
        adjusted_negative,
      
      N_all_503 =
        all(
          complete_54$N == 503
        ),
      
      pair_unique =
        length(
          unique(pair_id)
        ) == 54,
      
      core_genes_correct =
        setequal(
          unique(
            complete_54$Core_Gene
          ),
          expected_core_genes
        ),
      
      checkpoint_gene_count =
        length(
          unique(
            complete_54$Checkpoint_Gene
          )
        ) == 18,
      
      arithmetic_consistency =
        (
          LOCK_48 ==
            LOCK_34 + LOCK_14
          &
            LOCK_38 ==
            LOCK_34 + LOCK_4
        ),
      
      direction_consistency =
        (
          positive_to_positive == 21 &
            negative_to_negative == 1 &
            positive_to_negative == 32 &
            negative_to_positive == 0
        ),
      
      core_summary_consistency =
        core_summary_check
    )
)


# ============================================================
# 20. FINAL QC
# ============================================================

final_qc <- all(
  
  result_definition_lock$qc$total_pairs == 54,
  
  result_definition_lock$qc$unadjusted_significant == 48,
  
  result_definition_lock$qc$adjusted_significant == 38,
  
  result_definition_lock$qc$significant_both == 34,
  
  result_definition_lock$qc$significance_lost == 14,
  
  result_definition_lock$qc$significance_gained == 4,
  
  result_definition_lock$qc$direction_switches == 32,
  
  result_definition_lock$
    qc$positive_to_negative_adjusted_significant == 21,
  
  result_definition_lock$qc$adjusted_positive == 21,
  
  result_definition_lock$qc$adjusted_negative == 33,
  
  result_definition_lock$qc$N_all_503,
  
  result_definition_lock$qc$pair_unique,
  
  result_definition_lock$qc$core_genes_correct,
  
  result_definition_lock$qc$checkpoint_gene_count,
  
  result_definition_lock$qc$arithmetic_consistency,
  
  result_definition_lock$qc$direction_consistency,
  
  result_definition_lock$qc$core_summary_consistency
)


cat("============================================================\n")
cat("STEP 8C-7A FINAL QC\n")
cat("============================================================\n\n")

cat(
  "54 total pairs                        : ",
  LOCK_54,
  "\n",
  
  "48 unadjusted FDR-significant        : ",
  LOCK_48,
  "\n",
  
  "38 adjusted FDR-significant          : ",
  LOCK_38,
  "\n",
  
  "34 significant both                 : ",
  LOCK_34,
  "\n",
  
  "14 significance lost                : ",
  LOCK_14,
  "\n",
  
  " 4 significance gained              : ",
  LOCK_4,
  "\n",
  
  "32 direction switches               : ",
  LOCK_32,
  "\n",
  
  "21 Pos -> Neg & adjusted significant : ",
  LOCK_21,
  "\n",
  sep = ""
)

cat("\n")


if (!final_qc) {
  stop(
    "STEP 8C-7A FINAL QC FAILED."
  )
}

cat("QC STATUS: PASS\n")
cat(
  "STEP 8C-7A RESULT-DEFINITION LOCK: ",
  "SUCCESSFULLY LOCKED\n",
  sep = ""
)

cat(
  "No new correlation calculation was performed.\n"
)

cat("============================================================\n\n")


# ============================================================
# 21. SAVE MASTER RDS
# ============================================================

output_rds <- file.path(
  workdir,
  "Workflow-08-STEP8C7A_manuscript_result_definition_lock.rds"
)

saveRDS(
  result_definition_lock,
  output_rds
)

stopifnot(
  file.exists(output_rds)
)

cat("Master RDS saved:\n")
cat(output_rds, "\n\n")


# ============================================================
# 22. SAVE MANUSCRIPT LOCK TABLE
# ============================================================

output_lock_tsv <- file.path(
  workdir,
  "Workflow-08-STEP8C7A_manuscript_result_definition_lock.tsv"
)

write.table(
  manuscript_lock,
  output_lock_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

stopifnot(
  file.exists(output_lock_tsv)
)

cat("Manuscript lock TSV saved:\n")
cat(output_lock_tsv, "\n\n")


# ============================================================
# 23. SAVE CORE-GENE SUMMARY
# ============================================================

output_core_tsv <- file.path(
  workdir,
  "Workflow-08-STEP8C7A_core_gene_summary_locked.tsv"
)

write.table(
  core_summary_locked,
  output_core_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

stopifnot(
  file.exists(output_core_tsv)
)

cat("Core-gene summary TSV saved:\n")
cat(output_core_tsv, "\n\n")


# ============================================================
# 24. SAVE 21-PAIR SWITCH SET
# ============================================================

output_switch21_tsv <- file.path(
  workdir,
  "Workflow-08-STEP8C7A_positive_to_negative_significant_21.tsv"
)

write.table(
  switch_and_significant_21,
  output_switch21_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

stopifnot(
  file.exists(output_switch21_tsv)
)

cat("21-pair switch set TSV saved:\n")
cat(output_switch21_tsv, "\n\n")


# ============================================================
# 25. RELOAD TEST
# ============================================================

cat("Testing RDS reload...\n\n")

test_reload <- readRDS(
  output_rds
)

stopifnot(
  
  identical(
    test_reload$manuscript_lock,
    manuscript_lock
  ),
  
  identical(
    test_reload$core_gene_summary,
    core_summary_locked
  ),
  
  nrow(
    test_reload$switch_and_significant_21
  ) == 21
)

cat("RDS reload test: PASS\n\n")


# ============================================================
# 26. FINAL SUMMARY
# ============================================================

cat("============================================================\n")
cat("STEP 8C-7A COMPLETED\n")
cat("============================================================\n\n")

cat("Locked manuscript numbers:\n")
cat("  54 = total predefined associations\n")
cat("  48 = unadjusted FDR < 0.05\n")
cat("  38 = adjusted FDR < 0.05\n")
cat("  34 = significant both\n")
cat("  14 = significance lost\n")
cat("   4 = significance gained\n")
cat("  32 = direction switches\n")
cat("  21 = positive -> negative switches remaining significant\n\n")

cat("Locked direction strata:\n")
cat("  Positive -> Positive = 21\n")
cat("  Negative -> Negative = 1\n")
cat("  Positive -> Negative = 32\n")
cat("  Negative -> Positive = 0\n\n")

cat("Locked core-gene summary:\n")
print(
  core_summary_locked,
  row.names = FALSE
)

cat("\n")

cat("IMPORTANT:\n")
cat(
  "  Adjusted-positive total = 21\n",
  "  Positive-to-negative switches remaining significant = 21\n",
  "  These are different manuscript definitions.\n\n",
  sep = ""
)

cat("QC STATUS: PASS\n")
cat(
  "All manuscript-level values are locked from STEP 8C-6.\n"
)

cat("============================================================\n")


# ============================================================
# ============================================================
# Paper 1
# STEP 8C-7B
# Supplementary Table generation
#
# CORRECTED VERSION
#
# Source of truth:
#   step8c6$complete_54
#
# IMPORTANT:
#   No statistical recalculation.
#   No new correlation.
#   No new FDR correction.
#   No new significance testing.
#
# All manuscript-level definitions are inherited from
# STEP 8C-6 locked fields.
# ============================================================


cat("\n============================================================\n")
cat("STEP 8C-7B: Supplementary Table generation\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "readr",
  "openxlsx"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    stop(
      paste0(
        "Required package not installed: ",
        pkg
      )
    )
    
  }
  
}


# ------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------

analysis_dir <-
  "O:/Rproject/Paper1_tumor-only WGCNA"


checkpoint_file <- file.path(
  analysis_dir,
  "Workflow-08-STEP8C6_manuscript_checkpoint_extraction.rds"
)


output_xlsx <- file.path(
  analysis_dir,
  "Workflow-08-STEP8C7B_Supplementary_Table_Immune_Checkpoint_54.xlsx"
)


output_csv <- file.path(
  analysis_dir,
  "Workflow-08-STEP8C7B_Supplementary_Table_Immune_Checkpoint_54.csv"
)


output_rds <- file.path(
  analysis_dir,
  "Workflow-08-STEP8C7B_Supplementary_Table_Immune_Checkpoint_54.rds"
)


cat("Analysis directory:\n")
cat(analysis_dir, "\n\n")

cat("Checkpoint file:\n")
cat(checkpoint_file, "\n\n")


# ------------------------------------------------------------
# 3. Check checkpoint
# ------------------------------------------------------------

if (!file.exists(checkpoint_file)) {
  
  stop(
    paste0(
      "STEP 8C-6 checkpoint not found:\n",
      checkpoint_file
    )
  )
  
}


# ------------------------------------------------------------
# 4. Load locked STEP 8C-6 checkpoint
# ------------------------------------------------------------

step8c6 <- readRDS(
  checkpoint_file
)


cat(
  "STEP 8C-6 checkpoint loaded successfully.\n\n"
)


# ------------------------------------------------------------
# 5. Confirm complete_54
# ------------------------------------------------------------

if (!"complete_54" %in% names(step8c6)) {
  
  stop(
    "step8c6$complete_54 not found."
  )
  
}


complete_54 <- step8c6$complete_54


cat("Source object: step8c6$complete_54\n")

cat("Dimensions:\n")

print(
  dim(complete_54)
)

cat("\n")


# ------------------------------------------------------------
# 6. Structural QC
# ------------------------------------------------------------

if (!is.data.frame(complete_54)) {
  
  stop(
    "complete_54 is not a data.frame."
  )
  
}


if (nrow(complete_54) != 54) {
  
  stop(
    paste0(
      "Expected 54 rows; found ",
      nrow(complete_54)
    )
  )
  
}


if (ncol(complete_54) != 26) {
  
  warning(
    paste0(
      "Expected 26 columns from STEP 8C-6; found ",
      ncol(complete_54)
    )
  )
  
}


cat("QC PASS: complete_54 contains 54 associations.\n\n")


# ------------------------------------------------------------
# 7. Confirm exact locked columns
# ------------------------------------------------------------

required_columns <- c(
  
  "Core_Gene",
  "Checkpoint_Gene",
  "N",
  
  "Unadjusted_Spearman_rho",
  "Unadjusted_P",
  "Unadjusted_FDR_BH",
  "Unadjusted_Direction",
  
  "Adjusted_Partial_Spearman_rho",
  "Adjusted_P",
  "Adjusted_FDR",
  "Adjusted_Direction",
  
  "Adjustment",
  
  "Delta_rho",
  "Absolute_Delta_rho",
  "Absolute_rho_change",
  
  "Direction_Transition",
  "Direction_Switch",
  
  "Unadjusted_FDR_Significant",
  "Adjusted_FDR_Significant",
  
  "FDR_Transition",
  
  "Audit_Class",
  "Sig_Unadjusted",
  "Sig_Adjusted",
  
  "Association_Stratum",
  "Significance_Stratum",
  "Manuscript_Class"
  
)


missing_columns <- setdiff(
  required_columns,
  names(complete_54)
)


if (length(missing_columns) > 0) {
  
  stop(
    paste0(
      "Required STEP 8C-6 columns missing:\n",
      paste(
        missing_columns,
        collapse = "\n"
      )
    )
  )
  
}


cat(
  "QC PASS: all 26 STEP 8C-6 columns confirmed.\n\n"
)


# ------------------------------------------------------------
# 8. Preserve the locked source object
# ------------------------------------------------------------

supp_table <- complete_54


# ------------------------------------------------------------
# 9. Explicitly verify core genes
# ------------------------------------------------------------

core_gene_counts <- table(
  supp_table$Core_Gene
)


cat("Core-gene counts:\n")

print(
  core_gene_counts
)

cat("\n")


expected_core_counts <- c(
  AOC3 = 18,
  MS4A7 = 18,
  IL7R = 18
)


for (gene in names(expected_core_counts)) {
  
  observed <- as.numeric(
    core_gene_counts[gene]
  )
  
  expected <- expected_core_counts[gene]
  
  if (
    is.na(observed) ||
    observed != expected
  ) {
    
    stop(
      paste0(
        "Core-gene QC failed for ",
        gene,
        ": expected ",
        expected,
        ", observed ",
        observed
      )
    )
    
  }
  
}


cat(
  "QC PASS: AOC3 = 18, MS4A7 = 18, IL7R = 18.\n\n"
)


# ------------------------------------------------------------
# 10. Do NOT recalculate statistical fields
#
# These fields are already locked in STEP 8C-6:
#
#   Unadjusted_FDR_BH
#   Adjusted_FDR
#   Direction_Transition
#   Direction_Switch
#   Unadjusted_FDR_Significant
#   Adjusted_FDR_Significant
#   FDR_Transition
#   Manuscript_Class
#
# They are copied exactly.
# ------------------------------------------------------------


# ------------------------------------------------------------
# 11. Manuscript-friendly display table
#
# Preserve all 26 locked columns.
# Add human-readable significance labels only if needed.
# ------------------------------------------------------------

supp_table_export <- supp_table


# ------------------------------------------------------------
# 12. Create explicit display labels
#
# These are NOT new statistical calculations.
# They simply translate existing locked TRUE/FALSE
# fields into Yes/No for Excel.
# ------------------------------------------------------------

supp_table_export$Unadjusted_FDR_Significant_Display <-
  ifelse(
    supp_table$Unadjusted_FDR_Significant,
    "Yes",
    "No"
  )


supp_table_export$Adjusted_FDR_Significant_Display <-
  ifelse(
    supp_table$Adjusted_FDR_Significant,
    "Yes",
    "No"
  )


# ------------------------------------------------------------
# 13. Put the most important manuscript fields first
# ------------------------------------------------------------

first_columns <- c(
  
  "Core_Gene",
  "Checkpoint_Gene",
  "N",
  
  "Unadjusted_Spearman_rho",
  "Unadjusted_P",
  "Unadjusted_FDR_BH",
  "Unadjusted_Direction",
  
  "Adjusted_Partial_Spearman_rho",
  "Adjusted_P",
  "Adjusted_FDR",
  "Adjusted_Direction",
  
  "Direction_Transition",
  "Direction_Switch",
  
  "Unadjusted_FDR_Significant",
  "Adjusted_FDR_Significant",
  
  "FDR_Transition",
  "Manuscript_Class"
  
)


remaining_columns <- setdiff(
  names(supp_table_export),
  first_columns
)


supp_table_export <- supp_table_export[
  c(
    first_columns,
    remaining_columns
  )
]


# ------------------------------------------------------------
# 14. Order by core gene and checkpoint gene
# ------------------------------------------------------------

core_order <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


supp_table_export$Core_Gene <- factor(
  supp_table_export$Core_Gene,
  levels = core_order
)


supp_table_export <- supp_table_export[
  order(
    supp_table_export$Core_Gene,
    supp_table_export$Checkpoint_Gene
  ),
  ,
  drop = FALSE
]


supp_table_export$Core_Gene <-
  as.character(
    supp_table_export$Core_Gene
  )


# ------------------------------------------------------------
# 15. Reproduce STEP 8C-7A manuscript numbers
#
# IMPORTANT:
#   These are QC checks only.
#   Nothing is recalculated statistically.
# ------------------------------------------------------------


n_total <- nrow(
  supp_table
)


n_unadjusted <- sum(
  supp_table$Unadjusted_FDR_Significant,
  na.rm = TRUE
)


n_adjusted <- sum(
  supp_table$Adjusted_FDR_Significant,
  na.rm = TRUE
)


n_both <- sum(
  supp_table$Unadjusted_FDR_Significant &
    supp_table$Adjusted_FDR_Significant,
  na.rm = TRUE
)


n_lost <- sum(
  supp_table$FDR_Transition ==
    "Significance lost",
  na.rm = TRUE
)


n_gained <- sum(
  supp_table$FDR_Transition ==
    "Significance gained",
  na.rm = TRUE
)


n_direction_switch <- sum(
  supp_table$Direction_Switch,
  na.rm = TRUE
)


n_positive_negative_adjusted <- sum(
  supp_table$Direction_Transition ==
    "Positive -> Negative" &
    supp_table$Adjusted_FDR_Significant,
  na.rm = TRUE
)


# ------------------------------------------------------------
# 16. Print manuscript QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 8C-7A MANUSCRIPT NUMBER QC\n")
cat("============================================================\n\n")


cat(
  "54 = total predefined associations: ",
  n_total,
  "\n",
  sep = ""
)


cat(
  "48 = unadjusted FDR < 0.05: ",
  n_unadjusted,
  "\n",
  sep = ""
)


cat(
  "38 = adjusted FDR < 0.05: ",
  n_adjusted,
  "\n",
  sep = ""
)


cat(
  "34 = significant both: ",
  n_both,
  "\n",
  sep = ""
)


cat(
  "14 = significance lost: ",
  n_lost,
  "\n",
  sep = ""
)


cat(
  "4 = significance gained: ",
  n_gained,
  "\n",
  sep = ""
)


cat(
  "32 = direction switches: ",
  n_direction_switch,
  "\n",
  sep = ""
)


cat(
  "21 = positive -> negative switches remaining adjusted-significant: ",
  n_positive_negative_adjusted,
  "\n\n",
  sep = ""
)

# ============================================================
# STEP 8C-7B diagnostic
# Inspect the actual FDR_Transition labels
# ============================================================

cat("\n============================================================\n")
cat("DIAGNOSTIC: FDR_Transition\n")
cat("============================================================\n\n")

cat("Unique FDR_Transition values:\n")
print(
  unique(supp_table$FDR_Transition)
)

cat("\nFrequency table:\n")
print(
  table(
    supp_table$FDR_Transition,
    useNA = "ifany"
  )
)

cat("\nCross-tabulation with significance fields:\n")

print(
  table(
    Unadjusted = supp_table$Unadjusted_FDR_Significant,
    Adjusted   = supp_table$Adjusted_FDR_Significant,
    FDR_Transition = supp_table$FDR_Transition,
    useNA = "ifany"
  )
)

cat("\nRows where unadjusted significant = TRUE and adjusted = FALSE:\n")

print(
  supp_table[
    supp_table$Unadjusted_FDR_Significant &
      !supp_table$Adjusted_FDR_Significant,
    c(
      "Core_Gene",
      "Checkpoint_Gene",
      "Unadjusted_FDR_BH",
      "Adjusted_FDR",
      "Unadjusted_FDR_Significant",
      "Adjusted_FDR_Significant",
      "FDR_Transition"
    ),
    drop = FALSE
  ]
)

cat("\nRows where unadjusted significant = FALSE and adjusted = TRUE:\n")

print(
  supp_table[
    !supp_table$Unadjusted_FDR_Significant &
      supp_table$Adjusted_FDR_Significant,
    c(
      "Core_Gene",
      "Checkpoint_Gene",
      "Unadjusted_FDR_BH",
      "Adjusted_FDR",
      "Unadjusted_FDR_Significant",
      "Adjusted_FDR_Significant",
      "FDR_Transition"
    ),
    drop = FALSE
  ]
)


##STEP 8C-7B — Final Supplementary Table + QC Lock

# ============================================================
# Paper 1 — STEP 8C-7B
# Final Supplementary Table + Manuscript QC Lock
#
# Purpose:
#   1. Use ONLY locked STEP 8C-6 results
#   2. Generate the 54-association supplementary table
#   3. Reproduce all STEP 8C-7A manuscript-level locked values
#   4. Generate core-gene QC summaries
#   5. Export Excel / CSV / RDS
#
# IMPORTANT:
#   NO new statistical analysis
#   NO new Spearman correlation
#   NO new partial Spearman correlation
#   NO new FDR calculation
#
# Sole source of truth:
#   Workflow-08-STEP8C6_manuscript_checkpoint_extraction.rds
# ============================================================


cat("\n============================================================\n")
cat("STEP 8C-7B: FINAL SUPPLEMENTARY TABLE + QC LOCK\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Project directory
# ------------------------------------------------------------

project_dir <- "input_dir/Paper1_tumor-only WGCNA"

setwd(project_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 2. Load locked STEP 8C-6 object
# ------------------------------------------------------------

input_file <- file.path(
  project_dir,
  "Workflow-08-STEP8C6_manuscript_checkpoint_extraction.rds"
)

stopifnot(file.exists(input_file))

step8c6 <- readRDS(input_file)

cat("Loaded:\n")
cat(basename(input_file), "\n\n")


# ------------------------------------------------------------
# 3. Extract locked 54-row manuscript table
# ------------------------------------------------------------

stopifnot("complete_54" %in% names(step8c6))

supp_table <- step8c6$complete_54


# ------------------------------------------------------------
# 4. Basic structural QC
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("4. Structural QC\n")
cat("------------------------------------------------------------\n\n")

cat("Dimensions:\n")
print(dim(supp_table))

cat("\nColumn names:\n")
print(names(supp_table))


expected_n <- 54

if (nrow(supp_table) != expected_n) {
  stop(
    "QC FAILED: complete_54 does not contain exactly 54 associations."
  )
}


required_columns <- c(
  "Core_Gene",
  "Checkpoint_Gene",
  "N",
  "Unadjusted_Spearman_rho",
  "Unadjusted_P",
  "Unadjusted_FDR_BH",
  "Unadjusted_Direction",
  "Adjusted_Partial_Spearman_rho",
  "Adjusted_P",
  "Adjusted_FDR",
  "Adjusted_Direction",
  "Adjustment",
  "Delta_rho",
  "Absolute_Delta_rho",
  "Absolute_rho_change",
  "Direction_Transition",
  "Direction_Switch",
  "Unadjusted_FDR_Significant",
  "Adjusted_FDR_Significant",
  "FDR_Transition",
  "Audit_Class",
  "Sig_Unadjusted",
  "Sig_Adjusted",
  "Association_Stratum",
  "Significance_Stratum",
  "Manuscript_Class"
)

missing_columns <- setdiff(
  required_columns,
  names(supp_table)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "QC FAILED: Missing required columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

cat("\nStructural QC: PASS\n")


# ------------------------------------------------------------
# 5. Confirm core-gene composition
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("5. Core-gene composition\n")
cat("------------------------------------------------------------\n\n")

core_counts <- table(supp_table$Core_Gene)

print(core_counts)

expected_core_counts <- c(
  AOC3   = 18,
  MS4A7  = 18,
  IL7R   = 18
)

for (g in names(expected_core_counts)) {
  
  observed <- ifelse(
    g %in% names(core_counts),
    as.integer(core_counts[g]),
    0L
  )
  
  if (observed != expected_core_counts[g]) {
    
    stop(
      paste0(
        "QC FAILED: ",
        g,
        " expected ",
        expected_core_counts[g],
        " associations but observed ",
        observed,
        "."
      )
    )
  }
}

cat("\nCore-gene composition QC: PASS\n")


# ------------------------------------------------------------
# 6. Locked manuscript-level QC
#
# IMPORTANT:
#   lost/gained are derived from the already locked Boolean
#   significance fields, NOT from text labels.
# ------------------------------------------------------------

n_total <- nrow(supp_table)

n_unadjusted <- sum(
  supp_table$Unadjusted_FDR_Significant,
  na.rm = TRUE
)

n_adjusted <- sum(
  supp_table$Adjusted_FDR_Significant,
  na.rm = TRUE
)

n_both <- sum(
  supp_table$Unadjusted_FDR_Significant &
    supp_table$Adjusted_FDR_Significant,
  na.rm = TRUE
)

n_lost <- sum(
  supp_table$Unadjusted_FDR_Significant &
    !supp_table$Adjusted_FDR_Significant,
  na.rm = TRUE
)

n_gained <- sum(
  !supp_table$Unadjusted_FDR_Significant &
    supp_table$Adjusted_FDR_Significant,
  na.rm = TRUE
)

n_direction_switch <- sum(
  supp_table$Direction_Switch,
  na.rm = TRUE
)

n_pos_to_neg_adj_sig <- sum(
  supp_table$Direction_Transition ==
    "Positive -> Negative" &
    supp_table$Adjusted_FDR_Significant,
  na.rm = TRUE
)


manuscript_qc <- data.frame(
  
  Metric = c(
    "Total predefined associations",
    "Unadjusted FDR < 0.05",
    "Adjusted FDR < 0.05",
    "Significant both",
    "Significance lost",
    "Significance gained",
    "Direction switches",
    "Positive -> Negative switches remaining adjusted-significant"
  ),
  
  Locked = c(
    54,
    48,
    38,
    34,
    14,
    4,
    32,
    21
  ),
  
  Observed = c(
    n_total,
    n_unadjusted,
    n_adjusted,
    n_both,
    n_lost,
    n_gained,
    n_direction_switch,
    n_pos_to_neg_adj_sig
  ),
  
  stringsAsFactors = FALSE
)

manuscript_qc$QC <- manuscript_qc$Locked ==
  manuscript_qc$Observed


cat("\n------------------------------------------------------------\n")
cat("MANUSCRIPT-LEVEL QC\n")
cat("------------------------------------------------------------\n\n")

print(manuscript_qc)


if (!all(manuscript_qc$QC)) {
  
  print(
    manuscript_qc[
      !manuscript_qc$QC,
      ,
      drop = FALSE
    ]
  )
  
  stop(
    "QC FAILED: STEP 8C-7A locked manuscript values are not reproduced."
  )
}

cat("\nMANUSCRIPT-LEVEL QC: PASS\n")


# ------------------------------------------------------------
# 7. Direction-stratum QC
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("7. Direction-stratum QC\n")
cat("------------------------------------------------------------\n\n")

direction_table <- table(
  supp_table$Direction_Transition,
  useNA = "ifany"
)

print(direction_table)


expected_direction <- c(
  "Positive -> Positive" = 21,
  "Negative -> Negative" = 1,
  "Positive -> Negative" = 32,
  "Negative -> Positive" = 0
)


for (nm in names(expected_direction)) {
  
  observed <- ifelse(
    nm %in% names(direction_table),
    as.integer(direction_table[nm]),
    0L
  )
  
  if (observed != expected_direction[nm]) {
    
    stop(
      paste0(
        "QC FAILED: Direction stratum '",
        nm,
        "' expected ",
        expected_direction[nm],
        " but observed ",
        observed,
        "."
      )
    )
  }
}

cat("\nDirection-stratum QC: PASS\n")


# ------------------------------------------------------------
# 8. Core-gene manuscript summary
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


core_summary_list <- lapply(
  core_genes,
  function(g) {
    
    d <- supp_table[
      supp_table$Core_Gene == g,
      ,
      drop = FALSE
    ]
    
    data.frame(
      
      Core_Gene = g,
      
      N_pairs = nrow(d),
      
      Unadjusted_FDR_significant =
        sum(
          d$Unadjusted_FDR_Significant,
          na.rm = TRUE
        ),
      
      Adjusted_FDR_significant =
        sum(
          d$Adjusted_FDR_Significant,
          na.rm = TRUE
        ),
      
      Significant_both =
        sum(
          d$Unadjusted_FDR_Significant &
            d$Adjusted_FDR_Significant,
          na.rm = TRUE
        ),
      
      Significance_lost =
        sum(
          d$Unadjusted_FDR_Significant &
            !d$Adjusted_FDR_Significant,
          na.rm = TRUE
        ),
      
      Significance_gained =
        sum(
          !d$Unadjusted_FDR_Significant &
            d$Adjusted_FDR_Significant,
          na.rm = TRUE
        ),
      
      Direction_switches =
        sum(
          d$Direction_Switch,
          na.rm = TRUE
        ),
      
      Stable_positive =
        sum(
          d$Direction_Transition ==
            "Positive -> Positive",
          na.rm = TRUE
        ),
      
      Stable_negative =
        sum(
          d$Direction_Transition ==
            "Negative -> Negative",
          na.rm = TRUE
        ),
      
      Positive_to_negative =
        sum(
          d$Direction_Transition ==
            "Positive -> Negative",
          na.rm = TRUE
        ),
      
      Negative_to_positive =
        sum(
          d$Direction_Transition ==
            "Negative -> Positive",
          na.rm = TRUE
        ),
      
      Positive_to_negative_adjusted_significant =
        sum(
          d$Direction_Transition ==
            "Positive -> Negative" &
            d$Adjusted_FDR_Significant,
          na.rm = TRUE
        ),
      
      Adjusted_positive =
        sum(
          d$Adjusted_Direction == "Positive",
          na.rm = TRUE
        ),
      
      Adjusted_negative =
        sum(
          d$Adjusted_Direction == "Negative",
          na.rm = TRUE
        ),
      
      stringsAsFactors = FALSE
    )
  }
)

core_gene_qc <- do.call(
  rbind,
  core_summary_list
)

# ------------------------------------------------------------
# 9. Corrected locked core-gene summary
#
# Direction-related values are derived from the locked
# STEP 8C-6 Direction_Transition field.
#
# Internal consistency requirement:
#   Stable_positive +
#   Stable_negative +
#   Positive_to_negative +
#   Negative_to_positive
#   = N_pairs
#
# Direction_switches =
#   Positive_to_negative +
#   Negative_to_positive
# ------------------------------------------------------------

expected_core_summary <- data.frame(
  
  Core_Gene = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  
  N_pairs = c(
    18,
    18,
    18
  ),
  
  Unadjusted_FDR_significant = c(
    14,
    17,
    17
  ),
  
  Adjusted_FDR_significant = c(
    13,
    12,
    13
  ),
  
  Significant_both = c(
    9,
    12,
    13
  ),
  
  Significance_lost = c(
    5,
    5,
    4
  ),
  
  Significance_gained = c(
    4,
    0,
    0
  ),
  
  Direction_switches = c(
    13,
    12,
    7
  ),
  
  Stable_positive = c(
    4,
    6,
    11
  ),
  
  Stable_negative = c(
    1,
    0,
    0
  ),
  
  Positive_to_negative = c(
    13,
    12,
    7
  ),
  
  Negative_to_positive = c(
    0,
    0,
    0
  ),
  
  Positive_to_negative_adjusted_significant = c(
    9,
    8,
    4
  ),
  
  Adjusted_positive = c(
    4,
    6,
    11
  ),
  
  Adjusted_negative = c(
    14,
    12,
    7
  ),
  
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 9A. Core-gene internal direction consistency
# ------------------------------------------------------------

core_direction_check <- within(
  expected_core_summary,
  {
    Direction_total <-
      Stable_positive +
      Stable_negative +
      Positive_to_negative +
      Negative_to_positive
    
    Direction_switch_check <-
      Direction_switches ==
      (
        Positive_to_negative +
          Negative_to_positive
      )
    
    Direction_total_check <-
      Direction_total == N_pairs
  }
)

cat("\n------------------------------------------------------------\n")
cat("CORE-GENE INTERNAL DIRECTION CONSISTENCY\n")
cat("------------------------------------------------------------\n\n")

print(
  core_direction_check[
    ,
    c(
      "Core_Gene",
      "N_pairs",
      "Stable_positive",
      "Stable_negative",
      "Positive_to_negative",
      "Negative_to_positive",
      "Direction_switches",
      "Direction_total",
      "Direction_total_check",
      "Direction_switch_check"
    )
  ]
)

if (
  !all(core_direction_check$Direction_total_check) ||
  !all(core_direction_check$Direction_switch_check)
) {
  
  stop(
    "QC FAILED: Core-gene direction summary is internally inconsistent."
  )
}

cat("\nCore-gene internal direction QC: PASS\n")



# ------------------------------------------------------------
# 10. Core-gene QC comparison
# ------------------------------------------------------------

core_qc <- merge(
  expected_core_summary,
  core_gene_qc,
  by = "Core_Gene",
  suffixes = c("_Locked", "_Observed"),
  sort = FALSE
)


comparison_columns <- setdiff(
  names(expected_core_summary),
  "Core_Gene"
)


core_qc$QC <- apply(
  core_qc[
    ,
    paste0(comparison_columns, "_Locked"),
    drop = FALSE
  ],
  1,
  function(x) TRUE
)


# Recalculate QC explicitly column by column
core_qc$QC <- TRUE

for (v in comparison_columns) {
  
  locked_col <- paste0(v, "_Locked")
  observed_col <- paste0(v, "_Observed")
  
  core_qc$QC <- core_qc$QC &
    (
      core_qc[[locked_col]] ==
        core_qc[[observed_col]]
    )
}


cat("\n------------------------------------------------------------\n")
cat("CORE-GENE QC\n")
cat("------------------------------------------------------------\n\n")

print(core_qc)


if (!all(core_qc$QC)) {
  
  print(
    core_qc[
      !core_qc$QC,
      ,
      drop = FALSE
    ]
  )
  
  stop(
    "QC FAILED: Core-gene locked summaries are not reproduced."
  )
}

cat("\nCore-gene QC: PASS\n")


# ------------------------------------------------------------
# 11. Final Supplementary Table
#
# Keep the original 26 locked columns unchanged.
# ------------------------------------------------------------

supplementary_table <- supp_table


# ------------------------------------------------------------
# 12. Add a clean manuscript-facing row number
# ------------------------------------------------------------

supplementary_table$Association_ID <-
  seq_len(nrow(supplementary_table))


supplementary_table <- supplementary_table[
  ,
  c(
    "Association_ID",
    setdiff(
      names(supplementary_table),
      "Association_ID"
    )
  )
]


# ------------------------------------------------------------
# 13. Final table ordering
#
# Preserve original STEP 8C-6 ordering.
# No sorting by significance or effect size.
# ------------------------------------------------------------

stopifnot(
  identical(
    supplementary_table$Association_ID,
    seq_len(54)
  )
)


# ------------------------------------------------------------
# 14. Export filenames
# ------------------------------------------------------------

xlsx_file <- file.path(
  project_dir,
  "Workflow-08-STEP8C7B_Supplementary_Table_Immune_Checkpoint_54.xlsx"
)

csv_file <- file.path(
  project_dir,
  "Workflow-08-STEP8C7B_Supplementary_Table_Immune_Checkpoint_54.csv"
)

rds_file <- file.path(
  project_dir,
  "Workflow-08-STEP8C7B_Supplementary_Table_Immune_Checkpoint_54.rds"
)


# ------------------------------------------------------------
# 15. Export CSV
# ------------------------------------------------------------

write.csv(
  supplementary_table,
  csv_file,
  row.names = FALSE,
  na = ""
)


# ------------------------------------------------------------
# 16. Prepare Excel workbook
# ------------------------------------------------------------

if (!requireNamespace(
  "openxlsx",
  quietly = TRUE
)) {
  stop(
    "Package 'openxlsx' is required for Excel export."
  )
}


wb <- openxlsx::createWorkbook()


# ------------------------------------------------------------
# Sheet 1: Supplementary_Table
# ------------------------------------------------------------

openxlsx::addWorksheet(
  wb,
  "Supplementary_Table"
)

openxlsx::writeData(
  wb,
  "Supplementary_Table",
  supplementary_table,
  withFilter = TRUE
)

openxlsx::freezePane(
  wb,
  "Supplementary_Table",
  firstRow = TRUE
)


# ------------------------------------------------------------
# Sheet 2: Manuscript_QC
# ------------------------------------------------------------

openxlsx::addWorksheet(
  wb,
  "Manuscript_QC"
)

openxlsx::writeData(
  wb,
  "Manuscript_QC",
  manuscript_qc
)


# ------------------------------------------------------------
# Sheet 3: Core_Gene_QC
# ------------------------------------------------------------

openxlsx::addWorksheet(
  wb,
  "Core_Gene_QC"
)

openxlsx::writeData(
  wb,
  "Core_Gene_QC",
  core_gene_qc
)


# ------------------------------------------------------------
# Sheet 4: Locked_Core_QC
# ------------------------------------------------------------

openxlsx::addWorksheet(
  wb,
  "Locked_Core_QC"
)

openxlsx::writeData(
  wb,
  "Locked_Core_QC",
  core_qc
)


# ------------------------------------------------------------
# Basic Excel formatting
# ------------------------------------------------------------

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center"
)


for (sheet in c(
  "Supplementary_Table",
  "Manuscript_QC",
  "Core_Gene_QC",
  "Locked_Core_QC"
)) {
  
  openxlsx::addStyle(
    wb,
    sheet,
    header_style,
    rows = 1,
    cols = 1:ncol(
      openxlsx::read.xlsx(
        wb,
        sheet = sheet,
        rows = 1
      )
    ),
    gridExpand = TRUE
  )
  
  openxlsx::setColWidths(
    wb,
    sheet,
    cols = 1:100,
    widths = "auto"
  )
}


# ------------------------------------------------------------
# Save Excel
# ------------------------------------------------------------

openxlsx::saveWorkbook(
  wb,
  xlsx_file,
  overwrite = TRUE
)


# ------------------------------------------------------------
# 17. Save final RDS
#
# Store both table and QC objects for future reproducibility.
# ------------------------------------------------------------

step8c7b_final <- list(
  
  supplementary_table = supplementary_table,
  
  manuscript_qc = manuscript_qc,
  
  core_gene_qc = core_gene_qc,
  
  locked_core_qc = core_qc,
  
  source_file = input_file,
  
  analysis_status = "LOCKED — NO NEW STATISTICAL ANALYSIS",
  
  workflow = "STEP 8C-7B",
  
  manuscript_locked_values = list(
    total = 54,
    unadjusted_fdr_significant = 48,
    adjusted_fdr_significant = 38,
    significant_both = 34,
    significance_lost = 14,
    significance_gained = 4,
    direction_switches = 32,
    positive_to_negative_adjusted_significant = 21
  )
)


saveRDS(
  step8c7b_final,
  rds_file
)


# ------------------------------------------------------------
# 18. Final file existence check
# ------------------------------------------------------------

output_files <- c(
  xlsx_file,
  csv_file,
  rds_file
)

file_check <- file.exists(output_files)

names(file_check) <- basename(output_files)

cat("\n------------------------------------------------------------\n")
cat("OUTPUT FILE CHECK\n")
cat("------------------------------------------------------------\n\n")

print(file_check)


if (!all(file_check)) {
  
  stop(
    "QC FAILED: One or more output files were not created."
  )
}


# ------------------------------------------------------------
# 19. Final summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 8C-7B FINAL STATUS: PASS\n")
cat("============================================================\n\n")

cat("Locked manuscript values:\n")
print(manuscript_qc)

cat("\nCore-gene summary:\n")
print(core_gene_qc)

cat("\nOutput files:\n")

for (f in output_files) {
  cat(" - ", f, "\n", sep = "")
}

cat("\n============================================================\n")
cat("STEP 8C-7B COMPLETE\n")
cat("No new statistical analysis was performed.\n")
cat("============================================================\n")


# ============================================================

# Paper 1
# STEP 8C-8
# Supplementary Figure 4
#
# ACTUAL FIGURE GENERATION + NUMERICAL QC
#
# IMPORTANT:
#   1. STEP 8C-7B results are LOCKED.
#   2. No correlation is recalculated.
#   3. No P value is recalculated.
#   4. No FDR is recalculated.
#   5. Panel A / B / C are generated separately.
#   6. Final figure assembly will be performed in Illustrator.
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)


# ============================================================
# 1. Project paths
# ============================================================

project_dir <- "input_dir/Paper1_tumor-only WGCNA"

locked_rds <- file.path(
  project_dir,
  "Workflow-08-STEP8C7B_Supplementary_Table_Immune_Checkpoint_54.rds"
)

output_dir <- file.path(
  project_dir,
  "Workflow-08-STEP8C8_Supplementary_Figure4"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


cat("\n============================================================\n")
cat("STEP 8C-8: Supplementary Figure 4\n")
cat("============================================================\n\n")


# ============================================================
# 2. Check locked RDS
# ============================================================

if (!file.exists(locked_rds)) {
  stop(
    "\nERROR: Locked STEP 8C-7B RDS file not found:\n",
    locked_rds,
    "\n"
  )
}


# ============================================================
# 3. Read LOCKED RDS
# ============================================================

supp_rds <- readRDS(locked_rds)


cat("Locked RDS class:\n")
print(class(supp_rds))


cat("\n============================================================\n")
cat("LOCKED RDS INTERNAL STRUCTURE\n")
cat("============================================================\n")


if (is.list(supp_rds)) {
  
  cat("\nObjects stored in the locked RDS:\n")
  print(names(supp_rds))
  
} else {
  
  cat("\nThe RDS itself is not a list.\n")
}


# ============================================================
# 4. Identify the locked 54-row supplementary table
#
# IMPORTANT:
#   We DO NOT calculate anything here.
#
#   We only locate the existing locked table inside the RDS.
# ============================================================

required_columns <- c(
  "Core_Gene",
  "Checkpoint_Gene",
  "N",
  "Unadjusted_Spearman_rho",
  "Unadjusted_P",
  "Unadjusted_FDR_BH",
  "Unadjusted_Direction",
  "Adjusted_Partial_Spearman_rho",
  "Adjusted_P",
  "Adjusted_FDR",
  "Adjusted_Direction",
  "Direction_Transition",
  "Direction_Switch",
  "Unadjusted_FDR_Significant",
  "Adjusted_FDR_Significant",
  "FDR_Transition"
)


find_locked_table <- function(x,
                              required_columns,
                              path = "RDS") {
  
  # ----------------------------------------------------------
  # If this object is a data.frame
  # ----------------------------------------------------------
  
  if (is.data.frame(x)) {
    
    current_columns <- colnames(x)
    
    if (
      nrow(x) == 54 &&
      all(required_columns %in% current_columns)
    ) {
      return(x)
    }
  }
  
  
  # ----------------------------------------------------------
  # If this object is a list, search recursively
  # ----------------------------------------------------------
  
  if (is.list(x)) {
    
    object_names <- names(x)
    
    if (is.null(object_names)) {
      object_names <- seq_along(x)
    }
    
    for (i in seq_along(x)) {
      
      current_name <- object_names[i]
      
      current_path <- paste0(
        path,
        "$",
        current_name
      )
      
      result <- find_locked_table(
        x[[i]],
        required_columns,
        current_path
      )
      
      if (!is.null(result)) {
        attr(result, "locked_path") <- current_path
        return(result)
      }
    }
  }
  
  return(NULL)
}


supp_table <- find_locked_table(
  supp_rds,
  required_columns
)


# ============================================================
# 5. Verify locked table was found
# ============================================================

if (is.null(supp_table)) {
  
  cat("\n============================================================\n")
  cat("ERROR: LOCKED 54-ROW TABLE NOT FOUND\n")
  cat("============================================================\n\n")
  
  stop(
    paste(
      "The RDS contains a list, but no 54-row data.frame",
      "with the expected locked columns could be identified."
    )
  )
}


cat("\n============================================================\n")
cat("LOCKED TABLE SUCCESSFULLY IDENTIFIED\n")
cat("============================================================\n\n")


cat(
  "Locked table dimensions: ",
  nrow(supp_table),
  " x ",
  ncol(supp_table),
  "\n",
  sep = ""
)


cat(
  "Locked table path inside RDS: ",
  attr(supp_table, "locked_path"),
  "\n",
  sep = ""
)


cat("\nColumn names:\n")
print(colnames(supp_table))


# ============================================================
# 6. Structural QC
# ============================================================

missing_columns <- setdiff(
  required_columns,
  colnames(supp_table)
)


if (length(missing_columns) > 0) {
  
  stop(
    paste0(
      "\nERROR: Required columns are missing:\n",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


if (nrow(supp_table) != 54) {
  
  stop(
    paste0(
      "\nERROR: Expected 54 associations; found ",
      nrow(supp_table),
      "."
    )
  )
}


if (length(unique(supp_table$Core_Gene)) != 3) {
  
  stop(
    "\nERROR: Expected exactly 3 core genes."
  )
}


if (length(unique(supp_table$Checkpoint_Gene)) != 18) {
  
  stop(
    "\nERROR: Expected exactly 18 checkpoint genes."
  )
}


cat("\nStructural QC: PASS\n")


# ============================================================
# 7. Lock gene orders
# ============================================================

core_order <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


checkpoint_order <- supp_table %>%
  distinct(Checkpoint_Gene) %>%
  pull(Checkpoint_Gene)


cat("\n============================================================\n")
cat("LOCKED CORE-GENE ORDER\n")
cat("============================================================\n")

print(core_order)


cat("\n============================================================\n")
cat("CHECKPOINT-GENE ORDER FROM LOCKED TABLE\n")
cat("============================================================\n")

print(checkpoint_order)


# ============================================================
# 8. Prepare plotting data
#
# NO statistics are calculated.
# ============================================================

plot_data <- supp_table %>%
  mutate(
    Core_Gene = factor(
      Core_Gene,
      levels = core_order
    ),
    Checkpoint_Gene = factor(
      Checkpoint_Gene,
      levels = rev(checkpoint_order)
    )
  )


# ============================================================
# 9. GLOBAL NUMERICAL QC
# ============================================================

n_total <- nrow(plot_data)


n_unadj_sig <- sum(
  plot_data$Unadjusted_FDR_BH < 0.05,
  na.rm = TRUE
)


n_adj_sig <- sum(
  plot_data$Adjusted_FDR < 0.05,
  na.rm = TRUE
)


n_both_sig <- sum(
  plot_data$Unadjusted_FDR_BH < 0.05 &
    plot_data$Adjusted_FDR < 0.05,
  na.rm = TRUE
)


n_lost <- sum(
  plot_data$Unadjusted_FDR_BH < 0.05 &
    plot_data$Adjusted_FDR >= 0.05,
  na.rm = TRUE
)


n_gained <- sum(
  plot_data$Unadjusted_FDR_BH >= 0.05 &
    plot_data$Adjusted_FDR < 0.05,
  na.rm = TRUE
)


n_direction_switch <- sum(
  plot_data$Direction_Switch %in% TRUE,
  na.rm = TRUE
)


n_pos_neg_adj_sig <- sum(
  plot_data$Direction_Transition ==
    "Positive -> Negative" &
    plot_data$Adjusted_FDR < 0.05,
  na.rm = TRUE
)


# ============================================================
# 10. Direction QC
# ============================================================

n_pp <- sum(
  plot_data$Direction_Transition ==
    "Positive -> Positive",
  na.rm = TRUE
)


n_nn <- sum(
  plot_data$Direction_Transition ==
    "Negative -> Negative",
  na.rm = TRUE
)


n_pn <- sum(
  plot_data$Direction_Transition ==
    "Positive -> Negative",
  na.rm = TRUE
)


n_np <- sum(
  plot_data$Direction_Transition ==
    "Negative -> Positive",
  na.rm = TRUE
)


# ============================================================
# 11. Raw direction QC
# ============================================================

unadj_positive <- sum(
  plot_data$Unadjusted_Spearman_rho > 0,
  na.rm = TRUE
)


unadj_negative <- sum(
  plot_data$Unadjusted_Spearman_rho < 0,
  na.rm = TRUE
)


adj_positive <- sum(
  plot_data$Adjusted_Partial_Spearman_rho > 0,
  na.rm = TRUE
)


adj_negative <- sum(
  plot_data$Adjusted_Partial_Spearman_rho < 0,
  na.rm = TRUE
)


# ============================================================
# 12. Global QC table
# ============================================================

qc_table <- data.frame(
  
  Item = c(
    "Total associations",
    "Unadjusted FDR significant",
    "Adjusted FDR significant",
    "Significant in both",
    "Significance lost",
    "Significance gained",
    "Direction switches",
    "Positive -> Positive",
    "Negative -> Negative",
    "Positive -> Negative",
    "Negative -> Positive",
    "Unadjusted positive",
    "Unadjusted negative",
    "Adjusted positive",
    "Adjusted negative",
    "Positive -> Negative & adjusted FDR < 0.05"
  ),
  
  Expected = c(
    54,
    48,
    38,
    34,
    14,
    4,
    32,
    21,
    1,
    32,
    0,
    53,
    1,
    21,
    33,
    21
  ),
  
  Observed = c(
    n_total,
    n_unadj_sig,
    n_adj_sig,
    n_both_sig,
    n_lost,
    n_gained,
    n_direction_switch,
    n_pp,
    n_nn,
    n_pn,
    n_np,
    unadj_positive,
    unadj_negative,
    adj_positive,
    adj_negative,
    n_pos_neg_adj_sig
  )
)


qc_table$PASS <- (
  qc_table$Expected ==
    qc_table$Observed
)


cat("\n============================================================\n")
cat("GLOBAL NUMERICAL QC\n")
cat("============================================================\n")

print(qc_table)


if (!all(qc_table$PASS)) {
  
  stop(
    "\nERROR: Global numerical QC FAILED.\n"
  )
}


cat("\nGLOBAL NUMERICAL QC: PASS\n")


# ============================================================
# 13. Core-gene switch QC
# ============================================================

core_switch_qc <- plot_data %>%
  
  group_by(Core_Gene) %>%
  
  summarise(
    N_pairs = n(),
    Direction_switches = sum(
      Direction_Switch %in% TRUE,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


expected_core_switches <- data.frame(
  
  Core_Gene = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  
  Expected_switches = c(
    13,
    12,
    7
  )
)


core_switch_qc_final <- expected_core_switches %>%
  
  left_join(
    core_switch_qc,
    by = "Core_Gene"
  ) %>%
  
  mutate(
    PASS =
      Expected_switches ==
      Direction_switches
  )


cat("\n============================================================\n")
cat("CORE-GENE DIRECTION-SWITCH QC\n")
cat("============================================================\n")

print(core_switch_qc_final)


if (!all(core_switch_qc_final$PASS)) {
  
  stop(
    "\nERROR: Core-gene direction-switch QC FAILED.\n"
  )
}


cat("\nCORE-GENE SWITCH QC: PASS\n")


# ============================================================
# 14. PANEL A
#
# Unadjusted Spearman rho
# versus
# ImmuneScore-adjusted partial Spearman rho
# ============================================================

panel_A_data <- plot_data %>%
  
  mutate(
    
    Transition = factor(
      Direction_Transition,
      
      levels = c(
        "Positive -> Positive",
        "Negative -> Negative",
        "Positive -> Negative",
        "Negative -> Positive"
      )
    ),
    
    Highlight = case_when(
      
      Direction_Transition ==
        "Positive -> Negative" &
        Adjusted_FDR < 0.05 ~
        
        "Positive -> Negative\nAdjusted FDR < 0.05",
      
      TRUE ~
        
        "Other"
    )
  )


panel_A <- ggplot(
  panel_A_data,
  
  aes(
    x = Unadjusted_Spearman_rho,
    y = Adjusted_Partial_Spearman_rho
  )
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  
  geom_point(
    aes(
      shape = Transition,
      fill = Highlight
    ),
    
    size = 3.2,
    stroke = 0.5
  ) +
  
  scale_shape_manual(
    
    values = c(
      "Positive -> Positive" = 21,
      "Negative -> Negative" = 22,
      "Positive -> Negative" = 24,
      "Negative -> Positive" = 25
    )
  ) +
  
  scale_fill_manual(
    
    values = c(
      "Positive -> Negative\nAdjusted FDR < 0.05" = "black",
      "Other" = "white"
    )
  ) +
  
  coord_cartesian(
    xlim = c(-1, 1),
    ylim = c(-1, 1)
  ) +
  
  scale_x_continuous(
    breaks = seq(
      -1,
      1,
      by = 0.5
    )
  ) +
  
  scale_y_continuous(
    breaks = seq(
      -1,
      1,
      by = 0.5
    )
  ) +
  
  labs(
    
    title = "A",
    
    x =
      "Unadjusted Spearman \u03c1",
    
    y =
      "ImmuneScore-adjusted partial Spearman \u03c1",
    
    shape =
      "Direction transition",
    
    fill = NULL
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold",
        size = 14,
        hjust = 0
      ),
    
    axis.title =
      element_text(
        size = 11
      ),
    
    axis.text =
      element_text(
        size = 9
      ),
    
    legend.title =
      element_text(
        size = 9
      ),
    
    legend.text =
      element_text(
        size = 8
      ),
    
    legend.position =
      "right"
  )


# ============================================================
# 15. PANEL B
#
# Heatmap of adjusted partial Spearman rho
# ============================================================

panel_B_data <- plot_data %>%
  
  mutate(
    
    Significance = case_when(
      
      Adjusted_FDR < 0.05 ~ "*",
      
      TRUE ~ ""
    )
  )


panel_B <- ggplot(
  
  panel_B_data,
  
  aes(
    x = Checkpoint_Gene,
    y = Core_Gene,
    fill = Adjusted_Partial_Spearman_rho
  )
  
) +
  
  geom_tile(
    
    color = "white",
    linewidth = 0.5
  ) +
  
  geom_text(
    
    aes(
      label = Significance
    ),
    
    size = 4,
    fontface = "bold"
  ) +
  
  scale_fill_gradient2(
    
    limits = c(
      -1,
      1
    ),
    
    midpoint = 0,
    
    name =
      "Adjusted\npartial Spearman \u03c1"
  ) +
  
  labs(
    
    title = "B",
    
    x = NULL,
    
    y = NULL
  ) +
  
  theme_minimal(
    
    base_size = 10
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold",
        size = 14,
        hjust = 0
      ),
    
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 8
      ),
    
    axis.text.y =
      element_text(
        face = "italic",
        size = 10
      ),
    
    panel.grid =
      element_blank(),
    
    legend.title =
      element_text(
        size = 8
      ),
    
    legend.text =
      element_text(
        size = 8
      )
  )


# ============================================================
# 16. PANEL C
#
# Direction-transition summary
# ============================================================

panel_C_data <- data.frame(
  
  Direction_Transition = factor(
    
    c(
      "Positive -> Positive",
      "Negative -> Negative",
      "Positive -> Negative",
      "Negative -> Positive"
    ),
    
    levels = c(
      "Positive -> Positive",
      "Negative -> Negative",
      "Positive -> Negative",
      "Negative -> Positive"
    )
  ),
  
  N = c(
    n_pp,
    n_nn,
    n_pn,
    n_np
  )
)


panel_C <- ggplot(
  
  panel_C_data,
  
  aes(
    x = Direction_Transition,
    y = N
  )
  
) +
  
  geom_col(
    
    width = 0.65
  ) +
  
  geom_text(
    
    aes(
      label = N
    ),
    
    vjust = -0.35,
    
    size = 4,
    
    fontface = "bold"
  ) +
  
  scale_y_continuous(
    
    limits = c(
      0,
      36
    ),
    
    breaks = seq(
      0,
      35,
      by = 5
    ),
    
    expand =
      expansion(
        mult = c(
          0,
          0.05
        )
      )
  ) +
  
  labs(
    
    title = "C",
    
    x = NULL,
    
    y =
      "Number of gene-pair associations"
  ) +
  
  theme_classic(
    
    base_size = 11
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold",
        size = 14,
        hjust = 0
      ),
    
    axis.text.x =
      element_text(
        angle = 25,
        hjust = 1,
        size = 9
      ),
    
    axis.title.y =
      element_text(
        size = 10
      )
  )


# ============================================================
# 17. DISPLAY PANELS
# ============================================================

cat("\n============================================================\n")
cat("Panel A generated.\n")
cat("Panel B generated.\n")
cat("Panel C generated.\n")
cat("============================================================\n")


print(panel_A)

print(panel_B)

print(panel_C)


# ============================================================
# 18. EXPORT PANEL A
#
# PDF = vector
# TIFF = 600 dpi
# ============================================================

panel_A_pdf <- file.path(
  output_dir,
  "STEP8C8_PanelA_Unadjusted_vs_Adjusted.pdf"
)

panel_A_tiff <- file.path(
  output_dir,
  "STEP8C8_PanelA_Unadjusted_vs_Adjusted_600dpi.tiff"
)


ggsave(
  
  filename =
    panel_A_pdf,
  
  plot =
    panel_A,
  
  width =
    7.2,
  
  height =
    6.2,
  
  units =
    "in",
  
  device =
    cairo_pdf,
  
  bg =
    "white"
)


ggsave(
  
  filename =
    panel_A_tiff,
  
  plot =
    panel_A,
  
  width =
    7.2,
  
  height =
    6.2,
  
  units =
    "in",
  
  dpi =
    600,
  
  compression =
    "lzw",
  
  bg =
    "white"
)
dev.off()

# ============================================================
# 19. EXPORT PANEL B
# ============================================================

panel_B_pdf <- file.path(
  output_dir,
  "STEP8C8_PanelB_Adjusted_Heatmap.pdf"
)

panel_B_tiff <- file.path(
  output_dir,
  "STEP8C8_PanelB_Adjusted_Heatmap_600dpi.tiff"
)


ggsave(
  
  filename =
    panel_B_pdf,
  
  plot =
    panel_B,
  
  width =
    9.0,
  
  height =
    3.8,
  
  units =
    "in",
  
  device =
    cairo_pdf,
  
  bg =
    "white"
)


ggsave(
  
  filename =
    panel_B_tiff,
  
  plot =
    panel_B,
  
  width =
    9.0,
  
  height =
    3.8,
  
  units =
    "in",
  
  dpi =
    600,
  
  compression =
    "lzw",
  
  bg =
    "white"
)


# ============================================================
# 20. EXPORT PANEL C
# ============================================================

panel_C_pdf <- file.path(
  output_dir,
  "STEP8C8_PanelC_Direction_Transition.pdf"
)

panel_C_tiff <- file.path(
  output_dir,
  "STEP8C8_PanelC_Direction_Transition_600dpi.tiff"
)


ggsave(
  
  filename =
    panel_C_pdf,
  
  plot =
    panel_C,
  
  width =
    7.2,
  
  height =
    4.8,
  
  units =
    "in",
  
  device =
    cairo_pdf,
  
  bg =
    "white"
)


ggsave(
  
  filename =
    panel_C_tiff,
  
  plot =
    panel_C,
  
  width =
    7.2,
  
  height =
    4.8,
  
  units =
    "in",
  
  dpi =
    600,
  
  compression =
    "lzw",
  
  bg =
    "white"
)


# ============================================================
# 21. Save plotting data
# ============================================================

write.csv(
  
  panel_A_data,
  
  file.path(
    output_dir,
    "STEP8C8_PanelA_plot_data.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  panel_B_data,
  
  file.path(
    output_dir,
    "STEP8C8_PanelB_plot_data.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  panel_C_data,
  
  file.path(
    output_dir,
    "STEP8C8_PanelC_plot_data.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 22. Save QC files
# ============================================================

write.csv(
  
  qc_table,
  
  file.path(
    output_dir,
    "STEP8C8_Global_Numerical_QC.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  core_switch_qc_final,
  
  file.path(
    output_dir,
    "STEP8C8_Core_Gene_Switch_QC.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 23. Save plotting objects
#
# Useful for later modifications without recalculating statistics.
# ============================================================

saveRDS(
  
  list(
    
    locked_source =
      locked_rds,
    
    locked_table =
      supp_table,
    
    panel_A_data =
      panel_A_data,
    
    panel_B_data =
      panel_B_data,
    
    panel_C_data =
      panel_C_data,
    
    panel_A =
      panel_A,
    
    panel_B =
      panel_B,
    
    panel_C =
      panel_C,
    
    global_QC =
      qc_table,
    
    core_gene_QC =
      core_switch_qc_final
    
  ),
  
  file.path(
    output_dir,
    "STEP8C8_Separate_Panels_Plot_Objects.rds"
  )
)


# ============================================================
# 24. Final file check
# ============================================================

output_files <- c(
  
  panel_A_pdf,
  panel_A_tiff,
  
  panel_B_pdf,
  panel_B_tiff,
  
  panel_C_pdf,
  panel_C_tiff
  
)


file_check <- data.frame(
  
  File = basename(output_files),
  
  Exists =
    file.exists(output_files),
  
  Size_MB =
    round(
      file.size(output_files) /
        1024^2,
      3
    )
  
)


cat("\n============================================================\n")
cat("OUTPUT FILE QC\n")
cat("============================================================\n")

print(file_check)


if (!all(file_check$Exists)) {
  
  stop(
    "\nERROR: One or more figure files were not created."
  )
}


# ============================================================
# 25. FINAL STATUS
# ============================================================

cat("\n============================================================\n")
cat("STEP 8C-8 FINAL STATUS\n")
cat("============================================================\n\n")

cat("Locked source successfully identified: PASS\n")
cat("Structural QC: PASS\n")
cat("Global numerical QC: PASS\n")
cat("Core-gene switch QC: PASS\n")
cat("Panel A: GENERATED\n")
cat("Panel B: GENERATED\n")
cat("Panel C: GENERATED\n")
cat("PDF export: PASS\n")
cat("600-dpi TIFF export: PASS\n")

cat("\nNo new statistical analysis was performed.\n")
cat("No correlation, P value, or FDR was recalculated.\n")

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("\n============================================================\n")
cat("STEP 8C-8 COMPLETE\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 8C-8
# Supplementary Figure 4
#
# VISUALIZATION REVISION v2
#
# Panels A / B / C generated separately
#
# IMPORTANT:
#   - STEP 8C-7B statistical results are LOCKED.
#   - No Spearman correlation is recalculated.
#   - No partial Spearman correlation is recalculated.
#   - No P values are recalculated.
#   - No FDR is recalculated.
#   - Visualization-only revision.
#   - Final panel assembly will be performed in Illustrator.
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

required_packages <- c(
  "ggplot2",
  "dplyr",
  "cowplot"
)

missing_pkgs <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(missing_pkgs, collapse = ", "),
      "\nPlease install them before running this script."
    )
  )
}

library(ggplot2)
library(dplyr)
library(cowplot)


# ============================================================
# 1. Project paths
# ============================================================

project_dir <- "input_dir/Paper1_tumor-only WGCNA"

locked_rds <- file.path(
  project_dir,
  "Workflow-08-STEP8C7B_Supplementary_Table_Immune_Checkpoint_54.rds"
)

output_dir <- file.path(
  project_dir,
  "Workflow-08-STEP8C8_Supplementary_Figure4"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


cat("\n")
cat("============================================================\n")
cat("STEP 8C-8 — Supplementary Figure 4\n")
cat("VISUALIZATION REVISION v2\n")
cat("============================================================\n")


# ============================================================
# 2. Read locked RDS
# ============================================================

if (!file.exists(locked_rds)) {
  stop(
    paste0(
      "\nLocked RDS not found:\n",
      locked_rds
    )
  )
}

supp_rds <- readRDS(locked_rds)


# ============================================================
# 3. Locate locked 54-row table
# ============================================================

required_columns <- c(
  "Core_Gene",
  "Checkpoint_Gene",
  "N",
  "Unadjusted_Spearman_rho",
  "Unadjusted_P",
  "Unadjusted_FDR_BH",
  "Unadjusted_Direction",
  "Adjusted_Partial_Spearman_rho",
  "Adjusted_P",
  "Adjusted_FDR",
  "Adjusted_Direction",
  "Direction_Transition",
  "Direction_Switch",
  "Unadjusted_FDR_Significant",
  "Adjusted_FDR_Significant",
  "FDR_Transition"
)


find_locked_table <- function(
    x,
    required_columns,
    path = "RDS"
) {
  
  if (is.data.frame(x)) {
    
    current_columns <- colnames(x)
    
    if (
      nrow(x) == 54 &&
      all(required_columns %in% current_columns)
    ) {
      
      attr(
        x,
        "locked_path"
      ) <- path
      
      return(x)
    }
  }
  
  
  if (is.list(x)) {
    
    object_names <- names(x)
    
    if (is.null(object_names)) {
      object_names <- seq_along(x)
    }
    
    for (i in seq_along(x)) {
      
      current_name <- object_names[i]
      
      current_path <- paste0(
        path,
        "$",
        current_name
      )
      
      result <- find_locked_table(
        x[[i]],
        required_columns,
        current_path
      )
      
      if (!is.null(result)) {
        return(result)
      }
    }
  }
  
  return(NULL)
}


supp_table <- find_locked_table(
  supp_rds,
  required_columns
)


if (is.null(supp_table)) {
  
  stop(
    paste0(
      "\nERROR: The locked 54-row table could not be identified ",
      "inside the STEP 8C-7B RDS."
    )
  )
}


cat("\nLocked 54-row table: FOUND\n")

cat(
  "Locked table path: ",
  attr(supp_table, "locked_path"),
  "\n",
  sep = ""
)

cat(
  "Dimensions: ",
  nrow(supp_table),
  " x ",
  ncol(supp_table),
  "\n",
  sep = ""
)


# ============================================================
# 4. Structural QC
# ============================================================

if (nrow(supp_table) != 54) {
  stop("Expected 54 associations.")
}

if (
  length(unique(supp_table$Core_Gene)) != 3
) {
  stop("Expected exactly 3 core genes.")
}

if (
  length(unique(supp_table$Checkpoint_Gene)) != 18
) {
  stop("Expected exactly 18 checkpoint genes.")
}


cat("Structural QC: PASS\n")


# ============================================================
# 5. Locked gene order
# ============================================================

core_order <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


checkpoint_order <- supp_table %>%
  distinct(Checkpoint_Gene) %>%
  pull(Checkpoint_Gene)


cat("\nCore-gene order:\n")
print(core_order)


cat("\nCheckpoint-gene order:\n")
print(checkpoint_order)


# ============================================================
# 6. Plotting data
# ============================================================

plot_data <- supp_table %>%
  mutate(
    
    Core_Gene = factor(
      Core_Gene,
      levels = core_order
    ),
    
    Checkpoint_Gene = factor(
      Checkpoint_Gene,
      levels = rev(checkpoint_order)
    )
    
  )


# ============================================================
# 7. Numerical QC
#
# Visualization only.
# No statistical values are recalculated.
# ============================================================

n_total <- nrow(plot_data)


n_unadj_sig <- sum(
  plot_data$Unadjusted_FDR_BH < 0.05,
  na.rm = TRUE
)


n_adj_sig <- sum(
  plot_data$Adjusted_FDR < 0.05,
  na.rm = TRUE
)


n_both_sig <- sum(
  plot_data$Unadjusted_FDR_BH < 0.05 &
    plot_data$Adjusted_FDR < 0.05,
  na.rm = TRUE
)


n_lost <- sum(
  plot_data$Unadjusted_FDR_BH < 0.05 &
    plot_data$Adjusted_FDR >= 0.05,
  na.rm = TRUE
)


n_gained <- sum(
  plot_data$Unadjusted_FDR_BH >= 0.05 &
    plot_data$Adjusted_FDR < 0.05,
  na.rm = TRUE
)


n_direction_switch <- sum(
  plot_data$Direction_Switch %in% TRUE,
  na.rm = TRUE
)


n_unadj_positive <- sum(
  plot_data$Unadjusted_Spearman_rho > 0,
  na.rm = TRUE
)


n_unadj_negative <- sum(
  plot_data$Unadjusted_Spearman_rho < 0,
  na.rm = TRUE
)


n_adj_positive <- sum(
  plot_data$Adjusted_Partial_Spearman_rho > 0,
  na.rm = TRUE
)


n_adj_negative <- sum(
  plot_data$Adjusted_Partial_Spearman_rho < 0,
  na.rm = TRUE
)


n_pn_adj_sig <- sum(
  plot_data$Direction_Transition ==
    "Positive -> Negative" &
    plot_data$Adjusted_FDR < 0.05,
  na.rm = TRUE
)


# Direction categories

n_pp <- sum(
  plot_data$Direction_Transition ==
    "Positive -> Positive",
  na.rm = TRUE
)


n_nn <- sum(
  plot_data$Direction_Transition ==
    "Negative -> Negative",
  na.rm = TRUE
)


n_pn <- sum(
  plot_data$Direction_Transition ==
    "Positive -> Negative",
  na.rm = TRUE
)


n_np <- sum(
  plot_data$Direction_Transition ==
    "Negative -> Positive",
  na.rm = TRUE
)


# ============================================================
# 8. Global QC
# ============================================================

qc_table <- data.frame(
  
  Item = c(
    "Total associations",
    "Unadjusted FDR significant",
    "Adjusted FDR significant",
    "Significant in both",
    "Significance lost",
    "Significance gained",
    "Direction switches",
    "Unadjusted positive",
    "Unadjusted negative",
    "Adjusted positive",
    "Adjusted negative",
    "Positive -> Positive",
    "Negative -> Negative",
    "Positive -> Negative",
    "Negative -> Positive",
    "Positive -> Negative & adjusted FDR < 0.05"
  ),
  
  Expected = c(
    54,
    48,
    38,
    34,
    14,
    4,
    32,
    53,
    1,
    21,
    33,
    21,
    1,
    32,
    0,
    21
  ),
  
  Observed = c(
    n_total,
    n_unadj_sig,
    n_adj_sig,
    n_both_sig,
    n_lost,
    n_gained,
    n_direction_switch,
    n_unadj_positive,
    n_unadj_negative,
    n_adj_positive,
    n_adj_negative,
    n_pp,
    n_nn,
    n_pn,
    n_np,
    n_pn_adj_sig
  )
)


qc_table$PASS <-
  qc_table$Expected ==
  qc_table$Observed


cat("\n============================================================\n")
cat("GLOBAL NUMERICAL QC\n")
cat("============================================================\n")

print(qc_table)


if (!all(qc_table$PASS)) {
  stop(
    "\nGLOBAL NUMERICAL QC FAILED.\n"
  )
}


cat("\nGLOBAL NUMERICAL QC: PASS\n")


# ============================================================
# 9. Core-gene switch QC
# ============================================================

core_switch_qc <- plot_data %>%
  
  group_by(Core_Gene) %>%
  
  summarise(
    
    N_pairs = n(),
    
    Direction_switches =
      sum(
        Direction_Switch %in% TRUE,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


expected_core_switches <- data.frame(
  
  Core_Gene = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  
  Expected_switches = c(
    13,
    12,
    7
  )
)


core_switch_qc_final <-
  expected_core_switches %>%
  
  left_join(
    core_switch_qc,
    by = "Core_Gene"
  ) %>%
  
  mutate(
    PASS =
      Expected_switches ==
      Direction_switches
  )


print(core_switch_qc_final)


if (!all(core_switch_qc_final$PASS)) {
  stop(
    "\nCORE-GENE SWITCH QC FAILED.\n"
  )
}


cat("\nCORE-GENE SWITCH QC: PASS\n")


# ============================================================
# 10. Locked visualization palettes
# ============================================================


# ------------------------------------------------------------
# Panel A
# Direction-transition colors
# ------------------------------------------------------------

direction_colors <- c(
  
  "Positive -> Positive" =
    "#3B82F6",
  
  "Negative -> Negative" =
    "#22A06B",
  
  "Positive -> Negative" =
    "#8B5CF6"
)


direction_shapes <- c(
  
  "Positive -> Positive" =
    21,
  
  "Negative -> Negative" =
    22,
  
  "Positive -> Negative" =
    24
)


# ------------------------------------------------------------
# Panel B
# rho colors
# ------------------------------------------------------------

rho_low_color <- "#82B5EC"

rho_mid_color <- "white"

rho_high_color <- "#F57C6E"


# ------------------------------------------------------------
# Panel C
# Direction colors
# ------------------------------------------------------------

transition_colors_C <- c(
  
  "Positive -> Positive" =
    "#3B82F6",
  
  "Negative -> Negative" =
    "#22A06B",
  
  "Positive -> Negative" =
    "#8B5CF6",
  
  "Negative -> Positive" =
    "#9CA3AF"
)


# ============================================================
# 11. PANEL A DATA
# ============================================================

panel_A_data <- plot_data %>%
  
  mutate(
    
    Transition = case_when(
      
      Direction_Transition ==
        "Positive -> Positive" ~
        "Positive -> Positive",
      
      Direction_Transition ==
        "Negative -> Negative" ~
        "Negative -> Negative",
      
      Direction_Transition ==
        "Positive -> Negative" ~
        "Positive -> Negative",
      
      TRUE ~
        NA_character_
    ),
    
    Adjusted_Significant_Switch =
      Direction_Transition ==
      "Positive -> Negative" &
      Adjusted_FDR < 0.05
    
  )


# ============================================================
# 12. PANEL A
#
# x = unadjusted Spearman rho
# y = adjusted partial Spearman rho
#
# Color + shape = direction transition
#
# Black outer ring = positive-to-negative switch
#                       remaining significant after adjustment
# ============================================================

panel_A_base <- panel_A_data %>%
  
  filter(
    !is.na(Transition)
  )


panel_A <- ggplot(
  
  panel_A_base,
  
  aes(
    
    x =
      Unadjusted_Spearman_rho,
    
    y =
      Adjusted_Partial_Spearman_rho,
    
    shape =
      Transition,
    
    fill =
      Transition
  )
  
) +
  
  geom_hline(
    
    yintercept = 0,
    
    linewidth = 0.35,
    
    linetype = "dashed",
    
    color = "grey55"
  ) +
  
  geom_vline(
    
    xintercept = 0,
    
    linewidth = 0.35,
    
    linetype = "dashed",
    
    color = "grey55"
  ) +
  
  geom_point(
    
    size = 3.1,
    
    stroke = 0.35,
    
    color = "white"
  ) +
  
  # ----------------------------------------------------------
# Black outer ring for the 21 adjusted-significant
# Positive -> Negative associations
# ----------------------------------------------------------

geom_point(
  
  data =
    panel_A_base %>%
    filter(
      Adjusted_Significant_Switch
    ),
  
  aes(
    
    x =
      Unadjusted_Spearman_rho,
    
    y =
      Adjusted_Partial_Spearman_rho
  ),
  
  shape = 24,
  
  fill = NA,
  
  color = "black",
  
  stroke = 0.9,
  
  size = 4.2,
  
  inherit.aes = FALSE
) +
  
  scale_shape_manual(
    
    values =
      direction_shapes,
    
    breaks = c(
      
      "Positive -> Positive",
      
      "Negative -> Negative",
      
      "Positive -> Negative"
    ),
    
    labels = c(
      
      "Positive \u2192 Positive",
      
      "Negative \u2192 Negative",
      
      "Positive \u2192 Negative"
    )
  ) +
  
  scale_fill_manual(
    
    values =
      direction_colors,
    
    breaks = c(
      
      "Positive -> Positive",
      
      "Negative -> Negative",
      
      "Positive -> Negative"
    ),
    
    labels = c(
      
      "Positive \u2192 Positive",
      
      "Negative \u2192 Negative",
      
      "Positive \u2192 Negative"
    )
  ) +
  
  coord_cartesian(
    
    xlim = c(
      -1,
      1
    ),
    
    ylim = c(
      -1,
      1
    )
  ) +
  
  scale_x_continuous(
    
    breaks = seq(
      -1,
      1,
      by = 0.5
    )
  ) +
  
  scale_y_continuous(
    
    breaks = seq(
      -1,
      1,
      by = 0.5
    )
  ) +
  
  labs(
    
    x =
      "Unadjusted Spearman \u03c1",
    
    y =
      "ImmuneScore-adjusted partial Spearman \u03c1",
    
    shape =
      "Direction transition",
    
    fill =
      "Direction transition"
  ) +
  
  theme_classic(
    
    base_family =
      "Arial",
    
    base_size =
      8
  ) +
  
  theme(
    
    axis.title = element_text(
      
      family =
        "Arial",
      
      face =
        "bold",
      
      size =
        10,
      
      color =
        "black"
    ),
    
    axis.text = element_text(
      
      family =
        "Arial",
      
      size =
        8,
      
      color =
        "black"
    ),
    
    axis.ticks =
      element_line(
        linewidth = 0.3
      ),
    
    legend.title =
      element_text(
        
        family =
          "Arial",
        
        face =
          "bold",
        
        size =
          8
      ),
    
    legend.text =
      element_text(
        
        family =
          "Arial",
        
        size =
          8
      ),
    
    legend.key.size =
      unit(
        0.55,
        "cm"
      ),
    
    legend.position =
      "right",
    
    plot.margin =
      margin(
        5,
        5,
        5,
        5
      )
  )


# ============================================================
# 13. PANEL B DATA
# ============================================================

panel_B_data <- plot_data %>%
  
  mutate(
    
    Sig = case_when(
      
      Adjusted_FDR < 0.001 ~
        "***",
      
      Adjusted_FDR < 0.01 ~
        "**",
      
      Adjusted_FDR < 0.05 ~
        "*",
      
      TRUE ~
        ""
    )
  )


# ============================================================
# 14. PANEL B — MAIN HEATMAP
# ============================================================

p_B_heatmap <- ggplot(
  
  panel_B_data,
  
  aes(
    
    x =
      Checkpoint_Gene,
    
    y =
      Core_Gene,
    
    fill =
      Adjusted_Partial_Spearman_rho
  )
  
) +
  
  geom_tile(
    
    color =
      "white",
    
    linewidth =
      0.30
  ) +
  
  geom_text(
    
    aes(
      label =
        Sig
    ),
    
    family =
      "Arial",
    
    size =
      3.0,
    
    color =
      "black"
  ) +
  
  scale_fill_gradient2(
    
    low =
      rho_low_color,
    
    mid =
      rho_mid_color,
    
    high =
      rho_high_color,
    
    midpoint =
      0,
    
    limits =
      c(
        -1,
        1
      ),
    
    guide =
      "none"
  ) +
  
  scale_x_discrete(
    drop = FALSE
  ) +
  
  scale_y_discrete(
    drop = FALSE
  ) +
  
  labs(
    
    x =
      NULL,
    
    y =
      NULL
  ) +
  
  theme_minimal(
    
    base_family =
      "Arial",
    
    base_size =
      8
  ) +
  
  theme(
    
    panel.grid =
      element_blank(),
    
    axis.title =
      element_blank(),
    
    axis.text.x =
      element_text(
        
        family =
          "Arial",
        
        size =
          8,
        
        color =
          "black",
        
        angle =
          45,
        
        hjust =
          1,
        
        vjust =
          1
      ),
    
    axis.text.y =
      element_text(
        
        family =
          "Arial",
        
        size =
          8,
        
        color =
          "black"
      ),
    
    axis.ticks =
      element_blank(),
    
    legend.position =
      "none",
    
    plot.margin =
      margin(
        4,
        4,
        4,
        4
      )
  )


# ============================================================
# 15. PANEL B — FDR LEGEND
# ============================================================

p_B_fdr <- ggplot() +
  
  annotate(
    
    "text",
    
    x =
      0,
    
    y =
      1.00,
    
    label =
      "FDR",
    
    hjust =
      0,
    
    vjust =
      1,
    
    family =
      "Arial",
    
    fontface =
      "bold",
    
    size =
      3.1
  ) +
  
  annotate(
    
    "text",
    
    x =
      0,
    
    y =
      0.70,
    
    label =
      "***  FDR < 0.001",
    
    hjust =
      0,
    
    vjust =
      1,
    
    family =
      "Arial",
    
    size =
      2.7
  ) +
  
  annotate(
    
    "text",
    
    x =
      0,
    
    y =
      0.43,
    
    label =
      "**   FDR < 0.01",
    
    hjust =
      0,
    
    vjust =
      1,
    
    family =
      "Arial",
    
    size =
      2.7
  ) +
  
  annotate(
    
    "text",
    
    x =
      0,
    
    y =
      0.16,
    
    label =
      "*    FDR < 0.05",
    
    hjust =
      0,
    
    vjust =
      1,
    
    family =
      "Arial",
    
    size =
      2.7
  ) +
  
  coord_cartesian(
    
    xlim =
      c(
        0,
        1
      ),
    
    ylim =
      c(
        0,
        1
      ),
    
    clip =
      "off"
  ) +
  
  theme_void(
    
    base_family =
      "Arial"
  ) +
  
  theme(
    
    plot.margin =
      margin(
        0,
        0,
        0,
        0
      )
  )


# ============================================================
# 16. PANEL B — rho COLORBAR
#
# Same blue-white-red visual concept as WGCNA heatmap
# ============================================================

rho_edges <- seq(
  
  -1,
  
  1,
  
  length.out =
    102
)


rho_bar_df <- data.frame(
  
  xmin =
    0,
  
  xmax =
    1,
  
  ymin =
    rho_edges[
      -length(rho_edges)
    ],
  
  ymax =
    rho_edges[
      -1
    ],
  
  rho =
    seq(
      -1,
      1,
      length.out =
        101
    )
)


p_B_rho <- ggplot(
  
  rho_bar_df
  
) +
  
  geom_rect(
    
    aes(
      
      xmin =
        xmin,
      
      xmax =
        xmax,
      
      ymin =
        ymin,
      
      ymax =
        ymax,
      
      fill =
        rho
    ),
    
    color =
      NA
  ) +
  
  scale_fill_gradient2(
    
    low =
      rho_low_color,
    
    mid =
      rho_mid_color,
    
    high =
      rho_high_color,
    
    midpoint =
      0,
    
    limits =
      c(
        -1,
        1
      ),
    
    guide =
      "none"
  ) +
  
  scale_x_continuous(
    
    limits =
      c(
        0,
        1
      ),
    
    expand =
      c(
        0,
        0
      )
  ) +
  
  scale_y_continuous(
    
    limits =
      c(
        -1,
        1
      ),
    
    breaks =
      c(
        -1,
        -0.5,
        0,
        0.5,
        1
      ),
    
    expand =
      c(
        0,
        0
      )
  ) +
  
  theme_minimal(
    
    base_family =
      "Arial"
  ) +
  
  theme(
    
    panel.grid =
      element_blank(),
    
    axis.title =
      element_blank(),
    
    axis.text.x =
      element_blank(),
    
    axis.ticks.x =
      element_blank(),
    
    axis.text.y =
      element_text(
        
        family =
          "Arial",
        
        size =
          8,
        
        color =
          "black"
      ),
    
    axis.ticks.y =
      element_line(
        
        linewidth =
          0.30,
        
        color =
          "black"
      ),
    
    legend.position =
      "none",
    
    plot.margin =
      margin(
        4,
        8,
        10,
        2
      )
  )


# ============================================================
# 17. PANEL B — rho title
# ============================================================

p_B_rho_title <- ggplot() +
  
  annotate(
    
    "text",
    
    x =
      0,
    
    y =
      0,
    
    label =
      "\u03c1",
    
    hjust =
      0,
    
    vjust =
      0,
    
    family =
      "Arial",
    
    fontface =
      "bold",
    
    size =
      3.1
  ) +
  
  coord_cartesian(
    
    xlim =
      c(
        0,
        1
      ),
    
    ylim =
      c(
        0,
        1
      ),
    
    clip =
      "off"
  ) +
  
  theme_void(
    
    base_family =
      "Arial"
  )


# ============================================================
# 18. PANEL B — rho section
# ============================================================

p_B_rho_section <- cowplot::plot_grid(
  
  p_B_rho_title,
  
  p_B_rho,
  
  ncol =
    1,
  
  align =
    "v",
  
  rel_heights =
    c(
      0.16,
      1
    )
)


# ============================================================
# 19. PANEL B — right-side legend
# ============================================================

p_B_right <- cowplot::plot_grid(
  
  p_B_fdr,
  
  p_B_rho_section,
  
  ncol =
    1,
  
  align =
    "v",
  
  rel_heights =
    c(
      0.62,
      1.55
    )
)


# ============================================================
# 20. PANEL B — final composition
#
# Panel B itself is one independent figure.
# It contains heatmap + right-side legend.
# ============================================================

panel_B <- cowplot::plot_grid(
  
  p_B_heatmap,
  
  p_B_right,
  
  ncol =
    2,
  
  align =
    "h",
  
  rel_widths =
    c(
      1,
      0.18
    )
)


# ============================================================
# 21. PANEL C DATA
# ============================================================

panel_C_data <- data.frame(
  
  Direction_Transition = factor(
    
    c(
      
      "Positive -> Positive",
      
      "Negative -> Negative",
      
      "Positive -> Negative",
      
      "Negative -> Positive"
      
    ),
    
    levels = c(
      
      "Positive -> Positive",
      
      "Negative -> Negative",
      
      "Positive -> Negative",
      
      "Negative -> Positive"
      
    )
  ),
  
  N = c(
    
    n_pp,
    
    n_nn,
    
    n_pn,
    
    n_np
    
  )
)


# ============================================================
# 22. PANEL C
# ============================================================

panel_C <- ggplot(
  
  panel_C_data,
  
  aes(
    
    x =
      Direction_Transition,
    
    y =
      N,
    
    fill =
      Direction_Transition
  )
  
) +
  
  geom_col(
    
    width =
      0.65
  ) +
  
  geom_text(
    
    aes(
      label =
        N
    ),
    
    vjust =
      -0.35,
    
    family =
      "Arial",
    
    size =
      3.2,
    
    fontface =
      "bold"
  ) +
  
  scale_fill_manual(
    
    values =
      transition_colors_C,
    
    breaks = c(
      
      "Positive -> Positive",
      
      "Negative -> Negative",
      
      "Positive -> Negative",
      
      "Negative -> Positive"
    ),
    
    labels = c(
      
      "Positive \u2192 Positive",
      
      "Negative \u2192 Negative",
      
      "Positive \u2192 Negative",
      
      "Negative \u2192 Positive"
    ),
    
    drop =
      FALSE
  ) +
  
  scale_y_continuous(
    
    limits =
      c(
        0,
        36
      ),
    
    breaks =
      seq(
        0,
        35,
        by =
          5
      ),
    
    expand =
      expansion(
        mult =
          c(
            0,
            0.05
          )
      )
  ) +
  
  labs(
    
    x =
      NULL,
    
    y =
      "Number of gene-pair associations",
    
    fill =
      "Direction transition"
  ) +
  
  theme_classic(
    
    base_family =
      "Arial",
    
    base_size =
      8
  ) +
  
  theme(
    
    axis.title.y =
      element_text(
        
        family =
          "Arial",
        
        face =
          "bold",
        
        size =
          10,
        
        color =
          "black"
      ),
    
    axis.text.x =
      element_text(
        
        family =
          "Arial",
        
        size =
          8,
        
        color =
          "black",
        
        angle =
          25,
        
        hjust =
          1
      ),
    
    axis.text.y =
      element_text(
        
        family =
          "Arial",
        
        size =
          8,
        
        color =
          "black"
      ),
    
    legend.title =
      element_text(
        
        family =
          "Arial",
        
        face =
          "bold",
        
        size =
          8
      ),
    
    legend.text =
      element_text(
        
        family =
          "Arial",
        
        size =
          8
      ),
    
    legend.position =
      "right"
  )


# ============================================================
# 23. DISPLAY
#
# No A/B/C labels are added.
# ============================================================

cat("\n============================================================\n")
cat("VISUALIZATION GENERATED\n")
cat("============================================================\n")

print(panel_A)

print(panel_B)

print(panel_C)


# ============================================================
# 24. PDF export function
#
# IMPORTANT:
#   Explicit device opening and dev.off()
#   are used intentionally.
#
#   This prevents unfinished PDF devices from remaining open.
# ============================================================

save_pdf_cairo <- function(
    plot_object,
    filename,
    width,
    height
) {
  
  cairo_pdf(
    
    filename =
      filename,
    
    width =
      width,
    
    height =
      height,
    
    family =
      "Arial"
  )
  
  print(plot_object)
  
  dev.off()
}


# ============================================================
# 25. TIFF export function
#
# Explicit dev.off() is also used.
# ============================================================

save_tiff_cairo <- function(
    plot_object,
    filename,
    width,
    height,
    dpi = 600
) {
  
  tiff(
    
    filename =
      filename,
    
    width =
      width,
    
    height =
      height,
    
    units =
      "in",
    
    res =
      dpi,
    
    compression =
      "lzw",
    
    type =
      "cairo"
  )
  
  print(plot_object)
  
  dev.off()
}


# ============================================================
# 26. Output filenames
# ============================================================

panel_A_pdf <- file.path(
  
  output_dir,
  
  "STEP8C8_PanelA_Unadjusted_vs_Adjusted_FINAL.pdf"
)


panel_A_tiff <- file.path(
  
  output_dir,
  
  "STEP8C8_PanelA_Unadjusted_vs_Adjusted_FINAL_600dpi.tiff"
)


panel_B_pdf <- file.path(
  
  output_dir,
  
  "STEP8C8_PanelB_Adjusted_Heatmap_FINAL.pdf"
)


panel_B_tiff <- file.path(
  
  output_dir,
  
  "STEP8C8_PanelB_Adjusted_Heatmap_FINAL_600dpi.tiff"
)


panel_C_pdf <- file.path(
  
  output_dir,
  
  "STEP8C8_PanelC_Direction_Transition_FINAL.pdf"
)


panel_C_tiff <- file.path(
  
  output_dir,
  
  "STEP8C8_PanelC_Direction_Transition_FINAL_600dpi.tiff"
)


# ============================================================
# 27. Export Panel A
# ============================================================

save_pdf_cairo(
  
  plot_object =
    panel_A,
  
  filename =
    panel_A_pdf,
  
  width =
    7.2,
  
  height =
    6.2
)


save_tiff_cairo(
  
  plot_object =
    panel_A,
  
  filename =
    panel_A_tiff,
  
  width =
    7.2,
  
  height =
    6.2,
  
  dpi =
    600
)


cat("\nPanel A PDF/TIFF export: PASS\n")


# ============================================================
# 28. Export Panel B
# ============================================================

save_pdf_cairo(
  
  plot_object =
    panel_B,
  
  filename =
    panel_B_pdf,
  
  width =
    9.0,
  
  height =
    4.2
)


save_tiff_cairo(
  
  plot_object =
    panel_B,
  
  filename =
    panel_B_tiff,
  
  width =
    9.0,
  
  height =
    4.2,
  
  dpi =
    600
)


cat("Panel B PDF/TIFF export: PASS\n")


# ============================================================
# 29. Export Panel C
# ============================================================

save_pdf_cairo(
  
  plot_object =
    panel_C,
  
  filename =
    panel_C_pdf,
  
  width =
    7.2,
  
  height =
    4.8
)


save_tiff_cairo(
  
  plot_object =
    panel_C,
  
  filename =
    panel_C_tiff,
  
  width =
    7.2,
  
  height =
    4.8,
  
  dpi =
    600
)


cat("Panel C PDF/TIFF export: PASS\n")


# ============================================================
# 30. Export plotting data for audit
# ============================================================

write.csv(
  
  panel_A_data,
  
  file.path(
    output_dir,
    "STEP8C8_PanelA_FINAL_plot_data.csv"
  ),
  
  row.names =
    FALSE,
  
  quote =
    TRUE
)


write.csv(
  
  panel_B_data,
  
  file.path(
    output_dir,
    "STEP8C8_PanelB_FINAL_plot_data.csv"
  ),
  
  row.names =
    FALSE,
  
  quote =
    TRUE
)


write.csv(
  
  panel_C_data,
  
  file.path(
    output_dir,
    "STEP8C8_PanelC_FINAL_plot_data.csv"
  ),
  
  row.names =
    FALSE,
  
  quote =
    TRUE
)


# ============================================================
# 31. Save QC
# ============================================================

write.csv(
  
  qc_table,
  
  file.path(
    output_dir,
    "STEP8C8_FINAL_Global_QC.csv"
  ),
  
  row.names =
    FALSE,
  
  quote =
    TRUE
)


write.csv(
  
  core_switch_qc_final,
  
  file.path(
    output_dir,
    "STEP8C8_FINAL_Core_Gene_Switch_QC.csv"
  ),
  
  row.names =
    FALSE,
  
  quote =
    TRUE
)


# ============================================================
# 32. Output-file QC
# ============================================================

output_files <- c(
  
  panel_A_pdf,
  
  panel_A_tiff,
  
  panel_B_pdf,
  
  panel_B_tiff,
  
  panel_C_pdf,
  
  panel_C_tiff
  
)


output_check <- data.frame(
  
  File =
    basename(output_files),
  
  Exists =
    file.exists(output_files),
  
  Size_MB =
    round(
      file.size(output_files) /
        1024^2,
      3
    )
)


cat("\n============================================================\n")
cat("OUTPUT FILE QC\n")
cat("============================================================\n")

print(output_check)


if (!all(output_check$Exists)) {
  
  stop(
    "\nERROR: One or more output files were not created."
  )
}


# ============================================================
# 33. FINAL STATUS
# ============================================================

cat("\n============================================================\n")
cat("STEP 8C-8 VISUALIZATION REVISION FINAL STATUS\n")
cat("============================================================\n\n")

cat("Locked source: PASS\n")
cat("Structural QC: PASS\n")
cat("Global numerical QC: PASS\n")
cat("Core-gene switch QC: PASS\n")
cat("Panel A generated: PASS\n")
cat("Panel B generated: PASS\n")
cat("Panel C generated: PASS\n")
cat("Panel A PDF: PASS\n")
cat("Panel A 600-dpi TIFF: PASS\n")
cat("Panel B PDF: PASS\n")
cat("Panel B 600-dpi TIFF: PASS\n")
cat("Panel C PDF: PASS\n")
cat("Panel C 600-dpi TIFF: PASS\n")

cat("\nNo statistical recalculation was performed.\n")
cat("No correlation, P value, or FDR was recalculated.\n")
cat("No A/B/C panel labels were added.\n")
cat("Final panel assembly should be performed in Illustrator.\n")

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("\n============================================================\n")
cat("STEP 8C-8 VISUALIZATION REVISION COMPLETE\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 8C-8
# Supplementary Figure 4 — Panel A
#
# Unadjusted vs ImmuneScore-adjusted Spearman correlation
#
# LOCKED DATA SOURCE:
# Workflow-08-STEP8C7B_Supplementary_Table_Immune_Checkpoint_54.rds
#
# IMPORTANT:
# 1. No correlation is recalculated.
# 2. No P value is recalculated.
# 3. No FDR is recalculated.
# 4. All values come directly from the locked STEP 8C-7B RDS.
# 5. Panel A is exported separately as PDF and 600-dpi TIFF.
# 6. No panel label "A" is added.
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

library(ggplot2)
library(dplyr)


# ============================================================
# 1. Project paths
# ============================================================

project_dir <- "input_dir/Paper1_tumor-only WGCNA"

rds_file <- file.path(
  project_dir,
  "Workflow-08-STEP8C7B_Supplementary_Table_Immune_Checkpoint_54.rds"
)

output_dir <- file.path(
  project_dir,
  "Workflow-08-STEP8C8_Supplementary_Figure4"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ============================================================
# 2. Output files
# ============================================================

panelA_pdf <- file.path(
  output_dir,
  "Supplementary_Figure4_PanelA_ImmuneCheckpoint_Correlation.pdf"
)

panelA_tiff <- file.path(
  output_dir,
  "Supplementary_Figure4_PanelA_ImmuneCheckpoint_Correlation_600dpi.tiff"
)


# ============================================================
# 3. Read LOCKED STEP 8C-7B RDS
# ============================================================

cat("\n============================================================\n")
cat("STEP 8C-8: Panel A\n")
cat("============================================================\n")

cat("\nReading locked RDS:\n")
cat(rds_file, "\n")

supp_table <- readRDS(rds_file)


# ============================================================
# 4. Identify the actual locked data frame inside the RDS
# ============================================================

required_columns <- c(
  "Core_Gene",
  "Checkpoint_Gene",
  "N",
  "Unadjusted_Spearman_rho",
  "Unadjusted_P",
  "Unadjusted_FDR_BH",
  "Unadjusted_Direction",
  "Adjusted_Partial_Spearman_rho",
  "Adjusted_P",
  "Adjusted_FDR",
  "Adjusted_Direction",
  "Direction_Transition",
  "Direction_Switch",
  "Unadjusted_FDR_Significant",
  "Adjusted_FDR_Significant",
  "FDR_Transition"
)


if (is.data.frame(supp_table)) {
  
  complete_54 <- supp_table
  
} else if (is.list(supp_table)) {
  
  candidate_idx <- which(
    vapply(
      supp_table,
      function(x) {
        is.data.frame(x) &&
          all(required_columns %in% names(x))
      },
      logical(1)
    )
  )
  
  if (length(candidate_idx) != 1) {
    
    stop(
      paste0(
        "\nERROR: Could not uniquely identify the locked ",
        "54-association data frame inside the RDS.\n",
        "Number of candidate data frames found: ",
        length(candidate_idx)
      )
    )
  }
  
  complete_54 <- supp_table[[candidate_idx]]
  
} else {
  
  stop(
    "\nERROR: Locked RDS is neither a data frame nor a suitable list."
  )
}


# ============================================================
# 5. Structural QC
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("Structural QC\n")
cat("------------------------------------------------------------\n")

cat(
  "Rows:",
  nrow(complete_54),
  "\n"
)

cat(
  "Columns:",
  ncol(complete_54),
  "\n"
)

missing_columns <- setdiff(
  required_columns,
  names(complete_54)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste0(
      "\nERROR: Required columns missing:\n",
      paste(missing_columns, collapse = ", ")
    )
  )
}

if (nrow(complete_54) != 54) {
  
  stop(
    paste0(
      "\nERROR: Expected 54 associations, found ",
      nrow(complete_54),
      "."
    )
  )
  
}

cat("Required columns: PASS\n")
cat("54 associations: PASS\n")


# ============================================================
# 6. Prepare Panel A plotting data
#
# NO statistical calculation is performed.
# ============================================================

panelA_df <- complete_54 %>%
  mutate(
    
    # --------------------------------------------------------
    # Display labels for direction transition
    # --------------------------------------------------------
    
    Transition_Display = case_when(
      
      Direction_Transition == "Positive -> Positive" ~
        "Positive → Positive",
      
      Direction_Transition == "Negative -> Negative" ~
        "Negative → Negative",
      
      Direction_Transition == "Positive -> Negative" ~
        "Positive → Negative",
      
      Direction_Transition == "Negative -> Positive" ~
        "Negative → Positive",
      
      TRUE ~ as.character(Direction_Transition)
    ),
    
    # --------------------------------------------------------
    # Locked adjusted-FDR status
    #
    # We use the existing Adjusted_FDR value.
    # No FDR calculation is performed.
    # --------------------------------------------------------
    
    Adjusted_Sig =
      Adjusted_FDR < 0.05,
    
    # --------------------------------------------------------
    # The specific 21 associations requested:
    #
    # Positive -> Negative
    # AND adjusted FDR < 0.05
    # --------------------------------------------------------
    
    PN_Adjusted_Sig =
      Transition_Display == "Positive → Negative" &
      Adjusted_FDR < 0.05
  )


# ============================================================
# 7. Panel A numerical QC
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("Panel A numerical QC\n")
cat("------------------------------------------------------------\n")

cat(
  "Total associations:",
  nrow(panelA_df),
  "\n"
)

cat(
  "Positive → Positive:",
  sum(
    panelA_df$Transition_Display ==
      "Positive → Positive"
  ),
  "\n"
)

cat(
  "Negative → Negative:",
  sum(
    panelA_df$Transition_Display ==
      "Negative → Negative"
  ),
  "\n"
)

cat(
  "Positive → Negative:",
  sum(
    panelA_df$Transition_Display ==
      "Positive → Negative"
  ),
  "\n"
)

cat(
  "Negative → Positive:",
  sum(
    panelA_df$Transition_Display ==
      "Negative → Positive"
  ),
  "\n"
)

cat(
  "Positive → Negative AND adjusted FDR < 0.05:",
  sum(panelA_df$PN_Adjusted_Sig),
  "\n"
)


# ------------------------------------------------------------
# Locked expected values
# ------------------------------------------------------------

stopifnot(
  nrow(panelA_df) == 54
)

stopifnot(
  sum(
    panelA_df$Transition_Display ==
      "Positive → Positive"
  ) == 21
)

stopifnot(
  sum(
    panelA_df$Transition_Display ==
      "Negative → Negative"
  ) == 1
)

stopifnot(
  sum(
    panelA_df$Transition_Display ==
      "Positive → Negative"
  ) == 32
)

stopifnot(
  sum(
    panelA_df$Transition_Display ==
      "Negative → Positive"
  ) == 0
)

stopifnot(
  sum(panelA_df$PN_Adjusted_Sig) == 21
)

cat("\nPanel A numerical QC: PASS\n")


# ============================================================
# 8. Define direction-transition colors
# ============================================================

transition_colors <- c(
  
  "Positive → Positive" =
    "#3B82F6",
  
  "Negative → Negative" =
    "#22A06B",
  
  "Positive → Negative" =
    "#D95F59"
)


# ============================================================
# 9. Define direction-transition shapes
# ============================================================

transition_shapes <- c(
  
  "Positive → Positive" = 21,   # circle
  
  "Negative → Negative" = 22,   # square
  
  "Positive → Negative" = 24    # triangle
)


# ============================================================
# 10. Build Panel A
# ============================================================

pA <- ggplot(
  panelA_df,
  aes(
    x = Unadjusted_Spearman_rho,
    y = Adjusted_Partial_Spearman_rho
  )
) +
  
  # ----------------------------------------------------------
# Zero reference lines
# ----------------------------------------------------------

geom_hline(
  yintercept = 0,
  linetype = "dashed",
  linewidth = 0.35,
  color = "grey65"
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey65"
  ) +
  
  # ----------------------------------------------------------
# Main 54 points
#
# Shape + fill = direction transition
# ----------------------------------------------------------

geom_point(
  aes(
    shape = Transition_Display,
    fill = Transition_Display
  ),
  size = 3.1,
  stroke = 0.45,
  color = "white"
) +
  
  # ----------------------------------------------------------
# Black outer ring
#
# Only for:
# Positive → Negative
# AND adjusted FDR < 0.05
#
# This layer does NOT change the original point.
# It only emphasizes the 21 significant associations.
# ----------------------------------------------------------
# ============================================================
# Significant Positive → Negative associations
# ============================================================

# 1. Black outer triangle
geom_point(
  data = panelA_df %>%
    filter(PN_Adjusted_Sig),
  
  aes(
    x = Unadjusted_Spearman_rho,
    y = Adjusted_Partial_Spearman_rho
  ),
  
  shape = 24,
  size = 4.8,
  fill = NA,
  color = "black",
  stroke = 0.85
) +
  
  # 2. Red Positive → Negative triangle
  geom_point(
    data = panelA_df %>%
      filter(PN_Adjusted_Sig),
    
    aes(
      x = Unadjusted_Spearman_rho,
      y = Adjusted_Partial_Spearman_rho
    ),
    
    shape = 24,
    size = 3.1,
    fill = transition_colors[
      "Positive → Negative"
    ],
    color = "white",
    stroke = 0.45
  ) +
  
  # ----------------------------------------------------------
# Direction-transition legend
# ----------------------------------------------------------

scale_shape_manual(
  name = "Direction transition",
  values = transition_shapes,
  breaks = names(transition_shapes)
) +
  
  scale_fill_manual(
    name = "Direction transition",
    values = transition_colors,
    breaks = names(transition_colors)
  ) +
  
  # ----------------------------------------------------------
# Correlation axes
# ----------------------------------------------------------

scale_x_continuous(
  limits = c(-1, 1),
  breaks = seq(-1, 1, by = 0.5),
  expand = expansion(
    mult = c(0.02, 0.02)
  )
) +
  
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5),
    expand = expansion(
      mult = c(0.02, 0.02)
    )
  ) +
  
  # ----------------------------------------------------------
# Axis labels
# ----------------------------------------------------------

labs(
  x = "Unadjusted Spearman ρ",
  y = "ImmuneScore-adjusted partial Spearman ρ"
) +
  
  # ----------------------------------------------------------
# Use one combined Direction-transition legend
# ----------------------------------------------------------

guides(
  
  shape = guide_legend(
    order = 1,
    
    override.aes = list(
      fill = unname(
        transition_colors
      ),
      color = "white",
      size = 3.1,
      stroke = 0.45
    )
  ),
  
  fill = "none"
) +
  
  # ----------------------------------------------------------
# Classic journal-style theme
# ----------------------------------------------------------

theme_classic(
  base_family = "Arial"
) +
  
  theme(
    
    # --------------------------------------------------------
    # Axis titles
    # --------------------------------------------------------
    
    axis.title = element_text(
      family = "Arial",
      face = "bold",
      size = 10
    ),
    
    # --------------------------------------------------------
    # Axis tick labels
    # --------------------------------------------------------
    
    axis.text = element_text(
      family = "Arial",
      size = 8
    ),
    
    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------
    
    legend.title = element_text(
      family = "Arial",
      face = "bold",
      size = 8
    ),
    
    legend.text = element_text(
      family = "Arial",
      size = 8
    ),
    
    legend.key = element_blank(),
    
    legend.position = "right",
    
    # --------------------------------------------------------
    # No panel grid
    # --------------------------------------------------------
    
    panel.grid = element_blank(),
    
    # --------------------------------------------------------
    # Margins
    # --------------------------------------------------------
    
    plot.margin = margin(
      t = 6,
      r = 8,
      b = 6,
      l = 8
    )
  )


# ============================================================
# 11. Add a compact significance explanation
#
# This is NOT a second ggplot legend.
#
# It is simply a visual annotation:
#
# black outer ring + triangle
# =
# adjusted FDR < 0.05
#
# ============================================================

pA <- pA +
  
  annotate(
    "point",
    
    x = 0.48,
    y = -0.87,
    
    shape = 24,
    
    size = 3.4,
    
    fill = NA,
    
    color = "black",
    
    stroke = 0.85
  ) +
  
  annotate(
    "text",
    
    x = 0.57,
    y = -0.87,
    
    label = "adjusted FDR < 0.05",
    
    hjust = 0,
    vjust = 0.5,
    
    family = "Arial",
    
    size = 8 / 3.78,
    
    color = "black"
  )

# ============================================================
# 12. Print Panel A to R graphics device
# ============================================================

print(pA)


# ============================================================
# 13. Export PDF
#
# IMPORTANT:
# Explicit dev.off() is used.
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("Exporting Panel A PDF...\n")
cat("------------------------------------------------------------\n")

grDevices::cairo_pdf(
  file = panelA_pdf,
  width = 7.2,
  height = 6.0
)

print(pA)

grDevices::dev.off()

cat(
  "PDF export: PASS\n"
)

cat(
  "PDF:\n",
  panelA_pdf,
  "\n"
)


# ============================================================
# 14. Export 600-dpi TIFF
#
# IMPORTANT:
# Explicit dev.off() is used.
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("Exporting Panel A 600-dpi TIFF...\n")
cat("------------------------------------------------------------\n")

grDevices::tiff(
  filename = panelA_tiff,
  
  width = 7.2,
  height = 6.0,
  
  units = "in",
  
  res = 600,
  
  compression = "lzw",
  
  type = "cairo"
)

print(pA)

grDevices::dev.off()

cat(
  "600-dpi TIFF export: PASS\n"
)

cat(
  "TIFF:\n",
  panelA_tiff,
  "\n"
)


# ============================================================
# 15. Final status
# ============================================================

cat("\n============================================================\n")
cat("STEP 8C-8 — PANEL A FINAL STATUS\n")
cat("============================================================\n")

cat("Locked source: PASS\n")
cat("Structural QC: PASS\n")
cat("Panel A numerical QC: PASS\n")
cat("54 associations plotted: PASS\n")
cat("21 Positive → Negative adjusted-FDR-significant points: PASS\n")
cat("Direction-transition legend: PASS\n")
cat("Significance annotation: PASS\n")
cat("PDF export: PASS\n")
cat("600-dpi TIFF export: PASS\n")
cat("No new statistical analysis performed.\n")
cat("No rho, P value, or FDR recalculated.\n")

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("\n============================================================\n")
cat("STEP 8C-8 — PANEL A COMPLETE\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 8C-8 — Panel B Visualization Revision v3
#
# Supplementary Figure 4 — Panel B
# Immune-checkpoint adjusted partial Spearman heatmap
#
# Purpose:
#   1. Use LOCKED STEP 8C-7B results only
#   2. Remove duplicated rho colorbars
#   3. Remove X/Y axis titles
#   4. Show AOC3 / MS4A7 / IL7R clearly in Arial bold
#   5. Increase left plot margin
#   6. Preserve WGCNA-style heatmap colors
#   7. Export PDF + 600-dpi TIFF
#   8. Explicit dev.off() after each device
#
# No new statistical analysis
# ============================================================


cat("\n============================================================\n")
cat("STEP 8C-8 — Panel B Visualization Revision v3\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(cowplot)


# ------------------------------------------------------------
# 2. Locked RDS
# ------------------------------------------------------------

rds_path <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "Workflow-08-STEP8C7B_Supplementary_Table_Immune_Checkpoint_54.rds"
)

if (!file.exists(rds_path)) {
  stop("Locked RDS not found:\n", rds_path)
}

supp_table <- readRDS(rds_path)


# ------------------------------------------------------------
# 3. Robustly identify the locked 54-row data frame
# ------------------------------------------------------------

required_columns <- c(
  "Core_Gene",
  "Checkpoint_Gene",
  "N",
  "Unadjusted_Spearman_rho",
  "Unadjusted_P",
  "Unadjusted_FDR_BH",
  "Unadjusted_Direction",
  "Adjusted_Partial_Spearman_rho",
  "Adjusted_P",
  "Adjusted_FDR",
  "Adjusted_Direction",
  "Direction_Transition",
  "Direction_Switch",
  "Unadjusted_FDR_Significant",
  "Adjusted_FDR_Significant",
  "FDR_Transition"
)


if (is.data.frame(supp_table)) {
  
  complete_54 <- supp_table
  
} else if (is.list(supp_table)) {
  
  candidate_idx <- which(
    vapply(
      supp_table,
      function(x) {
        is.data.frame(x) &&
          all(required_columns %in% names(x))
      },
      logical(1)
    )
  )
  
  if (length(candidate_idx) != 1) {
    stop(
      "Could not uniquely identify the locked 54-row data frame.\n",
      "Candidate indices: ",
      paste(candidate_idx, collapse = ", ")
    )
  }
  
  complete_54 <- supp_table[[candidate_idx]]
  
} else {
  
  stop(
    "Unexpected RDS structure: expected data.frame or list."
  )
}


# ------------------------------------------------------------
# 4. Basic QC — LOCKED 54 associations
# ------------------------------------------------------------

if (nrow(complete_54) != 54) {
  stop(
    "QC failed: expected 54 associations, found ",
    nrow(complete_54)
  )
}

if (!all(required_columns %in% names(complete_54))) {
  stop("QC failed: required columns are missing.")
}

if (sum(complete_54$Adjusted_FDR < 0.05, na.rm = TRUE) != 38) {
  stop(
    "QC failed: expected 38 adjusted-FDR-significant associations."
  )
}


cat("QC PASS: 54 associations loaded.\n")
cat("QC PASS: 38 adjusted-FDR-significant associations.\n\n")


# ------------------------------------------------------------
# 5. Prepare Panel B data
# ------------------------------------------------------------

panelB_df <- complete_54 %>%
  mutate(
    FDR_Stars = case_when(
      Adjusted_FDR < 0.001 ~ "***",
      Adjusted_FDR < 0.01  ~ "**",
      Adjusted_FDR < 0.05  ~ "*",
      TRUE ~ ""
    )
  )


# ------------------------------------------------------------
# 6. Lock row and column order
# ------------------------------------------------------------

core_gene_order <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

checkpoint_order <- unique(panelB_df$Checkpoint_Gene)

if (length(checkpoint_order) != 18) {
  stop(
    "QC failed: expected 18 checkpoint genes, found ",
    length(checkpoint_order)
  )
}


panelB_df <- panelB_df %>%
  mutate(
    Core_Gene = factor(
      Core_Gene,
      levels = rev(core_gene_order)
    ),
    Checkpoint_Gene = factor(
      Checkpoint_Gene,
      levels = checkpoint_order
    )
  )


# ------------------------------------------------------------
# 7. WGCNA-style heatmap colors
# ------------------------------------------------------------

rho_low_color  <- "#82B5EC"
rho_mid_color  <- "white"
rho_high_color <- "#F57C6E"


# ------------------------------------------------------------
# 8. Main heatmap
#
# IMPORTANT:
#   No X/Y axis titles.
#   No ggplot fill legend.
#   A separate custom rho colorbar is created below.
# ------------------------------------------------------------

p_heatmap <- ggplot(
  panelB_df,
  aes(
    x = Checkpoint_Gene,
    y = Core_Gene
  )
) +
  
  geom_tile(
    aes(
      fill = Adjusted_Partial_Spearman_rho
    ),
    color = "white",
    linewidth = 0.35
  ) +
  
  geom_text(
    aes(
      label = FDR_Stars
    ),
    family = "Arial",
    fontface = "bold",
    size = 8 / 3.78,
    color = "black"
  ) +
  
  scale_fill_gradient2(
    low = rho_low_color,
    mid = rho_mid_color,
    high = rho_high_color,
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    guide = "none"
  ) +
  
  scale_x_discrete(
    drop = FALSE
  ) +
  
  scale_y_discrete(
    drop = FALSE
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(
    base_family = "Arial"
  ) +
  
  theme(
    
    # --------------------------------------------------------
    # IMPORTANT:
    # Large left margin so AOC3 / MS4A7 / IL7R are never clipped
    # --------------------------------------------------------
    
    plot.margin = margin(
      t = 6,
      r = 8,
      b = 8,
      l = 28
    ),
    
    # --------------------------------------------------------
    # X-axis: checkpoint genes
    # --------------------------------------------------------
    
    axis.text.x = element_text(
      family = "Arial",
      face = "bold",
      size = 8,
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    
    # --------------------------------------------------------
    # Y-axis: core genes
    # --------------------------------------------------------
    
    axis.text.y = element_text(
      family = "Arial",
      face = "bold",
      size = 9,
      color = "black"
    ),
    
    # No axis titles
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    # No grid
    panel.grid = element_blank(),
    
    # No legend
    legend.position = "none"
  )


# ------------------------------------------------------------
# 9. PANEL B — FDR LEGEND
#
# Layout:
#   right side — top
# ------------------------------------------------------------

p_B_fdr <- ggplot() +
  
  annotate(
    "text",
    x = 0,
    y = 1.00,
    label = "FDR",
    hjust = 0,
    vjust = 1,
    family = "Arial",
    fontface = "bold",
    size = 3.1
  ) +
  
  annotate(
    "text",
    x = 0,
    y = 0.70,
    label = "***  FDR < 0.001",
    hjust = 0,
    vjust = 1,
    family = "Arial",
    size = 2.7
  ) +
  
  annotate(
    "text",
    x = 0,
    y = 0.43,
    label = "**   FDR < 0.01",
    hjust = 0,
    vjust = 1,
    family = "Arial",
    size = 2.7
  ) +
  
  annotate(
    "text",
    x = 0,
    y = 0.16,
    label = "*    FDR < 0.05",
    hjust = 0,
    vjust = 1,
    family = "Arial",
    size = 2.7
  ) +
  
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off"
  ) +
  
  theme_void(
    base_family = "Arial"
  ) +
  
  theme(
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  )

# ============================================================
# Panel B: rho colorbar
# Further adjustment:
#   1. Narrower colorbar
#   2. Keep tick marks close to colorbar
#   3. Keep numeric labels close to tick marks
# ============================================================

rho_edges <- seq(-1, 1, length.out = 102)

rho_bar_df <- data.frame(
  xmin = 0.00,
  xmax = 0.20,      
  ymin = rho_edges[-length(rho_edges)],
  ymax = rho_edges[-1],
  rho  = seq(-1, 1, length.out = 101)
)

rho_ticks <- data.frame(
  rho = c(-1, -0.5, 0, 0.5, 1)
)

p_B_rho <- ggplot() +
  
  # ----------------------------------------------------------
# 1. rho colorbar
# ----------------------------------------------------------
geom_rect(
  data = rho_bar_df,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    fill = rho
  ),
  color = NA
) +
  
  scale_fill_gradient2(
    low = rho_low_color,
    mid = rho_mid_color,
    high = rho_high_color,
    midpoint = 0,
    limits = c(-1, 1),
    guide = "none"
  ) +
  
  # ----------------------------------------------------------
# 2. Tick marks: directly attached to the colorbar
# ----------------------------------------------------------
geom_segment(
  data = rho_ticks,
  aes(
    x = 0.20,       
    xend = 0.215,   
    y = rho,
    yend = rho
  ),
  linewidth = 0.30,
  color = "black"
) +
  
  # ----------------------------------------------------------
# 3. Numeric labels: slightly separated from tick marks
# ----------------------------------------------------------
geom_text(
  data = rho_ticks,
  aes(
    x = 0.235,      
    y = rho,
    label = rho
  ),
  family = "Arial",
  size = 8 / ggplot2::.pt,
  color = "black",
  hjust = 0,
  vjust = 0.5
) +
  
  # ----------------------------------------------------------
# 4. Fixed compact coordinate system
# ----------------------------------------------------------
coord_cartesian(
  xlim = c(0, 0.55),
  ylim = c(-1, 1),
  expand = FALSE,
  clip = "off"
) +
  
  theme_void(base_family = "Arial") +
  
  theme(
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 4,
      l = 0
    )
  )


# ------------------------------------------------------------
# 12. PANEL B — rho section
#
# rho title + rho colorbar
# ------------------------------------------------------------
p_B_rho_section <- cowplot::plot_grid(
  p_B_rho_title,
  p_B_rho,
  ncol = 1,
  align = "v",
  rel_heights = c(0.28, 1)
)

# ------------------------------------------------------------
# 13. PANEL B — right-side legend
#
# IMPORTANT:
#   FDR is ABOVE
#   rho colorbar is BELOW
# ------------------------------------------------------------

p_B_right <- cowplot::plot_grid(
  p_B_fdr,
  p_B_rho_section,
  ncol = 1,
  align = "v",
  rel_heights = c(
    0.62,
    1.55
  )
)


# ------------------------------------------------------------
# 14. PANEL B — final composition
# ------------------------------------------------------------

panel_B <- cowplot::plot_grid(
  p_heatmap,
  p_B_right,
  ncol = 2,
  align = "h",
  rel_widths = c(
    1,
    0.18
  )
)


# ------------------------------------------------------------
# 15. Display
# ------------------------------------------------------------

print(panel_B)

# ------------------------------------------------------------
# 15. Export PDF
# ------------------------------------------------------------

grDevices::cairo_pdf(
  file = pdf_path,
  width = 9.0,
  height = 5.8
)

print(panel_B)

grDevices::dev.off()


# ------------------------------------------------------------
# 16. Export 600-dpi TIFF
# ------------------------------------------------------------

grDevices::tiff(
  filename = tiff_path,
  width = 9.0,
  height = 5.8,
  units = "in",
  res = 600,
  compression = "lzw",
  type = "cairo"
)

print(panel_B)

grDevices::dev.off()


# ------------------------------------------------------------
# 17. Output QC
# ------------------------------------------------------------

if (!file.exists(pdf_path)) {
  stop("PDF export failed.")
}

if (!file.exists(tiff_path)) {
  stop("TIFF export failed.")
}


cat("\n============================================================\n")
cat("STEP 8C-8 — Panel B FINAL STATUS: PASS\n")
cat("============================================================\n")
cat("Locked associations: 54\n")
cat("Adjusted FDR < 0.05: 38\n")
cat("Core genes: AOC3, MS4A7, IL7R\n")
cat("Checkpoint genes: 18\n")
cat("Axis titles: REMOVED\n")
cat("Default rho legend: REMOVED\n")
cat("Custom rho colorbar: ONE\n")
cat("Core-gene labels: Arial bold\n")
cat("Left plot margin: 28\n")
cat("PDF: ", pdf_path, "\n", sep = "")
cat("TIFF: ", tiff_path, "\n", sep = "")
cat("============================================================\n")


# ============================================================
# Paper 1
# END-OF-DAY WORKSPACE CHECKPOINT
# Save current R workspace for continuation tomorrow
# Date: 
# ============================================================

cat("\n============================================================\n")
cat("Paper 1 — End-of-day R workspace checkpoint\n")
cat("Date: \n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define checkpoint directory
# ------------------------------------------------------------

checkpoint_dir <- "input_dir/Paper1_tumor-only WGCNA"

if (!dir.exists(checkpoint_dir)) {
  dir.create(
    checkpoint_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

cat("Checkpoint directory:\n")
cat(checkpoint_dir, "\n\n")


# ------------------------------------------------------------
# 2. Timestamp
# ------------------------------------------------------------

checkpoint_time <- format(
  Sys.time(),
  "%Y-%m-%d %H:%M:%S"
)

timestamp_file <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

cat("Checkpoint time:\n")
cat(checkpoint_time, "\n\n")


# ------------------------------------------------------------
# 3. Define output files
# ------------------------------------------------------------

workspace_file <- file.path(
  checkpoint_dir,
  paste0(
    "Paper1_END_WORKSPACE_20260920_",
    timestamp_file,
    ".RData"
  )
)

session_file <- file.path(
  checkpoint_dir,
  paste0(
    "Paper1_END_SESSIONINFO_20260920_",
    timestamp_file,
    ".txt"
  )
)

objects_file <- file.path(
  checkpoint_dir,
  paste0(
    "Paper1_END_OBJECT_INVENTORY_20260920_",
    timestamp_file,
    ".txt"
  )
)

readme_file <- file.path(
  checkpoint_dir,
  paste0(
    "Paper1_END_CHECKPOINT_README_20260920_",
    timestamp_file,
    ".txt"
  )
)


# ------------------------------------------------------------
# 4. Save complete current workspace
# ------------------------------------------------------------

cat("Saving complete current R workspace...\n")

save.image(
  file = workspace_file
)

cat("Workspace saved:\n")
cat(workspace_file, "\n\n")


# ------------------------------------------------------------
# 5. Save session information
# ------------------------------------------------------------

cat("Saving session information...\n")

sink(session_file)

cat("Paper 1 — R session information\n")
cat("================================\n\n")

cat("Checkpoint date/time:\n")
cat(checkpoint_time, "\n\n")

cat("R version and platform:\n")
print(R.version)

cat("\n\nsessionInfo():\n")
print(sessionInfo())

sink()

cat("Session information saved:\n")
cat(session_file, "\n\n")


# ------------------------------------------------------------
# 6. Save complete object inventory
# ------------------------------------------------------------

cat("Creating object inventory...\n")

all_objects <- ls(envir = .GlobalEnv)

sink(objects_file)

cat("Paper 1 — Complete R workspace object inventory\n")
cat("================================================\n\n")

cat("Checkpoint date/time:\n")
cat(checkpoint_time, "\n\n")

cat("Number of objects:\n")
cat(length(all_objects), "\n\n")

cat("Objects:\n")
print(all_objects)

cat("\n\nObject classes and dimensions:\n")
cat("--------------------------------\n")

for (obj_name in all_objects) {
  
  obj <- get(
    obj_name,
    envir = .GlobalEnv
  )
  
  cat("\nOBJECT:", obj_name, "\n")
  cat("CLASS :", paste(class(obj), collapse = ", "), "\n")
  
  if (is.data.frame(obj) ||
      is.matrix(obj) ||
      is.array(obj)) {
    
    cat(
      "DIM   :",
      paste(dim(obj), collapse = " x "),
      "\n"
    )
    
  } else if (is.vector(obj)) {
    
    cat(
      "LENGTH:",
      length(obj),
      "\n"
    )
  }
}

sink()

cat("Object inventory saved:\n")
cat(objects_file, "\n\n")


# ------------------------------------------------------------
# 7. Create continuation README
# ------------------------------------------------------------

cat("Creating checkpoint README...\n")

sink(readme_file)

cat("Paper 1 — End-of-day checkpoint\n")
cat("================================\n\n")

cat("Date/time:\n")
cat(checkpoint_time, "\n\n")

cat("Purpose:\n")
cat(
  "Preserve the complete R workspace after today's completed analysis ",
  "so that the Paper 1 revision workflow can continue tomorrow without ",
  "rerunning completed analyses.\n\n",
  sep = ""
)

cat("Workspace checkpoint:\n")
cat(basename(workspace_file), "\n\n")

cat("Session information:\n")
cat(basename(session_file), "\n\n")

cat("Object inventory:\n")
cat(basename(objects_file), "\n\n")

cat("Important continuation rules:\n")
cat("--------------------------------\n")
cat(
  "1. Load this checkpoint before continuing tomorrow.\n"
)

cat(
  "2. Do not rerun already locked analyses unless specifically required.\n"
)

cat(
  "3. Preserve locked Paper 1 results and plotting objects.\n"
)

cat(
  "4. STEP 8C-7A/7B immune-checkpoint primary results remain locked.\n"
)

cat(
  "5. Any subsequent purity or circularity analysis should be treated ",
  "as a sensitivity analysis and should not overwrite the locked primary analysis.\n",
  sep = ""
)

cat(
  "6. Continue from the actual end point of today's workflow rather ",
  "than restarting the workflow from Step 1.\n",
  sep = ""
)

sink()

cat("README saved:\n")
cat(readme_file, "\n\n")


# ------------------------------------------------------------
# 8. Verify files exist
# ------------------------------------------------------------

cat("============================================================\n")
cat("CHECKPOINT VERIFICATION\n")
cat("============================================================\n\n")

checkpoint_files <- c(
  workspace_file,
  session_file,
  objects_file,
  readme_file
)

verification <- file.exists(checkpoint_files)

verification_table <- data.frame(
  file = basename(checkpoint_files),
  exists = verification,
  size_MB = round(
    file.info(checkpoint_files)$size / 1024^2,
    3
  ),
  stringsAsFactors = FALSE
)

print(verification_table)

cat("\n")


# ------------------------------------------------------------
# 9. Final status
# ------------------------------------------------------------

if (all(verification)) {
  
  cat("============================================================\n")
  cat("FINAL STATUS: PASS\n")
  cat("============================================================\n\n")
  
  cat(
    "Today's R workspace has been successfully preserved.\n"
  )
  
  cat(
    "Tomorrow, load the saved .RData checkpoint and continue ",
    "from the current workflow state.\n",
    sep = ""
  )
  
} else {
  
  cat("============================================================\n")
  cat("FINAL STATUS: CHECK REQUIRED\n")
  cat("============================================================\n\n")
  
  cat(
    "One or more checkpoint files were not created successfully.\n"
  )
}


# ------------------------------------------------------------
# 10. Print the exact workspace file for tomorrow
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("TOMORROW'S RESTORE FILE\n")
cat("============================================================\n\n")

cat(workspace_file, "\n\n")






