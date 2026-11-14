############################################################
# Paper 1: LUAD transcriptomic analysis
#
# Script:
#   GitHub2_Paper1_tumor_only_WGCNA_1.R
#
# Workflow:
#   Workflow-04-1
#   Clinical cohort reconstruction and clinical data cleaning
#
# Research purpose:
#   1. Identify TCGA Primary Tumor samples (sample type 01)
#      from the original TCGA FPKM-UQ expression matrix.
#   2. Perform exact sample-level matching between expression
#      sample barcodes and the clinical data.
#   3. Reconstruct the Primary Tumor clinical cohort.
#   4. Audit clinical variable completeness and missing-value
#      encodings.
#   5. Standardize recognized missing-value encodings and
#      clean selected clinical variables.
#   6. Assess overall-survival data consistency without
#      automatically excluding samples.
#
# Input data:
#   1. data/Merge_RNA_seq_FPKM_UQ.txt
#      Original TCGA LUAD FPKM-UQ expression matrix.
#
#   2. data/<actual clinical data filename>
#      TCGA clinical data used in the original analysis.
#      The filename must be specified in the configuration
#      section below.
#
# Input data provenance:
#   Original TCGA LUAD expression and clinical data used in
#   the Paper 1 analysis.
#   The original input files are read-only and are not modified.
#
# Sample inclusion criteria:
#   - Standard TCGA sample barcode.
#   - TCGA sample type code = 01 (Primary Tumor).
#   - Exact sample-level match to the clinical Sample field.
#   - One Primary Tumor sample per patient in the locked cohort.
#
# Expected cohort:
#   513 Primary Tumor samples.
#
# Important:
#   - No expression filtering is performed here.
#   - No WGCNA is performed here.
#   - No clinical samples are removed because of missing
#     age, stage, OS.time, or OS.status.
#   - No undocumented recoding of survival status is allowed.
#   - No causal or diagnostic interpretation is made.
#
# Software environment:
#   R 4.5.1 (project target version).
#   Packages: here, readxl (versions recorded at runtime).
#
# Random seed:
#   Not applicable. This workflow uses deterministic operations
#   and does not perform random sampling or model fitting.
#
# Output directory:
#   results/02_tumor_only_WGCNA/Workflow-04-1/
#
# Downstream handoff:
#   Workflow-04-2A uses:
#   - Workflow-04-1E_clinical_clean.txt
#   - Workflow-04-1F-1_primary_tumor_sample_list.txt
#
#   Workflow-04-2A must independently align the clinical table
#   to the expression matrix sample-column order.
############################################################


############################################################
# 00. Project configuration and setup
############################################################

required_packages <- c(
  "here",
  "readxl"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "The following required packages are missing: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this analysis script."
  )
}

project_root <- here::here()

data_dir <- file.path(
  project_root,
  "data"
)

result_dir <- file.path(
  project_root,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-1"
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 00-1. Input file configuration
############################################################

expression_file <- file.path(
  data_dir,
  "Merge_RNA_seq_FPKM_UQ.txt"
)

# IMPORTANT:
# Replace this placeholder with the actual clinical filename
# located in the repository's data/ directory.
#
# Do not rename or modify the original clinical input file
# solely to satisfy this script.

clinical_filename <- "REPLACE_WITH_ACTUAL_CLINICAL_FILENAME.xlsx"

clinical_file <- file.path(
  data_dir,
  clinical_filename
)

# Locked project cohort size.
expected_primary_tumor_n <- 513L

# Required clinical fields.
required_clinical_vars <- c(
  "Sample",
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


############################################################
# 00-2. Input file checks
############################################################

if (!file.exists(expression_file)) {
  stop(
    "Expression input file was not found:\n",
    expression_file
  )
}

if (
  grepl(
    "^REPLACE_WITH_",
    clinical_filename
  )
) {
  stop(
    "Please replace clinical_filename with the actual ",
    "clinical input filename before running this script."
  )
}

if (!file.exists(clinical_file)) {
  stop(
    "Clinical input file was not found:\n",
    clinical_file,
    "\nPlease check the filename and data directory."
  )
}


############################################################
# 00-3. Software environment record
############################################################

r_version <- as.character(
  getRversion()
)

package_versions <- data.frame(
  Package = required_packages,
  Version = vapply(
    required_packages,
    function(pkg) {
      as.character(
        utils::packageVersion(pkg)
      )
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)

cat("\n========================================\n")
cat("SOFTWARE ENVIRONMENT\n")
cat("========================================\n\n")

cat("R version:", r_version, "\n\n")

print(
  package_versions,
  row.names = FALSE
)

if (r_version != "4.5.1") {
  warning(
    "The project target R version is 4.5.1. ",
    "Current version: ",
    r_version
  )
}

write.table(
  package_versions,
  file = file.path(
    result_dir,
    "Workflow-04-1_software_versions.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 01. Read expression matrix header
#
# Only the header is required in this workflow.
# The full expression matrix is read later in Workflow-04-2A.
############################################################

cat("\n========================================\n")
cat("EXPRESSION MATRIX HEADER AUDIT\n")
cat("========================================\n\n")

expr_header <- read.delim(
  expression_file,
  header = TRUE,
  sep = "\t",
  nrows = 0,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(expr_header) < 2L) {
  stop(
    "The expression matrix must contain one gene-ID column ",
    "and at least one sample column."
  )
}

expression_column_names <- colnames(
  expr_header
)

if (
  anyNA(expression_column_names) ||
  any(trimws(expression_column_names) == "")
) {
  stop(
    "The expression matrix contains missing or blank ",
    "column names."
  )
}

sample_ids <- expression_column_names[-1]

if (anyDuplicated(sample_ids)) {
  duplicated_samples <- unique(
    sample_ids[
      duplicated(sample_ids)
    ]
  )
  
  write.table(
    data.frame(
      Sample = duplicated_samples
    ),
    file = file.path(
      result_dir,
      "Workflow-04-1A_duplicate_expression_sample_IDs.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  stop(
    "Duplicated expression sample IDs were detected. ",
    "See the duplicate-ID audit file."
  )
}

cat(
  "Expression sample columns:",
  length(sample_ids),
  "\n"
)

cat(
  "First expression column (gene ID field):",
  expression_column_names[1],
  "\n"
)


############################################################
# 02. Validate TCGA sample barcodes
############################################################

cat("\n========================================\n")
cat("TCGA SAMPLE BARCODE AUDIT\n")
cat("========================================\n\n")

tcga_barcode_pattern <- paste0(
  "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}$"
)

valid_barcode <- grepl(
  tcga_barcode_pattern,
  sample_ids
)

invalid_expression_ids <- sample_ids[
  !valid_barcode
]

if (length(invalid_expression_ids) > 0) {
  
  write.table(
    data.frame(
      Sample = invalid_expression_ids
    ),
    file = file.path(
      result_dir,
      "Workflow-04-1A_invalid_expression_sample_IDs.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  stop(
    "Non-standard TCGA sample barcodes were detected. ",
    "Please inspect the invalid sample-ID audit file."
  )
}

sample_type_code <- substr(
  sample_ids,
  14,
  15
)

sample_type_table <- as.data.frame(
  table(
    sample_type_code,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

colnames(sample_type_table) <- c(
  "Sample_type_code",
  "N"
)

cat("Original sample composition:\n\n")

print(
  sample_type_table,
  row.names = FALSE
)

write.table(
  sample_type_table,
  file = file.path(
    result_dir,
    "Workflow-04-1A_expression_sample_composition.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 03. Identify Primary Tumor samples
#
# TCGA sample type:
#   01 = Primary Tumor
############################################################

primary_tumor_ids <- sample_ids[
  sample_type_code == "01"
]

if (length(primary_tumor_ids) == 0L) {
  stop(
    "No Primary Tumor samples (sample type 01) were detected."
  )
}

primary_tumor_patient_ids <- substr(
  primary_tumor_ids,
  1,
  12
)

patient_frequency <- table(
  primary_tumor_patient_ids
)

multiple_primary_tumor_patients <- names(
  patient_frequency[
    patient_frequency > 1
  ]
)

patient_duplication_audit <- data.frame(
  Patient_ID = names(patient_frequency),
  N_Primary_Tumor_Samples = as.integer(
    patient_frequency
  ),
  stringsAsFactors = FALSE
)

write.table(
  patient_duplication_audit,
  file = file.path(
    result_dir,
    "Workflow-04-1A_primary_tumor_patient_sample_count.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nPrimary Tumor samples:",
  length(primary_tumor_ids),
  "\n"
)

cat(
  "Unique Primary Tumor patients:",
  length(unique(primary_tumor_patient_ids)),
  "\n"
)

cat(
  "Patients with multiple Primary Tumor samples:",
  length(multiple_primary_tumor_patients),
  "\n"
)

if (length(multiple_primary_tumor_patients) > 0L) {
  stop(
    "More than one Primary Tumor sample was detected ",
    "for at least one patient. ",
    "The one-sample-per-patient cohort assumption ",
    "requires review before proceeding."
  )
}

if (
  length(primary_tumor_ids) !=
  expected_primary_tumor_n
) {
  stop(
    "The number of Primary Tumor samples does not match ",
    "the locked project cohort size.\n",
    "Expected: ",
    expected_primary_tumor_n,
    "\nObserved: ",
    length(primary_tumor_ids),
    "\nPlease investigate the input data and sample selection."
  )
}


############################################################
# 04. Read clinical data
############################################################

clinical_extension <- tolower(
  tools::file_ext(clinical_file)
)

if (clinical_extension %in% c("xlsx", "xls")) {
  
  clinical <- as.data.frame(
    readxl::read_excel(
      clinical_file
    ),
    stringsAsFactors = FALSE
  )
  
} else if (clinical_extension == "csv") {
  
  clinical <- read.csv(
    clinical_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
} else if (clinical_extension %in% c("txt", "tsv")) {
  
  clinical <- read.delim(
    clinical_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
} else {
  
  stop(
    "Unsupported clinical file extension: ",
    clinical_extension,
    "\nSupported formats: xlsx, xls, csv, txt, tsv."
  )
}

cat("\n========================================\n")
cat("CLINICAL DATA STRUCTURE\n")
cat("========================================\n\n")

cat("Clinical rows:", nrow(clinical), "\n")
cat("Clinical columns:", ncol(clinical), "\n\n")

print(
  colnames(clinical)
)

missing_required_columns <- setdiff(
  required_clinical_vars,
  colnames(clinical)
)

if (length(missing_required_columns) > 0L) {
  stop(
    "Required clinical columns are missing:\n",
    paste(
      missing_required_columns,
      collapse = ", "
    )
  )
}


############################################################
# 05. Clean and validate clinical sample IDs
############################################################

clinical$Sample <- trimws(
  as.character(
    clinical$Sample
  )
)

clinical$Sample[
  is.na(clinical$Sample) |
    clinical$Sample == ""
] <- NA_character_

clinical_missing_sample_id <- is.na(
  clinical$Sample
)

clinical_valid_ids <- clinical$Sample[
  !clinical_missing_sample_id
]

clinical_duplicate_ids <- unique(
  clinical_valid_ids[
    duplicated(clinical_valid_ids)
  ]
)

clinical_id_audit <- data.frame(
  Metric = c(
    "Clinical_rows",
    "Clinical_missing_or_blank_Sample_ID",
    "Clinical_nonmissing_Sample_ID",
    "Unique_clinical_Sample_ID",
    "Duplicated_clinical_Sample_ID"
  ),
  Value = c(
    nrow(clinical),
    sum(clinical_missing_sample_id),
    length(clinical_valid_ids),
    length(unique(clinical_valid_ids)),
    length(clinical_duplicate_ids)
  ),
  stringsAsFactors = FALSE
)

write.table(
  clinical_id_audit,
  file = file.path(
    result_dir,
    "Workflow-04-1B_clinical_sample_ID_audit.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (length(clinical_duplicate_ids) > 0L) {
  
  write.table(
    data.frame(
      Sample = clinical_duplicate_ids
    ),
    file = file.path(
      result_dir,
      "Workflow-04-1B_duplicate_clinical_sample_IDs.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  stop(
    "Duplicated clinical Sample IDs were detected. ",
    "Exact one-to-one matching cannot be guaranteed."
  )
}


############################################################
# 06. Bidirectional exact sample-level matching
############################################################

cat("\n========================================\n")
cat("BIDIRECTIONAL SAMPLE MATCHING\n")
cat("========================================\n\n")

# Expression -> clinical
expression_match_index <- match(
  sample_ids,
  clinical$Sample
)

expression_matched <- !is.na(
  expression_match_index
)

unmatched_expression <- sample_ids[
  !expression_matched
]

# Clinical -> expression
clinical_in_expression <- (
  !is.na(clinical$Sample) &
    clinical$Sample %in% sample_ids
)

unmatched_clinical <- clinical$Sample[
  !is.na(clinical$Sample) &
    !clinical_in_expression
]

unmatched_clinical <- unique(
  unmatched_clinical
)

# Primary Tumor matching
primary_match_index <- match(
  primary_tumor_ids,
  clinical$Sample
)

unmatched_primary_tumor <- primary_tumor_ids[
  is.na(primary_match_index)
]

matching_summary <- data.frame(
  Metric = c(
    "Expression_sample_count",
    "Expression_samples_matched_to_clinical",
    "Expression_samples_unmatched",
    "Clinical_rows",
    "Clinical_samples_matched_to_expression",
    "Clinical_samples_not_in_expression",
    "Primary_Tumor_sample_count",
    "Primary_Tumor_samples_unmatched"
  ),
  Value = c(
    length(sample_ids),
    sum(expression_matched),
    length(unmatched_expression),
    nrow(clinical),
    sum(clinical_in_expression),
    length(unmatched_clinical),
    length(primary_tumor_ids),
    length(unmatched_primary_tumor)
  ),
  stringsAsFactors = FALSE
)

print(
  matching_summary,
  row.names = FALSE
)

write.table(
  matching_summary,
  file = file.path(
    result_dir,
    "Workflow-04-1B_expression_clinical_matching_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    Sample = unmatched_expression
  ),
  file = file.path(
    result_dir,
    "Workflow-04-1B_unmatched_expression_samples.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    Sample = unmatched_clinical
  ),
  file = file.path(
    result_dir,
    "Workflow-04-1B_unmatched_clinical_samples.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (length(unmatched_primary_tumor) > 0L) {
  
  write.table(
    data.frame(
      Sample = unmatched_primary_tumor
    ),
    file = file.path(
      result_dir,
      "Workflow-04-1B_unmatched_primary_tumor_samples.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  stop(
    "At least one Primary Tumor sample has no exact ",
    "clinical match. The cohort cannot be finalized."
  )
}


############################################################
# 07. Construct matched clinical table
#
# Preserve expression sample order during matching.
# Only Primary Tumor samples are retained for this workflow.
############################################################

clinical_matched_sample <- clinical[
  primary_match_index,
  ,
  drop = FALSE
]

rownames(
  clinical_matched_sample
) <- primary_tumor_ids

if (
  nrow(clinical_matched_sample) !=
  length(primary_tumor_ids)
) {
  stop(
    "The matched clinical table has an unexpected row count."
  )
}

if (
  !identical(
    clinical_matched_sample$Sample,
    primary_tumor_ids
  )
) {
  stop(
    "Clinical sample order does not match the Primary Tumor ",
    "expression sample order."
  )
}

primary_tumor_clinical <- clinical_matched_sample

if (
  nrow(primary_tumor_clinical) !=
  expected_primary_tumor_n
) {
  stop(
    "The matched Primary Tumor clinical cohort does not ",
    "contain the expected 513 samples."
  )
}


############################################################
# 08. Verify tumor classification and patient IDs
############################################################

primary_tumor_clinical$tumor_normal <- tolower(
  trimws(
    as.character(
      primary_tumor_clinical$tumor_normal
    )
  )
)

primary_tumor_clinical$tumor_normal[
  primary_tumor_clinical$tumor_normal == ""
] <- NA_character_

tumor_classification_table <- as.data.frame(
  table(
    primary_tumor_clinical$tumor_normal,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

colnames(tumor_classification_table) <- c(
  "tumor_normal",
  "N"
)

write.table(
  tumor_classification_table,
  file = file.path(
    result_dir,
    "Workflow-04-1B_primary_tumor_classification.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (
  any(
    is.na(primary_tumor_clinical$tumor_normal) |
    primary_tumor_clinical$tumor_normal != "tumor"
  )
) {
  stop(
    "The clinical tumor_normal field is inconsistent ",
    "with the Primary Tumor sample barcode classification. ",
    "Please inspect the classification audit."
  )
}

expected_patient_ids <- substr(
  primary_tumor_clinical$Sample,
  1,
  12
)

observed_patient_ids <- trimws(
  as.character(
    primary_tumor_clinical$patient_id
  )
)

patient_id_missing <- (
  is.na(observed_patient_ids) |
    observed_patient_ids == ""
)

patient_id_mismatch <- (
  !patient_id_missing &
    observed_patient_ids != expected_patient_ids
)

patient_id_audit <- data.frame(
  Sample = primary_tumor_clinical$Sample,
  Expected_patient_ID = expected_patient_ids,
  Clinical_patient_ID = observed_patient_ids,
  Patient_ID_missing = patient_id_missing,
  Patient_ID_mismatch = patient_id_mismatch,
  stringsAsFactors = FALSE
)

write.table(
  patient_id_audit,
  file = file.path(
    result_dir,
    "Workflow-04-1B_patient_ID_consistency_audit.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 09. Missing-value definitions and helper functions
############################################################

missing_tokens <- c(
  "",
  "NA",
  "N/A",
  "NULL",
  "-",
  "."
)

is_missing_value <- function(x) {
  
  x_chr <- trimws(
    as.character(x)
  )
  
  is.na(x) |
    toupper(x_chr) %in% missing_tokens
}


clean_character_missing <- function(x) {
  
  x_chr <- trimws(
    as.character(x)
  )
  
  missing_idx <- (
    is.na(x) |
      toupper(x_chr) %in% missing_tokens
  )
  
  x_chr[missing_idx] <- NA_character_
  
  x_chr
}

parse_numeric_audited <- function(x, variable_name) {
  
  raw_chr <- trimws(
    as.character(x)
  )
  
  missing_idx <- (
    is.na(x) |
      toupper(raw_chr) %in% missing_tokens
  )
  
  parsed <- suppressWarnings(
    as.numeric(raw_chr)
  )
  
  nonfinite_idx <- (
    !missing_idx &
      !is.na(parsed) &
      !is.finite(parsed)
  )
  
  invalid_idx <- (
    !missing_idx &
      (
        is.na(parsed) |
          !is.finite(parsed)
      )
  )
  
  invalid_values <- unique(
    raw_chr[invalid_idx]
  )
  
  parsed[missing_idx | invalid_idx] <- NA_real_
  
  list(
    value = parsed,
    invalid = data.frame(
      Variable = variable_name,
      Invalid_value = invalid_values,
      stringsAsFactors = FALSE
    ),
    invalid_count = sum(invalid_idx),
    nonfinite = nonfinite_idx
  )
}

############################################################
# 10. Clinical completeness audit before cleaning
############################################################

clinical_vars_to_audit <- c(
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

clinical_missing_summary <- do.call(
  rbind,
  lapply(
    clinical_vars_to_audit,
    function(variable_name) {
      
      x <- primary_tumor_clinical[[variable_name]]
      
      n_missing <- sum(
        is_missing_value(x)
      )
      
      data.frame(
        Variable = variable_name,
        N = length(x),
        Missing_N = n_missing,
        Missing_Percent = (
          100 * n_missing / length(x)
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

cat("\n========================================\n")
cat("PRIMARY TUMOR CLINICAL COMPLETENESS\n")
cat("========================================\n\n")

print(
  clinical_missing_summary,
  row.names = FALSE
)

write.table(
  clinical_missing_summary,
  file = file.path(
    result_dir,
    "Workflow-04-1C_clinical_missing_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 11. Clinical variable distributions before cleaning
############################################################

distribution_vars <- c(
  "Gender",
  "OS.status",
  "pathologic_stage",
  "pathologic_stage_2",
  "pathologic_t",
  "pathologic_n",
  "pathologic_m"
)

for (variable_name in distribution_vars) {
  
  cat("\n========================================\n")
  cat(variable_name, " distribution\n")
  cat("========================================\n\n")
  
  print(
    table(
      primary_tumor_clinical[[variable_name]],
      useNA = "ifany"
    )
  )
}


############################################################
# 12. Detailed missing-value encoding audit
############################################################

audit_missing_encoding <- function(x) {
  
  x_chr <- trimws(
    as.character(x)
  )
  
  data.frame(
    True_NA = sum(is.na(x)),
    Empty = sum(
      !is.na(x) & x_chr == ""
    ),
    Literal_NA = sum(
      !is.na(x) & toupper(x_chr) == "NA"
    ),
    N_A = sum(
      !is.na(x) & toupper(x_chr) == "N/A"
    ),
    NULL = sum(
      !is.na(x) & toupper(x_chr) == "NULL"
    ),
    Dash = sum(
      !is.na(x) & x_chr == "-"
    ),
    Dot = sum(
      !is.na(x) & x_chr == "."
    ),
    stringsAsFactors = FALSE
  )
}

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

missing_encoding_audit <- do.call(
  rbind,
  lapply(
    audit_vars,
    function(variable_name) {
      
      result <- audit_missing_encoding(
        primary_tumor_clinical[[variable_name]]
      )
      
      result$Variable <- variable_name
      
      result
    }
  )
)

missing_encoding_audit <- missing_encoding_audit[
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

write.table(
  missing_encoding_audit,
  file = file.path(
    result_dir,
    "Workflow-04-1D_missing_value_encoding_audit.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 13. Create a cleaned clinical copy
#
# Do not overwrite primary_tumor_clinical.
# Do not delete samples because of missing clinical values.
############################################################

clinical_clean <- primary_tumor_clinical

# Clean character columns only.
character_columns <- names(clinical_clean)[
  vapply(
    clinical_clean,
    is.character,
    logical(1)
  )
]

for (variable_name in character_columns) {
  
  clinical_clean[[variable_name]] <- (
    clean_character_missing(
      clinical_clean[[variable_name]]
    )
  )
}


############################################################
# 14. Convert Age and OS.time with explicit audit
############################################################

age_result <- parse_numeric_audited(
  clinical_clean$Age,
  "Age"
)

os_time_result <- parse_numeric_audited(
  clinical_clean$OS.time,
  "OS.time"
)

clinical_clean$Age <- age_result$value
clinical_clean$OS.time <- os_time_result$value

numeric_conversion_audit <- rbind(
  age_result$invalid,
  os_time_result$invalid
)

write.table(
  numeric_conversion_audit,
  file = file.path(
    result_dir,
    "Workflow-04-1E_numeric_conversion_invalid_values.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

age_nonfinite_input <- age_result$nonfinite

os_time_nonfinite_input <- os_time_result$nonfinite

############################################################
# 15. Standardize recognized OS.status labels
#
# Only Alive and Dead are standardized.
# Other values are preserved and reported.
############################################################

clinical_clean$OS.status <- trimws(
  as.character(
    clinical_clean$OS.status
  )
)

os_status_lower <- tolower(
  clinical_clean$OS.status
)

clinical_clean$OS.status[
  os_status_lower == "alive"
] <- "Alive"

clinical_clean$OS.status[
  os_status_lower == "dead"
] <- "Dead"

recognized_os_status <- c(
  "Alive",
  "Dead"
)

unknown_os_status <- (
  !is.na(clinical_clean$OS.status) &
    !clinical_clean$OS.status %in%
    recognized_os_status
)

unknown_os_status_values <- unique(
  clinical_clean$OS.status[
    unknown_os_status
  ]
)

write.table(
  data.frame(
    Unknown_OS_status = unknown_os_status_values,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    result_dir,
    "Workflow-04-1E_unrecognized_OS_status_values.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 16. Age and OS.time numerical QC
############################################################

finite_age <- clinical_clean$Age[
  is.finite(clinical_clean$Age)
]

finite_os_time <- clinical_clean$OS.time[
  is.finite(clinical_clean$OS.time)
]

safe_summary <- function(x, statistic) {
  
  x <- x[
    is.finite(x)
  ]
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  switch(
    statistic,
    min = min(x),
    median = median(x),
    max = max(x),
    stop("Unsupported statistic.")
  )
}

numeric_qc_summary <- data.frame(
  Variable = c(
    "Age",
    "OS.time"
  ),
  N = c(
    length(clinical_clean$Age),
    length(clinical_clean$OS.time)
  ),
  Missing_N = c(
    sum(is.na(clinical_clean$Age)),
    sum(is.na(clinical_clean$OS.time))
  ),
  Nonfinite_N = c(
    sum(age_nonfinite_input),
    sum(os_time_nonfinite_input)
  ),
  Minimum = c(
    safe_summary(finite_age, "min"),
    safe_summary(finite_os_time, "min")
  ),
  Median = c(
    safe_summary(finite_age, "median"),
    safe_summary(finite_os_time, "median")
  ),
  Maximum = c(
    safe_summary(finite_age, "max"),
    safe_summary(finite_os_time, "max")
  ),
  stringsAsFactors = FALSE
)

write.table(
  numeric_qc_summary,
  file = file.path(
    result_dir,
    "Workflow-04-1E_numeric_variable_QC.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 17. Identify samples with missing Age, OS.time, or Stage
############################################################

age_missing <- is_missing_value(
  clinical_clean$Age
)

os_time_missing <- is_missing_value(
  clinical_clean$OS.time
)

stage2_missing <- is_missing_value(
  clinical_clean$pathologic_stage_2
)

write.table(
  data.frame(
    Sample = clinical_clean$Sample[
      age_missing
    ]
  ),
  file = file.path(
    result_dir,
    "Workflow-04-1E_samples_missing_Age.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    Sample = clinical_clean$Sample[
      os_time_missing
    ]
  ),
  file = file.path(
    result_dir,
    "Workflow-04-1E_samples_missing_OS_time.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    Sample = clinical_clean$Sample[
      stage2_missing
    ]
  ),
  file = file.path(
    result_dir,
    "Workflow-04-1E_samples_missing_pathologic_stage_2.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 18. Clinical variable distributions after cleaning
############################################################

clean_distribution_vars <- c(
  "Gender",
  "OS.status",
  "pathologic_stage",
  "pathologic_stage_2",
  "pathologic_t",
  "pathologic_n",
  "pathologic_m"
)

for (variable_name in clean_distribution_vars) {
  
  cat("\n========================================\n")
  cat(
    variable_name,
    " distribution after cleaning\n"
  )
  cat("========================================\n\n")
  
  print(
    table(
      clinical_clean[[variable_name]],
      useNA = "always"
    )
  )
}


############################################################
# 19. Final clinical sample-ID and cohort integrity checks
############################################################

if (anyNA(clinical_clean$Sample)) {
  stop(
    "Missing Sample IDs were detected in the cleaned ",
    "Primary Tumor clinical cohort."
  )
}

if (anyDuplicated(clinical_clean$Sample)) {
  stop(
    "Duplicated Sample IDs were detected in the cleaned ",
    "Primary Tumor clinical cohort."
  )
}

if (
  nrow(clinical_clean) !=
  expected_primary_tumor_n
) {
  stop(
    "The cleaned clinical cohort does not contain ",
    "the expected 513 samples."
  )
}

if (
  !identical(
    clinical_clean$Sample,
    primary_tumor_ids
  )
) {
  stop(
    "The cleaned clinical sample order differs from ",
    "the Primary Tumor expression sample order."
  )
}

if (
  any(
    substr(
      clinical_clean$Sample,
      14,
      15
    ) != "01"
  )
) {
  stop(
    "Non-Primary-Tumor samples were detected in ",
    "the cleaned clinical cohort."
  )
}


############################################################
# 20. Overall-survival consistency QC
#
# No sample is removed in this section.
############################################################

os_time <- clinical_clean$OS.time
os_status <- clinical_clean$OS.status

os_time_missing <- is.na(os_time)

os_time_nonfinite <- os_time_nonfinite_input

os_time_negative <- (
  is.finite(os_time) &
    os_time < 0
)

os_time_zero <- (
  is.finite(os_time) &
    os_time == 0
)

os_status_missing <- is.na(
  os_status
)

os_status_unrecognized <- (
  !os_status_missing &
    !os_status %in% recognized_os_status
)

os_valid_time <- (
  is.finite(os_time) &
    os_time > 0
)

os_valid_status <- (
  !os_status_missing &
    os_status %in% recognized_os_status
)

os_valid_complete <- (
  os_valid_time &
    os_valid_status
)

os_qc_summary <- data.frame(
  Metric = c(
    "Total_Primary_Tumor_samples",
    "Missing_OS_time",
    "Nonfinite_OS_time",
    "Negative_OS_time",
    "Zero_OS_time",
    "Missing_OS_status",
    "Unrecognized_OS_status",
    "Valid_positive_OS_time",
    "Valid_OS_status",
    "Valid_complete_OS_records"
  ),
  N = c(
    nrow(clinical_clean),
    sum(os_time_missing),
    sum(os_time_nonfinite),
    sum(os_time_negative),
    sum(os_time_zero),
    sum(os_status_missing),
    sum(os_status_unrecognized),
    sum(os_valid_time),
    sum(os_valid_status),
    sum(os_valid_complete)
  ),
  stringsAsFactors = FALSE
)

write.table(
  os_qc_summary,
  file = file.path(
    result_dir,
    "Workflow-04-1F-1_OS_consistency_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 21. Save sample-level OS QC records
############################################################

os_qc_records <- data.frame(
  Sample = clinical_clean$Sample,
  patient_id = clinical_clean$patient_id,
  OS.time = clinical_clean$OS.time,
  OS.status = clinical_clean$OS.status,
  Age = clinical_clean$Age,
  pathologic_stage_2 = clinical_clean$pathologic_stage_2,
  pathologic_stage = clinical_clean$pathologic_stage,
  OS_time_missing = os_time_missing,
  OS_time_nonfinite = os_time_nonfinite,
  OS_time_negative = os_time_negative,
  OS_time_zero = os_time_zero,
  OS_status_missing = os_status_missing,
  OS_status_unrecognized = os_status_unrecognized,
  OS_valid_complete = os_valid_complete,
  stringsAsFactors = FALSE
)

write.table(
  os_qc_records,
  file = file.path(
    result_dir,
    "Workflow-04-1F-1_OS_sample_level_QC.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 22. Save cleaned clinical table
#
# This is the primary clinical input for Workflow-04-2A.
############################################################

write.table(
  clinical_clean,
  file = file.path(
    result_dir,
    "Workflow-04-1E_clinical_clean.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)


############################################################
# 23. Save Primary Tumor sample manifest
#
# The order follows the original expression matrix columns.
############################################################

primary_tumor_sample_manifest <- data.frame(
  Sample = primary_tumor_ids,
  Patient_ID = substr(
    primary_tumor_ids,
    1,
    12
  ),
  Sample_type_code = substr(
    primary_tumor_ids,
    14,
    15
  ),
  stringsAsFactors = FALSE
)

write.table(
  primary_tumor_sample_manifest,
  file = file.path(
    result_dir,
    "Workflow-04-1F-1_primary_tumor_sample_list.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 24. Final workflow summary
############################################################

workflow_summary <- data.frame(
  Item = c(
    "Expression_input",
    "Clinical_input",
    "Total_expression_samples",
    "Primary_Tumor_samples",
    "Expected_Primary_Tumor_samples",
    "Unique_Primary_Tumor_patients",
    "Primary_Tumor_samples_with_clinical_match",
    "Cleaned_clinical_rows",
    "Valid_complete_OS_records",
    "Random_seed"
  ),
  Value = c(
    basename(expression_file),
    basename(clinical_file),
    length(sample_ids),
    length(primary_tumor_ids),
    expected_primary_tumor_n,
    length(unique(primary_tumor_patient_ids)),
    sum(!is.na(primary_match_index)),
    nrow(clinical_clean),
    sum(os_valid_complete),
    "Not applicable"
  ),
  stringsAsFactors = FALSE
)

write.table(
  workflow_summary,
  file = file.path(
    result_dir,
    "Workflow-04-1_workflow_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat("WORKFLOW-04-1 COMPLETED\n")
cat("========================================\n\n")

print(
  workflow_summary,
  row.names = FALSE
)

cat(
  "\nOutput directory:\n",
  result_dir,
  "\n"
)

cat(
  "\nPrimary downstream clinical input:\n",
  file.path(
    result_dir,
    "Workflow-04-1E_clinical_clean.txt"
  ),
  "\n"
)

cat(
  "\nPrimary Tumor sample manifest:\n",
  file.path(
    result_dir,
    "Workflow-04-1F-1_primary_tumor_sample_list.txt"
  ),
  "\n"
)

############################################################
# End of Workflow-04-1
############################################################

############################################################
# Paper 1: LUAD transcriptomic analysis
#
# Script:
#   GitHub2_Paper1_tumor_only_WGCNA_1.R
#
# Workflow:
#   Workflow-04-2A
#   Primary tumor expression matrix construction
#
# Research purpose:
#   1. Read the original TCGA LUAD FPKM-UQ expression matrix.
#   2. Identify Primary Tumor samples using TCGA sample type
#      code 01.
#   3. Extract the 513 Primary Tumor expression samples.
#   4. Verify exact sample-level consistency with the clinical
#      cohort reconstructed in Workflow-04-1.
#   5. Align the clinical table to the expression matrix
#      sample-column order.
#   6. Export the tumor-only expression matrix and its
#      corresponding clinical and sample manifests.
#
# Input data:
#
#   1. data/Merge_RNA_seq_FPKM_UQ.txt
#      Original TCGA LUAD FPKM-UQ expression matrix.
#
#   2. results/02_tumor_only_WGCNA/Workflow-04-1/
#      Workflow-04-1E_clinical_clean.txt
#      Cleaned Primary Tumor clinical table.
#
#   3. results/02_tumor_only_WGCNA/Workflow-04-1/
#      Workflow-04-1F-1_primary_tumor_sample_list.txt
#      Primary Tumor sample manifest.
#
# Sample inclusion criteria:
#
#   - Standard TCGA sample-level barcode.
#   - TCGA sample type code = 01 (Primary Tumor).
#   - Exact sample-level match with the Workflow-04-1
#     clinical cohort.
#   - One Primary Tumor sample per patient.
#
# Expected cohort:
#
#   513 Primary Tumor samples.
#
# Important:
#
#   - No gene filtering is performed.
#   - No Ensembl-to-Gene-Symbol conversion is performed.
#   - No log2 transformation is performed.
#   - No expression normalization is performed.
#   - No WGCNA is performed.
#   - Original expression values are retained.
#   - Original input files are not modified.
#
# Software environment:
#
#   R 4.5.1 (project target version).
#   Package: here (version recorded at runtime).
#
# Random seed:
#
#   Not applicable. This workflow uses deterministic
#   data extraction, matching, and file writing.
#
# Output directory:
#
#   results/02_tumor_only_WGCNA/Workflow-04-2A/
#
# Downstream workflow:
#
#   Workflow-04-2B reads:
#
#   Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt
#
#   The expression matrix contains:
#
#   - Rows: original gene identifiers.
#   - Columns: 513 Primary Tumor samples.
#   - First column: Tag (gene identifier).
#
#   The expression matrix remains in the original FPKM-UQ
#   scale at the end of this workflow.
############################################################


############################################################
# 00. Project configuration and setup
############################################################

required_packages <- c(
  "here"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "The following required packages are missing: ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\nInstall them before running this analysis script."
  )
}


############################################################
# 00-1. Project directories
############################################################

project_root <- here::here()

data_dir <- file.path(
  project_root,
  "data"
)

workflow04_1_dir <- file.path(
  project_root,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-1"
)

result_dir <- file.path(
  project_root,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2A"
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 00-2. Input file definitions
############################################################

expression_file <- file.path(
  data_dir,
  "Merge_RNA_seq_FPKM_UQ.txt"
)

clinical_file <- file.path(
  workflow04_1_dir,
  "Workflow-04-1E_clinical_clean.txt"
)

sample_manifest_file <- file.path(
  workflow04_1_dir,
  "Workflow-04-1F-1_primary_tumor_sample_list.txt"
)


############################################################
# 00-3. Expected cohort parameters
############################################################

expected_primary_tumor_n <- 513L

tcga_barcode_pattern <- (
  "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}$"
)


############################################################
# 00-4. Input file checks
############################################################

required_input_files <- c(
  expression_file,
  clinical_file,
  sample_manifest_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0L) {
  stop(
    "The following required input files were not found:\n",
    paste(
      missing_input_files,
      collapse = "\n"
    ),
    "\nPlease complete Workflow-04-1 before running Workflow-04-2A."
  )
}


############################################################
# 00-5. Software environment record
############################################################

r_version <- as.character(
  getRversion()
)

package_versions <- data.frame(
  Package = required_packages,
  Version = vapply(
    required_packages,
    function(pkg) {
      as.character(
        utils::packageVersion(pkg)
      )
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)

cat("\n========================================\n")
cat("WORKFLOW-04-2A SOFTWARE ENVIRONMENT\n")
cat("========================================\n\n")

cat(
  "R version:",
  r_version,
  "\n"
)

print(
  package_versions,
  row.names = FALSE
)

if (r_version != "4.5.1") {
  warning(
    "The project target R version is 4.5.1. ",
    "Current version: ",
    r_version
  )
}

write.table(
  package_versions,
  file = file.path(
    result_dir,
    "Workflow-04-2A_software_versions.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 01. Read original expression matrix
############################################################

cat("\n========================================\n")
cat("READING ORIGINAL EXPRESSION MATRIX\n")
cat("========================================\n\n")

cat(
  "Input file:\n",
  expression_file,
  "\n\n",
  sep = ""
)

expr_raw <- read.delim(
  expression_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

cat(
  "Original gene rows:",
  nrow(expr_raw),
  "\n"
)

cat(
  "Original total columns:",
  ncol(expr_raw),
  "\n"
)

cat(
  "Original sample columns:",
  ncol(expr_raw) - 1L,
  "\n"
)


############################################################
# 02. Basic expression matrix structure checks
############################################################

if (ncol(expr_raw) < 2L) {
  stop(
    "The expression matrix must contain one gene-ID column ",
    "and at least one sample column."
  )
}

if (nrow(expr_raw) == 0L) {
  stop(
    "The expression matrix contains no gene rows."
  )
}

expression_column_names <- colnames(
  expr_raw
)

if (
  anyNA(expression_column_names) ||
  any(trimws(expression_column_names) == "")
) {
  stop(
    "Missing or blank column names were detected ",
    "in the expression matrix."
  )
}

expr_sample_ids <- expression_column_names[-1]


############################################################
# 03. Check expression sample IDs
############################################################

cat("\n========================================\n")
cat("EXPRESSION SAMPLE ID QC\n")
cat("========================================\n\n")

if (anyDuplicated(expr_sample_ids)) {
  
  duplicated_expression_ids <- unique(
    expr_sample_ids[
      duplicated(expr_sample_ids)
    ]
  )
  
  write.table(
    data.frame(
      Sample = duplicated_expression_ids
    ),
    file = file.path(
      result_dir,
      "Workflow-04-2A_duplicate_expression_sample_IDs.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  stop(
    "Duplicated expression sample IDs were detected."
  )
}

barcode_valid <- grepl(
  tcga_barcode_pattern,
  expr_sample_ids
)

invalid_expression_ids <- expr_sample_ids[
  !barcode_valid
]

if (length(invalid_expression_ids) > 0L) {
  
  write.table(
    data.frame(
      Sample = invalid_expression_ids
    ),
    file = file.path(
      result_dir,
      "Workflow-04-2A_invalid_expression_sample_IDs.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  stop(
    "Non-standard TCGA sample barcodes were detected."
  )
}

cat(
  "Total expression sample IDs:",
  length(expr_sample_ids),
  "\n"
)

cat(
  "Valid TCGA sample IDs:",
  sum(barcode_valid),
  "\n"
)


############################################################
# 04. Determine TCGA sample type
############################################################

sample_type_code <- substr(
  expr_sample_ids,
  14,
  15
)

sample_type_table <- as.data.frame(
  table(
    sample_type_code,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

colnames(sample_type_table) <- c(
  "Sample_type_code",
  "N"
)

cat("\n========================================\n")
cat("TCGA SAMPLE-TYPE COMPOSITION\n")
cat("========================================\n\n")

print(
  sample_type_table,
  row.names = FALSE
)

write.table(
  sample_type_table,
  file = file.path(
    result_dir,
    "Workflow-04-2A_expression_sample_composition.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 05. Identify Primary Tumor samples
#
# TCGA sample type code:
#   01 = Primary Tumor
############################################################

primary_tumor_ids <- expr_sample_ids[
  sample_type_code == "01"
]

if (length(primary_tumor_ids) == 0L) {
  stop(
    "No Primary Tumor samples (sample type 01) were found."
  )
}

cat("\n========================================\n")
cat("PRIMARY TUMOR SAMPLE EXTRACTION\n")
cat("========================================\n\n")

cat(
  "Primary Tumor samples:",
  length(primary_tumor_ids),
  "\n"
)

cat(
  "Unique Primary Tumor samples:",
  length(unique(primary_tumor_ids)),
  "\n"
)

if (
  length(primary_tumor_ids) !=
  expected_primary_tumor_n
) {
  stop(
    "Unexpected Primary Tumor sample count.\n",
    "Expected: ",
    expected_primary_tumor_n,
    "\nObserved: ",
    length(primary_tumor_ids)
  )
}

if (anyDuplicated(primary_tumor_ids)) {
  stop(
    "Duplicated Primary Tumor sample IDs were detected."
  )
}


############################################################
# 06. Check patient-level duplication
############################################################

primary_tumor_patient_ids <- substr(
  primary_tumor_ids,
  1,
  12
)

patient_frequency <- table(
  primary_tumor_patient_ids
)

multiple_primary_tumor_patients <- names(
  patient_frequency[
    patient_frequency > 1
  ]
)

patient_duplication_audit <- data.frame(
  Patient_ID = names(patient_frequency),
  N_Primary_Tumor_Samples = as.integer(
    patient_frequency
  ),
  stringsAsFactors = FALSE
)

write.table(
  patient_duplication_audit,
  file = file.path(
    result_dir,
    "Workflow-04-2A_patient_sample_count.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "Unique Primary Tumor patients:",
  length(unique(primary_tumor_patient_ids)),
  "\n"
)

cat(
  "Patients with multiple Primary Tumor samples:",
  length(multiple_primary_tumor_patients),
  "\n"
)

if (length(multiple_primary_tumor_patients) > 0L) {
  stop(
    "Multiple Primary Tumor samples were detected ",
    "for at least one patient. Please investigate."
  )
}


############################################################
# 07. Extract Primary Tumor expression columns
#
# match() preserves the original expression-matrix order.
############################################################

primary_tumor_column_indices <- match(
  primary_tumor_ids,
  expr_sample_ids
)

if (anyNA(primary_tumor_column_indices)) {
  stop(
    "At least one Primary Tumor sample could not be ",
    "matched to the original expression matrix."
  )
}

if (
  length(primary_tumor_column_indices) !=
  expected_primary_tumor_n
) {
  stop(
    "The number of extracted expression columns ",
    "does not equal 513."
  )
}

expr_primary_tumor <- expr_raw[
  ,
  c(
    1,
    primary_tumor_column_indices + 1L
  ),
  drop = FALSE
]


############################################################
# 08. Standardize the gene identifier column name
#
# The original gene identifiers are retained without
# Ensembl-to-Symbol conversion.
#
# The standardized first-column name "Tag" is used by
# downstream Workflow-04-2B.
############################################################

gene_id_column_original <- colnames(
  expr_primary_tumor
)[1]

colnames(expr_primary_tumor)[1] <- "Tag"


############################################################
# 09. Gene identifier integrity checks
############################################################

gene_ids <- as.character(
  expr_primary_tumor$Tag
)

gene_ids <- trimws(
  gene_ids
)

gene_id_missing <- (
  is.na(gene_ids) |
    gene_ids == ""
)

if (any(gene_id_missing)) {
  
  write.table(
    data.frame(
      Original_row = which(gene_id_missing),
      Gene_ID = gene_ids[gene_id_missing],
      stringsAsFactors = FALSE
    ),
    file = file.path(
      result_dir,
      "Workflow-04-2A_missing_gene_IDs.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  stop(
    "Missing or blank gene identifiers were detected."
  )
}

duplicated_gene_ids <- unique(
  gene_ids[
    duplicated(gene_ids)
  ]
)

if (length(duplicated_gene_ids) > 0L) {
  
  write.table(
    data.frame(
      Gene_ID = duplicated_gene_ids
    ),
    file = file.path(
      result_dir,
      "Workflow-04-2A_duplicate_gene_IDs.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  stop(
    "Duplicated gene identifiers were detected. ",
    "Please inspect the original expression matrix."
  )
}

expr_primary_tumor$Tag <- gene_ids


############################################################
# 10. Check expression matrix numeric structure
#
# No expression transformation is performed here.
# Original FPKM-UQ values are retained.
############################################################

expression_sample_columns <- (
  colnames(expr_primary_tumor)[-1]
)

non_numeric_columns <- expression_sample_columns[
  !vapply(
    expr_primary_tumor[
      ,
      expression_sample_columns,
      drop = FALSE
    ],
    is.numeric,
    logical(1)
  )
]

if (length(non_numeric_columns) > 0L) {
  
  stop(
    "Non-numeric expression columns were detected:\n",
    paste(
      non_numeric_columns,
      collapse = ", "
    ),
    "\nPlease inspect the original expression input."
  )
}

expression_values <- as.matrix(
  expr_primary_tumor[
    ,
    expression_sample_columns,
    drop = FALSE
  ]
)

storage.mode(expression_values) <- "double"

if (nrow(expression_values) != length(gene_ids)) {
  stop(
    "Expression matrix row count does not match gene IDs."
  )
}

if (ncol(expression_values) != expected_primary_tumor_n) {
  stop(
    "Expression matrix does not contain exactly 513 ",
    "Primary Tumor sample columns."
  )
}

if (!identical(
  colnames(expression_values),
  primary_tumor_ids
)) {
  stop(
    "The extracted expression matrix sample order ",
    "differs from the original Primary Tumor sample order."
  )
}


############################################################
# 11. Expression matrix structure summary
############################################################

expression_matrix_summary <- data.frame(
  Metric = c(
    "Original_gene_rows",
    "Original_sample_columns",
    "Primary_Tumor_gene_rows",
    "Primary_Tumor_sample_columns",
    "Unique_gene_IDs",
    "Missing_expression_values",
    "Nonfinite_expression_values",
    "Negative_expression_values"
  ),
  Value = c(
    nrow(expr_raw),
    ncol(expr_raw) - 1L,
    nrow(expr_primary_tumor),
    ncol(expr_primary_tumor) - 1L,
    length(unique(gene_ids)),
    sum(is.na(expression_values)),
    sum(!is.na(expression_values) &
          !is.finite(expression_values)),
    sum(expression_values < 0, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

cat("\n========================================\n")
cat("PRIMARY TUMOR EXPRESSION MATRIX SUMMARY\n")
cat("========================================\n\n")

print(
  expression_matrix_summary,
  row.names = FALSE
)

write.table(
  expression_matrix_summary,
  file = file.path(
    result_dir,
    "Workflow-04-2A_expression_matrix_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 12. Read Workflow-04-1 clinical output
############################################################

cat("\n========================================\n")
cat("READING WORKFLOW-04-1 CLINICAL OUTPUT\n")
cat("========================================\n\n")

clinical_clean <- read.delim(
  clinical_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

sample_manifest <- read.delim(
  sample_manifest_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)


############################################################
# 13. Check required clinical columns
############################################################

if (!"Sample" %in% colnames(clinical_clean)) {
  stop(
    "The Workflow-04-1 clinical table does not contain ",
    "the required Sample column."
  )
}

if (!"tumor_normal" %in% colnames(clinical_clean)) {
  stop(
    "The Workflow-04-1 clinical table does not contain ",
    "the required tumor_normal column."
  )
}

if (!"Sample" %in% colnames(sample_manifest)) {
  stop(
    "The Workflow-04-1 sample manifest does not contain ",
    "the required Sample column."
  )
}


############################################################
# 14. Validate clinical sample IDs
############################################################

clinical_sample_ids <- trimws(
  as.character(
    clinical_clean$Sample
  )
)

clinical_sample_missing <- (
  is.na(clinical_sample_ids) |
    clinical_sample_ids == ""
)

if (any(clinical_sample_missing)) {
  stop(
    "Missing or blank Sample IDs were detected ",
    "in the Workflow-04-1 clinical table."
  )
}

if (anyDuplicated(clinical_sample_ids)) {
  
  duplicated_clinical_ids <- unique(
    clinical_sample_ids[
      duplicated(clinical_sample_ids)
    ]
  )
  
  write.table(
    data.frame(
      Sample = duplicated_clinical_ids
    ),
    file = file.path(
      result_dir,
      "Workflow-04-2A_duplicate_clinical_sample_IDs.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  stop(
    "Duplicated clinical sample IDs were detected."
  )
}

if (
  length(clinical_sample_ids) !=
  expected_primary_tumor_n
) {
  stop(
    "The Workflow-04-1 clinical table does not contain ",
    "exactly 513 samples."
  )
}


############################################################
# 15. Validate Workflow-04-1 sample manifest
############################################################

manifest_sample_ids <- trimws(
  as.character(
    sample_manifest$Sample
  )
)

if (
  anyNA(manifest_sample_ids) ||
  any(manifest_sample_ids == "")
) {
  stop(
    "Missing or blank Sample IDs were detected ",
    "in the Workflow-04-1 sample manifest."
  )
}

if (anyDuplicated(manifest_sample_ids)) {
  stop(
    "Duplicated Sample IDs were detected ",
    "in the Workflow-04-1 sample manifest."
  )
}

if (
  length(manifest_sample_ids) !=
  expected_primary_tumor_n
) {
  stop(
    "The Workflow-04-1 sample manifest does not ",
    "contain exactly 513 samples."
  )
}


############################################################
# 16. Bidirectional exact sample matching
############################################################

cat("\n========================================\n")
cat("BIDIRECTIONAL SAMPLE MATCHING\n")
cat("========================================\n\n")

# Expression Primary Tumor -> Clinical
expr_to_clinical_index <- match(
  primary_tumor_ids,
  clinical_sample_ids
)

expr_not_in_clinical <- primary_tumor_ids[
  is.na(expr_to_clinical_index)
]

# Clinical -> Expression Primary Tumor
clinical_to_expr_index <- match(
  clinical_sample_ids,
  primary_tumor_ids
)

clinical_not_in_expr <- clinical_sample_ids[
  is.na(clinical_to_expr_index)
]

# Manifest -> Expression Primary Tumor
manifest_to_expr_index <- match(
  manifest_sample_ids,
  primary_tumor_ids
)

manifest_not_in_expr <- manifest_sample_ids[
  is.na(manifest_to_expr_index)
]

matching_summary <- data.frame(
  Metric = c(
    "Expression_Primary_Tumor_samples",
    "Clinical_samples",
    "Expression_samples_not_in_clinical",
    "Clinical_samples_not_in_expression",
    "Workflow04_1_manifest_samples",
    "Manifest_samples_not_in_expression"
  ),
  Value = c(
    length(primary_tumor_ids),
    length(clinical_sample_ids),
    length(expr_not_in_clinical),
    length(clinical_not_in_expr),
    length(manifest_sample_ids),
    length(manifest_not_in_expr)
  ),
  stringsAsFactors = FALSE
)

print(
  matching_summary,
  row.names = FALSE
)

write.table(
  matching_summary,
  file = file.path(
    result_dir,
    "Workflow-04-2A_bidirectional_matching_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 17. Save unmatched sample audits
############################################################

write.table(
  data.frame(
    Sample = expr_not_in_clinical
  ),
  file = file.path(
    result_dir,
    "Workflow-04-2A_expression_not_in_clinical.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    Sample = clinical_not_in_expr
  ),
  file = file.path(
    result_dir,
    "Workflow-04-2A_clinical_not_in_expression.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    Sample = manifest_not_in_expr
  ),
  file = file.path(
    result_dir,
    "Workflow-04-2A_manifest_not_in_expression.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 18. Require exact sample-set consistency
############################################################

if (
  length(expr_not_in_clinical) > 0L ||
  length(clinical_not_in_expr) > 0L ||
  length(manifest_not_in_expr) > 0L
) {
  stop(
    "Bidirectional matching between the expression matrix, ",
    "clinical table, and Workflow-04-1 manifest is not exact."
  )
}

expression_clinical_sets_identical <- identical(
  sort(primary_tumor_ids),
  sort(clinical_sample_ids)
)

expression_manifest_sets_identical <- identical(
  sort(primary_tumor_ids),
  sort(manifest_sample_ids)
)

if (!expression_clinical_sets_identical) {
  stop(
    "Expression and clinical sample sets are not identical."
  )
}

if (!expression_manifest_sets_identical) {
  stop(
    "Expression and Workflow-04-1 manifest sample sets ",
    "are not identical."
  )
}


############################################################
# 19. Verify the original sample order against Workflow-04-1
############################################################

manifest_order_identical <- identical(
  primary_tumor_ids,
  manifest_sample_ids
)

cat(
  "Expression and Workflow-04-1 manifest order identical:",
  manifest_order_identical,
  "\n"
)

if (!manifest_order_identical) {
  warning(
    "The sample sets are identical, but their orders differ. ",
    "The expression matrix order will be retained as the ",
    "reference order for downstream analysis."
  )
}


############################################################
# 20. Align clinical data to expression sample order
############################################################

clinical_match_index <- match(
  primary_tumor_ids,
  clinical_sample_ids
)

if (anyNA(clinical_match_index)) {
  stop(
    "NA matching indices were detected when aligning ",
    "the clinical table."
  )
}

clinical_clean_aligned <- clinical_clean[
  clinical_match_index,
  ,
  drop = FALSE
]

rownames(
  clinical_clean_aligned
) <- NULL


############################################################
# 21. Verify expression-clinical sample alignment
############################################################

clinical_aligned_sample_ids <- as.character(
  clinical_clean_aligned$Sample
)

alignment_confirmed <- identical(
  primary_tumor_ids,
  clinical_aligned_sample_ids
)

cat("\n========================================\n")
cat("EXPRESSION-CLINICAL ALIGNMENT\n")
cat("========================================\n\n")

cat(
  "Expression Primary Tumor samples:",
  length(primary_tumor_ids),
  "\n"
)

cat(
  "Aligned clinical rows:",
  nrow(clinical_clean_aligned),
  "\n"
)

cat(
  "Sample order aligned:",
  alignment_confirmed,
  "\n"
)

if (!alignment_confirmed) {
  stop(
    "Expression and clinical sample order is not aligned."
  )
}


############################################################
# 22. Verify Primary Tumor classification
############################################################

clinical_clean_aligned$tumor_normal <- tolower(
  trimws(
    as.character(
      clinical_clean_aligned$tumor_normal
    )
  )
)

tumor_classification_table <- as.data.frame(
  table(
    clinical_clean_aligned$tumor_normal,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

colnames(tumor_classification_table) <- c(
  "tumor_normal",
  "N"
)

write.table(
  tumor_classification_table,
  file = file.path(
    result_dir,
    "Workflow-04-2A_tumor_classification_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (
  anyNA(clinical_clean_aligned$tumor_normal) ||
  any(
    clinical_clean_aligned$tumor_normal != "tumor"
  )
) {
  stop(
    "The aligned clinical table contains missing or ",
    "non-tumor classifications."
  )
}


############################################################
# 23. Final sample-order and matrix checks
############################################################

if (
  nrow(expr_primary_tumor) !=
  nrow(expression_values)
) {
  stop(
    "Unexpected expression matrix row count."
  )
}

if (
  ncol(expr_primary_tumor) - 1L !=
  expected_primary_tumor_n
) {
  stop(
    "The final expression matrix does not contain 513 samples."
  )
}

if (
  nrow(clinical_clean_aligned) !=
  expected_primary_tumor_n
) {
  stop(
    "The final aligned clinical table does not contain 513 rows."
  )
}

if (
  !identical(
    colnames(expr_primary_tumor)[-1],
    clinical_clean_aligned$Sample
  )
) {
  stop(
    "Final expression-column and clinical-row order check failed."
  )
}


############################################################
# 24. Final cohort summary
############################################################

workflow_summary <- data.frame(
  Item = c(
    "Expression_input_file",
    "Clinical_input_file",
    "Workflow04_1_sample_manifest",
    "Original_gene_rows",
    "Original_expression_sample_columns",
    "Primary_Tumor_samples_extracted",
    "Unique_Primary_Tumor_patients",
    "Clinical_samples_matched",
    "Expression_clinical_sample_sets_identical",
    "Expression_manifest_sample_sets_identical",
    "Expression_clinical_order_aligned",
    "Tumor_classification",
    "Gene_filtering_performed",
    "Log2_transformation_performed",
    "Normalization_performed",
    "WGCNA_performed",
    "Random_seed"
  ),
  Value = c(
    basename(expression_file),
    basename(clinical_file),
    basename(sample_manifest_file),
    nrow(expr_raw),
    ncol(expr_raw) - 1L,
    ncol(expr_primary_tumor) - 1L,
    length(unique(primary_tumor_patient_ids)),
    nrow(clinical_clean_aligned),
    expression_clinical_sets_identical,
    expression_manifest_sets_identical,
    alignment_confirmed,
    "tumor",
    "No",
    "No",
    "No",
    "No",
    "Not applicable"
  ),
  stringsAsFactors = FALSE
)

write.table(
  workflow_summary,
  file = file.path(
    result_dir,
    "Workflow-04-2A_workflow_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 25. Save Primary Tumor expression matrix
############################################################

expression_output_file <- file.path(
  result_dir,
  "Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt"
)

write.table(
  expr_primary_tumor,
  file = expression_output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)


############################################################
# 26. Save aligned clinical table
############################################################

clinical_output_file <- file.path(
  result_dir,
  "Workflow-04-2A_clinical_aligned.txt"
)

write.table(
  clinical_clean_aligned,
  file = clinical_output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)


############################################################
# 27. Save Primary Tumor sample list
############################################################

sample_list_output_file <- file.path(
  result_dir,
  "Workflow-04-2A_primary_tumor_sample_list.txt"
)

primary_tumor_sample_list <- data.frame(
  Sample = primary_tumor_ids,
  Patient_ID = primary_tumor_patient_ids,
  Sample_type_code = substr(
    primary_tumor_ids,
    14,
    15
  ),
  Expression_column_order = seq_along(
    primary_tumor_ids
  ),
  stringsAsFactors = FALSE
)

write.table(
  primary_tumor_sample_list,
  file = sample_list_output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 28. Verify saved output files
############################################################

output_files <- c(
  expression_output_file,
  clinical_output_file,
  sample_list_output_file,
  file.path(
    result_dir,
    "Workflow-04-2A_workflow_summary.txt"
  )
)

output_file_audit <- data.frame(
  File = basename(output_files),
  Exists = file.exists(output_files),
  File_size_bytes = ifelse(
    file.exists(output_files),
    file.info(output_files)$size,
    NA_real_
  ),
  stringsAsFactors = FALSE
)

write.table(
  output_file_audit,
  file = file.path(
    result_dir,
    "Workflow-04-2A_output_file_audit.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (!all(output_file_audit$Exists)) {
  stop(
    "At least one required output file was not generated."
  )
}


############################################################
# 29. Final workflow message
############################################################

cat("\n========================================\n")
cat("WORKFLOW-04-2A COMPLETED\n")
cat("========================================\n\n")

cat(
  "Original expression matrix:\n",
  nrow(expr_raw),
  " genes × ",
  ncol(expr_raw) - 1L,
  " samples\n\n",
  sep = ""
)

cat(
  "Primary Tumor expression matrix:\n",
  nrow(expr_primary_tumor),
  " genes × ",
  ncol(expr_primary_tumor) - 1L,
  " samples\n\n",
  sep = ""
)

cat(
  "Aligned clinical samples:",
  nrow(clinical_clean_aligned),
  "\n"
)

cat(
  "Expression-clinical alignment:",
  alignment_confirmed,
  "\n\n"
)

cat(
  "Output directory:\n",
  result_dir,
  "\n\n"
)

cat("Principal output files:\n")

print(
  output_files
)

cat(
  "\nNo expression filtering, log2 transformation, ",
  "normalization, or WGCNA was performed.\n",
  sep = ""
)

############################################################
# End of Workflow-04-2A
############################################################

############################################################
# Workflow-04-2B
# FPKM-UQ Expression Data Structure and Distribution Audit
#
# Project:
#   Paper 1 - LUAD transcriptomic analysis
#
# Code path:
#   scripts/02_tumor_only_WGCNA/
#   GitHub2_Paper1_tumor_only_WGCNA_1.R
#
# Research purpose:
#   To audit the structure, numeric integrity, missing/invalid
#   values, and expression distributions of the primary-tumor
#   FPKM-UQ matrix generated in Workflow-04-2A.
#
# Input data:
#   1. Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt
#   2. Workflow-04-2A_clinical_aligned.txt
#   3. Workflow-04-2A_primary_tumor_sample_list.txt
#
# Input directory:
#   results/02_tumor_only_WGCNA/Workflow-04-2A/
#
# Expected cohort:
#   513 TCGA primary LUAD tumor samples (sample type code 01).
#   The historical source matrix contained 60,483 gene IDs.
#   The gene count is reported for auditing, not hard-coded
#   as a mandatory inclusion criterion.
#
# Sample inclusion:
#   All samples retained in Workflow-04-2A.
#   Exact sample-level TCGA barcodes are required.
#   Expression, clinical data, and sample manifest must agree
#   in both sample identity and ordering.
#
# Analysis scope:
#   1. Expression matrix structure and identifier audit
#   2. Numeric data type and invalid-token audit
#   3. NA / NaN / Inf / negative-value assessment
#   4. Global expression distribution
#   5. Sample-level expression distribution
#   6. Gene-level expression distribution
#   7. Zero-expression and low-expression summaries
#   8. Descriptive identification of distribution tails
#   9. QC plots and machine-readable audit outputs
#
# Important restrictions:
#   No sample exclusion.
#   No gene exclusion.
#   No expression transformation.
#   No normalization.
#   No prevalence filtering.
#   No PCA or clustering.
#   No WGCNA network construction.
#
# Downstream:
#   Workflow-04-2B-1 will perform the prespecified expression
#   prevalence filtering and log2(FPKM-UQ + 1) transformation.
#
# Software:
#   R 4.5.1
#   Package: here
#   Base R functions are used for data import, QC, statistics,
#   plotting, and output.
#
# Random seed:
#   Not applicable. No stochastic procedure is performed.
#
# Reproducibility:
#   All paths are project-relative.
#   The script is designed to run in a clean R session.
############################################################


############################
# 0. Environment and setup
############################

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

if (!requireNamespace("here", quietly = TRUE)) {
  stop(
    "Required package 'here' is not installed. ",
    "Please install it before running this script."
  )
}

cat("\n============================================\n")
cat("Workflow-04-2B started\n")
cat("FPKM-UQ expression structure and distribution audit\n")
cat("============================================\n\n")

cat("R version: ", R.version.string, "\n", sep = "")
cat(
  "here version: ",
  as.character(utils::packageVersion("here")),
  "\n\n",
  sep = ""
)

project_dir <- here::here()

input_dir <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2A"
)

output_dir <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2B"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################
# 1. Input files
############################

expr_file <- file.path(
  input_dir,
  "Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt"
)

clinical_file <- file.path(
  input_dir,
  "Workflow-04-2A_clinical_aligned.txt"
)

sample_list_file <- file.path(
  input_dir,
  "Workflow-04-2A_primary_tumor_sample_list.txt"
)

input_files <- c(
  expr_file,
  clinical_file,
  sample_list_file
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {
  stop(
    "The following required input file(s) were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}


############################
# 2. Read input files
############################

cat("Reading Workflow-04-2A input files...\n")

expr_raw <- read.delim(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  na.strings = "NA"
)

clinical_aligned <- read.delim(
  clinical_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  na.strings = "NA"
)

sample_manifest <- read.delim(
  sample_list_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  na.strings = "NA"
)

if (ncol(expr_raw) < 2) {
  stop(
    "The expression input must contain one gene-ID column ",
    "and at least one sample column."
  )
}

if (!identical(colnames(expr_raw)[1], "Tag")) {
  stop(
    "The first expression-matrix column must be named 'Tag'. ",
    "Please check the Workflow-04-2A output."
  )
}

if (!"Sample" %in% colnames(clinical_aligned)) {
  stop(
    "The clinical-aligned input does not contain a 'Sample' column."
  )
}

if (!"Sample" %in% colnames(sample_manifest)) {
  stop(
    "The sample manifest does not contain a 'Sample' column."
  )
}

gene_id <- as.character(expr_raw[[1]])
sample_ids <- colnames(expr_raw)[-1]

clinical_ids <- as.character(
  clinical_aligned$Sample
)

manifest_ids <- as.character(
  sample_manifest$Sample
)


############################
# 3. Basic structure audit
############################

cat("\n--------------------------------------------\n")
cat("PART A: Basic expression matrix structure\n")
cat("--------------------------------------------\n")

n_genes <- nrow(expr_raw)
n_samples <- length(sample_ids)

cat("Number of gene rows: ", n_genes, "\n", sep = "")
cat("Number of sample columns: ", n_samples, "\n", sep = "")
cat("First column name: ", colnames(expr_raw)[1], "\n", sep = "")

cat("\nFirst 10 gene IDs:\n")
print(utils::head(gene_id, 10))

cat("\nFirst 10 sample IDs:\n")
print(utils::head(sample_ids, 10))

if (n_samples != 513) {
  stop(
    "Expected 513 primary tumor samples, but found ",
    n_samples,
    ". Please investigate the Workflow-04-2A output."
  )
}


############################
# 4. Gene ID audit
############################

cat("\n--------------------------------------------\n")
cat("PART B: Gene identifier audit\n")
cat("--------------------------------------------\n")

gene_id_missing <- is.na(gene_id)

gene_id_blank <- !gene_id_missing &
  trimws(gene_id) == ""

gene_id_duplicated <- duplicated(gene_id)

gene_id_audit <- data.frame(
  Metric = c(
    "Total_gene_IDs",
    "Missing_gene_IDs",
    "Blank_gene_IDs",
    "Duplicated_gene_IDs",
    "Unique_nonmissing_gene_IDs"
  ),
  Value = c(
    length(gene_id),
    sum(gene_id_missing),
    sum(gene_id_blank),
    sum(gene_id_duplicated),
    length(unique(gene_id[!gene_id_missing & !gene_id_blank]))
  )
)

print(gene_id_audit)

write.table(
  gene_id_audit,
  file.path(output_dir, "Workflow-04-2B_gene_ID_audit.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (
  any(gene_id_missing) ||
  any(gene_id_blank) ||
  any(gene_id_duplicated)
) {
  stop(
    "Gene identifiers contain missing, blank, or duplicated values. ",
    "No gene rows have been removed. Please investigate the input."
  )
}


############################
# 5. Sample ID and handoff audit
############################

cat("\n--------------------------------------------\n")
cat("PART C: Sample identifier and handoff audit\n")
cat("--------------------------------------------\n")

valid_tcga_barcode <- grepl(
  "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}$",
  sample_ids
)

sample_type_code <- substr(
  sample_ids,
  14,
  15
)

if (
  any(!valid_tcga_barcode) ||
  any(sample_type_code != "01")
) {
  stop(
    "Invalid TCGA sample barcode or non-primary-tumor ",
    "sample type detected. Expected sample type code 01."
  )
}

sample_id_missing <- is.na(sample_ids) |
  trimws(sample_ids) == ""

sample_id_duplicated <- duplicated(sample_ids)

clinical_id_missing <- is.na(clinical_ids) |
  trimws(clinical_ids) == ""

clinical_id_duplicated <- duplicated(clinical_ids)

manifest_id_missing <- is.na(manifest_ids) |
  trimws(manifest_ids) == ""

manifest_id_duplicated <- duplicated(manifest_ids)

sample_match_audit <- data.frame(
  Metric = c(
    "Expression_sample_count",
    "Clinical_row_count",
    "Manifest_sample_count",
    "Expression_missing_or_blank_IDs",
    "Expression_duplicated_IDs",
    "Clinical_missing_or_blank_IDs",
    "Clinical_duplicated_IDs",
    "Manifest_missing_or_blank_IDs",
    "Manifest_duplicated_IDs",
    "Expression_vs_clinical_order_identical",
    "Expression_vs_manifest_order_identical"
  ),
  Value = c(
    length(sample_ids),
    length(clinical_ids),
    length(manifest_ids),
    sum(sample_id_missing),
    sum(sample_id_duplicated),
    sum(clinical_id_missing),
    sum(clinical_id_duplicated),
    sum(manifest_id_missing),
    sum(manifest_id_duplicated),
    identical(sample_ids, clinical_ids),
    identical(sample_ids, manifest_ids)
  )
)

print(sample_match_audit)

write.table(
  sample_match_audit,
  file.path(output_dir, "Workflow-04-2B_sample_matching_audit.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sample_matching_ok <- (
  !any(sample_id_missing) &&
    !any(sample_id_duplicated) &&
    !any(clinical_id_missing) &&
    !any(clinical_id_duplicated) &&
    !any(manifest_id_missing) &&
    !any(manifest_id_duplicated) &&
    identical(sample_ids, clinical_ids) &&
    identical(sample_ids, manifest_ids)
)

if (!sample_matching_ok) {
  expression_not_in_clinical <- setdiff(
    sample_ids,
    clinical_ids
  )
  
  clinical_not_in_expression <- setdiff(
    clinical_ids,
    sample_ids
  )
  
  expression_not_in_manifest <- setdiff(
    sample_ids,
    manifest_ids
  )
  
  manifest_not_in_expression <- setdiff(
    manifest_ids,
    sample_ids
  )
  
  mismatch_details <- list(
    Expression_not_in_clinical = expression_not_in_clinical,
    Clinical_not_in_expression = clinical_not_in_expression,
    Expression_not_in_manifest = expression_not_in_manifest,
    Manifest_not_in_expression = manifest_not_in_expression
  )
  
  for (nm in names(mismatch_details)) {
    writeLines(
      c(
        paste0(nm, ":"),
        mismatch_details[[nm]]
      ),
      con = file.path(
        output_dir,
        paste0("Workflow-04-2B_", nm, ".txt")
      )
    )
  }
  
  stop(
    "Sample IDs or sample ordering are inconsistent across ",
    "the expression matrix, clinical table, and sample manifest. ",
    "No samples have been removed or reordered."
  )
}

cat(
  "PASS: Expression, clinical, and manifest sample IDs ",
  "are identical and in the same order.\n",
  sep = ""
)


############################
# 6. Clinical cohort check
############################

if (!"tumor_normal" %in% colnames(clinical_aligned)) {
  stop(
    "The clinical-aligned file does not contain 'tumor_normal'."
  )
}

tumor_status <- as.character(
  clinical_aligned$tumor_normal
)

if (
  any(is.na(tumor_status)) ||
  any(tumor_status != "tumor")
) {
  stop(
    "The clinical-aligned input contains missing or non-tumor ",
    "sample labels. No samples have been removed."
  )
}


############################
# 7. Strict numeric conversion
############################

cat("\n--------------------------------------------\n")
cat("PART D: Expression numeric-type audit\n")
cat("--------------------------------------------\n")

expression_columns <- expr_raw[, -1, drop = FALSE]

expr_mat <- matrix(
  NA_real_,
  nrow = nrow(expression_columns),
  ncol = ncol(expression_columns),
  dimnames = list(
    gene_id,
    sample_ids
  )
)

conversion_issues <- list()
issue_index <- 1L

for (j in seq_len(ncol(expression_columns))) {
  
  x <- expression_columns[[j]]
  
  if (is.numeric(x)) {
    
    expr_mat[, j] <- as.numeric(x)
    
  } else {
    
    x_character <- as.character(x)
    
    converted <- suppressWarnings(
      as.numeric(x_character)
    )
    
    recognized_nonfinite <- toupper(
      trimws(x_character)
    ) %in% c(
      "NAN",
      "INF",
      "+INF",
      "-INF",
      "INFINITY",
      "+INFINITY",
      "-INFINITY"
    )
    
    invalid <- (
      !is.na(x_character) &
        trimws(x_character) != "" &
        is.na(converted) &
        !recognized_nonfinite
    )
    
    if (any(invalid)) {
      
      invalid_rows <- which(invalid)
      
      conversion_issues[[issue_index]] <- data.frame(
        GeneID = gene_id[invalid_rows],
        Sample = sample_ids[j],
        OriginalValue = x_character[invalid_rows],
        stringsAsFactors = FALSE
      )
      
      issue_index <- issue_index + 1L
    }
    
    expr_mat[, j] <- converted
  }
}

if (length(conversion_issues) > 0) {
  
  conversion_issue_df <- do.call(
    rbind,
    conversion_issues
  )
  
  write.table(
    conversion_issue_df,
    file.path(
      output_dir,
      "Workflow-04-2B_numeric_conversion_issues.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  stop(
    "Non-numeric expression tokens were detected. ",
    "The offending values were saved to ",
    "Workflow-04-2B_numeric_conversion_issues.txt. ",
    "No expression values were silently discarded."
  )
}

if (!is.numeric(expr_mat)) {
  stop("Expression matrix conversion did not produce a numeric matrix.")
}

cat("Expression matrix storage mode: ", storage.mode(expr_mat), "\n", sep = "")
cat("Expression matrix class: ", paste(class(expr_mat), collapse = ", "), "\n", sep = "")
cat(
  "Numeric conversion completed without invalid text tokens.\n",
  sep = ""
)


############################
# 8. NA / NaN / Inf / negative audit
############################

cat("\n--------------------------------------------\n")
cat("PART E: Missing and invalid numeric values\n")
cat("--------------------------------------------\n")

n_total <- length(expr_mat)

n_nan <- sum(is.nan(expr_mat))

n_na_only <- sum(
  is.na(expr_mat) & !is.nan(expr_mat)
)

n_inf <- sum(is.infinite(expr_mat))

n_finite <- sum(is.finite(expr_mat))

n_negative <- sum(
  expr_mat < 0,
  na.rm = TRUE
)

n_zero <- sum(
  expr_mat == 0,
  na.rm = TRUE
)

n_positive <- sum(
  expr_mat > 0,
  na.rm = TRUE
)

n_nonfinite <- n_total - n_finite

global_value_audit <- data.frame(
  Metric = c(
    "Total_expression_values",
    "Finite_values",
    "NA_values_excluding_NaN",
    "NaN_values",
    "Positive_Infinity_values",
    "Negative_Infinity_values",
    "All_Infinite_values",
    "Nonfinite_values",
    "Negative_finite_values",
    "Zero_values",
    "Positive_finite_values",
    "Zero_fraction_among_finite_values"
  ),
  Value = c(
    n_total,
    n_finite,
    n_na_only,
    n_nan,
    sum(expr_mat == Inf, na.rm = TRUE),
    sum(expr_mat == -Inf, na.rm = TRUE),
    n_inf,
    n_nonfinite,
    n_negative,
    n_zero,
    n_positive,
    if (n_finite > 0) n_zero / n_finite else NA_real_
  )
)

print(global_value_audit)

write.table(
  global_value_audit,
  file.path(
    output_dir,
    "Workflow-04-2B_global_value_audit.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 9. Safe descriptive summary function
############################

safe_expression_summary <- function(x) {
  
  finite_x <- x[is.finite(x)]
  
  if (length(finite_x) == 0) {
    return(
      c(
        N_Finite = 0,
        N_NA = sum(is.na(x) & !is.nan(x)),
        N_NaN = sum(is.nan(x)),
        N_Inf = sum(is.infinite(x)),
        Min = NA_real_,
        Q01 = NA_real_,
        Q05 = NA_real_,
        Q25 = NA_real_,
        Median = NA_real_,
        Mean = NA_real_,
        SD = NA_real_,
        Q75 = NA_real_,
        Q90 = NA_real_,
        Q95 = NA_real_,
        Q99 = NA_real_,
        Q999 = NA_real_,
        Max = NA_real_,
        ZeroFraction = NA_real_,
        NegativeCount = 0
      )
    )
  }
  
  q <- stats::quantile(
    finite_x,
    probs = c(
      0,
      0.01,
      0.05,
      0.25,
      0.50,
      0.75,
      0.90,
      0.95,
      0.99,
      0.999,
      1
    ),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )
  
  c(
    N_Finite = length(finite_x),
    N_NA = sum(is.na(x) & !is.nan(x)),
    N_NaN = sum(is.nan(x)),
    N_Inf = sum(is.infinite(x)),
    Min = q[1],
    Q01 = q[2],
    Q05 = q[3],
    Q25 = q[4],
    Median = q[5],
    Mean = mean(finite_x),
    SD = if (length(finite_x) > 1) stats::sd(finite_x) else NA_real_,
    Q75 = q[6],
    Q90 = q[7],
    Q95 = q[8],
    Q99 = q[9],
    Q999 = q[10],
    Max = q[11],
    ZeroFraction = mean(finite_x == 0),
    NegativeCount = sum(finite_x < 0)
  )
}


############################
# 10. Global expression distribution
############################

cat("\n--------------------------------------------\n")
cat("PART F: Global expression distribution\n")
cat("--------------------------------------------\n")

finite_values <- as.numeric(
  expr_mat[is.finite(expr_mat)]
)

if (length(finite_values) == 0) {
  stop(
    "No finite expression values are available. ",
    "Distribution summaries and plots cannot be generated."
  )
}

global_summary_vector <- safe_expression_summary(
  as.numeric(expr_mat)
)

global_summary_df <- data.frame(
  Statistic = names(global_summary_vector),
  Value = as.numeric(global_summary_vector)
)

print(global_summary_df)

write.table(
  global_summary_df,
  file.path(
    output_dir,
    "Workflow-04-2B_global_expression_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 11. Per-sample expression summary
############################

cat("\n--------------------------------------------\n")
cat("PART G: Per-sample expression distribution\n")
cat("--------------------------------------------\n")

sample_summary_matrix <- t(
  vapply(
    seq_len(ncol(expr_mat)),
    function(j) {
      safe_expression_summary(expr_mat[, j])
    },
    FUN.VALUE = global_summary_vector
  )
)

sample_summary <- data.frame(
  Sample = sample_ids,
  sample_summary_matrix,
  check.names = FALSE
)

write.table(
  sample_summary,
  file.path(
    output_dir,
    "Workflow-04-2B_sample_expression_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nSummary of sample-level medians:\n")
print(summary(sample_summary$Median))

cat("\nSummary of sample-level means:\n")
print(summary(sample_summary$Mean))

cat("\nSummary of sample-level zero fractions:\n")
print(summary(sample_summary$ZeroFraction))


############################
# 12. Descriptive sample distribution tails
############################

get_tail_limits <- function(x) {
  
  valid_x <- x[is.finite(x)]
  
  if (length(valid_x) == 0) {
    return(c(Lower = NA_real_, Upper = NA_real_))
  }
  
  stats::quantile(
    valid_x,
    probs = c(0.01, 0.99),
    na.rm = TRUE,
    names = FALSE
  )
}

median_limits <- get_tail_limits(
  sample_summary$Median
)

zero_limits <- get_tail_limits(
  sample_summary$ZeroFraction
)

sample_summary$Median_Outside_01_99 <- (
  !is.na(sample_summary$Median) &
    (
      sample_summary$Median < median_limits[1] |
        sample_summary$Median > median_limits[2]
    )
)

sample_summary$ZeroFraction_Outside_01_99 <- (
  !is.na(sample_summary$ZeroFraction) &
    (
      sample_summary$ZeroFraction < zero_limits[1] |
        sample_summary$ZeroFraction > zero_limits[2]
    )
)

sample_median_tail <- sample_summary[
  sample_summary$Median_Outside_01_99,
  ,
  drop = FALSE
]

sample_zero_tail <- sample_summary[
  sample_summary$ZeroFraction_Outside_01_99,
  ,
  drop = FALSE
]

write.table(
  sample_median_tail,
  file.path(
    output_dir,
    "Workflow-04-2B_sample_median_distribution_tail.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sample_zero_tail,
  file.path(
    output_dir,
    "Workflow-04-2B_sample_zero_fraction_distribution_tail.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 13. Per-gene expression summary
############################

cat("\n--------------------------------------------\n")
cat("PART H: Per-gene expression distribution\n")
cat("--------------------------------------------\n")

gene_summary_matrix <- t(
  vapply(
    seq_len(nrow(expr_mat)),
    function(i) {
      safe_expression_summary(expr_mat[i, ])
    },
    FUN.VALUE = global_summary_vector
  )
)

gene_summary <- data.frame(
  GeneID = gene_id,
  gene_summary_matrix,
  check.names = FALSE
)

write.table(
  gene_summary,
  file.path(
    output_dir,
    "Workflow-04-2B_gene_expression_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nSummary of gene-level zero fractions:\n")
print(summary(gene_summary$ZeroFraction))


############################
# 14. Gene zero-fraction thresholds
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
  NumberOfGenes = vapply(
    zero_thresholds,
    function(threshold) {
      sum(
        gene_summary$ZeroFraction >= threshold,
        na.rm = TRUE
      )
    },
    numeric(1)
  )
)

write.table(
  gene_zero_counts,
  file.path(
    output_dir,
    "Workflow-04-2B_gene_zero_fraction_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 15. Highly expressed genes
############################

gene_order <- order(
  gene_summary$Median,
  decreasing = TRUE,
  na.last = NA
)

top_genes_median <- gene_summary[
  head(gene_order, 20),
  ,
  drop = FALSE
]

write.table(
  top_genes_median,
  file.path(
    output_dir,
    "Workflow-04-2B_top20_genes_by_median_expression.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 16. Low-expression gene summary
############################

n_gene_median_zero <- sum(
  gene_summary$Median == 0,
  na.rm = TRUE
)

low_expression_summary <- data.frame(
  Metric = c(
    "Total_genes",
    "Genes_with_median_FPKM_UQ_equal_to_zero",
    "Percentage_of_genes_with_median_equal_to_zero"
  ),
  Value = c(
    nrow(gene_summary),
    n_gene_median_zero,
    100 * n_gene_median_zero / nrow(gene_summary)
  )
)

write.table(
  low_expression_summary,
  file.path(
    output_dir,
    "Workflow-04-2B_low_expression_gene_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 17. Distribution plots
############################

cat("\n--------------------------------------------\n")
cat("PART I: Distribution plots\n")
cat("--------------------------------------------\n")

expr_mat_plot <- expr_mat
expr_mat_plot[!is.finite(expr_mat_plot)] <- NA_real_

plot_file <- file.path(
  output_dir,
  "Workflow-04-2B_expression_distribution_QC.pdf"
)

grDevices::pdf(
  plot_file,
  width = 10,
  height = 7,
  onefile = TRUE
)

old_par <- graphics::par(no.readonly = TRUE)

tryCatch({
  
  graphics::par(
    mar = c(5, 5, 4, 2),
    las = 1
  )
  
  # 17.1 Global raw FPKM-UQ distribution
  graphics::hist(
    finite_values,
    breaks = 100,
    main = "Global FPKM-UQ Expression Distribution",
    xlab = "Raw FPKM-UQ",
    ylab = "Frequency"
  )
  
  # 17.2 Distribution below the global 99th percentile
  global_99 <- stats::quantile(
    finite_values,
    probs = 0.99,
    names = FALSE
  )
  
  lower_range_values <- finite_values[
    finite_values <= global_99
  ]
  
  graphics::hist(
    lower_range_values,
    breaks = 100,
    main = "FPKM-UQ Distribution up to the 99th Percentile",
    xlab = "Raw FPKM-UQ",
    ylab = "Frequency"
  )
  
  # 17.3 Sample median distribution
  valid_sample_medians <- sample_summary$Median[
    is.finite(sample_summary$Median)
  ]
  
  graphics::hist(
    valid_sample_medians,
    breaks = 50,
    main = "Distribution of Sample Median FPKM-UQ",
    xlab = "Sample median FPKM-UQ",
    ylab = "Number of samples"
  )
  
  # 17.4 Sample zero fraction
  valid_sample_zero <- sample_summary$ZeroFraction[
    is.finite(sample_summary$ZeroFraction)
  ]
  
  graphics::hist(
    valid_sample_zero,
    breaks = 50,
    main = "Distribution of Sample Zero Fractions",
    xlab = "Fraction of zero expression values",
    ylab = "Number of samples"
  )
  
  # 17.5 Sample-level expression boxplots
  # Non-finite values are masked for plotting only.
  # The original expr_mat remains unchanged.
  
  graphics::boxplot(
    expr_mat_plot,
    outline = FALSE,
    las = 2,
    cex.axis = 0.35,
    main = "Raw FPKM-UQ Distribution Across Samples",
    xlab = "Primary tumor samples",
    ylab = "Raw FPKM-UQ"
  )
  
  # 17.6 Ranked sample medians
  sample_order <- order(
    sample_summary$Median,
    na.last = NA
  )
  
  graphics::plot(
    seq_along(sample_order),
    sample_summary$Median[sample_order],
    type = "p",
    pch = 16,
    main = "Sample Median FPKM-UQ (Ranked)",
    xlab = "Samples ranked by median expression",
    ylab = "Median FPKM-UQ"
  )
  
  # 17.7 Ranked sample zero fractions
  zero_order <- order(
    sample_summary$ZeroFraction,
    na.last = NA
  )
  
  graphics::plot(
    seq_along(zero_order),
    sample_summary$ZeroFraction[zero_order],
    type = "p",
    pch = 16,
    main = "Sample Zero Fraction (Ranked)",
    xlab = "Samples ranked by zero fraction",
    ylab = "Zero fraction"
  )
  
}, finally = {
  
  graphics::par(old_par)
  grDevices::dev.off()
  
})


############################
# 18. Final audit summary
############################

cat("\n============================================\n")
cat("Workflow-04-2B FINAL AUDIT SUMMARY\n")
cat("============================================\n")

audit_summary <- data.frame(
  Metric = c(
    "Number_of_genes",
    "Number_of_samples",
    "Total_expression_values",
    "Finite_values",
    "NA_values_excluding_NaN",
    "NaN_values",
    "Inf_values",
    "Negative_values",
    "Zero_values",
    "Positive_values",
    "Zero_fraction_among_finite_values",
    "Global_median_FPKM_UQ",
    "Global_mean_FPKM_UQ",
    "Global_99th_percentile",
    "Global_99.9th_percentile",
    "Global_maximum",
    "Genes_with_median_equal_to_zero",
    "Samples_with_median_in_distribution_tail",
    "Samples_with_zero_fraction_in_distribution_tail",
    "Sample_ID_and_order_check"
  ),
  Value = c(
    nrow(expr_mat),
    ncol(expr_mat),
    n_total,
    n_finite,
    n_na_only,
    n_nan,
    n_inf,
    n_negative,
    n_zero,
    n_positive,
    n_zero / n_finite,
    median(finite_values),
    mean(finite_values),
    stats::quantile(finite_values, 0.99, names = FALSE),
    stats::quantile(finite_values, 0.999, names = FALSE),
    max(finite_values),
    n_gene_median_zero,
    nrow(sample_median_tail),
    nrow(sample_zero_tail),
    "PASS"
  )
)

print(audit_summary)

write.table(
  audit_summary,
  file.path(
    output_dir,
    "Workflow-04-2B_audit_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 19. Final workflow report
############################

report_lines <- c(
  "Workflow-04-2B: FPKM-UQ Expression Data Structure and Distribution Audit",
  paste0("Date: ", Sys.Date()),
  paste0("R version: ", R.version.string),
  paste0(
    "here package version: ",
    as.character(utils::packageVersion("here"))
  ),
  "",
  paste0("Input expression file: ", basename(expr_file)),
  paste0("Input clinical file: ", basename(clinical_file)),
  paste0("Input sample manifest: ", basename(sample_list_file)),
  "",
  paste0("Number of genes: ", nrow(expr_mat)),
  paste0("Number of samples: ", ncol(expr_mat)),
  paste0("Finite expression values: ", n_finite),
  paste0("NA values excluding NaN: ", n_na_only),
  paste0("NaN values: ", n_nan),
  paste0("Infinite values: ", n_inf),
  paste0("Negative values: ", n_negative),
  paste0("Zero values: ", n_zero),
  paste0(
    "Zero fraction among finite values: ",
    signif(n_zero / n_finite, 6)
  ),
  "",
  "Sample matching: PASS",
  "Clinical sample order: PASS",
  "Sample manifest order: PASS",
  "",
  "No samples were excluded.",
  "No genes were excluded.",
  "No expression transformation was performed.",
  "No normalization was performed.",
  "No prevalence filtering was performed.",
  "No PCA or WGCNA was performed.",
  "",
  "Distribution-tail flags are descriptive only.",
  "They must not be interpreted as automatic exclusion criteria.",
  "",
  "Downstream handoff:",
  "Workflow-04-2B-1 may use the original expression matrix",
  "after reviewing this QC report and confirming that any",
  "non-finite or negative values do not require investigation."
)

writeLines(
  report_lines,
  con = file.path(
    output_dir,
    "Workflow-04-2B_workflow_report.txt"
  )
)

cat("\nWorkflow-04-2B completed.\n")
cat("No samples or genes were removed.\n")
cat("No transformation, normalization, or filtering was performed.\n")
cat(
  "All QC outputs were saved to:\n",
  output_dir,
  "\n",
  sep = ""
)

############################################################
# Workflow-04-2B-1
# Extreme Expression and Sample-Level Provenance Audit
#
# Project:
#   Paper 1 - LUAD transcriptomic analysis
#
# Code path:
#   scripts/02_tumor_only_WGCNA/
#   GitHub2_Paper1_tumor_only_WGCNA_1.R
#
# Research purpose:
#   To investigate extreme expression values and characterize
#   sample-level expression distributions in the 513 primary
#   LUAD tumor samples.
#
# Input files:
#   1. Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt
#   2. Workflow-04-2A_clinical_aligned.txt
#   3. Workflow-04-2A_primary_tumor_sample_list.txt
#
# Input directory:
#   results/02_tumor_only_WGCNA/Workflow-04-2A/
#
# Expected cohort:
#   513 TCGA primary LUAD tumor samples (sample type code 01).
#
# Analysis scope:
#   1. Global maximum expression value
#   2. Top 100 individual gene-sample expression observations
#   3. Sample-level maximum expression values
#   4. Counts above exploratory high-expression thresholds
#   5. TCGA-44 sample-level distribution audit
#   6. Prespecified sample-level distribution assessment
#   7. Top expressed genes in selected samples
#   8. Preliminary Ensembl ID pattern assessment
#   9. Relative sample distribution and barcode-prefix summary
#
# Sample inclusion:
#   All 513 samples retained in Workflow-04-2A.
#   Exact TCGA sample-level barcodes are required.
#   Expression, clinical, and manifest sample ordering must agree.
#
# Important restrictions:
#   No sample exclusion.
#   No gene exclusion.
#   No expression transformation.
#   No normalization.
#   No prevalence filtering.
#   No PCA or clustering.
#   No WGCNA.
#
# Downstream:
#   Workflow-04-2B-2 uses the Top 100 extreme observations
#   generated here for exact Ensembl feature annotation.
#
# Software:
#   R 4.5.1
#   Package: here
#   Base R
#
# Random seed:
#   Not applicable. No stochastic procedure is performed.
#
# Reproducibility:
#   Project-relative paths.
#   Designed to run in a clean R session.
############################################################


############################
# 0. Environment
############################

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

if (!requireNamespace("here", quietly = TRUE)) {
  stop(
    "Package 'here' is required. ",
    "Please install it before running this workflow."
  )
}

cat("\n============================================\n")
cat("Workflow-04-2B-1 started\n")
cat("Extreme expression and sample-level provenance audit\n")
cat("============================================\n\n")

cat("R version: ", R.version.string, "\n", sep = "")
cat(
  "here version: ",
  as.character(utils::packageVersion("here")),
  "\n\n",
  sep = ""
)


############################
# 1. Project directories
############################

project_dir <- here::here()

input_dir <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2A"
)

output_dir <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2B-1"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################
# 2. Input files
############################

expr_file <- file.path(
  input_dir,
  "Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt"
)

clinical_file <- file.path(
  input_dir,
  "Workflow-04-2A_clinical_aligned.txt"
)

sample_list_file <- file.path(
  input_dir,
  "Workflow-04-2A_primary_tumor_sample_list.txt"
)

input_files <- c(
  expr_file,
  clinical_file,
  sample_list_file
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing input file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}


############################
# 3. Read input data
############################

expr_raw <- read.delim(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  na.strings = "NA"
)

clinical_aligned <- read.delim(
  clinical_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  na.strings = "NA"
)

sample_manifest <- read.delim(
  sample_list_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = "",
  na.strings = "NA"
)

if (ncol(expr_raw) < 2) {
  stop("Expression matrix must contain gene IDs and sample columns.")
}

if (!identical(colnames(expr_raw)[1], "Tag")) {
  stop("The first expression matrix column must be named 'Tag'.")
}

if (!"Sample" %in% colnames(clinical_aligned)) {
  stop("Clinical-aligned file lacks the 'Sample' column.")
}

if (!"Sample" %in% colnames(sample_manifest)) {
  stop("Sample manifest lacks the 'Sample' column.")
}

gene_id <- as.character(expr_raw[[1]])
sample_ids <- colnames(expr_raw)[-1]

clinical_ids <- as.character(
  clinical_aligned$Sample
)

manifest_ids <- as.character(
  sample_manifest$Sample
)


############################
# 4. Input structure and ID checks
############################

cat("--------------------------------------------\n")
cat("PART A: Input structure and sample identity\n")
cat("--------------------------------------------\n")

if (nrow(expr_raw) == 0) {
  stop("Expression matrix contains no gene rows.")
}

if (length(sample_ids) != 513) {
  stop(
    "Expected 513 primary tumor samples; found ",
    length(sample_ids),
    "."
  )
}

if (
  any(is.na(gene_id)) ||
  any(trimws(gene_id) == "") ||
  anyDuplicated(gene_id) > 0
) {
  stop(
    "Gene IDs contain missing, blank, or duplicated values. ",
    "No genes have been removed."
  )
}

if (
  any(is.na(sample_ids)) ||
  any(trimws(sample_ids) == "") ||
  anyDuplicated(sample_ids) > 0
) {
  stop("Expression sample IDs are missing, blank, or duplicated.")
}

valid_barcode <- grepl(
  "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}$",
  sample_ids
)

sample_type <- substr(
  sample_ids,
  14,
  15
)

if (
  any(!valid_barcode) ||
  any(sample_type != "01")
) {
  stop(
    "Invalid TCGA barcode or non-primary-tumor sample detected. ",
    "Expected sample type code 01."
  )
}

if (
  !identical(sample_ids, clinical_ids) ||
  !identical(sample_ids, manifest_ids)
) {
  stop(
    "Expression, clinical, and sample manifest IDs or ordering ",
    "are inconsistent. No samples have been removed or reordered."
  )
}

if (!"tumor_normal" %in% colnames(clinical_aligned)) {
  stop("Clinical data lack the 'tumor_normal' column.")
}

if (
  any(is.na(clinical_aligned$tumor_normal)) ||
  any(as.character(clinical_aligned$tumor_normal) != "tumor")
) {
  stop("Clinical-aligned input contains missing or non-tumor labels.")
}

cat("Genes: ", length(gene_id), "\n", sep = "")
cat("Samples: ", length(sample_ids), "\n", sep = "")
cat("Sample ID/order check: PASS\n")


############################
# 5. Strict numeric conversion
############################

cat("\n--------------------------------------------\n")
cat("PART B: Numeric expression integrity\n")
cat("--------------------------------------------\n")

expression_df <- expr_raw[, -1, drop = FALSE]

expr_mat <- matrix(
  NA_real_,
  nrow = nrow(expression_df),
  ncol = ncol(expression_df),
  dimnames = list(gene_id, sample_ids)
)

conversion_issues <- list()
issue_index <- 1L

for (j in seq_len(ncol(expression_df))) {
  
  x <- expression_df[[j]]
  
  if (is.numeric(x)) {
    
    expr_mat[, j] <- as.numeric(x)
    
  } else {
    
    x_chr <- as.character(x)
    
    converted <- suppressWarnings(
      as.numeric(x_chr)
    )
    
    recognized_nonfinite <- toupper(
      trimws(x_chr)
    ) %in% c(
      "NAN",
      "INF",
      "+INF",
      "-INF",
      "INFINITY",
      "+INFINITY",
      "-INFINITY"
    )
    
    invalid <- (
      !is.na(x_chr) &
        trimws(x_chr) != "" &
        is.na(converted) &
        !recognized_nonfinite
    )
    
    if (any(invalid)) {
      
      idx <- which(invalid)
      
      conversion_issues[[issue_index]] <- data.frame(
        GeneID = gene_id[idx],
        Sample = sample_ids[j],
        OriginalValue = x_chr[idx],
        stringsAsFactors = FALSE
      )
      
      issue_index <- issue_index + 1L
    }
    
    expr_mat[, j] <- converted
  }
}

if (length(conversion_issues) > 0) {
  
  conversion_issue_df <- do.call(
    rbind,
    conversion_issues
  )
  
  write.table(
    conversion_issue_df,
    file.path(
      output_dir,
      "Workflow-04-2B-1_numeric_conversion_issues.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  stop(
    "Invalid non-numeric expression tokens detected. ",
    "See Workflow-04-2B-1_numeric_conversion_issues.txt."
  )
}

cat("Numeric conversion: PASS\n")
cat("Storage mode: ", storage.mode(expr_mat), "\n", sep = "")


############################
# 6. Finite-value assessment
############################

n_na <- sum(
  is.na(expr_mat) & !is.nan(expr_mat)
)

n_nan <- sum(is.nan(expr_mat))
n_inf <- sum(is.infinite(expr_mat))
n_finite <- sum(is.finite(expr_mat))

if (n_finite == 0) {
  stop("No finite expression values are available.")
}

finite_values <- as.numeric(
  expr_mat[is.finite(expr_mat)]
)

n_negative <- sum(
  finite_values < 0
)

value_integrity <- data.frame(
  Metric = c(
    "Total_values",
    "Finite_values",
    "NA_excluding_NaN",
    "NaN",
    "Inf",
    "Negative_finite_values"
  ),
  Value = c(
    length(expr_mat),
    n_finite,
    n_na,
    n_nan,
    n_inf,
    n_negative
  )
)

write.table(
  value_integrity,
  file.path(
    output_dir,
    "Workflow-04-2B-1_expression_value_integrity.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (n_negative > 0) {
  warning(
    "Negative finite expression values were detected. ",
    "They are retained and reported; investigate before ",
    "downstream expression filtering."
  )
}

if (n_na + n_nan + n_inf > 0) {
  warning(
    "Missing or non-finite values were detected. ",
    "They are retained and reported."
  )
}


############################
# 7. Safe summary function
############################

safe_summary <- function(x) {
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(
      c(
        N_Finite = 0,
        Min = NA_real_,
        Q1 = NA_real_,
        Median = NA_real_,
        Mean = NA_real_,
        Q3 = NA_real_,
        P90 = NA_real_,
        P95 = NA_real_,
        P99 = NA_real_,
        P999 = NA_real_,
        Max = NA_real_,
        ZeroFraction = NA_real_
      )
    )
  }
  
  q <- stats::quantile(
    x,
    probs = c(
      0,
      0.25,
      0.50,
      0.75,
      0.90,
      0.95,
      0.99,
      0.999,
      1
    ),
    names = FALSE,
    na.rm = TRUE
  )
  
  c(
    N_Finite = length(x),
    Min = q[1],
    Q1 = q[2],
    Median = q[3],
    Mean = mean(x),
    Q3 = q[4],
    P90 = q[5],
    P95 = q[6],
    P99 = q[7],
    P999 = q[8],
    Max = q[9],
    ZeroFraction = mean(x == 0)
  )
}


############################
# 8. Global maximum and Top 100
############################

cat("\n--------------------------------------------\n")
cat("PART C: Global extreme expression values\n")
cat("--------------------------------------------\n")

max_value <- max(finite_values)

max_index <- which(
  is.finite(expr_mat) &
    expr_mat == max_value,
  arr.ind = TRUE
)

max_info <- data.frame(
  GeneID = gene_id[max_index[, 1]],
  Sample = sample_ids[max_index[, 2]],
  Expression = expr_mat[max_index],
  stringsAsFactors = FALSE
)

write.table(
  max_info,
  file.path(
    output_dir,
    "Workflow-04-2B-1_global_maximum.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

finite_index <- which(
  is.finite(expr_mat),
  arr.ind = TRUE
)

finite_observations <- data.frame(
  GeneID = gene_id[finite_index[, 1]],
  Sample = sample_ids[finite_index[, 2]],
  Expression = expr_mat[finite_index],
  stringsAsFactors = FALSE
)

finite_observations <- finite_observations[
  order(
    finite_observations$Expression,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

n_top <- min(
  100L,
  nrow(finite_observations)
)

top100_extreme <- finite_observations[
  seq_len(n_top),
  ,
  drop = FALSE
]

top100_extreme$Rank <- seq_len(n_top)

top100_extreme <- top100_extreme[
  ,
  c("Rank", "GeneID", "Sample", "Expression")
]

write.table(
  top100_extreme,
  file.path(
    output_dir,
    "Workflow-04-2B-1_top100_extreme_values.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Global maximum: ", max_value, "\n", sep = "")
cat("Number of Top observations: ", n_top, "\n", sep = "")


############################
# 9. Sample-level expression summaries
############################

cat("\n--------------------------------------------\n")
cat("PART D: Sample-level expression summaries\n")
cat("--------------------------------------------\n")

sample_summary_matrix <- t(
  vapply(
    seq_len(ncol(expr_mat)),
    function(j) {
      safe_summary(expr_mat[, j])
    },
    FUN.VALUE = safe_summary(expr_mat[, 1])
  )
)

sample_summary <- data.frame(
  Sample = sample_ids,
  sample_summary_matrix,
  check.names = FALSE
)

sample_summary <- sample_summary[
  order(
    sample_summary$Max,
    decreasing = TRUE,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

write.table(
  sample_summary,
  file.path(
    output_dir,
    "Workflow-04-2B-1_sample_expression_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sample_max_df <- sample_summary[
  order(
    sample_summary$Max,
    decreasing = TRUE,
    na.last = TRUE
  ),
  c("Sample", "Max"),
  drop = FALSE
]

colnames(sample_max_df)[2] <- "MaxExpression"

write.table(
  sample_max_df,
  file.path(
    output_dir,
    "Workflow-04-2B-1_sample_max_expression.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 10. Exploratory high-expression thresholds
############################

cat("\n--------------------------------------------\n")
cat("PART E: Exploratory high-expression thresholds\n")
cat("--------------------------------------------\n")

thresholds <- c(
  1e6,
  1e7,
  1e8,
  1e9
)

sample_extreme_counts <- data.frame(
  Sample = sample_ids,
  stringsAsFactors = FALSE
)

for (thr in thresholds) {
  
  count_name <- paste0(
    "N_ge_",
    format(
      thr,
      scientific = FALSE,
      trim = TRUE
    )
  )
  
  sample_extreme_counts[[count_name]] <- vapply(
    seq_len(ncol(expr_mat)),
    function(j) {
      sum(
        is.finite(expr_mat[, j]) &
          expr_mat[, j] >= thr
      )
    },
    numeric(1)
  )
}

write.table(
  sample_extreme_counts,
  file.path(
    output_dir,
    "Workflow-04-2B-1_sample_extreme_expression_counts.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "Threshold counts are descriptive and are not used ",
  "for sample exclusion.\n",
  sep = ""
)


############################
# 11. TCGA-44 sample audit
############################

cat("\n--------------------------------------------\n")
cat("PART F: TCGA-44 sample audit\n")
cat("--------------------------------------------\n")

tcga44_samples <- grep(
  "^TCGA-44-",
  sample_ids,
  value = TRUE
)

tcga44_summary <- sample_summary[
  match(
    tcga44_samples,
    sample_summary$Sample
  ),
  ,
  drop = FALSE
]

tcga44_summary <- tcga44_summary[
  order(
    tcga44_summary$Median,
    decreasing = TRUE,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

write.table(
  data.frame(Sample = tcga44_samples),
  file.path(
    output_dir,
    "Workflow-04-2B-1_TCGA44_sample_list.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  tcga44_summary,
  file.path(
    output_dir,
    "Workflow-04-2B-1_TCGA44_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Number of TCGA-44 samples: ", length(tcga44_samples), "\n", sep = "")


############################
# 12. Prespecified sample audit
############################

cat("\n--------------------------------------------\n")
cat("PART G: Prespecified sample audit\n")
cat("--------------------------------------------\n")

selected_samples <- c(
  "TCGA-44-2662-01",
  "TCGA-44-2656-01",
  "TCGA-44-6147-01",
  "TCGA-44-4112-01",
  "TCGA-44-6146-01",
  "TCGA-44-3917-01"
)

selected_found <- selected_samples[
  selected_samples %in% sample_ids
]

selected_not_found <- selected_samples[
  !selected_samples %in% sample_ids
]

selected_sample_manifest <- data.frame(
  Sample = selected_samples,
  Found_in_expression_matrix = (
    selected_samples %in% sample_ids
  ),
  stringsAsFactors = FALSE
)

write.table(
  selected_sample_manifest,
  file.path(
    output_dir,
    "Workflow-04-2B-1_selected_sample_presence.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

selected_distribution <- sample_summary[
  match(
    selected_found,
    sample_summary$Sample
  ),
  ,
  drop = FALSE
]

selected_distribution <- selected_distribution[
  ,
  c(
    "Sample",
    "N_Finite",
    "Median",
    "Mean",
    "Q1",
    "Q3",
    "P90",
    "P95",
    "P99",
    "P999",
    "Max",
    "ZeroFraction"
  ),
  drop = FALSE
]

write.table(
  selected_distribution,
  file.path(
    output_dir,
    "Workflow-04-2B-1_selected_sample_distribution.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 13. Top genes in selected samples
############################

cat("\n--------------------------------------------\n")
cat("PART H: Top expressed genes in selected samples\n")
cat("--------------------------------------------\n")

top_genes_selected <- list()

for (s in selected_found) {
  
  x <- expr_mat[, s]
  
  valid_idx <- which(is.finite(x))
  
  valid_idx <- valid_idx[
    order(
      x[valid_idx],
      decreasing = TRUE
    )
  ]
  
  n_top_gene <- min(
    30L,
    length(valid_idx)
  )
  
  idx <- valid_idx[
    seq_len(n_top_gene)
  ]
  
  top_genes_selected[[s]] <- data.frame(
    Sample = s,
    Rank = seq_len(n_top_gene),
    GeneID = gene_id[idx],
    Expression = x[idx],
    stringsAsFactors = FALSE
  )
}

if (length(top_genes_selected) > 0) {
  
  top_genes_selected_df <- do.call(
    rbind,
    top_genes_selected
  )
  
} else {
  
  top_genes_selected_df <- data.frame(
    Sample = character(),
    Rank = integer(),
    GeneID = character(),
    Expression = numeric()
  )
}

write.table(
  top_genes_selected_df,
  file.path(
    output_dir,
    "Workflow-04-2B-1_selected_sample_top30_genes.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 14. Preliminary Ensembl ID pattern audit
############################

cat("\n--------------------------------------------\n")
cat("PART I: Preliminary Ensembl ID pattern audit\n")
cat("--------------------------------------------\n")

ensg_flag <- grepl(
  "^ENSG",
  gene_id
)

mito_pattern_flag <- grepl(
  "^ENSG000001988|^ENSG000002100",
  gene_id
)

gene_pattern_audit <- data.frame(
  Metric = c(
    "Total_gene_IDs",
    "IDs_matching_ENSG_prefix",
    "IDs_matching_preliminary_mitochondrial_pattern"
  ),
  Value = c(
    length(gene_id),
    sum(ensg_flag),
    sum(mito_pattern_flag)
  )
)

write.table(
  gene_pattern_audit,
  file.path(
    output_dir,
    "Workflow-04-2B-1_gene_ID_pattern_audit.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

top100_extreme$PreliminaryMitoPattern <- mito_pattern_flag[
  match(
    top100_extreme$GeneID,
    gene_id
  )
]

write.table(
  top100_extreme,
  file.path(
    output_dir,
    "Workflow-04-2B-1_top100_extreme_values.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

mito_top100_summary <- as.data.frame(
  table(
    top100_extreme$PreliminaryMitoPattern,
    useNA = "ifany"
  )
)

colnames(mito_top100_summary) <- c(
  "PreliminaryMitoPattern",
  "NumberOfObservations"
)

write.table(
  mito_top100_summary,
  file.path(
    output_dir,
    "Workflow-04-2B-1_top100_mito_pattern_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 15. Relative position among all samples
############################

cat("\n--------------------------------------------\n")
cat("PART J: Relative sample distribution\n")
cat("--------------------------------------------\n")

all_sample_summary <- sample_summary[
  match(
    sample_ids,
    sample_summary$Sample
  ),
  ,
  drop = FALSE
]

selected_relative <- all_sample_summary[
  match(
    selected_found,
    all_sample_summary$Sample
  ),
  ,
  drop = FALSE
]

selected_relative$MedianPercentile <- vapply(
  selected_relative$Median,
  function(x) {
    if (!is.finite(x)) {
      return(NA_real_)
    }
    
    mean(
      all_sample_summary$Median <= x,
      na.rm = TRUE
    )
  },
  numeric(1)
)

selected_relative$MeanPercentile <- vapply(
  selected_relative$Mean,
  function(x) {
    if (!is.finite(x)) {
      return(NA_real_)
    }
    
    mean(
      all_sample_summary$Mean <= x,
      na.rm = TRUE
    )
  },
  numeric(1)
)

write.table(
  selected_relative,
  file.path(
    output_dir,
    "Workflow-04-2B-1_selected_sample_relative_position.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 16. Barcode-prefix summary
############################

cat("\n--------------------------------------------\n")
cat("PART K: Barcode-prefix descriptive summary\n")
cat("--------------------------------------------\n")

patient_prefix <- sub(
  "^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})-[0-9]{2}$",
  "\\1",
  sample_ids
)

sample_prefix_summary <- data.frame(
  Sample = sample_ids,
  Prefix = patient_prefix,
  Median = all_sample_summary$Median,
  Max = all_sample_summary$Max,
  ZeroFraction = all_sample_summary$ZeroFraction,
  stringsAsFactors = FALSE
)

prefix_summary <- aggregate(
  cbind(
    Median,
    Max,
    ZeroFraction
  ) ~ Prefix,
  data = sample_prefix_summary,
  FUN = function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }
    
    median(x, na.rm = TRUE)
  }
)

prefix_sample_counts <- aggregate(
  Sample ~ Prefix,
  data = sample_prefix_summary,
  FUN = length
)

colnames(prefix_sample_counts)[2] <- "NumberOfSamples"

prefix_summary <- merge(
  prefix_summary,
  prefix_sample_counts,
  by = "Prefix",
  all.x = TRUE,
  sort = FALSE
)

prefix_summary <- prefix_summary[
  order(
    prefix_summary$Median,
    decreasing = TRUE,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

write.table(
  prefix_summary,
  file.path(
    output_dir,
    "Workflow-04-2B-1_patient_prefix_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 17. Final audit summary
############################

audit_summary <- data.frame(
  Metric = c(
    "Number_of_genes",
    "Number_of_samples",
    "Global_maximum_expression",
    "Number_of_top100_observations",
    "Number_of_TCGA44_samples",
    "Number_of_selected_samples_found",
    "Number_of_selected_samples_not_found",
    "NA_values_excluding_NaN",
    "NaN_values",
    "Inf_values",
    "Negative_finite_values",
    "Genes_matching_ENSG_prefix",
    "Genes_matching_preliminary_mito_pattern",
    "Sample_ID_and_order_check",
    "Sample_or_gene_exclusion"
  ),
  Value = c(
    nrow(expr_mat),
    ncol(expr_mat),
    max_value,
    n_top,
    length(tcga44_samples),
    length(selected_found),
    length(selected_not_found),
    n_na,
    n_nan,
    n_inf,
    n_negative,
    sum(ensg_flag),
    sum(mito_pattern_flag),
    "PASS",
    "NONE"
  )
)

write.table(
  audit_summary,
  file.path(
    output_dir,
    "Workflow-04-2B-1_audit_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 18. Workflow report
############################

report_lines <- c(
  "Workflow-04-2B-1: Extreme Expression and Sample-Level Provenance Audit",
  paste0("Date: ", Sys.Date()),
  paste0("R version: ", R.version.string),
  paste0(
    "here package version: ",
    as.character(utils::packageVersion("here"))
  ),
  "",
  paste0("Number of genes: ", nrow(expr_mat)),
  paste0("Number of samples: ", ncol(expr_mat)),
  paste0("Global maximum expression: ", signif(max_value, 8)),
  paste0("Top observations reported: ", n_top),
  paste0("TCGA-44 samples: ", length(tcga44_samples)),
  paste0("Selected samples found: ", length(selected_found)),
  paste0("Selected samples not found: ", length(selected_not_found)),
  "",
  paste0("NA excluding NaN: ", n_na),
  paste0("NaN: ", n_nan),
  paste0("Inf: ", n_inf),
  paste0("Negative finite values: ", n_negative),
  "",
  "All 513 primary tumor samples were retained.",
  "All input gene rows were retained.",
  "No expression transformation was performed.",
  "No normalization was performed.",
  "No prevalence filtering was performed.",
  "No PCA or WGCNA was performed.",
  "",
  "High-expression thresholds are exploratory descriptive summaries.",
  "Preliminary mitochondrial ID patterns are not formal gene annotations.",
  "No sample or gene exclusion decision was made in this workflow.",
  "",
  "Downstream handoff:",
  "Workflow-04-2B-2 uses Workflow-04-2B-1_top100_extreme_values.txt."
)

writeLines(
  report_lines,
  con = file.path(
    output_dir,
    "Workflow-04-2B-1_workflow_report.txt"
  )
)


############################
# 19. Completion
############################

cat("\n============================================\n")
cat("Workflow-04-2B-1 completed\n")
cat("============================================\n")

cat("Genes retained: ", nrow(expr_mat), "\n", sep = "")
cat("Samples retained: ", ncol(expr_mat), "\n", sep = "")
cat("No samples or genes were removed.\n")
cat("No transformation or normalization was performed.\n")
cat("No WGCNA was performed.\n")
cat("Output directory:\n", output_dir, "\n", sep = "")


############################################################
# Workflow-04-2B-2
# Annotation of extreme-expression features
#
# Research purpose:
#   Annotate the Top 100 extreme gene-expression observations
#   identified in Workflow-04-2B-1 and audit the annotation
#   completeness of the corresponding Ensembl gene IDs.
#
# Input:
#   results/02_tumor_only_WGCNA/Workflow-04-2B-1/
#   Workflow-04-2B-1_top100_extreme_values.txt
#
# Code path:
#   scripts/02_tumor_only_WGCNA/
#
# Output:
#   results/02_tumor_only_WGCNA/Workflow-04-2B-2/
#
# Software environment:
#   R 4.5.1
#   AnnotationDbi
#   org.Hs.eg.db
#   here
#
# Reproducibility:
#   No random procedure; random seed not applicable.
#   Package versions and annotation database metadata are recorded.
#
# Sample inclusion:
#   Uses the Top 100 observations from the 513 primary LUAD
#   samples audited in Workflow-04-2B-1.
#
# Important:
#   - No samples removed.
#   - No genes removed.
#   - No expression transformation.
#   - No normalization.
#   - No WGCNA.
#   - org.Hs.eg.db annotation is supplementary and is not
#     treated as a definitive gene-biotype annotation.
############################################################


############################
# 0. Initialize environment
############################

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install it before running this script.")
}

if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  stop("Package 'AnnotationDbi' is required. Install it before running this script.")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  stop("Package 'org.Hs.eg.db' is required. Install it before running this script.")
}

project_dir <- here::here()

input_file <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2B-1",
  "Workflow-04-2B-1_top100_extreme_values.txt"
)

output_dir <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2B-2"
)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


############################
# 1. Record software versions
############################

package_versions <- data.frame(
  Package = c(
    "R",
    "here",
    "AnnotationDbi",
    "org.Hs.eg.db"
  ),
  Version = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("here")),
    as.character(utils::packageVersion("AnnotationDbi")),
    as.character(utils::packageVersion("org.Hs.eg.db"))
  ),
  stringsAsFactors = FALSE
)

write.table(
  package_versions,
  file.path(output_dir, "Workflow-04-2B-2_software_versions.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 2. Load Top 100 observations
############################

top100 <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "Rank",
  "GeneID",
  "Sample",
  "Expression"
)

missing_columns <- setdiff(
  required_columns,
  colnames(top100)
)

if (length(missing_columns) > 0) {
  stop(
    "Required columns missing from Top 100 input: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (nrow(top100) == 0) {
  stop("The Top 100 input file contains no observations.")
}

if (any(is.na(top100$GeneID) | trimws(top100$GeneID) == "")) {
  stop("Missing or blank GeneID values detected.")
}

top100$Rank <- suppressWarnings(as.integer(top100$Rank))
top100$Expression <- suppressWarnings(as.numeric(top100$Expression))

if (anyNA(top100$Rank)) {
  stop("Invalid Rank values detected.")
}

if (any(!is.finite(top100$Expression))) {
  stop("Non-finite Expression values detected in Top 100 input.")
}

top100 <- top100[order(top100$Rank), , drop = FALSE]
rownames(top100) <- NULL

if (anyDuplicated(top100$Rank)) {
  stop("Duplicated Rank values detected.")
}


############################
# 3. Preserve raw and clean IDs
############################

top100$GeneID_raw <- top100$GeneID

top100$ENSEMBL_clean <- sub(
  "\\.[0-9]+$",
  "",
  top100$GeneID_raw
)

if (any(
  !grepl("^ENSG[0-9]+$", top100$ENSEMBL_clean)
)) {
  warning(
    "Some GeneID values do not follow the expected Ensembl gene-ID pattern."
  )
}

unique_ids <- unique(top100$ENSEMBL_clean)

if (length(unique_ids) == 0) {
  stop("No valid Ensembl IDs available for annotation.")
}


############################
# 4. Define 10 priority IDs
############################

priority_ids <- c(
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

annotation_ids <- unique(
  c(unique_ids, priority_ids)
)


############################
# 5. Inspect available annotation fields
############################

orgdb <- org.Hs.eg.db::org.Hs.eg.db

available_columns <- AnnotationDbi::columns(orgdb)

requested_columns <- c(
  "SYMBOL",
  "GENENAME",
  "ENTREZID"
)

usable_columns <- intersect(
  requested_columns,
  available_columns
)

if (length(usable_columns) == 0) {
  stop("No requested annotation columns are available.")
}

available_keytypes <- AnnotationDbi::keytypes(orgdb)

if (!"ENSEMBL" %in% available_keytypes) {
  stop("ENSEMBL is not a supported key type in org.Hs.eg.db.")
}


############################
# 6. Query org.Hs.eg.db
############################

valid_ensembl_keys <- AnnotationDbi::keys(
  orgdb,
  keytype = "ENSEMBL"
)

query_ids <- intersect(
  annotation_ids,
  valid_ensembl_keys
)

annotation_raw <- data.frame()

if (length(query_ids) > 0) {
  
  annotation_raw <- AnnotationDbi::select(
    orgdb,
    keys = query_ids,
    keytype = "ENSEMBL",
    columns = usable_columns
  )
  
}

if (nrow(annotation_raw) > 0) {
  
  annotation_raw$ENSEMBL <- as.character(
    annotation_raw$ENSEMBL
  )
  
  annotation_raw <- unique(annotation_raw)
  
} else {
  
  annotation_raw <- data.frame(
    ENSEMBL = character(),
    stringsAsFactors = FALSE
  )
  
}


############################
# 7. Save complete raw mapping
############################

write.table(
  annotation_raw,
  file.path(
    output_dir,
    "Workflow-04-2B-2_orgHsEgDB_raw_mapping.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)


############################
# 8. Create one-row-per-ID summary
############################

collapse_values <- function(x) {
  
  x <- unique(as.character(x))
  x <- x[!is.na(x) & trimws(x) != ""]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(sort(x), collapse = ";")
  
}

annotation_summary <- data.frame(
  ENSEMBL = annotation_ids,
  stringsAsFactors = FALSE
)

for (column_name in usable_columns) {
  
  annotation_summary[[column_name]] <- vapply(
    annotation_ids,
    function(id) {
      
      matched <- annotation_raw[
        annotation_raw$ENSEMBL == id,
        ,
        drop = FALSE
      ]
      
      if (nrow(matched) == 0 ||
          !column_name %in% colnames(matched)) {
        return(NA_character_)
      }
      
      collapse_values(matched[[column_name]])
      
    },
    character(1)
  )
  
}

mapping_count <- table(
  factor(
    annotation_raw$ENSEMBL,
    levels = annotation_ids
  )
)

annotation_summary$MappingRows <- as.integer(
  mapping_count[annotation_ids]
)

annotation_summary$AnnotationStatus <- ifelse(
  annotation_summary$MappingRows == 0,
  "No org.Hs.eg.db mapping",
  ifelse(
    annotation_summary$MappingRows > 1,
    "Multiple mapping rows",
    "One mapping row"
  )
)

annotation_summary$MitoPattern <- grepl(
  "^ENSG000001988|^ENSG000002100",
  annotation_summary$ENSEMBL
)

annotation_summary$MitoPatternInterpretation <-
  "Preliminary Ensembl-ID pattern only; not a formal biotype annotation"

write.table(
  annotation_summary,
  file.path(
    output_dir,
    "Workflow-04-2B-2_annotation_summary_by_ENSEMBL.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)


############################
# 9. Annotate Top 100 observations
############################

top100_annotated <- top100

match_idx <- match(
  top100_annotated$ENSEMBL_clean,
  annotation_summary$ENSEMBL
)

top100_annotated$orgDb_MappingRows <-
  annotation_summary$MappingRows[match_idx]

top100_annotated$orgDb_AnnotationStatus <-
  annotation_summary$AnnotationStatus[match_idx]

top100_annotated$MitoPattern <-
  annotation_summary$MitoPattern[match_idx]

for (column_name in usable_columns) {
  
  top100_annotated[[paste0("orgDb_", column_name)]] <-
    annotation_summary[[column_name]][match_idx]
  
}

top100_annotated <- top100_annotated[
  order(top100_annotated$Rank),
  ,
  drop = FALSE
]

rownames(top100_annotated) <- NULL

write.table(
  top100_annotated,
  file.path(
    output_dir,
    "Workflow-04-2B-2_top100_annotated.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)


############################
# 10. Audit all 10 priority IDs
############################

priority_audit <- annotation_summary[
  match(priority_ids, annotation_summary$ENSEMBL),
  ,
  drop = FALSE
]

priority_audit$InTop100 <- priority_audit$ENSEMBL %in% unique_ids

priority_audit$RawID_in_Top100 <- vapply(
  priority_audit$ENSEMBL,
  function(id) {
    
    x <- unique(top100$GeneID_raw[
      top100$ENSEMBL_clean == id
    ])
    
    if (length(x) == 0) {
      return(NA_character_)
    }
    
    paste(x, collapse = ";")
    
  },
  character(1)
)

write.table(
  priority_audit,
  file.path(
    output_dir,
    "Workflow-04-2B-2_priority_10_ID_audit.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)


############################
# 11. Annotation completeness summary
############################

annotation_audit <- data.frame(
  Metric = c(
    "Top100 observations",
    "Unique Ensembl IDs in Top100",
    "Total IDs submitted for annotation",
    "IDs with at least one org.Hs.eg.db mapping",
    "IDs without org.Hs.eg.db mapping",
    "IDs with multiple mapping rows",
    "Priority IDs",
    "Priority IDs found in Top100",
    "Priority IDs with at least one mapping"
  ),
  Value = c(
    nrow(top100),
    length(unique_ids),
    length(annotation_ids),
    sum(annotation_summary$MappingRows > 0),
    sum(annotation_summary$MappingRows == 0),
    sum(annotation_summary$MappingRows > 1),
    length(priority_ids),
    sum(priority_audit$InTop100),
    sum(priority_audit$MappingRows > 0)
  ),
  stringsAsFactors = FALSE
)

write.table(
  annotation_audit,
  file.path(
    output_dir,
    "Workflow-04-2B-2_annotation_audit_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 12. Save workflow report
############################

report_file <- file.path(
  output_dir,
  "Workflow-04-2B-2_workflow_report.txt"
)

report_lines <- c(
  "Workflow-04-2B-2: Annotation of extreme-expression features",
  "",
  paste("Input file:", input_file),
  paste("Input observations:", nrow(top100)),
  paste("Unique Top100 Ensembl IDs:", length(unique_ids)),
  paste("Total annotation IDs:", length(annotation_ids)),
  "",
  "Annotation source:",
  "org.Hs.eg.db Bioconductor annotation package",
  paste(
    "org.Hs.eg.db version:",
    as.character(utils::packageVersion("org.Hs.eg.db"))
  ),
  paste(
    "AnnotationDbi version:",
    as.character(utils::packageVersion("AnnotationDbi"))
  ),
  "",
  "Interpretation:",
  "org.Hs.eg.db mappings are supplementary gene-symbol/name/Entrez annotations.",
  "They are not treated as definitive GENCODE gene-biotype classifications.",
  "MitoPattern is a preliminary Ensembl-ID pattern only.",
  "",
  "No samples removed.",
  "No genes removed.",
  "No expression transformation.",
  "No normalization.",
  "No WGCNA."
)

writeLines(report_lines, report_file)


############################
# 13. Completion
############################

cat("\n========================================\n")
cat("Workflow-04-2B-2 completed.\n")
cat("Top 100 observations:", nrow(top100), "\n")
cat("Unique Top100 Ensembl IDs:", length(unique_ids), "\n")
cat("Priority IDs audited:", length(priority_ids), "\n")
cat("No samples removed.\n")
cat("No genes removed.\n")
cat("No expression transformation.\n")
cat("No normalization.\n")
cat("No WGCNA.\n")
cat("========================================\n")


############################################################
# Workflow-04-2B-2A
# GENCODE v22 annotation provenance and ID verification
#
# Research purpose:
#   Verify Ensembl gene IDs against the official GENCODE
#   human release 22 gene annotation (GRCh38.p2).
#
# Input:
#   1. Workflow-04-2B-2 annotation summary
#   2. Local GENCODE v22 comprehensive GTF file
#
# Input annotation file:
#   data/reference/GENCODE_v22/
#   gencode.v22.annotation.gtf.gz
#
# Code path:
#   scripts/02_tumor_only_WGCNA/
#
# Output:
#   results/02_tumor_only_WGCNA/Workflow-04-2B-2A/
#
# Software environment:
#   R 4.5.1
#   here
#
# Reproducibility:
#   No random procedure; random seed not applicable.
#   GTF file name, genome release, file size and MD5 checksum
#   are recorded.
#
# Annotation:
#   GENCODE human release 22, GRCh38.p2.
#
# Important:
#   - No samples removed.
#   - No genes removed from the expression matrix.
#   - No expression transformation.
#   - No normalization.
#   - No WGCNA.
############################################################


############################
# 0. Initialize environment
############################

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install it before running this script.")
}

project_dir <- here::here()

input_file <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2B-2",
  "Workflow-04-2B-2_annotation_summary_by_ENSEMBL.txt"
)

gtf_file <- file.path(
  project_dir,
  "data",
  "reference",
  "GENCODE_v22",
  "gencode.v22.annotation.gtf.gz"
)

output_dir <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2B-2A"
)

if (!file.exists(input_file)) {
  stop("Input annotation summary not found: ", input_file)
}

if (!file.exists(gtf_file)) {
  stop(
    "GENCODE v22 GTF file not found: ",
    gtf_file,
    "\nPlease download the official GENCODE v22 GTF file ",
    "and place it in data/reference/GENCODE_v22/."
  )
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


############################
# 1. Record software and file provenance
############################

software_info <- data.frame(
  Item = c(
    "R_version",
    "here_version",
    "GENCODE_release",
    "Genome_assembly",
    "GTF_filename",
    "GTF_size_bytes",
    "GTF_MD5"
  ),
  Value = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("here")),
    "GENCODE human release 22",
    "GRCh38.p2",
    basename(gtf_file),
    as.character(file.info(gtf_file)$size),
    unname(tools::md5sum(gtf_file))
  ),
  stringsAsFactors = FALSE
)

write.table(
  software_info,
  file.path(
    output_dir,
    "Workflow-04-2B-2A_annotation_provenance.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 2. Read IDs from Workflow-04-2B-2
############################

input_summary <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"ENSEMBL" %in% colnames(input_summary)) {
  stop("Column 'ENSEMBL' not found in Workflow-04-2B-2 input.")
}

top100_ids <- unique(
  input_summary$ENSEMBL[
    !is.na(input_summary$ENSEMBL) &
      trimws(input_summary$ENSEMBL) != ""
  ]
)


############################
# 3. Define the 10 priority IDs
############################

priority_ids <- c(
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

target_ids <- unique(
  c(top100_ids, priority_ids)
)

if (length(target_ids) == 0) {
  stop("No Ensembl IDs available for GENCODE v22 verification.")
}

if (any(!grepl("^ENSG[0-9]+$", target_ids))) {
  warning(
    "Some target IDs do not match the expected Ensembl gene-ID format."
  )
}


############################
# 4. Read GENCODE v22 GTF
############################

gtf <- read.delim(
  gzfile(gtf_file, open = "rt"),
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "#",
  stringsAsFactors = FALSE,
  colClasses = "character",
  check.names = FALSE
)

if (ncol(gtf) != 9) {
  stop(
    "Unexpected GTF structure: expected 9 columns, found ",
    ncol(gtf)
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


############################
# 5. Keep gene-level records
############################

gene_records <- gtf[
  gtf$feature == "gene",
  ,
  drop = FALSE
]

if (nrow(gene_records) == 0) {
  stop("No gene-level records found in the GENCODE GTF file.")
}


############################
# 6. Extract GTF attributes
############################

extract_gtf_attribute <- function(attribute_text, key) {
  
  pattern <- paste0(
    '(?:^|;[[:space:]]*)',
    key,
    '[[:space:]]+"([^"]*)"'
  )
  
  matched <- regmatches(
    attribute_text,
    regexec(pattern, attribute_text, perl = TRUE)
  )
  
  vapply(
    matched,
    function(x) {
      if (length(x) >= 2) x[2] else NA_character_
    },
    character(1)
  )
  
}


gene_records$gene_id_raw <- extract_gtf_attribute(
  gene_records$attribute,
  "gene_id"
)

gene_records$gene_name <- extract_gtf_attribute(
  gene_records$attribute,
  "gene_name"
)

gene_records$gene_type <- extract_gtf_attribute(
  gene_records$attribute,
  "gene_type"
)

gene_records$gene_biotype <- extract_gtf_attribute(
  gene_records$attribute,
  "gene_biotype"
)

gene_records$gene_status <- extract_gtf_attribute(
  gene_records$attribute,
  "gene_status"
)

gene_records$ENSEMBL_clean <- sub(
  "\\.[0-9]+$",
  "",
  gene_records$gene_id_raw
)

if (all(is.na(gene_records$gene_type)) &&
    all(is.na(gene_records$gene_biotype))) {
  stop(
    "Neither gene_type nor gene_biotype was found in GTF attributes."
  )
}

gene_records$GENE_BIOTYPE <- ifelse(
  !is.na(gene_records$gene_type) &
    gene_records$gene_type != "",
  gene_records$gene_type,
  gene_records$gene_biotype
)


############################
# 7. Select target gene records
############################

gencode_matches <- gene_records[
  gene_records$ENSEMBL_clean %in% target_ids,
  ,
  drop = FALSE
]

annotation_columns <- c(
  "gene_id_raw",
  "ENSEMBL_clean",
  "gene_name",
  "GENE_BIOTYPE",
  "gene_status",
  "seqname",
  "start",
  "end",
  "strand",
  "source"
)

gencode_annotation_raw <- gencode_matches[
  ,
  annotation_columns,
  drop = FALSE
]

colnames(gencode_annotation_raw) <- c(
  "GENCODE_GENE_ID_RAW",
  "ENSEMBL",
  "SYMBOL",
  "GENE_BIOTYPE",
  "GENE_STATUS",
  "CHROMOSOME",
  "START",
  "END",
  "STRAND",
  "ANNOTATION_SOURCE"
)

write.table(
  gencode_annotation_raw,
  file.path(
    output_dir,
    "Workflow-04-2B-2A_GENCODE_v22_raw_gene_mapping.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)


############################
# 8. Create one-row-per-ID audit
############################

collapse_values <- function(x) {
  
  x <- unique(as.character(x))
  x <- x[!is.na(x) & trimws(x) != ""]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(sort(x), collapse = ";")
  
}

gencode_summary <- data.frame(
  ENSEMBL = target_ids,
  stringsAsFactors = FALSE
)

summary_columns <- c(
  "SYMBOL",
  "GENE_BIOTYPE",
  "GENE_STATUS",
  "CHROMOSOME",
  "START",
  "END",
  "STRAND"
)

for (column_name in summary_columns) {
  
  gencode_summary[[column_name]] <- vapply(
    target_ids,
    function(id) {
      
      matched <- gencode_annotation_raw[
        gencode_annotation_raw$ENSEMBL == id,
        ,
        drop = FALSE
      ]
      
      if (nrow(matched) == 0) {
        return(NA_character_)
      }
      
      collapse_values(matched[[column_name]])
      
    },
    character(1)
  )
  
}


############################
# 9. Classify biotypes conservatively
############################

biotype <- gencode_summary$GENE_BIOTYPE

gencode_summary$BiotypeCategory <- ifelse(
  is.na(biotype),
  "Not annotated in GENCODE v22",
  ifelse(
    biotype == "protein_coding",
    "Protein-coding",
    ifelse(
      biotype == "Mt_rRNA",
      "Mitochondrial rRNA",
      ifelse(
        biotype == "Mt_tRNA",
        "Mitochondrial tRNA",
        ifelse(
          biotype == "rRNA",
          "Ribosomal RNA",
          ifelse(
            biotype == "lncRNA",
            "Long non-coding RNA",
            "Other annotated biotype"
          )
        )
      )
    )
  )
)

gencode_summary$MitochondrialChromosome <- grepl(
  "(^|;)MT($|;)",
  ifelse(
    is.na(gencode_summary$CHROMOSOME),
    "",
    gencode_summary$CHROMOSOME
  )
)

gencode_summary$MappingRows <- vapply(
  target_ids,
  function(id) {
    sum(gencode_annotation_raw$ENSEMBL == id)
  },
  integer(1)
)

gencode_summary$AnnotationStatus <- ifelse(
  gencode_summary$MappingRows == 0,
  "Not found in GENCODE v22",
  ifelse(
    gencode_summary$MappingRows > 1,
    "Multiple gene records",
    "One gene record"
  )
)

write.table(
  gencode_summary,
  file.path(
    output_dir,
    "Workflow-04-2B-2A_GENCODE_v22_annotation_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)


############################
# 10. Audit the 10 priority IDs
############################

priority_audit <- gencode_summary[
  match(priority_ids, gencode_summary$ENSEMBL),
  ,
  drop = FALSE
]

priority_audit$InTop100 <- priority_audit$ENSEMBL %in% top100_ids

write.table(
  priority_audit,
  file.path(
    output_dir,
    "Workflow-04-2B-2A_priority_10_ID_GENCODE_v22_audit.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)


############################
# 11. Identify unannotated IDs
############################

unannotated_ids <- gencode_summary[
  gencode_summary$MappingRows == 0,
  ,
  drop = FALSE
]

write.table(
  unannotated_ids,
  file.path(
    output_dir,
    "Workflow-04-2B-2A_GENCODE_v22_unannotated_IDs.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)


############################
# 12. Annotation completeness summary
############################

audit_summary <- data.frame(
  Metric = c(
    "Unique Top100 Ensembl IDs",
    "Priority IDs",
    "Total target IDs",
    "IDs found in GENCODE v22",
    "IDs not found in GENCODE v22",
    "IDs with multiple gene records",
    "Priority IDs found in Top100",
    "Priority IDs annotated in GENCODE v22"
  ),
  Value = c(
    length(top100_ids),
    length(priority_ids),
    length(target_ids),
    sum(gencode_summary$MappingRows > 0),
    sum(gencode_summary$MappingRows == 0),
    sum(gencode_summary$MappingRows > 1),
    sum(priority_audit$InTop100),
    sum(priority_audit$MappingRows > 0)
  ),
  stringsAsFactors = FALSE
)

write.table(
  audit_summary,
  file.path(
    output_dir,
    "Workflow-04-2B-2A_GENCODE_v22_audit_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 13. Save workflow report
############################

report_lines <- c(
  "Workflow-04-2B-2A",
  "GENCODE v22 annotation provenance and ID verification",
  "",
  "Annotation source: GENCODE human release 22",
  "Genome assembly: GRCh38.p2",
  paste("GTF file:", basename(gtf_file)),
  paste("GTF MD5:", unname(tools::md5sum(gtf_file))),
  "",
  paste("Unique Top100 IDs:", length(top100_ids)),
  paste("Priority IDs:", length(priority_ids)),
  paste("Total target IDs:", length(target_ids)),
  paste(
    "IDs found:",
    sum(gencode_summary$MappingRows > 0)
  ),
  paste(
    "IDs not found:",
    sum(gencode_summary$MappingRows == 0)
  ),
  "",
  "Biotype categories are based on the gene_type/gene_biotype",
  "attributes in the supplied GENCODE v22 GTF file.",
  "Mitochondrial chromosome location is recorded separately.",
  "",
  "No samples removed.",
  "No genes removed from the expression matrix.",
  "No expression transformation.",
  "No normalization.",
  "No WGCNA."
)

writeLines(
  report_lines,
  file.path(
    output_dir,
    "Workflow-04-2B-2A_workflow_report.txt"
  )
)


############################
# 14. Completion
############################

cat("\n========================================\n")
cat("Workflow-04-2B-2A completed.\n")
cat("GENCODE release: 22\n")
cat("Genome assembly: GRCh38.p2\n")
cat("Total target IDs:", length(target_ids), "\n")
cat(
  "Annotated IDs:",
  sum(gencode_summary$MappingRows > 0),
  "\n"
)
cat(
  "Unannotated IDs:",
  sum(gencode_summary$MappingRows == 0),
  "\n"
)
cat("No samples removed.\n")
cat("No genes removed from the expression matrix.\n")
cat("No expression transformation.\n")
cat("No normalization.\n")
cat("No WGCNA.\n")
cat("========================================\n")


############################################################
# Workflow-04-2C
# Expression prevalence filtering and log2 transformation
#
# Research purpose:
#   Identify genes with sufficient expression prevalence
#   across the 513 primary LUAD tumor samples and prepare
#   the retained expression matrix for downstream WGCNA.
#
# Input:
#   results/02_tumor_only_WGCNA/Workflow-04-2A/
#   Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt
#
# Code path:
#   scripts/02_tumor_only_WGCNA/
#
# Output:
#   results/02_tumor_only_WGCNA/Workflow-04-2C/
#
# Software environment:
#   R 4.5.1
#   Package: here
#   Other operations use base R.
#
# Random seed:
#   Not applicable; no stochastic procedure is used.
#
# Sample inclusion:
#   All 513 TCGA primary LUAD tumor samples.
#   TCGA sample type code = "01" (barcode positions 14-15).
#   No sample is removed during this workflow.
#
# Gene filtering:
#   Raw FPKM-UQ > 1 in at least 20% of the 513 samples.
#   Minimum required sample count = ceiling(513 * 0.20) = 103.
#
# Transformation:
#   log2(FPKM-UQ + 1), applied only after prevalence filtering.
#
# Important:
#   - Filtering is performed on the original FPKM-UQ scale.
#   - The inequality is strictly greater than 1.
#   - Gene annotation is not performed in this workflow.
#   - No gene is removed based on biotype or annotation status.
#   - No PCA, clustering, or WGCNA is performed.
############################################################


############################################################
# 0. Initialize environment and paths
############################################################

if (!requireNamespace("here", quietly = TRUE)) {
  stop(
    "Package 'here' is required. ",
    "Install it before running this script."
  )
}

project_dir <- here::here()

input_file <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2A",
  "Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt"
)

output_dir <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2C"
)

if (!file.exists(input_file)) {
  stop(
    "Input expression file not found:\n",
    input_file
  )
}

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


############################################################
# 1. Record software environment
############################################################

software_versions <- data.frame(
  Package = c(
    "R",
    "here"
  ),
  Version = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("here"))
  ),
  stringsAsFactors = FALSE
)

write.table(
  software_versions,
  file = file.path(
    output_dir,
    "Workflow-04-2C_software_versions.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 2. Define prespecified filtering parameters
############################################################

expected_samples <- 513L

expression_cutoff <- 1

prevalence_cutoff <- 0.20

minimum_samples_required <- ceiling(
  expected_samples * prevalence_cutoff
)

if (minimum_samples_required != 103L) {
  stop(
    "Unexpected minimum sample threshold. ",
    "Expected 103 samples for 513 samples at 20% prevalence."
  )
}

cat("\n========================================\n")
cat("Workflow-04-2C\n")
cat("Expression prevalence filtering and log2 transformation\n")
cat("========================================\n\n")

cat("Expected primary tumor samples:", expected_samples, "\n")
cat("Raw expression cutoff: FPKM-UQ >", expression_cutoff, "\n")
cat("Required prevalence:", prevalence_cutoff * 100, "%\n")
cat(
  "Minimum samples required:",
  minimum_samples_required,
  "\n\n"
)


############################################################
# 3. Read original FPKM-UQ expression matrix
############################################################

cat("----------------------------------------\n")
cat("PART A: Reading original expression matrix\n")
cat("----------------------------------------\n\n")

expr_raw <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

if (nrow(expr_raw) == 0) {
  stop("Input expression matrix contains no gene rows.")
}

if (ncol(expr_raw) < 2) {
  stop(
    "Input expression matrix must contain a gene ID column ",
    "and at least one sample column."
  )
}

if (!identical(
  colnames(expr_raw)[1],
  "Tag"
)) {
  stop(
    "Unexpected first column name: ",
    colnames(expr_raw)[1],
    ". Expected 'Tag'."
  )
}

cat("Input dimensions:\n")
cat("Genes:", nrow(expr_raw), "\n")
cat("Columns:", ncol(expr_raw), "\n")
cat("Samples:", ncol(expr_raw) - 1L, "\n\n")


############################################################
# 4. Validate Ensembl gene IDs
############################################################

cat("----------------------------------------\n")
cat("PART B: Gene ID quality control\n")
cat("----------------------------------------\n\n")

gene_id <- as.character(expr_raw$Tag)

gene_id_missing <- is.na(gene_id)

gene_id_blank <- !gene_id_missing &
  trimws(gene_id) == ""

gene_id_duplicate <- duplicated(gene_id) |
  duplicated(gene_id, fromLast = TRUE)

gene_id_audit <- data.frame(
  Metric = c(
    "Total gene rows",
    "Missing gene IDs",
    "Blank gene IDs",
    "Unique gene IDs",
    "Rows with duplicated gene IDs"
  ),
  Value = c(
    length(gene_id),
    sum(gene_id_missing),
    sum(gene_id_blank),
    length(unique(gene_id[!gene_id_missing & !gene_id_blank])),
    sum(gene_id_duplicate, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write.table(
  gene_id_audit,
  file = file.path(
    output_dir,
    "Workflow-04-2C_gene_ID_audit.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (any(gene_id_missing) ||
    any(gene_id_blank) ||
    any(gene_id_duplicate)) {
  stop(
    "Gene ID QC failed. ",
    "Missing, blank, or duplicated gene IDs were detected."
  )
}

cat("Gene ID QC passed.\n")
cat("Unique gene IDs:", length(gene_id), "\n\n")


############################################################
# 5. Validate sample IDs and sample inclusion
############################################################

cat("----------------------------------------\n")
cat("PART C: Sample ID and cohort QC\n")
cat("----------------------------------------\n\n")

sample_ids <- colnames(expr_raw)[-1]

if (is.null(sample_ids) ||
    anyNA(sample_ids) ||
    any(trimws(sample_ids) == "")) {
  stop("Missing or blank sample IDs detected.")
}

if (anyDuplicated(sample_ids)) {
  stop("Duplicated sample IDs detected.")
}

if (length(sample_ids) != expected_samples) {
  stop(
    "Unexpected sample count: ",
    length(sample_ids),
    ". Expected exactly ",
    expected_samples,
    " primary tumor samples."
  )
}

if (any(!grepl("^TCGA-", sample_ids))) {
  stop(
    "One or more sample IDs do not start with 'TCGA-'."
  )
}

sample_type_code <- substr(
  sample_ids,
  14,
  15
)

if (any(sample_type_code != "01")) {
  invalid_samples <- sample_ids[
    sample_type_code != "01"
  ]
  
  stop(
    "Non-primary-tumor sample IDs detected:\n",
    paste(invalid_samples, collapse = "\n")
  )
}

sample_manifest <- data.frame(
  Sample = sample_ids,
  Sample_Type_Code = sample_type_code,
  Sample_Order = seq_along(sample_ids),
  stringsAsFactors = FALSE
)

write.table(
  sample_manifest,
  file = file.path(
    output_dir,
    "Workflow-04-2C_primary_tumor_sample_manifest.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Sample ID QC passed.\n")
cat("Primary tumor samples:", length(sample_ids), "\n")
cat("All sample type codes are 01.\n")
cat("Original sample order will be preserved.\n\n")


############################################################
# 6. Convert expression columns to numeric safely
############################################################

cat("----------------------------------------\n")
cat("PART D: Numeric conversion and expression QC\n")
cat("----------------------------------------\n\n")

expression_data <- expr_raw[
  ,
  -1,
  drop = FALSE
]

conversion_issue_list <- vector(
  "list",
  length = ncol(expression_data)
)

conversion_issue_count <- integer(
  ncol(expression_data)
)

for (j in seq_len(ncol(expression_data))) {
  
  original_values <- as.character(
    expression_data[[j]]
  )
  
  converted_values <- suppressWarnings(
    as.numeric(original_values)
  )
  
  invalid_conversion <- !is.na(original_values) &
    trimws(original_values) != "" &
    is.na(converted_values)
  
  conversion_issue_count[j] <- sum(
    invalid_conversion
  )
  
  if (any(invalid_conversion)) {
    
    conversion_issue_list[[j]] <- data.frame(
      Sample = colnames(expression_data)[j],
      Original_Value = unique(
        original_values[invalid_conversion]
      ),
      stringsAsFactors = FALSE
    )
    
  }
  
  expression_data[[j]] <- converted_values
  
}

if (sum(conversion_issue_count) > 0) {
  
  conversion_issues <- do.call(
    rbind,
    conversion_issue_list[
      !vapply(
        conversion_issue_list,
        is.null,
        logical(1)
      )
    ]
  )
  
  write.table(
    conversion_issues,
    file = file.path(
      output_dir,
      "Workflow-04-2C_numeric_conversion_issues.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  stop(
    "Non-numeric expression values were detected. ",
    "See Workflow-04-2C_numeric_conversion_issues.txt."
  )
  
}

expr_mat <- as.matrix(expression_data)

storage.mode(expr_mat) <- "double"

rownames(expr_mat) <- gene_id

if (!identical(
  colnames(expr_mat),
  sample_ids
)) {
  stop("Sample order changed during numeric conversion.")
}

if (!identical(
  rownames(expr_mat),
  gene_id
)) {
  stop("Gene ID order changed during numeric conversion.")
}


############################################################
# 7. Check missing, infinite, and negative values
############################################################

n_na <- sum(is.na(expr_mat))

n_nan <- sum(is.nan(expr_mat))

n_inf <- sum(is.infinite(expr_mat))

n_negative <- sum(
  expr_mat < 0,
  na.rm = TRUE
)

expression_qc <- data.frame(
  Metric = c(
    "Expression matrix rows",
    "Expression matrix columns",
    "NA values",
    "NaN values",
    "Infinite values",
    "Negative values",
    "Numeric conversion issues"
  ),
  Value = c(
    nrow(expr_mat),
    ncol(expr_mat),
    n_na,
    n_nan,
    n_inf,
    n_negative,
    sum(conversion_issue_count)
  ),
  stringsAsFactors = FALSE
)

write.table(
  expression_qc,
  file = file.path(
    output_dir,
    "Workflow-04-2C_expression_integrity_audit.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (n_na > 0 ||
    n_nan > 0 ||
    n_inf > 0 ||
    n_negative > 0) {
  stop(
    "Expression integrity QC failed. ",
    "Missing, non-finite, or negative values were detected."
  )
}

if (any(!is.finite(expr_mat))) {
  stop("Non-finite expression values remain in the matrix.")
}

cat("Expression integrity QC passed.\n")
cat("NA:", n_na, "\n")
cat("NaN:", n_nan, "\n")
cat("Inf:", n_inf, "\n")
cat("Negative values:", n_negative, "\n\n")


############################################################
# 8. Record original matrix dimensions
############################################################

n_genes_original <- nrow(expr_mat)

n_samples <- ncol(expr_mat)

if (n_samples != expected_samples) {
  stop("The validated matrix does not contain 513 samples.")
}


############################################################
# 9. Calculate prevalence on raw FPKM-UQ scale
############################################################

cat("----------------------------------------\n")
cat("PART E: Raw FPKM-UQ prevalence calculation\n")
cat("----------------------------------------\n\n")

# Strict inequality: FPKM-UQ > 1
gene_detected_count <- rowSums(
  expr_mat > expression_cutoff
)

gene_prevalence <- gene_detected_count /
  n_samples

gene_filtering_table <- data.frame(
  ENSEMBL = rownames(expr_mat),
  Detected_Samples = as.integer(gene_detected_count),
  Total_Samples = n_samples,
  Prevalence = as.numeric(gene_prevalence),
  Median_FPKM_UQ = apply(
    expr_mat,
    1,
    stats::median
  ),
  Mean_FPKM_UQ = rowMeans(expr_mat),
  stringsAsFactors = FALSE
)

gene_filtering_table$Keep <-
  gene_filtering_table$Detected_Samples >=
  minimum_samples_required

gene_filtering_table$Filtering_Reason <- ifelse(
  gene_filtering_table$Keep,
  "Retained: prevalence threshold met",
  "Excluded: prevalence threshold not met"
)

if (any(
  gene_filtering_table$Detected_Samples < 0 |
  gene_filtering_table$Detected_Samples > n_samples
)) {
  stop("Invalid detected-sample counts.")
}

if (any(
  gene_filtering_table$Prevalence < 0 |
  gene_filtering_table$Prevalence > 1
)) {
  stop("Invalid gene prevalence values.")
}


############################################################
# 10. Determine retained genes
############################################################

genes_keep <- gene_filtering_table$ENSEMBL[
  gene_filtering_table$Keep
]

n_genes_kept <- length(genes_keep)

n_genes_removed <- n_genes_original -
  n_genes_kept

if (n_genes_kept == 0) {
  stop(
    "No genes passed the prespecified prevalence filter."
  )
}

retention_rate <- 100 *
  n_genes_kept /
  n_genes_original

cat("Original genes:", n_genes_original, "\n")
cat("Retained genes:", n_genes_kept, "\n")
cat("Removed genes:", n_genes_removed, "\n")
cat(
  "Retention rate:",
  round(retention_rate, 2),
  "%\n\n"
)


############################################################
# 11. Save complete gene filtering statistics
############################################################

write.table(
  gene_filtering_table,
  file = file.path(
    output_dir,
    "Workflow-04-2C_gene_prevalence_filtering_statistics.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

write.table(
  data.frame(
    ENSEMBL = genes_keep,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    output_dir,
    "Workflow-04-2C_retained_gene_IDs.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 12. Create prevalence-filtered raw FPKM-UQ matrix
############################################################

cat("----------------------------------------\n")
cat("PART F: Extracting prevalence-filtered raw matrix\n")
cat("----------------------------------------\n\n")

expr_raw_filtered <- expr_mat[
  genes_keep,
  ,
  drop = FALSE
]

if (nrow(expr_raw_filtered) != n_genes_kept ||
    ncol(expr_raw_filtered) != expected_samples) {
  stop("Unexpected dimensions in raw filtered expression matrix.")
}

if (!identical(
  colnames(expr_raw_filtered),
  sample_ids
)) {
  stop("Sample IDs or sample order changed during filtering.")
}

if (!identical(
  rownames(expr_raw_filtered),
  genes_keep
)) {
  stop("Retained gene order changed during filtering.")
}


############################################################
# 13. Apply log2 transformation after filtering
############################################################

cat("----------------------------------------\n")
cat("PART G: log2(FPKM-UQ + 1) transformation\n")
cat("----------------------------------------\n\n")

expr_log2_filtered <- log2(
  expr_raw_filtered + 1
)

if (any(!is.finite(expr_log2_filtered))) {
  stop(
    "Non-finite values detected after log2 transformation."
  )
}

if (!identical(
  dim(expr_log2_filtered),
  dim(expr_raw_filtered)
)) {
  stop("Matrix dimensions changed during log2 transformation.")
}

if (!identical(
  dimnames(expr_log2_filtered),
  dimnames(expr_raw_filtered)
)) {
  stop("Gene or sample names changed during transformation.")
}

cat(
  "Transformation applied to retained genes only:\n",
  "log2(FPKM-UQ + 1)\n\n",
  sep = ""
)


############################################################
# 14. Save filtered raw FPKM-UQ matrix
############################################################

raw_filtered_output <- file.path(
  output_dir,
  "Workflow-04-2C_expr_primary_tumor_FPKM_UQ_prevalence_filtered.txt"
)

raw_filtered_df <- data.frame(
  Tag = rownames(expr_raw_filtered),
  as.data.frame(
    expr_raw_filtered,
    check.names = FALSE
  ),
  check.names = FALSE
)

write.table(
  raw_filtered_df,
  file = raw_filtered_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)


############################################################
# 15. Save filtered log2 matrix: genes x samples
############################################################

log2_gene_by_sample_output <- file.path(
  output_dir,
  "Workflow-04-2C_expr_primary_tumor_log2_FPKM_UQ_prevalence_filtered.txt"
)

log2_gene_by_sample_df <- data.frame(
  Tag = rownames(expr_log2_filtered),
  as.data.frame(
    expr_log2_filtered,
    check.names = FALSE
  ),
  check.names = FALSE
)

write.table(
  log2_gene_by_sample_df,
  file = log2_gene_by_sample_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)


############################################################
# 16. Create WGCNA input matrix: samples x genes
############################################################

cat("----------------------------------------\n")
cat("PART H: Creating WGCNA candidate matrix\n")
cat("----------------------------------------\n\n")

expr_wgcna_samples_by_genes <- t(
  expr_log2_filtered
)

if (nrow(expr_wgcna_samples_by_genes) != expected_samples ||
    ncol(expr_wgcna_samples_by_genes) != n_genes_kept) {
  stop(
    "Unexpected dimensions in samples-by-genes WGCNA matrix."
  )
}

if (!identical(
  rownames(expr_wgcna_samples_by_genes),
  sample_ids
)) {
  stop("Sample order changed in WGCNA candidate matrix.")
}

if (!identical(
  colnames(expr_wgcna_samples_by_genes),
  genes_keep
)) {
  stop("Gene order changed in WGCNA candidate matrix.")
}


############################################################
# 17. Save WGCNA candidate matrix
############################################################

wgcna_candidate_output <- file.path(
  output_dir,
  "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ_samples_by_genes.txt"
)

wgcna_candidate_df <- data.frame(
  Sample = rownames(expr_wgcna_samples_by_genes),
  as.data.frame(
    expr_wgcna_samples_by_genes,
    check.names = FALSE
  ),
  check.names = FALSE
)

write.table(
  wgcna_candidate_df,
  file = wgcna_candidate_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)


############################################################
# 18. Save gene-level prevalence distribution
############################################################

prevalence_distribution <- data.frame(
  Metric = c(
    "Minimum detected samples",
    "Maximum detected samples",
    "Median detected samples",
    "Minimum prevalence",
    "Maximum prevalence",
    "Median prevalence"
  ),
  Value = c(
    min(gene_filtering_table$Detected_Samples),
    max(gene_filtering_table$Detected_Samples),
    stats::median(gene_filtering_table$Detected_Samples),
    min(gene_filtering_table$Prevalence),
    max(gene_filtering_table$Prevalence),
    stats::median(gene_filtering_table$Prevalence)
  ),
  stringsAsFactors = FALSE
)

write.table(
  prevalence_distribution,
  file = file.path(
    output_dir,
    "Workflow-04-2C_prevalence_distribution_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 19. Final matrix and filtering verification
############################################################

cat("----------------------------------------\n")
cat("PART I: Final verification\n")
cat("----------------------------------------\n\n")

final_verification <- data.frame(
  Check = c(
    "Original sample count is 513",
    "All sample IDs are unique",
    "All sample type codes are 01",
    "Original sample order preserved",
    "Minimum required sample count is 103",
    "All retained genes meet >=103 samples with FPKM-UQ > 1",
    "All excluded genes have <103 samples with FPKM-UQ > 1",
    "Retained raw matrix has 513 samples",
    "Filtered log2 matrix has 513 samples",
    "WGCNA matrix has 513 sample rows",
    "WGCNA matrix has retained genes as columns",
    "WGCNA matrix contains only finite values"
  ),
  PASS = c(
    n_samples == expected_samples,
    !anyDuplicated(sample_ids),
    all(sample_type_code == "01"),
    identical(colnames(expr_mat), sample_ids),
    minimum_samples_required == 103L,
    all(
      gene_filtering_table$Detected_Samples[
        gene_filtering_table$Keep
      ] >= minimum_samples_required
    ),
    all(
      gene_filtering_table$Detected_Samples[
        !gene_filtering_table$Keep
      ] < minimum_samples_required
    ),
    ncol(expr_raw_filtered) == expected_samples,
    ncol(expr_log2_filtered) == expected_samples,
    nrow(expr_wgcna_samples_by_genes) == expected_samples,
    ncol(expr_wgcna_samples_by_genes) == n_genes_kept,
    all(is.finite(expr_wgcna_samples_by_genes))
  ),
  stringsAsFactors = FALSE
)

write.table(
  final_verification,
  file = file.path(
    output_dir,
    "Workflow-04-2C_final_verification.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (!all(final_verification$PASS)) {
  stop(
    "At least one final verification check failed. ",
    "See Workflow-04-2C_final_verification.txt."
  )
}


############################################################
# 20. Save concise workflow summary
############################################################

workflow_summary <- data.frame(
  Item = c(
    "Input genes",
    "Input samples",
    "Sample type",
    "Expression scale used for filtering",
    "Expression cutoff",
    "Inequality",
    "Required prevalence",
    "Minimum samples required",
    "Genes retained",
    "Genes removed",
    "Retention rate (%)",
    "Filtered raw matrix dimensions",
    "Filtered log2 matrix dimensions (genes x samples)",
    "WGCNA matrix dimensions (samples x genes)",
    "Transformation",
    "Random seed",
    "Gene biotype filtering"
  ),
  Value = c(
    n_genes_original,
    n_samples,
    "TCGA primary tumor; sample type 01",
    "Original FPKM-UQ",
    expression_cutoff,
    "Strictly greater than",
    paste0(prevalence_cutoff * 100, "%"),
    minimum_samples_required,
    n_genes_kept,
    n_genes_removed,
    round(retention_rate, 2),
    paste(
      nrow(expr_raw_filtered),
      "x",
      ncol(expr_raw_filtered)
    ),
    paste(
      nrow(expr_log2_filtered),
      "x",
      ncol(expr_log2_filtered)
    ),
    paste(
      nrow(expr_wgcna_samples_by_genes),
      "x",
      ncol(expr_wgcna_samples_by_genes)
    ),
    "log2(FPKM-UQ + 1), after filtering",
    "Not applicable",
    "Not applied"
  ),
  stringsAsFactors = FALSE
)

write.table(
  workflow_summary,
  file = file.path(
    output_dir,
    "Workflow-04-2C_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 21. Save workflow report
############################################################

report_lines <- c(
  "Workflow-04-2C",
  "Expression prevalence filtering and log2 transformation",
  "",
  paste("Input file:", input_file),
  paste("Input genes:", n_genes_original),
  paste("Input samples:", n_samples),
  "Sample inclusion: TCGA primary tumor sample type 01",
  "",
  "Filtering was performed on the original FPKM-UQ scale.",
  "A gene was retained if FPKM-UQ was strictly greater than 1",
  "in at least 20% of the 513 primary tumor samples.",
  paste(
    "Minimum required sample count:",
    minimum_samples_required
  ),
  "",
  "The log2(FPKM-UQ + 1) transformation was applied",
  "only after prevalence filtering.",
  "",
  "No gene annotation was performed in this workflow.",
  "No genes were removed based on gene biotype.",
  "No samples were removed.",
  "No PCA or clustering was performed.",
  "No WGCNA was performed.",
  "",
  "The downstream WGCNA candidate matrix is stored",
  "with samples as rows and genes as columns."
)

writeLines(
  report_lines,
  file.path(
    output_dir,
    "Workflow-04-2C_workflow_report.txt"
  )
)


############################################################
# 22. Completion
############################################################

cat("\n========================================\n")
cat("Workflow-04-2C completed successfully.\n")
cat("========================================\n\n")

cat("Original genes:", n_genes_original, "\n")
cat("Primary tumor samples:", n_samples, "\n")
cat("Minimum samples required:", minimum_samples_required, "\n")
cat("Retained genes:", n_genes_kept, "\n")
cat("Removed genes:", n_genes_removed, "\n")
cat("Retention rate:", round(retention_rate, 2), "%\n\n")

cat(
  "Filtered log2 matrix (genes x samples):",
  nrow(expr_log2_filtered),
  "x",
  ncol(expr_log2_filtered),
  "\n"
)

cat(
  "WGCNA candidate matrix (samples x genes):",
  nrow(expr_wgcna_samples_by_genes),
  "x",
  ncol(expr_wgcna_samples_by_genes),
  "\n\n"
)

cat("All final verification checks passed.\n")
cat("No samples removed.\n")
cat("No gene-biotype filtering performed.\n")
cat("No PCA, clustering, or WGCNA performed.\n")
cat("========================================\n")

############################################################
# Workflow-04-2C-1
# GENCODE v22 Annotation Provenance Recovery and Cross-source Audit
#
# Purpose:
#   1. Recover gene-level annotation from the fixed GENCODE v22 GTF.
#   2. Audit the original Ensembl gene-ID universe in the source
#      FPKM-UQ expression file.
#   3. Extract the prevalence-filtered WGCNA candidate gene IDs
#      from Workflow-04-2C.
#   4. Cross-check BioMart and GENCODE v22 gene biotypes.
#   5. Identify unmatched IDs, missing biotypes and annotation
#      disagreements.
#   6. Produce annotation provenance tables for downstream review.
#
# Input data:
#   A. Workflow-04-2C:
#      results/02_tumor_only_WGCNA/Workflow-04-2C/
#      Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ_samples_by_genes.txt
#
#   B. Original source expression:
#      data/Merge_RNA_seq_FPKM_UQ.txt
#
#   C. Existing BioMart annotation table:
#      Set the exact local path in the configuration section below.
#
#   D. GENCODE v22:
#      data/reference/gencode.v22.annotation.gtf.gz
#
# Code location:
#   scripts/02_tumor_only_WGCNA/
#   Workflow-04-2C-1_GENCODE_v22_annotation_provenance.R
#
# Software environment:
#   R 4.5.1
#   Packages: here, data.table
#   Package versions are recorded in the output report.
#
# Random seed:
#   Not applicable; no stochastic procedure is used.
#
# Sample inclusion:
#   The candidate expression matrix must contain the 513 primary
#   LUAD tumor samples retained by Workflow-04-2A/04-2C.
#   Sample IDs and their original order are audited.
#
# Key analysis rules:
#   - Ensembl version suffixes are removed only for ID matching.
#   - Original IDs are retained in the audit tables.
#   - Only GTF feature == "gene" records are used.
#   - Duplicate annotation records are audited before collapsing.
#   - Conflicting annotation records are not silently resolved.
#   - No expression samples or genes are removed.
#   - No expression values are modified or transformed.
#   - No PCA, clustering or WGCNA is performed.
#   - Annotation status does not determine WGCNA gene inclusion.
############################################################


############################
# 0. Project paths and configuration
############################

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Please install it before running this script.")
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required. Please install it before running this script.")
}

project_root <- here::here()

candidate_file <- file.path(
  project_root,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2C",
  "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ_samples_by_genes.txt"
)

source_expression_file <- file.path(
  project_root,
  "data",
  "Merge_RNA_seq_FPKM_UQ.txt"
)

gencode_file <- file.path(
  project_root,
  "data",
  "reference",
  "gencode.v22.annotation.gtf.gz"
)

# IMPORTANT:
# Replace this path with the exact location and filename of the
# BioMart annotation table generated in the preceding annotation
# workflow. Do not rename or regenerate that file merely to match
# this example.
biomart_file <- file.path(
  project_root,
  "data",
  "reference",
  "biomart_gene_annotation.txt"
)

output_dir <- file.path(
  project_root,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2C-1"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################
# 1. Software and input checks
############################

cat("\n====================================================\n")
cat("Workflow-04-2C-1\n")
cat("GENCODE v22 annotation provenance recovery\n")
cat("====================================================\n\n")

software_info <- c(
  paste0("R version: ", R.version.string),
  paste0("here version: ", as.character(utils::packageVersion("here"))),
  paste0(
    "data.table version: ",
    as.character(utils::packageVersion("data.table"))
  ),
  paste0("Project root: ", project_root),
  "Random seed: Not applicable"
)

writeLines(
  software_info,
  file.path(output_dir, "Workflow-04-2C-1_software_versions.txt")
)

input_files <- c(
  candidate_expression = candidate_file,
  original_source_expression = source_expression_file,
  BioMart_annotation = biomart_file,
  GENCODE_v22_GTF = gencode_file
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "The following required input files were not found:\n",
      paste(names(missing_files), missing_files, sep = ": ", collapse = "\n"),
      "\n\nPlease correct the corresponding paths in Section 0."
    )
  )
}

cat("All four input files were found.\n\n")


############################
# 2. Helper functions
############################

normalize_ensembl_id <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- NA_character_
  
  # Remove version suffix for matching only.
  sub("\\.[0-9]+$", "", x)
}

find_column <- function(data, candidates, label) {
  matched <- candidates[candidates %in% colnames(data)]
  
  if (length(matched) == 0) {
    stop(
      paste0(
        "Could not identify the ", label, " column.\n",
        "Available columns: ",
        paste(colnames(data), collapse = ", ")
      )
    )
  }
  
  matched[1]
}

extract_gtf_attribute <- function(x, key) {
  pattern <- paste0(
    "(?:^|;[[:space:]]*)",
    key,
    "[[:space:]]+\"([^\"]*)\""
  )
  
  result <- regmatches(
    x,
    regexec(pattern, x, perl = TRUE)
  )
  
  vapply(
    result,
    function(z) {
      if (length(z) >= 2) z[2] else NA_character_
    },
    character(1)
  )
}

write_table <- function(x, filename) {
  utils::write.table(
    x,
    file = file.path(output_dir, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}


############################
# 3. Read Workflow-04-2C candidate matrix header
############################

cat("Reading Workflow-04-2C candidate matrix header...\n")

candidate_header <- data.table::fread(
  candidate_file,
  nrows = 0,
  data.table = FALSE,
  check.names = FALSE
)

candidate_columns <- colnames(candidate_header)

if (length(candidate_columns) < 2) {
  stop("Candidate matrix must contain a Sample column and at least one gene column.")
}

if (candidate_columns[1] != "Sample") {
  stop(
    paste0(
      "Unexpected first column in candidate matrix: ",
      candidate_columns[1],
      ". Expected 'Sample'."
    )
  )
}

sample_ids <- candidate_columns[1]
candidate_gene_ids_original <- candidate_columns[-1]
candidate_gene_ids <- normalize_ensembl_id(candidate_gene_ids_original)

if (anyNA(candidate_gene_ids)) {
  stop("Missing or empty gene IDs detected in candidate matrix column names.")
}

if (anyDuplicated(candidate_gene_ids)) {
  duplicated_ids <- unique(
    candidate_gene_ids[
      duplicated(candidate_gene_ids) |
        duplicated(candidate_gene_ids, fromLast = TRUE)
    ]
  )
  
  stop(
    paste0(
      "Duplicated gene IDs detected after removal of version suffixes:\n",
      paste(head(duplicated_ids, 20), collapse = ", ")
    )
  )
}

if (length(candidate_gene_ids) == 0) {
  stop("No candidate gene IDs were detected.")
}

candidate_sample_table <- data.table::fread(
  candidate_file,
  select = 1,
  data.table = FALSE,
  check.names = FALSE
)

if (nrow(candidate_sample_table) != 513) {
  stop(
    paste0(
      "Expected 513 primary tumor samples, but found ",
      nrow(candidate_sample_table),
      "."
    )
  )
}

if (anyNA(candidate_matrix$Sample) ||
    any(trimws(as.character(candidate_matrix$Sample)) == "")) {
  stop("Missing or empty sample IDs detected in the candidate matrix.")
}

if (anyDuplicated(candidate_matrix$Sample)) {
  stop("Duplicated sample IDs detected in the candidate matrix.")
}

candidate_sample_ids <- as.character(candidate_matrix$Sample)

if (!all(grepl("^TCGA-", candidate_sample_ids))) {
  stop("Some sample IDs do not have the expected TCGA barcode prefix.")
}

sample_type <- substr(candidate_sample_ids, 14, 15)

if (any(sample_type != "01")) {
  stop("Some candidate matrix samples are not TCGA sample type 01.")
}

sample_manifest <- data.frame(
  Row_order = seq_along(candidate_sample_ids),
  Sample = candidate_sample_ids,
  Sample_type = sample_type,
  stringsAsFactors = FALSE
)

write_table(
  sample_manifest,
  "Workflow-04-2C-1_candidate_sample_manifest.txt"
)

cat("Candidate matrix sample count:", nrow(candidate_matrix), "\n")
cat("Candidate gene count:", length(candidate_gene_ids), "\n\n")


############################
# 4. Audit original source Ensembl gene IDs
############################

cat("Reading original source expression gene IDs...\n")

source_header <- data.table::fread(
  source_expression_file,
  nrows = 0,
  data.table = FALSE,
  check.names = FALSE
)

source_id_column <- colnames(source_header)[1]

source_ids_table <- data.table::fread(
  source_expression_file,
  select = 1,
  data.table = FALSE,
  check.names = FALSE
)

source_gene_ids_original <- as.character(
  source_ids_table[[1]]
)

source_gene_ids <- normalize_ensembl_id(
  source_gene_ids_original
)

if (anyNA(source_gene_ids)) {
  stop("Missing or empty gene IDs detected in the original source expression file.")
}

source_id_audit <- data.frame(
  Original_gene_ID = source_gene_ids_original,
  Normalized_Ensembl_ID = source_gene_ids,
  Ensembl_format_valid = grepl("^ENSG[0-9]+$", source_gene_ids),
  stringsAsFactors = FALSE
)

if (anyDuplicated(source_gene_ids)) {
  duplicated_source_ids <- unique(
    source_gene_ids[
      duplicated(source_gene_ids) |
        duplicated(source_gene_ids, fromLast = TRUE)
    ]
  )
  
  write_table(
    data.frame(Ensembl_ID = duplicated_source_ids),
    "Workflow-04-2C-1_original_source_duplicate_gene_IDs.txt"
  )
  
  stop(
    paste0(
      "Duplicate normalized Ensembl IDs were detected in the original ",
      "source expression file. Resolve this before proceeding."
    )
  )
}

candidate_in_source <- candidate_gene_ids %in% source_gene_ids

if (!all(candidate_in_source)) {
  missing_from_source <- candidate_gene_ids[!candidate_in_source]
  
  write_table(
    data.frame(Ensembl_ID = missing_from_source),
    "Workflow-04-2C-1_candidate_IDs_not_found_in_original_source.txt"
  )
  
  stop(
    paste0(
      sum(!candidate_in_source),
      " candidate gene IDs were not found in the original source expression file."
    )
  )
}

write_table(
  source_id_audit,
  "Workflow-04-2C-1_original_source_gene_ID_audit.txt"
)

cat("Original source gene count:", length(source_gene_ids), "\n")
cat(
  "Candidate IDs recovered in original source:",
  sum(candidate_in_source),
  "/",
  length(candidate_gene_ids),
  "\n\n"
)


############################
# 5. Read BioMart annotation
############################

cat("Reading existing BioMart annotation...\n")

biomart <- data.table::fread(
  biomart_file,
  data.table = FALSE,
  check.names = FALSE
)

biomart_id_col <- find_column(
  biomart,
  c(
    "ensembl_gene_id",
    "GeneID",
    "gene_id",
    "Ensembl",
    "ENSEMBL",
    "Tag",
    "Ensembl_ID"
  ),
  "BioMart Ensembl gene ID"
)

biomart_biotype_col <- find_column(
  biomart,
  c(
    "gene_biotype",
    "Gene_Biotype",
    "biotype",
    "Biotype",
    "gene_type"
  ),
  "BioMart gene biotype"
)

biomart_name_col <- intersect(
  c("external_gene_name", "gene_name", "Gene_Name", "SYMBOL"),
  colnames(biomart)
)

biomart_clean <- data.frame(
  Ensembl_ID = normalize_ensembl_id(biomart[[biomart_id_col]]),
  BioMart_gene_biotype = trimws(as.character(biomart[[biomart_biotype_col]])),
  stringsAsFactors = FALSE
)

if (length(biomart_name_col) > 0) {
  biomart_clean$BioMart_gene_name <- as.character(
    biomart[[biomart_name_col[1]]]
  )
} else {
  biomart_clean$BioMart_gene_name <- NA_character_
}

biomart_clean$BioMart_gene_biotype[
  is.na(biomart_clean$BioMart_gene_biotype) |
    biomart_clean$BioMart_gene_biotype == ""
] <- NA_character_

biomart_clean <- biomart_clean[
  !is.na(biomart_clean$Ensembl_ID),
  ,
  drop = FALSE
]

# Audit duplicate BioMart IDs before reducing exact duplicate records.
biomart_duplicate_ids <- unique(
  biomart_clean$Ensembl_ID[
    duplicated(biomart_clean$Ensembl_ID) |
      duplicated(biomart_clean$Ensembl_ID, fromLast = TRUE)
  ]
)

if (length(biomart_duplicate_ids) > 0) {
  duplicate_biomart <- biomart_clean[
    biomart_clean$Ensembl_ID %in% biomart_duplicate_ids,
    ,
    drop = FALSE
  ]
  
  write_table(
    duplicate_biomart,
    "Workflow-04-2C-1_BioMart_duplicate_ID_records.txt"
  )
  
  # Multiple distinct biotypes for the same Ensembl ID are ambiguous.
  biotype_counts <- tapply(
    biomart_clean$BioMart_gene_biotype,
    biomart_clean$Ensembl_ID,
    function(x) length(unique(stats::na.omit(x)))
  )
  
  conflicting_ids <- names(
    biotype_counts[biotype_counts > 1]
  )
  
  if (length(conflicting_ids) > 0) {
    stop(
      paste0(
        "Conflicting BioMart biotypes were detected for ",
        length(conflicting_ids),
        " Ensembl IDs. Review Workflow-04-2C-1_BioMart_duplicate_ID_records.txt."
      )
    )
  }
}

# Collapse only records that do not conflict in biotype.
biomart_clean <- biomart_clean[
  order(biomart_clean$Ensembl_ID),
  ,
  drop = FALSE
]

biomart_clean <- biomart_clean[
  !duplicated(biomart_clean$Ensembl_ID),
  ,
  drop = FALSE
]

cat("BioMart unique Ensembl IDs:", nrow(biomart_clean), "\n\n")


############################
# 6. Read GENCODE v22 gene-level annotation
############################

cat("Reading GENCODE v22 GTF...\n")

gtf <- data.table::fread(
  gencode_file,
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  select = c(1, 3, 9),
  colClasses = "character",
  data.table = FALSE,
  quote = ""
)

if (ncol(gtf) != 3) {
  stop("GENCODE GTF parsing did not return the expected three selected columns.")
}

colnames(gtf) <- c(
  "seqname",
  "feature",
  "attribute"
)

gtf_gene <- gtf[
  gtf$feature == "gene",
  ,
  drop = FALSE
]

if (nrow(gtf_gene) == 0) {
  stop("No gene-level records were detected in the GENCODE v22 GTF.")
}

gtf_annotation <- data.frame(
  Ensembl_ID = normalize_ensembl_id(
    extract_gtf_attribute(gtf_gene$attribute, "gene_id")
  ),
  GENCODE_gene_name = extract_gtf_attribute(
    gtf_gene$attribute,
    "gene_name"
  ),
  GENCODE_gene_type = extract_gtf_attribute(
    gtf_gene$attribute,
    "gene_type"
  ),
  GENCODE_seqname = gtf_gene$seqname,
  stringsAsFactors = FALSE
)

gtf_annotation <- gtf_annotation[
  !is.na(gtf_annotation$Ensembl_ID),
  ,
  drop = FALSE
]

gtf_duplicate_ids <- unique(
  gtf_annotation$Ensembl_ID[
    duplicated(gtf_annotation$Ensembl_ID) |
      duplicated(gtf_annotation$Ensembl_ID, fromLast = TRUE)
  ]
)

if (length(gtf_duplicate_ids) > 0) {
  duplicate_gtf <- gtf_annotation[
    gtf_annotation$Ensembl_ID %in% gtf_duplicate_ids,
    ,
    drop = FALSE
  ]
  
  write_table(
    duplicate_gtf,
    "Workflow-04-2C-1_GENCODE_v22_duplicate_ID_records.txt"
  )
  
  gtf_biotype_counts <- tapply(
    gtf_annotation$GENCODE_gene_type,
    gtf_annotation$Ensembl_ID,
    function(x) length(unique(stats::na.omit(x)))
  )
  
  conflicting_gtf_ids <- names(
    gtf_biotype_counts[gtf_biotype_counts > 1]
  )
  
  if (length(conflicting_gtf_ids) > 0) {
    stop(
      paste0(
        "Conflicting GENCODE gene_type records were detected for ",
        length(conflicting_gtf_ids),
        " Ensembl IDs. Review the duplicate-record audit."
      )
    )
  }
}

gtf_annotation <- gtf_annotation[
  order(gtf_annotation$Ensembl_ID),
  ,
  drop = FALSE
]

gtf_annotation <- gtf_annotation[
  !duplicated(gtf_annotation$Ensembl_ID),
  ,
  drop = FALSE
]

cat("GENCODE v22 unique gene IDs:", nrow(gtf_annotation), "\n\n")


############################
# 7. Build full original-source annotation crosswalk
############################

cat("Building annotation crosswalk for the original source gene universe...\n")

source_crosswalk <- data.frame(
  Original_gene_ID = source_gene_ids_original,
  Ensembl_ID = source_gene_ids,
  stringsAsFactors = FALSE
)

source_crosswalk <- merge(
  source_crosswalk,
  biomart_clean,
  by = "Ensembl_ID",
  all.x = TRUE,
  sort = FALSE
)

source_crosswalk <- merge(
  source_crosswalk,
  gtf_annotation,
  by = "Ensembl_ID",
  all.x = TRUE,
  sort = FALSE
)

source_crosswalk$BioMart_biotype_available <- !is.na(
  source_crosswalk$BioMart_gene_biotype
)

source_crosswalk$GENCODE_biotype_available <- !is.na(
  source_crosswalk$GENCODE_gene_type
)

source_crosswalk$Biotype_comparison <- ifelse(
  !source_crosswalk$BioMart_biotype_available &
    !source_crosswalk$GENCODE_biotype_available,
  "Missing_in_both",
  ifelse(
    !source_crosswalk$BioMart_biotype_available &
      source_crosswalk$GENCODE_biotype_available,
    "GENCODE_only",
    ifelse(
      source_crosswalk$BioMart_biotype_available &
        !source_crosswalk$GENCODE_biotype_available,
      "BioMart_only",
      ifelse(
        source_crosswalk$BioMart_gene_biotype ==
          source_crosswalk$GENCODE_gene_type,
        "Agreement",
        "Disagreement"
      )
    )
  )
)

source_crosswalk <- source_crosswalk[
  match(source_gene_ids, source_crosswalk$Ensembl_ID),
  ,
  drop = FALSE
]

write_table(
  source_crosswalk,
  "Workflow-04-2C-1_original_source_BioMart_GENCODE_crosswalk.txt"
)


############################
# 8. Build candidate-gene annotation table
############################

cat("Building annotation table for Workflow-04-2C candidate genes...\n")

candidate_crosswalk <- data.frame(
  Candidate_order = seq_along(candidate_gene_ids),
  Original_gene_ID = candidate_gene_ids_original,
  Ensembl_ID = candidate_gene_ids,
  stringsAsFactors = FALSE
)

candidate_crosswalk <- merge(
  candidate_crosswalk,
  biomart_clean,
  by = "Ensembl_ID",
  all.x = TRUE,
  sort = FALSE
)

candidate_crosswalk <- merge(
  candidate_crosswalk,
  gtf_annotation,
  by = "Ensembl_ID",
  all.x = TRUE,
  sort = FALSE
)

candidate_crosswalk <- candidate_crosswalk[
  match(candidate_gene_ids, candidate_crosswalk$Ensembl_ID),
  ,
  drop = FALSE
]

candidate_crosswalk$BioMart_biotype_available <- !is.na(
  candidate_crosswalk$BioMart_gene_biotype
)

candidate_crosswalk$GENCODE_biotype_available <- !is.na(
  candidate_crosswalk$GENCODE_gene_type
)

candidate_crosswalk$Biotype_comparison <- ifelse(
  !candidate_crosswalk$BioMart_biotype_available &
    !candidate_crosswalk$GENCODE_biotype_available,
  "Missing_in_both",
  ifelse(
    !candidate_crosswalk$BioMart_biotype_available &
      candidate_crosswalk$GENCODE_biotype_available,
    "GENCODE_only",
    ifelse(
      candidate_crosswalk$BioMart_biotype_available &
        !candidate_crosswalk$GENCODE_biotype_available,
      "BioMart_only",
      ifelse(
        candidate_crosswalk$BioMart_gene_biotype ==
          candidate_crosswalk$GENCODE_gene_type,
        "Agreement",
        "Disagreement"
      )
    )
  )
)

candidate_crosswalk$Annotation_status <- ifelse(
  candidate_crosswalk$GENCODE_biotype_available,
  "GENCODE_v22_annotated",
  "GENCODE_v22_unresolved"
)

write_table(
  candidate_crosswalk,
  "Workflow-04-2C-1_candidate_gene_annotation_crosswalk.txt"
)


############################
# 9. Candidate-level unmatched and discrepancy lists
############################

write_table(
  candidate_crosswalk[
    !candidate_crosswalk$GENCODE_biotype_available,
    c("Original_gene_ID", "Ensembl_ID", "BioMart_gene_biotype"),
    drop = FALSE
  ],
  "Workflow-04-2C-1_candidate_genes_unresolved_by_GENCODE_v22.txt"
)

write_table(
  candidate_crosswalk[
    !candidate_crosswalk$BioMart_biotype_available,
    c("Original_gene_ID", "Ensembl_ID", "GENCODE_gene_type"),
    drop = FALSE
  ],
  "Workflow-04-2C-1_candidate_genes_missing_BioMart_biotype.txt"
)

write_table(
  candidate_crosswalk[
    candidate_crosswalk$Biotype_comparison == "Disagreement",
    ,
    drop = FALSE
  ],
  "Workflow-04-2C-1_candidate_BioMart_GENCODE_biotype_disagreements.txt"
)


############################
# 10. Annotation summaries
############################

make_summary <- function(data, label) {
  categories <- c(
    "Agreement",
    "GENCODE_only",
    "BioMart_only",
    "Missing_in_both",
    "Disagreement"
  )
  
  counts <- table(
    factor(data$Biotype_comparison, levels = categories)
  )
  
  data.frame(
    Gene_universe = label,
    Biotype_comparison = names(counts),
    N_genes = as.integer(counts),
    Percentage = round(
      100 * as.integer(counts) / nrow(data),
      3
    ),
    stringsAsFactors = FALSE
  )
}

source_summary <- make_summary(
  source_crosswalk,
  "Original_source_gene_universe"
)

candidate_summary <- make_summary(
  candidate_crosswalk,
  "Workflow-04-2C_candidate_gene_universe"
)

write_table(
  source_summary,
  "Workflow-04-2C-1_original_source_annotation_summary.txt"
)

write_table(
  candidate_summary,
  "Workflow-04-2C-1_candidate_annotation_summary.txt"
)

candidate_gene_type_summary <- as.data.frame(
  table(
    candidate_crosswalk$GENCODE_gene_type,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

colnames(candidate_gene_type_summary) <- c(
  "GENCODE_gene_type",
  "N_genes"
)

candidate_gene_type_summary$Percentage <- round(
  100 * candidate_gene_type_summary$N_genes /
    nrow(candidate_crosswalk),
  3
)

write_table(
  candidate_gene_type_summary,
  "Workflow-04-2C-1_candidate_GENCODE_v22_gene_type_distribution.txt"
)


############################
# 11. Audit the BioMart NA-biotype group
############################

source_biomart_na <- source_crosswalk[
  !source_crosswalk$BioMart_biotype_available,
  ,
  drop = FALSE
]

candidate_biomart_na <- candidate_crosswalk[
  !candidate_crosswalk$BioMart_biotype_available,
  ,
  drop = FALSE
]

write_table(
  source_biomart_na,
  "Workflow-04-2C-1_original_source_BioMart_NA_biotype_resolution.txt"
)

write_table(
  candidate_biomart_na,
  "Workflow-04-2C-1_candidate_BioMart_NA_biotype_resolution.txt"
)


############################
# 12. Final integrity checks
############################

if (nrow(candidate_crosswalk) != length(candidate_gene_ids)) {
  stop("Candidate annotation table row count does not match candidate gene count.")
}

if (!identical(
  candidate_crosswalk$Ensembl_ID,
  candidate_gene_ids
)) {
  stop("Candidate gene order was not preserved in the annotation crosswalk.")
}

if (nrow(candidate_matrix) != 513) {
  stop("Candidate expression matrix sample count changed unexpectedly.")
}

if (!identical(
  as.character(candidate_matrix$Sample),
  candidate_sample_ids
)) {
  stop("Candidate expression matrix sample order changed unexpectedly.")
}

final_checks <- c(
  "Workflow-04-2C-1 final integrity checks",
  "======================================",
  "",
  paste("Candidate samples:", nrow(candidate_matrix)),
  paste("Candidate genes:", length(candidate_gene_ids)),
  paste("Original source genes:", length(source_gene_ids)),
  paste("Candidate genes found in original source:", sum(candidate_in_source)),
  paste(
    "Candidate genes annotated by GENCODE v22:",
    sum(candidate_crosswalk$GENCODE_biotype_available)
  ),
  paste(
    "Candidate genes unresolved by GENCODE v22:",
    sum(!candidate_crosswalk$GENCODE_biotype_available)
  ),
  paste(
    "Candidate genes with BioMart/GENCODE agreement:",
    sum(candidate_crosswalk$Biotype_comparison == "Agreement")
  ),
  paste(
    "Candidate genes with BioMart/GENCODE disagreement:",
    sum(candidate_crosswalk$Biotype_comparison == "Disagreement")
  ),
  "",
  "Sample count and sample order preserved: PASS",
  "Candidate gene order preserved: PASS",
  "Candidate IDs found in original source: PASS",
  "No expression values modified: PASS",
  "No samples removed: PASS",
  "No candidate genes removed: PASS",
  "No expression transformation performed: PASS",
  "No WGCNA performed: PASS"
)

writeLines(
  final_checks,
  file.path(output_dir, "Workflow-04-2C-1_final_verification.txt")
)


############################
# 13. Workflow summary
############################

summary_lines <- c(
  "Workflow-04-2C-1",
  "GENCODE v22 Annotation Provenance Recovery",
  "==========================================",
  "",
  paste("R version:", R.version.string),
  paste("here version:", as.character(utils::packageVersion("here"))),
  paste("data.table version:", as.character(utils::packageVersion("data.table"))),
  "",
  paste("Original source gene count:", length(source_gene_ids)),
  paste("Workflow-04-2C candidate gene count:", length(candidate_gene_ids)),
  paste("Candidate sample count:", nrow(candidate_matrix)),
  "",
  paste(
    "Candidate genes annotated by GENCODE v22:",
    sum(candidate_crosswalk$GENCODE_biotype_available)
  ),
  paste(
    "Candidate genes unresolved by GENCODE v22:",
    sum(!candidate_crosswalk$GENCODE_biotype_available)
  ),
  paste(
    "BioMart/GENCODE biotype agreement:",
    sum(candidate_crosswalk$Biotype_comparison == "Agreement")
  ),
  paste(
    "BioMart/GENCODE biotype disagreement:",
    sum(candidate_crosswalk$Biotype_comparison == "Disagreement")
  ),
  "",
  "This workflow performs annotation provenance recovery only.",
  "No expression samples or genes were removed.",
  "No expression values were modified.",
  "No normalization or transformation was performed.",
  "No PCA, clustering or WGCNA was performed.",
  "",
  "Annotation status must not be used as an automatic gene-exclusion rule."
)

writeLines(
  summary_lines,
  file.path(output_dir, "Workflow-04-2C-1_summary.txt")
)

cat("\n====================================================\n")
cat("Workflow-04-2C-1 completed successfully.\n")
cat("====================================================\n")
cat("Output directory:\n", output_dir, "\n\n")
cat("Please inspect the annotation summaries and unresolved-ID lists.\n")
cat("No expression matrix was modified.\n\n")


############################################################
# Workflow-04-2C-2
# Gene Variability and MAD Assessment
#
# Purpose:
#   1. Describe gene-level expression variability across the
#      513 primary LUAD tumor samples.
#   2. Calculate mean, median, variance, SD and MAD for each
#      prevalence-filtered candidate gene.
#   3. Quantify MAD = 0 and MAD > 0 gene sets.
#   4. Compare variance = 0 and MAD = 0 gene sets.
#   5. Describe GENCODE v22 gene-type composition.
#   6. Compare annotation composition between MAD = 0 and
#      MAD > 0 genes.
#   7. Generate descriptive top-100 variance and MAD tables.
#
# Input data:
#   A. Workflow-04-2C:
#      results/02_tumor_only_WGCNA/Workflow-04-2C/
#      Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ_samples_by_genes.txt
#
#   B. Workflow-04-2C-1:
#      results/02_tumor_only_WGCNA/Workflow-04-2C-1/
#      Workflow-04-2C-1_candidate_gene_annotation_crosswalk.txt
#
# Code location:
#   scripts/02_tumor_only_WGCNA/
#   Workflow-04-2C-2_gene_variability_MAD_assessment.R
#
# Software environment:
#   R 4.5.1
#   Packages: here, data.table
#   Package versions are recorded in the output report.
#
# Random seed:
#   Not applicable; no stochastic procedure is used.
#
# Sample inclusion:
#   The 513 primary LUAD tumor samples retained by
#   Workflow-04-2A and Workflow-04-2C.
#   Sample IDs and their original order are preserved.
#
# Expression data:
#   log2(FPKM-UQ + 1), after the predefined prevalence filter:
#   raw FPKM-UQ > 1 in at least 20% of the 513 samples
#   (minimum 103 samples).
#
# Key analysis rules:
#   - Input matrix orientation: samples x genes.
#   - The first column must be Sample.
#   - Gene IDs are stored in expression matrix column names.
#   - Variance and SD use the sample-based R definitions.
#   - MAD uses stats::mad(constant = 1.4826).
#   - MAD = 0 and MAD > 0 are descriptive gene sets only.
#   - No genes or samples are removed.
#   - No additional expression filtering is performed.
#   - No PCA, clustering or WGCNA is performed.
############################################################


############################
# 0. Project paths
############################

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install it before running this script.")
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required. Install it before running this script.")
}

project_root <- here::here()

expression_file <- file.path(
  project_root,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2C",
  "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ_samples_by_genes.txt"
)

annotation_file <- file.path(
  project_root,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2C-1",
  "Workflow-04-2C-1_candidate_gene_annotation_crosswalk.txt"
)

output_dir <- file.path(
  project_root,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2C-2"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################
# 1. Software information
############################

software_info <- c(
  "Workflow-04-2C-2: Gene Variability and MAD Assessment",
  paste0("R version: ", R.version.string),
  paste0("here version: ", as.character(utils::packageVersion("here"))),
  paste0(
    "data.table version: ",
    as.character(utils::packageVersion("data.table"))
  ),
  "Random seed: Not applicable",
  "MAD definition: stats::mad(constant = 1.4826)",
  "Variance definition: stats::var()",
  "SD definition: stats::sd()",
  "Expression scale: log2(FPKM-UQ + 1)"
)

writeLines(
  software_info,
  file.path(
    output_dir,
    "Workflow-04-2C-2_software_versions.txt"
  )
)

cat("\n====================================================\n")
cat("Workflow-04-2C-2\n")
cat("Gene variability and MAD assessment\n")
cat("====================================================\n\n")


############################
# 2. Check input files
############################

required_files <- c(
  expression_file,
  annotation_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Required input file(s) not found:\n",
      paste(missing_files, collapse = "\n")
    )
  )
}

cat("Expression input:\n", expression_file, "\n\n")
cat("Annotation input:\n", annotation_file, "\n\n")


############################
# 3. Read expression matrix
############################

cat("Reading candidate expression matrix...\n")

expr_df <- data.table::fread(
  expression_file,
  data.table = FALSE,
  check.names = FALSE
)

if (ncol(expr_df) < 2) {
  stop("Expression matrix must contain Sample and at least one gene column.")
}

if (colnames(expr_df)[1] != "Sample") {
  stop(
    paste0(
      "The first expression matrix column must be 'Sample'. Found: ",
      colnames(expr_df)[1]
    )
  )
}

if (nrow(expr_df) != 513) {
  stop(
    paste0(
      "Expected 513 primary tumor samples; found ",
      nrow(expr_df),
      "."
    )
  )
}

sample_ids <- as.character(expr_df$Sample)

if (anyNA(sample_ids) ||
    any(trimws(sample_ids) == "")) {
  stop("Missing or empty sample IDs detected.")
}

if (anyDuplicated(sample_ids)) {
  stop("Duplicated sample IDs detected.")
}

if (!all(grepl("^TCGA-", sample_ids))) {
  stop("Some sample IDs do not have the expected TCGA barcode prefix.")
}

sample_type <- substr(sample_ids, 14, 15)

if (any(sample_type != "01")) {
  stop("Non-primary-tumor sample type detected.")
}

gene_ids_original <- colnames(expr_df)[-1]

gene_ids <- sub(
  "\\.[0-9]+$",
  "",
  gene_ids_original
)

if (anyNA(gene_ids) ||
    any(trimws(gene_ids) == "")) {
  stop("Missing or empty gene IDs detected in expression column names.")
}

if (anyDuplicated(gene_ids)) {
  stop("Duplicated gene IDs detected after removing version suffixes.")
}

if (!all(grepl("^ENSG[0-9]+$", gene_ids))) {
  stop("Some expression column names are not standard Ensembl gene IDs.")
}

cat("Sample count:", length(sample_ids), "\n")
cat("Candidate gene count:", length(gene_ids), "\n\n")


############################
# 4. Convert expression values safely
############################

cat("Checking numeric expression values...\n")

expression_columns <- expr_df[, -1, drop = FALSE]

numeric_columns <- lapply(
  seq_len(ncol(expression_columns)),
  function(i) {
    x <- as.character(expression_columns[[i]])
    suppressWarnings(as.numeric(x))
  }
)

expr_mat <- do.call(
  cbind,
  numeric_columns
)

rownames(expr_mat) <- sample_ids
colnames(expr_mat) <- gene_ids

storage.mode(expr_mat) <- "numeric"

if (anyNA(expr_mat)) {
  stop("Missing or non-numeric expression values detected.")
}

if (any(!is.finite(expr_mat))) {
  stop("Non-finite expression values detected.")
}

if (any(expr_mat < 0)) {
  stop("Negative log2(FPKM-UQ + 1) values detected.")
}

expression_integrity <- c(
  paste("Samples:", nrow(expr_mat)),
  paste("Genes:", ncol(expr_mat)),
  paste("NA values:", sum(is.na(expr_mat))),
  paste("NaN values:", sum(is.nan(expr_mat))),
  paste("Infinite values:", sum(!is.finite(expr_mat))),
  paste("Negative values:", sum(expr_mat < 0)),
  "Sample type 01 verification: PASS",
  "Sample IDs unique: PASS",
  "Gene IDs unique: PASS"
)

writeLines(
  expression_integrity,
  file.path(
    output_dir,
    "Workflow-04-2C-2_expression_integrity_audit.txt"
  )
)

cat("Expression matrix integrity checks passed.\n\n")


############################
# 5. Read GENCODE v22 annotation crosswalk
############################

cat("Reading Workflow-04-2C-1 annotation crosswalk...\n")

annotation <- data.table::fread(
  annotation_file,
  data.table = FALSE,
  check.names = FALSE
)

required_annotation_columns <- c(
  "Ensembl_ID",
  "GENCODE_gene_name",
  "GENCODE_gene_type"
)

missing_annotation_columns <- setdiff(
  required_annotation_columns,
  colnames(annotation)
)

if (length(missing_annotation_columns) > 0) {
  stop(
    paste0(
      "Required annotation columns are missing:\n",
      paste(missing_annotation_columns, collapse = ", ")
    )
  )
}

annotation$Ensembl_ID <- sub(
  "\\.[0-9]+$",
  "",
  as.character(annotation$Ensembl_ID)
)

if (anyNA(annotation$Ensembl_ID) ||
    any(trimws(annotation$Ensembl_ID) == "")) {
  stop("Missing Ensembl IDs detected in annotation crosswalk.")
}

if (anyDuplicated(annotation$Ensembl_ID)) {
  stop(
    "Duplicated Ensembl IDs detected in the candidate annotation crosswalk."
  )
}

annotation_index <- match(
  gene_ids,
  annotation$Ensembl_ID
)

annotation_matched <- !is.na(annotation_index)

cat(
  "Candidate genes matched to annotation:",
  sum(annotation_matched),
  "/",
  length(gene_ids),
  "\n"
)

cat(
  "Candidate genes without annotation:",
  sum(!annotation_matched),
  "\n\n"
)

gene_annotation <- annotation[
  annotation_index,
  ,
  drop = FALSE
]

if (!identical(
  as.character(gene_annotation$Ensembl_ID),
  gene_ids
)) {
  stop("Annotation order does not match expression gene order.")
}

gene_name <- as.character(
  gene_annotation$GENCODE_gene_name
)

gene_type <- as.character(
  gene_annotation$GENCODE_gene_type
)

gene_name[is.na(gene_name) | gene_name == ""] <- NA_character_
gene_type[is.na(gene_type) | gene_type == ""] <- NA_character_


############################
# 6. Calculate gene-level variability
############################

cat("Calculating gene-level variability...\n")

gene_mean <- colMeans(expr_mat)

gene_median <- apply(
  expr_mat,
  2,
  stats::median
)

gene_variance <- apply(
  expr_mat,
  2,
  stats::var
)

gene_sd <- apply(
  expr_mat,
  2,
  stats::sd
)

gene_mad <- apply(
  expr_mat,
  2,
  stats::mad,
  constant = 1.4826
)

if (any(!is.finite(gene_variance)) ||
    any(!is.finite(gene_sd)) ||
    any(!is.finite(gene_mad))) {
  stop("Non-finite variability statistics detected.")
}


############################
# 7. Descriptive detection prevalence
############################

# log2(FPKM-UQ + 1) > 1 is equivalent to raw FPKM-UQ > 1.
# This is descriptive only; no additional filtering is applied.

detection_threshold_log2 <- 1

gene_detected_n <- colSums(
  expr_mat > detection_threshold_log2
)

gene_detected_fraction <- gene_detected_n /
  nrow(expr_mat)


############################
# 8. Construct gene-level QC table
############################

gene_qc <- data.frame(
  Ensembl_ID = gene_ids,
  GENCODE_gene_name = gene_name,
  GENCODE_gene_type = gene_type,
  Annotation_available = !is.na(gene_type),
  Mean = as.numeric(gene_mean),
  Median = as.numeric(gene_median),
  Variance = as.numeric(gene_variance),
  SD = as.numeric(gene_sd),
  MAD = as.numeric(gene_mad),
  Detected_N = as.integer(gene_detected_n),
  Detected_fraction = as.numeric(gene_detected_fraction),
  stringsAsFactors = FALSE
)


############################
# 9. Define broad gene categories
############################

gene_qc$Broad_category <- "Other_annotated"

gene_qc$Broad_category[
  is.na(gene_qc$GENCODE_gene_type)
] <- "Unresolved_annotation"

gene_qc$Broad_category[
  gene_qc$GENCODE_gene_type == "protein_coding"
] <- "Protein_coding"

lncRNA_types <- c(
  "lncRNA",
  "lincRNA",
  "antisense",
  "sense_intronic",
  "sense_overlapping",
  "3prime_overlapping_ncRNA",
  "processed_transcript",
  "macro_lncRNA"
)

gene_qc$Broad_category[
  gene_qc$GENCODE_gene_type %in% lncRNA_types
] <- "Long_noncoding_RNA"

small_RNA_types <- c(
  "miRNA",
  "snRNA",
  "snoRNA",
  "scaRNA",
  "sRNA"
)

gene_qc$Broad_category[
  gene_qc$GENCODE_gene_type %in% small_RNA_types
] <- "Small_RNA"

gene_qc$Broad_category[
  gene_qc$GENCODE_gene_type %in% c("Mt_tRNA", "Mt_rRNA")
] <- "Mitochondrial_RNA"

gene_qc$Broad_category[
  gene_qc$GENCODE_gene_type == "rRNA"
] <- "Ribosomal_RNA"

gene_qc$Broad_category[
  grepl(
    "pseudogene",
    gene_qc$GENCODE_gene_type,
    ignore.case = TRUE
  ) &
    !is.na(gene_qc$GENCODE_gene_type)
] <- "Pseudogene"

gene_qc$Broad_category[
  grepl(
    "^(IG_|TR_)",
    gene_qc$GENCODE_gene_type
  ) &
    !is.na(gene_qc$GENCODE_gene_type)
] <- "Immune_receptor"

gene_qc$Broad_category[
  gene_qc$GENCODE_gene_type == "TEC"
] <- "TEC"


############################
# 10. Priority annotation flags
############################

gene_qc$Mitochondrial <- (
  gene_qc$GENCODE_gene_type %in%
    c("Mt_tRNA", "Mt_rRNA")
) |
  grepl(
    "^MT-",
    gene_qc$GENCODE_gene_name,
    ignore.case = TRUE
  )

gene_qc$Mitochondrial[is.na(gene_qc$Mitochondrial)] <- FALSE

gene_qc$rRNA_related <- gene_qc$GENCODE_gene_type %in%
  c("rRNA", "Mt_rRNA")

gene_qc$Protein_coding <- (
  gene_qc$GENCODE_gene_type == "protein_coding"
)

gene_qc$Protein_coding[
  is.na(gene_qc$Protein_coding)
] <- FALSE

gene_qc$Non_protein_coding <- (
  !is.na(gene_qc$GENCODE_gene_type) &
    gene_qc$GENCODE_gene_type != "protein_coding"
)


############################
# 11. MAD and variance groups
############################

gene_qc$MAD_group <- ifelse(
  gene_qc$MAD == 0,
  "MAD = 0",
  "MAD > 0"
)

gene_qc$Variance_group <- ifelse(
  gene_qc$Variance == 0,
  "Variance = 0",
  "Variance > 0"
)

gene_qc$MAD_Variance_relationship <- ifelse(
  gene_qc$Variance == 0 & gene_qc$MAD == 0,
  "Variance = 0 and MAD = 0",
  ifelse(
    gene_qc$Variance > 0 & gene_qc$MAD == 0,
    "Variance > 0 and MAD = 0",
    ifelse(
      gene_qc$Variance > 0 & gene_qc$MAD > 0,
      "Variance > 0 and MAD > 0",
      "Other_or_unexpected"
    )
  )
)


############################
# 12. Gene-type composition
############################

gene_type_summary <- as.data.frame(
  table(
    gene_qc$GENCODE_gene_type,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

colnames(gene_type_summary) <- c(
  "GENCODE_gene_type",
  "N_genes"
)

gene_type_summary$Percentage <- round(
  100 * gene_type_summary$N_genes / nrow(gene_qc),
  3
)

gene_type_summary <- gene_type_summary[
  order(-gene_type_summary$N_genes),
  ,
  drop = FALSE
]


############################
# 13. Broad-category composition
############################

broad_category_summary <- as.data.frame(
  table(
    gene_qc$Broad_category
  ),
  stringsAsFactors = FALSE
)

colnames(broad_category_summary) <- c(
  "Broad_category",
  "N_genes"
)

broad_category_summary$Percentage <- round(
  100 * broad_category_summary$N_genes / nrow(gene_qc),
  3
)

broad_category_summary <- broad_category_summary[
  order(-broad_category_summary$N_genes),
  ,
  drop = FALSE
]


############################
# 14. Special-category summary
############################

special_category_summary <- data.frame(
  Category = c(
    "All candidate genes",
    "GENCODE annotated",
    "Unresolved annotation",
    "Protein-coding",
    "Non-protein-coding",
    "Mitochondrial",
    "rRNA-related",
    "Long non-coding RNA",
    "Small RNA",
    "Pseudogene",
    "TEC"
  ),
  N_genes = c(
    nrow(gene_qc),
    sum(gene_qc$Annotation_available),
    sum(!gene_qc$Annotation_available),
    sum(gene_qc$Protein_coding),
    sum(gene_qc$Non_protein_coding),
    sum(gene_qc$Mitochondrial),
    sum(gene_qc$rRNA_related),
    sum(gene_qc$Broad_category == "Long_noncoding_RNA"),
    sum(gene_qc$Broad_category == "Small_RNA"),
    sum(gene_qc$Broad_category == "Pseudogene"),
    sum(gene_qc$Broad_category == "TEC")
  ),
  stringsAsFactors = FALSE
)

special_category_summary$Percentage <- round(
  100 * special_category_summary$N_genes / nrow(gene_qc),
  3
)


############################
# 15. MAD = 0 and MAD > 0 gene sets
############################

mad_zero_genes <- gene_qc[
  gene_qc$MAD == 0,
  ,
  drop = FALSE
]

mad_positive_genes <- gene_qc[
  gene_qc$MAD > 0,
  ,
  drop = FALSE
]

mad_group_summary <- data.frame(
  MAD_group = c("MAD = 0", "MAD > 0"),
  N_genes = c(
    nrow(mad_zero_genes),
    nrow(mad_positive_genes)
  ),
  Percentage = round(
    100 * c(
      nrow(mad_zero_genes),
      nrow(mad_positive_genes)
    ) / nrow(gene_qc),
    3
  ),
  stringsAsFactors = FALSE
)


############################
# 16. Compare gene-type composition by MAD group
############################

mad_gene_type_table <- as.data.frame(
  table(
    MAD_group = gene_qc$MAD_group,
    GENCODE_gene_type = ifelse(
      is.na(gene_qc$GENCODE_gene_type),
      "Unresolved_annotation",
      gene_qc$GENCODE_gene_type
    )
  ),
  stringsAsFactors = FALSE
)

colnames(mad_gene_type_table) <- c(
  "MAD_group",
  "GENCODE_gene_type",
  "N_genes"
)

mad_gene_type_table$Within_group_percentage <- ave(
  mad_gene_type_table$N_genes,
  mad_gene_type_table$MAD_group,
  FUN = function(x) {
    if (sum(x) == 0) rep(0, length(x)) else 100 * x / sum(x)
  }
)

mad_gene_type_table$Within_group_percentage <- round(
  mad_gene_type_table$Within_group_percentage,
  3
)


############################
# 17. Variance = 0 versus MAD = 0
############################

variance_mad_relationship <- as.data.frame(
  table(
    Variance_group = gene_qc$Variance_group,
    MAD_group = gene_qc$MAD_group
  ),
  stringsAsFactors = FALSE
)

colnames(variance_mad_relationship) <- c(
  "Variance_group",
  "MAD_group",
  "N_genes"
)

variance_mad_relationship$Percentage <- round(
  100 * variance_mad_relationship$N_genes / nrow(gene_qc),
  3
)


############################
# 18. Variance distribution
############################

variance_probabilities <- c(
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
)

variance_quantiles <- stats::quantile(
  gene_qc$Variance,
  probs = variance_probabilities,
  na.rm = TRUE,
  names = FALSE
)

variance_summary <- data.frame(
  Quantile = paste0(
    variance_probabilities * 100,
    "%"
  ),
  Value = as.numeric(variance_quantiles),
  stringsAsFactors = FALSE
)

variance_summary <- rbind(
  variance_summary,
  data.frame(
    Quantile = c("Mean", "SD"),
    Value = c(
      mean(gene_qc$Variance),
      stats::sd(gene_qc$Variance)
    )
  )
)


############################
# 19. MAD distribution
############################

mad_quantiles <- stats::quantile(
  gene_qc$MAD,
  probs = variance_probabilities,
  na.rm = TRUE,
  names = FALSE
)

mad_summary <- data.frame(
  Quantile = paste0(
    variance_probabilities * 100,
    "%"
  ),
  Value = as.numeric(mad_quantiles),
  stringsAsFactors = FALSE
)

mad_summary <- rbind(
  mad_summary,
  data.frame(
    Quantile = c("Mean", "SD"),
    Value = c(
      mean(gene_qc$MAD),
      stats::sd(gene_qc$MAD)
    )
  )
)


############################
# 20. Low-MAD descriptive thresholds
############################

low_mad_probabilities <- c(
  0.01,
  0.05,
  0.10,
  0.20,
  0.25
)

low_mad_values <- stats::quantile(
  gene_qc$MAD,
  probs = low_mad_probabilities,
  na.rm = TRUE,
  names = FALSE
)

low_mad_summary <- data.frame(
  Percentile = paste0(
    low_mad_probabilities * 100,
    "%"
  ),
  MAD_threshold = as.numeric(low_mad_values),
  N_genes_at_or_below = vapply(
    low_mad_values,
    function(x) sum(gene_qc$MAD <= x),
    integer(1)
  ),
  stringsAsFactors = FALSE
)


############################
# 21. Top 100 variable genes
############################

top100_variance_genes <- gene_qc[
  order(
    -gene_qc$Variance,
    gene_qc$Ensembl_ID
  ),
  ,
  drop = FALSE
]

top100_variance_genes <- head(
  top100_variance_genes,
  100
)


############################
# 22. Top 100 MAD genes
############################

top100_mad_genes <- gene_qc[
  order(
    -gene_qc$MAD,
    gene_qc$Ensembl_ID
  ),
  ,
  drop = FALSE
]

top100_mad_genes <- head(
  top100_mad_genes,
  100
)


############################
# 23. Top-gene type composition
############################

top100_variance_gene_type <- as.data.frame(
  table(
    top100_variance_genes$GENCODE_gene_type,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

colnames(top100_variance_gene_type) <- c(
  "GENCODE_gene_type",
  "N_genes"
)

top100_mad_gene_type <- as.data.frame(
  table(
    top100_mad_genes$GENCODE_gene_type,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

colnames(top100_mad_gene_type) <- c(
  "GENCODE_gene_type",
  "N_genes"
)


############################
# 24. Important gene audit
############################

important_gene_symbols <- c(
  "AOC3",
  "MS4A7",
  "IL7R",
  "GAPDH",
  "ACTB",
  "MT-ND1",
  "MT-ND2",
  "MT-ND3",
  "MT-ND4",
  "MT-CO1",
  "MT-CO2",
  "MT-CO3",
  "MT-ATP6",
  "MALAT1",
  "MT-RNR2",
  "RN7SL2",
  "RN7SK",
  "S100A9",
  "IGKC",
  "IGHA1",
  "SFTPA2",
  "SCGB3A1",
  "SLPI"
)

important_gene_audit <- gene_qc[
  !is.na(gene_qc$GENCODE_gene_name) &
    gene_qc$GENCODE_gene_name %in% important_gene_symbols,
  ,
  drop = FALSE
]

important_gene_audit$Requested_symbol <- important_gene_symbols[
  match(
    important_gene_audit$GENCODE_gene_name,
    important_gene_symbols
  )
]

missing_important_symbols <- setdiff(
  important_gene_symbols,
  important_gene_audit$GENCODE_gene_name
)

missing_important_gene_table <- data.frame(
  Requested_symbol = missing_important_symbols,
  Status = "Not_present_in_candidate_gene_set",
  stringsAsFactors = FALSE
)


############################
# 25. Historical extreme-feature ID audit
############################

important_ensembl_ids <- c(
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

important_ensembl_audit <- gene_qc[
  gene_qc$Ensembl_ID %in% important_ensembl_ids,
  ,
  drop = FALSE
]

missing_important_ensembl <- setdiff(
  important_ensembl_ids,
  important_ensembl_audit$Ensembl_ID
)

missing_important_ensembl_table <- data.frame(
  Ensembl_ID = missing_important_ensembl,
  Status = "Not_present_in_candidate_gene_set",
  stringsAsFactors = FALSE
)


############################
# 26. Candidate genes: AOC3, MS4A7 and IL7R
############################

candidate_symbols <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

candidate_gene_check <- gene_qc[
  !is.na(gene_qc$GENCODE_gene_name) &
    gene_qc$GENCODE_gene_name %in% candidate_symbols,
  ,
  drop = FALSE
]

candidate_gene_check <- candidate_gene_check[
  order(
    match(
      candidate_gene_check$GENCODE_gene_name,
      candidate_symbols
    )
  ),
  ,
  drop = FALSE
]

missing_candidate_symbols <- setdiff(
  candidate_symbols,
  candidate_gene_check$GENCODE_gene_name
)

missing_candidate_gene_table <- data.frame(
  Requested_symbol = missing_candidate_symbols,
  Status = "Not_present_in_candidate_gene_set",
  stringsAsFactors = FALSE
)


############################
# 27. Overall summary
############################

overall_summary <- data.frame(
  Metric = c(
    "Expression samples",
    "Candidate genes",
    "GENCODE annotated genes",
    "GENCODE unresolved genes",
    "MAD = 0 genes",
    "MAD > 0 genes",
    "Variance = 0 genes",
    "Variance > 0 genes",
    "Variance > 0 and MAD = 0 genes",
    "Variance > 0 and MAD > 0 genes",
    "Protein-coding genes",
    "Non-protein-coding genes",
    "Mitochondrial genes",
    "rRNA-related genes"
  ),
  Value = c(
    nrow(expr_mat),
    nrow(gene_qc),
    sum(gene_qc$Annotation_available),
    sum(!gene_qc$Annotation_available),
    nrow(mad_zero_genes),
    nrow(mad_positive_genes),
    sum(gene_qc$Variance == 0),
    sum(gene_qc$Variance > 0),
    sum(
      gene_qc$Variance > 0 &
        gene_qc$MAD == 0
    ),
    sum(
      gene_qc$Variance > 0 &
        gene_qc$MAD > 0
    ),
    sum(gene_qc$Protein_coding),
    sum(gene_qc$Non_protein_coding),
    sum(gene_qc$Mitochondrial),
    sum(gene_qc$rRNA_related)
  ),
  stringsAsFactors = FALSE
)


############################
# 28. Save output tables
############################

write_output <- function(x, filename) {
  utils::write.table(
    x,
    file = file.path(output_dir, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

write_output(
  gene_qc,
  "Workflow-04-2C-2_gene_variability_QC_full.txt"
)

write_output(
  gene_type_summary,
  "Workflow-04-2C-2_GENCODE_gene_type_summary.txt"
)

write_output(
  broad_category_summary,
  "Workflow-04-2C-2_broad_gene_category_summary.txt"
)

write_output(
  special_category_summary,
  "Workflow-04-2C-2_special_category_summary.txt"
)

write_output(
  mad_group_summary,
  "Workflow-04-2C-2_MAD_group_summary.txt"
)

write_output(
  mad_zero_genes,
  "Workflow-04-2C-2_MAD_zero_gene_set.txt"
)

write_output(
  mad_positive_genes,
  "Workflow-04-2C-2_MAD_positive_gene_set.txt"
)

write_output(
  mad_gene_type_table,
  "Workflow-04-2C-2_gene_type_by_MAD_group.txt"
)

write_output(
  variance_mad_relationship,
  "Workflow-04-2C-2_variance_vs_MAD_relationship.txt"
)

write_output(
  variance_summary,
  "Workflow-04-2C-2_variance_distribution_summary.txt"
)

write_output(
  mad_summary,
  "Workflow-04-2C-2_MAD_distribution_summary.txt"
)

write_output(
  low_mad_summary,
  "Workflow-04-2C-2_low_MAD_distribution_summary.txt"
)

write_output(
  top100_variance_genes,
  "Workflow-04-2C-2_top100_variance_genes.txt"
)

write_output(
  top100_mad_genes,
  "Workflow-04-2C-2_top100_MAD_genes.txt"
)

write_output(
  top100_variance_gene_type,
  "Workflow-04-2C-2_top100_variance_gene_type_summary.txt"
)

write_output(
  top100_mad_gene_type,
  "Workflow-04-2C-2_top100_MAD_gene_type_summary.txt"
)

write_output(
  important_gene_audit,
  "Workflow-04-2C-2_important_gene_audit.txt"
)

write_output(
  missing_important_gene_table,
  "Workflow-04-2C-2_important_symbols_not_in_candidate_set.txt"
)

write_output(
  important_ensembl_audit,
  "Workflow-04-2C-2_important_Ensembl_ID_audit.txt"
)

write_output(
  missing_important_ensembl_table,
  "Workflow-04-2C-2_important_Ensembl_IDs_not_in_candidate_set.txt"
)

write_output(
  candidate_gene_check,
  "Workflow-04-2C-2_AOC3_MS4A7_IL7R_QC.txt"
)

write_output(
  overall_summary,
  "Workflow-04-2C-2_overall_summary.txt"
)

write_output(
  missing_candidate_gene_table,
  "Workflow-04-2C-2_AOC3_MS4A7_IL7R_not_in_candidate_set.txt"
)

############################
# 29. Final verification
############################

if (nrow(expr_mat) != 513) {
  stop("Final verification failed: sample count is not 513.")
}

if (ncol(expr_mat) != nrow(gene_qc)) {
  stop("Final verification failed: gene count mismatch.")
}

if (sum(gene_qc$MAD_group == "MAD = 0") +
    sum(gene_qc$MAD_group == "MAD > 0") != nrow(gene_qc)) {
  stop("Final verification failed: MAD groups do not cover all candidate genes.")
}

if (sum(gene_qc$Variance_group == "Variance = 0") +
    sum(gene_qc$Variance_group == "Variance > 0") != nrow(gene_qc)) {
  stop("Final verification failed: variance groups do not cover all candidate genes.")
}

if (any(gene_qc$Variance == 0 & gene_qc$MAD > 0)) {
  stop("Unexpected result: variance = 0 but MAD > 0.")
}

if (!identical(
  rownames(expr_mat),
  sample_ids
)) {
  stop("Sample order changed unexpectedly.")
}

if (!identical(
  colnames(expr_mat),
  gene_ids
)) {
  stop("Gene order changed unexpectedly.")
}

final_verification <- c(
  "Workflow-04-2C-2 final verification",
  "====================================",
  "",
  paste("Number of samples:", nrow(expr_mat)),
  paste("Number of candidate genes:", ncol(expr_mat)),
  paste("MAD = 0 genes:", nrow(mad_zero_genes)),
  paste("MAD > 0 genes:", nrow(mad_positive_genes)),
  paste("Variance = 0 genes:", sum(gene_qc$Variance == 0)),
  paste("Variance > 0 genes:", sum(gene_qc$Variance > 0)),
  "",
  "Sample count = 513: PASS",
  "Sample IDs unique: PASS",
  "Sample type 01: PASS",
  "Sample order preserved: PASS",
  "Gene IDs unique: PASS",
  "Gene order preserved: PASS",
  "Annotation order aligned: PASS",
  "All expression values finite: PASS",
  "MAD groups cover all genes: PASS",
  "Variance groups cover all genes: PASS",
  "Variance = 0 implies MAD = 0: PASS",
  "",
  "No samples removed.",
  "No genes removed.",
  "No expression values modified.",
  "No additional prevalence filtering performed.",
  "No PCA or clustering performed.",
  "No WGCNA performed."
)

writeLines(
  final_verification,
  file.path(
    output_dir,
    "Workflow-04-2C-2_final_verification.txt"
  )
)


############################
# 30. Workflow summary
############################

summary_lines <- c(
  "Workflow-04-2C-2",
  "Gene Variability and MAD Assessment",
  "===================================",
  "",
  paste("R version:", R.version.string),
  paste("here version:", as.character(utils::packageVersion("here"))),
  paste(
    "data.table version:",
    as.character(utils::packageVersion("data.table"))
  ),
  "Random seed: Not applicable",
  "MAD definition: stats::mad(constant = 1.4826)",
  "",
  paste("Samples:", nrow(expr_mat)),
  paste("Candidate genes:", ncol(expr_mat)),
  paste("MAD = 0 genes:", nrow(mad_zero_genes)),
  paste("MAD > 0 genes:", nrow(mad_positive_genes)),
  paste("Variance = 0 genes:", sum(gene_qc$Variance == 0)),
  paste("Variance > 0 genes:", sum(gene_qc$Variance > 0)),
  "",
  "This workflow is descriptive only.",
  "No genes or samples were removed.",
  "No expression values were modified.",
  "No additional expression filtering was performed.",
  "No PCA, clustering or WGCNA was performed."
)

writeLines(
  summary_lines,
  file.path(
    output_dir,
    "Workflow-04-2C-2_summary.txt"
  )
)

cat("\n====================================================\n")
cat("Workflow-04-2C-2 completed successfully.\n")
cat("====================================================\n\n")

cat("Samples:", nrow(expr_mat), "\n")
cat("Candidate genes:", ncol(expr_mat), "\n")
cat("MAD = 0 genes:", nrow(mad_zero_genes), "\n")
cat("MAD > 0 genes:", nrow(mad_positive_genes), "\n")
cat("Variance = 0 genes:", sum(gene_qc$Variance == 0), "\n")
cat("Variance > 0 genes:", sum(gene_qc$Variance > 0), "\n\n")

cat("No genes or samples were removed.\n")
cat("All results saved to:\n")
cat(output_dir, "\n\n")

############################################################
# END OF WORKFLOW-04-2C-2
############################################################


############################################################
# Workflow-04-2C-3
# MAD Group Composition and Candidate Gene Comparison
#
# Research purpose:
# 1. Characterize the MAD = 0 and MAD > 0 gene groups
#    in the prevalence-filtered primary LUAD gene universe.
# 2. Compare gene-type and broad-category composition
#    between the complete candidate universe and MAD groups.
# 3. Audit the MAD status of AOC3, MS4A7, and IL7R.
# 4. Extract and verify the MAD > 0 expression matrix
#    for the subsequent WGCNA workflow.
#
# Input data:
# 1. Workflow-04-2C-2 full gene variability QC table.
# 2. Workflow-04-2C samples-by-genes log2 expression matrix.
#
# Input cohort:
# 513 TCGA primary LUAD tumor samples (sample type 01).
#
# MAD definition:
# stats::mad(x, constant = 1.4826)
#
# MAD groups:
# MAD = 0: MAD exactly equal to zero.
# MAD > 0: MAD strictly greater than zero.
#
# Important:
# - No samples are removed.
# - No additional expression-based gene filtering is performed.
# - The MAD > 0 matrix is extracted only as a downstream
#   WGCNA input; the original expression matrix is unchanged.
# - MAD = 0 does not necessarily mean Variance = 0.
# - Gene annotations are descriptive and do not determine
#   inclusion in the MAD > 0 expression matrix.
#
# Code location:
# scripts/02_tumor_only_WGCNA/
#
# Output location:
# results/02_tumor_only_WGCNA/Workflow-04-2C-3/
#
# Software:
# R 4.5.1
# Package: here
#
# Random seed:
# Not applicable; no stochastic procedure is used.
############################################################


############################################################
# 0. Setup
############################################################

options(stringsAsFactors = FALSE)

required_packages <- c("here")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace,
          logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Required package(s) not installed: ",
    paste(missing_packages, collapse = ", "),
    ". Please install them manually before running this script."
  )
}

project_dir <- here::here()

input_qc_file <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2C-2",
  "Workflow-04-2C-2_gene_variability_QC_full.txt"
)

input_expression_file <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2C",
  "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ_samples_by_genes.txt"
)

output_dir <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2C-3"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!file.exists(input_qc_file)) {
  stop("QC input file not found: ", input_qc_file)
}

if (!file.exists(input_expression_file)) {
  stop(
    "Expression input file not found: ",
    input_expression_file
  )
}

write_output <- function(x, filename) {
  output_file <- file.path(output_dir, filename)
  
  if (is.data.frame(x) || is.matrix(x)) {
    utils::write.table(
      x,
      file = output_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE,
      na = "NA"
    )
  } else {
    writeLines(as.character(x), output_file)
  }
  
  invisible(output_file)
}


############################################################
# 1. Record software versions and workflow settings
############################################################

software_info <- c(
  "Workflow: Workflow-04-2C-3",
  paste("R version:", R.version.string),
  paste("here version:",
        as.character(utils::packageVersion("here"))),
  paste("Project directory:", project_dir),
  paste("QC input:", input_qc_file),
  paste("Expression input:", input_expression_file),
  paste("Output directory:", output_dir),
  "Expected primary tumor sample number: 513",
  "TCGA sample type code: 01",
  "MAD definition: stats::mad(constant = 1.4826)",
  "MAD group definitions: MAD == 0; MAD > 0",
  "Random seed: Not applicable"
)

write_output(
  software_info,
  "Workflow-04-2C-3_software_versions_and_settings.txt"
)


############################################################
# 2. Read input files
############################################################

gene_qc <- utils::read.delim(
  input_qc_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "")
)

expression_df <- utils::read.delim(
  input_expression_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "")
)

if (nrow(gene_qc) == 0L || ncol(gene_qc) == 0L) {
  stop("The gene QC table is empty.")
}

if (nrow(expression_df) == 0L || ncol(expression_df) < 2L) {
  stop("The expression matrix is empty or has an invalid structure.")
}


############################################################
# 3. Standardize and validate QC column names
############################################################

# Accept only explicitly defined alternative column names.
rename_first_match <- function(data, target, alternatives) {
  if (target %in% colnames(data)) {
    return(data)
  }
  
  matched <- alternatives[alternatives %in% colnames(data)]
  
  if (length(matched) > 0L) {
    colnames(data)[match(matched[1], colnames(data))] <- target
  }
  
  data
}

gene_qc <- rename_first_match(
  gene_qc, "Ensembl_ID",
  c("gene_id", "Gene_ID", "ENSEMBL_ID")
)

gene_qc <- rename_first_match(
  gene_qc, "gene_name",
  c("GENCODE_gene_name", "gene_symbol", "Gene_name")
)

gene_qc <- rename_first_match(
  gene_qc, "gene_type",
  c("GENCODE_gene_type", "gene_biotype", "Gene_type")
)

gene_qc <- rename_first_match(
  gene_qc, "Broad_category",
  c("broad_category", "Gene_category")
)

required_qc_columns <- c(
  "Ensembl_ID",
  "Variance",
  "MAD"
)

missing_qc_columns <- setdiff(
  required_qc_columns,
  colnames(gene_qc)
)

if (length(missing_qc_columns) > 0L) {
  stop(
    "Required columns missing from gene QC table: ",
    paste(missing_qc_columns, collapse = ", ")
  )
}

if (!"gene_name" %in% colnames(gene_qc)) {
  gene_qc$gene_name <- NA_character_
}

if (!"gene_type" %in% colnames(gene_qc)) {
  gene_qc$gene_type <- NA_character_
}

if (!"Broad_category" %in% colnames(gene_qc)) {
  gene_qc$Broad_category <- NA_character_
}


############################################################
# 4. Validate Ensembl IDs and numerical QC variables
############################################################

gene_qc$Ensembl_ID <- trimws(
  as.character(gene_qc$Ensembl_ID)
)

if (anyNA(gene_qc$Ensembl_ID) ||
    any(gene_qc$Ensembl_ID == "")) {
  stop("Missing or blank Ensembl IDs detected in gene QC table.")
}

if (anyDuplicated(gene_qc$Ensembl_ID) > 0L) {
  duplicated_ids <- unique(
    gene_qc$Ensembl_ID[
      duplicated(gene_qc$Ensembl_ID) |
        duplicated(gene_qc$Ensembl_ID, fromLast = TRUE)
    ]
  )
  
  write_output(
    data.frame(Ensembl_ID = duplicated_ids),
    "Workflow-04-2C-3_duplicate_Ensembl_IDs.txt"
  )
  
  stop("Duplicated Ensembl IDs detected in gene QC table.")
}

for (column_name in c("Variance", "MAD")) {
  original_values <- gene_qc[[column_name]]
  
  converted_values <- suppressWarnings(
    as.numeric(as.character(original_values))
  )
  
  invalid_values <- is.na(converted_values) &
    !is.na(original_values)
  
  if (any(invalid_values)) {
    stop(
      "Non-numeric values detected in QC column: ",
      column_name
    )
  }
  
  if (anyNA(converted_values) ||
      any(!is.finite(converted_values))) {
    stop(
      "Missing or non-finite values detected in QC column: ",
      column_name
    )
  }
  
  if (any(converted_values < 0)) {
    stop("Negative values detected in QC column: ", column_name)
  }
  
  gene_qc[[column_name]] <- converted_values
}

if (any(gene_qc$MAD < 0)) {
  stop("MAD cannot be negative.")
}


############################################################
# 5. Validate the expression matrix structure
############################################################

if (!identical(colnames(expression_df)[1], "Sample")) {
  stop(
    "The first column of the expression matrix must be named 'Sample'."
  )
}

sample_ids <- trimws(as.character(expression_df$Sample))

if (length(sample_ids) != 513L) {
  stop(
    "Expected 513 primary tumor samples, but found ",
    length(sample_ids), "."
  )
}

if (anyNA(sample_ids) ||
    any(sample_ids == "") ||
    anyDuplicated(sample_ids) > 0L) {
  stop("Missing, blank, or duplicated sample IDs detected.")
}

if (any(!grepl("^TCGA-", sample_ids))) {
  stop("Non-TCGA sample identifiers detected.")
}

sample_type <- substr(sample_ids, 14, 15)

if (any(sample_type != "01")) {
  stop("One or more samples are not TCGA primary tumor type 01.")
}

expression_gene_ids <- colnames(expression_df)[-1]

if (anyNA(expression_gene_ids) ||
    any(expression_gene_ids == "") ||
    anyDuplicated(expression_gene_ids) > 0L) {
  stop("Missing, blank, or duplicated gene columns in expression matrix.")
}

expression_gene_ids <- trimws(expression_gene_ids)

if (any(expression_gene_ids == "")) {
  stop("Blank gene identifiers detected in expression matrix.")
}

# Strict numeric conversion for all expression columns.
expression_matrix <- matrix(
  NA_real_,
  nrow = nrow(expression_df),
  ncol = length(expression_gene_ids),
  dimnames = list(sample_ids, expression_gene_ids)
)

numeric_issue_count <- 0L

for (j in seq_along(expression_gene_ids)) {
  raw_values <- expression_df[[j + 1L]]
  
  converted_values <- suppressWarnings(
    as.numeric(as.character(raw_values))
  )
  
  invalid <- is.na(converted_values) & !is.na(raw_values)
  
  if (any(invalid)) {
    numeric_issue_count <- numeric_issue_count +
      sum(invalid)
  }
  
  if (any(invalid) ||
      anyNA(converted_values) ||
      any(!is.finite(converted_values))) {
    stop(
      "Invalid expression values detected in gene column: ",
      expression_gene_ids[j]
    )
  }
  
  expression_matrix[, j] <- converted_values
}

if (any(expression_matrix < 0)) {
  stop("Negative values detected in the log2 expression matrix.")
}

if (nrow(expression_matrix) != 513L) {
  stop("Expression matrix does not contain exactly 513 samples.")
}


############################################################
# 6. Check gene-universe agreement
############################################################

qc_gene_ids <- gene_qc$Ensembl_ID

qc_ids_not_in_expression <- setdiff(
  qc_gene_ids,
  expression_gene_ids
)

expression_ids_not_in_qc <- setdiff(
  expression_gene_ids,
  qc_gene_ids
)

gene_universe_audit <- data.frame(
  Metric = c(
    "Genes in QC table",
    "Genes in expression matrix",
    "QC genes absent from expression matrix",
    "Expression genes absent from QC table",
    "Shared genes"
  ),
  N = c(
    length(qc_gene_ids),
    length(expression_gene_ids),
    length(qc_ids_not_in_expression),
    length(expression_ids_not_in_qc),
    length(intersect(qc_gene_ids, expression_gene_ids))
  )
)

write_output(
  gene_universe_audit,
  "Workflow-04-2C-3_gene_universe_integrity_audit.txt"
)

if (length(qc_ids_not_in_expression) > 0L ||
    length(expression_ids_not_in_qc) > 0L) {
  
  write_output(
    data.frame(Ensembl_ID = qc_ids_not_in_expression),
    "Workflow-04-2C-3_QC_genes_absent_from_expression_matrix.txt"
  )
  
  write_output(
    data.frame(Ensembl_ID = expression_ids_not_in_qc),
    "Workflow-04-2C-3_expression_genes_absent_from_QC_table.txt"
  )
  
  stop(
    "The QC gene universe and expression matrix gene universe do not match. ",
    "Please resolve the upstream Workflow-04-2C-2 outputs before proceeding."
  )
}

# The same gene universe must be represented in both files.
if (length(qc_gene_ids) != length(expression_gene_ids)) {
  stop("Gene-universe size mismatch.")
}


############################################################
# 7. Define MAD groups
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

genes_mad0 <- gene_qc[
  gene_qc$MAD_group == "MAD = 0",
  ,
  drop = FALSE
]

genes_mad_gt0 <- gene_qc[
  gene_qc$MAD_group == "MAD > 0",
  ,
  drop = FALSE
]

if (nrow(genes_mad0) + nrow(genes_mad_gt0) != nrow(gene_qc)) {
  stop("MAD groups do not partition the complete gene QC table.")
}

if (length(intersect(
  genes_mad0$Ensembl_ID,
  genes_mad_gt0$Ensembl_ID
)) > 0L) {
  stop("An Ensembl ID occurs in both MAD groups.")
}

if (length(union(
  genes_mad0$Ensembl_ID,
  genes_mad_gt0$Ensembl_ID
)) != nrow(gene_qc)) {
  stop("MAD groups do not cover the complete gene universe.")
}


############################################################
# 8. MAD group composition
############################################################

mad_summary <- data.frame(
  MAD_group = c("MAD = 0", "MAD > 0"),
  Gene_N = c(nrow(genes_mad0), nrow(genes_mad_gt0)),
  Percentage = c(
    nrow(genes_mad0) / nrow(gene_qc) * 100,
    nrow(genes_mad_gt0) / nrow(gene_qc) * 100
  )
)

mad_summary$Percentage <- round(
  mad_summary$Percentage, 2
)

write_output(
  mad_summary,
  "Workflow-04-2C-3_MAD_group_summary.txt"
)


############################################################
# 9. Safe composition summary functions
############################################################

safe_category <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "Unresolved"
  x
}

composition_table <- function(data, category_column, group_label) {
  category <- safe_category(data[[category_column]])
  
  counts <- table(category)
  
  result <- data.frame(
    Category = names(counts),
    Gene_N = as.integer(counts),
    stringsAsFactors = FALSE
  )
  
  result$Group <- group_label
  result$Percentage <- if (nrow(data) > 0L) {
    round(result$Gene_N / nrow(data) * 100, 2)
  } else {
    NA_real_
  }
  
  result
}

compare_composition <- function(data_all, data_zero, data_positive,
                                category_column) {
  
  all_tab <- composition_table(
    data_all, category_column, "All genes"
  )
  
  zero_tab <- composition_table(
    data_zero, category_column, "MAD = 0"
  )
  
  positive_tab <- composition_table(
    data_positive, category_column, "MAD > 0"
  )
  
  categories <- sort(unique(c(
    all_tab$Category,
    zero_tab$Category,
    positive_tab$Category
  )))
  
  get_count <- function(tab, category) {
    index <- match(category, tab$Category)
    value <- tab$Gene_N[index]
    value[is.na(value)] <- 0L
    value
  }
  
  get_percentage <- function(tab, category) {
    index <- match(category, tab$Category)
    value <- tab$Percentage[index]
    value[is.na(value)] <- 0
    value
  }
  
  all_n <- get_count(all_tab, categories)
  zero_n <- get_count(zero_tab, categories)
  positive_n <- get_count(positive_tab, categories)
  
  all_pct <- get_percentage(all_tab, categories)
  zero_pct <- get_percentage(zero_tab, categories)
  positive_pct <- get_percentage(positive_tab, categories)
  
  result <- data.frame(
    Category = categories,
    All_N = all_n,
    MAD0_N = zero_n,
    MADgt0_N = positive_n,
    All_Percentage = all_pct,
    MAD0_Percentage = zero_pct,
    MADgt0_Percentage = positive_pct,
    Fraction_of_category_MAD0 = ifelse(
      all_n > 0,
      round(zero_n / all_n * 100, 2),
      NA_real_
    ),
    MAD0_composition_ratio = ifelse(
      all_pct > 0,
      round(zero_pct / all_pct, 3),
      NA_real_
    ),
    stringsAsFactors = FALSE
  )
  
  result <- result[
    order(-result$MAD0_N, result$Category),
    ,
    drop = FALSE
  ]
  
  result
}


############################################################
# 10. Gene-type and broad-category comparisons
############################################################

gene_type_comparison <- compare_composition(
  gene_qc,
  genes_mad0,
  genes_mad_gt0,
  "gene_type"
)

broad_category_comparison <- compare_composition(
  gene_qc,
  genes_mad0,
  genes_mad_gt0,
  "Broad_category"
)

write_output(
  gene_type_comparison,
  "Workflow-04-2C-3_gene_type_MAD_comparison.txt"
)

write_output(
  broad_category_comparison,
  "Workflow-04-2C-3_broad_category_MAD_comparison.txt"
)


############################################################
# 11. Special-category comparison
############################################################

broad_category <- safe_category(gene_qc$Broad_category)

special_category_names <- c(
  "Protein-coding",
  "Non-protein-coding",
  "Mitochondrial RNA",
  "Ribosomal RNA",
  "Long non-coding RNA",
  "Small RNA",
  "Pseudogene",
  "TEC",
  "Unresolved"
)

count_special <- function(data, category_name) {
  category_values <- safe_category(data$Broad_category)
  
  if (category_name == "Non-protein-coding") {
    return(sum(
      category_values != "Protein-coding" &
        category_values != "Unresolved"
    ))
  }
  
  sum(category_values == category_name)
}

special_categories <- data.frame(
  Category = special_category_names,
  All_N = vapply(
    special_category_names,
    function(x) count_special(gene_qc, x),
    numeric(1)
  ),
  MAD0_N = vapply(
    special_category_names,
    function(x) count_special(genes_mad0, x),
    numeric(1)
  ),
  MADgt0_N = vapply(
    special_category_names,
    function(x) count_special(genes_mad_gt0, x),
    numeric(1)
  ),
  stringsAsFactors = FALSE
)

special_categories$All_Percentage <- round(
  special_categories$All_N / nrow(gene_qc) * 100,
  2
)

special_categories$MAD0_Percentage <- round(
  special_categories$MAD0_N / nrow(genes_mad0) * 100,
  2
)

special_categories$MADgt0_Percentage <- round(
  special_categories$MADgt0_N / nrow(genes_mad_gt0) * 100,
  2
)

write_output(
  special_categories,
  "Workflow-04-2C-3_special_category_MAD_comparison.txt"
)


############################################################
# 12. Paper 1 candidate gene audit
############################################################

paper1_genes <- c("AOC3", "MS4A7", "IL7R")

candidate_rows <- gene_qc[
  gene_qc$gene_name %in% paper1_genes,
  ,
  drop = FALSE
]

candidate_rows$Requested_symbol <- candidate_rows$gene_name

candidate_rows <- candidate_rows[
  order(
    match(candidate_rows$Requested_symbol, paper1_genes),
    candidate_rows$Ensembl_ID
  ),
  ,
  drop = FALSE
]

candidate_status_list <- lapply(
  paper1_genes,
  function(symbol) {
    
    matched <- gene_qc[
      !is.na(gene_qc$gene_name) &
        gene_qc$gene_name == symbol,
      ,
      drop = FALSE
    ]
    
    if (nrow(matched) == 0L) {
      return(data.frame(
        Gene = symbol,
        Mapping_N = 0L,
        Ensembl_IDs = NA_character_,
        MAD_groups = NA_character_,
        Variance_values = NA_character_,
        MAD_values = NA_character_,
        Status = "Not found in QC gene universe",
        stringsAsFactors = FALSE
      ))
    }
    
    data.frame(
      Gene = symbol,
      Mapping_N = nrow(matched),
      Ensembl_IDs = paste(matched$Ensembl_ID, collapse = ";"),
      MAD_groups = paste(
        as.character(matched$MAD_group),
        collapse = ";"
      ),
      Variance_values = paste(
        signif(matched$Variance, 6),
        collapse = ";"
      ),
      MAD_values = paste(
        signif(matched$MAD, 6),
        collapse = ";"
      ),
      Status = if (nrow(matched) == 1L) {
        "One Ensembl mapping"
      } else {
        "Multiple Ensembl mappings; review required"
      },
      stringsAsFactors = FALSE
    )
  }
)

candidate_status <- do.call(
  rbind,
  candidate_status_list
)

write_output(
  candidate_rows,
  "Workflow-04-2C-3_Paper1_candidate_gene_QC_rows.txt"
)

write_output(
  candidate_status,
  "Workflow-04-2C-3_Paper1_candidate_gene_MAD_status.txt"
)


############################################################
# 13. Variance and MAD summaries by group
############################################################

summarize_numeric <- function(x) {
  c(
    N = length(x),
    Median = stats::median(x),
    Mean = mean(x),
    Minimum = min(x),
    Maximum = max(x)
  )
}

variance_mad_summary <- rbind(
  "MAD = 0" = summarize_numeric(genes_mad0$Variance),
  "MAD > 0" = summarize_numeric(genes_mad_gt0$Variance)
)

variance_mad_summary <- data.frame(
  MAD_group = rownames(variance_mad_summary),
  variance_mad_summary,
  row.names = NULL,
  check.names = FALSE
)

write_output(
  variance_mad_summary,
  "Workflow-04-2C-3_variance_summary_by_MAD_group.txt"
)


############################################################
# 14. Compare complete and MAD-positive gene universes
############################################################

count_broad <- function(data, category_name) {
  values <- safe_category(data$Broad_category)
  
  if (category_name == "Non-protein-coding") {
    return(sum(
      values != "Protein-coding" &
        values != "Unresolved"
    ))
  }
  
  sum(values == category_name)
}

universe_categories <- c(
  "Protein-coding",
  "Non-protein-coding",
  "Mitochondrial RNA",
  "Ribosomal RNA"
)

universe_comparison <- data.frame(
  Category = universe_categories,
  All_N = vapply(
    universe_categories,
    function(x) count_broad(gene_qc, x),
    numeric(1)
  ),
  MADgt0_N = vapply(
    universe_categories,
    function(x) count_broad(genes_mad_gt0, x),
    numeric(1)
  ),
  stringsAsFactors = FALSE
)

universe_comparison$All_Percentage <- round(
  universe_comparison$All_N / nrow(gene_qc) * 100,
  2
)

universe_comparison$MADgt0_Percentage <- round(
  universe_comparison$MADgt0_N /
    nrow(genes_mad_gt0) * 100,
  2
)

universe_summary <- data.frame(
  Universe = c(
    "Complete prevalence-filtered gene universe",
    "MAD > 0 gene universe"
  ),
  Gene_N = c(
    nrow(gene_qc),
    nrow(genes_mad_gt0)
  ),
  Protein_coding_N = c(
    count_broad(gene_qc, "Protein-coding"),
    count_broad(genes_mad_gt0, "Protein-coding")
  ),
  Non_protein_coding_N = c(
    count_broad(gene_qc, "Non-protein-coding"),
    count_broad(genes_mad_gt0, "Non-protein-coding")
  ),
  Mitochondrial_RNA_N = c(
    count_broad(gene_qc, "Mitochondrial RNA"),
    count_broad(genes_mad_gt0, "Mitochondrial RNA")
  ),
  Ribosomal_RNA_N = c(
    count_broad(gene_qc, "Ribosomal RNA"),
    count_broad(genes_mad_gt0, "Ribosomal RNA")
  ),
  stringsAsFactors = FALSE
)

universe_summary$Protein_coding_Percentage <- round(
  universe_summary$Protein_coding_N /
    universe_summary$Gene_N * 100,
  2
)

universe_summary$Non_protein_coding_Percentage <- round(
  universe_summary$Non_protein_coding_N /
    universe_summary$Gene_N * 100,
  2
)

write_output(
  universe_comparison,
  "Workflow-04-2C-3_gene_category_universe_comparison.txt"
)

write_output(
  universe_summary,
  "Workflow-04-2C-3_gene_universe_summary.txt"
)


############################################################
# 15. Extract MAD-positive gene IDs
############################################################

mad_positive_gene_ids <- genes_mad_gt0$Ensembl_ID

if (length(mad_positive_gene_ids) == 0L) {
  stop("No MAD > 0 genes were identified.")
}

if (anyDuplicated(mad_positive_gene_ids) > 0L) {
  stop("Duplicated gene IDs in MAD-positive gene set.")
}

# Preserve the original expression matrix gene-column order.
mad_positive_gene_ids_in_expression_order <-
  expression_gene_ids[
    expression_gene_ids %in% mad_positive_gene_ids
  ]

if (length(mad_positive_gene_ids_in_expression_order) !=
    length(mad_positive_gene_ids)) {
  stop("MAD-positive gene extraction does not match the QC gene set.")
}

write_output(
  data.frame(
    Ensembl_ID = mad_positive_gene_ids_in_expression_order
  ),
  "Workflow-04-2C-3_MAD_positive_gene_IDs.txt"
)

write_output(
  data.frame(
    Ensembl_ID = genes_mad0$Ensembl_ID
  ),
  "Workflow-04-2C-3_MAD_zero_gene_IDs.txt"
)


############################################################
# 16. Extract the MAD-positive expression matrix
############################################################

mad_positive_expression <- expression_matrix[
  ,
  mad_positive_gene_ids_in_expression_order,
  drop = FALSE
]

if (nrow(mad_positive_expression) != 513L) {
  stop("MAD-positive expression matrix does not retain 513 samples.")
}

if (ncol(mad_positive_expression) !=
    length(mad_positive_gene_ids)) {
  stop("MAD-positive expression matrix has an unexpected gene count.")
}

if (!identical(
  rownames(mad_positive_expression),
  sample_ids
)) {
  stop("Sample IDs or sample order changed during MAD-positive extraction.")
}

if (!identical(
  colnames(mad_positive_expression),
  mad_positive_gene_ids_in_expression_order
)) {
  stop("Gene-column order changed during MAD-positive extraction.")
}

# Write a samples-by-genes matrix with Sample as the first column.
mad_positive_output <- data.frame(
  Sample = sample_ids,
  as.data.frame(
    mad_positive_expression,
    check.names = FALSE
  ),
  check.names = FALSE
)

write_output(
  mad_positive_output,
  "Workflow-04-2C-3_WGCNA_input_MAD_positive_samples_by_genes.txt"
)


############################################################
# 17. Verify MAD-positive matrix values against source matrix
############################################################

source_subset <- expression_matrix[
  ,
  mad_positive_gene_ids_in_expression_order,
  drop = FALSE
]

values_identical <- isTRUE(all.equal(
  unname(mad_positive_expression),
  unname(source_subset),
  tolerance = 0,
  check.attributes = FALSE
))

if (!values_identical) {
  stop("Expression values changed during MAD-positive extraction.")
}

sample_manifest <- data.frame(
  Sample = sample_ids,
  Sample_type_code = sample_type,
  Sample_order = seq_along(sample_ids),
  stringsAsFactors = FALSE
)

write_output(
  sample_manifest,
  "Workflow-04-2C-3_primary_tumor_sample_manifest.txt"
)


############################################################
# 18. Final integrity and handoff audit
############################################################

final_audit <- data.frame(
  Check = c(
    "Expected primary tumor sample number",
    "Observed sample number",
    "All sample IDs are unique",
    "All samples are TCGA primary tumor type 01",
    "Sample order preserved",
    "QC and expression gene universes match",
    "MAD groups partition all genes",
    "MAD-positive gene count",
    "MAD-positive matrix sample number",
    "MAD-positive matrix gene number",
    "MAD-positive gene order follows source matrix",
    "MAD-positive expression values unchanged"
  ),
  Result = c(
    513,
    nrow(mad_positive_expression),
    !anyDuplicated(sample_ids),
    all(sample_type == "01"),
    identical(rownames(mad_positive_expression), sample_ids),
    length(qc_ids_not_in_expression) == 0L &&
      length(expression_ids_not_in_qc) == 0L,
    nrow(genes_mad0) + nrow(genes_mad_gt0) == nrow(gene_qc),
    ncol(mad_positive_expression),
    nrow(mad_positive_expression),
    ncol(mad_positive_expression),
    identical(
      colnames(mad_positive_expression),
      mad_positive_gene_ids_in_expression_order
    ),
    values_identical
  ),
  stringsAsFactors = FALSE
)

write_output(
  final_audit,
  "Workflow-04-2C-3_final_integrity_audit.txt"
)


############################################################
# 19. Workflow summary
############################################################

summary_lines <- c(
  "Workflow-04-2C-3 Summary",
  "========================",
  "",
  paste("Total genes in QC universe:", nrow(gene_qc)),
  paste("MAD = 0 genes:", nrow(genes_mad0)),
  paste("MAD > 0 genes:", nrow(genes_mad_gt0)),
  paste(
    "MAD = 0 percentage:",
    round(nrow(genes_mad0) / nrow(gene_qc) * 100, 2)
  ),
  paste(
    "MAD > 0 percentage:",
    round(nrow(genes_mad_gt0) / nrow(gene_qc) * 100, 2)
  ),
  "",
  paste("Primary tumor samples retained:",
        nrow(mad_positive_expression)),
  paste("MAD-positive genes extracted:",
        ncol(mad_positive_expression)),
  "",
  "Paper 1 candidate genes:",
  paste(paper1_genes, collapse = ", "),
  "",
  "Candidate gene MAD status:"
)

candidate_lines <- apply(
  candidate_status,
  1,
  function(x) {
    paste(
      x["Gene"],
      " | mappings=", x["Mapping_N"],
      " | MAD group=", x["MAD_groups"],
      " | status=", x["Status"],
      sep = ""
    )
  }
)

summary_lines <- c(
  summary_lines,
  candidate_lines,
  "",
  "No samples were removed.",
  "No additional genes were filtered in this workflow.",
  "The MAD-positive matrix preserves source sample order and expression values.",
  "The MAD-positive samples-by-genes matrix is the proposed input for Workflow-04-2D."
)

write_output(
  summary_lines,
  "Workflow-04-2C-3_summary.txt"
)


############################################################
# 20. Console summary
############################################################

cat("\n")
cat("====================================================\n")
cat("Workflow-04-2C-3 completed successfully.\n")
cat("====================================================\n")

cat("Total genes:", nrow(gene_qc), "\n")
cat("MAD = 0 genes:", nrow(genes_mad0), "\n")
cat("MAD > 0 genes:", nrow(genes_mad_gt0), "\n")
cat("Samples retained:", nrow(mad_positive_expression), "\n")
cat("MAD-positive genes extracted:",
    ncol(mad_positive_expression), "\n")

cat("\nPaper 1 candidate gene status:\n")
print(candidate_status, row.names = FALSE)

cat("\nFinal integrity audit:\n")
print(final_audit, row.names = FALSE)

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("====================================================\n")



############################################################
# Workflow-04-2D
# WGCNA Pre-analysis: Candidate Matrix and Sample-level QC
#
# Research purpose:
# 1. Import the MAD-positive expression matrix established in
#    Workflow-04-2C-3.
# 2. Verify the gene universe against the Workflow-04-2C-2
#    gene variability QC table.
# 3. Verify the retention and order of all 513 TCGA primary
#    LUAD tumor samples.
# 4. Run WGCNA goodSamplesGenes() as a data-integrity QC step.
# 5. Perform PCA and hierarchical clustering for descriptive
#    sample-level QC.
# 6. Calculate sample connectivity and examine its distribution
#    across TCGA Tissue Source Sites (TSS), including TSS 44.
# 7. Export the verified samples-by-genes expression matrix
#    for downstream WGCNA network construction.
#
# Important:
# - This workflow does NOT construct a WGCNA network.
# - PCA and hierarchical clustering use a QC-only subset of
#   highly variable genes; this does not change the WGCNA
#   candidate gene universe.
# - Connectivity and TSS results are descriptive only.
# - No sample is removed based on PCA, clustering, connectivity,
#   TSS, or survival information.
# - No additional gene filtering is performed in this workflow.
#
# Input data:
# 1. Workflow-04-2C-3 MAD-positive expression matrix:
#    results/02_tumor_only_WGCNA/Workflow-04-2C-3/
#    Workflow-04-2C-3_WGCNA_input_MAD_positive_samples_by_genes.txt
#
# 2. Workflow-04-2C-2 gene variability QC table:
#    results/02_tumor_only_WGCNA/Workflow-04-2C-2/
#    Workflow-04-2C-2_gene_variability_QC_full.txt
#
# 3. Workflow-04-2A aligned clinical data:
#    results/02_tumor_only_WGCNA/Workflow-04-2A/
#    Workflow-04-2A_clinical_aligned.txt
#
# Output:
#    results/02_tumor_only_WGCNA/Workflow-04-2D/
#
# Script location:
#    scripts/02_tumor_only_WGCNA/
#
# Sample inclusion:
# - 513 TCGA primary LUAD tumor samples.
# - TCGA sample type code 01 (barcode positions 14-15).
# - All sample IDs must be unique and retained in their
#   original Workflow-04-2C order.
#
# Expression:
# - log2(FPKM-UQ + 1).
# - Input orientation: samples in rows, genes in columns.
#
# Software environment:
# - R version: 4.5.1
# - Packages: here, WGCNA, ggplot2
# - Package versions are recorded in the output directory.
#
# Random seed:
# - Not applicable. No stochastic procedure is used.
#
# Network construction:
# - Not performed in this workflow.
############################################################


############################
# 0. Initialize
############################

options(stringsAsFactors = FALSE)

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Please install it before running this script.")
}

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop("Package 'WGCNA' is required. Please install it before running this script.")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Please install it before running this script.")
}

project_dir <- here::here()

input_expr_file <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2C-3",
  "Workflow-04-2C-3_WGCNA_input_MAD_positive_samples_by_genes.txt"
)

input_gene_qc_file <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2C-2",
  "Workflow-04-2C-2_gene_variability_QC_full.txt"
)

input_clinical_file <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2A",
  "Workflow-04-2A_clinical_aligned.txt"
)

output_dir <- file.path(
  project_dir,
  "results",
  "02_tumor_only_WGCNA",
  "Workflow-04-2D"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "02_tumor_only_WGCNA",
  "Workflow-04-2D"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  input_expr_file,
  input_gene_qc_file,
  input_clinical_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following required input files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}


############################
# 1. Record software versions
############################

package_versions <- data.frame(
  Package = c("R", "here", "WGCNA", "ggplot2"),
  Version = c(
    as.character(getRversion()),
    as.character(packageVersion("here")),
    as.character(packageVersion("WGCNA")),
    as.character(packageVersion("ggplot2"))
  ),
  stringsAsFactors = FALSE
)

write.table(
  package_versions,
  file = file.path(
    output_dir,
    "Workflow-04-2D_software_versions.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========== Workflow-04-2D ==========\n")
cat("Project directory:", project_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("R version:", as.character(getRversion()), "\n")
cat("WGCNA version:", as.character(packageVersion("WGCNA")), "\n")


############################
# 2. Read MAD-positive expression matrix
############################

expr_df <- read.delim(
  input_expr_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(expr_df) < 2) {
  stop("Expression input must contain a Sample column and at least one gene column.")
}

if (colnames(expr_df)[1] != "Sample") {
  stop(
    "The first column of the MAD-positive expression matrix must be named 'Sample'."
  )
}

sample_ids <- as.character(expr_df$Sample)
gene_ids <- colnames(expr_df)[-1]

if (anyNA(sample_ids) || any(trimws(sample_ids) == "")) {
  stop("Missing or blank sample IDs detected.")
}

if (anyDuplicated(sample_ids)) {
  stop("Duplicated sample IDs detected in the expression matrix.")
}

if (length(sample_ids) != 513) {
  stop(
    "Expected 513 primary tumor samples, but found ",
    length(sample_ids),
    "."
  )
}

if (any(substr(sample_ids, 14, 15) != "01")) {
  stop("At least one expression sample does not have TCGA sample type code 01.")
}

if (anyNA(gene_ids) || any(trimws(gene_ids) == "")) {
  stop("Missing or blank gene IDs detected.")
}

if (anyDuplicated(gene_ids)) {
  stop("Duplicated Ensembl gene IDs detected in the expression matrix.")
}

if (!all(grepl("^ENSG[0-9]+(\\.[0-9]+)?$", gene_ids))) {
  stop("Unexpected gene identifier format detected.")
}


############################
# 3. Strict numeric conversion
############################

expr_char <- as.matrix(
  expr_df[, -1, drop = FALSE]
)

expr_numeric <- matrix(
  NA_real_,
  nrow = nrow(expr_char),
  ncol = ncol(expr_char),
  dimnames = list(sample_ids, gene_ids)
)

conversion_issue_count <- 0L

for (j in seq_len(ncol(expr_char))) {
  
  original_values <- as.character(expr_char[, j])
  converted_values <- suppressWarnings(as.numeric(original_values))
  
  invalid <- is.na(converted_values) & !is.na(original_values)
  
  if (any(invalid)) {
    conversion_issue_count <- conversion_issue_count + sum(invalid)
  }
  
  expr_numeric[, j] <- converted_values
}

if (conversion_issue_count > 0) {
  stop(
    "Non-numeric expression values detected during conversion. ",
    "Please inspect the input matrix."
  )
}

if (anyNA(expr_numeric) ||
    any(is.nan(expr_numeric)) ||
    any(is.infinite(expr_numeric))) {
  stop("Expression matrix contains NA, NaN, or infinite values.")
}

if (any(expr_numeric < 0)) {
  stop("Negative expression values detected.")
}

if (!all(is.finite(expr_numeric))) {
  stop("Expression matrix contains non-finite values.")
}

datExpr_pre <- expr_numeric

if (!identical(rownames(datExpr_pre), sample_ids)) {
  stop("Sample row names were not preserved.")
}

if (!identical(colnames(datExpr_pre), gene_ids)) {
  stop("Gene column names were not preserved.")
}


############################
# 4. Read and validate gene QC table
############################

gene_qc <- read.delim(
  input_gene_qc_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_qc_columns <- c(
  "Ensembl_ID",
  "MAD"
)

missing_qc_columns <- setdiff(
  required_qc_columns,
  colnames(gene_qc)
)

if (length(missing_qc_columns) > 0) {
  stop(
    "Missing required columns in gene QC table: ",
    paste(missing_qc_columns, collapse = ", ")
  )
}

gene_qc$Ensembl_ID <- as.character(gene_qc$Ensembl_ID)
gene_qc$MAD <- suppressWarnings(as.numeric(gene_qc$MAD))

if (anyNA(gene_qc$Ensembl_ID) ||
    any(trimws(gene_qc$Ensembl_ID) == "")) {
  stop("Missing or blank Ensembl IDs in the gene QC table.")
}

if (anyDuplicated(gene_qc$Ensembl_ID)) {
  stop("Duplicated Ensembl IDs in the gene QC table.")
}

if (anyNA(gene_qc$MAD) || any(!is.finite(gene_qc$MAD))) {
  stop("Missing or non-finite MAD values in the gene QC table.")
}

if (any(gene_qc$MAD < 0)) {
  stop("Negative MAD values detected in the gene QC table.")
}

mad_positive_ids <- gene_qc$Ensembl_ID[gene_qc$MAD > 0]

if (length(mad_positive_ids) == 0) {
  stop("No MAD-positive genes were found in the gene QC table.")
}


############################
# 5. Verify gene universe
############################

expr_gene_ids <- colnames(datExpr_pre)

missing_from_expression <- setdiff(
  mad_positive_ids,
  expr_gene_ids
)

extra_in_expression <- setdiff(
  expr_gene_ids,
  mad_positive_ids
)

gene_universe_audit <- data.frame(
  Metric = c(
    "Total genes in gene QC table",
    "MAD = 0 genes in gene QC table",
    "MAD > 0 genes in gene QC table",
    "Genes in expression matrix",
    "MAD-positive genes missing from expression matrix",
    "Unexpected genes in expression matrix",
    "Gene universe exactly matched"
  ),
  Value = c(
    nrow(gene_qc),
    sum(gene_qc$MAD == 0),
    sum(gene_qc$MAD > 0),
    length(expr_gene_ids),
    length(missing_from_expression),
    length(extra_in_expression),
    as.character(
      length(missing_from_expression) == 0 &&
        length(extra_in_expression) == 0
    )
  ),
  stringsAsFactors = FALSE
)

write.table(
  gene_universe_audit,
  file = file.path(
    output_dir,
    "Workflow-04-2D_gene_universe_integrity_audit.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (length(missing_from_expression) > 0 ||
    length(extra_in_expression) > 0) {
  
  writeLines(
    c(
      "MAD-positive genes missing from expression matrix:",
      missing_from_expression,
      "",
      "Unexpected genes in expression matrix:",
      extra_in_expression
    ),
    con = file.path(
      output_dir,
      "Workflow-04-2D_gene_universe_mismatch_IDs.txt"
    )
  )
  
  stop(
    "The MAD-positive gene universe does not exactly match ",
    "the expression matrix. See gene_universe_mismatch_IDs.txt."
  )
}

# Reorder QC records to exactly match expression gene order.
gene_qc_matched <- gene_qc[
  match(expr_gene_ids, gene_qc$Ensembl_ID),
  ,
  drop = FALSE
]

if (!identical(gene_qc_matched$Ensembl_ID, expr_gene_ids)) {
  stop("Gene QC table could not be aligned to the expression matrix.")
}

if (any(gene_qc_matched$MAD <= 0)) {
  stop("The expression matrix contains genes that are not MAD-positive.")
}


############################
# 6. Read and align clinical data
############################

clinical_df <- read.delim(
  input_clinical_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"Sample" %in% colnames(clinical_df)) {
  stop("Clinical input must contain a column named 'Sample'.")
}

clinical_df$Sample <- as.character(clinical_df$Sample)

if (anyNA(clinical_df$Sample) ||
    any(trimws(clinical_df$Sample) == "")) {
  stop("Missing or blank clinical sample IDs detected.")
}

if (anyDuplicated(clinical_df$Sample)) {
  stop("Duplicated sample IDs detected in clinical data.")
}

clinical_missing <- setdiff(sample_ids, clinical_df$Sample)
clinical_extra <- setdiff(clinical_df$Sample, sample_ids)

if (length(clinical_missing) > 0 ||
    length(clinical_extra) > 0) {
  stop(
    "Clinical and expression sample universes do not match. ",
    "Missing clinical records: ", length(clinical_missing),
    "; extra clinical records: ", length(clinical_extra), "."
  )
}

clinical_aligned <- clinical_df[
  match(sample_ids, clinical_df$Sample),
  ,
  drop = FALSE
]

if (!identical(clinical_aligned$Sample, sample_ids)) {
  stop("Clinical data could not be aligned to expression sample order.")
}

clinical_alignment_audit <- data.frame(
  Metric = c(
    "Expression sample number",
    "Clinical sample number",
    "Clinical records missing for expression samples",
    "Extra clinical records",
    "Clinical sample order identical after alignment"
  ),
  Value = c(
    nrow(datExpr_pre),
    nrow(clinical_aligned),
    length(clinical_missing),
    length(clinical_extra),
    as.character(identical(clinical_aligned$Sample, sample_ids))
  ),
  stringsAsFactors = FALSE
)

write.table(
  clinical_alignment_audit,
  file = file.path(
    output_dir,
    "Workflow-04-2D_clinical_alignment_audit.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  clinical_aligned,
  file = file.path(
    output_dir,
    "Workflow-04-2D_clinical_aligned_to_expression.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 7. Sample manifest
############################

sample_manifest <- data.frame(
  Sample = sample_ids,
  Sample_order = seq_along(sample_ids),
  TSS = substr(sample_ids, 6, 7),
  Sample_type_code = substr(sample_ids, 14, 15),
  stringsAsFactors = FALSE
)

write.table(
  sample_manifest,
  file = file.path(
    output_dir,
    "Workflow-04-2D_primary_tumor_sample_manifest.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 8. WGCNA goodSamplesGenes QC
############################

cat("\n========== goodSamplesGenes QC ==========\n")

gsg <- WGCNA::goodSamplesGenes(
  datExpr_pre,
  verbose = 3
)

good_samples_df <- data.frame(
  Sample = rownames(datExpr_pre),
  GoodSample = gsg$goodSamples,
  stringsAsFactors = FALSE
)

good_genes_df <- data.frame(
  Ensembl_ID = colnames(datExpr_pre),
  GoodGene = gsg$goodGenes,
  stringsAsFactors = FALSE
)

write.table(
  good_samples_df,
  file = file.path(
    output_dir,
    "Workflow-04-2D_goodSamplesGenes_sample_status.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  good_genes_df,
  file = file.path(
    output_dir,
    "Workflow-04-2D_goodSamplesGenes_gene_status.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

gsg_summary <- data.frame(
  Metric = c(
    "allOK",
    "Total samples",
    "Good samples",
    "Bad samples",
    "Total genes",
    "Good genes",
    "Bad genes"
  ),
  Value = c(
    as.character(gsg$allOK),
    nrow(datExpr_pre),
    sum(gsg$goodSamples),
    sum(!gsg$goodSamples),
    ncol(datExpr_pre),
    sum(gsg$goodGenes),
    sum(!gsg$goodGenes)
  ),
  stringsAsFactors = FALSE
)

write.table(
  gsg_summary,
  file = file.path(
    output_dir,
    "Workflow-04-2D_goodSamplesGenes_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (!gsg$allOK) {
  
  writeLines(
    c(
      "Bad samples:",
      rownames(datExpr_pre)[!gsg$goodSamples],
      "",
      "Bad genes:",
      colnames(datExpr_pre)[!gsg$goodGenes]
    ),
    con = file.path(
      output_dir,
      "Workflow-04-2D_goodSamplesGenes_failed_items.txt"
    )
  )
  
  stop(
    "goodSamplesGenes() identified bad samples or genes. ",
    "No samples or genes were automatically removed. ",
    "Review the QC output before proceeding."
  )
}

# No samples or genes are removed.
datExpr <- datExpr_pre

if (nrow(datExpr) != 513) {
  stop("The verified WGCNA matrix does not contain 513 samples.")
}


############################
# 9. Select QC-only highly variable genes
############################

# This subset is used only for PCA, clustering and connectivity.
# It does not replace or modify the complete MAD-positive set.

gene_mad_qc <- apply(
  datExpr,
  2,
  stats::mad,
  constant = 1.4826,
  na.rm = FALSE
)

if (anyNA(gene_mad_qc) || any(!is.finite(gene_mad_qc))) {
  stop("Invalid MAD values encountered during QC gene selection.")
}

if (any(gene_mad_qc <= 0)) {
  stop("A non-positive MAD value was detected in the MAD-positive matrix.")
}

top_n <- min(5000L, length(gene_mad_qc))

top_genes <- names(
  sort(
    gene_mad_qc,
    decreasing = TRUE
  )
)[seq_len(top_n)]

datExpr_qc_top <- datExpr[
  ,
  top_genes,
  drop = FALSE
]

write.table(
  data.frame(
    Ensembl_ID = top_genes,
    MAD_in_expression_matrix = unname(gene_mad_qc[top_genes]),
    QC_gene_rank = seq_along(top_genes),
    stringsAsFactors = FALSE
  ),
  file = file.path(
    output_dir,
    "Workflow-04-2D_PCA_connectivity_QC_gene_list.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 10. PCA
############################

cat("\n========== PCA ==========\n")

pca <- stats::prcomp(
  datExpr_qc_top,
  center = TRUE,
  scale. = FALSE
)

pca_variance <- pca$sdev^2
pca_variance_percent <- 100 * pca_variance / sum(pca_variance)

pca_variance_df <- data.frame(
  PC = paste0("PC", seq_along(pca_variance_percent)),
  Standard_deviation = pca$sdev,
  Variance = pca_variance,
  Variance_percent = pca_variance_percent,
  Cumulative_variance_percent = cumsum(pca_variance_percent),
  stringsAsFactors = FALSE
)

write.table(
  pca_variance_df,
  file = file.path(
    output_dir,
    "Workflow-04-2D_PCA_variance_explained.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

n_pc_export <- min(5L, ncol(pca$x))

pca_scores <- data.frame(
  Sample = rownames(pca$x),
  pca$x[, seq_len(n_pc_export), drop = FALSE],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pca_scores <- merge(
  sample_manifest,
  pca_scores,
  by = "Sample",
  all.x = TRUE,
  sort = FALSE
)

pca_scores <- pca_scores[
  match(sample_ids, pca_scores$Sample),
  ,
  drop = FALSE
]

if (!identical(pca_scores$Sample, sample_ids)) {
  stop("PCA sample scores are not in the original sample order.")
}

write.table(
  pca_scores,
  file = file.path(
    output_dir,
    "Workflow-04-2D_PCA_sample_scores.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 11. Add clinical annotations to PCA
############################

pca_plot_df <- pca_scores

optional_clinical_columns <- c(
  "pathologic_stage",
  "tumor_normal"
)

for (column_name in optional_clinical_columns) {
  
  if (column_name %in% colnames(clinical_aligned)) {
    pca_plot_df[[column_name]] <- clinical_aligned[
      match(pca_plot_df$Sample, clinical_aligned$Sample),
      column_name
    ]
  }
}

if (ncol(pca_scores) < 4) {
  stop("Fewer than two principal components are available for PCA plotting.")
}

p_pca <- ggplot2::ggplot(
  pca_plot_df,
  ggplot2::aes(x = PC1, y = PC2)
) +
  ggplot2::geom_point(
    size = 2,
    alpha = 0.75
  ) +
  ggplot2::theme_classic() +
  ggplot2::labs(
    title = "PCA of Primary LUAD Samples",
    x = paste0(
      "PC1 (",
      round(pca_variance_percent[1], 2),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(pca_variance_percent[2], 2),
      "%)"
    )
  )

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "Workflow-04-2D_PCA_PC1_PC2.pdf"
  ),
  plot = p_pca,
  width = 7,
  height = 6
)

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "Workflow-04-2D_PCA_PC1_PC2.png"
  ),
  plot = p_pca,
  width = 7,
  height = 6,
  dpi = 300
)

if ("pathologic_stage" %in% colnames(pca_plot_df)) {
  
  p_pca_stage <- ggplot2::ggplot(
    pca_plot_df,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = pathologic_stage
    )
  ) +
    ggplot2::geom_point(
      size = 2,
      alpha = 0.75
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = "PCA of Primary LUAD Samples by Pathologic Stage",
      x = paste0(
        "PC1 (",
        round(pca_variance_percent[1], 2),
        "%)"
      ),
      y = paste0(
        "PC2 (",
        round(pca_variance_percent[2], 2),
        "%)"
      ),
      color = "Pathologic stage"
    )
  
  ggplot2::ggsave(
    filename = file.path(
      figure_dir,
      "Workflow-04-2D_PCA_PC1_PC2_by_stage.pdf"
    ),
    plot = p_pca_stage,
    width = 8,
    height = 6
  )
  
  ggplot2::ggsave(
    filename = file.path(
      figure_dir,
      "Workflow-04-2D_PCA_PC1_PC2_by_stage.png"
    ),
    plot = p_pca_stage,
    width = 8,
    height = 6,
    dpi = 300
  )
}


############################
# 12. Hierarchical clustering
############################

cat("\n========== Hierarchical clustering ==========\n")

sample_dist <- stats::dist(
  datExpr_qc_top,
  method = "euclidean"
)

sample_tree <- stats::hclust(
  sample_dist,
  method = "average"
)

cluster_order_df <- data.frame(
  Cluster_order = seq_along(sample_tree$order),
  Sample = sample_tree$labels[sample_tree$order],
  stringsAsFactors = FALSE
)

write.table(
  cluster_order_df,
  file = file.path(
    output_dir,
    "Workflow-04-2D_hierarchical_clustering_order.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

grDevices::pdf(
  file = file.path(
    figure_dir,
    "Workflow-04-2D_sample_hierarchical_clustering.pdf"
  ),
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

grDevices::dev.off()


############################
# 13. Sample connectivity
############################

cat("\n========== Sample connectivity ==========\n")

# datExpr_qc_top: samples x genes.
# Correlation of its transpose gives sample-by-sample correlations.

sample_cor <- stats::cor(
  t(datExpr_qc_top),
  method = "pearson",
  use = "everything"
)

if (!identical(rownames(sample_cor), sample_ids) ||
    !identical(colnames(sample_cor), sample_ids)) {
  stop("Sample correlation matrix is not aligned to the original sample order.")
}

if (anyNA(sample_cor) || any(!is.finite(sample_cor))) {
  stop("Invalid values detected in the sample correlation matrix.")
}

# Sum of Pearson correlations with all other samples.
sample_connectivity <- rowSums(sample_cor) - 1

connectivity_df <- data.frame(
  Sample = sample_ids,
  Connectivity = unname(sample_connectivity),
  stringsAsFactors = FALSE
)

connectivity_df$Connectivity_rank <- rank(
  connectivity_df$Connectivity,
  ties.method = "average"
)

connectivity_df$Connectivity_percentile <- 100 *
  connectivity_df$Connectivity_rank /
  length(sample_ids)

connectivity_df$TSS <- substr(
  connectivity_df$Sample,
  6,
  7
)

clinical_match_index <- match(
  connectivity_df$Sample,
  clinical_aligned$Sample
)

if ("pathologic_stage" %in% colnames(clinical_aligned)) {
  connectivity_df$Stage <- clinical_aligned$pathologic_stage[
    clinical_match_index
  ]
}

if ("OS.status" %in% colnames(clinical_aligned)) {
  connectivity_df$OS_status <- clinical_aligned$OS.status[
    clinical_match_index
  ]
}

if (!identical(connectivity_df$Sample, sample_ids)) {
  stop("Connectivity sample order differs from the original expression order.")
}

write.table(
  connectivity_df,
  file = file.path(
    output_dir,
    "Workflow-04-2D_sample_connectivity.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 14. Connectivity summary
############################

connectivity_summary <- data.frame(
  Metric = c(
    "N_samples",
    "Minimum",
    "First_quartile",
    "Median",
    "Mean",
    "Third_quartile",
    "Maximum",
    "10th_percentile"
  ),
  Value = c(
    length(sample_connectivity),
    min(sample_connectivity),
    unname(stats::quantile(sample_connectivity, 0.25)),
    stats::median(sample_connectivity),
    mean(sample_connectivity),
    unname(stats::quantile(sample_connectivity, 0.75)),
    max(sample_connectivity),
    unname(stats::quantile(sample_connectivity, 0.10))
  ),
  stringsAsFactors = FALSE
)

write.table(
  connectivity_summary,
  file = file.path(
    output_dir,
    "Workflow-04-2D_sample_connectivity_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 15. TSS-level connectivity audit
############################

tss_groups <- split(
  connectivity_df$Connectivity,
  connectivity_df$TSS
)

tss_summary_df <- do.call(
  rbind,
  lapply(
    names(tss_groups),
    function(tss_id) {
      
      values <- tss_groups[[tss_id]]
      
      data.frame(
        TSS = tss_id,
        N = length(values),
        Mean = mean(values),
        SD = stats::sd(values),
        Median = stats::median(values),
        Q1 = unname(stats::quantile(values, 0.25)),
        Q3 = unname(stats::quantile(values, 0.75)),
        Min = min(values),
        Max = max(values),
        stringsAsFactors = FALSE
      )
    }
  )
)

tss_summary_df <- tss_summary_df[
  order(tss_summary_df$Median),
  ,
  drop = FALSE
]

rownames(tss_summary_df) <- NULL

write.table(
  tss_summary_df,
  file = file.path(
    output_dir,
    "Workflow-04-2D_TSS_connectivity_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 16. Specific TCGA-44 audit
############################

tcga44 <- connectivity_df[
  connectivity_df$TSS == "44",
  ,
  drop = FALSE
]

tcga44 <- tcga44[
  order(tcga44$Connectivity),
  ,
  drop = FALSE
]

write.table(
  tcga44,
  file = file.path(
    output_dir,
    "Workflow-04-2D_TCGA44_connectivity.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 17. Lowest 10% connectivity audit
############################

connectivity_cutoff_10 <- unname(
  stats::quantile(
    connectivity_df$Connectivity,
    0.10
  )
)

low10 <- connectivity_df[
  connectivity_df$Connectivity <= connectivity_cutoff_10,
  ,
  drop = FALSE
]

low10 <- low10[
  order(low10$Connectivity),
  ,
  drop = FALSE
]

low10_summary <- data.frame(
  Metric = c(
    "10th_percentile_cutoff",
    "Number_of_samples_at_or_below_cutoff",
    "Number_of_TSS44_samples",
    "Percentage_of_low_connectivity_samples_from_TSS44"
  ),
  Value = c(
    connectivity_cutoff_10,
    nrow(low10),
    sum(low10$TSS == "44"),
    if (nrow(low10) > 0) {
      100 * mean(low10$TSS == "44")
    } else {
      NA_real_
    }
  ),
  stringsAsFactors = FALSE
)

write.table(
  low10,
  file = file.path(
    output_dir,
    "Workflow-04-2D_low10_connectivity_samples.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  low10_summary,
  file = file.path(
    output_dir,
    "Workflow-04-2D_low10_connectivity_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################
# 18. Connectivity plot
############################

p_connectivity <- ggplot2::ggplot(
  connectivity_df,
  ggplot2::aes(
    x = reorder(Sample, Connectivity),
    y = Connectivity
  )
) +
  ggplot2::geom_point(size = 1.5) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    title = "Sample Connectivity",
    x = "Primary LUAD samples",
    y = "Sum of Pearson correlations"
  )

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "Workflow-04-2D_sample_connectivity.pdf"
  ),
  plot = p_connectivity,
  width = 10,
  height = 5
)

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "Workflow-04-2D_sample_connectivity.png"
  ),
  plot = p_connectivity,
  width = 10,
  height = 5,
  dpi = 300
)


############################
# 19. Export verified WGCNA input matrix
############################

# Preserve all 513 samples and the complete MAD-positive
# gene set. No PCA, clustering, connectivity or TSS-based
# sample/gene filtering is applied.

wgcna_input_df <- data.frame(
  Sample = sample_ids,
  datExpr,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.table(
  wgcna_input_df,
  file = file.path(
    output_dir,
    "Workflow-04-2D_WGCNA_input_MAD_positive_samples_by_genes.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


############################
# 20. Final integrity verification
############################

wgcna_export <- read.delim(
  file.path(
    output_dir,
    "Workflow-04-2D_WGCNA_input_MAD_positive_samples_by_genes.txt"
  ),
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

export_sample_ids <- as.character(wgcna_export$Sample)
export_gene_ids <- colnames(wgcna_export)[-1]

final_checks <- data.frame(
  Check = c(
    "Exactly 513 samples retained",
    "All samples are TCGA primary tumor type 01",
    "Sample IDs unique",
    "Original sample order preserved",
    "Gene IDs unique",
    "Exported gene universe matches input",
    "Exported gene universe matches MAD-positive QC set",
    "No expression values changed during export",
    "No missing expression values",
    "No infinite expression values",
    "No samples removed",
    "No genes removed from MAD-positive input"
  ),
  PASS = c(
    nrow(wgcna_export) == 513,
    all(substr(export_sample_ids, 14, 15) == "01"),
    !anyDuplicated(export_sample_ids),
    identical(export_sample_ids, sample_ids),
    !anyDuplicated(export_gene_ids),
    identical(export_gene_ids, gene_ids),
    setequal(export_gene_ids, mad_positive_ids),
    isTRUE(all.equal(
      unname(as.matrix(wgcna_export[, -1, drop = FALSE])),
      unname(as.matrix(datExpr)),
      tolerance = 0
    )),
    !anyNA(wgcna_export[, -1, drop = FALSE]),
    all(is.finite(as.matrix(wgcna_export[, -1, drop = FALSE]))),
    nrow(wgcna_export) == nrow(datExpr_pre),
    ncol(wgcna_export) - 1 == ncol(datExpr_pre)
  ),
  stringsAsFactors = FALSE
)

write.table(
  final_checks,
  file = file.path(
    output_dir,
    "Workflow-04-2D_final_integrity_audit.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (!all(final_checks$PASS)) {
  stop("At least one final integrity check failed. Review final_integrity_audit.txt.")
}


############################
# 21. Workflow summary
############################

summary_lines <- c(
  "Workflow-04-2D summary",
  "======================",
  "",
  paste("R version:", as.character(getRversion())),
  paste("WGCNA version:", as.character(packageVersion("WGCNA"))),
  paste("ggplot2 version:", as.character(packageVersion("ggplot2"))),
  "",
  paste("Primary tumor samples:", nrow(datExpr)),
  paste("MAD-positive genes:", ncol(datExpr)),
  paste("Genes used for PCA/connectivity QC:", ncol(datExpr_qc_top)),
  "",
  paste("goodSamplesGenes allOK:", gsg$allOK),
  paste("Good samples:", sum(gsg$goodSamples)),
  paste("Good genes:", sum(gsg$goodGenes)),
  "",
  paste("TSS groups:", length(unique(connectivity_df$TSS))),
  paste("TCGA-44 sample number:", nrow(tcga44)),
  paste("Lowest 10% connectivity sample number:", nrow(low10)),
  "",
  "No samples were removed.",
  "No genes were removed from the MAD-positive candidate matrix.",
  "PCA, clustering, connectivity and TSS analyses are descriptive QC only.",
  "WGCNA network construction was not performed.",
  "",
  paste(
    "Final integrity checks passed:",
    all(final_checks$PASS)
  )
)

writeLines(
  summary_lines,
  con = file.path(
    output_dir,
    "Workflow-04-2D_summary.txt"
  )
)

cat("\n============================================\n")
cat("Workflow-04-2D completed successfully.\n")
cat("Samples:", nrow(datExpr), "\n")
cat("MAD-positive genes:", ncol(datExpr), "\n")
cat("QC genes:", ncol(datExpr_qc_top), "\n")
cat("goodSamplesGenes allOK:", gsg$allOK, "\n")
cat("Final checks passed:", all(final_checks$PASS), "\n")
cat("Output directory:", output_dir, "\n")
cat("Figure directory:", figure_dir, "\n")
cat("============================================\n")


