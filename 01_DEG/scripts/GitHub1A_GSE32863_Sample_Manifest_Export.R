
# ============================================================
# GitHub1A_GSE32863_Sample_Manifest_Export.R
#
# Purpose:
#   Export and audit sample manifests for the GSE32863 cohort.
#
# Outputs:
#   1. GSE32863_complete_sample_manifest.csv
#   2. GSE32863_57_complete_pairs.csv
#   3. GSE32863_singleton_samples_excluded_from_paired_analysis.csv
#
# Cohort characteristics established in the GSE32863 audit:
#   - 116 samples
#   - 59 unique subjects
#   - 58 tumor samples
#   - 58 adjacent-normal samples
#   - 57 complete tumor-normal pairs
#   - 1 tumor-only sample
#   - 1 normal-only sample
#
# Important:
#   This script only exports sample manifests.
#   It does NOT perform expression analysis, probe selection,
#   statistical testing, or figure generation.
#
# ============================================================


# ============================================================
# 01. Project configuration and package checks
# ============================================================

# Identify the project root using the here package.
project_root <- here::here()

# Input file
input_dir <- file.path(project_root, "data", "input_data")

series_matrix_file <- file.path(
  input_dir,
  "GSE32863_series_matrix.txt.gz"
)

# Output directories
manifest_dir <- file.path(
  project_root,
  "data",
  "sample_manifests"
)

audit_dir <- file.path(
  project_root,
  "results",
  "audit_reports"
)

# Create output directories if they do not exist.
dir.create(
  manifest_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Required packages
required_packages <- c(
  "here",
  "GEOquery",
  "Biobase"
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
    paste0(
      "The following required packages are missing: ",
      paste(missing_packages, collapse = ", "),
      "\nPlease install them before running this script."
    )
  )
}

# Check the input file.
if (!file.exists(series_matrix_file)) {
  stop(
    paste0(
      "GSE32863 Series Matrix file was not found:\n",
      series_matrix_file,
      "\nPlease check the input file name and location."
    )
  )
}


# ============================================================
# 02. Load the GSE32863 Series Matrix
# ============================================================

message("Reading GSE32863 Series Matrix...")

gse_object <- GEOquery::getGEO(
  filename = series_matrix_file
)

# Handle either a single ExpressionSet or a one-element list.
if (methods::is(gse_object, "ExpressionSet")) {
  
  eset <- gse_object
  
} else if (
  is.list(gse_object) &&
  length(gse_object) == 1L &&
  methods::is(gse_object[[1L]], "ExpressionSet")
) {
  
  eset <- gse_object[[1L]]
  
} else {
  
  stop(
    paste0(
      "Unexpected GEOquery output. ",
      "A single ExpressionSet was expected."
    )
  )
}

# Extract sample metadata.
pdata <- Biobase::pData(eset)

if (nrow(pdata) == 0L) {
  stop("The sample metadata table is empty.")
}

message(
  "Number of samples in metadata: ",
  nrow(pdata)
)


# ============================================================
# 03. Check required metadata fields
# ============================================================

required_metadata_fields <- c(
  "geo_accession",
  "subject:ch1",
  "tissue:ch1",
  "type:ch1"
)

missing_fields <- setdiff(
  required_metadata_fields,
  colnames(pdata)
)

if (length(missing_fields) > 0L) {
  stop(
    paste0(
      "Required metadata fields are missing: ",
      paste(missing_fields, collapse = ", ")
    )
  )
}

# Extract the required fields.
sample_metadata <- data.frame(
  GEO_ID = trimws(as.character(pdata[["geo_accession"]])),
  Subject = trimws(as.character(pdata[["subject:ch1"]])),
  Original_Tissue = trimws(as.character(pdata[["tissue:ch1"]])),
  Type = trimws(as.character(pdata[["type:ch1"]])),
  Title = if ("title" %in% colnames(pdata)) {
    as.character(pdata[["title"]])
  } else {
    NA_character_
  },
  Source = if ("source_name_ch1" %in% colnames(pdata)) {
    as.character(pdata[["source_name_ch1"]])
  } else {
    NA_character_
  },
  stringsAsFactors = FALSE
)

# Check missing or empty identifiers.
identifier_fields <- c(
  "GEO_ID",
  "Subject",
  "Original_Tissue"
)

for (field in identifier_fields) {
  
  invalid_values <- is.na(sample_metadata[[field]]) |
    sample_metadata[[field]] == ""
  
  if (any(invalid_values)) {
    stop(
      paste0(
        "Missing or empty values detected in metadata field: ",
        field
      )
    )
  }
}


# ============================================================
# 04. Validate GEO sample identifiers
# ============================================================

if (anyDuplicated(sample_metadata$GEO_ID) > 0L) {
  stop(
    "Duplicated GEO sample accessions were detected."
  )
}

# Record the original metadata order.
sample_metadata$Metadata_Order <- seq_len(
  nrow(sample_metadata)
)

# Confirm the expected cohort size.
if (nrow(sample_metadata) != 116L) {
  stop(
    paste0(
      "Expected 116 samples, but found ",
      nrow(sample_metadata),
      ". Please check the input Series Matrix."
    )
  )
}


# ============================================================
# 05. Define tumor and adjacent-normal groups
# ============================================================

# Use the audited tissue:ch1 values.
sample_metadata$Group <- NA_character_

sample_metadata$Group[
  sample_metadata$Original_Tissue == "Lung tumor"
] <- "Tumor"

sample_metadata$Group[
  sample_metadata$Original_Tissue == "Normal lung"
] <- "Normal"

# Stop if unexpected tissue labels are present.
if (anyNA(sample_metadata$Group)) {
  
  unexpected_tissues <- unique(
    sample_metadata$Original_Tissue[
      is.na(sample_metadata$Group)
    ]
  )
  
  stop(
    paste0(
      "Unexpected tissue labels were detected: ",
      paste(unexpected_tissues, collapse = "; ")
    )
  )
}

# Set a consistent group order.
sample_metadata$Group <- factor(
  sample_metadata$Group,
  levels = c("Normal", "Tumor")
)

# Check group sizes.
group_counts <- table(sample_metadata$Group)

print(group_counts)

if (
  unname(group_counts["Tumor"]) != 58L ||
  unname(group_counts["Normal"]) != 58L
) {
  stop(
    "The expected cohort composition (58 tumor and 58 normal samples) was not observed."
  )
}


# ============================================================
# 06. Generate the complete sample manifest
# ============================================================

complete_sample_manifest <- sample_metadata[
  order(sample_metadata$Metadata_Order),
  c(
    "GEO_ID",
    "Subject",
    "Group",
    "Original_Tissue",
    "Type",
    "Title",
    "Source"
  )
]

# Convert Group back to character for CSV export.
complete_sample_manifest$Group <- as.character(
  complete_sample_manifest$Group
)

# Export the complete 116-sample manifest.
complete_manifest_file <- file.path(
  manifest_dir,
  "GSE32863_complete_sample_manifest.csv"
)

utils::write.csv(
  complete_sample_manifest,
  file = complete_manifest_file,
  row.names = FALSE,
  na = ""
)

message(
  "Complete sample manifest exported: ",
  complete_manifest_file
)


# ============================================================
# 07. Identify complete tumor-normal pairs
# ============================================================

# Split sample records by subject.
subject_list <- split(
  sample_metadata,
  sample_metadata$Subject
)

pair_records <- lapply(
  names(subject_list),
  function(subject_id) {
    
    subject_data <- subject_list[[subject_id]]
    
    tumor_data <- subject_data[
      subject_data$Group == "Tumor",
      ,
      drop = FALSE
    ]
    
    normal_data <- subject_data[
      subject_data$Group == "Normal",
      ,
      drop = FALSE
    ]
    
    # A complete pair must contain exactly one tumor
    # and exactly one adjacent-normal sample.
    if (
      nrow(tumor_data) == 1L &&
      nrow(normal_data) == 1L
    ) {
      
      data.frame(
        Subject = subject_id,
        GEO_ID_Tumor = tumor_data$GEO_ID,
        GEO_ID_Normal = normal_data$GEO_ID,
        Title_Tumor = tumor_data$Title,
        Title_Normal = normal_data$Title,
        stringsAsFactors = FALSE
      )
      
    } else {
      
      NULL
      
    }
  }
)

# Combine all complete pairs.
complete_pairs <- do.call(
  rbind,
  pair_records
)

if (is.null(complete_pairs)) {
  stop("No complete tumor-normal pairs were identified.")
}

# Sort by subject ID for reproducible output.
complete_pairs <- complete_pairs[
  order(complete_pairs$Subject),
  ,
  drop = FALSE
]

rownames(complete_pairs) <- NULL


# ============================================================
# 08. Validate the complete-pair manifest
# ============================================================

if (nrow(complete_pairs) != 57L) {
  stop(
    paste0(
      "Expected 57 complete pairs, but found ",
      nrow(complete_pairs),
      "."
    )
  )
}

if (anyDuplicated(complete_pairs$Subject) > 0L) {
  stop(
    "Duplicated subject IDs were detected in the complete-pair manifest."
  )
}

if (
  anyDuplicated(complete_pairs$GEO_ID_Tumor) > 0L ||
  anyDuplicated(complete_pairs$GEO_ID_Normal) > 0L
) {
  stop(
    "Duplicated GEO sample IDs were detected in the pair manifest."
  )
}

# Confirm that each paired sample exists in the complete manifest.
if (
  !all(complete_pairs$GEO_ID_Tumor %in%
       complete_sample_manifest$GEO_ID) ||
  !all(complete_pairs$GEO_ID_Normal %in%
       complete_sample_manifest$GEO_ID)
) {
  stop(
    "Some paired GEO sample IDs are absent from the complete manifest."
  )
}

# Confirm that each subject has one tumor and one normal sample.
if (
  any(
    complete_pairs$GEO_ID_Tumor ==
    complete_pairs$GEO_ID_Normal
  )
) {
  stop(
    "A tumor and normal sample have the same GEO accession."
  )
}

# Export the 57 complete pairs.
complete_pairs_file <- file.path(
  manifest_dir,
  "GSE32863_57_complete_pairs.csv"
)

utils::write.csv(
  complete_pairs,
  file = complete_pairs_file,
  row.names = FALSE,
  na = ""
)

message(
  "Complete-pair manifest exported: ",
  complete_pairs_file
)


# ============================================================
# 09. Identify singleton samples excluded from paired analysis
# ============================================================

singleton_records <- lapply(
  names(subject_list),
  function(subject_id) {
    
    subject_data <- subject_list[[subject_id]]
    
    tumor_data <- subject_data[
      subject_data$Group == "Tumor",
      ,
      drop = FALSE
    ]
    
    normal_data <- subject_data[
      subject_data$Group == "Normal",
      ,
      drop = FALSE
    ]
    
    # Tumor-only subject.
    if (
      nrow(tumor_data) == 1L &&
      nrow(normal_data) == 0L
    ) {
      
      return(
        data.frame(
          Subject = subject_id,
          GEO_ID = tumor_data$GEO_ID,
          Group = "Tumor",
          Original_Tissue = tumor_data$Original_Tissue,
          Exclusion_Reason =
            "No matched adjacent-normal sample",
          stringsAsFactors = FALSE
        )
      )
    }
    
    # Normal-only subject.
    if (
      nrow(tumor_data) == 0L &&
      nrow(normal_data) == 1L
    ) {
      
      return(
        data.frame(
          Subject = subject_id,
          GEO_ID = normal_data$GEO_ID,
          Group = "Normal",
          Original_Tissue = normal_data$Original_Tissue,
          Exclusion_Reason =
            "No matched tumor sample",
          stringsAsFactors = FALSE
        )
      )
    }
    
    # Any unexpected subject structure is not silently ignored.
    if (
      nrow(tumor_data) != 1L ||
      nrow(normal_data) != 1L
    ) {
      
      stop(
        paste0(
          "Unexpected sample configuration for subject: ",
          subject_id
        )
      )
    }
    
    NULL
  }
)

singleton_samples <- do.call(
  rbind,
  singleton_records
)

if (is.null(singleton_samples)) {
  stop("No singleton samples were identified.")
}

singleton_samples <- singleton_samples[
  order(singleton_samples$Subject),
  ,
  drop = FALSE
]

rownames(singleton_samples) <- NULL


# ============================================================
# 10. Validate singleton samples
# ============================================================

if (nrow(singleton_samples) != 2L) {
  stop(
    paste0(
      "Expected 2 singleton samples, but found ",
      nrow(singleton_samples),
      "."
    )
  )
}

if (
  sum(singleton_samples$Group == "Tumor") != 1L ||
  sum(singleton_samples$Group == "Normal") != 1L
) {
  stop(
    "The expected one tumor-only and one normal-only sample were not identified."
  )
}

# Export singleton exclusion records.
singleton_file <- file.path(
  manifest_dir,
  "GSE32863_singleton_samples_excluded_from_paired_analysis.csv"
)

utils::write.csv(
  singleton_samples,
  file = singleton_file,
  row.names = FALSE,
  na = ""
)

message(
  "Singleton exclusion manifest exported: ",
  singleton_file
)


# ============================================================
# 11. Final cohort consistency checks
# ============================================================

n_total_samples <- nrow(complete_sample_manifest)

n_tumor_samples <- sum(
  complete_sample_manifest$Group == "Tumor"
)

n_normal_samples <- sum(
  complete_sample_manifest$Group == "Normal"
)

n_unique_subjects <- length(
  unique(complete_sample_manifest$Subject)
)

n_complete_pairs <- nrow(complete_pairs)

n_singleton_subjects <- nrow(singleton_samples)

n_paired_samples <- 2L * n_complete_pairs

# Check that all 116 samples are accounted for.
if (
  n_paired_samples + n_singleton_subjects !=
  n_total_samples
) {
  stop(
    "The number of paired and singleton samples does not equal the total cohort size."
  )
}

# Check the expected subject count.
if (n_unique_subjects != 59L) {
  stop(
    paste0(
      "Expected 59 unique subjects, but found ",
      n_unique_subjects,
      "."
    )
  )
}

# Check the expected cohort structure.
if (
  n_total_samples != 116L ||
  n_tumor_samples != 58L ||
  n_normal_samples != 58L ||
  n_complete_pairs != 57L ||
  n_singleton_subjects != 2L ||
  n_paired_samples != 114L
) {
  stop(
    "Final cohort consistency checks failed."
  )
}


# ============================================================
# 12. Save a cohort audit summary
# ============================================================

cohort_audit <- data.frame(
  Metric = c(
    "Total samples",
    "Tumor samples",
    "Adjacent-normal samples",
    "Unique subjects",
    "Complete tumor-normal pairs",
    "Samples included in paired analysis",
    "Singleton subjects",
    "Tumor-only samples",
    "Normal-only samples"
  ),
  Value = c(
    n_total_samples,
    n_tumor_samples,
    n_normal_samples,
    n_unique_subjects,
    n_complete_pairs,
    n_paired_samples,
    n_singleton_subjects,
    sum(singleton_samples$Group == "Tumor"),
    sum(singleton_samples$Group == "Normal")
  ),
  stringsAsFactors = FALSE
)

audit_summary_file <- file.path(
  audit_dir,
  "GSE32863_sample_manifest_audit_summary.csv"
)

utils::write.csv(
  cohort_audit,
  file = audit_summary_file,
  row.names = FALSE
)


# ============================================================
# 13. Print final audit results
# ============================================================

cat("\n")
cat("============================================================\n")
cat("GSE32863 SAMPLE MANIFEST EXPORT: COMPLETED\n")
cat("============================================================\n\n")

print(cohort_audit, row.names = FALSE)

cat("\nGenerated files:\n")

cat(
  "1. ",
  complete_manifest_file,
  "\n",
  sep = ""
)

cat(
  "2. ",
  complete_pairs_file,
  "\n",
  sep = ""
)

cat(
  "3. ",
  singleton_file,
  "\n",
  sep = ""
)

cat(
  "4. ",
  audit_summary_file,
  "\n",
  sep = ""
)

cat("\nAll cohort consistency checks passed.\n")

cat("============================================================\n")



# ============================================================
# 14. Cross-file sample ID consistency audit
# ============================================================

# Collect sample IDs from the three exported manifests.

all_manifest_ids <- complete_sample_manifest$GEO_ID

paired_tumor_ids <- complete_pairs$GEO_ID_Tumor
paired_normal_ids <- complete_pairs$GEO_ID_Normal

singleton_ids <- singleton_samples$GEO_ID

paired_ids <- c(
  paired_tumor_ids,
  paired_normal_ids
)

accounted_ids <- c(
  paired_ids,
  singleton_ids
)

# 1. Check uniqueness of all sample IDs.
check_unique_manifest_ids <-
  !anyDuplicated(all_manifest_ids)

# 2. Check that all paired samples are present.
check_paired_ids_present <-
  all(paired_ids %in% all_manifest_ids)

# 3. Check the group assignments of paired samples.
tumor_group_check <- all(
  complete_sample_manifest$Group[
    match(paired_tumor_ids, all_manifest_ids)
  ] == "Tumor"
)

normal_group_check <- all(
  complete_sample_manifest$Group[
    match(paired_normal_ids, all_manifest_ids)
  ] == "Normal"
)

# 4. Check that singleton samples are present.
check_singleton_ids_present <-
  all(singleton_ids %in% all_manifest_ids)

# 5. Check that no sample is duplicated across
#    paired and singleton records.
check_no_duplicate_accounting <-
  !anyDuplicated(accounted_ids)

# 6. Check exact coverage of the complete manifest.
check_exact_sample_coverage <-
  setequal(accounted_ids, all_manifest_ids) &&
  length(accounted_ids) == length(all_manifest_ids)

# Compile audit results.
cross_file_audit <- data.frame(
  Check = c(
    "Unique GEO IDs in complete manifest",
    "All paired IDs present in complete manifest",
    "Tumor IDs assigned to Tumor group",
    "Normal IDs assigned to Normal group",
    "All singleton IDs present in complete manifest",
    "No duplicated IDs across paired and singleton records",
    "Exact coverage of all 116 samples"
  ),
  PASS = c(
    check_unique_manifest_ids,
    check_paired_ids_present,
    tumor_group_check,
    normal_group_check,
    check_singleton_ids_present,
    check_no_duplicate_accounting,
    check_exact_sample_coverage
  ),
  stringsAsFactors = FALSE
)

print(cross_file_audit, row.names = FALSE)

# Stop if any check fails.
if (!all(cross_file_audit$PASS)) {
  stop(
    "Cross-file sample ID consistency audit failed."
  )
}

# Export the audit report.
cross_file_audit_file <- file.path(
  audit_dir,
  "GSE32863_cross_file_sample_ID_audit.csv"
)

utils::write.csv(
  cross_file_audit,
  cross_file_audit_file,
  row.names = FALSE
)

message(
  "Cross-file sample ID audit passed. Report saved to: ",
  cross_file_audit_file
)






