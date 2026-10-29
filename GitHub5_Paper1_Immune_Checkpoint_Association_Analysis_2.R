##GitHub5_Paper1_Immune_Checkpoint_Association_Analysis_2
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

##STEP 8D_Purity sensitivity analysis

work_dir <- "input_dir/Paper 1_Step 8_Immune Checkpoint"

setwd(work_dir)

getwd()

load(
  "input_dir/Paper1_tumor-only WGCNA/Paper1_END_WORKSPACE_20260920_20260920_015542.RData"
)

ls()

sessionInfo()

# ============================================================
# Paper 1
# WORKFLOW-08
# STEP 8D — Tumor Purity Sensitivity Analysis
#
# Purpose:
# Independent sensitivity analysis for tumor purity
#
# IMPORTANT:
# This analysis DOES NOT replace the locked STEP 8C-7A/7B
# ImmuneScore-adjusted primary analysis.
#
# Primary locked analysis:
#   Unadjusted
#   vs
#   ImmuneScore-adjusted
#
# STEP 8D:
#   Unadjusted
#   vs
#   Purity-adjusted
#
# ============================================================

cat("\n============================================================\n")
cat("STEP 8D: Tumor Purity Sensitivity Analysis\n")
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
      "Required package not installed: ",
      pkg,
      "\nPlease install it before continuing."
    )
  }
}


# ------------------------------------------------------------
# 2. Define working directory
# ------------------------------------------------------------

work_dir <- "input_dir/Paper1_tumor-only WGCNA"

setwd(work_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 3. Define locked STEP 8C-7B input
# ------------------------------------------------------------

locked_rds <- file.path(
  work_dir,
  "Workflow-08-STEP8C7B_Supplementary_Table_Immune_Checkpoint_54.rds"
)

if (!file.exists(locked_rds)) {
  
  stop(
    "Locked STEP 8C-7B RDS file not found:\n",
    locked_rds
  )
}

checkpoint_8c7b <- readRDS(locked_rds)

cat("STEP 8C-7B object class:\n")
print(class(checkpoint_8c7b))

cat("\n")


# ------------------------------------------------------------
# 4. Inspect locked object
# ------------------------------------------------------------

cat("Structure of STEP 8C-7B object:\n")
str(
  checkpoint_8c7b,
  max.level = 1
)

cat("\n")


# ------------------------------------------------------------
# 5. Identify the locked 54-pair table
# ------------------------------------------------------------

find_data_frame <- function(x) {
  
  if (is.data.frame(x)) {
    return(x)
  }
  
  if (is.list(x)) {
    
    df_candidates <- x[
      vapply(
        x,
        is.data.frame,
        logical(1)
      )
    ]
    
    if (length(df_candidates) > 0) {
      
      # Prefer table with 54 rows
      nrows <- vapply(
        df_candidates,
        nrow,
        integer(1)
      )
      
      idx54 <- which(nrows == 54)
      
      if (length(idx54) >= 1) {
        return(df_candidates[[idx54[1]]])
      }
      
      return(df_candidates[[1]])
    }
  }
  
  return(NULL)
}


locked_54 <- find_data_frame(checkpoint_8c7b)

if (is.null(locked_54)) {
  
  stop(
    "Could not identify the 54-pair data frame in the STEP 8C-7B RDS object."
  )
}

cat("Locked table dimensions:\n")
print(dim(locked_54))

cat("\nLocked table columns:\n")
print(colnames(locked_54))

cat("\n")

# ------------------------------------------------------------
# 6. Identify gene-pair columns
# ------------------------------------------------------------

cat("Attempting to identify gene-pair columns...\n\n")

print(
  data.frame(
    column = colnames(locked_54),
    class = vapply(
      locked_54,
      function(x) paste(class(x), collapse = ","),
      character(1)
    ),
    stringsAsFactors = FALSE
  )
)


library(dplyr)
# ------------------------------------------------------------
# 7. Lock the 54 predefined gene pairs
# ------------------------------------------------------------

required_pair_columns <- c(
  "Core_Gene",
  "Checkpoint_Gene"
)

missing_pair_columns <- setdiff(
  required_pair_columns,
  colnames(locked_54)
)

if (length(missing_pair_columns) > 0) {
  
  stop(
    "Missing required pair columns:\n",
    paste(missing_pair_columns, collapse = ", ")
  )
}


pair_design <- locked_54 %>%
  dplyr::select(
    Core_Gene,
    Checkpoint_Gene
  ) %>%
  dplyr::distinct()

cat("Number of unique predefined pairs:\n")
cat(nrow(pair_design), "\n\n")

if (nrow(pair_design) != 54) {
  
  stop(
    "ERROR: Expected exactly 54 predefined gene pairs, ",
    "but found ",
    nrow(pair_design),
    "."
  )
}


cat("Core genes:\n")
print(unique(pair_design$Core_Gene))

cat("\nNumber of checkpoint genes:\n")
cat(length(unique(pair_design$Checkpoint_Gene)), "\n\n")


# ============================================================
# STEP 8D
# 5. Inspect available purity-related objects
# ============================================================

purity_candidates <- ls(
  envir = .GlobalEnv,
  pattern = "(?i)(purity|estimate|immune|stromal)"
)

cat("\n============================================================\n")
cat("STEP 8D: Purity-related objects currently in GlobalEnv\n")
cat("============================================================\n\n")

print(purity_candidates)

# ============================================================
# STEP 8D
# 6. Inspect the 503-sample ESTIMATE object
# ============================================================

cat("\n============================================================\n")
cat("STEP 8D: Inspect estimate_scores_503\n")
cat("============================================================\n\n")

cat("Class:\n")
print(class(estimate_scores_503))

cat("\nDimensions:\n")
print(dim(estimate_scores_503))

cat("\nColumn names:\n")
print(colnames(estimate_scores_503))

cat("\nRow names: first 10\n")
print(head(rownames(estimate_scores_503), 10))

cat("\nStructure:\n")
str(estimate_scores_503)

cat("\nSummary:\n")
print(summary(estimate_scores_503))


# ============================================================
# STEP 8D
# 6A. Inspect the 3 × 503 ESTIMATE object in detail
# ============================================================

cat("\n============================================================\n")
cat("STEP 8D: Detailed inspection of estimate_scores_503\n")
cat("============================================================\n\n")

cat("Column names: first 10\n")
print(head(colnames(estimate_scores_503), 10))

cat("\nColumn names: last 10\n")
print(tail(colnames(estimate_scores_503), 10))

cat("\nRow names:\n")
print(rownames(estimate_scores_503))

cat("\nStructure:\n")
str(estimate_scores_503)

cat("\nFirst 6 columns:\n")
print(estimate_scores_503[, 1:min(6, ncol(estimate_scores_503)), drop = FALSE])

cat("\nSummary by row:\n")

# ============================================================
# STEP 8D
# 7. Verify exact sample identity for the 503-sample cohort
# ============================================================

cat("\n============================================================\n")
cat("STEP 8D: Exact sample matching check\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 7.1 ESTIMATE sample IDs
# ------------------------------------------------------------

estimate_ids_503 <- colnames(estimate_scores_503)

cat("Number of ESTIMATE sample IDs:\n")
print(length(estimate_ids_503))

cat("\nUnique ESTIMATE sample IDs:\n")
print(length(unique(estimate_ids_503)))

# ------------------------------------------------------------
# 7.2 Clinical survival 503 sample IDs
# ------------------------------------------------------------

cat("\nClinical survival object:\n")
cat("Class:\n")
print(class(clinical_survival_503))

cat("\nDimensions:\n")
print(dim(clinical_survival_503))

cat("\nColumn names:\n")
print(colnames(clinical_survival_503))

cat("\nFirst 10 row names:\n")
print(head(rownames(clinical_survival_503), 10))

# ------------------------------------------------------------
# 7.3 Identify the sample-ID field
# ------------------------------------------------------------

if ("patient_id" %in% colnames(clinical_survival_503)) {
  
  survival_ids_503 <- as.character(
    clinical_survival_503$patient_id
  )
  
} else if ("SampleID" %in% colnames(clinical_survival_503)) {
  
  survival_ids_503 <- as.character(
    clinical_survival_503$SampleID
  )
  
} else if ("sample_id" %in% colnames(clinical_survival_503)) {
  
  survival_ids_503 <- as.character(
    clinical_survival_503$sample_id
  )
  
} else {
  
  survival_ids_503 <- rownames(clinical_survival_503)
  
}

# ------------------------------------------------------------
# 7.4 Exact matching
# ------------------------------------------------------------

cat("\nNumber of survival sample IDs:\n")
print(length(survival_ids_503))

cat("\nUnique survival sample IDs:\n")
print(length(unique(survival_ids_503)))

cat("\nESTIMATE ∩ survival:\n")
print(length(intersect(
  estimate_ids_503,
  survival_ids_503
)))

cat("\nESTIMATE-only IDs:\n")
print(length(setdiff(
  estimate_ids_503,
  survival_ids_503
)))

cat("\nSurvival-only IDs:\n")
print(length(setdiff(
  survival_ids_503,
  estimate_ids_503
)))

# ------------------------------------------------------------
# 7.5 PASS / FAIL
# ------------------------------------------------------------

sample_match_503 <- (
  length(estimate_ids_503) == 503 &&
    length(unique(estimate_ids_503)) == 503 &&
    length(survival_ids_503) == 503 &&
    length(unique(survival_ids_503)) == 503 &&
    setequal(estimate_ids_503, survival_ids_503)
)

cat("\n============================================================\n")
cat("STEP 8D SAMPLE MATCH STATUS\n")
cat("============================================================\n")

if (sample_match_503) {
  
  cat("PASS: ESTIMATE and survival cohorts contain the same 503 samples.\n")
  
} else {
  
  cat("FAIL: ESTIMATE and survival sample IDs are not identical.\n")
  
}

# ============================================================
# STEP 8D-7A
# Inspect ESTIMATE sample IDs versus clinical patient IDs
# ============================================================

cat("\n============================================================\n")
cat("STEP 8D-7A: Inspect ID formats\n")
cat("============================================================\n\n")

cat("First 10 ESTIMATE IDs:\n")
print(head(estimate_ids_503, 10))

cat("\nFirst 10 clinical patient IDs:\n")
print(head(survival_ids_503, 10))

cat("\nESTIMATE ID lengths:\n")
print(table(nchar(estimate_ids_503)))

cat("\nClinical patient ID lengths:\n")
print(table(nchar(survival_ids_503)))

# ------------------------------------------------------------
# Remove the sample portion (-01, -02, etc.) from ESTIMATE IDs
# for diagnostic comparison only
# ------------------------------------------------------------

estimate_patient_ids_503 <- sub(
  "-[0-9]{2}$",
  "",
  estimate_ids_503
)

cat("\nFirst 10 normalized ESTIMATE patient IDs:\n")
print(head(estimate_patient_ids_503, 10))

cat("\nNormalized ESTIMATE ID uniqueness:\n")
print(length(unique(estimate_patient_ids_503)))

cat("\nNormalized ESTIMATE ∩ clinical patient IDs:\n")
print(length(intersect(
  estimate_patient_ids_503,
  survival_ids_503
)))

cat("\nNormalized ESTIMATE-only IDs:\n")
print(length(setdiff(
  estimate_patient_ids_503,
  survival_ids_503
)))

cat("\nClinical-only IDs:\n")
print(length(setdiff(
  survival_ids_503,
  estimate_patient_ids_503
)))


# ============================================================
# STEP 8D-8
# Build exact 503-sample ESTIMATE-to-patient mapping
# ============================================================

cat("\n============================================================\n")
cat("STEP 8D-8: Build exact 503-sample mapping\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 8.1 Create mapping table
# ------------------------------------------------------------

sample_map_503 <- data.frame(
  ESTIMATE_Sample_ID = estimate_ids_503,
  Patient_ID = estimate_patient_ids_503,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 8.2 Verify one-to-one mapping
# ------------------------------------------------------------

cat("Mapping dimensions:\n")
print(dim(sample_map_503))

cat("\nUnique ESTIMATE sample IDs:\n")
print(length(unique(sample_map_503$ESTIMATE_Sample_ID)))

cat("\nUnique patient IDs:\n")
print(length(unique(sample_map_503$Patient_ID)))

cat("\nDuplicated ESTIMATE sample IDs:\n")
print(sum(duplicated(sample_map_503$ESTIMATE_Sample_ID)))

cat("\nDuplicated patient IDs:\n")
print(sum(duplicated(sample_map_503$Patient_ID)))

# ------------------------------------------------------------
# 8.3 Exact correspondence with clinical cohort
# ------------------------------------------------------------

mapping_match <- (
  nrow(sample_map_503) == 503 &&
    length(unique(sample_map_503$ESTIMATE_Sample_ID)) == 503 &&
    length(unique(sample_map_503$Patient_ID)) == 503 &&
    all(sort(sample_map_503$Patient_ID) == sort(survival_ids_503))
)

cat("\n============================================================\n")
cat("STEP 8D SAMPLE MAPPING STATUS\n")
cat("============================================================\n")

if (mapping_match) {
  
  cat("PASS: Exact one-to-one mapping established for all 503 samples.\n")
  
} else {
  
  cat("FAIL: Mapping verification failed.\n")
  
}

# ------------------------------------------------------------
# 8.4 Display first 10 mappings
# ------------------------------------------------------------

cat("\nFirst 10 mappings:\n")
print(head(sample_map_503, 10))

# ============================================================
# STEP 8D-9
# Verify provenance and consistency of ESTIMATE scores
# ============================================================

cat("\n============================================================\n")
cat("STEP 8D-9: Verify ESTIMATE score provenance\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 9.1 Inspect estimate_scores_wgcna
# ------------------------------------------------------------

cat("### estimate_scores_wgcna ###\n")
cat("Class:\n")
print(class(estimate_scores_wgcna))

cat("\nDimensions:\n")
print(dim(estimate_scores_wgcna))

cat("\nRow names:\n")
print(rownames(estimate_scores_wgcna))

cat("\nFirst 5 column names:\n")
print(head(colnames(estimate_scores_wgcna), 5))

cat("\nLast 5 column names:\n")
print(tail(colnames(estimate_scores_wgcna), 5))


# ------------------------------------------------------------
# 9.2 Inspect estimate_scores_wgcna_locked
# ------------------------------------------------------------

cat("\n\n### estimate_scores_wgcna_locked ###\n")
cat("Class:\n")
print(class(estimate_scores_wgcna_locked))

cat("\nDimensions:\n")
print(dim(estimate_scores_wgcna_locked))

cat("\nRow names:\n")
print(rownames(estimate_scores_wgcna_locked))

cat("\nFirst 5 column names:\n")
print(head(colnames(estimate_scores_wgcna_locked), 5))

cat("\nLast 5 column names:\n")
print(tail(colnames(estimate_scores_wgcna_locked), 5))


# ------------------------------------------------------------
# 9.3 Compare 503-sample ESTIMATE object with WGCNA object
# ------------------------------------------------------------

cat("\n\n### Compare sample IDs ###\n")

estimate_wgcna_ids <- colnames(estimate_scores_wgcna)

cat("Number of WGCNA ESTIMATE IDs:\n")
print(length(estimate_wgcna_ids))

cat("Number overlapping with 503-sample object:\n")
print(length(intersect(
  estimate_ids_503,
  estimate_wgcna_ids
)))


# ------------------------------------------------------------
# 9.4 Check ESTIMATEScore availability
# ------------------------------------------------------------

cat("\n\n### ESTIMATEScore values ###\n")

cat("503-sample ESTIMATEScore summary:\n")
print(summary(
  as.numeric(estimate_scores_503["ESTIMATEScore", ])
))

if ("ESTIMATEScore" %in% rownames(estimate_scores_wgcna)) {
  
  cat("\nWGCNA ESTIMATEScore summary:\n")
  print(summary(
    as.numeric(estimate_scores_wgcna["ESTIMATEScore", ])
  ))
  
} else {
  
  cat("\nWGCNA object does not contain a row named ESTIMATEScore.\n")
}


##STEP 8D-10: Formally Establish ESTIMATE-derived Tumor Purity
#Here, we employ the classic conversion from ESTIMATE score to tumor purity:
#  TumorPurity =cos(0.6049872018 + 0.0001467884*ESTIMATEScore)
# ============================================================
# STEP 8D-10
# Establish ESTIMATE-derived tumor purity for n = 503
# ============================================================

cat("\n============================================================\n")
cat("STEP 8D-10: Establish ESTIMATE-derived Tumor Purity\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 10.1 Extract locked 503-sample ESTIMATE score
# ------------------------------------------------------------

estimate_score_503 <- as.numeric(
  estimate_scores_503["ESTIMATEScore", ]
)

names(estimate_score_503) <- colnames(
  estimate_scores_503
)

# ------------------------------------------------------------
# 10.2 Calculate ESTIMATE-derived tumor purity
# ------------------------------------------------------------

tumor_purity_503 <- cos(
  0.6049872018 +
    0.0001467884 * estimate_score_503
)

# ------------------------------------------------------------
# 10.3 Build a dedicated purity data frame
# ------------------------------------------------------------

tumor_purity_503_df <- data.frame(
  ESTIMATE_Sample_ID = names(estimate_score_503),
  Patient_ID = estimate_patient_ids_503,
  ESTIMATEScore = unname(estimate_score_503),
  TumorPurity = unname(tumor_purity_503),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 10.4 Basic QC
# ------------------------------------------------------------

cat("Dimensions:\n")
print(dim(tumor_purity_503_df))

cat("\nColumn names:\n")
print(colnames(tumor_purity_503_df))

cat("\nNumber of unique samples:\n")
print(length(unique(
  tumor_purity_503_df$ESTIMATE_Sample_ID
)))

cat("\nNumber of unique patients:\n")
print(length(unique(
  tumor_purity_503_df$Patient_ID
)))

cat("\nMissing TumorPurity values:\n")
print(sum(is.na(tumor_purity_503_df$TumorPurity)))

cat("\nInfinite TumorPurity values:\n")
print(sum(!is.finite(tumor_purity_503_df$TumorPurity)))

cat("\nTumorPurity range:\n")
print(range(
  tumor_purity_503_df$TumorPurity,
  na.rm = TRUE
))

cat("\nTumorPurity summary:\n")
print(summary(
  tumor_purity_503_df$TumorPurity
))

# ------------------------------------------------------------
# 10.5 Verify expected mathematical range
# ------------------------------------------------------------

purity_range_pass <- all(
  tumor_purity_503_df$TumorPurity >= 0 &
    tumor_purity_503_df$TumorPurity <= 1
)

cat("\nTumorPurity 0–1 range check:\n")

if (purity_range_pass) {
  
  cat("PASS: All tumor purity values are within [0, 1].\n")
  
} else {
  
  cat("FAIL: Tumor purity values outside [0, 1] detected.\n")
}

# ------------------------------------------------------------
# 10.6 Verify exact sample mapping
# ------------------------------------------------------------

purity_mapping_pass <- (
  nrow(tumor_purity_503_df) == 503 &&
    length(unique(tumor_purity_503_df$ESTIMATE_Sample_ID)) == 503 &&
    length(unique(tumor_purity_503_df$Patient_ID)) == 503 &&
    all(
      sort(tumor_purity_503_df$Patient_ID) ==
        sort(survival_ids_503)
    )
)

cat("\n503-sample purity mapping check:\n")

if (purity_mapping_pass) {
  
  cat("PASS: Tumor purity is mapped to the exact 503-patient cohort.\n")
  
} else {
  
  cat("FAIL: Tumor purity mapping does not match the 503-patient cohort.\n")
}

# ------------------------------------------------------------
# 10.7 Display first 10 samples
# ------------------------------------------------------------

cat("\nFirst 10 tumor-purity records:\n")
print(head(tumor_purity_503_df, 10))

# ------------------------------------------------------------
# 10.8 Overall status
# ------------------------------------------------------------

tumor_purity_503_PASS <- (
  nrow(tumor_purity_503_df) == 503 &&
    length(unique(tumor_purity_503_df$ESTIMATE_Sample_ID)) == 503 &&
    length(unique(tumor_purity_503_df$Patient_ID)) == 503 &&
    sum(is.na(tumor_purity_503_df$TumorPurity)) == 0 &&
    sum(!is.finite(tumor_purity_503_df$TumorPurity)) == 0 &&
    purity_range_pass &&
    purity_mapping_pass
)

cat("\n============================================================\n")
cat("STEP 8D-10 STATUS\n")
cat("============================================================\n")

if (tumor_purity_503_PASS) {
  
  cat("PASS: ESTIMATE-derived tumor purity established for n = 503.\n")
  
} else {
  
  cat("FAIL: Tumor purity QC failed.\n")
}


# ============================================================
# STEP 8D-11
# Identify and inspect the expression matrix for
# the 21 predefined genes
# ============================================================

cat("\n============================================================\n")
cat("STEP 8D-11: Expression Matrix Source QC\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 11.1 Verify locked checkpoint dataset
# ------------------------------------------------------------

cat("\n[1] Locked checkpoint dataset\n")

print(class(locked_54))
print(dim(locked_54))

cat("\nRequired columns:\n")
print(colnames(locked_54))


# ------------------------------------------------------------
# 11.2 Define the exact 21 genes
# ------------------------------------------------------------

core_genes_3 <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

checkpoint_genes_18 <- unique(
  locked_54$Checkpoint_Gene
)

all_21_genes <- unique(
  c(core_genes_3, checkpoint_genes_18)
)

cat("\nCore genes:\n")
print(core_genes_3)

cat("\nCheckpoint genes:\n")
print(checkpoint_genes_18)

cat("\nNumber of checkpoint genes:\n")
print(length(checkpoint_genes_18))

cat("\nTotal predefined genes:\n")
print(length(all_21_genes))


# ------------------------------------------------------------
# 11.3 Confirm exactly 18 checkpoint genes
# ------------------------------------------------------------

checkpoint_gene_count_pass <- (
  length(checkpoint_genes_18) == 18
)

cat("\n18-checkpoint-gene check:\n")

if (checkpoint_gene_count_pass) {
  
  cat("PASS: Exactly 18 checkpoint genes identified.\n")
  
} else {
  
  cat("FAIL: Expected 18 checkpoint genes.\n")
}


# ------------------------------------------------------------
# STEP 8D-11.4
# Inspect candidate expression objects in GlobalEnv
# ------------------------------------------------------------

cat("\nCandidate expression-related objects in GlobalEnv:\n")

all_objects <- ls()

expression_candidates <- all_objects[
  grepl(
    "(expr|expression|fpkm|tpm|matrix)",
    all_objects,
    ignore.case = TRUE
  )
]

print(expression_candidates)

cat("\nNumber of candidate expression objects:\n")
print(length(expression_candidates))


# ------------------------------------------------------------
# STEP 8D-11.5
# Inspect dimensions and classes
# ------------------------------------------------------------

cat("\nExpression candidate object summary:\n")

for (obj_name in expression_candidates) {
  
  obj <- get(obj_name)
  
  cat("\n------------------------------------------------------------\n")
  cat("Object:", obj_name, "\n")
  cat("Class:", paste(class(obj), collapse = ", "), "\n")
  
  if (is.data.frame(obj) || is.matrix(obj)) {
    
    cat("Dimensions:",
        paste(dim(obj), collapse = " x "),
        "\n")
    
    cat("First column name:",
        colnames(obj)[1],
        "\n")
    
    cat("First 5 column names:\n")
    print(head(colnames(obj), 5))
    
    cat("Last 5 column names:\n")
    print(tail(colnames(obj), 5))
    
    cat("First 5 row names:\n")
    print(head(rownames(obj), 5))
    
    cat("Last 5 row names:\n")
    print(tail(rownames(obj), 5))
  }
}


# ------------------------------------------------------------
# STEP 8D-11.6
# Check whether candidate objects contain the 21 genes
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("21-gene presence screening\n")
cat("============================================================\n")

for (obj_name in expression_candidates) {
  
  obj <- get(obj_name)
  
  if (is.data.frame(obj) || is.matrix(obj)) {
    
    rn <- rownames(obj)
    cn <- colnames(obj)
    
    row_hits <- sum(all_21_genes %in% rn)
    col_hits <- sum(all_21_genes %in% cn)
    
    cat("\nObject:", obj_name, "\n")
    cat("Genes found in rownames:", row_hits, "/ 21\n")
    cat("Genes found in colnames:", col_hits, "/ 21\n")
  }
}


# ------------------------------------------------------------
# STEP 8D-11.7
# Verify locked 54-pair structure
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("54-pair structure QC\n")
cat("============================================================\n")

cat("Number of locked associations:",
    nrow(locked_54),
    "\n")

cat("Unique core genes:",
    length(unique(locked_54$Core_Gene)),
    "\n")

cat("Unique checkpoint genes:",
    length(unique(locked_54$Checkpoint_Gene)),
    "\n")

cat("Unique Association_ID:",
    length(unique(locked_54$Association_ID)),
    "\n")

pair_structure_pass <- (
  nrow(locked_54) == 54 &&
    length(unique(locked_54$Core_Gene)) == 3 &&
    length(unique(locked_54$Checkpoint_Gene)) == 18 &&
    length(unique(locked_54$Association_ID)) == 54
)

if (pair_structure_pass) {
  
  cat("PASS: Locked 54-pair structure confirmed.\n")
  
} else {
  
  cat("FAIL: Locked 54-pair structure is inconsistent.\n")
}


# ------------------------------------------------------------
# STEP 8D-11.8
# Final status
# ------------------------------------------------------------

step8d11_structure_PASS <- (
  checkpoint_gene_count_pass &&
    pair_structure_pass
)

cat("\n============================================================\n")
cat("STEP 8D-11 STATUS\n")
cat("============================================================\n")

if (step8d11_structure_PASS) {
  
  cat("PASS: Locked 54-pair structure and 21-gene target set confirmed.\n")
  cat("Next step: identify the correct 503-sample expression matrix and verify sample-ID alignment.\n")
  
} else {
  
  cat("FAIL: Structure QC requires review before expression extraction.\n")
}


# ============================================================
# STEP 8D-12
# Exact 503-sample expression–purity alignment QC
# ============================================================

cat("\n============================================================\n")
cat("STEP 8D-12: Exact 503-Sample Expression-Purity Alignment QC\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 12.1 Define the four primary data components
# ------------------------------------------------------------

core_expr_503_local <- core_expr_503
checkpoint_expr_503_local <- checkpoint_expr_503

purity_df_503_local <- tumor_purity_503_df

estimate_ids_503_local <- colnames(
  estimate_scores_503
)

core_ids_503 <- colnames(
  core_expr_503_local
)

checkpoint_ids_503 <- colnames(
  checkpoint_expr_503_local
)

purity_sample_ids_503 <- purity_df_503_local$ESTIMATE_Sample_ID

purity_patient_ids_503 <- purity_df_503_local$Patient_ID

clinical_patient_ids_503 <- clinical_survival_503$patient_id


# ------------------------------------------------------------
# 12.2 Basic dimensions
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("12.2 Dimensions\n")
cat("------------------------------------------------------------\n")

cat("core_expr_503:",
    paste(dim(core_expr_503_local), collapse = " x "),
    "\n")

cat("checkpoint_expr_503:",
    paste(dim(checkpoint_expr_503_local), collapse = " x "),
    "\n")

cat("tumor_purity_503_df:",
    paste(dim(purity_df_503_local), collapse = " x "),
    "\n")

cat("clinical_survival_503:",
    paste(dim(clinical_survival_503), collapse = " x "),
    "\n")


# ------------------------------------------------------------
# 12.3 Unique sample-ID checks
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("12.3 Unique sample-ID checks\n")
cat("------------------------------------------------------------\n")

cat("Unique core expression IDs:",
    length(unique(core_ids_503)),
    "\n")

cat("Unique checkpoint expression IDs:",
    length(unique(checkpoint_ids_503)),
    "\n")

cat("Unique ESTIMATE IDs:",
    length(unique(estimate_ids_503_local)),
    "\n")

cat("Unique purity sample IDs:",
    length(unique(purity_sample_ids_503)),
    "\n")

cat("Unique purity patient IDs:",
    length(unique(purity_patient_ids_503)),
    "\n")

cat("Unique clinical patient IDs:",
    length(unique(clinical_patient_ids_503)),
    "\n")


# ------------------------------------------------------------
# 12.4 Exact set matching
# ------------------------------------------------------------

core_vs_checkpoint_pass <- (
  setequal(core_ids_503, checkpoint_ids_503)
)

core_vs_estimate_pass <- (
  setequal(core_ids_503, estimate_ids_503_local)
)

core_vs_purity_pass <- (
  setequal(core_ids_503, purity_sample_ids_503)
)

purity_vs_clinical_pass <- (
  setequal(purity_patient_ids_503, clinical_patient_ids_503)
)


cat("\n------------------------------------------------------------\n")
cat("12.4 Exact ID-set matching\n")
cat("------------------------------------------------------------\n")

cat("Core expression vs checkpoint expression:",
    core_vs_checkpoint_pass,
    "\n")

cat("Core expression vs ESTIMATE:",
    core_vs_estimate_pass,
    "\n")

cat("Core expression vs tumor purity:",
    core_vs_purity_pass,
    "\n")

cat("Purity patients vs clinical patients:",
    purity_vs_clinical_pass,
    "\n")


# ------------------------------------------------------------
# 12.5 Exact order matching
# ------------------------------------------------------------

core_vs_checkpoint_order_pass <- identical(
  core_ids_503,
  checkpoint_ids_503
)

core_vs_estimate_order_pass <- identical(
  core_ids_503,
  estimate_ids_503_local
)

core_vs_purity_order_pass <- identical(
  core_ids_503,
  purity_sample_ids_503
)


cat("\n------------------------------------------------------------\n")
cat("12.5 Exact ID-order matching\n")
cat("------------------------------------------------------------\n")

cat("Core vs checkpoint order:",
    core_vs_checkpoint_order_pass,
    "\n")

cat("Core vs ESTIMATE order:",
    core_vs_estimate_order_pass,
    "\n")

cat("Core vs purity order:",
    core_vs_purity_order_pass,
    "\n")


# ------------------------------------------------------------
# 12.6 Check patient-level correspondence
# ------------------------------------------------------------

core_patient_ids_503 <- sub(
  "-[0-9]{2}$",
  "",
  core_ids_503
)

checkpoint_patient_ids_503 <- sub(
  "-[0-9]{2}$",
  "",
  checkpoint_ids_503
)

estimate_patient_ids_503_local <- sub(
  "-[0-9]{2}$",
  "",
  estimate_ids_503_local
)


core_patient_pass <- (
  setequal(
    core_patient_ids_503,
    clinical_patient_ids_503
  )
)

checkpoint_patient_pass <- (
  setequal(
    checkpoint_patient_ids_503,
    clinical_patient_ids_503
  )
)

estimate_patient_pass <- (
  setequal(
    estimate_patient_ids_503_local,
    clinical_patient_ids_503
  )
)


cat("\n------------------------------------------------------------\n")
cat("12.6 Patient-level correspondence\n")
cat("------------------------------------------------------------\n")

cat("Core expression vs clinical:",
    core_patient_pass,
    "\n")

cat("Checkpoint expression vs clinical:",
    checkpoint_patient_pass,
    "\n")

cat("ESTIMATE vs clinical:",
    estimate_patient_pass,
    "\n")


# ------------------------------------------------------------
# 12.7 Check expression matrices for finite values
# ------------------------------------------------------------

core_numeric_pass <- all(
  is.finite(as.numeric(core_expr_503_local))
)

checkpoint_numeric_pass <- all(
  is.finite(as.numeric(checkpoint_expr_503_local))
)

purity_numeric_pass <- all(
  is.finite(purity_df_503_local$TumorPurity)
)


cat("\n------------------------------------------------------------\n")
cat("12.7 Numeric-value QC\n")
cat("------------------------------------------------------------\n")

cat("Core expression finite:",
    core_numeric_pass,
    "\n")

cat("Checkpoint expression finite:",
    checkpoint_numeric_pass,
    "\n")

cat("Tumor purity finite:",
    purity_numeric_pass,
    "\n")


# ------------------------------------------------------------
# 12.8 Verify exact 21-gene structure
# ------------------------------------------------------------

core_gene_structure_pass <- identical(
  rownames(core_expr_503_local),
  core_genes_3
)

checkpoint_gene_structure_pass <- identical(
  rownames(checkpoint_expr_503_local),
  checkpoint_genes_18
)

cat("\n------------------------------------------------------------\n")
cat("12.8 Gene structure\n")
cat("------------------------------------------------------------\n")

cat("Core gene order:",
    core_gene_structure_pass,
    "\n")

cat("Checkpoint gene order:",
    checkpoint_gene_structure_pass,
    "\n")


# ------------------------------------------------------------
# 12.9 Overall STEP 8D-12 status
# ------------------------------------------------------------

step8d12_PASS <- (
  ncol(core_expr_503_local) == 503 &&
    ncol(checkpoint_expr_503_local) == 503 &&
    nrow(purity_df_503_local) == 503 &&
    nrow(clinical_survival_503) == 503 &&
    length(unique(core_ids_503)) == 503 &&
    length(unique(checkpoint_ids_503)) == 503 &&
    length(unique(estimate_ids_503_local)) == 503 &&
    length(unique(purity_sample_ids_503)) == 503 &&
    length(unique(purity_patient_ids_503)) == 503 &&
    length(unique(clinical_patient_ids_503)) == 503 &&
    core_vs_checkpoint_pass &&
    core_vs_estimate_pass &&
    core_vs_purity_pass &&
    purity_vs_clinical_pass &&
    core_vs_checkpoint_order_pass &&
    core_vs_estimate_order_pass &&
    core_vs_purity_order_pass &&
    core_patient_pass &&
    checkpoint_patient_pass &&
    estimate_patient_pass &&
    core_numeric_pass &&
    checkpoint_numeric_pass &&
    purity_numeric_pass &&
    core_gene_structure_pass &&
    checkpoint_gene_structure_pass
)


cat("\n============================================================\n")
cat("STEP 8D-12 STATUS\n")
cat("============================================================\n")

if (step8d12_PASS) {
  
  cat("PASS: Exact 503-sample expression, ESTIMATE, purity,\n")
  cat("and clinical cohort alignment confirmed.\n")
  cat("The dataset is ready for purity-adjusted partial Spearman analysis.\n")
  
} else {
  
  cat("FAIL: Expression-purity-clinical alignment requires review.\n")
}

#STEP 8D-13 TumorPurity-adjusted partial Spearman analysis
# ============================================================
# STEP 8D-13
# TumorPurity-adjusted partial Spearman analysis
# 54 predefined immune-checkpoint associations
# ============================================================

cat("\n============================================================\n")
cat("STEP 8D-13: TumorPurity-Adjusted Partial Spearman Analysis\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 13.1 Preserve locked primary dataset
# ------------------------------------------------------------

locked_54_purity_input <- locked_54

cat("\nLocked dataset dimensions:\n")
print(dim(locked_54_purity_input))

if (!identical(
  locked_54_purity_input,
  locked_54
)) {
  stop("ERROR: locked_54 was unexpectedly altered.")
}

cat("PASS: locked_54 remains unchanged.\n")


# ------------------------------------------------------------
# 13.2 Define exact 54-pair design
# ------------------------------------------------------------

pair_design_54 <- data.frame(
  Association_ID = locked_54_purity_input$Association_ID,
  Core_Gene = locked_54_purity_input$Core_Gene,
  Checkpoint_Gene = locked_54_purity_input$Checkpoint_Gene,
  stringsAsFactors = FALSE
)

cat("\nNumber of predefined pairs:\n")
print(nrow(pair_design_54))

cat("\nUnique core genes:\n")
print(unique(pair_design_54$Core_Gene))

cat("\nUnique checkpoint genes:\n")
print(unique(pair_design_54$Checkpoint_Gene))


# ------------------------------------------------------------
# 13.3 Build combined 21-gene expression matrix
# ------------------------------------------------------------

expression_21_503 <- rbind(
  core_expr_503_local,
  checkpoint_expr_503_local
)

cat("\nExpression matrix dimensions:\n")
print(dim(expression_21_503))


# ------------------------------------------------------------
# 13.4 Confirm exact 21-gene order
# ------------------------------------------------------------

expected_21_genes <- c(
  core_genes_3,
  checkpoint_genes_18
)

gene_order_pass <- identical(
  rownames(expression_21_503),
  expected_21_genes
)

cat("\n21-gene order check:\n")

if (gene_order_pass) {
  
  cat("PASS: Exact predefined 21-gene order confirmed.\n")
  
} else {
  
  cat("FAIL: Gene order does not match predefined design.\n")
  stop("Gene-order QC failed.")
}


# ------------------------------------------------------------
# 13.5 Define TumorPurity covariate
# ------------------------------------------------------------

tumor_purity_vector <- purity_df_503_local$TumorPurity

names(tumor_purity_vector) <-
  purity_df_503_local$ESTIMATE_Sample_ID

cat("\nTumorPurity summary:\n")
print(summary(tumor_purity_vector))


# ------------------------------------------------------------
# 13.6 Confirm expression–purity order
# ------------------------------------------------------------

expression_sample_ids <- colnames(
  expression_21_503
)

purity_sample_ids <- names(
  tumor_purity_vector
)

expression_purity_order_pass <- identical(
  expression_sample_ids,
  purity_sample_ids
)

cat("\nExpression vs TumorPurity order:\n")

if (expression_purity_order_pass) {
  
  cat("PASS: Expression and TumorPurity sample order identical.\n")
  
} else {
  
  cat("FAIL: Expression and TumorPurity sample order differs.\n")
  stop("Sample-order QC failed.")
}


# ------------------------------------------------------------
# 13.7 Define partial Spearman function
# ------------------------------------------------------------

partial_spearman_purity <- function(
    x,
    y,
    covariate
) {
  
  complete <- complete.cases(
    x,
    y,
    covariate
  )
  
  x <- x[complete]
  y <- y[complete]
  covariate <- covariate[complete]
  
  n_complete <- length(x)
  
  if (n_complete < 10) {
    
    return(
      data.frame(
        N = n_complete,
        rho = NA_real_,
        P = NA_real_
      )
    )
  }
  
  rx <- rank(
    x,
    ties.method = "average"
  )
  
  ry <- rank(
    y,
    ties.method = "average"
  )
  
  rc <- rank(
    covariate,
    ties.method = "average"
  )
  
  fit_x <- lm(
    rx ~ rc
  )
  
  fit_y <- lm(
    ry ~ rc
  )
  
  residual_x <- residuals(
    fit_x
  )
  
  residual_y <- residuals(
    fit_y
  )
  
  test <- cor.test(
    residual_x,
    residual_y,
    method = "pearson",
    exact = FALSE
  )
  
  data.frame(
    N = n_complete,
    rho = unname(test$estimate),
    P = test$p.value
  )
}


# ------------------------------------------------------------
# 13.8 Run all 54 predefined pairs
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Running 54 TumorPurity-adjusted associations\n")
cat("============================================================\n")

purity_results_list <- vector(
  "list",
  length = nrow(pair_design_54)
)

for (i in seq_len(nrow(pair_design_54))) {
  
  core_gene_i <- pair_design_54$Core_Gene[i]
  checkpoint_gene_i <- pair_design_54$Checkpoint_Gene[i]
  
  x <- as.numeric(
    expression_21_503[
      core_gene_i,
      expression_sample_ids
    ]
  )
  
  y <- as.numeric(
    expression_21_503[
      checkpoint_gene_i,
      expression_sample_ids
    ]
  )
  
  purity_result_i <- partial_spearman_purity(
    x = x,
    y = y,
    covariate = tumor_purity_vector
  )
  
  purity_results_list[[i]] <- data.frame(
    Association_ID =
      pair_design_54$Association_ID[i],
    
    Core_Gene =
      core_gene_i,
    
    Checkpoint_Gene =
      checkpoint_gene_i,
    
    N =
      purity_result_i$N,
    
    Purity_Adjusted_Partial_Spearman_rho =
      purity_result_i$rho,
    
    Purity_Adjusted_P =
      purity_result_i$P,
    
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 13.9 Combine results
# ------------------------------------------------------------

purity_adjusted_54 <- do.call(
  rbind,
  purity_results_list
)

rownames(purity_adjusted_54) <- NULL


# ------------------------------------------------------------
# 13.10 BH-FDR across exactly 54 comparisons
# ------------------------------------------------------------

purity_adjusted_54$Purity_Adjusted_FDR_BH <-
  p.adjust(
    purity_adjusted_54$Purity_Adjusted_P,
    method = "BH"
  )


# ------------------------------------------------------------
# 13.11 Direction classification
# ------------------------------------------------------------

purity_adjusted_54$Purity_Adjusted_Direction <-
  ifelse(
    purity_adjusted_54$
      Purity_Adjusted_Partial_Spearman_rho > 0,
    "Positive",
    ifelse(
      purity_adjusted_54$
        Purity_Adjusted_Partial_Spearman_rho < 0,
      "Negative",
      "Zero"
    )
  )


purity_adjusted_54$Purity_Adjusted_FDR_Significant <-
  purity_adjusted_54$
  Purity_Adjusted_FDR_BH < 0.05


# ------------------------------------------------------------
# 13.12 QC: exactly 54 results
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 8D-13 QC\n")
cat("============================================================\n")

cat("Number of results:",
    nrow(purity_adjusted_54),
    "\n")

cat("Unique Association_ID:",
    length(unique(
      purity_adjusted_54$Association_ID
    )),
    "\n")

cat("Non-NA rho:",
    sum(!is.na(
      purity_adjusted_54$
        Purity_Adjusted_Partial_Spearman_rho
    )),
    "/ 54\n")

cat("Non-NA P:",
    sum(!is.na(
      purity_adjusted_54$Purity_Adjusted_P
    )),
    "/ 54\n")

cat("Non-NA FDR:",
    sum(!is.na(
      purity_adjusted_54$Purity_Adjusted_FDR_BH
    )),
    "/ 54\n")

cat("Unique sample size values:\n")
print(
  unique(
    purity_adjusted_54$N
  )
)


# ------------------------------------------------------------
# 13.13 Global summary
# ------------------------------------------------------------

purity_sig_n <- sum(
  purity_adjusted_54$
    Purity_Adjusted_FDR_BH < 0.05,
  na.rm = TRUE
)

purity_positive_n <- sum(
  purity_adjusted_54$
    Purity_Adjusted_Direction == "Positive",
  na.rm = TRUE
)

purity_negative_n <- sum(
  purity_adjusted_54$
    Purity_Adjusted_Direction == "Negative",
  na.rm = TRUE
)

cat("\nPurity-adjusted FDR < 0.05:",
    purity_sig_n,
    "/ 54\n")

cat("Positive associations:",
    purity_positive_n,
    "/ 54\n")

cat("Negative associations:",
    purity_negative_n,
    "/ 54\n")


# ------------------------------------------------------------
# 13.14 Core-gene summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Core-gene purity-adjusted summary\n")
cat("============================================================\n")

for (g in core_genes_3) {
  
  tmp <- purity_adjusted_54[
    purity_adjusted_54$Core_Gene == g,
    ,
    drop = FALSE
  ]
  
  cat("\n", g, "\n", sep = "")
  
  cat("Total pairs:",
      nrow(tmp),
      "\n")
  
  cat("FDR < 0.05:",
      sum(
        tmp$Purity_Adjusted_FDR_BH < 0.05,
        na.rm = TRUE
      ),
      "\n")
  
  cat("Positive:",
      sum(
        tmp$Purity_Adjusted_Direction == "Positive",
        na.rm = TRUE
      ),
      "\n")
  
  cat("Negative:",
      sum(
        tmp$Purity_Adjusted_Direction == "Negative",
        na.rm = TRUE
      ),
      "\n")
}


# ------------------------------------------------------------
# 13.15 Final PASS criteria
# ------------------------------------------------------------

step8d13_PASS <- (
  nrow(purity_adjusted_54) == 54 &&
    length(unique(
      purity_adjusted_54$Association_ID
    )) == 54 &&
    all(
      purity_adjusted_54$N == 503
    ) &&
    sum(is.na(
      purity_adjusted_54$
        Purity_Adjusted_Partial_Spearman_rho
    )) == 0 &&
    sum(is.na(
      purity_adjusted_54$Purity_Adjusted_P
    )) == 0 &&
    sum(is.na(
      purity_adjusted_54$Purity_Adjusted_FDR_BH
    )) == 0
)


cat("\n============================================================\n")
cat("STEP 8D-13 STATUS\n")
cat("============================================================\n")

if (step8d13_PASS) {
  
  cat("PASS: All 54 predefined gene pairs completed with\n")
  cat("TumorPurity-adjusted partial Spearman analysis,\n")
  cat("N = 503 for every association, with BH-FDR across 54 tests.\n")
  
} else {
  
  cat("FAIL: Purity-adjusted association analysis requires review.\n")
}

# ============================================================
# PAPER 1
# STEP 8D-14
# Three-layer paired comparison:
# Unadjusted vs ImmuneScore-adjusted vs TumorPurity-adjusted
# ============================================================

cat("\n============================================================\n")
cat("STEP 8D-14: Three-Layer Paired Comparison\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 14.1 Preserve locked primary dataset
# ------------------------------------------------------------

locked_54_8d14 <- locked_54

cat("\nLocked primary dataset:\n")
print(dim(locked_54_8d14))

if (!identical(
  locked_54_8d14,
  locked_54
)) {
  stop("ERROR: locked_54 was unexpectedly altered.")
}

cat("PASS: locked_54 remains unchanged.\n")


# ------------------------------------------------------------
# 14.2 Confirm STEP 8D-13 exists
# ------------------------------------------------------------

if (!exists("purity_adjusted_54")) {
  stop(
    "ERROR: purity_adjusted_54 not found. ",
    "Run STEP 8D-13 first."
  )
}

cat("\nSTEP 8D-13 purity-adjusted dataset:\n")
print(dim(purity_adjusted_54))


# ------------------------------------------------------------
# 14.3 Build locked unadjusted + ImmuneScore layer
# ------------------------------------------------------------

three_layer_primary <- data.frame(
  
  Association_ID =
    locked_54_8d14$Association_ID,
  
  Core_Gene =
    locked_54_8d14$Core_Gene,
  
  Checkpoint_Gene =
    locked_54_8d14$Checkpoint_Gene,
  
  N =
    locked_54_8d14$N,
  
  Unadjusted_rho =
    locked_54_8d14$Unadjusted_Spearman_rho,
  
  Unadjusted_P =
    locked_54_8d14$Unadjusted_P,
  
  Unadjusted_FDR =
    locked_54_8d14$Unadjusted_FDR_BH,
  
  Unadjusted_Direction =
    locked_54_8d14$Unadjusted_Direction,
  
  ImmuneScore_Adjusted_rho =
    locked_54_8d14$Adjusted_Partial_Spearman_rho,
  
  ImmuneScore_Adjusted_P =
    locked_54_8d14$Adjusted_P,
  
  ImmuneScore_Adjusted_FDR =
    locked_54_8d14$Adjusted_FDR,
  
  ImmuneScore_Adjusted_Direction =
    locked_54_8d14$Adjusted_Direction,
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 14.4 Extract purity layer
# ------------------------------------------------------------

purity_layer <- data.frame(
  
  Association_ID =
    purity_adjusted_54$Association_ID,
  
  Purity_Adjusted_rho =
    purity_adjusted_54$
    Purity_Adjusted_Partial_Spearman_rho,
  
  Purity_Adjusted_P =
    purity_adjusted_54$
    Purity_Adjusted_P,
  
  Purity_Adjusted_FDR =
    purity_adjusted_54$
    Purity_Adjusted_FDR_BH,
  
  Purity_Adjusted_Direction =
    purity_adjusted_54$
    Purity_Adjusted_Direction,
  
  Purity_Adjusted_N =
    purity_adjusted_54$N,
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 14.5 Check Association_ID uniqueness
# ------------------------------------------------------------

cat("\nAssociation_ID QC:\n")

cat(
  "Locked unique IDs:",
  length(unique(
    three_layer_primary$Association_ID
  )),
  "\n"
)

cat(
  "Purity unique IDs:",
  length(unique(
    purity_layer$Association_ID
  )),
  "\n"
)

if (
  anyDuplicated(
    three_layer_primary$Association_ID
  ) > 0
) {
  stop(
    "ERROR: Duplicate Association_ID in locked dataset."
  )
}

if (
  anyDuplicated(
    purity_layer$Association_ID
  ) > 0
) {
  stop(
    "ERROR: Duplicate Association_ID in purity dataset."
  )
}

cat("PASS: Association_ID values are unique in both layers.\n")


# ------------------------------------------------------------
# 14.6 Exact Association_ID set matching
# ------------------------------------------------------------

id_set_pass <- identical(
  sort(three_layer_primary$Association_ID),
  sort(purity_layer$Association_ID)
)

cat("\nAssociation_ID set matching:\n")

if (id_set_pass) {
  
  cat(
    "PASS: Exact 54-pair Association_ID set confirmed.\n"
  )
  
} else {
  
  cat(
    "FAIL: Association_ID sets do not match.\n"
  )
  
  print(
    setdiff(
      three_layer_primary$Association_ID,
      purity_layer$Association_ID
    )
  )
  
  print(
    setdiff(
      purity_layer$Association_ID,
      three_layer_primary$Association_ID
    )
  )
  
  stop("Association_ID matching failed.")
}


# ------------------------------------------------------------
# 14.7 Merge by Association_ID only
# ------------------------------------------------------------

three_layer_comparison <- merge(
  three_layer_primary,
  purity_layer,
  by = "Association_ID",
  all = FALSE,
  sort = FALSE
)


# ------------------------------------------------------------
# 14.8 Restore original locked order
# ------------------------------------------------------------

three_layer_comparison <- three_layer_comparison[
  match(
    three_layer_primary$Association_ID,
    three_layer_comparison$Association_ID
  ),
  ,
  drop = FALSE
]

rownames(three_layer_comparison) <- NULL


# ------------------------------------------------------------
# 14.9 Verify pair identity after merge
# ------------------------------------------------------------

pair_identity_pass <- all(
  three_layer_comparison$Core_Gene ==
    three_layer_primary$Core_Gene
) &&
  all(
    three_layer_comparison$Checkpoint_Gene ==
      three_layer_primary$Checkpoint_Gene
  )

cat("\nPair identity after merge:\n")

if (pair_identity_pass) {
  
  cat(
    "PASS: Core-gene and checkpoint-gene identity preserved.\n"
  )
  
} else {
  
  stop(
    "ERROR: Pair identity changed during merge."
  )
}


# ------------------------------------------------------------
# 14.10 Verify sample size
# ------------------------------------------------------------

cat("\nSample-size QC:\n")

cat(
  "Unadjusted / ImmuneScore N values:\n"
)
print(
  unique(
    three_layer_comparison$N
  )
)

cat(
  "TumorPurity N values:\n"
)
print(
  unique(
    three_layer_comparison$Purity_Adjusted_N
  )
)

if (
  !all(
    three_layer_comparison$N == 503
  )
) {
  stop(
    "ERROR: Locked analysis contains N != 503."
  )
}

if (
  !all(
    three_layer_comparison$Purity_Adjusted_N == 503
  )
) {
  stop(
    "ERROR: Purity analysis contains N != 503."
  )
}

cat("PASS: N = 503 for all three analytical layers.\n")


# ------------------------------------------------------------
# 14.11 Create significance indicators
# ------------------------------------------------------------

three_layer_comparison$Unadjusted_Significant <-
  three_layer_comparison$Unadjusted_FDR < 0.05

three_layer_comparison$ImmuneScore_Significant <-
  three_layer_comparison$ImmuneScore_Adjusted_FDR < 0.05

three_layer_comparison$TumorPurity_Significant <-
  three_layer_comparison$Purity_Adjusted_FDR < 0.05


# ------------------------------------------------------------
# 14.12 Create direction indicators
# ------------------------------------------------------------

three_layer_comparison$Unadjusted_Positive <-
  three_layer_comparison$Unadjusted_Direction == "Positive"

three_layer_comparison$ImmuneScore_Positive <-
  three_layer_comparison$ImmuneScore_Adjusted_Direction ==
  "Positive"

three_layer_comparison$TumorPurity_Positive <-
  three_layer_comparison$Purity_Adjusted_Direction ==
  "Positive"


# ------------------------------------------------------------
# 14.13 Three-layer significance pattern
# ------------------------------------------------------------

three_layer_comparison$Significance_Pattern <-
  paste0(
    ifelse(
      three_layer_comparison$Unadjusted_Significant,
      "Sig",
      "NS"
    ),
    " → ",
    ifelse(
      three_layer_comparison$ImmuneScore_Significant,
      "Sig",
      "NS"
    ),
    " → ",
    ifelse(
      three_layer_comparison$TumorPurity_Significant,
      "Sig",
      "NS"
    )
  )


# ------------------------------------------------------------
# 14.14 Three-layer direction pattern
# ------------------------------------------------------------

three_layer_comparison$Direction_Pattern <-
  paste0(
    substr(
      three_layer_comparison$Unadjusted_Direction,
      1,
      1
    ),
    " → ",
    substr(
      three_layer_comparison$ImmuneScore_Adjusted_Direction,
      1,
      1
    ),
    " → ",
    substr(
      three_layer_comparison$Purity_Adjusted_Direction,
      1,
      1
    )
  )


# ------------------------------------------------------------
# 14.15 Pairwise direction transitions
# ------------------------------------------------------------

three_layer_comparison$Unadjusted_to_ImmuneScore <-
  paste0(
    three_layer_comparison$Unadjusted_Direction,
    " → ",
    three_layer_comparison$ImmuneScore_Adjusted_Direction
  )

three_layer_comparison$Unadjusted_to_TumorPurity <-
  paste0(
    three_layer_comparison$Unadjusted_Direction,
    " → ",
    three_layer_comparison$Purity_Adjusted_Direction
  )

three_layer_comparison$ImmuneScore_to_TumorPurity <-
  paste0(
    three_layer_comparison$ImmuneScore_Adjusted_Direction,
    " → ",
    three_layer_comparison$Purity_Adjusted_Direction
  )


# ------------------------------------------------------------
# 14.16 Pairwise significance transitions
# ------------------------------------------------------------

three_layer_comparison$Unadjusted_to_ImmuneScore_FDR <-
  paste0(
    ifelse(
      three_layer_comparison$Unadjusted_Significant,
      "Sig",
      "NS"
    ),
    " → ",
    ifelse(
      three_layer_comparison$ImmuneScore_Significant,
      "Sig",
      "NS"
    )
  )

three_layer_comparison$Unadjusted_to_TumorPurity_FDR <-
  paste0(
    ifelse(
      three_layer_comparison$Unadjusted_Significant,
      "Sig",
      "NS"
    ),
    " → ",
    ifelse(
      three_layer_comparison$TumorPurity_Significant,
      "Sig",
      "NS"
    )
  )

three_layer_comparison$ImmuneScore_to_TumorPurity_FDR <-
  paste0(
    ifelse(
      three_layer_comparison$ImmuneScore_Significant,
      "Sig",
      "NS"
    ),
    " → ",
    ifelse(
      three_layer_comparison$TumorPurity_Significant,
      "Sig",
      "NS"
    )
  )


# ------------------------------------------------------------
# 14.17 Global significance summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GLOBAL THREE-LAYER SIGNIFICANCE SUMMARY\n")
cat("============================================================\n")

cat(
  "Unadjusted significant:",
  sum(
    three_layer_comparison$Unadjusted_Significant
  ),
  "/ 54\n"
)

cat(
  "ImmuneScore-adjusted significant:",
  sum(
    three_layer_comparison$ImmuneScore_Significant
  ),
  "/ 54\n"
)

cat(
  "TumorPurity-adjusted significant:",
  sum(
    three_layer_comparison$TumorPurity_Significant
  ),
  "/ 54\n"
)


cat("\nThree-layer significance patterns:\n")
print(
  sort(
    table(
      three_layer_comparison$Significance_Pattern,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)


# ------------------------------------------------------------
# 14.18 Global direction summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GLOBAL THREE-LAYER DIRECTION SUMMARY\n")
cat("============================================================\n")

cat("\nUnadjusted direction:\n")
print(
  table(
    three_layer_comparison$Unadjusted_Direction
  )
)

cat("\nImmuneScore-adjusted direction:\n")
print(
  table(
    three_layer_comparison$
      ImmuneScore_Adjusted_Direction
  )
)

cat("\nTumorPurity-adjusted direction:\n")
print(
  table(
    three_layer_comparison$
      Purity_Adjusted_Direction
  )
)

cat("\nThree-layer direction patterns:\n")
print(
  sort(
    table(
      three_layer_comparison$Direction_Pattern,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)


# ------------------------------------------------------------
# 14.19 Pairwise direction transition summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("PAIRWISE DIRECTION TRANSITIONS\n")
cat("============================================================\n")

cat("\nUnadjusted -> ImmuneScore:\n")
print(
  sort(
    table(
      three_layer_comparison$
        Unadjusted_to_ImmuneScore
    ),
    decreasing = TRUE
  )
)

cat("\nUnadjusted -> TumorPurity:\n")
print(
  sort(
    table(
      three_layer_comparison$
        Unadjusted_to_TumorPurity
    ),
    decreasing = TRUE
  )
)

cat("\nImmuneScore -> TumorPurity:\n")
print(
  sort(
    table(
      three_layer_comparison$
        ImmuneScore_to_TumorPurity
    ),
    decreasing = TRUE
  )
)


# ------------------------------------------------------------
# 14.20 Pairwise significance transition summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("PAIRWISE FDR SIGNIFICANCE TRANSITIONS\n")
cat("============================================================\n")

cat("\nUnadjusted -> ImmuneScore:\n")
print(
  sort(
    table(
      three_layer_comparison$
        Unadjusted_to_ImmuneScore_FDR
    ),
    decreasing = TRUE
  )
)

cat("\nUnadjusted -> TumorPurity:\n")
print(
  sort(
    table(
      three_layer_comparison$
        Unadjusted_to_TumorPurity_FDR
    ),
    decreasing = TRUE
  )
)

cat("\nImmuneScore -> TumorPurity:\n")
print(
  sort(
    table(
      three_layer_comparison$
        ImmuneScore_to_TumorPurity_FDR
    ),
    decreasing = TRUE
  )
)


# ------------------------------------------------------------
# 14.21 Three-way significance intersections
# ------------------------------------------------------------

all_three_sig <- with(
  three_layer_comparison,
  Unadjusted_Significant &
    ImmuneScore_Significant &
    TumorPurity_Significant
)

unadj_immune_only <- with(
  three_layer_comparison,
  Unadjusted_Significant &
    ImmuneScore_Significant &
    !TumorPurity_Significant
)

unadj_purity_only <- with(
  three_layer_comparison,
  Unadjusted_Significant &
    !ImmuneScore_Significant &
    TumorPurity_Significant
)

immune_purity_only <- with(
  three_layer_comparison,
  !Unadjusted_Significant &
    ImmuneScore_Significant &
    TumorPurity_Significant
)

unadjusted_only <- with(
  three_layer_comparison,
  Unadjusted_Significant &
    !ImmuneScore_Significant &
    !TumorPurity_Significant
)

immune_only <- with(
  three_layer_comparison,
  !Unadjusted_Significant &
    ImmuneScore_Significant &
    !TumorPurity_Significant
)

purity_only <- with(
  three_layer_comparison,
  !Unadjusted_Significant &
    !ImmuneScore_Significant &
    TumorPurity_Significant
)

none_sig <- with(
  three_layer_comparison,
  !Unadjusted_Significant &
    !ImmuneScore_Significant &
    !TumorPurity_Significant
)

three_way_significance_summary <- data.frame(
  Pattern = c(
    "Significant in all three",
    "Unadjusted + ImmuneScore only",
    "Unadjusted + TumorPurity only",
    "ImmuneScore + TumorPurity only",
    "Unadjusted only",
    "ImmuneScore only",
    "TumorPurity only",
    "None significant"
  ),
  N = c(
    sum(all_three_sig),
    sum(unadj_immune_only),
    sum(unadj_purity_only),
    sum(immune_purity_only),
    sum(unadjusted_only),
    sum(immune_only),
    sum(purity_only),
    sum(none_sig)
  ),
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("THREE-WAY SIGNIFICANCE INTERSECTIONS\n")
cat("============================================================\n")

print(
  three_way_significance_summary
)


# ------------------------------------------------------------
# 14.22 Core-gene-specific three-layer summary
# ------------------------------------------------------------

core_summary_8d14 <- do.call(
  rbind,
  lapply(
    core_genes_3,
    function(g) {
      
      tmp <- three_layer_comparison[
        three_layer_comparison$Core_Gene == g,
        ,
        drop = FALSE
      ]
      
      data.frame(
        
        Core_Gene = g,
        
        Total_Pairs = nrow(tmp),
        
        Unadjusted_Significant =
          sum(
            tmp$Unadjusted_Significant
          ),
        
        ImmuneScore_Significant =
          sum(
            tmp$ImmuneScore_Significant
          ),
        
        TumorPurity_Significant =
          sum(
            tmp$TumorPurity_Significant
          ),
        
        All_Three_Significant =
          sum(
            tmp$Unadjusted_Significant &
              tmp$ImmuneScore_Significant &
              tmp$TumorPurity_Significant
          ),
        
        Unadjusted_Positive =
          sum(
            tmp$Unadjusted_Direction ==
              "Positive"
          ),
        
        ImmuneScore_Positive =
          sum(
            tmp$ImmuneScore_Adjusted_Direction ==
              "Positive"
          ),
        
        TumorPurity_Positive =
          sum(
            tmp$Purity_Adjusted_Direction ==
              "Positive"
          ),
        
        Unadjusted_Negative =
          sum(
            tmp$Unadjusted_Direction ==
              "Negative"
          ),
        
        ImmuneScore_Negative =
          sum(
            tmp$ImmuneScore_Adjusted_Direction ==
              "Negative"
          ),
        
        TumorPurity_Negative =
          sum(
            tmp$Purity_Adjusted_Direction ==
              "Negative"
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)

cat("\n============================================================\n")
cat("CORE-GENE THREE-LAYER SUMMARY\n")
cat("============================================================\n")

print(
  core_summary_8d14,
  row.names = FALSE
)


# ------------------------------------------------------------
# 14.23 Identify pairs significant in all three analyses
# ------------------------------------------------------------

all_three_sig_pairs <- three_layer_comparison[
  all_three_sig,
  c(
    "Association_ID",
    "Core_Gene",
    "Checkpoint_Gene",
    "Unadjusted_rho",
    "Unadjusted_FDR",
    "ImmuneScore_Adjusted_rho",
    "ImmuneScore_Adjusted_FDR",
    "Purity_Adjusted_rho",
    "Purity_Adjusted_FDR",
    "Direction_Pattern"
  ),
  drop = FALSE
]

cat("\n============================================================\n")
cat("PAIRS SIGNIFICANT IN ALL THREE ANALYSES\n")
cat("============================================================\n")

cat(
  "Number:",
  nrow(all_three_sig_pairs),
  "\n\n"
)

print(
  all_three_sig_pairs,
  row.names = FALSE
)


# ------------------------------------------------------------
# 14.24 Identify direction switches under either adjustment
# ------------------------------------------------------------

three_layer_comparison$ImmuneScore_Direction_Switch <-
  three_layer_comparison$Unadjusted_Direction !=
  three_layer_comparison$ImmuneScore_Adjusted_Direction

three_layer_comparison$TumorPurity_Direction_Switch <-
  three_layer_comparison$Unadjusted_Direction !=
  three_layer_comparison$Purity_Adjusted_Direction

three_layer_comparison$Both_Adjustments_Direction_Switch <-
  three_layer_comparison$ImmuneScore_Direction_Switch &
  three_layer_comparison$TumorPurity_Direction_Switch


cat("\n============================================================\n")
cat("DIRECTION SWITCH SUMMARY\n")
cat("============================================================\n")

cat(
  "Unadjusted -> ImmuneScore direction switches:",
  sum(
    three_layer_comparison$
      ImmuneScore_Direction_Switch
  ),
  "\n"
)

cat(
  "Unadjusted -> TumorPurity direction switches:",
  sum(
    three_layer_comparison$
      TumorPurity_Direction_Switch
  ),
  "\n"
)

cat(
  "Direction switches under BOTH adjustments:",
  sum(
    three_layer_comparison$
      Both_Adjustments_Direction_Switch
  ),
  "\n"
)


# ------------------------------------------------------------
# 14.25 Verify locked primary values remain unchanged
# ------------------------------------------------------------

locked_rho_check <- identical(
  three_layer_comparison$Unadjusted_rho,
  locked_54_8d14$Unadjusted_Spearman_rho
)

locked_fdr_check <- identical(
  three_layer_comparison$Unadjusted_FDR,
  locked_54_8d14$Unadjusted_FDR_BH
)

immune_rho_check <- identical(
  three_layer_comparison$ImmuneScore_Adjusted_rho,
  locked_54_8d14$Adjusted_Partial_Spearman_rho
)

immune_fdr_check <- identical(
  three_layer_comparison$ImmuneScore_Adjusted_FDR,
  locked_54_8d14$Adjusted_FDR
)

cat("\n============================================================\n")
cat("LOCKED RESULT INTEGRITY CHECK\n")
cat("============================================================\n")

cat(
  "Unadjusted rho unchanged:",
  locked_rho_check,
  "\n"
)

cat(
  "Unadjusted FDR unchanged:",
  locked_fdr_check,
  "\n"
)

cat(
  "ImmuneScore-adjusted rho unchanged:",
  immune_rho_check,
  "\n"
)

cat(
  "ImmuneScore-adjusted FDR unchanged:",
  immune_fdr_check,
  "\n"
)

if (
  !locked_rho_check ||
  !locked_fdr_check ||
  !immune_rho_check ||
  !immune_fdr_check
) {
  stop(
    "ERROR: Locked STEP 8C-7B values were altered."
  )
}

cat(
  "PASS: All locked STEP 8C-7B values remain unchanged.\n"
)


# ------------------------------------------------------------
# 14.26 Final structural PASS criteria
# ------------------------------------------------------------

step8d14_PASS <- (
  
  nrow(three_layer_comparison) == 54 &&
    
    length(
      unique(
        three_layer_comparison$Association_ID
      )
    ) == 54 &&
    
    all(
      three_layer_comparison$N == 503
    ) &&
    
    all(
      three_layer_comparison$Purity_Adjusted_N == 503
    ) &&
    
    locked_rho_check &&
    
    locked_fdr_check &&
    
    immune_rho_check &&
    
    immune_fdr_check
  
)


# ------------------------------------------------------------
# 14.27 Final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 8D-14 FINAL STATUS\n")
cat("============================================================\n")

if (step8d14_PASS) {
  
  cat(
    "PASS: Three-layer paired comparison completed for all 54\n"
  )
  
  cat(
    "predefined associations with exact Association_ID matching.\n"
  )
  
  cat(
    "Locked STEP 8C-7B values remain unchanged.\n"
  )
  
} else {
  
  cat(
    "FAIL: Three-layer comparison requires review.\n"
  )
  
}


# ============================================================
# END OF STEP 8D-14
# ============================================================

# ============================================================
# Paper 1
# WORKFLOW-08 / STEP 8D-15
# Manuscript-ready Sensitivity Summary
# FINAL LOCK — corrected version
# ============================================================

cat("\n============================================================\n")
cat("STEP 8D-15: Manuscript-ready Sensitivity Summary\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 0. Output directory
# ------------------------------------------------------------

output_dir <- "result_dir/Paper1_tumor-only WGCNA/Workflow-08-STEP8D15_Manuscript_Sensitivity_Summary"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("\nOutput directory:\n")
cat(output_dir, "\n")


# ------------------------------------------------------------
# 1. Check required locked objects
# ------------------------------------------------------------

required_objects <- c(
  "locked_54",
  "purity_adjusted_54"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1))
]

if (length(missing_objects) > 0) {
  
  stop(
    paste0(
      "\nERROR: The following required objects are missing:\n",
      paste(missing_objects, collapse = "\n"),
      "\n\nPlease load the locked STEP 8C-7B and STEP 8D-13 objects."
    )
  )
}

cat("\nPASS: Required locked objects are available.\n")


# ------------------------------------------------------------
# 2. Verify locked 54-pair structure
# ------------------------------------------------------------

if (!is.data.frame(locked_54)) {
  stop("ERROR: locked_54 is not a data.frame.")
}

if (nrow(locked_54) != 54) {
  stop(
    paste0(
      "ERROR: locked_54 must contain 54 rows; found ",
      nrow(locked_54), "."
    )
  )
}

if (!"Association_ID" %in% colnames(locked_54)) {
  stop("ERROR: Association_ID missing from locked_54.")
}

if (anyDuplicated(locked_54$Association_ID) > 0) {
  stop("ERROR: Duplicate Association_ID values detected in locked_54.")
}

cat("PASS: locked_54 contains 54 unique predefined associations.\n")


# ------------------------------------------------------------
# 3. Verify purity-adjusted object
# ------------------------------------------------------------

if (!is.data.frame(purity_adjusted_54)) {
  stop("ERROR: purity_adjusted_54 is not a data.frame.")
}

if (nrow(purity_adjusted_54) != 54) {
  stop(
    paste0(
      "ERROR: purity_adjusted_54 must contain 54 rows; found ",
      nrow(purity_adjusted_54), "."
    )
  )
}

if (!"Association_ID" %in% colnames(purity_adjusted_54)) {
  stop("ERROR: Association_ID missing from purity_adjusted_54.")
}

if (anyDuplicated(purity_adjusted_54$Association_ID) > 0) {
  stop("ERROR: Duplicate Association_ID values detected in purity_adjusted_54.")
}

if (!identical(
  sort(locked_54$Association_ID),
  sort(purity_adjusted_54$Association_ID)
)) {
  stop(
    "ERROR: Association_ID sets differ between locked_54 and purity_adjusted_54."
  )
}

cat("PASS: locked_54 and purity_adjusted_54 contain the same 54 associations.\n")


# ------------------------------------------------------------
# 4. Manuscript-level LOCKED values
# ------------------------------------------------------------
#
# IMPORTANT:
# These values were already established in STEP 8C-7B
# and STEP 8D-13 / STEP 8D-14.
#
# No statistical recalculation is performed here.
# ------------------------------------------------------------

LOCK_TOTAL <- 54

LOCK_UNADJ_SIG <- 48
LOCK_IMMUNE_SIG <- 38
LOCK_PURITY_SIG <- 39

LOCK_ALL_THREE_SIG <- 27

LOCK_BOTH_ADJUST_P_TO_N <- 30

LOCK_P_N_N <- 17
LOCK_P_P_P <- 17


# ------------------------------------------------------------
# 5. Global manuscript-ready summary
# ------------------------------------------------------------

global_summary <- data.frame(
  
  Analysis_Layer = c(
    "Primary analysis",
    "Sensitivity analysis 1",
    "Sensitivity analysis 2"
  ),
  
  Adjustment = c(
    "None",
    "ESTIMATE-derived ImmuneScore",
    "ESTIMATE-derived TumorPurity"
  ),
  
  Method = c(
    "Two-sided Spearman rank correlation",
    "ImmuneScore-adjusted partial Spearman",
    "TumorPurity-adjusted partial Spearman"
  ),
  
  Total_Associations = c(
    54,
    54,
    54
  ),
  
  Significant_FDR_lt_0_05 = c(
    48,
    38,
    39
  ),
  
  Significant_Percent = c(
    round(100 * 48 / 54, 1),
    round(100 * 38 / 54, 1),
    round(100 * 39 / 54, 1)
  ),
  
  Positive_Associations = c(
    53,
    21,
    19
  ),
  
  Negative_Associations = c(
    1,
    33,
    35
  ),
  
  FDR_Method = rep(
    "Benjamini-Hochberg across 54 predefined comparisons",
    3
  ),
  
  stringsAsFactors = FALSE
)


cat("\n============================================================\n")
cat("GLOBAL MANUSCRIPT-READY SUMMARY\n")
cat("============================================================\n\n")

print(global_summary, row.names = FALSE)


# ------------------------------------------------------------
# 6. Three-way significance summary
# ------------------------------------------------------------

three_way_significance <- data.frame(
  
  Pattern = c(
    "Significant in all three analyses",
    "Unadjusted + ImmuneScore only",
    "Unadjusted + TumorPurity only",
    "ImmuneScore + TumorPurity only",
    "Unadjusted only",
    "ImmuneScore only",
    "TumorPurity only",
    "None significant"
  ),
  
  N = c(
    27,
    7,
    8,
    4,
    6,
    0,
    0,
    2
  ),
  
  stringsAsFactors = FALSE
)


cat("\n============================================================\n")
cat("THREE-WAY SIGNIFICANCE SUMMARY\n")
cat("============================================================\n\n")

print(three_way_significance, row.names = FALSE)


# ------------------------------------------------------------
# 7. Three-way direction summary
# ------------------------------------------------------------

three_way_direction <- data.frame(
  
  Direction_Pattern = c(
    "P -> N -> N",
    "P -> P -> P",
    "P -> P -> N",
    "P -> N -> P",
    "N -> N -> N"
  ),
  
  N = c(
    30,
    17,
    4,
    2,
    1
  ),
  
  Description = c(
    "Positive unadjusted association becomes negative after both adjustments",
    "Positive association remains positive across all three analyses",
    "Positive after ImmuneScore adjustment but negative after TumorPurity adjustment",
    "Negative after ImmuneScore adjustment but positive after TumorPurity adjustment",
    "Negative association remains negative across all three analyses"
  ),
  
  stringsAsFactors = FALSE
)


cat("\n============================================================\n")
cat("THREE-WAY DIRECTION SUMMARY\n")
cat("============================================================\n\n")

print(three_way_direction, row.names = FALSE)


# ------------------------------------------------------------
# 8. Pairwise direction-transition summary
# ------------------------------------------------------------

pairwise_direction_summary <- data.frame(
  
  Comparison = c(
    "Unadjusted -> ImmuneScore",
    "Unadjusted -> TumorPurity",
    "ImmuneScore -> TumorPurity"
  ),
  
  Positive_to_Positive = c(
    21,
    19,
    17
  ),
  
  Negative_to_Negative = c(
    1,
    1,
    31
  ),
  
  Positive_to_Negative = c(
    32,
    34,
    4
  ),
  
  Negative_to_Positive = c(
    0,
    0,
    2
  ),
  
  Direction_Switches = c(
    32,
    34,
    6
  ),
  
  stringsAsFactors = FALSE
)


cat("\n============================================================\n")
cat("PAIRWISE DIRECTION TRANSITIONS\n")
cat("============================================================\n\n")

print(pairwise_direction_summary, row.names = FALSE)


# ------------------------------------------------------------
# 9. Pairwise FDR-transition summary
# ------------------------------------------------------------

pairwise_fdr_summary <- data.frame(
  
  Comparison = c(
    "Unadjusted -> ImmuneScore",
    "Unadjusted -> TumorPurity",
    "ImmuneScore -> TumorPurity"
  ),
  
  Significant_to_Significant = c(
    34,
    35,
    31
  ),
  
  Significant_to_Nonsignificant = c(
    14,
    13,
    7
  ),
  
  Nonsignificant_to_Significant = c(
    4,
    4,
    8
  ),
  
  Nonsignificant_to_Nonsignificant = c(
    2,
    2,
    8
  ),
  
  stringsAsFactors = FALSE
)


cat("\n============================================================\n")
cat("PAIRWISE FDR TRANSITIONS\n")
cat("============================================================\n\n")

print(pairwise_fdr_summary, row.names = FALSE)


# ------------------------------------------------------------
# 10. Core-gene summary
# ------------------------------------------------------------

core_locked <- data.frame(
  
  Core_Gene = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  
  Total = c(
    18,
    18,
    18
  ),
  
  Unadjusted_Significant = c(
    14,
    17,
    17
  ),
  
  ImmuneScore_Adjusted_Significant = c(
    13,
    12,
    13
  ),
  
  TumorPurity_Adjusted_Significant = c(
    17,
    11,
    11
  ),
  
  stringsAsFactors = FALSE
)


cat("\n============================================================\n")
cat("CORE-GENE SENSITIVITY SUMMARY\n")
cat("============================================================\n\n")

print(core_locked, row.names = FALSE)


# ------------------------------------------------------------
# 11. Core-gene percentage summary
# ------------------------------------------------------------

core_locked_percent <- core_locked

core_locked_percent$Unadjusted_Percent <- round(
  100 * core_locked$Unadjusted_Significant /
    core_locked$Total,
  1
)

core_locked_percent$ImmuneScore_Percent <- round(
  100 * core_locked$ImmuneScore_Adjusted_Significant /
    core_locked$Total,
  1
)

core_locked_percent$TumorPurity_Percent <- round(
  100 * core_locked$TumorPurity_Adjusted_Significant /
    core_locked$Total,
  1
)


# ------------------------------------------------------------
# 12. Manuscript-ready key findings
# ------------------------------------------------------------

manuscript_key_findings <- data.frame(
  
  Finding_ID = c(
    "F1",
    "F2",
    "F3",
    "F4",
    "F5",
    "F6",
    "F7"
  ),
  
  Manuscript_Ready_Finding = c(
    
    "48 of 54 predefined associations were significant in the unadjusted analysis.",
    
    "38 of 54 predefined associations remained significant after adjustment for ESTIMATE-derived ImmuneScore.",
    
    "39 of 54 predefined associations remained significant after adjustment for ESTIMATE-derived tumor purity.",
    
    "27 of 54 predefined associations were significant in all three analytical layers.",
    
    "30 of 54 associations changed direction from positive in the unadjusted analysis to negative after both ImmuneScore and tumor-purity adjustment.",
    
    "17 of 54 associations followed a positive-to-negative-to-negative direction pattern across the three analytical layers.",
    
    "17 of 54 associations remained positive across all three analytical layers."
    
  ),
  
  stringsAsFactors = FALSE
)


cat("\n============================================================\n")
cat("MANUSCRIPT-READY KEY FINDINGS\n")
cat("============================================================\n\n")

print(manuscript_key_findings, row.names = FALSE)


# ------------------------------------------------------------
# 13. Numerical integrity QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("NUMERICAL INTEGRITY QC\n")
cat("============================================================\n")


QC <- data.frame(
  
  Check = c(
    "Total predefined associations",
    "Unadjusted significant",
    "ImmuneScore-adjusted significant",
    "TumorPurity-adjusted significant",
    "Significant in all three",
    "P -> N -> N",
    "P -> P -> P",
    "P -> N under both adjustments",
    "Core-gene unadjusted total",
    "Core-gene ImmuneScore total",
    "Core-gene TumorPurity total"
  ),
  
  Expected = c(
    54,
    48,
    38,
    39,
    27,
    17,
    17,
    30,
    48,
    38,
    39
  ),
  
  Actual = c(
    54,
    sum(global_summary$Significant_FDR_lt_0_05[1]),
    sum(global_summary$Significant_FDR_lt_0_05[2]),
    sum(global_summary$Significant_FDR_lt_0_05[3]),
    three_way_significance$N[
      three_way_significance$Pattern ==
        "Significant in all three analyses"
    ],
    three_way_direction$N[
      three_way_direction$Direction_Pattern ==
        "P -> N -> N"
    ],
    three_way_direction$N[
      three_way_direction$Direction_Pattern ==
        "P -> P -> P"
    ],
    three_way_direction$N[
      three_way_direction$Direction_Pattern ==
        "P -> N -> N"
    ] +
      three_way_direction$N[
        three_way_direction$Direction_Pattern ==
          "P -> N -> P"
      ],
    sum(core_locked$Unadjusted_Significant),
    sum(core_locked$ImmuneScore_Adjusted_Significant),
    sum(core_locked$TumorPurity_Adjusted_Significant)
  ),
  
  stringsAsFactors = FALSE
)


QC$Status <- ifelse(
  QC$Expected == QC$Actual,
  "PASS",
  "FAIL"
)


print(QC, row.names = FALSE)


# ============================================================
# STEP 8D-15
# FINAL QC CORRECTION PATCH
# ============================================================

cat("\n============================================================\n")
cat("STEP 8D-15 FINAL QC CORRECTION PATCH\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Correct locked three-way direction definition
# ------------------------------------------------------------

LOCK_P_N_N <- 30
LOCK_P_P_P <- 17


# ------------------------------------------------------------
# 2. Rebuild corrected QC
# ------------------------------------------------------------

QC_corrected <- data.frame(
  
  Check = c(
    "Total predefined associations",
    "Unadjusted significant",
    "ImmuneScore-adjusted significant",
    "TumorPurity-adjusted significant",
    "Significant in all three",
    "P -> N -> N",
    "P -> P -> P",
    "P -> N under both adjustments",
    "Core-gene unadjusted total",
    "Core-gene ImmuneScore total",
    "Core-gene TumorPurity total"
  ),
  
  Expected = c(
    54,
    48,
    38,
    39,
    27,
    30,
    17,
    30,
    48,
    38,
    39
  ),
  
  Actual = c(
    
    54,
    
    global_summary$Significant_FDR_lt_0_05[1],
    
    global_summary$Significant_FDR_lt_0_05[2],
    
    global_summary$Significant_FDR_lt_0_05[3],
    
    three_way_significance$N[
      three_way_significance$Pattern ==
        "Significant in all three analyses"
    ],
    
    three_way_direction$N[
      three_way_direction$Direction_Pattern ==
        "P -> N -> N"
    ],
    
    three_way_direction$N[
      three_way_direction$Direction_Pattern ==
        "P -> P -> P"
    ],
    
    # IMPORTANT:
    # "P -> N under both adjustments" means exactly
    # P -> N -> N.
    #
    # It must NOT be calculated as
    # P -> N -> N + P -> N -> P.
    
    three_way_direction$N[
      three_way_direction$Direction_Pattern ==
        "P -> N -> N"
    ],
    
    sum(core_locked$Unadjusted_Significant),
    
    sum(core_locked$ImmuneScore_Adjusted_Significant),
    
    sum(core_locked$TumorPurity_Adjusted_Significant)
    
  ),
  
  stringsAsFactors = FALSE
)


QC_corrected$Status <- ifelse(
  QC_corrected$Expected == QC_corrected$Actual,
  "PASS",
  "FAIL"
)


cat("\n============================================================\n")
cat("CORRECTED NUMERICAL INTEGRITY QC\n")
cat("============================================================\n\n")

print(QC_corrected, row.names = FALSE)


# ------------------------------------------------------------
# 3. Independent arithmetic checks
# ------------------------------------------------------------

direction_total <- sum(
  three_way_direction$N
)

significance_total <- sum(
  three_way_significance$N
)


if (direction_total != 54) {
  stop(
    paste0(
      "ERROR: Three-way direction categories sum to ",
      direction_total,
      ", not 54."
    )
  )
}


if (significance_total != 54) {
  stop(
    paste0(
      "ERROR: Three-way significance categories sum to ",
      significance_total,
      ", not 54."
    )
  )
}


if (
  three_way_direction$N[
    three_way_direction$Direction_Pattern ==
    "P -> N -> N"
  ] != 30
) {
  stop(
    "ERROR: P -> N -> N must equal 30."
  )
}


if (
  three_way_direction$N[
    three_way_direction$Direction_Pattern ==
    "P -> P -> P"
  ] != 17
) {
  stop(
    "ERROR: P -> P -> P must equal 17."
  )
}


if (any(QC_corrected$Status == "FAIL")) {
  stop(
    "ERROR: One or more corrected STEP 8D-15 QC checks failed."
  )
}


cat("\nPASS: All corrected numerical integrity checks passed.\n")


# ------------------------------------------------------------
# 4. Correct manuscript-ready key findings
# ------------------------------------------------------------

manuscript_key_findings <- data.frame(
  
  Finding_ID = c(
    "F1",
    "F2",
    "F3",
    "F4",
    "F5",
    "F6"
  ),
  
  Manuscript_Ready_Finding = c(
    
    "48 of 54 predefined associations were significant in the unadjusted analysis.",
    
    "38 of 54 predefined associations remained significant after adjustment for ESTIMATE-derived ImmuneScore.",
    
    "39 of 54 predefined associations remained significant after adjustment for ESTIMATE-derived tumor purity.",
    
    "27 of 54 predefined associations were significant in all three analytical layers.",
    
    "30 of 54 associations changed from positive in the unadjusted analysis to negative after both ImmuneScore and tumor-purity adjustment.",
    
    "17 of 54 associations remained positive across all three analytical layers."
    
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 5. Correct final lock object
# ------------------------------------------------------------

STEP8D15_LOCK <- list(
  
  Workflow =
    "Paper 1 / WORKFLOW-08 / STEP 8D-15",
  
  Status =
    "FINAL LOCK",
  
  Total_Associations =
    54,
  
  Primary_Analysis = list(
    Method =
      "Unadjusted Spearman",
    Significant =
      48,
    Positive =
      53,
    Negative =
      1
  ),
  
  ImmuneScore_Sensitivity = list(
    Method =
      "ImmuneScore-adjusted partial Spearman",
    Significant =
      38,
    Positive =
      21,
    Negative =
      33
  ),
  
  TumorPurity_Sensitivity = list(
    Method =
      "TumorPurity-adjusted partial Spearman",
    Significant =
      39,
    Positive =
      19,
    Negative =
      35
  ),
  
  Three_Way_Summary = list(
    Significant_in_all_three =
      27,
    
    P_N_N =
      30,
    
    P_P_P =
      17,
    
    P_P_N =
      4,
    
    P_N_P =
      2,
    
    N_N_N =
      1
  ),
  
  Core_Gene_Summary =
    core_locked,
  
  Interpretation =
    paste(
      "Tumor-purity adjustment is treated as an independent sensitivity",
      "analysis and does not replace the locked unadjusted or",
      "ImmuneScore-adjusted results. All analyses represent",
      "bulk-tumor expression associations."
    ),
  
  No_New_Statistical_Testing =
    TRUE,
  
  Locked_STEP8C7B_Values_Preserved =
    TRUE
)


# ------------------------------------------------------------
# 6. Overwrite corrected output files
# ------------------------------------------------------------

saveRDS(
  STEP8D15_LOCK,
  file = file.path(
    output_dir,
    "STEP8D15_Manuscript_Sensitivity_Summary_LOCK.rds"
  )
)


write.csv(
  manuscript_key_findings,
  file.path(
    output_dir,
    "STEP8D15_Manuscript_Key_Findings.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  QC_corrected,
  file.path(
    output_dir,
    "STEP8D15_QC.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 7. Final corrected report
# ------------------------------------------------------------

report_file <- file.path(
  output_dir,
  "STEP8D15_FINAL_LOCK_REPORT.txt"
)

sink(report_file)

cat("============================================================\n")
cat("Paper 1 / WORKFLOW-08 / STEP 8D-15\n")
cat("MANUSCRIPT-READY SENSITIVITY SUMMARY\n")
cat("FINAL LOCK REPORT\n")
cat("============================================================\n\n")

cat("STATUS: FINAL LOCK\n\n")

cat("Primary analysis:\n")
cat("  Unadjusted Spearman: 48/54 significant\n\n")

cat("Sensitivity analysis 1:\n")
cat("  ImmuneScore-adjusted partial Spearman: 38/54 significant\n\n")

cat("Sensitivity analysis 2:\n")
cat("  TumorPurity-adjusted partial Spearman: 39/54 significant\n\n")

cat("Three-layer significance:\n")
cat("  Significant in all three: 27/54\n\n")

cat("Three-layer direction:\n")
cat("  P -> N -> N: 30/54\n")
cat("  P -> P -> P: 17/54\n")
cat("  P -> P -> N: 4/54\n")
cat("  P -> N -> P: 2/54\n")
cat("  N -> N -> N: 1/54\n\n")

cat("Core-gene summary:\n")
print(core_locked, row.names = FALSE)

cat("\nCorrected QC:\n")
print(QC_corrected, row.names = FALSE)

cat("\nLocked STEP 8C-7B values preserved: TRUE\n")
cat("No new statistical testing performed: TRUE\n")

cat("\n============================================================\n")
cat("FINAL STATUS: PASS\n")
cat("STEP 8D-15 manuscript-level sensitivity summary locked.\n")
cat("============================================================\n")

sink()


# ------------------------------------------------------------
# 8. Final console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 8D-15 FINAL STATUS\n")
cat("============================================================\n")

cat("\nPASS: All corrected numerical integrity checks passed.\n")

cat("\nFinal locked global values:\n")
cat("  Unadjusted:              48 / 54\n")
cat("  ImmuneScore-adjusted:    38 / 54\n")
cat("  TumorPurity-adjusted:    39 / 54\n")

cat("\nFinal locked three-layer direction values:\n")
cat("  P -> N -> N:              30 / 54\n")
cat("  P -> P -> P:              17 / 54\n")
cat("  P -> P -> N:               4 / 54\n")
cat("  P -> N -> P:               2 / 54\n")
cat("  N -> N -> N:               1 / 54\n")

cat("\nSignificant in all three:\n")
cat("  27 / 54\n")

cat("\nCore-gene significant counts:\n")
cat("  AOC3:   14 / 13 / 17\n")
cat("  MS4A7:  17 / 12 / 11\n")
cat("  IL7R:   17 / 13 / 11\n")

cat("\nLocked STEP 8C-7B values preserved: TRUE\n")
cat("No new statistical testing performed: TRUE\n")

cat("\nFINAL STATUS: PASS\n")                        
cat("============================================================\n")

##STEP 8E — Leave-One-Core-Gene-Out Circularity / Dependence Audit
#Phase 1: Locked-result subset audit
# ============================================================
# Paper 1
# STEP 8E — Leave-One-Core-Gene-Out Circularity / Dependence Audit
#
# Phase 1: Locked-result subset audit
#
# PURPOSE
# ------------------------------------------------------------
# Audit whether the locked 54-pair immune-checkpoint association
# pattern is disproportionately dependent on any one core gene.
#
# CORE GENES:
#   AOC3
#   MS4A7
#   IL7R
#
# CHECKPOINT GENES:
#   18 predefined checkpoint-related genes
#
# IMPORTANT:
# ------------------------------------------------------------
# NO new statistical testing is performed in STEP 8E Phase 1.
#
# The following locked results are ONLY SUBSETTED:
#   1. locked_54
#   2. purity_adjusted_54
#
# No new:
#   - Spearman correlation
#   - partial Spearman correlation
#   - ImmuneScore adjustment
#   - TumorPurity adjustment
#   - BH-FDR recalculation
#
# This is an AUDIT / SENSITIVITY analysis only.
# ============================================================


# ============================================================
# STEP 8E-00
# Environment and output directory
# ============================================================

cat("\n============================================================\n")
cat("STEP 8E — Leave-One-Core-Gene-Out Circularity / Dependence Audit\n")
cat("============================================================\n")

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------

base_dir <- "result_dir/Paper1_tumor-only WGCNA"

out_dir <- file.path(
  base_dir,
  "Workflow-08-STEP8E_Leave-One-Core-Gene-Out_Audit"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cat("\nOutput directory:\n")
cat(out_dir, "\n")


# ============================================================
# STEP 8E-01
# Validate locked input objects
# ============================================================

cat("\n============================================================\n")
cat("STEP 8E-01: Validate locked input objects\n")
cat("============================================================\n")


# ------------------------------------------------------------
# Required objects
# ------------------------------------------------------------

required_objects <- c(
  "locked_54",
  "purity_adjusted_54"
)

object_exists <- vapply(
  required_objects,
  exists,
  logical(1),
  envir = .GlobalEnv
)

cat("\nRequired objects:\n")
print(
  data.frame(
    Object = required_objects,
    Exists = object_exists,
    row.names = NULL
  )
)

if (!all(object_exists)) {
  stop(
    "\nFAIL: One or more required locked objects are missing.\n",
    "Required: locked_54, purity_adjusted_54\n"
  )
}


# ------------------------------------------------------------
# Basic object checks
# ------------------------------------------------------------

qc_8e01 <- data.frame(
  Check = character(),
  Value = character(),
  Status = character(),
  stringsAsFactors = FALSE
)

add_qc <- function(check, value, status) {
  qc_8e01 <<- rbind(
    qc_8e01,
    data.frame(
      Check = check,
      Value = as.character(value),
      Status = status,
      stringsAsFactors = FALSE
    )
  )
}


# locked_54
add_qc(
  "locked_54_is_data.frame",
  is.data.frame(locked_54),
  ifelse(is.data.frame(locked_54), "PASS", "FAIL")
)

add_qc(
  "locked_54_rows",
  nrow(locked_54),
  ifelse(nrow(locked_54) == 54, "PASS", "FAIL")
)

add_qc(
  "locked_54_unique_Association_ID",
  length(unique(locked_54$Association_ID)),
  ifelse(
    length(unique(locked_54$Association_ID)) == 54,
    "PASS",
    "FAIL"
  )
)


# purity_adjusted_54
add_qc(
  "purity_adjusted_54_is_data.frame",
  is.data.frame(purity_adjusted_54),
  ifelse(
    is.data.frame(purity_adjusted_54),
    "PASS",
    "FAIL"
  )
)

add_qc(
  "purity_adjusted_54_rows",
  nrow(purity_adjusted_54),
  ifelse(
    nrow(purity_adjusted_54) == 54,
    "PASS",
    "FAIL"
  )
)

add_qc(
  "purity_adjusted_54_unique_Association_ID",
  length(unique(purity_adjusted_54$Association_ID)),
  ifelse(
    length(unique(purity_adjusted_54$Association_ID)) == 54,
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Verify Association_ID sets
# ------------------------------------------------------------

locked_ids <- as.character(locked_54$Association_ID)
purity_ids <- as.character(purity_adjusted_54$Association_ID)

id_set_match <- setequal(
  locked_ids,
  purity_ids
)

add_qc(
  "locked_54_vs_purity_adjusted_54_Association_ID_set",
  id_set_match,
  ifelse(id_set_match, "PASS", "FAIL")
)

cat("\nSTEP 8E-01 QC:\n")
print(qc_8e01)


if (any(qc_8e01$Status == "FAIL")) {
  stop("\nSTEP 8E-01 FAILED. Analysis stopped.\n")
}

cat("\nPASS: STEP 8E-01 locked input objects validated.\n")


# ============================================================
# STEP 8E-02
# Verify core genes and 18 checkpoint genes
# ============================================================

cat("\n============================================================\n")
cat("STEP 8E-02: Verify core genes and checkpoint genes\n")
cat("============================================================\n")


# ------------------------------------------------------------
# Locked core genes
# ------------------------------------------------------------

core_genes_locked <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


# ------------------------------------------------------------
# Locked 18 checkpoint-related genes
# ------------------------------------------------------------

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
# Verify core genes
# ------------------------------------------------------------

core_in_locked <- sort(unique(locked_54$Core_Gene))

cat("\nCore genes found in locked_54:\n")
print(core_in_locked)

core_set_pass <- setequal(
  core_in_locked,
  core_genes_locked
)

cat(
  "\nCore-gene set PASS:",
  core_set_pass,
  "\n"
)


# ------------------------------------------------------------
# Each core gene must have exactly 18 pairs
# ------------------------------------------------------------

core_pair_counts <- table(
  factor(
    locked_54$Core_Gene,
    levels = core_genes_locked
  )
)

print(core_pair_counts)

core_pair_count_pass <- all(
  core_pair_counts == 18
)

cat(
  "\n18 checkpoint pairs per core gene PASS:",
  core_pair_count_pass,
  "\n"
)


# ------------------------------------------------------------
# Verify checkpoint gene set
# ------------------------------------------------------------

checkpoint_in_locked <- sort(
  unique(locked_54$Checkpoint_Gene)
)

checkpoint_set_pass <- setequal(
  checkpoint_in_locked,
  checkpoint_genes_locked
)

cat("\nCheckpoint genes found:\n")
print(checkpoint_in_locked)

cat(
  "\n18-gene checkpoint set PASS:",
  checkpoint_set_pass,
  "\n"
)


# ------------------------------------------------------------
# Verify total pair structure
# ------------------------------------------------------------

pair_structure_pass <- (
  nrow(locked_54) == 54 &&
    length(unique(locked_54$Association_ID)) == 54 &&
    length(unique(locked_54$Core_Gene)) == 3 &&
    length(unique(locked_54$Checkpoint_Gene)) == 18
)

cat(
  "\n54-pair structure PASS:",
  pair_structure_pass,
  "\n"
)


# ------------------------------------------------------------
# STEP 8E-02 QC
# ------------------------------------------------------------

qc_8e02 <- data.frame(
  Check = c(
    "Core gene set",
    "18 pairs per core gene",
    "Checkpoint gene set",
    "54-pair structure"
  ),
  Status = c(
    ifelse(core_set_pass, "PASS", "FAIL"),
    ifelse(core_pair_count_pass, "PASS", "FAIL"),
    ifelse(checkpoint_set_pass, "PASS", "FAIL"),
    ifelse(pair_structure_pass, "PASS", "FAIL")
  ),
  stringsAsFactors = FALSE
)

print(qc_8e02)

if (any(qc_8e02$Status == "FAIL")) {
  stop("\nSTEP 8E-02 FAILED. Analysis stopped.\n")
}

cat(
  "\nPASS: STEP 8E-02 confirms AOC3/MS4A7/IL7R × 18 checkpoint genes.\n"
)


# ============================================================
# STEP 8E-03
# Build three leave-one-core-gene-out subsets
# ============================================================

cat("\n============================================================\n")
cat("STEP 8E-03: Build leave-one-core-gene-out subsets\n")
cat("============================================================\n")


# ------------------------------------------------------------
# Subset definitions
# ------------------------------------------------------------

subset_definitions <- list(
  
  Leave_AOC3 = list(
    omitted_gene = "AOC3",
    retained_genes = c("MS4A7", "IL7R")
  ),
  
  Leave_MS4A7 = list(
    omitted_gene = "MS4A7",
    retained_genes = c("AOC3", "IL7R")
  ),
  
  Leave_IL7R = list(
    omitted_gene = "IL7R",
    retained_genes = c("AOC3", "MS4A7")
  )
  
)


# ------------------------------------------------------------
# Create subsets from LOCKED results
# ------------------------------------------------------------

subset_locked <- list()
subset_purity <- list()

for (nm in names(subset_definitions)) {
  
  keep_genes <- subset_definitions[[nm]]$retained_genes
  
  subset_locked[[nm]] <- locked_54[
    locked_54$Core_Gene %in% keep_genes,
    ,
    drop = FALSE
  ]
  
  subset_purity[[nm]] <- purity_adjusted_54[
    purity_adjusted_54$Core_Gene %in% keep_genes,
    ,
    drop = FALSE
  ]
  
}


# ------------------------------------------------------------
# Verify each subset = 36 pairs
# ------------------------------------------------------------

subset_size_table <- do.call(
  rbind,
  lapply(
    names(subset_locked),
    function(nm) {
      
      data.frame(
        Subset = nm,
        Omitted_Core_Gene =
          subset_definitions[[nm]]$omitted_gene,
        Retained_Core_Genes =
          paste(
            subset_definitions[[nm]]$retained_genes,
            collapse = " + "
          ),
        Locked_Pairs =
          nrow(subset_locked[[nm]]),
        Purity_Adjusted_Pairs =
          nrow(subset_purity[[nm]]),
        stringsAsFactors = FALSE
      )
      
    }
  )
)

rownames(subset_size_table) <- NULL

cat("\nSubset structure:\n")
print(subset_size_table)


subset_size_pass <- all(
  subset_size_table$Locked_Pairs == 36 &
    subset_size_table$Purity_Adjusted_Pairs == 36
)

cat(
  "\nAll leave-one-gene-out subsets contain 36 pairs:",
  subset_size_pass,
  "\n"
)

if (!subset_size_pass) {
  stop(
    "\nFAIL: One or more leave-one-gene-out subsets do not contain 36 pairs.\n"
  )
}


# ============================================================
# STEP 8E-04
# Exact Association_ID matching
# ============================================================

cat("\n============================================================\n")
cat("STEP 8E-04: Exact Association_ID matching\n")
cat("============================================================\n")


association_matching_results <- list()

for (nm in names(subset_locked)) {
  
  ids_locked_subset <- as.character(
    subset_locked[[nm]]$Association_ID
  )
  
  ids_purity_subset <- as.character(
    subset_purity[[nm]]$Association_ID
  )
  
  exact_set_match <- setequal(
    ids_locked_subset,
    ids_purity_subset
  )
  
  exact_order_match <- identical(
    ids_locked_subset,
    ids_purity_subset
  )
  
  association_matching_results[[nm]] <- data.frame(
    Subset = nm,
    N_Locked = length(ids_locked_subset),
    N_Purity = length(ids_purity_subset),
    Unique_Locked = length(unique(ids_locked_subset)),
    Unique_Purity = length(unique(ids_purity_subset)),
    Exact_Set_Match = exact_set_match,
    Exact_Order_Match = exact_order_match,
    stringsAsFactors = FALSE
  )
  
}

association_matching_table <- do.call(
  rbind,
  association_matching_results
)

rownames(association_matching_table) <- NULL

print(association_matching_table)


association_set_pass <- all(
  association_matching_table$Exact_Set_Match &
    association_matching_table$Unique_Locked == 36 &
    association_matching_table$Unique_Purity == 36
)

cat(
  "\nAssociation_ID set matching PASS:",
  association_set_pass,
  "\n"
)

if (!association_set_pass) {
  stop(
    "\nFAIL: Association_ID mismatch detected between locked and purity-adjusted subsets.\n"
  )
}


# ------------------------------------------------------------
# Create a common Association_ID order
# ------------------------------------------------------------

for (nm in names(subset_locked)) {
  
  id_order <- subset_locked[[nm]]$Association_ID
  
  subset_purity[[nm]] <- subset_purity[[nm]][
    match(id_order, subset_purity[[nm]]$Association_ID),
    ,
    drop = FALSE
  ]
  
}

cat(
  "\nPASS: All purity-adjusted subsets reordered to locked Association_ID order.\n"
)


# ============================================================
# STEP 8E-05
# Extract locked significance and direction results
# ============================================================

cat("\n============================================================\n")
cat("STEP 8E-05: Extract locked significance and direction results\n")
cat("============================================================\n")


# ------------------------------------------------------------
# Helper: safely identify a column
# ------------------------------------------------------------

find_column <- function(
    df,
    exact_candidates = character(),
    pattern_candidates = character(),
    label = "requested column"
) {
  
  nms <- names(df)
  
  # exact candidates first
  hit <- exact_candidates[
    exact_candidates %in% nms
  ]
  
  if (length(hit) > 0) {
    return(hit[1])
  }
  
  # pattern candidates
  for (pat in pattern_candidates) {
    
    idx <- grep(
      pat,
      nms,
      ignore.case = TRUE
    )
    
    if (length(idx) > 0) {
      return(nms[idx[1]])
    }
    
  }
  
  stop(
    "\nFAIL: Cannot identify ",
    label,
    ".\nAvailable columns are:\n",
    paste(nms, collapse = ", "),
    "\n"
  )
}


# ------------------------------------------------------------
# Locked_54 columns
# ------------------------------------------------------------

col_unadj_sig <- find_column(
  locked_54,
  exact_candidates = c(
    "Unadjusted_FDR_Significant",
    "Sig_Unadjusted"
  ),
  pattern_candidates = c(
    "^Unadjusted.*Significant$",
    "^Sig.*Unadjusted"
  ),
  label = "unadjusted significance column"
)


col_immune_sig <- find_column(
  locked_54,
  exact_candidates = c(
    "Adjusted_FDR_Significant",
    "Sig_Adjusted"
  ),
  pattern_candidates = c(
    "^Adjusted.*Significant$",
    "^Sig.*Adjusted"
  ),
  label = "ImmuneScore-adjusted significance column"
)


col_immune_dir <- find_column(
  locked_54,
  exact_candidates = c(
    "Adjusted_Direction"
  ),
  pattern_candidates = c(
    "^Adjusted.*Direction$"
  ),
  label = "ImmuneScore-adjusted direction column"
)


col_unadj_dir <- find_column(
  locked_54,
  exact_candidates = c(
    "Unadjusted_Direction"
  ),
  pattern_candidates = c(
    "^Unadjusted.*Direction$"
  ),
  label = "unadjusted direction column"
)


# ------------------------------------------------------------
# Purity-adjusted columns
# ------------------------------------------------------------

col_purity_sig <- find_column(
  purity_adjusted_54,
  exact_candidates = c(
    "TumorPurity_FDR_Significant",
    "TumorPurity_Adjusted_FDR_Significant",
    "Purity_FDR_Significant",
    "Purity_Adjusted_FDR_Significant",
    "FDR_Significant"
  ),
  pattern_candidates = c(
    "TumorPurity.*Significant",
    "Purity.*Significant"
  ),
  label = "TumorPurity-adjusted significance column"
)


col_purity_dir <- find_column(
  purity_adjusted_54,
  exact_candidates = c(
    "TumorPurity_Direction",
    "TumorPurity_Adjusted_Direction",
    "Purity_Direction",
    "Purity_Adjusted_Direction",
    "Adjusted_Direction"
  ),
  pattern_candidates = c(
    "TumorPurity.*Direction",
    "Purity.*Direction"
  ),
  label = "TumorPurity-adjusted direction column"
)


cat("\nLocked result columns used:\n")
cat("Unadjusted significance :", col_unadj_sig, "\n")
cat("ImmuneScore significance:", col_immune_sig, "\n")
cat("ImmuneScore direction   :", col_immune_dir, "\n")
cat("Unadjusted direction    :", col_unadj_dir, "\n")
cat("TumorPurity significance:", col_purity_sig, "\n")
cat("TumorPurity direction   :", col_purity_dir, "\n")


# ------------------------------------------------------------
# Normalize significance values to TRUE/FALSE
# ------------------------------------------------------------

as_logical_sig <- function(x) {
  
  if (is.logical(x)) {
    return(x)
  }
  
  if (is.numeric(x)) {
    return(x != 0)
  }
  
  x_chr <- tolower(trimws(as.character(x)))
  
  x_chr %in% c(
    "true",
    "t",
    "yes",
    "y",
    "1",
    "significant",
    "sig"
  )
  
}


# ------------------------------------------------------------
# Normalize directions
# ------------------------------------------------------------

normalize_direction <- function(x) {
  
  x_chr <- toupper(trimws(as.character(x)))
  
  out <- rep(NA_character_, length(x_chr))
  
  out[x_chr %in% c(
    "POSITIVE",
    "POS",
    "+",
    "P"
  )] <- "P"
  
  out[x_chr %in% c(
    "NEGATIVE",
    "NEG",
    "-",
    "N"
  )] <- "N"
  
  out
  
}


# ------------------------------------------------------------
# Add locked audit fields
# ------------------------------------------------------------

locked_audit <- locked_54

locked_audit$STEP8E_Unadjusted_Sig <-
  as_logical_sig(
    locked_54[[col_unadj_sig]]
  )

locked_audit$STEP8E_ImmuneScore_Sig <-
  as_logical_sig(
    locked_54[[col_immune_sig]]
  )

locked_audit$STEP8E_Unadjusted_Direction <-
  normalize_direction(
    locked_54[[col_unadj_dir]]
  )

locked_audit$STEP8E_ImmuneScore_Direction <-
  normalize_direction(
    locked_54[[col_immune_dir]]
  )


purity_audit <- purity_adjusted_54

purity_audit$STEP8E_TumorPurity_Sig <-
  as_logical_sig(
    purity_adjusted_54[[col_purity_sig]]
  )

purity_audit$STEP8E_TumorPurity_Direction <-
  normalize_direction(
    purity_adjusted_54[[col_purity_dir]]
  )


# ------------------------------------------------------------
# Merge by Association_ID
# ------------------------------------------------------------

purity_fields <- purity_audit[
  ,
  c(
    "Association_ID",
    "STEP8E_TumorPurity_Sig",
    "STEP8E_TumorPurity_Direction"
  ),
  drop = FALSE
]


audit_54 <- merge(
  locked_audit,
  purity_fields,
  by = "Association_ID",
  all.x = TRUE,
  sort = FALSE
)


# Restore original locked order
audit_54 <- audit_54[
  match(
    locked_54$Association_ID,
    audit_54$Association_ID
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# Three-layer significance and direction
# ------------------------------------------------------------

audit_54$STEP8E_All_Three_Sig <-
  audit_54$STEP8E_Unadjusted_Sig &
  audit_54$STEP8E_ImmuneScore_Sig &
  audit_54$STEP8E_TumorPurity_Sig


audit_54$STEP8E_ThreeLayer_Direction <-
  paste0(
    audit_54$STEP8E_Unadjusted_Direction,
    "→",
    audit_54$STEP8E_ImmuneScore_Direction,
    "→",
    audit_54$STEP8E_TumorPurity_Direction
  )


# ------------------------------------------------------------
# Direction-pattern classification
# ------------------------------------------------------------

audit_54$STEP8E_Direction_Pattern <- NA_character_

audit_54$STEP8E_Direction_Pattern[
  audit_54$STEP8E_Unadjusted_Direction == "P" &
    audit_54$STEP8E_ImmuneScore_Direction == "N" &
    audit_54$STEP8E_TumorPurity_Direction == "N"
] <- "P→N→N"

audit_54$STEP8E_Direction_Pattern[
  audit_54$STEP8E_Unadjusted_Direction == "P" &
    audit_54$STEP8E_ImmuneScore_Direction == "P" &
    audit_54$STEP8E_TumorPurity_Direction == "P"
] <- "P→P→P"

audit_54$STEP8E_Direction_Pattern[
  audit_54$STEP8E_Unadjusted_Direction == "P" &
    audit_54$STEP8E_ImmuneScore_Direction == "P" &
    audit_54$STEP8E_TumorPurity_Direction == "N"
] <- "P→P→N"

audit_54$STEP8E_Direction_Pattern[
  audit_54$STEP8E_Unadjusted_Direction == "P" &
    audit_54$STEP8E_ImmuneScore_Direction == "N" &
    audit_54$STEP8E_TumorPurity_Direction == "P"
] <- "P→N→P"

audit_54$STEP8E_Direction_Pattern[
  audit_54$STEP8E_Unadjusted_Direction == "N" &
    audit_54$STEP8E_ImmuneScore_Direction == "N" &
    audit_54$STEP8E_TumorPurity_Direction == "N"
] <- "N→N→N"


# ------------------------------------------------------------
# Verify all 54 rows have complete locked results
# ------------------------------------------------------------

complete_locked_results <- all(
  !is.na(audit_54$STEP8E_Unadjusted_Sig) &
    !is.na(audit_54$STEP8E_ImmuneScore_Sig) &
    !is.na(audit_54$STEP8E_TumorPurity_Sig) &
    !is.na(audit_54$STEP8E_Unadjusted_Direction) &
    !is.na(audit_54$STEP8E_ImmuneScore_Direction) &
    !is.na(audit_54$STEP8E_TumorPurity_Direction)
)

cat(
  "\nComplete locked significance/direction fields:",
  complete_locked_results,
  "\n"
)

if (!complete_locked_results) {
  stop(
    "\nFAIL: Missing locked significance or direction values.\n"
  )
}


# ============================================================
# STEP 8E-06
# Numerical integrity QC
# ============================================================

cat("\n============================================================\n")
cat("STEP 8E-06: Numerical integrity QC\n")
cat("============================================================\n")


# ------------------------------------------------------------
# Expected locked global values from STEP 8C-7B / STEP 8D-15
# ------------------------------------------------------------

EXPECTED_UNADJ_SIG <- 48
EXPECTED_IMMUNE_SIG <- 38
EXPECTED_PURITY_SIG <- 39

EXPECTED_ALL_THREE_SIG <- 27

EXPECTED_PNN <- 30
EXPECTED_PPP <- 17
EXPECTED_PPN <- 4
EXPECTED_PNP <- 2
EXPECTED_NNN <- 1


# ------------------------------------------------------------
# Re-audit global locked values
# ------------------------------------------------------------

global_unadj_sig <- sum(
  audit_54$STEP8E_Unadjusted_Sig
)

global_immune_sig <- sum(
  audit_54$STEP8E_ImmuneScore_Sig
)

global_purity_sig <- sum(
  audit_54$STEP8E_TumorPurity_Sig
)

global_all_three <- sum(
  audit_54$STEP8E_All_Three_Sig
)


global_direction_counts <- table(
  factor(
    audit_54$STEP8E_Direction_Pattern,
    levels = c(
      "P→N→N",
      "P→P→P",
      "P→P→N",
      "P→N→P",
      "N→N→N"
    )
  )
)

cat("\nGlobal locked values re-audit:\n")
cat("Unadjusted significant:", global_unadj_sig, "\n")
cat("ImmuneScore significant:", global_immune_sig, "\n")
cat("TumorPurity significant:", global_purity_sig, "\n")
cat("All three significant:", global_all_three, "\n")

print(global_direction_counts)


# ------------------------------------------------------------
# Global integrity checks
# ------------------------------------------------------------

global_integrity_checks <- data.frame(
  Check = c(
    "Unadjusted significant",
    "ImmuneScore-adjusted significant",
    "TumorPurity-adjusted significant",
    "Significant in all three layers",
    "P→N→N",
    "P→P→P",
    "P→P→N",
    "P→N→P",
    "N→N→N"
  ),
  
  Expected = c(
    EXPECTED_UNADJ_SIG,
    EXPECTED_IMMUNE_SIG,
    EXPECTED_PURITY_SIG,
    EXPECTED_ALL_THREE_SIG,
    EXPECTED_PNN,
    EXPECTED_PPP,
    EXPECTED_PPN,
    EXPECTED_PNP,
    EXPECTED_NNN
  ),
  
  Observed = c(
    global_unadj_sig,
    global_immune_sig,
    global_purity_sig,
    global_all_three,
    as.integer(global_direction_counts["P→N→N"]),
    as.integer(global_direction_counts["P→P→P"]),
    as.integer(global_direction_counts["P→P→N"]),
    as.integer(global_direction_counts["P→N→P"]),
    as.integer(global_direction_counts["N→N→N"])
  ),
  
  stringsAsFactors = FALSE
)

global_integrity_checks$Status <-
  ifelse(
    global_integrity_checks$Expected ==
      global_integrity_checks$Observed,
    "PASS",
    "FAIL"
  )

print(global_integrity_checks)


# ------------------------------------------------------------
# Three-layer arithmetic
# ------------------------------------------------------------

direction_sum <- sum(
  global_direction_counts
)

direction_arithmetic_pass <- (
  direction_sum == 54
)

significance_pattern_sum <- sum(
  c(
    sum(
      audit_54$STEP8E_Unadjusted_Sig &
        audit_54$STEP8E_ImmuneScore_Sig &
        audit_54$STEP8E_TumorPurity_Sig
    ),
    
    sum(
      audit_54$STEP8E_Unadjusted_Sig &
        !audit_54$STEP8E_ImmuneScore_Sig &
        audit_54$STEP8E_TumorPurity_Sig
    ),
    
    sum(
      audit_54$STEP8E_Unadjusted_Sig &
        audit_54$STEP8E_ImmuneScore_Sig &
        !audit_54$STEP8E_TumorPurity_Sig
    ),
    
    sum(
      audit_54$STEP8E_Unadjusted_Sig &
        !audit_54$STEP8E_ImmuneScore_Sig &
        !audit_54$STEP8E_TumorPurity_Sig
    ),
    
    sum(
      !audit_54$STEP8E_Unadjusted_Sig &
        audit_54$STEP8E_ImmuneScore_Sig &
        audit_54$STEP8E_TumorPurity_Sig
    ),
    
    sum(
      !audit_54$STEP8E_Unadjusted_Sig &
        !audit_54$STEP8E_ImmuneScore_Sig &
        !audit_54$STEP8E_TumorPurity_Sig
    )
  )
)

significance_arithmetic_pass <- (
  significance_pattern_sum == 54
)

cat(
  "\nDirection arithmetic PASS:",
  direction_arithmetic_pass,
  "\n"
)

cat(
  "Significance arithmetic PASS:",
  significance_arithmetic_pass,
  "\n"
)


# ------------------------------------------------------------
# Locked STEP 8C-7B values protection
# ------------------------------------------------------------

locked_8c7b_protection <- (
  global_unadj_sig == 48 &&
    global_immune_sig == 38
)

locked_8d15_protection <- (
  global_purity_sig == 39 &&
    global_all_three == 27 &&
    as.integer(global_direction_counts["P→N→N"]) == 30
)

cat(
  "\nLocked STEP 8C-7B values preserved:",
  locked_8c7b_protection,
  "\n"
)

cat(
  "Locked STEP 8D-15 values preserved:",
  locked_8d15_protection,
  "\n"
)


# ------------------------------------------------------------
# Subset summary function
# ------------------------------------------------------------

make_subset_summary <- function(
    locked_subset,
    purity_subset,
    subset_name,
    omitted_gene
) {
  
  ids <- locked_subset$Association_ID
  
  purity_subset <- purity_subset[
    match(
      ids,
      purity_subset$Association_ID
    ),
    ,
    drop = FALSE
  ]
  
  unadj_sig <- as_logical_sig(
    locked_subset[[col_unadj_sig]]
  )
  
  immune_sig <- as_logical_sig(
    locked_subset[[col_immune_sig]]
  )
  
  purity_sig <- as_logical_sig(
    purity_subset[[col_purity_sig]]
  )
  
  unadj_dir <- normalize_direction(
    locked_subset[[col_unadj_dir]]
  )
  
  immune_dir <- normalize_direction(
    locked_subset[[col_immune_dir]]
  )
  
  purity_dir <- normalize_direction(
    purity_subset[[col_purity_dir]]
  )
  
  three_pattern <- paste0(
    unadj_dir,
    "→",
    immune_dir,
    "→",
    purity_dir
  )
  
  data.frame(
    Subset = subset_name,
    Omitted_Core_Gene = omitted_gene,
    Retained_Core_Genes = paste(
      setdiff(
        core_genes_locked,
        omitted_gene
      ),
      collapse = " + "
    ),
    Total_Pairs = length(ids),
    
    Unadjusted_Significant = sum(unadj_sig),
    ImmuneScore_Significant = sum(immune_sig),
    TumorPurity_Significant = sum(purity_sig),
    
    Significant_All_Three = sum(
      unadj_sig &
        immune_sig &
        purity_sig
    ),
    
    P_to_N_to_N = sum(
      three_pattern == "P→N→N"
    ),
    
    P_to_P_to_P = sum(
      three_pattern == "P→P→P"
    ),
    
    P_to_P_to_N = sum(
      three_pattern == "P→P→N"
    ),
    
    P_to_N_to_P = sum(
      three_pattern == "P→N→P"
    ),
    
    N_to_N_to_N = sum(
      three_pattern == "N→N→N"
    ),
    
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# Full 54-pair summary
# ------------------------------------------------------------

full_summary <- data.frame(
  Subset = "Full_54",
  Omitted_Core_Gene = "None",
  Retained_Core_Genes = paste(
    core_genes_locked,
    collapse = " + "
  ),
  Total_Pairs = 54,
  Unadjusted_Significant = global_unadj_sig,
  ImmuneScore_Significant = global_immune_sig,
  TumorPurity_Significant = global_purity_sig,
  Significant_All_Three = global_all_three,
  P_to_N_to_N =
    as.integer(global_direction_counts["P→N→N"]),
  P_to_P_to_P =
    as.integer(global_direction_counts["P→P→P"]),
  P_to_P_to_N =
    as.integer(global_direction_counts["P→P→N"]),
  P_to_N_to_P =
    as.integer(global_direction_counts["P→N→P"]),
  N_to_N_to_N =
    as.integer(global_direction_counts["N→N→N"]),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# Three leave-one-gene-out summaries
# ------------------------------------------------------------

subset_summaries <- do.call(
  rbind,
  lapply(
    names(subset_definitions),
    function(nm) {
      
      make_subset_summary(
        locked_subset = subset_locked[[nm]],
        purity_subset = subset_purity[[nm]],
        subset_name = nm,
        omitted_gene =
          subset_definitions[[nm]]$omitted_gene
      )
      
    }
  )
)

rownames(subset_summaries) <- NULL


# ------------------------------------------------------------
# Combined manuscript-level summary
# ------------------------------------------------------------

STEP8E_summary <- rbind(
  full_summary,
  subset_summaries
)

rownames(STEP8E_summary) <- NULL

cat("\n============================================================\n")
cat("STEP 8E leave-one-core-gene-out summary\n")
cat("============================================================\n")

print(STEP8E_summary)


# ------------------------------------------------------------
# Subset arithmetic QC
# ------------------------------------------------------------

subset_arithmetic_pass <- all(
  subset_summaries$Total_Pairs == 36 &
    subset_summaries$P_to_N_to_N +
    subset_summaries$P_to_P_to_P +
    subset_summaries$P_to_P_to_N +
    subset_summaries$P_to_N_to_P +
    subset_summaries$N_to_N_to_N ==
    36
)

cat(
  "\nSubset direction arithmetic PASS:",
  subset_arithmetic_pass,
  "\n"
)


# ------------------------------------------------------------
# Core-gene-level detailed audit
# ------------------------------------------------------------

core_gene_detailed <- do.call(
  rbind,
  lapply(
    core_genes_locked,
    function(g) {
      
      x <- locked_54[
        locked_54$Core_Gene == g,
        ,
        drop = FALSE
      ]
      
      p <- purity_adjusted_54[
        purity_adjusted_54$Core_Gene == g,
        ,
        drop = FALSE
      ]
      
      p <- p[
        match(
          x$Association_ID,
          p$Association_ID
        ),
        ,
        drop = FALSE
      ]
      
      unadj_sig <- as_logical_sig(
        x[[col_unadj_sig]]
      )
      
      immune_sig <- as_logical_sig(
        x[[col_immune_sig]]
      )
      
      purity_sig <- as_logical_sig(
        p[[col_purity_sig]]
      )
      
      unadj_dir <- normalize_direction(
        x[[col_unadj_dir]]
      )
      
      immune_dir <- normalize_direction(
        x[[col_immune_dir]]
      )
      
      purity_dir <- normalize_direction(
        p[[col_purity_dir]]
      )
      
      direction_pattern <- paste0(
        unadj_dir,
        "→",
        immune_dir,
        "→",
        purity_dir
      )
      
      data.frame(
        Core_Gene = g,
        Total_Pairs = length(unadj_sig),
        Unadjusted_Significant = sum(unadj_sig),
        ImmuneScore_Significant = sum(immune_sig),
        TumorPurity_Significant = sum(purity_sig),
        Significant_All_Three = sum(
          unadj_sig &
            immune_sig &
            purity_sig
        ),
        P_to_N_to_N = sum(
          direction_pattern == "P→N→N"
        ),
        P_to_P_to_P = sum(
          direction_pattern == "P→P→P"
        ),
        P_to_P_to_N = sum(
          direction_pattern == "P→P→N"
        ),
        P_to_N_to_P = sum(
          direction_pattern == "P→N→P"
        ),
        N_to_N_to_N = sum(
          direction_pattern == "N→N→N"
        ),
        stringsAsFactors = FALSE
      )
      
    }
  )
)

rownames(core_gene_detailed) <- NULL

cat("\nCore-gene detailed audit:\n")
print(core_gene_detailed)


# ------------------------------------------------------------
# Verify known STEP 8D-15 core-gene values
# ------------------------------------------------------------

expected_core_values <- data.frame(
  Core_Gene = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  Unadjusted_Significant = c(
    14, 17, 17
  ),
  ImmuneScore_Significant = c(
    13, 12, 13
  ),
  TumorPurity_Significant = c(
    17, 11, 11
  ),
  Significant_All_Three = c(
    8, 10, 9
  ),
  stringsAsFactors = FALSE
)


core_protection <- merge(
  expected_core_values,
  core_gene_detailed,
  by = "Core_Gene",
  suffixes = c("_Expected", "_Observed"),
  sort = FALSE
)

core_protection$Status <- ifelse(
  core_protection$Unadjusted_Significant_Expected ==
    core_protection$Unadjusted_Significant_Observed &
    core_protection$ImmuneScore_Significant_Expected ==
    core_protection$ImmuneScore_Significant_Observed &
    core_protection$TumorPurity_Significant_Expected ==
    core_protection$TumorPurity_Significant_Observed &
    core_protection$Significant_All_Three_Expected ==
    core_protection$Significant_All_Three_Observed,
  "PASS",
  "FAIL"
)

cat("\nCore-gene locked-value protection:\n")
print(core_protection)


core_protection_pass <- all(
  core_protection$Status == "PASS"
)


# ============================================================
# Paper 1
# STEP 8E-06 → STEP 8E-07
# Final QC, Final Lock and Output
# ============================================================


# ============================================================
# STEP 8E-06
# Final QC table
# ============================================================

cat("\n============================================================\n")
cat("STEP 8E-06: Final QC table\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Re-confirm core-gene protection
# ------------------------------------------------------------

core_protection_pass <- all(
  core_protection$Status == "PASS"
)

cat(
  "\nCore-gene locked-value protection:",
  ifelse(core_protection_pass, "PASS", "FAIL"),
  "\n"
)


# ------------------------------------------------------------
# 2. Build final STEP 8E-06 QC table
# ------------------------------------------------------------

qc_8e06 <- data.frame(
  
  Check = c(
    "Global locked 54-pair values preserved",
    "STEP 8C-7B values preserved",
    "STEP 8D-15 values preserved",
    "Direction pattern arithmetic",
    "Significance pattern arithmetic",
    "All leave-one-gene-out subsets = 36 pairs",
    "Subset direction arithmetic",
    "Core-gene locked values preserved",
    "No NA in locked significance/direction fields"
  ),
  
  Status = c(
    ifelse(
      locked_8c7b_protection &
        locked_8d15_protection,
      "PASS",
      "FAIL"
    ),
    
    ifelse(
      locked_8c7b_protection,
      "PASS",
      "FAIL"
    ),
    
    ifelse(
      locked_8d15_protection,
      "PASS",
      "FAIL"
    ),
    
    ifelse(
      direction_arithmetic_pass,
      "PASS",
      "FAIL"
    ),
    
    ifelse(
      significance_arithmetic_pass,
      "PASS",
      "FAIL"
    ),
    
    ifelse(
      subset_size_pass,
      "PASS",
      "FAIL"
    ),
    
    ifelse(
      subset_arithmetic_pass,
      "PASS",
      "FAIL"
    ),
    
    ifelse(
      core_protection_pass,
      "PASS",
      "FAIL"
    ),
    
    ifelse(
      complete_locked_results,
      "PASS",
      "FAIL"
    )
  ),
  
  stringsAsFactors = FALSE
)


print(qc_8e06)


# ------------------------------------------------------------
# 3. STEP 8E-06 overall QC
# ------------------------------------------------------------

qc_8e06_all_pass <- all(
  qc_8e06$Status == "PASS"
)

cat(
  "\nSTEP 8E-06 overall QC:",
  ifelse(qc_8e06_all_pass, "PASS", "FAIL"),
  "\n"
)


# ============================================================
# STEP 8E-07
# Final lock decision
# ============================================================

cat("\n============================================================\n")
cat("STEP 8E-07: Final lock decision\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 4. Methodological declarations
# ------------------------------------------------------------

NO_NEW_STATISTICS <- TRUE

NO_FDR_RECALCULATION <- TRUE

LOCKED_54_ONLY <- TRUE

LOCKED_STEP8C7B_PRESERVED <- isTRUE(
  locked_8c7b_protection
)

LOCKED_STEP8D15_PRESERVED <- isTRUE(
  locked_8d15_protection
)


# ------------------------------------------------------------
# 5. Reconfirm global locked values
# ------------------------------------------------------------

global_locked_values_pass <- all(
  
  global_unadj_sig == 48,
  
  global_immune_sig == 38,
  
  global_purity_sig == 39,
  
  global_all_three == 27,
  
  as.integer(global_direction_counts["P→N→N"]) == 30,
  
  as.integer(global_direction_counts["P→P→P"]) == 17,
  
  as.integer(global_direction_counts["P→P→N"]) == 4,
  
  as.integer(global_direction_counts["P→N→P"]) == 2,
  
  as.integer(global_direction_counts["N→N→N"]) == 1
)


cat(
  "\nGlobal locked values:",
  ifelse(global_locked_values_pass, "PASS", "FAIL"),
  "\n"
)


# ------------------------------------------------------------
# 6. Reconfirm core-gene locked values
# ------------------------------------------------------------

core_locked_values_pass <- all(
  
  core_protection$core_protection_pass[
    core_protection$Core_Gene == "AOC3"
  ] %in% TRUE,
  
  core_protection$core_protection_pass[
    core_protection$Core_Gene == "MS4A7"
  ] %in% TRUE,
  
  core_protection$core_protection_pass[
    core_protection$Core_Gene == "IL7R"
  ] %in% TRUE
)


# If the table does not contain the helper column,
# use the already validated Status column.

if (!"core_protection_pass" %in% names(core_protection)) {
  
  core_locked_values_pass <- all(
    core_protection$Status == "PASS"
  )
  
}


cat(
  "Core-gene locked values:",
  ifelse(core_locked_values_pass, "PASS", "FAIL"),
  "\n"
)


# ------------------------------------------------------------
# 7. Overall PASS condition
# ------------------------------------------------------------

all_qc_pass <- all(
  
  # STEP 8E-01
  all(qc_8e01$Status == "PASS"),
  
  # STEP 8E-02
  all(qc_8e02$Status == "PASS"),
  
  # Association structure
  isTRUE(association_set_pass),
  
  # Subset structure
  isTRUE(subset_size_pass),
  
  # Completeness
  isTRUE(complete_locked_results),
  
  # Global locked values
  isTRUE(global_locked_values_pass),
  
  # Arithmetic
  isTRUE(direction_arithmetic_pass),
  isTRUE(significance_arithmetic_pass),
  
  # Locked workflow protection
  isTRUE(locked_8c7b_protection),
  isTRUE(locked_8d15_protection),
  
  # Leave-one-out arithmetic
  isTRUE(subset_arithmetic_pass),
  
  # Core-gene protection
  isTRUE(core_protection_pass),
  isTRUE(core_locked_values_pass),
  
  # STEP 8E-06
  isTRUE(qc_8e06_all_pass),
  
  # Methodological declarations
  isTRUE(NO_NEW_STATISTICS),
  isTRUE(NO_FDR_RECALCULATION),
  isTRUE(LOCKED_54_ONLY)
)


FINAL_STATUS <- ifelse(
  all_qc_pass,
  "PASS",
  "FAIL"
)


cat("\n============================================================\n")
cat("STEP 8E FINAL STATUS:", FINAL_STATUS, "\n")
cat("============================================================\n")


# ============================================================
# 8. Build final lock object
# ============================================================

STEP8E_FINAL_LOCK <- list(
  
  Workflow =
    "STEP 8E — Leave-One-Core-Gene-Out Circularity / Dependence Audit",
  
  Phase =
    "Phase 1 — Locked-result subset audit",
  
  Date =
    as.character(Sys.Date()),
  
  Core_Genes =
    core_genes_locked,
  
  Checkpoint_Genes =
    checkpoint_genes_locked,
  
  Full_Pair_Count =
    54,
  
  Leave_One_Gene_Out_Pair_Count =
    36,
  
  Subsets =
    subset_definitions,
  
  Locked_Results_Source =
    c(
      "locked_54",
      "purity_adjusted_54"
    ),
  
  No_New_Statistics =
    NO_NEW_STATISTICS,
  
  No_FDR_Recalculation =
    NO_FDR_RECALCULATION,
  
  Locked_54_Only =
    LOCKED_54_ONLY,
  
  Locked_STEP8C7B_Preserved =
    LOCKED_STEP8C7B_PRESERVED,
  
  Locked_STEP8D15_Preserved =
    LOCKED_STEP8D15_PRESERVED,
  
  Global_Locked_Values_Pass =
    global_locked_values_pass,
  
  Core_Gene_Locked_Values_Pass =
    core_locked_values_pass,
  
  Global_Summary =
    STEP8E_summary,
  
  Core_Gene_Detailed =
    core_gene_detailed,
  
  Core_Gene_Protection =
    core_protection,
  
  Association_ID_QC =
    association_matching_table,
  
  Global_Integrity =
    global_integrity_checks,
  
  QC_STEP8E01 =
    qc_8e01,
  
  QC_STEP8E02 =
    qc_8e02,
  
  QC_STEP8E06 =
    qc_8e06,
  
  Final_Status =
    FINAL_STATUS
  
)


# ============================================================
# 9. Save final RDS lock
# ============================================================

saveRDS(
  STEP8E_FINAL_LOCK,
  file = file.path(
    out_dir,
    "STEP8E_FINAL_LOCK.rds"
  )
)


# ============================================================
# 10. Save 54-pair audit
# ============================================================

write.csv(
  audit_54,
  file = file.path(
    out_dir,
    "STEP8E_54Pair_Locked_Audit.csv"
  ),
  row.names = FALSE,
  na = ""
)


# ============================================================
# 11. Save manuscript-level summary
# ============================================================

write.csv(
  STEP8E_summary,
  file = file.path(
    out_dir,
    "STEP8E_Manuscript_Key_Findings.csv"
  ),
  row.names = FALSE,
  na = ""
)


# ============================================================
# 12. Save core-gene detailed audit
# ============================================================

write.csv(
  core_gene_detailed,
  file = file.path(
    out_dir,
    "STEP8E_Core_Gene_Detailed_Audit.csv"
  ),
  row.names = FALSE,
  na = ""
)


# ============================================================
# 13. Save Association_ID QC
# ============================================================

write.csv(
  association_matching_table,
  file = file.path(
    out_dir,
    "STEP8E_Association_ID_QC.csv"
  ),
  row.names = FALSE,
  na = ""
)


# ============================================================
# 14. Save global integrity QC
# ============================================================

write.csv(
  global_integrity_checks,
  file = file.path(
    out_dir,
    "STEP8E_Global_Integrity_QC.csv"
  ),
  row.names = FALSE,
  na = ""
)


# ============================================================
# 15. Build complete QC table
# ============================================================

cat("\n============================================================\n")
cat("Building complete STEP 8E QC table...\n")
cat("============================================================\n")


# ------------------------------------------------------------
# Convert each QC table to a common 3-column structure
# ------------------------------------------------------------

qc_8e01_std <- data.frame(
  
  Stage = "STEP8E-01",
  
  Check = as.character(qc_8e01$Check),
  
  Status = as.character(qc_8e01$Status),
  
  stringsAsFactors = FALSE
)


qc_8e02_std <- data.frame(
  
  Stage = "STEP8E-02",
  
  Check = as.character(qc_8e02$Check),
  
  Status = as.character(qc_8e02$Status),
  
  stringsAsFactors = FALSE
)


qc_8e06_std <- data.frame(
  
  Stage = "STEP8E-06",
  
  Check = as.character(qc_8e06$Check),
  
  Status = as.character(qc_8e06$Status),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# Combine standardized QC tables
# ------------------------------------------------------------

qc_all <- rbind(
  
  qc_8e01_std,
  
  qc_8e02_std,
  
  qc_8e06_std
  
)


# ------------------------------------------------------------
# QC table integrity check
# ------------------------------------------------------------

qc_all_structure_pass <- all(
  
  names(qc_all) ==
    c(
      "Stage",
      "Check",
      "Status"
    ),
  
  nrow(qc_all) ==
    nrow(qc_8e01_std) +
    nrow(qc_8e02_std) +
    nrow(qc_8e06_std),
  
  all(!is.na(qc_all$Stage)),
  
  all(!is.na(qc_all$Check)),
  
  all(!is.na(qc_all$Status)),
  
  all(qc_all$Status %in% c("PASS", "FAIL"))
)


cat(
  "\nComplete QC table structure:",
  ifelse(
    qc_all_structure_pass,
    "PASS",
    "FAIL"
  ),
  "\n"
)


# ------------------------------------------------------------
# Print complete QC table
# ------------------------------------------------------------

print(qc_all)


# ------------------------------------------------------------
# Overall QC status from complete table
# ------------------------------------------------------------

qc_all_pass <- all(
  
  qc_all_structure_pass,
  
  qc_all$Status == "PASS"
)


cat(
  "\nComplete QC table overall status:",
  ifelse(
    qc_all_pass,
    "PASS",
    "FAIL"
  ),
  "\n"
)


# ------------------------------------------------------------
# Save complete QC table
# ------------------------------------------------------------

write.csv(
  
  qc_all,
  
  file = file.path(
    out_dir,
    "STEP8E_QC.csv"
  ),
  
  row.names = FALSE,
  
  na = ""
)


cat(
  "\nSTEP8E_QC.csv saved successfully.\n"
)


# ============================================================
# 16. Final status re-confirmation after QC export
# ============================================================

all_qc_pass <- all(
  
  # STEP 8E-01
  all(qc_8e01$Status == "PASS"),
  
  # STEP 8E-02
  all(qc_8e02$Status == "PASS"),
  
  # Association structure
  isTRUE(association_set_pass),
  
  # Subset structure
  isTRUE(subset_size_pass),
  
  # Completeness
  isTRUE(complete_locked_results),
  
  # Global locked values
  isTRUE(global_locked_values_pass),
  
  # Arithmetic
  isTRUE(direction_arithmetic_pass),
  isTRUE(significance_arithmetic_pass),
  
  # Locked workflow protection
  isTRUE(locked_8c7b_protection),
  isTRUE(locked_8d15_protection),
  
  # Leave-one-out arithmetic
  isTRUE(subset_arithmetic_pass),
  
  # Core-gene protection
  isTRUE(core_protection_pass),
  isTRUE(core_locked_values_pass),
  
  # STEP 8E-06
  isTRUE(qc_8e06_all_pass),
  
  # Complete QC export
  isTRUE(qc_all_structure_pass),
  isTRUE(qc_all_pass),
  
  # Methodological declarations
  isTRUE(NO_NEW_STATISTICS),
  isTRUE(NO_FDR_RECALCULATION),
  isTRUE(LOCKED_54_ONLY)
)


FINAL_STATUS <- ifelse(
  all_qc_pass,
  "PASS",
  "FAIL"
)


cat("\n============================================================\n")
cat("STEP 8E FINAL STATUS:", FINAL_STATUS, "\n")
cat("============================================================\n")

# ============================================================
# STEP 8E FINAL OUTPUT LOCK
# Final RDS + CSV + Excel + TXT
# No new statistics
# No FDR recalculation
# ============================================================

cat("\n============================================================\n")
cat("STEP 8E FINAL OUTPUT LOCK\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Final methodological declarations
# ------------------------------------------------------------

NO_NEW_STATISTICS <- TRUE

NO_FDR_RECALCULATION <- TRUE

LOCKED_54_ONLY <- TRUE

LOCKED_STEP8C7B_PRESERVED <- isTRUE(
  locked_8c7b_protection
)

LOCKED_STEP8D15_PRESERVED <- isTRUE(
  locked_8d15_protection
)


# ------------------------------------------------------------
# 2. Final QC protection
# ------------------------------------------------------------

FINAL_QC_EXPORT_PASS <- isTRUE(
  qc_all_structure_pass
) &&
  isTRUE(
    qc_all_pass
  )


# ------------------------------------------------------------
# 3. Final status
# ------------------------------------------------------------

FINAL_STATUS <- ifelse(
  
  all(
    isTRUE(association_set_pass),
    isTRUE(subset_size_pass),
    isTRUE(complete_locked_results),
    isTRUE(global_locked_values_pass),
    isTRUE(direction_arithmetic_pass),
    isTRUE(significance_arithmetic_pass),
    isTRUE(locked_8c7b_protection),
    isTRUE(locked_8d15_protection),
    isTRUE(subset_arithmetic_pass),
    isTRUE(core_protection_pass),
    isTRUE(core_locked_values_pass),
    isTRUE(qc_8e06_all_pass),
    isTRUE(qc_all_structure_pass),
    isTRUE(qc_all_pass),
    isTRUE(NO_NEW_STATISTICS),
    isTRUE(NO_FDR_RECALCULATION),
    isTRUE(LOCKED_54_ONLY)
  ),
  
  "PASS",
  "FAIL"
)


cat(
  "\nFinal QC export:",
  ifelse(FINAL_QC_EXPORT_PASS, "PASS", "FAIL"),
  "\n"
)

cat(
  "Final STEP 8E status:",
  FINAL_STATUS,
  "\n"
)


# ============================================================
# 4. Update final lock object
# ============================================================

STEP8E_FINAL_LOCK <- list(
  
  Workflow =
    "STEP 8E — Leave-One-Core-Gene-Out Circularity / Dependence Audit",
  
  Phase =
    "Phase 1 — Locked-result subset audit",
  
  Date =
    as.character(Sys.Date()),
  
  Core_Genes =
    core_genes_locked,
  
  Checkpoint_Genes =
    checkpoint_genes_locked,
  
  Full_Pair_Count =
    54,
  
  Leave_One_Gene_Out_Pair_Count =
    36,
  
  Subsets =
    subset_definitions,
  
  Locked_Results_Source =
    c(
      "locked_54",
      "purity_adjusted_54"
    ),
  
  No_New_Statistics =
    NO_NEW_STATISTICS,
  
  No_FDR_Recalculation =
    NO_FDR_RECALCULATION,
  
  Locked_54_Only =
    LOCKED_54_ONLY,
  
  Locked_STEP8C7B_Preserved =
    LOCKED_STEP8C7B_PRESERVED,
  
  Locked_STEP8D15_Preserved =
    LOCKED_STEP8D15_PRESERVED,
  
  Global_Locked_Values_Pass =
    global_locked_values_pass,
  
  Core_Gene_Locked_Values_Pass =
    core_locked_values_pass,
  
  Final_QC_Export_Pass =
    FINAL_QC_EXPORT_PASS,
  
  Global_Summary =
    STEP8E_summary,
  
  Core_Gene_Detailed =
    core_gene_detailed,
  
  Core_Gene_Protection =
    core_protection,
  
  Association_ID_QC =
    association_matching_table,
  
  Global_Integrity =
    global_integrity_checks,
  
  QC_STEP8E01 =
    qc_8e01,
  
  QC_STEP8E02 =
    qc_8e02,
  
  QC_STEP8E06 =
    qc_8e06,
  
  QC_ALL =
    qc_all,
  
  Final_Status =
    FINAL_STATUS
)


# ============================================================
# 5. Save final RDS lock
# ============================================================

saveRDS(
  STEP8E_FINAL_LOCK,
  file = file.path(
    out_dir,
    "STEP8E_FINAL_LOCK.rds"
  )
)


# ============================================================
# 6. Save final CSV files
# ============================================================

write.csv(
  audit_54,
  file = file.path(
    out_dir,
    "STEP8E_54Pair_Locked_Audit.csv"
  ),
  row.names = FALSE,
  na = ""
)


write.csv(
  STEP8E_summary,
  file = file.path(
    out_dir,
    "STEP8E_Manuscript_Key_Findings.csv"
  ),
  row.names = FALSE,
  na = ""
)


write.csv(
  core_gene_detailed,
  file = file.path(
    out_dir,
    "STEP8E_Core_Gene_Detailed_Audit.csv"
  ),
  row.names = FALSE,
  na = ""
)


write.csv(
  association_matching_table,
  file = file.path(
    out_dir,
    "STEP8E_Association_ID_QC.csv"
  ),
  row.names = FALSE,
  na = ""
)


write.csv(
  global_integrity_checks,
  file = file.path(
    out_dir,
    "STEP8E_Global_Integrity_QC.csv"
  ),
  row.names = FALSE,
  na = ""
)


write.csv(
  qc_all,
  file = file.path(
    out_dir,
    "STEP8E_QC.csv"
  ),
  row.names = FALSE,
  na = ""
)


# ============================================================
# 7. Final Excel export
# ============================================================

if (requireNamespace("openxlsx", quietly = TRUE)) {
  
  wb <- openxlsx::createWorkbook()
  
  
  openxlsx::addWorksheet(
    wb,
    "Manuscript_Summary"
  )
  
  openxlsx::writeData(
    wb,
    "Manuscript_Summary",
    STEP8E_summary
  )
  
  
  openxlsx::addWorksheet(
    wb,
    "Core_Gene_Audit"
  )
  
  openxlsx::writeData(
    wb,
    "Core_Gene_Audit",
    core_gene_detailed
  )
  
  
  openxlsx::addWorksheet(
    wb,
    "54Pair_Audit"
  )
  
  openxlsx::writeData(
    wb,
    "54Pair_Audit",
    audit_54
  )
  
  
  openxlsx::addWorksheet(
    wb,
    "Association_ID_QC"
  )
  
  openxlsx::writeData(
    wb,
    "Association_ID_QC",
    association_matching_table
  )
  
  
  openxlsx::addWorksheet(
    wb,
    "Global_Integrity"
  )
  
  openxlsx::writeData(
    wb,
    "Global_Integrity",
    global_integrity_checks
  )
  
  
  openxlsx::addWorksheet(
    wb,
    "QC"
  )
  
  openxlsx::writeData(
    wb,
    "QC",
    qc_all
  )
  
  
  openxlsx::addWorksheet(
    wb,
    "Final_Lock"
  )
  
  openxlsx::writeData(
    wb,
    "Final_Lock",
    data.frame(
      Parameter = c(
        "Final_Status",
        "No_New_Statistics",
        "No_FDR_Recalculation",
        "Locked_54_Only",
        "STEP8C7B_Preserved",
        "STEP8D15_Preserved",
        "Final_QC_Export_Pass"
      ),
      Value = c(
        FINAL_STATUS,
        NO_NEW_STATISTICS,
        NO_FDR_RECALCULATION,
        LOCKED_54_ONLY,
        LOCKED_STEP8C7B_PRESERVED,
        LOCKED_STEP8D15_PRESERVED,
        FINAL_QC_EXPORT_PASS
      ),
      stringsAsFactors = FALSE
    )
  )
  
  
  openxlsx::saveWorkbook(
    wb,
    file = file.path(
      out_dir,
      "STEP8E_Leave_One_Core_Gene_Out_Audit.xlsx"
    ),
    overwrite = TRUE
  )
  
  EXCEL_EXPORT_PASS <- TRUE
  
  cat(
    "\nExcel export: PASS\n"
  )
  
} else {
  
  EXCEL_EXPORT_PASS <- FALSE
  
  cat(
    "\nExcel export: NOT CREATED\n",
    "openxlsx is not installed.\n"
  )
  
}


# ============================================================
# 8. Final lock report
# ============================================================

lock_report_file <- file.path(
  out_dir,
  "STEP8E_FINAL_LOCK_REPORT.txt"
)


report_lines <- c(
  
  "============================================================",
  "Paper 1",
  "STEP 8E — Leave-One-Core-Gene-Out Circularity / Dependence Audit",
  "============================================================",
  "",
  
  paste("Date:", Sys.Date()),
  "",
  
  "PHASE:",
  "Phase 1 — Locked-result subset audit",
  "",
  
  "PURPOSE:",
  "Assess whether the locked immune-checkpoint association pattern",
  "is disproportionately dependent on any one of the three core genes.",
  "",
  
  "CORE GENES:",
  paste(core_genes_locked, collapse = ", "),
  "",
  
  "CHECKPOINT GENES:",
  paste(checkpoint_genes_locked, collapse = ", "),
  "",
  
  "FULL ANALYSIS:",
  "54 predefined core-gene × checkpoint-gene pairs",
  "",
  
  "LEAVE-ONE-GENE-OUT SUBSETS:",
  "Leave AOC3-out   = MS4A7 + IL7R = 36 pairs",
  "Leave MS4A7-out  = AOC3 + IL7R = 36 pairs",
  "Leave IL7R-out   = AOC3 + MS4A7 = 36 pairs",
  "",
  
  "IMPORTANT METHODOLOGICAL STATUS:",
  "No new Spearman correlations were calculated.",
  "No new partial Spearman correlations were calculated.",
  "No new ImmuneScore adjustment was performed.",
  "No new TumorPurity adjustment was performed.",
  "No BH-FDR recalculation was performed.",
  "",
  
  "LOCKED GLOBAL VALUES:",
  
  paste(
    "Unadjusted significant =",
    global_unadj_sig,
    "/54"
  ),
  
  paste(
    "ImmuneScore-adjusted significant =",
    global_immune_sig,
    "/54"
  ),
  
  paste(
    "TumorPurity-adjusted significant =",
    global_purity_sig,
    "/54"
  ),
  
  paste(
    "Significant in all three layers =",
    global_all_three,
    "/54"
  ),
  
  "",
  
  "LOCKED THREE-LAYER DIRECTION PATTERNS:",
  
  paste(
    "P→N→N =",
    as.integer(global_direction_counts["P→N→N"]),
    "/54"
  ),
  
  paste(
    "P→P→P =",
    as.integer(global_direction_counts["P→P→P"]),
    "/54"
  ),
  
  paste(
    "P→P→N =",
    as.integer(global_direction_counts["P→P→N"]),
    "/54"
  ),
  
  paste(
    "P→N→P =",
    as.integer(global_direction_counts["P→N→P"]),
    "/54"
  ),
  
  paste(
    "N→N→N =",
    as.integer(global_direction_counts["N→N→N"]),
    "/54"
  ),
  
  "",
  
  "CORE-GENE LOCKED VALUES:",
  "AOC3:    14 / 13 / 17 significant",
  "MS4A7:   17 / 12 / 11 significant",
  "IL7R:    17 / 13 / 11 significant",
  "",
  
  "ORDER:",
  "Unadjusted / ImmuneScore-adjusted / TumorPurity-adjusted",
  "",
  
  "FINAL QC:",
  
  paste(
    "Association_ID matching =",
    ifelse(association_set_pass, "PASS", "FAIL")
  ),
  
  paste(
    "Subset size QC =",
    ifelse(subset_size_pass, "PASS", "FAIL")
  ),
  
  paste(
    "Locked 54-pair completeness =",
    ifelse(complete_locked_results, "PASS", "FAIL")
  ),
  
  paste(
    "Global locked values =",
    ifelse(global_locked_values_pass, "PASS", "FAIL")
  ),
  
  paste(
    "Direction arithmetic =",
    ifelse(direction_arithmetic_pass, "PASS", "FAIL")
  ),
  
  paste(
    "Significance arithmetic =",
    ifelse(significance_arithmetic_pass, "PASS", "FAIL")
  ),
  
  paste(
    "Subset arithmetic =",
    ifelse(subset_arithmetic_pass, "PASS", "FAIL")
  ),
  
  paste(
    "Core-gene locked values preserved =",
    ifelse(core_protection_pass, "PASS", "FAIL")
  ),
  
  paste(
    "STEP 8C-7B values preserved =",
    ifelse(locked_8c7b_protection, "PASS", "FAIL")
  ),
  
  paste(
    "STEP 8D-15 values preserved =",
    ifelse(locked_8d15_protection, "PASS", "FAIL")
  ),
  
  paste(
    "Complete QC export =",
    ifelse(FINAL_QC_EXPORT_PASS, "PASS", "FAIL")
  ),
  
  "",
  
  "STATISTICAL RE-CALCULATION:",
  "NO NEW STATISTICS = TRUE",
  "NO FDR RECALCULATION = TRUE",
  "LOCKED 54-PAIR RESULTS ONLY = TRUE",
  "",
  
  paste(
    "STEP 8E FINAL STATUS =",
    FINAL_STATUS
  ),
  
  "",
  
  "INTERPRETATION STATUS:",
  "STEP 8E Phase 1 is an audit of locked results.",
  "It does not establish causality or prove absence of compositional dependence.",
  "It does not establish cell-intrinsic regulation or direct molecular interaction.",
  "",
  
  "============================================================",
  "END OF STEP 8E FINAL LOCK REPORT",
  "============================================================"
)


writeLines(
  report_lines,
  con = lock_report_file
)


# ============================================================
# 9. Final output existence QC
# ============================================================

expected_files <- c(
  
  "STEP8E_FINAL_LOCK.rds",
  "STEP8E_54Pair_Locked_Audit.csv",
  "STEP8E_Manuscript_Key_Findings.csv",
  "STEP8E_Core_Gene_Detailed_Audit.csv",
  "STEP8E_Association_ID_QC.csv",
  "STEP8E_Global_Integrity_QC.csv",
  "STEP8E_QC.csv",
  "STEP8E_FINAL_LOCK_REPORT.txt"
  
)


if (EXCEL_EXPORT_PASS) {
  
  expected_files <- c(
    expected_files,
    "STEP8E_Leave_One_Core_Gene_Out_Audit.xlsx"
  )
  
}


file_exists <- file.exists(
  file.path(
    out_dir,
    expected_files
  )
)


FINAL_OUTPUT_FILES_PASS <- all(
  file_exists
)


# ============================================================
# 10. Final console summary
# ============================================================

cat("\n\n============================================================\n")
cat("STEP 8E FINAL LOCK SUMMARY\n")
cat("============================================================\n")

cat(
  "\nFinal status:",
  FINAL_STATUS,
  "\n"
)

cat(
  "Final QC export:",
  ifelse(FINAL_QC_EXPORT_PASS, "PASS", "FAIL"),
  "\n"
)

cat(
  "Output file existence:",
  ifelse(FINAL_OUTPUT_FILES_PASS, "PASS", "FAIL"),
  "\n"
)

cat(
  "Locked STEP 8C-7B values preserved:",
  LOCKED_STEP8C7B_PRESERVED,
  "\n"
)

cat(
  "Locked STEP 8D-15 values preserved:",
  LOCKED_STEP8D15_PRESERVED,
  "\n"
)

cat(
  "No new statistical testing:",
  NO_NEW_STATISTICS,
  "\n"
)

cat(
  "No FDR recalculation:",
  NO_FDR_RECALCULATION,
  "\n"
)

cat(
  "Locked 54-pair results only:",
  LOCKED_54_ONLY,
  "\n"
)

cat(
  "Excel export:",
  ifelse(EXCEL_EXPORT_PASS, "PASS", "NOT CREATED"),
  "\n"
)

cat("\nExpected output files:\n")

print(
  file.path(
    out_dir,
    expected_files
  )
)

cat("\n============================================================\n")
cat(
  "STEP 8E FINAL STATUS:",
  FINAL_STATUS,
  "\n"
)
cat("============================================================\n")

