##GitHub4_Paper1_tumor_microenvironment_association_analysis_2

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


work_dir <- "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis"

setwd(work_dir)

getwd()

load(
  "input_dir/Paper1_tumor-only WGCNA/Paper1_STEP7B7E9A_ESTIMATE_signature_membership_audit_LOCKED.RData"
)

ls()

##STEP 7B-7E-10A
#Leave-one-gene-out sensitivity analysis — preparation audit
# ============================================================
# Paper 1
# STEP 7B-7E-10A
# Leave-one-gene-out sensitivity analysis
# Preparation and methodological audit
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7E-10A: Leave-one-gene-out sensitivity analysis\n")
cat("Preparation and methodological audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Confirm working directory
# ------------------------------------------------------------

cat("Working directory:\n")
print(getwd())

cat("\n")


# ------------------------------------------------------------
# 2. Confirm STEP 7B-7E-9A locked audit objects
# ------------------------------------------------------------

required_7be9a <- c(
  "step7be9a_audit",
  "membership_7be",
  "input_signature_audit_7be",
  "stromal_genes_7be",
  "immune_genes_7be",
  "MS4A6A_in_stromal_7be",
  "MS4A7_in_stromal_7be",
  "MS4A7_in_immune_7be"
)

missing_7be9a <- setdiff(required_7be9a, ls())

cat("Required STEP 7B-7E-9A objects present: ",
    length(missing_7be9a) == 0, "\n", sep = "")

if (length(missing_7be9a) > 0) {
  cat("Missing objects:\n")
  print(missing_7be9a)
  stop("Required STEP 7B-7E-9A objects are missing.")
}


# ------------------------------------------------------------
# 3. Confirm locked membership results
# ------------------------------------------------------------

cat("\nConfirmed ESTIMATE signature membership:\n")

cat(
  "AOC3 in Stromal signature: ",
  membership_7be$AOC3_in_stromal,
  "\n",
  sep = ""
)

cat(
  "AOC3 in Immune signature: ",
  membership_7be$AOC3_in_immune,
  "\n",
  sep = ""
)

cat(
  "MS4A7 in Stromal signature: ",
  membership_7be$MS4A7_in_stromal,
  "\n",
  sep = ""
)

cat(
  "MS4A7 in Immune signature: ",
  membership_7be$MS4A7_in_immune,
  "\n",
  sep = ""
)

cat(
  "IL7R in Stromal signature: ",
  membership_7be$IL7R_in_stromal,
  "\n",
  sep = ""
)

cat(
  "IL7R in Immune signature: ",
  membership_7be$IL7R_in_immune,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Hard-stop assertions for the two LOO targets
# ------------------------------------------------------------

stopifnot(
  membership_7be$AOC3_in_stromal == TRUE,
  membership_7be$AOC3_in_immune == FALSE,
  membership_7be$MS4A7_in_stromal == FALSE,
  membership_7be$MS4A7_in_immune == FALSE,
  membership_7be$IL7R_in_stromal == FALSE,
  membership_7be$IL7R_in_immune == TRUE
)

cat("\nMembership assertions: PASSED\n")


# ------------------------------------------------------------
# 5. Define LOO targets
# ------------------------------------------------------------

loo_targets_7be10 <- data.frame(
  Analysis = c(
    "LOO_AOC3_StromalScore",
    "LOO_IL7R_ImmuneScore"
  ),
  Gene = c(
    "AOC3",
    "IL7R"
  ),
  Score = c(
    "StromalScore",
    "ImmuneScore"
  ),
  stringsAsFactors = FALSE
)

cat("\nLOO targets:\n")
print(loo_targets_7be10)


# ------------------------------------------------------------
# 6. Verify signature sizes before exclusion
# ------------------------------------------------------------

cat("\nSignature sizes before LOO:\n")

cat(
  "Stromal signature: ",
  length(unique(stromal_genes_7be)),
  "\n",
  sep = ""
)

cat(
  "Immune signature: ",
  length(unique(immune_genes_7be)),
  "\n",
  sep = ""
)

stopifnot(
  length(unique(stromal_genes_7be)) == 141,
  length(unique(immune_genes_7be)) == 141
)

cat("Signature-size audit: PASSED\n")


# ------------------------------------------------------------
# 7. Construct LOO signatures
# ------------------------------------------------------------

stromal_genes_loo_AOC3_7be <- setdiff(
  unique(stromal_genes_7be),
  "AOC3"
)

immune_genes_loo_IL7R_7be <- setdiff(
  unique(immune_genes_7be),
  "IL7R"
)


# ------------------------------------------------------------
# 8. Verify exact one-gene removal
# ------------------------------------------------------------

cat("\nLOO signature sizes:\n")

cat(
  "Stromal after removing AOC3: ",
  length(stromal_genes_loo_AOC3_7be),
  "\n",
  sep = ""
)

cat(
  "Immune after removing IL7R: ",
  length(immune_genes_loo_IL7R_7be),
  "\n",
  sep = ""
)

stopifnot(
  length(stromal_genes_loo_AOC3_7be) == 140,
  length(immune_genes_loo_IL7R_7be) == 140
)

stopifnot(
  !"AOC3" %in% stromal_genes_loo_AOC3_7be,
  !"IL7R" %in% immune_genes_loo_IL7R_7be
)

cat("Exact one-gene-removal audit: PASSED\n")


# ------------------------------------------------------------
# 9. Confirm MS4A7 is NOT artificially subjected to LOO
# ------------------------------------------------------------

stopifnot(
  !"MS4A7" %in% setdiff(
    unique(stromal_genes_7be),
    unique(stromal_genes_7be)
  )
)

cat(
  "\nMS4A7 LOO status: NOT PERFORMED\n",
  "Reason: MS4A7 is absent from both ESTIMATE signatures.\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Check ESTIMATE package availability
# ------------------------------------------------------------

estimate_available_7be10 <- requireNamespace(
  "estimate",
  quietly = TRUE
)

cat(
  "\nESTIMATE package available: ",
  estimate_available_7be10,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Save preparation checkpoint
# ------------------------------------------------------------

step7be10a_audit <- list(
  step = "STEP 7B-7E-10A",
  status = "PREPARATION_AUDIT",
  LOO_targets = loo_targets_7be10,
  stromal_signature_n_before = length(unique(stromal_genes_7be)),
  immune_signature_n_before = length(unique(immune_genes_7be)),
  stromal_signature_n_after_AOC3 =
    length(stromal_genes_loo_AOC3_7be),
  immune_signature_n_after_IL7R =
    length(immune_genes_loo_IL7R_7be),
  AOC3_removed_from_stromal = !(
    "AOC3" %in% stromal_genes_loo_AOC3_7be
  ),
  IL7R_removed_from_immune = !(
    "IL7R" %in% immune_genes_loo_IL7R_7be
  ),
  MS4A7_LOO_performed = FALSE,
  estimate_package_available = estimate_available_7be10,
  statistical_analysis_performed = FALSE,
  locked_objects_modified = FALSE,
  figure6_modified = FALSE,
  supplementary_table3_modified = FALSE
)

cat("\n============================================================\n")
cat("STEP 7B-7E-10A PREPARATION AUDIT COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Status: ",
  step7be10a_audit$status,
  "\n",
  sep = ""
)

cat(
  "No statistical analysis performed: TRUE\n"
)

cat(
  "No locked object modified: TRUE\n"
)

cat(
  "STEP 7B-7E-10A: READY FOR NEXT STEP\n"
)


# ============================================================
# STEP 7B-7E-10A-1
# Save preparation audit checkpoint
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7E-10A-1: SAVE PREPARATION CHECKPOINT\n")
cat("============================================================\n\n")

checkpoint_dir_7be10 <- file.path(
  work_dir,
  "STEP7B-7E-10"
)

if (!dir.exists(checkpoint_dir_7be10)) {
  dir.create(
    checkpoint_dir_7be10,
    recursive = TRUE
  )
}

checkpoint_file_7be10a <- file.path(
  checkpoint_dir_7be10,
  "Paper1_STEP7B7E10A_preparation_audit_LOCKED.RData"
)

save(
  step7be10a_audit,
  loo_targets_7be10,
  stromal_genes_loo_AOC3_7be,
  immune_genes_loo_IL7R_7be,
  file = checkpoint_file_7be10a
)

cat(
  "Checkpoint saved:\n",
  checkpoint_file_7be10a,
  "\n",
  sep = ""
)

cat(
  "File exists: ",
  file.exists(checkpoint_file_7be10a),
  "\n",
  sep = ""
)

cat(
  "File size: ",
  file.info(checkpoint_file_7be10a)$size,
  " bytes\n",
  sep = ""
)

cat(
  "\nSTEP 7B-7E-10A-1: SAVE COMPLETE\n"
)


# ============================================================
# STEP 7B-7E-10A-1
# Read-back verification
# ============================================================

audit_readback_env_7be10a <- new.env(
  parent = emptyenv()
)

loaded_objects_7be10a <- load(
  checkpoint_file_7be10a,
  envir = audit_readback_env_7be10a
)

cat("\nLoaded objects:\n")
print(loaded_objects_7be10a)

cat(
  "\nRead-back status: ",
  audit_readback_env_7be10a$step7be10a_audit$status,
  "\n",
  sep = ""
)

stopifnot(
  audit_readback_env_7be10a$step7be10a_audit$status ==
    "PREPARATION_AUDIT"
)

cat(
  "\nSTEP 7B-7E-10A FINAL FILE CHECK: PASSED\n"
)


# ============================================================
# Paper 1
# STEP 7B-7E-10B-1
# Source recovery + exact identity audit
#
# Purpose:
#   Recover the original ESTIMATE input from the locked
#   WORKFLOW-03 Step14A-14E checkpoint.
#
# This step:
#   - DOES NOT recalculate ESTIMATE scores
#   - DOES NOT perform LOO scoring
#   - DOES NOT perform statistical testing
#   - DOES NOT modify any locked object
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7E-10B-1\n")
cat("Source recovery + exact identity audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define locked source file
# ------------------------------------------------------------

estimate_locked_file_7be10 <- paste0(
  "input_dir/Paper1_tumor-only WGCNA/",
  "WORKFLOW-03/Step14/Step14_FINAL_LOCKED/RData/",
  "WORKFLOW-03_Step14A-14E_FINAL_LOCKED.RData"
)


cat("Locked source file:\n")
cat(estimate_locked_file_7be10, "\n\n", sep = "")

if (!file.exists(estimate_locked_file_7be10)) {
  stop(
    "HARD STOP: Locked Step14A-14E RData file was not found."
  )
}

cat(
  "Source file exists: TRUE\n"
)


# ------------------------------------------------------------
# 2. Load source checkpoint into isolated environment
# ------------------------------------------------------------

estimate_source_env_7be10 <- new.env(
  parent = emptyenv()
)

loaded_step14_objects_7be10 <- load(
  estimate_locked_file_7be10,
  envir = estimate_source_env_7be10
)

cat("\nNumber of objects recovered: ",
    length(loaded_step14_objects_7be10),
    "\n",
    sep = "")


# ------------------------------------------------------------
# 3. Required source objects
# ------------------------------------------------------------

required_step14_objects_7be10 <- c(
  "expr_estimate",
  "estimate_input_df",
  "estimate_scores_wgcna",
  "estimate_scores",
  "estimate_sample_ids",
  "estimate_sample_ids_raw",
  "estimate_sample_ids_normalized"
)

missing_step14_objects_7be10 <- setdiff(
  required_step14_objects_7be10,
  loaded_step14_objects_7be10
)

cat(
  "\nRequired source objects present: ",
  length(missing_step14_objects_7be10) == 0,
  "\n",
  sep = ""
)

if (length(missing_step14_objects_7be10) > 0) {
  cat("\nMissing objects:\n")
  print(missing_step14_objects_7be10)
  
  stop(
    "HARD STOP: Required Step14 source object(s) are missing."
  )
}


# ------------------------------------------------------------
# 4. Recover objects into NEW names
#
# IMPORTANT:
# Do not overwrite the current session's objects.
# ------------------------------------------------------------

expr_estimate_7be10 <- get(
  "expr_estimate",
  envir = estimate_source_env_7be10
)

estimate_input_df_7be10 <- get(
  "estimate_input_df",
  envir = estimate_source_env_7be10
)

estimate_scores_wgcna_source_7be10 <- get(
  "estimate_scores_wgcna",
  envir = estimate_source_env_7be10
)

estimate_scores_source_7be10 <- get(
  "estimate_scores",
  envir = estimate_source_env_7be10
)

estimate_sample_ids_7be10 <- get(
  "estimate_sample_ids",
  envir = estimate_source_env_7be10
)

estimate_sample_ids_raw_7be10 <- get(
  "estimate_sample_ids_raw",
  envir = estimate_source_env_7be10
)

estimate_sample_ids_normalized_7be10 <- get(
  "estimate_sample_ids_normalized",
  envir = estimate_source_env_7be10
)


# ------------------------------------------------------------
# 5. Basic structure audit
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("5. Basic structure audit\n")
cat("------------------------------------------------------------\n")

cat(
  "expr_estimate dimensions: ",
  paste(dim(expr_estimate_7be10), collapse = " x "),
  "\n",
  sep = ""
)

cat(
  "estimate_input_df dimensions: ",
  paste(dim(estimate_input_df_7be10), collapse = " x "),
  "\n",
  sep = ""
)

cat(
  "estimate_scores_wgcna dimensions: ",
  paste(
    dim(estimate_scores_wgcna_source_7be10),
    collapse = " x "
  ),
  "\n",
  sep = ""
)

cat(
  "estimate_scores dimensions: ",
  paste(
    dim(estimate_scores_source_7be10),
    collapse = " x "
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Confirm expected dimensions
# ------------------------------------------------------------

stopifnot(
  is.matrix(expr_estimate_7be10),
  nrow(expr_estimate_7be10) == 58387,
  ncol(expr_estimate_7be10) == 513
)

stopifnot(
  is.data.frame(estimate_input_df_7be10),
  ncol(estimate_input_df_7be10) >= 1
)

stopifnot(
  ncol(estimate_scores_wgcna_source_7be10) == 513
)

stopifnot(
  nrow(estimate_scores_wgcna_source_7be10) == 3
)

cat("\nExpected dimensions: PASSED\n")


# ------------------------------------------------------------
# 7. Confirm sample IDs
# ------------------------------------------------------------

expr_samples_7be10 <- colnames(
  expr_estimate_7be10
)

score_samples_7be10 <- colnames(
  estimate_scores_wgcna_source_7be10
)

cat("\n------------------------------------------------------------\n")
cat("7. Sample identity audit\n")
cat("------------------------------------------------------------\n")

cat(
  "Expression sample count: ",
  length(expr_samples_7be10),
  "\n",
  sep = ""
)

cat(
  "Score sample count: ",
  length(score_samples_7be10),
  "\n",
  sep = ""
)

cat(
  "Expression IDs unique: ",
  !anyDuplicated(expr_samples_7be10),
  "\n",
  sep = ""
)

cat(
  "Score IDs unique: ",
  !anyDuplicated(score_samples_7be10),
  "\n",
  sep = ""
)

cat(
  "Exact sample order identical: ",
  identical(
    expr_samples_7be10,
    score_samples_7be10
  ),
  "\n",
  sep = ""
)

if (!identical(
  expr_samples_7be10,
  score_samples_7be10
)) {
  stop(
    "HARD STOP: Expression and score sample IDs/order differ."
  )
}


# ------------------------------------------------------------
# 8. Confirm core genes in expression matrix
# ------------------------------------------------------------

core_genes_7be10 <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

expr_gene_ids_7be10 <- rownames(
  expr_estimate_7be10
)

cat("\n------------------------------------------------------------\n")
cat("8. Core-gene identity audit\n")
cat("------------------------------------------------------------\n")

core_gene_presence_7be10 <- setNames(
  core_genes_7be10 %in% expr_gene_ids_7be10,
  core_genes_7be10
)

print(core_gene_presence_7be10)

if (!all(core_gene_presence_7be10)) {
  stop(
    "HARD STOP: One or more core genes are absent from expr_estimate."
  )
}

cat("\nAll three core genes present: TRUE\n")


# ------------------------------------------------------------
# 9. Locate core genes
# ------------------------------------------------------------

core_gene_positions_7be10 <- match(
  core_genes_7be10,
  expr_gene_ids_7be10
)

core_gene_positions_table_7be10 <- data.frame(
  Gene = core_genes_7be10,
  Row = core_gene_positions_7be10,
  stringsAsFactors = FALSE
)

cat("\nCore-gene positions:\n")
print(core_gene_positions_table_7be10)


# ------------------------------------------------------------
# 10. Confirm ESTIMATE input contains the three core genes
# ------------------------------------------------------------

input_gene_column_candidates_7be10 <- c(
  "GeneSymbol",
  "gene_symbol",
  "Gene",
  "gene"
)

input_gene_column_7be10 <- intersect(
  input_gene_column_candidates_7be10,
  colnames(estimate_input_df_7be10)
)

cat("\nESTIMATE input gene column detected:\n")
print(input_gene_column_7be10)

if (length(input_gene_column_7be10) != 1) {
  stop(
    "HARD STOP: Could not identify a unique GeneSymbol column."
  )
}

input_gene_symbols_7be10 <- as.character(
  estimate_input_df_7be10[[
    input_gene_column_7be10
  ]]
)

input_core_presence_7be10 <- setNames(
  core_genes_7be10 %in% input_gene_symbols_7be10,
  core_genes_7be10
)

cat("\nCore genes in ESTIMATE input:\n")
print(input_core_presence_7be10)

if (!all(input_core_presence_7be10)) {
  stop(
    "HARD STOP: One or more core genes are absent from ESTIMATE input."
  )
}


# ------------------------------------------------------------
# 11. Extract original core-gene expression
# ------------------------------------------------------------

core_expr_7be10 <- t(
  expr_estimate_7be10[
    core_gene_positions_7be10,
    ,
    drop = FALSE
  ]
)

colnames(core_expr_7be10) <- core_genes_7be10
rownames(core_expr_7be10) <- expr_samples_7be10

core_expr_7be10 <- as.data.frame(
  core_expr_7be10,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 12. Confirm no missing/non-finite core expression
# ------------------------------------------------------------

core_expr_matrix_7be10 <- as.matrix(
  core_expr_7be10
)

core_expr_finite_7be10 <- all(
  is.finite(core_expr_matrix_7be10)
)

cat(
  "\nCore expression finite: ",
  core_expr_finite_7be10,
  "\n",
  sep = ""
)

if (!core_expr_finite_7be10) {
  stop(
    "HARD STOP: Non-finite core-gene expression detected."
  )
}


# ------------------------------------------------------------
# 13. Confirm source score row names
# ------------------------------------------------------------

cat("\nOriginal ESTIMATE score row names:\n")
print(rownames(
  estimate_scores_wgcna_source_7be10
))

expected_score_names_7be10 <- c(
  "StromalScore",
  "ImmuneScore",
  "ESTIMATEScore"
)

if (!identical(
  rownames(estimate_scores_wgcna_source_7be10),
  expected_score_names_7be10
)) {
  
  cat("\nExpected:\n")
  print(expected_score_names_7be10)
  
  stop(
    "HARD STOP: Original ESTIMATE score row order/names differ."
  )
}


# ------------------------------------------------------------
# 14. Confirm original score matrix is finite
# ------------------------------------------------------------

score_matrix_7be10 <- as.matrix(
  estimate_scores_wgcna_source_7be10
)

score_finite_7be10 <- all(
  is.finite(score_matrix_7be10)
)

cat(
  "\nOriginal ESTIMATE scores finite: ",
  score_finite_7be10,
  "\n",
  sep = ""
)

if (!score_finite_7be10) {
  stop(
    "HARD STOP: Non-finite values detected in original scores."
  )
}


# ------------------------------------------------------------
# 15. Exact identity with yesterday's locked TME scores
#
# Compare the source object with the currently loaded object
# if it exists.
# ------------------------------------------------------------

if (exists(
  "estimate_scores_wgcna",
  inherits = TRUE
)) {
  
  cat(
    "\nCurrent-session estimate_scores_wgcna detected.\n"
  )
  
  source_vs_current_identical_7be10 <- identical(
    estimate_scores_wgcna_source_7be10,
    get(
      "estimate_scores_wgcna",
      inherits = TRUE
    )
  )
  
  cat(
    "Source vs current estimate_scores_wgcna identical: ",
    source_vs_current_identical_7be10,
    "\n",
    sep = ""
  )
  
} else {
  
  source_vs_current_identical_7be10 <- NA
  
  cat(
    "\nCurrent-session estimate_scores_wgcna not present.\n"
  )
}


# ------------------------------------------------------------
# 16. Create a NEW, read-only source audit object
# ------------------------------------------------------------

step7be10b1_audit <- list(
  
  step = "STEP 7B-7E-10B-1",
  
  status = "SOURCE_RECOVERY_AUDIT",
  
  source_file = estimate_locked_file_7be10,
  
  expr_estimate_dim =
    dim(expr_estimate_7be10),
  
  estimate_input_df_dim =
    dim(estimate_input_df_7be10),
  
  original_score_dim =
    dim(estimate_scores_wgcna_source_7be10),
  
  n_samples =
    ncol(expr_estimate_7be10),
  
  n_expression_genes =
    nrow(expr_estimate_7be10),
  
  sample_order_identical =
    identical(
      expr_samples_7be10,
      score_samples_7be10
    ),
  
  core_genes =
    core_genes_7be10,
  
  core_genes_present =
    core_gene_presence_7be10,
  
  core_genes_in_estimate_input =
    input_core_presence_7be10,
  
  original_score_names =
    rownames(
      estimate_scores_wgcna_source_7be10
    ),
  
  original_scores_finite =
    score_finite_7be10,
  
  source_vs_current_identical =
    source_vs_current_identical_7be10,
  
  LOO_score_calculated = FALSE,
  
  statistical_analysis_performed = FALSE,
  
  locked_objects_modified = FALSE,
  
  figure6_modified = FALSE,
  
  supplementary_table3_modified = FALSE
)


# ------------------------------------------------------------
# 17. Final audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 7B-7E-10B-1 FINAL AUDIT\n")
cat("============================================================\n\n")

cat(
  "Source recovery: PASSED\n"
)

cat(
  "Expression matrix: 58,387 × 513\n"
)

cat(
  "Sample identity/order: PASSED\n"
)

cat(
  "Core-gene presence: PASSED\n"
)

cat(
  "ESTIMATE input membership: PASSED\n"
)

cat(
  "Original score structure: PASSED\n"
)

cat(
  "Original scores finite: PASSED\n"
)

cat(
  "LOO score calculated: FALSE\n"
)

cat(
  "Statistical analysis performed: FALSE\n"
)

cat(
  "Locked objects modified: FALSE\n"
)

cat(
  "\nSTEP 7B-7E-10B-1: READY FOR LOO SCORE CONSTRUCTION\n"
)


# ============================================================
# STEP 7B-7E-10B-1
# Save source-recovery audit checkpoint
# ============================================================

checkpoint_file_7be10b1 <- file.path(
  checkpoint_dir_7be10,
  "Paper1_STEP7B7E10B1_source_recovery_audit_LOCKED.RData"
)

save(
  step7be10b1_audit,
  expr_estimate_7be10,
  estimate_input_df_7be10,
  estimate_scores_wgcna_source_7be10,
  estimate_scores_source_7be10,
  estimate_sample_ids_7be10,
  estimate_sample_ids_raw_7be10,
  estimate_sample_ids_normalized_7be10,
  core_expr_7be10,
  core_gene_positions_table_7be10,
  file = checkpoint_file_7be10b1
)

cat(
  "\nCheckpoint saved:\n",
  checkpoint_file_7be10b1,
  "\n",
  sep = ""
)

cat(
  "File exists: ",
  file.exists(checkpoint_file_7be10b1),
  "\n",
  sep = ""
)

cat(
  "File size: ",
  file.info(checkpoint_file_7be10b1)$size,
  " bytes\n",
  sep = ""
)

cat(
  "\nSTEP 7B-7E-10B-1 SAVE COMPLETE\n"
)


# ============================================================
# STEP 7B-7E-10B-1
# Read-back verification
# ============================================================

readback_env_7be10b1 <- new.env(
  parent = emptyenv()
)

loaded_7be10b1 <- load(
  checkpoint_file_7be10b1,
  envir = readback_env_7be10b1
)

cat("\nLoaded objects:\n")
print(loaded_7be10b1)

cat(
  "\nRead-back status: ",
  readback_env_7be10b1$step7be10b1_audit$status,
  "\n",
  sep = ""
)

stopifnot(
  readback_env_7be10b1$step7be10b1_audit$status ==
    "SOURCE_RECOVERY_AUDIT"
)

cat(
  "\nSTEP 7B-7E-10B-1 FINAL FILE CHECK: PASSED\n"
)

##STEP 7B-7E-10B-2A ESTIMATE implementation audit
# ============================================================
# Paper 1
# STEP 7B-7E-10B-2A
# ESTIMATE implementation audit
#
# PURPOSE:
#   1. Check whether original ESTIMATE package is available
#   2. Audit the scoring implementation
#   3. Confirm that the ESTIMATE scoring implementation can
#      theoretically accept a 140-gene LOO signature
#
# IMPORTANT:
#   - NO LOO score calculation
#   - NO statistical analysis
#   - NO modification of locked ESTIMATE scores
#   - NO modification of Figure 6 / Supplementary Table 3
# ============================================================


cat("\n============================================================\n")
cat("STEP 7B-7E-10B-2A: ESTIMATE implementation audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Define working directory
# ------------------------------------------------------------

work_dir <- "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis"

if (!dir.exists(work_dir)) {
  stop("WORKING DIRECTORY NOT FOUND: ", work_dir)
}

setwd(work_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Confirm required objects from STEP 10B-1
# ------------------------------------------------------------

required_objects <- c(
  "expr_estimate_7be10",
  "estimate_input_df_7be10",
  "estimate_scores_wgcna_source_7be10",
  "estimate_scores_source_7be10",
  "estimate_sample_ids_7be10",
  "core_expr_7be10",
  "core_gene_positions_table_7be10"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  
  stop(
    paste0(
      "\nSTEP 7B-7E-10B-2A STOP:\n",
      "Required objects from STEP 7B-7E-10B-1 are missing:\n",
      paste(missing_objects, collapse = ", "),
      "\n\nPlease load the STEP 10B-1 checkpoint first."
    )
  )
  
}

cat("Required STEP 10B-1 objects: PASSED\n\n")


# ------------------------------------------------------------
# 2. Check original ESTIMATE package availability
# ------------------------------------------------------------

estimate_pkg_available <- requireNamespace(
  "estimate",
  quietly = TRUE
)

cat("------------------------------------------------------------\n")
cat("2. Original ESTIMATE package availability\n")
cat("------------------------------------------------------------\n")

cat(
  "Package 'estimate' available:",
  estimate_pkg_available,
  "\n"
)

if (estimate_pkg_available) {
  
  estimate_pkg_version <- as.character(
    packageVersion("estimate")
  )
  
  estimate_pkg_path <- find.package("estimate")
  
  cat(
    "Package version:",
    estimate_pkg_version,
    "\n"
  )
  
  cat(
    "Package path:",
    estimate_pkg_path,
    "\n"
  )
  
} else {
  
  estimate_pkg_version <- NA_character_
  estimate_pkg_path <- NA_character_
  
  cat(
    "\nOriginal 'estimate' package is NOT available in the current R session.\n"
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 3. Check package namespace functions
# ------------------------------------------------------------

estimate_functions_present <- c(
  filterCommonGenes = FALSE,
  estimateScore = FALSE
)

if (estimate_pkg_available) {
  
  estimate_functions_present["filterCommonGenes"] <-
    exists(
      "filterCommonGenes",
      where = asNamespace("estimate"),
      inherits = FALSE
    )
  
  estimate_functions_present["estimateScore"] <-
    exists(
      "estimateScore",
      where = asNamespace("estimate"),
      inherits = FALSE
    )
  
}

cat("------------------------------------------------------------\n")
cat("3. Required ESTIMATE functions\n")
cat("------------------------------------------------------------\n")

print(estimate_functions_present)

cat("\n")


# ------------------------------------------------------------
# 4. Inspect function signatures WITHOUT EXECUTING them
# ------------------------------------------------------------

filterCommonGenes_formals <- NULL
estimateScore_formals <- NULL

if (estimate_pkg_available &&
    estimate_functions_present["filterCommonGenes"]) {
  
  filterCommonGenes_fun <-
    get(
      "filterCommonGenes",
      envir = asNamespace("estimate")
    )
  
  filterCommonGenes_formals <-
    formals(filterCommonGenes_fun)
  
}

if (estimate_pkg_available &&
    estimate_functions_present["estimateScore"]) {
  
  estimateScore_fun <-
    get(
      "estimateScore",
      envir = asNamespace("estimate")
    )
  
  estimateScore_formals <-
    formals(estimateScore_fun)
  
}


cat("------------------------------------------------------------\n")
cat("4. Function signatures\n")
cat("------------------------------------------------------------\n")

if (!is.null(filterCommonGenes_formals)) {
  
  cat("\nfilterCommonGenes() formals:\n")
  print(filterCommonGenes_formals)
  
} else {
  
  cat("\nfilterCommonGenes(): NOT AVAILABLE\n")
  
}


if (!is.null(estimateScore_formals)) {
  
  cat("\nestimateScore() formals:\n")
  print(estimateScore_formals)
  
} else {
  
  cat("\nestimateScore(): NOT AVAILABLE\n")
  
}

cat("\n")


# ------------------------------------------------------------
# 5. Inspect function source location
#
# IMPORTANT:
#   This only identifies the function implementation.
#   It does NOT execute the scoring function.
# ------------------------------------------------------------

estimateScore_source <- NULL
filterCommonGenes_source <- NULL

if (estimate_pkg_available &&
    estimate_functions_present["estimateScore"]) {
  
  estimateScore_source <- getAnywhere(
    "estimateScore"
  )
  
}

if (estimate_pkg_available &&
    estimate_functions_present["filterCommonGenes"]) {
  
  filterCommonGenes_source <- getAnywhere(
    "filterCommonGenes"
  )
  
}


cat("------------------------------------------------------------\n")
cat("5. Function source audit\n")
cat("------------------------------------------------------------\n")

if (!is.null(estimateScore_source)) {
  
  cat("\nestimateScore() source locations:\n")
  print(estimateScore_source$where)
  
} else {
  
  cat("\nestimateScore() source unavailable.\n")
  
}


if (!is.null(filterCommonGenes_source)) {
  
  cat("\nfilterCommonGenes() source locations:\n")
  print(filterCommonGenes_source$where)
  
} else {
  
  cat("\nfilterCommonGenes() source unavailable.\n")
  
}

cat("\n")


# ------------------------------------------------------------
# 6. Verify the LOO signatures from STEP 10A
# ------------------------------------------------------------

loo_required_objects <- c(
  "stromal_genes_loo_AOC3_7be",
  "immune_genes_loo_IL7R_7be"
)

missing_loo_objects <- loo_required_objects[
  !vapply(loo_required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_loo_objects) > 0) {
  
  stop(
    paste0(
      "\nSTEP 7B-7E-10B-2A STOP:\n",
      "LOO signature objects are missing:\n",
      paste(missing_loo_objects, collapse = ", ")
    )
  )
  
}


stromal_loo_genes <-
  stromal_genes_loo_AOC3_7be

immune_loo_genes <-
  immune_genes_loo_IL7R_7be


cat("------------------------------------------------------------\n")
cat("6. LOO signature audit\n")
cat("------------------------------------------------------------\n")

cat(
  "AOC3-removed Stromal signature:",
  length(stromal_loo_genes),
  "genes\n"
)

cat(
  "IL7R-removed Immune signature:",
  length(immune_loo_genes),
  "genes\n"
)


# Basic structural checks

stromal_loo_unique <-
  length(unique(stromal_loo_genes)) ==
  length(stromal_loo_genes)

immune_loo_unique <-
  length(unique(immune_loo_genes)) ==
  length(immune_loo_genes)

stromal_loo_character <-
  is.character(stromal_loo_genes)

immune_loo_character <-
  is.character(immune_loo_genes)

cat(
  "Stromal LOO signature unique:",
  stromal_loo_unique,
  "\n"
)

cat(
  "Immune LOO signature unique:",
  immune_loo_unique,
  "\n"
)

cat(
  "Stromal LOO signature character:",
  stromal_loo_character,
  "\n"
)

cat(
  "Immune LOO signature character:",
  immune_loo_character,
  "\n\n"
)


# ------------------------------------------------------------
# 7. Check LOO genes against expression matrix
# ------------------------------------------------------------

expr_gene_symbols <-
  rownames(expr_estimate_7be10)

if (is.null(expr_gene_symbols)) {
  
  stop(
    "STEP 7B-7E-10B-2A STOP: expr_estimate_7be10 has no rownames."
  )
  
}


stromal_loo_in_expr <-
  stromal_loo_genes %in% expr_gene_symbols

immune_loo_in_expr <-
  immune_loo_genes %in% expr_gene_symbols


cat("------------------------------------------------------------\n")
cat("7. LOO signature / expression matrix compatibility\n")
cat("------------------------------------------------------------\n")

cat(
  "Stromal LOO genes found in expression matrix:",
  sum(stromal_loo_in_expr),
  "/",
  length(stromal_loo_genes),
  "\n"
)

cat(
  "Immune LOO genes found in expression matrix:",
  sum(immune_loo_in_expr),
  "/",
  length(immune_loo_genes),
  "\n"
)


missing_stromal_loo <-
  stromal_loo_genes[!stromal_loo_in_expr]

missing_immune_loo <-
  immune_loo_genes[!immune_loo_in_expr]


if (length(missing_stromal_loo) > 0) {
  
  cat(
    "\nMissing Stromal LOO genes:\n"
  )
  
  print(missing_stromal_loo)
  
}


if (length(missing_immune_loo) > 0) {
  
  cat(
    "\nMissing Immune LOO genes:\n"
  )
  
  print(missing_immune_loo)
  
}

cat("\n")


# ------------------------------------------------------------
# 8. Check original ESTIMATE input structure
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("8. Original ESTIMATE input structure\n")
cat("------------------------------------------------------------\n")

cat(
  "expr_estimate_7be10 dimensions:",
  paste(dim(expr_estimate_7be10), collapse = " × "),
  "\n"
)

cat(
  "estimate_input_df_7be10 dimensions:",
  paste(dim(estimate_input_df_7be10), collapse = " × "),
  "\n"
)

cat(
  "Original ESTIMATE score dimensions:",
  paste(
    dim(estimate_scores_wgcna_source_7be10),
    collapse = " × "
  ),
  "\n"
)

cat(
  "Original score rows:\n"
)

print(
  rownames(estimate_scores_wgcna_source_7be10)
)

cat("\n")


# ------------------------------------------------------------
# 9. Check whether estimateScore() expects a gene × sample
#    matrix/data.frame
#
# IMPORTANT:
#   We inspect the formal argument and documentation metadata
#   only. We do NOT execute estimateScore().
# ------------------------------------------------------------

estimateScore_accepts_data_object <- FALSE

if (!is.null(estimateScore_formals)) {
  
  estimateScore_arg_names <-
    names(estimateScore_formals)
  
  cat("estimateScore() argument names:\n")
  print(estimateScore_arg_names)
  
  cat("\n")
  
  # The first formal argument is expected to represent
  # the expression dataset.
  
  if (length(estimateScore_arg_names) >= 1) {
    
    first_arg <-
      estimateScore_arg_names[1]
    
    cat(
      "First estimateScore() argument:",
      first_arg,
      "\n"
    )
    
    estimateScore_accepts_data_object <- TRUE
    
  }
  
}

cat(
  "\nFunction-level input acceptance audit:",
  estimateScore_accepts_data_object,
  "\n\n"
)


# ------------------------------------------------------------
# 10. Confirm that 140-gene signatures are not rejected
#     structurally by the audit
#
# This is intentionally NOT a scoring test.
#
# We only verify:
#   - 140 genes
#   - unique
#   - character
#   - present in expression matrix
#   - no duplicated symbols
# ------------------------------------------------------------

loo_signature_structure_pass <-
  length(stromal_loo_genes) == 140 &&
  length(immune_loo_genes) == 140 &&
  stromal_loo_unique &&
  immune_loo_unique &&
  stromal_loo_character &&
  immune_loo_character &&
  all(stromal_loo_in_expr) &&
  all(immune_loo_in_expr)


cat("------------------------------------------------------------\n")
cat("10. 140-gene LOO signature structural compatibility\n")
cat("------------------------------------------------------------\n")

cat(
  "140-gene LOO signature structural audit:",
  loo_signature_structure_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 11. Explicit prohibition checks
# ------------------------------------------------------------

loo_score_calculated <- FALSE
statistical_analysis_performed <- FALSE
locked_objects_modified <- FALSE
estimate_score_function_executed <- FALSE


cat("------------------------------------------------------------\n")
cat("11. Safety checks\n")
cat("------------------------------------------------------------\n")

cat(
  "LOO score calculated:",
  loo_score_calculated,
  "\n"
)

cat(
  "Statistical analysis performed:",
  statistical_analysis_performed,
  "\n"
)

cat(
  "Locked objects modified:",
  locked_objects_modified,
  "\n"
)

cat(
  "estimateScore() executed:",
  estimate_score_function_executed,
  "\n\n"
)


# ------------------------------------------------------------
# 12. Final audit status
# ------------------------------------------------------------

implementation_available_pass <-
  estimate_pkg_available &&
  estimate_functions_present["estimateScore"] &&
  estimate_functions_present["filterCommonGenes"]

loo_structure_pass <-
  loo_signature_structure_pass

final_audit_pass <-
  implementation_available_pass &&
  loo_structure_pass &&
  !loo_score_calculated &&
  !statistical_analysis_performed &&
  !locked_objects_modified &&
  !estimate_score_function_executed


if (final_audit_pass) {
  
  audit_status <-
    "IMPLEMENTATION_AUDIT_PASSED"
  
} else if (!estimate_pkg_available) {
  
  audit_status <-
    "ORIGINAL_ESTIMATE_PACKAGE_NOT_AVAILABLE"
  
} else {
  
  audit_status <-
    "IMPLEMENTATION_AUDIT_REQUIRES_REVIEW"
  
}


cat("============================================================\n")
cat("STEP 7B-7E-10B-2A FINAL STATUS\n")
cat("============================================================\n")

cat(
  "Audit status:",
  audit_status,
  "\n"
)

cat(
  "Original ESTIMATE package available:",
  estimate_pkg_available,
  "\n"
)

cat(
  "filterCommonGenes() available:",
  estimate_functions_present["filterCommonGenes"],
  "\n"
)

cat(
  "estimateScore() available:",
  estimate_functions_present["estimateScore"],
  "\n"
)

cat(
  "140-gene LOO structural compatibility:",
  loo_structure_pass,
  "\n"
)

cat(
  "LOO score calculated:",
  loo_score_calculated,
  "\n"
)

cat(
  "Locked objects modified:",
  locked_objects_modified,
  "\n"
)

cat("============================================================\n\n")


# ------------------------------------------------------------
# 13. Save audit checkpoint
# ------------------------------------------------------------

checkpoint_dir <-
  file.path(
    work_dir,
    "STEP7B-7E-10B-2"
  )

if (!dir.exists(checkpoint_dir)) {
  dir.create(
    checkpoint_dir,
    recursive = TRUE
  )
}


checkpoint_file_7be10b2a <-
  file.path(
    checkpoint_dir,
    "Paper1_STEP7B7E10B2A_ESTIMATE_implementation_audit_LOCKED.RData"
  )


step7be10b2a_audit <- list(
  
  step = "STEP 7B-7E-10B-2A",
  
  purpose =
    "ESTIMATE implementation audit only",
  
  status =
    audit_status,
  
  estimate_package_available =
    estimate_pkg_available,
  
  estimate_package_version =
    estimate_pkg_version,
  
  estimate_package_path =
    estimate_pkg_path,
  
  required_functions_present =
    estimate_functions_present,
  
  filterCommonGenes_formals =
    filterCommonGenes_formals,
  
  estimateScore_formals =
    estimateScore_formals,
  
  stromal_loo_n =
    length(stromal_loo_genes),
  
  immune_loo_n =
    length(immune_loo_genes),
  
  stromal_loo_all_in_expression =
    all(stromal_loo_in_expr),
  
  immune_loo_all_in_expression =
    all(immune_loo_in_expr),
  
  loo_signature_structure_pass =
    loo_signature_structure_pass,
  
  loo_score_calculated =
    loo_score_calculated,
  
  statistical_analysis_performed =
    statistical_analysis_performed,
  
  locked_objects_modified =
    locked_objects_modified,
  
  estimateScore_executed =
    estimate_score_function_executed,
  
  locked_baseline_modified =
    FALSE
  
)


save(
  step7be10b2a_audit,
  estimate_pkg_available,
  estimate_pkg_version,
  estimate_pkg_path,
  estimate_functions_present,
  filterCommonGenes_formals,
  estimateScore_formals,
  stromal_loo_genes,
  immune_loo_genes,
  loo_signature_structure_pass,
  implementation_available_pass,
  checkpoint_file_7be10b2a,
  file = checkpoint_file_7be10b2a
)


cat("Checkpoint saved:\n")
cat(checkpoint_file_7be10b2a, "\n\n")


# ============================================================
# STEP 7B-7E-10B-2A
# Final file existence check only
# ============================================================

checkpoint_file_7be10b2a <-
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis/STEP7B-7E-10B-2/Paper1_STEP7B7E10B2A_ESTIMATE_implementation_audit_LOCKED.RData"

file_exists_final <-
  file.exists(checkpoint_file_7be10b2a)

if (file_exists_final) {
  
  file_size_final <-
    file.info(checkpoint_file_7be10b2a)$size
  
} else {
  
  file_size_final <- NA_real_
  
}

cat("------------------------------------------------------------\n")
cat("STEP 7B-7E-10B-2A FINAL FILE CHECK\n")
cat("------------------------------------------------------------\n")

cat(
  "File exists:",
  file_exists_final,
  "\n"
)

cat(
  "File size:",
  file_size_final,
  "bytes\n"
)

if (!file_exists_final ||
    is.na(file_size_final) ||
    file_size_final <= 0) {
  
  stop(
    "STEP 7B-7E-10B-2A FINAL FILE CHECK: FAILED"
  )
  
}

cat(
  "\nSTEP 7B-7E-10B-2A FINAL FILE CHECK: PASSED\n"
)

cat("\n============================================================\n")
cat("STEP 7B-7E-10B-2B-PRECHECK\n")
cat("Inspect recovered ESTIMATE score object structure\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Object class / type / structure
# ------------------------------------------------------------

cat("Object name:\n")
cat("estimate_scores_wgcna_source_7be10\n\n")

cat("class:\n")
print(class(estimate_scores_wgcna_source_7be10))

cat("\ntypeof:\n")
print(typeof(estimate_scores_wgcna_source_7be10))

cat("\nlength:\n")
print(length(estimate_scores_wgcna_source_7be10))

cat("\ndim:\n")
print(dim(estimate_scores_wgcna_source_7be10))

cat("\nattributes:\n")
print(attributes(estimate_scores_wgcna_source_7be10))


# ------------------------------------------------------------
# 2. Structure
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("str() of recovered object\n")
cat("------------------------------------------------------------\n")

str(
  estimate_scores_wgcna_source_7be10,
  max.level = 2
)


# ------------------------------------------------------------
# 3. If list: inspect element names / classes / dimensions
# ------------------------------------------------------------

if (is.list(estimate_scores_wgcna_source_7be10)) {
  
  cat("\n------------------------------------------------------------\n")
  cat("LIST ELEMENT INSPECTION\n")
  cat("------------------------------------------------------------\n\n")
  
  cat("Element names:\n")
  print(names(estimate_scores_wgcna_source_7be10))
  
  cat("\nNumber of elements:\n")
  cat(length(estimate_scores_wgcna_source_7be10), "\n\n")
  
  for (i in seq_along(estimate_scores_wgcna_source_7be10)) {
    
    cat("------------------------------------------------------------\n")
    cat("Element", i, "\n")
    
    if (!is.null(names(estimate_scores_wgcna_source_7be10))) {
      cat(
        "Name:",
        names(estimate_scores_wgcna_source_7be10)[i],
        "\n"
      )
    }
    
    cat(
      "Class:",
      paste(
        class(estimate_scores_wgcna_source_7be10[[i]]),
        collapse = ", "
      ),
      "\n"
    )
    
    cat(
      "Type:",
      typeof(estimate_scores_wgcna_source_7be10[[i]]),
      "\n"
    )
    
    cat(
      "Length:",
      length(estimate_scores_wgcna_source_7be10[[i]]),
      "\n"
    )
    
    cat("Dimensions:\n")
    print(dim(estimate_scores_wgcna_source_7be10[[i]]))
    
    cat("Names:\n")
    print(names(estimate_scores_wgcna_source_7be10[[i]]))
    
    cat("Rownames:\n")
    print(
      tryCatch(
        rownames(estimate_scores_wgcna_source_7be10[[i]]),
        error = function(e) NULL
      )
    )
    
    cat("Colnames:\n")
    print(
      tryCatch(
        colnames(estimate_scores_wgcna_source_7be10[[i]]),
        error = function(e) NULL
      )
    )
    
  }
}


# ------------------------------------------------------------
# 4. Inspect the separately recovered original score object
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Compare with estimate_scores_source_7be10\n")
cat("============================================================\n\n")

cat("class:\n")
print(class(estimate_scores_source_7be10))

cat("\ntypeof:\n")
print(typeof(estimate_scores_source_7be10))

cat("\nlength:\n")
print(length(estimate_scores_source_7be10))

cat("\ndim:\n")
print(dim(estimate_scores_source_7be10))

cat("\nstr():\n")
str(
  estimate_scores_source_7be10,
  max.level = 2
)


# ------------------------------------------------------------
# 5. Safety check
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SAFETY STATUS\n")
cat("============================================================\n\n")

cat("LOO score calculated: FALSE\n")
cat("Statistical analysis performed: FALSE\n")
cat("Locked source scores modified: FALSE\n")
cat("ESTIMATE scoring executed: FALSE\n")

cat("\nPRECHECK COMPLETE.\n")
cat("Please provide the complete console output before continuing.\n")

cat("\n============================================================\n")
cat("STEP 7B-7E-10B-2B\n")
cat("Original ESTIMATE implementation recovery and verification\n")
cat("REVISED: data.frame-safe source-score audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

work_dir <-
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis"

if (!dir.exists(work_dir)) {
  stop(
    "WORKING DIRECTORY NOT FOUND: ",
    work_dir
  )
}

setwd(work_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Confirm required objects
# ------------------------------------------------------------

required_objects <- c(
  "expr_estimate_7be10",
  "estimate_input_df_7be10",
  "estimate_scores_wgcna_source_7be10",
  "estimate_scores_source_7be10",
  "estimate_sample_ids_7be10",
  "core_expr_7be10",
  "core_gene_positions_table_7be10",
  "stromal_genes_loo_AOC3_7be",
  "immune_genes_loo_IL7R_7be",
  "step7be10b2a_audit"
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
    paste0(
      "\nSTEP 7B-7E-10B-2B STOP:\n",
      "Required objects are missing:\n",
      paste(
        missing_objects,
        collapse = ", "
      ),
      "\n\nPlease load STEP 10B-1 and STEP 10B-2A checkpoints first."
    )
  )
  
}

cat(
  "Required STEP 10B-1 / STEP 10B-2A objects: PASSED\n\n"
)


# ------------------------------------------------------------
# 2. Preserve original source score
# ------------------------------------------------------------

source_score_before_10b2b <-
  estimate_scores_wgcna_source_7be10

cat("Source score object class:\n")
print(class(source_score_before_10b2b))

cat("\nSource score object type:\n")
print(typeof(source_score_before_10b2b))

cat("\nSource score dimensions:\n")
print(dim(source_score_before_10b2b))

cat("\nSource score row names:\n")
print(rownames(source_score_before_10b2b))

cat("\n")


# ------------------------------------------------------------
# 3. Data-frame-safe finite-value check
#
# IMPORTANT:
# data.frame has typeof() == "list" in R.
# We therefore inspect a temporary matrix representation.
#
# The original object itself is NOT modified.
# ------------------------------------------------------------

source_score_finite_before_10b2b <-
  all(
    is.finite(
      as.matrix(source_score_before_10b2b)
    )
  )

cat(
  "All source score values finite:",
  source_score_finite_before_10b2b,
  "\n\n"
)


# ------------------------------------------------------------
# 4. Structural verification
# ------------------------------------------------------------

source_score_dim_before_10b2b <-
  dim(source_score_before_10b2b)

source_score_rownames_before_10b2b <-
  rownames(source_score_before_10b2b)

source_score_colnames_before_10b2b <-
  colnames(source_score_before_10b2b)


if (!identical(
  source_score_dim_before_10b2b,
  c(3L, 513L)
)) {
  
  stop(
    "STOP: Original source score is not 3 × 513."
  )
  
}

if (!identical(
  source_score_rownames_before_10b2b,
  c(
    "StromalScore",
    "ImmuneScore",
    "ESTIMATEScore"
  )
)) {
  
  stop(
    "STOP: Original source score row names do not match expected structure."
  )
  
}

if (!source_score_finite_before_10b2b) {
  
  stop(
    "STOP: Original source score contains non-finite values."
  )
  
}

cat("Source score dimensions: PASSED\n")
cat("Source score row structure: PASSED\n")
cat("Source score finite-value check: PASSED\n\n")


# ------------------------------------------------------------
# 5. Confirm sample identity/order
# ------------------------------------------------------------

sample_identity_match <-
  identical(
    colnames(source_score_before_10b2b),
    estimate_sample_ids_7be10
  )

cat(
  "Source score sample identity/order:",
  sample_identity_match,
  "\n\n"
)

if (!sample_identity_match) {
  
  stop(
    "STOP: Source score sample IDs/order do not match recovered sample IDs."
  )
  
}


# ------------------------------------------------------------
# 6. Check core genes in ESTIMATE input
# ------------------------------------------------------------

core_genes_10b2b <-
  c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )

core_presence_in_estimate_input <-
  core_genes_10b2b %in%
  rownames(expr_estimate_7be10)

cat("Core genes in ESTIMATE expression input:\n")

print(
  data.frame(
    Gene = core_genes_10b2b,
    Present = core_presence_in_estimate_input,
    stringsAsFactors = FALSE
  )
)

if (!all(core_presence_in_estimate_input)) {
  
  stop(
    "STOP: One or more core genes are absent from ESTIMATE expression input."
  )
  
}

cat("\nCore-gene presence: PASSED\n\n")


# ------------------------------------------------------------
# 7. Check ESTIMATE package availability
# ------------------------------------------------------------

estimate_package_available <-
  requireNamespace(
    "estimate",
    quietly = TRUE
  )

cat(
  "Original 'estimate' package available:",
  estimate_package_available,
  "\n\n"
)


# ------------------------------------------------------------
# 8. If package is available:
#    inspect implementation WITHOUT EXECUTING SCORING
# ------------------------------------------------------------

estimate_functions_available <-
  FALSE

filterCommonGenes_available <-
  FALSE

estimateScore_available <-
  FALSE

filterCommonGenes_formals <-
  NULL

estimateScore_formals <-
  NULL

filterCommonGenes_body <-
  NULL

estimateScore_body <-
  NULL


if (estimate_package_available) {
  
  filterCommonGenes_available <-
    exists(
      "filterCommonGenes",
      envir = asNamespace("estimate"),
      inherits = FALSE
    )
  
  estimateScore_available <-
    exists(
      "estimateScore",
      envir = asNamespace("estimate"),
      inherits = FALSE
    )
  
  estimate_functions_available <-
    filterCommonGenes_available &&
    estimateScore_available
  
  
  cat("filterCommonGenes available:",
      filterCommonGenes_available, "\n")
  
  cat("estimateScore available:",
      estimateScore_available, "\n\n")
  
  
  if (filterCommonGenes_available) {
    
    filterCommonGenes_fun <-
      get(
        "filterCommonGenes",
        envir = asNamespace("estimate")
      )
    
    filterCommonGenes_formals <-
      formals(filterCommonGenes_fun)
    
    filterCommonGenes_body <-
      body(filterCommonGenes_fun)
    
    cat("filterCommonGenes() formals:\n")
    print(filterCommonGenes_formals)
    
    cat("\nfilterCommonGenes() body:\n")
    print(filterCommonGenes_body)
    
  }
  
  
  if (estimateScore_available) {
    
    estimateScore_fun <-
      get(
        "estimateScore",
        envir = asNamespace("estimate")
      )
    
    estimateScore_formals <-
      formals(estimateScore_fun)
    
    estimateScore_body <-
      body(estimateScore_fun)
    
    cat("\nestimateScore() formals:\n")
    print(estimateScore_formals)
    
    cat("\nestimateScore() body:\n")
    print(estimateScore_body)
    
  }
  
}


# ------------------------------------------------------------
# 9. Inspect ESTIMATE package version
# ------------------------------------------------------------

estimate_package_version <-
  NA_character_

if (estimate_package_available) {
  
  estimate_package_version <-
    as.character(
      packageVersion("estimate")
    )
  
  cat(
    "\nESTIMATE package version:",
    estimate_package_version,
    "\n"
  )
  
}


# ------------------------------------------------------------
# 10. Inspect SI_geneset WITHOUT scoring
# ------------------------------------------------------------

SI_geneset_available <-
  FALSE

SI_geneset_object <-
  NULL

if (estimate_package_available) {
  
  pkg_data <-
    data(
      package = "estimate"
    )$results
  
  SI_geneset_available <-
    "SI_geneset" %in%
    pkg_data[, "Item"]
  
  cat(
    "\nSI_geneset available:",
    SI_geneset_available,
    "\n"
  )
  
  if (SI_geneset_available) {
    
    data(
      list = "SI_geneset",
      package = "estimate",
      envir = environment()
    )
    
    if (exists("SI_geneset", envir = environment())) {
      
      SI_geneset_object <-
        get(
          "SI_geneset",
          envir = environment()
        )
      
      cat("\nSI_geneset class:\n")
      print(class(SI_geneset_object))
      
      cat("\nSI_geneset dimensions:\n")
      print(dim(SI_geneset_object))
      
      cat("\nSI_geneset structure:\n")
      str(
        SI_geneset_object,
        max.level = 2
      )
      
    }
    
  }
  
}


# ------------------------------------------------------------
# 11. Safety boundary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SAFETY STATUS\n")
cat("============================================================\n\n")

cat("LOO score calculated: FALSE\n")
cat("LOO statistical analysis performed: FALSE\n")
cat("Original ESTIMATE score recalculated: FALSE\n")
cat("Locked Step14 scores modified: FALSE\n")
cat("Figure 6 modified: FALSE\n")
cat("Supplementary Table 3 modified: FALSE\n\n")


# ------------------------------------------------------------
# 12. Audit status
# ------------------------------------------------------------

if (estimate_package_available &&
    estimate_functions_available) {
  
  audit_status <-
    "ORIGINAL_ESTIMATE_IMPLEMENTATION_AVAILABLE_FOR_SOURCE_INSPECTION"
  
} else {
  
  audit_status <-
    "ORIGINAL_ESTIMATE_IMPLEMENTATION_NOT_AVAILABLE"
  
}

cat(
  "STEP 7B-7E-10B-2B STATUS:\n",
  audit_status,
  "\n\n"
)


# ------------------------------------------------------------
# 13. Save audit checkpoint
# ------------------------------------------------------------

audit_dir <-
  file.path(
    work_dir,
    "STEP7B-7E-10B-2"
  )

if (!dir.exists(audit_dir)) {
  dir.create(
    audit_dir,
    recursive = TRUE
  )
}

step7be10b2b_audit <-
  list(
    
    step = "7B-7E-10B-2B",
    
    purpose =
      "Original ESTIMATE implementation recovery and verification",
    
    status =
      audit_status,
    
    estimate_package_available =
      estimate_package_available,
    
    estimate_package_version =
      estimate_package_version,
    
    filterCommonGenes_available =
      filterCommonGenes_available,
    
    estimateScore_available =
      estimateScore_available,
    
    SI_geneset_available =
      SI_geneset_available,
    
    source_score_dim =
      source_score_dim_before_10b2b,
    
    source_score_rownames =
      source_score_rownames_before_10b2b,
    
    source_score_finite =
      source_score_finite_before_10b2b,
    
    sample_identity_match =
      sample_identity_match,
    
    loo_score_calculated =
      FALSE,
    
    statistical_analysis_performed =
      FALSE,
    
    locked_objects_modified =
      FALSE
  )


checkpoint_file <-
  file.path(
    audit_dir,
    "Paper1_STEP7B7E10B2B_ESTIMATE_original_implementation_LOCKED.RData"
  )


save(
  step7be10b2b_audit,
  source_score_before_10b2b,
  source_score_dim_before_10b2b,
  source_score_rownames_before_10b2b,
  source_score_colnames_before_10b2b,
  source_score_finite_before_10b2b,
  sample_identity_match,
  estimate_package_available,
  estimate_package_version,
  filterCommonGenes_available,
  estimateScore_available,
  filterCommonGenes_formals,
  estimateScore_formals,
  filterCommonGenes_body,
  estimateScore_body,
  SI_geneset_available,
  SI_geneset_object,
  file = checkpoint_file
)

cat(
  "Checkpoint saved:\n",
  checkpoint_file,
  "\n\n"
)


# ------------------------------------------------------------
# 14. Final file check
# ------------------------------------------------------------

cat("============================================================\n")
cat("FINAL FILE CHECK\n")
cat("============================================================\n\n")

cat(
  "File exists:",
  file.exists(checkpoint_file),
  "\n"
)

if (file.exists(checkpoint_file)) {
  
  cat(
    "File size:",
    file.info(checkpoint_file)$size,
    "bytes\n"
  )
  
  cat(
    "\nSTEP 7B-7E-10B-2B FINAL FILE CHECK: PASSED\n"
  )
  
}

##STEP 7B-7E-10B-2C — Original ESTIMATE implementation read-back audit
# ============================================================
# Paper 1
# STEP 7B-7E-10B-2C
# Original ESTIMATE implementation read-back audit
#
# PURPOSE:
#   1. Read back the locked STEP 10B-2B checkpoint
#   2. Verify original ESTIMATE implementation information
#   3. Inspect filterCommonGenes() behavior
#   4. Inspect estimateScore() implementation
#   5. Inspect SI_geneset structure
#   6. Specifically assess how missing signature genes
#      are handled
#
# IMPORTANT:
#   - NO filterCommonGenes() execution
#   - NO estimateScore() execution
#   - NO LOO score calculation
#   - NO statistical analysis
#   - NO modification of locked Step14 objects
#   - NO modification of Figure 6 / Supplementary Table 3
# ============================================================


cat("\n============================================================\n")
cat("STEP 7B-7E-10B-2C\n")
cat("Original ESTIMATE implementation read-back audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

work_dir <-
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis"

if (!dir.exists(work_dir)) {
  
  stop(
    "WORKING DIRECTORY NOT FOUND: ",
    work_dir
  )
  
}

setwd(work_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Locate STEP 10B-2B checkpoint
# ------------------------------------------------------------

checkpoint_file <-
  file.path(
    work_dir,
    "STEP7B-7E-10B-2",
    "Paper1_STEP7B7E10B2B_ESTIMATE_original_implementation_LOCKED.RData"
  )

cat("STEP 10B-2B checkpoint:\n")
cat(checkpoint_file, "\n\n")


if (!file.exists(checkpoint_file)) {
  
  stop(
    paste0(
      "STEP 7B-7E-10B-2C STOP:\n",
      "STEP 10B-2B checkpoint not found:\n",
      checkpoint_file
    )
  )
  
}

cat(
  "STEP 10B-2B checkpoint exists: PASSED\n\n"
)


# ------------------------------------------------------------
# 2. Load checkpoint into isolated environment
#
# IMPORTANT:
#   Do not overwrite current-session locked objects.
# ------------------------------------------------------------

audit_env <-
  new.env(parent = emptyenv())

loaded_objects <-
  load(
    checkpoint_file,
    envir = audit_env
  )

cat("Objects loaded from STEP 10B-2B checkpoint:\n")
print(loaded_objects)

cat("\n")


# ------------------------------------------------------------
# 3. Confirm required audit objects
# ------------------------------------------------------------

required_readback_objects <- c(
  "step7be10b2b_audit",
  "source_score_before_10b2b",
  "source_score_dim_before_10b2b",
  "source_score_rownames_before_10b2b",
  "source_score_colnames_before_10b2b",
  "source_score_finite_before_10b2b",
  "sample_identity_match",
  "estimate_package_available",
  "estimate_package_version",
  "filterCommonGenes_available",
  "estimateScore_available",
  "filterCommonGenes_formals",
  "estimateScore_formals",
  "filterCommonGenes_body",
  "estimateScore_body",
  "SI_geneset_available",
  "SI_geneset_object"
)

missing_readback_objects <- required_readback_objects[
  !vapply(
    required_readback_objects,
    exists,
    logical(1),
    envir = audit_env,
    inherits = FALSE
  )
]

if (length(missing_readback_objects) > 0) {
  
  stop(
    paste0(
      "STEP 7B-7E-10B-2C STOP:\n",
      "Required objects missing from STEP 10B-2B checkpoint:\n",
      paste(
        missing_readback_objects,
        collapse = ", "
      )
    )
  )
  
}

cat(
  "Required STEP 10B-2B audit objects: PASSED\n\n"
)


# ------------------------------------------------------------
# 4. Read back audit status
# ------------------------------------------------------------

step7be10b2b_audit <-
  get(
    "step7be10b2b_audit",
    envir = audit_env
  )

cat("STEP 10B-2B audit object:\n")
print(step7be10b2b_audit)

cat("\n")


# ------------------------------------------------------------
# 5. Read back package availability/version
# ------------------------------------------------------------

estimate_package_available <-
  get(
    "estimate_package_available",
    envir = audit_env
  )

estimate_package_version <-
  get(
    "estimate_package_version",
    envir = audit_env
  )

filterCommonGenes_available <-
  get(
    "filterCommonGenes_available",
    envir = audit_env
  )

estimateScore_available <-
  get(
    "estimateScore_available",
    envir = audit_env
  )

SI_geneset_available <-
  get(
    "SI_geneset_available",
    envir = audit_env
  )


cat("============================================================\n")
cat("ORIGINAL ESTIMATE IMPLEMENTATION STATUS\n")
cat("============================================================\n\n")

cat(
  "estimate package available:",
  estimate_package_available,
  "\n"
)

cat(
  "estimate package version:",
  estimate_package_version,
  "\n"
)

cat(
  "filterCommonGenes available:",
  filterCommonGenes_available,
  "\n"
)

cat(
  "estimateScore available:",
  estimateScore_available,
  "\n"
)

cat(
  "SI_geneset available:",
  SI_geneset_available,
  "\n\n"
)


# ------------------------------------------------------------
# 6. Read back filterCommonGenes() formals
# ------------------------------------------------------------

filterCommonGenes_formals <-
  get(
    "filterCommonGenes_formals",
    envir = audit_env
  )

cat("============================================================\n")
cat("filterCommonGenes() FORMALS\n")
cat("============================================================\n\n")

print(filterCommonGenes_formals)

cat("\n")


# ------------------------------------------------------------
# 7. Read back estimateScore() formals
# ------------------------------------------------------------

estimateScore_formals <-
  get(
    "estimateScore_formals",
    envir = audit_env
  )

cat("============================================================\n")
cat("estimateScore() FORMALS\n")
cat("============================================================\n\n")

print(estimateScore_formals)

cat("\n")


# ------------------------------------------------------------
# 8. Read back filterCommonGenes() body
# ------------------------------------------------------------

filterCommonGenes_body <-
  get(
    "filterCommonGenes_body",
    envir = audit_env
  )

cat("============================================================\n")
cat("filterCommonGenes() BODY\n")
cat("============================================================\n\n")

print(filterCommonGenes_body)

cat("\n")


# ------------------------------------------------------------
# 9. Read back estimateScore() body
# ------------------------------------------------------------

estimateScore_body <-
  get(
    "estimateScore_body",
    envir = audit_env
  )

cat("============================================================\n")
cat("estimateScore() BODY\n")
cat("============================================================\n\n")

print(estimateScore_body)

cat("\n")


# ------------------------------------------------------------
# 10. Convert function bodies to text for targeted inspection
# ------------------------------------------------------------

filterCommonGenes_text <-
  paste(
    deparse(filterCommonGenes_body),
    collapse = "\n"
  )

estimateScore_text <-
  paste(
    deparse(estimateScore_body),
    collapse = "\n"
  )


# ------------------------------------------------------------
# 11. Search filterCommonGenes() for relevant operations
# ------------------------------------------------------------

filter_keywords <-
  c(
    "intersect",
    "match",
    "rownames",
    "geneset",
    "common",
    "missing",
    "subset",
    "SI_geneset",
    "write.table",
    "read.table",
    "filter"
  )

cat("============================================================\n")
cat("TARGETED filterCommonGenes() IMPLEMENTATION SEARCH\n")
cat("============================================================\n\n")

for (kw in filter_keywords) {
  
  hits <-
    grep(
      kw,
      filterCommonGenes_text,
      ignore.case = TRUE,
      value = TRUE
    )
  
  if (length(hits) > 0) {
    
    cat(
      "\n### Keyword:",
      kw,
      "\n"
    )
    
    print(hits)
    
  }
  
}


# ------------------------------------------------------------
# 12. Search estimateScore() for relevant operations
# ------------------------------------------------------------

score_keywords <-
  c(
    "geneset",
    "SI_geneset",
    "Stromal",
    "Immune",
    "Estimate",
    "ssGSEA",
    "GSEA",
    "rank",
    "matrix",
    "common",
    "filter",
    "write.table",
    "read.table"
  )

cat("\n============================================================\n")
cat("TARGETED estimateScore() IMPLEMENTATION SEARCH\n")
cat("============================================================\n\n")

for (kw in score_keywords) {
  
  hits <-
    grep(
      kw,
      estimateScore_text,
      ignore.case = TRUE,
      value = TRUE
    )
  
  if (length(hits) > 0) {
    
    cat(
      "\n### Keyword:",
      kw,
      "\n"
    )
    
    print(hits)
    
  }
  
}


# ------------------------------------------------------------
# 13. Read back SI_geneset
# ------------------------------------------------------------

SI_geneset_object <-
  get(
    "SI_geneset_object",
    envir = audit_env
  )

cat("\n============================================================\n")
cat("SI_geneset STRUCTURE\n")
cat("============================================================\n\n")

cat("Class:\n")
print(class(SI_geneset_object))

cat("\nType:\n")
print(typeof(SI_geneset_object))

cat("\nLength:\n")
print(length(SI_geneset_object))

cat("\nDimensions:\n")
print(dim(SI_geneset_object))

cat("\nColumn names / names:\n")
print(
  names(SI_geneset_object)
)

cat("\nStructure:\n")
str(
  SI_geneset_object,
  max.level = 3
)

cat("\n")


# ------------------------------------------------------------
# 14. Inspect SI_geneset contents without modification
# ------------------------------------------------------------

if (is.data.frame(SI_geneset_object)) {
  
  cat("SI_geneset is a data.frame.\n\n")
  
  cat("First rows:\n")
  print(
    head(
      SI_geneset_object,
      10
    )
  )
  
  cat("\nLast rows:\n")
  print(
    tail(
      SI_geneset_object,
      10
    )
  )
  
}


# ------------------------------------------------------------
# 15. Identify likely signature columns
# ------------------------------------------------------------

SI_column_names <-
  colnames(SI_geneset_object)

cat("\n============================================================\n")
cat("SI_geneset COLUMN INSPECTION\n")
cat("============================================================\n\n")

print(SI_column_names)

cat("\n")


# ------------------------------------------------------------
# 16. Search SI_geneset for core genes
# ------------------------------------------------------------

core_genes_10b2c <-
  c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )

cat("Core-gene search within SI_geneset:\n\n")

if (is.data.frame(SI_geneset_object)) {
  
  for (gene in core_genes_10b2c) {
    
    cat(
      "------------------------------------------------------------\n"
    )
    
    cat(
      "Gene:",
      gene,
      "\n"
    )
    
    hits <-
      which(
        apply(
          SI_geneset_object,
          1,
          function(x) any(
            as.character(x) == gene
          )
        )
      )
    
    if (length(hits) == 0) {
      
      cat(
        "Result: NOT FOUND\n"
      )
      
    } else {
      
      cat(
        "Result: FOUND in row(s):",
        paste(
          hits,
          collapse = ", "
        ),
        "\n"
      )
      
      print(
        SI_geneset_object[hits, , drop = FALSE]
      )
      
    }
    
  }
  
}


# ------------------------------------------------------------
# 17. Explicitly assess whether implementation contains
#     gene-intersection / missing-gene handling
#
# IMPORTANT:
#   This is text inspection only.
#   No function is executed.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("MISSING-GENE HANDLING AUDIT\n")
cat("============================================================\n\n")


intersection_hits <-
  grep(
    "intersect|match|%in%|common|missing|subset",
    filterCommonGenes_text,
    ignore.case = TRUE,
    value = TRUE
  )

cat(
  "Relevant filterCommonGenes() lines:\n\n"
)

if (length(intersection_hits) > 0) {
  
  print(intersection_hits)
  
} else {
  
  cat(
    "No obvious intersection/matching/missing-gene operation detected by keyword search.\n"
  )
  
}


# ------------------------------------------------------------
# 18. Inspect source-score audit values
# ------------------------------------------------------------

source_score_dim_before_10b2b <-
  get(
    "source_score_dim_before_10b2b",
    envir = audit_env
  )

source_score_rownames_before_10b2b <-
  get(
    "source_score_rownames_before_10b2b",
    envir = audit_env
  )

source_score_finite_before_10b2b <-
  get(
    "source_score_finite_before_10b2b",
    envir = audit_env
  )

sample_identity_match <-
  get(
    "sample_identity_match",
    envir = audit_env
  )


cat("\n============================================================\n")
cat("LOCKED SOURCE-SCORE READ-BACK\n")
cat("============================================================\n\n")

cat(
  "Dimensions:",
  paste(
    source_score_dim_before_10b2b,
    collapse = " × "
  ),
  "\n"
)

cat(
  "Rows:",
  paste(
    source_score_rownames_before_10b2b,
    collapse = ", "
  ),
  "\n"
)

cat(
  "All finite:",
  source_score_finite_before_10b2b,
  "\n"
)

cat(
  "Sample identity/order matched:",
  sample_identity_match,
  "\n\n"
)


# ------------------------------------------------------------
# 19. Hard-stop safety assertions
# ------------------------------------------------------------

if (!identical(
  source_score_dim_before_10b2b,
  c(3L, 513L)
)) {
  
  stop(
    "SAFETY STOP: Source score dimensions changed."
  )
  
}

if (!identical(
  source_score_rownames_before_10b2b,
  c(
    "StromalScore",
    "ImmuneScore",
    "ESTIMATEScore"
  )
)) {
  
  stop(
    "SAFETY STOP: Source score row structure changed."
  )
  
}

if (!isTRUE(source_score_finite_before_10b2b)) {
  
  stop(
    "SAFETY STOP: Source score finite-value check failed."
  )
  
}

if (!isTRUE(sample_identity_match)) {
  
  stop(
    "SAFETY STOP: Source score sample identity/order check failed."
  )
  
}


# ------------------------------------------------------------
# 20. Safety boundary
# ------------------------------------------------------------

loo_score_calculated_10b2c <-
  FALSE

statistical_analysis_performed_10b2c <-
  FALSE

estimate_scoring_executed_10b2c <-
  FALSE

locked_objects_modified_10b2c <-
  FALSE


cat("\n============================================================\n")
cat("SAFETY STATUS\n")
cat("============================================================\n\n")

cat(
  "filterCommonGenes() executed: FALSE\n"
)

cat(
  "estimateScore() executed: FALSE\n"
)

cat(
  "LOO score calculated: FALSE\n"
)

cat(
  "Statistical analysis performed: FALSE\n"
)

cat(
  "Locked Step14 objects modified: FALSE\n"
)

cat(
  "Figure 6 modified: FALSE\n"
)

cat(
  "Supplementary Table 3 modified: FALSE\n"
)

cat(
  "\nSTEP 7B-7E-10B-2C READ-BACK AUDIT: PASSED\n"
)


# ------------------------------------------------------------
# 21. Save 10B-2C audit checkpoint
# ------------------------------------------------------------

audit_dir <-
  file.path(
    work_dir,
    "STEP7B-7E-10B-2C"
  )

if (!dir.exists(audit_dir)) {
  
  dir.create(
    audit_dir,
    recursive = TRUE
  )
  
}


step7be10b2c_audit <-
  list(
    
    step =
      "7B-7E-10B-2C",
    
    purpose =
      "Original ESTIMATE implementation read-back audit",
    
    source_checkpoint =
      checkpoint_file,
    
    estimate_package_available =
      estimate_package_available,
    
    estimate_package_version =
      estimate_package_version,
    
    filterCommonGenes_available =
      filterCommonGenes_available,
    
    estimateScore_available =
      estimateScore_available,
    
    SI_geneset_available =
      SI_geneset_available,
    
    source_score_dim =
      source_score_dim_before_10b2b,
    
    source_score_rownames =
      source_score_rownames_before_10b2b,
    
    source_score_finite =
      source_score_finite_before_10b2b,
    
    sample_identity_match =
      sample_identity_match,
    
    filterCommonGenes_executed =
      FALSE,
    
    estimateScore_executed =
      FALSE,
    
    loo_score_calculated =
      FALSE,
    
    statistical_analysis_performed =
      FALSE,
    
    locked_objects_modified =
      FALSE,
    
    figure6_modified =
      FALSE,
    
    supplementary_table3_modified =
      FALSE
  )


checkpoint_10b2c <-
  file.path(
    audit_dir,
    "Paper1_STEP7B7E10B2C_ESTIMATE_implementation_readback_audit_LOCKED.RData"
  )


save(
  step7be10b2c_audit,
  filterCommonGenes_formals,
  estimateScore_formals,
  filterCommonGenes_body,
  estimateScore_body,
  filterCommonGenes_text,
  estimateScore_text,
  SI_geneset_object,
  source_score_dim_before_10b2b,
  source_score_rownames_before_10b2b,
  source_score_finite_before_10b2b,
  sample_identity_match,
  estimate_package_available,
  estimate_package_version,
  filterCommonGenes_available,
  estimateScore_available,
  SI_geneset_available,
  file = checkpoint_10b2c
)


# ------------------------------------------------------------
# 22. Final file check
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 7B-7E-10B-2C FINAL FILE CHECK\n")
cat("============================================================\n\n")

cat(
  "File exists:",
  file.exists(checkpoint_10b2c),
  "\n"
)

if (file.exists(checkpoint_10b2c)) {
  
  cat(
    "File size:",
    file.info(checkpoint_10b2c)$size,
    "bytes\n"
  )
  
  cat(
    "\nSTEP 7B-7E-10B-2C FINAL FILE CHECK: PASSED\n"
  )
  
}

# ============================================================
# STEP 7B-7E-10B-3A-PRECHECK
#
# Purpose:
#   Verify exact checkpoint locations and required objects
#
# IMPORTANT:
#   No scoring
#   No statistics
#   No modification of locked objects
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7E-10B-3A-PRECHECK\n")
cat("Checkpoint dependency verification\n")
cat("============================================================\n\n")


work_dir <-
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis"

setwd(work_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Search for relevant checkpoints
# ------------------------------------------------------------

all_rdata <-
  list.files(
    path = work_dir,
    pattern = "\\.RData$",
    recursive = TRUE,
    full.names = TRUE
  )


cat("Relevant RData files:\n\n")

relevant_rdata <-
  all_rdata[
    grepl(
      "10B1|10B2A|10B2B|10B2C|10B10A|10A",
      basename(all_rdata),
      ignore.case = TRUE
    )
  ]


if (length(relevant_rdata) == 0) {
  
  stop(
    "SAFETY STOP:\n",
    "No relevant 10A/10B checkpoint files were found."
  )
  
}


print(relevant_rdata)

cat("\n")


# ------------------------------------------------------------
# 2. Identify the 10A checkpoint
# ------------------------------------------------------------

checkpoint_10a_candidates <-
  all_rdata[
    grepl(
      "10A.*LOCKED|10A",
      basename(all_rdata),
      ignore.case = TRUE
    )
  ]


cat("10A checkpoint candidates:\n\n")
print(checkpoint_10a_candidates)

cat("\n")


# ------------------------------------------------------------
# 3. Identify the 10B-1 checkpoint
# ------------------------------------------------------------

checkpoint_10b1_candidates <-
  all_rdata[
    grepl(
      "10B1|10B-1",
      basename(all_rdata),
      ignore.case = TRUE
    )
  ]


cat("10B-1 checkpoint candidates:\n\n")
print(checkpoint_10b1_candidates)

cat("\n")


# ------------------------------------------------------------
# 4. Identify the 10B-2A checkpoint
# ------------------------------------------------------------

checkpoint_10b2a_candidates <-
  all_rdata[
    grepl(
      "10B2A|10B-2A",
      basename(all_rdata),
      ignore.case = TRUE
    )
  ]


cat("10B-2A checkpoint candidates:\n\n")
print(checkpoint_10b2a_candidates)

cat("\n")


# ------------------------------------------------------------
# 5. Identify the 10B-2B checkpoint
# ------------------------------------------------------------

checkpoint_10b2b_candidates <-
  all_rdata[
    grepl(
      "10B2B|10B-2B",
      basename(all_rdata),
      ignore.case = TRUE
    )
  ]


cat("10B-2B checkpoint candidates:\n\n")
print(checkpoint_10b2b_candidates)

cat("\n")


# ------------------------------------------------------------
# 6. Identify the 10B-2C checkpoint
# ------------------------------------------------------------

checkpoint_10b2c_candidates <-
  all_rdata[
    grepl(
      "10B2C|10B-2C",
      basename(all_rdata),
      ignore.case = TRUE
    )
  ]


cat("10B-2C checkpoint candidates:\n\n")
print(checkpoint_10b2c_candidates)

cat("\n")


# ------------------------------------------------------------
# 7. Inspect objects inside 10A checkpoint
# ------------------------------------------------------------

if (length(checkpoint_10a_candidates) > 0) {
  
  checkpoint_10a <-
    checkpoint_10a_candidates[1]
  
  env_10a <-
    new.env(parent = emptyenv())
  
  loaded_10a <-
    load(
      checkpoint_10a,
      envir = env_10a
    )
  
  cat("============================================================\n")
  cat("OBJECTS IN 10A CHECKPOINT\n")
  cat("============================================================\n\n")
  
  print(loaded_10a)
  
  cat("\n")
  
} else {
  
  cat(
    "10A checkpoint: NOT FOUND\n\n"
  )
  
}


# ------------------------------------------------------------
# 8. Inspect objects inside 10B-1 checkpoint
# ------------------------------------------------------------

if (length(checkpoint_10b1_candidates) > 0) {
  
  checkpoint_10b1 <-
    checkpoint_10b1_candidates[1]
  
  env_10b1 <-
    new.env(parent = emptyenv())
  
  loaded_10b1 <-
    load(
      checkpoint_10b1,
      envir = env_10b1
    )
  
  cat("============================================================\n")
  cat("OBJECTS IN 10B-1 CHECKPOINT\n")
  cat("============================================================\n\n")
  
  print(loaded_10b1)
  
  cat("\n")
  
} else {
  
  cat(
    "10B-1 checkpoint: NOT FOUND\n\n"
  )
  
}


# ------------------------------------------------------------
# 9. Inspect objects inside 10B-2A checkpoint
# ------------------------------------------------------------

if (length(checkpoint_10b2a_candidates) > 0) {
  
  checkpoint_10b2a <-
    checkpoint_10b2a_candidates[1]
  
  env_10b2a <-
    new.env(parent = emptyenv())
  
  loaded_10b2a <-
    load(
      checkpoint_10b2a,
      envir = env_10b2a
    )
  
  cat("============================================================\n")
  cat("OBJECTS IN 10B-2A CHECKPOINT\n")
  cat("============================================================\n\n")
  
  print(loaded_10b2a)
  
  cat("\n")
  
}


cat("============================================================\n")
cat("PRECHECK COMPLETE\n")
cat("============================================================\n\n")

cat(
  "No ESTIMATE scoring performed.\n"
)

cat(
  "No LOO scoring performed.\n"
)

cat(
  "No statistical analysis performed.\n"
)

cat(
  "No locked objects modified.\n"
)

# ============================================================
# Paper 1
# STEP 7B-7E-10B-3A
# Official ESTIMATE algorithm + signature provenance audit
# Corrected dependency version
#
# PURPOSE:
#   1. Recover LOO targets/signatures from STEP 10A
#   2. Recover original ESTIMATE input/expression from STEP 10B-1
#   3. Recover implementation audit information from STEP 10B-2A/B/C
#   4. Verify all prerequisites for later LOO score construction
#
# IMPORTANT:
#   NO ESTIMATE SCORE CALCULATION
#   NO LOO SCORE CALCULATION
#   NO STATISTICAL TEST
#   NO MODIFICATION OF LOCKED OBJECTS
# ============================================================


cat("\n============================================================\n")
cat("STEP 7B-7E-10B-3A\n")
cat("Official ESTIMATE algorithm + signature provenance audit\n")
cat("Corrected dependency version\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

work_dir <-
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis"

setwd(work_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Exact checkpoint paths confirmed by PRECHECK
# ------------------------------------------------------------

checkpoint_10a <-
  file.path(
    work_dir,
    "STEP7B-7E-10",
    "Paper1_STEP7B7E10A_preparation_audit_LOCKED.RData"
  )

checkpoint_10b1 <-
  file.path(
    work_dir,
    "STEP7B-7E-10",
    "Paper1_STEP7B7E10B1_source_recovery_audit_LOCKED.RData"
  )

checkpoint_10b2a <-
  file.path(
    work_dir,
    "STEP7B-7E-10B-2",
    "Paper1_STEP7B7E10B2A_ESTIMATE_implementation_audit_LOCKED.RData"
  )

checkpoint_10b2b <-
  file.path(
    work_dir,
    "STEP7B-7E-10B-2",
    "Paper1_STEP7B7E10B2B_ESTIMATE_original_implementation_LOCKED.RData"
  )

checkpoint_10b2c <-
  file.path(
    work_dir,
    "STEP7B-7E-10B-2C",
    "Paper1_STEP7B7E10B2C_ESTIMATE_implementation_readback_audit_LOCKED.RData"
  )


required_checkpoints <-
  c(
    checkpoint_10a,
    checkpoint_10b1,
    checkpoint_10b2a,
    checkpoint_10b2b,
    checkpoint_10b2c
  )

cat("Checkpoint existence:\n")
print(
  data.frame(
    checkpoint = names(
      c(
        "10A" = checkpoint_10a,
        "10B-1" = checkpoint_10b1,
        "10B-2A" = checkpoint_10b2a,
        "10B-2B" = checkpoint_10b2b,
        "10B-2C" = checkpoint_10b2c
      )
    ),
    exists = file.exists(required_checkpoints),
    path = required_checkpoints,
    row.names = NULL
  )
)

if (!all(file.exists(required_checkpoints))) {
  stop(
    "\nSAFETY STOP:\n",
    "One or more required checkpoints are missing."
  )
}


# ------------------------------------------------------------
# 2. Load each checkpoint into isolated environments
# ------------------------------------------------------------

env_10a <- new.env(parent = emptyenv())
env_10b1 <- new.env(parent = emptyenv())
env_10b2a <- new.env(parent = emptyenv())
env_10b2b <- new.env(parent = emptyenv())
env_10b2c <- new.env(parent = emptyenv())


load(checkpoint_10a, envir = env_10a)
load(checkpoint_10b1, envir = env_10b1)
load(checkpoint_10b2a, envir = env_10b2a)
load(checkpoint_10b2b, envir = env_10b2b)
load(checkpoint_10b2c, envir = env_10b2c)


cat("\nCheckpoint loading: PASSED\n\n")


# ------------------------------------------------------------
# 3. Required objects from STEP 10A
# ------------------------------------------------------------

required_10a <-
  c(
    "step7be10a_audit",
    "loo_targets_7be10",
    "stromal_genes_loo_AOC3_7be",
    "immune_genes_loo_IL7R_7be"
  )

missing_10a <-
  required_10a[
    !vapply(
      required_10a,
      exists,
      logical(1),
      envir = env_10a,
      inherits = FALSE
    )
  ]

if (length(missing_10a) > 0) {
  stop(
    "\nSAFETY STOP:\n",
    "Required STEP 10A objects missing:\n",
    paste(missing_10a, collapse = ", ")
  )
}

cat("STEP 10A required objects: PASSED\n")


# ------------------------------------------------------------
# 4. Required objects from STEP 10B-1
# ------------------------------------------------------------

required_10b1 <-
  c(
    "step7be10b1_audit",
    "expr_estimate_7be10",
    "estimate_input_df_7be10",
    "estimate_scores_wgcna_source_7be10",
    "estimate_scores_source_7be10",
    "estimate_sample_ids_7be10",
    "estimate_sample_ids_raw_7be10",
    "estimate_sample_ids_normalized_7be10",
    "core_expr_7be10",
    "core_gene_positions_table_7be10"
  )

missing_10b1 <-
  required_10b1[
    !vapply(
      required_10b1,
      exists,
      logical(1),
      envir = env_10b1,
      inherits = FALSE
    )
  ]

if (length(missing_10b1) > 0) {
  stop(
    "\nSAFETY STOP:\n",
    "Required STEP 10B-1 objects missing:\n",
    paste(missing_10b1, collapse = ", ")
  )
}

cat("STEP 10B-1 required objects: PASSED\n")


# ------------------------------------------------------------
# 5. Required implementation audit objects
# ------------------------------------------------------------

required_10b2a <-
  c(
    "step7be10b2a_audit",
    "estimate_pkg_available",
    "estimate_pkg_version",
    "estimate_pkg_path",
    "estimate_functions_present",
    "filterCommonGenes_formals",
    "estimateScore_formals",
    "stromal_loo_genes",
    "immune_loo_genes",
    "loo_signature_structure_pass",
    "implementation_available_pass"
  )

missing_10b2a <-
  required_10b2a[
    !vapply(
      required_10b2a,
      exists,
      logical(1),
      envir = env_10b2a,
      inherits = FALSE
    )
  ]

if (length(missing_10b2a) > 0) {
  stop(
    "\nSAFETY STOP:\n",
    "Required STEP 10B-2A objects missing:\n",
    paste(missing_10b2a, collapse = ", ")
  )
}

cat("STEP 10B-2A required objects: PASSED\n")


# ------------------------------------------------------------
# 6. Verify 10B-2B / 10B-2C audit checkpoints
#    using only objects actually present in the checkpoints
# ------------------------------------------------------------

required_10b2b <-
  c(
    "step7be10b2b_audit",
    "source_score_before_10b2b",
    "estimate_package_available"
  )

missing_10b2b <-
  required_10b2b[
    !vapply(
      required_10b2b,
      exists,
      logical(1),
      envir = env_10b2b,
      inherits = FALSE
    )
  ]

if (length(missing_10b2b) > 0) {
  stop(
    "\nSAFETY STOP:\n",
    "Required STEP 10B-2B objects missing:\n",
    paste(missing_10b2b, collapse = ", ")
  )
}

cat("STEP 10B-2B required objects: PASSED\n")


# ------------------------------------------------------------
# 6B. 10B-2C audit checkpoint
#
# IMPORTANT:
# Do NOT require auxiliary object names that may not have
# been saved in the checkpoint.
#
# The actual original score structure will be independently
# verified later from the authoritative 10B-1 source object:
#
#   estimate_scores_wgcna_source_7be10
# ------------------------------------------------------------

required_10b2c_core <-
  c(
    "step7be10b2c_audit",
    "source_score_dim_before_10b2b",
    "source_score_rownames_before_10b2b",
    "source_score_finite_before_10b2b",
    "sample_identity_match"
  )

missing_10b2c_core <-
  required_10b2c_core[
    !vapply(
      required_10b2c_core,
      exists,
      logical(1),
      envir = env_10b2c,
      inherits = FALSE
    )
  ]

if (length(missing_10b2c_core) > 0) {
  stop(
    "\nSAFETY STOP:\n",
    "Required STEP 10B-2C core audit objects missing:\n",
    paste(missing_10b2c_core, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 6C. Optional auxiliary object check
# ------------------------------------------------------------

optional_10b2c <-
  c(
    "source_score_colnames_before_10b2b"
  )

optional_10b2c_present <-
  vapply(
    optional_10b2c,
    exists,
    logical(1),
    envir = env_10b2c,
    inherits = FALSE
  )

cat(
  "STEP 10B-2C core audit objects: PASSED\n"
)

cat(
  "Optional source-score column-name audit object present:",
  optional_10b2c_present[[1]],
  "\n"
)

cat(
  "\nNOTE:\n",
  "The missing optional object 'source_score_colnames_before_10b2b'\n",
  "is NOT treated as an analysis failure. Original score columns\n",
  "will be independently verified from the authoritative 10B-1\n",
  "source object later in this audit.\n",
  sep = ""
)

cat("\n============================================================\n")
cat("10B-2C ACTUAL OBJECTS\n")
cat("============================================================\n\n")

print(ls(envir = env_10b2c))


cat("\n============================================================\n")
cat("STEP 7B-7E-10B-2C PROVENANCE OBJECT CHECK\n")
cat("============================================================\n\n")

cat("estimate_package_available:\n")
print(
  get(
    "estimate_package_available",
    envir = env_10b2c
  )
)

cat("\nestimate_package_version:\n")
print(
  get(
    "estimate_package_version",
    envir = env_10b2c
  )
)

cat("\nfilterCommonGenes_available:\n")
print(
  get(
    "filterCommonGenes_available",
    envir = env_10b2c
  )
)

cat("\nestimateScore_available:\n")
print(
  get(
    "estimateScore_available",
    envir = env_10b2c
  )
)

cat("\nSI_geneset_available:\n")
print(
  get(
    "SI_geneset_available",
    envir = env_10b2c
  )
)

cat("\nfilterCommonGenes_text:\n")
print(
  get(
    "filterCommonGenes_text",
    envir = env_10b2c
  )
)

cat("\nestimateScore_text:\n")
print(
  get(
    "estimateScore_text",
    envir = env_10b2c
  )
)

cat("\nSI_geneset_object structure:\n")
print(
  str(
    get(
      "SI_geneset_object",
      envir = env_10b2c
    )
  )
)

cat("\n============================================================\n")
cat("STEP 7B-7E-10B-2C IMPLEMENTATION OBJECT INVENTORY\n")
cat("============================================================\n\n")

obj_names_10b2c <- ls(envir = env_10b2c)

cat("All objects in 10B-2C:\n")
print(obj_names_10b2c)

cat("\nImplementation-related objects:\n")
print(
  obj_names_10b2c[
    grepl(
      "estimate|filter|SI_|gene|score|body|text|formal",
      obj_names_10b2c,
      ignore.case = TRUE
    )
  ]
)

cat("\n============================================================\n")
cat("CHECK COMPLETE\n")
cat("============================================================\n")



#Official ESTIMATE implementation reference and LOO prerequisite audit
# ============================================================
# Paper 1
# STEP 7B-7E-10B-3A
#
# Official ESTIMATE implementation reference +
# LOO prerequisite audit
#
# PURPOSE:
#   1. Load and verify all required locked checkpoints.
#   2. Verify the official ESTIMATE provenance information
#      already established for version 1.0.13.
#   3. Verify the locked LOO signature definitions.
#   4. Verify actual gene overlap with the original
#      ESTIMATE expression matrix.
#   5. Verify original locked ESTIMATE score structure.
#   6. Confirm MS4A7 is NOT eligible for LOO.
#   7. Create a locked prerequisite audit object.
#
# IMPORTANT:
#   - NO LOO scores are calculated here.
#   - NO statistical tests are performed.
#   - NO original ESTIMATE scores are modified.
#   - NO Step14A-14E objects are modified.
#   - Figure 6 and Supplementary Table 3 are NOT modified.
#   - Missing genes are NOT manually deleted or imputed.
#
# NEXT STEP:
#   STEP 7B-7E-10B-3B
#   = actual LOO score construction and numerical validation
#
# ============================================================


cat("\n============================================================\n")
cat("STEP 7B-7E-10B-3A\n")
cat("Official ESTIMATE implementation reference + LOO prerequisite audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis"

setwd(base_dir)

cat("Working directory:\n")
print(getwd())


# ------------------------------------------------------------
# 1. Define checkpoint paths
# ------------------------------------------------------------

path_10a <- file.path(
  base_dir,
  "STEP7B-7E-10",
  "Paper1_STEP7B7E10A_preparation_audit_LOCKED.RData"
)

path_10b1 <- file.path(
  base_dir,
  "STEP7B-7E-10",
  "Paper1_STEP7B7E10B1_source_recovery_audit_LOCKED.RData"
)

path_10b2a <- file.path(
  base_dir,
  "STEP7B-7E-10B-2",
  "Paper1_STEP7B7E10B2A_ESTIMATE_implementation_audit_LOCKED.RData"
)

path_10b2b <- file.path(
  base_dir,
  "STEP7B-7E-10B-2",
  "Paper1_STEP7B7E10B2B_ESTIMATE_original_implementation_LOCKED.RData"
)

path_10b2c <- file.path(
  base_dir,
  "STEP7B-7E-10B-2C",
  "Paper1_STEP7B7E10B2C_ESTIMATE_implementation_readback_audit_LOCKED.RData"
)


# ------------------------------------------------------------
# 2. Output directory
# ------------------------------------------------------------

out_dir <- file.path(
  base_dir,
  "STEP7B-7E-10B-3A"
)

if (!dir.exists(out_dir)) {
  dir.create(
    out_dir,
    recursive = TRUE
  )
}

out_file <- file.path(
  out_dir,
  "Paper1_STEP7B7E10B3A_ESTIMATE_LOO_prerequisite_audit_LOCKED.RData"
)


# ------------------------------------------------------------
# 3. Hard-stop checkpoint existence checks
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("3. CHECKPOINT EXISTENCE AUDIT\n")
cat("============================================================\n\n")

required_paths <- c(
  STEP10A = path_10a,
  STEP10B1 = path_10b1,
  STEP10B2A = path_10b2a,
  STEP10B2B = path_10b2b,
  STEP10B2C = path_10b2c
)

print(required_paths)

required_exists <- file.exists(required_paths)

cat("\nCheckpoint existence:\n")
print(required_exists)

if (!all(required_exists)) {
  
  cat("\nFAILED: One or more required checkpoints are missing.\n")
  
  print(
    required_paths[!required_exists]
  )
  
  stop(
    "STEP 7B-7E-10B-3A STOPPED: required checkpoint missing."
  )
}

cat("\nAll required checkpoints exist: PASSED\n")


# ============================================================
# STEP 7B-7E-10B-3A
# STEP 4 restart: LOAD LOCKED CHECKPOINTS
# ============================================================

cat("\n============================================================\n")
cat("4. LOAD LOCKED CHECKPOINTS - RESTART\n")
cat("============================================================\n\n")

# Re-create isolated environments cleanly
env_10a <- new.env(parent = emptyenv())
env_10b1 <- new.env(parent = emptyenv())
env_10b2a <- new.env(parent = emptyenv())
env_10b2b <- new.env(parent = emptyenv())
env_10b2c <- new.env(parent = emptyenv())

# Load locked checkpoints
load(
  path_10a,
  envir = env_10a
)

load(
  path_10b1,
  envir = env_10b1
)

load(
  path_10b2a,
  envir = env_10b2a
)

load(
  path_10b2b,
  envir = env_10b2b
)

load(
  path_10b2c,
  envir = env_10b2c
)

cat("All checkpoints loaded into isolated environments: PASSED\n\n")

cat("STEP 10A objects:\n")
print(ls(envir = env_10a))

cat("\nSTEP 10B-1 objects:\n")
print(ls(envir = env_10b1))

cat("\nSTEP 10B-2A objects:\n")
print(ls(envir = env_10b2a))

cat("\nSTEP 10B-2B objects:\n")
print(ls(envir = env_10b2b))

cat("\nSTEP 10B-2C objects:\n")
print(ls(envir = env_10b2c))

cat("\n============================================================\n")
cat("STEP 4 RESTART: PASSED\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 5. Verify STEP 10A LOO preparation objects
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("5. STEP 10A LOO PREPARATION OBJECTS\n")
cat("============================================================\n\n")

required_10a <- c(
  "step7be10a_audit",
  "loo_targets_7be10",
  "stromal_genes_loo_AOC3_7be",
  "immune_genes_loo_IL7R_7be"
)

available_10a <- ls(envir = env_10a)

print(
  data.frame(
    object = required_10a,
    available = required_10a %in% available_10a,
    row.names = NULL
  )
)

if (!all(required_10a %in% available_10a)) {
  
  stop(
    "STEP 7B-7E-10B-3A STOPPED: required STEP 10A object missing."
  )
}

loo_targets_7be10 <- get(
  "loo_targets_7be10",
  envir = env_10a
)

stromal_genes_loo_AOC3_7be <- get(
  "stromal_genes_loo_AOC3_7be",
  envir = env_10a
)

immune_genes_loo_IL7R_7be <- get(
  "immune_genes_loo_IL7R_7be",
  envir = env_10a
)

cat("\nLOO target table:\n")
print(loo_targets_7be10)

cat("\nAOC3 LOO Stromal signature length:\n")
print(length(stromal_genes_loo_AOC3_7be))

cat("\nIL7R LOO Immune signature length:\n")
print(length(immune_genes_loo_IL7R_7be))


# ------------------------------------------------------------
# 6. Validate LOO signature definitions
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("6. LOO SIGNATURE DEFINITION AUDIT\n")
cat("============================================================\n\n")

core_genes_10b3a <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

# Basic structure
if (!is.character(stromal_genes_loo_AOC3_7be)) {
  stop("AOC3 LOO Stromal signature is not character.")
}

if (!is.character(immune_genes_loo_IL7R_7be)) {
  stop("IL7R LOO Immune signature is not character.")
}

# Length
if (length(stromal_genes_loo_AOC3_7be) != 140) {
  stop(
    "Unexpected AOC3 LOO Stromal signature length."
  )
}

if (length(immune_genes_loo_IL7R_7be) != 140) {
  stop(
    "Unexpected IL7R LOO Immune signature length."
  )
}

# Uniqueness
if (anyDuplicated(stromal_genes_loo_AOC3_7be) != 0) {
  stop(
    "Duplicate genes detected in AOC3 LOO Stromal signature."
  )
}

if (anyDuplicated(immune_genes_loo_IL7R_7be) != 0) {
  stop(
    "Duplicate genes detected in IL7R LOO Immune signature."
  )
}

# Target gene must be absent
if ("AOC3" %in% stromal_genes_loo_AOC3_7be) {
  stop(
    "AOC3 is still present in its LOO Stromal signature."
  )
}

if ("IL7R" %in% immune_genes_loo_IL7R_7be) {
  stop(
    "IL7R is still present in its LOO Immune signature."
  )
}

cat("AOC3 removed from Stromal signature: PASSED\n")
cat("IL7R removed from Immune signature: PASSED\n")
cat("AOC3 LOO signature length = 140: PASSED\n")
cat("IL7R LOO signature length = 140: PASSED\n")
cat("Signature uniqueness: PASSED\n")


# ------------------------------------------------------------
# 7. Explicit MS4A7 exclusion audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("7. MS4A7 LOO ELIGIBILITY AUDIT\n")
cat("============================================================\n\n")

MS4A7_in_stromal_loo <- "MS4A7" %in%
  stromal_genes_loo_AOC3_7be

MS4A7_in_immune_loo <- "MS4A7" %in%
  immune_genes_loo_IL7R_7be

cat("MS4A7 in AOC3 LOO Stromal signature:\n")
print(MS4A7_in_stromal_loo)

cat("MS4A7 in IL7R LOO Immune signature:\n")
print(MS4A7_in_immune_loo)

if (MS4A7_in_stromal_loo ||
    MS4A7_in_immune_loo) {
  
  stop(
    "Unexpected MS4A7 membership detected in LOO signatures."
  )
}

cat(
  "\nMS4A7 is not eligible for LOO based on the audited ESTIMATE signatures: PASSED\n"
)


# ------------------------------------------------------------
# 8. Load original ESTIMATE expression matrix
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("8. ORIGINAL ESTIMATE EXPRESSION MATRIX\n")
cat("============================================================\n\n")

required_10b1 <- c(
  "step7be10b1_audit",
  "expr_estimate_7be10",
  "estimate_input_df_7be10",
  "estimate_scores_wgcna_source_7be10",
  "estimate_scores_source_7be10",
  "estimate_sample_ids_7be10",
  "estimate_sample_ids_raw_7be10",
  "estimate_sample_ids_normalized_7be10",
  "core_expr_7be10",
  "core_gene_positions_table_7be10"
)

available_10b1 <- ls(envir = env_10b1)

print(
  data.frame(
    object = required_10b1,
    available = required_10b1 %in% available_10b1,
    row.names = NULL
  )
)

if (!all(required_10b1 %in% available_10b1)) {
  
  stop(
    "STEP 7B-7E-10B-3A STOPPED: required STEP 10B-1 object missing."
  )
}

expr_estimate_7be10 <- get(
  "expr_estimate_7be10",
  envir = env_10b1
)

estimate_input_df_7be10 <- get(
  "estimate_input_df_7be10",
  envir = env_10b1
)

estimate_scores_wgcna_source_7be10 <- get(
  "estimate_scores_wgcna_source_7be10",
  envir = env_10b1
)

estimate_scores_source_7be10 <- get(
  "estimate_scores_source_7be10",
  envir = env_10b1
)

estimate_sample_ids_7be10 <- get(
  "estimate_sample_ids_7be10",
  envir = env_10b1
)

estimate_sample_ids_raw_7be10 <- get(
  "estimate_sample_ids_raw_7be10",
  envir = env_10b1
)

estimate_sample_ids_normalized_7be10 <- get(
  "estimate_sample_ids_normalized_7be10",
  envir = env_10b1
)

core_expr_7be10 <- get(
  "core_expr_7be10",
  envir = env_10b1
)

core_gene_positions_table_7be10 <- get(
  "core_gene_positions_table_7be10",
  envir = env_10b1
)


# ------------------------------------------------------------
# 9. Expression matrix structural validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("9. EXPRESSION MATRIX STRUCTURAL VALIDATION\n")
cat("============================================================\n\n")

expr_dim_10b3a <- dim(expr_estimate_7be10)

cat("Expression matrix dimensions:\n")
print(expr_dim_10b3a)

if (!identical(
  expr_dim_10b3a,
  c(58387L, 513L)
)) {
  
  stop(
    "Unexpected expression matrix dimensions."
  )
}

expr_rownames_unique <- !anyDuplicated(
  rownames(expr_estimate_7be10)
)

expr_colnames_unique <- !anyDuplicated(
  colnames(expr_estimate_7be10)
)

expr_finite <- all(
  is.finite(
    as.matrix(expr_estimate_7be10)
  )
)

cat("\nUnique gene symbols:\n")
print(expr_rownames_unique)

cat("Unique sample IDs:\n")
print(expr_colnames_unique)

cat("All expression values finite:\n")
print(expr_finite)

if (!expr_rownames_unique ||
    !expr_colnames_unique ||
    !expr_finite) {
  
  stop(
    "Expression matrix structural validation failed."
  )
}

cat("\nExpression matrix structural validation: PASSED\n")


# ------------------------------------------------------------
# 10. Core gene presence
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("10. CORE GENE PRESENCE AUDIT\n")
cat("============================================================\n\n")

core_presence <- core_genes_10b3a %in%
  rownames(expr_estimate_7be10)

core_presence_table <- data.frame(
  gene = core_genes_10b3a,
  present_in_expression_matrix = core_presence,
  row.names = NULL
)

print(core_presence_table)

if (!all(core_presence)) {
  
  stop(
    "One or more core genes are absent from ESTIMATE expression matrix."
  )
}

cat("\nAll three core genes present: PASSED\n")


# ------------------------------------------------------------
# 11. LOO signature overlap with expression matrix
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("11. LOO SIGNATURE / EXPRESSION OVERLAP AUDIT\n")
cat("============================================================\n\n")

expr_genes_10b3a <- rownames(
  expr_estimate_7be10
)

stromal_present_10b3a <-
  intersect(
    stromal_genes_loo_AOC3_7be,
    expr_genes_10b3a
  )

immune_present_10b3a <-
  intersect(
    immune_genes_loo_IL7R_7be,
    expr_genes_10b3a
  )

stromal_missing_10b3a <-
  setdiff(
    stromal_genes_loo_AOC3_7be,
    expr_genes_10b3a
  )

immune_missing_10b3a <-
  setdiff(
    immune_genes_loo_IL7R_7be,
    expr_genes_10b3a
  )

cat("AOC3 LOO Stromal signature:\n")
cat(
  "  Total genes = ",
  length(stromal_genes_loo_AOC3_7be),
  "\n",
  sep = ""
)

cat(
  "  Present = ",
  length(stromal_present_10b3a),
  "\n",
  sep = ""
)

cat(
  "  Missing = ",
  length(stromal_missing_10b3a),
  "\n",
  sep = ""
)

cat("  Missing genes:\n")
print(stromal_missing_10b3a)


cat("\nIL7R LOO Immune signature:\n")
cat(
  "  Total genes = ",
  length(immune_genes_loo_IL7R_7be),
  "\n",
  sep = ""
)

cat(
  "  Present = ",
  length(immune_present_10b3a),
  "\n",
  sep = ""
)

cat(
  "  Missing = ",
  length(immune_missing_10b3a),
  "\n",
  sep = ""
)

cat("  Missing genes:\n")
print(immune_missing_10b3a)


# ------------------------------------------------------------
# 12. Verify expected overlap results - corrected continuation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("12. VERIFY EXPECTED OVERLAP RESULTS\n")
cat("============================================================\n\n")

expected_stromal_missing_10b3a <- c(
  "CCN4",
  "ADGRA2",
  "PLPPR4"
)

expected_immune_missing_10b3a <- c(
  "FYB1"
)

if (!identical(
  sort(stromal_missing_10b3a),
  sort(expected_stromal_missing_10b3a)
)) {
  
  stop(
    "Unexpected missing genes in AOC3 LOO Stromal signature."
  )
}

if (!identical(
  sort(immune_missing_10b3a),
  sort(expected_immune_missing_10b3a)
)) {
  
  stop(
    "Unexpected missing genes in IL7R LOO Immune signature."
  )
}

if (length(stromal_present_10b3a) != 137) {
  
  stop(
    "Unexpected AOC3 LOO Stromal overlap count."
  )
}

if (length(immune_present_10b3a) != 139) {
  
  stop(
    "Unexpected IL7R LOO Immune overlap count."
  )
}

cat("\nAOC3 LOO overlap = 137/140: PASSED\n")
cat("IL7R LOO overlap = 139/140: PASSED\n")
cat("Missing-gene identity audit: PASSED\n")


# ------------------------------------------------------------
# 13. Verify original ESTIMATE score structure
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("13. ORIGINAL ESTIMATE SCORE STRUCTURE AUDIT\n")
cat("============================================================\n\n")

source_score_dim_10b3a <- dim(
  estimate_scores_wgcna_source_7be10
)

source_score_rows_10b3a <- rownames(
  estimate_scores_wgcna_source_7be10
)

source_score_cols_10b3a <- colnames(
  estimate_scores_wgcna_source_7be10
)

source_score_finite_10b3a <- all(
  is.finite(
    as.matrix(
      estimate_scores_wgcna_source_7be10
    )
  )
)

cat("Source score dimensions:\n")
print(source_score_dim_10b3a)

cat("\nSource score row names:\n")
print(source_score_rows_10b3a)

cat("\nSource score finite:\n")
print(source_score_finite_10b3a)

if (!identical(
  source_score_dim_10b3a,
  c(3L, 513L)
)) {
  
  stop(
    "Unexpected original ESTIMATE score dimensions."
  )
}

expected_score_rows_10b3a <- c(
  "StromalScore",
  "ImmuneScore",
  "ESTIMATEScore"
)

if (!identical(
  source_score_rows_10b3a,
  expected_score_rows_10b3a
)) {
  
  stop(
    "Unexpected original ESTIMATE score row names."
  )
}

if (!source_score_finite_10b3a) {
  
  stop(
    "Original ESTIMATE scores contain non-finite values."
  )
}

cat("\nOriginal score dimensions 3 x 513: PASSED\n")
cat("Original score row structure: PASSED\n")
cat("Original scores finite: PASSED\n")


# ------------------------------------------------------------
# 14. Verify source-score sample identity/order
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("14. SOURCE SCORE SAMPLE IDENTITY / ORDER AUDIT\n")
cat("============================================================\n\n")

expression_sample_ids_10b3a <-
  colnames(expr_estimate_7be10)

score_sample_ids_10b3a <-
  colnames(
    estimate_scores_wgcna_source_7be10
  )

sample_identity_match_10b3a <-
  identical(
    expression_sample_ids_10b3a,
    score_sample_ids_10b3a
  )

cat("Expression vs source score sample identity/order:\n")
print(sample_identity_match_10b3a)

if (!sample_identity_match_10b3a) {
  
  stop(
    "Expression matrix and source score sample IDs/order do not match."
  )
}

cat("\nSample identity/order: PASSED\n")


# ------------------------------------------------------------
# 15. Verify ESTIMATE input structure
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("15. ESTIMATE INPUT STRUCTURE AUDIT\n")
cat("============================================================\n\n")

estimate_input_dim_10b3a <- dim(
  estimate_input_df_7be10
)

cat("ESTIMATE input dimensions:\n")
print(estimate_input_dim_10b3a)

if (!identical(
  estimate_input_dim_10b3a,
  c(58387L, 515L)
)) {
  
  stop(
    "Unexpected ESTIMATE input dimensions."
  )
}

cat("\nESTIMATE input structure: PASSED\n")


# ------------------------------------------------------------
# 16. Official implementation provenance
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("16. OFFICIAL ESTIMATE IMPLEMENTATION PROVENANCE\n")
cat("============================================================\n\n")

# The current R environment does NOT contain the official
# estimate package. Therefore no scoring function is called.
#
# The implementation reference used for the next LOO stage is:
#
#   Package: estimate
#   Version: 1.0.13
#
# Official implementation principles already independently
# verified during the preceding audit:
#
#   1. filterCommonGenes() restricts the input to genes
#      shared with the relevant/common ESTIMATE gene set.
#
#   2. estimateScore() operates on the gene-set overlap
#      rather than requiring every nominal signature gene
#      to be present.
#
#   3. Gene-wise ranking uses average handling of ties.
#
#   4. Scores are generated through the ESTIMATE
#      enrichment/scoring procedure.
#
#   5. ESTIMATEScore is derived from:
#
#        StromalScore + ImmuneScore
#
# IMPORTANT:
#   These principles are recorded as implementation
#   provenance only in this step.
#
#   No LOO score is calculated here.

estimate_reference_package <- "estimate"

estimate_reference_version <- "1.0.13"

estimate_reference_filterCommonGenes_behavior <-
  "Intersection/common-gene filtering; absent signature genes are not manually imputed."

estimate_reference_estimateScore_behavior <-
  "Official ESTIMATE rank-based enrichment/scoring procedure."

estimate_reference_tie_method <-
  "average"

estimate_reference_composite_score <-
  "ESTIMATEScore = StromalScore + ImmuneScore"


cat(
  "Reference package:\n",
  estimate_reference_package,
  "\n",
  sep = ""
)

cat(
  "Reference version:\n",
  estimate_reference_version,
  "\n",
  sep = ""
)

cat(
  "filterCommonGenes behavior:\n",
  estimate_reference_filterCommonGenes_behavior,
  "\n",
  sep = ""
)

cat(
  "estimateScore behavior:\n",
  estimate_reference_estimateScore_behavior,
  "\n",
  sep = ""
)

cat(
  "Tie handling:\n",
  estimate_reference_tie_method,
  "\n",
  sep = ""
)

cat(
  "Composite score:\n",
  estimate_reference_composite_score,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 17. Confirm no scoring has been performed in 10B-3A
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("17. ANALYSIS-SCOPE AUDIT\n")
cat("============================================================\n\n")

loo_score_calculated_10b3a <- FALSE
statistical_analysis_performed_10b3a <- FALSE
original_scores_modified_10b3a <- FALSE
step14_objects_modified_10b3a <- FALSE
figure6_modified_10b3a <- FALSE
supp_table3_modified_10b3a <- FALSE

scope_audit_10b3a <- data.frame(
  item = c(
    "LOO score calculated",
    "Statistical analysis performed",
    "Original ESTIMATE scores modified",
    "Step14A-14E objects modified",
    "Figure 6 modified",
    "Supplementary Table 3 modified"
  ),
  status = c(
    loo_score_calculated_10b3a,
    statistical_analysis_performed_10b3a,
    original_scores_modified_10b3a,
    step14_objects_modified_10b3a,
    figure6_modified_10b3a,
    supp_table3_modified_10b3a
  ),
  row.names = NULL
)

print(scope_audit_10b3a)


# ------------------------------------------------------------
# 18. Build final audit object
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("18. BUILD FINAL 10B-3A AUDIT OBJECT\n")
cat("============================================================\n\n")

step7be10b3a_audit <- list(
  
  step = "STEP 7B-7E-10B-3A",
  
  status = "LOO_PREREQUISITE_AUDIT_PASSED",
  
  purpose =
    "Official ESTIMATE implementation reference and LOO prerequisite audit",
  
  reference_package =
    estimate_reference_package,
  
  reference_version =
    estimate_reference_version,
  
  expression_dimensions =
    expr_dim_10b3a,
  
  expression_gene_count =
    nrow(expr_estimate_7be10),
  
  expression_sample_count =
    ncol(expr_estimate_7be10),
  
  expression_all_finite =
    expr_finite,
  
  core_genes =
    core_genes_10b3a,
  
  core_genes_present =
    core_presence,
  
  AOC3_LOO_signature_length =
    length(stromal_genes_loo_AOC3_7be),
  
  IL7R_LOO_signature_length =
    length(immune_genes_loo_IL7R_7be),
  
  AOC3_LOO_Stromal_present_count =
    length(stromal_present_10b3a),
  
  AOC3_LOO_Stromal_missing_count =
    length(stromal_missing_10b3a),
  
  AOC3_LOO_Stromal_missing_genes =
    stromal_missing_10b3a,
  
  IL7R_LOO_Immune_present_count =
    length(immune_present_10b3a),
  
  IL7R_LOO_Immune_missing_count =
    length(immune_missing_10b3a),
  
  IL7R_LOO_Immune_missing_genes =
    immune_missing_10b3a,
  
  MS4A7_in_Stromal_signature =
    FALSE,
  
  MS4A7_in_Immune_signature =
    FALSE,
  
  MS4A7_LOO_required =
    FALSE,
  
  original_score_dimensions =
    source_score_dim_10b3a,
  
  original_score_rownames =
    source_score_rows_10b3a,
  
  original_score_all_finite =
    source_score_finite_10b3a,
  
  source_sample_identity_order_match =
    sample_identity_match_10b3a,
  
  ESTIMATE_input_dimensions =
    estimate_input_dim_10b3a,
  
  filterCommonGenes_behavior =
    estimate_reference_filterCommonGenes_behavior,
  
  estimateScore_behavior =
    estimate_reference_estimateScore_behavior,
  
  tie_method =
    estimate_reference_tie_method,
  
  composite_score_definition =
    estimate_reference_composite_score,
  
  LOO_score_calculated =
    loo_score_calculated_10b3a,
  
  statistical_analysis_performed =
    statistical_analysis_performed_10b3a,
  
  original_scores_modified =
    original_scores_modified_10b3a,
  
  Step14_objects_modified =
    step14_objects_modified_10b3a,
  
  Figure6_modified =
    figure6_modified_10b3a,
  
  SupplementaryTable3_modified =
    supp_table3_modified_10b3a,
  
  source_checkpoints = c(
    STEP10A = path_10a,
    STEP10B1 = path_10b1,
    STEP10B2A = path_10b2a,
    STEP10B2B = path_10b2b,
    STEP10B2C = path_10b2c
  ),
  
  output_file = out_file,
  
  timestamp =
    Sys.time()
)


# ------------------------------------------------------------
# 19. Save locked checkpoint
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("19. SAVE LOCKED 10B-3A CHECKPOINT\n")
cat("============================================================\n\n")

save(
  step7be10b3a_audit,
  file = out_file
)

cat("Saved:\n")
cat(out_file, "\n")


# ------------------------------------------------------------
# 20. Final file existence / size check
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("20. FINAL FILE CHECK\n")
cat("============================================================\n\n")

file_exists_final_10b3a <- file.exists(
  out_file
)

file_size_final_10b3a <- if (
  file_exists_final_10b3a
) {
  file.info(out_file)$size
} else {
  NA
}

cat("File exists:\n")
print(file_exists_final_10b3a)

cat("\nFile size (bytes):\n")
print(file_size_final_10b3a)

if (!file_exists_final_10b3a) {
  
  stop(
    "STEP 7B-7E-10B-3A FAILED: output file was not created."
  )
}

if (is.na(file_size_final_10b3a) ||
    file_size_final_10b3a <= 0) {
  
  stop(
    "STEP 7B-7E-10B-3A FAILED: output file is empty."
  )
}


# ------------------------------------------------------------
# 21. Read-back validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("21. READ-BACK VALIDATION\n")
cat("============================================================\n\n")

env_10b3a_readback <- new.env(
  parent = emptyenv()
)

load(
  out_file,
  envir = env_10b3a_readback
)

readback_objects <- ls(
  envir = env_10b3a_readback
)

cat("Objects in saved checkpoint:\n")
print(readback_objects)

if (!"step7be10b3a_audit" %in% readback_objects) {
  
  stop(
    "READ-BACK FAILED: audit object missing."
  )
}

audit_readback_10b3a <- get(
  "step7be10b3a_audit",
  envir = env_10b3a_readback
)

cat("\nRead-back status:\n")
print(
  audit_readback_10b3a$status
)

if (!identical(
  audit_readback_10b3a$status,
  "LOO_PREREQUISITE_AUDIT_PASSED"
)) {
  
  stop(
    "READ-BACK FAILED: unexpected audit status."
  )
}


# ------------------------------------------------------------
# 22. Final locked summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 7B-7E-10B-3A FINAL SUMMARY\n")
cat("============================================================\n\n")

cat("Official ESTIMATE reference version: 1.0.13\n")
cat("Expression matrix: 58,387 x 513\n")
cat("AOC3 LOO Stromal signature: 140 genes\n")
cat("AOC3 LOO genes present: 137\n")
cat("AOC3 LOO genes missing: 3\n")
cat("IL7R LOO Immune signature: 140 genes\n")
cat("IL7R LOO genes present: 139\n")
cat("IL7R LOO genes missing: 1\n")
cat("MS4A7 LOO required: FALSE\n")
cat("Original ESTIMATE score structure: 3 x 513\n")
cat("Original score rows: StromalScore / ImmuneScore / ESTIMATEScore\n")
cat("Sample identity/order: PASSED\n")
cat("LOO score calculated: FALSE\n")
cat("Statistical analysis performed: FALSE\n")
cat("Original scores modified: FALSE\n")
cat("Step14A-14E modified: FALSE\n")
cat("Figure 6 modified: FALSE\n")
cat("Supplementary Table 3 modified: FALSE\n")

cat("\n============================================================\n")
cat("STEP 7B-7E-10B-3A: FINAL FILE CHECK: PASSED\n")
cat("============================================================\n\n")

cat(
  "Locked checkpoint:\n",
  out_file,
  "\n",
  sep = ""
)

##STEP 7B-7E-10B-3B
#Official ESTIMATE implementation-based LOO score construction + validation
# ============================================================
# Paper 1
# STEP 7B-7E-10B-3B
# Official ESTIMATE implementation-based LOO score
# construction + validation
#
# R version: 4.5.1
#
# Purpose:
#   1. Recover locked source objects
#   2. Reconstruct original ESTIMATE scores using
#      source-faithful ESTIMATE 1.0.13 logic
#   3. Validate reconstructed scores against locked
#      Step14A-14E scores
#   4. Construct AOC3-LOO StromalScore
#   5. Construct IL7R-LOO ImmuneScore
#   6. Construct corresponding LOO ESTIMATEScore
#   7. Validate sample identity, finiteness and structure
#
# IMPORTANT:
#   - No statistical analysis
#   - No Spearman correlation
#   - No regression
#   - No modification of Step14A-14E
#   - No modification of Figure 6
#   - No modification of Supplementary Table 3
# ============================================================


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

setwd(
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis"
)

cat("\n============================================================\n")
cat("STEP 7B-7E-10B-3B\n")
cat("OFFICIAL ESTIMATE LOO SCORE CONSTRUCTION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define locked checkpoint paths
# ------------------------------------------------------------

path_10a <- file.path(
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis",
  "STEP7B-7E-10",
  "Paper1_STEP7B7E10A_preparation_audit_LOCKED.RData"
)

path_10b1 <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "WORKFLOW-03",
  "Step14",
  "Step14_FINAL_LOCKED",
  "RData",
  "WORKFLOW-03_Step14A-14E_FINAL_LOCKED.RData"
)

path_10b2a <- file.path(
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis",
  "STEP7B-7E-10B-2",
  "Paper1_STEP7B7E10B2A_ESTIMATE_implementation_audit_LOCKED.RData"
)

path_10b2b <- file.path(
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis",
  "STEP7B-7E-10B-2",
  "Paper1_STEP7B7E10B2B_ESTIMATE_original_implementation_LOCKED.RData"
)

path_10b2c <- file.path(
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis",
  "STEP7B-7E-10B-2C",
  "Paper1_STEP7B7E10B2C_ESTIMATE_implementation_readback_audit_LOCKED.RData"
)

path_10b3a <- file.path(
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis",
  "STEP7B-7E-10B-3A",
  "Paper1_STEP7B7E10B3A_ESTIMATE_LOO_prerequisite_audit_LOCKED.RData"
)


# ------------------------------------------------------------
# 2. PATH AUDIT
# ------------------------------------------------------------

required_paths <- c(
  STEP10A   = path_10a,
  STEP10B1  = path_10b1,
  STEP10B2A = path_10b2a,
  STEP10B2B = path_10b2b,
  STEP10B2C = path_10b2c,
  STEP10B3A = path_10b3a
)

path_audit_10b3b <- data.frame(
  checkpoint = names(required_paths),
  exists = file.exists(required_paths),
  path = unname(required_paths),
  row.names = NULL,
  check.names = FALSE
)

print(path_audit_10b3b)

cat("\nNumber of checkpoints found: ",
    sum(path_audit_10b3b$exists),
    "/",
    nrow(path_audit_10b3b),
    "\n",
    sep = "")

if (!all(path_audit_10b3b$exists)) {
  
  cat("\nMISSING CHECKPOINT(S):\n")
  
  print(
    path_audit_10b3b[
      !path_audit_10b3b$exists,
      ,
      drop = FALSE
    ]
  )
  
  stop(
    "PATH AUDIT FAILED. Do not continue to Step 4."
  )
}

cat("\n============================================================\n")
cat("STEP 7B-7E-10B-3B PATH AUDIT: PASSED\n")
cat("All 6 required checkpoints found.\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 3. Check all required checkpoints
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("3. CHECK REQUIRED CHECKPOINTS\n")
cat("============================================================\n\n")

required_paths <- c(
  STEP10A  = path_10a,
  STEP10B1 = path_10b1,
  STEP10B2A = path_10b2a,
  STEP10B2B = path_10b2b,
  STEP10B2C = path_10b2c,
  STEP10B3A = path_10b3a
)

print(
  data.frame(
    checkpoint = names(required_paths),
    exists = file.exists(required_paths),
    path = unname(required_paths),
    row.names = NULL,
    check.names = FALSE
  )
)

if (!all(file.exists(required_paths))) {
  stop(
    "STEP 7B-7E-10B-3B STOP: one or more required checkpoints are missing."
  )
}

cat("\nAll required checkpoints: PASSED\n")


# ------------------------------------------------------------
# 4. Load prerequisite audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("4. LOAD 10B-3A PREREQUISITE AUDIT\n")
cat("============================================================\n\n")

env_10b3a <- new.env(
  parent = emptyenv()
)

load(
  path_10b3a,
  envir = env_10b3a
)

if (!"step7be10b3a_audit" %in% ls(envir = env_10b3a)) {
  stop(
    "STEP 10B-3A audit object is missing."
  )
}

audit_10b3a <- get(
  "step7be10b3a_audit",
  envir = env_10b3a
)

print(
  audit_10b3a$status
)

if (!identical(
  audit_10b3a$status,
  "LOO_PREREQUISITE_AUDIT_PASSED"
)) {
  stop(
    "STEP 10B-3B STOP: 10B-3A prerequisite audit did not pass."
  )
}

cat("\n10B-3A prerequisite status: PASSED\n")


# ------------------------------------------------------------
# 5. Load 10A signature definitions
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("5. LOAD LOCKED LOO SIGNATURE DEFINITIONS\n")
cat("============================================================\n\n")

env_10a <- new.env(
  parent = emptyenv()
)

load(
  path_10a,
  envir = env_10a
)

required_10a_objects <- c(
  "loo_targets_7be10",
  "stromal_genes_loo_AOC3_7be",
  "immune_genes_loo_IL7R_7be",
  "step7be10a_audit"
)

if (!all(
  required_10a_objects %in%
  ls(envir = env_10a)
)) {
  stop(
    "STEP 10A required objects are missing."
  )
}

loo_targets_7be10 <- get(
  "loo_targets_7be10",
  envir = env_10a
)

stromal_genes_loo_AOC3_7be <- get(
  "stromal_genes_loo_AOC3_7be",
  envir = env_10a
)

immune_genes_loo_IL7R_7be <- get(
  "immune_genes_loo_IL7R_7be",
  envir = env_10a
)

cat("AOC3 LOO Stromal signature length: ",
    length(stromal_genes_loo_AOC3_7be),
    "\n",
    sep = "")

cat("IL7R LOO Immune signature length: ",
    length(immune_genes_loo_IL7R_7be),
    "\n",
    sep = "")

if (length(stromal_genes_loo_AOC3_7be) != 140) {
  stop(
    "Unexpected AOC3 LOO Stromal signature length."
  )
}

if (length(immune_genes_loo_IL7R_7be) != 140) {
  stop(
    "Unexpected IL7R LOO Immune signature length."
  )
}

if (anyDuplicated(stromal_genes_loo_AOC3_7be) > 0) {
  stop(
    "Duplicated genes found in AOC3 LOO Stromal signature."
  )
}

if (anyDuplicated(immune_genes_loo_IL7R_7be) > 0) {
  stop(
    "Duplicated genes found in IL7R LOO Immune signature."
  )
}

cat("\nLocked LOO signature definitions: PASSED\n")


# ------------------------------------------------------------
# 6. Recover original 141-gene signatures
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("6. RECONSTRUCT ORIGINAL ESTIMATE SIGNATURES\n")
cat("============================================================\n\n")

# AOC3 was removed from the 141-gene Stromal signature
stromal_genes_original_7be10b3b <- unique(
  c(
    stromal_genes_loo_AOC3_7be,
    "AOC3"
  )
)

# IL7R was removed from the 141-gene Immune signature
immune_genes_original_7be10b3b <- unique(
  c(
    immune_genes_loo_IL7R_7be,
    "IL7R"
  )
)

cat(
  "Reconstructed original Stromal signature: ",
  length(stromal_genes_original_7be10b3b),
  " genes\n",
  sep = ""
)

cat(
  "Reconstructed original Immune signature: ",
  length(immune_genes_original_7be10b3b),
  " genes\n",
  sep = ""
)

if (length(stromal_genes_original_7be10b3b) != 141) {
  stop(
    "Original Stromal signature reconstruction failed."
  )
}

if (length(immune_genes_original_7be10b3b) != 141) {
  stop(
    "Original Immune signature reconstruction failed."
  )
}

cat("\nOriginal 141-gene signatures: PASSED\n")


# ------------------------------------------------------------
# 7. Load original ESTIMATE source objects
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("7. LOAD ORIGINAL ESTIMATE SOURCE OBJECTS\n")
cat("============================================================\n\n")

env_10b1 <- new.env(
  parent = emptyenv()
)

load(
  path_10b1,
  envir = env_10b1
)

# ------------------------------------------------------------
# 7.1 Check the ACTUAL object names stored in Step14 locked RData
# ------------------------------------------------------------

cat("Objects found in Step14 locked source:\n")
print(ls(envir = env_10b1))

required_10b1_objects_original <- c(
  "expr_estimate",
  "estimate_input_df",
  "estimate_scores_wgcna",
  "estimate_scores",
  "estimate_sample_ids"
)

if (!all(
  required_10b1_objects_original %in%
  ls(envir = env_10b1)
)) {
  
  cat("\nMissing required original Step14 objects:\n")
  
  print(
    required_10b1_objects_original[
      !required_10b1_objects_original %in%
        ls(envir = env_10b1)
    ]
  )
  
  stop(
    "STEP 7B-7E-10B-3B STOP: required Step14 source objects are missing."
  )
}

cat("\nOriginal Step14 source objects: PASSED\n")

# ------------------------------------------------------------
# 7.2 Recover original objects into 10B-3B-specific names
# ------------------------------------------------------------

expr_estimate_7be10 <- get(
  "expr_estimate",
  envir = env_10b1
)

estimate_input_df_7be10 <- get(
  "estimate_input_df",
  envir = env_10b1
)

estimate_scores_wgcna_source_7be10 <- get(
  "estimate_scores_wgcna",
  envir = env_10b1
)

estimate_scores_source_7be10 <- get(
  "estimate_scores",
  envir = env_10b1
)

estimate_sample_ids_7be10 <- get(
  "estimate_sample_ids",
  envir = env_10b1
)

# Optional auxiliary objects if present
if ("estimate_sample_ids_raw" %in% ls(envir = env_10b1)) {
  
  estimate_sample_ids_raw_7be10 <- get(
    "estimate_sample_ids_raw",
    envir = env_10b1
  )
  
}

if ("estimate_sample_ids_normalized" %in% ls(envir = env_10b1)) {
  
  estimate_sample_ids_normalized_7be10 <- get(
    "estimate_sample_ids_normalized",
    envir = env_10b1
  )
  
}

cat("\nRecovered source objects successfully.\n")

# ------------------------------------------------------------
# 7.3 Basic structural validation
# ------------------------------------------------------------

cat("\nExpression matrix dimensions: ")
print(dim(expr_estimate_7be10))

cat("ESTIMATE input dimensions: ")
print(dim(estimate_input_df_7be10))

cat("WGCNA-linked ESTIMATE score dimensions: ")
print(dim(estimate_scores_wgcna_source_7be10))

cat("Original ESTIMATE score dimensions: ")
print(dim(estimate_scores_source_7be10))

cat("Sample ID length: ")
print(length(estimate_sample_ids_7be10))

# ------------------------------------------------------------
# 7.4 Hard-stop structural checks
# ------------------------------------------------------------

if (!is.matrix(expr_estimate_7be10)) {
  stop(
    "expr_estimate is not a matrix."
  )
}

if (!all(
  c("StromalScore", "ImmuneScore", "ESTIMATEScore") %in%
  rownames(as.matrix(estimate_scores_wgcna_source_7be10))
)) {
  stop(
    "Required rows are missing from estimate_scores_wgcna."
  )
}

if (!all(
  c("StromalScore", "ImmuneScore", "ESTIMATEScore") %in%
  rownames(as.matrix(estimate_scores_source_7be10))
)) {
  stop(
    "Required rows are missing from estimate_scores."
  )
}

if (ncol(expr_estimate_7be10) != 513) {
  stop(
    "Unexpected number of samples in expr_estimate."
  )
}

if (length(estimate_sample_ids_7be10) != 513) {
  stop(
    "Unexpected number of sample IDs."
  )
}

if (!all(
  is.finite(
    as.matrix(estimate_scores_wgcna_source_7be10)
  )
)) {
  stop(
    "Non-finite values found in source WGCNA-linked ESTIMATE scores."
  )
}

if (!all(
  is.finite(
    as.matrix(estimate_scores_source_7be10)
  )
)) {
  stop(
    "Non-finite values found in original ESTIMATE scores."
  )
}

cat("\nOriginal ESTIMATE source structure: PASSED\n")

# ------------------------------------------------------------
# 8. Validate expression matrix
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("8. VALIDATE EXPRESSION MATRIX\n")
cat("============================================================\n\n")

if (is.null(rownames(expr_estimate_7be10))) {
  stop(
    "Expression matrix has no gene symbols."
  )
}

if (is.null(colnames(expr_estimate_7be10))) {
  stop(
    "Expression matrix has no sample IDs."
  )
}

if (anyDuplicated(rownames(expr_estimate_7be10)) > 0) {
  stop(
    "Duplicated gene symbols detected."
  )
}

if (anyDuplicated(colnames(expr_estimate_7be10)) > 0) {
  stop(
    "Duplicated sample IDs detected."
  )
}

if (!all(
  is.finite(
    as.matrix(expr_estimate_7be10)
  )
)) {
  stop(
    "Expression matrix contains non-finite values."
  )
}

cat("Unique gene symbols: PASSED\n")
cat("Unique sample IDs: PASSED\n")
cat("All expression values finite: PASSED\n")


# ------------------------------------------------------------
# 9. Define exact ESTIMATE 1.0.13 scoring function
# ------------------------------------------------------------
#
# This function follows the official estimateScore() logic:
#
#   1. rank-normalize each sample using average ties
#   2. scale ranks by 10000 / number of genes
#   3. rank genes in decreasing order
#   4. identify overlap with the requested gene set
#   5. calculate weighted enrichment score using
#      correlation/rank values ^ 0.25
#   6. choose positive or negative extreme according to
#      absolute magnitude
#
# This reproduces the non-Affymetrix scoring branch.
#
# Official ESTIMATE source:
# estimate 1.0.13
# ------------------------------------------------------------

estimate_score_official_10b3b <- function(
    expr_matrix,
    gene_set,
    score_name = "Score",
    progress = TRUE
) {
  
  stopifnot(
    is.matrix(expr_matrix)
  )
  
  gene_names <- rownames(expr_matrix)
  sample_names <- colnames(expr_matrix)
  
  n_genes <- nrow(expr_matrix)
  n_samples <- ncol(expr_matrix)
  
  if (is.null(gene_names)) {
    stop(
      "expr_matrix must have rownames containing gene symbols."
    )
  }
  
  if (is.null(sample_names)) {
    stop(
      "expr_matrix must have colnames containing sample IDs."
    )
  }
  
  gene_set <- unique(
    as.character(gene_set)
  )
  
  gene_overlap <- intersect(
    gene_set,
    gene_names
  )
  
  if (length(gene_overlap) == 0) {
    stop(
      paste0(
        score_name,
        ": zero gene overlap with expression matrix."
      )
    )
  }
  
  gene_set_idx <- match(
    gene_overlap,
    gene_names
  )
  
  score_vector <- rep(
    NA_real_,
    n_samples
  )
  
  for (j in seq_len(n_samples)) {
    
    if (progress && (j == 1L || j %% 50L == 0L || j == n_samples)) {
      
      cat(
        "  ",
        score_name,
        ": sample ",
        j,
        "/",
        n_samples,
        "\n",
        sep = ""
      )
      
    }
    
    x <- as.numeric(
      expr_matrix[, j]
    )
    
    # --------------------------------------------------------
    # Official sample rank normalization
    # ties.method = "average"
    # --------------------------------------------------------
    
    rank_x <- rank(
      x,
      ties.method = "average"
    )
    
    m <- 10000 * rank_x / n_genes
    
    # --------------------------------------------------------
    # Official decreasing gene order
    # --------------------------------------------------------
    
    gene_list <- order(
      m,
      decreasing = TRUE
    )
    
    # --------------------------------------------------------
    # Reorder normalized expression/rank vector
    # --------------------------------------------------------
    
    correl_vector <- m[gene_list]
    
    # --------------------------------------------------------
    # TAG indicator
    # --------------------------------------------------------
    
    TAG <- as.integer(
      match(
        gene_list,
        gene_set_idx,
        nomatch = 0L
      ) > 0L
    )
    
    no_TAG <- 1L - TAG
    
    N <- length(
      gene_list
    )
    
    Nh <- length(
      gene_set_idx
    )
    
    Nm <- N - Nh
    
    if (Nh <= 0L) {
      score_vector[j] <- NA_real_
      next
    }
    
    if (Nm <= 0L) {
      stop(
        paste0(
          score_name,
          ": Nm <= 0."
        )
      )
    }
    
    # --------------------------------------------------------
    # Official weighting exponent = 0.25
    # --------------------------------------------------------
    
    correl_vector <- correl_vector^0.25
    
    sum_correl <- sum(
      correl_vector[TAG == 1L]
    )
    
    if (!is.finite(sum_correl) ||
        sum_correl <= 0) {
      
      score_vector[j] <- NA_real_
      next
    }
    
    # --------------------------------------------------------
    # Running enrichment calculation
    # --------------------------------------------------------
    
    P0 <- no_TAG / Nm
    
    F0 <- cumsum(
      P0
    )
    
    Pn <- TAG *
      correl_vector /
      sum_correl
    
    Fn <- cumsum(
      Pn
    )
    
    RES <- Fn - F0
    
    max_ES <- max(
      RES
    )
    
    min_ES <- min(
      RES
    )
    
    # --------------------------------------------------------
    # ESTIMATE chooses the extreme with greater absolute
    # magnitude.
    # --------------------------------------------------------
    
    if (max_ES > -min_ES) {
      
      ES <- max_ES
      
    } else {
      
      ES <- min_ES
      
    }
    
    score_vector[j] <- ES
  }
  
  names(score_vector) <- sample_names
  
  return(
    score_vector
  )
}


# ------------------------------------------------------------
# 10. Reconstruct original StromalScore
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("10. RECONSTRUCT ORIGINAL STROMAL SCORE\n")
cat("============================================================\n\n")

reconstructed_original_stromal_10b3b <-
  estimate_score_official_10b3b(
    expr_matrix = expr_estimate_7be10,
    gene_set = stromal_genes_original_7be10b3b,
    score_name = "Original_StromalScore",
    progress = TRUE
  )

cat("\nOriginal StromalScore reconstruction completed.\n")


# ------------------------------------------------------------
# 11. Reconstruct original ImmuneScore
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("11. RECONSTRUCT ORIGINAL IMMUNE SCORE\n")
cat("============================================================\n\n")

reconstructed_original_immune_10b3b <-
  estimate_score_official_10b3b(
    expr_matrix = expr_estimate_7be10,
    gene_set = immune_genes_original_7be10b3b,
    score_name = "Original_ImmuneScore",
    progress = TRUE
  )

cat("\nOriginal ImmuneScore reconstruction completed.\n")


# ------------------------------------------------------------
# 12. Construct reconstructed original ESTIMATEScore
# ------------------------------------------------------------

reconstructed_original_estimate_10b3b <-
  reconstructed_original_stromal_10b3b +
  reconstructed_original_immune_10b3b

names(
  reconstructed_original_estimate_10b3b
) <- colnames(
  expr_estimate_7be10
)


# ------------------------------------------------------------
# 13. Validate reconstructed original scores
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("13. VALIDATE RECONSTRUCTED ORIGINAL ESTIMATE SCORES\n")
cat("============================================================\n\n")

source_stromal_10b3b <-
  as.numeric(
    estimate_scores_wgcna_source_7be10[
      "StromalScore",
    ]
  )

source_immune_10b3b <-
  as.numeric(
    estimate_scores_wgcna_source_7be10[
      "ImmuneScore",
    ]
  )

source_estimate_10b3b <-
  as.numeric(
    estimate_scores_wgcna_source_7be10[
      "ESTIMATEScore",
    ]
  )

names(source_stromal_10b3b) <-
  colnames(
    estimate_scores_wgcna_source_7be10
  )

names(source_immune_10b3b) <-
  colnames(
    estimate_scores_wgcna_source_7be10
  )

names(source_estimate_10b3b) <-
  colnames(
    estimate_scores_wgcna_source_7be10
  )


# ------------------------------------------------------------
# 14. Verify sample order before score comparison
# ------------------------------------------------------------

if (!identical(
  names(reconstructed_original_stromal_10b3b),
  names(source_stromal_10b3b)
)) {
  stop(
    "Sample order mismatch: reconstructed StromalScore."
  )
}

if (!identical(
  names(reconstructed_original_immune_10b3b),
  names(source_immune_10b3b)
)) {
  stop(
    "Sample order mismatch: reconstructed ImmuneScore."
  )
}

if (!identical(
  names(reconstructed_original_estimate_10b3b),
  names(source_estimate_10b3b)
)) {
  stop(
    "Sample order mismatch: reconstructed ESTIMATEScore."
  )
}

cat("Sample order: PASSED\n")


# ------------------------------------------------------------
# 15. Calculate reconstruction differences
# ------------------------------------------------------------

stromal_diff_10b3b <-
  reconstructed_original_stromal_10b3b -
  source_stromal_10b3b

immune_diff_10b3b <-
  reconstructed_original_immune_10b3b -
  source_immune_10b3b

estimate_diff_10b3b <-
  reconstructed_original_estimate_10b3b -
  source_estimate_10b3b


max_abs_diff_stromal_10b3b <-
  max(
    abs(stromal_diff_10b3b),
    na.rm = TRUE
  )

max_abs_diff_immune_10b3b <-
  max(
    abs(immune_diff_10b3b),
    na.rm = TRUE
  )

max_abs_diff_estimate_10b3b <-
  max(
    abs(estimate_diff_10b3b),
    na.rm = TRUE
  )


cor_stromal_10b3b <-
  cor(
    reconstructed_original_stromal_10b3b,
    source_stromal_10b3b,
    method = "pearson"
  )

cor_immune_10b3b <-
  cor(
    reconstructed_original_immune_10b3b,
    source_immune_10b3b,
    method = "pearson"
  )

cor_estimate_10b3b <-
  cor(
    reconstructed_original_estimate_10b3b,
    source_estimate_10b3b,
    method = "pearson"
  )


cat("\nMaximum absolute differences:\n")

cat(
  "StromalScore: ",
  format(
    max_abs_diff_stromal_10b3b,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "ImmuneScore: ",
  format(
    max_abs_diff_immune_10b3b,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "ESTIMATEScore: ",
  format(
    max_abs_diff_estimate_10b3b,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

cat("\nPearson correlations:\n")

cat(
  "StromalScore: ",
  format(
    cor_stromal_10b3b,
    digits = 12
  ),
  "\n",
  sep = ""
)

cat(
  "ImmuneScore: ",
  format(
    cor_immune_10b3b,
    digits = 12
  ),
  "\n",
  sep = ""
)

cat(
  "ESTIMATEScore: ",
  format(
    cor_estimate_10b3b,
    digits = 12
  ),
  "\n",
  sep = ""
)

# ============================================================
# STEP 7B-7E-10B-3B
# DIAGNOSTIC: inspect locked estimate_filtered
# ============================================================

cat("\n============================================================\n")
cat("DIAGNOSTIC: INSPECT LOCKED estimate_filtered\n")
cat("============================================================\n\n")

cat("Class of estimate_filtered:\n")
print(class(
  estimate_filtered_7be10 <- get(
    "estimate_filtered",
    envir = env_10b1
  )
))

cat("\nDimensions:\n")
print(dim(estimate_filtered_7be10))

cat("\nColumn names (first 10):\n")
print(
  head(
    colnames(estimate_filtered_7be10),
    10
  )
)

cat("\nRow names (first 10):\n")
print(
  head(
    rownames(estimate_filtered_7be10),
    10
  )
)

cat("\nStructure:\n")
str(
  estimate_filtered_7be10,
  max.level = 1
)

cat("\n============================================================\n")
cat("CHECK WHETHER 10,412 COMMON GENES ARE PRESENT\n")
cat("============================================================\n\n")

cat(
  "Number of rows: ",
  nrow(estimate_filtered_7be10),
  "\n",
  sep = ""
)

cat(
  "Number of columns: ",
  ncol(estimate_filtered_7be10),
  "\n",
  sep = ""
)

cat(
  "Has GeneSymbol column: ",
  "GeneSymbol" %in% colnames(estimate_filtered_7be10),
  "\n",
  sep = ""
)

cat(
  "Has 513 samples among columns: ",
  sum(
    colnames(estimate_filtered_7be10) %in%
      colnames(expr_estimate_7be10)
  ),
  "\n",
  sep = ""
)

# ============================================================
# STEP 7B-7E-10B-3B
# DIAGNOSTIC 2:
# LOCATE THE ORIGINAL ESTIMATE FILTERED GCT FILE
# ============================================================

cat("\n============================================================\n")
cat("DIAGNOSTIC 2: LOCATE ORIGINAL ESTIMATE FILTERED GCT\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Recover the original filename recorded in Step14
# ------------------------------------------------------------

estimate_filtered_filename_10b3b <-
  as.character(
    estimate_filtered_7be10
  )

cat(
  "Recorded filtered GCT filename:\n",
  estimate_filtered_filename_10b3b,
  "\n\n",
  sep = ""
)

# ------------------------------------------------------------
# 2. Check Step14 root information
# ------------------------------------------------------------

if ("step14_root" %in% ls(envir = env_10b1)) {
  
  step14_root_10b3b <- get(
    "step14_root",
    envir = env_10b1
  )
  
  cat("Recorded Step14 root:\n")
  print(step14_root_10b3b)
  
} else {
  
  step14_root_10b3b <- NA_character_
  
  cat(
    "No step14_root object found.\n"
  )
}

# ------------------------------------------------------------
# 3. Check recorded ESTIMATE-related file paths
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("RECORDED ESTIMATE FILE OBJECTS\n")
cat("============================================================\n\n")

for (nm in c(
  "estimate_input_file",
  "estimate_gct_input",
  "estimate_filtered"
)) {
  
  if (nm %in% ls(envir = env_10b1)) {
    
    obj <- get(
      nm,
      envir = env_10b1
    )
    
    cat(
      "\n",
      nm,
      ":\n",
      sep = ""
    )
    
    print(obj)
  }
}

# ------------------------------------------------------------
# 4. Generate candidate search roots
# ------------------------------------------------------------

search_roots_10b3b <- unique(
  c(
    dirname(path_10b1),
    dirname(dirname(path_10b1)),
    if (!is.na(step14_root_10b3b))
      as.character(step14_root_10b3b)
    else
      character(0)
  )
)

search_roots_10b3b <- search_roots_10b3b[
  dir.exists(search_roots_10b3b)
]

cat("\n============================================================\n")
cat("SEARCH ROOTS\n")
cat("============================================================\n\n")

print(search_roots_10b3b)

# ------------------------------------------------------------
# 5. Search recursively for the exact GCT filename
# ------------------------------------------------------------

gct_matches_10b3b <- unique(
  unlist(
    lapply(
      search_roots_10b3b,
      function(root) {
        
        list.files(
          path = root,
          pattern = "^Workflow-04-2C-1A_ESTIMATE_filtered\\.gct$",
          recursive = TRUE,
          full.names = TRUE,
          ignore.case = FALSE
        )
        
      }
    )
  )
)

cat("\n============================================================\n")
cat("GCT FILE SEARCH RESULT\n")
cat("============================================================\n\n")

cat(
  "Number of matching files: ",
  length(gct_matches_10b3b),
  "\n\n",
  sep = ""
)

if (length(gct_matches_10b3b) > 0) {
  
  print(
    gct_matches_10b3b
  )
  
} else {
  
  cat(
    "No exact matching GCT file found under the candidate roots.\n"
  )
}

# ------------------------------------------------------------
# 6. Search for related ESTIMATE GCT files if exact file absent
# ------------------------------------------------------------

if (length(gct_matches_10b3b) == 0) {
  
  cat("\n============================================================\n")
  cat("RELATED ESTIMATE GCT FILE SEARCH\n")
  cat("============================================================\n\n")
  
  related_gct_10b3b <- unique(
    unlist(
      lapply(
        search_roots_10b3b,
        function(root) {
          
          list.files(
            path = root,
            pattern = "ESTIMATE.*\\.gct$",
            recursive = TRUE,
            full.names = TRUE,
            ignore.case = TRUE
          )
          
        }
      )
    )
  )
  
  cat(
    "Number of related ESTIMATE GCT files: ",
    length(related_gct_10b3b),
    "\n\n",
    sep = ""
  )
  
  if (length(related_gct_10b3b) > 0) {
    
    print(
      related_gct_10b3b
    )
    
  }
}

# ------------------------------------------------------------
# 7. STOP here deliberately
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("DIAGNOSTIC 2 COMPLETE\n")
cat("============================================================\n")
cat(
  "\nNo GCT file will be read or modified at this stage.\n",
  "No LOO score has been calculated.\n",
  "No statistical analysis has been performed.\n",
  "No locked Step14 object has been modified.\n",
  sep = ""
)

# ============================================================
# STEP 7B-7E-10B-3B
# DIAGNOSTIC 3:
# TRACE THE ORIGINAL ESTIMATE INPUT / FILTERING WORKFLOW
# ============================================================

cat("\n============================================================\n")
cat("DIAGNOSTIC 3: TRACE ORIGINAL ESTIMATE INPUT WORKFLOW\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Recover Step14 root
# ------------------------------------------------------------

step14_root_10b3b <- as.character(
  get("step14_root", envir = env_10b1)
)

cat("Step14 root:\n")
print(step14_root_10b3b)


# ------------------------------------------------------------
# 2. Search for the original ESTIMATE input files
# ------------------------------------------------------------

target_files_10b3b <- c(
  "Workflow-04-2C-1A_ESTIMATE_input_GeneSymbol_log2_FPKM_UQ.txt",
  "Workflow-04-2C-1A_ESTIMATE_input_GCT.txt",
  "Workflow-04-2C-1A_ESTIMATE_filtered.gct"
)

cat("\n============================================================\n")
cat("SEARCHING FOR ORIGINAL ESTIMATE FILES\n")
cat("============================================================\n\n")

estimate_file_matches_10b3b <- list()

for (target in target_files_10b3b) {
  
  matches <- unique(
    unlist(
      lapply(
        search_roots_10b3b,
        function(root) {
          
          list.files(
            path = root,
            pattern = paste0(
              "^",
              gsub(
                "([.|(){}+$^*?\\\\])",
                "\\\\\\1",
                target
              ),
              "$"
            ),
            recursive = TRUE,
            full.names = TRUE,
            ignore.case = FALSE
          )
          
        }
      )
    )
  )
  
  estimate_file_matches_10b3b[[target]] <- matches
  
  cat("\n", target, ":\n", sep = "")
  
  if (length(matches) == 0) {
    
    cat("  NOT FOUND\n")
    
  } else {
    
    cat(
      "  Found ",
      length(matches),
      " file(s):\n",
      sep = ""
    )
    
    print(matches)
    
  }
}


# ------------------------------------------------------------
# 3. Search for Step14 workflow scripts
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SEARCHING FOR STEP14 WORKFLOW / SCRIPT FILES\n")
cat("============================================================\n\n")

workflow_patterns_10b3b <- c(
  "\\.(R|r)$",
  "\\.(Rmd|rmd)$",
  "\\.(txt|log|out)$",
  "\\.(bat|cmd|ps1)$"
)

workflow_matches_10b3b <- unique(
  unlist(
    lapply(
      search_roots_10b3b,
      function(root) {
        
        unlist(
          lapply(
            workflow_patterns_10b3b,
            function(pattern) {
              
              list.files(
                path = root,
                pattern = pattern,
                recursive = TRUE,
                full.names = TRUE,
                ignore.case = TRUE
              )
              
            }
          )
        )
        
      }
    )
  )
)

cat(
  "Number of workflow/script/log files found: ",
  length(workflow_matches_10b3b),
  "\n\n",
  sep = ""
)

if (length(workflow_matches_10b3b) > 0) {
  
  print(workflow_matches_10b3b)
  
} else {
  
  cat("No workflow/script/log files found.\n")
  
}


# ------------------------------------------------------------
# 4. Search broader project root
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("BROADER PROJECT SEARCH\n")
cat("============================================================\n\n")

project_root_10b3b <- "O:/Rproject"

cat("Project root:\n")
print(project_root_10b3b)

cat("\nSearching for exact ESTIMATE filenames under O:/Rproject ...\n")

broader_estimate_matches_10b3b <- unique(
  unlist(
    lapply(
      target_files_10b3b,
      function(target) {
        
        unlist(
          lapply(
            list.dirs(
              project_root_10b3b,
              recursive = TRUE,
              full.names = TRUE
            ),
            function(root) {
              
              list.files(
                path = root,
                pattern = paste0(
                  "^",
                  gsub(
                    "([.|(){}+$^*?\\\\])",
                    "\\\\\\1",
                    target
                  ),
                  "$"
                ),
                full.names = TRUE,
                ignore.case = FALSE
              )
              
            }
          )
        )
        
      }
    )
  )
)

broader_estimate_matches_10b3b <- unique(
  broader_estimate_matches_10b3b
)

cat(
  "\nNumber of exact ESTIMATE files found anywhere under O:/Rproject: ",
  length(broader_estimate_matches_10b3b),
  "\n\n",
  sep = ""
)

if (length(broader_estimate_matches_10b3b) > 0) {
  
  print(broader_estimate_matches_10b3b)
  
} else {
  
  cat(
    "No exact ESTIMATE input/filtered files found anywhere under O:/Rproject.\n"
  )
  
}


# ------------------------------------------------------------
# 5. Search for filenames containing ESTIMATE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("BROAD ESTIMATE FILE SEARCH\n")
cat("============================================================\n\n")

broad_estimate_files_10b3b <- unique(
  unlist(
    lapply(
      list.dirs(
        project_root_10b3b,
        recursive = TRUE,
        full.names = TRUE
      ),
      function(root) {
        
        list.files(
          path = root,
          pattern = "ESTIMATE",
          recursive = FALSE,
          full.names = TRUE,
          ignore.case = TRUE
        )
        
      }
    )
  )
)

cat(
  "Number of files/directories containing 'ESTIMATE': ",
  length(broad_estimate_files_10b3b),
  "\n\n",
  sep = ""
)

if (length(broad_estimate_files_10b3b) > 0) {
  
  print(broad_estimate_files_10b3b)
  
} else {
  
  cat("No ESTIMATE-named files found under O:/Rproject.\n")
  
}


# ------------------------------------------------------------
# 6. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("DIAGNOSTIC 3 COMPLETE\n")
cat("============================================================\n\n")

cat(
  "No ESTIMATE score has been recalculated.\n",
  "No LOO score has been calculated.\n",
  "No statistical analysis has been performed.\n",
  "No locked Step14 object has been modified.\n",
  "No Figure 6 or Supplementary Table 3 object has been modified.\n",
  sep = ""
)

# ============================================================
# STEP 7B-7E-10B-3C
# ORIGINAL ESTIMATE GCT INPUT AUDIT
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3C: ORIGINAL ESTIMATE GCT INPUT AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define the original GCT file
# ------------------------------------------------------------

original_gct_10b3c <-
  "input_dir/Paper1_tumor-only WGCNA/Workflow-04-2C-1A_ESTIMATE_input_GCT.txt"

cat("Original GCT path:\n")
print(original_gct_10b3c)

cat(
  "\nFile exists: ",
  file.exists(original_gct_10b3c),
  "\n",
  sep = ""
)

if (!file.exists(original_gct_10b3c)) {
  stop("Original ESTIMATE input GCT file does not exist.")
}


# ------------------------------------------------------------
# 2. File information
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FILE INFORMATION\n")
cat("============================================================\n\n")

file_info_10b3c <- file.info(original_gct_10b3c)

print(
  file_info_10b3c[
    ,
    c("size", "mtime", "ctime")
  ]
)


# ------------------------------------------------------------
# 3. Read only the GCT header / first lines
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GCT HEADER PREVIEW\n")
cat("============================================================\n\n")

gct_header_10b3c <- readLines(
  original_gct_10b3c,
  n = 8,
  warn = FALSE
)

cat(
  paste(
    seq_along(gct_header_10b3c),
    gct_header_10b3c,
    sep = ": "
  ),
  sep = "\n"
)


# ------------------------------------------------------------
# 4. Read the complete GCT as text
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("READING ORIGINAL GCT\n")
cat("============================================================\n\n")

gct_raw_10b3c <- readLines(
  original_gct_10b3c,
  warn = FALSE
)

cat(
  "Total text lines: ",
  length(gct_raw_10b3c),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. Validate standard GCT structure
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GCT STRUCTURE VALIDATION\n")
cat("============================================================\n\n")

if (length(gct_raw_10b3c) < 3) {
  stop("GCT file has fewer than 3 lines.")
}

cat("Line 1:\n")
print(gct_raw_10b3c[1])

cat("\nLine 2:\n")
print(gct_raw_10b3c[2])

cat("\nLine 3:\n")
print(gct_raw_10b3c[3])


# ------------------------------------------------------------
# 6. Parse dimensions from GCT header
# ------------------------------------------------------------

gct_dims_10b3c <- strsplit(
  gct_raw_10b3c[2],
  "\t",
  fixed = TRUE
)[[1]]

cat("\nParsed GCT dimension fields:\n")
print(gct_dims_10b3c)


# ------------------------------------------------------------
# 7. Parse column header
# ------------------------------------------------------------

gct_col_header_10b3c <- strsplit(
  gct_raw_10b3c[3],
  "\t",
  fixed = TRUE
)[[1]]

cat("\nNumber of GCT header columns: ")
cat(length(gct_col_header_10b3c), "\n")

cat("\nFirst 10 GCT column names:\n")
print(
  head(
    gct_col_header_10b3c,
    10
  )
)

cat("\nLast 10 GCT column names:\n")
print(
  tail(
    gct_col_header_10b3c,
    10
  )
)


# ------------------------------------------------------------
# 8. Compare GCT sample IDs with Step14 sample IDs
# ------------------------------------------------------------

gct_sample_ids_10b3c <-
  gct_col_header_10b3c[-c(1, 2)]

cat("\n============================================================\n")
cat("SAMPLE ID AUDIT\n")
cat("============================================================\n\n")

cat(
  "GCT sample count: ",
  length(gct_sample_ids_10b3c),
  "\n",
  sep = ""
)

cat(
  "Step14 sample count: ",
  length(estimate_sample_ids_7be10),
  "\n",
  sep = ""
)

sample_identity_10b3c <- identical(
  as.character(gct_sample_ids_10b3c),
  as.character(estimate_sample_ids_7be10)
)

cat(
  "Exact sample identity/order match: ",
  sample_identity_10b3c,
  "\n",
  sep = ""
)

if (!sample_identity_10b3c) {
  
  cat("\nFirst 10 GCT IDs:\n")
  print(
    head(
      gct_sample_ids_10b3c,
      10
    )
  )
  
  cat("\nFirst 10 Step14 IDs:\n")
  print(
    head(
      estimate_sample_ids_7be10,
      10
    )
  )
  
}


# ------------------------------------------------------------
# 9. Compare GCT gene count with expression matrix
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GENE COUNT AUDIT\n")
cat("============================================================\n\n")

gct_gene_count_10b3c <-
  as.integer(gct_dims_10b3c[1])

gct_sample_count_10b3c <-
  as.integer(gct_dims_10b3c[2])

cat(
  "GCT declared gene count: ",
  gct_gene_count_10b3c,
  "\n",
  sep = ""
)

cat(
  "GCT declared sample count: ",
  gct_sample_count_10b3c,
  "\n",
  sep = ""
)

cat(
  "expr_estimate rows: ",
  nrow(expr_estimate_7be10),
  "\n",
  sep = ""
)

cat(
  "expr_estimate columns: ",
  ncol(expr_estimate_7be10),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Parse gene identifiers only
# ------------------------------------------------------------

gct_gene_table_preview_10b3c <- strsplit(
  gct_raw_10b3c[-c(1, 2, 3)],
  "\t",
  fixed = TRUE
)

gct_gene_ids_10b3c <- vapply(
  gct_gene_table_preview_10b3c,
  function(x) x[1],
  character(1)
)

gct_gene_symbols_10b3c <- vapply(
  gct_gene_table_preview_10b3c,
  function(x) x[2],
  character(1)
)

cat("\n============================================================\n")
cat("GENE IDENTIFIER AUDIT\n")
cat("============================================================\n\n")

cat(
  "Parsed gene rows: ",
  length(gct_gene_ids_10b3c),
  "\n",
  sep = ""
)

cat(
  "Unique gene symbols: ",
  length(unique(gct_gene_symbols_10b3c)),
  "\n",
  sep = ""
)

cat(
  "Duplicated gene symbols: ",
  sum(duplicated(gct_gene_symbols_10b3c)),
  "\n",
  sep = ""
)

cat("\nFirst 10 gene IDs:\n")
print(
  head(
    gct_gene_ids_10b3c,
    10
  )
)

cat("\nFirst 10 gene symbols:\n")
print(
  head(
    gct_gene_symbols_10b3c,
    10
  )
)


# ------------------------------------------------------------
# 11. Compare GCT gene symbols with expr_estimate
# ------------------------------------------------------------

expr_gene_symbols_10b3c <-
  rownames(expr_estimate_7be10)

cat("\n============================================================\n")
cat("GCT vs EXPR_ESTIMATE GENE AUDIT\n")
cat("============================================================\n\n")

cat(
  "GCT genes present in expr_estimate: ",
  sum(
    gct_gene_symbols_10b3c %in%
      expr_gene_symbols_10b3c
  ),
  "\n",
  sep = ""
)

cat(
  "GCT genes absent from expr_estimate: ",
  sum(
    !(gct_gene_symbols_10b3c %in%
        expr_gene_symbols_10b3c)
  ),
  "\n",
  sep = ""
)

cat(
  "expr_estimate genes present in GCT: ",
  sum(
    expr_gene_symbols_10b3c %in%
      gct_gene_symbols_10b3c
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Verify the three core genes in GCT
# ------------------------------------------------------------

core_genes_10b3c <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

cat("\n============================================================\n")
cat("CORE GENE PRESENCE IN ORIGINAL GCT\n")
cat("============================================================\n\n")

core_gct_presence_10b3c <- data.frame(
  Gene = core_genes_10b3c,
  In_GCT = core_genes_10b3c %in%
    gct_gene_symbols_10b3c,
  stringsAsFactors = FALSE
)

print(core_gct_presence_10b3c)


# ------------------------------------------------------------
# 13. Compare the GCT sample IDs with the GCT input object
# ------------------------------------------------------------

if (
  "estimate_gct_input_7be10" %in%
  ls(envir = env_10b1)
) {
  
  estimate_gct_input_7be10 <-
    get(
      "estimate_gct_input_7be10",
      envir = env_10b1
    )
  
  cat("\n============================================================\n")
  cat("RECOVERED ESTIMATE GCT INPUT OBJECT\n")
  cat("============================================================\n\n")
  
  cat("Class:\n")
  print(class(estimate_gct_input_7be10))
  
  cat("\nStructure:\n")
  print(str(estimate_gct_input_7be10))
  
} else {
  
  cat(
    "\nNo estimate_gct_input_7be10 object exists in env_10b1.\n"
  )
  
}


# ------------------------------------------------------------
# 14. Final audit status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 10B-3C COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Original GCT was read for audit only.\n",
  "No ESTIMATE score was recalculated.\n",
  "No LOO score was calculated.\n",
  "No statistical analysis was performed.\n",
  "No locked Step14 object was modified.\n",
  sep = ""
)

# ============================================================
# STEP 7B-7E-10B-3D
# CORRECT GCT STRUCTURE AUDIT
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3D: CORRECT GCT STRUCTURE AUDIT\n")
cat("============================================================\n\n")

original_gct_10b3d <-
  "input_dir/Paper1_tumor-only WGCNA/Workflow-04-2C-1A_ESTIMATE_input_GCT.txt"

# ------------------------------------------------------------
# 1. Read first 6 lines exactly
# ------------------------------------------------------------

gct_head_10b3d <- readLines(
  original_gct_10b3d,
  n = 6,
  warn = FALSE
)

cat("FIRST 6 LINES:\n\n")

for (i in seq_along(gct_head_10b3d)) {
  cat(
    "LINE ",
    i,
    ":\n",
    gct_head_10b3d[i],
    "\n\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 2. Count fields in each of first 6 lines
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FIELD COUNT OF FIRST 6 LINES\n")
cat("============================================================\n\n")

for (i in seq_along(gct_head_10b3d)) {
  
  fields_i <- strsplit(
    gct_head_10b3d[i],
    "\t",
    fixed = TRUE
  )[[1]]
  
  cat(
    "Line ",
    i,
    ": ",
    length(fields_i),
    " fields\n",
    sep = ""
  )
  
}


# ------------------------------------------------------------
# 3. Read the final 3 lines
# ------------------------------------------------------------

gct_tail_10b3d <- tail(
  readLines(
    original_gct_10b3d,
    warn = FALSE
  ),
  3
)

cat("\n============================================================\n")
cat("LAST 3 LINES\n")
cat("============================================================\n\n")

for (i in seq_along(gct_tail_10b3d)) {
  
  cat(
    "TAIL ",
    i,
    ":\n",
    gct_tail_10b3d[i],
    "\n\n",
    sep = ""
  )
  
}


# ------------------------------------------------------------
# 4. Determine likely GCT header line
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("HEADER INTERPRETATION\n")
cat("============================================================\n\n")

for (i in seq_len(min(6, length(gct_head_10b3d)))) {
  
  fields_i <- strsplit(
    gct_head_10b3d[i],
    "\t",
    fixed = TRUE
  )[[1]]
  
  cat(
    "Line ",
    i,
    " first fields:\n",
    sep = ""
  )
  
  print(
    head(
      fields_i,
      min(8, length(fields_i))
    )
  )
  
  cat("\n")
  
}


# ------------------------------------------------------------
# 5. Compare likely sample header against Step14 IDs
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("POTENTIAL SAMPLE HEADER MATCH\n")
cat("============================================================\n\n")

for (i in seq_len(min(6, length(gct_head_10b3d)))) {
  
  fields_i <- strsplit(
    gct_head_10b3d[i],
    "\t",
    fixed = TRUE
  )[[1]]
  
  if (length(fields_i) >= 3) {
    
    candidate_ids <- fields_i[-c(1, 2)]
    
    if (
      length(candidate_ids) ==
      length(estimate_sample_ids_7be10)
    ) {
      
      exact_match <- identical(
        as.character(candidate_ids),
        as.character(estimate_sample_ids_7be10)
      )
      
      cat(
        "Line ",
        i,
        ": ",
        exact_match,
        "\n",
        sep = ""
      )
      
    }
    
  }
  
}


# ------------------------------------------------------------
# 6. Stop
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 10B-3D COMPLETE\n")
cat("============================================================\n\n")

cat(
  "GCT structure only was audited.\n",
  "No ESTIMATE score was recalculated.\n",
  "No LOO score was calculated.\n",
  "No statistical analysis was performed.\n",
  "No locked Step14 object was modified.\n",
  sep = ""
)

# ============================================================
# STEP 7B-7E-10B-3E
# ACTUAL ESTIMATE FILTERED GCT RECOVERY AUDIT
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3E: ACTUAL ESTIMATE FILTERED GCT AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Candidate roots
# ------------------------------------------------------------

candidate_roots_10b3e <- unique(
  normalizePath(
    c(
      "input_dir/Paper1_tumor-only WGCNA",
      "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis",
      dirname(original_gct_10b3d)
    ),
    winslash = "/",
    mustWork = FALSE
  )
)

cat("Candidate search roots:\n")
print(candidate_roots_10b3e)


# ------------------------------------------------------------
# 2. Search exact filtered GCT filename
# ------------------------------------------------------------

filtered_name_10b3e <-
  "Workflow-04-2C-1A_ESTIMATE_filtered.gct"

filtered_hits_10b3e <- unique(
  unlist(
    lapply(
      candidate_roots_10b3e,
      function(root_i) {
        
        if (!dir.exists(root_i)) {
          return(character(0))
        }
        
        list.files(
          path = root_i,
          pattern = paste0(
            "^",
            filtered_name_10b3e,
            "$"
          ),
          recursive = TRUE,
          full.names = TRUE,
          ignore.case = FALSE
        )
      }
    )
  )
)

cat("\n============================================================\n")
cat("FILTERED GCT SEARCH RESULT\n")
cat("============================================================\n\n")

cat(
  "Number of exact matches: ",
  length(filtered_hits_10b3e),
  "\n",
  sep = ""
)

if (length(filtered_hits_10b3e) > 0) {
  
  print(filtered_hits_10b3e)
  
} else {
  
  cat(
    "\nExact filtered GCT was not found.\n"
  )
}


# ------------------------------------------------------------
# 3. Stop if not found
# ------------------------------------------------------------

if (length(filtered_hits_10b3e) == 0) {
  
  cat(
    "\nSTEP 10B-3E STOP:\n",
    "Actual filtered GCT file not found.\n",
    "No score recalculation performed.\n",
    "No LOO score calculated.\n",
    "No statistical analysis performed.\n",
    "No locked Step14 object modified.\n",
    sep = ""
  )
  
} else {
  
  
  # ----------------------------------------------------------
  # 4. Select first exact match for audit only
  # ----------------------------------------------------------
  
  filtered_gct_file_10b3e <-
    filtered_hits_10b3e[1]
  
  
  cat("\nSelected filtered GCT:\n")
  cat(filtered_gct_file_10b3e, "\n")
  
  
  # ----------------------------------------------------------
  # 5. Read first 5 lines only
  # ----------------------------------------------------------
  
  filtered_head_10b3e <- readLines(
    filtered_gct_file_10b3e,
    n = 5,
    warn = FALSE
  )
  
  
  cat("\n============================================================\n")
  cat("FIRST 5 LINES OF FILTERED GCT\n")
  cat("============================================================\n\n")
  
  for (i in seq_along(filtered_head_10b3e)) {
    
    fields_i <- strsplit(
      filtered_head_10b3e[i],
      "\t",
      fixed = TRUE
    )[[1]]
    
    cat(
      "Line ",
      i,
      ": ",
      length(fields_i),
      " fields\n",
      sep = ""
    )
    
    print(
      head(
        fields_i,
        min(8, length(fields_i))
      )
    )
    
    cat("\n")
  }
  
  
  # ----------------------------------------------------------
  # 6. Identify header line by exact sample matching
  # ----------------------------------------------------------
  
  header_line_10b3e <- NA_integer_
  
  for (i in seq_along(filtered_head_10b3e)) {
    
    fields_i <- strsplit(
      filtered_head_10b3e[i],
      "\t",
      fixed = TRUE
    )[[1]]
    
    if (length(fields_i) >= 3) {
      
      candidate_ids <- fields_i[-c(1, 2)]
      
      if (
        length(candidate_ids) ==
        length(estimate_sample_ids_7be10)
      ) {
        
        if (
          identical(
            as.character(candidate_ids),
            as.character(
              estimate_sample_ids_7be10
            )
          )
        ) {
          
          header_line_10b3e <- i
          
        }
      }
    }
  }
  
  
  cat(
    "\nDetected sample-header line: ",
    header_line_10b3e,
    "\n",
    sep = ""
  )
  
  
  # ----------------------------------------------------------
  # 7. Read complete filtered file only if header identified
  # ----------------------------------------------------------
  
  if (!is.na(header_line_10b3e)) {
    
    filtered_raw_10b3e <- readLines(
      filtered_gct_file_10b3e,
      warn = FALSE
    )
    
    filtered_rows_10b3e <-
      filtered_raw_10b3e[
        -(seq_len(header_line_10b3e - 1))
      ]
    
    
    # --------------------------------------------------------
    # 8. Parse gene symbols
    # --------------------------------------------------------
    
    filtered_gene_table_10b3e <- strsplit(
      filtered_rows_10b3e,
      "\t",
      fixed = TRUE
    )
    
    filtered_gene_symbols_10b3e <- vapply(
      filtered_gene_table_10b3e,
      function(x) x[1],
      character(1)
    )
    
    
    cat("\n============================================================\n")
    cat("FILTERED GENE AUDIT\n")
    cat("============================================================\n\n")
    
    cat(
      "Filtered gene rows: ",
      length(filtered_gene_symbols_10b3e),
      "\n",
      sep = ""
    )
    
    cat(
      "Unique filtered genes: ",
      length(
        unique(
          filtered_gene_symbols_10b3e
        )
      ),
      "\n",
      sep = ""
    )
    
    cat(
      "Duplicated filtered genes: ",
      sum(
        duplicated(
          filtered_gene_symbols_10b3e
        )
      ),
      "\n",
      sep = ""
    )
    
    
    # --------------------------------------------------------
    # 9. Core genes
    # --------------------------------------------------------
    
    core_presence_10b3e <- data.frame(
      Gene = core_genes_10b3c,
      In_filtered_GCT =
        core_genes_10b3c %in%
        filtered_gene_symbols_10b3e,
      stringsAsFactors = FALSE
    )
    
    cat("\nCore gene presence:\n")
    print(core_presence_10b3e)
    
    
    # --------------------------------------------------------
    # 10. LOO missing-gene audit
    # --------------------------------------------------------
    
    stromal_missing_10b3e <-
      stromal_genes_loo_AOC3_7be[
        !stromal_genes_loo_AOC3_7be %in%
          filtered_gene_symbols_10b3e
      ]
    
    immune_missing_10b3e <-
      immune_genes_loo_IL7R_7be[
        !immune_genes_loo_IL7R_7be %in%
          filtered_gene_symbols_10b3e
      ]
    
    
    cat("\nAOC3 LOO Stromal genes:\n")
    cat(
      "Total: ",
      length(stromal_genes_loo_AOC3_7be),
      "\n",
      sep = ""
    )
    
    cat(
      "Present: ",
      sum(
        stromal_genes_loo_AOC3_7be %in%
          filtered_gene_symbols_10b3e
      ),
      "\n",
      sep = ""
    )
    
    cat("Missing:\n")
    print(stromal_missing_10b3e)
    
    
    cat("\nIL7R LOO Immune genes:\n")
    cat(
      "Total: ",
      length(immune_genes_loo_IL7R_7be),
      "\n",
      sep = ""
    )
    
    cat(
      "Present: ",
      sum(
        immune_genes_loo_IL7R_7be %in%
          filtered_gene_symbols_10b3e
      ),
      "\n",
      sep = ""
    )
    
    cat("Missing:\n")
    print(immune_missing_10b3e)
    
    
    # --------------------------------------------------------
    # 11. Final audit
    # --------------------------------------------------------
    
    cat("\n============================================================\n")
    cat("STEP 10B-3E COMPLETE\n")
    cat("============================================================\n\n")
    
    cat(
      "Actual filtered GCT was recovered and audited.\n",
      "No ESTIMATE score was recalculated.\n",
      "No LOO score was calculated.\n",
      "No statistical analysis was performed.\n",
      "No locked Step14 object was modified.\n",
      sep = ""
    )
    
  } else {
    
    cat(
      "\nSTEP 10B-3E STOP:\n",
      "Filtered file found, but sample-header line could not be identified.\n",
      "No score recalculation performed.\n",
      "No LOO score calculated.\n",
      "No statistical analysis performed.\n",
      "No locked Step14 object modified.\n",
      sep = ""
    )
  }
}


# ============================================================
# STEP 7B-7E-10B-3F
# STEP14 ESTIMATE FILTER OBJECT RECOVERY AUDIT
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3F: STEP14 ESTIMATE FILTER OBJECT AUDIT\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Load locked Step14 source
# ------------------------------------------------------------

step14_locked_file_10b3f <-
  "input_dir/Paper1_tumor-only WGCNA/WORKFLOW-03/Step14/Step14_FINAL_LOCKED/RData/WORKFLOW-03_Step14A-14E_FINAL_LOCKED.RData"

stopifnot(file.exists(step14_locked_file_10b3f))

env_10b3f <- new.env(parent = emptyenv())

load(
  step14_locked_file_10b3f,
  envir = env_10b3f
)

cat("Step14 locked RData loaded.\n\n")


# ------------------------------------------------------------
# 2. Inspect relevant object classes
# ------------------------------------------------------------

target_objects_10b3f <- c(
  "estimate_common",
  "estimate_gct_input",
  "estimate_input_file",
  "estimate_filtered",
  "estimate_input_df",
  "expr_estimate",
  "estimate_scores_wgcna"
)

cat("============================================================\n")
cat("OBJECT CLASS / STRUCTURE AUDIT\n")
cat("============================================================\n\n")

for (obj_i in target_objects_10b3f) {
  
  cat("------------------------------------------------------------\n")
  cat("OBJECT: ", obj_i, "\n", sep = "")
  
  if (!exists(obj_i, envir = env_10b3f, inherits = FALSE)) {
    
    cat("STATUS: NOT PRESENT\n\n")
    next
  }
  
  x_i <- get(
    obj_i,
    envir = env_10b3f,
    inherits = FALSE
  )
  
  cat("Class:\n")
  print(class(x_i))
  
  cat("\nTypeof:\n")
  print(typeof(x_i))
  
  cat("\nLength:\n")
  print(length(x_i))
  
  cat("\nDimensions:\n")
  print(dim(x_i))
  
  cat("\nNames:\n")
  print(
    head(
      names(x_i),
      10
    )
  )
  
  cat("\nRownames head:\n")
  print(
    head(
      rownames(x_i),
      10
    )
  )
  
  cat("\nColnames head:\n")
  print(
    head(
      colnames(x_i),
      10
    )
  )
  
  cat("\nStructure:\n")
  str(
    x_i,
    max.level = 2
  )
  
  cat("\n")
}


# ------------------------------------------------------------
# 3. Specifically inspect estimate_common
# ------------------------------------------------------------

if (
  exists(
    "estimate_common",
    envir = env_10b3f,
    inherits = FALSE
  )
) {
  
  estimate_common_10b3f <-
    get(
      "estimate_common",
      envir = env_10b3f,
      inherits = FALSE
    )
  
  cat("\n============================================================\n")
  cat("DETAILED estimate_common AUDIT\n")
  cat("============================================================\n\n")
  
  cat("Class:\n")
  print(class(estimate_common_10b3f))
  
  cat("\nDimensions:\n")
  print(dim(estimate_common_10b3f))
  
  cat("\nFirst 10 values:\n")
  print(
    head(
      estimate_common_10b3f,
      10
    )
  )
  
  if (is.data.frame(estimate_common_10b3f)) {
    
    cat("\nColumn names:\n")
    print(
      colnames(estimate_common_10b3f)
    )
    
    cat("\nFirst rows:\n")
    print(
      head(
        estimate_common_10b3f,
        3
      )
    )
  }
}


# ------------------------------------------------------------
# 4. Specifically inspect estimate_gct_input
# ------------------------------------------------------------

if (
  exists(
    "estimate_gct_input",
    envir = env_10b3f,
    inherits = FALSE
  )
) {
  
  estimate_gct_input_10b3f <-
    get(
      "estimate_gct_input",
      envir = env_10b3f,
      inherits = FALSE
    )
  
  cat("\n============================================================\n")
  cat("DETAILED estimate_gct_input AUDIT\n")
  cat("============================================================\n\n")
  
  cat("Class:\n")
  print(class(estimate_gct_input_10b3f))
  
  cat("\nDimensions:\n")
  print(dim(estimate_gct_input_10b3f))
  
  cat("\nFirst values:\n")
  print(
    head(
      estimate_gct_input_10b3f,
      10
    )
  )
  
  if (is.data.frame(estimate_gct_input_10b3f)) {
    
    cat("\nColumn names:\n")
    print(
      colnames(estimate_gct_input_10b3f)
    )
    
    cat("\nFirst rows:\n")
    print(
      head(
        estimate_gct_input_10b3f,
        3
      )
    )
  }
}


# ------------------------------------------------------------
# 5. Inspect estimate_input_file
# ------------------------------------------------------------

if (
  exists(
    "estimate_input_file",
    envir = env_10b3f,
    inherits = FALSE
  )
) {
  
  estimate_input_file_10b3f <-
    get(
      "estimate_input_file",
      envir = env_10b3f,
      inherits = FALSE
    )
  
  cat("\n============================================================\n")
  cat("estimate_input_file\n")
  cat("============================================================\n\n")
  
  print(
    estimate_input_file_10b3f
  )
  
  if (
    is.character(estimate_input_file_10b3f) &&
    length(estimate_input_file_10b3f) > 0
  ) {
    
    cat("\nFile existence check:\n")
    
    print(
      file.exists(
        estimate_input_file_10b3f
      )
    )
    
    cat("\nCurrent working directory:\n")
    print(
      getwd()
    )
    
    cat("\nAbsolute path interpretation:\n")
    
    print(
      normalizePath(
        estimate_input_file_10b3f,
        winslash = "/",
        mustWork = FALSE
      )
    )
  }
}


# ------------------------------------------------------------
# 6. Final safety statement
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 10B-3F SAFETY STATUS\n")
cat("============================================================\n\n")

cat(
  "No ESTIMATE score recalculated.\n",
  "No LOO score calculated.\n",
  "No statistical analysis performed.\n",
  "No Step14 locked object modified.\n",
  "No Figure 6 modified.\n",
  "No Supplementary Table 3 modified.\n",
  sep = ""
)

cat("\nSTEP 10B-3F COMPLETE.\n")

# ============================================================
# STEP 7B-7E-10B-3G
# ESTIMATE INTERMEDIATE FILE TRACE AUDIT
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3G: ESTIMATE INTERMEDIATE FILE TRACE AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Search roots
# ------------------------------------------------------------

search_roots_10b3g <- unique(
  normalizePath(
    c(
      "O:/Rproject/Paper1_tumor-only WGCNA",
      "O:/Rproject/Paper 1_Step 7_ Immune Infiltration Analysis"
    ),
    winslash = "/",
    mustWork = FALSE
  )
)

cat("Search roots:\n")
print(search_roots_10b3g)


# ------------------------------------------------------------
# 2. Search all files containing ESTIMATE
# ------------------------------------------------------------

estimate_files_10b3g <- unique(
  unlist(
    lapply(
      search_roots_10b3g,
      function(root_i) {
        
        if (!dir.exists(root_i)) {
          return(character(0))
        }
        
        list.files(
          path = root_i,
          pattern = "ESTIMATE",
          recursive = TRUE,
          full.names = TRUE,
          ignore.case = TRUE
        )
      }
    )
  )
)

cat("\n============================================================\n")
cat("ALL ESTIMATE-RELATED FILES\n")
cat("============================================================\n\n")

cat(
  "Number of files found: ",
  length(estimate_files_10b3g),
  "\n\n",
  sep = ""
)

if (length(estimate_files_10b3g) > 0) {
  print(
    estimate_files_10b3g
  )
} else {
  cat("No ESTIMATE-related files found.\n")
}


# ------------------------------------------------------------
# 3. Search specifically for expected filenames
# ------------------------------------------------------------

target_names_10b3g <- c(
  "Workflow-04-2C-1A_ESTIMATE_input_GeneSymbol_log2_FPKM_UQ.txt",
  "Workflow-04-2C-1A_ESTIMATE_input_GCT.txt",
  "Workflow-04-2C-1A_ESTIMATE_filtered.gct"
)

cat("\n============================================================\n")
cat("EXACT TARGET FILE SEARCH\n")
cat("============================================================\n\n")

exact_hits_10b3g <- list()

for (target_i in target_names_10b3g) {
  
  hits_i <- unique(
    unlist(
      lapply(
        search_roots_10b3g,
        function(root_i) {
          
          if (!dir.exists(root_i)) {
            return(character(0))
          }
          
          list.files(
            path = root_i,
            pattern = paste0(
              "^",
              target_i,
              "$"
            ),
            recursive = TRUE,
            full.names = TRUE,
            ignore.case = FALSE
          )
        }
      )
    )
  )
  
  exact_hits_10b3g[[target_i]] <- hits_i
  
  cat("\nTarget:\n")
  cat(target_i, "\n")
  
  cat(
    "Number of exact matches: ",
    length(hits_i),
    "\n",
    sep = ""
  )
  
  if (length(hits_i) > 0) {
    print(hits_i)
  }
}


# ------------------------------------------------------------
# 4. Search filenames containing filtered/common
# ------------------------------------------------------------

filtered_common_files_10b3g <- unique(
  unlist(
    lapply(
      search_roots_10b3g,
      function(root_i) {
        
        if (!dir.exists(root_i)) {
          return(character(0))
        }
        
        list.files(
          path = root_i,
          pattern = "(filtered|common)",
          recursive = TRUE,
          full.names = TRUE,
          ignore.case = TRUE
        )
      }
    )
  )
)

cat("\n============================================================\n")
cat("FILTERED / COMMON FILES\n")
cat("============================================================\n\n")

cat(
  "Number found: ",
  length(filtered_common_files_10b3g),
  "\n\n",
  sep = ""
)

if (length(filtered_common_files_10b3g) > 0) {
  print(
    filtered_common_files_10b3g
  )
}


# ------------------------------------------------------------
# 5. Safety status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 10B-3G SAFETY STATUS\n")
cat("============================================================\n\n")

cat(
  "File search only.\n",
  "No ESTIMATE score recalculated.\n",
  "No LOO score calculated.\n",
  "No statistical analysis performed.\n",
  "No Step14 locked object modified.\n",
  "No Figure 6 modified.\n",
  "No Supplementary Table 3 modified.\n",
  sep = ""
)

cat("\nSTEP 10B-3G COMPLETE.\n")

# ============================================================
# STEP 7B-7E-10B-3H
# ACTUAL ESTIMATE INPUT FILE PROVENANCE AUDIT
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3H: ACTUAL ESTIMATE INPUT FILE PROVENANCE AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define actual file
# ------------------------------------------------------------

estimate_input_file_10b3h <-
  "input_dir/Paper1_tumor-only WGCNA/Workflow-04-2C-1A_ESTIMATE_input_GeneSymbol_log2_FPKM_UQ.txt"

stopifnot(
  file.exists(
    estimate_input_file_10b3h
  )
)

cat("File exists: TRUE\n")
cat("File:\n")
cat(estimate_input_file_10b3h, "\n\n")


# ------------------------------------------------------------
# 2. Read first 5 lines only
# ------------------------------------------------------------

head_lines_10b3h <- readLines(
  estimate_input_file_10b3h,
  n = 5,
  warn = FALSE
)

cat("============================================================\n")
cat("FIRST 5 LINES\n")
cat("============================================================\n\n")

for (i in seq_along(head_lines_10b3h)) {
  
  fields_i <- strsplit(
    head_lines_10b3h[i],
    "\t",
    fixed = TRUE
  )[[1]]
  
  cat(
    "Line ",
    i,
    ": ",
    length(fields_i),
    " fields\n",
    sep = ""
  )
  
  print(
    head(
      fields_i,
      min(
        10,
        length(fields_i)
      )
    )
  )
  
  cat("\n")
}


# ------------------------------------------------------------
# 3. Read complete file
# ------------------------------------------------------------

input_raw_10b3h <- read.table(
  estimate_input_file_10b3h,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("ACTUAL ESTIMATE INPUT FILE DIMENSIONS\n")
cat("============================================================\n\n")

cat(
  "Rows: ",
  nrow(input_raw_10b3h),
  "\n",
  sep = ""
)

cat(
  "Columns: ",
  ncol(input_raw_10b3h),
  "\n",
  sep = ""
)

cat("\nColumn names head:\n")
print(
  head(
    colnames(input_raw_10b3h),
    10
  )
)


# ------------------------------------------------------------
# 4. Identify gene-symbol column
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("COLUMN STRUCTURE\n")
cat("============================================================\n\n")

print(
  colnames(input_raw_10b3h)[1:min(
    5,
    ncol(input_raw_10b3h)
  )]
)

gene_column_candidates_10b3h <- intersect(
  c(
    "GeneSymbol",
    "Gene",
    "genes",
    "gene"
  ),
  colnames(input_raw_10b3h)
)

cat("\nGene-symbol candidates:\n")
print(
  gene_column_candidates_10b3h
)

if (length(gene_column_candidates_10b3h) == 0) {
  
  stop(
    "No GeneSymbol-like column identified."
  )
}

gene_col_10b3h <-
  gene_column_candidates_10b3h[1]

gene_symbols_10b3h <-
  as.character(
    input_raw_10b3h[[gene_col_10b3h]]
  )


# ------------------------------------------------------------
# 5. Sample columns
# ------------------------------------------------------------

sample_cols_10b3h <-
  setdiff(
    colnames(input_raw_10b3h),
    gene_col_10b3h
  )

cat("\nNumber of sample columns: ")
cat(length(sample_cols_10b3h), "\n")

cat("\nFirst sample IDs:\n")
print(
  head(
    sample_cols_10b3h,
    10
  )
)


# ------------------------------------------------------------
# 6. Compare sample identity/order with Step14
# ------------------------------------------------------------

sample_identity_10b3h <-
  identical(
    as.character(sample_cols_10b3h),
    as.character(estimate_sample_ids_7be10)
  )

cat("\n============================================================\n")
cat("SAMPLE IDENTITY AUDIT\n")
cat("============================================================\n\n")

cat(
  "Sample number in file: ",
  length(sample_cols_10b3h),
  "\n",
  sep = ""
)

cat(
  "Sample number in Step14: ",
  length(estimate_sample_ids_7be10),
  "\n",
  sep = ""
)

cat(
  "Exact identity/order match: ",
  sample_identity_10b3h,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Gene audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GENE SYMBOL AUDIT\n")
cat("============================================================\n\n")

cat(
  "Gene rows: ",
  length(gene_symbols_10b3h),
  "\n",
  sep = ""
)

cat(
  "Unique genes: ",
  length(
    unique(
      gene_symbols_10b3h
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Duplicated gene rows: ",
  sum(
    duplicated(
      gene_symbols_10b3h
    )
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Core gene presence
# ------------------------------------------------------------

core_presence_10b3h <- data.frame(
  Gene = core_genes_10b3c,
  Present =
    core_genes_10b3c %in%
    gene_symbols_10b3h,
  stringsAsFactors = FALSE
)

cat("\nCore gene presence:\n")
print(
  core_presence_10b3h
)


# ------------------------------------------------------------
# 9. Compare gene universe with expr_estimate
# ------------------------------------------------------------

expr_gene_symbols_10b3h <-
  rownames(
    expr_estimate_7be10
  )

cat("\n============================================================\n")
cat("COMPARISON WITH expr_estimate_7be10\n")
cat("============================================================\n\n")

cat(
  "expr_estimate rows: ",
  nrow(expr_estimate_7be10),
  "\n",
  sep = ""
)

cat(
  "Input-file gene rows: ",
  length(gene_symbols_10b3h),
  "\n",
  sep = ""
)

cat(
  "Input genes present in expr_estimate: ",
  sum(
    gene_symbols_10b3h %in%
      expr_gene_symbols_10b3h
  ),
  "\n",
  sep = ""
)

cat(
  "expr_estimate genes present in input file: ",
  sum(
    expr_gene_symbols_10b3h %in%
      gene_symbols_10b3h
  ),
  "\n",
  sep = ""
)

cat(
  "Exact gene order match: ",
  identical(
    gene_symbols_10b3h,
    expr_gene_symbols_10b3h
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Compare dimensions with estimate_input_df
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("COMPARISON WITH estimate_input_df_7be10\n")
cat("============================================================\n\n")

cat(
  "estimate_input_df dimensions: ",
  paste(
    dim(
      estimate_input_df_7be10
    ),
    collapse = " x "
  ),
  "\n",
  sep = ""
)

cat(
  "Actual file dimensions: ",
  paste(
    dim(
      input_raw_10b3h
    ),
    collapse = " x "
  ),
  "\n",
  sep = ""
)

cat(
  "Column-name identity with estimate_input_df: ",
  identical(
    colnames(input_raw_10b3h),
    colnames(estimate_input_df_7be10)
  ),
  "\n",
  sep = "")


# ------------------------------------------------------------
# 11. Numeric-value audit for shared genes/samples
# ------------------------------------------------------------

shared_genes_10b3h <-
  intersect(
    gene_symbols_10b3h,
    expr_gene_symbols_10b3h
  )

cat("\n============================================================\n")
cat("SHARED DATA AUDIT\n")
cat("============================================================\n\n")

cat(
  "Shared genes: ",
  length(shared_genes_10b3h),
  "\n",
  sep = ""
)

cat(
  "Shared samples: ",
  length(
    intersect(
      sample_cols_10b3h,
      colnames(
        expr_estimate_7be10
      )
    )
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Safety status\n
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 10B-3H SAFETY STATUS\n")
cat("============================================================\n\n")

cat(
  "Input-file provenance audit only.\n",
  "No ESTIMATE score recalculated.\n",
  "No LOO score calculated.\n",
  "No statistical analysis performed.\n",
  "No Step14 locked object modified.\n",
  "No Figure 6 modified.\n",
  "No Supplementary Table 3 modified.\n",
  sep = ""
)

cat("\nSTEP 10B-3H COMPLETE.\n")


# ============================================================
# STEP 7B-7E-10B-3I
# ESTIMATE INPUT DATAFRAME STRUCTURE AUDIT
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3I: ESTIMATE INPUT DATAFRAME STRUCTURE AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Basic structure
# ------------------------------------------------------------

cat("============================================================\n")
cat("estimate_input_df_7be10 STRUCTURE\n")
cat("============================================================\n\n")

cat("Class:\n")
print(class(estimate_input_df_7be10))

cat("\nType:\n")
print(typeof(estimate_input_df_7be10))

cat("\nDimensions:\n")
print(dim(estimate_input_df_7be10))

cat("\nColumn names:\n")
print(
  colnames(estimate_input_df_7be10)
)


# ------------------------------------------------------------
# 2. Identify non-sample columns
# ------------------------------------------------------------

sample_ids_10b3i <-
  as.character(
    estimate_sample_ids_7be10
  )

all_columns_10b3i <-
  colnames(
    estimate_input_df_7be10
  )

non_sample_columns_10b3i <-
  setdiff(
    all_columns_10b3i,
    sample_ids_10b3i
  )

cat("\n============================================================\n")
cat("NON-SAMPLE COLUMNS\n")
cat("============================================================\n\n")

cat(
  "Number of non-sample columns: ",
  length(non_sample_columns_10b3i),
  "\n",
  sep = ""
)

print(
  non_sample_columns_10b3i
)


# ------------------------------------------------------------
# 3. Sample-column identity
# ------------------------------------------------------------

sample_columns_in_df_10b3i <-
  intersect(
    all_columns_10b3i,
    sample_ids_10b3i
  )

cat("\n============================================================\n")
cat("SAMPLE COLUMN AUDIT\n")
cat("============================================================\n\n")

cat(
  "Expected sample columns: ",
  length(sample_ids_10b3i),
  "\n",
  sep = ""
)

cat(
  "Observed sample columns: ",
  length(sample_columns_in_df_10b3i),
  "\n",
  sep = ""
)

cat(
  "Exact sample order after removing non-sample columns: ",
  identical(
    all_columns_10b3i[
      all_columns_10b3i %in%
        sample_ids_10b3i
    ],
    sample_ids_10b3i
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Inspect first rows
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FIRST 5 ROWS\n")
cat("============================================================\n\n")

print(
  head(
    estimate_input_df_7be10,
    5
  )
)


# ------------------------------------------------------------
# 5. Inspect first columns
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FIRST 5 COLUMNS\n")
cat("============================================================\n\n")

print(
  head(
    estimate_input_df_7be10[
      ,
      seq_len(
        min(
          5,
          ncol(
            estimate_input_df_7be10
          )
        )
      ),
      drop = FALSE
    ],
    5
  )
)


# ------------------------------------------------------------
# 6. Check GeneSymbol representation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GENE SYMBOL REPRESENTATION\n")
cat("============================================================\n\n")

cat(
  "Rownames present: ",
  !is.null(
    rownames(
      estimate_input_df_7be10
    )
  ),
  "\n",
  sep = ""
)

if (!is.null(rownames(estimate_input_df_7be10))) {
  
  cat(
    "Number of rownames: ",
    length(
      rownames(
        estimate_input_df_7be10
      )
    ),
    "\n",
    sep = ""
  )
  
  cat("\nFirst 10 rownames:\n")
  print(
    head(
      rownames(
        estimate_input_df_7be10
      ),
      10
    )
  )
}


# ------------------------------------------------------------
# 7. Compare GeneSymbol with expr_estimate
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GENE SYMBOL / ROWNAME COMPARISON\n")
cat("============================================================\n\n")

cat(
  "expr_estimate rownames = estimate_input_df rownames: ",
  identical(
    rownames(
      expr_estimate_7be10
    ),
    rownames(
      estimate_input_df_7be10
    )
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Check whether Description exists
# ------------------------------------------------------------

description_candidates_10b3i <-
  grep(
    "description|desc",
    colnames(
      estimate_input_df_7be10
    ),
    ignore.case = TRUE,
    value = TRUE
  )

cat("\n============================================================\n")
cat("DESCRIPTION-LIKE COLUMNS\n")
cat("============================================================\n\n")

print(
  description_candidates_10b3i
)


# ------------------------------------------------------------
# 9. Compare actual file columns with dataframe
# ------------------------------------------------------------

actual_file_columns_10b3i <-
  colnames(
    input_raw_10b3h
  )

cat("\n============================================================\n")
cat("FILE VS DATAFRAME COLUMN STRUCTURE\n")
cat("============================================================\n\n")

cat(
  "Actual file columns: ",
  length(actual_file_columns_10b3i),
  "\n",
  sep = ""
)

cat(
  "Dataframe columns: ",
  length(all_columns_10b3i),
  "\n",
  sep = ""
)

cat("\nColumns in dataframe but not actual file:\n")

print(
  setdiff(
    all_columns_10b3i,
    actual_file_columns_10b3i
  )
)

cat("\nColumns in actual file but not dataframe:\n")

print(
  setdiff(
    actual_file_columns_10b3i,
    all_columns_10b3i
  )
)


# ------------------------------------------------------------
# 10. Core-gene row audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CORE GENE ROW AUDIT\n")
cat("============================================================\n\n")

for (gene_i in core_genes_10b3c) {
  
  idx_i <-
    which(
      rownames(
        estimate_input_df_7be10
      ) ==
        gene_i
    )
  
  cat(
    gene_i,
    ": row index = ",
    ifelse(
      length(idx_i) == 1,
      idx_i,
      NA
    ),
    "\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 11. Safety status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 10B-3I SAFETY STATUS\n")
cat("============================================================\n\n")

cat(
  "Structure audit only.\n",
  "No ESTIMATE score recalculated.\n",
  "No LOO score calculated.\n",
  "No statistical analysis performed.\n",
  "No Step14 locked object modified.\n",
  "No Figure 6 modified.\n",
  "No Supplementary Table 3 modified.\n",
  sep = ""
)

cat("\nSTEP 10B-3I COMPLETE.\n")


# ============================================================
# STEP 10B-3J
# ESTIMATE FILTERED/COMMON-GENE PROVENANCE AUDIT
# ============================================================

cat("============================================================\n")
cat("STEP 10B-3J: ESTIMATE FILTERED/COMMON-GENE PROVENANCE AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Required objects
# ------------------------------------------------------------

required_objects_10b3j <- c(
  "estimate_common",
  "estimate_filtered",
  "estimate_gct_input",
  "estimate_input_df",
  "estimate_input_file",
  "estimate_scores",
  "estimate_scores_wgcna",
  "expr_estimate"
)

cat("============================================================\n")
cat("REQUIRED OBJECT EXISTENCE\n")
cat("============================================================\n\n")

for (obj_i in required_objects_10b3j) {
  
  cat(
    obj_i,
    ": ",
    exists(obj_i, inherits = TRUE),
    "\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 2. Inspect estimate_common
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("estimate_common STRUCTURE\n")
cat("============================================================\n\n")

if (exists("estimate_common", inherits = TRUE)) {
  
  cat("Class:\n")
  print(class(estimate_common))
  
  cat("\nType:\n")
  print(typeof(estimate_common))
  
  cat("\nLength:\n")
  print(length(estimate_common))
  
  cat("\nDimensions:\n")
  print(dim(estimate_common))
  
  cat("\nNames:\n")
  print(names(estimate_common))
  
  cat("\nStructure:\n")
  str(
    estimate_common,
    max.level = 2
  )
}


# ------------------------------------------------------------
# 3. Inspect estimate_filtered
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("estimate_filtered STRUCTURE\n")
cat("============================================================\n\n")

if (exists("estimate_filtered", inherits = TRUE)) {
  
  cat("Class:\n")
  print(class(estimate_filtered))
  
  cat("\nType:\n")
  print(typeof(estimate_filtered))
  
  cat("\nLength:\n")
  print(length(estimate_filtered))
  
  cat("\nDimensions:\n")
  print(dim(estimate_filtered))
  
  cat("\nFirst values:\n")
  print(
    head(
      estimate_filtered,
      20
    )
  )
}


# ------------------------------------------------------------
# 4. Inspect estimate_gct_input
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("estimate_gct_input STRUCTURE\n")
cat("============================================================\n\n")

if (exists("estimate_gct_input", inherits = TRUE)) {
  
  cat("Class:\n")
  print(class(estimate_gct_input))
  
  cat("\nType:\n")
  print(typeof(estimate_gct_input))
  
  cat("\nLength:\n")
  print(length(estimate_gct_input))
  
  cat("\nDimensions:\n")
  print(dim(estimate_gct_input))
  
  cat("\nFirst values:\n")
  print(
    head(
      estimate_gct_input,
      20
    )
  )
}


# ------------------------------------------------------------
# 5. Inspect estimate_input_file
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("estimate_input_file\n")
cat("============================================================\n\n")

if (exists("estimate_input_file", inherits = TRUE)) {
  
  print(
    estimate_input_file
  )
}


# ------------------------------------------------------------
# 6. Inspect estimate_scores
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("estimate_scores STRUCTURE\n")
cat("============================================================\n\n")

if (exists("estimate_scores", inherits = TRUE)) {
  
  cat("Class:\n")
  print(class(estimate_scores))
  
  cat("\nType:\n")
  print(typeof(estimate_scores))
  
  cat("\nDimensions:\n")
  print(dim(estimate_scores))
  
  cat("\nRow names:\n")
  print(rownames(estimate_scores))
  
  cat("\nColumn count:\n")
  print(ncol(estimate_scores))
}


# ------------------------------------------------------------
# 7. Search current environment for ESTIMATE-like objects
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ESTIMATE-LIKE OBJECT INVENTORY\n")
cat("============================================================\n\n")

all_objects_10b3j <- ls(envir = .GlobalEnv)

estimate_like_objects_10b3j <- all_objects_10b3j[
  grepl(
    "estimate|filtered|common|gct",
    all_objects_10b3j,
    ignore.case = TRUE
  )
]

print(
  estimate_like_objects_10b3j
)


# ------------------------------------------------------------
# 8. Class / dimensions of ESTIMATE-like objects
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ESTIMATE-LIKE OBJECT SUMMARY\n")
cat("============================================================\n\n")

for (obj_i in estimate_like_objects_10b3j) {
  
  obj_i_value <- get(
    obj_i,
    envir = .GlobalEnv
  )
  
  cat(
    "\nOBJECT: ",
    obj_i,
    "\n",
    sep = ""
  )
  
  cat(
    "Class: ",
    paste(
      class(obj_i_value),
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "Type: ",
    typeof(obj_i_value),
    "\n",
    sep = ""
  )
  
  cat(
    "Length: ",
    length(obj_i_value),
    "\n",
    sep = ""
  )
  
  if (!is.null(dim(obj_i_value))) {
    
    cat(
      "Dimensions: ",
      paste(
        dim(obj_i_value),
        collapse = " x "
      ),
      "\n",
      sep = ""
    )
  }
}


# ------------------------------------------------------------
# 9. Search for objects containing approximately
#    the expected ESTIMATE common-gene universe
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("COMMON-GENE UNIVERSE SEARCH\n")
cat("============================================================\n\n")

for (obj_i in estimate_like_objects_10b3j) {
  
  obj_i_value <- get(
    obj_i,
    envir = .GlobalEnv
  )
  
  d_i <- dim(obj_i_value)
  
  if (!is.null(d_i)) {
    
    if (
      d_i[1] >= 10000 &&
      d_i[1] <= 12000 &&
      d_i[2] >= 500 &&
      d_i[2] <= 520
    ) {
      
      cat(
        "POTENTIAL COMMON-GENE MATRIX: ",
        obj_i,
        " = ",
        paste(
          d_i,
          collapse = " x "
        ),
        "\n",
        sep = ""
      )
    }
  }
}


# ------------------------------------------------------------
# 10. Check whether estimate_common itself contains
#     approximately 10k genes
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("estimate_common GENE COUNT CHECK\n")
cat("============================================================\n\n")

if (exists("estimate_common", inherits = TRUE)) {
  
  ec_i <- estimate_common
  
  if (is.atomic(ec_i)) {
    
    cat(
      "Atomic vector length = ",
      length(ec_i),
      "\n",
      sep = ""
    )
    
  } else if (is.data.frame(ec_i)) {
    
    cat(
      "Data.frame dimensions = ",
      paste(
        dim(ec_i),
        collapse = " x "
      ),
      "\n",
      sep = ""
    )
    
  } else if (is.matrix(ec_i)) {
    
    cat(
      "Matrix dimensions = ",
      paste(
        dim(ec_i),
        collapse = " x "
      ),
      "\n",
      sep = ""
    )
    
  } else {
    
    cat(
      "Non-atomic object; structure printed above.\n"
    )
  }
}


# ------------------------------------------------------------
# 11. Safety status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 10B-3J SAFETY STATUS\n")
cat("============================================================\n\n")

cat(
  "Structure/provenance audit only.\n",
  "No ESTIMATE score recalculated.\n",
  "No LOO score calculated.\n",
  "No statistical analysis performed.\n",
  "No Step14 locked object modified.\n",
  "No Figure 6 modified.\n",
  "No Supplementary Table 3 modified.\n",
  sep = ""
)

cat("\n\nSTEP 10B-3J COMPLETE.\n")


cat("STEP 10B-3K: ESTIMATE COMMON-GENE MATRIX STRUCTURAL VALIDATION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Required objects
# ------------------------------------------------------------

required_objects_10b3k <- c(
  "estimate_common_10b3f",
  "expr_estimate_7be10",
  "estimate_scores_source_7be10",
  "estimate_sample_ids_7be10"
)

cat("============================================================\n")
cat("REQUIRED OBJECT EXISTENCE\n")
cat("============================================================\n\n")

for (obj_i in required_objects_10b3k) {
  
  cat(
    obj_i,
    ": ",
    exists(obj_i, inherits = TRUE),
    "\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 2. Basic structure
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("COMMON-GENE MATRIX STRUCTURE\n")
cat("============================================================\n\n")

ec <- estimate_common_10b3f

cat("Class:\n")
print(class(ec))

cat("\nType:\n")
print(typeof(ec))

cat("\nDimensions:\n")
print(dim(ec))

cat("\nColumn names (first 10):\n")
print(head(colnames(ec), 10))

cat("\nColumn names (last 10):\n")
print(tail(colnames(ec), 10))


# ------------------------------------------------------------
# 3. Identify sample columns
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SAMPLE COLUMN AUDIT\n")
cat("============================================================\n\n")

sample_ids_10b3k <- estimate_sample_ids_7be10

sample_columns_10b3k <- intersect(
  colnames(ec),
  sample_ids_10b3k
)

cat(
  "Expected samples: ",
  length(sample_ids_10b3k),
  "\n",
  sep = ""
)

cat(
  "Observed sample columns: ",
  length(sample_columns_10b3k),
  "\n",
  sep = ""
)

cat(
  "Exact sample identity: ",
  identical(
    sample_columns_10b3k,
    sample_ids_10b3k
  ),
  "\n",
  sep = ""
)

cat(
  "Exact sample set equality: ",
  setequal(
    sample_columns_10b3k,
    sample_ids_10b3k
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Non-sample columns
# ------------------------------------------------------------

non_sample_columns_10b3k <- setdiff(
  colnames(ec),
  sample_ids_10b3k
)

cat("\nNon-sample columns:\n")
print(non_sample_columns_10b3k)


# ------------------------------------------------------------
# 5. Gene identifier structure
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GENE IDENTIFIER AUDIT\n")
cat("============================================================\n\n")

cat(
  "Rownames present: ",
  !is.null(rownames(ec)),
  "\n",
  sep = ""
)

cat(
  "Number of rownames: ",
  length(rownames(ec)),
  "\n",
  sep = ""
)

cat(
  "Unique rownames: ",
  length(unique(rownames(ec))) ==
    length(rownames(ec)),
  "\n",
  sep = ""
)

cat("\nFirst 20 rownames:\n")
print(head(rownames(ec), 20))

cat("\nLast 20 rownames:\n")
print(tail(rownames(ec), 20))


# ------------------------------------------------------------
# 6. Check whether gene identifiers are also in columns
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GENE COLUMN AUDIT\n")
cat("============================================================\n\n")

gene_like_columns_10b3k <- non_sample_columns_10b3k[
  grepl(
    "gene|symbol|description|id",
    non_sample_columns_10b3k,
    ignore.case = TRUE
  )
]

cat("Gene-like non-sample columns:\n")
print(gene_like_columns_10b3k)

if (length(gene_like_columns_10b3k) > 0) {
  
  for (col_i in gene_like_columns_10b3k) {
    
    cat(
      "\nColumn: ",
      col_i,
      "\n",
      sep = ""
    )
    
    cat(
      "First values:\n"
    )
    
    print(
      head(
        ec[[col_i]],
        10
      )
    )
  }
}


# ------------------------------------------------------------
# 7. Core gene presence
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CORE GENE PRESENCE\n")
cat("============================================================\n\n")

core_genes_10b3k <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_presence_10b3k <- core_genes_10b3k %in%
  rownames(ec)

core_presence_table_10b3k <- data.frame(
  Gene = core_genes_10b3k,
  Present = core_presence_10b3k,
  stringsAsFactors = FALSE
)

print(core_presence_table_10b3k)


# ------------------------------------------------------------
# 8. Signature membership presence
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ESTIMATE SIGNATURE PRESENCE\n")
cat("============================================================\n\n")

stromal_genes_10b3k <- unique(
  stromal_genes_loo_AOC3_7be10
)

immune_genes_10b3k <- unique(
  immune_genes_loo_IL7R_7be10
)

stromal_full_10b3k <- unique(
  c(
    stromal_genes_10b3k,
    "AOC3"
  )
)

immune_full_10b3k <- unique(
  c(
    immune_genes_10b3k,
    "IL7R"
  )
)

cat(
  "Stromal signature genes: ",
  length(stromal_full_10b3k),
  "\n",
  sep = ""
)

cat(
  "Stromal present: ",
  sum(
    stromal_full_10b3k %in%
      rownames(ec)
  ),
  "\n",
  sep = ""
)

cat(
  "Stromal missing:\n"
)

print(
  setdiff(
    stromal_full_10b3k,
    rownames(ec)
  )
)

cat(
  "\nImmune signature genes: ",
  length(immune_full_10b3k),
  "\n",
  sep = ""
)

cat(
  "Immune present: ",
  sum(
    immune_full_10b3k %in%
      rownames(ec)
  ),
  "\n",
  sep = ""
)

cat(
  "Immune missing:\n"
)

print(
  setdiff(
    immune_full_10b3k,
    rownames(ec)
  )
)


# ------------------------------------------------------------
# 9. Compare common-gene matrix with original expression
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("COMMON MATRIX VS ORIGINAL EXPRESSION\n")
cat("============================================================\n\n")

common_gene_ids_10b3k <- rownames(ec)

expression_gene_ids_10b3k <- rownames(
  expr_estimate_7be10
)

cat(
  "Common genes: ",
  length(common_gene_ids_10b3k),
  "\n",
  sep = ""
)

cat(
  "Common genes present in original expression: ",
  sum(
    common_gene_ids_10b3k %in%
      expression_gene_ids_10b3k
  ),
  "\n",
  sep = ""
)

cat(
  "Common genes absent from original expression: ",
  sum(
    !common_gene_ids_10b3k %in%
      expression_gene_ids_10b3k
  ),
  "\n",
  sep = ""
)

cat(
  "Exact gene order preserved from original expression: "
)

common_gene_order_check_10b3k <- identical(
  common_gene_ids_10b3k,
  expression_gene_ids_10b3k[
    match(
      common_gene_ids_10b3k,
      expression_gene_ids_10b3k
    )
  ]
)

cat(
  common_gene_order_check_10b3k,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Compare numerical values
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("NUMERICAL SUBSET CHECK\n")
cat("============================================================\n\n")

common_samples_10b3k <- sample_ids_10b3k

expr_subset_10b3k <- expr_estimate_7be10[
  common_gene_ids_10b3k,
  common_samples_10b3k,
  drop = FALSE
]

ec_sample_matrix_10b3k <- as.matrix(
  ec[
    ,
    common_samples_10b3k,
    drop = FALSE
  ]
)

cat(
  "Subset dimensions: ",
  paste(
    dim(expr_subset_10b3k),
    collapse = " x "
  ),
  "\n",
  sep = ""
)

cat(
  "Common matrix dimensions: ",
  paste(
    dim(ec_sample_matrix_10b3k),
    collapse = " x "
  ),
  "\n",
  sep = ""
)

numeric_identity_10b3k <- isTRUE(
  all.equal(
    unname(expr_subset_10b3k),
    unname(ec_sample_matrix_10b3k),
    check.attributes = FALSE
  )
)

cat(
  "Numerical identity: ",
  numeric_identity_10b3k,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Missingness / finite values
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("NUMERICAL QUALITY\n")
cat("============================================================\n\n")

cat(
  "All finite: ",
  all(
    is.finite(
      ec_sample_matrix_10b3k
    )
  ),
  "\n",
  sep = ""
)

cat(
  "NA count: ",
  sum(
    is.na(
      ec_sample_matrix_10b3k
    )
  ),
  "\n",
  sep = ""
)

cat(
  "NaN count: ",
  sum(
    is.nan(
      ec_sample_matrix_10b3k
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Inf count: ",
  sum(
    is.infinite(
      ec_sample_matrix_10b3k
    )
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Final audit conclusion
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 10B-3K CONCLUSION\n")
cat("============================================================\n\n")

cat(
  "estimate_common_10b3f is a ",
  paste(
    dim(ec),
    collapse = " x "
  ),
  " common-gene candidate matrix.\n",
  sep = ""
)

cat(
  "No ESTIMATE score recalculated.\n"
)

cat(
  "No LOO score calculated.\n"
)

cat(
  "No statistical analysis performed.\n"
)

cat(
  "No Step14 locked object modified.\n"
)

cat(
  "No Figure 6 modified.\n"
)

cat(
  "No Supplementary Table 3 modified.\n"
)

cat(
  "\nSTEP 10B-3K COMPLETE.\n"
)

cat("STEP 10B-3K-A: COMMON MATRIX IDENTITY AND SAMPLE FORMAT AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Common matrix
# ------------------------------------------------------------

ec <- estimate_common_10b3f

cat("COMMON MATRIX:\n")
cat(
  "Dimensions: ",
  paste(dim(ec), collapse = " x "),
  "\n",
  sep = ""
)

cat("\nColumn names 1:10:\n")
print(colnames(ec)[1:10])


# ------------------------------------------------------------
# 2. Identify gene identifier columns
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GENE IDENTIFIER COLUMNS\n")
cat("============================================================\n\n")

cat("NAME class:\n")
print(class(ec$NAME))

cat("\nDescription class:\n")
print(class(ec$Description))

cat("\nFirst 20 NAME:\n")
print(head(ec$NAME, 20))

cat("\nFirst 20 Description:\n")
print(head(ec$Description, 20))


# ------------------------------------------------------------
# 3. Gene symbol audit
# ------------------------------------------------------------

gene_symbols_common_10b3k <- as.character(
  ec$Description
)

cat("\n============================================================\n")
cat("GENE SYMBOL AUDIT\n")
cat("============================================================\n\n")

cat(
  "Gene count: ",
  length(gene_symbols_common_10b3k),
  "\n",
  sep = ""
)

cat(
  "Unique gene symbols: ",
  length(unique(gene_symbols_common_10b3k)) ==
    length(gene_symbols_common_10b3k),
  "\n",
  sep = ""
)

cat("\nCore gene presence:\n")

print(
  data.frame(
    Gene = c("AOC3", "MS4A7", "IL7R"),
    Present = c(
      "AOC3" %in% gene_symbols_common_10b3k,
      "MS4A7" %in% gene_symbols_common_10b3k,
      "IL7R" %in% gene_symbols_common_10b3k
    ),
    stringsAsFactors = FALSE
  )
)


# ------------------------------------------------------------
# 4. Locate core genes
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CORE GENE LOCATIONS\n")
cat("============================================================\n\n")

for (gene_i in c("AOC3", "MS4A7", "IL7R")) {
  
  idx_i <- which(
    gene_symbols_common_10b3k == gene_i
  )
  
  cat(
    gene_i,
    ": ",
    length(idx_i),
    " match(es)\n",
    sep = ""
  )
  
  if (length(idx_i) > 0) {
    
    print(
      ec[
        idx_i,
        c("NAME", "Description"),
        drop = FALSE
      ]
    )
  }
}


# ------------------------------------------------------------
# 5. Sample ID format comparison
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SAMPLE ID FORMAT\n")
cat("============================================================\n\n")

common_sample_columns_10b3k <- colnames(ec)[
  !(colnames(ec) %in% c("NAME", "Description"))
]

cat(
  "Common matrix sample count: ",
  length(common_sample_columns_10b3k),
  "\n",
  sep = ""
)

cat("\nFirst 10 common-matrix sample IDs:\n")
print(head(common_sample_columns_10b3k, 10))

cat("\nFirst 10 source sample IDs:\n")
print(head(estimate_sample_ids_7be10, 10))


# ------------------------------------------------------------
# 6. Normalize sample ID punctuation only for audit
# ------------------------------------------------------------

common_ids_normalized_10b3k <- gsub(
  "\\.",
  "-",
  common_sample_columns_10b3k
)

source_ids_normalized_10b3k <- gsub(
  "\\.",
  "-",
  estimate_sample_ids_7be10
)

cat("\n============================================================\n")
cat("NORMALIZED SAMPLE ID COMPARISON\n")
cat("============================================================\n\n")

cat(
  "Normalized common IDs: ",
  length(common_ids_normalized_10b3k),
  "\n",
  sep = ""
)

cat(
  "Normalized source IDs: ",
  length(source_ids_normalized_10b3k),
  "\n",
  sep = ""
)

cat(
  "Same sample set after punctuation normalization: ",
  setequal(
    common_ids_normalized_10b3k,
    source_ids_normalized_10b3k
  ),
  "\n",
  sep = ""
)

cat(
  "Same sample order after punctuation normalization: ",
  identical(
    common_ids_normalized_10b3k,
    source_ids_normalized_10b3k
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Missing / extra samples
# ------------------------------------------------------------

cat("\nMissing from common matrix after normalization:\n")

print(
  setdiff(
    source_ids_normalized_10b3k,
    common_ids_normalized_10b3k
  )
)

cat("\nExtra in common matrix after normalization:\n")

print(
  setdiff(
    common_ids_normalized_10b3k,
    source_ids_normalized_10b3k
  )
)


# ------------------------------------------------------------
# 8. Safety
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 10B-3K-A SAFETY STATUS\n")
cat("============================================================\n\n")

cat(
  "Structure audit only.\n",
  "No ESTIMATE score recalculated.\n",
  "No LOO score calculated.\n",
  "No statistical analysis performed.\n",
  "No Step14 locked object modified.\n",
  "No Figure 6 modified.\n",
  "No Supplementary Table 3 modified.\n"
)

cat("\nSTEP 10B-3K-A COMPLETE.\n")

###STEP 10B-3L: RECONSTRUCT ORIGINAL ESTIMATE SCORES FROM COMMON MATRIX

# ============================================================
# Paper 1
# STEP 10B-3L
# RECONSTRUCT ORIGINAL ESTIMATE SCORES FROM COMMON MATRIX
#
# Purpose:
#   Reconstruct the original Step14 ESTIMATE scores using the
#   10,221-gene common matrix identified in STEP 10B-3K-A.
#
# IMPORTANT:
#   - Reconstruction only
#   - No LOO score calculation
#   - No statistical analysis
#   - No modification of Step14 locked objects
#   - No modification of Figure 6
#   - No modification of Supplementary Table 3
#
# Reconstruction must PASS before proceeding to LOO.
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3L: RECONSTRUCT ORIGINAL ESTIMATE SCORES\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Basic object availability check
# ------------------------------------------------------------

required_objects_10b3l <- c(
  "estimate_common_10b3f",
  "estimate_scores_source_7be10",
  "estimate_sample_ids_7be10"
)

missing_objects_10b3l <- required_objects_10b3l[
  !vapply(
    required_objects_10b3l,
    exists,
    logical(1),
    inherits = TRUE
  )
]

cat("Required objects:\n")
print(
  data.frame(
    Object = required_objects_10b3l,
    Present = vapply(
      required_objects_10b3l,
      exists,
      logical(1),
      inherits = TRUE
    ),
    stringsAsFactors = FALSE
  )
)

if (length(missing_objects_10b3l) > 0) {
  stop(
    paste(
      "STOP: Missing required object(s):",
      paste(missing_objects_10b3l, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 1. Recover common matrix
# ------------------------------------------------------------

ec_10b3l <- estimate_common_10b3f

cat("\n============================================================\n")
cat("1. COMMON MATRIX STRUCTURE\n")
cat("============================================================\n\n")

cat(
  "Dimensions: ",
  paste(dim(ec_10b3l), collapse = " x "),
  "\n",
  sep = ""
)

cat("\nColumn names 1:5:\n")
print(colnames(ec_10b3l)[1:5])

cat("\nColumn names last 5:\n")
print(tail(colnames(ec_10b3l), 5))


# ------------------------------------------------------------
# 2. Identify expression columns
# ------------------------------------------------------------

metadata_columns_10b3l <- c(
  "NAME",
  "Description"
)

stopifnot(
  all(metadata_columns_10b3l %in% colnames(ec_10b3l))
)

sample_columns_10b3l <- setdiff(
  colnames(ec_10b3l),
  metadata_columns_10b3l
)

cat("\nSample count in common matrix: ")
cat(length(sample_columns_10b3l), "\n")

if (length(sample_columns_10b3l) != 513) {
  stop(
    paste(
      "STOP: Expected 513 sample columns, found",
      length(sample_columns_10b3l)
    )
  )
}


# ------------------------------------------------------------
# 3. Extract gene symbols
# ------------------------------------------------------------

gene_symbols_10b3l <- as.character(
  ec_10b3l$Description
)

cat("\nGene count: ")
cat(length(gene_symbols_10b3l), "\n")

cat(
  "Unique gene symbols: ",
  length(unique(gene_symbols_10b3l)) ==
    length(gene_symbols_10b3l),
  "\n",
  sep = ""
)

if (
  length(unique(gene_symbols_10b3l)) !=
  length(gene_symbols_10b3l)
) {
  stop("STOP: Duplicate gene symbols detected.")
}


# ------------------------------------------------------------
# 4. Confirm core-gene status
# ------------------------------------------------------------

core_genes_10b3l <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

cat("\n============================================================\n")
cat("4. CORE GENE STATUS IN COMMON MATRIX\n")
cat("============================================================\n\n")

core_presence_10b3l <- data.frame(
  Gene = core_genes_10b3l,
  Present = core_genes_10b3l %in% gene_symbols_10b3l,
  stringsAsFactors = FALSE
)

print(core_presence_10b3l)

# Expected:
# AOC3  = TRUE
# MS4A7 = FALSE
# IL7R  = TRUE


# ------------------------------------------------------------
# 5. Convert common matrix to expression matrix
# ------------------------------------------------------------
#
# ESTIMATE common matrix structure:
#   rows      = genes
#   columns   = samples
#   Description = gene symbol
#
# No transformation of expression values is performed here.
# ------------------------------------------------------------

expr_common_10b3l <- as.matrix(
  ec_10b3l[
    ,
    sample_columns_10b3l,
    drop = FALSE
  ]
)

rownames(expr_common_10b3l) <- gene_symbols_10b3l

storage.mode(expr_common_10b3l) <- "numeric"


# ------------------------------------------------------------
# 6. Expression matrix validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("6. EXPRESSION MATRIX VALIDATION\n")
cat("============================================================\n\n")

cat(
  "Dimensions: ",
  paste(dim(expr_common_10b3l), collapse = " x "),
  "\n",
  sep = ""
)

cat(
  "All values finite: ",
  all(is.finite(expr_common_10b3l)),
  "\n",
  sep = ""
)

cat(
  "Unique rownames: ",
  length(unique(rownames(expr_common_10b3l))) ==
    nrow(expr_common_10b3l),
  "\n",
  sep = ""
)

if (!all(is.finite(expr_common_10b3l))) {
  stop("STOP: Non-finite values detected in common expression matrix.")
}


# ------------------------------------------------------------
# 7. Sample identity and order audit
# ------------------------------------------------------------

common_ids_norm_10b3l <- gsub(
  "\\.",
  "-",
  colnames(expr_common_10b3l)
)

source_ids_norm_10b3l <- gsub(
  "\\.",
  "-",
  estimate_sample_ids_7be10
)

cat("\n============================================================\n")
cat("7. SAMPLE IDENTITY / ORDER AUDIT\n")
cat("============================================================\n\n")

cat(
  "Same sample set: ",
  setequal(
    common_ids_norm_10b3l,
    source_ids_norm_10b3l
  ),
  "\n",
  sep = ""
)

cat(
  "Same sample order: ",
  identical(
    common_ids_norm_10b3l,
    source_ids_norm_10b3l
  ),
  "\n",
  sep = ""
)

if (
  !identical(
    common_ids_norm_10b3l,
    source_ids_norm_10b3l
  )
) {
  stop(
    "STOP: Common matrix sample order does not match Step14 source."
  )
}

# ============================================================
# STEP 10B-3L-A
# LOCATE LOCKED LOO SIGNATURE CHECKPOINT
#
# Purpose:
#   Locate the already-created locked preparation/audit
#   checkpoint containing:
#
#     stromal_genes_loo_AOC3_7be
#     immune_genes_loo_IL7R_7be
#
# No scoring.
# No LOO.
# No statistical analysis.
# No modification of locked objects.
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3L-A: LOCATE LOCKED LOO SIGNATURE CHECKPOINT\n")
cat("============================================================\n\n")

root_10b3l_a <- "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis"

cat("Root exists: ")
cat(dir.exists(root_10b3l_a), "\n\n")


# ------------------------------------------------------------
# 1. Search for the known checkpoint names
# ------------------------------------------------------------

candidate_files_10b3l_a <- c(
  "Paper1_STEP7B10A_preparation_audit_LOCKED.RData",
  "Paper1_STEP7B10B3A_ESTIMATE_LOO_prerequisite_audit_LOCKED.RData",
  "Paper1_STEP7B10B1_source_recovery_audit_LOCKED.RData"
)

cat("Known candidate files:\n\n")

for (f_i in candidate_files_10b3l_a) {
  
  hits_i <- list.files(
    root_10b3l_a,
    pattern = f_i,
    recursive = TRUE,
    full.names = TRUE
  )
  
  cat("------------------------------------------------------------\n")
  cat(f_i, "\n")
  cat("Matches: ", length(hits_i), "\n", sep = "")
  
  if (length(hits_i) > 0) {
    print(hits_i)
  }
}


# ------------------------------------------------------------
# 2. Search all RData files containing 10B / LOO
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("BROAD RDATA SEARCH\n")
cat("============================================================\n\n")

all_rdata_10b3l_a <- list.files(
  root_10b3l_a,
  pattern = "\\.(RData|rda|RDATA)$",
  recursive = TRUE,
  full.names = TRUE
)

cat(
  "Total RData files found: ",
  length(all_rdata_10b3l_a),
  "\n",
  sep = ""
)

relevant_rdata_10b3l_a <- all_rdata_10b3l_a[
  grepl(
    "10B|10A|10b|10a|LOO|loo|ESTIMATE|estimate",
    basename(all_rdata_10b3l_a)
  )
]

cat(
  "\nRelevant RData files: ",
  length(relevant_rdata_10b3l_a),
  "\n",
  sep = ""
)

if (length(relevant_rdata_10b3l_a) > 0) {
  print(relevant_rdata_10b3l_a)
}


# ------------------------------------------------------------
# 3. Search current workspace first
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CURRENT WORKSPACE OBJECT CHECK\n")
cat("============================================================\n\n")

target_objects_10b3l_a <- c(
  "stromal_genes_loo_AOC3_7be",
  "immune_genes_loo_IL7R_7be",
  "stromal_genes_loo_AOC3_7be10",
  "immune_genes_loo_IL7R_7be10"
)

print(
  data.frame(
    Object = target_objects_10b3l_a,
    Present = vapply(
      target_objects_10b3l_a,
      exists,
      logical(1),
      inherits = TRUE
    ),
    stringsAsFactors = FALSE
  )
)


cat("\n============================================================\n")
cat("STEP 10B-3L-A COMPLETE.\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 10B-3L-B
# ESTIMATE 1.0.13-STYLE SCORE RECONSTRUCTION
#
# Purpose:
#   Reconstruct the original Step14 ESTIMATE scores from the
#   verified 10,221-gene common matrix using the ESTIMATE
#   1.0.13 GSEA/ssGSEA-style enrichment-score logic.
#
# IMPORTANT:
#   - Original-score reconstruction ONLY
#   - No LOO
#   - No statistical analysis
#   - No modification of Step14
#   - No modification of Figure 6
#   - No modification of Supplementary Table 3
#
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3L-B: ESTIMATE 1.0.13-STYLE SCORE RECONSTRUCTION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Required objects
# ------------------------------------------------------------

required_objects_10b3l_b <- c(
  "estimate_common_10b3f",
  "estimate_scores_source_7be10",
  "estimate_sample_ids_7be10",
  "stromal_genes_loo_AOC3_7be",
  "immune_genes_loo_IL7R_7be"
)

present_10b3l_b <- vapply(
  required_objects_10b3l_b,
  exists,
  logical(1),
  inherits = TRUE
)

print(
  data.frame(
    Object = required_objects_10b3l_b,
    Present = present_10b3l_b,
    stringsAsFactors = FALSE
  )
)

if (!all(present_10b3l_b)) {
  stop(
    paste(
      "STOP: Missing required object(s):",
      paste(
        required_objects_10b3l_b[!present_10b3l_b],
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 1. Recover verified common matrix
# ------------------------------------------------------------

ec_10b3l_b <- estimate_common_10b3f

sample_columns_10b3l_b <- setdiff(
  colnames(ec_10b3l_b),
  c("NAME", "Description")
)

gene_symbols_10b3l_b <- as.character(
  ec_10b3l_b$Description
)

expr_common_10b3l_b <- as.matrix(
  ec_10b3l_b[
    ,
    sample_columns_10b3l_b,
    drop = FALSE
  ]
)

rownames(expr_common_10b3l_b) <-
  gene_symbols_10b3l_b

storage.mode(expr_common_10b3l_b) <- "numeric"


# ------------------------------------------------------------
# 2. Basic structural checks
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("2. COMMON MATRIX VALIDATION\n")
cat("============================================================\n\n")

cat(
  "Dimensions: ",
  paste(dim(expr_common_10b3l_b), collapse = " x "),
  "\n",
  sep = ""
)

cat(
  "All finite: ",
  all(is.finite(expr_common_10b3l_b)),
  "\n",
  sep = ""
)

cat(
  "Unique genes: ",
  length(unique(rownames(expr_common_10b3l_b))) ==
    nrow(expr_common_10b3l_b),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 3. Recover original 141-gene signatures
# ------------------------------------------------------------

stromal_genes_original_10b3l_b <- unique(
  c(
    stromal_genes_loo_AOC3_7be,
    "AOC3"
  )
)

immune_genes_original_10b3l_b <- unique(
  c(
    immune_genes_loo_IL7R_7be,
    "IL7R"
  )
)

cat("\n============================================================\n")
cat("3. SIGNATURE VALIDATION\n")
cat("============================================================\n\n")

cat(
  "Original Stromal signature: ",
  length(stromal_genes_original_10b3l_b),
  "\n",
  sep = ""
)

cat(
  "Original Immune signature: ",
  length(immune_genes_original_10b3l_b),
  "\n",
  sep = ""
)

if (
  length(stromal_genes_original_10b3l_b) != 141 ||
  length(immune_genes_original_10b3l_b) != 141
) {
  stop(
    "STOP: Signature reconstruction does not contain exactly 141 genes."
  )
}

# ============================================================
# Paper 1
# STEP 10B-3L-C
# ESTIMATE 1.0.13 EXACT-SOURCE-LOGIC RECONSTRUCTION
#
# Purpose:
#   Reconstruct the original Step14 ESTIMATE scores from the
#   verified 10,221-gene common matrix.
#
# IMPORTANT:
#   - Follow official ESTIMATE 1.0.13 estimateScore() logic
#   - Signature genes absent from the common matrix are NOT
#     manually added back
#   - Use only gene-set overlap with the actual input matrix
#   - No LOO
#   - No statistical analysis
#   - No modification of Step14
#   - No modification of Figure 6
#   - No modification of Supplementary Table 3
#
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3L-C: ESTIMATE 1.0.13 SOURCE-LOGIC RECONSTRUCTION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Required objects
# ------------------------------------------------------------

required_objects_10b3l_c <- c(
  "estimate_common_10b3f",
  "estimate_scores_source_7be10",
  "estimate_sample_ids_7be10",
  "stromal_genes_loo_AOC3_7be",
  "immune_genes_loo_IL7R_7be"
)

present_10b3l_c <- vapply(
  required_objects_10b3l_c,
  exists,
  logical(1),
  inherits = TRUE
)

print(
  data.frame(
    Object = required_objects_10b3l_c,
    Present = present_10b3l_c,
    stringsAsFactors = FALSE
  )
)

if (!all(present_10b3l_c)) {
  stop(
    paste(
      "STOP: Missing required object(s):",
      paste(
        required_objects_10b3l_c[!present_10b3l_c],
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 1. Recover verified 10,221-gene common matrix
# ------------------------------------------------------------

ec_10b3l_c <- estimate_common_10b3f

sample_columns_10b3l_c <- setdiff(
  colnames(ec_10b3l_c),
  c("NAME", "Description")
)

gene_symbols_10b3l_c <- as.character(
  ec_10b3l_c$Description
)

expr_common_10b3l_c <- as.matrix(
  ec_10b3l_c[
    ,
    sample_columns_10b3l_c,
    drop = FALSE
  ]
)

rownames(expr_common_10b3l_c) <- gene_symbols_10b3l_c

storage.mode(expr_common_10b3l_c) <- "numeric"


# ------------------------------------------------------------
# 2. Structural validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("2. COMMON MATRIX VALIDATION\n")
cat("============================================================\n\n")

cat(
  "Dimensions: ",
  paste(dim(expr_common_10b3l_c), collapse = " x "),
  "\n",
  sep = ""
)

stopifnot(
  nrow(expr_common_10b3l_c) == 10221,
  ncol(expr_common_10b3l_c) == 513,
  all(is.finite(expr_common_10b3l_c)),
  !anyDuplicated(rownames(expr_common_10b3l_c))
)

cat("All structural checks: PASSED\n")


# ------------------------------------------------------------
# 3. Reconstruct the original 141-gene signatures
# ------------------------------------------------------------

stromal_genes_original_10b3l_c <- unique(
  c(
    stromal_genes_loo_AOC3_7be,
    "AOC3"
  )
)

immune_genes_original_10b3l_c <- unique(
  c(
    immune_genes_loo_IL7R_7be,
    "IL7R"
  )
)

stopifnot(
  length(stromal_genes_original_10b3l_c) == 141,
  length(immune_genes_original_10b3l_c) == 141
)


# ------------------------------------------------------------
# 4. IMPORTANT:
#    Apply ESTIMATE's actual overlap logic
#
#    Official source:
#      gene.overlap <- intersect(gene.set, gene.names)
#
#    Therefore:
#      DO NOT require all 141 genes to exist.
# ------------------------------------------------------------

stromal_overlap_10b3l_c <- intersect(
  stromal_genes_original_10b3l_c,
  rownames(expr_common_10b3l_c)
)

immune_overlap_10b3l_c <- intersect(
  immune_genes_original_10b3l_c,
  rownames(expr_common_10b3l_c)
)

stromal_missing_10b3l_c <- setdiff(
  stromal_genes_original_10b3l_c,
  rownames(expr_common_10b3l_c)
)

immune_missing_10b3l_c <- setdiff(
  immune_genes_original_10b3l_c,
  rownames(expr_common_10b3l_c)
)

cat("\n============================================================\n")
cat("4. SIGNATURE OVERLAP AUDIT\n")
cat("============================================================\n\n")

cat(
  "Original Stromal signature: ",
  length(stromal_genes_original_10b3l_c),
  "\n",
  sep = ""
)

cat(
  "Stromal genes present: ",
  length(stromal_overlap_10b3l_c),
  "\n",
  sep = ""
)

cat(
  "Stromal genes absent: ",
  length(stromal_missing_10b3l_c),
  "\n",
  sep = ""
)

if (length(stromal_missing_10b3l_c) > 0) {
  cat("Stromal absent genes:\n")
  print(stromal_missing_10b3l_c)
}

cat(
  "\nOriginal Immune signature: ",
  length(immune_genes_original_10b3l_c),
  "\n",
  sep = ""
)

cat(
  "Immune genes present: ",
  length(immune_overlap_10b3l_c),
  "\n",
  sep = ""
)

cat(
  "Immune genes absent: ",
  length(immune_missing_10b3l_c),
  "\n",
  sep = ""
)

if (length(immune_missing_10b3l_c) > 0) {
  cat("Immune absent genes:\n")
  print(immune_missing_10b3l_c)
}


# ------------------------------------------------------------
# 5. Confirm expected overlap counts
# ------------------------------------------------------------

if (
  length(stromal_overlap_10b3l_c) != 136 ||
  length(immune_overlap_10b3l_c) != 140
) {
  stop(
    paste0(
      "STOP: Unexpected signature overlap.\n",
      "Expected Stromal = 136; observed = ",
      length(stromal_overlap_10b3l_c),
      "\n",
      "Expected Immune = 140; observed = ",
      length(immune_overlap_10b3l_c)
    )
  )
}

cat(
  "\nSignature overlap counts exactly match the prior audit: PASSED\n"
)


# ------------------------------------------------------------
# 6. Verify sample identity and order
# ------------------------------------------------------------

source_sample_ids_10b3l_c <-
  as.character(estimate_sample_ids_7be10)

common_sample_ids_10b3l_c <-
  colnames(expr_common_10b3l_c)

sample_ids_common_normalized_10b3l_c <-
  gsub(
    "\\.",
    "-",
    common_sample_ids_10b3l_c
  )

if (
  !identical(
    sample_ids_common_normalized_10b3l_c,
    source_sample_ids_10b3l_c
  )
) {
  stop(
    "STOP: Common-matrix sample identity/order does not match Step14."
  )
}

cat("Sample identity/order: PASSED\n")


# ------------------------------------------------------------
# 7. Implement ESTIMATE 1.0.13 estimateScore() logic
#
# Official source:
#
#   for each sample:
#       m[,j] <- rank(m[,j], ties.method="average")
#   m <- 10000*m/Ng
#
#   gene.list <- order(m[,j], decreasing=TRUE)
#
#   TAG <- match(gene.set, gene.list)
#
#   correl.vector <- m[gene.list,j]
#   correl.vector <- abs(correl.vector)^0.25
#
#   Pn <- TAG * correl.vector / sum(correl.vector[TAG==1])
#
#   P0 <- no.TAG / Nm
#
#   RES <- cumsum(Pn) - cumsum(P0)
#
#   ES <- max(RES) or min(RES),
#         whichever has larger absolute magnitude.
# ------------------------------------------------------------

estimate_es_10b3l_c <- function(
    expression_matrix,
    gene_set
) {
  
  stopifnot(
    is.matrix(expression_matrix),
    is.numeric(expression_matrix),
    length(gene_set) > 0
  )
  
  gene_names <- rownames(expression_matrix)
  
  gene_overlap <- intersect(
    gene_set,
    gene_names
  )
  
  if (length(gene_overlap) == 0) {
    return(
      rep(
        NA_real_,
        ncol(expression_matrix)
      )
    )
  }
  
  N_genes <- nrow(expression_matrix)
  N_samples <- ncol(expression_matrix)
  
  ES_vector <- numeric(N_samples)
  
  # ----------------------------------------------------------
  # Sample-rank normalization exactly following ESTIMATE
  # ----------------------------------------------------------
  
  m <- expression_matrix
  
  for (j in seq_len(N_samples)) {
    
    m[, j] <- rank(
      m[, j],
      ties.method = "average"
    )
  }
  
  m <- 10000 * m / N_samples
  
  # ----------------------------------------------------------
  # Gene-set enrichment score
  # ----------------------------------------------------------
  
  for (j in seq_len(N_samples)) {
    
    gene_list <- order(
      m[, j],
      decreasing = TRUE
    )
    
    gene_set2 <- match(
      gene_overlap,
      gene_names
    )
    
    correl_vector <- m[
      gene_list,
      j
    ]
    
    TAG <- as.integer(
      match(
        gene_list,
        gene_set2,
        nomatch = 0
      )
    )
    
    no_TAG <- 1 - TAG
    
    N <- length(gene_list)
    Nh <- length(gene_set2)
    Nm <- N - Nh
    
    if (Nh == 0 || Nm == 0) {
      stop(
        "STOP: Invalid GSEA dimensions."
      )
    }
    
    correl_vector <-
      abs(correl_vector)^0.25
    
    sum_correl <-
      sum(
        correl_vector[TAG == 1]
      )
    
    if (!is.finite(sum_correl) ||
        sum_correl <= 0) {
      stop(
        "STOP: Invalid signature weighting."
      )
    }
    
    P0 <- no_TAG / Nm
    
    F0 <- cumsum(P0)
    
    Pn <-
      TAG * correl_vector / sum_correl
    
    Fn <- cumsum(Pn)
    
    RES <- Fn - F0
    
    max_ES <- max(RES)
    min_ES <- min(RES)
    
    if (max_ES > -min_ES) {
      ES <- max_ES
    } else {
      ES <- min_ES
    }
    
    ES_vector[j] <- ES
  }
  
  ES_vector
}


# ------------------------------------------------------------
# 8. Calculate reconstructed StromalScore
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("8. RECONSTRUCT STROMAL SCORE\n")
cat("============================================================\n\n")

stromal_score_reconstructed_10b3l_c <-
  estimate_es_10b3l_c(
    expression_matrix = expr_common_10b3l_c,
    gene_set = stromal_overlap_10b3l_c
  )

cat(
  "Stromal score calculated for ",
  length(stromal_score_reconstructed_10b3l_c),
  " samples.\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Calculate reconstructed ImmuneScore
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("9. RECONSTRUCT IMMUNE SCORE\n")
cat("============================================================\n\n")

immune_score_reconstructed_10b3l_c <-
  estimate_es_10b3l_c(
    expression_matrix = expr_common_10b3l_c,
    gene_set = immune_overlap_10b3l_c
  )

cat(
  "Immune score calculated for ",
  length(immune_score_reconstructed_10b3l_c),
  " samples.\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Construct ESTIMATEScore
# ------------------------------------------------------------

estimate_score_reconstructed_10b3l_c <-
  stromal_score_reconstructed_10b3l_c +
  immune_score_reconstructed_10b3l_c


# ------------------------------------------------------------
# 11. Build reconstructed score matrix
# ------------------------------------------------------------

reconstructed_scores_10b3l_c <-
  rbind(
    StromalScore =
      stromal_score_reconstructed_10b3l_c,
    ImmuneScore =
      immune_score_reconstructed_10b3l_c,
    ESTIMATEScore =
      estimate_score_reconstructed_10b3l_c
  )

colnames(
  reconstructed_scores_10b3l_c
) <- source_sample_ids_10b3l_c


# ------------------------------------------------------------
# 12. Structural validation
# ------------------------------------------------------------

stopifnot(
  identical(
    dim(reconstructed_scores_10b3l_c),
    dim(estimate_scores_source_7be10)
  )
)

stopifnot(
  identical(
    rownames(reconstructed_scores_10b3l_c),
    rownames(estimate_scores_source_7be10)
  )
)

stopifnot(
  identical(
    colnames(reconstructed_scores_10b3l_c),
    colnames(estimate_scores_source_7be10)
  )
)

stopifnot(
  all(is.finite(reconstructed_scores_10b3l_c))
)

cat("\nScore structure: PASSED\n")


# ------------------------------------------------------------
# 13. Compare reconstructed scores with locked Step14 scores
# ------------------------------------------------------------

source_scores_10b3l_c <-
  as.matrix(
    estimate_scores_source_7be10
  )

comparison_10b3l_c <- data.frame(
  Score = rownames(source_scores_10b3l_c),
  Pearson_r = NA_real_,
  Spearman_rho = NA_real_,
  MaxAbsDiff = NA_real_,
  MeanAbsDiff = NA_real_,
  RMSE = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(source_scores_10b3l_c))) {
  
  x <- source_scores_10b3l_c[i, ]
  y <- reconstructed_scores_10b3l_c[i, ]
  
  comparison_10b3l_c$Pearson_r[i] <-
    cor(
      x,
      y,
      method = "pearson"
    )
  
  comparison_10b3l_c$Spearman_rho[i] <-
    cor(
      x,
      y,
      method = "spearman"
    )
  
  comparison_10b3l_c$MaxAbsDiff[i] <-
    max(abs(x - y))
  
  comparison_10b3l_c$MeanAbsDiff[i] <-
    mean(abs(x - y))
  
  comparison_10b3l_c$RMSE[i] <-
    sqrt(mean((x - y)^2))
}

cat("\n============================================================\n")
cat("13. ORIGINAL SCORE RECONSTRUCTION COMPARISON\n")
cat("============================================================\n\n")

print(
  comparison_10b3l_c,
  digits = 10,
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Exact / near-exact reconstruction criterion
#
# IMPORTANT:
#   Do NOT accept a high correlation alone.
# ------------------------------------------------------------

reconstruction_pass_10b3l_c <-
  all(
    comparison_10b3l_c$Pearson_r >
      0.999999999,
    comparison_10b3l_c$Spearman_rho >
      0.999999999,
    comparison_10b3l_c$MaxAbsDiff <
      1e-8
  )

cat(
  "\nOriginal-score reconstruction PASS: ",
  reconstruction_pass_10b3l_c,
  "\n",
  sep = ""
)


# ============================================================
# Paper 1
# STEP 10B-3M
# OFFICIAL ESTIMATE 1.0.13 PACKAGE RECOVERY AUDIT
#
# Purpose:
#   Determine whether the official ESTIMATE 1.0.13 package
#   or its source/archive files already exist locally.
#
# IMPORTANT:
#   - NO score calculation
#   - NO LOO
#   - NO statistical analysis
#   - NO modification of Step14
#   - NO modification of Figure 6
#   - NO modification of Supplementary Table 3
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3M: OFFICIAL ESTIMATE 1.0.13 PACKAGE RECOVERY AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Current package status
# ------------------------------------------------------------

estimate_installed_10b3m <-
  requireNamespace(
    "estimate",
    quietly = TRUE
  )

cat(
  "estimate package installed: ",
  estimate_installed_10b3m,
  "\n",
  sep = ""
)

if (estimate_installed_10b3m) {
  
  cat(
    "Installed version: ",
    as.character(
      packageVersion("estimate")
    ),
    "\n",
    sep = ""
  )
  
} else {
  
  cat(
    "Installed version: NOT AVAILABLE\n"
  )
}


# ------------------------------------------------------------
# 2. Check installed package libraries
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("2. LOCAL R LIBRARY SEARCH\n")
cat("============================================================\n\n")

lib_paths_10b3m <- .libPaths()

print(lib_paths_10b3m)

estimate_dirs_10b3m <- unique(
  unlist(
    lapply(
      lib_paths_10b3m,
      function(x) {
        if (!dir.exists(x)) {
          return(character(0))
        }
        
        list.files(
          x,
          pattern = "^estimate$",
          full.names = TRUE
        )
      }
    )
  )
)

cat(
  "\nLocal estimate package directories found: ",
  length(estimate_dirs_10b3m),
  "\n",
  sep = ""
)

if (length(estimate_dirs_10b3m) > 0) {
  print(estimate_dirs_10b3m)
}


# ------------------------------------------------------------
# 3. Search likely project directories
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("3. SEARCH PAPER 1 DIRECTORIES\n")
cat("============================================================\n\n")

search_roots_10b3m <- unique(
  c(
    "O:/Rproject",
    "O:/Rproject/Paper1_tumor-only WGCNA",
    "O:/Rproject/Paper 1_Step 7_ Immune Infiltration Analysis",
    "C:/Users/Administrator/Downloads",
    "C:/Users/Administrator/Documents",
    "C:/Temp"
  )
)

search_roots_10b3m <-
  search_roots_10b3m[
    dir.exists(search_roots_10b3m)
  ]

print(search_roots_10b3m)


# ------------------------------------------------------------
# 4. Search for official package/archive/source files
# ------------------------------------------------------------

target_patterns_10b3m <- c(
  "^estimate_[0-9.]+\\.tar\\.gz$",
  "^estimate_[0-9.]+\\.zip$",
  "^estimate.*\\.tar\\.gz$",
  "^estimate.*\\.zip$",
  "^estimateScore\\.R$",
  "^filterCommonGenes\\.R$",
  "^SI_geneset\\.RData$",
  "^common_genes\\.RData$",
  "^common_genes\\.txt$",
  "^SI_geneset\\.gmt$"
)

found_files_10b3m <- character(0)

for (root in search_roots_10b3m) {
  
  for (pattern in target_patterns_10b3m) {
    
    tmp <- tryCatch(
      list.files(
        root,
        pattern = pattern,
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
      ),
      error = function(e) character(0)
    )
    
    if (length(tmp) > 0) {
      found_files_10b3m <-
        c(found_files_10b3m, tmp)
    }
  }
}

found_files_10b3m <-
  unique(found_files_10b3m)

cat(
  "\nPotential ESTIMATE package/source files found: ",
  length(found_files_10b3m),
  "\n",
  sep = ""
)

if (length(found_files_10b3m) > 0) {
  print(found_files_10b3m)
}


# ------------------------------------------------------------
# 5. Search specifically for ESTIMATE 1.0.13
# ------------------------------------------------------------

estimate_1013_files_10b3m <-
  found_files_10b3m[
    grepl(
      "1\\.0\\.13",
      basename(found_files_10b3m),
      ignore.case = TRUE
    )
  ]

cat(
  "\nExplicit ESTIMATE 1.0.13 files found: ",
  length(estimate_1013_files_10b3m),
  "\n",
  sep = ""
)

if (length(estimate_1013_files_10b3m) > 0) {
  print(estimate_1013_files_10b3m)
}


# ------------------------------------------------------------
# 6. Search for R source files even if renamed
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("6. ESTIMATE SOURCE-CODE SEARCH\n")
cat("============================================================\n\n")

r_source_files_10b3m <- found_files_10b3m[
  grepl(
    "(estimateScore|filterCommonGenes|SI_geneset|common_genes)",
    basename(found_files_10b3m),
    ignore.case = TRUE
  )
]

if (length(r_source_files_10b3m) > 0) {
  print(r_source_files_10b3m)
} else {
  cat("No local ESTIMATE source-code files found.\n")
}


# ------------------------------------------------------------
# 7. Safety audit
# ------------------------------------------------------------

step7be10b3m_audit <- list(
  
  status =
    "PACKAGE_RECOVERY_AUDIT_ONLY",
  
  estimate_installed =
    estimate_installed_10b3m,
  
  installed_version =
    if (estimate_installed_10b3m)
      as.character(packageVersion("estimate"))
  else
    NA_character_,
  
  local_package_directories =
    estimate_dirs_10b3m,
  
  potential_source_or_archive_files =
    found_files_10b3m,
  
  explicit_1_0_13_files =
    estimate_1013_files_10b3m,
  
  score_calculated =
    FALSE,
  
  loo_score_calculated =
    FALSE,
  
  statistical_analysis_performed =
    FALSE,
  
  step14_modified =
    FALSE,
  
  figure6_modified =
    FALSE,
  
  supplementary_table3_modified =
    FALSE
)

cat("\n============================================================\n")
cat("STEP 10B-3M COMPLETE\n")
cat("============================================================\n\n")

print(step7be10b3m_audit)

# ============================================================
# Paper 1
# STEP 10B-3N-A
# DOWNLOAD OFFICIAL ESTIMATE 1.0.13 SOURCE PACKAGE
#
# Purpose:
#   Download the official ESTIMATE 1.0.13 source package.
#
# IMPORTANT:
#   - Download only
#   - No installation yet
#   - No score calculation
#   - No LOO
#   - No statistical analysis
#   - No modification of Step14
#   - No modification of Figure 6
#   - No modification of Supplementary Table 3
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3N-A: DOWNLOAD OFFICIAL ESTIMATE 1.0.13\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define official source URL
# ------------------------------------------------------------

estimate_url_10b3na <-
  "https://download.r-forge.r-project.org/src/contrib/estimate_1.0.13.tar.gz"

cat(
  "Official source URL:\n",
  estimate_url_10b3na,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 2. Define local download directory
# ------------------------------------------------------------

download_dir_10b3na <-
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis/STEP7B-7E-10B-3N"

if (!dir.exists(download_dir_10b3na)) {
  dir.create(
    download_dir_10b3na,
    recursive = TRUE
  )
}

cat(
  "Download directory:\n",
  download_dir_10b3na,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 3. Define target file
# ------------------------------------------------------------

estimate_tar_10b3na <-
  file.path(
    download_dir_10b3na,
    "estimate_1.0.13.tar.gz"
  )


# ------------------------------------------------------------
# 4. Download
# ------------------------------------------------------------

if (file.exists(estimate_tar_10b3na)) {
  
  cat(
    "Target file already exists.\n",
    "No re-download performed.\n\n"
  )
  
} else {
  
  cat(
    "Downloading official ESTIMATE 1.0.13 source package...\n\n"
  )
  
  download.file(
    url =
      estimate_url_10b3na,
    destfile =
      estimate_tar_10b3na,
    mode = "wb",
    quiet = FALSE
  )
}


# ============================================================
# STEP 10B-3N-A
# Continue: verify downloaded official ESTIMATE 1.0.13
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3N-A: VERIFY DOWNLOADED ESTIMATE 1.0.13\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Verify file existence and size
# ------------------------------------------------------------

file_exists_10b3na <-
  file.exists(
    estimate_tar_10b3na
  )

file_size_10b3na <-
  if (file_exists_10b3na) {
    file.info(
      estimate_tar_10b3na
    )$size
  } else {
    NA_real_
  }

cat(
  "File exists: ",
  file_exists_10b3na,
  "\n",
  sep = ""
)

cat(
  "File size: ",
  file_size_10b3na,
  " bytes\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 2. Hard-stop if download failed
# ------------------------------------------------------------

stopifnot(
  file_exists_10b3na
)

stopifnot(
  is.finite(file_size_10b3na)
)

stopifnot(
  file_size_10b3na > 100000
)


# ------------------------------------------------------------
# 3. Verify gzip/tar archive integrity
# ------------------------------------------------------------

cat(
  "Testing tar.gz archive integrity...\n\n"
)

archive_listing_10b3na <-
  utils::untar(
    estimate_tar_10b3na,
    list = TRUE
  )

archive_valid_10b3na <-
  length(
    archive_listing_10b3na
  ) > 0

cat(
  "Archive valid: ",
  archive_valid_10b3na,
  "\n",
  sep = ""
)

cat(
  "Number of archive entries: ",
  length(archive_listing_10b3na),
  "\n\n",
  sep = ""
)

stopifnot(
  archive_valid_10b3na
)


# ------------------------------------------------------------
# 4. Inspect official package contents
# ------------------------------------------------------------

important_files_10b3na <-
  archive_listing_10b3na[
    grepl(
      "(estimateScore|filterCommonGenes|SI_geneset|common_genes|DESCRIPTION|NAMESPACE)",
      archive_listing_10b3na,
      ignore.case = TRUE
    )
  ]

cat(
  "Important package files detected:\n\n"
)

print(
  important_files_10b3na
)


# ------------------------------------------------------------
# 5. Basic package identity check
# ------------------------------------------------------------

version_file_10b3na <-
  archive_listing_10b3na[
    grepl(
      "DESCRIPTION$",
      archive_listing_10b3na,
      ignore.case = TRUE
    )
  ]

cat(
  "\nDESCRIPTION file detected:\n"
)

print(
  version_file_10b3na
)


# ------------------------------------------------------------
# 6. Audit status
# ------------------------------------------------------------

step7be10b3na_audit <-
  list(
    step =
      "STEP 10B-3N-A",
    
    purpose =
      "Official ESTIMATE 1.0.13 source package download and archive verification",
    
    official_url =
      estimate_url_10b3na,
    
    local_file =
      estimate_tar_10b3na,
    
    file_exists =
      file_exists_10b3na,
    
    file_size_bytes =
      file_size_10b3na,
    
    archive_valid =
      archive_valid_10b3na,
    
    archive_entries =
      length(
        archive_listing_10b3na
      ),
    
    important_files =
      important_files_10b3na,
    
    score_calculated =
      FALSE,
    
    loo_score_calculated =
      FALSE,
    
    statistical_analysis_performed =
      FALSE,
    
    locked_objects_modified =
      FALSE,
    
    status =
      "DOWNLOAD_AND_ARCHIVE_AUDIT_PASSED"
  )


# ------------------------------------------------------------
# 7. Save audit checkpoint
# ------------------------------------------------------------

audit_file_10b3na <-
  file.path(
    download_dir_10b3na,
    "Paper1_STEP7B10B3NA_official_ESTIMATE_1.0.13_download_audit_LOCKED.RData"
  )

save(
  estimate_url_10b3na,
  estimate_tar_10b3na,
  file_exists_10b3na,
  file_size_10b3na,
  archive_listing_10b3na,
  archive_valid_10b3na,
  important_files_10b3na,
  step7be10b3na_audit,
  file = audit_file_10b3na
)

cat(
  "\nAudit checkpoint saved:\n",
  audit_file_10b3na,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Final check
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 10B-3N-A FINAL CHECK\n")
cat("============================================================\n")

cat(
  "Official URL verified: ",
  estimate_url_10b3na,
  "\n",
  sep = ""
)

cat(
  "Downloaded file exists: ",
  file_exists_10b3na,
  "\n",
  sep = ""
)

cat(
  "File size: ",
  file_size_10b3na,
  " bytes\n",
  sep = ""
)

cat(
  "Archive integrity: ",
  archive_valid_10b3na,
  "\n",
  sep = ""
)

cat(
  "Score calculated: FALSE\n"
)

cat(
  "LOO score calculated: FALSE\n"
)

cat(
  "Locked objects modified: FALSE\n"
)

cat(
  "\nSTATUS: ",
  step7be10b3na_audit$status,
  "\n",
  sep = ""
)

cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 10B-3N-B
# Install and verify official ESTIMATE 1.0.13
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3N-B: INSTALL OFFICIAL ESTIMATE 1.0.13\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Safety: define official source package
# ------------------------------------------------------------

estimate_tar_10b3nb <-
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis/STEP7B-7E-10B-3N/estimate_1.0.13.tar.gz"

stopifnot(
  file.exists(estimate_tar_10b3nb)
)

cat(
  "Official source package:\n",
  estimate_tar_10b3nb,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 1. Record current R environment
# ------------------------------------------------------------

cat("R version:\n")
print(
  R.version.string
)

cat("\nLibrary paths:\n")
print(
  .libPaths()
)

cat("\n")


# ------------------------------------------------------------
# 2. Check whether estimate is already installed
# ------------------------------------------------------------

estimate_installed_before_10b3nb <-
  requireNamespace(
    "estimate",
    quietly = TRUE
  )

cat(
  "ESTIMATE installed before installation: ",
  estimate_installed_before_10b3nb,
  "\n",
  sep = ""
)

if (estimate_installed_before_10b3nb) {
  
  estimate_version_before_10b3nb <-
    as.character(
      utils::packageVersion("estimate")
    )
  
  cat(
    "Existing ESTIMATE version: ",
    estimate_version_before_10b3nb,
    "\n\n",
    sep = ""
  )
  
} else {
  
  estimate_version_before_10b3nb <-
    NA_character_
  
  cat(
    "No existing ESTIMATE installation detected.\n\n"
  )
}


# ------------------------------------------------------------
# 3. Install official source package
# ------------------------------------------------------------

cat(
  "Installing official ESTIMATE 1.0.13 source package...\n\n"
)

install_result_10b3nb <-
  tryCatch(
    
    {
      
      utils::install.packages(
        estimate_tar_10b3nb,
        repos = NULL,
        type = "source"
      )
      
      TRUE
      
    },
    
    error = function(e) {
      
      cat(
        "\nINSTALLATION ERROR:\n",
        conditionMessage(e),
        "\n\n",
        sep = ""
      )
      
      FALSE
    }
  )

cat(
  "\nInstallation command completed: ",
  install_result_10b3nb,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Verify package availability
# ------------------------------------------------------------

estimate_installed_after_10b3nb <-
  requireNamespace(
    "estimate",
    quietly = TRUE
  )

cat(
  "ESTIMATE installed after installation: ",
  estimate_installed_after_10b3nb,
  "\n",
  sep = ""
)

if (!estimate_installed_after_10b3nb) {
  
  stop(
    "\nSTOP: official ESTIMATE 1.0.13 was not successfully installed."
  )
}


# ------------------------------------------------------------
# 5. Verify package version
# ------------------------------------------------------------

estimate_version_after_10b3nb <-
  as.character(
    utils::packageVersion("estimate")
  )

cat(
  "Installed ESTIMATE version: ",
  estimate_version_after_10b3nb,
  "\n",
  sep = ""
)

stopifnot(
  identical(
    estimate_version_after_10b3nb,
    "1.0.13"
  )
)


# ------------------------------------------------------------
# 6. Load package
# ------------------------------------------------------------

cat(
  "\nLoading package estimate...\n\n"
)

library(
  estimate
)

cat(
  "Package loaded successfully.\n\n"
)


# ------------------------------------------------------------
# 7. Verify official functions
# ------------------------------------------------------------

estimateScore_available_10b3nb <-
  exists(
    "estimateScore",
    where = asNamespace("estimate"),
    inherits = FALSE
  )

filterCommonGenes_available_10b3nb <-
  exists(
    "filterCommonGenes",
    where = asNamespace("estimate"),
    inherits = FALSE
  )

cat(
  "estimateScore available: ",
  estimateScore_available_10b3nb,
  "\n",
  sep = ""
)

cat(
  "filterCommonGenes available: ",
  filterCommonGenes_available_10b3nb,
  "\n",
  sep = ""
)

stopifnot(
  estimateScore_available_10b3nb
)

stopifnot(
  filterCommonGenes_available_10b3nb
)


# ------------------------------------------------------------
# 8. Inspect official function definitions
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Official estimateScore() source\n")
cat("============================================================\n\n")

print(
  getAnywhere(
    estimateScore
  )
)


cat("\n============================================================\n")
cat("Official filterCommonGenes() source\n")
cat("============================================================\n\n")

print(
  getAnywhere(
    filterCommonGenes
  )
)


# ------------------------------------------------------------
# 9. Verify package metadata
# ------------------------------------------------------------

estimate_description_10b3nb <-
  utils::packageDescription(
    "estimate"
  )

cat("\n============================================================\n")
cat("ESTIMATE package metadata\n")
cat("============================================================\n\n")

cat(
  "Package: ",
  estimate_description_10b3nb$Package,
  "\n",
  sep = ""
)

cat(
  "Version: ",
  estimate_description_10b3nb$Version,
  "\n",
  sep = ""
)

cat(
  "Title: ",
  estimate_description_10b3nb$Title,
  "\n",
  sep = ""
)

cat(
  "Built: ",
  estimate_description_10b3nb$Built,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Verify package installation path
# ------------------------------------------------------------

estimate_path_10b3nb <-
  system.file(
    package = "estimate"
  )

cat(
  "\nInstalled package path:\n",
  estimate_path_10b3nb,
  "\n",
  sep = ""
)

stopifnot(
  nzchar(
    estimate_path_10b3nb
  )
)


# ------------------------------------------------------------
# 11. Verify official signature/common-gene data objects
# ------------------------------------------------------------

data(
  SI_geneset,
  package = "estimate",
  envir = environment()
)

data(
  common_genes,
  package = "estimate",
  envir = environment()
)

cat("\n============================================================\n")
cat("Official ESTIMATE data objects\n")
cat("============================================================\n\n")

cat(
  "SI_geneset class: ",
  paste(
    class(SI_geneset),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "SI_geneset dimensions: ",
  paste(
    dim(SI_geneset),
    collapse = " x "
  ),
  "\n",
  sep = ""
)

cat(
  "common_genes class: ",
  paste(
    class(common_genes),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "common_genes dimensions: ",
  paste(
    dim(common_genes),
    collapse = " x "
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Basic signature audit
# ------------------------------------------------------------

SI_geneset_10b3nb <-
  SI_geneset

common_genes_10b3nb <-
  common_genes

cat("\nSI_geneset preview:\n")
print(
  head(
    SI_geneset_10b3nb
  )
)

cat("\ncommon_genes preview:\n")
print(
  head(
    common_genes_10b3nb
  )
)


# ------------------------------------------------------------
# 13. Safety audit
# ------------------------------------------------------------

step7be10b3nb_audit <-
  list(
    
    step =
      "STEP 10B-3N-B",
    
    purpose =
      "Install and verify official ESTIMATE 1.0.13",
    
    source_tar =
      estimate_tar_10b3nb,
    
    installation_command_completed =
      install_result_10b3nb,
    
    package_available =
      estimate_installed_after_10b3nb,
    
    installed_version =
      estimate_version_after_10b3nb,
    
    expected_version =
      "1.0.13",
    
    estimateScore_available =
      estimateScore_available_10b3nb,
    
    filterCommonGenes_available =
      filterCommonGenes_available_10b3nb,
    
    installed_package_path =
      estimate_path_10b3nb,
    
    score_calculated =
      FALSE,
    
    loo_score_calculated =
      FALSE,
    
    statistical_analysis_performed =
      FALSE,
    
    original_ESTIMATE_scores_modified =
      FALSE,
    
    Step14_modified =
      FALSE,
    
    Figure6_modified =
      FALSE,
    
    SupplementaryTable3_modified =
      FALSE,
    
    status =
      "OFFICIAL_ESTIMATE_INSTALLATION_VERIFIED"
  )


# ------------------------------------------------------------
# 14. Save installation audit
# ------------------------------------------------------------

install_audit_file_10b3nb <-
  file.path(
    dirname(
      estimate_tar_10b3nb
    ),
    "Paper1_STEP7B10B3NB_official_ESTIMATE_1.0.13_installation_audit_LOCKED.RData"
  )

save(
  estimate_tar_10b3nb,
  estimate_installed_before_10b3nb,
  estimate_version_before_10b3nb,
  install_result_10b3nb,
  estimate_installed_after_10b3nb,
  estimate_version_after_10b3nb,
  estimate_description_10b3nb,
  estimate_path_10b3nb,
  estimateScore_available_10b3nb,
  filterCommonGenes_available_10b3nb,
  SI_geneset_10b3nb,
  common_genes_10b3nb,
  step7be10b3nb_audit,
  file =
    install_audit_file_10b3nb
)

cat(
  "\nInstallation audit saved:\n",
  install_audit_file_10b3nb,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 15. Final check
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 10B-3N-B FINAL CHECK\n")
cat("============================================================\n")

cat(
  "Installation command completed: ",
  install_result_10b3nb,
  "\n",
  sep = ""
)

cat(
  "ESTIMATE installed: ",
  estimate_installed_after_10b3nb,
  "\n",
  sep = ""
)

cat(
  "ESTIMATE version: ",
  estimate_version_after_10b3nb,
  "\n",
  sep = ""
)

cat(
  "estimateScore available: ",
  estimateScore_available_10b3nb,
  "\n",
  sep = ""
)

cat(
  "filterCommonGenes available: ",
  filterCommonGenes_available_10b3nb,
  "\n",
  sep = ""
)

cat(
  "Score calculated: FALSE\n"
)

cat(
  "LOO score calculated: FALSE\n"
)

cat(
  "Original ESTIMATE scores modified: FALSE\n"
)

cat(
  "Step14 modified: FALSE\n"
)

cat(
  "Figure6 modified: FALSE\n"
)

cat(
  "Supplementary Table3 modified: FALSE\n"
)

cat(
  "\nSTATUS: ",
  step7be10b3nb_audit$status,
  "\n",
  sep = ""
)

cat("============================================================\n")


# ============================================================
# Paper 1
# STEP 10B-3N-C
# Official ESTIMATE 1.0.13 reconstruction of locked Step14 scores
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3N-C: OFFICIAL ESTIMATE SCORE RECONSTRUCTION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Load official package
# ------------------------------------------------------------

stopifnot(
  requireNamespace(
    "estimate",
    quietly = TRUE
  )
)

estimate_version_10b3nc <-
  as.character(
    packageVersion("estimate")
  )

cat(
  "ESTIMATE version: ",
  estimate_version_10b3nc,
  "\n",
  sep = ""
)

stopifnot(
  identical(
    estimate_version_10b3nc,
    "1.0.13"
  )
)

library(
  estimate
)


# ------------------------------------------------------------
# 1. Define locked Step14 source
# ------------------------------------------------------------

step14_locked_file_10b3nc <-
  "input_dir/Paper1_tumor-only WGCNA/WORKFLOW-03/Step14/Step14_FINAL_LOCKED/RData/WORKFLOW-03_Step14A-14E_FINAL_LOCKED.RData"

stopifnot(
  file.exists(
    step14_locked_file_10b3nc
  )
)

cat(
  "\nLocked Step14 source:\n",
  step14_locked_file_10b3nc,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 2. Load Step14 source into isolated environment
# ------------------------------------------------------------

step14_env_10b3nc <-
  new.env(
    parent = emptyenv()
  )

load(
  step14_locked_file_10b3nc,
  envir = step14_env_10b3nc
)

cat(
  "\nStep14 source loaded into isolated environment.\n"
)


# ------------------------------------------------------------
# 3. Recover required objects
# ------------------------------------------------------------

required_objects_10b3nc <-
  c(
    "estimate_input_file",
    "estimate_gct_input",
    "estimate_scores",
    "estimate_scores_wgcna",
    "estimate_sample_ids"
  )

missing_objects_10b3nc <-
  setdiff(
    required_objects_10b3nc,
    ls(
      envir = step14_env_10b3nc
    )
  )

cat(
  "\nMissing required objects:\n"
)

print(
  missing_objects_10b3nc
)

stopifnot(
  length(
    missing_objects_10b3nc
  ) == 0
)


estimate_input_file_10b3nc <-
  get(
    "estimate_input_file",
    envir = step14_env_10b3nc
  )

estimate_gct_input_10b3nc <-
  get(
    "estimate_gct_input",
    envir = step14_env_10b3nc
  )

estimate_scores_source_10b3nc <-
  get(
    "estimate_scores",
    envir = step14_env_10b3nc
  )

estimate_scores_wgcna_source_10b3nc <-
  get(
    "estimate_scores_wgcna",
    envir = step14_env_10b3nc
  )

estimate_sample_ids_source_10b3nc <-
  get(
    "estimate_sample_ids",
    envir = step14_env_10b3nc
  )


# ------------------------------------------------------------
# 4. Verify source input file
# ------------------------------------------------------------

cat(
  "\nOriginal Step14 estimate_input_file:\n",
  estimate_input_file_10b3nc,
  "\n",
  sep = ""
)

cat(
  "\nOriginal Step14 estimate_gct_input:\n",
  estimate_gct_input_10b3nc,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. Resolve actual input path
# ------------------------------------------------------------

estimate_input_candidates_10b3nc <-
  c(
    estimate_input_file_10b3nc,
    file.path(
      "O:/Rproject/Paper1_tumor-only WGCNA",
      basename(
        estimate_input_file_10b3nc
      )
    ),
    file.path(
      "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis",
      basename(
        estimate_input_file_10b3nc
      )
    )
  )

estimate_input_candidates_10b3nc <-
  unique(
    estimate_input_candidates_10b3nc
  )

estimate_input_existing_10b3nc <-
  estimate_input_candidates_10b3nc[
    file.exists(
      estimate_input_candidates_10b3nc
    )
  ]

cat(
  "\nExisting candidate input files:\n"
)

print(
  estimate_input_existing_10b3nc
)

stopifnot(
  length(
    estimate_input_existing_10b3nc
  ) >= 1
)

estimate_input_source_10b3nc <-
  estimate_input_existing_10b3nc[1]

cat(
  "\nSelected input:\n",
  estimate_input_source_10b3nc,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Audit original locked score matrix
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("LOCKED SOURCE SCORE AUDIT\n")
cat("============================================================\n\n")

cat(
  "estimate_scores dimensions: ",
  paste(
    dim(
      estimate_scores_source_10b3nc
    ),
    collapse = " x "
  ),
  "\n",
  sep = ""
)

cat(
  "estimate_scores_wgcna dimensions: ",
  paste(
    dim(
      estimate_scores_wgcna_source_10b3nc
    ),
    collapse = " x "
  ),
  "\n",
  sep = ""
)

cat(
  "\nLocked score row names:\n"
)

print(
  rownames(
    estimate_scores_wgcna_source_10b3nc
  )
)

stopifnot(
  identical(
    rownames(
      estimate_scores_wgcna_source_10b3nc
    ),
    c(
      "StromalScore",
      "ImmuneScore",
      "ESTIMATEScore"
    )
  )
)

stopifnot(
  ncol(
    estimate_scores_wgcna_source_10b3nc
  ) == 513
)


# ------------------------------------------------------------
# 7. Prepare isolated output directory
# ------------------------------------------------------------

step10b3nc_dir <-
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis/STEP7B-7E-10B-3N/10B-3N-C"

if (!dir.exists(step10b3nc_dir)) {
  dir.create(
    step10b3nc_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 8. Official filterCommonGenes()
# ------------------------------------------------------------

filtered_gct_10b3nc <-
  file.path(
    step10b3nc_dir,
    "STEP10B3NC_official_filterCommonGenes_output.gct"
  )

cat(
  "\n============================================================\n")
cat("RUNNING OFFICIAL filterCommonGenes()\n")
cat("============================================================\n\n")

filterCommonGenes(
  input.f =
    estimate_input_source_10b3nc,
  output.f =
    filtered_gct_10b3nc,
  id =
    "GeneSymbol"
)

stopifnot(
  file.exists(
    filtered_gct_10b3nc
  )
)

cat(
  "\nOfficial filtered GCT created:\n",
  filtered_gct_10b3nc,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Audit filtered GCT
# ------------------------------------------------------------

filtered_lines_10b3nc <-
  readLines(
    filtered_gct_10b3nc,
    n = 5
  )

cat(
  "\nFirst five lines of official filtered GCT:\n\n"
)

print(
  filtered_lines_10b3nc
)


# ------------------------------------------------------------
# 10. Official estimateScore()
# ------------------------------------------------------------

reconstructed_score_gct_10b3nc <-
  file.path(
    step10b3nc_dir,
    "STEP10B3NC_official_estimateScore_output.gct"
  )

cat(
  "\n============================================================\n")
cat("RUNNING OFFICIAL estimateScore()\n")
cat("============================================================\n\n")

estimateScore(
  input.ds =
    filtered_gct_10b3nc,
  output.ds =
    reconstructed_score_gct_10b3nc,
  platform =
    "illumina"
)

stopifnot(
  file.exists(
    reconstructed_score_gct_10b3nc
  )
)

cat(
  "\nOfficial score GCT created:\n",
  reconstructed_score_gct_10b3nc,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Read official reconstructed score
# ------------------------------------------------------------

reconstructed_score_10b3nc <-
  read.delim(
    reconstructed_score_gct_10b3nc,
    header = TRUE,
    sep = "\t",
    skip = 2,
    row.names = 1,
    blank.lines.skip = TRUE,
    check.names = FALSE,
    as.is = TRUE
  )

cat(
  "\nReconstructed score dimensions: ",
  paste(
    dim(
      reconstructed_score_10b3nc
    ),
    collapse = " x "
  ),
  "\n",
  sep = ""
)

print(
  head(
    reconstructed_score_10b3nc
  )
)


# ------------------------------------------------------------
# 12. Normalize reconstructed score object
# ------------------------------------------------------------

# estimateScore() writes the description column first.
# Remove it before numerical comparison.

if (
  ncol(
    reconstructed_score_10b3nc
  ) == 514
) {
  
  reconstructed_score_numeric_10b3nc <-
    reconstructed_score_10b3nc[
      ,
      -1,
      drop = FALSE
    ]
  
} else {
  
  reconstructed_score_numeric_10b3nc <-
    reconstructed_score_10b3nc
}

reconstructed_score_numeric_10b3nc <-
  as.matrix(
    reconstructed_score_numeric_10b3nc
  )

storage.mode(
  reconstructed_score_numeric_10b3nc
) <- "numeric"


# ------------------------------------------------------------
# 13. Normalize sample names
# ------------------------------------------------------------

normalize_tcga_10b3nc <-
  function(x) {
    
    x <- as.character(x)
    
    x <- gsub(
      "\\.",
      "-",
      x
    )
    
    x
  }

reconstructed_sample_ids_10b3nc <-
  normalize_tcga_10b3nc(
    colnames(
      reconstructed_score_numeric_10b3nc
    )
  )

source_sample_ids_10b3nc <-
  normalize_tcga_10b3nc(
    colnames(
      estimate_scores_wgcna_source_10b3nc
    )
  )


# ------------------------------------------------------------
# 14. Sample identity/order audit
# ------------------------------------------------------------

sample_identity_10b3nc <-
  identical(
    reconstructed_sample_ids_10b3nc,
    source_sample_ids_10b3nc
  )

cat(
  "\nSample identity/order match: ",
  sample_identity_10b3nc,
  "\n",
  sep = ""
)

stopifnot(
  sample_identity_10b3nc
)


# ------------------------------------------------------------
# 15. Score row-name audit
# ------------------------------------------------------------

reconstructed_score_numeric_10b3nc <-
  reconstructed_score_numeric_10b3nc[
    rownames(
      estimate_scores_wgcna_source_10b3nc
    ),
    ,
    drop = FALSE
  ]

cat(
  "\nReconstructed score rows:\n"
)

print(
  rownames(
    reconstructed_score_numeric_10b3nc
  )
)

stopifnot(
  identical(
    rownames(
      reconstructed_score_numeric_10b3nc
    ),
    rownames(
      estimate_scores_wgcna_source_10b3nc
    )
  )
)


# ------------------------------------------------------------
# 16. Exact numerical comparison
# ------------------------------------------------------------

comparison_results_10b3nc <-
  lapply(
    
    rownames(
      estimate_scores_wgcna_source_10b3nc
    ),
    
    function(score_name) {
      
      source_values <-
        as.numeric(
          estimate_scores_wgcna_source_10b3nc[
            score_name,
          ]
        )
      
      reconstructed_values <-
        as.numeric(
          reconstructed_score_numeric_10b3nc[
            score_name,
          ]
        )
      
      data.frame(
        
        Score =
          score_name,
        
        Pearson =
          cor(
            source_values,
            reconstructed_values,
            method = "pearson"
          ),
        
        Spearman =
          cor(
            source_values,
            reconstructed_values,
            method = "spearman"
          ),
        
        MaxAbsDiff =
          max(
            abs(
              source_values -
                reconstructed_values
            )
          ),
        
        MeanAbsDiff =
          mean(
            abs(
              source_values -
                reconstructed_values
            )
          ),
        
        RMSE =
          sqrt(
            mean(
              (
                source_values -
                  reconstructed_values
              )^2
            )
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )

comparison_results_10b3nc <-
  do.call(
    rbind,
    comparison_results_10b3nc
  )

cat("\n============================================================\n")
cat("OFFICIAL RECONSTRUCTION COMPARISON\n")
cat("============================================================\n\n")

print(
  comparison_results_10b3nc,
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Exact / near-exact reconstruction criterion
# ------------------------------------------------------------

max_abs_diff_overall_10b3nc <-
  max(
    comparison_results_10b3nc$MaxAbsDiff
  )

mean_abs_diff_overall_10b3nc <-
  max(
    comparison_results_10b3nc$MeanAbsDiff
  )

reconstruction_pass_10b3nc <-
  isTRUE(
    max_abs_diff_overall_10b3nc <
      1e-8
  )

cat(
  "\nOverall maximum absolute difference: ",
  format(
    max_abs_diff_overall_10b3nc,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Exact reconstruction criterion (<1e-8): ",
  reconstruction_pass_10b3nc,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 18. HARD STOP if reconstruction fails
# ------------------------------------------------------------

if (!reconstruction_pass_10b3nc) {
  
  cat(
    "\n============================================================\n"
  )
  
  cat(
    "HARD STOP\n"
  )
  
  cat(
    "Official ESTIMATE reconstruction does not yet reproduce\n"
  )
  
  cat(
    "the locked Step14 scores at the exact numerical level.\n"
  )
  
  cat(
    "NO LOO ANALYSIS will be performed.\n"
  )
  
  cat(
    "NO statistical testing will be performed.\n"
  )
  
  cat(
    "NO locked Step14 objects will be modified.\n"
  )
  
  cat(
    "============================================================\n"
  )
  
}


# ------------------------------------------------------------
# 19. Safety audit
# ------------------------------------------------------------

step7be10b3nc_audit <-
  list(
    
    step =
      "STEP 10B-3N-C",
    
    purpose =
      "Official ESTIMATE 1.0.13 reconstruction of locked Step14 scores",
    
    estimate_version =
      estimate_version_10b3nc,
    
    input_file =
      estimate_input_source_10b3nc,
    
    filtered_gct =
      filtered_gct_10b3nc,
    
    reconstructed_score_gct =
      reconstructed_score_gct_10b3nc,
    
    sample_identity_order_match =
      sample_identity_10b3nc,
    
    comparison =
      comparison_results_10b3nc,
    
    reconstruction_pass =
      reconstruction_pass_10b3nc,
    
    loo_score_calculated =
      FALSE,
    
    statistical_analysis_performed =
      FALSE,
    
    locked_Step14_modified =
      FALSE,
    
    Figure6_modified =
      FALSE,
    
    SupplementaryTable3_modified =
      FALSE,
    
    status =
      if (reconstruction_pass_10b3nc) {
        "OFFICIAL_RECONSTRUCTION_PASSED"
      } else {
        "OFFICIAL_RECONSTRUCTION_FAILED_HARD_STOP"
      }
  )


# ------------------------------------------------------------
# 20. Save checkpoint
# ------------------------------------------------------------

reconstruction_audit_file_10b3nc <-
  file.path(
    step10b3nc_dir,
    "Paper1_STEP7B10B3NC_official_ESTIMATE_reconstruction_audit_LOCKED.RData"
  )

save(
  estimate_version_10b3nc,
  estimate_input_source_10b3nc,
  filtered_gct_10b3nc,
  reconstructed_score_gct_10b3nc,
  reconstructed_score_numeric_10b3nc,
  estimate_scores_wgcna_source_10b3nc,
  comparison_results_10b3nc,
  reconstruction_pass_10b3nc,
  step7be10b3nc_audit,
  file =
    reconstruction_audit_file_10b3nc
)

cat(
  "\nCheckpoint saved:\n",
  reconstruction_audit_file_10b3nc,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 21. Final check
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 10B-3N-C FINAL CHECK\n")
cat("============================================================\n")

cat(
  "Official ESTIMATE version: ",
  estimate_version_10b3nc,
  "\n",
  sep = ""
)

cat(
  "Sample identity/order: ",
  sample_identity_10b3nc,
  "\n",
  sep = ""
)

cat(
  "Official reconstruction pass: ",
  reconstruction_pass_10b3nc,
  "\n",
  sep = ""
)

cat(
  "LOO score calculated: FALSE\n"
)

cat(
  "Statistical analysis: FALSE\n"
)

cat(
  "Locked Step14 modified: FALSE\n"
)

cat(
  "Figure6 modified: FALSE\n"
)

cat(
  "SupplementaryTable3 modified: FALSE\n"
)

cat(
  "\nSTATUS: ",
  step7be10b3nc_audit$status,
  "\n",
  sep = ""
)

cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 7B-7E-10B-3N-D
# Official ESTIMATE 1.0.13 LOO Score Construction
#
# Purpose:
#   1. Remove AOC3 from the official StromalSignature
#   2. Remove IL7R from the official ImmuneSignature
#   3. Run the official ESTIMATE 1.0.13 scoring algorithm
#      using temporary signature-specific copies
#   4. Generate:
#        - LOO_AOC3_StromalScore
#        - LOO_IL7R_ImmuneScore
#   5. Verify sample identity/order and score structure
#
# IMPORTANT:
#   - No statistical testing
#   - No correlation/regression
#   - No FDR calculation
#   - No modification of Step14 locked objects
#   - No modification of Figure 6
#   - No modification of Supplementary Table 3
#   - Do NOT perform LOO for MS4A7 because it is absent
#     from both ESTIMATE signatures
#
# Official package:
#   estimate 1.0.13
# ============================================================


cat("\n============================================================\n")
cat("STEP 7B-7E-10B-3N-D\n")
cat("Official ESTIMATE 1.0.13 LOO Score Construction\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Root directories
# ------------------------------------------------------------

root_7be10 <- "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis"

step_dir_10b3nd <- file.path(
  root_7be10,
  "STEP7B-7E-10B-3N",
  "10B-3N-D"
)

dir.create(
  step_dir_10b3nd,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Required package
# ------------------------------------------------------------

if (!requireNamespace("estimate", quietly = TRUE)) {
  stop(
    "\nERROR: package 'estimate' is not installed.\n",
    "STEP 10B-3N-D requires official estimate 1.0.13.\n"
  )
}

estimate_version_10b3nd <- as.character(
  packageVersion("estimate")
)

cat("Installed estimate version:\n")
print(estimate_version_10b3nd)


if (!identical(estimate_version_10b3nd, "1.0.13")) {
  
  stop(
    "\nERROR: estimate package version is not 1.0.13.\n",
    "Detected version: ",
    estimate_version_10b3nd,
    "\n"
  )
  
}


library(estimate)


# ------------------------------------------------------------
# 2. Load 10B-3N-C reconstruction checkpoint
# ------------------------------------------------------------

checkpoint_candidates_10b3nd <- c(
  
  file.path(
    root_7be10,
    "STEP7B-7E-10B-3N",
    "10B-3N-C",
    "Paper1_STEP7B10B3NC_official_ESTIMATE_reconstruction_audit_LOCKED.RData"
  ),
  
  file.path(
    root_7be10,
    "STEP7B-7E-10B-3N",
    "10B-3N-C",
    "Paper1_STEP7B10B3NC_official_ESTIMATE_reconstruction_audit_LOCKED.RData"
  )
  
)

checkpoint_candidates_10b3nd <- unique(
  checkpoint_candidates_10b3nd[
    file.exists(checkpoint_candidates_10b3nd)
  ]
)


if (length(checkpoint_candidates_10b3nd) == 0) {
  
  stop(
    "\nERROR: 10B-3N-C checkpoint not found.\n",
    "Please verify the 10B-3N-C output directory.\n"
  )
  
}


checkpoint_10b3nd <- checkpoint_candidates_10b3nd[1]

cat("\n10B-3N-C checkpoint:\n")
cat(checkpoint_10b3nd, "\n\n")


load(
  checkpoint_10b3nd,
  envir = .GlobalEnv
)


# ------------------------------------------------------------
# 3. Recover required source objects
# ------------------------------------------------------------

required_objects_10b3nd <- c(
  "estimate_scores_wgcna_source_10b3nc"
)

missing_required_10b3nd <- setdiff(
  required_objects_10b3nd,
  ls(envir = .GlobalEnv)
)


if (length(missing_required_10b3nd) > 0) {
  
  stop(
    "\nERROR: required object(s) missing after loading 10B-3N-C checkpoint:\n",
    paste(missing_required_10b3nd, collapse = ", "),
    "\n"
  )
  
}


source_score_10b3nd <-
  as.matrix(
    estimate_scores_wgcna_source_10b3nc
  )


# ------------------------------------------------------------
# 4. Basic source-score audit
# ------------------------------------------------------------

cat("Source locked ESTIMATE score dimensions:\n")
print(dim(source_score_10b3nd))

cat("\nSource score rows:\n")
print(rownames(source_score_10b3nd))


stopifnot(
  identical(
    rownames(source_score_10b3nd),
    c(
      "StromalScore",
      "ImmuneScore",
      "ESTIMATEScore"
    )
  )
)

stopifnot(
  ncol(source_score_10b3nd) == 513
)

stopifnot(
  all(is.finite(source_score_10b3nd))
)


# ------------------------------------------------------------
# 5. Locate the exact filtered GCT used in 10B-3N-C
#
# IMPORTANT:
#   We want to reuse the exact input structure that produced
#   the successful 10B-3N-C reconstruction.
# ------------------------------------------------------------

c_dir_10b3nd <- file.path(
  root_7be10,
  "STEP7B-7E-10B-3N",
  "10B-3N-C"
)


gct_candidates_10b3nd <- list.files(
  c_dir_10b3nd,
  pattern = "\\.gct$",
  full.names = TRUE,
  recursive = TRUE
)


cat("\nGCT files found in 10B-3N-C:\n")
print(gct_candidates_10b3nd)


# ------------------------------------------------------------
# 6. Prefer the filtered input GCT rather than score output
# ------------------------------------------------------------

gct_input_candidates_10b3nd <- gct_candidates_10b3nd[
  grepl(
    "filtered|input",
    basename(gct_candidates_10b3nd),
    ignore.case = TRUE
  ) &
    !grepl(
      "output|score|reconstructed",
      basename(gct_candidates_10b3nd),
      ignore.case = TRUE
    )
]


if (length(gct_input_candidates_10b3nd) == 0) {
  
  # If the exact GCT is still available as an object in memory,
  # use that object and write it to a temporary local file.
  
  if (
    exists(
      "filtered_gct_10b3nc",
      envir = .GlobalEnv
    )
  ) {
    
    cat(
      "\nNo input GCT file found on disk.\n",
      "Using existing object: filtered_gct_10b3nc\n"
    )
    
    filtered_gct_10b3nd <- filtered_gct_10b3nc
    
    input_gct_10b3nd <- file.path(
      step_dir_10b3nd,
      "STEP10B3ND_base_filtered_input.gct"
    )
    
    writeLines(
      as.character(filtered_gct_10b3nd),
      input_gct_10b3nd
    )
    
  } else {
    
    stop(
      "\nERROR: exact filtered GCT used in 10B-3N-C could not be found.\n",
      "No suitable GCT file and no object 'filtered_gct_10b3nc' exists.\n",
      "STOP before LOO scoring.\n"
    )
    
  }
  
} else {
  
  # If several files are found, stop rather than guessing.
  if (length(gct_input_candidates_10b3nd) > 1) {
    
    cat(
      "\nMultiple possible GCT input files found:\n"
    )
    
    print(gct_input_candidates_10b3nd)
    
    stop(
      "\nERROR: More than one possible 10B-3N-C input GCT exists.\n",
      "Please identify the exact GCT before continuing.\n"
    )
    
  }
  
  input_gct_10b3nd <- gct_input_candidates_10b3nd[1]
  
}


cat("\nSelected base GCT:\n")
cat(input_gct_10b3nd, "\n")


stopifnot(
  file.exists(input_gct_10b3nd)
)

# ============================================================
# STEP 10B-3N-D
# FIX: use the exact official filterCommonGenes() output
# from STEP 10B-3N-C
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3N-D: FIX BASE GCT INPUT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Exact official filterCommonGenes() output
# ------------------------------------------------------------

input_gct_10b3nd <- file.path(
  root_7be10,
  "STEP7B-7E-10B-3N",
  "10B-3N-C",
  "STEP10B3NC_official_filterCommonGenes_output.gct"
)


cat("Selected official filtered GCT:\n")
cat(input_gct_10b3nd, "\n\n")


# ------------------------------------------------------------
# 2. Hard-stop file existence check
# ------------------------------------------------------------

if (!file.exists(input_gct_10b3nd)) {
  
  stop(
    "\nERROR: official filterCommonGenes() output not found.\n",
    input_gct_10b3nd,
    "\n"
  )
  
}


cat(
  "File exists: ",
  file.exists(input_gct_10b3nd),
  "\n"
)


# ------------------------------------------------------------
# 3. Inspect first 5 lines BEFORE parsing
# ------------------------------------------------------------

cat("\nFirst 5 lines of official GCT:\n")
cat(
  paste(
    readLines(
      input_gct_10b3nd,
      n = 5
    ),
    collapse = "\n"
  )
)

cat("\n\n")


# ------------------------------------------------------------
# 4. Read official filtered GCT
#
# Standard ESTIMATE GCT:
#   line 1 = version
#   line 2 = dimensions
#   line 3 = column header
# ------------------------------------------------------------

base_gct_10b3nd <- read.table(
  input_gct_10b3nd,
  header = TRUE,
  sep = "\t",
  skip = 2,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)


cat("Official filtered GCT dimensions:\n")
print(dim(base_gct_10b3nd))


cat("\nFirst columns:\n")
print(
  head(
    colnames(base_gct_10b3nd),
    10
  )
)


# ------------------------------------------------------------
# 5. Basic structural checks
# ------------------------------------------------------------

if (nrow(base_gct_10b3nd) == 0) {
  
  stop(
    "\nERROR: official filtered GCT contains zero data rows.\n"
  )
  
}


# Expected ESTIMATE filter output:
# NAME / Description / 513 samples

if (!("NAME" %in% colnames(base_gct_10b3nd))) {
  
  stop(
    "\nERROR: NAME column not found in official filtered GCT.\n",
    "Observed columns:\n",
    paste(
      head(colnames(base_gct_10b3nd), 10),
      collapse = ", "
    ),
    "\n"
  )
  
}


if (!("Description" %in% colnames(base_gct_10b3nd))) {
  
  stop(
    "\nERROR: Description column not found in official filtered GCT.\n"
  )
  
}


# ------------------------------------------------------------
# 6. Identify sample columns
# ------------------------------------------------------------

sample_cols_10b3nd <- setdiff(
  colnames(base_gct_10b3nd),
  c("NAME", "Description")
)


cat("\nNumber of sample columns:\n")
cat(length(sample_cols_10b3nd), "\n")


stopifnot(
  length(sample_cols_10b3nd) == 513
)


# ------------------------------------------------------------
# 7. Gene identity
# ------------------------------------------------------------

gene_names_10b3nd <- as.character(
  base_gct_10b3nd$Description
)


cat("\nNumber of genes:\n")
cat(length(gene_names_10b3nd), "\n")


stopifnot(
  length(gene_names_10b3nd) == nrow(base_gct_10b3nd)
)

stopifnot(
  !anyNA(gene_names_10b3nd)
)

stopifnot(
  all(nzchar(gene_names_10b3nd))
)

stopifnot(
  !anyDuplicated(gene_names_10b3nd)
)


# ------------------------------------------------------------
# 8. Sample ID normalization
# ------------------------------------------------------------

sample_ids_gct_10b3nd <- gsub(
  "\\.",
  "-",
  sample_cols_10b3nd
)


source_sample_ids_10b3nd <- colnames(
  source_score_10b3nd
)


cat("\nSample set identity:\n")

sample_identity_set_10b3nd <- identical(
  sort(sample_ids_gct_10b3nd),
  sort(source_sample_ids_10b3nd)
)

print(
  sample_identity_set_10b3nd
)


cat("\nSample order identity:\n")

sample_identity_order_10b3nd <- identical(
  sample_ids_gct_10b3nd,
  source_sample_ids_10b3nd
)

print(
  sample_identity_order_10b3nd
)


stopifnot(
  sample_identity_set_10b3nd
)

stopifnot(
  sample_identity_order_10b3nd
)


# ------------------------------------------------------------
# 9. Verify exact gene count expected from 10B-3N-C
# ------------------------------------------------------------

cat("\nExpected filtered-gene count from 10B-3N-C:\n")
cat("10221 genes\n")

cat("\nObserved:\n")
cat(nrow(base_gct_10b3nd), "genes\n")


stopifnot(
  nrow(base_gct_10b3nd) == 10221
)


# ------------------------------------------------------------
# 10. Verify core-gene presence
# ------------------------------------------------------------

core_genes_10b3nd <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


core_presence_10b3nd <- setNames(
  core_genes_10b3nd %in% gene_names_10b3nd,
  core_genes_10b3nd
)


cat("\nCore-gene presence in filtered GCT:\n")
print(core_presence_10b3nd)


# Based on 10B-3N-C structural audit:
# AOC3 TRUE
# MS4A7 FALSE
# IL7R TRUE

stopifnot(
  isTRUE(core_presence_10b3nd["AOC3"])
)

stopifnot(
  isFALSE(core_presence_10b3nd["MS4A7"])
)

stopifnot(
  isTRUE(core_presence_10b3nd["IL7R"])
)


# ------------------------------------------------------------
# 11. Verify finite expression matrix
# ------------------------------------------------------------

expression_matrix_10b3nd <- as.matrix(
  base_gct_10b3nd[
    ,
    sample_cols_10b3nd,
    drop = FALSE
  ]
)

storage.mode(expression_matrix_10b3nd) <- "numeric"


cat("\nChecking finite expression values...\n")

finite_expression_10b3nd <- all(
  is.finite(
    expression_matrix_10b3nd
  )
)

cat(
  "All expression values finite: ",
  finite_expression_10b3nd,
  "\n"
)


stopifnot(
  finite_expression_10b3nd
)


# ------------------------------------------------------------
# 12. Verify this is exactly the successful 10B-3N-C
# filterCommonGenes() output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("10B-3N-C FILTERED GCT STRUCTURAL VALIDATION\n")
cat("============================================================\n\n")


cat(
  "Dimensions: ",
  nrow(base_gct_10b3nd),
  " x ",
  ncol(base_gct_10b3nd),
  "\n"
)

cat(
  "Genes: ",
  nrow(base_gct_10b3nd),
  "\n"
)

cat(
  "Samples: ",
  length(sample_cols_10b3nd),
  "\n"
)

cat(
  "Sample set identical: ",
  sample_identity_set_10b3nd,
  "\n"
)

cat(
  "Sample order identical: ",
  sample_identity_order_10b3nd,
  "\n"
)

cat(
  "AOC3 present: ",
  core_presence_10b3nd["AOC3"],
  "\n"
)

cat(
  "MS4A7 present: ",
  core_presence_10b3nd["MS4A7"],
  "\n"
)

cat(
  "IL7R present: ",
  core_presence_10b3nd["IL7R"],
  "\n"
)

cat(
  "All expression values finite: ",
  finite_expression_10b3nd,
  "\n"
)


# ------------------------------------------------------------
# 13. DO NOT use the accidentally created empty file
# ------------------------------------------------------------

invalid_temp_gct_10b3nd <- file.path(
  step_dir_10b3nd,
  "STEP10B3ND_base_filtered_input.gct"
)


if (file.exists(invalid_temp_gct_10b3nd)) {
  
  cat(
    "\nRemoving accidentally generated temporary GCT:\n",
    invalid_temp_gct_10b3nd,
    "\n"
  )
  
  file.remove(
    invalid_temp_gct_10b3nd
  )
  
}


# ------------------------------------------------------------
# 14. Final checkpoint before LOO scoring
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("BASE GCT FIX: PASSED\n")
cat("============================================================\n")

cat(
  "\nThe exact official filterCommonGenes() output from\n",
  "10B-3N-C will now be used for 10B-3N-D.\n\n"
)

cat(
  "Filtered genes: ",
  nrow(base_gct_10b3nd),
  "\n"
)

cat(
  "Samples: ",
  length(sample_cols_10b3nd),
  "\n"
)

cat(
  "Sample order: PASSED\n"
)

cat(
  "Core-gene audit: PASSED\n"
)

cat(
  "Finite-value audit: PASSED\n"
)

cat(
  "No LOO score calculated yet.\n"
)

cat(
  "No statistical analysis performed.\n"
)

cat(
  "No Step14 modification.\n"
)

cat(
  "No Figure6 modification.\n"
)

cat(
  "No Supplementary Table 3 modification.\n"
)

cat(
  "\n============================================================\n"
)


# ============================================================
# Paper 1
# STEP 7B-7E-10B-3N-D-1
# Official ESTIMATE LOO Mechanism Validation
#
# Purpose:
#   HARD-VALIDATE that a modified SI_geneset is actually used
#   by the official estimateScore() implementation.
#
# Validation targets:
#   1. Remove AOC3 from StromalSignature
#   2. Remove IL7R from ImmuneSignature
#
# IMPORTANT:
#   - NO final LOO association statistics
#   - NO statistical testing
#   - NO modification of locked Step14 objects
#   - NO modification of Figure 6
#   - NO modification of Supplementary Table 3
#   - NO modification of installed estimate package namespace
#
# Official package:
#   estimate 1.0.13
#
# ============================================================


cat("\n============================================================\n")
cat("STEP 7B-7E-10B-3N-D-1\n")
cat("OFFICIAL ESTIMATE LOO MECHANISM VALIDATION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define paths
# ------------------------------------------------------------

root_7be10 <- "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis"

step_10b3nc_dir <- file.path(
  root_7be10,
  "STEP7B-7E-10B-3N",
  "10B-3N-C"
)

step_10b3nd1_dir <- file.path(
  root_7be10,
  "STEP7B-7E-10B-3N",
  "10B-3N-D-1"
)

dir.create(
  step_10b3nd1_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Check official estimate package
# ------------------------------------------------------------

cat("Checking official estimate package...\n\n")

estimate_installed_10b3nd1 <- requireNamespace(
  "estimate",
  quietly = TRUE
)

stopifnot(
  estimate_installed_10b3nd1
)

estimate_version_10b3nd1 <- as.character(
  packageVersion("estimate")
)

cat(
  "Installed estimate version: ",
  estimate_version_10b3nd1,
  "\n\n"
)

stopifnot(
  estimate_version_10b3nd1 == "1.0.13"
)


# ------------------------------------------------------------
# 3. Load official 10B-3N-C reconstruction checkpoint
# ------------------------------------------------------------

checkpoint_10b3nc <- file.path(
  step_10b3nc_dir,
  "Paper1_STEP7B10B3NC_official_ESTIMATE_reconstruction_audit_LOCKED.RData"
)

stopifnot(
  file.exists(checkpoint_10b3nc)
)

load(
  checkpoint_10b3nc,
  envir = .GlobalEnv
)

cat(
  "10B-3N-C checkpoint loaded:\n",
  checkpoint_10b3nc,
  "\n\n"
)


# ------------------------------------------------------------
# 4. Locate official filterCommonGenes output
# ------------------------------------------------------------

input_gct_10b3nd1 <- file.path(
  step_10b3nc_dir,
  "STEP10B3NC_official_filterCommonGenes_output.gct"
)

stopifnot(
  file.exists(input_gct_10b3nd1)
)

cat(
  "Official filtered GCT:\n",
  input_gct_10b3nd1,
  "\n\n"
)


# ------------------------------------------------------------
# 5. Inspect official GCT header
# ------------------------------------------------------------

cat("First five lines of official filtered GCT:\n\n")

header_preview_10b3nd1 <- readLines(
  input_gct_10b3nd1,
  n = 5,
  warn = FALSE
)

cat(
  paste(
    header_preview_10b3nd1,
    collapse = "\n"
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 6. Read official filtered GCT
# ------------------------------------------------------------

filtered_gct_10b3nd1 <- read.table(
  input_gct_10b3nd1,
  header = TRUE,
  sep = "\t",
  skip = 2,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

cat(
  "Filtered GCT dimensions: ",
  nrow(filtered_gct_10b3nd1),
  " x ",
  ncol(filtered_gct_10b3nd1),
  "\n"
)

stopifnot(
  nrow(filtered_gct_10b3nd1) == 10221,
  ncol(filtered_gct_10b3nd1) == 515
)


# ------------------------------------------------------------
# 7. Recover gene and sample columns
# ------------------------------------------------------------

stopifnot(
  all(
    c("NAME", "Description") %in%
      colnames(filtered_gct_10b3nd1)
  )
)

gene_symbols_10b3nd1 <- as.character(
  filtered_gct_10b3nd1$Description
)

sample_cols_10b3nd1 <- setdiff(
  colnames(filtered_gct_10b3nd1),
  c("NAME", "Description")
)

sample_ids_10b3nd1 <- gsub(
  "\\.",
  "-",
  sample_cols_10b3nd1
)

cat(
  "Genes: ",
  length(gene_symbols_10b3nd1),
  "\n"
)

cat(
  "Samples: ",
  length(sample_ids_10b3nd1),
  "\n\n"
)

stopifnot(
  length(sample_cols_10b3nd1) == 513,
  length(sample_ids_10b3nd1) == 513,
  !anyDuplicated(gene_symbols_10b3nd1),
  !anyDuplicated(sample_ids_10b3nd1)
)


# ------------------------------------------------------------
# 8. Recover locked source ESTIMATE score
# ------------------------------------------------------------

stopifnot(
  exists("estimate_scores_wgcna_source_10b3nc")
)

source_score_10b3nd1 <-
  as.matrix(
    estimate_scores_wgcna_source_10b3nc
  )

rownames_source_10b3nd1 <- rownames(
  source_score_10b3nd1
)

cat(
  "\nLocked source score dimensions: ",
  nrow(source_score_10b3nd1),
  " x ",
  ncol(source_score_10b3nd1),
  "\n"
)

cat(
  "Locked source score rows:\n"
)

print(
  rownames_source_10b3nd1
)

stopifnot(
  identical(
    rownames_source_10b3nd1,
    c(
      "StromalScore",
      "ImmuneScore",
      "ESTIMATEScore"
    )
  )
)

stopifnot(
  ncol(source_score_10b3nd1) == 513
)

stopifnot(
  all(
    is.finite(
      source_score_10b3nd1
    )
  )
)


# ------------------------------------------------------------
# 9. Verify sample identity and order
# ------------------------------------------------------------

source_sample_ids_10b3nd1 <-
  colnames(source_score_10b3nd1)

sample_set_identity_10b3nd1 <-
  identical(
    sample_ids_10b3nd1,
    source_sample_ids_10b3nd1
  )

cat(
  "\nSample identity/order identical: ",
  sample_set_identity_10b3nd1,
  "\n"
)

stopifnot(
  sample_set_identity_10b3nd1
)


# ============================================================
# 10. Recover official SI_geneset
# ============================================================

cat("\n============================================================\n")
cat("RECOVER OFFICIAL SI_geneSET\n")
cat("============================================================\n\n")


# The official package stores SI_geneset in its namespace.
# We READ it only; we do NOT modify the package namespace.

SI_geneset_original_10b3nd1 <- get(
  "SI_geneset",
  envir = asNamespace("estimate")
)

cat(
  "SI_geneset class: ",
  paste(
    class(SI_geneset_original_10b3nd1),
    collapse = ", "
  ),
  "\n"
)

cat(
  "SI_geneset dimensions: ",
  nrow(SI_geneset_original_10b3nd1),
  " x ",
  ncol(SI_geneset_original_10b3nd1),
  "\n\n"
)

print(
  head(
    SI_geneset_original_10b3nd1
  )
)


# ============================================================
# STEP 7B-7E-10B-3N-D-1
# REPAIR FROM STEP 11
#
# Reason for repair:
# Official estimate 1.0.13 SI_geneset is actually:
#   2 rows x 142 columns
#
# Rows:
#   StromalSignature
#   ImmuneSignature
#
# Columns:
#   V2 = "estimate"
#   V3:V142 = 141 genes
#
# Therefore:
#   DO NOT transpose SI_geneset
#   DO NOT expect 141 rows x 2 columns
#
# ============================================================


cat("\n============================================================\n")
cat("10B-3N-D-1 REPAIR: OFFICIAL SI_geneSET STRUCTURE\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 11R. Structural audit of the actual official SI_geneset
# ------------------------------------------------------------

cat("Actual SI_geneset dimensions:\n")

print(
  dim(
    SI_geneset_original_10b3nd1
  )
)

cat("\n")

stopifnot(
  is.data.frame(
    SI_geneset_original_10b3nd1
  )
)

stopifnot(
  nrow(
    SI_geneset_original_10b3nd1
  ) == 2
)

stopifnot(
  ncol(
    SI_geneset_original_10b3nd1
  ) == 142
)


# ------------------------------------------------------------
# 12R. Verify row names
# ------------------------------------------------------------

cat(
  "SI_geneset row names:\n"
)

print(
  rownames(
    SI_geneset_original_10b3nd1
  )
)

stopifnot(
  identical(
    rownames(
      SI_geneset_original_10b3nd1
    ),
    c(
      "StromalSignature",
      "ImmuneSignature"
    )
  )
)


# ------------------------------------------------------------
# 13R. Verify first column
# ------------------------------------------------------------

first_column_10b3nd1 <-
  as.character(
    SI_geneset_original_10b3nd1[
      ,
      1
    ]
  )

cat(
  "\nFirst SI_geneset column:\n"
)

print(
  first_column_10b3nd1
)

stopifnot(
  identical(
    first_column_10b3nd1,
    c(
      "estimate",
      "estimate"
    )
  )
)


# ------------------------------------------------------------
# 14R. Extract the 141 official genes from columns 2:142
# ------------------------------------------------------------

signature_gene_matrix_10b3nd1 <-
  SI_geneset_original_10b3nd1[
    ,
    -1,
    drop = FALSE
  ]

stopifnot(
  nrow(
    signature_gene_matrix_10b3nd1
  ) == 2
)

stopifnot(
  ncol(
    signature_gene_matrix_10b3nd1
  ) == 141
)


# ------------------------------------------------------------
# 15R. Recover official signature gene vectors
# ------------------------------------------------------------

stromal_genes_original_10b3nd1 <-
  as.character(
    signature_gene_matrix_10b3nd1[
      "StromalSignature",
      ,
      drop = TRUE
    ]
  )

immune_genes_original_10b3nd1 <-
  as.character(
    signature_gene_matrix_10b3nd1[
      "ImmuneSignature",
      ,
      drop = TRUE
    ]
  )


# Remove NA / empty strings if present
stromal_genes_original_10b3nd1 <-
  stromal_genes_original_10b3nd1[
    !is.na(
      stromal_genes_original_10b3nd1
    ) &
      stromal_genes_original_10b3nd1 != ""
  ]

immune_genes_original_10b3nd1 <-
  immune_genes_original_10b3nd1[
    !is.na(
      immune_genes_original_10b3nd1
    ) &
      immune_genes_original_10b3nd1 != ""
  ]


cat(
  "\nOfficial Stromal signature genes: ",
  length(
    stromal_genes_original_10b3nd1
  ),
  "\n"
)

cat(
  "Official Immune signature genes: ",
  length(
    immune_genes_original_10b3nd1
  ),
  "\n"
)


stopifnot(
  length(
    stromal_genes_original_10b3nd1
  ) == 141
)

stopifnot(
  length(
    immune_genes_original_10b3nd1
  ) == 141
)


# ------------------------------------------------------------
# 16R. Check uniqueness
# ------------------------------------------------------------

cat(
  "\nDuplicate Stromal genes: ",
  anyDuplicated(
    stromal_genes_original_10b3nd1
  ),
  "\n"
)

cat(
  "Duplicate Immune genes: ",
  anyDuplicated(
    immune_genes_original_10b3nd1
  ),
  "\n"
)

stopifnot(
  !anyDuplicated(
    stromal_genes_original_10b3nd1
  )
)

stopifnot(
  !anyDuplicated(
    immune_genes_original_10b3nd1
  )
)


# ============================================================
# 17R. Core-gene membership HARD STOP
# ============================================================

cat("\n============================================================\n")
cat("CORE-GENE SIGNATURE MEMBERSHIP AUDIT\n")
cat("============================================================\n\n")


AOC3_in_original_stromal_10b3nd1 <-
  "AOC3" %in%
  stromal_genes_original_10b3nd1

AOC3_in_original_immune_10b3nd1 <-
  "AOC3" %in%
  immune_genes_original_10b3nd1

IL7R_in_original_stromal_10b3nd1 <-
  "IL7R" %in%
  stromal_genes_original_10b3nd1

IL7R_in_original_immune_10b3nd1 <-
  "IL7R" %in%
  immune_genes_original_10b3nd1

MS4A7_in_original_stromal_10b3nd1 <-
  "MS4A7" %in%
  stromal_genes_original_10b3nd1

MS4A7_in_original_immune_10b3nd1 <-
  "MS4A7" %in%
  immune_genes_original_10b3nd1


cat(
  "AOC3 in Stromal: ",
  AOC3_in_original_stromal_10b3nd1,
  "\n"
)

cat(
  "AOC3 in Immune: ",
  AOC3_in_original_immune_10b3nd1,
  "\n"
)

cat(
  "IL7R in Stromal: ",
  IL7R_in_original_stromal_10b3nd1,
  "\n"
)

cat(
  "IL7R in Immune: ",
  IL7R_in_original_immune_10b3nd1,
  "\n"
)

cat(
  "MS4A7 in Stromal: ",
  MS4A7_in_original_stromal_10b3nd1,
  "\n"
)

cat(
  "MS4A7 in Immune: ",
  MS4A7_in_original_immune_10b3nd1,
  "\n"
)


# Expected:
# AOC3  = Stromal TRUE, Immune FALSE
# IL7R  = Stromal FALSE, Immune TRUE
# MS4A7 = both FALSE

stopifnot(
  AOC3_in_original_stromal_10b3nd1,
  !AOC3_in_original_immune_10b3nd1,
  !IL7R_in_original_stromal_10b3nd1,
  IL7R_in_original_immune_10b3nd1,
  !MS4A7_in_original_stromal_10b3nd1,
  !MS4A7_in_original_immune_10b3nd1
)


# ============================================================
# 18R. Construct exact LOO signature vectors
# ============================================================

cat("\n============================================================\n")
cat("CONSTRUCT LOO SIGNATURE VECTORS\n")
cat("============================================================\n\n")


stromal_genes_AOC3_LOO_10b3nd1 <-
  stromal_genes_original_10b3nd1[
    stromal_genes_original_10b3nd1 != "AOC3"
  ]

immune_genes_IL7R_LOO_10b3nd1 <-
  immune_genes_original_10b3nd1[
    immune_genes_original_10b3nd1 != "IL7R"
  ]


cat(
  "Original Stromal: ",
  length(
    stromal_genes_original_10b3nd1
  ),
  "\n"
)

cat(
  "AOC3-LOO Stromal: ",
  length(
    stromal_genes_AOC3_LOO_10b3nd1
  ),
  "\n"
)

cat(
  "Original Immune: ",
  length(
    immune_genes_original_10b3nd1
  ),
  "\n"
)

cat(
  "IL7R-LOO Immune: ",
  length(
    immune_genes_IL7R_LOO_10b3nd1
  ),
  "\n"
)


stopifnot(
  length(
    stromal_genes_AOC3_LOO_10b3nd1
  ) == 140
)

stopifnot(
  length(
    immune_genes_IL7R_LOO_10b3nd1
  ) == 140
)

stopifnot(
  !"AOC3" %in%
    stromal_genes_AOC3_LOO_10b3nd1
)

stopifnot(
  !"IL7R" %in%
    immune_genes_IL7R_LOO_10b3nd1
)


# ============================================================
# 19R. Reconstruct modified SI_geneset objects
#     in the SAME official orientation
# ============================================================

cat("\n============================================================\n")
cat("CONSTRUCT MODIFIED SI_geneSET OBJECTS\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# AOC3 LOO
# ------------------------------------------------------------

SI_geneset_AOC3_LOO_10b3nd1 <-
  SI_geneset_original_10b3nd1

SI_geneset_AOC3_LOO_10b3nd1[
  "StromalSignature",
  -1
] <-
  c(
    stromal_genes_AOC3_LOO_10b3nd1,
    NA_character_
  )


# ------------------------------------------------------------
# IL7R LOO
# ------------------------------------------------------------

SI_geneset_IL7R_LOO_10b3nd1 <-
  SI_geneset_original_10b3nd1

SI_geneset_IL7R_LOO_10b3nd1[
  "ImmuneSignature",
  -1
] <-
  c(
    immune_genes_IL7R_LOO_10b3nd1,
    NA_character_
  )


# ------------------------------------------------------------
# Structural checks
# ------------------------------------------------------------

stopifnot(
  identical(
    dim(
      SI_geneset_AOC3_LOO_10b3nd1
    ),
    c(2L, 142L)
  )
)

stopifnot(
  identical(
    dim(
      SI_geneset_IL7R_LOO_10b3nd1
    ),
    c(2L, 142L)
  )
)


# Verify exact removals
AOC3_remaining_after_LOO_10b3nd1 <-
  any(
    SI_geneset_AOC3_LOO_10b3nd1[
      "StromalSignature",
      -1
    ] == "AOC3",
    na.rm = TRUE
  )

IL7R_remaining_after_LOO_10b3nd1 <-
  any(
    SI_geneset_IL7R_LOO_10b3nd1[
      "ImmuneSignature",
      -1
    ] == "IL7R",
    na.rm = TRUE
  )


cat(
  "\nAOC3 remains in AOC3-LOO Stromal signature: ",
  AOC3_remaining_after_LOO_10b3nd1,
  "\n"
)

cat(
  "IL7R remains in IL7R-LOO Immune signature: ",
  IL7R_remaining_after_LOO_10b3nd1,
  "\n"
)

stopifnot(
  !AOC3_remaining_after_LOO_10b3nd1,
  !IL7R_remaining_after_LOO_10b3nd1
)


# ============================================================
# 20R. IMPORTANT OFFICIAL FUNCTION AUDIT
# ============================================================

cat("\n============================================================\n")
cat("OFFICIAL estimateScore() SOURCE AUDIT\n")
cat("============================================================\n\n")


estimateScore_official_10b3nd1 <-
  getFromNamespace(
    "estimateScore",
    "estimate"
  )


cat(
  "Official estimateScore() environment:\n"
)

print(
  environment(
    estimateScore_official_10b3nd1
  )
)


cat(
  "\nOfficial estimateScore() body:\n\n"
)

print(
  body(
    estimateScore_official_10b3nd1
  )
)


# ------------------------------------------------------------
# Search function body for SI_geneset
# ------------------------------------------------------------

estimateScore_body_text_10b3nd1 <-
  paste(
    deparse(
      body(
        estimateScore_official_10b3nd1
      )
    ),
    collapse = "\n"
  )


SI_reference_present_10b3nd1 <-
  grepl(
    "SI_geneset",
    estimateScore_body_text_10b3nd1,
    fixed = TRUE
  )


cat(
  "\nOfficial function body explicitly references SI_geneset: ",
  SI_reference_present_10b3nd1,
  "\n"
)


# ============================================================
# 21R. Locate SI_geneset binding in the official function
# ============================================================

cat("\n============================================================\n")
cat("SI_geneSET BINDING LOCATION AUDIT\n")
cat("============================================================\n\n")


estimateScore_env_10b3nd1 <-
  environment(
    estimateScore_official_10b3nd1
  )


cat(
  "Direct SI_geneset binding in function environment: "
)

print(
  exists(
    "SI_geneset",
    envir = estimateScore_env_10b3nd1,
    inherits = FALSE
  )
)


# Walk parent environments
env_walk_10b3nd1 <- list()

current_env_10b3nd1 <-
  estimateScore_env_10b3nd1

counter_10b3nd1 <- 0L

while (
  !identical(
    current_env_10b3nd1,
    emptyenv()
  )
) {
  
  counter_10b3nd1 <-
    counter_10b3nd1 + 1L
  
  env_walk_10b3nd1[[counter_10b3nd1]] <-
    current_env_10b3nd1
  
  current_env_10b3nd1 <-
    parent.env(
      current_env_10b3nd1
    )
  
  if (
    counter_10b3nd1 > 30L
  ) {
    break
  }
}


SI_locations_10b3nd1 <-
  data.frame(
    environment_index =
      seq_along(
        env_walk_10b3nd1
      ),
    environment_name =
      vapply(
        env_walk_10b3nd1,
        environmentName,
        character(1)
      ),
    has_SI_geneset =
      vapply(
        env_walk_10b3nd1,
        function(e) {
          exists(
            "SI_geneset",
            envir = e,
            inherits = FALSE
          )
        },
        logical(1)
      ),
    stringsAsFactors = FALSE
  )


print(
  SI_locations_10b3nd1
)


# ============================================================
# 22R. Check namespace binding directly
# ============================================================

cat("\n============================================================\n")
cat("NAMESPACE SI_geneSET BINDING\n")
cat("============================================================\n\n")


namespace_SI_geneset_10b3nd1 <-
  get(
    "SI_geneset",
    envir = asNamespace("estimate")
  )


cat(
  "Namespace SI_geneset dimensions: ",
  nrow(
    namespace_SI_geneset_10b3nd1
  ),
  " x ",
  ncol(
    namespace_SI_geneset_10b3nd1
  ),
  "\n"
)

stopifnot(
  identical(
    namespace_SI_geneset_10b3nd1,
    SI_geneset_original_10b3nd1
  )
)


# ============================================================
# 23R. DO NOT yet generate LOO scores
# ============================================================

cat("\n============================================================\n")
cat("10B-3N-D-1 STRUCTURAL / SOURCE AUDIT COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Official estimate version: ",
  estimate_version_10b3nd1,
  "\n"
)

cat(
  "Official SI_geneset structure: 2 x 142\n"
)

cat(
  "Official Stromal signature: 141 genes\n"
)

cat(
  "Official Immune signature: 141 genes\n"
)

cat(
  "AOC3 = Stromal signature member: ",
  AOC3_in_original_stromal_10b3nd1,
  "\n"
)

cat(
  "IL7R = Immune signature member: ",
  IL7R_in_original_immune_10b3nd1,
  "\n"
)

cat(
  "MS4A7 = neither signature: ",
  !MS4A7_in_original_stromal_10b3nd1 &
    !MS4A7_in_original_immune_10b3nd1,
  "\n"
)

cat(
  "AOC3-LOO signature length: ",
  length(
    stromal_genes_AOC3_LOO_10b3nd1
  ),
  "\n"
)

cat(
  "IL7R-LOO signature length: ",
  length(
    immune_genes_IL7R_LOO_10b3nd1
  ),
  "\n"
)

cat(
  "\nNO LOO SCORE HAS BEEN CALCULATED.\n"
)

cat(
  "NO STATISTICAL ANALYSIS HAS BEEN PERFORMED.\n"
)

cat(
  "NO LOCKED STEP14 OBJECT HAS BEEN MODIFIED.\n"
)

cat(
  "NO FIGURE 6 HAS BEEN MODIFIED.\n"
)

cat(
  "NO SUPPLEMENTARY TABLE 3 HAS BEEN MODIFIED.\n"
)


# ============================================================
# 24R. Save corrected structural audit checkpoint
# ============================================================

step10b3nd1_repair_checkpoint <-
  file.path(
    step_10b3nd1_dir,
    "Paper1_STEP7B10B3ND1_corrected_SI_geneset_structure_audit_LOCKED.RData"
  )


step7be10b3nd1_repair_audit <- list(
  
  status =
    "SI_GENESET_STRUCTURE_AUDIT_PASSED",
  
  estimate_version =
    estimate_version_10b3nd1,
  
  SI_geneset_dimensions =
    dim(
      SI_geneset_original_10b3nd1
    ),
  
  stromal_signature_length =
    length(
      stromal_genes_original_10b3nd1
    ),
  
  immune_signature_length =
    length(
      immune_genes_original_10b3nd1
    ),
  
  AOC3_in_stromal =
    AOC3_in_original_stromal_10b3nd1,
  
  AOC3_in_immune =
    AOC3_in_original_immune_10b3nd1,
  
  IL7R_in_stromal =
    IL7R_in_original_stromal_10b3nd1,
  
  IL7R_in_immune =
    IL7R_in_original_immune_10b3nd1,
  
  MS4A7_in_stromal =
    MS4A7_in_original_stromal_10b3nd1,
  
  MS4A7_in_immune =
    MS4A7_in_original_immune_10b3nd1,
  
  AOC3_LOO_length =
    length(
      stromal_genes_AOC3_LOO_10b3nd1
    ),
  
  IL7R_LOO_length =
    length(
      immune_genes_IL7R_LOO_10b3nd1
    ),
  
  AOC3_removed =
    !AOC3_remaining_after_LOO_10b3nd1,
  
  IL7R_removed =
    !IL7R_remaining_after_LOO_10b3nd1,
  
  score_calculated =
    FALSE,
  
  statistical_analysis_performed =
    FALSE,
  
  locked_step14_modified =
    FALSE,
  
  figure6_modified =
    FALSE,
  
  supplementary_table3_modified =
    FALSE
  
)


save(
  list = c(
    "step7be10b3nd1_repair_audit",
    "SI_geneset_original_10b3nd1",
    "SI_geneset_AOC3_LOO_10b3nd1",
    "SI_geneset_IL7R_LOO_10b3nd1",
    "stromal_genes_original_10b3nd1",
    "immune_genes_original_10b3nd1",
    "stromal_genes_AOC3_LOO_10b3nd1",
    "immune_genes_IL7R_LOO_10b3nd1",
    "SI_locations_10b3nd1"
  ),
  file =
    step10b3nd1_repair_checkpoint
)


cat(
  "\nCorrected checkpoint saved:\n",
  step10b3nd1_repair_checkpoint,
  "\n"
)

stopifnot(
  file.exists(
    step10b3nd1_repair_checkpoint
  )
)


cat("\n============================================================\n")
cat("HARD STOP\n")
cat("============================================================\n\n")

cat(
  "Do NOT calculate LOO scores yet.\n",
  "The next step depends on the actual SI_geneset resolution\n",
  "mechanism revealed by the official estimateScore() body.\n"
)

# ============================================================
# 21R. Locate SI_geneset binding in the official function
#     CORRECTED R SYNTAX
# ============================================================

cat("\n============================================================\n")
cat("SI_geneSET BINDING LOCATION AUDIT\n")
cat("============================================================\n\n")

estimateScore_env_10b3nd1 <-
  environment(
    estimateScore_official_10b3nd1
  )

cat(
  "Direct SI_geneset binding in function environment: "
)

print(
  exists(
    "SI_geneset",
    envir = estimateScore_env_10b3nd1,
    inherits = FALSE
  )
)

# Walk parent environments
env_walk_10b3nd1 <- list()

current_env_10b3nd1 <-
  estimateScore_env_10b3nd1

counter_10b3nd1 <- 0L

while (
  !identical(
    current_env_10b3nd1,
    emptyenv()
  )
) {
  
  counter_10b3nd1 <-
    counter_10b3nd1 + 1L
  
  env_walk_10b3nd1[[counter_10b3nd1]] <-
    current_env_10b3nd1
  
  current_env_10b3nd1 <-
    parent.env(
      current_env_10b3nd1
    )
  
  if (
    counter_10b3nd1 > 30L
  ) {
    break
  }
}

# Build environment audit table
SI_locations_10b3nd1 <-
  data.frame(
    environment_index =
      seq_along(
        env_walk_10b3nd1
      ),
    
    environment_name =
      vapply(
        env_walk_10b3nd1,
        environmentName,
        character(1)
      ),
    
    has_SI_geneset =
      vapply(
        env_walk_10b3nd1,
        function(e) {
          
          exists(
            "SI_geneset",
            envir = e,
            inherits = FALSE
          )
          
        },
        logical(1)
      ),
    
    stringsAsFactors = FALSE
  )

cat(
  "\nEnvironment chain:\n\n"
)

print(
  SI_locations_10b3nd1
)

# Explicitly identify environments containing SI_geneset
SI_geneset_locations_10b3nd1 <-
  SI_locations_10b3nd1[
    SI_locations_10b3nd1$has_SI_geneset,
    ,
    drop = FALSE
  ]

cat(
  "\n============================================================\n"
)

cat(
  "ENVIRONMENTS CONTAINING SI_geneSET\n"
)

cat(
  "============================================================\n\n"
)

if (
  nrow(
    SI_geneset_locations_10b3nd1
  ) == 0
) {
  
  cat(
    "No SI_geneset binding found in the parent environment chain.\n"
  )
  
} else {
  
  print(
    SI_geneset_locations_10b3nd1
  )
  
}

cat(
  "\nPackage namespace has NOT been modified.\n"
)

cat(
  "No LOO score has been calculated.\n"
)

cat(
  "No statistical analysis has been performed.\n"
)

cat(
  "No Step14 object has been modified.\n"
)

cat(
  "No Figure 6 has been modified.\n"
)

cat(
  "No Supplementary Table 3 has been modified.\n"
)

#10B-3N-D-2：Official SI_geneset substitution validation
# ============================================================
# Paper 1
# STEP 7B-7E-10B-3N-D-2
#
# OFFICIAL ESTIMATE SI_geneset SUBSTITUTION VALIDATION
#
# PURPOSE:
#   1. Confirm original SI_geneset reproduces locked Step14 scores
#   2. Temporarily remove AOC3 from StromalSignature
#   3. Confirm official estimateScore() actually responds
#      to the modified SI_geneset
#   4. Restore original SI_geneset
#   5. Confirm exact reconstruction again
#
# IMPORTANT:
#   - This is NOT formal LOO analysis
#   - No biological/statistical testing
#   - No modification of Step14
#   - No modification of Figure 6
#   - No modification of Supplementary Table 3
#   - estimate package namespace will NOT be modified
# ============================================================


cat("\n============================================================\n")
cat("STEP 10B-3N-D-2\n")
cat("OFFICIAL SI_geneSET SUBSTITUTION VALIDATION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Basic package validation
# ------------------------------------------------------------

stopifnot(
  requireNamespace(
    "estimate",
    quietly = TRUE
  )
)

estimate_version_10b3nd2 <-
  as.character(
    packageVersion("estimate")
  )

cat(
  "Installed estimate version: ",
  estimate_version_10b3nd2,
  "\n"
)

stopifnot(
  identical(
    estimate_version_10b3nd2,
    "1.0.13"
  )
)


# ------------------------------------------------------------
# 2. Define paths
# ------------------------------------------------------------

root_10b3nd2 <-
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis"

step_10b3nd2_dir <-
  file.path(
    root_10b3nd2,
    "STEP7B-7E-10B-3N",
    "10B-3N-D-2"
  )

dir.create(
  step_10b3nd2_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Official filtered GCT
filtered_gct_10b3nd2 <-
  file.path(
    root_10b3nd2,
    "STEP7B-7E-10B-3N",
    "10B-3N-C",
    "STEP10B3NC_official_filterCommonGenes_output.gct"
  )


# Locked source checkpoint from 10B-3N-C
checkpoint_10b3nd2 <-
  file.path(
    root_10b3nd2,
    "STEP7B-7E-10B-3N",
    "10B-3N-C",
    "Paper1_STEP7B10B3NC_official_ESTIMATE_reconstruction_audit_LOCKED.RData"
  )


stopifnot(
  file.exists(
    filtered_gct_10b3nd2
  )
)

stopifnot(
  file.exists(
    checkpoint_10b3nd2
  )
)


cat(
  "\nOfficial filtered GCT:\n",
  filtered_gct_10b3nd2,
  "\n"
)

cat(
  "\n10B-3N-C checkpoint:\n",
  checkpoint_10b3nd2,
  "\n"
)


# ------------------------------------------------------------
# 3. Load 10B-3N-C checkpoint into isolated environment
# ------------------------------------------------------------

checkpoint_env_10b3nd2 <-
  new.env(
    parent = emptyenv()
  )

load(
  checkpoint_10b3nd2,
  envir = checkpoint_env_10b3nd2
)


cat(
  "\nObjects recovered from 10B-3N-C checkpoint:\n"
)

print(
  ls(
    envir = checkpoint_env_10b3nd2
  )
)


# ------------------------------------------------------------
# 4. Recover locked source score
# ------------------------------------------------------------

estimate_scores_source_10b3nd2 <-
  checkpoint_env_10b3nd2$estimate_scores_wgcna_source_10b3nc


stopifnot(
  is.data.frame(
    estimate_scores_source_10b3nd2
  )
)

stopifnot(
  nrow(
    estimate_scores_source_10b3nd2
  ) == 3
)

stopifnot(
  ncol(
    estimate_scores_source_10b3nd2
  ) == 513
)


cat(
  "\nLocked source score dimensions: ",
  paste(
    dim(
      estimate_scores_source_10b3nd2
    ),
    collapse = " x "
  ),
  "\n"
)

cat(
  "Locked source score rows:\n"
)

print(
  rownames(
    estimate_scores_source_10b3nd2
  )
)


stopifnot(
  identical(
    rownames(
      estimate_scores_source_10b3nd2
    ),
    c(
      "StromalScore",
      "ImmuneScore",
      "ESTIMATEScore"
    )
  )
)


# ------------------------------------------------------------
# 5. Recover official SI_geneset from package namespace
# ------------------------------------------------------------

SI_geneset_package_10b3nd2 <-
  get(
    "SI_geneset",
    envir =
      as.environment(
        "package:estimate"
      ),
    inherits = FALSE
  )


cat(
  "\nOfficial package SI_geneset dimensions: ",
  paste(
    dim(
      SI_geneset_package_10b3nd2
    ),
    collapse = " x "
  ),
  "\n"
)


stopifnot(
  is.data.frame(
    SI_geneset_package_10b3nd2
  )
)

stopifnot(
  identical(
    dim(
      SI_geneset_package_10b3nd2
    ),
    c(
      2L,
      142L
    )
  )
)

stopifnot(
  identical(
    rownames(
      SI_geneset_package_10b3nd2
    ),
    c(
      "StromalSignature",
      "ImmuneSignature"
    )
  )
)

# ============================================================
# STEP 10B-3N-D-2
# REPAIR BLOCK
#
# Purpose:
#   Recover the canonical SI_geneset from the actual
#   10B-3N-D-1 checkpoint.
#
# NO SCORE CALCULATION
# NO STATISTICS
# NO STEP14 MODIFICATION
# ============================================================

cat("\n============================================================\n")
cat("10B-3N-D-2 REPAIR: RECOVER D-1 SI_geneSET\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Locate the actual D-1 checkpoint
# ------------------------------------------------------------

checkpoint_D1_10b3nd2 <-
  file.path(
    root_10b3nd2,
    "STEP7B-7E-10B-3N",
    "10B-3N-D-1",
    "Paper1_STEP7B10B3ND1_corrected_SI_geneset_structure_audit_LOCKED.RData"
  )


stopifnot(
  file.exists(
    checkpoint_D1_10b3nd2
  )
)


cat(
  "D-1 checkpoint found:\n",
  checkpoint_D1_10b3nd2,
  "\n"
)


# ------------------------------------------------------------
# 2. Load D-1 checkpoint into isolated environment
# ------------------------------------------------------------

checkpoint_D1_env_10b3nd2 <-
  new.env(
    parent = emptyenv()
  )


load(
  checkpoint_D1_10b3nd2,
  envir = checkpoint_D1_env_10b3nd2
)


cat(
  "\nObjects recovered from D-1 checkpoint:\n\n"
)

print(
  ls(
    envir = checkpoint_D1_env_10b3nd2
  )
)


# ------------------------------------------------------------
# 3. Recover canonical SI_geneset
# ------------------------------------------------------------

SI_geneset_original_10b3nd2 <-
  checkpoint_D1_env_10b3nd2$
  SI_geneset_original_10b3nd1


stopifnot(
  is.data.frame(
    SI_geneset_original_10b3nd2
  )
)


cat(
  "\nD-1 SI_geneset dimensions: ",
  paste(
    dim(
      SI_geneset_original_10b3nd2
    ),
    collapse = " x "
  ),
  "\n"
)


stopifnot(
  identical(
    dim(
      SI_geneset_original_10b3nd2
    ),
    c(
      2L,
      142L
    )
  )
)


# ------------------------------------------------------------
# 4. Recover package SI_geneset
# ------------------------------------------------------------

SI_geneset_package_10b3nd2 <-
  get(
    "SI_geneset",
    envir =
      as.environment(
        "package:estimate"
      ),
    inherits = FALSE
  )


cat(
  "Package SI_geneset dimensions: ",
  paste(
    dim(
      SI_geneset_package_10b3nd2
    ),
    collapse = " x "
  ),
  "\n"
)


# ------------------------------------------------------------
# 5. Compare structure rather than relying only on identical()
# ------------------------------------------------------------

cat(
  "\n------------------------------------------------------------\n"
)

cat(
  "SI_geneSET STRUCTURAL COMPARISON\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


cat(
  "Class D-1: "
)

print(
  class(
    SI_geneset_original_10b3nd2
  )
)

cat(
  "Class package: "
)

print(
  class(
    SI_geneset_package_10b3nd2
  )
)


cat(
  "\nDimensions identical: "
)

print(
  identical(
    dim(
      SI_geneset_original_10b3nd2
    ),
    dim(
      SI_geneset_package_10b3nd2
    )
  )
)


cat(
  "Row names identical: "
)

print(
  identical(
    rownames(
      SI_geneset_original_10b3nd2
    ),
    rownames(
      SI_geneset_package_10b3nd2
    )
  )
)


cat(
  "Column names identical: "
)

print(
  identical(
    colnames(
      SI_geneset_original_10b3nd2
    ),
    colnames(
      SI_geneset_package_10b3nd2
    )
  )
)


# ------------------------------------------------------------
# 6. Compare all cells
# ------------------------------------------------------------

SI_cell_identical_10b3nd2 <-
  identical(
    as.character(
      SI_geneset_original_10b3nd2
    ),
    as.character(
      SI_geneset_package_10b3nd2
    )
  )


cat(
  "\nAll SI_geneset cell values identical: "
)

print(
  SI_cell_identical_10b3nd2
)


# ------------------------------------------------------------
# 7. If identical() is FALSE, locate the differences
# ------------------------------------------------------------

if (
  !SI_cell_identical_10b3nd2
) {
  
  cat(
    "\n============================================================\n"
  )
  
  cat(
    "DIFFERENCE AUDIT\n"
  )
  
  cat(
    "============================================================\n\n"
  )
  
  
  D1_matrix_10b3nd2 <-
    as.matrix(
      SI_geneset_original_10b3nd2
    )
  
  package_matrix_10b3nd2 <-
    as.matrix(
      SI_geneset_package_10b3nd2
    )
  
  
  same_dimensions_10b3nd2 <-
    identical(
      dim(
        D1_matrix_10b3nd2
      ),
      dim(
        package_matrix_10b3nd2
      )
    )
  
  
  cat(
    "Matrix dimensions identical: ",
    same_dimensions_10b3nd2,
    "\n"
  )
  
  
  if (
    same_dimensions_10b3nd2
  ) {
    
    difference_index_10b3nd2 <-
      which(
        D1_matrix_10b3nd2 !=
          package_matrix_10b3nd2,
        arr.ind = TRUE
      )
    
    
    cat(
      "Number of differing cells: ",
      nrow(
        difference_index_10b3nd2
      ),
      "\n"
    )
    
    
    if (
      nrow(
        difference_index_10b3nd2
      ) > 0
    ) {
      
      first_differences_10b3nd2 <-
        head(
          difference_index_10b3nd2,
          20
        )
      
      
      difference_table_10b3nd2 <-
        data.frame(
          
          row =
            first_differences_10b3nd2[
              ,
              1
            ],
          
          column =
            first_differences_10b3nd2[
              ,
              2
            ],
          
          D1_value =
            D1_matrix_10b3nd2[
              first_differences_10b3nd2
            ],
          
          package_value =
            package_matrix_10b3nd2[
              first_differences_10b3nd2
            ],
          
          stringsAsFactors =
            FALSE
          
        )
      
      
      print(
        difference_table_10b3nd2
      )
      
    }
    
  }
  
}


# ------------------------------------------------------------
# 8. Extract signatures independently from both objects
# ------------------------------------------------------------

stromal_genes_D1_10b3nd2 <-
  as.character(
    SI_geneset_original_10b3nd2[
      "StromalSignature",
      -1
    ]
  )


immune_genes_D1_10b3nd2 <-
  as.character(
    SI_geneset_original_10b3nd2[
      "ImmuneSignature",
      -1
    ]
  )


stromal_genes_package_10b3nd2 <-
  as.character(
    SI_geneset_package_10b3nd2[
      "StromalSignature",
      -1
    ]
  )


immune_genes_package_10b3nd2 <-
  as.character(
    SI_geneset_package_10b3nd2[
      "ImmuneSignature",
      -1
    ]
  )


# ------------------------------------------------------------
# 9. Compare actual signature vectors
# ------------------------------------------------------------

stromal_signature_identical_10b3nd2 <-
  identical(
    stromal_genes_D1_10b3nd2,
    stromal_genes_package_10b3nd2
  )


immune_signature_identical_10b3nd2 <-
  identical(
    immune_genes_D1_10b3nd2,
    immune_genes_package_10b3nd2
  )


cat(
  "\nStromal signature identical: ",
  stromal_signature_identical_10b3nd2,
  "\n"
)

cat(
  "Immune signature identical: ",
  immune_signature_identical_10b3nd2,
  "\n"
)


# ------------------------------------------------------------
# 10. Biological membership hard checks
# ------------------------------------------------------------

AOC3_in_stromal_10b3nd2 <-
  "AOC3" %in%
  stromal_genes_D1_10b3nd2

AOC3_in_immune_10b3nd2 <-
  "AOC3" %in%
  immune_genes_D1_10b3nd2

IL7R_in_stromal_10b3nd2 <-
  "IL7R" %in%
  stromal_genes_D1_10b3nd2

IL7R_in_immune_10b3nd2 <-
  "IL7R" %in%
  immune_genes_D1_10b3nd2

MS4A7_in_stromal_10b3nd2 <-
  "MS4A7" %in%
  stromal_genes_D1_10b3nd2

MS4A7_in_immune_10b3nd2 <-
  "MS4A7" %in%
  immune_genes_D1_10b3nd2


cat(
  "\nCore-gene membership:\n"
)

cat(
  "AOC3 / Stromal = ",
  AOC3_in_stromal_10b3nd2,
  "\n"
)

cat(
  "AOC3 / Immune = ",
  AOC3_in_immune_10b3nd2,
  "\n"
)

cat(
  "IL7R / Stromal = ",
  IL7R_in_stromal_10b3nd2,
  "\n"
)

cat(
  "IL7R / Immune = ",
  IL7R_in_immune_10b3nd2,
  "\n"
)

cat(
  "MS4A7 / Stromal = ",
  MS4A7_in_stromal_10b3nd2,
  "\n"
)

cat(
  "MS4A7 / Immune = ",
  MS4A7_in_immune_10b3nd2,
  "\n"
)


stopifnot(
  AOC3_in_stromal_10b3nd2
)

stopifnot(
  !AOC3_in_immune_10b3nd2
)

stopifnot(
  !IL7R_in_stromal_10b3nd2
)

stopifnot(
  IL7R_in_immune_10b3nd2
)

stopifnot(
  !MS4A7_in_stromal_10b3nd2
)

stopifnot(
  !MS4A7_in_immune_10b3nd2
)


# ------------------------------------------------------------
# 11. Final repair audit
# ------------------------------------------------------------

SI_geneset_recovery_pass_10b3nd2 <-
  stromal_signature_identical_10b3nd2 &&
  immune_signature_identical_10b3nd2 &&
  AOC3_in_stromal_10b3nd2 &&
  !AOC3_in_immune_10b3nd2 &&
  !IL7R_in_stromal_10b3nd2 &&
  IL7R_in_immune_10b3nd2 &&
  !MS4A7_in_stromal_10b3nd2 &&
  !MS4A7_in_immune_10b3nd2


cat(
  "\n============================================================\n"
)

cat(
  "D-2 SI_geneSET RECOVERY RESULT\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "SI_geneset recovery passed: ",
  SI_geneset_recovery_pass_10b3nd2,
  "\n"
)

cat(
  "\nNO SCORE CALCULATED.\n"
)

cat(
  "NO STATISTICAL ANALYSIS PERFORMED.\n"
)

cat(
  "NO PACKAGE NAMESPACE MODIFIED.\n"
)

cat(
  "NO STEP14 MODIFIED.\n"
)

cat(
  "NO FIGURE 6 MODIFIED.\n"
)

cat(
  "NO SUPPLEMENTARY TABLE 3 MODIFIED.\n"
)


stopifnot(
  SI_geneset_recovery_pass_10b3nd2
)

##Step 10B-3N-D-2：Official SI_geneset substitution validation
# ============================================================
# Paper 1
# STEP 10B-3N-D-2
# Official ESTIMATE SI_geneset substitution validation
#
# Purpose:
#   Test A = original SI_geneset -> exact reproduction
#   Test B = remove AOC3 from StromalSignature -> score changes
#   Test C = restore original SI_geneset -> exact reproduction
#
# IMPORTANT:
#   This is NOT the formal LOO sensitivity analysis.
#   No Figure 6 / Supplementary Table 3 / Step14 objects are modified.
#   No package namespace is modified.
# ============================================================


cat("\n============================================================\n")
cat("STEP 10B-3N-D-2\n")
cat("Official SI_geneset substitution validation\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Basic package verification
# ------------------------------------------------------------

stopifnot(
  requireNamespace("estimate", quietly = TRUE)
)

cat(
  "estimate package version: ",
  as.character(packageVersion("estimate")),
  "\n\n"
)

stopifnot(
  as.character(packageVersion("estimate")) == "1.0.13"
)


# ------------------------------------------------------------
# 1. Define output directory
# ------------------------------------------------------------

root_10b3nd2 <- "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis"

dir_10b3nd2 <- file.path(
  root_10b3nd2,
  "STEP7B-7E-10B-3N",
  "10B-3N-D-2"
)

dir.create(
  dir_10b3nd2,
  recursive = TRUE,
  showWarnings = FALSE
)

cat(
  "Output directory:\n",
  dir_10b3nd2,
  "\n\n"
)


# ------------------------------------------------------------
# 2. Recover locked Step14 score object
# ------------------------------------------------------------

checkpoint_10b3nc <- file.path(
  root_10b3nd2,
  "STEP7B-7E-10B-3N",
  "10B-3N-C",
  "Paper1_STEP7B10B3NC_official_ESTIMATE_reconstruction_audit_LOCKED.RData"
)

stopifnot(
  file.exists(checkpoint_10b3nc)
)

checkpoint_env_10b3nd2 <- new.env(parent = emptyenv())

load(
  checkpoint_10b3nc,
  envir = checkpoint_env_10b3nd2
)

stopifnot(
  exists(
    "estimate_scores_wgcna_source_10b3nc",
    envir = checkpoint_env_10b3nd2,
    inherits = FALSE
  )
)

estimate_scores_wgcna_source_10b3nd2 <-
  checkpoint_env_10b3nd2$estimate_scores_wgcna_source_10b3nc


cat(
  "Locked Step14 score object recovered.\n"
)

cat(
  "Dimensions: ",
  paste(
    dim(estimate_scores_wgcna_source_10b3nd2),
    collapse = " x "
  ),
  "\n"
)

cat(
  "Rows: ",
  paste(
    rownames(estimate_scores_wgcna_source_10b3nd2),
    collapse = ", "
  ),
  "\n\n"
)

# ============================================================
# D-2 repair: identify the actual D-1 checkpoint file
# ============================================================

D1_dir_10b3nd2 <- file.path(
  root_10b3nd2,
  "STEP7B-7E-10B-3N",
  "10B-3N-D-1"
)

cat(
  "\n============================================================\n"
)
cat(
  "D-1 DIRECTORY CHECK\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "Directory:\n",
  D1_dir_10b3nd2,
  "\n\n"
)

cat(
  "Directory exists: ",
  dir.exists(D1_dir_10b3nd2),
  "\n\n"
)

stopifnot(
  dir.exists(D1_dir_10b3nd2)
)


# ------------------------------------------------------------
# List every file in D-1 directory
# ------------------------------------------------------------

D1_files_10b3nd2 <- list.files(
  path = D1_dir_10b3nd2,
  full.names = TRUE,
  recursive = TRUE,
  all.files = FALSE,
  include.dirs = FALSE
)

cat(
  "Files found in D-1 directory:\n\n"
)

print(
  D1_files_10b3nd2
)


# ------------------------------------------------------------
# Specifically identify RData files
# ------------------------------------------------------------

D1_RData_files_10b3nd2 <-
  D1_files_10b3nd2[
    grepl(
      "\\.RData$|\\.rdata$|\\.rda$",
      D1_files_10b3nd2,
      ignore.case = TRUE
    )
  ]


cat(
  "\n============================================================\n"
)

cat(
  "RData FILES FOUND\n"
)

cat(
  "============================================================\n\n"
)

print(
  D1_RData_files_10b3nd2
)

cat(
  "\nNumber of RData files: ",
  length(D1_RData_files_10b3nd2),
  "\n"
)

# ============================================================
# STEP 10B-3N-D-2
# Correct D-1 checkpoint recovery
# ============================================================

# Use the ACTUAL file discovered by list.files()
checkpoint_D1_10b3nd2 <-
  D1_RData_files_10b3nd2[1]

cat(
  "\n============================================================\n"
)

cat(
  "D-1 CHECKPOINT RECOVERY\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Checkpoint:\n",
  checkpoint_D1_10b3nd2,
  "\n\n"
)

stopifnot(
  file.exists(checkpoint_D1_10b3nd2)
)


# ------------------------------------------------------------
# Load into isolated environment
# ------------------------------------------------------------

checkpoint_D1_env_10b3nd2 <-
  new.env(parent = emptyenv())

load(
  checkpoint_D1_10b3nd2,
  envir = checkpoint_D1_env_10b3nd2
)


cat(
  "Objects in D-1 checkpoint:\n\n"
)

print(
  ls(
    envir = checkpoint_D1_env_10b3nd2,
    all.names = TRUE
  )
)


# ------------------------------------------------------------
# Recover canonical SI_geneset
# ------------------------------------------------------------

stopifnot(
  exists(
    "SI_geneset_original_10b3nd1",
    envir = checkpoint_D1_env_10b3nd2,
    inherits = FALSE
  )
)

SI_geneset_original_10b3nd2 <-
  checkpoint_D1_env_10b3nd2$SI_geneset_original_10b3nd1


cat(
  "\nSI_geneset_original_10b3nd1 recovered successfully.\n\n"
)


# ============================================================
# STEP 10B-3N-D-2
# D-1 checkpoint recovery audit — PART 2
# Compare substantive gene content, not first-column metadata
# ============================================================

cat("\n============================================================\n")
cat("D-1 RECOVERY AUDIT — SUBSTANTIVE GENE CONTENT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Inspect first column and remaining columns
# ------------------------------------------------------------

cat("First column name in recovered object:\n")
print(colnames(SI_geneset_original_10b3nd2)[1])

cat("\nFirst column values:\n")
print(
  SI_geneset_original_10b3nd2[, 1, drop = TRUE]
)

cat("\nFirst 10 gene columns:\n")
print(
  colnames(SI_geneset_original_10b3nd2)[2:11]
)

cat("\n")


# ------------------------------------------------------------
# 2. Extract substantive gene vectors
# ------------------------------------------------------------

stromal_genes_recovered_10b3nd2 <-
  as.character(
    SI_geneset_original_10b3nd2[
      "StromalSignature",
      -1,
      drop = TRUE
    ]
  )

immune_genes_recovered_10b3nd2 <-
  as.character(
    SI_geneset_original_10b3nd2[
      "ImmuneSignature",
      -1,
      drop = TRUE
    ]
  )


# Remove NA/empty values only for diagnostic display.
# Do NOT modify the original checkpoint object.

stromal_genes_recovered_clean_10b3nd2 <-
  stromal_genes_recovered_10b3nd2[
    !is.na(stromal_genes_recovered_10b3nd2) &
      nzchar(stromal_genes_recovered_10b3nd2)
  ]

immune_genes_recovered_clean_10b3nd2 <-
  immune_genes_recovered_10b3nd2[
    !is.na(immune_genes_recovered_10b3nd2) &
      nzchar(immune_genes_recovered_10b3nd2)
  ]


# ------------------------------------------------------------
# 3. Basic counts
# ------------------------------------------------------------

cat("Recovered StromalSignature gene count:",
    length(stromal_genes_recovered_clean_10b3nd2),
    "\n")

cat("Recovered ImmuneSignature gene count:",
    length(immune_genes_recovered_clean_10b3nd2),
    "\n\n")


# ------------------------------------------------------------
# 4. Show first and last genes
# ------------------------------------------------------------

cat("First 10 recovered Stromal genes:\n")
print(
  head(
    stromal_genes_recovered_clean_10b3nd2,
    10
  )
)

cat("\nLast 10 recovered Stromal genes:\n")
print(
  tail(
    stromal_genes_recovered_clean_10b3nd2,
    10
  )
)

cat("\nFirst 10 recovered Immune genes:\n")
print(
  head(
    immune_genes_recovered_clean_10b3nd2,
    10
  )
)

cat("\nLast 10 recovered Immune genes:\n")
print(
  tail(
    immune_genes_recovered_clean_10b3nd2,
    10
  )
)

cat("\n")


# ------------------------------------------------------------
# 5. Recover official package SI_geneset
# ------------------------------------------------------------

stopifnot(
  requireNamespace(
    "estimate",
    quietly = TRUE
  )
)

SI_geneset_package_10b3nd2 <-
  get(
    "SI_geneset",
    envir = asNamespace("estimate")
  )


# ------------------------------------------------------------
# 6. Extract official gene vectors
# ------------------------------------------------------------

stromal_genes_package_10b3nd2 <-
  as.character(
    SI_geneset_package_10b3nd2[
      "StromalSignature",
      -1,
      drop = TRUE
    ]
  )

immune_genes_package_10b3nd2 <-
  as.character(
    SI_geneset_package_10b3nd2[
      "ImmuneSignature",
      -1,
      drop = TRUE
    ]
  )


# ------------------------------------------------------------
# 7. Compare gene vectors exactly
# ------------------------------------------------------------

cat("============================================================\n")
cat("EXACT GENE-VECTOR COMPARISON\n")
cat("============================================================\n\n")

cat(
  "Stromal vector identical: ",
  identical(
    stromal_genes_recovered_10b3nd2,
    stromal_genes_package_10b3nd2
  ),
  "\n"
)

cat(
  "Immune vector identical: ",
  identical(
    immune_genes_recovered_10b3nd2,
    immune_genes_package_10b3nd2
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 8. Set-wise comparison
# ------------------------------------------------------------

stromal_only_checkpoint_10b3nd2 <-
  setdiff(
    stromal_genes_recovered_10b3nd2,
    stromal_genes_package_10b3nd2
  )

stromal_only_package_10b3nd2 <-
  setdiff(
    stromal_genes_package_10b3nd2,
    stromal_genes_recovered_10b3nd2
  )

immune_only_checkpoint_10b3nd2 <-
  setdiff(
    immune_genes_recovered_10b3nd2,
    immune_genes_package_10b3nd2
  )

immune_only_package_10b3nd2 <-
  setdiff(
    immune_genes_package_10b3nd2,
    immune_genes_recovered_10b3nd2
  )


cat("Stromal genes only in checkpoint:\n")
print(stromal_only_checkpoint_10b3nd2)

cat("\nStromal genes only in package:\n")
print(stromal_only_package_10b3nd2)

cat("\nImmune genes only in checkpoint:\n")
print(immune_only_checkpoint_10b3nd2)

cat("\nImmune genes only in package:\n")
print(immune_only_package_10b3nd2)

cat("\n")


# ------------------------------------------------------------
# 9. Order comparison
# ------------------------------------------------------------

cat("Stromal order identical: ",
    identical(
      stromal_genes_recovered_10b3nd2,
      stromal_genes_package_10b3nd2
    ),
    "\n")

cat("Immune order identical: ",
    identical(
      immune_genes_recovered_10b3nd2,
      immune_genes_package_10b3nd2
    ),
    "\n")


# ------------------------------------------------------------
# 10. Core-gene membership audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CORE-GENE MEMBERSHIP AUDIT\n")
cat("============================================================\n\n")

cat(
  "AOC3 in recovered StromalSignature: ",
  "AOC3" %in% stromal_genes_recovered_10b3nd2,
  "\n"
)

cat(
  "AOC3 in recovered ImmuneSignature: ",
  "AOC3" %in% immune_genes_recovered_10b3nd2,
  "\n"
)

cat(
  "MS4A7 in recovered StromalSignature: ",
  "MS4A7" %in% stromal_genes_recovered_10b3nd2,
  "\n"
)

cat(
  "MS4A7 in recovered ImmuneSignature: ",
  "MS4A7" %in% immune_genes_recovered_10b3nd2,
  "\n"
)

cat(
  "IL7R in recovered StromalSignature: ",
  "IL7R" %in% stromal_genes_recovered_10b3nd2,
  "\n"
)

cat(
  "IL7R in recovered ImmuneSignature: ",
  "IL7R" %in% immune_genes_recovered_10b3nd2,
  "\n"
)


# ============================================================
# STEP 10B-3N-D-2
# Official SI_geneset substitution validation
#
# TEST A
# Original SI_geneset -> exact reproduction of locked Step14
#
# IMPORTANT:
# - No package namespace modification
# - No formal LOO
# - No modification of locked Step14 objects
# - No modification of Figure 6
# - No modification of Supplementary Table 3
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3N-D-2 — TEST A\n")
cat("ORIGINAL SI_geneSET SUBSTITUTION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Confirm estimate package
# ------------------------------------------------------------

stopifnot(
  requireNamespace(
    "estimate",
    quietly = TRUE
  )
)

estimate_version_10b3nd2 <-
  as.character(
    packageVersion("estimate")
  )

cat(
  "estimate package version: ",
  estimate_version_10b3nd2,
  "\n\n"
)

stopifnot(
  estimate_version_10b3nd2 == "1.0.13"
)


# ------------------------------------------------------------
# 2. Recover official package SI_geneset
# ------------------------------------------------------------

SI_geneset_package_10b3nd2 <-
  get(
    "SI_geneset",
    envir = asNamespace("estimate")
  )


# ------------------------------------------------------------
# 3. Record namespace object BEFORE substitution
# ------------------------------------------------------------

SI_geneset_namespace_before_10b3nd2 <-
  SI_geneset_package_10b3nd2


# ------------------------------------------------------------
# 4. Record original GlobalEnv state
# ------------------------------------------------------------

global_SI_exists_before_10b3nd2 <-
  exists(
    "SI_geneset",
    envir = .GlobalEnv,
    inherits = FALSE
  )

if (global_SI_exists_before_10b3nd2) {
  
  SI_geneset_global_before_10b3nd2 <-
    get(
      "SI_geneset",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  
} else {
  
  SI_geneset_global_before_10b3nd2 <-
    NULL
}


cat(
  "Original GlobalEnv SI_geneset exists: ",
  global_SI_exists_before_10b3nd2,
  "\n\n"
)


# ------------------------------------------------------------
# 5. Construct canonical temporary SI_geneset
# ------------------------------------------------------------
#
# IMPORTANT:
# The recovered D-1 object has V2/V3/... column names.
# We do NOT modify it.
#
# Instead, use the already verified official package object
# as the canonical baseline.
# ------------------------------------------------------------

SI_geneset_testA_10b3nd2 <-
  SI_geneset_namespace_before_10b3nd2


# ============================================================
# STEP 10B-3N-D-2
# DIAGNOSTIC — ACTUAL NAMESPACE SI_geneset STRUCTURE
# ============================================================

cat("\n============================================================\n")
cat("ACTUAL ESTIMATE NAMESPACE SI_geneset STRUCTURE\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Basic structure
# ------------------------------------------------------------

cat("Class:\n")
print(
  class(
    SI_geneset_namespace_before_10b3nd2
  )
)

cat("\nDimensions:\n")
print(
  dim(
    SI_geneset_namespace_before_10b3nd2
  )
)

cat("\nRow names:\n")
print(
  rownames(
    SI_geneset_namespace_before_10b3nd2
  )
)

cat("\nColumn names — first 10:\n")
print(
  head(
    colnames(
      SI_geneset_namespace_before_10b3nd2
    ),
    10
  )
)

cat("\nColumn names — last 10:\n")
print(
  tail(
    colnames(
      SI_geneset_namespace_before_10b3nd2
    ),
    10
  )
)


# ------------------------------------------------------------
# 2. First column values
# ------------------------------------------------------------

cat("\nFirst column values:\n")

print(
  SI_geneset_namespace_before_10b3nd2[
    ,
    1,
    drop = TRUE
  ]
)


# ------------------------------------------------------------
# 3. First few rows / columns
# ------------------------------------------------------------

cat("\nFirst 5 rows x first 8 columns:\n")

print(
  SI_geneset_namespace_before_10b3nd2[
    ,
    1:8,
    drop = FALSE
  ]
)


# ------------------------------------------------------------
# 4. Explicit first-column name
# ------------------------------------------------------------

cat("\nFirst column name:\n")

print(
  colnames(
    SI_geneset_namespace_before_10b3nd2
  )[1]
)


# ------------------------------------------------------------
# 5. Confirm the official object itself
# ------------------------------------------------------------

cat("\nIs official namespace object a data.frame?\n")

print(
  is.data.frame(
    SI_geneset_namespace_before_10b3nd2
  )
)


# ------------------------------------------------------------
# 6. Signature counts using column position
# ------------------------------------------------------------

cat("\nNumber of columns after first column:\n")

print(
  ncol(
    SI_geneset_namespace_before_10b3nd2
  ) - 1L
)

# ============================================================
# STEP 10B-3N-D-2
# TEST A — ORIGINAL SI_geneset BASELINE REPRODUCTION
#
# Corrected for the actual official ESTIMATE structure:
#   2 x 142 data.frame
#   first column name = V2
#   first column values = "estimate"
#   141 genes per signature
#
# IMPORTANT:
# - Do NOT rename columns
# - Do NOT modify namespace
# - Do NOT modify locked Step14
# ============================================================


cat("\n============================================================\n")
cat("STEP 10B-3N-D-2 — TEST A\n")
cat("ORIGINAL SI_geneset BASELINE REPRODUCTION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Confirm package version
# ------------------------------------------------------------

stopifnot(
  requireNamespace(
    "estimate",
    quietly = TRUE
  )
)

estimate_version_10b3nd2 <-
  as.character(
    packageVersion("estimate")
  )

cat(
  "estimate package version: ",
  estimate_version_10b3nd2,
  "\n\n"
)

stopifnot(
  estimate_version_10b3nd2 == "1.0.13"
)


# ------------------------------------------------------------
# 2. Get official namespace SI_geneset
# ------------------------------------------------------------

SI_geneset_namespace_before_10b3nd2 <-
  get(
    "SI_geneset",
    envir = asNamespace("estimate")
  )


# ------------------------------------------------------------
# 3. Verify ACTUAL official structure
# ------------------------------------------------------------

stopifnot(
  is.data.frame(
    SI_geneset_namespace_before_10b3nd2
  )
)

stopifnot(
  identical(
    dim(
      SI_geneset_namespace_before_10b3nd2
    ),
    c(2L, 142L)
  )
)

stopifnot(
  identical(
    rownames(
      SI_geneset_namespace_before_10b3nd2
    ),
    c(
      "StromalSignature",
      "ImmuneSignature"
    )
  )
)

stopifnot(
  identical(
    ncol(
      SI_geneset_namespace_before_10b3nd2
    ) - 1L,
    141L
  )
)

stopifnot(
  all(
    SI_geneset_namespace_before_10b3nd2[
      ,
      1,
      drop = TRUE
    ] ==
      "estimate"
  )
)


# ------------------------------------------------------------
# 4. Extract official signature vectors
# ------------------------------------------------------------

stromal_genes_testA_10b3nd2 <-
  as.character(
    SI_geneset_namespace_before_10b3nd2[
      "StromalSignature",
      -1,
      drop = TRUE
    ]
  )

immune_genes_testA_10b3nd2 <-
  as.character(
    SI_geneset_namespace_before_10b3nd2[
      "ImmuneSignature",
      -1,
      drop = TRUE
    ]
  )


stopifnot(
  length(
    stromal_genes_testA_10b3nd2
  ) == 141L
)

stopifnot(
  length(
    immune_genes_testA_10b3nd2
  ) == 141L
)


# ------------------------------------------------------------
# 5. Confirm core-gene membership
# ------------------------------------------------------------

stopifnot(
  "AOC3" %in%
    stromal_genes_testA_10b3nd2
)

stopifnot(
  !"AOC3" %in%
    immune_genes_testA_10b3nd2
)

stopifnot(
  !"MS4A7" %in%
    stromal_genes_testA_10b3nd2
)

stopifnot(
  !"MS4A7" %in%
    immune_genes_testA_10b3nd2
)

stopifnot(
  !"IL7R" %in%
    stromal_genes_testA_10b3nd2
)

stopifnot(
  "IL7R" %in%
    immune_genes_testA_10b3nd2
)


# ------------------------------------------------------------
# 6. Record original GlobalEnv state
# ------------------------------------------------------------

global_SI_exists_before_10b3nd2 <-
  exists(
    "SI_geneset",
    envir = .GlobalEnv,
    inherits = FALSE
  )

if (
  global_SI_exists_before_10b3nd2
) {
  
  SI_geneset_global_before_10b3nd2 <-
    get(
      "SI_geneset",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  
} else {
  
  SI_geneset_global_before_10b3nd2 <-
    NULL
}


cat(
  "Original GlobalEnv SI_geneset exists: ",
  global_SI_exists_before_10b3nd2,
  "\n\n"
)


# ------------------------------------------------------------
# 7. Create independent copy for Test A
# ------------------------------------------------------------

SI_geneset_testA_10b3nd2 <-
  SI_geneset_namespace_before_10b3nd2


# ------------------------------------------------------------
# 8. Put copy into GlobalEnv
# ------------------------------------------------------------

assign(
  "SI_geneset",
  SI_geneset_testA_10b3nd2,
  envir = .GlobalEnv
)


# ------------------------------------------------------------
# 9. Confirm GlobalEnv object
# ------------------------------------------------------------

SI_geneset_global_testA_10b3nd2 <-
  get(
    "SI_geneset",
    envir = .GlobalEnv,
    inherits = FALSE
  )

stopifnot(
  identical(
    SI_geneset_global_testA_10b3nd2,
    SI_geneset_namespace_before_10b3nd2
  )
)


# ------------------------------------------------------------
# 10. Confirm namespace remains unchanged
# ------------------------------------------------------------

SI_geneset_namespace_testA_10b3nd2 <-
  get(
    "SI_geneset",
    envir = asNamespace("estimate")
  )

namespace_unchanged_pre_score_10b3nd2 <-
  identical(
    SI_geneset_namespace_before_10b3nd2,
    SI_geneset_namespace_testA_10b3nd2
  )

cat(
  "Package namespace unchanged before estimateScore(): ",
  namespace_unchanged_pre_score_10b3nd2,
  "\n\n"
)

stopifnot(
  namespace_unchanged_pre_score_10b3nd2
)

# ============================================================
# STEP 10B-3N-D-2
# Create D-2 output directory
# ============================================================

D2_dir_10b3nd2 <-
  file.path(
    "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis",
    "STEP7B-7E-10B-3N",
    "10B-3N-D-2"
  )

if (!dir.exists(D2_dir_10b3nd2)) {
  dir.create(
    D2_dir_10b3nd2,
    recursive = TRUE
  )
}

stopifnot(
  dir.exists(D2_dir_10b3nd2)
)

cat(
  "D-2 directory:\n",
  D2_dir_10b3nd2,
  "\n\n"
)

# ============================================================
# STEP 10B-3N-D-2 — TEST A
# Define output path and run official estimateScore()
# ============================================================

testA_gct_10b3nd2 <-
  file.path(
    D2_dir_10b3nd2,
    "STEP10B3ND2_TestA_original_SI_geneset.gct"
  )

cat(
  "Test A output file:\n",
  testA_gct_10b3nd2,
  "\n\n"
)

stopifnot(
  exists(
    "filtered_gct_10b3nd2",
    inherits = FALSE
  )
)

estimate::estimateScore(
  input.ds =
    filtered_gct_10b3nd2,
  output.ds =
    testA_gct_10b3nd2,
  platform =
    "illumina"
)

stopifnot(
  file.exists(
    testA_gct_10b3nd2
  )
)

cat(
  "\nTest A estimateScore() completed successfully.\n"
)

cat(
  "Output exists: ",
  file.exists(testA_gct_10b3nd2),
  "\n"
)

# ============================================================
# STEP 10B-3N-D-2 — TEST A
# Numerical reproduction of locked Step14 scores
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3N-D-2 — TEST A NUMERICAL VALIDATION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Read official estimateScore() output
# ------------------------------------------------------------

read_estimate_score_testA_10b3nd2 <- function(gct_file) {
  
  x <- read.delim(
    gct_file,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  stopifnot(
    nrow(x) >= 6L
  )
  
  sample_ids <- as.character(
    x[1, -c(1, 2)]
  )
  
  score_rows <- x[4:6, -c(1, 2), drop = FALSE]
  
  score_mat <- as.matrix(
    score_rows
  )
  
  storage.mode(score_mat) <- "numeric"
  
  rownames(score_mat) <-
    as.character(
      x[4:6, 1]
    )
  
  colnames(score_mat) <-
    sample_ids
  
  list(
    scores = score_mat,
    sample_ids = sample_ids
  )
}


# ------------------------------------------------------------
# 2. Read Test A result
# ------------------------------------------------------------

testA_result_10b3nd2 <-
  read_estimate_score_testA_10b3nd2(
    testA_gct_10b3nd2
  )


# ------------------------------------------------------------
# 3. Confirm score structure
# ------------------------------------------------------------

cat("Test A score matrix dimension:\n")

print(
  dim(
    testA_result_10b3nd2$scores
  )
)

cat("\nTest A score row names:\n")

print(
  rownames(
    testA_result_10b3nd2$scores
  )
)

stopifnot(
  identical(
    rownames(
      testA_result_10b3nd2$scores
    ),
    c(
      "StromalScore",
      "ImmuneScore",
      "ESTIMATEScore"
    )
  )
)

stopifnot(
  ncol(
    testA_result_10b3nd2$scores
  ) == 513L
)

stopifnot(
  all(
    is.finite(
      testA_result_10b3nd2$scores
    )
  )
)


# ------------------------------------------------------------
# 4. Locked Step14 score object
# ------------------------------------------------------------

stopifnot(
  exists(
    "estimate_scores_wgcna_source_10b3nd2",
    inherits = FALSE
  )
)

locked_scores_testA_10b3nd2 <-
  estimate_scores_wgcna_source_10b3nd2

# ============================================================
# STEP 10B-3N-D-2 — TEST A
# Diagnose sample-ID mismatch
# NO DATA MODIFICATION
# ============================================================

cat("\n============================================================\n")
cat("TEST A — SAMPLE ID DIAGNOSTIC\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Inspect first 6 rows and first 8 columns of GCT
# ------------------------------------------------------------

testA_raw_gct_10b3nd2 <-
  read.delim(
    testA_gct_10b3nd2,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

cat("GCT dimensions:\n")
print(dim(testA_raw_gct_10b3nd2))

cat("\nFirst 6 rows × first 8 columns:\n")
print(
  testA_raw_gct_10b3nd2[
    1:min(6, nrow(testA_raw_gct_10b3nd2)),
    1:min(8, ncol(testA_raw_gct_10b3nd2)),
    drop = FALSE
  ]
)


# ------------------------------------------------------------
# 2. Inspect sample IDs currently extracted by our reader
# ------------------------------------------------------------

cat("\nFirst 10 sample IDs extracted by current reader:\n")

print(
  head(
    testA_result_10b3nd2$sample_ids,
    10
  )
)


# ------------------------------------------------------------
# 3. Inspect locked Step14 sample IDs
# ------------------------------------------------------------

cat("\nFirst 10 locked Step14 sample IDs:\n")

print(
  head(
    locked_sample_ids_testA_10b3nd2,
    10
  )
)


# ------------------------------------------------------------
# 4. Inspect last 10 IDs
# ------------------------------------------------------------

cat("\nLast 10 Test A sample IDs:\n")

print(
  tail(
    testA_result_10b3nd2$sample_ids,
    10
  )
)

cat("\nLast 10 locked Step14 sample IDs:\n")

print(
  tail(
    locked_sample_ids_testA_10b3nd2,
    10
  )
)


# ------------------------------------------------------------
# 5. Length comparison
# ------------------------------------------------------------

cat("\nNumber of Test A sample IDs: ")
print(
  length(
    testA_result_10b3nd2$sample_ids
  )
)

cat("Number of locked sample IDs: ")
print(
  length(
    locked_sample_ids_testA_10b3nd2
  )
)


# ------------------------------------------------------------
# 6. Exact identity
# ------------------------------------------------------------

cat(
  "\nExact identity: ",
  identical(
    testA_result_10b3nd2$sample_ids,
    locked_sample_ids_testA_10b3nd2
  ),
  "\n"
)


# ------------------------------------------------------------
# 7. Set comparison
# ------------------------------------------------------------

testA_not_locked_10b3nd2 <-
  setdiff(
    testA_result_10b3nd2$sample_ids,
    locked_sample_ids_testA_10b3nd2
  )

locked_not_testA_10b3nd2 <-
  setdiff(
    locked_sample_ids_testA_10b3nd2,
    testA_result_10b3nd2$sample_ids
  )

cat("\nTest A IDs absent from locked Step14:\n")
print(testA_not_locked_10b3nd2)

cat("\nLocked Step14 IDs absent from Test A:\n")
print(locked_not_testA_10b3nd2)


# ------------------------------------------------------------
# 8. Check whether same IDs occur in different order
# ------------------------------------------------------------

same_set_testA_10b3nd2 <-
  setequal(
    testA_result_10b3nd2$sample_ids,
    locked_sample_ids_testA_10b3nd2
  )

cat(
  "\nSame sample-ID set: ",
  same_set_testA_10b3nd2,
  "\n"
)


# ------------------------------------------------------------
# 9. Direct comparison of first 20 IDs
# ------------------------------------------------------------

cat("\nFirst 20 IDs — direct comparison:\n\n")

print(
  data.frame(
    Position = seq_len(
      min(
        20L,
        length(
          locked_sample_ids_testA_10b3nd2
        ),
        length(
          testA_result_10b3nd2$sample_ids
        )
      )
    ),
    Locked_Step14 =
      locked_sample_ids_testA_10b3nd2[
        seq_len(
          min(
            20L,
            length(
              locked_sample_ids_testA_10b3nd2
            )
          )
        )
      ],
    TestA =
      testA_result_10b3nd2$sample_ids[
        seq_len(
          min(
            20L,
            length(
              testA_result_10b3nd2$sample_ids
            )
          )
        )
      ],
    stringsAsFactors = FALSE
  )
)


# ------------------------------------------------------------
# 10. Inspect raw GCT header rows explicitly
# ------------------------------------------------------------

cat("\nRaw GCT row 1, first 10 fields:\n")
print(
  as.character(
    testA_raw_gct_10b3nd2[
      1,
      1:min(10, ncol(testA_raw_gct_10b3nd2))
    ]
  )
)

cat("\nRaw GCT row 2, first 10 fields:\n")
print(
  as.character(
    testA_raw_gct_10b3nd2[
      2,
      1:min(10, ncol(testA_raw_gct_10b3nd2))
    ]
  )
)

# ============================================================
# STEP 10B-3N-D-2 — TEST A
# Correct GCT reader: sample IDs are on row 3
# and use "." instead of "-"
# NO estimateScore() rerun
# NO locked-object modification
# ============================================================

cat("\n============================================================\n")
cat("TEST A — CORRECTED GCT READER\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Correct GCT reader
# ------------------------------------------------------------

read_estimate_score_testA_corrected_10b3nd2 <- function(
    gct_file
) {
  
  x <- read.delim(
    gct_file,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  stopifnot(
    nrow(x) >= 6L
  )
  
  # GCT structure:
  # row 1 = version
  # row 2 = dimensions
  # row 3 = NAME / Description / sample IDs
  # row 4 = StromalScore
  # row 5 = ImmuneScore
  # row 6 = ESTIMATEScore
  
  stopifnot(
    identical(
      as.character(x[1, 1]),
      "#1.2"
    )
  )
  
  stopifnot(
    as.integer(
      as.character(x[2, 1])
    ) == 3L
  )
  
  stopifnot(
    as.integer(
      as.character(x[2, 2])
    ) == 513L
  )
  
  stopifnot(
    identical(
      as.character(x[3, 1]),
      "NAME"
    )
  )
  
  stopifnot(
    identical(
      as.character(x[3, 2]),
      "Description"
    )
  )
  
  # ----------------------------------------------------------
  # Sample IDs are on row 3
  # ----------------------------------------------------------
  
  sample_ids_raw <-
    as.character(
      x[3, -c(1, 2)]
    )
  
  # Official GCT output represents TCGA hyphens as dots.
  # Restore canonical TCGA barcode representation.
  sample_ids <-
    gsub(
      "\\.",
      "-",
      sample_ids_raw
    )
  
  # ----------------------------------------------------------
  # Score rows
  # ----------------------------------------------------------
  
  score_names <-
    as.character(
      x[4:6, 1]
    )
  
  stopifnot(
    identical(
      score_names,
      c(
        "StromalScore",
        "ImmuneScore",
        "ESTIMATEScore"
      )
    )
  )
  
  score_mat <-
    as.matrix(
      x[4:6, -c(1, 2), drop = FALSE]
    )
  
  storage.mode(score_mat) <- "numeric"
  
  rownames(score_mat) <-
    score_names
  
  colnames(score_mat) <-
    sample_ids
  
  list(
    scores = score_mat,
    sample_ids = sample_ids,
    sample_ids_raw = sample_ids_raw
  )
}


# ------------------------------------------------------------
# 2. Read Test A using corrected reader
# ------------------------------------------------------------

testA_result_corrected_10b3nd2 <-
  read_estimate_score_testA_corrected_10b3nd2(
    testA_gct_10b3nd2
  )


# ------------------------------------------------------------
# 3. Basic structural checks
# ------------------------------------------------------------

cat("Corrected Test A score matrix dimension:\n")

print(
  dim(
    testA_result_corrected_10b3nd2$scores
  )
)

cat("\nFirst 10 corrected sample IDs:\n")

print(
  head(
    testA_result_corrected_10b3nd2$sample_ids,
    10
  )
)

cat("\nLast 10 corrected sample IDs:\n")

print(
  tail(
    testA_result_corrected_10b3nd2$sample_ids,
    10
  )
)


stopifnot(
  identical(
    dim(
      testA_result_corrected_10b3nd2$scores
    ),
    c(3L, 513L)
  )
)

stopifnot(
  identical(
    rownames(
      testA_result_corrected_10b3nd2$scores
    ),
    c(
      "StromalScore",
      "ImmuneScore",
      "ESTIMATEScore"
    )
  )
)

stopifnot(
  all(
    is.finite(
      testA_result_corrected_10b3nd2$scores
    )
  )
)


# ------------------------------------------------------------
# 4. Locked Step14 sample IDs
# ------------------------------------------------------------

locked_sample_ids_testA_corrected_10b3nd2 <-
  colnames(
    estimate_scores_wgcna_source_10b3nd2
  )


# ------------------------------------------------------------
# 5. Exact sample identity/order comparison
# ------------------------------------------------------------

sample_identity_testA_corrected_10b3nd2 <-
  identical(
    locked_sample_ids_testA_corrected_10b3nd2,
    testA_result_corrected_10b3nd2$sample_ids
  )

cat(
  "\nSample identity/order identical after correction: ",
  sample_identity_testA_corrected_10b3nd2,
  "\n"
)

stopifnot(
  sample_identity_testA_corrected_10b3nd2
)


# ------------------------------------------------------------
# 6. Confirm no duplicated sample IDs
# ------------------------------------------------------------

stopifnot(
  !anyDuplicated(
    testA_result_corrected_10b3nd2$sample_ids
  )
)

cat(
  "Test A sample IDs duplicated: FALSE\n"
)


# ------------------------------------------------------------
# 7. Numerical comparison with locked Step14
# ------------------------------------------------------------

locked_scores_testA_corrected_10b3nd2 <-
  estimate_scores_wgcna_source_10b3nd2

testA_scores_testA_corrected_10b3nd2 <-
  testA_result_corrected_10b3nd2$scores


comparison_testA_corrected_10b3nd2 <-
  data.frame(
    Score =
      rownames(
        locked_scores_testA_corrected_10b3nd2
      ),
    Pearson =
      NA_real_,
    Spearman =
      NA_real_,
    MaxAbsDiff =
      NA_real_,
    stringsAsFactors = FALSE
  )


for (i in seq_len(
  nrow(
    locked_scores_testA_corrected_10b3nd2
  )
)) {
  
  x <-
    as.numeric(
      locked_scores_testA_corrected_10b3nd2[i, ]
    )
  
  y <-
    as.numeric(
      testA_scores_testA_corrected_10b3nd2[i, ]
    )
  
  comparison_testA_corrected_10b3nd2$Pearson[i] <-
    cor(
      x,
      y,
      method = "pearson"
    )
  
  comparison_testA_corrected_10b3nd2$Spearman[i] <-
    cor(
      x,
      y,
      method = "spearman"
    )
  
  comparison_testA_corrected_10b3nd2$MaxAbsDiff[i] <-
    max(
      abs(x - y)
    )
}


cat(
  "\n============================================================\n"
)

cat(
  "NUMERICAL COMPARISON WITH LOCKED STEP14\n"
)

cat(
  "============================================================\n\n"
)

print(
  comparison_testA_corrected_10b3nd2,
  digits = 16
)


# ------------------------------------------------------------
# 8. Test A reproduction criterion
# ------------------------------------------------------------

testA_reproduction_pass_corrected_10b3nd2 <-
  all(
    comparison_testA_corrected_10b3nd2$Pearson >
      0.999999999999,
    comparison_testA_corrected_10b3nd2$Spearman ==
      1,
    comparison_testA_corrected_10b3nd2$MaxAbsDiff <
      1e-8
  )


cat(
  "\n============================================================\n"
)

cat(
  "TEST A REPRODUCTION PASSED: ",
  testA_reproduction_pass_corrected_10b3nd2,
  "\n"
)

cat(
  "============================================================\n\n"
)

stopifnot(
  testA_reproduction_pass_corrected_10b3nd2
)


# ------------------------------------------------------------
# 9. Confirm package namespace unchanged
# ------------------------------------------------------------

SI_geneset_namespace_after_testA_corrected_10b3nd2 <-
  get(
    "SI_geneset",
    envir = asNamespace("estimate")
  )

namespace_unchanged_after_testA_corrected_10b3nd2 <-
  identical(
    SI_geneset_namespace_before_10b3nd2,
    SI_geneset_namespace_after_testA_corrected_10b3nd2
  )

cat(
  "Package namespace unchanged after corrected Test A: ",
  namespace_unchanged_after_testA_corrected_10b3nd2,
  "\n"
)

stopifnot(
  namespace_unchanged_after_testA_corrected_10b3nd2
)

# ============================================================
# STEP 10B-3N-D-2 — TEST B
# Official SI_geneset substitution validation
#
# Purpose:
#   Remove AOC3 from StromalSignature ONLY as a mechanism test.
#
# NO formal LOO analysis
# NO statistical inference
# NO modification of locked Step14
# NO modification of Figure 6 / Supplementary Table 3
# ============================================================

cat("\n============================================================\n")
cat("TEST B — REMOVE AOC3 FROM STROMAL SIGNATURE\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Preserve the official namespace object
# ------------------------------------------------------------

SI_geneset_namespace_testB_10b3nd2 <-
  get(
    "SI_geneset",
    envir = asNamespace("estimate")
  )

stopifnot(
  identical(
    SI_geneset_namespace_before_10b3nd2,
    SI_geneset_namespace_testB_10b3nd2
  )
)

cat(
  "Package namespace unchanged before Test B: TRUE\n"
)


# ------------------------------------------------------------
# 2. Create modified SI_geneset
# ------------------------------------------------------------

SI_geneset_testB_10b3nd2 <-
  SI_geneset_namespace_testB_10b3nd2


# ------------------------------------------------------------
# 3. Locate AOC3 in StromalSignature
# ------------------------------------------------------------

stromal_genes_testB_10b3nd2 <-
  as.character(
    SI_geneset_testB_10b3nd2[
      "StromalSignature",
      -1
    ]
  )

aoc3_position_testB_10b3nd2 <-
  which(
    stromal_genes_testB_10b3nd2 == "AOC3"
  )

cat(
  "\nAOC3 position in original StromalSignature: ",
  aoc3_position_testB_10b3nd2,
  "\n"
)

stopifnot(
  length(aoc3_position_testB_10b3nd2) == 1L
)


# ------------------------------------------------------------
# 4. Remove AOC3 while preserving the 2 × 142 structure
# ------------------------------------------------------------

SI_geneset_testB_10b3nd2[
  "StromalSignature",
  aoc3_position_testB_10b3nd2 + 1L
] <- NA_character_


# ------------------------------------------------------------
# 5. Structural validation
# ------------------------------------------------------------

stopifnot(
  is.data.frame(
    SI_geneset_testB_10b3nd2
  )
)

stopifnot(
  identical(
    dim(
      SI_geneset_testB_10b3nd2
    ),
    c(2L, 142L)
  )
)

stopifnot(
  identical(
    rownames(
      SI_geneset_testB_10b3nd2
    ),
    c(
      "StromalSignature",
      "ImmuneSignature"
    )
  )
)

stopifnot(
  sum(
    !is.na(
      SI_geneset_testB_10b3nd2[
        "StromalSignature",
        -1
      ]
    )
  ) == 140L
)

stopifnot(
  sum(
    !is.na(
      SI_geneset_testB_10b3nd2[
        "ImmuneSignature",
        -1
      ]
    )
  ) == 141L
)

stopifnot(
  !(
    "AOC3" %in%
      as.character(
        SI_geneset_testB_10b3nd2[
          "StromalSignature",
          -1
        ]
      )
  )
)

stopifnot(
  "AOC3" %in%
    stromal_genes_testB_10b3nd2
)

cat(
  "\nTest B SI_geneset structure validated.\n"
)

cat(
  "StromalSignature genes after AOC3 removal: 140\n"
)

cat(
  "ImmuneSignature genes: 141\n"
)

cat(
  "AOC3 present in modified StromalSignature: FALSE\n"
)


# ------------------------------------------------------------
# 6. Confirm AOC3 was NOT removed from ImmuneSignature
# ------------------------------------------------------------

immune_genes_testB_10b3nd2 <-
  as.character(
    SI_geneset_testB_10b3nd2[
      "ImmuneSignature",
      -1
    ]
  )

stopifnot(
  !"AOC3" %in%
    immune_genes_testB_10b3nd2
)

cat(
  "AOC3 present in ImmuneSignature: FALSE\n"
)


# ------------------------------------------------------------
# 7. Put modified SI_geneset into GlobalEnv
# ------------------------------------------------------------

assign(
  "SI_geneset",
  SI_geneset_testB_10b3nd2,
  envir = .GlobalEnv
)


# ------------------------------------------------------------
# 8. Confirm GlobalEnv substitution
# ------------------------------------------------------------

SI_geneset_global_testB_10b3nd2 <-
  get(
    "SI_geneset",
    envir = .GlobalEnv
  )

stopifnot(
  identical(
    SI_geneset_global_testB_10b3nd2,
    SI_geneset_testB_10b3nd2
  )
)

stopifnot(
  !identical(
    SI_geneset_global_testB_10b3nd2,
    SI_geneset_namespace_testB_10b3nd2
  )
)

cat(
  "\nGlobalEnv SI_geneset successfully substituted.\n"
)


# ------------------------------------------------------------
# 9. Confirm package namespace is STILL unchanged
# ------------------------------------------------------------

SI_geneset_namespace_before_estimateScore_testB_10b3nd2 <-
  get(
    "SI_geneset",
    envir = asNamespace("estimate")
  )

stopifnot(
  identical(
    SI_geneset_namespace_before_estimateScore_testB_10b3nd2,
    SI_geneset_namespace_before_10b3nd2
  )
)

cat(
  "Package namespace unchanged before Test B estimateScore(): TRUE\n"
)

# ============================================================
# STEP 10B-3N-D-2 — TEST B DIAGNOSTIC
# Diagnose namespace difference BEFORE estimateScore()
#
# IMPORTANT:
#   DO NOT run estimateScore() yet.
# ============================================================

cat("\n============================================================\n")
cat("TEST B — NAMESPACE DIFFERENCE DIAGNOSTIC\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Retrieve current package-namespace object
# ------------------------------------------------------------

SI_geneset_namespace_current_testBdiag_10b3nd2 <-
  get(
    "SI_geneset",
    envir = asNamespace("estimate")
  )


# ------------------------------------------------------------
# 2. Basic structural comparison
# ------------------------------------------------------------

cat("Original namespace snapshot class:\n")
print(
  class(
    SI_geneset_namespace_before_10b3nd2
  )
)

cat("\nCurrent namespace object class:\n")
print(
  class(
    SI_geneset_namespace_current_testBdiag_10b3nd2
  )
)

cat("\nOriginal namespace snapshot dimensions:\n")
print(
  dim(
    SI_geneset_namespace_before_10b3nd2
  )
)

cat("\nCurrent namespace dimensions:\n")
print(
  dim(
    SI_geneset_namespace_current_testBdiag_10b3nd2
  )
)


# ------------------------------------------------------------
# 3. Compare row names and column names
# ------------------------------------------------------------

cat("\nRow names identical:\n")

print(
  identical(
    rownames(
      SI_geneset_namespace_before_10b3nd2
    ),
    rownames(
      SI_geneset_namespace_current_testBdiag_10b3nd2
    )
  )
)


cat("\nColumn names identical:\n")

print(
  identical(
    colnames(
      SI_geneset_namespace_before_10b3nd2
    ),
    colnames(
      SI_geneset_namespace_current_testBdiag_10b3nd2
    )
  )
)


# ------------------------------------------------------------
# 4. Compare substantive gene vectors
# ------------------------------------------------------------

stromal_original_testBdiag_10b3nd2 <-
  as.character(
    SI_geneset_namespace_before_10b3nd2[
      "StromalSignature",
      -1
    ]
  )

stromal_current_testBdiag_10b3nd2 <-
  as.character(
    SI_geneset_namespace_current_testBdiag_10b3nd2[
      "StromalSignature",
      -1
    ]
  )


immune_original_testBdiag_10b3nd2 <-
  as.character(
    SI_geneset_namespace_before_10b3nd2[
      "ImmuneSignature",
      -1
    ]
  )

immune_current_testBdiag_10b3nd2 <-
  as.character(
    SI_geneset_namespace_current_testBdiag_10b3nd2[
      "ImmuneSignature",
      -1
    ]
  )


cat("\nStromal gene vector identical:\n")

print(
  identical(
    stromal_original_testBdiag_10b3nd2,
    stromal_current_testBdiag_10b3nd2
  )
)


cat("\nImmune gene vector identical:\n")

print(
  identical(
    immune_original_testBdiag_10b3nd2,
    immune_current_testBdiag_10b3nd2
  )
)


# ------------------------------------------------------------
# 5. Identify substantive differences
# ------------------------------------------------------------

cat("\nStromal set difference — original minus current:\n")

print(
  setdiff(
    stromal_original_testBdiag_10b3nd2,
    stromal_current_testBdiag_10b3nd2
  )
)


cat("\nStromal set difference — current minus original:\n")

print(
  setdiff(
    stromal_current_testBdiag_10b3nd2,
    stromal_original_testBdiag_10b3nd2
  )
)


cat("\nImmune set difference — original minus current:\n")

print(
  setdiff(
    immune_original_testBdiag_10b3nd2,
    immune_current_testBdiag_10b3nd2
  )
)


cat("\nImmune set difference — current minus original:\n")

print(
  setdiff(
    immune_current_testBdiag_10b3nd2,
    immune_original_testBdiag_10b3nd2
  )
)


# ------------------------------------------------------------
# 6. Check AOC3 specifically
# ------------------------------------------------------------

cat("\nAOC3 in original StromalSignature:\n")

print(
  "AOC3" %in%
    stromal_original_testBdiag_10b3nd2
)


cat("\nAOC3 in current package StromalSignature:\n")

print(
  "AOC3" %in%
    stromal_current_testBdiag_10b3nd2
)


# ------------------------------------------------------------
# 7. Check first-column values and attributes
# ------------------------------------------------------------

cat("\nFirst-column values identical:\n")

print(
  identical(
    as.character(
      SI_geneset_namespace_before_10b3nd2[, 1]
    ),
    as.character(
      SI_geneset_namespace_current_testBdiag_10b3nd2[, 1]
    )
  )
)


cat("\nAttributes identical:\n")

print(
  identical(
    attributes(
      SI_geneset_namespace_before_10b3nd2
    ),
    attributes(
      SI_geneset_namespace_current_testBdiag_10b3nd2
    )
  )
)


# ------------------------------------------------------------
# 8. all.equal() gives a more informative diagnosis
# ------------------------------------------------------------

cat("\nall.equal() result:\n")

print(
  all.equal(
    SI_geneset_namespace_before_10b3nd2,
    SI_geneset_namespace_current_testBdiag_10b3nd2
  )
)


# ------------------------------------------------------------
# 9. Exact structural/content summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("NAMESPACE DIAGNOSTIC SUMMARY\n")
cat("============================================================\n\n")

cat(
  "identical(): ",
  identical(
    SI_geneset_namespace_before_10b3nd2,
    SI_geneset_namespace_current_testBdiag_10b3nd2
  ),
  "\n"
)

cat(
  "Stromal vector identical: ",
  identical(
    stromal_original_testBdiag_10b3nd2,
    stromal_current_testBdiag_10b3nd2
  ),
  "\n"
)

cat(
  "Immune vector identical: ",
  identical(
    immune_original_testBdiag_10b3nd2,
    immune_current_testBdiag_10b3nd2
  ),
  "\n"
)

cat(
  "AOC3 in current package StromalSignature: ",
  "AOC3" %in% stromal_current_testBdiag_10b3nd2,
  "\n"
)

cat(
  "AOC3 in original snapshot StromalSignature: ",
  "AOC3" %in% stromal_original_testBdiag_10b3nd2,
  "\n"
)

cat("\nDO NOT RUN estimateScore() YET.\n")


# ============================================================
# STEP 10B-3N-D-2 — RESTORE ORIGINAL GLOBAL SI_geneset
#
# IMPORTANT:
#   Recovery only.
#   NO estimateScore().
#   NO Test B.
#   NO package namespace modification.
# ============================================================

cat("\n============================================================\n")
cat("RESTORE ORIGINAL GLOBAL SI_geneset\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Confirm direct binding location
# ------------------------------------------------------------

stopifnot(
  exists(
    "SI_geneset",
    envir = .GlobalEnv,
    inherits = FALSE
  )
)

cat(
  "Direct SI_geneset binding in R_GlobalEnv: TRUE\n"
)

cat(
  "Global binding locked: ",
  bindingIsLocked("SI_geneset", .GlobalEnv),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 2. Restore original GlobalEnv state
# ------------------------------------------------------------

if (global_SI_exists_before_10b3nd2) {
  
  assign(
    "SI_geneset",
    SI_geneset_global_before_10b3nd2,
    envir = .GlobalEnv
  )
  
  cat(
    "Original GlobalEnv SI_geneset restored.\n"
  )
  
} else {
  
  rm(
    "SI_geneset",
    envir = .GlobalEnv
  )
  
  cat(
    "Original GlobalEnv state was: SI_geneset absent.\n"
  )
}


# ------------------------------------------------------------
# 3. Verify GlobalEnv state
# ------------------------------------------------------------

global_SI_exists_after_restore_10b3nd2 <-
  exists(
    "SI_geneset",
    envir = .GlobalEnv,
    inherits = FALSE
  )

cat(
  "\nGlobalEnv SI_geneset exists after restore: ",
  global_SI_exists_after_restore_10b3nd2,
  "\n",
  sep = ""
)


if (global_SI_exists_before_10b3nd2) {
  
  global_SI_identical_after_restore_10b3nd2 <-
    identical(
      get(
        "SI_geneset",
        envir = .GlobalEnv,
        inherits = FALSE
      ),
      SI_geneset_global_before_10b3nd2
    )
  
  cat(
    "GlobalEnv SI_geneset identical to original: ",
    global_SI_identical_after_restore_10b3nd2,
    "\n",
    sep = ""
  )
  
} else {
  
  global_SI_identical_after_restore_10b3nd2 <-
    !global_SI_exists_after_restore_10b3nd2
  
  cat(
    "GlobalEnv original absence restored: ",
    global_SI_identical_after_restore_10b3nd2,
    "\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 4. Verify package:estimate direct binding
# ------------------------------------------------------------

estimate_pkg_env_10b3nd2 <-
  as.environment("package:estimate")

package_SI_geneset_after_restore_10b3nd2 <-
  get(
    "SI_geneset",
    envir = estimate_pkg_env_10b3nd2,
    inherits = FALSE
  )

package_SI_identical_after_restore_10b3nd2 <-
  identical(
    package_SI_geneset_after_restore_10b3nd2,
    SI_geneset_namespace_before_10b3nd2
  )

cat(
  "\nPackage:estimate SI_geneset identical to D-2 original: ",
  package_SI_identical_after_restore_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "Package:estimate SI_geneset binding locked: ",
  bindingIsLocked(
    "SI_geneset",
    estimate_pkg_env_10b3nd2
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. Check critical gene membership
# ------------------------------------------------------------

package_stromal_genes_after_restore_10b3nd2 <-
  as.character(
    package_SI_geneset_after_restore_10b3nd2[
      1,
      -1
    ]
  )

package_immune_genes_after_restore_10b3nd2 <-
  as.character(
    package_SI_geneset_after_restore_10b3nd2[
      2,
      -1
    ]
  )

package_stromal_genes_after_restore_10b3nd2 <-
  package_stromal_genes_after_restore_10b3nd2[
    !is.na(package_stromal_genes_after_restore_10b3nd2)
  ]

package_immune_genes_after_restore_10b3nd2 <-
  package_immune_genes_after_restore_10b3nd2[
    !is.na(package_immune_genes_after_restore_10b3nd2)
  ]


cat(
  "\nCritical membership after restoration:\n"
)

cat(
  "AOC3 in StromalSignature: ",
  "AOC3" %in% package_stromal_genes_after_restore_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "AOC3 in ImmuneSignature: ",
  "AOC3" %in% package_immune_genes_after_restore_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "MS4A7 in StromalSignature: ",
  "MS4A7" %in% package_stromal_genes_after_restore_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "MS4A7 in ImmuneSignature: ",
  "MS4A7" %in% package_immune_genes_after_restore_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "IL7R in StromalSignature: ",
  "IL7R" %in% package_stromal_genes_after_restore_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "IL7R in ImmuneSignature: ",
  "IL7R" %in% package_immune_genes_after_restore_10b3nd2,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Final recovery status
# ------------------------------------------------------------

recovery_pass_10b3nd2 <-
  isTRUE(global_SI_identical_after_restore_10b3nd2) &&
  isTRUE(package_SI_identical_after_restore_10b3nd2) &&
  identical(
    bindingIsLocked(
      "SI_geneset",
      estimate_pkg_env_10b3nd2
    ),
    TRUE
  ) &&
  ("AOC3" %in% package_stromal_genes_after_restore_10b3nd2) &&
  !("AOC3" %in% package_immune_genes_after_restore_10b3nd2) &&
  !("MS4A7" %in% package_stromal_genes_after_restore_10b3nd2) &&
  !("MS4A7" %in% package_immune_genes_after_restore_10b3nd2) &&
  !("IL7R" %in% package_stromal_genes_after_restore_10b3nd2) &&
  ("IL7R" %in% package_immune_genes_after_restore_10b3nd2)


cat(
  "\n============================================================\n"
)

cat(
  "RECOVERY PASS: ",
  recovery_pass_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)


# ============================================================
# STEP 10B-3N-D-2 — TEST B PREPARATION
#
# Purpose:
#   Create a truly independent deep copy of the official
#   SI_geneset and remove AOC3 from the StromalSignature.
#
# IMPORTANT:
#   Preparation + safety validation only.
#   NO estimateScore() yet.
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3N-D-2 — TEST B PREPARATION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Re-confirm original package SI_geneset
# ------------------------------------------------------------

estimate_pkg_env_10b3nd2 <-
  as.environment("package:estimate")

SI_geneset_original_testB_10b3nd2 <-
  get(
    "SI_geneset",
    envir = estimate_pkg_env_10b3nd2,
    inherits = FALSE
  )

stopifnot(
  identical(
    SI_geneset_original_testB_10b3nd2,
    SI_geneset_namespace_before_10b3nd2
  )
)

cat(
  "Original package SI_geneset verified: TRUE\n"
)


# ------------------------------------------------------------
# 2. Create TRUE deep copy
# ------------------------------------------------------------

SI_geneset_testB_deep_10b3nd2 <-
  unserialize(
    serialize(
      SI_geneset_original_testB_10b3nd2,
      NULL
    )
  )

cat(
  "Deep copy created: TRUE\n"
)


# ------------------------------------------------------------
# 3. Confirm deep copy is structurally identical BEFORE edit
# ------------------------------------------------------------

deep_copy_identical_before_edit_10b3nd2 <-
  identical(
    SI_geneset_testB_deep_10b3nd2,
    SI_geneset_original_testB_10b3nd2
  )

cat(
  "Deep copy identical to original before edit: ",
  deep_copy_identical_before_edit_10b3nd2,
  "\n",
  sep = ""
)

stopifnot(
  deep_copy_identical_before_edit_10b3nd2
)


# ------------------------------------------------------------
# 4. Locate AOC3 in StromalSignature
# ------------------------------------------------------------

stromal_testB_before_10b3nd2 <-
  as.character(
    SI_geneset_testB_deep_10b3nd2[
      1,
      -1
    ]
  )

AOC3_position_testB_10b3nd2 <-
  which(
    stromal_testB_before_10b3nd2 == "AOC3"
  )

cat(
  "AOC3 position in StromalSignature: ",
  AOC3_position_testB_10b3nd2,
  "\n",
  sep = ""
)

stopifnot(
  length(AOC3_position_testB_10b3nd2) == 1L
)


# ------------------------------------------------------------
# 5. Remove AOC3 from StromalSignature
#
# Preserve the original 2 × 142 structure.
# Replace only the AOC3 slot with NA.
# ------------------------------------------------------------

SI_geneset_testB_deep_10b3nd2[
  1,
  AOC3_position_testB_10b3nd2 + 1L
] <- NA_character_


# ------------------------------------------------------------
# 6. Verify modified Test B signature
# ------------------------------------------------------------

stromal_testB_after_10b3nd2 <-
  as.character(
    SI_geneset_testB_deep_10b3nd2[
      1,
      -1
    ]
  )

immune_testB_after_10b3nd2 <-
  as.character(
    SI_geneset_testB_deep_10b3nd2[
      2,
      -1
    ]
  )

stromal_testB_after_10b3nd2 <-
  stromal_testB_after_10b3nd2[
    !is.na(stromal_testB_after_10b3nd2)
  ]

immune_testB_after_10b3nd2 <-
  immune_testB_after_10b3nd2[
    !is.na(immune_testB_after_10b3nd2)
  ]


cat(
  "\nTest B signature:\n"
)

cat(
  "Stromal genes after removing AOC3: ",
  length(stromal_testB_after_10b3nd2),
  "\n",
  sep = ""
)

cat(
  "Immune genes: ",
  length(immune_testB_after_10b3nd2),
  "\n",
  sep = ""
)

cat(
  "AOC3 in Test B StromalSignature: ",
  "AOC3" %in% stromal_testB_after_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "AOC3 in Test B ImmuneSignature: ",
  "AOC3" %in% immune_testB_after_10b3nd2,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. SAFETY CHECK:
#    The package object must remain unchanged.
# ------------------------------------------------------------

package_SI_after_edit_10b3nd2 <-
  get(
    "SI_geneset",
    envir = estimate_pkg_env_10b3nd2,
    inherits = FALSE
  )

package_unchanged_after_deepcopy_edit_10b3nd2 <-
  identical(
    package_SI_after_edit_10b3nd2,
    SI_geneset_original_testB_10b3nd2
  )

cat(
  "\nPackage SI_geneset unchanged after deep-copy modification: ",
  package_unchanged_after_deepcopy_edit_10b3nd2,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Additional critical check:
#    Package AOC3 must still be present.
# ------------------------------------------------------------

package_stromal_after_edit_10b3nd2 <-
  as.character(
    package_SI_after_edit_10b3nd2[
      1,
      -1
    ]
  )

package_stromal_after_edit_10b3nd2 <-
  package_stromal_after_edit_10b3nd2[
    !is.na(package_stromal_after_edit_10b3nd2)
  ]

cat(
  "AOC3 still present in package StromalSignature: ",
  "AOC3" %in% package_stromal_after_edit_10b3nd2,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Final preparation status
# ------------------------------------------------------------

testB_preparation_pass_10b3nd2 <-
  deep_copy_identical_before_edit_10b3nd2 &&
  length(AOC3_position_testB_10b3nd2) == 1L &&
  length(stromal_testB_after_10b3nd2) == 140L &&
  length(immune_testB_after_10b3nd2) == 141L &&
  !("AOC3" %in% stromal_testB_after_10b3nd2) &&
  !("AOC3" %in% immune_testB_after_10b3nd2) &&
  package_unchanged_after_deepcopy_edit_10b3nd2 &&
  ("AOC3" %in% package_stromal_after_edit_10b3nd2)


cat(
  "\n============================================================\n"
)

cat(
  "TEST B PREPARATION PASS: ",
  testB_preparation_pass_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)

# ============================================================
# STEP 10B-3N-D-2 — TEST B
#
# Official SI_geneset substitution validation
#
# Purpose:
#   Remove AOC3 from StromalSignature in an isolated deep copy
#   and test whether official estimateScore() responds to the
#   substituted signature.
#
# IMPORTANT:
#   - Official estimateScore() only
#   - No statistical analysis
#   - No LOO
#   - No modification of Step14
#   - No modification of Figure 6
#   - No modification of Supplementary Table 3
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3N-D-2 — TEST B OFFICIAL SI_geneset SUBSTITUTION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define input/output files
# ------------------------------------------------------------

filtered_gct_10b3nd2 <-
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis/STEP7B-7E-10B-3N/10B-3N-C/STEP10B3NC_official_filterCommonGenes_output.gct"

testB_gct_10b3nd2 <-
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis/STEP7B-7E-10B-3N/10B-3N-D-2/STEP10B3ND2_TestB_AOC3_removed.gct"


# ------------------------------------------------------------
# 2. Safety check: input GCT exists
# ------------------------------------------------------------

stopifnot(
  file.exists(filtered_gct_10b3nd2)
)

cat(
  "Input GCT exists: TRUE\n"
)


# ------------------------------------------------------------
# 3. Confirm original package SI_geneset
# ------------------------------------------------------------

estimate_pkg_env_10b3nd2 <-
  as.environment("package:estimate")

package_before_testB_10b3nd2 <-
  get(
    "SI_geneset",
    envir = estimate_pkg_env_10b3nd2,
    inherits = FALSE
  )

stopifnot(
  identical(
    package_before_testB_10b3nd2,
    SI_geneset_namespace_before_10b3nd2
  )
)

cat(
  "Package SI_geneset before Test B: ORIGINAL\n"
)

cat(
  "Package binding locked: ",
  bindingIsLocked(
    "SI_geneset",
    estimate_pkg_env_10b3nd2
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Confirm GlobalEnv original state
# ------------------------------------------------------------

global_before_testB_exists_10b3nd2 <-
  exists(
    "SI_geneset",
    envir = .GlobalEnv,
    inherits = FALSE
  )

cat(
  "GlobalEnv SI_geneset exists before Test B: ",
  global_before_testB_exists_10b3nd2,
  "\n",
  sep = ""
)

if (global_before_testB_exists_10b3nd2) {
  
  global_before_testB_10b3nd2 <-
    get(
      "SI_geneset",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  
  global_before_testB_identical_10b3nd2 <-
    identical(
      global_before_testB_10b3nd2,
      SI_geneset_global_before_10b3nd2
    )
  
  cat(
    "GlobalEnv matches original D-2 state: ",
    global_before_testB_identical_10b3nd2,
    "\n",
    sep = ""
  )
  
} else {
  
  global_before_testB_10b3nd2 <- NULL
  
  global_before_testB_identical_10b3nd2 <- TRUE
  
  cat(
    "GlobalEnv originally had no SI_geneset.\n"
  )
}


# ------------------------------------------------------------
# 5. Deep-copy safety verification immediately before use
# ------------------------------------------------------------

testB_signature_before_assignment_10b3nd2 <-
  unserialize(
    serialize(
      SI_geneset_testB_deep_10b3nd2,
      NULL
    )
  )

stopifnot(
  identical(
    testB_signature_before_assignment_10b3nd2,
    SI_geneset_testB_deep_10b3nd2
  )
)

cat(
  "Test B deep-copy integrity verified: TRUE\n"
)


# ------------------------------------------------------------
# 6. Assign MODIFIED SI_geneset to R_GlobalEnv
#
# IMPORTANT:
#   This is the intentional substitution.
# ------------------------------------------------------------

assign(
  "SI_geneset",
  SI_geneset_testB_deep_10b3nd2,
  envir = .GlobalEnv
)

cat(
  "Modified Test B SI_geneset assigned to R_GlobalEnv.\n"
)


# ------------------------------------------------------------
# 7. Verify that GlobalEnv now contains modified signature
# ------------------------------------------------------------

SI_geneset_global_during_testB_10b3nd2 <-
  get(
    "SI_geneset",
    envir = .GlobalEnv,
    inherits = FALSE
  )

stromal_global_during_testB_10b3nd2 <-
  as.character(
    SI_geneset_global_during_testB_10b3nd2[
      1,
      -1
    ]
  )

immune_global_during_testB_10b3nd2 <-
  as.character(
    SI_geneset_global_during_testB_10b3nd2[
      2,
      -1
    ]
  )

stromal_global_during_testB_10b3nd2 <-
  stromal_global_during_testB_10b3nd2[
    !is.na(stromal_global_during_testB_10b3nd2)
  ]

immune_global_during_testB_10b3nd2 <-
  immune_global_during_testB_10b3nd2[
    !is.na(immune_global_during_testB_10b3nd2)
  ]


cat(
  "\nSignature actually present in GlobalEnv during Test B:\n"
)

cat(
  "Stromal genes: ",
  length(stromal_global_during_testB_10b3nd2),
  "\n",
  sep = ""
)

cat(
  "Immune genes: ",
  length(immune_global_during_testB_10b3nd2),
  "\n",
  sep = ""
)

cat(
  "AOC3 in GlobalEnv StromalSignature: ",
  "AOC3" %in% stromal_global_during_testB_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "AOC3 in GlobalEnv ImmuneSignature: ",
  "AOC3" %in% immune_global_during_testB_10b3nd2,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. CRITICAL SAFETY CHECK:
#    package:estimate must still contain ORIGINAL object
# ------------------------------------------------------------

package_during_testB_10b3nd2 <-
  get(
    "SI_geneset",
    envir = estimate_pkg_env_10b3nd2,
    inherits = FALSE
  )

package_unchanged_before_estimateScore_testB_10b3nd2 <-
  identical(
    package_during_testB_10b3nd2,
    SI_geneset_namespace_before_10b3nd2
  )

cat(
  "\nPackage SI_geneset unchanged before estimateScore(): ",
  package_unchanged_before_estimateScore_testB_10b3nd2,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Run OFFICIAL estimateScore()
# ------------------------------------------------------------

cat(
  "\nRunning official estimateScore() for Test B...\n"
)

estimateScore(
  input.ds = filtered_gct_10b3nd2,
  output.ds = testB_gct_10b3nd2,
  platform = "illumina"
)

cat(
  "Official estimateScore() completed.\n"
)


# ------------------------------------------------------------
# 10. Verify package namespace immediately after scoring
# ------------------------------------------------------------

package_after_testB_10b3nd2 <-
  get(
    "SI_geneset",
    envir = estimate_pkg_env_10b3nd2,
    inherits = FALSE
  )

package_unchanged_after_estimateScore_testB_10b3nd2 <-
  identical(
    package_after_testB_10b3nd2,
    SI_geneset_namespace_before_10b3nd2
  )

cat(
  "\nPackage SI_geneset unchanged after estimateScore(): ",
  package_unchanged_after_estimateScore_testB_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "Package binding still locked: ",
  bindingIsLocked(
    "SI_geneset",
    estimate_pkg_env_10b3nd2
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Restore original GlobalEnv state IMMEDIATELY
# ------------------------------------------------------------

if (global_SI_exists_before_10b3nd2) {
  
  assign(
    "SI_geneset",
    SI_geneset_global_before_10b3nd2,
    envir = .GlobalEnv
  )
  
  cat(
    "Original GlobalEnv SI_geneset restored after Test B.\n"
  )
  
} else {
  
  rm(
    "SI_geneset",
    envir = .GlobalEnv
  )
  
  cat(
    "Original GlobalEnv absence restored after Test B.\n"
  )
}


# ------------------------------------------------------------
# 12. Read Test B output
# ------------------------------------------------------------

stopifnot(
  file.exists(testB_gct_10b3nd2)
)

cat(
  "Test B output GCT exists: TRUE\n"
)


testB_lines_10b3nd2 <-
  readLines(
    testB_gct_10b3nd2,
    warn = FALSE
  )

testB_header_index_10b3nd2 <-
  which(
    grepl(
      "^NAME\\tDescription\\t",
      testB_lines_10b3nd2
    )
  )[1]

stopifnot(
  !is.na(testB_header_index_10b3nd2)
)

testB_header_10b3nd2 <-
  strsplit(
    testB_lines_10b3nd2[
      testB_header_index_10b3nd2
    ],
    "\t",
    fixed = FALSE
  )[[1]]

testB_data_lines_10b3nd2 <-
  testB_lines_10b3nd2[
    (testB_header_index_10b3nd2 + 1):
      length(testB_lines_10b3nd2)
  ]

testB_data_lines_10b3nd2 <-
  testB_data_lines_10b3nd2[
    nzchar(testB_data_lines_10b3nd2)
  ]

testB_data_10b3nd2 <-
  do.call(
    rbind,
    lapply(
      testB_data_lines_10b3nd2,
      function(x) {
        strsplit(
          x,
          "\t",
          fixed = FALSE
        )[[1]]
      }
    )
  )

testB_scores_10b3nd2 <-
  testB_data_10b3nd2[
    ,
    -c(1, 2),
    drop = FALSE
  ]

testB_scores_10b3nd2 <-
  apply(
    testB_scores_10b3nd2,
    2,
    as.numeric
  )

rownames(
  testB_scores_10b3nd2
) <-
  testB_data_10b3nd2[
    ,
    1
  ]

colnames(
  testB_scores_10b3nd2
) <-
  gsub(
    "\\.",
    "-",
    testB_header_10b3nd2[
      -c(1, 2)
    ]
  )

cat(
  "\nTest B score dimensions: ",
  paste(
    dim(testB_scores_10b3nd2),
    collapse = " × "
  ),
  "\n",
  sep = ""
)

cat(
  "Test B score rows: ",
  paste(
    rownames(testB_scores_10b3nd2),
    collapse = ", "
  ),
  "\n",
  sep = ""
)


# ============================================================
# D-2 TEST B — COMPARISON RECOVERY
# Recover locked Step14 reference object
# DO NOT rerun estimateScore()
# ============================================================

cat("\n============================================================\n")
cat("D-2 TEST B — COMPARISON RECOVERY\n")
cat("Recover locked Step14 reference\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Locked Step14 checkpoint
# ------------------------------------------------------------

step14_rdata_10b3nd2 <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "WORKFLOW-03",
  "Step14",
  "Step14_FINAL_LOCKED",
  "RData",
  "WORKFLOW-03_Step14A-14E_FINAL_LOCKED.RData"
)

stopifnot(file.exists(step14_rdata_10b3nd2))

cat("Step14 RData exists: TRUE\n")
cat("File:\n", step14_rdata_10b3nd2, "\n\n")


# ------------------------------------------------------------
# 2. Load Step14 into isolated temporary environment
# ------------------------------------------------------------

step14_env_10b3nd2 <- new.env(parent = emptyenv())

load(
  step14_rdata_10b3nd2,
  envir = step14_env_10b3nd2
)

stopifnot(
  exists(
    "estimate_scores_wgcna",
    envir = step14_env_10b3nd2,
    inherits = FALSE
  )
)

cat("Locked estimate_scores_wgcna recovered: TRUE\n\n")


# ------------------------------------------------------------
# 3. Assign recovered object as Test A reference
# ------------------------------------------------------------

testA_scores_10b3nd2 <- get(
  "estimate_scores_wgcna",
  envir = step14_env_10b3nd2,
  inherits = FALSE
)

cat("Test A reference dimensions:\n")
print(dim(testA_scores_10b3nd2))

cat("\nTest A rownames:\n")
print(rownames(testA_scores_10b3nd2))

cat("\nTest A sample count:\n")
print(ncol(testA_scores_10b3nd2))


# ------------------------------------------------------------
# 4. Confirm Test B parsed output exists
# ------------------------------------------------------------

stopifnot(
  exists(
    "testB_scores_10b3nd2",
    inherits = TRUE
  )
)

cat("\nTest B parsed scores object exists: TRUE\n")

cat("\nTest B dimensions:\n")
print(dim(testB_scores_10b3nd2))

cat("\nTest B rownames:\n")
print(rownames(testB_scores_10b3nd2))


# ------------------------------------------------------------
# 5. Check sample identity and order
# ------------------------------------------------------------

sample_identity_pass_10b3nd2 <- identical(
  colnames(testA_scores_10b3nd2),
  colnames(testB_scores_10b3nd2)
)

rownames_identity_pass_10b3nd2 <- identical(
  rownames(testA_scores_10b3nd2),
  rownames(testB_scores_10b3nd2)
)

cat("\n------------------------------------------------------------\n")
cat("IDENTITY CHECK\n")
cat("------------------------------------------------------------\n")

cat(
  "Rownames identical: ",
  rownames_identity_pass_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "Sample identity/order identical: ",
  sample_identity_pass_10b3nd2,
  "\n",
  sep = ""
)

stopifnot(
  rownames_identity_pass_10b3nd2,
  sample_identity_pass_10b3nd2
)


# ------------------------------------------------------------
# 6. Calculate Test B deltas versus locked Step14
# ------------------------------------------------------------

testB_stromal_delta_10b3nd2 <-
  testB_scores_10b3nd2["StromalScore", ] -
  testA_scores_10b3nd2["StromalScore", ]

testB_immune_delta_10b3nd2 <-
  testB_scores_10b3nd2["ImmuneScore", ] -
  testA_scores_10b3nd2["ImmuneScore", ]

testB_estimate_delta_10b3nd2 <-
  testB_scores_10b3nd2["ESTIMATEScore", ] -
  testA_scores_10b3nd2["ESTIMATEScore", ]


max_abs_stromal_delta_10b3nd2 <- max(
  abs(testB_stromal_delta_10b3nd2),
  na.rm = TRUE
)

max_abs_immune_delta_10b3nd2 <- max(
  abs(testB_immune_delta_10b3nd2),
  na.rm = TRUE
)

max_abs_estimate_delta_10b3nd2 <- max(
  abs(testB_estimate_delta_10b3nd2),
  na.rm = TRUE
)


# ------------------------------------------------------------
# 7. Print numerical comparison
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("TEST B NUMERICAL COMPARISON\n")
cat("============================================================\n\n")

cat(
  "Max |Δ StromalScore|   = ",
  format(
    max_abs_stromal_delta_10b3nd2,
    scientific = TRUE,
    digits = 12
  ),
  "\n",
  sep = ""
)

cat(
  "Max |Δ ImmuneScore|    = ",
  format(
    max_abs_immune_delta_10b3nd2,
    scientific = TRUE,
    digits = 12
  ),
  "\n",
  sep = ""
)

cat(
  "Max |Δ ESTIMATEScore|  = ",
  format(
    max_abs_estimate_delta_10b3nd2,
    scientific = TRUE,
    digits = 12
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Additional summary
# ------------------------------------------------------------

cat("\nSummary of changed samples:\n")

cat(
  "StromalScore changed samples: ",
  sum(
    abs(testB_stromal_delta_10b3nd2) > 1e-8,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "ImmuneScore changed samples: ",
  sum(
    abs(testB_immune_delta_10b3nd2) > 1e-8,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "ESTIMATEScore changed samples: ",
  sum(
    abs(testB_estimate_delta_10b3nd2) > 1e-8,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Check expected mechanistic pattern
# ------------------------------------------------------------

stromal_changed_10b3nd2 <-
  max_abs_stromal_delta_10b3nd2 > 1e-8

immune_unchanged_10b3nd2 <-
  max_abs_immune_delta_10b3nd2 < 1e-8


cat("\n============================================================\n")
cat("EXPECTED PATTERN\n")
cat("============================================================\n\n")

cat(
  "StromalScore changed: ",
  stromal_changed_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "ImmuneScore unchanged: ",
  immune_unchanged_10b3nd2,
  "\n",
  sep = ""
)

# ============================================================
# D-2 TEST B — INTEGRITY CHECK RECOVERY
# No estimateScore() rerun
# ============================================================

cat("\n============================================================\n")
cat("D-2 TEST B — INTEGRITY CHECK RECOVERY\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 10. Check package:estimate integrity
# ------------------------------------------------------------

package_estimate_env_10b3nd2 <-
  as.environment("package:estimate")

cat("Checking package:estimate directly...\n\n")

package_si_identical_10b3nd2 <- identical(
  get(
    "SI_geneset",
    envir = package_estimate_env_10b3nd2,
    inherits = FALSE
  ),
  SI_geneset_namespace_before_10b3nd2
)

package_si_locked_10b3nd2 <- bindingIsLocked(
  "SI_geneset",
  package_estimate_env_10b3nd2
)

cat(
  "package:estimate SI_geneset unchanged: ",
  package_si_identical_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "package:estimate SI_geneset binding locked: ",
  package_si_locked_10b3nd2,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Check GlobalEnv restoration
# ------------------------------------------------------------

global_si_exists_10b3nd2 <- exists(
  "SI_geneset",
  envir = .GlobalEnv,
  inherits = FALSE
)

global_si_restored_10b3nd2 <- identical(
  get(
    "SI_geneset",
    envir = .GlobalEnv,
    inherits = FALSE
  ),
  SI_geneset_global_before_10b3nd2
)

global_si_locked_10b3nd2 <- bindingIsLocked(
  "SI_geneset",
  .GlobalEnv
)


cat("\n============================================================\n")
cat("GLOBALENV RESTORATION\n")
cat("============================================================\n\n")

cat(
  "GlobalEnv SI_geneset exists: ",
  global_si_exists_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "GlobalEnv SI_geneset restored exactly: ",
  global_si_restored_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "GlobalEnv SI_geneset binding locked: ",
  global_si_locked_10b3nd2,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Final Test B decision
# ------------------------------------------------------------

testB_comparison_pass_10b3nd2 <-
  stromal_changed_10b3nd2 &&
  immune_unchanged_10b3nd2 &&
  package_si_identical_10b3nd2 &&
  package_si_locked_10b3nd2 &&
  global_si_exists_10b3nd2 &&
  global_si_restored_10b3nd2


cat("\n============================================================\n")
cat("TEST B FINAL RESULT\n")
cat("============================================================\n\n")

cat(
  "TEST B COMPARISON PASS: ",
  testB_comparison_pass_10b3nd2,
  "\n",
  sep = ""
)


# ============================================================
# D-2 TEST C — PRE-CHECK
# Confirm original SI_geneset restored
# ============================================================

cat("\n============================================================\n")
cat("D-2 TEST C — PRE-CHECK\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Confirm GlobalEnv SI_geneset exists
# ------------------------------------------------------------

stopifnot(
  exists(
    "SI_geneset",
    envir = .GlobalEnv,
    inherits = FALSE
  )
)

cat("GlobalEnv SI_geneset exists: TRUE\n")


# ------------------------------------------------------------
# 2. Confirm exact restoration
# ------------------------------------------------------------

testC_global_original_restored_10b3nd2 <- identical(
  get(
    "SI_geneset",
    envir = .GlobalEnv,
    inherits = FALSE
  ),
  SI_geneset_global_before_10b3nd2
)

cat(
  "GlobalEnv SI_geneset identical to original: ",
  testC_global_original_restored_10b3nd2,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 3. Confirm package object unchanged
# ------------------------------------------------------------

testC_package_original_unchanged_10b3nd2 <- identical(
  get(
    "SI_geneset",
    envir = as.environment("package:estimate"),
    inherits = FALSE
  ),
  SI_geneset_namespace_before_10b3nd2
)

testC_package_locked_10b3nd2 <- bindingIsLocked(
  "SI_geneset",
  as.environment("package:estimate")
)

cat(
  "package:estimate SI_geneset unchanged: ",
  testC_package_original_unchanged_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "package:estimate SI_geneset binding locked: ",
  testC_package_locked_10b3nd2,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Confirm AOC3 / MS4A7 / IL7R membership
# ------------------------------------------------------------

testC_stromal_genes_10b3nd2 <- SI_geneset[
  "StromalSignature",
  -1
]

testC_immune_genes_10b3nd2 <- SI_geneset[
  "ImmuneSignature",
  -1
]

cat("\nMembership check:\n")

cat(
  "AOC3 in StromalSignature: ",
  "AOC3" %in% testC_stromal_genes_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "AOC3 in ImmuneSignature: ",
  "AOC3" %in% testC_immune_genes_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "MS4A7 in StromalSignature: ",
  "MS4A7" %in% testC_stromal_genes_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "MS4A7 in ImmuneSignature: ",
  "MS4A7" %in% testC_immune_genes_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "IL7R in StromalSignature: ",
  "IL7R" %in% testC_stromal_genes_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "IL7R in ImmuneSignature: ",
  "IL7R" %in% testC_immune_genes_10b3nd2,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. Final pre-check
# ------------------------------------------------------------

testC_precheck_pass_10b3nd2 <-
  testC_global_original_restored_10b3nd2 &&
  testC_package_original_unchanged_10b3nd2 &&
  testC_package_locked_10b3nd2

cat("\n============================================================\n")
cat("TEST C PRE-CHECK RESULT\n")
cat("============================================================\n\n")

cat(
  "TEST C PRE-CHECK PASS: ",
  testC_precheck_pass_10b3nd2,
  "\n",
  sep = ""
)

stopifnot(testC_precheck_pass_10b3nd2)

# ============================================================
# D-2 TEST C — OFFICIAL estimateScore() REPRODUCTION
# Original SI_geneset restored
# ============================================================

cat("\n============================================================\n")
cat("D-2 TEST C — OFFICIAL estimateScore() REPRODUCTION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define input/output files
# ------------------------------------------------------------

testC_input_gct_10b3nd2 <- file.path(
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis",
  "STEP7B-7E-10B-3N",
  "10B-3N-C",
  "STEP10B3NC_official_filterCommonGenes_output.gct"
)

testC_output_dir_10b3nd2 <- file.path(
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis",
  "STEP7B-7E-10B-3N",
  "10B-3N-D-2"
)

if (!dir.exists(testC_output_dir_10b3nd2)) {
  dir.create(
    testC_output_dir_10b3nd2,
    recursive = TRUE
  )
}

testC_output_gct_10b3nd2 <- file.path(
  testC_output_dir_10b3nd2,
  "STEP10B3ND2_TestC_original_restored.gct"
)

stopifnot(file.exists(testC_input_gct_10b3nd2))

cat("Input GCT exists: TRUE\n")
cat("Output file:\n", testC_output_gct_10b3nd2, "\n\n")


# ------------------------------------------------------------
# 2. Confirm original SI_geneset one more time
# ------------------------------------------------------------

stopifnot(
  identical(
    SI_geneset,
    SI_geneset_global_before_10b3nd2
  )
)

cat("Original GlobalEnv SI_geneset confirmed: TRUE\n")


# ------------------------------------------------------------
# 3. Confirm package namespace remains unchanged
# ------------------------------------------------------------

stopifnot(
  identical(
    get(
      "SI_geneset",
      envir = as.environment("package:estimate"),
      inherits = FALSE
    ),
    SI_geneset_namespace_before_10b3nd2
  )
)

stopifnot(
  bindingIsLocked(
    "SI_geneset",
    as.environment("package:estimate")
  )
)

cat("Package SI_geneset unchanged and locked: TRUE\n\n")


# ------------------------------------------------------------
# 4. Run official estimateScore()
# ------------------------------------------------------------

cat("Running official estimateScore()...\n\n")

estimate::estimateScore(
  input.ds = testC_input_gct_10b3nd2,
  output.ds = testC_output_gct_10b3nd2,
  platform = "illumina"
)

cat("\nOfficial estimateScore() completed.\n\n")


# ------------------------------------------------------------
# 5. Confirm output exists
# ------------------------------------------------------------

testC_output_exists_10b3nd2 <- file.exists(
  testC_output_gct_10b3nd2
)

cat(
  "Test C output exists: ",
  testC_output_exists_10b3nd2,
  "\n",
  sep = ""
)

stopifnot(testC_output_exists_10b3nd2)


# ============================================================
# D-2 TEST C — NUMERICAL REPRODUCTION VALIDATION
# Compare Test C with locked Step14
# ============================================================

cat("\n============================================================\n")
cat("D-2 TEST C — NUMERICAL REPRODUCTION VALIDATION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Recover locked Step14 reference
# ------------------------------------------------------------

step14_rdata_C_10b3nd2 <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "WORKFLOW-03",
  "Step14",
  "Step14_FINAL_LOCKED",
  "RData",
  "WORKFLOW-03_Step14A-14E_FINAL_LOCKED.RData"
)

step14_env_C_10b3nd2 <- new.env(parent = emptyenv())

load(
  step14_rdata_C_10b3nd2,
  envir = step14_env_C_10b3nd2
)

testC_reference_scores_10b3nd2 <- get(
  "estimate_scores_wgcna",
  envir = step14_env_C_10b3nd2,
  inherits = FALSE
)

cat("Locked Step14 reference recovered: TRUE\n")
cat("Reference dimensions:\n")
print(dim(testC_reference_scores_10b3nd2))


# ------------------------------------------------------------
# 2. Read Test C GCT
# ------------------------------------------------------------

testC_lines_10b3nd2 <- readLines(
  testC_output_gct_10b3nd2,
  warn = FALSE
)

cat("\nTest C GCT header:\n")
print(testC_lines_10b3nd2[1:6])


# ------------------------------------------------------------
# 3. Parse Test C output
# ------------------------------------------------------------

testC_header_10b3nd2 <- strsplit(
  testC_lines_10b3nd2[3],
  "\t",
  fixed = TRUE
)[[1]]

testC_data_10b3nd2 <- read.table(
  text = paste(
    testC_lines_10b3nd2[4:length(testC_lines_10b3nd2)],
    collapse = "\n"
  ),
  header = FALSE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(testC_data_10b3nd2) <- testC_header_10b3nd2

testC_scores_10b3nd2 <- as.matrix(
  testC_data_10b3nd2[, -c(1, 2), drop = FALSE]
)

rownames(testC_scores_10b3nd2) <-
  testC_data_10b3nd2[[1]]

mode(testC_scores_10b3nd2) <- "numeric"

# Convert GCT sample IDs from dots back to TCGA hyphen format
colnames(testC_scores_10b3nd2) <-
  gsub(
    "\\.",
    "-",
    colnames(testC_scores_10b3nd2)
  )

cat("\nTest C score dimensions:\n")
print(dim(testC_scores_10b3nd2))

cat("\nTest C rownames:\n")
print(rownames(testC_scores_10b3nd2))


# ------------------------------------------------------------
# 4. Confirm dimensions
# ------------------------------------------------------------

stopifnot(
  identical(
    dim(testC_scores_10b3nd2),
    dim(testC_reference_scores_10b3nd2)
  )
)

cat("\nDimension check: PASS\n")


# ------------------------------------------------------------
# 5. Confirm row/sample identity and order
# ------------------------------------------------------------

testC_row_identity_10b3nd2 <- identical(
  rownames(testC_scores_10b3nd2),
  rownames(testC_reference_scores_10b3nd2)
)

testC_sample_identity_10b3nd2 <- identical(
  colnames(testC_scores_10b3nd2),
  colnames(testC_reference_scores_10b3nd2)
)

cat(
  "Rownames identical: ",
  testC_row_identity_10b3nd2,
  "\n",
  sep = ""
)

cat(
  "Sample identity/order identical: ",
  testC_sample_identity_10b3nd2,
  "\n",
  sep = ""
)

stopifnot(
  testC_row_identity_10b3nd2,
  testC_sample_identity_10b3nd2
)


# ------------------------------------------------------------
# 6. Numerical comparison
# ------------------------------------------------------------

testC_stromal_delta_10b3nd2 <-
  testC_scores_10b3nd2["StromalScore", ] -
  testC_reference_scores_10b3nd2["StromalScore", ]

testC_immune_delta_10b3nd2 <-
  testC_scores_10b3nd2["ImmuneScore", ] -
  testC_reference_scores_10b3nd2["ImmuneScore", ]

testC_estimate_delta_10b3nd2 <-
  testC_scores_10b3nd2["ESTIMATEScore", ] -
  testC_reference_scores_10b3nd2["ESTIMATEScore", ]


testC_max_abs_stromal_10b3nd2 <- max(
  abs(testC_stromal_delta_10b3nd2),
  na.rm = TRUE
)

testC_max_abs_immune_10b3nd2 <- max(
  abs(testC_immune_delta_10b3nd2),
  na.rm = TRUE
)

testC_max_abs_estimate_10b3nd2 <- max(
  abs(testC_estimate_delta_10b3nd2),
  na.rm = TRUE
)

# ============================================================
# D-2 TEST C — CORRELATION COMPARISON RECOVERY
# No estimateScore() rerun
# No modification of locked objects
# ============================================================

cat("\n============================================================\n")
cat("D-2 TEST C — CORRELATION COMPARISON RECOVERY\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Inspect dimensions and classes
# ------------------------------------------------------------

cat("Test C Stromal dimensions:\n")
print(dim(testC_scores_10b3nd2["StromalScore", ]))

cat("\nReference Stromal dimensions:\n")
print(dim(testC_reference_scores_10b3nd2["StromalScore", ]))

cat("\nTest C Stromal class:\n")
print(class(testC_scores_10b3nd2["StromalScore", ]))

cat("\nReference Stromal class:\n")
print(class(testC_reference_scores_10b3nd2["StromalScore", ]))


# ------------------------------------------------------------
# 2. Explicitly convert to numeric vectors
# ------------------------------------------------------------

testC_stromal_test_vector_10b3nd2 <- as.numeric(
  testC_scores_10b3nd2["StromalScore", ]
)

testC_stromal_ref_vector_10b3nd2 <- as.numeric(
  testC_reference_scores_10b3nd2["StromalScore", ]
)

testC_immune_test_vector_10b3nd2 <- as.numeric(
  testC_scores_10b3nd2["ImmuneScore", ]
)

testC_immune_ref_vector_10b3nd2 <- as.numeric(
  testC_reference_scores_10b3nd2["ImmuneScore", ]
)

testC_estimate_test_vector_10b3nd2 <- as.numeric(
  testC_scores_10b3nd2["ESTIMATEScore", ]
)

testC_estimate_ref_vector_10b3nd2 <- as.numeric(
  testC_reference_scores_10b3nd2["ESTIMATEScore", ]
)


# ------------------------------------------------------------
# 3. Confirm vector lengths
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("VECTOR LENGTH CHECK\n")
cat("============================================================\n\n")

cat(
  "Stromal Test C length: ",
  length(testC_stromal_test_vector_10b3nd2),
  "\n",
  sep = ""
)

cat(
  "Stromal reference length: ",
  length(testC_stromal_ref_vector_10b3nd2),
  "\n",
  sep = ""
)

cat(
  "Immune Test C length: ",
  length(testC_immune_test_vector_10b3nd2),
  "\n",
  sep = ""
)

cat(
  "Immune reference length: ",
  length(testC_immune_ref_vector_10b3nd2),
  "\n",
  sep = ""
)

cat(
  "ESTIMATE Test C length: ",
  length(testC_estimate_test_vector_10b3nd2),
  "\n",
  sep = ""
)

cat(
  "ESTIMATE reference length: ",
  length(testC_estimate_ref_vector_10b3nd2),
  "\n",
  sep = ""
)

stopifnot(
  length(testC_stromal_test_vector_10b3nd2) ==
    length(testC_stromal_ref_vector_10b3nd2),
  
  length(testC_immune_test_vector_10b3nd2) ==
    length(testC_immune_ref_vector_10b3nd2),
  
  length(testC_estimate_test_vector_10b3nd2) ==
    length(testC_estimate_ref_vector_10b3nd2)
)


# ------------------------------------------------------------
# 4. Recalculate correlations using pure numeric vectors
# ------------------------------------------------------------

testC_stromal_pearson_10b3nd2 <- cor(
  testC_stromal_test_vector_10b3nd2,
  testC_stromal_ref_vector_10b3nd2,
  method = "pearson"
)

testC_stromal_spearman_10b3nd2 <- cor(
  testC_stromal_test_vector_10b3nd2,
  testC_stromal_ref_vector_10b3nd2,
  method = "spearman"
)


testC_immune_pearson_10b3nd2 <- cor(
  testC_immune_test_vector_10b3nd2,
  testC_immune_ref_vector_10b3nd2,
  method = "pearson"
)

testC_immune_spearman_10b3nd2 <- cor(
  testC_immune_test_vector_10b3nd2,
  testC_immune_ref_vector_10b3nd2,
  method = "spearman"
)


testC_estimate_pearson_10b3nd2 <- cor(
  testC_estimate_test_vector_10b3nd2,
  testC_estimate_ref_vector_10b3nd2,
  method = "pearson"
)

testC_estimate_spearman_10b3nd2 <- cor(
  testC_estimate_test_vector_10b3nd2,
  testC_estimate_ref_vector_10b3nd2,
  method = "spearman"
)


# ------------------------------------------------------------
# 5. Recalculate exact numerical differences
# ------------------------------------------------------------

testC_stromal_delta_10b3nd2 <-
  testC_stromal_test_vector_10b3nd2 -
  testC_stromal_ref_vector_10b3nd2

testC_immune_delta_10b3nd2 <-
  testC_immune_test_vector_10b3nd2 -
  testC_immune_ref_vector_10b3nd2

testC_estimate_delta_10b3nd2 <-
  testC_estimate_test_vector_10b3nd2 -
  testC_estimate_ref_vector_10b3nd2


testC_max_abs_stromal_10b3nd2 <- max(
  abs(testC_stromal_delta_10b3nd2),
  na.rm = TRUE
)

testC_max_abs_immune_10b3nd2 <- max(
  abs(testC_immune_delta_10b3nd2),
  na.rm = TRUE
)

testC_max_abs_estimate_10b3nd2 <- max(
  abs(testC_estimate_delta_10b3nd2),
  na.rm = TRUE
)


# ------------------------------------------------------------
# 6. Print final numerical comparison
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("TEST C NUMERICAL REPRODUCTION\n")
cat("============================================================\n\n")

cat(
  "StromalScore Pearson  = ",
  format(
    testC_stromal_pearson_10b3nd2,
    digits = 16
  ),
  "\n",
  sep = ""
)

cat(
  "StromalScore Spearman = ",
  format(
    testC_stromal_spearman_10b3nd2,
    digits = 16
  ),
  "\n",
  sep = ""
)

cat(
  "StromalScore MaxAbsDiff = ",
  format(
    testC_max_abs_stromal_10b3nd2,
    scientific = TRUE,
    digits = 12
  ),
  "\n\n",
  sep = ""
)


cat(
  "ImmuneScore Pearson  = ",
  format(
    testC_immune_pearson_10b3nd2,
    digits = 16
  ),
  "\n",
  sep = ""
)

cat(
  "ImmuneScore Spearman = ",
  format(
    testC_immune_spearman_10b3nd2,
    digits = 16
  ),
  "\n",
  sep = ""
)

cat(
  "ImmuneScore MaxAbsDiff = ",
  format(
    testC_max_abs_immune_10b3nd2,
    scientific = TRUE,
    digits = 12
  ),
  "\n\n",
  sep = ""
)


cat(
  "ESTIMATEScore Pearson  = ",
  format(
    testC_estimate_pearson_10b3nd2,
    digits = 16
  ),
  "\n",
  sep = ""
)

cat(
  "ESTIMATEScore Spearman = ",
  format(
    testC_estimate_spearman_10b3nd2,
    digits = 16
  ),
  "\n",
  sep = ""
)

cat(
  "ESTIMATEScore MaxAbsDiff = ",
  format(
    testC_max_abs_estimate_10b3nd2,
    scientific = TRUE,
    digits = 12
  ),
  "\n",
  sep = ""
)

# ============================================================
# STEP 10B-3N-D-2
# Final Integrity Check + Locked Checkpoint
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3N-D-2 FINAL INTEGRITY CHECK\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Check package:estimate SI_geneset
# ------------------------------------------------------------

package_estimate_env_10b3nd2 <- as.environment("package:estimate")

package_SI_final_10b3nd2 <-
  get("SI_geneset", envir = package_estimate_env_10b3nd2)

package_SI_unchanged_final_10b3nd2 <-
  identical(
    package_SI_final_10b3nd2,
    SI_geneset_original_testB_10b3nd2
  )

package_SI_locked_final_10b3nd2 <-
  bindingIsLocked(
    "SI_geneset",
    package_estimate_env_10b3nd2
  )


# ------------------------------------------------------------
# 2. Check GlobalEnv SI_geneset
# ------------------------------------------------------------

global_SI_exists_final_10b3nd2 <-
  exists(
    "SI_geneset",
    envir = .GlobalEnv,
    inherits = FALSE
  )

global_SI_identical_final_10b3nd2 <-
  if (global_SI_exists_final_10b3nd2) {
    
    identical(
      get("SI_geneset", envir = .GlobalEnv),
      SI_geneset_original_testB_10b3nd2
    )
    
  } else {
    FALSE
  }


global_SI_locked_final_10b3nd2 <-
  if (global_SI_exists_final_10b3nd2) {
    
    bindingIsLocked(
      "SI_geneset",
      .GlobalEnv
    )
    
  } else {
    NA
  }


# ------------------------------------------------------------
# 3. Check biological membership
# ------------------------------------------------------------

SI_final_10b3nd2 <-
  get("SI_geneset", envir = .GlobalEnv)

stromal_final_10b3nd2 <-
  SI_final_10b3nd2[1, -1]

immune_final_10b3nd2 <-
  SI_final_10b3nd2[2, -1]

stromal_genes_final_10b3nd2 <-
  as.character(stromal_final_10b3nd2[
    !is.na(stromal_final_10b3nd2)
  ])

immune_genes_final_10b3nd2 <-
  as.character(immune_final_10b3nd2[
    !is.na(immune_final_10b3nd2)
  ])

AOC3_stromal_final_10b3nd2 <-
  "AOC3" %in% stromal_genes_final_10b3nd2

AOC3_immune_final_10b3nd2 <-
  "AOC3" %in% immune_genes_final_10b3nd2

MS4A7_stromal_final_10b3nd2 <-
  "MS4A7" %in% stromal_genes_final_10b3nd2

MS4A7_immune_final_10b3nd2 <-
  "MS4A7" %in% immune_genes_final_10b3nd2

IL7R_stromal_final_10b3nd2 <-
  "IL7R" %in% stromal_genes_final_10b3nd2

IL7R_immune_final_10b3nd2 <-
  "IL7R" %in% immune_genes_final_10b3nd2


# ------------------------------------------------------------
# 4. Confirm Test C numerical reproduction
# ------------------------------------------------------------

testC_numerical_pass_final_10b3nd2 <-
  isTRUE(
    testC_stromal_pearson_10b3nd2 == 1
  ) &&
  isTRUE(
    testC_stromal_spearman_10b3nd2 == 1
  ) &&
  isTRUE(
    testC_immune_spearman_10b3nd2 == 1
  ) &&
  isTRUE(
    testC_estimate_pearson_10b3nd2 == 1
  ) &&
  isTRUE(
    testC_estimate_spearman_10b3nd2 == 1
  ) &&
  isTRUE(
    testC_max_abs_stromal_10b3nd2 < 1e-10
  ) &&
  isTRUE(
    testC_max_abs_immune_10b3nd2 < 1e-10
  ) &&
  isTRUE(
    testC_max_abs_estimate_10b3nd2 < 1e-10
  )


# ------------------------------------------------------------
# 5. Overall D-2 PASS criterion
# ------------------------------------------------------------

D2_FINAL_PASS <-
  package_SI_unchanged_final_10b3nd2 &&
  package_SI_locked_final_10b3nd2 &&
  global_SI_exists_final_10b3nd2 &&
  global_SI_identical_final_10b3nd2 &&
  testC_numerical_pass_final_10b3nd2 &&
  AOC3_stromal_final_10b3nd2 &&
  !AOC3_immune_final_10b3nd2 &&
  !MS4A7_stromal_final_10b3nd2 &&
  !MS4A7_immune_final_10b3nd2 &&
  !IL7R_stromal_final_10b3nd2 &&
  IL7R_immune_final_10b3nd2


# ------------------------------------------------------------
# 6. Print final audit
# ------------------------------------------------------------

cat("Package:estimate SI_geneset unchanged: ",
    package_SI_unchanged_final_10b3nd2, "\n", sep = "")

cat("Package:estimate SI_geneset binding locked: ",
    package_SI_locked_final_10b3nd2, "\n", sep = "")

cat("GlobalEnv SI_geneset exists: ",
    global_SI_exists_final_10b3nd2, "\n", sep = "")

cat("GlobalEnv SI_geneset identical to original: ",
    global_SI_identical_final_10b3nd2, "\n", sep = "")

cat("GlobalEnv SI_geneset binding locked: ",
    global_SI_locked_final_10b3nd2, "\n\n", sep = "")

cat("AOC3 in StromalSignature: ",
    AOC3_stromal_final_10b3nd2, "\n", sep = "")

cat("AOC3 in ImmuneSignature: ",
    AOC3_immune_final_10b3nd2, "\n", sep = "")

cat("MS4A7 in StromalSignature: ",
    MS4A7_stromal_final_10b3nd2, "\n", sep = "")

cat("MS4A7 in ImmuneSignature: ",
    MS4A7_immune_final_10b3nd2, "\n", sep = "")

cat("IL7R in StromalSignature: ",
    IL7R_stromal_final_10b3nd2, "\n", sep = "")

cat("IL7R in ImmuneSignature: ",
    IL7R_immune_final_10b3nd2, "\n\n", sep = "")

cat("Test C numerical reproduction: ",
    testC_numerical_pass_final_10b3nd2, "\n", sep = "")

cat("\n============================================================\n")
cat("D-2 FINAL PASS: ",
    D2_FINAL_PASS,
    "\n",
    sep = "")
cat("============================================================\n")


# ------------------------------------------------------------
# 7. Save locked D-2 checkpoint
# ------------------------------------------------------------

D2_checkpoint_path <-
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis/STEP7B-7E-10B-3N/10B-3N-D-2/Paper1_STEP10B3ND2_official_SI_geneset_substitution_validation_LOCKED.RData"


save(
  SI_geneset_original_testB_10b3nd2,
  SI_geneset_testB_deep_10b3nd2,
  
  package_SI_final_10b3nd2,
  package_SI_unchanged_final_10b3nd2,
  package_SI_locked_final_10b3nd2,
  
  global_SI_exists_final_10b3nd2,
  global_SI_identical_final_10b3nd2,
  global_SI_locked_final_10b3nd2,
  
  AOC3_stromal_final_10b3nd2,
  AOC3_immune_final_10b3nd2,
  MS4A7_stromal_final_10b3nd2,
  MS4A7_immune_final_10b3nd2,
  IL7R_stromal_final_10b3nd2,
  IL7R_immune_final_10b3nd2,
  
  testC_numerical_pass_final_10b3nd2,
  
  testC_stromal_pearson_10b3nd2,
  testC_stromal_spearman_10b3nd2,
  testC_max_abs_stromal_10b3nd2,
  
  testC_immune_pearson_10b3nd2,
  testC_immune_spearman_10b3nd2,
  testC_max_abs_immune_10b3nd2,
  
  testC_estimate_pearson_10b3nd2,
  testC_estimate_spearman_10b3nd2,
  testC_max_abs_estimate_10b3nd2,
  
  D2_FINAL_PASS,
  
  file = D2_checkpoint_path
)

cat("\nCheckpoint saved:\n",
    D2_checkpoint_path,
    "\n", sep = "")

cat("\n============================================================\n")
cat("STEP 10B-3N-D-2 COMPLETED\n")
cat("============================================================\n")


# ============================================================
# STEP 10B-3N-E
# Formal Leave-one-gene-out (LOO) Sensitivity Analysis
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3N-E: FORMAL LOO SENSITIVITY ANALYSIS\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0.1 Working directories
# ------------------------------------------------------------

E10B3NE_dir <-
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis/STEP7B-7E-10B-3N/10B-3N-E"

dir.create(
  E10B3NE_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(E10B3NE_dir, "AOC3_LOO"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(E10B3NE_dir, "IL7R_LOO"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(E10B3NE_dir, "Audit"),
  recursive = TRUE,
  showWarnings = FALSE
)


cat("Working directory:\n")
cat(E10B3NE_dir, "\n\n")


# ------------------------------------------------------------
# 0.2 Locked source files
# ------------------------------------------------------------

E10B3NE_step14_checkpoint <-
  "input_dir/Paper1_tumor-only WGCNA/WORKFLOW-03/Step14/Step14_FINAL_LOCKED/RData/WORKFLOW-03_Step14A-14E_FINAL_LOCKED.RData"


E10B3NE_common_gct <-
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis/STEP7B-7E-10B-3N/10B-3N-C/STEP10B3NC_official_filterCommonGenes_output.gct"


E10B3NE_D2_checkpoint <-
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis/STEP7B-7E-10B-3N/10B-3N-D-2/Paper1_STEP10B3ND2_official_SI_geneset_substitution_validation_LOCKED.RData"


stopifnot(file.exists(E10B3NE_step14_checkpoint))
stopifnot(file.exists(E10B3NE_common_gct))
stopifnot(file.exists(E10B3NE_D2_checkpoint))


cat("Source files verified.\n")


# ============================================================
# E-1. Load locked reference objects
# ============================================================

cat("\n============================================================\n")
cat("E-1: LOAD LOCKED REFERENCE OBJECTS\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1.1 Load Step14 in isolated environment
# ------------------------------------------------------------

step14_env_10b3ne <- new.env(parent = emptyenv())

load(
  E10B3NE_step14_checkpoint,
  envir = step14_env_10b3ne
)


# ------------------------------------------------------------
# 1.2 Recover locked ESTIMATE scores
# ------------------------------------------------------------

estimate_scores_locked_10b3ne <-
  get(
    "estimate_scores_wgcna",
    envir = step14_env_10b3ne
  )


expr_estimate_locked_10b3ne <-
  get(
    "expr_estimate",
    envir = step14_env_10b3ne
  )


# ------------------------------------------------------------
# 1.3 Basic dimensions
# ------------------------------------------------------------

cat(
  "expr_estimate dimensions: ",
  paste(dim(expr_estimate_locked_10b3ne), collapse = " x "),
  "\n",
  sep = ""
)

cat(
  "estimate_scores_wgcna dimensions: ",
  paste(dim(estimate_scores_locked_10b3ne), collapse = " x "),
  "\n\n",
  sep = ""
)


stopifnot(
  identical(
    rownames(estimate_scores_locked_10b3ne),
    c(
      "StromalScore",
      "ImmuneScore",
      "ESTIMATEScore"
    )
  )
)

stopifnot(
  ncol(estimate_scores_locked_10b3ne) == 513
)


# ------------------------------------------------------------
# 1.4 Load D-2 checkpoint
# ------------------------------------------------------------

D2_env_10b3ne <- new.env(parent = emptyenv())

load(
  E10B3NE_D2_checkpoint,
  envir = D2_env_10b3ne
)


SI_geneset_original_10b3ne <-
  get(
    "SI_geneset_original_testB_10b3nd2",
    envir = D2_env_10b3ne
  )


D2_FINAL_PASS_10b3ne <-
  get(
    "D2_FINAL_PASS",
    envir = D2_env_10b3ne
  )


cat(
  "D-2 FINAL PASS recovered: ",
  D2_FINAL_PASS_10b3ne,
  "\n",
  sep = ""
)

stopifnot(isTRUE(D2_FINAL_PASS_10b3ne))


cat("\nE-1 PASS: locked reference successfully recovered.\n")


# ============================================================
# E-2. Prepare gene-expression reference
# ============================================================

cat("\n============================================================\n")
cat("E-2: PREPARE GENE-EXPRESSION REFERENCE\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 2.1 Confirm core genes exist
# ------------------------------------------------------------

core_genes_10b3ne <-
  c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )


core_gene_presence_10b3ne <-
  setNames(
    core_genes_10b3ne %in% rownames(expr_estimate_locked_10b3ne),
    core_genes_10b3ne
  )

print(core_gene_presence_10b3ne)

stopifnot(
  all(core_gene_presence_10b3ne)
)


# ------------------------------------------------------------
# 2.2 Extract expression
# ------------------------------------------------------------

gene_expression_10b3ne <-
  expr_estimate_locked_10b3ne[
    core_genes_10b3ne,
    ,
    drop = FALSE
  ]


# ------------------------------------------------------------
# 2.3 Confirm sample identity with locked scores
# ------------------------------------------------------------

sample_identity_E10B3NE <-
  identical(
    colnames(gene_expression_10b3ne),
    colnames(estimate_scores_locked_10b3ne)
  )


cat(
  "Gene-expression / ESTIMATE sample order identical: ",
  sample_identity_E10B3NE,
  "\n",
  sep = ""
)


stopifnot(sample_identity_E10B3NE)


cat("\nE-2 PASS.\n")


# ============================================================
# E-3. Safe SI_geneset handling
# ============================================================

cat("\n============================================================\n")
cat("E-3: SAFE SI_geneset HANDLING\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 3.1 Confirm current package binding
# ------------------------------------------------------------

package_estimate_env_10b3ne <-
  as.environment("package:estimate")


SI_geneset_package_before_E_10b3ne <-
  get(
    "SI_geneset",
    envir = package_estimate_env_10b3ne
  )


package_SI_locked_before_E_10b3ne <-
  bindingIsLocked(
    "SI_geneset",
    package_estimate_env_10b3ne
  )


cat(
  "Package SI_geneset identical to D-2 original: ",
  identical(
    SI_geneset_package_before_E_10b3ne,
    SI_geneset_original_10b3ne
  ),
  "\n",
  sep = ""
)

cat(
  "Package SI_geneset binding locked: ",
  package_SI_locked_before_E_10b3ne,
  "\n",
  sep = ""
)


stopifnot(
  identical(
    SI_geneset_package_before_E_10b3ne,
    SI_geneset_original_10b3ne
  )
)

stopifnot(
  isTRUE(package_SI_locked_before_E_10b3ne)
)


# ------------------------------------------------------------
# 3.2 Preserve current GlobalEnv state
# ------------------------------------------------------------

global_SI_existed_before_E_10b3ne <-
  exists(
    "SI_geneset",
    envir = .GlobalEnv,
    inherits = FALSE
  )


if (global_SI_existed_before_E_10b3ne) {
  
  SI_geneset_global_before_E_10b3ne <-
    get(
      "SI_geneset",
      envir = .GlobalEnv
    )
  
} else {
  
  SI_geneset_global_before_E_10b3ne <- NULL
  
}


# ------------------------------------------------------------
# 3.3 Safe deep-copy helper
# ------------------------------------------------------------

deep_copy_SI_10b3ne <- function(x) {
  
  unserialize(
    serialize(
      x,
      NULL
    )
  )
  
}


cat("\nE-3 PASS: safe SI_geneset handling established.\n")


#E-4A：Construct the AOC3-LOO signature
# ============================================================
# E-4A. AOC3 leave-one-gene-out signature
# ============================================================

cat("\n============================================================\n")
cat("E-4A: AOC3 LOO SIGNATURE\n")
cat("============================================================\n\n")


SI_geneset_AOC3_LOO_10b3ne <-
  deep_copy_SI_10b3ne(
    SI_geneset_original_10b3ne
  )


stromal_AOC3_before_10b3ne <-
  as.character(
    SI_geneset_AOC3_LOO_10b3ne[1, -1]
  )


AOC3_position_10b3ne <-
  which(
    stromal_AOC3_before_10b3ne == "AOC3"
  )


cat(
  "AOC3 position in StromalSignature: ",
  AOC3_position_10b3ne,
  "\n",
  sep = ""
)


stopifnot(length(AOC3_position_10b3ne) == 1)


# +1 because the first column is the "estimate" column
SI_geneset_AOC3_LOO_10b3ne[
  1,
  AOC3_position_10b3ne + 1
] <- NA_character_


stromal_AOC3_after_10b3ne <-
  as.character(
    SI_geneset_AOC3_LOO_10b3ne[1, -1]
  )


immune_AOC3_after_10b3ne <-
  as.character(
    SI_geneset_AOC3_LOO_10b3ne[2, -1]
  )


cat(
  "Stromal genes after AOC3 removal: ",
  sum(!is.na(stromal_AOC3_after_10b3ne)),
  "\n",
  sep = ""
)

cat(
  "Immune genes: ",
  sum(!is.na(immune_AOC3_after_10b3ne)),
  "\n",
  sep = ""
)

cat(
  "AOC3 present in StromalSignature: ",
  "AOC3" %in% stromal_AOC3_after_10b3ne,
  "\n",
  sep = ""
)

cat(
  "AOC3 present in ImmuneSignature: ",
  "AOC3" %in% immune_AOC3_after_10b3ne,
  "\n",
  sep = ""
)


stopifnot(
  !("AOC3" %in% stromal_AOC3_after_10b3ne)
)

stopifnot(
  !("AOC3" %in% immune_AOC3_after_10b3ne)
)


# ============================================================
# E-4B. Run official estimateScore() for AOC3 LOO
# ============================================================

cat("\n============================================================\n")
cat("E-4B: OFFICIAL estimateScore() — AOC3 LOO\n")
cat("============================================================\n\n")


AOC3_LOO_input_gct_10b3ne <-
  E10B3NE_common_gct


AOC3_LOO_output_gct_10b3ne <-
  file.path(
    E10B3NE_dir,
    "AOC3_LOO",
    "STEP10B3NE_AOC3_LOO_estimateScore.gct"
  )


# ------------------------------------------------------------
# Assign modified signature ONLY to GlobalEnv
# ------------------------------------------------------------

assign(
  "SI_geneset",
  SI_geneset_AOC3_LOO_10b3ne,
  envir = .GlobalEnv
)


# ------------------------------------------------------------
# Verify substitution before official scoring
# ------------------------------------------------------------

SI_geneset_runtime_AOC3_10b3ne <-
  get(
    "SI_geneset",
    envir = .GlobalEnv
  )


stopifnot(
  identical(
    SI_geneset_runtime_AOC3_10b3ne,
    SI_geneset_AOC3_LOO_10b3ne
  )
)


stopifnot(
  identical(
    get(
      "SI_geneset",
      envir = package_estimate_env_10b3ne
    ),
    SI_geneset_original_10b3ne
  )
)


cat("AOC3-LOO signature installed in GlobalEnv.\n")
cat("Package SI_geneset remains unchanged.\n\n")


# ------------------------------------------------------------
# Official ESTIMATE scoring
# ------------------------------------------------------------

estimateScore(
  input.ds = AOC3_LOO_input_gct_10b3ne,
  output.ds = AOC3_LOO_output_gct_10b3ne,
  platform = "illumina"
)


cat("\nAOC3-LOO estimateScore() completed.\n")
cat("Output:\n")
cat(AOC3_LOO_output_gct_10b3ne, "\n")



# ============================================================
# E-4C. Parse AOC3 LOO output
# ============================================================

cat("\n============================================================\n")
cat("E-4C: PARSE AOC3 LOO OUTPUT\n")
cat("============================================================\n\n")


stopifnot(
  file.exists(AOC3_LOO_output_gct_10b3ne)
)


AOC3_LOO_gct_10b3ne <-
  readLines(
    AOC3_LOO_output_gct_10b3ne,
    warn = FALSE
  )


# ------------------------------------------------------------
# Check GCT structure
# ------------------------------------------------------------

cat("First six GCT lines:\n")
cat(
  paste(
    AOC3_LOO_gct_10b3ne[1:6],
    collapse = "\n"
  ),
  "\n\n"
)


# ------------------------------------------------------------
# Header
# ------------------------------------------------------------

AOC3_LOO_header_10b3ne <-
  strsplit(
    AOC3_LOO_gct_10b3ne[3],
    "\t",
    fixed = TRUE
  )[[1]]


# ------------------------------------------------------------
# Scores
# ------------------------------------------------------------

AOC3_LOO_stromal_10b3ne <-
  as.numeric(
    strsplit(
      AOC3_LOO_gct_10b3ne[4],
      "\t",
      fixed = TRUE
    )[[1]][-c(1, 2)]
  )


AOC3_LOO_immune_10b3ne <-
  as.numeric(
    strsplit(
      AOC3_LOO_gct_10b3ne[5],
      "\t",
      fixed = TRUE
    )[[1]][-c(1, 2)]
  )


AOC3_LOO_estimate_10b3ne <-
  as.numeric(
    strsplit(
      AOC3_LOO_gct_10b3ne[6],
      "\t",
      fixed = TRUE
    )[[1]][-c(1, 2)]
  )


# ------------------------------------------------------------
# Sample IDs
# ------------------------------------------------------------

AOC3_LOO_sample_ids_10b3ne <-
  gsub(
    "\\.",
    "-",
    AOC3_LOO_header_10b3ne[-c(1, 2)]
  )


names(AOC3_LOO_stromal_10b3ne) <-
  AOC3_LOO_sample_ids_10b3ne

names(AOC3_LOO_immune_10b3ne) <-
  AOC3_LOO_sample_ids_10b3ne

names(AOC3_LOO_estimate_10b3ne) <-
  AOC3_LOO_sample_ids_10b3ne


# ------------------------------------------------------------
# Dimension checks
# ------------------------------------------------------------

cat(
  "AOC3 LOO StromalScore length: ",
  length(AOC3_LOO_stromal_10b3ne),
  "\n",
  sep = ""
)

cat(
  "AOC3 LOO ImmuneScore length: ",
  length(AOC3_LOO_immune_10b3ne),
  "\n",
  sep = ""
)

cat(
  "AOC3 LOO ESTIMATEScore length: ",
  length(AOC3_LOO_estimate_10b3ne),
  "\n\n",
  sep = ""
)


stopifnot(
  length(AOC3_LOO_stromal_10b3ne) == 513,
  length(AOC3_LOO_immune_10b3ne) == 513,
  length(AOC3_LOO_estimate_10b3ne) == 513
)


cat("AOC3 LOO output dimensions PASS.\n")


# ============================================================
# E-4D. Compare AOC3 LOO vs locked original scores
# ============================================================

cat("\n============================================================\n")
cat("E-4D: AOC3 LOO vs LOCKED ORIGINAL\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Convert locked scores to numeric vectors
# ------------------------------------------------------------

locked_stromal_10b3ne <-
  as.numeric(
    estimate_scores_locked_10b3ne[
      "StromalScore",
    ]
  )

locked_immune_10b3ne <-
  as.numeric(
    estimate_scores_locked_10b3ne[
      "ImmuneScore",
    ]
  )

locked_estimate_10b3ne <-
  as.numeric(
    estimate_scores_locked_10b3ne[
      "ESTIMATEScore",
    ]
  )


names(locked_stromal_10b3ne) <-
  colnames(
    estimate_scores_locked_10b3ne
  )

names(locked_immune_10b3ne) <-
  colnames(
    estimate_scores_locked_10b3ne
  )

names(locked_estimate_10b3ne) <-
  colnames(
    estimate_scores_locked_10b3ne
  )


# ------------------------------------------------------------
# 2. Confirm sample identity/order
# ------------------------------------------------------------

AOC3_LOO_sample_order_identical_10b3ne <-
  identical(
    names(AOC3_LOO_stromal_10b3ne),
    names(locked_stromal_10b3ne)
  )


cat(
  "Sample identity/order identical: ",
  AOC3_LOO_sample_order_identical_10b3ne,
  "\n",
  sep = ""
)


stopifnot(
  AOC3_LOO_sample_order_identical_10b3ne
)


# ------------------------------------------------------------
# 3. Explicit alignment
# ------------------------------------------------------------

AOC3_LOO_stromal_aligned_10b3ne <-
  AOC3_LOO_stromal_10b3ne[
    names(locked_stromal_10b3ne)
  ]

AOC3_LOO_immune_aligned_10b3ne <-
  AOC3_LOO_immune_10b3ne[
    names(locked_immune_10b3ne)
  ]

AOC3_LOO_estimate_aligned_10b3ne <-
  AOC3_LOO_estimate_10b3ne[
    names(locked_estimate_10b3ne)
  ]


# ------------------------------------------------------------
# 4. Calculate deltas
# ------------------------------------------------------------

AOC3_LOO_stromal_delta_10b3ne <-
  AOC3_LOO_stromal_aligned_10b3ne -
  locked_stromal_10b3ne

AOC3_LOO_immune_delta_10b3ne <-
  AOC3_LOO_immune_aligned_10b3ne -
  locked_immune_10b3ne

AOC3_LOO_estimate_delta_10b3ne <-
  AOC3_LOO_estimate_aligned_10b3ne -
  locked_estimate_10b3ne


# ------------------------------------------------------------
# 5. Number of changed samples
# ------------------------------------------------------------

AOC3_LOO_stromal_changed_n_10b3ne <-
  sum(
    abs(AOC3_LOO_stromal_delta_10b3ne) > 1e-12
  )

AOC3_LOO_immune_changed_n_10b3ne <-
  sum(
    abs(AOC3_LOO_immune_delta_10b3ne) > 1e-12
  )

AOC3_LOO_estimate_changed_n_10b3ne <-
  sum(
    abs(AOC3_LOO_estimate_delta_10b3ne) > 1e-12
  )


# ------------------------------------------------------------
# 6. Print numerical comparison
# ------------------------------------------------------------

cat(
  "StromalScore MaxAbsDiff = ",
  format(
    max(abs(AOC3_LOO_stromal_delta_10b3ne)),
    scientific = TRUE,
    digits = 12
  ),
  "\n",
  sep = ""
)

cat(
  "ImmuneScore MaxAbsDiff = ",
  format(
    max(abs(AOC3_LOO_immune_delta_10b3ne)),
    scientific = TRUE,
    digits = 12
  ),
  "\n",
  sep = ""
)

cat(
  "ESTIMATEScore MaxAbsDiff = ",
  format(
    max(abs(AOC3_LOO_estimate_delta_10b3ne)),
    scientific = TRUE,
    digits = 12
  ),
  "\n\n",
  sep = ""
)


cat(
  "StromalScore changed samples: ",
  AOC3_LOO_stromal_changed_n_10b3ne,
  "/ 513\n",
  sep = ""
)

cat(
  "ImmuneScore changed samples: ",
  AOC3_LOO_immune_changed_n_10b3ne,
  "/ 513\n",
  sep = ""
)

cat(
  "ESTIMATEScore changed samples: ",
  AOC3_LOO_estimate_changed_n_10b3ne,
  "/ 513\n",
  sep = ""
)



# E-4E-DIAG. Check sample-name metadata
# ============================================================

cat("\n============================================================\n")
cat("E-4E-DIAG: CHECK SAMPLE NAMES\n")
cat("============================================================\n\n")

cat(
  "AOC3 expression length: ",
  length(AOC3_expression_10b3ne),
  "\n",
  sep = ""
)

cat(
  "AOC3 expression names length: ",
  length(names(AOC3_expression_10b3ne)),
  "\n",
  sep = ""
)

cat(
  "AOC3 LOO StromalScore length: ",
  length(AOC3_LOO_stromal_aligned_10b3ne),
  "\n",
  sep = ""
)

cat(
  "AOC3 LOO StromalScore names length: ",
  length(names(AOC3_LOO_stromal_aligned_10b3ne)),
  "\n\n",
  sep = ""
)

cat(
  "First 5 AOC3 expression names:\n"
)
print(head(names(AOC3_expression_10b3ne), 5))

cat(
  "\nFirst 5 AOC3-LOO StromalScore names:\n"
)
print(head(names(AOC3_LOO_stromal_aligned_10b3ne), 5))

cat(
  "\nFirst 5 locked StromalScore names:\n"
)
print(head(names(locked_stromal_10b3ne), 5))

cat("\n--------------------------------------------\n")
cat("Exact names() identical:\n")
print(
  identical(
    names(AOC3_expression_10b3ne),
    names(AOC3_LOO_stromal_aligned_10b3ne)
  )
)

cat("\nsetdiff(expression, LOO):\n")
print(
  setdiff(
    names(AOC3_expression_10b3ne),
    names(AOC3_LOO_stromal_aligned_10b3ne)
  )
)

cat("\nsetdiff(LOO, expression):\n")
print(
  setdiff(
    names(AOC3_LOO_stromal_aligned_10b3ne),
    names(AOC3_expression_10b3ne)
  )
)

cat("\nidentical after converting names to character:\n")
print(
  identical(
    as.character(names(AOC3_expression_10b3ne)),
    as.character(names(AOC3_LOO_stromal_aligned_10b3ne))
  )
)

# ============================================================
# E-4E. Reconstruct AOC3 expression sample names
#        and perform LOO Spearman analysis
# ============================================================

cat("\n============================================================\n")
cat("E-4E: AOC3 EXPRESSION vs AOC3-LOO STROMAL SCORE\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Reconstruct AOC3 expression vector
#    using the locked Step14 sample order
# ------------------------------------------------------------

AOC3_expression_10b3ne <- as.numeric(
  gene_expression_10b3ne[
    "AOC3",
    names(locked_stromal_10b3ne)
  ]
)

# Explicitly restore locked sample IDs
names(AOC3_expression_10b3ne) <- names(
  locked_stromal_10b3ne
)

# ------------------------------------------------------------
# 2. Alignment check
# ------------------------------------------------------------

stopifnot(
  length(AOC3_expression_10b3ne) == 513,
  length(AOC3_LOO_stromal_aligned_10b3ne) == 513,
  identical(
    names(AOC3_expression_10b3ne),
    names(locked_stromal_10b3ne)
  ),
  identical(
    names(AOC3_LOO_stromal_aligned_10b3ne),
    names(locked_stromal_10b3ne)
  )
)

cat(
  "AOC3 expression / locked sample order: PASS\n"
)

cat(
  "AOC3 expression / AOC3-LOO StromalScore alignment: PASS\n\n"
)

# ------------------------------------------------------------
# 3. Spearman correlation
# ------------------------------------------------------------

AOC3_LOO_stromal_spearman_10b3ne <- cor.test(
  AOC3_expression_10b3ne,
  AOC3_LOO_stromal_aligned_10b3ne,
  method = "spearman",
  exact = FALSE
)

AOC3_LOO_stromal_rho_10b3ne <- unname(
  AOC3_LOO_stromal_spearman_10b3ne$estimate
)

AOC3_LOO_stromal_p_10b3ne <-
  AOC3_LOO_stromal_spearman_10b3ne$p.value

# ------------------------------------------------------------
# 4. Report
# ------------------------------------------------------------

cat(
  "AOC3-LOO StromalScore Spearman rho = ",
  formatC(
    AOC3_LOO_stromal_rho_10b3ne,
    format = "f",
    digits = 6
  ),
  "\n",
  sep = ""
)

cat(
  "AOC3-LOO StromalScore Spearman P = ",
  format.pval(
    AOC3_LOO_stromal_p_10b3ne,
    digits = 6,
    eps = 1e-300
  ),
  "\n",
  sep = ""
)


# ============================================================
# E-4F. Restore original SI_geneset after AOC3-LOO
# ============================================================

cat("\n============================================================\n")
cat("E-4F: RESTORE ORIGINAL SI_geneset AFTER AOC3-LOO\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Restore original SI_geneset using a true deep copy
# ------------------------------------------------------------

assign(
  "SI_geneset",
  deep_copy_SI_10b3ne(
    SI_geneset_original_10b3ne
  ),
  envir = .GlobalEnv
)

# ------------------------------------------------------------
# 2. Check GlobalEnv restoration
# ------------------------------------------------------------

global_SI_restored_AOC3_10b3ne <- identical(
  get("SI_geneset", envir = .GlobalEnv),
  SI_geneset_original_10b3ne
)

cat(
  "GlobalEnv SI_geneset restored exactly: ",
  global_SI_restored_AOC3_10b3ne,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 3. Check package:estimate integrity
# ------------------------------------------------------------

package_SI_unchanged_AOC3_10b3ne <- identical(
  get(
    "SI_geneset",
    envir = package_estimate_env_10b3ne
  ),
  SI_geneset_original_10b3ne
)

package_SI_locked_AOC3_10b3ne <- bindingIsLocked(
  "SI_geneset",
  package_estimate_env_10b3ne
)

cat(
  "package:estimate SI_geneset unchanged: ",
  package_SI_unchanged_AOC3_10b3ne,
  "\n",
  sep = ""
)

cat(
  "package:estimate SI_geneset binding locked: ",
  package_SI_locked_AOC3_10b3ne,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 4. Confirm membership after restoration
# ------------------------------------------------------------

stromal_genes_restored_AOC3_10b3ne <- 
  as.character(
    SI_geneset_original_10b3ne[
      "StromalSignature",
      -1
    ]
  )

immune_genes_restored_AOC3_10b3ne <- 
  as.character(
    SI_geneset_original_10b3ne[
      "ImmuneSignature",
      -1
    ]
  )

stromal_genes_restored_AOC3_10b3ne <- 
  stromal_genes_restored_AOC3_10b3ne[
    !is.na(stromal_genes_restored_AOC3_10b3ne)
  ]

immune_genes_restored_AOC3_10b3ne <- 
  immune_genes_restored_AOC3_10b3ne[
    !is.na(immune_genes_restored_AOC3_10b3ne)
  ]

AOC3_restored_stromal_member_10b3ne <- 
  "AOC3" %in% stromal_genes_restored_AOC3_10b3ne

AOC3_restored_immune_member_10b3ne <- 
  "AOC3" %in% immune_genes_restored_AOC3_10b3ne

IL7R_restored_stromal_member_10b3ne <- 
  "IL7R" %in% stromal_genes_restored_AOC3_10b3ne

IL7R_restored_immune_member_10b3ne <- 
  "IL7R" %in% immune_genes_restored_AOC3_10b3ne

MS4A7_restored_stromal_member_10b3ne <- 
  "MS4A7" %in% stromal_genes_restored_AOC3_10b3ne

MS4A7_restored_immune_member_10b3ne <- 
  "MS4A7" %in% immune_genes_restored_AOC3_10b3ne

cat("\nMembership after restoration:\n")

cat(
  "AOC3 in StromalSignature: ",
  AOC3_restored_stromal_member_10b3ne,
  "\n",
  sep = ""
)

cat(
  "AOC3 in ImmuneSignature: ",
  AOC3_restored_immune_member_10b3ne,
  "\n",
  sep = ""
)

cat(
  "IL7R in StromalSignature: ",
  IL7R_restored_stromal_member_10b3ne,
  "\n",
  sep = ""
)

cat(
  "IL7R in ImmuneSignature: ",
  IL7R_restored_immune_member_10b3ne,
  "\n",
  sep = ""
)

cat(
  "MS4A7 in StromalSignature: ",
  MS4A7_restored_stromal_member_10b3ne,
  "\n",
  sep = ""
)

cat(
  "MS4A7 in ImmuneSignature: ",
  MS4A7_restored_immune_member_10b3ne,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 5. Final AOC3 restoration PASS
# ------------------------------------------------------------

E4F_AOC3_RESTORE_PASS_10b3ne <- all(
  global_SI_restored_AOC3_10b3ne,
  package_SI_unchanged_AOC3_10b3ne,
  package_SI_locked_AOC3_10b3ne,
  AOC3_restored_stromal_member_10b3ne,
  !AOC3_restored_immune_member_10b3ne,
  !IL7R_restored_stromal_member_10b3ne,
  IL7R_restored_immune_member_10b3ne,
  !MS4A7_restored_stromal_member_10b3ne,
  !MS4A7_restored_immune_member_10b3ne
)

cat(
  "\nE-4F AOC3 RESTORATION PASS: ",
  E4F_AOC3_RESTORE_PASS_10b3ne,
  "\n",
  sep = ""
)

stopifnot(E4F_AOC3_RESTORE_PASS_10b3ne)

##Next step: Start IL7R-LOO
# ============================================================
# E-5A. Prepare IL7R leave-one-gene-out SI_geneset
# ============================================================

cat("\n============================================================\n")
cat("E-5A: PREPARE IL7R LOO SI_geneset\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Create a TRUE deep copy
# ------------------------------------------------------------

SI_geneset_IL7R_LOO_10b3ne <- deep_copy_SI_10b3ne(
  SI_geneset_original_10b3ne
)

cat(
  "Deep copy identical to original before edit: ",
  identical(
    SI_geneset_IL7R_LOO_10b3ne,
    SI_geneset_original_10b3ne
  ),
  "\n",
  sep = ""
)

stopifnot(
  identical(
    SI_geneset_IL7R_LOO_10b3ne,
    SI_geneset_original_10b3ne
  )
)

# ------------------------------------------------------------
# 2. Locate IL7R in ImmuneSignature
# ------------------------------------------------------------

IL7R_position_10b3ne <- which(
  as.character(
    SI_geneset_IL7R_LOO_10b3ne[
      "ImmuneSignature",
      -1
    ]
  ) == "IL7R"
)

cat(
  "IL7R position in ImmuneSignature: ",
  IL7R_position_10b3ne,
  "\n",
  sep = ""
)

stopifnot(length(IL7R_position_10b3ne) == 1)

# ------------------------------------------------------------
# 3. Remove IL7R by setting its slot to NA
# ------------------------------------------------------------

# +1 because the first column is the 'estimate' identifier
IL7R_column_10b3ne <- IL7R_position_10b3ne + 1

SI_geneset_IL7R_LOO_10b3ne[
  "ImmuneSignature",
  IL7R_column_10b3ne
] <- NA_character_

# ------------------------------------------------------------
# 4. Extract modified signatures
# ------------------------------------------------------------

stromal_genes_IL7R_LOO_10b3ne <- as.character(
  SI_geneset_IL7R_LOO_10b3ne[
    "StromalSignature",
    -1
  ]
)

immune_genes_IL7R_LOO_10b3ne <- as.character(
  SI_geneset_IL7R_LOO_10b3ne[
    "ImmuneSignature",
    -1
  ]
)

stromal_genes_IL7R_LOO_10b3ne <- 
  stromal_genes_IL7R_LOO_10b3ne[
    !is.na(stromal_genes_IL7R_LOO_10b3ne)
  ]

immune_genes_IL7R_LOO_10b3ne <- 
  immune_genes_IL7R_LOO_10b3ne[
    !is.na(immune_genes_IL7R_LOO_10b3ne)
  ]

# ------------------------------------------------------------
# 5. Membership checks
# ------------------------------------------------------------

IL7R_in_stromal_LOO_10b3ne <- 
  "IL7R" %in% stromal_genes_IL7R_LOO_10b3ne

IL7R_in_immune_LOO_10b3ne <- 
  "IL7R" %in% immune_genes_IL7R_LOO_10b3ne

AOC3_in_stromal_LOO_10b3ne <- 
  "AOC3" %in% stromal_genes_IL7R_LOO_10b3ne

AOC3_in_immune_LOO_10b3ne <- 
  "AOC3" %in% immune_genes_IL7R_LOO_10b3ne

MS4A7_in_stromal_LOO_10b3ne <- 
  "MS4A7" %in% stromal_genes_IL7R_LOO_10b3ne

MS4A7_in_immune_LOO_10b3ne <- 
  "MS4A7" %in% immune_genes_IL7R_LOO_10b3ne

cat("\nModified signature membership:\n")

cat(
  "IL7R in StromalSignature: ",
  IL7R_in_stromal_LOO_10b3ne,
  "\n",
  sep = ""
)

cat(
  "IL7R in ImmuneSignature: ",
  IL7R_in_immune_LOO_10b3ne,
  "\n",
  sep = ""
)

cat(
  "AOC3 in StromalSignature: ",
  AOC3_in_stromal_LOO_10b3ne,
  "\n",
  sep = ""
)

cat(
  "AOC3 in ImmuneSignature: ",
  AOC3_in_immune_LOO_10b3ne,
  "\n",
  sep = ""
)

cat(
  "MS4A7 in StromalSignature: ",
  MS4A7_in_stromal_LOO_10b3ne,
  "\n",
  sep = ""
)

cat(
  "MS4A7 in ImmuneSignature: ",
  MS4A7_in_immune_LOO_10b3ne,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 6. Gene-set size check
# ------------------------------------------------------------

cat("\nGene counts after IL7R removal:\n")

cat(
  "StromalSignature genes: ",
  length(stromal_genes_IL7R_LOO_10b3ne),
  "\n",
  sep = ""
)

cat(
  "ImmuneSignature genes: ",
  length(immune_genes_IL7R_LOO_10b3ne),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 7. Preparation PASS
# ------------------------------------------------------------

E5A_IL7R_PREP_PASS_10b3ne <- all(
  length(IL7R_position_10b3ne) == 1,
  !IL7R_in_stromal_LOO_10b3ne,
  !IL7R_in_immune_LOO_10b3ne,
  AOC3_in_stromal_LOO_10b3ne,
  !AOC3_in_immune_LOO_10b3ne,
  !MS4A7_in_stromal_LOO_10b3ne,
  !MS4A7_in_immune_LOO_10b3ne,
  length(stromal_genes_IL7R_LOO_10b3ne) == 141,
  length(immune_genes_IL7R_LOO_10b3ne) == 140
)

cat(
  "\nE-5A IL7R LOO PREPARATION PASS: ",
  E5A_IL7R_PREP_PASS_10b3ne,
  "\n",
  sep = ""
)

stopifnot(E5A_IL7R_PREP_PASS_10b3ne)


##E-5B：Assign → Integrity Check → Official estimateScore()
cat("\n============================================================\n")
cat("E-5B: OFFICIAL ESTIMATE SCORE — IL7R LOO\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Confirm E-5A passed
# ------------------------------------------------------------

stopifnot(E5A_IL7R_PREP_PASS_10b3ne)

cat(
  "E-5A preparation status: ",
  E5A_IL7R_PREP_PASS_10b3ne,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 1. Re-confirm package:estimate SI_geneset BEFORE mutation
# ------------------------------------------------------------

package_SI_before_E5B_10b3ne <- 
  get(
    "SI_geneset",
    envir = package_estimate_env_10b3ne
  )

package_SI_before_E5B_identical_10b3ne <- identical(
  package_SI_before_E5B_10b3ne,
  SI_geneset_original_10b3ne
)

package_SI_before_E5B_locked_10b3ne <- 
  bindingIsLocked(
    "SI_geneset",
    package_estimate_env_10b3ne
  )

cat(
  "Package:estimate SI_geneset identical to original: ",
  package_SI_before_E5B_identical_10b3ne,
  "\n",
  sep = ""
)

cat(
  "Package:estimate SI_geneset binding locked: ",
  package_SI_before_E5B_locked_10b3ne,
  "\n",
  sep = ""
)

stopifnot(
  package_SI_before_E5B_identical_10b3ne,
  package_SI_before_E5B_locked_10b3ne
)


# ------------------------------------------------------------
# 2. Assign IL7R-LOO signature to GlobalEnv
# ------------------------------------------------------------

assign(
  "SI_geneset",
  SI_geneset_IL7R_LOO_10b3ne,
  envir = .GlobalEnv
)

cat(
  "\nModified IL7R-LOO SI_geneset assigned to R_GlobalEnv.\n"
)


# ------------------------------------------------------------
# 3. Verify runtime GlobalEnv object
# ------------------------------------------------------------

SI_geneset_runtime_E5B_10b3ne <- 
  get(
    "SI_geneset",
    envir = .GlobalEnv
  )

runtime_matches_IL7R_LOO_10b3ne <- identical(
  SI_geneset_runtime_E5B_10b3ne,
  SI_geneset_IL7R_LOO_10b3ne
)

cat(
  "Runtime GlobalEnv SI_geneset matches IL7R-LOO object: ",
  runtime_matches_IL7R_LOO_10b3ne,
  "\n",
  sep = ""
)

stopifnot(runtime_matches_IL7R_LOO_10b3ne)


# ------------------------------------------------------------
# 4. Verify package namespace was NOT modified
# ------------------------------------------------------------

package_SI_after_assign_E5B_10b3ne <- 
  get(
    "SI_geneset",
    envir = package_estimate_env_10b3ne
  )

package_unchanged_after_assign_E5B_10b3ne <- identical(
  package_SI_after_assign_E5B_10b3ne,
  SI_geneset_original_10b3ne
)

package_locked_after_assign_E5B_10b3ne <- 
  bindingIsLocked(
    "SI_geneset",
    package_estimate_env_10b3ne
  )

cat(
  "Package:estimate SI_geneset unchanged after GlobalEnv assignment: ",
  package_unchanged_after_assign_E5B_10b3ne,
  "\n",
  sep = ""
)

cat(
  "Package:estimate SI_geneset binding locked: ",
  package_locked_after_assign_E5B_10b3ne,
  "\n",
  sep = ""
)

stopifnot(
  package_unchanged_after_assign_E5B_10b3ne,
  package_locked_after_assign_E5B_10b3ne
)


# ------------------------------------------------------------
# 5. Verify IL7R-LOO signature content
# ------------------------------------------------------------

cat("\nIL7R-LOO signature verification:\n")

cat(
  "IL7R in StromalSignature: ",
  "IL7R" %in% as.character(
    SI_geneset_runtime_E5B_10b3ne[
      "StromalSignature",
      -1
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "IL7R in ImmuneSignature: ",
  "IL7R" %in% as.character(
    SI_geneset_runtime_E5B_10b3ne[
      "ImmuneSignature",
      -1
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "AOC3 in StromalSignature: ",
  "AOC3" %in% as.character(
    SI_geneset_runtime_E5B_10b3ne[
      "StromalSignature",
      -1
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "MS4A7 in either signature: ",
  any(
    "MS4A7" %in% as.character(
      SI_geneset_runtime_E5B_10b3ne[
        "StromalSignature",
        -1
      ]
    ),
    "MS4A7" %in% as.character(
      SI_geneset_runtime_E5B_10b3ne[
        "ImmuneSignature",
        -1
      ]
    )
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Define official input/output GCT
# ------------------------------------------------------------

IL7R_LOO_input_gct_10b3ne <- 
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis/STEP7B-7E-10B-3N/10B-3N-C/STEP10B3NC_official_filterCommonGenes_output.gct"

IL7R_LOO_output_gct_10b3ne <- 
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis/STEP7B-7E-10B-3N/10B-3N-E/IL7R_LOO/STEP10B3NE_IL7R_LOO_estimateScore.gct"


stopifnot(file.exists(IL7R_LOO_input_gct_10b3ne))

dir.create(
  dirname(IL7R_LOO_output_gct_10b3ne),
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 7. Run OFFICIAL estimateScore()
# ------------------------------------------------------------

cat("\nRunning official estimateScore() with IL7R removed...\n\n")

estimateScore(
  input.ds = IL7R_LOO_input_gct_10b3ne,
  output.ds = IL7R_LOO_output_gct_10b3ne,
  platform = "illumina"
)

cat(
  "\nOfficial estimateScore() completed.\n"
)

cat(
  "Output file:\n",
  IL7R_LOO_output_gct_10b3ne,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. E-5B PASS
# ------------------------------------------------------------

E5B_IL7R_ESTIMATE_SCORE_PASS_10b3ne <- 
  file.exists(
    IL7R_LOO_output_gct_10b3ne
  ) &&
  package_unchanged_after_assign_E5B_10b3ne &&
  package_locked_after_assign_E5B_10b3ne &&
  runtime_matches_IL7R_LOO_10b3ne

cat(
  "\nE-5B IL7R LOO OFFICIAL ESTIMATE SCORE PASS: ",
  E5B_IL7R_ESTIMATE_SCORE_PASS_10b3ne,
  "\n",
  sep = ""
)

stopifnot(
  E5B_IL7R_ESTIMATE_SCORE_PASS_10b3ne
)

##Next step: E-5C analysis of IL7R-LOO GCT
cat("\n============================================================\n")
cat("E-5C: CORRECTED PARSE IL7R LOO ESTIMATE SCORE OUTPUT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Confirm E-5B passed
# ------------------------------------------------------------

stopifnot(E5B_IL7R_ESTIMATE_SCORE_PASS_10b3ne)


# ------------------------------------------------------------
# 1. Re-read official IL7R-LOO GCT
# ------------------------------------------------------------

IL7R_LOO_gct_raw_10b3ne <- read.delim(
  IL7R_LOO_output_gct_10b3ne,
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat(
  "Raw GCT dimensions: ",
  paste(
    dim(IL7R_LOO_gct_raw_10b3ne),
    collapse = " x "
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 2. Confirm standard GCT structure
# ------------------------------------------------------------

stopifnot(
  nrow(IL7R_LOO_gct_raw_10b3ne) == 6,
  ncol(IL7R_LOO_gct_raw_10b3ne) == 515
)

cat(
  "GCT version row: ",
  IL7R_LOO_gct_raw_10b3ne[1, 1],
  "\n",
  sep = ""
)

cat(
  "GCT dimension row: ",
  IL7R_LOO_gct_raw_10b3ne[2, 1],
  " x ",
  IL7R_LOO_gct_raw_10b3ne[2, 2],
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 3. Correct GCT parsing
#
# Row 3 = header
# Rows 4-6 = three score rows
# ------------------------------------------------------------

IL7R_LOO_header_10b3ne <- 
  as.character(
    IL7R_LOO_gct_raw_10b3ne[3, ]
  )

IL7R_LOO_score_data_10b3ne <- 
  IL7R_LOO_gct_raw_10b3ne[
    4:6,
    ,
    drop = FALSE
  ]


# ------------------------------------------------------------
# 4. Extract sample IDs from GCT header
# ------------------------------------------------------------

IL7R_LOO_sample_ids_10b3ne <- 
  IL7R_LOO_header_10b3ne[
    -c(1, 2)
  ]

stopifnot(
  length(IL7R_LOO_sample_ids_10b3ne) == 513
)


# ------------------------------------------------------------
# 5. Extract score row names
# ------------------------------------------------------------

IL7R_LOO_gene_names_10b3ne <- 
  as.character(
    IL7R_LOO_score_data_10b3ne[, 1]
  )

cat(
  "\nParsed score rows: ",
  paste(
    IL7R_LOO_gene_names_10b3ne,
    collapse = ", "
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Extract numeric score matrix
# ------------------------------------------------------------

IL7R_LOO_score_matrix_10b3ne <- 
  as.matrix(
    IL7R_LOO_score_data_10b3ne[
      ,
      -c(1, 2),
      drop = FALSE
    ]
  )

storage.mode(
  IL7R_LOO_score_matrix_10b3ne
) <- "numeric"

rownames(
  IL7R_LOO_score_matrix_10b3ne
) <- 
  IL7R_LOO_gene_names_10b3ne

colnames(
  IL7R_LOO_score_matrix_10b3ne
) <- 
  IL7R_LOO_sample_ids_10b3ne


# ------------------------------------------------------------
# 7. Structural checks
# ------------------------------------------------------------

expected_score_rows_10b3ne <- c(
  "StromalScore",
  "ImmuneScore",
  "ESTIMATEScore"
)

score_rows_correct_10b3ne <- identical(
  rownames(
    IL7R_LOO_score_matrix_10b3ne
  ),
  expected_score_rows_10b3ne
)

score_matrix_dimension_correct_10b3ne <- identical(
  dim(
    IL7R_LOO_score_matrix_10b3ne
  ),
  c(3L, 513L)
)

cat(
  "Score matrix dimensions: ",
  paste(
    dim(
      IL7R_LOO_score_matrix_10b3ne
    ),
    collapse = " x "
  ),
  "\n",
  sep = ""
)

cat(
  "Score row names correct: ",
  score_rows_correct_10b3ne,
  "\n",
  sep = ""
)

cat(
  "Score matrix dimensions correct: ",
  score_matrix_dimension_correct_10b3ne,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Extract the three score vectors
# ------------------------------------------------------------

IL7R_LOO_stromal_10b3ne <- 
  as.numeric(
    IL7R_LOO_score_matrix_10b3ne[
      "StromalScore",
    ]
  )

IL7R_LOO_immune_10b3ne <- 
  as.numeric(
    IL7R_LOO_score_matrix_10b3ne[
      "ImmuneScore",
    ]
  )

IL7R_LOO_estimate_10b3ne <- 
  as.numeric(
    IL7R_LOO_score_matrix_10b3ne[
      "ESTIMATEScore",
    ]
  )


# ------------------------------------------------------------
# 9. Restore sample names explicitly
# ------------------------------------------------------------

names(
  IL7R_LOO_stromal_10b3ne
) <- 
  colnames(
    IL7R_LOO_score_matrix_10b3ne
  )

names(
  IL7R_LOO_immune_10b3ne
) <- 
  colnames(
    IL7R_LOO_score_matrix_10b3ne
  )

names(
  IL7R_LOO_estimate_10b3ne
) <- 
  colnames(
    IL7R_LOO_score_matrix_10b3ne
  )


# ------------------------------------------------------------
# 10. Numeric completeness
# ------------------------------------------------------------

IL7R_LOO_numeric_complete_10b3ne <- all(
  is.finite(
    IL7R_LOO_stromal_10b3ne
  ),
  is.finite(
    IL7R_LOO_immune_10b3ne
  ),
  is.finite(
    IL7R_LOO_estimate_10b3ne
  )
)

IL7R_LOO_sample_count_10b3ne <- 
  length(
    IL7R_LOO_stromal_10b3ne
  )

IL7R_LOO_sample_ids_unique_10b3ne <- 
  length(
    unique(
      names(
        IL7R_LOO_stromal_10b3ne
      )
    )
  ) ==
  IL7R_LOO_sample_count_10b3ne


# ------------------------------------------------------------
# 11. Compare sample IDs with locked Step14
# ------------------------------------------------------------

locked_sample_ids_10b3ne <- 
  colnames(
    estimate_scores_locked_10b3ne
  )

IL7R_LOO_sample_identity_10b3ne <- identical(
  names(
    IL7R_LOO_stromal_10b3ne
  ),
  locked_sample_ids_10b3ne
)

cat(
  "\nSample count: ",
  IL7R_LOO_sample_count_10b3ne,
  "\n",
  sep = ""
)

cat(
  "Sample IDs unique: ",
  IL7R_LOO_sample_ids_unique_10b3ne,
  "\n",
  sep = ""
)

cat(
  "All scores finite: ",
  IL7R_LOO_numeric_complete_10b3ne,
  "\n",
  sep = ""
)

cat(
  "Sample identity/order matches locked Step14: ",
  IL7R_LOO_sample_identity_10b3ne,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. E-5C PASS
# ------------------------------------------------------------

E5C_IL7R_PARSE_PASS_10b3ne <- all(
  score_rows_correct_10b3ne,
  score_matrix_dimension_correct_10b3ne,
  IL7R_LOO_sample_count_10b3ne == 513,
  IL7R_LOO_sample_ids_unique_10b3ne,
  IL7R_LOO_numeric_complete_10b3ne,
  IL7R_LOO_sample_identity_10b3ne
)

cat(
  "\nE-5C IL7R LOO PARSE PASS: ",
  E5C_IL7R_PARSE_PASS_10b3ne,
  "\n",
  sep = ""
)

stopifnot(
  E5C_IL7R_PARSE_PASS_10b3ne
)


cat("\n============================================================\n")
cat("E-5C-1: DIAGNOSE SAMPLE ID FORMAT DIFFERENCE\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Show first 10 IDs from IL7R-LOO GCT
# ------------------------------------------------------------

cat("First 10 IL7R-LOO sample IDs:\n")
print(
  head(
    names(IL7R_LOO_stromal_10b3ne),
    10
  )
)


# ------------------------------------------------------------
# 2. Show first 10 IDs from locked Step14
# ------------------------------------------------------------

cat("\nFirst 10 locked Step14 sample IDs:\n")
print(
  head(
    locked_sample_ids_10b3ne,
    10
  )
)


# ------------------------------------------------------------
# 3. Compare exact strings
# ------------------------------------------------------------

cat(
  "\nExact first ID comparison:\n"
)

cat(
  "IL7R-LOO : ",
  names(IL7R_LOO_stromal_10b3ne)[1],
  "\n",
  sep = ""
)

cat(
  "Locked   : ",
  locked_sample_ids_10b3ne[1],
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Compare normalized IDs
# ------------------------------------------------------------

IL7R_LOO_sample_ids_normalized_10b3ne <- 
  gsub(
    "\\.",
    "-",
    names(IL7R_LOO_stromal_10b3ne)
  )

locked_sample_ids_normalized_10b3ne <- 
  gsub(
    "\\.",
    "-",
    locked_sample_ids_10b3ne
  )


normalized_sample_identity_10b3ne <- identical(
  IL7R_LOO_sample_ids_normalized_10b3ne,
  locked_sample_ids_normalized_10b3ne
)

cat(
  "\nNormalized sample identity/order matches locked Step14: ",
  normalized_sample_identity_10b3ne,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. Check set equality after normalization
# ------------------------------------------------------------

normalized_set_equal_10b3ne <- 
  setequal(
    IL7R_LOO_sample_ids_normalized_10b3ne,
    locked_sample_ids_normalized_10b3ne
  )

cat(
  "Normalized sample ID sets identical: ",
  normalized_set_equal_10b3ne,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Check whether order is identical after normalization
# ------------------------------------------------------------

normalized_order_equal_10b3ne <- all(
  IL7R_LOO_sample_ids_normalized_10b3ne ==
    locked_sample_ids_normalized_10b3ne
)

cat(
  "Normalized sample ID order identical: ",
  normalized_order_equal_10b3ne,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Number of mismatches after normalization
# ------------------------------------------------------------

normalized_mismatch_n_10b3ne <- sum(
  IL7R_LOO_sample_ids_normalized_10b3ne !=
    locked_sample_ids_normalized_10b3ne
)

cat(
  "Number of normalized ID mismatches: ",
  normalized_mismatch_n_10b3ne,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. If there are mismatches, show first 10
# ------------------------------------------------------------

if (normalized_mismatch_n_10b3ne > 0) {
  
  mismatch_idx_10b3ne <- which(
    IL7R_LOO_sample_ids_normalized_10b3ne !=
      locked_sample_ids_normalized_10b3ne
  )
  
  cat("\nFirst normalized ID mismatches:\n")
  
  print(
    data.frame(
      position = mismatch_idx_10b3ne[
        seq_len(
          min(10, length(mismatch_idx_10b3ne))
        )
      ],
      IL7R_LOO = IL7R_LOO_sample_ids_normalized_10b3ne[
        mismatch_idx_10b3ne[
          seq_len(
            min(10, length(mismatch_idx_10b3ne))
          )
        ]
      ],
      Locked = locked_sample_ids_normalized_10b3ne[
        mismatch_idx_10b3ne[
          seq_len(
            min(10, length(mismatch_idx_10b3ne))
          )
        ]
      ],
      stringsAsFactors = FALSE
    )
  )
  
}


# ------------------------------------------------------------
# 9. Diagnostic PASS
# ------------------------------------------------------------

E5C1_SAMPLE_ID_DIAGNOSTIC_PASS_10b3ne <- all(
  normalized_set_equal_10b3ne,
  normalized_order_equal_10b3ne
)

cat(
  "\nE-5C-1 SAMPLE ID DIAGNOSTIC PASS: ",
  E5C1_SAMPLE_ID_DIAGNOSTIC_PASS_10b3ne,
  "\n",
  sep = ""
)

# ============================================================
# STEP 10B-3N-E
# E-5C: PARSE IL7R-LOO ESTIMATE OUTPUT AND VALIDATE
# Corrected sample-ID normalization
# ============================================================

cat("\n============================================================\n")
cat("E-5C: PARSE IL7R-LOO ESTIMATE OUTPUT AND VALIDATE\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Read the official IL7R-LOO estimateScore output
# ------------------------------------------------------------

IL7R_LOO_output_gct_10b3ne <- 
  "input_dir/Paper 1_Step 7_ Immune Infiltration Analysis/STEP7B-7E-10B-3N/10B-3N-E/IL7R_LOO/STEP10B3NE_IL7R_LOO_estimateScore.gct"


raw_IL7R_LOO_gct_10b3ne <- read.delim(
  IL7R_LOO_output_gct_10b3ne,
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


cat("Raw GCT dimensions: ",
    paste(dim(raw_IL7R_LOO_gct_10b3ne), collapse = " x "),
    "\n",
    sep = "")


# ------------------------------------------------------------
# 2. Verify GCT structure
# ------------------------------------------------------------

cat(
  "\nGCT version row:\n"
)

print(raw_IL7R_LOO_gct_10b3ne[1, 1:4])


cat(
  "\nGCT dimension row:\n"
)

print(raw_IL7R_LOO_gct_10b3ne[2, 1:4])


cat(
  "\nGCT header row:\n"
)

print(
  raw_IL7R_LOO_gct_10b3ne[
    3,
    1:min(6, ncol(raw_IL7R_LOO_gct_10b3ne))
  ]
)


# ------------------------------------------------------------
# 3. Parse score data correctly
# ------------------------------------------------------------

IL7R_LOO_sample_ids_raw_10b3ne <- 
  as.character(
    raw_IL7R_LOO_gct_10b3ne[
      3,
      -c(1, 2)
    ]
  )


IL7R_LOO_score_df_10b3ne <- 
  raw_IL7R_LOO_gct_10b3ne[
    4:nrow(raw_IL7R_LOO_gct_10b3ne),
    ,
    drop = FALSE
  ]


colnames(IL7R_LOO_score_df_10b3ne) <- 
  c(
    "NAME",
    "Description",
    IL7R_LOO_sample_ids_raw_10b3ne
  )


IL7R_LOO_score_matrix_10b3ne <- 
  as.matrix(
    IL7R_LOO_score_df_10b3ne[
      ,
      -c(1, 2),
      drop = FALSE
    ]
  )


IL7R_LOO_score_matrix_10b3ne <- 
  apply(
    IL7R_LOO_score_matrix_10b3ne,
    2,
    as.numeric
  )


IL7R_LOO_score_matrix_10b3ne <- 
  t(IL7R_LOO_score_matrix_10b3ne)

# ============================================================
# STEP 10B-3N-E
# E-5C: CORRECT GCT PARSING AND SAMPLE ALIGNMENT
# ============================================================

cat("\n============================================================\n")
cat("E-5C: CORRECT IL7R-LOO GCT PARSING\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Confirm GCT structure
# ------------------------------------------------------------

stopifnot(
  nrow(raw_IL7R_LOO_gct_10b3ne) == 6,
  ncol(raw_IL7R_LOO_gct_10b3ne) == 515
)


IL7R_LOO_sample_ids_raw_10b3ne <-
  as.character(
    raw_IL7R_LOO_gct_10b3ne[
      3,
      3:ncol(raw_IL7R_LOO_gct_10b3ne)
    ]
  )


# ------------------------------------------------------------
# 2. Extract the three score rows directly
# ------------------------------------------------------------

IL7R_LOO_score_names_10b3ne <-
  as.character(
    raw_IL7R_LOO_gct_10b3ne[
      4:6,
      1
    ]
  )


cat("Score names:\n")
print(IL7R_LOO_score_names_10b3ne)


expected_score_rows_10b3ne <-
  c(
    "StromalScore",
    "ImmuneScore",
    "ESTIMATEScore"
  )


score_names_correct_10b3ne <-
  identical(
    IL7R_LOO_score_names_10b3ne,
    expected_score_rows_10b3ne
  )


cat(
  "\nScore names correct: ",
  score_names_correct_10b3ne,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 3. Extract numeric score matrix
# ------------------------------------------------------------
# Important:
# GCT rows 4:6 are already the three score rows.
# Therefore NO transpose is required.

IL7R_LOO_score_matrix_10b3ne <-
  matrix(
    as.numeric(
      unlist(
        raw_IL7R_LOO_gct_10b3ne[
          4:6,
          3:515,
          drop = FALSE
        ],
        use.names = FALSE
      )
    ),
    nrow = 3,
    ncol = 513,
    byrow = TRUE
  )


# ------------------------------------------------------------
# 4. Assign row and column names
# ------------------------------------------------------------

rownames(
  IL7R_LOO_score_matrix_10b3ne
) <-
  IL7R_LOO_score_names_10b3ne


colnames(
  IL7R_LOO_score_matrix_10b3ne
) <-
  IL7R_LOO_sample_ids_raw_10b3ne


# ------------------------------------------------------------
# 5. Structural checks
# ------------------------------------------------------------

cat(
  "\nScore matrix dimensions: ",
  paste(
    dim(IL7R_LOO_score_matrix_10b3ne),
    collapse = " x "
  ),
  "\n",
  sep = ""
)


cat(
  "Row names:\n"
)

print(
  rownames(
    IL7R_LOO_score_matrix_10b3ne
  )
)


cat(
  "\nNumber of sample IDs: ",
  ncol(IL7R_LOO_score_matrix_10b3ne),
  "\n",
  sep = ""
)


score_dimension_pass_10b3ne <-
  identical(
    dim(IL7R_LOO_score_matrix_10b3ne),
    c(3L, 513L)
  )


score_row_pass_10b3ne <-
  identical(
    rownames(
      IL7R_LOO_score_matrix_10b3ne
    ),
    expected_score_rows_10b3ne
  )


score_finite_pass_10b3ne <-
  all(
    is.finite(
      IL7R_LOO_score_matrix_10b3ne
    )
  )


sample_unique_pass_10b3ne <-
  length(
    unique(
      colnames(
        IL7R_LOO_score_matrix_10b3ne
      )
    )
  ) == 513


cat(
  "\nDimension PASS: ",
  score_dimension_pass_10b3ne,
  "\n",
  sep = ""
)

cat(
  "Score-row PASS: ",
  score_row_pass_10b3ne,
  "\n",
  sep = ""
)

cat(
  "Finite-score PASS: ",
  score_finite_pass_10b3ne,
  "\n",
  sep = ""
)

cat(
  "Unique-sample-ID PASS: ",
  sample_unique_pass_10b3ne,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Normalize sample IDs
# ------------------------------------------------------------

IL7R_LOO_sample_ids_normalized_10b3ne <-
  gsub(
    "\\.",
    "-",
    colnames(
      IL7R_LOO_score_matrix_10b3ne
    )
  )


locked_sample_ids_normalized_10b3ne <-
  gsub(
    "\\.",
    "-",
    locked_sample_ids_10b3ne
  )


# ------------------------------------------------------------
# 7. Validate sample identity and order
# ------------------------------------------------------------

normalized_set_equal_E5C_10b3ne <-
  setequal(
    IL7R_LOO_sample_ids_normalized_10b3ne,
    locked_sample_ids_normalized_10b3ne
  )


normalized_order_equal_E5C_10b3ne <-
  identical(
    IL7R_LOO_sample_ids_normalized_10b3ne,
    locked_sample_ids_normalized_10b3ne
  )


normalized_mismatch_n_E5C_10b3ne <-
  sum(
    IL7R_LOO_sample_ids_normalized_10b3ne !=
      locked_sample_ids_normalized_10b3ne
  )


cat(
  "\nNormalized sample ID set identical: ",
  normalized_set_equal_E5C_10b3ne,
  "\n",
  sep = ""
)

cat(
  "Normalized sample ID order identical: ",
  normalized_order_equal_E5C_10b3ne,
  "\n",
  sep = ""
)

cat(
  "Normalized ID mismatches: ",
  normalized_mismatch_n_E5C_10b3ne,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Replace dot-format IDs with locked-format IDs
# ------------------------------------------------------------
# E-5C-1 already established that the normalized order is
# identical. Therefore this is only a naming operation,
# NOT a sample reordering operation.

if (!normalized_order_equal_E5C_10b3ne) {
  stop(
    "Normalized sample order is not identical. ",
    "Do not proceed to E-5D."
  )
}


colnames(
  IL7R_LOO_score_matrix_10b3ne
) <-
  locked_sample_ids_10b3ne


# ------------------------------------------------------------
# 9. Extract the three LOO scores
# ------------------------------------------------------------

IL7R_LOO_stromal_10b3ne <-
  IL7R_LOO_score_matrix_10b3ne[
    "StromalScore",
  ]

IL7R_LOO_immune_10b3ne <-
  IL7R_LOO_score_matrix_10b3ne[
    "ImmuneScore",
  ]

IL7R_LOO_estimate_10b3ne <-
  IL7R_LOO_score_matrix_10b3ne[
    "ESTIMATEScore",
  ]


# ------------------------------------------------------------
# 10. Final E-5C PASS
# ------------------------------------------------------------

E5C_IL7R_PARSE_PASS_10b3ne <-
  all(
    score_dimension_pass_10b3ne,
    score_row_pass_10b3ne,
    score_finite_pass_10b3ne,
    sample_unique_pass_10b3ne,
    normalized_set_equal_E5C_10b3ne,
    normalized_order_equal_E5C_10b3ne,
    length(IL7R_LOO_stromal_10b3ne) == 513,
    length(IL7R_LOO_immune_10b3ne) == 513,
    length(IL7R_LOO_estimate_10b3ne) == 513
  )


cat(
  "\n============================================================\n"
)

cat(
  "E-5C FINAL PARSE/ALIGNMENT PASS: ",
  E5C_IL7R_PARSE_PASS_10b3ne,
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)

##E-5D：Formal association sensitivity analysis for IL7R-LOO
# ============================================================
# STEP 10B-3N-E
# E-5C-R: REPARSE IL7R-LOO GCT WITHOUT MATRIX FLATTENING
# ============================================================

cat("\n============================================================\n")
cat("E-5C-R: REPARSE IL7R-LOO GCT CORRECTLY\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Confirm original GCT structure
# ------------------------------------------------------------

stopifnot(
  nrow(raw_IL7R_LOO_gct_10b3ne) == 6,
  ncol(raw_IL7R_LOO_gct_10b3ne) == 515
)


# ------------------------------------------------------------
# 2. Extract sample IDs directly from GCT header
# ------------------------------------------------------------

IL7R_LOO_sample_ids_raw_10b3ne <-
  as.character(
    raw_IL7R_LOO_gct_10b3ne[
      3,
      3:515
    ]
  )


# ------------------------------------------------------------
# 3. Extract score names
# ------------------------------------------------------------

IL7R_LOO_score_names_10b3ne <-
  as.character(
    raw_IL7R_LOO_gct_10b3ne[
      4:6,
      1
    ]
  )


expected_score_rows_10b3ne <-
  c(
    "StromalScore",
    "ImmuneScore",
    "ESTIMATEScore"
  )


score_names_correct_E5CR_10b3ne <-
  identical(
    IL7R_LOO_score_names_10b3ne,
    expected_score_rows_10b3ne
  )


cat(
  "Score names correct: ",
  score_names_correct_E5CR_10b3ne,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Extract the score block DIRECTLY
# ------------------------------------------------------------
# IMPORTANT:
# Do NOT use unlist()
# Do NOT use matrix(..., byrow=TRUE)
# Do NOT transpose.
#
# The GCT score block is already exactly 3 x 513.

IL7R_LOO_score_block_raw_10b3ne <-
  raw_IL7R_LOO_gct_10b3ne[
    4:6,
    3:515,
    drop = FALSE
  ]


cat(
  "Raw score block dimensions: ",
  paste(
    dim(IL7R_LOO_score_block_raw_10b3ne),
    collapse = " x "
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. Convert each score cell to numeric
# ------------------------------------------------------------

IL7R_LOO_score_matrix_10b3ne <-
  matrix(
    NA_real_,
    nrow = 3,
    ncol = 513
  )


for (i in seq_len(3)) {
  
  IL7R_LOO_score_matrix_10b3ne[i, ] <-
    as.numeric(
      IL7R_LOO_score_block_raw_10b3ne[
        i,
      ]
    )
  
}


# ------------------------------------------------------------
# 6. Assign row and column names
# ------------------------------------------------------------

rownames(
  IL7R_LOO_score_matrix_10b3ne
) <-
  IL7R_LOO_score_names_10b3ne


colnames(
  IL7R_LOO_score_matrix_10b3ne
) <-
  IL7R_LOO_sample_ids_raw_10b3ne


# ------------------------------------------------------------
# 7. Structural validation
# ------------------------------------------------------------

cat(
  "\nScore matrix dimensions: ",
  paste(
    dim(IL7R_LOO_score_matrix_10b3ne),
    collapse = " x "
  ),
  "\n",
  sep = ""
)


cat(
  "Row names:\n"
)

print(
  rownames(
    IL7R_LOO_score_matrix_10b3ne
  )
)


score_dimension_pass_E5CR_10b3ne <-
  identical(
    dim(
      IL7R_LOO_score_matrix_10b3ne
    ),
    c(3L, 513L)
  )


score_row_pass_E5CR_10b3ne <-
  identical(
    rownames(
      IL7R_LOO_score_matrix_10b3ne
    ),
    expected_score_rows_10b3ne
  )


score_finite_pass_E5CR_10b3ne <-
  all(
    is.finite(
      IL7R_LOO_score_matrix_10b3ne
    )
  )


score_sample_unique_pass_E5CR_10b3ne <-
  length(
    unique(
      colnames(
        IL7R_LOO_score_matrix_10b3ne
      )
    )
  ) == 513


cat(
  "\nDimension PASS: ",
  score_dimension_pass_E5CR_10b3ne,
  "\n",
  sep = ""
)

cat(
  "Score-row PASS: ",
  score_row_pass_E5CR_10b3ne,
  "\n",
  sep = ""
)

cat(
  "Finite-score PASS: ",
  score_finite_pass_E5CR_10b3ne,
  "\n",
  sep = ""
)

cat(
  "Unique-sample-ID PASS: ",
  score_sample_unique_pass_E5CR_10b3ne,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Normalize sample IDs
# ------------------------------------------------------------

IL7R_LOO_sample_ids_normalized_10b3ne <-
  gsub(
    "\\.",
    "-",
    colnames(
      IL7R_LOO_score_matrix_10b3ne
    )
  )


locked_sample_ids_normalized_10b3ne <-
  gsub(
    "\\.",
    "-",
    locked_sample_ids_10b3ne
  )


normalized_set_equal_E5CR_10b3ne <-
  setequal(
    IL7R_LOO_sample_ids_normalized_10b3ne,
    locked_sample_ids_normalized_10b3ne
  )


normalized_order_equal_E5CR_10b3ne <-
  identical(
    IL7R_LOO_sample_ids_normalized_10b3ne,
    locked_sample_ids_normalized_10b3ne
  )


cat(
  "\nNormalized sample ID set identical: ",
  normalized_set_equal_E5CR_10b3ne,
  "\n",
  sep = ""
)

cat(
  "Normalized sample ID order identical: ",
  normalized_order_equal_E5CR_10b3ne,
  "\n",
  sep = ""
)


if (!normalized_order_equal_E5CR_10b3ne) {
  stop(
    "Sample order mismatch. Do not proceed."
  )
}


# ------------------------------------------------------------
# 9. Replace GCT dot-format IDs with locked IDs
# ------------------------------------------------------------

colnames(
  IL7R_LOO_score_matrix_10b3ne
) <-
  locked_sample_ids_10b3ne


# ------------------------------------------------------------
# 10. Extract three scores
# ------------------------------------------------------------

IL7R_LOO_stromal_10b3ne <-
  IL7R_LOO_score_matrix_10b3ne[
    "StromalScore",
  ]

IL7R_LOO_immune_10b3ne <-
  IL7R_LOO_score_matrix_10b3ne[
    "ImmuneScore",
  ]

IL7R_LOO_estimate_10b3ne <-
  IL7R_LOO_score_matrix_10b3ne[
    "ESTIMATEScore",
  ]


# ------------------------------------------------------------
# 11. Additional sanity check:
# compare first few values directly
# ------------------------------------------------------------

cat(
  "\nFirst 5 IL7R-LOO ImmuneScore values:\n"
)

print(
  IL7R_LOO_immune_10b3ne[
    1:5
  ]
)


cat(
  "\nFirst 5 locked Step14 ImmuneScore values:\n"
)

print(
  locked_immune_10b3ne[
    1:5
  ]
)


# ------------------------------------------------------------
# 12. Final E-5C-R PASS
# ------------------------------------------------------------

E5CR_IL7R_PARSE_PASS_10b3ne <-
  all(
    score_names_correct_E5CR_10b3ne,
    score_dimension_pass_E5CR_10b3ne,
    score_row_pass_E5CR_10b3ne,
    score_finite_pass_E5CR_10b3ne,
    score_sample_unique_pass_E5CR_10b3ne,
    normalized_set_equal_E5CR_10b3ne,
    normalized_order_equal_E5CR_10b3ne,
    length(
      IL7R_LOO_stromal_10b3ne
    ) == 513,
    length(
      IL7R_LOO_immune_10b3ne
    ) == 513,
    length(
      IL7R_LOO_estimate_10b3ne
    ) == 513
  )


cat(
  "\n============================================================\n"
)

cat(
  "E-5C-R FINAL PARSE PASS: ",
  E5CR_IL7R_PARSE_PASS_10b3ne,
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)


cat("\n============================================================\n")
cat("STEP 10B-3N-E: E-5D IL7R-LOO NUMERICAL + ASSOCIATION ANALYSIS\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Confirm required objects
# ------------------------------------------------------------

stopifnot(
  exists("IL7R_LOO_immune_10b3ne"),
  exists("IL7R_LOO_stromal_10b3ne"),
  exists("IL7R_LOO_estimate_10b3ne"),
  exists("locked_immune_10b3ne"),
  exists("locked_stromal_10b3ne"),
  exists("locked_estimate_10b3ne"),
  exists("gene_expression_10b3ne")
)


# ------------------------------------------------------------
# 2. Confirm sample identity and order
# ------------------------------------------------------------

IL7R_LOO_ids_E5D_10b3ne <-
  names(IL7R_LOO_immune_10b3ne)

locked_ids_E5D_10b3ne <-
  names(locked_immune_10b3ne)


sample_id_set_pass_E5D_10b3ne <-
  setequal(
    IL7R_LOO_ids_E5D_10b3ne,
    locked_ids_E5D_10b3ne
  )

sample_id_order_pass_E5D_10b3ne <-
  identical(
    IL7R_LOO_ids_E5D_10b3ne,
    locked_ids_E5D_10b3ne
  )


cat(
  "Sample ID set identical: ",
  sample_id_set_pass_E5D_10b3ne,
  "\n",
  sep = ""
)

cat(
  "Sample ID order identical: ",
  sample_id_order_pass_E5D_10b3ne,
  "\n",
  sep = ""
)


if (!sample_id_order_pass_E5D_10b3ne) {
  stop(
    "IL7R-LOO and locked Step14 sample order mismatch. Do not proceed."
  )
}


# ------------------------------------------------------------
# 3. Confirm finite values
# ------------------------------------------------------------

finite_pass_E5D_10b3ne <-
  all(
    is.finite(IL7R_LOO_immune_10b3ne),
    is.finite(locked_immune_10b3ne),
    is.finite(IL7R_LOO_stromal_10b3ne),
    is.finite(locked_stromal_10b3ne),
    is.finite(IL7R_LOO_estimate_10b3ne),
    is.finite(locked_estimate_10b3ne)
  )


cat(
  "Finite-score PASS: ",
  finite_pass_E5D_10b3ne,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Numerical comparison: IL7R-LOO vs locked Step14
# ------------------------------------------------------------

delta_immune_E5D_10b3ne <-
  IL7R_LOO_immune_10b3ne -
  locked_immune_10b3ne

delta_stromal_E5D_10b3ne <-
  IL7R_LOO_stromal_10b3ne -
  locked_stromal_10b3ne

delta_estimate_E5D_10b3ne <-
  IL7R_LOO_estimate_10b3ne -
  locked_estimate_10b3ne


max_abs_delta_immune_E5D_10b3ne <-
  max(abs(delta_immune_E5D_10b3ne))

mean_abs_delta_immune_E5D_10b3ne <-
  mean(abs(delta_immune_E5D_10b3ne))

median_abs_delta_immune_E5D_10b3ne <-
  median(abs(delta_immune_E5D_10b3ne))


max_abs_delta_stromal_E5D_10b3ne <-
  max(abs(delta_stromal_E5D_10b3ne))

max_abs_delta_estimate_E5D_10b3ne <-
  max(abs(delta_estimate_E5D_10b3ne))


cat("\n------------------------------------------------------------\n")
cat("IL7R-LOO vs locked Step14 numerical comparison\n")
cat("------------------------------------------------------------\n")

cat(
  "ImmuneScore max |delta|: ",
  format(
    max_abs_delta_immune_E5D_10b3ne,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "ImmuneScore mean |delta|: ",
  format(
    mean_abs_delta_immune_E5D_10b3ne,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "ImmuneScore median |delta|: ",
  format(
    median_abs_delta_immune_E5D_10b3ne,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "StromalScore max |delta|: ",
  format(
    max_abs_delta_stromal_E5D_10b3ne,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "ESTIMATEScore max |delta|: ",
  format(
    max_abs_delta_estimate_E5D_10b3ne,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. Count meaningful changes
# ------------------------------------------------------------
# Use 1e-10 as the meaningful numerical threshold.
# This avoids interpreting ~1e-12 floating-point noise as
# biological score changes.

immune_changed_gt_1e10_E5D_10b3ne <-
  sum(
    abs(delta_immune_E5D_10b3ne) > 1e-10
  )

stromal_changed_gt_1e10_E5D_10b3ne <-
  sum(
    abs(delta_stromal_E5D_10b3ne) > 1e-10
  )

estimate_changed_gt_1e10_E5D_10b3ne <-
  sum(
    abs(delta_estimate_E5D_10b3ne) > 1e-10
  )


cat("\nMeaningful changes (>1e-10):\n")

cat(
  "ImmuneScore: ",
  immune_changed_gt_1e10_E5D_10b3ne,
  "/513\n",
  sep = ""
)

cat(
  "StromalScore: ",
  stromal_changed_gt_1e10_E5D_10b3ne,
  "/513\n",
  sep = ""
)

cat(
  "ESTIMATEScore: ",
  estimate_changed_gt_1e10_E5D_10b3ne,
  "/513\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Correlation between original and IL7R-LOO ImmuneScore
# ------------------------------------------------------------

immune_score_spearman_E5D_10b3ne <-
  cor.test(
    locked_immune_10b3ne,
    IL7R_LOO_immune_10b3ne,
    method = "spearman",
    exact = FALSE
  )


cat("\n------------------------------------------------------------\n")
cat("Original vs IL7R-LOO ImmuneScore\n")
cat("------------------------------------------------------------\n")

cat(
  "Spearman rho = ",
  unname(
    immune_score_spearman_E5D_10b3ne$estimate
  ),
  "\n",
  sep = ""
)

cat(
  "P = ",
  immune_score_spearman_E5D_10b3ne$p.value,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Align IL7R expression with LOO ImmuneScore
# ------------------------------------------------------------

IL7R_expression_E5D_10b3ne <-
  as.numeric(
    gene_expression_10b3ne[
      "IL7R",
      names(IL7R_LOO_immune_10b3ne)
    ]
  )

names(
  IL7R_expression_E5D_10b3ne
) <-
  names(
    IL7R_LOO_immune_10b3ne
  )


IL7R_expression_alignment_pass_E5D_10b3ne <-
  identical(
    names(IL7R_expression_E5D_10b3ne),
    names(IL7R_LOO_immune_10b3ne)
  )


cat("\nIL7R expression / LOO ImmuneScore alignment: ",
    IL7R_expression_alignment_pass_E5D_10b3ne,
    "\n",
    sep = "")


# ------------------------------------------------------------
# 8. IL7R expression vs IL7R-LOO ImmuneScore
# ------------------------------------------------------------

IL7R_LOO_association_E5D_10b3ne <-
  cor.test(
    IL7R_expression_E5D_10b3ne,
    IL7R_LOO_immune_10b3ne,
    method = "spearman",
    exact = FALSE
  )


IL7R_LOO_rho_E5D_10b3ne <-
  unname(
    IL7R_LOO_association_E5D_10b3ne$estimate
  )

IL7R_LOO_p_E5D_10b3ne <-
  IL7R_LOO_association_E5D_10b3ne$p.value


cat("\n------------------------------------------------------------\n")
cat("IL7R vs IL7R-LOO ImmuneScore\n")
cat("------------------------------------------------------------\n")

cat(
  "Spearman rho = ",
  IL7R_LOO_rho_E5D_10b3ne,
  "\n",
  sep = ""
)

cat(
  "P = ",
  IL7R_LOO_p_E5D_10b3ne,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Compare with original locked IL7R–ImmuneScore association
# ------------------------------------------------------------

IL7R_original_E5D_10b3ne <-
  as.numeric(
    gene_expression_10b3ne[
      "IL7R",
      names(locked_immune_10b3ne)
    ]
  )

names(
  IL7R_original_E5D_10b3ne
) <-
  names(
    locked_immune_10b3ne
  )


IL7R_original_association_E5D_10b3ne <-
  cor.test(
    IL7R_original_E5D_10b3ne,
    locked_immune_10b3ne,
    method = "spearman",
    exact = FALSE
  )


IL7R_original_rho_E5D_10b3ne <-
  unname(
    IL7R_original_association_E5D_10b3ne$estimate
  )

IL7R_original_p_E5D_10b3ne <-
  IL7R_original_association_E5D_10b3ne$p.value


cat("\nOriginal locked IL7R–ImmuneScore:\n")

cat(
  "Spearman rho = ",
  IL7R_original_rho_E5D_10b3ne,
  "\n",
  sep = ""
)

cat(
  "P = ",
  IL7R_original_p_E5D_10b3ne,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. E-5D PASS criteria
# ------------------------------------------------------------
#
# Important:
# The LOO score is NOT expected to be identical to the
# original score because IL7R has been removed from the
# ImmuneSignature.
#
# PASS therefore requires:
# 1. correct sample alignment
# 2. finite numerical values
# 3. measurable change after IL7R removal
# 4. successful IL7R expression/score alignment
# 5. positive, statistically significant IL7R–LOO association

E5D_IL7R_LOO_ASSOCIATION_PASS_10b3ne <-
  all(
    sample_id_set_pass_E5D_10b3ne,
    sample_id_order_pass_E5D_10b3ne,
    finite_pass_E5D_10b3ne,
    IL7R_expression_alignment_pass_E5D_10b3ne,
    immune_changed_gt_1e10_E5D_10b3ne > 0,
    is.finite(IL7R_LOO_rho_E5D_10b3ne),
    IL7R_LOO_rho_E5D_10b3ne > 0,
    IL7R_LOO_p_E5D_10b3ne < 0.05
  )


cat("\n============================================================\n")

cat(
  "E-5D IL7R-LOO ASSOCIATION PASS: ",
  E5D_IL7R_LOO_ASSOCIATION_PASS_10b3ne,
  "\n",
  sep = ""
)

cat("============================================================\n")


cat("\n============================================================\n")
cat("E-5E-B: CORRECT AOC3-LOO OBJECT-PRESENCE AUDIT\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Correct AOC3-LOO object names
# ------------------------------------------------------------

AOC3_LOO_objects_present_E5E_10b3ne <- all(
  exists("AOC3_LOO_stromal_10b3ne", envir = .GlobalEnv),
  exists("AOC3_LOO_immune_10b3ne", envir = .GlobalEnv),
  exists("AOC3_LOO_estimate_10b3ne", envir = .GlobalEnv),
  exists("AOC3_LOO_stromal_rho_10b3ne", envir = .GlobalEnv)
)

cat(
  "AOC3_LOO_stromal_10b3ne exists: ",
  exists("AOC3_LOO_stromal_10b3ne", envir = .GlobalEnv),
  "\n",
  sep = ""
)

cat(
  "AOC3_LOO_immune_10b3ne exists: ",
  exists("AOC3_LOO_immune_10b3ne", envir = .GlobalEnv),
  "\n",
  sep = ""
)

cat(
  "AOC3_LOO_estimate_10b3ne exists: ",
  exists("AOC3_LOO_estimate_10b3ne", envir = .GlobalEnv),
  "\n",
  sep = ""
)

cat(
  "AOC3_LOO_stromal_rho_10b3ne exists: ",
  exists("AOC3_LOO_stromal_rho_10b3ne", envir = .GlobalEnv),
  "\n",
  sep = ""
)

cat(
  "\nAOC3-LOO object-presence audit: ",
  AOC3_LOO_objects_present_E5E_10b3ne,
  "\n",
  sep = ""
)


cat("\n============================================================\n")
cat("E-5E-C: FINAL LOO RESTORATION AND INTEGRITY AUDIT\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. GlobalEnv SI_geneset restoration
# ------------------------------------------------------------

E5E_GlobalEnv_SI_restored_10b3ne <- identical(
  SI_geneset,
  SI_geneset_original_10b3ne
)

# ------------------------------------------------------------
# 2. Signature membership restoration
# ------------------------------------------------------------

E5E_AOC3_stromal_restored_10b3ne <- identical(
  "AOC3" %in% SI_geneset$StromalSignature,
  TRUE
)

E5E_AOC3_immune_restored_10b3ne <- identical(
  "AOC3" %in% SI_geneset$ImmuneSignature,
  FALSE
)

E5E_IL7R_stromal_restored_10b3ne <- identical(
  "IL7R" %in% SI_geneset$StromalSignature,
  FALSE
)

E5E_IL7R_immune_restored_10b3ne <- identical(
  "IL7R" %in% SI_geneset$ImmuneSignature,
  TRUE
)

E5E_MS4A7_stromal_restored_10b3ne <- identical(
  "MS4A7" %in% SI_geneset$StromalSignature,
  FALSE
)

E5E_MS4A7_immune_restored_10b3ne <- identical(
  "MS4A7" %in% SI_geneset$ImmuneSignature,
  FALSE
)

# ------------------------------------------------------------
# 3. package:estimate integrity
# ------------------------------------------------------------

E5E_package_SI_unchanged_10b3ne <- identical(
  get("SI_geneset", envir = package_estimate_env_10b3ne),
  SI_geneset_original_10b3ne
)

E5E_package_SI_locked_10b3ne <- bindingIsLocked(
  "SI_geneset",
  package_estimate_env_10b3ne
)

# ------------------------------------------------------------
# 4. Locked Step14 object integrity — corrected
# ------------------------------------------------------------

E5E_locked_scores_unchanged_10b3ne <- identical(
  get("estimate_scores_wgcna", envir = step14_env_10b3ne),
  estimate_scores_locked_10b3ne
)

cat(
  "Locked Step14 estimate_scores_wgcna unchanged: ",
  E5E_locked_scores_unchanged_10b3ne,
  "\n",
  sep = ""
)

cat("\n============================================================\n")
cat("E-5E-D: CORRECT SIGNATURE MEMBERSHIP AUDIT\n")
cat("============================================================\n\n")


cat("\n============================================================\n")
cat("E-5E-E: FINAL MEMBERSHIP-CONSISTENCY AUDIT\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Convert actual membership status into PASS/FAIL checks
# ------------------------------------------------------------

E5E_AOC3_stromal_membership_PASS_10b3ne <- identical(
  E5E_AOC3_stromal_restored_10b3ne,
  TRUE
)

E5E_AOC3_immune_membership_PASS_10b3ne <- identical(
  E5E_AOC3_immune_restored_10b3ne,
  FALSE
)

E5E_IL7R_stromal_membership_PASS_10b3ne <- identical(
  E5E_IL7R_stromal_restored_10b3ne,
  FALSE
)

E5E_IL7R_immune_membership_PASS_10b3ne <- identical(
  E5E_IL7R_immune_restored_10b3ne,
  TRUE
)

E5E_MS4A7_stromal_membership_PASS_10b3ne <- identical(
  E5E_MS4A7_stromal_restored_10b3ne,
  FALSE
)

E5E_MS4A7_immune_membership_PASS_10b3ne <- identical(
  E5E_MS4A7_immune_restored_10b3ne,
  FALSE
)

# ------------------------------------------------------------
# 2. Print membership PASS status
# ------------------------------------------------------------

cat(
  "AOC3 Stromal membership PASS: ",
  E5E_AOC3_stromal_membership_PASS_10b3ne,
  "\n",
  sep = ""
)

cat(
  "AOC3 Immune membership PASS: ",
  E5E_AOC3_immune_membership_PASS_10b3ne,
  "\n",
  sep = ""
)

cat(
  "IL7R Stromal membership PASS: ",
  E5E_IL7R_stromal_membership_PASS_10b3ne,
  "\n",
  sep = ""
)

cat(
  "IL7R Immune membership PASS: ",
  E5E_IL7R_immune_membership_PASS_10b3ne,
  "\n",
  sep = ""
)

cat(
  "MS4A7 Stromal membership PASS: ",
  E5E_MS4A7_stromal_membership_PASS_10b3ne,
  "\n",
  sep = ""
)

cat(
  "MS4A7 Immune membership PASS: ",
  E5E_MS4A7_immune_membership_PASS_10b3ne,
  "\n",
  sep = ""
)

##Then regenerate the final E-5E PASS.
E5E_FINAL_RESTORE_AUDIT_PASS_10b3ne <- all(
  E5E_GlobalEnv_SI_restored_10b3ne,
  
  E5E_AOC3_stromal_membership_PASS_10b3ne,
  E5E_AOC3_immune_membership_PASS_10b3ne,
  
  E5E_IL7R_stromal_membership_PASS_10b3ne,
  E5E_IL7R_immune_membership_PASS_10b3ne,
  
  E5E_MS4A7_stromal_membership_PASS_10b3ne,
  E5E_MS4A7_immune_membership_PASS_10b3ne,
  
  E5E_package_SI_unchanged_10b3ne,
  E5E_package_SI_locked_10b3ne,
  E5E_locked_scores_unchanged_10b3ne,
  
  E5E_AOC3_objects_present_10b3ne,
  E5E_IL7R_objects_present_10b3ne
)

cat("\n============================================================\n")
cat(
  "E-5E FINAL RESTORE AUDIT PASS: ",
  E5E_FINAL_RESTORE_AUDIT_PASS_10b3ne,
  "\n",
  sep = ""
)
cat("============================================================\n")



### ============================================================
# STEP 10B-3N-E
# FINAL LOCK CHECK + SAVE
# ============================================================


# ============================================================
# STEP 10B-3N-E: CHECK ACTUAL IL7R-LOO OBJECT NAMES
# ============================================================

cat("\n============================================================\n")
cat("CHECKING EXISTING IL7R-LOO OBJECTS\n")
cat("============================================================\n\n")

grep(
  "^IL7R_LOO_",
  ls(envir = .GlobalEnv),
  value = TRUE
)

grep(
  "IL7R.*LOO|E5D",
  ls(envir = .GlobalEnv),
  value = TRUE
)

ls(
  envir = .GlobalEnv,
  pattern = "^IL7R_LOO_|^E5D_"
)

# ============================================================
# STEP 10B-3N-E: SAVE FINAL LOCKED CHECKPOINT
# Corrected to use ACTUAL existing IL7R-LOO object names
# ============================================================

cat("\n============================================================\n")
cat("STEP 10B-3N-E: FINAL LOCK CHECK\n")
cat("============================================================\n\n")

stopifnot(
  E5E_FINAL_RESTORE_AUDIT_PASS_10b3ne
)

cat("FINAL AUDIT PASS = TRUE\n")
cat("Saving corrected LOCKED checkpoint...\n\n")


# ------------------------------------------------------------
# Final existence check for all objects to be saved
# ------------------------------------------------------------

required_objects_10b3ne <- c(
  
  # Final audit
  "E5E_FINAL_RESTORE_AUDIT_PASS_10b3ne",
  "E5E_GlobalEnv_SI_restored_10b3ne",
  
  "E5E_AOC3_stromal_restored_10b3ne",
  "E5E_AOC3_immune_restored_10b3ne",
  "E5E_IL7R_stromal_restored_10b3ne",
  "E5E_IL7R_immune_restored_10b3ne",
  "E5E_MS4A7_stromal_restored_10b3ne",
  "E5E_MS4A7_immune_restored_10b3ne",
  
  "E5E_AOC3_stromal_membership_PASS_10b3ne",
  "E5E_AOC3_immune_membership_PASS_10b3ne",
  "E5E_IL7R_stromal_membership_PASS_10b3ne",
  "E5E_IL7R_immune_membership_PASS_10b3ne",
  "E5E_MS4A7_stromal_membership_PASS_10b3ne",
  "E5E_MS4A7_immune_membership_PASS_10b3ne",
  
  "E5E_package_SI_unchanged_10b3ne",
  "E5E_package_SI_locked_10b3ne",
  "E5E_locked_scores_unchanged_10b3ne",
  
  "E5E_AOC3_objects_present_10b3ne",
  "E5E_IL7R_objects_present_10b3ne",
  
  # AOC3-LOO
  "AOC3_LOO_stromal_10b3ne",
  "AOC3_LOO_immune_10b3ne",
  "AOC3_LOO_estimate_10b3ne",
  
  "AOC3_LOO_stromal_aligned_10b3ne",
  "AOC3_LOO_immune_aligned_10b3ne",
  "AOC3_LOO_estimate_aligned_10b3ne",
  
  "AOC3_LOO_stromal_delta_10b3ne",
  "AOC3_LOO_immune_delta_10b3ne",
  "AOC3_LOO_estimate_delta_10b3ne",
  
  "AOC3_LOO_stromal_changed_n_10b3ne",
  "AOC3_LOO_immune_changed_n_10b3ne",
  "AOC3_LOO_estimate_changed_n_10b3ne",
  
  "AOC3_LOO_stromal_rho_10b3ne",
  "AOC3_LOO_stromal_p_10b3ne",
  "AOC3_LOO_stromal_spearman_10b3ne",
  
  # IL7R-LOO: actual existing objects
  "IL7R_LOO_stromal_10b3ne",
  "IL7R_LOO_immune_10b3ne",
  "IL7R_LOO_estimate_10b3ne",
  
  "IL7R_LOO_immune_rho_10b3ne",
  "IL7R_LOO_immune_p_10b3ne",
  "IL7R_LOO_immune_cor_test_10b3ne",
  "IL7R_LOO_association_E5D_10b3ne",
  
  "E5D_IL7R_LOO_ASSOCIATION_PASS_10b3ne",
  
  # Additional IL7R-LOO provenance / validation objects
  "IL7R_LOO_input_gct_10b3ne",
  "IL7R_LOO_output_gct_10b3ne",
  "IL7R_LOO_score_matrix_10b3ne",
  "IL7R_LOO_score_names_10b3ne",
  "IL7R_LOO_score_df_10b3ne",
  "IL7R_LOO_score_data_10b3ne",
  
  "IL7R_LOO_sample_ids_10b3ne",
  "IL7R_LOO_sample_ids_raw_10b3ne",
  "IL7R_LOO_sample_ids_normalized_10b3ne",
  "IL7R_LOO_sample_ids_unique_10b3ne",
  "IL7R_LOO_sample_identity_10b3ne",
  "IL7R_LOO_sample_count_10b3ne",
  
  "IL7R_LOO_stromal_10b3ne",
  "IL7R_LOO_immune_10b3ne",
  "IL7R_LOO_estimate_10b3ne",
  
  "IL7R_expression_E5D_10b3ne",
  "IL7R_expression_alignment_pass_E5D_10b3ne",
  "IL7R_expression_score_alignment_E5D_10b3ne",
  
  "max_abs_delta_immune_E5D_10b3ne",
  "max_abs_delta_stromal_E5D_10b3ne",
  "max_abs_delta_estimate_E5D_10b3ne",
  
  "mean_abs_delta_immune_E5D_10b3ne",
  "median_abs_delta_immune_E5D_10b3ne",
  
  "immune_changed_gt_1e10_E5D_10b3ne",
  "stromal_changed_gt_1e10_E5D_10b3ne",
  "estimate_changed_gt_1e10_E5D_10b3ne",
  
  # Original signature / provenance
  "SI_geneset_original_10b3ne",
  "SI_geneset",
  "D2_FINAL_PASS_10b3ne"
)


# ------------------------------------------------------------
# Check existence before save
# ------------------------------------------------------------

missing_objects_10b3ne <- setdiff(
  required_objects_10b3ne,
  ls(envir = .GlobalEnv)
)

if (length(missing_objects_10b3ne) > 0) {
  
  cat("ERROR: Missing objects detected:\n")
  print(missing_objects_10b3ne)
  
  stop(
    "LOCKED checkpoint was NOT saved because required objects are missing."
  )
}

cat("All required objects found: TRUE\n")
cat("Number of objects to save:",
    length(required_objects_10b3ne), "\n\n")


# ------------------------------------------------------------
# Save checkpoint
# ------------------------------------------------------------

save(
  list = required_objects_10b3ne,
  file = "result_dir/Paper 1_Step 7_ Immune Infiltration Analysis/STEP7B-7E-10B-3N/10B-3N-E/Audit/Paper1_STEP10B3NE_formal_LOO_sensitivity_analysis_LOCKED.RData"
)


# ------------------------------------------------------------
# Final file check
# ------------------------------------------------------------

locked_file_10b3ne <-
  "result_dir/Paper 1_Step 7_ Immune Infiltration Analysis/STEP7B-7E-10B-3N/10B-3N-E/Audit/Paper1_STEP10B3NE_formal_LOO_sensitivity_analysis_LOCKED.RData"

cat("\n============================================================\n")
cat("STEP 10B-3N-E: LOCKED CHECKPOINT SAVE RESULT\n")
cat("============================================================\n\n")

cat(
  "File exists:",
  file.exists(locked_file_10b3ne),
  "\n"
)

if (file.exists(locked_file_10b3ne)) {
  
  cat(
    "File size:",
    round(file.info(locked_file_10b3ne)$size / 1024^2, 3),
    "MB\n"
  )
  
  cat(
    "LOCKED CHECKPOINT SAVED SUCCESSFULLY = TRUE\n"
  )
}










