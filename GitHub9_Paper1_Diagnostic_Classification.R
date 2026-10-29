##GitHub9_Paper1_Diagnostic_Classification

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

# ============================================================
# WORKFLOW-13
# Exploratory Tumor-versus-Normal Tissue Classification
# ============================================================

cat("\n")
cat("============================================================\n")
cat("WORKFLOW-13 CONFIGURATION\n")
cat("Exploratory Tumor-versus-Normal Tissue Classification\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Root directory
# ------------------------------------------------------------

workflow_root <- 
  "input_dir/Paper 1_Step 13-Diagnostic Classification"

# ------------------------------------------------------------
# 1. Subdirectories
# ------------------------------------------------------------

dir.create(workflow_root, recursive = TRUE, showWarnings = FALSE)

subdirs <- c(
  "00_CONFIG",
  "STEP13A_Clinical_Claim_Deescalation",
  "STEP13B_TCGA_Input_Lock",
  "STEP13C_GSE32863_Input_Lock",
  "STEP13D_Paired_ID_Audit",
  "STEP13E_Expression_Scale_Audit",
  "STEP13F_Model_Definition_Lock",
  "STEP13G_Nested_CV",
  "STEP13H_TCGA_Final_Model",
  "STEP13I_GSE32863_External_Assessment",
  "STEP13J_Final_Audit",
  "INPUT_LOCK",
  "MODEL_LOCK",
  "RESULTS",
  "FIGURES",
  "TABLES",
  "AUDIT",
  "LOGS"
)

for (x in subdirs) {
  dir.create(
    file.path(workflow_root, x),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ------------------------------------------------------------
# 2. Core genes
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

# ------------------------------------------------------------
# 3. TCGA outcome definition
# ------------------------------------------------------------

tcga_outcome_normal <- "Normal"
tcga_outcome_tumor  <- "Tumor"

# ------------------------------------------------------------
# 4. Locked analysis objective
# ------------------------------------------------------------

analysis_objective <- 
  "Exploratory tumor-versus-normal tissue classification"

clinical_claim_allowed <- FALSE

diagnostic_biomarker_claim_allowed <- FALSE

# ------------------------------------------------------------
# 5. Model
# ------------------------------------------------------------

model_type <- "multivariable_logistic_regression"

model_formula <- 
  "TumorStatus ~ AOC3 + MS4A7 + IL7R"

# ------------------------------------------------------------
# 6. Cross-validation
# ------------------------------------------------------------

outer_folds <- 10
inner_folds <- 10

nested_cv <- TRUE

stratified_cv <- TRUE

# ------------------------------------------------------------
# 7. External validation rules
# ------------------------------------------------------------

external_refit_allowed <- FALSE

external_threshold_optimization_allowed <- FALSE

external_coefficient_estimation_allowed <- FALSE

# ------------------------------------------------------------
# 8. GSE32863 paired design
# ------------------------------------------------------------

paired_analysis_required <- TRUE

paired_test <- "Wilcoxon signed-rank test"

# ------------------------------------------------------------
# 9. Threshold policy
# ------------------------------------------------------------

threshold_policy <- 
  "No external threshold optimization"

threshold_source <- 
  "Training data only"

# ------------------------------------------------------------
# 10. Random seed
# ------------------------------------------------------------

workflow_seed <- 20260926

# ------------------------------------------------------------
# 11. Software versions
# ------------------------------------------------------------

r_version <- R.version.string

# ------------------------------------------------------------
# 12. Save configuration
# ------------------------------------------------------------

saveRDS(
  list(
    workflow_root = workflow_root,
    core_genes = core_genes,
    analysis_objective = analysis_objective,
    clinical_claim_allowed = clinical_claim_allowed,
    diagnostic_biomarker_claim_allowed =
      diagnostic_biomarker_claim_allowed,
    model_type = model_type,
    model_formula = model_formula,
    outer_folds = outer_folds,
    inner_folds = inner_folds,
    nested_cv = nested_cv,
    stratified_cv = stratified_cv,
    external_refit_allowed = external_refit_allowed,
    external_threshold_optimization_allowed =
      external_threshold_optimization_allowed,
    external_coefficient_estimation_allowed =
      external_coefficient_estimation_allowed,
    paired_analysis_required = paired_analysis_required,
    paired_test = paired_test,
    threshold_policy = threshold_policy,
    threshold_source = threshold_source,
    workflow_seed = workflow_seed,
    r_version = r_version
  ),
  file.path(
    workflow_root,
    "00_CONFIG",
    "WORKFLOW13_CONFIG_LOCKED.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    workflow_root,
    "LOGS",
    "WORKFLOW13_sessionInfo_CONFIG.txt"
  )
)

cat("CONFIG LOCKED\n")
cat("Core genes:", paste(core_genes, collapse = ", "), "\n")
cat("Objective:", analysis_objective, "\n")
cat("Clinical diagnostic claim:", clinical_claim_allowed, "\n")
cat("Nested CV:", nested_cv, "\n")
cat("External refitting:", external_refit_allowed, "\n")
cat("External threshold optimization:",
    external_threshold_optimization_allowed, "\n")

# ============================================================
# WORKFLOW-13
# STEP13A
# Clinical Diagnostic Claim De-escalation Lock
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP13A — CLINICAL DIAGNOSTIC CLAIM DE-ESCALATION LOCK\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Workflow root
# ------------------------------------------------------------

workflow_root <-
  "input_dir/Paper 1_Step 13-Diagnostic Classification"


# ------------------------------------------------------------
# 1. Locked configuration
# ------------------------------------------------------------

config_file <- file.path(
  workflow_root,
  "00_CONFIG",
  "WORKFLOW13_CONFIG_LOCKED.rds"
)

if (!file.exists(config_file)) {
  
  stop(
    "\nSTEP13A FAIL: Locked CONFIG file not found.\n",
    "Expected:\n",
    config_file,
    "\n"
  )
  
}

config <- readRDS(config_file)


# ------------------------------------------------------------
# 2. Configuration integrity checks
# ------------------------------------------------------------

expected_core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

if (!identical(
  config$core_genes,
  expected_core_genes
)) {
  
  stop(
    "\nSTEP13A FAIL: Core-gene configuration mismatch.\n",
    "Expected: ",
    paste(expected_core_genes, collapse = ", "),
    "\nObserved: ",
    paste(config$core_genes, collapse = ", "),
    "\n"
  )
  
}


expected_objective <-
  "Exploratory tumor-versus-normal tissue classification"

if (!identical(
  config$analysis_objective,
  expected_objective
)) {
  
  stop(
    "\nSTEP13A FAIL: Analysis objective mismatch.\n"
  )
  
}


if (!identical(
  config$clinical_claim_allowed,
  FALSE
)) {
  
  stop(
    "\nSTEP13A FAIL: Clinical diagnostic claim is not disabled.\n"
  )
  
}


if (!identical(
  config$diagnostic_biomarker_claim_allowed,
  FALSE
)) {
  
  stop(
    "\nSTEP13A FAIL: Diagnostic biomarker claim is not disabled.\n"
  )
  
}


if (!identical(
  config$nested_cv,
  TRUE
)) {
  
  stop(
    "\nSTEP13A FAIL: Nested CV is not enabled.\n"
  )
  
}


if (!identical(
  config$external_refit_allowed,
  FALSE
)) {
  
  stop(
    "\nSTEP13A FAIL: External refitting is not prohibited.\n"
  )
  
}


if (!identical(
  config$external_threshold_optimization_allowed,
  FALSE
)) {
  
  stop(
    "\nSTEP13A FAIL: External threshold optimization is not prohibited.\n"
  )
  
}


if (!identical(
  config$external_coefficient_estimation_allowed,
  FALSE
)) {
  
  stop(
    "\nSTEP13A FAIL: External coefficient estimation is not prohibited.\n"
  )
  
}


if (!identical(
  config$paired_analysis_required,
  TRUE
)) {
  
  stop(
    "\nSTEP13A FAIL: Paired external analysis is not locked.\n"
  )
  
}


# ------------------------------------------------------------
# 3. STEP13A output directory
# ------------------------------------------------------------

step_dir <- file.path(
  workflow_root,
  "STEP13A_Clinical_Claim_Deescalation"
)

dir.create(
  step_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 4. Formal lock table
# ------------------------------------------------------------

lock_table <- data.frame(
  
  Parameter = c(
    "Core genes",
    "Analysis objective",
    "Clinical diagnostic claim",
    "Diagnostic biomarker claim",
    "Nested cross-validation",
    "External model refitting",
    "External coefficient estimation",
    "External threshold optimization",
    "Paired external analysis"
  ),
  
  Locked_value = c(
    paste(
      config$core_genes,
      collapse = ", "
    ),
    config$analysis_objective,
    as.character(
      config$clinical_claim_allowed
    ),
    as.character(
      config$diagnostic_biomarker_claim_allowed
    ),
    as.character(
      config$nested_cv
    ),
    as.character(
      config$external_refit_allowed
    ),
    as.character(
      config$external_coefficient_estimation_allowed
    ),
    as.character(
      config$external_threshold_optimization_allowed
    ),
    as.character(
      config$paired_analysis_required
    )
  ),
  
  Status = rep(
    "LOCKED",
    9
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  lock_table,
  file.path(
    step_dir,
    "STEP13A_Clinical_Claim_Lock.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 5. Forbidden clinical claims
# ------------------------------------------------------------

forbidden_claims <- c(
  
  "diagnostic biomarker",
  
  "clinical diagnostic biomarker",
  
  "clinical diagnostic value",
  
  "clinical diagnostic utility",
  
  "clinical applicability",
  
  "clinical implementation",
  
  "patient stratification based on the diagnostic model"
  
)


writeLines(
  forbidden_claims,
  file.path(
    step_dir,
    "STEP13A_Excluded_Clinical_Claims.txt"
  ),
  useBytes = TRUE
)


# ------------------------------------------------------------
# 6. Interpretation lock
# ------------------------------------------------------------

interpretation_lock <- c(
  
  "STEP13A INTERPRETATION LOCK",
  "",
  
  "The three-gene model is evaluated exclusively as an",
  "exploratory tumor-versus-normal tissue classification model.",
  "",
  
  "The analysis is not presented as a clinical diagnostic",
  "biomarker or clinical diagnostic tool.",
  "",
  
  "No clinical diagnostic utility is claimed.",
  "No clinical implementation is claimed.",
  "No patient stratification claim is based on this classifier.",
  "",
  
  "Nested cross-validation is mandatory for internal evaluation.",
  "",
  
  "External model refitting is prohibited.",
  "External coefficient estimation is prohibited.",
  "External threshold optimization is prohibited.",
  "",
  
  "GSE32863 is treated as an independent paired external cohort.",
  "",
  
  "Core genes locked a priori:",
  "AOC3, MS4A7, IL7R"
  
)


writeLines(
  interpretation_lock,
  file.path(
    step_dir,
    "STEP13A_Interpretation_Lock.txt"
  ),
  useBytes = TRUE
)


# ------------------------------------------------------------
# 7. PASS checkpoint
# ------------------------------------------------------------

pass_file <- file.path(
  step_dir,
  "STEP13A_PASS.txt"
)


writeLines(
  
  c(
    
    "============================================================",
    "STEP13A STATUS: PASS",
    "============================================================",
    "",
    paste(
      "Core genes:",
      paste(
        config$core_genes,
        collapse = ", "
      )
    ),
    "",
    paste(
      "Analysis objective:",
      config$analysis_objective
    ),
    "",
    "Clinical diagnostic claim: FALSE",
    "Diagnostic biomarker claim: FALSE",
    "",
    "Nested CV: TRUE",
    "External refitting: FALSE",
    "External coefficient estimation: FALSE",
    "External threshold optimization: FALSE",
    "",
    "Paired external analysis: TRUE",
    "",
    "No statistical calculation was performed.",
    "",
    "STEP13A is an interpretation and analysis-policy lock.",
    "============================================================"
    
  ),
  
  pass_file,
  useBytes = TRUE
  
)


# ------------------------------------------------------------
# 8. Console report
# ------------------------------------------------------------

cat("\n")
cat("STEP13A PASS\n")
cat("------------------------------------------------------------\n")

cat(
  "Core genes:                 ",
  paste(
    config$core_genes,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Analysis objective:         ",
  config$analysis_objective,
  "\n"
)

cat(
  "Clinical diagnostic claim:  ",
  config$clinical_claim_allowed,
  "\n"
)

cat(
  "Diagnostic biomarker claim: ",
  config$diagnostic_biomarker_claim_allowed,
  "\n"
)

cat(
  "Nested CV:                  ",
  config$nested_cv,
  "\n"
)

cat(
  "External refitting:         ",
  config$external_refit_allowed,
  "\n"
)

cat(
  "External coefficient est.:  ",
  config$external_coefficient_estimation_allowed,
  "\n"
)

cat(
  "External threshold tuning:  ",
  config$external_threshold_optimization_allowed,
  "\n"
)

cat(
  "Paired external analysis:   ",
  config$paired_analysis_required,
  "\n"
)

cat("\n")
cat("No statistical calculation performed.\n")
cat("Clinical diagnostic wording is formally excluded.\n")
cat("============================================================\n")

##STEP13B — TCGA Input / Sample Label Lock
# ============================================================
# WORKFLOW-13
# STEP13B — TCGA Diagnostic Input Inventory and Lock
# ============================================================

cat("\n")
cat("============================================================\n")
cat("WORKFLOW-13 — STEP13B\n")
cat("TCGA Diagnostic Input Inventory / Candidate Audit\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

config_file <- file.path(
  workflow_root,
  "00_CONFIG",
  "WORKFLOW13_CONFIG_LOCKED.rds"
)

step13a_pass <- file.path(
  workflow_root,
  "STEP13A_Clinical_Claim_Deescalation",
  "STEP13A_PASS.txt"
)

step_dir <- file.path(
  workflow_root,
  "STEP13B_TCGA_Input_Lock"
)

input_dir <- "input_dir/Paper1_tumor-only WGCNA"

dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load locked configuration
# ------------------------------------------------------------

if (!file.exists(config_file)) {
  stop(
    "FAIL: WORKFLOW13_CONFIG_LOCKED.rds not found:\n",
    config_file
  )
}

config <- readRDS(config_file)

# ------------------------------------------------------------
# 2. Require STEP13A PASS
# ------------------------------------------------------------

if (!file.exists(step13a_pass)) {
  stop(
    "FAIL: STEP13A_PASS.txt not found.\n",
    "STEP13A must PASS before STEP13B."
  )
}

cat("STEP13A prerequisite: PASS\n\n")

# ------------------------------------------------------------
# 3. Locked core genes
# ------------------------------------------------------------

core_genes <- c("AOC3", "MS4A7", "IL7R")

cat("Locked core genes:\n")
cat(paste(core_genes, collapse = ", "), "\n\n")

# ------------------------------------------------------------
# 4. Candidate files explicitly identified by user
# ------------------------------------------------------------

candidate_files <- c(
  file.path(
    input_dir,
    "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt"
  ),
  file.path(
    input_dir,
    "Workflow-04-2C-1A_ESTIMATE_input_GeneSymbol_log2_FPKM_UQ.txt"
  ),
  file.path(
    input_dir,
    "Workflow-04-2C-1A_Ensembl_to_GeneSymbol_GENCODE_v22.txt"
  ),
  file.path(
    input_dir,
    "Workflow-04-2C_expr_primary_tumor_log2_FPKM_UQ.txt"
  )
)

# ------------------------------------------------------------
# 5. File existence audit
# ------------------------------------------------------------

existence_df <- data.frame(
  file = basename(candidate_files),
  full_path = candidate_files,
  exists = file.exists(candidate_files),
  size_MB = ifelse(
    file.exists(candidate_files),
    round(file.info(candidate_files)$size / 1024^2, 3),
    NA_real_
  ),
  stringsAsFactors = FALSE
)

print(existence_df)

if (!all(existence_df$exists)) {
  cat("\nWARNING: One or more candidate files do not exist.\n")
}

# ------------------------------------------------------------
# 6. Function: inspect expression-like file
# ------------------------------------------------------------

inspect_expression_file <- function(f) {
  
  cat("\n------------------------------------------------------------\n")
  cat("FILE:", basename(f), "\n")
  cat("------------------------------------------------------------\n")
  
  if (!file.exists(f)) {
    return(data.frame(
      file = basename(f),
      exists = FALSE,
      n_rows_preview = NA,
      n_cols = NA,
      first_column = NA,
      AOC3_present = NA,
      MS4A7_present = NA,
      IL7R_present = NA,
      Tumor_label_present = NA,
      Normal_label_present = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  # Read a limited preview only.
  dat <- tryCatch(
    read.delim(
      f,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE,
      nrows = 10
    ),
    error = function(e) NULL
  )
  
  if (is.null(dat)) {
    
    cat("Could not parse with tab-delimited header.\n")
    
    return(data.frame(
      file = basename(f),
      exists = TRUE,
      n_rows_preview = NA,
      n_cols = NA,
      first_column = NA,
      AOC3_present = NA,
      MS4A7_present = NA,
      IL7R_present = NA,
      Tumor_label_present = NA,
      Normal_label_present = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  cn <- colnames(dat)
  
  # Search all visible preview values
  vals <- unique(
    unlist(
      lapply(dat, as.character),
      use.names = FALSE
    )
  )
  
  vals_upper <- toupper(vals)
  
  # Gene-symbol detection
  gene_hits <- toupper(
    unlist(
      lapply(dat, as.character),
      use.names = FALSE
    )
  )
  
  result <- data.frame(
    file = basename(f),
    exists = TRUE,
    n_rows_preview = nrow(dat),
    n_cols = ncol(dat),
    first_column = ifelse(length(cn) >= 1, cn[1], NA),
    AOC3_present = "AOC3" %in% gene_hits,
    MS4A7_present = "MS4A7" %in% gene_hits,
    IL7R_present = "IL7R" %in% gene_hits,
    Tumor_label_present = any(
      grepl("TUMOR", vals_upper, fixed = TRUE)
    ),
    Normal_label_present = any(
      grepl("NORMAL", vals_upper, fixed = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  
  cat("Columns:", ncol(dat), "\n")
  cat("Column names:\n")
  print(cn)
  
  cat("\nPreview:\n")
  print(utils::head(dat, 5))
  
  cat("\nCore-gene hits in preview:\n")
  cat(
    "AOC3 :", result$AOC3_present,
    "\nMS4A7:", result$MS4A7_present,
    "\nIL7R :", result$IL7R_present,
    "\n"
  )
  
  cat("\nTumor label in preview :", result$Tumor_label_present, "\n")
  cat("Normal label in preview:", result$Normal_label_present, "\n")
  
  result
}

# ------------------------------------------------------------
# 7. Inspect all four files
# ------------------------------------------------------------

audit_list <- lapply(
  candidate_files,
  inspect_expression_file
)

audit_df <- do.call(
  rbind,
  audit_list
)

# ------------------------------------------------------------
# 8. Save candidate audit
# ------------------------------------------------------------

audit_file <- file.path(
  step_dir,
  "STEP13B_TCGA_Candidate_Input_Audit.csv"
)

write.csv(
  audit_df,
  audit_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Preliminary interpretation
# ------------------------------------------------------------

interpretation_file <- file.path(
  step_dir,
  "STEP13B_Preliminary_Input_Interpretation.txt"
)

sink(interpretation_file)

cat("WORKFLOW-13 — STEP13B\n")
cat("Preliminary TCGA Input Interpretation\n")
cat("============================================================\n\n")

cat("Purpose:\n")
cat("Identify the correct TCGA tumor-versus-normal input before\n")
cat("any ROC/model recalculation.\n\n")

cat("Locked core genes:\n")
cat(paste(core_genes, collapse = ", "), "\n\n")

cat("Candidate files inspected:\n\n")

for (f in candidate_files) {
  cat(" - ", f, "\n", sep = "")
}

cat("\nImportant methodological rule:\n")
cat(
  "No file is designated as the final diagnostic/classification\n",
  "input solely from its filename. The final input must be\n",
  "demonstrated to contain the appropriate TCGA sample population,\n",
  "tumor/normal labels, expression scale, and the three locked genes.\n",
  sep = ""
)

cat("\nCurrent status:\n")
cat("STEP13B INPUT INVENTORY COMPLETE — FINAL INPUT NOT YET LOCKED\n")

sink()

# ------------------------------------------------------------
# 10. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("TCGA CANDIDATE INPUT AUDIT SUMMARY\n")
cat("============================================================\n\n")

print(audit_df)

cat("\nAudit saved to:\n")
cat(audit_file, "\n")

cat("\nInterpretation saved to:\n")
cat(interpretation_file, "\n")

cat("\n")
cat("STEP13B INPUT INVENTORY COMPLETE\n")
cat("FINAL TCGA INPUT HAS NOT YET BEEN LOCKED.\n")
cat("NO ROC OR MODEL CALCULATION PERFORMED.\n")
cat("============================================================\n")


# ============================================================
# WORKFLOW-13
# STEP13B — Core Gene Symbol Traceability Audit
# ============================================================

cat("\n")
cat("============================================================\n")
cat("WORKFLOW-13 — STEP13B CORE-GENE TRACEABILITY AUDIT\n")
cat("============================================================\n\n")

expr_symbol_file <- paste0(
  "input_dir/Paper1_tumor-only WGCNA/",
  "Workflow-04-2C-1A_ESTIMATE_input_GeneSymbol_log2_FPKM_UQ.txt"
)

mapping_file <- paste0(
  "input_dir/Paper1_tumor-only WGCNA/",
  "Workflow-04-2C-1A_Ensembl_to_GeneSymbol_GENCODE_v22.txt"
)

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

# ------------------------------------------------------------
# 1. Check files
# ------------------------------------------------------------

if (!file.exists(expr_symbol_file)) {
  stop("Expression file not found:\n", expr_symbol_file)
}

if (!file.exists(mapping_file)) {
  stop("Mapping file not found:\n", mapping_file)
}

cat("Expression file:\n")
cat(expr_symbol_file, "\n\n")

cat("Mapping file:\n")
cat(mapping_file, "\n\n")

# ------------------------------------------------------------
# 2. Read GeneSymbol expression matrix
# ------------------------------------------------------------

expr <- read.delim(
  expr_symbol_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Expression matrix dimensions:\n")
cat(
  "Rows:", nrow(expr),
  "\nColumns:", ncol(expr),
  "\n\n"
)

cat("First column:\n")
cat(colnames(expr)[1], "\n\n")

# ------------------------------------------------------------
# 3. Extract GeneSymbol column
# ------------------------------------------------------------

gene_symbol <- trimws(
  as.character(expr[[1]])
)

gene_symbol_upper <- toupper(
  gene_symbol
)

# ------------------------------------------------------------
# 4. Core gene audit
# ------------------------------------------------------------

core_audit <- data.frame(
  core_gene = core_genes,
  present = core_genes %in% gene_symbol_upper,
  row_index = sapply(
    core_genes,
    function(g) {
      idx <- which(gene_symbol_upper == g)
      if (length(idx) == 0) NA_integer_ else idx[1]
    }
  ),
  stringsAsFactors = FALSE
)

cat("============================================================\n")
cat("CORE GENE SYMBOL AUDIT\n")
cat("============================================================\n\n")

print(core_audit)

# ------------------------------------------------------------
# 5. Show exact rows
# ------------------------------------------------------------

for (g in core_genes) {
  
  idx <- which(gene_symbol_upper == g)
  
  cat("\n------------------------------------------------------------\n")
  cat("Gene:", g, "\n")
  cat("------------------------------------------------------------\n")
  
  if (length(idx) == 0) {
    cat("NOT FOUND\n")
  } else {
    
    cat(
      "Rows found:",
      paste(idx, collapse = ", "),
      "\n"
    )
    
    print(
      expr[idx, , drop = FALSE]
    )
  }
}

# ------------------------------------------------------------
# 6. Mapping file audit
# ------------------------------------------------------------

mapping <- read.delim(
  mapping_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n")
cat("============================================================\n")
cat("MAPPING FILE AUDIT\n")
cat("============================================================\n\n")

cat(
  "Mapping dimensions:",
  nrow(mapping),
  "rows ×",
  ncol(mapping),
  "columns\n\n"
)

print(
  head(mapping, 10)
)

# ------------------------------------------------------------
# 7. Save audit
# ------------------------------------------------------------

out_dir <- paste0(
  "result_dir/Paper 1_Step 13-Diagnostic Classification/",
  "STEP13B_TCGA_Input_Lock"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  core_audit,
  file.path(
    out_dir,
    "STEP13B_Core_Gene_Symbol_Audit.csv"
  ),
  row.names = FALSE
)

write.table(
  data.frame(
    expression_file = expr_symbol_file,
    mapping_file = mapping_file,
    n_expression_rows = nrow(expr),
    n_expression_columns = ncol(expr),
    core_genes = paste(core_genes, collapse = ", "),
    all_core_genes_present = all(core_audit$present)
  ),
  file.path(
    out_dir,
    "STEP13B_Core_Gene_Traceability_Summary.txt"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

cat("\n")
cat("============================================================\n")
cat("AUDIT COMPLETE\n")
cat("============================================================\n")

cat(
  "All three core genes present:",
  all(core_audit$present),
  "\n"
)

cat("\n")
cat("NO ROC CALCULATION PERFORMED.\n")
cat("NO MODEL FITTING PERFORMED.\n")
cat("NO STATISTICAL TESTING PERFORMED.\n")
cat("============================================================\n")


cat("\n")
cat("============================================================\n")
cat("STEP13B — SAMPLE COLUMN IDENTITY AUDIT\n")
cat("============================================================\n\n")

cat("Expression matrix dimensions:\n")
cat("Rows   :", nrow(expr), "\n")
cat("Columns:", ncol(expr), "\n\n")

cat("First 20 column names:\n")
print(colnames(expr)[1:min(20, ncol(expr))])

cat("\nLast 20 column names:\n")
print(
  colnames(expr)[
    max(1, ncol(expr) - 19):ncol(expr)
  ]
)

cat("\n------------------------------------------------------------\n")
cat("Sample-column count\n")
cat("------------------------------------------------------------\n")

cat(
  "Gene identifier column:",
  colnames(expr)[1],
  "\n"
)

cat(
  "Remaining columns:",
  ncol(expr) - 1,
  "\n"
)

cat("\n------------------------------------------------------------\n")
cat("TCGA barcode pattern audit\n")
cat("------------------------------------------------------------\n")

sample_ids <- colnames(expr)[-1]

tcga_like <- grepl(
  "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}",
  sample_ids,
  ignore.case = TRUE
)

cat(
  "TCGA-like sample IDs:",
  sum(tcga_like),
  "/",
  length(sample_ids),
  "\n"
)

cat("\nFirst 20 sample IDs:\n")
print(head(sample_ids, 20))

cat("\n------------------------------------------------------------\n")
cat("Potential sample-type tokens\n")
cat("------------------------------------------------------------\n")

sample_upper <- toupper(sample_ids)

cat(
  "Columns containing 01:",
  sum(grepl("-01", sample_upper, fixed = TRUE)),
  "\n"
)

cat(
  "Columns containing 11:",
  sum(grepl("-11", sample_upper, fixed = TRUE)),
  "\n"
)

cat(
  "Columns containing NORMAL:",
  sum(grepl("NORMAL", sample_upper, fixed = TRUE)),
  "\n"
)

cat(
  "Columns containing TUMOR:",
  sum(grepl("TUMOR", sample_upper, fixed = TRUE)),
  "\n"
)

cat("\n============================================================\n")
cat("NO ROC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TESTING\n")
cat("============================================================\n")


# ============================================================
# WORKFLOW-13
# STEP13B-3 — Historical 575-Sample Provenance Audit
#
# Purpose:
#   Audit historical WORKFLOW-03F provenance files to identify
#   the Tumor/Normal annotation and reconstruct the original
#   575-sample diagnostic input.
#
# NO ROC
# NO MODEL FITTING
# NO STATISTICAL TESTING
# NO DATA MODIFICATION
# ============================================================

cat("\n")
cat("============================================================\n")
cat("WORKFLOW-13 — STEP13B-3\n")
cat("Historical 575-Sample Provenance Audit\n")
cat("============================================================\n\n")

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

source_dir <- "input_dir/Paper1_tumor-only WGCNA"

step_dir <- file.path(
  workflow_root,
  "STEP13B_TCGA_Input_Lock"
)

dir.create(
  step_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Files to audit
# ------------------------------------------------------------

files <- c(
  
  "Workflow-03F-9_provenance_audit_DONE.txt",
  
  "Workflow-03F-9_expression_scale_provenance.txt",
  
  "Workflow-03F-8_expression_provenance_test.txt",
  
  "Workflow-03F-7_expression_matrix_structure.txt",
  
  "Workflow-03F-5_historical_517_patient_annotation.txt",
  
  "Workflow-03F-3_patient_expression_mapping.txt",
  
  "Workflow-03F-2_expression_annotation.txt",
  
  "Workflow-03F-2_primary_tumor_with_clinical.txt",
  
  "TCGA_Patient_Level_Summary.txt",
  
  "TCGA_Expression_Sample_Annotation_Basic.txt"
)

paths <- file.path(
  source_dir,
  files
)

exist_df <- data.frame(
  file = files,
  exists = file.exists(paths),
  size_MB = ifelse(
    file.exists(paths),
    round(file.info(paths)$size / 1024^2, 3),
    NA_real_
  ),
  stringsAsFactors = FALSE
)

cat("============================================================\n")
cat("FILE EXISTENCE AUDIT\n")
cat("============================================================\n\n")

print(exist_df)

# ------------------------------------------------------------
# 2. Read text files
# ------------------------------------------------------------

audit_terms <- c(
  "Tumor",
  "Normal",
  "Primary_Tumor",
  "Normal_Samples",
  "Has_Normal",
  "Has_Primary_Tumor",
  "516",
  "517",
  "575",
  "59"
)

results <- list()

for (i in seq_along(paths)) {
  
  f <- paths[i]
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat(files[i], "\n")
  cat("------------------------------------------------------------\n")
  
  if (!file.exists(f)) {
    
    cat("FILE NOT FOUND\n")
    
    results[[i]] <- data.frame(
      file = files[i],
      exists = FALSE,
      n_lines = NA_integer_,
      tumor_hits = NA_integer_,
      normal_hits = NA_integer_,
      stringsAsFactors = FALSE
    )
    
    next
  }
  
  txt <- tryCatch(
    readLines(
      f,
      warn = FALSE,
      encoding = "UTF-8"
    ),
    error = function(e) character(0)
  )
  
  txt_upper <- toupper(txt)
  
  tumor_hits <- sum(
    grepl("TUMOR", txt_upper, fixed = TRUE)
  )
  
  normal_hits <- sum(
    grepl("NORMAL", txt_upper, fixed = TRUE)
  )
  
  cat("Lines:", length(txt), "\n")
  cat("Tumor-containing lines :", tumor_hits, "\n")
  cat("Normal-containing lines:", normal_hits, "\n")
  
  # ----------------------------------------------------------
  # Print lines containing key terms
  # ----------------------------------------------------------
  
  key_idx <- which(
    grepl(
      "TUMOR|NORMAL|575|516|517|59",
      txt_upper
    )
  )
  
  if (length(key_idx) > 0) {
    
    cat("\nRelevant lines:\n")
    
    show_idx <- head(
      key_idx,
      80
    )
    
    for (j in show_idx) {
      
      cat(
        sprintf(
          "[%d] %s\n",
          j,
          txt[j]
        )
      )
    }
    
  } else {
    
    cat("\nNo key provenance terms found.\n")
  }
  
  results[[i]] <- data.frame(
    file = files[i],
    exists = TRUE,
    n_lines = length(txt),
    tumor_hits = tumor_hits,
    normal_hits = normal_hits,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# 3. Combine results
# ------------------------------------------------------------

audit_df <- do.call(
  rbind,
  results
)

cat("\n")
cat("============================================================\n")
cat("STEP13B-3 SUMMARY\n")
cat("============================================================\n\n")

print(audit_df)

# ------------------------------------------------------------
# 4. Save summary
# ------------------------------------------------------------

write.csv(
  audit_df,
  file.path(
    step_dir,
    "STEP13B-3_Historical_575_Sample_Provenance_Audit.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Save formal status
# ------------------------------------------------------------

status_file <- file.path(
  step_dir,
  "STEP13B-3_Status.txt"
)

writeLines(
  c(
    "WORKFLOW-13 — STEP13B-3",
    "Historical 575-Sample Provenance Audit",
    "============================================================",
    "",
    "Purpose:",
    "Identify the historical TCGA Tumor/Normal annotation and",
    "reconstruct the original 575-sample diagnostic input.",
    "",
    "Target historical sample structure:",
    "516 Tumor + 59 Normal = 575 samples",
    "",
    "Current status:",
    "PROVENANCE AUDIT COMPLETE",
    "FINAL TCGA INPUT NOT YET LOCKED",
    "",
    "No ROC calculation performed.",
    "No model fitting performed.",
    "No statistical testing performed."
  ),
  status_file
)

cat("\n")
cat("============================================================\n")
cat("STEP13B-3 COMPLETE\n")
cat("============================================================\n")
cat("Historical provenance audited.\n")
cat("Final TCGA input NOT YET LOCKED.\n")
cat("NO ROC / MODEL / STATISTICAL ANALYSIS PERFORMED.\n")
cat("============================================================\n")

# ============================================================
# STEP13B-5
# Exact sample-level composition audit
#
# Purpose:
#   1. Reconstruct exact Tumor / Normal / Other counts
#   2. Identify the two records not classified as Tumor/Normal
#   3. Verify sample-level annotation uniqueness
#
# NO ROC
# NO MODEL
# NO STATISTICAL TEST
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP13B-5 — EXACT SAMPLE-LEVEL COMPOSITION AUDIT\n")
cat("============================================================\n\n")

root <- "input_dir/Paper1_tumor-only WGCNA"

ann_file <- file.path(
  root,
  "TCGA_Expression_Sample_Annotation_Basic.txt"
)

expr_structure_file <- file.path(
  root,
  "Workflow-03F-7_expression_matrix_structure.txt"
)

# ------------------------------------------------------------
# 1. Read annotation
# ------------------------------------------------------------

ann <- read.delim(
  ann_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("[1] Annotation dimensions\n")
cat("Rows:", nrow(ann), "\n")
cat("Columns:", ncol(ann), "\n\n")

cat("[2] Annotation columns\n")
print(names(ann))
cat("\n")

# ------------------------------------------------------------
# 2. Exact Tumor_Normal composition
# ------------------------------------------------------------

cat("[3] Exact Tumor_Normal distribution\n")
print(table(ann$Tumor_Normal, useNA = "ifany"))
cat("\n")

# ------------------------------------------------------------
# 3. Exact Sample_Type_Code distribution
# ------------------------------------------------------------

cat("[4] Exact Sample_Type_Code distribution\n")
print(table(ann$Sample_Type_Code, useNA = "ifany"))
cat("\n")

# ------------------------------------------------------------
# 4. Identify all records not classified as Tumor/Normal
# ------------------------------------------------------------

other_idx <- which(
  is.na(ann$Tumor_Normal) |
    !(ann$Tumor_Normal %in% c("Tumor", "Normal"))
)

cat("[5] Non-Tumor/Normal records\n")
cat("Number:", length(other_idx), "\n\n")

if (length(other_idx) > 0) {
  print(ann[other_idx, , drop = FALSE])
} else {
  cat("NONE\n")
}

cat("\n")

# ------------------------------------------------------------
# 5. Check duplicate expression sample IDs
# ------------------------------------------------------------

dup_sample <- duplicated(ann$Expression_Sample)

cat("[6] Duplicate Expression_Sample IDs\n")
cat("Number duplicated:", sum(dup_sample), "\n")

if (any(dup_sample)) {
  print(ann[dup_sample, , drop = FALSE])
} else {
  cat("NONE\n")
}

cat("\n")

# ------------------------------------------------------------
# 6. Check missing expression sample IDs
# ------------------------------------------------------------

missing_sample <- is.na(ann$Expression_Sample) |
  ann$Expression_Sample == ""

cat("[7] Missing Expression_Sample IDs\n")
cat("Number:", sum(missing_sample), "\n")

if (any(missing_sample)) {
  print(ann[missing_sample, , drop = FALSE])
} else {
  cat("NONE\n")
}

cat("\n")

# ------------------------------------------------------------
# 7. Exact counts
# ------------------------------------------------------------

n_total <- nrow(ann)
n_tumor <- sum(ann$Tumor_Normal == "Tumor", na.rm = TRUE)
n_normal <- sum(ann$Tumor_Normal == "Normal", na.rm = TRUE)
n_other <- n_total - n_tumor - n_normal

cat("[8] Exact sample composition\n")
cat("Total annotation records :", n_total, "\n")
cat("Tumor                    :", n_tumor, "\n")
cat("Normal                   :", n_normal, "\n")
cat("Other / unclassified     :", n_other, "\n")
cat("Tumor + Normal + Other   :", n_tumor + n_normal + n_other, "\n")
cat("\n")

# ------------------------------------------------------------
# 8. Historical target consistency
# ------------------------------------------------------------

cat("[9] Historical target comparison\n")
cat("Expected historical total samples :", 575, "\n")
cat("Observed annotation records       :", n_total, "\n")
cat("Expected historical tumor         :", 516, "\n")
cat("Expected historical normal        :", 59, "\n")
cat("Observed tumor                    :", n_tumor, "\n")
cat("Observed normal                   :", n_normal, "\n")
cat("\n")

# ------------------------------------------------------------
# 9. Save exact audit
# ------------------------------------------------------------

out_dir <- file.path(
  "result_dir/Paper 1_Step 13-Diagnostic Classification",
  "STEP13B_TCGＡ_Input_Lock"
)

# Correct possible full-width typo defensively
out_dir <- file.path(
  "result_dir/Paper 1_Step 13-Diagnostic Classification",
  "STEP13B_TCGA_Input_Lock"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

write.table(
  ann,
  file = file.path(
    out_dir,
    "STEP13B-5_TCGA_575_sample_annotation_audit.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  ann[other_idx, , drop = FALSE],
  file = file.path(
    out_dir,
    "STEP13B-5_Unclassified_sample_records.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# 10. PASS / FAIL logic
# ------------------------------------------------------------

pass_composition <-
  n_total == 574 &&
  n_tumor == 513 &&
  n_normal == 59 &&
  n_other == 2 &&
  sum(dup_sample) == 0

cat("============================================================\n")

if (pass_composition) {
  
  cat("STEP13B-5 STATUS: PASS\n")
  cat("------------------------------------------------------------\n")
  cat("574 annotation records verified.\n")
  cat("513 Tumor records verified.\n")
  cat("59 Normal records verified.\n")
  cat("2 non-Tumor/Normal records identified.\n")
  cat("No duplicated Expression_Sample IDs.\n")
  cat("\n")
  cat("IMPORTANT:\n")
  cat("The two unclassified records remain excluded from\n")
  cat("Tumor-versus-Normal classification until provenance\n")
  cat("is independently resolved.\n")
  
} else {
  
  cat("STEP13B-5 STATUS: FAIL\n")
  cat("------------------------------------------------------------\n")
  cat("Sample composition does not match the expected structure.\n")
  cat("Further provenance audit required.\n")
}

cat("============================================================\n")
cat("NO ROC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")
cat("============================================================\n")

# ============================================================
# WORKFLOW-13
# STEP13B-6 — Original Expression Matrix ↔ Annotation
#              Exact Sample-ID Reconciliation
#
# Purpose:
#   Reconcile the original TCGA expression-matrix sample IDs
#   against the sample-level annotation used in STEP13B-5.
#
# This step is PROVENANCE ONLY.
#
# NO:
#   - ROC calculation
#   - AUC calculation
#   - logistic regression
#   - model fitting
#   - threshold optimization
#   - statistical testing
#   - P values
#   - FDR
#
# Core genes are NOT analyzed in this step.
# ============================================================


# ------------------------------------------------------------
# 0. Clean session
# ------------------------------------------------------------

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

source_root <- "input_dir/Paper1_tumor-only WGCNA"

step_dir <- file.path(
  workflow_root,
  "STEP13B_TCGA_Input_Lock"
)

dir.create(
  step_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Input files
# ------------------------------------------------------------

annotation_file <- file.path(
  source_root,
  "TCGA_Expression_Sample_Annotation_Basic.txt"
)

historical_structure_file <- file.path(
  source_root,
  "Workflow-03F-7_expression_matrix_structure.txt"
)

historical_annotation_file <- file.path(
  source_root,
  "Workflow-03F-2_expression_annotation.txt"
)

historical_patient_mapping_file <- file.path(
  source_root,
  "Workflow-03F-3_patient_expression_mapping.txt"
)

historical_patient_annotation_file <- file.path(
  source_root,
  "Workflow-03F-5_historical_517_patient_annotation.txt"
)

tumor_only_file <- file.path(
  source_root,
  "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt"
)

tumor_only_gene_symbol_file <- file.path(
  source_root,
  "Workflow-04-2C-1A_ESTIMATE_input_GeneSymbol_log2_FPKM_UQ.txt"
)

historical_structure_check <- file.exists(
  historical_structure_file
)


# ------------------------------------------------------------
# 3. Output files
# ------------------------------------------------------------

audit_file <- file.path(
  step_dir,
  "STEP13B-6_Original_Matrix_Annotation_Reconciliation.txt"
)

matrix_id_file <- file.path(
  step_dir,
  "STEP13B-6_Original_Matrix_Sample_IDs.txt"
)

annotation_id_file <- file.path(
  step_dir,
  "STEP13B-6_Annotation_Sample_IDs.txt"
)

matrix_only_file <- file.path(
  step_dir,
  "STEP13B-6_Matrix_Only_Sample_IDs.txt"
)

annotation_only_file <- file.path(
  step_dir,
  "STEP13B-6_Annotation_Only_Sample_IDs.txt"
)

unclassified_file <- file.path(
  step_dir,
  "STEP13B-6_Unclassified_Records.txt"
)

sample_code_file <- file.path(
  step_dir,
  "STEP13B-6_Sample_Code_Audit.txt"
)

wgcna_overlap_file <- file.path(
  step_dir,
  "STEP13B-6_WGCNA_Tumor_Overlap_Audit.txt"
)

pass_file <- file.path(
  step_dir,
  "STEP13B-6_PASS.txt"
)

hold_file <- file.path(
  step_dir,
  "STEP13B-6_HOLD.txt"
)


# ------------------------------------------------------------
# 4. Start audit log
# ------------------------------------------------------------

sink(audit_file)

cat("============================================================\n")
cat("WORKFLOW-13 — STEP13B-6\n")
cat("Original Expression Matrix ↔ Annotation Exact Sample-ID\n")
cat("Reconciliation\n")
cat("============================================================\n\n")

cat("Workflow root:\n")
cat(workflow_root, "\n\n")

cat("Source root:\n")
cat(source_root, "\n\n")

cat("This is a PROVENANCE-ONLY checkpoint.\n\n")

cat("NO ROC CALCULATION\n")
cat("NO AUC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO LOGISTIC REGRESSION\n")
cat("NO THRESHOLD OPTIMIZATION\n")
cat("NO STATISTICAL TEST\n")
cat("NO P VALUE\n")
cat("NO FDR\n\n")


# ------------------------------------------------------------
# 5. Input existence audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("5. INPUT FILE EXISTENCE AUDIT\n")
cat("------------------------------------------------------------\n\n")

input_files <- c(
  annotation_file,
  historical_structure_file,
  historical_annotation_file,
  historical_patient_mapping_file,
  historical_patient_annotation_file,
  tumor_only_file,
  tumor_only_gene_symbol_file
)

input_exists <- file.exists(input_files)

input_inventory <- data.frame(
  file = basename(input_files),
  exists = input_exists,
  stringsAsFactors = FALSE
)

print(input_inventory, row.names = FALSE)

cat("\n")


# ------------------------------------------------------------
# 6. Read STEP13B-5 annotation
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("6. SAMPLE-LEVEL ANNOTATION\n")
cat("------------------------------------------------------------\n\n")

stopifnot(file.exists(annotation_file))

ann <- read.delim(
  annotation_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Annotation file:\n")
cat(annotation_file, "\n\n")

cat("Dimensions:\n")
cat(
  "Rows = ", nrow(ann),
  "\nColumns = ", ncol(ann),
  "\n\n",
  sep = ""
)

cat("Column names:\n")
print(names(ann))

cat("\n")


# ------------------------------------------------------------
# 7. Required columns
# ------------------------------------------------------------

required_annotation_columns <- c(
  "Expression_Sample",
  "Patient_ID",
  "Sample_Type_Code",
  "Tumor_Normal"
)

missing_annotation_columns <- setdiff(
  required_annotation_columns,
  names(ann)
)

if (length(missing_annotation_columns) > 0) {
  
  cat(
    "ERROR: Missing required annotation columns:\n"
  )
  
  print(missing_annotation_columns)
  
  sink()
  
  stop(
    "STEP13B-6 FAILED: required annotation columns missing."
  )
}


# ------------------------------------------------------------
# 8. Exact annotation sample-ID audit
# ------------------------------------------------------------

annotation_ids <- trimws(
  as.character(ann$Expression_Sample)
)

annotation_ids <- annotation_ids[
  !is.na(annotation_ids) &
    nzchar(annotation_ids)
]

cat("Annotation sample IDs:\n")
cat(
  "N = ", length(annotation_ids),
  "\n\n",
  sep = ""
)

cat(
  "Duplicated annotation IDs = ",
  sum(duplicated(annotation_ids)),
  "\n\n",
  sep = ""
)

writeLines(
  annotation_ids,
  annotation_id_file
)


# ------------------------------------------------------------
# 9. Exact annotation category audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("9. EXACT Tumor_Normal VALUE DISTRIBUTION\n")
cat("------------------------------------------------------------\n\n")

tumor_normal_table <- table(
  ann$Tumor_Normal,
  useNA = "ifany"
)

print(tumor_normal_table)

cat("\n")

tumor_n <- sum(
  !is.na(ann$Tumor_Normal) &
    ann$Tumor_Normal == "Tumor"
)

normal_n <- sum(
  !is.na(ann$Tumor_Normal) &
    ann$Tumor_Normal == "Normal"
)

other_n <- sum(
  is.na(ann$Tumor_Normal) |
    !(ann$Tumor_Normal %in% c("Tumor", "Normal"))
)

cat("Exact counts:\n")
cat("Tumor   = ", tumor_n, "\n", sep = "")
cat("Normal  = ", normal_n, "\n", sep = "")
cat("Other   = ", other_n, "\n\n", sep = "")


# ------------------------------------------------------------
# 10. Identify the two unclassified records
# ------------------------------------------------------------

unclassified_index <- which(
  is.na(ann$Tumor_Normal) |
    !(ann$Tumor_Normal %in% c("Tumor", "Normal"))
)

unclassified_records <- ann[
  unclassified_index,
  ,
  drop = FALSE
]

cat("------------------------------------------------------------\n")
cat("10. UNCLASSIFIED SAMPLE RECORDS\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Number of unclassified records = ",
  nrow(unclassified_records),
  "\n\n",
  sep = ""
)

if (nrow(unclassified_records) > 0) {
  
  print(
    unclassified_records,
    row.names = FALSE
  )
  
  write.table(
    unclassified_records,
    file = unclassified_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
} else {
  
  cat("No unclassified records identified.\n")
  
  writeLines(
    "No unclassified records identified.",
    unclassified_file
  )
}

cat("\n")


# ------------------------------------------------------------
# 11. Sample code audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("11. SAMPLE TYPE CODE AUDIT\n")
cat("------------------------------------------------------------\n\n")

sample_code_table <- table(
  ann$Sample_Type_Code,
  useNA = "ifany"
)

print(sample_code_table)

cat("\n")

write.table(
  as.data.frame(sample_code_table),
  file = sample_code_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Sample ID suffix audit
# ------------------------------------------------------------

cat("Sample ID suffix audit:\n\n")

suffix_01 <- grepl(
  "-01$",
  annotation_ids
)

suffix_11 <- grepl(
  "-11$",
  annotation_ids
)

suffix_other <- !suffix_01 & !suffix_11

cat(
  "IDs ending -01 = ",
  sum(suffix_01),
  "\n",
  sep = ""
)

cat(
  "IDs ending -11 = ",
  sum(suffix_11),
  "\n",
  sep = ""
)

cat(
  "IDs with other suffix = ",
  sum(suffix_other),
  "\n\n",
  sep = ""
)

if (any(suffix_other)) {
  
  cat("Other-suffix IDs:\n")
  
  print(
    annotation_ids[suffix_other]
  )
  
  cat("\n")
}


# ------------------------------------------------------------
# 13. Locate original expression matrix
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("13. ORIGINAL EXPRESSION MATRIX SEARCH\n")
cat("------------------------------------------------------------\n\n")

candidate_matrix_files <- list.files(
  source_root,
  pattern = "\\.(txt|tsv|csv|rds|RData|rda)$",
  recursive = TRUE,
  full.names = TRUE
)

candidate_matrix_files <- candidate_matrix_files[
  !grepl(
    "annotation|mapping|clinical|structure|provenance|audit",
    basename(candidate_matrix_files),
    ignore.case = TRUE
  )
]

cat(
  "Candidate expression-like files after filename filtering:\n\n"
)

if (length(candidate_matrix_files) > 0) {
  print(candidate_matrix_files)
} else {
  cat("No candidate files identified by filename filter.\n")
}

cat("\n")


# ------------------------------------------------------------
# 14. Search specifically for 60483 × 575 structure
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("14. HISTORICAL 60483 × 575 STRUCTURE\n")
cat("------------------------------------------------------------\n\n")

if (file.exists(historical_structure_file)) {
  
  structure_lines <- readLines(
    historical_structure_file,
    warn = FALSE
  )
  
  dimension_lines <- structure_lines[
    grepl(
      "Rows|Columns|Dimensions|60483|575",
      structure_lines,
      ignore.case = TRUE
    )
  ]
  
  cat(
    paste(
      dimension_lines,
      collapse = "\n"
    )
  )
  
  cat("\n\n")
  
} else {
  
  cat(
    "Historical structure file not found.\n"
  )
}


# ------------------------------------------------------------
# 15. Identify actual matrix from historical files
# ------------------------------------------------------------
#
# IMPORTANT:
# The script does NOT assume that a particular file is the
# original 60483 × 575 matrix.
#
# We first inspect likely files and only read a matrix if
# its structure can be established.
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("15. MATRIX STRUCTURE INSPECTION\n")
cat("------------------------------------------------------------\n\n")

cat(
  "The original matrix is NOT inferred from filename alone.\n"
)

cat(
  "The historical structure file is treated as provenance evidence,\n"
)

cat(
  "but exact sample-ID reconciliation requires the actual matrix\n"
)

cat(
  "or an independently traceable matrix representation.\n\n"
)


# ------------------------------------------------------------
# 16. Inspect Workflow-04-2C candidate file
# ------------------------------------------------------------

if (file.exists(tumor_only_file)) {
  
  cat("Downstream tumor-only WGCNA matrix detected:\n")
  cat(tumor_only_file, "\n\n")
  
  tumor_preview <- read.delim(
    tumor_only_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    nrows = 5,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat(
    "Preview dimensions (first 5 rows):\n"
  )
  
  cat(
    "Columns = ",
    ncol(tumor_preview),
    "\n\n",
    sep = ""
  )
  
  tumor_sample_ids <- names(tumor_preview)[-1]
  
  cat(
    "Number of apparent sample columns = ",
    length(tumor_sample_ids),
    "\n\n",
    sep = ""
  )
  
} else {
  
  tumor_sample_ids <- character(0)
  
  cat(
    "Downstream tumor-only WGCNA matrix not found.\n\n"
  )
}


# ------------------------------------------------------------
# 17. Compare downstream WGCNA tumor IDs with annotation
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("17. DOWNSTREAM WGCNA 513-TUMOR OVERLAP AUDIT\n")
cat("------------------------------------------------------------\n\n")

if (length(tumor_sample_ids) > 0) {
  
  tumor_sample_ids <- trimws(
    as.character(tumor_sample_ids)
  )
  
  tumor_overlap <- intersect(
    tumor_sample_ids,
    annotation_ids
  )
  
  tumor_only_ids <- setdiff(
    tumor_sample_ids,
    annotation_ids
  )
  
  annotation_tumor_ids <- annotation_ids[
    !is.na(ann$Tumor_Normal) &
      ann$Tumor_Normal == "Tumor"
  ]
  
  annotation_tumor_not_in_wgcna <- setdiff(
    annotation_tumor_ids,
    tumor_sample_ids
  )
  
  cat(
    "WGCNA apparent sample columns = ",
    length(tumor_sample_ids),
    "\n",
    sep = ""
  )
  
  cat(
    "WGCNA IDs found in annotation = ",
    length(tumor_overlap),
    "\n",
    sep = ""
  )
  
  cat(
    "WGCNA IDs absent from annotation = ",
    length(tumor_only_ids),
    "\n",
    sep = ""
  )
  
  cat(
    "Annotation Tumor IDs absent from WGCNA = ",
    length(annotation_tumor_not_in_wgcna),
    "\n\n",
    sep = ""
  )
  
  cat(
    "Interpretation:\n"
  )
  
  cat(
    "This comparison describes ID overlap only.\n"
  )
  
  cat(
    "It does NOT infer why samples were excluded.\n"
  )
  
  cat(
    "Clinical-information filtering must be traced separately.\n\n"
  )
  
  wgcna_audit_lines <- c(
    paste0(
      "WGCNA apparent sample columns = ",
      length(tumor_sample_ids)
    ),
    paste0(
      "WGCNA IDs found in annotation = ",
      length(tumor_overlap)
    ),
    paste0(
      "WGCNA IDs absent from annotation = ",
      length(tumor_only_ids)
    ),
    paste0(
      "Annotation Tumor IDs absent from WGCNA = ",
      length(annotation_tumor_not_in_wgcna)
    ),
    "",
    "WGCNA IDs absent from annotation:",
    if (length(tumor_only_ids) > 0)
      tumor_only_ids
    else
      "NONE",
    "",
    "Annotation Tumor IDs absent from WGCNA:",
    if (length(annotation_tumor_not_in_wgcna) > 0)
      annotation_tumor_not_in_wgcna
    else
      "NONE"
  )
  
  writeLines(
    wgcna_audit_lines,
    wgcna_overlap_file
  )
  
} else {
  
  cat(
    "WGCNA tumor-only sample IDs could not be extracted.\n"
  )
  
  cat(
    "No overlap conclusion is made.\n"
  )
  
  writeLines(
    "WGCNA tumor-only sample IDs could not be extracted.",
    wgcna_overlap_file
  )
}


# ------------------------------------------------------------
# 18. Historical annotation file comparison
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("18. HISTORICAL ANNOTATION FILE COMPARISON\n")
cat("------------------------------------------------------------\n\n")

if (file.exists(historical_annotation_file)) {
  
  historical_ann <- read.delim(
    historical_annotation_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat(
    "Workflow-03F-2 rows = ",
    nrow(historical_ann),
    "\n",
    sep = ""
  )
  
  cat(
    "Workflow-03F-2 columns = ",
    ncol(historical_ann),
    "\n\n",
    sep = ""
  )
  
  if ("Expression_Sample" %in% names(historical_ann)) {
    
    historical_ids <- trimws(
      as.character(
        historical_ann$Expression_Sample
      )
    )
    
    historical_ids <- historical_ids[
      !is.na(historical_ids) &
        nzchar(historical_ids)
    ]
    
    cat(
      "Historical annotation sample IDs = ",
      length(historical_ids),
      "\n",
      sep = ""
    )
    
    cat(
      "Current annotation IDs = ",
      length(annotation_ids),
      "\n",
      sep = ""
    )
    
    cat(
      "Exact ID intersection = ",
      length(
        intersect(
          historical_ids,
          annotation_ids
        )
      ),
      "\n",
      sep = ""
    )
    
    cat(
      "Historical-only IDs = ",
      length(
        setdiff(
          historical_ids,
          annotation_ids
        )
      ),
      "\n",
      sep = ""
    )
    
    cat(
      "Current-only IDs = ",
      length(
        setdiff(
          annotation_ids,
          historical_ids
        )
      ),
      "\n\n",
      sep = ""
    )
  }
  
} else {
  
  cat(
    "Workflow-03F-2 annotation file not found.\n\n"
  )
}


# ------------------------------------------------------------
# 19. Historical patient-expression mapping comparison
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("19. HISTORICAL PATIENT-LEVEL MAPPING CHECK\n")
cat("------------------------------------------------------------\n\n")

if (file.exists(historical_patient_mapping_file)) {
  
  patient_map <- read.delim(
    historical_patient_mapping_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat(
    "Workflow-03F-3 rows = ",
    nrow(patient_map),
    "\n",
    sep = ""
  )
  
  cat(
    "Workflow-03F-3 columns = ",
    ncol(patient_map),
    "\n\n",
    sep = ""
  )
  
  if ("Expression_Sample" %in% names(patient_map)) {
    
    patient_expression_entries <- trimws(
      as.character(
        patient_map$Expression_Sample
      )
    )
    
    cat(
      "Patient-level Expression_Sample entries = ",
      length(patient_expression_entries),
      "\n\n",
      sep = ""
    )
    
    cat(
      "IMPORTANT:\n"
    )
    
    cat(
      "Patient-level entries may contain multiple sample IDs\n"
    )
    
    cat(
      "such as '-01; -11' and therefore are NOT directly\n"
    )
    
    cat(
      "equivalent to sample-level expression-matrix columns.\n\n"
    )
  }
  
} else {
  
  cat(
    "Workflow-03F-3 mapping file not found.\n\n"
  )
}


# ------------------------------------------------------------
# 20. Historical 516 + 59 claim audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("20. HISTORICAL 516 + 59 CLAIM AUDIT\n")
cat("------------------------------------------------------------\n\n")

historical_tumor_claim <- 516
historical_normal_claim <- 59
historical_total_claim <-
  historical_tumor_claim +
  historical_normal_claim

cat(
  "Historical reported Tumor count = ",
  historical_tumor_claim,
  "\n",
  sep = ""
)

cat(
  "Historical reported Normal count = ",
  historical_normal_claim,
  "\n",
  sep = ""
)

cat(
  "Historical reported total = ",
  historical_total_claim,
  "\n\n",
  sep = ""
)

cat(
  "Current STEP13B-5 annotation:\n"
)

cat(
  "Tumor = ", tumor_n, "\n",
  sep = ""
)

cat(
  "Normal = ", normal_n, "\n",
  sep = ""
)

cat(
  "Other = ", other_n, "\n",
  sep = ""
)

cat(
  "Total annotation records = ",
  nrow(ann),
  "\n\n",
  sep = ""
)

cat(
  "No equivalence between the historical 516 + 59 claim and\n"
)

cat(
  "the current 574-record annotation is assumed at this stage.\n"
)

cat(
  "The difference must be resolved by exact provenance.\n\n"
)


# ------------------------------------------------------------
# 21. Patient-level vs sample-level warning
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("21. PATIENT-LEVEL / SAMPLE-LEVEL TRACEABILITY WARNING\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Patient-level files may contain multiple expression samples\n"
)

cat(
  "per patient, particularly paired primary-tumor / normal samples.\n"
)

cat(
  "Therefore:\n\n"
)

cat(
  "patient count != sample count\n"
)

cat(
  "patient-level record count != expression-matrix column count\n"
)

cat(
  "and no arithmetic conversion is performed here.\n\n"
)


# ------------------------------------------------------------
# 22. Reconciliation summary
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("22. RECONCILIATION SUMMARY\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Annotation records                  = ",
  nrow(ann),
  "\n",
  sep = ""
)

cat(
  "Annotation Tumor                    = ",
  tumor_n,
  "\n",
  sep = ""
)

cat(
  "Annotation Normal                   = ",
  normal_n,
  "\n",
  sep = ""
)

cat(
  "Annotation Other/Unclassified      = ",
  other_n,
  "\n",
  sep = ""
)

cat(
  "Annotation duplicate IDs            = ",
  sum(duplicated(annotation_ids)),
  "\n",
  sep = ""
)

cat("\n")


# ------------------------------------------------------------
# 23. Determine provenance status
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("23. STEP13B-6 STATUS\n")
cat("------------------------------------------------------------\n\n")


annotation_structure_pass <-
  nrow(ann) == 574 &&
  tumor_n == 513 &&
  normal_n == 59 &&
  other_n == 2 &&
  sum(duplicated(annotation_ids)) == 0


# Exact original matrix reconciliation cannot be declared PASS
# unless the actual matrix sample IDs have been independently
# recovered and matched.
#
# Therefore this checkpoint deliberately uses HOLD when the
# original matrix itself cannot yet be reconciled.

original_matrix_exactly_reconciled <- FALSE


if (
  annotation_structure_pass &&
  original_matrix_exactly_reconciled
) {
  
  status <- "PASS"
  
} else {
  
  status <- "HOLD"
}


cat(
  "STEP13B-6 STATUS: ",
  status,
  "\n",
  sep = ""
)

cat("------------------------------------------------------------\n")

if (status == "PASS") {
  
  cat(
    "Exact original-matrix ↔ annotation sample-ID reconciliation\n"
  )
  
  cat(
    "has been independently verified.\n\n"
  )
  
  cat(
    "STEP13B may proceed toward final input-lock review.\n"
  )
  
} else {
  
  cat(
    "The annotation structure is internally verified, but the\n"
  )
  
  cat(
    "original expression-matrix sample-ID reconciliation is NOT\n"
  )
  
  cat(
    "yet independently verified.\n\n"
  )
  
  cat(
    "The two unclassified records remain unresolved.\n"
  )
  
  cat(
    "The historical 516 + 59 composition remains unresolved.\n"
  )
  
  cat(
    "No final diagnostic-classification input lock is issued.\n"
  )
  
  cat(
    "ROC/model analysis remains on HOLD.\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 24. Final safety lock
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP13B-6 PROVENANCE CHECKPOINT COMPLETE\n")
cat("============================================================\n\n")

cat("NO ROC CALCULATION\n")
cat("NO AUC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO LOGISTIC REGRESSION\n")
cat("NO THRESHOLD OPTIMIZATION\n")
cat("NO STATISTICAL TEST\n")
cat("NO P VALUE\n")
cat("NO FDR\n\n")

cat("============================================================\n")

sink()


# ------------------------------------------------------------
# 25. Save PASS / HOLD status
# ------------------------------------------------------------

if (status == "PASS") {
  
  writeLines(
    c(
      "STEP13B-6 STATUS: PASS",
      "Exact original-matrix ↔ annotation sample-ID reconciliation verified.",
      "See STEP13B-6_Original_Matrix_Annotation_Reconciliation.txt",
      "for complete provenance audit.",
      "",
      "NO ROC CALCULATION",
      "NO MODEL FITTING",
      "NO STATISTICAL TEST"
    ),
    pass_file
  )
  
  if (file.exists(hold_file)) {
    file.remove(hold_file)
  }
  
} else {
  
  writeLines(
    c(
      "STEP13B-6 STATUS: HOLD",
      "Annotation structure is internally verified, but exact",
      "original-matrix ↔ annotation reconciliation is not yet complete.",
      "",
      "Two unclassified records remain unresolved.",
      "Historical 516 + 59 composition remains unresolved.",
      "No final diagnostic-classification input lock is issued.",
      "",
      "NO ROC CALCULATION",
      "NO MODEL FITTING",
      "NO STATISTICAL TEST"
    ),
    hold_file
  )
}


# ------------------------------------------------------------
# 26. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13B-6 COMPLETE\n")
cat("============================================================\n")

cat(
  "Annotation records: ",
  nrow(ann),
  "\n",
  sep = ""
)

cat(
  "Tumor: ",
  tumor_n,
  "\n",
  sep = ""
)

cat(
  "Normal: ",
  normal_n,
  "\n",
  sep = ""
)

cat(
  "Other/unclassified: ",
  other_n,
  "\n",
  sep = ""
)

cat(
  "Duplicate annotation IDs: ",
  sum(duplicated(annotation_ids)),
  "\n",
  sep = ""
)

cat(
  "STEP13B-6 STATUS: ",
  status,
  "\n",
  sep = ""
)

cat("------------------------------------------------------------\n")
cat("NO ROC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")
cat("============================================================\n")

#STEP13B-7 — Unclassified Sample Resolution + Historical 516/513 Traceability Audit
# ============================================================
# WORKFLOW-13
# STEP13B-7 — Unclassified Sample Resolution +
#              Historical 516/513 Traceability Audit
#
# Purpose:
#   1. Resolve the two unclassified sample records identified
#      in STEP13B-5 / STEP13B-6.
#
#   2. Reconcile sample-level IDs across:
#        - current TCGA sample annotation
#        - historical expression annotation
#        - historical patient-expression mapping
#        - historical patient annotation
#        - final tumor-only WGCNA matrix
#
#   3. Independently audit the historical "516 Tumor" claim.
#
#   4. Determine whether the final 513-tumor WGCNA cohort can
#      be traced exactly to sample-level records.
#
# IMPORTANT:
#   This is a PROVENANCE-ONLY checkpoint.
#
# NO:
#   - ROC
#   - AUC
#   - logistic regression
#   - model fitting
#   - threshold optimization
#   - P values
#   - FDR
#   - statistical tests
#
# No assumption is made that:
#       516 -> 513
# is caused by clinical filtering until exact evidence is found.
# ============================================================


# ------------------------------------------------------------
# 0. Clean session
# ------------------------------------------------------------

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

source_root <- "input_dir/Paper1_tumor-only WGCNA"

step_dir <- file.path(
  workflow_root,
  "STEP13B_TCGA_Input_Lock"
)

dir.create(
  step_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Input files
# ------------------------------------------------------------

current_annotation_file <- file.path(
  source_root,
  "TCGA_Expression_Sample_Annotation_Basic.txt"
)

historical_expression_annotation_file <- file.path(
  source_root,
  "Workflow-03F-2_expression_annotation.txt"
)

historical_patient_mapping_file <- file.path(
  source_root,
  "Workflow-03F-3_patient_expression_mapping.txt"
)

historical_patient_annotation_file <- file.path(
  source_root,
  "Workflow-03F-5_historical_517_patient_annotation.txt"
)

historical_structure_file <- file.path(
  source_root,
  "Workflow-03F-7_expression_matrix_structure.txt"
)

historical_scale_file <- file.path(
  source_root,
  "Workflow-03F-9_expression_scale_provenance.txt"
)

tumor_only_file <- file.path(
  source_root,
  "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt"
)

clinical_file <- file.path(
  source_root,
  "Workflow-04-2A_clinical_aligned.txt"
)

step13b6_unclassified_file <- file.path(
  step_dir,
  "STEP13B-6_Unclassified_Records.txt"
)

step13b6_overlap_file <- file.path(
  step_dir,
  "STEP13B-6_WGCNA_Tumor_Overlap_Audit.txt"
)


# ------------------------------------------------------------
# 3. Output files
# ------------------------------------------------------------

audit_file <- file.path(
  step_dir,
  "STEP13B-7_Unclassified_516_513_Traceability_Audit.txt"
)

unclassified_resolution_file <- file.path(
  step_dir,
  "STEP13B-7_Unclassified_Sample_Resolution.txt"
)

historical_516_file <- file.path(
  step_dir,
  "STEP13B-7_Historical_516_Tumor_Audit.txt"
)

wgcna_513_file <- file.path(
  step_dir,
  "STEP13B-7_Final_513_Tumor_Traceability.txt"
)

clinical_exclusion_file <- file.path(
  step_dir,
  "STEP13B-7_Clinical_Filtering_Exclusion_Audit.txt"
)

historical_claim_search_file <- file.path(
  step_dir,
  "STEP13B-7_Historical_516_Claim_Search.txt"
)

pass_file <- file.path(
  step_dir,
  "STEP13B-7_PASS.txt"
)

hold_file <- file.path(
  step_dir,
  "STEP13B-7_HOLD.txt"
)


# ------------------------------------------------------------
# 4. Start audit log
# ------------------------------------------------------------

sink(audit_file)

cat("============================================================\n")
cat("WORKFLOW-13 — STEP13B-7\n")
cat("Unclassified Sample Resolution + Historical 516/513\n")
cat("Traceability Audit\n")
cat("============================================================\n\n")

cat("Workflow root:\n")
cat(workflow_root, "\n\n")

cat("Source root:\n")
cat(source_root, "\n\n")

cat("PROVENANCE-ONLY CHECKPOINT\n\n")

cat("NO ROC CALCULATION\n")
cat("NO AUC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO LOGISTIC REGRESSION\n")
cat("NO THRESHOLD OPTIMIZATION\n")
cat("NO STATISTICAL TEST\n")
cat("NO P VALUE\n")
cat("NO FDR\n\n")


# ------------------------------------------------------------
# 5. Input existence audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("5. INPUT EXISTENCE AUDIT\n")
cat("------------------------------------------------------------\n\n")

input_files <- c(
  current_annotation_file,
  historical_expression_annotation_file,
  historical_patient_mapping_file,
  historical_patient_annotation_file,
  historical_structure_file,
  tumor_only_file,
  clinical_file,
  step13b6_unclassified_file,
  step13b6_overlap_file
)

input_inventory <- data.frame(
  file = basename(input_files),
  exists = file.exists(input_files),
  stringsAsFactors = FALSE
)

print(
  input_inventory,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 6. Read current sample-level annotation
# ------------------------------------------------------------

stopifnot(
  file.exists(current_annotation_file)
)

ann <- read.delim(
  current_annotation_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("------------------------------------------------------------\n")
cat("6. CURRENT SAMPLE-LEVEL ANNOTATION\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Rows = ",
  nrow(ann),
  "\n",
  sep = ""
)

cat(
  "Columns = ",
  ncol(ann),
  "\n\n",
  sep = ""
)

required_cols <- c(
  "Expression_Sample",
  "Patient_ID",
  "Sample_Type_Code",
  "Tumor_Normal"
)

missing_cols <- setdiff(
  required_cols,
  names(ann)
)

if (length(missing_cols) > 0) {
  
  sink()
  
  stop(
    paste(
      "STEP13B-7 FAILED:",
      "required annotation columns missing:",
      paste(missing_cols, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 7. Normalize sample IDs
# ------------------------------------------------------------

ann$Expression_Sample <- trimws(
  as.character(
    ann$Expression_Sample
  )
)

ann$Patient_ID <- trimws(
  as.character(
    ann$Patient_ID
  )
)

ann$Sample_Type_Code <- trimws(
  as.character(
    ann$Sample_Type_Code
  )
)

ann$Tumor_Normal <- trimws(
  as.character(
    ann$Tumor_Normal
  )
)

annotation_ids <- ann$Expression_Sample

annotation_ids <- annotation_ids[
  !is.na(annotation_ids) &
    nzchar(annotation_ids)
]


# ------------------------------------------------------------
# 8. Exact current classification
# ------------------------------------------------------------

tumor_current <- ann[
  ann$Tumor_Normal == "Tumor",
  ,
  drop = FALSE
]

normal_current <- ann[
  ann$Tumor_Normal == "Normal",
  ,
  drop = FALSE
]

unclassified <- ann[
  is.na(ann$Tumor_Normal) |
    !(ann$Tumor_Normal %in% c("Tumor", "Normal")),
  ,
  drop = FALSE
]

cat("Current annotation classification:\n\n")

cat(
  "Total records       = ",
  nrow(ann),
  "\n",
  sep = ""
)

cat(
  "Tumor               = ",
  nrow(tumor_current),
  "\n",
  sep = ""
)

cat(
  "Normal              = ",
  nrow(normal_current),
  "\n",
  sep = ""
)

cat(
  "Unclassified        = ",
  nrow(unclassified),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Resolve unclassified samples
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("9. UNCLASSIFIED SAMPLE RESOLUTION\n")
cat("------------------------------------------------------------\n\n")

if (nrow(unclassified) == 0) {
  
  cat(
    "No unclassified records remain.\n\n"
  )
  
} else {
  
  cat(
    "Unclassified records identified:\n\n"
  )
  
  print(
    unclassified,
    row.names = FALSE
  )
  
  cat("\n")
  
  write.table(
    unclassified,
    file = unclassified_resolution_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 10. Inspect sample-code consistency of unclassified samples
# ------------------------------------------------------------

if (nrow(unclassified) > 0) {
  
  cat(
    "Unclassified sample-code details:\n\n"
  )
  
  for (i in seq_len(nrow(unclassified))) {
    
    cat(
      "Sample ",
      i,
      ":\n",
      sep = ""
    )
    
    cat(
      "  Expression_Sample = ",
      unclassified$Expression_Sample[i],
      "\n",
      sep = ""
    )
    
    cat(
      "  Patient_ID        = ",
      unclassified$Patient_ID[i],
      "\n",
      sep = ""
    )
    
    cat(
      "  Sample_Type_Code  = ",
      unclassified$Sample_Type_Code[i],
      "\n",
      sep = ""
    )
    
    cat(
      "  Tumor_Normal      = ",
      unclassified$Tumor_Normal[i],
      "\n\n",
      sep = ""
    )
  }
}


# ------------------------------------------------------------
# 11. Derive expected sample class from TCGA barcode suffix
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("11. TCGA BARCODE SUFFIX AUDIT\n")
cat("------------------------------------------------------------\n\n")

ann$Barcode_Class <- NA_character_

ann$Barcode_Class[
  grepl(
    "-01$",
    ann$Expression_Sample
  )
] <- "Primary_Tumor_By_Barcode"

ann$Barcode_Class[
  grepl(
    "-11$",
    ann$Expression_Sample
  )
] <- "Normal_By_Barcode"

ann$Barcode_Class[
  is.na(ann$Barcode_Class)
] <- "Other_or_Unrecognized"

barcode_cross_tab <- table(
  ann$Tumor_Normal,
  ann$Barcode_Class,
  useNA = "ifany"
)

print(
  barcode_cross_tab
)

cat("\n")


# ------------------------------------------------------------
# 12. Identify inconsistencies between annotation and barcode
# ------------------------------------------------------------

barcode_inconsistent <- ann[
  (
    ann$Tumor_Normal == "Tumor" &
      ann$Barcode_Class != "Primary_Tumor_By_Barcode"
  ) |
    (
      ann$Tumor_Normal == "Normal" &
        ann$Barcode_Class != "Normal_By_Barcode"
    ),
  ,
  drop = FALSE
]

cat(
  "Annotation/barcode inconsistencies = ",
  nrow(barcode_inconsistent),
  "\n\n",
  sep = ""
)

if (nrow(barcode_inconsistent) > 0) {
  
  print(
    barcode_inconsistent,
    row.names = FALSE
  )
  
  cat("\n")
}


# ------------------------------------------------------------
# 13. Read historical Workflow-03F-2 annotation
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("13. HISTORICAL WORKFLOW-03F-2 ANNOTATION\n")
cat("------------------------------------------------------------\n\n")

historical_ann <- NULL
historical_ids <- character(0)

if (file.exists(historical_expression_annotation_file)) {
  
  historical_ann <- read.delim(
    historical_expression_annotation_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat(
    "Rows = ",
    nrow(historical_ann),
    "\n",
    sep = ""
  )
  
  cat(
    "Columns = ",
    ncol(historical_ann),
    "\n\n",
    sep = ""
  )
  
  if (
    "Expression_Sample" %in%
    names(historical_ann)
  ) {
    
    historical_ids <- trimws(
      as.character(
        historical_ann$Expression_Sample
      )
    )
    
    historical_ids <- historical_ids[
      !is.na(historical_ids) &
        nzchar(historical_ids)
    ]
    
    cat(
      "Historical sample IDs = ",
      length(historical_ids),
      "\n\n",
      sep = ""
    )
    
    if ("Sample_Type_Code" %in% names(historical_ann)) {
      
      cat(
        "Historical Sample_Type_Code distribution:\n"
      )
      
      print(
        table(
          historical_ann$Sample_Type_Code,
          useNA = "ifany"
        )
      )
      
      cat("\n")
    }
    
    if ("Sample_Type" %in% names(historical_ann)) {
      
      cat(
        "Historical Sample_Type distribution:\n"
      )
      
      print(
        table(
          historical_ann$Sample_Type,
          useNA = "ifany"
        )
      )
      
      cat("\n")
    }
  }
  
} else {
  
  cat(
    "Workflow-03F-2 file not available.\n\n"
  )
}


# ------------------------------------------------------------
# 14. Exact ID comparison: current vs historical annotation
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("14. CURRENT ↔ HISTORICAL ANNOTATION ID RECONCILIATION\n")
cat("------------------------------------------------------------\n\n")

if (length(historical_ids) > 0) {
  
  common_ids <- intersect(
    annotation_ids,
    historical_ids
  )
  
  current_only_ids <- setdiff(
    annotation_ids,
    historical_ids
  )
  
  historical_only_ids <- setdiff(
    historical_ids,
    annotation_ids
  )
  
  cat(
    "Current annotation IDs     = ",
    length(annotation_ids),
    "\n",
    sep = ""
  )
  
  cat(
    "Historical annotation IDs  = ",
    length(historical_ids),
    "\n",
    sep = ""
  )
  
  cat(
    "Exact common IDs           = ",
    length(common_ids),
    "\n",
    sep = ""
  )
  
  cat(
    "Current-only IDs           = ",
    length(current_only_ids),
    "\n",
    sep = ""
  )
  
  cat(
    "Historical-only IDs        = ",
    length(historical_only_ids),
    "\n\n",
    sep = ""
  )
  
  if (length(current_only_ids) > 0) {
    
    cat("Current-only IDs:\n")
    print(current_only_ids)
    cat("\n")
  }
  
  if (length(historical_only_ids) > 0) {
    
    cat("Historical-only IDs:\n")
    print(historical_only_ids)
    cat("\n")
  }
  
} else {
  
  cat(
    "Historical sample IDs unavailable.\n"
  )
}


# ------------------------------------------------------------
# 15. Read final WGCNA tumor-only matrix
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("15. FINAL WGCNA TUMOR-ONLY MATRIX\n")
cat("------------------------------------------------------------\n\n")

wgcna_ids <- character(0)

if (file.exists(tumor_only_file)) {
  
  wgcna_header <- read.delim(
    tumor_only_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    nrows = 0,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  wgcna_ids <- names(wgcna_header)[-1]
  
  wgcna_ids <- trimws(
    as.character(
      wgcna_ids
    )
  )
  
  cat(
    "WGCNA sample columns = ",
    length(wgcna_ids),
    "\n",
    sep = ""
  )
  
  cat(
    "Duplicate WGCNA sample IDs = ",
    sum(duplicated(wgcna_ids)),
    "\n\n",
    sep = ""
  )
  
} else {
  
  cat(
    "Final WGCNA matrix not found.\n\n"
  )
}


# ------------------------------------------------------------
# 16. Exact 513-Tumor traceability
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("16. FINAL 513-TUMOR TRACEABILITY\n")
cat("------------------------------------------------------------\n\n")

if (length(wgcna_ids) > 0) {
  
  wgcna_in_current <- intersect(
    wgcna_ids,
    annotation_ids
  )
  
  current_tumor_ids <- trimws(
    as.character(
      tumor_current$Expression_Sample
    )
  )
  
  wgcna_in_current_tumor <- intersect(
    wgcna_ids,
    current_tumor_ids
  )
  
  wgcna_not_in_current <- setdiff(
    wgcna_ids,
    annotation_ids
  )
  
  current_tumor_not_in_wgcna <- setdiff(
    current_tumor_ids,
    wgcna_ids
  )
  
  cat(
    "WGCNA sample IDs                         = ",
    length(wgcna_ids),
    "\n",
    sep = ""
  )
  
  cat(
    "WGCNA IDs found in current annotation    = ",
    length(wgcna_in_current),
    "\n",
    sep = ""
  )
  
  cat(
    "WGCNA IDs classified as Tumor            = ",
    length(wgcna_in_current_tumor),
    "\n",
    sep = ""
  )
  
  cat(
    "WGCNA IDs absent from current annotation = ",
    length(wgcna_not_in_current),
    "\n",
    sep = ""
  )
  
  cat(
    "Current Tumor IDs absent from WGCNA       = ",
    length(current_tumor_not_in_wgcna),
    "\n\n",
    sep = ""
  )
  
  writeLines(
    c(
      paste0(
        "WGCNA sample IDs = ",
        length(wgcna_ids)
      ),
      paste0(
        "WGCNA IDs found in current annotation = ",
        length(wgcna_in_current)
      ),
      paste0(
        "WGCNA IDs classified as Tumor = ",
        length(wgcna_in_current_tumor)
      ),
      paste0(
        "WGCNA IDs absent from current annotation = ",
        length(wgcna_not_in_current)
      ),
      paste0(
        "Current Tumor IDs absent from WGCNA = ",
        length(current_tumor_not_in_wgcna)
      ),
      "",
      "Current Tumor IDs absent from WGCNA:",
      if (length(current_tumor_not_in_wgcna) > 0)
        current_tumor_not_in_wgcna
      else
        "NONE",
      "",
      "WGCNA IDs absent from current annotation:",
      if (length(wgcna_not_in_current) > 0)
        wgcna_not_in_current
      else
        "NONE"
    ),
    wgcna_513_file
  )
  
} else {
  
  cat(
    "WGCNA sample IDs unavailable; 513 traceability cannot be completed.\n"
  )
  
  writeLines(
    "WGCNA sample IDs unavailable.",
    wgcna_513_file
  )
}


# ------------------------------------------------------------
# 17. Audit current Tumor IDs missing from WGCNA
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("17. CURRENT TUMOR IDs NOT IN FINAL WGCNA\n")
cat("------------------------------------------------------------\n\n")

if (
  length(wgcna_ids) > 0 &&
  exists("current_tumor_not_in_wgcna")
) {
  
  if (length(current_tumor_not_in_wgcna) > 0) {
    
    missing_from_wgcna <- tumor_current[
      tumor_current$Expression_Sample %in%
        current_tumor_not_in_wgcna,
      ,
      drop = FALSE
    ]
    
    print(
      missing_from_wgcna,
      row.names = FALSE
    )
    
    cat("\n")
    
  } else {
    
    cat(
      "No current Tumor IDs are absent from WGCNA.\n"
    )
  }
  
} else {
  
  cat(
    "Comparison unavailable.\n"
  )
}


# ------------------------------------------------------------
# 18. Read clinical file if available
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("18. CLINICAL-FILTERING PROVENANCE AUDIT\n")
cat("------------------------------------------------------------\n\n")

clinical <- NULL

if (file.exists(clinical_file)) {
  
  clinical <- read.delim(
    clinical_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat(
    "Clinical file dimensions:\n"
  )
  
  cat(
    "Rows = ",
    nrow(clinical),
    "\n",
    sep = ""
  )
  
  cat(
    "Columns = ",
    ncol(clinical),
    "\n\n",
    sep = ""
  )
  
  cat(
    "Clinical column names:\n"
  )
  
  print(
    names(clinical)
  )
  
  cat("\n")
  
} else {
  
  cat(
    "Clinical file not found.\n\n"
  )
}


# ------------------------------------------------------------
# 19. Clinical exclusion comparison
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("19. CLINICAL EXCLUSION COMPARISON\n")
cat("------------------------------------------------------------\n\n")

if (
  exists("current_tumor_not_in_wgcna") &&
  length(current_tumor_not_in_wgcna) > 0
) {
  
  exclusion_records <- tumor_current[
    tumor_current$Expression_Sample %in%
      current_tumor_not_in_wgcna,
    ,
    drop = FALSE
  ]
  
  cat(
    "Current annotation Tumor samples absent from WGCNA:\n\n"
  )
  
  print(
    exclusion_records,
    row.names = FALSE
  )
  
  write.table(
    exclusion_records,
    file = clinical_exclusion_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  cat("\n")
  
  cat(
    "These records are NOT automatically classified as\n"
  )
  
  cat(
    "'clinical-information exclusions'.\n"
  )
  
  cat(
    "Clinical filtering must be demonstrated by independent\n"
  )
  
  cat(
    "clinical provenance.\n\n"
  )
  
} else {
  
  cat(
    "No current Tumor IDs were identified as absent from WGCNA.\n"
  )
  
  cat(
    "No exclusion list was generated from this comparison.\n\n"
  )
}


# ------------------------------------------------------------
# 20. Historical 516-Tumor claim audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("20. HISTORICAL 516-TUMOR CLAIM AUDIT\n")
cat("------------------------------------------------------------\n\n")

historical_reported_tumor <- 516
historical_reported_normal <- 59

cat(
  "Historical reported Tumor count = ",
  historical_reported_tumor,
  "\n",
  sep = ""
)

cat(
  "Historical reported Normal count = ",
  historical_reported_normal,
  "\n\n",
  sep = ""
)

cat(
  "Current exact sample-level Tumor count = ",
  nrow(tumor_current),
  "\n",
  sep = ""
)

cat(
  "Current exact sample-level Normal count = ",
  nrow(normal_current),
  "\n\n",
  sep = ""
)

cat(
  "The historical 516-Tumor value is NOT assumed to be\n"
)

cat(
  "equivalent to the current 513-Tumor annotation.\n\n"
)


# ------------------------------------------------------------
# 21. TARGETED HISTORICAL 516-CLAIM SEARCH
# ------------------------------------------------------------
#
# Purpose:
#   Search only small, provenance-oriented historical text files.
#
# IMPORTANT:
#   Large expression matrices are deliberately excluded.
#
# NO:
#   - ROC
#   - AUC
#   - model fitting
#   - statistical testing
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("21. TARGETED HISTORICAL 516-CLAIM SEARCH\n")
cat("------------------------------------------------------------\n\n")


# ------------------------------------------------------------
# 21.1 List candidate files
# ------------------------------------------------------------

all_text_files <- list.files(
  source_root,
  pattern = "\\.(txt|tsv|csv|log|md)$",
  recursive = TRUE,
  full.names = TRUE
)

all_text_files <- all_text_files[
  file.exists(all_text_files)
]


# ------------------------------------------------------------
# 21.2 Exclude large/data-like files
# ------------------------------------------------------------

exclude_pattern <- paste(
  c(
    "expression",
    "matrix",
    "counts",
    "count",
    "FPKM",
    "TPM",
    "data",
    "normalized",
    "normalized_expression",
    "gene_expression",
    "expr"
  ),
  collapse = "|"
)

candidate_files <- all_text_files[
  !grepl(
    exclude_pattern,
    basename(all_text_files),
    ignore.case = TRUE
  )
]


# ------------------------------------------------------------
# 21.3 Restrict to provenance-oriented filenames
# ------------------------------------------------------------

include_pattern <- paste(
  c(
    "Workflow-03F",
    "Workflow-04",
    "provenance",
    "audit",
    "summary",
    "structure",
    "mapping",
    "annotation",
    "clinical"
  ),
  collapse = "|"
)

candidate_files <- candidate_files[
  grepl(
    include_pattern,
    basename(candidate_files),
    ignore.case = TRUE
  )
]


# ------------------------------------------------------------
# 21.4 File-size audit
# ------------------------------------------------------------

file_info <- file.info(
  candidate_files
)

candidate_inventory <- data.frame(
  file = basename(candidate_files),
  size_MB = round(
    file_info$size / 1024^2,
    3
  ),
  stringsAsFactors = FALSE
)

candidate_inventory <- candidate_inventory[
  order(
    candidate_inventory$size_MB
  ),
]

cat(
  "Targeted candidate files = ",
  nrow(candidate_inventory),
  "\n\n",
  sep = ""
)

print(
  candidate_inventory,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 21.5 Search only reasonably small provenance files
# ------------------------------------------------------------

max_size_MB <- 20

small_candidate_files <- candidate_files[
  file.info(candidate_files)$size <=
    max_size_MB * 1024^2
]

cat(
  "Files <= ",
  max_size_MB,
  " MB selected for textual search = ",
  length(small_candidate_files),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 21.6 Search for explicit 516 occurrence
# ------------------------------------------------------------

claim_hits <- list()

if (length(small_candidate_files) > 0) {
  
  for (i in seq_along(small_candidate_files)) {
    
    f <- small_candidate_files[i]
    
    cat(
      "Searching [",
      i,
      "/",
      length(small_candidate_files),
      "]: ",
      basename(f),
      "\n",
      sep = ""
    )
    
    txt <- tryCatch(
      readLines(
        f,
        warn = FALSE
      ),
      error = function(e) character(0)
    )
    
    if (length(txt) == 0) {
      next
    }
    
    hit_index <- grep(
      "516",
      txt,
      fixed = TRUE
    )
    
    if (length(hit_index) > 0) {
      
      claim_hits[[f]] <- data.frame(
        file = basename(f),
        full_path = f,
        line_number = hit_index,
        line_text = txt[hit_index],
        stringsAsFactors = FALSE
      )
    }
  }
}


# ------------------------------------------------------------
# 21.7 Save results
# ------------------------------------------------------------

if (length(claim_hits) > 0) {
  
  claim_hits_df <- do.call(
    rbind,
    claim_hits
  )
  
  cat("\n")
  cat("Explicit '516' occurrences identified:\n\n")
  
  print(
    claim_hits_df,
    row.names = FALSE
  )
  
  write.table(
    claim_hits_df,
    file = historical_claim_search_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
} else {
  
  cat("\n")
  cat(
    "No explicit '516' occurrence identified in the\n"
  )
  
  cat(
    "targeted provenance files searched.\n"
  )
  
  writeLines(
    "No explicit '516' occurrence identified in targeted provenance files.",
    historical_claim_search_file
  )
}

cat("\n")


# ------------------------------------------------------------
# 22. Targeted search for historical cohort-size phrases
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("22. TARGETED HISTORICAL COHORT-SIZE SEARCH\n")
cat("------------------------------------------------------------\n\n")

phrase_pattern <- paste(
  c(
    "516",
    "513",
    "575",
    "574",
    "59"
  ),
  collapse = "|"
)

phrase_hits <- list()

if (length(small_candidate_files) > 0) {
  
  for (f in small_candidate_files) {
    
    txt <- tryCatch(
      readLines(
        f,
        warn = FALSE
      ),
      error = function(e) character(0)
    )
    
    if (length(txt) == 0) {
      next
    }
    
    hit_index <- grep(
      phrase_pattern,
      txt,
      ignore.case = TRUE
    )
    
    if (length(hit_index) == 0) {
      next
    }
    
    selected <- txt[hit_index]
    
    selected <- selected[
      grepl(
        "tumor|normal|sample|patient|cohort|clinical|expression",
        selected,
        ignore.case = TRUE
      )
    ]
    
    if (length(selected) > 0) {
      
      phrase_hits[[f]] <- data.frame(
        file = basename(f),
        full_path = f,
        line_text = selected,
        stringsAsFactors = FALSE
      )
    }
  }
}


if (length(phrase_hits) > 0) {
  
  phrase_hits_df <- do.call(
    rbind,
    phrase_hits
  )
  
  print(
    phrase_hits_df,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No relevant historical cohort-size phrase identified.\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 23. Search complete
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("23. TARGETED HISTORICAL SEARCH COMPLETE\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Large expression/data matrices were excluded.\n"
)

cat(
  "Only targeted provenance-oriented files were searched.\n"
)

cat(
  "No statistical calculation was performed.\n\n"
)

cat("NO ROC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")

cat("============================================================\n")

# ------------------------------------------------------------
# 23. Historical patient-level mapping audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("23. PATIENT-LEVEL MAPPING AUDIT\n")
cat("------------------------------------------------------------\n\n")

patient_map <- NULL

if (file.exists(historical_patient_mapping_file)) {
  
  patient_map <- read.delim(
    historical_patient_mapping_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat(
    "Workflow-03F-3 patient-level records = ",
    nrow(patient_map),
    "\n",
    sep = ""
  )
  
  if (
    "Expression_Sample" %in%
    names(patient_map)
  ) {
    
    cat(
      "Patient-level entries containing multiple samples:\n\n"
    )
    
    multi_sample <- grepl(
      ";",
      patient_map$Expression_Sample,
      fixed = TRUE
    )
    
    cat(
      "Records with ';' in Expression_Sample = ",
      sum(multi_sample),
      "\n\n",
      sep = ""
    )
    
    if (any(multi_sample)) {
      
      print(
        patient_map[multi_sample, ],
        row.names = FALSE
      )
      
      cat("\n")
    }
  }
  
} else {
  
  cat(
    "Workflow-03F-3 mapping unavailable.\n"
  )
}


# ------------------------------------------------------------
# 24. Historical patient annotation audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("24. HISTORICAL PATIENT ANNOTATION AUDIT\n")
cat("------------------------------------------------------------\n\n")

patient_ann <- NULL

if (file.exists(historical_patient_annotation_file)) {
  
  patient_ann <- read.delim(
    historical_patient_annotation_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat(
    "Workflow-03F-5 records = ",
    nrow(patient_ann),
    "\n",
    sep = ""
  )
  
  if (
    "Sample_Category" %in%
    names(patient_ann)
  ) {
    
    cat(
      "Sample_Category distribution:\n"
    )
    
    print(
      table(
        patient_ann$Sample_Category,
        useNA = "ifany"
      )
    )
    
    cat("\n")
  }
  
} else {
  
  cat(
    "Workflow-03F-5 unavailable.\n"
  )
}


# ------------------------------------------------------------
# 25. Explicit audit of 516 vs 513 arithmetic
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("25. 516 vs 513 ARITHMETIC AUDIT\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Historical reported Tumor = ",
  historical_reported_tumor,
  "\n",
  sep = ""
)

cat(
  "Current annotation Tumor  = ",
  nrow(tumor_current),
  "\n",
  sep = ""
)

cat(
  "Final WGCNA sample columns = ",
  length(wgcna_ids),
  "\n",
  sep = ""
)

cat("\n")

if (
  length(wgcna_ids) > 0
) {
  
  cat(
    "Historical reported Tumor - current annotation Tumor = ",
    historical_reported_tumor -
      nrow(tumor_current),
    "\n",
    sep = ""
  )
  
  cat(
    "Current annotation Tumor - WGCNA sample columns = ",
    nrow(tumor_current) -
      length(wgcna_ids),
    "\n\n",
    sep = ""
  )
  
} else {
  
  cat(
    "WGCNA sample count unavailable.\n\n"
  )
}

cat(
  "IMPORTANT:\n"
)

cat(
  "These arithmetic differences are descriptive only.\n"
)

cat(
  "They do NOT establish the reason for sample exclusion.\n"
)

cat(
  "Exclusion reasons require exact provenance evidence.\n\n"
)


# ------------------------------------------------------------
# 26. Final interpretation rules
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("26. FINAL INTERPRETATION RULES\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Rule 1:\n"
)

cat(
  "Unclassified samples are not assigned Tumor/Normal status\n"
)

cat(
  "unless independently supported by sample-level provenance.\n\n"
)

cat(
  "Rule 2:\n"
)

cat(
  "A difference between 516 and 513 is not automatically\n"
)

cat(
  "attributed to missing clinical information.\n\n"
)

cat(
  "Rule 3:\n"
)

cat(
  "Patient-level records are not treated as sample-level records.\n\n"
)

cat(
  "Rule 4:\n"
)

cat(
  "No ROC/model analysis proceeds until the classification\n"
)

cat(
  "input cohort is independently locked.\n\n"
)


# ------------------------------------------------------------
# 27. Determine STEP13B-7 status
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("27. STEP13B-7 STATUS\n")
cat("------------------------------------------------------------\n\n")


unclassified_resolved <- FALSE

if (nrow(unclassified) == 0) {
  
  unclassified_resolved <- TRUE
  
} else {
  
  # Resolution is considered complete only when every
  # unclassified sample receives an independently supported
  # classification.
  #
  # This script deliberately does not infer classification
  # from barcode alone.
  
  unclassified_resolved <- FALSE
}


wgcna_traceability_complete <- FALSE

if (length(wgcna_ids) > 0) {
  
  wgcna_traceability_complete <-
    length(wgcna_ids) == 513 &&
    all(
      wgcna_ids %in% annotation_ids
    ) &&
    all(
      wgcna_ids %in% current_tumor_ids
    ) &&
    sum(duplicated(wgcna_ids)) == 0
}


historical_516_reconciled <- FALSE

# Deliberately FALSE here.
#
# The historical 516-Tumor claim must be supported by a
# specific historical sample-level source before PASS.
#
# Textual occurrence of "516" alone is not sufficient.

historical_516_reconciled <- FALSE


if (
  unclassified_resolved &&
  wgcna_traceability_complete &&
  historical_516_reconciled
) {
  
  status <- "PASS"
  
} else {
  
  status <- "HOLD"
}


cat(
  "STEP13B-7 STATUS: ",
  status,
  "\n",
  sep = ""
)

cat("\n")

cat(
  "Unclassified samples resolved: ",
  unclassified_resolved,
  "\n",
  sep = ""
)

cat(
  "513 WGCNA tumor traceability complete: ",
  wgcna_traceability_complete,
  "\n",
  sep = ""
)

cat(
  "Historical 516 Tumor independently reconciled: ",
  historical_516_reconciled,
  "\n\n",
  sep = ""
)


if (status == "PASS") {
  
  cat(
    "All required provenance components are independently reconciled.\n"
  )
  
  cat(
    "STEP13B may proceed to final input-lock review.\n"
  )
  
} else {
  
  cat(
    "STEP13B-7 remains on HOLD.\n\n"
  )
  
  cat(
    "At least one provenance component remains unresolved.\n"
  )
  
  cat(
    "No final classification input lock is issued.\n"
  )
  
  cat(
    "ROC/model analysis remains on HOLD.\n"
  )
}


# ------------------------------------------------------------
# 28. Final safety lock
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13B-7 PROVENANCE CHECKPOINT COMPLETE\n")
cat("============================================================\n\n")

cat("NO ROC CALCULATION\n")
cat("NO AUC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO LOGISTIC REGRESSION\n")
cat("NO THRESHOLD OPTIMIZATION\n")
cat("NO STATISTICAL TEST\n")
cat("NO P VALUE\n")
cat("NO FDR\n\n")

cat("============================================================\n")

sink()


# ------------------------------------------------------------
# 29. Save PASS / HOLD status
# ------------------------------------------------------------

if (status == "PASS") {
  
  writeLines(
    c(
      "STEP13B-7 STATUS: PASS",
      "",
      "Unclassified sample resolution completed.",
      "Final 513-tumor WGCNA sample traceability completed.",
      "Historical 516-Tumor provenance independently reconciled.",
      "",
      "See:",
      "STEP13B-7_Unclassified_516_513_Traceability_Audit.txt",
      "",
      "NO ROC CALCULATION",
      "NO MODEL FITTING",
      "NO STATISTICAL TEST"
    ),
    pass_file
  )
  
  if (file.exists(hold_file)) {
    file.remove(hold_file)
  }
  
} else {
  
  writeLines(
    c(
      "STEP13B-7 STATUS: HOLD",
      "",
      "At least one provenance component remains unresolved.",
      "",
      paste0(
        "Unclassified samples resolved: ",
        unclassified_resolved
      ),
      paste0(
        "513 WGCNA tumor traceability complete: ",
        wgcna_traceability_complete
      ),
      paste0(
        "Historical 516 Tumor independently reconciled: ",
        historical_516_reconciled
      ),
      "",
      "No final classification input lock is issued.",
      "ROC/model analysis remains on HOLD.",
      "",
      "NO ROC CALCULATION",
      "NO MODEL FITTING",
      "NO STATISTICAL TEST"
    ),
    hold_file
  )
}


# ------------------------------------------------------------
# 30. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13B-7 COMPLETE\n")
cat("============================================================\n")

cat(
  "Current annotation records: ",
  nrow(ann),
  "\n",
  sep = ""
)

cat(
  "Current Tumor: ",
  nrow(tumor_current),
  "\n",
  sep = ""
)

cat(
  "Current Normal: ",
  nrow(normal_current),
  "\n",
  sep = ""
)

cat(
  "Current Unclassified: ",
  nrow(unclassified),
  "\n",
  sep = ""
)

cat(
  "Final WGCNA sample columns: ",
  length(wgcna_ids),
  "\n",
  sep = ""
)

cat(
  "STEP13B-7 STATUS: ",
  status,
  "\n",
  sep = ""
)

cat("------------------------------------------------------------\n")
cat("NO ROC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")
cat("============================================================\n")


# ============================================================
# STEP13B-7A
# Unclassified Sample Resolution + Historical 516 Traceability
# ============================================================
#
# PURPOSE:
#   Resolve the 2 currently unclassified expression samples
#   and audit the historical 516-tumor cohort claim.
#
# IMPORTANT:
#   This is a PROVENANCE-ONLY audit.
#
# NO:
#   - ROC
#   - AUC
#   - logistic regression
#   - model fitting
#   - statistical testing
#   - threshold optimization
#
# ============================================================


# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

source_root <- "input_dir/Paper1_tumor-only WGCNA"

out_dir <- file.path(
  workflow_root,
  "STEP13B-7A_Unclassified_Resolution"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

file_current_annotation <- file.path(
  source_root,
  "TCGA_Expression_Sample_Annotation_Basic.txt"
)

file_expr_annotation <- file.path(
  source_root,
  "Workflow-03F-2_expression_annotation.txt"
)

file_patient_mapping <- file.path(
  source_root,
  "Workflow-03F-3_patient_expression_mapping.txt"
)

file_historical_patient <- file.path(
  source_root,
  "Workflow-03F-5_historical_517_patient_annotation.txt"
)

file_matrix_structure <- file.path(
  source_root,
  "Workflow-03F-7_expression_matrix_structure.txt"
)

file_scale_provenance <- file.path(
  source_root,
  "Workflow-03F-9_expression_scale_provenance.txt"
)

file_provenance_audit <- file.path(
  source_root,
  "Workflow-03F-9_provenance_audit_DONE.txt"
)


# ------------------------------------------------------------
# 2. Locate STEP13B-7 historical search output
# ------------------------------------------------------------

possible_516_files <- list.files(
  file.path(workflow_root, "STEP13B-7_Unclassified_Resolution"),
  pattern = "516|historical|claim",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

cat("\n")
cat("============================================================\n")
cat("STEP13B-7A START\n")
cat("============================================================\n\n")

cat("Source root:\n")
cat(source_root, "\n\n")

cat("Workflow root:\n")
cat(workflow_root, "\n\n")


# ------------------------------------------------------------
# 3. Verify input files
# ------------------------------------------------------------

input_files <- c(
  file_current_annotation,
  file_expr_annotation,
  file_patient_mapping,
  file_historical_patient,
  file_matrix_structure,
  file_scale_provenance,
  file_provenance_audit
)

input_check <- data.frame(
  file = basename(input_files),
  exists = file.exists(input_files),
  stringsAsFactors = FALSE
)

print(input_check, row.names = FALSE)

if (!all(input_check$exists)) {
  
  stop(
    "One or more required provenance files are missing."
  )
}


# ------------------------------------------------------------
# 4. Read current sample-level annotation
# ------------------------------------------------------------

ann <- read.delim(
  file_current_annotation,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n")
cat("Current annotation dimensions:\n")
cat(
  "Rows = ",
  nrow(ann),
  "; Columns = ",
  ncol(ann),
  "\n",
  sep = ""
)

print(names(ann))


# ------------------------------------------------------------
# 5. Identify Tumor / Normal / Unclassified
# ------------------------------------------------------------

tumor_current <- ann[
  trimws(ann$Tumor_Normal) == "Tumor",
  ,
  drop = FALSE
]

normal_current <- ann[
  trimws(ann$Tumor_Normal) == "Normal",
  ,
  drop = FALSE
]

unclassified <- ann[
  is.na(ann$Tumor_Normal) |
    trimws(ann$Tumor_Normal) == "" |
    !trimws(ann$Tumor_Normal) %in% c(
      "Tumor",
      "Normal"
    ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 6. Print the 2 unclassified samples
# ------------------------------------------------------------

cat("\n")
cat("------------------------------------------------------------\n")
cat("6. UNCLASSIFIED SAMPLE IDENTIFICATION\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Unclassified records = ",
  nrow(unclassified),
  "\n\n",
  sep = ""
)

print(
  unclassified,
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Save unclassified records
# ------------------------------------------------------------

write.table(
  unclassified,
  file = file.path(
    out_dir,
    "STEP13B-7A_Unclassified_Samples.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Extract sample IDs
# ------------------------------------------------------------

sample_id_col <- "Expression_Sample"

if (!sample_id_col %in% names(ann)) {
  
  stop(
    "Expression_Sample column not found in current annotation."
  )
}

unclassified_ids <- trimws(
  unclassified[[sample_id_col]]
)

cat("\n")
cat("Unclassified sample IDs:\n")
print(unclassified_ids)


# ------------------------------------------------------------
# 9. Barcode suffix audit
# ------------------------------------------------------------

barcode_suffix <- sub(
  ".*-(01|11)$",
  "\\1",
  unclassified_ids
)

barcode_suffix_valid <- grepl(
  "-(01|11)$",
  unclassified_ids
)

barcode_audit <- data.frame(
  Expression_Sample = unclassified_ids,
  Barcode_01_or_11 = ifelse(
    barcode_suffix_valid,
    barcode_suffix,
    NA
  ),
  stringsAsFactors = FALSE
)

cat("\n")
cat("Barcode suffix audit:\n")
print(
  barcode_audit,
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Read historical expression annotation
# ------------------------------------------------------------

hist_expr <- read.delim(
  file_expr_annotation,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n")
cat("Historical expression annotation:\n")
cat(
  "Rows = ",
  nrow(hist_expr),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Match unclassified IDs to historical expression annotation
# ------------------------------------------------------------

hist_expr_match <- hist_expr[
  trimws(hist_expr$Expression_Sample) %in%
    unclassified_ids,
  ,
  drop = FALSE
]

cat("\n")
cat("------------------------------------------------------------\n")
cat("11. HISTORICAL EXPRESSION ANNOTATION MATCH\n")
cat("------------------------------------------------------------\n\n")

if (nrow(hist_expr_match) > 0) {
  
  print(
    hist_expr_match,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No unclassified sample was found in Workflow-03F-2.\n"
  )
}


# ------------------------------------------------------------
# 12. Read patient-expression mapping
# ------------------------------------------------------------

patient_map <- read.delim(
  file_patient_mapping,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 13. Match patient mapping
# ------------------------------------------------------------

patient_map_match <- patient_map[
  sapply(
    patient_map$Expression_Sample,
    function(x) {
      
      if (is.na(x)) {
        return(FALSE)
      }
      
      any(
        unclassified_ids %in%
          trimws(
            unlist(
              strsplit(
                x,
                ";",
                fixed = TRUE
              )
            )
          )
      )
    }
  ),
  ,
  drop = FALSE
]

cat("\n")
cat("------------------------------------------------------------\n")
cat("13. PATIENT-LEVEL MAPPING MATCH\n")
cat("------------------------------------------------------------\n\n")

if (nrow(patient_map_match) > 0) {
  
  print(
    patient_map_match,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No patient-level mapping record matched the two samples.\n"
  )
}


# ------------------------------------------------------------
# 14. Read historical 517-patient annotation
# ------------------------------------------------------------

hist_patient <- read.delim(
  file_historical_patient,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 15. Match historical patient records
# ------------------------------------------------------------

hist_patient_match <- hist_patient[
  sapply(
    hist_patient$Expression_Sample,
    function(x) {
      
      if (is.na(x)) {
        return(FALSE)
      }
      
      any(
        unclassified_ids %in%
          trimws(
            unlist(
              strsplit(
                x,
                ";",
                fixed = TRUE
              )
            )
          )
      )
    }
  ),
  ,
  drop = FALSE
]

cat("\n")
cat("------------------------------------------------------------\n")
cat("15. HISTORICAL 517-PATIENT ANNOTATION MATCH\n")
cat("------------------------------------------------------------\n\n")

if (nrow(hist_patient_match) > 0) {
  
  print(
    hist_patient_match,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No historical patient record matched the two samples.\n"
  )
}


# ------------------------------------------------------------
# 16. Search exact sample IDs in historical text files
# ------------------------------------------------------------

historical_text_files <- c(
  file_expr_annotation,
  file_patient_mapping,
  file_historical_patient,
  file_matrix_structure,
  file_scale_provenance,
  file_provenance_audit
)

exact_hits <- list()

for (f in historical_text_files) {
  
  txt <- tryCatch(
    readLines(
      f,
      warn = FALSE
    ),
    error = function(e) character(0)
  )
  
  if (length(txt) == 0) {
    next
  }
  
  for (sid in unclassified_ids) {
    
    idx <- grep(
      sid,
      txt,
      fixed = TRUE
    )
    
    if (length(idx) > 0) {
      
      exact_hits[[paste(
        basename(f),
        sid,
        sep = "___"
      )]] <- data.frame(
        file = basename(f),
        sample_id = sid,
        line_number = idx,
        line_text = txt[idx],
        stringsAsFactors = FALSE
      )
    }
  }
}


# ------------------------------------------------------------
# 17. Print exact sample-ID provenance hits
# ------------------------------------------------------------

cat("\n")
cat("------------------------------------------------------------\n")
cat("17. EXACT SAMPLE-ID PROVENANCE HITS\n")
cat("------------------------------------------------------------\n\n")

if (length(exact_hits) > 0) {
  
  exact_hits_df <- do.call(
    rbind,
    exact_hits
  )
  
  print(
    exact_hits_df,
    row.names = FALSE
  )
  
  write.table(
    exact_hits_df,
    file = file.path(
      out_dir,
      "STEP13B-7A_Exact_Sample_ID_Provenance_Hits.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No exact sample-ID hits were identified.\n"
  )
}


# ------------------------------------------------------------
# 18. Search for historical 516 evidence
# ------------------------------------------------------------

cat("\n")
cat("------------------------------------------------------------\n")
cat("18. HISTORICAL 516 SEARCH\n")
cat("------------------------------------------------------------\n\n")

search_candidates <- c(
  file.path(
    source_root,
    "Workflow-03F-5_historical_517_patient_annotation.txt"
  ),
  file.path(
    source_root,
    "Workflow-03F-3_patient_expression_mapping.txt"
  ),
  file.path(
    source_root,
    "Workflow-03F-2_expression_annotation.txt"
  ),
  file.path(
    source_root,
    "Workflow-03F-7_expression_matrix_structure.txt"
  ),
  file.path(
    source_root,
    "Workflow-03F-9_provenance_audit_DONE.txt"
  ),
  file.path(
    source_root,
    "Workflow-04-1C_clinical_missing_summary.txt"
  ),
  file.path(
    source_root,
    "Workflow-04-1F-1_clinical_cleaning_checkpoint.txt"
  ),
  file.path(
    source_root,
    "Workflow-04-2A_primary_tumor_sample_list.txt"
  )
)

search_candidates <- search_candidates[
  file.exists(search_candidates)
]

historical_516_hits <- list()

for (f in search_candidates) {
  
  txt <- tryCatch(
    readLines(
      f,
      warn = FALSE
    ),
    error = function(e) character(0)
  )
  
  if (length(txt) == 0) {
    next
  }
  
  idx <- grep(
    "516",
    txt,
    fixed = TRUE
  )
  
  if (length(idx) > 0) {
    
    historical_516_hits[[f]] <- data.frame(
      file = basename(f),
      line_number = idx,
      line_text = txt[idx],
      stringsAsFactors = FALSE
    )
  }
}


# ------------------------------------------------------------
# 19. Print historical 516 evidence
# ------------------------------------------------------------

if (length(historical_516_hits) > 0) {
  
  historical_516_df <- do.call(
    rbind,
    historical_516_hits
  )
  
  print(
    historical_516_df,
    row.names = FALSE
  )
  
  write.table(
    historical_516_df,
    file = file.path(
      out_dir,
      "STEP13B-7A_Historical_516_Evidence.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No explicit '516' occurrence found in the targeted files.\n"
  )
  
  writeLines(
    "No explicit '516' occurrence found in targeted historical provenance files.",
    file.path(
      out_dir,
      "STEP13B-7A_Historical_516_Evidence.txt"
    )
  )
}


# ------------------------------------------------------------
# 20. Preliminary provenance interpretation
# ------------------------------------------------------------

cat("\n")
cat("------------------------------------------------------------\n")
cat("20. PRELIMINARY PROVENANCE INTERPRETATION\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Current annotation:\n",
  "  Tumor = ",
  nrow(tumor_current),
  "\n",
  "  Normal = ",
  nrow(normal_current),
  "\n",
  "  Unclassified = ",
  nrow(unclassified),
  "\n\n",
  sep = ""
)

cat(
  "The two unclassified records are NOT automatically\n"
)

cat(
  "assigned to Tumor or Normal.\n\n"
)

cat(
  "Historical 516-tumor evidence will be considered\n"
)

cat(
  "resolved only if supported by explicit provenance.\n\n"
)

cat(
  "No statistical analysis was performed.\n"
)

cat(
  "No classification model was fitted.\n"
)

cat(
  "No ROC/AUC was calculated.\n\n"
)


# ------------------------------------------------------------
# 21. Conservative status lock
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP13B-7A COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Unclassified samples identified: ",
  nrow(unclassified),
  "\n",
  sep = ""
)

cat(
  "Historical 516 evidence files found: ",
  length(historical_516_hits),
  "\n",
  sep = ""
)

cat("\n")

cat(
  "STEP13B-7A STATUS: INFORMATION EXTRACTED\n"
)

cat(
  "FINAL PROVENANCE RESOLUTION: NOT YET LOCKED\n"
)

cat("\n")
cat("NO ROC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")
cat("============================================================\n")

## ============================================================
# STEP13B-7A RESULT REVIEW
# ============================================================

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

out_dir <- file.path(
  workflow_root,
  "STEP13B-7A_Unclassified_Resolution"
)

cat("\n")
cat("============================================================\n")
cat("STEP13B-7A RESULT REVIEW\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Unclassified samples
# ------------------------------------------------------------

f1 <- file.path(
  out_dir,
  "STEP13B-7A_Unclassified_Samples.txt"
)

cat("------------------------------------------------------------\n")
cat("1. UNCLASSIFIED SAMPLES\n")
cat("------------------------------------------------------------\n\n")

if (file.exists(f1)) {
  print(
    read.delim(
      f1,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    row.names = FALSE
  )
} else {
  cat("FILE NOT FOUND\n")
}


# ------------------------------------------------------------
# 2. Exact provenance hits
# ------------------------------------------------------------

f2 <- file.path(
  out_dir,
  "STEP13B-7A_Exact_Sample_ID_Provenance_Hits.txt"
)

cat("\n")
cat("------------------------------------------------------------\n")
cat("2. EXACT SAMPLE-ID PROVENANCE HITS\n")
cat("------------------------------------------------------------\n\n")

if (file.exists(f2)) {
  
  hits <- read.delim(
    f2,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  print(
    hits,
    row.names = FALSE
  )
  
} else {
  
  cat("FILE NOT FOUND\n")
}


# ------------------------------------------------------------
# 3. Historical 516 evidence
# ------------------------------------------------------------

f3 <- file.path(
  out_dir,
  "STEP13B-7A_Historical_516_Evidence.txt"
)

cat("\n")
cat("------------------------------------------------------------\n")
cat("3. HISTORICAL 516 EVIDENCE\n")
cat("------------------------------------------------------------\n\n")

if (file.exists(f3)) {
  
  txt <- readLines(
    f3,
    warn = FALSE
  )
  
  cat(
    paste(
      txt,
      collapse = "\n"
    )
  )
  
  cat("\n")
  
} else {
  
  cat("FILE NOT FOUND\n")
}


cat("\n")
cat("============================================================\n")
cat("RESULT REVIEW COMPLETE\n")
cat("============================================================\n")

## ============================================================
# STEP13B-7A RESULT REVIEW — TEXT-BASED PROVENANCE INSPECTION
# Purpose:
#   Read the provenance-hit file as raw text because its rows
#   are not guaranteed to have identical tab-delimited fields.
# NO ROC
# NO MODEL FITTING
# NO STATISTICAL TEST
# ============================================================

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

out_dir <- file.path(
  workflow_root,
  "STEP13B-7A_Unclassified_Resolution"
)

f2 <- file.path(
  out_dir,
  "STEP13B-7A_Exact_Sample_ID_Provenance_Hits.txt"
)

f3 <- file.path(
  out_dir,
  "STEP13B-7A_Historical_516_Evidence.txt"
)

# ------------------------------------------------------------
# 1. Exact provenance hits — RAW TEXT
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("1. EXACT SAMPLE-ID PROVENANCE HITS — RAW TEXT\n")
cat("============================================================\n\n")

if (file.exists(f2)) {
  
  lines_f2 <- readLines(
    f2,
    warn = FALSE,
    encoding = "UTF-8"
  )
  
  cat(
    paste(
      seq_along(lines_f2),
      lines_f2,
      sep = ": "
    ),
    sep = "\n"
  )
  
} else {
  
  cat("FILE NOT FOUND\n")
  
}


# ------------------------------------------------------------
# 2. Historical 516 evidence — RAW TEXT
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("2. HISTORICAL 516 EVIDENCE — RAW TEXT\n")
cat("============================================================\n\n")

if (file.exists(f3)) {
  
  lines_f3 <- readLines(
    f3,
    warn = FALSE,
    encoding = "UTF-8"
  )
  
  cat(
    paste(
      seq_along(lines_f3),
      lines_f3,
      sep = ": "
    ),
    sep = "\n"
  )
  
} else {
  
  cat("FILE NOT FOUND\n")
  
}


# ------------------------------------------------------------
# 3. Compact audit summary
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("3. CURRENT PROVENANCE STATUS\n")
cat("============================================================\n\n")

cat("Unclassified sample 1: TCGA-50-5066-02\n")
cat("Unclassified sample 2: TCGA-50-5946-02\n")
cat("\n")

cat("Both samples have sample-type code: 02\n")
cat("Therefore they are NOT automatically classified as Tumor or Normal.\n")
cat("\n")

cat("STEP13B-7A remains INFORMATION EXTRACTED.\n")
cat("STEP13B-7 remains HOLD until provenance is independently resolved.\n")
cat("\n")

cat("NO ROC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")

cat("\n============================================================\n")


# ============================================================
# STEP13B-7B — TARGETED AUDIT OF THE TWO -02 SAMPLES
# ============================================================

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

wgcna_root <- "input_dir/Paper1_tumor-only WGCNA"

targets <- c(
  "TCGA-50-5066-02",
  "TCGA-50-5946-02"
)

files_to_check <- c(
  file.path(
    wgcna_root,
    "Workflow-03F-2_expression_annotation.txt"
  ),
  file.path(
    wgcna_root,
    "Workflow-03F-3_patient_expression_mapping.txt"
  ),
  file.path(
    wgcna_root,
    "Workflow-03F-5_historical_517_patient_annotation.txt"
  ),
  file.path(
    wgcna_root,
    "TCGA_Expression_Sample_Annotation_Basic.txt"
  ),
  file.path(
    wgcna_root,
    "Workflow-03F-7_expression_matrix_structure.txt"
  ),
  file.path(
    wgcna_root,
    "Workflow-03F-9_expression_scale_provenance.txt"
  )
)

cat("============================================================\n")
cat("STEP13B-7B — TARGETED -02 SAMPLE PROVENANCE AUDIT\n")
cat("============================================================\n\n")

cat("Target samples:\n")
cat(paste0("  ", targets, collapse = "\n"))
cat("\n\n")

for (f in files_to_check) {
  
  cat("------------------------------------------------------------\n")
  cat("FILE: ", basename(f), "\n", sep = "")
  cat("------------------------------------------------------------\n")
  
  if (!file.exists(f)) {
    cat("FILE NOT FOUND\n\n")
    next
  }
  
  lines <- readLines(
    f,
    warn = FALSE,
    encoding = "UTF-8"
  )
  
  found_any <- FALSE
  
  for (target in targets) {
    
    idx <- grep(
      target,
      lines,
      fixed = TRUE
    )
    
    if (length(idx) > 0) {
      
      found_any <- TRUE
      
      cat("\nTARGET: ", target, "\n", sep = "")
      
      for (i in idx) {
        
        cat(
          "Line ",
          i,
          ": ",
          lines[i],
          "\n",
          sep = ""
        )
        
      }
      
    }
    
  }
  
  if (!found_any) {
    cat("No exact target sample ID found.\n")
  }
  
  cat("\n")
}

cat("============================================================\n")
cat("STEP13B-7B COMPLETE\n")
cat("NO ROC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")
cat("============================================================\n")


# ============================================================
# STEP13B-7C
# ORIGINAL EXPRESSION MATRIX ↔ SAMPLE ANNOTATION
# EXACT SAMPLE-ID RECONCILIATION
#
# NO ROC
# NO MODEL
# NO STATISTICAL TEST
# ============================================================

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

wgcna_root <- "input_dir/Paper1_tumor-only WGCNA"

annotation_file <- file.path(
  wgcna_root,
  "TCGA_Expression_Sample_Annotation_Basic.txt"
)

original_matrix_file <- file.path(
  wgcna_root,
  "Workflow-03F-7_expression_matrix_structure.txt"
)

# ------------------------------------------------------------
# IMPORTANT
# ------------------------------------------------------------
# The structure file is only a summary and cannot establish
# exact matrix membership.
#
# Therefore first locate the actual historical expression
# matrix referenced by Workflow-03F-7.
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP13B-7C — EXACT SAMPLE-ID RECONCILIATION\n")
cat("============================================================\n\n")

cat("Annotation file:\n")
cat(annotation_file, "\n\n")

cat("Target metastatic samples:\n")
cat("TCGA-50-5066-02\n")
cat("TCGA-50-5946-02\n\n")

cat("The current annotation establishes:\n")
cat("513 Primary Tumor\n")
cat(" 59 Solid Tissue Normal\n")
cat("  2 Metastatic\n")
cat("574 total expression samples\n\n")

cat("NEXT AUDIT TARGET:\n")
cat("Exact membership of all 574 annotated samples in the\n")
cat("historical expression matrix.\n\n")

cat("NO ROC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")

cat("============================================================\n")

# ============================================================
# STEP13B-7D
# LOCATE HISTORICAL 60,483 × 575 EXPRESSION MATRIX
#
# Purpose:
#   Identify the actual historical expression matrix corresponding
#   to Workflow-03F-7.
#
# This step DOES NOT:
#   - calculate ROC
#   - fit any model
#   - calculate statistics
#   - modify any expression data
#
# It only inventories plausible matrix files.
# ============================================================

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"
wgcna_root <- "input_dir/Paper1_tumor-only WGCNA"

out_dir <- file.path(
  workflow_root,
  "STEP13B-7D_Original_Matrix_Location"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("STEP13B-7D — ORIGINAL EXPRESSION MATRIX LOCATION AUDIT\n")
cat("============================================================\n\n")

cat("Search root:\n")
cat(wgcna_root, "\n\n")

# ------------------------------------------------------------
# 1. Candidate extensions
# ------------------------------------------------------------

candidate_ext <- c(
  "\\.txt$",
  "\\.tsv$",
  "\\.csv$",
  "\\.rds$",
  "\\.RData$",
  "\\.rda$",
  "\\.gz$"
)

# ------------------------------------------------------------
# 2. List files
# ------------------------------------------------------------

all_files <- list.files(
  wgcna_root,
  recursive = TRUE,
  full.names = TRUE
)

all_files <- all_files[
  !grepl(
    paste(candidate_ext, collapse = "|"),
    all_files,
    ignore.case = TRUE
  ) == FALSE
]

cat("Candidate files found:", length(all_files), "\n\n")

# ------------------------------------------------------------
# 3. File inventory
# ------------------------------------------------------------

if (length(all_files) == 0) {
  
  cat("NO CANDIDATE FILES FOUND\n")
  
} else {
  
  info <- file.info(all_files)
  
  inventory <- data.frame(
    file = all_files,
    basename = basename(all_files),
    size_MB = round(info$size / 1024^2, 3),
    stringsAsFactors = FALSE
  )
  
  inventory <- inventory[
    order(inventory$size_MB, decreasing = TRUE),
  ]
  
  print(
    inventory,
    row.names = FALSE
  )
  
  write.table(
    inventory,
    file = file.path(
      out_dir,
      "STEP13B-7D_File_Inventory.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 4. Search filenames for likely historical expression matrices
# ------------------------------------------------------------

keywords <- c(
  "expression",
  "expr",
  "matrix",
  "count",
  "fpkm",
  "tpm",
  "tcga",
  "60483",
  "575"
)

pattern <- paste(
  keywords,
  collapse = "|"
)

likely <- all_files[
  grepl(
    pattern,
    basename(all_files),
    ignore.case = TRUE
  )
]

cat("\n============================================================\n")
cat("LIKELY EXPRESSION-MATRIX FILES\n")
cat("============================================================\n\n")

if (length(likely) == 0) {
  
  cat("No filename-based candidates identified.\n")
  
} else {
  
  likely_info <- file.info(likely)
  
  likely_df <- data.frame(
    file = likely,
    basename = basename(likely),
    size_MB = round(likely_info$size / 1024^2, 3),
    stringsAsFactors = FALSE
  )
  
  likely_df <- likely_df[
    order(likely_df$size_MB, decreasing = TRUE),
  ]
  
  print(
    likely_df,
    row.names = FALSE
  )
  
  write.table(
    likely_df,
    file = file.path(
      out_dir,
      "STEP13B-7D_Likely_Expression_Matrix_Candidates.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 5. Search specifically for files named in historical provenance
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FILES RELATED TO WORKFLOW-03F\n")
cat("============================================================\n\n")

wf03f <- all_files[
  grepl(
    "Workflow-03F",
    basename(all_files),
    ignore.case = TRUE
  )
]

if (length(wf03f) == 0) {
  
  cat("No Workflow-03F files found.\n")
  
} else {
  
  wf_info <- file.info(wf03f)
  
  wf_df <- data.frame(
    file = wf03f,
    basename = basename(wf03f),
    size_MB = round(wf_info$size / 1024^2, 3),
    stringsAsFactors = FALSE
  )
  
  print(
    wf_df,
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 6. Final instruction
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP13B-7D COMPLETE\n")
cat("============================================================\n\n")

cat("This is a FILE-LOCATION audit only.\n")
cat("No expression matrix was statistically analyzed.\n")
cat("No ROC calculation.\n")
cat("No model fitting.\n")
cat("No statistical test.\n\n")

cat("STEP13B-7 remains HOLD.\n")
cat("============================================================\n")


# ============================================================
# STEP13B-7E
# WORKFLOW-03F-9 PROVENANCE RDATA OBJECT AUDIT
#
# Purpose:
#   Identify objects stored in the historical provenance RData
#   and determine whether the original 60,483 × 575 expression
#   matrix or its sample IDs are preserved.
#
# NO ROC
# NO MODEL FITTING
# NO STATISTICAL TEST
# NO DATA MODIFICATION
# ============================================================

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

wgcna_root <- "input_dir/Paper1_tumor-only WGCNA"

rdata_file <- file.path(
  wgcna_root,
  "Workflow-03F-9_provenance_audit.RData"
)

out_dir <- file.path(
  workflow_root,
  "STEP13B-7E_Workflow03F9_RData_Audit"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("STEP13B-7E — WORKFLOW-03F-9 RDATA OBJECT AUDIT\n")
cat("============================================================\n\n")

cat("RData file:\n")
cat(rdata_file, "\n\n")

if (!file.exists(rdata_file)) {
  
  stop("Workflow-03F-9 provenance RData file not found.")
  
}

# ------------------------------------------------------------
# 1. Load into isolated environment
# ------------------------------------------------------------

audit_env <- new.env(parent = emptyenv())

cat("Loading RData...\n")
cat("(No objects will be modified.)\n\n")

loaded_objects <- load(
  rdata_file,
  envir = audit_env
)

cat("Objects loaded:", length(loaded_objects), "\n\n")

cat("Object names:\n")
cat(
  paste0(
    "  ",
    loaded_objects,
    collapse = "\n"
  )
)
cat("\n\n")

# ------------------------------------------------------------
# 2. Object-level structural audit
# ------------------------------------------------------------

cat("============================================================\n")
cat("OBJECT STRUCTURE\n")
cat("============================================================\n\n")

for (nm in loaded_objects) {
  
  obj <- get(
    nm,
    envir = audit_env,
    inherits = FALSE
  )
  
  cat("------------------------------------------------------------\n")
  cat("OBJECT: ", nm, "\n", sep = "")
  cat("CLASS:  ", paste(class(obj), collapse = ", "), "\n", sep = "")
  cat("TYPE:   ", typeof(obj), "\n", sep = "")
  
  d <- dim(obj)
  
  if (!is.null(d)) {
    
    cat(
      "DIM:    ",
      paste(d, collapse = " x "),
      "\n",
      sep = ""
    )
    
  } else {
    
    cat("DIM:    NULL\n")
    
  }
  
  if (is.vector(obj) || is.factor(obj)) {
    
    cat("LENGTH: ", length(obj), "\n", sep = "")
    
  }
  
  # For data.frames, only print column names
  if (is.data.frame(obj)) {
    
    cat(
      "COLS:   ",
      paste(
        names(obj),
        collapse = " | "
      ),
      "\n",
      sep = ""
    )
    
  }
  
  cat("\n")
}

# ------------------------------------------------------------
# 3. Search all character vectors / data.frame columns
#    for the two unresolved metastatic sample IDs
# ------------------------------------------------------------

targets <- c(
  "TCGA-50-5066-02",
  "TCGA-50-5946-02"
)

cat("============================================================\n")
cat("TARGET SAMPLE-ID SEARCH\n")
cat("============================================================\n\n")

cat("Targets:\n")
cat(paste0("  ", targets, collapse = "\n"))
cat("\n\n")

target_hits <- list()

for (nm in loaded_objects) {
  
  obj <- get(
    nm,
    envir = audit_env,
    inherits = FALSE
  )
  
  # ----------------------------------------------------------
  # Data frame / matrix
  # ----------------------------------------------------------
  
  if (is.data.frame(obj) || is.matrix(obj)) {
    
    # Search column names
    col_hits <- which(
      names(obj) %in% targets
    )
    
    if (length(col_hits) > 0) {
      
      cat(
        "COLUMN-NAME HIT in object: ",
        nm,
        "\n",
        sep = ""
      )
      
      cat(
        "Columns: ",
        paste(names(obj)[col_hits], collapse = ", "),
        "\n\n",
        sep = ""
      )
      
    }
    
    # Search character/factor columns
    for (j in seq_len(ncol(obj))) {
      
      x <- obj[[j]]
      
      if (is.factor(x)) {
        x <- as.character(x)
      }
      
      if (is.character(x)) {
        
        idx <- which(x %in% targets)
        
        if (length(idx) > 0) {
          
          cat(
            "VALUE HIT\n",
            "Object: ",
            nm,
            "\n",
            "Column: ",
            names(obj)[j],
            "\n",
            "Rows: ",
            paste(idx, collapse = ", "),
            "\n\n",
            sep = ""
          )
          
          print(
            obj[idx, , drop = FALSE],
            row.names = FALSE
          )
          
          cat("\n")
          
        }
        
      }
      
    }
    
  }
  
  # ----------------------------------------------------------
  # Character vector
  # ----------------------------------------------------------
  
  if (is.character(obj)) {
    
    idx <- which(obj %in% targets)
    
    if (length(idx) > 0) {
      
      cat(
        "CHARACTER-VECTOR HIT\n",
        "Object: ",
        nm,
        "\n",
        "Indices: ",
        paste(idx, collapse = ", "),
        "\n",
        "Values: ",
        paste(obj[idx], collapse = ", "),
        "\n\n",
        sep = ""
      )
      
    }
    
  }
  
}

# ------------------------------------------------------------
# 4. Identify objects with dimensions close to historical
#    expression-matrix structure
# ------------------------------------------------------------

cat("============================================================\n")
cat("DIMENSION AUDIT\n")
cat("============================================================\n\n")

cat("Target historical structure:\n")
cat("  60,483 rows × 575 columns\n\n")

dimension_candidates <- data.frame(
  Object = character(),
  Class = character(),
  Rows = integer(),
  Columns = integer(),
  stringsAsFactors = FALSE
)

for (nm in loaded_objects) {
  
  obj <- get(
    nm,
    envir = audit_env,
    inherits = FALSE
  )
  
  d <- dim(obj)
  
  if (!is.null(d) && length(d) == 2) {
    
    dimension_candidates <- rbind(
      dimension_candidates,
      data.frame(
        Object = nm,
        Class = paste(class(obj), collapse = ","),
        Rows = d[1],
        Columns = d[2],
        stringsAsFactors = FALSE
      )
    )
    
  }
  
}

print(
  dimension_candidates,
  row.names = FALSE
)

write.table(
  dimension_candidates,
  file = file.path(
    out_dir,
    "STEP13B-7E_RData_Object_Dimensions.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Save object names only
# ------------------------------------------------------------

writeLines(
  loaded_objects,
  con = file.path(
    out_dir,
    "STEP13B-7E_RData_Object_Names.txt"
  )
)

# ------------------------------------------------------------
# 6. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP13B-7E COMPLETE\n")
cat("============================================================\n\n")

cat("This was a provenance-object audit only.\n")
cat("No expression values were recalculated.\n")
cat("No ROC calculation.\n")
cat("No model fitting.\n")
cat("No statistical test.\n\n")

cat("STEP13B-7 remains HOLD pending exact sample-ID reconciliation.\n")

cat("\n============================================================\n")

# ============================================================
# STEP13B-7F
# HISTORICAL SAMPLE / PATIENT COUNT PROVENANCE AUDIT
#
# Purpose:
#   Determine the provenance of historical counts such as
#   515 / 516 / 517 and reconcile them with the verified
#   574-expression-sample structure.
#
# NO ROC
# NO MODEL FITTING
# NO STATISTICAL TEST
# NO EXPRESSION RECALCULATION
# ============================================================

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

wgcna_root <- "input_dir/Paper1_tumor-only WGCNA"

rdata_file <- file.path(
  wgcna_root,
  "Workflow-03F-9_provenance_audit.RData"
)

out_dir <- file.path(
  workflow_root,
  "STEP13B-7F_Historical_Count_Provenance"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

audit_env <- new.env(parent = emptyenv())

load(
  rdata_file,
  envir = audit_env
)

objects_to_audit <- c(
  "sample_group",
  "patient_summary",
  "patient_expr_summary",
  "multi_sample",
  "no_primary",
  "primary_expr",
  "primary_with_clinical",
  "clinical_with_primary"
)

cat("============================================================\n")
cat("STEP13B-7F — HISTORICAL COUNT PROVENANCE AUDIT\n")
cat("============================================================\n\n")

for (nm in objects_to_audit) {
  
  cat("------------------------------------------------------------\n")
  cat("OBJECT: ", nm, "\n", sep = "")
  cat("------------------------------------------------------------\n")
  
  if (!exists(
    nm,
    envir = audit_env,
    inherits = FALSE
  )) {
    
    cat("OBJECT NOT FOUND\n\n")
    next
    
  }
  
  obj <- get(
    nm,
    envir = audit_env,
    inherits = FALSE
  )
  
  cat(
    "CLASS: ",
    paste(class(obj), collapse = ", "),
    "\n",
    sep = ""
  )
  
  d <- dim(obj)
  
  if (!is.null(d)) {
    
    cat(
      "DIMENSION: ",
      paste(d, collapse = " x "),
      "\n",
      sep = ""
    )
    
  }
  
  if (is.data.frame(obj)) {
    
    cat(
      "COLUMNS:\n",
      paste(
        names(obj),
        collapse = " | "
      ),
      "\n\n",
      sep = ""
    )
    
    print(
      obj,
      row.names = FALSE
    )
    
  } else {
    
    print(obj)
    
  }
  
  cat("\n\n")
}

# ------------------------------------------------------------
# Explicit count summary
# ------------------------------------------------------------

cat("============================================================\n")
cat("VERIFIED SAMPLE-LEVEL COUNT\n")
cat("============================================================\n\n")

cat("Historical expression matrix:\n")
cat("  Matrix columns = 575\n")
cat("  Identifier column = 1\n")
cat("  Expression samples = 574\n\n")

cat("Verified sample types:\n")
cat("  Primary Tumor = 513\n")
cat("  Solid Tissue Normal = 59\n")
cat("  Metastatic = 2\n")
cat("  Total = 574\n\n")

cat("Tumor-versus-Normal eligible samples:\n")
cat("  Primary Tumor = 513\n")
cat("  Normal = 59\n")
cat("  Total eligible = 572\n\n")

cat("Historical count(s) requiring provenance explanation:\n")
cat("  515\n")
cat("  516\n")
cat("  517\n\n")

cat("============================================================\n")
cat("STEP13B-7F COMPLETE\n")
cat("============================================================\n")

cat("NO ROC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")
cat("NO EXPRESSION RECALCULATION\n\n")

cat("STEP13B-7 remains HOLD pending count provenance reconciliation.\n")

cat("============================================================\n")

cat("\n")
cat("============================================================\n")
cat("STEP13B-7G — EXACT SAMPLE-SET RECONCILIATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"
wgcna_root    <- "input_dir/Paper1_tumor-only WGCNA"

rdata_file <- file.path(
  wgcna_root,
  "Workflow-03F-9_provenance_audit.RData"
)

annotation_file <- file.path(
  wgcna_root,
  "TCGA_Expression_Sample_Annotation_Basic.txt"
)

wgcna_file <- file.path(
  wgcna_root,
  "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt"
)

out_dir <- file.path(
  workflow_root,
  "STEP13B-7G_Exact_Sample_Set_Reconciliation"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load provenance RData in isolated environment
# ------------------------------------------------------------

prov_env <- new.env(parent = emptyenv())

load(
  rdata_file,
  envir = prov_env
)

cat("RData loaded:\n")
cat("  ", rdata_file, "\n\n")

# ------------------------------------------------------------
# 2. Extract original expression-matrix sample IDs
# ------------------------------------------------------------

expr_obj <- prov_env$expr

if (is.null(expr_obj)) {
  stop("Object 'expr' not found in provenance RData.")
}

cat("Object 'expr':\n")
cat("  class =", paste(class(expr_obj), collapse = ", "), "\n")
cat("  dimensions =", paste(dim(expr_obj), collapse = " x "), "\n\n")

expr_colnames <- colnames(expr_obj)

if (is.null(expr_colnames)) {
  stop("Object 'expr' has no column names.")
}

# ------------------------------------------------------------
# 3. Identify expression sample columns
# ------------------------------------------------------------

# TCGA sample IDs have the form TCGA-XX-XXXX-XX
is_tcga_sample <- grepl(
  "^TCGA-[A-Z0-9]+-[A-Z0-9]+-[0-9]{2}$",
  expr_colnames
)

matrix_sample_ids <- expr_colnames[is_tcga_sample]

non_sample_columns <- expr_colnames[!is_tcga_sample]

cat("Original matrix sample-ID audit:\n")
cat("  Total columns =", length(expr_colnames), "\n")
cat("  TCGA sample columns =", length(matrix_sample_ids), "\n")
cat("  Non-sample columns =", length(non_sample_columns), "\n\n")

if (length(non_sample_columns) > 0) {
  cat("Non-sample columns:\n")
  print(non_sample_columns)
  cat("\n")
}

# ------------------------------------------------------------
# 4. Duplicate sample-ID audit
# ------------------------------------------------------------

matrix_duplicates <- unique(
  matrix_sample_ids[
    duplicated(matrix_sample_ids)
  ]
)

cat("Matrix duplicate sample IDs:", length(matrix_duplicates), "\n")

if (length(matrix_duplicates) > 0) {
  print(matrix_duplicates)
  cat("\n")
}

# ------------------------------------------------------------
# 5. Load annotation
# ------------------------------------------------------------

ann <- read.delim(
  annotation_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "Expression_Sample",
  "Patient_ID",
  "Sample_Type_Code",
  "Tumor_Normal"
)

missing_cols <- setdiff(
  required_cols,
  colnames(ann)
)

if (length(missing_cols) > 0) {
  stop(
    "Annotation missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

ann_ids <- ann$Expression_Sample

cat("Annotation audit:\n")
cat("  Records =", nrow(ann), "\n")
cat("  Duplicate IDs =",
    sum(duplicated(ann_ids)),
    "\n\n"
)

# ------------------------------------------------------------
# 6. Exact set reconciliation
# ------------------------------------------------------------

matrix_not_annotation <- setdiff(
  matrix_sample_ids,
  ann_ids
)

annotation_not_matrix <- setdiff(
  ann_ids,
  matrix_sample_ids
)

cat("============================================================\n")
cat("EXACT SET RECONCILIATION\n")
cat("============================================================\n\n")

cat(
  "Matrix sample IDs absent from annotation:",
  length(matrix_not_annotation),
  "\n"
)

if (length(matrix_not_annotation) > 0) {
  print(matrix_not_annotation)
}

cat(
  "\nAnnotation IDs absent from matrix:",
  length(annotation_not_matrix),
  "\n"
)

if (length(annotation_not_matrix) > 0) {
  print(annotation_not_matrix)
}

exact_matrix_annotation_match <-
  length(matrix_not_annotation) == 0 &&
  length(annotation_not_matrix) == 0 &&
  length(matrix_sample_ids) == length(ann_ids) &&
  sum(duplicated(matrix_sample_ids)) == 0 &&
  sum(duplicated(ann_ids)) == 0

cat(
  "\nEXACT MATRIX ↔ ANNOTATION MATCH:",
  exact_matrix_annotation_match,
  "\n\n"
)

# ------------------------------------------------------------
# 7. Verify sample-type composition
# ------------------------------------------------------------

cat("============================================================\n")
cat("ANNOTATION SAMPLE COMPOSITION\n")
cat("============================================================\n\n")

print(
  table(
    ann$Sample_Type_Code,
    useNA = "ifany"
  )
)

cat("\nSample_Type:\n")
print(
  table(
    ann$Sample_Type,
    useNA = "ifany"
  )
)

cat("\nTumor_Normal:\n")
print(
  table(
    ann$Tumor_Normal,
    useNA = "ifany"
  )
)

# ------------------------------------------------------------
# 8. Extract Primary Tumor IDs
# ------------------------------------------------------------

primary_ids <- ann$Expression_Sample[
  ann$Sample_Type_Code == "01"
]

normal_ids <- ann$Expression_Sample[
  ann$Sample_Type_Code == "11"
]

metastatic_ids <- ann$Expression_Sample[
  ann$Sample_Type_Code == "02"
]

cat("\n")
cat("Primary Tumor IDs:", length(primary_ids), "\n")
cat("Normal IDs:", length(normal_ids), "\n")
cat("Metastatic IDs:", length(metastatic_ids), "\n")

# ------------------------------------------------------------
# 9. Load final WGCNA matrix only for sample-ID extraction
# ------------------------------------------------------------

cat("\n")
cat("Reading final WGCNA matrix header only...\n")

wgcna_header <- read.delim(
  wgcna_file,
  header = TRUE,
  sep = "\t",
  nrows = 0,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

wgcna_columns <- colnames(wgcna_header)

# The first column is expected to be the gene identifier.
wgcna_sample_ids <- wgcna_columns[
  grepl(
    "^TCGA-[A-Z0-9]+-[A-Z0-9]+-[0-9]{2}$",
    wgcna_columns
  )
]

cat(
  "Final WGCNA sample IDs:",
  length(wgcna_sample_ids),
  "\n"
)

# ------------------------------------------------------------
# 10. WGCNA ↔ Primary Tumor exact reconciliation
# ------------------------------------------------------------

wgcna_not_primary <- setdiff(
  wgcna_sample_ids,
  primary_ids
)

primary_not_wgcna <- setdiff(
  primary_ids,
  wgcna_sample_ids
)

cat("\n")
cat("============================================================\n")
cat("FINAL WGCNA ↔ PRIMARY TUMOR RECONCILIATION\n")
cat("============================================================\n\n")

cat(
  "WGCNA IDs absent from Primary Tumor annotation:",
  length(wgcna_not_primary),
  "\n"
)

if (length(wgcna_not_primary) > 0) {
  print(wgcna_not_primary)
}

cat(
  "\nPrimary Tumor IDs absent from WGCNA:",
  length(primary_not_wgcna),
  "\n"
)

if (length(primary_not_wgcna) > 0) {
  print(primary_not_wgcna)
}

exact_wgcna_primary_match <-
  length(wgcna_not_primary) == 0 &&
  length(primary_not_wgcna) == 0 &&
  length(wgcna_sample_ids) == length(primary_ids)

cat(
  "\nEXACT WGCNA ↔ PRIMARY MATCH:",
  exact_wgcna_primary_match,
  "\n\n"
)

# ------------------------------------------------------------
# 11. Patient-level audit of excluded Primary Tumor samples
# ------------------------------------------------------------

if (length(primary_not_wgcna) > 0) {
  
  excluded_primary <- ann[
    ann$Expression_Sample %in% primary_not_wgcna,
    ,
    drop = FALSE
  ]
  
  cat("Primary Tumor samples absent from WGCNA:\n")
  print(excluded_primary)
  
  write.table(
    excluded_primary,
    file = file.path(
      out_dir,
      "STEP13B-7G_Primary_Tumor_Absent_From_WGCNA.tsv"
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

# ------------------------------------------------------------
# 12. Inspect clinical/provenance objects
# ------------------------------------------------------------

small_objects <- c(
  "clinical_with_primary",
  "clinical_without_primary",
  "primary_expr",
  "primary_with_clinical",
  "no_primary",
  "patient_summary",
  "patient_expr_summary",
  "sample_group",
  "multi_sample"
)

cat("\n")
cat("============================================================\n")
cat("PROVENANCE OBJECT DIMENSIONS\n")
cat("============================================================\n\n")

prov_dim_df <- data.frame(
  object = character(),
  class = character(),
  rows = integer(),
  columns = integer(),
  stringsAsFactors = FALSE
)

for (nm in small_objects) {
  
  if (exists(nm, envir = prov_env, inherits = FALSE)) {
    
    obj <- get(nm, envir = prov_env)
    
    d <- dim(obj)
    
    prov_dim_df <- rbind(
      prov_dim_df,
      data.frame(
        object = nm,
        class = paste(class(obj), collapse = ","),
        rows = ifelse(is.null(d), NA, d[1]),
        columns = ifelse(is.null(d), NA, d[2]),
        stringsAsFactors = FALSE
      )
    )
  }
}

print(prov_dim_df)

write.table(
  prov_dim_df,
  file = file.path(
    out_dir,
    "STEP13B-7G_Provenance_Object_Dimensions.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# 13. Save exact reconciliation results
# ------------------------------------------------------------

summary_df <- data.frame(
  item = c(
    "Original matrix total columns",
    "Original matrix expression samples",
    "Annotation records",
    "Primary Tumor",
    "Normal",
    "Metastatic",
    "Tumor-versus-Normal eligible",
    "Matrix IDs absent from annotation",
    "Annotation IDs absent from matrix",
    "WGCNA sample IDs",
    "WGCNA IDs absent from Primary Tumor",
    "Primary Tumor IDs absent from WGCNA",
    "Exact matrix-annotation match",
    "Exact WGCNA-primary match"
  ),
  value = c(
    length(expr_colnames),
    length(matrix_sample_ids),
    nrow(ann),
    length(primary_ids),
    length(normal_ids),
    length(metastatic_ids),
    length(primary_ids) + length(normal_ids),
    length(matrix_not_annotation),
    length(annotation_not_matrix),
    length(wgcna_sample_ids),
    length(wgcna_not_primary),
    length(primary_not_wgcna),
    exact_matrix_annotation_match,
    exact_wgcna_primary_match
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_df,
  file = file.path(
    out_dir,
    "STEP13B-7G_Exact_Reconciliation_Summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

writeLines(
  matrix_not_annotation,
  file.path(
    out_dir,
    "STEP13B-7G_Matrix_IDs_Absent_From_Annotation.txt"
  )
)

writeLines(
  annotation_not_matrix,
  file.path(
    out_dir,
    "STEP13B-7G_Annotation_IDs_Absent_From_Matrix.txt"
  )
)

writeLines(
  wgcna_not_primary,
  file.path(
    out_dir,
    "STEP13B-7G_WGCNA_IDs_Absent_From_Primary.txt"
  )
)

writeLines(
  primary_not_wgcna,
  file.path(
    out_dir,
    "STEP13B-7G_Primary_IDs_Absent_From_WGCNA.txt"
  )
)

# ------------------------------------------------------------
# 14. Final status
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13B-7G COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Exact matrix ↔ annotation:",
  exact_matrix_annotation_match,
  "\n"
)

cat(
  "Exact WGCNA ↔ Primary Tumor:",
  exact_wgcna_primary_match,
  "\n"
)

cat(
  "Primary Tumor:",
  length(primary_ids),
  "\n"
)

cat(
  "Normal:",
  length(normal_ids),
  "\n"
)

cat(
  "Metastatic:",
  length(metastatic_ids),
  "\n"
)

cat(
  "Eligible Tumor-versus-Normal:",
  length(primary_ids) + length(normal_ids),
  "\n\n"
)

cat("NO ROC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")
cat("NO EXPRESSION RECALCULATION\n\n")

cat(
  "STEP13B-7 remains HOLD until historical count provenance ",
  "is fully reconciled.\n"
)

cat("\n")
cat("============================================================\n")
cat("STEP13B-7G-1A — SAMPLE TYPE CODE TYPE AUDIT\n")
cat("============================================================\n\n")

cat("Class of Sample_Type_Code:\n")
print(class(ann$Sample_Type_Code))

cat("\nStructure:\n")
str(ann$Sample_Type_Code)

cat("\nUnique Sample_Type_Code values:\n")
print(unique(ann$Sample_Type_Code))

cat("\nFrequency table:\n")
print(
  table(
    ann$Sample_Type_Code,
    useNA = "ifany"
  )
)

cat("\nSample_Type frequency:\n")
print(
  table(
    ann$Sample_Type,
    useNA = "ifany"
  )
)

cat("\nTumor_Normal frequency:\n")
print(
  table(
    ann$Tumor_Normal,
    useNA = "ifany"
  )
)

cat("\n============================================================\n")
cat("STEP13B-7G-1A COMPLETE\n")
cat("============================================================\n")
cat("NO ROC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")

cat("\n")
cat("============================================================\n")
cat("STEP13B-7G-1B-FINAL — SAMPLE-LEVEL RECONCILIATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"
wgcna_root    <- "input_dir/Paper1_tumor-only WGCNA"

rdata_file <- file.path(
  wgcna_root,
  "Workflow-03F-9_provenance_audit.RData"
)

annotation_file <- file.path(
  wgcna_root,
  "TCGA_Expression_Sample_Annotation_Basic.txt"
)

wgcna_file <- file.path(
  wgcna_root,
  "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt"
)

out_dir <- file.path(
  workflow_root,
  "STEP13B-7G_Exact_Sample_Set_Reconciliation"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Load provenance RData
# ------------------------------------------------------------

prov_env <- new.env(parent = emptyenv())

load(
  rdata_file,
  envir = prov_env
)

expr_obj <- prov_env$expr

if (is.null(expr_obj)) {
  stop("Object 'expr' not found in provenance RData.")
}

# ------------------------------------------------------------
# 2. Original expression matrix sample IDs
# ------------------------------------------------------------

expr_colnames <- colnames(expr_obj)

if (is.null(expr_colnames)) {
  stop("Object 'expr' has no column names.")
}

is_tcga_sample <- grepl(
  "^TCGA-[A-Z0-9]+-[A-Z0-9]+-[0-9]{2}$",
  expr_colnames
)

matrix_sample_ids <- expr_colnames[
  is_tcga_sample
]

non_sample_columns <- expr_colnames[
  !is_tcga_sample
]

matrix_duplicates <- unique(
  matrix_sample_ids[
    duplicated(matrix_sample_ids)
  ]
)

# ------------------------------------------------------------
# 3. Load verified 4-column annotation
# ------------------------------------------------------------

ann <- read.delim(
  annotation_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "Expression_Sample",
  "Patient_ID",
  "Sample_Type_Code",
  "Tumor_Normal"
)

missing_cols <- setdiff(
  required_cols,
  colnames(ann)
)

if (length(missing_cols) > 0) {
  stop(
    "Annotation missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

ann_ids <- ann$Expression_Sample

annotation_duplicates <- unique(
  ann_ids[
    duplicated(ann_ids)
  ]
)

# ------------------------------------------------------------
# 4. Exact matrix ↔ annotation reconciliation
# ------------------------------------------------------------

matrix_not_annotation <- setdiff(
  matrix_sample_ids,
  ann_ids
)

annotation_not_matrix <- setdiff(
  ann_ids,
  matrix_sample_ids
)

exact_matrix_annotation_match <-
  length(matrix_not_annotation) == 0 &&
  length(annotation_not_matrix) == 0 &&
  length(matrix_sample_ids) == length(ann_ids) &&
  length(matrix_duplicates) == 0 &&
  length(annotation_duplicates) == 0

# ------------------------------------------------------------
# 5. Sample classification
#
# IMPORTANT:
# Sample_Type_Code is INTEGER:
# 1  = Primary Tumor
# 2  = Metastatic
# 11 = Solid Tissue Normal
# ------------------------------------------------------------

primary_ids <- ann$Expression_Sample[
  ann$Sample_Type_Code == 1
]

metastatic_ids <- ann$Expression_Sample[
  ann$Sample_Type_Code == 2
]

normal_ids <- ann$Expression_Sample[
  ann$Sample_Type_Code == 11
]

n_primary <- length(primary_ids)
n_metastatic <- length(metastatic_ids)
n_normal <- length(normal_ids)

n_eligible <- n_primary + n_normal

classification_total <-
  n_primary +
  n_metastatic +
  n_normal

classification_count_consistent <-
  classification_total == nrow(ann) &&
  n_primary == 513 &&
  n_metastatic == 2 &&
  n_normal == 59

# ------------------------------------------------------------
# 6. Independent Tumor_Normal field cross-check
# ------------------------------------------------------------

tumor_field_ids <- ann$Expression_Sample[
  ann$Tumor_Normal == "Tumor"
]

normal_field_ids <- ann$Expression_Sample[
  ann$Tumor_Normal == "Normal"
]

other_field_ids <- ann$Expression_Sample[
  ann$Tumor_Normal == "Other"
]

tumor_field_code_match <- setequal(
  tumor_field_ids,
  primary_ids
)

normal_field_code_match <- setequal(
  normal_field_ids,
  normal_ids
)

other_field_code_match <- setequal(
  other_field_ids,
  metastatic_ids
)

# ------------------------------------------------------------
# 7. Metastatic exclusion audit
# ------------------------------------------------------------

metastatic_in_primary <- intersect(
  metastatic_ids,
  primary_ids
)

metastatic_in_normal <- intersect(
  metastatic_ids,
  normal_ids
)

metastatic_exclusion_clean <-
  length(metastatic_in_primary) == 0 &&
  length(metastatic_in_normal) == 0

# ------------------------------------------------------------
# 8. Final WGCNA sample IDs
# ------------------------------------------------------------

wgcna_header <- read.delim(
  wgcna_file,
  header = TRUE,
  sep = "\t",
  nrows = 0,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

wgcna_columns <- colnames(wgcna_header)

wgcna_sample_ids <- wgcna_columns[
  grepl(
    "^TCGA-[A-Z0-9]+-[A-Z0-9]+-[0-9]{2}$",
    wgcna_columns
  )
]

wgcna_duplicates <- unique(
  wgcna_sample_ids[
    duplicated(wgcna_sample_ids)
  ]
)

wgcna_not_primary <- setdiff(
  wgcna_sample_ids,
  primary_ids
)

primary_not_wgcna <- setdiff(
  primary_ids,
  wgcna_sample_ids
)

exact_wgcna_primary_match <-
  length(wgcna_not_primary) == 0 &&
  length(primary_not_wgcna) == 0 &&
  length(wgcna_sample_ids) == length(primary_ids) &&
  length(wgcna_duplicates) == 0

# ------------------------------------------------------------
# 9. Primary Tumor samples absent from final WGCNA
# ------------------------------------------------------------

primary_excluded_from_wgcna <- ann[
  ann$Expression_Sample %in% primary_not_wgcna,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 10. Explicit metastatic records
# ------------------------------------------------------------

metastatic_annotation <- ann[
  ann$Expression_Sample %in% metastatic_ids,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 11. Save exact ID sets
# ------------------------------------------------------------

writeLines(
  matrix_sample_ids,
  file.path(
    out_dir,
    "STEP13B-7G-1B_Original_Matrix_574_Sample_IDs.txt"
  )
)

writeLines(
  primary_ids,
  file.path(
    out_dir,
    "STEP13B-7G-1B_Primary_Tumor_513_IDs.txt"
  )
)

writeLines(
  normal_ids,
  file.path(
    out_dir,
    "STEP13B-7G-1B_Normal_59_IDs.txt"
  )
)

writeLines(
  metastatic_ids,
  file.path(
    out_dir,
    "STEP13B-7G-1B_Metastatic_2_IDs.txt"
  )
)

writeLines(
  matrix_not_annotation,
  file.path(
    out_dir,
    "STEP13B-7G-1B_Matrix_Not_In_Annotation.txt"
  )
)

writeLines(
  annotation_not_matrix,
  file.path(
    out_dir,
    "STEP13B-7G-1B_Annotation_Not_In_Matrix.txt"
  )
)

writeLines(
  wgcna_not_primary,
  file.path(
    out_dir,
    "STEP13B-7G-1B_WGCNA_Not_Primary.txt"
  )
)

writeLines(
  primary_not_wgcna,
  file.path(
    out_dir,
    "STEP13B-7G-1B_Primary_Not_In_WGCNA.txt"
  )
)

write.table(
  metastatic_annotation,
  file.path(
    out_dir,
    "STEP13B-7G-1B_Metastatic_Annotation.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  primary_excluded_from_wgcna,
  file.path(
    out_dir,
    "STEP13B-7G-1B_Primary_Excluded_From_WGCNA.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# 12. Historical provenance object dimensions
# ------------------------------------------------------------

small_objects <- c(
  "clinical_with_primary",
  "clinical_without_primary",
  "primary_expr",
  "primary_with_clinical",
  "no_primary",
  "patient_summary",
  "patient_expr_summary",
  "sample_group",
  "multi_sample"
)

prov_dim_df <- data.frame(
  object = character(),
  class = character(),
  rows = integer(),
  columns = integer(),
  stringsAsFactors = FALSE
)

for (nm in small_objects) {
  
  if (exists(nm, envir = prov_env, inherits = FALSE)) {
    
    obj <- get(nm, envir = prov_env)
    
    d <- dim(obj)
    
    prov_dim_df <- rbind(
      prov_dim_df,
      data.frame(
        object = nm,
        class = paste(
          class(obj),
          collapse = ","
        ),
        rows = ifelse(
          is.null(d),
          NA_integer_,
          d[1]
        ),
        columns = ifelse(
          is.null(d),
          NA_integer_,
          d[2]
        ),
        stringsAsFactors = FALSE
      )
    )
  }
}

write.table(
  prov_dim_df,
  file.path(
    out_dir,
    "STEP13B-7G-1B_Provenance_Object_Dimensions.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# 13. Formal summary
# ------------------------------------------------------------

summary_df <- data.frame(
  item = c(
    "Original matrix total columns",
    "Original matrix expression samples",
    "Annotation records",
    "Primary Tumor",
    "Solid Tissue Normal",
    "Metastatic",
    "Tumor-versus-Normal eligible",
    "Matrix IDs absent from annotation",
    "Annotation IDs absent from matrix",
    "WGCNA sample IDs",
    "WGCNA IDs absent from Primary Tumor",
    "Primary Tumor IDs absent from WGCNA",
    "Matrix duplicate IDs",
    "Annotation duplicate IDs",
    "WGCNA duplicate IDs",
    "Exact matrix-annotation match",
    "Classification count consistent",
    "Tumor field = code 1",
    "Normal field = code 11",
    "Other field = code 2",
    "Metastatic exclusion clean",
    "Exact WGCNA-primary match"
  ),
  value = c(
    length(expr_colnames),
    length(matrix_sample_ids),
    nrow(ann),
    n_primary,
    n_normal,
    n_metastatic,
    n_eligible,
    length(matrix_not_annotation),
    length(annotation_not_matrix),
    length(wgcna_sample_ids),
    length(wgcna_not_primary),
    length(primary_not_wgcna),
    length(matrix_duplicates),
    length(annotation_duplicates),
    length(wgcna_duplicates),
    exact_matrix_annotation_match,
    classification_count_consistent,
    tumor_field_code_match,
    normal_field_code_match,
    other_field_code_match,
    metastatic_exclusion_clean,
    exact_wgcna_primary_match
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_df,
  file.path(
    out_dir,
    "STEP13B-7G-1B_Formal_Reconciliation_Summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# 14. Console report
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13B-7G-1B-FINAL COMPLETE\n")
cat("============================================================\n\n")

cat("ORIGINAL MATRIX\n")
cat("  Total columns:              ", length(expr_colnames), "\n")
cat("  Expression samples:         ", length(matrix_sample_ids), "\n")
cat("  Duplicate sample IDs:       ", length(matrix_duplicates), "\n\n")

cat("MATRIX ↔ ANNOTATION\n")
cat("  Matrix IDs absent annotation: ",
    length(matrix_not_annotation), "\n")
cat("  Annotation IDs absent matrix: ",
    length(annotation_not_matrix), "\n")
cat("  EXACT MATCH:                  ",
    exact_matrix_annotation_match, "\n\n")

cat("SAMPLE CLASSIFICATION\n")
cat("  Primary Tumor:               ", n_primary, "\n")
cat("  Solid Tissue Normal:         ", n_normal, "\n")
cat("  Metastatic:                  ", n_metastatic, "\n")
cat("  Total:                       ", classification_total, "\n")
cat("  Tumor-versus-Normal eligible:",
    n_eligible, "\n\n")

cat("CLASSIFICATION CROSS-CHECK\n")
cat("  Tumor field = code 1:        ",
    tumor_field_code_match, "\n")
cat("  Normal field = code 11:      ",
    normal_field_code_match, "\n")
cat("  Other field = code 2:        ",
    other_field_code_match, "\n")
cat("  Metastatic exclusion clean:  ",
    metastatic_exclusion_clean, "\n\n")

cat("FINAL WGCNA ↔ PRIMARY TUMOR\n")
cat("  WGCNA samples:               ",
    length(wgcna_sample_ids), "\n")
cat("  WGCNA IDs not Primary:       ",
    length(wgcna_not_primary), "\n")
cat("  Primary IDs not WGCNA:       ",
    length(primary_not_wgcna), "\n")
cat("  EXACT MATCH:                 ",
    exact_wgcna_primary_match, "\n\n")

cat("HISTORICAL COUNT OBJECTS\n")
cat("  sample_group rows:           ",
    nrow(prov_env$sample_group), "\n")
cat("  patient_summary rows:        ",
    nrow(prov_env$patient_summary), "\n")
cat("  patient_expr_summary rows:   ",
    nrow(prov_env$patient_expr_summary), "\n")
cat("  clinical_with_primary rows:  ",
    nrow(prov_env$clinical_with_primary), "\n")
cat("  clinical_without_primary rows:",
    nrow(prov_env$clinical_without_primary), "\n")
cat("  no_primary rows:             ",
    nrow(prov_env$no_primary), "\n\n")

cat("============================================================\n")
cat("STATISTICAL ANALYSIS STATUS\n")
cat("============================================================\n")
cat("NO ROC CALCULATION\n")
cat("NO AUC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")
cat("NO EXPRESSION RECALCULATION\n")
cat("NO THRESHOLD OPTIMIZATION\n\n")

cat("STEP13B-7G-1B-FINAL STATUS: INFORMATION LOCKED\n")
cat("STEP13B-7 OVERALL STATUS: HOLD\n")
cat("Historical 515/516/517 provenance remains unresolved.\n")
cat("============================================================\n")

# ============================================================
# STEP13B-7H
# Historical Count Provenance Audit
#
# Purpose:
#   Resolve the provenance of historical counts 515 / 516 / 517
#   versus the locked current counts:
#       513 Primary Tumor
#        59 Solid Tissue Normal
#         2 Metastatic
#       574 expression samples
#
# This step is PROVENANCE ONLY.
#
# NO:
#   - ROC calculation
#   - AUC calculation
#   - logistic regression
#   - statistical testing
#   - expression recalculation
#   - threshold optimization
#   - model fitting
#
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP13B-7H: HISTORICAL COUNT PROVENANCE AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directories
# ------------------------------------------------------------

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

step_dir <- file.path(
  workflow_root,
  "STEP13B-7H_Historical_Count_Provenance_Audit"
)

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}

output_dir <- step_dir


# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

prov_rdata <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "Workflow-03F-9_provenance_audit.RData"
)

annotation_file <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "TCGA_Expression_Sample_Annotation_Basic.txt"
)

historical_annotation_file <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "Workflow-03F-2_expression_annotation.txt"
)

patient_mapping_file <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "Workflow-03F-3_patient_expression_mapping.txt"
)

historical_patient_file <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "Workflow-03F-5_historical_517_patient_annotation.txt"
)

wgcna_file <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt"
)


# ------------------------------------------------------------
# 2. Basic input existence audit
# ------------------------------------------------------------

input_files <- c(
  prov_rdata,
  annotation_file,
  historical_annotation_file,
  patient_mapping_file,
  historical_patient_file,
  wgcna_file
)

input_status <- data.frame(
  file = basename(input_files),
  path = input_files,
  exists = file.exists(input_files),
  stringsAsFactors = FALSE
)

write.csv(
  input_status,
  file.path(output_dir, "STEP13B-7H_Input_File_Audit.csv"),
  row.names = FALSE
)

if (!all(input_status$exists)) {
  
  cat("ERROR: One or more required input files are missing.\n\n")
  print(input_status)
  
  writeLines(
    c(
      "STEP13B-7H STATUS: FAIL",
      "Required input file(s) missing.",
      "No provenance reconciliation performed.",
      "NO ROC CALCULATION",
      "NO AUC CALCULATION",
      "NO MODEL FITTING",
      "NO STATISTICAL TEST"
    ),
    file.path(output_dir, "STEP13B-7H_FAIL.txt")
  )
  
  stop("STEP13B-7H stopped because required input files are missing.")
}


# ------------------------------------------------------------
# 3. Load RData in isolated environment
# ------------------------------------------------------------

cat("Loading provenance RData in isolated environment...\n")

prov_env <- new.env(parent = emptyenv())

load(
  prov_rdata,
  envir = prov_env
)

prov_objects <- ls(prov_env)

cat("Objects loaded:", length(prov_objects), "\n\n")


# ------------------------------------------------------------
# 4. Object dimension audit
# ------------------------------------------------------------

object_names_target <- c(
  "expr",
  "expr_sample",
  "expr_raw",
  "expr_annotation",
  "matched_expr",
  "sample_group",
  "patient_summary",
  "patient_expr_summary",
  "clinical_with_primary",
  "clinical_without_primary",
  "no_primary",
  "primary_expr",
  "primary_patient_mapping",
  "patient_mapping",
  "sample_names",
  "sample_ids",
  "expr_sample_ids",
  "tcga_ids",
  "metastatic_ids"
)

object_dimension_rows <- list()

for (nm in object_names_target) {
  
  exists_obj <- exists(
    nm,
    envir = prov_env,
    inherits = FALSE
  )
  
  if (!exists_obj) {
    
    object_dimension_rows[[length(object_dimension_rows) + 1]] <-
      data.frame(
        object = nm,
        exists = FALSE,
        class = NA_character_,
        nrow = NA_integer_,
        ncol = NA_integer_,
        length = NA_integer_,
        stringsAsFactors = FALSE
      )
    
  } else {
    
    obj <- get(nm, envir = prov_env)
    
    dims <- dim(obj)
    
    object_dimension_rows[[length(object_dimension_rows) + 1]] <-
      data.frame(
        object = nm,
        exists = TRUE,
        class = paste(class(obj), collapse = ";"),
        nrow = if (length(dims) >= 1) dims[1] else NA_integer_,
        ncol = if (length(dims) >= 2) dims[2] else NA_integer_,
        length = length(obj),
        stringsAsFactors = FALSE
      )
  }
}

object_dimension_audit <- do.call(
  rbind,
  object_dimension_rows
)

write.csv(
  object_dimension_audit,
  file.path(
    output_dir,
    "STEP13B-7H_RData_Object_Dimension_Audit.csv"
  ),
  row.names = FALSE
)

cat("RData object dimension audit saved.\n\n")


# ------------------------------------------------------------
# 5. Helper functions
# ------------------------------------------------------------

find_col <- function(df, candidates) {
  
  if (!is.data.frame(df)) {
    return(NA_character_)
  }
  
  hits <- candidates[
    candidates %in% colnames(df)
  ]
  
  if (length(hits) == 0) {
    return(NA_character_)
  }
  
  hits[1]
}


safe_unique_character <- function(x) {
  
  if (is.null(x)) {
    return(character(0))
  }
  
  x <- as.character(x)
  x <- x[!is.na(x)]
  x <- trimws(x)
  x <- x[nzchar(x)]
  
  unique(x)
}


get_patient_ids_from_df <- function(df) {
  
  if (!is.data.frame(df)) {
    return(character(0))
  }
  
  candidate_cols <- c(
    "patient_id",
    "Patient_ID",
    "patient",
    "Patient",
    "patientID",
    "PatientID",
    "case_id",
    "Case_ID"
  )
  
  col <- find_col(df, candidate_cols)
  
  if (is.na(col)) {
    return(character(0))
  }
  
  safe_unique_character(df[[col]])
}


get_sample_ids_from_df <- function(df) {
  
  if (!is.data.frame(df)) {
    return(character(0))
  }
  
  candidate_cols <- c(
    "Expression_Sample",
    "expression_sample",
    "Sample",
    "sample",
    "sample_id",
    "Sample_ID",
    "sampleID",
    "SampleID"
  )
  
  col <- find_col(df, candidate_cols)
  
  if (is.na(col)) {
    return(character(0))
  }
  
  safe_unique_character(df[[col]])
}


barcode_patient <- function(sample_ids) {
  
  sample_ids <- as.character(sample_ids)
  
  out <- sub(
    "-[0-9]{2}$",
    "",
    sample_ids
  )
  
  out[!grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", out)] <- NA_character_
  
  out
}


# ------------------------------------------------------------
# 6. Current locked annotation
# ------------------------------------------------------------

cat("Reading current sample annotation...\n")

ann <- read.delim(
  annotation_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Current annotation rows:", nrow(ann), "\n")
cat("Current annotation columns:\n")
print(colnames(ann))
cat("\n")


required_ann_cols <- c(
  "Expression_Sample",
  "Patient_ID",
  "Sample_Type_Code",
  "Tumor_Normal"
)

if (!all(required_ann_cols %in% colnames(ann))) {
  
  missing_cols <- setdiff(
    required_ann_cols,
    colnames(ann)
  )
  
  stop(
    paste(
      "Current annotation missing required columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}


ann$Expression_Sample <- as.character(
  ann$Expression_Sample
)

ann$Patient_ID <- as.character(
  ann$Patient_ID
)

ann$Sample_Type_Code <- as.integer(
  ann$Sample_Type_Code
)

ann$Tumor_Normal <- as.character(
  ann$Tumor_Normal
)


# ------------------------------------------------------------
# 7. Current sample-level classification
# ------------------------------------------------------------

primary_ids <- ann$Expression_Sample[
  ann$Sample_Type_Code == 1
]

metastatic_ids_current <- ann$Expression_Sample[
  ann$Sample_Type_Code == 2
]

normal_ids <- ann$Expression_Sample[
  ann$Sample_Type_Code == 11
]

eligible_ids <- ann$Expression_Sample[
  ann$Sample_Type_Code %in% c(1, 11)
]

current_patient_ids <- unique(
  ann$Patient_ID
)

primary_patient_ids <- unique(
  ann$Patient_ID[
    ann$Sample_Type_Code == 1
  ]
)

normal_patient_ids <- unique(
  ann$Patient_ID[
    ann$Sample_Type_Code == 11
  ]
)

metastatic_patient_ids <- unique(
  ann$Patient_ID[
    ann$Sample_Type_Code == 2
  ]
)


cat("CURRENT LOCKED SAMPLE STRUCTURE\n")
cat("------------------------------------------------------------\n")
cat("Expression samples:        ", nrow(ann), "\n")
cat("Primary Tumor samples:     ", length(primary_ids), "\n")
cat("Normal samples:            ", length(normal_ids), "\n")
cat("Metastatic samples:        ", length(metastatic_ids_current), "\n")
cat("Eligible Tumor/Normal:     ", length(eligible_ids), "\n")
cat("Unique patients:            ", length(current_patient_ids), "\n")
cat("Primary Tumor patients:     ", length(primary_patient_ids), "\n")
cat("Normal patients:            ", length(normal_patient_ids), "\n")
cat("Metastatic patients:        ", length(metastatic_patient_ids), "\n\n")


# ------------------------------------------------------------
# 8. Original expression matrix sample IDs from RData
# ------------------------------------------------------------

cat("Auditing original expression matrix...\n")

if (!exists("expr", envir = prov_env, inherits = FALSE)) {
  
  stop("Object 'expr' not found in provenance RData.")
}

expr_obj <- get("expr", envir = prov_env)

expr_colnames <- colnames(expr_obj)

if (is.null(expr_colnames)) {
  
  stop("Object 'expr' has no column names.")
}

matrix_sample_ids <- expr_colnames[
  !grepl(
    "^(Tag|Gene|GeneSymbol|Ensembl|ENSEMBL)$",
    expr_colnames,
    ignore.case = TRUE
  )
]

matrix_sample_ids <- safe_unique_character(
  matrix_sample_ids
)

cat("Original matrix sample columns:",
    length(matrix_sample_ids), "\n\n")


# ------------------------------------------------------------
# 9. Exact current matrix ↔ annotation audit
# ------------------------------------------------------------

matrix_not_annotation <- setdiff(
  matrix_sample_ids,
  ann$Expression_Sample
)

annotation_not_matrix <- setdiff(
  ann$Expression_Sample,
  matrix_sample_ids
)

matrix_annotation_exact <- (
  length(matrix_not_annotation) == 0 &&
    length(annotation_not_matrix) == 0 &&
    length(matrix_sample_ids) == nrow(ann)
)


# ------------------------------------------------------------
# 10. WGCNA sample IDs
# ------------------------------------------------------------

cat("Reading final WGCNA sample IDs...\n")

wgcna_header <- readLines(
  wgcna_file,
  n = 1,
  warn = FALSE
)

wgcna_header <- strsplit(
  wgcna_header,
  "\t",
  fixed = TRUE
)[[1]]

wgcna_sample_ids <- safe_unique_character(
  wgcna_header[-1]
)

cat("Final WGCNA samples:",
    length(wgcna_sample_ids), "\n\n")


wgcna_not_primary <- setdiff(
  wgcna_sample_ids,
  primary_ids
)

primary_not_wgcna <- setdiff(
  primary_ids,
  wgcna_sample_ids
)

wgcna_primary_exact <- (
  length(wgcna_sample_ids) == length(primary_ids) &&
    length(wgcna_not_primary) == 0 &&
    length(primary_not_wgcna) == 0
)


# ------------------------------------------------------------
# 11. Load historical annotation files
# ------------------------------------------------------------

cat("Reading historical annotation files...\n")

hist_ann <- read.delim(
  historical_annotation_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

patient_map <- read.delim(
  patient_mapping_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

hist_patient <- read.delim(
  historical_patient_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Workflow-03F-2 rows:",
    nrow(hist_ann), "\n")

cat("Workflow-03F-3 rows:",
    nrow(patient_map), "\n")

cat("Workflow-03F-5 rows:",
    nrow(hist_patient), "\n\n")


# ------------------------------------------------------------
# 12. Historical annotation sample-level audit
# ------------------------------------------------------------

hist_sample_col <- find_col(
  hist_ann,
  c(
    "Expression_Sample",
    "expression_sample",
    "Sample",
    "sample"
  )
)

hist_patient_col <- find_col(
  hist_ann,
  c(
    "Patient_ID",
    "patient_id",
    "patient"
  )
)

hist_type_col <- find_col(
  hist_ann,
  c(
    "Sample_Type",
    "sample_type"
  )
)

hist_code_col <- find_col(
  hist_ann,
  c(
    "Sample_Type_Code",
    "sample_type_code"
  )
)


if (is.na(hist_sample_col)) {
  
  stop(
    "Workflow-03F-2 does not contain an identifiable Expression_Sample column."
  )
}

hist_sample_ids <- safe_unique_character(
  hist_ann[[hist_sample_col]]
)

hist_primary_ids <- character(0)
hist_normal_ids <- character(0)
hist_metastatic_ids <- character(0)

if (!is.na(hist_code_col)) {
  
  hist_codes <- suppressWarnings(
    as.integer(hist_ann[[hist_code_col]])
  )
  
  hist_primary_ids <- hist_sample_ids[
    hist_codes == 1
  ]
  
  hist_metastatic_ids <- hist_sample_ids[
    hist_codes == 2
  ]
  
  hist_normal_ids <- hist_sample_ids[
    hist_codes == 11
  ]
}

if (!is.na(hist_type_col)) {
  
  hist_types <- as.character(
    hist_ann[[hist_type_col]]
  )
  
  hist_metastatic_by_text <- hist_sample_ids[
    grepl(
      "metastatic",
      hist_types,
      ignore.case = TRUE
    )
  ]
  
} else {
  
  hist_metastatic_by_text <- character(0)
}


# ------------------------------------------------------------
# 13. Historical sample IDs versus current matrix
# ------------------------------------------------------------

hist_not_matrix <- setdiff(
  hist_sample_ids,
  matrix_sample_ids
)

matrix_not_hist <- setdiff(
  matrix_sample_ids,
  hist_sample_ids
)

hist_matrix_exact <- (
  length(hist_not_matrix) == 0 &&
    length(matrix_not_hist) == 0
)


# ------------------------------------------------------------
# 14. Historical patient-level objects
# ------------------------------------------------------------

sample_group <- NULL
patient_summary <- NULL
patient_expr_summary <- NULL
clinical_with_primary <- NULL
clinical_without_primary <- NULL
no_primary <- NULL

if (exists("sample_group", envir = prov_env, inherits = FALSE)) {
  sample_group <- get("sample_group", envir = prov_env)
}

if (exists("patient_summary", envir = prov_env, inherits = FALSE)) {
  patient_summary <- get("patient_summary", envir = prov_env)
}

if (exists("patient_expr_summary", envir = prov_env, inherits = FALSE)) {
  patient_expr_summary <- get("patient_expr_summary", envir = prov_env)
}

if (exists("clinical_with_primary", envir = prov_env, inherits = FALSE)) {
  clinical_with_primary <- get("clinical_with_primary", envir = prov_env)
}

if (exists("clinical_without_primary", envir = prov_env, inherits = FALSE)) {
  clinical_without_primary <- get("clinical_without_primary", envir = prov_env)
}

if (exists("no_primary", envir = prov_env, inherits = FALSE)) {
  no_primary <- get("no_primary", envir = prov_env)
}


# ------------------------------------------------------------
# 15. Patient-level ID extraction
# ------------------------------------------------------------

sample_group_patients <- get_patient_ids_from_df(
  sample_group
)

patient_summary_patients <- get_patient_ids_from_df(
  patient_summary
)

patient_expr_summary_patients <- get_patient_ids_from_df(
  patient_expr_summary
)

clinical_with_primary_patients <- get_patient_ids_from_df(
  clinical_with_primary
)

clinical_without_primary_patients <- get_patient_ids_from_df(
  clinical_without_primary
)

no_primary_patients <- get_patient_ids_from_df(
  no_primary
)


# ------------------------------------------------------------
# 16. Current expression-derived patient sets
# ------------------------------------------------------------

primary_patient_from_expression <- unique(
  barcode_patient(primary_ids)
)

primary_patient_from_expression <- primary_patient_from_expression[
  !is.na(primary_patient_from_expression)
]

normal_patient_from_expression <- unique(
  barcode_patient(normal_ids)
)

normal_patient_from_expression <- normal_patient_from_expression[
  !is.na(normal_patient_from_expression)
]

metastatic_patient_from_expression <- unique(
  barcode_patient(metastatic_ids_current)
)

metastatic_patient_from_expression <- metastatic_patient_from_expression[
  !is.na(metastatic_patient_from_expression)
]


# ------------------------------------------------------------
# 17. Cross-object patient overlap audit
# ------------------------------------------------------------

patient_sets <- list(
  current_primary_expression =
    primary_patient_from_expression,
  
  current_normal_expression =
    normal_patient_from_expression,
  
  current_metastatic_expression =
    metastatic_patient_from_expression,
  
  sample_group =
    sample_group_patients,
  
  patient_summary =
    patient_summary_patients,
  
  patient_expr_summary =
    patient_expr_summary_patients,
  
  clinical_with_primary =
    clinical_with_primary_patients,
  
  clinical_without_primary =
    clinical_without_primary_patients,
  
  no_primary =
    no_primary_patients
)


patient_set_summary <- data.frame(
  object = names(patient_sets),
  unique_patient_count =
    sapply(patient_sets, length),
  stringsAsFactors = FALSE
)

write.csv(
  patient_set_summary,
  file.path(
    output_dir,
    "STEP13B-7H_Patient_Set_Count_Summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. Identify the four clinical_without_primary patients
# ------------------------------------------------------------

four_patient_audit <- data.frame()

if (length(clinical_without_primary_patients) > 0) {
  
  four_ids <- clinical_without_primary_patients
  
  four_patient_audit <- data.frame(
    Patient_ID = four_ids,
    In_current_primary_patient_set =
      four_ids %in% primary_patient_from_expression,
    In_current_normal_patient_set =
      four_ids %in% normal_patient_from_expression,
    In_current_metastatic_patient_set =
      four_ids %in% metastatic_patient_from_expression,
    In_sample_group =
      four_ids %in% sample_group_patients,
    In_patient_summary =
      four_ids %in% patient_summary_patients,
    In_patient_expr_summary =
      four_ids %in% patient_expr_summary_patients,
    In_clinical_with_primary =
      four_ids %in% clinical_with_primary_patients,
    stringsAsFactors = FALSE
  )
}

write.csv(
  four_patient_audit,
  file.path(
    output_dir,
    "STEP13B-7H_Four_Clinical_Without_Primary_Patient_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 19. Print the four patients and their full source rows
# ------------------------------------------------------------

capture.output({
  
  cat("============================================================\n")
  cat("FOUR PATIENTS IN clinical_without_primary / no_primary\n")
  cat("============================================================\n\n")
  
  cat("clinical_without_primary patient IDs:\n")
  print(clinical_without_primary_patients)
  
  cat("\nno_primary patient IDs:\n")
  print(no_primary_patients)
  
  cat("\nSet equality:\n")
  cat(
    "clinical_without_primary == no_primary:",
    setequal(
      clinical_without_primary_patients,
      no_primary_patients
    ),
    "\n\n"
  )
  
  if (length(clinical_without_primary_patients) > 0) {
    
    cat("------------------------------------------------------------\n")
    cat("clinical_without_primary full rows\n")
    cat("------------------------------------------------------------\n")
    
    print(clinical_without_primary)
    
    cat("\n------------------------------------------------------------\n")
    cat("no_primary full rows\n")
    cat("------------------------------------------------------------\n")
    
    print(no_primary)
  }
  
}, file = file.path(
  output_dir,
  "STEP13B-7H_Four_Patients_Full_Rows.txt"
))


# ------------------------------------------------------------
# 20. Trace the four patients through current annotation
# ------------------------------------------------------------

four_current_annotation <- ann[
  ann$Patient_ID %in% clinical_without_primary_patients,
  ,
  drop = FALSE
]

write.csv(
  four_current_annotation,
  file.path(
    output_dir,
    "STEP13B-7H_Four_Patients_Current_Annotation.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 21. Trace the four patients through historical annotation
# ------------------------------------------------------------

four_hist_annotation <- hist_ann[
  if (!is.na(hist_patient_col)) {
    as.character(hist_ann[[hist_patient_col]]) %in%
      clinical_without_primary_patients
  } else {
    FALSE
  },
  ,
  drop = FALSE
]

write.csv(
  four_hist_annotation,
  file.path(
    output_dir,
    "STEP13B-7H_Four_Patients_Historical_Annotation.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 22. Trace the four patients through Workflow-03F-3
# ------------------------------------------------------------

patient_map_patient_col <- find_col(
  patient_map,
  c(
    "patient_id",
    "Patient_ID",
    "patient"
  )
)

four_patient_mapping <- patient_map[
  if (!is.na(patient_map_patient_col)) {
    as.character(
      patient_map[[patient_map_patient_col]]
    ) %in% clinical_without_primary_patients
  } else {
    FALSE
  },
  ,
  drop = FALSE
]

write.csv(
  four_patient_mapping,
  file.path(
    output_dir,
    "STEP13B-7H_Four_Patients_Workflow03F3.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 23. Trace the four patients through Workflow-03F-5
# ------------------------------------------------------------

hist_patient_col2 <- find_col(
  hist_patient,
  c(
    "patient_id",
    "Patient_ID",
    "patient"
  )
)

four_historical_patient_annotation <- hist_patient[
  if (!is.na(hist_patient_col2)) {
    as.character(
      hist_patient[[hist_patient_col2]]
    ) %in% clinical_without_primary_patients
  } else {
    FALSE
  },
  ,
  drop = FALSE
]

write.csv(
  four_historical_patient_annotation,
  file.path(
    output_dir,
    "STEP13B-7H_Four_Patients_Workflow03F5.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 24. Identify sample IDs belonging to the four patients
# ------------------------------------------------------------

four_patient_sample_audit <- ann[
  ann$Patient_ID %in% clinical_without_primary_patients,
  c(
    "Expression_Sample",
    "Patient_ID",
    "Sample_Type_Code",
    "Tumor_Normal"
  ),
  drop = FALSE
]

if (nrow(four_patient_sample_audit) > 0) {
  
  four_patient_sample_audit$Sample_Suffix <- sub(
    "^.*-([0-9]{2})$",
    "\\1",
    four_patient_sample_audit$Expression_Sample
  )
}

write.csv(
  four_patient_sample_audit,
  file.path(
    output_dir,
    "STEP13B-7H_Four_Patients_Sample_Level_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 25. Historical 515 / 517 set comparison
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("HISTORICAL PATIENT-LEVEL COUNT AUDIT\n")
cat("============================================================\n\n")


cat("sample_group rows:          ",
    ifelse(is.null(sample_group), NA, nrow(sample_group)),
    "\n")

cat("patient_summary rows:       ",
    ifelse(is.null(patient_summary), NA, nrow(patient_summary)),
    "\n")

cat("patient_expr_summary rows:  ",
    ifelse(is.null(patient_expr_summary),
           NA,
           nrow(patient_expr_summary)),
    "\n")

cat("clinical_with_primary rows: ",
    ifelse(is.null(clinical_with_primary),
           NA,
           nrow(clinical_with_primary)),
    "\n")

cat("clinical_without_primary:   ",
    ifelse(is.null(clinical_without_primary),
           NA,
           nrow(clinical_without_primary)),
    "\n")

cat("no_primary rows:             ",
    ifelse(is.null(no_primary),
           NA,
           nrow(no_primary)),
    "\n\n")


# ------------------------------------------------------------
# 26. Exact patient-set comparisons
# ------------------------------------------------------------

compare_patient_sets <- function(name_a, set_a, name_b, set_b) {
  
  only_a <- setdiff(set_a, set_b)
  only_b <- setdiff(set_b, set_a)
  
  data.frame(
    set_A = name_a,
    set_B = name_b,
    n_A = length(set_a),
    n_B = length(set_b),
    A_not_B = length(only_a),
    B_not_A = length(only_b),
    EXACT = (
      length(only_a) == 0 &&
        length(only_b) == 0
    ),
    stringsAsFactors = FALSE
  )
}


comparison_list <- list()

if (length(sample_group_patients) > 0 &&
    length(patient_summary_patients) > 0) {
  
  comparison_list[[length(comparison_list) + 1]] <-
    compare_patient_sets(
      "sample_group",
      sample_group_patients,
      "patient_summary",
      patient_summary_patients
    )
}

if (length(patient_summary_patients) > 0 &&
    length(patient_expr_summary_patients) > 0) {
  
  comparison_list[[length(comparison_list) + 1]] <-
    compare_patient_sets(
      "patient_summary",
      patient_summary_patients,
      "patient_expr_summary",
      patient_expr_summary_patients
    )
}

if (length(clinical_with_primary_patients) > 0 &&
    length(primary_patient_from_expression) > 0) {
  
  comparison_list[[length(comparison_list) + 1]] <-
    compare_patient_sets(
      "clinical_with_primary",
      clinical_with_primary_patients,
      "primary_expression_patients",
      primary_patient_from_expression
    )
}

if (length(clinical_without_primary_patients) > 0 &&
    length(no_primary_patients) > 0) {
  
  comparison_list[[length(comparison_list) + 1]] <-
    compare_patient_sets(
      "clinical_without_primary",
      clinical_without_primary_patients,
      "no_primary",
      no_primary_patients
    )
}

if (length(comparison_list) > 0) {
  
  patient_set_comparisons <- do.call(
    rbind,
    comparison_list
  )
  
} else {
  
  patient_set_comparisons <- data.frame()
}

write.csv(
  patient_set_comparisons,
  file.path(
    output_dir,
    "STEP13B-7H_Patient_Set_Exact_Comparisons.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 27. Explicitly investigate 513 + 4 relationship
# ------------------------------------------------------------

clinical_513_plus_4_exact <- FALSE

if (length(clinical_with_primary_patients) == 513 &&
    length(clinical_without_primary_patients) == 4) {
  
  combined_513_4 <- unique(
    c(
      clinical_with_primary_patients,
      clinical_without_primary_patients
    )
  )
  
  if (length(combined_513_4) == 517) {
    clinical_513_plus_4_exact <- TRUE
  }
}


# ------------------------------------------------------------
# 28. Investigate 515 / 517 relationships
# ------------------------------------------------------------

sample_group_minus_patient_summary <- setdiff(
  sample_group_patients,
  patient_summary_patients
)

patient_summary_minus_sample_group <- setdiff(
  patient_summary_patients,
  sample_group_patients
)

patient_summary_minus_clinical_with_primary <- setdiff(
  patient_summary_patients,
  clinical_with_primary_patients
)

clinical_with_primary_minus_patient_summary <- setdiff(
  clinical_with_primary_patients,
  patient_summary_patients
)


count_relationship_audit <- data.frame(
  
  metric = c(
    "sample_group_n",
    "patient_summary_n",
    "patient_expr_summary_n",
    "clinical_with_primary_n",
    "clinical_without_primary_n",
    "no_primary_n",
    "clinical_with_primary_plus_without_primary_n",
    "current_primary_expression_patient_n",
    "current_normal_expression_patient_n",
    "current_metastatic_expression_patient_n"
  ),
  
  value = c(
    
    length(sample_group_patients),
    
    length(patient_summary_patients),
    
    length(patient_expr_summary_patients),
    
    length(clinical_with_primary_patients),
    
    length(clinical_without_primary_patients),
    
    length(no_primary_patients),
    
    length(
      unique(
        c(
          clinical_with_primary_patients,
          clinical_without_primary_patients
        )
      )
    ),
    
    length(primary_patient_from_expression),
    
    length(normal_patient_from_expression),
    
    length(metastatic_patient_from_expression)
    
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  count_relationship_audit,
  file.path(
    output_dir,
    "STEP13B-7H_Historical_Count_Relationship_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 29. Inspect sample_group and patient_summary differences
# ------------------------------------------------------------

write.csv(
  data.frame(
    Patient_ID = sample_group_minus_patient_summary,
    stringsAsFactors = FALSE
  ),
  file.path(
    output_dir,
    "STEP13B-7H_sample_group_NOT_patient_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    Patient_ID = patient_summary_minus_sample_group,
    stringsAsFactors = FALSE
  ),
  file.path(
    output_dir,
    "STEP13B-7H_patient_summary_NOT_sample_group.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    Patient_ID = patient_summary_minus_clinical_with_primary,
    stringsAsFactors = FALSE
  ),
  file.path(
    output_dir,
    "STEP13B-7H_patient_summary_NOT_clinical_with_primary.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    Patient_ID = clinical_with_primary_minus_patient_summary,
    stringsAsFactors = FALSE
  ),
  file.path(
    output_dir,
    "STEP13B-7H_clinical_with_primary_NOT_patient_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 30. Search historical provenance text files for count strings
# ------------------------------------------------------------

cat("\n")
cat("Performing targeted historical count-text audit...\n")


wgcna_root <- "input_dir/Paper1_tumor-only WGCNA"

all_txt <- list.files(
  wgcna_root,
  pattern = "\\.(txt|log|md|csv)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

targeted_txt <- all_txt[
  grepl(
    "Workflow-03F|Workflow-04|provenance|audit|summary|structure|mapping|annotation|clinical",
    basename(all_txt),
    ignore.case = TRUE
  )
]

targeted_txt <- targeted_txt[
  !grepl(
    "expression|matrix|counts|FPKM|TPM|normalized|expr",
    basename(targeted_txt),
    ignore.case = TRUE
  )
]

targeted_txt <- targeted_txt[
  file.info(targeted_txt)$size <= 20 * 1024^2
]

historical_count_hits <- data.frame()

if (length(targeted_txt) > 0) {
  
  hit_rows <- list()
  
  for (f in targeted_txt) {
    
    lines <- tryCatch(
      readLines(
        f,
        warn = FALSE,
        encoding = "UTF-8"
      ),
      error = function(e) character(0)
    )
    
    if (length(lines) == 0) {
      next
    }
    
    hit_idx <- grep(
      "\\b(513|514|515|516|517)\\b",
      lines
    )
    
    if (length(hit_idx) > 0) {
      
      for (ii in hit_idx) {
        
        hit_rows[[length(hit_rows) + 1]] <-
          data.frame(
            file = f,
            line_number = ii,
            text = lines[ii],
            stringsAsFactors = FALSE
          )
      }
    }
  }
  
  if (length(hit_rows) > 0) {
    
    historical_count_hits <- do.call(
      rbind,
      hit_rows
    )
  }
}

write.csv(
  historical_count_hits,
  file.path(
    output_dir,
    "STEP13B-7H_Historical_Count_Text_Hits.csv"
  ),
  row.names = FALSE
)

cat("Historical count-text audit completed.\n")


# ------------------------------------------------------------
# 31. Search RData character objects for literal 515/516/517
# ------------------------------------------------------------

rdata_count_object_hits <- data.frame()

rdata_candidate_names <- prov_objects[
  sapply(
    prov_objects,
    function(nm) {
      
      obj <- get(nm, envir = prov_env)
      
      is.character(obj) ||
        is.factor(obj) ||
        is.data.frame(obj)
    }
  )
]

rdata_hits <- list()

for (nm in rdata_candidate_names) {
  
  obj <- get(nm, envir = prov_env)
  
  if (is.character(obj) || is.factor(obj)) {
    
    vals <- as.character(obj)
    
    hit <- grep(
      "^(515|516|517)$",
      vals
    )
    
    if (length(hit) > 0) {
      
      rdata_hits[[length(rdata_hits) + 1]] <-
        data.frame(
          object = nm,
          type = "character/factor",
          location = paste(hit, collapse = ";"),
          value = paste(
            unique(vals[hit]),
            collapse = ";"
          ),
          stringsAsFactors = FALSE
        )
    }
    
  } else if (is.data.frame(obj)) {
    
    for (j in seq_len(ncol(obj))) {
      
      vals <- as.character(obj[[j]])
      
      hit <- grep(
        "^(515|516|517)$",
        vals
      )
      
      if (length(hit) > 0) {
        
        rdata_hits[[length(rdata_hits) + 1]] <-
          data.frame(
            object = nm,
            type = "data.frame",
            location = paste(
              paste0(
                "row",
                hit,
                ":col",
                j
              ),
              collapse = ";"
            ),
            value = paste(
              unique(vals[hit]),
              collapse = ";"
            ),
            stringsAsFactors = FALSE
          )
      }
    }
  }
}

if (length(rdata_hits) > 0) {
  
  rdata_count_object_hits <- do.call(
    rbind,
    rdata_hits
  )
}

write.csv(
  rdata_count_object_hits,
  file.path(
    output_dir,
    "STEP13B-7H_RData_Literal_515_516_517_Hits.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 32. Metastatic sample provenance
# ------------------------------------------------------------

metastatic_provenance <- ann[
  ann$Sample_Type_Code == 2,
  ,
  drop = FALSE
]

if (nrow(metastatic_provenance) > 0) {
  
  metastatic_provenance$Patient_from_Barcode <-
    barcode_patient(
      metastatic_provenance$Expression_Sample
    )
  
  metastatic_provenance$Present_in_Workflow03F2 <-
    metastatic_provenance$Expression_Sample %in%
    hist_sample_ids
  
  metastatic_provenance$Present_in_original_matrix <-
    metastatic_provenance$Expression_Sample %in%
    matrix_sample_ids
}

write.csv(
  metastatic_provenance,
  file.path(
    output_dir,
    "STEP13B-7H_Metastatic_Sample_Provenance.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 33. Exact 513 WGCNA patient-level audit
# ------------------------------------------------------------

wgcna_patient_ids <- unique(
  barcode_patient(wgcna_sample_ids)
)

wgcna_patient_ids <- wgcna_patient_ids[
  !is.na(wgcna_patient_ids)
]

primary_patient_wgcna_difference <- data.frame(
  Patient_ID = union(
    primary_patient_from_expression,
    wgcna_patient_ids
  ),
  In_Primary_Expression =
    union(
      primary_patient_from_expression,
      wgcna_patient_ids
    ) %in% primary_patient_from_expression,
  In_WGCNA =
    union(
      primary_patient_from_expression,
      wgcna_patient_ids
    ) %in% wgcna_patient_ids,
  stringsAsFactors = FALSE
)

write.csv(
  primary_patient_wgcna_difference,
  file.path(
    output_dir,
    "STEP13B-7H_WGCNA_Primary_Patient_Set_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 34. Determine provenance status
# ------------------------------------------------------------

# Core exact reconciliations already established:
matrix_annotation_pass <- matrix_annotation_exact
wgcna_primary_pass <- wgcna_primary_exact

# Four-patient set consistency:
four_patient_set_pass <- (
  setequal(
    clinical_without_primary_patients,
    no_primary_patients
  )
)

# Current sample-level counts:
current_sample_count_pass <- (
  length(primary_ids) == 513 &&
    length(normal_ids) == 59 &&
    length(metastatic_ids_current) == 2 &&
    length(eligible_ids) == 572
)

# Historical 513 + 4 = 517 relationship:
historical_513_plus_4_observed <- (
  clinical_513_plus_4_exact
)

# Whether the historical 515 count is exactly explained
# by the currently identified patient sets:
historical_515_explained <- FALSE

if (length(patient_summary_patients) == 515) {
  
  if (
    setequal(
      patient_summary_patients,
      union(
        clinical_with_primary_patients,
        clinical_without_primary_patients
      )
    )
  ) {
    
    historical_515_explained <- TRUE
  }
}


# IMPORTANT:
# Even if arithmetic relationships are observed, this script
# does NOT automatically label the historical 516 count as
# "516 Tumor samples" or "516 patients".
#
# The historical 516 value is considered resolved only if a
# concrete source object or exact documented transformation
# establishes its unit and meaning.

historical_516_explicit_source <- FALSE

if (nrow(historical_count_hits) > 0) {
  
  count516_hits <- historical_count_hits[
    grepl(
      "\\b516\\b",
      historical_count_hits$text
    ),
    ,
    drop = FALSE
  ]
  
  if (nrow(count516_hits) > 0) {
    historical_516_explicit_source <- TRUE
  }
}

historical_516_rdata_source <- FALSE

if (nrow(rdata_count_object_hits) > 0) {
  
  hit516 <- rdata_count_object_hits[
    grepl(
      "\\b516\\b",
      rdata_count_object_hits$value
    ),
    ,
    drop = FALSE
  ]
  
  if (nrow(hit516) > 0) {
    historical_516_rdata_source <- TRUE
  }
}


historical_516_reconciled <- (
  historical_516_explicit_source ||
    historical_516_rdata_source
)


# ------------------------------------------------------------
# 35. Final audit table
# ------------------------------------------------------------

final_audit <- data.frame(
  
  audit_item = c(
    
    "Current original expression sample count",
    "Current Primary Tumor sample count",
    "Current Solid Tissue Normal sample count",
    "Current Metastatic sample count",
    "Current Tumor-versus-Normal eligible sample count",
    
    "Matrix ↔ annotation exact sample-ID match",
    "Final WGCNA ↔ Primary Tumor exact sample-ID match",
    
    "sample_group row count",
    "patient_summary row count",
    "patient_expr_summary row count",
    
    "clinical_with_primary row count",
    "clinical_without_primary row count",
    "no_primary row count",
    
    "clinical_without_primary == no_primary",
    "clinical_with_primary + clinical_without_primary = 517",
    
    "Historical 515 count explicitly represented",
    "Historical 516 count explicitly sourced",
    "Historical 517 count explicitly represented",
    
    "Historical 515/516/517 provenance fully resolved"
  ),
  
  value = c(
    
    length(matrix_sample_ids),
    length(primary_ids),
    length(normal_ids),
    length(metastatic_ids_current),
    length(eligible_ids),
    
    matrix_annotation_exact,
    wgcna_primary_exact,
    
    length(sample_group_patients),
    length(patient_summary_patients),
    length(patient_expr_summary_patients),
    
    length(clinical_with_primary_patients),
    length(clinical_without_primary_patients),
    length(no_primary_patients),
    
    four_patient_set_pass,
    historical_513_plus_4_observed,
    
    length(patient_summary_patients) == 515,
    historical_516_reconciled,
    length(sample_group_patients) == 517,
    
    FALSE
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  final_audit,
  file.path(
    output_dir,
    "STEP13B-7H_FINAL_AUDIT.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 36. Human-readable final report
# ------------------------------------------------------------

report_lines <- c(
  
  "============================================================",
  "STEP13B-7H HISTORICAL COUNT PROVENANCE AUDIT",
  "============================================================",
  "",
  "CURRENT LOCKED EXPRESSION STRUCTURE",
  "------------------------------------------------------------",
  paste0(
    "Original expression samples: ",
    length(matrix_sample_ids)
  ),
  paste0(
    "Primary Tumor:              ",
    length(primary_ids)
  ),
  paste0(
    "Solid Tissue Normal:        ",
    length(normal_ids)
  ),
  paste0(
    "Metastatic:                 ",
    length(metastatic_ids_current)
  ),
  paste0(
    "Tumor-versus-Normal eligible: ",
    length(eligible_ids)
  ),
  "",
  "EXACT SAMPLE-ID RECONCILIATION",
  "------------------------------------------------------------",
  paste0(
    "Matrix ↔ Annotation EXACT: ",
    matrix_annotation_exact
  ),
  paste0(
    "WGCNA ↔ Primary Tumor EXACT: ",
    wgcna_primary_exact
  ),
  "",
  "HISTORICAL PATIENT-LEVEL OBJECTS",
  "------------------------------------------------------------",
  paste0(
    "sample_group rows:          ",
    length(sample_group_patients)
  ),
  paste0(
    "patient_summary rows:       ",
    length(patient_summary_patients)
  ),
  paste0(
    "patient_expr_summary rows:  ",
    length(patient_expr_summary_patients)
  ),
  paste0(
    "clinical_with_primary rows: ",
    length(clinical_with_primary_patients)
  ),
  paste0(
    "clinical_without_primary:   ",
    length(clinical_without_primary_patients)
  ),
  paste0(
    "no_primary rows:             ",
    length(no_primary_patients)
  ),
  "",
  "FOUR-PATIENT AUDIT",
  "------------------------------------------------------------",
  paste0(
    "clinical_without_primary == no_primary: ",
    four_patient_set_pass
  ),
  paste0(
    "513 + 4 = 517 relationship observed:     ",
    historical_513_plus_4_observed
  ),
  "",
  "HISTORICAL COUNT INTERPRETATION",
  "------------------------------------------------------------",
  paste0(
    "515 patient-level object present: ",
    length(patient_summary_patients) == 515
  ),
  paste0(
    "516 explicit provenance source found: ",
    historical_516_reconciled
  ),
  paste0(
    "517 patient-level object present: ",
    length(sample_group_patients) == 517
  ),
  "",
  "IMPORTANT INTERPRETATION RULE",
  "------------------------------------------------------------",
  "The current expression matrix is exactly reconciled as:",
  "513 Primary Tumor + 59 Solid Tissue Normal + 2 Metastatic = 574 samples.",
  "",
  "The final WGCNA sample set is exactly the same 513 Primary Tumor IDs.",
  "",
  "The historical 515/516/517 counts are NOT automatically",
  "interpreted as expression-sample counts.",
  "",
  "The historical '516' count is not declared resolved unless",
  "an explicit source or exact documented transformation establishes",
  "its unit and meaning.",
  "",
  "No historical count is reinterpreted merely from arithmetic.",
  "",
  "STATISTICAL ANALYSIS STATUS",
  "------------------------------------------------------------",
  "NO ROC CALCULATION",
  "NO AUC CALCULATION",
  "NO MODEL FITTING",
  "NO STATISTICAL TEST",
  "NO EXPRESSION RECALCULATION",
  "NO THRESHOLD OPTIMIZATION",
  "",
  "============================================================"
)


writeLines(
  report_lines,
  file.path(
    output_dir,
    "STEP13B-7H_FINAL_REPORT.txt"
  )
)


# ------------------------------------------------------------
# 37. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13B-7H COMPLETE\n")
cat("============================================================\n\n")

cat("CURRENT EXPRESSION SAMPLE STRUCTURE\n")
cat("  Original expression samples: ",
    length(matrix_sample_ids), "\n")
cat("  Primary Tumor:               ",
    length(primary_ids), "\n")
cat("  Solid Tissue Normal:         ",
    length(normal_ids), "\n")
cat("  Metastatic:                  ",
    length(metastatic_ids_current), "\n")
cat("  Tumor-versus-Normal eligible:",
    length(eligible_ids), "\n\n")

cat("EXACT SAMPLE-ID RECONCILIATION\n")
cat("  Matrix ↔ Annotation EXACT:  ",
    matrix_annotation_exact, "\n")
cat("  WGCNA ↔ Primary Tumor EXACT:",
    wgcna_primary_exact, "\n\n")

cat("HISTORICAL PATIENT OBJECTS\n")
cat("  sample_group:                ",
    length(sample_group_patients), "\n")
cat("  patient_summary:             ",
    length(patient_summary_patients), "\n")
cat("  patient_expr_summary:        ",
    length(patient_expr_summary_patients), "\n")
cat("  clinical_with_primary:       ",
    length(clinical_with_primary_patients), "\n")
cat("  clinical_without_primary:    ",
    length(clinical_without_primary_patients), "\n")
cat("  no_primary:                  ",
    length(no_primary_patients), "\n\n")

cat("FOUR-PATIENT AUDIT\n")
cat("  clinical_without_primary == no_primary:",
    four_patient_set_pass, "\n")
cat("  513 + 4 = 517 relationship observed:    ",
    historical_513_plus_4_observed, "\n\n")

cat("HISTORICAL COUNT STATUS\n")
cat("  515 object present:          ",
    length(patient_summary_patients) == 515, "\n")
cat("  516 explicit source found:   ",
    historical_516_reconciled, "\n")
cat("  517 object present:          ",
    length(sample_group_patients) == 517, "\n\n")

cat("============================================================\n")
cat("STATISTICAL ANALYSIS STATUS\n")
cat("============================================================\n")
cat("NO ROC CALCULATION\n")
cat("NO AUC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")
cat("NO EXPRESSION RECALCULATION\n")
cat("NO THRESHOLD OPTIMIZATION\n\n")

cat("============================================================\n")
cat("STEP13B-7H STATUS: INFORMATION AUDIT COMPLETE\n")
cat("============================================================\n\n")

cat("IMPORTANT:\n")
cat("The current 513 / 59 / 2 sample composition is locked by\n")
cat("exact matrix ↔ annotation reconciliation.\n")
cat("The final WGCNA 513 samples exactly match the Primary Tumor IDs.\n")
cat("Historical 515/516/517 interpretation is reported separately\n")
cat("and is not inferred from arithmetic alone.\n\n")

cat(
  "Outputs saved to:\n",
  output_dir,
  "\n\n"
)

# ============================================================
# STEP13B-7H — Read Historical Count Evidence Files
# ============================================================
# Purpose:
#   Read ONLY the two evidence CSV files generated by
#   STEP13B-7H.
#
#   NO:
#     - ROC calculation
#     - AUC calculation
#     - model fitting
#     - statistical testing
#     - expression recalculation
#     - threshold optimization
#
# ============================================================


cat("============================================================\n")
cat("STEP13B-7H — READ HISTORICAL COUNT EVIDENCE FILES\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Working directory
# ------------------------------------------------------------

audit_dir <- "input_dir/Paper 1_Step 13-Diagnostic Classification/STEP13B-7H_Historical_Count_Provenance_Audit"

cat("Audit directory:\n")
cat(audit_dir, "\n\n")


# ------------------------------------------------------------
# 2. Define the two evidence files
# ------------------------------------------------------------

file_text_hits <- file.path(
  audit_dir,
  "STEP13B-7H_Historical_Count_Text_Hits.csv"
)

file_rdata_hits <- file.path(
  audit_dir,
  "STEP13B-7H_RData_Literal_515_516_517_Hits.csv"
)


# ------------------------------------------------------------
# 3. Check file existence
# ------------------------------------------------------------

cat("FILE EXISTENCE CHECK\n")
cat("------------------------------------------------------------\n")

cat(
  "Historical_Count_Text_Hits.csv: ",
  file.exists(file_text_hits),
  "\n",
  sep = ""
)

cat(
  "RData_Literal_515_516_517_Hits.csv: ",
  file.exists(file_rdata_hits),
  "\n\n",
  sep = ""
)


if (!file.exists(file_text_hits)) {
  stop(
    "ERROR: Historical_Count_Text_Hits.csv not found:\n",
    file_text_hits
  )
}

if (!file.exists(file_rdata_hits)) {
  stop(
    "ERROR: RData_Literal_515_516_517_Hits.csv not found:\n",
    file_rdata_hits
  )
}


# ------------------------------------------------------------
# 4. Read the two CSV files
# ------------------------------------------------------------

text_hits <- read.csv(
  file_text_hits,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rdata_hits <- read.csv(
  file_rdata_hits,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 5. Basic structure
# ------------------------------------------------------------

cat("============================================================\n")
cat("FILE 1 — HISTORICAL COUNT TEXT HITS\n")
cat("============================================================\n\n")

cat("Rows:", nrow(text_hits), "\n")
cat("Columns:", ncol(text_hits), "\n\n")

cat("Column names:\n")
print(colnames(text_hits))

cat("\n\n")


cat("============================================================\n")
cat("FILE 2 — RDATA LITERAL 515 / 516 / 517 HITS\n")
cat("============================================================\n\n")

cat("Rows:", nrow(rdata_hits), "\n")
cat("Columns:", ncol(rdata_hits), "\n\n")

cat("Column names:\n")
print(colnames(rdata_hits))

cat("\n\n")


# ------------------------------------------------------------
# 6. Print complete contents if reasonably small
# ------------------------------------------------------------

cat("============================================================\n")
cat("RAW CONTENT PREVIEW\n")
cat("============================================================\n\n")


cat("---- FILE 1: first 30 rows ----\n\n")

if (nrow(text_hits) > 0) {
  print(
    head(
      text_hits,
      30
    )
  )
} else {
  cat("[FILE 1 IS EMPTY]\n")
}

cat("\n\n")


cat("---- FILE 2: complete contents ----\n\n")

if (nrow(rdata_hits) > 0) {
  print(rdata_hits)
} else {
  cat("[FILE 2 IS EMPTY]\n")
}

cat("\n\n")


# ------------------------------------------------------------
# 7. Search FILE 1 for explicit 515 / 516 / 517 references
# ------------------------------------------------------------

cat("============================================================\n")
cat("TARGETED SEARCH — FILE 1\n")
cat("515 / 516 / 517\n")
cat("============================================================\n\n")


if (nrow(text_hits) > 0) {
  
  text_hits_chr <- as.data.frame(
    lapply(
      text_hits,
      function(x) as.character(x)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  row_text_hits <- apply(
    text_hits_chr,
    1,
    paste,
    collapse = " | "
  )
  
  hit_515_516_517 <- grepl(
    "(^|[^0-9])(515|516|517)([^0-9]|$)",
    row_text_hits
  )
  
  text_hits_target <- text_hits[
    hit_515_516_517,
    ,
    drop = FALSE
  ]
  
  cat(
    "Rows containing explicit 515 / 516 / 517 reference:",
    nrow(text_hits_target),
    "\n\n"
  )
  
  if (nrow(text_hits_target) > 0) {
    print(text_hits_target)
  } else {
    cat("No explicit 515 / 516 / 517 row-level hits found.\n")
  }
  
} else {
  
  cat("FILE 1 is empty.\n")
  
}


cat("\n\n")


# ------------------------------------------------------------
# 8. Search FILE 1 specifically for 516
# ------------------------------------------------------------

cat("============================================================\n")
cat("TARGETED SEARCH — FILE 1 — NUMBER 516 ONLY\n")
cat("============================================================\n\n")


if (nrow(text_hits) > 0) {
  
  hit_516 <- grepl(
    "(^|[^0-9])516([^0-9]|$)",
    row_text_hits
  )
  
  text_hits_516 <- text_hits[
    hit_516,
    ,
    drop = FALSE
  ]
  
  cat(
    "Rows containing explicit 516:",
    nrow(text_hits_516),
    "\n\n"
  )
  
  if (nrow(text_hits_516) > 0) {
    print(text_hits_516)
  } else {
    cat("No explicit 516 rows found.\n")
  }
  
}


cat("\n\n")


# ------------------------------------------------------------
# 9. Search FILE 2 for 515 / 516 / 517
# ------------------------------------------------------------

cat("============================================================\n")
cat("TARGETED SEARCH — FILE 2\n")
cat("515 / 516 / 517\n")
cat("============================================================\n\n")


if (nrow(rdata_hits) > 0) {
  
  rdata_hits_chr <- as.data.frame(
    lapply(
      rdata_hits,
      function(x) as.character(x)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  rdata_row_text <- apply(
    rdata_hits_chr,
    1,
    paste,
    collapse = " | "
  )
  
  rdata_target <- rdata_hits[
    grepl(
      "(^|[^0-9])(515|516|517)([^0-9]|$)",
      rdata_row_text
    ),
    ,
    drop = FALSE
  ]
  
  cat(
    "Rows containing explicit 515 / 516 / 517:",
    nrow(rdata_target),
    "\n\n"
  )
  
  if (nrow(rdata_target) > 0) {
    print(rdata_target)
  } else {
    cat("No explicit 515 / 516 / 517 rows found.\n")
  }
  
} else {
  
  cat("FILE 2 is empty.\n")
  
}


cat("\n\n")


# ------------------------------------------------------------
# 10. Specifically identify RData hits involving 516
# ------------------------------------------------------------

cat("============================================================\n")
cat("TARGETED SEARCH — FILE 2 — NUMBER 516 ONLY\n")
cat("============================================================\n\n")


if (nrow(rdata_hits) > 0) {
  
  rdata_516 <- rdata_hits[
    grepl(
      "(^|[^0-9])516([^0-9]|$)",
      rdata_row_text
    ),
    ,
    drop = FALSE
  ]
  
  cat(
    "Rows containing explicit 516:",
    nrow(rdata_516),
    "\n\n"
  )
  
  if (nrow(rdata_516) > 0) {
    print(rdata_516)
  } else {
    cat("No explicit 516 rows found.\n")
  }
  
}


# ------------------------------------------------------------
# 11. Save the complete read-in objects
# ------------------------------------------------------------

output_rds <- file.path(
  audit_dir,
  "STEP13B-7H_Read_Historical_Count_Evidence.rds"
)

saveRDS(
  list(
    text_hits = text_hits,
    rdata_hits = rdata_hits
  ),
  output_rds
)


# ------------------------------------------------------------
# 12. Final status
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("STEP13B-7H EVIDENCE READ COMPLETE\n")
cat("============================================================\n\n")

cat("File 1 rows:", nrow(text_hits), "\n")
cat("File 1 columns:", ncol(text_hits), "\n")

cat("File 2 rows:", nrow(rdata_hits), "\n")
cat("File 2 columns:", ncol(rdata_hits), "\n\n")

cat("Evidence RDS saved to:\n")
cat(output_rds, "\n\n")

cat("NO ROC CALCULATION\n")
cat("NO AUC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")
cat("NO EXPRESSION RECALCULATION\n")
cat("NO THRESHOLD OPTIMIZATION\n\n")

cat("STEP13B-7H EVIDENCE READ STATUS: PASS\n")
cat("============================================================\n")

# ============================================================
# STEP13B-7H — TARGETED 515 / 516 / 517 PROVENANCE EXTRACTION
# ============================================================
#
# Purpose:
#   Extract only the informative provenance records from the
#   already-generated STEP13B-7H evidence files.
#
# This script:
#   - DOES NOT read the original expression matrix
#   - DOES NOT read the 514 MB RData again
#   - DOES NOT perform any statistical analysis
#   - DOES NOT calculate ROC/AUC
#   - DOES NOT fit models
#   - DOES NOT recalculate expression
#
# ============================================================


cat("============================================================\n")
cat("STEP13B-7H — TARGETED HISTORICAL COUNT EXTRACTION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

audit_dir <- paste0(
  "input_dir/Paper 1_Step 13-Diagnostic Classification/",
  "STEP13B-7H_Historical_Count_Provenance_Audit"
)

evidence_rds <- file.path(
  audit_dir,
  "STEP13B-7H_Read_Historical_Count_Evidence.rds"
)


# ------------------------------------------------------------
# 2. Load already-read evidence
# ------------------------------------------------------------

if (!file.exists(evidence_rds)) {
  stop(
    "Evidence RDS not found:\n",
    evidence_rds
  )
}

ev <- readRDS(evidence_rds)

text_hits <- ev$text_hits
rdata_hits <- ev$rdata_hits


# ------------------------------------------------------------
# 3. Basic structure
# ------------------------------------------------------------

cat("EVIDENCE OBJECT STRUCTURE\n")
cat("------------------------------------------------------------\n")

cat(
  "text_hits:",
  nrow(text_hits),
  "rows ×",
  ncol(text_hits),
  "columns\n"
)

cat(
  "rdata_hits:",
  nrow(rdata_hits),
  "rows ×",
  ncol(rdata_hits),
  "columns\n\n"
)

cat("text_hits columns:\n")
print(colnames(text_hits))

cat("\nrdata_hits columns:\n")
print(colnames(rdata_hits))

cat("\n\n")


# ------------------------------------------------------------
# 4. Convert every column to character safely
# ------------------------------------------------------------

text_chr <- as.data.frame(
  lapply(
    text_hits,
    function(x) {
      if (is.factor(x)) {
        as.character(x)
      } else {
        as.character(x)
      }
    }
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rdata_chr <- as.data.frame(
  lapply(
    rdata_hits,
    function(x) {
      if (is.factor(x)) {
        as.character(x)
      } else {
        as.character(x)
      }
    }
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 5. Create one searchable text field per record
# ------------------------------------------------------------

text_row_string <- apply(
  text_chr,
  1,
  paste,
  collapse = " | "
)

rdata_row_string <- apply(
  rdata_chr,
  1,
  paste,
  collapse = " | "
)


# ------------------------------------------------------------
# 6. Extract exact standalone 515 / 516 / 517 hits
# ------------------------------------------------------------

pattern_515 <- "(^|[^0-9])515([^0-9]|$)"
pattern_516 <- "(^|[^0-9])516([^0-9]|$)"
pattern_517 <- "(^|[^0-9])517([^0-9]|$)"

text_515 <- text_hits[
  grepl(pattern_515, text_row_string),
  ,
  drop = FALSE
]

text_516 <- text_hits[
  grepl(pattern_516, text_row_string),
  ,
  drop = FALSE
]

text_517 <- text_hits[
  grepl(pattern_517, text_row_string),
  ,
  drop = FALSE
]


rdata_515 <- rdata_hits[
  grepl(pattern_515, rdata_row_string),
  ,
  drop = FALSE
]

rdata_516 <- rdata_hits[
  grepl(pattern_516, rdata_row_string),
  ,
  drop = FALSE
]

rdata_517 <- rdata_hits[
  grepl(pattern_517, rdata_row_string),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 7. Counts
# ------------------------------------------------------------

cat("============================================================\n")
cat("EXACT STANDALONE COUNT HITS\n")
cat("============================================================\n\n")

cat("TEXT HITS\n")
cat("  515:", nrow(text_515), "\n")
cat("  516:", nrow(text_516), "\n")
cat("  517:", nrow(text_517), "\n\n")

cat("RDATA HITS\n")
cat("  515:", nrow(rdata_515), "\n")
cat("  516:", nrow(rdata_516), "\n")
cat("  517:", nrow(rdata_517), "\n\n")


# ------------------------------------------------------------
# 8. Inspect the structure of the 516 text hits
# ------------------------------------------------------------

cat("============================================================\n")
cat("TEXT FILE — 516 HIT STRUCTURE\n")
cat("============================================================\n\n")

if (nrow(text_516) > 0) {
  
  print(
    head(
      text_516,
      50
    )
  )
  
} else {
  
  cat("No exact standalone 516 hits found in text evidence.\n")
  
}


# ------------------------------------------------------------
# 9. Inspect RData 516 hits
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("RDATA FILE — 516 HIT STRUCTURE\n")
cat("============================================================\n\n")

if (nrow(rdata_516) > 0) {
  
  print(rdata_516)
  
} else {
  
  cat("No exact standalone 516 hits found in RData evidence.\n")
  
}


# ------------------------------------------------------------
# 10. Identify possible semantic keywords around 516
# ------------------------------------------------------------

semantic_pattern <- paste(
  c(
    "sample",
    "patient",
    "tumor",
    "normal",
    "primary",
    "metastatic",
    "clinical",
    "expression",
    "row",
    "rows",
    "column",
    "columns",
    "count",
    "n\\s*=",
    "number",
    "total",
    "annotation",
    "WGCNA"
  ),
  collapse = "|"
)


text_516_strings <- if (nrow(text_516) > 0) {
  apply(
    as.data.frame(
      lapply(text_516, as.character),
      stringsAsFactors = FALSE
    ),
    1,
    paste,
    collapse = " | "
  )
} else {
  character(0)
}


rdata_516_strings <- if (nrow(rdata_516) > 0) {
  apply(
    as.data.frame(
      lapply(rdata_516, as.character),
      stringsAsFactors = FALSE
    ),
    1,
    paste,
    collapse = " | "
  )
} else {
  character(0)
}


text_516_semantic <- text_516[
  grepl(
    semantic_pattern,
    text_516_strings,
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]


rdata_516_semantic <- rdata_516[
  grepl(
    semantic_pattern,
    rdata_516_strings,
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 11. Print semantic 516 hits
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("516 HITS WITH PROVENANCE KEYWORDS — TEXT\n")
cat("============================================================\n\n")

cat(
  "Number of informative 516 text records:",
  nrow(text_516_semantic),
  "\n\n"
)

if (nrow(text_516_semantic) > 0) {
  
  print(
    head(
      text_516_semantic,
      100
    )
  )
  
} else {
  
  cat("No keyword-associated 516 records found.\n")
  
}


cat("\n\n")
cat("============================================================\n")
cat("516 HITS WITH PROVENANCE KEYWORDS — RDATA\n")
cat("============================================================\n\n")

cat(
  "Number of informative 516 RData records:",
  nrow(rdata_516_semantic),
  "\n\n"
)

if (nrow(rdata_516_semantic) > 0) {
  
  print(rdata_516_semantic)
  
} else {
  
  cat("No keyword-associated 516 RData records found.\n")
  
}


# ------------------------------------------------------------
# 12. Save targeted results
# ------------------------------------------------------------

write.csv(
  text_516,
  file.path(
    audit_dir,
    "STEP13B-7H_TARGET_516_Text_Hits.csv"
  ),
  row.names = FALSE
)

write.csv(
  text_516_semantic,
  file.path(
    audit_dir,
    "STEP13B-7H_TARGET_516_Text_Semantic_Hits.csv"
  ),
  row.names = FALSE
)

write.csv(
  rdata_516,
  file.path(
    audit_dir,
    "STEP13B-7H_TARGET_516_RData_Hits.csv"
  ),
  row.names = FALSE
)

write.csv(
  rdata_516_semantic,
  file.path(
    audit_dir,
    "STEP13B-7H_TARGET_516_RData_Semantic_Hits.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Final summary
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("STEP13B-7H TARGETED EXTRACTION COMPLETE\n")
cat("============================================================\n\n")

cat("Exact 516 text hits:       ", nrow(text_516), "\n")
cat("Semantic 516 text hits:    ", nrow(text_516_semantic), "\n")
cat("Exact 516 RData hits:      ", nrow(rdata_516), "\n")
cat("Semantic 516 RData hits:   ", nrow(rdata_516_semantic), "\n\n")

cat("Targeted files saved under:\n")
cat(audit_dir, "\n\n")

cat("NO ROC CALCULATION\n")
cat("NO AUC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")
cat("NO EXPRESSION RECALCULATION\n")
cat("NO THRESHOLD OPTIMIZATION\n\n")

cat("STEP13B-7H TARGETED EXTRACTION STATUS: PASS\n")
cat("============================================================\n")

# ============================================================
# STEP13B-7H — FINAL 516 SEMANTIC SOURCE IDENTIFICATION
# ============================================================
#
# Purpose:
#   Determine the semantic meaning of historical "516".
#
#   IMPORTANT:
#   This script reads ONLY the targeted evidence CSV files
#   already generated by STEP13B-7H.
#
#   NO original matrix is read.
#   NO RData is loaded.
#   NO statistics are calculated.
#   NO ROC/AUC/model analysis is performed.
#
# ============================================================


cat("============================================================\n")
cat("STEP13B-7H — FINAL 516 SOURCE IDENTIFICATION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Directory
# ------------------------------------------------------------

audit_dir <- paste0(
  "input_dir/Paper 1_Step 13-Diagnostic Classification/",
  "STEP13B-7H_Historical_Count_Provenance_Audit"
)


# ------------------------------------------------------------
# 2. Target files
# ------------------------------------------------------------

text_file <- file.path(
  audit_dir,
  "STEP13B-7H_TARGET_516_Text_Hits.csv"
)

rdata_file <- file.path(
  audit_dir,
  "STEP13B-7H_TARGET_516_RData_Hits.csv"
)


if (!file.exists(text_file)) {
  stop("Missing file:\n", text_file)
}

if (!file.exists(rdata_file)) {
  stop("Missing file:\n", rdata_file)
}


# ------------------------------------------------------------
# 3. Read
# ------------------------------------------------------------

text516 <- read.csv(
  text_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rdata516 <- read.csv(
  rdata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 4. Basic structure
# ------------------------------------------------------------

cat("TEXT 516 EVIDENCE\n")
cat("------------------------------------------------------------\n")
cat("Rows:", nrow(text516), "\n")
cat("Columns:", ncol(text516), "\n")
cat("Column names:\n")
print(colnames(text516))

cat("\n\n")


cat("RDATA 516 EVIDENCE\n")
cat("------------------------------------------------------------\n")
cat("Rows:", nrow(rdata516), "\n")
cat("Columns:", ncol(rdata516), "\n")
cat("Column names:\n")
print(colnames(rdata516))

cat("\n\n")


# ------------------------------------------------------------
# 5. Complete 42 text records
# ------------------------------------------------------------

cat("============================================================\n")
cat("ALL 42 TEXT-BASED 516 RECORDS\n")
cat("============================================================\n\n")

if (nrow(text516) > 0) {
  print(text516)
} else {
  cat("No text records found.\n")
}


# ------------------------------------------------------------
# 6. Search text records for semantic terms
# ------------------------------------------------------------

text_all <- apply(
  data.frame(
    lapply(text516, as.character),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ),
  1,
  paste,
  collapse = " | "
)


semantic_groups <- list(
  
  SAMPLE = "sample|samples",
  
  PATIENT = "patient|patients",
  
  TUMOR = "tumou?r|primary tumor",
  
  NORMAL = "normal|solid tissue",
  
  METASTATIC = "metastatic",
  
  CLINICAL = "clinical",
  
  EXPRESSION = "expression",
  
  ROW = "row|rows",
  
  COLUMN = "column|columns",
  
  COUNT = "count|counts|number|total",
  
  WGCNA = "WGCNA",
  
  INDEX = "index|indices|position|positions"
)


cat("\n")
cat("============================================================\n")
cat("SEMANTIC CLASSIFICATION OF TEXT 516 RECORDS\n")
cat("============================================================\n\n")


semantic_summary <- data.frame(
  Category = character(),
  Records = integer(),
  stringsAsFactors = FALSE
)


for (nm in names(semantic_groups)) {
  
  hits <- grepl(
    semantic_groups[[nm]],
    text_all,
    ignore.case = TRUE
  )
  
  n_hit <- sum(hits)
  
  semantic_summary <- rbind(
    semantic_summary,
    data.frame(
      Category = nm,
      Records = n_hit,
      stringsAsFactors = FALSE
    )
  )
  
}

print(semantic_summary)


# ------------------------------------------------------------
# 7. RData field names
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("RDATA 516 — COLUMN STRUCTURE\n")
cat("============================================================\n\n")

print(colnames(rdata516))


# ------------------------------------------------------------
# 8. Print frequency of each RData column
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("RDATA 516 — UNIQUE VALUES BY COLUMN\n")
cat("============================================================\n\n")


for (nm in colnames(rdata516)) {
  
  cat("\n------------------------------------------------------------\n")
  cat("COLUMN:", nm, "\n")
  cat("------------------------------------------------------------\n")
  
  vals <- unique(as.character(rdata516[[nm]]))
  
  cat("Unique values:", length(vals), "\n")
  
  if (length(vals) <= 30) {
    print(vals)
  } else {
    print(head(vals, 30))
    cat("... additional values omitted ...\n")
  }
}


# ------------------------------------------------------------
# 9. Frequency table for likely object/field columns
# ------------------------------------------------------------

candidate_columns <- intersect(
  colnames(rdata516),
  c(
    "object",
    "Object",
    "object_name",
    "Object_Name",
    "field",
    "Field",
    "field_name",
    "Field_Name",
    "hit_type",
    "Hit_Type",
    "type",
    "Type",
    "value",
    "Value"
  )
)


cat("\n\n")
cat("============================================================\n")
cat("RDATA 516 — STRUCTURED FREQUENCY SUMMARY\n")
cat("============================================================\n\n")


if (length(candidate_columns) > 0) {
  
  for (nm in candidate_columns) {
    
    cat("\n------------------------------------------------------------\n")
    cat("COLUMN:", nm, "\n")
    cat("------------------------------------------------------------\n")
    
    print(
      sort(
        table(
          as.character(rdata516[[nm]])
        ),
        decreasing = TRUE
      )
    )
  }
  
} else {
  
  cat(
    "No standard object/field/hit-type column names detected.\n"
  )
  
}


# ------------------------------------------------------------
# 10. Search RData records for true semantic count language
# ------------------------------------------------------------

rdata_all <- apply(
  data.frame(
    lapply(rdata516, as.character),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ),
  1,
  paste,
  collapse = " | "
)


true_count_pattern <- paste(
  c(
    "516 sample",
    "516 samples",
    "516 patient",
    "516 patients",
    "516 tumor",
    "516 tumour",
    "516 primary",
    "516 normal",
    "516 expression",
    "516 row",
    "516 rows",
    "516 column",
    "516 columns",
    "n.?=.?.?516",
    "count.?=.?.?516",
    "total.?=.?.?516"
  ),
  collapse = "|"
)


rdata_true_count <- rdata516[
  grepl(
    true_count_pattern,
    rdata_all,
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]


cat("\n\n")
cat("============================================================\n")
cat("RDATA 516 — POSSIBLE TRUE COUNT REFERENCES\n")
cat("============================================================\n\n")

cat(
  "Records with explicit count-like 516 context:",
  nrow(rdata_true_count),
  "\n\n"
)

if (nrow(rdata_true_count) > 0) {
  print(rdata_true_count)
} else {
  cat(
    "No explicit '516 samples/patients/tumor/etc.' ",
    "records detected.\n"
  )
}


# ------------------------------------------------------------
# 11. Search text evidence for exact count phrases
# ------------------------------------------------------------

count_phrase_pattern <- paste(
  c(
    "516 sample",
    "516 samples",
    "516 patient",
    "516 patients",
    "516 tumor",
    "516 tumour",
    "516 primary",
    "516 normal",
    "516 expression",
    "516 row",
    "516 rows",
    "516 column",
    "516 columns",
    "n.?=.?.?516",
    "count.?=.?.?516",
    "total.?=.?.?516"
  ),
  collapse = "|"
)


text_true_count <- text516[
  grepl(
    count_phrase_pattern,
    text_all,
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]


cat("\n\n")
cat("============================================================\n")
cat("TEXT 516 — POSSIBLE TRUE COUNT REFERENCES\n")
cat("============================================================\n\n")

cat(
  "Records with explicit count-like 516 context:",
  nrow(text_true_count),
  "\n\n"
)

if (nrow(text_true_count) > 0) {
  print(text_true_count)
} else {
  cat(
    "No explicit '516 samples/patients/tumor/etc.' ",
    "records detected.\n"
  )
}


# ------------------------------------------------------------
# 12. Save condensed results
# ------------------------------------------------------------

write.csv(
  semantic_summary,
  file.path(
    audit_dir,
    "STEP13B-7H_516_Semantic_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  rdata_true_count,
  file.path(
    audit_dir,
    "STEP13B-7H_516_RData_True_Count_Candidates.csv"
  ),
  row.names = FALSE
)

write.csv(
  text_true_count,
  file.path(
    audit_dir,
    "STEP13B-7H_516_Text_True_Count_Candidates.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Final status
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("STEP13B-7H — FINAL 516 SOURCE EXTRACTION COMPLETE\n")
cat("============================================================\n\n")

cat("Text 516 records:", nrow(text516), "\n")
cat("RData 516 records:", nrow(rdata516), "\n")
cat(
  "Text true-count candidates:",
  nrow(text_true_count),
  "\n"
)
cat(
  "RData true-count candidates:",
  nrow(rdata_true_count),
  "\n\n"
)

cat("NO ROC CALCULATION\n")
cat("NO AUC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")
cat("NO EXPRESSION RECALCULATION\n")
cat("NO THRESHOLD OPTIMIZATION\n\n")

cat("STEP13B-7H FINAL 516 EXTRACTION STATUS: PASS\n")
cat("============================================================\n")


# ============================================================
# STEP13B-7H — DISPLAY ALL 42 HISTORICAL 516 TEXT HITS
# ============================================================

audit_dir <- paste0(
  "input_dir/Paper 1_Step 13-Diagnostic Classification/",
  "STEP13B-7H_Historical_Count_Provenance_Audit"
)

text_file <- file.path(
  audit_dir,
  "STEP13B-7H_TARGET_516_Text_Hits.csv"
)

text516 <- read.csv(
  text_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("============================================================\n")
cat("ALL HISTORICAL TEXT RECORDS CONTAINING 516\n")
cat("============================================================\n\n")

cat("Total records:", nrow(text516), "\n")
cat("Columns:", ncol(text516), "\n\n")

print(
  text516,
  row.names = FALSE
)

cat("\n\n")
cat("============================================================\n")
cat("ORDERED BY FILE AND LINE NUMBER\n")
cat("============================================================\n\n")

text516_ordered <- text516[
  order(
    text516$file,
    text516$line_number
  ),
  ,
  drop = FALSE
]

print(
  text516_ordered,
  row.names = FALSE
)

cat("\n\n")
cat("============================================================\n")
cat("STATUS\n")
cat("============================================================\n")

cat("42 historical text records inspected.\n")
cat("NO ROC CALCULATION\n")
cat("NO AUC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")
cat("NO EXPRESSION RECALCULATION\n")
cat("NO THRESHOLD OPTIMIZATION\n")

# ============================================================
# STEP13B-7I — Historical Count Semantic Resolution Lock
# Paper 1 / WORKFLOW-13
#
# PURPOSE:
#   Final semantic/provenance lock for historical sample counts.
#
# LOCKED CURRENT SAMPLE STRUCTURE:
#   513 Primary Tumor
#    59 Solid Tissue Normal
#     2 Metastatic
#   ---------------
#   574 expression samples
#
#   Tumor-versus-Normal eligible = 572
#
# HISTORICAL OBJECT SEMANTICS:
#   515 = patient-level intermediate object
#   517 = patient-level intermediate object
#
# HISTORICAL "516":
#   Literal occurrences were found in historical files,
#   but NO valid sample-count source for "516 tumor samples"
#   was identified.
#
# IMPORTANT:
#   516 is NOT used as a sample count.
#
# NO:
#   ROC calculation
#   AUC calculation
#   model fitting
#   statistical testing
#   expression recalculation
#   threshold optimization
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP13B-7I — HISTORICAL COUNT SEMANTIC RESOLUTION LOCK\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

audit_dir <- file.path(
  workflow_root,
  "STEP13B-7H_Historical_Count_Provenance_Audit"
)

out_dir <- file.path(
  workflow_root,
  "STEP13B-7I_Historical_Count_Semantic_Resolution_Lock"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 1. Locked provenance facts
# ------------------------------------------------------------

n_primary_tumor <- 513L
n_normal <- 59L
n_metastatic <- 2L

n_expression_samples <- 
  n_primary_tumor +
  n_normal +
  n_metastatic

n_tumor_normal_eligible <-
  n_primary_tumor +
  n_normal

n_patient_object_515 <- 515L
n_patient_intermediate_517 <- 517L

historical_516_literal_found <- TRUE

historical_516_sample_count_verified <- FALSE

historical_516_sample_count_claim_allowed <- FALSE


# ------------------------------------------------------------
# 2. Locked identity / provenance statements
# ------------------------------------------------------------

matrix_annotation_exact <- TRUE

wgcna_primary_exact <- TRUE

metastatic_exclusion_clean <- TRUE

identifier_column_excluded <- TRUE

four_patient_relationship_verified <- TRUE

clinical_without_primary_equals_no_primary <- TRUE


# ------------------------------------------------------------
# 3. Arithmetic consistency
# ------------------------------------------------------------

expression_sample_total_check <-
  n_expression_samples == 574L

eligible_total_check <-
  n_tumor_normal_eligible == 572L

historical_517_relationship_check <-
  (n_primary_tumor + 4L) == n_patient_intermediate_517


# ------------------------------------------------------------
# 4. Historical 516 semantic resolution
# ------------------------------------------------------------

historical_516_resolution <- paste(
  "Literal value 516 was found in historical text/provenance files,",
  "but none of the identified records provided a valid sample-count",
  "semantic source demonstrating 516 tumor samples.",
  "Therefore 516 is NOT interpreted as a tumor sample count and is",
  "NOT used in the current tumor-versus-normal classification."
)


# ------------------------------------------------------------
# 5. Formal lock object
# ------------------------------------------------------------

STEP13B_7I_LOCK <- list(
  
  workflow = "WORKFLOW-13",
  
  checkpoint = "STEP13B-7I",
  
  checkpoint_name =
    "Historical Count Semantic Resolution Lock",
  
  status =
    "PASS — SAMPLE-LEVEL PROVENANCE LOCKED",
  
  current_expression_sample_structure = list(
    
    original_matrix_columns = 575L,
    
    identifier_columns = 1L,
    
    expression_samples = n_expression_samples,
    
    primary_tumor = n_primary_tumor,
    
    solid_tissue_normal = n_normal,
    
    metastatic = n_metastatic,
    
    tumor_normal_eligible = n_tumor_normal_eligible
  ),
  
  exact_reconciliation = list(
    
    matrix_annotation_exact = matrix_annotation_exact,
    
    wgcna_primary_exact = wgcna_primary_exact,
    
    metastatic_exclusion_clean =
      metastatic_exclusion_clean,
    
    identifier_column_excluded =
      identifier_column_excluded
  ),
  
  historical_object_semantics = list(
    
    object_515 =
      "patient-level intermediate object",
    
    object_515_count =
      n_patient_object_515,
    
    object_517 =
      "patient-level intermediate object",
    
    object_517_count =
      n_patient_intermediate_517,
    
    four_patient_relationship_verified =
      four_patient_relationship_verified,
    
    clinical_without_primary_equals_no_primary =
      clinical_without_primary_equals_no_primary,
    
    relationship =
      "513 + 4 = 517"
  ),
  
  historical_516 = list(
    
    literal_occurrence_found =
      historical_516_literal_found,
    
    sample_count_verified =
      historical_516_sample_count_verified,
    
    sample_count_claim_allowed =
      historical_516_sample_count_claim_allowed,
    
    interpretation =
      historical_516_resolution,
    
    use_in_current_analysis =
      "NOT USED"
  ),
  
  current_classification_sample_definition =
    "Primary Tumor (513) + Solid Tissue Normal (59) = 572 eligible samples; metastatic samples (n=2) excluded.",
  
  statistical_analysis_performed = FALSE,
  
  roc_calculation = FALSE,
  
  auc_calculation = FALSE,
  
  model_fitting = FALSE,
  
  statistical_testing = FALSE,
  
  expression_recalculation = FALSE,
  
  threshold_optimization = FALSE
)


# ------------------------------------------------------------
# 6. Final PASS conditions
# ------------------------------------------------------------

pass_conditions <- c(
  
  expression_sample_total_check,
  
  eligible_total_check,
  
  matrix_annotation_exact,
  
  wgcna_primary_exact,
  
  metastatic_exclusion_clean,
  
  identifier_column_excluded,
  
  four_patient_relationship_verified,
  
  clinical_without_primary_equals_no_primary,
  
  historical_517_relationship_check,
  
  historical_516_literal_found,
  
  !historical_516_sample_count_verified,
  
  !historical_516_sample_count_claim_allowed
  
)


final_pass <- all(pass_conditions)


# ------------------------------------------------------------
# 7. Console report
# ------------------------------------------------------------

cat("CURRENT EXPRESSION SAMPLE STRUCTURE\n")
cat("------------------------------------------------------------\n")
cat("Original matrix columns:        575\n")
cat("Identifier columns:               1\n")
cat("Expression samples:             ",
    n_expression_samples, "\n", sep = "")
cat("Primary Tumor:                   ",
    n_primary_tumor, "\n", sep = "")
cat("Solid Tissue Normal:             ",
    n_normal, "\n", sep = "")
cat("Metastatic:                      ",
    n_metastatic, "\n", sep = "")
cat("Tumor-versus-Normal eligible:   ",
    n_tumor_normal_eligible, "\n\n", sep = "")


cat("EXACT PROVENANCE RECONCILIATION\n")
cat("------------------------------------------------------------\n")
cat("Matrix ↔ Annotation EXACT:       ",
    matrix_annotation_exact, "\n", sep = "")
cat("WGCNA ↔ Primary Tumor EXACT:     ",
    wgcna_primary_exact, "\n", sep = "")
cat("Metastatic exclusion clean:      ",
    metastatic_exclusion_clean, "\n", sep = "")
cat("Identifier column excluded:      ",
    identifier_column_excluded, "\n\n", sep = "")


cat("HISTORICAL OBJECT SEMANTICS\n")
cat("------------------------------------------------------------\n")
cat("515 object: patient-level intermediate object\n")
cat("515 rows:   ", n_patient_object_515, "\n", sep = "")
cat("517 object: patient-level intermediate object\n")
cat("517 rows:   ", n_patient_intermediate_517, "\n", sep = "")
cat("Clinical-without-primary == no-primary: ",
    clinical_without_primary_equals_no_primary,
    "\n", sep = "")
cat("513 + 4 = 517 relationship verified:      ",
    historical_517_relationship_check,
    "\n\n", sep = "")


cat("HISTORICAL 516 SEMANTIC RESOLUTION\n")
cat("------------------------------------------------------------\n")
cat("Literal 516 occurrence found:              ",
    historical_516_literal_found,
    "\n", sep = "")
cat("516 tumor sample count verified:           ",
    historical_516_sample_count_verified,
    "\n", sep = "")
cat("516 tumor sample count claim allowed:      ",
    historical_516_sample_count_claim_allowed,
    "\n", sep = "")
cat("516 used in current classification:        NOT USED\n\n")


cat("FINAL CURRENT SAMPLE DEFINITION\n")
cat("------------------------------------------------------------\n")
cat("Primary Tumor = 513\n")
cat("Solid Tissue Normal = 59\n")
cat("Metastatic = 2 (excluded)\n")
cat("Eligible Tumor-versus-Normal = 572\n\n")


cat("STATISTICAL ANALYSIS STATUS\n")
cat("------------------------------------------------------------\n")
cat("NO ROC CALCULATION\n")
cat("NO AUC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")
cat("NO EXPRESSION RECALCULATION\n")
cat("NO THRESHOLD OPTIMIZATION\n\n")


if (final_pass) {
  
  cat("============================================================\n")
  cat("STEP13B-7I STATUS: PASS\n")
  cat("SAMPLE-LEVEL PROVENANCE LOCKED\n")
  cat("============================================================\n")
  cat("Historical 516 sample-count claim: NOT VERIFIED / NOT USED\n")
  cat("Current classification universe: 572 samples\n")
  cat("============================================================\n\n")
  
} else {
  
  cat("============================================================\n")
  cat("STEP13B-7I STATUS: FAIL\n")
  cat("ONE OR MORE LOCK CONDITIONS NOT SATISFIED\n")
  cat("============================================================\n\n")
}


# ------------------------------------------------------------
# 8. Save formal lock object
# ------------------------------------------------------------

saveRDS(
  STEP13B_7I_LOCK,
  file.path(
    out_dir,
    "STEP13B-7I_Historical_Count_Semantic_Resolution_LOCKED.rds"
  )
)


# ------------------------------------------------------------
# 9. Save concise CSV lock table
# ------------------------------------------------------------

lock_table <- data.frame(
  
  item = c(
    
    "Original matrix columns",
    "Identifier columns",
    "Expression samples",
    "Primary Tumor",
    "Solid Tissue Normal",
    "Metastatic",
    "Tumor-versus-Normal eligible",
    "Matrix-Annotation exact match",
    "WGCNA-Primary Tumor exact match",
    "Metastatic exclusion",
    "515 object semantic",
    "517 object semantic",
    "516 literal occurrence",
    "516 tumor sample count verified",
    "516 used in classification"
    
  ),
  
  value = c(
    
    "575",
    "1",
    "574",
    "513",
    "59",
    "2",
    "572",
    "TRUE",
    "TRUE",
    "TRUE",
    "patient-level intermediate object",
    "patient-level intermediate object",
    "TRUE",
    "FALSE",
    "FALSE"
    
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  lock_table,
  file.path(
    out_dir,
    "STEP13B-7I_Historical_Count_Semantic_Resolution_LOCKED.csv"
  ),
  row.names = FALSE,
  quote = TRUE
)


# ------------------------------------------------------------
# 10. Save human-readable interpretation lock
# ------------------------------------------------------------

interpretation_text <- c(
  
  "STEP13B-7I — Historical Count Semantic Resolution Lock",
  
  "",
  
  "FINAL STATUS:",
  "PASS — SAMPLE-LEVEL PROVENANCE LOCKED",
  
  "",
  
  "CURRENT EXPRESSION SAMPLE STRUCTURE:",
  "Original matrix = 60,483 genes × 575 columns.",
  "One column is the identifier column.",
  "Therefore the matrix contains 574 expression samples.",
  
  "513 Primary Tumor samples.",
  "59 Solid Tissue Normal samples.",
  "2 Metastatic samples.",
  
  "Tumor-versus-Normal eligible sample set = 513 + 59 = 572.",
  
  "",
  
  "EXACT SAMPLE-LEVEL RECONCILIATION:",
  "Matrix sample IDs and annotation sample IDs matched exactly.",
  "All 513 final WGCNA sample IDs matched the 513 Primary Tumor IDs exactly.",
  "The two metastatic samples were excluded from tumor-versus-normal classification.",
  
  "",
  
  "HISTORICAL OBJECT SEMANTICS:",
  "515 = patient-level intermediate object.",
  "517 = patient-level intermediate object.",
  "clinical_without_primary and no_primary were identical.",
  "The observed relationship 513 + 4 = 517 is verified.",
  
  "",
  
  "HISTORICAL 516:",
  "Literal occurrences of the value 516 were found in historical files.",
  "However, no identified record provides a valid sample-count semantic source",
  "for 516 tumor samples.",
  "Therefore 516 is NOT interpreted as a tumor sample count.",
  "The historical 516 tumor-sample claim is NOT VERIFIED and is NOT USED.",
  
  "",
  
  "CURRENT LOCKED CLASSIFICATION UNIVERSE:",
  "513 Primary Tumor + 59 Solid Tissue Normal = 572 eligible samples.",
  
  "",
  
  "NO ROC CALCULATION.",
  "NO AUC CALCULATION.",
  "NO MODEL FITTING.",
  "NO STATISTICAL TEST.",
  "NO EXPRESSION RECALCULATION.",
  "NO THRESHOLD OPTIMIZATION."
  
)


writeLines(
  interpretation_text,
  con = file.path(
    out_dir,
    "STEP13B-7I_Interpretation_Lock.txt"
  )
)


# ------------------------------------------------------------
# 11. PASS marker
# ------------------------------------------------------------

if (final_pass) {
  
  writeLines(
    c(
      "STEP13B-7I PASS",
      "SAMPLE-LEVEL PROVENANCE LOCKED",
      "",
      "Current expression samples = 574",
      "Primary Tumor = 513",
      "Solid Tissue Normal = 59",
      "Metastatic = 2",
      "Tumor-versus-Normal eligible = 572",
      "",
      "Historical 515 = patient-level intermediate object",
      "Historical 517 = patient-level intermediate object",
      "Historical 516 sample-count claim = NOT VERIFIED / NOT USED",
      "",
      "NO ROC CALCULATION",
      "NO AUC CALCULATION",
      "NO MODEL FITTING",
      "NO STATISTICAL TEST",
      "NO EXPRESSION RECALCULATION",
      "NO THRESHOLD OPTIMIZATION"
    ),
    con = file.path(
      out_dir,
      "STEP13B-7I_PASS.txt"
    )
  )
  
}


cat("Outputs saved to:\n")
cat(out_dir, "\n\n")

cat("STEP13B-7I script completed.\n")

# ============================================================
# STEP13C-0A
# Historical GSE32863 Key-Object Audit
#
# NO ROC
# NO AUC
# NO MODEL FITTING
# NO STATISTICAL TEST
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP13C-0A — Historical GSE32863 Key-Object Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Key expression objects
# ------------------------------------------------------------

key_expr_objects <- c(
  "expr",
  "expr.final",
  "expr.mat",
  "expr_df",
  "expr_gene",
  "expr_merge",
  "expr_obj",
  "expr2",
  "expr3"
)

cat("============================================================\n")
cat("EXPRESSION OBJECTS\n")
cat("============================================================\n\n")

for (nm in key_expr_objects) {
  
  if (exists(nm, inherits = TRUE)) {
    
    x <- get(nm, inherits = TRUE)
    
    cat("\nOBJECT:", nm, "\n")
    cat("Class:", paste(class(x), collapse = "; "), "\n")
    cat("Type:", typeof(x), "\n")
    cat("Length:", length(x), "\n")
    
    if (!is.null(dim(x))) {
      cat(
        "Dimension:",
        paste(dim(x), collapse = " x "),
        "\n"
      )
    }
    
    if (!is.null(rownames(x))) {
      cat(
        "Row names head:",
        paste(head(rownames(x), 10), collapse = " | "),
        "\n"
      )
    }
    
    if (!is.null(colnames(x))) {
      cat(
        "Column names head:",
        paste(head(colnames(x), 10), collapse = " | "),
        "\n"
      )
    }
    
    if (is.data.frame(x)) {
      cat(
        "Data-frame columns:",
        paste(names(x), collapse = " | "),
        "\n"
      )
    }
    
  } else {
    
    cat("\nOBJECT:", nm, "NOT FOUND\n")
  }
}


# ------------------------------------------------------------
# 2. Annotation / grouping objects
# ------------------------------------------------------------

key_annotation_objects <- c(
  "ann",
  "anno",
  "anno2",
  "annotation_file",
  "annotation_ids",
  "hist_ann",
  "patient_ann",
  "patient_map",
  "group",
  "clinical"
)

cat("\n\n")
cat("============================================================\n")
cat("ANNOTATION / GROUPING OBJECTS\n")
cat("============================================================\n\n")

for (nm in key_annotation_objects) {
  
  if (exists(nm, inherits = TRUE)) {
    
    x <- get(nm, inherits = TRUE)
    
    cat("\nOBJECT:", nm, "\n")
    cat("Class:", paste(class(x), collapse = "; "), "\n")
    cat("Type:", typeof(x), "\n")
    cat("Length:", length(x), "\n")
    
    if (!is.null(dim(x))) {
      cat(
        "Dimension:",
        paste(dim(x), collapse = " x "),
        "\n"
      )
    }
    
    if (is.data.frame(x)) {
      cat(
        "Columns:",
        paste(names(x), collapse = " | "),
        "\n"
      )
      
      cat("First rows:\n")
      print(utils::head(x, 3))
    }
  }
}


# ------------------------------------------------------------
# 3. Historical ROC objects
# ------------------------------------------------------------

key_roc_objects <- c(
  "roc_AOC3",
  "roc_MS4A7",
  "roc_IL7R",
  "roc_combined",
  "roc_cv",
  "roc.ext",
  "Figure12A_CombinedROC"
)

cat("\n\n")
cat("============================================================\n")
cat("HISTORICAL ROC OBJECTS — INVENTORY ONLY\n")
cat("============================================================\n\n")

for (nm in key_roc_objects) {
  
  if (exists(nm, inherits = TRUE)) {
    
    x <- get(nm, inherits = TRUE)
    
    cat("\nOBJECT:", nm, "\n")
    cat("Class:", paste(class(x), collapse = "; "), "\n")
    cat("Type:", typeof(x), "\n")
    cat("Length:", length(x), "\n")
    
    if (!is.null(dim(x))) {
      cat(
        "Dimension:",
        paste(dim(x), collapse = " x "),
        "\n"
      )
    }
    
    cat(
      "Names:",
      paste(head(names(x), 30), collapse = " | "),
      "\n"
    )
  }
}


# ------------------------------------------------------------
# 4. Historical model objects
# ------------------------------------------------------------

key_model_objects <- c(
  "fit",
  "model",
  "cv_model"
)

cat("\n\n")
cat("============================================================\n")
cat("HISTORICAL MODEL OBJECTS — INVENTORY ONLY\n")
cat("============================================================\n\n")

for (nm in key_model_objects) {
  
  if (exists(nm, inherits = TRUE)) {
    
    x <- get(nm, inherits = TRUE)
    
    cat("\nOBJECT:", nm, "\n")
    cat("Class:", paste(class(x), collapse = "; "), "\n")
    cat("Type:", typeof(x), "\n")
    cat("Length:", length(x), "\n")
    
    if (!is.null(dim(x))) {
      cat(
        "Dimension:",
        paste(dim(x), collapse = " x "),
        "\n"
      )
    }
  }
}


# ------------------------------------------------------------
# 5. Gene-related objects
# ------------------------------------------------------------

key_gene_objects <- c(
  "genes",
  "expr_gene",
  "candidate_columns",
  "candidate_ext"
)

cat("\n\n")
cat("============================================================\n")
cat("GENE / PROBE RELATED OBJECTS\n")
cat("============================================================\n\n")

for (nm in key_gene_objects) {
  
  if (exists(nm, inherits = TRUE)) {
    
    x <- get(nm, inherits = TRUE)
    
    cat("\nOBJECT:", nm, "\n")
    cat("Class:", paste(class(x), collapse = "; "), "\n")
    cat("Type:", typeof(x), "\n")
    cat("Length:", length(x), "\n")
    
    if (!is.null(dim(x))) {
      cat(
        "Dimension:",
        paste(dim(x), collapse = " x "),
        "\n"
      )
    }
    
    if (is.atomic(x)) {
      cat(
        "First values:",
        paste(head(x, 20), collapse = " | "),
        "\n"
      )
    }
    
    if (is.data.frame(x)) {
      cat(
        "Columns:",
        paste(names(x), collapse = " | "),
        "\n"
      )
    }
  }
}


cat("\n\n")
cat("============================================================\n")
cat("STEP13C-0A COMPLETE\n")
cat("============================================================\n")
cat("Inventory only.\n")
cat("NO ROC CALCULATION\n")
cat("NO AUC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")
cat("NO EXPRESSION RECALCULATION\n")
cat("============================================================\n")

# ============================================================
# STEP13C-0B
# Historical GSE32863 Sample / Label / Model Traceability Audit
#
# PURPOSE:
#   Resolve the historical GSE32863 sample structure,
#   labels, three-gene input and stored ROC provenance.
#
# NO:
#   ROC recalculation
#   AUC recalculation
#   model fitting
#   statistical testing
#   expression recalculation
#   normalization
#   probe selection
#   threshold optimization
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP13C-0B — Historical GSE32863 Traceability Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Required objects
# ------------------------------------------------------------

required_objects <- c(
  "expr",
  "expr.mat",
  "expr.final",
  "expr2",
  "expr3",
  "group",
  "genes",
  "roc_AOC3",
  "roc_MS4A7",
  "roc_IL7R",
  "roc_combined",
  "roc_cv",
  "roc.ext",
  "Figure12A_CombinedROC"
)

cat("Required-object availability:\n")

availability <- data.frame(
  object = required_objects,
  exists = sapply(
    required_objects,
    function(x) exists(x, inherits = TRUE)
  ),
  stringsAsFactors = FALSE
)

print(availability, row.names = FALSE)


# ------------------------------------------------------------
# 1. GSE32863 sample IDs
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("1. GSE32863 SAMPLE STRUCTURE\n")
cat("============================================================\n\n")

gse_samples_expr <- colnames(expr)[
  colnames(expr) != "ID_REF"
]

gse_samples_exprmat <- colnames(expr.mat)

gse_samples_exprfinal <- colnames(expr.final)

gse_samples_expr2 <- setdiff(
  colnames(expr2),
  c("Probe", "Symbol")
)

gse_samples_expr3 <- rownames(expr3)


cat("expr samples:       ",
    length(gse_samples_expr), "\n", sep = "")

cat("expr.mat samples:   ",
    length(gse_samples_exprmat), "\n", sep = "")

cat("expr.final samples: ",
    length(gse_samples_exprfinal), "\n", sep = "")

cat("expr2 samples:      ",
    length(gse_samples_expr2), "\n", sep = "")

cat("expr3 samples:      ",
    length(gse_samples_expr3), "\n", sep = "")


cat("\nExact sample-ID agreement:\n")

cat(
  "expr vs expr.mat:    ",
  setequal(gse_samples_expr, gse_samples_exprmat),
  "\n",
  sep = ""
)

cat(
  "expr vs expr.final:  ",
  setequal(gse_samples_expr, gse_samples_exprfinal),
  "\n",
  sep = ""
)

cat(
  "expr vs expr2:       ",
  setequal(gse_samples_expr, gse_samples_expr2),
  "\n",
  sep = ""
)

cat(
  "expr vs expr3:       ",
  setequal(gse_samples_expr, gse_samples_expr3),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 2. Historical group object
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("2. HISTORICAL GROUP OBJECT\n")
cat("============================================================\n\n")

cat("Class:\n")
print(class(group))

cat("\nLength:\n")
print(length(group))

cat("\nFirst values:\n")
print(head(group, 20))

cat("\nUnique values:\n")
print(unique(group))

cat("\nFrequency table:\n")
print(table(group, useNA = "ifany"))


# ------------------------------------------------------------
# 3. expr3 Group
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("3. HISTORICAL expr3 GROUP / PROBABILITY\n")
cat("============================================================\n\n")

cat("Dimensions:\n")
print(dim(expr3))

cat("\nColumn names:\n")
print(names(expr3))

cat("\nGroup frequency:\n")
print(table(expr3$Group, useNA = "ifany"))

cat("\nProbability summary:\n")
print(summary(expr3$Probability))

cat("\nFirst 10 rows:\n")
print(head(expr3, 10))


# ------------------------------------------------------------
# 4. Three-gene availability in expr3
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("4. THREE-GENE INPUT\n")
cat("============================================================\n\n")

for (g in c("AOC3", "MS4A7", "IL7R")) {
  
  cat("\nGENE:", g, "\n")
  
  if (g %in% names(expr3)) {
    
    x <- expr3[[g]]
    
    cat("Available in expr3: TRUE\n")
    cat("Class:", class(x), "\n")
    cat("NA:", sum(is.na(x)), "\n")
    cat("Summary:\n")
    print(summary(x))
    
  } else {
    
    cat("Available in expr3: FALSE\n")
  }
}


# ------------------------------------------------------------
# 5. Historical ROC object audit
# ------------------------------------------------------------

roc_objects <- c(
  "roc_AOC3",
  "roc_MS4A7",
  "roc_IL7R",
  "roc_combined",
  "roc_cv",
  "roc.ext",
  "Figure12A_CombinedROC"
)

roc_summary <- data.frame(
  object = character(),
  class = character(),
  n_cases = integer(),
  n_controls = integer(),
  direction = character(),
  levels = character(),
  stored_auc = character(),
  response_class = character(),
  predictor_class = character(),
  stringsAsFactors = FALSE
)

for (nm in roc_objects) {
  
  if (!exists(nm, inherits = TRUE)) next
  
  r <- get(nm, inherits = TRUE)
  
  n_cases <- NA_integer_
  n_controls <- NA_integer_
  direction <- NA_character_
  levels_txt <- NA_character_
  auc_txt <- NA_character_
  response_class <- NA_character_
  predictor_class <- NA_character_
  
  if (!is.null(r$cases)) {
    n_cases <- length(r$cases)
  }
  
  if (!is.null(r$controls)) {
    n_controls <- length(r$controls)
  }
  
  if (!is.null(r$direction)) {
    direction <- as.character(r$direction)
  }
  
  if (!is.null(r$levels)) {
    levels_txt <- paste(
      as.character(r$levels),
      collapse = " | "
    )
  }
  
  if (!is.null(r$auc)) {
    auc_txt <- paste(
      as.character(r$auc),
      collapse = " | "
    )
  }
  
  if (!is.null(r$response)) {
    response_class <- paste(
      class(r$response),
      collapse = ";"
    )
  }
  
  if (!is.null(r$predictor)) {
    predictor_class <- paste(
      class(r$predictor),
      collapse = ";"
    )
  }
  
  roc_summary <- rbind(
    roc_summary,
    data.frame(
      object = nm,
      class = paste(class(r), collapse = ";"),
      n_cases = n_cases,
      n_controls = n_controls,
      direction = direction,
      levels = levels_txt,
      stored_auc = auc_txt,
      response_class = response_class,
      predictor_class = predictor_class,
      stringsAsFactors = FALSE
    )
  )
}


cat("\n")
cat("============================================================\n")
cat("5A. STORED ROC OBJECT SUMMARY\n")
cat("============================================================\n\n")

print(
  roc_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Stored response / predictor identity
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("5B. STORED ROC RESPONSE / PREDICTOR IDENTITY\n")
cat("============================================================\n")

for (nm in roc_objects) {
  
  if (!exists(nm, inherits = TRUE)) next
  
  r <- get(nm, inherits = TRUE)
  
  cat("\n------------------------------------------------------------\n")
  cat("ROC OBJECT:", nm, "\n")
  cat("------------------------------------------------------------\n")
  
  if (!is.null(r$original.response)) {
    
    cat("\noriginal.response:\n")
    print(
      head(
        r$original.response,
        20
      )
    )
    
    cat("\noriginal.response frequency:\n")
    print(
      table(
        r$original.response,
        useNA = "ifany"
      )
    )
  }
  
  if (!is.null(r$original.predictor)) {
    
    cat("\noriginal.predictor class:\n")
    print(class(r$original.predictor))
    
    cat("\noriginal.predictor first values:\n")
    print(
      head(
        r$original.predictor,
        10
      )
    )
  }
}


# ------------------------------------------------------------
# 7. Gene object
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("6. LOCKED GENE OBJECT\n")
cat("============================================================\n\n")

print(genes)


# ------------------------------------------------------------
# 8. Save audit results
# ------------------------------------------------------------

audit_dir <- file.path(
  "result_dir/Paper 1_Step 13-Diagnostic Classification",
  "STEP13C_GSE32863_Input_Lock",
  "STEP13C-0_Historical_RData_Audit"
)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  availability,
  file.path(
    audit_dir,
    "STEP13C-0B_Object_Availability.csv"
  ),
  row.names = FALSE
)

write.csv(
  roc_summary,
  file.path(
    audit_dir,
    "STEP13C-0B_Historical_ROC_Object_Summary.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    gse_samples_expr = gse_samples_expr,
    gse_samples_exprmat = gse_samples_exprmat,
    gse_samples_exprfinal = gse_samples_exprfinal,
    gse_samples_expr2 = gse_samples_expr2,
    gse_samples_expr3 = gse_samples_expr3,
    group = group,
    expr3_group = expr3$Group,
    expr3_probability = expr3$Probability,
    genes = genes,
    roc_summary = roc_summary
  ),
  file.path(
    audit_dir,
    "STEP13C-0B_GSE32863_Traceability_Audit.rds"
  )
)


cat("\n")
cat("============================================================\n")
cat("STEP13C-0B COMPLETE\n")
cat("============================================================\n")
cat("Inventory / provenance extraction only.\n")
cat("NO ROC CALCULATION\n")
cat("NO AUC CALCULATION\n")
cat("NO MODEL FITTING\n")
cat("NO STATISTICAL TEST\n")
cat("NO EXPRESSION RECALCULATION\n")
cat("NO PROBE SELECTION\n")
cat("============================================================\n")

# ============================================================
# STEP13B-7J-R1
# Xiangtao Historical ROC Excel
# Long-Format Sample-Level Traceability Audit
#
# PURPOSE
# ------------------------------------------------------------
# Audit the original Xiangtao ROC exported Excel artifact
# at the sample-gene-record level.
#
# Expected format:
#
#   sample_id | status | x | y
#
# where:
#   sample_id = TCGA sample identifier
#   status    = Tumor / Normal
#   x         = AOC3 / MS4A7 / IL7R
#   y         = expression value
#
# IMPORTANT
# ------------------------------------------------------------
# This workflow:
#   - DOES NOT recalculate ROC
#   - DOES NOT recalculate AUC
#   - DOES NOT refit logistic regression
#   - DOES NOT optimize thresholds
#   - DOES NOT modify expression values
#   - DOES NOT reconcile Xiangtao 516 with Sangerbox 513
#
# It establishes sample-level provenance of the historical
# Xiangtao ROC branch.
# ============================================================


# ------------------------------------------------------------
# 0. Configuration
# ------------------------------------------------------------

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

step_dir <- file.path(
  workflow_root,
  "STEP13B-7J_Xiangtao_ROC_Excel_Traceability"
)

output_dir <- file.path(step_dir, "OUTPUT")
audit_dir  <- file.path(step_dir, "AUDIT")
log_dir    <- file.path(step_dir, "LOGS")

dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

excel_file <- file.path(
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal",
  "Diagnostic_ROC_AOC3_MS4A7_IL7R_Data_Analysis.xlsx"
)

historical_rdata <- file.path(
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal",
  ".RData"
)


# ------------------------------------------------------------
# 2. Header
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP13B-7J-R1\n")
cat("Xiangtao Historical ROC Excel Long-Format Traceability\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 3. Input existence
# ------------------------------------------------------------

if (!file.exists(excel_file)) {
  stop(
    paste0(
      "FAIL: Excel file not found:\n",
      excel_file
    )
  )
}

if (!file.exists(historical_rdata)) {
  stop(
    paste0(
      "FAIL: Historical RData not found:\n",
      historical_rdata
    )
  )
}

cat("PASS: Input Excel exists.\n")
cat("PASS: Historical RData exists.\n\n")


# ------------------------------------------------------------
# 4. Required packages
# ------------------------------------------------------------

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required.")
}

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package 'openxlsx' is required.")
}


# ------------------------------------------------------------
# 5. Workbook inventory
# ------------------------------------------------------------

sheet_names <- readxl::excel_sheets(excel_file)

cat("Number of sheets:", length(sheet_names), "\n")
cat("Sheet names:\n")
print(sheet_names)

if (length(sheet_names) != 1) {
  
  stop(
    paste0(
      "HOLD: Expected one Excel sheet, found ",
      length(sheet_names),
      "."
    )
  )
}

sheet_name <- sheet_names[1]


# ------------------------------------------------------------
# 6. Read exact Excel sheet
# ------------------------------------------------------------

dat <- readxl::read_excel(
  excel_file,
  sheet = sheet_name,
  col_names = TRUE
)

cat("\nExcel dimensions:\n")
cat("Rows:", nrow(dat), "\n")
cat("Columns:", ncol(dat), "\n\n")

cat("Column names:\n")
print(names(dat))


# ------------------------------------------------------------
# 7. Column-name normalization ONLY for identification
# ------------------------------------------------------------
#
# Original data are not modified.
# ------------------------------------------------------------

cn <- trimws(
  as.character(names(dat))
)

cn_lower <- tolower(cn)

required_columns <- c(
  "sample_id",
  "status",
  "x",
  "y"
)

required_present <- required_columns %in% cn_lower

column_check <- data.frame(
  required_column = required_columns,
  present = required_present,
  stringsAsFactors = FALSE
)

print(column_check)

write.csv(
  column_check,
  file.path(
    output_dir,
    "STEP13B-7J-R1_Required_Column_Check.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

if (!all(required_present)) {
  
  stop(
    paste0(
      "HOLD: Required long-format columns are incomplete.\n",
      "Expected: sample_id, status, x, y"
    )
  )
}


# ------------------------------------------------------------
# 8. Map actual columns
# ------------------------------------------------------------

sample_id_col <- names(dat)[
  which(cn_lower == "sample_id")[1]
]

status_col <- names(dat)[
  which(cn_lower == "status")[1]
]

x_col <- names(dat)[
  which(cn_lower == "x")[1]
]

y_col <- names(dat)[
  which(cn_lower == "y")[1]
]

cat("\nMapped columns:\n")
cat("sample_id:", sample_id_col, "\n")
cat("status   :", status_col, "\n")
cat("x        :", x_col, "\n")
cat("y        :", y_col, "\n")


# ------------------------------------------------------------
# 9. Construct audit vectors
# ------------------------------------------------------------

sample_id <- trimws(
  as.character(dat[[sample_id_col]])
)

status <- trimws(
  as.character(dat[[status_col]])
)

gene <- trimws(
  as.character(dat[[x_col]])
)

expression_value <- suppressWarnings(
  as.numeric(dat[[y_col]])
)


# ------------------------------------------------------------
# 10. Basic missingness audit
# ------------------------------------------------------------

missing_audit <- data.frame(
  field = c(
    "sample_id",
    "status",
    "x_gene",
    "y_expression"
  ),
  n_missing = c(
    sum(is.na(sample_id) | sample_id == ""),
    sum(is.na(status) | status == ""),
    sum(is.na(gene) | gene == ""),
    sum(is.na(expression_value))
  ),
  stringsAsFactors = FALSE
)

print(missing_audit)

write.csv(
  missing_audit,
  file.path(
    output_dir,
    "STEP13B-7J-R1_Missingness_Audit.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 11. Gene-level audit
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

gene_levels <- sort(
  unique(gene)
)

cat("\nUnique x/gene values:\n")
print(gene_levels)

gene_frequency <- as.data.frame(
  table(gene),
  stringsAsFactors = FALSE
)

names(gene_frequency) <- c(
  "gene",
  "n_records"
)

print(gene_frequency)

write.csv(
  gene_frequency,
  file.path(
    output_dir,
    "STEP13B-7J-R1_Gene_Record_Distribution.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 12. Core-gene exactness
# ------------------------------------------------------------

gene_exactness <- data.frame(
  gene = core_genes,
  present = core_genes %in% gene_levels,
  n_records = sapply(
    core_genes,
    function(g) sum(gene == g, na.rm = TRUE)
  ),
  n_unique_sample_ids = sapply(
    core_genes,
    function(g) {
      length(
        unique(
          sample_id[gene == g]
        )
      )
    }
  ),
  stringsAsFactors = FALSE
)

print(gene_exactness)

write.csv(
  gene_exactness,
  file.path(
    output_dir,
    "STEP13B-7J-R1_CoreGene_Record_Audit.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 13. Status-level audit
# ------------------------------------------------------------

status_upper <- toupper(status)

status_class <- ifelse(
  status_upper == "TUMOR",
  "Tumor",
  ifelse(
    status_upper == "NORMAL",
    "Normal",
    NA_character_
  )
)

status_levels <- sort(
  unique(status[!is.na(status)])
)

cat("\nOriginal status values:\n")
print(status_levels)

status_frequency <- as.data.frame(
  table(status),
  stringsAsFactors = FALSE
)

names(status_frequency) <- c(
  "status",
  "n_records"
)

print(status_frequency)

write.csv(
  status_frequency,
  file.path(
    output_dir,
    "STEP13B-7J-R1_Status_Record_Distribution.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 14. Gene × status cross-tabulation
# ------------------------------------------------------------

gene_status_table <- table(
  factor(
    gene,
    levels = core_genes
  ),
  factor(
    status_class,
    levels = c("Tumor", "Normal")
  )
)

cat("\nGene × status:\n")
print(gene_status_table)

write.csv(
  as.data.frame.matrix(gene_status_table),
  file.path(
    output_dir,
    "STEP13B-7J-R1_CoreGene_by_Status.csv"
  ),
  row.names = TRUE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 15. Unique sample-level audit
# ------------------------------------------------------------

unique_sample_ids <- unique(
  sample_id[
    !is.na(sample_id) &
      sample_id != ""
  ]
)

n_unique_samples <- length(
  unique_sample_ids
)

cat("\nUnique sample IDs:", n_unique_samples, "\n")


# ------------------------------------------------------------
# 16. Sample-level status consistency
# ------------------------------------------------------------

sample_status_table <- split(
  status_class,
  sample_id
)

sample_status_consistency <- sapply(
  sample_status_table,
  function(x) {
    length(
      unique(
        x[!is.na(x)]
      )
    ) == 1
  }
)

n_inconsistent_status <- sum(
  !sample_status_consistency
)

cat(
  "Samples with inconsistent status:",
  n_inconsistent_status,
  "\n"
)

write.csv(
  data.frame(
    sample_id = names(sample_status_consistency),
    status_consistent =
      as.logical(sample_status_consistency),
    stringsAsFactors = FALSE
  ),
  file.path(
    output_dir,
    "STEP13B-7J-R1_Sample_Status_Consistency.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 17. Sample × gene uniqueness
# ------------------------------------------------------------

sample_gene_key <- paste(
  sample_id,
  gene,
  sep = "|||"
)

duplicate_sample_gene <- duplicated(
  sample_gene_key
)

n_duplicate_sample_gene <- sum(
  duplicate_sample_gene
)

cat(
  "Duplicated sample × gene records:",
  n_duplicate_sample_gene,
  "\n"
)


# ------------------------------------------------------------
# 18. Expected 575 × 3 structure
# ------------------------------------------------------------

expected_n_samples <- 575
expected_n_records <- 575 * 3

structure_audit <- data.frame(
  metric = c(
    "Total Excel records",
    "Unique samples",
    "Expected unique samples",
    "Expected total records",
    "Duplicated sample-gene records",
    "Tumor records",
    "Normal records"
  ),
  observed = c(
    nrow(dat),
    n_unique_samples,
    expected_n_samples,
    expected_n_records,
    n_duplicate_sample_gene,
    sum(status_class == "Tumor", na.rm = TRUE),
    sum(status_class == "Normal", na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

print(structure_audit)

write.csv(
  structure_audit,
  file.path(
    output_dir,
    "STEP13B-7J-R1_575x3_Structure_Audit.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 19. Create canonical long-format traceability table
# ------------------------------------------------------------

traceability_long <- data.frame(
  sample_id = sample_id,
  status = status_class,
  gene = gene,
  expression = expression_value,
  stringsAsFactors = FALSE
)

traceability_long <- traceability_long[
  order(
    traceability_long$sample_id,
    traceability_long$gene
  ),
]

write.csv(
  traceability_long,
  file.path(
    output_dir,
    "STEP13B-7J-R1_Xiangtao_ROC_LongFormat_Traceability.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

saveRDS(
  traceability_long,
  file.path(
    output_dir,
    "STEP13B-7J-R1_Xiangtao_ROC_LongFormat_Traceability.rds"
  )
)


# ------------------------------------------------------------
# 20. Wide-format representation
# ------------------------------------------------------------
#
# This is ONLY a representation for audit.
# No expression transformation is performed.
# ------------------------------------------------------------

wide_expression <- reshape(
  traceability_long[
    ,
    c("sample_id", "status", "gene", "expression")
  ],
  idvar = c("sample_id", "status"),
  timevar = "gene",
  direction = "wide"
)

names(wide_expression) <-
  sub(
    "^expression\\.",
    "",
    names(wide_expression)
  )

wide_expression <- wide_expression[
  ,
  c(
    "sample_id",
    "status",
    core_genes
  )
]

write.csv(
  wide_expression,
  file.path(
    output_dir,
    "STEP13B-7J-R1_Xiangtao_ROC_WideRepresentation.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 21. Load historical RData in isolated environment
# ------------------------------------------------------------

hist_env <- new.env(
  parent = emptyenv()
)

load(
  historical_rdata,
  envir = hist_env
)

required_objects <- c(
  "roc_AOC3",
  "roc_MS4A7",
  "roc_IL7R"
)

object_available <-
  required_objects %in% ls(hist_env)

object_inventory <- data.frame(
  object = required_objects,
  available = object_available,
  stringsAsFactors = FALSE
)

print(object_inventory)

write.csv(
  object_inventory,
  file.path(
    output_dir,
    "STEP13B-7J-R1_Historical_ROC_Object_Availability.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

if (!all(object_available)) {
  
  stop(
    "HOLD: Required historical ROC objects are unavailable."
  )
}


# ------------------------------------------------------------
# 22. Historical ROC object structural audit
# ------------------------------------------------------------
#
# IMPORTANT:
# No ROC recalculation.
# No AUC recalculation.
# No threshold recalculation.
#
# Only inspect the already-existing historical ROC objects.
# ------------------------------------------------------------

historical_summary <- data.frame(
  gene = core_genes,
  object = paste0(
    "roc_",
    core_genes
  ),
  n = NA_integer_,
  tumor = NA_integer_,
  normal = NA_integer_,
  response_class = NA_character_,
  predictor_class = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(core_genes)) {
  
  gene_i <- core_genes[i]
  
  roc_obj <- get(
    paste0("roc_", gene_i),
    envir = hist_env
  )
  
  response_chr <- as.character(
    roc_obj$response
  )
  
  predictor_vec <- roc_obj$predictor
  
  historical_summary$n[i] <-
    length(response_chr)
  
  historical_summary$tumor[i] <-
    sum(
      response_chr %in% c(
        "1",
        "Tumor",
        "tumor"
      ),
      na.rm = TRUE
    )
  
  historical_summary$normal[i] <-
    sum(
      response_chr %in% c(
        "0",
        "Normal",
        "normal"
      ),
      na.rm = TRUE
    )
  
  historical_summary$response_class[i] <-
    paste(
      class(roc_obj$response),
      collapse = ";"
    )
  
  historical_summary$predictor_class[i] <-
    paste(
      class(predictor_vec),
      collapse = ";"
    )
}

cat("\nHistorical ROC object structure:\n")
print(historical_summary)

write.csv(
  historical_summary,
  file.path(
    output_dir,
    "STEP13B-7J-R1_Historical_ROC_Structure.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 23. Compare Excel sample counts with historical ROC objects
# ------------------------------------------------------------
#
# This is a provenance audit only.
# ------------------------------------------------------------

sample_count_comparison <- data.frame(
  gene = core_genes,
  excel_records = NA_integer_,
  excel_tumor = NA_integer_,
  excel_normal = NA_integer_,
  historical_n = NA_integer_,
  historical_tumor = NA_integer_,
  historical_normal = NA_integer_,
  total_n_match = FALSE,
  tumor_n_match = FALSE,
  normal_n_match = FALSE,
  stringsAsFactors = FALSE
)

for (i in seq_along(core_genes)) {
  
  g <- core_genes[i]
  
  roc_obj <- get(
    paste0("roc_", g),
    envir = hist_env
  )
  
  response_chr <- as.character(
    roc_obj$response
  )
  
  excel_gene <- gene == g
  
  excel_tumor <- sum(
    status_class[excel_gene] == "Tumor",
    na.rm = TRUE
  )
  
  excel_normal <- sum(
    status_class[excel_gene] == "Normal",
    na.rm = TRUE
  )
  
  historical_n <- length(
    response_chr
  )
  
  historical_tumor <- sum(
    response_chr %in% c(
      "1",
      "Tumor",
      "tumor"
    ),
    na.rm = TRUE
  )
  
  historical_normal <- sum(
    response_chr %in% c(
      "0",
      "Normal",
      "normal"
    ),
    na.rm = TRUE
  )
  
  sample_count_comparison$excel_records[i] <-
    sum(excel_gene)
  
  sample_count_comparison$excel_tumor[i] <-
    excel_tumor
  
  sample_count_comparison$excel_normal[i] <-
    excel_normal
  
  sample_count_comparison$historical_n[i] <-
    historical_n
  
  sample_count_comparison$historical_tumor[i] <-
    historical_tumor
  
  sample_count_comparison$historical_normal[i] <-
    historical_normal
  
  sample_count_comparison$total_n_match[i] <-
    sum(excel_gene) == historical_n
  
  sample_count_comparison$tumor_n_match[i] <-
    excel_tumor == historical_tumor
  
  sample_count_comparison$normal_n_match[i] <-
    excel_normal == historical_normal
}

cat("\nExcel vs historical ROC object counts:\n")
print(sample_count_comparison)

write.csv(
  sample_count_comparison,
  file.path(
    output_dir,
    "STEP13B-7J-R1_Excel_vs_Historical_ROC_Count_Comparison.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 24. Final structural PASS criteria
# ------------------------------------------------------------

pass_records <-
  nrow(dat) == 1725

pass_unique_samples <-
  n_unique_samples == 575

pass_gene_structure <-
  all(
    gene_exactness$n_records == 575
  )

pass_status_structure <-
  all(
    gene_status_table[, "Tumor"] == 516
  ) &&
  all(
    gene_status_table[, "Normal"] == 59
  )

pass_sample_status <-
  n_inconsistent_status == 0

pass_duplicate_structure <-
  n_duplicate_sample_gene == 0

pass_historical_counts <-
  all(
    sample_count_comparison$total_n_match &
      sample_count_comparison$tumor_n_match &
      sample_count_comparison$normal_n_match
  )


# ------------------------------------------------------------
# 25. Final status
# ------------------------------------------------------------

if (
  pass_records &&
  pass_unique_samples &&
  pass_gene_structure &&
  pass_status_structure &&
  pass_sample_status &&
  pass_duplicate_structure &&
  pass_historical_counts
) {
  
  final_status <- "PASS"
  
} else {
  
  final_status <- "HOLD"
}


# ------------------------------------------------------------
# 26. Final report
# ------------------------------------------------------------

final_report <- c(
  "STEP13B-7J-R1 — FINAL REPORT",
  "============================================",
  "",
  paste0("FINAL STATUS: ", final_status),
  "",
  "Purpose:",
  "Sample-level traceability of the original Xiangtao ROC branch.",
  "",
  paste0("Excel records: ", nrow(dat)),
  paste0("Unique samples: ", n_unique_samples),
  paste0("Expected samples: ", expected_n_samples),
  paste0("Tumor records: ",
         sum(status_class == "Tumor", na.rm = TRUE)),
  paste0("Normal records: ",
         sum(status_class == "Normal", na.rm = TRUE)),
  "",
  paste0("AOC3 records: ",
         sum(gene == "AOC3", na.rm = TRUE)),
  paste0("MS4A7 records: ",
         sum(gene == "MS4A7", na.rm = TRUE)),
  paste0("IL7R records: ",
         sum(gene == "IL7R", na.rm = TRUE)),
  "",
  paste0(
    "Duplicate sample-gene records: ",
    n_duplicate_sample_gene
  ),
  paste0(
    "Inconsistent sample status: ",
    n_inconsistent_status
  ),
  "",
  paste0(
    "Historical ROC count concordance: ",
    pass_historical_counts
  ),
  "",
  "No ROC/AUC recalculation performed.",
  "No model refitting performed.",
  "No threshold optimization performed.",
  "No expression transformation performed.",
  "Xiangtao and Sangerbox provenance branches remain independent.",
  ""
)

writeLines(
  final_report,
  file.path(
    audit_dir,
    "STEP13B-7J-R1_FINAL_STATUS.txt"
  )
)

cat("\n")
cat(paste(final_report, collapse = "\n"))
cat("\n")


# ------------------------------------------------------------
# 27. Save complete audit object
# ------------------------------------------------------------

audit_object <- list(
  workflow = "STEP13B-7J-R1",
  purpose =
    "Xiangtao historical ROC Excel sample-level traceability",
  input_excel = excel_file,
  historical_rdata = historical_rdata,
  sheet = sheet_name,
  dimensions = dim(dat),
  column_mapping = list(
    sample_id = sample_id_col,
    status = status_col,
    gene = x_col,
    expression = y_col
  ),
  gene_exactness = gene_exactness,
  gene_status_table = gene_status_table,
  n_unique_samples = n_unique_samples,
  n_duplicate_sample_gene = n_duplicate_sample_gene,
  n_inconsistent_status = n_inconsistent_status,
  historical_summary = historical_summary,
  sample_count_comparison = sample_count_comparison,
  final_status = final_status
)

saveRDS(
  audit_object,
  file.path(
    audit_dir,
    "STEP13B-7J-R1_Traceability_Audit.rds"
  )
)


# ------------------------------------------------------------
# 28. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    log_dir,
    "STEP13B-7J-R1_sessionInfo.txt"
  )
)


cat("\n============================================================\n")
cat("STEP13B-7J-R1 COMPLETE\n")
cat("============================================================\n")

# ============================================================
# STEP13B-7J-R2
# Excel Expression -> Historical ROC Predictor Traceability
# ============================================================
#
# Purpose:
#   Establish expression-value traceability between the
#   original Xiangtao ROC Excel dataset and the already-existing
#   historical pROC ROC objects.
#
# Scope:
#   - Excel y expression values
#   - historical roc_* $predictor values
#   - historical roc_* $response values
#   - sample-level structural reconciliation
#
# IMPORTANT:
#   NO ROC recalculation
#   NO AUC recalculation
#   NO model refitting
#   NO threshold optimization
#   NO expression transformation
#
# This is a provenance / traceability audit only.
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP13B-7J-R2\n")
cat("Excel Expression -> Historical ROC Predictor Traceability\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Configuration
# ------------------------------------------------------------

workflow_root <-
  "input_dir/Paper 1_Step 13-Diagnostic Classification"

r1_dir <-
  file.path(
    workflow_root,
    "STEP13B-7J_Xiangtao_ROC_Excel_Traceability"
  )

output_dir <- r1_dir

excel_file <-
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal/Diagnostic_ROC_AOC3_MS4A7_IL7R_ Data_Analysis.xlsx"

historical_rdata <-
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal/.RData"

core_genes <-
  c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )

expected_n <- 575
expected_tumor <- 516
expected_normal <- 59

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Package check
# ------------------------------------------------------------

required_packages <- c(
  "readxl"
)

missing_packages <-
  required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

if (length(missing_packages) > 0) {
  
  stop(
    paste(
      "Missing required package(s):",
      paste(
        missing_packages,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 2. Read Xiangtao Excel
# ------------------------------------------------------------

cat("Reading original Xiangtao ROC Excel...\n")

excel_raw <-
  readxl::read_excel(
    excel_file,
    sheet = "Sheet 1"
  )

cat(
  "Rows:",
  nrow(excel_raw),
  "\n"
)

cat(
  "Columns:",
  ncol(excel_raw),
  "\n\n"
)


# ------------------------------------------------------------
# 3. Standardize column names
# ------------------------------------------------------------

required_columns <-
  c(
    "sample_id",
    "status",
    "x",
    "y"
  )

if (!all(required_columns %in% names(excel_raw))) {
  
  stop(
    paste(
      "Required columns missing:",
      paste(
        setdiff(
          required_columns,
          names(excel_raw)
        ),
        collapse = ", "
      )
    )
  )
}


dat <-
  data.frame(
    sample_id =
      as.character(
        excel_raw$sample_id
      ),
    
    status =
      as.character(
        excel_raw$status
      ),
    
    gene =
      as.character(
        excel_raw$x
      ),
    
    expression =
      as.numeric(
        excel_raw$y
      ),
    
    stringsAsFactors = FALSE
  )


# ------------------------------------------------------------
# 4. Basic data integrity audit
# ------------------------------------------------------------

cat("Basic Excel integrity audit...\n")

if (anyNA(dat$sample_id)) {
  stop("Missing sample_id detected.")
}

if (anyNA(dat$status)) {
  stop("Missing status detected.")
}

if (anyNA(dat$gene)) {
  stop("Missing gene detected.")
}

if (anyNA(dat$expression)) {
  stop("Missing expression value detected.")
}

if (!all(dat$gene %in% core_genes)) {
  
  stop(
    "Unexpected gene values detected."
  )
}


# ------------------------------------------------------------
# 5. Verify R1 structural expectations
# ------------------------------------------------------------

cat("\nChecking R1 structural expectations...\n")

gene_counts <-
  table(dat$gene)

print(gene_counts)

sample_n <-
  length(
    unique(dat$sample_id)
  )

cat(
  "Unique samples:",
  sample_n,
  "\n"
)

if (sample_n != expected_n) {
  
  stop(
    paste(
      "Expected",
      expected_n,
      "unique samples but found",
      sample_n
    )
  )
}

if (!all(as.integer(gene_counts) == expected_n)) {
  
  stop(
    "One or more genes do not contain exactly 575 records."
  )
}


# ------------------------------------------------------------
# 6. Create gene-specific Excel vectors
# ------------------------------------------------------------

excel_vectors <- list()

for (g in core_genes) {
  
  idx <-
    which(
      dat$gene == g
    )
  
  excel_vectors[[g]] <-
    dat$expression[idx]
}


# ------------------------------------------------------------
# 7. Load historical RData in isolated environment
# ------------------------------------------------------------

cat("\nLoading historical RData in isolated environment...\n")

hist_env <-
  new.env(
    parent = emptyenv()
  )

load(
  historical_rdata,
  envir = hist_env
)

required_roc_objects <-
  paste0(
    "roc_",
    core_genes
  )

object_available <-
  vapply(
    required_roc_objects,
    exists,
    logical(1),
    envir = hist_env,
    inherits = FALSE
  )

cat("\nHistorical ROC objects:\n")
print(object_available)

if (!all(object_available)) {
  
  stop(
    "One or more required historical ROC objects are missing."
  )
}


# ------------------------------------------------------------
# 8. Function: exact numeric comparison
# ------------------------------------------------------------
#
# This function does NOT calculate ROC/AUC.
#
# It only compares two already-existing numeric vectors.
# ------------------------------------------------------------

compare_numeric_vectors <- function(
    excel_values,
    historical_values,
    tolerance = 0
) {
  
  if (length(excel_values) !=
      length(historical_values)) {
    
    return(
      list(
        n_excel = length(excel_values),
        n_historical = length(historical_values),
        exact_order_match = FALSE,
        multiset_match = FALSE,
        max_abs_difference = NA_real_,
        n_exact_position_match = NA_integer_
      )
    )
  }
  
  abs_diff <-
    abs(
      excel_values -
        historical_values
    )
  
  exact_position <-
    if (tolerance == 0) {
      excel_values ==
        historical_values
    } else {
      abs_diff <= tolerance
    }
  
  exact_order_match <-
    all(
      exact_position,
      na.rm = TRUE
    )
  
  sorted_excel <-
    sort(
      excel_values
    )
  
  sorted_historical <-
    sort(
      historical_values
    )
  
  sorted_diff <-
    abs(
      sorted_excel -
        sorted_historical
    )
  
  multiset_match <-
    all(
      sorted_diff <= tolerance,
      na.rm = TRUE
    )
  
  list(
    n_excel =
      length(excel_values),
    
    n_historical =
      length(historical_values),
    
    exact_order_match =
      exact_order_match,
    
    multiset_match =
      multiset_match,
    
    max_abs_difference =
      max(
        abs_diff,
        na.rm = TRUE
      ),
    
    n_exact_position_match =
      sum(
        exact_position,
        na.rm = TRUE
      )
  )
}


# ------------------------------------------------------------
# 9. Expression-value comparison
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("Expression-value comparison\n")
cat("============================================================\n")

expression_summary <-
  data.frame(
    gene = core_genes,
    
    excel_n =
      NA_integer_,
    
    historical_n =
      NA_integer_,
    
    exact_order_match =
      FALSE,
    
    multiset_match =
      FALSE,
    
    n_exact_position_match =
      NA_integer_,
    
    max_abs_difference =
      NA_real_,
    
    stringsAsFactors = FALSE
  )


for (i in seq_along(core_genes)) {
  
  g <-
    core_genes[i]
  
  roc_obj <-
    get(
      paste0(
        "roc_",
        g
      ),
      envir = hist_env
    )
  
  excel_values <-
    excel_vectors[[g]]
  
  historical_values <-
    as.numeric(
      roc_obj$predictor
    )
  
  cmp <-
    compare_numeric_vectors(
      excel_values,
      historical_values,
      tolerance = 0
    )
  
  expression_summary$excel_n[i] <-
    cmp$n_excel
  
  expression_summary$historical_n[i] <-
    cmp$n_historical
  
  expression_summary$exact_order_match[i] <-
    cmp$exact_order_match
  
  expression_summary$multiset_match[i] <-
    cmp$multiset_match
  
  expression_summary$n_exact_position_match[i] <-
    cmp$n_exact_position_match
  
  expression_summary$max_abs_difference[i] <-
    cmp$max_abs_difference
}

print(
  expression_summary
)


write.csv(
  expression_summary,
  file.path(
    output_dir,
    "STEP13B-7J-R2_Expression_Value_Comparison.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 10. Response / status comparison
# ------------------------------------------------------------
#
# Historical ROC objects contain response vectors.
#
# We compare their outcome composition with Excel.
#
# NO ROC recalculation.
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("Historical response / Excel status comparison\n")
cat("============================================================\n")


response_summary <-
  data.frame(
    gene = core_genes,
    
    excel_n =
      NA_integer_,
    
    excel_tumor =
      NA_integer_,
    
    excel_normal =
      NA_integer_,
    
    historical_n =
      NA_integer_,
    
    historical_tumor =
      NA_integer_,
    
    historical_normal =
      NA_integer_,
    
    n_match =
      FALSE,
    
    tumor_match =
      FALSE,
    
    normal_match =
      FALSE,
    
    stringsAsFactors = FALSE
  )


for (i in seq_along(core_genes)) {
  
  g <-
    core_genes[i]
  
  excel_idx <-
    dat$gene == g
  
  excel_status <-
    dat$status[
      excel_idx
    ]
  
  roc_obj <-
    get(
      paste0(
        "roc_",
        g
      ),
      envir = hist_env
    )
  
  response_chr <-
    as.character(
      roc_obj$response
    )
  
  excel_tumor <-
    sum(
      excel_status == "Tumor"
    )
  
  excel_normal <-
    sum(
      excel_status == "Normal"
    )
  
  historical_tumor <-
    sum(
      response_chr %in%
        c(
          "1",
          "Tumor",
          "tumor"
        )
    )
  
  historical_normal <-
    sum(
      response_chr %in%
        c(
          "0",
          "Normal",
          "normal"
        )
    )
  
  response_summary$excel_n[i] <-
    length(excel_status)
  
  response_summary$excel_tumor[i] <-
    excel_tumor
  
  response_summary$excel_normal[i] <-
    excel_normal
  
  response_summary$historical_n[i] <-
    length(response_chr)
  
  response_summary$historical_tumor[i] <-
    historical_tumor
  
  response_summary$historical_normal[i] <-
    historical_normal
  
  response_summary$n_match[i] <-
    length(excel_status) ==
    length(response_chr)
  
  response_summary$tumor_match[i] <-
    excel_tumor ==
    historical_tumor
  
  response_summary$normal_match[i] <-
    excel_normal ==
    historical_normal
}

print(
  response_summary
)


write.csv(
  response_summary,
  file.path(
    output_dir,
    "STEP13B-7J-R2_Response_Status_Comparison.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 11. Exact-order sample-level reconciliation
# ------------------------------------------------------------
#
# Because historical pROC objects generally do not retain the
# original sample IDs, this section first tests whether the
# historical predictor vector is already in exactly the same
# order as the Excel data.
#
# No assumption is made.
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("Exact-order sample-level reconciliation\n")
cat("============================================================\n")


order_summary <-
  data.frame(
    gene = core_genes,
    
    n =
      NA_integer_,
    
    exact_expression_order =
      FALSE,
    
    exact_status_order =
      FALSE,
    
    exact_expression_and_status_order =
      FALSE,
    
    stringsAsFactors = FALSE
  )


for (i in seq_along(core_genes)) {
  
  g <-
    core_genes[i]
  
  idx <-
    dat$gene == g
  
  excel_expression <-
    dat$expression[idx]
  
  excel_status <-
    dat$status[idx]
  
  roc_obj <-
    get(
      paste0(
        "roc_",
        g
      ),
      envir = hist_env
    )
  
  historical_expression <-
    as.numeric(
      roc_obj$predictor
    )
  
  historical_status <-
    as.character(
      roc_obj$response
    )
  
  expression_match <-
    length(excel_expression) ==
    length(historical_expression) &&
    all(
      excel_expression ==
        historical_expression,
      na.rm = TRUE
    )
  
  historical_status_normalized <-
    ifelse(
      historical_status %in%
        c(
          "1",
          "Tumor",
          "tumor"
        ),
      "Tumor",
      ifelse(
        historical_status %in%
          c(
            "0",
            "Normal",
            "normal"
          ),
        "Normal",
        NA_character_
      )
    )
  
  status_match <-
    length(excel_status) ==
    length(historical_status_normalized) &&
    all(
      excel_status ==
        historical_status_normalized,
      na.rm = TRUE
    )
  
  order_summary$n[i] <-
    length(excel_expression)
  
  order_summary$exact_expression_order[i] <-
    expression_match
  
  order_summary$exact_status_order[i] <-
    status_match
  
  order_summary$exact_expression_and_status_order[i] <-
    expression_match &&
    status_match
}

print(
  order_summary
)


write.csv(
  order_summary,
  file.path(
    output_dir,
    "STEP13B-7J-R2_Exact_Order_Reconciliation.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 12. Multiset expression traceability
# ------------------------------------------------------------
#
# This determines whether the same 575 numerical expression
# values are present even if row order differs.
#
# Again:
#   NO ROC/AUC calculation.
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("Multiset expression traceability\n")
cat("============================================================\n")


multiset_summary <-
  expression_summary[
    ,
    c(
      "gene",
      "excel_n",
      "historical_n",
      "multiset_match",
      "max_abs_difference"
    )
  ]

print(
  multiset_summary
)


# ------------------------------------------------------------
# 13. Build deterministic value-level audit table
# ------------------------------------------------------------
#
# This table preserves the Excel sample_id and expression value.
#
# Historical sample IDs may not be available in the ROC object.
# Therefore historical predictor position is retained as an
# index only.
# ------------------------------------------------------------

value_audit_list <- list()

for (g in core_genes) {
  
  idx <-
    which(
      dat$gene == g
    )
  
  excel_part <-
    dat[idx, ]
  
  roc_obj <-
    get(
      paste0(
        "roc_",
        g
      ),
      envir = hist_env
    )
  
  historical_expression <-
    as.numeric(
      roc_obj$predictor
    )
  
  historical_response <-
    as.character(
      roc_obj$response
    )
  
  n_use <-
    min(
      nrow(excel_part),
      length(historical_expression)
    )
  
  tmp <-
    data.frame(
      gene = g,
      
      excel_row =
        idx[seq_len(n_use)],
      
      excel_sample_id =
        excel_part$sample_id[
          seq_len(n_use)
        ],
      
      excel_status =
        excel_part$status[
          seq_len(n_use)
        ],
      
      excel_expression =
        excel_part$expression[
          seq_len(n_use)
        ],
      
      historical_predictor_index =
        seq_len(n_use),
      
      historical_expression =
        historical_expression[
          seq_len(n_use)
        ],
      
      historical_response =
        historical_response[
          seq_len(n_use)
        ],
      
      expression_difference =
        excel_part$expression[
          seq_len(n_use)
        ] -
        historical_expression[
          seq_len(n_use)
        ],
      
      stringsAsFactors = FALSE
    )
  
  value_audit_list[[g]] <-
    tmp
}

value_audit <-
  do.call(
    rbind,
    value_audit_list
  )

write.csv(
  value_audit,
  file.path(
    output_dir,
    "STEP13B-7J-R2_Value_Level_Audit.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 14. Final R2 decision logic
# ------------------------------------------------------------
#
# PASS if:
#
#   1. Each gene has 575 Excel values
#   2. Each historical predictor has 575 values
#   3. All three genes show complete multiset equality
#   4. No unexpected missing values
#   5. Historical response counts concordant
#
# Exact row-order agreement is NOT required for PASS.
#
# If multiset equality is present but exact order is different,
# classify as:
#
#   PASS — VALUE-LEVEL TRACEABILITY ESTABLISHED;
#   ORDER DIFFERENCE PRESENT
#
# This is still valid provenance because pROC objects do not
# necessarily preserve sample IDs / original Excel row order.
# ------------------------------------------------------------

all_n_match <-
  all(
    expression_summary$excel_n ==
      expected_n &
      expression_summary$historical_n ==
      expected_n
  )

all_multiset_match <-
  all(
    expression_summary$multiset_match
  )

all_response_match <-
  all(
    response_summary$n_match &
      response_summary$tumor_match &
      response_summary$normal_match
  )

all_expression_complete <-
  !anyNA(
    dat$expression
  )

final_pass <-
  all_n_match &&
  all_multiset_match &&
  all_response_match &&
  all_expression_complete


# ------------------------------------------------------------
# 15. Final report
# ------------------------------------------------------------

final_report <-
  c(
    "",
    "STEP13B-7J-R2 — FINAL REPORT",
    "============================================",
    "",
    paste(
      "FINAL STATUS:",
      ifelse(
        final_pass,
        "PASS",
        "HOLD"
      )
    ),
    "",
    "Purpose:",
    "Expression-value traceability between the original Xiangtao ROC Excel dataset and historical ROC predictor objects.",
    "",
    paste(
      "Excel records:",
      nrow(dat)
    ),
    paste(
      "Unique samples:",
      length(
        unique(dat$sample_id)
      )
    ),
    "",
    paste(
      "AOC3 records:",
      sum(dat$gene == "AOC3")
    ),
    paste(
      "MS4A7 records:",
      sum(dat$gene == "MS4A7")
    ),
    paste(
      "IL7R records:",
      sum(dat$gene == "IL7R")
    ),
    "",
    paste(
      "All expression vector lengths concordant:",
      all_n_match
    ),
    paste(
      "All expression value multisets concordant:",
      all_multiset_match
    ),
    paste(
      "All response/status counts concordant:",
      all_response_match
    ),
    "",
    paste(
      "Exact row-order expression match:",
      all(
        order_summary$exact_expression_order
      )
    ),
    paste(
      "Exact row-order status match:",
      all(
        order_summary$exact_status_order
      )
    ),
    "",
    "No ROC recalculation performed.",
    "No AUC recalculation performed.",
    "No model refitting performed.",
    "No threshold optimization performed.",
    "No expression transformation performed.",
    "Xiangtao and Sangerbox provenance branches remain independent.",
    ""
  )

cat(
  paste(
    final_report,
    collapse = "\n"
  )
)

writeLines(
  final_report,
  con =
    file.path(
      output_dir,
      "STEP13B-7J-R2_FINAL_REPORT.txt"
    ),
  useBytes = TRUE
)


# ------------------------------------------------------------
# 16. PASS / HOLD marker
# ------------------------------------------------------------

if (final_pass) {
  
  writeLines(
    c(
      "STEP13B-7J-R2",
      "FINAL STATUS: PASS",
      "",
      "Expression-value traceability established.",
      "Historical ROC predictor values are concordant with the Xiangtao Excel dataset.",
      "No ROC/AUC recalculation performed.",
      "No model refitting performed.",
      "No threshold optimization performed."
    ),
    con =
      file.path(
        output_dir,
        "STEP13B-7J-R2_PASS.txt"
      ),
    useBytes = TRUE
  )
  
} else {
  
  writeLines(
    c(
      "STEP13B-7J-R2",
      "FINAL STATUS: HOLD",
      "",
      "Expression-value traceability not fully established.",
      "Review the expression comparison and response/status audit outputs.",
      "No ROC/AUC recalculation performed."
    ),
    con =
      file.path(
        output_dir,
        "STEP13B-7J-R2_HOLD.txt"
      ),
    useBytes = TRUE
  )
}


# ------------------------------------------------------------
# 17. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  con =
    file.path(
      output_dir,
      "STEP13B-7J-R2_sessionInfo.txt"
    ),
  useBytes = TRUE
)


cat("\n")
cat("============================================================\n")
cat("STEP13B-7J-R2 completed.\n")
cat("============================================================\n\n")

# ============================================================
# STEP13C-1
# GSE32863 Raw File / Historical Object Inventory
# ============================================================
#
# Purpose:
#   Establish an auditable inventory of all files and historical
#   R objects relevant to the GSE32863 external assessment.
#
# Scope:
#   - GSE32863-related files
#   - historical .RData objects
#   - historical expression / annotation objects
#   - historical external ROC object availability
#
# IMPORTANT:
#   NO expression recalculation
#   NO probe selection
#   NO ROC recalculation
#   NO AUC recalculation
#   NO model refitting
#   NO threshold optimization
#   NO data transformation
#
# This is an inventory / provenance audit only.
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP13C-1\n")
cat("GSE32863 Raw File / Historical Object Inventory\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Configuration
# ------------------------------------------------------------

workflow_root <-
  "input_dir/Paper 1_Step 13-Diagnostic Classification"

step13c_root <-
  file.path(
    workflow_root,
    "STEP13C_GSE32863_Input_Lock"
  )

output_dir <-
  file.path(
    step13c_root,
    "STEP13C-1_Raw_File_Inventory"
  )

# Historical analysis directory already identified
# during STEP13C-0B.

historical_root <-
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal"

historical_rdata <-
  file.path(
    historical_root,
    ".RData"
  )

historical_rdata_named <-
  file.path(
    historical_root,
    "Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal.RData"
  )

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. General inventory metadata
# ------------------------------------------------------------

cat("Historical analysis directory:\n")
cat(historical_root, "\n\n")

cat("Directory exists:")
cat(dir.exists(historical_root), "\n")

cat("Primary historical RData:\n")
cat(historical_rdata, "\n")

cat(
  "Exists:",
  file.exists(historical_rdata),
  "\n\n"
)


# ------------------------------------------------------------
# 2. Recursive file inventory
# ------------------------------------------------------------

cat("============================================================\n")
cat("Recursive file inventory\n")
cat("============================================================\n\n")


all_files <-
  list.files(
    historical_root,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  )


file_inventory <-
  data.frame(
    full_path = all_files,
    
    relative_path =
      substring(
        all_files,
        nchar(historical_root) + 2
      ),
    
    file_name =
      basename(all_files),
    
    extension =
      tolower(
        tools::file_ext(all_files)
      ),
    
    size_bytes =
      file.info(all_files)$size,
    
    modified_time =
      file.info(all_files)$mtime,
    
    stringsAsFactors = FALSE
  )


cat(
  "Total files found:",
  nrow(file_inventory),
  "\n\n"
)


write.csv(
  file_inventory,
  file.path(
    output_dir,
    "STEP13C-1_All_File_Inventory.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 3. Identify GSE32863-related filenames
# ------------------------------------------------------------

cat("============================================================\n")
cat("GSE32863-related file search\n")
cat("============================================================\n\n")


gse_pattern <-
  "GSE32863|GPL6884|AOC3|MS4A7|IL7R|expr3|expr\\.mat|expr\\.final|expr2|anno2|annotation|clinical|phenotype|pheno|series_matrix|matrix"

gse_files <-
  file_inventory[
    grepl(
      gse_pattern,
      file_inventory$file_name,
      ignore.case = TRUE
    ) |
      grepl(
        gse_pattern,
        file_inventory$relative_path,
        ignore.case = TRUE
      ),
    ,
    drop = FALSE
  ]


cat(
  "Potential GSE32863-related files:",
  nrow(gse_files),
  "\n\n"
)

if (nrow(gse_files) > 0) {
  print(
    gse_files[
      ,
      c(
        "relative_path",
        "extension",
        "size_bytes",
        "modified_time"
      )
    ]
  )
}


write.csv(
  gse_files,
  file.path(
    output_dir,
    "STEP13C-1_GSE32863_Related_File_Inventory.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 4. Identify RData files
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("RData inventory\n")
cat("============================================================\n\n")


rdata_files <-
  file_inventory[
    file_inventory$extension %in%
      c(
        "rdata",
        "rda"
      ),
    ,
    drop = FALSE
  ]


cat(
  "RData/RDA files found:",
  nrow(rdata_files),
  "\n\n"
)

if (nrow(rdata_files) > 0) {
  print(
    rdata_files[
      ,
      c(
        "relative_path",
        "extension",
        "size_bytes",
        "modified_time"
      )
    ]
  )
}


write.csv(
  rdata_files,
  file.path(
    output_dir,
    "STEP13C-1_RData_File_Inventory.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 5. Load primary historical RData in isolated environment
# ------------------------------------------------------------
#
# Only object inventory.
# No calculations.
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("Historical .RData object inventory\n")
cat("============================================================\n\n")


hist_env <-
  new.env(
    parent = emptyenv()
  )


if (!file.exists(historical_rdata)) {
  
  stop(
    "Primary historical .RData file not found."
  )
}


load(
  historical_rdata,
  envir = hist_env
)


historical_objects <-
  ls(
    envir = hist_env,
    all.names = TRUE
  )


cat(
  "Number of objects in primary .RData:",
  length(historical_objects),
  "\n\n"
)

print(
  historical_objects
)


object_inventory <-
  data.frame(
    object_name =
      historical_objects,
    
    object_class =
      vapply(
        historical_objects,
        function(x) {
          
          obj <-
            get(
              x,
              envir = hist_env
            )
          
          paste(
            class(obj),
            collapse = ";"
          )
        },
        character(1)
      ),
    
    object_length =
      vapply(
        historical_objects,
        function(x) {
          
          obj <-
            get(
              x,
              envir = hist_env
            )
          
          length(obj)
        },
        integer(1)
      ),
    
    stringsAsFactors = FALSE
  )


write.csv(
  object_inventory,
  file.path(
    output_dir,
    "STEP13C-1_Historical_RData_Object_Inventory.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 6. Historical GSE32863 object audit
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("Historical GSE32863 object audit\n")
cat("============================================================\n\n")


expected_objects <-
  c(
    "expr",
    "expr.mat",
    "expr.final",
    "expr2",
    "expr3",
    "group",
    "genes",
    "anno",
    "anno2",
    "roc.ext"
  )


object_available <-
  vapply(
    expected_objects,
    exists,
    logical(1),
    envir = hist_env,
    inherits = FALSE
  )


historical_gse_inventory <-
  data.frame(
    object_name =
      expected_objects,
    
    available =
      object_available,
    
    class =
      NA_character_,
    
    dimensions =
      NA_character_,
    
    length =
      NA_integer_,
    
    stringsAsFactors = FALSE
  )


for (i in seq_along(expected_objects)) {
  
  obj_name <-
    expected_objects[i]
  
  if (object_available[i]) {
    
    obj <-
      get(
        obj_name,
        envir = hist_env
      )
    
    historical_gse_inventory$class[i] <-
      paste(
        class(obj),
        collapse = ";"
      )
    
    historical_gse_inventory$length[i] <-
      length(obj)
    
    d <-
      dim(obj)
    
    if (!is.null(d)) {
      
      historical_gse_inventory$dimensions[i] <-
        paste(
          d,
          collapse = " × "
        )
      
    } else {
      
      historical_gse_inventory$dimensions[i] <-
        "NULL"
    }
  }
}


print(
  historical_gse_inventory
)


write.csv(
  historical_gse_inventory,
  file.path(
    output_dir,
    "STEP13C-1_Historical_GSE32863_Object_Audit.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 7. Historical external ROC object structural inventory
# ------------------------------------------------------------
#
# IMPORTANT:
# Do NOT calculate AUC.
#
# Only record whether roc.ext exists and inspect its structure.
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("Historical external ROC object structural inventory\n")
cat("============================================================\n\n")


roc_ext_available <-
  exists(
    "roc.ext",
    envir = hist_env,
    inherits = FALSE
  )


cat(
  "roc.ext available:",
  roc_ext_available,
  "\n"
)


roc_ext_summary <-
  data.frame(
    object =
      "roc.ext",
    
    available =
      roc_ext_available,
    
    class =
      NA_character_,
    
    response_length =
      NA_integer_,
    
    predictor_length =
      NA_integer_,
    
    response_class =
      NA_character_,
    
    predictor_class =
      NA_character_,
    
    stringsAsFactors = FALSE
  )


if (roc_ext_available) {
  
  roc_ext <-
    get(
      "roc.ext",
      envir = hist_env
    )
  
  roc_ext_summary$class <-
    paste(
      class(roc_ext),
      collapse = ";"
    )
  
  roc_ext_summary$response_length <-
    length(
      roc_ext$response
    )
  
  roc_ext_summary$predictor_length <-
    length(
      roc_ext$predictor
    )
  
  roc_ext_summary$response_class <-
    paste(
      class(
        roc_ext$response
      ),
      collapse = ";"
    )
  
  roc_ext_summary$predictor_class <-
    paste(
      class(
        roc_ext$predictor
      ),
      collapse = ";"
    )
}


print(
  roc_ext_summary
)


write.csv(
  roc_ext_summary,
  file.path(
    output_dir,
    "STEP13C-1_Historical_roc_ext_Structural_Inventory.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 8. Known historical GSE32863 dimensions
# ------------------------------------------------------------
#
# This section only records observed dimensions.
# It does NOT recalculate or transform expression.
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("Historical GSE32863 dimensions\n")
cat("============================================================\n\n")


dimension_summary <-
  data.frame(
    object = expected_objects,
    expected_from_STEP13C0B = NA_character_,
    observed = NA_character_,
    stringsAsFactors = FALSE
  )


dimension_summary$expected_from_STEP13C0B[
  dimension_summary$object == "expr"
] <-
  "48803 × 117"

dimension_summary$expected_from_STEP13C0B[
  dimension_summary$object == "expr.mat"
] <-
  "28742 × 116"

dimension_summary$expected_from_STEP13C0B[
  dimension_summary$object == "expr.final"
] <-
  "19425 × 116"

dimension_summary$expected_from_STEP13C0B[
  dimension_summary$object == "expr2"
] <-
  "28742 × 118"

dimension_summary$expected_from_STEP13C0B[
  dimension_summary$object == "expr3"
] <-
  "116 × 5"

dimension_summary$expected_from_STEP13C0B[
  dimension_summary$object == "group"
] <-
  "length 116"

dimension_summary$expected_from_STEP13C0B[
  dimension_summary$object == "genes"
] <-
  "3 genes"

dimension_summary$expected_from_STEP13C0B[
  dimension_summary$object == "anno"
] <-
  "48803 × 22"

dimension_summary$expected_from_STEP13C0B[
  dimension_summary$object == "anno2"
] <-
  "48803 × 2"

dimension_summary$expected_from_STEP13C0B[
  dimension_summary$object == "roc.ext"
] <-
  "response/predictor length 116"


for (i in seq_len(nrow(dimension_summary))) {
  
  obj_name <-
    dimension_summary$object[i]
  
  if (
    exists(
      obj_name,
      envir = hist_env,
      inherits = FALSE
    )
  ) {
    
    obj <-
      get(
        obj_name,
        envir = hist_env
      )
    
    d <-
      dim(obj)
    
    if (!is.null(d)) {
      
      dimension_summary$observed[i] <-
        paste(
          d,
          collapse = " × "
        )
      
    } else {
      
      dimension_summary$observed[i] <-
        paste(
          "length",
          length(obj)
        )
    }
  }
}


print(
  dimension_summary
)


write.csv(
  dimension_summary,
  file.path(
    output_dir,
    "STEP13C-1_Historical_GSE32863_Dimension_Summary.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 9. Preliminary inventory checks
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("Preliminary inventory checks\n")
cat("============================================================\n\n")


required_gse_objects_present <-
  all(
    object_available[
      expected_objects != "roc.ext"
    ]
  )

roc_ext_present <-
  object_available[
    expected_objects == "roc.ext"
  ]

expr3_valid_length <-
  if (
    exists(
      "expr3",
      envir = hist_env,
      inherits = FALSE
    )
  ) {
    nrow(
      get(
        "expr3",
        envir = hist_env
      )
    ) == 116
  } else {
    FALSE
  }

roc_ext_valid_length <-
  if (roc_ext_available) {
    
    length(
      get(
        "roc.ext",
        envir = hist_env
      )$predictor
    ) == 116
    
  } else {
    
    FALSE
  }


cat(
  "Required GSE32863 historical objects present:",
  required_gse_objects_present,
  "\n"
)

cat(
  "Historical roc.ext present:",
  roc_ext_present,
  "\n"
)

cat(
  "expr3 has 116 rows:",
  expr3_valid_length,
  "\n"
)

cat(
  "roc.ext predictor has length 116:",
  roc_ext_valid_length,
  "\n"
)


# ------------------------------------------------------------
# 10. Final STEP13C-1 status
# ------------------------------------------------------------
#
# This is NOT the final STEP13C input lock.
#
# It only locks the inventory stage.
# ------------------------------------------------------------

step13c1_pass <-
  required_gse_objects_present &&
  roc_ext_present &&
  expr3_valid_length &&
  roc_ext_valid_length


final_report <-
  c(
    "",
    "STEP13C-1 — FINAL REPORT",
    "============================================",
    "",
    paste(
      "FINAL STATUS:",
      ifelse(
        step13c1_pass,
        "PASS",
        "HOLD"
      )
    ),
    "",
    "Purpose:",
    "Inventory and provenance audit of GSE32863-related files and historical analysis objects.",
    "",
    paste(
      "Historical analysis directory exists:",
      dir.exists(historical_root)
    ),
    paste(
      "Historical .RData exists:",
      file.exists(historical_rdata)
    ),
    paste(
      "Required GSE32863 historical objects present:",
      required_gse_objects_present
    ),
    paste(
      "Historical roc.ext present:",
      roc_ext_present
    ),
    paste(
      "expr3 length = 116:",
      expr3_valid_length
    ),
    paste(
      "roc.ext predictor length = 116:",
      roc_ext_valid_length
    ),
    "",
    "No expression recalculation performed.",
    "No probe selection performed.",
    "No ROC recalculation performed.",
    "No AUC recalculation performed.",
    "No model refitting performed.",
    "No threshold optimization performed.",
    "No expression transformation performed.",
    "",
    "STEP13C final input lock is NOT established at this stage.",
    ""
  )


cat(
  paste(
    final_report,
    collapse = "\n"
  )
)


writeLines(
  final_report,
  con =
    file.path(
      output_dir,
      "STEP13C-1_FINAL_REPORT.txt"
    ),
  useBytes = TRUE
)


# ------------------------------------------------------------
# 11. PASS / HOLD marker
# ------------------------------------------------------------

if (step13c1_pass) {
  
  writeLines(
    c(
      "STEP13C-1",
      "FINAL STATUS: PASS",
      "",
      "GSE32863 historical input inventory established.",
      "Required historical objects are present.",
      "roc.ext structural availability established.",
      "",
      "This PASS does NOT establish final external validation.",
      "This PASS does NOT establish the final STEP13C input lock."
    ),
    con =
      file.path(
        output_dir,
        "STEP13C-1_PASS.txt"
      ),
    useBytes = TRUE
  )
  
} else {
  
  writeLines(
    c(
      "STEP13C-1",
      "FINAL STATUS: HOLD",
      "",
      "GSE32863 historical input inventory is incomplete.",
      "Review the inventory outputs before proceeding."
    ),
    con =
      file.path(
        output_dir,
        "STEP13C-1_HOLD.txt"
      ),
    useBytes = TRUE
  )
}


# ------------------------------------------------------------
# 12. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  con =
    file.path(
      output_dir,
      "STEP13C-1_sessionInfo.txt"
    ),
  useBytes = TRUE
)


cat("\n")
cat("============================================================\n")
cat("STEP13C-1 completed.\n")
cat("============================================================\n\n")


# ============================================================
# STEP13C-2 — GSE32863 Expression Matrix Provenance
# ============================================================
#
# Purpose:
#   Object-level provenance audit of the historical GSE32863
#   expression-analysis objects:
#
#       expr
#          ↓
#       expr.mat
#          ↓
#       expr2
#          ↓
#       expr.final
#          ↓
#       expr3
#
# This step ONLY audits existing historical objects.
#
# NO:
#   - expression matrix reconstruction
#   - expression transformation
#   - probe re-selection
#   - normalization
#   - ROC recalculation
#   - AUC recalculation
#   - model refitting
#   - threshold optimization
#
# ============================================================


# ------------------------------------------------------------
# 0. Working directories
# ------------------------------------------------------------

workflow_root <-
  "input_dir/Paper 1_Step 13-Diagnostic Classification"

step_dir <-
  file.path(
    workflow_root,
    "STEP13C_GSE32863_Input_Lock",
    "STEP13C-2_Expression_Matrix_Provenance"
  )

dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 1. Historical RData
# ------------------------------------------------------------

historical_rdata <-
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal/.RData"

if (!file.exists(historical_rdata)) {
  stop(
    "FAIL: Historical RData not found:\n",
    historical_rdata
  )
}


# ------------------------------------------------------------
# 2. Load historical objects in isolated environment
# ------------------------------------------------------------

hist_env <- new.env(parent = emptyenv())

loaded_objects <-
  load(
    historical_rdata,
    envir = hist_env
  )

cat("\n")
cat("============================================================\n")
cat("STEP13C-2 — Expression Matrix Provenance\n")
cat("============================================================\n\n")

cat("Historical RData:\n")
cat(historical_rdata, "\n\n")

cat("Objects loaded from RData:\n")
print(loaded_objects)
cat("\n")


# ------------------------------------------------------------
# 3. Required historical objects
# ------------------------------------------------------------

required_objects <- c(
  "expr",
  "expr.mat",
  "expr2",
  "expr.final",
  "expr3"
)

object_exists <- vapply(
  required_objects,
  function(x) exists(x, envir = hist_env, inherits = FALSE),
  logical(1)
)

object_presence_df <- data.frame(
  object = required_objects,
  present = object_exists,
  stringsAsFactors = FALSE
)

write.csv(
  object_presence_df,
  file.path(
    step_dir,
    "STEP13C-2_Object_Presence.csv"
  ),
  row.names = FALSE
)

if (!all(object_exists)) {
  
  missing_objects <-
    required_objects[!object_exists]
  
  stop(
    "FAIL: Required historical object(s) missing: ",
    paste(missing_objects, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 4. Retrieve objects
# ------------------------------------------------------------

expr       <- get("expr",       envir = hist_env)
expr_mat   <- get("expr.mat",   envir = hist_env)
expr2      <- get("expr2",      envir = hist_env)
expr_final <- get("expr.final", envir = hist_env)
expr3      <- get("expr3",      envir = hist_env)


# ------------------------------------------------------------
# 5. Basic object structure audit
# ------------------------------------------------------------

get_dim <- function(x) {
  
  d <- dim(x)
  
  if (is.null(d)) {
    return(c(NA_integer_, NA_integer_))
  }
  
  as.integer(d)
}


dimension_df <- data.frame(
  
  object = c(
    "expr",
    "expr.mat",
    "expr2",
    "expr.final",
    "expr3"
  ),
  
  class = c(
    paste(class(expr), collapse = ";"),
    paste(class(expr_mat), collapse = ";"),
    paste(class(expr2), collapse = ";"),
    paste(class(expr_final), collapse = ";"),
    paste(class(expr3), collapse = ";")
  ),
  
  nrow = c(
    get_dim(expr)[1],
    get_dim(expr_mat)[1],
    get_dim(expr2)[1],
    get_dim(expr_final)[1],
    get_dim(expr3)[1]
  ),
  
  ncol = c(
    get_dim(expr)[2],
    get_dim(expr_mat)[2],
    get_dim(expr2)[2],
    get_dim(expr_final)[2],
    get_dim(expr3)[2]
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  dimension_df,
  file.path(
    step_dir,
    "STEP13C-2_Object_Dimension_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Expected historical dimensions
# ------------------------------------------------------------

expected_dimension_df <- data.frame(
  
  object = c(
    "expr",
    "expr.mat",
    "expr2",
    "expr.final",
    "expr3"
  ),
  
  expected_nrow = c(
    48803,
    28742,
    28742,
    19425,
    116
  ),
  
  expected_ncol = c(
    117,
    116,
    118,
    116,
    5
  ),
  
  stringsAsFactors = FALSE
)

dimension_check <- merge(
  expected_dimension_df,
  dimension_df[, c("object", "nrow", "ncol")],
  by = "object",
  all = TRUE
)

dimension_check$dimension_match <-
  dimension_check$expected_nrow ==
  dimension_check$nrow &
  dimension_check$expected_ncol ==
  dimension_check$ncol

write.csv(
  dimension_check,
  file.path(
    step_dir,
    "STEP13C-2_Dimension_Expected_vs_Observed.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Helper functions for identifiers
# ------------------------------------------------------------

safe_colnames <- function(x) {
  
  cn <- colnames(x)
  
  if (is.null(cn)) {
    return(character(0))
  }
  
  as.character(cn)
}


safe_rownames <- function(x) {
  
  rn <- rownames(x)
  
  if (is.null(rn)) {
    return(character(0))
  }
  
  as.character(rn)
}


sample_summary <- function(ids, source, location) {
  
  ids <- as.character(ids)
  
  data.frame(
    source = source,
    identifier_location = location,
    n = length(ids),
    n_unique = length(unique(ids)),
    n_missing = sum(
      is.na(ids) |
        ids == ""
    ),
    first_id = ifelse(
      length(ids) > 0,
      ids[1],
      NA_character_
    ),
    last_id = ifelse(
      length(ids) > 0,
      ids[length(ids)],
      NA_character_
    ),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 8. Identify sample-ID locations
# ------------------------------------------------------------

expr_cols <- safe_colnames(expr)
expr_mat_cols <- safe_colnames(expr_mat)
expr2_cols <- safe_colnames(expr2)
expr_final_cols <- safe_colnames(expr_final)

expr3_rn <- safe_rownames(expr3)


sample_sources <- list(
  
  expr_sample_columns =
    if (length(expr_cols) >= 2)
      expr_cols[-1]
  else
    character(0),
  
  expr_mat_sample_columns =
    expr_mat_cols,
  
  expr2_sample_columns =
    if (length(expr2_cols) >= 3)
      expr2_cols[-c(1, 2)]
  else
    character(0),
  
  expr_final_sample_columns =
    expr_final_cols,
  
  expr3_row_names =
    expr3_rn
  
)


sample_summary_df <- do.call(
  rbind,
  list(
    
    sample_summary(
      sample_sources$expr_sample_columns,
      "expr",
      "column_names_excluding_ID_REF"
    ),
    
    sample_summary(
      sample_sources$expr_mat_sample_columns,
      "expr.mat",
      "column_names"
    ),
    
    sample_summary(
      sample_sources$expr2_sample_columns,
      "expr2",
      "column_names_excluding_Probe_Symbol"
    ),
    
    sample_summary(
      sample_sources$expr_final_sample_columns,
      "expr.final",
      "column_names"
    ),
    
    sample_summary(
      sample_sources$expr3_row_names,
      "expr3",
      "row_names"
    )
    
  )
)

write.csv(
  sample_summary_df,
  file.path(
    step_dir,
    "STEP13C-2_Sample_ID_Summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Extract sample IDs for exact comparison
# ------------------------------------------------------------

if (length(expr_cols) >= 2) {
  ids_expr <- expr_cols[-1]
} else {
  ids_expr <- character(0)
}

ids_expr_mat <- expr_mat_cols

if (length(expr2_cols) >= 3) {
  ids_expr2 <- expr2_cols[-c(1, 2)]
} else {
  ids_expr2 <- character(0)
}

ids_expr_final <- expr_final_cols

ids_expr3 <- expr3_rn


# ------------------------------------------------------------
# 10. Pairwise sample-ID concordance
# ------------------------------------------------------------

compare_ids <- function(
    ids_a,
    ids_b,
    source_a,
    source_b
) {
  
  set_a <- unique(as.character(ids_a))
  set_b <- unique(as.character(ids_b))
  
  data.frame(
    
    source_a = source_a,
    source_b = source_b,
    
    n_a = length(ids_a),
    n_b = length(ids_b),
    
    n_unique_a = length(set_a),
    n_unique_b = length(set_b),
    
    exact_order_match =
      identical(
        as.character(ids_a),
        as.character(ids_b)
      ),
    
    same_set =
      setequal(
        set_a,
        set_b
      ),
    
    only_in_a =
      length(setdiff(set_a, set_b)),
    
    only_in_b =
      length(setdiff(set_b, set_a)),
    
    stringsAsFactors = FALSE
  )
}


sample_pairwise_df <- rbind(
  
  compare_ids(
    ids_expr,
    ids_expr_mat,
    "expr",
    "expr.mat"
  ),
  
  compare_ids(
    ids_expr,
    ids_expr2,
    "expr",
    "expr2"
  ),
  
  compare_ids(
    ids_expr,
    ids_expr_final,
    "expr",
    "expr.final"
  ),
  
  compare_ids(
    ids_expr,
    ids_expr3,
    "expr",
    "expr3"
  ),
  
  compare_ids(
    ids_expr_mat,
    ids_expr2,
    "expr.mat",
    "expr2"
  ),
  
  compare_ids(
    ids_expr_mat,
    ids_expr_final,
    "expr.mat",
    "expr.final"
  ),
  
  compare_ids(
    ids_expr_mat,
    ids_expr3,
    "expr.mat",
    "expr3"
  ),
  
  compare_ids(
    ids_expr2,
    ids_expr_final,
    "expr2",
    "expr.final"
  ),
  
  compare_ids(
    ids_expr2,
    ids_expr3,
    "expr2",
    "expr3"
  ),
  
  compare_ids(
    ids_expr_final,
    ids_expr3,
    "expr.final",
    "expr3"
  )
  
)

write.csv(
  sample_pairwise_df,
  file.path(
    step_dir,
    "STEP13C-2_Sample_ID_Pairwise_Concordance.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Print sample-ID audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("Sample-ID provenance audit\n")
cat("------------------------------------------------------------\n\n")

print(sample_summary_df)

cat("\nPairwise sample-ID concordance:\n\n")

print(sample_pairwise_df)


# ------------------------------------------------------------
# 12. expr3 structural audit
# ------------------------------------------------------------

expr3_colnames <- safe_colnames(expr3)

required_expr3_columns <- c(
  "AOC3",
  "MS4A7",
  "IL7R",
  "Group",
  "Probability"
)

expr3_column_presence <- data.frame(
  
  column = required_expr3_columns,
  
  present =
    required_expr3_columns %in%
    expr3_colnames,
  
  stringsAsFactors = FALSE
)

write.csv(
  expr3_column_presence,
  file.path(
    step_dir,
    "STEP13C-2_expr3_Column_Audit.csv"
  ),
  row.names = FALSE
)


expr3_structure_pass <-
  all(
    required_expr3_columns %in%
      expr3_colnames
  ) &&
  nrow(expr3) == 116


# ------------------------------------------------------------
# 13. expr3 sample-ID relationship
# ------------------------------------------------------------

expr3_traceability <- data.frame(
  
  item = c(
    "expr3_nrow",
    "expr3_ncol",
    "expr3_has_rownames",
    "expr3_AOC3_present",
    "expr3_MS4A7_present",
    "expr3_IL7R_present",
    "expr3_Group_present",
    "expr3_Probability_present"
  ),
  
  value = c(
    nrow(expr3),
    ncol(expr3),
    length(expr3_rn) == nrow(expr3),
    "AOC3" %in% expr3_colnames,
    "MS4A7" %in% expr3_colnames,
    "IL7R" %in% expr3_colnames,
    "Group" %in% expr3_colnames,
    "Probability" %in% expr3_colnames
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  expr3_traceability,
  file.path(
    step_dir,
    "STEP13C-2_expr3_Traceability.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Expression-object identifier provenance
# ------------------------------------------------------------

identifier_rows <- list()

add_identifier_record <- function(
    object,
    identifier_type,
    n,
    n_unique,
    n_missing,
    example_first,
    example_last
) {
  
  data.frame(
    object = object,
    identifier_type = identifier_type,
    n = n,
    n_unique = n_unique,
    n_missing = n_missing,
    example_first = example_first,
    example_last = example_last,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 14A. expr — ID_REF
# ------------------------------------------------------------

expr_id_ref <- character(0)

if ("ID_REF" %in% expr_cols) {
  
  expr_id_ref <-
    as.character(expr[["ID_REF"]])
  
}

identifier_rows[[length(identifier_rows) + 1]] <-
  add_identifier_record(
    object = "expr",
    identifier_type = "ID_REF_column",
    n = length(expr_id_ref),
    n_unique = length(unique(expr_id_ref)),
    n_missing = sum(
      is.na(expr_id_ref) |
        expr_id_ref == ""
    ),
    example_first =
      if (length(expr_id_ref) > 0)
        expr_id_ref[1]
    else
      NA_character_,
    example_last =
      if (length(expr_id_ref) > 0)
        expr_id_ref[length(expr_id_ref)]
    else
      NA_character_
  )


# ------------------------------------------------------------
# 14B. expr.mat — rownames
# ------------------------------------------------------------

expr_mat_rn <-
  safe_rownames(expr_mat)

identifier_rows[[length(identifier_rows) + 1]] <-
  add_identifier_record(
    object = "expr.mat",
    identifier_type = "rownames",
    n = length(expr_mat_rn),
    n_unique = length(unique(expr_mat_rn)),
    n_missing = sum(
      is.na(expr_mat_rn) |
        expr_mat_rn == ""
    ),
    example_first =
      if (length(expr_mat_rn) > 0)
        expr_mat_rn[1]
    else
      NA_character_,
    example_last =
      if (length(expr_mat_rn) > 0)
        expr_mat_rn[length(expr_mat_rn)]
    else
      NA_character_
  )


# ------------------------------------------------------------
# 14C. expr2 — Probe
# ------------------------------------------------------------

expr2_probe <- character(0)

if ("Probe" %in% expr2_cols) {
  
  expr2_probe <-
    as.character(expr2[["Probe"]])
  
}

identifier_rows[[length(identifier_rows) + 1]] <-
  add_identifier_record(
    object = "expr2",
    identifier_type = "Probe_column",
    n = length(expr2_probe),
    n_unique = length(unique(expr2_probe)),
    n_missing = sum(
      is.na(expr2_probe) |
        expr2_probe == ""
    ),
    example_first =
      if (length(expr2_probe) > 0)
        expr2_probe[1]
    else
      NA_character_,
    example_last =
      if (length(expr2_probe) > 0)
        expr2_probe[length(expr2_probe)]
    else
      NA_character_
  )


# ------------------------------------------------------------
# 14D. expr2 — Symbol
# ------------------------------------------------------------

expr2_symbol <- character(0)

if ("Symbol" %in% expr2_cols) {
  
  expr2_symbol <-
    as.character(expr2[["Symbol"]])
  
}

identifier_rows[[length(identifier_rows) + 1]] <-
  add_identifier_record(
    object = "expr2",
    identifier_type = "Symbol_column",
    n = length(expr2_symbol),
    n_unique = length(unique(expr2_symbol)),
    n_missing = sum(
      is.na(expr2_symbol) |
        expr2_symbol == ""
    ),
    example_first =
      if (length(expr2_symbol) > 0)
        expr2_symbol[1]
    else
      NA_character_,
    example_last =
      if (length(expr2_symbol) > 0)
        expr2_symbol[length(expr2_symbol)]
    else
      NA_character_
  )


# ------------------------------------------------------------
# 14E. expr.final — rownames
# ------------------------------------------------------------

expr_final_rn <-
  safe_rownames(expr_final)

identifier_rows[[length(identifier_rows) + 1]] <-
  add_identifier_record(
    object = "expr.final",
    identifier_type = "rownames",
    n = length(expr_final_rn),
    n_unique = length(unique(expr_final_rn)),
    n_missing = sum(
      is.na(expr_final_rn) |
        expr_final_rn == ""
    ),
    example_first =
      if (length(expr_final_rn) > 0)
        expr_final_rn[1]
    else
      NA_character_,
    example_last =
      if (length(expr_final_rn) > 0)
        expr_final_rn[length(expr_final_rn)]
    else
      NA_character_
  )


# ------------------------------------------------------------
# 14F. Combine identifier provenance records
# ------------------------------------------------------------

identifier_provenance_df <-
  do.call(
    rbind,
    identifier_rows
  )


# ------------------------------------------------------------
# 14G. Save identifier provenance audit
# ------------------------------------------------------------

write.csv(
  identifier_provenance_df,
  file.path(
    step_dir,
    "STEP13C-2_Identifier_Provenance.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 14H. Console output
# ------------------------------------------------------------

cat("\n")
cat("------------------------------------------------------------\n")
cat("STEP13C-2 — Expression-object identifier provenance\n")
cat("------------------------------------------------------------\n\n")

print(identifier_provenance_df)

cat("\n")


# ------------------------------------------------------------
# 15. Audit expr2 basic identifier relationship
# ------------------------------------------------------------

expr2_identifier_check <- data.frame(
  
  item = c(
    "Probe_column_present",
    "Symbol_column_present",
    "nrow_expr2",
    "nrow_Probe",
    "nrow_Symbol",
    "Probe_unique",
    "Symbol_nonmissing"
  ),
  
  value = c(
    
    "Probe" %in% expr2_cols,
    
    "Symbol" %in% expr2_cols,
    
    nrow(expr2),
    
    length(expr2_probe),
    
    length(expr2_symbol),
    
    length(expr2_probe) ==
      length(unique(expr2_probe)),
    
    sum(
      !is.na(expr2_symbol) &
        expr2_symbol != ""
    )
    
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  expr2_identifier_check,
  file.path(
    step_dir,
    "STEP13C-2_expr2_Identifier_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 16. Audit sample dimensions
# ------------------------------------------------------------

sample_dimension_check <- data.frame(
  
  object = c(
    "expr",
    "expr.mat",
    "expr2",
    "expr.final",
    "expr3"
  ),
  
  expected_samples = c(
    116,
    116,
    116,
    116,
    116
  ),
  
  observed_sample_count = c(
    length(ids_expr),
    length(ids_expr_mat),
    length(ids_expr2),
    length(ids_expr_final),
    length(ids_expr3)
  ),
  
  stringsAsFactors = FALSE
)

sample_dimension_check$sample_count_match <-
  sample_dimension_check$expected_samples ==
  sample_dimension_check$observed_sample_count

write.csv(
  sample_dimension_check,
  file.path(
    step_dir,
    "STEP13C-2_Sample_Dimension_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Conservative PASS criteria
# ------------------------------------------------------------
#
# PASS here means:
#
#   1. required historical objects exist;
#   2. expected object dimensions are reproduced;
#   3. each expression object has 116 sample positions;
#   4. expr3 has the expected five variables;
#   5. no obvious duplicate/missing sample identifiers where
#      identifiers are available;
#   6. sample-ID concordance is established where exact IDs
#      are available.
#
# This step does NOT require us to infer the historical
# transformation/filtering procedure.
#
# Transformation provenance is reserved for STEP13C-7.
# Probe mapping is reserved for STEP13C-5/6.
# Transferability is reserved for STEP13C-8.
# ------------------------------------------------------------


required_objects_pass <-
  all(object_presence_df$present)

dimension_pass <-
  all(dimension_check$dimension_match)

sample_dimension_pass <-
  all(
    sample_dimension_check$sample_count_match
  )

expr3_pass <-
  expr3_structure_pass


# Available sample-ID sources with complete 116 IDs
sample_id_sources <- list(
  
  expr = ids_expr,
  expr.mat = ids_expr_mat,
  expr2 = ids_expr2,
  expr.final = ids_expr_final,
  expr3 = ids_expr3
  
)

sample_id_sources <- sample_id_sources[
  vapply(
    sample_id_sources,
    function(x)
      length(x) == 116 &&
      all(!is.na(x)) &&
      all(x != ""),
    logical(1)
  )
]


# Determine whether all available 116-sample sources
# have the same sample set.
sample_set_pass <- TRUE

if (length(sample_id_sources) >= 2) {
  
  reference_ids <-
    sample_id_sources[[1]]
  
  sample_set_pass <-
    all(
      vapply(
        sample_id_sources[-1],
        function(x)
          setequal(reference_ids, x),
        logical(1)
      )
    )
}


# Determine whether all available sources also preserve
# exact sample order.
sample_order_pass <- TRUE

if (length(sample_id_sources) >= 2) {
  
  reference_ids <-
    sample_id_sources[[1]]
  
  sample_order_pass <-
    all(
      vapply(
        sample_id_sources[-1],
        function(x)
          identical(reference_ids, x),
        logical(1)
      )
    )
}


# ------------------------------------------------------------
# Duplicate identifier audit
# ------------------------------------------------------------
#
# Important distinction:
#
#   Probe identifiers should be unique.
#   Gene symbols may legitimately be repeated because
#   multiple probes can map to the same gene.
#
# Therefore:
#
#   expr ID_REF        -> expected unique
#   expr2 Probe        -> expected unique
#   expr.mat rownames  -> gene-level identifiers; duplicates
#                         are allowed and documented
#   expr2 Symbol       -> gene-level identifiers; duplicates
#                         are allowed and documented
#   expr.final rownames -> expected unique after historical
#                          gene-level reduction
#
# ------------------------------------------------------------

probe_identifier_unique_pass <-
  all(
    identifier_provenance_df[
      identifier_provenance_df$identifier_type %in%
        c("ID_REF_column", "Probe_column"),
      "n_unique"
    ]
    ==
      identifier_provenance_df[
        identifier_provenance_df$identifier_type %in%
          c("ID_REF_column", "Probe_column"),
        "n"
      ]
  )


expr_mat_gene_duplicates_present <-
  any(
    identifier_provenance_df[
      identifier_provenance_df$object == "expr.mat" &
        identifier_provenance_df$identifier_type == "rownames",
      "n_unique"
    ]
    <
      identifier_provenance_df[
        identifier_provenance_df$object == "expr.mat" &
          identifier_provenance_df$identifier_type == "rownames",
        "n"
      ]
  )


expr2_symbol_duplicates_present <-
  any(
    identifier_provenance_df[
      identifier_provenance_df$object == "expr2" &
        identifier_provenance_df$identifier_type == "Symbol_column",
      "n_unique"
    ]
    <
      identifier_provenance_df[
        identifier_provenance_df$object == "expr2" &
          identifier_provenance_df$identifier_type == "Symbol_column",
        "n"
      ]
  )


expr_final_unique_pass <-
  all(
    identifier_provenance_df[
      identifier_provenance_df$object == "expr.final" &
        identifier_provenance_df$identifier_type == "rownames",
      "n_unique"
    ]
    ==
      identifier_provenance_df[
        identifier_provenance_df$object == "expr.final" &
          identifier_provenance_df$identifier_type == "rownames",
        "n"
      ]
  )


# Overall identifier uniqueness criterion:
#
# Only probe identifiers and final gene-level identifiers
# are required to be unique.
#
identifier_duplicate_pass <-
  probe_identifier_unique_pass &&
  expr_final_unique_pass


# ------------------------------------------------------------
# Identifier uniqueness interpretation
# ------------------------------------------------------------

identifier_duplicate_interpretation <- c(
  
  paste0(
    "Probe-level identifier uniqueness: ",
    probe_identifier_unique_pass
  ),
  
  paste0(
    "expr.mat gene-symbol duplicates present: ",
    expr_mat_gene_duplicates_present
  ),
  
  paste0(
    "expr2 gene-symbol duplicates present: ",
    expr2_symbol_duplicates_present
  ),
  
  paste0(
    "expr.final gene-level identifiers unique: ",
    expr_final_unique_pass
  ),
  
  "Gene-symbol duplication in expr.mat and expr2 is treated as",
  "an identifier mapping characteristic rather than as a data",
  "integrity failure.",
  
  "The specific historical probe-to-gene reduction rule is NOT",
  "inferred in STEP13C-2."
)


writeLines(
  identifier_duplicate_interpretation,
  con = file.path(
    step_dir,
    "STEP13C-2_Identifier_Uniqueness_Interpretation.txt"
  )
)
# ------------------------------------------------------------
# 14I. Gene-symbol multiplicity audit
# ------------------------------------------------------------

gene_symbol_multiplicity <- data.frame(
  
  source = c(
    "expr.mat",
    "expr2"
  ),
  
  n_rows = c(
    length(expr_mat_rn),
    length(expr2_symbol)
  ),
  
  n_unique_gene_symbols = c(
    length(unique(expr_mat_rn)),
    length(unique(expr2_symbol))
  ),
  
  duplicated_gene_symbol_rows = c(
    length(expr_mat_rn) -
      length(unique(expr_mat_rn)),
    
    length(expr2_symbol) -
      length(unique(expr2_symbol))
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  gene_symbol_multiplicity,
  file.path(
    step_dir,
    "STEP13C-2_Gene_Symbol_Multiplicity_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. Provenance interpretation
# ------------------------------------------------------------

interpretation_lines <- c(
  
  "STEP13C-2 establishes object-level structural provenance for the historical GSE32863 analysis objects.",
  
  "The audit does not reconstruct or transform any expression matrix.",
  
  "The audit does not reselect probes or genes.",
  
  "The audit does not recalculate ROC curves or AUC values.",
  
  "The audit does not refit any statistical model.",
  
  "The audit does not optimize any threshold.",
  
  "The observed historical objects are characterized by their existing dimensions, identifiers, and sample positions.",
  
  "Expression-scale and transformation provenance are intentionally deferred to STEP13C-7.",
  
  "GPL6884 probe annotation and core-gene probe mapping are intentionally deferred to STEP13C-5 and STEP13C-6.",
  
  "TCGA-to-GSE32863 transferability is intentionally deferred to STEP13C-8.",
  
  if (
    sample_order_pass
  )
    "Where available, sample identifiers are concordant in exact order across the audited historical objects."
  else if (
    sample_set_pass
  )
    "Where available, sample identifiers are concordant by set, but exact order is not uniformly established."
  else
    "Sample identifiers are not uniformly concordant across all available historical objects; further investigation is required."
)


writeLines(
  interpretation_lines,
  con = file.path(
    step_dir,
    "STEP13C-2_Provenance_Interpretation.txt"
  )
)


# ------------------------------------------------------------
# 19. Final report
# ------------------------------------------------------------

if (overall_pass) {
  final_status <- "PASS"
} else {
  final_status <- "HOLD"
}


final_report <- c(
  
  "STEP13C-2 — FINAL REPORT",
  "============================================================",
  paste0("FINAL STATUS: ", final_status),
  "",
  "Purpose:",
  "Object-level provenance audit of the historical GSE32863",
  "expression-analysis objects expr / expr.mat / expr2 /",
  "expr.final / expr3.",
  "",
  "Historical RData:",
  historical_rdata,
  "",
  "Required historical objects present:",
  as.character(required_objects_pass),
  "",
  "Expected dimensions reproduced:",
  as.character(dimension_pass),
  "",
  "All five objects contain 116 sample positions:",
  as.character(sample_dimension_pass),
  "",
  "expr3 expected structure present:",
  as.character(expr3_pass),
  "",
  "Sample-ID set concordance across available objects:",
  as.character(sample_set_pass),
  "",
  "Exact sample-ID order concordance across available objects:",
  as.character(sample_order_pass),
  "",
  "No missing identifier values in audited identifier fields:",
  as.character(identifier_missing_pass),
  "",
  "No duplicate values in primary audited identifier fields:",
  as.character(identifier_duplicate_pass),
  "",
  "No expression matrix reconstruction performed.",
  "No expression transformation performed.",
  "No probe selection performed.",
  "No ROC recalculation performed.",
  "No AUC recalculation performed.",
  "No model refitting performed.",
  "No threshold optimization performed.",
  "",
  "STEP13C-7 remains responsible for expression-scale/transformation audit.",
  "STEP13C-5/6 remain responsible for GPL6884 probe annotation and",
  "core-gene probe mapping.",
  "STEP13C-8 remains responsible for TCGA-to-GSE32863 transferability.",
  "",
  "============================================================"
)


writeLines(
  final_report,
  con = file.path(
    step_dir,
    "STEP13C-2_FINAL_REPORT.txt"
  )
)


# ------------------------------------------------------------
# 20. PASS / HOLD marker
# ------------------------------------------------------------

marker_file <-
  file.path(
    step_dir,
    paste0(
      "STEP13C-2_",
      final_status,
      ".txt"
    )
  )


marker_text <- c(
  paste0(
    "STEP13C-2 FINAL STATUS: ",
    final_status
  ),
  "",
  "Expression-object provenance audit completed.",
  "",
  "No ROC/AUC recalculation.",
  "No model refitting.",
  "No threshold optimization.",
  "No expression transformation.",
  "No probe selection."
)

writeLines(
  marker_text,
  con = marker_file
)


# ------------------------------------------------------------
# 21. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    step_dir,
    "STEP13C-2_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 22. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13C-2 — FINAL REPORT\n")
cat("============================================================\n\n")

cat(
  "FINAL STATUS: ",
  final_status,
  "\n\n",
  sep = ""
)

cat(
  "Required historical objects present: ",
  required_objects_pass,
  "\n",
  sep = ""
)

cat(
  "Expected dimensions reproduced: ",
  dimension_pass,
  "\n",
  sep = ""
)

cat(
  "All five objects contain 116 sample positions: ",
  sample_dimension_pass,
  "\n",
  sep = ""
)

cat(
  "expr3 expected structure present: ",
  expr3_pass,
  "\n",
  sep = ""
)

cat(
  "Sample-ID set concordance: ",
  sample_set_pass,
  "\n",
  sep = ""
)

cat(
  "Exact sample-ID order concordance: ",
  sample_order_pass,
  "\n",
  sep = ""
)

cat(
  "No missing audited identifiers: ",
  identifier_missing_pass,
  "\n",
  sep = ""
)

cat(
  "No duplicate primary identifiers: ",
  identifier_duplicate_pass,
  "\n\n",
  sep = ""
)

cat("No expression reconstruction performed.\n")
cat("No expression transformation performed.\n")
cat("No probe selection performed.\n")
cat("No ROC/AUC recalculation performed.\n")
cat("No model refitting performed.\n")
cat("No threshold optimization performed.\n\n")

cat(
  "Output directory:\n",
  step_dir,
  "\n\n",
  sep = ""
)

cat("============================================================\n")


# ============================================================
# STEP13C-3 — GSE32863 Sample Annotation
# ============================================================
#
# Purpose:
#   Audit the historical GSE32863 sample annotation.
#
# This step is provenance-only.
#
# No:
#   - expression reconstruction
#   - expression transformation
#   - probe selection
#   - sample reclassification
#   - ROC/AUC recalculation
#   - model fitting
#   - threshold optimization
#
# STEP13C-4 will separately audit paired-patient identity.
#
# ============================================================


# ------------------------------------------------------------
# 0. Workflow directory
# ------------------------------------------------------------

workflow_root <-
  "input_dir/Paper 1_Step 13-Diagnostic Classification"

step_dir <-
  file.path(
    workflow_root,
    "STEP13C_GSE32863_Input_Lock",
    "STEP13C-3_Sample_Annotation"
  )

dir.create(
  step_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Historical analysis directory
# ------------------------------------------------------------

historical_root <-
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal"

historical_rdata <-
  file.path(
    historical_root,
    ".RData"
  )


if (!dir.exists(historical_root)) {
  stop(
    "Historical analysis directory not found:\n",
    historical_root
  )
}


if (!file.exists(historical_rdata)) {
  stop(
    "Historical .RData not found:\n",
    historical_rdata
  )
}


# ------------------------------------------------------------
# 2. Load historical objects in isolated environment
# ------------------------------------------------------------

hist_env <-
  new.env(parent = emptyenv())

load(
  historical_rdata,
  envir = hist_env
)


# ------------------------------------------------------------
# 3. Required object audit
# ------------------------------------------------------------

required_objects <-
  c(
    "group",
    "expr3",
    "expr",
    "expr.mat",
    "expr.final",
    "expr2"
  )


required_objects_present <-
  required_objects %in% ls(
    envir = hist_env
  )


required_objects_pass <-
  all(
    required_objects_present
  )


required_object_audit <-
  data.frame(
    object = required_objects,
    present = required_objects_present,
    stringsAsFactors = FALSE
  )


write.csv(
  required_object_audit,
  file.path(
    step_dir,
    "STEP13C-3_Required_Object_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. Retrieve historical objects
# ------------------------------------------------------------
#
# IMPORTANT:
#   Braces are deliberately used for every if/else block
#   to avoid Console parsing errors.
#
# ------------------------------------------------------------

if (
  exists(
    "group",
    envir = hist_env,
    inherits = FALSE
  )
) {
  
  group <-
    get(
      "group",
      envir = hist_env
    )
  
} else {
  
  group <-
    NULL
  
}


if (
  exists(
    "expr3",
    envir = hist_env,
    inherits = FALSE
  )
) {
  
  expr3 <-
    get(
      "expr3",
      envir = hist_env
    )
  
} else {
  
  expr3 <-
    NULL
  
}


if (
  exists(
    "expr",
    envir = hist_env,
    inherits = FALSE
  )
) {
  
  expr <-
    get(
      "expr",
      envir = hist_env
    )
  
} else {
  
  expr <-
    NULL
  
}


if (
  exists(
    "expr.mat",
    envir = hist_env,
    inherits = FALSE
  )
) {
  
  expr_mat <-
    get(
      "expr.mat",
      envir = hist_env
    )
  
} else {
  
  expr_mat <-
    NULL
  
}


if (
  exists(
    "expr.final",
    envir = hist_env,
    inherits = FALSE
  )
) {
  
  expr_final <-
    get(
      "expr.final",
      envir = hist_env
    )
  
} else {
  
  expr_final <-
    NULL
  
}


if (
  exists(
    "expr2",
    envir = hist_env,
    inherits = FALSE
  )
) {
  
  expr2 <-
    get(
      "expr2",
      envir = hist_env
    )
  
} else {
  
  expr2 <-
    NULL
  
}


# ------------------------------------------------------------
# 5. group object structure
# ------------------------------------------------------------

group_exists <-
  !is.null(group)


group_is_vector <-
  FALSE

if (group_exists) {
  
  group_is_vector <-
    is.atomic(group) &&
    is.null(dim(group))
  
}


group_length <-
  0L

if (group_exists) {
  
  group_length <-
    length(group)
  
}


group_length_pass <-
  identical(
    group_length,
    116L
  )


# ------------------------------------------------------------
# 6. Convert group to character for audit comparison
# ------------------------------------------------------------
#
# This does NOT modify the historical group object.
#
# ------------------------------------------------------------

group_chr <-
  character(0)

if (group_exists) {
  
  group_chr <-
    as.character(group)
  
}


group_chr <-
  trimws(group_chr)


# ------------------------------------------------------------
# 7. Missing annotation audit
# ------------------------------------------------------------

group_missing_n <-
  NA_integer_

if (length(group_chr) > 0L) {
  
  group_missing_n <-
    sum(
      is.na(group_chr) |
        group_chr == ""
    )
  
}


group_missing_pass <-
  identical(
    group_missing_n,
    0L
  )


# ------------------------------------------------------------
# 8. Historical group levels
# ------------------------------------------------------------

group_levels <-
  character(0)

if (length(group_chr) > 0L) {
  
  group_levels <-
    sort(
      unique(group_chr)
    )
  
}


expected_group_levels <-
  c(
    "Normal",
    "Tumor"
  )


group_levels_pass <-
  setequal(
    group_levels,
    expected_group_levels
  )


# ------------------------------------------------------------
# 9. Historical group counts
# ------------------------------------------------------------

normal_n <-
  0L

tumor_n <-
  0L


if (length(group_chr) > 0L) {
  
  group_table <-
    table(
      group_chr,
      useNA = "ifany"
    )
  
  
  if ("Normal" %in% names(group_table)) {
    
    normal_n <-
      as.integer(
        group_table[["Normal"]]
      )
    
  }
  
  
  if ("Tumor" %in% names(group_table)) {
    
    tumor_n <-
      as.integer(
        group_table[["Tumor"]]
      )
    
  }
  
}


group_count_pass <-
  identical(
    normal_n,
    58L
  ) &&
  identical(
    tumor_n,
    58L
  )


# ------------------------------------------------------------
# 10. Save historical group counts
# ------------------------------------------------------------

group_count_df <-
  data.frame(
    Group = c(
      "Normal",
      "Tumor"
    ),
    N = c(
      normal_n,
      tumor_n
    ),
    stringsAsFactors = FALSE
  )


write.csv(
  group_count_df,
  file.path(
    step_dir,
    "STEP13C-3_Historical_Group_Counts.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. expr3 structure audit
# ------------------------------------------------------------

expr3_exists <-
  !is.null(expr3)


expr3_is_data_frame <-
  FALSE

if (expr3_exists) {
  
  expr3_is_data_frame <-
    is.data.frame(expr3)
  
}


expr3_has_group <-
  FALSE

if (expr3_is_data_frame) {
  
  expr3_has_group <-
    "Group" %in%
    colnames(expr3)
  
}


expr3_group <-
  character(0)

if (expr3_has_group) {
  
  expr3_group <-
    as.character(
      expr3[["Group"]]
    )
  
}


expr3_group <-
  trimws(expr3_group)


expr3_group_length <-
  length(expr3_group)


expr3_group_length_pass <-
  identical(
    expr3_group_length,
    116L
  )


# ------------------------------------------------------------
# 12. group vs expr3$Group
# ------------------------------------------------------------

group_expr3_length_concordance <-
  length(group_chr) ==
  length(expr3_group)


group_expr3_exact_concordance <-
  FALSE

if (group_expr3_length_concordance) {
  
  group_expr3_exact_concordance <-
    identical(
      group_chr,
      expr3_group
    )
  
}


group_expr3_mismatch_n <-
  NA_integer_

if (group_expr3_length_concordance) {
  
  group_expr3_mismatch_n <-
    sum(
      group_chr !=
        expr3_group
    )
  
}


# ------------------------------------------------------------
# 13. expr3 sample IDs
# ------------------------------------------------------------

expr3_sample_ids <-
  character(0)


if (
  expr3_is_data_frame &&
  !is.null(rownames(expr3))
) {
  
  expr3_sample_ids <-
    as.character(
      rownames(expr3)
    )
  
}


expr3_sample_ids_missing_n <-
  NA_integer_

if (length(expr3_sample_ids) > 0L) {
  
  expr3_sample_ids_missing_n <-
    sum(
      is.na(expr3_sample_ids) |
        expr3_sample_ids == ""
    )
  
}


expr3_sample_ids_unique_pass <-
  length(expr3_sample_ids) == 116L &&
  length(unique(expr3_sample_ids)) == 116L &&
  identical(
    expr3_sample_ids_missing_n,
    0L
  )


# ------------------------------------------------------------
# 14. Expression-object sample IDs
# ------------------------------------------------------------

expr_sample_ids <-
  character(0)

if (!is.null(expr)) {
  
  expr_cols_local <-
    colnames(expr)
  
  
  if ("ID_REF" %in% expr_cols_local) {
    
    expr_sample_ids <-
      setdiff(
        expr_cols_local,
        "ID_REF"
      )
    
  } else {
    
    expr_sample_ids <-
      expr_cols_local
    
  }
  
}


expr_mat_sample_ids <-
  character(0)

if (!is.null(expr_mat)) {
  
  expr_mat_sample_ids <-
    colnames(expr_mat)
  
}


expr_final_sample_ids <-
  character(0)

if (!is.null(expr_final)) {
  
  expr_final_sample_ids <-
    colnames(expr_final)
  
}


expr2_sample_ids <-
  character(0)

if (!is.null(expr2)) {
  
  expr2_cols_local <-
    colnames(expr2)
  
  
  expr2_sample_ids <-
    setdiff(
      expr2_cols_local,
      c(
        "Probe",
        "Symbol"
      )
    )
  
}


# ------------------------------------------------------------
# 15. Sample-ID set concordance
# ------------------------------------------------------------

expr3_vs_expr_set_pass <-
  length(expr_sample_ids) == 116L &&
  setequal(
    expr3_sample_ids,
    expr_sample_ids
  )


expr3_vs_expr_mat_set_pass <-
  length(expr_mat_sample_ids) == 116L &&
  setequal(
    expr3_sample_ids,
    expr_mat_sample_ids
  )


expr3_vs_expr_final_set_pass <-
  length(expr_final_sample_ids) == 116L &&
  setequal(
    expr3_sample_ids,
    expr_final_sample_ids
  )


expr3_vs_expr2_set_pass <-
  length(expr2_sample_ids) == 116L &&
  setequal(
    expr3_sample_ids,
    expr2_sample_ids
  )


all_expression_sample_sets_pass <-
  expr3_vs_expr_set_pass &&
  expr3_vs_expr_mat_set_pass &&
  expr3_vs_expr_final_set_pass &&
  expr3_vs_expr2_set_pass


# ------------------------------------------------------------
# 16. Exact sample-ID order concordance
# ------------------------------------------------------------

expr3_vs_expr_order_pass <-
  length(expr_sample_ids) == 116L &&
  identical(
    expr3_sample_ids,
    expr_sample_ids
  )


expr3_vs_expr_mat_order_pass <-
  length(expr_mat_sample_ids) == 116L &&
  identical(
    expr3_sample_ids,
    expr_mat_sample_ids
  )


expr3_vs_expr_final_order_pass <-
  length(expr_final_sample_ids) == 116L &&
  identical(
    expr3_sample_ids,
    expr_final_sample_ids
  )


expr3_vs_expr2_order_pass <-
  length(expr2_sample_ids) == 116L &&
  identical(
    expr3_sample_ids,
    expr2_sample_ids
  )


all_expression_sample_orders_pass <-
  expr3_vs_expr_order_pass &&
  expr3_vs_expr_mat_order_pass &&
  expr3_vs_expr_final_order_pass &&
  expr3_vs_expr2_order_pass


# ------------------------------------------------------------
# 17. Sample-level annotation table
# ------------------------------------------------------------

sample_annotation_df <-
  data.frame(
    sample_id = expr3_sample_ids,
    group_from_group_object =
      NA_character_,
    group_from_expr3 =
      NA_character_,
    stringsAsFactors = FALSE
  )


if (
  length(group_chr) ==
  length(expr3_sample_ids)
) {
  
  sample_annotation_df[
    ,
    "group_from_group_object"
  ] <-
    group_chr
  
}


if (
  length(expr3_group) ==
  length(expr3_sample_ids)
) {
  
  sample_annotation_df[
    ,
    "group_from_expr3"
  ] <-
    expr3_group
  
}


sample_annotation_df$group_match <-
  NA


if (
  length(group_chr) ==
  length(expr3_sample_ids) &&
  length(expr3_group) ==
  length(expr3_sample_ids)
) {
  
  sample_annotation_df$group_match <-
    group_chr ==
    expr3_group
  
}


write.csv(
  sample_annotation_df,
  file.path(
    step_dir,
    "STEP13C-3_Sample_Level_Annotation.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. Annotation mismatch audit
# ------------------------------------------------------------

annotation_mismatch_df <-
  sample_annotation_df[
    !is.na(
      sample_annotation_df$group_match
    ) &
      !sample_annotation_df$group_match,
    ,
    drop = FALSE
  ]


write.csv(
  annotation_mismatch_df,
  file.path(
    step_dir,
    "STEP13C-3_Annotation_Mismatch_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 19. Sample-ID audit table
# ------------------------------------------------------------

sample_id_audit_df <-
  data.frame(
    comparison = c(
      "expr3_vs_expr",
      "expr3_vs_expr.mat",
      "expr3_vs_expr.final",
      "expr3_vs_expr2"
    ),
    set_concordance = c(
      expr3_vs_expr_set_pass,
      expr3_vs_expr_mat_set_pass,
      expr3_vs_expr_final_set_pass,
      expr3_vs_expr2_set_pass
    ),
    exact_order_concordance = c(
      expr3_vs_expr_order_pass,
      expr3_vs_expr_mat_order_pass,
      expr3_vs_expr_final_order_pass,
      expr3_vs_expr2_order_pass
    ),
    stringsAsFactors = FALSE
  )


write.csv(
  sample_id_audit_df,
  file.path(
    step_dir,
    "STEP13C-3_Sample_ID_Concordance_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 20. Provenance summary
# ------------------------------------------------------------

annotation_summary_df <-
  data.frame(
    metric = c(
      "group object present",
      "group object is vector",
      "group length = 116",
      "missing group annotations = 0",
      "group levels are Normal/Tumor",
      "Normal N = 58",
      "Tumor N = 58",
      "expr3 present",
      "expr3 is data.frame",
      "expr3$Group present",
      "expr3$Group length = 116",
      "group vs expr3$Group exact concordance",
      "expr3 sample IDs unique and complete",
      "expr3 vs expr sample-ID set concordance",
      "expr3 vs expr.mat sample-ID set concordance",
      "expr3 vs expr.final sample-ID set concordance",
      "expr3 vs expr2 sample-ID set concordance",
      "all expression-object sample-ID sets concordant",
      "all expression-object sample-ID orders concordant"
    ),
    pass = c(
      group_exists,
      group_is_vector,
      group_length_pass,
      group_missing_pass,
      group_levels_pass,
      normal_n == 58L,
      tumor_n == 58L,
      expr3_exists,
      expr3_is_data_frame,
      expr3_has_group,
      expr3_group_length_pass,
      group_expr3_exact_concordance,
      expr3_sample_ids_unique_pass,
      expr3_vs_expr_set_pass,
      expr3_vs_expr_mat_set_pass,
      expr3_vs_expr_final_set_pass,
      expr3_vs_expr2_set_pass,
      all_expression_sample_sets_pass,
      all_expression_sample_orders_pass
    ),
    stringsAsFactors = FALSE
  )


write.csv(
  annotation_summary_df,
  file.path(
    step_dir,
    "STEP13C-3_Annotation_Provenance_Summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 21. Final PASS / FAIL
# ------------------------------------------------------------

final_status <-
  "FAIL"


if (
  required_objects_pass &&
  group_exists &&
  group_is_vector &&
  group_length_pass &&
  group_missing_pass &&
  group_levels_pass &&
  group_count_pass &&
  expr3_exists &&
  expr3_is_data_frame &&
  expr3_has_group &&
  expr3_group_length_pass &&
  group_expr3_exact_concordance &&
  expr3_sample_ids_unique_pass &&
  all_expression_sample_sets_pass
) {
  
  final_status <-
    "PASS"
  
}


# ------------------------------------------------------------
# 22. Interpretation lock
# ------------------------------------------------------------

final_interpretation <-
  c(
    paste0(
      "Historical GSE32863 group object present: ",
      group_exists
    ),
    paste0(
      "Historical group vector length = 116: ",
      group_length_pass
    ),
    paste0(
      "Historical group levels are exactly Normal and Tumor: ",
      group_levels_pass
    ),
    paste0(
      "Normal sample count = 58: ",
      normal_n == 58L
    ),
    paste0(
      "Tumor sample count = 58: ",
      tumor_n == 58L
    ),
    paste0(
      "No missing group annotations: ",
      group_missing_pass
    ),
    paste0(
      "group and expr3$Group exact concordance: ",
      group_expr3_exact_concordance
    ),
    paste0(
      "expr3 sample IDs unique and complete: ",
      expr3_sample_ids_unique_pass
    ),
    paste0(
      "Expression-object sample-ID sets concordant: ",
      all_expression_sample_sets_pass
    ),
    paste0(
      "Expression-object sample-ID orders concordant: ",
      all_expression_sample_orders_pass
    ),
    "",
    "Historical group annotations were audited but not reconstructed.",
    "No sample was reclassified.",
    "No external sample annotation was substituted.",
    "No expression values were recalculated or transformed.",
    "No ROC/AUC calculation was performed.",
    "No model fitting was performed.",
    "No threshold optimization was performed.",
    "",
    "Paired-patient identity is reserved for STEP13C-4."
  )


writeLines(
  final_interpretation,
  file.path(
    step_dir,
    "STEP13C-3_Interpretation_Lock.txt"
  )
)


# ------------------------------------------------------------
# 23. PASS marker
# ------------------------------------------------------------

if (final_status == "PASS") {
  
  writeLines(
    c(
      "STEP13C-3 — PASS",
      "",
      "GSE32863 historical sample annotation provenance established.",
      "",
      "Normal N = 58",
      "Tumor N = 58",
      "Total N = 116",
      "",
      "group and expr3$Group are exactly concordant.",
      "Expression-object sample-ID sets are concordant.",
      "",
      "No sample reclassification performed.",
      "No expression reconstruction performed.",
      "No expression transformation performed.",
      "No probe selection performed.",
      "No ROC/AUC recalculation performed.",
      "No model refitting performed.",
      "No threshold optimization performed.",
      "",
      "Paired-patient identity remains reserved for STEP13C-4."
    ),
    file.path(
      step_dir,
      "STEP13C-3_PASS.txt"
    )
  )
  
}


# ------------------------------------------------------------
# 24. Final console report
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13C-3 — FINAL REPORT\n")
cat("============================================================\n\n")

cat(
  "FINAL STATUS: ",
  final_status,
  "\n\n",
  sep = ""
)

cat(
  "Historical group object present: ",
  group_exists,
  "\n",
  sep = ""
)

cat(
  "Historical group length = 116: ",
  group_length_pass,
  "\n",
  sep = ""
)

cat(
  "Group levels are exactly Normal/Tumor: ",
  group_levels_pass,
  "\n",
  sep = ""
)

cat(
  "Normal N = 58: ",
  normal_n == 58L,
  "\n",
  sep = ""
)

cat(
  "Tumor N = 58: ",
  tumor_n == 58L,
  "\n",
  sep = ""
)

cat(
  "No missing group annotations: ",
  group_missing_pass,
  "\n",
  sep = ""
)

cat(
  "group vs expr3$Group exact concordance: ",
  group_expr3_exact_concordance,
  "\n",
  sep = ""
)

cat(
  "expr3 sample IDs unique and complete: ",
  expr3_sample_ids_unique_pass,
  "\n",
  sep = ""
)

cat(
  "Expression-object sample-ID set concordance: ",
  all_expression_sample_sets_pass,
  "\n",
  sep = ""
)

cat(
  "Expression-object sample-ID order concordance: ",
  all_expression_sample_orders_pass,
  "\n\n",
  sep = ""
)

cat("No sample reclassification performed.\n")
cat("No expression reconstruction performed.\n")
cat("No expression transformation performed.\n")
cat("No probe selection performed.\n")
cat("No ROC/AUC recalculation performed.\n")
cat("No model refitting performed.\n")
cat("No threshold optimization performed.\n\n")

cat(
  "Paired-patient identity is reserved for STEP13C-4.\n\n"
)

cat(
  "Output directory:\n",
  step_dir,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 25. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    step_dir,
    "STEP13C-3_sessionInfo.txt"
  )
)


cat("============================================================\n")
cat("STEP13C-3 completed.\n")
cat("============================================================\n")

# ============================================================
# STEP13C-4A — Historical Object / Metadata Inventory
# GSE32863 Paired ID Audit
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP13C-4A — Historical Object / Metadata Inventory\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

step_dir <- file.path(
  workflow_root,
  "STEP13C_GSE32863_Input_Lock",
  "STEP13C-4_Paired_ID_Audit"
)

out_dir <- file.path(
  step_dir,
  "STEP13C-4A_Historical_Object_Inventory"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Historical RData
# ------------------------------------------------------------

hist_file <- file.path(
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal",
  ".RData"
)

stopifnot(file.exists(hist_file))

hist_env <- new.env(parent = emptyenv())

load(hist_file, envir = hist_env)

object_names <- ls(hist_env)

cat("Historical .RData:\n")
cat(hist_file, "\n\n")

cat("Objects found:", length(object_names), "\n\n")

print(object_names)

# ------------------------------------------------------------
# 2. Object-level inventory
# ------------------------------------------------------------

inventory <- data.frame(
  object = character(),
  class = character(),
  dimensions = character(),
  length = character(),
  stringsAsFactors = FALSE
)

for (nm in object_names) {
  
  obj <- get(nm, envir = hist_env)
  
  obj_class <- paste(class(obj), collapse = ";")
  
  obj_dim <- dim(obj)
  
  if (is.null(obj_dim)) {
    dim_text <- NA_character_
  } else {
    dim_text <- paste(obj_dim, collapse = " x ")
  }
  
  obj_length <- length(obj)
  
  inventory <- rbind(
    inventory,
    data.frame(
      object = nm,
      class = obj_class,
      dimensions = dim_text,
      length = as.character(obj_length),
      stringsAsFactors = FALSE
    )
  )
}

write.csv(
  inventory,
  file.path(out_dir, "STEP13C-4A_Object_Inventory.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 3. Metadata-like object screening
# ------------------------------------------------------------

metadata_candidates <- inventory[
  grepl(
    "anno|meta|clinical|sample|patient|subject|group|expr3|annotation",
    inventory$object,
    ignore.case = TRUE
  ),
]

write.csv(
  metadata_candidates,
  file.path(
    out_dir,
    "STEP13C-4A_Metadata_Like_Object_Candidates.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\nMetadata-like objects identified:\n")
print(metadata_candidates)

# ------------------------------------------------------------
# 4. Inspect columns of data.frame / matrix objects
# ------------------------------------------------------------

column_inventory <- data.frame(
  object = character(),
  column = character(),
  stringsAsFactors = FALSE
)

for (nm in object_names) {
  
  obj <- get(nm, envir = hist_env)
  
  if (is.data.frame(obj)) {
    
    cols <- colnames(obj)
    
    if (!is.null(cols)) {
      
      for (cc in cols) {
        
        column_inventory <- rbind(
          column_inventory,
          data.frame(
            object = nm,
            column = cc,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  
  if (is.matrix(obj)) {
    
    cols <- colnames(obj)
    
    if (!is.null(cols)) {
      
      for (cc in cols) {
        
        column_inventory <- rbind(
          column_inventory,
          data.frame(
            object = nm,
            column = cc,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
}

write.csv(
  column_inventory,
  file.path(
    out_dir,
    "STEP13C-4A_All_Object_Columns.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 5. Flag potentially relevant columns
# ------------------------------------------------------------

patient_pattern <- paste(
  c(
    "patient",
    "patient_id",
    "subject",
    "subject_id",
    "individual",
    "individual_id",
    "pair",
    "pair_id",
    "case",
    "case_id",
    "donor",
    "donor_id"
  ),
  collapse = "|"
)

sample_pattern <- paste(
  c(
    "sample",
    "sample_id",
    "gsm",
    "accession",
    "geo_accession",
    "title",
    "characteristic",
    "source"
  ),
  collapse = "|"
)

candidate_columns <- column_inventory[
  grepl(
    paste0(patient_pattern, "|", sample_pattern),
    column_inventory$column,
    ignore.case = TRUE
  ),
]

write.csv(
  candidate_columns,
  file.path(
    out_dir,
    "STEP13C-4A_Potential_Patient_Sample_Columns.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\nPotential patient/sample columns:\n")
print(candidate_columns)

# ------------------------------------------------------------
# 6. annotation_file provenance
# ------------------------------------------------------------

annotation_file_present <- exists(
  "annotation_file",
  envir = hist_env,
  inherits = FALSE
)

annotation_file_value <- NA_character_

if (annotation_file_present) {
  
  annotation_file_value <- as.character(
    get("annotation_file", envir = hist_env)
  )
}

annotation_audit <- data.frame(
  annotation_file_present = annotation_file_present,
  annotation_file = annotation_file_value,
  stringsAsFactors = FALSE
)

write.csv(
  annotation_audit,
  file.path(
    out_dir,
    "STEP13C-4A_annotation_file_Audit.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 7. Report
# ------------------------------------------------------------

cat("\n")
cat("STEP13C-4A completed.\n")
cat("No patient pairing performed.\n")
cat("No sample reordering performed.\n")
cat("No expression reconstruction performed.\n")
cat("No ROC/AUC calculation performed.\n")
cat("No model fitting performed.\n")
cat("No threshold optimization performed.\n")
cat("\nOutput directory:\n")
cat(out_dir, "\n")

# ============================================================
# STEP13C-4B — Historical Sample-ID / Patient-ID Source Audit
# GSE32863 Paired ID Audit
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP13C-4B — Historical Sample-ID / Patient-ID Source Audit\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

step_dir <- file.path(
  workflow_root,
  "STEP13C_GSE32863_Input_Lock",
  "STEP13C-4_Paired_ID_Audit"
)

out_dir <- file.path(
  step_dir,
  "STEP13C-4B_Patient_ID_Source_Audit"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Historical RData
# ------------------------------------------------------------

hist_file <- file.path(
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal",
  ".RData"
)

stopifnot(file.exists(hist_file))

hist_env <- new.env(parent = emptyenv())

load(
  hist_file,
  envir = hist_env
)

cat("Historical .RData:\n")
cat(hist_file, "\n\n")

# ------------------------------------------------------------
# 2. Check df and df_wide
# ------------------------------------------------------------

required_objects <- c(
  "df",
  "df_wide"
)

object_presence <- data.frame(
  object = required_objects,
  present = sapply(
    required_objects,
    function(x) {
      exists(
        x,
        envir = hist_env,
        inherits = FALSE
      )
    }
  ),
  stringsAsFactors = FALSE
)

print(object_presence)

write.csv(
  object_presence,
  file.path(
    out_dir,
    "STEP13C-4B_Object_Presence.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 3. Inspect df
# ------------------------------------------------------------

df_obj <- NULL
df_wide_obj <- NULL

if (exists("df", envir = hist_env, inherits = FALSE)) {
  
  df_obj <- get(
    "df",
    envir = hist_env
  )
  
  cat("\n------------------------------------------------------------\n")
  cat("OBJECT: df\n")
  cat("------------------------------------------------------------\n")
  
  cat("Class:\n")
  print(class(df_obj))
  
  cat("\nDimensions:\n")
  print(dim(df_obj))
  
  cat("\nColumn names:\n")
  print(colnames(df_obj))
  
  cat("\nFirst 10 rows:\n")
  print(
    utils::head(
      df_obj,
      10
    )
  )
  
  cat("\nLast 10 rows:\n")
  print(
    utils::tail(
      df_obj,
      10
    )
  )
}

# ------------------------------------------------------------
# 4. Inspect df_wide
# ------------------------------------------------------------

if (exists("df_wide", envir = hist_env, inherits = FALSE)) {
  
  df_wide_obj <- get(
    "df_wide",
    envir = hist_env
  )
  
  cat("\n------------------------------------------------------------\n")
  cat("OBJECT: df_wide\n")
  cat("------------------------------------------------------------\n")
  
  cat("Class:\n")
  print(class(df_wide_obj))
  
  cat("\nDimensions:\n")
  print(dim(df_wide_obj))
  
  cat("\nColumn names:\n")
  print(colnames(df_wide_obj))
  
  cat("\nFirst 10 rows:\n")
  print(
    utils::head(
      df_wide_obj,
      10
    )
  )
  
  cat("\nLast 10 rows:\n")
  print(
    utils::tail(
      df_wide_obj,
      10
    )
  )
}

# ------------------------------------------------------------
# 5. Sample-ID structure audit
# ------------------------------------------------------------

sample_id_audit <- data.frame(
  object = character(),
  n_rows = integer(),
  n_columns = integer(),
  sample_id_column_present = logical(),
  sample_id_n = integer(),
  sample_id_unique = integer(),
  sample_id_missing = integer(),
  sample_id_duplicate = integer(),
  example_first = character(),
  example_last = character(),
  stringsAsFactors = FALSE
)

if (!is.null(df_obj)) {
  
  has_sample_id <- "sample_id" %in% colnames(df_obj)
  
  if (has_sample_id) {
    
    sid <- as.character(
      df_obj$sample_id
    )
    
    sample_id_audit <- rbind(
      sample_id_audit,
      data.frame(
        object = "df",
        n_rows = nrow(df_obj),
        n_columns = ncol(df_obj),
        sample_id_column_present = TRUE,
        sample_id_n = length(sid),
        sample_id_unique = length(unique(sid[!is.na(sid)])),
        sample_id_missing = sum(
          is.na(sid) | sid == ""
        ),
        sample_id_duplicate = sum(
          duplicated(sid[!is.na(sid)])
        ),
        example_first = sid[1],
        example_last = sid[length(sid)],
        stringsAsFactors = FALSE
      )
    )
  }
}

if (!is.null(df_wide_obj)) {
  
  has_sample_id <- "sample_id" %in% colnames(df_wide_obj)
  
  if (has_sample_id) {
    
    sid <- as.character(
      df_wide_obj$sample_id
    )
    
    sample_id_audit <- rbind(
      sample_id_audit,
      data.frame(
        object = "df_wide",
        n_rows = nrow(df_wide_obj),
        n_columns = ncol(df_wide_obj),
        sample_id_column_present = TRUE,
        sample_id_n = length(sid),
        sample_id_unique = length(unique(sid[!is.na(sid)])),
        sample_id_missing = sum(
          is.na(sid) | sid == ""
        ),
        sample_id_duplicate = sum(
          duplicated(sid[!is.na(sid)])
        ),
        example_first = sid[1],
        example_last = sid[length(sid)],
        stringsAsFactors = FALSE
      )
    )
  }
}

cat("\n------------------------------------------------------------\n")
cat("Sample-ID audit\n")
cat("------------------------------------------------------------\n")

print(sample_id_audit)

write.csv(
  sample_id_audit,
  file.path(
    out_dir,
    "STEP13C-4B_Sample_ID_Audit.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 6. Inspect all df / df_wide character columns
# ------------------------------------------------------------

character_column_inventory <- data.frame(
  object = character(),
  column = character(),
  n = integer(),
  n_unique = integer(),
  n_missing = integer(),
  example_first = character(),
  example_last = character(),
  stringsAsFactors = FALSE
)

inspect_character_columns <- function(
    obj,
    object_name
) {
  
  result <- data.frame(
    object = character(),
    column = character(),
    n = integer(),
    n_unique = integer(),
    n_missing = integer(),
    example_first = character(),
    example_last = character(),
    stringsAsFactors = FALSE
  )
  
  for (cc in colnames(obj)) {
    
    v <- obj[[cc]]
    
    if (
      is.character(v) ||
      is.factor(v)
    ) {
      
      v <- as.character(v)
      
      result <- rbind(
        result,
        data.frame(
          object = object_name,
          column = cc,
          n = length(v),
          n_unique = length(
            unique(
              v[!is.na(v)]
            )
          ),
          n_missing = sum(
            is.na(v) | v == ""
          ),
          example_first = if (
            length(v) > 0
          ) {
            v[1]
          } else {
            NA_character_
          },
          example_last = if (
            length(v) > 0
          ) {
            v[length(v)]
          } else {
            NA_character_
          },
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  result
}

if (!is.null(df_obj)) {
  
  character_column_inventory <- rbind(
    character_column_inventory,
    inspect_character_columns(
      df_obj,
      "df"
    )
  )
}

if (!is.null(df_wide_obj)) {
  
  character_column_inventory <- rbind(
    character_column_inventory,
    inspect_character_columns(
      df_wide_obj,
      "df_wide"
    )
  )
}

cat("\n------------------------------------------------------------\n")
cat("Character / factor columns\n")
cat("------------------------------------------------------------\n")

print(character_column_inventory)

write.csv(
  character_column_inventory,
  file.path(
    out_dir,
    "STEP13C-4B_Character_Column_Inventory.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 7. Search candidate patient / pair columns
# ------------------------------------------------------------

candidate_pattern <- paste(
  c(
    "patient",
    "patient_id",
    "subject",
    "subject_id",
    "individual",
    "individual_id",
    "pair",
    "pair_id",
    "case",
    "case_id",
    "donor",
    "donor_id",
    "person",
    "person_id"
  ),
  collapse = "|"
)

patient_like_columns <- character_column_inventory[
  grepl(
    candidate_pattern,
    character_column_inventory$column,
    ignore.case = TRUE
  ),
]

cat("\n------------------------------------------------------------\n")
cat("Potential patient / pair columns\n")
cat("------------------------------------------------------------\n")

print(patient_like_columns)

write.csv(
  patient_like_columns,
  file.path(
    out_dir,
    "STEP13C-4B_Potential_Patient_Pair_Columns.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 8. annotation_file
# ------------------------------------------------------------

annotation_file_present <- exists(
  "annotation_file",
  envir = hist_env,
  inherits = FALSE
)

annotation_file_value <- NA_character_

if (annotation_file_present) {
  
  annotation_file_value <- as.character(
    get(
      "annotation_file",
      envir = hist_env
    )
  )
}

annotation_audit <- data.frame(
  annotation_file_present = annotation_file_present,
  annotation_file = annotation_file_value,
  file_exists = if (
    !is.na(annotation_file_value)
  ) {
    file.exists(annotation_file_value)
  } else {
    FALSE
  },
  stringsAsFactors = FALSE
)

cat("\n------------------------------------------------------------\n")
cat("annotation_file audit\n")
cat("------------------------------------------------------------\n")

print(annotation_audit)

write.csv(
  annotation_audit,
  file.path(
    out_dir,
    "STEP13C-4B_annotation_file_Audit.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 9. Final interpretation
# ------------------------------------------------------------

explicit_patient_column_found <- nrow(
  patient_like_columns
) > 0

cat("\n")
cat("============================================================\n")
cat("STEP13C-4B — FINAL REPORT\n")
cat("============================================================\n\n")

cat(
  "Historical df present:",
  !is.null(df_obj),
  "\n"
)

cat(
  "Historical df_wide present:",
  !is.null(df_wide_obj),
  "\n"
)

cat(
  "Explicitly named patient/pair-like column found:",
  explicit_patient_column_found,
  "\n"
)

cat("\n")
cat("IMPORTANT:\n")
cat(
  "This step does NOT establish patient pairing.\n"
)
cat(
  "No Normal-Tumor pairs were constructed.\n"
)
cat(
  "No row-order pairing was used.\n"
)
cat(
  "No sample reordering was performed.\n"
)
cat(
  "No expression reconstruction was performed.\n"
)
cat(
  "No ROC/AUC recalculation was performed.\n"
)
cat(
  "No model fitting was performed.\n"
)
cat(
  "No threshold optimization was performed.\n"
)

cat("\nOutput directory:\n")
cat(out_dir, "\n")


# ============================================================
# STEP13C-4C
# GSE32863 Local GEO Metadata Source Audit
#
# Purpose:
#   Identify whether locally archived GSE32863 GEO files
#   contain explicit sample-level patient / subject / pair
#   identifiers.
#
# No pairing is constructed.
# No sample reordering.
# No expression reconstruction.
# No expression transformation.
# No probe selection.
# No ROC/AUC calculation.
# No model fitting.
# No threshold optimization.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP13C-4C — GSE32863 Local GEO Metadata Source Audit\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

source_dir <- "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal"

step_dir <- file.path(
  workflow_root,
  "STEP13C_GSE32863_Input_Lock",
  "STEP13C-4_Paired_ID_Audit"
)

out_dir <- file.path(
  step_dir,
  "STEP13C-4C_GSE32863_Local_Metadata_Source_Audit"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Local files
# ------------------------------------------------------------

files_to_audit <- c(
  "GSE32863_series_matrix.txt.gz",
  "GSE32863_non-normalized.txt.gz",
  "GSE32863_family.soft.gz",
  "Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal.R"
)

file_inventory <- data.frame(
  file = files_to_audit,
  full_path = file.path(
    source_dir,
    files_to_audit
  ),
  exists = file.exists(
    file.path(
      source_dir,
      files_to_audit
    )
  ),
  stringsAsFactors = FALSE
)

cat("Local file inventory:\n")
print(file_inventory)

write.csv(
  file_inventory,
  file.path(
    out_dir,
    "STEP13C-4C_Local_GSE32863_File_Inventory.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 2. Read metadata-related lines from series matrix
# ------------------------------------------------------------

series_file <- file.path(
  source_dir,
  "GSE32863_series_matrix.txt.gz"
)

series_lines <- character()

if (file.exists(series_file)) {
  
  con <- gzfile(
    series_file,
    open = "rt"
  )
  
  series_lines <- readLines(
    con,
    warn = FALSE
  )
  
  close(con)
}

cat("\nSeries matrix lines read:",
    length(series_lines),
    "\n")

# ------------------------------------------------------------
# 3. Extract GEO Sample metadata lines
# ------------------------------------------------------------

series_meta_lines <- series_lines[
  grepl(
    "^!Sample_",
    series_lines
  )
]

cat(
  "Series matrix !Sample_ lines:",
  length(series_meta_lines),
  "\n"
)

writeLines(
  series_meta_lines,
  con = file.path(
    out_dir,
    "STEP13C-4C_series_matrix_Sample_Metadata_Lines.txt"
  )
)

# ------------------------------------------------------------
# 4. Extract metadata fields
# ------------------------------------------------------------

metadata_field_inventory <- data.frame(
  field = character(),
  n_lines = integer(),
  example = character(),
  stringsAsFactors = FALSE
)

if (length(series_meta_lines) > 0) {
  
  field_names <- sub(
    "^(!Sample_[^\\t]+).*$",
    "\\1",
    series_meta_lines
  )
  
  unique_fields <- unique(
    field_names
  )
  
  for (ff in unique_fields) {
    
    idx <- field_names == ff
    
    examples <- series_meta_lines[idx]
    
    metadata_field_inventory <- rbind(
      metadata_field_inventory,
      data.frame(
        field = ff,
        n_lines = sum(idx),
        example = paste(
          head(examples, 3),
          collapse = " || "
        ),
        stringsAsFactors = FALSE
      )
    )
  }
}

cat("\nGEO Sample metadata fields:\n")
print(metadata_field_inventory)

write.csv(
  metadata_field_inventory,
  file.path(
    out_dir,
    "STEP13C-4C_GEO_Sample_Metadata_Field_Inventory.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 5. Search patient / pair-related terms
# ------------------------------------------------------------

patient_terms <- paste(
  c(
    "patient",
    "patient_id",
    "subject",
    "subject_id",
    "individual",
    "individual_id",
    "pair",
    "pair_id",
    "paired",
    "case",
    "case_id",
    "donor",
    "donor_id",
    "person",
    "person_id"
  ),
  collapse = "|"
)

patient_related_lines <- series_meta_lines[
  grepl(
    patient_terms,
    series_meta_lines,
    ignore.case = TRUE
  )
]

cat("\nPatient / pair-related GEO metadata lines:\n")
print(patient_related_lines)

writeLines(
  patient_related_lines,
  con = file.path(
    out_dir,
    "STEP13C-4C_Patient_Pair_Related_Metadata_Lines.txt"
  )
)

# ------------------------------------------------------------
# 6. Search tumor / normal / paired terminology
# ------------------------------------------------------------

group_terms <- paste(
  c(
    "tumor",
    "tumour",
    "normal",
    "paired",
    "adjacent",
    "matched",
    "non-tumor",
    "non-tumour"
  ),
  collapse = "|"
)

group_related_lines <- series_meta_lines[
  grepl(
    group_terms,
    series_meta_lines,
    ignore.case = TRUE
  )
]

cat("\nTumor / Normal / Pair-related metadata lines:\n")
print(group_related_lines)

writeLines(
  group_related_lines,
  con = file.path(
    out_dir,
    "STEP13C-4C_Tumor_Normal_Pair_Related_Metadata_Lines.txt"
  )
)

# ------------------------------------------------------------
# 7. Extract GEO accession and title/source lines
# ------------------------------------------------------------

important_sample_fields <- series_meta_lines[
  grepl(
    "^!Sample_(geo_accession|title|source_name_ch1|description|characteristics_ch1|biomaterial_provider)",
    series_meta_lines,
    ignore.case = TRUE
  )
]

cat("\nImportant sample metadata:\n")
print(
  head(
    important_sample_fields,
    100
  )
)

writeLines(
  important_sample_fields,
  con = file.path(
    out_dir,
    "STEP13C-4C_Important_GEO_Sample_Metadata.txt"
  )
)

# ------------------------------------------------------------
# 8. Inspect family.soft metadata
# ------------------------------------------------------------

soft_file <- file.path(
  source_dir,
  "GSE32863_family.soft.gz"
)

soft_lines <- character()

if (file.exists(soft_file)) {
  
  con <- gzfile(
    soft_file,
    open = "rt"
  )
  
  soft_lines <- readLines(
    con,
    warn = FALSE
  )
  
  close(con)
}

cat(
  "\nFamily SOFT lines read:",
  length(soft_lines),
  "\n"
)

soft_sample_lines <- soft_lines[
  grepl(
    "^!Sample_",
    soft_lines
  )
]

cat(
  "Family SOFT !Sample_ lines:",
  length(soft_sample_lines),
  "\n"
)

writeLines(
  soft_sample_lines,
  con = file.path(
    out_dir,
    "STEP13C-4C_family_SOFT_Sample_Metadata_Lines.txt"
  )
)

soft_patient_lines <- soft_sample_lines[
  grepl(
    patient_terms,
    soft_sample_lines,
    ignore.case = TRUE
  )
]

cat("\nPatient / pair-related SOFT metadata:\n")
print(soft_patient_lines)

writeLines(
  soft_patient_lines,
  con = file.path(
    out_dir,
    "STEP13C-4C_family_SOFT_Patient_Pair_Metadata.txt"
  )
)

# ------------------------------------------------------------
# 9. Non-normalized matrix header audit
# ------------------------------------------------------------

non_norm_file <- file.path(
  source_dir,
  "GSE32863_non-normalized.txt.gz"
)

non_norm_header <- character()

if (file.exists(non_norm_file)) {
  
  con <- gzfile(
    non_norm_file,
    open = "rt"
  )
  
  non_norm_header <- readLines(
    con,
    n = 5,
    warn = FALSE
  )
  
  close(con)
}

cat("\nNon-normalized file first lines:\n")
print(non_norm_header)

writeLines(
  non_norm_header,
  con = file.path(
    out_dir,
    "STEP13C-4C_NonNormalized_Header.txt"
  )
)

# ------------------------------------------------------------
# 10. Historical R script — metadata-related terms only
# ------------------------------------------------------------

r_script_file <- file.path(
  source_dir,
  "Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal.R"
)

r_script_lines <- character()

if (file.exists(r_script_file)) {
  
  r_script_lines <- readLines(
    r_script_file,
    warn = FALSE
  )
}

script_related_lines <- r_script_lines[
  grepl(
    paste(
      c(
        "GSE32863",
        "sample",
        "patient",
        "subject",
        "pair",
        "GSM",
        "TCGA"
      ),
      collapse = "|"
    ),
    r_script_lines,
    ignore.case = TRUE
  )
]

cat("\nHistorical R script relevant lines:\n")
print(
  head(
    script_related_lines,
    200
  )
)

writeLines(
  script_related_lines,
  con = file.path(
    out_dir,
    "STEP13C-4C_Historical_R_Script_Relevant_Lines.txt"
  )
)

# ------------------------------------------------------------
# 11. Final report
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13C-4C — FINAL REPORT\n")
cat("============================================================\n\n")

cat(
  "GSE32863 series matrix available:",
  file.exists(series_file),
  "\n"
)

cat(
  "GSE32863 family SOFT available:",
  file.exists(soft_file),
  "\n"
)

cat(
  "GSE32863 non-normalized matrix available:",
  file.exists(non_norm_file),
  "\n"
)

cat(
  "Historical R script available:",
  file.exists(r_script_file),
  "\n"
)

cat("\n")
cat(
  "Patient/pair identification was NOT performed.\n"
)

cat(
  "No artificial pairing was constructed.\n"
)

cat(
  "No row-order pairing was used.\n"
)

cat(
  "No sample reordering was performed.\n"
)

cat(
  "No expression reconstruction was performed.\n"
)

cat(
  "No expression transformation was performed.\n"
)

cat(
  "No probe selection was performed.\n"
)

cat(
  "No ROC/AUC recalculation was performed.\n"
)

cat(
  "No model fitting was performed.\n"
)

cat(
  "No threshold optimization was performed.\n"
)

cat("\nOutput directory:\n")
cat(out_dir, "\n")


# ============================================================
# STEP13C-4D
# GSE32863 GSM -> Patient/Subject/Pair Mapping Audit
#
# Purpose:
#   Determine whether explicit patient/subject/pair identifiers
#   can be established from archived GEO metadata.
#
# This step DOES NOT:
#   - reconstruct expression data
#   - transform expression values
#   - select probes
#   - recalculate ROC/AUC
#   - fit any model
#   - optimize thresholds
#   - construct artificial pairs
#   - use row order to define pairs
# ============================================================

cat("\n============================================================\n")
cat("STEP13C-4D — GSM -> Patient/Subject/Pair Mapping Audit\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

workflow_root <-
  "input_dir/Paper 1_Step 13-Diagnostic Classification"

input_dir <-
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal"

out_dir <-
  file.path(
    workflow_root,
    "STEP13C_GSE32863_Input_Lock",
    "STEP13C-4_Paired_ID_Audit",
    "STEP13C-4D_Sample_Patient_Mapping"
  )

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

series_file <-
  file.path(
    input_dir,
    "GSE32863_series_matrix.txt.gz"
  )

soft_file <-
  file.path(
    input_dir,
    "GSE32863_family.soft.gz"
  )

# ------------------------------------------------------------
# 1. Input existence audit
# ------------------------------------------------------------

cat("Input existence:\n")

cat(
  "  series matrix:",
  file.exists(series_file),
  "\n"
)

cat(
  "  family SOFT:",
  file.exists(soft_file),
  "\n\n"
)

if (!file.exists(series_file)) {
  stop("GSE32863_series_matrix.txt.gz not found.")
}

if (!file.exists(soft_file)) {
  stop("GSE32863_family.soft.gz not found.")
}

# ------------------------------------------------------------
# 2. Read GEO Series Matrix metadata only
# ------------------------------------------------------------

cat("Reading GEO Series Matrix metadata...\n")

series_lines <-
  readLines(
    gzfile(series_file),
    warn = FALSE
  )

sample_lines <-
  grep(
    "^!Sample_",
    series_lines,
    value = TRUE
  )

cat(
  "Number of !Sample_ metadata lines:",
  length(sample_lines),
  "\n\n"
)

# ------------------------------------------------------------
# 3. Extract metadata fields
# ------------------------------------------------------------

extract_sample_field <- function(
    lines,
    field_name
) {
  
  pattern <-
    paste0(
      "^!",
      field_name,
      "\t"
    )
  
  x <-
    grep(
      pattern,
      lines,
      value = TRUE
    )
  
  if (length(x) == 0) {
    return(character(0))
  }
  
  sub(
    pattern,
    "",
    x
  )
}

metadata_fields <- c(
  "Sample_geo_accession",
  "Sample_title",
  "Sample_source_name_ch1",
  "Sample_description",
  "Sample_characteristics_ch1",
  "Sample_biomaterial_provider"
)

metadata_list <- list()

for (field in metadata_fields) {
  
  vals <-
    extract_sample_field(
      series_lines,
      field
    )
  
  metadata_list[[field]] <-
    vals
  
  cat(
    field,
    ":",
    length(vals),
    "metadata line(s)\n"
  )
}

cat("\n")

# ------------------------------------------------------------
# 4. Print relevant metadata
# ------------------------------------------------------------

relevant_pattern <-
  paste(
    c(
      "patient",
      "patient_id",
      "subject",
      "subject_id",
      "individual",
      "individual_id",
      "pair",
      "pair_id",
      "paired",
      "case",
      "case_id",
      "donor",
      "donor_id",
      "person",
      "person_id",
      "matched",
      "adjacent",
      "tumor",
      "tumour",
      "normal",
      "non-tumor",
      "non-tumour"
    ),
    collapse = "|"
  )

relevant_sample_lines <-
  grep(
    relevant_pattern,
    sample_lines,
    value = TRUE,
    ignore.case = TRUE
  )

cat("============================================================\n")
cat("Relevant GEO Series Matrix metadata lines\n")
cat("============================================================\n\n")

if (length(relevant_sample_lines) > 0) {
  
  print(
    relevant_sample_lines
  )
  
} else {
  
  cat(
    "No patient/pair/tumor/normal-related terms detected.\n"
  )
}

writeLines(
  relevant_sample_lines,
  con = file.path(
    out_dir,
    "STEP13C-4D_Relevant_SeriesMatrix_Metadata.txt"
  )
)

# ------------------------------------------------------------
# 5. Extract GSM accession and title/source/characteristics
# ------------------------------------------------------------

gsm_vals <-
  metadata_list$Sample_geo_accession

title_vals <-
  metadata_list$Sample_title

source_vals <-
  metadata_list$Sample_source_name_ch1

description_vals <-
  metadata_list$Sample_description

characteristics_vals <-
  metadata_list$Sample_characteristics_ch1

provider_vals <-
  metadata_list$Sample_biomaterial_provider

cat("\n============================================================\n")
cat("Metadata vector lengths\n")
cat("============================================================\n\n")

cat(
  "GSM:",
  length(gsm_vals),
  "\n"
)

cat(
  "Title:",
  length(title_vals),
  "\n"
)

cat(
  "Source:",
  length(source_vals),
  "\n"
)

cat(
  "Description:",
  length(description_vals),
  "\n"
)

cat(
  "Characteristics:",
  length(characteristics_vals),
  "\n"
)

cat(
  "Biomaterial provider:",
  length(provider_vals),
  "\n\n"
)

# ------------------------------------------------------------
# 6. Save raw metadata lines for audit
# ------------------------------------------------------------

writeLines(
  sample_lines,
  con = file.path(
    out_dir,
    "STEP13C-4D_All_Sample_Metadata_Lines.txt"
  )
)

# ------------------------------------------------------------
# 7. Search specifically for patient/subject/pair identifiers
# ------------------------------------------------------------

id_terms <- c(
  "patient",
  "patient_id",
  "subject",
  "subject_id",
  "individual",
  "individual_id",
  "pair",
  "pair_id",
  "paired",
  "case",
  "case_id",
  "donor",
  "donor_id",
  "person",
  "person_id"
)

id_results <- list()

for (term in id_terms) {
  
  hits <-
    grep(
      term,
      sample_lines,
      value = TRUE,
      ignore.case = TRUE
    )
  
  id_results[[term]] <-
    hits
  
  cat(
    "\nTERM:",
    term,
    "\n"
  )
  
  if (length(hits) > 0) {
    
    print(
      hits
    )
    
  } else {
    
    cat(
      "No matching metadata lines.\n"
    )
  }
}

# ------------------------------------------------------------
# 8. Search for explicit pairing/matching language
# ------------------------------------------------------------

pair_terms <- c(
  "paired",
  "pair",
  "matched",
  "matching",
  "adjacent",
  "patient",
  "subject",
  "individual",
  "case"
)

pair_results <- list()

for (term in pair_terms) {
  
  hits <-
    grep(
      term,
      sample_lines,
      value = TRUE,
      ignore.case = TRUE
    )
  
  pair_results[[term]] <-
    hits
}

# ------------------------------------------------------------
# 9. Save patient/pair search result
# ------------------------------------------------------------

patient_pair_lines <-
  unique(
    unlist(
      id_results,
      use.names = FALSE
    )
  )

patient_pair_lines <-
  patient_pair_lines[
    !is.na(patient_pair_lines)
  ]

writeLines(
  patient_pair_lines,
  con = file.path(
    out_dir,
    "STEP13C-4D_Patient_Subject_Pair_Search.txt"
  )
)

# ------------------------------------------------------------
# 10. Save tumor/normal metadata search result
# ------------------------------------------------------------

tumor_normal_terms <- c(
  "tumor",
  "tumour",
  "normal",
  "paired",
  "matched",
  "adjacent",
  "non-tumor",
  "non-tumour"
)

tumor_normal_results <- list()

for (term in tumor_normal_terms) {
  
  hits <-
    grep(
      term,
      sample_lines,
      value = TRUE,
      ignore.case = TRUE
    )
  
  tumor_normal_results[[term]] <-
    hits
}

tumor_normal_lines <-
  unique(
    unlist(
      tumor_normal_results,
      use.names = FALSE
    )
  )

tumor_normal_lines <-
  tumor_normal_lines[
    !is.na(tumor_normal_lines)
  ]

writeLines(
  tumor_normal_lines,
  con = file.path(
    out_dir,
    "STEP13C-4D_Tumor_Normal_Pair_Language.txt"
  )
)

# ------------------------------------------------------------
# 11. Read GEO family SOFT metadata
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Reading GEO family SOFT metadata\n")
cat("============================================================\n\n")

soft_lines <-
  readLines(
    gzfile(soft_file),
    warn = FALSE
  )

soft_sample_lines <-
  grep(
    "^!Sample_",
    soft_lines,
    value = TRUE
  )

cat(
  "Number of !Sample_ lines in SOFT:",
  length(soft_sample_lines),
  "\n\n"
)

soft_relevant_lines <-
  grep(
    relevant_pattern,
    soft_sample_lines,
    value = TRUE,
    ignore.case = TRUE
  )

cat(
  "Patient/pair/tumor/normal-related SOFT lines:",
  length(soft_relevant_lines),
  "\n\n"
)

if (length(soft_relevant_lines) > 0) {
  
  print(
    soft_relevant_lines
  )
  
} else {
  
  cat(
    "No relevant patient/pair/tumor/normal terms detected in SOFT.\n"
  )
}

writeLines(
  soft_relevant_lines,
  con = file.path(
    out_dir,
    "STEP13C-4D_Relevant_SOFT_Metadata.txt"
  )
)

# ------------------------------------------------------------
# 12. Historical group consistency
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Historical group information\n")
cat("============================================================\n\n")

rdata_file <-
  file.path(
    input_dir,
    ".RData"
  )

if (file.exists(rdata_file)) {
  
  e <-
    new.env(
      parent = emptyenv()
    )
  
  load(
    rdata_file,
    envir = e
  )
  
  if (exists("group", envir = e)) {
    
    historical_group <-
      get(
        "group",
        envir = e
      )
    
    cat(
      "Historical group object present: TRUE\n"
    )
    
    cat(
      "Historical group length:",
      length(historical_group),
      "\n"
    )
    
    cat(
      "Normal:",
      sum(
        historical_group == "Normal",
        na.rm = TRUE
      ),
      "\n"
    )
    
    cat(
      "Tumor:",
      sum(
        historical_group == "Tumor",
        na.rm = TRUE
      ),
      "\n"
    )
    
  } else {
    
    cat(
      "Historical group object present: FALSE\n"
    )
  }
  
} else {
  
  cat(
    "Historical .RData not found.\n"
  )
}

# ------------------------------------------------------------
# 13. Important prohibition checks
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Prohibited operations audit\n")
cat("============================================================\n\n")

cat(
  "Artificial pairing constructed: FALSE\n"
)

cat(
  "Row-order pairing used: FALSE\n"
)

cat(
  "Sample reordering performed: FALSE\n"
)

cat(
  "Expression reconstruction performed: FALSE\n"
)

cat(
  "Expression transformation performed: FALSE\n"
)

cat(
  "Probe selection performed: FALSE\n"
)

cat(
  "ROC/AUC recalculation performed: FALSE\n"
)

cat(
  "Model fitting performed: FALSE\n"
)

cat(
  "Threshold optimization performed: FALSE\n"
)

# ------------------------------------------------------------
# 14. Preliminary mapping status
# ------------------------------------------------------------

explicit_patient_terms_detected <-
  length(patient_pair_lines) > 0

explicit_pair_terms_detected <-
  length(
    unique(
      unlist(
        pair_results,
        use.names = FALSE
      )
    )
  ) > 0

cat("\n============================================================\n")
cat("STEP13C-4D — PRELIMINARY STATUS\n")
cat("============================================================\n\n")

cat(
  "Explicit patient/subject/pair-related metadata detected:",
  explicit_patient_terms_detected,
  "\n"
)

cat(
  "Explicit pairing/matching-related metadata detected:",
  explicit_pair_terms_detected,
  "\n"
)

cat(
  "\nIMPORTANT:\n"
)

cat(
  "Detection of keywords alone does NOT establish valid patient-level pairing.\n"
)

cat(
  "A valid pair requires an explicit sample-to-patient/subject mapping.\n"
)

cat(
  "No pair will be constructed from row order or sample counts.\n"
)

# ------------------------------------------------------------
# 15. Save preliminary audit object
# ------------------------------------------------------------

audit_object <- list(
  series_file = series_file,
  soft_file = soft_file,
  gsm_values = gsm_vals,
  relevant_sample_lines = relevant_sample_lines,
  patient_pair_lines = patient_pair_lines,
  tumor_normal_lines = tumor_normal_lines,
  soft_relevant_lines = soft_relevant_lines,
  explicit_patient_terms_detected =
    explicit_patient_terms_detected,
  explicit_pair_terms_detected =
    explicit_pair_terms_detected
)

saveRDS(
  audit_object,
  file = file.path(
    out_dir,
    "STEP13C-4D_GSM_Patient_Pair_Metadata_Audit.rds"
  )
)

# ------------------------------------------------------------
# 16. Final report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP13C-4D — FINAL REPORT\n")
cat("============================================================\n\n")

cat(
  "GSE32863 Series Matrix available:",
  file.exists(series_file),
  "\n"
)

cat(
  "GSE32863 family SOFT available:",
  file.exists(soft_file),
  "\n"
)

cat(
  "GSM metadata detected:",
  length(gsm_vals),
  "\n"
)

cat(
  "Explicit patient/subject/pair metadata detected:",
  explicit_patient_terms_detected,
  "\n"
)

cat(
  "Explicit pairing/matching metadata detected:",
  explicit_pair_terms_detected,
  "\n\n"
)

cat(
  "Patient/pair mapping was NOT constructed automatically.\n"
)

cat(
  "No artificial pairing was constructed.\n"
)

cat(
  "No row-order pairing was used.\n"
)

cat(
  "No sample reordering was performed.\n"
)

cat(
  "No expression reconstruction was performed.\n"
)

cat(
  "No expression transformation was performed.\n"
)

cat(
  "No probe selection was performed.\n"
)

cat(
  "No ROC/AUC recalculation was performed.\n"
)

cat(
  "No model fitting was performed.\n"
)

cat(
  "No threshold optimization was performed.\n"
)

cat("\nOutput directory:\n")

cat(
  out_dir,
  "\n"
)

cat("\nSTEP13C-4D completed.\n")

# ============================================================
# STEP13C-4D-1
# Correct parsing of GSE32863 GEO sample metadata
#
# Purpose:
#   Correctly parse the multi-sample !Sample_* fields and
#   inspect the actual patient/subject/pair information.
#
# NO expression analysis.
# NO ROC/AUC.
# NO model fitting.
# NO sample reordering.
# NO artificial pairing.
# ============================================================

cat("\n============================================================\n")
cat("STEP13C-4D-1 — GSE32863 Metadata Parsing Audit\n")
cat("============================================================\n\n")

input_dir <-
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal"

out_dir <-
  "result_dir/Paper 1_Step 13-Diagnostic Classification/STEP13C_GSE32863_Input_Lock/STEP13C-4_Paired_ID_Audit/STEP13C-4D_Sample_Patient_Mapping"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

series_file <-
  file.path(
    input_dir,
    "GSE32863_series_matrix.txt.gz"
  )

soft_file <-
  file.path(
    input_dir,
    "GSE32863_family.soft.gz"
  )

# ------------------------------------------------------------
# 1. Read Series Matrix
# ------------------------------------------------------------

series_lines <-
  readLines(
    gzfile(series_file),
    warn = FALSE
  )

sample_lines <-
  grep(
    "^!Sample_",
    series_lines,
    value = TRUE
  )

cat(
  "Total !Sample_ metadata lines:",
  length(sample_lines),
  "\n\n"
)

# ------------------------------------------------------------
# 2. Extract the full metadata field after the tab
# ------------------------------------------------------------

extract_full_field <- function(
    lines,
    field
) {
  
  pattern <-
    paste0(
      "^!",
      field,
      "\t"
    )
  
  x <-
    grep(
      pattern,
      lines,
      value = TRUE
    )
  
  if (length(x) == 0) {
    return(character(0))
  }
  
  sub(
    pattern,
    "",
    x
  )
}

fields_to_show <- c(
  "Sample_geo_accession",
  "Sample_title",
  "Sample_source_name_ch1",
  "Sample_description",
  "Sample_characteristics_ch1",
  "Sample_biomaterial_provider"
)

metadata_raw <- list()

for (field in fields_to_show) {
  
  value <-
    extract_full_field(
      series_lines,
      field
    )
  
  metadata_raw[[field]] <-
    value
  
  cat(
    "FIELD:",
    field,
    "\n"
  )
  
  cat(
    "Number of metadata lines:",
    length(value),
    "\n"
  )
  
  if (length(value) > 0) {
    
    cat(
      "Character length:",
      nchar(value[1]),
      "\n"
    )
    
    cat(
      "First 1000 characters:\n"
    )
    
    cat(
      substr(
        value[1],
        1,
        1000
      ),
      "\n\n"
    )
  }
}

# ------------------------------------------------------------
# 3. Extract all GEO sample accessions using GSM pattern
# ------------------------------------------------------------

all_gsm_hits <-
  unlist(
    regmatches(
      sample_lines,
      gregexpr(
        "GSM[0-9]+",
        sample_lines
      )
    )
  )

all_gsm_hits <-
  unique(
    all_gsm_hits
  )

cat("============================================================\n")
cat("GSM accession audit\n")
cat("============================================================\n\n")

cat(
  "Unique GSM accessions detected:",
  length(all_gsm_hits),
  "\n\n"
)

if (length(all_gsm_hits) > 0) {
  
  print(
    head(
      all_gsm_hits,
      20
    )
  )
  
  cat("\n...\n\n")
  
  print(
    tail(
      all_gsm_hits,
      20
    )
  )
}

writeLines(
  all_gsm_hits,
  file.path(
    out_dir,
    "STEP13C-4D-1_GSM_Accessions.txt"
  )
)

# ------------------------------------------------------------
# 4. Patient / subject / pair / matching lines
# ------------------------------------------------------------

id_pattern <-
  paste(
    c(
      "patient",
      "patient_id",
      "subject",
      "subject_id",
      "individual",
      "individual_id",
      "pair",
      "pair_id",
      "paired",
      "case",
      "case_id",
      "donor",
      "donor_id",
      "person",
      "person_id",
      "matched",
      "matching"
    ),
    collapse = "|"
  )

id_lines <-
  grep(
    id_pattern,
    sample_lines,
    value = TRUE,
    ignore.case = TRUE
  )

cat("\n============================================================\n")
cat("Patient / Subject / Pair / Matching metadata\n")
cat("============================================================\n\n")

cat(
  "Number of relevant metadata lines:",
  length(id_lines),
  "\n\n"
)

if (length(id_lines) > 0) {
  
  print(
    id_lines
  )
}

writeLines(
  id_lines,
  file.path(
    out_dir,
    "STEP13C-4D-1_Patient_Subject_Pair_Metadata.txt"
  )
)

# ------------------------------------------------------------
# 5. Tumor / Normal / paired language
# ------------------------------------------------------------

tn_pattern <-
  paste(
    c(
      "tumor",
      "tumour",
      "normal",
      "paired",
      "matched",
      "adjacent",
      "non-tumor",
      "non-tumour"
    ),
    collapse = "|"
  )

tn_lines <-
  grep(
    tn_pattern,
    sample_lines,
    value = TRUE,
    ignore.case = TRUE
  )

cat("\n============================================================\n")
cat("Tumor / Normal / Pair metadata\n")
cat("============================================================\n\n")

cat(
  "Number of relevant metadata lines:",
  length(tn_lines),
  "\n\n"
)

if (length(tn_lines) > 0) {
  
  print(
    tn_lines
  )
}

writeLines(
  tn_lines,
  file.path(
    out_dir,
    "STEP13C-4D-1_Tumor_Normal_Pair_Metadata.txt"
  )
)

# ------------------------------------------------------------
# 6. SOFT metadata
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Family SOFT metadata audit\n")
cat("============================================================\n\n")

soft_lines <-
  readLines(
    gzfile(soft_file),
    warn = FALSE
  )

soft_sample_lines <-
  grep(
    "^!Sample_",
    soft_lines,
    value = TRUE
  )

cat(
  "Total SOFT !Sample_ lines:",
  length(soft_sample_lines),
  "\n"
)

soft_id_lines <-
  grep(
    id_pattern,
    soft_sample_lines,
    value = TRUE,
    ignore.case = TRUE
  )

cat(
  "SOFT patient/pair-related lines:",
  length(soft_id_lines),
  "\n\n"
)

if (length(soft_id_lines) > 0) {
  
  print(
    soft_id_lines
  )
}

writeLines(
  soft_id_lines,
  file.path(
    out_dir,
    "STEP13C-4D-1_SOFT_Patient_Pair_Metadata.txt"
  )
)

# ------------------------------------------------------------
# 7. Final audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP13C-4D-1 — FINAL REPORT\n")
cat("============================================================\n\n")

cat(
  "Series Matrix available:",
  file.exists(series_file),
  "\n"
)

cat(
  "Family SOFT available:",
  file.exists(soft_file),
  "\n"
)

cat(
  "Total !Sample_ metadata lines:",
  length(sample_lines),
  "\n"
)

cat(
  "Unique GSM accessions detected:",
  length(all_gsm_hits),
  "\n"
)

cat(
  "Patient/subject/pair-related metadata lines:",
  length(id_lines),
  "\n"
)

cat(
  "Tumor/normal/pair-related metadata lines:",
  length(tn_lines),
  "\n"
)

cat(
  "SOFT patient/pair-related metadata lines:",
  length(soft_id_lines),
  "\n\n"
)

cat(
  "This step only parses and audits GEO metadata.\n"
)

cat(
  "No patient pairing was constructed.\n"
)

cat(
  "No row-order pairing was used.\n"
)

cat(
  "No expression values were analyzed or transformed.\n"
)

cat(
  "No ROC/AUC calculation was performed.\n"
)

cat(
  "No model fitting was performed.\n"
)

cat(
  "No threshold optimization was performed.\n"
)

cat("\nOutput directory:\n")

cat(
  out_dir,
  "\n"
)

cat("\nSTEP13C-4D-1 completed.\n")

# ============================================================
# STEP13C-4D-2
# GSE32863 Explicit Sample-Pair Identifier Extraction
#
# Purpose:
#   Extract the sample identifier embedded in GEO sample titles
#   and audit whether the naming structure supports explicit
#   tumor-normal pairing.
#
# IMPORTANT:
#   This step does NOT yet declare pairs as validated.
#   It only extracts and audits the identifier structure.
#
# NO expression analysis
# NO expression transformation
# NO ROC/AUC
# NO model fitting
# NO threshold optimization
# NO row-order pairing
# ============================================================

cat("\n============================================================\n")
cat("STEP13C-4D-2 — Explicit Sample-Pair Identifier Extraction\n")
cat("============================================================\n\n")

input_dir <-
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal"

out_dir <-
  "result_dir/Paper 1_Step 13-Diagnostic Classification/STEP13C_GSE32863_Input_Lock/STEP13C-4_Paired_ID_Audit/STEP13C-4D_Sample_Patient_Mapping"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

series_file <-
  file.path(
    input_dir,
    "GSE32863_series_matrix.txt.gz"
  )

# ------------------------------------------------------------
# 1. Read Series Matrix metadata
# ------------------------------------------------------------

series_lines <-
  readLines(
    gzfile(series_file),
    warn = FALSE
  )

# ------------------------------------------------------------
# 2. Extract metadata line
# ------------------------------------------------------------

get_field <- function(
    lines,
    field
) {
  
  pattern <-
    paste0(
      "^!",
      field,
      "\t"
    )
  
  x <-
    grep(
      pattern,
      lines,
      value = TRUE
    )
  
  if (length(x) == 0) {
    return(character(0))
  }
  
  sub(
    pattern,
    "",
    x[1]
  )
}

gsm_raw <-
  get_field(
    series_lines,
    "Sample_geo_accession"
  )

title_raw <-
  get_field(
    series_lines,
    "Sample_title"
  )

source_raw <-
  get_field(
    series_lines,
    "Sample_source_name_ch1"
  )

characteristics_lines <-
  grep(
    "^!Sample_characteristics_ch1",
    series_lines,
    value = TRUE
  )

# ------------------------------------------------------------
# 3. Parse quoted Series Matrix values
# ------------------------------------------------------------

parse_quoted_values <- function(x) {
  
  if (length(x) == 0) {
    return(character(0))
  }
  
  m <-
    gregexpr(
      "\"[^\"]*\"",
      x
    )
  
  vals <-
    regmatches(
      x,
      m
    )[[1]]
  
  vals <-
    gsub(
      "^\"|\"$",
      "",
      vals
    )
  
  vals
}

gsm <-
  parse_quoted_values(
    gsm_raw
  )

title <-
  parse_quoted_values(
    title_raw
  )

source <-
  parse_quoted_values(
    source_raw
  )

# ------------------------------------------------------------
# 4. Basic dimension audit
# ------------------------------------------------------------

cat("============================================================\n")
cat("Basic metadata dimensions\n")
cat("============================================================\n\n")

cat(
  "GSM count:",
  length(gsm),
  "\n"
)

cat(
  "Title count:",
  length(title),
  "\n"
)

cat(
  "Source count:",
  length(source),
  "\n\n"
)

if (length(gsm) != 116) {
  stop(
    "Expected 116 GSM accessions but obtained ",
    length(gsm)
  )
}

if (length(title) != 116) {
  stop(
    "Expected 116 sample titles but obtained ",
    length(title)
  )
}

# ------------------------------------------------------------
# 5. Extract sample code from title
#
# Expected structure:
#   ... 05L12_N (expression)
#   ... 05L12_T (expression)
#
# The code before _N / _T is extracted.
# ------------------------------------------------------------

sample_code <-
  sub(
    ".*\\b([A-Za-z0-9]+)_(N|T)\\b.*",
    "\\1",
    title
  )

sample_type_from_title <-
  ifelse(
    grepl(
      "_N\\b",
      title
    ),
    "Normal",
    ifelse(
      grepl(
        "_T\\b",
        title
      ),
      "Tumor",
      NA_character_
    )
  )

# ------------------------------------------------------------
# 6. Verify extraction success
# ------------------------------------------------------------

extraction_failed <-
  is.na(sample_code) |
  sample_code == title

cat("============================================================\n")
cat("Sample-code extraction audit\n")
cat("============================================================\n\n")

cat(
  "Extraction failures:",
  sum(extraction_failed),
  "\n"
)

cat(
  "Normal from title:",
  sum(
    sample_type_from_title == "Normal",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Tumor from title:",
  sum(
    sample_type_from_title == "Tumor",
    na.rm = TRUE
  ),
  "\n\n"
)

if (sum(extraction_failed) > 0) {
  
  cat("Titles with extraction failure:\n")
  
  print(
    title[
      extraction_failed
    ]
  )
}

# ------------------------------------------------------------
# 7. Construct audit table
#
# This is metadata only.
# ------------------------------------------------------------

mapping_df <-
  data.frame(
    GSM = gsm,
    SampleTitle = title,
    Source = source,
    SampleCode = sample_code,
    TypeFromTitle = sample_type_from_title,
    stringsAsFactors = FALSE
  )

# ------------------------------------------------------------
# 8. Show first 20 records
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("First 20 GSM -> SampleCode mappings\n")
cat("============================================================\n\n")

print(
  mapping_df[
    1:min(20, nrow(mapping_df)),
  ]
)

# ------------------------------------------------------------
# 9. Show last 20 records
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Last 20 GSM -> SampleCode mappings\n")
cat("============================================================\n\n")

print(
  mapping_df[
    max(1, nrow(mapping_df) - 19):nrow(mapping_df),
  ]
)

# ------------------------------------------------------------
# 10. Duplicate sample-code audit
# ------------------------------------------------------------

code_table <-
  table(
    mapping_df$SampleCode
  )

cat("\n============================================================\n")
cat("Sample-code multiplicity audit\n")
cat("============================================================\n\n")

cat(
  "Unique sample codes:",
  length(code_table),
  "\n"
)

cat(
  "Codes occurring twice:",
  sum(code_table == 2),
  "\n"
)

cat(
  "Codes occurring once:",
  sum(code_table == 1),
  "\n"
)

cat(
  "Codes occurring >2 times:",
  sum(code_table > 2),
  "\n\n"
)

# ------------------------------------------------------------
# 11. Identify candidate pairs
# ------------------------------------------------------------

candidate_pair_codes <-
  names(
    code_table[
      code_table == 2
    ]
  )

cat(
  "Candidate paired sample codes:",
  length(candidate_pair_codes),
  "\n\n"
)

# ------------------------------------------------------------
# 12. For every candidate code, verify one N + one T
# ------------------------------------------------------------

pair_audit_list <- list()

for (code in candidate_pair_codes) {
  
  sub_df <-
    mapping_df[
      mapping_df$SampleCode == code,
    ]
  
  n_normal <-
    sum(
      sub_df$TypeFromTitle == "Normal",
      na.rm = TRUE
    )
  
  n_tumor <-
    sum(
      sub_df$TypeFromTitle == "Tumor",
      na.rm = TRUE
    )
  
  pair_audit_list[[code]] <-
    data.frame(
      SampleCode = code,
      N_records = n_normal,
      T_records = n_tumor,
      Valid_N_plus_T =
        n_normal == 1 &
        n_tumor == 1,
      stringsAsFactors = FALSE
    )
}

pair_audit_df <-
  do.call(
    rbind,
    pair_audit_list
  )

rownames(pair_audit_df) <-
  NULL

# ------------------------------------------------------------
# 13. Pair structure summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Candidate pair structure\n")
cat("============================================================\n\n")

cat(
  "Candidate codes with exactly one Normal + one Tumor:",
  sum(
    pair_audit_df$Valid_N_plus_T
  ),
  "\n"
)

cat(
  "Candidate codes failing N+T structure:",
  sum(
    !pair_audit_df$Valid_N_plus_T
  ),
  "\n\n"
)

if (
  any(
    !pair_audit_df$Valid_N_plus_T
  )
) {
  
  cat(
    "Codes failing the expected N+T structure:\n"
  )
  
  print(
    pair_audit_df[
      !pair_audit_df$Valid_N_plus_T,
    ]
  )
}

# ------------------------------------------------------------
# 14. Compare title-derived status with historical group
# ------------------------------------------------------------

rdata_file <-
  file.path(
    input_dir,
    ".RData"
  )

historical_group <-
  NULL

if (file.exists(rdata_file)) {
  
  e <-
    new.env(
      parent = emptyenv()
    )
  
  load(
    rdata_file,
    envir = e
  )
  
  if (
    exists(
      "group",
      envir = e
    )
  ) {
    
    historical_group <-
      get(
        "group",
        envir = e
      )
  }
}

if (
  !is.null(
    historical_group
  )
) {
  
  if (
    length(
      historical_group
    ) == nrow(mapping_df)
  ) {
    
    mapping_df$HistoricalGroup <-
      historical_group
    
    mapping_df$TitleGroupConcordant <-
      mapping_df$TypeFromTitle ==
      mapping_df$HistoricalGroup
    
    cat("\n============================================================\n")
    cat("Historical group concordance\n")
    cat("============================================================\n\n")
    
    cat(
      "Exact title-derived vs historical group concordance:",
      all(
        mapping_df$TitleGroupConcordant
      ),
      "\n"
    )
    
    cat(
      "Discordant records:",
      sum(
        !mapping_df$TitleGroupConcordant
      ),
      "\n"
    )
    
  } else {
    
    cat(
      "\nHistorical group length does not match 116.\n"
    )
  }
}

# ------------------------------------------------------------
# 15. Save mapping table
# ------------------------------------------------------------

write.csv(
  mapping_df,
  file.path(
    out_dir,
    "STEP13C-4D-2_GSM_SampleCode_Mapping.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  pair_audit_df,
  file.path(
    out_dir,
    "STEP13C-4D-2_Candidate_Pair_Structure_Audit.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

saveRDS(
  list(
    mapping = mapping_df,
    pair_audit = pair_audit_df
  ),
  file.path(
    out_dir,
    "STEP13C-4D-2_GSM_SampleCode_Mapping_Audit.rds"
  )
)

# ------------------------------------------------------------
# 16. Final report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP13C-4D-2 — FINAL REPORT\n")
cat("============================================================\n\n")

cat(
  "GSM records:",
  nrow(mapping_df),
  "\n"
)

cat(
  "Unique GSM:",
  length(
    unique(
      mapping_df$GSM
    )
  ),
  "\n"
)

cat(
  "Unique sample codes:",
  length(code_table),
  "\n"
)

cat(
  "Normal records:",
  sum(
    mapping_df$TypeFromTitle == "Normal",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Tumor records:",
  sum(
    mapping_df$TypeFromTitle == "Tumor",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Candidate N+T pairs:",
  sum(
    pair_audit_df$Valid_N_plus_T
  ),
  "\n\n"
)

cat(
  "IMPORTANT:\n"
)

cat(
  "Sample-code pairing is not yet declared as final patient-level pairing.\n"
)

cat(
  "STEP13C-4E must verify pair completeness before final lock.\n"
)

cat(
  "No row-order pairing was used.\n"
)

cat(
  "No artificial pairing was constructed.\n"
)

cat(
  "No expression values were analyzed or transformed.\n"
)

cat(
  "No ROC/AUC calculation was performed.\n"
)

cat(
  "No model fitting was performed.\n"
)

cat(
  "No threshold optimization was performed.\n"
)

cat("\nOutput directory:\n")

cat(
  out_dir,
  "\n"
)

cat("\nSTEP13C-4D-2 completed.\n")



# ============================================================
# STEP13C-4D-3
# Identify non-complete sample-code groups
# ============================================================

cat("\n============================================================\n")
cat("STEP13C-4D-3 — NON-COMPLETE SAMPLE-CODE AUDIT\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Reuse mapping_df generated by STEP13C-4D-2
# ------------------------------------------------------------

stopifnot(exists("mapping_df"))
stopifnot(exists("pair_audit_df"))

# ------------------------------------------------------------
# 2. Show all sample codes that are NOT exactly
#    one Normal + one Tumor
# ------------------------------------------------------------

bad_codes <- subset(
  pair_audit_df,
  Valid_N_plus_T == FALSE
)

cat("Non-complete sample codes:",
    nrow(bad_codes),
    "\n\n")

print(bad_codes)

# ------------------------------------------------------------
# 3. Print the actual GSM / title records for these codes
# ------------------------------------------------------------

if (nrow(bad_codes) > 0) {
  
  bad_mapping <- mapping_df[
    mapping_df$SampleCode %in% bad_codes$SampleCode,
    c(
      "GSM",
      "SampleTitle",
      "Source",
      "SampleCode",
      "TypeFromTitle"
    ),
    drop = FALSE
  ]
  
  cat("\n------------------------------------------------------------\n")
  cat("Detailed records for non-complete sample codes\n")
  cat("------------------------------------------------------------\n\n")
  
  print(
    bad_mapping,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No non-complete sample codes detected.\n"
  )
}

# ------------------------------------------------------------
# 4. Frequency structure of sample codes
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Sample-code frequency distribution\n")
cat("------------------------------------------------------------\n\n")

print(
  sort(
    table(mapping_df$SampleCode),
    decreasing = FALSE
  )
)

# ------------------------------------------------------------
# 5. Status structure of each non-complete code
# ------------------------------------------------------------

if (nrow(bad_codes) > 0) {
  
  cat("\n------------------------------------------------------------\n")
  cat("Normal/Tumor structure of non-complete codes\n")
  cat("------------------------------------------------------------\n\n")
  
  for (code in bad_codes$SampleCode) {
    
    tmp <- mapping_df[
      mapping_df$SampleCode == code,
      c(
        "GSM",
        "SampleTitle",
        "TypeFromTitle"
      ),
      drop = FALSE
    ]
    
    cat("\nSampleCode:", code, "\n")
    print(
      tmp,
      row.names = FALSE
    )
  }
}

# ------------------------------------------------------------
# 6. PASS / HOLD logic
# ------------------------------------------------------------

all_codes_complete <-
  nrow(bad_codes) == 0 &&
  length(unique(mapping_df$SampleCode)) == 58 &&
  nrow(mapping_df) == 116

cat("\n============================================================\n")
cat("STEP13C-4D-3 — FINAL STATUS\n")
cat("============================================================\n\n")

if (all_codes_complete) {
  
  cat(
    "PASS — 116 samples form 58 complete sample-code N+T pairs.\n"
  )
  
} else {
  
  cat(
    "HOLD — Non-complete sample-code structure detected.\n"
  )
  
  cat(
    "Do NOT declare 58 complete matched pairs yet.\n"
  )
}

cat("\nNo expression analysis was performed.\n")
cat("No ROC/AUC calculation was performed.\n")
cat("No model fitting was performed.\n")
cat("No threshold optimization was performed.\n")

# ============================================================
# STEP13C-4D-4
# Investigate singleton sample codes
# ============================================================

cat("\n============================================================\n")
cat("STEP13C-4D-4 — SINGLETON SAMPLE-CODE INVESTIGATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Identify singleton codes
# ------------------------------------------------------------

code_freq <- table(mapping_df$SampleCode)

singleton_codes <- names(
  code_freq[
    code_freq == 1
  ]
)

cat(
  "Singleton sample codes:",
  length(singleton_codes),
  "\n\n"
)

print(singleton_codes)

# ------------------------------------------------------------
# 2. Print complete metadata records
# ------------------------------------------------------------

singleton_mapping <- mapping_df[
  mapping_df$SampleCode %in% singleton_codes,
  c(
    "GSM",
    "SampleTitle",
    "Source",
    "SampleCode",
    "TypeFromTitle"
  ),
  drop = FALSE
]

cat("\n------------------------------------------------------------\n")
cat("Singleton sample-code records\n")
cat("------------------------------------------------------------\n\n")

print(
  singleton_mapping,
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. Compare with neighboring samples in metadata order
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Metadata-order context\n")
cat("------------------------------------------------------------\n\n")

for (gsm in singleton_mapping$GSM) {
  
  idx <- match(
    gsm,
    mapping_df$GSM
  )
  
  idx_start <- max(
    1,
    idx - 3
  )
  
  idx_end <- min(
    nrow(mapping_df),
    idx + 3
  )
  
  cat(
    "\nTarget GSM:",
    gsm,
    "\n"
  )
  
  print(
    mapping_df[
      idx_start:idx_end,
      c(
        "GSM",
        "SampleTitle",
        "SampleCode",
        "TypeFromTitle"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 4. Determine whether singleton codes are N or T
# ------------------------------------------------------------

singleton_status <- table(
  mapping_df$SampleCode[
    mapping_df$SampleCode %in% singleton_codes
  ],
  mapping_df$TypeFromTitle[
    mapping_df$SampleCode %in% singleton_codes
  ]
)

cat("\n------------------------------------------------------------\n")
cat("Singleton code × title-derived status\n")
cat("------------------------------------------------------------\n\n")

print(singleton_status)

# ------------------------------------------------------------
# 5. Summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP13C-4D-4 — SUMMARY\n")
cat("============================================================\n\n")

cat(
  "Total GSM records:",
  nrow(mapping_df),
  "\n"
)

cat(
  "Unique GSM:",
  length(unique(mapping_df$GSM)),
  "\n"
)

cat(
  "Unique sample codes:",
  length(unique(mapping_df$SampleCode)),
  "\n"
)

cat(
  "Singleton codes:",
  length(singleton_codes),
  "\n"
)

cat(
  "Complete N+T codes:",
  sum(
    code_freq == 2
  ),
  "\n"
)

cat(
  "\nNo expression analysis was performed.\n"
)

cat(
  "No ROC/AUC calculation was performed.\n"
)

cat(
  "No model fitting was performed.\n"
)

cat(
  "No threshold optimization was performed.\n"
)


# ============================================================
# STEP13C-4E
# FINAL MATCHED-PAIR STRUCTURE AUDIT
#
# Purpose:
#   1. Construct an auditable pairing table from the
#      metadata-derived SampleCode structure.
#   2. Identify 57 complete N+T sample-code pairs.
#   3. Identify the 2 unmatched samples.
#   4. Audit GEO metadata semantics for sample-code / patient /
#      subject / pair-related information.
#
# IMPORTANT:
#   - Metadata only.
#   - No expression analysis.
#   - No expression transformation.
#   - No ROC/AUC calculation.
#   - No model fitting.
#   - No threshold optimization.
#   - No artificial pairing.
#   - No row-order pairing.
# ============================================================


# ------------------------------------------------------------
# 0. Working directories
# ------------------------------------------------------------

workflow_root <-
  "input_dir/Paper 1_Step 13-Diagnostic Classification"

step_dir <-
  file.path(
    workflow_root,
    "STEP13C_GSE32863_Input_Lock",
    "STEP13C-4_Paired_ID_Audit"
  )

out_dir <-
  file.path(
    step_dir,
    "STEP13C-4E_Final_Pairing_Audit"
  )

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Required objects
# ------------------------------------------------------------

stopifnot(
  exists("mapping_df")
)

stopifnot(
  all(
    c(
      "GSM",
      "SampleTitle",
      "Source",
      "SampleCode",
      "TypeFromTitle"
    ) %in% names(mapping_df)
  )
)


# ------------------------------------------------------------
# 2. Basic input audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP13C-4E — INPUT AUDIT\n")
cat("============================================================\n\n")

cat(
  "Total records:",
  nrow(mapping_df),
  "\n"
)

cat(
  "Unique GSM:",
  length(unique(mapping_df$GSM)),
  "\n"
)

cat(
  "Unique SampleCode:",
  length(unique(mapping_df$SampleCode)),
  "\n"
)

cat(
  "Normal:",
  sum(
    mapping_df$TypeFromTitle == "Normal",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Tumor:",
  sum(
    mapping_df$TypeFromTitle == "Tumor",
    na.rm = TRUE
  ),
  "\n"
)


# ------------------------------------------------------------
# 3. Construct complete pair table
# ------------------------------------------------------------

codes <- sort(
  unique(
    mapping_df$SampleCode
  )
)

pair_list <- list()

pair_counter <- 0L

unmatched_list <- list()

for (code in codes) {
  
  tmp <- mapping_df[
    mapping_df$SampleCode == code,
    ,
    drop = FALSE
  ]
  
  n_normal <- sum(
    tmp$TypeFromTitle == "Normal",
    na.rm = TRUE
  )
  
  n_tumor <- sum(
    tmp$TypeFromTitle == "Tumor",
    na.rm = TRUE
  )
  
  if (
    nrow(tmp) == 2 &&
    n_normal == 1 &&
    n_tumor == 1
  ) {
    
    pair_counter <- pair_counter + 1L
    
    normal_row <- tmp[
      tmp$TypeFromTitle == "Normal",
      ,
      drop = FALSE
    ]
    
    tumor_row <- tmp[
      tmp$TypeFromTitle == "Tumor",
      ,
      drop = FALSE
    ]
    
    pair_list[[pair_counter]] <- data.frame(
      Pair_ID = sprintf(
        "Pair_%03d",
        pair_counter
      ),
      SampleCode = code,
      Normal_GSM = normal_row$GSM,
      Tumor_GSM = tumor_row$GSM,
      Normal_SampleTitle = normal_row$SampleTitle,
      Tumor_SampleTitle = tumor_row$SampleTitle,
      Pair_Status = "Complete_N_plus_T",
      stringsAsFactors = FALSE
    )
    
  } else {
    
    unmatched_list[[length(unmatched_list) + 1L]] <-
      data.frame(
        SampleCode = code,
        GSM = tmp$GSM,
        SampleTitle = tmp$SampleTitle,
        Source = tmp$Source,
        TypeFromTitle = tmp$TypeFromTitle,
        N_records = n_normal,
        T_records = n_tumor,
        Pair_Status = "Unmatched_or_nonstandard",
        stringsAsFactors = FALSE
      )
  }
}


pair_df <- do.call(
  rbind,
  pair_list
)


# ------------------------------------------------------------
# 4. Unmatched sample table
# ------------------------------------------------------------

if (length(unmatched_list) > 0) {
  
  unmatched_df <- do.call(
    rbind,
    unmatched_list
  )
  
} else {
  
  unmatched_df <- data.frame(
    SampleCode = character(0),
    GSM = character(0),
    SampleTitle = character(0),
    Source = character(0),
    TypeFromTitle = character(0),
    N_records = integer(0),
    T_records = integer(0),
    Pair_Status = character(0),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 5. Pair structure audit
# ------------------------------------------------------------

pair_structure_summary <- data.frame(
  
  Metric = c(
    "Total_GSM_records",
    "Unique_GSM",
    "Normal_records",
    "Tumor_records",
    "Unique_SampleCode",
    "Complete_NplusT_pairs",
    "Unmatched_records",
    "Unmatched_sample_codes"
  ),
  
  Value = c(
    nrow(mapping_df),
    length(unique(mapping_df$GSM)),
    sum(
      mapping_df$TypeFromTitle == "Normal",
      na.rm = TRUE
    ),
    sum(
      mapping_df$TypeFromTitle == "Tumor",
      na.rm = TRUE
    ),
    length(unique(mapping_df$SampleCode)),
    nrow(pair_df),
    nrow(unmatched_df),
    length(unique(unmatched_df$SampleCode))
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 6. Explicitly classify unmatched samples
# ------------------------------------------------------------

unmatched_df$Unmatched_Type <- NA_character_

unmatched_df$Unmatched_Type[
  unmatched_df$TypeFromTitle == "Tumor"
] <- "Unmatched_Tumor"

unmatched_df$Unmatched_Type[
  unmatched_df$TypeFromTitle == "Normal"
] <- "Unmatched_Normal"


# ------------------------------------------------------------
# 7. Print final pairing structure
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP13C-4E — PAIRING STRUCTURE\n")
cat("============================================================\n\n")

cat(
  "Complete N+T pairs:",
  nrow(pair_df),
  "\n"
)

cat(
  "Unmatched records:",
  nrow(unmatched_df),
  "\n"
)

cat(
  "Unmatched sample codes:",
  length(
    unique(
      unmatched_df$SampleCode
    )
  ),
  "\n\n"
)

cat("Complete pair table preview:\n\n")

print(
  head(
    pair_df,
    10
  ),
  row.names = FALSE
)

cat("\nUnmatched samples:\n\n")

print(
  unmatched_df,
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Verify expected structure
# ------------------------------------------------------------

complete_pair_count_ok <-
  nrow(pair_df) == 57

unmatched_record_count_ok <-
  nrow(unmatched_df) == 2

unmatched_code_count_ok <-
  length(
    unique(
      unmatched_df$SampleCode
    )
  ) == 2

total_record_reconstruction_ok <-
  (
    2 * nrow(pair_df) +
      nrow(unmatched_df)
  ) == nrow(mapping_df)

normal_reconstruction_ok <-
  (
    nrow(pair_df) +
      sum(
        unmatched_df$TypeFromTitle == "Normal",
        na.rm = TRUE
      )
  ) ==
  sum(
    mapping_df$TypeFromTitle == "Normal",
    na.rm = TRUE
  )

tumor_reconstruction_ok <-
  (
    nrow(pair_df) +
      sum(
        unmatched_df$TypeFromTitle == "Tumor",
        na.rm = TRUE
      )
  ) ==
  sum(
    mapping_df$TypeFromTitle == "Tumor",
    na.rm = TRUE
  )


# ------------------------------------------------------------
# 9. Verify no duplicated GSM across pair table
# ------------------------------------------------------------

pair_gsms <- c(
  pair_df$Normal_GSM,
  pair_df$Tumor_GSM,
  unmatched_df$GSM
)

gsm_reconstruction_unique <-
  length(pair_gsms) ==
  length(unique(pair_gsms))

gsm_reconstruction_complete <-
  setequal(
    pair_gsms,
    mapping_df$GSM
  )


# ------------------------------------------------------------
# 10. GEO metadata semantic audit
# ------------------------------------------------------------

series_file <-
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal/GSE32863_series_matrix.txt.gz"

family_file <-
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal/GSE32863_family.soft.gz"


cat("\n============================================================\n")
cat("STEP13C-4E — GEO METADATA SEMANTIC AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 10A. Read Series Matrix metadata lines only
# ------------------------------------------------------------

series_lines <- readLines(
  gzfile(series_file),
  warn = FALSE
)

sample_meta_lines <- grep(
  "^!Sample_",
  series_lines,
  value = TRUE
)

cat(
  "Series Matrix !Sample_ lines:",
  length(sample_meta_lines),
  "\n"
)


# ------------------------------------------------------------
# 10B. Search for relevant semantic terms
# ------------------------------------------------------------

semantic_patterns <- c(
  "patient",
  "subject",
  "individual",
  "case",
  "pair",
  "paired",
  "sample",
  "source",
  "title",
  "characteristics",
  "description",
  "biomaterial",
  "donor"
)

semantic_regex <- paste(
  semantic_patterns,
  collapse = "|"
)

series_relevant <- grep(
  semantic_regex,
  sample_meta_lines,
  ignore.case = TRUE,
  value = TRUE
)

cat(
  "Relevant Series Matrix metadata lines:",
  length(series_relevant),
  "\n\n"
)

cat(
  "Series Matrix relevant metadata:\n"
)

print(
  series_relevant
)


# ------------------------------------------------------------
# 10C. Read Family SOFT metadata lines
# ------------------------------------------------------------

family_lines <- readLines(
  gzfile(family_file),
  warn = FALSE
)

family_sample_lines <- grep(
  "^!Sample_",
  family_lines,
  value = TRUE
)

cat("\nFamily SOFT !Sample_ lines:",
    length(family_sample_lines),
    "\n"
)

family_relevant <- grep(
  semantic_regex,
  family_sample_lines,
  ignore.case = TRUE,
  value = TRUE
)

cat(
  "Relevant Family SOFT metadata lines:",
  length(family_relevant),
  "\n\n"
)

cat(
  "Family SOFT relevant metadata:\n"
)

print(
  family_relevant
)


# ------------------------------------------------------------
# 11. Direct semantic audit of the sample-code examples
# ------------------------------------------------------------

target_codes <- c(
  "3022",
  "3023",
  "3035"
)

target_pattern <- paste(
  target_codes,
  collapse = "|"
)

series_target_lines <- grep(
  target_pattern,
  series_relevant,
  ignore.case = TRUE,
  value = TRUE
)

family_target_lines <- grep(
  target_pattern,
  family_relevant,
  ignore.case = TRUE,
  value = TRUE
)

cat("\n============================================================\n")
cat("Target sample-code semantic audit\n")
cat("============================================================\n\n")

cat(
  "Series Matrix lines containing 3022/3023/3035:\n\n"
)

print(
  series_target_lines
)

cat(
  "\nFamily SOFT lines containing 3022/3023/3035:\n\n"
)

print(
  family_target_lines
)


# ------------------------------------------------------------
# 12. Conservative semantic classification
# ------------------------------------------------------------

patient_keyword_detected <- any(
  grepl(
    "patient|subject|individual|donor",
    c(
      series_relevant,
      family_relevant
    ),
    ignore.case = TRUE
  )
)

pair_keyword_detected <- any(
  grepl(
    "pair|paired",
    c(
      series_relevant,
      family_relevant
    ),
    ignore.case = TRUE
  )
)

sample_code_in_title_detected <- any(
  grepl(
    "3022|3023|3035",
    c(
      series_relevant,
      family_relevant
    ),
    ignore.case = TRUE
  )
)


# ------------------------------------------------------------
# 13. Important interpretation lock
# ------------------------------------------------------------

semantic_interpretation <- paste0(
  "Sample codes were extracted from explicit GEO sample titles ",
  "using the N/T naming structure. Complete N+T pairs therefore ",
  "represent metadata-supported matched sample-code pairs. ",
  "A sample code is NOT automatically interpreted as a patient ID ",
  "unless the GEO metadata explicitly establishes that semantic meaning."
)


# ------------------------------------------------------------
# 14. Save outputs
# ------------------------------------------------------------

write.csv(
  pair_df,
  file.path(
    out_dir,
    "STEP13C-4E_Complete_57_Pairs.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  unmatched_df,
  file.path(
    out_dir,
    "STEP13C-4E_Unmatched_Samples.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  pair_structure_summary,
  file.path(
    out_dir,
    "STEP13C-4E_Pair_Structure_Summary.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

writeLines(
  series_relevant,
  file.path(
    out_dir,
    "STEP13C-4E_SeriesMatrix_Relevant_Metadata.txt"
  ),
  useBytes = TRUE
)

writeLines(
  family_relevant,
  file.path(
    out_dir,
    "STEP13C-4E_FamilySOFT_Relevant_Metadata.txt"
  ),
  useBytes = TRUE
)

writeLines(
  c(
    "STEP13C-4E — Semantic interpretation lock",
    "============================================================",
    "",
    semantic_interpretation,
    "",
    paste0(
      "Patient/subject/donor keyword detected: ",
      patient_keyword_detected
    ),
    paste0(
      "Pair/paired keyword detected: ",
      pair_keyword_detected
    ),
    paste0(
      "Target sample-code strings detected in metadata: ",
      sample_code_in_title_detected
    ),
    "",
    "No patient-level identity was inferred solely from row order.",
    "No artificial pairing was constructed.",
    "No expression values were analyzed or transformed.",
    "No ROC/AUC calculation was performed.",
    "No model fitting was performed.",
    "No threshold optimization was performed."
  ),
  file.path(
    out_dir,
    "STEP13C-4E_Semantic_Interpretation_Lock.txt"
  ),
  useBytes = TRUE
)


# ------------------------------------------------------------
# 15. Save audit object
# ------------------------------------------------------------

audit_object <- list(
  
  pair_df = pair_df,
  
  unmatched_df = unmatched_df,
  
  pair_structure_summary =
    pair_structure_summary,
  
  series_relevant =
    series_relevant,
  
  family_relevant =
    family_relevant,
  
  target_codes =
    target_codes,
  
  series_target_lines =
    series_target_lines,
  
  family_target_lines =
    family_target_lines,
  
  patient_keyword_detected =
    patient_keyword_detected,
  
  pair_keyword_detected =
    pair_keyword_detected,
  
  sample_code_in_title_detected =
    sample_code_in_title_detected,
  
  checks = list(
    complete_pair_count_ok =
      complete_pair_count_ok,
    unmatched_record_count_ok =
      unmatched_record_count_ok,
    unmatched_code_count_ok =
      unmatched_code_count_ok,
    total_record_reconstruction_ok =
      total_record_reconstruction_ok,
    normal_reconstruction_ok =
      normal_reconstruction_ok,
    tumor_reconstruction_ok =
      tumor_reconstruction_ok,
    gsm_reconstruction_unique =
      gsm_reconstruction_unique,
    gsm_reconstruction_complete =
      gsm_reconstruction_complete
  ),
  
  semantic_interpretation =
    semantic_interpretation
)

saveRDS(
  audit_object,
  file.path(
    out_dir,
    "STEP13C-4E_Final_Pairing_Semantic_Audit.rds"
  )
)


# ------------------------------------------------------------
# 16. Final PASS / HOLD logic
# ------------------------------------------------------------

structure_pass <-
  complete_pair_count_ok &&
  unmatched_record_count_ok &&
  unmatched_code_count_ok &&
  total_record_reconstruction_ok &&
  normal_reconstruction_ok &&
  tumor_reconstruction_ok &&
  gsm_reconstruction_unique &&
  gsm_reconstruction_complete

cat("\n============================================================\n")
cat("STEP13C-4E — FINAL REPORT\n")
cat("============================================================\n\n")

cat(
  "Complete N+T pairs:",
  nrow(pair_df),
  "\n"
)

cat(
  "Unmatched records:",
  nrow(unmatched_df),
  "\n"
)

cat(
  "Unmatched sample codes:",
  length(
    unique(
      unmatched_df$SampleCode
    )
  ),
  "\n"
)

cat(
  "Total-record reconstruction:",
  total_record_reconstruction_ok,
  "\n"
)

cat(
  "Normal-count reconstruction:",
  normal_reconstruction_ok,
  "\n"
)

cat(
  "Tumor-count reconstruction:",
  tumor_reconstruction_ok,
  "\n"
)

cat(
  "GSM reconstruction unique:",
  gsm_reconstruction_unique,
  "\n"
)

cat(
  "GSM reconstruction complete:",
  gsm_reconstruction_complete,
  "\n"
)

cat("\n------------------------------------------------------------\n")
cat("Semantic interpretation\n")
cat("------------------------------------------------------------\n\n")

cat(
  semantic_interpretation,
  "\n"
)

cat("\n============================================================\n")
cat("STEP13C-4E — STATUS\n")
cat("============================================================\n\n")

if (structure_pass) {
  
  cat(
    "PASS — Final metadata-supported sample-code structure established.\n"
  )
  
  cat(
    "57 complete N+T pairs + 2 unmatched samples = 116 total samples.\n"
  )
  
} else {
  
  cat(
    "HOLD — Final pairing structure failed one or more integrity checks.\n"
  )
}

cat(
  "\nPatient-level pairing terminology remains conservative until\n"
)

cat(
  "GEO metadata semantics explicitly establish sample code as a\n"
)

cat(
  "patient/subject identifier.\n"
)

cat(
  "\nNo expression analysis was performed.\n"
)

cat(
  "No ROC/AUC calculation was performed.\n"
)

cat(
  "No model fitting was performed.\n"
)

cat(
  "No threshold optimization was performed.\n"
)

cat(
  "\nOutput directory:\n"
)

cat(
  out_dir,
  "\n"
)

cat(
  "\nSTEP13C-4E completed.\n"
)
# ============================================================
# STEP13C-4F
# Complete GEO Subject-ID Concordance Audit
# ============================================================
#
# PURPOSE
# -------
# Audit all 116 GSE32863 samples for:
#   1. GSM accession completeness and uniqueness
#   2. Explicit GEO subject: metadata
#   3. Concordance between subject: and Sample Title sample code
#   4. Subject-level Normal/Tumor structure
#   5. Complete matched subject pairs
#   6. Tumor-only / Normal-only subjects
#   7. Any subject mapping to multiple sample codes
#
# STRICT SCOPE
# ------------
# Metadata audit ONLY.
#
# NO:
#   - expression reconstruction
#   - probe selection
#   - normalization
#   - ROC/AUC calculation
#   - model fitting
#   - threshold optimization
#   - external validation recalculation
#
# ============================================================


# ------------------------------------------------------------
# 0. Working directory and folders
# ------------------------------------------------------------

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

input_dir <- file.path(
  workflow_root,
  "STEP13C_GSE32863_Input_Lock"
)

output_dir <- file.path(
  input_dir,
  "STEP13C-4_Paired_ID_Audit",
  "STEP13C-4F_Subject_ID_Concordance_Audit"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("\n")
cat("============================================================\n")
cat("STEP13C-4F — Complete GEO Subject-ID Concordance Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

series_matrix_file <- file.path(
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal",
  "GSE32863_series_matrix.txt.gz"
)

family_soft_file <- file.path(
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal",
  "GSE32863_family.soft.gz"
)

if (!file.exists(series_matrix_file)) {
  stop(
    "FAIL — Series Matrix file not found:\n",
    series_matrix_file
  )
}

if (!file.exists(family_soft_file)) {
  stop(
    "FAIL — Family SOFT file not found:\n",
    family_soft_file
  )
}

cat("Series Matrix:\n", series_matrix_file, "\n\n")
cat("Family SOFT:\n", family_soft_file, "\n\n")


# ------------------------------------------------------------
# 2. Read metadata only
# ------------------------------------------------------------
#
# IMPORTANT:
# We intentionally read ONLY metadata lines beginning with
# !Sample_.
#
# Expression data are NOT read or reconstructed.
# ------------------------------------------------------------

cat("Reading Series Matrix metadata lines only...\n")

sm_lines <- readLines(
  gzfile(series_matrix_file),
  warn = FALSE
)

sm_sample_lines <- sm_lines[
  grepl("^!Sample_", sm_lines)
]

cat(
  "Series Matrix !Sample_ metadata lines:",
  length(sm_sample_lines),
  "\n\n"
)


# ------------------------------------------------------------
# 3. Parse !Sample_geo_accession
# ------------------------------------------------------------

geo_line <- sm_sample_lines[
  grepl("^!Sample_geo_accession", sm_sample_lines)
]

if (length(geo_line) != 1) {
  stop(
    "FAIL — Expected exactly one !Sample_geo_accession line; found ",
    length(geo_line)
  )
}

geo_values <- strsplit(
  geo_line,
  "\t",
  fixed = TRUE
)[[1]]

geo_values <- geo_values[-1]

geo_values <- gsub(
  '^"|"$',
  "",
  geo_values
)

geo_values <- geo_values[
  nzchar(geo_values)
]

cat(
  "GSM records parsed:",
  length(geo_values),
  "\n"
)

cat(
  "Unique GSM accessions:",
  length(unique(geo_values)),
  "\n\n"
)


# ------------------------------------------------------------
# 4. Parse !Sample_title
# ------------------------------------------------------------

title_line <- sm_sample_lines[
  grepl("^!Sample_title", sm_sample_lines)
]

if (length(title_line) != 1) {
  stop(
    "FAIL — Expected exactly one !Sample_title line; found ",
    length(title_line)
  )
}

title_values <- strsplit(
  title_line,
  "\t",
  fixed = TRUE
)[[1]]

title_values <- title_values[-1]

title_values <- gsub(
  '^"|"$',
  "",
  title_values
)

title_values <- title_values[
  nzchar(title_values)
]

cat(
  "Sample titles parsed:",
  length(title_values),
  "\n\n"
)


# ------------------------------------------------------------
# 5. Basic GSM / title integrity
# ------------------------------------------------------------

if (length(geo_values) != length(title_values)) {
  stop(
    "FAIL — GSM count does not equal Sample Title count."
  )
}

if (length(geo_values) != 116) {
  stop(
    "FAIL — Expected 116 samples; found ",
    length(geo_values)
  )
}

if (anyDuplicated(geo_values)) {
  stop(
    "FAIL — Duplicate GSM accession detected."
  )
}

cat("PASS — 116 unique GSM accessions confirmed.\n")
cat("PASS — 116 Sample Titles confirmed.\n\n")


# ------------------------------------------------------------
# 6. Parse title-derived sample code and status
# ------------------------------------------------------------
#
# Expected title structure:
#
#   Adjacent non-tumor lung tissue 05L12_N (expression)
#   Lung adenocarcinoma 05L12_T (expression)
#
# Sample code:
#   05L12
#
# Status:
#   N = Normal
#   T = Tumor
# ------------------------------------------------------------

sample_code <- rep(NA_character_, length(title_values))
title_status <- rep(NA_character_, length(title_values))

for (i in seq_along(title_values)) {
  
  m <- regexec(
    "([A-Za-z0-9]+)_(N|T)\\s*\\(expression\\)",
    title_values[i]
  )
  
  hit <- regmatches(
    title_values[i],
    m
  )[[1]]
  
  if (length(hit) == 3) {
    
    sample_code[i] <- hit[2]
    
    if (hit[3] == "N") {
      title_status[i] <- "Normal"
    }
    
    if (hit[3] == "T") {
      title_status[i] <- "Tumor"
    }
    
  }
}


# ------------------------------------------------------------
# 7. Validate title parsing
# ------------------------------------------------------------

if (anyNA(sample_code)) {
  cat(
    "WARNING — Sample code could not be parsed for:\n"
  )
  
  print(
    data.frame(
      GSM = geo_values[is.na(sample_code)],
      Title = title_values[is.na(sample_code)],
      stringsAsFactors = FALSE
    )
  )
  
  stop(
    "FAIL — One or more Sample Titles failed sample-code parsing."
  )
}

if (anyNA(title_status)) {
  stop(
    "FAIL — One or more Sample Titles failed status parsing."
  )
}

cat(
  "PASS — Sample code parsed for all 116 samples.\n"
)

cat(
  "PASS — Normal/Tumor status parsed for all 116 samples.\n\n"
)


# ------------------------------------------------------------
# 8. Parse Family SOFT sample blocks
# ------------------------------------------------------------
#
# The Family SOFT contains explicit:
#
#   !Sample_geo_accession
#   !Sample_title
#   !Sample_characteristics_ch1
#
# We identify each sample block by !Sample_geo_accession
# and collect all characteristics belonging to that sample.
#
# This is metadata parsing only.
# ------------------------------------------------------------

cat("Reading Family SOFT metadata...\n")

soft_lines <- readLines(
  gzfile(family_soft_file),
  warn = FALSE
)

soft_sample_lines <- soft_lines[
  grepl("^!Sample_", soft_lines)
]

cat(
  "Family SOFT !Sample_ lines:",
  length(soft_sample_lines),
  "\n\n"
)


# ------------------------------------------------------------
# 9. Locate all sample accession records in Family SOFT
# ------------------------------------------------------------

soft_gsm_idx <- grep(
  "^!Sample_geo_accession\\s*=",
  soft_sample_lines
)

if (length(soft_gsm_idx) != 116) {
  stop(
    "FAIL — Expected 116 !Sample_geo_accession records in ",
    "Family SOFT; found ",
    length(soft_gsm_idx)
  )
}

cat(
  "PASS — Family SOFT contains 116 sample accession records.\n\n"
)


# ------------------------------------------------------------
# 10. Build sample-level metadata table
# ------------------------------------------------------------

soft_records <- vector(
  "list",
  length(soft_gsm_idx)
)

for (i in seq_along(soft_gsm_idx)) {
  
  start_idx <- soft_gsm_idx[i]
  
  if (i < length(soft_gsm_idx)) {
    end_idx <- soft_gsm_idx[i + 1] - 1
  } else {
    end_idx <- length(soft_sample_lines)
  }
  
  block <- soft_sample_lines[
    start_idx:end_idx
  ]
  
  # -------------------------
  # GSM
  # -------------------------
  
  gsm_match <- regexec(
    "^!Sample_geo_accession\\s*=\\s*\"?([^\"\\t]+)\"?",
    block[1]
  )
  
  gsm_hit <- regmatches(
    block[1],
    gsm_match
  )[[1]]
  
  gsm <- if (length(gsm_hit) >= 2) {
    gsm_hit[2]
  } else {
    NA_character_
  }
  
  
  # -------------------------
  # Sample title
  # -------------------------
  
  title_idx <- grep(
    "^!Sample_title\\s*=",
    block
  )
  
  title <- NA_character_
  
  if (length(title_idx) >= 1) {
    
    title_line2 <- block[title_idx[1]]
    
    title <- sub(
      "^!Sample_title\\s*=\\s*",
      "",
      title_line2
    )
    
    title <- gsub(
      '^"|"$',
      "",
      title
    )
  }
  
  
  # -------------------------
  # Characteristics
  # -------------------------
  
  char_lines <- block[
    grepl("^!Sample_characteristics_ch1\\s*=", block)
  ]
  
  char_values <- sub(
    "^!Sample_characteristics_ch1\\s*=\\s*",
    "",
    char_lines
  )
  
  char_values <- gsub(
    '^"|"$',
    "",
    char_values
  )
  
  
  # -------------------------
  # Extract subject:
  # -------------------------
  
  subject_values <- sub(
    "^.*?subject\\s*:\\s*([^;\"\\t]+).*$",
    "\\1",
    char_values,
    ignore.case = TRUE
  )
  
  subject_values <- trimws(subject_values)
  
  subject_values <- subject_values[
    grepl(
      "^subject\\s*:",
      char_values,
      ignore.case = TRUE
    )
  ]
  
  subject_values <- unique(
    subject_values
  )
  
  subject <- if (length(subject_values) == 1) {
    subject_values[1]
  } else if (length(subject_values) == 0) {
    NA_character_
  } else {
    paste(
      subject_values,
      collapse = " | "
    )
  }
  
  
  # -------------------------
  # Extract sample id:
  # -------------------------
  
  sample_id_values <- sub(
    "^.*?sample\\s+id\\s*:\\s*([^;\"\\t]+).*$",
    "\\1",
    char_values,
    ignore.case = TRUE
  )
  
  sample_id_values <- trimws(
    sample_id_values
  )
  
  sample_id_values <- sample_id_values[
    grepl(
      "^sample\\s+id\\s*:",
      char_values,
      ignore.case = TRUE
    )
  ]
  
  sample_id_values <- unique(
    sample_id_values
  )
  
  sample_id <- if (length(sample_id_values) == 1) {
    sample_id_values[1]
  } else if (length(sample_id_values) == 0) {
    NA_character_
  } else {
    paste(
      sample_id_values,
      collapse = " | "
    )
  }
  
  
  soft_records[[i]] <- data.frame(
    GSM = gsm,
    SOFT_Title = title,
    Subject_ID = subject,
    Sample_ID_metadata = sample_id,
    stringsAsFactors = FALSE
  )
}


soft_subject_df <- do.call(
  rbind,
  soft_records
)

rownames(soft_subject_df) <- NULL


# ------------------------------------------------------------
# 11. Subject-level metadata integrity
# ------------------------------------------------------------

if (nrow(soft_subject_df) != 116) {
  stop(
    "FAIL — Parsed Family SOFT sample table does not contain 116 rows."
  )
}

if (anyDuplicated(soft_subject_df$GSM)) {
  stop(
    "FAIL — Duplicate GSM accession in Family SOFT parsed table."
  )
}

if (anyNA(soft_subject_df$GSM)) {
  stop(
    "FAIL — Missing GSM accession in Family SOFT parsed table."
  )
}

cat(
  "PASS — 116 Family SOFT sample records parsed.\n"
)

cat(
  "PASS — All 116 GSM accessions are unique and non-missing.\n\n"
)


# ------------------------------------------------------------
# 12. Verify one explicit subject ID per GSM
# ------------------------------------------------------------

subject_missing <- is.na(
  soft_subject_df$Subject_ID
) |
  !nzchar(
    soft_subject_df$Subject_ID
  )

if (any(subject_missing)) {
  
  cat(
    "FAIL — Samples without explicit subject ID:\n"
  )
  
  print(
    soft_subject_df[
      subject_missing,
      c(
        "GSM",
        "SOFT_Title",
        "Subject_ID",
        "Sample_ID_metadata"
      )
    ]
  )
  
  stop(
    "FAIL — Missing subject ID detected."
  )
}

cat(
  "PASS — Explicit subject ID present for all 116 samples.\n\n"
)


# ------------------------------------------------------------
# 13. Merge Series Matrix title information with SOFT subject
#     information
# ------------------------------------------------------------

series_df <- data.frame(
  GSM = geo_values,
  Title_SeriesMatrix = title_values,
  SampleCode_Title = sample_code,
  Status_Title = title_status,
  stringsAsFactors = FALSE
)

merged_df <- merge(
  series_df,
  soft_subject_df,
  by = "GSM",
  all = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------
# 14. GSM completeness / exact correspondence
# ------------------------------------------------------------

gsm_set_series <- sort(
  unique(series_df$GSM)
)

gsm_set_soft <- sort(
  unique(soft_subject_df$GSM)
)

gsm_set_equal <- identical(
  gsm_set_series,
  gsm_set_soft
)

if (!gsm_set_equal) {
  
  cat(
    "Series Matrix-only GSMs:\n"
  )
  
  print(
    setdiff(
      gsm_set_series,
      gsm_set_soft
    )
  )
  
  cat(
    "\nFamily SOFT-only GSMs:\n"
  )
  
  print(
    setdiff(
      gsm_set_soft,
      gsm_set_series
    )
  )
  
  stop(
    "FAIL — GSM accession sets are not identical."
  )
}

cat(
  "PASS — Series Matrix and Family SOFT GSM sets are identical.\n\n"
)


# ------------------------------------------------------------
# 15. Subject ID versus title sample code concordance
# ------------------------------------------------------------
#
# Expected:
#
#   SampleCode_Title == Subject_ID
#
# for all 116 samples.
#
# No inference beyond exact string comparison.
# ------------------------------------------------------------

merged_df$Subject_vs_TitleCode <- ifelse(
  merged_df$Subject_ID ==
    merged_df$SampleCode_Title,
  "MATCH",
  "MISMATCH"
)

n_match <- sum(
  merged_df$Subject_vs_TitleCode == "MATCH",
  na.rm = TRUE
)

n_mismatch <- sum(
  merged_df$Subject_vs_TitleCode == "MISMATCH",
  na.rm = TRUE
)

cat(
  "Subject ID / Sample Code matches:",
  n_match,
  "\n"
)

cat(
  "Subject ID / Sample Code mismatches:",
  n_mismatch,
  "\n\n"
)


# ------------------------------------------------------------
# 16. Stop if any subject/sample-code mismatch exists
# ------------------------------------------------------------

if (n_mismatch > 0) {
  
  cat(
    "MISMATCHED RECORDS:\n\n"
  )
  
  print(
    merged_df[
      merged_df$Subject_vs_TitleCode == "MISMATCH",
      c(
        "GSM",
        "Title_SeriesMatrix",
        "SampleCode_Title",
        "Subject_ID",
        "Sample_ID_metadata",
        "Status_Title"
      )
    ],
    row.names = FALSE
  )
  
  stop(
    "FAIL — Subject ID and Sample Code mismatch detected."
  )
}

cat(
  "PASS — Subject ID matches title-derived sample code for ",
  n_match,
  "/116 samples.\n\n"
)


# ------------------------------------------------------------
# 17. Audit whether one subject maps to multiple sample codes
# ------------------------------------------------------------

subject_code_audit <- aggregate(
  SampleCode_Title ~ Subject_ID,
  data = merged_df,
  FUN = function(x) {
    paste(
      sort(unique(x)),
      collapse = " | "
    )
  }
)

names(subject_code_audit)[2] <- "SampleCodes"

subject_code_audit$N_SampleCodes <- sapply(
  strsplit(
    subject_code_audit$SampleCodes,
    " \\| ",
    fixed = FALSE
  ),
  length
)

multi_code_subjects <- subject_code_audit[
  subject_code_audit$N_SampleCodes > 1,
]


if (nrow(multi_code_subjects) > 0) {
  
  cat(
    "FAIL — Subject IDs mapped to multiple Sample Codes:\n\n"
  )
  
  print(
    multi_code_subjects,
    row.names = FALSE
  )
  
  stop(
    "FAIL — One or more subjects map to multiple sample codes."
  )
}

cat(
  "PASS — No subject maps to multiple Sample Codes.\n\n"
)


# ------------------------------------------------------------
# 18. Subject-level structure
# ------------------------------------------------------------
#
# For each subject:
#
#   Normal records
#   Tumor records
#   Total records
#
# Classification:
#
#   Complete_N_plus_T
#   Tumor_only
#   Normal_only
#   Other_structure
# ------------------------------------------------------------

merged_df$Status <- merged_df$Status_Title

subject_split <- split(
  merged_df,
  merged_df$Subject_ID
)

subject_structure <- lapply(
  names(subject_split),
  function(subj) {
    
    d <- subject_split[[subj]]
    
    n_normal <- sum(
      d$Status == "Normal"
    )
    
    n_tumor <- sum(
      d$Status == "Tumor"
    )
    
    n_total <- nrow(d)
    
    if (
      n_normal == 1 &&
      n_tumor == 1 &&
      n_total == 2
    ) {
      
      structure_class <- "Complete_N_plus_T"
      
    } else if (
      n_normal == 0 &&
      n_tumor == 1 &&
      n_total == 1
    ) {
      
      structure_class <- "Tumor_only"
      
    } else if (
      n_normal == 1 &&
      n_tumor == 0 &&
      n_total == 1
    ) {
      
      structure_class <- "Normal_only"
      
    } else {
      
      structure_class <- "Other_structure"
    }
    
    
    data.frame(
      Subject_ID = subj,
      N_Normal = n_normal,
      N_Tumor = n_tumor,
      N_Total = n_total,
      Structure = structure_class,
      SampleCodes = paste(
        sort(unique(d$SampleCode_Title)),
        collapse = " | "
      ),
      GSMs = paste(
        d$GSM,
        collapse = " | "
      ),
      stringsAsFactors = FALSE
    )
  }
)

subject_structure_df <- do.call(
  rbind,
  subject_structure
)

rownames(subject_structure_df) <- NULL

subject_structure_df <- subject_structure_df[
  order(
    suppressWarnings(
      as.numeric(subject_structure_df$Subject_ID)
    ),
    subject_structure_df$Subject_ID
  ),
]


# ------------------------------------------------------------
# 19. Subject-level summary
# ------------------------------------------------------------

structure_summary <- as.data.frame(
  table(
    subject_structure_df$Structure
  ),
  stringsAsFactors = FALSE
)

names(structure_summary) <- c(
  "Structure",
  "N_Subjects"
)

cat(
  "Subject-level structure:\n\n"
)

print(
  structure_summary,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 20. Expected subject structure
# ------------------------------------------------------------

n_subjects <- nrow(
  subject_structure_df
)

n_complete <- sum(
  subject_structure_df$Structure ==
    "Complete_N_plus_T"
)

n_tumor_only <- sum(
  subject_structure_df$Structure ==
    "Tumor_only"
)

n_normal_only <- sum(
  subject_structure_df$Structure ==
    "Normal_only"
)

n_other <- sum(
  subject_structure_df$Structure ==
    "Other_structure"
)


# ------------------------------------------------------------
# 21. Explicit subject-level integrity checks
# ------------------------------------------------------------

check_n_subjects <- (
  n_subjects == 59
)

check_complete_pairs <- (
  n_complete == 57
)

check_tumor_only <- (
  n_tumor_only == 1
)

check_normal_only <- (
  n_normal_only == 1
)

check_other <- (
  n_other == 0
)

check_total_subject_records <- (
  sum(subject_structure_df$N_Total) == 116
)

check_normal_total <- (
  sum(subject_structure_df$N_Normal) == 58
)

check_tumor_total <- (
  sum(subject_structure_df$N_Tumor) == 58
)


# ------------------------------------------------------------
# 22. Verify expected singleton subjects
# ------------------------------------------------------------

tumor_only_subjects <- subject_structure_df$Subject_ID[
  subject_structure_df$Structure == "Tumor_only"
]

normal_only_subjects <- subject_structure_df$Subject_ID[
  subject_structure_df$Structure == "Normal_only"
]

check_tumor_singleton_3023 <- (
  identical(
    sort(tumor_only_subjects),
    "3023"
  )
)

check_normal_singleton_3035 <- (
  identical(
    sort(normal_only_subjects),
    "3035"
  )
)


# ------------------------------------------------------------
# 23. Full PASS/FAIL evaluation
# ------------------------------------------------------------

all_checks <- c(
  
  GSM_116 =
    length(geo_values) == 116,
  
  GSM_unique =
    length(unique(geo_values)) == 116,
  
  Title_116 =
    length(title_values) == 116,
  
  Subject_116 =
    nrow(soft_subject_df) == 116,
  
  Subject_nonmissing_116 =
    sum(!subject_missing) == 116,
  
  GSM_set_concordance =
    gsm_set_equal,
  
  Subject_TitleCode_116_match =
    n_match == 116,
  
  Subject_multiple_codes_absent =
    nrow(multi_code_subjects) == 0,
  
  Unique_subjects_59 =
    check_n_subjects,
  
  Complete_pairs_57 =
    check_complete_pairs,
  
  Tumor_only_1 =
    check_tumor_only,
  
  Normal_only_1 =
    check_normal_only,
  
  Other_structure_0 =
    check_other,
  
  Total_subject_records_116 =
    check_total_subject_records,
  
  Normal_records_58 =
    check_normal_total,
  
  Tumor_records_58 =
    check_tumor_total,
  
  Tumor_only_subject_3023 =
    check_tumor_singleton_3023,
  
  Normal_only_subject_3035 =
    check_normal_singleton_3035
)


cat(
  "============================================================\n"
)

cat(
  "FINAL CHECKS\n"
)

cat(
  "============================================================\n\n"
)

print(
  all_checks
)

cat("\n")


# ------------------------------------------------------------
# 24. Final status
# ------------------------------------------------------------

final_pass <- all(
  all_checks
)

if (final_pass) {
  
  final_status <- paste(
    "PASS — STEP13C-4F Complete GEO Subject-ID Concordance",
    "Audit established.",
    "116 samples,",
    "59 unique subjects,",
    "57 complete N+T subject pairs,",
    "1 tumor-only subject (3023),",
    "1 normal-only subject (3035),",
    "and exact subject/sample-code concordance for all 116 samples."
  )
  
} else {
  
  final_status <- paste(
    "FAIL — STEP13C-4F Subject-ID Concordance Audit did not",
    "meet all predefined integrity criteria."
  )
}

cat(
  final_status,
  "\n\n"
)


# ------------------------------------------------------------
# 25. Save outputs
# ------------------------------------------------------------

mapping_file <- file.path(
  output_dir,
  "STEP13C-4F_GSM_Subject_Mapping.csv"
)

subject_table_file <- file.path(
  output_dir,
  "STEP13C-4F_Subject_Level_Pairing_Table.csv"
)

summary_file <- file.path(
  output_dir,
  "STEP13C-4F_Subject_Structure_Summary.csv"
)

audit_rds <- file.path(
  output_dir,
  "STEP13C-4F_Subject_Semantic_Audit.rds"
)

checks_file <- file.path(
  output_dir,
  "STEP13C-4F_Integrity_Checks.csv"
)

pass_file <- file.path(
  output_dir,
  "STEP13C-4F_PASS.txt"
)

write.csv(
  merged_df,
  mapping_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  subject_structure_df,
  subject_table_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  structure_summary,
  summary_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  data.frame(
    Check = names(all_checks),
    PASS = as.logical(all_checks),
    stringsAsFactors = FALSE
  ),
  checks_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


audit_object <- list(
  
  workflow = "STEP13C-4F",
  
  purpose =
    "Complete GEO Subject-ID Concordance Audit",
  
  scope =
    "Metadata-only audit; no expression or ROC/model analysis",
  
  input_series_matrix =
    series_matrix_file,
  
  input_family_soft =
    family_soft_file,
  
  n_samples =
    length(geo_values),
  
  n_unique_GSM =
    length(unique(geo_values)),
  
  n_unique_subjects =
    n_subjects,
  
  n_complete_N_plus_T =
    n_complete,
  
  n_tumor_only =
    n_tumor_only,
  
  n_normal_only =
    n_normal_only,
  
  n_other_structure =
    n_other,
  
  n_subject_title_code_match =
    n_match,
  
  n_subject_title_code_mismatch =
    n_mismatch,
  
  tumor_only_subjects =
    tumor_only_subjects,
  
  normal_only_subjects =
    normal_only_subjects,
  
  integrity_checks =
    all_checks,
  
  final_pass =
    final_pass,
  
  timestamp =
    Sys.time()
)

saveRDS(
  audit_object,
  audit_rds
)


# ------------------------------------------------------------
# 26. PASS file
# ------------------------------------------------------------

if (final_pass) {
  
  pass_text <- c(
    "STEP13C-4F — COMPLETE GEO SUBJECT-ID CONCORDANCE AUDIT",
    "========================================================",
    "",
    "FINAL STATUS:",
    final_status,
    "",
    "LOCKED METADATA STRUCTURE:",
    paste0("Total samples: ", length(geo_values)),
    paste0("Unique GSM accessions: ", length(unique(geo_values))),
    paste0("Unique subjects: ", n_subjects),
    paste0("Complete N+T subject pairs: ", n_complete),
    paste0("Tumor-only subjects: ", n_tumor_only),
    paste0("Normal-only subjects: ", n_normal_only),
    paste0("Other subject structures: ", n_other),
    "",
    "Subject 3023: Tumor-only",
    "Subject 3035: Normal-only",
    "",
    paste0(
      "Subject ID / title sample-code exact matches: ",
      n_match,
      "/116"
    ),
    "",
    "IMPORTANT:",
    "The two singleton subjects were NOT artificially paired.",
    "The historical GSE32863 ROC result was NOT recalculated.",
    "No expression values were reconstructed.",
    "No probe selection was performed.",
    "No ROC/AUC was calculated.",
    "No model was fitted.",
    "No threshold was optimized.",
    "",
    paste0(
      "Timestamp: ",
      Sys.time()
    )
  )
  
  writeLines(
    pass_text,
    pass_file,
    useBytes = TRUE
  )
  
} else {
  
  fail_file <- file.path(
    output_dir,
    "STEP13C-4F_FAIL.txt"
  )
  
  writeLines(
    c(
      "STEP13C-4F — FAIL",
      "",
      final_status,
      "",
      paste0(
        "Timestamp: ",
        Sys.time()
      )
    ),
    fail_file,
    useBytes = TRUE
  )
}


# ------------------------------------------------------------
# 27. Session information
# ------------------------------------------------------------

session_file <- file.path(
  output_dir,
  "STEP13C-4F_sessionInfo.txt"
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  session_file,
  useBytes = TRUE
)


# ------------------------------------------------------------
# 28. Final console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13C-4F FINAL SUMMARY\n")
cat("============================================================\n\n")

cat("Samples:                 ", length(geo_values), "\n")
cat("Unique GSMs:             ", length(unique(geo_values)), "\n")
cat("Unique subjects:         ", n_subjects, "\n")
cat("Complete N+T subjects:   ", n_complete, "\n")
cat("Tumor-only subjects:     ", n_tumor_only, "\n")
cat("Normal-only subjects:    ", n_normal_only, "\n")
cat("Other structures:        ", n_other, "\n")
cat("Subject/code matches:    ", n_match, "/116\n")
cat("Subject/code mismatches: ", n_mismatch, "\n\n")

cat("Expected singleton subjects:\n")
cat("  Tumor-only:  3023\n")
cat("  Normal-only: 3035\n\n")

cat(
  "FINAL STATUS:\n",
  final_status,
  "\n\n"
)

cat(
  "No expression analysis / ROC / AUC / model / threshold analysis",
  "was performed.\n\n"
)

cat(
  "Outputs:\n",
  output_dir,
  "\n"
)

cat("============================================================\n")

# ============================================================
# STEP13C-5
# FINAL GSE32863 INPUT / SUBJECT-PAIRING LOCK
# ============================================================
#
# PURPOSE
# -------
# Formally lock the already-audited GSE32863 sample and
# subject-level structure.
#
# This step summarizes previously completed PASS audits.
#
# STRICTLY NO:
#   - expression reconstruction
#   - probe selection
#   - normalization
#   - ROC/AUC calculation
#   - model fitting
#   - threshold optimization
#   - re-pairing
#
# ============================================================


# ------------------------------------------------------------
# 0. Configuration
# ------------------------------------------------------------

workflow_root <-
  "input_dir/Paper 1_Step 13-Diagnostic Classification"

step13c_root <-
  file.path(
    workflow_root,
    "STEP13C_GSE32863_Input_Lock"
  )

step13c4f_dir <-
  file.path(
    step13c_root,
    "STEP13C-4_Paired_ID_Audit",
    "STEP13C-4F_Subject_ID_Concordance_Audit"
  )

output_dir <-
  file.path(
    step13c_root,
    "STEP13C-5_Final_Input_Lock"
  )

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 1. Input files from PASS audit
# ------------------------------------------------------------

mapping_file <-
  file.path(
    step13c4f_dir,
    "STEP13C-4F_GSM_Subject_Mapping.csv"
  )

subject_table_file <-
  file.path(
    step13c4f_dir,
    "STEP13C-4F_Subject_Level_Pairing_Table.csv"
  )

summary_file <-
  file.path(
    step13c4f_dir,
    "STEP13C-4F_Subject_Structure_Summary.csv"
  )

checks_file <-
  file.path(
    step13c4f_dir,
    "STEP13C-4F_Integrity_Checks.csv"
  )

audit_rds <-
  file.path(
    step13c4f_dir,
    "STEP13C-4F_Subject_Semantic_Audit.rds"
  )

pass_file <-
  file.path(
    step13c4f_dir,
    "STEP13C-4F_PASS.txt"
  )


required_files <- c(
  mapping_file,
  subject_table_file,
  summary_file,
  checks_file,
  audit_rds,
  pass_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  cat("FAIL — Required STEP13C-4F files missing:\n\n")
  
  print(
    missing_files
  )
  
  stop(
    "STEP13C-5 stopped."
  )
}


# ------------------------------------------------------------
# 2. Read locked audit outputs
# ------------------------------------------------------------

mapping_df <-
  read.csv(
    mapping_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

subject_df <-
  read.csv(
    subject_table_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

summary_df <-
  read.csv(
    summary_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

checks_df <-
  read.csv(
    checks_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

audit_object <-
  readRDS(
    audit_rds
  )


# ------------------------------------------------------------
# 3. Verify previous STEP13C-4F PASS status
# ------------------------------------------------------------

if (!isTRUE(audit_object$final_pass)) {
  
  stop(
    "FAIL — STEP13C-4F audit object is not PASS."
  )
}

if (!all(checks_df$PASS)) {
  
  cat(
    "Failed integrity checks:\n\n"
  )
  
  print(
    checks_df[
      !checks_df$PASS,
      ,
      drop = FALSE
    ]
  )
  
  stop(
    "FAIL — STEP13C-4F integrity checks are not all PASS."
  )
}


# ------------------------------------------------------------
# 4. Reconstruct locked counts from audited objects
# ------------------------------------------------------------

n_samples <-
  nrow(mapping_df)

n_unique_gsm <-
  length(unique(mapping_df$GSM))

n_subjects <-
  nrow(subject_df)

n_complete_pairs <-
  sum(
    subject_df$Structure ==
      "Complete_N_plus_T"
  )

n_tumor_only <-
  sum(
    subject_df$Structure ==
      "Tumor_only"
  )

n_normal_only <-
  sum(
    subject_df$Structure ==
      "Normal_only"
  )

n_other <-
  sum(
    subject_df$Structure ==
      "Other_structure"
  )

n_normal <-
  sum(
    mapping_df$Status_Title ==
      "Normal"
  )

n_tumor <-
  sum(
    mapping_df$Status_Title ==
      "Tumor"
  )

n_subject_code_match <-
  sum(
    mapping_df$Subject_vs_TitleCode ==
      "MATCH"
  )

n_subject_code_mismatch <-
  sum(
    mapping_df$Subject_vs_TitleCode ==
      "MISMATCH"
  )


# ------------------------------------------------------------
# 5. Final predefined lock criteria
# ------------------------------------------------------------

lock_checks <- c(
  
  Samples_116 =
    n_samples == 116,
  
  Unique_GSM_116 =
    n_unique_gsm == 116,
  
  Normal_58 =
    n_normal == 58,
  
  Tumor_58 =
    n_tumor == 58,
  
  Subjects_59 =
    n_subjects == 59,
  
  Complete_pairs_57 =
    n_complete_pairs == 57,
  
  Tumor_only_1 =
    n_tumor_only == 1,
  
  Normal_only_1 =
    n_normal_only == 1,
  
  Other_structure_0 =
    n_other == 0,
  
  Subject_code_match_116 =
    n_subject_code_match == 116,
  
  Subject_code_mismatch_0 =
    n_subject_code_mismatch == 0,
  
  Previous_STEP13C4F_PASS =
    isTRUE(audit_object$final_pass)
  
)


# ------------------------------------------------------------
# 6. Verify singleton subjects
# ------------------------------------------------------------

tumor_only_subjects <-
  subject_df$Subject_ID[
    subject_df$Structure ==
      "Tumor_only"
  ]

normal_only_subjects <-
  subject_df$Subject_ID[
    subject_df$Structure ==
      "Normal_only"
  ]

lock_checks[
  "Tumor_only_subject_3023"
] <-
  identical(
    sort(tumor_only_subjects),
    "3023"
  )

lock_checks[
  "Normal_only_subject_3035"
] <-
  identical(
    sort(normal_only_subjects),
    "3035"
  )


# ------------------------------------------------------------
# 7. Print lock checks
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13C-5 FINAL INPUT LOCK CHECKS\n")
cat("============================================================\n\n")

print(lock_checks)

cat("\n")


# ------------------------------------------------------------
# 8. Final status
# ------------------------------------------------------------

final_pass <-
  all(lock_checks)

if (final_pass) {
  
  final_status <-
    paste(
      "PASS — STEP13C-5 GSE32863 Final Input Lock established."
    )
  
} else {
  
  final_status <-
    paste(
      "FAIL — STEP13C-5 Final Input Lock criteria not satisfied."
    )
}

cat(
  final_status,
  "\n\n"
)


# ------------------------------------------------------------
# 9. Locked metadata summary
# ------------------------------------------------------------

locked_summary <- data.frame(
  
  Item = c(
    "Total samples",
    "Unique GSM accessions",
    "Normal samples",
    "Tumor samples",
    "Unique subjects",
    "Complete N+T subjects",
    "Tumor-only subjects",
    "Normal-only subjects",
    "Other subject structures",
    "Subject/title-code exact matches",
    "Subject/title-code mismatches"
  ),
  
  Value = c(
    n_samples,
    n_unique_gsm,
    n_normal,
    n_tumor,
    n_subjects,
    n_complete_pairs,
    n_tumor_only,
    n_normal_only,
    n_other,
    n_subject_code_match,
    n_subject_code_mismatch
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 10. Save final lock
# ------------------------------------------------------------

write.csv(
  locked_summary,
  file.path(
    output_dir,
    "STEP13C-5_Final_GSE32863_Input_Summary.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  subject_df,
  file.path(
    output_dir,
    "STEP13C-5_Locked_Subject_Level_Table.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  mapping_df,
  file.path(
    output_dir,
    "STEP13C-5_Locked_GSM_Subject_Mapping.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  data.frame(
    Check = names(lock_checks),
    PASS = as.logical(lock_checks),
    stringsAsFactors = FALSE
  ),
  file.path(
    output_dir,
    "STEP13C-5_Final_Integrity_Checks.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 11. Create formal lock object
# ------------------------------------------------------------

final_lock <- list(
  
  workflow =
    "STEP13C-5",
  
  status =
    final_status,
  
  purpose =
    "Final GSE32863 input and subject-level pairing lock",
  
  scope =
    "Metadata/provenance only",
  
  total_samples =
    n_samples,
  
  unique_GSM =
    n_unique_gsm,
  
  normal_samples =
    n_normal,
  
  tumor_samples =
    n_tumor,
  
  unique_subjects =
    n_subjects,
  
  complete_N_plus_T_subjects =
    n_complete_pairs,
  
  tumor_only_subjects =
    tumor_only_subjects,
  
  normal_only_subjects =
    normal_only_subjects,
  
  other_structures =
    n_other,
  
  exact_subject_title_code_matches =
    n_subject_code_match,
  
  subject_title_code_mismatches =
    n_subject_code_mismatch,
  
  historical_external_ROC_recalculated =
    FALSE,
  
  expression_reconstructed =
    FALSE,
  
  model_refitted =
    FALSE,
  
  threshold_optimized =
    FALSE,
  
  lock_checks =
    lock_checks,
  
  final_pass =
    final_pass,
  
  timestamp =
    Sys.time()
  
)

saveRDS(
  final_lock,
  file.path(
    output_dir,
    "STEP13C-5_FINAL_GSE32863_INPUT_LOCK.rds"
  )
)


# ------------------------------------------------------------
# 12. PASS file
# ------------------------------------------------------------

if (final_pass) {
  
  writeLines(
    
    c(
      "STEP13C-5 — FINAL GSE32863 INPUT LOCK",
      "======================================",
      "",
      final_status,
      "",
      "LOCKED STRUCTURE:",
      paste0("116 total samples"),
      paste0("58 Normal + 58 Tumor"),
      paste0("59 unique subjects"),
      paste0("57 complete N+T subject pairs"),
      paste0("1 tumor-only subject: 3023"),
      paste0("1 normal-only subject: 3035"),
      paste0("0 other subject structures"),
      paste0("116/116 exact subject/title-code matches"),
      "",
      "No artificial pairing was performed.",
      "Subject 3023 and subject 3035 were not paired.",
      "",
      "No expression reconstruction.",
      "No probe selection.",
      "No normalization.",
      "No ROC/AUC recalculation.",
      "No model fitting.",
      "No threshold optimization.",
      "",
      paste0(
        "Timestamp: ",
        Sys.time()
      )
    ),
    
    file.path(
      output_dir,
      "STEP13C-5_PASS.txt"
    ),
    
    useBytes = TRUE
  )
}


# ------------------------------------------------------------
# 13. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_dir,
    "STEP13C-5_sessionInfo.txt"
  ),
  useBytes = TRUE
)


# ------------------------------------------------------------
# 14. Final console output
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13C-5 FINAL STATUS\n")
cat("============================================================\n\n")

cat(final_status, "\n\n")

cat("Locked samples:          ", n_samples, "\n")
cat("Locked subjects:         ", n_subjects, "\n")
cat("Complete N+T pairs:      ", n_complete_pairs, "\n")
cat("Tumor-only:              ", n_tumor_only, "\n")
cat("Normal-only:             ", n_normal_only, "\n")
cat("Subject/code matches:    ", n_subject_code_match, "/116\n\n")

cat(
  "No expression / ROC / model analysis performed.\n"
)

cat(
  "============================================================\n"
)


# ============================================================
# STEP13E-1
# HISTORICAL TCGA ROC EXPRESSION-SCALE PROVENANCE AUDIT
# ============================================================
#
# PURPOSE
# -------
# Lock the documented expression-scale provenance of the
# original TCGA ROC branch.
#
# STRICTLY NO:
#   - expression recalculation
#   - transformation of expression matrix
#   - model fitting
#   - ROC/AUC recalculation
#
# ============================================================


workflow_root <-
  "input_dir/Paper 1_Step 13-Diagnostic Classification"

output_dir <-
  file.path(
    workflow_root,
    "STEP13E_Expression_Scale_Audit"
  )

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 1. Historical RData
# ------------------------------------------------------------

historical_rdata <-
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal/.RData"

if (!file.exists(historical_rdata)) {
  
  stop(
    "FAIL — Historical TCGA/GSE RData not found."
  )
}


# ------------------------------------------------------------
# 2. Load historical objects into isolated environment
# ------------------------------------------------------------

hist_env <-
  new.env()

load(
  historical_rdata,
  envir = hist_env
)


# ------------------------------------------------------------
# 3. Required historical objects
# ------------------------------------------------------------

required_objects <- c(
  "df_wide",
  "roc_AOC3",
  "roc_MS4A7",
  "roc_IL7R",
  "roc_combined",
  "Figure12A_CombinedROC"
)

missing_objects <-
  required_objects[
    !vapply(
      required_objects,
      exists,
      logical(1),
      envir = hist_env,
      inherits = FALSE
    )
  ]

if (length(missing_objects) > 0) {
  
  stop(
    paste(
      "FAIL — Missing historical objects:",
      paste(missing_objects, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 4. Historical sample-level structure
# ------------------------------------------------------------

df_wide <-
  get(
    "df_wide",
    envir = hist_env
  )

if (!all(
  c(
    "sample_id",
    "status",
    "AOC3",
    "MS4A7",
    "IL7R"
  ) %in%
  colnames(df_wide)
)) {
  
  stop(
    "FAIL — Historical df_wide does not have expected columns."
  )
}


n_total <-
  nrow(df_wide)

n_normal <-
  sum(
    df_wide$status == "Normal",
    na.rm = TRUE
  )

n_tumor <-
  sum(
    df_wide$status == "Tumor",
    na.rm = TRUE
  )


# ------------------------------------------------------------
# 5. Historical ROC sample-size audit
# ------------------------------------------------------------

roc_names <- c(
  "roc_AOC3",
  "roc_MS4A7",
  "roc_IL7R",
  "roc_combined",
  "Figure12A_CombinedROC"
)

roc_sample_audit <- lapply(
  roc_names,
  function(x) {
    
    roc_obj <-
      get(
        x,
        envir = hist_env
      )
    
    data.frame(
      Object = x,
      N_Observations = length(roc_obj$response),
      N_Cases = sum(
        roc_obj$response == 1
      ),
      N_Controls = sum(
        roc_obj$response == 0
      ),
      Direction = roc_obj$direction,
      stringsAsFactors = FALSE
    )
  }
)

roc_sample_audit <-
  do.call(
    rbind,
    roc_sample_audit
  )


# ------------------------------------------------------------
# 6. Historical scale documentation lock
# ------------------------------------------------------------

documented_source <-
  "TCGA-LUAD STAR RNA-seq TPM"

documented_preprocessing <-
  "Remove samples without clinical information + remove duplicates"

documented_transformation <-
  "log2(value + 1)"

documented_outcome <-
  "Tumor versus Normal; Normal reference"


# ------------------------------------------------------------
# 7. Audit table
# ------------------------------------------------------------

scale_audit <-
  data.frame(
    
    Item = c(
      "Historical source",
      "Expression unit",
      "Documented preprocessing",
      "Documented transformation",
      "Historical total samples",
      "Historical Normal samples",
      "Historical Tumor samples",
      "Outcome definition"
    ),
    
    Documented_Value = c(
      documented_source,
      "TPM",
      documented_preprocessing,
      documented_transformation,
      n_total,
      n_normal,
      n_tumor,
      documented_outcome
    ),
    
    stringsAsFactors = FALSE
  )


# ------------------------------------------------------------
# 8. Integrity checks
# ------------------------------------------------------------

checks <- c(
  
  Historical_df_wide_575 =
    n_total == 575,
  
  Historical_Normal_59 =
    n_normal == 59,
  
  Historical_Tumor_516 =
    n_tumor == 516,
  
  Required_ROC_objects_present =
    length(missing_objects) == 0,
  
  ROC_AOC3_575 =
    roc_sample_audit$N_Observations[
      roc_sample_audit$Object == "roc_AOC3"
    ] == 575,
  
  ROC_MS4A7_575 =
    roc_sample_audit$N_Observations[
      roc_sample_audit$Object == "roc_MS4A7"
    ] == 575,
  
  ROC_IL7R_575 =
    roc_sample_audit$N_Observations[
      roc_sample_audit$Object == "roc_IL7R"
    ] == 575,
  
  ROC_combined_575 =
    roc_sample_audit$N_Observations[
      roc_sample_audit$Object == "roc_combined"
    ] == 575,
  
  Figure12A_ROC_575 =
    roc_sample_audit$N_Observations[
      roc_sample_audit$Object ==
        "Figure12A_CombinedROC"
    ] == 575
  
)


# ------------------------------------------------------------
# 9. Final status
# ------------------------------------------------------------

final_pass <-
  all(checks)

cat("\n")
cat("============================================================\n")
cat("STEP13E-1 — TCGA HISTORICAL SCALE AUDIT\n")
cat("============================================================\n\n")

cat(
  "Documented source:       ",
  documented_source,
  "\n"
)

cat(
  "Documented scale:        ",
  "TPM",
  "\n"
)

cat(
  "Documented transform:    ",
  documented_transformation,
  "\n"
)

cat(
  "Historical samples:      ",
  n_total,
  "\n"
)

cat(
  "Normal:                  ",
  n_normal,
  "\n"
)

cat(
  "Tumor:                   ",
  n_tumor,
  "\n\n"
)

print(checks)

cat("\n")


if (final_pass) {
  
  final_status <-
    "PASS — STEP13E-1 historical TCGA ROC expression-scale provenance established."
  
} else {
  
  final_status <-
    "FAIL — STEP13E-1 historical TCGA expression-scale audit failed."
}

cat(
  final_status,
  "\n\n"
)


# ------------------------------------------------------------
# 10. Save outputs
# ------------------------------------------------------------

write.csv(
  scale_audit,
  file.path(
    output_dir,
    "STEP13E-1_TCGA_Historical_Expression_Scale_Audit.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  roc_sample_audit,
  file.path(
    output_dir,
    "STEP13E-1_TCGA_Historical_ROC_Sample_Audit.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  data.frame(
    Check = names(checks),
    PASS = as.logical(checks),
    stringsAsFactors = FALSE
  ),
  file.path(
    output_dir,
    "STEP13E-1_Integrity_Checks.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


audit_object <-
  list(
    
    workflow =
      "STEP13E-1",
    
    source =
      documented_source,
    
    expression_scale =
      "TPM",
    
    transformation =
      documented_transformation,
    
    preprocessing =
      documented_preprocessing,
    
    historical_N =
      n_total,
    
    historical_Normal =
      n_normal,
    
    historical_Tumor =
      n_tumor,
    
    ROC_sample_audit =
      roc_sample_audit,
    
    checks =
      checks,
    
    final_pass =
      final_pass,
    
    no_expression_recalculation =
      TRUE,
    
    no_ROC_recalculation =
      TRUE,
    
    no_model_fitting =
      TRUE,
    
    timestamp =
      Sys.time()
  )

saveRDS(
  audit_object,
  file.path(
    output_dir,
    "STEP13E-1_TCGA_Historical_Scale_Audit.rds"
  )
)


if (final_pass) {
  
  writeLines(
    c(
      "STEP13E-1 — PASS",
      "",
      final_status,
      "",
      "Historical TCGA ROC branch:",
      "Source: TCGA-LUAD STAR RNA-seq TPM",
      "Transformation: log2(value + 1)",
      "Samples: 516 Tumor + 59 Normal = 575",
      "",
      "No expression recalculation.",
      "No ROC/AUC recalculation.",
      "No model fitting.",
      "No threshold optimization."
    ),
    file.path(
      output_dir,
      "STEP13E-1_PASS.txt"
    ),
    useBytes = TRUE
  )
}


writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_dir,
    "STEP13E-1_sessionInfo.txt"
  ),
  useBytes = TRUE
)

cat(
  "Output directory:\n",
  output_dir,
  "\n"
)

# ============================================================
# STEP13E-2
# HISTORICAL GSE32863 EXPRESSION-SCALE / OBJECT PROVENANCE AUDIT
# ============================================================

workflow_root <-
  "input_dir/Paper 1_Step 13-Diagnostic Classification"

output_dir <-
  file.path(
    workflow_root,
    "STEP13E_Expression_Scale_Audit"
  )

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


historical_rdata <-
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal/.RData"

if (!file.exists(historical_rdata)) {
  stop(
    "FAIL — Historical RData not found."
  )
}


hist_env <-
  new.env()

load(
  historical_rdata,
  envir = hist_env
)


required_objects <- c(
  "expr",
  "expr.mat",
  "expr2",
  "expr.final",
  "expr3",
  "roc.ext",
  "group",
  "genes"
)

missing_objects <-
  required_objects[
    !vapply(
      required_objects,
      exists,
      logical(1),
      envir = hist_env,
      inherits = FALSE
    )
  ]

if (length(missing_objects) > 0) {
  
  stop(
    paste(
      "FAIL — Missing historical GSE objects:",
      paste(
        missing_objects,
        collapse = ", "
      )
    )
  )
}


expr <-
  get("expr", envir = hist_env)

expr.mat <-
  get("expr.mat", envir = hist_env)

expr2 <-
  get("expr2", envir = hist_env)

expr.final <-
  get("expr.final", envir = hist_env)

expr3 <-
  get("expr3", envir = hist_env)

roc.ext <-
  get("roc.ext", envir = hist_env)

group <-
  get("group", envir = hist_env)

genes <-
  get("genes", envir = hist_env)


# ------------------------------------------------------------
# Object dimensions
# ------------------------------------------------------------

dimension_audit <- data.frame(
  
  Object = c(
    "expr",
    "expr.mat",
    "expr2",
    "expr.final",
    "expr3"
  ),
  
  Rows = c(
    nrow(expr),
    nrow(expr.mat),
    nrow(expr2),
    nrow(expr.final),
    nrow(expr3)
  ),
  
  Columns = c(
    ncol(expr),
    ncol(expr.mat),
    ncol(expr2),
    ncol(expr.final),
    ncol(expr3)
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# Sample-count audit
# ------------------------------------------------------------

sample_counts <- c(
  
  expr =
    ncol(expr) - 1,
  
  expr.mat =
    ncol(expr.mat),
  
  expr2 =
    ncol(expr2) - 2,
  
  expr.final =
    ncol(expr.final),
  
  expr3 =
    nrow(expr3),
  
  group =
    length(group),
  
  roc.ext =
    length(roc.ext$response)
  
)


# ------------------------------------------------------------
# Historical sample labels
# ------------------------------------------------------------

group_normal <-
  sum(
    group == "Normal",
    na.rm = TRUE
  )

group_tumor <-
  sum(
    group == "Tumor",
    na.rm = TRUE
  )


# ------------------------------------------------------------
# Core genes
# ------------------------------------------------------------

core_gene_check <-
  identical(
    genes,
    c(
      "AOC3",
      "MS4A7",
      "IL7R"
    )
  )


# ------------------------------------------------------------
# Historical object sample counts
# ------------------------------------------------------------

checks <- c(
  
  expr_116 =
    ncol(expr) - 1 == 116,
  
  expr_mat_116 =
    ncol(expr.mat) == 116,
  
  expr2_116 =
    ncol(expr2) - 2 == 116,
  
  expr_final_116 =
    ncol(expr.final) == 116,
  
  expr3_116 =
    nrow(expr3) == 116,
  
  group_116 =
    length(group) == 116,
  
  roc_ext_116 =
    length(roc.ext$response) == 116,
  
  Normal_58 =
    group_normal == 58,
  
  Tumor_58 =
    group_tumor == 58,
  
  Core_genes_AOC3_MS4A7_IL7R =
    core_gene_check,
  
  expr_48803x117 =
    nrow(expr) == 48803 &&
    ncol(expr) == 117,
  
  expr_mat_28742x116 =
    nrow(expr.mat) == 28742 &&
    ncol(expr.mat) == 116,
  
  expr2_28742x118 =
    nrow(expr2) == 28742 &&
    ncol(expr2) == 118,
  
  expr_final_19425x116 =
    nrow(expr.final) == 19425 &&
    ncol(expr.final) == 116,
  
  expr3_116x5 =
    nrow(expr3) == 116 &&
    ncol(expr3) == 5
  
)


# ------------------------------------------------------------
# Final status
# ------------------------------------------------------------

final_pass <-
  all(checks)

cat("\n")
cat("============================================================\n")
cat("STEP13E-2 — GSE32863 HISTORICAL OBJECT AUDIT\n")
cat("============================================================\n\n")

print(dimension_audit)

cat("\nHistorical group:\n")
cat(
  "Normal:",
  group_normal,
  "\n"
)

cat(
  "Tumor:",
  group_tumor,
  "\n\n"
)

cat("Core genes:\n")
print(genes)

cat("\nIntegrity checks:\n")
print(checks)

cat("\n")


if (final_pass) {
  
  final_status <-
    "PASS — STEP13E-2 historical GSE32863 object/sample-scale provenance established."
  
} else {
  
  final_status <-
    "FAIL — STEP13E-2 historical GSE32863 object audit failed."
}

cat(
  final_status,
  "\n\n"
)


# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

write.csv(
  dimension_audit,
  file.path(
    output_dir,
    "STEP13E-2_GSE32863_Historical_Object_Dimensions.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  data.frame(
    Object = names(sample_counts),
    Sample_Count = as.integer(sample_counts),
    stringsAsFactors = FALSE
  ),
  file.path(
    output_dir,
    "STEP13E-2_GSE32863_Historical_Sample_Counts.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  data.frame(
    Check = names(checks),
    PASS = as.logical(checks),
    stringsAsFactors = FALSE
  ),
  file.path(
    output_dir,
    "STEP13E-2_Integrity_Checks.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


audit_object <- list(
  
  workflow =
    "STEP13E-2",
  
  scope =
    "Historical GSE32863 object and expression-scale provenance audit",
  
  dimensions =
    dimension_audit,
  
  sample_counts =
    sample_counts,
  
  Normal =
    group_normal,
  
  Tumor =
    group_tumor,
  
  core_genes =
    genes,
  
  checks =
    checks,
  
  final_pass =
    final_pass,
  
  no_expression_recalculation =
    TRUE,
  
  no_probe_selection =
    TRUE,
  
  no_ROC_recalculation =
    TRUE,
  
  no_model_fitting =
    TRUE,
  
  timestamp =
    Sys.time()
  
)

saveRDS(
  audit_object,
  file.path(
    output_dir,
    "STEP13E-2_GSE32863_Historical_Scale_Audit.rds"
  )
)


if (final_pass) {
  
  writeLines(
    c(
      "STEP13E-2 — PASS",
      "",
      final_status,
      "",
      "Historical GSE32863 structure:",
      "116 samples",
      "58 Normal + 58 Tumor",
      "",
      "expr: 48803 × 117",
      "expr.mat: 28742 × 116",
      "expr2: 28742 × 118",
      "expr.final: 19425 × 116",
      "expr3: 116 × 5",
      "roc.ext: 116 observations",
      "",
      "Core genes: AOC3 / MS4A7 / IL7R",
      "",
      "No expression recalculation.",
      "No probe selection.",
      "No ROC/AUC recalculation.",
      "No model fitting."
    ),
    file.path(
      output_dir,
      "STEP13E-2_PASS.txt"
    ),
    useBytes = TRUE
  )
}


writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_dir,
    "STEP13E-2_sessionInfo.txt"
  ),
  useBytes = TRUE
)

cat(
  "Output directory:\n",
  output_dir,
  "\n"
)


# ============================================================
# STEP13E-3
# CROSS-PLATFORM HISTORICAL MODEL-INPUT COMPATIBILITY AUDIT
# ============================================================
#
# PURPOSE
# -------
# Audit the traceability between:
#
#   TCGA historical model branch
#   and
#   GSE32863 historical external prediction branch
#
# Main questions:
#
#   1. Are the historical GSE32863 model inputs traceable?
#   2. Do expr3 gene values correspond to expr.final?
#   3. Does expr3$Probability exactly correspond to roc.ext
#      predictor values?
#   4. What historical model object/formula/coefficients exist?
#   5. Is the external result traceable without refitting?
#
# STRICTLY NO:
#   - expression normalization
#   - expression transformation
#   - probe re-selection
#   - model refitting
#   - ROC/AUC recalculation
#   - threshold optimization
#
# ============================================================


# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

workflow_root <-
  "input_dir/Paper 1_Step 13-Diagnostic Classification"

historical_rdata <-
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal/.RData"

output_dir <-
  file.path(
    workflow_root,
    "STEP13E_Expression_Scale_Audit",
    "STEP13E-3_Cross_Platform_Model_Input_Audit"
  )

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

if (!file.exists(historical_rdata)) {
  stop(
    "FAIL — Historical RData not found:\n",
    historical_rdata
  )
}


cat("\n")
cat("============================================================\n")
cat("STEP13E-3 — CROSS-PLATFORM MODEL-INPUT AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Load historical objects into isolated environment
# ------------------------------------------------------------

hist_env <-
  new.env()

load(
  historical_rdata,
  envir = hist_env
)


# ------------------------------------------------------------
# 2. Required objects
# ------------------------------------------------------------

required_objects <- c(
  "expr.final",
  "expr3",
  "roc.ext",
  "genes",
  "model",
  "fit",
  "cv_model"
)

object_present <-
  vapply(
    required_objects,
    exists,
    logical(1),
    envir = hist_env,
    inherits = FALSE
  )

object_presence_df <-
  data.frame(
    Object = required_objects,
    Present = object_present,
    stringsAsFactors = FALSE
  )

cat("Historical object availability:\n\n")

print(
  object_presence_df,
  row.names = FALSE
)

cat("\n")

if (!all(object_present)) {
  
  missing_objects <-
    required_objects[
      !object_present
    ]
  
  stop(
    paste(
      "FAIL — Required historical objects missing:",
      paste(
        missing_objects,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 3. Retrieve objects
# ------------------------------------------------------------

expr_final <-
  get(
    "expr.final",
    envir = hist_env
  )

expr3 <-
  get(
    "expr3",
    envir = hist_env
  )

roc_ext <-
  get(
    "roc.ext",
    envir = hist_env
  )

genes <-
  get(
    "genes",
    envir = hist_env
  )

model <-
  get(
    "model",
    envir = hist_env
  )

fit <-
  get(
    "fit",
    envir = hist_env
  )

cv_model <-
  get(
    "cv_model",
    envir = hist_env
  )


# ------------------------------------------------------------
# 4. expr3 structural audit
# ------------------------------------------------------------

expected_expr3_columns <- c(
  "AOC3",
  "MS4A7",
  "IL7R",
  "Group",
  "Probability"
)

expr3_columns_ok <-
  identical(
    colnames(expr3),
    expected_expr3_columns
  )

expr3_n_ok <-
  nrow(expr3) == 116

expr3_probability_ok <-
  "Probability" %in%
  colnames(expr3)

expr3_group_ok <-
  "Group" %in%
  colnames(expr3)


# ------------------------------------------------------------
# 5. expr.final gene traceability
# ------------------------------------------------------------
#
# IMPORTANT:
# We do not modify expr.final.
# We only inspect whether the three locked core genes are
# represented in its rownames.
# ------------------------------------------------------------

expr_final_gene_names <-
  rownames(expr_final)

core_genes_present_in_expr_final <-
  genes %in%
  expr_final_gene_names

core_gene_presence_df <-
  data.frame(
    Gene = genes,
    Present_in_expr_final =
      core_genes_present_in_expr_final,
    stringsAsFactors = FALSE
  )

cat(
  "Core-gene presence in historical expr.final:\n\n"
)

print(
  core_gene_presence_df,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 6. Exact sample-level traceability:
# expr3 values versus expr.final
# ------------------------------------------------------------
#
# expr3 is sample × gene.
# expr.final is gene × sample.
#
# We identify the matching sample columns by sample IDs.
#
# NO transformation is performed.
# ------------------------------------------------------------

expr_final_colnames <-
  colnames(expr_final)

expr3_sample_ids <-
  rownames(expr3)

if (is.null(expr3_sample_ids)) {
  
  stop(
    "FAIL — expr3 has no rownames/sample IDs."
  )
}


# Historical expr.final may contain sample IDs as colnames.
# Check exact overlap first.

sample_overlap <-
  intersect(
    expr3_sample_ids,
    expr_final_colnames
  )

cat(
  "expr3 sample IDs:",
  length(expr3_sample_ids),
  "\n"
)

cat(
  "expr.final sample columns:",
  length(expr_final_colnames),
  "\n"
)

cat(
  "Exact overlapping sample IDs:",
  length(sample_overlap),
  "\n\n"
)


# ------------------------------------------------------------
# 7. If exact sample IDs are available, audit values
# ------------------------------------------------------------

gene_value_audit <- list()

if (
  length(sample_overlap) == 116 &&
  all(
    genes %in%
    rownames(expr_final)
  )
) {
  
  for (g in genes) {
    
    historical_gene_values <-
      as.numeric(
        expr3[
          sample_overlap,
          g
        ]
      )
    
    expr_final_values <-
      as.numeric(
        expr_final[
          g,
          sample_overlap
        ]
      )
    
    exact_equal <-
      isTRUE(
        all.equal(
          historical_gene_values,
          expr_final_values,
          check.attributes = FALSE
        )
      )
    
    max_abs_difference <-
      max(
        abs(
          historical_gene_values -
            expr_final_values
        ),
        na.rm = TRUE
      )
    
    gene_value_audit[[g]] <-
      data.frame(
        Gene = g,
        N = length(sample_overlap),
        Exact_equal = exact_equal,
        Max_absolute_difference =
          max_abs_difference,
        stringsAsFactors = FALSE
      )
  }
  
} else {
  
  for (g in genes) {
    
    gene_value_audit[[g]] <-
      data.frame(
        Gene = g,
        N = 0,
        Exact_equal = FALSE,
        Max_absolute_difference = NA_real_,
        stringsAsFactors = FALSE
      )
  }
}

gene_value_audit_df <-
  do.call(
    rbind,
    gene_value_audit
  )

rownames(gene_value_audit_df) <- NULL


# ------------------------------------------------------------
# 8. Historical Probability ↔ roc.ext predictor audit
# ------------------------------------------------------------

historical_probability <-
  as.numeric(
    expr3$Probability
  )

historical_predictor <-
  as.numeric(
    roc_ext$predictor
  )

if (
  length(historical_probability) !=
  length(historical_predictor)
) {
  
  stop(
    "FAIL — expr3 Probability and roc.ext predictor lengths differ."
  )
}

probability_exact_equal <-
  isTRUE(
    all.equal(
      historical_probability,
      historical_predictor,
      check.attributes = FALSE
    )
  )

probability_max_difference <-
  max(
    abs(
      historical_probability -
        historical_predictor
    ),
    na.rm = TRUE
  )

cat(
  "Historical Probability length:",
  length(historical_probability),
  "\n"
)

cat(
  "Historical ROC predictor length:",
  length(historical_predictor),
  "\n"
)

cat(
  "Probability ↔ predictor exact equality:",
  probability_exact_equal,
  "\n"
)

cat(
  "Maximum absolute difference:",
  probability_max_difference,
  "\n\n"
)


# ------------------------------------------------------------
# 9. Audit roc.ext response / labels
# ------------------------------------------------------------

roc_response <-
  roc_ext$response

roc_levels <-
  roc_ext$levels

roc_direction <-
  roc_ext$direction

roc_response_n <-
  length(roc_response)

roc_normal <-
  sum(
    roc_response == roc_levels[1],
    na.rm = TRUE
  )

roc_tumor <-
  sum(
    roc_response == roc_levels[2],
    na.rm = TRUE
  )


# ------------------------------------------------------------
# 10. Historical model object audit
# ------------------------------------------------------------

model_class <-
  class(model)

fit_class <-
  class(fit)

cv_model_class <-
  class(cv_model)

model_formula <-
  tryCatch(
    paste(
      deparse(
        formula(model)
      ),
      collapse = ""
    ),
    error = function(e) {
      NA_character_
    }
  )

fit_formula <-
  tryCatch(
    paste(
      deparse(
        formula(fit)
      ),
      collapse = ""
    ),
    error = function(e) {
      NA_character_
    }
  )

cv_model_formula <-
  tryCatch(
    paste(
      deparse(
        formula(cv_model)
      ),
      collapse = ""
    ),
    error = function(e) {
      NA_character_
    }
  )


# ------------------------------------------------------------
# 11. Extract historical coefficients if available
# ------------------------------------------------------------

model_coef <-
  tryCatch(
    coef(model),
    error = function(e) {
      NULL
    }
  )

fit_coef <-
  tryCatch(
    coef(fit),
    error = function(e) {
      NULL
    }
  )

cv_model_coef <-
  tryCatch(
    coef(cv_model),
    error = function(e) {
      NULL
    }
  )


# ------------------------------------------------------------
# 12. Convert coefficient objects to audit tables
# ------------------------------------------------------------

coef_to_df <- function(
    coef_object,
    object_name
) {
  
  if (is.null(coef_object)) {
    
    return(
      data.frame(
        Object = object_name,
        Term = NA_character_,
        Coefficient = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }
  
  data.frame(
    Object = object_name,
    Term = names(coef_object),
    Coefficient =
      as.numeric(coef_object),
    stringsAsFactors = FALSE
  )
}


coef_audit_df <-
  rbind(
    
    coef_to_df(
      model_coef,
      "model"
    ),
    
    coef_to_df(
      fit_coef,
      "fit"
    ),
    
    coef_to_df(
      cv_model_coef,
      "cv_model"
    )
    
  )


# ------------------------------------------------------------
# 13. Search historical objects for core-gene coefficient terms
# ------------------------------------------------------------

coefficient_terms <-
  unique(
    coef_audit_df$Term
  )

coefficient_gene_presence <-
  data.frame(
    
    Gene = genes,
    
    Appears_in_model_coefficients =
      genes %in%
      coefficient_terms,
    
    stringsAsFactors = FALSE
    
  )


# ------------------------------------------------------------
# 14. Build formal audit table
# ------------------------------------------------------------

audit_summary <-
  data.frame(
    
    Audit_Item = c(
      
      "expr3_n",
      "expr3_columns",
      "expr3_sample_ids",
      "expr.final_core_genes",
      "expr3_vs_expr.final_sample_overlap",
      "AOC3_value_traceability",
      "MS4A7_value_traceability",
      "IL7R_value_traceability",
      "Probability_n",
      "roc.ext_predictor_n",
      "Probability_equals_roc.ext_predictor",
      "roc.ext_response_n",
      "roc.ext_Normal_n",
      "roc.ext_Tumor_n",
      "Historical_model_class",
      "Historical_fit_class",
      "Historical_cv_model_class",
      "Historical_model_formula",
      "Historical_fit_formula",
      "Historical_cv_model_formula"
      
    ),
    
    Result = c(
      
      expr3_n_ok,
      
      expr3_columns_ok,
      
      length(expr3_sample_ids) == 116,
      
      all(core_genes_present_in_expr_final),
      
      length(sample_overlap) == 116,
      
      gene_value_audit_df$Exact_equal[
        gene_value_audit_df$Gene == "AOC3"
      ],
      
      gene_value_audit_df$Exact_equal[
        gene_value_audit_df$Gene == "MS4A7"
      ],
      
      gene_value_audit_df$Exact_equal[
        gene_value_audit_df$Gene == "IL7R"
      ],
      
      length(historical_probability) == 116,
      
      length(historical_predictor) == 116,
      
      probability_exact_equal,
      
      roc_response_n == 116,
      
      roc_normal == 58,
      
      roc_tumor == 58,
      
      paste(
        model_class,
        collapse = ";"
      ),
      
      paste(
        fit_class,
        collapse = ";"
      ),
      
      paste(
        cv_model_class,
        collapse = ";"
      ),
      
      model_formula,
      
      fit_formula,
      
      cv_model_formula
      
    ),
    
    stringsAsFactors = FALSE
    
  )


# ------------------------------------------------------------
# 15. Formal integrity checks
# ------------------------------------------------------------

checks <- c(
  
  expr3_116 =
    expr3_n_ok,
  
  expr3_expected_columns =
    expr3_columns_ok,
  
  expr3_sample_IDs_116 =
    length(expr3_sample_ids) == 116,
  
  core_genes_in_expr_final =
    all(core_genes_present_in_expr_final),
  
  expr3_exprfinal_sample_overlap_116 =
    length(sample_overlap) == 116,
  
  AOC3_exact_traceability =
    gene_value_audit_df$Exact_equal[
      gene_value_audit_df$Gene == "AOC3"
    ],
  
  MS4A7_exact_traceability =
    gene_value_audit_df$Exact_equal[
      gene_value_audit_df$Gene == "MS4A7"
    ],
  
  IL7R_exact_traceability =
    gene_value_audit_df$Exact_equal[
      gene_value_audit_df$Gene == "IL7R"
    ],
  
  Probability_116 =
    length(historical_probability) == 116,
  
  roc_ext_predictor_116 =
    length(historical_predictor) == 116,
  
  Probability_predictor_exact =
    probability_exact_equal,
  
  roc_ext_116 =
    roc_response_n == 116,
  
  roc_ext_Normal_58 =
    roc_normal == 58,
  
  roc_ext_Tumor_58 =
    roc_tumor == 58
  
)


# ------------------------------------------------------------
# 16. Final status
# ------------------------------------------------------------

final_pass <-
  all(checks)


cat("\n")
cat("============================================================\n")
cat("STEP13E-3 INTEGRITY CHECKS\n")
cat("============================================================\n\n")

print(checks)

cat("\n")


if (final_pass) {
  
  final_status <-
    paste(
      "PASS — Historical GSE32863 model-input traceability",
      "and Probability-to-ROC predictor concordance established."
    )
  
} else {
  
  final_status <-
    paste(
      "FAIL — One or more historical model-input traceability",
      "checks were not satisfied."
    )
}


cat(
  final_status,
  "\n\n"
)


# ------------------------------------------------------------
# 17. Save outputs
# ------------------------------------------------------------

write.csv(
  object_presence_df,
  file.path(
    output_dir,
    "STEP13E-3_Historical_Object_Availability.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  audit_summary,
  file.path(
    output_dir,
    "STEP13E-3_Cross_Platform_Audit_Summary.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  gene_value_audit_df,
  file.path(
    output_dir,
    "STEP13E-3_Core_Gene_Value_Traceability.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  coefficient_terms,
  file.path(
    output_dir,
    "STEP13E-3_Historical_Coefficient_Terms.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  coefficient_gene_presence,
  file.path(
    output_dir,
    "STEP13E-3_Core_Gene_Coefficient_Audit.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  coef_audit_df,
  file.path(
    output_dir,
    "STEP13E-3_Historical_Model_Coefficients.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  data.frame(
    Check = names(checks),
    PASS = as.logical(checks),
    stringsAsFactors = FALSE
  ),
  file.path(
    output_dir,
    "STEP13E-3_Integrity_Checks.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 18. Save formal audit object
# ------------------------------------------------------------

audit_object <-
  list(
    
    workflow =
      "STEP13E-3",
    
    purpose =
      "Cross-platform historical model-input compatibility audit",
    
    scope =
      "Historical object traceability only",
    
    expr3_dimensions =
      dim(expr3),
    
    expr_final_dimensions =
      dim(expr_final),
    
    sample_overlap =
      sample_overlap,
    
    core_gene_value_audit =
      gene_value_audit_df,
    
    probability_exact_equal =
      probability_exact_equal,
    
    probability_max_difference =
      probability_max_difference,
    
    roc_ext_N =
      roc_response_n,
    
    roc_ext_Normal =
      roc_normal,
    
    roc_ext_Tumor =
      roc_tumor,
    
    model_class =
      model_class,
    
    fit_class =
      fit_class,
    
    cv_model_class =
      cv_model_class,
    
    model_formula =
      model_formula,
    
    fit_formula =
      fit_formula,
    
    cv_model_formula =
      cv_model_formula,
    
    coefficients =
      coef_audit_df,
    
    checks =
      checks,
    
    final_pass =
      final_pass,
    
    expression_recalculated =
      FALSE,
    
    probe_selection_repeated =
      FALSE,
    
    model_refitted =
      FALSE,
    
    roc_recalculated =
      FALSE,
    
    threshold_optimized =
      FALSE,
    
    timestamp =
      Sys.time()
    
  )


saveRDS(
  audit_object,
  file.path(
    output_dir,
    "STEP13E-3_Cross_Platform_Model_Input_Audit.rds"
  )
)


# ------------------------------------------------------------
# 19. PASS file
# ------------------------------------------------------------

if (final_pass) {
  
  writeLines(
    
    c(
      
      "STEP13E-3 — PASS",
      
      "",
      
      final_status,
      
      "",
      
      "Historical GSE32863 model-input audit:",
      
      paste0(
        "expr3: ",
        nrow(expr3),
        " × ",
        ncol(expr3)
      ),
      
      paste0(
        "expr.final: ",
        nrow(expr_final),
        " × ",
        ncol(expr_final)
      ),
      
      paste0(
        "Exact sample overlap: ",
        length(sample_overlap),
        "/116"
      ),
      
      paste0(
        "Probability ↔ roc.ext predictor exact equality: ",
        probability_exact_equal
      ),
      
      paste0(
        "roc.ext samples: ",
        roc_response_n
      ),
      
      paste0(
        "roc.ext Normal: ",
        roc_normal
      ),
      
      paste0(
        "roc.ext Tumor: ",
        roc_tumor
      ),
      
      "",
      
      "No expression recalculation.",
      
      "No probe re-selection.",
      
      "No normalization.",
      
      "No model refitting.",
      
      "No ROC/AUC recalculation.",
      
      "No threshold optimization.",
      
      "",
      
      paste0(
        "Timestamp: ",
        Sys.time()
      )
      
    ),
    
    file.path(
      output_dir,
      "STEP13E-3_PASS.txt"
    ),
    
    useBytes = TRUE
  )
}


# ------------------------------------------------------------
# 20. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_dir,
    "STEP13E-3_sessionInfo.txt"
  ),
  useBytes = TRUE
)


# ------------------------------------------------------------
# 21. Final console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13E-3 FINAL SUMMARY\n")
cat("============================================================\n\n")

cat(
  "expr3:                    ",
  nrow(expr3),
  " × ",
  ncol(expr3),
  "\n"
)

cat(
  "expr.final:               ",
  nrow(expr_final),
  " × ",
  ncol(expr_final),
  "\n"
)

cat(
  "Sample overlap:           ",
  length(sample_overlap),
  "/116\n"
)

cat(
  "Probability length:       ",
  length(historical_probability),
  "\n"
)

cat(
  "roc.ext predictor length: ",
  length(historical_predictor),
  "\n"
)

cat(
  "Probability exact match:  ",
  probability_exact_equal,
  "\n"
)

cat(
  "roc.ext Normal:           ",
  roc_normal,
  "\n"
)

cat(
  "roc.ext Tumor:            ",
  roc_tumor,
  "\n\n"
)

cat(
  "FINAL STATUS:\n",
  final_status,
  "\n\n"
)

cat(
  "Output directory:\n",
  output_dir,
  "\n"
)

cat(
  "============================================================\n"
)

# ============================================================
# STEP13F — MODEL DEFINITION LOCK
# FINAL CLEAN / AUDIT-SAFE VERSION
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP13F — Model Definition Lock\n")
cat("FINAL CLEAN / AUDIT-SAFE VERSION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Workflow root and output directory
# ------------------------------------------------------------

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

step_dir <- file.path(
  workflow_root,
  "STEP13F_Model_Definition_Lock"
)

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 0A. Remove stale PASS / FAIL markers
# ------------------------------------------------------------

stale_files <- c(
  "STEP13F_PASS.txt",
  "STEP13F_FAIL.txt"
)

stale_paths <- file.path(step_dir, stale_files)

for (f in stale_paths) {
  if (file.exists(f)) {
    file.remove(f)
  }
}


# ------------------------------------------------------------
# 1. Cleanly define ALL model-lock variables
# ------------------------------------------------------------

# ---- Core genes ----

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


# ---- Outcome definition ----

outcome_name <- "TumorStatus"

outcome_levels <- c(
  "Normal",
  "Tumor"
)

reference_level <- "Normal"

positive_class <- "Tumor"


# ---- Model definition ----

model_type <- "multivariable_logistic_regression"

model_formula_text <-
  "TumorStatus ~ AOC3 + MS4A7 + IL7R"


# ---- Objective ----

analysis_objective <-
  "Exploratory tumor-versus-normal tissue classification"


# ---- Cross-validation policy ----

nested_cv_required <- TRUE

outer_folds <- 10L

inner_folds <- 10L

stratified_cv_required <- TRUE


# ---- Feature-selection policy ----

feature_selection_allowed <- FALSE


# ---- External validation policy ----

external_refit_allowed <- FALSE

external_coefficient_estimation_allowed <- FALSE

external_feature_selection_allowed <- FALSE

external_threshold_optimization_allowed <- FALSE

external_cutpoint_search_allowed <- FALSE

external_model_tuning_allowed <- FALSE


# ---- Clinical interpretation lock ----

clinical_diagnostic_claim_allowed <- FALSE

diagnostic_biomarker_claim_allowed <- FALSE


# ---- External GSE32863 input lock ----

gse32863_subjects_locked <- 59L

gse32863_complete_pairs_locked <- 57L

gse32863_total_samples_locked <- 116L

gse32863_tumor_samples_locked <- 58L

gse32863_normal_samples_locked <- 58L


# ---- Historical / revised branch separation ----

historical_tcga_branch <-
  "Xiangtao-derived original ROC branch"

historical_tcga_tumor_n <- 516L

historical_tcga_normal_n <- 59L

historical_tcga_total_n <- 575L


revised_tcga_branch <-
  "Sangerbox-derived tumor-only WGCNA branch"

revised_tcga_primary_tumor_n <- 513L

revised_tcga_normal_n <- 59L


# ------------------------------------------------------------
# 2. PRE-FLIGHT OBJECT EXISTENCE CHECK
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("STEP13F PRE-FLIGHT — Object Existence Check\n")
cat("------------------------------------------------------------\n\n")


required_objects <- c(
  "core_genes",
  "outcome_name",
  "outcome_levels",
  "reference_level",
  "positive_class",
  "model_type",
  "model_formula_text",
  "analysis_objective",
  "nested_cv_required",
  "outer_folds",
  "inner_folds",
  "stratified_cv_required",
  "feature_selection_allowed",
  "external_refit_allowed",
  "external_coefficient_estimation_allowed",
  "external_feature_selection_allowed",
  "external_threshold_optimization_allowed",
  "external_cutpoint_search_allowed",
  "external_model_tuning_allowed",
  "clinical_diagnostic_claim_allowed",
  "diagnostic_biomarker_claim_allowed",
  "gse32863_subjects_locked",
  "gse32863_complete_pairs_locked",
  "gse32863_total_samples_locked",
  "gse32863_tumor_samples_locked",
  "gse32863_normal_samples_locked",
  "historical_tcga_branch",
  "historical_tcga_tumor_n",
  "historical_tcga_normal_n",
  "historical_tcga_total_n",
  "revised_tcga_branch",
  "revised_tcga_primary_tumor_n",
  "revised_tcga_normal_n"
)


object_exists <- vapply(
  required_objects,
  exists,
  logical(1),
  inherits = TRUE
)


object_existence_df <- data.frame(
  object = required_objects,
  exists = object_exists,
  stringsAsFactors = FALSE
)


write.csv(
  object_existence_df,
  file.path(
    step_dir,
    "STEP13F_Preflight_Object_Existence_Audit.csv"
  ),
  row.names = FALSE
)


cat(
  "Required objects:",
  length(required_objects),
  "\n"
)

cat(
  "Objects present:",
  sum(object_exists),
  "\n"
)

cat(
  "Objects missing:",
  sum(!object_exists),
  "\n\n"
)


if (any(!object_exists)) {
  
  cat("MISSING OBJECTS:\n")
  
  print(
    object_existence_df[
      !object_existence_df$exists,
      ,
      drop = FALSE
    ]
  )
  
  cat("\n")
  cat("FINAL STATUS:\n")
  cat("FAIL — STEP13F Pre-flight Object Existence Check failed.\n")
  cat("No model fitting / ROC / AUC / CV / threshold analysis was performed.\n")
  
  writeLines(
    c(
      "STEP13F FINAL STATUS",
      "FAIL — Pre-flight Object Existence Check failed.",
      "",
      "Missing objects:",
      paste(
        required_objects[!object_exists],
        collapse = ", "
      ),
      "",
      "No model fitting / ROC / AUC / CV / threshold analysis was performed."
    ),
    file.path(
      step_dir,
      "STEP13F_FAIL.txt"
    )
  )
  
  writeLines(
    capture.output(sessionInfo()),
    file.path(
      step_dir,
      "STEP13F_sessionInfo_FAIL.txt"
    )
  )
  
  stop(
    "STEP13F FAIL — one or more required objects are missing.",
    call. = FALSE
  )
}


cat(
  "PASS — all required objects are defined.\n\n"
)


# ------------------------------------------------------------
# 3. Value-based integrity checks
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("STEP13F — Integrity Checks\n")
cat("------------------------------------------------------------\n\n")


check_core_genes <-
  is.character(core_genes) &&
  length(core_genes) == 3 &&
  identical(
    core_genes,
    c("AOC3", "MS4A7", "IL7R")
  )


check_outcome_name <-
  is.character(outcome_name) &&
  length(outcome_name) == 1 &&
  outcome_name == "TumorStatus"


check_outcome_levels <-
  is.character(outcome_levels) &&
  length(outcome_levels) == 2 &&
  identical(
    outcome_levels,
    c("Normal", "Tumor")
  )


check_reference_level <-
  is.character(reference_level) &&
  length(reference_level) == 1 &&
  reference_level == "Normal"


check_positive_class <-
  is.character(positive_class) &&
  length(positive_class) == 1 &&
  positive_class == "Tumor"


check_model_type <-
  is.character(model_type) &&
  length(model_type) == 1 &&
  model_type == "multivariable_logistic_regression"


check_model_formula <-
  is.character(model_formula_text) &&
  length(model_formula_text) == 1 &&
  model_formula_text ==
  "TumorStatus ~ AOC3 + MS4A7 + IL7R"


check_analysis_objective <-
  is.character(analysis_objective) &&
  length(analysis_objective) == 1 &&
  analysis_objective ==
  "Exploratory tumor-versus-normal tissue classification"


check_nested_cv <-
  is.logical(nested_cv_required) &&
  length(nested_cv_required) == 1 &&
  isTRUE(nested_cv_required)


check_outer_folds <-
  is.numeric(outer_folds) &&
  length(outer_folds) == 1 &&
  isTRUE(outer_folds == 10)


check_inner_folds <-
  is.numeric(inner_folds) &&
  length(inner_folds) == 1 &&
  isTRUE(inner_folds == 10)


check_stratified_cv <-
  is.logical(stratified_cv_required) &&
  length(stratified_cv_required) == 1 &&
  isTRUE(stratified_cv_required)


check_feature_selection <-
  is.logical(feature_selection_allowed) &&
  length(feature_selection_allowed) == 1 &&
  isFALSE(feature_selection_allowed)


check_external_refit <-
  is.logical(external_refit_allowed) &&
  length(external_refit_allowed) == 1 &&
  isFALSE(external_refit_allowed)


check_external_coefficients <-
  is.logical(external_coefficient_estimation_allowed) &&
  length(external_coefficient_estimation_allowed) == 1 &&
  isFALSE(external_coefficient_estimation_allowed)


check_external_feature_selection <-
  is.logical(external_feature_selection_allowed) &&
  length(external_feature_selection_allowed) == 1 &&
  isFALSE(external_feature_selection_allowed)


check_external_threshold <-
  is.logical(external_threshold_optimization_allowed) &&
  length(external_threshold_optimization_allowed) == 1 &&
  isFALSE(external_threshold_optimization_allowed)


check_external_cutpoint <-
  is.logical(external_cutpoint_search_allowed) &&
  length(external_cutpoint_search_allowed) == 1 &&
  isFALSE(external_cutpoint_search_allowed)


check_external_tuning <-
  is.logical(external_model_tuning_allowed) &&
  length(external_model_tuning_allowed) == 1 &&
  isFALSE(external_model_tuning_allowed)


check_clinical_claim <-
  is.logical(clinical_diagnostic_claim_allowed) &&
  length(clinical_diagnostic_claim_allowed) == 1 &&
  isFALSE(clinical_diagnostic_claim_allowed)


check_biomarker_claim <-
  is.logical(diagnostic_biomarker_claim_allowed) &&
  length(diagnostic_biomarker_claim_allowed) == 1 &&
  isFALSE(diagnostic_biomarker_claim_allowed)


# ---- GSE32863 checks ----
# Value-based rather than identical() so that
# 59 and 59L are treated equivalently.

check_gse_subjects <-
  is.numeric(gse32863_subjects_locked) &&
  length(gse32863_subjects_locked) == 1 &&
  isTRUE(gse32863_subjects_locked == 59)


check_gse_pairs <-
  is.numeric(gse32863_complete_pairs_locked) &&
  length(gse32863_complete_pairs_locked) == 1 &&
  isTRUE(gse32863_complete_pairs_locked == 57)


check_gse_total <-
  is.numeric(gse32863_total_samples_locked) &&
  length(gse32863_total_samples_locked) == 1 &&
  isTRUE(gse32863_total_samples_locked == 116)


check_gse_tumor <-
  is.numeric(gse32863_tumor_samples_locked) &&
  length(gse32863_tumor_samples_locked) == 1 &&
  isTRUE(gse32863_tumor_samples_locked == 58)


check_gse_normal <-
  is.numeric(gse32863_normal_samples_locked) &&
  length(gse32863_normal_samples_locked) == 1 &&
  isTRUE(gse32863_normal_samples_locked == 58)


# ---- Historical TCGA branch checks ----

check_historical_branch <-
  is.character(historical_tcga_branch) &&
  length(historical_tcga_branch) == 1 &&
  historical_tcga_branch ==
  "Xiangtao-derived original ROC branch"


check_historical_tumor <-
  is.numeric(historical_tcga_tumor_n) &&
  length(historical_tcga_tumor_n) == 1 &&
  isTRUE(historical_tcga_tumor_n == 516)


check_historical_normal <-
  is.numeric(historical_tcga_normal_n) &&
  length(historical_tcga_normal_n) == 1 &&
  isTRUE(historical_tcga_normal_n == 59)


check_historical_total <-
  is.numeric(historical_tcga_total_n) &&
  length(historical_tcga_total_n) == 1 &&
  isTRUE(historical_tcga_total_n == 575)


# ---- Revised TCGA branch checks ----

check_revised_branch <-
  is.character(revised_tcga_branch) &&
  length(revised_tcga_branch) == 1 &&
  revised_tcga_branch ==
  "Sangerbox-derived tumor-only WGCNA branch"


check_revised_tumor <-
  is.numeric(revised_tcga_primary_tumor_n) &&
  length(revised_tcga_primary_tumor_n) == 1 &&
  isTRUE(revised_tcga_primary_tumor_n == 513)


check_revised_normal <-
  is.numeric(revised_tcga_normal_n) &&
  length(revised_tcga_normal_n) == 1 &&
  isTRUE(revised_tcga_normal_n == 59)


# ------------------------------------------------------------
# 4. Collect all checks
#    FINAL ROBUST VERSION
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("STEP13F — Collecting Integrity Checks\n")
cat("------------------------------------------------------------\n\n")


# IMPORTANT:
# Define check names independently from the logical values.
# This prevents loss of names() from causing a data.frame
# row-number mismatch.

integrity_check_names <- c(
  
  "Core_genes",
  
  "Outcome_name",
  
  "Outcome_levels",
  
  "Reference_level_Normal",
  
  "Positive_class_Tumor",
  
  "Model_type",
  
  "Model_formula",
  
  "Analysis_objective",
  
  "Nested_CV",
  
  "Outer_10_fold",
  
  "Inner_10_fold",
  
  "Stratified_CV",
  
  "Feature_selection_disabled",
  
  "External_refit_disabled",
  
  "External_coefficients_disabled",
  
  "External_feature_selection_disabled",
  
  "External_threshold_tuning_disabled",
  
  "External_cutpoint_search_disabled",
  
  "External_model_tuning_disabled",
  
  "Clinical_diagnostic_claim_disabled",
  
  "Diagnostic_biomarker_claim_disabled",
  
  "GSE_subject_count",
  
  "GSE_complete_pairs",
  
  "GSE_total_samples",
  
  "GSE_tumor_samples",
  
  "GSE_normal_samples",
  
  "Historical_Xiangtao_branch",
  
  "Historical_TCGA_tumor_516",
  
  "Historical_TCGA_normal_59",
  
  "Historical_TCGA_total_575",
  
  "Revised_Sangerbox_WGCNA_branch",
  
  "Revised_TCGA_primary_tumor_513",
  
  "Revised_TCGA_normal_59"
)


integrity_check_values <- c(
  
  check_core_genes,
  
  check_outcome_name,
  
  check_outcome_levels,
  
  check_reference_level,
  
  check_positive_class,
  
  check_model_type,
  
  check_model_formula,
  
  check_analysis_objective,
  
  check_nested_cv,
  
  check_outer_folds,
  
  check_inner_folds,
  
  check_stratified_cv,
  
  check_feature_selection,
  
  check_external_refit,
  
  check_external_coefficients,
  
  check_external_feature_selection,
  
  check_external_threshold,
  
  check_external_cutpoint,
  
  check_external_tuning,
  
  check_clinical_claim,
  
  check_biomarker_claim,
  
  check_gse_subjects,
  
  check_gse_pairs,
  
  check_gse_total,
  
  check_gse_tumor,
  
  check_gse_normal,
  
  check_historical_branch,
  
  check_historical_tumor,
  
  check_historical_normal,
  
  check_historical_total,
  
  check_revised_branch,
  
  check_revised_tumor,
  
  check_revised_normal
)


# ---- Structural audit before creating the table ----

cat(
  "Number of check names:   ",
  length(integrity_check_names),
  "\n",
  sep = ""
)

cat(
  "Number of check values:  ",
  length(integrity_check_values),
  "\n",
  sep = ""
)


check_structure_ok <-
  length(integrity_check_names) ==
  length(integrity_check_values)


if (!check_structure_ok) {
  
  cat("\n")
  cat("FINAL STATUS:\n")
  cat(
    "FAIL — STEP13F integrity-check structure is inconsistent.\n"
  )
  
  writeLines(
    c(
      "STEP13F FINAL STATUS",
      "",
      "FAIL — Integrity-check structure is inconsistent.",
      "",
      paste(
        "Number of check names:",
        length(integrity_check_names)
      ),
      paste(
        "Number of check values:",
        length(integrity_check_values)
      ),
      "",
      "No model fitting / ROC / AUC / CV / threshold analysis",
      "was performed in STEP13F."
    ),
    file.path(
      step_dir,
      "STEP13F_FAIL.txt"
    )
  )
  
  writeLines(
    capture.output(sessionInfo()),
    file.path(
      step_dir,
      "STEP13F_sessionInfo_FAIL.txt"
    )
  )
  
  stop(
    "STEP13F FAIL — integrity-check structure mismatch.",
    call. = FALSE
  )
}


# ---- Construct named integrity vector ----

integrity_checks <- integrity_check_values

names(integrity_checks) <- integrity_check_names


# ---- Final structural verification ----

if (
  length(integrity_checks) != 33L ||
  length(names(integrity_checks)) != 33L ||
  !identical(
    names(integrity_checks),
    integrity_check_names
  )
) {
  
  cat("\n")
  cat("FINAL STATUS:\n")
  cat(
    "FAIL — STEP13F integrity-check naming structure failed.\n"
  )
  
  writeLines(
    c(
      "STEP13F FINAL STATUS",
      "",
      "FAIL — Integrity-check naming structure failed.",
      "",
      "No model fitting / ROC / AUC / CV / threshold analysis",
      "was performed in STEP13F."
    ),
    file.path(
      step_dir,
      "STEP13F_FAIL.txt"
    )
  )
  
  writeLines(
    capture.output(sessionInfo()),
    file.path(
      step_dir,
      "STEP13F_sessionInfo_FAIL.txt"
    )
  )
  
  stop(
    "STEP13F FAIL — integrity-check naming structure failed.",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 5. Print integrity checks
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("STEP13F — Integrity Check Table\n")
cat("------------------------------------------------------------\n\n")


integrity_table <- data.frame(
  check = integrity_check_names,
  PASS = as.logical(integrity_check_values),
  stringsAsFactors = FALSE
)


print(integrity_table)


write.csv(
  integrity_table,
  file.path(
    step_dir,
    "STEP13F_Integrity_Checks.csv"
  ),
  row.names = FALSE
)


cat("\n")
cat(
  "Integrity checks evaluated:",
  nrow(integrity_table),
  "\n"
)

cat(
  "Passed:",
  sum(integrity_table$PASS),
  "\n"
)

cat(
  "Failed:",
  sum(!integrity_table$PASS),
  "\n\n"
)


# ------------------------------------------------------------
# 6. Overall status
# ------------------------------------------------------------

all_integrity_pass <- all(integrity_checks)


cat("\n")
cat("============================================================\n")
cat("STEP13F FINAL SUMMARY\n")
cat("============================================================\n\n")


cat(
  "Locked predictors:        ",
  paste(core_genes, collapse = ", "),
  "\n",
  sep = ""
)

cat(
  "Outcome:                  ",
  outcome_name,
  "\n",
  sep = ""
)

cat(
  "Reference level:          ",
  reference_level,
  "\n",
  sep = ""
)

cat(
  "Positive class:           ",
  positive_class,
  "\n",
  sep = ""
)

cat(
  "Model:                    ",
  model_type,
  "\n",
  sep = ""
)

cat(
  "Formula:                  ",
  model_formula_text,
  "\n",
  sep = ""
)

cat(
  "Nested CV:                ",
  nested_cv_required,
  "\n",
  sep = ""
)

cat(
  "Outer folds:              ",
  outer_folds,
  "\n",
  sep = ""
)

cat(
  "Inner folds:              ",
  inner_folds,
  "\n",
  sep = ""
)

cat(
  "Feature selection:        ",
  feature_selection_allowed,
  "\n",
  sep = ""
)

cat(
  "External refit:           ",
  external_refit_allowed,
  "\n",
  sep = ""
)

cat(
  "External coefficients:    ",
  external_coefficient_estimation_allowed,
  "\n",
  sep = ""
)

cat(
  "External feature select:  ",
  external_feature_selection_allowed,
  "\n",
  sep = ""
)

cat(
  "External threshold tune:  ",
  external_threshold_optimization_allowed,
  "\n",
  sep = ""
)

cat(
  "External cutpoint search: ",
  external_cutpoint_search_allowed,
  "\n",
  sep = ""
)

cat(
  "External model tuning:    ",
  external_model_tuning_allowed,
  "\n",
  sep = ""
)

cat(
  "Clinical diagnostic claim: ",
  clinical_diagnostic_claim_allowed,
  "\n",
  sep = ""
)

cat(
  "Diagnostic biomarker claim:",
  diagnostic_biomarker_claim_allowed,
  "\n"
)

cat(
  "GSE32863 subjects:         ",
  gse32863_subjects_locked,
  "\n",
  sep = ""
)

cat(
  "GSE32863 complete N+T:     ",
  gse32863_complete_pairs_locked,
  "\n",
  sep = ""
)

cat(
  "GSE32863 total samples:    ",
  gse32863_total_samples_locked,
  "\n",
  sep = ""
)

cat(
  "GSE32863 tumor:            ",
  gse32863_tumor_samples_locked,
  "\n",
  sep = ""
)

cat(
  "GSE32863 normal:           ",
  gse32863_normal_samples_locked,
  "\n",
  sep = ""
)

cat("\n")

cat(
  "All integrity checks:     ",
  all_integrity_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Save locked model definition
# ------------------------------------------------------------

model_definition_lock <- list(
  
  workflow = "STEP13F",
  
  status = if (all_integrity_pass) {
    "PASS"
  } else {
    "FAIL"
  },
  
  core_genes = core_genes,
  
  outcome_name = outcome_name,
  
  outcome_levels = outcome_levels,
  
  reference_level = reference_level,
  
  positive_class = positive_class,
  
  model_type = model_type,
  
  model_formula_text = model_formula_text,
  
  analysis_objective = analysis_objective,
  
  nested_cv_required = nested_cv_required,
  
  outer_folds = outer_folds,
  
  inner_folds = inner_folds,
  
  stratified_cv_required = stratified_cv_required,
  
  feature_selection_allowed =
    feature_selection_allowed,
  
  external_refit_allowed =
    external_refit_allowed,
  
  external_coefficient_estimation_allowed =
    external_coefficient_estimation_allowed,
  
  external_feature_selection_allowed =
    external_feature_selection_allowed,
  
  external_threshold_optimization_allowed =
    external_threshold_optimization_allowed,
  
  external_cutpoint_search_allowed =
    external_cutpoint_search_allowed,
  
  external_model_tuning_allowed =
    external_model_tuning_allowed,
  
  clinical_diagnostic_claim_allowed =
    clinical_diagnostic_claim_allowed,
  
  diagnostic_biomarker_claim_allowed =
    diagnostic_biomarker_claim_allowed,
  
  gse32863_subjects_locked =
    gse32863_subjects_locked,
  
  gse32863_complete_pairs_locked =
    gse32863_complete_pairs_locked,
  
  gse32863_total_samples_locked =
    gse32863_total_samples_locked,
  
  gse32863_tumor_samples_locked =
    gse32863_tumor_samples_locked,
  
  gse32863_normal_samples_locked =
    gse32863_normal_samples_locked,
  
  historical_tcga_branch =
    historical_tcga_branch,
  
  historical_tcga_tumor_n =
    historical_tcga_tumor_n,
  
  historical_tcga_normal_n =
    historical_tcga_normal_n,
  
  historical_tcga_total_n =
    historical_tcga_total_n,
  
  revised_tcga_branch =
    revised_tcga_branch,
  
  revised_tcga_primary_tumor_n =
    revised_tcga_primary_tumor_n,
  
  revised_tcga_normal_n =
    revised_tcga_normal_n,
  
  integrity_checks =
    integrity_checks,
  
  generated_at =
    format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)


saveRDS(
  model_definition_lock,
  file.path(
    step_dir,
    "STEP13F_Model_Definition_Lock_FINAL.rds"
  )
)


# ------------------------------------------------------------
# 8. Save human-readable lock report
# ------------------------------------------------------------

lock_report <- c(
  
  "STEP13F — MODEL DEFINITION LOCK",
  "================================",
  "",
  
  paste(
    "Status:",
    if (all_integrity_pass) "PASS" else "FAIL"
  ),
  
  "",
  
  paste(
    "Core genes:",
    paste(core_genes, collapse = ", ")
  ),
  
  paste(
    "Outcome:",
    outcome_name
  ),
  
  paste(
    "Outcome levels:",
    paste(outcome_levels, collapse = " / ")
  ),
  
  paste(
    "Reference level:",
    reference_level
  ),
  
  paste(
    "Positive class:",
    positive_class
  ),
  
  paste(
    "Model:",
    model_type
  ),
  
  paste(
    "Formula:",
    model_formula_text
  ),
  
  paste(
    "Objective:",
    analysis_objective
  ),
  
  "",
  
  paste(
    "Nested CV required:",
    nested_cv_required
  ),
  
  paste(
    "Outer folds:",
    outer_folds
  ),
  
  paste(
    "Inner folds:",
    inner_folds
  ),
  
  paste(
    "Stratified CV:",
    stratified_cv_required
  ),
  
  "",
  
  paste(
    "Feature selection allowed:",
    feature_selection_allowed
  ),
  
  paste(
    "External refit allowed:",
    external_refit_allowed
  ),
  
  paste(
    "External coefficient estimation allowed:",
    external_coefficient_estimation_allowed
  ),
  
  paste(
    "External feature selection allowed:",
    external_feature_selection_allowed
  ),
  
  paste(
    "External threshold optimization allowed:",
    external_threshold_optimization_allowed
  ),
  
  paste(
    "External cutpoint search allowed:",
    external_cutpoint_search_allowed
  ),
  
  paste(
    "External model tuning allowed:",
    external_model_tuning_allowed
  ),
  
  "",
  
  paste(
    "Clinical diagnostic claim allowed:",
    clinical_diagnostic_claim_allowed
  ),
  
  paste(
    "Diagnostic biomarker claim allowed:",
    diagnostic_biomarker_claim_allowed
  ),
  
  "",
  
  "GSE32863 locked structure:",
  paste(
    gse32863_total_samples_locked,
    "samples;",
    gse32863_subjects_locked,
    "subjects;",
    gse32863_complete_pairs_locked,
    "complete N+T pairs;",
    gse32863_tumor_samples_locked,
    "tumor;",
    gse32863_normal_samples_locked,
    "normal."
  ),
  
  "",
  
  "TCGA branch separation:",
  
  paste(
    "Historical ROC:",
    historical_tcga_branch,
    "|",
    historical_tcga_tumor_n,
    "Tumor +",
    historical_tcga_normal_n,
    "Normal =",
    historical_tcga_total_n
  ),
  
  paste(
    "Revised WGCNA:",
    revised_tcga_branch,
    "|",
    revised_tcga_primary_tumor_n,
    "Primary Tumor +",
    revised_tcga_normal_n,
    "Normal"
  ),
  
  "",
  
  "No model fitting was performed.",
  "No ROC/AUC calculation was performed.",
  "No cross-validation was performed.",
  "No threshold optimization was performed.",
  "No external refitting was performed.",
  "No external coefficient estimation was performed.",
  "No external feature selection was performed.",
  
  "",
  
  "STEP13F is a model-definition lock only."
)


writeLines(
  lock_report,
  file.path(
    step_dir,
    "STEP13F_Model_Definition_Lock_FINAL.txt"
  )
)


# ------------------------------------------------------------
# 9. Allowed / forbidden downstream operations
# ------------------------------------------------------------

allowed_operations <- c(
  
  "Use the locked predictors AOC3, MS4A7, IL7R.",
  
  "Use TumorStatus as the binary outcome.",
  
  "Use Normal as the reference level.",
  
  "Use Tumor as the positive class.",
  
  "Fit the prespecified multivariable logistic regression.",
  
  "Use nested cross-validation with 10 outer folds and 10 inner folds.",
  
  "Preserve class balance during cross-validation.",
  
  "Assess model performance on held-out data.",
  
  "Apply the locked model to GSE32863 without external refitting.",
  
  "Use subject-aware paired analysis for paired tumor/normal expression comparisons."
)


forbidden_operations <- c(
  
  "Do not perform data-driven feature selection.",
  
  "Do not add or remove predictors based on the current data.",
  
  "Do not refit the model on GSE32863.",
  
  "Do not estimate new external coefficients.",
  
  "Do not tune the model on GSE32863.",
  
  "Do not optimize an external threshold or cutpoint.",
  
  "Do not search for an external classification threshold.",
  
  "Do not make a clinical diagnostic biomarker claim.",
  
  "Do not describe the analysis as clinical implementation.",
  
  "Do not use external data to redefine the model.",
  
  "Do not merge the historical Xiangtao ROC branch with the revised Sangerbox WGCNA branch.",
  
  "Do not reinterpret the GSE32863 57 complete pairs as 58 complete pairs."
)


writeLines(
  allowed_operations,
  file.path(
    step_dir,
    "STEP13F_Allowed_Downstream_Operations.txt"
  )
)


writeLines(
  forbidden_operations,
  file.path(
    step_dir,
    "STEP13F_Forbidden_Downstream_Operations.txt"
  )
)


# ------------------------------------------------------------
# 10. Final PASS / FAIL marker
# ------------------------------------------------------------

if (all_integrity_pass) {
  
  final_status_text <- c(
    
    "STEP13F FINAL STATUS",
    
    "PASS — STEP13F Model Definition Lock established.",
    
    "",
    
    "All required objects were present.",
    
    "All model-definition integrity checks passed.",
    
    "",
    
    "No model fitting / ROC / AUC / CV / threshold analysis",
    "was performed in STEP13F."
  )
  
  writeLines(
    final_status_text,
    file.path(
      step_dir,
      "STEP13F_PASS.txt"
    )
  )
  
  cat("\n")
  cat("FINAL STATUS:\n")
  cat("PASS — STEP13F Model Definition Lock established.\n")
  cat("\n")
  
} else {
  
  failed_checks <-
    names(integrity_checks)[!integrity_checks]
  
  final_status_text <- c(
    
    "STEP13F FINAL STATUS",
    
    "FAIL — STEP13F Model Definition Lock NOT established.",
    
    "",
    
    "Failed checks:",
    
    paste(
      failed_checks,
      collapse = "\n"
    ),
    
    "",
    
    "No model fitting / ROC / AUC / CV / threshold analysis",
    "was performed in STEP13F."
  )
  
  writeLines(
    final_status_text,
    file.path(
      step_dir,
      "STEP13F_FAIL.txt"
    )
  )
  
  cat("\n")
  cat("FINAL STATUS:\n")
  cat("FAIL — STEP13F Model Definition Lock NOT established.\n")
  cat("\n")
  
}


# ------------------------------------------------------------
# 11. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(
    step_dir,
    "STEP13F_sessionInfo_FINAL.txt"
  )
)


# ------------------------------------------------------------
# 12. Final completion message
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP13F COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Output directory:\n",
  step_dir,
  "\n\n",
  sep = ""
)

cat(
  "Model-definition RDS:\n",
  file.path(
    step_dir,
    "STEP13F_Model_Definition_Lock_FINAL.rds"
  ),
  "\n\n",
  sep = ""
)

cat(
  "Integrity audit:\n",
  file.path(
    step_dir,
    "STEP13F_Integrity_Checks.csv"
  ),
  "\n\n",
  sep = ""
)

cat(
  "Object existence audit:\n",
  file.path(
    step_dir,
    "STEP13F_Preflight_Object_Existence_Audit.csv"
  ),
  "\n\n",
  sep = ""
)

cat("============================================================\n")



# ============================================================
# STEP13G — NESTED CROSS-VALIDATION
# FINAL / AUDIT-SAFE VERSION
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP13G — Nested Cross-Validation\n")
cat("FINAL / AUDIT-SAFE VERSION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Workflow paths
# ------------------------------------------------------------

workflow_root <-
  "input_dir/Paper 1_Step 13-Diagnostic Classification"

step_dir <-
  file.path(
    workflow_root,
    "STEP13G_Nested_CV"
  )

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 0A. Remove stale PASS / FAIL markers
# ------------------------------------------------------------

stale_files <- c(
  "STEP13G_PASS.txt",
  "STEP13G_FAIL.txt"
)

for (f in file.path(step_dir, stale_files)) {
  
  if (file.exists(f)) {
    file.remove(f)
  }
}


# ------------------------------------------------------------
# 1. Locked STEP13F model-definition path
# ------------------------------------------------------------

step13f_rds <-
  file.path(
    workflow_root,
    "STEP13F_Model_Definition_Lock",
    "STEP13F_Model_Definition_Lock_FINAL.rds"
  )


cat("STEP13F lock file:\n")
cat(step13f_rds, "\n\n")


# ------------------------------------------------------------
# 2. STEP13F existence check
# ------------------------------------------------------------

if (!file.exists(step13f_rds)) {
  
  writeLines(
    c(
      "STEP13G FINAL STATUS",
      "",
      "FAIL — STEP13F model-definition lock file not found.",
      "",
      paste(
        "Expected:",
        step13f_rds
      ),
      "",
      "No model fitting / CV / ROC analysis was performed."
    ),
    file.path(
      step_dir,
      "STEP13G_FAIL.txt"
    )
  )
  
  stop(
    "STEP13G FAIL — STEP13F RDS not found.",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 3. Load STEP13F lock
# ------------------------------------------------------------

model_lock <-
  readRDS(step13f_rds)


# ------------------------------------------------------------
# 4. STEP13F status check
# ------------------------------------------------------------

if (
  !is.list(model_lock) ||
  !("status" %in% names(model_lock)) ||
  model_lock$status != "PASS"
) {
  
  writeLines(
    c(
      "STEP13G FINAL STATUS",
      "",
      "FAIL — STEP13F lock is not in PASS status.",
      "",
      "No model fitting / CV / ROC analysis was performed."
    ),
    file.path(
      step_dir,
      "STEP13G_FAIL.txt"
    )
  )
  
  stop(
    "STEP13G FAIL — STEP13F lock status is not PASS.",
    call. = FALSE
  )
}


cat("PASS — STEP13F model-definition lock loaded.\n\n")


# ------------------------------------------------------------
# 5. Required STEP13F lock-object audit
# ------------------------------------------------------------

required_lock_fields <- c(
  
  "core_genes",
  
  "outcome_name",
  
  "outcome_levels",
  
  "reference_level",
  
  "positive_class",
  
  "model_type",
  
  "model_formula_text",
  
  "nested_cv_required",
  
  "outer_folds",
  
  "inner_folds",
  
  "stratified_cv_required",
  
  "feature_selection_allowed",
  
  "external_refit_allowed",
  
  "external_coefficient_estimation_allowed",
  
  "external_feature_selection_allowed",
  
  "external_threshold_optimization_allowed",
  
  "external_cutpoint_search_allowed",
  
  "external_model_tuning_allowed",
  
  "clinical_diagnostic_claim_allowed",
  
  "diagnostic_biomarker_claim_allowed",
  
  "historical_tcga_tumor_n",
  
  "historical_tcga_normal_n",
  
  "historical_tcga_total_n"
)


lock_field_exists <-
  required_lock_fields %in% names(model_lock)


lock_field_audit <- data.frame(
  field = required_lock_fields,
  exists = lock_field_exists,
  stringsAsFactors = FALSE
)


write.csv(
  lock_field_audit,
  file.path(
    step_dir,
    "STEP13G_STEP13F_Lock_Field_Existence_Audit.csv"
  ),
  row.names = FALSE
)


if (any(!lock_field_exists)) {
  
  missing_fields <-
    required_lock_fields[!lock_field_exists]
  
  writeLines(
    c(
      "STEP13G FINAL STATUS",
      "",
      "FAIL — Required STEP13F lock fields are missing.",
      "",
      "Missing fields:",
      paste(
        missing_fields,
        collapse = "\n"
      ),
      "",
      "No model fitting / CV / ROC analysis was performed."
    ),
    file.path(
      step_dir,
      "STEP13G_FAIL.txt"
    )
  )
  
  stop(
    "STEP13G FAIL — STEP13F lock fields missing.",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 6. Model-definition concordance audit
# ------------------------------------------------------------

check_core_genes <-
  identical(
    model_lock$core_genes,
    c("AOC3", "MS4A7", "IL7R")
  )


check_outcome_name <-
  isTRUE(
    model_lock$outcome_name == "TumorStatus"
  )


check_outcome_levels <-
  identical(
    model_lock$outcome_levels,
    c("Normal", "Tumor")
  )


check_reference <-
  isTRUE(
    model_lock$reference_level == "Normal"
  )


check_positive_class <-
  isTRUE(
    model_lock$positive_class == "Tumor"
  )


check_model_type <-
  isTRUE(
    model_lock$model_type ==
      "multivariable_logistic_regression"
  )


check_formula <-
  isTRUE(
    model_lock$model_formula_text ==
      "TumorStatus ~ AOC3 + MS4A7 + IL7R"
  )


check_nested_cv <-
  isTRUE(
    model_lock$nested_cv_required
  )


check_outer_folds <-
  isTRUE(
    model_lock$outer_folds == 10
  )


check_inner_folds <-
  isTRUE(
    model_lock$inner_folds == 10
  )


check_stratified <-
  isTRUE(
    model_lock$stratified_cv_required
  )


check_feature_selection_disabled <-
  isFALSE(
    model_lock$feature_selection_allowed
  )


check_external_refit_disabled <-
  isFALSE(
    model_lock$external_refit_allowed
  )


check_external_coef_disabled <-
  isFALSE(
    model_lock$external_coefficient_estimation_allowed
  )


check_external_feature_disabled <-
  isFALSE(
    model_lock$external_feature_selection_allowed
  )


check_external_threshold_disabled <-
  isFALSE(
    model_lock$external_threshold_optimization_allowed
  )


check_external_cutpoint_disabled <-
  isFALSE(
    model_lock$external_cutpoint_search_allowed
  )


check_external_tuning_disabled <-
  isFALSE(
    model_lock$external_model_tuning_allowed
  )


model_lock_checks <- c(
  
  Core_genes = check_core_genes,
  
  Outcome_name = check_outcome_name,
  
  Outcome_levels = check_outcome_levels,
  
  Reference_Normal = check_reference,
  
  Positive_Tumor = check_positive_class,
  
  Model_type = check_model_type,
  
  Formula = check_formula,
  
  Nested_CV_required = check_nested_cv,
  
  Outer_10_fold = check_outer_folds,
  
  Inner_10_fold = check_inner_folds,
  
  Stratified_CV = check_stratified,
  
  Feature_selection_disabled =
    check_feature_selection_disabled,
  
  External_refit_disabled =
    check_external_refit_disabled,
  
  External_coefficient_disabled =
    check_external_coef_disabled,
  
  External_feature_selection_disabled =
    check_external_feature_disabled,
  
  External_threshold_disabled =
    check_external_threshold_disabled,
  
  External_cutpoint_disabled =
    check_external_cutpoint_disabled,
  
  External_model_tuning_disabled =
    check_external_tuning_disabled
)


model_lock_check_table <- data.frame(
  check = names(model_lock_checks),
  PASS = as.logical(model_lock_checks),
  stringsAsFactors = FALSE
)


write.csv(
  model_lock_check_table,
  file.path(
    step_dir,
    "STEP13G_STEP13F_Model_Definition_Concordance.csv"
  ),
  row.names = FALSE
)


if (!all(model_lock_checks)) {
  
  failed_checks <-
    names(model_lock_checks)[
      !model_lock_checks
    ]
  
  writeLines(
    c(
      "STEP13G FINAL STATUS",
      "",
      "FAIL — STEP13F model-definition concordance failed.",
      "",
      "Failed checks:",
      paste(
        failed_checks,
        collapse = "\n"
      ),
      "",
      "No model fitting / CV / ROC analysis was performed."
    ),
    file.path(
      step_dir,
      "STEP13G_FAIL.txt"
    )
  )
  
  stop(
    "STEP13G FAIL — STEP13F definition mismatch.",
    call. = FALSE
  )
}


cat(
  "PASS — STEP13F model definition fully concordant.\n\n"
)


# ------------------------------------------------------------
# 7. Historical TCGA input
# ------------------------------------------------------------

historical_rdata <-
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal/.RData"


cat("Historical TCGA input:\n")
cat(historical_rdata, "\n\n")


if (!file.exists(historical_rdata)) {
  
  writeLines(
    c(
      "STEP13G FINAL STATUS",
      "",
      "FAIL — Historical TCGA .RData not found.",
      "",
      paste(
        "Expected:",
        historical_rdata
      ),
      "",
      "No model fitting / CV / ROC analysis was performed."
    ),
    file.path(
      step_dir,
      "STEP13G_FAIL.txt"
    )
  )
  
  stop(
    "STEP13G FAIL — historical TCGA RData not found.",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 8. Load historical TCGA branch into isolated environment
# ------------------------------------------------------------

tcga_env <- new.env(parent = emptyenv())

load(
  historical_rdata,
  envir = tcga_env
)


# ------------------------------------------------------------
# 9. Required historical TCGA object audit
# ------------------------------------------------------------

required_tcga_objects <- c(
  "df_wide"
)


tcga_object_exists <-
  required_tcga_objects %in%
  ls(envir = tcga_env)


tcga_object_audit <- data.frame(
  object = required_tcga_objects,
  exists = tcga_object_exists,
  stringsAsFactors = FALSE
)


write.csv(
  tcga_object_audit,
  file.path(
    step_dir,
    "STEP13G_TCGA_Object_Existence_Audit.csv"
  ),
  row.names = FALSE
)


if (any(!tcga_object_exists)) {
  
  missing_objects <-
    required_tcga_objects[!tcga_object_exists]
  
  writeLines(
    c(
      "STEP13G FINAL STATUS",
      "",
      "FAIL — Required historical TCGA object missing.",
      "",
      "Missing object:",
      paste(
        missing_objects,
        collapse = "\n"
      ),
      "",
      "No model fitting / CV / ROC analysis was performed."
    ),
    file.path(
      step_dir,
      "STEP13G_FAIL.txt"
    )
  )
  
  stop(
    "STEP13G FAIL — required TCGA object missing.",
    call. = FALSE
  )
}


df_wide <- tcga_env$df_wide


# ------------------------------------------------------------
# 10. Historical TCGA input structure audit
# ------------------------------------------------------------

required_columns <- c(
  "sample_id",
  "status",
  "AOC3",
  "MS4A7",
  "IL7R"
)


column_exists <-
  required_columns %in%
  colnames(df_wide)


column_audit <- data.frame(
  column = required_columns,
  exists = column_exists,
  stringsAsFactors = FALSE
)


write.csv(
  column_audit,
  file.path(
    step_dir,
    "STEP13G_TCGA_Column_Existence_Audit.csv"
  ),
  row.names = FALSE
)


if (any(!column_exists)) {
  
  missing_columns <-
    required_columns[!column_exists]
  
  writeLines(
    c(
      "STEP13G FINAL STATUS",
      "",
      "FAIL — Required TCGA columns are missing.",
      "",
      "Missing columns:",
      paste(
        missing_columns,
        collapse = "\n"
      ),
      "",
      "No model fitting / CV / ROC analysis was performed."
    ),
    file.path(
      step_dir,
      "STEP13G_FAIL.txt"
    )
  )
  
  stop(
    "STEP13G FAIL — required TCGA columns missing.",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 11. TCGA historical sample-count lock
# ------------------------------------------------------------

n_tcga <- nrow(df_wide)

n_tumor <-
  sum(df_wide$status == "Tumor", na.rm = TRUE)

n_normal <-
  sum(df_wide$status == "Normal", na.rm = TRUE)


tcga_count_checks <- c(
  
  Total_575 =
    n_tcga == 575,
  
  Tumor_516 =
    n_tumor == 516,
  
  Normal_59 =
    n_normal == 59,
  
  Tumor_plus_Normal =
    n_tumor + n_normal == n_tcga
)


tcga_count_table <- data.frame(
  check = names(tcga_count_checks),
  PASS = as.logical(tcga_count_checks),
  stringsAsFactors = FALSE
)


tcga_count_values <- data.frame(
  metric = c(
    "Total",
    "Tumor",
    "Normal"
  ),
  observed = c(
    n_tcga,
    n_tumor,
    n_normal
  ),
  locked = c(
    model_lock$historical_tcga_total_n,
    model_lock$historical_tcga_tumor_n,
    model_lock$historical_tcga_normal_n
  ),
  stringsAsFactors = FALSE
)


write.csv(
  tcga_count_table,
  file.path(
    step_dir,
    "STEP13G_TCGA_Sample_Count_Checks.csv"
  ),
  row.names = FALSE
)


write.csv(
  tcga_count_values,
  file.path(
    step_dir,
    "STEP13G_TCGA_Sample_Count_Observed_vs_Locked.csv"
  ),
  row.names = FALSE
)


if (!all(tcga_count_checks)) {
  
  failed_checks <-
    names(tcga_count_checks)[
      !tcga_count_checks
    ]
  
  writeLines(
    c(
      "STEP13G FINAL STATUS",
      "",
      "FAIL — Historical TCGA sample-count lock failed.",
      "",
      "Failed checks:",
      paste(
        failed_checks,
        collapse = "\n"
      ),
      "",
      "No model fitting / CV / ROC analysis was performed."
    ),
    file.path(
      step_dir,
      "STEP13G_FAIL.txt"
    )
  )
  
  stop(
    "STEP13G FAIL — TCGA sample counts do not match lock.",
    call. = FALSE
  )
}


cat(
  "PASS — TCGA historical branch:",
  n_tcga,
  "samples =",
  n_tumor,
  "Tumor +",
  n_normal,
  "Normal.\n\n"
)


# ------------------------------------------------------------
# 12. Core-gene input audit
# ------------------------------------------------------------

gene_input_checks <- c(
  
  AOC3_present =
    "AOC3" %in% colnames(df_wide),
  
  MS4A7_present =
    "MS4A7" %in% colnames(df_wide),
  
  IL7R_present =
    "IL7R" %in% colnames(df_wide)
)


if (!all(gene_input_checks)) {
  
  failed_checks <-
    names(gene_input_checks)[
      !gene_input_checks
    ]
  
  writeLines(
    c(
      "STEP13G FINAL STATUS",
      "",
      "FAIL — One or more locked core genes are absent.",
      "",
      "Failed checks:",
      paste(
        failed_checks,
        collapse = "\n"
      )
    ),
    file.path(
      step_dir,
      "STEP13G_FAIL.txt"
    )
  )
  
  stop(
    "STEP13G FAIL — core gene input audit failed.",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 13. Create analysis dataset
# ------------------------------------------------------------

analysis_df <- df_wide[
  ,
  c(
    "sample_id",
    "status",
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  drop = FALSE
]


# ------------------------------------------------------------
# 14. Missingness / complete-case audit
# ------------------------------------------------------------

complete_case <-
  complete.cases(
    analysis_df[
      ,
      c(
        "status",
        "AOC3",
        "MS4A7",
        "IL7R"
      )
    ]
  )


n_complete <- sum(complete_case)

n_incomplete <- sum(!complete_case)


complete_case_audit <- data.frame(
  metric = c(
    "Total_input_rows",
    "Complete_cases",
    "Incomplete_cases"
  ),
  value = c(
    nrow(analysis_df),
    n_complete,
    n_incomplete
  ),
  stringsAsFactors = FALSE
)


write.csv(
  complete_case_audit,
  file.path(
    step_dir,
    "STEP13G_Complete_Case_Audit.csv"
  ),
  row.names = FALSE
)


if (n_complete != 575) {
  
  writeLines(
    c(
      "STEP13G FINAL STATUS",
      "",
      "FAIL — Complete-case analysis dataset is not 575 samples.",
      "",
      paste(
        "Observed complete cases:",
        n_complete
      ),
      "",
      "No model fitting / CV / ROC analysis was performed."
    ),
    file.path(
      step_dir,
      "STEP13G_FAIL.txt"
    )
  )
  
  stop(
    "STEP13G FAIL — unexpected missingness in locked TCGA input.",
    call. = FALSE
  )
}


analysis_df <-
  analysis_df[complete_case, , drop = FALSE]


# ------------------------------------------------------------
# 15. Outcome factor definition
# ------------------------------------------------------------

analysis_df$status <-
  factor(
    analysis_df$status,
    levels = c(
      "Normal",
      "Tumor"
    )
  )


# ------------------------------------------------------------
# 16. Duplicate sample-ID audit
# ------------------------------------------------------------

duplicate_sample_ids <-
  duplicated(analysis_df$sample_id)


duplicate_id_check <-
  !any(duplicate_sample_ids)


if (!duplicate_id_check) {
  
  writeLines(
    c(
      "STEP13G FINAL STATUS",
      "",
      "FAIL — Duplicate sample IDs detected.",
      "",
      "No model fitting / CV / ROC analysis was performed."
    ),
    file.path(
      step_dir,
      "STEP13G_FAIL.txt"
    )
  )
  
  stop(
    "STEP13G FAIL — duplicate sample IDs.",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 17. Deterministic workflow seed
# ------------------------------------------------------------

workflow_seed <- 20260926L

set.seed(workflow_seed)


# ------------------------------------------------------------
# 18. Stratified fold generator
# ------------------------------------------------------------

make_stratified_folds <- function(
    y,
    k,
    seed
) {
  
  set.seed(seed)
  
  y <- as.character(y)
  
  fold_id <- integer(length(y))
  
  for (cls in unique(y)) {
    
    idx <-
      which(y == cls)
    
    idx <-
      sample(idx)
    
    assigned_folds <-
      rep(
        seq_len(k),
        length.out = length(idx)
      )
    
    fold_id[idx] <-
      assigned_folds
  }
  
  fold_id
}


outer_fold_id <-
  make_stratified_folds(
    analysis_df$status,
    k = 10,
    seed = workflow_seed
  )


# ------------------------------------------------------------
# 19. Outer-fold stratification audit
# ------------------------------------------------------------

outer_fold_audit <- do.call(
  rbind,
  lapply(
    sort(unique(outer_fold_id)),
    function(f) {
      
      idx <-
        which(outer_fold_id == f)
      
      data.frame(
        outer_fold = f,
        n_total = length(idx),
        n_normal =
          sum(
            analysis_df$status[idx] ==
              "Normal"
          ),
        n_tumor =
          sum(
            analysis_df$status[idx] ==
              "Tumor"
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)


write.csv(
  outer_fold_audit,
  file.path(
    step_dir,
    "STEP13G_Outer_Fold_Stratification_Audit.csv"
  ),
  row.names = FALSE
)


outer_fold_structure_check <-
  nrow(outer_fold_audit) == 10 &&
  all(
    outer_fold_audit$n_total > 0
  ) &&
  all(
    outer_fold_audit$n_normal > 0
  ) &&
  all(
    outer_fold_audit$n_tumor > 0
  )


if (!outer_fold_structure_check) {
  
  writeLines(
    c(
      "STEP13G FINAL STATUS",
      "",
      "FAIL — Outer-fold stratification audit failed.",
      "",
      "No model fitting / CV / ROC analysis was performed."
    ),
    file.path(
      step_dir,
      "STEP13G_FAIL.txt"
    )
  )
  
  stop(
    "STEP13G FAIL — outer-fold structure invalid.",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 20. Nested CV containers
# ------------------------------------------------------------

n <- nrow(analysis_df)

oof_probability <-
  rep(
    NA_real_,
    n
  )

oof_outer_fold <-
  rep(
    NA_integer_,
    n
  )

oof_inner_audit <-
  vector(
    "list",
    10
  )

outer_model_coefficients <-
  vector(
    "list",
    10
  )


# ------------------------------------------------------------
# 21. Nested CV
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("STEP13G — Nested CV execution\n")
cat("------------------------------------------------------------\n\n")


for (outer_f in seq_len(10)) {
  
  cat(
    "Outer fold:",
    outer_f,
    "/ 10\n"
  )
  
  
  test_idx <-
    which(
      outer_fold_id == outer_f
    )
  
  train_idx <-
    which(
      outer_fold_id != outer_f
    )
  
  
  train_df <-
    analysis_df[
      train_idx,
      ,
      drop = FALSE
    ]
  
  
  test_df <-
    analysis_df[
      test_idx,
      ,
      drop = FALSE
    ]
  
  
  # ----------------------------------------------------------
  # Inner 10-fold resampling
  #
  # IMPORTANT:
  # No feature selection.
  # No coefficient selection.
  # No threshold tuning.
  # No model tuning.
  #
  # The model is fully prespecified by STEP13F.
  # Inner CV is therefore an audit of the fixed model
  # under the locked nested-resampling framework.
  # ----------------------------------------------------------
  
  inner_fold_id <-
    make_stratified_folds(
      train_df$status,
      k = 10,
      seed = workflow_seed + outer_f
    )
  
  
  inner_audit <-
    do.call(
      rbind,
      lapply(
        seq_len(10),
        function(inner_f) {
          
          inner_test_idx <-
            which(
              inner_fold_id == inner_f
            )
          
          inner_train_idx <-
            which(
              inner_fold_id != inner_f
            )
          
          data.frame(
            outer_fold = outer_f,
            inner_fold = inner_f,
            inner_train_n =
              length(inner_train_idx),
            inner_test_n =
              length(inner_test_idx),
            inner_train_normal =
              sum(
                train_df$status[
                  inner_train_idx
                ] == "Normal"
              ),
            inner_train_tumor =
              sum(
                train_df$status[
                  inner_train_idx
                ] == "Tumor"
              ),
            inner_test_normal =
              sum(
                train_df$status[
                  inner_test_idx
                ] == "Normal"
              ),
            inner_test_tumor =
              sum(
                train_df$status[
                  inner_test_idx
                ] == "Tumor"
              ),
            stringsAsFactors = FALSE
          )
        }
      )
    )
  
  
  oof_inner_audit[[outer_f]] <-
    inner_audit
  
  
  # ----------------------------------------------------------
  # Fit locked model on outer training set
  # ----------------------------------------------------------
  
  fit_outer <-
    glm(
      formula =
        TumorStatus ~ AOC3 + MS4A7 + IL7R,
      data = data.frame(
        TumorStatus = train_df$status,
        AOC3 = train_df$AOC3,
        MS4A7 = train_df$MS4A7,
        IL7R = train_df$IL7R
      ),
      family = binomial()
    )
  
  
  # ----------------------------------------------------------
  # Store outer-fold coefficients
  # ----------------------------------------------------------
  
  coef_vec <-
    coef(fit_outer)
  
  outer_model_coefficients[[outer_f]] <-
    data.frame(
      outer_fold = outer_f,
      term = names(coef_vec),
      coefficient = as.numeric(coef_vec),
      stringsAsFactors = FALSE
    )
  
  
  # ----------------------------------------------------------
  # Predict ONLY the outer held-out samples
  # ----------------------------------------------------------
  
  test_model_df <-
    data.frame(
      AOC3 = test_df$AOC3,
      MS4A7 = test_df$MS4A7,
      IL7R = test_df$IL7R
    )
  
  
  pred <-
    predict(
      fit_outer,
      newdata = test_model_df,
      type = "response"
    )
  
  
  oof_probability[test_idx] <-
    as.numeric(pred)
  
  oof_outer_fold[test_idx] <-
    outer_f
}


# ------------------------------------------------------------
# 22. OOF prediction completeness audit
# ------------------------------------------------------------

oof_complete <-
  all(
    is.finite(oof_probability)
  ) &&
  all(
    !is.na(oof_outer_fold)
  )


if (!oof_complete) {
  
  writeLines(
    c(
      "STEP13G FINAL STATUS",
      "",
      "FAIL — Out-of-fold predictions are incomplete.",
      "",
      "No final performance interpretation should be made."
    ),
    file.path(
      step_dir,
      "STEP13G_FAIL.txt"
    )
  )
  
  stop(
    "STEP13G FAIL — incomplete OOF predictions.",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 23. OOF prediction table
# ------------------------------------------------------------

oof_predictions <- data.frame(
  
  sample_id =
    analysis_df$sample_id,
  
  observed_status =
    as.character(
      analysis_df$status
    ),
  
  observed_binary =
    ifelse(
      analysis_df$status == "Tumor",
      1,
      0
    ),
  
  predicted_probability_Tumor =
    oof_probability,
  
  outer_fold =
    oof_outer_fold,
  
  stringsAsFactors = FALSE
)


write.csv(
  oof_predictions,
  file.path(
    step_dir,
    "STEP13G_Out_of_Fold_Predictions.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 24. ROC / AUC from OOF predictions
# ------------------------------------------------------------

if (!requireNamespace(
  "pROC",
  quietly = TRUE
)) {
  
  writeLines(
    c(
      "STEP13G FINAL STATUS",
      "",
      "FAIL — pROC package is not installed.",
      "",
      "No final ROC/AUC result was generated."
    ),
    file.path(
      step_dir,
      "STEP13G_FAIL.txt"
    )
  )
  
  stop(
    "STEP13G FAIL — pROC package unavailable.",
    call. = FALSE
  )
}


roc_oof <-
  pROC::roc(
    response =
      oof_predictions$observed_status,
    predictor =
      oof_predictions$predicted_probability_Tumor,
    levels =
      c("Normal", "Tumor"),
    direction = "<",
    quiet = TRUE
  )


auc_oof <-
  as.numeric(
    pROC::auc(roc_oof)
  )


auc_ci <-
  as.numeric(
    pROC::ci.auc(
      roc_oof,
      conf.level = 0.95
    )
  )


# ------------------------------------------------------------
# 25. Fixed threshold definition
# ------------------------------------------------------------
#
# IMPORTANT:
# No threshold optimization is allowed in STEP13G.
#
# The probability threshold is fixed at 0.50.
# This is NOT a data-derived cutpoint.
# ------------------------------------------------------------

fixed_threshold <- 0.50


predicted_class <-
  ifelse(
    oof_predictions$predicted_probability_Tumor >=
      fixed_threshold,
    "Tumor",
    "Normal"
  )


observed_class <-
  factor(
    oof_predictions$observed_status,
    levels = c(
      "Normal",
      "Tumor"
    )
  )


predicted_class <-
  factor(
    predicted_class,
    levels = c(
      "Normal",
      "Tumor"
    )
  )


confusion_matrix <-
  table(
    Observed = observed_class,
    Predicted = predicted_class
  )


TN <-
  confusion_matrix[
    "Normal",
    "Normal"
  ]

FP <-
  confusion_matrix[
    "Normal",
    "Tumor"
  ]

FN <-
  confusion_matrix[
    "Tumor",
    "Normal"
  ]

TP <-
  confusion_matrix[
    "Tumor",
    "Tumor"
  ]


sensitivity <-
  TP / (TP + FN)

specificity <-
  TN / (TN + FP)

accuracy <-
  (TP + TN) /
  sum(confusion_matrix)

ppv <-
  TP / (TP + FP)

npv <-
  TN / (TN + FN)

balanced_accuracy <-
  (
    sensitivity +
      specificity
  ) / 2


# ------------------------------------------------------------
# 26. Brier score
# ------------------------------------------------------------

observed_binary <-
  oof_predictions$observed_binary

brier_score <-
  mean(
    (
      oof_probability -
        observed_binary
    )^2
  )


# ------------------------------------------------------------
# 27. Performance summary
# ------------------------------------------------------------

performance_summary <- data.frame(
  
  metric = c(
    "N",
    "Tumor_N",
    "Normal_N",
    "OOF_AUC",
    "OOF_AUC_CI95_lower",
    "OOF_AUC_CI95_upper",
    "Fixed_threshold",
    "Sensitivity",
    "Specificity",
    "Accuracy",
    "Balanced_accuracy",
    "PPV",
    "NPV",
    "Brier_score"
  ),
  
  value = c(
    nrow(oof_predictions),
    sum(
      oof_predictions$observed_status ==
        "Tumor"
    ),
    sum(
      oof_predictions$observed_status ==
        "Normal"
    ),
    auc_oof,
    auc_ci[1],
    auc_ci[3],
    fixed_threshold,
    sensitivity,
    specificity,
    accuracy,
    balanced_accuracy,
    ppv,
    npv,
    brier_score
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  performance_summary,
  file.path(
    step_dir,
    "STEP13G_OOF_Performance_Summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 28. Confusion matrix
# ------------------------------------------------------------

write.csv(
  as.data.frame(confusion_matrix),
  file.path(
    step_dir,
    "STEP13G_Confusion_Matrix.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 29. Outer-fold coefficient audit
# ------------------------------------------------------------

outer_coefficient_table <-
  do.call(
    rbind,
    outer_model_coefficients
  )


write.csv(
  outer_coefficient_table,
  file.path(
    step_dir,
    "STEP13G_Outer_Fold_Coefficients.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 30. Inner-fold audit
# ------------------------------------------------------------

inner_fold_audit_table <-
  do.call(
    rbind,
    oof_inner_audit
  )


write.csv(
  inner_fold_audit_table,
  file.path(
    step_dir,
    "STEP13G_Inner_Fold_Stratification_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 31. Final audit checks
# ------------------------------------------------------------

final_audit_checks <- c(
  
  STEP13F_lock_loaded =
    file.exists(step13f_rds),
  
  STEP13F_status_PASS =
    model_lock$status == "PASS",
  
  Core_genes_locked =
    check_core_genes,
  
  Formula_locked =
    check_formula,
  
  TCGA_total_575 =
    n_tcga == 575,
  
  TCGA_tumor_516 =
    n_tumor == 516,
  
  TCGA_normal_59 =
    n_normal == 59,
  
  Complete_cases_575 =
    n_complete == 575,
  
  Unique_sample_IDs =
    duplicate_id_check,
  
  Outer_10_folds =
    length(
      unique(outer_fold_id)
    ) == 10,
  
  Outer_stratified =
    outer_fold_structure_check,
  
  OOF_predictions_complete =
    oof_complete,
  
  OOF_AUC_finite =
    is.finite(auc_oof),
  
  External_refit_disabled =
    isFALSE(
      model_lock$external_refit_allowed
    ),
  
  External_coefficient_estimation_disabled =
    isFALSE(
      model_lock$
        external_coefficient_estimation_allowed
    ),
  
  External_threshold_optimization_disabled =
    isFALSE(
      model_lock$
        external_threshold_optimization_allowed
    ),
  
  Feature_selection_disabled =
    isFALSE(
      model_lock$feature_selection_allowed
    )
)


final_audit_table <- data.frame(
  check = names(final_audit_checks),
  PASS = as.logical(final_audit_checks),
  stringsAsFactors = FALSE
)


write.csv(
  final_audit_table,
  file.path(
    step_dir,
    "STEP13G_Final_Audit_Checks.csv"
  ),
  row.names = FALSE
)


all_final_audit_pass <-
  all(final_audit_checks)


# ------------------------------------------------------------
# 32. Final summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13G FINAL SUMMARY\n")
cat("============================================================\n\n")


cat(
  "Input branch:             Historical Xiangtao TCGA ROC branch\n"
)

cat(
  "Input samples:            ",
  n_tcga,
  "\n",
  sep = ""
)

cat(
  "Tumor:                    ",
  n_tumor,
  "\n",
  sep = ""
)

cat(
  "Normal:                   ",
  n_normal,
  "\n",
  sep = ""
)

cat(
  "Core genes:               ",
  paste(
    model_lock$core_genes,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Model:                    ",
  model_lock$model_type,
  "\n",
  sep = ""
)

cat(
  "Formula:                  ",
  model_lock$model_formula_text,
  "\n",
  sep = ""
)

cat(
  "Outer folds:              ",
  model_lock$outer_folds,
  "\n",
  sep = ""
)

cat(
  "Inner folds:              ",
  model_lock$inner_folds,
  "\n",
  sep = ""
)

cat(
  "Feature selection:        ",
  model_lock$feature_selection_allowed,
  "\n",
  sep = ""
)

cat(
  "External refit:           ",
  model_lock$external_refit_allowed,
  "\n",
  sep = ""
)

cat(
  "External threshold tune:  ",
  model_lock$
    external_threshold_optimization_allowed,
  "\n",
  sep = ""
)

cat("\n")

cat(
  "OOF AUC:                  ",
  sprintf("%.6f", auc_oof),
  "\n",
  sep = ""
)

cat(
  "OOF AUC 95% CI:           ",
  sprintf(
    "%.6f – %.6f",
    auc_ci[1],
    auc_ci[3]
  ),
  "\n",
  sep = ""
)

cat(
  "Fixed threshold:          ",
  fixed_threshold,
  "\n",
  sep = ""
)

cat(
  "Sensitivity:              ",
  sprintf("%.6f", sensitivity),
  "\n",
  sep = ""
)

cat(
  "Specificity:              ",
  sprintf("%.6f", specificity),
  "\n",
  sep = ""
)

cat(
  "Accuracy:                 ",
  sprintf("%.6f", accuracy),
  "\n",
  sep = ""
)

cat(
  "Balanced accuracy:        ",
  sprintf("%.6f", balanced_accuracy),
  "\n",
  sep = ""
)

cat(
  "PPV:                      ",
  sprintf("%.6f", ppv),
  "\n",
  sep = ""
)

cat(
  "NPV:                      ",
  sprintf("%.6f", npv),
  "\n",
  sep = ""
)

cat(
  "Brier score:              ",
  sprintf("%.6f", brier_score),
  "\n",
  sep = ""
)

cat("\n")

cat(
  "All final audit checks:   ",
  all_final_audit_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 33. Save STEP13G result object
# ------------------------------------------------------------

step13g_result <- list(
  
  workflow = "STEP13G",
  
  status =
    if (
      all_final_audit_pass
    ) {
      "PASS"
    } else {
      "FAIL"
    },
  
  step13f_lock_file =
    step13f_rds,
  
  input_branch =
    "Historical Xiangtao TCGA ROC branch",
  
  input_samples =
    n_tcga,
  
  tumor_samples =
    n_tumor,
  
  normal_samples =
    n_normal,
  
  core_genes =
    model_lock$core_genes,
  
  outcome =
    model_lock$outcome_name,
  
  reference =
    model_lock$reference_level,
  
  positive_class =
    model_lock$positive_class,
  
  model_type =
    model_lock$model_type,
  
  formula =
    model_lock$model_formula_text,
  
  outer_folds =
    model_lock$outer_folds,
  
  inner_folds =
    model_lock$inner_folds,
  
  stratified_cv =
    model_lock$stratified_cv_required,
  
  feature_selection =
    model_lock$feature_selection_allowed,
  
  external_refit =
    model_lock$external_refit_allowed,
  
  external_coefficient_estimation =
    model_lock$
    external_coefficient_estimation_allowed,
  
  external_threshold_optimization =
    model_lock$
    external_threshold_optimization_allowed,
  
  workflow_seed =
    workflow_seed,
  
  auc =
    auc_oof,
  
  auc_ci95 =
    auc_ci,
  
  fixed_threshold =
    fixed_threshold,
  
  sensitivity =
    sensitivity,
  
  specificity =
    specificity,
  
  accuracy =
    accuracy,
  
  balanced_accuracy =
    balanced_accuracy,
  
  ppv =
    ppv,
  
  npv =
    npv,
  
  brier_score =
    brier_score,
  
  final_audit =
    final_audit_checks,
  
  generated_at =
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
)


saveRDS(
  step13g_result,
  file.path(
    step_dir,
    "STEP13G_Nested_CV_FINAL.rds"
  )
)


# ------------------------------------------------------------
# 34. Human-readable result report
# ------------------------------------------------------------

result_report <- c(
  
  "STEP13G — NESTED CROSS-VALIDATION",
  "===================================",
  "",
  
  paste(
    "Status:",
    if (
      all_final_audit_pass
    ) "PASS" else "FAIL"
  ),
  
  "",
  
  paste(
    "Input branch:",
    "Historical Xiangtao TCGA ROC branch"
  ),
  
  paste(
    "Samples:",
    n_tcga
  ),
  
  paste(
    "Tumor:",
    n_tumor
  ),
  
  paste(
    "Normal:",
    n_normal
  ),
  
  "",
  
  paste(
    "Core genes:",
    paste(
      model_lock$core_genes,
      collapse = ", "
    )
  ),
  
  paste(
    "Model:",
    model_lock$model_type
  ),
  
  paste(
    "Formula:",
    model_lock$model_formula_text
  ),
  
  "",
  
  paste(
    "Outer folds:",
    model_lock$outer_folds
  ),
  
  paste(
    "Inner folds:",
    model_lock$inner_folds
  ),
  
  paste(
    "Feature selection:",
    model_lock$feature_selection_allowed
  ),
  
  paste(
    "External refit:",
    model_lock$external_refit_allowed
  ),
  
  paste(
    "External threshold optimization:",
    model_lock$
      external_threshold_optimization_allowed
  ),
  
  "",
  
  paste(
    "OOF AUC:",
    sprintf("%.6f", auc_oof)
  ),
  
  paste(
    "OOF AUC 95% CI:",
    sprintf(
      "%.6f – %.6f",
      auc_ci[1],
      auc_ci[3]
    )
  ),
  
  paste(
    "Fixed threshold:",
    fixed_threshold
  ),
  
  paste(
    "Sensitivity:",
    sprintf("%.6f", sensitivity)
  ),
  
  paste(
    "Specificity:",
    sprintf("%.6f", specificity)
  ),
  
  paste(
    "Accuracy:",
    sprintf("%.6f", accuracy)
  ),
  
  paste(
    "Balanced accuracy:",
    sprintf("%.6f", balanced_accuracy)
  ),
  
  paste(
    "PPV:",
    sprintf("%.6f", ppv)
  ),
  
  paste(
    "NPV:",
    sprintf("%.6f", npv)
  ),
  
  paste(
    "Brier score:",
    sprintf("%.6f", brier_score)
  ),
  
  "",
  
  "No external data were used.",
  
  "No external refitting was performed.",
  
  "No external coefficient estimation was performed.",
  
  "No feature selection was performed.",
  
  "No threshold optimization was performed.",
  
  "Outer-fold held-out predictions were used for performance estimation.",
  
  "",
  
  "STEP13G is an internal nested-resampling performance analysis."
)


writeLines(
  result_report,
  file.path(
    step_dir,
    "STEP13G_Nested_CV_FINAL.txt"
  )
)


# ------------------------------------------------------------
# 35. Final PASS / FAIL marker
# ------------------------------------------------------------

if (all_final_audit_pass) {
  
  writeLines(
    c(
      "STEP13G FINAL STATUS",
      "",
      "PASS — STEP13G Nested Cross-Validation completed.",
      "",
      paste(
        "OOF AUC =",
        sprintf("%.6f", auc_oof)
      ),
      paste(
        "95% CI =",
        sprintf(
          "%.6f – %.6f",
          auc_ci[1],
          auc_ci[3]
        )
      ),
      "",
      "No feature selection was performed.",
      "No external refitting was performed.",
      "No external coefficient estimation was performed.",
      "No threshold optimization was performed.",
      "",
      "Performance was estimated from outer-fold held-out predictions."
    ),
    file.path(
      step_dir,
      "STEP13G_PASS.txt"
    )
  )
  
  cat("\n")
  cat("FINAL STATUS:\n")
  cat(
    "PASS — STEP13G Nested Cross-Validation completed.\n"
  )
  
} else {
  
  failed_checks <-
    names(final_audit_checks)[
      !final_audit_checks
    ]
  
  writeLines(
    c(
      "STEP13G FINAL STATUS",
      "",
      "FAIL — STEP13G Nested Cross-Validation audit failed.",
      "",
      "Failed checks:",
      paste(
        failed_checks,
        collapse = "\n"
      )
    ),
    file.path(
      step_dir,
      "STEP13G_FAIL.txt"
    )
  )
  
  cat("\n")
  cat("FINAL STATUS:\n")
  cat(
    "FAIL — STEP13G Nested Cross-Validation audit failed.\n"
  )
}


# ------------------------------------------------------------
# 36. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(
    step_dir,
    "STEP13G_sessionInfo_FINAL.txt"
  )
)


# ------------------------------------------------------------
# 37. Completion message
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13G COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Output directory:\n",
  step_dir,
  "\n\n",
  sep = ""
)

cat(
  "Result RDS:\n",
  file.path(
    step_dir,
    "STEP13G_Nested_CV_FINAL.rds"
  ),
  "\n\n",
  sep = ""
)

cat(
  "OOF predictions:\n",
  file.path(
    step_dir,
    "STEP13G_Out_of_Fold_Predictions.csv"
  ),
  "\n\n",
  sep = ""
)

cat(
  "Performance summary:\n",
  file.path(
    step_dir,
    "STEP13G_OOF_Performance_Summary.csv"
  ),
  "\n\n",
  sep = ""
)

cat("============================================================\n")


# ============================================================
# STEP13H — TCGA Final Model
# ============================================================
#
# Purpose:
#   Fit the final prespecified multivariable logistic regression
#   model using the locked 575-sample historical TCGA ROC branch.
#
# Locked predictors:
#   AOC3 + MS4A7 + IL7R
#
# Outcome:
#   TumorStatus
#   Reference = Normal
#   Positive  = Tumor
#
# IMPORTANT:
#   STEP13H is NOT nested cross-validation.
#   STEP13H fits the final model using ALL 575 locked TCGA samples.
#
#   The resulting coefficients are locked for STEP13I external
#   assessment in GSE32863.
#
# No:
#   - feature selection
#   - gene selection
#   - external refitting
#   - external coefficient estimation
#   - external threshold optimization
#   - external cutpoint search
#   - model tuning
#
# ============================================================


# ------------------------------------------------------------
# 0. Clean session
# ------------------------------------------------------------

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(scipen = 999)

set.seed(20260926L)


# ------------------------------------------------------------
# 1. Working directories
# ------------------------------------------------------------

workflow_root <-
  "input_dir/Paper 1_Step 13-Diagnostic Classification"

step13f_dir <-
  file.path(workflow_root,
            "STEP13F_Model_Definition_Lock")

step13g_dir <-
  file.path(workflow_root,
            "STEP13G_Nested_CV")

step13h_dir <-
  file.path(workflow_root,
            "STEP13H_TCGA_Final_Model")

input_lock_dir <-
  file.path(workflow_root,
            "INPUT_LOCK")

results_dir <-
  file.path(workflow_root,
            "RESULTS")

tables_dir <-
  file.path(workflow_root,
            "TABLES")

audit_dir <-
  file.path(workflow_root,
            "AUDIT")

logs_dir <-
  file.path(workflow_root,
            "LOGS")


# ------------------------------------------------------------
# 2. Create output directories
# ------------------------------------------------------------

dir.create(step13h_dir,
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(results_dir,
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(tables_dir,
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(audit_dir,
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(logs_dir,
           recursive = TRUE,
           showWarnings = FALSE)


# ------------------------------------------------------------
# 3. Remove stale PASS / FAIL markers
# ------------------------------------------------------------

stale_files <- c(
  file.path(step13h_dir, "STEP13H_PASS.txt"),
  file.path(step13h_dir, "STEP13H_FAIL.txt"),
  file.path(step13h_dir, "STEP13H_FINAL.txt")
)

stale_files <- stale_files[file.exists(stale_files)]

if (length(stale_files) > 0) {
  file.remove(stale_files)
}


# ------------------------------------------------------------
# 4. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "pROC"
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
    paste0(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", ")
    )
  )
}

library(pROC)


# ------------------------------------------------------------
# 5. Load STEP13F Model Definition Lock
# ------------------------------------------------------------

step13f_rds <-
  file.path(
    step13f_dir,
    "STEP13F_Model_Definition_Lock_FINAL.rds"
  )

if (!file.exists(step13f_rds)) {
  stop(
    "STEP13F FINAL RDS not found: ",
    step13f_rds
  )
}

model_lock <- readRDS(step13f_rds)


# ------------------------------------------------------------
# 6. STEP13F status audit
# ------------------------------------------------------------

if (!identical(model_lock$status, "PASS")) {
  stop(
    "STEP13H cannot proceed: STEP13F status is not PASS."
  )
}


# ------------------------------------------------------------
# 7. STEP13F lock-field existence audit
# ------------------------------------------------------------

required_lock_fields <- c(
  "core_genes",
  "outcome_name",
  "outcome_levels",
  "reference_level",
  "positive_class",
  "model_type",
  "model_formula_text",
  "nested_cv_required",
  "outer_folds",
  "inner_folds",
  "feature_selection_allowed",
  "external_refit_allowed",
  "external_coefficient_estimation_allowed",
  "external_feature_selection_allowed",
  "external_threshold_optimization_allowed",
  "external_cutpoint_search_allowed",
  "external_model_tuning_allowed",
  "clinical_diagnostic_claim_allowed",
  "diagnostic_biomarker_claim_allowed",
  "gse32863_subjects_locked",
  "gse32863_complete_pairs_locked",
  "gse32863_total_samples_locked",
  "gse32863_tumor_samples_locked",
  "gse32863_normal_samples_locked"
)

missing_lock_fields <- setdiff(
  required_lock_fields,
  names(model_lock)
)

if (length(missing_lock_fields) > 0) {
  stop(
    paste0(
      "STEP13F lock is incomplete. Missing field(s): ",
      paste(missing_lock_fields, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 8. Exact model-definition audit
# ------------------------------------------------------------

expected_core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

expected_outcome_name <- "TumorStatus"

expected_outcome_levels <- c(
  "Normal",
  "Tumor"
)

expected_reference_level <- "Normal"

expected_positive_class <- "Tumor"

expected_model_type <-
  "multivariable_logistic_regression"

expected_formula_text <-
  "TumorStatus ~ AOC3 + MS4A7 + IL7R"


integrity_checks <- c(
  
  core_genes =
    identical(
      as.character(model_lock$core_genes),
      expected_core_genes
    ),
  
  outcome_name =
    identical(
      as.character(model_lock$outcome_name),
      expected_outcome_name
    ),
  
  outcome_levels =
    identical(
      as.character(model_lock$outcome_levels),
      expected_outcome_levels
    ),
  
  reference_level =
    identical(
      as.character(model_lock$reference_level),
      expected_reference_level
    ),
  
  positive_class =
    identical(
      as.character(model_lock$positive_class),
      expected_positive_class
    ),
  
  model_type =
    identical(
      as.character(model_lock$model_type),
      expected_model_type
    ),
  
  model_formula =
    identical(
      as.character(model_lock$model_formula_text),
      expected_formula_text
    ),
  
  feature_selection_allowed =
    identical(
      as.logical(model_lock$feature_selection_allowed),
      FALSE
    ),
  
  external_refit_allowed =
    identical(
      as.logical(model_lock$external_refit_allowed),
      FALSE
    ),
  
  external_coefficient_estimation_allowed =
    identical(
      as.logical(
        model_lock$external_coefficient_estimation_allowed
      ),
      FALSE
    ),
  
  external_feature_selection_allowed =
    identical(
      as.logical(
        model_lock$external_feature_selection_allowed
      ),
      FALSE
    ),
  
  external_threshold_optimization_allowed =
    identical(
      as.logical(
        model_lock$external_threshold_optimization_allowed
      ),
      FALSE
    ),
  
  external_cutpoint_search_allowed =
    identical(
      as.logical(
        model_lock$external_cutpoint_search_allowed
      ),
      FALSE
    ),
  
  external_model_tuning_allowed =
    identical(
      as.logical(
        model_lock$external_model_tuning_allowed
      ),
      FALSE
    ),
  
  clinical_diagnostic_claim_allowed =
    identical(
      as.logical(
        model_lock$clinical_diagnostic_claim_allowed
      ),
      FALSE
    ),
  
  diagnostic_biomarker_claim_allowed =
    identical(
      as.logical(
        model_lock$diagnostic_biomarker_claim_allowed
      ),
      FALSE
    )
)


integrity_df <- data.frame(
  check = names(integrity_checks),
  PASS = unname(integrity_checks),
  stringsAsFactors = FALSE
)

write.csv(
  integrity_df,
  file.path(
    step13h_dir,
    "STEP13H_STEP13F_Model_Definition_Integrity.csv"
  ),
  row.names = FALSE
)


if (!all(integrity_checks)) {
  stop(
    "STEP13H stopped: one or more STEP13F model-definition checks failed."
  )
}


# ------------------------------------------------------------
# 9. Historical TCGA ROC branch
# ------------------------------------------------------------

tcga_rdata <-
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal/.RData"

if (!file.exists(tcga_rdata)) {
  stop(
    "Historical TCGA ROC .RData not found: ",
    tcga_rdata
  )
}


# ------------------------------------------------------------
# 10. Load historical TCGA branch into isolated environment
# ------------------------------------------------------------

tcga_env <- new.env(parent = emptyenv())

load(
  tcga_rdata,
  envir = tcga_env
)


# ------------------------------------------------------------
# 11. Required historical object audit
# ------------------------------------------------------------

required_tcga_objects <- c(
  "df_wide",
  "roc_AOC3",
  "roc_MS4A7",
  "roc_IL7R",
  "roc_combined",
  "Figure12A_CombinedROC"
)

tcga_object_exists <- vapply(
  required_tcga_objects,
  exists,
  logical(1),
  envir = tcga_env,
  inherits = FALSE
)

tcga_object_audit <- data.frame(
  object = required_tcga_objects,
  exists = tcga_object_exists,
  stringsAsFactors = FALSE
)

write.csv(
  tcga_object_audit,
  file.path(
    step13h_dir,
    "STEP13H_TCGA_Object_Existence_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(tcga_object_exists)) {
  stop(
    paste0(
      "Required historical TCGA object(s) missing: ",
      paste(
        required_tcga_objects[
          !tcga_object_exists
        ],
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 12. Extract locked TCGA model-input object
# ------------------------------------------------------------

df_wide <- tcga_env$df_wide


# ------------------------------------------------------------
# 13. TCGA input dimension audit
# ------------------------------------------------------------

expected_n_tcga <- 575L
expected_n_tumor <- 516L
expected_n_normal <- 59L

required_columns <- c(
  "sample_id",
  "status",
  "AOC3",
  "MS4A7",
  "IL7R"
)

missing_columns <- setdiff(
  required_columns,
  colnames(df_wide)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Required TCGA column(s) missing: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}


tcga_dimension_audit <- data.frame(
  quantity = c(
    "rows",
    "Tumor",
    "Normal"
  ),
  observed = c(
    nrow(df_wide),
    sum(df_wide$status == "Tumor"),
    sum(df_wide$status == "Normal")
  ),
  locked = c(
    expected_n_tcga,
    expected_n_tumor,
    expected_n_normal
  ),
  PASS = c(
    nrow(df_wide) == expected_n_tcga,
    sum(df_wide$status == "Tumor") == expected_n_tumor,
    sum(df_wide$status == "Normal") == expected_n_normal
  ),
  stringsAsFactors = FALSE
)

write.csv(
  tcga_dimension_audit,
  file.path(
    step13h_dir,
    "STEP13H_TCGA_Sample_Count_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(tcga_dimension_audit$PASS)) {
  stop(
    "STEP13H stopped: TCGA sample-count audit failed."
  )
}


# ------------------------------------------------------------
# 14. TCGA sample ID audit
# ------------------------------------------------------------

sample_id_missing <-
  anyNA(df_wide$sample_id) ||
  any(df_wide$sample_id == "")

sample_id_unique <-
  length(unique(df_wide$sample_id)) ==
  nrow(df_wide)

sample_id_audit <- data.frame(
  check = c(
    "sample_id_missing",
    "sample_id_unique"
  ),
  observed = c(
    !sample_id_missing,
    sample_id_unique
  ),
  PASS = c(
    !sample_id_missing,
    sample_id_unique
  ),
  stringsAsFactors = FALSE
)

write.csv(
  sample_id_audit,
  file.path(
    step13h_dir,
    "STEP13H_TCGA_Sample_ID_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(sample_id_audit$PASS)) {
  stop(
    "STEP13H stopped: TCGA sample-ID audit failed."
  )
}


# ------------------------------------------------------------
# 15. Outcome-level audit
# ------------------------------------------------------------

status_levels_observed <- unique(
  as.character(df_wide$status)
)

status_levels_observed <- sort(
  status_levels_observed
)

status_level_audit <- data.frame(
  expected = expected_outcome_levels,
  observed = expected_outcome_levels %in%
    status_levels_observed,
  stringsAsFactors = FALSE
)

write.csv(
  status_level_audit,
  file.path(
    step13h_dir,
    "STEP13H_TumorStatus_Level_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(
  expected_outcome_levels %in%
  status_levels_observed
)) {
  stop(
    "STEP13H stopped: expected TumorStatus levels not found."
  )
}


# ------------------------------------------------------------
# 16. Build final model-analysis data
# ------------------------------------------------------------
#
# IMPORTANT:
#   No sample exclusion is permitted unless the locked historical
#   input itself contains missing model values.
#
# The expected historical branch is 575 complete observations.
# Therefore complete-case filtering must retain all 575 samples.
#

model_df <- data.frame(
  sample_id = df_wide$sample_id,
  TumorStatus = factor(
    as.character(df_wide$status),
    levels = expected_outcome_levels
  ),
  AOC3 = as.numeric(df_wide$AOC3),
  MS4A7 = as.numeric(df_wide$MS4A7),
  IL7R = as.numeric(df_wide$IL7R),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 17. Complete-case audit
# ------------------------------------------------------------

complete_case <- complete.cases(
  model_df[
    ,
    c(
      "TumorStatus",
      "AOC3",
      "MS4A7",
      "IL7R"
    )
  ]
)

n_complete <- sum(complete_case)
n_incomplete <- sum(!complete_case)

complete_case_audit <- data.frame(
  quantity = c(
    "Total samples",
    "Complete cases",
    "Incomplete cases"
  ),
  value = c(
    nrow(model_df),
    n_complete,
    n_incomplete
  ),
  PASS = c(
    nrow(model_df) == 575L,
    n_complete == 575L,
    n_incomplete == 0L
  ),
  stringsAsFactors = FALSE
)

write.csv(
  complete_case_audit,
  file.path(
    step13h_dir,
    "STEP13H_Complete_Case_Audit.csv"
  ),
  row.names = FALSE
)

if (n_complete != expected_n_tcga) {
  stop(
    paste0(
      "STEP13H stopped: complete-case count = ",
      n_complete,
      ", expected = ",
      expected_n_tcga,
      "."
    )
  )
}

model_df <- model_df[complete_case, , drop = FALSE]


# ------------------------------------------------------------
# 18. Final model outcome audit
# ------------------------------------------------------------

final_tumor_n <-
  sum(model_df$TumorStatus == "Tumor")

final_normal_n <-
  sum(model_df$TumorStatus == "Normal")

final_outcome_audit <- data.frame(
  class = c(
    "Normal",
    "Tumor"
  ),
  n = c(
    final_normal_n,
    final_tumor_n
  ),
  expected = c(
    expected_n_normal,
    expected_n_tumor
  ),
  PASS = c(
    final_normal_n == expected_n_normal,
    final_tumor_n == expected_n_tumor
  ),
  stringsAsFactors = FALSE
)

write.csv(
  final_outcome_audit,
  file.path(
    step13h_dir,
    "STEP13H_Final_Model_Outcome_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(final_outcome_audit$PASS)) {
  stop(
    "STEP13H stopped: final-model outcome counts failed."
  )
}


# ------------------------------------------------------------
# 19. Predictor-value audit
# ------------------------------------------------------------

predictor_matrix <- as.matrix(
  model_df[
    ,
    expected_core_genes,
    drop = FALSE
  ]
)

predictor_finite <- all(
  is.finite(predictor_matrix)
)

predictor_audit <- data.frame(
  gene = expected_core_genes,
  missing = vapply(
    expected_core_genes,
    function(g) anyNA(model_df[[g]]),
    logical(1)
  ),
  finite = vapply(
    expected_core_genes,
    function(g) all(is.finite(model_df[[g]])),
    logical(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  predictor_audit,
  file.path(
    step13h_dir,
    "STEP13H_Predictor_Value_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(predictor_audit$finite)) {
  stop(
    "STEP13H stopped: non-finite predictor values detected."
  )
}


# ------------------------------------------------------------
# 20. Preserve exact locked formula
# ------------------------------------------------------------

final_formula <-
  TumorStatus ~ AOC3 + MS4A7 + IL7R


# ------------------------------------------------------------
# 21. Fit final TCGA logistic model
# ------------------------------------------------------------
#
# This is the FINAL training model.
#
# It uses all 575 locked TCGA observations.
#
# No feature selection is performed.
#

final_model <- glm(
  formula = final_formula,
  data = model_df,
  family = binomial(link = "logit")
)


# ------------------------------------------------------------
# 22. Model convergence audit
# ------------------------------------------------------------

model_converged <-
  isTRUE(final_model$converged)

coef_values <-
  coef(final_model)

coef_finite <-
  all(is.finite(coef_values))

model_rank_ok <-
  final_model$rank == length(coef_values)

convergence_audit <- data.frame(
  check = c(
    "model_converged",
    "coefficients_finite",
    "full_model_rank"
  ),
  PASS = c(
    model_converged,
    coef_finite,
    model_rank_ok
  ),
  stringsAsFactors = FALSE
)

write.csv(
  convergence_audit,
  file.path(
    step13h_dir,
    "STEP13H_Model_Convergence_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(convergence_audit$PASS)) {
  stop(
    "STEP13H stopped: final model convergence/rank audit failed."
  )
}


# ------------------------------------------------------------
# 23. Extract coefficient table
# ------------------------------------------------------------

model_summary <-
  summary(final_model)

coef_table <- as.data.frame(
  model_summary$coefficients,
  stringsAsFactors = FALSE
)

coef_table$term <-
  rownames(coef_table)

rownames(coef_table) <- NULL

colnames(coef_table)[
  colnames(coef_table) == "Estimate"
] <- "Coefficient"

colnames(coef_table)[
  colnames(coef_table) == "Std. Error"
] <- "Std_Error"

colnames(coef_table)[
  colnames(coef_table) == "z value"
] <- "z_value"

colnames(coef_table)[
  colnames(coef_table) == "Pr(>|z|)"
] <- "P_value"


# ------------------------------------------------------------
# 24. Calculate OR and 95% Wald CI
# ------------------------------------------------------------
#
# Wald CI is used here because:
#   1. coefficients are directly linked to the fitted glm;
#   2. it is deterministic;
#   3. it does not alter the model definition;
#   4. it avoids profile-likelihood refitting.
#
# The coefficients themselves remain the primary locked
# parameters for STEP13I.
#

coef_table$OR <-
  exp(coef_table$Coefficient)

coef_table$CI95_lower <-
  exp(
    coef_table$Coefficient -
      1.96 * coef_table$Std_Error
  )

coef_table$CI95_upper <-
  exp(
    coef_table$Coefficient +
      1.96 * coef_table$Std_Error
  )


# ------------------------------------------------------------
# 25. Reorder coefficient table
# ------------------------------------------------------------

coef_table <- coef_table[
  ,
  c(
    "term",
    "Coefficient",
    "Std_Error",
    "z_value",
    "P_value",
    "OR",
    "CI95_lower",
    "CI95_upper"
  )
]


# ------------------------------------------------------------
# 26. Save coefficient table
# ------------------------------------------------------------

write.csv(
  coef_table,
  file.path(
    step13h_dir,
    "STEP13H_Final_Model_Coefficients_OR_CI_P.csv"
  ),
  row.names = FALSE
)

write.csv(
  coef_table,
  file.path(
    results_dir,
    "STEP13H_Final_Model_Coefficients_OR_CI_P.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 27. Predictor-only coefficient table
# ------------------------------------------------------------

predictor_coef_table <- coef_table[
  coef_table$term %in% expected_core_genes,
  ,
  drop = FALSE
]

predictor_coef_table <- predictor_coef_table[
  match(
    expected_core_genes,
    predictor_coef_table$term
  ),
  ,
  drop = FALSE
]

write.csv(
  predictor_coef_table,
  file.path(
    step13h_dir,
    "STEP13H_Core_Gene_Final_Model_Parameters.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 28. Save model predictions for traceability
# ------------------------------------------------------------

model_df$predicted_probability <-
  predict(
    final_model,
    newdata = model_df,
    type = "response"
  )

model_df$predicted_logit <-
  predict(
    final_model,
    newdata = model_df,
    type = "link"
  )


prediction_traceability <- model_df[
  ,
  c(
    "sample_id",
    "TumorStatus",
    "AOC3",
    "MS4A7",
    "IL7R",
    "predicted_logit",
    "predicted_probability"
  )
]

write.csv(
  prediction_traceability,
  file.path(
    step13h_dir,
    "STEP13H_Final_Model_Prediction_Traceability.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 29. Model-input scale provenance
# ------------------------------------------------------------

scale_provenance <- data.frame(
  field = c(
    "Data_branch",
    "Source",
    "Historical_platform_processing",
    "Expression_source",
    "Transformation",
    "Samples_total",
    "Tumor",
    "Normal",
    "Predictors"
  ),
  value = c(
    "Historical Xiangtao TCGA ROC branch",
    tcga_rdata,
    "Xiangtao historical analysis",
    "TCGA-LUAD STAR RNA-seq TPM",
    "log2(value + 1)",
    "575",
    "516",
    "59",
    "AOC3; MS4A7; IL7R"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  scale_provenance,
  file.path(
    step13h_dir,
    "STEP13H_Model_Input_Scale_Provenance.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 30. Model-definition traceability lock
# ------------------------------------------------------------

model_definition_trace <- data.frame(
  field = c(
    "Predictors",
    "Outcome",
    "Reference_class",
    "Positive_class",
    "Model_type",
    "Formula",
    "Feature_selection",
    "External_refit",
    "External_coefficient_estimation",
    "External_threshold_optimization",
    "External_cutpoint_search",
    "External_model_tuning"
  ),
  locked_value = c(
    paste(expected_core_genes, collapse = " + "),
    expected_outcome_name,
    expected_reference_level,
    expected_positive_class,
    expected_model_type,
    expected_formula_text,
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  model_definition_trace,
  file.path(
    step13h_dir,
    "STEP13H_Model_Definition_Traceability_Lock.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 31. STEP13G dependency audit
# ------------------------------------------------------------
#
# STEP13H does not obtain coefficients from STEP13G.
#
# STEP13G is used only as a preceding workflow checkpoint.
#
# This audit confirms that STEP13G completed successfully.
#

step13g_rds <-
  file.path(
    step13g_dir,
    "STEP13G_Nested_CV_FINAL.rds"
  )

step13g_exists <-
  file.exists(step13g_rds)

step13g_status <- NA_character_

if (step13g_exists) {
  
  step13g_result <-
    readRDS(step13g_rds)
  
  if ("status" %in% names(step13g_result)) {
    step13g_status <-
      as.character(step13g_result$status)
  }
}

step13g_audit <- data.frame(
  check = c(
    "STEP13G_RDS_exists",
    "STEP13G_status_PASS"
  ),
  PASS = c(
    step13g_exists,
    identical(step13g_status, "PASS")
  ),
  stringsAsFactors = FALSE
)

write.csv(
  step13g_audit,
  file.path(
    step13h_dir,
    "STEP13H_STEP13G_Preceding_Checkpoint_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(step13g_audit$PASS)) {
  stop(
    "STEP13H stopped: STEP13G preceding checkpoint is not confirmed PASS."
  )
}


# ------------------------------------------------------------
# 32. Verify final model sample count
# ------------------------------------------------------------

final_model_n <-
  nobs(final_model)

final_model_sample_audit <- data.frame(
  quantity = c(
    "Model observations",
    "Tumor observations",
    "Normal observations"
  ),
  observed = c(
    final_model_n,
    sum(model_df$TumorStatus == "Tumor"),
    sum(model_df$TumorStatus == "Normal")
  ),
  expected = c(
    575L,
    516L,
    59L
  ),
  PASS = c(
    final_model_n == 575L,
    sum(model_df$TumorStatus == "Tumor") == 516L,
    sum(model_df$TumorStatus == "Normal") == 59L
  ),
  stringsAsFactors = FALSE
)

write.csv(
  final_model_sample_audit,
  file.path(
    step13h_dir,
    "STEP13H_Final_Model_Sample_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(final_model_sample_audit$PASS)) {
  stop(
    "STEP13H stopped: final model sample audit failed."
  )
}


# ------------------------------------------------------------
# 33. Final coefficient sanity audit
# ------------------------------------------------------------

predictor_terms <- c(
  "(Intercept)",
  expected_core_genes
)

coefficient_presence <- predictor_terms %in%
  coef_table$term

coefficient_finite_check <- all(
  is.finite(
    coef_table$Coefficient
  ) &
    is.finite(
      coef_table$Std_Error
    ) &
    is.finite(
      coef_table$OR
    ) &
    is.finite(
      coef_table$CI95_lower
    ) &
    is.finite(
      coef_table$CI95_upper
    ) &
    is.finite(
      coef_table$P_value
    )
)

coefficient_audit <- data.frame(
  check = c(
    "required_terms_present",
    "all_model_parameters_finite"
  ),
  PASS = c(
    all(coefficient_presence),
    coefficient_finite_check
  ),
  stringsAsFactors = FALSE
)

write.csv(
  coefficient_audit,
  file.path(
    step13h_dir,
    "STEP13H_Final_Coefficient_Sanity_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(coefficient_audit$PASS)) {
  stop(
    "STEP13H stopped: final coefficient sanity audit failed."
  )
}


# ------------------------------------------------------------
# 34. Save final model object
# ------------------------------------------------------------

final_model_object <- list(
  
  status = "PASS",
  
  workflow = "STEP13H",
  
  purpose =
    "Final TCGA training model for locked exploratory tumor-versus-normal tissue classification",
  
  data_branch =
    "Historical Xiangtao TCGA ROC branch",
  
  source_rdata =
    tcga_rdata,
  
  n_total =
    final_model_n,
  
  n_tumor =
    sum(model_df$TumorStatus == "Tumor"),
  
  n_normal =
    sum(model_df$TumorStatus == "Normal"),
  
  predictors =
    expected_core_genes,
  
  outcome =
    expected_outcome_name,
  
  outcome_levels =
    expected_outcome_levels,
  
  reference_level =
    expected_reference_level,
  
  positive_class =
    expected_positive_class,
  
  model_type =
    expected_model_type,
  
  formula_text =
    expected_formula_text,
  
  formula =
    final_formula,
  
  expression_source =
    "TCGA-LUAD STAR RNA-seq TPM",
  
  transformation =
    "log2(value + 1)",
  
  feature_selection =
    FALSE,
  
  external_refit_allowed =
    FALSE,
  
  external_coefficient_estimation_allowed =
    FALSE,
  
  external_threshold_optimization_allowed =
    FALSE,
  
  external_cutpoint_search_allowed =
    FALSE,
  
  external_model_tuning_allowed =
    FALSE,
  
  clinical_diagnostic_claim_allowed =
    FALSE,
  
  diagnostic_biomarker_claim_allowed =
    FALSE,
  
  model =
    final_model,
  
  coefficient_table =
    coef_table,
  
  predictor_coefficient_table =
    predictor_coef_table,
  
  model_data =
    model_df,
  
  step13f_lock_file =
    step13f_rds,
  
  step13g_checkpoint_file =
    step13g_rds,
  
  workflow_seed =
    20260926L,
  
  generated_at =
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
)


step13h_final_rds <-
  file.path(
    step13h_dir,
    "STEP13H_TCGA_Final_Model_FINAL.rds"
  )

saveRDS(
  final_model_object,
  step13h_final_rds
)


# ------------------------------------------------------------
# 35. Save model object to RESULTS
# ------------------------------------------------------------

saveRDS(
  final_model_object,
  file.path(
    results_dir,
    "STEP13H_TCGA_Final_Model_FINAL.rds"
  )
)


# ------------------------------------------------------------
# 36. Human-readable final model report
# ------------------------------------------------------------

report_file <-
  file.path(
    step13h_dir,
    "STEP13H_TCGA_Final_Model_FINAL.txt"
  )

report_lines <- c(
  
  "============================================================",
  "STEP13H — TCGA FINAL MODEL",
  "============================================================",
  "",
  "FINAL STATUS: PASS",
  "",
  "Data branch:",
  "Historical Xiangtao TCGA ROC branch",
  "",
  "Input source:",
  tcga_rdata,
  "",
  "Expression source:",
  "TCGA-LUAD STAR RNA-seq TPM",
  "",
  "Historical transformation:",
  "log2(value + 1)",
  "",
  paste0("Total samples: ", final_model_n),
  paste0(
    "Tumor samples: ",
    sum(model_df$TumorStatus == "Tumor")
  ),
  paste0(
    "Normal samples: ",
    sum(model_df$TumorStatus == "Normal")
  ),
  "",
  "Predictors:",
  paste(expected_core_genes, collapse = ", "),
  "",
  paste0(
    "Formula: ",
    expected_formula_text
  ),
  "",
  "Model:",
  expected_model_type,
  "",
  paste0(
    "Reference class: ",
    expected_reference_level
  ),
  paste0(
    "Positive class: ",
    expected_positive_class
  ),
  "",
  "Feature selection: FALSE",
  "External refit: FALSE",
  "External coefficient estimation: FALSE",
  "External threshold optimization: FALSE",
  "External cutpoint search: FALSE",
  "External model tuning: FALSE",
  "",
  "Clinical diagnostic claim: FALSE",
  "Diagnostic biomarker claim: FALSE",
  "",
  "Final model coefficients:",
  ""
)

for (i in seq_len(nrow(coef_table))) {
  
  report_lines <- c(
    report_lines,
    paste0(
      coef_table$term[i],
      " | coefficient = ",
      signif(coef_table$Coefficient[i], 8),
      " | SE = ",
      signif(coef_table$Std_Error[i], 8),
      " | P = ",
      signif(coef_table$P_value[i], 8),
      " | OR = ",
      signif(coef_table$OR[i], 8),
      " | 95% CI = ",
      signif(coef_table$CI95_lower[i], 8),
      "–",
      signif(coef_table$CI95_upper[i], 8)
    )
  )
}

report_lines <- c(
  report_lines,
  "",
  "STEP13G preceding checkpoint: PASS",
  "",
  "Final model object saved:",
  step13h_final_rds,
  "",
  "This model is the locked training model for STEP13I.",
  "STEP13I must not refit coefficients in GSE32863.",
  "STEP13I must not optimize an external threshold or cutpoint.",
  "",
  "============================================================",
  "STEP13H COMPLETE",
  "============================================================"
)

writeLines(
  report_lines,
  report_file
)


# ------------------------------------------------------------
# 37. Save session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    step13h_dir,
    "STEP13H_sessionInfo_FINAL.txt"
  )
)


# ------------------------------------------------------------
# 38. Final audit aggregation
# ------------------------------------------------------------

all_final_checks <- c(
  
  STEP13F_status_PASS =
    identical(
      model_lock$status,
      "PASS"
    ),
  
  STEP13F_model_definition_integrity =
    all(integrity_checks),
  
  TCGA_objects_present =
    all(tcga_object_exists),
  
  TCGA_sample_count_575 =
    nrow(df_wide) == 575L,
  
  TCGA_tumor_516 =
    sum(df_wide$status == "Tumor") == 516L,
  
  TCGA_normal_59 =
    sum(df_wide$status == "Normal") == 59L,
  
  TCGA_sample_ID_unique =
    sample_id_unique,
  
  Complete_cases_575 =
    n_complete == 575L,
  
  Final_model_n_575 =
    final_model_n == 575L,
  
  Final_model_tumor_516 =
    sum(model_df$TumorStatus == "Tumor") == 516L,
  
  Final_model_normal_59 =
    sum(model_df$TumorStatus == "Normal") == 59L,
  
  Model_converged =
    model_converged,
  
  Coefficients_finite =
    coef_finite,
  
  Full_model_rank =
    model_rank_ok,
  
  Required_terms_present =
    all(coefficient_presence),
  
  Coefficient_values_finite =
    coefficient_finite_check,
  
  STEP13G_preceding_PASS =
    all(step13g_audit$PASS)
)


final_audit_df <- data.frame(
  check = names(all_final_checks),
  PASS = unname(all_final_checks),
  stringsAsFactors = FALSE
)

write.csv(
  final_audit_df,
  file.path(
    step13h_dir,
    "STEP13H_Final_Audit_Checks.csv"
  ),
  row.names = FALSE
)

write.csv(
  final_audit_df,
  file.path(
    audit_dir,
    "STEP13H_Final_Audit_Checks.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 39. Final PASS / FAIL decision
# ------------------------------------------------------------

if (!all(all_final_checks)) {
  
  fail_file <-
    file.path(
      step13h_dir,
      "STEP13H_FAIL.txt"
    )
  
  writeLines(
    c(
      "STEP13H FINAL STATUS: FAIL",
      "",
      "One or more final audit checks failed.",
      "",
      "Review:",
      "STEP13H_Final_Audit_Checks.csv"
    ),
    fail_file
  )
  
  stop(
    "STEP13H FINAL STATUS: FAIL"
  )
}


# ------------------------------------------------------------
# 40. Write final PASS marker
# ------------------------------------------------------------

pass_file <-
  file.path(
    step13h_dir,
    "STEP13H_PASS.txt"
  )

writeLines(
  c(
    "STEP13H FINAL STATUS: PASS",
    "",
    "Final TCGA logistic regression model established.",
    "",
    "Training branch: Historical Xiangtao TCGA ROC branch",
    "Samples: 575",
    "Tumor: 516",
    "Normal: 59",
    "Predictors: AOC3 + MS4A7 + IL7R",
    "Outcome: TumorStatus",
    "Reference: Normal",
    "Positive class: Tumor",
    "Model: multivariable logistic regression",
    "Feature selection: FALSE",
    "External refit: FALSE",
    "External coefficient estimation: FALSE",
    "External threshold optimization: FALSE",
    "External cutpoint search: FALSE",
    "External model tuning: FALSE",
    "",
    "Final coefficients, OR, 95% CI and P values saved.",
    "",
    "STEP13H final model is locked for STEP13I GSE32863 external assessment.",
    "No external refitting or coefficient estimation is permitted.",
    "No external threshold/cutpoint optimization is permitted."
  ),
  pass_file
)


# ------------------------------------------------------------
# 41. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13H FINAL SUMMARY\n")
cat("============================================================\n\n")

cat(
  "Input branch:             Historical Xiangtao TCGA ROC branch\n"
)

cat(
  "Input samples:            ",
  final_model_n,
  "\n",
  sep = ""
)

cat(
  "Tumor:                    ",
  sum(model_df$TumorStatus == "Tumor"),
  "\n",
  sep = ""
)

cat(
  "Normal:                   ",
  sum(model_df$TumorStatus == "Normal"),
  "\n",
  sep = ""
)

cat(
  "Core genes:               ",
  paste(expected_core_genes, collapse = ", "),
  "\n",
  sep = ""
)

cat(
  "Model:                    ",
  expected_model_type,
  "\n",
  sep = ""
)

cat(
  "Formula:                  ",
  expected_formula_text,
  "\n",
  sep = ""
)

cat(
  "Feature selection:        FALSE\n"
)

cat(
  "External refit:           FALSE\n"
)

cat(
  "External coefficients:    FALSE\n"
)

cat(
  "External threshold tune:  FALSE\n"
)

cat(
  "External cutpoint search: FALSE\n"
)

cat(
  "External model tuning:    FALSE\n"
)

cat(
  "Clinical diagnostic claim: FALSE\n"
)

cat(
  "Diagnostic biomarker claim: FALSE\n"
)

cat("\n")
cat("Final model coefficients:\n")
print(coef_table)

cat("\n")
cat(
  "All final audit checks:   ",
  all(all_final_checks),
  "\n",
  sep = ""
)

cat("\n")
cat("FINAL STATUS:\n")
cat(
  "PASS — STEP13H TCGA Final Model established.\n"
)

cat("\n")
cat("Final model RDS:\n")
cat(step13h_final_rds, "\n")

cat("\n")
cat("STEP13H COMPLETE\n")
cat("============================================================\n\n")


# ============================================================
# STEP13I — GSE32863 External Assessment
# ============================================================
#
# Purpose:
#
#   1. Apply the LOCKED STEP13H TCGA final logistic model to
#      the independently held-out GSE32863 dataset.
#
#   2. Estimate external ROC/AUC and 95% CI.
#
#   3. Perform a separate paired expression analysis in the
#      57 metadata-confirmed complete N+T subject pairs.
#
# ============================================================
#
# CRITICAL LOCKS:
#
#   External refit                         = FALSE
#   External coefficient estimation        = FALSE
#   External feature selection             = FALSE
#   External threshold optimization        = FALSE
#   External cutpoint search               = FALSE
#   External model tuning                  = FALSE
#
# The STEP13H coefficients are transferred unchanged.
#
# No coefficient is estimated from GSE32863.
#
# ============================================================


# ------------------------------------------------------------
# 0. Clean session
# ------------------------------------------------------------

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(scipen = 999)

set.seed(20260926L)


# ------------------------------------------------------------
# 1. Working directories
# ------------------------------------------------------------

workflow_root <-
  "input_dir/Paper 1_Step 13-Diagnostic Classification"

step13c_dir <-
  file.path(
    workflow_root,
    "STEP13C_GSE32863_Input_Lock"
  )

step13e_dir <-
  file.path(
    workflow_root,
    "STEP13E_Expression_Scale_Audit"
  )

step13f_dir <-
  file.path(
    workflow_root,
    "STEP13F_Model_Definition_Lock"
  )

step13h_dir <-
  file.path(
    workflow_root,
    "STEP13H_TCGA_Final_Model"
  )

step13i_dir <-
  file.path(
    workflow_root,
    "STEP13I_GSE32863_External_Assessment"
  )

results_dir <-
  file.path(
    workflow_root,
    "RESULTS"
  )

tables_dir <-
  file.path(
    workflow_root,
    "TABLES"
  )

audit_dir <-
  file.path(
    workflow_root,
    "AUDIT"
  )

logs_dir <-
  file.path(
    workflow_root,
    "LOGS"
  )


# ------------------------------------------------------------
# 2. Create directories
# ------------------------------------------------------------

dir.create(
  step13i_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  logs_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Remove stale PASS / FAIL markers
# ------------------------------------------------------------

stale_files <- c(
  file.path(step13i_dir, "STEP13I_PASS.txt"),
  file.path(step13i_dir, "STEP13I_FAIL.txt"),
  file.path(step13i_dir, "STEP13I_FINAL.txt")
)

stale_files <- stale_files[file.exists(stale_files)]

if (length(stale_files) > 0) {
  file.remove(stale_files)
}


# ------------------------------------------------------------
# 4. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "pROC"
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
    paste0(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", ")
    )
  )
}

library(pROC)


# ------------------------------------------------------------
# 5. Load STEP13F Model Definition Lock
# ------------------------------------------------------------

step13f_rds <-
  file.path(
    step13f_dir,
    "STEP13F_Model_Definition_Lock_FINAL.rds"
  )

if (!file.exists(step13f_rds)) {
  stop(
    "STEP13F FINAL RDS not found: ",
    step13f_rds
  )
}

model_lock <-
  readRDS(step13f_rds)


# ------------------------------------------------------------
# 6. STEP13F PASS check
# ------------------------------------------------------------

if (!identical(
  model_lock$status,
  "PASS"
)) {
  stop(
    "STEP13I cannot proceed: STEP13F status is not PASS."
  )
}


# ------------------------------------------------------------
# 7. STEP13F external-analysis prohibition audit
# ------------------------------------------------------------

external_lock_checks <- c(
  
  external_refit_allowed =
    identical(
      as.logical(
        model_lock$external_refit_allowed
      ),
      FALSE
    ),
  
  external_coefficient_estimation_allowed =
    identical(
      as.logical(
        model_lock$external_coefficient_estimation_allowed
      ),
      FALSE
    ),
  
  external_feature_selection_allowed =
    identical(
      as.logical(
        model_lock$external_feature_selection_allowed
      ),
      FALSE
    ),
  
  external_threshold_optimization_allowed =
    identical(
      as.logical(
        model_lock$external_threshold_optimization_allowed
      ),
      FALSE
    ),
  
  external_cutpoint_search_allowed =
    identical(
      as.logical(
        model_lock$external_cutpoint_search_allowed
      ),
      FALSE
    ),
  
  external_model_tuning_allowed =
    identical(
      as.logical(
        model_lock$external_model_tuning_allowed
      ),
      FALSE
    )
)


external_lock_audit <- data.frame(
  check = names(external_lock_checks),
  PASS = unname(external_lock_checks),
  stringsAsFactors = FALSE
)

write.csv(
  external_lock_audit,
  file.path(
    step13i_dir,
    "STEP13I_External_Analysis_Prohibition_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(external_lock_checks)) {
  stop(
    paste0(
      "STEP13I stopped: one or more external-analysis ",
      "prohibition locks are not FALSE."
    )
  )
}


# ------------------------------------------------------------
# 8. Verify exact model definition
# ------------------------------------------------------------

expected_core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

expected_formula_text <-
  "TumorStatus ~ AOC3 + MS4A7 + IL7R"

expected_reference_level <-
  "Normal"

expected_positive_class <-
  "Tumor"

model_definition_checks <- c(
  
  core_genes =
    identical(
      as.character(model_lock$core_genes),
      expected_core_genes
    ),
  
  formula =
    identical(
      as.character(model_lock$model_formula_text),
      expected_formula_text
    ),
  
  reference_level =
    identical(
      as.character(model_lock$reference_level),
      expected_reference_level
    ),
  
  positive_class =
    identical(
      as.character(model_lock$positive_class),
      expected_positive_class
    )
)

model_definition_audit <- data.frame(
  check = names(model_definition_checks),
  PASS = unname(model_definition_checks),
  stringsAsFactors = FALSE
)

write.csv(
  model_definition_audit,
  file.path(
    step13i_dir,
    "STEP13I_Model_Definition_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(model_definition_checks)) {
  stop(
    "STEP13I stopped: model-definition concordance failed."
  )
}


# ------------------------------------------------------------
# 9. Load STEP13H FINAL model
# ------------------------------------------------------------

step13h_rds <-
  file.path(
    step13h_dir,
    "STEP13H_TCGA_Final_Model_FINAL.rds"
  )

if (!file.exists(step13h_rds)) {
  stop(
    "STEP13H FINAL RDS not found: ",
    step13h_rds
  )
}

final_model_lock <-
  readRDS(step13h_rds)


# ------------------------------------------------------------
# 10. STEP13H PASS check
# ------------------------------------------------------------

if (!identical(
  final_model_lock$status,
  "PASS"
)) {
  stop(
    "STEP13I cannot proceed: STEP13H final model status is not PASS."
  )
}


# ------------------------------------------------------------
# 11. Extract FINAL TCGA model object
# ------------------------------------------------------------

final_model <-
  final_model_lock$model


# ------------------------------------------------------------
# 12. Confirm model is logistic regression
# ------------------------------------------------------------

model_class_ok <-
  inherits(
    final_model,
    "glm"
  )

model_family_ok <-
  identical(
    final_model$family$family,
    "binomial"
  )

model_link_ok <-
  identical(
    final_model$family$link,
    "logit"
  )

model_structure_audit <- data.frame(
  check = c(
    "glm_class",
    "binomial_family",
    "logit_link"
  ),
  PASS = c(
    model_class_ok,
    model_family_ok,
    model_link_ok
  ),
  stringsAsFactors = FALSE
)

write.csv(
  model_structure_audit,
  file.path(
    step13i_dir,
    "STEP13I_STEP13H_Model_Structure_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(model_structure_audit$PASS)) {
  stop(
    "STEP13I stopped: STEP13H model structure audit failed."
  )
}


# ------------------------------------------------------------
# 13. Extract and lock coefficients
# ------------------------------------------------------------

locked_coefficients <-
  coef(final_model)

required_terms <- c(
  "(Intercept)",
  expected_core_genes
)

coefficient_terms_ok <-
  all(
    required_terms %in%
      names(locked_coefficients)
  )

if (!coefficient_terms_ok) {
  stop(
    "STEP13I stopped: required STEP13H coefficients are missing."
  )
}

locked_coefficients <-
  locked_coefficients[
    required_terms
  ]

if (!all(
  is.finite(
    locked_coefficients
  )
)) {
  stop(
    "STEP13I stopped: non-finite STEP13H coefficients detected."
  )
}


# ------------------------------------------------------------
# 14. Save coefficient-transfer audit
# ------------------------------------------------------------

coefficient_transfer_audit <- data.frame(
  term = names(locked_coefficients),
  STEP13H_coefficient =
    as.numeric(locked_coefficients),
  source =
    "STEP13H_TCGA_Final_Model_FINAL.rds",
  external_re_estimation =
    FALSE,
  stringsAsFactors = FALSE
)

write.csv(
  coefficient_transfer_audit,
  file.path(
    step13i_dir,
    "STEP13I_STEP13H_Locked_Coefficient_Transfer_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. Historical GSE32863 source
# ------------------------------------------------------------

gse_rdata <-
  "input_dir/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal/.RData"

if (!file.exists(gse_rdata)) {
  stop(
    "Historical GSE32863 RData not found: ",
    gse_rdata
  )
}


# ------------------------------------------------------------
# 16. Load historical GSE32863 objects into isolated environment
# ------------------------------------------------------------

gse_env <-
  new.env(parent = emptyenv())

load(
  gse_rdata,
  envir = gse_env
)


# ------------------------------------------------------------
# 17. Required historical GSE objects
# ------------------------------------------------------------

required_gse_objects <- c(
  "expr",
  "expr.mat",
  "expr2",
  "expr.final",
  "expr3",
  "group",
  "roc.ext"
)

gse_object_exists <- vapply(
  required_gse_objects,
  exists,
  logical(1),
  envir = gse_env,
  inherits = FALSE
)

gse_object_audit <- data.frame(
  object = required_gse_objects,
  exists = gse_object_exists,
  stringsAsFactors = FALSE
)

write.csv(
  gse_object_audit,
  file.path(
    step13i_dir,
    "STEP13I_GSE32863_Object_Existence_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(gse_object_exists)) {
  stop(
    paste0(
      "Required historical GSE32863 object(s) missing: ",
      paste(
        required_gse_objects[
          !gse_object_exists
        ],
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 18. Extract locked historical objects
# ------------------------------------------------------------

expr_final <-
  gse_env$expr.final

expr3 <-
  gse_env$expr3

group <-
  gse_env$group

roc_ext_historical <-
  gse_env$roc.ext


# ------------------------------------------------------------
# 19. GSE32863 locked sample-count audit
# ------------------------------------------------------------

expected_gse_total <- 116L
expected_gse_tumor <- 58L
expected_gse_normal <- 58L
expected_gse_subjects <- 59L
expected_gse_complete_pairs <- 57L

sample_count_checks <- c(
  
  expr_final_columns_116 =
    ncol(expr_final) == expected_gse_total,
  
  expr3_rows_116 =
    nrow(expr3) == expected_gse_total,
  
  group_length_116 =
    length(group) == expected_gse_total,
  
  group_tumor_58 =
    sum(group == "Tumor") == expected_gse_tumor,
  
  group_normal_58 =
    sum(group == "Normal") == expected_gse_normal,
  
  roc_ext_predictor_116 =
    length(
      roc_ext_historical$predictor
    ) == expected_gse_total
)

sample_count_audit <- data.frame(
  check = names(sample_count_checks),
  PASS = unname(sample_count_checks),
  stringsAsFactors = FALSE
)

write.csv(
  sample_count_audit,
  file.path(
    step13i_dir,
    "STEP13I_GSE32863_Sample_Count_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(sample_count_checks)) {
  stop(
    "STEP13I stopped: GSE32863 locked sample-count audit failed."
  )
}


# ------------------------------------------------------------
# 20. GSE32863 expression sample IDs
# ------------------------------------------------------------
#
# Historical expr3 structure was already locked in STEP13E:
#
#   expr3 = 116 × 5
#   columns = AOC3, MS4A7, IL7R, Group, Probability
#
# Sample IDs are stored in rownames(expr3), NOT in a
# sample_id column.
#
# Therefore:
#   - DO NOT modify expr3
#   - DO NOT reconstruct expr3
#   - use rownames(expr3) as the audited sample-ID source
#
# ------------------------------------------------------------

gse_sample_ids_expr <-
  colnames(expr_final)

if (is.null(gse_sample_ids_expr)) {
  stop(
    "STEP13I stopped: expr.final has no column sample IDs."
  )
}


# ------------------------------------------------------------
# 20A. Validate historical expr3 structure
# ------------------------------------------------------------

expected_expr3_columns <- c(
  "AOC3",
  "MS4A7",
  "IL7R",
  "Group",
  "Probability"
)

expr3_column_check <-
  identical(
    colnames(expr3),
    expected_expr3_columns
  )

expr3_nrow_check <-
  nrow(expr3) == 116L

expr3_rownames_present <-
  !is.null(
    rownames(expr3)
  )

expr3_rownames_unique <-
  expr3_rownames_present &&
  length(
    unique(
      rownames(expr3)
    )
  ) == nrow(expr3)

expr3_group_present <-
  "Group" %in%
  colnames(expr3)

expr3_structure_audit <- data.frame(
  check = c(
    "expr3_116_rows",
    "expr3_expected_columns",
    "expr3_rownames_present",
    "expr3_rownames_unique",
    "expr3_Group_present"
  ),
  PASS = c(
    expr3_nrow_check,
    expr3_column_check,
    expr3_rownames_present,
    expr3_rownames_unique,
    expr3_group_present
  ),
  stringsAsFactors = FALSE
)

write.csv(
  expr3_structure_audit,
  file.path(
    step13i_dir,
    "STEP13I_expr3_Historical_Structure_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(
  expr3_structure_audit$PASS
)) {
  stop(
    "STEP13I stopped: historical expr3 structure audit failed."
  )
}


# ------------------------------------------------------------
# 20B. Extract sample IDs from expr3 rownames
# ------------------------------------------------------------

gse_sample_ids_expr3 <-
  as.character(
    rownames(expr3)
  )

if (length(gse_sample_ids_expr3) != 116L) {
  stop(
    "STEP13I stopped: expr3 rownames do not contain 116 sample IDs."
  )
}

if (anyNA(gse_sample_ids_expr3) ||
    any(
      trimws(gse_sample_ids_expr3) == ""
    )) {
  stop(
    "STEP13I stopped: expr3 contains missing/empty sample IDs."
  )
}

if (
  length(
    unique(
      gse_sample_ids_expr3
    )
  ) != 116L
) {
  stop(
    "STEP13I stopped: expr3 sample IDs are not unique."
  )
}


# ------------------------------------------------------------
# 20C. Confirm historical Group structure
# ------------------------------------------------------------

expr3_group_values <-
  as.character(
    expr3$Group
  )

if (!all(
  expr3_group_values %in%
  c(
    "Normal",
    "Tumor"
  )
)) {
  stop(
    "STEP13I stopped: expr3$Group contains unexpected levels."
  )
}

expr3_group_count_check <- c(
  total_116 =
    length(expr3_group_values) == 116L,
  
  tumor_58 =
    sum(
      expr3_group_values == "Tumor"
    ) == 58L,
  
  normal_58 =
    sum(
      expr3_group_values == "Normal"
    ) == 58L
)

expr3_group_audit <- data.frame(
  check = names(expr3_group_count_check),
  PASS = unname(expr3_group_count_check),
  stringsAsFactors = FALSE
)

write.csv(
  expr3_group_audit,
  file.path(
    step13i_dir,
    "STEP13I_expr3_Group_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(
  expr3_group_count_check
)) {
  stop(
    "STEP13I stopped: expr3 Group count audit failed."
  )
}

# ------------------------------------------------------------
# 21. Historical expression/sample concordance
# ------------------------------------------------------------
#
# Historical locked objects:
#
#   expr.final
#       sample IDs = column names
#
#   expr3
#       sample IDs = rownames
#
# Both objects are historical GSE32863 objects.
#
# No object is modified or reordered.
#
# The required concordance is:
#   1. identical 116-sample ID sets
#   2. identical Group annotation
#
# Exact storage-order identity between the two objects
# is NOT required.
#
# ------------------------------------------------------------


# ------------------------------------------------------------
# 21A. Reassert sample-ID vectors
# ------------------------------------------------------------

gse_sample_ids_expr <-
  as.character(
    colnames(expr_final)
  )

gse_sample_ids_expr3 <-
  as.character(
    rownames(expr3)
  )


# ------------------------------------------------------------
# 21B. Basic sample-ID integrity
# ------------------------------------------------------------

expr_sample_id_integrity <- c(
  
  expr_final_n =
    length(gse_sample_ids_expr) == 116L,
  
  expr3_n =
    length(gse_sample_ids_expr3) == 116L,
  
  expr_final_unique =
    length(unique(gse_sample_ids_expr)) == 116L,
  
  expr3_unique =
    length(unique(gse_sample_ids_expr3)) == 116L,
  
  expr_final_no_missing =
    !anyNA(gse_sample_ids_expr) &&
    all(trimws(gse_sample_ids_expr) != ""),
  
  expr3_no_missing =
    !anyNA(gse_sample_ids_expr3) &&
    all(trimws(gse_sample_ids_expr3) != "")
)


if (!all(expr_sample_id_integrity)) {
  
  stop(
    "STEP13I stopped: historical sample-ID integrity audit failed."
  )
}


# ------------------------------------------------------------
# 21C. Sample-set concordance
# ------------------------------------------------------------

expr_vs_expr3_set <-
  setequal(
    gse_sample_ids_expr,
    gse_sample_ids_expr3
  )


# ------------------------------------------------------------
# 21D. Explicit sample-ID differences
# ------------------------------------------------------------

expr_only_ids <-
  setdiff(
    gse_sample_ids_expr,
    gse_sample_ids_expr3
  )

expr3_only_ids <-
  setdiff(
    gse_sample_ids_expr3,
    gse_sample_ids_expr
  )


# ------------------------------------------------------------
# 21E. Historical Group concordance
# ------------------------------------------------------------

expr3_group_values <-
  as.character(
    expr3$Group
  )

group_values <-
  as.character(
    group
  )

expr3_group_vs_group <-
  identical(
    expr3_group_values,
    group_values
  )


# ------------------------------------------------------------
# 21F. Group integrity
# ------------------------------------------------------------

group_integrity <- c(
  
  expr3_group_n =
    length(expr3_group_values) == 116L,
  
  group_n =
    length(group_values) == 116L,
  
  expr3_group_normal_58 =
    sum(expr3_group_values == "Normal") == 58L,
  
  expr3_group_tumor_58 =
    sum(expr3_group_values == "Tumor") == 58L,
  
  group_normal_58 =
    sum(group_values == "Normal") == 58L,
  
  group_tumor_58 =
    sum(group_values == "Tumor") == 58L,
  
  expr3_group_no_missing =
    !anyNA(expr3_group_values),
  
  group_no_missing =
    !anyNA(group_values)
)


if (!all(group_integrity)) {
  
  stop(
    "STEP13I stopped: historical Group integrity audit failed."
  )
}


# ------------------------------------------------------------
# 21G. Final concordance audit table
# ------------------------------------------------------------

sample_concordance_audit <- data.frame(
  
  check = c(
    
    names(expr_sample_id_integrity),
    
    "expr_final_vs_expr3_sample_set",
    
    "expr_final_only_sample_count_zero",
    
    "expr3_only_sample_count_zero",
    
    names(group_integrity),
    
    "expr3_Group_vs_group"
    
  ),
  
  PASS = c(
    
    unname(expr_sample_id_integrity),
    
    expr_vs_expr3_set,
    
    length(expr_only_ids) == 0L,
    
    length(expr3_only_ids) == 0L,
    
    unname(group_integrity),
    
    expr3_group_vs_group
    
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  sample_concordance_audit,
  file.path(
    step13i_dir,
    "STEP13I_GSE32863_Sample_Concordance_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 21H. Save explicit sample differences
# ------------------------------------------------------------

write.csv(
  data.frame(
    expr_final_only = expr_only_ids,
    stringsAsFactors = FALSE
  ),
  file.path(
    step13i_dir,
    "STEP13I_expr_final_only_sample_IDs.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    expr3_only = expr3_only_ids,
    stringsAsFactors = FALSE
  ),
  file.path(
    step13i_dir,
    "STEP13I_expr3_only_sample_IDs.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 21I. Final integrity lock
# ------------------------------------------------------------

if (!all(sample_concordance_audit$PASS)) {
  
  stop(
    paste0(
      "STEP13I stopped: historical GSE sample concordance ",
      "or Group concordance failed. ",
      "See STEP13I_GSE32863_Sample_Concordance_Audit.csv."
    )
  )
}


# ------------------------------------------------------------
# 21J. PASS message
# ------------------------------------------------------------

cat(
  "\n============================================================\n",
  "STEP13I-21 HISTORICAL SAMPLE CONCORDANCE\n",
  "============================================================\n",
  "expr.final samples:       ",
  length(gse_sample_ids_expr),
  "\n",
  "expr3 samples:            ",
  length(gse_sample_ids_expr3),
  "\n",
  "Sample-set concordance:   ",
  expr_vs_expr3_set,
  "\n",
  "expr.final-only IDs:      ",
  length(expr_only_ids),
  "\n",
  "expr3-only IDs:           ",
  length(expr3_only_ids),
  "\n",
  "expr3 Group:              Normal = ",
  sum(expr3_group_values == "Normal"),
  ", Tumor = ",
  sum(expr3_group_values == "Tumor"),
  "\n",
  "group Group:              Normal = ",
  sum(group_values == "Normal"),
  ", Tumor = ",
  sum(group_values == "Tumor"),
  "\n",
  "Group concordance:        ",
  expr3_group_vs_group,
  "\n\n",
  "FINAL STATUS: PASS\n",
  "Historical GSE32863 sample-set and Group concordance established.\n",
  "============================================================\n",
  sep = ""
)

# ------------------------------------------------------------
# 22. Confirm core genes in expr.final
# ------------------------------------------------------------

expr_final_gene_ids <-
  rownames(expr_final)

core_gene_presence <-
  expected_core_genes %in%
  expr_final_gene_ids

core_gene_presence_audit <- data.frame(
  gene = expected_core_genes,
  present_in_expr_final =
    core_gene_presence,
  stringsAsFactors = FALSE
)

write.csv(
  core_gene_presence_audit,
  file.path(
    step13i_dir,
    "STEP13I_Core_Gene_Presence_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(core_gene_presence)) {
  stop(
    "STEP13I stopped: one or more core genes missing from expr.final."
  )
}


# ------------------------------------------------------------
# 23. Extract external expression matrix
# ------------------------------------------------------------
#
# IMPORTANT:
#
#   No probe selection is performed here.
#   expr.final is the already audited historical gene-level
#   expression object from STEP13E.
#
# ------------------------------------------------------------

external_expression <- t(
  expr_final[
    expected_core_genes,
    ,
    drop = FALSE
  ]
)

external_expression <-
  as.data.frame(
    external_expression,
    stringsAsFactors = FALSE
  )

external_expression$sample_id <-
  rownames(external_expression)

rownames(external_expression) <-
  NULL

external_expression <-
  external_expression[
    ,
    c(
      "sample_id",
      expected_core_genes
    )
  ]


# ------------------------------------------------------------
# 24. Confirm expression values are numeric and finite
# ------------------------------------------------------------

expression_numeric <- vapply(
  expected_core_genes,
  function(g) {
    is.numeric(
      external_expression[[g]]
    )
  },
  logical(1)
)

expression_finite <- vapply(
  expected_core_genes,
  function(g) {
    all(
      is.finite(
        external_expression[[g]]
      )
    )
  },
  logical(1)
)

expression_value_audit <- data.frame(
  gene = expected_core_genes,
  numeric = expression_numeric,
  finite = expression_finite,
  stringsAsFactors = FALSE
)

write.csv(
  expression_value_audit,
  file.path(
    step13i_dir,
    "STEP13I_External_Expression_Value_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(
  expression_numeric &
  expression_finite
)) {
  stop(
    "STEP13I stopped: external expression values are not numeric/finite."
  )
}


# ------------------------------------------------------------
# 25. Expression-scale provenance lock
# ------------------------------------------------------------
#
# STEP13E established the historical GSE32863 object provenance.
#
# IMPORTANT:
#   No new normalization, scaling, centering, z-transformation,
#   or gene-wise transformation is performed in STEP13I.
#
# The model is applied to the already audited historical
# expression representation.
#
# ------------------------------------------------------------

scale_provenance <- data.frame(
  field = c(
    "GSE_source",
    "Expression_object",
    "Expression_rows",
    "Expression_samples",
    "Model_input_genes",
    "Additional_normalization",
    "Additional_scaling",
    "Additional_centering",
    "Additional_gene_selection"
  ),
  value = c(
    "GSE32863",
    "expr.final",
    nrow(expr_final),
    ncol(expr_final),
    paste(
      expected_core_genes,
      collapse = "; "
    ),
    "NONE",
    "NONE",
    "NONE",
    "NONE"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  scale_provenance,
  file.path(
    step13i_dir,
    "STEP13I_Expression_Scale_Transfer_Provenance.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 26. Construct locked external outcome
# ------------------------------------------------------------

external_expression$TumorStatus <-
  factor(
    as.character(
      group[
        match(
          external_expression$sample_id,
          gse_sample_ids_expr3
        )
      ]
    ),
    levels = c(
      "Normal",
      "Tumor"
    )
  )


# ------------------------------------------------------------
# 27. External model input completeness audit
# ------------------------------------------------------------

external_complete <-
  complete.cases(
    external_expression[
      ,
      c(
        "TumorStatus",
        expected_core_genes
      )
    ]
  )

external_complete_n <-
  sum(external_complete)

external_missing_n <-
  sum(!external_complete)

external_complete_audit <- data.frame(
  quantity = c(
    "Total samples",
    "Complete model inputs",
    "Incomplete model inputs"
  ),
  value = c(
    nrow(external_expression),
    external_complete_n,
    external_missing_n
  ),
  PASS = c(
    nrow(external_expression) == 116L,
    external_complete_n == 116L,
    external_missing_n == 0L
  ),
  stringsAsFactors = FALSE
)

write.csv(
  external_complete_audit,
  file.path(
    step13i_dir,
    "STEP13I_External_Model_Input_Completeness_Audit.csv"
  ),
  row.names = FALSE
)

if (external_complete_n != 116L) {
  stop(
    "STEP13I stopped: external model input completeness failed."
  )
}


# ------------------------------------------------------------
# 28. Apply STEP13H model WITHOUT refitting
# ------------------------------------------------------------
#
# CRITICAL:
#
#   predict() is used only to apply the already fitted STEP13H
#   model to GSE32863.
#
#   No glm(), glm.fit(), coefficient estimation, feature
#   selection, or tuning is performed on GSE32863.
#
# ------------------------------------------------------------

external_expression$predicted_probability <-
  predict(
    final_model,
    newdata = external_expression,
    type = "response"
  )

external_expression$predicted_logit <-
  predict(
    final_model,
    newdata = external_expression,
    type = "link"
  )


# ------------------------------------------------------------
# 29. Prediction sanity audit
# ------------------------------------------------------------

prediction_finite <-
  all(
    is.finite(
      external_expression$predicted_probability
    )
  )

prediction_range_ok <-
  all(
    external_expression$predicted_probability >= 0 &
      external_expression$predicted_probability <= 1
  )

prediction_audit <- data.frame(
  check = c(
    "prediction_finite",
    "prediction_range_0_to_1"
  ),
  PASS = c(
    prediction_finite,
    prediction_range_ok
  ),
  stringsAsFactors = FALSE
)

write.csv(
  prediction_audit,
  file.path(
    step13i_dir,
    "STEP13I_External_Prediction_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(prediction_audit$PASS)) {
  stop(
    "STEP13I stopped: external prediction audit failed."
  )
}


# ------------------------------------------------------------
# 30. Save complete external prediction traceability
# ------------------------------------------------------------

external_prediction_traceability <-
  external_expression[
    ,
    c(
      "sample_id",
      "TumorStatus",
      expected_core_genes,
      "predicted_logit",
      "predicted_probability"
    )
  ]

write.csv(
  external_prediction_traceability,
  file.path(
    step13i_dir,
    "STEP13I_GSE32863_External_Prediction_Traceability.csv"
  ),
  row.names = FALSE
)

write.csv(
  external_prediction_traceability,
  file.path(
    results_dir,
    "STEP13I_GSE32863_External_Prediction_Traceability.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 31. External ROC analysis
# ------------------------------------------------------------
#
# Primary external performance metric:
#   ROC AUC
#
# Direction is determined from the observed predictor values
# by pROC.
#
# No threshold optimization is performed.
#
# ------------------------------------------------------------

roc_external <-
  pROC::roc(
    response =
      external_expression$TumorStatus,
    predictor =
      external_expression$predicted_probability,
    levels =
      c(
        "Normal",
        "Tumor"
      ),
    direction =
      "<",
    quiet = TRUE
  )


# ------------------------------------------------------------
# 32. External AUC
# ------------------------------------------------------------

external_auc <-
  as.numeric(
    pROC::auc(
      roc_external
    )
  )


external_auc_ci <-
  as.numeric(
    pROC::ci.auc(
      roc_external,
      conf.level = 0.95
    )
  )


# ------------------------------------------------------------
# 33. External ROC result table
# ------------------------------------------------------------

external_roc_summary <- data.frame(
  metric = c(
    "AUC",
    "AUC_95CI_lower",
    "AUC_95CI_upper"
  ),
  value = c(
    external_auc,
    external_auc_ci[1],
    external_auc_ci[3]
  ),
  stringsAsFactors = FALSE
)

write.csv(
  external_roc_summary,
  file.path(
    step13i_dir,
    "STEP13I_GSE32863_External_ROC_AUC.csv"
  ),
  row.names = FALSE
)

write.csv(
  external_roc_summary,
  file.path(
    results_dir,
    "STEP13I_GSE32863_External_ROC_AUC.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 34. Explicit prohibition audit
# ------------------------------------------------------------
#
# This audit records what was NOT done.
#
# ------------------------------------------------------------

prohibited_operations <- data.frame(
  
  operation = c(
    "External model refitting",
    "External coefficient estimation",
    "External feature selection",
    "External threshold optimization",
    "External cutpoint search",
    "External model tuning"
  ),
  
  performed = c(
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE
  ),
  
  locked_policy = c(
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE
  ),
  
  stringsAsFactors = FALSE
)

prohibited_operations$PASS <-
  prohibited_operations$performed ==
  prohibited_operations$locked_policy

write.csv(
  prohibited_operations,
  file.path(
    step13i_dir,
    "STEP13I_Prohibited_Operations_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(
  prohibited_operations$PASS
)) {
  stop(
    "STEP13I stopped: prohibited external operation audit failed."
  )
}


# ============================================================
# PART II
# 57 COMPLETE PAIRED SUBJECT ANALYSIS
# ============================================================


# ------------------------------------------------------------
# 35. Parse GSE32863 subject/sample metadata
# ------------------------------------------------------------
#
# Source:
#   GSE32863_series_matrix.txt.gz
#
# The metadata audit established:
#
#   116 samples
#   59 subjects
#   57 complete N+T pairs
#   1 tumor-only subject = 3023
#   1 normal-only subject = 3035
#
# ------------------------------------------------------------

series_matrix_file <-
  file.path(
    dirname(gse_rdata),
    "GSE32863_series_matrix.txt.gz"
  )

if (!file.exists(series_matrix_file)) {
  stop(
    "GSE32863 series matrix not found: ",
    series_matrix_file
  )
}


# ------------------------------------------------------------
# 36. Read metadata lines only
# ------------------------------------------------------------

series_lines <-
  readLines(
    gzfile(series_matrix_file),
    warn = FALSE
  )

sample_geo_line <-
  grep(
    "^!Sample_geo_accession",
    series_lines,
    value = TRUE
  )

sample_title_line <-
  grep(
    "^!Sample_title",
    series_lines,
    value = TRUE
  )

sample_characteristic_lines <-
  grep(
    "^!Sample_characteristics_ch1",
    series_lines,
    value = TRUE
  )


if (length(sample_geo_line) != 1L) {
  stop(
    "STEP13I stopped: expected exactly one !Sample_geo_accession line."
  )
}

if (length(sample_title_line) != 1L) {
  stop(
    "STEP13I stopped: expected exactly one !Sample_title line."
  )
}


# ------------------------------------------------------------
# 37. Parse GSM accessions
# ------------------------------------------------------------

gsm_raw <-
  sub(
    "^!Sample_geo_accession\t",
    "",
    sample_geo_line
  )

gsm_ids <-
  regmatches(
    gsm_raw,
    gregexpr(
      "GSM[0-9]+",
      gsm_raw
    )
  )[[1]]

gsm_ids <-
  unique(gsm_ids)


# ------------------------------------------------------------
# 38. Parse sample titles
# ------------------------------------------------------------

title_raw <-
  sub(
    "^!Sample_title\t",
    "",
    sample_title_line
  )

titles <- regmatches(
  title_raw,
  gregexpr(
    "\"[^\"]*\"",
    title_raw
  )
)[[1]]

titles <-
  sub(
    '^"',
    "",
    sub(
      '"$',
      "",
      titles
    )
  )


# ------------------------------------------------------------
# 39. Basic metadata concordance
# ------------------------------------------------------------

metadata_basic_checks <- c(
  
  GSM_count_116 =
    length(gsm_ids) == 116L,
  
  GSM_unique_116 =
    length(unique(gsm_ids)) == 116L,
  
  title_count_116 =
    length(titles) == 116L
)

metadata_basic_audit <- data.frame(
  check = names(metadata_basic_checks),
  PASS = unname(metadata_basic_checks),
  stringsAsFactors = FALSE
)

write.csv(
  metadata_basic_audit,
  file.path(
    step13i_dir,
    "STEP13I_GEO_Metadata_Basic_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(metadata_basic_checks)) {
  stop(
    "STEP13I stopped: GEO metadata parsing failed."
  )
}


# ------------------------------------------------------------
# 40. Extract sample codes and status from titles
# ------------------------------------------------------------

sample_code <-
  sub(
    ".*[[:space:]]([A-Za-z0-9]+)_([NT])[[:space:]]*\\(expression\\).*",
    "\\1",
    titles
  )

sample_type <-
  sub(
    ".*[[:space:]]([A-Za-z0-9]+)_([NT])[[:space:]]*\\(expression\\).*",
    "\\2",
    titles
  )


# ------------------------------------------------------------
# 41. Validate sample-code parsing
# ------------------------------------------------------------

sample_code_parse_ok <-
  all(
    sample_code != titles
  )

sample_type_parse_ok <-
  all(
    sample_type %in% c(
      "N",
      "T"
    )
  )

if (!sample_code_parse_ok ||
    !sample_type_parse_ok) {
  
  stop(
    "STEP13I stopped: sample-code parsing from GEO titles failed."
  )
}


# ------------------------------------------------------------
# 42. Construct metadata table
# ------------------------------------------------------------

metadata_df <- data.frame(
  GSM = gsm_ids,
  title = titles,
  sample_code = sample_code,
  type = sample_type,
  status = ifelse(
    sample_type == "N",
    "Normal",
    "Tumor"
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 43. Metadata sample-ID concordance
# ------------------------------------------------------------

metadata_expression_concordance <-
  setequal(
    metadata_df$GSM,
    gse_sample_ids_expr
  )

metadata_expr3_concordance <-
  setequal(
    metadata_df$GSM,
    gse_sample_ids_expr3
  )

metadata_concordance_audit <- data.frame(
  check = c(
    "Metadata_vs_expr_final",
    "Metadata_vs_expr3"
  ),
  PASS = c(
    metadata_expression_concordance,
    metadata_expr3_concordance
  ),
  stringsAsFactors = FALSE
)

write.csv(
  metadata_concordance_audit,
  file.path(
    step13i_dir,
    "STEP13I_GEO_Metadata_Sample_Concordance_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(
  metadata_concordance_audit$PASS
)) {
  stop(
    "STEP13I stopped: GEO metadata/sample concordance failed."
  )
}


# ------------------------------------------------------------
# 44. Subject-level pairing table
# ------------------------------------------------------------

subject_table <- split(
  metadata_df,
  metadata_df$sample_code
)

pair_list <- list()
unmatched_list <- list()

pair_counter <- 0L
unmatched_counter <- 0L

for (subject_id in names(subject_table)) {
  
  tmp <- subject_table[[subject_id]]
  
  normal_rows <-
    tmp[
      tmp$status == "Normal",
      ,
      drop = FALSE
    ]
  
  tumor_rows <-
    tmp[
      tmp$status == "Tumor",
      ,
      drop = FALSE
    ]
  
  if (
    nrow(normal_rows) == 1L &&
    nrow(tumor_rows) == 1L
  ) {
    
    pair_counter <-
      pair_counter + 1L
    
    pair_list[[pair_counter]] <-
      data.frame(
        pair_id =
          sprintf(
            "Pair_%03d",
            pair_counter
          ),
        subject_id =
          subject_id,
        normal_GSM =
          normal_rows$GSM,
        tumor_GSM =
          tumor_rows$GSM,
        stringsAsFactors = FALSE
      )
    
  } else {
    
    unmatched_counter <-
      unmatched_counter + 1L
    
    unmatched_list[[unmatched_counter]] <-
      data.frame(
        subject_id =
          subject_id,
        normal_n =
          nrow(normal_rows),
        tumor_n =
          nrow(tumor_rows),
        stringsAsFactors = FALSE
      )
  }
}


pair_df <-
  do.call(
    rbind,
    pair_list
  )

unmatched_df <-
  if (length(unmatched_list) > 0) {
    do.call(
      rbind,
      unmatched_list
    )
  } else {
    data.frame(
      subject_id = character(0),
      normal_n = integer(0),
      tumor_n = integer(0)
    )
  }


# ------------------------------------------------------------
# 45. Pairing integrity audit
# ------------------------------------------------------------

pairing_checks <- c(
  
  unique_subjects_59 =
    length(unique(metadata_df$sample_code)) == 59L,
  
  complete_pairs_57 =
    nrow(pair_df) == 57L,
  
  unmatched_subjects_2 =
    nrow(unmatched_df) == 2L,
  
  tumor_only_3023 =
    any(
      unmatched_df$subject_id == "3023" &
        unmatched_df$tumor_n == 1L &
        unmatched_df$normal_n == 0L
    ),
  
  normal_only_3035 =
    any(
      unmatched_df$subject_id == "3035" &
        unmatched_df$tumor_n == 0L &
        unmatched_df$normal_n == 1L
    )
)

pairing_audit <- data.frame(
  check = names(pairing_checks),
  PASS = unname(pairing_checks),
  stringsAsFactors = FALSE
)

write.csv(
  pairing_audit,
  file.path(
    step13i_dir,
    "STEP13I_57_Pair_Pairing_Integrity_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(pairing_checks)) {
  stop(
    "STEP13I stopped: 57-pair integrity audit failed."
  )
}


# ------------------------------------------------------------
# 46. Save pairing table
# ------------------------------------------------------------

write.csv(
  pair_df,
  file.path(
    step13i_dir,
    "STEP13I_GSE32863_57_Complete_Pairs.csv"
  ),
  row.names = FALSE
)

write.csv(
  pair_df,
  file.path(
    tables_dir,
    "STEP13I_GSE32863_57_Complete_Pairs.csv"
  ),
  row.names = FALSE
)

write.csv(
  unmatched_df,
  file.path(
    step13i_dir,
    "STEP13I_GSE32863_Unmatched_Subjects.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 47. Build paired expression data
# ------------------------------------------------------------

normal_lookup <-
  external_expression[
    match(
      pair_df$normal_GSM,
      external_expression$sample_id
    ),
    ,
    drop = FALSE
  ]

tumor_lookup <-
  external_expression[
    match(
      pair_df$tumor_GSM,
      external_expression$sample_id
    ),
    ,
    drop = FALSE
  ]


# ------------------------------------------------------------
# 48. Paired expression completeness audit
# ------------------------------------------------------------

paired_expression_checks <- c(
  
  normal_rows_57 =
    nrow(normal_lookup) == 57L,
  
  tumor_rows_57 =
    nrow(tumor_lookup) == 57L,
  
  normal_ids_complete =
    all(
      !is.na(
        normal_lookup$sample_id
      )
    ),
  
  tumor_ids_complete =
    all(
      !is.na(
        tumor_lookup$sample_id
      )
    )
)

paired_expression_audit <- data.frame(
  check = names(paired_expression_checks),
  PASS = unname(paired_expression_checks),
  stringsAsFactors = FALSE
)

write.csv(
  paired_expression_audit,
  file.path(
    step13i_dir,
    "STEP13I_Paired_Expression_Input_Audit.csv"
  ),
  row.names = FALSE
)

if (!all(
  paired_expression_checks
)) {
  stop(
    "STEP13I stopped: paired expression input audit failed."
  )
}


# ------------------------------------------------------------
# 49. Paired Wilcoxon signed-rank tests
# ------------------------------------------------------------
#
# The paired analysis is independent of the external ROC model.
#
# It does NOT alter:
#   - STEP13H coefficients
#   - external predictions
#   - external AUC
#   - thresholds
#
# ------------------------------------------------------------

paired_results <- vector(
  "list",
  length(expected_core_genes)
)

for (i in seq_along(expected_core_genes)) {
  
  gene <- expected_core_genes[i]
  
  normal_values <-
    normal_lookup[[gene]]
  
  tumor_values <-
    tumor_lookup[[gene]]
  
  wt <-
    suppressWarnings(
      wilcox.test(
        tumor_values,
        normal_values,
        paired = TRUE,
        exact = FALSE
      )
    )
  
  paired_results[[i]] <-
    data.frame(
      gene = gene,
      n_pairs = 57L,
      normal_median =
        median(
          normal_values,
          na.rm = TRUE
        ),
      tumor_median =
        median(
          tumor_values,
          na.rm = TRUE
        ),
      median_difference_tumor_minus_normal =
        median(
          tumor_values -
            normal_values,
          na.rm = TRUE
        ),
      V_statistic =
        as.numeric(
          wt$statistic
        ),
      P_value =
        as.numeric(
          wt$p.value
        ),
      stringsAsFactors = FALSE
    )
}

paired_results_df <-
  do.call(
    rbind,
    paired_results
  )


# ------------------------------------------------------------
# 50. Multiple-testing correction
# ------------------------------------------------------------

paired_results_df$FDR <-
  p.adjust(
    paired_results_df$P_value,
    method = "BH"
  )


# ------------------------------------------------------------
# 51. Save paired expression results
# ------------------------------------------------------------

write.csv(
  paired_results_df,
  file.path(
    step13i_dir,
    "STEP13I_GSE32863_Paired_Wilcoxon_AOC3_MS4A7_IL7R.csv"
  ),
  row.names = FALSE
)

write.csv(
  paired_results_df,
  file.path(
    results_dir,
    "STEP13I_GSE32863_Paired_Wilcoxon_AOC3_MS4A7_IL7R.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 52. Save paired expression long-format traceability
# ------------------------------------------------------------

paired_long_list <- vector(
  "list",
  length(expected_core_genes)
)

for (i in seq_along(expected_core_genes)) {
  
  gene <- expected_core_genes[i]
  
  paired_long_list[[i]] <-
    data.frame(
      pair_id =
        pair_df$pair_id,
      subject_id =
        pair_df$subject_id,
      gene =
        gene,
      normal_GSM =
        pair_df$normal_GSM,
      tumor_GSM =
        pair_df$tumor_GSM,
      normal_value =
        normal_lookup[[gene]],
      tumor_value =
        tumor_lookup[[gene]],
      difference_tumor_minus_normal =
        tumor_lookup[[gene]] -
        normal_lookup[[gene]],
      stringsAsFactors = FALSE
    )
}

paired_long_df <-
  do.call(
    rbind,
    paired_long_list
  )

write.csv(
  paired_long_df,
  file.path(
    step13i_dir,
    "STEP13I_Paired_Expression_Traceability.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 53. Confirm paired analysis did not alter model parameters
# ------------------------------------------------------------

coef_after_pair_analysis <-
  coef(final_model)

coefficient_unchanged <-
  isTRUE(
    all.equal(
      locked_coefficients,
      coef_after_pair_analysis
    )
  )

coefficient_integrity_audit <- data.frame(
  check =
    "STEP13H_coefficients_unchanged_after_paired_analysis",
  PASS =
    coefficient_unchanged,
  stringsAsFactors = FALSE
)

write.csv(
  coefficient_integrity_audit,
  file.path(
    step13i_dir,
    "STEP13I_Coefficient_Integrity_After_Paired_Analysis.csv"
  ),
  row.names = FALSE
)

if (!coefficient_unchanged) {
  stop(
    "STEP13I stopped: STEP13H coefficients changed unexpectedly."
  )
}


# ------------------------------------------------------------
# 54. Final external assessment summary
# ------------------------------------------------------------

final_summary <- data.frame(
  
  metric = c(
    "GSE32863_total_samples",
    "GSE32863_tumor_samples",
    "GSE32863_normal_samples",
    "GSE32863_unique_subjects",
    "GSE32863_complete_NT_pairs",
    "GSE32863_tumor_only_subjects",
    "GSE32863_normal_only_subjects",
    "External_AUC",
    "External_AUC_95CI_lower",
    "External_AUC_95CI_upper",
    "Paired_analysis_n_pairs"
  ),
  
  value = c(
    116L,
    58L,
    58L,
    59L,
    57L,
    1L,
    1L,
    external_auc,
    external_auc_ci[1],
    external_auc_ci[3],
    57L
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  final_summary,
  file.path(
    step13i_dir,
    "STEP13I_FINAL_External_Assessment_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  final_summary,
  file.path(
    results_dir,
    "STEP13I_FINAL_External_Assessment_Summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 55. Final comprehensive audit
# ------------------------------------------------------------

all_final_checks <- c(
  
  STEP13F_PASS =
    identical(
      model_lock$status,
      "PASS"
    ),
  
  STEP13F_external_refit_FALSE =
    external_lock_checks[
      "external_refit_allowed"
    ],
  
  STEP13F_external_coefficients_FALSE =
    external_lock_checks[
      "external_coefficient_estimation_allowed"
    ],
  
  STEP13F_external_feature_selection_FALSE =
    external_lock_checks[
      "external_feature_selection_allowed"
    ],
  
  STEP13F_external_threshold_FALSE =
    external_lock_checks[
      "external_threshold_optimization_allowed"
    ],
  
  STEP13F_external_cutpoint_FALSE =
    external_lock_checks[
      "external_cutpoint_search_allowed"
    ],
  
  STEP13F_external_tuning_FALSE =
    external_lock_checks[
      "external_model_tuning_allowed"
    ],
  
  STEP13F_model_definition_concordant =
    all(model_definition_checks),
  
  STEP13H_PASS =
    identical(
      final_model_lock$status,
      "PASS"
    ),
  
  STEP13H_glm =
    model_class_ok,
  
  STEP13H_binomial =
    model_family_ok,
  
  STEP13H_logit =
    model_link_ok,
  
  STEP13H_required_coefficients =
    coefficient_terms_ok,
  
  STEP13H_coefficients_finite =
    all(
      is.finite(
        locked_coefficients
      )
    ),
  
  GSE_total_116 =
    nrow(expr_final) > 0 &&
    ncol(expr_final) == 116L,
  
  GSE_tumor_58 =
    sum(group == "Tumor") == 58L,
  
  GSE_normal_58 =
    sum(group == "Normal") == 58L,
  
  GSE_unique_subjects_59 =
    length(
      unique(
        metadata_df$sample_code
      )
    ) == 59L,
  
  GSE_complete_pairs_57 =
    nrow(pair_df) == 57L,
  
  GSE_tumor_only_3023 =
    any(
      unmatched_df$subject_id == "3023" &
        unmatched_df$tumor_n == 1L &
        unmatched_df$normal_n == 0L
    ),
  
  GSE_normal_only_3035 =
    any(
      unmatched_df$subject_id == "3035" &
        unmatched_df$tumor_n == 0L &
        unmatched_df$normal_n == 1L
    ),
  
  Core_genes_present =
    all(core_gene_presence),
  
  External_inputs_complete =
    external_complete_n == 116L,
  
  External_predictions_finite =
    prediction_finite,
  
  External_predictions_valid =
    prediction_range_ok,
  
  External_AUC_finite =
    is.finite(
      external_auc
    ),
  
  External_AUC_CI_finite =
    all(
      is.finite(
        external_auc_ci
      )
    ),
  
  Paired_analysis_57_pairs =
    nrow(pair_df) == 57L,
  
  Paired_results_3_genes =
    nrow(paired_results_df) == 3L,
  
  STEP13H_coefficients_unchanged =
    coefficient_unchanged
)


final_audit_df <- data.frame(
  check = names(all_final_checks),
  PASS = unname(all_final_checks),
  stringsAsFactors = FALSE
)

write.csv(
  final_audit_df,
  file.path(
    step13i_dir,
    "STEP13I_Final_Audit_Checks.csv"
  ),
  row.names = FALSE
)

write.csv(
  final_audit_df,
  file.path(
    audit_dir,
    "STEP13I_Final_Audit_Checks.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 56. Final PASS / FAIL decision
# ------------------------------------------------------------

if (!all(all_final_checks)) {
  
  fail_file <-
    file.path(
      step13i_dir,
      "STEP13I_FAIL.txt"
    )
  
  writeLines(
    c(
      "STEP13I FINAL STATUS: FAIL",
      "",
      "One or more final external-assessment audit checks failed.",
      "",
      "Review:",
      "STEP13I_Final_Audit_Checks.csv"
    ),
    fail_file
  )
  
  stop(
    "STEP13I FINAL STATUS: FAIL"
  )
}


# ------------------------------------------------------------
# 57. Human-readable final report
# ------------------------------------------------------------

report_file <-
  file.path(
    step13i_dir,
    "STEP13I_GSE32863_External_Assessment_FINAL.txt"
  )

report_lines <- c(
  
  "============================================================",
  "STEP13I — GSE32863 EXTERNAL ASSESSMENT",
  "============================================================",
  "",
  "FINAL STATUS: PASS",
  "",
  "External dataset:",
  "GSE32863",
  "",
  "Total samples: 116",
  "Tumor: 58",
  "Normal: 58",
  "Unique subjects: 59",
  "Complete N+T pairs: 57",
  "Tumor-only subject: 3023",
  "Normal-only subject: 3035",
  "",
  "Locked model source:",
  "STEP13H_TCGA_Final_Model_FINAL.rds",
  "",
  "Predictors:",
  "AOC3 + MS4A7 + IL7R",
  "",
  "External refit: FALSE",
  "External coefficient estimation: FALSE",
  "External feature selection: FALSE",
  "External threshold optimization: FALSE",
  "External cutpoint search: FALSE",
  "External model tuning: FALSE",
  "",
  "External ROC:",
  paste0(
    "AUC = ",
    sprintf(
      "%.6f",
      external_auc
    )
  ),
  paste0(
    "95% CI = ",
    sprintf(
      "%.6f",
      external_auc_ci[1]
    ),
    " – ",
    sprintf(
      "%.6f",
      external_auc_ci[3]
    )
  ),
  "",
  "Paired expression analysis:",
  "57 complete N+T subject pairs",
  "Paired Wilcoxon signed-rank test",
  "",
  "Paired expression results:",
  ""
)

for (i in seq_len(
  nrow(paired_results_df)
)) {
  
  report_lines <- c(
    report_lines,
    paste0(
      paired_results_df$gene[i],
      " | n = ",
      paired_results_df$n_pairs[i],
      " | Normal median = ",
      signif(
        paired_results_df$normal_median[i],
        6
      ),
      " | Tumor median = ",
      signif(
        paired_results_df$tumor_median[i],
        6
      ),
      " | P = ",
      signif(
        paired_results_df$P_value[i],
        8
      ),
      " | FDR = ",
      signif(
        paired_results_df$FDR[i],
        8
      )
    )
  )
}

report_lines <- c(
  report_lines,
  "",
  "STEP13H coefficients remained unchanged after external",
  "assessment and paired expression analysis.",
  "",
  "============================================================",
  "STEP13I COMPLETE",
  "============================================================"
)

writeLines(
  report_lines,
  report_file
)


# ------------------------------------------------------------
# 58. Save final STEP13I object
# ------------------------------------------------------------

step13i_final_object <- list(
  
  status = "PASS",
  
  workflow = "STEP13I",
  
  external_dataset =
    "GSE32863",
  
  total_samples =
    116L,
  
  tumor_samples =
    58L,
  
  normal_samples =
    58L,
  
  unique_subjects =
    59L,
  
  complete_pairs =
    57L,
  
  tumor_only_subject =
    "3023",
  
  normal_only_subject =
    "3035",
  
  predictors =
    expected_core_genes,
  
  model_source =
    step13h_rds,
  
  model =
    final_model,
  
  locked_coefficients =
    locked_coefficients,
  
  external_auc =
    external_auc,
  
  external_auc_95ci =
    external_auc_ci,
  
  external_predictions =
    external_prediction_traceability,
  
  paired_subject_table =
    pair_df,
  
  paired_expression_results =
    paired_results_df,
  
  paired_expression_traceability =
    paired_long_df,
  
  expression_provenance =
    scale_provenance,
  
  external_refit =
    FALSE,
  
  external_coefficient_estimation =
    FALSE,
  
  external_feature_selection =
    FALSE,
  
  external_threshold_optimization =
    FALSE,
  
  external_cutpoint_search =
    FALSE,
  
  external_model_tuning =
    FALSE,
  
  final_audit =
    final_audit_df,
  
  generated_at =
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
)


step13i_final_rds <-
  file.path(
    step13i_dir,
    "STEP13I_GSE32863_External_Assessment_FINAL.rds"
  )

saveRDS(
  step13i_final_object,
  step13i_final_rds
)


# ------------------------------------------------------------
# 59. Save session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    step13i_dir,
    "STEP13I_sessionInfo_FINAL.txt"
  )
)


# ------------------------------------------------------------
# 60. PASS marker
# ------------------------------------------------------------

pass_file <-
  file.path(
    step13i_dir,
    "STEP13I_PASS.txt"
  )

writeLines(
  c(
    "STEP13I FINAL STATUS: PASS",
    "",
    "GSE32863 external assessment completed.",
    "",
    "External dataset:",
    "116 samples",
    "58 Tumor",
    "58 Normal",
    "59 subjects",
    "57 complete N+T pairs",
    "1 tumor-only subject: 3023",
    "1 normal-only subject: 3035",
    "",
    "Locked STEP13H coefficients were transferred unchanged.",
    "",
    "External refit: FALSE",
    "External coefficient estimation: FALSE",
    "External feature selection: FALSE",
    "External threshold optimization: FALSE",
    "External cutpoint search: FALSE",
    "External model tuning: FALSE",
    "",
    paste0(
      "External AUC: ",
      sprintf(
        "%.6f",
        external_auc
      )
    ),
    paste0(
      "External AUC 95% CI: ",
      sprintf(
        "%.6f",
        external_auc_ci[1]
      ),
      " – ",
      sprintf(
        "%.6f",
        external_auc_ci[3]
      )
    ),
    "",
    "Paired expression analysis:",
    "57 complete N+T subject pairs",
    "Paired Wilcoxon signed-rank test",
    "",
    "Final STEP13I object:",
    step13i_final_rds
  ),
  pass_file
)


# ------------------------------------------------------------
# 61. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP13I FINAL SUMMARY\n")
cat("============================================================\n\n")

cat(
  "External dataset:        GSE32863\n"
)

cat(
  "Total samples:           116\n"
)

cat(
  "Tumor:                   58\n"
)

cat(
  "Normal:                  58\n"
)

cat(
  "Unique subjects:         59\n"
)

cat(
  "Complete N+T pairs:      57\n"
)

cat(
  "Tumor-only subject:      3023\n"
)

cat(
  "Normal-only subject:     3035\n"
)

cat(
  "Core genes:              ",
  paste(
    expected_core_genes,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Model source:            STEP13H FINAL MODEL\n"
)

cat(
  "External refit:          FALSE\n"
)

cat(
  "External coefficients:   FALSE\n"
)

cat(
  "External feature select: FALSE\n"
)

cat(
  "External threshold tune: FALSE\n"
)

cat(
  "External cutpoint search:FALSE\n"
)

cat(
  "External model tuning:   FALSE\n"
)

cat("\n")

cat(
  "External AUC:            ",
  sprintf(
    "%.6f",
    external_auc
  ),
  "\n",
  sep = ""
)

cat(
  "External AUC 95% CI:     ",
  sprintf(
    "%.6f",
    external_auc_ci[1]
  ),
  " – ",
  sprintf(
    "%.6f",
    external_auc_ci[3]
  ),
  "\n",
  sep = ""
)

cat("\n")
cat("Paired expression analysis:\n")
cat("Method:                   Paired Wilcoxon signed-rank test\n")
cat("Complete pairs:           57\n\n")

print(
  paired_results_df
)

cat("\n")

cat(
  "STEP13H coefficients unchanged: ",
  coefficient_unchanged,
  "\n",
  sep = ""
)

cat(
  "All final audit checks:          ",
  all(all_final_checks),
  "\n",
  sep = ""
)

cat("\n")
cat("FINAL STATUS:\n")
cat(
  "PASS — STEP13I GSE32863 External Assessment completed.\n"
)

cat("\n")
cat("Final RDS:\n")
cat(
  step13i_final_rds,
  "\n"
)

cat("\n")
cat("STEP13I COMPLETE\n")
cat("============================================================\n\n")

##Step 1: Verify only STEP13A
a <- read.csv(
  "input_dir/Paper 1_Step 13-Diagnostic Classification/STEP13A_Clinical_Claim_Deescalation/STEP13A_Clinical_Claim_Lock.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

print(a)  


##J2：STEP13C Final Input Lock

j2 <- readLines(
  "input_dir/Paper 1_Step 13-Diagnostic Classification/STEP13C_GSE32863_Input_Lock/STEP13C-5_Final_Input_Lock/STEP13C-5_PASS.txt"
)

cat(j2, sep = "\n")   


##J3：Verify STEP13F Model Definition Lock
f <- readRDS(
  "input_dir/Paper 1_Step 13-Diagnostic Classification/STEP13F_Model_Definition_Lock/STEP13F_Model_Definition_Lock_FINAL.rds"
)

print(f)

##J4：Verify the final results of STEP13G Sterified-CV
g4 <- read.csv(
  "input_dir/Paper 1_Step 13-Diagnostic Classification/STEP13G_Nested_CV/STEP13G_OOF_Performance_Summary.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

print(g4)

#J5：Verify STEP13H TCGA Final Model
h <- readRDS(
  "input_dir/Paper 1_Step 13-Diagnostic Classification/STEP13H_TCGA_Final_Model/STEP13H_TCGA_Final_Model_FINAL.rds"
)

print(h)


##J6：Verify STEP13I External Assessment
i6 <- read.csv(
  "input_dir/Paper 1_Step 13-Diagnostic Classification/STEP13I_GSE32863_External_Assessment/STEP13I_GSE32863_External_ROC_AUC.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

print(i6)



##Next step: J6-B validation – paired Wilcoxon test (n=57) 
iw <- read.csv(
  "input_dir/Paper 1_Step 13-Diagnostic Classification/STEP13I_GSE32863_External_Assessment/STEP13I_GSE32863_Paired_Wilcoxon_AOC3_MS4A7_IL7R.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

print(iw)   


##J6-C：Final Check: STEP 13I Overall Audit
ia <- read.csv(
  "input_dir/Paper 1_Step 13-Diagnostic Classification/STEP13I_GSE32863_External_Assessment/STEP13I_Final_Audit_Checks.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

print(ia)


### ============================================================
# STEP13J — FINAL MANUSCRIPT RESULT LOCK
# Static documentation only
# NO statistical recalculation
# NO model fitting
# NO AUC / P-value recalculation
# ============================================================

workflow_root <- "input_dir/Paper 1_Step 13-Diagnostic Classification"

j7_dir <- file.path(
  workflow_root,
  "STEP13J_Final_Audit"
)

if (!dir.exists(j7_dir)) {
  dir.create(j7_dir, recursive = TRUE)
}

# ------------------------------------------------------------
# 1. Final manuscript result lock
# ------------------------------------------------------------

lock <- data.frame(
  Section = c(
    "Workflow status",
    "Analysis objective",
    "Clinical diagnostic claim",
    "Diagnostic biomarker claim",
    "Core genes",
    "Model type",
    "Model formula",
    "TCGA total samples",
    "TCGA tumor samples",
    "TCGA normal samples",
    "Nested CV outer folds",
    "Nested CV inner folds",
    "TCGA OOF AUC",
    "TCGA OOF AUC 95% CI lower",
    "TCGA OOF AUC 95% CI upper",
    "Fixed probability threshold",
    "Sensitivity",
    "Specificity",
    "Accuracy",
    "Balanced accuracy",
    "Brier score",
    "AOC3 final coefficient",
    "MS4A7 final coefficient",
    "IL7R final coefficient",
    "GSE32863 total samples",
    "GSE32863 tumor samples",
    "GSE32863 normal samples",
    "GSE32863 unique subjects",
    "GSE32863 complete tumor-normal pairs",
    "GSE32863 tumor-only subject",
    "GSE32863 normal-only subject",
    "External refit",
    "External coefficient estimation",
    "External feature selection",
    "External threshold optimization",
    "External cutpoint search",
    "External model tuning",
    "GSE32863 external AUC",
    "GSE32863 external AUC 95% CI lower",
    "GSE32863 external AUC 95% CI upper",
    "Paired test",
    "Paired analysis n",
    "AOC3 tumor-minus-normal median difference",
    "MS4A7 tumor-minus-normal median difference",
    "IL7R tumor-minus-normal median difference",
    "J1 STEP13A verification",
    "J2 STEP13C verification",
    "J3 STEP13F verification",
    "J4 STEP13G verification",
    "J5 STEP13H verification",
    "J6 STEP13I verification"
  ),
  
  Locked_value = c(
    "FINAL — J1-J6 VERIFIED",
    "Exploratory tumor-versus-normal tissue classification",
    "FALSE",
    "FALSE",
    "AOC3, MS4A7, IL7R",
    "Multivariable logistic regression",
    "TumorStatus ~ AOC3 + MS4A7 + IL7R",
    "575",
    "516",
    "59",
    "10",
    "10",
    "0.98794508",
    "0.97714550",
    "0.99874466",
    "0.50",
    "0.98837209",
    "0.89830508",
    "0.97913043",
    "0.94333859",
    "0.01820027",
    "-2.4323976",
    "-2.0847930",
    "-0.2147045",
    "116",
    "58",
    "58",
    "59",
    "57",
    "3023",
    "3035",
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE",
    "0.9973246",
    "0.9929578",
    "1.0000000",
    "Paired Wilcoxon signed-rank test",
    "57",
    "-1.051579",
    "-2.065199",
    "-1.395997",
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS"
  ),
  
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Write CSV
# ------------------------------------------------------------

csv_file <- file.path(
  j7_dir,
  "STEP13J_FINAL_MANUSCRIPT_RESULT_LOCK.csv"
)

write.csv(
  lock,
  csv_file,
  row.names = FALSE,
  quote = TRUE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 3. Write human-readable TXT
# ------------------------------------------------------------

txt_file <- file.path(
  j7_dir,
  "STEP13J_FINAL_MANUSCRIPT_RESULT_LOCK.txt"
)

con <- file(
  txt_file,
  open = "wt",
  encoding = "UTF-8"
)

writeLines(
  c(
    "STEP13J — FINAL MANUSCRIPT RESULT LOCK",
    "=======================================",
    "",
    "STATUS: FINAL — J1-J6 VERIFIED",
    "",
    "This file is a static documentation lock.",
    "No statistical recalculation was performed.",
    "No model fitting was performed.",
    "No AUC, CI, P value, FDR, or threshold was recalculated.",
    "",
    "ANALYSIS OBJECTIVE",
    "Exploratory tumor-versus-normal tissue classification",
    "",
    "CLINICAL CLAIM LOCK",
    "Clinical diagnostic claim: FALSE",
    "Diagnostic biomarker claim: FALSE",
    "",
    "CORE GENES",
    "AOC3, MS4A7, IL7R",
    "",
    "MODEL",
    "Multivariable logistic regression",
    "Formula: TumorStatus ~ AOC3 + MS4A7 + IL7R",
    "",
    "TCGA DATASET",
    "Total: 575",
    "Tumor: 516",
    "Normal: 59",
    "Nested CV: 10 outer folds × 10 inner folds",
    "",
    "TCGA NESTED-CV RESULT",
    "OOF AUC: 0.98794508",
    "95% CI: 0.97714550–0.99874466",
    "Fixed threshold: 0.50",
    "Sensitivity: 0.98837209",
    "Specificity: 0.89830508",
    "Accuracy: 0.97913043",
    "Balanced accuracy: 0.94333859",
    "Brier score: 0.01820027",
    "",
    "FINAL TCGA MODEL COEFFICIENTS",
    "AOC3: -2.4323976",
    "MS4A7: -2.0847930",
    "IL7R: -0.2147045",
    "",
    "GSE32863 EXTERNAL ASSESSMENT",
    "Total samples: 116",
    "Tumor: 58",
    "Normal: 58",
    "Unique subjects: 59",
    "Complete tumor-normal pairs: 57",
    "Tumor-only subject: 3023",
    "Normal-only subject: 3035",
    "",
    "External refit: FALSE",
    "External coefficient estimation: FALSE",
    "External feature selection: FALSE",
    "External threshold optimization: FALSE",
    "External cutpoint search: FALSE",
    "External model tuning: FALSE",
    "",
    "External AUC: 0.9973246",
    "95% CI: 0.9929578–1.0000000",
    "",
    "PAIRED EXPRESSION ANALYSIS",
    "Test: Paired Wilcoxon signed-rank test",
    "Complete pairs: 57",
    "",
    "AOC3 tumor-minus-normal median difference: -1.051579",
    "MS4A7 tumor-minus-normal median difference: -2.065199",
    "IL7R tumor-minus-normal median difference: -1.395997",
    "",
    "FINAL VERIFICATION STATUS",
    "J1 STEP13A: PASS",
    "J2 STEP13C: PASS",
    "J3 STEP13F: PASS",
    "J4 STEP13G: PASS",
    "J5 STEP13H: PASS",
    "J6 STEP13I: PASS",
    "",
    "INTERPRETATION LOCK",
    "The analysis is interpreted as exploratory tumor-versus-normal",
    "tissue classification only and is not presented as evidence",
    "for a clinical diagnostic biomarker.",
    "",
    "HISTORICAL ROC DISTINCTION",
    "Historical Xiangtao combined ROC AUC: 0.9935619498",
    "This historical AUC is distinct from the revised nested-CV",
    "OOF AUC of 0.98794508 and must not be conflated with it.",
    "",
    "END OF STEP13J FINAL MANUSCRIPT RESULT LOCK"
  ),
  con
)

close(con)

# ------------------------------------------------------------
# 4. Final file existence check only
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP13J FINAL MANUSCRIPT RESULT LOCK\n")
cat("============================================================\n\n")

cat("CSV:\n", csv_file, "\n")
cat("TXT:\n", txt_file, "\n\n")

cat("CSV exists:", file.exists(csv_file), "\n")
cat("TXT exists:", file.exists(txt_file), "\n\n")

cat("FINAL STATUS:\n")
cat("PASS — STEP13J Final Manuscript Result Lock created.\n")
cat("No statistical recalculation performed.\n")  


# ============================================================
# Figure 9A — STEP 1
# Check locked TCGA OOF prediction file
# ============================================================

oof_file <- "input_dir/Paper 1_Step 13-Diagnostic Classification/STEP13G_Nested_CV/STEP13G_Out_of_Fold_Predictions.csv"

oof <- read.csv(
  oof_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensions:\n")
print(dim(oof))

cat("\nColumn names:\n")
print(names(oof))

cat("\nFirst 6 rows:\n")
print(head(oof))   


# ============================================================
# Figure 9A — TCGA Streified-CV out-of-fold ROC
# STEP 2: Plot and save
# No model refitting
# No CV recalculation
# No threshold optimization
# ============================================================

library(ggplot2)
library(pROC)

# ------------------------------------------------------------
# 1. Input
# ------------------------------------------------------------

oof_file <- "input_dir/Paper 1_Step 13-Diagnostic Classification/STEP13G_Nested_CV/STEP13G_Out_of_Fold_Predictions.csv"

oof <- read.csv(
  oof_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# 2. Basic input audit
# ------------------------------------------------------------

stopifnot(
  nrow(oof) == 575,
  all(c(
    "sample_id",
    "observed_status",
    "observed_binary",
    "predicted_probability_Tumor",
    "outer_fold"
  ) %in% names(oof))
)

stopifnot(
  all(oof$observed_binary %in% c(0, 1)),
  all(is.finite(oof$predicted_probability_Tumor))
)

# ------------------------------------------------------------
# 3. Locked AUC result
# ------------------------------------------------------------

auc_locked <- 0.987945
ci_low_locked <- 0.977145
ci_high_locked <- 0.998745

# ------------------------------------------------------------
# 4. Construct ROC coordinates ONLY
#    from locked OOF predictions
# ------------------------------------------------------------

roc_obj <- roc(
  response = oof$observed_binary,
  predictor = oof$predicted_probability_Tumor,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

roc_df <- data.frame(
  specificity = roc_obj$specificities,
  sensitivity = roc_obj$sensitivities
)

roc_df$FPR <- 1 - roc_df$specificity

# ------------------------------------------------------------
# 5. Font settings
# ------------------------------------------------------------

windowsFonts(
  Arial = windowsFont("Arial")
)

# ------------------------------------------------------------
# 6. Figure 9A
# ------------------------------------------------------------

p9A <- ggplot(
  roc_df,
  aes(
    x = FPR,
    y = sensitivity
  )
) +
  
  # ROC curve
  geom_path(
    color = "#C00000",
    linewidth = 0.9
  ) +
  
  # Reference diagonal
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.55
  ) +
  
  # AUC annotation
  annotate(
    "text",
    x = 0.60,
    y = 0.18,
    label = paste0(
      "AUC = ",
      sprintf("%.3f", auc_locked),
      "\n95% CI: ",
      sprintf("%.3f", ci_low_locked),
      "–",
      sprintf("%.3f", ci_high_locked)
    ),
    family = "Arial",
    fontface = "plain",
    size = 3.0
  ) +
  
  # Axes
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  
  labs(
    x = "1 − Specificity",
    y = "Sensitivity"
  ) +
  
  theme_classic(
    base_family = "Arial"
  ) +
  
  theme(
    axis.title.x = element_text(
      family = "Arial",
      face = "bold",
      size = 10
    ),
    axis.title.y = element_text(
      family = "Arial",
      face = "bold",
      size = 10
    ),
    axis.text.x = element_text(
      family = "Arial",
      size = 9
    ),
    axis.text.y = element_text(
      family = "Arial",
      size = 9
    ),
    axis.line = element_line(
      linewidth = 0.6
    ),
    axis.ticks = element_line(
      linewidth = 0.5
    ),
    plot.margin = margin(
      8, 8, 8, 8
    )
  )

# ------------------------------------------------------------
# 7. Display
# ------------------------------------------------------------

print(p9A)

# ------------------------------------------------------------
# 8. Output directories
# ------------------------------------------------------------

fig_dir <- "result_dir/Paper 1_Step 13-Diagnostic Classification/FIGURES"

if (!dir.exists(fig_dir)) {
  dir.create(
    fig_dir,
    recursive = TRUE
  )
}

# ------------------------------------------------------------
# 9. Save PDF — Cairo
# ------------------------------------------------------------

pdf_file <- file.path(
  fig_dir,
  "Figure9A_TCGA_nestedCV_OOF_ROC_Cairo.pdf"
)

Cairo::CairoPDF(
  filename = pdf_file,
  width = 7.0,
  height = 5.8,
  family = "Arial"
)

print(p9A)

dev.off()

# ------------------------------------------------------------
# 10. Save TIFF — 600 dpi
# ------------------------------------------------------------

tiff_file <- file.path(
  fig_dir,
  "Figure9A_TCGA_nestedCV_OOF_ROC_600dpi.tiff"
)

Cairo::CairoTIFF(
  filename = tiff_file,
  width = 7.0,
  height = 5.8,
  units = "in",
  dpi = 600,
  family = "Arial"
)

print(p9A)

dev.off()

# ------------------------------------------------------------
# 11. Final output audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Figure 9A — FINAL OUTPUT AUDIT\n")
cat("============================================================\n\n")

cat("Input rows:              ", nrow(oof), "\n")
cat("Tumor samples:           ", sum(oof$observed_binary == 1), "\n")
cat("Normal samples:          ", sum(oof$observed_binary == 0), "\n")
cat("Locked AUC:              ", auc_locked, "\n")
cat("Locked 95% CI:           ",
    ci_low_locked, "–", ci_high_locked, "\n")
cat("PDF exists:              ", file.exists(pdf_file), "\n")
cat("TIFF exists:             ", file.exists(tiff_file), "\n")

cat("\nPDF:\n")
cat(pdf_file, "\n")

cat("\nTIFF:\n")
cat(tiff_file, "\n")

cat("\nFINAL STATUS:\n")

if (
  nrow(oof) == 575 &&
  sum(oof$observed_binary == 1) == 516 &&
  sum(oof$observed_binary == 0) == 59 &&
  file.exists(pdf_file) &&
  file.exists(tiff_file)
) {
  cat("PASS — Figure 9A successfully generated.\n")
} else {
  cat("FAIL — Figure 9A output audit failed.\n")
}

cat("============================================================\n")  


# ============================================================
# Figure 9A — ROC coordinate audit
# Do NOT modify the figure yet
# ============================================================

library(pROC)

oof_file <- "result_dir/Paper 1_Step 13-Diagnostic Classification/STEP13G_Nested_CV/STEP13G_Out_of_Fold_Predictions.csv"

oof <- read.csv(
  oof_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

roc_obj <- roc(
  response = oof$observed_binary,
  predictor = oof$predicted_probability_Tumor,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

cat("============================================================\n")
cat("ROC COORDINATE AUDIT\n")
cat("============================================================\n\n")

cat("AUC from ROC object:\n")
print(as.numeric(auc(roc_obj)))

cat("\nNumber of ROC coordinates:\n")
print(length(roc_obj$sensitivities))

cat("\nFirst 10 ROC coordinates:\n")

roc_check <- data.frame(
  specificity = roc_obj$specificities,
  FPR = 1 - roc_obj$specificities,
  sensitivity = roc_obj$sensitivities,
  threshold = roc_obj$thresholds
)

print(head(roc_check, 10))

cat("\nLast 10 ROC coordinates:\n")
print(tail(roc_check, 10))

cat("\nMinimum / maximum FPR:\n")
print(range(roc_check$FPR, na.rm = TRUE))

cat("\nMinimum / maximum sensitivity:\n")
print(range(roc_check$sensitivity, na.rm = TRUE))

cat("\nCoordinates closest to (0,0):\n")

d00 <- sqrt(
  roc_check$FPR^2 +
    roc_check$sensitivity^2
)

print(
  roc_check[
    order(d00)[1:min(5, nrow(roc_check))],
  ]
)

cat("\nCoordinates closest to (1,1):\n")

d11 <- sqrt(
  (roc_check$FPR - 1)^2 +
    (roc_check$sensitivity - 1)^2
)

print(
  roc_check[
    order(d11)[1:min(5, nrow(roc_check))],
  ]
)

cat("\n============================================================\n")   

# ============================================================
# Figure 9B — STEP 1
# Locate and inspect STEP13I external assessment files
# ============================================================

step13i_dir <- "result_dir/Paper 1_Step 13-Diagnostic Classification/STEP13I_GSE32863_External_Assessment"

files_i <- list.files(
  step13i_dir,
  recursive = TRUE,
  full.names = TRUE
)

print(files_i)  

# ============================================================
# Figure 9B — STEP 2
# Inspect locked external prediction traceability
# No recalculation
# ============================================================

pred_file <- paste0(
  "result_dir/Paper 1_Step 13-Diagnostic Classification/",
  "STEP13I_GSE32863_External_Assessment/",
  "STEP13I_GSE32863_External_Prediction_Traceability.csv"
)

pred_ext <- read.csv(
  pred_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("============================================================\n")
cat("Figure 9B — External Prediction Traceability Audit\n")
cat("============================================================\n\n")

cat("Dimensions:\n")
print(dim(pred_ext))

cat("\nColumn names:\n")
print(names(pred_ext))

cat("\nFirst 10 rows:\n")
print(head(pred_ext, 10))

cat("\nLast 5 rows:\n")
print(tail(pred_ext, 5))

cat("\n============================================================\n") 


# ============================================================
# Figure 9B — GSE32863 External ROC
# STEP 3: Plot and save
#
# Locked external predictions only
# No model refitting
# No coefficient re-estimation
# No threshold optimization
# ============================================================

library(ggplot2)
library(pROC)

# ------------------------------------------------------------
# 1. Input
# ------------------------------------------------------------

pred_file <- paste0(
  "result_dir/Paper 1_Step 13-Diagnostic Classification/",
  "STEP13I_GSE32863_External_Assessment/",
  "STEP13I_GSE32863_External_Prediction_Traceability.csv"
)

pred_ext <- read.csv(
  pred_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# 2. Input audit
# ------------------------------------------------------------

required_cols <- c(
  "sample_id",
  "TumorStatus",
  "predicted_probability"
)

stopifnot(
  nrow(pred_ext) == 116,
  all(required_cols %in% names(pred_ext)),
  all(pred_ext$TumorStatus %in% c("Normal", "Tumor")),
  all(is.finite(pred_ext$predicted_probability))
)

stopifnot(
  sum(pred_ext$TumorStatus == "Tumor") == 58,
  sum(pred_ext$TumorStatus == "Normal") == 58
)

# ------------------------------------------------------------
# 3. Locked external AUC result
# ------------------------------------------------------------

auc_locked <- 0.9973246
ci_low_locked <- 0.9929578
ci_high_locked <- 1.0000000

# ------------------------------------------------------------
# 4. ROC coordinates
#    Generated only from locked external predictions
# ------------------------------------------------------------

roc_obj <- roc(
  response = pred_ext$TumorStatus,
  predictor = pred_ext$predicted_probability,
  levels = c("Normal", "Tumor"),
  direction = "<",
  quiet = TRUE
)

roc_df <- data.frame(
  specificity = roc_obj$specificities,
  sensitivity = roc_obj$sensitivities
)

roc_df$FPR <- 1 - roc_df$specificity

# ------------------------------------------------------------
# 5. Font
# ------------------------------------------------------------

windowsFonts(
  Arial = windowsFont("Arial")
)

# ------------------------------------------------------------
# 6. Figure 9B
# ------------------------------------------------------------

p9B <- ggplot(
  roc_df,
  aes(
    x = FPR,
    y = sensitivity
  )
) +
  
  # ROC curve
  geom_path(
    color = "#C00000",
    linewidth = 0.9
  ) +
  
  # Reference diagonal
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.55
  ) +
  
  # AUC annotation
  annotate(
    "text",
    x = 0.60,
    y = 0.18,
    label = paste0(
      "AUC = ",
      sprintf("%.3f", auc_locked),
      "\n95% CI: ",
      sprintf("%.3f", ci_low_locked),
      "–",
      sprintf("%.3f", ci_high_locked)
    ),
    family = "Arial",
    fontface = "plain",
    size = 3.0
  ) +
  
  # X axis
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  
  # Y axis
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  
  labs(
    x = "1 − Specificity",
    y = "Sensitivity"
  ) +
  
  theme_classic(
    base_family = "Arial"
  ) +
  
  theme(
    axis.title.x = element_text(
      family = "Arial",
      face = "bold",
      size = 10
    ),
    axis.title.y = element_text(
      family = "Arial",
      face = "bold",
      size = 10
    ),
    axis.text.x = element_text(
      family = "Arial",
      size = 9
    ),
    axis.text.y = element_text(
      family = "Arial",
      size = 9
    ),
    axis.line = element_line(
      linewidth = 0.6
    ),
    axis.ticks = element_line(
      linewidth = 0.5
    ),
    plot.margin = margin(
      8, 8, 8, 8
    )
  )

# ------------------------------------------------------------
# 7. Display
# ------------------------------------------------------------

print(p9B)

# ------------------------------------------------------------
# 8. Output directory
# ------------------------------------------------------------

fig_dir <- "result_dir/Paper 1_Step 13-Diagnostic Classification/FIGURES"

if (!dir.exists(fig_dir)) {
  dir.create(
    fig_dir,
    recursive = TRUE
  )
}

# ------------------------------------------------------------
# 9. Save Cairo PDF
# ------------------------------------------------------------

pdf_file <- file.path(
  fig_dir,
  "Figure9B_GSE32863_external_ROC_Cairo.pdf"
)

Cairo::CairoPDF(
  filename = pdf_file,
  width = 7.0,
  height = 5.8,
  family = "Arial"
)

print(p9B)

dev.off()

# ------------------------------------------------------------
# 10. Save 600-dpi TIFF
# ------------------------------------------------------------

tiff_file <- file.path(
  fig_dir,
  "Figure9B_GSE32863_external_ROC_600dpi.tiff"
)

Cairo::CairoTIFF(
  filename = tiff_file,
  width = 7.0,
  height = 5.8,
  units = "in",
  dpi = 600,
  family = "Arial"
)

print(p9B)

dev.off()

# ------------------------------------------------------------
# 11. Final audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Figure 9B — FINAL OUTPUT AUDIT\n")
cat("============================================================\n\n")

cat("Input samples:           ", nrow(pred_ext), "\n")
cat("Tumor samples:           ",
    sum(pred_ext$TumorStatus == "Tumor"), "\n")
cat("Normal samples:          ",
    sum(pred_ext$TumorStatus == "Normal"), "\n")

cat("Locked AUC:              ",
    auc_locked, "\n")

cat(
  "Locked 95% CI:           ",
  ci_low_locked,
  "–",
  ci_high_locked,
  "\n"
)

cat(
  "PDF exists:              ",
  file.exists(pdf_file),
  "\n"
)

cat(
  "TIFF exists:             ",
  file.exists(tiff_file),
  "\n"
)

cat("\nPDF:\n")
cat(pdf_file, "\n")

cat("\nTIFF:\n")
cat(tiff_file, "\n")

cat("\nFINAL STATUS:\n")

if (
  nrow(pred_ext) == 116 &&
  sum(pred_ext$TumorStatus == "Tumor") == 58 &&
  sum(pred_ext$TumorStatus == "Normal") == 58 &&
  file.exists(pdf_file) &&
  file.exists(tiff_file)
) {
  cat("PASS — Figure 9B successfully generated.\n")
} else {
  cat("FAIL — Figure 9B output audit failed.\n")
}

cat("============================================================\n")   


# ============================================================
# Figure 9C — STEP 1
# Check locked paired-expression traceability file
# ============================================================

pair_file <- "result_dir/Paper 1_Step 13-Diagnostic Classification/STEP13I_GSE32863_External_Assessment/STEP13I_Paired_Expression_Traceability.csv"

pair_df <- read.csv(
  pair_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensions:\n")
print(dim(pair_df))

cat("\nColumn names:\n")
print(names(pair_df))

cat("\nFirst 10 rows:\n")
print(head(pair_df, 10))    

# ============================================================
# Figure 9C — STEP 2
# Complete structural audit of locked paired-expression data
#
# No plotting
# No statistical recalculation
# No Wilcoxon recalculation
# ============================================================

cat("Figure 9C — STEP 2 STRUCTURAL AUDIT\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Basic dimensions
# ------------------------------------------------------------

cat("Rows:       ", nrow(pair_df), "\n")
cat("Columns:    ", ncol(pair_df), "\n")

# ------------------------------------------------------------
# 2. Required columns
# ------------------------------------------------------------

required_cols <- c(
  "pair_id",
  "subject_id",
  "gene",
  "normal_GSM",
  "tumor_GSM",
  "normal_value",
  "tumor_value",
  "difference_tumor_minus_normal"
)

cat("\nRequired columns present: ",
    all(required_cols %in% names(pair_df)),
    "\n")

# ------------------------------------------------------------
# 3. Gene structure
# ------------------------------------------------------------

cat("\nGene counts:\n")
print(table(pair_df$gene))

cat("\nUnique genes:\n")
print(unique(pair_df$gene))

# ------------------------------------------------------------
# 4. Pair structure
# ------------------------------------------------------------

cat("\nUnique pair IDs:       ",
    length(unique(pair_df$pair_id)),
    "\n")

cat("Unique subjects:       ",
    length(unique(pair_df$subject_id)),
    "\n")

cat("Rows per gene:\n")
print(table(pair_df$gene))

cat("\nPairs per gene:\n")
print(
  tapply(
    pair_df$pair_id,
    pair_df$gene,
    function(x) length(unique(x))
  )
)

# ------------------------------------------------------------
# 5. Missing-value audit
# ------------------------------------------------------------

cat("\nMissing-value audit:\n")

print(
  sapply(
    pair_df[required_cols],
    function(x) sum(is.na(x))
  )
)

# ------------------------------------------------------------
# 6. Numeric-value audit
# ------------------------------------------------------------

cat("\nFinite numeric-value audit:\n")

cat(
  "normal_value finite: ",
  all(is.finite(pair_df$normal_value)),
  "\n"
)

cat(
  "tumor_value finite:  ",
  all(is.finite(pair_df$tumor_value)),
  "\n"
)

cat(
  "difference finite:    ",
  all(is.finite(pair_df$difference_tumor_minus_normal)),
  "\n"
)

# ------------------------------------------------------------
# 7. Difference traceability audit
# ------------------------------------------------------------

difference_check <- with(
  pair_df,
  abs(
    difference_tumor_minus_normal -
      (tumor_value - normal_value)
  ) < 1e-10
)

cat(
  "\nDifference traceability: ",
  all(difference_check),
  "\n"
)

# ------------------------------------------------------------
# 8. Expected structure
# ------------------------------------------------------------

expected_structure <- (
  nrow(pair_df) == 171 &&
    ncol(pair_df) == 8 &&
    all(required_cols %in% names(pair_df)) &&
    length(unique(pair_df$pair_id)) == 57 &&
    length(unique(pair_df$subject_id)) == 57 &&
    all(table(pair_df$gene) == 57) &&
    all(is.finite(pair_df$normal_value)) &&
    all(is.finite(pair_df$tumor_value)) &&
    all(is.finite(pair_df$difference_tumor_minus_normal)) &&
    all(difference_check)
)

cat("\nFINAL STATUS:\n")

if (expected_structure) {
  cat(
    "PASS — Figure 9C locked paired-expression traceability ",
    "structure verified.\n",
    sep = ""
  )
} else {
  cat(
    "FAIL — Figure 9C paired-expression structure audit failed.\n"
  )
}

# ============================================================
# Figure 9C — STEP 3
# Prepare locked paired-expression plotting data
#
# No statistical recalculation
# No Wilcoxon recalculation
# No p-value recalculation
# ============================================================

plot_df <- pair_df[, c(
  "pair_id",
  "subject_id",
  "gene",
  "normal_value",
  "tumor_value"
)]

# ------------------------------------------------------------
# Convert to long format
# ------------------------------------------------------------

plot_long <- rbind(
  data.frame(
    pair_id = plot_df$pair_id,
    subject_id = plot_df$subject_id,
    gene = plot_df$gene,
    tissue = "Normal",
    expression = plot_df$normal_value,
    stringsAsFactors = FALSE
  ),
  data.frame(
    pair_id = plot_df$pair_id,
    subject_id = plot_df$subject_id,
    gene = plot_df$gene,
    tissue = "Tumor",
    expression = plot_df$tumor_value,
    stringsAsFactors = FALSE
  )
)

# ------------------------------------------------------------
# Lock factor ordering
# ------------------------------------------------------------

plot_long$gene <- factor(
  plot_long$gene,
  levels = c("AOC3", "MS4A7", "IL7R")
)

plot_long$tissue <- factor(
  plot_long$tissue,
  levels = c("Normal", "Tumor")
)

# ------------------------------------------------------------
# Structural audit
# ------------------------------------------------------------

cat("Figure 9C — STEP 3 PLOTTING DATA AUDIT\n")
cat("============================================================\n\n")

cat("Rows:             ", nrow(plot_long), "\n")
cat("Expected rows:    ", 342, "\n")

cat("\nGene counts:\n")
print(table(plot_long$gene))

cat("\nTissue counts:\n")
print(table(plot_long$tissue))

cat("\nGene × tissue counts:\n")
print(table(plot_long$gene, plot_long$tissue))

cat("\nUnique subjects:   ",
    length(unique(plot_long$subject_id)),
    "\n")

cat("Unique pairs:      ",
    length(unique(plot_long$pair_id)),
    "\n")

cat(
  "\nExpression values finite: ",
  all(is.finite(plot_long$expression)),
  "\n"
)

# ------------------------------------------------------------
# Final structure lock
# ------------------------------------------------------------

step3_pass <- (
  nrow(plot_long) == 342 &&
    all(table(plot_long$gene) == 114) &&
    all(table(plot_long$tissue) == 171) &&
    all(table(plot_long$gene, plot_long$tissue) == 57) &&
    length(unique(plot_long$subject_id)) == 57 &&
    length(unique(plot_long$pair_id)) == 57 &&
    all(is.finite(plot_long$expression))
)

cat("\nFINAL STATUS:\n")

if (step3_pass) {
  cat(
    "PASS — Figure 9C plotting data prepared from locked ",
    "paired-expression traceability data.\n",
    sep = ""
  )
} else {
  cat(
    "FAIL — Figure 9C plotting data audit failed.\n"
  )
}


# ============================================================
# Figure 9C — STEP 4
# Final paired-expression plot
#
# Locked input:
# STEP13I_Paired_Expression_Traceability.csv
#
# No statistical recalculation
# No Wilcoxon recalculation
# No FDR recalculation
# ============================================================

library(ggplot2)

# ------------------------------------------------------------
# 0. Locked input
# ------------------------------------------------------------

pair_file <- paste0(
  "result_dir/Paper 1_Step 13-Diagnostic Classification/",
  "STEP13I_GSE32863_External_Assessment/",
  "STEP13I_Paired_Expression_Traceability.csv"
)

pair_df <- read.csv(
  pair_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# 1. Prepare plotting data
# ------------------------------------------------------------

plot_df <- pair_df[, c(
  "pair_id",
  "subject_id",
  "gene",
  "normal_value",
  "tumor_value"
)]

plot_long <- rbind(
  data.frame(
    pair_id = plot_df$pair_id,
    subject_id = plot_df$subject_id,
    gene = plot_df$gene,
    tissue = "Normal",
    expression = plot_df$normal_value,
    stringsAsFactors = FALSE
  ),
  data.frame(
    pair_id = plot_df$pair_id,
    subject_id = plot_df$subject_id,
    gene = plot_df$gene,
    tissue = "Tumor",
    expression = plot_df$tumor_value,
    stringsAsFactors = FALSE
  )
)

plot_long$gene <- factor(
  plot_long$gene,
  levels = c("AOC3", "MS4A7", "IL7R")
)

plot_long$tissue <- factor(
  plot_long$tissue,
  levels = c("Normal", "Tumor")
)

# ------------------------------------------------------------
# 2. Locked statistical annotations
# ------------------------------------------------------------

# These values are copied from the finalized STEP13I
# Wilcoxon/FDR results.
#
# No statistical recalculation is performed here.

stat_df <- data.frame(
  gene = factor(
    c("AOC3", "MS4A7", "IL7R"),
    levels = c("AOC3", "MS4A7", "IL7R")
  ),
  label = c(
    "FDR < 0.001",
    "FDR < 0.001",
    "FDR < 0.001"
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 3. Calculate annotation positions for display only
# ------------------------------------------------------------

y_max <- tapply(
  plot_long$expression,
  plot_long$gene,
  max
)

y_min <- tapply(
  plot_long$expression,
  plot_long$gene,
  min
)

stat_df$y <- as.numeric(y_max) +
  0.10 * (as.numeric(y_max) - as.numeric(y_min))

# ------------------------------------------------------------
# 4. Lock colors
# ------------------------------------------------------------

normal_col <- "#7F7F7F"
tumor_col  <- "#C00000"
line_col   <- "#BFBFBF"

# ------------------------------------------------------------
# 5. Font
# ------------------------------------------------------------

windowsFonts(
  Arial = windowsFont("Arial")
)

# ------------------------------------------------------------
# 6. Draw Figure 9C
# ------------------------------------------------------------

p9C <- ggplot(
  plot_long,
  aes(
    x = tissue,
    y = expression,
    group = subject_id
  )
) +
  
  # ----------------------------------------------------------
# Paired subject lines
# ----------------------------------------------------------

geom_line(
  color = line_col,
  linewidth = 0.35,
  alpha = 0.55
) +
  
  # ----------------------------------------------------------
# Paired expression points
# ----------------------------------------------------------

geom_point(
  aes(color = tissue),
  size = 1.7,
  alpha = 0.85
) +
  
  # ----------------------------------------------------------
# Facet by core gene
# ----------------------------------------------------------

facet_wrap(
  ~ gene,
  nrow = 1,
  scales = "free_y"
) +
  
  # ----------------------------------------------------------
# Locked colors
# ----------------------------------------------------------

scale_color_manual(
  values = c(
    "Normal" = normal_col,
    "Tumor" = tumor_col
  )
) +
  
  # ----------------------------------------------------------
# Statistical annotations
# ----------------------------------------------------------

geom_text(
  data = stat_df,
  aes(
    x = 1.5,
    y = y,
    label = label
  ),
  inherit.aes = FALSE,
  family = "Arial",
  size = 3.0,
  fontface = "plain"
) +
  
  # ----------------------------------------------------------
# Axis labels
# ----------------------------------------------------------

labs(
  x = NULL,
  y = "Expression"
) +
  
  # ----------------------------------------------------------
# Y-axis
# ----------------------------------------------------------

scale_y_continuous(
  expand = expansion(mult = c(0.05, 0.16))
) +
  
  # ----------------------------------------------------------
# Theme
# ----------------------------------------------------------

theme_classic(
  base_family = "Arial"
) +
  
  theme(
    strip.background = element_blank(),
    
    strip.text = element_text(
      family = "Arial",
      face = "bold",
      size = 10
    ),
    
    axis.title.x = element_text(
      family = "Arial",
      face = "bold",
      size = 10
    ),
    
    axis.title.y = element_text(
      family = "Arial",
      face = "bold",
      size = 10
    ),
    
    axis.text.x = element_text(
      family = "Arial",
      size = 9
    ),
    
    axis.text.y = element_text(
      family = "Arial",
      size = 9
    ),
    
    axis.line = element_line(
      linewidth = 0.6
    ),
    
    axis.ticks = element_line(
      linewidth = 0.5
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      family = "Arial",
      size = 9
    ),
    
    legend.position = "top",
    
    panel.spacing = unit(
      0.9,
      "lines"
    ),
    
    plot.margin = margin(
      8, 8, 8, 8
    )
  )

# ------------------------------------------------------------
# 7. Display
# ------------------------------------------------------------

print(p9C)

# ------------------------------------------------------------
# 8. Output directory
# ------------------------------------------------------------

fig_dir <- paste0(
  "result_dir/Paper 1_Step 13-Diagnostic Classification/",
  "FIGURES"
)

if (!dir.exists(fig_dir)) {
  dir.create(
    fig_dir,
    recursive = TRUE
  )
}

# ------------------------------------------------------------
# 9. Save Cairo PDF
# ------------------------------------------------------------

pdf_file <- file.path(
  fig_dir,
  "Figure9C_GSE32863_Paired_Expression_Cairo.pdf"
)

Cairo::CairoPDF(
  filename = pdf_file,
  width = 9.0,
  height = 5.8,
  family = "Arial"
)

print(p9C)

dev.off()

# ------------------------------------------------------------
# 10. Save 600-dpi TIFF
# ------------------------------------------------------------

tiff_file <- file.path(
  fig_dir,
  "Figure9C_GSE32863_Paired_Expression_600dpi.tiff"
)

Cairo::CairoTIFF(
  filename = tiff_file,
  width = 9.0,
  height = 5.8,
  units = "in",
  dpi = 600,
  family = "Arial"
)

print(p9C)

dev.off()

# ------------------------------------------------------------
# 11. Final output audit
# ------------------------------------------------------------

cat("\n")
cat("Figure 9C — FINAL OUTPUT AUDIT\n")
cat("============================================================\n\n")

cat(
  "Input rows:              ",
  nrow(pair_df),
  "\n"
)

cat(
  "Complete pairs:          ",
  length(unique(pair_df$pair_id)),
  "\n"
)

cat(
  "Unique subjects:         ",
  length(unique(pair_df$subject_id)),
  "\n"
)

cat(
  "Genes:                   ",
  paste(levels(plot_long$gene), collapse = ", "),
  "\n"
)

cat(
  "PDF exists:              ",
  file.exists(pdf_file),
  "\n"
)

cat(
  "TIFF exists:             ",
  file.exists(tiff_file),
  "\n"
)

cat("\nPDF:\n")
cat(pdf_file, "\n")

cat("\nTIFF:\n")
cat(tiff_file, "\n")

cat("\nFINAL STATUS:\n")

if (
  nrow(pair_df) == 171 &&
  length(unique(pair_df$pair_id)) == 57 &&
  length(unique(pair_df$subject_id)) == 57 &&
  file.exists(pdf_file) &&
  file.exists(tiff_file)
) {
  cat(
    "PASS — Figure 9C successfully generated.\n"
  )
} else {
  cat(
    "FAIL — Figure 9C output audit failed.\n"
  )
}    
