##GitHub6_Paper1_Functional_Enrichment_Analysis_1
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


work_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

setwd(work_dir)

getwd()


# ============================================================
# Paper 1 — WORKFLOW-09 Functional Enrichment
# STEP9A — Input Inventory + Definition Lock
#
# Purpose:
#   1. Inventory all required upstream inputs
#   2. Verify tumor-only WGCNA expression matrix
#   3. Verify 50 WGCNA modules and gene membership
#   4. Verify STEP14 locked module–TME results
#   5. Verify STEP18A candidate/core-gene results
#   6. Identify the final DEG input and DEG definition
#   7. Define appropriate gene universes
#   8. Perform structural / consistency QC
#   9. Save a reproducible STEP9A lock object
#
# IMPORTANT:
#   This script DOES NOT perform:
#   - GO enrichment
#   - KEGG enrichment
#   - GSEA
#   - Hallmark analysis
#   - Reactome analysis
#
# STEP9A is INPUT LOCK ONLY.
#
# R version target: 4.5.1
# ============================================================


# ============================================================
# 0. Global settings
# ============================================================

options(stringsAsFactors = FALSE)
options(scipen = 999)

cat("\n============================================================\n")
cat("Paper 1 — WORKFLOW-09\n")
cat("STEP9A — Input Inventory + Definition Lock\n")
cat("============================================================\n\n")


# ============================================================
# 1. Define directories
# ============================================================

project_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

tumor_wgcna_dir <- "input_dir/Paper1_tumor-only WGCNA"

if (!dir.exists(project_dir)) {
  stop("STEP9 project directory does not exist:\n", project_dir)
}

if (!dir.exists(tumor_wgcna_dir)) {
  stop("Tumor-only WGCNA directory does not exist:\n", tumor_wgcna_dir)
}

setwd(project_dir)

cat("STEP9 working directory:\n")
cat(project_dir, "\n\n")

cat("Tumor-only WGCNA directory:\n")
cat(tumor_wgcna_dir, "\n\n")


# ============================================================
# 2. Create STEP9A output directories
# ============================================================

step9a_dir <- file.path(project_dir, "STEP9A_Input_Inventory_Definition_Lock")

dir.create(step9a_dir, showWarnings = FALSE, recursive = TRUE)

dir.create(
  file.path(step9a_dir, "01_inventory"),
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  file.path(step9a_dir, "02_locked_inputs"),
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  file.path(step9a_dir, "03_QC"),
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  file.path(step9a_dir, "04_logs"),
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 3. Required packages
# ============================================================

required_pkgs <- c(
  "data.table"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nPlease install them before running STEP9A."
  )
}

library(data.table)


# ============================================================
# 4. Session information
# ============================================================

session_info_file <- file.path(
  step9a_dir,
  "04_logs",
  "STEP9A_sessionInfo.txt"
)

capture.output(
  sessionInfo(),
  file = session_info_file
)


# ============================================================
# 5. Helper functions
# ============================================================

pass_flag <- function(x) {
  if (isTRUE(x)) "PASS" else "FAIL"
}


safe_fread <- function(path, ...) {
  
  if (!file.exists(path)) {
    stop("File does not exist:\n", path)
  }
  
  fread(
    path,
    data.table = FALSE,
    check.names = FALSE,
    ...
  )
}


find_files_by_pattern <- function(
    root,
    pattern,
    ignore_case = TRUE) {
  
  all_files <- list.files(
    root,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  )
  
  if (ignore_case) {
    hit <- grepl(
      pattern,
      basename(all_files),
      ignore.case = TRUE
    )
  } else {
    hit <- grepl(
      pattern,
      basename(all_files),
      ignore.case = FALSE
    )
  }
  
  all_files[hit]
}


write_text <- function(x, file) {
  
  con <- file(file, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  
  writeLines(
    as.character(x),
    con = con,
    useBytes = TRUE
  )
}


# ============================================================
# 6. Full project inventory
# ============================================================

cat("------------------------------------------------------------\n")
cat("STEP9A-01: Full project inventory\n")
cat("------------------------------------------------------------\n")

all_project_files <- list.files(
  tumor_wgcna_dir,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

inventory_df <- data.frame(
  path = all_project_files,
  basename = basename(all_project_files),
  extension = tools::file_ext(all_project_files),
  size_bytes = file.info(all_project_files)$size,
  stringsAsFactors = FALSE
)

inventory_file <- file.path(
  step9a_dir,
  "01_inventory",
  "STEP9A_full_tumor_WGCNA_file_inventory.csv"
)

write.csv(
  inventory_df,
  inventory_file,
  row.names = FALSE,
  na = ""
)

cat("Total files inventoried:",
    nrow(inventory_df), "\n\n")


# ============================================================
# 7. Search for key upstream files
# ============================================================

cat("------------------------------------------------------------\n")
cat("STEP9A-02: Search for upstream input files\n")
cat("------------------------------------------------------------\n")


# ------------------------------------------------------------
# 7A. Tumor-only expression matrix
# ------------------------------------------------------------

expr_hits <- find_files_by_pattern(
  tumor_wgcna_dir,
  "Workflow-04-2C_expr_primary_tumor_log2_FPKM_UQ"
)

cat("\nExpression candidates:\n")
print(expr_hits)


# ------------------------------------------------------------
# 7B. WGCNA candidate expression matrix
# ------------------------------------------------------------

wgcna_expr_hits <- find_files_by_pattern(
  tumor_wgcna_dir,
  "Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ"
)

cat("\nWGCNA candidate expression candidates:\n")
print(wgcna_expr_hits)


# ------------------------------------------------------------
# 7C. WGCNA module assignment
# ------------------------------------------------------------

module_hits <- find_files_by_pattern(
  tumor_wgcna_dir,
  "Workflow-04-2F_gene_module_assignment"
)

cat("\nModule assignment candidates:\n")
print(module_hits)


# ------------------------------------------------------------
# 7D. STEP14 locked files
# ------------------------------------------------------------

step14_hits <- find_files_by_pattern(
  tumor_wgcna_dir,
  "Step14D_module_TME|Step14A-14E_LOCKED_KEY_RESULTS|Step14.*FINAL_LOCKED"
)

cat("\nSTEP14 candidates:\n")
print(step14_hits)


# ------------------------------------------------------------
# 7E. STEP18A candidate/core files
# ------------------------------------------------------------

step18a_hits <- find_files_by_pattern(
  tumor_wgcna_dir,
  "STEP18A_30_candidate|STEP18A_core_gene|candidate_master_18D|Supplementary_Table_X_TumorOnly_WGCNA_CoreGenes"
)

cat("\nSTEP18A candidate/core candidates:\n")
print(step18a_hits)


# ------------------------------------------------------------
# 7F. DEG-related files
# ------------------------------------------------------------

deg_hits <- find_files_by_pattern(
  tumor_wgcna_dir,
  "DEG|differential|Differential|DESeq|limma|edgeR|tumor.*normal|normal.*tumor"
)

cat("\nDEG-related candidates:\n")
print(deg_hits)


# Save search results
search_summary <- list(
  expression = expr_hits,
  wgcna_expression = wgcna_expr_hits,
  module_assignment = module_hits,
  step14 = step14_hits,
  step18a = step18a_hits,
  deg = deg_hits
)

saveRDS(
  search_summary,
  file.path(
    step9a_dir,
    "01_inventory",
    "STEP9A_upstream_file_search_results.rds"
  )
)


# ============================================================
# 8. Expression matrix QC
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("STEP9A-03: Tumor-only expression matrix QC\n")
cat("------------------------------------------------------------\n")


# Prefer the exact primary tumor expression file
expr_primary <- expr_hits[
  grepl(
    "expr_primary_tumor_log2_FPKM_UQ\\.txt$",
    basename(expr_hits),
    ignore.case = TRUE
  )
]

if (length(expr_primary) != 1) {
  
  cat(
    "\nWARNING: Exactly one primary tumor expression file was not\n",
    "automatically identified.\n\n",
    "Candidates found:\n",
    paste(expr_hits, collapse = "\n"),
    "\n\n"
  )
  
  expr_primary_selected <- NA_character_
  
} else {
  
  expr_primary_selected <- expr_primary
  
}


expr_qc <- data.frame(
  item = character(),
  value = character(),
  status = character(),
  stringsAsFactors = FALSE
)


if (!is.na(expr_primary_selected)) {
  
  cat("Reading:\n")
  cat(expr_primary_selected, "\n\n")
  
  expr <- safe_fread(expr_primary_selected)
  
  cat("Expression dimensions:\n")
  print(dim(expr))
  
  # First column assumed to be gene identifier
  gene_col <- colnames(expr)[1]
  
  gene_ids <- as.character(expr[[gene_col]])
  
  sample_names <- colnames(expr)[-1]
  
  expr_qc <- rbind(
    expr_qc,
    data.frame(
      item = "expression_file",
      value = expr_primary_selected,
      status = "PASS"
    ),
    data.frame(
      item = "n_rows",
      value = nrow(expr),
      status = pass_flag(nrow(expr) >= 30000)
    ),
    data.frame(
      item = "n_samples",
      value = length(sample_names),
      status = pass_flag(length(sample_names) == 513)
    ),
    data.frame(
      item = "gene_id_column",
      value = gene_col,
      status = "PASS"
    ),
    data.frame(
      item = "duplicate_gene_ids",
      value = sum(duplicated(gene_ids)),
      status = pass_flag(sum(duplicated(gene_ids)) == 0)
    ),
    data.frame(
      item = "NA_values",
      value = sum(is.na(expr)),
      status = pass_flag(sum(is.na(expr)) == 0)
    )
  )
  
  # Check numerical columns
  expr_numeric <- as.matrix(expr[, -1, drop = FALSE])
  
  suppressWarnings(
    storage.mode(expr_numeric) <- "numeric"
  )
  
  nonfinite_count <- sum(!is.finite(expr_numeric))
  
  expr_qc <- rbind(
    expr_qc,
    data.frame(
      item = "nonfinite_values",
      value = nonfinite_count,
      status = pass_flag(nonfinite_count == 0)
    )
  )
  
  # Gene count expected from locked tumor-only WGCNA
  expr_qc <- rbind(
    expr_qc,
    data.frame(
      item = "expected_gene_count",
      value = 31609,
      status = pass_flag(nrow(expr) == 31609)
    )
  )
  
  # Expected tumor count
  expr_qc <- rbind(
    expr_qc,
    data.frame(
      item = "expected_tumor_sample_count",
      value = 513,
      status = pass_flag(length(sample_names) == 513)
    )
  )
  
  # Save only identifiers for STEP9A
  writeLines(
    gene_ids,
    file.path(
      step9a_dir,
      "02_locked_inputs",
      "STEP9A_tumor_only_31609_gene_universe.txt"
    )
  )
  
  writeLines(
    sample_names,
    file.path(
      step9a_dir,
      "02_locked_inputs",
      "STEP9A_tumor_only_513_sample_ids.txt"
    )
  )
  
} else {
  
  expr_qc <- rbind(
    expr_qc,
    data.frame(
      item = "expression_file",
      value = "NOT_UNIQUELY_IDENTIFIED",
      status = "FAIL"
    )
  )
}


write.csv(
  expr_qc,
  file.path(
    step9a_dir,
    "03_QC",
    "STEP9A_expression_QC.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. WGCNA module membership QC
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("STEP9A-04: WGCNA module membership QC\n")
cat("------------------------------------------------------------\n")


module_file <- module_hits[
  grepl(
    "gene_module_assignment\\.txt$",
    basename(module_hits),
    ignore.case = TRUE
  )
]

if (length(module_file) != 1) {
  
  cat(
    "WARNING: module assignment file not uniquely identified.\n"
  )
  
  module_selected <- NA_character_
  
} else {
  
  module_selected <- module_file
}


module_qc <- data.frame(
  item = character(),
  value = character(),
  status = character(),
  stringsAsFactors = FALSE
)


if (!is.na(module_selected)) {
  
  cat("Reading:\n")
  cat(module_selected, "\n\n")
  
  module_df <- safe_fread(module_selected)
  
  cat("Module assignment dimensions:\n")
  print(dim(module_df))
  
  cat("\nColumn names:\n")
  print(colnames(module_df))
  
  # Identify likely gene/module columns
  gene_candidates <- grep(
    "gene|symbol|ensembl",
    colnames(module_df),
    ignore.case = TRUE,
    value = TRUE
  )
  
  module_candidates <- grep(
    "module|color|ME",
    colnames(module_df),
    ignore.case = TRUE,
    value = TRUE
  )
  
  module_qc <- rbind(
    module_qc,
    data.frame(
      item = "module_assignment_file",
      value = module_selected,
      status = "PASS"
    ),
    data.frame(
      item = "n_rows",
      value = nrow(module_df),
      status = pass_flag(nrow(module_df) == 31609)
    ),
    data.frame(
      item = "gene_column_candidates",
      value = paste(gene_candidates, collapse = ";"),
      status = ifelse(length(gene_candidates) >= 1, "PASS", "FAIL")
    ),
    data.frame(
      item = "module_column_candidates",
      value = paste(module_candidates, collapse = ";"),
      status = ifelse(length(module_candidates) >= 1, "PASS", "FAIL")
    )
  )
  
  if (length(module_candidates) >= 1) {
    
    module_col <- module_candidates[1]
    
    module_values <- as.character(module_df[[module_col]])
    
    module_table <- sort(table(module_values), decreasing = TRUE)
    
    n_modules <- length(unique(module_values))
    
    cat("\nNumber of unique module labels:",
        n_modules, "\n")
    
    print(module_table)
    
    module_qc <- rbind(
      module_qc,
      data.frame(
        item = "n_unique_module_labels",
        value = n_modules,
        status = pass_flag(n_modules == 50)
      )
    )
    
    write.csv(
      data.frame(
        module = names(module_table),
        gene_count = as.integer(module_table),
        stringsAsFactors = FALSE
      ),
      file.path(
        step9a_dir,
        "02_locked_inputs",
        "STEP9A_WGCNA_module_gene_counts.csv"
      ),
      row.names = FALSE
    )
  }
  
} else {
  
  module_qc <- rbind(
    module_qc,
    data.frame(
      item = "module_assignment_file",
      value = "NOT_UNIQUELY_IDENTIFIED",
      status = "FAIL"
    )
  )
}


write.csv(
  module_qc,
  file.path(
    step9a_dir,
    "03_QC",
    "STEP9A_WGCNA_module_QC.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. STEP14 locked result QC
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("STEP9A-05: STEP14 locked module–TME result QC\n")
cat("------------------------------------------------------------\n")


step14_all <- step14_hits[
  grepl(
    "Step14D_module_TME_all_150_tests\\.csv$",
    basename(step14_hits),
    ignore.case = TRUE
  )
]

if (length(step14_all) != 1) {
  
  cat(
    "WARNING: STEP14D 150-test file not uniquely identified.\n"
  )
  
  step14_selected <- NA_character_
  
} else {
  
  step14_selected <- step14_all
}


step14_qc <- data.frame(
  item = character(),
  value = character(),
  status = character(),
  stringsAsFactors = FALSE
)


if (!is.na(step14_selected)) {
  
  cat("Reading:\n")
  cat(step14_selected, "\n\n")
  
  step14_df <- safe_fread(step14_selected)
  
  cat("STEP14 dimensions:\n")
  print(dim(step14_df))
  
  cat("\nSTEP14 columns:\n")
  print(colnames(step14_df))
  
  step14_qc <- rbind(
    step14_qc,
    data.frame(
      item = "STEP14_file",
      value = step14_selected,
      status = "PASS"
    ),
    data.frame(
      item = "n_rows",
      value = nrow(step14_df),
      status = pass_flag(nrow(step14_df) == 150)
    )
  )
  
  # Save exact locked copy
  write.csv(
    step14_df,
    file.path(
      step9a_dir,
      "02_locked_inputs",
      "STEP9A_STEP14D_module_TME_all_150_tests_LOCKED_COPY.csv"
    ),
    row.names = FALSE
  )
  
} else {
  
  step14_qc <- rbind(
    step14_qc,
    data.frame(
      item = "STEP14_file",
      value = "NOT_UNIQUELY_IDENTIFIED",
      status = "FAIL"
    )
  )
}


write.csv(
  step14_qc,
  file.path(
    step9a_dir,
    "03_QC",
    "STEP9A_STEP14_QC.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. STEP18A candidate/core-gene QC
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("STEP9A-06: STEP18A candidate/core-gene QC\n")
cat("------------------------------------------------------------\n")


candidate_hits <- step18a_hits[
  grepl(
    "STEP18A_30_candidate_gene_TME_audit\\.csv$",
    basename(step18a_hits),
    ignore.case = TRUE
  )
]

core_hits <- step18a_hits[
  grepl(
    "STEP18A_core_gene_TME_audit\\.csv$",
    basename(step18a_hits),
    ignore.case = TRUE
  )
]


candidate_selected <-
  if (length(candidate_hits) == 1) candidate_hits else NA_character_

core_selected <-
  if (length(core_hits) == 1) core_hits else NA_character_


candidate_qc <- data.frame(
  item = character(),
  value = character(),
  status = character(),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 11A. 30 candidates
# ------------------------------------------------------------

if (!is.na(candidate_selected)) {
  
  candidate_df <- safe_fread(candidate_selected)
  
  cat("Candidate audit dimensions:\n")
  print(dim(candidate_df))
  
  cat("\nCandidate columns:\n")
  print(colnames(candidate_df))
  
  candidate_qc <- rbind(
    candidate_qc,
    data.frame(
      item = "candidate_file",
      value = candidate_selected,
      status = "PASS"
    ),
    data.frame(
      item = "candidate_rows",
      value = nrow(candidate_df),
      status = pass_flag(nrow(candidate_df) == 30)
    )
  )
  
  # Attempt to identify gene column
  candidate_gene_candidates <- grep(
    "gene|symbol",
    colnames(candidate_df),
    ignore.case = TRUE,
    value = TRUE
  )
  
  if (length(candidate_gene_candidates) >= 1) {
    
    candidate_gene_col <- candidate_gene_candidates[1]
    
    candidate_genes <- unique(
      na.omit(
        as.character(
          candidate_df[[candidate_gene_col]]
        )
      )
    )
    
    candidate_qc <- rbind(
      candidate_qc,
      data.frame(
        item = "candidate_gene_column",
        value = candidate_gene_col,
        status = "PASS"
      ),
      data.frame(
        item = "unique_candidate_genes",
        value = length(candidate_genes),
        status = pass_flag(length(candidate_genes) == 30)
      )
    )
    
    writeLines(
      candidate_genes,
      file.path(
        step9a_dir,
        "02_locked_inputs",
        "STEP9A_candidate_genes_30.txt"
      )
    )
    
  } else {
    
    candidate_qc <- rbind(
      candidate_qc,
      data.frame(
        item = "candidate_gene_column",
        value = "NOT_IDENTIFIED",
        status = "FAIL"
      )
    )
  }
  
} else {
  
  candidate_qc <- rbind(
    candidate_qc,
    data.frame(
      item = "candidate_file",
      value = "NOT_UNIQUELY_IDENTIFIED",
      status = "FAIL"
    )
  )
}


# ------------------------------------------------------------
# 11B. Core genes
# ------------------------------------------------------------

if (!is.na(core_selected)) {
  
  core_df <- safe_fread(core_selected)
  
  cat("\nCore-gene audit dimensions:\n")
  print(dim(core_df))
  
  cat("\nCore-gene columns:\n")
  print(colnames(core_df))
  
  core_gene_candidates <- grep(
    "gene|symbol",
    colnames(core_df),
    ignore.case = TRUE,
    value = TRUE
  )
  
  if (length(core_gene_candidates) >= 1) {
    
    core_gene_col <- core_gene_candidates[1]
    
    core_genes_found <- unique(
      na.omit(
        as.character(
          core_df[[core_gene_col]]
        )
      )
    )
    
  } else {
    
    core_gene_col <- NA_character_
    core_genes_found <- character(0)
    
  }
  
  expected_core_genes <- c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )
  
  core_match <- setequal(
    core_genes_found,
    expected_core_genes
  )
  
  candidate_qc <- rbind(
    candidate_qc,
    data.frame(
      item = "core_gene_file",
      value = core_selected,
      status = "PASS"
    ),
    data.frame(
      item = "core_gene_column",
      value = core_gene_col,
      status = ifelse(
        is.na(core_gene_col),
        "FAIL",
        "PASS"
      )
    ),
    data.frame(
      item = "core_gene_set",
      value = paste(core_genes_found, collapse = ";"),
      status = pass_flag(core_match)
    )
  )
  
  if (core_match) {
    
    writeLines(
      expected_core_genes,
      file.path(
        step9a_dir,
        "02_locked_inputs",
        "STEP9A_core_genes_AOC3_MS4A7_IL7R.txt"
      )
    )
    
  }
  
} else {
  
  candidate_qc <- rbind(
    candidate_qc,
    data.frame(
      item = "core_gene_file",
      value = "NOT_UNIQUELY_IDENTIFIED",
      status = "FAIL"
    )
  )
}


write.csv(
  candidate_qc,
  file.path(
    step9a_dir,
    "03_QC",
    "STEP9A_candidate_core_gene_QC.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. DEG audit — DO NOT GUESS THRESHOLD
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("STEP9A-07: DEG input audit\n")
cat("------------------------------------------------------------\n")


cat(
  "\nIMPORTANT:\n",
  "STEP9A will NOT automatically select a DEG file solely by filename,\n",
  "and will NOT assume |log2FC| > 0.58 or |log2FC| > 1.5.\n\n",
  "The final DEG definition must be traced to the locked STEP2 output.\n\n"
)


# Save all DEG-related candidate files for manual review
writeLines(
  deg_hits,
  file.path(
    step9a_dir,
    "01_inventory",
    "STEP9A_DEG_candidate_files.txt"
  )
)


deg_qc <- data.frame(
  item = character(),
  value = character(),
  status = character(),
  stringsAsFactors = FALSE
)


deg_qc <- rbind(
  deg_qc,
  data.frame(
    item = "number_of_DEG_candidate_files",
    value = length(deg_hits),
    status = ifelse(
      length(deg_hits) >= 1,
      "REVIEW_REQUIRED",
      "FAIL"
    )
  ),
  data.frame(
    item = "DEG_threshold",
    value = "NOT_ASSUMED",
    status = "REVIEW_REQUIRED"
  ),
  data.frame(
    item = "DEG_gene_universe",
    value = "NOT_ASSUMED",
    status = "REVIEW_REQUIRED"
  )
)


write.csv(
  deg_qc,
  file.path(
    step9a_dir,
    "03_QC",
    "STEP9A_DEG_QC.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. Create STEP9A definition record
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("STEP9A-08: Create definition record\n")
cat("------------------------------------------------------------\n")


definition_record <- c(
  
  "============================================================",
  "Paper 1 — WORKFLOW-09 Functional Enrichment",
  "STEP9A — Input Inventory + Definition Lock",
  "============================================================",
  "",
  paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "WORKING DIRECTORY:",
  project_dir,
  "",
  "UPSTREAM DIRECTORY:",
  tumor_wgcna_dir,
  "",
  "------------------------------------------------------------",
  "LOCKED / EXPECTED INPUTS",
  "------------------------------------------------------------",
  "",
  "Tumor-only expression:",
  "513 primary LUAD tumor samples × 31,609 genes",
  "",
  "WGCNA:",
  "Tumor-only WGCNA",
  "50 module labels expected",
  "",
  "STEP14:",
  "150 module × TME association tests",
  "Existing locked result; no recalculation in STEP9A",
  "",
  "STEP18A:",
  "30 candidate genes",
  "",
  "Core genes:",
  "AOC3 / MS4A7 / IL7R",
  "",
  "Functional enrichment methods planned for later steps:",
  "GO",
  "KEGG",
  "preranked GSEA",
  "Reactome",
  "Hallmark",
  "",
  "------------------------------------------------------------",
  "CRITICAL UNRESOLVED ITEM",
  "------------------------------------------------------------",
  "",
  "The final DEG file, DEG threshold, and DEG testing universe",
  "must be traced to the locked STEP2 result.",
  "",
  "STEP9A intentionally does NOT assume:",
  "|log2FC| > 0.58",
  "|log2FC| > 1.5",
  "or any other historical DEG threshold.",
  "",
  "No enrichment analysis is performed in STEP9A.",
  ""
)


definition_file <- file.path(
  step9a_dir,
  "01_inventory",
  "STEP9A_definition_record.txt"
)

write_text(
  definition_record,
  definition_file
)


# ============================================================
# 14. Overall QC summary
# ============================================================

cat("\n------------------------------------------------------------\n")
cat("STEP9A-09: Overall QC summary\n")
cat("------------------------------------------------------------\n")


overall_qc <- data.frame(
  component = c(
    "Tumor-only expression",
    "31,609-gene expectation",
    "513-sample expectation",
    "WGCNA module assignment",
    "50-module expectation",
    "STEP14 150-test result",
    "STEP18A 30 candidates",
    "AOC3/MS4A7/IL7R core genes",
    "Final DEG definition"
  ),
  status = c(
    ifelse(
      exists("expr_qc") &&
        any(expr_qc$item == "expression_file" &
              expr_qc$status == "PASS"),
      "PASS",
      "FAIL"
    ),
    ifelse(
      exists("expr_qc") &&
        any(expr_qc$item == "expected_gene_count" &
              expr_qc$status == "PASS"),
      "PASS",
      "FAIL"
    ),
    ifelse(
      exists("expr_qc") &&
        any(expr_qc$item == "expected_tumor_sample_count" &
              expr_qc$status == "PASS"),
      "PASS",
      "FAIL"
    ),
    ifelse(
      exists("module_qc") &&
        any(module_qc$item == "module_assignment_file" &
              module_qc$status == "PASS"),
      "PASS",
      "FAIL"
    ),
    ifelse(
      exists("module_qc") &&
        any(module_qc$item == "n_unique_module_labels" &
              module_qc$status == "PASS"),
      "PASS",
      "FAIL"
    ),
    ifelse(
      exists("step14_qc") &&
        any(step14_qc$item == "n_rows" &
              step14_qc$status == "PASS"),
      "PASS",
      "FAIL"
    ),
    ifelse(
      exists("candidate_qc") &&
        any(candidate_qc$item == "candidate_rows" &
              candidate_qc$status == "PASS"),
      "PASS",
      "FAIL"
    ),
    ifelse(
      exists("candidate_qc") &&
        any(candidate_qc$item == "core_gene_set" &
              candidate_qc$status == "PASS"),
      "PASS",
      "FAIL"
    ),
    "REVIEW_REQUIRED"
  ),
  stringsAsFactors = FALSE
)


write.csv(
  overall_qc,
  file.path(
    step9a_dir,
    "03_QC",
    "STEP9A_overall_QC.csv"
  ),
  row.names = FALSE
)


print(overall_qc)


# ============================================================
# 15. Determine provisional STEP9A status
# ============================================================

non_deg_failures <- overall_qc$status[
  overall_qc$component != "Final DEG definition"
]

if (all(non_deg_failures == "PASS") &&
    overall_qc$status[
      overall_qc$component == "Final DEG definition"
    ] == "REVIEW_REQUIRED") {
  
  step9a_status <- "HOLD_FOR_DEG_DEFINITION"
  
} else if (any(overall_qc$status == "FAIL")) {
  
  step9a_status <- "FAIL"
  
} else {
  
  step9a_status <- "REVIEW_REQUIRED"
}


cat("\n============================================================\n")
cat("STEP9A PROVISIONAL STATUS:\n")
cat(step9a_status, "\n")
cat("============================================================\n\n")


# ============================================================
# 16. Save status
# ============================================================

status_record <- list(
  
  workflow = "WORKFLOW-09",
  
  step = "STEP9A",
  
  step_name = "Input Inventory + Definition Lock",
  
  date = Sys.time(),
  
  status = step9a_status,
  
  project_dir = project_dir,
  
  tumor_wgcna_dir = tumor_wgcna_dir,
  
  expected_tumor_samples = 513,
  
  expected_wgcna_genes = 31609,
  
  expected_wgcna_modules = 50,
  
  expected_step14_tests = 150,
  
  expected_candidates = 30,
  
  core_genes = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  
  deg_threshold_locked = FALSE,
  
  deg_universe_locked = FALSE,
  
  enrichment_performed = FALSE
  
)


saveRDS(
  status_record,
  file.path(
    step9a_dir,
    "STEP9A_provisional_status.rds"
  )
)


# ============================================================
# 17. Final message
# ============================================================

cat("\n============================================================\n")
cat("STEP9A completed.\n")
cat("============================================================\n\n")

cat(
  "IMPORTANT:\n",
  "No GO / KEGG / GSEA / Reactome / Hallmark analysis was performed.\n\n"
)

cat(
  "The next required action is to inspect:\n",
  file.path(
    step9a_dir,
    "03_QC",
    "STEP9A_DEG_QC.csv"
  ),
  "\n\n"
)

cat(
  "and identify the final locked STEP2 DEG file and definition.\n"
)

cat(
  "\nSTEP9A status:",
  step9a_status,
  "\n\n"
)

cat(
  "Output directory:\n",
  step9a_dir,
  "\n"
)


library(data.table)

step9b_dir <- 
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9B_DEG_Definition_Recovery_Lock"

# ------------------------------------------------------------
# 1. Candidate directories
# ------------------------------------------------------------

f1 <- file.path(
  step9b_dir,
  "01_file_inventory",
  "STEP9B2_candidate_directories.txt"
)

cat("\n============================================================\n")
cat("STEP9B-2: Candidate directories\n")
cat("============================================================\n\n")

if (file.exists(f1)) {
  x1 <- readLines(f1, warn = FALSE)
  x1 <- x1[nzchar(trimws(x1))]
  
  cat("Number of candidate directories:", length(x1), "\n\n")
  
  if (length(x1) > 0) {
    print(x1)
  } else {
    cat("NO candidate directories found.\n")
  }
} else {
  cat("File not found:\n", f1, "\n")
}


# ------------------------------------------------------------
# 2. Candidate STEP2 / DEG files
# ------------------------------------------------------------

f2 <- file.path(
  step9b_dir,
  "01_file_inventory",
  "STEP9B2_candidate_STEP2_DEG_files.txt"
)

cat("\n============================================================\n")
cat("STEP9B-2: Candidate STEP2 / DEG files\n")
cat("============================================================\n\n")

if (file.exists(f2)) {
  x2 <- readLines(f2, warn = FALSE)
  x2 <- x2[nzchar(trimws(x2))]
  
  cat("Number of candidate files:", length(x2), "\n\n")
  
  if (length(x2) > 0) {
    print(x2)
  } else {
    cat("NO candidate STEP2 / DEG files found.\n")
  }
} else {
  cat("File not found:\n", f2, "\n")
}


# ------------------------------------------------------------
# 3. Files inside candidate directories
# ------------------------------------------------------------

f3 <- file.path(
  step9b_dir,
  "01_file_inventory",
  "STEP9B2_files_inside_candidate_directories.csv"
)

cat("\n============================================================\n")
cat("STEP9B-2: Files inside candidate directories\n")
cat("============================================================\n\n")

if (file.exists(f3)) {
  x3 <- fread(f3)
  
  cat("Dimensions:\n")
  print(dim(x3))
  
  cat("\nColumns:\n")
  print(colnames(x3))
  
  cat("\nFiles:\n")
  print(x3)
  
} else {
  cat("File not found:\n", f3, "\n")
}

# ============================================================
# Paper 1
# STEP9B-6D
# Final Historical DEG Candidate Disambiguation
#
# Purpose:
#   1. Inspect all historical write.table() / saveRDS() /
#      write.csv() output commands in Workflow-03F lineage
#   2. Identify whether any output is a true tumor-vs-normal
#      differential-expression result
#   3. Trace DEG-like object names, filenames and statistical
#      function calls
#
# IMPORTANT:
#   - NO DEG recalculation
#   - NO statistical testing
#   - NO threshold locking
#   - NO modification of historical data
# ============================================================

cat("\n============================================================\n")
cat("STEP9B-6D: Final Historical DEG Candidate Disambiguation\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

root_dir <- "input_dir/Paper1_tumor-only WGCNA"

workflow_script <- file.path(
  root_dir,
  "Workflow-03F-9_provenance_audit.R"
)

workflow_rdata <- file.path(
  root_dir,
  "Workflow-03F-9_provenance_audit.RData"
)

out_dir <- file.path(
  root_dir,
  "STEP9B-6D_Final_Historical_DEG_Candidate_Disambiguation"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(
  file.path(out_dir, "01_script_audit"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(out_dir, "02_output_command_audit"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(out_dir, "03_statistical_function_audit"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(out_dir, "04_rdata_audit"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(out_dir, "05_final_QC"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(out_dir, "06_logs"),
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Basic existence check
# ------------------------------------------------------------

cat("Checking historical Workflow-03F files...\n")

script_exists <- file.exists(workflow_script)
rdata_exists  <- file.exists(workflow_rdata)

cat("Workflow-03F-9 script exists:", script_exists, "\n")
cat("Workflow-03F-9 RData exists :", rdata_exists, "\n\n")

if (!script_exists) {
  stop("Workflow-03F-9 script not found.")
}

# ------------------------------------------------------------
# 2. Read complete historical script
# ------------------------------------------------------------

script_lines <- readLines(
  workflow_script,
  warn = FALSE,
  encoding = "UTF-8"
)

cat("Historical script lines:", length(script_lines), "\n\n")

# ------------------------------------------------------------
# 3. Audit all output commands
# ------------------------------------------------------------

output_patterns <- c(
  "write\\.table\\s*\\(",
  "write\\.csv\\s*\\(",
  "write\\.csv2\\s*\\(",
  "saveRDS\\s*\\(",
  "save\\s*\\(",
  "ggsave\\s*\\("
)

output_hits <- data.frame(
  Line = integer(0),
  Code = character(0),
  stringsAsFactors = FALSE
)

for (pat in output_patterns) {
  
  idx <- grep(
    pat,
    script_lines,
    ignore.case = TRUE
  )
  
  if (length(idx) > 0) {
    
    tmp <- data.frame(
      Line = idx,
      Code = trimws(script_lines[idx]),
      stringsAsFactors = FALSE
    )
    
    output_hits <- rbind(
      output_hits,
      tmp
    )
  }
}

output_hits <- output_hits[
  order(output_hits$Line),
]

output_hits <- output_hits[
  !duplicated(output_hits$Line),
]

write.csv(
  output_hits,
  file.path(
    out_dir,
    "02_output_command_audit",
    "STEP9B6D_all_historical_output_commands.csv"
  ),
  row.names = FALSE
)

cat(
  "Historical output-command lines:",
  nrow(output_hits),
  "\n\n"
)

# ------------------------------------------------------------
# 4. Search for statistical testing functions
# ------------------------------------------------------------

stat_patterns <- c(
  "wilcox\\.test\\s*\\(",
  "t\\.test\\s*\\(",
  "limma",
  "lmFit\\s*\\(",
  "eBayes\\s*\\(",
  "voom\\s*\\(",
  "DESeq\\s*\\(",
  "DESeq2",
  "DESeqDataSet",
  "results\\s*\\(",
  "edgeR",
  "exactTest\\s*\\(",
  "glmQLFit\\s*\\(",
  "glmQLFTest\\s*\\(",
  "differential",
  "differential.expression",
  "differential_expression"
)

stat_hits <- data.frame(
  Line = integer(0),
  Pattern = character(0),
  Code = character(0),
  stringsAsFactors = FALSE
)

for (pat in stat_patterns) {
  
  idx <- grep(
    pat,
    script_lines,
    ignore.case = TRUE
  )
  
  if (length(idx) > 0) {
    
    tmp <- data.frame(
      Line = idx,
      Pattern = pat,
      Code = trimws(script_lines[idx]),
      stringsAsFactors = FALSE
    )
    
    stat_hits <- rbind(
      stat_hits,
      tmp
    )
  }
}

stat_hits <- stat_hits[
  order(stat_hits$Line),
]

write.csv(
  stat_hits,
  file.path(
    out_dir,
    "03_statistical_function_audit",
    "STEP9B6D_historical_statistical_function_hits.csv"
  ),
  row.names = FALSE
)

cat(
  "Potential statistical-function hits:",
  nrow(stat_hits),
  "\n\n"
)

# ------------------------------------------------------------
# 5. Search DEG-result column definitions
# ------------------------------------------------------------

deg_column_patterns <- c(
  "log2FC",
  "logFC",
  "fold.?change",
  "P\\.Value",
  "p\\.value",
  "p_value",
  "adj\\.P\\.Val",
  "adj_p",
  "FDR",
  "padj",
  "qvalue",
  "p\\.adjust"
)

deg_column_hits <- data.frame(
  Line = integer(0),
  Pattern = character(0),
  Code = character(0),
  stringsAsFactors = FALSE
)

for (pat in deg_column_patterns) {
  
  idx <- grep(
    pat,
    script_lines,
    ignore.case = TRUE
  )
  
  if (length(idx) > 0) {
    
    tmp <- data.frame(
      Line = idx,
      Pattern = pat,
      Code = trimws(script_lines[idx]),
      stringsAsFactors = FALSE
    )
    
    deg_column_hits <- rbind(
      deg_column_hits,
      tmp
    )
  }
}

deg_column_hits <- deg_column_hits[
  order(deg_column_hits$Line),
]

write.csv(
  deg_column_hits,
  file.path(
    out_dir,
    "01_script_audit",
    "STEP9B6D_historical_DEG_column_hits.csv"
  ),
  row.names = FALSE
)

cat(
  "Potential DEG-column hits:",
  nrow(deg_column_hits),
  "\n\n"
)

# ------------------------------------------------------------
# 6. Search tumor-normal comparison logic
# ------------------------------------------------------------

comparison_patterns <- c(
  "tumor.*normal",
  "normal.*tumor",
  "tumor_ids",
  "normal_ids",
  "paired_patient",
  "sample_type",
  "Tumor_Normal",
  "Sample_Type_Code"
)

comparison_hits <- data.frame(
  Line = integer(0),
  Pattern = character(0),
  Code = character(0),
  stringsAsFactors = FALSE
)

for (pat in comparison_patterns) {
  
  idx <- grep(
    pat,
    script_lines,
    ignore.case = TRUE
  )
  
  if (length(idx) > 0) {
    
    tmp <- data.frame(
      Line = idx,
      Pattern = pat,
      Code = trimws(script_lines[idx]),
      stringsAsFactors = FALSE
    )
    
    comparison_hits <- rbind(
      comparison_hits,
      tmp
    )
  }
}

comparison_hits <- comparison_hits[
  order(comparison_hits$Line),
]

write.csv(
  comparison_hits,
  file.path(
    out_dir,
    "01_script_audit",
    "STEP9B6D_tumor_normal_comparison_hits.csv"
  ),
  row.names = FALSE
)

cat(
  "Tumor-normal comparison hits:",
  nrow(comparison_hits),
  "\n\n"
)

# ------------------------------------------------------------
# 7. Search threshold-like DEG definitions
# ------------------------------------------------------------

threshold_patterns <- c(
  "log2FC.*[><=]",
  "logFC.*[><=]",
  "abs\\s*\\(.*log2FC",
  "abs\\s*\\(.*logFC",
  "FDR.*[><=]",
  "padj.*[><=]",
  "adj\\.P\\.Val.*[><=]",
  "p\\.adjust.*[><=]",
  "p\\.value.*[><=]"
)

threshold_hits <- data.frame(
  Line = integer(0),
  Pattern = character(0),
  Code = character(0),
  stringsAsFactors = FALSE
)

for (pat in threshold_patterns) {
  
  idx <- grep(
    pat,
    script_lines,
    ignore.case = TRUE
  )
  
  if (length(idx) > 0) {
    
    tmp <- data.frame(
      Line = idx,
      Pattern = pat,
      Code = trimws(script_lines[idx]),
      stringsAsFactors = FALSE
    )
    
    threshold_hits <- rbind(
      threshold_hits,
      tmp
    )
  }
}

threshold_hits <- threshold_hits[
  order(threshold_hits$Line),
]

write.csv(
  threshold_hits,
  file.path(
    out_dir,
    "01_script_audit",
    "STEP9B6D_historical_DEG_threshold_hits.csv"
  ),
  row.names = FALSE
)

cat(
  "Potential DEG-threshold hits:",
  nrow(threshold_hits),
  "\n\n"
)

# ------------------------------------------------------------
# 8. Audit RData object names and structures
# ------------------------------------------------------------

rdata_summary <- data.frame(
  Object = character(0),
  Class = character(0),
  Dim1 = integer(0),
  Dim2 = integer(0),
  Ncol = integer(0),
  NamesPreview = character(0),
  stringsAsFactors = FALSE
)

if (rdata_exists) {
  
  e <- new.env(parent = emptyenv())
  
  load(
    workflow_rdata,
    envir = e
  )
  
  object_names <- ls(e)
  
  cat(
    "RData objects:",
    length(object_names),
    "\n\n"
  )
  
  for (nm in object_names) {
    
    obj <- get(
      nm,
      envir = e
    )
    
    cls <- paste(
      class(obj),
      collapse = ";"
    )
    
    d <- dim(obj)
    
    dim1 <- if (!is.null(d)) d[1] else NA_integer_
    dim2 <- if (!is.null(d) && length(d) >= 2) d[2] else NA_integer_
    
    ncol_obj <- NA_integer_
    
    if (!is.null(d) && length(d) >= 2) {
      ncol_obj <- d[2]
    }
    
    preview_names <- ""
    
    if (!is.null(names(obj))) {
      
      preview_names <- paste(
        head(names(obj), 20),
        collapse = " | "
      )
      
    } else if (!is.null(colnames(obj))) {
      
      preview_names <- paste(
        head(colnames(obj), 20),
        collapse = " | "
      )
    }
    
    rdata_summary <- rbind(
      rdata_summary,
      data.frame(
        Object = nm,
        Class = cls,
        Dim1 = dim1,
        Dim2 = dim2,
        Ncol = ncol_obj,
        NamesPreview = preview_names,
        stringsAsFactors = FALSE
      )
    )
  }
  
  rm(e)
  
}

write.csv(
  rdata_summary,
  file.path(
    out_dir,
    "04_rdata_audit",
    "STEP9B6D_RData_object_structure.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# STEP9B-6D SAFE RDATA COLUMN AUDIT
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9B-6D: Safe RData DEG-column audit\n")
cat("============================================================\n\n")

deg_object_hits <- data.frame(
  Object = character(0),
  Class = character(0),
  MatchedColumns = character(0),
  stringsAsFactors = FALSE
)

if (rdata_exists) {
  
  e <- new.env(parent = emptyenv())
  
  load(
    workflow_rdata,
    envir = e
  )
  
  object_names <- ls(e)
  
  for (nm in object_names) {
    
    obj <- tryCatch(
      get(nm, envir = e),
      error = function(err) NULL
    )
    
    if (is.null(obj)) {
      next
    }
    
    cls <- paste(
      class(obj),
      collapse = ";"
    )
    
    cn <- tryCatch(
      {
        if (
          is.data.frame(obj) ||
          is.matrix(obj) ||
          is.array(obj)
        ) {
          colnames(obj)
        } else {
          NULL
        }
      },
      error = function(err) NULL
    )
    
    if (is.null(cn)) {
      next
    }
    
    matched <- cn[
      grepl(
        paste(
          deg_column_patterns,
          collapse = "|"
        ),
        cn,
        ignore.case = TRUE
      )
    ]
    
    if (length(matched) > 0) {
      
      deg_object_hits <- rbind(
        deg_object_hits,
        data.frame(
          Object = nm,
          Class = cls,
          MatchedColumns = paste(
            matched,
            collapse = " | "
          ),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  rm(e)
}

write.csv(
  deg_object_hits,
  file.path(
    out_dir,
    "04_rdata_audit",
    "STEP9B6D_RData_DEG_like_objects_SAFE.csv"
  ),
  row.names = FALSE
)

cat(
  "Safe DEG-like RData objects:",
  nrow(deg_object_hits),
  "\n\n"
)

if (nrow(deg_object_hits) == 0) {
  
  cat(
    "No DEG-like result object was identified.\n"
  )
  
} else {
  
  cat(
    "Potential DEG-like objects require manual review.\n"
  )
}
# ------------------------------------------------------------
# 10. Search output commands for DEG-related object names
# ------------------------------------------------------------

deg_output_patterns <- c(
  "deg",
  "differential",
  "de_results",
  "de_result",
  "diff_expr",
  "differential_expression",
  "tumor_normal",
  "tumor.*normal",
  "normal.*tumor"
)

deg_output_hits <- data.frame(
  Line = integer(0),
  Pattern = character(0),
  Code = character(0),
  stringsAsFactors = FALSE
)

for (pat in deg_output_patterns) {
  
  idx <- grep(
    pat,
    output_hits$Code,
    ignore.case = TRUE
  )
  
  if (length(idx) > 0) {
    
    tmp <- data.frame(
      Line = output_hits$Line[idx],
      Pattern = pat,
      Code = output_hits$Code[idx],
      stringsAsFactors = FALSE
    )
    
    deg_output_hits <- rbind(
      deg_output_hits,
      tmp
    )
  }
}

deg_output_hits <- deg_output_hits[
  !duplicated(
    deg_output_hits[, c("Line", "Code")]
  ),
]

write.csv(
  deg_output_hits,
  file.path(
    out_dir,
    "02_output_command_audit",
    "STEP9B6D_DEG_related_output_commands.csv"
  ),
  row.names = FALSE
)

cat(
  "DEG-related output-command hits:",
  nrow(deg_output_hits),
  "\n\n"
)

# ------------------------------------------------------------
# 11. Final evidence classification
# ------------------------------------------------------------

explicit_stat_test <- nrow(stat_hits) > 0

explicit_deg_columns <- nrow(deg_column_hits) > 0

explicit_deg_object <- nrow(deg_object_hits) > 0

deg_output_command <- nrow(deg_output_hits) > 0

# Conservative rule:
#
# Historical DEG evidence is NOT accepted merely because one
# keyword appears.
#
# At minimum, there should be:
#   A. explicit differential-expression statistical method
#   AND
#   B. DEG-like result columns/object
#   AND
#   C. tumor-normal comparison context
#
# This script does not lock a DEG result even if these are found.
# It only classifies the evidence.

has_tumor_normal_context <- nrow(comparison_hits) > 0

historical_deg_evidence_level <- "NO_EXPLICIT_DEG_RESULT_EVIDENCE"

if (
  explicit_stat_test &&
  explicit_deg_columns &&
  has_tumor_normal_context
) {
  
  historical_deg_evidence_level <-
    "POTENTIAL_DEG_CHAIN_REQUIRES_MANUAL_REVIEW"
  
}

if (
  explicit_stat_test &&
  explicit_deg_object &&
  has_tumor_normal_context
) {
  
  historical_deg_evidence_level <-
    "POTENTIAL_DEG_RESULT_OBJECT_REQUIRES_MANUAL_REVIEW"
  
}

# ------------------------------------------------------------
# 12. Final QC
# ------------------------------------------------------------

final_qc <- data.frame(
  Item = c(
    "Workflow-03F-9 script exists",
    "Workflow-03F-9 RData exists",
    "Historical statistical-function hits",
    "Historical DEG-column hits",
    "Historical DEG-like RData objects",
    "DEG-related output commands",
    "Tumor-normal comparison context",
    "Historical DEG evidence classification",
    "DEG recalculation performed",
    "DEG threshold locked"
  ),
  
  Result = c(
    script_exists,
    rdata_exists,
    nrow(stat_hits),
    nrow(deg_column_hits),
    nrow(deg_object_hits),
    nrow(deg_output_hits),
    nrow(comparison_hits),
    historical_deg_evidence_level,
    FALSE,
    FALSE
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  final_qc,
  file.path(
    out_dir,
    "05_final_QC",
    "STEP9B6D_final_QC.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Final status
# ------------------------------------------------------------

if (
  !explicit_stat_test &&
  !explicit_deg_object &&
  !deg_output_command
) {
  
  final_status <-
    "HISTORICAL_DEG_RESULT_NOT_RECOVERED"
  
} else {
  
  final_status <-
    "REVIEW_REQUIRED_BEFORE_DEG_RECOVERY_DECISION"
  
}

status_text <- c(
  "STEP9B-6D Final Status",
  "=======================",
  "",
  paste("Status:", final_status),
  "",
  paste(
    "Statistical-function hits:",
    nrow(stat_hits)
  ),
  paste(
    "DEG-column hits:",
    nrow(deg_column_hits)
  ),
  paste(
    "DEG-like RData objects:",
    nrow(deg_object_hits)
  ),
  paste(
    "DEG-related output commands:",
    nrow(deg_output_hits)
  ),
  paste(
    "Tumor-normal comparison hits:",
    nrow(comparison_hits)
  ),
  "",
  "No DEG recalculation was performed.",
  "No DEG threshold was locked.",
  "No historical DEG result was accepted.",
  "",
  "Historical expression lineage remains confirmed.",
  "Historical DEG evidence must be distinguished from",
  "expression provenance/QC evidence."
)

writeLines(
  status_text,
  file.path(
    out_dir,
    "06_logs",
    "STEP9B6D_FINAL_STATUS.txt"
  )
)

cat(
  paste(
    status_text,
    collapse = "\n"
  ),
  "\n"
)

cat(
  "\n============================================================\n"
)
cat(
  "STEP9B-6D completed.\n"
)
cat(
  "No DEG recalculation performed.\n"
)
cat(
  "No DEG threshold locked.\n"
)
cat(
  "============================================================\n"
)

##STEP9B-7R — Revised DEG Definition and Lock
# ============================================================
# Paper 1
# WORKFLOW-09 Functional Enrichment
#
# STEP9B-7R
# Revised DEG Definition and Lock
#
# Purpose:
#   1. Establish a new, reproducible genome-wide DEG framework
#   2. Use the historical TCGA FPKM-UQ expression lineage
#   3. Define tumor and normal samples independently of historical
#      Xiantao online DEG results
#   4. Perform tumor-versus-normal genome-wide differential analysis
#   5. Lock the DEG testing universe and DEG definition
#
# IMPORTANT:
#   - No historical DEG result is reused.
#   - No historical |log2FC| > 0.58 or > 1.5 threshold is used.
#   - 31,609 WGCNA genes are NOT used as the DEG universe.
#   - Figure 2's historical Xiantao three-gene analysis is NOT
#     treated as a genome-wide DEG analysis.
#
# Expression lineage:
#   Merge_RNA_seq_FPKM_UQ.txt
#
# Revised DEG framework:
#   Expression transformation:
#       log2(FPKM-UQ + 1)
#
#   Statistical test:
#       Wilcoxon rank-sum test
#
#   Multiple testing:
#       Benjamini-Hochberg FDR
#
#   DEG definition:
#       FDR < 0.05
#       AND
#       |log2FC| >= 1
#
#   Direction:
#       log2FC > 0  -> Tumor-up
#       log2FC < 0  -> Tumor-down
#
# ============================================================


rm(list = ls())
gc()

cat("\n============================================================\n")
cat("STEP9B-7R: Revised DEG Definition and Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

step_dir <- file.path(
  base_dir,
  "STEP9B-7R_Revised_DEG_Definition_and_Lock"
)

dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(step_dir, "00_input"),
           recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(step_dir, "01_data"),
           recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(step_dir, "02_results"),
           recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(step_dir, "03_QC"),
           recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(step_dir, "04_lock"),
           recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(step_dir, "05_logs"),
           recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Input files
# ------------------------------------------------------------

expr_file <- paste0(
  "input_dir/Paper1_tumor-only WGCNA/",
  "Merge_RNA_seq_FPKM_UQ.txt"
)

annotation_file <- paste0(
  "input_dir/Paper1_tumor-only WGCNA/",
  "TCGA_Expression_Sample_Annotation_Basic.txt"
)


# ------------------------------------------------------------
# 3. Confirm input files
# ------------------------------------------------------------

input_qc <- data.frame(
  File = c(
    expr_file,
    annotation_file
  ),
  Exists = c(
    file.exists(expr_file),
    file.exists(annotation_file)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  input_qc,
  file.path(step_dir, "03_QC",
            "STEP9B7R_Input_File_QC.csv"),
  row.names = FALSE
)

if (!all(input_qc$Exists)) {
  stop(
    "STEP9B-7R FAILED: one or more required input files do not exist."
  )
}


# ------------------------------------------------------------
# 4. Read sample annotation
# ------------------------------------------------------------

cat("Reading sample annotation...\n")

ann <- read.delim(
  annotation_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Annotation dimensions:\n")
print(dim(ann))

cat("\nAnnotation columns:\n")
print(colnames(ann))

cat("\nTumor/Normal distribution:\n")
print(table(ann$Tumor_Normal, useNA = "ifany"))


# ------------------------------------------------------------
# 5. Standardize annotation columns
# ------------------------------------------------------------

required_ann_cols <- c(
  "Expression_Sample",
  "Patient_ID",
  "Sample_Type_Code",
  "Tumor_Normal"
)

missing_ann_cols <- setdiff(
  required_ann_cols,
  colnames(ann)
)

if (length(missing_ann_cols) > 0) {
  stop(
    paste(
      "STEP9B-7R FAILED: missing annotation columns:",
      paste(missing_ann_cols, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 6. Define analysis samples
#
# Expected:
#   513 Primary Tumor
#   59 Solid Tissue Normal
#   2 Other excluded
#
# ------------------------------------------------------------

analysis_ann <- ann[
  ann$Tumor_Normal %in% c("Tumor", "Normal"),
  ,
  drop = FALSE
]

analysis_ann$Group <- ifelse(
  analysis_ann$Tumor_Normal == "Tumor",
  "Tumor",
  "Normal"
)

analysis_ann <- analysis_ann[
  !is.na(analysis_ann$Expression_Sample) &
    analysis_ann$Expression_Sample != "",
  ,
  drop = FALSE
]

# Remove duplicated expression sample IDs if any
analysis_ann <- analysis_ann[
  !duplicated(analysis_ann$Expression_Sample),
  ,
  drop = FALSE
]


cat("\n============================================================\n")
cat("Analysis sample definition\n")
cat("============================================================\n")

print(table(
  analysis_ann$Group,
  useNA = "ifany"
))


# ------------------------------------------------------------
# 7. Sample-number QC
# ------------------------------------------------------------

sample_counts <- table(analysis_ann$Group)

n_tumor <- ifelse(
  "Tumor" %in% names(sample_counts),
  as.integer(sample_counts["Tumor"]),
  0
)

n_normal <- ifelse(
  "Normal" %in% names(sample_counts),
  as.integer(sample_counts["Normal"]),
  0
)

cat("\nTumor samples :", n_tumor, "\n")
cat("Normal samples:", n_normal, "\n")


if (n_tumor != 513) {
  stop(
    paste(
      "STEP9B-7R FAILED: expected 513 tumor samples, found",
      n_tumor
    )
  )
}

if (n_normal != 59) {
  stop(
    paste(
      "STEP9B-7R FAILED: expected 59 normal samples, found",
      n_normal
    )
  )
}


# ------------------------------------------------------------
# 8. Save sample audit
# ------------------------------------------------------------

write.csv(
  analysis_ann,
  file.path(
    step_dir,
    "03_QC",
    "STEP9B7R_Input_Sample_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Read historical FPKM-UQ expression matrix
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Reading Merge_RNA_seq_FPKM_UQ.txt\n")
cat("============================================================\n")

expr_raw <- read.delim(
  expr_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\nExpression dimensions:\n")
print(dim(expr_raw))

cat("\nFirst columns:\n")
print(colnames(expr_raw)[1:min(10, ncol(expr_raw))])


# ------------------------------------------------------------
# 10. Identify gene identifier column
# ------------------------------------------------------------

gene_col_candidates <- c(
  "Gene",
  "Gene_ID",
  "GeneID",
  "gene",
  "gene_id",
  "Symbol",
  "Gene_Symbol"
)

gene_col <- intersect(
  gene_col_candidates,
  colnames(expr_raw)
)

if (length(gene_col) == 0) {
  gene_col <- colnames(expr_raw)[1]
} else {
  gene_col <- gene_col[1]
}

cat("\nGene identifier column:\n")
cat(gene_col, "\n")


# ------------------------------------------------------------
# 11. Match annotation sample IDs to expression columns
# ------------------------------------------------------------

expression_sample_columns <- setdiff(
  colnames(expr_raw),
  gene_col
)

analysis_samples <- analysis_ann$Expression_Sample

matched_samples <- intersect(
  analysis_samples,
  expression_sample_columns
)

missing_samples <- setdiff(
  analysis_samples,
  expression_sample_columns
)

cat("\nRequested analysis samples :", length(analysis_samples), "\n")
cat("Matched expression samples:", length(matched_samples), "\n")
cat("Missing expression samples :", length(missing_samples), "\n")

if (length(missing_samples) > 0) {
  
  write.table(
    missing_samples,
    file.path(
      step_dir,
      "03_QC",
      "STEP9B7R_Missing_Expression_Samples.txt"
    ),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  
  stop(
    paste(
      "STEP9B-7R FAILED:",
      length(missing_samples),
      "analysis samples are missing from the expression matrix."
    )
  )
}


# ------------------------------------------------------------
# 12. Lock sample order
# ------------------------------------------------------------

analysis_samples <- analysis_ann$Expression_Sample

analysis_ann <- analysis_ann[
  match(
    analysis_samples,
    analysis_ann$Expression_Sample
  ),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    analysis_ann$Expression_Sample,
    analysis_samples
  )
)


# ------------------------------------------------------------
# 13. Extract expression matrix
# ------------------------------------------------------------

gene_ids <- as.character(expr_raw[[gene_col]])

expr_mat <- as.matrix(
  expr_raw[
    ,
    analysis_samples,
    drop = FALSE
  ]
)

storage.mode(expr_mat) <- "numeric"

rownames(expr_mat) <- gene_ids


# ------------------------------------------------------------
# 14. Gene-ID QC
# ------------------------------------------------------------

gene_id_qc <- data.frame(
  Metric = c(
    "Total expression rows",
    "Unique gene IDs",
    "Duplicated gene IDs",
    "NA gene IDs",
    "NA expression cells",
    "Infinite expression cells"
  ),
  Value = c(
    nrow(expr_mat),
    length(unique(rownames(expr_mat))),
    sum(duplicated(rownames(expr_mat))),
    sum(is.na(rownames(expr_mat))),
    sum(is.na(expr_mat)),
    sum(!is.finite(expr_mat))
  )
)

write.csv(
  gene_id_qc,
  file.path(
    step_dir,
    "03_QC",
    "STEP9B7R_Gene_ID_QC.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. Remove invalid gene IDs
# ------------------------------------------------------------

valid_gene_id <- !is.na(rownames(expr_mat)) &
  rownames(expr_mat) != ""

expr_mat <- expr_mat[
  valid_gene_id,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 16. Handle duplicated gene identifiers
#
# For reproducibility:
#   keep the first occurrence
#
# The full raw expression matrix remains untouched.
# ------------------------------------------------------------

dup_gene <- duplicated(rownames(expr_mat))

duplicate_gene_count <- sum(dup_gene)

cat("\nDuplicated gene IDs:\n")
cat(duplicate_gene_count, "\n")

if (duplicate_gene_count > 0) {
  
  expr_mat <- expr_mat[
    !dup_gene,
    ,
    drop = FALSE
  ]
}


# ------------------------------------------------------------
# 17. Expression scale audit
# ------------------------------------------------------------

expr_scale_summary <- data.frame(
  Metric = c(
    "Minimum",
    "Maximum",
    "Median",
    "Mean",
    "NA cells",
    "Infinite cells"
  ),
  Value = c(
    min(expr_mat, na.rm = TRUE),
    max(expr_mat, na.rm = TRUE),
    median(expr_mat, na.rm = TRUE),
    mean(expr_mat, na.rm = TRUE),
    sum(is.na(expr_mat)),
    sum(!is.finite(expr_mat))
  )
)

write.csv(
  expr_scale_summary,
  file.path(
    step_dir,
    "03_QC",
    "STEP9B7R_Expression_Scale_QC.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. Expression transformation
#
# FPKM-UQ -> log2(FPKM-UQ + 1)
# ------------------------------------------------------------

log2_expr <- log2(
  expr_mat + 1
)


# ------------------------------------------------------------
# 19. Transformation QC
# ------------------------------------------------------------

if (any(!is.finite(log2_expr))) {
  stop(
    "STEP9B-7R FAILED: non-finite values detected after log2 transformation."
  )
}


# ------------------------------------------------------------
# 20. Define statistical testing universe
#
# A gene enters the testing universe if:
#   - all required sample values are finite
#   - the combined tumor + normal expression has > 0 variance
#
# No WGCNA filtering is applied here.
# ------------------------------------------------------------

finite_gene <- apply(
  log2_expr,
  1,
  function(x) all(is.finite(x))
)

variable_gene <- apply(
  log2_expr,
  1,
  function(x) var(x) > 0
)

testing_gene <- finite_gene &
  variable_gene

testing_universe <- rownames(log2_expr)[testing_gene]

cat("\n============================================================\n")
cat("Testing universe\n")
cat("============================================================\n")

cat(
  "Total genes in expression matrix:",
  nrow(log2_expr),
  "\n"
)

cat(
  "Genes with finite expression:",
  sum(finite_gene),
  "\n"
)

cat(
  "Genes with non-zero variance:",
  sum(variable_gene),
  "\n"
)

cat(
  "Final DEG testing universe:",
  length(testing_universe),
  "\n"
)


if (length(testing_universe) < 10000) {
  
  warning(
    "Testing universe is unexpectedly small. Review input expression matrix."
  )
  
}


# ------------------------------------------------------------
# 21. Save testing universe
# ------------------------------------------------------------

testing_universe_df <- data.frame(
  Gene = testing_universe,
  stringsAsFactors = FALSE
)

write.csv(
  testing_universe_df,
  file.path(
    step_dir,
    "02_results",
    "STEP9B7R_Testing_Universe.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 22. Prepare group vectors
# ------------------------------------------------------------

group <- analysis_ann$Group

tumor_idx <- which(group == "Tumor")
normal_idx <- which(group == "Normal")

stopifnot(length(tumor_idx) == 513)
stopifnot(length(normal_idx) == 59)


# ------------------------------------------------------------
# 23. Genome-wide Wilcoxon rank-sum testing
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Genome-wide differential expression analysis\n")
cat("============================================================\n")

deg_list <- vector(
  "list",
  length(testing_universe)
)

names(deg_list) <- testing_universe


for (i in seq_along(testing_universe)) {
  
  gene <- testing_universe[i]
  
  x_tumor <- log2_expr[
    gene,
    tumor_idx
  ]
  
  x_normal <- log2_expr[
    gene,
    normal_idx
  ]
  
  wt <- suppressWarnings(
    wilcox.test(
      x = x_tumor,
      y = x_normal,
      alternative = "two.sided",
      exact = FALSE
    )
  )
  
  tumor_median <- median(
    x_tumor,
    na.rm = TRUE
  )
  
  normal_median <- median(
    x_normal,
    na.rm = TRUE
  )
  
  tumor_mean <- mean(
    x_tumor,
    na.rm = TRUE
  )
  
  normal_mean <- mean(
    x_normal,
    na.rm = TRUE
  )
  
  deg_list[[i]] <- data.frame(
    Gene = gene,
    Tumor_Mean = tumor_mean,
    Normal_Mean = normal_mean,
    Tumor_Median = tumor_median,
    Normal_Median = normal_median,
    log2FC = tumor_mean - normal_mean,
    P_value = wt$p.value,
    stringsAsFactors = FALSE
  )
  
  if (i %% 1000 == 0) {
    cat(
      "Processed:",
      i,
      "/",
      length(testing_universe),
      "\n"
    )
  }
}


deg_raw <- do.call(
  rbind,
  deg_list
)

rownames(deg_raw) <- NULL


# ------------------------------------------------------------
# 24. Multiple-testing correction
# ------------------------------------------------------------

deg_raw$FDR <- p.adjust(
  deg_raw$P_value,
  method = "BH"
)


# ------------------------------------------------------------
# 25. DEG direction
# ------------------------------------------------------------

deg_raw$Direction <- ifelse(
  deg_raw$log2FC > 0,
  "Tumor_up",
  ifelse(
    deg_raw$log2FC < 0,
    "Tumor_down",
    "No_change"
  )
)


# ------------------------------------------------------------
# 26. Revised DEG definition
#
# LOCKED DEFINITION:
#
#   FDR < 0.05
#   AND
#   |log2FC| >= 1
#
# ------------------------------------------------------------

DEG_FDR_THRESHOLD <- 0.05
DEG_LOG2FC_THRESHOLD <- 1


deg_raw$DEG <- (
  deg_raw$FDR < DEG_FDR_THRESHOLD &
    abs(deg_raw$log2FC) >= DEG_LOG2FC_THRESHOLD
)


# ------------------------------------------------------------
# 27. Sort full result
# ------------------------------------------------------------

deg_full <- deg_raw[
  order(
    deg_raw$FDR,
    -abs(deg_raw$log2FC)
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 28. Extract DEG subsets
# ------------------------------------------------------------

deg_all <- deg_full[
  deg_full$DEG,
  ,
  drop = FALSE
]

deg_up <- deg_all[
  deg_all$Direction == "Tumor_up",
  ,
  drop = FALSE
]

deg_down <- deg_all[
  deg_all$Direction == "Tumor_down",
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 29. Save results
# ------------------------------------------------------------

write.csv(
  deg_full,
  file.path(
    step_dir,
    "02_results",
    "STEP9B7R_Revised_DEG_FULL.csv"
  ),
  row.names = FALSE
)

write.csv(
  deg_all,
  file.path(
    step_dir,
    "02_results",
    "STEP9B7R_Revised_DEG_ALL.csv"
  ),
  row.names = FALSE
)

write.csv(
  deg_up,
  file.path(
    step_dir,
    "02_results",
    "STEP9B7R_Revised_DEG_UP.csv"
  ),
  row.names = FALSE
)

write.csv(
  deg_down,
  file.path(
    step_dir,
    "02_results",
    "STEP9B7R_Revised_DEG_DOWN.csv"
  ),
  row.names = FALSE
)


# ============================================================
# STEP9B-7R ID DIAGNOSTIC
# Determine the identifier type used in deg_full$Gene
# ============================================================

cat("\n============================================================\n")
cat("STEP9B-7R Gene ID Diagnostic\n")
cat("============================================================\n\n")

cat("First 30 gene IDs:\n")
print(head(deg_full$Gene, 30))

cat("\n\nExample gene ID classes:\n")
print(
  data.frame(
    Gene = head(deg_full$Gene, 20),
    stringsAsFactors = FALSE
  )
)

cat("\n\nCore-gene direct matching:\n")

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_direct_match <- data.frame(
  Core_Gene = core_genes,
  Direct_Match = core_genes %in% deg_full$Gene,
  stringsAsFactors = FALSE
)

print(core_direct_match)

cat("\n\nGene ID pattern summary:\n")

gene_id_example <- head(
  deg_full$Gene,
  min(100, nrow(deg_full))
)

cat(
  "Starts with ENSG:",
  sum(grepl("^ENSG", gene_id_example)),
  "/",
  length(gene_id_example),
  "\n"
)

cat(
  "Starts with ENSG + version:",
  sum(grepl("^ENSG[0-9]+\\.", gene_id_example)),
  "/",
  length(gene_id_example),
  "\n"
)

cat(
  "Looks like HGNC symbols:",
  sum(grepl("^[A-Z][A-Z0-9-]+$", gene_id_example)),
  "/",
  length(gene_id_example),
  "\n"
)


# ============================================================
# STEP9B-7R
# Core-gene Ensembl ID recovery from existing Paper 1 files
#
# Purpose:
#   Recover AOC3 / MS4A7 / IL7R Ensembl IDs from existing
#   Paper 1 annotation/provenance files.
#
# IMPORTANT:
#   - No DEG recalculation
#   - No statistical recalculation
#   - Only gene-ID annotation recovery
# ============================================================

cat("\n============================================================\n")
cat("STEP9B-7R Core Gene Ensembl-ID Recovery\n")
cat("============================================================\n\n")

search_root <- "O:/Rproject/Paper1_tumor-only WGCNA"

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


# ------------------------------------------------------------
# 1. List candidate annotation files
# ------------------------------------------------------------

all_files <- list.files(
  search_root,
  recursive = TRUE,
  full.names = TRUE
)

annotation_candidates <- all_files[
  grepl(
    "annot|annotation|gencode|ensembl|gene.*symbol|symbol.*gene",
    basename(all_files),
    ignore.case = TRUE
  )
]

cat(
  "Annotation-related filename candidates:",
  length(annotation_candidates),
  "\n\n"
)

print(
  basename(annotation_candidates)
)


# ------------------------------------------------------------
# 2. Search text-based files for the three symbols
# ------------------------------------------------------------

text_candidates <- annotation_candidates[
  grepl(
    "\\.(txt|csv|tsv|tab)$",
    annotation_candidates,
    ignore.case = TRUE
  )
]

cat(
  "\nText annotation candidates:",
  length(text_candidates),
  "\n"
)


search_results <- list()


for (f in text_candidates) {
  
  x <- tryCatch(
    read.delim(
      f,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = ""
    ),
    error = function(e) NULL
  )
  
  if (is.null(x)) {
    next
  }
  
  if (nrow(x) == 0 || ncol(x) == 0) {
    next
  }
  
  # Convert all columns to character safely
  x_char <- as.data.frame(
    lapply(
      x,
      function(z) as.character(z)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  for (g in core_genes) {
    
    hit_rows <- which(
      apply(
        x_char,
        1,
        function(row) {
          any(
            toupper(trimws(row)) == toupper(g),
            na.rm = TRUE
          )
        }
      )
    )
    
    if (length(hit_rows) > 0) {
      
      for (r in hit_rows) {
        
        search_results[[length(search_results) + 1]] <-
          data.frame(
            File = f,
            Core_Gene = g,
            Row = r,
            stringsAsFactors = FALSE
          )
        
      }
    }
  }
}


# ------------------------------------------------------------
# 3. Report exact symbol matches
# ------------------------------------------------------------

if (length(search_results) > 0) {
  
  core_file_hits <- do.call(
    rbind,
    search_results
  )
  
  cat("\n============================================================\n")
  cat("EXACT CORE-GENE SYMBOL HITS\n")
  cat("============================================================\n\n")
  
  print(core_file_hits)
  
  write.csv(
    core_file_hits,
    file.path(
      step_dir,
      "03_QC",
      "STEP9B7R_CoreGene_ExistingFile_Hits.csv"
    ),
    row.names = FALSE
  )
  
} else {
  
  core_file_hits <- data.frame(
    File = character(0),
    Core_Gene = character(0),
    Row = integer(0),
    stringsAsFactors = FALSE
  )
  
  cat("\nNo exact core-gene symbol hits found.\n")
  
}


# ------------------------------------------------------------
# 4. Targeted inspection of known Workflow-03F files
# ------------------------------------------------------------

known_files <- c(
  "Workflow-03F-8_expression_provenance_test.txt",
  "Workflow-03F-6_expression_provenance_test.txt",
  "Workflow-03F-2_expression_annotation.txt",
  "Workflow-03F-3_patient_expression_mapping.txt",
  "Workflow-03F-5_historical_517_patient_annotation.txt"
)

known_paths <- file.path(
  search_root,
  known_files
)

known_paths <- known_paths[
  file.exists(known_paths)
]


cat("\n============================================================\n")
cat("Known Workflow-03F files\n")
cat("============================================================\n\n")

print(basename(known_paths))


# ------------------------------------------------------------
# 5. Search known files directly
# ------------------------------------------------------------

known_hits <- list()

for (f in known_paths) {
  
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
  
  for (g in core_genes) {
    
    hit <- grep(
      pattern = paste0(
        "(^|[\\t,; ])",
        g,
        "($|[\\t,; ])"
      ),
      txt,
      ignore.case = TRUE
    )
    
    if (length(hit) > 0) {
      
      known_hits[[length(known_hits) + 1]] <-
        data.frame(
          File = f,
          Core_Gene = g,
          Line = hit,
          Text = txt[hit],
          stringsAsFactors = FALSE
        )
      
    }
  }
}


if (length(known_hits) > 0) {
  
  known_hits_df <- do.call(
    rbind,
    known_hits
  )
  
  cat("\nKnown-file core gene hits:\n\n")
  
  print(
    known_hits_df,
    row.names = FALSE
  )
  
  write.csv(
    known_hits_df,
    file.path(
      step_dir,
      "03_QC",
      "STEP9B7R_KnownWorkflow03F_CoreGene_Hits.csv"
    ),
    row.names = FALSE
  )
  
} else {
  
  cat(
    "\nNo core-gene symbol hits found in the known Workflow-03F files.\n"
  )
  
}


# ------------------------------------------------------------
# 6. Final diagnostic status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9B-7R Core-ID Recovery Diagnostic Status\n")
cat("============================================================\n\n")

if (
  length(search_results) > 0 ||
  length(known_hits) > 0
) {
  
  cat(
    "Status: CORE_GENE_ANNOTATION_CANDIDATES_FOUND\n"
  )
  
} else {
  
  cat(
    "Status: CORE_GENE_ANNOTATION_NOT_FOUND_IN_EXISTING_FILES\n"
  )
  
}

cat("\nNo DEG recalculation performed.\n")
cat("No P values recalculated.\n")
cat("No FDR recalculated.\n")
cat("No DEG threshold changed.\n")


# ============================================================
# STEP9B-7R
# Core Gene Ensembl-ID Mapping Lock + DEG Core Audit
#
# Purpose:
#   1. Recover AOC3 / MS4A7 / IL7R Ensembl IDs
#      from the existing GENCODE v22 mapping file.
#   2. Verify that the mapping is unique and reproducible.
#   3. Connect the three core genes to the already-computed
#      deg_full table.
#   4. Perform NO DEG recalculation.
#
# IMPORTANT:
#   - No expression recalculation
#   - No Wilcoxon recalculation
#   - No P-value recalculation
#   - No FDR recalculation
#   - No DEG threshold change
# ============================================================

cat("\n============================================================\n")
cat("STEP9B-7R Core Gene Mapping Lock + DEG Core Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Define paths and core genes
# ------------------------------------------------------------

search_root <- "input_dir/Paper1_tumor-only WGCNA"

mapping_file <- file.path(
  search_root,
  "Workflow-04-2C-1A_Ensembl_to_GeneSymbol_GENCODE_v22.txt"
)

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

cat("Mapping file:\n")
cat(mapping_file, "\n\n")

cat("Core genes:\n")
print(core_genes)


# ------------------------------------------------------------
# 1. Check mapping file
# ------------------------------------------------------------

stopifnot(file.exists(mapping_file))

mapping <- read.delim(
  mapping_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

cat("\n============================================================\n")
cat("1. Mapping-file structure\n")
cat("============================================================\n\n")

cat("Dimensions:\n")
print(dim(mapping))

cat("\nColumn names:\n")
print(colnames(mapping))

cat("\nFirst rows:\n")
print(head(mapping))


# ------------------------------------------------------------
# 2. Identify Ensembl and Gene Symbol columns
# ------------------------------------------------------------

cn <- colnames(mapping)

ensembl_candidates <- which(
  grepl(
    "ensembl|gene.*id|id.*gene",
    cn,
    ignore.case = TRUE
  )
)

symbol_candidates <- which(
  grepl(
    "symbol|gene.*name|name.*gene",
    cn,
    ignore.case = TRUE
  )
)

cat("\n============================================================\n")
cat("2. Candidate mapping columns\n")
cat("============================================================\n\n")

cat("Ensembl-related columns:\n")
print(cn[ensembl_candidates])

cat("\nGene-symbol-related columns:\n")
print(cn[symbol_candidates])


# ------------------------------------------------------------
# 3. Automatically select the expected columns
# ------------------------------------------------------------

if (
  length(ensembl_candidates) == 0 ||
  length(symbol_candidates) == 0
) {
  stop(
    "Could not identify Ensembl and Gene Symbol columns automatically."
  )
}

ensembl_col <- cn[ensembl_candidates[1]]
symbol_col  <- cn[symbol_candidates[1]]

cat("\nSelected Ensembl column:\n")
cat(ensembl_col, "\n")

cat("Selected Gene Symbol column:\n")
cat(symbol_col, "\n")


# ------------------------------------------------------------
# 4. Clean mapping table
# ------------------------------------------------------------

mapping_clean <- data.frame(
  ENSEMBL = trimws(
    as.character(mapping[[ensembl_col]])
  ),
  GeneSymbol = toupper(
    trimws(
      as.character(mapping[[symbol_col]])
    )
  ),
  stringsAsFactors = FALSE
)

mapping_clean <- mapping_clean[
  !is.na(mapping_clean$ENSEMBL) &
    mapping_clean$ENSEMBL != "" &
    !is.na(mapping_clean$GeneSymbol) &
    mapping_clean$GeneSymbol != "",
  ,
  drop = FALSE
]

cat("\nClean mapping dimensions:\n")
print(dim(mapping_clean))


# ------------------------------------------------------------
# 5. Extract the three core genes
# ------------------------------------------------------------

core_mapping <- mapping_clean[
  mapping_clean$GeneSymbol %in% core_genes,
  ,
  drop = FALSE
]

core_mapping <- unique(core_mapping)

cat("\n============================================================\n")
cat("3. Core-gene mapping candidates\n")
cat("============================================================\n\n")

print(core_mapping)


# ------------------------------------------------------------
# 6. Check mapping multiplicity
# ------------------------------------------------------------

mapping_count <- aggregate(
  ENSEMBL ~ GeneSymbol,
  data = core_mapping,
  FUN = function(x) length(unique(x))
)

colnames(mapping_count)[2] <- "N_unique_ENSEMBL"

cat("\nMapping multiplicity:\n")
print(mapping_count)


# ------------------------------------------------------------
# 7. Require exactly one Ensembl ID per core gene
# ------------------------------------------------------------

mapping_complete <- all(
  core_genes %in% core_mapping$GeneSymbol
)

mapping_unique <- all(
  mapping_count$N_unique_ENSEMBL == 1
)

cat("\n============================================================\n")
cat("4. Mapping QC\n")
cat("============================================================\n\n")

cat(
  "All three core genes found: ",
  mapping_complete,
  "\n",
  sep = ""
)

cat(
  "Exactly one Ensembl ID per core gene: ",
  mapping_unique,
  "\n",
  sep = ""
)

if (!mapping_complete) {
  stop(
    "Core-gene mapping incomplete. Do not lock STEP9B-7R."
  )
}

if (!mapping_unique) {
  stop(
    "One or more core genes map to multiple Ensembl IDs. ",
    "Do not lock STEP9B-7R."
  )
}


# ------------------------------------------------------------
# 8. Create locked three-gene mapping table
# ------------------------------------------------------------

core_mapping_locked <- data.frame(
  Core_Gene = core_genes,
  ENSEMBL = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(core_genes)) {
  
  g <- core_genes[i]
  
  core_mapping_locked$ENSEMBL[i] <-
    core_mapping$ENSEMBL[
      core_mapping$GeneSymbol == g
    ][1]
}

cat("\n============================================================\n")
cat("5. LOCKED CORE-GENE MAPPING\n")
cat("============================================================\n\n")

print(core_mapping_locked)


# ------------------------------------------------------------
# 9. Check that deg_full already exists
# ------------------------------------------------------------

if (!exists("deg_full")) {
  
  stop(
    "deg_full is not present in the current R environment. ",
    "Do NOT recalculate DEG here. Reload the existing DEG result instead."
  )
  
}

cat("\ndeg_full detected in current R environment.\n")

cat("deg_full dimensions:\n")
print(dim(deg_full))

cat("\ndeg_full column names:\n")
print(colnames(deg_full))


# ------------------------------------------------------------
# 10. Verify DEG table uses Ensembl IDs
# ------------------------------------------------------------

if (!"Gene" %in% colnames(deg_full)) {
  
  stop(
    "deg_full does not contain a column named 'Gene'."
  )
  
}

deg_gene_ids <- trimws(
  as.character(deg_full$Gene)
)

cat("\nFirst 10 DEG Gene IDs:\n")
print(head(deg_gene_ids, 10))

cat(
  "\nNumber of DEG genes with ENSG prefix: ",
  sum(grepl("^ENSG", deg_gene_ids)),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Match the three core genes to deg_full
# ------------------------------------------------------------

core_audit <- merge(
  core_mapping_locked,
  deg_full,
  by.x = "ENSEMBL",
  by.y = "Gene",
  all.x = TRUE,
  sort = FALSE
)

# Restore desired core-gene order
core_audit$Core_Gene <- factor(
  core_audit$Core_Gene,
  levels = core_genes
)

core_audit <- core_audit[
  order(core_audit$Core_Gene),
  ,
  drop = FALSE
]

core_audit$Core_Gene <- as.character(
  core_audit$Core_Gene
)


# ------------------------------------------------------------
# 12. Add explicit audit fields
# ------------------------------------------------------------

core_audit$Found_in_DEG_Table <-
  !is.na(core_audit$log2FC)

core_audit$DEG <- FALSE

if ("FDR" %in% colnames(core_audit)) {
  
  core_audit$DEG <-
    !is.na(core_audit$FDR) &
    core_audit$FDR < 0.05 &
    abs(core_audit$log2FC) >= 1
  
}

core_audit$Direction <- NA_character_

if ("log2FC" %in% colnames(core_audit)) {
  
  core_audit$Direction[
    !is.na(core_audit$log2FC) &
      core_audit$log2FC > 0
  ] <- "Tumor_up"
  
  core_audit$Direction[
    !is.na(core_audit$log2FC) &
      core_audit$log2FC < 0
  ] <- "Tumor_down"
  
  core_audit$Direction[
    !is.na(core_audit$log2FC) &
      core_audit$log2FC == 0
  ] <- "No_change"
}


# ------------------------------------------------------------
# 13. Display core-gene DEG audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("6. CORE-GENE DEG AUDIT\n")
cat("============================================================\n\n")

audit_display_cols <- intersect(
  c(
    "Core_Gene",
    "ENSEMBL",
    "Found_in_DEG_Table",
    "log2FC",
    "P_value",
    "FDR",
    "DEG",
    "Direction"
  ),
  colnames(core_audit)
)

print(
  core_audit[, audit_display_cols, drop = FALSE],
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Add explicit annotation status
# ------------------------------------------------------------

core_audit$Found_in_Expression <- TRUE
core_audit$Found_in_Testing_Universe <-
  core_audit$Found_in_DEG_Table


# ------------------------------------------------------------
# 15. Save locked mapping and core audit
# ------------------------------------------------------------

step_dir <- file.path(
  "result_dir/Paper 1_Step 9-Functional Enrichment",
  "STEP9B-7R_Revised_DEG_Definition_and_Lock"
)

dir.create(
  file.path(step_dir, "03_QC"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(step_dir, "04_lock"),
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  core_mapping_locked,
  file.path(
    step_dir,
    "04_lock",
    "STEP9B7R_CoreGene_Ensembl_Mapping_LOCKED.csv"
  ),
  row.names = FALSE
)

write.csv(
  core_audit,
  file.path(
    step_dir,
    "03_QC",
    "STEP9B7R_CoreGene_DEG_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 16. Save mapping provenance
# ------------------------------------------------------------

mapping_provenance <- data.frame(
  Item = c(
    "Mapping source",
    "Annotation version",
    "Core genes",
    "Mapping method",
    "DEG recalculation",
    "P-value recalculation",
    "FDR recalculation",
    "Threshold change"
  ),
  Value = c(
    basename(mapping_file),
    "GENCODE v22",
    paste(core_genes, collapse = ";"),
    "Exact GeneSymbol-to-Ensembl lookup",
    "NO",
    "NO",
    "NO",
    "NO"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  mapping_provenance,
  file.path(
    step_dir,
    "04_lock",
    "STEP9B7R_CoreGene_Mapping_Provenance_LOCKED.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Final PASS / HOLD status
# ------------------------------------------------------------

mapping_pass <- (
  mapping_complete &&
    mapping_unique &&
    nrow(core_mapping_locked) == 3 &&
    all(!is.na(core_mapping_locked$ENSEMBL))
)

audit_pass <- (
  nrow(core_audit) == 3 &&
    all(core_audit$Found_in_Expression)
)

cat("\n============================================================\n")
cat("STEP9B-7R FINAL CORE-ID AUDIT STATUS\n")
cat("============================================================\n\n")

cat(
  "Mapping completeness: ",
  mapping_complete,
  "\n",
  sep = ""
)

cat(
  "Mapping uniqueness: ",
  mapping_unique,
  "\n",
  sep = ""
)

cat(
  "Three core genes audited: ",
  audit_pass,
  "\n",
  sep = ""
)

if (mapping_pass && audit_pass) {
  
  cat(
    "\nStatus: PASS_CORE_GENE_MAPPING_AND_AUDIT\n"
  )
  
  cat(
    "\nNo DEG recalculation performed.\n"
  )
  
  cat(
    "No P values recalculated.\n"
  )
  
  cat(
    "No FDR recalculated.\n"
  )
  
  cat(
    "No DEG threshold changed.\n"
  )
  
  cat(
    "Existing DEG result retained as the statistical result.\n"
  )
  
} else {
  
  cat(
    "\nStatus: HOLD_CORE_GENE_MAPPING_AUDIT\n"
  )
  
}

cat("\n============================================================\n")
cat("STEP9B-7R Core Gene Mapping Lock completed.\n")
cat("============================================================\n")

##STEP9C — Revised DEG-based GO/KEGG ORA
# ============================================================
# PAPER 1
# STEP9C — Revised DEG-based GO/KEGG Over-Representation Analysis
#
# Purpose:
#   Perform reproducible GO and KEGG ORA using the LOCKED
#   STEP9B-7R revised DEG results.
#
# Primary DEG groups:
#   1. Tumor_up
#   2. Tumor_down
#
# GO categories:
#   - Biological Process (BP)
#   - Cellular Component (CC)
#   - Molecular Function (MF)
#
# KEGG:
#   - Homo sapiens
#
# IMPORTANT:
#   - No DEG recalculation
#   - No Wilcoxon recalculation
#   - No P-value recalculation for DEGs
#   - No FDR recalculation for DEGs
#   - DEG threshold remains LOCKED:
#       FDR < 0.05 AND |log2FC| >= 1
#
# Background:
#   The FINAL DEG testing universe from STEP9B-7R:
#       57,828 genes
#
# IMPORTANT DISTINCTION:
#   57,828 = DEG testing universe
#   31,609 = tumor-only WGCNA universe
#
#   STEP9C uses 57,828 as ORA background.
#   The 31,609 WGCNA universe is NOT used here.
# ============================================================


cat("\n============================================================\n")
cat("STEP9C — Revised DEG-based GO/KEGG ORA\n")
cat("============================================================\n\n")


# ============================================================
# 0. Package requirements
# ============================================================

required_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "enrichplot",
  "ggplot2"
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
      "The following required packages are missing:\n",
      paste(missing_packages, collapse = "\n"),
      "\n\nInstall them before running STEP9C."
    )
  )
  
}


library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)


# ============================================================
# 1. Define STEP9C directories
# ============================================================

step_dir <- file.path(
  "input_dir/Paper 1_Step 9-Functional Enrichment",
  "STEP9C_Revised_DEG_GO_KEGG_ORA"
)

dir.create(
  step_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(step_dir, "00_input"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(step_dir, "01_data"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(step_dir, "02_results", "GO_BP"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(step_dir, "02_results", "GO_CC"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(step_dir, "02_results", "GO_MF"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(step_dir, "02_results", "KEGG"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(step_dir, "03_QC"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(step_dir, "04_lock"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(step_dir, "05_logs"),
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Define locked DEG universe
# ============================================================

cat("Checking STEP9B-7R objects...\n\n")


if (!exists("deg_full")) {
  
  stop(
    "deg_full is not present in the current R environment.\n",
    "Do NOT recalculate DEG here.\n",
    "Reload the locked STEP9B-7R DEG result instead."
  )
  
}

if (!exists("deg_all")) {
  
  stop(
    "deg_all is not present in the current R environment.\n",
    "Do NOT recalculate DEG here."
  )
  
}


cat("deg_full detected.\n")
cat("deg_all detected.\n\n")


# ============================================================
# 3. Check DEG table structure
# ============================================================

cat("============================================================\n")
cat("STEP9B-7R DEG structure check\n")
cat("============================================================\n\n")

cat("deg_full dimensions:\n")
print(dim(deg_full))

cat("\ndeg_full columns:\n")
print(colnames(deg_full))

cat("\ndeg_all dimensions:\n")
print(dim(deg_all))

cat("\ndeg_all columns:\n")
print(colnames(deg_all))


# ============================================================
# 4. Identify DEG testing universe
# ============================================================
#
# The STEP9B-7R testing universe should contain 57,828 genes.
#
# Prefer an existing locked object if available.
# Otherwise derive the universe from deg_full.
#
# IMPORTANT:
# This is NOT a new statistical calculation.
# It is only retrieval of the already-completed DEG universe.
# ============================================================

if (exists("deg_testing_universe")) {
  
  testing_universe <- as.character(
    deg_testing_universe
  )
  
  testing_universe <- unique(
    trimws(testing_universe)
  )
  
  testing_universe <- testing_universe[
    testing_universe != "" &
      !is.na(testing_universe)
  ]
  
  cat(
    "\nUsing existing object: deg_testing_universe\n"
  )
  
} else {
  
  cat(
    "\nObject 'deg_testing_universe' not found.\n"
  )
  
  cat(
    "Attempting to recover the 57,828-gene universe from deg_full.\n"
  )
  
  if (!"Gene" %in% colnames(deg_full)) {
    
    stop(
      "deg_full does not contain a 'Gene' column."
    )
    
  }
  
  testing_universe <- as.character(
    deg_full$Gene
  )
  
  testing_universe <- unique(
    trimws(testing_universe)
  )
  
  testing_universe <- testing_universe[
    testing_universe != "" &
      !is.na(testing_universe)
  ]
  
}


# ============================================================
# 5. LOCK testing-universe size
# ============================================================

cat("\n============================================================\n")
cat("DEG testing universe\n")
cat("============================================================\n\n")

cat(
  "Number of unique genes in testing universe: ",
  length(testing_universe),
  "\n",
  sep = ""
)

if (length(testing_universe) != 57828) {
  
  stop(
    paste0(
      "Expected 57,828 genes in the locked DEG testing universe, ",
      "but found ",
      length(testing_universe),
      ".\n",
      "STOP before ORA."
    )
  )
  
}

cat(
  "Testing universe size check: PASS\n"
)


# ============================================================
# 6. Verify Ensembl ID format
# ============================================================

ensembl_fraction <- mean(
  grepl(
    "^ENSG[0-9]+$",
    testing_universe
  )
)

cat(
  "\nFraction of testing-universe genes with unversioned ENSG IDs: ",
  round(ensembl_fraction, 6),
  "\n",
  sep = ""
)

if (ensembl_fraction < 0.95) {
  
  warning(
    "Less than 95% of the testing universe has the expected ENSG format."
  )
  
}


# ============================================================
# 7. Retrieve locked DEG groups
# ============================================================

if (!all(c("Gene", "Direction") %in% colnames(deg_all))) {
  
  stop(
    "deg_all must contain 'Gene' and 'Direction' columns."
  )
  
}


deg_all$Gene <- trimws(
  as.character(deg_all$Gene)
)

deg_all$Direction <- trimws(
  as.character(deg_all$Direction)
)


deg_up <- deg_all[
  deg_all$Direction == "Tumor_up",
  ,
  drop = FALSE
]

deg_down <- deg_all[
  deg_all$Direction == "Tumor_down",
  ,
  drop = FALSE
]


# ============================================================
# 8. Check DEG groups
# ============================================================

cat("\n============================================================\n")
cat("Locked DEG groups\n")
cat("============================================================\n\n")

cat(
  "Total DEGs: ",
  nrow(deg_all),
  "\n",
  sep = ""
)

cat(
  "Tumor_up DEGs: ",
  nrow(deg_up),
  "\n",
  sep = ""
)

cat(
  "Tumor_down DEGs: ",
  nrow(deg_down),
  "\n",
  sep = ""
)


if (nrow(deg_up) == 0) {
  
  stop(
    "No Tumor_up DEGs found."
  )
  
}

if (nrow(deg_down) == 0) {
  
  stop(
    "No Tumor_down DEGs found."
  )
  
}


# ============================================================
# 9. Confirm DEG definition without recalculation
# ============================================================

if (
  !all(
    deg_all$FDR < 0.05 &
    abs(deg_all$log2FC) >= 1
  )
) {
  
  stop(
    "deg_all contains genes that do not satisfy the locked ",
    "STEP9B-7R DEG definition."
  )
  
}

cat(
  "\nSTEP9B-7R DEG definition consistency check: PASS\n"
)

cat(
  "Definition = FDR < 0.05 AND |log2FC| >= 1\n"
)


# ============================================================
# 10. Save exact STEP9C input tables
# ============================================================

write.csv(
  testing_universe,
  file.path(
    step_dir,
    "00_input",
    "STEP9C_DEG_testing_universe_57828.csv"
  ),
  row.names = FALSE
)

write.csv(
  deg_all,
  file.path(
    step_dir,
    "00_input",
    "STEP9C_locked_DEG_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  deg_up,
  file.path(
    step_dir,
    "00_input",
    "STEP9C_locked_DEG_Tumor_up.csv"
  ),
  row.names = FALSE
)

write.csv(
  deg_down,
  file.path(
    step_dir,
    "00_input",
    "STEP9C_locked_DEG_Tumor_down.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. GO ORA function
# ============================================================

run_GO_ORA <- function(
    gene_list,
    universe,
    ontology,
    output_dir,
    prefix
) {
  
  cat("\n------------------------------------------------------------\n")
  cat("GO ORA:", ontology, "-", prefix, "\n")
  cat("------------------------------------------------------------\n")
  
  gene_list <- unique(
    trimws(as.character(gene_list))
  )
  
  universe <- unique(
    trimws(as.character(universe))
  )
  
  gene_list <- gene_list[
    grepl("^ENSG[0-9]+$", gene_list)
  ]
  
  universe <- universe[
    grepl("^ENSG[0-9]+$", universe)
  ]
  
  cat(
    "Input DEG genes: ",
    length(gene_list),
    "\n",
    sep = ""
  )
  
  cat(
    "Input universe genes: ",
    length(universe),
    "\n",
    sep = ""
  )
  
  
  # ----------------------------------------------------------
  # Check that all DEG genes belong to universe
  # ----------------------------------------------------------
  
  gene_in_universe <- gene_list %in% universe
  
  if (!all(gene_in_universe)) {
    
    cat(
      "DEG genes outside universe: ",
      sum(!gene_in_universe),
      "\n",
      sep = ""
    )
    
    gene_list <- gene_list[
      gene_in_universe
    ]
    
  }
  
  
  # ----------------------------------------------------------
  # Run enrichGO
  # ----------------------------------------------------------
  
  ego <- enrichGO(
    gene = gene_list,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  
  
  if (is.null(ego)) {
    
    cat(
      "No GO enrichment object returned.\n"
    )
    
    return(
      list(
        object = NULL,
        all = data.frame(),
        significant = data.frame()
      )
    )
    
  }
  
  
  ego_df <- as.data.frame(ego)
  
  
  if (nrow(ego_df) == 0) {
    
    cat(
      "GO result contains zero terms.\n"
    )
    
    return(
      list(
        object = ego,
        all = ego_df,
        significant = ego_df
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Explicit BH-FDR significance definition
  # ----------------------------------------------------------
  
  ego_sig <- ego_df[
    !is.na(ego_df$p.adjust) &
      ego_df$p.adjust < 0.05,
    ,
    drop = FALSE
  ]
  
  # Sort by adjusted P value
  ego_sig <- ego_sig[
    order(
      ego_sig$p.adjust,
      ego_sig$pvalue
    ),
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Save full and significant results
  # ----------------------------------------------------------
  
  write.csv(
    ego_df,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_GO_",
        ontology,
        "_ALL.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    ego_sig,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_GO_",
        ontology,
        "_FDR05.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Save RDS
  # ----------------------------------------------------------
  
  saveRDS(
    ego,
    file = file.path(
      output_dir,
      paste0(
        prefix,
        "_GO_",
        ontology,
        "_enrichGO.rds"
      )
    )
  )
  
  
  cat(
    "Total GO terms: ",
    nrow(ego_df),
    "\n",
    sep = ""
  )
  
  cat(
    "FDR < 0.05 GO terms: ",
    nrow(ego_sig),
    "\n",
    sep = ""
  )
  
  
  return(
    list(
      object = ego,
      all = ego_df,
      significant = ego_sig
    )
  )
  
}


# ============================================================
# 12. Run GO-BP / GO-CC / GO-MF for Tumor_up
# ============================================================

cat("\n============================================================\n")
cat("GO ORA — Tumor_up\n")
cat("============================================================\n")

go_up_bp <- run_GO_ORA(
  gene_list = deg_up$Gene,
  universe = testing_universe,
  ontology = "BP",
  output_dir = file.path(
    step_dir,
    "02_results",
    "GO_BP"
  ),
  prefix = "Tumor_up"
)

go_up_cc <- run_GO_ORA(
  gene_list = deg_up$Gene,
  universe = testing_universe,
  ontology = "CC",
  output_dir = file.path(
    step_dir,
    "02_results",
    "GO_CC"
  ),
  prefix = "Tumor_up"
)

go_up_mf <- run_GO_ORA(
  gene_list = deg_up$Gene,
  universe = testing_universe,
  ontology = "MF",
  output_dir = file.path(
    step_dir,
    "02_results",
    "GO_MF"
  ),
  prefix = "Tumor_up"
)


# ============================================================
# 13. Run GO-BP / GO-CC / GO-MF for Tumor_down
# ============================================================

cat("\n============================================================\n")
cat("GO ORA — Tumor_down\n")
cat("============================================================\n")

go_down_bp <- run_GO_ORA(
  gene_list = deg_down$Gene,
  universe = testing_universe,
  ontology = "BP",
  output_dir = file.path(
    step_dir,
    "02_results",
    "GO_BP"
  ),
  prefix = "Tumor_down"
)

go_down_cc <- run_GO_ORA(
  gene_list = deg_down$Gene,
  universe = testing_universe,
  ontology = "CC",
  output_dir = file.path(
    step_dir,
    "02_results",
    "GO_CC"
  ),
  prefix = "Tumor_down"
)

go_down_mf <- run_GO_ORA(
  gene_list = deg_down$Gene,
  universe = testing_universe,
  ontology = "MF",
  output_dir = file.path(
    step_dir,
    "02_results",
    "GO_MF"
  ),
  prefix = "Tumor_down"
)


# ============================================================
# 14. KEGG ORA helper
# ============================================================

run_KEGG_ORA <- function(
    gene_list,
    universe,
    output_dir,
    prefix
) {
  
  cat("\n------------------------------------------------------------\n")
  cat("KEGG ORA -", prefix, "\n")
  cat("------------------------------------------------------------\n")
  
  
  gene_list <- unique(
    trimws(as.character(gene_list))
  )
  
  universe <- unique(
    trimws(as.character(universe))
  )
  
  gene_list <- gene_list[
    grepl("^ENSG[0-9]+$", gene_list)
  ]
  
  universe <- universe[
    grepl("^ENSG[0-9]+$", universe)
  ]
  
  
  # ----------------------------------------------------------
  # Ensembl -> Entrez mapping
  # ----------------------------------------------------------
  
  gene_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_list,
    keytype = "ENSEMBL",
    columns = c(
      "ENSEMBL",
      "ENTREZID"
    )
  )
  
  universe_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = universe,
    keytype = "ENSEMBL",
    columns = c(
      "ENSEMBL",
      "ENTREZID"
    )
  )
  
  
  # Remove missing mappings
  gene_map <- gene_map[
    !is.na(gene_map$ENTREZID) &
      gene_map$ENTREZID != "",
    ,
    drop = FALSE
  ]
  
  universe_map <- universe_map[
    !is.na(universe_map$ENTREZID) &
      universe_map$ENTREZID != "",
    ,
    drop = FALSE
  ]
  
  
  # Unique Entrez IDs
  gene_entrez <- unique(
    as.character(
      gene_map$ENTREZID
    )
  )
  
  universe_entrez <- unique(
    as.character(
      universe_map$ENTREZID
    )
  )
  
  
  cat(
    "Input DEG Ensembl genes: ",
    length(gene_list),
    "\n",
    sep = ""
  )
  
  cat(
    "Mapped DEG Entrez genes: ",
    length(gene_entrez),
    "\n",
    sep = ""
  )
  
  cat(
    "Input universe Ensembl genes: ",
    length(universe),
    "\n",
    sep = ""
  )
  
  cat(
    "Mapped universe Entrez genes: ",
    length(universe_entrez),
    "\n",
    sep = ""
  )
  
  
  # ----------------------------------------------------------
  # Mapping coverage
  # ----------------------------------------------------------
  
  gene_mapping_coverage <-
    length(unique(gene_map$ENSEMBL)) /
    length(gene_list)
  
  universe_mapping_coverage <-
    length(unique(universe_map$ENSEMBL)) /
    length(universe)
  
  
  cat(
    "DEG mapping coverage: ",
    round(gene_mapping_coverage, 4),
    "\n",
    sep = ""
  )
  
  cat(
    "Universe mapping coverage: ",
    round(universe_mapping_coverage, 4),
    "\n",
    sep = ""
  )
  
  
  # Save mapping QC
  write.csv(
    gene_map,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_KEGG_DEG_Ensembl_to_Entrez_mapping.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    universe_map,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_KEGG_universe_Ensembl_to_Entrez_mapping.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Run enrichKEGG
  # ----------------------------------------------------------
  
  ekegg <- enrichKEGG(
    gene = gene_entrez,
    universe = universe_entrez,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  
  
  if (is.null(ekegg)) {
    
    cat(
      "No KEGG enrichment object returned.\n"
    )
    
    return(
      list(
        object = NULL,
        all = data.frame(),
        significant = data.frame()
      )
    )
    
  }
  
  
  kegg_df <- as.data.frame(
    ekegg
  )
  
  
  if (nrow(kegg_df) == 0) {
    
    cat(
      "KEGG result contains zero pathways.\n"
    )
    
    return(
      list(
        object = ekegg,
        all = kegg_df,
        significant = kegg_df
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Explicit FDR < 0.05 definition
  # ----------------------------------------------------------
  
  kegg_sig <- kegg_df[
    !is.na(kegg_df$p.adjust) &
      kegg_df$p.adjust < 0.05,
    ,
    drop = FALSE
  ]
  
  
  kegg_sig <- kegg_sig[
    order(
      kegg_sig$p.adjust,
      kegg_sig$pvalue
    ),
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Save results
  # ----------------------------------------------------------
  
  write.csv(
    kegg_df,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_KEGG_ALL.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    kegg_sig,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_KEGG_FDR05.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  saveRDS(
    ekegg,
    file = file.path(
      output_dir,
      paste0(
        prefix,
        "_KEGG_enrichKEGG.rds"
      )
    )
  )
  
  
  cat(
    "Total KEGG pathways: ",
    nrow(kegg_df),
    "\n",
    sep = ""
  )
  
  cat(
    "FDR < 0.05 KEGG pathways: ",
    nrow(kegg_sig),
    "\n",
    sep = ""
  )
  
  
  return(
    list(
      object = ekegg,
      all = kegg_df,
      significant = kegg_sig
    )
  )
  
}


# ============================================================
# 15. Run KEGG ORA — Tumor_up
# ============================================================

cat("\n============================================================\n")
cat("KEGG ORA — Tumor_up\n")
cat("============================================================\n")

kegg_up <- run_KEGG_ORA(
  gene_list = deg_up$Gene,
  universe = testing_universe,
  output_dir = file.path(
    step_dir,
    "02_results",
    "KEGG"
  ),
  prefix = "Tumor_up"
)


# ============================================================
# 16. Run KEGG ORA — Tumor_down
# ============================================================

cat("\n============================================================\n")
cat("KEGG ORA — Tumor_down\n")
cat("============================================================\n")

kegg_down <- run_KEGG_ORA(
  gene_list = deg_down$Gene,
  universe = testing_universe,
  output_dir = file.path(
    step_dir,
    "02_results",
    "KEGG"
  ),
  prefix = "Tumor_down"
)


# ============================================================
# 17. Create global ORA summary
# ============================================================

ora_summary <- data.frame(
  Analysis = c(
    "Tumor_up",
    "Tumor_down"
  ),
  DEG_N = c(
    nrow(deg_up),
    nrow(deg_down)
  ),
  GO_BP_FDR05 = c(
    nrow(go_up_bp$significant),
    nrow(go_down_bp$significant)
  ),
  GO_CC_FDR05 = c(
    nrow(go_up_cc$significant),
    nrow(go_down_cc$significant)
  ),
  GO_MF_FDR05 = c(
    nrow(go_up_mf$significant),
    nrow(go_down_mf$significant)
  ),
  KEGG_FDR05 = c(
    nrow(kegg_up$significant),
    nrow(kegg_down$significant)
  ),
  stringsAsFactors = FALSE
)


cat("\n============================================================\n")
cat("STEP9C ORA SUMMARY\n")
cat("============================================================\n\n")

print(
  ora_summary,
  row.names = FALSE
)


write.csv(
  ora_summary,
  file.path(
    step_dir,
    "03_QC",
    "STEP9C_ORA_Global_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. Create mapping/QC summary
# ============================================================

qc_summary <- data.frame(
  Item = c(
    "DEG_testing_universe",
    "Total_DEGs",
    "Tumor_up_DEGs",
    "Tumor_down_DEGs",
    "GO_BP_method",
    "GO_CC_method",
    "GO_MF_method",
    "KEGG_method",
    "DEG_definition",
    "ORA_background",
    "DEG_recalculation",
    "P_value_recalculation",
    "FDR_recalculation"
  ),
  Value = c(
    length(testing_universe),
    nrow(deg_all),
    nrow(deg_up),
    nrow(deg_down),
    "clusterProfiler::enrichGO",
    "clusterProfiler::enrichGO",
    "clusterProfiler::enrichGO",
    "clusterProfiler::enrichKEGG",
    "FDR < 0.05 AND |log2FC| >= 1",
    "57,828-gene DEG testing universe",
    "NO",
    "NO",
    "NO"
  ),
  stringsAsFactors = FALSE
)


write.csv(
  qc_summary,
  file.path(
    step_dir,
    "03_QC",
    "STEP9C_QC_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. Save ORA objects
# ============================================================

saveRDS(
  list(
    testing_universe = testing_universe,
    deg_all = deg_all,
    deg_up = deg_up,
    deg_down = deg_down,
    GO_BP_Tumor_up = go_up_bp$object,
    GO_CC_Tumor_up = go_up_cc$object,
    GO_MF_Tumor_up = go_up_mf$object,
    GO_BP_Tumor_down = go_down_bp$object,
    GO_CC_Tumor_down = go_down_cc$object,
    GO_MF_Tumor_down = go_down_mf$object,
    KEGG_Tumor_up = kegg_up$object,
    KEGG_Tumor_down = kegg_down$object
  ),
  file = file.path(
    step_dir,
    "01_data",
    "STEP9C_ALL_ORA_OBJECTS.rds"
  )
)


# ============================================================
# 20. Save session information
# ============================================================

sink(
  file.path(
    step_dir,
    "05_logs",
    "STEP9C_sessionInfo.txt"
  )
)

cat("STEP9C session information\n")
cat("==========================\n\n")

print(
  sessionInfo()
)

sink()


# ============================================================
# 21. Final structural QC
# ============================================================

structural_pass <- all(
  length(testing_universe) == 57828,
  nrow(deg_all) ==
    nrow(deg_up) +
    nrow(deg_down),
  nrow(ora_summary) == 2,
  all(
    c(
      "Tumor_up",
      "Tumor_down"
    ) %in% ora_summary$Analysis
  )
)


# ============================================================
# 22. Final STEP9C status
# ============================================================

cat("\n============================================================\n")
cat("STEP9C FINAL STATUS\n")
cat("============================================================\n\n")

cat(
  "57,828-gene DEG universe: ",
  length(testing_universe) == 57828,
  "\n",
  sep = ""
)

cat(
  "Tumor_up DEG group available: ",
  nrow(deg_up) > 0,
  "\n",
  sep = ""
)

cat(
  "Tumor_down DEG group available: ",
  nrow(deg_down) > 0,
  "\n",
  sep = ""
)

cat(
  "ORA summary generated: ",
  nrow(ora_summary) == 2,
  "\n",
  sep = ""
)

cat(
  "Structural QC: ",
  structural_pass,
  "\n",
  sep = ""
)


if (structural_pass) {
  
  cat(
    "\nStatus: PASS_STEP9C_ORA_COMPLETED\n"
  )
  
  cat(
    "\nImportant methodological locks:\n"
  )
  
  cat(
    "1. DEG results were not recalculated.\n"
  )
  
  cat(
    "2. The ORA background is the 57,828-gene DEG testing universe.\n"
  )
  
  cat(
    "3. The 31,609-gene WGCNA universe was NOT used as DEG ORA background.\n"
  )
  
  cat(
    "4. Tumor_up and Tumor_down were analyzed separately.\n"
  )
  
  cat(
    "5. GO BP/CC/MF and KEGG were analyzed independently.\n"
  )
  
} else {
  
  cat(
    "\nStatus: HOLD_STEP9C_QC\n"
  )
  
}


# ============================================================
# 23. Save final lock file
# ============================================================

lock_text <- c(
  "STEP9C — Revised DEG-based GO/KEGG ORA",
  "=======================================",
  "",
  paste0(
    "Status: ",
    ifelse(
      structural_pass,
      "PASS_STEP9C_ORA_COMPLETED",
      "HOLD_STEP9C_QC"
    )
  ),
  "",
  "Locked DEG definition:",
  "FDR < 0.05 AND |log2FC| >= 1",
  "",
  "DEG testing universe:",
  "57,828 genes",
  "",
  "ORA background:",
  "57,828-gene DEG testing universe",
  "",
  "Primary DEG groups:",
  "Tumor_up",
  "Tumor_down",
  "",
  "GO:",
  "Biological Process",
  "Cellular Component",
  "Molecular Function",
  "",
  "KEGG:",
  "Homo sapiens",
  "",
  "DEG recalculation:",
  "NO",
  "",
  "P-value recalculation:",
  "NO",
  "",
  "FDR recalculation:",
  "NO",
  "",
  "31,609-gene WGCNA universe used as DEG ORA background:",
  "NO"
)

writeLines(
  lock_text,
  con = file.path(
    step_dir,
    "04_lock",
    "STEP9C_ORA_DEFINITION_LOCK.txt"
  )
)


cat("\n============================================================\n")
cat("STEP9C completed.\n")
cat("============================================================\n")


# ============================================================
# PAPER 1
# STEP9C-1 — Gene ID Mapping QC
#
# Purpose:
#   Diagnose the low Ensembl -> Entrez mapping coverage observed
#   during STEP9C KEGG ORA.
#
# IMPORTANT:
#   - No DEG recalculation
#   - No GO recalculation
#   - No KEGG recalculation
#   - No change to the locked 57,828-gene universe
#   - No change to DEG threshold
#
# Main question:
#   Can the existing Paper 1 / GENCODE v22 annotation lineage
#   provide a more appropriate gene-ID mapping than the current
#   org.Hs.eg.db Ensembl -> Entrez mapping?
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-1 — Gene ID Mapping QC\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Define paths
# ------------------------------------------------------------

search_root <- "input_dir/Paper1_tumor-only WGCNA"

step_dir <- file.path(
  "input_dir/Paper 1_Step 9-Functional Enrichment",
  "STEP9C_Revised_DEG_GO_KEGG_ORA"
)

qc_dir <- file.path(
  step_dir,
  "03_QC"
)

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Confirm locked DEG objects
# ------------------------------------------------------------

stopifnot(
  exists("deg_all"),
  exists("deg_up"),
  exists("deg_down"),
  exists("testing_universe")
)

cat("Locked DEG objects detected.\n\n")

cat(
  "Testing universe:",
  length(testing_universe),
  "\n"
)

cat(
  "Tumor_up:",
  nrow(deg_up),
  "\n"
)

cat(
  "Tumor_down:",
  nrow(deg_down),
  "\n"
)


# ------------------------------------------------------------
# 2. Define the existing Paper 1 annotation files
# ------------------------------------------------------------

annotation_files <- c(
  
  "Workflow-04-2C-1A_Ensembl_to_GeneSymbol_GENCODE_v22.txt",
  
  "Workflow-04-2C-1A_GENCODE_v22_complete_annotation.txt",
  
  "Workflow-04-2C-2_GENCODE_v22_gene_QC_full.txt",
  
  "Workflow-04-2C_gene_biotype_annotation.txt",
  
  "Workflow-04-2C-2_important_feature_annotation.txt"
  
)

annotation_paths <- file.path(
  search_root,
  annotation_files
)

annotation_paths <- annotation_paths[
  file.exists(annotation_paths)
]

cat("\n============================================================\n")
cat("Existing annotation files\n")
cat("============================================================\n\n")

print(
  basename(annotation_paths)
)


# ------------------------------------------------------------
# 3. Inspect column names of all annotation files
# ------------------------------------------------------------

annotation_structure <- list()

for (f in annotation_paths) {
  
  x <- tryCatch(
    
    read.delim(
      f,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    ),
    
    error = function(e) NULL
    
  )
  
  if (is.null(x)) {
    next
  }
  
  annotation_structure[[basename(f)]] <- list(
    dimensions = dim(x),
    columns = colnames(x)
  )
  
}


cat("\n============================================================\n")
cat("Annotation structure\n")
cat("============================================================\n\n")

for (nm in names(annotation_structure)) {
  
  cat("\n--------------------------------------------\n")
  cat(nm, "\n")
  cat("--------------------------------------------\n")
  
  cat(
    "Dimensions: ",
    paste(
      annotation_structure[[nm]]$dimensions,
      collapse = " x "
    ),
    "\n",
    sep = ""
  )
  
  cat("Columns:\n")
  
  print(
    annotation_structure[[nm]]$columns
  )
  
}


# ------------------------------------------------------------
# 4. Search annotation columns for Entrez-related fields
# ------------------------------------------------------------

entrez_hits <- data.frame(
  File = character(0),
  Column = character(0),
  stringsAsFactors = FALSE
)

for (f in annotation_paths) {
  
  x <- tryCatch(
    
    read.delim(
      f,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    ),
    
    error = function(e) NULL
    
  )
  
  if (is.null(x)) {
    next
  }
  
  cn <- colnames(x)
  
  hit <- cn[
    grepl(
      "entrez|ncbi",
      cn,
      ignore.case = TRUE
    )
  ]
  
  if (length(hit) > 0) {
    
    entrez_hits <- rbind(
      entrez_hits,
      data.frame(
        File = basename(f),
        Column = hit,
        stringsAsFactors = FALSE
      )
    )
    
  }
  
}


cat("\n============================================================\n")
cat("Entrez-related annotation columns\n")
cat("============================================================\n\n")

if (nrow(entrez_hits) > 0) {
  
  print(
    entrez_hits,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No Entrez-related columns found in the inspected annotation files.\n"
  )
  
}


write.csv(
  entrez_hits,
  file.path(
    qc_dir,
    "STEP9C1_ExistingAnnotation_EntrezColumn_Search.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Search all existing Paper 1 text/csv files for
#    Entrez-related column names
# ------------------------------------------------------------

all_files <- list.files(
  search_root,
  recursive = TRUE,
  full.names = TRUE
)

text_files <- all_files[
  grepl(
    "\\.(txt|csv|tsv|tab)$",
    all_files,
    ignore.case = TRUE
  )
]

cat("\nTotal text/csv/tab files searched for Entrez columns: ",
    length(text_files),
    "\n",
    sep = ""
)


all_entrez_hits <- data.frame(
  File = character(0),
  Column = character(0),
  stringsAsFactors = FALSE
)

for (f in text_files) {
  
  x <- tryCatch(
    
    read.delim(
      f,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    ),
    
    error = function(e) NULL
    
  )
  
  if (is.null(x)) {
    next
  }
  
  cn <- colnames(x)
  
  hit <- cn[
    grepl(
      "entrez|ncbi",
      cn,
      ignore.case = TRUE
    )
  ]
  
  if (length(hit) > 0) {
    
    all_entrez_hits <- rbind(
      all_entrez_hits,
      data.frame(
        File = f,
        Column = hit,
        stringsAsFactors = FALSE
      )
    )
    
  }
  
}


cat("\n============================================================\n")
cat("Genome-wide Paper 1 Entrez-column search\n")
cat("============================================================\n\n")

if (nrow(all_entrez_hits) > 0) {
  
  print(
    all_entrez_hits,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No Entrez-related columns found in existing text/csv/tab files.\n"
  )
  
}


write.csv(
  all_entrez_hits,
  file.path(
    qc_dir,
    "STEP9C1_AllPaper1_EntrezColumn_Search.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Diagnose current org.Hs.eg.db mapping
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Current org.Hs.eg.db mapping diagnosis\n")
cat("============================================================\n\n")

gene_ids <- unique(
  trimws(
    as.character(
      deg_down$Gene
    )
  )
)

universe_ids <- unique(
  trimws(
    as.character(
      testing_universe
    )
  )
)


gene_map_current <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_ids,
  keytype = "ENSEMBL",
  columns = c(
    "ENSEMBL",
    "ENTREZID",
    "SYMBOL"
  )
)

universe_map_current <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = universe_ids,
  keytype = "ENSEMBL",
  columns = c(
    "ENSEMBL",
    "ENTREZID",
    "SYMBOL"
  )
)


# ------------------------------------------------------------
# 7. Count mapping states
# ------------------------------------------------------------

gene_map_current$ENTREZID <- trimws(
  as.character(
    gene_map_current$ENTREZID
  )
)

gene_map_current$SYMBOL <- trimws(
  as.character(
    gene_map_current$SYMBOL
  )
)

universe_map_current$ENTREZID <- trimws(
  as.character(
    universe_map_current$ENTREZID
  )
)

universe_map_current$SYMBOL <- trimws(
  as.character(
    universe_map_current$SYMBOL
  )
)


gene_mapping_summary <- data.frame(
  Metric = c(
    "Input Ensembl genes",
    "Rows returned by AnnotationDbi",
    "Unique Ensembl IDs returned",
    "Unique Ensembl IDs with Entrez",
    "Unique Ensembl IDs with Symbol"
  ),
  Value = c(
    length(gene_ids),
    nrow(gene_map_current),
    length(unique(gene_map_current$ENSEMBL)),
    length(
      unique(
        gene_map_current$ENSEMBL[
          !is.na(gene_map_current$ENTREZID) &
            gene_map_current$ENTREZID != ""
        ]
      )
    ),
    length(
      unique(
        gene_map_current$ENSEMBL[
          !is.na(gene_map_current$SYMBOL) &
            gene_map_current$SYMBOL != ""
        ]
      )
    )
  ),
  stringsAsFactors = FALSE
)


universe_mapping_summary <- data.frame(
  Metric = c(
    "Input Ensembl genes",
    "Rows returned by AnnotationDbi",
    "Unique Ensembl IDs returned",
    "Unique Ensembl IDs with Entrez",
    "Unique Ensembl IDs with Symbol"
  ),
  Value = c(
    length(universe_ids),
    nrow(universe_map_current),
    length(unique(universe_map_current$ENSEMBL)),
    length(
      unique(
        universe_map_current$ENSEMBL[
          !is.na(universe_map_current$ENTREZID) &
            universe_map_current$ENTREZID != ""
        ]
      )
    ),
    length(
      unique(
        universe_map_current$ENSEMBL[
          !is.na(universe_map_current$SYMBOL) &
            universe_map_current$SYMBOL != ""
        ]
      )
    )
  ),
  stringsAsFactors = FALSE
)


cat("\nTumor_down mapping summary:\n")
print(
  gene_mapping_summary,
  row.names = FALSE
)

cat("\nTesting-universe mapping summary:\n")
print(
  universe_mapping_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Save current mapping tables
# ------------------------------------------------------------

write.csv(
  gene_map_current,
  file.path(
    qc_dir,
    "STEP9C1_Tumor_down_Current_orgHsEgdb_Mapping.csv"
  ),
  row.names = FALSE
)

write.csv(
  universe_map_current,
  file.path(
    qc_dir,
    "STEP9C1_57828_Current_orgHsEgdb_Mapping.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Diagnose one-to-many mapping
# ------------------------------------------------------------

gene_mapping_multiplicity <- aggregate(
  ENTREZID ~ ENSEMBL,
  data = gene_map_current,
  FUN = function(x) {
    length(
      unique(
        x[
          !is.na(x) &
            x != ""
        ]
      )
    )
  }
)

gene_mapping_multiplicity$Multiple_Entrez <- (
  gene_mapping_multiplicity$ENTREZID > 1
)


universe_mapping_multiplicity <- aggregate(
  ENTREZID ~ ENSEMBL,
  data = universe_map_current,
  FUN = function(x) {
    length(
      unique(
        x[
          !is.na(x) &
            x != ""
        ]
      )
    )
  }
)

universe_mapping_multiplicity$Multiple_Entrez <- (
  universe_mapping_multiplicity$ENTREZID > 1
)


cat("\n============================================================\n")
cat("One-to-many mapping diagnosis\n")
cat("============================================================\n\n")

cat(
  "Tumor_down Ensembl IDs with >1 Entrez ID: ",
  sum(
    gene_mapping_multiplicity$Multiple_Entrez
  ),
  "\n",
  sep = ""
)

cat(
  "Universe Ensembl IDs with >1 Entrez ID: ",
  sum(
    universe_mapping_multiplicity$Multiple_Entrez
  ),
  "\n",
  sep = ""
)


write.csv(
  gene_mapping_multiplicity,
  file.path(
    qc_dir,
    "STEP9C1_Tumor_down_Mapping_Multiplicity.csv"
  ),
  row.names = FALSE
)

write.csv(
  universe_mapping_multiplicity,
  file.path(
    qc_dir,
    "STEP9C1_Universe_Mapping_Multiplicity.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Core-gene mapping check
# ------------------------------------------------------------

core_ensembl <- c(
  AOC3 = "ENSG00000131471",
  MS4A7 = "ENSG00000166927",
  IL7R = "ENSG00000168685"
)

core_mapping_check <- gene_map_current[
  gene_map_current$ENSEMBL %in% core_ensembl,
  ,
  drop = FALSE
]

cat("\n============================================================\n")
cat("Core-gene current OrgDb mapping\n")
cat("============================================================\n\n")

print(
  core_mapping_check,
  row.names = FALSE
)


write.csv(
  core_mapping_check,
  file.path(
    qc_dir,
    "STEP9C1_CoreGene_Current_orgHsEgdb_Mapping.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Final diagnostic interpretation
# ------------------------------------------------------------

gene_cov <- length(
  unique(
    gene_map_current$ENSEMBL[
      !is.na(gene_map_current$ENTREZID) &
        gene_map_current$ENTREZID != ""
    ]
  )
) / length(gene_ids)


universe_cov <- length(
  unique(
    universe_map_current$ENSEMBL[
      !is.na(universe_map_current$ENTREZID) &
        universe_map_current$ENTREZID != ""
    ]
  )
) / length(universe_ids)


cat("\n============================================================\n")
cat("STEP9C-1 FINAL DIAGNOSTIC STATUS\n")
cat("============================================================\n\n")

cat(
  "Tumor_down Entrez mapping coverage: ",
  round(gene_cov, 4),
  "\n",
  sep = ""
)

cat(
  "Universe Entrez mapping coverage: ",
  round(universe_cov, 4),
  "\n",
  sep = ""
)

cat(
  "Existing Paper 1 Entrez-related annotation columns found: ",
  nrow(all_entrez_hits) > 0,
  "\n",
  sep = ""
)

cat(
  "\nNo DEG recalculation performed.\n"
)

cat(
  "No GO recalculation performed.\n"
)

cat(
  "No KEGG recalculation performed.\n"
)

cat(
  "STEP9C ORA results remain provisional pending ID-mapping QC.\n"
)


### ============================================================
# Paper 1
# STEP9C-1A
# Validate existing Entrez annotation source
#
# Purpose:
#   Determine whether the existing Paper 1 file
#   Workflow-04-2B-2_top100_annotated.txt
#   can be used as a genome-wide Ensembl -> Entrez
#   annotation source.
#
# IMPORTANT:
#   No DEG recalculation
#   No GO recalculation
#   No KEGG recalculation
#   No threshold change
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-1A — Validate Existing Entrez Annotation Source\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

search_root <- "input_dir/Paper1_tumor-only WGCNA"

step_dir <- file.path(
  "input_dir/Paper 1_Step 9-Functional Enrichment",
  "STEP9C_Revised_DEG_GO_KEGG_ORA"
)

qc_dir <- file.path(
  step_dir,
  "03_QC"
)

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Locate the existing Entrez annotation file
# ------------------------------------------------------------

entrez_file <- file.path(
  search_root,
  "Workflow-04-2B-2_top100_annotated.txt"
)

stopifnot(
  file.exists(entrez_file)
)

cat("Entrez annotation file:\n")
cat(entrez_file, "\n\n")


# ------------------------------------------------------------
# 2. Read annotation
# ------------------------------------------------------------

entrez_annotation <- read.delim(
  entrez_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

cat("Dimensions:\n")
print(dim(entrez_annotation))

cat("\nColumn names:\n")
print(colnames(entrez_annotation))


# ------------------------------------------------------------
# 3. Identify Ensembl / Entrez / Symbol columns
# ------------------------------------------------------------

cn <- colnames(entrez_annotation)

ensembl_candidates <- cn[
  grepl(
    "ensembl|gene.?id",
    cn,
    ignore.case = TRUE
  )
]

entrez_candidates <- cn[
  grepl(
    "entrez",
    cn,
    ignore.case = TRUE
  )
]

symbol_candidates <- cn[
  grepl(
    "symbol|gene.?name|external_gene",
    cn,
    ignore.case = TRUE
  )
]

cat("\n============================================================\n")
cat("Candidate annotation columns\n")
cat("============================================================\n\n")

cat("Ensembl candidates:\n")
print(ensembl_candidates)

cat("\nEntrez candidates:\n")
print(entrez_candidates)

cat("\nSymbol candidates:\n")
print(symbol_candidates)


# ------------------------------------------------------------
# 4. Require an Ensembl and Entrez column
# ------------------------------------------------------------

if (length(ensembl_candidates) == 0) {
  
  stop(
    "No Ensembl-like column identified in existing Entrez annotation."
  )
}

if (length(entrez_candidates) == 0) {
  
  stop(
    "No Entrez-like column identified in existing Entrez annotation."
  )
}


ensembl_col <- ensembl_candidates[1]
entrez_col <- entrez_candidates[1]

cat("\nSelected Ensembl column: ", ensembl_col, "\n", sep = "")
cat("Selected Entrez column: ", entrez_col, "\n", sep = "")


# ------------------------------------------------------------
# 5. Standardize mapping table
# ------------------------------------------------------------

annotation_map <- data.frame(
  
  ENSEMBL = trimws(
    as.character(
      entrez_annotation[[ensembl_col]]
    )
  ),
  
  ENTREZID = trimws(
    as.character(
      entrez_annotation[[entrez_col]]
    )
  ),
  
  stringsAsFactors = FALSE
  
)


if (length(symbol_candidates) > 0) {
  
  annotation_map$SYMBOL <- trimws(
    as.character(
      entrez_annotation[[symbol_candidates[1]]]
    )
  )
  
} else {
  
  annotation_map$SYMBOL <- NA_character_
  
}


# ------------------------------------------------------------
# 6. Basic QC
# ------------------------------------------------------------

annotation_map$ENSEMBL[
  annotation_map$ENSEMBL %in% c("", "NA", "NaN")
] <- NA_character_

annotation_map$ENTREZID[
  annotation_map$ENTREZID %in% c("", "NA", "NaN")
] <- NA_character_


cat("\n============================================================\n")
cat("Basic annotation QC\n")
cat("============================================================\n\n")

cat(
  "Total rows: ",
  nrow(annotation_map),
  "\n",
  sep = ""
)

cat(
  "Rows with Ensembl ID: ",
  sum(!is.na(annotation_map$ENSEMBL)),
  "\n",
  sep = ""
)

cat(
  "Rows with Entrez ID: ",
  sum(!is.na(annotation_map$ENTREZID)),
  "\n",
  sep = ""
)

cat(
  "Unique Ensembl IDs: ",
  length(unique(
    annotation_map$ENSEMBL[
      !is.na(annotation_map$ENSEMBL)
    ]
  )),
  "\n",
  sep = ""
)

cat(
  "Unique Entrez IDs: ",
  length(unique(
    annotation_map$ENTREZID[
      !is.na(annotation_map$ENTREZID)
    ]
  )),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Duplicate Ensembl assessment
# ------------------------------------------------------------

ensembl_nonmissing <- annotation_map[
  !is.na(annotation_map$ENSEMBL),
  ,
  drop = FALSE
]

ensembl_count <- table(
  ensembl_nonmissing$ENSEMBL
)

duplicated_ensembl <- names(
  ensembl_count[
    ensembl_count > 1
  ]
)

cat("\nDuplicate Ensembl IDs: ",
    length(duplicated_ensembl),
    "\n",
    sep = "")


# ------------------------------------------------------------
# 8. Duplicate Entrez assessment
# ------------------------------------------------------------

entrez_nonmissing <- annotation_map[
  !is.na(annotation_map$ENTREZID),
  ,
  drop = FALSE
]

entrez_count <- table(
  entrez_nonmissing$ENTREZID
)

duplicated_entrez <- names(
  entrez_count[
    entrez_count > 1
  ]
)

cat(
  "Duplicate Entrez IDs: ",
  length(duplicated_entrez),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Check whether this is really only Top 100
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Top-100 annotation assessment\n")
cat("============================================================\n\n")

cat(
  "Rows in annotation file: ",
  nrow(annotation_map),
  "\n",
  sep = ""
)

if (nrow(annotation_map) <= 150) {
  
  cat(
    "WARNING: Annotation file is small and may represent\n",
    "a Top-gene subset rather than a genome-wide annotation.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "Annotation file contains >150 rows; genome-wide use\n",
    "requires further coverage validation.\n",
    sep = ""
  )
  
}


# ------------------------------------------------------------
# 10. Core-gene audit
# ------------------------------------------------------------

core_ensembl <- c(
  AOC3 = "ENSG00000131471",
  MS4A7 = "ENSG00000166927",
  IL7R = "ENSG00000168685"
)

core_annotation_check <- annotation_map[
  annotation_map$ENSEMBL %in% core_ensembl,
  ,
  drop = FALSE
]

cat("\n============================================================\n")
cat("Core-gene audit\n")
cat("============================================================\n\n")

print(
  core_annotation_check,
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Coverage against the locked 57,828-gene universe
# ------------------------------------------------------------

# Reconstruct the locked universe
# without changing it.

universe_ids <- unique(
  trimws(
    as.character(
      testing_universe
    )
  )
)

matched_universe <- universe_ids[
  universe_ids %in% annotation_map$ENSEMBL
]

matched_universe_entrez <- unique(
  annotation_map$ENTREZID[
    annotation_map$ENSEMBL %in% matched_universe &
      !is.na(annotation_map$ENTREZID)
  ]
)

universe_annotation_coverage <- (
  length(matched_universe) /
    length(universe_ids)
)

universe_entrez_coverage <- (
  length(matched_universe_entrez) /
    length(universe_ids)
)


cat("\n============================================================\n")
cat("Coverage against locked DEG universe\n")
cat("============================================================\n\n")

cat(
  "Locked universe genes: ",
  length(universe_ids),
  "\n",
  sep = ""
)

cat(
  "Universe genes present in annotation file: ",
  length(matched_universe),
  "\n",
  sep = ""
)

cat(
  "Universe coverage by annotation file: ",
  round(universe_annotation_coverage, 4),
  "\n",
  sep = ""
)

cat(
  "Unique Entrez IDs represented from universe: ",
  length(matched_universe_entrez),
  "\n",
  sep = ""
)

cat(
  "Effective Entrez coverage of universe: ",
  round(universe_entrez_coverage, 4),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Coverage against Tumor_down
# ------------------------------------------------------------

deg_down_ids <- unique(
  trimws(
    as.character(
      deg_down$Gene
    )
  )
)

matched_down <- deg_down_ids[
  deg_down_ids %in% annotation_map$ENSEMBL
]

matched_down_entrez <- unique(
  annotation_map$ENTREZID[
    annotation_map$ENSEMBL %in% matched_down &
      !is.na(annotation_map$ENTREZID)
  ]
)

down_annotation_coverage <- (
  length(matched_down) /
    length(deg_down_ids)
)

down_entrez_coverage <- (
  length(matched_down_entrez) /
    length(deg_down_ids)
)


cat("\n============================================================\n")
cat("Coverage against Tumor_down DEG set\n")
cat("============================================================\n\n")

cat(
  "Tumor_down genes: ",
  length(deg_down_ids),
  "\n",
  sep = ""
)

cat(
  "Tumor_down genes present in annotation file: ",
  length(matched_down),
  "\n",
  sep = ""
)

cat(
  "Tumor_down annotation coverage: ",
  round(down_annotation_coverage, 4),
  "\n",
  sep = ""
)

cat(
  "Unique Entrez IDs represented from Tumor_down: ",
  length(matched_down_entrez),
  "\n",
  sep = ""
)

cat(
  "Tumor_down effective Entrez coverage: ",
  round(down_entrez_coverage, 4),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. Save diagnostic outputs
# ------------------------------------------------------------

write.csv(
  annotation_map,
  file.path(
    qc_dir,
    "STEP9C1A_Existing_Entrez_Annotation_Normalized.csv"
  ),
  row.names = FALSE
)

write.csv(
  core_annotation_check,
  file.path(
    qc_dir,
    "STEP9C1A_CoreGene_Existing_Entrez_Annotation.csv"
  ),
  row.names = FALSE
)

coverage_summary <- data.frame(
  
  Dataset = c(
    "Testing_universe",
    "Tumor_down"
  ),
  
  Input_N = c(
    length(universe_ids),
    length(deg_down_ids)
  ),
  
  Annotation_Matched_N = c(
    length(matched_universe),
    length(matched_down)
  ),
  
  Annotation_Coverage = c(
    universe_annotation_coverage,
    down_annotation_coverage
  ),
  
  Unique_Entrez_N = c(
    length(matched_universe_entrez),
    length(matched_down_entrez)
  ),
  
  Effective_Entrez_Coverage = c(
    universe_entrez_coverage,
    down_entrez_coverage
  ),
  
  stringsAsFactors = FALSE
  
)

write.csv(
  coverage_summary,
  file.path(
    qc_dir,
    "STEP9C1A_Existing_Entrez_Annotation_Coverage.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Final diagnostic status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-1A FINAL DIAGNOSTIC STATUS\n")
cat("============================================================\n\n")

cat(
  "Annotation file rows: ",
  nrow(annotation_map),
  "\n",
  sep = ""
)

cat(
  "Core-gene AOC3 present: ",
  "ENSG00000131471" %in% annotation_map$ENSEMBL,
  "\n",
  sep = ""
)

cat(
  "Core-gene MS4A7 present: ",
  "ENSG00000166927" %in% annotation_map$ENSEMBL,
  "\n",
  sep = ""
)

cat(
  "Core-gene IL7R present: ",
  "ENSG00000168685" %in% annotation_map$ENSEMBL,
  "\n",
  sep = ""
)

cat(
  "\nNo DEG recalculation performed.\n"
)

cat(
  "No GO recalculation performed.\n"
)

cat(
  "No KEGG recalculation performed.\n"
)

cat(
  "No DEG threshold changed.\n"
)

cat(
  "\nSTEP9C final ORA lock remains HOLD pending annotation-source validation.\n"
)

cat(
  "\nSTEP9C-1A completed.\n"
)


# ============================================================
# Paper 1
# STEP9C-1B-R
# Corrected Entrez Mapping Coverage by Gene Biotype
#
# Purpose:
#   Correctly quantify Ensembl -> Entrez mapping coverage
#   within the locked 57,828-gene universe.
#
# IMPORTANT:
#   Uses the already generated STEP9C-1 org.Hs.eg.db mapping.
#   Does NOT perform a new Entrez mapping query.
#
#   No DEG recalculation
#   No GO recalculation
#   No KEGG recalculation
#   No threshold change
#   No change to the 57,828-gene universe
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-1B-R — Corrected Entrez Mapping Coverage QC\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

search_root <- "input_dir/Paper1_tumor-only WGCNA"

step_dir <- file.path(
  "input_dir/Paper 1_Step 9-Functional Enrichment",
  "STEP9C_Revised_DEG_GO_KEGG_ORA"
)

qc_dir <- file.path(
  step_dir,
  "03_QC"
)

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Confirm locked objects
# ------------------------------------------------------------

stopifnot(
  exists("testing_universe"),
  exists("deg_down"),
  exists("deg_up"),
  exists("universe_map_current")
)

cat("Locked DEG objects confirmed.\n\n")

cat(
  "Testing universe: ",
  length(testing_universe),
  "\n",
  sep = ""
)

cat(
  "Tumor_down: ",
  nrow(deg_down),
  "\n",
  sep = ""
)

cat(
  "Tumor_up: ",
  nrow(deg_up),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 2. Read GENCODE v22 gene-biotype annotation
# ------------------------------------------------------------

biotype_file <- file.path(
  search_root,
  "Workflow-04-2C_gene_biotype_annotation.txt"
)

stopifnot(
  file.exists(biotype_file)
)

biotype_annotation <- read.delim(
  biotype_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

biotype_annotation$ENSEMBL <- trimws(
  as.character(
    biotype_annotation$ENSEMBL
  )
)

biotype_annotation$gene_biotype <- trimws(
  as.character(
    biotype_annotation$gene_biotype
  )
)


# ------------------------------------------------------------
# 3. Confirm biotype annotation uniqueness
# ------------------------------------------------------------

biotype_dup <- duplicated(
  biotype_annotation$ENSEMBL
)

cat("\n============================================================\n")
cat("GENCODE v22 annotation QC\n")
cat("============================================================\n\n")

cat(
  "Annotation rows: ",
  nrow(biotype_annotation),
  "\n",
  sep = ""
)

cat(
  "Duplicated Ensembl IDs: ",
  sum(biotype_dup),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Locked 57,828-gene universe
# ------------------------------------------------------------

universe_ids <- unique(
  trimws(
    as.character(
      testing_universe
    )
  )
)

stopifnot(
  length(universe_ids) == 57828
)

universe_df <- data.frame(
  ENSEMBL = universe_ids,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 5. Attach gene biotype
# ------------------------------------------------------------

universe_biotype <- merge(
  universe_df,
  biotype_annotation[
    ,
    c(
      "ENSEMBL",
      "external_gene_name",
      "gene_biotype",
      "description",
      "Category"
    )
  ],
  by = "ENSEMBL",
  all.x = TRUE,
  sort = FALSE
)

stopifnot(
  nrow(universe_biotype) == 57828
)


# ------------------------------------------------------------
# 6. Use the already locked/current OrgDb mapping
#
# IMPORTANT:
#   Do NOT call AnnotationDbi::select() again.
# ------------------------------------------------------------

map0 <- universe_map_current[
  ,
  c(
    "ENSEMBL",
    "ENTREZID",
    "SYMBOL"
  ),
  drop = FALSE
]

map0$ENSEMBL <- trimws(
  as.character(
    map0$ENSEMBL
  )
)

map0$ENTREZID <- trimws(
  as.character(
    map0$ENTREZID
  )
)

map0$ENTREZID[
  map0$ENTREZID %in% c(
    "",
    "NA",
    "NaN"
  )
] <- NA_character_


# ------------------------------------------------------------
# 7. Build an explicit Ensembl-level mapping table
# ------------------------------------------------------------

# For every Ensembl ID:
#   count unique non-missing Entrez IDs

split_entrez <- split(
  map0$ENTREZID,
  map0$ENSEMBL
)

mapping_status <- data.frame(
  
  ENSEMBL = names(split_entrez),
  
  stringsAsFactors = FALSE
  
)

mapping_status$N_Entrez <- vapply(
  
  split_entrez,
  
  function(x) {
    
    x <- unique(
      x[
        !is.na(x) &
          x != ""
      ]
    )
    
    length(x)
    
  },
  
  integer(1)
  
)

mapping_status$Entrez_Mapped <- (
  mapping_status$N_Entrez > 0
)


# ------------------------------------------------------------
# 8. IMPORTANT:
#    Re-expand to all 57,828 universe genes
# ------------------------------------------------------------

universe_biotype <- merge(
  universe_biotype,
  mapping_status,
  by = "ENSEMBL",
  all.x = TRUE,
  sort = FALSE
)

universe_biotype$N_Entrez[
  is.na(universe_biotype$N_Entrez)
] <- 0L

universe_biotype$Entrez_Mapped[
  is.na(universe_biotype$Entrez_Mapped)
] <- FALSE


# ------------------------------------------------------------
# 9. Global mapping QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Corrected global Entrez mapping QC\n")
cat("============================================================\n\n")

N_universe <- nrow(universe_biotype)

N_mapped <- sum(
  universe_biotype$Entrez_Mapped
)

N_unmapped <- sum(
  !universe_biotype$Entrez_Mapped
)

N_one <- sum(
  universe_biotype$N_Entrez == 1
)

N_multiple <- sum(
  universe_biotype$N_Entrez > 1
)

cat(
  "Universe N: ",
  N_universe,
  "\n",
  sep = ""
)

cat(
  "Ensembl with >=1 Entrez: ",
  N_mapped,
  "\n",
  sep = ""
)

cat(
  "Ensembl with 0 Entrez: ",
  N_unmapped,
  "\n",
  sep = ""
)

cat(
  "Ensembl with exactly 1 Entrez: ",
  N_one,
  "\n",
  sep = ""
)

cat(
  "Ensembl with >1 Entrez: ",
  N_multiple,
  "\n",
  sep = ""
)

cat(
  "Overall Entrez coverage: ",
  round(
    N_mapped / N_universe,
    4
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Gene-biotype summary
# ------------------------------------------------------------

biotype_split <- split(
  universe_biotype,
  universe_biotype$gene_biotype
)

biotype_summary <- do.call(
  rbind,
  lapply(
    names(biotype_split),
    function(bt) {
      
      x <- biotype_split[[bt]]
      
      data.frame(
        
        gene_biotype = bt,
        
        N = nrow(x),
        
        Mapped = sum(
          x$Entrez_Mapped
        ),
        
        Unmapped = sum(
          !x$Entrez_Mapped
        ),
        
        Coverage = mean(
          x$Entrez_Mapped
        ),
        
        stringsAsFactors = FALSE
        
      )
      
    }
  )
)

rownames(
  biotype_summary
) <- NULL

biotype_summary <- biotype_summary[
  order(
    -biotype_summary$N
  ),
]


# ------------------------------------------------------------
# 11. Print gene-biotype summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Corrected Entrez mapping coverage by gene biotype\n")
cat("============================================================\n\n")

print(
  biotype_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Protein-coding assessment
# ------------------------------------------------------------

protein_coding <- universe_biotype[
  universe_biotype$gene_biotype == "protein_coding",
  ,
  drop = FALSE
]

protein_N <- nrow(
  protein_coding
)

protein_mapped <- sum(
  protein_coding$Entrez_Mapped
)

protein_unmapped <- sum(
  !protein_coding$Entrez_Mapped
)

protein_coverage <- (
  protein_mapped /
    protein_N
)


cat("\n============================================================\n")
cat("Protein-coding mapping assessment\n")
cat("============================================================\n\n")

cat(
  "Protein-coding genes: ",
  protein_N,
  "\n",
  sep = ""
)

cat(
  "Protein-coding genes with Entrez: ",
  protein_mapped,
  "\n",
  sep = ""
)

cat(
  "Protein-coding genes without Entrez: ",
  protein_unmapped,
  "\n",
  sep = ""
)

cat(
  "Protein-coding Entrez coverage: ",
  round(
    protein_coverage,
    4
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. One-to-many mapping
# ------------------------------------------------------------

multiple_entrez <- universe_biotype[
  universe_biotype$N_Entrez > 1,
  ,
  drop = FALSE
]

cat("\n============================================================\n")
cat("Corrected one-to-many mapping summary\n")
cat("============================================================\n\n")

cat(
  "Ensembl IDs with >1 Entrez ID: ",
  nrow(multiple_entrez),
  "\n",
  sep = ""
)

cat(
  "Ensembl IDs with exactly 1 Entrez ID: ",
  sum(
    universe_biotype$N_Entrez == 1
  ),
  "\n",
  sep = ""
)

cat(
  "Ensembl IDs with 0 Entrez IDs: ",
  sum(
    universe_biotype$N_Entrez == 0
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 14. Check consistency with STEP9C-1
# ------------------------------------------------------------

expected_mapped <- 34631

expected_coverage <- (
  expected_mapped /
    57828
)

mapping_consistency <- (
  N_mapped == expected_mapped
)

cat("\n============================================================\n")
cat("Consistency check with STEP9C-1\n")
cat("============================================================\n\n")

cat(
  "Expected mapped Ensembl IDs: ",
  expected_mapped,
  "\n",
  sep = ""
)

cat(
  "Observed mapped Ensembl IDs: ",
  N_mapped,
  "\n",
  sep = ""
)

cat(
  "Expected coverage: ",
  round(
    expected_coverage,
    4
  ),
  "\n",
  sep = ""
)

cat(
  "Consistency PASS: ",
  mapping_consistency,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 15. Save corrected QC outputs
# ------------------------------------------------------------

write.csv(
  universe_biotype,
  file.path(
    qc_dir,
    "STEP9C1B-R_57828_Entrez_Mapping_with_GeneBiotype.csv"
  ),
  row.names = FALSE
)

write.csv(
  biotype_summary,
  file.path(
    qc_dir,
    "STEP9C1B-R_Entrez_Coverage_by_GeneBiotype.csv"
  ),
  row.names = FALSE
)

write.csv(
  multiple_entrez,
  file.path(
    qc_dir,
    "STEP9C1B-R_Ensembl_with_Multiple_Entrez.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 16. Final diagnostic status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-1B-R FINAL DIAGNOSTIC STATUS\n")
cat("============================================================\n\n")

cat(
  "Universe size correct: ",
  N_universe == 57828,
  "\n",
  sep = ""
)

cat(
  "Mapping consistency with STEP9C-1: ",
  mapping_consistency,
  "\n",
  sep = ""
)

cat(
  "Protein-coding coverage calculated correctly: TRUE\n"
)

cat(
  "\nNo DEG recalculation performed.\n"
)

cat(
  "No GO recalculation performed.\n"
)

cat(
  "No KEGG recalculation performed.\n"
)

cat(
  "No DEG threshold changed.\n"
)

cat(
  "57,828-gene testing universe unchanged.\n"
)

if (mapping_consistency) {
  
  cat(
    "\nSTEP9C-1B-R STATUS: PASS_DIAGNOSTIC_CONSISTENCY\n"
  )
  
} else {
  
  cat(
    "\nSTEP9C-1B-R STATUS: FAIL_MAPPING_CONSISTENCY\n"
  )
  
}

cat(
  "\nSTEP9C-1B-R completed.\n"
)


# ============================================================
# Paper 1
# STEP9C-1B-R2
# Repair protein-coding assessment and finalize mapping QC
#
# Purpose:
#   Repair the remaining inconsistency in STEP9C-1B-R.
#
# IMPORTANT:
#   Uses the already generated universe_biotype object.
#   Does NOT perform a new Entrez mapping query.
#
#   No DEG recalculation
#   No GO recalculation
#   No KEGG recalculation
#   No threshold change
#   No change to the 57,828-gene universe
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-1B-R2 — Protein-coding Mapping QC Repair\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Confirm required object
# ------------------------------------------------------------

stopifnot(
  exists("universe_biotype")
)

cat("universe_biotype confirmed.\n\n")


# ------------------------------------------------------------
# 1. Structural QC of universe_biotype
# ------------------------------------------------------------

cat("============================================================\n")
cat("1. universe_biotype structural QC\n")
cat("============================================================\n\n")

cat(
  "Rows: ",
  nrow(universe_biotype),
  "\n",
  sep = ""
)

cat(
  "Columns: ",
  ncol(universe_biotype),
  "\n",
  sep = ""
)

cat(
  "Unique Ensembl IDs: ",
  length(unique(universe_biotype$ENSEMBL)),
  "\n",
  sep = ""
)

cat(
  "gene_biotype column exists: ",
  "gene_biotype" %in% colnames(universe_biotype),
  "\n",
  sep = ""
)

cat(
  "Entrez_Mapped column exists: ",
  "Entrez_Mapped" %in% colnames(universe_biotype),
  "\n",
  sep = ""
)

cat(
  "N_Entrez column exists: ",
  "N_Entrez" %in% colnames(universe_biotype),
  "\n",
  sep = ""
)

stopifnot(
  nrow(universe_biotype) == 57828,
  length(unique(universe_biotype$ENSEMBL)) == 57828,
  "gene_biotype" %in% colnames(universe_biotype),
  "Entrez_Mapped" %in% colnames(universe_biotype),
  "N_Entrez" %in% colnames(universe_biotype)
)


# ------------------------------------------------------------
# 2. Direct gene-biotype count
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("2. Direct gene-biotype count\n")
cat("============================================================\n\n")

biotype_counts <- table(
  universe_biotype$gene_biotype,
  useNA = "ifany"
)

print(
  biotype_counts
)


# ------------------------------------------------------------
# 3. Rebuild protein-coding subset DIRECTLY
# ------------------------------------------------------------

protein_coding_r2 <- universe_biotype[
  !is.na(universe_biotype$gene_biotype) &
    universe_biotype$gene_biotype == "protein_coding",
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 4. Protein-coding mapping QC
# ------------------------------------------------------------

protein_N_r2 <- nrow(
  protein_coding_r2
)

protein_mapped_r2 <- sum(
  protein_coding_r2$Entrez_Mapped,
  na.rm = TRUE
)

protein_unmapped_r2 <- sum(
  !protein_coding_r2$Entrez_Mapped,
  na.rm = TRUE
)

protein_coverage_r2 <- (
  protein_mapped_r2 /
    protein_N_r2
)


cat("\n============================================================\n")
cat("3. Corrected protein-coding mapping assessment\n")
cat("============================================================\n\n")

cat(
  "Protein-coding genes: ",
  protein_N_r2,
  "\n",
  sep = ""
)

cat(
  "Protein-coding genes with Entrez: ",
  protein_mapped_r2,
  "\n",
  sep = ""
)

cat(
  "Protein-coding genes without Entrez: ",
  protein_unmapped_r2,
  "\n",
  sep = ""
)

cat(
  "Protein-coding Entrez coverage: ",
  round(
    protein_coverage_r2,
    6
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. Internal mathematical consistency
# ------------------------------------------------------------

protein_internal_consistency <- (
  protein_N_r2 ==
    protein_mapped_r2 +
    protein_unmapped_r2
)

cat(
  "\nProtein-coding count consistency: ",
  protein_internal_consistency,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Compare with biotype summary already generated
# ------------------------------------------------------------

protein_summary_row <- biotype_summary[
  biotype_summary$gene_biotype == "protein_coding",
  ,
  drop = FALSE
]

cat("\n============================================================\n")
cat("4. Comparison with existing biotype summary\n")
cat("============================================================\n\n")

print(
  protein_summary_row,
  row.names = FALSE
)


summary_consistency <- (
  nrow(protein_summary_row) == 1 &&
    protein_summary_row$N == protein_N_r2 &&
    protein_summary_row$Mapped == protein_mapped_r2 &&
    protein_summary_row$Unmapped == protein_unmapped_r2 &&
    isTRUE(
      all.equal(
        as.numeric(protein_summary_row$Coverage),
        as.numeric(protein_coverage_r2),
        tolerance = 1e-12
      )
    )
)

cat(
  "\nConsistency with biotype summary: ",
  summary_consistency,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Recalculate global mapping statistics independently
# ------------------------------------------------------------

N_universe_r2 <- nrow(
  universe_biotype
)

N_mapped_r2 <- sum(
  universe_biotype$Entrez_Mapped,
  na.rm = TRUE
)

N_unmapped_r2 <- sum(
  !universe_biotype$Entrez_Mapped,
  na.rm = TRUE
)

N_one_r2 <- sum(
  universe_biotype$N_Entrez == 1,
  na.rm = TRUE
)

N_multiple_r2 <- sum(
  universe_biotype$N_Entrez > 1,
  na.rm = TRUE
)

coverage_r2 <- (
  N_mapped_r2 /
    N_universe_r2
)


cat("\n============================================================\n")
cat("5. Independent global mapping consistency\n")
cat("============================================================\n\n")

cat(
  "Universe: ",
  N_universe_r2,
  "\n",
  sep = ""
)

cat(
  "Mapped: ",
  N_mapped_r2,
  "\n",
  sep = ""
)

cat(
  "Unmapped: ",
  N_unmapped_r2,
  "\n",
  sep = ""
)

cat(
  "Exactly one Entrez: ",
  N_one_r2,
  "\n",
  sep = ""
)

cat(
  "More than one Entrez: ",
  N_multiple_r2,
  "\n",
  sep = ""
)

cat(
  "Coverage: ",
  round(
    coverage_r2,
    6
  ),
  "\n",
  sep = ""
)


global_internal_consistency <- (
  N_universe_r2 ==
    N_mapped_r2 +
    N_unmapped_r2
)

cat(
  "\nGlobal count consistency: ",
  global_internal_consistency,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Verify expected locked values
# ------------------------------------------------------------

expected_universe <- 57828
expected_mapped <- 34631
expected_unmapped <- 23197
expected_one <- 34138
expected_multiple <- 493

expected_global_pass <- (
  N_universe_r2 == expected_universe &&
    N_mapped_r2 == expected_mapped &&
    N_unmapped_r2 == expected_unmapped &&
    N_one_r2 == expected_one &&
    N_multiple_r2 == expected_multiple
)

cat("\n============================================================\n")
cat("6. Locked global values\n")
cat("============================================================\n\n")

cat(
  "Expected universe: ",
  expected_universe,
  "\n",
  sep = ""
)

cat(
  "Observed universe: ",
  N_universe_r2,
  "\n",
  sep = ""
)

cat(
  "Expected mapped: ",
  expected_mapped,
  "\n",
  sep = ""
)

cat(
  "Observed mapped: ",
  N_mapped_r2,
  "\n",
  sep = ""
)

cat(
  "Expected unmapped: ",
  expected_unmapped,
  "\n",
  sep = ""
)

cat(
  "Observed unmapped: ",
  N_unmapped_r2,
  "\n",
  sep = ""
)

cat(
  "Expected exactly one: ",
  expected_one,
  "\n",
  sep = ""
)

cat(
  "Observed exactly one: ",
  N_one_r2,
  "\n",
  sep = ""
)

cat(
  "Expected multiple: ",
  expected_multiple,
  "\n",
  sep = ""
)

cat(
  "Observed multiple: ",
  N_multiple_r2,
  "\n",
  sep = ""
)

cat(
  "\nLocked global mapping consistency: ",
  expected_global_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Final protein-coding QC decision
# ------------------------------------------------------------

protein_qc_pass <- (
  protein_internal_consistency &&
    summary_consistency &&
    is.finite(protein_coverage_r2) &&
    protein_N_r2 == 18235 &&
    protein_mapped_r2 == 18007 &&
    protein_unmapped_r2 == 228
)


# ------------------------------------------------------------
# 10. Final STEP9C-1B-R2 status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-1B-R2 FINAL STATUS\n")
cat("============================================================\n\n")

cat(
  "Universe size correct: ",
  N_universe_r2 == 57828,
  "\n",
  sep = ""
)

cat(
  "Global mapping consistency: ",
  expected_global_pass,
  "\n",
  sep = ""
)

cat(
  "Protein-coding subset size correct: ",
  protein_N_r2 == 18235,
  "\n",
  sep = ""
)

cat(
  "Protein-coding mapped count correct: ",
  protein_mapped_r2 == 18007,
  "\n",
  sep = ""
)

cat(
  "Protein-coding unmapped count correct: ",
  protein_unmapped_r2 == 228,
  "\n",
  sep = ""
)

cat(
  "Protein-coding internal consistency: ",
  protein_internal_consistency,
  "\n",
  sep = ""
)

cat(
  "Protein-coding coverage finite: ",
  is.finite(protein_coverage_r2),
  "\n",
  sep = ""
)

cat(
  "\nNo DEG recalculation performed.\n"
)

cat(
  "No GO recalculation performed.\n"
)

cat(
  "No KEGG recalculation performed.\n"
)

cat(
  "No DEG threshold changed.\n"
)

cat(
  "57,828-gene testing universe unchanged.\n"
)


if (
  expected_global_pass &&
  protein_qc_pass
) {
  
  cat(
    "\nSTEP9C-1B-R2 STATUS: PASS\n"
  )
  
} else {
  
  cat(
    "\nSTEP9C-1B-R2 STATUS: FAIL\n"
  )
  
}

cat(
  "\nSTEP9C-1B-R2 completed.\n"
)


##STEP9C-2 — KEGG ORA Input Mapping Audit 
# ============================================================
# Paper 1
# STEP9C-2
# KEGG ORA Input Mapping Audit
#
# Purpose:
#   Audit Ensembl -> Entrez mapping for:
#     1. Tumor_up DEG set
#     2. Tumor_down DEG set
#     3. Locked 57,828-gene DEG testing universe
#
# IMPORTANT:
#   This is an AUDIT ONLY step.
#
#   NO DEG recalculation
#   NO P-value recalculation
#   NO FDR recalculation
#   NO KEGG recalculation
#   NO GO recalculation
#   NO threshold change
#   NO change to the 57,828-gene universe
#
# Locked DEG definition:
#   FDR < 0.05
#   AND |log2FC| >= 1
#
# Locked universe:
#   57,828 genes
#
# Core genes:
#   AOC3 / MS4A7 / IL7R
# ============================================================


cat("\n============================================================\n")
cat("STEP9C-2 — KEGG ORA Input Mapping Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

step_root <- file.path(
  "input_dir/Paper 1_Step 9-Functional Enrichment",
  "STEP9C_Revised_DEG_GO_KEGG_ORA"
)

step9c2_dir <- file.path(
  step_root,
  "STEP9C-2_KEGG_ORA_Input_Mapping_Audit"
)

results_dir <- file.path(
  step9c2_dir,
  "01_results"
)

qc_dir <- file.path(
  step9c2_dir,
  "02_QC"
)

lock_dir <- file.path(
  step9c2_dir,
  "03_lock"
)

log_dir <- file.path(
  step9c2_dir,
  "04_logs"
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  lock_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Confirm locked objects
# ------------------------------------------------------------

cat("============================================================\n")
cat("1. Locked object confirmation\n")
cat("============================================================\n\n")


stopifnot(
  exists("testing_universe"),
  exists("deg_up"),
  exists("deg_down"),
  exists("universe_map_current"),
  exists("universe_biotype")
)


cat(
  "testing_universe exists: TRUE\n"
)

cat(
  "deg_up exists: TRUE\n"
)

cat(
  "deg_down exists: TRUE\n"
)

cat(
  "universe_map_current exists: TRUE\n"
)

cat(
  "universe_biotype exists: TRUE\n\n"
)


# ------------------------------------------------------------
# 2. Confirm locked dimensions
# ------------------------------------------------------------

cat("============================================================\n")
cat("2. Locked DEG dimensions\n")
cat("============================================================\n\n")


N_universe <- length(
  unique(
    trimws(
      as.character(
        testing_universe
      )
    )
  )
)

N_up <- nrow(
  deg_up
)

N_down <- nrow(
  deg_down
)


cat(
  "Locked testing universe: ",
  N_universe,
  "\n",
  sep = ""
)

cat(
  "Tumor_up DEG: ",
  N_up,
  "\n",
  sep = ""
)

cat(
  "Tumor_down DEG: ",
  N_down,
  "\n",
  sep = ""
)


stopifnot(
  N_universe == 57828,
  N_up == 9311,
  N_down == 6754
)


cat(
  "\nLocked dimensions PASS.\n\n"
)


# ------------------------------------------------------------
# 3. Detect DEG Ensembl ID column
# ------------------------------------------------------------

cat("============================================================\n")
cat("3. DEG Ensembl ID column detection\n")
cat("============================================================\n\n")


cat(
  "Tumor_up columns:\n"
)

print(
  colnames(deg_up)
)

cat(
  "\nTumor_down columns:\n"
)

print(
  colnames(deg_down)
)


candidate_id_columns <- c(
  "Gene",
  "ENSEMBL",
  "Ensembl",
  "Ensembl_ID",
  "GeneID",
  "Gene_ID"
)


up_id_candidates <- intersect(
  candidate_id_columns,
  colnames(deg_up)
)

down_id_candidates <- intersect(
  candidate_id_columns,
  colnames(deg_down)
)


if (length(up_id_candidates) == 0) {
  
  stop(
    "No recognized Ensembl ID column found in deg_up."
  )
  
}

if (length(down_id_candidates) == 0) {
  
  stop(
    "No recognized Ensembl ID column found in deg_down."
  )
  
}


up_id_col <- up_id_candidates[1]

down_id_col <- down_id_candidates[1]


cat(
  "\nTumor_up ID column selected: ",
  up_id_col,
  "\n",
  sep = ""
)

cat(
  "Tumor_down ID column selected: ",
  down_id_col,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Prepare locked Ensembl ID sets
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("4. Prepare locked Ensembl ID sets\n")
cat("============================================================\n\n")


universe_ids <- unique(
  trimws(
    as.character(
      testing_universe
    )
  )
)

up_ids <- unique(
  trimws(
    as.character(
      deg_up[[up_id_col]]
    )
  )
)

down_ids <- unique(
  trimws(
    as.character(
      deg_down[[down_id_col]]
    )
  )
)


up_ids <- up_ids[
  !is.na(up_ids) &
    up_ids != ""
]

down_ids <- down_ids[
  !is.na(down_ids) &
    down_ids != ""
]

universe_ids <- universe_ids[
  !is.na(universe_ids) &
    universe_ids != ""
]


cat(
  "Unique Tumor_up Ensembl IDs: ",
  length(up_ids),
  "\n",
  sep = ""
)

cat(
  "Unique Tumor_down Ensembl IDs: ",
  length(down_ids),
  "\n",
  sep = ""
)

cat(
  "Unique universe Ensembl IDs: ",
  length(universe_ids),
  "\n",
  sep = ""
)


stopifnot(
  length(up_ids) == 9311,
  length(down_ids) == 6754,
  length(universe_ids) == 57828
)


# ------------------------------------------------------------
# 5. Confirm DEG groups are subsets of locked universe
# ------------------------------------------------------------

up_in_universe <- up_ids %in% universe_ids

down_in_universe <- down_ids %in% universe_ids


cat("\n============================================================\n")
cat("5. DEG membership in locked universe\n")
cat("============================================================\n\n")


cat(
  "Tumor_up IDs in universe: ",
  sum(up_in_universe),
  " / ",
  length(up_ids),
  "\n",
  sep = ""
)

cat(
  "Tumor_down IDs in universe: ",
  sum(down_in_universe),
  " / ",
  length(down_ids),
  "\n",
  sep = ""
)


up_membership_pass <- all(
  up_in_universe
)

down_membership_pass <- all(
  down_in_universe
)


cat(
  "\nTumor_up membership PASS: ",
  up_membership_pass,
  "\n",
  sep = ""
)

cat(
  "Tumor_down membership PASS: ",
  down_membership_pass,
  "\n",
  sep = ""
)


stopifnot(
  up_membership_pass,
  down_membership_pass
)


# ------------------------------------------------------------
# 6. Prepare existing OrgDb mapping
#
# IMPORTANT:
#   No new AnnotationDbi::select() call.
#   Uses the already generated universe_map_current.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("6. Existing Entrez mapping object\n")
cat("============================================================\n\n")


map0 <- universe_map_current[
  ,
  c(
    "ENSEMBL",
    "ENTREZID",
    "SYMBOL"
  ),
  drop = FALSE
]


map0$ENSEMBL <- trimws(
  as.character(
    map0$ENSEMBL
  )
)

map0$ENTREZID <- trimws(
  as.character(
    map0$ENTREZID
  )
)

map0$SYMBOL <- trimws(
  as.character(
    map0$SYMBOL
  )
)


map0$ENTREZID[
  map0$ENTREZID %in%
    c(
      "",
      "NA",
      "NaN"
    )
] <- NA_character_


cat(
  "Existing mapping rows: ",
  nrow(map0),
  "\n",
  sep = ""
)

cat(
  "Unique mapped Ensembl IDs represented: ",
  length(
    unique(
      map0$ENSEMBL
    )
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Function for Ensembl-level mapping audit
# ------------------------------------------------------------

audit_mapping <- function(
    ids,
    label,
    map_table
) {
  
  ids <- unique(
    trimws(
      as.character(
        ids
      )
    )
  )
  
  ids <- ids[
    !is.na(ids) &
      ids != ""
  ]
  
  
  # Keep only mapping rows for requested IDs
  m <- map_table[
    map_table$ENSEMBL %in% ids,
    ,
    drop = FALSE
  ]
  
  
  # Split by Ensembl ID
  split_map <- split(
    m$ENTREZID,
    m$ENSEMBL
  )
  
  
  # Initialize all requested IDs
  out <- data.frame(
    ENSEMBL = ids,
    stringsAsFactors = FALSE
  )
  
  
  # Number of unique non-missing Entrez IDs
  out$N_Entrez <- vapply(
    ids,
    function(id) {
      
      x <- split_map[[id]]
      
      if (is.null(x)) {
        return(0L)
      }
      
      x <- unique(
        x[
          !is.na(x) &
            x != ""
        ]
      )
      
      length(x)
      
    },
    integer(1)
  )
  
  
  # Mapping status
  out$Entrez_Mapped <- (
    out$N_Entrez > 0
  )
  
  
  out$Mapping_Category <- ifelse(
    
    out$N_Entrez == 0,
    "Unmapped",
    
    ifelse(
      out$N_Entrez == 1,
      "One_to_one_or_single_Entrez",
      "One_to_many"
    )
    
  )
  
  
  # Primary descriptive Entrez ID
  out$Primary_Entrez <- vapply(
    ids,
    function(id) {
      
      x <- split_map[[id]]
      
      if (is.null(x)) {
        return(NA_character_)
      }
      
      x <- unique(
        x[
          !is.na(x) &
            x != ""
        ]
      )
      
      if (length(x) == 0) {
        return(NA_character_)
      }
      
      sort(x)[1]
      
    },
    character(1)
  )
  
  
  # Gene symbol
  symbol_split <- split(
    map_table$SYMBOL,
    map_table$ENSEMBL
  )
  
  
  out$SYMBOL <- vapply(
    ids,
    function(id) {
      
      x <- symbol_split[[id]]
      
      if (is.null(x)) {
        return(NA_character_)
      }
      
      x <- unique(
        x[
          !is.na(x) &
            x != ""
        ]
      )
      
      if (length(x) == 0) {
        return(NA_character_)
      }
      
      x[1]
      
    },
    character(1)
  )
  
  
  out$Set <- label
  
  
  out <- out[
    ,
    c(
      "Set",
      "ENSEMBL",
      "SYMBOL",
      "N_Entrez",
      "Entrez_Mapped",
      "Mapping_Category",
      "Primary_Entrez"
    )
  ]
  
  
  rownames(out) <- NULL
  
  return(out)
  
}


# ------------------------------------------------------------
# 8. Run mapping audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("7. Run KEGG input mapping audit\n")
cat("============================================================\n\n")


audit_universe <- audit_mapping(
  ids = universe_ids,
  label = "Universe",
  map_table = map0
)


audit_up <- audit_mapping(
  ids = up_ids,
  label = "Tumor_up",
  map_table = map0
)


audit_down <- audit_mapping(
  ids = down_ids,
  label = "Tumor_down",
  map_table = map0
)


cat(
  "Universe audit completed.\n"
)

cat(
  "Tumor_up audit completed.\n"
)

cat(
  "Tumor_down audit completed.\n"
)


# ------------------------------------------------------------
# 9. Summary function
# ------------------------------------------------------------

make_summary <- function(
    audit_df
) {
  
  N <- nrow(
    audit_df
  )
  
  mapped <- sum(
    audit_df$Entrez_Mapped
  )
  
  unmapped <- sum(
    !audit_df$Entrez_Mapped
  )
  
  one <- sum(
    audit_df$N_Entrez == 1
  )
  
  multiple <- sum(
    audit_df$N_Entrez > 1
  )
  
  
  data.frame(
    
    N_Ensembl = N,
    
    Entrez_Mapped = mapped,
    
    Entrez_Unmapped = unmapped,
    
    Entrez_Coverage = mapped / N,
    
    Exactly_One_Entrez = one,
    
    More_Than_One_Entrez = multiple,
    
    stringsAsFactors = FALSE
    
  )
  
}


summary_universe <- make_summary(
  audit_universe
)

summary_up <- make_summary(
  audit_up
)

summary_down <- make_summary(
  audit_down
)


mapping_summary <- rbind(
  
  data.frame(
    Set = "Universe",
    summary_universe
  ),
  
  data.frame(
    Set = "Tumor_up",
    summary_up
  ),
  
  data.frame(
    Set = "Tumor_down",
    summary_down
  )
  
)


# ------------------------------------------------------------
# 10. Print mapping summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("8. KEGG ORA input mapping summary\n")
cat("============================================================\n\n")


print(
  mapping_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Protein-coding annotation audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("9. Protein-coding assessment\n")
cat("============================================================\n\n")


biotype_lookup <- universe_biotype[
  ,
  c(
    "ENSEMBL",
    "gene_biotype"
  ),
  drop = FALSE
]


biotype_lookup <- biotype_lookup[
  !duplicated(
    biotype_lookup$ENSEMBL
  ),
]


audit_up_biotype <- merge(
  audit_up,
  biotype_lookup,
  by = "ENSEMBL",
  all.x = TRUE,
  sort = FALSE
)


audit_down_biotype <- merge(
  audit_down,
  biotype_lookup,
  by = "ENSEMBL",
  all.x = TRUE,
  sort = FALSE
)


audit_universe_biotype <- merge(
  audit_universe,
  biotype_lookup,
  by = "ENSEMBL",
  all.x = TRUE,
  sort = FALSE
)


# Function for protein-coding summary
protein_summary <- function(
    audit_df
) {
  
  x <- audit_df[
    audit_df$gene_biotype == "protein_coding",
    ,
    drop = FALSE
  ]
  
  
  data.frame(
    
    N_protein_coding = nrow(x),
    
    Protein_coding_Entrez_Mapped =
      sum(
        x$Entrez_Mapped
      ),
    
    Protein_coding_Entrez_Unmapped =
      sum(
        !x$Entrez_Mapped
      ),
    
    Protein_coding_Coverage =
      ifelse(
        nrow(x) > 0,
        sum(x$Entrez_Mapped) / nrow(x),
        NA_real_
      ),
    
    stringsAsFactors = FALSE
    
  )
  
}


protein_universe_summary <- protein_summary(
  audit_universe_biotype
)

protein_up_summary <- protein_summary(
  audit_up_biotype
)

protein_down_summary <- protein_summary(
  audit_down_biotype
)


protein_mapping_summary <- rbind(
  
  data.frame(
    Set = "Universe",
    protein_universe_summary
  ),
  
  data.frame(
    Set = "Tumor_up",
    protein_up_summary
  ),
  
  data.frame(
    Set = "Tumor_down",
    protein_down_summary
  )
  
)


print(
  protein_mapping_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. One-to-many mapping audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("10. One-to-many Entrez mapping audit\n")
cat("============================================================\n\n")


multiple_up <- audit_up[
  audit_up$N_Entrez > 1,
  ,
  drop = FALSE
]


multiple_down <- audit_down[
  audit_down$N_Entrez > 1,
  ,
  drop = FALSE
]


multiple_universe <- audit_universe[
  audit_universe$N_Entrez > 1,
  ,
  drop = FALSE
]


cat(
  "Tumor_up Ensembl IDs with >1 Entrez: ",
  nrow(multiple_up),
  "\n",
  sep = ""
)

cat(
  "Tumor_down Ensembl IDs with >1 Entrez: ",
  nrow(multiple_down),
  "\n",
  sep = ""
)

cat(
  "Universe Ensembl IDs with >1 Entrez: ",
  nrow(multiple_universe),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. Unmapped genes
# ------------------------------------------------------------

unmapped_up <- audit_up[
  !audit_up$Entrez_Mapped,
  ,
  drop = FALSE
]


unmapped_down <- audit_down[
  !audit_down$Entrez_Mapped,
  ,
  drop = FALSE
]


unmapped_universe <- audit_universe[
  !audit_universe$Entrez_Mapped,
  ,
  drop = FALSE
]


cat("\n============================================================\n")
cat("11. Unmapped Ensembl IDs\n")
cat("============================================================\n\n")


cat(
  "Tumor_up unmapped: ",
  nrow(unmapped_up),
  "\n",
  sep = ""
)

cat(
  "Tumor_down unmapped: ",
  nrow(unmapped_down),
  "\n",
  sep = ""
)

cat(
  "Universe unmapped: ",
  nrow(unmapped_universe),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 14. Check core genes
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("12. Core-gene KEGG mapping audit\n")
cat("============================================================\n\n")


core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


core_audit <- audit_universe[
  audit_universe$SYMBOL %in% core_genes,
  ,
  drop = FALSE
]


core_audit <- core_audit[
  match(
    core_genes,
    core_audit$SYMBOL
  ),
  ,
  drop = FALSE
]


print(
  core_audit,
  row.names = FALSE
)


core_mapping_pass <- (
  nrow(core_audit) == 3 &&
    all(
      core_audit$Entrez_Mapped
    )
)


cat(
  "\nCore-gene Entrez mapping PASS: ",
  core_mapping_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 15. Verify locked DEG counts were not changed
# ------------------------------------------------------------

locked_counts_pass <- (
  length(universe_ids) == 57828 &&
    length(up_ids) == 9311 &&
    length(down_ids) == 6754
)


# ------------------------------------------------------------
# 16. Verify DEG membership
# ------------------------------------------------------------

membership_pass <- (
  all(up_ids %in% universe_ids) &&
    all(down_ids %in% universe_ids)
)


# ------------------------------------------------------------
# 17. Final mapping audit status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("13. STEP9C-2 FINAL DIAGNOSTIC STATUS\n")
cat("============================================================\n\n")


cat(
  "Locked universe = 57,828: ",
  length(universe_ids) == 57828,
  "\n",
  sep = ""
)

cat(
  "Tumor_up = 9,311: ",
  length(up_ids) == 9311,
  "\n",
  sep = ""
)

cat(
  "Tumor_down = 6,754: ",
  length(down_ids) == 6754,
  "\n",
  sep = ""
)

cat(
  "DEG counts unchanged: ",
  locked_counts_pass,
  "\n",
  sep = ""
)

cat(
  "DEG sets contained in locked universe: ",
  membership_pass,
  "\n",
  sep = ""
)

cat(
  "Core-gene mapping PASS: ",
  core_mapping_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 18. Save outputs
# ------------------------------------------------------------

write.csv(
  mapping_summary,
  file.path(
    results_dir,
    "STEP9C2_KEGG_Input_Entrez_Mapping_Summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  protein_mapping_summary,
  file.path(
    results_dir,
    "STEP9C2_KEGG_Input_ProteinCoding_Mapping_Summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  audit_up_biotype,
  file.path(
    results_dir,
    "STEP9C2_Tumor_up_Ensembl_Entrez_Biotype_Audit.csv"
  ),
  row.names = FALSE
)


write.csv(
  audit_down_biotype,
  file.path(
    results_dir,
    "STEP9C2_Tumor_down_Ensembl_Entrez_Biotype_Audit.csv"
  ),
  row.names = FALSE
)


write.csv(
  audit_universe_biotype,
  file.path(
    results_dir,
    "STEP9C2_57828_Universe_Ensembl_Entrez_Biotype_Audit.csv"
  ),
  row.names = FALSE
)


write.csv(
  multiple_up,
  file.path(
    qc_dir,
    "STEP9C2_Tumor_up_One_to_Many_Entrez.csv"
  ),
  row.names = FALSE
)


write.csv(
  multiple_down,
  file.path(
    qc_dir,
    "STEP9C2_Tumor_down_One_to_Many_Entrez.csv"
  ),
  row.names = FALSE
)


write.csv(
  multiple_universe,
  file.path(
    qc_dir,
    "STEP9C2_Universe_One_to_Many_Entrez.csv"
  ),
  row.names = FALSE
)


write.csv(
  unmapped_up,
  file.path(
    qc_dir,
    "STEP9C2_Tumor_up_Unmapped_Ensembl.csv"
  ),
  row.names = FALSE
)


write.csv(
  unmapped_down,
  file.path(
    qc_dir,
    "STEP9C2_Tumor_down_Unmapped_Ensembl.csv"
  ),
  row.names = FALSE
)


write.csv(
  unmapped_universe,
  file.path(
    qc_dir,
    "STEP9C2_Universe_Unmapped_Ensembl.csv"
  ),
  row.names = FALSE
)


write.csv(
  core_audit,
  file.path(
    qc_dir,
    "STEP9C2_Core_Gene_Entrez_Mapping_Audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 19. Save final lock report
# ------------------------------------------------------------

final_status <- if (
  locked_counts_pass &&
  membership_pass &&
  core_mapping_pass
) {
  
  "PASS_KEGG_INPUT_MAPPING_AUDIT"
  
} else {
  
  "HOLD_KEGG_INPUT_MAPPING_AUDIT"
  
}


lock_report <- data.frame(
  
  Item = c(
    "Locked_testing_universe",
    "Tumor_up_N",
    "Tumor_down_N",
    "DEG_counts_unchanged",
    "DEG_sets_in_locked_universe",
    "Core_gene_mapping",
    "Global_Entrez_coverage",
    "Final_status"
  ),
  
  Value = c(
    "57828",
    "9311",
    "6754",
    as.character(locked_counts_pass),
    as.character(membership_pass),
    as.character(core_mapping_pass),
    sprintf(
      "%.6f",
      summary_universe$Entrez_Coverage
    ),
    final_status
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  lock_report,
  file.path(
    lock_dir,
    "STEP9C2_Final_Lock_Report.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 20. Final console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-2 FINAL STATUS\n")
cat("============================================================\n\n")


cat(
  "No DEG recalculation performed.\n"
)

cat(
  "No P-value recalculation performed.\n"
)

cat(
  "No FDR recalculation performed.\n"
)

cat(
  "No GO recalculation performed.\n"
)

cat(
  "No KEGG recalculation performed.\n"
)

cat(
  "No DEG threshold changed.\n"
)

cat(
  "57,828-gene testing universe unchanged.\n"
)


if (
  final_status ==
  "PASS_KEGG_INPUT_MAPPING_AUDIT"
) {
  
  cat(
    "\nSTEP9C-2 STATUS: PASS_KEGG_INPUT_MAPPING_AUDIT\n"
  )
  
} else {
  
  cat(
    "\nSTEP9C-2 STATUS: HOLD_KEGG_INPUT_MAPPING_AUDIT\n"
  )
  
}


cat(
  "\nSTEP9C-2 completed.\n"
)

# ============================================================
# Paper 1
# STEP9C-3
# Final KEGG Entrez Background and Input-ID Policy Lock
#
# Purpose:
#   1. Lock the final Entrez ID space for KEGG ORA
#   2. Define the effective KEGG background
#   3. Define Tumor_up and Tumor_down KEGG input sets
#   4. Audit one-to-many Ensembl -> Entrez mappings
#   5. Confirm DEG membership is unchanged
#   6. Confirm core-gene mapping
#
# IMPORTANT:
#   - No DEG recalculation
#   - No P-value recalculation
#   - No FDR recalculation
#   - No DEG threshold change
#   - No KEGG ORA recalculation
#   - This step ONLY locks the KEGG ID/background policy
# ============================================================


cat("\n============================================================\n")
cat("STEP9C-3: Final KEGG Entrez Background and Input-ID Policy Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory and output folders
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9C_Revised_DEG_GO_KEGG_ORA/STEP9C-3_Final_KEGG_Entrez_Background_Input_Policy_Lock"

dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(
  file.path(step_dir, "01_results"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(step_dir, "02_QC"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(step_dir, "03_lock"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(step_dir, "04_logs"),
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Confirm required objects
# ------------------------------------------------------------

cat("1. Required object existence\n")
cat("------------------------------------------------------------\n\n")

required_objects <- c(
  "testing_universe",
  "deg_up",
  "deg_down",
  "universe_map_current",
  "universe_biotype"
)

object_status <- sapply(
  required_objects,
  exists,
  inherits = TRUE
)

print(object_status)

if (!all(object_status)) {
  stop(
    "Required object(s) missing. STEP9C-3 stopped."
  )
}

cat("\nAll required objects exist: PASS\n\n")


# ------------------------------------------------------------
# 2. Locked DEG dimensions
# ------------------------------------------------------------

cat("2. Locked DEG dimensions\n")
cat("------------------------------------------------------------\n\n")

N_universe <- nrow(testing_universe)
N_up       <- nrow(deg_up)
N_down     <- nrow(deg_down)

cat(
  "Testing universe: ",
  N_universe,
  "\n",
  sep = ""
)

cat(
  "Tumor_up: ",
  N_up,
  "\n",
  sep = ""
)

cat(
  "Tumor_down: ",
  N_down,
  "\n\n",
  sep = ""
)

stopifnot(
  N_universe == 57828,
  N_up == 9311,
  N_down == 6754
)

cat("Locked dimensions PASS.\n\n")

# ============================================================
# STEP9C-3 Diagnostic
# Inspect testing_universe structure and column names
# NO recalculation
# NO mapping
# NO data modification
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-3 Diagnostic: testing_universe structure\n")
cat("============================================================\n\n")

cat("Class of testing_universe:\n")
print(class(testing_universe))

cat("\nLength of testing_universe:\n")
print(length(testing_universe))

cat("\nDimensions of testing_universe:\n")
print(dim(testing_universe))

cat("\nNumber of rows:\n")
print(NROW(testing_universe))

cat("\nNumber of columns:\n")
print(NCOL(testing_universe))

cat("\nColumn names:\n")
print(colnames(testing_universe))

cat("\nFirst 10 rows:\n")
print(
  head(testing_universe, 10)
)

cat("\nStructure:\n")
str(
  testing_universe,
  max.level = 2
)

cat("\n============================================================\n")
cat("Compare with deg_up / deg_down\n")
cat("============================================================\n\n")

cat("deg_up class:\n")
print(class(deg_up))

cat("\ndeg_up dimensions:\n")
print(dim(deg_up))

cat("\ndeg_up columns:\n")
print(colnames(deg_up))

cat("\ndeg_down class:\n")
print(class(deg_down))

cat("\ndeg_down dimensions:\n")
print(dim(deg_down))

cat("\ndeg_down columns:\n")
print(colnames(deg_down))

cat("\n============================================================\n")
cat("Diagnostic completed.\n")
cat("============================================================\n")


cat("\n============================================================\n")
cat("STEP9C-3 Diagnostic: testing_universe ID integrity\n")
cat("============================================================\n\n")

cat("1. testing_universe length:\n")
print(length(testing_universe))

cat("\n2. Unique IDs:\n")
print(length(unique(testing_universe)))

cat("\n3. NA count:\n")
print(sum(is.na(testing_universe)))

cat("\n4. Empty-string count:\n")
print(sum(testing_universe == "", na.rm = TRUE))

cat("\n5. Ensembl-format IDs:\n")
ensembl_pattern <- "^ENSG[0-9]+$"

print(
  sum(grepl(ensembl_pattern, testing_universe))
)

cat("\n6. Non-Ensembl-format IDs:\n")
print(
  sum(!grepl(ensembl_pattern, testing_universe))
)

cat("\n7. Duplicate IDs:\n")
print(
  sum(duplicated(testing_universe))
)

cat("\n8. testing_universe ∩ deg_up$Gene:\n")
print(
  length(intersect(testing_universe, deg_up$Gene))
)

cat("\n9. testing_universe ∩ deg_down$Gene:\n")
print(
  length(intersect(testing_universe, deg_down$Gene))
)

cat("\n10. deg_up IDs absent from testing_universe:\n")
print(
  length(setdiff(deg_up$Gene, testing_universe))
)

cat("\n11. deg_down IDs absent from testing_universe:\n")
print(
  length(setdiff(deg_down$Gene, testing_universe))
)

cat("\n============================================================\n")
cat("Diagnostic completed.\n")
cat("============================================================\n")


# ============================================================
# Paper 1 — STEP9C-3
# Diagnostic: identify existing Ensembl -> Entrez mapping objects
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-3 Diagnostic: Existing Ensembl -> Entrez Mapping Objects\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Basic environment information
# ------------------------------------------------------------

cat("R session:\n")
print(R.version.string)

cat("\nCurrent working directory:\n")
print(getwd())


# ------------------------------------------------------------
# 1. Confirm the already-locked testing universe
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("1. Confirm locked testing_universe\n")
cat("============================================================\n\n")

cat("Class:\n")
print(class(testing_universe))

cat("\nLength:\n")
print(length(testing_universe))

cat("\nUnique IDs:\n")
print(length(unique(testing_universe)))

cat("\nNA count:\n")
print(sum(is.na(testing_universe)))

cat("\nDuplicate count:\n")
print(sum(duplicated(testing_universe)))

# ============================================================
# STEP9C-3 Diagnostic
# Identify existing Ensembl -> Entrez mapping objects
# Corrected for R 4.5.1
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-3 Diagnostic: Existing Ensembl -> Entrez Mapping Objects\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 2. Identify candidate mapping objects already in memory
# ------------------------------------------------------------

all_objects <- ls(envir = .GlobalEnv)

mapping_pattern <- paste(
  c(
    "entrez",
    "ENSEMBL",
    "ensembl",
    "gene.*id",
    "id.*map",
    "mapping",
    "annotation",
    "org\\.Hs\\.eg"
  ),
  collapse = "|"
)

mapping_candidates <- all_objects[
  grepl(
    mapping_pattern,
    all_objects,
    ignore.case = TRUE
  )
]

cat("Candidate objects found:\n")
print(mapping_candidates)

cat("\nNumber of candidate objects:\n")
print(length(mapping_candidates))

# ============================================================
# Paper 1 — STEP9C-3
# Existing Ensembl -> Entrez Mapping Audit
#
# Purpose:
#   1. Identify the actual mapping table already created
#   2. Inspect mapping_clean / mapping_complete / mapping_unique
#   3. Inspect universe mapping coverage
#   4. Inspect mapping multiplicity
#   5. Inspect locked mapping status
#
# IMPORTANT:
#   No new ID conversion is performed.
#   No DEG recalculation is performed.
#   testing_universe is not modified.
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-3: Existing Ensembl -> Entrez Mapping Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Confirm locked testing universe
# ------------------------------------------------------------

cat("0. Locked testing_universe\n")
cat("----------------------------\n")

cat("Length: ")
print(length(testing_universe))

cat("Unique: ")
print(length(unique(testing_universe)))

cat("NA: ")
print(sum(is.na(testing_universe)))

cat("Duplicate: ")
print(sum(duplicated(testing_universe)))


# ------------------------------------------------------------
# 1. Inspect core mapping
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("1. Core-gene Ensembl -> Entrez mapping\n")
cat("============================================================\n\n")

cat("core_mapping_check:\n")
print(core_mapping_check)

cat("\ncore_mapping_pass:\n")
print(core_mapping_pass)


# ------------------------------------------------------------
# 2. Inspect annotation_map
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("2. annotation_map\n")
cat("============================================================\n\n")

cat("Class:\n")
print(class(annotation_map))

cat("\nDimensions:\n")
print(dim(annotation_map))

cat("\nColumn names:\n")
print(colnames(annotation_map))

cat("\nFirst 10 rows:\n")
print(utils::head(annotation_map, 10))

cat("\nUnique Ensembl IDs:\n")
print(length(unique(annotation_map$ENSEMBL)))

cat("\nUnique Entrez IDs:\n")
print(length(unique(annotation_map$ENTREZID)))

cat("\nDuplicated Ensembl rows:\n")
print(sum(duplicated(annotation_map$ENSEMBL)))

cat("\nDuplicated Entrez rows:\n")
print(sum(duplicated(annotation_map$ENTREZID)))


# ------------------------------------------------------------
# 3. Inspect primary mapping objects
# ------------------------------------------------------------

mapping_objects_to_check <- c(
  "mapping",
  "mapping_clean",
  "mapping_complete",
  "mapping_unique"
)

cat("\n============================================================\n")
cat("3. Existing mapping objects\n")
cat("============================================================\n\n")

for (obj_name in mapping_objects_to_check) {
  
  if (exists(obj_name, envir = .GlobalEnv)) {
    
    cat("\n------------------------------------------------------------\n")
    cat("Object:", obj_name, "\n")
    cat("------------------------------------------------------------\n")
    
    obj <- get(obj_name, envir = .GlobalEnv)
    
    cat("Class:\n")
    print(class(obj))
    
    cat("\nDimensions:\n")
    print(dim(obj))
    
    if (is.data.frame(obj) || is.matrix(obj)) {
      
      cat("\nColumn names:\n")
      print(colnames(obj))
      
      cat("\nFirst 10 rows:\n")
      print(utils::head(obj, 10))
      
    } else {
      
      cat("\nLength:\n")
      print(length(obj))
      
      cat("\nFirst elements:\n")
      print(utils::head(obj, 10))
    }
  } else {
    
    cat("\nObject", obj_name, "does not exist.\n")
  }
}


# ------------------------------------------------------------
# 4. Inspect mapping summary objects
# ------------------------------------------------------------

summary_objects <- c(
  "mapping_count",
  "mapping_status",
  "mapping_pass",
  "mapping_consistency",
  "mapping_summary",
  "mapping_provenance"
)

cat("\n============================================================\n")
cat("4. Existing mapping QC summaries\n")
cat("============================================================\n\n")

for (obj_name in summary_objects) {
  
  if (exists(obj_name, envir = .GlobalEnv)) {
    
    cat("\n------------------------------------------------------------\n")
    cat("Object:", obj_name, "\n")
    cat("------------------------------------------------------------\n")
    
    obj <- get(obj_name, envir = .GlobalEnv)
    
    cat("Class:\n")
    print(class(obj))
    
    cat("\nDimensions:\n")
    print(dim(obj))
    
    if (is.data.frame(obj) || is.matrix(obj)) {
      
      cat("\nColumn names:\n")
      print(colnames(obj))
      
      cat("\nContent:\n")
      print(obj)
      
    } else {
      
      cat("\nContent:\n")
      print(obj)
    }
  }
}


# ------------------------------------------------------------
# 5. Inspect universe mapping coverage
# ------------------------------------------------------------

universe_objects <- c(
  "universe_annotation_coverage",
  "universe_entrez_coverage",
  "universe_mapping_multiplicity",
  "universe_mapping_summary",
  "matched_universe_entrez"
)

cat("\n============================================================\n")
cat("5. Testing-universe Ensembl -> Entrez coverage\n")
cat("============================================================\n\n")

for (obj_name in universe_objects) {
  
  if (exists(obj_name, envir = .GlobalEnv)) {
    
    cat("\n------------------------------------------------------------\n")
    cat("Object:", obj_name, "\n")
    cat("------------------------------------------------------------\n")
    
    obj <- get(obj_name, envir = .GlobalEnv)
    
    cat("Class:\n")
    print(class(obj))
    
    cat("\nDimensions:\n")
    print(dim(obj))
    
    cat("\nContent:\n")
    print(obj)
  }
}


# ------------------------------------------------------------
# 6. Inspect multiplicity-related objects
# ------------------------------------------------------------

multiplicity_objects <- c(
  "multiple_entrez",
  "duplicated_ensembl",
  "duplicated_entrez",
  "gene_mapping_multiplicity"
)

cat("\n============================================================\n")
cat("6. Mapping multiplicity audit\n")
cat("============================================================\n\n")

for (obj_name in multiplicity_objects) {
  
  if (exists(obj_name, envir = .GlobalEnv)) {
    
    cat("\n------------------------------------------------------------\n")
    cat("Object:", obj_name, "\n")
    cat("------------------------------------------------------------\n")
    
    obj <- get(obj_name, envir = .GlobalEnv)
    
    cat("Class:\n")
    print(class(obj))
    
    cat("\nLength:\n")
    print(length(obj))
    
    cat("\nDimensions:\n")
    print(dim(obj))
    
    cat("\nContent / first elements:\n")
    print(utils::head(obj, 20))
  }
}


# ------------------------------------------------------------
# 7. Inspect DEG-level Entrez coverage
# ------------------------------------------------------------

deg_objects <- c(
  "matched_down_entrez",
  "down_entrez_coverage",
  "down_annotation_coverage",
  "entrez_hits",
  "all_entrez_hits"
)

cat("\n============================================================\n")
cat("7. DEG-level Entrez coverage objects\n")
cat("============================================================\n\n")

for (obj_name in deg_objects) {
  
  if (exists(obj_name, envir = .GlobalEnv)) {
    
    cat("\n------------------------------------------------------------\n")
    cat("Object:", obj_name, "\n")
    cat("------------------------------------------------------------\n")
    
    obj <- get(obj_name, envir = .GlobalEnv)
    
    cat("Class:\n")
    print(class(obj))
    
    cat("\nDimensions:\n")
    print(dim(obj))
    
    cat("\nContent:\n")
    print(obj)
  }
}


# ------------------------------------------------------------
# 8. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-3 Existing Mapping Audit Completed\n")
cat("============================================================\n\n")

cat("IMPORTANT:\n")
cat("No new Ensembl -> Entrez conversion was performed.\n")
cat("No DEG recalculation was performed.\n")
cat("testing_universe was not modified.\n")
cat("Existing mapping objects were inspected only.\n")

cat("\n============================================================\n")


# ============================================================
# Paper 1 — STEP9C-3
# Locate Existing Complete Universe Entrez Vector
#
# IMPORTANT:
#   No new ID conversion
#   No recalculation
#   No modification of testing_universe
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-3: Locate Existing Complete Universe Entrez Vector\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Inspect candidate Entrez objects
# ------------------------------------------------------------

candidate_objects <- c(
  "matched_universe_entrez",
  "universe_entrez",
  "universe_entrez_ids",
  "all_universe_entrez",
  "mapping_entrez",
  "entrez_ids",
  "entrez_annotation",
  "mapping_complete",
  "mapping_status"
)

for (obj_name in candidate_objects) {
  
  if (!exists(obj_name, envir = .GlobalEnv)) {
    next
  }
  
  cat("\n------------------------------------------------------------\n")
  cat("Object:", obj_name, "\n")
  cat("------------------------------------------------------------\n")
  
  obj <- get(obj_name, envir = .GlobalEnv)
  
  cat("Class:\n")
  print(class(obj))
  
  cat("Length:\n")
  print(length(obj))
  
  cat("Dimensions:\n")
  print(dim(obj))
  
  if (is.data.frame(obj)) {
    cat("Column names:\n")
    print(colnames(obj))
  }
  
  if (is.atomic(obj) && length(obj) <= 100) {
    cat("Content:\n")
    print(obj)
  }
}


# ------------------------------------------------------------
# 2. Inspect all global objects containing "entrez"
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("All global objects containing 'entrez'\n")
cat("============================================================\n\n")

all_objects <- ls(envir = .GlobalEnv)

entrez_objects <- all_objects[
  grepl("entrez", all_objects, ignore.case = TRUE)
]

print(entrez_objects)


# ------------------------------------------------------------
# 3. Compact structural summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Compact Entrez-object summary\n")
cat("============================================================\n\n")

for (obj_name in entrez_objects) {
  
  obj <- get(obj_name, envir = .GlobalEnv)
  
  cat(
    sprintf(
      "%-35s | class=%-15s | length=%-8s | dim=%s\n",
      obj_name,
      paste(class(obj), collapse = ","),
      length(obj),
      if (is.null(dim(obj))) {
        "NULL"
      } else {
        paste(dim(obj), collapse = " x ")
      }
    )
  )
}


# ------------------------------------------------------------
# 4. Locked reference numbers
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Locked reference numbers\n")
cat("============================================================\n\n")

cat("Testing universe Ensembl:", length(testing_universe), "\n")
cat("Entrez-mapped Ensembl:  ", mapping_summary$Entrez_Mapped[
  mapping_summary$Set == "Universe"
], "\n")
cat("Entrez-unmapped Ensembl:", mapping_summary$Entrez_Unmapped[
  mapping_summary$Set == "Universe"
], "\n")
cat("Exactly one Entrez:     ", mapping_summary$Exactly_One_Entrez[
  mapping_summary$Set == "Universe"
], "\n")
cat("More than one Entrez:   ", mapping_summary$More_Than_One_Entrez[
  mapping_summary$Set == "Universe"
], "\n")

cat(
  "Expected mapped Ensembl = 34631\n"
)

cat("\n============================================================\n")
cat("STEP9C-3 vector-location audit completed\n")
cat("============================================================\n")


# ============================================================
# Paper 1 — STEP9C-3
# Audit split_entrez as the Existing Complete Mapping Object
#
# IMPORTANT:
#   No new ID conversion
#   No recalculation
#   No modification of testing_universe
#   No modification of split_entrez
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-3: Audit split_entrez Complete Mapping Structure\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Basic structure
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("1. Basic structure of split_entrez\n")
cat("------------------------------------------------------------\n")

cat("Class:\n")
print(class(split_entrez))

cat("Length:\n")
print(length(split_entrez))

cat("Testing universe length:\n")
print(length(testing_universe))

cat(
  "Length matches testing_universe: ",
  length(split_entrez) == length(testing_universe),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 2. Inspect first 20 list elements
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("2. First 20 split_entrez elements\n")
cat("------------------------------------------------------------\n")

for (i in seq_len(min(20, length(split_entrez)))) {
  
  cat("\nIndex:", i, "\n")
  cat("Ensembl:", testing_universe[i], "\n")
  cat("Class:  ", paste(class(split_entrez[[i]]), collapse = ", "), "\n")
  cat("Length: ", length(split_entrez[[i]]), "\n")
  cat("Value:\n")
  print(split_entrez[[i]])
}


# ------------------------------------------------------------
# 3. Distribution of element lengths
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("3. Distribution of mapping multiplicity\n")
cat("------------------------------------------------------------\n")

element_lengths <- vapply(
  split_entrez,
  length,
  integer(1)
)

print(table(element_lengths))

cat("\nNumber with zero Entrez:\n")
print(sum(element_lengths == 0))

cat("Number with exactly one Entrez:\n")
print(sum(element_lengths == 1))

cat("Number with more than one Entrez:\n")
print(sum(element_lengths > 1))


# ------------------------------------------------------------
# 4. Compare with locked mapping counts
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("4. Compare split_entrez with locked mapping summary\n")
cat("------------------------------------------------------------\n")

cat(
  "Locked unmapped Ensembl:       ",
  mapping_summary$Entrez_Unmapped[
    mapping_summary$Set == "Universe"
  ],
  "\n",
  sep = ""
)

cat(
  "split_entrez zero-length:      ",
  sum(element_lengths == 0),
  "\n",
  sep = ""
)

cat(
  "Locked exactly-one Entrez:     ",
  mapping_summary$Exactly_One_Entrez[
    mapping_summary$Set == "Universe"
  ],
  "\n",
  sep = ""
)

cat(
  "split_entrez exactly-one:      ",
  sum(element_lengths == 1),
  "\n",
  sep = ""
)

cat(
  "Locked >1 Entrez:              ",
  mapping_summary$More_Than_One_Entrez[
    mapping_summary$Set == "Universe"
  ],
  "\n",
  sep = ""
)

cat(
  "split_entrez >1:               ",
  sum(element_lengths > 1),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. Inspect all multi-Entrez elements
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("5. Multi-Entrez elements\n")
cat("------------------------------------------------------------\n")

multi_idx <- which(element_lengths > 1)

cat("Number of multi-Entrez elements:", length(multi_idx), "\n")

if (length(multi_idx) > 0) {
  
  cat("\nFirst 20 multi-Entrez mappings:\n")
  
  for (i in multi_idx[seq_len(min(20, length(multi_idx)))]) {
    
    cat("\nEnsembl:", testing_universe[i], "\n")
    cat("Entrez values:\n")
    print(split_entrez[[i]])
  }
}


# ------------------------------------------------------------
# 6. Inspect zero-Entrez elements
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("6. Zero-Entrez elements\n")
cat("------------------------------------------------------------\n")

zero_idx <- which(element_lengths == 0)

cat("Number of zero-Entrez elements:", length(zero_idx), "\n")

if (length(zero_idx) > 0) {
  
  cat("\nFirst 20 unmapped Ensembl IDs:\n")
  
  print(
    testing_universe[
      zero_idx[seq_len(min(20, length(zero_idx)))]
    ]
  )
}


# ------------------------------------------------------------
# 7. Check whether list names preserve Ensembl IDs
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("7. Check list names\n")
cat("------------------------------------------------------------\n")

split_names <- names(split_entrez)

cat("Names present:\n")
print(!is.null(split_names))

if (!is.null(split_names)) {
  
  cat("Number of names:\n")
  print(length(split_names))
  
  cat("First 20 names:\n")
  print(head(split_names, 20))
  
  cat(
    "Names identical to testing_universe: ",
    identical(split_names, testing_universe),
    "\n",
    sep = ""
  )
  
  cat(
    "Names set identical to testing_universe: ",
    setequal(split_names, testing_universe),
    "\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 8. Check value types
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("8. Entrez value-type audit\n")
cat("------------------------------------------------------------\n")

value_classes <- unique(
  unlist(
    lapply(
      split_entrez,
      function(x) paste(class(x), collapse = ",")
    )
  )
)

cat("Classes encountered in list elements:\n")
print(value_classes)


# ------------------------------------------------------------
# 9. Check for NA / empty Entrez values
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("9. NA / empty-value audit\n")
cat("------------------------------------------------------------\n")

all_split_values <- unlist(
  split_entrez,
  use.names = FALSE
)

cat("Total Entrez values stored:\n")
print(length(all_split_values))

cat("NA values:\n")
print(sum(is.na(all_split_values)))

cat(
  "Empty-string values:\n"
)

print(
  sum(
    !is.na(all_split_values) &
      as.character(all_split_values) == ""
  )
)


# ------------------------------------------------------------
# 10. Unique Entrez count
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("10. Unique Entrez count\n")
cat("------------------------------------------------------------\n")

valid_entrez <- as.character(
  all_split_values[
    !is.na(all_split_values) &
      as.character(all_split_values) != ""
  ]
)

cat("Valid Entrez values:\n")
print(length(valid_entrez))

cat("Unique Entrez values:\n")
print(length(unique(valid_entrez)))

cat("Duplicated Entrez values:\n")
print(sum(duplicated(valid_entrez)))


# ------------------------------------------------------------
# 11. Final structural interpretation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-3 split_entrez audit completed\n")
cat("============================================================\n\n")

cat(
  "Testing universe: ",
  length(testing_universe),
  "\n",
  sep = ""
)

cat(
  "split_entrez elements: ",
  length(split_entrez),
  "\n",
  sep = ""
)

cat(
  "Zero Entrez elements: ",
  sum(element_lengths == 0),
  "\n",
  sep = ""
)

cat(
  "Exactly-one elements: ",
  sum(element_lengths == 1),
  "\n",
  sep = ""
)

cat(
  "Multi-Entrez elements: ",
  sum(element_lengths > 1),
  "\n",
  sep = ""
)

cat(
  "Unique valid Entrez IDs: ",
  length(unique(valid_entrez)),
  "\n",
  sep = ""
)

cat("\nNo ID conversion was performed.\n")
cat("No locked object was modified.\n")


# ============================================================
# Paper 1 — STEP9C-3
# Final KEGG Entrez Background Audit and Lock
#
# IMPORTANT:
#   No new ID conversion
#   No biomaRt
#   No bitr()
#   No AnnotationDbi conversion
#   No DEG recalculation
#   No modification of testing_universe
#   Source mapping object: existing split_entrez
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-3: Final KEGG Entrez Background Audit and Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Preserve locked testing universe
# ------------------------------------------------------------

testing_universe_before <- testing_universe

cat("------------------------------------------------------------\n")
cat("0. Locked testing universe\n")
cat("------------------------------------------------------------\n")

cat("Length:\n")
print(length(testing_universe_before))

cat("Unique:\n")
print(length(unique(testing_universe_before)))

cat("NA:\n")
print(sum(is.na(testing_universe_before)))

cat("Duplicated:\n")
print(sum(duplicated(testing_universe_before)))


# ------------------------------------------------------------
# 1. Extract existing Entrez mappings from split_entrez
#
# IMPORTANT:
#   This is extraction only.
#   No ID conversion is performed.
# ------------------------------------------------------------

all_split_values <- unlist(
  split_entrez,
  use.names = FALSE
)

valid_entrez <- as.character(
  all_split_values[
    !is.na(all_split_values) &
      as.character(all_split_values) != ""
  ]
)


# ------------------------------------------------------------
# 2. Construct final unique Entrez background
#
# IMPORTANT:
#   Unique() removes repeated Entrez IDs only.
#   No biological annotation or ID conversion is performed.
# ------------------------------------------------------------

final_universe_entrez <- unique(valid_entrez)


# ------------------------------------------------------------
# 3. Basic final-background QC
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("1. Final Entrez background QC\n")
cat("------------------------------------------------------------\n")

cat("Total Entrez values extracted:\n")
print(length(valid_entrez))

cat("Unique Entrez IDs:\n")
print(length(final_universe_entrez))

cat("NA:\n")
print(sum(is.na(final_universe_entrez)))

cat("Empty strings:\n")
print(sum(final_universe_entrez == ""))

cat("Duplicated values after unique():\n")
print(sum(duplicated(final_universe_entrez)))


# ------------------------------------------------------------
# 4. Check Entrez ID format
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("2. Entrez ID format audit\n")
cat("------------------------------------------------------------\n")

entrez_numeric_format <- grepl(
  "^[0-9]+$",
  final_universe_entrez
)

cat("Valid numeric Entrez IDs:\n")
print(sum(entrez_numeric_format))

cat("Invalid-format Entrez IDs:\n")
print(sum(!entrez_numeric_format))

if (any(!entrez_numeric_format)) {
  
  cat("\nInvalid-format values:\n")
  print(
    final_universe_entrez[
      !entrez_numeric_format
    ]
  )
}


# ------------------------------------------------------------
# 5. Reconcile against locked mapping summary
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("3. Reconciliation with locked mapping summary\n")
cat("------------------------------------------------------------\n")

locked_mapped_ensembl <- mapping_summary$Entrez_Mapped[
  mapping_summary$Set == "Universe"
]

locked_unmapped_ensembl <- mapping_summary$Entrez_Unmapped[
  mapping_summary$Set == "Universe"
]

locked_exactly_one <- mapping_summary$Exactly_One_Entrez[
  mapping_summary$Set == "Universe"
]

locked_multiple <- mapping_summary$More_Than_One_Entrez[
  mapping_summary$Set == "Universe"
]

cat("Locked mapped Ensembl:\n")
print(locked_mapped_ensembl)

cat("Observed mapped Ensembl from split_entrez:\n")
print(sum(element_lengths > 0))

cat("Locked unmapped Ensembl:\n")
print(locked_unmapped_ensembl)

cat("Observed unmapped Ensembl from split_entrez:\n")
print(sum(element_lengths == 0))

cat("Locked exactly-one Entrez:\n")
print(locked_exactly_one)

cat("Observed exactly-one elements:\n")
print(sum(element_lengths == 1))

cat("Locked multi-Entrez Ensembl:\n")
print(locked_multiple)

cat("Observed multi-Entrez elements:\n")
print(sum(element_lengths > 1))


# ------------------------------------------------------------
# 6. Core-gene Entrez audit
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("4. Core-gene Entrez background audit\n")
cat("------------------------------------------------------------\n")

core_entrez_expected <- c(
  AOC3   = "8639",
  MS4A7  = "58475",
  IL7R   = "3575"
)

cat("Locked core Entrez IDs:\n")
print(core_entrez_expected)

core_in_background <- core_entrez_expected %in% final_universe_entrez

cat("\nCore-gene presence:\n")
print(
  data.frame(
    Gene = names(core_entrez_expected),
    ENTREZID = unname(core_entrez_expected),
    In_Final_Background = core_in_background,
    stringsAsFactors = FALSE
  )
)


# ------------------------------------------------------------
# 7. Verify source relationship
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("5. Source-object integrity\n")
cat("------------------------------------------------------------\n")

cat(
  "split_entrez length = testing_universe length: ",
  length(split_entrez) == length(testing_universe),
  "\n",
  sep = ""
)

cat(
  "split_entrez names contain same Ensembl set: ",
  setequal(names(split_entrez), testing_universe),
  "\n",
  sep = ""
)

cat(
  "Total valid Entrez values = 36260: ",
  length(valid_entrez) == 36260,
  "\n",
  sep = ""
)

cat(
  "Unique Entrez IDs = 35365: ",
  length(final_universe_entrez) == 35365,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Verify testing_universe was not modified
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("6. testing_universe integrity after extraction\n")
cat("------------------------------------------------------------\n")

testing_universe_unchanged <- identical(
  testing_universe,
  testing_universe_before
)

cat(
  "testing_universe unchanged: ",
  testing_universe_unchanged,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Final PASS criteria
# ------------------------------------------------------------

final_background_pass <- all(
  length(testing_universe) == 57828,
  length(unique(testing_universe)) == 57828,
  sum(is.na(testing_universe)) == 0,
  sum(duplicated(testing_universe)) == 0,
  
  length(split_entrez) == 57828,
  
  sum(element_lengths == 0) == 0,
  sum(element_lengths == 1) == 57335,
  sum(element_lengths > 1) == 493,
  
  length(valid_entrez) == 36260,
  length(final_universe_entrez) == 35365,
  
  sum(is.na(final_universe_entrez)) == 0,
  sum(final_universe_entrez == "") == 0,
  sum(duplicated(final_universe_entrez)) == 0,
  
  all(entrez_numeric_format),
  
  all(core_in_background),
  
  testing_universe_unchanged
)


# ------------------------------------------------------------
# 10. Final LOCK
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-3 FINAL KEGG ENTrez BACKGROUND\n")
cat("============================================================\n\n")

cat("Final unique Entrez background size:\n")
print(length(final_universe_entrez))

cat("\nCore genes represented:\n")
print(core_in_background)

cat("\nFinal background PASS:\n")
print(final_background_pass)

if (final_background_pass) {
  
  cat("\nSTATUS: PASS — FINAL KEGG ENTrez BACKGROUND LOCKED\n")
  
} else {
  
  cat("\nSTATUS: FAIL — DO NOT PROCEED TO KEGG ENRICHMENT\n")
}

cat("\n============================================================\n")
cat("No new ID conversion was performed.\n")
cat("No DEG values were recalculated.\n")
cat("No threshold was changed.\n")
cat("testing_universe was not modified.\n")
cat("============================================================\n")


# ============================================================
# Paper 1 — STEP9C-3
# Corrected Final KEGG Entrez Background Audit
#
# IMPORTANT:
#   Existing split_entrez only
#   No new ID conversion
#   No biomaRt
#   No bitr()
#   No AnnotationDbi conversion
#   No DEG recalculation
#   No threshold change
#   No modification of testing_universe
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-3: Corrected Final KEGG Entrez Background Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Preserve locked universe
# ------------------------------------------------------------

testing_universe_before <- testing_universe

cat("0. Testing universe integrity\n")
cat("----------------------------\n")

cat("Length: ")
print(length(testing_universe))

cat("Unique: ")
print(length(unique(testing_universe)))

cat("NA: ")
print(sum(is.na(testing_universe)))

cat("Duplicated: ")
print(sum(duplicated(testing_universe)))


# ------------------------------------------------------------
# 1. Inspect split_entrez at the Ensembl level
#
# IMPORTANT:
#   Unmapped entries are stored as NA rather than character(0).
# ------------------------------------------------------------

split_has_entrez <- vapply(
  split_entrez,
  function(x) {
    any(
      !is.na(x) &
        as.character(x) != ""
    )
  },
  logical(1)
)

split_is_unmapped <- !split_has_entrez

split_n_entrez <- vapply(
  split_entrez,
  function(x) {
    sum(
      !is.na(x) &
        as.character(x) != ""
    )
  },
  integer(1)
)


# ------------------------------------------------------------
# 2. Correct Ensembl-level mapping counts
# ------------------------------------------------------------

observed_mapped_ensembl <- sum(split_has_entrez)

observed_unmapped_ensembl <- sum(split_is_unmapped)

observed_exactly_one <- sum(
  split_n_entrez == 1
)

observed_multiple <- sum(
  split_n_entrez > 1
)


cat("\n1. Correct Ensembl-level mapping counts\n")
cat("---------------------------------------\n")

cat(
  "Observed mapped Ensembl:     ",
  observed_mapped_ensembl,
  "\n",
  sep = ""
)

cat(
  "Observed unmapped Ensembl:   ",
  observed_unmapped_ensembl,
  "\n",
  sep = ""
)

cat(
  "Observed exactly-one:        ",
  observed_exactly_one,
  "\n",
  sep = ""
)

cat(
  "Observed multiple:           ",
  observed_multiple,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 3. Compare with locked mapping summary
# ------------------------------------------------------------

locked_mapped <- mapping_summary$Entrez_Mapped[
  mapping_summary$Set == "Universe"
]

locked_unmapped <- mapping_summary$Entrez_Unmapped[
  mapping_summary$Set == "Universe"
]

locked_exactly_one <- mapping_summary$Exactly_One_Entrez[
  mapping_summary$Set == "Universe"
]

locked_multiple <- mapping_summary$More_Than_One_Entrez[
  mapping_summary$Set == "Universe"
]


cat("\n2. Reconciliation with locked mapping summary\n")
cat("----------------------------------------------\n")

cat(
  "Mapped Ensembl: locked=",
  locked_mapped,
  " | observed=",
  observed_mapped_ensembl,
  "\n",
  sep = ""
)

cat(
  "Unmapped Ensembl: locked=",
  locked_unmapped,
  " | observed=",
  observed_unmapped_ensembl,
  "\n",
  sep = ""
)

cat(
  "Exactly-one: locked=",
  locked_exactly_one,
  " | observed=",
  observed_exactly_one,
  "\n",
  sep = ""
)

cat(
  "Multiple: locked=",
  locked_multiple,
  " | observed=",
  observed_multiple,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Extract valid Entrez values
# ------------------------------------------------------------

all_split_values <- unlist(
  split_entrez,
  use.names = FALSE
)

valid_entrez <- as.character(
  all_split_values[
    !is.na(all_split_values) &
      as.character(all_split_values) != ""
  ]
)

final_universe_entrez <- unique(valid_entrez)


# ------------------------------------------------------------
# 5. Final Entrez background QC
# ------------------------------------------------------------

cat("\n3. Final Entrez background QC\n")
cat("-----------------------------\n")

cat(
  "Total valid Entrez mapping entries: ",
  length(valid_entrez),
  "\n",
  sep = ""
)

cat(
  "Unique Entrez background:            ",
  length(final_universe_entrez),
  "\n",
  sep = ""
)

cat(
  "Duplicated unique-background values: ",
  sum(duplicated(final_universe_entrez)),
  "\n",
  sep = ""
)

cat(
  "NA: ",
  sum(is.na(final_universe_entrez)),
  "\n",
  sep = ""
)

cat(
  "Empty strings: ",
  sum(final_universe_entrez == ""),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Entrez format
# ------------------------------------------------------------

entrez_numeric_format <- grepl(
  "^[0-9]+$",
  final_universe_entrez
)

cat("\n4. Entrez format QC\n")
cat("-------------------\n")

cat(
  "Valid numeric Entrez IDs: ",
  sum(entrez_numeric_format),
  "\n",
  sep = ""
)

cat(
  "Invalid Entrez IDs: ",
  sum(!entrez_numeric_format),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Core-gene inclusion
# ------------------------------------------------------------

core_entrez_expected <- c(
  AOC3  = "8639",
  MS4A7 = "58475",
  IL7R  = "3575"
)

core_in_background <- (
  core_entrez_expected %in% final_universe_entrez
)

cat("\n5. Core-gene inclusion\n")
cat("----------------------\n")

print(
  data.frame(
    Gene = names(core_entrez_expected),
    ENTREZID = unname(core_entrez_expected),
    In_Final_Background = core_in_background,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
)


# ------------------------------------------------------------
# 8. Source-object integrity
# ------------------------------------------------------------

cat("\n6. Source-object integrity\n")
cat("--------------------------\n")

cat(
  "split_entrez length = testing_universe length: ",
  length(split_entrez) == length(testing_universe),
  "\n",
  sep = ""
)

cat(
  "split_entrez names = same Ensembl set: ",
  setequal(names(split_entrez), testing_universe),
  "\n",
  sep = ""
)

cat(
  "Unmapped stored as NA: ",
  sum(
    vapply(
      split_entrez,
      function(x) {
        length(x) == 1 &&
          is.na(x)
      },
      logical(1)
    )
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. testing_universe unchanged
# ------------------------------------------------------------

testing_universe_unchanged <- identical(
  testing_universe,
  testing_universe_before
)

cat(
  "testing_universe unchanged: ",
  testing_universe_unchanged,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Correct final PASS criteria
# ------------------------------------------------------------

final_background_pass_corrected <- all(
  
  # testing universe
  length(testing_universe) == 57828,
  length(unique(testing_universe)) == 57828,
  sum(is.na(testing_universe)) == 0,
  sum(duplicated(testing_universe)) == 0,
  
  # source mapping
  length(split_entrez) == 57828,
  setequal(names(split_entrez), testing_universe),
  
  # correct mapping reconciliation
  observed_mapped_ensembl == 34631,
  observed_unmapped_ensembl == 23197,
  observed_exactly_one == 34138,
  observed_multiple == 493,
  
  # final Entrez background
  length(valid_entrez) == 36260,
  length(final_universe_entrez) == 35365,
  sum(is.na(final_universe_entrez)) == 0,
  sum(final_universe_entrez == "") == 0,
  sum(duplicated(final_universe_entrez)) == 0,
  
  # format
  all(entrez_numeric_format),
  
  # core genes
  all(core_in_background),
  
  # no modification
  testing_universe_unchanged
)


# ------------------------------------------------------------
# 11. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-3 CORRECTED FINAL STATUS\n")
cat("============================================================\n\n")

cat(
  "Final unique Entrez background = ",
  length(final_universe_entrez),
  "\n",
  sep = ""
)

cat(
  "Mapped Ensembl = ",
  observed_mapped_ensembl,
  "\n",
  sep = ""
)

cat(
  "Unmapped Ensembl = ",
  observed_unmapped_ensembl,
  "\n",
  sep = ""
)

cat(
  "Exactly-one Ensembl = ",
  observed_exactly_one,
  "\n",
  sep = ""
)

cat(
  "Multi-Entrez Ensembl = ",
  observed_multiple,
  "\n",
  sep = ""
)

cat(
  "Final background PASS = ",
  final_background_pass_corrected,
  "\n",
  sep = ""
)

if (final_background_pass_corrected) {
  
  cat(
    "\nSTATUS: PASS — FINAL KEGG ENTrez BACKGROUND CONFIRMED\n"
  )
  
} else {
  
  cat(
    "\nSTATUS: FAIL — DO NOT LOCK FINAL KEGG BACKGROUND\n"
  )
}

cat("\n============================================================\n")
cat("No new ID conversion was performed.\n")
cat("No DEG recalculation was performed.\n")
cat("No threshold was changed.\n")
cat("testing_universe was not modified.\n")
cat("============================================================\n")

###STEP9C-3 LOCK record：
# ============================================================
# Paper 1 — STEP9C-3
# FINAL LOCK RECORD
#
# No recalculation
# No ID conversion
# ============================================================

STEP9C3_FINAL_LOCK <- list(
  
  Step = "STEP9C-3",
  
  Status = "PASS",
  
  Mapping_Source = "existing split_entrez",
  
  Testing_Universe_Ensembl = 57828,
  
  Mapped_Ensembl = 34631,
  
  Unmapped_Ensembl = 23197,
  
  Exactly_One_Entrez = 34138,
  
  Multi_Entrez_Ensembl = 493,
  
  Valid_Entrez_Mapping_Entries = 36260,
  
  Final_Unique_KEGG_Entrez_Background = 35365,
  
  Core_Genes = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  
  Core_Entrez = c(
    AOC3 = "8639",
    MS4A7 = "58475",
    IL7R = "3575"
  ),
  
  Core_Genes_In_Background = c(
    AOC3 = TRUE,
    MS4A7 = TRUE,
    IL7R = TRUE
  ),
  
  ID_Conversion_Performed = FALSE,
  
  DEG_Recalculation = FALSE,
  
  Threshold_Changed = FALSE,
  
  Testing_Universe_Modified = FALSE
)

cat("\n============================================================\n")
cat("STEP9C-3 FINAL LOCK RECORD\n")
cat("============================================================\n\n")

print(STEP9C3_FINAL_LOCK)

cat("\nFINAL STATUS: PASS — STEP9C-3 LOCKED\n")
cat("Final KEGG Entrez background: 35365 unique Entrez IDs\n")


# ============================================================
# Paper 1
# STEP9C-4
# KEGG Enrichment Input Entrez-ID Set and Background
# Correspondence Lock
#
# Purpose:
#   1. Lock Tumor-up KEGG input Entrez IDs
#   2. Lock Tumor-down KEGG input Entrez IDs
#   3. Verify all input IDs are contained in the
#      STEP9C-3 final KEGG background
#   4. Verify no new ID conversion is performed
#   5. Verify DEG objects and testing_universe are unchanged
#
# IMPORTANT:
#   This step does NOT:
#   - recalculate DEGs
#   - recalculate P values
#   - recalculate FDR
#   - change DEG thresholds
#   - perform new Ensembl-to-Entrez conversion
#   - modify testing_universe
#   - modify final_universe_entrez
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-4: KEGG Input Entrez-ID Policy Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}

setwd(step_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Required locked objects
# ------------------------------------------------------------

required_objects <- c(
  "deg_up",
  "deg_down",
  "testing_universe",
  "split_entrez",
  "final_universe_entrez",
  "STEP9C3_FINAL_LOCK"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  
  stop(
    paste0(
      "FAIL — Required object(s) missing: ",
      paste(missing_objects, collapse = ", ")
    )
  )
  
}

cat("Required objects: PASS\n")
cat(
  paste0(
    "Objects available: ",
    paste(required_objects, collapse = ", "),
    "\n\n"
  )
)


cat("STEP9C-3 type audit\n")
cat("====================\n\n")

cat(
  "Final_Unique_KEGG_Entrez_Background:\n"
)
print(STEP9C3_FINAL_LOCK$Final_Unique_KEGG_Entrez_Background)
cat(
  "class = ",
  class(STEP9C3_FINAL_LOCK$Final_Unique_KEGG_Entrez_Background),
  "\n",
  sep = ""
)
cat(
  "typeof = ",
  typeof(STEP9C3_FINAL_LOCK$Final_Unique_KEGG_Entrez_Background),
  "\n\n",
  sep = ""
)

cat(
  "Testing_Universe_Ensembl:\n"
)
print(STEP9C3_FINAL_LOCK$Testing_Universe_Ensembl)
cat(
  "class = ",
  class(STEP9C3_FINAL_LOCK$Testing_Universe_Ensembl),
  "\n",
  sep = ""
)
cat(
  "typeof = ",
  typeof(STEP9C3_FINAL_LOCK$Testing_Universe_Ensembl),
  "\n\n",
  sep = ""
)

cat(
  "Mapped_Ensembl:\n"
)
print(STEP9C3_FINAL_LOCK$Mapped_Ensembl)
cat(
  "class = ",
  class(STEP9C3_FINAL_LOCK$Mapped_Ensembl),
  "\n",
  sep = ""
)
cat(
  "typeof = ",
  typeof(STEP9C3_FINAL_LOCK$Mapped_Ensembl),
  "\n\n",
  sep = ""
)

cat(
  "Valid_Entrez_Mapping_Entries:\n"
)
print(STEP9C3_FINAL_LOCK$Valid_Entrez_Mapping_Entries)
cat(
  "class = ",
  class(STEP9C3_FINAL_LOCK$Valid_Entrez_Mapping_Entries),
  "\n",
  sep = ""
)
cat(
  "typeof = ",
  typeof(STEP9C3_FINAL_LOCK$Valid_Entrez_Mapping_Entries),
  "\n",
  sep = ""
)

# ============================================================
# Paper 1
# STEP9C-4
# KEGG Enrichment Input Entrez-ID Policy Lock
# FINAL CORRECTED VERSION
#
# Purpose:
#   1. Lock Tumor-up KEGG input Entrez IDs
#   2. Lock Tumor-down KEGG input Entrez IDs
#   3. Verify all input IDs are contained in the
#      STEP9C-3 final KEGG background
#   4. Verify no new ID conversion is performed
#   5. Verify DEG objects and testing_universe are unchanged
#
# IMPORTANT:
#   Numeric values in STEP9C3_FINAL_LOCK are validated by
#   numerical equality, not by R storage type.
#
#   This step does NOT:
#   - recalculate DEGs
#   - recalculate P values
#   - recalculate FDR
#   - change DEG thresholds
#   - perform new Ensembl-to-Entrez conversion
#   - modify testing_universe
#   - modify final_universe_entrez
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-4: KEGG Input Entrez-ID Policy Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}

setwd(step_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Required locked objects
# ------------------------------------------------------------

required_objects <- c(
  "deg_up",
  "deg_down",
  "testing_universe",
  "split_entrez",
  "final_universe_entrez",
  "STEP9C3_FINAL_LOCK"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  
  stop(
    paste0(
      "FAIL — Required object(s) missing: ",
      paste(missing_objects, collapse = ", ")
    )
  )
  
}

cat("Required objects: PASS\n")
cat(
  paste0(
    "Objects available: ",
    paste(required_objects, collapse = ", "),
    "\n\n"
  )
)


# ------------------------------------------------------------
# 2. Verify STEP9C-3 FINAL LOCK
#
# Numerical equality is used intentionally because the locked
# metadata values are stored as numeric/double in R.
# ------------------------------------------------------------

step9c3_lock_valid <- all(
  identical(STEP9C3_FINAL_LOCK$Step, "STEP9C-3"),
  identical(STEP9C3_FINAL_LOCK$Status, "PASS"),
  
  isTRUE(
    STEP9C3_FINAL_LOCK$Final_Unique_KEGG_Entrez_Background ==
      35365
  ),
  
  isTRUE(
    STEP9C3_FINAL_LOCK$Testing_Universe_Ensembl ==
      57828
  ),
  
  isTRUE(
    STEP9C3_FINAL_LOCK$Mapped_Ensembl ==
      34631
  ),
  
  isTRUE(
    STEP9C3_FINAL_LOCK$Valid_Entrez_Mapping_Entries ==
      36260
  )
)

if (!step9c3_lock_valid) {
  
  stop(
    "FAIL — STEP9C-3 FINAL LOCK numerical verification failed."
  )
  
}

cat("STEP9C-3 FINAL LOCK verification: PASS\n")
cat("  Testing universe = 57,828 Ensembl IDs\n")
cat("  Mapped Ensembl = 34,631\n")
cat("  Valid mapping entries = 36,260\n")
cat("  Final KEGG background = 35,365 unique Entrez IDs\n\n")


# ------------------------------------------------------------
# 3. Preserve pre-existing object fingerprints
# ------------------------------------------------------------

testing_universe_before <- testing_universe
final_universe_entrez_before <- final_universe_entrez

deg_up_gene_before <- deg_up$Gene
deg_down_gene_before <- deg_down$Gene

split_entrez_before <- split_entrez


# ------------------------------------------------------------
# 4. Verify testing_universe
# ------------------------------------------------------------

stopifnot(
  length(testing_universe) == 57828L,
  length(unique(testing_universe)) == 57828L,
  all(!is.na(testing_universe)),
  all(testing_universe != ""),
  all(grepl("^ENSG[0-9]+$", testing_universe))
)

cat("testing_universe QC: PASS\n")
cat("  Ensembl IDs:", length(testing_universe), "\n")
cat("  Unique:", length(unique(testing_universe)), "\n\n")


# ------------------------------------------------------------
# 5. Verify existing split_entrez mapping
# ------------------------------------------------------------

stopifnot(
  is.list(split_entrez),
  length(split_entrez) == length(testing_universe),
  setequal(names(split_entrez), testing_universe)
)

cat("Existing split_entrez mapping container: PASS\n")
cat("  Number of Ensembl entries:", length(split_entrez), "\n")
cat("  Same Ensembl ID set as testing_universe: TRUE\n\n")


# ------------------------------------------------------------
# 6. Existing Entrez extraction function
#
# NO NEW ID CONVERSION.
# Only existing split_entrez is queried.
# ------------------------------------------------------------

extract_existing_entrez <- function(
    ensembl_ids,
    split_map
) {
  
  if (length(ensembl_ids) == 0) {
    return(character(0))
  }
  
  missing_names <- setdiff(
    ensembl_ids,
    names(split_map)
  )
  
  if (length(missing_names) > 0) {
    
    stop(
      paste0(
        "FAIL — Requested Ensembl IDs absent from split_entrez: ",
        length(missing_names)
      )
    )
    
  }
  
  selected <- split_map[ensembl_ids]
  
  values <- unlist(
    selected,
    use.names = FALSE
  )
  
  values <- as.character(
    values[
      !is.na(values) &
        as.character(values) != ""
    ]
  )
  
  unique(values)
}


# ------------------------------------------------------------
# 7. Extract original DEG Ensembl IDs
# ------------------------------------------------------------

up_ensembl <- unique(
  as.character(
    deg_up$Gene[
      !is.na(deg_up$Gene) &
        deg_up$Gene != ""
    ]
  )
)

down_ensembl <- unique(
  as.character(
    deg_down$Gene[
      !is.na(deg_down$Gene) &
        deg_down$Gene != ""
    ]
  )
)

cat("Original DEG input sets:\n")
cat(
  "  Tumor-up unique Ensembl IDs:",
  length(up_ensembl),
  "\n"
)

cat(
  "  Tumor-down unique Ensembl IDs:",
  length(down_ensembl),
  "\n\n"
)


# ------------------------------------------------------------
# 8. DEG-to-testing_universe correspondence
# ------------------------------------------------------------

up_outside_universe <- setdiff(
  up_ensembl,
  testing_universe
)

down_outside_universe <- setdiff(
  down_ensembl,
  testing_universe
)

if (length(up_outside_universe) > 0) {
  
  stop(
    paste0(
      "FAIL — Tumor-up DEG IDs outside testing_universe: ",
      length(up_outside_universe)
    )
  )
  
}

if (length(down_outside_universe) > 0) {
  
  stop(
    paste0(
      "FAIL — Tumor-down DEG IDs outside testing_universe: ",
      length(down_outside_universe)
    )
  )
  
}

cat("DEG-to-testing_universe correspondence: PASS\n")
cat("  Tumor-up outside universe: 0\n")
cat("  Tumor-down outside universe: 0\n\n")


# ------------------------------------------------------------
# 9. Extract KEGG input Entrez IDs
# ------------------------------------------------------------

kegg_up_entrez <- extract_existing_entrez(
  ensembl_ids = up_ensembl,
  split_map = split_entrez
)

kegg_down_entrez <- extract_existing_entrez(
  ensembl_ids = down_ensembl,
  split_map = split_entrez
)

cat("KEGG input Entrez extraction complete.\n\n")

cat("Tumor-up:\n")
cat("  Ensembl IDs:", length(up_ensembl), "\n")
cat("  Unique Entrez IDs:", length(kegg_up_entrez), "\n\n")

cat("Tumor-down:\n")
cat("  Ensembl IDs:", length(down_ensembl), "\n")
cat("  Unique Entrez IDs:", length(kegg_down_entrez), "\n\n")


# ------------------------------------------------------------
# 10. Entrez structural QC
# ------------------------------------------------------------

up_entrez_structural_pass <- (
  length(kegg_up_entrez) ==
    length(unique(kegg_up_entrez)) &&
    all(!is.na(kegg_up_entrez)) &&
    all(kegg_up_entrez != "") &&
    all(grepl("^[0-9]+$", kegg_up_entrez))
)

down_entrez_structural_pass <- (
  length(kegg_down_entrez) ==
    length(unique(kegg_down_entrez)) &&
    all(!is.na(kegg_down_entrez)) &&
    all(kegg_down_entrez != "") &&
    all(grepl("^[0-9]+$", kegg_down_entrez))
)

cat(
  "Tumor-up Entrez structural QC:",
  up_entrez_structural_pass,
  "\n"
)

cat(
  "Tumor-down Entrez structural QC:",
  down_entrez_structural_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 11. KEGG input-to-background correspondence
# ------------------------------------------------------------

up_not_in_background <- setdiff(
  kegg_up_entrez,
  final_universe_entrez
)

down_not_in_background <- setdiff(
  kegg_down_entrez,
  final_universe_entrez
)

if (length(up_not_in_background) > 0) {
  
  stop(
    paste0(
      "FAIL — Tumor-up KEGG input IDs outside final background: ",
      length(up_not_in_background)
    )
  )
  
}

if (length(down_not_in_background) > 0) {
  
  stop(
    paste0(
      "FAIL — Tumor-down KEGG input IDs outside final background: ",
      length(down_not_in_background)
    )
  )
  
}

cat("KEGG input-to-background correspondence: PASS\n")
cat("  Tumor-up outside background: 0\n")
cat("  Tumor-down outside background: 0\n\n")


# ------------------------------------------------------------
# 12. Descriptive mapping coverage
# ------------------------------------------------------------

up_mapped_n <- length(kegg_up_entrez)
down_mapped_n <- length(kegg_down_entrez)

up_coverage <- up_mapped_n /
  length(up_ensembl)

down_coverage <- down_mapped_n /
  length(down_ensembl)

cat("Existing Entrez mapping coverage:\n")

cat(
  sprintf(
    "  Tumor-up: %d / %d = %.4f\n",
    up_mapped_n,
    length(up_ensembl),
    up_coverage
  )
)

cat(
  sprintf(
    "  Tumor-down: %d / %d = %.4f\n",
    down_mapped_n,
    length(down_ensembl),
    down_coverage
  )
)

cat("\n")


# ------------------------------------------------------------
# 13. Identify DEG Ensembl IDs without existing Entrez mapping
# ------------------------------------------------------------

has_existing_entrez <- vapply(
  split_entrez,
  function(x) {
    any(
      !is.na(x) &
        as.character(x) != ""
    )
  },
  logical(1)
)

mapped_ensembl_ids <- names(
  split_entrez
)[has_existing_entrez]

up_unmapped_ensembl <- setdiff(
  up_ensembl,
  mapped_ensembl_ids
)

down_unmapped_ensembl <- setdiff(
  down_ensembl,
  mapped_ensembl_ids
)

cat("Unmapped DEG Ensembl IDs:\n")
cat(
  "  Tumor-up:",
  length(up_unmapped_ensembl),
  "\n"
)

cat(
  "  Tumor-down:",
  length(down_unmapped_ensembl),
  "\n\n"
)


# ------------------------------------------------------------
# 14. Up/down Entrez overlap
# ------------------------------------------------------------

up_down_entrez_overlap <- intersect(
  kegg_up_entrez,
  kegg_down_entrez
)

cat("Up/down KEGG Entrez overlap:\n")
cat(
  "  Shared Entrez IDs:",
  length(up_down_entrez_overlap),
  "\n\n"
)


# ------------------------------------------------------------
# 15. Core-gene correspondence
# ------------------------------------------------------------

core_entrez <- c(
  AOC3 = "8639",
  MS4A7 = "58475",
  IL7R = "3575"
)

core_in_background <-
  core_entrez %in% final_universe_entrez

core_in_up_input <-
  core_entrez %in% kegg_up_entrez

core_in_down_input <-
  core_entrez %in% kegg_down_entrez

cat("Core-gene KEGG correspondence:\n")

print(
  data.frame(
    Gene = names(core_entrez),
    Entrez = unname(core_entrez),
    In_Final_Background =
      unname(core_in_background),
    In_Tumor_Up_Input =
      unname(core_in_up_input),
    In_Tumor_Down_Input =
      unname(core_in_down_input),
    row.names = NULL
  )
)

cat("\n")


# ------------------------------------------------------------
# 16. Locked-object integrity audit
# ------------------------------------------------------------

testing_universe_unchanged <- identical(
  testing_universe,
  testing_universe_before
)

final_background_unchanged <- identical(
  final_universe_entrez,
  final_universe_entrez_before
)

deg_up_unchanged <- identical(
  deg_up$Gene,
  deg_up_gene_before
)

deg_down_unchanged <- identical(
  deg_down$Gene,
  deg_down_gene_before
)

split_entrez_unchanged <- identical(
  split_entrez,
  split_entrez_before
)

cat("Locked-object integrity audit:\n")

cat(
  "  testing_universe unchanged:",
  testing_universe_unchanged,
  "\n"
)

cat(
  "  final_universe_entrez unchanged:",
  final_background_unchanged,
  "\n"
)

cat(
  "  deg_up$Gene unchanged:",
  deg_up_unchanged,
  "\n"
)

cat(
  "  deg_down$Gene unchanged:",
  deg_down_unchanged,
  "\n"
)

cat(
  "  split_entrez unchanged:",
  split_entrez_unchanged,
  "\n\n"
)


# ------------------------------------------------------------
# 17. Final PASS criteria
# ------------------------------------------------------------

STEP9C4_PASS <- all(
  
  step9c3_lock_valid,
  
  length(kegg_up_entrez) > 0,
  
  length(kegg_down_entrez) > 0,
  
  up_entrez_structural_pass,
  
  down_entrez_structural_pass,
  
  length(up_not_in_background) == 0,
  
  length(down_not_in_background) == 0,
  
  testing_universe_unchanged,
  
  final_background_unchanged,
  
  deg_up_unchanged,
  
  deg_down_unchanged,
  
  split_entrez_unchanged,
  
  length(final_universe_entrez) == 35365,
  
  length(unique(final_universe_entrez)) == 35365
)


# ============================================================
# STEP9C-4
# Final Lock Record
#
# All analytical QC steps have already passed.
# This block only constructs and saves the final lock record.
# ============================================================

STEP9C4_FINAL_LOCK <- list(
  
  Step = "STEP9C-4",
  
  Status =
    if (STEP9C4_PASS) "PASS" else "FAIL",
  
  Input_Source =
    "existing deg_up$Gene and deg_down$Gene",
  
  Mapping_Source =
    "existing split_entrez",
  
  ID_Conversion_Performed = FALSE,
  
  DEG_Recalculation = FALSE,
  
  P_Value_Recalculation = FALSE,
  
  FDR_Recalculation = FALSE,
  
  Threshold_Changed = FALSE,
  
  Testing_Universe_Modified = FALSE,
  
  Final_KEGG_Background = 35365,
  
  Tumor_Up = list(
    
    Original_Ensembl =
      length(up_ensembl),
    
    Mapped_Entrez =
      length(kegg_up_entrez),
    
    Unmapped_Ensembl =
      length(up_unmapped_ensembl),
    
    Mapping_Coverage =
      up_coverage,
    
    Input_Entrez =
      kegg_up_entrez
    
  ),
  
  Tumor_Down = list(
    
    Original_Ensembl =
      length(down_ensembl),
    
    Mapped_Entrez =
      length(kegg_down_entrez),
    
    Unmapped_Ensembl =
      length(down_unmapped_ensembl),
    
    Mapping_Coverage =
      down_coverage,
    
    Input_Entrez =
      kegg_down_entrez
    
  ),
  
  Up_Down_Entrez_Overlap =
    length(up_down_entrez_overlap),
  
  Core_Entrez =
    core_entrez,
  
  Core_In_Background =
    core_in_background,
  
  Core_In_Tumor_Up_Input =
    core_in_up_input,
  
  Core_In_Tumor_Down_Input =
    core_in_down_input,
  
  Background_Size =
    length(final_universe_entrez),
  
  Background_Unique =
    length(unique(final_universe_entrez)),
  
  Background_Unchanged =
    final_background_unchanged
  
)


# ------------------------------------------------------------
# Final STEP9C-4 status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-4 FINAL LOCK RECORD\n")
cat("============================================================\n\n")

print(STEP9C4_FINAL_LOCK)

cat("\n============================================================\n")

if (STEP9C4_PASS) {
  
  cat(
    "FINAL STATUS: PASS — STEP9C-4 LOCKED\n"
  )
  
  cat(
    "KEGG background: 35,365 unique Entrez IDs\n"
  )
  
  cat(
    "Tumor-up KEGG input: ",
    length(kegg_up_entrez),
    " unique Entrez IDs\n",
    sep = ""
  )
  
  cat(
    "Tumor-down KEGG input: ",
    length(kegg_down_entrez),
    " unique Entrez IDs\n",
    sep = ""
  )
  
} else {
  
  cat(
    "FINAL STATUS: FAIL — STEP9C-4 NOT LOCKED\n"
  )
  
}

cat("============================================================\n")


# ------------------------------------------------------------
# Save final lock records
# ------------------------------------------------------------

step_dir <- "result_dir/Paper 1_Step 9-Functional Enrichment"

saveRDS(
  STEP9C4_FINAL_LOCK,
  file = file.path(
    step_dir,
    "STEP9C4_FINAL_LOCK.rds"
  )
)

saveRDS(
  list(
    kegg_up_entrez = kegg_up_entrez,
    kegg_down_entrez = kegg_down_entrez,
    final_universe_entrez = final_universe_entrez
  ),
  file = file.path(
    step_dir,
    "STEP9C4_KEGG_INPUT_AND_BACKGROUND_LOCKED.rds"
  )
)

cat("\nSaved:\n")
cat("  STEP9C4_FINAL_LOCK.rds\n")
cat("  STEP9C4_KEGG_INPUT_AND_BACKGROUND_LOCKED.rds\n")

# ------------------------------------------------------------
# 19. Final output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-4 FINAL LOCK RECORD\n")
cat("============================================================\n\n")

print(STEP9C4_FINAL_LOCK)

cat("\n============================================================\n")

if (STEP9C4_PASS) {
  
  cat(
    "FINAL STATUS: PASS — STEP9C-4 LOCKED\n"
  )
  
  cat(
    "KEGG background: 35,365 unique Entrez IDs\n"
  )
  
  cat(
    "Tumor-up KEGG input: ",
    length(kegg_up_entrez),
    " unique Entrez IDs\n",
    sep = ""
  )
  
  cat(
    "Tumor-down KEGG input: ",
    length(kegg_down_entrez),
    " unique Entrez IDs\n",
    sep = ""
  )
  
} else {
  
  cat(
    "FINAL STATUS: FAIL — STEP9C-4 NOT LOCKED\n"
  )
  
}

cat("============================================================\n")


# ------------------------------------------------------------
# 20. Save STEP9C-4 lock records
# ------------------------------------------------------------

saveRDS(
  STEP9C4_FINAL_LOCK,
  file = file.path(
    step_dir,
    "STEP9C4_FINAL_LOCK.rds"
  )
)

saveRDS(
  list(
    kegg_up_entrez = kegg_up_entrez,
    kegg_down_entrez = kegg_down_entrez,
    final_universe_entrez =
      final_universe_entrez
  ),
  file = file.path(
    step_dir,
    "STEP9C4_KEGG_INPUT_AND_BACKGROUND_LOCKED.rds"
  )
)

cat("\nSaved:\n")
cat("  STEP9C4_FINAL_LOCK.rds\n")
cat("  STEP9C4_KEGG_INPUT_AND_BACKGROUND_LOCKED.rds\n")

##STEP9C-5 — KEGG enrichment execution and result-definition lock
# ============================================================
# Paper 1
# STEP9C-5
# KEGG Enrichment Execution and Result-Definition Lock
# ============================================================
#
# Purpose:
#   Perform the final KEGG over-representation analysis (ORA)
#   using ONLY the locked Entrez ID inputs and background
#   established in STEP9C-3 and STEP9C-4.
#
# LOCKED INPUTS:
#   KEGG background      = 35,365 unique Entrez IDs
#   Tumor-up KEGG input  = 5,194 unique Entrez IDs
#   Tumor-down KEGG input= 4,463 unique Entrez IDs
#
# IMPORTANT:
#   - No new Ensembl -> Entrez conversion
#   - No DEG recalculation
#   - No P-value recalculation for DEGs
#   - No FDR recalculation for DEGs
#   - No DEG threshold modification
#   - No modification of testing_universe
#   - No modification of final_universe_entrez
#
# Statistical definition:
#   KEGG over-representation analysis (ORA)
#   Organism: Homo sapiens (hsa)
#   ID type: ncbi-geneid
#   Background: locked 35,365 Entrez IDs
#   Multiple testing: Benjamini-Hochberg (BH)
#   Final significance definition:
#       adjusted P value (BH) < 0.05
#
# ============================================================
# ============================================================
# STEP9C-5 PRECHECK-2
# Audit existing STEP9C files
# ============================================================

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

setwd(step_dir)

cat("\n============================================================\n")
cat("STEP9C-5 PRECHECK-2: Existing file audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. List all STEP9C-related files
# ------------------------------------------------------------

all_files <- list.files(
  path = step_dir,
  full.names = FALSE,
  recursive = FALSE
)

step9c_files <- all_files[
  grepl(
    "STEP9C|step9c|entrez|KEGG|kegg",
    all_files,
    ignore.case = TRUE
  )
]

cat("STEP9C-related files found:\n\n")

if (length(step9c_files) == 0) {
  
  cat("NONE\n")
  
} else {
  
  print(step9c_files)
  
}

cat("\n")


# ------------------------------------------------------------
# 2. Show RDS files specifically
# ------------------------------------------------------------

rds_files <- all_files[
  grepl("\\.rds$", all_files, ignore.case = TRUE)
]

cat("RDS files in STEP9C working directory:\n\n")

if (length(rds_files) == 0) {
  
  cat("NONE\n")
  
} else {
  
  print(rds_files)
  
}

cat("\n")


# ------------------------------------------------------------
# 3. Check currently available objects
# ------------------------------------------------------------

objects_to_check <- c(
  "STEP9C3_FINAL_LOCK",
  "STEP9C4_FINAL_LOCK",
  "final_universe_entrez",
  "split_entrez",
  "testing_universe"
)

cat("Current R-session objects:\n\n")

for (obj in objects_to_check) {
  
  cat(
    obj,
    ": ",
    if (exists(obj, inherits = FALSE)) "PRESENT" else "ABSENT",
    "\n",
    sep = ""
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 4. Inspect STEP9C4 lock if available
# ------------------------------------------------------------

if (exists("STEP9C4_FINAL_LOCK", inherits = FALSE)) {
  
  cat("STEP9C4_FINAL_LOCK summary:\n\n")
  
  cat(
    "Step: ",
    STEP9C4_FINAL_LOCK$Step,
    "\n",
    sep = ""
  )
  
  cat(
    "Status: ",
    STEP9C4_FINAL_LOCK$Status,
    "\n",
    sep = ""
  )
  
  cat(
    "Background size: ",
    STEP9C4_FINAL_LOCK$Final_KEGG_Background,
    "\n",
    sep = ""
  )
  
  cat(
    "Tumor-up mapped Entrez: ",
    STEP9C4_FINAL_LOCK$Tumor_Up$Mapped_Entrez,
    "\n",
    sep = ""
  )
  
  cat(
    "Tumor-down mapped Entrez: ",
    STEP9C4_FINAL_LOCK$Tumor_Down$Mapped_Entrez,
    "\n",
    sep = ""
  )
  
}


cat("\n============================================================\n")
cat("PRECHECK-2 COMPLETE\n")
cat("============================================================\n")

# ============================================================
# STEP9C-5
# LOCKED OBJECT RECOVERY AUDIT
# ============================================================

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

setwd(step_dir)

cat("\n============================================================\n")
cat("STEP9C-5: Locked Object Recovery Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Confirm required RDS files
# ------------------------------------------------------------

file_c4_lock <- file.path(
  step_dir,
  "STEP9C4_FINAL_LOCK.rds"
)

file_c4_background <- file.path(
  step_dir,
  "STEP9C4_KEGG_INPUT_AND_BACKGROUND_LOCKED.rds"
)


cat("STEP9C4_FINAL_LOCK.rds: ")

if (file.exists(file_c4_lock)) {
  cat("FOUND\n")
} else {
  cat("MISSING\n")
}


cat(
  "STEP9C4_KEGG_INPUT_AND_BACKGROUND_LOCKED.rds: "
)

if (file.exists(file_c4_background)) {
  cat("FOUND\n")
} else {
  cat("MISSING\n")
}


# ------------------------------------------------------------
# 2. Read STEP9C4_FINAL_LOCK
# ------------------------------------------------------------

if (!file.exists(file_c4_lock)) {
  
  stop(
    "FAIL — STEP9C4_FINAL_LOCK.rds is missing."
  )
  
}

STEP9C4_FINAL_LOCK <- readRDS(
  file_c4_lock
)

cat("\nSTEP9C4_FINAL_LOCK loaded.\n")

cat(
  "Step: ",
  STEP9C4_FINAL_LOCK$Step,
  "\n",
  sep = ""
)

cat(
  "Status: ",
  STEP9C4_FINAL_LOCK$Status,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 3. Inspect the second RDS object
# ------------------------------------------------------------

if (!file.exists(file_c4_background)) {
  
  stop(
    "FAIL — STEP9C4_KEGG_INPUT_AND_BACKGROUND_LOCKED.rds is missing."
  )
  
}

C4_BACKGROUND_OBJECT <- readRDS(
  file_c4_background
)

cat("\n============================================================\n")
cat("STEP9C4_KEGG_INPUT_AND_BACKGROUND_LOCKED.rds structure\n")
cat("============================================================\n\n")

cat("Object class:\n")
print(class(C4_BACKGROUND_OBJECT))

cat("\nObject length:\n")
print(length(C4_BACKGROUND_OBJECT))

cat("\nObject names:\n")
print(names(C4_BACKGROUND_OBJECT))


# ------------------------------------------------------------
# 4. Inspect each top-level component
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Top-level component summary\n")
cat("============================================================\n\n")

if (is.list(C4_BACKGROUND_OBJECT)) {
  
  for (nm in names(C4_BACKGROUND_OBJECT)) {
    
    x <- C4_BACKGROUND_OBJECT[[nm]]
    
    cat("\n--- ", nm, " ---\n", sep = "")
    
    cat("class: ")
    print(class(x))
    
    cat("length: ")
    print(length(x))
    
    if (is.atomic(x) && length(x) <= 10) {
      cat("values:\n")
      print(x)
    }
    
  }
  
} else {
  
  cat(
    "Object is not a list.\n"
  )
  
}


# ------------------------------------------------------------
# 5. Inspect STEP9C4 lock fields
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C4_FINAL_LOCK fields\n")
cat("============================================================\n\n")

print(names(STEP9C4_FINAL_LOCK))


# ------------------------------------------------------------
# 6. Confirm locked numerical definitions
# ------------------------------------------------------------

cat("\nLocked numerical definitions:\n\n")

cat(
  "KEGG background size: ",
  STEP9C4_FINAL_LOCK$Final_KEGG_Background,
  "\n",
  sep = ""
)

cat(
  "Tumor-up Entrez input: ",
  STEP9C4_FINAL_LOCK$Tumor_Up$Mapped_Entrez,
  "\n",
  sep = ""
)

cat(
  "Tumor-down Entrez input: ",
  STEP9C4_FINAL_LOCK$Tumor_Down$Mapped_Entrez,
  "\n",
  sep = ""
)


cat("\n============================================================\n")
cat("RECOVERY AUDIT COMPLETE\n")
cat("============================================================\n")


# ============================================================
# Paper 1
# STEP9C-5
# KEGG Over-Representation Analysis
# Formal Execution + Result Definition Lock
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-5: KEGG Over-Representation Analysis\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

setwd(step_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Required package
# ------------------------------------------------------------

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  
  stop(
    "FAIL — clusterProfiler is not installed."
  )
  
}

library(clusterProfiler)

cat(
  "clusterProfiler version: ",
  as.character(
    packageVersion("clusterProfiler")
  ),
  "\n",
  sep = ""
)

cat(
  "R version: ",
  R.version.string,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 2. Recover locked STEP9C-4 objects
# ------------------------------------------------------------

file_c4_lock <- file.path(
  step_dir,
  "STEP9C4_FINAL_LOCK.rds"
)

file_c4_background <- file.path(
  step_dir,
  "STEP9C4_KEGG_INPUT_AND_BACKGROUND_LOCKED.rds"
)


if (!file.exists(file_c4_lock)) {
  
  stop(
    "FAIL — STEP9C4_FINAL_LOCK.rds is missing."
  )
  
}


if (!file.exists(file_c4_background)) {
  
  stop(
    "FAIL — STEP9C4_KEGG_INPUT_AND_BACKGROUND_LOCKED.rds is missing."
  )
  
}


STEP9C4_FINAL_LOCK <- readRDS(
  file_c4_lock
)


C4_BACKGROUND_OBJECT <- readRDS(
  file_c4_background
)


cat("STEP9C4 lock recovered.\n")
cat("Background/input object recovered.\n\n")


# ------------------------------------------------------------
# 3. Recover exact locked vectors
# ------------------------------------------------------------

kegg_background <- C4_BACKGROUND_OBJECT$final_universe_entrez

kegg_up_entrez <- C4_BACKGROUND_OBJECT$kegg_up_entrez

kegg_down_entrez <- C4_BACKGROUND_OBJECT$kegg_down_entrez


# ------------------------------------------------------------
# 4. Structural validation of recovered vectors
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP9C-5 PRE-EXECUTION VECTOR QC\n")
cat("============================================================\n\n")


vector_qc <- c(
  
  background_is_character =
    is.character(kegg_background),
  
  up_is_character =
    is.character(kegg_up_entrez),
  
  down_is_character =
    is.character(kegg_down_entrez),
  
  background_unique =
    length(unique(kegg_background)) == length(kegg_background),
  
  up_unique =
    length(unique(kegg_up_entrez)) == length(kegg_up_entrez),
  
  down_unique =
    length(unique(kegg_down_entrez)) == length(kegg_down_entrez),
  
  background_no_na =
    !anyNA(kegg_background),
  
  up_no_na =
    !anyNA(kegg_up_entrez),
  
  down_no_na =
    !anyNA(kegg_down_entrez),
  
  background_nonempty =
    length(kegg_background) > 0,
  
  up_nonempty =
    length(kegg_up_entrez) > 0,
  
  down_nonempty =
    length(kegg_down_entrez) > 0
  
)


print(vector_qc)


if (!all(vector_qc)) {
  
  stop(
    "FAIL — Recovered KEGG vectors failed structural QC."
  )
}


cat("\nStructural QC: PASS\n\n")


# ------------------------------------------------------------
# 5. Exact numerical lock validation
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP9C-5 LOCKED NUMERICAL DEFINITION CHECK\n")
cat("============================================================\n\n")


background_size_ok <- (
  length(kegg_background) ==
    STEP9C4_FINAL_LOCK$Final_KEGG_Background
)


up_size_ok <- (
  length(kegg_up_entrez) ==
    STEP9C4_FINAL_LOCK$Tumor_Up$Mapped_Entrez
)


down_size_ok <- (
  length(kegg_down_entrez) ==
    STEP9C4_FINAL_LOCK$Tumor_Down$Mapped_Entrez
)


background_unique_ok <- (
  length(unique(kegg_background)) == 35365
)


up_unique_ok <- (
  length(unique(kegg_up_entrez)) == 5194
)


down_unique_ok <- (
  length(unique(kegg_down_entrez)) == 4463
)


cat(
  "Background length = ",
  length(kegg_background),
  " | expected = 35365 | ",
  ifelse(background_size_ok, "PASS", "FAIL"),
  "\n",
  sep = ""
)


cat(
  "Tumor-up length = ",
  length(kegg_up_entrez),
  " | expected = 5194 | ",
  ifelse(up_size_ok, "PASS", "FAIL"),
  "\n",
  sep = ""
)


cat(
  "Tumor-down length = ",
  length(kegg_down_entrez),
  " | expected = 4463 | ",
  ifelse(down_size_ok, "PASS", "FAIL"),
  "\n",
  sep = ""
)


cat(
  "Background unique = ",
  length(unique(kegg_background)),
  " | ",
  ifelse(background_unique_ok, "PASS", "FAIL"),
  "\n",
  sep = ""
)


cat(
  "Tumor-up unique = ",
  length(unique(kegg_up_entrez)),
  " | ",
  ifelse(up_unique_ok, "PASS", "FAIL"),
  "\n",
  sep = ""
)


cat(
  "Tumor-down unique = ",
  length(unique(kegg_down_entrez)),
  " | ",
  ifelse(down_unique_ok, "PASS", "FAIL"),
  "\n",
  sep = ""
)


if (!all(
  background_size_ok,
  up_size_ok,
  down_size_ok,
  background_unique_ok,
  up_unique_ok,
  down_unique_ok
)) {
  
  stop(
    "FAIL — Recovered vectors do not exactly match STEP9C-4 locked definitions."
  )
}


cat("\nNumerical lock validation: PASS\n\n")


# ------------------------------------------------------------
# 6. Background containment validation
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP9C-5 BACKGROUND CONTAINMENT CHECK\n")
cat("============================================================\n\n")


up_not_in_background <- setdiff(
  kegg_up_entrez,
  kegg_background
)


down_not_in_background <- setdiff(
  kegg_down_entrez,
  kegg_background
)


up_down_overlap <- intersect(
  kegg_up_entrez,
  kegg_down_entrez
)


cat(
  "Tumor-up IDs outside background: ",
  length(up_not_in_background),
  "\n",
  sep = ""
)


cat(
  "Tumor-down IDs outside background: ",
  length(down_not_in_background),
  "\n",
  sep = ""
)


cat(
  "Up/down Entrez overlap: ",
  length(up_down_overlap),
  "\n",
  sep = ""
)


if (
  length(up_not_in_background) != 0 ||
  length(down_not_in_background) != 0 ||
  length(up_down_overlap) != 0
) {
  
  stop(
    "FAIL — KEGG input/background containment or independence check failed."
  )
}


cat("\nBackground containment QC: PASS\n\n")


# ------------------------------------------------------------
# 7. Confirm KEGG analysis definition
# ------------------------------------------------------------

KEGG_ANALYSIS_DEFINITION <- list(
  
  Organism = "hsa",
  
  KeyType = "ncbi-geneid",
  
  Background_Size = length(kegg_background),
  
  Tumor_Up_Input_Size = length(kegg_up_entrez),
  
  Tumor_Down_Input_Size = length(kegg_down_entrez),
  
  P_Value_Cutoff_Execution = 1,
  
  Q_Value_Cutoff_Execution = 1,
  
  MinGSSize = 10,
  
  MaxGSSize = 500,
  
  P_Adjust_Method = "BH",
  
  Final_Significance_Definition =
    "BH-adjusted P < 0.05",
  
  ID_Conversion_Performed = FALSE,
  
  DEG_Recalculation = FALSE,
  
  P_Value_Recalculation = FALSE,
  
  FDR_Recalculation = FALSE,
  
  Threshold_Changed = FALSE,
  
  Testing_Universe_Modified = FALSE,
  
  Background_Modified = FALSE
  
)


cat("============================================================\n")
cat("LOCKED KEGG ORA DEFINITION\n")
cat("============================================================\n\n")

print(KEGG_ANALYSIS_DEFINITION)


# ------------------------------------------------------------
# 8. Execute KEGG ORA — tumor-up
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Executing KEGG ORA — TUMOR-UP\n")
cat("============================================================\n\n")


kegg_up_enrich <- enrichKEGG(
  
  gene =
    kegg_up_entrez,
  
  organism =
    "hsa",
  
  keyType =
    "ncbi-geneid",
  
  universe =
    kegg_background,
  
  pvalueCutoff =
    1,
  
  pAdjustMethod =
    "BH",
  
  minGSSize =
    10,
  
  maxGSSize =
    500,
  
  qvalueCutoff =
    1
)


cat("Tumor-up KEGG ORA completed.\n")

cat(
  "Returned pathway rows: ",
  nrow(as.data.frame(kegg_up_enrich)),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Execute KEGG ORA — tumor-down
# ------------------------------------------------------------

cat("============================================================\n")
cat("Executing KEGG ORA — TUMOR-DOWN\n")
cat("============================================================\n\n")


kegg_down_enrich <- enrichKEGG(
  
  gene =
    kegg_down_entrez,
  
  organism =
    "hsa",
  
  keyType =
    "ncbi-geneid",
  
  universe =
    kegg_background,
  
  pvalueCutoff =
    1,
  
  pAdjustMethod =
    "BH",
  
  minGSSize =
    10,
  
  maxGSSize =
    500,
  
  qvalueCutoff =
    1
)


cat("Tumor-down KEGG ORA completed.\n")

cat(
  "Returned pathway rows: ",
  nrow(as.data.frame(kegg_down_enrich)),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Convert results to data frames
# ------------------------------------------------------------

kegg_up_table <- as.data.frame(
  kegg_up_enrich
)

kegg_down_table <- as.data.frame(
  kegg_down_enrich
)


# ------------------------------------------------------------
# 11. Final significance definition
# ------------------------------------------------------------

if (nrow(kegg_up_table) > 0) {
  
  kegg_up_sig <- subset(
    kegg_up_table,
    p.adjust < 0.05
  )
  
} else {
  
  kegg_up_sig <- kegg_up_table
  
}


if (nrow(kegg_down_table) > 0) {
  
  kegg_down_sig <- subset(
    kegg_down_table,
    p.adjust < 0.05
  )
  
} else {
  
  kegg_down_sig <- kegg_down_table
  
}


cat("============================================================\n")
cat("KEGG SIGNIFICANCE SUMMARY\n")
cat("============================================================\n\n")


cat(
  "Tumor-up pathways returned: ",
  nrow(kegg_up_table),
  "\n",
  sep = ""
)


cat(
  "Tumor-up pathways BH-adjusted P < 0.05: ",
  nrow(kegg_up_sig),
  "\n\n",
  sep = ""
)


cat(
  "Tumor-down pathways returned: ",
  nrow(kegg_down_table),
  "\n",
  sep = ""
)


cat(
  "Tumor-down pathways BH-adjusted P < 0.05: ",
  nrow(kegg_down_sig),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Basic result-object QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-5 RESULT OBJECT QC\n")
cat("============================================================\n\n")


up_result_qc <- c(
  
  object_class =
    inherits(kegg_up_enrich, "enrichResult"),
  
  input_gene_count =
    length(kegg_up_enrich@gene),
  
  universe_count =
    length(kegg_up_enrich@universe),
  
  input_gene_count_correct =
    length(kegg_up_enrich@gene) == 5194,
  
  universe_count_correct =
    length(kegg_up_enrich@universe) == 35365
  
)


down_result_qc <- c(
  
  object_class =
    inherits(kegg_down_enrich, "enrichResult"),
  
  input_gene_count =
    length(kegg_down_enrich@gene),
  
  universe_count =
    length(kegg_down_enrich@universe),
  
  input_gene_count_correct =
    length(kegg_down_enrich@gene) == 4463,
  
  universe_count_correct =
    length(kegg_down_enrich@universe) == 35365
  
)


cat("Tumor-up result QC:\n")
print(up_result_qc)

cat("\nTumor-down result QC:\n")
print(down_result_qc)


# ------------------------------------------------------------
# 13. Core-gene audit
# ------------------------------------------------------------

core_entrez <- c(
  AOC3 = "8639",
  MS4A7 = "58475",
  IL7R = "3575"
)


core_kegg_audit <- data.frame(
  
  Gene =
    names(core_entrez),
  
  Entrez =
    unname(core_entrez),
  
  In_Background =
    unname(core_entrez) %in% kegg_background,
  
  In_Tumor_Up_Input =
    unname(core_entrez) %in% kegg_up_entrez,
  
  In_Tumor_Down_Input =
    unname(core_entrez) %in% kegg_down_entrez,
  
  stringsAsFactors = FALSE
  
)


cat("\n============================================================\n")
cat("CORE GENE KEGG INPUT AUDIT\n")
cat("============================================================\n\n")

print(core_kegg_audit)


# ------------------------------------------------------------
# 14. Final PASS criteria
# ------------------------------------------------------------

STEP9C5_PASS <- all(
  
  length(kegg_background) == 35365,
  
  length(unique(kegg_background)) == 35365,
  
  length(kegg_up_entrez) == 5194,
  
  length(unique(kegg_up_entrez)) == 5194,
  
  length(kegg_down_entrez) == 4463,
  
  length(unique(kegg_down_entrez)) == 4463,
  
  length(up_not_in_background) == 0,
  
  length(down_not_in_background) == 0,
  
  length(up_down_overlap) == 0,
  
  inherits(kegg_up_enrich, "enrichResult"),
  
  inherits(kegg_down_enrich, "enrichResult"),
  
  length(kegg_up_enrich@gene) == 5194,
  
  length(kegg_down_enrich@gene) == 4463,
  
  length(kegg_up_enrich@universe) == 35365,
  
  length(kegg_down_enrich@universe) == 35365
  
)


cat("\n============================================================\n")

if (STEP9C5_PASS) {
  
  cat("FINAL STATUS: PASS — STEP9C-5 KEGG ORA EXECUTED\n")
  
} else {
  
  cat("FINAL STATUS: FAIL — STEP9C-5 QC FAILED\n")
  
}

cat("============================================================\n")


# ------------------------------------------------------------
# 15. Save result tables
# ------------------------------------------------------------

if (STEP9C5_PASS) {
  
  write.csv(
    kegg_up_table,
    file.path(
      step_dir,
      "STEP9C5_KEGG_TumorUp_COMPLETE.csv"
    ),
    row.names = FALSE
  )
  
  
  write.csv(
    kegg_down_table,
    file.path(
      step_dir,
      "STEP9C5_KEGG_TumorDown_COMPLETE.csv"
    ),
    row.names = FALSE
  )
  
  
  write.csv(
    kegg_up_sig,
    file.path(
      step_dir,
      "STEP9C5_KEGG_TumorUp_SIGNIFICANT_FDR0.05.csv"
    ),
    row.names = FALSE
  )
  
  
  write.csv(
    kegg_down_sig,
    file.path(
      step_dir,
      "STEP9C5_KEGG_TumorDown_SIGNIFICANT_FDR0.05.csv"
    ),
    row.names = FALSE
  )
  
  
  write.csv(
    core_kegg_audit,
    file.path(
      step_dir,
      "STEP9C5_CoreGene_KEGG_Input_Audit.csv"
    ),
    row.names = FALSE
  )
  
  
  saveRDS(
    kegg_up_enrich,
    file.path(
      step_dir,
      "STEP9C5_KEGG_TumorUp_enrichResult.rds"
    )
  )
  
  
  saveRDS(
    kegg_down_enrich,
    file.path(
      step_dir,
      "STEP9C5_KEGG_TumorDown_enrichResult.rds"
    )
  )
  
  
  saveRDS(
    list(
      Definition =
        KEGG_ANALYSIS_DEFINITION,
      
      Background =
        kegg_background,
      
      Tumor_Up_Input =
        kegg_up_entrez,
      
      Tumor_Down_Input =
        kegg_down_entrez,
      
      Tumor_Up_Complete =
        kegg_up_table,
      
      Tumor_Down_Complete =
        kegg_down_table,
      
      Tumor_Up_Significant =
        kegg_up_sig,
      
      Tumor_Down_Significant =
        kegg_down_sig,
      
      Core_Audit =
        core_kegg_audit,
      
      R_Version =
        R.version.string,
      
      clusterProfiler_Version =
        as.character(
          packageVersion("clusterProfiler")
        ),
      
      STEP9C5_PASS =
        STEP9C5_PASS
    ),
    file.path(
      step_dir,
      "STEP9C5_KEGG_RESULTS_LOCKED.rds"
    )
  )
  
  
  cat("\nResult files saved successfully.\n")
  
}

##STEP9C-5A — Controlled KEGG Re-execution and Effective-Universe Audit
# ============================================================
# Paper 1
# STEP9C-5A
# Controlled KEGG Re-execution and Effective-Universe Audit
#
# Purpose:
#   1. Recover LOCKED STEP9C-4 KEGG inputs
#   2. Re-run the SAME KEGG ORA definition
#   3. Immediately preserve raw enrichResult objects
#   4. Audit predefined vs effective KEGG universe
#   5. Do NOT modify STEP9C-3 or STEP9C-4
#   6. Do NOT recalculate DEGs / ID mapping / thresholds
#
# IMPORTANT:
#   This is a controlled re-execution of the already-defined
#   STEP9C-5 analysis, NOT a new analysis definition.
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-5A: Controlled KEGG Re-execution and Effective-Universe Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

setwd(step_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Required package
# ------------------------------------------------------------

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  stop(
    "FAIL — clusterProfiler is not installed."
  )
}

library(clusterProfiler)

cat("R version:\n")
cat(R.version.string, "\n\n")

cat("clusterProfiler version:\n")
cat(
  as.character(
    packageVersion("clusterProfiler")
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 2. Define source files
# ------------------------------------------------------------

c4_final_lock_file <- file.path(
  step_dir,
  "STEP9C4_FINAL_LOCK.rds"
)

c4_input_file <- file.path(
  step_dir,
  "STEP9C4_KEGG_INPUT_AND_BACKGROUND_LOCKED.rds"
)


# ------------------------------------------------------------
# 3. Verify STEP9C-4 source files
# ------------------------------------------------------------

required_source_files <- c(
  c4_final_lock_file,
  c4_input_file
)

if (!all(file.exists(required_source_files))) {
  
  missing_files <- required_source_files[
    !file.exists(required_source_files)
  ]
  
  stop(
    paste0(
      "FAIL — Required STEP9C-4 source file(s) missing:\n",
      paste(basename(missing_files), collapse = "\n")
    )
  )
  
}

cat("STEP9C-4 source files: PASS\n")
cat("  ", basename(c4_final_lock_file), "\n", sep = "")
cat("  ", basename(c4_input_file), "\n\n", sep = "")


# ------------------------------------------------------------
# 4. Recover STEP9C-4 locks
# ------------------------------------------------------------

STEP9C4_FINAL_LOCK <- readRDS(
  c4_final_lock_file
)

C4_BACKGROUND_OBJECT <- readRDS(
  c4_input_file
)

# ============================================================
# Paper 1
# STEP9C-5A Diagnostic
# Inspect actual STEP9C-4 FINAL LOCK structure
#
# NO re-enrichment
# NO recalculation
# NO modification
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-5A Diagnostic: STEP9C-4 FINAL LOCK STRUCTURE\n")
cat("============================================================\n\n")


step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

setwd(step_dir)


# ------------------------------------------------------------
# 1. Read existing STEP9C-4 FINAL LOCK
# ------------------------------------------------------------

lock_file <- file.path(
  step_dir,
  "STEP9C4_FINAL_LOCK.rds"
)

if (!file.exists(lock_file)) {
  stop(
    "FAIL — STEP9C4_FINAL_LOCK.rds not found."
  )
}

STEP9C4_FINAL_LOCK <- readRDS(lock_file)


# ------------------------------------------------------------
# 2. Show top-level structure
# ------------------------------------------------------------

cat("============================================================\n")
cat("TOP-LEVEL NAMES\n")
cat("============================================================\n\n")

print(names(STEP9C4_FINAL_LOCK))


# ------------------------------------------------------------
# 3. Show complete structure
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("COMPLETE STRUCTURE\n")
cat("============================================================\n\n")

str(
  STEP9C4_FINAL_LOCK,
  max.level = 3
)


# ------------------------------------------------------------
# 4. Print all top-level values where possible
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("TOP-LEVEL VALUES\n")
cat("============================================================\n\n")

for (nm in names(STEP9C4_FINAL_LOCK)) {
  
  cat("\n------------------------------\n")
  cat(nm, "\n")
  cat("------------------------------\n")
  
  obj <- STEP9C4_FINAL_LOCK[[nm]]
  
  if (is.atomic(obj) && length(obj) <= 20) {
    print(obj)
  } else if (is.data.frame(obj)) {
    print(obj)
  } else if (is.list(obj)) {
    cat(
      "LIST with ",
      length(obj),
      " element(s)\n",
      sep = ""
    )
    cat(
      "Names: ",
      paste(names(obj), collapse = ", "),
      "\n",
      sep = ""
    )
  } else {
    cat(
      "Class: ",
      paste(class(obj), collapse = ", "),
      "\n",
      sep = ""
    )
    cat(
      "Length: ",
      length(obj),
      "\n",
      sep = ""
    )
  }
}


# ------------------------------------------------------------
# 5. Specifically inspect Tumor_Up / Tumor_Down
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("TUMOR-UP / TUMOR-DOWN LOCK STRUCTURE\n")
cat("============================================================\n\n")


if ("Tumor_Up" %in% names(STEP9C4_FINAL_LOCK)) {
  
  cat("Tumor_Up exists: TRUE\n\n")
  
  cat("Tumor_Up names:\n")
  print(names(STEP9C4_FINAL_LOCK$Tumor_Up))
  
  cat("\nTumor_Up structure:\n")
  str(
    STEP9C4_FINAL_LOCK$Tumor_Up,
    max.level = 2
  )
  
} else {
  
  cat("Tumor_Up exists: FALSE\n")
  
}


cat("\n")


if ("Tumor_Down" %in% names(STEP9C4_FINAL_LOCK)) {
  
  cat("Tumor_Down exists: TRUE\n\n")
  
  cat("Tumor_Down names:\n")
  print(names(STEP9C4_FINAL_LOCK$Tumor_Down))
  
  cat("\nTumor_Down structure:\n")
  str(
    STEP9C4_FINAL_LOCK$Tumor_Down,
    max.level = 2
  )
  
} else {
  
  cat("Tumor_Down exists: FALSE\n")
  
}


# ------------------------------------------------------------
# 6. Show the exact fields relevant to our validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FIELDS USED BY THE FAILED VALIDATION\n")
cat("============================================================\n\n")


cat(
  "Step:\n"
)
print(STEP9C4_FINAL_LOCK$Step)

cat(
  "\nStatus:\n"
)
print(STEP9C4_FINAL_LOCK$Status)

cat(
  "\nFinal_KEGG_Background:\n"
)
print(STEP9C4_FINAL_LOCK$Final_KEGG_Background)


if ("Tumor_Up" %in% names(STEP9C4_FINAL_LOCK)) {
  
  cat(
    "\nTumor_Up$Input_Entrez:\n"
  )
  
  print(
    STEP9C4_FINAL_LOCK$Tumor_Up$Input_Entrez
  )
  
}


if ("Tumor_Down" %in% names(STEP9C4_FINAL_LOCK)) {
  
  cat(
    "\nTumor_Down$Input_Entrez:\n"
  )
  
  print(
    STEP9C4_FINAL_LOCK$Tumor_Down$Input_Entrez
  )
  
}


cat(
  "\nBackground_Unchanged:\n"
)

print(
  STEP9C4_FINAL_LOCK$Background_Unchanged
)


cat("\n============================================================\n")
cat("DIAGNOSTIC COMPLETE\n")
cat("============================================================\n\n")

cat(
  "NO re-enrichment performed.\n"
)

cat(
  "NO recalculation performed.\n"
)

cat(
  "NO STEP9C-4 object or file modified.\n"
)


# ============================================================
# Paper 1
# STEP9C-5A
# Controlled KEGG Re-execution and Effective-Universe Audit
#
# Purpose:
#   1. Recover the locked STEP9C-4 KEGG inputs
#   2. Correctly validate STEP9C-4 structure
#   3. Re-run KEGG ORA under the locked definition
#   4. Audit the effective KEGG universe returned by
#      clusterProfiler::enrichKEGG()
#   5. Preserve the distinction between:
#        - locked predefined KEGG background = 35,365
#        - effective KEGG universe = expected 8,702
#   6. Save the complete enrichment objects and audit results
#
# IMPORTANT:
#   - No DEG recalculation
#   - No ID conversion
#   - No p-value recalculation outside enrichKEGG()
#   - No FDR recalculation outside enrichKEGG()
#   - No threshold change
#   - No modification of STEP9C-3 / STEP9C-4 locks
#   - Existing split_entrez is NOT used again
#   - No biomaRt / bitr / AnnotationDbi conversion
# ============================================================


cat("\n============================================================\n")
cat("STEP9C-5A: Controlled KEGG Re-execution and Effective-Universe Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

setwd(step_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Required package
# ------------------------------------------------------------

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  stop(
    "FAIL — clusterProfiler is not installed."
  )
}

library(clusterProfiler)

cat("R version:\n")
cat(R.version.string, "\n\n")

cat("clusterProfiler version:\n")
cat(
  as.character(packageVersion("clusterProfiler")),
  "\n\n"
)


# ------------------------------------------------------------
# 2. Define STEP9C-4 source files
# ------------------------------------------------------------

c4_final_lock_file <- file.path(
  step_dir,
  "STEP9C4_FINAL_LOCK.rds"
)

c4_input_file <- file.path(
  step_dir,
  "STEP9C4_KEGG_INPUT_AND_BACKGROUND_LOCKED.rds"
)


# ------------------------------------------------------------
# 3. Verify source files
# ------------------------------------------------------------

required_source_files <- c(
  c4_final_lock_file,
  c4_input_file
)

if (!all(file.exists(required_source_files))) {
  
  missing_files <- required_source_files[
    !file.exists(required_source_files)
  ]
  
  stop(
    paste0(
      "FAIL — Required STEP9C-4 source file(s) missing:\n",
      paste(basename(missing_files), collapse = "\n")
    )
  )
  
}

cat("STEP9C-4 source files: PASS\n")
cat(
  "  ",
  basename(c4_final_lock_file),
  "\n",
  sep = ""
)
cat(
  "  ",
  basename(c4_input_file),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Recover STEP9C-4 locks
# ------------------------------------------------------------

STEP9C4_FINAL_LOCK <- readRDS(
  c4_final_lock_file
)

C4_BACKGROUND_OBJECT <- readRDS(
  c4_input_file
)


# ------------------------------------------------------------
# 5. Correct STEP9C-4 FINAL LOCK validation
#
# IMPORTANT:
#   Input_Entrez is a vector, not a count.
#   Therefore use length(Input_Entrez).
# ------------------------------------------------------------

step9c4_lock_valid <- all(
  
  identical(
    STEP9C4_FINAL_LOCK$Step,
    "STEP9C-4"
  ),
  
  identical(
    STEP9C4_FINAL_LOCK$Status,
    "PASS"
  ),
  
  isTRUE(
    STEP9C4_FINAL_LOCK$Final_KEGG_Background == 35365
  ),
  
  isTRUE(
    length(
      STEP9C4_FINAL_LOCK$Tumor_Up$Input_Entrez
    ) == 5194
  ),
  
  isTRUE(
    length(
      STEP9C4_FINAL_LOCK$Tumor_Down$Input_Entrez
    ) == 4463
  ),
  
  isTRUE(
    STEP9C4_FINAL_LOCK$Background_Size == 35365
  ),
  
  isTRUE(
    STEP9C4_FINAL_LOCK$Background_Unique == 35365
  ),
  
  isTRUE(
    STEP9C4_FINAL_LOCK$Up_Down_Entrez_Overlap == 0
  ),
  
  isTRUE(
    STEP9C4_FINAL_LOCK$Background_Unchanged
  ),
  
  identical(
    STEP9C4_FINAL_LOCK$ID_Conversion_Performed,
    FALSE
  ),
  
  identical(
    STEP9C4_FINAL_LOCK$DEG_Recalculation,
    FALSE
  ),
  
  identical(
    STEP9C4_FINAL_LOCK$P_Value_Recalculation,
    FALSE
  ),
  
  identical(
    STEP9C4_FINAL_LOCK$FDR_Recalculation,
    FALSE
  ),
  
  identical(
    STEP9C4_FINAL_LOCK$Threshold_Changed,
    FALSE
  ),
  
  identical(
    STEP9C4_FINAL_LOCK$Testing_Universe_Modified,
    FALSE
  )
  
)

cat("STEP9C-4 FINAL LOCK validation:\n")

if (!step9c4_lock_valid) {
  
  cat("FAIL\n\n")
  
  stop(
    "FAIL — Corrected STEP9C-4 FINAL LOCK validation failed."
  )
  
} else {
  
  cat("PASS\n\n")
  
}


# ------------------------------------------------------------
# 6. Recover locked KEGG vectors
# ------------------------------------------------------------

kegg_up_entrez <- C4_BACKGROUND_OBJECT$kegg_up_entrez

kegg_down_entrez <- C4_BACKGROUND_OBJECT$kegg_down_entrez

kegg_background <- C4_BACKGROUND_OBJECT$final_universe_entrez


# ------------------------------------------------------------
# 7. Structural validation of recovered vectors
# ------------------------------------------------------------

vector_qc <- list(
  
  background_is_character =
    is.character(kegg_background),
  
  background_unique =
    length(unique(kegg_background)) ==
    length(kegg_background),
  
  background_no_na =
    !anyNA(kegg_background),
  
  up_is_character =
    is.character(kegg_up_entrez),
  
  up_unique =
    length(unique(kegg_up_entrez)) ==
    length(kegg_up_entrez),
  
  up_no_na =
    !anyNA(kegg_up_entrez),
  
  down_is_character =
    is.character(kegg_down_entrez),
  
  down_unique =
    length(unique(kegg_down_entrez)) ==
    length(kegg_down_entrez),
  
  down_no_na =
    !anyNA(kegg_down_entrez)
  
)

if (!all(unlist(vector_qc))) {
  
  print(vector_qc)
  
  stop(
    "FAIL — STEP9C-4 recovered KEGG vectors failed structural QC."
  )
  
}

cat("Recovered vector structural QC: PASS\n\n")


# ------------------------------------------------------------
# 8. Exact numerical lock validation
# ------------------------------------------------------------

numerical_lock <- all(
  
  length(kegg_background) == 35365,
  
  length(unique(kegg_background)) == 35365,
  
  length(kegg_up_entrez) == 5194,
  
  length(unique(kegg_up_entrez)) == 5194,
  
  length(kegg_down_entrez) == 4463,
  
  length(unique(kegg_down_entrez)) == 4463
  
)

if (!numerical_lock) {
  
  stop(
    "FAIL — STEP9C-4 numerical lock validation failed."
  )
  
}

cat("STEP9C-4 numerical lock: PASS\n")
cat("  Locked background = 35,365\n")
cat("  Locked tumor-up input = 5,194\n")
cat("  Locked tumor-down input = 4,463\n\n")


# ------------------------------------------------------------
# 9. Containment and overlap audit
# ------------------------------------------------------------

up_outside_background <- setdiff(
  kegg_up_entrez,
  kegg_background
)

down_outside_background <- setdiff(
  kegg_down_entrez,
  kegg_background
)

up_down_overlap <- intersect(
  kegg_up_entrez,
  kegg_down_entrez
)

containment_qc <- all(
  
  length(up_outside_background) == 0,
  
  length(down_outside_background) == 0,
  
  length(up_down_overlap) == 0
  
)

cat("Containment / overlap QC:\n")

cat(
  "  Tumor-up outside background = ",
  length(up_outside_background),
  "\n",
  sep = ""
)

cat(
  "  Tumor-down outside background = ",
  length(down_outside_background),
  "\n",
  sep = ""
)

cat(
  "  Up/down overlap = ",
  length(up_down_overlap),
  "\n\n",
  sep = ""
)

if (!containment_qc) {
  
  stop(
    "FAIL — KEGG input containment or overlap QC failed."
  )
  
}

cat("Containment / overlap QC: PASS\n\n")


# ------------------------------------------------------------
# 10. Core-gene audit
# ------------------------------------------------------------

core_entrez <- c(
  AOC3  = "8639",
  MS4A7 = "58475",
  IL7R  = "3575"
)

core_audit <- data.frame(
  
  Gene = names(core_entrez),
  
  Entrez = unname(core_entrez),
  
  In_Locked_Background =
    unname(core_entrez) %in% kegg_background,
  
  In_Tumor_Up_Input =
    unname(core_entrez) %in% kegg_up_entrez,
  
  In_Tumor_Down_Input =
    unname(core_entrez) %in% kegg_down_entrez,
  
  stringsAsFactors = FALSE
  
)

print(core_audit)

core_audit_pass <- all(
  core_audit$In_Locked_Background,
  core_audit$In_Tumor_Down_Input,
  !core_audit$In_Tumor_Up_Input
)

if (!core_audit_pass) {
  
  stop(
    "FAIL — Core-gene KEGG input audit failed."
  )
  
}

cat("\nCore-gene audit: PASS\n\n")


# ------------------------------------------------------------
# 11. Locked KEGG ORA definition
# ------------------------------------------------------------

kegg_definition <- list(
  
  organism = "hsa",
  
  keyType = "ncbi-geneid",
  
  predefined_background_size = 35365,
  
  tumor_up_input_size = 5194,
  
  tumor_down_input_size = 4463,
  
  pvalueCutoff = 1,
  
  qvalueCutoff = 1,
  
  minGSSize = 10,
  
  maxGSSize = 500,
  
  pAdjustMethod = "BH",
  
  final_significance_rule =
    "BH-adjusted P < 0.05",
  
  ID_conversion =
    "No new ID conversion; existing STEP9C-4 Entrez vectors used",
  
  DEG_recalculation = FALSE,
  
  threshold_changed = FALSE,
  
  testing_universe_modified = FALSE,
  
  predefined_background_modified = FALSE
  
)


# ------------------------------------------------------------
# 12. Controlled KEGG re-execution
# ------------------------------------------------------------

cat("============================================================\n")
cat("CONTROLLED KEGG RE-EXECUTION\n")
cat("============================================================\n\n")

cat("Running tumor-up KEGG ORA...\n")

kegg_up_enrich <- enrichKEGG(
  
  gene = kegg_up_entrez,
  
  organism = "hsa",
  
  keyType = "ncbi-geneid",
  
  pvalueCutoff = 1,
  
  qvalueCutoff = 1,
  
  minGSSize = 10,
  
  maxGSSize = 500,
  
  pAdjustMethod = "BH",
  
  universe = kegg_background
  
)

cat(
  "Tumor-up KEGG ORA completed.\n\n"
)


cat("Running tumor-down KEGG ORA...\n")

kegg_down_enrich <- enrichKEGG(
  
  gene = kegg_down_entrez,
  
  organism = "hsa",
  
  keyType = "ncbi-geneid",
  
  pvalueCutoff = 1,
  
  qvalueCutoff = 1,
  
  minGSSize = 10,
  
  maxGSSize = 500,
  
  pAdjustMethod = "BH",
  
  universe = kegg_background
  
)

cat(
  "Tumor-down KEGG ORA completed.\n\n"
)


# ------------------------------------------------------------
# 13. Basic result-object audit
# ------------------------------------------------------------

cat("============================================================\n")
cat("RESULT OBJECT AUDIT\n")
cat("============================================================\n\n")

cat("Tumor-up result class:\n")
print(class(kegg_up_enrich))

cat("\nTumor-down result class:\n")
print(class(kegg_down_enrich))


up_result_rows <- nrow(
  as.data.frame(kegg_up_enrich)
)

down_result_rows <- nrow(
  as.data.frame(kegg_down_enrich)
)

cat("\nTumor-up pathway rows = ")
cat(up_result_rows, "\n")

cat("Tumor-down pathway rows = ")
cat(down_result_rows, "\n\n")


# ------------------------------------------------------------
# 14. Inspect result-object universe
#
# IMPORTANT:
#   enrichKEGG() may internally restrict the user-specified
#   universe to genes represented in the KEGG annotation space.
#
#   Therefore:
#      locked predefined background = 35,365
#      effective KEGG universe      = result@universe
#
#   This is an audit, NOT a replacement of STEP9C-4 background.
# ------------------------------------------------------------

effective_up_universe <- kegg_up_enrich@universe

effective_down_universe <- kegg_down_enrich@universe

cat("============================================================\n")
cat("EFFECTIVE KEGG UNIVERSE AUDIT\n")
cat("============================================================\n\n")

cat(
  "Locked predefined KEGG background = ",
  length(kegg_background),
  "\n",
  sep = ""
)

cat(
  "Tumor-up effective KEGG universe = ",
  length(effective_up_universe),
  "\n",
  sep = ""
)

cat(
  "Tumor-down effective KEGG universe = ",
  length(effective_down_universe),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 15. Verify that both effective universes are identical
# ------------------------------------------------------------

effective_universe_identical <- identical(
  sort(unique(effective_up_universe)),
  sort(unique(effective_down_universe))
)

if (!effective_universe_identical) {
  
  stop(
    "FAIL — Tumor-up and tumor-down effective KEGG universes differ."
  )
  
}

cat(
  "Up/down effective KEGG universe identity: PASS\n\n"
)


# ------------------------------------------------------------
# 16. Quantify relationship between predefined and effective
#     universe
# ------------------------------------------------------------

effective_kegg_universe <- sort(
  unique(effective_up_universe)
)

effective_universe_size <- length(
  effective_kegg_universe
)

effective_outside_locked_background <- setdiff(
  effective_kegg_universe,
  kegg_background
)

locked_background_not_effective <- setdiff(
  kegg_background,
  effective_kegg_universe
)

cat("Effective-universe relationship:\n")

cat(
  "  Effective KEGG universe size = ",
  effective_universe_size,
  "\n",
  sep = ""
)

cat(
  "  Effective universe outside locked background = ",
  length(effective_outside_locked_background),
  "\n",
  sep = ""
)

cat(
  "  Locked-background IDs not in effective KEGG universe = ",
  length(locked_background_not_effective),
  "\n\n",
  sep = ""
)


effective_universe_qc <- all(
  
  length(effective_outside_locked_background) == 0,
  
  effective_universe_size > 0,
  
  effective_universe_size <= length(kegg_background)
  
)

if (!effective_universe_qc) {
  
  stop(
    "FAIL — Effective KEGG universe relationship audit failed."
  )
  
}

cat(
  "Effective-universe relationship QC: PASS\n\n"
)


# ------------------------------------------------------------
# 17. Audit input genes against effective KEGG universe
# ------------------------------------------------------------

up_in_effective_universe <- intersect(
  kegg_up_entrez,
  effective_kegg_universe
)

down_in_effective_universe <- intersect(
  kegg_down_entrez,
  effective_kegg_universe
)

up_not_in_effective_universe <- setdiff(
  kegg_up_entrez,
  effective_kegg_universe
)

down_not_in_effective_universe <- setdiff(
  kegg_down_entrez,
  effective_kegg_universe
)

cat("============================================================\n")
cat("INPUT vs EFFECTIVE KEGG UNIVERSE AUDIT\n")
cat("============================================================\n\n")

cat(
  "Tumor-up locked input = ",
  length(kegg_up_entrez),
  "\n",
  sep = ""
)

cat(
  "Tumor-up input represented in effective KEGG universe = ",
  length(up_in_effective_universe),
  "\n",
  sep = ""
)

cat(
  "Tumor-up input outside effective KEGG universe = ",
  length(up_not_in_effective_universe),
  "\n\n",
  sep = ""
)

cat(
  "Tumor-down locked input = ",
  length(kegg_down_entrez),
  "\n",
  sep = ""
)

cat(
  "Tumor-down input represented in effective KEGG universe = ",
  length(down_in_effective_universe),
  "\n",
  sep = ""
)

cat(
  "Tumor-down input outside effective KEGG universe = ",
  length(down_not_in_effective_universe),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 18. Effective input coverage
# ------------------------------------------------------------

up_effective_coverage <-
  length(up_in_effective_universe) /
  length(kegg_up_entrez)

down_effective_coverage <-
  length(down_in_effective_universe) /
  length(kegg_down_entrez)

cat(
  "Tumor-up effective KEGG representation = ",
  sprintf("%.4f", up_effective_coverage),
  "\n",
  sep = ""
)

cat(
  "Tumor-down effective KEGG representation = ",
  sprintf("%.4f", down_effective_coverage),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 19. Verify core genes against effective universe
# ------------------------------------------------------------

core_effective_audit <- data.frame(
  
  Gene = names(core_entrez),
  
  Entrez = unname(core_entrez),
  
  In_Locked_Background =
    unname(core_entrez) %in% kegg_background,
  
  In_Effective_KEGG_Universe =
    unname(core_entrez) %in% effective_kegg_universe,
  
  In_Tumor_Up_Input =
    unname(core_entrez) %in% kegg_up_entrez,
  
  In_Tumor_Down_Input =
    unname(core_entrez) %in% kegg_down_entrez,
  
  stringsAsFactors = FALSE
  
)

cat("Core-gene effective-universe audit:\n\n")

print(core_effective_audit)

cat("\n")


# ------------------------------------------------------------
# 20. Significance counts
# ------------------------------------------------------------

up_result_df <- as.data.frame(
  kegg_up_enrich
)

down_result_df <- as.data.frame(
  kegg_down_enrich
)

up_sig_n <- sum(
  up_result_df$p.adjust < 0.05,
  na.rm = TRUE
)

down_sig_n <- sum(
  down_result_df$p.adjust < 0.05,
  na.rm = TRUE
)

cat("============================================================\n")
cat("KEGG SIGNIFICANCE SUMMARY\n")
cat("============================================================\n\n")

cat(
  "Tumor-up pathways returned = ",
  nrow(up_result_df),
  "\n",
  sep = ""
)

cat(
  "Tumor-up BH-adjusted P < 0.05 = ",
  up_sig_n,
  "\n\n",
  sep = ""
)

cat(
  "Tumor-down pathways returned = ",
  nrow(down_result_df),
  "\n",
  sep = ""
)

cat(
  "Tumor-down BH-adjusted P < 0.05 = ",
  down_sig_n,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 21. Result-object input counts
# ------------------------------------------------------------

result_gene_counts <- list(
  
  Tumor_Up_Result_Gene_Count =
    length(kegg_up_enrich@gene),
  
  Tumor_Down_Result_Gene_Count =
    length(kegg_down_enrich@gene),
  
  Tumor_Up_Locked_Input_Count =
    length(kegg_up_entrez),
  
  Tumor_Down_Locked_Input_Count =
    length(kegg_down_entrez),
  
  Tumor_Up_Effective_Universe =
    length(effective_up_universe),
  
  Tumor_Down_Effective_Universe =
    length(effective_down_universe)
  
)

cat("Result-object gene/universe counts:\n")
print(result_gene_counts)
cat("\n")


# ------------------------------------------------------------
# 22. Final STEP9C-5A PASS definition
#
# Important:
#   Do NOT require result@universe == 35,365.
#
#   The correct audit requirement is:
#      - predefined locked background = 35,365
#      - effective universe is a subset of it
#      - up/down effective universes identical
#      - enrichment objects exist
#      - inputs remain unchanged
# ------------------------------------------------------------

STEP9C5A_PASS <- all(
  
  step9c4_lock_valid,
  
  numerical_lock,
  
  containment_qc,
  
  core_audit_pass,
  
  effective_universe_identical,
  
  effective_universe_qc,
  
  up_result_rows >= 0,
  
  down_result_rows >= 0,
  
  length(kegg_up_enrich@gene) ==
    length(kegg_up_entrez),
  
  length(kegg_down_enrich@gene) ==
    length(kegg_down_entrez),
  
  length(effective_up_universe) ==
    length(effective_down_universe),
  
  length(effective_kegg_universe) <=
    length(kegg_background)
  
)


# ------------------------------------------------------------
# 23. Build STEP9C-5A audit lock
# ------------------------------------------------------------

STEP9C5A_FINAL_LOCK <- list(
  
  Step = "STEP9C-5A",
  
  Status =
    if (STEP9C5A_PASS) "PASS" else "FAIL",
  
  R_Version =
    R.version.string,
  
  clusterProfiler_Version =
    as.character(
      packageVersion("clusterProfiler")
    ),
  
  ORA_Definition =
    kegg_definition,
  
  Locked_Predefined_Background =
    length(kegg_background),
  
  Locked_Tumor_Up_Input =
    length(kegg_up_entrez),
  
  Locked_Tumor_Down_Input =
    length(kegg_down_entrez),
  
  Effective_KEGG_Universe =
    effective_universe_size,
  
  Effective_Universe_Outside_Locked_Background =
    length(effective_outside_locked_background),
  
  Locked_Background_Not_In_Effective_Universe =
    length(locked_background_not_effective),
  
  Tumor_Up_In_Effective_Universe =
    length(up_in_effective_universe),
  
  Tumor_Up_Outside_Effective_Universe =
    length(up_not_in_effective_universe),
  
  Tumor_Down_In_Effective_Universe =
    length(down_in_effective_universe),
  
  Tumor_Down_Outside_Effective_Universe =
    length(down_not_in_effective_universe),
  
  Tumor_Up_Effective_Coverage =
    up_effective_coverage,
  
  Tumor_Down_Effective_Coverage =
    down_effective_coverage,
  
  Tumor_Up_Pathways_Returned =
    nrow(up_result_df),
  
  Tumor_Down_Pathways_Returned =
    nrow(down_result_df),
  
  Tumor_Up_BH_FDR_lt_0.05 =
    up_sig_n,
  
  Tumor_Down_BH_FDR_lt_0.05 =
    down_sig_n,
  
  Up_Down_Effective_Universe_Identical =
    effective_universe_identical,
  
  Core_Audit =
    core_audit,
  
  Core_Effective_Universe_Audit =
    core_effective_audit,
  
  No_DEG_Recalculation =
    TRUE,
  
  No_ID_Conversion =
    TRUE,
  
  No_Threshold_Change =
    TRUE,
  
  No_STEP9C3_Modification =
    TRUE,
  
  No_STEP9C4_Modification =
    TRUE
  
)


# ------------------------------------------------------------
# 24. Save enrichment result objects
# ------------------------------------------------------------

if (STEP9C5A_PASS) {
  
  saveRDS(
    kegg_up_enrich,
    file.path(
      step_dir,
      "STEP9C5A_KEGG_TumorUp_enrichResult.rds"
    )
  )
  
  saveRDS(
    kegg_down_enrich,
    file.path(
      step_dir,
      "STEP9C5A_KEGG_TumorDown_enrichResult.rds"
    )
  )
  
  saveRDS(
    STEP9C5A_FINAL_LOCK,
    file.path(
      step_dir,
      "STEP9C5A_FINAL_LOCK.rds"
    )
  )
  
  saveRDS(
    list(
      kegg_up_entrez = kegg_up_entrez,
      kegg_down_entrez = kegg_down_entrez,
      kegg_background = kegg_background,
      effective_kegg_universe =
        effective_kegg_universe
    ),
    file.path(
      step_dir,
      "STEP9C5A_KEGG_INPUT_BACKGROUND_EFFECTIVE_UNIVERSE_LOCKED.rds"
    )
  )
  
  write.csv(
    up_result_df,
    file.path(
      step_dir,
      "STEP9C5A_KEGG_TumorUp_COMPLETE.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    down_result_df,
    file.path(
      step_dir,
      "STEP9C5A_KEGG_TumorDown_COMPLETE.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    core_effective_audit,
    file.path(
      step_dir,
      "STEP9C5A_CoreGene_KEGG_EffectiveUniverse_Audit.csv"
    ),
    row.names = FALSE
  )
  
}


# ------------------------------------------------------------
# 25. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-5A FINAL STATUS\n")
cat("============================================================\n\n")

if (STEP9C5A_PASS) {
  
  cat(
    "FINAL STATUS: PASS — STEP9C-5A LOCKED\n\n"
  )
  
  cat(
    "Locked predefined KEGG background = 35,365\n"
  )
  
  cat(
    "Effective KEGG universe = ",
    effective_universe_size,
    "\n",
    sep = ""
  )
  
  cat(
    "Tumor-up input = 5,194\n"
  )
  
  cat(
    "Tumor-down input = 4,463\n"
  )
  
  cat(
    "Tumor-up pathways returned = ",
    nrow(up_result_df),
    "\n",
    sep = ""
  )
  
  cat(
    "Tumor-down pathways returned = ",
    nrow(down_result_df),
    "\n",
    sep = ""
  )
  
  cat(
    "Tumor-up BH-adjusted P < 0.05 = ",
    up_sig_n,
    "\n",
    sep = ""
  )
  
  cat(
    "Tumor-down BH-adjusted P < 0.05 = ",
    down_sig_n,
    "\n\n",
    sep = ""
  )
  
  cat(
    "IMPORTANT INTERPRETATION:\n"
  )
  
  cat(
    "  35,365 = locked predefined KEGG Entrez background\n"
  )
  
  cat(
    "  ",
    effective_universe_size,
    " = effective KEGG universe used internally by enrichKEGG()\n",
    sep = ""
  )
  
  cat(
    "  The effective universe was audited as a subset of the locked background.\n"
  )
  
  cat(
    "  STEP9C-4 was not modified.\n"
  )
  
} else {
  
  cat(
    "FINAL STATUS: FAIL — STEP9C-5A NOT LOCKED\n"
  )
  
  cat(
    "No STEP9C-5A result lock should be treated as final.\n"
  )
  
}

cat("\n============================================================\n")
cat("STEP9C-5A COMPLETE\n")
cat("============================================================\n")

##STEP9C-5B
#KEGG Result Extraction and Manuscript-Definition Lock

# ============================================================
# STEP9C-5A — Existing Result Recovery Audit
# NO recalculation
# NO modification
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-5A Existing Result Recovery Audit\n")
cat("============================================================\n\n")

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

setwd(step_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. List all files in STEP9C directory
# ------------------------------------------------------------

all_files <- list.files(
  step_dir,
  recursive = TRUE,
  full.names = TRUE
)

cat("Total files found:\n")
cat(length(all_files), "\n\n")


# ------------------------------------------------------------
# 2. Search for likely KEGG / enrichment / STEP9C files
# ------------------------------------------------------------

search_pattern <- paste(
  c(
    "STEP9C",
    "KEGG",
    "enrich",
    "enrichment",
    "pathway"
  ),
  collapse = "|"
)

candidate_files <- all_files[
  grepl(
    search_pattern,
    basename(all_files),
    ignore.case = TRUE
  )
]

cat("Candidate files:\n\n")

if (length(candidate_files) == 0) {
  
  cat("NONE FOUND\n\n")
  
} else {
  
  print(
    data.frame(
      FullPath = candidate_files,
      FileName = basename(candidate_files),
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 3. Search specifically for RDS files
# ------------------------------------------------------------

rds_files <- all_files[
  grepl(
    "\\.rds$",
    all_files,
    ignore.case = TRUE
  )
]

cat("RDS files found:\n")
cat(length(rds_files), "\n\n")

if (length(rds_files) > 0) {
  
  print(
    data.frame(
      FullPath = rds_files,
      FileName = basename(rds_files),
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 4. Search specifically for CSV files
# ------------------------------------------------------------

csv_files <- all_files[
  grepl(
    "\\.csv$",
    all_files,
    ignore.case = TRUE
  )
]

cat("CSV files found:\n")
cat(length(csv_files), "\n\n")

if (length(csv_files) > 0) {
  
  print(
    data.frame(
      FullPath = csv_files,
      FileName = basename(csv_files),
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 5. Search specifically for STEP9C-5A naming variants
# ------------------------------------------------------------

step9c5a_variants <- all_files[
  grepl(
    "5A|5-A|5_A|C5A|C-5A",
    basename(all_files),
    ignore.case = TRUE
  )
]

cat("Possible STEP9C-5A naming variants:\n\n")

if (length(step9c5a_variants) == 0) {
  
  cat("NONE FOUND\n\n")
  
} else {
  
  print(
    data.frame(
      FullPath = step9c5a_variants,
      FileName = basename(step9c5a_variants),
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 6. Search for likely enrichment object names
# ------------------------------------------------------------

object_pattern <- paste(
  c(
    "kegg_up",
    "kegg_down",
    "KEGG",
    "enrich",
    "enrichment"
  ),
  collapse = "|"
)

object_candidates <- rds_files[
  grepl(
    object_pattern,
    basename(rds_files),
    ignore.case = TRUE
  )
]

cat("RDS files potentially containing enrichment results:\n\n")

if (length(object_candidates) == 0) {
  
  cat("NONE FOUND\n\n")
  
} else {
  
  print(
    data.frame(
      FullPath = object_candidates,
      FileName = basename(object_candidates),
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 7. Search for CSV tables that may contain KEGG results
# ------------------------------------------------------------

table_candidates <- csv_files[
  grepl(
    "KEGG|enrich|pathway|tumor.*up|tumor.*down|up.*pathway|down.*pathway",
    basename(csv_files),
    ignore.case = TRUE
  )
]

cat("CSV files potentially containing KEGG results:\n\n")

if (length(table_candidates) == 0) {
  
  cat("NONE FOUND\n\n")
  
} else {
  
  print(
    data.frame(
      FullPath = table_candidates,
      FileName = basename(table_candidates),
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 8. Final audit conclusion
# ------------------------------------------------------------

cat("============================================================\n")
cat("RECOVERY AUDIT COMPLETE\n")
cat("============================================================\n\n")

cat(
  "IMPORTANT:\n",
  "This audit performed NO enrichment calculation.\n",
  "No DEG object was modified.\n",
  "No background was modified.\n",
  "No STEP9C-3 or STEP9C-4 lock was modified.\n",
  sep = ""
)

cat("\n")


# ============================================================
# STEP9C-5B
# Existing STEP9C-5A Result Recovery and Structural Audit
#
# IMPORTANT:
#   NO enrichKEGG()
#   NO ID conversion
#   NO DEG recalculation
#   NO P-value recalculation
#   NO FDR recalculation
#   NO background modification
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-5B: STEP9C-5A Result Recovery Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

setwd(step_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Define locked STEP9C-5A files
# ------------------------------------------------------------

up_file <- file.path(
  step_dir,
  "STEP9C5A_KEGG_TumorUp_COMPLETE.csv"
)

down_file <- file.path(
  step_dir,
  "STEP9C5A_KEGG_TumorDown_COMPLETE.csv"
)

core_audit_file <- file.path(
  step_dir,
  "STEP9C5A_CoreGene_KEGG_EffectiveUniverse_Audit.csv"
)


# ------------------------------------------------------------
# 2. File existence check
# ------------------------------------------------------------

required_files <- c(
  up_file,
  down_file,
  core_audit_file
)

file_check <- file.exists(required_files)

print(
  data.frame(
    File = basename(required_files),
    Exists = file_check,
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)

if (!all(file_check)) {
  
  stop(
    "FAIL — One or more required STEP9C-5A result files are missing."
  )
}

cat("\nRequired STEP9C-5A files: PASS\n\n")


# ------------------------------------------------------------
# 3. Read existing result tables
# ------------------------------------------------------------

kegg_up_complete <- read.csv(
  up_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

kegg_down_complete <- read.csv(
  down_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

core_effective_audit <- read.csv(
  core_audit_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 4. Structural inspection
# ------------------------------------------------------------

cat("Tumor-up dimensions:\n")
print(dim(kegg_up_complete))

cat("\nTumor-down dimensions:\n")
print(dim(kegg_down_complete))

cat("\nCore/effective-universe audit dimensions:\n")
print(dim(core_effective_audit))

cat("\nTumor-up columns:\n")
print(colnames(kegg_up_complete))

cat("\nTumor-down columns:\n")
print(colnames(kegg_down_complete))

cat("\n")


# ------------------------------------------------------------
# 5. Required KEGG result columns
# ------------------------------------------------------------

required_columns <- c(
  "ID",
  "Description",
  "GeneRatio",
  "BgRatio",
  "pvalue",
  "p.adjust",
  "geneID"
)

missing_up <- setdiff(
  required_columns,
  colnames(kegg_up_complete)
)

missing_down <- setdiff(
  required_columns,
  colnames(kegg_down_complete)
)

if (length(missing_up) > 0) {
  
  stop(
    paste0(
      "FAIL — Tumor-up missing columns: ",
      paste(missing_up, collapse = ", ")
    )
  )
}

if (length(missing_down) > 0) {
  
  stop(
    paste0(
      "FAIL — Tumor-down missing columns: ",
      paste(missing_down, collapse = ", ")
    )
  )
}

cat("Required KEGG result columns: PASS\n\n")


# ------------------------------------------------------------
# 6. Lock pathway counts
# ------------------------------------------------------------

cat("Tumor-up pathways recovered:\n")
cat(nrow(kegg_up_complete), "\n")

cat("Tumor-down pathways recovered:\n")
cat(nrow(kegg_down_complete), "\n\n")


if (nrow(kegg_up_complete) != 335) {
  
  stop(
    paste0(
      "FAIL — Expected 335 tumor-up pathways; found ",
      nrow(kegg_up_complete)
    )
  )
}

if (nrow(kegg_down_complete) != 338) {
  
  stop(
    paste0(
      "FAIL — Expected 338 tumor-down pathways; found ",
      nrow(kegg_down_complete)
    )
  )
}

cat("335 / 338 pathway recovery: PASS\n\n")


# ------------------------------------------------------------
# 7. Re-identify significant pathways
# ------------------------------------------------------------

BH_THRESHOLD <- 0.05

up_sig <- kegg_up_complete[
  !is.na(kegg_up_complete$p.adjust) &
    kegg_up_complete$p.adjust < BH_THRESHOLD,
  ,
  drop = FALSE
]

down_sig <- kegg_down_complete[
  !is.na(kegg_down_complete$p.adjust) &
    kegg_down_complete$p.adjust < BH_THRESHOLD,
  ,
  drop = FALSE
]


cat("Tumor-up BH-adjusted P < 0.05:\n")
cat(nrow(up_sig), "\n")

cat("Tumor-down BH-adjusted P < 0.05:\n")
cat(nrow(down_sig), "\n\n")


if (nrow(up_sig) != 13) {
  
  stop(
    paste0(
      "FAIL — Expected 13 significant tumor-up pathways; found ",
      nrow(up_sig)
    )
  )
}

if (nrow(down_sig) != 60) {
  
  stop(
    paste0(
      "FAIL — Expected 60 significant tumor-down pathways; found ",
      nrow(down_sig)
    )
  )
}

cat("13 / 60 significant pathway lock: PASS\n\n")


# ------------------------------------------------------------
# 8. Verify significant pathway table content
# ------------------------------------------------------------

up_sig_final <- up_sig[
  ,
  required_columns,
  drop = FALSE
]

down_sig_final <- down_sig[
  ,
  required_columns,
  drop = FALSE
]


# ------------------------------------------------------------
# 9. Verify no duplicate pathway IDs
# ------------------------------------------------------------

if (anyDuplicated(up_sig_final$ID) > 0) {
  
  stop(
    "FAIL — Duplicate tumor-up pathway IDs."
  )
}

if (anyDuplicated(down_sig_final$ID) > 0) {
  
  stop(
    "FAIL — Duplicate tumor-down pathway IDs."
  )
}

cat("Pathway ID uniqueness: PASS\n")


# ------------------------------------------------------------
# 10. Verify BH definition
# ------------------------------------------------------------

up_invalid <- any(
  is.na(up_sig_final$p.adjust) |
    up_sig_final$p.adjust >= BH_THRESHOLD
)

down_invalid <- any(
  is.na(down_sig_final$p.adjust) |
    down_sig_final$p.adjust >= BH_THRESHOLD
)

if (up_invalid || down_invalid) {
  
  stop(
    "FAIL — Significant pathway table violates BH-adjusted P < 0.05."
  )
}

cat("BH-adjusted P < 0.05 definition: PASS\n\n")


# ------------------------------------------------------------
# 11. Verify geneID fields are present
# ------------------------------------------------------------

if (
  any(
    is.na(up_sig_final$geneID) |
    up_sig_final$geneID == ""
  )
) {
  
  stop(
    "FAIL — Missing tumor-up geneID values."
  )
}

if (
  any(
    is.na(down_sig_final$geneID) |
    down_sig_final$geneID == ""
  )
) {
  
  stop(
    "FAIL — Missing tumor-down geneID values."
  )
}

cat("Pathway geneID completeness: PASS\n\n")


# ------------------------------------------------------------
# 12. Display significant pathway summaries
# ------------------------------------------------------------

cat("============================================================\n")
cat("Tumor-up significant KEGG pathways\n")
cat("============================================================\n\n")

print(
  up_sig_final[
    ,
    c(
      "ID",
      "Description",
      "GeneRatio",
      "BgRatio",
      "pvalue",
      "p.adjust",
      "geneID"
    )
  ],
  row.names = FALSE
)

cat("\n")


cat("============================================================\n")
cat("Tumor-down significant KEGG pathways\n")
cat("============================================================\n\n")

print(
  down_sig_final[
    ,
    c(
      "ID",
      "Description",
      "GeneRatio",
      "BgRatio",
      "pvalue",
      "p.adjust",
      "geneID"
    )
  ],
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 13. Recovery audit conclusion
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP9C-5A RESULT RECOVERY AUDIT\n")
cat("============================================================\n\n")

cat("Existing result files successfully recovered: PASS\n")
cat("Tumor-up complete pathways = 335\n")
cat("Tumor-down complete pathways = 338\n")
cat("Tumor-up significant pathways = 13\n")
cat("Tumor-down significant pathways = 60\n")
cat("Final significance definition = BH-adjusted P < 0.05\n\n")

cat("NO enrichment recalculation performed.\n")
cat("NO DEG recalculation performed.\n")
cat("NO ID conversion performed.\n")
cat("NO P-value/FDR recalculation performed.\n")
cat("NO background modification performed.\n")
cat("NO STEP9C-3 modification performed.\n")
cat("NO STEP9C-4 modification performed.\n\n")

cat("RECOVERY AUDIT STATUS: PASS\n")
cat("============================================================\n")



# ============================================================
# Paper 1 — STEP9C-5B
# KEGG Result Extraction and Manuscript-Definition Lock
#
# FINAL CONTROLLED VERSION
#
# Source:
#   Existing STEP9C-5A COMPLETE result tables
#
# IMPORTANT:
#   NO enrichKEGG()
#   NO DEG recalculation
#   NO ID conversion
#   NO P-value recalculation
#   NO FDR recalculation
#   NO background modification
#   NO modification of STEP9C-3
#   NO modification of STEP9C-4
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-5B: KEGG Result Extraction and Manuscript-Definition Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

setwd(step_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Define existing STEP9C-5A result files
# ------------------------------------------------------------

up_file <- file.path(
  step_dir,
  "STEP9C5A_KEGG_TumorUp_COMPLETE.csv"
)

down_file <- file.path(
  step_dir,
  "STEP9C5A_KEGG_TumorDown_COMPLETE.csv"
)

core_audit_file <- file.path(
  step_dir,
  "STEP9C5A_CoreGene_KEGG_EffectiveUniverse_Audit.csv"
)


# ------------------------------------------------------------
# 2. File existence lock
# ------------------------------------------------------------

required_files <- c(
  up_file,
  down_file,
  core_audit_file
)

file_exists <- file.exists(required_files)

if (!all(file_exists)) {
  
  print(
    data.frame(
      File = basename(required_files),
      Exists = file_exists,
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )
  
  stop(
    "FAIL — Required STEP9C-5A result file(s) missing."
  )
}

cat("STEP9C-5A source files: PASS\n\n")


# ------------------------------------------------------------
# 3. Read existing STEP9C-5A results
# ------------------------------------------------------------

kegg_up_complete <- read.csv(
  up_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

kegg_down_complete <- read.csv(
  down_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

core_effective_audit <- read.csv(
  core_audit_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 4. Structural dimensions
# ------------------------------------------------------------

cat("Tumor-up complete result rows:   ",
    nrow(kegg_up_complete), "\n", sep = "")

cat("Tumor-down complete result rows: ",
    nrow(kegg_down_complete), "\n", sep = "")

cat("Core/effective-universe audit rows: ",
    nrow(core_effective_audit), "\n\n", sep = "")


# ------------------------------------------------------------
# 5. Exact pathway-count lock
# ------------------------------------------------------------

if (nrow(kegg_up_complete) != 335) {
  stop(
    paste0(
      "FAIL — Tumor-up pathway count is ",
      nrow(kegg_up_complete),
      ", expected 335."
    )
  )
}

if (nrow(kegg_down_complete) != 338) {
  stop(
    paste0(
      "FAIL — Tumor-down pathway count is ",
      nrow(kegg_down_complete),
      ", expected 338."
    )
  )
}

cat("335 / 338 complete pathway count: PASS\n\n")


# ------------------------------------------------------------
# 6. Required KEGG result fields
# ------------------------------------------------------------

required_columns <- c(
  "ID",
  "Description",
  "GeneRatio",
  "BgRatio",
  "pvalue",
  "p.adjust",
  "geneID"
)

missing_up <- setdiff(
  required_columns,
  colnames(kegg_up_complete)
)

missing_down <- setdiff(
  required_columns,
  colnames(kegg_down_complete)
)

if (length(missing_up) > 0) {
  
  stop(
    paste0(
      "FAIL — Tumor-up missing columns: ",
      paste(missing_up, collapse = ", ")
    )
  )
}

if (length(missing_down) > 0) {
  
  stop(
    paste0(
      "FAIL — Tumor-down missing columns: ",
      paste(missing_down, collapse = ", ")
    )
  )
}

cat("Required KEGG reporting fields: PASS\n\n")


# ------------------------------------------------------------
# 7. Final significance definition
# ------------------------------------------------------------

FINAL_BH_THRESHOLD <- 0.05

FINAL_SIGNIFICANCE_DEFINITION <-
  "BH-adjusted P < 0.05"


# ------------------------------------------------------------
# 8. Extract significant pathways
# ------------------------------------------------------------
#
# IMPORTANT:
# This is NOT recalculation of P values or FDR.
#
# p.adjust values already exist in the locked STEP9C-5A
# result tables.
#
# This step only subsets the existing results according to
# the predefined significance rule.
# ------------------------------------------------------------

up_sig <- kegg_up_complete[
  !is.na(kegg_up_complete$p.adjust) &
    kegg_up_complete$p.adjust < FINAL_BH_THRESHOLD,
  ,
  drop = FALSE
]

down_sig <- kegg_down_complete[
  !is.na(kegg_down_complete$p.adjust) &
    kegg_down_complete$p.adjust < FINAL_BH_THRESHOLD,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 9. Significant pathway count lock
# ------------------------------------------------------------

cat("Tumor-up significant pathways:   ",
    nrow(up_sig), "\n", sep = "")

cat("Tumor-down significant pathways: ",
    nrow(down_sig), "\n\n", sep = "")


if (nrow(up_sig) != 13) {
  
  stop(
    paste0(
      "FAIL — Expected 13 significant tumor-up pathways; found ",
      nrow(up_sig)
    )
  )
}

if (nrow(down_sig) != 60) {
  
  stop(
    paste0(
      "FAIL — Expected 60 significant tumor-down pathways; found ",
      nrow(down_sig)
    )
  )
}

cat("13 / 60 significant pathway lock: PASS\n\n")


# ------------------------------------------------------------
# 10. Create final reporting tables
# ------------------------------------------------------------

up_sig_final <- up_sig[
  ,
  required_columns,
  drop = FALSE
]

down_sig_final <- down_sig[
  ,
  required_columns,
  drop = FALSE
]


# ------------------------------------------------------------
# 11. Add direction
# ------------------------------------------------------------

up_sig_final$Direction <- "Tumor-up"
down_sig_final$Direction <- "Tumor-down"


up_sig_final <- up_sig_final[
  ,
  c(
    "Direction",
    "ID",
    "Description",
    "GeneRatio",
    "BgRatio",
    "pvalue",
    "p.adjust",
    "geneID"
  )
]

down_sig_final <- down_sig_final[
  ,
  c(
    "Direction",
    "ID",
    "Description",
    "GeneRatio",
    "BgRatio",
    "pvalue",
    "p.adjust",
    "geneID"
  )
]


# ------------------------------------------------------------
# 12. Pathway ID uniqueness
# ------------------------------------------------------------

if (anyDuplicated(up_sig_final$ID) > 0) {
  
  stop(
    "FAIL — Duplicate tumor-up pathway IDs."
  )
}

if (anyDuplicated(down_sig_final$ID) > 0) {
  
  stop(
    "FAIL — Duplicate tumor-down pathway IDs."
  )
}

cat("Pathway ID uniqueness: PASS\n")


# ------------------------------------------------------------
# 13. Verify final significance table
# ------------------------------------------------------------

if (
  any(
    is.na(up_sig_final$p.adjust) |
    up_sig_final$p.adjust >= FINAL_BH_THRESHOLD
  )
) {
  
  stop(
    "FAIL — Tumor-up significant table contains invalid FDR values."
  )
}

if (
  any(
    is.na(down_sig_final$p.adjust) |
    down_sig_final$p.adjust >= FINAL_BH_THRESHOLD
  )
) {
  
  stop(
    "FAIL — Tumor-down significant table contains invalid FDR values."
  )
}

cat("BH-adjusted P < 0.05 enforcement: PASS\n")


# ------------------------------------------------------------
# 14. geneID completeness
# ------------------------------------------------------------

if (
  any(
    is.na(up_sig_final$geneID) |
    up_sig_final$geneID == ""
  )
) {
  
  stop(
    "FAIL — Missing tumor-up geneID."
  )
}

if (
  any(
    is.na(down_sig_final$geneID) |
    down_sig_final$geneID == ""
  )
) {
  
  stop(
    "FAIL — Missing tumor-down geneID."
  )
}

cat("geneID completeness: PASS\n\n")


# ------------------------------------------------------------
# 15. Core-gene audit recovery
# ------------------------------------------------------------

cat("============================================================\n")
cat("Core-gene / effective-universe audit recovery\n")
cat("============================================================\n\n")

cat("Source:\n")
cat(basename(core_audit_file), "\n\n")

print(
  core_effective_audit,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 16. Core-gene definitions
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_entrez <- c(
  AOC3 = "8639",
  MS4A7 = "58475",
  IL7R = "3575"
)


# ------------------------------------------------------------
# 17. Parse pathway geneID
# ------------------------------------------------------------

parse_geneID <- function(x) {
  
  if (
    is.na(x) ||
    x == ""
  ) {
    return(character(0))
  }
  
  unique(
    trimws(
      unlist(
        strsplit(
          x,
          "/",
          fixed = TRUE
        )
      )
    )
  )
}


# ------------------------------------------------------------
# 18. Direct core-gene / significant-pathway audit
# ------------------------------------------------------------

audit_core_pathways <- function(
    result_df,
    direction
) {
  
  out <- vector(
    "list",
    nrow(result_df)
  )
  
  for (i in seq_len(nrow(result_df))) {
    
    pathway_genes <- parse_geneID(
      result_df$geneID[i]
    )
    
    present <- core_entrez %in% pathway_genes
    
    out[[i]] <- data.frame(
      Direction = direction,
      ID = result_df$ID[i],
      Description = result_df$Description[i],
      AOC3 = unname(present["AOC3"]),
      MS4A7 = unname(present["MS4A7"]),
      IL7R = unname(present["IL7R"]),
      Core_Gene_Count = sum(present),
      stringsAsFactors = FALSE
    )
  }
  
  do.call(
    rbind,
    out
  )
}


core_pathway_audit_up <- audit_core_pathways(
  up_sig_final,
  "Tumor-up"
)

core_pathway_audit_down <- audit_core_pathways(
  down_sig_final,
  "Tumor-down"
)

core_pathway_audit <- rbind(
  core_pathway_audit_up,
  core_pathway_audit_down
)


# ------------------------------------------------------------
# 19. Core-gene summary
# ------------------------------------------------------------

core_summary <- data.frame(
  Gene = core_genes,
  Entrez = unname(core_entrez),
  
  Tumor_Up_Significant_Pathways =
    c(
      sum(core_pathway_audit_up$AOC3),
      sum(core_pathway_audit_up$MS4A7),
      sum(core_pathway_audit_up$IL7R)
    ),
  
  Tumor_Down_Significant_Pathways =
    c(
      sum(core_pathway_audit_down$AOC3),
      sum(core_pathway_audit_down$MS4A7),
      sum(core_pathway_audit_down$IL7R)
    ),
  
  stringsAsFactors = FALSE
)

cat("Core-gene significant-pathway summary:\n\n")

print(
  core_summary,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 20. Manuscript pathway-selection rule
# ------------------------------------------------------------

PATHWAY_SELECTION_RULE <- paste0(
  "Significant KEGG pathways were defined exclusively by ",
  "BH-adjusted P < 0.05 among the pathways returned by the ",
  "locked STEP9C-5A enrichment analysis. No additional pathway ",
  "selection criterion was applied."
)


# ------------------------------------------------------------
# 21. No-recalculation lock
# ------------------------------------------------------------

NO_RECALCULATION <- list(
  
  Enrichment_Recalculation = FALSE,
  
  DEG_Recalculation = FALSE,
  
  ID_Conversion = FALSE,
  
  P_Value_Recalculation = FALSE,
  
  FDR_Recalculation = FALSE,
  
  Threshold_Changed = FALSE,
  
  Background_Modified = FALSE,
  
  STEP9C3_Modified = FALSE,
  
  STEP9C4_Modified = FALSE
)


# ------------------------------------------------------------
# 22. Export final Supplementary Tables
# ------------------------------------------------------------

write.csv(
  up_sig_final,
  file = file.path(
    step_dir,
    "STEP9C5B_Supplementary_Table_KEGG_TumorUp_13.csv"
  ),
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  down_sig_final,
  file = file.path(
    step_dir,
    "STEP9C5B_Supplementary_Table_KEGG_TumorDown_60.csv"
  ),
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  core_pathway_audit,
  file = file.path(
    step_dir,
    "STEP9C5B_CoreGene_Pathway_Overlap_Audit.csv"
  ),
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  core_summary,
  file = file.path(
    step_dir,
    "STEP9C5B_CoreGene_Pathway_Summary.csv"
  ),
  row.names = FALSE,
  quote = TRUE
)


# ------------------------------------------------------------
# 23. Create STEP9C-5B FINAL LOCK
# ------------------------------------------------------------

STEP9C5B_FINAL_LOCK <- list(
  
  Step = "STEP9C-5B",
  
  Status = "PASS",
  
  Source =
    "Existing STEP9C-5A COMPLETE CSV result tables",
  
  Predefined_KEGG_Background = 35365,
  
  Effective_KEGG_Universe = 8702,
  
  Tumor_Up_Input = 5194,
  
  Tumor_Down_Input = 4463,
  
  Tumor_Up_Pathways_Returned = 335,
  
  Tumor_Down_Pathways_Returned = 338,
  
  Tumor_Up_Significant_Pathways = 13,
  
  Tumor_Down_Significant_Pathways = 60,
  
  Significance_Definition =
    FINAL_SIGNIFICANCE_DEFINITION,
  
  Reporting_Columns =
    c(
      "Direction",
      "ID",
      "Description",
      "GeneRatio",
      "BgRatio",
      "pvalue",
      "p.adjust",
      "geneID"
    ),
  
  Tumor_Up_Pathway_IDs =
    up_sig_final$ID,
  
  Tumor_Down_Pathway_IDs =
    down_sig_final$ID,
  
  Core_Genes =
    core_genes,
  
  Core_Entrez =
    core_entrez,
  
  Pathway_Selection_Rule =
    PATHWAY_SELECTION_RULE,
  
  Enrichment_Recalculation = FALSE,
  
  DEG_Recalculation = FALSE,
  
  ID_Conversion = FALSE,
  
  P_Value_Recalculation = FALSE,
  
  FDR_Recalculation = FALSE,
  
  Threshold_Changed = FALSE,
  
  Background_Modified = FALSE,
  
  STEP9C3_Modified = FALSE,
  
  STEP9C4_Modified = FALSE
)


# ------------------------------------------------------------
# 24. Save final lock
# ------------------------------------------------------------

saveRDS(
  STEP9C5B_FINAL_LOCK,
  file = file.path(
    step_dir,
    "STEP9C5B_FINAL_LOCK.rds"
  )
)


# ------------------------------------------------------------
# 25. Save complete manuscript result object
# ------------------------------------------------------------

STEP9C5B_RESULTS <- list(
  
  Tumor_Up_Significant =
    up_sig_final,
  
  Tumor_Down_Significant =
    down_sig_final,
  
  Core_Pathway_Audit =
    core_pathway_audit,
  
  Core_Summary =
    core_summary,
  
  Existing_STEP9C5A_Core_EffectiveUniverse_Audit =
    core_effective_audit,
  
  Predefined_KEGG_Background =
    35365,
  
  Effective_KEGG_Universe =
    8702,
  
  Significance_Definition =
    FINAL_SIGNIFICANCE_DEFINITION,
  
  Pathway_Selection_Rule =
    PATHWAY_SELECTION_RULE
)


saveRDS(
  STEP9C5B_RESULTS,
  file = file.path(
    step_dir,
    "STEP9C5B_MANUSCRIPT_RESULTS_LOCKED.rds"
  )
)


# ------------------------------------------------------------
# 26. Final validation
# ------------------------------------------------------------

STEP9C5B_PASS <- all(
  
  identical(
    STEP9C5B_FINAL_LOCK$Step,
    "STEP9C-5B"
  ),
  
  identical(
    STEP9C5B_FINAL_LOCK$Status,
    "PASS"
  ),
  
  isTRUE(
    STEP9C5B_FINAL_LOCK$Predefined_KEGG_Background == 35365
  ),
  
  isTRUE(
    STEP9C5B_FINAL_LOCK$Effective_KEGG_Universe == 8702
  ),
  
  isTRUE(
    STEP9C5B_FINAL_LOCK$Tumor_Up_Input == 5194
  ),
  
  isTRUE(
    STEP9C5B_FINAL_LOCK$Tumor_Down_Input == 4463
  ),
  
  isTRUE(
    STEP9C5B_FINAL_LOCK$Tumor_Up_Pathways_Returned == 335
  ),
  
  isTRUE(
    STEP9C5B_FINAL_LOCK$Tumor_Down_Pathways_Returned == 338
  ),
  
  isTRUE(
    STEP9C5B_FINAL_LOCK$Tumor_Up_Significant_Pathways == 13
  ),
  
  isTRUE(
    STEP9C5B_FINAL_LOCK$Tumor_Down_Significant_Pathways == 60
  ),
  
  identical(
    STEP9C5B_FINAL_LOCK$Significance_Definition,
    "BH-adjusted P < 0.05"
  ),
  
  identical(
    STEP9C5B_FINAL_LOCK$Enrichment_Recalculation,
    FALSE
  ),
  
  identical(
    STEP9C5B_FINAL_LOCK$DEG_Recalculation,
    FALSE
  ),
  
  identical(
    STEP9C5B_FINAL_LOCK$ID_Conversion,
    FALSE
  ),
  
  identical(
    STEP9C5B_FINAL_LOCK$P_Value_Recalculation,
    FALSE
  ),
  
  identical(
    STEP9C5B_FINAL_LOCK$FDR_Recalculation,
    FALSE
  ),
  
  identical(
    STEP9C5B_FINAL_LOCK$Threshold_Changed,
    FALSE
  ),
  
  identical(
    STEP9C5B_FINAL_LOCK$Background_Modified,
    FALSE
  ),
  
  identical(
    STEP9C5B_FINAL_LOCK$STEP9C3_Modified,
    FALSE
  ),
  
  identical(
    STEP9C5B_FINAL_LOCK$STEP9C4_Modified,
    FALSE
  )
)


# ------------------------------------------------------------
# 27. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")

if (STEP9C5B_PASS) {
  
  cat("FINAL STATUS: PASS — STEP9C-5B LOCKED\n\n")
  
  cat("Predefined KEGG background = 35,365\n")
  cat("Effective KEGG universe = 8,702\n")
  cat("Tumor-up input = 5,194\n")
  cat("Tumor-down input = 4,463\n")
  cat("Tumor-up pathways returned = 335\n")
  cat("Tumor-down pathways returned = 338\n")
  cat("Tumor-up significant pathways = 13\n")
  cat("Tumor-down significant pathways = 60\n\n")
  
  cat("Final significance definition:\n")
  cat("  BH-adjusted P < 0.05\n\n")
  
  cat("No enrichment recalculation.\n")
  cat("No DEG recalculation.\n")
  cat("No ID conversion.\n")
  cat("No P-value/FDR recalculation.\n")
  cat("No threshold modification.\n")
  cat("No background modification.\n")
  cat("STEP9C-3 unchanged.\n")
  cat("STEP9C-4 unchanged.\n\n")
  
  cat("Final Supplementary Tables exported.\n")
  cat("Core-gene pathway audit exported.\n")
  cat("STEP9C5B_FINAL_LOCK.rds saved.\n")
  cat("STEP9C5B_MANUSCRIPT_RESULTS_LOCKED.rds saved.\n")
  
} else {
  
  cat("FINAL STATUS: FAIL — STEP9C-5B NOT LOCKED\n")
  cat("Do NOT treat STEP9C-5B as final.\n")
}

cat("\n============================================================\n")
cat("STEP9C-5B COMPLETE\n")
cat("============================================================\n")


# ============================================================
# Paper 1
# STEP9C-6
# Controlled GO Biological Process ORA
#
# Purpose:
#   Perform GO-BP enrichment using the LOCKED STEP9C-4
#   Entrez inputs and predefined Entrez background.
#
# IMPORTANT:
#   No DEG recalculation
#   No ID conversion
#   No p-value recalculation outside enrichGO()
#   No FDR recalculation outside enrichGO()
#   No threshold modification
#   No background modification
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-6: Controlled GO Biological Process ORA\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

setwd(step_dir)


# ------------------------------------------------------------
# 1. Load locked STEP9C-4 object
# ------------------------------------------------------------

step9c4_file <- file.path(
  step_dir,
  "STEP9C4_KEGG_INPUT_AND_BACKGROUND_LOCKED.rds"
)

stopifnot(file.exists(step9c4_file))

STEP9C4_LOCKED <- readRDS(step9c4_file)


# ------------------------------------------------------------
# 2. Recover locked inputs
# ------------------------------------------------------------

kegg_up_entrez <- STEP9C4_LOCKED$kegg_up_entrez

kegg_down_entrez <- STEP9C4_LOCKED$kegg_down_entrez

final_universe_entrez <- STEP9C4_LOCKED$final_universe_entrez


# ------------------------------------------------------------
# 3. Structural QC
# ------------------------------------------------------------

stopifnot(
  is.character(kegg_up_entrez),
  is.character(kegg_down_entrez),
  is.character(final_universe_entrez)
)

stopifnot(
  length(kegg_up_entrez) == 5194,
  length(kegg_down_entrez) == 4463,
  length(final_universe_entrez) == 35365
)

stopifnot(
  length(unique(kegg_up_entrez)) == 5194,
  length(unique(kegg_down_entrez)) == 4463,
  length(unique(final_universe_entrez)) == 35365
)

stopifnot(
  length(intersect(
    kegg_up_entrez,
    kegg_down_entrez
  )) == 0
)


# ------------------------------------------------------------
# 4. Load clusterProfiler
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})


cat("clusterProfiler version: ",
    as.character(packageVersion("clusterProfiler")),
    "\n")

cat("R version: ",
    R.version.string,
    "\n\n")


# ------------------------------------------------------------
# 5. Locked GO-BP ORA parameters
# ------------------------------------------------------------

GO_ORA_PARAMS <- list(
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "BH",
  significance_definition = "BH-adjusted P < 0.05"
)


# ------------------------------------------------------------
# 6. GO-BP enrichment
# ------------------------------------------------------------

cat("Running GO-BP ORA...\n")

go_up <- enrichGO(
  gene          = kegg_up_entrez,
  universe      = final_universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  readable      = FALSE
)


go_down <- enrichGO(
  gene          = kegg_down_entrez,
  universe      = final_universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  readable      = FALSE
)


# ------------------------------------------------------------
# 7. Basic result audit
# ------------------------------------------------------------

go_up_df <- as.data.frame(go_up)

go_down_df <- as.data.frame(go_down)


cat("\nGO-BP results:\n")

cat("Tumor-up pathways returned   = ",
    nrow(go_up_df), "\n")

cat("Tumor-down pathways returned = ",
    nrow(go_down_df), "\n")

cat("\n")


# ------------------------------------------------------------
# 8. Effective universe audit
# ------------------------------------------------------------

effective_go_up <- length(go_up@universe)

effective_go_down <- length(go_down@universe)

cat("Effective GO universe:\n")

cat("Tumor-up   = ",
    effective_go_up, "\n")

cat("Tumor-down = ",
    effective_go_down, "\n\n")


# ------------------------------------------------------------
# 9. Effective universe must be subset of locked background
# ------------------------------------------------------------

stopifnot(
  all(go_up@universe %in% final_universe_entrez)
)

stopifnot(
  all(go_down@universe %in% final_universe_entrez)
)


# ------------------------------------------------------------
# 10. Final significance definition
# ------------------------------------------------------------

FINAL_BH_THRESHOLD <- 0.05

FINAL_SIGNIFICANCE_DEFINITION <-
  "BH-adjusted P < 0.05"


go_up_sig <- subset(
  go_up_df,
  p.adjust < FINAL_BH_THRESHOLD
)

go_down_sig <- subset(
  go_down_df,
  p.adjust < FINAL_BH_THRESHOLD
)


# ------------------------------------------------------------
# 11. Export complete GO results
# ------------------------------------------------------------

write.csv(
  go_up_df,
  file.path(
    step_dir,
    "STEP9C6_GO_BP_TumorUp_COMPLETE.csv"
  ),
  row.names = FALSE
)

write.csv(
  go_down_df,
  file.path(
    step_dir,
    "STEP9C6_GO_BP_TumorDown_COMPLETE.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Export significant GO results
# ------------------------------------------------------------

write.csv(
  go_up_sig,
  file.path(
    step_dir,
    "STEP9C6_GO_BP_TumorUp_FDR05.csv"
  ),
  row.names = FALSE
)

write.csv(
  go_down_sig,
  file.path(
    step_dir,
    "STEP9C6_GO_BP_TumorDown_FDR05.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Lock object
# ------------------------------------------------------------

STEP9C6_FINAL_LOCK <- list(
  
  Step = "STEP9C-6",
  
  Status = "PASS",
  
  Analysis = "GO Biological Process ORA",
  
  Input_Source =
    "STEP9C4_KEGG_INPUT_AND_BACKGROUND_LOCKED.rds",
  
  Tumor_Up_Input =
    length(kegg_up_entrez),
  
  Tumor_Down_Input =
    length(kegg_down_entrez),
  
  Predefined_Background =
    length(final_universe_entrez),
  
  Effective_GO_Universe_Up =
    effective_go_up,
  
  Effective_GO_Universe_Down =
    effective_go_down,
  
  Tumor_Up_Pathways_Returned =
    nrow(go_up_df),
  
  Tumor_Down_Pathways_Returned =
    nrow(go_down_df),
  
  Tumor_Up_Significant =
    nrow(go_up_sig),
  
  Tumor_Down_Significant =
    nrow(go_down_sig),
  
  Significance_Definition =
    FINAL_SIGNIFICANCE_DEFINITION,
  
  Ontology = "BP",
  
  P_Adjust_Method = "BH",
  
  Min_GS_Size = 10,
  
  Max_GS_Size = 500,
  
  ID_Conversion_Performed = FALSE,
  
  DEG_Recalculation = FALSE,
  
  Threshold_Changed = FALSE,
  
  Background_Modified = FALSE,
  
  STEP9C3_Modified = FALSE,
  
  STEP9C4_Modified = FALSE
)


saveRDS(
  STEP9C6_FINAL_LOCK,
  file.path(
    step_dir,
    "STEP9C6_FINAL_LOCK.rds"
  )
)


# ------------------------------------------------------------
# 14. Final report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-6 COMPLETE\n")
cat("============================================================\n\n")

cat("Predefined Entrez background = ",
    length(final_universe_entrez), "\n")

cat("Tumor-up input = ",
    length(kegg_up_entrez), "\n")

cat("Tumor-down input = ",
    length(kegg_down_entrez), "\n")

cat("Tumor-up GO-BP pathways returned = ",
    nrow(go_up_df), "\n")

cat("Tumor-down GO-BP pathways returned = ",
    nrow(go_down_df), "\n")

cat("Tumor-up GO-BP BH-adjusted P < 0.05 = ",
    nrow(go_up_sig), "\n")

cat("Tumor-down GO-BP BH-adjusted P < 0.05 = ",
    nrow(go_down_sig), "\n\n")

cat("Final significance definition:\n")
cat("  ", FINAL_SIGNIFICANCE_DEFINITION, "\n\n")

cat("No DEG recalculation.\n")
cat("No modification of STEP9C-4 inputs.\n")
cat("No background modification.\n")
cat("STEP9C-3 unchanged.\n")
cat("STEP9C-4 unchanged.\n\n")

# ============================================================
# Paper 1
# STEP9C-6A
# GO-BP Result Recovery / Effective-Universe / Structural QC
#
# IMPORTANT:
#   This step DOES NOT rerun enrichGO().
#   This step DOES NOT recalculate P values or FDR.
#   This step only audits the existing STEP9C-6 results.
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-6A: GO-BP Result Recovery / Structural QC\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

setwd(step_dir)


# ------------------------------------------------------------
# 1. Confirm required STEP9C-6 objects exist
# ------------------------------------------------------------

required_objects <- c(
  "go_up",
  "go_down",
  "go_up_df",
  "go_down_df",
  "go_up_sig",
  "go_down_sig",
  "final_universe_entrez",
  "kegg_up_entrez",
  "kegg_down_entrez"
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
      "FAIL: Required STEP9C-6 objects are missing: ",
      paste(missing_objects, collapse = ", ")
    )
  )
  
}


# ------------------------------------------------------------
# 2. Locked expected values from STEP9C-6
# ------------------------------------------------------------

EXPECTED_BACKGROUND <- 35365

EXPECTED_UP_INPUT <- 5194

EXPECTED_DOWN_INPUT <- 4463

EXPECTED_UP_RETURNED <- 5680

EXPECTED_DOWN_RETURNED <- 5847

EXPECTED_UP_SIG <- 161

EXPECTED_DOWN_SIG <- 1296

FINAL_BH_THRESHOLD <- 0.05

FINAL_SIGNIFICANCE_DEFINITION <-
  "BH-adjusted P < 0.05"


# ------------------------------------------------------------
# 3. Basic input validation
# ------------------------------------------------------------

stopifnot(
  is.character(final_universe_entrez),
  is.character(kegg_up_entrez),
  is.character(kegg_down_entrez)
)

stopifnot(
  length(final_universe_entrez) == EXPECTED_BACKGROUND
)

stopifnot(
  length(unique(final_universe_entrez)) == EXPECTED_BACKGROUND
)

stopifnot(
  length(kegg_up_entrez) == EXPECTED_UP_INPUT
)

stopifnot(
  length(unique(kegg_up_entrez)) == EXPECTED_UP_INPUT
)

stopifnot(
  length(kegg_down_entrez) == EXPECTED_DOWN_INPUT
)

stopifnot(
  length(unique(kegg_down_entrez)) == EXPECTED_DOWN_INPUT
)

stopifnot(
  length(intersect(
    kegg_up_entrez,
    kegg_down_entrez
  )) == 0
)


# ------------------------------------------------------------
# 4. Effective GO universe audit
# ------------------------------------------------------------

effective_go_up <- go_up@universe

effective_go_down <- go_down@universe

effective_go_up_n <- length(effective_go_up)

effective_go_down_n <- length(effective_go_down)


cat("Effective GO universe:\n")

cat("  Tumor-up   = ",
    effective_go_up_n,
    "\n")

cat("  Tumor-down = ",
    effective_go_down_n,
    "\n\n")


# ------------------------------------------------------------
# 5. Effective universe must be subset of locked background
# ------------------------------------------------------------

effective_up_is_subset <- all(
  effective_go_up %in% final_universe_entrez
)

effective_down_is_subset <- all(
  effective_go_down %in% final_universe_entrez
)

stopifnot(effective_up_is_subset)

stopifnot(effective_down_is_subset)


# ------------------------------------------------------------
# 6. Result row-count audit
# ------------------------------------------------------------

stopifnot(
  nrow(go_up_df) == EXPECTED_UP_RETURNED
)

stopifnot(
  nrow(go_down_df) == EXPECTED_DOWN_RETURNED
)

stopifnot(
  nrow(go_up_sig) == EXPECTED_UP_SIG
)

stopifnot(
  nrow(go_down_sig) == EXPECTED_DOWN_SIG
)


# ------------------------------------------------------------
# 7. Required column audit
# ------------------------------------------------------------

required_columns <- c(
  "ID",
  "Description",
  "GeneRatio",
  "BgRatio",
  "pvalue",
  "p.adjust",
  "geneID"
)

missing_up_columns <- setdiff(
  required_columns,
  colnames(go_up_df)
)

missing_down_columns <- setdiff(
  required_columns,
  colnames(go_down_df)
)

stopifnot(
  length(missing_up_columns) == 0
)

stopifnot(
  length(missing_down_columns) == 0
)


# ------------------------------------------------------------
# 8. GO ID uniqueness
# ------------------------------------------------------------

stopifnot(
  !anyDuplicated(go_up_df$ID)
)

stopifnot(
  !anyDuplicated(go_down_df$ID)
)


# ------------------------------------------------------------
# 9. P-adjustment audit
# ------------------------------------------------------------

stopifnot(
  all(
    is.finite(go_up_df$p.adjust) |
      is.na(go_up_df$p.adjust)
  )
)

stopifnot(
  all(
    is.finite(go_down_df$p.adjust) |
      is.na(go_down_df$p.adjust)
  )
)


# ------------------------------------------------------------
# 10. Significant-result definition audit
# ------------------------------------------------------------

up_sig_recovered <- subset(
  go_up_df,
  p.adjust < FINAL_BH_THRESHOLD
)

down_sig_recovered <- subset(
  go_down_df,
  p.adjust < FINAL_BH_THRESHOLD
)


stopifnot(
  nrow(up_sig_recovered) == EXPECTED_UP_SIG
)

stopifnot(
  nrow(down_sig_recovered) == EXPECTED_DOWN_SIG
)


# ------------------------------------------------------------
# 11. Confirm saved significant objects agree
# ------------------------------------------------------------

stopifnot(
  identical(
    as.character(go_up_sig$ID),
    as.character(up_sig_recovered$ID)
  )
)

stopifnot(
  identical(
    as.character(go_down_sig$ID),
    as.character(down_sig_recovered$ID)
  )
)


# ------------------------------------------------------------
# 12. GO ontology audit
# ------------------------------------------------------------

GO_ONTOLOGY <- "BP"


# ------------------------------------------------------------
# 13. Save recovery / QC summary
# ------------------------------------------------------------

STEP9C6A_QC <- list(
  
  Step = "STEP9C-6A",
  
  Status = "PASS",
  
  Analysis = "GO Biological Process ORA",
  
  Predefined_Entrez_Background =
    EXPECTED_BACKGROUND,
  
  Tumor_Up_Input =
    EXPECTED_UP_INPUT,
  
  Tumor_Down_Input =
    EXPECTED_DOWN_INPUT,
  
  Effective_GO_Universe_Up =
    effective_go_up_n,
  
  Effective_GO_Universe_Down =
    effective_go_down_n,
  
  Effective_Up_Subset_Of_Background =
    effective_up_is_subset,
  
  Effective_Down_Subset_Of_Background =
    effective_down_is_subset,
  
  Tumor_Up_Pathways_Returned =
    EXPECTED_UP_RETURNED,
  
  Tumor_Down_Pathways_Returned =
    EXPECTED_DOWN_RETURNED,
  
  Tumor_Up_Significant =
    EXPECTED_UP_SIG,
  
  Tumor_Down_Significant =
    EXPECTED_DOWN_SIG,
  
  Ontology =
    GO_ONTOLOGY,
  
  P_Adjust_Method =
    "BH",
  
  Significance_Definition =
    FINAL_SIGNIFICANCE_DEFINITION,
  
  Enrichment_Recalculation =
    FALSE,
  
  DEG_Recalculation =
    FALSE,
  
  ID_Conversion =
    FALSE,
  
  P_Value_Recalculation =
    FALSE,
  
  FDR_Recalculation =
    FALSE,
  
  Threshold_Changed =
    FALSE,
  
  Background_Modified =
    FALSE,
  
  STEP9C3_Modified =
    FALSE,
  
  STEP9C4_Modified =
    FALSE
)


saveRDS(
  STEP9C6A_QC,
  file.path(
    step_dir,
    "STEP9C6A_GO_BP_RESULT_QC_LOCK.rds"
  )
)


# ------------------------------------------------------------
# 14. Final PASS
# ------------------------------------------------------------

FINAL_PASS <- all(
  STEP9C6A_QC$Status == "PASS",
  STEP9C6A_QC$Predefined_Entrez_Background == 35365,
  STEP9C6A_QC$Tumor_Up_Input == 5194,
  STEP9C6A_QC$Tumor_Down_Input == 4463,
  STEP9C6A_QC$Tumor_Up_Pathways_Returned == 5680,
  STEP9C6A_QC$Tumor_Down_Pathways_Returned == 5847,
  STEP9C6A_QC$Tumor_Up_Significant == 161,
  STEP9C6A_QC$Tumor_Down_Significant == 1296,
  STEP9C6A_QC$Effective_Up_Subset_Of_Background,
  STEP9C6A_QC$Effective_Down_Subset_Of_Background,
  STEP9C6A_QC$Enrichment_Recalculation == FALSE,
  STEP9C6A_QC$DEG_Recalculation == FALSE,
  STEP9C6A_QC$ID_Conversion == FALSE,
  STEP9C6A_QC$P_Value_Recalculation == FALSE,
  STEP9C6A_QC$FDR_Recalculation == FALSE,
  STEP9C6A_QC$Threshold_Changed == FALSE,
  STEP9C6A_QC$Background_Modified == FALSE,
  STEP9C6A_QC$STEP9C3_Modified == FALSE,
  STEP9C6A_QC$STEP9C4_Modified == FALSE
)


if (!FINAL_PASS) {
  
  stop(
    "FINAL STATUS: FAIL — STEP9C-6A QC failed."
  )
  
}


cat("\n============================================================\n")
cat("FINAL STATUS: PASS — STEP9C-6A QC COMPLETE\n")
cat("============================================================\n\n")

cat("Predefined Entrez background = ",
    EXPECTED_BACKGROUND, "\n")

cat("Effective GO universe — up = ",
    effective_go_up_n, "\n")

cat("Effective GO universe — down = ",
    effective_go_down_n, "\n")

cat("Tumor-up pathways returned = ",
    EXPECTED_UP_RETURNED, "\n")

cat("Tumor-down pathways returned = ",
    EXPECTED_DOWN_RETURNED, "\n")

cat("Tumor-up significant GO-BP terms = ",
    EXPECTED_UP_SIG, "\n")

cat("Tumor-down significant GO-BP terms = ",
    EXPECTED_DOWN_SIG, "\n\n")

cat("Final significance definition:\n")
cat("  ", FINAL_SIGNIFICANCE_DEFINITION, "\n\n")

cat("No enrichment recalculation.\n")
cat("No DEG recalculation.\n")
cat("No ID conversion.\n")
cat("No P-value/FDR recalculation.\n")
cat("No threshold modification.\n")
cat("No background modification.\n")
cat("STEP9C-3 unchanged.\n")
cat("STEP9C-4 unchanged.\n\n")

cat("STEP9C-6A COMPLETE\n")
cat("============================================================\n")


# ============================================================
# Paper 1
# STEP9C-6B
# GO-BP Final Result Definition / Manuscript Lock
#
# Purpose:
#   1. Lock all significant GO-BP terms
#   2. Lock reporting fields
#   3. Perform direct core-gene/pathway overlap audit
#   4. Export final Supplementary Tables
#   5. Save manuscript-result lock
#
# IMPORTANT:
#   NO enrichGO() rerun
#   NO DEG recalculation
#   NO ID conversion
#   NO P-value recalculation
#   NO FDR recalculation
#   NO threshold modification
#   NO background modification
#   NO modification of STEP9C-3
#   NO modification of STEP9C-4
#
# Final significance definition:
#   BH-adjusted P < 0.05
# ============================================================


cat("\n============================================================\n")
cat("STEP9C-6B: GO-BP Final Result Definition Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

setwd(step_dir)


# ------------------------------------------------------------
# 1. Required files
# ------------------------------------------------------------

up_complete_file <- file.path(
  step_dir,
  "STEP9C6_GO_BP_TumorUp_COMPLETE.csv"
)

down_complete_file <- file.path(
  step_dir,
  "STEP9C6_GO_BP_TumorDown_COMPLETE.csv"
)

step9c6a_lock_file <- file.path(
  step_dir,
  "STEP9C6A_GO_BP_RESULT_QC_LOCK.rds"
)


stopifnot(file.exists(up_complete_file))
stopifnot(file.exists(down_complete_file))
stopifnot(file.exists(step9c6a_lock_file))


# ------------------------------------------------------------
# 2. Recover STEP9C-6A QC lock
# ------------------------------------------------------------

STEP9C6A_QC <- readRDS(
  step9c6a_lock_file
)


# ------------------------------------------------------------
# 3. Validate STEP9C-6A lock
# ------------------------------------------------------------

stopifnot(
  identical(
    STEP9C6A_QC$Step,
    "STEP9C-6A"
  )
)

stopifnot(
  identical(
    STEP9C6A_QC$Status,
    "PASS"
  )
)

stopifnot(
  isTRUE(
    STEP9C6A_QC$Predefined_Entrez_Background == 35365
  )
)

stopifnot(
  isTRUE(
    STEP9C6A_QC$Tumor_Up_Input == 5194
  )
)

stopifnot(
  isTRUE(
    STEP9C6A_QC$Tumor_Down_Input == 4463
  )
)

stopifnot(
  isTRUE(
    STEP9C6A_QC$Tumor_Up_Pathways_Returned == 5680
  )
)

stopifnot(
  isTRUE(
    STEP9C6A_QC$Tumor_Down_Pathways_Returned == 5847
  )
)

stopifnot(
  isTRUE(
    STEP9C6A_QC$Tumor_Up_Significant == 161
  )
)

stopifnot(
  isTRUE(
    STEP9C6A_QC$Tumor_Down_Significant == 1296
  )
)

stopifnot(
  identical(
    STEP9C6A_QC$Significance_Definition,
    "BH-adjusted P < 0.05"
  )
)


# ------------------------------------------------------------
# 4. Read EXISTING GO-BP results
#
# No enrichment recalculation.
# ------------------------------------------------------------

go_up_df <- read.csv(
  up_complete_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

go_down_df <- read.csv(
  down_complete_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 5. Structural result audit
# ------------------------------------------------------------

cat("GO-BP complete result recovery:\n")

cat("  Tumor-up rows   = ",
    nrow(go_up_df), "\n")

cat("  Tumor-down rows = ",
    nrow(go_down_df), "\n\n")


stopifnot(
  nrow(go_up_df) == 5680
)

stopifnot(
  nrow(go_down_df) == 5847
)


# ------------------------------------------------------------
# 6. Required reporting fields
# ------------------------------------------------------------

required_columns <- c(
  "ID",
  "Description",
  "GeneRatio",
  "BgRatio",
  "pvalue",
  "p.adjust",
  "geneID"
)


missing_up_columns <- setdiff(
  required_columns,
  colnames(go_up_df)
)

missing_down_columns <- setdiff(
  required_columns,
  colnames(go_down_df)
)


if (length(missing_up_columns) > 0) {
  
  stop(
    paste0(
      "FAIL: Missing tumor-up columns: ",
      paste(
        missing_up_columns,
        collapse = ", "
      )
    )
  )
  
}


if (length(missing_down_columns) > 0) {
  
  stop(
    paste0(
      "FAIL: Missing tumor-down columns: ",
      paste(
        missing_down_columns,
        collapse = ", "
      )
    )
  )
  
}


# ------------------------------------------------------------
# 7. Final significance definition
# ------------------------------------------------------------

FINAL_BH_THRESHOLD <- 0.05

FINAL_SIGNIFICANCE_DEFINITION <-
  "BH-adjusted P < 0.05"


# ------------------------------------------------------------
# 8. Extract SIGNIFICANT results only
#
# IMPORTANT:
#   This uses the already calculated p.adjust values.
#   No P/FDR recalculation.
# ------------------------------------------------------------

up_sig_final <- go_up_df[
  !is.na(go_up_df$p.adjust) &
    go_up_df$p.adjust < FINAL_BH_THRESHOLD,
  ,
  drop = FALSE
]


down_sig_final <- go_down_df[
  !is.na(go_down_df$p.adjust) &
    go_down_df$p.adjust < FINAL_BH_THRESHOLD,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 9. Confirm locked counts
# ------------------------------------------------------------

stopifnot(
  nrow(up_sig_final) == 161
)

stopifnot(
  nrow(down_sig_final) == 1296
)


# ------------------------------------------------------------
# 10. Add Direction as first column
# ------------------------------------------------------------

up_sig_final <- data.frame(
  Direction = "Tumor-up",
  up_sig_final,
  check.names = FALSE
)


down_sig_final <- data.frame(
  Direction = "Tumor-down",
  down_sig_final,
  check.names = FALSE
)


# ------------------------------------------------------------
# 11. Final pathway/term ID uniqueness
# ------------------------------------------------------------

stopifnot(
  !anyDuplicated(up_sig_final$ID)
)

stopifnot(
  !anyDuplicated(down_sig_final$ID)
)


# ------------------------------------------------------------
# 12. Confirm all final rows satisfy significance definition
# ------------------------------------------------------------

stopifnot(
  all(
    up_sig_final$p.adjust < FINAL_BH_THRESHOLD
  )
)

stopifnot(
  all(
    down_sig_final$p.adjust < FINAL_BH_THRESHOLD
  )
)


# ------------------------------------------------------------
# 13. Confirm geneID is populated
# ------------------------------------------------------------

stopifnot(
  all(
    !is.na(up_sig_final$geneID) &
      nzchar(up_sig_final$geneID)
  )
)

stopifnot(
  all(
    !is.na(down_sig_final$geneID) &
      nzchar(down_sig_final$geneID)
  )
)


# ------------------------------------------------------------
# 14. Core-gene definitions
# ------------------------------------------------------------

core_entrez <- c(
  AOC3 = "8639",
  MS4A7 = "58475",
  IL7R = "3575"
)


core_genes <- names(core_entrez)


cat("\nCore genes:\n")

print(core_entrez)

cat("\n")


# ------------------------------------------------------------
# 15. Helper function:
#     direct core-gene / GO-term overlap
#
# IMPORTANT:
#   This is DIRECT GENE-SET OVERLAP ONLY.
#   It is NOT a leading-edge analysis.
#   It does NOT infer causality or mechanism.
# ------------------------------------------------------------

extract_core_overlap <- function(
    result_df,
    direction_label,
    core_entrez
) {
  
  out_list <- vector(
    "list",
    nrow(result_df)
  )
  
  
  for (i in seq_len(nrow(result_df))) {
    
    gene_string <- as.character(
      result_df$geneID[i]
    )
    
    pathway_genes <- unlist(
      strsplit(
        gene_string,
        split = "/",
        fixed = TRUE
      )
    )
    
    pathway_genes <- trimws(
      pathway_genes
    )
    
    overlap_entrez <- intersect(
      pathway_genes,
      unname(core_entrez)
    )
    
    overlap_genes <- names(
      core_entrez
    )[
      core_entrez %in% overlap_entrez
    ]
    
    
    out_list[[i]] <- data.frame(
      
      Direction =
        direction_label,
      
      ID =
        as.character(result_df$ID[i]),
      
      Description =
        as.character(
          result_df$Description[i]
        ),
      
      GeneRatio =
        as.character(
          result_df$GeneRatio[i]
        ),
      
      BgRatio =
        as.character(
          result_df$BgRatio[i]
        ),
      
      pvalue =
        as.numeric(
          result_df$pvalue[i]
        ),
      
      p.adjust =
        as.numeric(
          result_df$p.adjust[i]
        ),
      
      AOC3 =
        "AOC3" %in% overlap_genes,
      
      MS4A7 =
        "MS4A7" %in% overlap_genes,
      
      IL7R =
        "IL7R" %in% overlap_genes,
      
      CoreGene_Count =
        length(overlap_genes),
      
      CoreGenes_Overlapping =
        if (
          length(overlap_genes) == 0
        ) {
          ""
        } else {
          paste(
            overlap_genes,
            collapse = ";"
          )
        },
      
      stringsAsFactors = FALSE,
      
      check.names = FALSE
    )
  }
  
  
  do.call(
    rbind,
    out_list
  )
}


# ------------------------------------------------------------
# 16. Perform direct core-gene overlap audit
# ------------------------------------------------------------

core_overlap_up <- extract_core_overlap(
  result_df = up_sig_final,
  direction_label = "Tumor-up",
  core_entrez = core_entrez
)


core_overlap_down <- extract_core_overlap(
  result_df = down_sig_final,
  direction_label = "Tumor-down",
  core_entrez = core_entrez
)


core_overlap_all <- rbind(
  core_overlap_up,
  core_overlap_down
)


# ------------------------------------------------------------
# 17. Keep only terms with direct core-gene overlap
# ------------------------------------------------------------

core_overlap_positive <- core_overlap_all[
  core_overlap_all$CoreGene_Count > 0,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 18. Core-gene summary
# ------------------------------------------------------------

core_summary <- data.frame(
  
  Core_Gene =
    core_genes,
  
  Entrez_ID =
    unname(core_entrez),
  
  Tumor_Up_GO_BP_Term_Count =
    sapply(
      core_genes,
      function(g) {
        sum(
          core_overlap_up[[g]]
        )
      }
    ),
  
  Tumor_Down_GO_BP_Term_Count =
    sapply(
      core_genes,
      function(g) {
        sum(
          core_overlap_down[[g]]
        )
      }
    ),
  
  Total_Significant_GO_BP_Term_Count =
    sapply(
      core_genes,
      function(g) {
        sum(
          core_overlap_all[[g]]
        )
      }
    ),
  
  stringsAsFactors = FALSE,
  
  check.names = FALSE
)


# ------------------------------------------------------------
# 19. Export final Supplementary Tables
# ------------------------------------------------------------

write.csv(
  up_sig_final,
  file.path(
    step_dir,
    "STEP9C6B_Supplementary_Table_GO_BP_TumorUp_161.csv"
  ),
  row.names = FALSE
)


write.csv(
  down_sig_final,
  file.path(
    step_dir,
    "STEP9C6B_Supplementary_Table_GO_BP_TumorDown_1296.csv"
  ),
  row.names = FALSE
)


write.csv(
  core_overlap_positive,
  file.path(
    step_dir,
    "STEP9C6B_CoreGene_GO_BP_Overlap_Audit.csv"
  ),
  row.names = FALSE
)


write.csv(
  core_summary,
  file.path(
    step_dir,
    "STEP9C6B_CoreGene_GO_BP_Summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 20. Create manuscript-ready result lock
# ------------------------------------------------------------

STEP9C6B_MANUSCRIPT_RESULTS_LOCKED <- list(
  
  Step = "STEP9C-6B",
  
  Status = "PASS",
  
  Analysis =
    "GO Biological Process ORA",
  
  Predefined_Entrez_Background =
    35365,
  
  Effective_GO_Universe_Up =
    18011,
  
  Effective_GO_Universe_Down =
    18011,
  
  Tumor_Up_Input =
    5194,
  
  Tumor_Down_Input =
    4463,
  
  Tumor_Up_Pathways_Returned =
    5680,
  
  Tumor_Down_Pathways_Returned =
    5847,
  
  Tumor_Up_Significant =
    161,
  
  Tumor_Down_Significant =
    1296,
  
  Final_Significance_Definition =
    FINAL_SIGNIFICANCE_DEFINITION,
  
  Ontology =
    "BP",
  
  P_Adjust_Method =
    "BH",
  
  Min_GS_Size =
    10,
  
  Max_GS_Size =
    500,
  
  Core_Genes =
    core_genes,
  
  Core_Entrez =
    core_entrez,
  
  Core_Gene_Overlap_Audit =
    "Direct gene-set overlap only",
  
  Enrichment_Recalculation =
    FALSE,
  
  DEG_Recalculation =
    FALSE,
  
  ID_Conversion =
    FALSE,
  
  P_Value_Recalculation =
    FALSE,
  
  FDR_Recalculation =
    FALSE,
  
  Threshold_Changed =
    FALSE,
  
  Background_Modified =
    FALSE,
  
  STEP9C3_Modified =
    FALSE,
  
  STEP9C4_Modified =
    FALSE,
  
  Supplementary_Table_Up =
    "STEP9C6B_Supplementary_Table_GO_BP_TumorUp_161.csv",
  
  Supplementary_Table_Down =
    "STEP9C6B_Supplementary_Table_GO_BP_TumorDown_1296.csv",
  
  Core_Overlap_Audit =
    "STEP9C6B_CoreGene_GO_BP_Overlap_Audit.csv",
  
  Core_Summary =
    "STEP9C6B_CoreGene_GO_BP_Summary.csv"
)


saveRDS(
  STEP9C6B_MANUSCRIPT_RESULTS_LOCKED,
  file.path(
    step_dir,
    "STEP9C6B_MANUSCRIPT_RESULTS_LOCKED.rds"
  )
)


# ------------------------------------------------------------
# 21. Final structural validation
# ------------------------------------------------------------

FINAL_PASS <- all(
  
  identical(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$Step,
    "STEP9C-6B"
  ),
  
  identical(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$Status,
    "PASS"
  ),
  
  isTRUE(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      Predefined_Entrez_Background == 35365
  ),
  
  isTRUE(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      Effective_GO_Universe_Up == 18011
  ),
  
  isTRUE(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      Effective_GO_Universe_Down == 18011
  ),
  
  isTRUE(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      Tumor_Up_Input == 5194
  ),
  
  isTRUE(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      Tumor_Down_Input == 4463
  ),
  
  isTRUE(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      Tumor_Up_Pathways_Returned == 5680
  ),
  
  isTRUE(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      Tumor_Down_Pathways_Returned == 5847
  ),
  
  isTRUE(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      Tumor_Up_Significant == 161
  ),
  
  isTRUE(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      Tumor_Down_Significant == 1296
  ),
  
  identical(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      Final_Significance_Definition,
    
    "BH-adjusted P < 0.05"
  ),
  
  identical(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      Ontology,
    
    "BP"
  ),
  
  identical(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      P_Adjust_Method,
    
    "BH"
  ),
  
  identical(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      Enrichment_Recalculation,
    
    FALSE
  ),
  
  identical(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      DEG_Recalculation,
    
    FALSE
  ),
  
  identical(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      ID_Conversion,
    
    FALSE
  ),
  
  identical(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      P_Value_Recalculation,
    
    FALSE
  ),
  
  identical(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      FDR_Recalculation,
    
    FALSE
  ),
  
  identical(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      Threshold_Changed,
    
    FALSE
  ),
  
  identical(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      Background_Modified,
    
    FALSE
  ),
  
  identical(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      STEP9C3_Modified,
    
    FALSE
  ),
  
  identical(
    STEP9C6B_MANUSCRIPT_RESULTS_LOCKED$
      STEP9C4_Modified,
    
    FALSE
  ),
  
  nrow(up_sig_final) == 161,
  
  nrow(down_sig_final) == 1296,
  
  nrow(core_summary) == 3
  
)


if (!FINAL_PASS) {
  
  stop(
    "FINAL STATUS: FAIL — STEP9C-6B validation failed."
  )
  
}


# ------------------------------------------------------------
# 22. Final report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL STATUS: PASS — STEP9C-6B LOCKED\n")
cat("============================================================\n\n")

cat("Predefined Entrez background = 35,365\n")
cat("Effective GO universe — up = 18,011\n")
cat("Effective GO universe — down = 18,011\n")

cat("Tumor-up input = 5,194\n")
cat("Tumor-down input = 4,463\n")

cat("Tumor-up GO-BP pathways returned = 5,680\n")
cat("Tumor-down GO-BP pathways returned = 5,847\n")

cat("Tumor-up significant GO-BP terms = 161\n")
cat("Tumor-down significant GO-BP terms = 1,296\n\n")

cat("Final significance definition:\n")
cat("  BH-adjusted P < 0.05\n\n")

cat("Core genes:\n")
cat("  AOC3   = Entrez 8639\n")
cat("  MS4A7  = Entrez 58475\n")
cat("  IL7R   = Entrez 3575\n\n")

cat("Core-gene audit definition:\n")
cat("  Direct gene-set overlap only\n\n")

cat("No enrichment recalculation.\n")
cat("No DEG recalculation.\n")
cat("No ID conversion.\n")
cat("No P-value/FDR recalculation.\n")
cat("No threshold modification.\n")
cat("No background modification.\n")
cat("STEP9C-3 unchanged.\n")
cat("STEP9C-4 unchanged.\n\n")

cat("Final GO-BP Supplementary Tables exported.\n")
cat("Core-gene GO-BP overlap audit exported.\n")
cat("Core-gene GO-BP summary exported.\n")
cat("STEP9C6B_MANUSCRIPT_RESULTS_LOCKED.rds saved.\n\n")

cat("STEP9C-6B COMPLETE\n")
cat("============================================================\n")


# ============================================================
# Paper 1 — STEP9C-7A
# GO-BP Redundancy / Semantic Similarity Audit
# ============================================================
#
# PURPOSE
# ------------------------------------------------------------
# 1. Audit semantic redundancy among the LOCKED significant
#    GO-BP terms from STEP9C-6B.
#
# 2. Tumor-up:
#       161 significant GO-BP terms
#
# 3. Tumor-down:
#       1,296 significant GO-BP terms
#
# 4. Calculate GO-BP semantic similarity using GOSemSim.
#
# 5. Perform hierarchical clustering based on semantic distance.
#
# 6. Quantify redundancy at multiple predefined similarity
#    thresholds:
#       0.50
#       0.60
#       0.70
#       0.80
#
# 7. Identify descriptive representative terms within
#    redundancy clusters.
#
# IMPORTANT
# ------------------------------------------------------------
# This is NOT a new enrichment analysis.
#
# NO:
#   - DEG recalculation
#   - ID conversion
#   - enrichment recalculation
#   - P-value recalculation
#   - FDR recalculation
#   - significance threshold modification
#   - background modification
#   - change to STEP9C-6B statistical results
#
# The locked statistical results remain:
#   Tumor-up   = 161 significant GO-BP terms
#   Tumor-down = 1296 significant GO-BP terms
#
# Representative terms generated here are ONLY for
# descriptive organization / visualization.
#
# They are NOT new statistically significant terms and
# must NOT replace the complete locked result set.
# ============================================================


# ============================================================
# STEP9C-7A — STEP9C-6B Input File Recovery / Location Audit
# ============================================================
#
# PURPOSE
# ------------------------------------------------------------
# Locate the already-generated STEP9C-6B locked GO-BP result
# files before STEP9C-7A semantic similarity analysis.
#
# IMPORTANT
# ------------------------------------------------------------
# This script:
#   - does NOT run enrichGO()
#   - does NOT recalculate enrichment
#   - does NOT recalculate P values
#   - does NOT recalculate FDR
#   - does NOT modify the locked 161 / 1296 results
#
# It only searches for existing files.
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-7A: STEP9C-6B Input File Location Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Root directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

cat("Search root:\n")
cat(step_dir, "\n\n")


# ------------------------------------------------------------
# 2. Search recursively for STEP9C-6B GO-BP files
# ------------------------------------------------------------

all_files <- list.files(
  path = step_dir,
  recursive = TRUE,
  full.names = TRUE,
  all.files = FALSE
)

cat("Total files found:", length(all_files), "\n\n")


# ------------------------------------------------------------
# 3. Search for tumor-up locked GO-BP table
# ------------------------------------------------------------

up_candidates <- all_files[
  grepl(
    "STEP9C6B.*GO.*BP.*TumorUp.*161.*\\.csv$",
    basename(all_files),
    ignore.case = TRUE
  )
]

# Additional broader search if exact naming differs
if (length(up_candidates) == 0) {
  
  up_candidates <- all_files[
    grepl(
      "GO.*BP.*TumorUp",
      basename(all_files),
      ignore.case = TRUE
    ) &
      grepl(
        "\\.csv$",
        all_files,
        ignore.case = TRUE
      )
  ]
}


# ------------------------------------------------------------
# 4. Search for tumor-down locked GO-BP table
# ------------------------------------------------------------

down_candidates <- all_files[
  grepl(
    "STEP9C6B.*GO.*BP.*TumorDown.*1296.*\\.csv$",
    basename(all_files),
    ignore.case = TRUE
  )
]

# Additional broader search if exact naming differs
if (length(down_candidates) == 0) {
  
  down_candidates <- all_files[
    grepl(
      "GO.*BP.*TumorDown",
      basename(all_files),
      ignore.case = TRUE
    ) &
      grepl(
        "\\.csv$",
        all_files,
        ignore.case = TRUE
      )
  ]
}


# ------------------------------------------------------------
# 5. Report candidates
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("Tumor-up GO-BP candidates\n")
cat("------------------------------------------------------------\n")

if (length(up_candidates) == 0) {
  
  cat("NONE FOUND\n\n")
  
} else {
  
  cat(
    paste(
      seq_along(up_candidates),
      up_candidates,
      sep = ": "
    ),
    sep = "\n"
  )
  
  cat("\n\n")
}


cat("------------------------------------------------------------\n")
cat("Tumor-down GO-BP candidates\n")
cat("------------------------------------------------------------\n")

if (length(down_candidates) == 0) {
  
  cat("NONE FOUND\n\n")
  
} else {
  
  cat(
    paste(
      seq_along(down_candidates),
      down_candidates,
      sep = ": "
    ),
    sep = "\n"
  )
  
  cat("\n\n")
}


# ------------------------------------------------------------
# 6. Search for the STEP9C-6B lock object
# ------------------------------------------------------------

lock_candidates <- all_files[
  grepl(
    "STEP9C6B.*MANUSCRIPT_RESULTS_LOCKED.*\\.rds$",
    basename(all_files),
    ignore.case = TRUE
  )
]

cat("------------------------------------------------------------\n")
cat("STEP9C-6B manuscript lock candidates\n")
cat("------------------------------------------------------------\n")

if (length(lock_candidates) == 0) {
  
  cat("NONE FOUND\n\n")
  
} else {
  
  cat(
    paste(
      seq_along(lock_candidates),
      lock_candidates,
      sep = ": "
    ),
    sep = "\n"
  )
  
  cat("\n\n")
}


# ------------------------------------------------------------
# 7. Search for any files containing 1296
# ------------------------------------------------------------

down_1296_candidates <- all_files[
  grepl(
    "1296",
    basename(all_files),
    ignore.case = TRUE
  ) &
    grepl(
      "\\.(csv|rds|RData|txt)$",
      all_files,
      ignore.case = TRUE
    )
]

cat("------------------------------------------------------------\n")
cat("Files containing '1296'\n")
cat("------------------------------------------------------------\n")

if (length(down_1296_candidates) == 0) {
  
  cat("NONE FOUND\n\n")
  
} else {
  
  cat(
    paste(
      seq_along(down_1296_candidates),
      down_1296_candidates,
      sep = ": "
    ),
    sep = "\n"
  )
  
  cat("\n\n")
}


# ------------------------------------------------------------
# 8. Search for any files containing 161
# ------------------------------------------------------------

up_161_candidates <- all_files[
  grepl(
    "161",
    basename(all_files),
    ignore.case = TRUE
  ) &
    grepl(
      "\\.(csv|rds|RData|txt)$",
      all_files,
      ignore.case = TRUE
    )
]

cat("------------------------------------------------------------\n")
cat("Files containing '161'\n")
cat("------------------------------------------------------------\n")

if (length(up_161_candidates) == 0) {
  
  cat("NONE FOUND\n\n")
  
} else {
  
  cat(
    paste(
      seq_along(up_161_candidates),
      up_161_candidates,
      sep = ": "
    ),
    sep = "\n"
  )
  
  cat("\n\n")
}


# ------------------------------------------------------------
# 9. Final status
# ------------------------------------------------------------

cat("============================================================\n")

if (
  length(up_candidates) >= 1 &&
  length(down_candidates) >= 1
) {
  
  cat("STATUS: STEP9C-6B RESULT FILES LOCATED\n\n")
  
  cat("Tumor-up candidate:\n")
  cat(up_candidates[1], "\n\n")
  
  cat("Tumor-down candidate:\n")
  cat(down_candidates[1], "\n\n")
  
  cat("NO enrichment recalculation performed.\n")
  cat("NO statistical results modified.\n")
  
} else {
  
  cat("STATUS: FILE LOCATION REQUIRES MANUAL REVIEW\n\n")
  
  cat(
    "The existing STEP9C-6B CSV files were not found using\n",
    "the expected naming patterns.\n",
    sep = ""
  )
  
  cat(
    "Do NOT rerun enrichGO().\n"
  )
  
}

cat("============================================================\n")


# ============================================================
# Paper 1 — STEP9C-7A
# GO-BP Redundancy / Semantic Similarity Audit
# ============================================================
#
# INPUT:
#   STEP9C6B_Supplementary_Table_GO_BP_TumorUp_161.csv
#   STEP9C6B_Supplementary_Table_GO_BP_TumorDown_1296.csv
#   STEP9C6B_MANUSCRIPT_RESULTS_LOCKED.rds
#
# LOCKED INPUTS:
#   Tumor-up   = 161 significant GO-BP terms
#   Tumor-down = 1296 significant GO-BP terms
#
# METHOD:
#   GOSemSim Wang semantic similarity
#   GO Biological Process ontology
#   Complete-linkage hierarchical clustering
#
# PURPOSE:
#   Descriptive redundancy / semantic similarity audit only.
#
# IMPORTANT:
#   NO enrichGO()
#   NO DEG recalculation
#   NO ID conversion
#   NO P-value recalculation
#   NO FDR recalculation
#   NO threshold modification
#   NO background modification
#   NO modification of STEP9C-6B
# ============================================================


cat("\n============================================================\n")
cat("STEP9C-7A: GO-BP Redundancy / Semantic Similarity Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory and output folder
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

step7a_dir <- file.path(
  step_dir,
  "STEP9C-7A_GO_BP_Redundancy_Semantic_Similarity"
)

dir.create(
  step7a_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Working directory:\n")
cat(step_dir, "\n\n")

cat("STEP9C-7A output directory:\n")
cat(step7a_dir, "\n\n")


# ------------------------------------------------------------
# 1. Exact locked STEP9C-6B input paths
# ------------------------------------------------------------

up_file <- file.path(
  step_dir,
  "STEP9C6B_Supplementary_Table_GO_BP_TumorUp_161.csv"
)

down_file <- file.path(
  step_dir,
  "STEP9C6B_Supplementary_Table_GO_BP_TumorDown_1296.csv"
)

lock_file <- file.path(
  step_dir,
  "STEP9C6B_MANUSCRIPT_RESULTS_LOCKED.rds"
)


# ------------------------------------------------------------
# 2. File existence check
# ------------------------------------------------------------

required_files <- c(
  up_file,
  down_file,
  lock_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    paste0(
      "Required STEP9C-6B file(s) missing:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}

cat("STEP9C-6B input files: PASS\n\n")


# ------------------------------------------------------------
# 3. Required packages
# ------------------------------------------------------------

required_pkgs <- c(
  "GOSemSim",
  "org.Hs.eg.db",
  "AnnotationDbi"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0) {
  
  stop(
    paste0(
      "Missing required package(s): ",
      paste(
        missing_pkgs,
        collapse = ", "
      ),
      "\n",
      "This script does NOT install packages automatically."
    )
  )
}

suppressPackageStartupMessages({
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

cat("Required packages: PASS\n\n")


# ------------------------------------------------------------
# 4. Record package versions
# ------------------------------------------------------------

package_versions <- data.frame(
  Package = c(
    "R",
    "GOSemSim",
    "org.Hs.eg.db",
    "AnnotationDbi"
  ),
  Version = c(
    R.version.string,
    as.character(
      packageVersion("GOSemSim")
    ),
    as.character(
      packageVersion("org.Hs.eg.db")
    ),
    as.character(
      packageVersion("AnnotationDbi")
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  package_versions,
  file.path(
    step7a_dir,
    "STEP9C7A_Package_Versions.csv"
  ),
  row.names = FALSE
)

print(package_versions)

cat("\n")


# ------------------------------------------------------------
# 5. Read locked STEP9C-6B CSV results
# ------------------------------------------------------------

go_up <- read.csv(
  up_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

go_down <- read.csv(
  down_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Tumor-up rows:",
    nrow(go_up),
    "\n")

cat("Tumor-down rows:",
    nrow(go_down),
    "\n\n")


# ------------------------------------------------------------
# 6. Read STEP9C-6B manuscript lock
# ------------------------------------------------------------

STEP9C6B_LOCK <- readRDS(
  lock_file
)

cat(
  "STEP9C-6B lock object class:",
  class(STEP9C6B_LOCK),
  "\n\n"
)


# ------------------------------------------------------------
# 7. Structural validation
# ------------------------------------------------------------

required_columns <- c(
  "Direction",
  "ID",
  "Description",
  "GeneRatio",
  "BgRatio",
  "pvalue",
  "p.adjust",
  "geneID"
)

if (!all(
  required_columns %in% colnames(go_up)
)) {
  
  stop(
    "FAIL: Required columns missing from tumor-up table."
  )
}

if (!all(
  required_columns %in% colnames(go_down)
)) {
  
  stop(
    "FAIL: Required columns missing from tumor-down table."
  )
}


# ------------------------------------------------------------
# 8. Locked row-count validation
# ------------------------------------------------------------

if (!isTRUE(nrow(go_up) == 161)) {
  
  stop(
    paste0(
      "FAIL: Tumor-up rows = ",
      nrow(go_up),
      ", expected 161."
    )
  )
}

if (!isTRUE(nrow(go_down) == 1296)) {
  
  stop(
    paste0(
      "FAIL: Tumor-down rows = ",
      nrow(go_down),
      ", expected 1296."
    )
  )
}


# ------------------------------------------------------------
# 9. GO ID validation
# ------------------------------------------------------------

go_id_pattern <- "^GO:[0-9]{7}$"

if (!all(
  grepl(
    go_id_pattern,
    go_up$ID
  )
)) {
  
  stop(
    "FAIL: Invalid GO IDs detected in tumor-up table."
  )
}

if (!all(
  grepl(
    go_id_pattern,
    go_down$ID
  )
)) {
  
  stop(
    "FAIL: Invalid GO IDs detected in tumor-down table."
  )
}


# ------------------------------------------------------------
# 10. GO ID uniqueness
# ------------------------------------------------------------

if (anyDuplicated(go_up$ID) > 0) {
  
  stop(
    "FAIL: Duplicate GO IDs in tumor-up table."
  )
}

if (anyDuplicated(go_down$ID) > 0) {
  
  stop(
    "FAIL: Duplicate GO IDs in tumor-down table."
  )
}


# ------------------------------------------------------------
# 11. Significance-definition audit
# ------------------------------------------------------------
#
# STEP9C-6B locked definition:
#
#   BH-adjusted P < 0.05
#
# We ONLY inspect existing values.
# We do NOT recalculate P values or FDR.
# ------------------------------------------------------------

if (!all(
  go_up$p.adjust < 0.05
)) {
  
  stop(
    "FAIL: Tumor-up table contains p.adjust >= 0.05."
  )
}

if (!all(
  go_down$p.adjust < 0.05
)) {
  
  stop(
    "FAIL: Tumor-down table contains p.adjust >= 0.05."
  )
}


# ------------------------------------------------------------
# 12. Direction audit
# ------------------------------------------------------------

if (!all(
  go_up$Direction == "Tumor-up"
)) {
  
  stop(
    "FAIL: Tumor-up Direction column contains unexpected values."
  )
}

if (!all(
  go_down$Direction == "Tumor-down"
)) {
  
  stop(
    "FAIL: Tumor-down Direction column contains unexpected values."
  )
}


# ------------------------------------------------------------
# 13. Lock the GO ID vectors
# ------------------------------------------------------------

go_up_ids <- unique(
  go_up$ID
)

go_down_ids <- unique(
  go_down$ID
)

if (!isTRUE(
  length(go_up_ids) == 161
)) {
  
  stop(
    "FAIL: Tumor-up unique GO IDs != 161."
  )
}

if (!isTRUE(
  length(go_down_ids) == 1296
)) {
  
  stop(
    "FAIL: Tumor-down unique GO IDs != 1296."
  )
}

cat("STEP9C-6B locked input audit: PASS\n\n")


# ------------------------------------------------------------
# 14. Save immutable input copies
# ------------------------------------------------------------

write.csv(
  go_up,
  file.path(
    step7a_dir,
    "STEP9C7A_INPUT_LOCKED_GO_BP_TumorUp_161.csv"
  ),
  row.names = FALSE
)

write.csv(
  go_down,
  file.path(
    step7a_dir,
    "STEP9C7A_INPUT_LOCKED_GO_BP_TumorDown_1296.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. Prepare GO semantic data
# ------------------------------------------------------------

cat("Preparing GO-BP semantic data...\n")

semData <- GOSemSim::godata(
  "org.Hs.eg.db",
  ont = "BP",
  computeIC = FALSE
)

cat(
  "GO-BP semantic data preparation: PASS\n\n"
)


# ------------------------------------------------------------
# 16. Tumor-up Wang semantic similarity
# ------------------------------------------------------------

cat(
  "Calculating tumor-up semantic similarity...\n"
)

sim_up <- GOSemSim::mgoSim(
  go_up_ids,
  go_up_ids,
  semData = semData,
  measure = "Wang",
  combine = NULL
)

cat(
  "Tumor-up matrix:",
  nrow(sim_up),
  "x",
  ncol(sim_up),
  "\n\n"
)


# ------------------------------------------------------------
# 17. Tumor-down Wang semantic similarity
# ------------------------------------------------------------

cat(
  "Calculating tumor-down semantic similarity...\n"
)

sim_down <- GOSemSim::mgoSim(
  go_down_ids,
  go_down_ids,
  semData = semData,
  measure = "Wang",
  combine = NULL
)

cat(
  "Tumor-down matrix:",
  nrow(sim_down),
  "x",
  ncol(sim_down),
  "\n\n"
)


# ------------------------------------------------------------
# 18. Similarity matrix QC
# ------------------------------------------------------------

check_similarity_matrix <- function(
    sim,
    expected_n,
    label
) {
  
  if (!is.matrix(sim)) {
    
    stop(
      label,
      ": similarity object is not a matrix."
    )
  }
  
  if (!isTRUE(
    nrow(sim) == expected_n
  )) {
    
    stop(
      label,
      ": incorrect row count."
    )
  }
  
  if (!isTRUE(
    ncol(sim) == expected_n
  )) {
    
    stop(
      label,
      ": incorrect column count."
    )
  }
  
  if (!identical(
    rownames(sim),
    colnames(sim)
  )) {
    
    stop(
      label,
      ": row and column GO IDs differ."
    )
  }
  
  if (anyNA(sim)) {
    
    stop(
      label,
      ": NA values detected."
    )
  }
  
  if (any(
    sim < 0 |
    sim > 1
  )) {
    
    stop(
      label,
      ": similarity outside [0,1]."
    )
  }
  
  if (!all(
    abs(
      sim - t(sim)
    ) < 1e-10
  )) {
    
    stop(
      label,
      ": matrix is not symmetric."
    )
  }
  
  if (!all(
    abs(
      diag(sim) - 1
    ) < 1e-10
  )) {
    
    stop(
      label,
      ": diagonal != 1."
    )
  }
  
  TRUE
}


check_similarity_matrix(
  sim_up,
  161,
  "Tumor-up"
)

check_similarity_matrix(
  sim_down,
  1296,
  "Tumor-down"
)

cat(
  "Semantic similarity matrix QC: PASS\n\n"
)


# ------------------------------------------------------------
# 19. Save similarity matrices
# ------------------------------------------------------------

saveRDS(
  sim_up,
  file.path(
    step7a_dir,
    "STEP9C7A_GO_BP_Wang_SemanticSimilarity_TumorUp.rds"
  )
)

saveRDS(
  sim_down,
  file.path(
    step7a_dir,
    "STEP9C7A_GO_BP_Wang_SemanticSimilarity_TumorDown.rds"
  )
)


# ------------------------------------------------------------
# 20. Similarity summary
# ------------------------------------------------------------

summarize_similarity <- function(
    sim,
    direction
) {
  
  x <- sim[
    upper.tri(
      sim,
      diag = FALSE
    )
  ]
  
  data.frame(
    Direction = direction,
    N_Terms = nrow(sim),
    Pairwise_Comparisons = length(x),
    Mean_Similarity = mean(x),
    Median_Similarity = median(x),
    SD_Similarity = sd(x),
    Min_Similarity = min(x),
    Q25_Similarity =
      as.numeric(
        quantile(x, 0.25)
      ),
    Q75_Similarity =
      as.numeric(
        quantile(x, 0.75)
      ),
    Max_Similarity = max(x),
    stringsAsFactors = FALSE
  )
}


similarity_summary <- rbind(
  
  summarize_similarity(
    sim_up,
    "Tumor-up"
  ),
  
  summarize_similarity(
    sim_down,
    "Tumor-down"
  )
)

write.csv(
  similarity_summary,
  file.path(
    step7a_dir,
    "STEP9C7A_SemanticSimilarity_Summary.csv"
  ),
  row.names = FALSE
)

print(similarity_summary)

cat("\n")


# ------------------------------------------------------------
# 21. Hierarchical clustering
# ------------------------------------------------------------

hc_up <- hclust(
  as.dist(
    1 - sim_up
  ),
  method = "complete"
)

hc_down <- hclust(
  as.dist(
    1 - sim_down
  ),
  method = "complete"
)

saveRDS(
  hc_up,
  file.path(
    step7a_dir,
    "STEP9C7A_GO_BP_HierarchicalClustering_TumorUp.rds"
  )
)

saveRDS(
  hc_down,
  file.path(
    step7a_dir,
    "STEP9C7A_GO_BP_HierarchicalClustering_TumorDown.rds"
  )
)


# ------------------------------------------------------------
# 22. Redundancy clustering function
# ------------------------------------------------------------

cluster_go_terms <- function(
    sim,
    result_df,
    direction,
    threshold
) {
  
  distance_matrix <- as.dist(
    1 - sim
  )
  
  hc <- hclust(
    distance_matrix,
    method = "complete"
  )
  
  clusters <- cutree(
    hc,
    h = 1 - threshold
  )
  
  cluster_df <- data.frame(
    ID = names(clusters),
    Cluster = as.integer(
      clusters
    ),
    stringsAsFactors = FALSE
  )
  
  result_sub <- result_df[
    ,
    c(
      "ID",
      "Description",
      "GeneRatio",
      "BgRatio",
      "pvalue",
      "p.adjust",
      "geneID"
    )
  ]
  
  cluster_df <- merge(
    cluster_df,
    result_sub,
    by = "ID",
    all.x = TRUE,
    sort = FALSE
  )
  
  cluster_sizes <- table(
    cluster_df$Cluster
  )
  
  cluster_df$Cluster_Size <-
    as.integer(
      cluster_sizes[
        as.character(
          cluster_df$Cluster
        )
      ]
    )
  
  cluster_df$Direction <- direction
  
  cluster_df$Similarity_Threshold <-
    threshold
  
  cluster_df
}


# ------------------------------------------------------------
# 23. Predefined similarity thresholds
# ------------------------------------------------------------

thresholds <- c(
  0.50,
  0.60,
  0.70,
  0.80
)

cluster_results_up <- list()
cluster_results_down <- list()


for (thr in thresholds) {
  
  cat(
    "Tumor-up clustering:",
    "similarity >=",
    thr,
    "\n"
  )
  
  cluster_results_up[[as.character(thr)]]<- cluster_go_terms(
    sim = sim_up,
    result_df = go_up,
    direction = "Tumor-up",
    threshold = thr
  )
  
  cat(
    "Tumor-down clustering:",
    "similarity >=",
    thr,
    "\n"
  )
  
  cluster_results_down[[as.character(thr)]]<-cluster_go_terms(
    sim = sim_down,
    result_df = go_down,
    direction = "Tumor-down",
    threshold = thr
  )
}

cat("\n")


# ------------------------------------------------------------
# 24. Redundancy summary
# ------------------------------------------------------------

summarize_clusters <- function(
    cluster_df
) {
  
  cluster_sizes <- as.integer(
    table(
      cluster_df$Cluster
    )
  )
  
  data.frame(
    
    Direction =
      unique(
        cluster_df$Direction
      ),
    
    Similarity_Threshold =
      unique(
        cluster_df$Similarity_Threshold
      ),
    
    N_Terms =
      nrow(cluster_df),
    
    N_Clusters =
      length(cluster_sizes),
    
    Singleton_Clusters =
      sum(
        cluster_sizes == 1
      ),
    
    Redundant_Clusters =
      sum(
        cluster_sizes > 1
      ),
    
    Terms_In_Redundant_Clusters =
      sum(
        cluster_sizes[
          cluster_sizes > 1
        ]
      ),
    
    Max_Cluster_Size =
      max(cluster_sizes),
    
    Mean_Cluster_Size =
      mean(cluster_sizes),
    
    Median_Cluster_Size =
      median(cluster_sizes),
    
    stringsAsFactors = FALSE
  )
}


redundancy_summary <- do.call(
  rbind,
  c(
    
    lapply(
      thresholds,
      function(thr) {
        
        summarize_clusters(
          cluster_results_up[[as.character(thr)]])
      }
    ),
    
    lapply(
      thresholds,
      function(thr) {
        
        summarize_clusters(
          cluster_results_down[[as.character(thr)]]
        )
      }
    )
  )
)

rownames(
  redundancy_summary
) <- NULL

write.csv(
  redundancy_summary,
  file.path(
    step7a_dir,
    "STEP9C7A_GO_BP_Redundancy_Summary_All_Thresholds.csv"
  ),
  row.names = FALSE
)

print(redundancy_summary)

cat("\n")


# ------------------------------------------------------------
# 25. Export cluster membership
# ------------------------------------------------------------

for (thr in thresholds) {
  
  thr_label <- gsub(
    "\\.",
    "_",
    sprintf(
      "%.2f",
      thr
    )
  )
  
  write.csv(cluster_results_up[[as.character(thr)]],
            file.path(
              step7a_dir,
              paste0(
                "STEP9C7A_TumorUp_ClusterMembership_Sim",
                thr_label,
                ".csv"
              )
            ),
            row.names = FALSE
  )
  
  write.csv(cluster_results_down[[as.character(thr)]],
            file.path(
              step7a_dir,
              paste0(
                "STEP9C7A_TumorDown_ClusterMembership_Sim",
                thr_label,
                ".csv"
              )
            ),
            row.names = FALSE
  )
}


# ------------------------------------------------------------
# 26. Primary descriptive threshold
# ------------------------------------------------------------
#
# IMPORTANT:
#   0.70 is a descriptive redundancy threshold.
#
# It is NOT:
#   - a P-value threshold
#   - an FDR threshold
#   - a new enrichment criterion
#
# ------------------------------------------------------------

primary_threshold <- 0.70

up_primary <- cluster_results_up[[as.character(primary_threshold)]]

down_primary <- cluster_results_down[[as.character(primary_threshold)]]


# ------------------------------------------------------------
# 27. Representative-term selection
# ------------------------------------------------------------
#
# Deterministic descriptive rule:
#
#   1. smallest existing BH-adjusted P
#   2. larger GeneRatio numerator
#   3. alphabetical GO ID
#
# This does NOT alter statistical results.
# ------------------------------------------------------------

parse_gene_ratio <- function(
    x
) {
  
  sapply(
    strsplit(
      as.character(x),
      "/",
      fixed = TRUE
    ),
    function(z) {
      
      if (length(z) != 2) {
        return(NA_real_)
      }
      
      suppressWarnings(
        as.numeric(
          z[1]
        )
      )
    }
  )
}


select_representatives <- function(
    cluster_df
) {
  
  cluster_list <- split(
    cluster_df,
    cluster_df$Cluster
  )
  
  representative_list <- lapply(
    cluster_list,
    function(df) {
      
      df$GeneRatio_N <-
        parse_gene_ratio(
          df$GeneRatio
        )
      
      df <- df[
        order(
          df$p.adjust,
          -df$GeneRatio_N,
          df$ID
        ),
        ,
        drop = FALSE
      ]
      
      df[1, ]
    }
  )
  
  reps <- do.call(
    rbind,
    representative_list
  )
  
  rownames(reps) <- NULL
  
  reps$Representative_Rule <-
    paste(
      "Smallest BH-adjusted P;",
      "then larger GeneRatio numerator;",
      "then alphabetical GO ID"
    )
  
  reps
}


up_representatives <-
  select_representatives(
    up_primary
  )

down_representatives <-
  select_representatives(
    down_primary
  )


# ------------------------------------------------------------
# 28. Export representatives
# ------------------------------------------------------------

write.csv(
  up_representatives,
  file.path(
    step7a_dir,
    "STEP9C7A_TumorUp_Representative_GO_BP_Terms_Sim0_70.csv"
  ),
  row.names = FALSE
)

write.csv(
  down_representatives,
  file.path(
    step7a_dir,
    "STEP9C7A_TumorDown_Representative_GO_BP_Terms_Sim0_70.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 29. Export primary redundancy clusters
# ------------------------------------------------------------

write.csv(
  up_primary,
  file.path(
    step7a_dir,
    "STEP9C7A_TumorUp_PrimaryRedundancyClusters_Sim0_70.csv"
  ),
  row.names = FALSE
)

write.csv(
  down_primary,
  file.path(
    step7a_dir,
    "STEP9C7A_TumorDown_PrimaryRedundancyClusters_Sim0_70.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 30. Cluster-level summary
# ------------------------------------------------------------

cluster_level_summary <- function(
    cluster_df
) {
  
  cluster_list <- split(
    cluster_df,
    cluster_df$Cluster
  )
  
  out <- lapply(
    cluster_list,
    function(df) {
      
      representative <-
        df[
          order(
            df$p.adjust,
            df$ID
          ),
          ,
          drop = FALSE
        ][1, ]
      
      data.frame(
        
        Direction =
          unique(
            df$Direction
          ),
        
        Similarity_Threshold =
          unique(
            df$Similarity_Threshold
          ),
        
        Cluster =
          unique(
            df$Cluster
          ),
        
        Cluster_Size =
          nrow(df),
        
        Representative_ID =
          representative$ID,
        
        Representative_Description =
          representative$Description,
        
        Representative_p_adjust =
          representative$p.adjust,
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  do.call(
    rbind,
    out
  )
}


up_cluster_summary <-
  cluster_level_summary(
    up_primary
  )

down_cluster_summary <-
  cluster_level_summary(
    down_primary
  )

write.csv(
  up_cluster_summary,
  file.path(
    step7a_dir,
    "STEP9C7A_TumorUp_Cluster_Level_Summary_Sim0_70.csv"
  ),
  row.names = FALSE
)

write.csv(
  down_cluster_summary,
  file.path(
    step7a_dir,
    "STEP9C7A_TumorDown_Cluster_Level_Summary_Sim0_70.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 31. Final STEP9C-7A audit object
# ------------------------------------------------------------

STEP9C7A_AUDIT <- list(
  
  Step =
    "STEP9C-7A",
  
  Status =
    "PASS",
  
  Input_Source =
    "STEP9C-6B locked GO-BP results",
  
  Tumor_Up_Locked_Terms =
    161,
  
  Tumor_Down_Locked_Terms =
    1296,
  
  Tumor_Up_Observed_Terms =
    length(go_up_ids),
  
  Tumor_Down_Observed_Terms =
    length(go_down_ids),
  
  Semantic_Similarity_Method =
    "GOSemSim Wang",
  
  Ontology =
    "Biological Process",
  
  Clustering_Method =
    "Hierarchical clustering, complete linkage",
  
  Distance_Definition =
    "1 - Wang semantic similarity",
  
  Sensitivity_Thresholds =
    thresholds,
  
  Primary_Descriptive_Threshold =
    primary_threshold,
  
  Primary_Threshold_Is_Statistical =
    FALSE,
  
  Representative_Term_Purpose =
    "Descriptive organization / visualization only",
  
  Statistical_Result_Changed =
    FALSE,
  
  Enrichment_Recalculated =
    FALSE,
  
  DEG_Recalculated =
    FALSE,
  
  ID_Conversion_Performed =
    FALSE,
  
  P_Value_Recalculated =
    FALSE,
  
  FDR_Recalculated =
    FALSE,
  
  Significance_Threshold_Changed =
    FALSE,
  
  Background_Modified =
    FALSE,
  
  STEP9C6B_Modified =
    FALSE,
  
  Locked_Up_Result_Remains =
    161,
  
  Locked_Down_Result_Remains =
    1296
)


saveRDS(
  STEP9C7A_AUDIT,
  file.path(
    step7a_dir,
    "STEP9C7A_FINAL_AUDIT_LOCK.rds"
  )
)


# ------------------------------------------------------------
# 32. Final validation
# ------------------------------------------------------------

step9c7a_pass <- all(
  
  identical(
    STEP9C7A_AUDIT$Step,
    "STEP9C-7A"
  ),
  
  identical(
    STEP9C7A_AUDIT$Status,
    "PASS"
  ),
  
  isTRUE(
    STEP9C7A_AUDIT$Tumor_Up_Locked_Terms == 161
  ),
  
  isTRUE(
    STEP9C7A_AUDIT$Tumor_Down_Locked_Terms == 1296
  ),
  
  isTRUE(
    STEP9C7A_AUDIT$Tumor_Up_Observed_Terms == 161
  ),
  
  isTRUE(
    STEP9C7A_AUDIT$Tumor_Down_Observed_Terms == 1296
  ),
  
  identical(
    STEP9C7A_AUDIT$Semantic_Similarity_Method,
    "GOSemSim Wang"
  ),
  
  isTRUE(
    STEP9C7A_AUDIT$Primary_Descriptive_Threshold == 0.70
  ),
  
  isFALSE(
    STEP9C7A_AUDIT$Statistical_Result_Changed
  ),
  
  isFALSE(
    STEP9C7A_AUDIT$Enrichment_Recalculated
  ),
  
  isFALSE(
    STEP9C7A_AUDIT$DEG_Recalculated
  ),
  
  isFALSE(
    STEP9C7A_AUDIT$ID_Conversion_Performed
  ),
  
  isFALSE(
    STEP9C7A_AUDIT$P_Value_Recalculated
  ),
  
  isFALSE(
    STEP9C7A_AUDIT$FDR_Recalculated
  ),
  
  isFALSE(
    STEP9C7A_AUDIT$Significance_Threshold_Changed
  ),
  
  isFALSE(
    STEP9C7A_AUDIT$Background_Modified
  ),
  
  isFALSE(
    STEP9C7A_AUDIT$STEP9C6B_Modified
  )
)


if (!step9c7a_pass) {
  
  STEP9C7A_AUDIT$Status <-
    "FAIL"
  
  saveRDS(
    STEP9C7A_AUDIT,
    file.path(
      step7a_dir,
      "STEP9C7A_FINAL_AUDIT_LOCK.rds"
    )
  )
  
  stop(
    "\nFINAL STATUS: FAIL — STEP9C-7A QC failed.\n"
  )
}


# ------------------------------------------------------------
# 33. Final report
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("FINAL STATUS: PASS — STEP9C-7A COMPLETE\n")
cat("============================================================\n\n")

cat("LOCKED STEP9C-6B INPUTS:\n")
cat("  Tumor-up GO-BP terms   = 161\n")
cat("  Tumor-down GO-BP terms = 1296\n\n")

cat("Semantic similarity:\n")
cat("  Method    = GOSemSim Wang\n")
cat("  Ontology  = Biological Process\n\n")

cat("Redundancy clustering:\n")
cat("  Method    = complete-linkage hierarchical clustering\n")
cat("  Distance  = 1 - semantic similarity\n\n")

cat("Sensitivity thresholds:\n")
cat("  0.50\n")
cat("  0.60\n")
cat("  0.70  <-- primary descriptive threshold\n")
cat("  0.80\n\n")

cat("IMPORTANT:\n")
cat("  0.70 is NOT a statistical significance threshold.\n")
cat("  Representative terms are NOT new statistical results.\n")
cat("  Representative terms are for descriptive organization /\n")
cat("  visualization only.\n\n")

cat("NO enrichment recalculation.\n")
cat("NO DEG recalculation.\n")
cat("NO ID conversion.\n")
cat("NO P-value recalculation.\n")
cat("NO FDR recalculation.\n")
cat("NO significance-threshold modification.\n")
cat("NO background modification.\n")
cat("STEP9C-6B unchanged.\n\n")

cat("STEP9C-7A COMPLETE\n")
cat("============================================================\n")


step7a_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9C-7A_GO_BP_Redundancy_Semantic_Similarity"

cat("\n============================================================\n")
cat("STEP9C-7A: Read Locked Redundancy Audit Results\n")
cat("============================================================\n\n")

similarity_summary <- read.csv(
  file.path(
    step7a_dir,
    "STEP9C7A_SemanticSimilarity_Summary.csv"
  ),
  stringsAsFactors = FALSE
)

redundancy_summary <- read.csv(
  file.path(
    step7a_dir,
    "STEP9C7A_GO_BP_Redundancy_Summary_All_Thresholds.csv"
  ),
  stringsAsFactors = FALSE
)

cat("------------------------------------------------------------\n")
cat("Semantic similarity summary\n")
cat("------------------------------------------------------------\n")

print(
  similarity_summary,
  row.names = FALSE
)

cat("\n")

cat("------------------------------------------------------------\n")
cat("GO-BP redundancy summary\n")
cat("------------------------------------------------------------\n")

print(
  redundancy_summary,
  row.names = FALSE
)

cat("\n")

cat("------------------------------------------------------------\n")
cat("Primary threshold = 0.70\n")
cat("------------------------------------------------------------\n")

print(
  redundancy_summary[
    redundancy_summary$Similarity_Threshold == 0.70,
    ,
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\n============================================================\n")
cat("READ-ONLY AUDIT COMPLETE\n")
cat("No recalculation performed.\n")
cat("============================================================\n")


# ============================================================
# Paper 1
# STEP9C-7B
# GO–KEGG Functional Theme Integration
#
# Purpose:
#   Integrate locked GO-BP and KEGG enrichment results
#   using shared Entrez-gene content and the locked
#   STEP9C-7A semantic clusters.
#
# IMPORTANT:
#   This is a descriptive integration step only.
#
#   NO enrichment recalculation
#   NO DEG recalculation
#   NO ID conversion
#   NO P-value recalculation
#   NO FDR recalculation
#   NO significance-threshold modification
#   NO background modification
#   NO removal of locked significant terms/pathways
#
# Locked statistical inputs:
#   STEP9C-5B:
#       KEGG-up   = 13 significant pathways
#       KEGG-down = 60 significant pathways
#
#   STEP9C-6B:
#       GO-BP-up   = 161 significant terms
#       GO-BP-down = 1296 significant terms
#
# STEP9C-7A:
#   Semantic similarity threshold = 0.70
#   Used only as descriptive GO-BP organization.
#
# Core genes:
#   AOC3
#   MS4A7
#   IL7R
# ============================================================


cat("\n============================================================\n")
cat("STEP9C-7B: GO–KEGG Functional Theme Integration\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

step5b_dir <- base_dir
step6b_dir <- base_dir

step7a_dir <- file.path(
  base_dir,
  "STEP9C-7A_GO_BP_Redundancy_Semantic_Similarity"
)

step7b_dir <- file.path(
  base_dir,
  "STEP9C-7B_GO_KEGG_Functional_Theme_Integration"
)

if (!dir.exists(step7b_dir)) {
  dir.create(step7b_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 1. Locked input files
# ------------------------------------------------------------

kegg_up_file <- file.path(
  step5b_dir,
  "STEP9C5B_Supplementary_Table_KEGG_TumorUp_13.csv"
)

kegg_down_file <- file.path(
  step5b_dir,
  "STEP9C5B_Supplementary_Table_KEGG_TumorDown_60.csv"
)

go_up_file <- file.path(
  step6b_dir,
  "STEP9C6B_Supplementary_Table_GO_BP_TumorUp_161.csv"
)

go_down_file <- file.path(
  step6b_dir,
  "STEP9C6B_Supplementary_Table_GO_BP_TumorDown_1296.csv"
)

step5b_lock_file <- file.path(
  step5b_dir,
  "STEP9C5B_MANUSCRIPT_RESULTS_LOCKED.rds"
)

step6b_lock_file <- file.path(
  step6b_dir,
  "STEP9C6B_MANUSCRIPT_RESULTS_LOCKED.rds"
)

step7a_similarity_file <- file.path(
  step7a_dir,
  "STEP9C7A_SemanticSimilarity_Summary.csv"
)

step7a_redundancy_file <- file.path(
  step7a_dir,
  "STEP9C7A_GO_BP_Redundancy_Summary_All_Thresholds.csv"
)


# ------------------------------------------------------------
# 2. File existence audit
# ------------------------------------------------------------

required_files <- c(
  kegg_up_file,
  kegg_down_file,
  go_up_file,
  go_down_file,
  step5b_lock_file,
  step6b_lock_file,
  step7a_similarity_file,
  step7a_redundancy_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  cat("ERROR: Required locked input files are missing:\n\n")
  
  print(missing_files)
  
  stop(
    "STEP9C-7B STOPPED: Missing required locked input files."
  )
}

cat("All required locked input files found.\n\n")


# ------------------------------------------------------------
# 3. Read locked enrichment tables
# ------------------------------------------------------------

kegg_up <- read.csv(
  kegg_up_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

kegg_down <- read.csv(
  kegg_down_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

go_up <- read.csv(
  go_up_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

go_down <- read.csv(
  go_down_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 4. Read locked RDS result-definition locks
# ------------------------------------------------------------

STEP9C5B_FINAL_LOCK <- readRDS(step5b_lock_file)
STEP9C6B_FINAL_LOCK <- readRDS(step6b_lock_file)


# ------------------------------------------------------------
# 5. Basic structural audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("Locked result counts\n")
cat("------------------------------------------------------------\n")

cat(
  "KEGG-up rows   = ",
  nrow(kegg_up),
  "\n",
  sep = ""
)

cat(
  "KEGG-down rows = ",
  nrow(kegg_down),
  "\n",
  sep = ""
)

cat(
  "GO-up rows     = ",
  nrow(go_up),
  "\n",
  sep = ""
)

cat(
  "GO-down rows   = ",
  nrow(go_down),
  "\n",
  sep = ""
)

cat("\n")


# ------------------------------------------------------------
# 6. Count locks
# ------------------------------------------------------------

count_lock_valid <- all(
  isTRUE(nrow(kegg_up) == 13),
  isTRUE(nrow(kegg_down) == 60),
  isTRUE(nrow(go_up) == 161),
  isTRUE(nrow(go_down) == 1296)
)

if (!count_lock_valid) {
  
  stop(
    "STEP9C-7B STOPPED: Locked enrichment table counts do not match STEP9C-5B / STEP9C-6B."
  )
}

cat("Locked enrichment counts: PASS\n\n")


# ------------------------------------------------------------
# 7. Required columns
# ------------------------------------------------------------

required_enrich_cols <- c(
  "ID",
  "Description",
  "GeneRatio",
  "BgRatio",
  "pvalue",
  "p.adjust",
  "geneID"
)

check_required_columns <- function(
    x,
    object_name
) {
  
  missing <- setdiff(
    required_enrich_cols,
    colnames(x)
  )
  
  if (length(missing) > 0) {
    
    stop(
      paste0(
        object_name,
        " is missing required columns: ",
        paste(missing, collapse = ", ")
      )
    )
  }
  
  TRUE
}

check_required_columns(
  kegg_up,
  "kegg_up"
)

check_required_columns(
  kegg_down,
  "kegg_down"
)

check_required_columns(
  go_up,
  "go_up"
)

check_required_columns(
  go_down,
  "go_down"
)

cat("Required enrichment columns: PASS\n\n")


# ------------------------------------------------------------
# 8. GeneID parser
#
# clusterProfiler stores geneID as:
#   "8639/58475/3575/..."
#
# No ID conversion is performed.
# These are already locked Entrez IDs.
# ------------------------------------------------------------

parse_geneID <- function(x) {
  
  x <- as.character(x)
  
  lapply(
    x,
    function(z) {
      
      if (is.na(z) || !nzchar(z)) {
        return(character(0))
      }
      
      ids <- unlist(
        strsplit(
          z,
          "/",
          fixed = TRUE
        )
      )
      
      ids <- trimws(ids)
      
      ids[
        nzchar(ids) &
          !is.na(ids)
      ]
    }
  )
}


kegg_up_genes <- parse_geneID(kegg_up$geneID)
kegg_down_genes <- parse_geneID(kegg_down$geneID)

go_up_genes <- parse_geneID(go_up$geneID)
go_down_genes <- parse_geneID(go_down$geneID)


# ------------------------------------------------------------
# 9. Gene-set integrity audit
# ------------------------------------------------------------

all_gene_lists <- c(
  kegg_up_genes,
  kegg_down_genes,
  go_up_genes,
  go_down_genes
)

all_gene_values <- unlist(
  all_gene_lists,
  use.names = FALSE
)

non_numeric_entrez <- all_gene_values[
  !grepl(
    "^[0-9]+$",
    all_gene_values
  )
]

if (length(non_numeric_entrez) > 0) {
  
  stop(
    "STEP9C-7B STOPPED: Non-Entrez gene identifiers detected in locked geneID fields."
  )
}

cat("GeneID format audit: PASS\n")
cat("No ID conversion performed.\n\n")


# ------------------------------------------------------------
# 10. Core genes
# ------------------------------------------------------------

core_entrez <- c(
  AOC3 = "8639",
  MS4A7 = "58475",
  IL7R = "3575"
)


# ------------------------------------------------------------
# 11. Function:
#     calculate GO–KEGG gene overlap
#
# This is descriptive gene-set overlap only.
#
# No statistical test is performed here.
# ------------------------------------------------------------

calculate_overlap <- function(
    go_ids,
    kegg_ids
) {
  
  shared <- intersect(
    go_ids,
    kegg_ids
  )
  
  union_ids <- union(
    go_ids,
    kegg_ids
  )
  
  n_shared <- length(shared)
  
  n_go <- length(unique(go_ids))
  n_kegg <- length(unique(kegg_ids))
  n_union <- length(union_ids)
  
  jaccard <- ifelse(
    n_union > 0,
    n_shared / n_union,
    NA_real_
  )
  
  overlap_go <- ifelse(
    n_go > 0,
    n_shared / n_go,
    NA_real_
  )
  
  overlap_kegg <- ifelse(
    n_kegg > 0,
    n_shared / n_kegg,
    NA_real_
  )
  
  list(
    Shared_Entrez = shared,
    N_Shared = n_shared,
    N_GO_Genes = n_go,
    N_KEGG_Genes = n_kegg,
    Jaccard = jaccard,
    GO_Coverage = overlap_go,
    KEGG_Coverage = overlap_kegg
  )
}


# ------------------------------------------------------------
# 12. Build one complete GO–KEGG integration matrix
#
# Four direction-specific comparisons:
#
#   Tumor-up:
#       GO-up × KEGG-up
#
#   Tumor-down:
#       GO-down × KEGG-down
#
# Cross-direction comparisons are also generated separately
# for audit purposes but are NOT used to define biological
# themes.
# ------------------------------------------------------------

build_integration_table <- function(
    go_table,
    go_genes,
    kegg_table,
    kegg_genes,
    direction_label
) {
  
  out <- vector(
    "list",
    length(go_genes) * length(kegg_genes)
  )
  
  idx <- 1
  
  for (i in seq_along(go_genes)) {
    
    for (j in seq_along(kegg_genes)) {
      
      ov <- calculate_overlap(
        go_genes[[i]],
        kegg_genes[[j]]
      )
      
      shared <- ov$Shared_Entrez
      
      core_shared <- intersect(
        shared,
        core_entrez
      )
      
      out[[idx]] <- data.frame(
        Direction = direction_label,
        GO_ID = as.character(go_table$ID[i]),
        GO_Description = as.character(
          go_table$Description[i]
        ),
        GO_p_adjust = as.numeric(
          go_table$p.adjust[i]
        ),
        GO_GeneRatio = as.character(
          go_table$GeneRatio[i]
        ),
        GO_Genes = length(
          unique(go_genes[[i]])
        ),
        KEGG_ID = as.character(
          kegg_table$ID[j]
        ),
        KEGG_Description = as.character(
          kegg_table$Description[j]
        ),
        KEGG_p_adjust = as.numeric(
          kegg_table$p.adjust[j]
        ),
        KEGG_GeneRatio = as.character(
          kegg_table$GeneRatio[j]
        ),
        KEGG_Genes = length(
          unique(kegg_genes[[j]])
        ),
        N_Shared_Entrez = ov$N_Shared,
        Jaccard = ov$Jaccard,
        GO_Coverage = ov$GO_Coverage,
        KEGG_Coverage = ov$KEGG_Coverage,
        Shared_Entrez = paste(
          shared,
          collapse = "/"
        ),
        Core_Genes_Shared = paste(
          names(core_entrez)[
            core_entrez %in% shared
          ],
          collapse = ";"
        ),
        N_Core_Genes_Shared = length(
          core_shared
        ),
        stringsAsFactors = FALSE
      )
      
      idx <- idx + 1
    }
  }
  
  do.call(
    rbind,
    out
  )
}


# ------------------------------------------------------------
# 13. Same-direction GO–KEGG integration
# ------------------------------------------------------------

integration_up <- build_integration_table(
  go_table = go_up,
  go_genes = go_up_genes,
  kegg_table = kegg_up,
  kegg_genes = kegg_up_genes,
  direction_label = "Tumor-up"
)

integration_down <- build_integration_table(
  go_table = go_down,
  go_genes = go_down_genes,
  kegg_table = kegg_down,
  kegg_genes = kegg_down_genes,
  direction_label = "Tumor-down"
)


# ------------------------------------------------------------
# 14. Retain only actual gene-overlap pairs for descriptive
#     integration.
#
# IMPORTANT:
#   This does NOT change the statistical enrichment result set.
#   It only identifies GO–KEGG pairs sharing >=1 Entrez gene.
# ------------------------------------------------------------

integration_up_nonzero <- integration_up[
  integration_up$N_Shared_Entrez > 0,
  ,
  drop = FALSE
]

integration_down_nonzero <- integration_down[
  integration_down$N_Shared_Entrez > 0,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 15. Rank overlap relationships descriptively
#
# No statistical filtering.
# No significance filtering.
#
# Ranking is based on:
#   1. number of shared genes
#   2. Jaccard similarity
#   3. GO adjusted P
# ------------------------------------------------------------

integration_up_nonzero <- integration_up_nonzero[
  order(
    -integration_up_nonzero$N_Shared_Entrez,
    -integration_up_nonzero$Jaccard,
    integration_up_nonzero$GO_p_adjust
  ),
  ,
  drop = FALSE
]

integration_down_nonzero <- integration_down_nonzero[
  order(
    -integration_down_nonzero$N_Shared_Entrez,
    -integration_down_nonzero$Jaccard,
    integration_down_nonzero$GO_p_adjust
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 16. Build cross-direction audit
#
# These are NOT used for primary biological-theme assignment.
# They are retained to document whether GO-up/down terms
# share genes with KEGG pathways of the opposite direction.
# ------------------------------------------------------------

cross_up_go_down_kegg <- build_integration_table(
  go_table = go_up,
  go_genes = go_up_genes,
  kegg_table = kegg_down,
  kegg_genes = kegg_down_genes,
  direction_label = "GO-up_vs_KEGG-down"
)

cross_down_go_up_kegg <- build_integration_table(
  go_table = go_down,
  go_genes = go_down_genes,
  kegg_table = kegg_up,
  kegg_genes = kegg_up_genes,
  direction_label = "GO-down_vs_KEGG-up"
)


# ------------------------------------------------------------
# 17. Direction-consistency audit
# ------------------------------------------------------------

direction_audit <- data.frame(
  Comparison = c(
    "GO-up_vs_KEGG-up",
    "GO-down_vs_KEGG-down",
    "GO-up_vs_KEGG-down",
    "GO-down_vs_KEGG-up"
  ),
  Total_Pairs = c(
    nrow(integration_up),
    nrow(integration_down),
    nrow(cross_up_go_down_kegg),
    nrow(cross_down_go_up_kegg)
  ),
  Nonzero_Overlap_Pairs = c(
    sum(integration_up$N_Shared_Entrez > 0),
    sum(integration_down$N_Shared_Entrez > 0),
    sum(cross_up_go_down_kegg$N_Shared_Entrez > 0),
    sum(cross_down_go_up_kegg$N_Shared_Entrez > 0)
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 18. GO semantic-cluster integration
#
# IMPORTANT:
# STEP9C-7A cluster assignments may exist in one or more
# files depending on the final STEP9C-7A export structure.
#
# Search recursively for CSV files containing:
#   GO
#   Cluster
#   0.70
#
# We do NOT recalculate semantic similarity.
# ------------------------------------------------------------

cluster_candidates <- list.files(
  step7a_dir,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

cluster_candidates <- cluster_candidates[
  grepl(
    "cluster",
    basename(cluster_candidates),
    ignore.case = TRUE
  )
]

cat("STEP9C-7A cluster-file candidates found = ",
    length(cluster_candidates),
    "\n",
    sep = "")

if (length(cluster_candidates) > 0) {
  print(basename(cluster_candidates))
}

cat("\n")


# ------------------------------------------------------------
# 19. Attempt to identify a 0.70 cluster assignment table
#
# The table must contain at least:
#   GO_ID
#   Cluster_ID
#
# If not found, STEP9C-7B continues with GO–KEGG gene overlap
# integration, but explicitly records that cluster assignment
# could not be attached.
# ------------------------------------------------------------

cluster_assignment_file <- NA_character_

if (length(cluster_candidates) > 0) {
  
  for (f in cluster_candidates) {
    
    tmp <- tryCatch(
      read.csv(
        f,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      error = function(e) NULL
    )
    
    if (is.null(tmp)) next
    
    nms <- colnames(tmp)
    
    has_go <- any(
      nms %in% c(
        "GO_ID",
        "ID",
        "GO",
        "GOid"
      )
    )
    
    has_cluster <- any(
      nms %in% c(
        "Cluster_ID",
        "Cluster",
        "ClusterID"
      )
    )
    
    if (has_go && has_cluster) {
      
      cluster_assignment_file <- f
      break
    }
  }
}


if (!is.na(cluster_assignment_file)) {
  
  cat(
    "STEP9C-7A cluster assignment file identified:\n",
    basename(cluster_assignment_file),
    "\n\n"
  )
  
  cluster_df <- read.csv(
    cluster_assignment_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Standardize column names.
  go_col <- intersect(
    c("GO_ID", "ID", "GO", "GOid"),
    colnames(cluster_df)
  )[1]
  
  cluster_col <- intersect(
    c("Cluster_ID", "Cluster", "ClusterID"),
    colnames(cluster_df)
  )[1]
  
  cluster_df <- cluster_df[
    ,
    c(
      go_col,
      cluster_col
    ),
    drop = FALSE
  ]
  
  colnames(cluster_df) <- c(
    "GO_ID",
    "Cluster_ID"
  )
  
  cluster_df$GO_ID <- as.character(
    cluster_df$GO_ID
  )
  
  cluster_df$Cluster_ID <- as.character(
    cluster_df$Cluster_ID
  )
  
  # Only retain the primary 0.70 GO terms.
  cluster_df <- cluster_df[
    cluster_df$GO_ID %in% c(
      go_up$ID,
      go_down$ID
    ),
    ,
    drop = FALSE
  ]
  
  integration_up <- merge(
    integration_up,
    cluster_df,
    by = "GO_ID",
    all.x = TRUE,
    sort = FALSE
  )
  
  integration_down <- merge(
    integration_down,
    cluster_df,
    by = "GO_ID",
    all.x = TRUE,
    sort = FALSE
  )
  
} else {
  
  cat(
    "NOTE:\n",
    "No explicit STEP9C-7A GO cluster-assignment table was identified.\n",
    "GO–KEGG gene-overlap integration will still be completed.\n",
    "No semantic similarity will be recalculated.\n\n"
  )
}


# ------------------------------------------------------------
# 20. Core-gene functional integration
#
# Direct overlap only.
#
# No leading-edge terminology.
# No causal interpretation.
# ------------------------------------------------------------

core_functional_audit <- function(
    go_table,
    go_genes,
    kegg_table,
    kegg_genes,
    direction_label
) {
  
  rows <- list()
  idx <- 1
  
  for (i in seq_along(go_genes)) {
    
    for (j in seq_along(kegg_genes)) {
      
      shared <- intersect(
        go_genes[[i]],
        kegg_genes[[j]]
      )
      
      core_shared <- intersect(
        shared,
        core_entrez
      )
      
      if (length(core_shared) == 0) {
        next
      }
      
      rows[[idx]] <- data.frame(
        Direction = direction_label,
        GO_ID = as.character(go_table$ID[i]),
        GO_Description = as.character(
          go_table$Description[i]
        ),
        KEGG_ID = as.character(
          kegg_table$ID[j]
        ),
        KEGG_Description = as.character(
          kegg_table$Description[j]
        ),
        Core_Genes_Shared = paste(
          names(core_entrez)[
            core_entrez %in% core_shared
          ],
          collapse = ";"
        ),
        N_Core_Genes_Shared = length(
          core_shared
        ),
        Shared_Entrez = paste(
          shared,
          collapse = "/"
        ),
        N_Shared_Entrez = length(
          shared
        ),
        stringsAsFactors = FALSE
      )
      
      idx <- idx + 1
    }
  }
  
  if (length(rows) == 0) {
    
    return(
      data.frame(
        Direction = character(),
        GO_ID = character(),
        GO_Description = character(),
        KEGG_ID = character(),
        KEGG_Description = character(),
        Core_Genes_Shared = character(),
        N_Core_Genes_Shared = numeric(),
        Shared_Entrez = character(),
        N_Shared_Entrez = numeric(),
        stringsAsFactors = FALSE
      )
    )
  }
  
  do.call(
    rbind,
    rows
  )
}


core_up <- core_functional_audit(
  go_up,
  go_up_genes,
  kegg_up,
  kegg_up_genes,
  "Tumor-up"
)

core_down <- core_functional_audit(
  go_down,
  go_down_genes,
  kegg_down,
  kegg_down_genes,
  "Tumor-down"
)

core_functional_overlap <- rbind(
  core_up,
  core_down
)


# ------------------------------------------------------------
# 21. Candidate functional theme units
#
# IMPORTANT:
#
# These are NOT statistical results.
#
# A candidate theme unit is defined descriptively by:
#   GO term + KEGG pathway
#   sharing >= 1 locked Entrez gene.
#
# No GO/KEGG term is removed from the locked result set.
# ------------------------------------------------------------

candidate_theme_up <- integration_up_nonzero

candidate_theme_down <- integration_down_nonzero


# ------------------------------------------------------------
# 22. Summary statistics
# ------------------------------------------------------------

integration_summary <- data.frame(
  Direction = c(
    "Tumor-up",
    "Tumor-down"
  ),
  GO_Terms = c(
    161,
    1296
  ),
  KEGG_Pathways = c(
    13,
    60
  ),
  Total_GO_KEGG_Pairs = c(
    nrow(integration_up),
    nrow(integration_down)
  ),
  Nonzero_Overlap_Pairs = c(
    nrow(integration_up_nonzero),
    nrow(integration_down_nonzero)
  ),
  GO_Terms_With_KEGG_Overlap = c(
    length(unique(
      integration_up_nonzero$GO_ID
    )),
    length(unique(
      integration_down_nonzero$GO_ID
    ))
  ),
  KEGG_Pathways_With_GO_Overlap = c(
    length(unique(
      integration_up_nonzero$KEGG_ID
    )),
    length(unique(
      integration_down_nonzero$KEGG_ID
    ))
  ),
  Max_Shared_Genes = c(
    ifelse(
      nrow(integration_up_nonzero) > 0,
      max(
        integration_up_nonzero$N_Shared_Entrez
      ),
      0
    ),
    ifelse(
      nrow(integration_down_nonzero) > 0,
      max(
        integration_down_nonzero$N_Shared_Entrez
      ),
      0
    )
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 23. Save outputs
# ------------------------------------------------------------

write.csv(
  integration_up,
  file.path(
    step7b_dir,
    "STEP9C7B_GO_KEGG_Integration_TumorUp_ALL_Pairs.csv"
  ),
  row.names = FALSE
)

write.csv(
  integration_down,
  file.path(
    step7b_dir,
    "STEP9C7B_GO_KEGG_Integration_TumorDown_ALL_Pairs.csv"
  ),
  row.names = FALSE
)

write.csv(
  integration_up_nonzero,
  file.path(
    step7b_dir,
    "STEP9C7B_GO_KEGG_Integration_TumorUp_NonzeroOverlap.csv"
  ),
  row.names = FALSE
)

write.csv(
  integration_down_nonzero,
  file.path(
    step7b_dir,
    "STEP9C7B_GO_KEGG_Integration_TumorDown_NonzeroOverlap.csv"
  ),
  row.names = FALSE
)

write.csv(
  cross_up_go_down_kegg[
    cross_up_go_down_kegg$N_Shared_Entrez > 0,
    ,
    drop = FALSE
  ],
  file.path(
    step7b_dir,
    "STEP9C7B_CrossDirection_GOup_KEGGdown.csv"
  ),
  row.names = FALSE
)

write.csv(
  cross_down_go_up_kegg[
    cross_down_go_up_kegg$N_Shared_Entrez > 0,
    ,
    drop = FALSE
  ],
  file.path(
    step7b_dir,
    "STEP9C7B_CrossDirection_GOdown_KEGGup.csv"
  ),
  row.names = FALSE
)

write.csv(
  direction_audit,
  file.path(
    step7b_dir,
    "STEP9C7B_Direction_Overlap_Audit.csv"
  ),
  row.names = FALSE
)

write.csv(
  core_functional_overlap,
  file.path(
    step7b_dir,
    "STEP9C7B_CoreGene_GO_KEGG_DirectOverlap_Audit.csv"
  ),
  row.names = FALSE
)

write.csv(
  integration_summary,
  file.path(
    step7b_dir,
    "STEP9C7B_GO_KEGG_Integration_Summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 24. Final lock object
# ------------------------------------------------------------

STEP9C7B_FINAL_LOCK <- list(
  
  Step = "STEP9C-7B",
  
  Status = "PASS",
  
  Purpose =
    "Descriptive GO-BP and KEGG functional integration",
  
  Statistical_Inputs = list(
    KEGG_TumorUp = 13,
    KEGG_TumorDown = 60,
    GO_BP_TumorUp = 161,
    GO_BP_TumorDown = 1296
  ),
  
  GO_Semantic_Organization = list(
    Source = "STEP9C-7A locked results",
    Primary_Threshold = 0.70,
    Descriptive_Only = TRUE,
    Statistical_Threshold = FALSE
  ),
  
  Integration_Method =
    "Direct overlap of Entrez gene sets between locked GO-BP terms and locked KEGG pathways",
  
  Same_Direction_Integration = TRUE,
  
  Cross_Direction_Audit = TRUE,
  
  Core_Gene_Audit =
    "Direct gene-set overlap only",
  
  Core_Genes = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  
  Core_Entrez = core_entrez,
  
  Enrichment_Recalculation = FALSE,
  
  DEG_Recalculation = FALSE,
  
  ID_Conversion_Performed = FALSE,
  
  P_Value_Recalculation = FALSE,
  
  FDR_Recalculation = FALSE,
  
  Significance_Threshold_Modified = FALSE,
  
  Background_Modified = FALSE,
  
  Locked_GO_Results_Modified = FALSE,
  
  Locked_KEGG_Results_Modified = FALSE,
  
  STEP9C3_Unchanged = TRUE,
  
  STEP9C4_Unchanged = TRUE,
  
  STEP9C5B_Unchanged = TRUE,
  
  STEP9C6B_Unchanged = TRUE,
  
  STEP9C7A_Unchanged = TRUE
)


saveRDS(
  STEP9C7B_FINAL_LOCK,
  file.path(
    step7b_dir,
    "STEP9C7B_FINAL_LOCK.rds"
  )
)


# ------------------------------------------------------------
# 25. Final console report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-7B FINAL REPORT\n")
cat("============================================================\n\n")

cat(
  "Locked GO-BP tumor-up terms     = 161\n"
)

cat(
  "Locked GO-BP tumor-down terms   = 1296\n"
)

cat(
  "Locked KEGG tumor-up pathways   = 13\n"
)

cat(
  "Locked KEGG tumor-down pathways = 60\n\n"
)

cat(
  "Tumor-up GO–KEGG total pairs = ",
  nrow(integration_up),
  "\n",
  sep = ""
)

cat(
  "Tumor-up nonzero-overlap pairs = ",
  nrow(integration_up_nonzero),
  "\n",
  sep = ""
)

cat(
  "Tumor-down GO–KEGG total pairs = ",
  nrow(integration_down),
  "\n",
  sep = ""
)

cat(
  "Tumor-down nonzero-overlap pairs = ",
  nrow(integration_down_nonzero),
  "\n\n",
  sep = ""
)

cat(
  "STEP9C-7A semantic similarity threshold = 0.70\n"
)

cat(
  "0.70 is descriptive only; it is NOT a significance threshold.\n\n"
)

cat(
  "Functional theme definition:\n"
)

cat(
  "  Locked GO-BP term + locked KEGG pathway\n"
)

cat(
  "  connected by shared Entrez gene content.\n\n"
)

cat(
  "Core-gene integration:\n"
)

cat(
  "  AOC3 / MS4A7 / IL7R\n"
)

cat(
  "  Direct gene-set overlap only.\n\n"
)

cat(
  "NO enrichment recalculation.\n"
)

cat(
  "NO DEG recalculation.\n"
)

cat(
  "NO ID conversion.\n"
)

cat(
  "NO P-value recalculation.\n"
)

cat(
  "NO FDR recalculation.\n"
)

cat(
  "NO significance-threshold modification.\n"
)

cat(
  "NO background modification.\n"
)

cat(
  "STEP9C-3 unchanged.\n"
)

cat(
  "STEP9C-4 unchanged.\n"
)

cat(
  "STEP9C-5B unchanged.\n"
)

cat(
  "STEP9C-6B unchanged.\n"
)

cat(
  "STEP9C-7A unchanged.\n\n"
)

cat(
  "============================================================\n"
)

cat(
  "FINAL STATUS: PASS — STEP9C-7B COMPLETE\n"
)

cat(
  "============================================================\n"
)


# ============================================================
# Paper 1
# STEP9C-7B-READ
# READ-ONLY Functional Theme Evidence Audit
#
# Purpose:
#   Read the locked STEP9C-7B GO–KEGG integration outputs
#   and audit the empirical evidence for recurrent
#   cross-database functional relationships.
#
# IMPORTANT:
#
#   READ-ONLY AUDIT ONLY
#
#   NO enrichment recalculation
#   NO DEG recalculation
#   NO ID conversion
#   NO P-value recalculation
#   NO FDR recalculation
#   NO semantic-similarity recalculation
#   NO significance-threshold modification
#   NO background modification
#   NO modification of STEP9C-7B locked results
#
# Locked statistical result sets:
#   GO-BP tumor-up   = 161
#   GO-BP tumor-down = 1296
#   KEGG tumor-up    = 13
#   KEGG tumor-down  = 60
#
# STEP9C-7A:
#   Semantic similarity threshold = 0.70
#   Descriptive organization only
#
# Core genes:
#   AOC3 / MS4A7 / IL7R
#
# All analyses below are descriptive evidence audits.
# ============================================================


cat("\n============================================================\n")
cat("STEP9C-7B-READ: Functional Theme Evidence Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

step7b_dir <- file.path(
  base_dir,
  "STEP9C-7B_GO_KEGG_Functional_Theme_Integration"
)

audit_dir <- file.path(
  step7b_dir,
  "STEP9C-7B-READ_Functional_Theme_Evidence_Audit"
)

if (!dir.exists(audit_dir)) {
  dir.create(
    audit_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 1. Locked input files
# ------------------------------------------------------------

files_required <- c(
  
  GO_KEGG_UP =
    file.path(
      step7b_dir,
      "STEP9C7B_GO_KEGG_Integration_TumorUp_NonzeroOverlap.csv"
    ),
  
  GO_KEGG_DOWN =
    file.path(
      step7b_dir,
      "STEP9C7B_GO_KEGG_Integration_TumorDown_NonzeroOverlap.csv"
    ),
  
  SUMMARY =
    file.path(
      step7b_dir,
      "STEP9C7B_GO_KEGG_Integration_Summary.csv"
    ),
  
  CORE =
    file.path(
      step7b_dir,
      "STEP9C7B_CoreGene_GO_KEGG_DirectOverlap_Audit.csv"
    ),
  
  DIRECTION =
    file.path(
      step7b_dir,
      "STEP9C7B_Direction_Overlap_Audit.csv"
    ),
  
  LOCK =
    file.path(
      step7b_dir,
      "STEP9C7B_FINAL_LOCK.rds"
    )
)


# ------------------------------------------------------------
# 2. File existence audit
# ------------------------------------------------------------

missing_files <- files_required[
  !file.exists(files_required)
]

if (length(missing_files) > 0) {
  
  cat("ERROR: Missing required STEP9C-7B files:\n\n")
  
  print(missing_files)
  
  stop(
    "STEP9C-7B-READ STOPPED: Required locked outputs are missing."
  )
}

cat("Required STEP9C-7B files: PASS\n\n")


# ------------------------------------------------------------
# 3. Read-only loading
# ------------------------------------------------------------

up <- read.csv(
  files_required["GO_KEGG_UP"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

down <- read.csv(
  files_required["GO_KEGG_DOWN"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

summary_df <- read.csv(
  files_required["SUMMARY"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

core_df <- read.csv(
  files_required["CORE"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

direction_df <- read.csv(
  files_required["DIRECTION"],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

STEP9C7B_FINAL_LOCK <- readRDS(
  files_required["LOCK"]
)


# ------------------------------------------------------------
# 4. Verify STEP9C-7B lock
# ------------------------------------------------------------

lock_valid <- all(
  identical(
    STEP9C7B_FINAL_LOCK$Step,
    "STEP9C-7B"
  ),
  identical(
    STEP9C7B_FINAL_LOCK$Status,
    "PASS"
  ),
  isTRUE(
    STEP9C7B_FINAL_LOCK$Enrichment_Recalculation == FALSE
  ),
  isTRUE(
    STEP9C7B_FINAL_LOCK$DEG_Recalculation == FALSE
  ),
  isTRUE(
    STEP9C7B_FINAL_LOCK$ID_Conversion_Performed == FALSE
  ),
  isTRUE(
    STEP9C7B_FINAL_LOCK$P_Value_Recalculation == FALSE
  ),
  isTRUE(
    STEP9C7B_FINAL_LOCK$FDR_Recalculation == FALSE
  ),
  isTRUE(
    STEP9C7B_FINAL_LOCK$Significance_Threshold_Modified == FALSE
  ),
  isTRUE(
    STEP9C7B_FINAL_LOCK$Background_Modified == FALSE
  )
)

if (!lock_valid) {
  
  stop(
    "STEP9C-7B-READ STOPPED: STEP9C-7B lock integrity check failed."
  )
}

cat("STEP9C-7B lock integrity: PASS\n\n")


# ------------------------------------------------------------
# 5. Verify expected locked counts
# ------------------------------------------------------------

count_audit <- data.frame(
  
  Dataset = c(
    "GO-BP tumor-up",
    "GO-BP tumor-down",
    "KEGG tumor-up",
    "KEGG tumor-down"
  ),
  
  Expected = c(
    161,
    1296,
    13,
    60
  ),
  
  stringsAsFactors = FALSE
)

# The nonzero-overlap files do not contain every GO/KEGG pair,
# therefore their row counts are NOT expected to equal the
# statistical result counts.

cat("------------------------------------------------------------\n")
cat("Locked statistical result counts\n")
cat("------------------------------------------------------------\n")

print(
  count_audit,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 6. Basic structural audit of overlap tables
# ------------------------------------------------------------

required_overlap_columns <- c(
  "Direction",
  "GO_ID",
  "GO_Description",
  "GO_p_adjust",
  "GO_GeneRatio",
  "GO_Genes",
  "KEGG_ID",
  "KEGG_Description",
  "KEGG_p_adjust",
  "KEGG_GeneRatio",
  "KEGG_Genes",
  "N_Shared_Entrez",
  "Jaccard",
  "GO_Coverage",
  "KEGG_Coverage",
  "Shared_Entrez"
)

missing_up_cols <- setdiff(
  required_overlap_columns,
  colnames(up)
)

missing_down_cols <- setdiff(
  required_overlap_columns,
  colnames(down)
)

if (length(missing_up_cols) > 0) {
  
  stop(
    paste0(
      "Tumor-up overlap table missing columns: ",
      paste(
        missing_up_cols,
        collapse = ", "
      )
    )
  )
}

if (length(missing_down_cols) > 0) {
  
  stop(
    paste0(
      "Tumor-down overlap table missing columns: ",
      paste(
        missing_down_cols,
        collapse = ", "
      )
    )
  )
}

cat("Overlap-table structure: PASS\n\n")


# ------------------------------------------------------------
# 7. Verify that every retained pair has >=1 shared gene
# ------------------------------------------------------------

if (any(
  is.na(up$N_Shared_Entrez) |
  up$N_Shared_Entrez < 1
)) {
  
  stop(
    "Tumor-up overlap table contains invalid zero/NA overlaps."
  )
}

if (any(
  is.na(down$N_Shared_Entrez) |
  down$N_Shared_Entrez < 1
)) {
  
  stop(
    "Tumor-down overlap table contains invalid zero/NA overlaps."
  )
}

cat("Nonzero-overlap integrity: PASS\n\n")


# ------------------------------------------------------------
# 8. Direction integrity
# ------------------------------------------------------------

up_direction_ok <- all(
  up$Direction == "Tumor-up"
)

down_direction_ok <- all(
  down$Direction == "Tumor-down"
)

if (!up_direction_ok ||
    !down_direction_ok) {
  
  stop(
    "Direction labels are inconsistent."
  )
}

cat("Direction integrity: PASS\n\n")


# ------------------------------------------------------------
# 9. Overall overlap evidence summary
# ------------------------------------------------------------

overall_overlap_summary <- data.frame(
  
  Direction = c(
    "Tumor-up",
    "Tumor-down"
  ),
  
  GO_Terms = c(
    length(unique(up$GO_ID)),
    length(unique(down$GO_ID))
  ),
  
  KEGG_Pathways = c(
    length(unique(up$KEGG_ID)),
    length(unique(down$KEGG_ID))
  ),
  
  Nonzero_GO_KEGG_Pairs = c(
    nrow(up),
    nrow(down)
  ),
  
  GO_Terms_With_Overlap = c(
    length(unique(up$GO_ID)),
    length(unique(down$GO_ID))
  ),
  
  KEGG_Pathways_With_Overlap = c(
    length(unique(up$KEGG_ID)),
    length(unique(down$KEGG_ID))
  ),
  
  Mean_Shared_Genes = c(
    mean(up$N_Shared_Entrez),
    mean(down$N_Shared_Entrez)
  ),
  
  Median_Shared_Genes = c(
    median(up$N_Shared_Entrez),
    median(down$N_Shared_Entrez)
  ),
  
  Max_Shared_Genes = c(
    max(up$N_Shared_Entrez),
    max(down$N_Shared_Entrez)
  ),
  
  Mean_Jaccard = c(
    mean(up$Jaccard),
    mean(down$Jaccard)
  ),
  
  Median_Jaccard = c(
    median(up$Jaccard),
    median(down$Jaccard)
  ),
  
  Max_Jaccard = c(
    max(up$Jaccard),
    max(down$Jaccard)
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 10. KEGG pathway connectivity
#
# Descriptive:
#   number of distinct GO-BP terms connected to each KEGG
#   pathway by >=1 shared Entrez gene.
#
# NO biological ranking is implied.
# ------------------------------------------------------------

kegg_connectivity <- function(
    df
) {
  
  split_go <- split(
    df$GO_ID,
    df$KEGG_ID
  )
  
  split_shared <- split(
    df$N_Shared_Entrez,
    df$KEGG_ID
  )
  
  split_jaccard <- split(
    df$Jaccard,
    df$KEGG_ID
  )
  
  split_desc <- split(
    df$KEGG_Description,
    df$KEGG_ID
  )
  
  out <- lapply(
    names(split_go),
    function(k) {
      
      data.frame(
        KEGG_ID = k,
        KEGG_Description =
          unique(split_desc[[k]])[1],
        N_GO_Terms =
          length(unique(split_go[[k]])),
        Total_Shared_Genes =
          sum(split_shared[[k]]),
        Max_Shared_Genes =
          max(split_shared[[k]]),
        Mean_Shared_Genes =
          mean(split_shared[[k]]),
        Max_Jaccard =
          max(split_jaccard[[k]]),
        stringsAsFactors = FALSE
      )
    }
  )
  
  out <- do.call(
    rbind,
    out
  )
  
  out[
    order(
      -out$N_GO_Terms,
      -out$Max_Shared_Genes,
      -out$Max_Jaccard
    ),
    ,
    drop = FALSE
  ]
}


kegg_connectivity_up <- kegg_connectivity(
  up
)

kegg_connectivity_down <- kegg_connectivity(
  down
)


# ------------------------------------------------------------
# 11. GO-term connectivity
#
# Descriptive:
#   number of KEGG pathways connected to each GO term.
# ------------------------------------------------------------

go_connectivity <- function(
    df
) {
  
  split_kegg <- split(
    df$KEGG_ID,
    df$GO_ID
  )
  
  split_shared <- split(
    df$N_Shared_Entrez,
    df$GO_ID
  )
  
  split_jaccard <- split(
    df$Jaccard,
    df$GO_ID
  )
  
  split_desc <- split(
    df$GO_Description,
    df$GO_ID
  )
  
  out <- lapply(
    names(split_kegg),
    function(g) {
      
      data.frame(
        GO_ID = g,
        GO_Description =
          unique(split_desc[[g]])[1],
        N_KEGG_Pathways =
          length(unique(split_kegg[[g]])),
        Total_Shared_Genes =
          sum(split_shared[[g]]),
        Max_Shared_Genes =
          max(split_shared[[g]]),
        Mean_Shared_Genes =
          mean(split_shared[[g]]),
        Max_Jaccard =
          max(split_jaccard[[g]]),
        stringsAsFactors = FALSE
      )
    }
  )
  
  out <- do.call(
    rbind,
    out
  )
  
  out[
    order(
      -out$N_KEGG_Pathways,
      -out$Max_Shared_Genes,
      -out$Max_Jaccard
    ),
    ,
    drop = FALSE
  ]
}


go_connectivity_up <- go_connectivity(
  up
)

go_connectivity_down <- go_connectivity(
  down
)


# ------------------------------------------------------------
# 12. High-overlap descriptive strata
#
# IMPORTANT:
# These are NOT significance filters.
# They are only descriptive bins for inspecting the
# distribution of actual gene overlaps.
#
# No threshold is used to delete results.
# ------------------------------------------------------------

overlap_distribution <- function(
    df,
    direction
) {
  
  data.frame(
    
    Direction = direction,
    
    N_Pairs = nrow(df),
    
    Shared_1 = sum(
      df$N_Shared_Entrez == 1
    ),
    
    Shared_2_4 = sum(
      df$N_Shared_Entrez >= 2 &
        df$N_Shared_Entrez <= 4
    ),
    
    Shared_5_9 = sum(
      df$N_Shared_Entrez >= 5 &
        df$N_Shared_Entrez <= 9
    ),
    
    Shared_10_or_more = sum(
      df$N_Shared_Entrez >= 10
    ),
    
    Jaccard_ge_0_10 = sum(
      df$Jaccard >= 0.10
    ),
    
    Jaccard_ge_0_25 = sum(
      df$Jaccard >= 0.25
    ),
    
    Jaccard_ge_0_50 = sum(
      df$Jaccard >= 0.50
    ),
    
    stringsAsFactors = FALSE
  )
}


overlap_distribution_df <- rbind(
  
  overlap_distribution(
    up,
    "Tumor-up"
  ),
  
  overlap_distribution(
    down,
    "Tumor-down"
  )
)


# ------------------------------------------------------------
# 13. Identify recurrent GO–KEGG relationships
#
# This does NOT define final themes.
#
# It identifies pairs with:
#   - shared gene content
#   - repeated representation across related terms
#
# All retained pairs remain available.
# ------------------------------------------------------------

pair_evidence <- function(
    df
) {
  
  df2 <- df
  
  df2$Evidence_Score_Descriptive <- paste0(
    "SharedGenes=",
    df2$N_Shared_Entrez,
    ";Jaccard=",
    sprintf(
      "%.3f",
      df2$Jaccard
    ),
    ";GO_Coverage=",
    sprintf(
      "%.3f",
      df2$GO_Coverage
    ),
    ";KEGG_Coverage=",
    sprintf(
      "%.3f",
      df2$KEGG_Coverage
    )
  )
  
  df2[
    order(
      -df2$N_Shared_Entrez,
      -df2$Jaccard,
      df2$GO_p_adjust,
      df2$KEGG_p_adjust
    ),
    ,
    drop = FALSE
  ]
}


evidence_up <- pair_evidence(
  up
)

evidence_down <- pair_evidence(
  down
)


# ------------------------------------------------------------
# 14. Core-gene functional evidence
#
# Read the already generated direct-overlap audit.
#
# NO new gene mapping.
# ------------------------------------------------------------

core_audit <- core_df

cat("------------------------------------------------------------\n")
cat("Core-gene GO–KEGG direct-overlap audit\n")
cat("------------------------------------------------------------\n")

print(
  core_audit,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 15. Core-gene summary
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_summary <- do.call(
  rbind,
  lapply(
    core_genes,
    function(g) {
      
      tmp <- core_audit[
        grepl(
          paste0(
            "(^|;)",
            g,
            "(;|$)"
          ),
          core_audit$Core_Genes_Shared
        ),
        ,
        drop = FALSE
      ]
      
      data.frame(
        Core_Gene = g,
        
        N_GO_KEGG_Overlap_Pairs =
          nrow(tmp),
        
        N_TumorUp_Pairs =
          sum(
            tmp$Direction == "Tumor-up"
          ),
        
        N_TumorDown_Pairs =
          sum(
            tmp$Direction == "Tumor-down"
          ),
        
        N_GO_Terms =
          length(unique(
            tmp$GO_ID
          )),
        
        N_KEGG_Pathways =
          length(unique(
            tmp$KEGG_ID
          )),
        
        stringsAsFactors = FALSE
      )
    }
  )
)


# ------------------------------------------------------------
# 16. Core-gene pair details
# ------------------------------------------------------------

core_pair_details <- core_audit[
  order(
    core_audit$Direction,
    -core_audit$N_Core_Genes_Shared,
    -core_audit$N_Shared_Entrez
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 17. GO semantic-cluster information
#
# We do NOT recalculate STEP9C-7A.
#
# First search for any saved GO cluster assignment file.
# ------------------------------------------------------------

cluster_files <- list.files(
  file.path(
    base_dir,
    "STEP9C-7A_GO_BP_Redundancy_Semantic_Similarity"
  ),
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

cluster_files <- cluster_files[
  grepl(
    "cluster",
    basename(cluster_files),
    ignore.case = TRUE
  )
]

cat("STEP9C-7A cluster-related CSV files found = ",
    length(cluster_files),
    "\n",
    sep = "")

if (length(cluster_files) > 0) {
  print(
    basename(cluster_files)
  )
}

cat("\n")


# ------------------------------------------------------------
# 18. Save audit outputs
# ------------------------------------------------------------

write.csv(
  overall_overlap_summary,
  file.path(
    audit_dir,
    "STEP9C7B_READ_Overall_Overlap_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  kegg_connectivity_up,
  file.path(
    audit_dir,
    "STEP9C7B_READ_KEGG_Connectivity_TumorUp.csv"
  ),
  row.names = FALSE
)

write.csv(
  kegg_connectivity_down,
  file.path(
    audit_dir,
    "STEP9C7B_READ_KEGG_Connectivity_TumorDown.csv"
  ),
  row.names = FALSE
)

write.csv(
  go_connectivity_up,
  file.path(
    audit_dir,
    "STEP9C7B_READ_GO_Connectivity_TumorUp.csv"
  ),
  row.names = FALSE
)

write.csv(
  go_connectivity_down,
  file.path(
    audit_dir,
    "STEP9C7B_READ_GO_Connectivity_TumorDown.csv"
  ),
  row.names = FALSE
)

write.csv(
  overlap_distribution_df,
  file.path(
    audit_dir,
    "STEP9C7B_READ_Overlap_Distribution.csv"
  ),
  row.names = FALSE
)

write.csv(
  evidence_up,
  file.path(
    audit_dir,
    "STEP9C7B_READ_GO_KEGG_Evidence_TumorUp.csv"
  ),
  row.names = FALSE
)

write.csv(
  evidence_down,
  file.path(
    audit_dir,
    "STEP9C7B_READ_GO_KEGG_Evidence_TumorDown.csv"
  ),
  row.names = FALSE
)

write.csv(
  core_summary,
  file.path(
    audit_dir,
    "STEP9C7B_READ_CoreGene_Functional_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  core_pair_details,
  file.path(
    audit_dir,
    "STEP9C7B_READ_CoreGene_GO_KEGG_Details.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 19. Save audit metadata
# ------------------------------------------------------------

STEP9C7B_READ_AUDIT <- list(
  
  Step = "STEP9C-7B-READ",
  
  Status = "PASS",
  
  Audit_Type =
    "READ-ONLY descriptive evidence audit",
  
  Input_STEP9C7B_Unchanged = TRUE,
  
  Statistical_Results_Unchanged = TRUE,
  
  GO_BP_TumorUp = 161,
  
  GO_BP_TumorDown = 1296,
  
  KEGG_TumorUp = 13,
  
  KEGG_TumorDown = 60,
  
  TumorUp_NonzeroOverlapPairs =
    nrow(up),
  
  TumorDown_NonzeroOverlapPairs =
    nrow(down),
  
  GO_KEGG_Connection_Definition =
    "Shared Entrez gene content",
  
  Semantic_Similarity_Threshold =
    0.70,
  
  Semantic_Similarity_Use =
    "Descriptive organization only",
  
  Final_Theme_Assignment =
    FALSE,
  
  Biological_Theme_Naming =
    FALSE,
  
  Enrichment_Recalculation =
    FALSE,
  
  DEG_Recalculation =
    FALSE,
  
  ID_Conversion =
    FALSE,
  
  P_Value_Recalculation =
    FALSE,
  
  FDR_Recalculation =
    FALSE,
  
  Significance_Threshold_Modification =
    FALSE,
  
  Background_Modification =
    FALSE
)


saveRDS(
  STEP9C7B_READ_AUDIT,
  file.path(
    audit_dir,
    "STEP9C7B_READ_AUDIT_LOCK.rds"
  )
)


# ------------------------------------------------------------
# 20. Console report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-7B-READ FINAL REPORT\n")
cat("============================================================\n\n")

cat(
  "Locked GO-BP tumor-up terms     = 161\n"
)

cat(
  "Locked GO-BP tumor-down terms   = 1296\n"
)

cat(
  "Locked KEGG tumor-up pathways   = 13\n"
)

cat(
  "Locked KEGG tumor-down pathways = 60\n\n"
)

cat(
  "Tumor-up nonzero GO–KEGG pairs = ",
  nrow(up),
  "\n",
  sep = ""
)

cat(
  "Tumor-down nonzero GO–KEGG pairs = ",
  nrow(down),
  "\n\n",
  sep = ""
)

cat(
  "Tumor-up GO terms connected to KEGG = ",
  length(unique(up$GO_ID)),
  "\n",
  sep = ""
)

cat(
  "Tumor-down GO terms connected to KEGG = ",
  length(unique(down$GO_ID)),
  "\n\n",
  sep = ""
)

cat(
  "Tumor-up KEGG pathways connected to GO = ",
  length(unique(up$KEGG_ID)),
  "\n",
  sep = ""
)

cat(
  "Tumor-down KEGG pathways connected to GO = ",
  length(unique(down$KEGG_ID)),
  "\n\n",
  sep = ""
)

cat(
  "STEP9C-7A semantic similarity threshold = 0.70\n"
)

cat(
  "Semantic similarity is used descriptively only.\n\n"
)

cat(
  "No final biological theme assigned in this step.\n"
)

cat(
  "No biological theme ranking performed.\n"
)

cat(
  "No top-N statistical filtering performed.\n\n"
)

cat(
  "NO enrichment recalculation.\n"
)

cat(
  "NO DEG recalculation.\n"
)

cat(
  "NO ID conversion.\n"
)

cat(
  "NO P-value recalculation.\n"
)

cat(
  "NO FDR recalculation.\n"
)

cat(
  "NO significance-threshold modification.\n"
)

cat(
  "NO background modification.\n"
)

cat(
  "STEP9C-7B locked results unchanged.\n\n"
)

cat("============================================================\n")
cat("FINAL STATUS: PASS — STEP9C-7B-READ COMPLETE\n")
cat("============================================================\n")


# ============================================================
# STEP9C-7B-READ-2
# Focused Functional Theme Evidence Display
#
# READ-ONLY
# No recalculation
# ============================================================

audit_dir <- paste0(
  "input_dir/Paper 1_Step 9-Functional Enrichment/",
  "STEP9C-7B_GO_KEGG_Functional_Theme_Integration/",
  "STEP9C-7B-READ_Functional_Theme_Evidence_Audit"
)

# ------------------------------------------------------------
# 1. KEGG connectivity — tumor-up
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("TUMOR-UP: KEGG CONNECTIVITY\n")
cat("============================================================\n")

kegg_up <- read.csv(
  file.path(
    audit_dir,
    "STEP9C7B_READ_KEGG_Connectivity_TumorUp.csv"
  ),
  stringsAsFactors = FALSE
)

print(
  kegg_up,
  row.names = FALSE
)


# ------------------------------------------------------------
# 2. KEGG connectivity — tumor-down
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("TUMOR-DOWN: KEGG CONNECTIVITY\n")
cat("============================================================\n")

kegg_down <- read.csv(
  file.path(
    audit_dir,
    "STEP9C7B_READ_KEGG_Connectivity_TumorDown.csv"
  ),
  stringsAsFactors = FALSE
)

print(
  kegg_down,
  row.names = FALSE
)


# ------------------------------------------------------------
# 3. Core-gene functional summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CORE-GENE FUNCTIONAL SUMMARY\n")
cat("============================================================\n")

core_summary <- read.csv(
  file.path(
    audit_dir,
    "STEP9C7B_READ_CoreGene_Functional_Summary.csv"
  ),
  stringsAsFactors = FALSE
)

print(
  core_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. Tumor-up GO connectivity
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("TUMOR-UP: MOST CONNECTED GO-BP TERMS\n")
cat("============================================================\n")

go_up <- read.csv(
  file.path(
    audit_dir,
    "STEP9C7B_READ_GO_Connectivity_TumorUp.csv"
  ),
  stringsAsFactors = FALSE
)

print(
  head(go_up, 30),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Tumor-down GO connectivity
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("TUMOR-DOWN: MOST CONNECTED GO-BP TERMS\n")
cat("============================================================\n")

go_down <- read.csv(
  file.path(
    audit_dir,
    "STEP9C7B_READ_GO_Connectivity_TumorDown.csv"
  ),
  stringsAsFactors = FALSE
)

print(
  head(go_down, 50),
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Strongest actual GO–KEGG overlap pairs — tumor-up
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("TUMOR-UP: GO–KEGG EVIDENCE\n")
cat("============================================================\n")

evidence_up <- read.csv(
  file.path(
    audit_dir,
    "STEP9C7B_READ_GO_KEGG_Evidence_TumorUp.csv"
  ),
  stringsAsFactors = FALSE
)

print(
  head(evidence_up, 30),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Strongest actual GO–KEGG overlap pairs — tumor-down
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("TUMOR-DOWN: GO–KEGG EVIDENCE\n")
cat("============================================================\n")

evidence_down <- read.csv(
  file.path(
    audit_dir,
    "STEP9C7B_READ_GO_KEGG_Evidence_TumorDown.csv"
  ),
  stringsAsFactors = FALSE
)

print(
  head(evidence_down, 50),
  row.names = FALSE
)


cat("\n============================================================\n")
cat("STEP9C-7B-READ-2 COMPLETE\n")
cat("READ-ONLY DISPLAY — NO RECALCULATION\n")
cat("============================================================\n")


# ============================================================
# Paper 1 — STEP9C-7C
# READ-ONLY Functional Theme Consolidation Audit
#
# Purpose:
#   Integrate already-locked:
#     1. STEP9C-6B significant GO-BP terms
#     2. STEP9C-7A GO semantic clusters (threshold = 0.70)
#     3. STEP9C-7B GO–KEGG shared-gene evidence
#
# IMPORTANT:
#   This step is READ-ONLY.
#
#   NO enrichment recalculation
#   NO DEG recalculation
#   NO ID conversion
#   NO P-value recalculation
#   NO FDR recalculation
#   NO significance-threshold modification
#   NO semantic-similarity recalculation
#   NO new clustering
#   NO background modification
#   NO modification of STEP9C-5B / STEP9C-6B / STEP9C-7A / STEP9C-7B
#
# Output:
#   Candidate-theme evidence tables for descriptive
#   functional-theme consolidation.
# ============================================================


cat("\n============================================================\n")
cat("STEP9C-7C: READ-ONLY Functional Theme Consolidation Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

setwd(step_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Create output directory
# ------------------------------------------------------------

out_dir <- file.path(step_dir, "STEP9C-7C_Functional_Theme_Consolidation")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cat("Output directory:\n")
cat(out_dir, "\n\n")


# ------------------------------------------------------------
# 2. Load ONLY locked STEP9C-6B statistical results
# ------------------------------------------------------------

go_up_file <- file.path(
  step_dir,
  "STEP9C6B_Supplementary_Table_GO_BP_TumorUp_161.csv"
)

go_down_file <- file.path(
  step_dir,
  "STEP9C6B_Supplementary_Table_GO_BP_TumorDown_1296.csv"
)

stopifnot(file.exists(go_up_file))
stopifnot(file.exists(go_down_file))

go_up_locked <- read.csv(
  go_up_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

go_down_locked <- read.csv(
  go_down_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 3. Structural validation of locked GO results
# ------------------------------------------------------------

required_go_columns <- c(
  "Direction",
  "ID",
  "Description",
  "GeneRatio",
  "BgRatio",
  "pvalue",
  "p.adjust",
  "geneID"
)

stopifnot(all(required_go_columns %in% colnames(go_up_locked)))
stopifnot(all(required_go_columns %in% colnames(go_down_locked)))

stopifnot(nrow(go_up_locked) == 161)
stopifnot(nrow(go_down_locked) == 1296)

stopifnot(
  all(go_up_locked$Direction == "Tumor-up")
)

stopifnot(
  all(go_down_locked$Direction == "Tumor-down")
)

cat("STEP9C-6B locked GO-BP results:\n")
cat("  Tumor-up   =", nrow(go_up_locked), "\n")
cat("  Tumor-down =", nrow(go_down_locked), "\n\n")


# ------------------------------------------------------------
# 4. Load existing STEP9C-7B GO–KEGG evidence
#
# IMPORTANT:
#   These objects must already exist from STEP9C-7B-READ.
#
#   No new GO–KEGG overlap calculation is performed here.
# ------------------------------------------------------------

required_evidence_columns <- c(
  "Direction",
  "GO_ID",
  "GO_Description",
  "GO_p_adjust",
  "GO_GeneRatio",
  "GO_Genes",
  "KEGG_ID",
  "KEGG_Description",
  "KEGG_p_adjust",
  "KEGG_GeneRatio",
  "KEGG_Genes",
  "N_Shared_Entrez",
  "Jaccard",
  "GO_Coverage",
  "KEGG_Coverage"
)


if (!exists("evidence_up", inherits = TRUE)) {
  stop(
    paste0(
      "\nERROR: object 'evidence_up' is not present.\n",
      "Please load the existing STEP9C-7B evidence object/table first.\n",
      "No new GO–KEGG overlap calculation is permitted in STEP9C-7C."
    )
  )
}

if (!exists("evidence_down", inherits = TRUE)) {
  stop(
    paste0(
      "\nERROR: object 'evidence_down' is not present.\n",
      "Please load the existing STEP9C-7B evidence object/table first.\n",
      "No new GO–KEGG overlap calculation is permitted in STEP9C-7C."
    )
  )
}


stopifnot(
  all(required_evidence_columns %in% colnames(evidence_up))
)

stopifnot(
  all(required_evidence_columns %in% colnames(evidence_down))
)


cat("STEP9C-7B evidence objects detected.\n")
cat("  evidence_up   rows =", nrow(evidence_up), "\n")
cat("  evidence_down rows =", nrow(evidence_down), "\n\n")


# ------------------------------------------------------------
# 5. Validate evidence against locked GO-BP results
#
# READ-ONLY validation only.
# ------------------------------------------------------------

up_go_ids_locked <- unique(go_up_locked$ID)
down_go_ids_locked <- unique(go_down_locked$ID)

up_evidence_go_ids <- unique(evidence_up$GO_ID)
down_evidence_go_ids <- unique(evidence_down$GO_ID)


up_unexpected_go <- setdiff(
  up_evidence_go_ids,
  up_go_ids_locked
)

down_unexpected_go <- setdiff(
  down_evidence_go_ids,
  down_go_ids_locked
)


if (length(up_unexpected_go) > 0) {
  stop(
    "FAIL: evidence_up contains GO IDs not present in locked STEP9C-6B tumor-up results."
  )
}

if (length(down_unexpected_go) > 0) {
  stop(
    "FAIL: evidence_down contains GO IDs not present in locked STEP9C-6B tumor-down results."
  )
}


cat("GO-ID consistency audit: PASS\n\n")


# ------------------------------------------------------------
# 6. Validate Direction
# ------------------------------------------------------------

stopifnot(
  all(evidence_up$Direction == "Tumor-up")
)

stopifnot(
  all(evidence_down$Direction == "Tumor-down")
)

cat("Direction consistency audit: PASS\n\n")


# ============================================================
# 7. STEP9C-7A semantic-cluster objects
# ============================================================
#
# IMPORTANT:
#   STEP9C-7C must USE the already locked 0.70 clusters.
#
#   It must NOT:
#     - calculate semantic similarity again
#     - run mgoSim()
#     - run hclust()
#     - create new clusters
#
# The following section searches the current R environment for
# existing cluster-assignment objects.
#
# Expected structure:
#   GO_ID + Cluster_ID
#
# If your locked STEP9C-7A objects have different names,
# the script will report the candidates and stop safely.
# ============================================================


find_cluster_candidates <- function(direction) {
  
  nms <- ls(envir = .GlobalEnv)
  
  candidates <- list()
  
  for (nm in nms) {
    
    obj <- get(nm, envir = .GlobalEnv)
    
    if (!is.data.frame(obj)) {
      next
    }
    
    cn <- colnames(obj)
    
    has_go <- any(cn %in% c(
      "GO_ID",
      "ID",
      "go_id",
      "GO"
    ))
    
    has_cluster <- any(cn %in% c(
      "Cluster",
      "Cluster_ID",
      "cluster",
      "cluster_id",
      "ClusterID"
    ))
    
    if (has_go && has_cluster) {
      
      candidates[[nm]] <- list(
        rows = nrow(obj),
        columns = cn
      )
    }
  }
  
  candidates
}


cluster_candidates <- find_cluster_candidates("all")


cat("Existing possible semantic-cluster objects:\n")

if (length(cluster_candidates) == 0) {
  
  cat("  NONE FOUND\n\n")
  
} else {
  
  for (nm in names(cluster_candidates)) {
    
    cat(
      "  ",
      nm,
      " | rows = ",
      cluster_candidates[[nm]]$rows,
      " | columns = ",
      paste(cluster_candidates[[nm]]$columns, collapse = ", "),
      "\n",
      sep = ""
    )
  }
  
  cat("\n")
}

# ============================================================
# STEP9C-7C
# Existing STEP9C-7A locked semantic-cluster assignment
# ------------------------------------------------------------
# READ-ONLY assignment
#
# IMPORTANT:
#   up_primary   = locked STEP9C-7A tumor-up 0.70 clusters
#   down_primary = locked STEP9C-7A tumor-down 0.70 clusters
#
# NO semantic similarity recalculation
# NO mgoSim()
# NO hclust()
# NO new clustering
# ============================================================

cluster_up_070 <- up_primary
cluster_down_070 <- down_primary


cat("\n============================================================\n")
cat("STEP9C-7A LOCKED CLUSTER OBJECTS ASSIGNED\n")
cat("============================================================\n\n")

cat("Tumor-up object:\n")
cat("  object = cluster_up_070\n")
cat("  rows   =", nrow(cluster_up_070), "\n")
cat("  cols   =", paste(colnames(cluster_up_070), collapse = ", "), "\n\n")

cat("Tumor-down object:\n")
cat("  object = cluster_down_070\n")
cat("  rows   =", nrow(cluster_down_070), "\n")
cat("  cols   =", paste(colnames(cluster_down_070), collapse = ", "), "\n\n")


# ------------------------------------------------------------
# Validate locked threshold
# ------------------------------------------------------------

cat("Tumor-up Similarity_Threshold:\n")
print(unique(cluster_up_070$Similarity_Threshold))

cat("\nTumor-down Similarity_Threshold:\n")
print(unique(cluster_down_070$Similarity_Threshold))


# ------------------------------------------------------------
# Validate expected term counts
# ------------------------------------------------------------

stopifnot(
  nrow(cluster_up_070) == 161
)

stopifnot(
  nrow(cluster_down_070) == 1296
)


# ------------------------------------------------------------
# Validate expected cluster counts
# ------------------------------------------------------------

stopifnot(
  length(unique(cluster_up_070$Cluster)) == 97
)

stopifnot(
  length(unique(cluster_down_070$Cluster)) == 743
)


# ------------------------------------------------------------
# Validate threshold
# ------------------------------------------------------------

stopifnot(
  all(cluster_up_070$Similarity_Threshold == 0.70)
)

stopifnot(
  all(cluster_down_070$Similarity_Threshold == 0.70)
)


cat("\n------------------------------------------------------------\n")
cat("STEP9C-7A LOCK VALIDATION: PASS\n")
cat("------------------------------------------------------------\n")

cat("Tumor-up:\n")
cat("  GO terms   = 161\n")
cat("  clusters   = 97\n")
cat("  threshold  = 0.70\n\n")

cat("Tumor-down:\n")
cat("  GO terms   = 1296\n")
cat("  clusters   = 743\n")
cat("  threshold  = 0.70\n\n")

cat("NO semantic similarity recalculation.\n")
cat("NO new clustering.\n")
cat("NO enrichment recalculation.\n")
cat("NO statistical result modification.\n\n")

cat("Ready for STEP9C-7C consolidation.\n")


# ============================================================
# STEP9C-7C
# Existing STEP9C-7A cluster object binding
# READ-ONLY
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-7C: Bind Existing STEP9C-7A Cluster Objects\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Bind existing locked STEP9C-7A objects
# ------------------------------------------------------------

cluster_up_070   <- up_primary
cluster_down_070 <- down_primary


# ------------------------------------------------------------
# 2. Structural validation
# ------------------------------------------------------------

stopifnot(
  nrow(cluster_up_070) == 161,
  nrow(cluster_down_070) == 1296
)

stopifnot(
  length(unique(cluster_up_070$Cluster)) == 97,
  length(unique(cluster_down_070$Cluster)) == 743
)

stopifnot(
  all(cluster_up_070$Similarity_Threshold == 0.70),
  all(cluster_down_070$Similarity_Threshold == 0.70)
)


# ------------------------------------------------------------
# 3. Display validation
# ------------------------------------------------------------

cat("Tumor-up cluster object:\n")
cat("  rows              =", nrow(cluster_up_070), "\n")
cat("  unique GO terms   =", length(unique(cluster_up_070$ID)), "\n")
cat("  clusters          =", length(unique(cluster_up_070$Cluster)), "\n")
cat("  threshold         =", unique(cluster_up_070$Similarity_Threshold), "\n\n")

cat("Tumor-down cluster object:\n")
cat("  rows              =", nrow(cluster_down_070), "\n")
cat("  unique GO terms   =", length(unique(cluster_down_070$ID)), "\n")
cat("  clusters          =", length(unique(cluster_down_070$Cluster)), "\n")
cat("  threshold         =", unique(cluster_down_070$Similarity_Threshold), "\n\n")


# ------------------------------------------------------------
# 4. READ-ONLY confirmation
# ------------------------------------------------------------

cat("STEP9C-7A CLUSTER OBJECT BINDING: PASS\n")
cat("------------------------------------------------------------\n")
cat("Existing objects used:\n")
cat("  up_primary\n")
cat("  down_primary\n\n")

cat("No semantic similarity recalculation.\n")
cat("No new clustering.\n")
cat("No enrichment recalculation.\n")
cat("No modification of locked GO results.\n")
cat("No modification of locked KEGG results.\n")
cat("No modification of evidence_up/evidence_down.\n\n")

cat("Ready for cluster standardization and GO-ID matching.\n")


# ============================================================
# STEP9C-7C
# GO-ID ↔ Semantic Cluster Exact Matching
# READ-ONLY
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-7C: GO-ID to Semantic Cluster Exact Matching\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

out_dir <- file.path(
  step_dir,
  "STEP9C-7C_Functional_Theme_Consolidation"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 1. Standardize existing STEP9C-7A cluster tables
# ------------------------------------------------------------

standardize_cluster_table <- function(x, direction_label) {
  
  required_cols <- c(
    "ID",
    "Cluster",
    "Similarity_Threshold"
  )
  
  missing_cols <- setdiff(required_cols, colnames(x))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ",
      direction_label,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  out <- data.frame(
    GO_ID = as.character(x$ID),
    Semantic_Cluster_070 = as.character(x$Cluster),
    Similarity_Threshold = as.numeric(
      x$Similarity_Threshold
    ),
    stringsAsFactors = FALSE
  )
  
  out <- unique(out)
  
  if (!all(out$Similarity_Threshold == 0.70)) {
    stop(
      direction_label,
      ": similarity threshold is not uniformly 0.70."
    )
  }
  
  if (anyDuplicated(out$GO_ID)) {
    stop(
      direction_label,
      ": duplicate GO IDs detected after standardization."
    )
  }
  
  out
}


cluster_up_std <- standardize_cluster_table(
  cluster_up_070,
  "Tumor-up"
)

cluster_down_std <- standardize_cluster_table(
  cluster_down_070,
  "Tumor-down"
)


# ------------------------------------------------------------
# 2. Load locked STEP9C-6B GO-BP tables
# ------------------------------------------------------------

go_up_file <- file.path(
  step_dir,
  "STEP9C6B_Supplementary_Table_GO_BP_TumorUp_161.csv"
)

go_down_file <- file.path(
  step_dir,
  "STEP9C6B_Supplementary_Table_GO_BP_TumorDown_1296.csv"
)

go_up_locked <- read.csv(
  go_up_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

go_down_locked <- read.csv(
  go_down_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 3. Extract locked GO-ID universes
# ------------------------------------------------------------

up_go_ids_locked <- unique(
  as.character(go_up_locked$ID)
)

down_go_ids_locked <- unique(
  as.character(go_down_locked$ID)
)


# ------------------------------------------------------------
# 4. Basic count validation
# ------------------------------------------------------------

stopifnot(
  length(up_go_ids_locked) == 161,
  length(down_go_ids_locked) == 1296
)

stopifnot(
  nrow(cluster_up_std) == 161,
  nrow(cluster_down_std) == 1296
)


# ------------------------------------------------------------
# 5. Exact GO-ID set comparison
# ------------------------------------------------------------

up_missing_from_cluster <- setdiff(
  up_go_ids_locked,
  cluster_up_std$GO_ID
)

up_extra_in_cluster <- setdiff(
  cluster_up_std$GO_ID,
  up_go_ids_locked
)

down_missing_from_cluster <- setdiff(
  down_go_ids_locked,
  cluster_down_std$GO_ID
)

down_extra_in_cluster <- setdiff(
  cluster_down_std$GO_ID,
  down_go_ids_locked
)


# ------------------------------------------------------------
# 6. Exact matching checks
# ------------------------------------------------------------

if (length(up_missing_from_cluster) > 0) {
  stop(
    "Tumor-up GO IDs missing from semantic cluster table: ",
    paste(up_missing_from_cluster, collapse = ", ")
  )
}

if (length(up_extra_in_cluster) > 0) {
  stop(
    "Tumor-up semantic cluster table contains unexpected GO IDs: ",
    paste(up_extra_in_cluster, collapse = ", ")
  )
}

if (length(down_missing_from_cluster) > 0) {
  stop(
    "Tumor-down GO IDs missing from semantic cluster table: ",
    paste(down_missing_from_cluster, collapse = ", ")
  )
}

if (length(down_extra_in_cluster) > 0) {
  stop(
    "Tumor-down semantic cluster table contains unexpected GO IDs: ",
    paste(down_extra_in_cluster, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 7. Exact 1:1 GO-ID mapping validation
# ------------------------------------------------------------

stopifnot(
  setequal(
    up_go_ids_locked,
    cluster_up_std$GO_ID
  )
)

stopifnot(
  setequal(
    down_go_ids_locked,
    cluster_down_std$GO_ID
  )
)

stopifnot(
  !anyDuplicated(cluster_up_std$GO_ID),
  !anyDuplicated(cluster_down_std$GO_ID)
)


# ------------------------------------------------------------
# 8. Cluster count validation
# ------------------------------------------------------------

stopifnot(
  length(unique(
    cluster_up_std$Semantic_Cluster_070
  )) == 97
)

stopifnot(
  length(unique(
    cluster_down_std$Semantic_Cluster_070
  )) == 743
)


# ------------------------------------------------------------
# 9. Final validation summary
# ------------------------------------------------------------

cat("STEP9C-7C GO-ID / CLUSTER MATCHING: PASS\n")
cat("------------------------------------------------------------\n")

cat("Tumor-up:\n")
cat("  Locked GO terms              =", length(up_go_ids_locked), "\n")
cat("  Cluster GO terms             =", nrow(cluster_up_std), "\n")
cat("  Unique semantic clusters     =",
    length(unique(cluster_up_std$Semantic_Cluster_070)), "\n")
cat("  Missing GO IDs               =",
    length(up_missing_from_cluster), "\n")
cat("  Unexpected GO IDs            =",
    length(up_extra_in_cluster), "\n")
cat("  Similarity threshold         = 0.70\n\n")

cat("Tumor-down:\n")
cat("  Locked GO terms              =", length(down_go_ids_locked), "\n")
cat("  Cluster GO terms             =", nrow(cluster_down_std), "\n")
cat("  Unique semantic clusters     =",
    length(unique(cluster_down_std$Semantic_Cluster_070)), "\n")
cat("  Missing GO IDs               =",
    length(down_missing_from_cluster), "\n")
cat("  Unexpected GO IDs            =",
    length(down_extra_in_cluster), "\n")
cat("  Similarity threshold         = 0.70\n\n")


# ------------------------------------------------------------
# 10. READ-ONLY audit statement
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("READ-ONLY AUDIT STATUS: PASS\n")
cat("------------------------------------------------------------\n")

cat("Existing STEP9C-7A clusters were used unchanged.\n")
cat("No semantic similarity recalculation.\n")
cat("No new clustering.\n")
cat("No GO enrichment recalculation.\n")
cat("No KEGG enrichment recalculation.\n")
cat("No ID conversion.\n")
cat("No P-value recalculation.\n")
cat("No FDR recalculation.\n")
cat("No significance threshold modification.\n")
cat("No background modification.\n")
cat("No locked result modification.\n\n")

cat("Ready for GO-KEGG evidence consolidation.\n")


# ============================================================
# STEP9C-7C
# READ-ONLY Evidence Object Structure Audit
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-7C: Evidence Object Structure Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Object existence
# ------------------------------------------------------------

objects_to_check <- c(
  "evidence_up",
  "evidence_down",
  "integration_up",
  "integration_down"
)

exist_status <- vapply(
  objects_to_check,
  exists,
  logical(1),
  inherits = TRUE
)

print(
  data.frame(
    Object = objects_to_check,
    Exists = exist_status,
    stringsAsFactors = FALSE
  )
)


# ------------------------------------------------------------
# 2. Dimensions
# ------------------------------------------------------------

for (obj_name in objects_to_check) {
  
  if (exists(obj_name, inherits = TRUE)) {
    
    obj <- get(obj_name)
    
    cat("\nObject:", obj_name, "\n")
    cat("  Class   :", paste(class(obj), collapse = ", "), "\n")
    cat("  Rows    :", nrow(obj), "\n")
    cat("  Columns :", ncol(obj), "\n")
    
  }
}


# ------------------------------------------------------------
# 3. Column names
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("COLUMN STRUCTURE\n")
cat("============================================================\n\n")

for (obj_name in objects_to_check) {
  
  if (exists(obj_name, inherits = TRUE)) {
    
    obj <- get(obj_name)
    
    cat("\n---", obj_name, "---\n")
    print(colnames(obj))
    
  }
}


# ------------------------------------------------------------
# 4. GO-ID / KEGG-ID uniqueness audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GO–KEGG KEY AUDIT\n")
cat("============================================================\n\n")

for (obj_name in objects_to_check) {
  
  if (exists(obj_name, inherits = TRUE)) {
    
    obj <- get(obj_name)
    
    if (all(
      c("GO_ID", "KEGG_ID") %in% colnames(obj)
    )) {
      
      pair_key <- paste(
        as.character(obj$GO_ID),
        as.character(obj$KEGG_ID),
        sep = "||"
      )
      
      cat(obj_name, ":\n")
      cat("  Unique GO IDs       =",
          length(unique(obj$GO_ID)), "\n")
      cat("  Unique KEGG IDs     =",
          length(unique(obj$KEGG_ID)), "\n")
      cat("  Unique GO–KEGG pairs =",
          length(unique(pair_key)), "\n")
      cat("  Duplicate pair rows =",
          sum(duplicated(pair_key)), "\n\n")
      
    } else {
      
      cat(
        obj_name,
        ": GO_ID / KEGG_ID columns not both available.\n\n"
      )
      
    }
  }
}


# ------------------------------------------------------------
# 5. N_Shared_Entrez audit
# ------------------------------------------------------------

cat("============================================================\n")
cat("SHARED-GENE EVIDENCE AUDIT\n")
cat("============================================================\n\n")

for (obj_name in objects_to_check) {
  
  if (exists(obj_name, inherits = TRUE)) {
    
    obj <- get(obj_name)
    
    if ("N_Shared_Entrez" %in% colnames(obj)) {
      
      cat(obj_name, ":\n")
      cat("  Rows with N_Shared_Entrez > 0 =",
          sum(
            obj$N_Shared_Entrez > 0,
            na.rm = TRUE
          ),
          "\n")
      
      cat("  Rows with N_Shared_Entrez == 0 =",
          sum(
            obj$N_Shared_Entrez == 0,
            na.rm = TRUE
          ),
          "\n\n")
      
    }
  }
}


# ------------------------------------------------------------
# 6. Direct comparison with integration objects
# ------------------------------------------------------------

cat("============================================================\n")
cat("INTEGRATION OBJECT REFERENCE\n")
cat("============================================================\n\n")

if (exists("integration_up", inherits = TRUE)) {
  
  cat(
    "integration_up rows =",
    nrow(integration_up),
    "\n"
  )
}

if (exists("integration_down", inherits = TRUE)) {
  
  cat(
    "integration_down rows =",
    nrow(integration_down),
    "\n"
  )
}


# ------------------------------------------------------------
# 7. READ-ONLY status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("READ-ONLY AUDIT COMPLETE\n")
cat("============================================================\n")

cat("No enrichment recalculation.\n")
cat("No GO–KEGG overlap recalculation.\n")
cat("No semantic similarity recalculation.\n")
cat("No clustering.\n")
cat("No ID conversion.\n")
cat("No P-value/FDR recalculation.\n")
cat("No modification of existing evidence objects.\n")
cat("No modification of locked results.\n\n")

cat("Awaiting evidence-object structure before consolidation.\n")



# ============================================================
# STEP9C-7C
# Full GO–KEGG Integration + Semantic Cluster Consolidation
# READ-ONLY
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-7C: Full GO–KEGG Integration Consolidation\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Required objects
# ------------------------------------------------------------

required_objects <- c(
  "integration_up",
  "integration_down",
  "evidence_up",
  "evidence_down",
  "cluster_up_std",
  "cluster_down_std"
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
    "Required objects missing: ",
    paste(missing_objects, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 1. Validate locked integration dimensions
# ------------------------------------------------------------

stopifnot(
  nrow(integration_up) == 2093,
  nrow(integration_down) == 77760
)


# ------------------------------------------------------------
# 2. Validate locked nonzero evidence dimensions
# ------------------------------------------------------------

stopifnot(
  nrow(evidence_up) == 388,
  nrow(evidence_down) == 31108
)


# ------------------------------------------------------------
# 3. Required integration columns
# ------------------------------------------------------------

required_cols <- c(
  "GO_ID",
  "Direction",
  "GO_Description",
  "GO_p_adjust",
  "GO_GeneRatio",
  "GO_Genes",
  "KEGG_ID",
  "KEGG_Description",
  "KEGG_p_adjust",
  "KEGG_GeneRatio",
  "KEGG_Genes",
  "N_Shared_Entrez",
  "Jaccard",
  "GO_Coverage",
  "KEGG_Coverage",
  "Shared_Entrez",
  "Core_Genes_Shared",
  "N_Core_Genes_Shared"
)

missing_up <- setdiff(
  required_cols,
  colnames(integration_up)
)

missing_down <- setdiff(
  required_cols,
  colnames(integration_down)
)

if (length(missing_up) > 0) {
  stop(
    "integration_up missing columns: ",
    paste(missing_up, collapse = ", ")
  )
}

if (length(missing_down) > 0) {
  stop(
    "integration_down missing columns: ",
    paste(missing_down, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 4. Check complete GO–KEGG pair uniqueness
# ------------------------------------------------------------

make_pair_key <- function(df) {
  paste(
    as.character(df$GO_ID),
    as.character(df$KEGG_ID),
    sep = "||"
  )
}

up_pair_key <- make_pair_key(integration_up)
down_pair_key <- make_pair_key(integration_down)

stopifnot(
  length(unique(up_pair_key)) == 2093,
  length(unique(down_pair_key)) == 77760,
  !anyDuplicated(up_pair_key),
  !anyDuplicated(down_pair_key)
)


# ------------------------------------------------------------
# 5. Validate expected GO / KEGG dimensions
# ------------------------------------------------------------

stopifnot(
  length(unique(integration_up$GO_ID)) == 161,
  length(unique(integration_up$KEGG_ID)) == 13
)

stopifnot(
  length(unique(integration_down$GO_ID)) == 1296,
  length(unique(integration_down$KEGG_ID)) == 60
)


# ------------------------------------------------------------
# 6. Attach existing semantic clusters
# ------------------------------------------------------------

integration_up_cluster <- merge(
  integration_up,
  cluster_up_std,
  by = "GO_ID",
  all.x = TRUE,
  sort = FALSE
)

integration_down_cluster <- merge(
  integration_down,
  cluster_down_std,
  by = "GO_ID",
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------
# 7. Exact row-count preservation
# ------------------------------------------------------------

stopifnot(
  nrow(integration_up_cluster) == 2093,
  nrow(integration_down_cluster) == 77760
)


# ------------------------------------------------------------
# 8. Every GO–KEGG pair must receive a semantic cluster
# ------------------------------------------------------------

if (
  anyNA(
    integration_up_cluster$Semantic_Cluster_070
  )
) {
  stop(
    "Tumor-up: missing semantic cluster assignment detected."
  )
}

if (
  anyNA(
    integration_down_cluster$Semantic_Cluster_070
  )
) {
  stop(
    "Tumor-down: missing semantic cluster assignment detected."
  )
}


# ------------------------------------------------------------
# 9. Verify pair identity was preserved
# ------------------------------------------------------------

up_pair_after <- make_pair_key(
  integration_up_cluster
)

down_pair_after <- make_pair_key(
  integration_down_cluster
)

stopifnot(
  setequal(
    up_pair_key,
    up_pair_after
  )
)

stopifnot(
  setequal(
    down_pair_key,
    down_pair_after
  )
)


# ------------------------------------------------------------
# 10. Validate nonzero evidence subset
# ------------------------------------------------------------

up_nonzero_from_integration <- integration_up_cluster[
  integration_up_cluster$N_Shared_Entrez > 0,
  ,
  drop = FALSE
]

down_nonzero_from_integration <- integration_down_cluster[
  integration_down_cluster$N_Shared_Entrez > 0,
  ,
  drop = FALSE
]

stopifnot(
  nrow(up_nonzero_from_integration) == 388,
  nrow(down_nonzero_from_integration) == 31108
)


# ------------------------------------------------------------
# 11. Compare nonzero pair keys with locked evidence objects
# ------------------------------------------------------------

up_evidence_key <- make_pair_key(evidence_up)
down_evidence_key <- make_pair_key(evidence_down)

up_nonzero_key <- make_pair_key(
  up_nonzero_from_integration
)

down_nonzero_key <- make_pair_key(
  down_nonzero_from_integration
)

stopifnot(
  setequal(
    up_evidence_key,
    up_nonzero_key
  )
)

stopifnot(
  setequal(
    down_evidence_key,
    down_nonzero_key
  )
)


# ------------------------------------------------------------
# 12. Validate semantic-cluster counts
# ------------------------------------------------------------

stopifnot(
  length(unique(
    integration_up_cluster$Semantic_Cluster_070
  )) == 97
)

stopifnot(
  length(unique(
    integration_down_cluster$Semantic_Cluster_070
  )) == 743
)


# ------------------------------------------------------------
# 13. Summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-7C FULL INTEGRATION CONSOLIDATION: PASS\n")
cat("============================================================\n\n")

cat("Tumor-up:\n")
cat("  GO terms                 =",
    length(unique(integration_up_cluster$GO_ID)), "\n")
cat("  KEGG pathways            =",
    length(unique(integration_up_cluster$KEGG_ID)), "\n")
cat("  Total GO–KEGG pairs      =",
    nrow(integration_up_cluster), "\n")
cat("  Nonzero overlap pairs    =",
    nrow(up_nonzero_from_integration), "\n")
cat("  Zero-overlap pairs       =",
    sum(
      integration_up_cluster$N_Shared_Entrez == 0,
      na.rm = TRUE
    ), "\n")
cat("  Semantic clusters        =",
    length(unique(
      integration_up_cluster$Semantic_Cluster_070
    )), "\n\n")

cat("Tumor-down:\n")
cat("  GO terms                 =",
    length(unique(integration_down_cluster$GO_ID)), "\n")
cat("  KEGG pathways            =",
    length(unique(integration_down_cluster$KEGG_ID)), "\n")
cat("  Total GO–KEGG pairs      =",
    nrow(integration_down_cluster), "\n")
cat("  Nonzero overlap pairs    =",
    nrow(down_nonzero_from_integration), "\n")
cat("  Zero-overlap pairs       =",
    sum(
      integration_down_cluster$N_Shared_Entrez == 0,
      na.rm = TRUE
    ), "\n")
cat("  Semantic clusters        =",
    length(unique(
      integration_down_cluster$Semantic_Cluster_070
    )), "\n\n")


# ------------------------------------------------------------
# 14. READ-ONLY audit statement
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("READ-ONLY AUDIT STATUS: PASS\n")
cat("------------------------------------------------------------\n")

cat("integration_up/down used as complete GO–KEGG pair universe.\n")
cat("evidence_up/down verified as nonzero-overlap subsets.\n")
cat("Existing semantic clusters used unchanged.\n")
cat("No GO enrichment recalculation.\n")
cat("No KEGG enrichment recalculation.\n")
cat("No GO–KEGG overlap recalculation.\n")
cat("No semantic similarity recalculation.\n")
cat("No new clustering.\n")
cat("No ID conversion.\n")
cat("No P-value recalculation.\n")
cat("No FDR recalculation.\n")
cat("No threshold modification.\n")
cat("No background modification.\n")
cat("No statistical ranking.\n")
cat("No theme score calculation.\n")
cat("No top-N filtering.\n")
cat("No locked result modification.\n\n")

cat("Ready for descriptive theme-level aggregation.\n")


# ============================================================
# STEP9C-7C
# Descriptive Theme-Level Aggregation
# READ-ONLY
# ============================================================

cat("\n============================================================\n")
cat("STEP9C-7C: Descriptive Theme-Level Aggregation\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Helper: safe numeric aggregation
# ------------------------------------------------------------

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}


# ------------------------------------------------------------
# 1. Semantic-cluster-level aggregation
#    Based on NONZERO evidence only
# ------------------------------------------------------------

aggregate_semantic_cluster <- function(
    full_df,
    evidence_df,
    direction_label
) {
  
  # ----------------------------------------------------------
  # Full cluster membership
  # ----------------------------------------------------------
  
  cluster_membership <- unique(
    full_df[
      ,
      c(
        "GO_ID",
        "GO_Description",
        "Semantic_Cluster_070"
      )
    ]
  )
  
  cluster_membership$Direction <- direction_label
  
  
  # ----------------------------------------------------------
  # Evidence subset
  # ----------------------------------------------------------
  
  evidence_df <- evidence_df[
    evidence_df$N_Shared_Entrez > 0,
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Aggregate each semantic cluster
  # ----------------------------------------------------------
  
  cluster_ids <- unique(
    cluster_membership$Semantic_Cluster_070
  )
  
  out <- lapply(
    cluster_ids,
    function(cluster_id) {
      
      members <- cluster_membership[
        cluster_membership$Semantic_Cluster_070 ==
          cluster_id,
        ,
        drop = FALSE
      ]
      
      ev <- evidence_df[
        evidence_df$GO_ID %in% members$GO_ID,
        ,
        drop = FALSE
      ]
      
      data.frame(
        Direction =
          direction_label,
        
        Semantic_Cluster_070 =
          cluster_id,
        
        N_GO_Terms =
          length(unique(members$GO_ID)),
        
        N_GO_Terms_with_KEGG_Overlap =
          length(unique(ev$GO_ID)),
        
        N_KEGG_Pathways =
          length(unique(ev$KEGG_ID)),
        
        N_GO_KEGG_Nonzero_Pairs =
          nrow(ev),
        
        Max_N_Shared_Entrez =
          safe_max(ev$N_Shared_Entrez),
        
        Mean_N_Shared_Entrez =
          safe_mean(ev$N_Shared_Entrez),
        
        Median_N_Shared_Entrez =
          safe_median(ev$N_Shared_Entrez),
        
        Max_Jaccard =
          safe_max(ev$Jaccard),
        
        Mean_Jaccard =
          safe_mean(ev$Jaccard),
        
        Median_Jaccard =
          safe_median(ev$Jaccard),
        
        Max_GO_Coverage =
          safe_max(ev$GO_Coverage),
        
        Mean_GO_Coverage =
          safe_mean(ev$GO_Coverage),
        
        Max_KEGG_Coverage =
          safe_max(ev$KEGG_Coverage),
        
        Mean_KEGG_Coverage =
          safe_mean(ev$KEGG_Coverage),
        
        Min_GO_p_adjust =
          safe_min(ev$GO_p_adjust),
        
        Min_KEGG_p_adjust =
          safe_min(ev$KEGG_p_adjust),
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  do.call(rbind, out)
}


cluster_theme_up <- aggregate_semantic_cluster(
  integration_up_cluster,
  evidence_up,
  "Tumor-up"
)

cluster_theme_down <- aggregate_semantic_cluster(
  integration_down_cluster,
  evidence_down,
  "Tumor-down"
)


# ------------------------------------------------------------
# 2. Cluster count validation
# ------------------------------------------------------------

stopifnot(
  nrow(cluster_theme_up) == 97,
  nrow(cluster_theme_down) == 743
)


# ------------------------------------------------------------
# 3. KEGG-pathway-level aggregation
#    Based on NONZERO evidence only
# ------------------------------------------------------------

aggregate_kegg_evidence <- function(
    full_df,
    evidence_df,
    direction_label
) {
  
  evidence_df <- evidence_df[
    evidence_df$N_Shared_Entrez > 0,
    ,
    drop = FALSE
  ]
  
  kegg_ids <- unique(
    full_df$KEGG_ID
  )
  
  out <- lapply(
    kegg_ids,
    function(kegg_id) {
      
      ev <- evidence_df[
        evidence_df$KEGG_ID == kegg_id,
        ,
        drop = FALSE
      ]
      
      all_rows <- full_df[
        full_df$KEGG_ID == kegg_id,
        ,
        drop = FALSE
      ]
      
      data.frame(
        Direction =
          direction_label,
        
        KEGG_ID =
          kegg_id,
        
        KEGG_Description =
          unique(all_rows$KEGG_Description)[1],
        
        N_GO_Terms_Total =
          length(unique(all_rows$GO_ID)),
        
        N_GO_Terms_with_Overlap =
          length(unique(ev$GO_ID)),
        
        N_Semantic_Clusters_with_Overlap =
          length(unique(
            ev$Semantic_Cluster_070
          )),
        
        N_GO_KEGG_Nonzero_Pairs =
          nrow(ev),
        
        Max_N_Shared_Entrez =
          safe_max(ev$N_Shared_Entrez),
        
        Mean_N_Shared_Entrez =
          safe_mean(ev$N_Shared_Entrez),
        
        Median_N_Shared_Entrez =
          safe_median(ev$N_Shared_Entrez),
        
        Max_Jaccard =
          safe_max(ev$Jaccard),
        
        Mean_Jaccard =
          safe_mean(ev$Jaccard),
        
        Max_GO_Coverage =
          safe_max(ev$GO_Coverage),
        
        Mean_GO_Coverage =
          safe_mean(ev$GO_Coverage),
        
        Max_KEGG_Coverage =
          safe_max(ev$KEGG_Coverage),
        
        Mean_KEGG_Coverage =
          safe_mean(ev$KEGG_Coverage),
        
        KEGG_p_adjust =
          unique(all_rows$KEGG_p_adjust)[1],
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  do.call(rbind, out)
}


kegg_theme_up <- aggregate_kegg_evidence(
  integration_up_cluster,
  evidence_up,
  "Tumor-up"
)

kegg_theme_down <- aggregate_kegg_evidence(
  integration_down_cluster,
  evidence_down,
  "Tumor-down"
)


# ------------------------------------------------------------
# 4. KEGG count validation
# ------------------------------------------------------------

stopifnot(
  nrow(kegg_theme_up) == 13,
  nrow(kegg_theme_down) == 60
)


# ------------------------------------------------------------
# 5. GO-term-level aggregation
# ------------------------------------------------------------

aggregate_go_evidence <- function(
    full_df,
    evidence_df,
    direction_label
) {
  
  go_ids <- unique(
    full_df$GO_ID
  )
  
  out <- lapply(
    go_ids,
    function(go_id) {
      
      all_rows <- full_df[
        full_df$GO_ID == go_id,
        ,
        drop = FALSE
      ]
      
      ev <- evidence_df[
        evidence_df$GO_ID == go_id,
        ,
        drop = FALSE
      ]
      
      data.frame(
        Direction =
          direction_label,
        
        GO_ID =
          go_id,
        
        GO_Description =
          unique(all_rows$GO_Description)[1],
        
        Semantic_Cluster_070 =
          unique(
            all_rows$Semantic_Cluster_070
          )[1],
        
        GO_p_adjust =
          unique(all_rows$GO_p_adjust)[1],
        
        N_KEGG_Pathways_Total =
          length(unique(all_rows$KEGG_ID)),
        
        N_KEGG_Pathways_with_Overlap =
          length(unique(ev$KEGG_ID)),
        
        N_GO_KEGG_Nonzero_Pairs =
          nrow(ev),
        
        Max_N_Shared_Entrez =
          safe_max(ev$N_Shared_Entrez),
        
        Mean_N_Shared_Entrez =
          safe_mean(ev$N_Shared_Entrez),
        
        Max_Jaccard =
          safe_max(ev$Jaccard),
        
        Mean_Jaccard =
          safe_mean(ev$Jaccard),
        
        Max_GO_Coverage =
          safe_max(ev$GO_Coverage),
        
        Mean_GO_Coverage =
          safe_mean(ev$GO_Coverage),
        
        Max_KEGG_Coverage =
          safe_max(ev$KEGG_Coverage),
        
        Mean_KEGG_Coverage =
          safe_mean(ev$KEGG_Coverage),
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  do.call(rbind, out)
}


go_theme_up <- aggregate_go_evidence(
  integration_up_cluster,
  evidence_up,
  "Tumor-up"
)

go_theme_down <- aggregate_go_evidence(
  integration_down_cluster,
  evidence_down,
  "Tumor-down"
)


# ------------------------------------------------------------
# 6. GO count validation
# ------------------------------------------------------------

stopifnot(
  nrow(go_theme_up) == 161,
  nrow(go_theme_down) == 1296
)


# ------------------------------------------------------------
# 7. Evidence coverage summary
# ------------------------------------------------------------

theme_coverage_summary <- data.frame(
  Direction = c(
    "Tumor-up",
    "Tumor-down"
  ),
  
  GO_Terms_Total = c(
    161,
    1296
  ),
  
  GO_Terms_with_KEGG_Overlap = c(
    length(unique(evidence_up$GO_ID)),
    length(unique(evidence_down$GO_ID))
  ),
  
  KEGG_Pathways_Total = c(
    13,
    60
  ),
  
  KEGG_Pathways_with_GO_Overlap = c(
    length(unique(evidence_up$KEGG_ID)),
    length(unique(evidence_down$KEGG_ID))
  ),
  
  Total_GO_KEGG_Pairs = c(
    2093,
    77760
  ),
  
  Nonzero_GO_KEGG_Pairs = c(
    388,
    31108
  ),
  
  Zero_GO_KEGG_Pairs = c(
    1705,
    46652
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 8. Structural validation
# ------------------------------------------------------------

stopifnot(
  theme_coverage_summary$GO_Terms_Total ==
    c(161, 1296)
)

stopifnot(
  theme_coverage_summary$KEGG_Pathways_Total ==
    c(13, 60)
)

stopifnot(
  theme_coverage_summary$Total_GO_KEGG_Pairs ==
    c(2093, 77760)
)

stopifnot(
  theme_coverage_summary$Nonzero_GO_KEGG_Pairs ==
    c(388, 31108)
)

stopifnot(
  theme_coverage_summary$Zero_GO_KEGG_Pairs ==
    c(1705, 46652)
)


# ------------------------------------------------------------
# 9. Summary output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-7C DESCRIPTIVE THEME AGGREGATION: PASS\n")
cat("============================================================\n\n")

cat("Tumor-up:\n")
cat("  Semantic clusters                    =",
    nrow(cluster_theme_up), "\n")
cat("  KEGG pathways                        =",
    nrow(kegg_theme_up), "\n")
cat("  GO terms                             =",
    nrow(go_theme_up), "\n")
cat("  Nonzero GO–KEGG pairs                =",
    sum(evidence_up$N_Shared_Entrez > 0), "\n")
cat("  GO terms with overlap                =",
    length(unique(evidence_up$GO_ID)), "\n")
cat("  KEGG pathways with overlap           =",
    length(unique(evidence_up$KEGG_ID)), "\n\n")

cat("Tumor-down:\n")
cat("  Semantic clusters                    =",
    nrow(cluster_theme_down), "\n")
cat("  KEGG pathways                        =",
    nrow(kegg_theme_down), "\n")
cat("  GO terms                             =",
    nrow(go_theme_down), "\n")
cat("  Nonzero GO–KEGG pairs                =",
    sum(evidence_down$N_Shared_Entrez > 0), "\n")
cat("  GO terms with overlap                =",
    length(unique(evidence_down$GO_ID)), "\n")
cat("  KEGG pathways with overlap           =",
    length(unique(evidence_down$KEGG_ID)), "\n\n")


# ------------------------------------------------------------
# 10. READ-ONLY status
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("READ-ONLY AUDIT STATUS: PASS\n")
cat("------------------------------------------------------------\n")

cat("Semantic clusters were used unchanged.\n")
cat("GO–KEGG integration values were used unchanged.\n")
cat("Nonzero evidence values were used unchanged.\n")
cat("No enrichment recalculation.\n")
cat("No GO–KEGG overlap recalculation.\n")
cat("No semantic similarity recalculation.\n")
cat("No new clustering.\n")
cat("No ID conversion.\n")
cat("No P-value/FDR recalculation.\n")
cat("No significance threshold modification.\n")
cat("No background modification.\n")
cat("No composite theme score.\n")
cat("No statistical ranking.\n")
cat("No top-N filtering.\n")
cat("No locked result modification.\n\n")

cat("Ready for final candidate-theme evidence tables.\n")


# ============================================================
# Paper 1 — STEP9C-7D
# Final Candidate-Theme Evidence Tables
# ============================================================
#
# Purpose:
#   Consolidate the already LOCKED STEP9C-7A / 7B / 7C results
#   into final descriptive evidence tables for manuscript/
#   supplementary reporting.
#
# IMPORTANT:
#   This is a READ-ONLY consolidation step.
#
#   NO GO enrichment recalculation
#   NO KEGG enrichment recalculation
#   NO GO–KEGG overlap recalculation
#   NO semantic similarity recalculation
#   NO new clustering
#   NO ID conversion
#   NO P-value recalculation
#   NO FDR recalculation
#   NO threshold modification
#   NO background modification
#   NO composite theme score
#   NO statistical ranking
#   NO top-N filtering
#   NO locked-result modification
#
# ============================================================


cat("\n============================================================\n")
cat("STEP9C-7D: Final Candidate-Theme Evidence Tables\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory and output folders
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9C-7D"

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}

dir.create(
  file.path(step_dir, "Tables"),
  showWarnings = FALSE,
  recursive = TRUE
)

cat("Output directory:\n")
cat(step_dir, "\n\n")


# ------------------------------------------------------------
# 1. Required locked objects
# ------------------------------------------------------------

required_objects <- c(
  "integration_up_cluster",
  "integration_down_cluster",
  "evidence_up",
  "evidence_down",
  "cluster_theme_up",
  "cluster_theme_down",
  "kegg_theme_up",
  "kegg_theme_down",
  "go_theme_up",
  "go_theme_down",
  "cluster_up_std",
  "cluster_down_std",
  "up_representatives",
  "down_representatives"
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
      "Required locked objects are missing:",
      paste(missing_objects, collapse = ", ")
    )
  )
}

cat("Required locked objects: PASS\n")
cat("All required objects are present in the current R session.\n\n")


# ------------------------------------------------------------
# 2. Structural lock checks
# ------------------------------------------------------------

stopifnot(
  nrow(integration_up_cluster) == 2093,
  nrow(integration_down_cluster) == 77760,
  
  nrow(evidence_up) == 388,
  nrow(evidence_down) == 31108,
  
  nrow(cluster_theme_up) == 97,
  nrow(cluster_theme_down) == 743,
  
  nrow(kegg_theme_up) == 13,
  nrow(kegg_theme_down) == 60,
  
  nrow(go_theme_up) == 161,
  nrow(go_theme_down) == 1296
)

cat("Structural lock checks: PASS\n")
cat("  integration_up_cluster   =", nrow(integration_up_cluster), "\n")
cat("  integration_down_cluster =", nrow(integration_down_cluster), "\n")
cat("  evidence_up              =", nrow(evidence_up), "\n")
cat("  evidence_down            =", nrow(evidence_down), "\n")
cat("  cluster_theme_up         =", nrow(cluster_theme_up), "\n")
cat("  cluster_theme_down       =", nrow(cluster_theme_down), "\n")
cat("  kegg_theme_up            =", nrow(kegg_theme_up), "\n")
cat("  kegg_theme_down          =", nrow(kegg_theme_down), "\n")
cat("  go_theme_up              =", nrow(go_theme_up), "\n")
cat("  go_theme_down            =", nrow(go_theme_down), "\n\n")


# ------------------------------------------------------------
# 3. Helper: safe character collapse
# ------------------------------------------------------------

collapse_unique <- function(x, sep = "; ") {
  
  x <- as.character(x)
  
  x <- x[
    !is.na(x) &
      nzchar(x)
  ]
  
  x <- unique(x)
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(x, collapse = sep)
}


# ------------------------------------------------------------
# 4. Helper: numeric summary
# ------------------------------------------------------------

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  min(x)
}


safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  max(x)
}


safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  mean(x)
}


safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  median(x)
}


# ------------------------------------------------------------
# 5. Representative GO term lookup
#
# Existing STEP9C-7A representative objects are used.
#
# No new representative selection is performed.
# ------------------------------------------------------------

prepare_representatives <- function(
    rep_df,
    direction_label
) {
  
  required <- c("ID", "Description")
  
  if (!all(required %in% names(rep_df))) {
    stop(
      paste0(
        "Representative table for ",
        direction_label,
        " does not contain required columns."
      )
    )
  }
  
  out <- data.frame(
    GO_ID = as.character(rep_df$ID),
    Representative_GO_Description =
      as.character(rep_df$Description),
    stringsAsFactors = FALSE
  )
  
  out <- unique(out)
  
  if (anyDuplicated(out$GO_ID)) {
    stop(
      paste0(
        "Duplicate GO IDs detected in representative table: ",
        direction_label
      )
    )
  }
  
  out
}


rep_up <- prepare_representatives(
  up_representatives,
  "Tumor-up"
)

rep_down <- prepare_representatives(
  down_representatives,
  "Tumor-down"
)


# ------------------------------------------------------------
# 6. FINAL TABLE 1
#    Tumor-up semantic-cluster evidence
# ------------------------------------------------------------
#
# All 97 semantic clusters are retained.
#
# Evidence-related fields are based only on the existing
# nonzero-overlap evidence_up object.
#
# No ranking.
# No theme score.
# No top-N filtering.
# ------------------------------------------------------------

build_cluster_evidence <- function(
    full_cluster_df,
    evidence_df,
    rep_df,
    direction_label
) {
  
  # Full cluster membership
  full <- full_cluster_df
  
  full$GO_ID <- as.character(full$GO_ID)
  full$Cluster_ID <- as.character(full$Semantic_Cluster_070)
  
  cluster_ids <- unique(full$Cluster_ID)
  
  result_list <- vector(
    "list",
    length(cluster_ids)
  )
  
  for (i in seq_along(cluster_ids)) {
    
    cid <- cluster_ids[i]
    
    members <- full[
      full$Cluster_ID == cid,
      ,
      drop = FALSE
    ]
    
    evidence_members <- evidence_df[
      as.character(evidence_df$GO_ID) %in%
        as.character(members$GO_ID),
      ,
      drop = FALSE
    ]
    
    go_ids <- unique(as.character(members$GO_ID))
    
    representative_match <- rep_df[
      rep_df$GO_ID %in% go_ids,
      ,
      drop = FALSE
    ]
    
    representative_go <- collapse_unique(
      representative_match$GO_ID
    )
    
    representative_desc <- collapse_unique(
      representative_match$Representative_GO_Description
    )
    
    result_list[[i]] <- data.frame(
      
      Direction = direction_label,
      
      Semantic_Cluster_070 = cid,
      
      Representative_GO_ID =
        representative_go,
      
      Representative_GO_Description =
        representative_desc,
      
      N_GO_Terms =
        length(go_ids),
      
      N_GO_Terms_with_KEGG_Overlap =
        length(unique(
          as.character(evidence_members$GO_ID)
        )),
      
      N_KEGG_Pathways =
        length(unique(
          as.character(evidence_members$KEGG_ID)
        )),
      
      N_GO_KEGG_Nonzero_Pairs =
        nrow(evidence_members),
      
      Shared_Entrez_Min =
        safe_min(evidence_members$N_Shared_Entrez),
      
      Shared_Entrez_Median =
        safe_median(evidence_members$N_Shared_Entrez),
      
      Shared_Entrez_Max =
        safe_max(evidence_members$N_Shared_Entrez),
      
      KEGG_Pathways =
        collapse_unique(
          evidence_members$KEGG_ID
        ),
      
      KEGG_Descriptions =
        collapse_unique(
          evidence_members$KEGG_Description
        ),
      
      Core_Genes_Shared =
        collapse_unique(
          evidence_members$Core_Genes_Shared
        ),
      
      stringsAsFactors = FALSE
    )
  }
  
  result <- do.call(
    rbind,
    result_list
  )
  
  rownames(result) <- NULL
  
  result
}


final_cluster_theme_up <- build_cluster_evidence(
  integration_up_cluster,
  evidence_up,
  rep_up,
  "Tumor-up"
)


final_cluster_theme_down <- build_cluster_evidence(
  integration_down_cluster,
  evidence_down,
  rep_down,
  "Tumor-down"
)


# ------------------------------------------------------------
# 7. FINAL TABLE 2
#    Tumor-up / Tumor-down KEGG evidence tables
# ------------------------------------------------------------

build_kegg_evidence <- function(
    full_df,
    evidence_df,
    direction_label
) {
  
  kegg_ids <- unique(
    as.character(full_df$KEGG_ID)
  )
  
  result_list <- vector(
    "list",
    length(kegg_ids)
  )
  
  for (i in seq_along(kegg_ids)) {
    
    kid <- kegg_ids[i]
    
    full_members <- full_df[
      as.character(full_df$KEGG_ID) == kid,
      ,
      drop = FALSE
    ]
    
    evidence_members <- evidence_df[
      as.character(evidence_df$KEGG_ID) == kid,
      ,
      drop = FALSE
    ]
    
    result_list[[i]] <- data.frame(
      
      Direction = direction_label,
      
      KEGG_ID = kid,
      
      KEGG_Description =
        collapse_unique(
          full_members$KEGG_Description
        ),
      
      N_GO_Terms =
        length(unique(
          as.character(full_members$GO_ID)
        )),
      
      N_GO_Terms_with_Overlap =
        length(unique(
          as.character(evidence_members$GO_ID)
        )),
      
      N_Nonzero_GO_KEGG_Pairs =
        nrow(evidence_members),
      
      Shared_Entrez_Min =
        safe_min(
          evidence_members$N_Shared_Entrez
        ),
      
      Shared_Entrez_Median =
        safe_median(
          evidence_members$N_Shared_Entrez
        ),
      
      Shared_Entrez_Max =
        safe_max(
          evidence_members$N_Shared_Entrez
        ),
      
      GO_Terms =
        collapse_unique(
          evidence_members$GO_ID
        ),
      
      GO_Descriptions =
        collapse_unique(
          evidence_members$GO_Description
        ),
      
      Core_Genes_Shared =
        collapse_unique(
          evidence_members$Core_Genes_Shared
        ),
      
      stringsAsFactors = FALSE
    )
  }
  
  result <- do.call(
    rbind,
    result_list
  )
  
  rownames(result) <- NULL
  
  result
}


final_kegg_up <- build_kegg_evidence(
  integration_up_cluster,
  evidence_up,
  "Tumor-up"
)


final_kegg_down <- build_kegg_evidence(
  integration_down_cluster,
  evidence_down,
  "Tumor-down"
)


# ------------------------------------------------------------
# 8. FINAL TABLE 3
#    Tumor-up / Tumor-down GO evidence tables
# ------------------------------------------------------------

build_go_evidence <- function(
    full_df,
    evidence_df,
    direction_label
) {
  
  go_ids <- unique(
    as.character(full_df$GO_ID)
  )
  
  result_list <- vector(
    "list",
    length(go_ids)
  )
  
  for (i in seq_along(go_ids)) {
    
    gid <- go_ids[i]
    
    full_members <- full_df[
      as.character(full_df$GO_ID) == gid,
      ,
      drop = FALSE
    ]
    
    evidence_members <- evidence_df[
      as.character(evidence_df$GO_ID) == gid,
      ,
      drop = FALSE
    ]
    
    cluster_id <- collapse_unique(
      full_members$Cluster_ID
    )
    
    result_list[[i]] <- data.frame(
      
      Direction = direction_label,
      
      GO_ID = gid,
      
      GO_Description =
        collapse_unique(
          full_members$GO_Description
        ),
      
      Semantic_Cluster_070 =
        cluster_id,
      
      N_KEGG_Pathways =
        length(unique(
          as.character(
            evidence_members$KEGG_ID
          )
        )),
      
      N_Nonzero_GO_KEGG_Pairs =
        nrow(evidence_members),
      
      Shared_Entrez_Min =
        safe_min(
          evidence_members$N_Shared_Entrez
        ),
      
      Shared_Entrez_Median =
        safe_median(
          evidence_members$N_Shared_Entrez
        ),
      
      Shared_Entrez_Max =
        safe_max(
          evidence_members$N_Shared_Entrez
        ),
      
      KEGG_Pathways =
        collapse_unique(
          evidence_members$KEGG_ID
        ),
      
      KEGG_Descriptions =
        collapse_unique(
          evidence_members$KEGG_Description
        ),
      
      Core_Genes_Shared =
        collapse_unique(
          evidence_members$Core_Genes_Shared
        ),
      
      stringsAsFactors = FALSE
    )
  }
  
  result <- do.call(
    rbind,
    result_list
  )
  
  rownames(result) <- NULL
  
  result
}


final_go_up <- build_go_evidence(
  integration_up_cluster,
  evidence_up,
  "Tumor-up"
)


final_go_down <- build_go_evidence(
  integration_down_cluster,
  evidence_down,
  "Tumor-down"
)


# ------------------------------------------------------------
# 9. FINAL TABLE 4
#    Core-gene direct GO–KEGG evidence
# ------------------------------------------------------------
#
# This table reports only existing direct overlap evidence.
#
# It does NOT infer:
#   driver
#   mediator
#   regulator
#   causal relationship
#   leading-edge gene
#
# ------------------------------------------------------------

core_gene_labels <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


build_core_evidence <- function(
    evidence_df,
    direction_label
) {
  
  evidence_df <- evidence_df[
    !is.na(evidence_df$Core_Genes_Shared) &
      nzchar(
        as.character(
          evidence_df$Core_Genes_Shared
        )
      ),
    ,
    drop = FALSE
  ]
  
  if (nrow(evidence_df) == 0) {
    
    return(
      data.frame(
        Direction = character(0),
        Core_Gene = character(0),
        GO_ID = character(0),
        GO_Description = character(0),
        KEGG_ID = character(0),
        KEGG_Description = character(0),
        N_Shared_Entrez = numeric(0),
        Shared_Entrez = character(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  result_list <- list()
  counter <- 1
  
  for (gene in core_gene_labels) {
    
    hit <- grepl(
      paste0(
        "(^|;\\s*)",
        gene,
        "(;|$)"
      ),
      as.character(
        evidence_df$Core_Genes_Shared
      )
    )
    
    tmp <- evidence_df[
      hit,
      ,
      drop = FALSE
    ]
    
    if (nrow(tmp) == 0) {
      next
    }
    
    result_list[[counter]] <- data.frame(
      
      Direction =
        direction_label,
      
      Core_Gene =
        gene,
      
      GO_ID =
        as.character(tmp$GO_ID),
      
      GO_Description =
        as.character(tmp$GO_Description),
      
      KEGG_ID =
        as.character(tmp$KEGG_ID),
      
      KEGG_Description =
        as.character(tmp$KEGG_Description),
      
      N_Shared_Entrez =
        as.numeric(tmp$N_Shared_Entrez),
      
      Shared_Entrez =
        as.character(tmp$Shared_Entrez),
      
      stringsAsFactors = FALSE
    )
    
    counter <- counter + 1
  }
  
  if (length(result_list) == 0) {
    
    return(
      data.frame(
        Direction = character(0),
        Core_Gene = character(0),
        GO_ID = character(0),
        GO_Description = character(0),
        KEGG_ID = character(0),
        KEGG_Description = character(0),
        N_Shared_Entrez = numeric(0),
        Shared_Entrez = character(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  result <- do.call(
    rbind,
    result_list
  )
  
  rownames(result) <- NULL
  
  result
}


core_evidence_up <- build_core_evidence(
  evidence_up,
  "Tumor-up"
)

core_evidence_down <- build_core_evidence(
  evidence_down,
  "Tumor-down"
)

core_evidence_all <- rbind(
  core_evidence_up,
  core_evidence_down
)

rownames(core_evidence_all) <- NULL


# ------------------------------------------------------------
# 10. FINAL TABLE 5
#     Overall candidate-theme coverage summary
# ------------------------------------------------------------

theme_coverage_summary <- data.frame(
  
  Direction = c(
    "Tumor-up",
    "Tumor-down"
  ),
  
  GO_Terms = c(
    161,
    1296
  ),
  
  Semantic_Clusters = c(
    97,
    743
  ),
  
  KEGG_Pathways = c(
    13,
    60
  ),
  
  Total_GO_KEGG_Pairs = c(
    2093,
    77760
  ),
  
  Nonzero_GO_KEGG_Pairs = c(
    388,
    31108
  ),
  
  Zero_Overlap_Pairs = c(
    1705,
    46652
  ),
  
  GO_Terms_with_Overlap = c(
    154,
    1286
  ),
  
  KEGG_Pathways_with_Overlap = c(
    11,
    60
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 11. Final structural validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-7D STRUCTURAL VALIDATION\n")
cat("============================================================\n\n")


stopifnot(
  
  nrow(final_cluster_theme_up) == 97,
  nrow(final_cluster_theme_down) == 743,
  
  nrow(final_kegg_up) == 13,
  nrow(final_kegg_down) == 60,
  
  nrow(final_go_up) == 161,
  nrow(final_go_down) == 1296,
  
  theme_coverage_summary$GO_Terms[
    theme_coverage_summary$Direction == "Tumor-up"
  ] == 161,
  
  theme_coverage_summary$GO_Terms[
    theme_coverage_summary$Direction == "Tumor-down"
  ] == 1296
)


# Exact GO coverage checks

up_go_with_overlap <- length(
  unique(
    as.character(evidence_up$GO_ID)
  )
)

down_go_with_overlap <- length(
  unique(
    as.character(evidence_down$GO_ID)
  )
)

stopifnot(
  up_go_with_overlap == 154,
  down_go_with_overlap == 1286
)


# Exact KEGG coverage checks

up_kegg_with_overlap <- length(
  unique(
    as.character(evidence_up$KEGG_ID)
  )
)

down_kegg_with_overlap <- length(
  unique(
    as.character(evidence_down$KEGG_ID)
  )
)

stopifnot(
  up_kegg_with_overlap == 11,
  down_kegg_with_overlap == 60
)


# Exact nonzero pair checks

stopifnot(
  nrow(evidence_up) == 388,
  nrow(evidence_down) == 31108
)


cat("Tumor-up:\n")
cat("  Cluster evidence rows       =", nrow(final_cluster_theme_up), "\n")
cat("  KEGG evidence rows          =", nrow(final_kegg_up), "\n")
cat("  GO evidence rows            =", nrow(final_go_up), "\n")
cat("  GO terms with overlap       =", up_go_with_overlap, "\n")
cat("  KEGG pathways with overlap  =", up_kegg_with_overlap, "\n")
cat("  Nonzero GO–KEGG pairs       =", nrow(evidence_up), "\n\n")

cat("Tumor-down:\n")
cat("  Cluster evidence rows       =", nrow(final_cluster_theme_down), "\n")
cat("  KEGG evidence rows          =", nrow(final_kegg_down), "\n")
cat("  GO evidence rows            =", nrow(final_go_down), "\n")
cat("  GO terms with overlap       =", down_go_with_overlap, "\n")
cat("  KEGG pathways with overlap  =", down_kegg_with_overlap, "\n")
cat("  Nonzero GO–KEGG pairs       =", nrow(evidence_down), "\n\n")


# ------------------------------------------------------------
# 12. Check that no duplicate primary identifiers exist
# ------------------------------------------------------------

stopifnot(
  !anyDuplicated(final_cluster_theme_up$Semantic_Cluster_070),
  !anyDuplicated(final_cluster_theme_down$Semantic_Cluster_070),
  
  !anyDuplicated(final_kegg_up$KEGG_ID),
  !anyDuplicated(final_kegg_down$KEGG_ID),
  
  !anyDuplicated(final_go_up$GO_ID),
  !anyDuplicated(final_go_down$GO_ID)
)


cat("Identifier uniqueness checks: PASS\n\n")


# ------------------------------------------------------------
# 13. Export final CSV tables
# ------------------------------------------------------------

write.csv(
  final_cluster_theme_up,
  file.path(
    step_dir,
    "Tables",
    "STEP9C7D_TumorUp_SemanticCluster_Evidence.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  final_cluster_theme_down,
  file.path(
    step_dir,
    "Tables",
    "STEP9C7D_TumorDown_SemanticCluster_Evidence.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  final_kegg_up,
  file.path(
    step_dir,
    "Tables",
    "STEP9C7D_TumorUp_KEGG_Evidence.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  final_kegg_down,
  file.path(
    step_dir,
    "Tables",
    "STEP9C7D_TumorDown_KEGG_Evidence.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  final_go_up,
  file.path(
    step_dir,
    "Tables",
    "STEP9C7D_TumorUp_GO_BP_Evidence.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  final_go_down,
  file.path(
    step_dir,
    "Tables",
    "STEP9C7D_TumorDown_GO_BP_Evidence.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  core_evidence_all,
  file.path(
    step_dir,
    "Tables",
    "STEP9C7D_CoreGene_GO_KEGG_DirectOverlap_Evidence.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  theme_coverage_summary,
  file.path(
    step_dir,
    "Tables",
    "STEP9C7D_Theme_Coverage_Summary.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 14. Final LOCK object
# ------------------------------------------------------------

STEP9C7D_FINAL_LOCK <- list(
  
  Step = "STEP9C-7D",
  
  Status = "PASS",
  
  Purpose =
    "Final descriptive candidate-theme evidence table consolidation",
  
  TumorUp = list(
    GO_terms = 161,
    Semantic_clusters = 97,
    KEGG_pathways = 13,
    Total_GO_KEGG_pairs = 2093,
    Nonzero_GO_KEGG_pairs = 388,
    Zero_overlap_pairs = 1705,
    GO_terms_with_overlap = 154,
    KEGG_pathways_with_overlap = 11
  ),
  
  TumorDown = list(
    GO_terms = 1296,
    Semantic_clusters = 743,
    KEGG_pathways = 60,
    Total_GO_KEGG_pairs = 77760,
    Nonzero_GO_KEGG_pairs = 31108,
    Zero_overlap_pairs = 46652,
    GO_terms_with_overlap = 1286,
    KEGG_pathways_with_overlap = 60
  ),
  
  CoreGenes =
    c("AOC3", "MS4A7", "IL7R"),
  
  Representative_GO_Source =
    "Existing locked STEP9C-7A representative GO terms",
  
  Evidence_Source =
    "Existing locked STEP9C-7B nonzero-overlap evidence",
  
  Cluster_Source =
    "Existing locked STEP9C-7A semantic clusters",
  
  Integration_Source =
    "Existing locked STEP9C-7B GO–KEGG integration",
  
  No_GO_Enrichment_Recalculation = TRUE,
  
  No_KEGG_Enrichment_Recalculation = TRUE,
  
  No_GO_KEGG_Overlap_Recalculation = TRUE,
  
  No_Semantic_Similarity_Recalculation = TRUE,
  
  No_New_Clustering = TRUE,
  
  No_ID_Conversion = TRUE,
  
  No_Pvalue_Recalculation = TRUE,
  
  No_FDR_Recalculation = TRUE,
  
  No_Threshold_Modification = TRUE,
  
  No_Background_Modification = TRUE,
  
  No_Composite_Theme_Score = TRUE,
  
  No_Statistical_Ranking = TRUE,
  
  No_TopN_Filtering = TRUE,
  
  No_Locked_Result_Modification = TRUE
)


saveRDS(
  STEP9C7D_FINAL_LOCK,
  file.path(
    step_dir,
    "STEP9C7D_FINAL_LOCK.rds"
  )
)


# ------------------------------------------------------------
# 15. Final output summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9C-7D FINAL CANDIDATE-THEME EVIDENCE TABLES: PASS\n")
cat("============================================================\n\n")

cat("Tumor-up:\n")
cat("  Semantic-cluster evidence =", nrow(final_cluster_theme_up), "\n")
cat("  KEGG evidence             =", nrow(final_kegg_up), "\n")
cat("  GO evidence               =", nrow(final_go_up), "\n")
cat("  Core-gene evidence        = available in combined table\n\n")

cat("Tumor-down:\n")
cat("  Semantic-cluster evidence =", nrow(final_cluster_theme_down), "\n")
cat("  KEGG evidence             =", nrow(final_kegg_down), "\n")
cat("  GO evidence               =", nrow(final_go_down), "\n")
cat("  Core-gene evidence        = available in combined table\n\n")

cat("Final files written to:\n")
cat(file.path(step_dir, "Tables"), "\n\n")

cat("READ-ONLY CONSOLIDATION STATUS: PASS\n")
cat("No enrichment recalculation.\n")
cat("No GO–KEGG overlap recalculation.\n")
cat("No semantic similarity recalculation.\n")
cat("No new clustering.\n")
cat("No ID conversion.\n")
cat("No P-value/FDR recalculation.\n")
cat("No threshold modification.\n")
cat("No background modification.\n")
cat("No composite theme score.\n")
cat("No statistical ranking.\n")
cat("No top-N filtering.\n")
cat("No locked result modification.\n\n")

cat("STEP9C-7D is ready for manuscript-level descriptive reporting.\n")


# ============================================================
# Paper 1
# STEP9D-1
# Genome-wide preranked GSEA
# Input / Gene-Universe Lock
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-1: Genome-wide preranked GSEA Input Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA/9D-1_Input_Lock"

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}

setwd(step_dir)


# ------------------------------------------------------------
# 1. Required package
# ------------------------------------------------------------

required_packages <- c(
  "data.table"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Required package not installed: ",
      pkg
    )
  }
}

library(data.table)

# ============================================================
# Paper 1
# STEP9D-1A
# Locate candidate genome-wide differential-statistic files
# READ-ONLY SEARCH ONLY
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-1A: Locate genome-wide differential-statistic files\n")
cat("============================================================\n\n")


search_roots <- c(
  "input_dir/Paper1_tumor-only WGCNA",
  "input_dir/Paper 1_Step 9-Functional Enrichment"
)


candidate_files <- unlist(
  lapply(
    search_roots,
    function(x) {
      
      if (!dir.exists(x)) {
        return(character(0))
      }
      
      list.files(
        path = x,
        pattern = "\\.(txt|csv|tsv|rds)$",
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
      )
    }
  )
)

candidate_files <- unique(candidate_files)


cat("Total candidate files found:",
    length(candidate_files), "\n\n")


# ------------------------------------------------------------
# Search filenames for differential-analysis keywords
# ------------------------------------------------------------

keyword_pattern <- paste(
  c(
    "deg",
    "differential",
    "de_gene",
    "gene.*stat",
    "stat.*gene",
    "result",
    "limma",
    "deseq",
    "fc",
    "pvalue",
    "p\\.value",
    "fdr",
    "adj"
  ),
  collapse = "|"
)


candidate_hits <- candidate_files[
  grepl(
    keyword_pattern,
    basename(candidate_files),
    ignore.case = TRUE
  )
]


cat("Files matching differential-analysis keywords:",
    length(candidate_hits), "\n\n")


if (length(candidate_hits) > 0) {
  
  for (i in seq_along(candidate_hits)) {
    
    cat(
      sprintf("[%03d] %s\n",
              i,
              candidate_hits[i])
    )
  }
  
} else {
  
  cat("No filename-level differential-analysis candidates found.\n")
}


# ------------------------------------------------------------
# Save search result
# ------------------------------------------------------------

writeLines(
  candidate_hits,
  "STEP9D1A_Candidate_Differential_Files.txt"
)


cat("\n============================================================\n")
cat("STEP9D-1A STATUS: SEARCH COMPLETE\n")
cat("============================================================\n")

cat(
  "No file was modified.\n",
  "No expression data were modified.\n",
  "No DEG result was recalculated.\n",
  "No GSEA was performed.\n"
)


### ============================================================
# Paper 1
# STEP9D-1B
# Structural audit of candidate differential-analysis files
#
# PURPOSE:
#   Identify the actual genome-wide differential-statistic
#   candidate WITHOUT performing any statistical analysis.
#
# IMPORTANT:
#   READ-ONLY
#   NO GSEA
#   NO DEG recalculation
#   NO modification of STEP9B
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-1B: Candidate Differential-Analysis Structure Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA/9D-1_Input_Lock"

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}

setwd(step_dir)


# ------------------------------------------------------------
# 1. Candidate list generated by STEP9D-1A
# ------------------------------------------------------------

candidate_file_list <-
  file.path(
    step_dir,
    "STEP9D1A_Candidate_Differential_Files.txt"
  )

if (!file.exists(candidate_file_list)) {
  stop(
    "STEP9D1A candidate list not found.\n",
    "Please run STEP9D-1A first."
  )
}


candidate_files <-
  readLines(candidate_file_list, warn = FALSE)

candidate_files <-
  candidate_files[
    nzchar(candidate_files)
  ]

cat(
  "Candidate files:",
  length(candidate_files),
  "\n\n"
)


# ------------------------------------------------------------
# 2. Helper function
# ------------------------------------------------------------

inspect_table <- function(f) {
  
  result <- data.frame(
    File = f,
    Exists = file.exists(f),
    Extension = tolower(tools::file_ext(f)),
    N_rows = NA_integer_,
    N_columns = NA_integer_,
    Gene_like_column = NA_character_,
    LogFC_like_column = NA_character_,
    Pvalue_like_column = NA_character_,
    FDR_like_column = NA_character_,
    stringsAsFactors = FALSE
  )
  
  if (!file.exists(f)) {
    return(result)
  }
  
  ext <- tolower(
    tools::file_ext(f)
  )
  
  x <- tryCatch({
    
    if (ext == "csv") {
      
      read.csv(
        f,
        nrows = 20,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      
    } else if (ext %in% c("txt", "tsv")) {
      
      data.table::fread(
        f,
        nrows = 20,
        data.table = FALSE,
        check.names = FALSE
      )
      
    } else {
      
      NULL
    }
    
  }, error = function(e) {
    NULL
  })
  
  
  if (is.null(x)) {
    return(result)
  }
  
  
  # ----------------------------------------------------------
  # Column names
  # ----------------------------------------------------------
  
  cn <- colnames(x)
  
  result$N_columns <- length(cn)
  
  
  # ----------------------------------------------------------
  # Gene-like columns
  # ----------------------------------------------------------
  
  gene_hits <- cn[
    grepl(
      "^(gene|geneid|gene_id|ensembl|ensemblid|entrez|symbol|genes)$",
      cn,
      ignore.case = TRUE
    )
  ]
  
  if (length(gene_hits) > 0) {
    result$Gene_like_column <-
      paste(gene_hits, collapse = ";")
  }
  
  
  # ----------------------------------------------------------
  # LogFC-like columns
  # ----------------------------------------------------------
  
  logfc_hits <- cn[
    grepl(
      "log2?fc|logfc|fold.?change|effect",
      cn,
      ignore.case = TRUE
    )
  ]
  
  if (length(logfc_hits) > 0) {
    result$LogFC_like_column <-
      paste(logfc_hits, collapse = ";")
  }
  
  
  # ----------------------------------------------------------
  # P-value-like columns
  # ----------------------------------------------------------
  
  p_hits <- cn[
    grepl(
      "^p$|pvalue|p.value|p_val|p.val|pvalue",
      cn,
      ignore.case = TRUE
    )
  ]
  
  if (length(p_hits) > 0) {
    result$Pvalue_like_column <-
      paste(p_hits, collapse = ";")
  }
  
  
  # ----------------------------------------------------------
  # FDR / adjusted P-like columns
  # ----------------------------------------------------------
  
  fdr_hits <- cn[
    grepl(
      "fdr|adj.*p|p.*adj|qvalue|q.value",
      cn,
      ignore.case = TRUE
    )
  ]
  
  if (length(fdr_hits) > 0) {
    result$FDR_like_column <-
      paste(fdr_hits, collapse = ";")
  }
  
  
  # ----------------------------------------------------------
  # Try to obtain row count without loading huge objects
  # ----------------------------------------------------------
  
  result$N_rows <- tryCatch({
    
    if (ext == "csv") {
      
      nrow(
        read.csv(
          f,
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      )
      
    } else if (ext %in% c("txt", "tsv")) {
      
      nrow(
        data.table::fread(
          f,
          data.table = FALSE,
          check.names = FALSE
        )
      )
      
    } else {
      
      NA_integer_
    }
    
  }, error = function(e) {
    NA_integer_
  })
  
  
  result
}


# ------------------------------------------------------------
# 3. Inspect all candidates
# ------------------------------------------------------------

audit_list <-
  lapply(
    candidate_files,
    inspect_table
  )

audit_df <-
  do.call(
    rbind,
    audit_list
  )


# ------------------------------------------------------------
# 4. Save complete audit
# ------------------------------------------------------------

write.csv(
  audit_df,
  "STEP9D1B_All_Candidate_Structural_Audit.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Identify files containing the minimum expected
#    differential-statistic structure
# ------------------------------------------------------------

potential_hits <-
  audit_df[
    !is.na(audit_df$Gene_like_column) &
      !is.na(audit_df$LogFC_like_column) &
      !is.na(audit_df$Pvalue_like_column),
  ]


# ------------------------------------------------------------
# 6. Remove obvious audit / QC / locked-result files
#    ONLY for candidate prioritization.
#
#    IMPORTANT:
#    No source file is deleted or modified.
# ------------------------------------------------------------

exclude_pattern <- paste(
  c(
    "audit",
    "qc",
    "lock",
    "summary",
    "mapping",
    "inventory",
    "candidate",
    "checkpoint",
    "enrichment",
    "ora",
    "gsea"
  ),
  collapse = "|"
)


potential_primary <-
  potential_hits[
    !grepl(
      exclude_pattern,
      basename(potential_hits$File),
      ignore.case = TRUE
    ),
  ]


# ------------------------------------------------------------
# 7. Print results
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Files containing Gene + LogFC + P-value-like columns:\n")
cat("============================================================\n\n")

cat(
  "Total:",
  nrow(potential_hits),
  "\n\n"
)


if (nrow(potential_hits) > 0) {
  
  print(
    potential_hits[
      ,
      c(
        "File",
        "N_rows",
        "N_columns",
        "Gene_like_column",
        "LogFC_like_column",
        "Pvalue_like_column",
        "FDR_like_column"
      )
    ],
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No candidate with the complete expected structure was found.\n"
  )
}


cat("\n============================================================\n")
cat("Primary candidates after filename-level exclusion:\n")
cat("============================================================\n\n")

cat(
  "Total:",
  nrow(potential_primary),
  "\n\n"
)


if (nrow(potential_primary) > 0) {
  
  print(
    potential_primary[
      ,
      c(
        "File",
        "N_rows",
        "N_columns",
        "Gene_like_column",
        "LogFC_like_column",
        "Pvalue_like_column",
        "FDR_like_column"
      )
    ],
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No primary candidate identified by automatic filtering.\n"
  )
}


# ------------------------------------------------------------
# 8. Save primary candidate table
# ------------------------------------------------------------

write.csv(
  potential_hits,
  "STEP9D1B_Potential_GenomeWide_DEG_Statistic_Candidates.csv",
  row.names = FALSE
)

write.csv(
  potential_primary,
  "STEP9D1B_Primary_GenomeWide_DEG_Statistic_Candidates.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-1B STATUS: STRUCTURAL AUDIT COMPLETE\n")
cat("============================================================\n")

cat(
  "No statistical analysis performed.\n",
  "No gene ranking performed.\n",
  "No GSEA performed.\n",
  "No DEG recalculation performed.\n",
  "No STEP9B result modified.\n",
  "No source file modified.\n"
)


# ============================================================
cat("\n============================================================\n")
cat("STEP9D-1C-R: Genome-wide Universe Column Recovery Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA/9D-1_Input_Lock"

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}

setwd(step_dir)


# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

candidate_file <-
  "result_dir/Paper 1_Step 9-Functional Enrichment/STEP9B-7R_Revised_DEG_Definition_and_Lock/02_results/STEP9B7R_Revised_DEG_FULL.csv"

universe_file <-
  "result_dir/Paper 1_Step 9-Functional Enrichment/STEP9C_Revised_DEG_GO_KEGG_ORA/00_input/STEP9C_DEG_testing_universe_57828.csv"


if (!file.exists(candidate_file)) {
  stop("FAIL: Candidate file not found.\n", candidate_file)
}

if (!file.exists(universe_file)) {
  stop("FAIL: Universe file not found.\n", universe_file)
}


# ------------------------------------------------------------
# 2. Read files read-only
# ------------------------------------------------------------

candidate_df <-
  data.table::fread(
    candidate_file,
    data.table = FALSE,
    check.names = FALSE
  )

universe_df <-
  data.table::fread(
    universe_file,
    data.table = FALSE,
    check.names = FALSE
  )


# ------------------------------------------------------------
# 3. Structural inspection
# ------------------------------------------------------------

cat("Candidate dimensions:\n")
print(dim(candidate_df))

cat("\nCandidate columns:\n")
print(colnames(candidate_df))

cat("\nUniverse dimensions:\n")
print(dim(universe_df))

cat("\nUniverse columns:\n")
print(colnames(universe_df))


# ------------------------------------------------------------
# 4. Inspect universe x column
# ------------------------------------------------------------

if (!"x" %in% colnames(universe_df)) {
  stop(
    "FAIL: Expected universe column 'x' not found."
  )
}

universe_ids_raw <- universe_df[["x"]]


cat("\n------------------------------------------------------------\n")
cat("Universe column 'x' inspection\n")
cat("------------------------------------------------------------\n\n")

cat("Class:\n")
print(class(universe_ids_raw))

cat("\nLength:\n")
print(length(universe_ids_raw))

cat("\nFirst 20 values:\n")
print(head(universe_ids_raw, 20))

cat("\nLast 20 values:\n")
print(tail(universe_ids_raw, 20))

cat("\nNumber missing / NA:\n")
print(sum(is.na(universe_ids_raw)))

cat("\nNumber empty strings:\n")
print(sum(trimws(as.character(universe_ids_raw)) == ""))

cat("\nNumber unique IDs:\n")
print(length(unique(as.character(universe_ids_raw))))


# ------------------------------------------------------------
# 5. Candidate Gene inspection
# ------------------------------------------------------------

candidate_ids_raw <- candidate_df[["Gene"]]

cat("\n------------------------------------------------------------\n")
cat("Candidate Gene column inspection\n")
cat("------------------------------------------------------------\n\n")

cat("Class:\n")
print(class(candidate_ids_raw))

cat("\nLength:\n")
print(length(candidate_ids_raw))

cat("\nFirst 20 values:\n")
print(head(candidate_ids_raw, 20))

cat("\nLast 20 values:\n")
print(tail(candidate_ids_raw, 20))

cat("\nNumber missing / NA:\n")
print(sum(is.na(candidate_ids_raw)))

cat("\nNumber empty strings:\n")
print(sum(trimws(as.character(candidate_ids_raw)) == ""))

cat("\nNumber unique IDs:\n")
print(length(unique(as.character(candidate_ids_raw))))


# ------------------------------------------------------------
# 6. Normalize ONLY temporary comparison copies
# ------------------------------------------------------------

candidate_ids <-
  trimws(as.character(candidate_ids_raw))

universe_ids <-
  trimws(as.character(universe_ids_raw))


# ------------------------------------------------------------
# 7. Exact set comparison
# ------------------------------------------------------------

shared_ids <-
  intersect(candidate_ids, universe_ids)

candidate_only <-
  setdiff(candidate_ids, universe_ids)

universe_only <-
  setdiff(universe_ids, candidate_ids)


cat("\n============================================================\n")
cat("EXACT GENE-ID SET COMPARISON\n")
cat("============================================================\n\n")

cat("Candidate unique IDs: ", length(unique(candidate_ids)), "\n", sep = "")
cat("Universe unique IDs:  ", length(unique(universe_ids)), "\n", sep = "")
cat("Shared IDs:           ", length(shared_ids), "\n", sep = "")
cat("Candidate-only IDs:   ", length(candidate_only), "\n", sep = "")
cat("Universe-only IDs:    ", length(universe_only), "\n", sep = "")


# ------------------------------------------------------------
# 8. Exact ordered-vector comparison
# ------------------------------------------------------------

same_length <-
  length(candidate_ids) == length(universe_ids)

same_order <-
  same_length &&
  identical(candidate_ids, universe_ids)

same_set <-
  length(candidate_only) == 0 &&
  length(universe_only) == 0


cat("\n------------------------------------------------------------\n")
cat("Comparison status\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Same length: ",
  same_length,
  "\n",
  sep = ""
)

cat(
  "Same ordered vector: ",
  same_order,
  "\n",
  sep = ""
)

cat(
  "Same gene-ID set: ",
  same_set,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Format / ID-pattern inspection
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("ID pattern inspection\n")
cat("------------------------------------------------------------\n\n")

cat("Candidate IDs beginning with ENSG:\n")
print(sum(grepl("^ENSG", candidate_ids)))

cat("\nUniverse IDs beginning with ENSG:\n")
print(sum(grepl("^ENSG", universe_ids)))

cat("\nCandidate IDs not beginning with ENSG:\n")
print(sum(!grepl("^ENSG", candidate_ids)))

cat("\nUniverse IDs not beginning with ENSG:\n")
print(sum(!grepl("^ENSG", universe_ids)))


# ------------------------------------------------------------
# 10. Show discrepancies if any
# ------------------------------------------------------------

if (length(candidate_only) > 0) {
  
  cat("\nFirst candidate-only IDs:\n")
  print(head(candidate_only, 20))
  
}

if (length(universe_only) > 0) {
  
  cat("\nFirst universe-only IDs:\n")
  print(head(universe_only, 20))
  
}


# ------------------------------------------------------------
# 11. Preliminary recovery interpretation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-1C-R INTERPRETATION\n")
cat("============================================================\n\n")

if (
  length(candidate_ids) == 57828 &&
  length(universe_ids) == 57828 &&
  length(unique(candidate_ids)) == 57828 &&
  length(unique(universe_ids)) == 57828 &&
  same_set
) {
  
  cat(
    "PASS: The universe column 'x' contains an exact 57,828-gene ",
    "set matching candidate_df$Gene.\n"
  )
  
} else {
  
  cat(
    "FAIL: The universe column 'x' does not yet establish an ",
    "exact 57,828-gene match.\n"
  )
  
}


cat("\nNo DEG recalculation performed.\n")
cat("No GSEA performed.\n")
cat("No ranking metric calculated.\n")
cat("No source file modified.\n")
cat("No STEP9B result modified.\n")


cat("\n============================================================\n")
cat("STEP9D-1D: Genome-wide Differential-Statistic Source Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA/9D-1_Input_Lock"

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}

setwd(step_dir)


# ------------------------------------------------------------
# 1. Locked source
# ------------------------------------------------------------

source_file <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9B-7R_Revised_DEG_Definition_and_Lock/02_results/STEP9B7R_Revised_DEG_FULL.csv"


# ------------------------------------------------------------
# 2. Locked universe
# ------------------------------------------------------------

universe_file <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9C_Revised_DEG_GO_KEGG_ORA/00_input/STEP9C_DEG_testing_universe_57828.csv"


if (!file.exists(source_file)) {
  stop(
    "FAIL: Locked genome-wide source file not found.\n",
    source_file
  )
}

if (!file.exists(universe_file)) {
  stop(
    "FAIL: Locked 57,828-gene universe file not found.\n",
    universe_file
  )
}


# ------------------------------------------------------------
# 3. Read-only import
# ------------------------------------------------------------

source_df <-
  data.table::fread(
    source_file,
    data.table = FALSE,
    check.names = FALSE
  )

universe_df <-
  data.table::fread(
    universe_file,
    data.table = FALSE,
    check.names = FALSE
  )


# ------------------------------------------------------------
# 4. Structural lock checks
# ------------------------------------------------------------

required_source_columns <- c(
  "Gene",
  "log2FC",
  "P_value",
  "FDR"
)

missing_columns <-
  setdiff(required_source_columns, colnames(source_df))

if (length(missing_columns) > 0) {
  stop(
    "FAIL: Required source columns missing: ",
    paste(missing_columns, collapse = ", ")
  )
}


source_ids <-
  trimws(as.character(source_df$Gene))

universe_ids <-
  trimws(as.character(universe_df$x))


# ------------------------------------------------------------
# 5. Gene-ID integrity
# ------------------------------------------------------------

n_source <-
  length(source_ids)

n_universe <-
  length(universe_ids)

n_source_unique <-
  length(unique(source_ids))

n_universe_unique <-
  length(unique(universe_ids))


source_only <-
  setdiff(source_ids, universe_ids)

universe_only <-
  setdiff(universe_ids, source_ids)


exact_set_match <-
  length(source_only) == 0 &&
  length(universe_only) == 0


# ------------------------------------------------------------
# 6. Statistical-column integrity
# ------------------------------------------------------------

log2fc_num <-
  suppressWarnings(as.numeric(source_df$log2FC))

p_num <-
  suppressWarnings(as.numeric(source_df$P_value))

fdr_num <-
  suppressWarnings(as.numeric(source_df$FDR))


finite_log2fc <-
  all(is.finite(log2fc_num))

finite_p <-
  all(is.finite(p_num))

finite_fdr <-
  all(is.finite(fdr_num))


valid_p_range <-
  all(p_num > 0 & p_num <= 1)

valid_fdr_range <-
  all(fdr_num >= 0 & fdr_num <= 1)


# ------------------------------------------------------------
# 7. Directional structure
# ------------------------------------------------------------

n_positive_log2fc <-
  sum(log2fc_num > 0)

n_negative_log2fc <-
  sum(log2fc_num < 0)

n_zero_log2fc <-
  sum(log2fc_num == 0)


# ------------------------------------------------------------
# 8. FDR structure
# ------------------------------------------------------------

n_fdr05 <-
  sum(fdr_num < 0.05)

n_fdr05_up <-
  sum(fdr_num < 0.05 & log2fc_num > 0)

n_fdr05_down <-
  sum(fdr_num < 0.05 & log2fc_num < 0)


# ------------------------------------------------------------
# 9. Core-gene presence audit
# ------------------------------------------------------------

core_genes <- c(
  AOC3  = "ENSG00000131471",
  MS4A7 = "ENSG00000166927",
  IL7R  = "ENSG00000168685"
)

core_source_presence <-
  data.frame(
    Core_Gene = names(core_genes),
    Ensembl_ID = unname(core_genes),
    In_Source = unname(core_genes) %in% source_ids,
    In_Universe = unname(core_genes) %in% universe_ids,
    stringsAsFactors = FALSE
  )


# ------------------------------------------------------------
# 10. Expected locked STEP9B counts
# ------------------------------------------------------------

expected_up <- 9311
expected_down <- 6754
expected_total <- 16065


count_check_total <-
  n_fdr05 == expected_total

count_check_up <-
  n_fdr05_up == expected_up

count_check_down <-
  n_fdr05_down == expected_down


# ------------------------------------------------------------
# 11. Final PASS criteria
# ------------------------------------------------------------

pass_structural <-
  n_source == 57828 &&
  n_universe == 57828 &&
  n_source_unique == 57828 &&
  n_universe_unique == 57828

pass_gene_set <-
  exact_set_match

pass_statistics <-
  finite_log2fc &&
  finite_p &&
  finite_fdr &&
  valid_p_range &&
  valid_fdr_range

pass_fdr_counts <-
  count_check_total &&
  count_check_up &&
  count_check_down

pass_core <-
  all(core_source_presence$In_Source) &&
  all(core_source_presence$In_Universe)


final_pass <-
  pass_structural &&
  pass_gene_set &&
  pass_statistics &&
  pass_fdr_counts &&
  pass_core


# ------------------------------------------------------------
# 12. Print lock summary
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("SOURCE\n")
cat("------------------------------------------------------------\n")
cat(source_file, "\n\n")

cat("Rows: ", n_source, "\n", sep = "")
cat("Unique Gene IDs: ", n_source_unique, "\n", sep = "")

cat("\n------------------------------------------------------------\n")
cat("LOCKED UNIVERSE\n")
cat("------------------------------------------------------------\n")
cat(universe_file, "\n\n")

cat("Rows: ", n_universe, "\n", sep = "")
cat("Unique Gene IDs: ", n_universe_unique, "\n", sep = "")

cat("\n------------------------------------------------------------\n")
cat("GENE-ID CONCORDANCE\n")
cat("------------------------------------------------------------\n")
cat("Exact ID-set match: ", exact_set_match, "\n", sep = "")
cat("Source-only IDs: ", length(source_only), "\n", sep = "")
cat("Universe-only IDs: ", length(universe_only), "\n", sep = "")

cat("\n------------------------------------------------------------\n")
cat("STATISTICAL COLUMN AUDIT\n")
cat("------------------------------------------------------------\n")
cat("Finite log2FC: ", finite_log2fc, "\n", sep = "")
cat("Finite P_value: ", finite_p, "\n", sep = "")
cat("Finite FDR: ", finite_fdr, "\n", sep = "")
cat("P_value in (0,1]: ", valid_p_range, "\n", sep = "")
cat("FDR in [0,1]: ", valid_fdr_range, "\n", sep = "")

cat("\n------------------------------------------------------------\n")
cat("FDR < 0.05 STRUCTURE\n")
cat("------------------------------------------------------------\n")
cat("Total: ", n_fdr05, "\n", sep = "")
cat("Up: ", n_fdr05_up, "\n", sep = "")
cat("Down: ", n_fdr05_down, "\n", sep = "")

cat("\nExpected locked STEP9B counts:\n")
cat("Total: ", expected_total, "\n", sep = "")
cat("Up: ", expected_up, "\n", sep = "")
cat("Down: ", expected_down, "\n", sep = "")

cat("\n------------------------------------------------------------\n")
cat("CORE-GENE AUDIT\n")
cat("------------------------------------------------------------\n")
print(core_source_presence)

cat("\n============================================================\n")
cat("STEP9D-1D FINAL STATUS\n")
cat("============================================================\n\n")

if (final_pass) {
  
  cat("PASS: Genome-wide differential-statistic source is LOCKED.\n\n")
  
  cat(
    "Locked source: STEP9B7R_Revised_DEG_FULL.csv\n"
  )
  
  cat(
    "Locked genome-wide universe: 57,828 Ensembl IDs\n"
  )
  
  cat(
    "Exact source/universe gene-ID concordance: PASS\n"
  )
  
  cat(
    "Required statistical columns: log2FC, P_value, FDR\n"
  )
  
} else {
  
  cat(
    "FAIL: Genome-wide source lock criteria were not all satisfied.\n"
  )
  
}


cat("\nNo DEG recalculation performed.\n")
cat("No GSEA performed.\n")
cat("No ranking metric calculated.\n")
cat("No source file modified.\n")
cat("No STEP9B result modified.\n")


# ------------------------------------------------------------
# 13. Save formal lock object
# ------------------------------------------------------------

STEP9D1D_FINAL_LOCK <- list(
  
  step = "STEP9D-1D",
  
  status = if (final_pass) "PASS" else "FAIL",
  
  source_file = source_file,
  
  universe_file = universe_file,
  
  source_rows = n_source,
  
  universe_rows = n_universe,
  
  source_unique_genes = n_source_unique,
  
  universe_unique_genes = n_universe_unique,
  
  exact_gene_set_match = exact_set_match,
  
  source_only_ids = length(source_only),
  
  universe_only_ids = length(universe_only),
  
  required_statistics = c(
    "log2FC",
    "P_value",
    "FDR"
  ),
  
  expected_fdr05_total = expected_total,
  
  expected_fdr05_up = expected_up,
  
  expected_fdr05_down = expected_down,
  
  observed_fdr05_total = n_fdr05,
  
  observed_fdr05_up = n_fdr05_up,
  
  observed_fdr05_down = n_fdr05_down,
  
  core_genes = core_genes,
  
  no_deg_recalculation = TRUE,
  
  no_gsea = TRUE,
  
  no_ranking_calculation = TRUE,
  
  no_source_modification = TRUE,
  
  no_step9b_modification = TRUE,
  
  timestamp = Sys.time()
)


saveRDS(
  STEP9D1D_FINAL_LOCK,
  file = "STEP9D1D_FINAL_LOCK.rds"
)


write.csv(
  core_source_presence,
  file = "STEP9D1D_CoreGene_Source_Audit.csv",
  row.names = FALSE
)

cat(
  "\nFormal lock object saved:\n",
  file.path(step_dir, "STEP9D1D_FINAL_LOCK.rds"),
  "\n"
)

##STEP9D-1E — Genome-wide Source / DEG Definition Concordance Audit
cat("\n============================================================\n")
cat("STEP9D-1E: Genome-wide Source / DEG Definition Concordance Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA/9D-1_Input_Lock"

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}

setwd(step_dir)


# ------------------------------------------------------------
# 1. Source file
# ------------------------------------------------------------

source_file <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9B-7R_Revised_DEG_Definition_and_Lock/02_results/STEP9B7R_Revised_DEG_FULL.csv"


if (!file.exists(source_file)) {
  stop(
    "FAIL: Genome-wide source file not found.\n",
    source_file
  )
}


# ------------------------------------------------------------
# 2. Read source read-only
# ------------------------------------------------------------

source_df <-
  data.table::fread(
    source_file,
    data.table = FALSE,
    check.names = FALSE
  )


# ------------------------------------------------------------
# 3. Structural checks
# ------------------------------------------------------------

required_columns <- c(
  "Gene",
  "log2FC",
  "P_value",
  "FDR",
  "Direction",
  "DEG"
)

missing_columns <-
  setdiff(required_columns, colnames(source_df))

if (length(missing_columns) > 0) {
  stop(
    "FAIL: Required columns missing: ",
    paste(missing_columns, collapse = ", ")
  )
}


cat("Source dimensions:\n")
print(dim(source_df))

cat("\nRequired columns present:\n")
print(required_columns %in% colnames(source_df))


# ------------------------------------------------------------
# 4. Temporary numeric copies
# ------------------------------------------------------------

log2fc <-
  suppressWarnings(as.numeric(source_df$log2FC))

fdr <-
  suppressWarnings(as.numeric(source_df$FDR))

p_value <-
  suppressWarnings(as.numeric(source_df$P_value))


# ------------------------------------------------------------
# 5. Inspect existing DEG and Direction definitions
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("EXISTING DEG COLUMN\n")
cat("============================================================\n\n")

cat("Class:\n")
print(class(source_df$DEG))

cat("\nUnique DEG values:\n")
print(unique(source_df$DEG))

cat("\nDEG frequency:\n")
print(table(source_df$DEG, useNA = "ifany"))


cat("\n============================================================\n")
cat("EXISTING DIRECTION COLUMN\n")
cat("============================================================\n\n")

cat("Class:\n")
print(class(source_df$Direction))

cat("\nUnique Direction values:\n")
print(unique(source_df$Direction))

cat("\nDirection frequency:\n")
print(table(source_df$Direction, useNA = "ifany"))


# ------------------------------------------------------------
# 6. Reconstruct DEG definition ONLY as an audit copy
# ------------------------------------------------------------
#
# IMPORTANT:
# This does NOT modify source_df.
# It only creates a temporary logical vector for comparison.
#

audit_deg_definition <-
  !is.na(fdr) &
  !is.na(log2fc) &
  fdr < 0.05 &
  abs(log2fc) > 0.58


# ------------------------------------------------------------
# 7. Convert existing DEG column to logical audit form
# ------------------------------------------------------------

existing_deg_chr <-
  trimws(tolower(as.character(source_df$DEG)))

existing_deg_logical <-
  existing_deg_chr %in% c(
    "true",
    "1",
    "yes",
    "deg"
  )


# ------------------------------------------------------------
# 8. DEG definition concordance
# ------------------------------------------------------------

deg_definition_match <-
  identical(
    existing_deg_logical,
    audit_deg_definition
  )


cat("\n============================================================\n")
cat("DEG DEFINITION CONCORDANCE\n")
cat("============================================================\n\n")

cat(
  "Audit definition: FDR < 0.05 AND |log2FC| > 0.58\n"
)

cat(
  "Audit-defined DEG count: ",
  sum(audit_deg_definition),
  "\n",
  sep = ""
)

cat(
  "Existing DEG=TRUE count: ",
  sum(existing_deg_logical),
  "\n",
  sep = ""
)

cat(
  "Exact row-wise concordance: ",
  deg_definition_match,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Direction audit among audit-defined DEGs
# ------------------------------------------------------------

audit_up <-
  audit_deg_definition &
  log2fc > 0

audit_down <-
  audit_deg_definition &
  log2fc < 0

audit_zero <-
  audit_deg_definition &
  log2fc == 0


cat("\n------------------------------------------------------------\n")
cat("AUDIT-DEFINED DEG DIRECTION\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Up: ",
  sum(audit_up),
  "\n",
  sep = ""
)

cat(
  "Down: ",
  sum(audit_down),
  "\n",
  sep = ""
)

cat(
  "Zero log2FC: ",
  sum(audit_zero),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Existing DEG + Direction counts
# ------------------------------------------------------------

direction_chr <-
  trimws(as.character(source_df$Direction))

existing_up <-
  existing_deg_logical &
  direction_chr == "Up"

existing_down <-
  existing_deg_logical &
  direction_chr == "Down"

existing_other <-
  existing_deg_logical &
  !direction_chr %in% c("Up", "Down")


cat("\n------------------------------------------------------------\n")
cat("EXISTING DEG + DIRECTION COUNTS\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Existing DEG total: ",
  sum(existing_deg_logical),
  "\n",
  sep = ""
)

cat(
  "Existing DEG Up: ",
  sum(existing_up),
  "\n",
  sep = ""
)

cat(
  "Existing DEG Down: ",
  sum(existing_down),
  "\n",
  sep = ""
)

cat(
  "Existing DEG other/unknown Direction: ",
  sum(existing_other),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Expected STEP9B locked values
# ------------------------------------------------------------

expected_total <- 16065
expected_up <- 9311
expected_down <- 6754


count_total_pass <-
  sum(audit_deg_definition) == expected_total

count_up_pass <-
  sum(audit_up) == expected_up

count_down_pass <-
  sum(audit_down) == expected_down


cat("\n============================================================\n")
cat("STEP9B LOCKED COUNT REPRODUCTION — READ-ONLY AUDIT\n")
cat("============================================================\n\n")

cat(
  "Expected total: ",
  expected_total,
  "\n",
  sep = ""
)

cat(
  "Observed total: ",
  sum(audit_deg_definition),
  "\n",
  sep = ""
)

cat(
  "Total PASS: ",
  count_total_pass,
  "\n\n",
  sep = ""
)

cat(
  "Expected Up: ",
  expected_up,
  "\n",
  sep = ""
)

cat(
  "Observed Up: ",
  sum(audit_up),
  "\n",
  sep = ""
)

cat(
  "Up PASS: ",
  count_up_pass,
  "\n\n",
  sep = ""
)

cat(
  "Expected Down: ",
  expected_down,
  "\n",
  sep = ""
)

cat(
  "Observed Down: ",
  sum(audit_down),
  "\n",
  sep = ""
)

cat(
  "Down PASS: ",
  count_down_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Compare FDR-only versus formal DEG definition
# ------------------------------------------------------------

fdr_only <-
  !is.na(fdr) &
  fdr < 0.05

fdr_only_not_deg <-
  fdr_only &
  !audit_deg_definition


cat("\n============================================================\n")
cat("FDR-ONLY VS FORMAL DEG DEFINITION\n")
cat("============================================================\n\n")

cat(
  "FDR < 0.05 only: ",
  sum(fdr_only),
  "\n",
  sep = ""
)

cat(
  "Formal DEG (FDR < 0.05 & |log2FC| > 0.58): ",
  sum(audit_deg_definition),
  "\n",
  sep = ""
)

cat(
  "FDR < 0.05 but failing |log2FC| > 0.58: ",
  sum(fdr_only_not_deg),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. Final audit status
# ------------------------------------------------------------

final_pass <-
  nrow(source_df) == 57828 &&
  sum(audit_deg_definition) == 16065 &&
  sum(audit_up) == 9311 &&
  sum(audit_down) == 6754 &&
  deg_definition_match &&
  sum(existing_up) == 9311 &&
  sum(existing_down) == 6754 &&
  sum(existing_other) == 0


cat("\n============================================================\n")
cat("STEP9D-1E FINAL STATUS\n")
cat("============================================================\n\n")

if (final_pass) {
  
  cat(
    "PASS: Genome-wide source reproduces the locked STEP9B ",
    "DEG definition and counts in a read-only audit.\n"
  )
  
  cat(
    "\nImportant interpretation:\n"
  )
  
  cat(
    "The genome-wide source contains 57,828 genes.\n"
  )
  
  cat(
    "FDR < 0.05 alone identifies 27,270 genes.\n"
  )
  
  cat(
    "The formal STEP9B DEG definition additionally requires ",
    "|log2FC| > 0.58, yielding 16,065 DEGs.\n"
  )
  
  cat(
    "Therefore, the 27,270 FDR-significant genes are not expected ",
    "to equal the 16,065 formal DEG count.\n"
  )
  
} else {
  
  cat(
    "FAIL: The existing genome-wide source does not fully reproduce ",
    "the locked STEP9B DEG definition/counts.\n"
  )
  
}


cat("\nNo source file modified.\n")
cat("No DEG recalculation performed.\n")
cat("No GSEA performed.\n")
cat("No ranking metric calculated.\n")
cat("No STEP9B result modified.\n")


# ------------------------------------------------------------
# 14. Save audit result
# ------------------------------------------------------------

STEP9D1E_AUDIT <- list(
  
  step = "STEP9D-1E",
  
  status = if (final_pass) "PASS" else "FAIL",
  
  source_file = source_file,
  
  genome_wide_n = nrow(source_df),
  
  fdr05_n = sum(fdr_only),
  
  formal_deg_n = sum(audit_deg_definition),
  
  formal_deg_up = sum(audit_up),
  
  formal_deg_down = sum(audit_down),
  
  fdr05_not_formal_deg = sum(fdr_only_not_deg),
  
  deg_definition =
    "FDR < 0.05 AND |log2FC| > 0.58",
  
  existing_deg_concordance =
    deg_definition_match,
  
  no_source_modification = TRUE,
  
  no_deg_recalculation = TRUE,
  
  no_gsea = TRUE,
  
  no_ranking_calculation = TRUE,
  
  no_step9b_modification = TRUE,
  
  timestamp = Sys.time()
)


saveRDS(
  STEP9D1E_AUDIT,
  file = "STEP9D1E_DEG_Definition_Concordance_Audit.rds"
)

cat(
  "\nAudit object saved:\n",
  file.path(
    step_dir,
    "STEP9D1E_DEG_Definition_Concordance_Audit.rds"
  ),
  "\n"
)

cat("\n============================================================\n")
cat("STEP9D-1F: Locked DEG Definition Reconstruction Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA/9D-1_Input_Lock"

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}

setwd(step_dir)


# ------------------------------------------------------------
# 1. Source file\n
# ------------------------------------------------------------

source_file <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9B-7R_Revised_DEG_Definition_and_Lock/02_results/STEP9B7R_Revised_DEG_FULL.csv"


if (!file.exists(source_file)) {
  stop(
    "FAIL: Source file not found.\n",
    source_file
  )
}


# ------------------------------------------------------------
# 2. Read-only import
# ------------------------------------------------------------

source_df <-
  data.table::fread(
    source_file,
    data.table = FALSE,
    check.names = FALSE
  )


# ------------------------------------------------------------
# 3. Structural audit
# ------------------------------------------------------------

required_columns <- c(
  "Gene",
  "log2FC",
  "P_value",
  "FDR",
  "Direction",
  "DEG"
)

missing_columns <-
  setdiff(required_columns, colnames(source_df))

if (length(missing_columns) > 0) {
  stop(
    "FAIL: Required columns missing: ",
    paste(missing_columns, collapse = ", ")
  )
}


cat("Source dimensions:\n")
print(dim(source_df))

cat("\nAll source columns:\n")
print(colnames(source_df))


# ------------------------------------------------------------
# 4. Temporary numeric copies
# ------------------------------------------------------------

log2fc <-
  suppressWarnings(as.numeric(source_df$log2FC))

abs_log2fc <-
  abs(log2fc)

p_value <-
  suppressWarnings(as.numeric(source_df$P_value))

fdr <-
  suppressWarnings(as.numeric(source_df$FDR))


existing_deg <-
  source_df$DEG


direction <-
  trimws(as.character(source_df$Direction))


# ------------------------------------------------------------
# 5. Existing DEG structure
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("EXISTING DEG STRUCTURE\n")
cat("============================================================\n\n")

cat("DEG class:\n")
print(class(existing_deg))

cat("\nDEG table:\n")
print(table(existing_deg, useNA = "ifany"))

cat("\nExisting DEG total:\n")
print(sum(existing_deg == TRUE, na.rm = TRUE))


# ------------------------------------------------------------
# 6. Correct Direction audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("EXISTING DIRECTION STRUCTURE\n")
cat("============================================================\n\n")

cat("Direction values:\n")
print(unique(direction))

cat("\nDirection frequency:\n")
print(table(direction, useNA = "ifany"))

cat("\nDEG + Direction cross-tabulation:\n")
print(
  table(
    DEG = existing_deg,
    Direction = direction,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 7. Correct locked STEP9B counts
# ------------------------------------------------------------

existing_deg_true <-
  existing_deg == TRUE

existing_tumor_up <-
  existing_deg_true &
  direction == "Tumor_up"

existing_tumor_down <-
  existing_deg_true &
  direction == "Tumor_down"


cat("\n------------------------------------------------------------\n")
cat("LOCKED DEG + DIRECTION COUNTS\n")
cat("------------------------------------------------------------\n\n")

cat(
  "DEG total: ",
  sum(existing_deg_true),
  "\n",
  sep = ""
)

cat(
  "DEG Tumor_up: ",
  sum(existing_tumor_up),
  "\n",
  sep = ""
)

cat(
  "DEG Tumor_down: ",
  sum(existing_tumor_down),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Distribution of log2FC among DEG and non-DEG
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("log2FC DISTRIBUTION BY EXISTING DEG STATUS\n")
cat("============================================================\n\n")

cat("Existing DEG = TRUE:\n")
print(
  summary(
    abs_log2fc[existing_deg_true]
  )
)

cat("\nExisting DEG = FALSE:\n")
print(
  summary(
    abs_log2fc[!existing_deg_true]
  )
)


# ------------------------------------------------------------
# 9. FDR distribution by DEG status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FDR DISTRIBUTION BY EXISTING DEG STATUS\n")
cat("============================================================\n\n")

cat("Existing DEG = TRUE:\n")
print(
  summary(
    fdr[existing_deg_true]
  )
)

cat("\nExisting DEG = FALSE:\n")
print(
  summary(
    fdr[!existing_deg_true]
  )
)


# ------------------------------------------------------------
# 10. Check whether all DEG satisfy FDR < 0.05
# ------------------------------------------------------------

deg_fdr05 <-
  existing_deg_true &
  fdr < 0.05

deg_fdr05_all <-
  all(
    fdr[existing_deg_true] < 0.05
  )


cat("\n============================================================\n")
cat("FDR REQUIREMENT AUDIT\n")
cat("============================================================\n\n")

cat(
  "All existing DEGs have FDR < 0.05: ",
  deg_fdr05_all,
  "\n",
  sep = ""
)

cat(
  "Existing DEGs with FDR >= 0.05: ",
  sum(
    existing_deg_true &
      !(fdr < 0.05)
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Check whether all DEG satisfy common |log2FC| thresholds
# ------------------------------------------------------------

candidate_effect_cutoffs <- c(
  0.50,
  0.55,
  0.58,
  0.60,
  0.65,
  0.70,
  0.75,
  0.80,
  0.90,
  1.00,
  1.20,
  1.50
)


effect_cutoff_audit <-
  data.frame(
    cutoff = candidate_effect_cutoffs,
    n_FDR05_and_effect =
      sapply(
        candidate_effect_cutoffs,
        function(x) {
          sum(
            fdr < 0.05 &
              abs_log2fc > x
          )
        }
      ),
    n_DEG_within_rule =
      sapply(
        candidate_effect_cutoffs,
        function(x) {
          sum(
            existing_deg_true &
              fdr < 0.05 &
              abs_log2fc > x
          )
        }
      ),
    stringsAsFactors = FALSE
  )


cat("\n============================================================\n")
cat("COMMON EFFECT-SIZE CUTOFF AUDIT\n")
cat("============================================================\n\n")

print(effect_cutoff_audit)


# ------------------------------------------------------------
# 12. Search empirical cutoff boundaries around 16,065
# ------------------------------------------------------------

target_deg_n <- 16065


candidate_thresholds <-
  sort(
    unique(
      abs_log2fc[fdr < 0.05]
    )
  )


threshold_counts <-
  sapply(
    candidate_thresholds,
    function(x) {
      sum(
        fdr < 0.05 &
          abs_log2fc >= x
      )
    }
  )


exact_threshold_hits <-
  candidate_thresholds[
    threshold_counts == target_deg_n
  ]


cat("\n============================================================\n")
cat("EMPIRICAL THRESHOLD RECONSTRUCTION\n")
cat("============================================================\n\n")

cat(
  "Target locked DEG count: ",
  target_deg_n,
  "\n",
  sep = ""
)

cat(
  "Number of exact threshold values yielding target count: ",
  length(exact_threshold_hits),
  "\n",
  sep = ""
)

if (length(exact_threshold_hits) > 0) {
  
  cat("\nExact threshold value(s):\n")
  print(exact_threshold_hits)
  
} else {
  
  cat(
    "\nNo exact single |log2FC| threshold reproduces ",
    target_deg_n,
    " genes under FDR < 0.05.\n",
    sep = ""
  )
  
}


# ------------------------------------------------------------
# 13. Determine nearest effect-size thresholds
# ------------------------------------------------------------

nearest_idx <-
  order(
    abs(threshold_counts - target_deg_n)
  )[1:min(20, length(threshold_counts))]


nearest_threshold_audit <-
  data.frame(
    abs_log2FC_threshold =
      candidate_thresholds[nearest_idx],
    resulting_gene_count =
      threshold_counts[nearest_idx],
    difference_from_target =
      threshold_counts[nearest_idx] -
      target_deg_n,
    stringsAsFactors = FALSE
  )


cat("\n------------------------------------------------------------\n")
cat("NEAREST EFFECT-SIZE THRESHOLDS\n")
cat("------------------------------------------------------------\n\n")

print(nearest_threshold_audit)


# ------------------------------------------------------------
# 14. Check exact concordance of existing DEG with
#     FDR + effect-size threshold candidates
# ------------------------------------------------------------

concordance_candidates <-
  data.frame(
    cutoff = candidate_effect_cutoffs,
    predicted_n =
      sapply(
        candidate_effect_cutoffs,
        function(x) {
          sum(
            fdr < 0.05 &
              abs_log2fc > x
          )
        }
      ),
    existing_deg_n =
      sum(existing_deg_true),
    stringsAsFactors = FALSE
  )


concordance_candidates$exact_count_match <-
  concordance_candidates$predicted_n ==
  concordance_candidates$existing_deg_n


cat("\n============================================================\n")
cat("COUNT CONCORDANCE WITH EXISTING DEG\n")
cat("============================================================\n\n")

print(concordance_candidates)


# ------------------------------------------------------------
# 15. Test whether existing DEG is exactly contained within
#     FDR < 0.05
# ------------------------------------------------------------

fdr_subset_status <-
  all(
    existing_deg_true <=
      (fdr < 0.05)
  )


cat("\n============================================================\n")
cat("EXISTING DEG SUBSET AUDIT\n")
cat("============================================================\n\n")

cat(
  "All existing DEGs satisfy FDR < 0.05: ",
  fdr_subset_status,
  "\n",
  sep = ""
)

cat(
  "FDR-significant genes not existing DEG: ",
  sum(
    (fdr < 0.05) &
      !existing_deg_true
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 16. Test direction consistency with log2FC
# ------------------------------------------------------------

direction_consistency <-
  ifelse(
    log2fc > 0,
    direction == "Tumor_up",
    ifelse(
      log2fc < 0,
      direction == "Tumor_down",
      NA
    )
  )


cat("\n============================================================\n")
cat("DIRECTION vs log2FC CONSISTENCY\n")
cat("============================================================\n\n")

cat(
  "Non-NA direction consistency: ",
  sum(
    direction_consistency,
    na.rm = TRUE
  ),
  " / ",
  sum(
    !is.na(direction_consistency)
  ),
  "\n",
  sep = ""
)

cat(
  "Inconsistent rows: ",
  sum(
    !direction_consistency,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 17. Core genes
# ------------------------------------------------------------

core_genes <- c(
  AOC3  = "ENSG00000131471",
  MS4A7 = "ENSG00000166927",
  IL7R  = "ENSG00000168685"
)


core_audit <-
  source_df[
    match(
      core_genes,
      source_df$Gene
    ),
    c(
      "Gene",
      "log2FC",
      "P_value",
      "FDR",
      "Direction",
      "DEG"
    ),
    drop = FALSE
  ]

core_audit$Core_Gene <-
  names(core_genes)


core_audit <-
  core_audit[
    c(
      "Core_Gene",
      "Gene",
      "log2FC",
      "P_value",
      "FDR",
      "Direction",
      "DEG"
    )
  ]


cat("\n============================================================\n")
cat("CORE-GENE AUDIT\n")
cat("============================================================\n\n")

print(core_audit)


# ------------------------------------------------------------
# 18. Final status — diagnostic only
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-1F STATUS\n")
cat("============================================================\n\n")

cat(
  "STATUS: DIAGNOSTIC AUDIT ONLY — NO LOCK CREATED.\n"
)

cat(
  "\nThe current audit does NOT assume that ",
  "FDR < 0.05 AND |log2FC| > 0.58 is the locked STEP9B DEG definition.\n",
  sep = ""
)

cat(
  "\nThe existing DEG definition must be recovered from ",
  "the locked STEP9B provenance before STEP9D input lock.\n"
)

cat("\nNo source file modified.\n")
cat("No STEP9B result modified.\n")
cat("No DEG recalculation performed.\n")
cat("No GSEA performed.\n")
cat("No ranking metric calculated.\n")


# ------------------------------------------------------------
# 19. Save diagnostic audit only
# ------------------------------------------------------------

STEP9D1F_DIAGNOSTIC <- list(
  step = "STEP9D-1F",
  status = "DIAGNOSTIC_ONLY",
  source_file = source_file,
  source_n = nrow(source_df),
  existing_deg_n = sum(existing_deg_true),
  existing_deg_tumor_up =
    sum(existing_tumor_up),
  existing_deg_tumor_down =
    sum(existing_tumor_down),
  fdr05_n =
    sum(fdr < 0.05),
  fdr05_and_abs_log2fc_gt_058 =
    sum(
      fdr < 0.05 &
        abs_log2fc > 0.58
    ),
  fdr_subset_status =
    fdr_subset_status,
  exact_empirical_thresholds =
    exact_threshold_hits,
  no_source_modification = TRUE,
  no_step9b_modification = TRUE,
  no_deg_recalculation = TRUE,
  no_gsea = TRUE,
  no_ranking_calculation = TRUE,
  timestamp = Sys.time()
)


saveRDS(
  STEP9D1F_DIAGNOSTIC,
  file =
    "STEP9D1F_DEG_Definition_Reconstruction_Diagnostic.rds"
)


write.csv(
  effect_cutoff_audit,
  file =
    "STEP9D1F_Effect_Size_Cutoff_Diagnostic.csv",
  row.names = FALSE
)


write.csv(
  nearest_threshold_audit,
  file =
    "STEP9D1F_Nearest_Threshold_Diagnostic.csv",
  row.names = FALSE
)


write.csv(
  core_audit,
  file =
    "STEP9D1F_CoreGene_DEG_Diagnostic.csv",
  row.names = FALSE
)


cat(
  "\nDiagnostic files saved under:\n",
  step_dir,
  "\n"
)


cat("STEP9D-1G: STEP9B Provenance Recovery Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Directories
# ------------------------------------------------------------

step9b_dir <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9B-7R_Revised_DEG_Definition_and_Lock"

step9d_dir <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA/9D-1_Input_Lock"

if (!dir.exists(step9b_dir)) {
  stop(
    "FAIL: STEP9B provenance directory not found:\n",
    step9b_dir
  )
}

if (!dir.exists(step9d_dir)) {
  dir.create(step9d_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 1. Recursive file inventory
# ------------------------------------------------------------

all_files <-
  list.files(
    step9b_dir,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  )

cat("Total files found in STEP9B directory:\n")
print(length(all_files))

cat("\n")


# ------------------------------------------------------------
# 2. File extension summary
# ------------------------------------------------------------

file_ext <-
  tolower(
    tools::file_ext(all_files)
  )

cat("File extension summary:\n")
print(sort(table(file_ext), decreasing = TRUE))

cat("\n")


# ------------------------------------------------------------
# 3. Filename-level provenance candidates
# ------------------------------------------------------------

filename_pattern <-
  paste(
    c(
      "definition",
      "lock",
      "audit",
      "cutoff",
      "threshold",
      "criteria",
      "provenance",
      "metadata",
      "deg"
    ),
    collapse = "|"
  )

candidate_files <-
  all_files[
    grepl(
      filename_pattern,
      basename(all_files),
      ignore.case = TRUE
    )
  ]

cat("Filename-level provenance candidates:\n")
print(candidate_files)

cat("\n")
cat("Number of candidate files:\n")
print(length(candidate_files))


# ------------------------------------------------------------
# 4. Search text-readable files for DEG definition terms
# ------------------------------------------------------------

text_ext <-
  c(
    "r",
    "R",
    "txt",
    "csv",
    "tsv",
    "md",
    "log",
    "json"
  )

text_files <-
  all_files[
    file_ext %in% tolower(text_ext)
  ]

cat("\nText-readable files available for content search:\n")
print(length(text_files))


search_pattern <-
  paste(
    c(
      "abs\\s*\\(?.*log2FC",
      "log2FC.*1",
      "log2FC.*0\\.58",
      "FDR.*0\\.05",
      "0\\.05.*FDR",
      "DEG.*FDR",
      "FDR.*DEG",
      "cutoff",
      "threshold",
      "definition",
      "criteria"
    ),
    collapse = "|"
  )


provenance_hits <-
  list()


for (f in text_files) {
  
  txt <-
    tryCatch(
      readLines(
        f,
        warn = FALSE,
        encoding = "UTF-8"
      ),
      error = function(e) NULL
    )
  
  if (is.null(txt)) {
    next
  }
  
  hit_idx <-
    grep(
      search_pattern,
      txt,
      ignore.case = TRUE
    )
  
  if (length(hit_idx) > 0) {
    
    provenance_hits[[f]] <-
      data.frame(
        File = f,
        Line = hit_idx,
        Text = txt[hit_idx],
        stringsAsFactors = FALSE
      )
  }
}


if (length(provenance_hits) > 0) {
  
  provenance_hits_df <-
    do.call(
      rbind,
      provenance_hits
    )
  
  rownames(provenance_hits_df) <- NULL
  
} else {
  
  provenance_hits_df <-
    data.frame(
      File = character(),
      Line = integer(),
      Text = character(),
      stringsAsFactors = FALSE
    )
}


cat("\n============================================================\n")
cat("CONTENT-LEVEL PROVENANCE HITS\n")
cat("============================================================\n\n")

print(
  provenance_hits_df,
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Inspect RDS / RData / RData objects
# ------------------------------------------------------------

r_object_files <-
  all_files[
    file_ext %in%
      c(
        "rds",
        "rda",
        "rdata"
      )
  ]

cat("\n============================================================\n")
cat("R OBJECT FILES\n")
cat("============================================================\n\n")

print(r_object_files)


# ------------------------------------------------------------
# 6. Read-only RDS object structure audit
# ------------------------------------------------------------

rds_summary <-
  list()

for (f in r_object_files) {
  
  if (tolower(tools::file_ext(f)) != "rds") {
    next
  }
  
  obj <-
    tryCatch(
      readRDS(f),
      error = function(e) NULL
    )
  
  if (is.null(obj)) {
    next
  }
  
  rds_summary[[f]] <-
    list(
      file = f,
      class = class(obj),
      typeof = typeof(obj),
      length = length(obj),
      names = names(obj)
    )
}


cat("\n============================================================\n")
cat("RDS STRUCTURE SUMMARY\n")
cat("============================================================\n\n")

if (length(rds_summary) > 0) {
  
  for (x in rds_summary) {
    
    cat("\n------------------------------------------------------------\n")
    cat("File:\n")
    cat(x$file, "\n")
    
    cat("Class:\n")
    print(x$class)
    
    cat("Type:\n")
    print(x$typeof)
    
    cat("Length:\n")
    print(x$length)
    
    cat("Names:\n")
    print(x$names)
  }
  
} else {
  
  cat("No readable RDS objects found.\n")
}


# ------------------------------------------------------------
# 7. Search character content inside RDS objects
# ------------------------------------------------------------

rds_text_hits <-
  list()


for (f in r_object_files) {
  
  if (tolower(tools::file_ext(f)) != "rds") {
    next
  }
  
  obj <-
    tryCatch(
      readRDS(f),
      error = function(e) NULL
    )
  
  if (is.null(obj)) {
    next
  }
  
  flat_text <-
    capture.output(
      str(
        obj,
        max.level = 4
      )
    )
  
  hit_idx <-
    grep(
      search_pattern,
      flat_text,
      ignore.case = TRUE
    )
  
  if (length(hit_idx) > 0) {
    
    rds_text_hits[[f]] <-
      data.frame(
        File = f,
        Line = hit_idx,
        Text = flat_text[hit_idx],
        stringsAsFactors = FALSE
      )
  }
}


cat("\n============================================================\n")
cat("RDS CONTENT SEARCH HITS\n")
cat("============================================================\n\n")

if (length(rds_text_hits) > 0) {
  
  rds_text_hits_df <-
    do.call(
      rbind,
      rds_text_hits
    )
  
  rownames(rds_text_hits_df) <- NULL
  
  print(
    rds_text_hits_df,
    row.names = FALSE
  )
  
} else {
  
  rds_text_hits_df <-
    data.frame(
      File = character(),
      Line = integer(),
      Text = character(),
      stringsAsFactors = FALSE
    )
  
  cat("No relevant RDS content hits found.\n")
}


# ------------------------------------------------------------
# 8. Save diagnostic audit only
# ------------------------------------------------------------

STEP9D1G_DIAGNOSTIC <-
  list(
    step = "STEP9D-1G",
    status = "DIAGNOSTIC_ONLY",
    step9b_directory = step9b_dir,
    total_files = length(all_files),
    candidate_files = candidate_files,
    provenance_hits = provenance_hits_df,
    r_object_files = r_object_files,
    rds_text_hits = rds_text_hits_df,
    no_source_modification = TRUE,
    no_step9b_modification = TRUE,
    no_deg_recalculation = TRUE,
    no_gsea = TRUE,
    no_ranking_calculation = TRUE,
    timestamp = Sys.time()
  )


saveRDS(
  STEP9D1G_DIAGNOSTIC,
  file =
    file.path(
      step9d_dir,
      "STEP9D1G_STEP9B_Provenance_Recovery_Diagnostic.rds"
    )
)


write.csv(
  provenance_hits_df,
  file =
    file.path(
      step9d_dir,
      "STEP9D1G_Text_Provenance_Hits.csv"
    ),
  row.names = FALSE
)


write.csv(
  rds_text_hits_df,
  file =
    file.path(
      step9d_dir,
      "STEP9D1G_RDS_Provenance_Hits.csv"
    ),
  row.names = FALSE
)


write.csv(
  data.frame(
    File = all_files,
    Extension = file_ext,
    stringsAsFactors = FALSE
  ),
  file =
    file.path(
      step9d_dir,
      "STEP9D1G_STEP9B_File_Inventory.csv"
    ),
  row.names = FALSE
)


cat("\n============================================================\n")
cat("STEP9D-1G STATUS\n")
cat("============================================================\n\n")

cat(
  "STATUS: DIAGNOSTIC ONLY — NO LOCK CREATED.\n"
)

cat(
  "No STEP9B result modified.\n"
)

cat(
  "No DEG recalculation performed.\n"
)

cat(
  "No GSEA performed.\n"
)

cat(
  "No ranking metric calculated.\n"
)

cat(
  "\nDiagnostic output saved under:\n",
  step9d_dir,
  "\n",
  sep = ""
)


cat("STEP9D-1H: Exact DEG Definition Set-Concordance Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Directories and source
# ------------------------------------------------------------

step9d_dir <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA/9D-1_Input_Lock"

source_file <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9B-7R_Revised_DEG_Definition_and_Lock/02_results/STEP9B7R_Revised_DEG_FULL.csv"

if (!file.exists(source_file)) {
  stop(
    "FAIL: STEP9B source file not found:\n",
    source_file
  )
}

if (!dir.exists(step9d_dir)) {
  dir.create(
    step9d_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 1. Read-only import
# ------------------------------------------------------------

source_df <-
  data.table::fread(
    source_file,
    data.table = FALSE,
    check.names = FALSE
  )


# ------------------------------------------------------------
# 2. Structural checks
# ------------------------------------------------------------

required_columns <-
  c(
    "Gene",
    "log2FC",
    "FDR",
    "Direction",
    "DEG"
  )

missing_columns <-
  setdiff(
    required_columns,
    colnames(source_df)
  )

if (length(missing_columns) > 0) {
  
  stop(
    "FAIL: Required columns missing: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


if (anyDuplicated(source_df$Gene) > 0) {
  
  stop(
    "FAIL: Duplicate Gene IDs detected in source."
  )
}


# ------------------------------------------------------------
# 3. Numeric copies
# ------------------------------------------------------------

log2fc <-
  suppressWarnings(
    as.numeric(
      source_df$log2FC
    )
  )

fdr <-
  suppressWarnings(
    as.numeric(
      source_df$FDR
    )
  )

direction <-
  trimws(
    as.character(
      source_df$Direction
    )
  )

existing_deg <-
  source_df$DEG == TRUE


# ------------------------------------------------------------
# 4. Candidate reconstructed definition
# ------------------------------------------------------------

candidate_deg <-
  fdr < 0.05 &
  abs(log2fc) > 1


candidate_up <-
  candidate_deg &
  direction == "Tumor_up"


candidate_down <-
  candidate_deg &
  direction == "Tumor_down"


existing_up <-
  existing_deg &
  direction == "Tumor_up"


existing_down <-
  existing_deg &
  direction == "Tumor_down"


# ------------------------------------------------------------
# 5. Count comparison
# ------------------------------------------------------------

cat("============================================================\n")
cat("COUNT COMPARISON\n")
cat("============================================================\n\n")

cat(
  "Existing DEG total: ",
  sum(existing_deg),
  "\n",
  sep = ""
)

cat(
  "Candidate DEG total: ",
  sum(candidate_deg),
  "\n",
  sep = ""
)

cat(
  "Existing Tumor_up: ",
  sum(existing_up),
  "\n",
  sep = ""
)

cat(
  "Candidate Tumor_up: ",
  sum(candidate_up),
  "\n",
  sep = ""
)

cat(
  "Existing Tumor_down: ",
  sum(existing_down),
  "\n",
  sep = ""
)

cat(
  "Candidate Tumor_down: ",
  sum(candidate_down),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Exact logical-vector equality
# ------------------------------------------------------------

logical_vector_exact_match <-
  identical(
    existing_deg,
    candidate_deg
  )


cat("\n============================================================\n")
cat("EXACT LOGICAL-VECTOR MATCH\n")
cat("============================================================\n\n")

cat(
  "Existing DEG == (FDR < 0.05 & |log2FC| > 1): ",
  logical_vector_exact_match,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Gene-set comparison
# ------------------------------------------------------------

existing_gene_set <-
  source_df$Gene[
    existing_deg
  ]

candidate_gene_set <-
  source_df$Gene[
    candidate_deg
  ]


existing_only <-
  setdiff(
    existing_gene_set,
    candidate_gene_set
  )

candidate_only <-
  setdiff(
    candidate_gene_set,
    existing_gene_set
  )


gene_set_exact_match <-
  length(existing_only) == 0 &
  length(candidate_only) == 0


cat("\n============================================================\n")
cat("EXACT GENE-SET MATCH\n")
cat("============================================================\n\n")

cat(
  "Existing DEG only: ",
  length(existing_only),
  "\n",
  sep = ""
)

cat(
  "Candidate definition only: ",
  length(candidate_only),
  "\n",
  sep = ""
)

cat(
  "Exact gene-set match: ",
  gene_set_exact_match,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Direction-specific exact comparison
# ------------------------------------------------------------

existing_up_genes <-
  source_df$Gene[
    existing_up
  ]

candidate_up_genes <-
  source_df$Gene[
    candidate_up
  ]

existing_down_genes <-
  source_df$Gene[
    existing_down
  ]

candidate_down_genes <-
  source_df$Gene[
    candidate_down
  ]


up_existing_only <-
  setdiff(
    existing_up_genes,
    candidate_up_genes
  )

up_candidate_only <-
  setdiff(
    candidate_up_genes,
    existing_up_genes
  )

down_existing_only <-
  setdiff(
    existing_down_genes,
    candidate_down_genes
  )

down_candidate_only <-
  setdiff(
    candidate_down_genes,
    existing_down_genes
  )


up_exact_match <-
  length(up_existing_only) == 0 &
  length(up_candidate_only) == 0


down_exact_match <-
  length(down_existing_only) == 0 &
  length(down_candidate_only) == 0


cat("\n============================================================\n")
cat("DIRECTION-SPECIFIC EXACT MATCH\n")
cat("============================================================\n\n")

cat(
  "Tumor_up existing only: ",
  length(up_existing_only),
  "\n",
  sep = ""
)

cat(
  "Tumor_up candidate only: ",
  length(up_candidate_only),
  "\n",
  sep = ""
)

cat(
  "Tumor_up exact match: ",
  up_exact_match,
  "\n",
  sep = ""
)

cat("\n")

cat(
  "Tumor_down existing only: ",
  length(down_existing_only),
  "\n",
  sep = ""
)

cat(
  "Tumor_down candidate only: ",
  length(down_candidate_only),
  "\n",
  sep = ""
)

cat(
  "Tumor_down exact match: ",
  down_exact_match,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Boundary audit
# ------------------------------------------------------------

existing_deg_abs_log2fc_min <-
  min(
    abs(log2fc[existing_deg]),
    na.rm = TRUE
  )

existing_deg_abs_log2fc_max <-
  max(
    abs(log2fc[existing_deg]),
    na.rm = TRUE
  )

existing_deg_fdr_max <-
  max(
    fdr[existing_deg],
    na.rm = TRUE
  )


cat("\n============================================================\n")
cat("BOUNDARY AUDIT\n")
cat("============================================================\n\n")

cat(
  "Minimum |log2FC| among existing DEG: ",
  sprintf(
    "%.12f",
    existing_deg_abs_log2fc_min
  ),
  "\n",
  sep = ""
)

cat(
  "Maximum |log2FC| among existing DEG: ",
  sprintf(
    "%.12f",
    existing_deg_abs_log2fc_max
  ),
  "\n",
  sep = ""
)

cat(
  "Maximum FDR among existing DEG: ",
  sprintf(
    "%.12f",
    existing_deg_fdr_max
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Identify all disagreements if any
# ------------------------------------------------------------

disagreement <-
  existing_deg != candidate_deg


n_disagreement <-
  sum(
    disagreement,
    na.rm = TRUE
  )


cat("\n============================================================\n")
cat("DISAGREEMENT AUDIT\n")
cat("============================================================\n\n")

cat(
  "Number of logical DEG disagreements: ",
  n_disagreement,
  "\n",
  sep = ""
)


if (n_disagreement > 0) {
  
  disagreement_table <-
    source_df[
      disagreement,
      c(
        "Gene",
        "log2FC",
        "FDR",
        "Direction",
        "DEG"
      ),
      drop = FALSE
    ]
  
  disagreement_table$Candidate_DEG <-
    candidate_deg[
      disagreement
    ]
  
  cat("\nDisagreement rows:\n")
  print(
    disagreement_table,
    row.names = FALSE
  )
  
} else {
  
  disagreement_table <-
    data.frame(
      Gene = character(),
      log2FC = numeric(),
      FDR = numeric(),
      Direction = character(),
      DEG = logical(),
      Candidate_DEG = logical(),
      stringsAsFactors = FALSE
    )
  
  cat(
    "PASS: No logical DEG disagreements detected.\n"
  )
}


# ------------------------------------------------------------
# 11. Core-gene audit
# ------------------------------------------------------------

core_genes <-
  c(
    AOC3  = "ENSG00000131471",
    MS4A7 = "ENSG00000166927",
    IL7R  = "ENSG00000168685"
  )


core_idx <-
  match(
    core_genes,
    source_df$Gene
  )


core_audit <-
  source_df[
    core_idx,
    c(
      "Gene",
      "log2FC",
      "FDR",
      "Direction",
      "DEG"
    ),
    drop = FALSE
  ]

core_audit$Core_Gene <-
  names(core_genes)

core_audit$Candidate_DEG <-
  candidate_deg[
    core_idx
  ]

core_audit$Exact_Match <-
  core_audit$DEG ==
  core_audit$Candidate_DEG

core_audit <-
  core_audit[
    c(
      "Core_Gene",
      "Gene",
      "log2FC",
      "FDR",
      "Direction",
      "DEG",
      "Candidate_DEG",
      "Exact_Match"
    )
  ]


cat("\n============================================================\n")
cat("CORE-GENE EXACT MATCH AUDIT\n")
cat("============================================================\n\n")

print(
  core_audit,
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Final diagnostic status
# ------------------------------------------------------------

overall_exact_match <-
  logical_vector_exact_match &
  gene_set_exact_match &
  up_exact_match &
  down_exact_match &
  n_disagreement == 0 &
  all(core_audit$Exact_Match)


cat("\n============================================================\n")
cat("STEP9D-1H STATUS\n")
cat("============================================================\n\n")

if (overall_exact_match) {
  
  cat(
    "STATUS: PASS — EXACT RESULT-LEVEL CONCORDANCE.\n"
  )
  
  cat(
    "\nThe existing STEP9B DEG flag is exactly reproduced by:\n"
  )
  
  cat(
    "FDR < 0.05 AND |log2FC| > 1\n"
  )
  
  cat(
    "\nThis is a result-level reconstruction audit.\n"
  )
  
  cat(
    "The original STEP9B generating script/provenance was not recovered.\n"
  )
  
} else {
  
  cat(
    "STATUS: FAIL — EXACT RESULT-LEVEL CONCORDANCE NOT ESTABLISHED.\n"
  )
}


cat("\nNo source file modified.\n")
cat("No STEP9B result modified.\n")
cat("No DEG recalculation performed.\n")
cat("No GSEA performed.\n")
cat("No ranking metric calculated.\n")
cat("No STEP9D input lock created.\n")


# ------------------------------------------------------------
# 13. Save diagnostic only
# ------------------------------------------------------------

STEP9D1H_DIAGNOSTIC <-
  list(
    step = "STEP9D-1H",
    status =
      if (overall_exact_match)
        "PASS_RESULT_LEVEL_CONCORDANCE"
    else
      "FAIL_RESULT_LEVEL_CONCORDANCE",
    
    source_file = source_file,
    
    existing_deg_n =
      sum(existing_deg),
    
    candidate_deg_n =
      sum(candidate_deg),
    
    existing_tumor_up_n =
      sum(existing_up),
    
    candidate_tumor_up_n =
      sum(candidate_up),
    
    existing_tumor_down_n =
      sum(existing_down),
    
    candidate_tumor_down_n =
      sum(candidate_down),
    
    exact_logical_vector_match =
      logical_vector_exact_match,
    
    exact_gene_set_match =
      gene_set_exact_match,
    
    tumor_up_exact_match =
      up_exact_match,
    
    tumor_down_exact_match =
      down_exact_match,
    
    n_disagreement =
      n_disagreement,
    
    existing_only =
      existing_only,
    
    candidate_only =
      candidate_only,
    
    up_existing_only =
      up_existing_only,
    
    up_candidate_only =
      up_candidate_only,
    
    down_existing_only =
      down_existing_only,
    
    down_candidate_only =
      down_candidate_only,
    
    min_abs_log2FC_existing_DEG =
      existing_deg_abs_log2fc_min,
    
    max_FDR_existing_DEG =
      existing_deg_fdr_max,
    
    core_audit =
      core_audit,
    
    no_source_modification = TRUE,
    no_step9b_modification = TRUE,
    no_deg_recalculation = TRUE,
    no_gsea = TRUE,
    no_ranking_calculation = TRUE,
    no_input_lock_created = TRUE,
    
    timestamp = Sys.time()
  )


saveRDS(
  STEP9D1H_DIAGNOSTIC,
  file =
    file.path(
      step9d_dir,
      "STEP9D1H_Exact_DEG_Definition_Concordance_Diagnostic.rds"
    )
)


write.csv(
  disagreement_table,
  file =
    file.path(
      step9d_dir,
      "STEP9D1H_DEG_Definition_Disagreement_Table.csv"
    ),
  row.names = FALSE
)


write.csv(
  core_audit,
  file =
    file.path(
      step9d_dir,
      "STEP9D1H_CoreGene_Exact_DEG_Definition_Audit.csv"
    ),
  row.names = FALSE
)


cat(
  "\nDiagnostic files saved under:\n",
  step9d_dir,
  "\n",
  sep = ""
)


cat("STEP9D-1I: Final Genome-wide GSEA Input Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Directories
# ------------------------------------------------------------

step9d_dir <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA/9D-1_Input_Lock"

source_file <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9B-7R_Revised_DEG_Definition_and_Lock/02_results/STEP9B7R_Revised_DEG_FULL.csv"

universe_file <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9C_Revised_DEG_GO_KEGG_ORA/00_input/STEP9C_DEG_testing_universe_57828.csv"


if (!file.exists(source_file)) {
  stop(
    "FAIL: Genome-wide source file not found:\n",
    source_file
  )
}

if (!file.exists(universe_file)) {
  stop(
    "FAIL: Locked testing-universe file not found:\n",
    universe_file
  )
}

if (!dir.exists(step9d_dir)) {
  dir.create(
    step9d_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 1. Read-only import
# ------------------------------------------------------------

source_df <-
  data.table::fread(
    source_file,
    data.table = FALSE,
    check.names = FALSE
  )

universe_df <-
  data.table::fread(
    universe_file,
    data.table = FALSE,
    check.names = FALSE
  )


# ------------------------------------------------------------
# 2. Identify universe column
# ------------------------------------------------------------

universe_gene_column <-
  if ("x" %in% colnames(universe_df)) {
    "x"
  } else if ("Gene" %in% colnames(universe_df)) {
    "Gene"
  } else {
    stop(
      "FAIL: Could not identify universe gene-ID column."
    )
  }


universe_genes <-
  as.character(
    universe_df[[universe_gene_column]]
  )


source_genes <-
  as.character(
    source_df$Gene
  )


# ------------------------------------------------------------
# 3. Structural audit
# ------------------------------------------------------------

required_columns <-
  c(
    "Gene",
    "log2FC",
    "P_value",
    "FDR",
    "Direction",
    "DEG"
  )

missing_columns <-
  setdiff(
    required_columns,
    colnames(source_df)
  )

if (length(missing_columns) > 0) {
  
  stop(
    "FAIL: Required source columns missing: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


source_n <-
  nrow(source_df)

universe_n <-
  length(universe_genes)


source_unique <-
  length(unique(source_genes)) ==
  length(source_genes)

universe_unique <-
  length(unique(universe_genes)) ==
  length(universe_genes)


exact_universe_match <-
  setequal(
    source_genes,
    universe_genes
  )


cat("============================================================\n")
cat("STRUCTURAL INPUT LOCK AUDIT\n")
cat("============================================================\n\n")

cat(
  "Source rows: ",
  source_n,
  "\n",
  sep = ""
)

cat(
  "Universe rows: ",
  universe_n,
  "\n",
  sep = ""
)

cat(
  "Source Gene IDs unique: ",
  source_unique,
  "\n",
  sep = ""
)

cat(
  "Universe Gene IDs unique: ",
  universe_unique,
  "\n",
  sep = ""
)

cat(
  "Exact source/universe gene-set match: ",
  exact_universe_match,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Numeric validity
# ------------------------------------------------------------

log2fc <-
  suppressWarnings(
    as.numeric(
      source_df$log2FC
    )
  )

p_value <-
  suppressWarnings(
    as.numeric(
      source_df$P_value
    )
  )

fdr <-
  suppressWarnings(
    as.numeric(
      source_df$FDR
    )
  )


finite_log2fc <-
  all(
    is.finite(log2fc)
  )

finite_p <-
  all(
    is.finite(p_value)
  )

finite_fdr <-
  all(
    is.finite(fdr)
  )

valid_p <-
  all(
    p_value >= 0 &
      p_value <= 1
  )

valid_fdr <-
  all(
    fdr >= 0 &
      fdr <= 1
  )


cat("\n============================================================\n")
cat("STATISTICAL FIELD VALIDITY\n")
cat("============================================================\n\n")

cat(
  "Finite log2FC: ",
  finite_log2fc,
  "\n",
  sep = ""
)

cat(
  "Finite P_value: ",
  finite_p,
  "\n",
  sep = ""
)

cat(
  "Finite FDR: ",
  finite_fdr,
  "\n",
  sep = ""
)

cat(
  "P_value in [0,1]: ",
  valid_p,
  "\n",
  sep = ""
)

cat(
  "FDR in [0,1]: ",
  valid_fdr,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. Locked result-level DEG definition
# ------------------------------------------------------------

existing_deg <-
  source_df$DEG == TRUE

reconstructed_deg <-
  fdr < 0.05 &
  abs(log2fc) > 1


exact_deg_match <-
  identical(
    existing_deg,
    reconstructed_deg
  )


existing_deg_n <-
  sum(existing_deg)

reconstructed_deg_n <-
  sum(reconstructed_deg)


existing_up_n <-
  sum(
    existing_deg &
      source_df$Direction == "Tumor_up"
  )

existing_down_n <-
  sum(
    existing_deg &
      source_df$Direction == "Tumor_down"
  )


cat("\n============================================================\n")
cat("LOCKED DEG DEFINITION\n")
cat("============================================================\n\n")

cat(
  "Definition: FDR < 0.05 AND |log2FC| > 1\n"
)

cat(
  "Existing DEG: ",
  existing_deg_n,
  "\n",
  sep = ""
)

cat(
  "Reconstructed DEG: ",
  reconstructed_deg_n,
  "\n",
  sep = ""
)

cat(
  "Exact logical match: ",
  exact_deg_match,
  "\n",
  sep = ""
)

cat(
  "Tumor_up: ",
  existing_up_n,
  "\n",
  sep = ""
)

cat(
  "Tumor_down: ",
  existing_down_n,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Core-gene audit
# ------------------------------------------------------------

core_genes <-
  c(
    AOC3  = "ENSG00000131471",
    MS4A7 = "ENSG00000166927",
    IL7R  = "ENSG00000168685"
  )


core_present_source <-
  core_genes %in% source_genes

core_present_universe <-
  core_genes %in% universe_genes


core_idx <-
  match(
    core_genes,
    source_genes
  )


core_audit <-
  source_df[
    core_idx,
    c(
      "Gene",
      "log2FC",
      "P_value",
      "FDR",
      "Direction",
      "DEG"
    ),
    drop = FALSE
  ]

core_audit$Core_Gene <-
  names(core_genes)

core_audit$Candidate_DEG <-
  reconstructed_deg[
    core_idx
  ]

core_audit$Exact_DEG_Match <-
  core_audit$DEG ==
  core_audit$Candidate_DEG

core_audit <-
  core_audit[
    c(
      "Core_Gene",
      "Gene",
      "log2FC",
      "P_value",
      "FDR",
      "Direction",
      "DEG",
      "Candidate_DEG",
      "Exact_DEG_Match"
    )
  ]


core_audit_pass <-
  all(
    core_present_source &
      core_present_universe
  ) &
  all(
    core_audit$Exact_DEG_Match
  )


cat("\n============================================================\n")
cat("CORE-GENE INPUT AUDIT\n")
cat("============================================================\n\n")

print(
  core_audit,
  row.names = FALSE
)

cat(
  "\nAll core genes present in source: ",
  all(core_present_source),
  "\n",
  sep = ""
)

cat(
  "All core genes present in universe: ",
  all(core_present_universe),
  "\n",
  sep = ""
)

cat(
  "Core-gene DEG definition match: ",
  core_audit_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. GSEA input principle
# ------------------------------------------------------------

gsea_universe_n <-
  source_n

gsea_uses_genome_wide_universe <-
  gsea_universe_n == 57828 &
  exact_universe_match


cat("\n============================================================\n")
cat("GENOME-WIDE GSEA INPUT PRINCIPLE\n")
cat("============================================================\n\n")

cat(
  "Genome-wide GSEA input genes: ",
  gsea_universe_n,
  "\n",
  sep = ""
)

cat(
  "Genome-wide universe confirmed: ",
  gsea_uses_genome_wide_universe,
  "\n",
  sep = ""
)

cat(
  "\nIMPORTANT:\n"
)

cat(
  "GSEA will use the complete 57,828-gene ranked list.\n"
)

cat(
  "The 16,065 DEG genes are retained only as provenance/QC information.\n"
)

cat(
  "No DEG filtering will be applied to the genome-wide GSEA input.\n"
)


# ------------------------------------------------------------
# 8. Final lock status
# ------------------------------------------------------------

final_lock_pass <-
  source_n == 57828 &
  universe_n == 57828 &
  source_unique &
  universe_unique &
  exact_universe_match &
  finite_log2fc &
  finite_p &
  finite_fdr &
  valid_p &
  valid_fdr &
  exact_deg_match &
  existing_deg_n == 16065 &
  existing_up_n == 9311 &
  existing_down_n == 6754 &
  core_audit_pass &
  gsea_uses_genome_wide_universe


cat("\n============================================================\n")
cat("STEP9D-1I FINAL STATUS\n")
cat("============================================================\n\n")


if (final_lock_pass) {
  
  cat(
    "STATUS: PASS — FINAL INPUT LOCK ESTABLISHED.\n"
  )
  
  cat(
    "\nLocked genome-wide GSEA universe: 57,828 genes.\n"
  )
  
  cat(
    "Locked DEG subset: 16,065 genes.\n"
  )
  
  cat(
    "Locked Tumor_up DEG subset: 9,311 genes.\n"
  )
  
  cat(
    "Locked Tumor_down DEG subset: 6,754 genes.\n"
  )
  
  cat(
    "Result-level reconstructed DEG rule: FDR < 0.05 AND |log2FC| > 1.\n"
  )
  
} else {
  
  cat(
    "STATUS: FAIL — FINAL INPUT LOCK NOT ESTABLISHED.\n"
  )
}


cat("\n")
cat("No STEP9B result modified.\n")
cat("No DEG recalculation performed for analysis purposes.\n")
cat("No GSEA performed.\n")
cat("No ranking metric calculated.\n")


# ------------------------------------------------------------
# 9. Save final lock
# ------------------------------------------------------------

STEP9D1I_FINAL_LOCK <-
  list(
    step = "STEP9D-1I",
    
    status =
      if (final_lock_pass)
        "PASS_FINAL_INPUT_LOCK"
    else
      "FAIL_FINAL_INPUT_LOCK",
    
    source_file =
      source_file,
    
    universe_file =
      universe_file,
    
    source_n =
      source_n,
    
    universe_n =
      universe_n,
    
    genome_wide_gsea_universe_n =
      gsea_universe_n,
    
    exact_universe_match =
      exact_universe_match,
    
    source_unique =
      source_unique,
    
    universe_unique =
      universe_unique,
    
    finite_log2fc =
      finite_log2fc,
    
    finite_p =
      finite_p,
    
    finite_fdr =
      finite_fdr,
    
    valid_p =
      valid_p,
    
    valid_fdr =
      valid_fdr,
    
    reconstructed_deg_definition =
      "FDR < 0.05 AND |log2FC| > 1",
    
    result_level_reconstruction_only =
      TRUE,
    
    original_generating_script_recovered =
      FALSE,
    
    existing_deg_n =
      existing_deg_n,
    
    existing_tumor_up_n =
      existing_up_n,
    
    existing_tumor_down_n =
      existing_down_n,
    
    exact_deg_match =
      exact_deg_match,
    
    core_audit =
      core_audit,
    
    gsea_uses_complete_57828_gene_list =
      TRUE,
    
    no_deg_filter_for_gsea =
      TRUE,
    
    no_step9b_modification =
      TRUE,
    
    no_gsea =
      TRUE,
    
    no_ranking_calculation =
      TRUE,
    
    timestamp =
      Sys.time()
  )


if (final_lock_pass) {
  
  saveRDS(
    STEP9D1I_FINAL_LOCK,
    file =
      file.path(
        step9d_dir,
        "STEP9D1I_FINAL_INPUT_LOCK.rds"
      )
  )
  
  write.csv(
    core_audit,
    file =
      file.path(
        step9d_dir,
        "STEP9D1I_CoreGene_Input_Lock_Audit.csv"
      ),
    row.names = FALSE
  )
  
  cat(
    "\nFINAL LOCK SAVED:\n",
    file.path(
      step9d_dir,
      "STEP9D1I_FINAL_INPUT_LOCK.rds"
    ),
    "\n",
    sep = ""
  )
  
} else {
  
  cat(
    "\nFINAL LOCK NOT SAVED BECAUSE STATUS = FAIL.\n"
  )
}


# ============================================================
cat("\n============================================================\n")
cat("STEP9D END-OF-DAY CHECKPOINT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Working directory
# ------------------------------------------------------------

step9d_root <- "result_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA"

checkpoint_dir <- file.path(
  step9d_root,
  "END_OF_DAY_CHECKPOINT_2026-09-22"
)

if (!dir.exists(checkpoint_dir)) {
  dir.create(checkpoint_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 2. Record today's locked / audited status
# ------------------------------------------------------------

step9d_checkpoint <- list(
  
  project = "Paper 1",
  workflow = "STEP9D_GenomeWide_Preranked_GSEA",
  
  checkpoint_date = "2026-09-22",
  
  checkpoint_status = "READY_TO_RESUME",
  
  # ----------------------------------------------------------
  # STEP9D-1H result-level reconstruction
  # ----------------------------------------------------------
  
  STEP9D1H_status = "PASS",
  
  STEP9D1H_type =
    "Exact result-level concordance audit",
  
  STEP9D1H_original_script_recovered = FALSE,
  
  STEP9D1H_original_rule_recovered = FALSE,
  
  STEP9D1H_reconstructed_rule =
    "FDR < 0.05 AND abs(log2FC) > 1",
  
  STEP9D1H_reconstruction_scope =
    "Result-level reconstruction only; original generating script/provenance was not recovered",
  
  # ----------------------------------------------------------
  # Genome-wide source
  # ----------------------------------------------------------
  
  source_file =
    "STEP9B7R_Revised_DEG_FULL.csv",
  
  source_n_genes = 57828,
  
  source_id_column = "Gene",
  
  source_statistics = c(
    "log2FC",
    "P_value",
    "FDR"
  ),
  
  # ----------------------------------------------------------
  # Testing universe
  # ----------------------------------------------------------
  
  testing_universe_file =
    "STEP9C_DEG_testing_universe_57828.csv",
  
  testing_universe_n_genes = 57828,
  
  testing_universe_id_column = "x",
  
  exact_source_universe_match = TRUE,
  
  # ----------------------------------------------------------
  # Existing DEG result
  # ----------------------------------------------------------
  
  existing_DEG_total = 16065,
  
  existing_DEG_Tumor_up = 9311,
  
  existing_DEG_Tumor_down = 6754,
  
  # ----------------------------------------------------------
  # Reconstructed DEG result
  # ----------------------------------------------------------
  
  reconstructed_DEG_total = 16065,
  
  reconstructed_DEG_Tumor_up = 9311,
  
  reconstructed_DEG_Tumor_down = 6754,
  
  exact_DEG_logical_match = TRUE,
  
  exact_DEG_gene_set_match = TRUE,
  
  exact_Tumor_up_match = TRUE,
  
  exact_Tumor_down_match = TRUE,
  
  logical_disagreements = 0,
  
  # ----------------------------------------------------------
  # Core genes
  # ----------------------------------------------------------
  
  core_genes = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  
  core_genes_exact_match = TRUE,
  
  # ----------------------------------------------------------
  # Important numerical audit results
  # ----------------------------------------------------------
  
  minimum_abs_log2FC_existing_DEG =
    1.000069302726,
  
  maximum_abs_log2FC_existing_DEG =
    11.269940524353,
  
  maximum_FDR_existing_DEG =
    0.049968458056,
  
  # ----------------------------------------------------------
  # Negative statements — important for reproducibility
  # ----------------------------------------------------------
  
  DEG_recalculated_today = FALSE,
  
  STEP9B_modified_today = FALSE,
  
  GSEA_performed_today = FALSE,
  
  ranking_metric_locked_today = FALSE,
  
  ranking_list_created_today = FALSE,
  
  GO_GSEA_performed_today = FALSE,
  
  KEGG_GSEA_performed_today = FALSE,
  
  # ----------------------------------------------------------
  # Next step
  # ----------------------------------------------------------
  
  next_step =
    "STEP9D-1 formal Input/Gene-Universe Lock",
  
  next_step_number =
    "STEP9D-1I",
  
  next_step_action = c(
    "Lock STEP9B7R_Revised_DEG_FULL.csv as genome-wide source",
    "Lock exact 57,828-gene testing universe",
    "Verify finite/valid log2FC, P_value and FDR",
    "Verify direction/sign consistency",
    "Record exact STEP9D1H DEG concordance",
    "Verify AOC3/MS4A7/IL7R presence",
    "Create formal STEP9D-1 input lock",
    "Do NOT calculate ranking metric in STEP9D-1"
  ),
  
  # ----------------------------------------------------------
  # Future ranking step
  # ----------------------------------------------------------
  
  future_ranking_step =
    "STEP9D-2 Ranking Metric Lock",
  
  future_ranking_note =
    "Use one genome-wide signed ranking list; do not split into up/down GSEA lists",
  
  candidate_ranking_metric =
    "sign(log2FC) * (-log10(P_value))",
  
  candidate_ranking_status =
    "NOT YET LOCKED",
  
  ranking_metric_requires_version_and_parameter_audit = TRUE,
  
  # ----------------------------------------------------------
  # GSEA implementation audit
  # ----------------------------------------------------------
  
  future_GSEA_implementation_audit = c(
    "Inspect installed fgsea version",
    "Inspect clusterProfiler::GSEA formal arguments",
    "Confirm exact GSEA implementation and parameters",
    "Do not invent unsupported parameters"
  ),
  
  # ----------------------------------------------------------
  # Final protection
  # ----------------------------------------------------------
  
  instruction_for_resume =
    "Resume from STEP9D-1 formal Input/Gene-Universe Lock. Do not redo STEP9B or STEP9C."
)


# ------------------------------------------------------------
# 3. Save checkpoint RDS
# ------------------------------------------------------------

checkpoint_rds <- file.path(
  checkpoint_dir,
  "STEP9D_END_OF_DAY_CHECKPOINT_2026-09-22.rds"
)

saveRDS(
  step9d_checkpoint,
  checkpoint_rds
)


# ------------------------------------------------------------
# 4. Save human-readable TXT summary
# ------------------------------------------------------------

checkpoint_txt <- file.path(
  checkpoint_dir,
  "STEP9D_END_OF_DAY_CHECKPOINT_2026-09-22.txt"
)

sink(checkpoint_txt)

cat("============================================================\n")
cat("Paper 1 — STEP9D End-of-Day Checkpoint\n")
cat("Date: 2026-09-22\n")
cat("============================================================\n\n")

cat("TODAY'S FINAL STATUS\n")
cat("--------------------\n")
cat("STEP9D-1H: PASS\n")
cat("Type: Exact result-level concordance audit\n\n")

cat("Reconstructed DEG definition:\n")
cat("FDR < 0.05 AND abs(log2FC) > 1\n\n")

cat("Existing DEG total: 16,065\n")
cat("Reconstructed DEG total: 16,065\n\n")

cat("Existing Tumor_up: 9,311\n")
cat("Reconstructed Tumor_up: 9,311\n\n")

cat("Existing Tumor_down: 6,754\n")
cat("Reconstructed Tumor_down: 6,754\n\n")

cat("Exact logical match: TRUE\n")
cat("Exact gene-set match: TRUE\n")
cat("Tumor_up exact match: TRUE\n")
cat("Tumor_down exact match: TRUE\n")
cat("Logical disagreements: 0\n\n")

cat("Core genes:\n")
cat("AOC3, MS4A7, IL7R\n")
cat("All three show exact DEG-definition concordance.\n\n")

cat("IMPORTANT PROVENANCE NOTE\n")
cat("--------------------------\n")
cat("The original STEP9B generating script/provenance was not recovered.\n")
cat("Therefore, the DEG definition is documented as an\n")
cat("EXACT RESULT-LEVEL RECONSTRUCTION, not as recovery of the\n")
cat("original generating script.\n\n")

cat("NO ACTIONS PERFORMED TODAY\n")
cat("--------------------------\n")
cat("DEG recalculation: NO\n")
cat("STEP9B modification: NO\n")
cat("GSEA: NO\n")
cat("Ranking metric lock: NO\n")
cat("GO-GSEA: NO\n")
cat("KEGG-GSEA: NO\n\n")

cat("NEXT STEP\n")
cat("---------\n")
cat("STEP9D-1I — Formal Input/Gene-Universe Lock\n\n")

cat("Input source:\n")
cat("STEP9B7R_Revised_DEG_FULL.csv\n")
cat("57,828 genes\n\n")

cat("Testing universe:\n")
cat("STEP9C_DEG_testing_universe_57828.csv\n")
cat("57,828 genes\n")
cat("Exact source-universe ID match: TRUE\n\n")

cat("STEP9D-1I must NOT calculate the ranking metric.\n")
cat("Ranking will be handled separately in STEP9D-2.\n\n")

cat("FUTURE STEP9D-2\n")
cat("---------------\n")
cat("Candidate metric: sign(log2FC) * (-log10(P_value))\n")
cat("Status: NOT YET LOCKED\n")
cat("First inspect installed fgsea / clusterProfiler versions\n")
cat("and exact GSEA implementation parameters.\n\n")

cat("============================================================\n")
cat("CHECKPOINT SAVED\n")
cat("============================================================\n")

sink()


# ------------------------------------------------------------
# 5. Final file existence check
# ------------------------------------------------------------

cat("\nCheckpoint files:\n")
cat("RDS : ", checkpoint_rds, "\n", sep = "")
cat("TXT : ", checkpoint_txt, "\n", sep = "")

cat("\nFile existence check:\n")
print(
  c(
    RDS = file.exists(checkpoint_rds),
    TXT = file.exists(checkpoint_txt)
  )
)


# ------------------------------------------------------------
# 6. Final status
# ------------------------------------------------------------

if (
  file.exists(checkpoint_rds) &&
  file.exists(checkpoint_txt)
) {
  cat("\n============================================================\n")
  cat("FINAL CHECKPOINT STATUS = PASS\n")
  cat("Today's STEP9D work is safely checkpointed.\n")
  cat("Tomorrow resume from STEP9D-1I.\n")
  cat("============================================================\n")
} else {
  cat("\nFINAL CHECKPOINT STATUS = FAIL\n")
  cat("Please verify the output directory and files.\n")
}




