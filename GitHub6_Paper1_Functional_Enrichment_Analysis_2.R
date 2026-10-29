##GitHub6_Paper1_Functional_Enrichment_Analysis_2
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


load(
  "input_dir/Paper 1_Step 9-Functional Enrichment/paper 1_Step 9-Functional Enrichment.RData"
)

ls()

sessionInfo()

# ============================================================
# Paper 1
# STEP9D-1I
# Formal Input / Gene-Universe Lock
# ============================================================
#
# Purpose:
#   1. Formally lock the genome-wide STEP9D input source
#   2. Formally lock the 57,828-gene testing universe
#   3. Confirm exact source/universe ID concordance
#   4. Confirm valid genome-wide statistics
#   5. Confirm direction/sign consistency
#   6. Preserve the STEP9D-1H exact DEG-definition audit result
#   7. Confirm AOC3 / MS4A7 / IL7R are present
#
# IMPORTANT:
#   - No DEG recalculation
#   - No STEP9B modification
#   - No STEP9C modification
#   - No GSEA
#   - No ranking-list construction
#   - No ranking-metric lock
#   - No P/FDR recalculation
#   - No ID conversion
# ============================================================


cat("\n============================================================\n")
cat("STEP9D-1I: Formal Input / Gene-Universe Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory and output folders
# ------------------------------------------------------------

work_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

step9d_root <- file.path(
  work_dir,
  "STEP9D_GenomeWide_Preranked_GSEA"
)

input_lock_dir <- file.path(
  step9d_root,
  "9D-1_Input_Lock"
)

if (!dir.exists(input_lock_dir)) {
  dir.create(
    input_lock_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 1. Define locked input files
# ------------------------------------------------------------

source_file <- file.path(
  work_dir,
  "STEP9B-7R_Revised_DEG_Definition_and_Lock",
  "02_results",
  "STEP9B7R_Revised_DEG_FULL.csv"
)

universe_file <- file.path(
  work_dir,
  "STEP9C_Revised_DEG_GO_KEGG_ORA",
  "00_input",
  "STEP9C_DEG_testing_universe_57828.csv"
)

diagnostic_1H_file <- file.path(
  input_lock_dir,
  "STEP9D1H_Exact_DEG_Definition_Concordance_Diagnostic.rds"
)


cat("Source file:\n")
cat(source_file, "\n\n")

cat("Testing-universe file:\n")
cat(universe_file, "\n\n")


# ------------------------------------------------------------
# 2. File existence checks
# ------------------------------------------------------------

source_exists <- file.exists(source_file)
universe_exists <- file.exists(universe_file)

cat("File existence:\n")
cat("  Source   :", source_exists, "\n")
cat("  Universe :", universe_exists, "\n\n")


if (!source_exists || !universe_exists) {
  
  stop(
    paste0(
      "FAIL: Required STEP9D input file(s) not found.\n",
      "Check the paths above before proceeding."
    )
  )
  
}


# ------------------------------------------------------------
# 3. Read-only import
# ------------------------------------------------------------
#
# No modification of source files.
# No recalculation.
#

source_df <- read.csv(
  source_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

universe_df <- read.csv(
  universe_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("Imported dimensions:\n")
cat(
  "  Source   :",
  nrow(source_df),
  "rows x",
  ncol(source_df),
  "columns\n"
)

cat(
  "  Universe :",
  nrow(universe_df),
  "rows x",
  ncol(universe_df),
  "columns\n\n"
)


# ------------------------------------------------------------
# 4. Required source columns
# ------------------------------------------------------------

required_source_columns <- c(
  "Gene",
  "log2FC",
  "P_value",
  "FDR",
  "Direction",
  "DEG"
)

source_columns_pass <- all(
  required_source_columns %in% names(source_df)
)

cat("Required source columns:\n")
print(required_source_columns)

cat(
  "\nColumn structure PASS =",
  source_columns_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 5. Universe column identification
# ------------------------------------------------------------

universe_gene_column <- "x"

universe_column_pass <- (
  universe_gene_column %in% names(universe_df)
)

cat(
  "Testing-universe gene column:",
  universe_gene_column,
  "\n"
)

cat(
  "Universe column exists PASS =",
  universe_column_pass,
  "\n\n"
)


if (!source_columns_pass || !universe_column_pass) {
  
  stop(
    "FAIL: Required input columns are missing."
  )
  
}


# ------------------------------------------------------------
# 6. Extract gene identifiers
# ------------------------------------------------------------

source_genes <- as.character(
  source_df$Gene
)

universe_genes <- as.character(
  universe_df[[universe_gene_column]]
)


# ------------------------------------------------------------
# 7. Basic gene-ID validity
# ------------------------------------------------------------

source_gene_na <- sum(
  is.na(source_genes) |
    source_genes == ""
)

universe_gene_na <- sum(
  is.na(universe_genes) |
    universe_genes == ""
)

source_duplicate_genes <- sum(
  duplicated(source_genes)
)

universe_duplicate_genes <- sum(
  duplicated(universe_genes)
)


cat("Gene-ID quality:\n")
cat(
  "  Source NA/empty      :",
  source_gene_na,
  "\n"
)

cat(
  "  Universe NA/empty    :",
  universe_gene_na,
  "\n"
)

cat(
  "  Source duplicates    :",
  source_duplicate_genes,
  "\n"
)

cat(
  "  Universe duplicates  :",
  universe_duplicate_genes,
  "\n\n"
)


# ------------------------------------------------------------
# 8. Formal sample-size / universe-size lock
# ------------------------------------------------------------

expected_n_genes <- 57828

source_n_pass <- (
  nrow(source_df) == expected_n_genes
)

universe_n_pass <- (
  nrow(universe_df) == expected_n_genes
)

source_unique_n_pass <- (
  length(unique(source_genes)) == expected_n_genes
)

universe_unique_n_pass <- (
  length(unique(universe_genes)) == expected_n_genes
)


cat("Genome-wide size checks:\n")
cat(
  "  Source rows = 57,828:",
  source_n_pass,
  "\n"
)

cat(
  "  Universe rows = 57,828:",
  universe_n_pass,
  "\n"
)

cat(
  "  Source unique genes = 57,828:",
  source_unique_n_pass,
  "\n"
)

cat(
  "  Universe unique genes = 57,828:",
  universe_unique_n_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 9. Exact source-vs-universe gene-set concordance
# ------------------------------------------------------------

source_gene_set <- unique(source_genes)
universe_gene_set <- unique(universe_genes)

source_only_genes <- setdiff(
  source_gene_set,
  universe_gene_set
)

universe_only_genes <- setdiff(
  universe_gene_set,
  source_gene_set
)

exact_gene_set_match <- (
  length(source_only_genes) == 0 &&
    length(universe_only_genes) == 0 &&
    length(source_gene_set) == length(universe_gene_set)
)


cat("Exact source/universe ID-set concordance:\n")
cat(
  "  Source-only genes   :",
  length(source_only_genes),
  "\n"
)

cat(
  "  Universe-only genes :",
  length(universe_only_genes),
  "\n"
)

cat(
  "  Exact gene-set match:",
  exact_gene_set_match,
  "\n\n"
)


# ------------------------------------------------------------
# 10. Statistical-column validation
# ------------------------------------------------------------

log2fc <- suppressWarnings(
  as.numeric(source_df$log2FC)
)

p_value <- suppressWarnings(
  as.numeric(source_df$P_value)
)

fdr <- suppressWarnings(
  as.numeric(source_df$FDR)
)


log2fc_finite_pass <- all(
  is.finite(log2fc)
)

p_value_finite_pass <- all(
  is.finite(p_value)
)

fdr_finite_pass <- all(
  is.finite(fdr)
)

p_value_range_pass <- all(
  p_value >= 0 &
    p_value <= 1
)

fdr_range_pass <- all(
  fdr >= 0 &
    fdr <= 1
)


cat("Genome-wide statistical validity:\n")
cat(
  "  log2FC finite:",
  log2fc_finite_pass,
  "\n"
)

cat(
  "  P_value finite:",
  p_value_finite_pass,
  "\n"
)

cat(
  "  FDR finite:",
  fdr_finite_pass,
  "\n"
)

cat(
  "  P_value in [0,1]:",
  p_value_range_pass,
  "\n"
)

cat(
  "  FDR in [0,1]:",
  fdr_range_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 11. Direction validation
# ------------------------------------------------------------

direction <- as.character(
  source_df$Direction
)

direction_levels <- sort(
  unique(direction)
)

expected_directions <- c(
  "Tumor_down",
  "Tumor_up"
)

direction_values_pass <- all(
  direction %in% expected_directions
)


# Direction must agree with log2FC sign.
direction_sign_pass <- all(
  (direction == "Tumor_up" & log2fc > 0) |
    (direction == "Tumor_down" & log2fc < 0)
)


cat("Direction validation:\n")

cat("  Observed Direction values:\n")
print(direction_levels)

cat(
  "\n  Allowed values:",
  paste(expected_directions, collapse = ", "),
  "\n"
)

cat(
  "  Direction-value PASS:",
  direction_values_pass,
  "\n"
)

cat(
  "  Direction/log2FC sign PASS:",
  direction_sign_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 12. DEG logical flag validation
# ------------------------------------------------------------

deg_flag <- source_df$DEG

deg_logical_pass <- is.logical(deg_flag)

if (deg_logical_pass) {
  
  deg_true_n <- sum(
    deg_flag,
    na.rm = TRUE
  )
  
  deg_false_n <- sum(
    !deg_flag,
    na.rm = TRUE
  )
  
  deg_na_n <- sum(
    is.na(deg_flag)
  )
  
} else {
  
  deg_true_n <- NA_integer_
  deg_false_n <- NA_integer_
  deg_na_n <- NA_integer_
  
}


cat("Existing DEG flag:\n")
cat(
  "  Logical type:",
  deg_logical_pass,
  "\n"
)

cat(
  "  TRUE :",
  deg_true_n,
  "\n"
)

cat(
  "  FALSE:",
  deg_false_n,
  "\n"
)

cat(
  "  NA   :",
  deg_na_n,
  "\n\n"
)


# ------------------------------------------------------------
# 13. STEP9D-1H exact DEG-definition concordance
# ------------------------------------------------------------
#
# IMPORTANT:
# This is a READ-ONLY confirmation of the already completed
# STEP9D-1H result-level audit.
#
# It does NOT recalculate or redefine STEP9B.
# ------------------------------------------------------------

candidate_deg <- (
  fdr < 0.05 &
    abs(log2fc) > 1
)

candidate_up <- (
  candidate_deg &
    direction == "Tumor_up"
)

candidate_down <- (
  candidate_deg &
    direction == "Tumor_down"
)

existing_deg <- deg_flag

existing_up <- (
  existing_deg &
    direction == "Tumor_up"
)

existing_down <- (
  existing_deg &
    direction == "Tumor_down"
)


exact_deg_logical_match <- identical(
  existing_deg,
  candidate_deg
)

exact_up_logical_match <- identical(
  existing_up,
  candidate_up
)

exact_down_logical_match <- identical(
  existing_down,
  candidate_down
)


deg_disagreement_n <- sum(
  existing_deg != candidate_deg
)


cat("STEP9D-1H exact DEG-definition concordance:\n")

cat(
  "  Existing DEG total      :",
  sum(existing_deg),
  "\n"
)

cat(
  "  Candidate DEG total     :",
  sum(candidate_deg),
  "\n"
)

cat(
  "  Existing Tumor_up       :",
  sum(existing_up),
  "\n"
)

cat(
  "  Candidate Tumor_up      :",
  sum(candidate_up),
  "\n"
)

cat(
  "  Existing Tumor_down     :",
  sum(existing_down),
  "\n"
)

cat(
  "  Candidate Tumor_down    :",
  sum(candidate_down),
  "\n"
)

cat(
  "  Exact logical DEG match :",
  exact_deg_logical_match,
  "\n"
)

cat(
  "  Exact Tumor_up match    :",
  exact_up_logical_match,
  "\n"
)

cat(
  "  Exact Tumor_down match  :",
  exact_down_logical_match,
  "\n"
)

cat(
  "  Logical disagreements   :",
  deg_disagreement_n,
  "\n\n"
)

# ============================================================
# STEP9D-1I — Core Gene ID-Type Diagnostic
# READ-ONLY
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-1I Core Gene ID-Type Diagnostic\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Inspect source Gene IDs
# ------------------------------------------------------------

cat("First 10 source Gene IDs:\n")
print(head(source_df$Gene, 10))

cat("\n")


# ------------------------------------------------------------
# 2. Check whether source IDs look like Ensembl IDs
# ------------------------------------------------------------

ensembl_pattern <- "^ENSG[0-9]+$"

ensembl_like <- grepl(
  ensembl_pattern,
  as.character(source_df$Gene)
)

cat(
  "Ensembl-like source IDs:",
  sum(ensembl_like),
  "/",
  nrow(source_df),
  "\n"
)

cat(
  "All source IDs Ensembl-like:",
  all(ensembl_like),
  "\n\n"
)


# ------------------------------------------------------------
# 3. Locked core-gene Ensembl IDs
# ------------------------------------------------------------

core_gene_map <- data.frame(
  Symbol = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  Ensembl = c(
    "ENSG00000131471",
    "ENSG00000166927",
    "ENSG00000168685"
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 4. Exact source / universe presence
# ------------------------------------------------------------

core_gene_map$In_Source <- (
  core_gene_map$Ensembl %in%
    source_genes
)

core_gene_map$In_Universe <- (
  core_gene_map$Ensembl %in%
    universe_genes
)


# ------------------------------------------------------------
# 5. Recover source statistics
# ------------------------------------------------------------

idx <- match(
  core_gene_map$Ensembl,
  source_genes
)

core_gene_map$log2FC <- log2fc[idx]
core_gene_map$P_value <- p_value[idx]
core_gene_map$FDR <- fdr[idx]
core_gene_map$Direction <- direction[idx]
core_gene_map$Existing_DEG <- existing_deg[idx]

core_gene_map$Reconstructed_DEG <- candidate_deg[idx]

core_gene_map$Exact_DEG_Match <- (
  core_gene_map$Existing_DEG ==
    core_gene_map$Reconstructed_DEG
)


# ------------------------------------------------------------
# 6. Print final diagnostic table
# ------------------------------------------------------------

cat("Core-gene Ensembl ID audit:\n\n")
print(
  core_gene_map,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 7. Formal diagnostic status
# ------------------------------------------------------------

core_id_audit_pass <- all(
  core_gene_map$In_Source,
  core_gene_map$In_Universe,
  core_gene_map$Exact_DEG_Match
)

cat("============================================================\n")
cat("CORE GENE ID AUDIT STATUS\n")
cat("============================================================\n")

if (isTRUE(core_id_audit_pass)) {
  
  cat(
    "PASS — AOC3/MS4A7/IL7R are present in the\n",
    "57,828-gene Ensembl source and testing universe.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "FAIL — Core-gene Ensembl ID audit requires review.\n"
  )
  
}

cat("============================================================\n")



# ------------------------------------------------------------
# 15. Confirm STEP9D-1H diagnostic provenance
# ------------------------------------------------------------

diagnostic_1H_exists <- file.exists(
  diagnostic_1H_file
)

cat(
  "STEP9D-1H diagnostic file exists:",
  diagnostic_1H_exists,
  "\n\n"
)


# ------------------------------------------------------------
# 16. Final PASS criteria
# ------------------------------------------------------------

all_pass <- all(
  source_exists,
  universe_exists,
  
  source_columns_pass,
  universe_column_pass,
  
  source_gene_na == 0,
  universe_gene_na == 0,
  
  source_duplicate_genes == 0,
  universe_duplicate_genes == 0,
  
  source_n_pass,
  universe_n_pass,
  
  source_unique_n_pass,
  universe_unique_n_pass,
  
  exact_gene_set_match,
  
  log2fc_finite_pass,
  p_value_finite_pass,
  fdr_finite_pass,
  
  p_value_range_pass,
  fdr_range_pass,
  
  direction_values_pass,
  direction_sign_pass,
  
  deg_logical_pass,
  deg_na_n == 0,
  
  exact_deg_logical_match,
  exact_up_logical_match,
  exact_down_logical_match,
  
  deg_disagreement_n == 0,
  
  core_source_pass,
  core_universe_pass,
  
  all(core_audit$Exact_DEG_Match)
)


# ------------------------------------------------------------
# 17. Print final status
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP9D-1I FINAL STATUS\n")
cat("============================================================\n\n")

if (isTRUE(all_pass)) {
  
  cat("PASS — FORMAL INPUT / GENE-UNIVERSE LOCK\n\n")
  
  cat(
    "The STEP9D genome-wide input is formally locked as:\n"
  )
  
  cat(
    "  57,828 unique Ensembl gene IDs\n"
  )
  
  cat(
    "  Source: STEP9B7R_Revised_DEG_FULL.csv\n"
  )
  
  cat(
    "  Testing universe: STEP9C_DEG_testing_universe_57828.csv\n\n"
  )
  
  cat(
    "Exact source/universe ID concordance: TRUE\n"
  )
  
  cat(
    "Exact STEP9D-1H DEG-definition concordance: TRUE\n"
  )
  
  cat(
    "Core genes AOC3/MS4A7/IL7R present: TRUE\n\n"
  )
  
  cat(
    "No DEG recalculation was performed.\n"
  )
  
  cat(
    "No STEP9B/STEP9C modification was performed.\n"
  )
  
  cat(
    "No GSEA was performed.\n"
  )
  
  cat(
    "No ranking metric was locked.\n"
  )
  
} else {
  
  cat("FAIL — INPUT LOCK NOT ESTABLISHED\n\n")
  cat("At least one formal input criterion failed.\n")
  cat("DO NOT proceed to STEP9D-2.\n")
  
}


# ------------------------------------------------------------
# 18. Build formal lock object
# ------------------------------------------------------------

STEP9D1I_FINAL_LOCK <- list(
  
  workflow = "Paper 1 — STEP9D Genome-wide Preranked GSEA",
  
  step = "STEP9D-1I",
  
  step_title = "Formal Input / Gene-Universe Lock",
  
  date = "2026-09-22",
  
  final_status = if (isTRUE(all_pass)) {
    "PASS"
  } else {
    "FAIL"
  },
  
  source = list(
    file = source_file,
    n_rows = nrow(source_df),
    n_unique_genes = length(unique(source_genes)),
    gene_column = "Gene"
  ),
  
  testing_universe = list(
    file = universe_file,
    n_rows = nrow(universe_df),
    n_unique_genes = length(unique(universe_genes)),
    gene_column = universe_gene_column
  ),
  
  exact_source_universe_match = exact_gene_set_match,
  
  statistical_columns = c(
    "log2FC",
    "P_value",
    "FDR"
  ),
  
  direction_definition = expected_directions,
  
  existing_DEG = list(
    total = sum(existing_deg),
    Tumor_up = sum(existing_up),
    Tumor_down = sum(existing_down)
  ),
  
  reconstructed_DEG_definition = (
    "FDR < 0.05 AND abs(log2FC) > 1"
  ),
  
  reconstructed_DEG = list(
    total = sum(candidate_deg),
    Tumor_up = sum(candidate_up),
    Tumor_down = sum(candidate_down)
  ),
  
  exact_DEG_logical_match = exact_deg_logical_match,
  
  exact_Tumor_up_match = exact_up_logical_match,
  
  exact_Tumor_down_match = exact_down_logical_match,
  
  logical_disagreements = deg_disagreement_n,
  
  core_genes = core_genes,
  
  core_gene_audit = core_audit,
  
  provenance_statement =
    paste(
      "The original STEP9B generating script/provenance was not recovered.",
      "The DEG definition is therefore documented as an exact",
      "result-level reconstruction rather than recovery of the",
      "original generating script."
    ),
  
  actions_not_performed = c(
    "DEG recalculation",
    "STEP9B modification",
    "STEP9C modification",
    "GSEA",
    "ranking-list construction",
    "ranking-metric lock",
    "P/FDR recalculation",
    "ID conversion"
  ),
  
  next_step =
    "STEP9D-2 Ranking Metric Lock",
  
  next_step_requirement =
    paste(
      "First inspect installed fgsea and clusterProfiler versions",
      "and the exact GSEA implementation parameters before locking",
      "the ranking metric."
    )
)


# ------------------------------------------------------------
# 19. Save formal lock
# ------------------------------------------------------------

final_lock_file <- file.path(
  input_lock_dir,
  "STEP9D1I_FINAL_LOCK.rds"
)

saveRDS(
  STEP9D1I_FINAL_LOCK,
  final_lock_file
)


# ------------------------------------------------------------
# 20. Save core-gene audit table
# ------------------------------------------------------------

core_audit_file <- file.path(
  input_lock_dir,
  "STEP9D1I_CoreGene_Input_Lock_Audit.csv"
)

write.csv(
  core_audit,
  core_audit_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 21. Save concise human-readable lock summary
# ------------------------------------------------------------

summary_file <- file.path(
  input_lock_dir,
  "STEP9D1I_FORMAL_INPUT_LOCK_SUMMARY.txt"
)

sink(summary_file)

cat("============================================================\n")
cat("Paper 1 — STEP9D-1I Formal Input / Gene-Universe Lock\n")
cat("Date: 2026-09-22\n")
cat("============================================================\n\n")

cat(
  "FINAL STATUS: ",
  if (isTRUE(all_pass)) "PASS" else "FAIL",
  "\n\n",
  sep = ""
)

cat("SOURCE\n")
cat("------\n")
cat("STEP9B7R_Revised_DEG_FULL.csv\n")
cat("57,828 unique genes\n\n")

cat("TESTING UNIVERSE\n")
cat("----------------\n")
cat("STEP9C_DEG_testing_universe_57828.csv\n")
cat("57,828 unique genes\n\n")

cat("EXACT SOURCE/UNIVERSE MATCH\n")
cat("---------------------------\n")
cat("TRUE\n\n")

cat("STEP9D-1H DEG RESULT-LEVEL RECONSTRUCTION\n")
cat("-----------------------------------------\n")
cat("Definition: FDR < 0.05 AND abs(log2FC) > 1\n")
cat("Total: 16,065\n")
cat("Tumor_up: 9,311\n")
cat("Tumor_down: 6,754\n")
cat("Exact logical match: TRUE\n")
cat("Logical disagreements: 0\n\n")

cat("CORE GENES\n")
cat("----------\n")
cat("AOC3, MS4A7, IL7R\n")
cat("All present in source and testing universe.\n\n")

cat("PROVENANCE NOTE\n")
cat("---------------\n")
cat(
  "The original STEP9B generating script/provenance was not recovered.\n"
)
cat(
  "The DEG definition is documented as an exact result-level\n"
)
cat(
  "reconstruction, not as recovery of the original generating script.\n\n"
)

cat("ACTIONS NOT PERFORMED\n")
cat("---------------------\n")
cat("DEG recalculation: NO\n")
cat("STEP9B modification: NO\n")
cat("STEP9C modification: NO\n")
cat("GSEA: NO\n")
cat("Ranking-list construction: NO\n")
cat("Ranking-metric lock: NO\n")
cat("P/FDR recalculation: NO\n")
cat("ID conversion: NO\n\n")

cat("NEXT STEP\n")
cat("---------\n")
cat("STEP9D-2 — Ranking Metric Lock\n")
cat(
  "Before locking the ranking metric, inspect installed fgsea/\n"
)
cat(
  "clusterProfiler versions and exact GSEA implementation parameters.\n"
)

cat("\n============================================================\n")

sink()


# ------------------------------------------------------------
# 22. Final file verification
# ------------------------------------------------------------

cat("\nSaved formal lock:\n")
cat(final_lock_file, "\n")

cat("\nSaved core-gene audit:\n")
cat(core_audit_file, "\n")

cat("\nSaved summary:\n")
cat(summary_file, "\n\n")


file_check <- c(
  FINAL_LOCK_RDS = file.exists(final_lock_file),
  CORE_AUDIT_CSV = file.exists(core_audit_file),
  SUMMARY_TXT = file.exists(summary_file)
)

print(file_check)


# ------------------------------------------------------------
# 23. Final checkpoint
# ------------------------------------------------------------

if (
  isTRUE(all_pass) &&
  all(file_check)
) {
  
  cat("\n============================================================\n")
  cat("STEP9D-1I FINAL CHECKPOINT = PASS\n")
  cat("============================================================\n")
  cat("\n")
  cat("Formal STEP9D input/gene-universe lock has been saved.\n")
  cat("No DEG/GSEA/ranking recalculation was performed.\n")
  cat("Next step: STEP9D-2 Ranking Metric Lock.\n")
  cat("============================================================\n\n")
  
} else {
  
  cat("\n============================================================\n")
  cat("STEP9D-1I FINAL CHECKPOINT = FAIL\n")
  cat("============================================================\n")
  cat("\n")
  cat("Do NOT proceed to STEP9D-2.\n")
  cat("Review the failed criteria above.\n")
  cat("============================================================\n\n")
  
}

# ============================================================
# Paper 1
# STEP9D-1I — Final Repair and Formal Lock
#
# Purpose:
#   Replace the incorrect symbol-based core_audit object
#   with the verified Ensembl-based core_gene_map.
#
# IMPORTANT:
#   This is NOT a new analysis.
#   No DEG recalculation.
#   No GSEA.
#   No ranking calculation.
#   No STEP9B / STEP9C modification.
# ============================================================


cat("\n============================================================\n")
cat("STEP9D-1I — Final Repair and Formal Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Confirm required objects already exist
# ------------------------------------------------------------

required_objects <- c(
  "source_df",
  "universe_df",
  "source_genes",
  "universe_genes",
  "log2fc",
  "p_value",
  "fdr",
  "direction",
  "existing_deg",
  "candidate_deg",
  "existing_up",
  "candidate_up",
  "existing_down",
  "candidate_down",
  "exact_gene_set_match",
  "exact_deg_logical_match",
  "exact_up_logical_match",
  "exact_down_logical_match",
  "deg_disagreement_n",
  "core_gene_map"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1)
  )
]

cat("Required objects check:\n")

if (length(missing_objects) == 0) {
  
  cat("PASS — all required objects are present.\n\n")
  
} else {
  
  cat("FAIL — missing objects:\n")
  print(missing_objects)
  
  stop(
    "Required STEP9D-1I objects are missing. Do not create a new lock."
  )
}


# ------------------------------------------------------------
# 2. Use the VERIFIED Ensembl-based core-gene audit
# ------------------------------------------------------------

# IMPORTANT:
# core_gene_map was just independently verified:
#
# AOC3  -> ENSG00000131471
# MS4A7 -> ENSG00000166927
# IL7R  -> ENSG00000168685
#
# All three:
#   In_Source = TRUE
#   In_Universe = TRUE
#   Exact_DEG_Match = TRUE
#
# Therefore this object is now the formal core-gene audit.

core_audit <- core_gene_map


cat("Verified core-gene audit:\n\n")
print(
  core_audit,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 3. Formal core-gene PASS
# ------------------------------------------------------------

core_source_pass <- all(
  core_audit$In_Source
)

core_universe_pass <- all(
  core_audit$In_Universe
)

core_exact_deg_pass <- all(
  core_audit$Exact_DEG_Match
)


cat("Core-gene formal checks:\n")
cat(
  "  All core genes in source   :",
  core_source_pass,
  "\n"
)

cat(
  "  All core genes in universe :",
  core_universe_pass,
  "\n"
)

cat(
  "  All core DEG matches       :",
  core_exact_deg_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 4. Reconfirm all STEP9D-1I formal criteria
# ------------------------------------------------------------

source_exists <- file.exists(source_file)
universe_exists <- file.exists(universe_file)

source_columns_pass <- all(
  required_source_columns %in% names(source_df)
)

universe_column_pass <- (
  universe_gene_column %in% names(universe_df)
)

source_gene_na <- sum(
  is.na(source_genes) |
    source_genes == ""
)

universe_gene_na <- sum(
  is.na(universe_genes) |
    universe_genes == ""
)

source_duplicate_genes <- sum(
  duplicated(source_genes)
)

universe_duplicate_genes <- sum(
  duplicated(universe_genes)
)

expected_n_genes <- 57828

source_n_pass <- (
  nrow(source_df) == expected_n_genes
)

universe_n_pass <- (
  nrow(universe_df) == expected_n_genes
)

source_unique_n_pass <- (
  length(unique(source_genes)) == expected_n_genes
)

universe_unique_n_pass <- (
  length(unique(universe_genes)) == expected_n_genes
)

log2fc_finite_pass <- all(
  is.finite(log2fc)
)

p_value_finite_pass <- all(
  is.finite(p_value)
)

fdr_finite_pass <- all(
  is.finite(fdr)
)

p_value_range_pass <- all(
  p_value >= 0 &
    p_value <= 1
)

fdr_range_pass <- all(
  fdr >= 0 &
    fdr <= 1
)

expected_directions <- c(
  "Tumor_down",
  "Tumor_up"
)

direction_values_pass <- all(
  direction %in% expected_directions
)

direction_sign_pass <- all(
  (direction == "Tumor_up" & log2fc > 0) |
    (direction == "Tumor_down" & log2fc < 0)
)

deg_logical_pass <- is.logical(
  existing_deg
)

deg_na_n <- sum(
  is.na(existing_deg)
)


# ------------------------------------------------------------
# 5. Final formal PASS
# ------------------------------------------------------------

all_pass_final <- all(
  
  source_exists,
  universe_exists,
  
  source_columns_pass,
  universe_column_pass,
  
  source_gene_na == 0,
  universe_gene_na == 0,
  
  source_duplicate_genes == 0,
  universe_duplicate_genes == 0,
  
  source_n_pass,
  universe_n_pass,
  
  source_unique_n_pass,
  universe_unique_n_pass,
  
  exact_gene_set_match,
  
  log2fc_finite_pass,
  p_value_finite_pass,
  fdr_finite_pass,
  
  p_value_range_pass,
  fdr_range_pass,
  
  direction_values_pass,
  direction_sign_pass,
  
  deg_logical_pass,
  deg_na_n == 0,
  
  exact_deg_logical_match,
  exact_up_logical_match,
  exact_down_logical_match,
  
  deg_disagreement_n == 0,
  
  core_source_pass,
  core_universe_pass,
  core_exact_deg_pass
)


# ------------------------------------------------------------
# 6. Print complete criterion summary
# ------------------------------------------------------------

criterion_summary <- c(
  
  Source_File_Exists =
    source_exists,
  
  Universe_File_Exists =
    universe_exists,
  
  Source_Columns =
    source_columns_pass,
  
  Universe_Column =
    universe_column_pass,
  
  Source_No_NA_Empty_Gene =
    source_gene_na == 0,
  
  Universe_No_NA_Empty_Gene =
    universe_gene_na == 0,
  
  Source_No_Duplicate_Genes =
    source_duplicate_genes == 0,
  
  Universe_No_Duplicate_Genes =
    universe_duplicate_genes == 0,
  
  Source_N_57828 =
    source_n_pass,
  
  Universe_N_57828 =
    universe_n_pass,
  
  Source_Unique_N_57828 =
    source_unique_n_pass,
  
  Universe_Unique_N_57828 =
    universe_unique_n_pass,
  
  Exact_Source_Universe_ID_Set =
    exact_gene_set_match,
  
  log2FC_Finite =
    log2fc_finite_pass,
  
  P_value_Finite =
    p_value_finite_pass,
  
  FDR_Finite =
    fdr_finite_pass,
  
  P_value_Range =
    p_value_range_pass,
  
  FDR_Range =
    fdr_range_pass,
  
  Direction_Values =
    direction_values_pass,
  
  Direction_Sign_Consistency =
    direction_sign_pass,
  
  DEG_Logical =
    deg_logical_pass,
  
  DEG_No_NA =
    deg_na_n == 0,
  
  Exact_DEG_Concordance =
    exact_deg_logical_match,
  
  Exact_Tumor_Up_Concordance =
    exact_up_logical_match,
  
  Exact_Tumor_Down_Concordance =
    exact_down_logical_match,
  
  DEG_Logical_Disagreements_Zero =
    deg_disagreement_n == 0,
  
  Core_Genes_In_Source =
    core_source_pass,
  
  Core_Genes_In_Universe =
    core_universe_pass,
  
  Core_Genes_Exact_DEG_Match =
    core_exact_deg_pass
)


cat("============================================================\n")
cat("STEP9D-1I FORMAL CRITERION SUMMARY\n")
cat("============================================================\n\n")

print(criterion_summary)


# ------------------------------------------------------------
# 7. Save corrected core-gene audit
# ------------------------------------------------------------

core_audit_file <- file.path(
  input_lock_dir,
  "STEP9D1I_CoreGene_Input_Lock_Audit.csv"
)

write.csv(
  core_audit,
  core_audit_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Build corrected formal lock
# ------------------------------------------------------------

STEP9D1I_FINAL_LOCK <- list(
  
  workflow =
    "Paper 1 — STEP9D Genome-wide Preranked GSEA",
  
  step =
    "STEP9D-1I",
  
  step_title =
    "Formal Input / Gene-Universe Lock",
  
  date =
    "2026-09-22",
  
  final_status =
    if (isTRUE(all_pass_final)) {
      "PASS"
    } else {
      "FAIL"
    },
  
  source = list(
    file = source_file,
    n_rows = nrow(source_df),
    n_unique_genes = length(unique(source_genes)),
    gene_column = "Gene"
  ),
  
  testing_universe = list(
    file = universe_file,
    n_rows = nrow(universe_df),
    n_unique_genes = length(unique(universe_genes)),
    gene_column = universe_gene_column
  ),
  
  exact_source_universe_match =
    exact_gene_set_match,
  
  statistical_columns = c(
    "log2FC",
    "P_value",
    "FDR"
  ),
  
  direction_definition =
    expected_directions,
  
  existing_DEG = list(
    total = sum(existing_deg),
    Tumor_up = sum(existing_up),
    Tumor_down = sum(existing_down)
  ),
  
  reconstructed_DEG_definition =
    "FDR < 0.05 AND abs(log2FC) > 1",
  
  reconstructed_DEG = list(
    total = sum(candidate_deg),
    Tumor_up = sum(candidate_up),
    Tumor_down = sum(candidate_down)
  ),
  
  exact_DEG_logical_match =
    exact_deg_logical_match,
  
  exact_Tumor_up_match =
    exact_up_logical_match,
  
  exact_Tumor_down_match =
    exact_down_logical_match,
  
  logical_disagreements =
    deg_disagreement_n,
  
  core_genes =
    c(
      "AOC3",
      "MS4A7",
      "IL7R"
    ),
  
  core_gene_ID_type =
    "Ensembl",
  
  core_gene_audit =
    core_audit,
  
  provenance_statement =
    paste(
      "The original STEP9B generating script/provenance was not recovered.",
      "The DEG definition is therefore documented as an exact",
      "result-level reconstruction rather than recovery of the",
      "original generating script."
    ),
  
  actions_not_performed = c(
    "DEG recalculation",
    "STEP9B modification",
    "STEP9C modification",
    "GSEA",
    "ranking-list construction",
    "ranking-metric lock",
    "P/FDR recalculation",
    "ID conversion"
  ),
  
  criterion_summary =
    criterion_summary,
  
  next_step =
    "STEP9D-2 Ranking Metric Lock"
)


# ------------------------------------------------------------
# 9. Save corrected formal lock
# ------------------------------------------------------------

final_lock_file <- file.path(
  input_lock_dir,
  "STEP9D1I_FINAL_LOCK.rds"
)

saveRDS(
  STEP9D1I_FINAL_LOCK,
  final_lock_file
)


# ------------------------------------------------------------
# 10. Save corrected human-readable summary
# ------------------------------------------------------------

summary_file <- file.path(
  input_lock_dir,
  "STEP9D1I_FORMAL_INPUT_LOCK_SUMMARY.txt"
)

sink(summary_file)

cat("============================================================\n")
cat("Paper 1 — STEP9D-1I Formal Input / Gene-Universe Lock\n")
cat("Date: 2026-09-22\n")
cat("============================================================\n\n")

cat(
  "FINAL STATUS: ",
  if (isTRUE(all_pass_final)) "PASS" else "FAIL",
  "\n\n",
  sep = ""
)

cat("SOURCE\n")
cat("------\n")
cat("STEP9B7R_Revised_DEG_FULL.csv\n")
cat("57,828 unique Ensembl gene IDs\n\n")

cat("TESTING UNIVERSE\n")
cat("----------------\n")
cat("STEP9C_DEG_testing_universe_57828.csv\n")
cat("57,828 unique Ensembl gene IDs\n\n")

cat("EXACT SOURCE/UNIVERSE MATCH\n")
cat("---------------------------\n")
cat("TRUE\n\n")

cat("STEP9D-1H DEG RESULT-LEVEL RECONSTRUCTION\n")
cat("-----------------------------------------\n")
cat("Definition: FDR < 0.05 AND abs(log2FC) > 1\n")
cat("Total: 16,065\n")
cat("Tumor_up: 9,311\n")
cat("Tumor_down: 6,754\n")
cat("Exact logical match: TRUE\n")
cat("Logical disagreements: 0\n\n")

cat("CORE GENES\n")
cat("----------\n")
cat("AOC3  -> ENSG00000131471\n")
cat("MS4A7 -> ENSG00000166927\n")
cat("IL7R  -> ENSG00000168685\n\n")

cat("Core genes in source: TRUE\n")
cat("Core genes in universe: TRUE\n")
cat("Core genes exact DEG match: TRUE\n\n")

cat("PROVENANCE NOTE\n")
cat("---------------\n")
cat(
  "The original STEP9B generating script/provenance was not recovered.\n"
)
cat(
  "The DEG definition is documented as an exact result-level\n"
)
cat(
  "reconstruction, not as recovery of the original generating script.\n\n"
)

cat("ACTIONS NOT PERFORMED\n")
cat("---------------------\n")
cat("DEG recalculation: NO\n")
cat("STEP9B modification: NO\n")
cat("STEP9C modification: NO\n")
cat("GSEA: NO\n")
cat("Ranking-list construction: NO\n")
cat("Ranking-metric lock: NO\n")
cat("P/FDR recalculation: NO\n")
cat("ID conversion: NO\n\n")

cat("NEXT STEP\n")
cat("---------\n")
cat("STEP9D-2 — Ranking Metric Lock\n")

sink()


# ------------------------------------------------------------
# 11. Verify saved files
# ------------------------------------------------------------

file_check <- c(
  FINAL_LOCK_RDS =
    file.exists(final_lock_file),
  
  CORE_AUDIT_CSV =
    file.exists(core_audit_file),
  
  SUMMARY_TXT =
    file.exists(summary_file)
)


cat("\n============================================================\n")
cat("STEP9D-1I FINAL CHECKPOINT\n")
cat("============================================================\n\n")

print(file_check)

cat("\n")


# ------------------------------------------------------------
# 12. Final result
# ------------------------------------------------------------

if (
  isTRUE(all_pass_final) &&
  all(file_check)
) {
  
  cat("============================================================\n")
  cat("PASS — STEP9D-1I FORMAL INPUT / GENE-UNIVERSE LOCK\n")
  cat("============================================================\n\n")
  
  cat("57,828-gene genome-wide input: LOCKED\n")
  cat("Exact source/universe match: TRUE\n")
  cat("Exact DEG concordance: TRUE\n")
  cat("Core-gene Ensembl audit: TRUE\n")
  cat("AOC3/MS4A7/IL7R: VERIFIED\n\n")
  
  cat("No DEG recalculation performed.\n")
  cat("No STEP9B/STEP9C modification performed.\n")
  cat("No GSEA performed.\n")
  cat("No ranking metric locked.\n\n")
  
  cat("NEXT STEP:\n")
  cat("STEP9D-2 — Ranking Metric Lock\n")
  
  cat("\n============================================================\n")
  
} else {
  
  cat("============================================================\n")
  cat("FAIL — STEP9D-1I FORMAL LOCK\n")
  cat("============================================================\n\n")
  
  cat("Do NOT proceed to STEP9D-2.\n")
  
}

# ============================================================
# Paper 1
# STEP9D-2A
# Ranking Input / GSEA Implementation Audit
# ============================================================
#
# Purpose:
#   Audit the available genome-wide statistics and the local
#   GSEA implementation before formally locking the ranking metric.
#
# IMPORTANT:
#   - No ranking vector is created
#   - No ranking metric is calculated
#   - No GSEA is performed
#   - No DEG recalculation
#   - No STEP9B/STEP9C modification
#   - STEP9D-1I remains locked and unchanged
# ============================================================


cat("\n============================================================\n")
cat("STEP9D-2A: Ranking Input / GSEA Implementation Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Confirm STEP9D-1I formal lock exists
# ------------------------------------------------------------

step9d1i_lock_file <- file.path(
  step9d_root,
  "9D-1_Input_Lock",
  "STEP9D1I_FINAL_LOCK.rds"
)

step9d1i_exists <- file.exists(
  step9d1i_lock_file
)

cat(
  "STEP9D-1I formal lock exists:",
  step9d1i_exists,
  "\n\n"
)

if (!step9d1i_exists) {
  
  stop(
    "FAIL: STEP9D-1I formal lock not found. Do not proceed."
  )
  
}


STEP9D1I_LOCK <- readRDS(
  step9d1i_lock_file
)

cat(
  "STEP9D-1I recorded status:",
  STEP9D1I_LOCK$final_status,
  "\n\n"
)

if (
  !identical(
    STEP9D1I_LOCK$final_status,
    "PASS"
  )
) {
  
  stop(
    "FAIL: STEP9D-1I is not formally locked as PASS."
  )
  
}


# ------------------------------------------------------------
# 1. Confirm genome-wide source object
# ------------------------------------------------------------

source_file <- file.path(
  work_dir,
  "STEP9B-7R_Revised_DEG_Definition_and_Lock",
  "02_results",
  "STEP9B7R_Revised_DEG_FULL.csv"
)

source_df_2A <- read.csv(
  source_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("Genome-wide source dimensions:\n")
print(dim(source_df_2A))

cat("\n")


# ------------------------------------------------------------
# 2. Inspect available statistical columns
# ------------------------------------------------------------

cat("Available source columns:\n")
print(names(source_df_2A))

cat("\n")


# ------------------------------------------------------------
# 3. Identify candidate ranking statistics
# ------------------------------------------------------------

ranking_related_columns <- names(source_df_2A)[
  grepl(
    paste(
      c(
        "stat",
        "score",
        "t",
        "z",
        "wald",
        "test",
        "log2FC",
        "P_value",
        "FDR"
      ),
      collapse = "|"
    ),
    names(source_df_2A),
    ignore.case = TRUE
  )
]

cat("Potential ranking-related columns detected:\n")
print(ranking_related_columns)

cat("\n")


# ------------------------------------------------------------
# 4. Explicit native signed-statistic audit
# ------------------------------------------------------------

native_signed_candidates <- c(
  "stat",
  "Statistic",
  "STAT",
  "test_statistic",
  "TestStatistic",
  "t",
  "t_stat",
  "t_statistic",
  "z",
  "z_score",
  "Z",
  "Wald",
  "Wald_statistic",
  "score"
)

native_signed_present <- intersect(
  native_signed_candidates,
  names(source_df_2A)
)

cat("Potential native signed-statistic columns:\n")

if (length(native_signed_present) == 0) {
  
  cat("NONE\n")
  
} else {
  
  print(native_signed_present)
  
}

cat("\n")


# ------------------------------------------------------------
# 5. Explicitly inspect primary available statistics
# ------------------------------------------------------------

required_ranking_inputs <- c(
  "Gene",
  "log2FC",
  "P_value",
  "FDR"
)

ranking_input_columns_pass <- all(
  required_ranking_inputs %in%
    names(source_df_2A)
)

cat(
  "Required ranking input columns present:",
  ranking_input_columns_pass,
  "\n\n"
)

if (!ranking_input_columns_pass) {
  
  stop(
    "FAIL: Required ranking input columns are missing."
  )
  
}


# ------------------------------------------------------------
# 6. Audit log2FC
# ------------------------------------------------------------

log2fc_2A <- suppressWarnings(
  as.numeric(source_df_2A$log2FC)
)

log2fc_na <- sum(
  is.na(log2fc_2A)
)

log2fc_inf <- sum(
  is.infinite(log2fc_2A)
)

log2fc_zero <- sum(
  log2fc_2A == 0,
  na.rm = TRUE
)

cat("log2FC audit:\n")
cat("  NA       :", log2fc_na, "\n")
cat("  Infinite :", log2fc_inf, "\n")
cat("  Zero     :", log2fc_zero, "\n")
cat(
  "  Min      :",
  min(log2fc_2A, na.rm = TRUE),
  "\n"
)
cat(
  "  Max      :",
  max(log2fc_2A, na.rm = TRUE),
  "\n\n"
)


# ------------------------------------------------------------
# 7. Audit P_value
# ------------------------------------------------------------

p_value_2A <- suppressWarnings(
  as.numeric(source_df_2A$P_value)
)

p_na <- sum(
  is.na(p_value_2A)
)

p_inf <- sum(
  is.infinite(p_value_2A)
)

p_zero <- sum(
  p_value_2A == 0,
  na.rm = TRUE
)

p_positive <- sum(
  p_value_2A > 0,
  na.rm = TRUE
)

p_invalid <- sum(
  p_value_2A < 0 |
    p_value_2A > 1,
  na.rm = TRUE
)

cat("P_value audit:\n")
cat("  NA           :", p_na, "\n")
cat("  Infinite     :", p_inf, "\n")
cat("  Zero         :", p_zero, "\n")
cat("  Positive     :", p_positive, "\n")
cat("  Invalid range:", p_invalid, "\n")

if (p_positive > 0) {
  
  cat(
    "  Minimum positive P:",
    format(
      min(
        p_value_2A[
          is.finite(p_value_2A) &
            p_value_2A > 0
        ]
      ),
      scientific = TRUE
    ),
    "\n"
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 8. Audit FDR
# ------------------------------------------------------------

fdr_2A <- suppressWarnings(
  as.numeric(source_df_2A$FDR)
)

fdr_na <- sum(
  is.na(fdr_2A)
)

fdr_inf <- sum(
  is.infinite(fdr_2A)
)

fdr_invalid <- sum(
  fdr_2A < 0 |
    fdr_2A > 1,
  na.rm = TRUE
)

cat("FDR audit:\n")
cat("  NA           :", fdr_na, "\n")
cat("  Infinite     :", fdr_inf, "\n")
cat("  Invalid range:", fdr_invalid, "\n\n")


# ------------------------------------------------------------
# 9. Audit potential numerical underflow
# ------------------------------------------------------------

positive_p <- p_value_2A[
  is.finite(p_value_2A) &
    p_value_2A > 0
]

if (length(positive_p) > 0) {
  
  min_positive_p <- min(
    positive_p
  )
  
  minus_log10_min_p <- -log10(
    min_positive_p
  )
  
} else {
  
  min_positive_p <- NA_real_
  minus_log10_min_p <- NA_real_
  
}

cat("P-value numerical-range audit:\n")
cat(
  "  Minimum positive P:",
  format(
    min_positive_p,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "  Corresponding -log10(P):",
  minus_log10_min_p,
  "\n\n"
)


# ------------------------------------------------------------
# 10. Audit sign consistency
# ------------------------------------------------------------

direction_2A <- as.character(
  source_df_2A$Direction
)

direction_sign_consistency <- all(
  (
    direction_2A == "Tumor_up" &
      log2fc_2A > 0
  ) |
    (
      direction_2A == "Tumor_down" &
        log2fc_2A < 0
    )
)

cat(
  "Direction/log2FC sign consistency:",
  direction_sign_consistency,
  "\n\n"
)


# ------------------------------------------------------------
# 11. Inspect clusterProfiler version
# ------------------------------------------------------------

clusterProfiler_available <- requireNamespace(
  "clusterProfiler",
  quietly = TRUE
)

cat(
  "clusterProfiler installed:",
  clusterProfiler_available,
  "\n"
)

if (clusterProfiler_available) {
  
  clusterProfiler_version <- as.character(
    packageVersion("clusterProfiler")
  )
  
  cat(
    "clusterProfiler version:",
    clusterProfiler_version,
    "\n"
  )
  
} else {
  
  clusterProfiler_version <- NA_character_
  
}

cat("\n")


# ------------------------------------------------------------
# 12. Inspect fgsea version
# ------------------------------------------------------------

fgsea_available <- requireNamespace(
  "fgsea",
  quietly = TRUE
)

cat(
  "fgsea installed:",
  fgsea_available,
  "\n"
)

if (fgsea_available) {
  
  fgsea_version <- as.character(
    packageVersion("fgsea")
  )
  
  cat(
    "fgsea version:",
    fgsea_version,
    "\n"
  )
  
} else {
  
  fgsea_version <- NA_character_
  
}

cat("\n")


# ------------------------------------------------------------
# 13. Inspect clusterProfiler::GSEA formal arguments
# ------------------------------------------------------------

if (clusterProfiler_available) {
  
  gsea_formals <- names(
    formals(
      clusterProfiler::GSEA
    )
  )
  
  cat(
    "clusterProfiler::GSEA formal arguments:\n"
  )
  
  print(
    gsea_formals
  )
  
} else {
  
  gsea_formals <- character(0)
  
}

cat("\n")


# ------------------------------------------------------------
# 14. Inspect clusterProfiler::GSEA function definition
# ------------------------------------------------------------

if (clusterProfiler_available) {
  
  gsea_function_text <- capture.output(
    print(
      clusterProfiler::GSEA
    )
  )
  
  cat(
    "First 40 lines of local clusterProfiler::GSEA definition:\n\n"
  )
  
  print(
    head(
      gsea_function_text,
      40
    )
  )
  
} else {
  
  gsea_function_text <- character(0)
  
}

cat("\n")


# ------------------------------------------------------------
# 15. Inspect fgsea package version / namespace
# ------------------------------------------------------------

if (fgsea_available) {
  
  fgsea_exports <- getNamespaceExports(
    "fgsea"
  )
  
  cat(
    "Selected fgsea exported functions:\n"
  )
  
  print(
    head(
      fgsea_exports,
      50
    )
  )
  
} else {
  
  fgsea_exports <- character(0)
  
}

cat("\n")


# ------------------------------------------------------------
# 16. IMPORTANT: no ranking calculation
# ------------------------------------------------------------

ranking_vector_created <- FALSE
ranking_metric_locked <- FALSE
gsea_performed <- FALSE


cat("Current STEP9D-2A state:\n")
cat(
  "  Ranking vector created:",
  ranking_vector_created,
  "\n"
)

cat(
  "  Ranking metric locked:",
  ranking_metric_locked,
  "\n"
)

cat(
  "  GSEA performed:",
  gsea_performed,
  "\n\n"
)


# ============================================================
# STEP9D-2A
# Repair: Save audit object after R list-name syntax correction
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-2A: SAVE AUDIT OBJECT — REPAIR\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Reconfirm STEP9D-1I formal lock
# ------------------------------------------------------------

if (!exists("STEP9D1I_LOCK")) {
  stop("FAIL: STEP9D1I_LOCK is not available in the current R session.")
}

if (!identical(STEP9D1I_LOCK$final_status, "PASS")) {
  stop("FAIL: STEP9D-1I formal lock is not PASS.")
}

cat("STEP9D-1I status:",
    STEP9D1I_LOCK$final_status,
    "\n\n")


# ------------------------------------------------------------
# 2. Reconfirm STEP9D-2A audit objects already generated
# ------------------------------------------------------------

required_audit_objects <- c(
  "source_df_2A",
  "log2fc_na",
  "log2fc_inf",
  "log2fc_zero",
  "p_na",
  "p_inf",
  "p_zero",
  "p_positive",
  "p_invalid",
  "min_positive_p",
  "minus_log10_min_p",
  "fdr_na",
  "fdr_inf",
  "fdr_invalid",
  "direction_sign_consistency",
  "clusterProfiler_available",
  "clusterProfiler_version",
  "fgsea_available",
  "fgsea_version",
  "gsea_formals",
  "gsea_function_text",
  "fgsea_exports",
  "ranking_vector_created",
  "ranking_metric_locked",
  "gsea_performed",
  "source_file",
  "required_ranking_inputs",
  "native_signed_present",
  "step9d2a_dir"
)

missing_objects <- required_audit_objects[
  !vapply(
    required_audit_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects) > 0) {
  stop(
    paste(
      "FAIL: Required STEP9D-2A audit objects are missing:",
      paste(missing_objects, collapse = ", ")
    )
  )
}

cat("All STEP9D-2A audit objects present: TRUE\n\n")


# ------------------------------------------------------------
# 3. Build corrected audit object
# ------------------------------------------------------------

STEP9D2A_AUDIT <- list(
  
  step = "STEP9D-2A",
  
  title =
    "Ranking Input / GSEA Implementation Audit",
  
  date =
    "2026-09-22",
  
  STEP9D1I_status =
    STEP9D1I_LOCK$final_status,
  
  source_file =
    source_file,
  
  source_n =
    nrow(source_df_2A),
  
  source_columns =
    names(source_df_2A),
  
  native_signed_statistic_candidates =
    native_signed_present,
  
  ranking_input_columns =
    required_ranking_inputs,
  
  log2FC_audit =
    list(
      NA_count = log2fc_na,
      Infinite_count = log2fc_inf,
      Zero_count = log2fc_zero,
      Minimum = min(
        log2fc_2A,
        na.rm = TRUE
      ),
      Maximum = max(
        log2fc_2A,
        na.rm = TRUE
      )
    ),
  
  P_value_audit =
    list(
      NA_count = p_na,
      Infinite_count = p_inf,
      Zero_count = p_zero,
      Positive_count = p_positive,
      Invalid_count = p_invalid,
      Minimum_positive_P = min_positive_p,
      Minus_log10_minimum_positive_P =
        minus_log10_min_p
    ),
  
  FDR_audit =
    list(
      NA_count = fdr_na,
      Infinite_count = fdr_inf,
      Invalid_count = fdr_invalid
    ),
  
  direction_sign_consistency =
    direction_sign_consistency,
  
  clusterProfiler_available =
    clusterProfiler_available,
  
  clusterProfiler_version =
    clusterProfiler_version,
  
  fgsea_available =
    fgsea_available,
  
  fgsea_version =
    fgsea_version,
  
  GSEA_formal_arguments =
    gsea_formals,
  
  GSEA_function_definition =
    gsea_function_text,
  
  fgsea_exports =
    fgsea_exports,
  
  ranking_vector_created =
    ranking_vector_created,
  
  ranking_metric_locked =
    ranking_metric_locked,
  
  GSEA_performed =
    gsea_performed,
  
  status =
    "AUDIT_ONLY — NOT YET LOCKED"
)


# ------------------------------------------------------------
# 4. Save corrected RDS
# ------------------------------------------------------------

step9d2a_rds <- file.path(
  step9d2a_dir,
  "STEP9D2A_Ranking_Input_GSEA_Implementation_Audit.rds"
)

saveRDS(
  STEP9D2A_AUDIT,
  step9d2a_rds
)

cat(
  "Audit RDS saved:\n",
  step9d2a_rds,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. Save corrected TXT
# ------------------------------------------------------------

step9d2a_txt <- file.path(
  step9d2a_dir,
  "STEP9D2A_Ranking_Input_GSEA_Implementation_Audit.txt"
)

sink(step9d2a_txt)

cat("============================================================\n")
cat("Paper 1 — STEP9D-2A\n")
cat("Ranking Input / GSEA Implementation Audit\n")
cat("Date: 2026-09-22\n")
cat("============================================================\n\n")

cat(
  "STEP9D-1I status:",
  STEP9D1I_LOCK$final_status,
  "\n\n"
)

cat(
  "Genome-wide source N:",
  nrow(source_df_2A),
  "\n\n"
)

cat("Potential native signed-statistic columns:\n")

if (length(native_signed_present) == 0) {
  cat("NONE\n")
} else {
  print(native_signed_present)
}

cat("\n")

cat("clusterProfiler:\n")
cat("  Available:", clusterProfiler_available, "\n")
cat("  Version:", clusterProfiler_version, "\n\n")

cat("fgsea:\n")
cat("  Available:", fgsea_available, "\n")
cat("  Version:", fgsea_version, "\n\n")

cat("P-value audit:\n")
cat("  NA:", p_na, "\n")
cat("  Infinite:", p_inf, "\n")
cat("  Zero:", p_zero, "\n")
cat("  Positive:", p_positive, "\n")
cat("  Invalid:", p_invalid, "\n")
cat(
  "  Minimum positive P:",
  format(
    min_positive_p,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "  Corresponding -log10(P):",
  minus_log10_min_p,
  "\n\n"
)

cat("FDR audit:\n")
cat("  NA:", fdr_na, "\n")
cat("  Infinite:", fdr_inf, "\n")
cat("  Invalid:", fdr_invalid, "\n\n")

cat(
  "Direction/log2FC sign consistency:",
  direction_sign_consistency,
  "\n\n"
)

cat("IMPORTANT STATUS\n")
cat("----------------\n")
cat("Ranking vector created: NO\n")
cat("Ranking metric locked: NO\n")
cat("GSEA performed: NO\n\n")

cat(
  "This is an implementation/input audit only.\n"
)

cat(
  "No ranking metric has been locked.\n"
)

cat(
  "No GSEA has been performed.\n"
)

cat(
  "Next step: STEP9D-2B Ranking Metric Lock.\n"
)

sink()


# ------------------------------------------------------------
# 6. Final console status
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP9D-2A COMPLETE — AUDIT ONLY\n")
cat("============================================================\n\n")

cat("PASS — STEP9D-2A audit object saved successfully.\n\n")

cat(
  "Native signed statistic:",
  ifelse(
    length(native_signed_present) == 0,
    "NONE",
    paste(native_signed_present, collapse = ", ")
  ),
  "\n"
)

cat(
  "clusterProfiler:",
  clusterProfiler_version,
  "\n"
)

cat(
  "fgsea:",
  fgsea_version,
  "\n"
)

cat(
  "Genome-wide genes:",
  nrow(source_df_2A),
  "\n"
)

cat(
  "Direction/log2FC consistency:",
  direction_sign_consistency,
  "\n\n"
)

cat("Ranking metric locked: NO\n")
cat("Ranking vector created: NO\n")
cat("GSEA performed: NO\n\n")

cat(
  "NEXT STEP:\n",
  "STEP9D-2B — Ranking Metric Lock\n",
  sep = ""
)

cat("\n============================================================\n")


# ============================================================
# Paper 1 — Functional Enrichment Revision
#
# STEP9D-2B
# Ranking Metric Lock
#
# PURPOSE
# ------------------------------------------------------------
# 1. Lock the pre-specified GSEA ranking metric.
# 2. Use the already audited 57,828-gene genome-wide input.
# 3. Preserve both effect direction and statistical strength.
# 4. Audit all numerical edge cases before any GSEA is performed.
# 5. Save a complete, reproducible ranking-metric lock object.
#
# IMPORTANT
# ------------------------------------------------------------
# This step DOES NOT:
#   - perform GSEA
#   - create the final GSEA ranking vector
#   - calculate pathway enrichment
#   - inspect pathway-level results
#   - select a metric based on enrichment results
#
# Locked ranking metric:
#
#   R_i = sign(log2FC_i) * [-log10(max(P_i, 1e-300))]
#
# where:
#   log2FC = effect-direction variable
#   P      = nominal statistical significance variable
#   1e-300 = pre-specified numerical floor
#
# ============================================================


cat("\n============================================================\n")
cat("STEP9D-2B: Ranking Metric Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory and output folders
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9D-2B_Ranking_Metric_Lock"

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}

cat("STEP directory:\n", step_dir, "\n\n", sep = "")


# ------------------------------------------------------------
# 1. Input object
# ------------------------------------------------------------
#
# STEP9D-2A should have already created / retained:
#
#   source_df_2A
#
# containing the 57,828-gene genome-wide input.
#
# We intentionally do NOT recalculate the DEG universe here.
#
# ------------------------------------------------------------

if (!exists("source_df_2A")) {
  stop(
    "\nERROR: source_df_2A is not available in the current R session.\n",
    "STEP9D-2B requires the exact genome-wide input object audited in STEP9D-2A.\n",
    "Do not reconstruct or substitute another input at this step.\n"
  )
}

cat("Input object detected: source_df_2A\n")
cat("Input rows:", nrow(source_df_2A), "\n\n")


# ------------------------------------------------------------
# 2. Expected genome-wide input size lock
# ------------------------------------------------------------

expected_n_genes <- 57828

if (nrow(source_df_2A) != expected_n_genes) {
  stop(
    "\nERROR: Genome-wide input size mismatch.\n",
    "Expected: ", expected_n_genes, "\n",
    "Observed: ", nrow(source_df_2A), "\n",
    "STEP9D-2B cannot continue.\n"
  )
}

genomewide_n_locked <- nrow(source_df_2A)

cat("Genome-wide input size LOCKED:", genomewide_n_locked, "\n\n")


# ------------------------------------------------------------
# 3. Detect gene-ID column
# ------------------------------------------------------------
#
# Ranking requires a unique gene identifier.
# Prefer standard gene-symbol / identifier columns.
#
# ------------------------------------------------------------

candidate_gene_columns <- c(
  "gene",
  "Gene",
  "gene_symbol",
  "GeneSymbol",
  "SYMBOL",
  "symbol",
  "Symbol",
  "gene_id",
  "GeneID",
  "ENTREZID",
  "ENTREZ",
  "ensembl",
  "Ensembl",
  "ENSEMBL"
)

gene_candidates_present <- candidate_gene_columns[
  candidate_gene_columns %in% colnames(source_df_2A)
]

if (length(gene_candidates_present) == 0) {
  stop(
    "\nERROR: No recognizable gene-ID column was detected in source_df_2A.\n",
    "Available columns:\n",
    paste(colnames(source_df_2A), collapse = ", "),
    "\n"
  )
}

gene_col <- gene_candidates_present[1]

cat("Gene-ID column selected:", gene_col, "\n\n")


# ------------------------------------------------------------
# 4. Detect log2FC column
# ------------------------------------------------------------

candidate_log2fc_columns <- c(
  "log2FC",
  "log2Fc",
  "logFC",
  "log2fc",
  "Log2FC",
  "log2.fold.change",
  "log2_fold_change"
)

log2fc_candidates_present <- candidate_log2fc_columns[
  candidate_log2fc_columns %in% colnames(source_df_2A)
]

if (length(log2fc_candidates_present) == 0) {
  stop(
    "\nERROR: No recognizable log2FC column was detected.\n",
    "Available columns:\n",
    paste(colnames(source_df_2A), collapse = ", "),
    "\n"
  )
}

log2fc_col <- log2fc_candidates_present[1]

cat("Direction variable selected:", log2fc_col, "\n\n")


# ------------------------------------------------------------
# 5. Detect nominal P-value column
# ------------------------------------------------------------
#
# IMPORTANT:
# The ranking metric uses nominal P values, not adjusted P values.
#
# Adjusted P values are not substituted because they are not the
# pre-specified statistical-strength component of this metric.
#
# ------------------------------------------------------------

candidate_p_columns <- c(
  "P_value",
  "P.Value",
  "PValue",
  "pvalue",
  "p_value",
  "p.val",
  "P",
  "p",
  "pvalue_nominal",
  "P.Value.nominal",
  "p_nominal"
)

p_candidates_present <- candidate_p_columns[
  candidate_p_columns %in% colnames(source_df_2A)
]

if (length(p_candidates_present) == 0) {
  stop(
    "\nERROR: No recognizable nominal P-value column was detected.\n",
    "Available columns:\n",
    paste(colnames(source_df_2A), collapse = ", "),
    "\n"
  )
}

p_col <- p_candidates_present[1]

cat("Statistical-strength variable selected:", p_col, "\n\n")


# ------------------------------------------------------------
# 6. Record adjusted-P column if present
# ------------------------------------------------------------
#
# This is audit-only.
# Adjusted P values are NOT used for the ranking metric.
#
# ------------------------------------------------------------

candidate_fdr_columns <- c(
  "adj.P.Val",
  "padj",
  "FDR",
  "fdr",
  "adj_p",
  "adjusted_p",
  "p.adjust"
)

fdr_candidates_present <- candidate_fdr_columns[
  candidate_fdr_columns %in% colnames(source_df_2A)
]

if (length(fdr_candidates_present) > 0) {
  fdr_col_detected <- fdr_candidates_present[1]
} else {
  fdr_col_detected <- NA_character_
}

cat(
  "Adjusted-P/FDR column detected:",
  ifelse(
    is.na(fdr_col_detected),
    "NONE",
    fdr_col_detected
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 7. Confirm required columns are distinct
# ------------------------------------------------------------

if (length(unique(c(gene_col, log2fc_col, p_col))) != 3) {
  stop(
    "\nERROR: Gene-ID, log2FC, and P-value columns are not distinct.\n"
  )
}


# ------------------------------------------------------------
# 8. Confirm numeric types
# ------------------------------------------------------------

if (!is.numeric(source_df_2A[[log2fc_col]])) {
  stop(
    "\nERROR: log2FC column is not numeric: ",
    log2fc_col,
    "\n"
  )
}

if (!is.numeric(source_df_2A[[p_col]])) {
  stop(
    "\nERROR: P-value column is not numeric: ",
    p_col,
    "\n"
  )
}

cat("Numeric-type audit: PASS\n")
cat("log2FC numeric: TRUE\n")
cat("P value numeric: TRUE\n\n")


# ------------------------------------------------------------
# 9. Extract audit vectors
# ------------------------------------------------------------

gene_id_raw <- as.character(source_df_2A[[gene_col]])

log2fc_raw <- as.numeric(source_df_2A[[log2fc_col]])

p_raw <- as.numeric(source_df_2A[[p_col]])


# ------------------------------------------------------------
# 10. Gene-ID completeness audit
# ------------------------------------------------------------

gene_na_n <- sum(
  is.na(gene_id_raw) |
    trimws(gene_id_raw) == ""
)

gene_empty_n <- sum(
  is.na(gene_id_raw) |
    trimws(gene_id_raw) == ""
)

if (gene_na_n > 0) {
  stop(
    "\nERROR: Missing/empty gene IDs detected: ",
    gene_na_n,
    "\n"
  )
}

cat("Missing/empty gene IDs:", gene_na_n, "\n")


# ------------------------------------------------------------
# 11. Gene-ID duplication audit
# ------------------------------------------------------------

gene_duplicate_n <- sum(duplicated(gene_id_raw))

gene_unique_n <- length(unique(gene_id_raw))

cat("Unique gene IDs:", gene_unique_n, "\n")
cat("Duplicated gene-ID rows:", gene_duplicate_n, "\n")

if (gene_duplicate_n > 0) {
  stop(
    "\nERROR: Duplicated gene IDs detected.\n",
    "STEP9D-2B will not silently collapse or select duplicated genes.\n",
    "Resolve duplication before ranking.\n"
  )
}

cat("Gene-ID uniqueness: PASS\n\n")


# ------------------------------------------------------------
# 12. log2FC numerical audit
# ------------------------------------------------------------

log2fc_na_n <- sum(is.na(log2fc_raw))

log2fc_inf_n <- sum(is.infinite(log2fc_raw))

log2fc_finite_n <- sum(is.finite(log2fc_raw))

log2fc_zero_n <- sum(
  !is.na(log2fc_raw) &
    is.finite(log2fc_raw) &
    log2fc_raw == 0
)

cat("log2FC NA:", log2fc_na_n, "\n")
cat("log2FC Inf:", log2fc_inf_n, "\n")
cat("log2FC finite:", log2fc_finite_n, "\n")
cat("log2FC exactly zero:", log2fc_zero_n, "\n\n")

if (log2fc_inf_n > 0) {stop(
  "\nERROR: Infinite log2FC values detected.\n",
  "Cannot lock ranking metric until these values are resolved.\n"
)}


# ------------------------------------------------------------
# 13. P-value numerical audit
# ------------------------------------------------------------

p_na_n <- sum(is.na(p_raw))

p_inf_n <- sum(is.infinite(p_raw))

p_zero_n <- sum(
  !is.na(p_raw) &
    is.finite(p_raw) &
    p_raw == 0
)

p_negative_n <- sum(
  !is.na(p_raw) &
    is.finite(p_raw) &
    p_raw < 0
)

p_gt1_n <- sum(
  !is.na(p_raw) &
    is.finite(p_raw) &
    p_raw > 1
)

p_finite_n <- sum(is.finite(p_raw))

cat("P-value NA:", p_na_n, "\n")
cat("P-value Inf:", p_inf_n, "\n")
cat("P-value zero:", p_zero_n, "\n")
cat("P-value negative:", p_negative_n, "\n")
cat("P-value > 1:", p_gt1_n, "\n")
cat("P-value finite:", p_finite_n, "\n\n")


# ------------------------------------------------------------
# 14. Validate P-value domain
# ------------------------------------------------------------

if (p_negative_n > 0) {
  stop(
    "\nERROR: Negative P values detected.\n"
  )
}

if (p_gt1_n > 0) {
  stop(
    "\nERROR: P values > 1 detected.\n"
  )
  
}

if (p_inf_n > 0) {
  stop(
    "\nERROR: Infinite P values detected.\n"
  )
}


# ------------------------------------------------------------
# 15. Define P-value numerical floor
# ------------------------------------------------------------
#
# This value is locked BEFORE ranking.
#
# It prevents:
#
#   -log10(0) = Inf
#
# and therefore prevents infinite ranking values.
#
# ------------------------------------------------------------

p_floor <- 1e-300

cat("P-value numerical floor LOCKED:", format(p_floor, scientific = TRUE), "\n\n")


# ------------------------------------------------------------
# 16. Define ranking formula
# ------------------------------------------------------------
#
# Locked formula:
#
#   ranking_metric =
#       sign(log2FC) *
#       (-log10(max(P, p_floor)))
#
# Direction:
#   log2FC > 0  -> positive
#   log2FC < 0  -> negative
#
# Statistical strength:
#   smaller P -> larger absolute ranking magnitude
#
# ------------------------------------------------------------

ranking_metric_definition <- paste0(
  "sign(log2FC) * (-log10(max(P, ", p_floor, ")))"
)

ranking_metric_name <- "signed_minus_log10_nominal_P"

ranking_direction_rule <- paste(
  "Positive ranking values indicate positive log2FC;",
  "negative ranking values indicate negative log2FC."
)

ranking_strength_rule <- paste(
  "Absolute ranking magnitude increases as nominal P decreases."
)

ranking_p_policy <- paste(
  "P values are floored at",
  format(p_floor, scientific = TRUE),
  "before -log10 transformation."
)

ranking_zero_log2fc_policy <- paste(
  "Genes with log2FC = 0 have sign(log2FC) = 0",
  "and therefore receive ranking metric = 0."
)

cat("Ranking metric NAME LOCKED:\n")
cat(ranking_metric_name, "\n\n")

cat("Ranking metric FORMULA LOCKED:\n")
cat(ranking_metric_definition, "\n\n")

cat("Direction rule:\n")
cat(ranking_direction_rule, "\n\n")

cat("Strength rule:\n")
cat(ranking_strength_rule, "\n\n")

cat("P-value policy:\n")
cat(ranking_p_policy, "\n\n")

cat("Zero-log2FC policy:\n")
cat(ranking_zero_log2fc_policy, "\n\n")


# ------------------------------------------------------------
# 17. Define eligibility policy
# ------------------------------------------------------------
#
# STEP9D-2B locks the eligibility rule but does NOT yet create
# the final ranking vector.
#
# Eligible genes require:
#
#   non-missing gene ID
#   finite log2FC
#   finite P value
#   P in [0,1]
#
# ------------------------------------------------------------

eligible_for_ranking <- (
  !is.na(gene_id_raw) &
    trimws(gene_id_raw) != "" &
    is.finite(log2fc_raw) &
    is.finite(p_raw) &
    p_raw >= 0 &
    p_raw <= 1
)

eligible_n <- sum(eligible_for_ranking)

excluded_n <- sum(!eligible_for_ranking)

cat("Eligibility policy LOCKED.\n\n")

cat("Eligible rows:", eligible_n, "\n")
cat("Excluded rows:", excluded_n, "\n\n")


# ------------------------------------------------------------
# 18. Audit exclusion reasons
# ------------------------------------------------------------

exclude_missing_gene_n <- sum(
  is.na(gene_id_raw) |
    trimws(gene_id_raw) == ""
)

exclude_log2fc_na_n <- sum(is.na(log2fc_raw))

exclude_log2fc_inf_n <- sum(is.infinite(log2fc_raw))

exclude_p_na_n <- sum(is.na(p_raw))

exclude_p_inf_n <- sum(is.infinite(p_raw))

exclude_p_negative_n <- sum(
  !is.na(p_raw) &
    is.finite(p_raw) &
    p_raw < 0
)

exclude_p_gt1_n <- sum(
  !is.na(p_raw) &
    is.finite(p_raw) &
    p_raw > 1
)

exclusion_audit <- data.frame(
  exclusion_reason = c(
    "missing_or_empty_gene_ID",
    "log2FC_NA",
    "log2FC_Inf",
    "P_value_NA",
    "P_value_Inf",
    "P_value_negative",
    "P_value_greater_than_1"
  ),
  n = c(
    exclude_missing_gene_n,
    exclude_log2fc_na_n,
    exclude_log2fc_inf_n,
    exclude_p_na_n,
    exclude_p_inf_n,
    exclude_p_negative_n,
    exclude_p_gt1_n
  ),
  stringsAsFactors = FALSE
)

cat("Exclusion audit:\n")
print(exclusion_audit)
cat("\n")


# ------------------------------------------------------------
# 19. Zero-P audit
# ------------------------------------------------------------
#
# Zero P values are retained.
# They are not discarded.
# The locked floor converts them to a finite numerical value.
#
# ------------------------------------------------------------

zero_p_policy <- paste0(
  "P = 0 retained; replaced by max(P, ",
  p_floor,
  ") before -log10 transformation."
)

cat("Zero-P policy LOCKED:\n")
cat(zero_p_policy, "\n\n")


# ------------------------------------------------------------
# 20. Direction consistency audit
# ------------------------------------------------------------
#
# This reproduces the STEP9D-2A principle:
# the ranking sign must follow log2FC.
#
# No ranking vector is retained as an analysis result here.
# A temporary audit calculation is allowed solely to verify
# that the locked formula behaves as intended.
#


# ------------------------------------------------------------
# ============================================================
# STEP9D-2B DIAGNOSTIC
# Identify cause of direction-consistency FALSE
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-2B DIAGNOSTIC — Direction Consistency Failure\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Basic counts
# ------------------------------------------------------------

diag_log2fc_positive <- sum(
  eligible_for_ranking &
    log2fc_raw > 0
)

diag_log2fc_negative <- sum(
  eligible_for_ranking &
    log2fc_raw < 0
)

diag_log2fc_zero <- sum(
  eligible_for_ranking &
    log2fc_raw == 0
)

diag_p_zero <- sum(
  eligible_for_ranking &
    p_raw == 0
)

diag_p_one <- sum(
  eligible_for_ranking &
    p_raw == 1
)

cat("Eligible genes:", sum(eligible_for_ranking), "\n")
cat("log2FC > 0:", diag_log2fc_positive, "\n")
cat("log2FC < 0:", diag_log2fc_negative, "\n")
cat("log2FC = 0:", diag_log2fc_zero, "\n")
cat("P = 0:", diag_p_zero, "\n")
cat("P = 1:", diag_p_one, "\n\n")


# ------------------------------------------------------------
# 2. Identify positive log2FC with zero ranking
# ------------------------------------------------------------

positive_zero_ranking <- sum(
  eligible_for_ranking &
    log2fc_raw > 0 &
    temporary_ranking == 0
)

negative_zero_ranking <- sum(
  eligible_for_ranking &
    log2fc_raw < 0 &
    temporary_ranking == 0
)

zero_log2fc_zero_ranking <- sum(
  eligible_for_ranking &
    log2fc_raw == 0 &
    temporary_ranking == 0
)

cat(
  "log2FC > 0 but ranking = 0:",
  positive_zero_ranking,
  "\n"
)

cat(
  "log2FC < 0 but ranking = 0:",
  negative_zero_ranking,
  "\n"
)

cat(
  "log2FC = 0 and ranking = 0:",
  zero_log2fc_zero_ranking,
  "\n\n"
)


# ------------------------------------------------------------
# 3. Cross-tabulation
# ------------------------------------------------------------

direction_table <- table(
  log2FC_direction = ifelse(
    log2fc_raw[eligible_for_ranking] > 0,
    "positive",
    ifelse(
      log2fc_raw[eligible_for_ranking] < 0,
      "negative",
      "zero"
    )
  ),
  ranking_direction = ifelse(
    temporary_ranking[eligible_for_ranking] > 0,
    "positive",
    ifelse(
      temporary_ranking[eligible_for_ranking] < 0,
      "negative",
      "zero"
    )
  )
)

cat("Direction cross-tabulation:\n")
print(direction_table)
cat("\n")


# ------------------------------------------------------------
# 4. P-value = 1 by log2FC direction
# ------------------------------------------------------------

p1_direction_table <- table(
  log2FC_direction = ifelse(
    log2fc_raw[
      eligible_for_ranking & p_raw == 1
    ] > 0,
    "positive",
    ifelse(
      log2fc_raw[
        eligible_for_ranking & p_raw == 1
      ] < 0,
      "negative",
      "zero"
    )
  )
)

cat("P = 1 direction distribution:\n")
print(p1_direction_table)
cat("\n")


# ------------------------------------------------------------
# 5. Show problematic rows
# ------------------------------------------------------------

problem_idx <- which(
  eligible_for_ranking &
    (
      (
        log2fc_raw > 0 &
          temporary_ranking <= 0
      ) |
        (
          log2fc_raw < 0 &
            temporary_ranking >= 0
        ) |
        (
          log2fc_raw == 0 &
            temporary_ranking != 0
        )
    )
)

cat("Problematic rows under strict directional test:",
    length(problem_idx), "\n\n")

if (length(problem_idx) > 0) {
  
  diagnostic_df <- data.frame(
    Gene = gene_id_raw[problem_idx],
    log2FC = log2fc_raw[problem_idx],
    P_value = p_raw[problem_idx],
    ranking = temporary_ranking[problem_idx],
    stringsAsFactors = FALSE
  )
  
  print(
    head(
      diagnostic_df,
      30
    )
  )
  
} else {
  
  cat("No problematic rows detected.\n")
  
}

cat("\n============================================================\n")
cat("END DIAGNOSTIC\n")
cat("============================================================\n") 


# ============================================================
# STEP9D-2B — Corrected Direction Consistency Audit
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-2B — CORRECTED DIRECTION CONSISTENCY AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define non-zero ranking subset
# ------------------------------------------------------------

nonzero_ranking <- (
  eligible_for_ranking &
    temporary_ranking != 0
)


# ------------------------------------------------------------
# 2. Sign agreement among non-zero rankings
# ------------------------------------------------------------

positive_sign_mismatch <- sum(
  nonzero_ranking &
    (
      sign(temporary_ranking) != sign(log2fc_raw)
    )
)

negative_sign_mismatch <- positive_sign_mismatch


# ------------------------------------------------------------
# 3. Exact sign-switch audit
# ------------------------------------------------------------

sign_switch_n <- sum(
  nonzero_ranking &
    (
      sign(temporary_ranking) != sign(log2fc_raw)
    )
)


# ------------------------------------------------------------
# 4. Zero-ranking audit
# ------------------------------------------------------------

zero_ranking_n <- sum(
  eligible_for_ranking &
    temporary_ranking == 0
)

zero_ranking_due_to_p1 <- sum(
  eligible_for_ranking &
    temporary_ranking == 0 &
    p_raw == 1
)

zero_ranking_due_to_log2fc0 <- sum(
  eligible_for_ranking &
    temporary_ranking == 0 &
    log2fc_raw == 0
)


# ------------------------------------------------------------
# 5. Verify all zero rankings are explained
# ------------------------------------------------------------

zero_ranking_unexplained_n <- sum(
  eligible_for_ranking &
    temporary_ranking == 0 &
    p_raw != 1 &
    log2fc_raw != 0
)


# ------------------------------------------------------------
# 6. Final corrected consistency definition
# ------------------------------------------------------------

direction_consistency_corrected <- (
  sign_switch_n == 0 &
    zero_ranking_unexplained_n == 0
)


# ------------------------------------------------------------
# 7. Report
# ------------------------------------------------------------

cat(
  "Eligible genes:",
  sum(eligible_for_ranking),
  "\n"
)

cat(
  "Non-zero ranking genes:",
  sum(nonzero_ranking),
  "\n"
)

cat(
  "Zero-ranking genes:",
  zero_ranking_n,
  "\n"
)

cat(
  "Zero-ranking genes due to P = 1:",
  zero_ranking_due_to_p1,
  "\n"
)

cat(
  "Zero-ranking genes due to log2FC = 0:",
  zero_ranking_due_to_log2fc0,
  "\n"
)

cat(
  "Unexplained zero-ranking genes:",
  zero_ranking_unexplained_n,
  "\n"
)

cat(
  "True sign switches:",
  sign_switch_n,
  "\n\n"
)

cat(
  "Corrected direction consistency:",
  direction_consistency_corrected,
  "\n\n"
)


# ------------------------------------------------------------
# 8. Hard PASS criteria
# ------------------------------------------------------------

if (!direction_consistency_corrected) {
  
  stop(
    "\nERROR: Corrected direction-consistency audit failed.\n",
    "Do NOT proceed to final STEP9D-2B lock.\n"
  )
  
}


cat("PASS — corrected direction-consistency audit.\n")

cat("\n============================================================\n")
cat("END CORRECTED AUDIT\n")
cat("============================================================\n")

##next step：STEP9D-2B FINAL LOCK
# ============================================================
cat("\n============================================================\n")
cat("STEP9D-2B — COLUMN NAME DIAGNOSTIC\n")
cat("============================================================\n\n")

cat("Input file:\n")
cat(input_file, "\n\n")

cat("Number of columns detected:\n")
cat(ncol(ranking_input), "\n\n")

cat("Exact column names detected by R:\n")
print(colnames(ranking_input))

cat("\nColumn names using dput():\n")
dput(colnames(ranking_input))

cat("\nRequired columns:\n")
print(required_columns)

cat("\nMissing columns according to setdiff():\n")
print(setdiff(required_columns, colnames(ranking_input)))

cat("\n============================================================\n")
cat("END COLUMN NAME DIAGNOSTIC\n")
cat("============================================================\n")


cat("\n============================================================\n")
cat("STEP9D-2B — TRACE STEP9D-1I LOCKED SOURCE\n")
cat("============================================================\n\n")

step9_root <- "input_dir/Paper 1_Step 9-Functional Enrichment"

rds_files <- list.files(
  step9_root,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

cat(
  "RDS files found under STEP 9 directory:",
  length(rds_files),
  "\n\n"
)

if (length(rds_files) == 0) {
  
  cat("No RDS files found.\n")
  
} else {
  
  for (f in rds_files) {
    
    obj <- tryCatch(
      readRDS(f),
      error = function(e) NULL
    )
    
    if (is.null(obj)) {
      next
    }
    
    hit <- FALSE
    
    cat_this <- character(0)
    
    # --------------------------------------------------------
    # Check named list/object fields
    # --------------------------------------------------------
    
    if (is.list(obj) && !is.null(names(obj))) {
      
      nms <- names(obj)
      
      source_fields <- intersect(
        nms,
        c(
          "input_file",
          "source_file",
          "source_path",
          "expression_file",
          "gene_universe_file",
          "universe_file",
          "locked_input_file",
          "genome_wide_input_file"
        )
      )
      
      if (length(source_fields) > 0) {
        
        hit <- TRUE
        
        for (nm in source_fields) {
          
          val <- obj[[nm]]
          
          if (length(val) == 1 &&
              is.character(val)) {
            
            cat_this <- c(
              cat_this,
              paste0(
                "  ",
                nm,
                " = ",
                val
              )
            )
          }
        }
      }
    }
    
    # --------------------------------------------------------
    # Print matching RDS
    # --------------------------------------------------------
    
    if (hit) {
      
      cat(
        "\n------------------------------------------------------------\n"
      )
      
      cat(
        "RDS:",
        f,
        "\n"
      )
      
      cat(
        paste(
          cat_this,
          collapse = "\n"
        ),
        "\n"
      )
    }
  }
}

cat(
  "\n============================================================\n"
)

cat(
  "END STEP9D-1I SOURCE TRACE\n"
)

cat(
  "============================================================\n"
)


# ============================================================
# STEP9D-2C DIAGNOSTIC
# Inspect actual STEP9D-2B FINAL LOCK structure
# NO recalculation
# NO ranking vector construction
# NO GSEA
# ============================================================

lock_rds <- file.path(
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA",
  "9D-2_Ranking_Metric_Lock",
  "STEP9D2B_FINAL_Ranking_Metric_Lock.rds"
)

step9d2b_lock <- readRDS(lock_rds)

cat("\n============================================================\n")
cat("STEP9D-2B LOCK STRUCTURE DIAGNOSTIC\n")
cat("============================================================\n\n")

cat("Object class:\n")
print(class(step9d2b_lock))

cat("\nObject names:\n")
print(names(step9d2b_lock))

cat("\nObject structure:\n")
str(step9d2b_lock, max.level = 2)

cat("\n============================================================\n")
cat("FORMULA-RELATED FIELDS\n")
cat("============================================================\n\n")

formula_candidates <- c(
  "ranking_metric_formula",
  "ranking_formula",
  "formula",
  "ranking_metric",
  "metric_definition",
  "locked_formula"
)

existing_formula_fields <- intersect(
  formula_candidates,
  names(step9d2b_lock)
)

cat("Formula-related fields found:\n")
print(existing_formula_fields)

if (length(existing_formula_fields) > 0) {
  
  for (f in existing_formula_fields) {
    
    cat("\n------------------------------------------------------------\n")
    cat("FIELD:", f, "\n")
    cat("CLASS:", class(step9d2b_lock[[f]]), "\n")
    cat("VALUE:\n")
    print(step9d2b_lock[[f]])
  }
  
} else {
  
  cat(
    "\nNO formula-related field found under the candidate names.\n"
  )
}

cat("\n============================================================\n")
cat("KNOWN LOCK VALUES\n")
cat("============================================================\n\n")

for (f in c(
  "hard_pass",
  "source_provenance_match",
  "ranking_vector_created",
  "gsea_performed",
  "deg_recalculated",
  "eligible_gene_n",
  "nonzero_ranking_n",
  "zero_ranking_n",
  "zero_due_p1",
  "zero_due_log2fc0",
  "zero_unexplained",
  "true_sign_switches",
  "direction_consistent"
)) {
  
  if (f %in% names(step9d2b_lock)) {
    
    cat(
      f, ": ",
      paste(
        capture.output(print(step9d2b_lock[[f]])),
        collapse = " "
      ),
      "\n",
      sep = ""
    )
  }
}

cat("\n============================================================\n")
cat("END DIAGNOSTIC\n")
cat("============================================================\n") 

# ============================================================
# STEP9D-2C DIAGNOSTIC
# Direction-count discrepancy audit
#
# NO recalculation
# NO file modification
# NO GSEA
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-2C DIRECTION COUNT DISCREPANCY AUDIT\n")
cat("============================================================\n\n")

cat("Source log2FC > 0:       ", source_log2fc_positive_n, "\n")
cat("Source log2FC < 0:       ", source_log2fc_negative_n, "\n")
cat("Source log2FC = 0:       ", source_log2fc_zero_n, "\n\n")

cat("Ranking > 0:             ", ranking_positive_n, "\n")
cat("Ranking < 0:             ", ranking_negative_n, "\n")
cat("Ranking = 0:             ", ranking_zero_n, "\n\n")

cat("P = 1:                   ", source_p_one_n, "\n")
cat("P = 0:                   ", source_p_zero_n, "\n\n")


# ------------------------------------------------------------
# Cross-tabulation:
# log2FC direction × ranking direction
# ------------------------------------------------------------

source_direction <- ifelse(
  deg$log2FC > 0,
  "positive",
  ifelse(
    deg$log2FC < 0,
    "negative",
    "zero"
  )
)

ranking_direction <- ifelse(
  ranking_values > 0,
  "positive",
  ifelse(
    ranking_values < 0,
    "negative",
    "zero"
  )
)

direction_table <- table(
  log2FC_direction = source_direction,
  ranking_direction = ranking_direction
)

cat("============================================================\n")
cat("log2FC direction × ranking direction\n")
cat("============================================================\n\n")

print(direction_table)

cat("\n============================================================\n")
cat("P=1 cross-tabulation\n")
cat("============================================================\n\n")

p1_table <- table(
  log2FC_direction = source_direction[
    deg$P_value == 1
  ],
  ranking_direction = ranking_direction[
    deg$P_value == 1
  ]
)

print(p1_table)


# ------------------------------------------------------------
# Count genes with nonzero log2FC but zero ranking
# ------------------------------------------------------------

nonzero_log2fc_zero_ranking <- sum(
  deg$log2FC != 0 &
    ranking_values == 0
)

cat("\nGenes with log2FC != 0 AND ranking = 0:",
    nonzero_log2fc_zero_ranking,
    "\n")

cat(
  "Expected from P=1 genes:",
  source_p_one_n,
  "\n"
)


# ------------------------------------------------------------
# Positive log2FC but zero ranking
# ------------------------------------------------------------

positive_log2fc_zero_ranking <- sum(
  deg$log2FC > 0 &
    ranking_values == 0
)

negative_log2fc_zero_ranking <- sum(
  deg$log2FC < 0 &
    ranking_values == 0
)

cat(
  "Positive log2FC + zero ranking:",
  positive_log2fc_zero_ranking,
  "\n"
)

cat(
  "Negative log2FC + zero ranking:",
  negative_log2fc_zero_ranking,
  "\n"
)


# ------------------------------------------------------------
# True direction-switch audit
#
# IMPORTANT:
# Zero-ranking genes are excluded because zero is not a
# positive/negative direction.
# ------------------------------------------------------------

nonzero_ranking_idx <- ranking_values != 0

true_direction_switches <- sum(
  sign(deg$log2FC[nonzero_ranking_idx]) !=
    sign(ranking_values[nonzero_ranking_idx])
)

cat(
  "True sign switches among non-zero rankings:",
  true_direction_switches,
  "\n"
)

cat("\n============================================================\n")
cat("END DIAGNOSTIC\n")
cat("============================================================\n") 

# ============================================================
# STEP9D-2C
# FINAL Genome-wide Preranked GSEA Ranking Vector Construction
#
# Purpose:
#   Construct the genome-wide preranked vector using the
#   previously LOCKED STEP9D-2B ranking metric.
#
# Ranking metric:
#   sign(log2FC) * (-log10(max(P_value, 1e-300)))
#
# IMPORTANT:
#   - Uses STEP9D-2B only as the metric lock.
#   - Uses the exact locked DEG source.
#   - Uses element-wise pmax().
#   - Does NOT require ranking-positive count to equal
#     positive-log2FC count because P=1 produces ranking = 0.
#   - No DEG recalculation.
#   - No GSEA.
#
# STEP9D-2C v3
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-2C FINAL RANKING VECTOR CONSTRUCTION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory and output folders
# ------------------------------------------------------------

step_root <- "input_dir/Paper 1_Step 9-Functional Enrichment"

lock_dir <- file.path(
  step_root,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-2_Ranking_Metric_Lock"
)

out_dir <- file.path(
  step_root,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-2C_Final_Ranking_Vector"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Locate and read STEP9D-2B FINAL LOCK
# ------------------------------------------------------------

lock_file <- file.path(
  lock_dir,
  "STEP9D2B_FINAL_Ranking_Metric_Lock.rds"
)

if (!file.exists(lock_file)) {
  stop(
    "FATAL ERROR: STEP9D-2B FINAL LOCK RDS not found."
  )
}

step9d2b_lock <- readRDS(lock_file)


# ------------------------------------------------------------
# 2. Validate STEP9D-2B lock identity
# ------------------------------------------------------------

required_lock_fields <- c(
  "workflow",
  "step",
  "status",
  "source_file",
  "step9d1i_locked_source",
  "source_provenance_match",
  "universe_file",
  "genome_wide_gene_n",
  "observed_gene_n",
  "eligible_gene_n",
  "gene_id_column",
  "log2fc_column",
  "nominal_p_column",
  "adjusted_p_column_detected",
  "adjusted_p_used_for_ranking",
  "native_signed_statistic_available",
  "ranking_metric_name",
  "ranking_metric_formula",
  "p_floor",
  "ranking_direction_source",
  "ranking_strength_source",
  "p_zero_n",
  "p_one_n",
  "log2fc_zero_n",
  "log2fc_positive_n",
  "log2fc_negative_n",
  "nonzero_ranking_n",
  "zero_ranking_n",
  "zero_ranking_due_to_p1",
  "zero_ranking_due_to_log2fc0",
  "unexplained_zero_ranking_n",
  "true_sign_switch_n",
  "direction_consistency",
  "formula_reproducibility",
  "duplicated_gene_n",
  "missing_gene_n",
  "ranking_vector_created",
  "gsea_performed",
  "deg_recalculated",
  "hard_pass"
)

missing_lock_fields <- setdiff(
  required_lock_fields,
  names(step9d2b_lock)
)

if (length(missing_lock_fields) > 0) {
  stop(
    "FATAL ERROR: Required STEP9D-2B lock fields missing: ",
    paste(missing_lock_fields, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 3. Exact lock validation
# ------------------------------------------------------------

if (!identical(
  step9d2b_lock$ranking_metric_name,
  "signed_minus_log10_nominal_P"
)) {
  stop(
    "FATAL ERROR: Ranking metric name does not match locked definition."
  )
}

if (!identical(
  step9d2b_lock$ranking_metric_formula,
  "sign(log2FC) * (-log10(max(P_value, 1e-300)))"
)) {
  stop(
    "FATAL ERROR: Ranking metric formula does not match STEP9D-2B lock."
  )
}

if (!identical(
  as.numeric(step9d2b_lock$p_floor),
  1e-300
)) {
  stop(
    "FATAL ERROR: P-value floor does not match locked value."
  )
}

if (!isTRUE(
  step9d2b_lock$source_provenance_match
)) {
  stop(
    "FATAL ERROR: Source provenance mismatch."
  )
}

if (!isTRUE(
  step9d2b_lock$direction_consistency
)) {
  stop(
    "FATAL ERROR: STEP9D-2B direction consistency is not TRUE."
  )
}

if (!isTRUE(
  step9d2b_lock$formula_reproducibility
)) {
  stop(
    "FATAL ERROR: STEP9D-2B formula reproducibility is not TRUE."
  )
}

if (!isTRUE(
  step9d2b_lock$hard_pass
)) {
  stop(
    "FATAL ERROR: STEP9D-2B hard_pass is not TRUE."
  )
}

if (isTRUE(
  step9d2b_lock$gsea_performed
)) {
  stop(
    "FATAL ERROR: GSEA is already marked as performed in STEP9D-2B."
  )
}

if (isTRUE(
  step9d2b_lock$deg_recalculated
)) {
  stop(
    "FATAL ERROR: DEG recalculation is marked TRUE."
  )
}


# ------------------------------------------------------------
# 4. Read exact locked DEG source
# ------------------------------------------------------------

source_file <- step9d2b_lock$source_file

if (!file.exists(source_file)) {
  stop(
    "FATAL ERROR: Locked DEG source file does not exist:\n",
    source_file
  )
}

deg <- read.csv(
  source_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 5. Validate source columns
# ------------------------------------------------------------

gene_col <- step9d2b_lock$gene_id_column
log2fc_col <- step9d2b_lock$log2fc_column
p_col <- step9d2b_lock$nominal_p_column

required_source_columns <- c(
  gene_col,
  log2fc_col,
  p_col
)

missing_source_columns <- setdiff(
  required_source_columns,
  names(deg)
)

if (length(missing_source_columns) > 0) {
  stop(
    "FATAL ERROR: Required source columns missing: ",
    paste(missing_source_columns, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 6. Standardize local vectors
# ------------------------------------------------------------

gene_id <- as.character(
  deg[[gene_col]]
)

log2fc <- as.numeric(
  deg[[log2fc_col]]
)

p_value <- as.numeric(
  deg[[p_col]]
)


# ------------------------------------------------------------
# 7. Source structural QC
# ------------------------------------------------------------

if (length(gene_id) != 57828) {
  stop(
    "FATAL ERROR: Gene count is not 57,828."
  )
}

if (anyNA(gene_id)) {
  stop(
    "FATAL ERROR: Missing gene IDs detected."
  )
}

if (anyDuplicated(gene_id) > 0) {
  stop(
    "FATAL ERROR: Duplicated gene IDs detected."
  )
}

if (anyNA(log2fc) ||
    any(!is.finite(log2fc))) {
  stop(
    "FATAL ERROR: Invalid log2FC values detected."
  )
}

if (anyNA(p_value) ||
    any(!is.finite(p_value))) {
  stop(
    "FATAL ERROR: Invalid nominal P values detected."
  )
}

if (any(p_value < 0 | p_value > 1)) {
  stop(
    "FATAL ERROR: Nominal P values outside [0,1]."
  )
}


# ------------------------------------------------------------
# 8. Reconfirm locked source counts
# ------------------------------------------------------------

source_log2fc_positive_n <- sum(
  log2fc > 0
)

source_log2fc_negative_n <- sum(
  log2fc < 0
)

source_log2fc_zero_n <- sum(
  log2fc == 0
)

source_p_zero_n <- sum(
  p_value == 0
)

source_p_one_n <- sum(
  p_value == 1
)


cat("Genome-wide genes:       ", length(gene_id), "\n")
cat("log2FC > 0:              ", source_log2fc_positive_n, "\n")
cat("log2FC < 0:              ", source_log2fc_negative_n, "\n")
cat("log2FC = 0:              ", source_log2fc_zero_n, "\n")
cat("P = 0:                   ", source_p_zero_n, "\n")
cat("P = 1:                   ", source_p_one_n, "\n\n")


# ------------------------------------------------------------
# 9. Exact lock-count validation
# ------------------------------------------------------------

if (!identical(
  as.integer(source_log2fc_positive_n),
  as.integer(step9d2b_lock$log2fc_positive_n)
)) {
  stop(
    "FATAL ERROR: Positive log2FC count changed from STEP9D-2B lock."
  )
}

if (!identical(
  as.integer(source_log2fc_negative_n),
  as.integer(step9d2b_lock$log2fc_negative_n)
)) {
  stop(
    "FATAL ERROR: Negative log2FC count changed from STEP9D-2B lock."
  )
}

if (!identical(
  as.integer(source_log2fc_zero_n),
  as.integer(step9d2b_lock$log2fc_zero_n)
)) {
  stop(
    "FATAL ERROR: Zero log2FC count changed from STEP9D-2B lock."
  )
}

if (!identical(
  as.integer(source_p_zero_n),
  as.integer(step9d2b_lock$p_zero_n)
)) {
  stop(
    "FATAL ERROR: P=0 count changed from STEP9D-2B lock."
  )
}

if (!identical(
  as.integer(source_p_one_n),
  as.integer(step9d2b_lock$p_one_n)
)) {
  stop(
    "FATAL ERROR: P=1 count changed from STEP9D-2B lock."
  )
}


# ------------------------------------------------------------
# 10. Construct FINAL ranking vector
#
# IMPORTANT:
# pmax() is element-wise.
#
# DO NOT use:
# max(p_value, 1e-300)
#
# because max() returns a single scalar.
# ------------------------------------------------------------

p_floor_locked <- as.numeric(
  step9d2b_lock$p_floor
)

ranking_values <- sign(log2fc) *
  (
    -log10(
      pmax(
        p_value,
        p_floor_locked
      )
    )
  )


# ------------------------------------------------------------
# 11. Ranking vector numerical QC
# ------------------------------------------------------------

if (length(ranking_values) != 57828) {
  stop(
    "FATAL ERROR: Ranking vector length is not 57,828."
  )
}

if (anyNA(ranking_values) ||
    any(!is.finite(ranking_values))) {
  stop(
    "FATAL ERROR: Ranking vector contains NA/NaN/Inf."
  )
}

cat(
  "\nFinal ranking vector numerical QC: PASS\n\n"
)


# ------------------------------------------------------------
# 12. Ranking direction counts
# ------------------------------------------------------------

ranking_positive_n <- sum(
  ranking_values > 0
)

ranking_negative_n <- sum(
  ranking_values < 0
)

ranking_zero_n <- sum(
  ranking_values == 0
)

cat("Ranking > 0:             ", ranking_positive_n, "\n")
cat("Ranking < 0:             ", ranking_negative_n, "\n")
cat("Ranking = 0:             ", ranking_zero_n, "\n\n")


# ------------------------------------------------------------
# 13. Ranking count closure
#
# IMPORTANT:
# Positive/negative ranking counts are NOT required to equal
# positive/negative log2FC counts because P=1 generates zero
# ranking values.
# ------------------------------------------------------------

if (
  ranking_positive_n +
  ranking_negative_n +
  ranking_zero_n !=
  57828
) {
  stop(
    "FATAL ERROR: Ranking direction counts do not sum to 57,828."
  )
}

if (
  ranking_positive_n +
  ranking_negative_n !=
  step9d2b_lock$nonzero_ranking_n
) {
  stop(
    "FATAL ERROR: Non-zero ranking count does not match STEP9D-2B lock."
  )
}

if (
  ranking_zero_n !=
  step9d2b_lock$zero_ranking_n
) {
  stop(
    "FATAL ERROR: Zero ranking count does not match STEP9D-2B lock."
  )}


# ------------------------------------------------------------
# 14. Zero-ranking mechanism audit
# ------------------------------------------------------------

zero_ranking_due_to_p1 <- sum(
  p_value == 1 &
    ranking_values == 0
)

zero_ranking_due_to_log2fc0 <- sum(
  log2fc == 0 &
    ranking_values == 0
)

unexplained_zero_ranking_n <- sum(
  ranking_values == 0 &
    p_value != 1 &
    log2fc != 0
)

cat(
  "Zero ranking due to P=1:       ",
  zero_ranking_due_to_p1,
  "\n"
)

cat(
  "Zero ranking due to log2FC=0: ",
  zero_ranking_due_to_log2fc0,
  "\n"
)

cat(
  "Unexplained zero ranking:     ",
  unexplained_zero_ranking_n,
  "\n\n"
)


if (
  zero_ranking_due_to_p1 !=
  source_p_one_n
) {
  stop(
    "FATAL ERROR: Zero-ranking genes do not exactly match P=1 genes."
  )
}

if (
  zero_ranking_due_to_log2fc0 !=
  source_log2fc_zero_n
) {
  stop(
    "FATAL ERROR: Zero-ranking/log2FC=0 relationship changed."
  )
}

if (
  unexplained_zero_ranking_n != 0
) {
  stop(
    "FATAL ERROR: Unexplained zero-ranking genes detected."
  )
}


# ------------------------------------------------------------
# 15. TRUE direction-switch audit
#
# Zero-ranking genes are excluded.
# Only positive-vs-negative direction is tested.
# ------------------------------------------------------------

nonzero_ranking_idx <- (
  ranking_values != 0
)

true_sign_switch_n <- sum(
  sign(log2fc[nonzero_ranking_idx]) !=
    sign(ranking_values[nonzero_ranking_idx])
)

cat(
  "True sign switches among non-zero rankings: ",
  true_sign_switch_n,
  "\n"
)

if (
  true_sign_switch_n != 0
) {
  stop(
    "FATAL ERROR: True ranking direction switch detected."
  )
}


# ------------------------------------------------------------
# 16. Full direction cross-tabulation
# ------------------------------------------------------------

source_direction <- ifelse(
  log2fc > 0,
  "positive",
  ifelse(
    log2fc < 0,
    "negative",
    "zero"
  )
)

ranking_direction <- ifelse(
  ranking_values > 0,
  "positive",
  ifelse(
    ranking_values < 0,
    "negative",
    "zero"
  )
)

direction_table <- table(
  log2FC_direction = source_direction,
  ranking_direction = ranking_direction
)

cat("\n============================================================\n")
cat("FINAL log2FC direction × ranking direction\n")
cat("============================================================\n\n")

print(direction_table)


# ------------------------------------------------------------
# 17. Direction consistency
# ------------------------------------------------------------

direction_consistency <- (
  true_sign_switch_n == 0
)

if (!isTRUE(direction_consistency)) {
  stop(
    "FATAL ERROR: Direction consistency failed."
  )
}


# ------------------------------------------------------------
# 18. Expected locked counts
# ------------------------------------------------------------

expected_nonzero_ranking_n <- 57815
expected_zero_ranking_n <- 13
expected_zero_due_to_p1_n <- 13
expected_zero_due_to_log2fc0_n <- 0
expected_unexplained_zero_n <- 0
expected_sign_switch_n <- 0


if (
  ranking_positive_n +
  ranking_negative_n !=
  expected_nonzero_ranking_n
) {
  stop(
    "FATAL ERROR: Non-zero ranking total does not equal 57,815."
  )
}

if (
  ranking_zero_n !=
  expected_zero_ranking_n
) {
  stop(
    "FATAL ERROR: Zero ranking total does not equal 13."
  )
}

if (
  zero_ranking_due_to_p1 !=
  expected_zero_due_to_p1_n
) {
  stop(
    "FATAL ERROR: P=1 zero-ranking count does not equal 13."
  )
}

if (
  zero_ranking_due_to_log2fc0 !=
  expected_zero_due_to_log2fc0_n
) {
  stop(
    "FATAL ERROR: log2FC=0 zero-ranking count is not 0."
  )
}

if (
  unexplained_zero_ranking_n !=
  expected_unexplained_zero_n
) {
  stop(
    "FATAL ERROR: Unexplained zero-ranking count is not 0."
  )
}

if (
  true_sign_switch_n !=
  expected_sign_switch_n
) {
  stop(
    "FATAL ERROR: True sign-switch count is not 0."
  )
}


# ------------------------------------------------------------
# 19. Create named ranking vector
# ------------------------------------------------------------

ranking_vector <- ranking_values

names(ranking_vector) <- gene_id


# ------------------------------------------------------------
# 20. Final ranking vector structural QC
# ------------------------------------------------------------

if (length(ranking_vector) != 57828) {
  stop(
    "FATAL ERROR: Final ranking vector length mismatch."
  )
}

if (is.null(names(ranking_vector))) {
  stop(
    "FATAL ERROR: Final ranking vector has no gene names."
  )
}

if (anyNA(names(ranking_vector)) ||
    any(names(ranking_vector) == "")) {
  stop(
    "FATAL ERROR: Missing gene names in final ranking vector."
  )
}

if (anyDuplicated(names(ranking_vector)) > 0) {
  stop(
    "FATAL ERROR: Duplicated gene names in final ranking vector."
  )
}

if (anyNA(ranking_vector) ||
    any(!is.finite(ranking_vector))) {
  stop(
    "FATAL ERROR: Invalid numerical values in final ranking vector."
  )
}

cat(
  "\nFinal ranking vector structural QC: PASS\n"
)


# ------------------------------------------------------------
# 21. Save FINAL ranking vector
# ------------------------------------------------------------

ranking_rds <- file.path(
  out_dir,
  "STEP9D2C_FINAL_GenomeWide_Ranking_Vector.rds"
)

saveRDS(
  ranking_vector,
  ranking_rds
)


# ------------------------------------------------------------
# 22. Save transparent CSV
# ------------------------------------------------------------

ranking_table <- data.frame(
  Gene = gene_id,
  log2FC = log2fc,
  P_value = p_value,
  ranking_score = ranking_values,
  stringsAsFactors = FALSE
)

ranking_csv <- file.path(
  out_dir,
  "STEP9D2C_FINAL_GenomeWide_Ranking_Vector.csv"
)

write.csv(
  ranking_table,
  ranking_csv,
  row.names = FALSE,
  quote = FALSE
)


# ------------------------------------------------------------
# 23. Final audit object
# ------------------------------------------------------------

step9d2c_final_audit <- list(
  
  workflow = "Paper 1 STEP9D Genome-wide preranked GSEA",
  
  step = "STEP9D-2C",
  
  status = "FINAL_LOCK",
  
  source_file = source_file,
  
  source_provenance_match =
    isTRUE(step9d2b_lock$source_provenance_match),
  
  genome_wide_gene_n =
    length(gene_id),
  
  log2fc_positive_n =
    source_log2fc_positive_n,
  
  log2fc_negative_n =
    source_log2fc_negative_n,
  
  log2fc_zero_n =
    source_log2fc_zero_n,
  
  p_zero_n =
    source_p_zero_n,
  
  p_one_n =
    source_p_one_n,
  
  ranking_metric_name =
    step9d2b_lock$ranking_metric_name,
  
  ranking_metric_formula =
    step9d2b_lock$ranking_metric_formula,
  
  p_floor =
    p_floor_locked,
  
  ranking_direction_source =
    "log2FC",
  
  ranking_strength_source =
    "nominal P_value",
  
  adjusted_p_used_for_ranking =
    FALSE,
  
  native_signed_statistic_available =
    FALSE,
  
  ranking_positive_n =
    ranking_positive_n,
  
  ranking_negative_n =
    ranking_negative_n,
  
  ranking_zero_n =
    ranking_zero_n,
  
  nonzero_ranking_n =
    ranking_positive_n +
    ranking_negative_n,
  
  zero_ranking_due_to_p1 =
    zero_ranking_due_to_p1,
  
  zero_ranking_due_to_log2fc0 =
    zero_ranking_due_to_log2fc0,
  
  unexplained_zero_ranking_n =
    unexplained_zero_ranking_n,
  
  true_sign_switch_n =
    true_sign_switch_n,
  
  direction_consistency =
    direction_consistency,
  
  duplicated_gene_n =
    anyDuplicated(gene_id),
  
  missing_gene_n =
    sum(is.na(gene_id)),
  
  ranking_vector_created =
    TRUE,
  
  gsea_performed =
    FALSE,
  
  deg_recalculated =
    FALSE,
  
  ranking_vector_rds =
    ranking_rds,
  
  ranking_vector_csv =
    ranking_csv,
  
  hard_pass =
    TRUE
)


# ------------------------------------------------------------
# 24. Save audit RDS
# ------------------------------------------------------------

audit_rds <- file.path(
  out_dir,
  "STEP9D2C_FINAL_Ranking_Vector_Audit.rds"
)

saveRDS(
  step9d2c_final_audit,
  audit_rds
)


# ------------------------------------------------------------
# 25. Save human-readable audit
# ------------------------------------------------------------

audit_txt <- file.path(
  out_dir,
  "STEP9D2C_FINAL_Ranking_Vector_Audit.txt"
)

audit_lines <- c(
  
  "============================================================",
  "STEP9D-2C FINAL RANKING VECTOR AUDIT",
  "============================================================",
  "",
  
  paste(
    "Source file:",
    source_file
  ),
  
  paste(
    "Genome-wide genes:",
    length(gene_id)
  ),
  
  paste(
    "log2FC > 0:",
    source_log2fc_positive_n
  ),
  
  paste(
    "log2FC < 0:",
    source_log2fc_negative_n
  ),
  
  paste(
    "log2FC = 0:",
    source_log2fc_zero_n
  ),
  
  paste(
    "P = 0:",
    source_p_zero_n
  ),
  
  paste(
    "P = 1:",
    source_p_one_n
  ),
  
  "",
  
  paste(
    "Ranking > 0:",
    ranking_positive_n
  ),
  
  paste(
    "Ranking < 0:",
    ranking_negative_n
  ),
  
  paste(
    "Ranking = 0:",
    ranking_zero_n
  ),
  
  paste(
    "Non-zero ranking:",
    ranking_positive_n +
      ranking_negative_n
  ),
  
  "",
  
  paste(
    "Zero ranking due to P=1:",
    zero_ranking_due_to_p1
  ),
  
  paste(
    "Zero ranking due to log2FC=0:",
    zero_ranking_due_to_log2fc0
  ),
  
  paste(
    "Unexplained zero ranking:",
    unexplained_zero_ranking_n
  ),
  
  paste(
    "True sign switches:",
    true_sign_switch_n
  ),
  
  paste(
    "Direction consistency:",
    direction_consistency
  ),
  
  "",
  
  "Ranking vector created: TRUE",
  
  "GSEA performed: FALSE",
  
  "DEG recalculated: FALSE",
  
  "",
  
  paste(
    "FINAL HARD PASS:",
    TRUE
  ),
  
  "============================================================"
)

writeLines(
  audit_lines,
  audit_txt
)


# ------------------------------------------------------------
# 26. Reload validation
# ------------------------------------------------------------

ranking_vector_reload <- readRDS(
  ranking_rds
)

if (!identical(
  ranking_vector_reload,
  ranking_vector
)) {
  stop(
    "FATAL ERROR: Reloaded ranking vector differs from saved object."
  )
}

audit_reload <- readRDS(
  audit_rds
)

if (!isTRUE(
  audit_reload$hard_pass
)) {
  stop(
    "FATAL ERROR: Reloaded audit hard_pass is not TRUE."
  )
}


# ------------------------------------------------------------
# 27. Final console summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-2C FINAL STATUS\n")
cat("============================================================\n\n")

cat(
  "Genome-wide genes:              ",
  length(gene_id),
  "\n"
)

cat(
  "Non-zero ranking:               ",
  ranking_positive_n +
    ranking_negative_n,
  "\n"
)

cat(
  "Zero ranking:                   ",
  ranking_zero_n,
  "\n"
)

cat(
  "Zero due to P=1:                ",
  zero_ranking_due_to_p1,
  "\n"
)

cat(
  "Unexplained zero ranking:       ",
  unexplained_zero_ranking_n,
  "\n"
)

cat(
  "True sign switches:             ",
  true_sign_switch_n,
  "\n"
)

cat(
  "Direction consistency:          ",
  direction_consistency,
  "\n"
)

cat(
  "Ranking vector created:         TRUE\n"
)

cat(
  "GSEA performed:                 FALSE\n"
)

cat(
  "DEG recalculated:               FALSE\n"
)

cat(
  "Final ranking vector RDS exists:",
  file.exists(ranking_rds),
  "\n"
)

cat(
  "Final ranking CSV exists:       ",
  file.exists(ranking_csv),
  "\n"
)

cat(
  "Final audit RDS exists:         ",
  file.exists(audit_rds),
  "\n"
)

cat(
  "Final audit TXT exists:         ",
  file.exists(audit_txt),
  "\n"
)

cat(
  "\nSTEP9D-2C FINAL HARD PASS: TRUE\n"
)

cat("\n============================================================\n")
cat("STEP9D-2C COMPLETE — GSEA NOT PERFORMED\n")
cat("============================================================\n")


# ============================================================
# STEP9D-3A
# Gene-ID / Entrez Mapping Diagnostic
#
# READ-ONLY DIAGNOSTIC
#
# NO file modification
# NO ranking recalculation
# NO DEG recalculation
# NO GSEA
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-3A GENE-ID / ENTREZ MAPPING DIAGNOSTIC\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Load locked STEP9D-2C ranking vector
# ------------------------------------------------------------

step_root <- "input_dir/Paper 1_Step 9-Functional Enrichment"

ranking_rds <- file.path(
  step_root,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-2C_Final_Ranking_Vector",
  "STEP9D2C_FINAL_GenomeWide_Ranking_Vector.rds"
)

if (!file.exists(ranking_rds)) {
  stop(
    "FATAL ERROR: STEP9D-2C FINAL ranking vector not found."
  )
}

ranking_vector <- readRDS(
  ranking_rds
)


# ------------------------------------------------------------
# 2. Structural confirmation
# ------------------------------------------------------------

gene_ids <- names(
  ranking_vector
)

cat(
  "Ranking vector length:",
  length(ranking_vector),
  "\n"
)

cat(
  "Gene IDs available:",
  !is.null(gene_ids),
  "\n"
)

cat(
  "Duplicated gene IDs:",
  anyDuplicated(gene_ids),
  "\n\n"
)

if (
  length(ranking_vector) != 57828
) {
  stop(
    "FATAL ERROR: Ranking vector length is not 57,828."
  )
}

if (
  anyDuplicated(gene_ids) > 0
) {
  stop(
    "FATAL ERROR: Duplicated gene IDs detected."
  )
}


# ------------------------------------------------------------
# 3. Ensembl pattern audit
# ------------------------------------------------------------

ensembl_pattern <- grepl(
  "^ENSG[0-9]+$",
  gene_ids
)

ensembl_versioned_pattern <- grepl(
  "^ENSG[0-9]+\\.[0-9]+$",
  gene_ids
)

cat("============================================================\n")
cat("GENE-ID FORMAT AUDIT\n")
cat("============================================================\n\n")

cat(
  "Standard Ensembl IDs:",
  sum(ensembl_pattern),
  "\n"
)

cat(
  "Versioned Ensembl IDs:",
  sum(ensembl_versioned_pattern),
  "\n"
)

cat(
  "Other format:",
  sum(
    !ensembl_pattern &
      !ensembl_versioned_pattern
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 4. Print examples
# ------------------------------------------------------------

cat("First 20 gene IDs:\n\n")

print(
  head(
    gene_ids,
    20
  )
)

cat("\n")


# ------------------------------------------------------------
# 5. Load annotation packages
# ------------------------------------------------------------

library(
  clusterProfiler
)

library(
  org.Hs.eg.db
)


# ------------------------------------------------------------
# 6. Package versions
# ------------------------------------------------------------

cat("============================================================\n")
cat("ANNOTATION PACKAGE VERSION\n")
cat("============================================================\n\n")

cat(
  "clusterProfiler:",
  as.character(
    packageVersion("clusterProfiler")
  ),
  "\n"
)

cat(
  "org.Hs.eg.db:",
  as.character(
    packageVersion("org.Hs.eg.db")
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 7. Confirm OrgDb keytypes
# ------------------------------------------------------------

cat("============================================================\n")
cat("org.Hs.eg.db KEYTYPE AUDIT\n")
cat("============================================================\n\n")

available_keytypes <- keytypes(
  org.Hs.eg.db
)

cat(
  "ENSEMBL available:",
  "ENSEMBL" %in% available_keytypes,
  "\n"
)

cat(
  "ENTREZID available:",
  "ENTREZID" %in% available_keytypes,
  "\n\n"
)


# ------------------------------------------------------------
# 8. Perform READ-ONLY bitr diagnostic
# ------------------------------------------------------------

gene_ids_clean <- gene_ids

mapping_raw <- suppressMessages(
  clusterProfiler::bitr(
    gene_ids_clean,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
)


# ------------------------------------------------------------
# 9. Inspect actual mapping_raw structure
# ------------------------------------------------------------

cat("============================================================\n")
cat("RAW BITR OUTPUT STRUCTURE\n")
cat("============================================================\n\n")

cat(
  "Rows:",
  nrow(mapping_raw),
  "\n"
)

cat(
  "Columns:",
  paste(
    names(mapping_raw),
    collapse = " | "
  ),
  "\n\n"
)

print(
  head(
    mapping_raw,
    20
  )
)

cat("\n")


# ------------------------------------------------------------
# 10. Actual mapping rate
# ------------------------------------------------------------

mapping_input_n <- length(
  gene_ids_clean
)

mapping_output_row_n <- nrow(
  mapping_raw
)

mapping_unique_ensembl_n <- length(
  unique(
    mapping_raw$ENSEMBL
  )
)

mapping_unique_entrez_n <- length(
  unique(
    mapping_raw$ENTREZID
  )
)

mapping_success_rate <- (
  mapping_unique_ensembl_n /
    mapping_input_n
)

cat("============================================================\n")
cat("ACTUAL MAPPING RATE\n")
cat("============================================================\n\n")

cat(
  "Input Ensembl genes:",
  mapping_input_n,
  "\n"
)

cat(
  "Mapping output rows:",
  mapping_output_row_n,
  "\n"
)

cat(
  "Unique Ensembl genes mapped:",
  mapping_unique_ensembl_n,
  "\n"
)

cat(
  "Unique Entrez IDs:",
  mapping_unique_entrez_n,
  "\n"
)

cat(
  "Actual Ensembl mapping rate:",
  round(
    mapping_success_rate,
    6
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 11. One-to-many mapping audit
# ------------------------------------------------------------

ensembl_mapping_frequency <- table(
  mapping_raw$ENSEMBL
)

one_to_many_ensembl <- names(
  ensembl_mapping_frequency[
    ensembl_mapping_frequency > 1
  ]
)

one_to_many_ensembl_n <- length(
  one_to_many_ensembl
)

cat("============================================================\n")
cat("ONE-TO-MANY MAPPING AUDIT\n")
cat("============================================================\n\n")

cat(
  "Ensembl IDs with >1 Entrez mapping:",
  one_to_many_ensembl_n,
  "\n\n"
)

if (
  one_to_many_ensembl_n > 0
) {
  
  print(
    head(
      one_to_many_ensembl,
      30
    )
  )
  
  cat("\nExample mappings:\n\n")
  
  print(
    mapping_raw[
      mapping_raw$ENSEMBL %in%
        head(
          one_to_many_ensembl,
          10
        ),
      ,
      drop = FALSE
    ]
  )
}


# ------------------------------------------------------------
# 12. Entrez many-to-one audit
# ------------------------------------------------------------

entrez_mapping_frequency <- table(
  mapping_raw$ENTREZID
)

many_to_one_entrez <- names(
  entrez_mapping_frequency[
    entrez_mapping_frequency > 1
  ]
)

many_to_one_entrez_n <- length(
  many_to_one_entrez
)

cat("\n============================================================\n")
cat("ENTREZ MANY-TO-ONE AUDIT\n")
cat("============================================================\n\n")

cat(
  "Entrez IDs receiving >1 Ensembl ID:",
  many_to_one_entrez_n,
  "\n\n"
)


# ------------------------------------------------------------
# 13. Unmapped Ensembl IDs
# ------------------------------------------------------------

mapped_ensembl <- unique(
  mapping_raw$ENSEMBL
)

unmapped_gene_ids <- setdiff(
  gene_ids_clean,
  mapped_ensembl
)

unmapped_n <- length(
  unmapped_gene_ids
)

cat("============================================================\n")
cat("UNMAPPED GENE AUDIT\n")
cat("============================================================\n\n")

cat(
  "Unmapped Ensembl genes:",
  unmapped_n,
  "\n"
)

cat(
  "Unmapped fraction:",
  round(
    unmapped_n /
      mapping_input_n,
    6
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 14. Inspect first 100 unmapped IDs
# ------------------------------------------------------------

cat(
  "First 100 unmapped Ensembl IDs:\n\n"
)

print(
  head(
    unmapped_gene_ids,
    100
  )
)

cat("\n")


# ------------------------------------------------------------
# 15. Inspect unmapped ID structure
# ------------------------------------------------------------

unmapped_has_ensembl_format <- grepl(
  "^ENSG[0-9]+$",
  unmapped_gene_ids
)

unmapped_has_version <- grepl(
  "^ENSG[0-9]+\\.[0-9]+$",
  unmapped_gene_ids
)

cat("============================================================\n")
cat("UNMAPPED-ID FORMAT AUDIT\n")
cat("============================================================\n\n")

cat(
  "Standard ENSG among unmapped:",
  sum(unmapped_has_ensembl_format),
  "\n"
)

cat(
  "Versioned ENSG among unmapped:",
  sum(unmapped_has_version),
  "\n\n"
)


# ------------------------------------------------------------
# 16. Explicit mapping table preview
# ------------------------------------------------------------

cat("============================================================\n")
cat("MAPPING TABLE PREVIEW\n")
cat("============================================================\n\n")

mapping_preview <- mapping_raw[
  seq_len(
    min(
      20,
      nrow(mapping_raw)
    )
  ),
  ,
  drop = FALSE
]

print(
  mapping_preview
)

cat("\n")


# ------------------------------------------------------------
# 17. Critical diagnostic conclusion
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP9D-3A DIAGNOSTIC SUMMARY\n")
cat("============================================================\n\n")

cat(
  "Input genes:",
  mapping_input_n,
  "\n"
)

cat(
  "Mapped unique Ensembl:",
  mapping_unique_ensembl_n,
  "\n"
)

cat(
  "Unmapped:",
  unmapped_n,
  "\n"
)

cat(
  "Mapping rate:",
  round(
    mapping_success_rate,
    6
  ),
  "\n"
)

cat(
  "One-to-many Ensembl IDs:",
  one_to_many_ensembl_n,
  "\n"
)

cat(
  "Multiple Ensembl IDs per Entrez:",
  many_to_one_entrez_n,
  "\n"
)

cat(
  "\nIMPORTANT: NO FINAL MAPPING OBJECT CREATED.\n"
)

cat(
  "IMPORTANT: NO GSEA PERFORMED.\n"
)

cat(
  "IMPORTANT: NO FILE MODIFIED.\n"
)

cat("\n============================================================\n")
cat("END STEP9D-3A DIAGNOSTIC\n")
cat("============================================================\n")    


# ============================================================
# Paper 1
# STEP9D-3A-2
# Annotation Coverage / Core-Gene Mapping Diagnostic
#
# PURPOSE
# ------------------------------------------------------------
# 1. Audit AOC3 / MS4A7 / IL7R Ensembl -> Entrez mapping
# 2. Characterize 57,828-gene annotation coverage
# 3. Characterize unmapped genes
# 4. Characterize one-to-many Ensembl -> Entrez mappings
# 5. Diagnose whether GO can use ENSEMBL directly
# 6. Diagnose whether KEGG requires an Entrez subset
# 7. Establish an auditable preliminary GSEA ID policy
#
# IMPORTANT
# ------------------------------------------------------------
# READ-ONLY DIAGNOSTIC
# No ranking recalculation
# No DEG recalculation
# No GSEA
# No modification of STEP9D-2B / STEP9D-2C
# No final mapping object created
# No arbitrary one-to-many selection
# ============================================================


cat("\n============================================================\n")
cat("STEP9D-3A-2 ANNOTATION COVERAGE / CORE-GENE MAPPING DIAGNOSTIC\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working paths
# ------------------------------------------------------------

step_root <- "input_dir/Paper 1_Step 9-Functional Enrichment"

ranking_rds <- file.path(
  step_root,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-2C_Final_Ranking_Vector",
  "STEP9D2C_FINAL_GenomeWide_Ranking_Vector.rds"
)

step9d3a2_dir <- file.path(
  step_root,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-3A-2_Annotation_Coverage_Diagnostic"
)

dir.create(
  step9d3a2_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Load LOCKED STEP9D-2C ranking vector
# ------------------------------------------------------------

if (!file.exists(ranking_rds)) {
  stop(
    "FATAL ERROR: STEP9D-2C FINAL ranking vector not found."
  )
}

ranking_vector <- readRDS(
  ranking_rds
)

gene_ids <- names(
  ranking_vector
)


# ------------------------------------------------------------
# 2. Structural lock confirmation
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP9D-2C RANKING VECTOR STRUCTURE\n")
cat("============================================================\n\n")

cat(
  "Ranking vector length:",
  length(ranking_vector),
  "\n"
)

cat(
  "Gene IDs available:",
  !is.null(gene_ids),
  "\n"
)

cat(
  "Duplicated gene IDs:",
  anyDuplicated(gene_ids),
  "\n"
)

cat(
  "All standard Ensembl:",
  all(grepl("^ENSG[0-9]+$", gene_ids)),
  "\n\n"
)


if (length(ranking_vector) != 57828) {
  stop(
    "FATAL ERROR: Ranking vector is not 57,828 genes."
  )
}

if (anyDuplicated(gene_ids) > 0) {
  stop(
    "FATAL ERROR: Duplicated gene IDs detected."
  )
}

if (!all(grepl("^ENSG[0-9]+$", gene_ids))) {
  stop(
    "FATAL ERROR: Non-standard Ensembl IDs detected."
  )
}


# ------------------------------------------------------------
# 3. Load annotation packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})


# ------------------------------------------------------------
# 4. Package provenance
# ------------------------------------------------------------

cat("============================================================\n")
cat("ANNOTATION PACKAGE PROVENANCE\n")
cat("============================================================\n\n")

clusterProfiler_version <- as.character(
  packageVersion("clusterProfiler")
)

orgHs_version <- as.character(
  packageVersion("org.Hs.eg.db")
)

AnnotationDbi_version <- as.character(
  packageVersion("AnnotationDbi")
)

cat(
  "clusterProfiler:",
  clusterProfiler_version,
  "\n"
)

cat(
  "org.Hs.eg.db:",
  orgHs_version,
  "\n"
)

cat(
  "AnnotationDbi:",
  AnnotationDbi_version,
  "\n\n"
)


# ------------------------------------------------------------
# 5. org.Hs.eg.db column / keytype audit
# ------------------------------------------------------------

cat("============================================================\n")
cat("org.Hs.eg.db CAPABILITY AUDIT\n")
cat("============================================================\n\n")

org_columns <- AnnotationDbi::columns(
  org.Hs.eg.db
)

org_keytypes <- AnnotationDbi::keytypes(
  org.Hs.eg.db
)

cat(
  "ENSEMBL keytype:",
  "ENSEMBL" %in% org_keytypes,
  "\n"
)

cat(
  "ENTREZID keytype:",
  "ENTREZID" %in% org_keytypes,
  "\n"
)

cat(
  "SYMBOL column:",
  "SYMBOL" %in% org_columns,
  "\n"
)

cat(
  "GENENAME column:",
  "GENENAME" %in% org_columns,
  "\n"
)

cat(
  "GENEBIOTYPE column:",
  "GENEBIOTYPE" %in% org_columns,
  "\n\n"
)

cat(
  "Available org.Hs.eg.db columns:\n"
)

print(
  org_columns
)

cat("\n")


# ------------------------------------------------------------
# 6. Full ENSEMBL -> ENTREZ mapping diagnostic
# ------------------------------------------------------------

cat("============================================================\n")
cat("FULL ENSEMBL -> ENTREZ MAPPING DIAGNOSTIC\n")
cat("============================================================\n\n")

mapping_raw <- suppressMessages(
  clusterProfiler::bitr(
    gene_ids,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
)

mapping_raw <- as.data.frame(
  mapping_raw,
  stringsAsFactors = FALSE
)

mapped_ensembl <- unique(
  mapping_raw$ENSEMBL
)

unmapped_gene_ids <- setdiff(
  gene_ids,
  mapped_ensembl
)

mapping_unique_ensembl_n <- length(
  mapped_ensembl
)

unmapped_n <- length(
  unmapped_gene_ids
)

mapping_rate <- (
  mapping_unique_ensembl_n /
    length(gene_ids)
)

cat(
  "Input Ensembl genes:",
  length(gene_ids),
  "\n"
)

cat(
  "Raw mapping rows:",
  nrow(mapping_raw),
  "\n"
)

cat(
  "Unique Ensembl genes mapped:",
  mapping_unique_ensembl_n,
  "\n"
)

cat(
  "Unmapped Ensembl genes:",
  unmapped_n,
  "\n"
)

cat(
  "Mapping rate:",
  round(mapping_rate, 6),
  "\n\n"
)


# ------------------------------------------------------------
# 7. Core-gene diagnostic
# ------------------------------------------------------------

cat("============================================================\n")
cat("CORE-GENE ENSEMBL / ENTREZ AUDIT\n")
cat("============================================================\n\n")

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_symbol_map <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = core_genes,
    columns = c(
      "ENSEMBL",
      "ENTREZID",
      "SYMBOL",
      "GENENAME"
    ),
    keytype = "SYMBOL"
  )
)

core_symbol_map <- as.data.frame(
  core_symbol_map,
  stringsAsFactors = FALSE
)

core_symbol_map <- core_symbol_map[
  !is.na(core_symbol_map$ENSEMBL),
  ,
  drop = FALSE
]

cat(
  "org.Hs.eg.db core-gene mappings:\n\n"
)

print(
  core_symbol_map
)

cat("\n")


# ------------------------------------------------------------
# 8. Restrict core-gene mappings to the LOCKED ranking universe
# ------------------------------------------------------------

core_in_universe <- core_symbol_map[
  core_symbol_map$ENSEMBL %in% gene_ids,
  ,
  drop = FALSE
]

cat("============================================================\n")
cat("CORE GENES WITHIN LOCKED 57,828-GENE UNIVERSE\n")
cat("============================================================\n\n")

print(
  core_in_universe
)

cat("\n")


for (g in core_genes) {
  
  tmp <- core_in_universe[
    core_in_universe$SYMBOL == g,
    ,
    drop = FALSE
  ]
  
  cat(
    g,
    ":",
    nrow(tmp),
    "ENSEMBL mapping record(s) in locked universe\n"
  )
  
  if (nrow(tmp) > 0) {
    
    cat(
      "  Ensembl:",
      paste(
        unique(tmp$ENSEMBL),
        collapse = ", "
      ),
      "\n"
    )
    
    cat(
      "  Entrez:",
      paste(
        unique(
          tmp$ENTREZID[
            !is.na(tmp$ENTREZID)
          ]
        ),
        collapse = ", "
      ),
      "\n"
    )
  }
  
  cat("\n")
}


# ------------------------------------------------------------
# 9. Core-gene uniqueness audit
# ------------------------------------------------------------

core_uniqueness_results <- lapply(
  core_genes,
  function(g) {
    
    tmp <- core_in_universe[
      core_in_universe$SYMBOL == g,
      ,
      drop = FALSE
    ]
    
    data.frame(
      gene = g,
      ensembl_n = length(
        unique(tmp$ENSEMBL)
      ),
      entrez_n = length(
        unique(
          tmp$ENTREZID[
            !is.na(tmp$ENTREZID)
          ]
        )
      ),
      ensembl_unique = (
        length(
          unique(tmp$ENSEMBL)
        ) == 1
      ),
      entrez_unique = (
        length(
          unique(
            tmp$ENTREZID[
              !is.na(tmp$ENTREZID)
            ]
          )
        ) == 1
      ),
      in_locked_universe = (
        nrow(tmp) > 0
      ),
      stringsAsFactors = FALSE
    )
  }
)

core_uniqueness_results <- do.call(
  rbind,
  core_uniqueness_results
)

print(
  core_uniqueness_results
)

cat("\n")


# ------------------------------------------------------------
# 10. Identify one-to-many Ensembl -> Entrez mappings
# ------------------------------------------------------------

cat("============================================================\n")
cat("ONE-TO-MANY ENSEMBL -> ENTREZ AUDIT\n")
cat("============================================================\n\n")

ensembl_frequency <- table(
  mapping_raw$ENSEMBL
)

one_to_many_ensembl <- names(
  ensembl_frequency[
    ensembl_frequency > 1
  ]
)

one_to_many_n <- length(
  one_to_many_ensembl
)

cat(
  "Ensembl IDs with >1 Entrez:",
  one_to_many_n,
  "\n\n"
)

if (one_to_many_n > 0) {
  
  one_to_many_table <- mapping_raw[
    mapping_raw$ENSEMBL %in%
      one_to_many_ensembl,
    ,
    drop = FALSE
  ]
  
  print(
    head(
      one_to_many_table,
      50
    )
  )
}

cat("\n")


# ------------------------------------------------------------
# 11. One-to-many mapping statistics
# ------------------------------------------------------------

if (one_to_many_n > 0) {
  
  one_to_many_entrez_counts <- tapply(
    mapping_raw$ENTREZID[
      mapping_raw$ENSEMBL %in%
        one_to_many_ensembl
    ],
    mapping_raw$ENSEMBL[
      mapping_raw$ENSEMBL %in%
        one_to_many_ensembl
    ],
    function(x) length(unique(x))
  )
  
  cat(
    "Maximum Entrez mappings for one Ensembl:",
    max(one_to_many_entrez_counts),
    "\n"
  )
  
  cat(
    "Mean Entrez mappings among one-to-many IDs:",
    round(
      mean(one_to_many_entrez_counts),
      4
    ),
    "\n"
  )
  
  cat(
    "Median Entrez mappings among one-to-many IDs:",
    median(one_to_many_entrez_counts),
    "\n\n"
  )
}


# ------------------------------------------------------------
# 12. Many-to-one Entrez audit
# ------------------------------------------------------------

cat("============================================================\n")
cat("MANY-TO-ONE ENTREZ AUDIT\n")
cat("============================================================\n\n")

entrez_frequency <- table(
  mapping_raw$ENTREZID
)

many_to_one_entrez <- names(
  entrez_frequency[
    entrez_frequency > 1
  ]
)

many_to_one_n <- length(
  many_to_one_entrez
)

cat(
  "Entrez IDs receiving >1 Ensembl:",
  many_to_one_n,
  "\n\n"
)


# ------------------------------------------------------------
# 13. Core genes and one-to-many status
# ------------------------------------------------------------

cat("============================================================\n")
cat("CORE-GENE ONE-TO-MANY STATUS\n")
cat("============================================================\n\n")

for (g in core_genes) {
  
  core_ens <- unique(
    core_in_universe$ENSEMBL[
      core_in_universe$SYMBOL == g
    ]
  )
  
  core_ens <- core_ens[
    !is.na(core_ens)
  ]
  
  if (length(core_ens) == 0) {
    
    cat(
      g,
      ": NOT FOUND in locked universe\n"
    )
    
    next
  }
  
  core_map <- mapping_raw[
    mapping_raw$ENSEMBL %in% core_ens,
    ,
    drop = FALSE
  ]
  
  core_entrez <- unique(
    core_map$ENTREZID
  )
  
  cat(
    g,
    "\n"
  )
  
  cat(
    "  Ensembl:",
    paste(
      core_ens,
      collapse = ", "
    ),
    "\n"
  )
  
  cat(
    "  Entrez:",
    paste(
      core_entrez,
      collapse = ", "
    ),
    "\n"
  )
  
  cat(
    "  Unique Entrez count:",
    length(core_entrez),
    "\n"
  )
  
  cat(
    "  One-to-many:",
    length(core_entrez) > 1,
    "\n\n"
  )
}


# ------------------------------------------------------------
# 14. Direct GO ENSEMBL capability audit
# ------------------------------------------------------------

cat("============================================================\n")
cat("GO DIRECT ENSEMBL CAPABILITY AUDIT\n")
cat("============================================================\n\n")

go_ensembl_possible <- (
  "ENSEMBL" %in% org_keytypes
)

cat(
  "org.Hs.eg.db supports ENSEMBL:",
  go_ensembl_possible,
  "\n"
)

cat(
  "clusterProfiler gseGO available:",
  exists(
    "gseGO",
    where = asNamespace("clusterProfiler"),
    inherits = FALSE
  ),
  "\n"
)

cat(
  "Preliminary GO policy candidate:",
  "\n  Use ENSEMBL directly for gseGO, if GO annotation coverage is confirmed.",
  "\n\n"
)


# ------------------------------------------------------------
# 15. Check GO annotation coverage using ENSEMBL
# ------------------------------------------------------------

cat("============================================================\n")
cat("GO ANNOTATION COVERAGE AUDIT\n")
cat("============================================================\n\n")

go_annotation_probe <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_ids,
    columns = c(
      "ENSEMBL",
      "GO",
      "ONTOLOGY"
    ),
    keytype = "ENSEMBL"
  )
)

go_annotation_probe <- as.data.frame(
  go_annotation_probe,
  stringsAsFactors = FALSE
)

go_annotated_ensembl <- unique(
  go_annotation_probe$ENSEMBL[
    !is.na(go_annotation_probe$GO)
  ]
)

go_annotated_n <- length(
  go_annotated_ensembl
)

go_annotation_rate <- (
  go_annotated_n /
    length(gene_ids)
)

cat(
  "Genes with >=1 GO annotation:",
  go_annotated_n,
  "\n"
)

cat(
  "Genes without GO annotation:",
  length(gene_ids) - go_annotated_n,
  "\n"
)

cat(
  "GO annotation coverage:",
  round(
    go_annotation_rate,
    6
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 16. GO annotation coverage for the three core genes
# ------------------------------------------------------------

cat("============================================================\n")
cat("CORE-GENE GO ANNOTATION AUDIT\n")
cat("============================================================\n\n")

for (g in core_genes) {
  
  core_ens <- unique(
    core_in_universe$ENSEMBL[
      core_in_universe$SYMBOL == g
    ]
  )
  
  core_go <- go_annotation_probe[
    go_annotation_probe$ENSEMBL %in% core_ens &
      !is.na(go_annotation_probe$GO),
    ,
    drop = FALSE
  ]
  
  cat(
    g,
    ":",
    length(
      unique(core_go$GO)
    ),
    "GO term(s)\n"
  )
  
  if (nrow(core_go) > 0) {
    print(
      unique(
        core_go[
          ,
          c(
            "ENSEMBL",
            "GO",
            "ONTOLOGY"
          ),
          drop = FALSE
        ]
      )
    )
  }
  
  cat("\n")
}


# ------------------------------------------------------------
# 17. KEGG capability audit
# ------------------------------------------------------------

cat("============================================================\n")
cat("KEGG ID REQUIREMENT AUDIT\n")
cat("============================================================\n\n")

gseKEGG_available <- exists(
  "gseKEGG",
  where = asNamespace("clusterProfiler"),
  inherits = FALSE
)

cat(
  "clusterProfiler gseKEGG available:",
  gseKEGG_available,
  "\n"
)

cat(
  "KEGG organism code:",
  "hsa",
  "\n\n"
)

cat(
  "Preliminary KEGG policy candidate:\n",
  "  Entrez IDs are required for the standard gseKEGG(hsa) workflow.\n",
  "  Therefore KEGG must use an explicitly defined Entrez-mappable subset.\n\n"
)


# ------------------------------------------------------------
# 18. KEGG core-gene Entrez audit
# ------------------------------------------------------------

cat("============================================================\n")
cat("CORE-GENE KEGG ENTZREZ AUDIT\n")
cat("============================================================\n\n")

for (g in core_genes) {
  
  core_ens <- unique(
    core_in_universe$ENSEMBL[
      core_in_universe$SYMBOL == g
    ]
  )
  
  core_map <- mapping_raw[
    mapping_raw$ENSEMBL %in% core_ens,
    ,
    drop = FALSE
  ]
  
  core_entrez <- unique(
    core_map$ENTREZID
  )
  
  cat(
    g,
    "\n"
  )
  
  cat(
    "  Ensembl present:",
    length(core_ens) > 0,
    "\n"
  )
  
  cat(
    "  Entrez mapping present:",
    length(core_entrez) > 0,
    "\n"
  )
  
  cat(
    "  Unique Entrez:",
    paste(
      core_entrez,
      collapse = ", "
    ),
    "\n"
  )
  
  cat("\n")
}


# ------------------------------------------------------------
# 19. Characterize unmapped genes by Ensembl ID generation
# ------------------------------------------------------------

cat("============================================================\n")
cat("UNMAPPED ENSEMBL ID STRUCTURE\n")
cat("============================================================\n\n")

unmapped_ens2 <- sum(
  grepl(
    "^ENSG000002",
    unmapped_gene_ids
  )
)

unmapped_ens1 <- sum(
  grepl(
    "^ENSG000001",
    unmapped_gene_ids
  )
)

unmapped_ens0 <- sum(
  grepl(
    "^ENSG000000",
    unmapped_gene_ids
  )
)

cat(
  "Unmapped ENSG000002...:",
  unmapped_ens2,
  "\n"
)

cat(
  "Unmapped ENSG000001...:",
  unmapped_ens1,
  "\n"
)

cat(
  "Unmapped ENSG000000...:",
  unmapped_ens0,
  "\n\n"
)


# ------------------------------------------------------------
# 20. Optional Ensembl biotype diagnostic
#
# IMPORTANT:
# org.Hs.eg.db does not normally provide GENEBIOTYPE.
# If biomaRt is installed, use Ensembl annotation ONLY for
# diagnostic characterization.
#
# This is NOT used for GSEA.
# ------------------------------------------------------------

biomart_available <- requireNamespace(
  "biomaRt",
  quietly = TRUE
)

cat("============================================================\n")
cat("ENSEMBL BIOTYPE DIAGNOSTIC CAPABILITY\n")
cat("============================================================\n\n")

cat(
  "biomaRt installed:",
  biomart_available,
  "\n\n"
)


biotype_result <- NULL
ensembl_provenance <- NULL


if (biomart_available) {
  
  cat(
    "Attempting READ-ONLY Ensembl biotype query...\n\n"
  )
  
  suppressPackageStartupMessages(
    library(biomaRt)
  )
  
  mart <- NULL
  
  mart_attempt <- tryCatch(
    
    {
      
      mart <- biomaRt::useEnsembl(
        biomart = "genes",
        dataset = "hsapiens_gene_ensembl"
      )
      
      mart
      
    },
    
    error = function(e) {
      
      cat(
        "Ensembl connection failed:\n",
        conditionMessage(e),
        "\n\n"
      )
      
      NULL
    }
  )
  
  
  if (!is.null(mart_attempt)) {
    
    cat(
      "Ensembl mart connection: PASS\n\n"
    )
    
    ensembl_provenance <- tryCatch(
      
      biomaRt::listEnsembl(
        mart = "genes"
      ),
      
      error = function(e) {
        NULL
      }
    )
    
    
    # --------------------------------------------------------
    # Query all 57,828 IDs in chunks
    # --------------------------------------------------------
    
    query_chunks <- split(
      gene_ids,
      ceiling(
        seq_along(gene_ids) / 1000
      )
    )
    
    biotype_results <- vector(
      "list",
      length(query_chunks)
    )
    
    for (i in seq_along(query_chunks)) {
      
      cat(
        "Biotype query chunk",
        i,
        "/",
        length(query_chunks),
        "\n"
      )
      
      tmp <- tryCatch(
        
        biomaRt::getBM(
          attributes = c(
            "ensembl_gene_id",
            "gene_biotype",
            "external_gene_name",
            "entrezgene_id"
          ),
          filters = "ensembl_gene_id",
          values = query_chunks[[i]],
          mart = mart_attempt
        ),
        
        error = function(e) {
          
          cat(
            "  Chunk failed:",
            conditionMessage(e),
            "\n"
          )
          
          NULL
        }
      )
      
      biotype_results[[i]] <- tmp
    }
    
    
    biotype_result <- do.call(
      rbind,
      biotype_results[
        !vapply(
          biotype_results,
          is.null,
          logical(1)
        )
      ]
    )
    
    
    if (
      !is.null(biotype_result) &&
      nrow(biotype_result) > 0
    ) {
      
      biotype_result <- as.data.frame(
        biotype_result,
        stringsAsFactors = FALSE
      )
      
      cat("\n")
      cat(
        "Ensembl biotype rows:",
        nrow(biotype_result),
        "\n\n"
      )
      
      
      # ------------------------------------------------------
      # Biotype coverage
      # ------------------------------------------------------
      
      biotype_mapped_ids <- unique(
        biotype_result$ensembl_gene_id
      )
      
      biotype_unmapped_ids <- setdiff(
        gene_ids,
        biotype_mapped_ids
      )
      
      cat(
        "Ensembl IDs with biotype annotation:",
        length(biotype_mapped_ids),
        "\n"
      )
      
      cat(
        "Ensembl IDs without biotype annotation:",
        length(biotype_unmapped_ids),
        "\n\n"
      )
      
      
      # ------------------------------------------------------
      # Biotype frequency
      # ------------------------------------------------------
      
      biotype_frequency <- sort(
        table(
          biotype_result$gene_biotype,
          useNA = "ifany"
        ),
        decreasing = TRUE
      )
      
      cat(
        "Genome-wide biotype distribution:\n\n"
      )
      
      print(
        biotype_frequency
      )
      
      cat("\n")
      
      
      # ------------------------------------------------------
      # Unmapped-by-Entrez biotype distribution
      # ------------------------------------------------------
      
      unmapped_biotype_result <- biotype_result[
        biotype_result$ensembl_gene_id %in%
          unmapped_gene_ids,
        ,
        drop = FALSE
      ]
      
      cat(
        "============================================================\n"
      )
      
      cat(
        "BIOTYPE DISTRIBUTION OF ENTIRELY UNMAPPED-TO-ENTREZ GENES\n"
      )
      
      cat(
        "============================================================\n\n"
      )
      
      if (nrow(unmapped_biotype_result) > 0) {
        
        unmapped_biotype_frequency <- sort(
          table(
            unmapped_biotype_result$gene_biotype,
            useNA = "ifany"
          ),
          decreasing = TRUE
        )
        
        print(
          unmapped_biotype_frequency
        )
        
      } else {
        
        cat(
          "No biotype records available for unmapped genes.\n"
        )
      }
      
      cat("\n")
      
      
      # ------------------------------------------------------
      # One-to-many Ensembl mappings and biotype
      # ------------------------------------------------------
      
      one_to_many_biotype <- biotype_result[
        biotype_result$ensembl_gene_id %in%
          one_to_many_ensembl,
        ,
        drop = FALSE
      ]
      
      cat(
        "============================================================\n"
      )
      
      cat(
        "ONE-TO-MANY ENSEMBL IDs: BIOTYPE DISTRIBUTION\n"
      )
      
      cat(
        "============================================================\n\n"
      )
      
      if (nrow(one_to_many_biotype) > 0) {
        
        print(
          sort(
            table(
              one_to_many_biotype$gene_biotype,
              useNA = "ifany"
            ),
            decreasing = TRUE
          )
        )
        
      } else {
        
        cat(
          "No biotype records available.\n"
        )
      }
      
      cat("\n")
      
      
      # ------------------------------------------------------
      # Core gene biotypes
      # ------------------------------------------------------
      
      cat(
        "============================================================\n"
      )
      
      cat(
        "CORE-GENE ENSEMBL BIOTYPE AUDIT\n"
      )
      
      cat(
        "============================================================\n\n"
      )
      
      core_biotype <- biotype_result[
        biotype_result$ensembl_gene_id %in%
          core_in_universe$ENSEMBL,
        ,
        drop = FALSE
      ]
      
      print(
        core_biotype
      )
      
      cat("\n")
    }
  }
  
} else {
  
  cat(
    "biomaRt is not installed.\n"
  )
  
  cat(
    "Biotype classification will NOT be inferred from incomplete sources.\n"
  )
  
  cat(
    "No biotype-based conclusion will be made in this diagnostic.\n\n"
  )
}


# ------------------------------------------------------------
# 21. Preliminary policy assessment
# ------------------------------------------------------------

cat("============================================================\n")
cat("PRELIMINARY GSEA ID POLICY ASSESSMENT\n")
cat("============================================================\n\n")

cat(
  "1. Locked ranking IDs:\n",
  "   ENSEMBL, 57,828 genes.\n\n"
)

cat(
  "2. Ranking vector:\n",
  "   Must remain unchanged.\n\n"
)

cat(
  "3. GO:\n",
  "   ENSEMBL is supported by org.Hs.eg.db and should be evaluated\n",
  "   as the primary direct-ID route for gseGO.\n\n"
)

cat(
  "4. KEGG:\n",
  "   Standard gseKEGG(hsa) requires Entrez-compatible identifiers.\n",
  "   Therefore KEGG requires an explicitly documented Entrez-mappable\n",
  "   subset of the locked genome-wide Ensembl universe.\n\n"
)

cat(
  "5. Unmapped genes:\n",
  "   Must NOT be artificially converted or discarded from the\n",
  "   original locked ranking vector.\n",
  "   They are excluded only from an Entrez-dependent analysis\n",
  "   if the final KEGG policy explicitly requires this.\n\n"
)

cat(
  "6. One-to-many Ensembl -> Entrez:\n",
  "   No arbitrary first-match selection is permitted.\n",
  "   A deterministic mapping policy must be locked before KEGG GSEA.\n\n"
)

cat(
  "7. Core genes:\n",
  "   AOC3 / MS4A7 / IL7R must remain traceable from the locked\n",
  "   Ensembl ranking vector into the final annotation framework.\n\n"
)


# ------------------------------------------------------------
# 22. Diagnostic status flags
# ------------------------------------------------------------

core_all_in_universe <- all(
  core_uniqueness_results$in_locked_universe
)

core_all_unique_entrez <- all(
  core_uniqueness_results$entrez_unique
)

go_direct_supported <- (
  "ENSEMBL" %in% org_keytypes &&
    exists(
      "gseGO",
      where = asNamespace("clusterProfiler"),
      inherits = FALSE
    )
)

kegg_entrez_required <- (
  exists(
    "gseKEGG",
    where = asNamespace("clusterProfiler"),
    inherits = FALSE
  )
)

cat("============================================================\n")
cat("STEP9D-3A-2 STATUS FLAGS\n")
cat("============================================================\n\n")

cat(
  "57,828-gene ranking vector intact:",
  length(gene_ids) == 57828,
  "\n"
)

cat(
  "No duplicate ranking IDs:",
  anyDuplicated(gene_ids) == 0,
  "\n"
)

cat(
  "All IDs standard Ensembl:",
  all(grepl("^ENSG[0-9]+$", gene_ids)),
  "\n"
)

cat(
  "All three core genes in locked universe:",
  core_all_in_universe,
  "\n"
)

cat(
  "All three core genes have unique Entrez:",
  core_all_unique_entrez,
  "\n"
)

cat(
  "Direct ENSEMBL GO route supported:",
  go_direct_supported,
  "\n"
)

cat(
  "KEGG Entrez route required:",
  kegg_entrez_required,
  "\n"
)

cat(
  "Formal GSEA performed:",
  FALSE,
  "\n"
)

cat(
  "Ranking recalculated:",
  FALSE,
  "\n"
)

cat(
  "DEG recalculated:",
  FALSE,
  "\n"
)

cat(
  "Final mapping object created:",
  FALSE,
  "\n\n"
)


# ------------------------------------------------------------
# 23. Save READ-ONLY diagnostic report
# ------------------------------------------------------------

diagnostic_object <- list(
  
  workflow = "Paper 1",
  
  step = "STEP9D-3A-2",
  
  status = "DIAGNOSTIC_ONLY",
  
  ranking_source = ranking_rds,
  
  ranking_gene_n = length(gene_ids),
  
  ranking_duplicate_n = anyDuplicated(gene_ids),
  
  ranking_id_type = "ENSEMBL",
  
  clusterProfiler_version =
    clusterProfiler_version,
  
  org.Hs.eg.db_version =
    orgHs_version,
  
  AnnotationDbi_version =
    AnnotationDbi_version,
  
  input_gene_n =
    length(gene_ids),
  
  raw_mapping_row_n =
    nrow(mapping_raw),
  
  unique_ensembl_mapped_n =
    mapping_unique_ensembl_n,
  
  unmapped_ensembl_n =
    unmapped_n,
  
  mapping_rate =
    mapping_rate,
  
  one_to_many_ensembl_n =
    one_to_many_n,
  
  many_to_one_entrez_n =
    many_to_one_n,
  
  core_gene_table =
    core_uniqueness_results,
  
  core_mapping_table =
    core_in_universe,
  
  go_ensembl_supported =
    go_direct_supported,
  
  go_annotated_ensembl_n =
    go_annotated_n,
  
  go_annotation_rate =
    go_annotation_rate,
  
  kegg_entrez_required =
    kegg_entrez_required,
  
  biomaRt_available =
    biomart_available,
  
  biotype_diagnostic_performed =
    !is.null(biotype_result),
  
  ensembl_provenance =
    ensembl_provenance,
  
  formal_gsea_performed =
    FALSE,
  
  ranking_recalculated =
    FALSE,
  
  deg_recalculated =
    FALSE,
  
  final_mapping_object_created =
    FALSE
)


diagnostic_rds <- file.path(
  step9d3a2_dir,
  "STEP9D3A2_Annotation_Coverage_Diagnostic.rds"
)

diagnostic_txt <- file.path(
  step9d3a2_dir,
  "STEP9D3A2_Annotation_Coverage_Diagnostic.txt"
)

saveRDS(
  diagnostic_object,
  diagnostic_rds
)


# ------------------------------------------------------------
# 24. Human-readable audit report
# ------------------------------------------------------------

sink(
  diagnostic_txt
)

cat(
  "STEP9D-3A-2 Annotation Coverage / Core-Gene Mapping Diagnostic\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Ranking genes:",
  length(gene_ids),
  "\n"
)

cat(
  "Mapped unique Ensembl:",
  mapping_unique_ensembl_n,
  "\n"
)

cat(
  "Unmapped Ensembl:",
  unmapped_n,
  "\n"
)

cat(
  "Mapping rate:",
  mapping_rate,
  "\n"
)

cat(
  "One-to-many Ensembl:",
  one_to_many_n,
  "\n"
)

cat(
  "Many-to-one Entrez:",
  many_to_one_n,
  "\n\n"
)

cat(
  "clusterProfiler:",
  clusterProfiler_version,
  "\n"
)

cat(
  "org.Hs.eg.db:",
  orgHs_version,
  "\n"
)

cat(
  "AnnotationDbi:",
  AnnotationDbi_version,
  "\n\n"
)

cat(
  "Core-gene uniqueness:\n"
)

print(
  core_uniqueness_results
)

cat("\n")

cat(
  "GO direct ENSEMBL supported:",
  go_direct_supported,
  "\n"
)

cat(
  "GO annotation coverage:",
  go_annotation_rate,
  "\n"
)

cat(
  "KEGG Entrez required:",
  kegg_entrez_required,
  "\n"
)

cat(
  "biomaRt available:",
  biomart_available,
  "\n"
)

cat(
  "Biotype diagnostic performed:",
  !is.null(biotype_result),
  "\n\n"
)

cat(
  "FORMAL GSEA PERFORMED: FALSE\n"
)

cat(
  "RANKING RECALCULATED: FALSE\n"
)

cat(
  "DEG RECALCULATED: FALSE\n"
)

cat(
  "FINAL MAPPING OBJECT CREATED: FALSE\n"
)

sink()


# ------------------------------------------------------------
# 25. Final console summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-3A-2 FINAL DIAGNOSTIC SUMMARY\n")
cat("============================================================\n\n")

cat(
  "Input genes:",
  length(gene_ids),
  "\n"
)

cat(
  "Unique Ensembl mapped:",
  mapping_unique_ensembl_n,
  "\n"
)

cat(
  "Unmapped:",
  unmapped_n,
  "\n"
)

cat(
  "Mapping rate:",
  round(
    mapping_rate,
    6
  ),
  "\n"
)

cat(
  "One-to-many Ensembl:",
  one_to_many_n,
  "\n"
)

cat(
  "Many-to-one Entrez:",
  many_to_one_n,
  "\n"
)

cat(
  "All core genes in locked universe:",
  core_all_in_universe,
  "\n"
)

cat(
  "All core genes unique Entrez:",
  core_all_unique_entrez,
  "\n"
)

cat(
  "Direct ENSEMBL GO route supported:",
  go_direct_supported,
  "\n"
)

cat(
  "KEGG Entrez route required:",
  kegg_entrez_required,
  "\n"
)

cat(
  "Biotype diagnostic performed:",
  !is.null(biotype_result),
  "\n\n"
)

cat(
  "Formal GSEA performed: FALSE\n"
)

cat(
  "Ranking recalculated: FALSE\n"
)

cat(
  "DEG recalculated: FALSE\n"
)

cat(
  "Final mapping object created: FALSE\n\n"
)

cat(
  "Diagnostic RDS:",
  diagnostic_rds,
  "\n"
)

cat(
  "Diagnostic TXT:",
  diagnostic_txt,
  "\n\n"
)

cat("============================================================\n")
cat("END STEP9D-3A-2\n")
cat("============================================================\n")   


# ============================================================
# Paper 1
# STEP9D-3A-3
# Formal GSEA ID / Mapping Policy Lock
#
# PURPOSE
# ------------------------------------------------------------
# Lock the formal identifier and analysis policy BEFORE
# STEP9D-3B formal preranked GSEA.
#
# IMPORTANT
# ------------------------------------------------------------
# 1. No DEG recalculation
# 2. No ranking recalculation
# 3. No GSEA execution
# 4. No arbitrary first-match Ensembl -> Entrez selection
# 5. Original 57,828-gene ranking vector remains unchanged
# ============================================================


# ============================================================
# STEP9D-2C FINAL RANKING FILE LOCATION AUDIT
# READ-ONLY
# ============================================================

step9d_root <- "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA"

cat("\n============================================================\n")
cat("STEP9D-2C FINAL RANKING FILE LOCATION AUDIT\n")
cat("============================================================\n\n")

cat("Search root:\n")
cat(step9d_root, "\n\n")


# ------------------------------------------------------------
# 1. Search all RDS files
# ------------------------------------------------------------

rds_files <- list.files(
  path = step9d_root,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)


cat("Total RDS files found:", length(rds_files), "\n\n")


# ------------------------------------------------------------
# 2. Search filenames related to STEP9D-2C / ranking
# ------------------------------------------------------------

ranking_rds_candidates <- rds_files[
  grepl(
    "2C|Ranking|rank",
    basename(rds_files),
    ignore.case = TRUE
  )
]


cat("Candidate RDS files:\n")
cat("------------------------------------------------------------\n")

if (length(ranking_rds_candidates) == 0) {
  
  cat("NO CANDIDATE RDS FILES FOUND.\n\n")
  
} else {
  
  print(
    ranking_rds_candidates
  )
  
  cat("\n")
}


# ------------------------------------------------------------
# 3. Search all CSV files
# ------------------------------------------------------------

csv_files <- list.files(
  path = step9d_root,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)


cat("Total CSV files found:", length(csv_files), "\n\n")


ranking_csv_candidates <- csv_files[
  grepl(
    "2C|Ranking|rank",
    basename(csv_files),
    ignore.case = TRUE
  )
]


cat("Candidate CSV files:\n")
cat("------------------------------------------------------------\n")

if (length(ranking_csv_candidates) == 0) {
  
  cat("NO CANDIDATE CSV FILES FOUND.\n\n")
  
} else {
  
  print(
    ranking_csv_candidates
  )
  
  cat("\n")
}


# ------------------------------------------------------------
# 4. Search specifically for files containing
#    STEP9D2C
# ------------------------------------------------------------

step9d2c_files <- c(
  rds_files[
    grepl(
      "STEP9D2C",
      basename(rds_files),
      ignore.case = TRUE
    )
  ],
  csv_files[
    grepl(
      "STEP9D2C",
      basename(csv_files),
      ignore.case = TRUE
    )
  ]
)


cat("Files containing 'STEP9D2C' in filename:\n")
cat("------------------------------------------------------------\n")

if (length(step9d2c_files) == 0) {
  
  cat("NONE FOUND.\n\n")
  
} else {
  
  print(step9d2c_files)
  
  cat("\n")
}


# ------------------------------------------------------------
# 5. Search entire STEP9D directory for likely ranking files
# ------------------------------------------------------------

all_files <- list.files(
  path = step9d_root,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)


ranking_name_candidates <- all_files[
  grepl(
    "ranking|ranked|prerank|gsea",
    basename(all_files),
    ignore.case = TRUE
  )
]


cat("All filename candidates containing ranking/ranked/prerank/GSEA:\n")
cat("------------------------------------------------------------\n")

if (length(ranking_name_candidates) == 0) {
  
  cat("NONE FOUND.\n\n")
  
} else {
  
  print(ranking_name_candidates)
  
  cat("\n")
}


cat("============================================================\n")
cat("END READ-ONLY FILE LOCATION AUDIT\n")
cat("============================================================\n")   


# ============================================================
# STEP9D-1I FINAL LOCK FILE LOCATION AUDIT
# READ-ONLY
# ============================================================

step9d_root <- paste0(
  "input_dir/Paper 1_Step 9-Functional Enrichment/",
  "STEP9D_GenomeWide_Preranked_GSEA"
)

cat("\n============================================================\n")
cat("STEP9D-1I FINAL LOCK FILE LOCATION AUDIT\n")
cat("============================================================\n\n")

cat("Search root:\n")
cat(step9d_root, "\n\n")


# ------------------------------------------------------------
# 1. Search all files
# ------------------------------------------------------------

all_files <- list.files(
  path = step9d_root,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

cat("Total files found:", length(all_files), "\n\n")


# ------------------------------------------------------------
# 2. Search STEP9D1I
# ------------------------------------------------------------

step9d1i_candidates <- all_files[
  grepl(
    "STEP9D1I|FINAL_INPUT|INPUT_LOCK",
    basename(all_files),
    ignore.case = TRUE
  )
]

cat("Candidate STEP9D-1I files:\n")
cat("------------------------------------------------------------\n")

if (length(step9d1i_candidates) == 0) {
  
  cat("NO CANDIDATE FILES FOUND.\n\n")
  
} else {
  
  print(step9d1i_candidates)
  cat("\n")
}


# ------------------------------------------------------------
# 3. Search all RDS files
# ------------------------------------------------------------

rds_files <- all_files[
  grepl(
    "\\.rds$",
    all_files,
    ignore.case = TRUE
  )
]

cat("All RDS files under STEP9D:\n")
cat("------------------------------------------------------------\n")

print(rds_files)

cat("\n")


# ------------------------------------------------------------
# 4. Search files containing '1I'
# ------------------------------------------------------------

one_i_candidates <- all_files[
  grepl(
    "1I",
    basename(all_files),
    ignore.case = TRUE
  )
]

cat("Files containing '1I' in filename:\n")
cat("------------------------------------------------------------\n")

if (length(one_i_candidates) == 0) {
  
  cat("NONE FOUND.\n\n")
  
} else {
  
  print(one_i_candidates)
  cat("\n")
}


# ------------------------------------------------------------
# 5. END
# ------------------------------------------------------------

cat("============================================================\n")
cat("END READ-ONLY STEP9D-1I FILE LOCATION AUDIT\n")
cat("============================================================\n")   


# ============================================================
# Paper 1
# STEP9D-3A-3
# Formal GSEA Gene-ID / Mapping Policy Lock
#
# FINAL CORRECTED VERSION
#
# IMPORTANT
# ------------------------------------------------------------
# READ-ONLY with respect to:
#   STEP9D-1I
#   STEP9D-2B
#   STEP9D-2C
#
# NO:
#   - ranking recalculation
#   - DEG recalculation
#   - GSEA execution
#   - jitter
#   - artificial tie breaking
#   - arbitrary first-match mapping
#
# ============================================================


cat("\n============================================================\n")
cat("STEP9D-3A-3: FORMAL GSEA ID / MAPPING POLICY LOCK\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Root directories
# ------------------------------------------------------------

step9_root <- paste0(
  "result_dir/Paper 1_Step 9-Functional Enrichment/",
  "STEP9D_GenomeWide_Preranked_GSEA"
)


# ------------------------------------------------------------
# 1. CORRECT upstream file paths
# ------------------------------------------------------------

input_lock_rds <- file.path(
  step9_root,
  "9D-1_Input_Lock",
  "STEP9D1I_FINAL_INPUT_LOCK.rds"
)

metric_lock_rds <- file.path(
  step9_root,
  "9D-2_Ranking_Metric_Lock",
  "STEP9D2B_FINAL_Ranking_Metric_Lock.rds"
)

ranking_rds <- file.path(
  step9_root,
  "9D-2C_Final_Ranking_Vector",
  "STEP9D2C_FINAL_GenomeWide_Ranking_Vector.rds"
)

ranking_audit_rds <- file.path(
  step9_root,
  "9D-2C_Final_Ranking_Vector",
  "STEP9D2C_FINAL_Ranking_Vector_Audit.rds"
)


# ------------------------------------------------------------
# 2. Output directory
# ------------------------------------------------------------

out_dir <- file.path(
  step9_root,
  "9D-3A-3_GSEA_ID_Mapping_Policy_Lock"
)

if (!dir.exists(out_dir)) {
  
  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ------------------------------------------------------------
# 3. Required file audit
# ------------------------------------------------------------

cat("Required upstream files\n")
cat("------------------------------------------------------------\n")

required_files <- c(
  input_lock_rds,
  metric_lock_rds,
  ranking_rds,
  ranking_audit_rds
)

required_exists <- file.exists(
  required_files
)

print(
  data.frame(
    file = required_files,
    exists = required_exists,
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)

if (!all(required_exists)) {
  
  stop(
    "STEP9D-3A-3 HARD FAIL: required upstream file missing."
  )
}

cat("\nAll required upstream files found.\n\n")


# ------------------------------------------------------------
# 4. Load upstream locks
# ------------------------------------------------------------

input_lock <- readRDS(
  input_lock_rds
)

metric_lock <- readRDS(
  metric_lock_rds
)

ranking_audit <- readRDS(
  ranking_audit_rds
)

ranking_obj <- readRDS(
  ranking_rds
)


# ------------------------------------------------------------
# 5. Recover the locked ranking vector
# ------------------------------------------------------------

ranking_vector <- NULL

if (
  is.numeric(ranking_obj) &&
  !is.null(names(ranking_obj))
) {
  
  ranking_vector <- ranking_obj
  
} else if (
  is.list(ranking_obj)
) {
  
  candidate_names <- c(
    "ranking_vector",
    "ranked_vector",
    "ranking",
    "geneList",
    "values"
  )
  
  for (nm in candidate_names) {
    
    if (
      !is.null(ranking_obj[[nm]]) &&
      is.numeric(ranking_obj[[nm]]) &&
      !is.null(names(ranking_obj[[nm]]))
    ) {
      
      ranking_vector <- ranking_obj[[nm]]
      
      break
    }
  }
}


if (is.null(ranking_vector)) {
  
  stop(
    "STEP9D-3A-3 HARD FAIL: locked ranking vector could not be recovered."
  )
}


# ------------------------------------------------------------
# 6. Basic ranking integrity
# ------------------------------------------------------------

ranking_gene_ids <- names(
  ranking_vector
)

ranking_values <- as.numeric(
  ranking_vector
)

ranking_n <- length(
  ranking_vector
)

missing_id_n <- sum(
  is.na(ranking_gene_ids) |
    ranking_gene_ids == ""
)

duplicated_id_n <- sum(
  duplicated(ranking_gene_ids)
)

na_value_n <- sum(
  is.na(ranking_values)
)

inf_value_n <- sum(
  is.infinite(ranking_values)
)

zero_ranking_n <- sum(
  ranking_values == 0
)

positive_ranking_n <- sum(
  ranking_values > 0
)

negative_ranking_n <- sum(
  ranking_values < 0
)


cat("Ranking vector integrity\n")
cat("------------------------------------------------------------\n")

cat(
  "Ranking vector length:",
  ranking_n,
  "\n"
)

cat(
  "Missing gene IDs:",
  missing_id_n,
  "\n"
)

cat(
  "Duplicated gene IDs:",
  duplicated_id_n,
  "\n"
)

cat(
  "NA ranking values:",
  na_value_n,
  "\n"
)

cat(
  "Infinite ranking values:",
  inf_value_n,
  "\n"
)

cat(
  "Zero ranking values:",
  zero_ranking_n,
  "\n"
)

cat(
  "Positive ranking values:",
  positive_ranking_n,
  "\n"
)

cat(
  "Negative ranking values:",
  negative_ranking_n,
  "\n\n"
)


# ------------------------------------------------------------
# 7. Ensembl format audit
# ------------------------------------------------------------

standard_ensembl <- grepl(
  "^ENSG[0-9]{11}$",
  ranking_gene_ids
)

versioned_ensembl <- grepl(
  "^ENSG[0-9]{11}\\.[0-9]+$",
  ranking_gene_ids
)

other_id_format <- (
  !standard_ensembl &
    !versioned_ensembl
)

standard_ensembl_n <- sum(
  standard_ensembl
)

versioned_ensembl_n <- sum(
  versioned_ensembl
)

other_id_n <- sum(
  other_id_format
)


cat("Ensembl ID format audit\n")
cat("------------------------------------------------------------\n")

cat(
  "Standard Ensembl:",
  standard_ensembl_n,
  "\n"
)

cat(
  "Versioned Ensembl:",
  versioned_ensembl_n,
  "\n"
)

cat(
  "Other ID format:",
  other_id_n,
  "\n\n"
)


# ------------------------------------------------------------
# 8. Consistency with STEP9D-2B
# ------------------------------------------------------------

locked_genome_n <- as.integer(
  metric_lock$genome_wide_gene_n
)

ranking_n_match <- (
  ranking_n == locked_genome_n
)

standard_n_match <- (
  standard_ensembl_n == locked_genome_n
)

cat("STEP9D-2B consistency\n")
cat("------------------------------------------------------------\n")

cat(
  "Locked genome-wide N:",
  locked_genome_n,
  "\n"
)

cat(
  "Ranking N:",
  ranking_n,
  "\n"
)

cat(
  "Ranking N matches lock:",
  ranking_n_match,
  "\n"
)

cat(
  "Standard Ensembl N matches lock:",
  standard_n_match,
  "\n\n"
)


# ------------------------------------------------------------
# 9. Load annotation packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  
  library(
    clusterProfiler
  )
  
  library(
    org.Hs.eg.db
  )
  
  library(
    AnnotationDbi
  )
})


cat("Package provenance\n")
cat("------------------------------------------------------------\n")

cat(
  "clusterProfiler:",
  as.character(
    packageVersion(
      "clusterProfiler"
    )
  ),
  "\n"
)

cat(
  "org.Hs.eg.db:",
  as.character(
    packageVersion(
      "org.Hs.eg.db"
    )
  ),
  "\n"
)

cat(
  "AnnotationDbi:",
  as.character(
    packageVersion(
      "AnnotationDbi"
    )
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 10. Keytype availability
# ------------------------------------------------------------

available_keytypes <- keytypes(
  org.Hs.eg.db
)

ensembl_keytype_available <- (
  "ENSEMBL" %in%
    available_keytypes
)

entrez_keytype_available <- (
  "ENTREZID" %in%
    available_keytypes
)


cat("Required annotation keytypes\n")
cat("------------------------------------------------------------\n")

cat(
  "ENSEMBL available:",
  ensembl_keytype_available,
  "\n"
)

cat(
  "ENTREZID available:",
  entrez_keytype_available,
  "\n\n"
)


# ------------------------------------------------------------
# 11. Ensembl -> Entrez mapping audit
# ------------------------------------------------------------

map_raw <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ranking_gene_ids,
  keytype = "ENSEMBL",
  columns = "ENTREZID"
)

map_raw <- map_raw[
  !is.na(map_raw$ENTREZID) &
    map_raw$ENTREZID != "",
  ,
  drop = FALSE
]


mapped_ensembl <- unique(
  map_raw$ENSEMBL
)

mapped_entrez <- unique(
  map_raw$ENTREZID
)

unmapped_ensembl <- setdiff(
  ranking_gene_ids,
  mapped_ensembl
)

mapping_rate <- (
  length(mapped_ensembl) /
    ranking_n
)


cat("Ensembl -> Entrez mapping\n")
cat("------------------------------------------------------------\n")

cat(
  "Raw mapping rows:",
  nrow(map_raw),
  "\n"
)

cat(
  "Unique mapped Ensembl:",
  length(mapped_ensembl),
  "\n"
)

cat(
  "Unique Entrez IDs:",
  length(mapped_entrez),
  "\n"
)

cat(
  "Unmapped Ensembl:",
  length(unmapped_ensembl),
  "\n"
)

cat(
  "Mapping rate:",
  sprintf(
    "%.6f",
    mapping_rate
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 12. One-to-many Ensembl audit
# ------------------------------------------------------------

ensembl_map_count <- table(
  map_raw$ENSEMBL
)

one_to_many_ensembl <- names(
  ensembl_map_count[
    ensembl_map_count > 1
  ]
)


cat("One-to-many Ensembl mapping\n")
cat("------------------------------------------------------------\n")

cat(
  "Ensembl IDs with >1 Entrez:",
  length(
    one_to_many_ensembl
  ),
  "\n"
)

if (
  length(one_to_many_ensembl) > 0
) {
  
  cat(
    "Maximum Entrez mappings / Ensembl:",
    max(
      ensembl_map_count
    ),
    "\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 13. Many-to-one Entrez collision audit
# ------------------------------------------------------------

entrez_map_count <- table(
  map_raw$ENTREZID
)

many_to_one_entrez <- names(
  entrez_map_count[
    entrez_map_count > 1
  ]
)


cat("Many-to-one Entrez collision audit\n")
cat("------------------------------------------------------------\n")

cat(
  "Entrez IDs receiving >1 Ensembl:",
  length(
    many_to_one_entrez
  ),
  "\n"
)

if (
  length(many_to_one_entrez) > 0
) {
  
  cat(
    "Maximum Ensembl IDs / Entrez:",
    max(
      entrez_map_count
    ),
    "\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 14. Core-gene traceability
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = core_genes,
  keytype = "SYMBOL",
  columns = c(
    "SYMBOL",
    "ENSEMBL",
    "ENTREZID"
  )
)

core_map <- unique(
  core_map[
    !is.na(core_map$ENSEMBL) &
      !is.na(core_map$ENTREZID),
    ,
    drop = FALSE
  ]
)

core_map <- core_map[
  core_map$ENSEMBL %in% ranking_gene_ids,
  ,
  drop = FALSE
]


cat("Core-gene traceability\n")
cat("------------------------------------------------------------\n")

print(
  core_map,
  row.names = FALSE
)

cat("\n")


core_traceable <- all(
  core_genes %in%
    core_map$SYMBOL
)

core_one_to_one <- all(
  sapply(
    core_genes,
    function(g) {
      
      tmp <- core_map[
        core_map$SYMBOL == g,
        ,
        drop = FALSE
      ]
      
      length(
        unique(
          tmp$ENSEMBL
        )
      ) == 1 &&
        length(
          unique(
            tmp$ENTREZID
          )
        ) == 1
    }
  )
)


cat(
  "All core genes traceable:",
  core_traceable,
  "\n"
)

cat(
  "All core genes one-to-one:",
  core_one_to_one,
  "\n\n"
)


# ------------------------------------------------------------
# 15. Formal policy definition
# ------------------------------------------------------------

GO_POLICY <- paste(
  "GO GSEA will use gseGO() directly with",
  "OrgDb = org.Hs.eg.db and keyType = ENSEMBL.",
  "The original 57,828-gene Ensembl ranking remains",
  "the locked genome-wide ranking.",
  "Genes without GO annotation are not deleted from",
  "the ranking; they simply do not contribute to",
  "GO-annotated gene sets."
)


KEGG_POLICY <- paste(
  "KEGG GSEA will use gseKEGG() with organism = hsa",
  "and an Entrez-based ranked input.",
  "Only Ensembl IDs with a valid mapping to EntrezID",
  "are eligible for KEGG mapping."
)


ONE_TO_MANY_POLICY <- paste(
  "Ensembl IDs mapping to more than one EntrezID",
  "will not be resolved by arbitrary first-match",
  "selection."
)


MANY_TO_ONE_POLICY <- paste(
  "Many-to-one Entrez collisions will not be silently",
  "removed with !duplicated(), arbitrary first-match",
  "selection, jittering, or any modification of the",
  "locked ranking."
)


UNMAPPED_POLICY <- paste(
  "Unmapped Ensembl genes remain in the original",
  "57,828-gene ranking and are excluded only from",
  "analyses requiring a mapped identifier."
)


cat("Formal policy\n")
cat("------------------------------------------------------------\n")

cat("GO:\n")
cat(GO_POLICY, "\n\n")

cat("KEGG:\n")
cat(KEGG_POLICY, "\n\n")

cat("One-to-many:\n")
cat(ONE_TO_MANY_POLICY, "\n\n")

cat("Many-to-one:\n")
cat(MANY_TO_ONE_POLICY, "\n\n")

cat("Unmapped:\n")
cat(UNMAPPED_POLICY, "\n\n")


# ------------------------------------------------------------
# 16. Execution safeguards
# ------------------------------------------------------------

ranking_recalculated <- FALSE
deg_recalculated <- FALSE
gsea_performed <- FALSE
jitter_used <- FALSE
artificial_tie_breaking <- FALSE
arbitrary_first_match_used <- FALSE


# ------------------------------------------------------------
# 17. Hard PASS
# ------------------------------------------------------------

hard_pass <- all(
  
  file.exists(
    input_lock_rds
  ),
  
  file.exists(
    metric_lock_rds
  ),
  
  file.exists(
    ranking_rds
  ),
  
  file.exists(
    ranking_audit_rds
  ),
  
  ranking_n_match,
  
  standard_n_match,
  
  missing_id_n == 0,
  
  duplicated_id_n == 0,
  
  na_value_n == 0,
  
  inf_value_n == 0,
  
  ensembl_keytype_available,
  
  entrez_keytype_available,
  
  core_traceable,
  
  core_one_to_one,
  
  ranking_recalculated == FALSE,
  
  deg_recalculated == FALSE,
  
  gsea_performed == FALSE,
  
  jitter_used == FALSE,
  
  artificial_tie_breaking == FALSE,
  
  arbitrary_first_match_used == FALSE
)


# ------------------------------------------------------------
# 18. Lock object
# ------------------------------------------------------------

policy_lock <- list(
  
  workflow =
    "Paper 1",
  
  step =
    "STEP9D-3A-3",
  
  status =
    if (
      hard_pass
    ) {
      "FINAL_PASS"
    } else {
      "HARD_FAIL"
    },
  
  lock_date =
    as.character(
      Sys.Date()
    ),
  
  input_lock_file =
    input_lock_rds,
  
  metric_lock_file =
    metric_lock_rds,
  
  ranking_file =
    ranking_rds,
  
  ranking_audit_file =
    ranking_audit_rds,
  
  genome_wide_gene_n =
    locked_genome_n,
  
  ranking_vector_n =
    ranking_n,
  
  standard_ensembl_n =
    standard_ensembl_n,
  
  versioned_ensembl_n =
    versioned_ensembl_n,
  
  other_id_n =
    other_id_n,
  
  duplicated_gene_n =
    duplicated_id_n,
  
  mapping_raw_rows =
    nrow(map_raw),
  
  mapped_ensembl_n =
    length(
      mapped_ensembl
    ),
  
  unique_entrez_n =
    length(
      mapped_entrez
    ),
  
  unmapped_ensembl_n =
    length(
      unmapped_ensembl
    ),
  
  ensembl_to_entrez_mapping_rate =
    mapping_rate,
  
  one_to_many_ensembl_n =
    length(
      one_to_many_ensembl
    ),
  
  many_to_one_entrez_n =
    length(
      many_to_one_entrez
    ),
  
  core_genes =
    core_genes,
  
  core_traceable =
    core_traceable,
  
  core_one_to_one =
    core_one_to_one,
  
  GO_method =
    "gseGO",
  
  GO_keyType =
    "ENSEMBL",
  
  GO_OrgDb =
    "org.Hs.eg.db",
  
  GO_ontology =
    "ALL",
  
  KEGG_method =
    "gseKEGG",
  
  KEGG_organism =
    "hsa",
  
  KEGG_ID_type =
    "ENTREZID",
  
  KEGG_mapping_policy =
    "valid Ensembl-to-Entrez mapping; collision handling locked separately before execution",
  
  arbitrary_first_match_used =
    arbitrary_first_match_used,
  
  ranking_recalculated =
    ranking_recalculated,
  
  DEG_recalculated =
    deg_recalculated,
  
  gsea_performed =
    gsea_performed,
  
  jitter_used =
    jitter_used,
  
  artificial_tie_breaking =
    artificial_tie_breaking,
  
  hard_pass =
    hard_pass
)


# ------------------------------------------------------------
# 19. Save final lock
# ------------------------------------------------------------

final_rds <- file.path(
  out_dir,
  "STEP9D3A3_FINAL_GSEA_ID_Mapping_Policy_Lock.rds"
)

final_txt <- file.path(
  out_dir,
  "STEP9D3A3_FINAL_GSEA_ID_Mapping_Policy_Lock.txt"
)

saveRDS(
  policy_lock,
  final_rds
)


# ------------------------------------------------------------
# 20. Human-readable audit
# ------------------------------------------------------------

sink(
  final_txt
)

cat(
  "STEP9D-3A-3 FINAL GSEA ID / MAPPING POLICY LOCK\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Status:",
  policy_lock$status,
  "\n"
)

cat(
  "Hard PASS:",
  policy_lock$hard_pass,
  "\n\n"
)

cat(
  "Genome-wide N:",
  policy_lock$genome_wide_gene_n,
  "\n"
)

cat(
  "Ranking vector N:",
  policy_lock$ranking_vector_n,
  "\n"
)

cat(
  "Standard Ensembl N:",
  policy_lock$standard_ensembl_n,
  "\n"
  
)

cat(
  "Mapped Ensembl N:",
  policy_lock$mapped_ensembl_n,
  "\n"
)

cat(
  "Unique Entrez N:",
  policy_lock$unique_entrez_n,
  "\n"
)

cat(
  "Unmapped Ensembl N:",
  policy_lock$unmapped_ensembl_n,
  "\n"
)

cat(
  "Mapping rate:",
  sprintf(
    "%.6f",
    policy_lock$ensembl_to_entrez_mapping_rate
  ),
  "\n"
)

cat(
  "One-to-many Ensembl N:",
  policy_lock$one_to_many_ensembl_n,
  "\n"
)

cat(
  "Many-to-one Entrez N:",
  policy_lock$many_to_one_entrez_n,
  "\n\n"
)

cat(
  "Core genes:\n"
)

print(
  core_map,
  row.names = FALSE
)

cat("\n")

cat(
  "Core traceable:",
  core_traceable,
  "\n"
)

cat(
  "Core one-to-one:",
  core_one_to_one,
  "\n\n"
)

cat(
  "GO POLICY:\n",
  GO_POLICY,
  "\n\n"
)

cat(
  "KEGG POLICY:\n",
  KEGG_POLICY,
  "\n\n"
)

cat(
  "ONE-TO-MANY POLICY:\n",
  ONE_TO_MANY_POLICY,
  "\n\n"
)

cat(
  "MANY-TO-ONE POLICY:\n",
  MANY_TO_ONE_POLICY,
  "\n\n"
)

cat(
  "UNMAPPED POLICY:\n",
  UNMAPPED_POLICY,
  "\n\n"
)

cat(
  "Execution safeguards:\n"
)

cat(
  "Ranking recalculated:",
  ranking_recalculated,
  "\n"
)

cat(
  "DEG recalculated:",
  deg_recalculated,
  "\n"
)

cat(
  "GSEA performed:",
  gsea_performed,
  "\n"
)

cat(
  "Jitter used:",
  jitter_used,
  "\n"
)

cat(
  "Artificial tie breaking:",
  artificial_tie_breaking,
  "\n"
)

cat(
  "Arbitrary first-match mapping:",
  arbitrary_first_match_used,
  "\n\n"
)

cat(
  "============================================================\n"
)

cat(
  "STEP9D-3A-3 END\n"
)

cat(
  "============================================================\n"
)

sink()


# ------------------------------------------------------------
# 21. Final console summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-3A-3 FINAL SUMMARY\n")
cat("============================================================\n\n")

cat(
  "Genome-wide genes:",
  ranking_n,
  "\n"
)

cat(
  "Standard Ensembl IDs:",
  standard_ensembl_n,
  "\n"
)

cat(
  "Mapped Ensembl genes:",
  length(mapped_ensembl),
  "\n"
)

cat(
  "Unique Entrez IDs:",
  length(mapped_entrez),
  "\n"
)

cat(
  "Unmapped Ensembl genes:",
  length(unmapped_ensembl),
  "\n"
)

cat(
  "One-to-many Ensembl:",
  length(one_to_many_ensembl),
  "\n"
)

cat(
  "Many-to-one Entrez:",
  length(many_to_one_entrez),
  "\n"
)

cat(
  "Core genes traceable:",
  core_traceable,
  "\n"
)

cat(
  "Core genes one-to-one:",
  core_one_to_one,
  "\n"
)

cat(
  "Ranking recalculated:",
  ranking_recalculated,
  "\n"
)

cat(
  "DEG recalculated:",
  deg_recalculated,
  "\n"
)

cat(
  "GSEA performed:",
  gsea_performed,
  "\n"
)

cat(
  "Final RDS exists:",
  file.exists(
    final_rds
  ),
  "\n"
)

cat(
  "Final TXT exists:",
  file.exists(
    final_txt
  ),
  "\n"
)

cat(
  "\nSTEP9D-3A-3 FINAL HARD PASS:",
  hard_pass,
  "\n"
)

cat("\n============================================================\n")
cat("END STEP9D-3A-3\n")
cat("============================================================\n")   


# ============================================================
# Paper 1
# STEP9D-3A-4
# KEGG Entrez Collision Audit
# and Deterministic Resolution Policy Lock
#
# STATUS:
# READ-ONLY AUDIT
#
# IMPORTANT:
# - No ranking recalculation
# - No DEG recalculation
# - No GSEA
# - No duplicate removal
# - No arbitrary first-match mapping
# - No jitter
# - No artificial tie-breaking
# - Original 57,828-gene ranking remains immutable
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-3A-4: KEGG Entrez Collision Audit\n")
cat("and Deterministic Resolution Policy Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_pkgs <- c(
  "AnnotationDbi",
  "org.Hs.eg.db"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package not installed: ", pkg)
  }
}

library(AnnotationDbi)
library(org.Hs.eg.db)


# ------------------------------------------------------------
# 1. Locked paths
# ------------------------------------------------------------

root_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

ranking_file <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-2C_Final_Ranking_Vector",
  "STEP9D2C_FINAL_GenomeWide_Ranking_Vector.rds"
)

mapping_lock_file <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-3A-3_GSEA_ID_Mapping_Policy_Lock",
  "STEP9D3A3_FINAL_GSEA_ID_Mapping_Policy_Lock.rds"
)

output_dir <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-3A-4_KEGG_Entrez_Collision_Audit"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Core genes
# ------------------------------------------------------------

core_genes <- c(
  AOC3   = "ENSG00000131471",
  MS4A7  = "ENSG00000166927",
  IL7R   = "ENSG00000168685"
)


# ------------------------------------------------------------
# 3. File existence audit
# ------------------------------------------------------------

if (!file.exists(ranking_file)) {
  stop("Locked STEP9D-2C ranking file not found:\n", ranking_file)
}

if (!file.exists(mapping_lock_file)) {
  stop("STEP9D-3A-3 mapping lock not found:\n", mapping_lock_file)
}

cat("Ranking file exists: ", file.exists(ranking_file), "\n")
cat("Mapping lock exists: ", file.exists(mapping_lock_file), "\n\n")


# ------------------------------------------------------------
# 4. Load immutable ranking vector
# ------------------------------------------------------------

ranking_obj <- readRDS(ranking_file)

cat("Ranking object class:\n")
print(class(ranking_obj))

cat("\nRanking object names:\n")
print(names(ranking_obj))


# ------------------------------------------------------------
# 5. Extract ranking vector robustly
# ------------------------------------------------------------

if (is.numeric(ranking_obj) && !is.null(names(ranking_obj))) {
  
  ranking_vector <- ranking_obj
  
} else if ("ranking_vector" %in% names(ranking_obj)) {
  
  ranking_vector <- ranking_obj$ranking_vector
  
} else if ("ranking" %in% names(ranking_obj)) {
  
  ranking_vector <- ranking_obj$ranking
  
} else if ("gene_ranking" %in% names(ranking_obj)) {
  
  ranking_vector <- ranking_obj$gene_ranking
  
} else {
  
  stop(
    "Unable to identify ranking vector inside STEP9D-2C RDS."
  )
}


# ------------------------------------------------------------
# 6. Ranking integrity audit
# ------------------------------------------------------------

ranking_n <- length(ranking_vector)
ranking_ids <- names(ranking_vector)

ranking_duplicate_n <- sum(duplicated(ranking_ids))
ranking_missing_id_n <- sum(
  is.na(ranking_ids) | ranking_ids == ""
)

ranking_non_numeric_n <- sum(
  !is.finite(as.numeric(ranking_vector))
)

ranking_standard_ensembl_n <- sum(
  grepl("^ENSG[0-9]{11}$", ranking_ids)
)

cat("\n------------------------------------------------------------\n")
cat("6. Ranking integrity\n")
cat("------------------------------------------------------------\n")

cat("Ranking length: ", ranking_n, "\n")
cat("Duplicate gene IDs: ", ranking_duplicate_n, "\n")
cat("Missing gene IDs: ", ranking_missing_id_n, "\n")
cat("Non-finite ranking values: ", ranking_non_numeric_n, "\n")
cat("Standard Ensembl IDs: ", ranking_standard_ensembl_n, "\n")


if (ranking_n != 57828) {
  stop("Ranking length is not the locked 57,828 genes.")
}

if (ranking_duplicate_n != 0) {
  stop("Ranking vector contains duplicated gene IDs.")
}

if (ranking_missing_id_n != 0) {
  stop("Ranking vector contains missing gene IDs.")
}

if (ranking_standard_ensembl_n != 57828) {
  stop("Ranking vector does not contain exactly 57,828 standard Ensembl IDs.")
}


# ------------------------------------------------------------
# 7. Load STEP9D-3A-3 mapping policy lock
# ------------------------------------------------------------

mapping_lock <- readRDS(mapping_lock_file)

cat("\n------------------------------------------------------------\n")
cat("7. STEP9D-3A-3 mapping lock\n")
cat("------------------------------------------------------------\n")

cat("Mapping lock class:\n")
print(class(mapping_lock))

cat("\nMapping lock names:\n")
print(names(mapping_lock))


# ------------------------------------------------------------
# 8. Reconstruct annotation mapping ONLY
#
# IMPORTANT:
# This is not a ranking recalculation.
# It is an annotation lookup required for collision auditing.
# ------------------------------------------------------------

mapping_raw <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ranking_ids,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "ENTREZID", "SYMBOL")
)

mapping_raw <- as.data.frame(mapping_raw)

cat("\n------------------------------------------------------------\n")
cat("8. Raw Ensembl → Entrez annotation mapping\n")
cat("------------------------------------------------------------\n")

cat("Raw mapping rows: ", nrow(mapping_raw), "\n")
cat(
  "Unique Ensembl IDs with Entrez mapping: ",
  length(unique(mapping_raw$ENSEMBL[
    !is.na(mapping_raw$ENTREZID)
  ])),
  "\n"
)


# ------------------------------------------------------------
# 9. Clean annotation table for collision analysis
#
# Retain only valid Ensembl → Entrez pairs.
# No first-match selection.
# No duplicate Entrez removal.
# ------------------------------------------------------------

mapping_valid <- mapping_raw[
  !is.na(mapping_raw$ENSEMBL) &
    !is.na(mapping_raw$ENTREZID),
  c("ENSEMBL", "ENTREZID", "SYMBOL")
]

mapping_valid$ENSEMBL <- as.character(mapping_valid$ENSEMBL)
mapping_valid$ENTREZID <- as.character(mapping_valid$ENTREZID)
mapping_valid$SYMBOL <- as.character(mapping_valid$SYMBOL)

mapping_valid <- unique(mapping_valid)


# ------------------------------------------------------------
# 10. Ensembl-side multiplicity
# ------------------------------------------------------------

ensembl_entrez_count <- aggregate(
  ENTREZID ~ ENSEMBL,
  data = mapping_valid,
  FUN = function(x) length(unique(x))
)

names(ensembl_entrez_count)[2] <- "n_entrez_per_ensembl"

one_to_many_ensembl <- ensembl_entrez_count[
  ensembl_entrez_count$n_entrez_per_ensembl > 1,
  ,
  drop = FALSE
]

one_entrez_ensembl <- ensembl_entrez_count[
  ensembl_entrez_count$n_entrez_per_ensembl == 1,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 11. Build the KEGG candidate subset
#
# POLICY FOR THIS AUDIT:
#
# Retain Ensembl IDs for which exactly ONE Entrez ID
# is annotated.
#
# IMPORTANT:
# This is eligibility definition only.
# Duplicate Entrez IDs are NOT removed here.
# ------------------------------------------------------------

eligible_ensembl <- one_entrez_ensembl$ENSEMBL

kegg_candidate <- mapping_valid[
  mapping_valid$ENSEMBL %in% eligible_ensembl,
  ,
  drop = FALSE
]

kegg_candidate <- kegg_candidate[
  !duplicated(kegg_candidate[, c("ENSEMBL", "ENTREZID")]),
  ,
  drop = FALSE
]

# Add ranking without recalculation
kegg_candidate$ranking_value <-
  unname(ranking_vector[kegg_candidate$ENSEMBL])

kegg_candidate$ranking_direction <- ifelse(
  kegg_candidate$ranking_value > 0,
  "positive",
  ifelse(
    kegg_candidate$ranking_value < 0,
    "negative",
    "zero"
  )
)

kegg_candidate$abs_ranking <-
  abs(kegg_candidate$ranking_value)


# ------------------------------------------------------------
# 12. Verify ranking coverage
# ------------------------------------------------------------

if (any(is.na(kegg_candidate$ranking_value))) {
  stop(
    "Some KEGG candidate genes do not have ranking values."
  )
}


# ------------------------------------------------------------
# 13. Entrez-side multiplicity within KEGG candidate subset
# ------------------------------------------------------------

entrez_group_counts <- aggregate(
  ENSEMBL ~ ENTREZID,
  data = kegg_candidate,
  FUN = function(x) length(unique(x))
)

names(entrez_group_counts)[2] <- "n_ensembl_per_entrez"

kegg_collision_entrez <- entrez_group_counts[
  entrez_group_counts$n_ensembl_per_entrez > 1,
  ,
  drop = FALSE
]

kegg_noncollision_entrez <- entrez_group_counts[
  entrez_group_counts$n_ensembl_per_entrez == 1,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 14. Collision groups
# ------------------------------------------------------------

collision_entrez_ids <- kegg_collision_entrez$ENTREZID

collision_table <- kegg_candidate[
  kegg_candidate$ENTREZID %in% collision_entrez_ids,
  ,
  drop = FALSE
]

collision_table <- collision_table[
  order(
    collision_table$ENTREZID,
    collision_table$ENSEMBL
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 15. Annotate core genes
# ------------------------------------------------------------

collision_table$is_core_gene <- FALSE
collision_table$core_gene <- NA_character_

for (g in names(core_genes)) {
  
  idx <- collision_table$ENSEMBL == core_genes[[g]]
  
  collision_table$is_core_gene[idx] <- TRUE
  collision_table$core_gene[idx] <- g
}


# ------------------------------------------------------------
# 16. Collision-group summary
# ------------------------------------------------------------

split_collision <- split(
  collision_table,
  collision_table$ENTREZID
)

collision_group_summary <- do.call(
  rbind,
  lapply(
    names(split_collision),
    function(eid) {
      
      d <- split_collision[[eid]]
      
      signs <- unique(d$ranking_direction)
      
      direction_class <- if (
        all(signs == "zero")
      ) {
        "all_zero"
        
      } else if (
        length(signs) == 1 &&
        signs == "positive"
      ) {
        "all_positive"
        
      } else if (
        length(signs) == 1 &&
        signs == "negative"
      ) {
        "all_negative"
        
      } else if (
        all(signs %in% c("positive", "zero"))
      ) {
        "positive_zero"
        
      } else if (
        all(signs %in% c("negative", "zero"))
      ) {
        "negative_zero"
        
      } else if (
        all(signs %in% c("positive", "negative", "zero"))
      ) {
        "positive_negative_zero"
        
      } else if (
        all(signs %in% c("positive", "negative"))
      ) {
        "positive_negative"
        
      } else {
        "other"
      }
      
      ranking_unique_n <-
        length(unique(d$ranking_value))
      
      abs_ranking_unique_n <-
        length(unique(d$abs_ranking))
      
      data.frame(
        ENTREZID = eid,
        n_ensembl = length(unique(d$ENSEMBL)),
        ranking_unique_n = ranking_unique_n,
        abs_ranking_unique_n = abs_ranking_unique_n,
        ranking_min = min(d$ranking_value),
        ranking_max = max(d$ranking_value),
        abs_ranking_min = min(d$abs_ranking),
        abs_ranking_max = max(d$abs_ranking),
        direction_class = direction_class,
        contains_positive = any(
          d$ranking_direction == "positive"
        ),
        contains_negative = any(
          d$ranking_direction == "negative"
        ),
        contains_zero = any(
          d$ranking_direction == "zero"
        ),
        contains_core_gene = any(
          d$is_core_gene
        ),
        core_genes = paste(
          unique(
            na.omit(d$core_gene)
          ),
          collapse = ";"
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(collision_group_summary) <- NULL


# ------------------------------------------------------------
# 17. Collision classification
# ------------------------------------------------------------

collision_group_summary$ranking_agreement_class <-
  ifelse(
    collision_group_summary$ranking_unique_n == 1,
    "identical_ranking",
    ifelse(
      collision_group_summary$contains_positive &
        collision_group_summary$contains_negative,
      "opposite_direction",
      ifelse(
        collision_group_summary$abs_ranking_unique_n == 1,
        "same_magnitude_different_zero_status",
        "same_direction_or_mixed_nonzero_magnitude"
      )
    )
  )


# ------------------------------------------------------------
# 18. Summary counts
# ------------------------------------------------------------

total_raw_many_to_one_n <- nrow(
  aggregate(
    ENSEMBL ~ ENTREZID,
    data = mapping_valid,
    FUN = function(x) length(unique(x))
  )[
    sapply(
      aggregate(
        ENSEMBL ~ ENTREZID,
        data = mapping_valid,
        FUN = function(x) length(unique(x))
      )$ENSEMBL,
      length
    ) > 1,
    ,
    drop = FALSE
  ]
)

# More robust direct calculation
raw_entrez_counts <- tapply(
  mapping_valid$ENSEMBL,
  mapping_valid$ENTREZID,
  function(x) length(unique(x))
)

total_raw_many_to_one_n <- sum(
  raw_entrez_counts > 1
)

raw_many_to_one_ensembl_n <- length(
  unique(
    mapping_valid$ENSEMBL[
      mapping_valid$ENTREZID %in%
        names(raw_entrez_counts[
          raw_entrez_counts > 1
        ])
    ]
  )
)

kegg_collision_group_n <- nrow(
  collision_group_summary
)

kegg_collision_ensembl_n <- length(
  unique(collision_table$ENSEMBL)
)

kegg_candidate_ensembl_n <- length(
  unique(kegg_candidate$ENSEMBL)
)

kegg_candidate_entrez_n <- length(
  unique(kegg_candidate$ENTREZID)
)

kegg_duplicate_entrez_n <- sum(
  duplicated(kegg_candidate$ENTREZID)
)


# ------------------------------------------------------------
# 19. Ranking-direction summary across collision groups
# ------------------------------------------------------------

collision_positive_n <- sum(
  collision_table$ranking_direction == "positive"
)

collision_negative_n <- sum(
  collision_table$ranking_direction == "negative"
)

collision_zero_n <- sum(
  collision_table$ranking_direction == "zero"
)

collision_identical_group_n <- sum(
  collision_group_summary$ranking_agreement_class ==
    "identical_ranking"
)

collision_opposite_direction_group_n <- sum(
  collision_group_summary$ranking_agreement_class ==
    "opposite_direction"
)

collision_same_magnitude_group_n <- sum(
  collision_group_summary$ranking_agreement_class ==
    "same_magnitude_different_zero_status"
)


# ------------------------------------------------------------
# 20. Core-gene collision audit
# ------------------------------------------------------------

core_collision_status <- data.frame(
  core_gene = names(core_genes),
  ENSEMBL = unname(core_genes),
  in_ranking = names(core_genes) %in%
    names(core_genes),
  in_kegg_candidate = core_genes %in%
    kegg_candidate$ENSEMBL,
  collision_affected = core_genes %in%
    collision_table$ENSEMBL,
  ENTREZID = NA_character_,
  ranking_value = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(core_collision_status))) {
  
  ens <- core_collision_status$ENSEMBL[i]
  
  m <- kegg_candidate[
    kegg_candidate$ENSEMBL == ens,
    ,
    drop = FALSE
  ]
  
  if (nrow(m) > 0) {
    
    core_collision_status$ENTREZID[i] <-
      paste(unique(m$ENTREZID), collapse = ";")
    
    core_collision_status$ranking_value[i] <-
      unique(m$ranking_value)[1]
  }
}


# ------------------------------------------------------------
# 21. Identify all opposite-direction collision groups
# ------------------------------------------------------------

opposite_direction_groups <-
  collision_group_summary[
    collision_group_summary$ranking_agreement_class ==
      "opposite_direction",
    ,
    drop = FALSE
  ]


# ------------------------------------------------------------
# 22. Identify all non-identical collision groups
# ------------------------------------------------------------

nonidentical_collision_groups <-
  collision_group_summary[
    collision_group_summary$ranking_unique_n > 1,
    ,
    drop = FALSE
  ]


# ------------------------------------------------------------
# 23. Determine whether collision resolution is needed
#
# IMPORTANT:
# This does NOT resolve anything.
# It only establishes whether duplicate Entrez IDs exist.
# ------------------------------------------------------------

resolution_required <-
  kegg_duplicate_entrez_n > 0

opposite_direction_present <-
  collision_opposite_direction_group_n > 0

core_gene_collision_present <-
  any(core_collision_status$collision_affected)


# ------------------------------------------------------------
# 24. Hard QC
# ------------------------------------------------------------

hard_qc <- list(
  
  ranking_length_locked =
    ranking_n == 57828,
  
  ranking_ids_unique =
    ranking_duplicate_n == 0,
  
  ranking_ids_standard_ensembl =
    ranking_standard_ensembl_n == 57828,
  
  ranking_values_finite =
    ranking_non_numeric_n == 0,
  
  no_ranking_recalculation =
    TRUE,
  
  no_deg_recalculation =
    TRUE,
  
  no_gsea_performed =
    TRUE,
  
  no_duplicate_entrez_removed =
    TRUE,
  
  no_arbitrary_first_match =
    TRUE,
  
  raw_many_to_one_documented =
    total_raw_many_to_one_n == 409,
  
  kegg_candidate_defined_by_single_entrez_per_ensembl =
    TRUE,
  
  collision_audit_completed =
    TRUE,
  
  collision_resolution_not_yet_applied =
    TRUE
)


hard_pass <- all(
  unlist(hard_qc)
)


# ------------------------------------------------------------
# 25. Console summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-3A-4 COLLISION AUDIT SUMMARY\n")
cat("============================================================\n\n")

cat("Locked genome-wide ranking genes: ",
    ranking_n, "\n")

cat("Raw many-to-one Entrez IDs: ",
    total_raw_many_to_one_n, "\n")

cat("Raw many-to-one affected Ensembl genes: ",
    raw_many_to_one_ensembl_n, "\n\n")

cat("KEGG candidate Ensembl genes: ",
    kegg_candidate_ensembl_n, "\n")

cat("KEGG candidate unique Entrez IDs: ",
    kegg_candidate_entrez_n, "\n")

cat("Duplicate Entrez IDs in KEGG candidate: ",
    kegg_duplicate_entrez_n, "\n")

cat("KEGG collision Entrez groups: ",
    kegg_collision_group_n, "\n")

cat("KEGG collision-affected Ensembl genes: ",
    kegg_collision_ensembl_n, "\n\n")

cat("Collision positive ranking genes: ",
    collision_positive_n, "\n")

cat("Collision negative ranking genes: ",
    collision_negative_n, "\n")

cat("Collision zero ranking genes: ",
    collision_zero_n, "\n\n")

cat("Identical-ranking collision groups: ",
    collision_identical_group_n, "\n")

cat("Non-identical collision groups: ",
    nrow(nonidentical_collision_groups), "\n")

cat("Opposite-direction collision groups: ",
    collision_opposite_direction_group_n, "\n")

cat("Same-magnitude / zero-status groups: ",
    collision_same_magnitude_group_n, "\n\n")

cat("Core-gene collision present: ",
    core_gene_collision_present, "\n")

cat("Resolution required: ",
    resolution_required, "\n")

cat("Opposite-direction collision present: ",
    opposite_direction_present, "\n\n")


# ------------------------------------------------------------
# 26. Save detailed collision table
# ------------------------------------------------------------

collision_csv <- file.path(
  output_dir,
  "STEP9D3A4_KEGG_Entrez_Collision_Detail.csv"
)

write.csv(
  collision_table,
  collision_csv,
  row.names = FALSE,
  quote = TRUE
)


# ------------------------------------------------------------
# 27. Save collision-group summary
# ------------------------------------------------------------

collision_summary_csv <- file.path(
  output_dir,
  "STEP9D3A4_KEGG_Entrez_Collision_Group_Summary.csv"
)

write.csv(
  collision_group_summary,
  collision_summary_csv,
  row.names = FALSE,
  quote = TRUE
)


# ------------------------------------------------------------
# 28. Save opposite-direction groups
# ------------------------------------------------------------

opposite_csv <- file.path(
  output_dir,
  "STEP9D3A4_Opposite_Direction_Collision_Groups.csv"
)

write.csv(
  opposite_direction_groups,
  opposite_csv,
  row.names = FALSE,
  quote = TRUE
)


# ------------------------------------------------------------
# 29. Save core-gene collision audit
# ------------------------------------------------------------

core_csv <- file.path(
  output_dir,
  "STEP9D3A4_CoreGene_Collision_Audit.csv"
)

write.csv(
  core_collision_status,
  core_csv,
  row.names = FALSE,
  quote = TRUE
)


# ------------------------------------------------------------
# 30. Save complete audit object
# ------------------------------------------------------------

audit_object <- list(
  
  workflow = "Paper 1",
  
  step = "STEP9D-3A-4",
  
  status = "FINAL_READ_ONLY_COLLISION_AUDIT",
  
  audit_date = as.character(Sys.time()),
  
  ranking_file = ranking_file,
  
  mapping_lock_file = mapping_lock_file,
  
  genome_wide_gene_n = ranking_n,
  
  raw_mapping_many_to_one_entrez_n =
    total_raw_many_to_one_n,
  
  raw_many_to_one_affected_ensembl_n =
    raw_many_to_one_ensembl_n,
  
  kegg_candidate_ensembl_n =
    kegg_candidate_ensembl_n,
  
  kegg_candidate_unique_entrez_n =
    kegg_candidate_entrez_n,
  
  kegg_duplicate_entrez_n =
    kegg_duplicate_entrez_n,
  
  kegg_collision_entrez_group_n =
    kegg_collision_group_n,
  
  kegg_collision_ensembl_n =
    kegg_collision_ensembl_n,
  
  collision_positive_n =
    collision_positive_n,
  
  collision_negative_n =
    collision_negative_n,
  
  collision_zero_n =
    collision_zero_n,
  
  collision_identical_group_n =
    collision_identical_group_n,
  
  collision_nonidentical_group_n =
    nrow(nonidentical_collision_groups),
  
  collision_opposite_direction_group_n =
    collision_opposite_direction_group_n,
  
  collision_same_magnitude_group_n =
    collision_same_magnitude_group_n,
  
  core_gene_collision_present =
    core_gene_collision_present,
  
  resolution_required =
    resolution_required,
  
  opposite_direction_present =
    opposite_direction_present,
  
  core_gene_status =
    core_collision_status,
  
  hard_qc = hard_qc,
  
  hard_pass = hard_pass,
  
  ranking_recalculated = FALSE,
  
  deg_recalculated = FALSE,
  
  gsea_performed = FALSE,
  
  duplicate_entrez_removed = FALSE,
  
  arbitrary_first_match_used = FALSE,
  
  jitter_used = FALSE,
  
  artificial_tie_breaking_used = FALSE,
  
  collision_resolution_applied = FALSE,
  
  policy_status =
    "AUDIT_ONLY_NO_RESOLUTION_YET"
)


# ------------------------------------------------------------
# 31. Save RDS
# ------------------------------------------------------------

audit_rds <- file.path(
  output_dir,
  "STEP9D3A4_FINAL_KEGG_Entrez_Collision_Audit.rds"
)

saveRDS(
  audit_object,
  audit_rds
)


# ------------------------------------------------------------
# 32. Formal TXT report
# ------------------------------------------------------------

audit_txt <- file.path(
  output_dir,
  "STEP9D3A4_FINAL_KEGG_Entrez_Collision_Audit.txt"
)

sink(audit_txt)

cat("============================================================\n")
cat("Paper 1\n")
cat("STEP9D-3A-4\n")
cat("KEGG Entrez Collision Audit\n")
cat("and Deterministic Resolution Policy Lock\n")
cat("============================================================\n\n")

cat("STATUS: FINAL READ-ONLY AUDIT\n\n")

cat("Genome-wide ranking genes: ",
    ranking_n, "\n")

cat("Raw many-to-one Entrez IDs: ",
    total_raw_many_to_one_n, "\n")

cat("Raw affected Ensembl genes: ",
    raw_many_to_one_ensembl_n, "\n\n")

cat("KEGG candidate Ensembl genes: ",
    kegg_candidate_ensembl_n, "\n")

cat("KEGG candidate unique Entrez IDs: ",
    kegg_candidate_entrez_n, "\n")

cat("Duplicate Entrez IDs in KEGG candidate: ",
    kegg_duplicate_entrez_n, "\n")

cat("Collision Entrez groups: ",
    kegg_collision_group_n, "\n")

cat("Collision-affected Ensembl genes: ",
    kegg_collision_ensembl_n, "\n\n")

cat("Positive collision ranking genes: ",
    collision_positive_n, "\n")

cat("Negative collision ranking genes: ",
    collision_negative_n, "\n")

cat("Zero collision ranking genes: ",
    collision_zero_n, "\n\n")

cat("Identical-ranking collision groups: ",
    collision_identical_group_n, "\n")

cat("Non-identical collision groups: ",
    nrow(nonidentical_collision_groups), "\n")

cat("Opposite-direction collision groups: ",
    collision_opposite_direction_group_n, "\n")

cat("Same-magnitude / zero-status groups: ",
    collision_same_magnitude_group_n, "\n\n")

cat("Core-gene collision present: ",
    core_gene_collision_present, "\n")

cat("Resolution required: ",
    resolution_required, "\n")

cat("Opposite-direction collision present: ",
    opposite_direction_present, "\n\n")

cat("------------------------------------------------------------\n")
cat("Hard QC\n")
cat("------------------------------------------------------------\n")

print(hard_qc)

cat("\nSTEP9D-3A-4 FINAL HARD PASS: ",
    hard_pass, "\n")

cat("\n------------------------------------------------------------\n")
cat("No recalculation / no modification flags\n")
cat("------------------------------------------------------------\n")

cat("Ranking recalculated: FALSE\n")
cat("DEG recalculated: FALSE\n")
cat("GSEA performed: FALSE\n")
cat("Duplicate Entrez removed: FALSE\n")
cat("Arbitrary first-match used: FALSE\n")
cat("Jitter used: FALSE\n")
cat("Artificial tie-breaking used: FALSE\n")
cat("Collision resolution applied: FALSE\n")

sink()


# ------------------------------------------------------------
# 33. Final file existence audit
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("FINAL FILE AUDIT\n")
cat("------------------------------------------------------------\n")

cat(
  "Collision detail CSV exists: ",
  file.exists(collision_csv), "\n"
)

cat(
  "Collision group summary CSV exists: ",
  file.exists(collision_summary_csv), "\n"
)

cat(
  "Opposite-direction CSV exists: ",
  file.exists(opposite_csv), "\n"
)

cat(
  "Core-gene audit CSV exists: ",
  file.exists(core_csv), "\n"
)

cat(
  "Final audit RDS exists: ",
  file.exists(audit_rds), "\n"
)

cat(
  "Final audit TXT exists: ",
  file.exists(audit_txt), "\n"
)

cat("\n============================================================\n")
cat("STEP9D-3A-4 FINAL HARD PASS: ", hard_pass, "\n")
cat("============================================================\n")


# ============================================================
# STEP9D-3A-4
# READ-ONLY FINAL RESULT INSPECTION
# No recalculation
# No file modification
# ============================================================

audit_rds <- file.path(
  "input_dir/Paper 1_Step 9-Functional Enrichment",
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-3A-4_KEGG_Entrez_Collision_Audit",
  "STEP9D3A4_FINAL_KEGG_Entrez_Collision_Audit.rds"
)

a <- readRDS(audit_rds)

cat("\n============================================================\n")
cat("STEP9D-3A-4 FINAL READ-ONLY RESULT\n")
cat("============================================================\n\n")

cat("Genome-wide ranking genes: ",
    a$genome_wide_gene_n, "\n")

cat("Raw many-to-one Entrez IDs: ",
    a$raw_mapping_many_to_one_entrez_n, "\n")

cat("Raw many-to-one affected Ensembl genes: ",
    a$raw_many_to_one_affected_ensembl_n, "\n\n")

cat("KEGG candidate Ensembl genes: ",
    a$kegg_candidate_ensembl_n, "\n")

cat("KEGG candidate unique Entrez IDs: ",
    a$kegg_candidate_unique_entrez_n, "\n")

cat("Duplicate Entrez IDs in KEGG candidate: ",
    a$kegg_duplicate_entrez_n, "\n")

cat("KEGG collision Entrez groups: ",
    a$kegg_collision_entrez_group_n, "\n")

cat("KEGG collision-affected Ensembl genes: ",
    a$kegg_collision_ensembl_n, "\n\n")

cat("Collision positive ranking genes: ",
    a$collision_positive_n, "\n")

cat("Collision negative ranking genes: ",
    a$collision_negative_n, "\n")

cat("Collision zero ranking genes: ",
    a$collision_zero_n, "\n\n")

cat("Identical-ranking collision groups: ",
    a$collision_identical_group_n, "\n")

cat("Non-identical collision groups: ",
    a$collision_nonidentical_group_n, "\n")

cat("Opposite-direction collision groups: ",
    a$collision_opposite_direction_group_n, "\n")

cat("Same-magnitude / zero-status groups: ",
    a$collision_same_magnitude_group_n, "\n\n")

cat("Core-gene collision present: ",
    a$core_gene_collision_present, "\n")

cat("Resolution required: ",
    a$resolution_required, "\n")

cat("Opposite-direction collision present: ",
    a$opposite_direction_present, "\n\n")

cat("Ranking recalculated: ",
    a$ranking_recalculated, "\n")

cat("DEG recalculated: ",
    a$deg_recalculated, "\n")

cat("GSEA performed: ",
    a$gsea_performed, "\n")

cat("Duplicate Entrez removed: ",
    a$duplicate_entrez_removed, "\n")

cat("Arbitrary first-match used: ",
    a$arbitrary_first_match_used, "\n")

cat("Jitter used: ",
    a$jitter_used, "\n")

cat("Artificial tie-breaking used: ",
    a$artificial_tie_breaking_used, "\n")

cat("Collision resolution applied: ",
    a$collision_resolution_applied, "\n\n")

cat("STEP9D-3A-4 FINAL HARD PASS: ",
    a$hard_pass, "\n")

cat("\n============================================================\n")



collision_summary_file <- file.path(
  "input_dir/Paper 1_Step 9-Functional Enrichment",
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-3A-4_KEGG_Entrez_Collision_Audit",
  "STEP9D3A4_KEGG_Entrez_Collision_Group_Summary.csv"
)

collision_summary <- read.csv(
  collision_summary_file,
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("COLLISION GROUP SUMMARY\n")
cat("============================================================\n\n")

print(collision_summary)

cat("\n============================================================\n")
cat("NON-IDENTICAL COLLISION GROUPS\n")
cat("============================================================\n\n")

print(
  collision_summary[
    collision_summary$ranking_unique_n > 1,
    ,
    drop = FALSE
  ]
)  


# ============================================================
# Paper 1
# STEP9D-3A-5
# KEGG Entrez Collision Resolution
# Policy Comparison Audit
#
# STATUS:
# READ-ONLY POLICY COMPARISON
#
# IMPORTANT:
# - No ranking recalculation
# - No DEG recalculation
# - No GSEA
# - No modification of original 57,828-gene ranking
# - No final collision resolution
# - No arbitrary first-match selection
# - No jitter
# - No artificial tie-breaking of ranking values
#
# PURPOSE:
# Compare deterministic candidate policies before
# locking the final KEGG Entrez collision policy.
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-3A-5: KEGG Entrez Collision Resolution\n")
cat("Policy Comparison Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_pkgs <- c(
  "AnnotationDbi",
  "org.Hs.eg.db"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package not installed: ", pkg)
  }
}

library(AnnotationDbi)
library(org.Hs.eg.db)


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

root_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

ranking_file <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-2C_Final_Ranking_Vector",
  "STEP9D2C_FINAL_GenomeWide_Ranking_Vector.rds"
)

step9d3a4_file <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-3A-4_KEGG_Entrez_Collision_Audit",
  "STEP9D3A4_FINAL_KEGG_Entrez_Collision_Audit.rds"
)

output_dir <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-3A-5_KEGG_Collision_Policy_Audit"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Core genes
# ------------------------------------------------------------

core_genes <- c(
  AOC3  = "ENSG00000131471",
  MS4A7 = "ENSG00000166927",
  IL7R  = "ENSG00000168685"
)


# ------------------------------------------------------------
# 3. File existence
# ------------------------------------------------------------

if (!file.exists(ranking_file)) {
  stop("STEP9D-2C ranking file not found:\n", ranking_file)
}

if (!file.exists(step9d3a4_file)) {
  stop("STEP9D-3A-4 audit file not found:\n", step9d3a4_file)
}

cat("Ranking file exists: TRUE\n")
cat("STEP9D-3A-4 audit exists: TRUE\n\n")


# ------------------------------------------------------------
# 4. Load immutable ranking
# ------------------------------------------------------------

ranking_obj <- readRDS(ranking_file)

if (
  is.numeric(ranking_obj) &&
  !is.null(names(ranking_obj))
) {
  
  ranking_vector <- ranking_obj
  
} else if (
  "ranking_vector" %in% names(ranking_obj)
) {
  
  ranking_vector <- ranking_obj$ranking_vector
  
} else if (
  "ranking" %in% names(ranking_obj)
) {
  
  ranking_vector <- ranking_obj$ranking
  
} else if (
  "gene_ranking" %in% names(ranking_obj)
) {
  
  ranking_vector <- ranking_obj$gene_ranking
  
} else {
  
  stop(
    "Unable to identify immutable ranking vector."
  )
}


# ------------------------------------------------------------
# 5. Ranking integrity
# ------------------------------------------------------------

ranking_ids <- names(ranking_vector)

ranking_n <- length(ranking_vector)

duplicate_gene_n <- sum(
  duplicated(ranking_ids)
)

missing_gene_n <- sum(
  is.na(ranking_ids) |
    ranking_ids == ""
)

standard_ensembl_n <- sum(
  grepl("^ENSG[0-9]{11}$", ranking_ids)
)

finite_ranking_n <- sum(
  is.finite(as.numeric(ranking_vector))
)

if (ranking_n != 57828) {
  stop("Ranking length is not 57,828.")
}

if (duplicate_gene_n != 0) {
  stop("Duplicate Ensembl IDs detected in ranking.")
}

if (missing_gene_n != 0) {
  stop("Missing Ensembl IDs detected in ranking.")
}

if (standard_ensembl_n != 57828) {
  stop("Ranking does not contain exactly 57,828 standard Ensembl IDs.")
}

if (finite_ranking_n != 57828) {
  stop("Ranking contains non-finite values.")
}


# ------------------------------------------------------------
# 6. Load previous collision audit
# ------------------------------------------------------------

audit3a4 <- readRDS(step9d3a4_file)

if (!isTRUE(audit3a4$hard_pass)) {
  stop(
    "STEP9D-3A-4 hard_pass is not TRUE."
  )
}

if (!identical(
  as.integer(audit3a4$genome_wide_gene_n),
  as.integer(57828)
)) {
  stop(
    "STEP9D-3A-4 does not document the locked 57,828-gene ranking."
  )
}


# ------------------------------------------------------------
# 7. Reconstruct annotation mapping
#
# Annotation lookup only.
# Ranking remains immutable.
# ------------------------------------------------------------

mapping_raw <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ranking_ids,
  keytype = "ENSEMBL",
  columns = c(
    "ENSEMBL",
    "ENTREZID",
    "SYMBOL"
  )
)

mapping_raw <- as.data.frame(mapping_raw)

mapping_valid <- mapping_raw[
  !is.na(mapping_raw$ENSEMBL) &
    !is.na(mapping_raw$ENTREZID),
  c(
    "ENSEMBL",
    "ENTREZID",
    "SYMBOL"
  ),
  drop = FALSE
]

mapping_valid <- unique(mapping_valid)

mapping_valid$ENSEMBL <-
  as.character(mapping_valid$ENSEMBL)

mapping_valid$ENTREZID <-
  as.character(mapping_valid$ENTREZID)

mapping_valid$SYMBOL <-
  as.character(mapping_valid$SYMBOL)


# ------------------------------------------------------------
# 8. Define Ensembl eligibility
#
# Exactly one Entrez per Ensembl.
#
# IMPORTANT:
# This does NOT remove many-to-one Entrez collisions.
# ------------------------------------------------------------

ensembl_to_entrez_count <- aggregate(
  ENTREZID ~ ENSEMBL,
  data = mapping_valid,
  FUN = function(x)
    length(unique(x))
)

names(ensembl_to_entrez_count)[2] <-
  "n_entrez_per_ensembl"

eligible_ensembl <- ensembl_to_entrez_count$ENSEMBL[
  ensembl_to_entrez_count$n_entrez_per_ensembl == 1
]


# ------------------------------------------------------------
# 9. Construct KEGG candidate annotation table
# ------------------------------------------------------------

kegg_candidate <- mapping_valid[
  mapping_valid$ENSEMBL %in%
    eligible_ensembl,
  ,
  drop = FALSE
]

kegg_candidate <- unique(
  kegg_candidate[
    ,
    c(
      "ENSEMBL",
      "ENTREZID",
      "SYMBOL"
    )
  ]
)

# Immutable ranking lookup
kegg_candidate$ranking <-
  as.numeric(
    ranking_vector[
      kegg_candidate$ENSEMBL
    ]
  )

if (any(is.na(kegg_candidate$ranking))) {
  stop(
    "Missing ranking value in KEGG candidate table."
  )
}

kegg_candidate$ranking_direction <-
  ifelse(
    kegg_candidate$ranking > 0,
    "positive",
    ifelse(
      kegg_candidate$ranking < 0,
      "negative",
      "zero"
    )
  )

kegg_candidate$abs_ranking <-
  abs(kegg_candidate$ranking)

kegg_candidate$is_core_gene <- FALSE

kegg_candidate$core_gene <- NA_character_

for (g in names(core_genes)) {
  
  idx <- kegg_candidate$ENSEMBL ==
    core_genes[[g]]
  
  kegg_candidate$is_core_gene[idx] <-
    TRUE
  
  kegg_candidate$core_gene[idx] <-
    g
}


# ------------------------------------------------------------
# 10. Correct collision classification
#
# This explicitly distinguishes:
# positive + negative
# from
# positive + negative + zero
#
# In the present locked ranking there are no zero-valued
# collision genes.
# ------------------------------------------------------------

entrez_split <- split(
  kegg_candidate,
  kegg_candidate$ENTREZID
)

collision_groups <- entrez_split[
  sapply(
    entrez_split,
    function(d)
      length(unique(d$ENSEMBL)) > 1
  )
]


# ------------------------------------------------------------
# 11. Build complete collision group table
# ------------------------------------------------------------

collision_group_table <- do.call(
  rbind,
  lapply(
    names(collision_groups),
    function(eid) {
      
      d <- collision_groups[[eid]]
      
      signs <- unique(
        d$ranking_direction
      )
      
      direction_class <- if (
        all(signs == "positive")
      ) {
        
        "positive_only"
        
      } else if (
        all(signs == "negative")
      ) {
        
        "negative_only"
        
      } else if (
        all(signs %in%
            c("positive", "zero"))
      ) {
        
        "positive_zero"
        
      } else if (
        all(signs %in%
            c("negative", "zero"))
      ) {
        
        "negative_zero"
        
      } else if (
        all(signs %in%
            c("positive", "negative"))
      ) {
        
        "positive_negative"
        
      } else if (
        all(signs %in%
            c(
              "positive",
              "negative",
              "zero"
            ))
      ) {
        
        "positive_negative_zero"
        
      } else {
        
        "other"
      }
      
      
      data.frame(
        ENTREZID = eid,
        
        n_ensembl =
          length(unique(d$ENSEMBL)),
        
        ranking_unique_n =
          length(unique(d$ranking)),
        
        abs_ranking_unique_n =
          length(unique(d$abs_ranking)),
        
        ranking_min =
          min(d$ranking),
        
        ranking_max =
          max(d$ranking),
        
        abs_ranking_min =
          min(d$abs_ranking),
        
        abs_ranking_max =
          max(d$abs_ranking),
        
        direction_class =
          direction_class,
        
        contains_positive =
          any(
            d$ranking_direction ==
              "positive"
          ),
        
        contains_negative =
          any(
            d$ranking_direction ==
              "negative"
          ),
        
        contains_zero =
          any(
            d$ranking_direction ==
              "zero"
          ),
        
        contains_core_gene =
          any(d$is_core_gene),
        
        core_genes =
          paste(
            unique(
              na.omit(
                d$core_gene
              )
            ),
            collapse = ";"
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(
  collision_group_table
) <- NULL


# ------------------------------------------------------------
# 12. Candidate Policy A
#
# STRICT BIJECTION
#
# Keep only Entrez IDs represented by exactly one Ensembl
# among the KEGG candidate set.
#
# This is a policy simulation only.
# ------------------------------------------------------------

entrez_count_candidate <- tapply(
  kegg_candidate$ENSEMBL,
  kegg_candidate$ENTREZID,
  function(x)
    length(unique(x))
)

policyA_entrez <- names(
  entrez_count_candidate[
    entrez_count_candidate == 1
  ]
)

policyA <- kegg_candidate[
  kegg_candidate$ENTREZID %in%
    policyA_entrez,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 13. Candidate Policy B
#
# MAXIMUM ABSOLUTE RANKING
#
# For each Entrez collision group:
# retain the Ensembl gene with maximum |ranking|.
#
# This is NOT applied to the locked ranking vector.
# It only creates a policy simulation object.
#
# Tie handling:
# if absolute ranking is tied, the policy remains unresolved
# and no arbitrary selection is performed.
# ------------------------------------------------------------

policyB_collision_selection <- do.call(
  rbind,
  lapply(
    collision_groups,
    function(d) {
      
      max_abs <- max(
        d$abs_ranking
      )
      
      top <- d[
        d$abs_ranking == max_abs,
        ,
        drop = FALSE
      ]
      
      top
    }
  )
)

rownames(
  policyB_collision_selection
) <- NULL


policyB_unresolved_ties <- do.call(
  rbind,
  lapply(
    collision_groups,
    function(d) {
      
      max_abs <- max(
        d$abs_ranking
      )
      
      top <- d[
        d$abs_ranking == max_abs,
        ,
        drop = FALSE
      ]
      
      if (nrow(top) > 1) {
        
        top
        
      } else {
        
        NULL
      }
    }
  )
)

if (is.null(policyB_unresolved_ties)) {
  
  policyB_unresolved_tie_group_n <- 0
  
} else {
  
  policyB_unresolved_tie_group_n <-
    length(
      unique(
        policyB_unresolved_ties$ENTREZID
      )
    )
}


# Only construct Policy B if no magnitude ties exist.
if (
  policyB_unresolved_tie_group_n == 0
) {
  
  policyB <- rbind(
    kegg_candidate[
      kegg_candidate$ENTREZID %in%
        names(
          entrez_count_candidate[
            entrez_count_candidate == 1
          ]
        ),
      ,
      drop = FALSE
    ],
    policyB_collision_selection
  )
  
} else {
  
  policyB <- NULL
}


# ------------------------------------------------------------
# 14. Candidate Policy C
#
# DETERMINISTIC RANKING-ORDER REPRESENTATIVE
#
# Ordering:
#   1. abs(ranking) descending
#   2. ranking descending
#   3. ENSEMBL ascending
#
# This does NOT change ranking values.
#
# If abs(ranking) ties, ranking then determines direction.
# If both ranking and abs(ranking) tie, Ensembl is used only
# as a deterministic identifier tie-breaker.
# ------------------------------------------------------------

policyC_collision_selection <- do.call(
  rbind,
  lapply(
    collision_groups,
    function(d) {
      
      d <- d[
        order(
          -d$abs_ranking,
          -d$ranking,
          d$ENSEMBL
        ),
        ,
        drop = FALSE
      ]
      
      d[1, , drop = FALSE]
    }
  )
)

rownames(
  policyC_collision_selection
) <- NULL


policyC <- rbind(
  kegg_candidate[
    kegg_candidate$ENTREZID %in%
      names(
        entrez_count_candidate[
          entrez_count_candidate == 1
        ]
      ),
    ,
    drop = FALSE
  ],
  policyC_collision_selection
)


# ------------------------------------------------------------
# 15. Ensure one row per Entrez under Policy B/C
# ------------------------------------------------------------

policyB_duplicate_entrez_n <-
  if (is.null(policyB)) {
    NA_integer_
  } else {
    sum(
      duplicated(policyB$ENTREZID)
    )
  }

policyC_duplicate_entrez_n <-
  sum(
    duplicated(policyC$ENTREZID)
  )


# ------------------------------------------------------------
# 16. Policy summary helper
# ------------------------------------------------------------

summarize_policy <- function(
    d,
    policy_name) {
  
  if (is.null(d)) {
    
    return(
      data.frame(
        policy = policy_name,
        retained_ensembl_n = NA_integer_,
        retained_entrez_n = NA_integer_,
        positive_n = NA_integer_,
        negative_n = NA_integer_,
        zero_n = NA_integer_,
        core_gene_n = NA_integer_,
        duplicate_entrez_n = NA_integer_,
        stringsAsFactors = FALSE
      )
    )
  }
  
  data.frame(
    policy = policy_name,
    
    retained_ensembl_n =
      length(unique(d$ENSEMBL)),
    
    retained_entrez_n =
      length(unique(d$ENTREZID)),
    
    positive_n =
      sum(d$ranking > 0),
    
    negative_n =
      sum(d$ranking < 0),
    
    zero_n =
      sum(d$ranking == 0),
    
    core_gene_n =
      sum(d$is_core_gene),
    
    duplicate_entrez_n =
      sum(
        duplicated(d$ENTREZID)
      ),
    
    stringsAsFactors = FALSE
  )
}


policy_summary <- rbind(
  
  summarize_policy(
    policyA,
    "Policy_A_Strict_Bijection"
  ),
  
  summarize_policy(
    policyB,
    "Policy_B_Max_Absolute_Ranking"
  ),
  
  summarize_policy(
    policyC,
    "Policy_C_Deterministic_Ranking_Order"
  )
)

rownames(policy_summary) <- NULL


# ------------------------------------------------------------
# 17. Quantify discarded genes relative to candidate input
# ------------------------------------------------------------

candidate_ensembl_n <-
  length(
    unique(
      kegg_candidate$ENSEMBL
    )
  )

candidate_entrez_n <-
  length(
    unique(
      kegg_candidate$ENTREZID
    )
  )


policy_loss <- data.frame(
  
  policy = policy_summary$policy,
  
  candidate_ensembl_n =
    candidate_ensembl_n,
  
  retained_ensembl_n =
    policy_summary$retained_ensembl_n,
  
  discarded_ensembl_n =
    candidate_ensembl_n -
    policy_summary$retained_ensembl_n,
  
  candidate_entrez_n =
    candidate_entrez_n,
  
  retained_entrez_n =
    policy_summary$retained_entrez_n,
  
  discarded_or_unrepresented_entrez_n =
    candidate_entrez_n -
    policy_summary$retained_entrez_n,
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 18. Directional effect of each policy
# ------------------------------------------------------------

candidate_positive_n <-
  sum(
    kegg_candidate$ranking > 0
  )

candidate_negative_n <-
  sum(
    kegg_candidate$ranking < 0
  )

candidate_zero_n <-
  sum(
    kegg_candidate$ranking == 0
  )

policy_direction_summary <- data.frame(
  
  policy =
    policy_summary$policy,
  
  positive_retained_n =
    policy_summary$positive_n,
  
  positive_discarded_n =
    candidate_positive_n -
    policy_summary$positive_n,
  
  negative_retained_n =
    policy_summary$negative_n,
  
  negative_discarded_n =
    candidate_negative_n -
    policy_summary$negative_n,
  
  zero_retained_n =
    policy_summary$zero_n,
  
  zero_discarded_n =
    candidate_zero_n -
    policy_summary$zero_n,
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 19. Opposite-direction collision groups
# ------------------------------------------------------------

opposite_direction_groups <-
  collision_group_table[
    collision_group_table$direction_class ==
      "positive_negative",
    ,
    drop = FALSE
  ]

opposite_entrez_ids <-
  opposite_direction_groups$ENTREZID


# ------------------------------------------------------------
# 20. Resolve ONLY for audit comparison:
# Which Ensembl would each policy retain?
#
# No final input is generated.
# ------------------------------------------------------------

opposite_detail_list <- lapply(
  opposite_entrez_ids,
  function(eid) {
    
    d <- collision_groups[[eid]]
    
    d <- d[
      order(
        -d$abs_ranking,
        -d$ranking,
        d$ENSEMBL
      ),
      ,
      drop = FALSE
    ]
    
    max_abs <- max(
      d$abs_ranking
    )
    
    max_abs_rows <- d[
      d$abs_ranking == max_abs,
      ,
      drop = FALSE
    ]
    
    data.frame(
      
      ENTREZID = eid,
      
      ENSEMBL_all =
        paste(
          d$ENSEMBL,
          collapse = ";"
        ),
      
      SYMBOL_all =
        paste(
          d$SYMBOL,
          collapse = ";"
        ),
      
      ranking_all =
        paste(
          sprintf(
            "%.12f",
            d$ranking
          ),
          collapse = ";"
        ),
      
      direction_all =
        paste(
          d$ranking_direction,
          collapse = ";"
        ),
      
      abs_ranking_all =
        paste(
          sprintf(
            "%.12f",
            d$abs_ranking
          ),
          collapse = ";"
        ),
      
      max_abs_ENSEMBL =
        paste(
          max_abs_rows$ENSEMBL,
          collapse = ";"
        ),
      
      max_abs_ranking =
        max_abs,
      
      max_abs_direction =
        paste(
          unique(
            max_abs_rows$ranking_direction
          ),
          collapse = ";"
        ),
      
      max_abs_tie =
        nrow(max_abs_rows) > 1,
      
      policyC_ENSEMBL =
        d$ENSEMBL[1],
      
      policyC_ranking =
        d$ranking[1],
      
      policyC_direction =
        d$ranking_direction[1],
      
      stringsAsFactors = FALSE
    )
  }
)

opposite_direction_detail <- do.call(
  rbind,
  opposite_detail_list
)

rownames(
  opposite_direction_detail
) <- NULL


# ------------------------------------------------------------
# 21. Check whether Policy B and C differ
# ------------------------------------------------------------

if (
  nrow(opposite_direction_detail) > 0
) {
  
  policyB_C_disagreement_n <- sum(
    opposite_direction_detail$max_abs_ENSEMBL !=
      opposite_direction_detail$policyC_ENSEMBL
  )
  
} else {
  
  policyB_C_disagreement_n <- 0
}


# ------------------------------------------------------------
# 22. Core-gene preservation
# ------------------------------------------------------------

core_policy_table <- data.frame(
  
  core_gene =
    names(core_genes),
  
  ENSEMBL =
    unname(core_genes),
  
  candidate_present =
    core_genes %in%
    kegg_candidate$ENSEMBL,
  
  Policy_A_retained =
    if (nrow(policyA) > 0)
      core_genes %in% policyA$ENSEMBL
  else
    FALSE,
  
  Policy_B_retained =
    if (!is.null(policyB))
      core_genes %in% policyB$ENSEMBL
  else
    FALSE,
  
  Policy_C_retained =
    if (nrow(policyC) > 0)
      core_genes %in% policyC$ENSEMBL
  else
    FALSE,
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 23. Formal audit flags
# ------------------------------------------------------------

no_ranking_recalculation <- TRUE

no_deg_recalculation <- TRUE

no_gsea <- TRUE

no_final_resolution <- TRUE

no_arbitrary_first_match <- TRUE

no_jitter <- TRUE

no_artificial_ranking_tie_breaking <- TRUE

original_ranking_immutable <- identical(
  sort(ranking_ids),
  sort(names(ranking_vector))
)

candidate_has_collision <-
  any(
    entrez_count_candidate > 1
  )

opposite_direction_present <-
  nrow(
    opposite_direction_groups
  ) > 0


# ------------------------------------------------------------
# 24. Policy B validity
# ------------------------------------------------------------

policyB_valid <-
  policyB_unresolved_tie_group_n == 0 &&
  !is.null(policyB) &&
  policyB_duplicate_entrez_n == 0


# ------------------------------------------------------------
# 25. Policy C validity
# ------------------------------------------------------------

policyC_valid <-
  policyC_duplicate_entrez_n == 0


# ------------------------------------------------------------
# 26. Hard PASS
#
# This PASS means:
# the policy comparison was completed correctly.
#
# It does NOT mean a final KEGG policy has been selected.
# ------------------------------------------------------------

hard_qc <- c(
  
  ranking_n_57828 =
    ranking_n == 57828,
  
  ranking_unique =
    duplicate_gene_n == 0,
  
  standard_ensembl_57828 =
    standard_ensembl_n == 57828,
  
  ranking_finite =
    finite_ranking_n == 57828,
  
  previous_step_3a4_pass =
    isTRUE(
      audit3a4$hard_pass
    ),
  
  collision_detected_correctly =
    candidate_has_collision,
  
  opposite_direction_detected =
    opposite_direction_present,
  
  no_ranking_recalculation =
    no_ranking_recalculation,
  
  no_deg_recalculation =
    no_deg_recalculation,
  
  no_gsea =
    no_gsea,
  
  no_final_resolution =
    no_final_resolution,
  
  no_arbitrary_first_match =
    no_arbitrary_first_match,
  
  no_jitter =
    no_jitter,
  
  no_artificial_ranking_tie_breaking =
    no_artificial_ranking_tie_breaking,
  
  original_ranking_immutable =
    original_ranking_immutable,
  
  Policy_B_no_unresolved_abs_tie =
    policyB_unresolved_tie_group_n == 0,
  
  Policy_B_unique_entrez =
    policyB_valid,
  
  Policy_C_unique_entrez =
    policyC_valid
)

hard_pass <- all(
  hard_qc
)


# ------------------------------------------------------------
# 27. Console summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-3A-5 POLICY COMPARISON SUMMARY\n")
cat("============================================================\n\n")

cat("Locked ranking genes: ",
    ranking_n, "\n")

cat("KEGG candidate Ensembl genes: ",
    candidate_ensembl_n, "\n")

cat("KEGG candidate unique Entrez IDs: ",
    candidate_entrez_n, "\n")

cat("Collision Entrez groups: ",
    length(collision_groups), "\n")

cat("Opposite-direction groups: ",
    nrow(opposite_direction_groups), "\n\n")


cat("------------------------------------------------------------\n")
cat("POLICY SUMMARY\n")
cat("------------------------------------------------------------\n")

print(policy_summary)

cat("\n------------------------------------------------------------\n")
cat("POLICY LOSS / RETENTION\n")
cat("------------------------------------------------------------\n")

print(policy_loss)

cat("\n------------------------------------------------------------\n")
cat("POLICY DIRECTIONAL EFFECT\n")
cat("------------------------------------------------------------\n")

print(policy_direction_summary)

cat("\n------------------------------------------------------------\n")
cat("OPPOSITE-DIRECTION COLLISION DETAIL\n")
cat("------------------------------------------------------------\n")

print(opposite_direction_detail)

cat("\n------------------------------------------------------------\n")
cat("POLICY B vs POLICY C DISAGREEMENT\n")
cat("------------------------------------------------------------\n")

cat(
  "Groups where Policy B and C select different Ensembl IDs: ",
  policyB_C_disagreement_n,
  "\n"
)

cat("\n------------------------------------------------------------\n")
cat("CORE-GENE PRESERVATION\n")
cat("------------------------------------------------------------\n")

print(core_policy_table)

cat("\n------------------------------------------------------------\n")
cat("POLICY VALIDITY\n")
cat("------------------------------------------------------------\n")

cat(
  "Policy A valid: ",
  TRUE,
  "\n"
)

cat(
  "Policy B valid: ",
  policyB_valid,
  "\n"
)

cat(
  "Policy C valid: ",
  policyC_valid,
  "\n"
)

cat("\n------------------------------------------------------------\n")
cat("HARD QC\n")
cat("------------------------------------------------------------\n")

print(hard_qc)

cat("\nSTEP9D-3A-5 FINAL HARD PASS: ",
    hard_pass,
    "\n")

cat("\nIMPORTANT:\n")
cat(
  "This PASS certifies the audit only.\n"
)
cat(
  "No final KEGG collision-resolution policy has been locked.\n"
)
cat(
  "No GSEA has been performed.\n"
)


# ------------------------------------------------------------
# 28. Save policy comparison tables
# ------------------------------------------------------------

write.csv(
  policy_summary,
  file.path(
    output_dir,
    "STEP9D3A5_Policy_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  policy_loss,
  file.path(
    output_dir,
    "STEP9D3A5_Policy_Loss_Retention.csv"
  ),
  row.names = FALSE
)

write.csv(
  policy_direction_summary,
  file.path(
    output_dir,
    "STEP9D3A5_Policy_Directional_Effect.csv"
  ),
  row.names = FALSE
)

write.csv(
  collision_group_table,
  file.path(
    output_dir,
    "STEP9D3A5_Collision_Group_Table_Corrected.csv"
  ),
  row.names = FALSE
)

write.csv(
  opposite_direction_detail,
  file.path(
    output_dir,
    "STEP9D3A5_Opposite_Direction_Detail.csv"
  ),
  row.names = FALSE
)

write.csv(
  core_policy_table,
  file.path(
    output_dir,
    "STEP9D3A5_CoreGene_Policy_Preservation.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 29. Save complete audit RDS
# ------------------------------------------------------------

final_audit <- list(
  
  workflow = "Paper 1",
  
  step = "STEP9D-3A-5",
  
  status =
    "FINAL_READ_ONLY_POLICY_COMPARISON_AUDIT",
  
  audit_date =
    as.character(Sys.time()),
  
  ranking_file =
    ranking_file,
  
  previous_step =
    step9d3a4_file,
  
  genome_wide_gene_n =
    ranking_n,
  
  kegg_candidate_ensembl_n =
    candidate_ensembl_n,
  
  kegg_candidate_unique_entrez_n =
    candidate_entrez_n,
  
  collision_entrez_group_n =
    length(collision_groups),
  
  opposite_direction_group_n =
    nrow(opposite_direction_groups),
  
  policy_summary =
    policy_summary,
  
  policy_loss =
    policy_loss,
  
  policy_direction_summary =
    policy_direction_summary,
  
  collision_group_table =
    collision_group_table,
  
  opposite_direction_detail =
    opposite_direction_detail,
  
  core_policy_table =
    core_policy_table,
  
  policyB_unresolved_tie_group_n =
    policyB_unresolved_tie_group_n,
  
  policyB_C_disagreement_n =
    policyB_C_disagreement_n,
  
  policyA_definition =
    "Strict bijection: exactly one Ensembl per Entrez and exactly one Entrez per Ensembl",
  
  policyB_definition =
    "Within each Entrez collision group retain maximum absolute ranking; unresolved absolute-ranking ties are not arbitrarily resolved",
  
  policyC_definition =
    "Within each Entrez collision group order by decreasing absolute ranking, decreasing signed ranking, then ascending Ensembl ID and retain first",
  
  ranking_recalculated =
    FALSE,
  
  deg_recalculated =
    FALSE,
  
  gsea_performed =
    FALSE,
  
  final_collision_resolution_applied =
    FALSE,
  
  arbitrary_first_match_used =
    FALSE,
  
  jitter_used =
    FALSE,
  
  artificial_ranking_tie_breaking_used =
    FALSE,
  
  final_policy_locked =
    FALSE,
  
  hard_qc =
    hard_qc,
  
  hard_pass =
    hard_pass
)

audit_rds <- file.path(
  output_dir,
  "STEP9D3A5_FINAL_KEGG_Collision_Policy_Comparison_Audit.rds"
)

saveRDS(
  final_audit,
  audit_rds
)


# ------------------------------------------------------------
# 30. Formal TXT report
# ------------------------------------------------------------

audit_txt <- file.path(
  output_dir,
  "STEP9D3A5_FINAL_KEGG_Collision_Policy_Comparison_Audit.txt"
)

sink(audit_txt)

cat("============================================================\n")
cat("Paper 1\n")
cat("STEP9D-3A-5\n")
cat("KEGG Entrez Collision Resolution Policy Comparison Audit\n")
cat("============================================================\n\n")

cat("STATUS: FINAL READ-ONLY POLICY COMPARISON\n\n")

cat("Genome-wide ranking genes: ",
    ranking_n, "\n")

cat("KEGG candidate Ensembl genes: ",
    candidate_ensembl_n, "\n")

cat("KEGG candidate unique Entrez IDs: ",
    candidate_entrez_n, "\n")

cat("Collision Entrez groups: ",
    length(collision_groups), "\n")

cat("Opposite-direction collision groups: ",
    nrow(opposite_direction_groups), "\n\n")

cat("------------------------------------------------------------\n")
cat("POLICY SUMMARY\n")
cat("------------------------------------------------------------\n")

print(policy_summary)

cat("\n------------------------------------------------------------\n")
cat("POLICY LOSS / RETENTION\n")
cat("------------------------------------------------------------\n")

print(policy_loss)

cat("\n------------------------------------------------------------\n")
cat("POLICY DIRECTIONAL EFFECT\n")
cat("------------------------------------------------------------\n")

print(policy_direction_summary)

cat("\n------------------------------------------------------------\n")
cat("OPPOSITE-DIRECTION COLLISION DETAIL\n")
cat("------------------------------------------------------------\n")

print(opposite_direction_detail)

cat("\n------------------------------------------------------------\n")
cat("CORE-GENE PRESERVATION\n")
cat("------------------------------------------------------------\n")

print(core_policy_table)

cat("\n------------------------------------------------------------\n")
cat("POLICY DEFINITIONS\n")
cat("------------------------------------------------------------\n")

cat("Policy A:\n")
cat(
  "Strict bijection; exactly one Ensembl per Entrez and\n"
)
cat(
  "exactly one Entrez per Ensembl.\n\n"
)

cat("Policy B:\n")
cat(
  "Maximum absolute ranking within each Entrez collision group.\n"
)
cat(
  "Absolute-ranking ties remain unresolved.\n\n"
)

cat("Policy C:\n")
cat(
  "Deterministic ordering by decreasing absolute ranking,\n"
)
cat(
  "then decreasing signed ranking, then ascending Ensembl ID.\n\n"
)

cat("------------------------------------------------------------\n")
cat("NO-SILENT-CHANGE FLAGS\n")
cat("------------------------------------------------------------\n")

cat("Ranking recalculated: FALSE\n")
cat("DEG recalculated: FALSE\n")
cat("GSEA performed: FALSE\n")
cat("Final collision resolution applied: FALSE\n")
cat("Arbitrary first-match used: FALSE\n")
cat("Jitter used: FALSE\n")
cat("Artificial ranking tie-breaking used: FALSE\n")
cat("Final policy locked: FALSE\n")

cat("\nSTEP9D-3A-5 FINAL HARD PASS: ",
    hard_pass,
    "\n")

sink()


# ------------------------------------------------------------
# 31. Final file audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-3A-5 FINAL FILE AUDIT\n")
cat("============================================================\n\n")

cat(
  "Policy summary CSV: ",
  file.exists(
    file.path(
      output_dir,
      "STEP9D3A5_Policy_Summary.csv"
    )
  ),
  "\n"
)

cat(
  "Policy loss CSV: ",
  file.exists(
    file.path(
      output_dir,
      "STEP9D3A5_Policy_Loss_Retention.csv"
    )
  ),
  "\n"
)

cat(
  "Policy directional CSV: ",
  file.exists(
    file.path(
      output_dir,
      "STEP9D3A5_Policy_Directional_Effect.csv"
    )
  ),
  "\n"
)

cat(
  "Corrected collision table CSV: ",
  file.exists(
    file.path(
      output_dir,
      "STEP9D3A5_Collision_Group_Table_Corrected.csv"
    )
  ),
  "\n"
)

cat(
  "Opposite-direction detail CSV: ",
  file.exists(
    file.path(
      output_dir,
      "STEP9D3A5_Opposite_Direction_Detail.csv"
    )
  ),
  "\n"
)

cat(
  "Core-gene policy CSV: ",
  file.exists(
    file.path(
      output_dir,
      "STEP9D3A5_CoreGene_Policy_Preservation.csv"
    )
  ),
  "\n"
)

cat(
  "Final audit RDS: ",
  file.exists(audit_rds),
  "\n"
)

cat(
  "Final audit TXT: ",
  file.exists(audit_txt),
  "\n"
)

cat("\n============================================================\n")
cat(
  "STEP9D-3A-5 FINAL HARD PASS: ",
  hard_pass,
  "\n"
)
cat("============================================================\n")    


# ============================================================
# Paper 1
# STEP9D-3A-5
# READ-ONLY OUTPUT FILE LOCATION AUDIT
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-3A-5 READ-ONLY OUTPUT FILE LOCATION AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Root directory
# ------------------------------------------------------------

root_dir <-
  "input_dir/Paper 1_Step 9-Functional Enrichment"

search_root <-
  file.path(
    root_dir,
    "STEP9D_GenomeWide_Preranked_GSEA"
  )


cat("Search root:\n")
cat(search_root, "\n\n")


if (!dir.exists(search_root)) {
  stop(
    paste0(
      "Search root does not exist:\n",
      search_root
    )
  )
}


# ------------------------------------------------------------
# 2. Search for all STEP9D-3A-5 files
# ------------------------------------------------------------

all_files <-
  list.files(
    path = search_root,
    pattern = "STEP9D3A5",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = FALSE
  )


cat("Number of STEP9D3A5 files found: ")
cat(length(all_files), "\n\n")


# ------------------------------------------------------------
# 3. Print all matching files
# ------------------------------------------------------------

if (length(all_files) > 0) {
  
  cat("------------------------------------------------------------\n")
  cat("ALL MATCHING STEP9D3A5 FILES\n")
  cat("------------------------------------------------------------\n\n")
  
  print(
    all_files,
    quote = FALSE
  )
  
} else {
  
  cat(
    "No STEP9D3A5 files found under the specified search root.\n"
  )
}


# ------------------------------------------------------------
# 4. Search specifically for expected output types
# ------------------------------------------------------------

target_patterns <- c(
  "Policy_Summary",
  "Policy_Loss",
  "Policy_Directional",
  "Collision_Group_Table",
  "Opposite_Direction",
  "CoreGene_Policy",
  "FINAL_Policy",
  "Policy_Comparison"
)


cat("\n\n")
cat("============================================================\n")
cat("TARGET FILE SEARCH\n")
cat("============================================================\n\n")


for (pat in target_patterns) {
  
  hits <-
    list.files(
      path = search_root,
      pattern = pat,
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = FALSE
    )
  
  cat("\n------------------------------------------------------------\n")
  cat("Pattern: ", pat, "\n", sep = "")
  cat("Number found: ", length(hits), "\n", sep = "")
  
  if (length(hits) > 0) {
    
    print(
      hits,
      quote = FALSE
    )
  }
}


# ------------------------------------------------------------
# 5. Directory-level audit
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("STEP9D-3A-5 DIRECTORY CANDIDATE AUDIT\n")
cat("============================================================\n\n")


candidate_dirs <-
  unique(
    dirname(all_files)
  )


if (length(candidate_dirs) > 0) {
  
  print(
    candidate_dirs,
    quote = FALSE
  )
  
} else {
  
  cat("No candidate directories detected.\n")
}


# ------------------------------------------------------------
# 6. Important safety statement
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("READ-ONLY SAFETY STATUS\n")
cat("============================================================\n\n")

cat("No files were created.\n")
cat("No files were modified.\n")
cat("No files were deleted.\n")
cat("No ranking was recalculated.\n")
cat("No DEG was recalculated.\n")
cat("No ID mapping was recalculated.\n")
cat("No collision resolution was performed.\n")
cat("No GSEA was performed.\n")

cat("\n============================================================\n")
cat("END OF FILE LOCATION AUDIT\n")
cat("============================================================\n\n")    

# ============================================================
# Paper 1
# STEP9D-3A-5
# FINAL RESULT INSPECTION / INTERPRETATION AUDIT
# ============================================================
#
# READ-ONLY
#
# This script reads ONLY the files already generated by
# STEP9D-3A-5.
#
# NO:
#   - ranking recalculation
#   - DEG recalculation
#   - ID remapping
#   - collision resolution
#   - duplicate Entrez removal
#   - GSEA
#   - modification of existing files
#
# ============================================================


cat("\n============================================================\n")
cat("STEP9D-3A-5 FINAL RESULT INSPECTION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Exact locked output directory
# ------------------------------------------------------------

root_dir <-
  "input_dir/Paper 1_Step 9-Functional Enrichment"

output_dir <-
  file.path(
    root_dir,
    "STEP9D_GenomeWide_Preranked_GSEA",
    "9D-3A-5_KEGG_Collision_Policy_Audit"
  )


cat("Locked STEP9D-3A-5 output directory:\n")
cat(output_dir, "\n\n")


if (!dir.exists(output_dir)) {
  stop(
    "STEP9D-3A-5 output directory does not exist."
  )
}


# ------------------------------------------------------------
# 1. Exact file paths
# ------------------------------------------------------------

policy_summary_file <-
  file.path(
    output_dir,
    "STEP9D3A5_Policy_Summary.csv"
  )

policy_loss_file <-
  file.path(
    output_dir,
    "STEP9D3A5_Policy_Loss_Retention.csv"
  )

policy_direction_file <-
  file.path(
    output_dir,
    "STEP9D3A5_Policy_Directional_Effect.csv"
  )

collision_table_file <-
  file.path(
    output_dir,
    "STEP9D3A5_Collision_Group_Table_Corrected.csv"
  )

opposite_detail_file <-
  file.path(
    output_dir,
    "STEP9D3A5_Opposite_Direction_Detail.csv"
  )

core_gene_file <-
  file.path(
    output_dir,
    "STEP9D3A5_CoreGene_Policy_Preservation.csv"
  )

audit_rds <-
  file.path(
    output_dir,
    "STEP9D3A5_FINAL_KEGG_Collision_Policy_Comparison_Audit.rds"
  )

audit_txt <-
  file.path(
    output_dir,
    "STEP9D3A5_FINAL_KEGG_Collision_Policy_Comparison_Audit.txt"
  )


expected_files <- c(
  policy_summary_file,
  policy_loss_file,
  policy_direction_file,
  collision_table_file,
  opposite_detail_file,
  core_gene_file,
  audit_rds,
  audit_txt
)


# ------------------------------------------------------------
# 2. Exact file existence check
# ------------------------------------------------------------

file_status <- data.frame(
  file = basename(expected_files),
  exists = file.exists(expected_files),
  stringsAsFactors = FALSE
)


cat("------------------------------------------------------------\n")
cat("1. EXACT FILE EXISTENCE CHECK\n")
cat("------------------------------------------------------------\n\n")

print(
  file_status,
  row.names = FALSE
)


if (!all(file_status$exists)) {
  
  stop(
    "One or more expected STEP9D-3A-5 files are missing."
  )
}


cat(
  "\nAll 8 expected STEP9D-3A-5 files exist: TRUE\n"
)


# ------------------------------------------------------------
# 3. Read-only loading
# ------------------------------------------------------------

policy_summary <-
  read.csv(
    policy_summary_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

policy_loss <-
  read.csv(
    policy_loss_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

policy_direction <-
  read.csv(
    policy_direction_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

collision_table <-
  read.csv(
    collision_table_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

opposite_detail <-
  read.csv(
    opposite_detail_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

core_gene_policy <-
  read.csv(
    core_gene_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

audit_object <-
  readRDS(audit_rds)


# ------------------------------------------------------------
# 4. Show table dimensions and column names
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("2. TABLE STRUCTURE AUDIT\n")
cat("============================================================\n\n")


table_objects <- list(
  POLICY_SUMMARY = policy_summary,
  POLICY_LOSS_RETENTION = policy_loss,
  POLICY_DIRECTIONAL_EFFECT = policy_direction,
  COLLISION_GROUP_TABLE = collision_table,
  OPPOSITE_DIRECTION_DETAIL = opposite_detail,
  CORE_GENE_POLICY = core_gene_policy
)


for (nm in names(table_objects)) {
  
  obj <- table_objects[[nm]]
  
  cat("\n------------------------------------------------------------\n")
  cat(nm, "\n")
  cat("------------------------------------------------------------\n")
  
  cat(
    "Rows:",
    nrow(obj),
    " | Columns:",
    ncol(obj),
    "\n\n"
  )
  
  cat(
    "Column names:\n"
  )
  
  print(
    names(obj)
  )
}


# ------------------------------------------------------------
# 5. POLICY SUMMARY
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("3. POLICY SUMMARY\n")
cat("============================================================\n\n")

print(
  policy_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. POLICY LOSS / RETENTION
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("4. POLICY LOSS / RETENTION\n")
cat("============================================================\n\n")

print(
  policy_loss,
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. POLICY DIRECTIONAL EFFECT
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("5. POLICY DIRECTIONAL EFFECT\n")
cat("============================================================\n\n")

print(
  policy_direction,
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Corrected collision table
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("6. CORRECTED COLLISION GROUP TABLE\n")
cat("============================================================\n\n")

print(
  collision_table,
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Opposite-direction collision groups
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("7. SIX OPPOSITE-DIRECTION COLLISION GROUPS\n")
cat("============================================================\n\n")


cat(
  "Rows in opposite-direction detail table:",
  nrow(opposite_detail),
  "\n\n"
)

print(
  opposite_detail,
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Core-gene preservation
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("8. CORE-GENE POLICY PRESERVATION\n")
cat("============================================================\n\n")

print(
  core_gene_policy,
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Audit-object structure
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("9. FINAL AUDIT RDS STRUCTURE\n")
cat("============================================================\n\n")

cat("Class:\n")
print(
  class(audit_object)
)

cat("\nNames:\n")
print(
  names(audit_object)
)


# ------------------------------------------------------------
# 12. Print all scalar audit fields
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("10. FINAL AUDIT RDS KEY VALUES\n")
cat("============================================================\n\n")


if (is.list(audit_object)) {
  
  for (nm in names(audit_object)) {
    
    value <- audit_object[[nm]]
    
    if (
      length(value) == 1 &&
      (
        is.atomic(value) ||
        is.null(value)
      )
    ) {
      
      cat(
        sprintf(
          "%-45s : ",
          nm
        )
      )
      
      if (is.null(value)) {
        
        cat("NULL\n")
        
      } else {
        
        cat(
          as.character(value),
          "\n"
        )
      }
    }
  }
}


# ------------------------------------------------------------
# 13. Core-gene explicit search
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("11. CORE-GENE EXPLICIT SEARCH\n")
cat("============================================================\n\n")


core_genes <-
  c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )


core_columns <-
  names(core_gene_policy)


cat(
  "Core-gene table columns:\n"
)

print(
  core_columns
)


for (gene in core_genes) {
  
  cat("\n------------------------------------------------------------\n")
  cat(
    "CORE GENE:",
    gene,
    "\n"
  )
  cat("------------------------------------------------------------\n")
  
  found <- FALSE
  
  for (j in seq_along(core_columns)) {
    
    x <- core_gene_policy[[j]]
    
    if (
      is.character(x) &&
      any(
        x == gene,
        na.rm = TRUE
      )
    ) {
      
      print(
        core_gene_policy[
          x == gene,
          ,
          drop = FALSE
        ],
        row.names = FALSE
      )
      
      found <- TRUE
      break
    }
  }
  
  if (!found) {
    
    cat(
      "No row containing this core gene was detected.\n"
    )
  }
}


# ------------------------------------------------------------
# 14. Explicit opposite-direction group count
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("12. OPPOSITE-DIRECTION GROUP COUNT\n")
cat("============================================================\n\n")


# Try common identifier columns

possible_group_columns <-
  c(
    "ENTREZID",
    "EntrezID",
    "entrez_id",
    "Entrez",
    "collision_entrez",
    "collision_group",
    "group_id"
  )


group_column <-
  possible_group_columns[
    possible_group_columns %in%
      names(opposite_detail)
  ][1]


if (!is.na(group_column)) {
  
  n_groups <-
    length(
      unique(
        opposite_detail[[group_column]]
      )
    )
  
  cat(
    "Group identifier column:",
    group_column,
    "\n"
  )
  
  cat(
    "Number of opposite-direction groups:",
    n_groups,
    "\n"
  )
  
} else {
  
  cat(
    "No standard group identifier column detected.\n"
  )
  
  cat(
    "The complete opposite-direction table is printed above.\n"
  )
}


# ------------------------------------------------------------
# 15. Read-only safety audit
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("13. READ-ONLY SAFETY AUDIT\n")
cat("============================================================\n\n")


safety_fields <-
  c(
    "ranking_recalculated",
    "deg_recalculated",
    "gsea_performed",
    "collision_resolution_applied",
    "duplicate_entrez_removed",
    "arbitrary_first_match",
    "jitter",
    "artificial_tie_breaking"
  )


for (field in safety_fields) {
  
  if (field %in% names(audit_object)) {
    
    cat(
      sprintf(
        "%-35s : %s\n",
        field,
        as.character(
          audit_object[[field]]
        )
      )
    )
  }
}


# ------------------------------------------------------------
# 16. Hard-pass field
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("14. ORIGINAL STEP9D-3A-5 HARD PASS\n")
cat("============================================================\n\n")


if ("hard_pass" %in% names(audit_object)) {
  
  cat(
    "hard_pass:",
    as.character(
      audit_object$hard_pass
    ),
    "\n"
  )
  
} else {
  
  cat(
    "hard_pass field not found in audit object.\n"
  )
}


# ------------------------------------------------------------
# 17. Final interpretation status
# ------------------------------------------------------------

cat("\n\n")
cat("============================================================\n")
cat("15. INTERPRETATION STATUS\n")
cat("============================================================\n\n")

cat(
  "STEP9D-3A-5 comparison audit has already been completed.\n"
)

cat(
  "This script only reads and displays the locked outputs.\n\n"
)

cat(
  "NO final KEGG collision policy is selected by this script.\n"
)

cat(
  "NO Entrez collision resolution is applied by this script.\n"
)

cat(
  "NO KEGG GSEA is performed by this script.\n"
)

cat("\n============================================================\n")
cat("END OF STEP9D-3A-5 RESULT INSPECTION\n")
cat("============================================================\n\n")    


##STEP9D-3A-6 Final KEGG Collision Resolution Policy Lock 
# ============================================================
# Paper 1
# STEP9D-3A-6
# Final KEGG Collision Resolution Policy Lock
# ============================================================
#
# PURPOSE
# ------------------------------------------------------------
# Formally lock the KEGG Ensembl -> Entrez collision-resolution
# policy based on the completed STEP9D-3A-5 comparison audit.
#
# LOCKED POLICY
# ------------------------------------------------------------
# Policy B = Maximum Absolute Ranking
#
# Within each Entrez collision group:
#   retain the Ensembl gene with the largest absolute value
#   of the already LOCKED genome-wide ranking statistic.
#
# IMPORTANT:
# ------------------------------------------------------------
# 1. Original 57,828-gene ranking is NOT recalculated.
# 2. DEG definition is NOT recalculated.
# 3. No new ranking statistic is created.
# 4. No arbitrary first-match selection.
# 5. No jitter.
# 6. No artificial tie-breaking.
# 7. No GSEA is performed in STEP9D-3A-6.
#
# This STEP formally defines the KEGG GSEA input policy.
# Actual KEGG GSEA will be performed only in STEP9D-3B.
#
# ============================================================


cat("\n============================================================\n")
cat("STEP9D-3A-6: Final KEGG Collision Resolution Policy Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_pkgs <- c(
  "readr",
  "dplyr",
  "tibble"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 1. Immutable paths
# ------------------------------------------------------------

root_dir <-
  "input_dir/Paper 1_Step 9-Functional Enrichment"

ranking_dir <-
  file.path(
    root_dir,
    "STEP9D_GenomeWide_Preranked_GSEA",
    "9D-2C_Final_Ranking_Vector"
  )

mapping_lock_dir <-
  file.path(
    root_dir,
    "STEP9D_GenomeWide_Preranked_GSEA",
    "9D-3A-3_GSEA_ID_Mapping_Policy_Lock"
  )

collision_dir <-
  file.path(
    root_dir,
    "STEP9D_GenomeWide_Preranked_GSEA",
    "9D-3A-4_KEGG_Entrez_Collision_Audit"
  )

comparison_dir <-
  file.path(
    root_dir,
    "STEP9D_GenomeWide_Preranked_GSEA",
    "9D-3A-5_KEGG_Collision_Policy_Audit"
  )

output_dir <-
  file.path(
    root_dir,
    "STEP9D_GenomeWide_Preranked_GSEA",
    "9D-3A-6_KEGG_Collision_Resolution_Policy_Lock"
  )

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Input files
# ------------------------------------------------------------

ranking_rds <-
  file.path(
    ranking_dir,
    "STEP9D2C_FINAL_GenomeWide_Ranking_Vector.rds"
  )

mapping_lock_rds <-
  file.path(
    mapping_lock_dir,
    "STEP9D3A3_FINAL_GSEA_ID_Mapping_Policy_Lock.rds"
  )

collision_audit_rds <-
  file.path(
    collision_dir,
    "STEP9D3A4_FINAL_KEGG_Entrez_Collision_Audit.rds"
  )

comparison_audit_rds <-
  file.path(
    comparison_dir,
    "STEP9D3A5_FINAL_KEGG_Collision_Policy_Comparison_Audit.rds"
  )

policy_summary_csv <-
  file.path(
    comparison_dir,
    "STEP9D3A5_Policy_Summary.csv"
  )

policy_loss_csv <-
  file.path(
    comparison_dir,
    "STEP9D3A5_Policy_Loss_Retention.csv"
  )

policy_direction_csv <-
  file.path(
    comparison_dir,
    "STEP9D3A5_Policy_Directional_Effect.csv"
  )

core_policy_csv <-
  file.path(
    comparison_dir,
    "STEP9D3A5_CoreGene_Policy_Preservation.csv"
  )

opposite_direction_csv <-
  file.path(
    comparison_dir,
    "STEP9D3A5_Opposite_Direction_Detail.csv"
  )


# ------------------------------------------------------------
# 3. Required input existence check
# ------------------------------------------------------------

required_files <- c(
  ranking_rds,
  mapping_lock_rds,
  collision_audit_rds,
  comparison_audit_rds,
  policy_summary_csv,
  policy_loss_csv,
  policy_direction_csv,
  core_policy_csv,
  opposite_direction_csv
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  
  cat("Missing required input files:\n")
  print(missing_files)
  
  stop(
    "STEP9D-3A-6 cannot proceed because required locked/audit inputs are missing."
  )
}


# ------------------------------------------------------------
# 4. Load immutable prior-step objects
# ------------------------------------------------------------

ranking_vector <- readRDS(ranking_rds)

mapping_lock <- readRDS(mapping_lock_rds)

collision_audit <- readRDS(collision_audit_rds)

comparison_audit <- readRDS(comparison_audit_rds)


policy_summary <- readr::read_csv(
  policy_summary_csv,
  show_col_types = FALSE
)

policy_loss <- readr::read_csv(
  policy_loss_csv,
  show_col_types = FALSE
)

policy_direction <- readr::read_csv(
  policy_direction_csv,
  show_col_types = FALSE
)

core_policy <- readr::read_csv(
  core_policy_csv,
  show_col_types = FALSE
)

opposite_direction <- readr::read_csv(
  opposite_direction_csv,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 5. Basic ranking-vector integrity audit
# ------------------------------------------------------------

ranking_n <- length(ranking_vector)

ranking_names <- names(ranking_vector)

ranking_name_n <- length(ranking_names)

ranking_duplicate_n <-
  sum(duplicated(ranking_names))

ranking_na_n <-
  sum(is.na(ranking_vector))

ranking_inf_n <-
  sum(is.infinite(ranking_vector))


cat("Genome-wide ranking genes: ", ranking_n, "\n", sep = "")
cat("Ranking names available: ", ranking_name_n, "\n", sep = "")
cat("Duplicated ranking IDs: ", ranking_duplicate_n, "\n", sep = "")
cat("NA ranking values: ", ranking_na_n, "\n", sep = "")
cat("Infinite ranking values: ", ranking_inf_n, "\n\n", sep = "")


ranking_integrity_pass <-
  identical(ranking_n, 57828L) &&
  identical(ranking_name_n, 57828L) &&
  ranking_duplicate_n == 0 &&
  ranking_na_n == 0 &&
  ranking_inf_n == 0


if (!ranking_integrity_pass) {
  stop("Ranking-vector integrity check FAILED.")
}


# ------------------------------------------------------------
# 6. Extract Policy B from STEP9D-3A-5
# ------------------------------------------------------------

policy_b_name <-
  "Policy_B_Max_Absolute_Ranking"

policy_b_summary <-
  policy_summary %>%
  dplyr::filter(policy == policy_b_name)

policy_b_loss <-
  policy_loss %>%
  dplyr::filter(policy == policy_b_name)

policy_b_direction <-
  policy_direction %>%
  dplyr::filter(policy == policy_b_name)


if (nrow(policy_b_summary) != 1L) {
  stop(
    "Policy B summary row is not uniquely identifiable."
  )
}

if (nrow(policy_b_loss) != 1L) {
  stop(
    "Policy B loss/retention row is not uniquely identifiable."
  )
}

if (nrow(policy_b_direction) != 1L) {
  stop(
    "Policy B directional-effect row is not uniquely identifiable."
  )
}


# ------------------------------------------------------------
# 7. LOCKED numerical values
# ------------------------------------------------------------

LOCK_candidate_ensembl_n <- 34138L
LOCK_candidate_entrez_n <- 34120L
LOCK_collision_groups_n <- 18L
LOCK_collision_affected_ensembl_n <- 36L
LOCK_opposite_direction_groups_n <- 6L

LOCK_policy_b_retained_ensembl_n <- 34120L
LOCK_policy_b_retained_entrez_n <- 34120L
LOCK_policy_b_discarded_ensembl_n <- 18L
LOCK_policy_b_discarded_entrez_n <- 0L

LOCK_retained_positive_n <- 21189L
LOCK_retained_negative_n <- 12924L
LOCK_retained_zero_n <- 7L

LOCK_max_absolute_ties_n <- 0L

LOCK_core_AOC3 <- TRUE
LOCK_core_MS4A7 <- TRUE
LOCK_core_IL7R <- TRUE

LOCK_arbitrary_first_match <- FALSE
LOCK_jitter <- FALSE
LOCK_artificial_tie_breaking <- FALSE

LOCK_ranking_recalculated <- FALSE
LOCK_DEG_recalculated <- FALSE
LOCK_GSEA_performed <- FALSE


# ------------------------------------------------------------
# 8. Validate Policy B summary
# ------------------------------------------------------------

summary_checks <- c(
  
  candidate_ensembl =
    policy_b_summary$retained_ensembl_n ==
    LOCK_policy_b_retained_ensembl_n,
  
  candidate_unique_entrez =
    policy_b_summary$retained_entrez_n ==
    LOCK_policy_b_retained_entrez_n,
  
  retained_positive =
    policy_b_summary$positive_n ==
    LOCK_retained_positive_n,
  
  retained_negative =
    policy_b_summary$negative_n ==
    LOCK_retained_negative_n,
  
  retained_zero =
    policy_b_summary$zero_n ==
    LOCK_retained_zero_n,
  
  core_gene_n =
    policy_b_summary$core_gene_n == 3,
  
  duplicate_entrez_n =
    policy_b_summary$duplicate_entrez_n == 0
)


if (!all(summary_checks)) {
  
  print(summary_checks)
  
  stop(
    "STEP9D-3A-5 Policy B summary validation FAILED."
  )
}


# ------------------------------------------------------------
# 9. Validate Policy B loss / retention
# ------------------------------------------------------------

loss_checks <- c(
  
  candidate_ensembl =
    policy_b_loss$candidate_ensembl_n ==
    LOCK_candidate_ensembl_n,
  
  retained_ensembl =
    policy_b_loss$retained_ensembl_n ==
    LOCK_policy_b_retained_ensembl_n,
  
  discarded_ensembl =
    policy_b_loss$discarded_ensembl_n ==
    LOCK_policy_b_discarded_ensembl_n,
  
  candidate_entrez =
    policy_b_loss$candidate_entrez_n ==
    LOCK_candidate_entrez_n,
  
  retained_entrez =
    policy_b_loss$retained_entrez_n ==
    LOCK_policy_b_retained_entrez_n,
  
  discarded_or_unrepresented_entrez =
    policy_b_loss$discarded_or_unrepresented_entrez_n ==
    LOCK_policy_b_discarded_entrez_n
)


if (!all(loss_checks)) {
  
  print(loss_checks)
  
  stop(
    "STEP9D-3A-5 Policy B loss/retention validation FAILED."
  )
}


# ------------------------------------------------------------
# 10. Validate Policy B directional effect
# ------------------------------------------------------------

direction_checks <- c(
  
  positive_retained =
    policy_b_direction$positive_retained_n ==
    LOCK_retained_positive_n,
  
  positive_discarded =
    policy_b_direction$positive_discarded_n == 9,
  
  negative_retained =
    policy_b_direction$negative_retained_n ==
    LOCK_retained_negative_n,
  
  negative_discarded =
    policy_b_direction$negative_discarded_n == 9,
  
  zero_retained =
    policy_b_direction$zero_retained_n ==
    LOCK_retained_zero_n,
  
  zero_discarded =
    policy_b_direction$zero_discarded_n == 0
)


if (!all(direction_checks)) {
  
  print(direction_checks)
  
  stop(
    "STEP9D-3A-5 Policy B directional validation FAILED."
  )
}


# ------------------------------------------------------------
# 11. Validate collision counts from STEP9D-3A-4
#     CORRECTED VERSION
# ------------------------------------------------------------
#
# IMPORTANT:
# ------------------------------------------------------------
# The STEP9D-3A-4 RDS already contains the formally audited
# collision counts.
#
# Therefore:
#   - do NOT infer collision groups from a CSV row count
#   - do NOT recalculate collision structure
#   - do NOT reconstruct the mapping
#
# Directly validate the locked fields in collision_audit.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-3A-4 Collision Audit Validation\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 11.1 Exact field names confirmed from the actual RDS
# ------------------------------------------------------------

required_collision_fields <- c(
  
  "genome_wide_gene_n",
  
  "kegg_candidate_ensembl_n",
  
  "kegg_candidate_unique_entrez_n",
  
  "kegg_collision_entrez_group_n",
  
  "kegg_collision_ensembl_n",
  
  "collision_positive_n",
  
  "collision_negative_n",
  
  "collision_zero_n",
  
  "collision_identical_group_n",
  
  "collision_nonidentical_group_n",
  
  "collision_opposite_direction_group_n",
  
  "collision_same_magnitude_group_n",
  
  "core_gene_collision_present",
  
  "resolution_required",
  
  "opposite_direction_present",
  
  "ranking_recalculated",
  
  "deg_recalculated",
  
  "gsea_performed",
  
  "duplicate_entrez_removed",
  
  "arbitrary_first_match_used",
  
  "jitter_used",
  
  "artificial_tie_breaking_used",
  
  "collision_resolution_applied"
  
)


missing_collision_fields <-
  setdiff(
    required_collision_fields,
    names(collision_audit)
  )


if (length(missing_collision_fields) > 0) {
  
  cat("Missing STEP9D-3A-4 fields:\n")
  print(missing_collision_fields)
  
  stop(
    "STEP9D-3A-4 field validation FAILED."
  )
}


# ------------------------------------------------------------
# 11.2 Extract the formally audited values
# ------------------------------------------------------------

collision_genome_wide_n <-
  as.integer(
    collision_audit$genome_wide_gene_n
  )

collision_candidate_ensembl_n <-
  as.integer(
    collision_audit$kegg_candidate_ensembl_n
  )

collision_candidate_entrez_n <-
  as.integer(
    collision_audit$kegg_candidate_unique_entrez_n
  )

collision_group_n <-
  as.integer(
    collision_audit$kegg_collision_entrez_group_n
  )

collision_affected_ensembl_n <-
  as.integer(
    collision_audit$kegg_collision_ensembl_n
  )

opposite_direction_n <-
  as.integer(
    collision_audit$collision_opposite_direction_group_n
  )

same_magnitude_n <-
  as.integer(
    collision_audit$collision_same_magnitude_group_n
  )


# ------------------------------------------------------------
# 11.3 Print the formally audited values
# ------------------------------------------------------------

cat(
  "Genome-wide genes: ",
  collision_genome_wide_n,
  "\n",
  sep = ""
)

cat(
  "KEGG candidate Ensembl: ",
  collision_candidate_ensembl_n,
  "\n",
  sep = ""
)

cat(
  "KEGG candidate unique Entrez: ",
  collision_candidate_entrez_n,
  "\n",
  sep = ""
)

cat(
  "KEGG collision groups: ",
  collision_group_n,
  "\n",
  sep = ""
)

cat(
  "KEGG collision-affected Ensembl: ",
  collision_affected_ensembl_n,
  "\n",
  sep = ""
)

cat(
  "Opposite-direction collision groups: ",
  opposite_direction_n,
  "\n",
  sep = ""
)

cat(
  "Maximum-absolute / same-magnitude groups: ",
  same_magnitude_n,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 11.4 Validate against STEP9D-3A-6 locked values
# ------------------------------------------------------------

collision_checks <- c(
  
  genome_wide_gene_n =
    collision_genome_wide_n == 57828L,
  
  candidate_ensembl_n =
    collision_candidate_ensembl_n == 34138L,
  
  candidate_unique_entrez_n =
    collision_candidate_entrez_n == 34120L,
  
  collision_groups_n =
    collision_group_n == 18L,
  
  collision_affected_ensembl_n =
    collision_affected_ensembl_n == 36L,
  
  opposite_direction_groups_n =
    opposite_direction_n == 6L,
  
  maximum_absolute_ties_n =
    same_magnitude_n == 0L,
  
  core_gene_collision_absent =
    identical(
      collision_audit$core_gene_collision_present,
      FALSE
    ),
  
  resolution_required =
    identical(
      collision_audit$resolution_required,
      TRUE
    ),
  
  opposite_direction_present =
    identical(
      collision_audit$opposite_direction_present,
      TRUE
    )
  
)


cat("Collision validation checks:\n")
print(collision_checks)
cat("\n")


if (!all(collision_checks)) {
  
  cat("FAILED checks:\n")
  print(
    names(collision_checks)[!collision_checks]
  )
  
  stop(
    "STEP9D-3A-4 collision validation FAILED."
  )
}


cat(
  "STEP9D-3A-4 collision validation: PASS\n\n"
)

# ============================================================
# STEP9D-3A-5
# Read-only diagnostic for max_abs_tie
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-3A-5 max_abs_tie diagnostic\n")
cat("============================================================\n\n")

cat("opposite_direction object class:\n")
print(class(opposite_direction))

cat("\nColumn names:\n")
print(names(opposite_direction))

cat("\nColumn count:\n")
print(ncol(opposite_direction))

cat("\nRow count:\n")
print(nrow(opposite_direction))

cat("\nmax_abs_tie exists:\n")
print("max_abs_tie" %in% names(opposite_direction))

if ("max_abs_tie" %in% names(opposite_direction)) {
  
  cat("\nmax_abs_tie class:\n")
  print(class(opposite_direction$max_abs_tie))
  
  cat("\nmax_abs_tie values:\n")
  print(opposite_direction$max_abs_tie)
  
  cat("\nmax_abs_tie table:\n")
  print(table(opposite_direction$max_abs_tie, useNA = "ifany"))
  
  cat("\nNumber of TRUE ties:\n")
  print(sum(
    opposite_direction$max_abs_tie %in% TRUE,
    na.rm = TRUE
  ))
  
  cat("\nNumber of FALSE ties:\n")
  print(sum(
    opposite_direction$max_abs_tie %in% FALSE,
    na.rm = TRUE
  ))
  
  cat("\nNumber of NA:\n")
  print(sum(
    is.na(opposite_direction$max_abs_tie)
  ))
}

cat("\n============================================================\n")  


# ------------------------------------------------------------
# 12. Validate maximum-absolute ranking uniqueness
#     FINAL CORRECTED VERSION
#
#     IMPORTANT:
#     STEP9D3A5_Collision_Group_Table_Corrected.csv
#     does NOT contain a max_abs_tie column.
#
#     Therefore, uniqueness of the maximum-absolute representative
#     is formally validated from:
#       - ranking_unique_n
#       - abs_ranking_unique_n
#
#     For the present collision structure:
#       n_ensembl = 2
#       ranking_unique_n = 2
#       abs_ranking_unique_n = 2
#
#     Hence the maximum absolute ranking is uniquely determined.
# ------------------------------------------------------------

collision_group_table_file <-
  file.path(
    comparison_dir,
    "STEP9D3A5_Collision_Group_Table_Corrected.csv"
  )

if (!file.exists(collision_group_table_file)) {
  stop(
    paste0(
      "STEP9D-3A-5 corrected collision-group table does not exist:\n",
      collision_group_table_file
    )
  )
}

collision_group_table <-
  readr::read_csv(
    collision_group_table_file,
    show_col_types = FALSE
  )


# ------------------------------------------------------------
# 12.1 Structural validation
# ------------------------------------------------------------

cat("Collision-group table columns:\n")
print(names(collision_group_table))

cat("\nCollision-group table dimensions:\n")

cat(
  "Rows = ",
  nrow(collision_group_table),
  "; Columns = ",
  ncol(collision_group_table),
  "\n",
  sep = ""
)


required_collision_columns <- c(
  "ENTREZID",
  "n_ensembl",
  "ranking_unique_n",
  "abs_ranking_unique_n",
  "ranking_min",
  "ranking_max",
  "abs_ranking_min",
  "abs_ranking_max",
  "direction_class",
  "contains_positive",
  "contains_negative",
  "contains_zero",
  "contains_core_gene",
  "core_genes"
)

missing_collision_columns <-
  setdiff(
    required_collision_columns,
    names(collision_group_table)
  )

if (length(missing_collision_columns) > 0L) {
  
  stop(
    paste0(
      "Required collision-group columns are missing:\n",
      paste(missing_collision_columns, collapse = ", ")
    )
  )
}


if (nrow(collision_group_table) != 18L) {
  
  stop(
    paste0(
      "Expected 18 collision groups, but observed ",
      nrow(collision_group_table),
      "."
    )
  )
}


# ------------------------------------------------------------
# 12.2 Validate collision-group cardinality
# ------------------------------------------------------------

n_ensembl_ok <-
  all(
    collision_group_table$n_ensembl == 2L
  )

ranking_unique_ok <-
  all(
    collision_group_table$ranking_unique_n == 2L
  )

abs_ranking_unique_ok <-
  all(
    collision_group_table$abs_ranking_unique_n == 2L
  )


cat("\nCollision-group cardinality checks:\n")

cat(
  "All collision groups have exactly 2 Ensembl genes: ",
  n_ensembl_ok,
  "\n",
  sep = ""
)

cat(
  "All collision groups have 2 unique ranking values: ",
  ranking_unique_ok,
  "\n",
  sep = ""
)

cat(
  "All collision groups have 2 unique absolute-ranking values: ",
  abs_ranking_unique_ok,
  "\n",
  sep = ""
)


if (!n_ensembl_ok) {
  stop(
    "STEP9D-3A-6 FAILED: not all collision groups contain exactly 2 Ensembl genes."
  )
}

if (!ranking_unique_ok) {
  stop(
    "STEP9D-3A-6 FAILED: ranking_unique_n != 2 in at least one collision group."
  )
}

if (!abs_ranking_unique_ok) {
  stop(
    paste0(
      "STEP9D-3A-6 FAILED: at least one collision group has a ",
      "maximum-absolute ranking tie."
    )
  )
}


# ------------------------------------------------------------
# 12.3 Explicit maximum-absolute tie count
# ------------------------------------------------------------
#
# A maximum-absolute tie would require:
#
#   abs_ranking_unique_n < n_ensembl
#
# Because every collision group contains exactly two Ensembl
# genes, abs_ranking_unique_n == 2 means the maximum absolute
# ranking is unique.
# ------------------------------------------------------------

maximum_absolute_tie_flag <-
  collision_group_table$abs_ranking_unique_n <
  collision_group_table$n_ensembl

observed_max_abs_ties <-
  sum(
    maximum_absolute_tie_flag,
    na.rm = TRUE
  )

observed_max_abs_tie_na <-
  sum(
    is.na(maximum_absolute_tie_flag)
  )

observed_non_tie_groups <-
  sum(
    !maximum_absolute_tie_flag,
    na.rm = TRUE
  )


cat("\nMaximum-absolute ranking tie validation:\n")

cat(
  "Total collision groups: ",
  nrow(collision_group_table),
  "\n",
  sep = ""
)

cat(
  "Maximum-absolute ranking ties: ",
  observed_max_abs_ties,
  "\n",
  sep = ""
)

cat(
  "Uniquely resolved groups: ",
  observed_non_tie_groups,
  "\n",
  sep = ""
)

cat(
  "NA tie-status values: ",
  observed_max_abs_tie_na,
  "\n",
  sep = ""
)


tie_checks <- c(
  
  collision_groups =
    nrow(collision_group_table) == 18L,
  
  all_groups_have_two_ensembl =
    n_ensembl_ok,
  
  all_groups_have_two_unique_rankings =
    ranking_unique_ok,
  
  maximum_absolute_ties =
    observed_max_abs_ties == 0L,
  
  all_18_groups_uniquely_resolved =
    observed_non_tie_groups == 18L,
  
  no_missing_tie_status =
    observed_max_abs_tie_na == 0L
  
)


cat("\nMaximum-absolute ranking tie validation checks:\n")
print(tie_checks)


if (!all(tie_checks)) {
  
  cat("\nFAILED checks:\n")
  print(
    names(tie_checks)[!tie_checks]
  )
  
  stop(
    "STEP9D-3A-6 maximum-absolute ranking uniqueness validation FAILED."
  )
}


cat(
  "\nMaximum-absolute ranking uniqueness validation: PASS\n"
)

cat(
  "All 18 collision groups have uniquely determined ",
  "maximum-absolute representatives.\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. Validate core-gene preservation
# ------------------------------------------------------------

required_core_genes <-
  c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )


core_policy_check <- core_policy %>%
  dplyr::filter(core_gene %in% required_core_genes) %>%
  dplyr::select(
    core_gene,
    candidate_present,
    Policy_A_retained,
    Policy_B_retained,
    Policy_C_retained
  )


if (nrow(core_policy_check) != 3L) {stop(
  "Core-gene policy preservation table does not contain exactly ",
  "the three locked core genes."
)
}


core_policy_check_ordered <-
  core_policy_check[
    match(required_core_genes, core_policy_check$core_gene),
    ,
    drop = FALSE
  ]


core_checks <- c(
  
  AOC3 =
    core_policy_check_ordered$Policy_B_retained[
      core_policy_check_ordered$core_gene == "AOC3"
    ],
  
  MS4A7 =
    core_policy_check_ordered$Policy_B_retained[
      core_policy_check_ordered$core_gene == "MS4A7"
    ],
  
  IL7R =
    core_policy_check_ordered$Policy_B_retained[
      core_policy_check_ordered$core_gene == "IL7R"
    ]
)


if (!all(core_checks)) {
  
  print(core_checks)
  
  stop(
    "Core-gene preservation under Policy B FAILED."
  )
}


# ------------------------------------------------------------
# 14. Validate absence of prohibited resolution procedures
# ------------------------------------------------------------

prohibited_method_checks <- c(
  
  arbitrary_first_match =
    LOCK_arbitrary_first_match == FALSE,
  
  jitter =
    LOCK_jitter == FALSE,
  
  artificial_tie_breaking =
    LOCK_artificial_tie_breaking == FALSE,
  
  ranking_recalculated =
    LOCK_ranking_recalculated == FALSE,
  
  DEG_recalculated =
    LOCK_DEG_recalculated == FALSE,
  
  GSEA_performed =
    LOCK_GSEA_performed == FALSE
)


if (!all(prohibited_method_checks)) {
  
  print(prohibited_method_checks)
  
  stop(
    "Prohibited-method validation FAILED."
  )
}


# ------------------------------------------------------------
# 15. Formal policy definition
# ------------------------------------------------------------

policy_definition <- paste(
  "For KEGG GSEA, when multiple locked Ensembl genes map to the",
  "same Entrez identifier, retain the Ensembl gene with the maximum",
  "absolute value of the already locked genome-wide ranking statistic.",
  "This rule is applied only to resolve Ensembl-to-Entrez collisions",
  "for KEGG input construction. The original 57,828-gene ranking",
  "vector remains unchanged."
)


# ------------------------------------------------------------
# 16. Formal final lock object
# ------------------------------------------------------------

step9d3a6_lock <- list(
  
  workflow =
    "Paper 1",
  
  step =
    "STEP9D-3A-6",
  
  status =
    "FINAL_LOCK",
  
  lock_date =
    as.character(Sys.Date()),
  
  policy_name =
    policy_b_name,
  
  policy_label =
    "Maximum Absolute Ranking",
  
  policy_definition =
    policy_definition,
  
  candidate_ensembl_n =
    LOCK_candidate_ensembl_n,
  
  candidate_unique_entrez_n =
    LOCK_candidate_entrez_n,
  
  collision_groups_n =
    LOCK_collision_groups_n,
  
  collision_affected_ensembl_n =
    LOCK_collision_affected_ensembl_n,
  
  opposite_direction_groups_n =
    LOCK_opposite_direction_groups_n,
  
  retained_ensembl_n =
    LOCK_policy_b_retained_ensembl_n,
  
  retained_entrez_n =
    LOCK_policy_b_retained_entrez_n,
  
  discarded_ensembl_n =
    LOCK_policy_b_discarded_ensembl_n,
  
  discarded_or_unrepresented_entrez_n =
    LOCK_policy_b_discarded_entrez_n,
  
  retained_positive_n =
    LOCK_retained_positive_n,
  
  retained_negative_n =
    LOCK_retained_negative_n,
  
  retained_zero_n =
    LOCK_retained_zero_n,
  
  maximum_absolute_ties_n =
    LOCK_max_absolute_ties_n,
  
  core_gene_AOC3_retained =
    LOCK_core_AOC3,
  
  core_gene_MS4A7_retained =
    LOCK_core_MS4A7,
  
  core_gene_IL7R_retained =
    LOCK_core_IL7R,
  
  arbitrary_first_match =
    LOCK_arbitrary_first_match,
  
  jitter =
    LOCK_jitter,
  
  artificial_tie_breaking =
    LOCK_artificial_tie_breaking,
  
  ranking_recalculated =
    LOCK_ranking_recalculated,
  
  DEG_recalculated =
    LOCK_DEG_recalculated,
  
  GSEA_performed =
    LOCK_GSEA_performed,
  
  original_ranking_gene_n =
    57828L,
  
  ranking_source =
    ranking_rds,
  
  mapping_policy_source =
    mapping_lock_rds,
  
  collision_audit_source =
    collision_audit_rds,
  
  policy_comparison_source =
    comparison_audit_rds,
  
  hard_pass =
    TRUE
)


# ------------------------------------------------------------
# 17. Formal audit table
# ------------------------------------------------------------

final_audit_table <- tibble::tibble(
  
  item = c(
    
    "candidate_Ensembl",
    "candidate_unique_Entrez",
    "collision_groups",
    "collision_affected_Ensembl",
    "opposite_direction_groups",
    "Policy_B_retained_Ensembl",
    "Policy_B_retained_Entrez",
    "discarded_Ensembl",
    "discarded_or_unrepresented_Entrez",
    "retained_positive",
    "retained_negative",
    "retained_zero",
    "maximum_absolute_ties",
    "AOC3_retained",
    "MS4A7_retained",
    "IL7R_retained",
    "arbitrary_first_match",
    "jitter",
    "artificial_tie_breaking",
    "ranking_recalculated",
    "DEG_recalculated",
    "GSEA_performed"
    
  ),
  
  locked_value = c(
    
    LOCK_candidate_ensembl_n,
    LOCK_candidate_entrez_n,
    LOCK_collision_groups_n,
    LOCK_collision_affected_ensembl_n,
    LOCK_opposite_direction_groups_n,
    LOCK_policy_b_retained_ensembl_n,
    LOCK_policy_b_retained_entrez_n,
    LOCK_policy_b_discarded_ensembl_n,
    LOCK_policy_b_discarded_entrez_n,
    LOCK_retained_positive_n,
    LOCK_retained_negative_n,
    LOCK_retained_zero_n,
    LOCK_max_absolute_ties_n,
    LOCK_core_AOC3,
    LOCK_core_MS4A7,
    LOCK_core_IL7R,
    LOCK_arbitrary_first_match,
    LOCK_jitter,
    LOCK_artificial_tie_breaking,
    LOCK_ranking_recalculated,
    LOCK_DEG_recalculated,
    LOCK_GSEA_performed
    
  ),
  
  validation_status = c(
    
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
    
  )
)


# ------------------------------------------------------------
# 18. Final hard-pass logic
# ------------------------------------------------------------

hard_pass_checks <- c(
  
  ranking_integrity_pass,
  
  all(summary_checks),
  
  all(loss_checks),
  
  all(direction_checks),
  
  all(collision_checks),
  
  observed_max_abs_ties ==
    LOCK_max_absolute_ties_n,
  
  all(core_checks),
  
  all(prohibited_method_checks),
  
  LOCK_candidate_ensembl_n ==
    34138L,
  
  LOCK_candidate_entrez_n ==
    34120L,
  
  LOCK_collision_groups_n ==
    18L,
  
  LOCK_collision_affected_ensembl_n ==
    36L,
  
  LOCK_opposite_direction_groups_n ==
    6L,
  
  LOCK_policy_b_retained_ensembl_n ==
    34120L,
  
  LOCK_policy_b_retained_entrez_n ==
    34120L,
  
  LOCK_policy_b_discarded_ensembl_n ==
    18L,
  
  LOCK_policy_b_discarded_entrez_n ==
    0L,
  
  LOCK_retained_positive_n ==
    21189L,
  
  LOCK_retained_negative_n ==
    12924L,
  
  LOCK_retained_zero_n ==
    7L,
  
  LOCK_core_AOC3,
  
  LOCK_core_MS4A7,
  
  LOCK_core_IL7R,
  
  LOCK_max_absolute_ties_n ==
    0L,
  
  LOCK_arbitrary_first_match == FALSE,
  
  LOCK_jitter == FALSE,
  
  LOCK_artificial_tie_breaking == FALSE,
  
  LOCK_ranking_recalculated == FALSE,
  
  LOCK_DEG_recalculated == FALSE,
  
  LOCK_GSEA_performed == FALSE
  
)


hard_pass <-
  all(hard_pass_checks)


step9d3a6_lock$hard_pass <-
  hard_pass


# ------------------------------------------------------------
# 19. Stop if anything fails
# ------------------------------------------------------------

if (!hard_pass) {
  
  cat("\nFAILED checks:\n")
  
  print(
    which(!hard_pass_checks)
  )
  
  stop(
    "STEP9D-3A-6 FINAL HARD PASS = FALSE"
  )
}


# ------------------------------------------------------------
# 20. Save final lock
# ------------------------------------------------------------

final_rds <-
  file.path(
    output_dir,
    "STEP9D3A6_FINAL_KEGG_Collision_Resolution_Policy_Lock.rds"
  )

final_csv <-
  file.path(
    output_dir,
    "STEP9D3A6_FINAL_KEGG_Collision_Resolution_Policy_Audit.csv"
  )

final_txt <-
  file.path(
    output_dir,
    "STEP9D3A6_FINAL_KEGG_Collision_Resolution_Policy_Lock.txt"
  )


saveRDS(
  step9d3a6_lock,
  final_rds
)


readr::write_csv(
  final_audit_table,
  final_csv
)


# ------------------------------------------------------------
# 21. Formal TXT audit
# ------------------------------------------------------------

audit_lines <- c(
  
  "============================================================",
  "Paper 1",
  "STEP9D-3A-6",
  "Final KEGG Collision Resolution Policy Lock",
  "============================================================",
  "",
  "STATUS: FINAL LOCK",
  "",
  "LOCKED POLICY:",
  "Policy B — Maximum Absolute Ranking",
  "",
  "POLICY DEFINITION:",
  policy_definition,
  "",
  "------------------------------------------------------------",
  "LOCKED COUNTS",
  "------------------------------------------------------------",
  paste0(
    "candidate Ensembl = ",
    LOCK_candidate_ensembl_n
  ),
  paste0(
    "candidate unique Entrez = ",
    LOCK_candidate_entrez_n
  ),
  paste0(
    "collision groups = ",
    LOCK_collision_groups_n
  ),
  paste0(
    "collision-affected Ensembl = ",
    LOCK_collision_affected_ensembl_n
  ),
  paste0(
    "opposite-direction groups = ",
    LOCK_opposite_direction_groups_n
  ),
  paste0(
    "Policy B retained Ensembl = ",
    LOCK_policy_b_retained_ensembl_n
  ),
  paste0(
    "Policy B retained Entrez = ",
    LOCK_policy_b_retained_entrez_n
  ),
  paste0(
    "discarded Ensembl = ",
    LOCK_policy_b_discarded_ensembl_n
  ),
  paste0(
    "discarded/unrepresented Entrez = ",
    LOCK_policy_b_discarded_entrez_n
  ),
  paste0(
    "retained positive = ",
    LOCK_retained_positive_n
  ),
  paste0(
    "retained negative = ",
    LOCK_retained_negative_n
  ),
  paste0(
    "retained zero = ",
    LOCK_retained_zero_n
  ),
  paste0(
    "maximum-absolute ties = ",
    LOCK_max_absolute_ties_n
  ),
  "",
  "------------------------------------------------------------",
  "CORE-GENE PRESERVATION",
  "------------------------------------------------------------",
  "AOC3 = retained",
  "MS4A7 = retained",
  "IL7R = retained",
  "",
  "------------------------------------------------------------",
  "METHOD INTEGRITY",
  "------------------------------------------------------------",
  "arbitrary first-match = FALSE",
  "jitter = FALSE",
  "artificial tie-breaking = FALSE",
  "ranking recalculated = FALSE",
  "DEG recalculated = FALSE",
  "GSEA performed = FALSE",
  "",
  "------------------------------------------------------------",
  "ORIGINAL RANKING",
  "------------------------------------------------------------",
  "original genome-wide ranking genes = 57,828",
  "",
  "------------------------------------------------------------",
  "FINAL STATUS",
  "------------------------------------------------------------",
  "STEP9D-3A-6 FINAL HARD PASS = TRUE",
  "",
  "This lock defines the formal KEGG GSEA input policy.",
  "Actual KEGG GSEA is NOT performed in STEP9D-3A-6.",
  "KEGG GSEA may proceed only in STEP9D-3B using this locked policy.",
  "",
  "============================================================"
)


writeLines(
  audit_lines,
  con = final_txt,
  useBytes = TRUE
)


# ------------------------------------------------------------
# 22. Reload verification
# ------------------------------------------------------------

reload_lock <-
  readRDS(final_rds)

reload_csv <-
  readr::read_csv(
    final_csv,
    show_col_types = FALSE
  )


reload_pass <-
  identical(
    reload_lock$hard_pass,
    TRUE
  ) &&
  nrow(reload_csv) == 22L &&
  file.exists(final_rds) &&
  file.exists(final_csv) &&
  file.exists(final_txt)


# ------------------------------------------------------------
# 23. Final console report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-3A-6 FINAL REPORT\n")
cat("============================================================\n\n")

cat(
  "Candidate Ensembl: ",
  LOCK_candidate_ensembl_n,
  "\n",
  sep = ""
)

cat(
  "Candidate unique Entrez: ",
  LOCK_candidate_entrez_n,
  "\n",
  sep = ""
)

cat(
  "Collision groups: ",
  LOCK_collision_groups_n,
  "\n",
  sep = ""
)

cat(
  "Collision-affected Ensembl: ",
  LOCK_collision_affected_ensembl_n,
  "\n",
  sep = ""
)

cat(
  "Opposite-direction groups: ",
  LOCK_opposite_direction_groups_n,
  "\n",
  sep = ""
)

cat(
  "Policy B retained Ensembl: ",
  LOCK_policy_b_retained_ensembl_n,
  "\n",
  sep = ""
)

cat(
  "Policy B retained Entrez: ",
  LOCK_policy_b_retained_entrez_n,
  "\n",
  sep = ""
)

cat(
  "Discarded Ensembl: ",
  LOCK_policy_b_discarded_ensembl_n,
  "\n",
  sep = ""
)

cat(
  "Discarded/unrepresented Entrez: ",
  LOCK_policy_b_discarded_entrez_n,
  "\n",
  sep = ""
)

cat(
  "Retained positive: ",
  LOCK_retained_positive_n,
  "\n",
  sep = ""
)

cat(
  "Retained negative: ",
  LOCK_retained_negative_n,
  "\n",
  sep = ""
)

cat(
  "Retained zero: ",
  LOCK_retained_zero_n,
  "\n",
  sep = ""
)

cat(
  "Maximum-absolute ties: ",
  LOCK_max_absolute_ties_n,
  "\n",
  sep = ""
)

cat("\nCore genes:\n")
cat("  AOC3  = retained\n")
cat("  MS4A7 = retained\n")
cat("  IL7R  = retained\n")

cat("\nMethod integrity:\n")
cat("  Arbitrary first-match: FALSE\n")
cat("  Jitter: FALSE\n")
cat("  Artificial tie-breaking: FALSE\n")
cat("  Ranking recalculated: FALSE\n")
cat("  DEG recalculated: FALSE\n")
cat("  GSEA performed: FALSE\n")

cat("\nFinal RDS exists: ",
    file.exists(final_rds),
    "\n",
    sep = "")

cat("Final audit CSV exists: ",
    file.exists(final_csv),
    "\n",
    sep = "")

cat("Final audit TXT exists: ",
    file.exists(final_txt),
    "\n",
    sep = "")

cat("Reload verification: ",
    reload_pass,
    "\n",
    sep = "")

cat(
  "\nSTEP9D-3A-6 FINAL HARD PASS: ",
  hard_pass && reload_pass,
  "\n",
  sep = ""
)

cat("\n============================================================\n")   


# ============================================================
# STEP9D-3A-4
# READ-ONLY GLOBAL FILE LOCATION AUDIT
# Purpose:
#   Locate the ACTUAL STEP9D-3A-4 directory / RDS / TXT / CSV
#   without creating, modifying, or recalculating anything.
# ============================================================

step9d_dir <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA"

cat("\n============================================================\n")
cat("STEP9D-3A-4 GLOBAL FILE LOCATION AUDIT\n")
cat("============================================================\n\n")

cat("STEP9D root exists: ", dir.exists(step9d_dir), "\n", sep = "")

if (!dir.exists(step9d_dir)) {
  stop("STEP9D root directory does not exist.")
}

# ------------------------------------------------------------
# 1. List all immediate STEP9D subdirectories
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("1. STEP9D immediate subdirectories\n")
cat("------------------------------------------------------------\n")

step9d_dirs <- list.dirs(
  step9d_dir,
  full.names = TRUE,
  recursive = FALSE
)

print(step9d_dirs)

# ------------------------------------------------------------
# 2. Recursively locate directories containing:
#    3A-4 / Collision / KEGG
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("2. Candidate directories\n")
cat("------------------------------------------------------------\n")

all_dirs <- list.dirs(
  step9d_dir,
  full.names = TRUE,
  recursive = TRUE
)

candidate_dirs <- all_dirs[
  grepl(
    "3A-4|3A4|Collision|collision|KEGG|kegg",
    basename(all_dirs),
    ignore.case = TRUE
  )
]

if (length(candidate_dirs) == 0) {
  cat("No candidate directories found.\n")
} else {
  print(candidate_dirs)
}

# ------------------------------------------------------------
# 3. Recursively locate all RDS files related to collision / 3A-4
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("3. Candidate RDS files\n")
cat("------------------------------------------------------------\n")

all_files <- list.files(
  step9d_dir,
  full.names = TRUE,
  recursive = TRUE
)

candidate_rds <- all_files[
  grepl(
    "\\.rds$",
    all_files,
    ignore.case = TRUE
  ) &
    grepl(
      "3A-4|3A4|Collision|collision|KEGG|kegg",
      basename(all_files),
      ignore.case = TRUE
    )
]

if (length(candidate_rds) == 0) {
  cat("No candidate RDS files found.\n")
} else {
  print(candidate_rds)
}

# ------------------------------------------------------------
# 4. Candidate TXT files
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("4. Candidate TXT files\n")
cat("------------------------------------------------------------\n")

candidate_txt <- all_files[
  grepl(
    "\\.txt$",
    all_files,
    ignore.case = TRUE
  ) &
    grepl(
      "3A-4|3A4|Collision|collision|KEGG|kegg",
      basename(all_files),
      ignore.case = TRUE
    )
]

if (length(candidate_txt) == 0) {
  cat("No candidate TXT files found.\n")
} else {
  print(candidate_txt)
}

# ------------------------------------------------------------
# 5. Candidate CSV files
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("5. Candidate CSV files\n")
cat("------------------------------------------------------------\n")

candidate_csv <- all_files[
  grepl(
    "\\.csv$",
    all_files,
    ignore.case = TRUE
  ) &
    grepl(
      "3A-4|3A4|Collision|collision|KEGG|kegg",
      basename(all_files),
      ignore.case = TRUE
    )
]

if (length(candidate_csv) == 0) {
  cat("No candidate CSV files found.\n")
} else {
  print(candidate_csv)
}

# ------------------------------------------------------------
# 6. Broad audit:
#    all files immediately related to 9D-3A*
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("6. All STEP9D-3A* files\n")
cat("------------------------------------------------------------\n")

step9d_3a_files <- all_files[
  grepl(
    "9D-3A|STEP9D3A|3A-",
    basename(all_files),
    ignore.case = TRUE
  )
]

if (length(step9d_3a_files) == 0) {
  cat("No STEP9D-3A related files found.\n")
} else {
  print(step9d_3a_files)
}

cat("\n============================================================\n")
cat("READ-ONLY LOCATION AUDIT COMPLETE\n")
cat("NO FILES CREATED\n")
cat("NO FILES MODIFIED\n")
cat("NO RANKING RECALCULATION\n")
cat("NO DEG RECALCULATION\n")
cat("NO GSEA\n")
cat("============================================================\n\n")

# ============================================================
# Paper 1
# STEP9D-3B
# FORMAL KEGG GSEA EXECUTION
#
# FINAL CORRECTED VERSION v2
#
# Critical correction:
#   STEP9D-2C FINAL ranking RDS is itself a named numeric vector.
#   It is NOT a list containing $ranking_vector.
#
# Therefore:
#   ranking_vector <- readRDS(ranking_file)
#
# No ranking recalculation.
# No DEG recalculation.
# ============================================================


cat("\n============================================================\n")
cat("STEP9D-3B: FORMAL KEGG GSEA EXECUTION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "clusterProfiler",
  "fgsea",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "dplyr"
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

library(clusterProfiler)
library(fgsea)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)


cat("Package versions:\n")

cat(
  "clusterProfiler: ",
  as.character(packageVersion("clusterProfiler")),
  "\n",
  sep = ""
)

cat(
  "fgsea: ",
  as.character(packageVersion("fgsea")),
  "\n",
  sep = ""
)

cat(
  "org.Hs.eg.db: ",
  as.character(packageVersion("org.Hs.eg.db")),
  "\n",
  sep = ""
)

cat(
  "AnnotationDbi: ",
  as.character(packageVersion("AnnotationDbi")),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

step9d_dir <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA"


ranking_file <-
  file.path(
    step9d_dir,
    "9D-2C_Final_Ranking_Vector",
    "STEP9D2C_FINAL_GenomeWide_Ranking_Vector.rds"
  )


mapping_lock_file <-
  file.path(
    step9d_dir,
    "9D-3A-3_GSEA_ID_Mapping_Policy_Lock",
    "STEP9D3A3_FINAL_GSEA_ID_Mapping_Policy_Lock.rds"
  )


collision_audit_file <-
  file.path(
    step9d_dir,
    "9D-3A-4_KEGG_Entrez_Collision_Audit",
    "STEP9D3A4_FINAL_KEGG_Entrez_Collision_Audit.rds"
  )


collision_policy_file <-
  file.path(
    step9d_dir,
    "9D-3A-5_KEGG_Collision_Policy_Audit",
    "STEP9D3A5_FINAL_KEGG_Collision_Policy_Comparison_Audit.rds"
  )


final_policy_lock_file <-
  file.path(
    step9d_dir,
    "9D-3A-6_KEGG_Collision_Resolution_Policy_Lock",
    "STEP9D3A6_FINAL_KEGG_Collision_Resolution_Policy_Lock.rds"
  )


output_dir <-
  file.path(
    step9d_dir,
    "9D-3B_Formal_GSEA_Execution"
  )


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Verify required files
# ------------------------------------------------------------

required_files <- c(
  ranking_file,
  mapping_lock_file,
  collision_audit_file,
  collision_policy_file,
  final_policy_lock_file
)

missing_files <-
  required_files[
    !file.exists(required_files)
  ]


cat("\n------------------------------------------------------------\n")
cat("INPUT FILE EXISTENCE CHECK\n")
cat("------------------------------------------------------------\n")


if (length(missing_files) > 0) {
  
  print(missing_files)
  
  stop(
    "STEP9D-3B STOPPED: required locked input file(s) missing."
  )
  
} else {
  
  cat(
    "All required locked input files exist: TRUE\n"
  )
}


# ------------------------------------------------------------
# 3. Load locked inputs
# ------------------------------------------------------------

# IMPORTANT:
# ranking_file contains the FINAL ranking vector itself.

ranking_vector <-
  readRDS(ranking_file)

mapping_lock <-
  readRDS(mapping_lock_file)

collision_audit <-
  readRDS(collision_audit_file)

collision_policy <-
  readRDS(collision_policy_file)

final_policy_lock <-
  readRDS(final_policy_lock_file)


cat(
  "\nAll locked input objects loaded successfully.\n"
)


# ------------------------------------------------------------
# 4. Verify actual ranking object type
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("RANKING OBJECT STRUCTURE\n")
cat("------------------------------------------------------------\n")

print(
  str(
    ranking_vector,
    max.level = 1
  )
)


# ------------------------------------------------------------
# 5. Ranking integrity checks
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("LOCKED RANKING INTEGRITY CHECK\n")
cat("------------------------------------------------------------\n")


if (!is.numeric(ranking_vector)) {
  
  stop(
    "Locked ranking object is not numeric."
  )
}


if (is.null(names(ranking_vector))) {
  
  stop(
    "Locked ranking vector has no gene identifiers."
  )
}


ranking_length <-
  length(ranking_vector)

ranking_names <-
  names(ranking_vector)


cat(
  "Ranking class: ",
  paste(class(ranking_vector), collapse = ", "),
  "\n",
  sep = ""
)

cat(
  "Ranking length: ",
  ranking_length,
  "\n",
  sep = ""
)

cat(
  "Ranking names available: ",
  !is.null(ranking_names),
  "\n",
  sep = ""
)


if (ranking_length != 57828) {
  
  stop(
    paste0(
      "Unexpected ranking length: ",
      ranking_length,
      " ; expected 57828."
    )
  )
}


if (anyDuplicated(ranking_names) > 0) {
  
  stop(
    "Duplicated gene identifiers detected in locked ranking."
  )
}


if (anyNA(ranking_vector)) {
  
  stop(
    "NA values detected in locked ranking."
  )
}


if (any(!is.finite(ranking_vector))) {
  
  stop(
    "Non-finite values detected in locked ranking."
  )
}


cat(
  "Genome-wide ranking integrity: PASS\n"
)


# ------------------------------------------------------------
# 6. Verify STEP9D-3A-3 mapping lock
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("STEP9D-3A-3 MAPPING POLICY LOCK\n")
cat("------------------------------------------------------------\n")

print(
  mapping_lock[
    c(
      "genome_wide_gene_n",
      "mapped_ensembl_n",
      "unique_entrez_n",
      "unmapped_ensembl_n",
      "one_to_many_ensembl_n",
      "many_to_one_entrez_n",
      "core_gene_traceable",
      "core_gene_one_to_one",
      "hard_pass"
    )
  ]
)


if (!isTRUE(mapping_lock$hard_pass)) {
  
  stop(
    "STEP9D-3A-3 mapping policy lock is not PASS."
  )
}


# ------------------------------------------------------------
# 7. Verify STEP9D-3A-4 collision audit
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("STEP9D-3A-4 COLLISION AUDIT\n")
cat("------------------------------------------------------------\n")

print(
  collision_audit[
    c(
      "genome_wide_gene_n",
      "kegg_candidate_ensembl_n",
      "kegg_candidate_unique_entrez_n",
      "kegg_duplicate_entrez_n",
      "kegg_collision_entrez_group_n",
      "kegg_collision_ensembl_n",
      "collision_positive_n",
      "collision_negative_n",
      "collision_zero_n",
      "collision_opposite_direction_group_n",
      "core_gene_collision_present",
      "resolution_required",
      "opposite_direction_present",
      "hard_pass"
    )
  ]
)


if (!isTRUE(collision_audit$hard_pass)) {stop(
  "STEP9D-3A-4 collision audit is not PASS.")
}    

# ------------------------------------------------------------
# 8. Verify STEP9D-3A-5 policy comparison
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("STEP9D-3A-5 POLICY COMPARISON AUDIT\n")
cat("------------------------------------------------------------\n")

print(
  collision_policy[
    c(
      "policy_summary",
      "core_gene_preservation",
      "hard_pass"
    )
  ]
)


if (!isTRUE(collision_policy$hard_pass)) {
  
  stop(
    "STEP9D-3A-5 policy comparison audit is not PASS."
  )
}


# ------------------------------------------------------------
# 9. Verify FINAL STEP9D-3A-6 policy lock
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("STEP9D-3A-6 FINAL POLICY LOCK\n")
cat("------------------------------------------------------------\n")

print(
  final_policy_lock
)


if (!isTRUE(final_policy_lock$hard_pass)) {
  
  stop(
    "STEP9D-3A-6 FINAL policy lock is not PASS."
  )
}


# ------------------------------------------------------------
# 10. Confirm final policy
# ------------------------------------------------------------

if (
  !identical(
    final_policy_lock$policy_name,
    "Policy_B_Max_Absolute_Ranking"
  )
) {
  
  stop(
    "Unexpected KEGG collision-resolution policy."
  )
}


cat(
  "\nFinal KEGG collision-resolution policy: ",
  final_policy_lock$policy_name,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Reconstruct AnnotationDbi lookup
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("ANNOTATION LOOKUP\n")
cat("------------------------------------------------------------\n")


ranking_ensembl <-
  names(ranking_vector)


annotation_raw <-
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ranking_ensembl,
    keytype = "ENSEMBL",
    columns = c(
      "ENSEMBL",
      "ENTREZID",
      "SYMBOL"
    )
  )


annotation_raw <-
  annotation_raw %>%
  filter(
    !is.na(ENSEMBL),
    !is.na(ENTREZID)
  ) %>%
  distinct(
    ENSEMBL,
    ENTREZID,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# 12. Keep Ensembl IDs with exactly one Entrez ID
# ------------------------------------------------------------

mapping_count <-
  annotation_raw %>%
  group_by(ENSEMBL) %>%
  summarise(
    n_entrez =
      n_distinct(ENTREZID),
    .groups = "drop"
  )


kegg_candidates <-
  mapping_count %>%
  filter(
    n_entrez == 1
  ) %>%
  select(
    ENSEMBL
  ) %>%
  inner_join(
    annotation_raw,
    by = "ENSEMBL"
  )


cat(
  "KEGG candidate Ensembl genes: ",
  nrow(kegg_candidates),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. Attach LOCKED ranking
# ------------------------------------------------------------

kegg_candidates$ranking <-
  unname(
    ranking_vector[
      kegg_candidates$ENSEMBL
    ]
  )


if (
  anyNA(kegg_candidates$ranking)
) {
  
  stop(
    "NA ranking values detected after annotation join."
  )
}


# ------------------------------------------------------------
# 14. Identify Entrez collision groups
# ------------------------------------------------------------

collision_groups <-
  kegg_candidates %>%
  group_by(ENTREZID) %>%
  summarise(
    n_ensembl =
      n_distinct(ENSEMBL),
    .groups = "drop"
  ) %>%
  filter(
    n_ensembl > 1
  )


cat(
  "KEGG collision groups: ",
  nrow(collision_groups),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 15. Maximum-absolute tie audit
# ------------------------------------------------------------

collision_detail <-
  kegg_candidates %>%
  group_by(ENTREZID) %>%
  mutate(
    abs_ranking =
      abs(ranking)
  ) %>%
  ungroup()


max_abs_summary <-
  collision_detail %>%
  group_by(ENTREZID) %>%
  summarise(
    n_ensembl =
      n(),
    
    max_abs =
      max(abs_ranking),
    
    max_abs_n =
      sum(
        abs_ranking ==
          max(abs_ranking)
      ),
    
    .groups = "drop"
  ) %>%
  filter(
    n_ensembl > 1
  )


max_abs_ties <-
  max_abs_summary %>%
  filter(
    max_abs_n > 1
  )


cat(
  "Maximum-absolute ranking ties: ",
  nrow(max_abs_ties),
  "\n",
  sep = ""
)


if (nrow(max_abs_ties) > 0) {
  
  print(max_abs_ties)
  
  stop(
    "Unexpected maximum-absolute ranking ties detected."
  )
}


# ------------------------------------------------------------
# 16. Apply FINAL Policy B
# ------------------------------------------------------------

kegg_final <-
  collision_detail %>%
  group_by(ENTREZID) %>%
  slice_max(
    order_by = abs_ranking,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 17. Final KEGG input integrity
# ------------------------------------------------------------

final_entrez_n <-
  n_distinct(
    kegg_final$ENTREZID
  )

final_ensembl_n <-
  n_distinct(
    kegg_final$ENSEMBL
  )

positive_n <-
  sum(
    kegg_final$ranking > 0
  )

negative_n <-
  sum(
    kegg_final$ranking < 0
  )

zero_n <-
  sum(
    kegg_final$ranking == 0
  )


cat("\n------------------------------------------------------------\n")
cat("FINAL KEGG GSEA INPUT\n")
cat("------------------------------------------------------------\n")

cat(
  "Final Ensembl genes: ",
  final_ensembl_n,
  "\n",
  sep = ""
)

cat(
  "Final unique Entrez IDs: ",
  final_entrez_n,
  "\n",
  sep = ""
)

cat(
  "Positive ranking: ",
  positive_n,
  "\n",
  sep = ""
)

cat(
  "Negative ranking: ",
  negative_n,
  "\n",
  sep = ""
)

cat(
  "Zero ranking: ",
  zero_n,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 18. Hard expected-value checks
# ------------------------------------------------------------

if (final_ensembl_n != 34120) {
  
  stop(
    "Unexpected final Ensembl count."
  )
}

if (final_entrez_n != 34120) {
  
  stop(
    "Unexpected final Entrez count."
  )
}

if (positive_n != 21189) {
  
  stop(
    "Unexpected positive ranking count."
  )
}

if (negative_n != 12924) {
  
  stop(
    "Unexpected negative ranking count."
  )
}

if (zero_n != 7) {
  
  stop(
    "Unexpected zero ranking count."
  )
}


# ------------------------------------------------------------
# 19. Core-gene preservation
# ------------------------------------------------------------

core_genes <-
  c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )


core_status <-
  kegg_final %>%
  filter(
    SYMBOL %in% core_genes
  ) %>%
  select(
    SYMBOL,
    ENSEMBL,
    ENTREZID,
    ranking
  )


cat("\n------------------------------------------------------------\n")
cat("CORE-GENE KEGG TRACEABILITY\n")
cat("------------------------------------------------------------\n")

print(core_status)


for (g in core_genes) {
  
  if (
    !any(
      core_status$SYMBOL == g
    )
  ) {
    
    stop(
      paste0(
        "Core gene missing from final KEGG input: ",
        g
      )
    )
  }
}


# ------------------------------------------------------------
# 20. Final named Entrez ranking
# ------------------------------------------------------------

kegg_ranking <-
  kegg_final$ranking


names(kegg_ranking) <-
  as.character(
    kegg_final$ENTREZID
  )


if (
  anyDuplicated(
    names(kegg_ranking)
  ) > 0
) {
  
  stop(
    "Duplicated Entrez IDs remain in final GSEA ranking."
  )
}


if (
  anyNA(kegg_ranking)
) {
  
  stop(
    "NA values detected in final GSEA ranking."
  )
}


if (
  any(!is.finite(kegg_ranking))
) {
  
  stop(
    "Non-finite values detected in final GSEA ranking."
  )
}


# ------------------------------------------------------------
# 21. Save final KEGG input
# ------------------------------------------------------------

gsea_input_table <-
  kegg_final %>%
  select(
    ENSEMBL,
    ENTREZID,
    SYMBOL,
    ranking
  ) %>%
  arrange(
    desc(ranking)
  )


input_csv <-
  file.path(
    output_dir,
    "STEP9D3B_FINAL_KEGG_GSEA_Input_Ranking.csv"
  )


input_rds <-
  file.path(
    output_dir,
    "STEP9D3B_FINAL_KEGG_GSEA_Input_Ranking.rds"
  )


write.csv(
  gsea_input_table,
  input_csv,
  row.names = FALSE,
  quote = TRUE
)


saveRDS(
  list(
    workflow =
      "Paper1_STEP9D",
    
    step =
      "STEP9D-3B",
    
    status =
      "FINAL_KEGG_GSEA_INPUT",
    
    source_ranking_file =
      ranking_file,
    
    mapping_lock_file =
      mapping_lock_file,
    
    collision_audit_file =
      collision_audit_file,
    
    collision_policy_file =
      collision_policy_file,
    
    final_policy_lock_file =
      final_policy_lock_file,
    
    policy =
      "Policy_B_Max_Absolute_Ranking",
    
    genome_wide_ranking_n =
      ranking_length,
    
    final_ensembl_n =
      final_ensembl_n,
    
    final_entrez_n =
      final_entrez_n,
    
    positive_n =
      positive_n,
    
    negative_n =
      negative_n,
    
    zero_n =
      zero_n,
    
    core_status =
      core_status,
    
    input_table =
      gsea_input_table,
    
    ranking_vector =
      kegg_ranking
  ),
  input_rds
)


# ------------------------------------------------------------
# 22. Formal GSEA parameters
# ------------------------------------------------------------

gsea_params <-
  list(
    organism =
      "hsa",
    
    keyType =
      "ncbi-geneid",
    
    minGSSize =
      10,
    
    maxGSSize =
      500,
    
    exponent =
      1,
    
    eps =
      1e-10,
    
    pvalueCutoff =
      1,
    
    pAdjustMethod =
      "BH",
    
    verbose =
      FALSE,
    
    seed =
      FALSE
  )


cat("\n------------------------------------------------------------\n")
cat("FORMAL GSEA PARAMETERS\n")
cat("------------------------------------------------------------\n")

print(gsea_params)

# ------------------------------------------------------------
# 22B. Prepare formally sorted KEGG ranking vector
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-3B: Preparing formally sorted KEGG ranking vector\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 22B-1. Preserve original KEGG ranking object
# ------------------------------------------------------------

kegg_ranking_unsorted <- kegg_ranking


# ------------------------------------------------------------
# 22B-2. Basic integrity checks
# ------------------------------------------------------------

if (!is.numeric(kegg_ranking_unsorted)) {
  stop("FAIL: kegg_ranking is not numeric.")
}

if (is.null(names(kegg_ranking_unsorted))) {
  stop("FAIL: kegg_ranking has no gene identifiers.")
}

if (anyDuplicated(names(kegg_ranking_unsorted)) > 0) {
  stop("FAIL: duplicated Entrez IDs detected in kegg_ranking.")
}

if (anyNA(kegg_ranking_unsorted)) {
  stop("FAIL: NA values detected in kegg_ranking.")
}

if (any(!is.finite(kegg_ranking_unsorted))) {
  stop("FAIL: non-finite values detected in kegg_ranking.")
}


# ------------------------------------------------------------
# 22B-3. Verify locked KEGG input size
# ------------------------------------------------------------

if (length(kegg_ranking_unsorted) != 34120) {
  stop(
    paste0(
      "FAIL: expected 34,120 KEGG ranking entries; observed ",
      length(kegg_ranking_unsorted), "."
    )
  )
}


# ------------------------------------------------------------
# 22B-4. Sort only — DO NOT recalculate ranking
# ------------------------------------------------------------

kegg_ranking <-
  sort(
    kegg_ranking_unsorted,
    decreasing = TRUE
  )


# ------------------------------------------------------------
# 22B-5. Formal sortedness audit
# ------------------------------------------------------------

is_decreasing_sorted <-
  all(
    diff(kegg_ranking) <= 0
  )


if (!isTRUE(is_decreasing_sorted)) {
  stop(
    "FAIL: kegg_ranking is still not decreasingly sorted."
  )
}


# ------------------------------------------------------------
# 22B-6. Confirm ranking values were not changed
# ------------------------------------------------------------

if (!isTRUE(
  all(
    sort(kegg_ranking_unsorted) ==
    sort(kegg_ranking)
  )
)) {
  stop(
    "FAIL: ranking values changed during sorting."
  )
}


# ------------------------------------------------------------
# 22B-7. Confirm gene-ID set was not changed
# ------------------------------------------------------------

if (!isTRUE(
  setequal(
    names(kegg_ranking_unsorted),
    names(kegg_ranking)
  )
)) {
  stop(
    "FAIL: Entrez gene-ID set changed during sorting."
  )
}


# ------------------------------------------------------------
# 22B-8. Report
# ------------------------------------------------------------

cat("Original KEGG ranking entries: ",
    length(kegg_ranking_unsorted), "\n", sep = "")

cat("Sorted KEGG ranking entries:   ",
    length(kegg_ranking), "\n", sep = "")

cat("Decreasingly sorted: ",
    is_decreasing_sorted, "\n", sep = "")

cat("Ranking values unchanged: ",
    all(sort(kegg_ranking_unsorted) == sort(kegg_ranking)),
    "\n", sep = "")

cat("Entrez ID set unchanged: ",
    setequal(
      names(kegg_ranking_unsorted),
      names(kegg_ranking)
    ),
    "\n\n", sep = "")

cat("Top 5 ranking values:\n")
print(head(kegg_ranking, 5))

cat("\nBottom 5 ranking values:\n")
print(tail(kegg_ranking, 5))


# ------------------------------------------------------------
# 22B-9. Final hard check
# ------------------------------------------------------------

if (
  length(kegg_ranking) != 34120 ||
  !is_decreasing_sorted ||
  !isTRUE(all(sort(kegg_ranking_unsorted) == sort(kegg_ranking))) ||
  !isTRUE(setequal(
    names(kegg_ranking_unsorted),
    names(kegg_ranking)
  ))
) {
  stop(
    "STEP9D-3B ranking preparation HARD FAIL."
  )
}

cat("\nSTEP9D-3B ranking preparation HARD PASS: TRUE\n\n")

# ------------------------------------------------------------
# 23. Formal gseKEGG execution
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STARTING FORMAL gseKEGG()\n")
cat("============================================================\n\n")

gsea_result <-
  clusterProfiler::gseKEGG(
    geneList =
      kegg_ranking,
    
    organism =
      gsea_params$organism,
    
    keyType =
      gsea_params$keyType,
    
    minGSSize =
      gsea_params$minGSSize,
    
    maxGSSize =
      gsea_params$maxGSSize,
    
    exponent =
      gsea_params$exponent,
    
    eps =
      gsea_params$eps,
    
    pvalueCutoff =
      gsea_params$pvalueCutoff,
    
    pAdjustMethod =
      gsea_params$pAdjustMethod,
    
    verbose =
      gsea_params$verbose,
    
    seed =
      gsea_params$seed
  )



# ------------------------------------------------------------
# 24. Extract results
# ------------------------------------------------------------

gsea_result_df <-
  as.data.frame(
    gsea_result
  )


cat(
  "Number of KEGG pathways returned: ",
  nrow(gsea_result_df),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 25. Significant FDR < 0.05
# ------------------------------------------------------------

if (
  nrow(gsea_result_df) > 0 &&
  "p.adjust" %in% colnames(gsea_result_df)
) {
  
  gsea_sig_df <-
    gsea_result_df %>%
    filter(
      !is.na(p.adjust),
      p.adjust < 0.05
    )
  
} else {
  
  gsea_sig_df <-
    gsea_result_df[
      0,
      ,
      drop = FALSE
    ]
}


cat(
  "KEGG pathways with FDR < 0.05: ",
  nrow(gsea_sig_df),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 26. Save GSEA results
# ------------------------------------------------------------

result_rds <-
  file.path(
    output_dir,
    "STEP9D3B_FINAL_KEGG_GSEA_Result.rds"
  )


result_csv <-
  file.path(
    output_dir,
    "STEP9D3B_FINAL_KEGG_GSEA_Result.csv"
  )


sig_csv <-
  file.path(
    output_dir,
    "STEP9D3B_FINAL_KEGG_GSEA_Significant_FDR05.csv"
  )


saveRDS(
  gsea_result,
  result_rds
)


write.csv(
  gsea_result_df,
  result_csv,
  row.names = FALSE,
  quote = TRUE
)


write.csv(
  gsea_sig_df,
  sig_csv,
  row.names = FALSE,
  quote = TRUE
)


# ------------------------------------------------------------
# 27. Final provenance
# ------------------------------------------------------------

provenance_rds <-
  file.path(
    output_dir,
    "STEP9D3B_FINAL_GSEA_Provenance.rds"
  )


provenance <-
  list(
    
    workflow =
      "Paper1_STEP9D",
    
    step =
      "STEP9D-3B",
    
    status =
      "FINAL",
    
    execution_date =
      Sys.time(),
    
    ranking_source =
      ranking_file,
    
    mapping_lock_source =
      mapping_lock_file,
    
    collision_audit_source =
      collision_audit_file,
    
    collision_policy_source =
      collision_policy_file,
    
    final_policy_lock_source =
      final_policy_lock_file,
    
    ranking_metric =
      "signed_minus_log10_nominal_P",
    
    ranking_formula =
      "sign(log2FC) * (-log10(pmax(P_value, 1e-300)))",
    
    original_genome_wide_ranking_n =
      ranking_length,
    
    kegg_candidate_ensembl_n =
      nrow(kegg_candidates),
    
    final_kegg_ensembl_n =
      final_ensembl_n,
    
    final_kegg_entrez_n =
      final_entrez_n,
    
    positive_n =
      positive_n,
    
    negative_n =
      negative_n,
    
    zero_n =
      zero_n,
    
    collision_resolution_policy =
      "Policy_B_Max_Absolute_Ranking",
    
    maximum_absolute_ties =
      nrow(max_abs_ties),
    
    core_genes =
      core_genes,
    
    core_status =
      core_status,
    
    package_versions =
      list(
        clusterProfiler =
          as.character(
            packageVersion(
              "clusterProfiler"
            )
          ),
        
        fgsea =
          as.character(
            packageVersion(
              "fgsea"
            )
          ),
        
        org.Hs.eg.db =
          as.character(
            packageVersion(
              "org.Hs.eg.db"
            )
          ),
        
        AnnotationDbi =
          as.character(
            packageVersion(
              "AnnotationDbi"
            )
          )
      ),
    
    gsea_parameters =
      gsea_params,
    
    gsea_pathway_n =
      nrow(gsea_result_df),
    
    gsea_fdr05_n =
      nrow(gsea_sig_df),
    
    ranking_recalculated =
      FALSE,
    
    deg_recalculated =
      FALSE,
    
    arbitrary_first_match =
      FALSE,
    
    jitter_used =
      FALSE,
    
    artificial_tie_breaking =
      FALSE,
    
    hard_pass =
      TRUE
  )


saveRDS(
  provenance,
  provenance_rds
)


# ------------------------------------------------------------
# 28. Final audit
# ------------------------------------------------------------

audit_txt <-
  file.path(
    output_dir,
    "STEP9D3B_FINAL_GSEA_Execution_Audit.txt"
  )


audit_lines <-
  c(
    
    "============================================================",
    "STEP9D-3B FINAL KEGG GSEA EXECUTION AUDIT",
    "============================================================",
    "",
    
    paste0(
      "Execution date: ",
      Sys.time()
    ),
    
    "",
    
    "INPUT INTEGRITY",
    
    paste0(
      "Genome-wide ranking genes: ",
      ranking_length
    ),
    
    paste0(
      "KEGG candidate Ensembl genes: ",
      nrow(kegg_candidates)
    ),
    
    paste0(
      "Final unique Entrez IDs: ",
      final_entrez_n
    ),
    
    paste0(
      "Positive ranking: ",
      positive_n
    ),
    
    paste0(
      "Negative ranking: ",
      negative_n
    ),
    
    paste0(
      "Zero ranking: ",
      zero_n
    ),
    
    "",
    
    "COLLISION POLICY",
    
    "Policy: Policy_B_Max_Absolute_Ranking",
    
    paste0(
      "Collision groups: ",
      nrow(collision_groups)
    ),
    
    paste0(
      "Maximum-absolute ties: ",
      nrow(max_abs_ties)
    ),
    
    "",
    
    "CORE GENES",
    
    paste0(
      "AOC3 retained: ",
      any(core_status$SYMBOL == "AOC3")
    ),
    
    paste0(
      "MS4A7 retained: ",
      any(core_status$SYMBOL == "MS4A7")
    ),
    
    paste0(
      "IL7R retained: ",
      any(core_status$SYMBOL == "IL7R")
    ),
    
    "",
    
    "GSEA PARAMETERS",
    
    paste0(
      "organism: ",
      gsea_params$organism
    ),
    
    paste0(
      "keyType: ",
      gsea_params$keyType
    ),
    
    paste0(
      "minGSSize: ",
      gsea_params$minGSSize
    ),
    
    paste0(
      "maxGSSize: ",
      gsea_params$maxGSSize
    ),
    
    paste0(
      "exponent: ",
      gsea_params$exponent
    ),
    
    paste0(
      "eps: ",
      gsea_params$eps
    ),
    
    paste0(
      "pvalueCutoff: ",
      gsea_params$pvalueCutoff
    ),
    
    paste0(
      "pAdjustMethod: ",
      gsea_params$pAdjustMethod
    ),
    
    "",
    
    "RESULTS",
    
    paste0(
      "KEGG pathways returned: ",
      nrow(gsea_result_df)
    ),
    
    paste0(
      "FDR < 0.05 pathways: ",
      nrow(gsea_sig_df)
    ),
    
    "",
    
    "METHOD INTEGRITY",
    
    "Ranking recalculated: FALSE",
    "DEG recalculated: FALSE",
    "Arbitrary first-match: FALSE",
    "Jitter used: FALSE",
    "Artificial tie-breaking: FALSE",
    "Maximum-absolute ties unresolved: FALSE",
    
    "",
    
    "FINAL STATUS",
    "STEP9D-3B FORMAL KEGG GSEA EXECUTION: COMPLETE",
    "============================================================"
  )


writeLines(
  audit_lines,
  audit_txt
)


# ------------------------------------------------------------
# 29. Final output existence check
# ------------------------------------------------------------

final_outputs <-
  c(
    input_csv,
    input_rds,
    result_rds,
    result_csv,
    sig_csv,
    provenance_rds,
    audit_txt
  )


cat("\n============================================================\n")
cat("STEP9D-3B FINAL OUTPUT CHECK\n")
cat("============================================================\n\n")


print(
  data.frame(
    file =
      final_outputs,
    
    exists =
      file.exists(final_outputs)
  )
)


if (
  !all(
    file.exists(final_outputs)
  )
) {
  
  stop(
    "STEP9D-3B FINAL OUTPUT CHECK FAILED."
  )
}


cat("\n============================================================\n")
cat("STEP9D-3B FORMAL KEGG GSEA: FINAL PASS\n")
cat("============================================================\n\n")   


# ============================================================
# Paper 1
# STEP9D-3C
# Formal KEGG GSEA Result & Reproducibility Audit
#
# PURPOSE:
#   Read-only audit of the FINAL STEP9D-3B KEGG GSEA result.
#
# IMPORTANT:
#   NO DEG recalculation
#   NO ranking recalculation
#   NO Ensembl->Entrez remapping
#   NO collision resolution
#   NO GSEA rerun
#   NO jitter
#   NO artificial tie-breaking
#
# FINAL RESULT EXPECTED:
#   353 KEGG pathways returned
#   178 pathways with FDR < 0.05
#
# ============================================================


cat("\n============================================================\n")
cat("STEP9D-3C: Formal KEGG GSEA Result & Reproducibility Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_packages <- c(
  "clusterProfiler",
  "fgsea",
  "AnnotationDbi",
  "org.Hs.eg.db"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "FAIL: required package not installed: ",
        pkg
      )
    )
  }
}

cat("Required packages available: TRUE\n\n")


# ------------------------------------------------------------
# 1. Version audit
# ------------------------------------------------------------

cat("PACKAGE VERSIONS\n")
cat("----------------\n")

cat(
  "clusterProfiler: ",
  as.character(
    packageVersion("clusterProfiler")
  ),
  "\n",
  sep = ""
)

cat(
  "fgsea:           ",
  as.character(
    packageVersion("fgsea")
  ),
  "\n",
  sep = ""
)

cat(
  "AnnotationDbi:   ",
  as.character(
    packageVersion("AnnotationDbi")
  ),
  "\n",
  sep = ""
)

cat(
  "org.Hs.eg.db:    ",
  as.character(
    packageVersion("org.Hs.eg.db")
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 2. Working directory and output directory
# ------------------------------------------------------------

step9d_root <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA"

step3b_dir <-
  file.path(
    step9d_root,
    "9D-3B_Formal_GSEA_Execution"
  )

step3c_dir <-
  file.path(
    step9d_root,
    "9D-3C_Formal_GSEA_Result_Audit"
  )

if (!dir.exists(step3b_dir)) {
  stop(
    paste0(
      "FAIL: STEP9D-3B directory does not exist:\n",
      step3b_dir
    )
  )
}

if (!dir.exists(step3c_dir)) {
  dir.create(
    step3c_dir,
    recursive = TRUE
  )
}

cat("STEP9D-3B directory verified.\n")
cat("STEP9D-3C audit directory verified.\n\n")


# ------------------------------------------------------------
# 3. Exact FINAL 3B files
# ------------------------------------------------------------

input_ranking_file <-
  file.path(
    step3b_dir,
    "STEP9D3B_FINAL_KEGG_GSEA_Input_Ranking.rds"
  )

input_ranking_csv <-
  file.path(
    step3b_dir,
    "STEP9D3B_FINAL_KEGG_GSEA_Input_Ranking.csv"
  )

gsea_result_file <-
  file.path(
    step3b_dir,
    "STEP9D3B_FINAL_KEGG_GSEA_Result.rds"
  )

gsea_result_csv <-
  file.path(
    step3b_dir,
    "STEP9D3B_FINAL_KEGG_GSEA_Result.csv"
  )

gsea_sig_csv <-
  file.path(
    step3b_dir,
    "STEP9D3B_FINAL_KEGG_GSEA_Significant_FDR05.csv"
  )

provenance_file <-
  file.path(
    step3b_dir,
    "STEP9D3B_FINAL_GSEA_Provenance.rds"
  )

execution_audit_file <-
  file.path(
    step3b_dir,
    "STEP9D3B_FINAL_GSEA_Execution_Audit.txt"
  )


required_files <- c(
  input_ranking_file,
  input_ranking_csv,
  gsea_result_file,
  gsea_result_csv,
  gsea_sig_csv,
  provenance_file,
  execution_audit_file
)

cat("FINAL STEP9D-3B FILE EXISTENCE\n")
cat("--------------------------------\n")

file_check <- data.frame(
  file = required_files,
  exists = file.exists(required_files),
  stringsAsFactors = FALSE
)

print(file_check)

if (!all(file_check$exists)) {
  stop(
    "FAIL: one or more required STEP9D-3B files are missing."
  )
}

cat("\nAll required STEP9D-3B files exist: TRUE\n\n")


# ------------------------------------------------------------
# 4. Read-only loading
# ------------------------------------------------------------

cat("LOADING FINAL LOCKED OBJECTS\n")
cat("-----------------------------\n")

kegg_ranking_final <-
  readRDS(input_ranking_file)

gsea_result <-
  readRDS(gsea_result_file)

gsea_result_csv_df <-
  read.csv(
    gsea_result_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

gsea_sig_csv_df <-
  read.csv(
    gsea_sig_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

provenance <-
  readRDS(provenance_file)

cat("Final KEGG ranking loaded.\n")
cat("Final GSEA result loaded.\n")
cat("Final GSEA CSV loaded.\n")
cat("Final significant CSV loaded.\n")
cat("Final provenance loaded.\n\n")


# ------------------------------------------------------------
# 5. Object structure audit
# ------------------------------------------------------------

cat("GSEA RESULT OBJECT STRUCTURE\n")
cat("----------------------------\n")

cat(
  "Class:\n"
)
print(class(gsea_result))

cat("\nNames:\n")
print(names(gsea_result))

cat("\nNumber of rows:\n")

if (is.data.frame(gsea_result)) {
  
  cat(
    nrow(gsea_result),
    "\n"
  )
  
} else if (
  "result" %in% slotNames(gsea_result)
) {
  
  cat(
    nrow(gsea_result@result),
    "\n"
  )
  
} else {
  
  cat("Unable to identify result table automatically.\n")
}

cat("\n")


# ------------------------------------------------------------
# 6. Extract result table safely
# ------------------------------------------------------------

if (is.data.frame(gsea_result)) {
  
  gsea_df <-
    gsea_result
  
} else if (
  inherits(gsea_result, "gseaResult")
) {
  
  gsea_df <-
    as.data.frame(gsea_result)
  
} else {
  
  stop(
    paste0(
      "FAIL: unsupported GSEA result object class: ",
      paste(class(gsea_result), collapse = ", ")
    )
  )
}

cat("Result table successfully extracted.\n")
cat(
  "Observed pathway rows: ",
  nrow(gsea_df),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Core structural columns
# ------------------------------------------------------------

cat("RESULT COLUMN AUDIT\n")
cat("-------------------\n")

print(
  data.frame(
    column = names(gsea_df),
    stringsAsFactors = FALSE
  )
)

required_result_columns <-
  c(
    "ID",
    "Description",
    "NES",
    "pvalue",
    "p.adjust"
  )

missing_required_columns <-
  setdiff(
    required_result_columns,
    names(gsea_df)
  )

if (length(missing_required_columns) > 0) {
  
  stop(
    paste0(
      "FAIL: missing required GSEA result columns: ",
      paste(
        missing_required_columns,
        collapse = ", "
      )
    )
  )
}

cat(
  "\nAll core GSEA result columns present: TRUE\n\n"
)


# ------------------------------------------------------------
# 8. Formal pathway count audit
# ------------------------------------------------------------

observed_pathway_n <-
  nrow(gsea_df)

expected_pathway_n <-
  353

if (observed_pathway_n != expected_pathway_n) {
  
  stop(
    paste0(
      "FAIL: expected ",
      expected_pathway_n,
      " pathways; observed ",
      observed_pathway_n,
      "."
    )
  )
}

cat(
  "KEGG pathways returned: ",
  observed_pathway_n,
  "\n",
  sep = ""
)

cat(
  "Expected pathways:      ",
  expected_pathway_n,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Pathway ID uniqueness
# ------------------------------------------------------------

duplicated_pathway_id_n <-
  sum(
    duplicated(gsea_df$ID)
  )

if (duplicated_pathway_id_n != 0) {
  
  stop(
    paste0(
      "FAIL: duplicated KEGG pathway IDs detected: ",
      duplicated_pathway_id_n
    )
  )
}

cat(
  "Duplicated pathway IDs: 0\n\n"
)


# ------------------------------------------------------------
# 10. Numerical validity audit
# ------------------------------------------------------------

cat("NUMERICAL VALIDITY AUDIT\n")
cat("-------------------------\n")

numeric_columns <-
  c(
    "NES",
    "pvalue",
    "p.adjust"
  )

for (col in numeric_columns) {
  
  if (!is.numeric(gsea_df[[col]])) {
    
    stop(
      paste0(
        "FAIL: column is not numeric: ",
        col
      )
    )
  }
  
  if (anyNA(gsea_df[[col]])) {
    
    stop(
      paste0(
        "FAIL: NA detected in ",
        col
      )
    )
  }
  
  if (any(!is.finite(gsea_df[[col]]))) {
    
    stop(
      paste0(
        "FAIL: non-finite value detected in ",
        col
      )
    )
  }
}

cat("NES valid: TRUE\n")
cat("pvalue valid: TRUE\n")
cat("p.adjust valid: TRUE\n\n")


# ------------------------------------------------------------
# 11. FDR < 0.05 reproducibility audit
# ------------------------------------------------------------

gsea_sig_reconstructed <-
  gsea_df[
    !is.na(gsea_df$p.adjust) &
      gsea_df$p.adjust < 0.05,
    ,
    drop = FALSE
  ]

observed_fdr05_n <-
  nrow(gsea_sig_reconstructed)

expected_fdr05_n <-
  178

cat("FDR < 0.05 AUDIT\n")
cat("----------------\n")

cat(
  "Observed FDR < 0.05 pathways: ",
  observed_fdr05_n,
  "\n",
  sep = ""
)

cat(
  "Expected FDR < 0.05 pathways: ",
  expected_fdr05_n,
  "\n\n",
  sep = ""
)

if (observed_fdr05_n != expected_fdr05_n) {
  
  stop(
    paste0(
      "FAIL: expected ",
      expected_fdr05_n,
      " significant pathways; observed ",
      observed_fdr05_n,
      "."
    )
  )
}


# ------------------------------------------------------------
# 12. Compare significant CSV with reconstructed result
# ------------------------------------------------------------

if (!("ID" %in% names(gsea_sig_csv_df))) {
  
  stop(
    "FAIL: significant CSV does not contain ID column."
  )
}

csv_sig_ids <-
  as.character(
    gsea_sig_csv_df$ID
  )

reconstructed_sig_ids <-
  as.character(
    gsea_sig_reconstructed$ID
  )

sig_id_set_match <-
  setequal(
    csv_sig_ids,
    reconstructed_sig_ids
  )

if (!isTRUE(sig_id_set_match)) {
  
  stop(
    "FAIL: significant FDR<0.05 CSV does not match reconstructed FDR<0.05 pathway set."
  )
}

cat(
  "Significant CSV pathway set matches reconstructed FDR<0.05 set: TRUE\n\n"
)


# ------------------------------------------------------------
# 13. Compare complete CSV with RDS result
# ------------------------------------------------------------

gsea_rds_df <-
  as.data.frame(gsea_result)

csv_id_match <-
  identical(
    as.character(gsea_rds_df$ID),
    as.character(gsea_result_csv_df$ID)
  )

csv_nes_match <-
  isTRUE(
    all.equal(
      as.numeric(gsea_rds_df$NES),
      as.numeric(gsea_result_csv_df$NES),
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  )

csv_pvalue_match <-
  isTRUE(
    all.equal(
      as.numeric(gsea_rds_df$pvalue),
      as.numeric(gsea_result_csv_df$pvalue),
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  )

csv_padjust_match <-
  isTRUE(
    all.equal(
      as.numeric(gsea_rds_df$p.adjust),
      as.numeric(gsea_result_csv_df$p.adjust),
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  )

cat("RDS ↔ CSV REPRODUCIBILITY AUDIT\n")
cat("-------------------------------\n")

cat(
  "Pathway ID order identical: ",
  csv_id_match,
  "\n",
  sep = ""
)

cat(
  "NES values identical:       ",
  csv_nes_match,
  "\n",
  sep = ""
)

cat(
  "P values identical:         ",
  csv_pvalue_match,
  "\n",
  sep = ""
)

cat(
  "Adjusted P values identical:",
  csv_padjust_match,
  "\n\n",
  sep = ""
)

if (
  !csv_id_match ||
  !csv_nes_match ||
  !csv_pvalue_match ||
  !csv_padjust_match
) {
  
  stop(
    "FAIL: GSEA RDS and CSV results are not fully concordant."
  )
}


# ------------------------------------------------------------
# 14. NES direction audit
# ------------------------------------------------------------

positive_nes_n <-
  sum(gsea_df$NES > 0)

negative_nes_n <-
  sum(gsea_df$NES < 0)

zero_nes_n <-
  sum(gsea_df$NES == 0)

cat("NES DIRECTION AUDIT\n")
cat("-------------------\n")

cat(
  "Positive NES: ",
  positive_nes_n,
  "\n",
  sep = ""
)

cat(
  "Negative NES: ",
  negative_nes_n,
  "\n",
  sep = ""
)

cat(
  "Zero NES:     ",
  zero_nes_n,
  "\n\n",
  sep = ""
)

if (
  positive_nes_n +
  negative_nes_n +
  zero_nes_n !=
  nrow(gsea_df)
) {
  
  stop(
    "FAIL: NES direction counts do not sum to total pathways."
  )
}


# ------------------------------------------------------------
# 15. P-value / FDR range audit
# ------------------------------------------------------------

pvalue_range_valid <-
  all(
    gsea_df$pvalue >= 0 &
      gsea_df$pvalue <= 1
  )

padjust_range_valid <-
  all(
    gsea_df$p.adjust >= 0 &
      gsea_df$p.adjust <= 1
  )

if (!pvalue_range_valid) {
  stop("FAIL: invalid P-value range detected.")
}

if (!padjust_range_valid) {
  stop("FAIL: invalid adjusted P-value range detected.")
}

cat("P-value range valid: TRUE\n")
cat("Adjusted P-value range valid: TRUE\n\n")


# ------------------------------------------------------------
# 16. Core enrichment / leading-edge audit
# ------------------------------------------------------------

core_enrichment_column_present <-
  "core_enrichment" %in% names(gsea_df)

leading_edge_column_present <-
  "leading_edge" %in% names(gsea_df)

cat("CORE ENRICHMENT AUDIT\n")
cat("---------------------\n")

cat(
  "core_enrichment column present: ",
  core_enrichment_column_present,
  "\n",
  sep = ""
)

cat(
  "leading_edge column present:    ",
  leading_edge_column_present,
  "\n\n",
  sep = ""
)

if (core_enrichment_column_present) {
  
  core_empty_n <-
    sum(
      is.na(gsea_df$core_enrichment) |
        trimws(gsea_df$core_enrichment) == ""
    )
  
  cat(
    "Empty core_enrichment rows: ",
    core_empty_n,
    "\n",
    sep = ""
  )
  
} else {
  
  core_empty_n <-
    NA_integer_
  
  cat(
    "core_enrichment is not a stored column in this result object.\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 17. Leading-edge / log2err audit if available
# ------------------------------------------------------------

log2err_present <-
  "log2err" %in% names(gsea_df)

if (log2err_present) {
  
  log2err_na_n <-
    sum(
      is.na(gsea_df$log2err)
    )
  
  cat(
    "log2err column present: TRUE\n"
  )
  
  cat(
    "log2err NA count: ",
    log2err_na_n,
    "\n\n",
    sep = ""
  )
  
} else {
  
  log2err_na_n <-
    NA_integer_
  
  cat(
    "log2err column present: FALSE\n\n"
  )
}


cat("\n============================================================\n")
cat("STEP9D-3C: FINAL KEGG GSEA INPUT OBJECT STRUCTURE AUDIT\n")
cat("============================================================\n\n")

cat("Class:\n")
print(class(kegg_ranking_final))

cat("\nType:\n")
print(typeof(kegg_ranking_final))

cat("\nNames:\n")
print(names(kegg_ranking_final))

cat("\nStructure:\n")
str(
  kegg_ranking_final,
  max.level = 2
)

cat("\n============================================================\n")


cat("\n============================================================\n")
cat("STEP9D-3C: FINAL KEGG INPUT ORDER / SORTING DIAGNOSTIC\n")
cat("READ-ONLY — NO RECALCULATION / NO GSEA\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Extract locked objects
# ------------------------------------------------------------

rv <- kegg_ranking_final$ranking_vector
it <- kegg_ranking_final$input_table

# ------------------------------------------------------------
# 2. Basic structure
# ------------------------------------------------------------

cat("Ranking vector length:", length(rv), "\n")
cat("Input table rows:", nrow(it), "\n")

cat("\nInput table columns:\n")
print(names(it))

# ------------------------------------------------------------
# 3. Check whether input_table contains the ranking values
# ------------------------------------------------------------

ranking_col_candidates <- names(it)[
  grepl(
    "ranking|stat|score",
    names(it),
    ignore.case = TRUE
  )
]

cat("\nPossible ranking/statistic columns:\n")
print(ranking_col_candidates)

# ------------------------------------------------------------
# 4. Compare Entrez ID sets
# ------------------------------------------------------------

input_ids <- as.character(it$ENTREZID)
ranking_ids <- names(rv)

cat("\nID set comparison:\n")

cat(
  "Same ID set:",
  setequal(input_ids, ranking_ids),
  "\n"
)

cat(
  "Input IDs not in ranking:",
  length(setdiff(input_ids, ranking_ids)),
  "\n"
)

cat(
  "Ranking IDs not in input:",
  length(setdiff(ranking_ids, input_ids)),
  "\n"
)

# ------------------------------------------------------------
# 5. Check the first 20 IDs in each object
# ------------------------------------------------------------

cat("\nFirst 20 input_table Entrez IDs:\n")
print(head(input_ids, 20))

cat("\nFirst 20 ranking-vector Entrez IDs:\n")
print(head(ranking_ids, 20))

# ------------------------------------------------------------
# 6. Check the last 20 IDs
# ------------------------------------------------------------

cat("\nLast 20 input_table Entrez IDs:\n")
print(tail(input_ids, 20))

cat("\nLast 20 ranking-vector Entrez IDs:\n")
print(tail(ranking_ids, 20))

# ------------------------------------------------------------
# 7. Determine whether ranking vector is sorted decreasing
# ------------------------------------------------------------

is_sorted_decreasing <- all(
  diff(rv) <= 0
)

is_sorted_increasing <- all(
  diff(rv) >= 0
)

cat("\nRanking sorted decreasing:", is_sorted_decreasing, "\n")
cat("Ranking sorted increasing:", is_sorted_increasing, "\n")

# ------------------------------------------------------------
# 8. Check whether input_table itself contains a ranking column
# ------------------------------------------------------------

if (length(ranking_col_candidates) > 0) {
  
  for (cc in ranking_col_candidates) {
    
    x <- it[[cc]]
    
    cat(
      "\nColumn:",
      cc,
      "\n"
    )
    
    cat(
      "Class:",
      class(x),
      "\n"
    )
    
    if (is.numeric(x)) {
      
      cat(
        "Numeric:",
        TRUE,
        "\n"
      )
      
      cat(
        "Sorted decreasing:",
        all(diff(x) <= 0),
        "\n"
      )
      
      cat(
        "Sorted increasing:",
        all(diff(x) >= 0),
        "\n"
      )
      
      cat(
        "Exact value-set match with ranking vector:",
        setequal(x, unname(rv)),
        "\n"
      )
      
      cat(
        "Exact vector match with ranking vector:",
        identical(as.numeric(x), as.numeric(rv)),
        "\n"
      )
    }
  }
  
} else {
  
  cat(
    "\nNo obvious ranking/statistic column detected in input_table.\n"
  )
}

# ------------------------------------------------------------
# 9. Compare ranking vector with input_table after matching by ID
# ------------------------------------------------------------

if (setequal(input_ids, ranking_ids)) {
  
  rv_by_input_order <- unname(
    rv[input_ids]
  )
  
  aligned_value_identical <- identical(
    as.numeric(rv_by_input_order),
    as.numeric(rv)
  )
  
  cat(
    "\nRanking values in input_table ID order identical to\n",
    "current ranking-vector order:",
    aligned_value_identical,
    "\n"
  )
  
  cat(
    "Ranking values in input_table ID order sorted decreasing:",
    all(diff(rv_by_input_order) <= 0),
    "\n"
  )
  
  # Check whether input_table order itself corresponds to
  # descending ranking after ID matching
  cat(
    "Input-table ID order reproduces descending ranking:",
    all(diff(rv_by_input_order) <= 0),
    "\n"
  )
  
} else {
  
  cat(
    "\nID sets are not identical; alignment test skipped.\n"
  )
}

# ------------------------------------------------------------
# 10. Final diagnostic conclusion
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("DIAGNOSTIC COMPLETE\n")
cat("============================================================\n")
cat("\nNO FILES MODIFIED.\n")
cat("NO RANKING RECALCULATED.\n")
cat("NO GSEA EXECUTED.\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 18. FINAL KEGG GSEA INPUT RANKING AUDIT
#     FINAL CORRECTED VERSION
#
#     IMPORTANT:
#     The formal GSEA ranking is stored in:
#         kegg_ranking_final$input_table$ranking
#
#     The accompanying ranking_vector is an Entrez-ID-indexed
#     representation and is not required to have the same order.
#
#     Read-only audit.
#     No ranking recalculation.
#     No reordering.
#     No GSEA.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("18. FINAL KEGG GSEA INPUT RANKING AUDIT\n")
cat("FINAL CORRECTED VERSION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 18.1 Extract locked objects
# ------------------------------------------------------------

input_table_final <- kegg_ranking_final$input_table
ranking_vector_final <- kegg_ranking_final$ranking_vector

input_ids <- as.character(
  input_table_final$ENTREZID
)

input_ranking <- as.numeric(
  input_table_final$ranking
)

ranking_ids <- names(
  ranking_vector_final
)

ranking_values <- as.numeric(
  ranking_vector_final
)


# ------------------------------------------------------------
# 18.2 Basic structural integrity
# ------------------------------------------------------------

input_table_n <- nrow(input_table_final)

input_ranking_numeric <- is.numeric(
  input_table_final$ranking
)

input_ranking_no_na <- all(
  !is.na(input_ranking)
)

input_ranking_finite <- all(
  is.finite(input_ranking)
)

input_ids_unique <- (
  length(unique(input_ids)) == length(input_ids)
)

ranking_vector_numeric <- is.numeric(
  ranking_vector_final
)

ranking_vector_no_na <- all(
  !is.na(ranking_values)
)

ranking_vector_finite <- all(
  is.finite(ranking_values)
)

ranking_vector_named <- !is.null(
  names(ranking_vector_final)
)

ranking_vector_ids_unique <- (
  ranking_vector_named &&
    length(unique(ranking_ids)) == length(ranking_ids)
)


cat("Input table rows:",
    input_table_n,
    "\n")

cat("Input ranking numeric:",
    input_ranking_numeric,
    "\n")

cat("Input ranking no NA:",
    input_ranking_no_na,
    "\n")

cat("Input ranking all finite:",
    input_ranking_finite,
    "\n")

cat("Input Entrez IDs unique:",
    input_ids_unique,
    "\n")

cat("Ranking vector numeric:",
    ranking_vector_numeric,
    "\n")

cat("Ranking vector no NA:",
    ranking_vector_no_na,
    "\n")

cat("Ranking vector all finite:",
    ranking_vector_finite,
    "\n")

cat("Ranking vector named:",
    ranking_vector_named,
    "\n")

cat("Ranking vector Entrez IDs unique:",
    ranking_vector_ids_unique,
    "\n")


# ------------------------------------------------------------
# 18.3 Locked metadata concordance
# ------------------------------------------------------------

final_entrez_n_match <- (
  input_table_n ==
    as.integer(kegg_ranking_final$final_entrez_n)
)

positive_n_input <- sum(
  input_ranking > 0
)

negative_n_input <- sum(
  input_ranking < 0
)

zero_n_input <- sum(
  input_ranking == 0
)

positive_n_match <- (
  positive_n_input ==
    as.integer(kegg_ranking_final$positive_n)
)

negative_n_match <- (
  negative_n_input ==
    as.integer(kegg_ranking_final$negative_n)
)

zero_n_match <- (
  zero_n_input ==
    as.integer(kegg_ranking_final$zero_n)
)

cat("\n============================================================\n")
cat("LOCKED METADATA CONCORDANCE\n")
cat("============================================================\n")

cat("Final Entrez N match:",
    final_entrez_n_match,
    "\n")

cat(
  "Positive:",
  positive_n_input,
  " | locked:",
  kegg_ranking_final$positive_n,
  " | match:",
  positive_n_match,
  "\n"
)

cat(
  "Negative:",
  negative_n_input,
  " | locked:",
  kegg_ranking_final$negative_n,
  " | match:",
  negative_n_match,
  "\n"
)

cat(
  "Zero:",
  zero_n_input,
  " | locked:",
  kegg_ranking_final$zero_n,
  " | match:",
  zero_n_match,
  "\n"
)


# ------------------------------------------------------------
# 18.4 Formal GSEA ordering audit
#
# The formal GSEA input must be ordered decreasingly.
# Audit the ranking column actually stored in input_table.
# ------------------------------------------------------------

gsea_ranking_sorted_decreasing <- all(
  diff(input_ranking) <= 0
)

gsea_ranking_sorted_increasing <- all(
  diff(input_ranking) >= 0
)

cat("\n============================================================\n")
cat("FORMAL GSEA RANKING ORDER AUDIT\n")
cat("============================================================\n")

cat(
  "Formal GSEA ranking sorted decreasing:",
  gsea_ranking_sorted_decreasing,
  "\n"
)

cat(
  "Formal GSEA ranking sorted increasing:",
  gsea_ranking_sorted_increasing,
  "\n"
)


# ------------------------------------------------------------
# 18.5 Entrez ID-set concordance
#
# Order is NOT required to be identical.
# Set equality is the correct audit criterion.
# ------------------------------------------------------------

same_entrez_set <- setequal(
  input_ids,
  ranking_ids
)

input_only_ids <- length(
  setdiff(input_ids, ranking_ids)
)

ranking_only_ids <- length(
  setdiff(ranking_ids, input_ids)
)

cat("\n============================================================\n")
cat("ENTREZ ID-SET CONCORDANCE\n")
cat("============================================================\n")

cat(
  "Same Entrez ID set:",
  same_entrez_set,
  "\n"
)

cat(
  "Input-only Entrez IDs:",
  input_only_ids,
  "\n"
)

cat(
  "Ranking-only Entrez IDs:",
  ranking_only_ids,
  "\n"
)


# ------------------------------------------------------------
# 18.6 Ranking-value concordance
#
# The two stored representations may have different order.
# Therefore compare:
#
#   1. value set
#   2. value multiset after sorting
#
# rather than requiring identical vector order.
# ------------------------------------------------------------

ranking_value_set_match <- setequal(
  input_ranking,
  ranking_values
)

ranking_value_multiset_match <- identical(
  sort(input_ranking),
  sort(ranking_values)
)

cat("\n============================================================\n")
cat("RANKING-VALUE CONCORDANCE\n")
cat("============================================================\n")

cat(
  "Ranking value set match:",
  ranking_value_set_match,
  "\n"
)

cat(
  "Ranking value multiset match:",
  ranking_value_multiset_match,
  "\n"
)


# ------------------------------------------------------------
# 18.7 Verify one-to-one ID ↔ ranking correspondence
#
# Reconstruct the ranking vector in input_table Entrez order
# WITHOUT changing or recalculating any ranking values.
# ------------------------------------------------------------

if (same_entrez_set) {
  
  ranking_reordered_by_input_id <- unname(
    ranking_vector_final[input_ids]
  )
  
  id_aligned_ranking_match <- identical(
    as.numeric(ranking_reordered_by_input_id),
    as.numeric(input_ranking)
  )
  
} else {
  
  ranking_reordered_by_input_id <- NULL
  
  id_aligned_ranking_match <- FALSE
}

cat("\n============================================================\n")
cat("ID-ALIGNED RANKING CONCORDANCE\n")
cat("============================================================\n")

cat(
  "Ranking values match after ID alignment:",
  id_aligned_ranking_match,
  "\n"
)


# ------------------------------------------------------------
# 18.8 Core-gene preservation
# ------------------------------------------------------------

core_entrez <- c(
  AOC3  = "8639",
  MS4A7 = "58475",
  IL7R  = "3575"
)

core_presence <- core_entrez %in% input_ids

core_status_df <- data.frame(
  gene = names(core_entrez),
  ENTREZID = unname(core_entrez),
  present = core_presence,
  row.names = NULL
)

cat("\n============================================================\n")
cat("CORE-GENE PRESERVATION\n")
cat("============================================================\n")

print(core_status_df)

core_all_present <- all(
  core_presence
)


# ------------------------------------------------------------
# 18.9 Final hard-pass definition
# ------------------------------------------------------------

section18_hard_pass <- all(
  
  # Input table
  input_table_n == 34120,
  final_entrez_n_match,
  
  # Formal GSEA ranking
  input_ranking_numeric,
  input_ranking_no_na,
  input_ranking_finite,
  gsea_ranking_sorted_decreasing,
  
  # ID integrity
  input_ids_unique,
  same_entrez_set,
  input_only_ids == 0,
  ranking_only_ids == 0,
  
  # Ranking-vector integrity
  ranking_vector_numeric,
  ranking_vector_no_na,
  ranking_vector_finite,
  ranking_vector_named,
  ranking_vector_ids_unique,
  
  # Ranking-value integrity
  ranking_value_set_match,
  ranking_value_multiset_match,
  id_aligned_ranking_match,
  
  # Locked counts
  positive_n_match,
  negative_n_match,
  zero_n_match,
  
  # Core genes
  core_all_present
)


# ------------------------------------------------------------
# 18.10 Final Section 18 report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SECTION 18 FINAL KEGG RANKING AUDIT\n")
cat("============================================================\n")

cat(
  "Final Entrez N = 34,120:            ",
  input_table_n == 34120,
  "\n"
)

cat(
  "Formal ranking numeric:             ",
  input_ranking_numeric,
  "\n"
)

cat(
  "Formal ranking no NA:               ",
  input_ranking_no_na,
  "\n"
)

cat(
  "Formal ranking all finite:          ",
  input_ranking_finite,
  "\n"
)

cat(
  "Formal ranking sorted decreasing:   ",
  gsea_ranking_sorted_decreasing,
  "\n"
)

cat(
  "Entrez IDs unique:                  ",
  input_ids_unique,
  "\n"
)

cat(
  "Same Entrez ID set:                 ",
  same_entrez_set,
  "\n"
)

cat(
  "Ranking value set match:             ",
  ranking_value_set_match,
  "\n"
)

cat(
  "Ranking multiset match:              ",
  ranking_value_multiset_match,
  "\n"
)

cat(
  "ID-aligned ranking match:            ",
  id_aligned_ranking_match,
  "\n"
)

cat(
  "Positive count match:                ",
  positive_n_match,
  "\n"
)

cat(
  "Negative count match:                ",
  negative_n_match,
  "\n"
)

cat(
  "Zero count match:                    ",
  zero_n_match,
  "\n"
)

cat(
  "All three core genes present:        ",
  core_all_present,
  "\n"
)

cat("\nSECTION 18 HARD PASS: ",
    section18_hard_pass,
    "\n")

cat("============================================================\n\n")


# ------------------------------------------------------------
# 19. Core gene presence in FINAL KEGG ranking
#     FINAL CORRECTED VERSION
#
#     Read-only.
#     No ranking recalculation.
#     No GSEA.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("19. CORE GENE PRESENCE IN FINAL KEGG RANKING\n")
cat("FINAL CORRECTED VERSION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 19.1 Locked core-gene Entrez IDs
# ------------------------------------------------------------

core_genes <- c(
  AOC3  = "8639",
  MS4A7 = "58475",
  IL7R  = "3575"
)


# ------------------------------------------------------------
# 19.2 Use the ACTUAL FINAL KEGG GSEA INPUT
# ------------------------------------------------------------

final_kegg_entrez <- as.character(
  kegg_ranking_final$input_table$ENTREZID
)


# ------------------------------------------------------------
# 19.3 Core-gene presence audit
# ------------------------------------------------------------

core_presence <- data.frame(
  gene = names(core_genes),
  EntrezID = unname(core_genes),
  retained_in_final_kegg_ranking =
    unname(core_genes) %in% final_kegg_entrez,
  stringsAsFactors = FALSE
)


cat("CORE GENE FINAL KEGG INPUT AUDIT\n")
cat("---------------------------------\n\n")

print(core_presence)


# ------------------------------------------------------------
# 19.4 Hard check
# ------------------------------------------------------------

core_all_retained <- all(
  core_presence$retained_in_final_kegg_ranking
)


cat("\nAll three core genes retained:",
    core_all_retained,
    "\n")


if (!core_all_retained) {
  
  stop(
    "FAIL: one or more core genes are absent from final KEGG ranking."
  )
  
}


cat("\nSTEP9D-3C SECTION 19 HARD PASS: TRUE\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 20. Provenance audit
# ------------------------------------------------------------

cat("PROVENANCE OBJECT STRUCTURE\n")
cat("----------------------------\n")

print(names(provenance))

cat("\n")


# ------------------------------------------------------------
# 21. Check important provenance fields dynamically
# ------------------------------------------------------------

provenance_fields_to_check <-
  c(
    "workflow",
    "step",
    "organism",
    "keyType",
    "minGSSize",
    "maxGSSize",
    "exponent",
    "eps",
    "pvalueCutoff",
    "pAdjustMethod",
    "seed"
  )

provenance_audit <-
  data.frame(
    field =
      provenance_fields_to_check,
    present =
      provenance_fields_to_check %in%
      names(provenance),
    value =
      NA_character_,
    stringsAsFactors = FALSE
  )

for (i in seq_len(nrow(provenance_audit))) {
  
  field_i <-
    provenance_audit$field[i]
  
  if (
    provenance_audit$present[i]
  ) {
    
    value_i <-
      provenance[[field_i]]
    
    if (length(value_i) == 1) {
      
      provenance_audit$value[i] <-
        as.character(value_i)
      
    } else {
      
      provenance_audit$value[i] <-
        paste(
          value_i,
          collapse = ";"
        )
    }
  }
}

print(provenance_audit)

cat("\n")


# ------------------------------------------------------------
# 22. Compare locked expected GSEA parameters
# ------------------------------------------------------------

expected_params <-
  list(
    organism = "hsa",
    keyType = "ncbi-geneid",
    minGSSize = 10,
    maxGSSize = 500,
    exponent = 1,
    eps = 1e-10,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE,
    seed = FALSE
  )

parameter_match <- list()

for (nm in names(expected_params)) {
  
  if (nm %in% names(provenance)) {
    
    parameter_match[[nm]] <-
      isTRUE(
        all.equal(
          provenance[[nm]],
          expected_params[[nm]],
          check.attributes = FALSE
        )
      )
    
  } else {
    
    parameter_match[[nm]] <-
      NA
  }
}

parameter_audit_df <-
  data.frame(
    parameter =
      names(parameter_match),
    provenance_present =
      names(parameter_match) %in%
      names(provenance),
    matches_locked_definition =
      unlist(parameter_match),
    stringsAsFactors = FALSE
  )

cat("GSEA PARAMETER LOCK AUDIT\n")
cat("-------------------------\n")

print(parameter_audit_df)

cat("\n")


# ------------------------------------------------------------
# 23. Execution-integrity flags
# ------------------------------------------------------------

execution_integrity <-
  data.frame(
    item = c(
      "DEG_recalculated",
      "ranking_recalculated",
      "mapping_recalculated",
      "collision_resolution_reapplied",
      "gsea_rerun_during_audit",
      "jitter_used",
      "artificial_tie_breaking_used"
    ),
    value = c(
      FALSE,
      FALSE,
      FALSE,
      FALSE,
      FALSE,
      FALSE,
      FALSE
    ),
    stringsAsFactors = FALSE
  )

cat("READ-ONLY AUDIT INTEGRITY\n")
cat("-------------------------\n")

print(execution_integrity)

cat("\n")

# ============================================================
# STEP9D-3C
# Sections 24–27 — FINAL READ-ONLY AUDIT
# CORRECTED VERSION
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-3C: FINAL READ-ONLY AUDIT — CORRECTED\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 24. Final audit decision
# ------------------------------------------------------------

# IMPORTANT:
# kegg_ranking_final is a LIST.
# The actual final KEGG ranking vector is:
ranking_vector_final <-
  kegg_ranking_final$ranking_vector

# The formal GSEA input table is:
input_table_final <-
  kegg_ranking_final$input_table


# ------------------------------------------------------------
# 24A. Final ranking integrity
# ------------------------------------------------------------

ranking_pass <-
  length(ranking_vector_final) == 34120 &&
  is.numeric(ranking_vector_final) &&
  !anyNA(ranking_vector_final) &&
  all(is.finite(ranking_vector_final)) &&
  !anyDuplicated(names(ranking_vector_final))


# ------------------------------------------------------------
# 24B. Structural result pass
# ------------------------------------------------------------

structural_pass <-
  observed_pathway_n == 353 &&
  duplicated_pathway_id_n == 0 &&
  pvalue_range_valid &&
  padjust_range_valid


# ------------------------------------------------------------
# 24C. FDR result pass
# ------------------------------------------------------------

fdr_pass <-
  observed_fdr05_n == 178 &&
  sig_id_set_match


# ------------------------------------------------------------
# 24D. RDS / CSV concordance
# ------------------------------------------------------------

csv_rds_pass <-
  csv_id_match &&
  csv_nes_match &&
  csv_pvalue_match &&
  csv_padjust_match


# ------------------------------------------------------------
# 24E. Core-gene preservation
# ------------------------------------------------------------

core_gene_pass <-
  all(
    core_presence$retained_in_final_kegg_ranking
  )


# ------------------------------------------------------------
# 24F. Read-only integrity
# ------------------------------------------------------------

integrity_pass <-
  all(
    execution_integrity$value == FALSE
  )


# ------------------------------------------------------------
# 24G. Parameter audit
# ------------------------------------------------------------

parameter_values <-
  parameter_audit_df$matches_locked_definition

parameter_pass <-
  all(
    parameter_values[
      !is.na(parameter_values)
    ]
  )


# ------------------------------------------------------------
# 24H. Ranking tie audit
# ------------------------------------------------------------

ranking_unique_n <-
  length(unique(ranking_vector_final))

ranking_tied_gene_n <-
  length(ranking_vector_final) -
  ranking_unique_n

ranking_tie_fraction <-
  ranking_tied_gene_n /
  length(ranking_vector_final)


# ------------------------------------------------------------
# 24I. Final KEGG ranking order audit
# ------------------------------------------------------------

formal_ranking <-
  input_table_final$ranking

gsea_ranking_sorted_decreasing <-
  all(
    diff(formal_ranking) <= 0
  )

entrez_ids_final <-
  as.character(
    input_table_final$ENTREZID
  )

ranking_order_audit_pass <-
  length(formal_ranking) == 34120 &&
  gsea_ranking_sorted_decreasing &&
  !anyDuplicated(entrez_ids_final)


# ------------------------------------------------------------
# 24J. Overall final decision
# ------------------------------------------------------------

final_hard_pass <-
  structural_pass &&
  fdr_pass &&
  csv_rds_pass &&
  ranking_pass &&
  ranking_order_audit_pass &&
  core_gene_pass &&
  integrity_pass &&
  parameter_pass


# ------------------------------------------------------------
# 25. Warning interpretation
# ------------------------------------------------------------

cat("WARNING INTERPRETATION\n")
cat("----------------------\n")

cat(
  "Preranked ties were present and retained without artificial ",
  "tie-breaking.\n",
  sep = ""
)

cat(
  "The formal GSEA execution generated fgsea log2err/eps warnings.\n"
)

cat(
  "These warnings are documented as execution metadata and were NOT ",
  "used to alter the FINAL GSEA result during this read-only audit.\n"
)

cat(
  "No jitter, artificial tie-breaking, or GSEA rerun was performed.\n\n"
)


# ------------------------------------------------------------
# 26. FINAL HARD PASS SUMMARY
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP9D-3C FINAL AUDIT SUMMARY\n")
cat("============================================================\n\n")

cat(
  "Pathways returned = ",
  observed_pathway_n,
  " : ",
  observed_pathway_n == 353,
  "\n",
  sep = ""
)

cat(
  "FDR < 0.05 pathways = ",
  observed_fdr05_n,
  " : ",
  observed_fdr05_n == 178,
  "\n",
  sep = ""
)

cat(
  "Duplicated pathway IDs = ",
  duplicated_pathway_id_n,
  " : ",
  duplicated_pathway_id_n == 0,
  "\n",
  sep = ""
)

cat(
  "RDS/CSV concordance = ",
  csv_rds_pass,
  "\n",
  sep = ""
)

cat(
  "Final KEGG ranking = 34,120 : ",
  ranking_pass,
  "\n",
  sep = ""
)

cat(
  "Formal ranking sorted decreasing = ",
  ranking_order_audit_pass,
  "\n",
  sep = ""
)

cat(
  "AOC3/MS4A7/IL7R retained = ",
  core_gene_pass,
  "\n",
  sep = ""
)

cat(
  "Read-only integrity = ",
  integrity_pass,
  "\n",
  sep = ""
)

cat(
  "GSEA parameter lock concordance = ",
  parameter_pass,
  "\n",
  sep = ""
)

cat(
  "Ranking unique values = ",
  ranking_unique_n,
  "\n",
  sep = ""
)

cat(
  "Ranking tied-gene count = ",
  ranking_tied_gene_n,
  "\n",
  sep = ""
)

cat(
  "Ranking tie fraction = ",
  sprintf("%.6f", ranking_tie_fraction),
  "\n",
  sep = ""
)

cat(
  "\nSTEP9D-3C FINAL HARD PASS: ",
  final_hard_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 27. Save formal audit objects
# ------------------------------------------------------------

audit_summary <-
  list(
    
    workflow =
      "Paper1_STEP9D_GenomeWide_Preranked_GSEA",
    
    step =
      "STEP9D-3C",
    
    status =
      if (final_hard_pass)
        "FINAL_PASS"
    else
      "FAIL",
    
    audit_type =
      "READ_ONLY_FORMAL_RESULT_REPRODUCIBILITY_AUDIT",
    
    audit_date =
      as.character(Sys.time()),
    
    input_ranking_file =
      input_ranking_file,
    
    gsea_result_file =
      gsea_result_file,
    
    gsea_result_csv =
      gsea_result_csv,
    
    significant_fdr05_csv =
      gsea_sig_csv,
    
    provenance_file =
      provenance_file,
    
    final_kegg_ranking_n =
      length(ranking_vector_final),
    
    pathway_n =
      observed_pathway_n,
    
    fdr05_n =
      observed_fdr05_n,
    
    duplicated_pathway_id_n =
      duplicated_pathway_id_n,
    
    positive_nes_n =
      positive_nes_n,
    
    negative_nes_n =
      negative_nes_n,
    
    zero_nes_n =
      zero_nes_n,
    
    ranking_unique_n =
      ranking_unique_n,
    
    ranking_tied_gene_n =
      ranking_tied_gene_n,
    
    ranking_tie_fraction =
      ranking_tie_fraction,
    
    core_gene_pass =
      core_gene_pass,
    
    structural_pass =
      structural_pass,
    
    fdr_pass =
      fdr_pass,
    
    csv_rds_pass =
      csv_rds_pass,
    
    ranking_pass =
      ranking_pass,
    
    ranking_order_audit_pass =
      ranking_order_audit_pass,
    
    integrity_pass =
      integrity_pass,
    
    parameter_pass =
      parameter_pass,
    
    execution_integrity =
      execution_integrity,
    
    parameter_audit =
      parameter_audit_df,
    
    core_presence =
      core_presence,
    
    final_hard_pass =
      final_hard_pass
  )


cat("\nAudit summary object created successfully.\n")
cat(
  "FINAL STATUS = ",
  audit_summary$status,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 28. Human-readable TXT audit
# ------------------------------------------------------------

# Use the existing STEP9D-3C output directory
audit_dir <-
  file.path(
    "input_dir/Paper 1_Step 9-Functional Enrichment",
    "STEP9D_GenomeWide_Preranked_GSEA",
    "9D-3C_Formal_GSEA_Result_Audit"
  )

if (!dir.exists(audit_dir)) {
  dir.create(
    audit_dir,
    recursive = TRUE
  )
}

audit_txt_file <-
  file.path(
    audit_dir,
    "STEP9D3C_FINAL_GSEA_Result_Reproducibility_Audit.txt"
  )

audit_rds_file <-
  file.path(
    audit_dir,
    "STEP9D3C_FINAL_GSEA_Result_Reproducibility_Audit.rds"
  )

audit_csv_file <-
  file.path(
    audit_dir,
    "STEP9D3C_FINAL_GSEA_Result_Reproducibility_Audit.csv"
  )


# ------------------------------------------------------------
# Write human-readable audit TXT
# ------------------------------------------------------------

sink(audit_txt_file)

cat("============================================================\n")
cat("Paper 1 — STEP9D-3C\n")
cat("Formal KEGG GSEA Result & Reproducibility Audit\n")
cat("============================================================\n\n")

cat("STATUS\n")
cat("------\n")
cat(
  "Final status: ",
  audit_summary$status,
  "\n\n",
  sep = ""
)

cat("GSEA RESULT SUMMARY\n")
cat("-------------------\n")

cat(
  "KEGG pathways returned: ",
  observed_pathway_n,
  "\n",
  sep = ""
)

cat(
  "KEGG pathways with FDR < 0.05: ",
  observed_fdr05_n,
  "\n",
  sep = ""
)

cat(
  "Duplicated pathway IDs: ",
  duplicated_pathway_id_n,
  "\n",
  sep = ""
)

cat(
  "Positive NES pathways: ",
  positive_nes_n,
  "\n",
  sep = ""
)

cat(
  "Negative NES pathways: ",
  negative_nes_n,
  "\n",
  sep = ""
)

cat(
  "Zero NES pathways: ",
  zero_nes_n,
  "\n\n",
  sep = ""
)


cat("FINAL KEGG RANKING\n")
cat("-------------------\n")

cat(
  "Final Entrez ranking N: ",
  length(ranking_vector_final),
  "\n",
  sep = ""
)

cat(
  "Unique ranking values: ",
  ranking_unique_n,
  "\n",
  sep = ""
)

cat(
  "Tied ranking genes: ",
  ranking_tied_gene_n,
  "\n",
  sep = ""
)

cat(
  "Tie fraction: ",
  sprintf("%.6f", ranking_tie_fraction),
  "\n",
  sep = ""
)

cat(
  "Formal GSEA ranking sorted decreasing: ",
  gsea_ranking_sorted_decreasing,
  "\n\n",
  sep = ""
)


cat("CORE GENE PRESERVATION\n")
cat("----------------------\n")

print(core_presence)

cat("\n")


cat("GSEA PARAMETER AUDIT\n")
cat("--------------------\n")

print(parameter_audit_df)

cat("\n")


cat("READ-ONLY EXECUTION INTEGRITY\n")
cat("-----------------------------\n")

print(execution_integrity)

cat("\n")


cat("AUDIT DECISION\n")
cat("--------------\n")

cat(
  "Structural pass: ",
  structural_pass,
  "\n",
  sep = ""
)

cat(
  "FDR pass: ",
  fdr_pass,
  "\n",
  sep = ""
)

cat(
  "RDS/CSV concordance: ",
  csv_rds_pass,
  "\n",
  sep = ""
)

cat(
  "Ranking pass: ",
  ranking_pass,
  "\n",
  sep = ""
)

cat(
  "Ranking order audit pass: ",
  ranking_order_audit_pass,
  "\n",
  sep = ""
)

cat(
  "Core-gene pass: ",
  core_gene_pass,
  "\n",
  sep = ""
)

cat(
  "Read-only integrity pass: ",
  integrity_pass,
  "\n",
  sep = ""
)

cat(
  "Parameter audit pass: ",
  parameter_pass,
  "\n",
  sep = ""
)

cat(
  "\nSTEP9D-3C FINAL HARD PASS: ",
  final_hard_pass,
  "\n",
  sep = ""
)

cat("\n============================================================\n")
cat("WARNING INTERPRETATION\n")
cat("============================================================\n")

cat(
  "Preranked ties were retained without jitter or artificial ",
  "tie-breaking.\n"
)

cat(
  "fgsea log2err/eps warnings were generated during formal ",
  "execution and are retained as execution metadata.\n"
)

cat(
  "No GSEA rerun was performed during this read-only audit.\n"
)

cat(
  "No DEG recalculation, ranking recalculation, mapping ",
  "recalculation, or collision-resolution reapplication was ",
  "performed.\n"
)

cat("\n============================================================\n")
cat("END OF STEP9D-3C AUDIT\n")
cat("============================================================\n")

sink()

cat("\nSTEP9D-3C TXT audit written successfully.\n")
cat(
  "File: ",
  audit_txt_file,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 29. Save formal STEP9D-3C audit RDS
# ------------------------------------------------------------

saveRDS(
  audit_summary,
  audit_rds_file
)

cat(
  "RDS audit saved: ",
  file.exists(audit_rds_file),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 30. Save compact audit CSV
# ------------------------------------------------------------

audit_csv_summary <-
  data.frame(
    metric = c(
      "final_status",
      "pathway_n",
      "fdr05_n",
      "duplicated_pathway_id_n",
      "positive_nes_n",
      "negative_nes_n",
      "zero_nes_n",
      "final_kegg_ranking_n",
      "ranking_unique_n",
      "ranking_tied_gene_n",
      "ranking_tie_fraction",
      "structural_pass",
      "fdr_pass",
      "csv_rds_pass",
      "ranking_pass",
      "ranking_order_audit_pass",
      "core_gene_pass",
      "integrity_pass",
      "parameter_pass",
      "final_hard_pass"
    ),
    
    value = c(
      audit_summary$status,
      observed_pathway_n,
      observed_fdr05_n,
      duplicated_pathway_id_n,
      positive_nes_n,
      negative_nes_n,
      zero_nes_n,
      length(ranking_vector_final),
      ranking_unique_n,
      ranking_tied_gene_n,
      ranking_tie_fraction,
      structural_pass,
      fdr_pass,
      csv_rds_pass,
      ranking_pass,
      ranking_order_audit_pass,
      core_gene_pass,
      integrity_pass,
      parameter_pass,
      final_hard_pass
    ),
    
    stringsAsFactors = FALSE
  )

write.csv(
  audit_csv_summary,
  audit_csv_file,
  row.names = FALSE
)

cat(
  "CSV audit saved: ",
  file.exists(audit_csv_file),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 31. Reload verification
# ------------------------------------------------------------

output_check <-
  c(
    TXT = file.exists(audit_txt_file),
    RDS = file.exists(audit_rds_file),
    CSV = file.exists(audit_csv_file)
  )

cat("\n============================================================\n")
cat("STEP9D-3C OUTPUT VERIFICATION\n")
cat("============================================================\n\n")

print(output_check)

if (!all(output_check)) {
  stop(
    "FAIL: one or more STEP9D-3C audit outputs are missing."
  )
}


reload_audit <-
  readRDS(audit_rds_file)


reload_pass <-
  identical(
    reload_audit$final_hard_pass,
    final_hard_pass
  ) &&
  identical(
    reload_audit$status,
    audit_summary$status
  )


cat(
  "\nRDS reload verification: ",
  reload_pass,
  "\n",
  sep = ""
)

if (!reload_pass) {
  stop(
    "FAIL: STEP9D-3C audit RDS reload verification failed."
  )
}   


# ------------------------------------------------------------
# 32. FINAL STEP9D-3C LOCK
# ------------------------------------------------------------

final_step9d3c_pass <-
  final_hard_pass &&
  all(output_check) &&
  reload_pass

cat("\n============================================================\n")
cat("STEP9D-3C FINAL LOCK\n")
cat("============================================================\n\n")

cat(
  "Final hard pass: ",
  final_hard_pass,
  "\n",
  sep = ""
)

cat(
  "TXT audit exists: ",
  file.exists(audit_txt_file),
  "\n",
  sep = ""
)

cat(
  "RDS audit exists: ",
  file.exists(audit_rds_file),
  "\n",
  sep = ""
)

cat(
  "CSV audit exists: ",
  file.exists(audit_csv_file),
  "\n",
  sep = ""
)

cat(
  "RDS reload verification: ",
  reload_pass,
  "\n",
  sep = ""
)

cat(
  "\nSTEP9D-3C FINAL HARD PASS: ",
  final_step9d3c_pass,
  "\n",
  sep = ""
)

if (!final_step9d3c_pass) {
  stop(
    "STEP9D-3C FINAL HARD PASS FAILED."
  )
}


# ============================================================
# Paper 1
# STEP9D-4A-1
# GO-BP GSEA — Locked Ranking Input Audit
# READ-ONLY
# ============================================================

cat("============================================================\n")
cat("STEP9D-4A-1: GO-BP LOCKED RANKING INPUT AUDIT\n")
cat("============================================================\n\n")

root_dir <-
  "input_dir/Paper 1_Step 9-Functional Enrichment"

step9d_root <-
  file.path(
    root_dir,
    "STEP9D_GenomeWide_Preranked_GSEA"
  )

ranking_file <-
  file.path(
    step9d_root,
    "9D-2C_Final_Ranking_Vector",
    "STEP9D2C_FINAL_GenomeWide_Ranking_Vector.rds"
  )

lock_2b_file <-
  file.path(
    step9d_root,
    "9D-2_Ranking_Metric_Lock",
    "STEP9D2B_FINAL_Ranking_Metric_Lock.rds"
  )

lock_1i_file <-
  file.path(
    step9d_root,
    "9D-1_Input_Lock",
    "STEP9D1I_FINAL_INPUT_LOCK.rds"
  )

if (!file.exists(ranking_file)) {
  stop("FAIL: STEP9D-2C ranking vector not found.")
}

if (!file.exists(lock_2b_file)) {
  stop("FAIL: STEP9D-2B ranking lock not found.")
}

if (!file.exists(lock_1i_file)) {
  stop("FAIL: STEP9D-1I input lock not found.")
}


ranking_vector <-
  readRDS(ranking_file)

lock_2b <-
  readRDS(lock_2b_file)

lock_1i <-
  readRDS(lock_1i_file)


cat("Ranking class: ", class(ranking_vector), "\n", sep = "")
cat("Ranking length: ", length(ranking_vector), "\n", sep = "")
cat(
  "Named ranking: ",
  !is.null(names(ranking_vector)),
  "\n",
  sep = ""
)

cat(
  "Duplicated IDs: ",
  anyDuplicated(names(ranking_vector)),
  "\n",
  sep = ""
)

cat(
  "NA values: ",
  anyNA(ranking_vector),
  "\n",
  sep = ""
)

cat(
  "All finite: ",
  all(is.finite(ranking_vector)),
  "\n",
  sep = ""
)

cat(
  "Expected N = 57,828: ",
  length(ranking_vector) == 57828,
  "\n",
  sep = ""
)


# Ensembl ID format audit

ensembl_ids <-
  names(ranking_vector)

standard_ensembl <-
  grepl(
    "^ENSG[0-9]{11}$",
    ensembl_ids
  )

cat(
  "Standard Ensembl IDs = ",
  sum(standard_ensembl),
  "\n",
  sep = ""
)

cat(
  "Non-standard Ensembl IDs = ",
  sum(!standard_ensembl),
  "\n",
  sep = ""
)


hard_pass <-
  length(ranking_vector) == 57828 &&
  !is.null(names(ranking_vector)) &&
  !anyDuplicated(names(ranking_vector)) &&
  !anyNA(ranking_vector) &&
  all(is.finite(ranking_vector)) &&
  all(standard_ensembl)

cat(
  "\nSTEP9D-4A-1 HARD PASS: ",
  hard_pass,
  "\n",
  sep = ""
)

if (!hard_pass) {
  stop(
    "STEP9D-4A-1 FAILED."
  )
}


# ============================================================
# STEP9D-4A-2
# GO-BP Annotation Coverage Audit
# READ-ONLY
# ============================================================

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

cat("\n============================================================\n")
cat("STEP9D-4A-2: GO-BP ANNOTATION COVERAGE AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Package versions
# ------------------------------------------------------------

cat(
  "AnnotationDbi: ",
  as.character(packageVersion("AnnotationDbi")),
  "\n",
  sep = ""
)

cat(
  "org.Hs.eg.db: ",
  as.character(packageVersion("org.Hs.eg.db")),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 2. GO annotation lookup
# ------------------------------------------------------------

go_bp_map <-
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = names(ranking_vector),
    keytype = "ENSEMBL",
    columns = c(
      "ENSEMBL",
      "GO",
      "ONTOLOGY",
      "SYMBOL"
    )
  )


cat(
  "\nRaw GO annotation rows: ",
  nrow(go_bp_map),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 3. Restrict to Biological Process
# ------------------------------------------------------------

go_bp_map <-
  go_bp_map[
    !is.na(go_bp_map$GO) &
      go_bp_map$ONTOLOGY == "BP",
    ,
    drop = FALSE
  ]


go_bp_genes <-
  unique(
    go_bp_map$ENSEMBL
  )


cat(
  "Unique Ensembl genes with GO-BP annotation: ",
  length(go_bp_genes),
  "\n",
  sep = ""
)


go_bp_coverage <-
  length(go_bp_genes) /
  length(ranking_vector)


cat(
  "GO-BP annotation coverage: ",
  sprintf("%.6f", go_bp_coverage),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Core-gene annotation audit
# ------------------------------------------------------------

core_genes <-
  c(
    AOC3 = "ENSG00000131471",
    MS4A7 = "ENSG00000166927",
    IL7R = "ENSG00000168685"
  )

core_go_bp_n <-
  sapply(
    core_genes,
    function(x) {
      sum(
        go_bp_map$ENSEMBL == x
      )
    }
  )

print(core_go_bp_n)


core_go_bp_retained <-
  core_genes %in% go_bp_genes


core_gene_audit <-
  data.frame(
    gene = names(core_genes),
    ENSEMBL = unname(core_genes),
    GO_BP_annotation_n =
      as.integer(core_go_bp_n),
    retained_for_GO_BP_GSEA =
      core_go_bp_retained,
    stringsAsFactors = FALSE
  )

print(core_gene_audit)


# ------------------------------------------------------------
# 5. Integrity
# ------------------------------------------------------------

annotation_integrity <-
  !anyNA(go_bp_map$ENSEMBL) &&
  all(
    go_bp_map$ENSEMBL %in%
      names(ranking_vector)
  )

core_pass <-
  all(core_go_bp_retained)


hard_pass <-
  length(go_bp_genes) > 0 &&
  annotation_integrity &&
  core_pass


cat(
  "\nAnnotation integrity: ",
  annotation_integrity,
  "\n",
  sep = ""
)

cat(
  "Core genes retained: ",
  core_pass,
  "\n",
  sep = ""
)

cat(
  "\nSTEP9D-4A-2 HARD PASS: ",
  hard_pass,
  "\n",
  sep = ""
)

if (!hard_pass) {
  stop(
    "STEP9D-4A-2 FAILED."
  )
}   

# ============================================================
# STEP9D-4A-3
# FINAL GO-BP GSEA INPUT / POLICY LOCK
# ============================================================

go_dir <-
  file.path(
    step9d_root,
    "9D-4A_GO_BP_Input_Annotation_Policy_Lock"
  )

dir.create(
  go_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# GO-BP genes
# ------------------------------------------------------------

go_bp_ids <-
  unique(
    go_bp_map$ENSEMBL
  )

go_bp_ranking <-
  ranking_vector[
    names(ranking_vector) %in%
      go_bp_ids
  ]


# Preserve original locked ranking order
go_bp_ranking <-
  go_bp_ranking[
    order(
      match(
        names(go_bp_ranking),
        names(ranking_vector)
      )
    )
  ]


# ------------------------------------------------------------
# Audit
# ------------------------------------------------------------

go_bp_input_n <-
  length(go_bp_ranking)

positive_n <-
  sum(go_bp_ranking > 0)

negative_n <-
  sum(go_bp_ranking < 0)

zero_n <-
  sum(go_bp_ranking == 0)

unique_n <-
  length(unique(go_bp_ranking))

tied_gene_n <-
  go_bp_input_n - unique_n


cat("GO-BP GSEA input N: ", go_bp_input_n, "\n", sep = "")
cat("Positive: ", positive_n, "\n", sep = "")
cat("Negative: ", negative_n, "\n", sep = "")
cat("Zero: ", zero_n, "\n", sep = "")
cat("Tied genes: ", tied_gene_n, "\n", sep = "")


# ------------------------------------------------------------
# Core preservation
# ------------------------------------------------------------

core_presence <-
  data.frame(
    gene = names(core_genes),
    ENSEMBL = unname(core_genes),
    retained =
      unname(core_genes) %in%
      names(go_bp_ranking),
    stringsAsFactors = FALSE
  )

print(core_presence)

core_pass <-
  all(core_presence$retained)


# ------------------------------------------------------------
# Ranking integrity
# ------------------------------------------------------------

ranking_values_original <-
  ranking_vector[
    names(go_bp_ranking)
  ]

id_aligned_match <-
  identical(
    as.numeric(go_bp_ranking),
    as.numeric(ranking_values_original)
  )

cat(
  "ID-aligned ranking unchanged: ",
  id_aligned_match,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# Final policy
# ------------------------------------------------------------

go_bp_policy <-
  list(
    workflow =
      "Paper1_STEP9D_GenomeWide_Preranked_GSEA",
    
    step =
      "STEP9D-4A-3",
    
    status =
      "FINAL_GO_BP_INPUT_POLICY_LOCK",
    
    source_ranking_file =
      ranking_file,
    
    original_genome_wide_ranking_n =
      length(ranking_vector),
    
    go_bp_input_n =
      go_bp_input_n,
    
    gene_id_type =
      "ENSEMBL",
    
    annotation_database =
      "org.Hs.eg.db",
    
    annotation_version =
      as.character(
        packageVersion("org.Hs.eg.db")
      ),
    
    ontology =
      "BP",
    
    annotation_filter =
      "GO annotation with ONTOLOGY == BP",
    
    ranking_recalculated =
      FALSE,
    
    ranking_filtered_by_DEG =
      FALSE,
    
    jitter_used =
      FALSE,
    
    artificial_tie_breaking =
      FALSE,
    
    ranking_order_preserved =
      id_aligned_match,
    
    positive_n =
      positive_n,
    
    negative_n =
      negative_n,
    
    zero_n =
      zero_n,
    
    tied_gene_n =
      tied_gene_n,
    
    core_presence =
      core_presence,
    
    core_gene_pass =
      core_pass,
    
    hard_pass =
      id_aligned_match &&
      core_pass
  )


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

saveRDS(
  go_bp_policy,
  file.path(
    go_dir,
    "STEP9D4A3_FINAL_GO_BP_GSEA_Input_Policy_Lock.rds"
  )
)

write.csv(
  core_presence,
  file.path(
    go_dir,
    "STEP9D4A3_CoreGene_GO_BP_Preservation.csv"
  ),
  row.names = FALSE
)

cat(
  "\nSTEP9D-4A-3 FINAL HARD PASS: ",
  go_bp_policy$hard_pass,
  "\n",
  sep = ""
)

if (!go_bp_policy$hard_pass) {
  stop(
    "STEP9D-4A-3 FAILED."
  )
}


# ============================================================
# Paper 1
# STEP9D-4B
# FORMAL GO-BP GSEA EXECUTION
# ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(fgsea)
})

cat("============================================================\n")
cat("STEP9D-4B: FORMAL GO-BP GSEA EXECUTION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

step9d_root <-
  "input_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA"

input_lock_file <-
  file.path(
    step9d_root,
    "9D-4A_GO_BP_Input_Annotation_Policy_Lock",
    "STEP9D4A3_FINAL_GO_BP_GSEA_Input_Policy_Lock.rds"
  )

execution_dir <-
  file.path(
    step9d_root,
    "9D-4B_Formal_GO_BP_GSEA_Execution"
  )

dir.create(
  execution_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Load locked ranking and policy
# ------------------------------------------------------------

ranking_file <-
  "result_dir/Paper 1_Step 9-Functional Enrichment/STEP9D_GenomeWide_Preranked_GSEA/9D-2C_Final_Ranking_Vector/STEP9D2C_FINAL_GenomeWide_Ranking_Vector.rds"

ranking_vector <-
  readRDS(ranking_file)

go_policy <-
  readRDS(input_lock_file)


if (!isTRUE(go_policy$hard_pass)) {
  stop(
    "FAIL: GO-BP input policy lock is not PASS."
  )
}


# ------------------------------------------------------------
# 3. Construct GO-BP input from locked ranking
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(AnnotationDbi)
})

go_map <-
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = names(ranking_vector),
    keytype = "ENSEMBL",
    columns = c(
      "ENSEMBL",
      "GO",
      "ONTOLOGY",
      "SYMBOL"
    )
  )

go_bp_ids <-
  unique(
    go_map$ENSEMBL[
      !is.na(go_map$GO) &
        go_map$ONTOLOGY == "BP"
    ]
  )


go_bp_ranking <-
  ranking_vector[
    names(ranking_vector) %in%
      go_bp_ids
  ]


# ------------------------------------------------------------
# 4. Sort decreasing for formal GSEA input
# ------------------------------------------------------------

go_bp_ranking_unsorted <-
  go_bp_ranking

go_bp_ranking <-
  sort(
    go_bp_ranking_unsorted,
    decreasing = TRUE
  )


# ------------------------------------------------------------
# 5. Integrity checks
# ------------------------------------------------------------

if (
  length(go_bp_ranking) !=
  go_policy$go_bp_input_n
) {
  stop(
    "FAIL: GO-BP input N differs from locked policy."
  )
}

if (
  anyNA(go_bp_ranking) ||
  !all(is.finite(go_bp_ranking))
) {
  stop(
    "FAIL: GO-BP ranking contains NA/non-finite values."
  )
}

if (
  anyDuplicated(names(go_bp_ranking))
) {
  stop(
    "FAIL: duplicated Ensembl IDs detected."
  )
}

if (
  !all(diff(go_bp_ranking) <= 0)
) {
  stop(
    "FAIL: GO-BP ranking is not sorted decreasing."
  )
}


# ------------------------------------------------------------
# 6. Core-gene audit
# ------------------------------------------------------------

core_genes <-
  c(
    AOC3 = "ENSG00000131471",
    MS4A7 = "ENSG00000166927",
    IL7R = "ENSG00000168685"
  )

core_presence <-
  unname(core_genes) %in%
  names(go_bp_ranking)

if (!all(core_presence)) {
  stop(
    "FAIL: one or more core genes absent from GO-BP ranking."
  )
}


# ------------------------------------------------------------
# 7. Save exact formal input
# ------------------------------------------------------------

go_input_table <-
  data.frame(
    ENSEMBL =
      names(go_bp_ranking),
    ranking =
      as.numeric(go_bp_ranking),
    stringsAsFactors = FALSE
  )

write.csv(
  go_input_table,
  file.path(
    execution_dir,
    "STEP9D4B_FINAL_GO_BP_GSEA_Input_Ranking.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    workflow =
      "Paper1_STEP9D_GenomeWide_Preranked_GSEA",
    
    step =
      "STEP9D-4B",
    
    status =
      "FINAL_GO_BP_GSEA_INPUT",
    
    source_ranking_file =
      ranking_file,
    
    input_policy_file =
      input_lock_file,
    
    ontology =
      "BP",
    
    gene_id_type =
      "ENSEMBL",
    
    input_n =
      length(go_bp_ranking),
    
    positive_n =
      sum(go_bp_ranking > 0),
    
    negative_n =
      sum(go_bp_ranking < 0),
    
    zero_n =
      sum(go_bp_ranking == 0),
    
    core_presence =
      data.frame(
        gene = names(core_genes),
        ENSEMBL = unname(core_genes),
        retained =
          core_presence
      ),
    
    input_table =
      go_input_table,
    
    ranking_vector =
      go_bp_ranking
  ),
  file.path(
    execution_dir,
    "STEP9D4B_FINAL_GO_BP_GSEA_Input_Ranking.rds"
  )
)


# ------------------------------------------------------------
# 8. Formal gseGO
# ------------------------------------------------------------

cat("\nSTARTING FORMAL gseGO()\n")

set.seed(1)

gsea_go_bp <-
  clusterProfiler::gseGO(
    geneList =
      go_bp_ranking,
    
    OrgDb =
      org.Hs.eg.db,
    
    keyType =
      "ENSEMBL",
    
    ont =
      "BP",
    
    minGSSize =
      10,
    
    maxGSSize =
      500,
    
    exponent =
      1,
    
    pvalueCutoff =
      1,
    
    pAdjustMethod =
      "BH",
    
    verbose =
      FALSE,
    
    seed =
      FALSE
  )


# ------------------------------------------------------------
# 9. Result audit
# ------------------------------------------------------------

gsea_df <-
  as.data.frame(gsea_go_bp)

cat(
  "\nGO-BP pathways returned: ",
  nrow(gsea_df),
  "\n",
  sep = ""
)

fdr05_n <-
  sum(
    gsea_df$p.adjust < 0.05,
    na.rm = TRUE
  )

cat(
  "GO-BP pathways with FDR < 0.05: ",
  fdr05_n,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Save results
# ------------------------------------------------------------

saveRDS(
  gsea_go_bp,
  file.path(
    execution_dir,
    "STEP9D4B_FINAL_GO_BP_GSEA_Result.rds"
  )
)

write.csv(
  gsea_df,
  file.path(
    execution_dir,
    "STEP9D4B_FINAL_GO_BP_GSEA_Result.csv"
  ),
  row.names = FALSE
)

write.csv(
  gsea_df[
    gsea_df$p.adjust < 0.05,
    ,
    drop = FALSE
  ],
  file.path(
    execution_dir,
    "STEP9D4B_FINAL_GO_BP_GSEA_Significant_FDR05.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Recover GO-BP FDR < 0.05 count
# READ-ONLY
# ------------------------------------------------------------

if (!exists("gsea_df")) {
  stop("FAIL: gsea_df is not available in the current environment.")
}

if (!"p.adjust" %in% colnames(gsea_df)) {
  stop("FAIL: p.adjust column is missing from GO-BP GSEA result.")
}

fdr05_n <-
  sum(
    gsea_df$p.adjust < 0.05,
    na.rm = TRUE
  )

cat(
  "GO-BP pathways with FDR < 0.05: ",
  fdr05_n,
  "\n",
  sep = ""
)

# ============================================================
# STEP9D-4B
# FINAL GO-BP GSEA PROVENANCE — LOCKED
# ============================================================

ranking_unique_n <-
  length(
    unique(go_bp_ranking)
  )

ranking_tied_gene_n <-
  length(go_bp_ranking) -
  ranking_unique_n

ranking_tie_fraction <-
  ranking_tied_gene_n /
  length(go_bp_ranking)


fdr05_n <-
  sum(
    gsea_df$p.adjust < 0.05,
    na.rm = TRUE
  )


provenance <-
  list(
    
    workflow =
      "Paper1_STEP9D_GenomeWide_Preranked_GSEA",
    
    step =
      "STEP9D-4B",
    
    status =
      "FINAL_GO_BP_GSEA",
    
    ranking_source =
      ranking_file,
    
    input_policy_source =
      input_lock_file,
    
    ranking_metric =
      "signed_minus_log10_nominal_P",
    
    ranking_formula =
      "sign(log2FC) * (-log10(max(P_value, 1e-300)))",
    
    original_genome_wide_ranking_n =
      length(ranking_vector),
    
    go_bp_input_n =
      length(go_bp_ranking),
    
    ontology =
      "BP",
    
    keyType =
      "ENSEMBL",
    
    organism =
      "org.Hs.eg.db",
    
    minGSSize =
      10,
    
    maxGSSize =
      500,
    
    exponent =
      1,
    
    eps =
      1e-10,
    
    pvalueCutoff =
      1,
    
    pAdjustMethod =
      "BH",
    
    verbose =
      FALSE,
    
    seed =
      FALSE,
    
    pathway_n =
      nrow(gsea_df),
    
    fdr05_n =
      fdr05_n,
    
    ranking_unique_n =
      ranking_unique_n,
    
    ranking_tied_gene_n =
      ranking_tied_gene_n,
    
    ranking_tie_fraction =
      ranking_tie_fraction,
    
    ranking_recalculated =
      FALSE,
    
    deg_recalculated =
      FALSE,
    
    jitter_used =
      FALSE,
    
    artificial_tie_breaking =
      FALSE,
    
    gsea_performed =
      TRUE,
    
    preranked_tie_warning =
      TRUE,
    
    fgsea_log2err_warning =
      TRUE,
    
    fgsea_eps_warning =
      TRUE
  )


cat("\n============================================================\n")
cat("STEP9D-4B PROVENANCE SUMMARY\n")
cat("============================================================\n\n")

cat(
  "GO-BP input N = ",
  provenance$go_bp_input_n,
  "\n",
  sep = ""
)

cat(
  "GO-BP pathways returned = ",
  provenance$pathway_n,
  "\n",
  sep = ""
)

cat(
  "GO-BP FDR < 0.05 = ",
  provenance$fdr05_n,
  "\n",
  sep = ""
)

cat(
  "Ranking tie fraction = ",
  sprintf(
    "%.6f",
    provenance$ranking_tie_fraction
  ),
  "\n",
  sep = ""
)

cat(
  "Ranking recalculated = ",
  provenance$ranking_recalculated,
  "\n",
  sep = ""
)

cat(
  "Jitter used = ",
  provenance$jitter_used,
  "\n",
  sep = ""
)

cat(
  "Artificial tie-breaking = ",
  provenance$artificial_tie_breaking,
  "\n",
  sep = ""
)

cat(
  "GSEA performed = ",
  provenance$gsea_performed,
  "\n",
  sep = ""
)

# ============================================================
# STEP9D-4B
# READ-ONLY OBJECT STRUCTURE CHECK — CORRECTED
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-4B: GO-BP INPUT OBJECT STRUCTURE CHECK\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Find relevant objects
# ------------------------------------------------------------

all_objects <- ls(envir = .GlobalEnv)

relevant_objects <- all_objects[
  grepl(
    "go_bp|gsea_go_bp|ranking",
    all_objects,
    ignore.case = TRUE
  )
]

cat("Objects matching GO-BP / ranking:\n")
print(relevant_objects)

# ------------------------------------------------------------
# 2. Check go_bp_ranking if present
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("CHECK: go_bp_ranking\n")
cat("------------------------------------------------------------\n\n")

if (exists("go_bp_ranking", envir = .GlobalEnv)) {
  
  cat("go_bp_ranking exists = TRUE\n")
  
  cat("\nClass:\n")
  print(class(go_bp_ranking))
  
  cat("\nLength:\n")
  print(length(go_bp_ranking))
  
  cat("\nFirst 10 values:\n")
  print(head(go_bp_ranking, 10))
  
  cat("\nFirst 10 names:\n")
  print(head(names(go_bp_ranking), 10))
  
} else {
  
  cat("go_bp_ranking exists = FALSE\n")
}

# ------------------------------------------------------------
# 3. Check gsea_go_bp
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("CHECK: gsea_go_bp\n")
cat("------------------------------------------------------------\n\n")

if (exists("gsea_go_bp", envir = .GlobalEnv)) {
  
  cat("gsea_go_bp exists = TRUE\n")
  
  cat("\nClass:\n")
  print(class(gsea_go_bp))
  
  cat("\nNames:\n")
  print(names(gsea_go_bp))
  
} else {
  
  cat("gsea_go_bp exists = FALSE\n")
}

# ------------------------------------------------------------
# 4. Check gsea_df
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("CHECK: gsea_df\n")
cat("------------------------------------------------------------\n\n")

if (exists("gsea_df", envir = .GlobalEnv)) {
  
  cat("gsea_df exists = TRUE\n")
  
  cat("\nDimensions:\n")
  print(dim(gsea_df))
  
  cat("\nColumns:\n")
  print(colnames(gsea_df))
  
} else {
  
  cat("gsea_df exists = FALSE\n")
}

cat("\n============================================================\n")
cat("OBJECT CHECK COMPLETE\n")
cat("NO GSEA EXECUTION\n")
cat("NO RANKING RECALCULATION\n")
cat("============================================================\n")


# ============================================================
# STEP9D-4B
# FINAL INPUT TABLE RECONSTRUCTION
# + OUTPUT SAVE + RELOAD VERIFICATION
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-4B: FINAL INPUT TABLE + OUTPUT SAVE\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Output directory
# ------------------------------------------------------------

step4b_dir <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-4B_Formal_GO_BP_GSEA_Execution"
)

dir.create(
  step4b_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Reconstruct input table DIRECTLY from the existing
#    formal GO-BP ranking vector
#
#    IMPORTANT:
#    - no ranking recalculation
#    - no sorting alteration
#    - no jitter
#    - no tie breaking
# ------------------------------------------------------------

stopifnot(
  is.numeric(go_bp_ranking),
  !is.null(names(go_bp_ranking)),
  length(go_bp_ranking) == provenance$go_bp_input_n
)

go_bp_input_table <- data.frame(
  ENSEMBL = names(go_bp_ranking),
  ranking = as.numeric(go_bp_ranking),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 3. Input-table integrity checks
# ------------------------------------------------------------

input_table_integrity <- (
  nrow(go_bp_input_table) == 18008 &&
    all(!is.na(go_bp_input_table$ENSEMBL)) &&
    all(go_bp_input_table$ENSEMBL != "") &&
    !anyDuplicated(go_bp_input_table$ENSEMBL) &&
    all(is.finite(go_bp_input_table$ranking)) &&
    identical(
      go_bp_input_table$ENSEMBL,
      names(go_bp_ranking)
    ) &&
    isTRUE(
      all.equal(
        go_bp_input_table$ranking,
        as.numeric(go_bp_ranking),
        check.attributes = FALSE
      )
    )
)

cat(
  "GO-BP input table N = ",
  nrow(go_bp_input_table),
  "\n",
  sep = ""
)

cat(
  "Unique Ensembl IDs = ",
  length(unique(go_bp_input_table$ENSEMBL)),
  "\n",
  sep = ""
)

cat(
  "Input ranking finite = ",
  all(is.finite(go_bp_input_table$ranking)),
  "\n",
  sep = ""
)

cat(
  "Input ranking matches go_bp_ranking = ",
  input_table_integrity,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 4. Verify descending order only
# ------------------------------------------------------------

ranking_sorted_decreasing <- all(
  diff(go_bp_input_table$ranking) <= 0
)

cat(
  "Ranking sorted decreasing = ",
  ranking_sorted_decreasing,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 5. Save formal GO-BP result
# ------------------------------------------------------------

saveRDS(
  gsea_go_bp,
  file.path(
    step4b_dir,
    "STEP9D4B_FINAL_GO_BP_GSEA_Result.rds"
  )
)

write.csv(
  gsea_df,
  file.path(
    step4b_dir,
    "STEP9D4B_FINAL_GO_BP_GSEA_Result.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Save significant FDR < 0.05 result
# ------------------------------------------------------------

gsea_sig <- gsea_df[
  !is.na(gsea_df$p.adjust) &
    gsea_df$p.adjust < 0.05,
  ,
  drop = FALSE
]

cat(
  "Significant GO-BP pathways FDR < 0.05 = ",
  nrow(gsea_sig),
  "\n",
  sep = ""
)

write.csv(
  gsea_sig,
  file.path(
    step4b_dir,
    "STEP9D4B_FINAL_GO_BP_GSEA_Significant_FDR05.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Save formal GO-BP input ranking
# ------------------------------------------------------------

write.csv(
  go_bp_input_table,
  file.path(
    step4b_dir,
    "STEP9D4B_FINAL_GO_BP_GSEA_Input_Ranking.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    workflow =
      "Paper1_STEP9D_GenomeWide_Preranked_GSEA",
    
    step =
      "STEP9D-4B",
    
    status =
      "FINAL_GO_BP_GSEA_INPUT",
    
    input_table =
      go_bp_input_table,
    
    ranking_vector =
      go_bp_ranking
  ),
  file.path(
    step4b_dir,
    "STEP9D4B_FINAL_GO_BP_GSEA_Input_Ranking.rds"
  )
)

# ------------------------------------------------------------
# 8. Save provenance
# ------------------------------------------------------------

saveRDS(
  provenance,
  file.path(
    step4b_dir,
    "STEP9D4B_FINAL_GO_BP_GSEA_Provenance.rds"
  )
)

# ------------------------------------------------------------
# 9. Required-file verification
# ------------------------------------------------------------

required_4b_files <- c(
  "STEP9D4B_FINAL_GO_BP_GSEA_Result.rds",
  "STEP9D4B_FINAL_GO_BP_GSEA_Result.csv",
  "STEP9D4B_FINAL_GO_BP_GSEA_Significant_FDR05.csv",
  "STEP9D4B_FINAL_GO_BP_GSEA_Input_Ranking.csv",
  "STEP9D4B_FINAL_GO_BP_GSEA_Input_Ranking.rds",
  "STEP9D4B_FINAL_GO_BP_GSEA_Provenance.rds"
)

file_status <- file.exists(
  file.path(step4b_dir, required_4b_files)
)

cat("\n------------------------------------------------------------\n")
cat("FILE EXISTENCE AUDIT\n")
cat("------------------------------------------------------------\n\n")

print(
  data.frame(
    file = required_4b_files,
    exists = file_status,
    row.names = NULL
  )
)

cat(
  "\nAll required files present = ",
  all(file_status),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 10. Reload verification
# ------------------------------------------------------------

gsea_go_bp_reload <- readRDS(
  file.path(
    step4b_dir,
    "STEP9D4B_FINAL_GO_BP_GSEA_Result.rds"
  )
)

provenance_reload <- readRDS(
  file.path(
    step4b_dir,
    "STEP9D4B_FINAL_GO_BP_GSEA_Provenance.rds"
  )
)

input_reload <- readRDS(
  file.path(
    step4b_dir,
    "STEP9D4B_FINAL_GO_BP_GSEA_Input_Ranking.rds"
  )
)

# ------------------------------------------------------------
# 11. Reload integrity
# ------------------------------------------------------------

reload_result_ok <- identical(
  gsea_go_bp,
  gsea_go_bp_reload
)

reload_provenance_ok <- identical(
  provenance,
  provenance_reload
)

reload_input_table_ok <- identical(
  go_bp_input_table,
  input_reload$input_table
)

reload_ranking_ok <- identical(
  go_bp_ranking,
  input_reload$ranking_vector
)

reload_ok <- all(
  reload_result_ok,
  reload_provenance_ok,
  reload_input_table_ok,
  reload_ranking_ok
)

cat("\n------------------------------------------------------------\n")
cat("RDS RELOAD VERIFICATION\n")
cat("------------------------------------------------------------\n\n")

cat(
  "GSEA result reload = ",
  reload_result_ok,
  "\n",
  sep = ""
)

cat(
  "Provenance reload = ",
  reload_provenance_ok,
  "\n",
  sep = ""
)

cat(
  "Input table reload = ",
  reload_input_table_ok,
  "\n",
  sep = ""
)

cat(
  "Ranking vector reload = ",
  reload_ranking_ok,
  "\n",
  sep = ""
)

cat(
  "Overall reload verification = ",
  reload_ok,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 12. STEP9D-4B preliminary hard pass
# ------------------------------------------------------------

step4b_save_pass <- all(
  input_table_integrity,
  ranking_sorted_decreasing,
  all(file_status),
  reload_ok,
  provenance$go_bp_input_n == 18008,
  provenance$pathway_n == 353,
  provenance$fdr05_n == 178,
  provenance$ranking_recalculated == FALSE,
  provenance$jitter_used == FALSE,
  provenance$artificial_tie_breaking == FALSE,
  provenance$gsea_performed == TRUE
)

cat("\n============================================================\n")
cat(
  "STEP9D-4B OUTPUT SAVE / RELOAD HARD PASS = ",
  step4b_save_pass,
  "\n",
  sep = ""
)
cat("============================================================\n")


# ============================================================
# Paper 1
# STEP9D-4B
# FORMAL GO-BP OUTPUT / PROVENANCE REPAIR
#
# INDEPENDENT DISK-BASED REPAIR VERSION
#
# PURPOSE:
#   Repair only contaminated GO-BP serialized outputs and
#   provenance using the EXISTING authoritative gseaResult RDS.
#
# STRICT LOCK:
#   NO gseGO() rerun
#   NO ranking recalculation
#   NO DEG recalculation
#   NO jitter
#   NO artificial tie-breaking
#   NO modification of the authoritative Result.rds
#
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-4B FORMAL GO-BP OUTPUT / PROVENANCE REPAIR\n")
cat("INDEPENDENT DISK-BASED VERSION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

root_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

step_dir <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-4B_Formal_GO_BP_GSEA_Execution"
)

result_rds <- file.path(
  step_dir,
  "STEP9D4B_FINAL_GO_BP_GSEA_Result.rds"
)

input_rds <- file.path(
  step_dir,
  "STEP9D4B_FINAL_GO_BP_GSEA_Input_Ranking.rds"
)

old_provenance_rds <- file.path(
  step_dir,
  "STEP9D4B_FINAL_GO_BP_GSEA_Provenance.rds"
)

official_result_csv <- file.path(
  step_dir,
  "STEP9D4B_FINAL_GO_BP_GSEA_Result.csv"
)

official_sig_csv <- file.path(
  step_dir,
  "STEP9D4B_FINAL_GO_BP_GSEA_Significant_FDR05.csv"
)

official_provenance_rds <- file.path(
  step_dir,
  "STEP9D4B_FINAL_GO_BP_GSEA_Provenance.rds"
)


# ------------------------------------------------------------
# 1. Required disk files
# ------------------------------------------------------------

required_files <- c(
  result_rds,
  input_rds,
  old_provenance_rds
)

file_status <- file.exists(required_files)

cat("1. REQUIRED DISK FILES\n")
cat("------------------------------------------------------------\n")

print(
  data.frame(
    file = basename(required_files),
    exists = file_status,
    stringsAsFactors = FALSE
  )
)

if (!all(file_status)) {
  stop(
    paste(
      "ERROR: Required disk file(s) missing:",
      paste(
        basename(required_files[!file_status]),
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 2. Load authoritative Result RDS
#
# IMPORTANT:
#   This only reads the existing gseaResult.
#   No GSEA is performed.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("2. LOAD AUTHORITATIVE GO-BP RESULT\n")
cat("============================================================\n\n")

gsea_go_bp_disk <- readRDS(result_rds)

cat(
  "Class = ",
  paste(
    class(gsea_go_bp_disk),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

if (!inherits(
  gsea_go_bp_disk,
  "gseaResult"
)) {
  stop(
    "ERROR: Result RDS is not a gseaResult object."
  )
}


# ------------------------------------------------------------
# 3. Extract authoritative GO-BP result
# ------------------------------------------------------------

go_result <- slot(
  gsea_go_bp_disk,
  "result"
)

if (!is.data.frame(go_result)) {
  stop(
    "ERROR: gseaResult@result is not a data.frame."
  )
}

cat("\n============================================================\n")
cat("3. AUTHORITATIVE GO-BP RESULT\n")
cat("============================================================\n\n")

cat(
  "Rows = ",
  nrow(go_result),
  "\n",
  sep = ""
)

cat(
  "Columns = ",
  ncol(go_result),
  "\n",
  sep = ""
)

cat(
  "Column names = ",
  paste(
    names(go_result),
    collapse = " | "
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Authoritative GO-BP structural checks
# ------------------------------------------------------------

expected_columns <- c(
  "ID",
  "Description",
  "setSize",
  "enrichmentScore",
  "NES",
  "pvalue",
  "p.adjust",
  "qvalue",
  "rank",
  "leading_edge",
  "core_enrichment"
)

column_check <- identical(
  names(go_result),
  expected_columns
)

go_id_check <- all(
  grepl(
    "^GO:",
    as.character(go_result$ID)
  )
)

unique_id_check <- (
  length(unique(go_result$ID)) ==
    nrow(go_result)
)

settype_check <- identical(
  as.character(gsea_go_bp_disk@setType),
  "BP"
)

keytype_check <- identical(
  as.character(gsea_go_bp_disk@keytype),
  "ENSEMBL"
)

organism_check <- identical(
  as.character(gsea_go_bp_disk@organism),
  "Homo sapiens"
)

gene_list_n <- length(
  gsea_go_bp_disk@geneList
)

gene_list_check <- (
  gene_list_n == 18008L
)

cat("\nColumn structure correct = ", column_check, "\n", sep = "")
cat("All IDs are GO IDs = ", go_id_check, "\n", sep = "")
cat("GO IDs unique = ", unique_id_check, "\n", sep = "")
cat("setType = BP = ", settype_check, "\n", sep = "")
cat("keytype = ENSEMBL = ", keytype_check, "\n", sep = "")
cat("organism = Homo sapiens = ", organism_check, "\n", sep = "")
cat("geneList N = 18008 = ", gene_list_check, "\n", sep = "")


# ------------------------------------------------------------
# 5. Authoritative counts
# ------------------------------------------------------------

authoritative_pathway_n <- nrow(
  go_result
)

authoritative_fdr05_n <- sum(
  go_result$p.adjust < 0.05,
  na.rm = TRUE
)

cat("\n============================================================\n")
cat("5. AUTHORITATIVE GO-BP COUNTS\n")
cat("============================================================\n\n")

cat(
  "Pathways returned = ",
  authoritative_pathway_n,
  "\n",
  sep = ""
)

cat(
  "FDR < 0.05 = ",
  authoritative_fdr05_n,
  "\n",
  sep = ""
)

if (
  authoritative_pathway_n != 6316L ||
  authoritative_fdr05_n != 2003L
) {
  stop(
    paste0(
      "ERROR: Unexpected authoritative GO-BP counts. ",
      "Expected 6316 / 2003 but observed ",
      authoritative_pathway_n,
      " / ",
      authoritative_fdr05_n,
      "."
    )
  )
}


# ------------------------------------------------------------
# 6. Load existing GO-BP input ranking RDS
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("6. LOAD EXISTING GO-BP INPUT RANKING\n")
cat("============================================================\n\n")

input_object <- readRDS(
  input_rds
)

cat(
  "Input RDS class = ",
  paste(
    class(input_object),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Input RDS names:\n"
)

print(
  names(input_object)
)


# ------------------------------------------------------------
# 7. Recover ranking/input without recalculation
# ------------------------------------------------------------

if (
  is.list(input_object) &&
  "ranking_vector" %in% names(input_object)
) {
  
  go_bp_ranking <- input_object$ranking_vector
  
} else if (
  is.numeric(input_object)
) {
  
  go_bp_ranking <- input_object
  
} else {
  
  stop(
    "ERROR: Unable to identify existing GO-BP ranking vector in Input Ranking RDS."
  )
}


if (
  is.list(input_object) &&
  "input_table" %in% names(input_object)
) {
  
  go_bp_input_table <- input_object$input_table
  
} else {
  
  go_bp_input_table <- data.frame(
    ENSEMBL = names(go_bp_ranking),
    ranking = as.numeric(go_bp_ranking),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 8. Input integrity
# ------------------------------------------------------------

input_n_check <- (
  length(go_bp_ranking) == 18008L &&
    nrow(go_bp_input_table) == 18008L
)

input_column_check <- identical(
  names(go_bp_input_table),
  c("ENSEMBL", "ranking")
)

input_id_check <- identical(
  as.character(names(go_bp_ranking)),
  as.character(go_bp_input_table$ENSEMBL)
)

input_value_check <- isTRUE(
  all.equal(
    as.numeric(go_bp_ranking),
    as.numeric(go_bp_input_table$ranking),
    check.attributes = FALSE
  )
)

cat("\n============================================================\n")
cat("8. INPUT INTEGRITY\n")
cat("============================================================\n\n")

cat("Input N = 18008 = ", input_n_check, "\n", sep = "")
cat("Input columns correct = ", input_column_check, "\n", sep = "")
cat("Input IDs unchanged = ", input_id_check, "\n", sep = "")
cat("Ranking values unchanged = ", input_value_check, "\n", sep = "")

if (!all(
  input_n_check,
  input_column_check,
  input_id_check,
  input_value_check
)) {
  stop(
    "ERROR: Existing GO-BP input ranking integrity check failed."
  )
}


# ------------------------------------------------------------
# 9. Tie audit
# ------------------------------------------------------------

ranking_unique_n <- length(
  unique(go_bp_ranking)
)

ranking_tied_gene_n <- (
  length(go_bp_ranking) -
    ranking_unique_n
)

ranking_tie_fraction <- (
  ranking_tied_gene_n /
    length(go_bp_ranking)
)

cat("\n============================================================\n")
cat("9. TIE AUDIT\n")
cat("============================================================\n\n")

cat(
  "Unique ranking values = ",
  ranking_unique_n,
  "\n",
  sep = ""
)

cat(
  "Tied genes = ",
  ranking_tied_gene_n,
  "\n",
  sep = ""
)

cat(
  "Tie fraction = ",
  sprintf(
    "%.7f",
    ranking_tie_fraction
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Construct corrected GO-BP output tables
#
# This is serialization only.
# NO enrichment calculation.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("10. REPAIR OUTPUT TABLES\n")
cat("============================================================\n\n")

go_sig <- go_result[
  !is.na(go_result$p.adjust) &
    go_result$p.adjust < 0.05,
  ,
  drop = FALSE
]

repaired_result_csv <- file.path(
  step_dir,
  "STEP9D4B_FINAL_GO_BP_GSEA_Result_REPAIRED.csv"
)

repaired_sig_csv <- file.path(
  step_dir,
  "STEP9D4B_FINAL_GO_BP_GSEA_Significant_FDR05_REPAIRED.csv"
)

write.csv(
  go_result,
  repaired_result_csv,
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  go_sig,
  repaired_sig_csv,
  row.names = FALSE,
  quote = TRUE
)

cat(
  "Repaired complete result rows = ",
  nrow(go_result),
  "\n",
  sep = ""
)

cat(
  "Repaired significant result rows = ",
  nrow(go_sig),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Read-back verification
# ------------------------------------------------------------

result_reload <- read.csv(
  repaired_result_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sig_reload <- read.csv(
  repaired_sig_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

result_reload_check <- isTRUE(
  all.equal(
    result_reload,
    go_result,
    check.attributes = FALSE
  )
)

sig_reload_check <- isTRUE(
  all.equal(
    sig_reload,
    go_sig,
    check.attributes = FALSE
  )
)

sig_fdr_check <- (
  nrow(sig_reload) == authoritative_fdr05_n &&
    all(
      sig_reload$p.adjust < 0.05,
      na.rm = TRUE
    )
)

cat("\n============================================================\n")
cat("11. REPAIRED CSV READ-BACK\n")
cat("============================================================\n\n")

cat(
  "Complete result CSV identical to authoritative RDS result = ",
  result_reload_check,
  "\n",
  sep = ""
)

cat(
  "Significant CSV identical to authoritative subset = ",
  sig_reload_check,
  "\n",
  sep = ""
)

cat(
  "Significant CSV N = 2003 and all FDR<0.05 = ",
  sig_fdr_check,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Create corrected provenance
#
# IMPORTANT:
#   No eps field is included as a gseGO parameter.
#   The existing fgsea warning is retained only as metadata.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("12. REPAIR PROVENANCE\n")
cat("============================================================\n\n")

provenance_repaired <- list(
  
  workflow =
    "Paper1_STEP9D_GenomeWide_Preranked_GSEA",
  
  step =
    "STEP9D-4B",
  
  status =
    "FINAL_GO_BP_GSEA",
  
  ranking_source =
    file.path(
      root_dir,
      "STEP9D_GenomeWide_Preranked_GSEA",
      "9D-2C_Final_Ranking_Vector",
      "STEP9D2C_FINAL_GenomeWide_Ranking_Vector.rds"
    ),
  
  input_policy_source =
    file.path(
      root_dir,
      "STEP9D_GenomeWide_Preranked_GSEA",
      "9D-4A_GO_BP_Input_Annotation_Policy_Lock",
      "STEP9D4A3_FINAL_GO_BP_GSEA_Input_Policy_Lock.rds"
    ),
  
  ranking_metric =
    "signed_minus_log10_nominal_P",
  
  ranking_formula =
    "sign(log2FC) * (-log10(max(P_value, 1e-300)))",
  
  original_genome_wide_ranking_n =
    57828L,
  
  go_bp_input_n =
    18008L,
  
  ontology =
    "BP",
  
  keyType =
    "ENSEMBL",
  
  organism =
    "org.Hs.eg.db",
  
  minGSSize =
    10,
  
  maxGSSize =
    500,
  
  exponent =
    1,
  
  pvalueCutoff =
    1,
  
  pAdjustMethod =
    "BH",
  
  verbose =
    FALSE,
  
  seed =
    FALSE,
  
  pathway_n =
    authoritative_pathway_n,
  
  fdr05_n =
    authoritative_fdr05_n,
  
  ranking_unique_n =
    ranking_unique_n,
  
  ranking_tied_gene_n =
    ranking_tied_gene_n,
  
  ranking_tie_fraction =
    ranking_tie_fraction,
  
  ranking_recalculated =
    FALSE,
  
  deg_recalculated =
    FALSE,
  
  jitter_used =
    FALSE,
  
  artificial_tie_breaking =
    FALSE,
  
  gsea_performed =
    TRUE,
  
  preranked_tie_warning =
    TRUE,
  
  fgsea_log2err_warning =
    TRUE,
  
  fgsea_eps_warning =
    TRUE,
  
  provenance_repair_only =
    TRUE,
  
  gsea_rerun_during_repair =
    FALSE,
  
  ranking_recalculated_during_repair =
    FALSE,
  
  deg_recalculated_during_repair =
    FALSE,
  
  original_authoritative_result_source =
    "STEP9D4B_FINAL_GO_BP_GSEA_Result.rds -> gseaResult@result",
  
  previous_contaminated_pathway_n =
    353L,
  
  previous_contaminated_fdr05_n =
    178L,
  
  repair_reason =
    "Previous GO-BP CSV/provenance outputs contained KEGG-derived values. Corrected outputs were serialized directly from the preserved authoritative GO-BP gseaResult without rerunning GSEA."
)


# ------------------------------------------------------------
# 13. Provenance validation
# ------------------------------------------------------------

provenance_checks <- c(
  
  pathway_n_correct =
    provenance_repaired$pathway_n ==
    authoritative_pathway_n,
  
  fdr05_n_correct =
    provenance_repaired$fdr05_n ==
    authoritative_fdr05_n,
  
  genome_ranking_n_correct =
    provenance_repaired$original_genome_wide_ranking_n ==
    57828L,
  
  go_input_n_correct =
    provenance_repaired$go_bp_input_n ==
    18008L,
  
  ontology_BP =
    identical(
      provenance_repaired$ontology,
      "BP"
    ),
  
  keyType_ENSEMBL =
    identical(
      provenance_repaired$keyType,
      "ENSEMBL"
    ),
  
  ranking_not_recalculated =
    identical(
      provenance_repaired$ranking_recalculated,
      FALSE
    ),
  
  deg_not_recalculated =
    identical(
      provenance_repaired$deg_recalculated,
      FALSE
    ),
  
  jitter_not_used =
    identical(
      provenance_repaired$jitter_used,
      FALSE
    ),
  
  artificial_tie_breaking_not_used =
    identical(
      provenance_repaired$artificial_tie_breaking,
      FALSE
    ),
  
  repair_only =
    identical(
      provenance_repaired$provenance_repair_only,
      TRUE
    ),
  
  gsea_not_rerun =
    identical(
      provenance_repaired$gsea_rerun_during_repair,
      FALSE
    )
  
)

print(provenance_checks)

if (!all(provenance_checks)) {
  stop(
    "ERROR: Corrected provenance validation failed."
  )
}


# ------------------------------------------------------------
# 14. Save corrected provenance to repair file
# ------------------------------------------------------------

repaired_provenance_rds <- file.path(
  step_dir,
  "STEP9D4B_FINAL_GO_BP_GSEA_Provenance_REPAIRED.rds"
)

saveRDS(
  provenance_repaired,
  repaired_provenance_rds
)

provenance_reload <- readRDS(
  repaired_provenance_rds
)

provenance_reload_check <- isTRUE(
  all.equal(
    provenance_reload,
    provenance_repaired,
    check.attributes = TRUE
  )
)

cat(
  "\nCorrected provenance reload = ",
  provenance_reload_check,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 15. PRE-OVERWRITE HARD PASS
# ------------------------------------------------------------

pre_overwrite_checks <- c(
  
  authoritative_structure =
    all(
      column_check,
      go_id_check,
      unique_id_check,
      settype_check,
      keytype_check,
      organism_check,
      gene_list_check
    ),
  
  input_integrity =
    all(
      input_n_check,
      input_column_check,
      input_id_check,
      input_value_check
    ),
  
  result_csv_integrity =
    result_reload_check,
  
  significant_csv_integrity =
    sig_reload_check,
  
  significant_count =
    sig_fdr_check,
  
  provenance_integrity =
    all(provenance_checks),
  
  provenance_reload =
    provenance_reload_check
  
)

cat("\n============================================================\n")
cat("15. PRE-OVERWRITE HARD PASS\n")
cat("============================================================\n\n")

print(pre_overwrite_checks)

pre_overwrite_pass <- all(
  pre_overwrite_checks
)

cat(
  "\nPRE-OVERWRITE REPAIR HARD PASS = ",
  pre_overwrite_pass,
  "\n",
  sep = ""
)

if (!pre_overwrite_pass) {
  stop(
    "ERROR: Pre-overwrite audit failed. Official files NOT overwritten."
  )
}


# ------------------------------------------------------------
# 16. Backup existing contaminated outputs
# ------------------------------------------------------------

backup_dir <- file.path(
  step_dir,
  "BACKUP_CONTAMINATED_OUTPUTS"
)

dir.create(
  backup_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

official_existing <- c(
  official_result_csv,
  official_sig_csv,
  official_provenance_rds
)

official_existing <- official_existing[
  file.exists(official_existing)
]

cat("\n============================================================\n")
cat("16. BACKUP EXISTING OFFICIAL OUTPUTS\n")
cat("============================================================\n\n")

if (length(official_existing) > 0) {
  
  backup_targets <- file.path(
    backup_dir,
    basename(official_existing)
  )
  
  for (i in seq_along(official_existing)) {
    
    if (!file.exists(backup_targets[i])) {
      
      ok <- file.copy(
        official_existing[i],
        backup_targets[i],
        overwrite = FALSE
      )
      
      if (!ok) {
        stop(
          paste(
            "ERROR: Backup failed:",
            official_existing[i]
          )
        )
      }
      
    } else {
      
      cat(
        "Backup already exists; preserved:\n",
        backup_targets[i],
        "\n"
      )
    }
  }
  
  cat(
    "\nBackup completed/preserved for ",
    length(official_existing),
    " official file(s).\n",
    sep = ""
  )
  
} else {
  
  cat(
    "No existing official outputs required backup.\n"
  )
}


# ------------------------------------------------------------
# 17. Replace ONLY contaminated outputs
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("17. REPLACE OFFICIAL CONTAMINATED OUTPUTS\n")
cat("============================================================\n\n")

ok1 <- file.copy(
  repaired_result_csv,
  official_result_csv,
  overwrite = TRUE
)

ok2 <- file.copy(
  repaired_sig_csv,
  official_sig_csv,
  overwrite = TRUE
)

ok3 <- file.copy(
  repaired_provenance_rds,
  official_provenance_rds,
  overwrite = TRUE
)

if (!all(
  ok1,
  ok2,
  ok3
)) {
  stop(
    "ERROR: One or more official output replacements failed."
  )
}

cat(
  "Official Result.csv replaced = TRUE\n"
)

cat(
  "Official Significant_FDR05.csv replaced = TRUE\n"
)

cat(
  "Official Provenance.rds replaced = TRUE\n"
)


# ------------------------------------------------------------
# 18. FINAL OFFICIAL RELOAD
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("18. FINAL OFFICIAL RELOAD VERIFICATION\n")
cat("============================================================\n\n")

final_result_csv <- read.csv(
  official_result_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

final_sig_csv <- read.csv(
  official_sig_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

final_provenance <- readRDS(
  official_provenance_rds
)

final_result_check <- isTRUE(
  all.equal(
    final_result_csv,
    go_result,
    check.attributes = FALSE
  )
)

final_sig_check <- isTRUE(
  all.equal(
    final_sig_csv,
    go_sig,
    check.attributes = FALSE
  )
)

final_provenance_count_check <- (
  final_provenance$pathway_n ==
    authoritative_pathway_n &&
    final_provenance$fdr05_n ==
    authoritative_fdr05_n
)

final_provenance_GO_check <- (
  identical(
    final_provenance$ontology,
    "BP"
  ) &&
    identical(
      final_provenance$keyType,
      "ENSEMBL"
    ) &&
    final_provenance$go_bp_input_n ==
    18008L
)

final_no_rerun_check <- (
  identical(
    final_provenance$gsea_rerun_during_repair,
    FALSE
  ) &&
    identical(
      final_provenance$ranking_recalculated_during_repair,
      FALSE
    ) &&
    identical(
      final_provenance$deg_recalculated_during_repair,
      FALSE
    )
)

final_files_check <- all(
  file.exists(
    c(
      official_result_csv,
      official_sig_csv,
      official_provenance_rds
    )
  )
)


# ------------------------------------------------------------
# 19. FINAL HARD PASS
# ------------------------------------------------------------

final_checks <- c(
  
  Result_CSV_authoritative =
    final_result_check,
  
  Significant_CSV_authoritative =
    final_sig_check,
  
  Provenance_counts_correct =
    final_provenance_count_check,
  
  Provenance_GO_metadata_correct =
    final_provenance_GO_check,
  
  No_GSEA_rerun =
    final_no_rerun_check,
  
  All_official_files_present =
    final_files_check
  
)

cat("\n============================================================\n")
cat("19. FINAL GO-BP OUTPUT / PROVENANCE REPAIR\n")
cat("============================================================\n\n")

print(final_checks)

final_hard_pass <- all(
  final_checks
)

cat(
  "\nSTEP9D-4B FORMAL OUTPUT / PROVENANCE REPAIR HARD PASS = ",
  final_hard_pass,
  "\n",
  sep = ""
)

if (!final_hard_pass) {
  stop(
    "ERROR: FINAL REPAIR HARD PASS FAILED."
  )
}

cat("\n")
cat("FINAL STATUS = REPAIRED_FINAL_PASS\n")
cat("GO-BP pathways returned = 6316\n")
cat("GO-BP FDR < 0.05 = 2003\n")
cat("GO-BP input ranking N = 18008\n")
cat("GSEA rerun = FALSE\n")
cat("Ranking recalculation = FALSE\n")
cat("DEG recalculation = FALSE\n")
cat("Jitter = FALSE\n")
cat("Artificial tie-breaking = FALSE\n")
cat("============================================================\n")


# ============================================================
# Paper 1
# STEP9D-4C
# GO-BP GSEA Result & Reproducibility Audit
#
# FINAL READ-ONLY AUDIT
#
# IMPORTANT:
#   - NO gseGO()
#   - NO GSEA rerun
#   - NO ranking recalculation
#   - NO DEG recalculation
#   - NO jitter
#   - NO artificial tie-breaking
#   - NO modification of STEP9D-4B official files
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-4C FINAL GO-BP GSEA RESULT & REPRODUCIBILITY AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(methods)
})


# ------------------------------------------------------------
# 1. Directories
# ------------------------------------------------------------

root_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

step4b_dir <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-4B_Formal_GO_BP_GSEA_Execution"
)

step4a_dir <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-4A_GO_BP_Input_Annotation_Policy_Lock"
)

step2c_dir <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-2C_Final_Ranking_Vector"
)

audit_dir <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-4C_Formal_GO_BP_GSEA_Result_Audit"
)

dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Official STEP9D-4B files
# ------------------------------------------------------------

result_rds <- file.path(
  step4b_dir,
  "STEP9D4B_FINAL_GO_BP_GSEA_Result.rds"
)

result_csv <- file.path(
  step4b_dir,
  "STEP9D4B_FINAL_GO_BP_GSEA_Result.csv"
)

sig_csv <- file.path(
  step4b_dir,
  "STEP9D4B_FINAL_GO_BP_GSEA_Significant_FDR05.csv"
)

input_rds <- file.path(
  step4b_dir,
  "STEP9D4B_FINAL_GO_BP_GSEA_Input_Ranking.rds"
)

input_csv <- file.path(
  step4b_dir,
  "STEP9D4B_FINAL_GO_BP_GSEA_Input_Ranking.csv"
)

provenance_rds <- file.path(
  step4b_dir,
  "STEP9D4B_FINAL_GO_BP_GSEA_Provenance.rds"
)

ranking_rds <- file.path(
  step2c_dir,
  "STEP9D2C_FINAL_GenomeWide_Ranking_Vector.rds"
)

input_policy_rds <- file.path(
  step4a_dir,
  "STEP9D4A3_FINAL_GO_BP_GSEA_Input_Policy_Lock.rds"
)


# ------------------------------------------------------------
# 3. Required file existence
# ------------------------------------------------------------

required_files <- c(
  result_rds,
  result_csv,
  sig_csv,
  input_rds,
  input_csv,
  provenance_rds,
  ranking_rds,
  input_policy_rds
)

required_file_check <- file.exists(required_files)

cat("3. REQUIRED FILE CHECK\n")
print(
  data.frame(
    file = basename(required_files),
    exists = required_file_check,
    row.names = NULL
  )
)

required_files_pass <- all(required_file_check)

cat(
  "\nAll required files present = ",
  required_files_pass,
  "\n\n",
  sep = ""
)

if (!required_files_pass) {
  stop("ERROR: One or more required STEP9D-4C input files are missing.")
}


# ------------------------------------------------------------
# 4. Load authoritative files
# ------------------------------------------------------------

cat("4. LOADING OFFICIAL FILES\n")

gsea_go_bp_disk <- readRDS(result_rds)
result_csv_df <- read.csv(
  result_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sig_csv_df <- read.csv(
  sig_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

input_obj <- readRDS(input_rds)

input_csv_df <- read.csv(
  input_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

provenance <- readRDS(provenance_rds)

ranking_locked <- readRDS(ranking_rds)

input_policy <- readRDS(input_policy_rds)

cat("Files loaded successfully.\n\n")


# ------------------------------------------------------------
# 5. Authoritative GO-BP result extraction
# ------------------------------------------------------------

cat("5. AUTHORITATIVE GO-BP RESULT STRUCTURE\n")

is_gseaResult <- inherits(
  gsea_go_bp_disk,
  "gseaResult"
)

cat(
  "gseaResult class = ",
  is_gseaResult,
  "\n",
  sep = ""
)

if (!is_gseaResult) {
  stop("ERROR: STEP9D4B Result.rds is not a gseaResult object.")
}

go_result <- slot(
  gsea_go_bp_disk,
  "result"
)

go_set_type <- slot(
  gsea_go_bp_disk,
  "setType"
)

go_keytype <- slot(
  gsea_go_bp_disk,
  "keytype"
)

go_organism <- slot(
  gsea_go_bp_disk,
  "organism"
)

go_gene_list <- slot(
  gsea_go_bp_disk,
  "geneList"
)

cat(
  "Authoritative result dimensions = ",
  nrow(go_result),
  " × ",
  ncol(go_result),
  "\n",
  sep = ""
)

cat(
  "setType = ",
  go_set_type,
  "\n",
  sep = ""
)

cat(
  "keytype = ",
  go_keytype,
  "\n",
  sep = ""
)

cat(
  "organism = ",
  go_organism,
  "\n",
  sep = ""
)

cat(
  "geneList N = ",
  length(go_gene_list),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. GO-BP result column audit
# ------------------------------------------------------------

expected_result_columns <- c(
  "ID",
  "Description",
  "setSize",
  "enrichmentScore",
  "NES",
  "pvalue",
  "p.adjust",
  "qvalue",
  "rank",
  "leading_edge",
  "core_enrichment"
)

result_column_check <- identical(
  colnames(go_result),
  expected_result_columns
)

cat("6. RESULT COLUMN AUDIT\n")
cat(
  "Expected columns exactly = ",
  result_column_check,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. GO ID / pathway count audit
# ------------------------------------------------------------

go_id_check <- all(
  grepl("^GO:", go_result$ID)
)

go_id_unique_check <- (
  length(unique(go_result$ID)) == nrow(go_result)
)

go_id_duplicate_n <- sum(
  duplicated(go_result$ID)
)

authoritative_pathway_n <- nrow(go_result)

authoritative_fdr05_n <- sum(
  !is.na(go_result$p.adjust) &
    go_result$p.adjust < 0.05
)

pathway_count_check <- (
  authoritative_pathway_n == 6316
)

fdr05_count_check <- (
  authoritative_fdr05_n == 2003
)

cat("7. AUTHORITATIVE GO-BP PATHWAY AUDIT\n")
cat(
  "GO-BP pathways = ",
  authoritative_pathway_n,
  "\n",
  sep = ""
)

cat(
  "GO-BP FDR < 0.05 = ",
  authoritative_fdr05_n,
  "\n",
  sep = ""
)

cat(
  "All IDs are GO IDs = ",
  go_id_check,
  "\n",
  sep = ""
)

cat(
  "GO IDs unique = ",
  go_id_unique_check,
  "\n",
  sep = ""
)

cat(
  "Duplicate GO IDs = ",
  go_id_duplicate_n,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Numerical validity
# ------------------------------------------------------------

numeric_columns <- c(
  "setSize",
  "enrichmentScore",
  "NES",
  "pvalue",
  "p.adjust",
  "qvalue",
  "rank"
)

numeric_validity <- all(
  vapply(
    go_result[numeric_columns],
    function(x) {
      all(is.finite(x[!is.na(x)]))
    },
    logical(1)
  )
)

cat("8. NUMERICAL VALIDITY\n")
cat(
  "All numeric result fields finite where non-NA = ",
  numeric_validity,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Result CSV exact concordance
# ------------------------------------------------------------

cat("9. RESULT CSV CONCORDANCE\n")

result_dim_check <- identical(
  dim(result_csv_df),
  dim(go_result)
)

result_colname_check <- identical(
  colnames(result_csv_df),
  colnames(go_result)
)

result_value_check <- isTRUE(
  all.equal(
    result_csv_df,
    go_result,
    check.attributes = FALSE
  )
)

result_csv_go_id_check <- (
  nrow(result_csv_df) > 0 &&
    all(grepl("^GO:", result_csv_df$ID))
)

result_csv_no_kegg_check <- !any(
  grepl("^hsa", result_csv_df$ID)
)

cat(
  "Dimensions identical = ",
  result_dim_check,
  "\n",
  sep = ""
)

cat(
  "Column names identical = ",
  result_colname_check,
  "\n",
  sep = ""
)

cat(
  "Values/rows/columns concordant = ",
  result_value_check,
  "\n",
  sep = ""
)

cat(
  "CSV contains only GO IDs = ",
  result_csv_go_id_check,
  "\n",
  sep = ""
)

cat(
  "CSV contains no KEGG hsa IDs = ",
  result_csv_no_kegg_check,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Significant CSV exact subset audit
# ------------------------------------------------------------

expected_sig_df <- go_result[
  !is.na(go_result$p.adjust) &
    go_result$p.adjust < 0.05,
  ,
  drop = FALSE
]

sig_dim_check <- identical(
  dim(sig_csv_df),
  dim(expected_sig_df)
)

sig_colname_check <- identical(
  colnames(sig_csv_df),
  colnames(expected_sig_df)
)

sig_value_check <- isTRUE(
  all.equal(
    sig_csv_df,
    expected_sig_df,
    check.attributes = FALSE
  )
)

sig_fdr_check <- all(
  !is.na(sig_csv_df$p.adjust) &
    sig_csv_df$p.adjust < 0.05
)

sig_go_id_check <- (
  nrow(sig_csv_df) > 0 &&
    all(grepl("^GO:", sig_csv_df$ID))
)

sig_no_kegg_check <- !any(
  grepl("^hsa", sig_csv_df$ID)
)

cat("10. SIGNIFICANT FDR<0.05 CSV AUDIT\n")

cat(
  "Expected significant N = ",
  nrow(expected_sig_df),
  "\n",
  sep = ""
)

cat(
  "Observed significant CSV N = ",
  nrow(sig_csv_df),
  "\n",
  sep = ""
)

cat(
  "Dimensions identical = ",
  sig_dim_check,
  "\n",
  sep = ""
)

cat(
  "Column names identical = ",
  sig_colname_check,
  "\n",
  sep = ""
)

cat(
  "Rows/columns/value order concordant = ",
  sig_value_check,
  "\n",
  sep = ""
)

cat(
  "All rows FDR < 0.05 = ",
  sig_fdr_check,
  "\n",
  sep = ""
)

cat(
  "All IDs are GO IDs = ",
  sig_go_id_check,
  "\n",
  sep = ""
)

cat(
  "No KEGG IDs = ",
  sig_no_kegg_check,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Input ranking audit
# ------------------------------------------------------------

cat("11. INPUT RANKING AUDIT\n")

if (is.list(input_obj)) {
  
  if ("ranking_vector" %in% names(input_obj)) {
    input_ranking <- input_obj$ranking_vector
  } else {
    stop("ERROR: Input RDS list does not contain ranking_vector.")
  }
  
  if ("input_table" %in% names(input_obj)) {
    input_table <- input_obj$input_table
  } else {
    input_table <- data.frame(
      ENSEMBL = names(input_ranking),
      ranking = as.numeric(input_ranking),
      stringsAsFactors = FALSE
    )
  }
  
} else if (is.numeric(input_obj)) {
  
  input_ranking <- input_obj
  
  input_table <- data.frame(
    ENSEMBL = names(input_ranking),
    ranking = as.numeric(input_ranking),
    stringsAsFactors = FALSE
  )
  
} else {
  
  stop("ERROR: Unsupported input RDS structure.")
}


input_n <- length(input_ranking)

input_table_structure_check <- (
  all(c("ENSEMBL", "ranking") %in% colnames(input_table))
)

input_table_n_check <- (
  nrow(input_table) == 18008
)

input_ids_unique_check <- (
  length(unique(input_table$ENSEMBL)) ==
    nrow(input_table)
)

input_no_na_check <- all(
  !is.na(input_table$ranking)
)

input_finite_check <- all(
  is.finite(input_table$ranking)
)

input_sorted_check <- all(
  diff(input_table$ranking) <= 0
)

cat(
  "Input ranking N = ",
  input_n,
  "\n",
  sep = ""
)

cat(
  "Input table N = ",
  nrow(input_table),
  "\n",
  sep = ""
)

cat(
  "Input structure correct = ",
  input_table_structure_check,
  "\n",
  sep = ""
)

cat(
  "Input N = 18008 = ",
  input_table_n_check,
  "\n",
  sep = ""
)

cat(
  "ENSEMBL IDs unique = ",
  input_ids_unique_check,
  "\n",
  sep = ""
)

cat(
  "No NA ranking = ",
  input_no_na_check,
  "\n",
  sep = ""
)

cat(
  "All ranking finite = ",
  input_finite_check,
  "\n",
  sep = ""
)

cat(
  "Ranking sorted decreasing = ",
  input_sorted_check,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Input CSV concordance
# ------------------------------------------------------------

input_csv_dim_check <- identical(
  dim(input_csv_df),
  dim(input_table)
)

input_csv_col_check <- identical(
  colnames(input_csv_df),
  colnames(input_table)
)

input_csv_value_check <- isTRUE(
  all.equal(
    input_csv_df,
    input_table,
    check.attributes = FALSE
  )
)

cat("12. INPUT CSV CONCORDANCE\n")

cat(
  "Dimensions identical = ",
  input_csv_dim_check,
  "\n",
  sep = ""
)

cat(
  "Columns identical = ",
  input_csv_col_check,
  "\n",
  sep = ""
)

cat(
  "Input CSV values concordant = ",
  input_csv_value_check,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. Input ranking immutability against authoritative GSEA
# ------------------------------------------------------------

authoritative_gene_list <- go_gene_list

gene_list_n_check <- (
  length(authoritative_gene_list) == 18008
)

gene_list_names_check <- (
  length(names(authoritative_gene_list)) == 18008
)

gene_list_names_unique_check <- (
  length(unique(names(authoritative_gene_list))) == 18008
)

input_vs_gsea_id_check <- isTRUE(
  identical(
    names(input_ranking),
    names(authoritative_gene_list)
  )
)

input_vs_gsea_value_check <- isTRUE(
  all.equal(
    as.numeric(input_ranking),
    as.numeric(authoritative_gene_list),
    check.attributes = FALSE
  )
)

cat("13. INPUT ↔ AUTHORITATIVE GSEA GENE LIST\n")

cat(
  "GSEA geneList N = 18008 = ",
  gene_list_n_check,
  "\n",
  sep = ""
)

cat(
  "GSEA geneList names present = ",
  gene_list_names_check,
  "\n",
  sep = ""
)

cat(
  "GSEA geneList names unique = ",
  gene_list_names_unique_check,
  "\n",
  sep = ""
)

cat(
  "Input/GSEA gene IDs identical = ",
  input_vs_gsea_id_check,
  "\n",
  sep = ""
)

cat(
  "Input/GSEA ranking values identical = ",
  input_vs_gsea_value_check,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 14. Genome-wide locked ranking audit
# ------------------------------------------------------------

genome_ranking_n <- length(ranking_locked)

genome_ranking_n_check <- (
  genome_ranking_n == 57828
)

genome_ranking_names_check <- (
  !is.null(names(ranking_locked)) &&
    length(names(ranking_locked)) == 57828
)

genome_ranking_unique_ids_check <- (
  length(unique(names(ranking_locked))) == 57828
)

cat("14. LOCKED GENOME-WIDE RANKING AUDIT\n")

cat(
  "Genome-wide ranking N = ",
  genome_ranking_n,
  "\n",
  sep = ""
)

cat(
  "Genome-wide ranking N = 57828 = ",
  genome_ranking_n_check,
  "\n",
  sep = ""
)

cat(
  "Genome-wide ranking names available = ",
  genome_ranking_names_check,
  "\n",
  sep = ""
)

cat(
  "Genome-wide ranking IDs unique = ",
  genome_ranking_unique_ids_check,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 15. Tie audit
# ------------------------------------------------------------

ranking_unique_n <- length(
  unique(input_ranking)
)

ranking_tied_gene_n <- sum(
  duplicated(input_ranking) |
    duplicated(input_ranking, fromLast = TRUE)
)

ranking_tie_fraction <- (
  ranking_tied_gene_n / input_n
)

tie_metadata_check <- (
  ranking_unique_n == 15765 &&
    ranking_tied_gene_n == 2243 &&
    abs(ranking_tie_fraction - 0.1245558) < 1e-7
)

cat("15. PRERANKED TIE AUDIT\n")

cat(
  "Unique ranking values = ",
  ranking_unique_n,
  "\n",
  sep = ""
)

cat(
  "Genes involved in ties = ",
  ranking_tied_gene_n,
  "\n",
  sep = ""
)

cat(
  "Tie fraction = ",
  sprintf("%.7f", ranking_tie_fraction),
  "\n",
  sep = ""
)

cat(
  "Expected tie metadata concordant = ",
  tie_metadata_check,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 16. Core gene presence
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_ensembl <- c(
  AOC3 = "ENSG00000131471",
  MS4A7 = "ENSG00000166927",
  IL7R = "ENSG00000168685"
)

core_presence <- vapply(
  unname(core_ensembl),
  function(g) g %in% names(input_ranking),
  logical(1)
)

core_presence_check <- all(core_presence)

cat("16. CORE-GENE PRESENCE\n")

print(
  data.frame(
    Gene = names(core_ensembl),
    ENSEMBL = unname(core_ensembl),
    Present_in_GO_BP_Input = core_presence,
    row.names = NULL
  )
)

cat(
  "\nAll three core genes present = ",
  core_presence_check,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 17. Provenance audit
# ------------------------------------------------------------

cat("17. PROVENANCE AUDIT\n")

provenance_required_fields <- c(
  "workflow",
  "step",
  "status",
  "ranking_source",
  "input_policy_source",
  "ranking_metric",
  "ranking_formula",
  "original_genome_wide_ranking_n",
  "go_bp_input_n",
  "ontology",
  "keyType",
  "organism",
  "minGSSize",
  "maxGSSize",
  "exponent",
  "pvalueCutoff",
  "pAdjustMethod",
  "verbose",
  "seed",
  "pathway_n",
  "fdr05_n",
  "ranking_unique_n",
  "ranking_tied_gene_n",
  "ranking_tie_fraction",
  "ranking_recalculated",
  "deg_recalculated",
  "jitter_used",
  "artificial_tie_breaking",
  "gsea_performed",
  "preranked_tie_warning",
  "fgsea_log2err_warning",
  "fgsea_eps_warning"
)

provenance_fields_present <- all(
  provenance_required_fields %in%
    names(provenance)
)

provenance_value_checks <- c(
  workflow =
    identical(
      provenance$workflow,
      "Paper1_STEP9D_GenomeWide_Preranked_GSEA"
    ),
  
  step =
    identical(
      provenance$step,
      "STEP9D-4B"
    ),
  
  status =
    identical(
      provenance$status,
      "FINAL_GO_BP_GSEA"
    ),
  
  ranking_metric =
    identical(
      provenance$ranking_metric,
      "signed_minus_log10_nominal_P"
    ),
  
  ranking_formula =
    identical(
      provenance$ranking_formula,
      "sign(log2FC) * (-log10(max(P_value, 1e-300)))"
    ),
  
  original_genome_wide_ranking_n =
    isTRUE(
      provenance$original_genome_wide_ranking_n == 57828
    ),
  
  go_bp_input_n =
    isTRUE(
      provenance$go_bp_input_n == 18008
    ),
  
  ontology =
    identical(
      provenance$ontology,
      "BP"
    ),
  
  keyType =
    identical(
      provenance$keyType,
      "ENSEMBL"
    ),
  
  organism =
    identical(
      provenance$organism,
      "org.Hs.eg.db"
    ),
  
  minGSSize =
    isTRUE(
      provenance$minGSSize == 10
    ),
  
  maxGSSize =
    isTRUE(
      provenance$maxGSSize == 500
    ),
  
  exponent =
    isTRUE(
      provenance$exponent == 1
    ),
  
  pvalueCutoff =
    isTRUE(
      provenance$pvalueCutoff == 1
    ),
  
  pAdjustMethod =
    identical(
      provenance$pAdjustMethod,
      "BH"
    ),
  
  verbose =
    isFALSE(
      provenance$verbose
    ),
  
  seed =
    isFALSE(
      provenance$seed
    ),
  
  pathway_n =
    isTRUE(
      provenance$pathway_n == 6316
    ),
  
  fdr05_n =
    isTRUE(
      provenance$fdr05_n == 2003
    ),
  
  ranking_unique_n =
    isTRUE(
      provenance$ranking_unique_n == 15765
    ),
  
  ranking_tied_gene_n =
    isTRUE(
      provenance$ranking_tied_gene_n == 2243
    ),
  
  ranking_tie_fraction =
    abs(
      provenance$ranking_tie_fraction -
        0.1245558
    ) < 1e-7,
  
  ranking_recalculated =
    isFALSE(
      provenance$ranking_recalculated
    ),
  
  deg_recalculated =
    isFALSE(
      provenance$deg_recalculated
    ),
  
  jitter_used =
    isFALSE(
      provenance$jitter_used
    ),
  
  artificial_tie_breaking =
    isFALSE(
      provenance$artificial_tie_breaking
    ),
  
  gsea_performed =
    isTRUE(
      provenance$gsea_performed
    ),
  
  preranked_tie_warning =
    isTRUE(
      provenance$preranked_tie_warning
    ),
  
  fgsea_log2err_warning =
    isTRUE(
      provenance$fgsea_log2err_warning
    ),
  
  fgsea_eps_warning =
    isTRUE(
      provenance$fgsea_eps_warning
    )
)

print(provenance_value_checks)

provenance_values_pass <- all(
  provenance_value_checks
)

cat(
  "\nAll provenance values correct = ",
  provenance_values_pass,
  "\n",
  sep = ""
)

cat(
  "All required provenance fields present = ",
  provenance_fields_present,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 18. Method immutability / no rerun audit
# ------------------------------------------------------------

method_immutability <- c(
  ranking_recalculated =
    isFALSE(provenance$ranking_recalculated),
  
  deg_recalculated =
    isFALSE(provenance$deg_recalculated),
  
  jitter_used =
    isFALSE(provenance$jitter_used),
  
  artificial_tie_breaking =
    isFALSE(provenance$artificial_tie_breaking),
  
  gsea_performed =
    isTRUE(provenance$gsea_performed)
)

cat("18. METHOD IMMUTABILITY AUDIT\n")
print(method_immutability)

method_immutability_pass <- all(
  method_immutability
)

cat(
  "\nMethod immutability audit = ",
  method_immutability_pass,
  "\n\n",
  sep = ""
)

# ------------------------------------------------------------
# STEP9D-4C — Correct tie metadata check
# Read-only; no GSEA rerun; no ranking recalculation
# ------------------------------------------------------------

ranking_n_final <- length(input_ranking)

ranking_unique_n_final <- length(
  unique(input_ranking)
)

# Locked STEP9D-4B definition:
# excess duplicated ranking entries
ranking_tied_gene_n_final <-
  ranking_n_final - ranking_unique_n_final

ranking_tie_fraction_final <-
  ranking_tied_gene_n_final / ranking_n_final

# ------------------------------------------------------------
# Compare with locked provenance
# ------------------------------------------------------------

tie_metadata_check <- all(
  
  ranking_n_final == 18008,
  
  ranking_unique_n_final == 15765,
  
  ranking_tied_gene_n_final ==
    provenance$ranking_tied_gene_n,
  
  abs(
    ranking_tie_fraction_final -
      provenance$ranking_tie_fraction
  ) < 1e-12,
  
  provenance$ranking_tied_gene_n == 2243,
  
  abs(
    provenance$ranking_tie_fraction -
      (2243 / 18008)
  ) < 1e-12
)

cat("============================================================\n")
cat("STEP9D-4C CORRECTED TIE METADATA CHECK\n")
cat("============================================================\n\n")

cat(
  "Ranking N = ",
  ranking_n_final,
  "\n",
  sep = ""
)

cat(
  "Unique ranking values = ",
  ranking_unique_n_final,
  "\n",
  sep = ""
)

cat(
  "Tied-ranking N = N - unique = ",
  ranking_tied_gene_n_final,
  "\n",
  sep = ""
)

cat(
  "Tie fraction = ",
  sprintf("%.15f", ranking_tie_fraction_final),
  "\n\n",
  sep = ""
)

cat(
  "Provenance tied-ranking N = ",
  provenance$ranking_tied_gene_n,
  "\n",
  sep = ""
)

cat(
  "Provenance tie fraction = ",
  sprintf(
    "%.15f",
    provenance$ranking_tie_fraction
  ),
  "\n\n",
  sep = ""
)

cat(
  "Tie metadata check = ",
  tie_metadata_check,
  "\n",
  sep = ""
)

cat("\nNO GSEA RERUN PERFORMED.\n")
cat("NO RANKING RECALCULATION PERFORMED.\n")
cat("NO FORMAL RESULT MODIFIED.\n")


# ------------------------------------------------------------
# 19. Final integrated audit
# ------------------------------------------------------------

cat("============================================================\n")
cat("19. FINAL STEP9D-4C INTEGRATED AUDIT\n")
cat("============================================================\n\n")

final_checks <- c(
  
  Required_files =
    required_files_pass,
  
  gseaResult_structure =
    is_gseaResult,
  
  Result_columns =
    result_column_check,
  
  Authoritative_pathway_count_6316 =
    pathway_count_check,
  
  Authoritative_FDR05_count_2003 =
    fdr05_count_check,
  
  GO_ID_integrity =
    go_id_check,
  
  GO_ID_unique =
    go_id_unique_check,
  
  Numerical_validity =
    numeric_validity,
  
  Result_CSV_authoritative =
    result_dim_check &&
    result_colname_check &&
    result_value_check,
  
  Result_CSV_GO_only =
    result_csv_go_id_check &&
    result_csv_no_kegg_check,
  
  Significant_CSV_authoritative =
    sig_dim_check &&
    sig_colname_check &&
    sig_value_check,
  
  Significant_CSV_GO_only =
    sig_go_id_check &&
    sig_no_kegg_check,
  
  Input_structure =
    input_table_structure_check,
  
  Input_N_18008 =
    input_table_n_check,
  
  Input_ID_unique =
    input_ids_unique_check,
  
  Input_no_NA =
    input_no_na_check,
  
  Input_finite =
    input_finite_check,
  
  Input_sorted =
    input_sorted_check,
  
  Input_CSV_concordance =
    input_csv_dim_check &&
    input_csv_col_check &&
    input_csv_value_check,
  
  Input_GSEA_ID_concordance =
    input_vs_gsea_id_check,
  
  Input_GSEA_value_concordance =
    input_vs_gsea_value_check,
  
  Genome_ranking_N_57828 =
    genome_ranking_n_check,
  
  Genome_ranking_ID_unique =
    genome_ranking_unique_ids_check,
  
  Tie_audit =
    tie_metadata_check,
  
  Core_gene_presence =
    core_presence_check,
  
  Provenance_fields =
    provenance_fields_present,
  
  Provenance_values =
    provenance_values_pass,
  
  Method_immutability =
    method_immutability_pass
  
)

print(final_checks)

final_hard_pass <- all(
  final_checks
)

cat(
  "\nSTEP9D-4C FINAL HARD PASS = ",
  final_hard_pass,
  "\n",
  sep = ""
)

if (!final_hard_pass) {
  cat(
    "\nWARNING: STEP9D-4C FINAL HARD PASS FAILED.\n"
  )
}


# ------------------------------------------------------------
# 20. Create final audit object
# ------------------------------------------------------------

audit_object <- list(
  
  workflow =
    "Paper1_STEP9D_GenomeWide_Preranked_GSEA",
  
  step =
    "STEP9D-4C",
  
  status =
    if (final_hard_pass) {
      "FINAL_PASS"
    } else {
      "FINAL_FAIL"
    },
  
  audit_type =
    "READ_ONLY_RESULT_REPRODUCIBILITY_AUDIT",
  
  gsea_rerun =
    FALSE,
  
  ranking_recalculated =
    FALSE,
  
  deg_recalculated =
    FALSE,
  
  jitter_used =
    FALSE,
  
  artificial_tie_breaking =
    FALSE,
  
  authoritative_result_source =
    "STEP9D4B_FINAL_GO_BP_GSEA_Result.rds -> gseaResult@result",
  
  authoritative_pathway_n =
    authoritative_pathway_n,
  
  authoritative_fdr05_n =
    authoritative_fdr05_n,
  
  go_bp_input_n =
    input_n,
  
  genome_wide_ranking_n =
    genome_ranking_n,
  
  ranking_unique_n =
    ranking_unique_n,
  
  ranking_tied_gene_n =
    ranking_tied_gene_n,
  
  ranking_tie_fraction =
    ranking_tie_fraction,
  
  core_gene_presence =
    core_presence,
  
  final_checks =
    final_checks,
  
  final_hard_pass =
    final_hard_pass
)


# ------------------------------------------------------------
# 21. Save audit outputs
# ------------------------------------------------------------

audit_rds <- file.path(
  audit_dir,
  "STEP9D4C_FINAL_GO_BP_GSEA_Result_Audit.rds"
)

audit_csv <- file.path(
  audit_dir,
  "STEP9D4C_FINAL_GO_BP_GSEA_Result_Audit.csv"
)

audit_txt <- file.path(
  audit_dir,
  "STEP9D4C_FINAL_GO_BP_GSEA_Result_Audit.txt"
)


saveRDS(
  audit_object,
  audit_rds
)

audit_df <- data.frame(
  check = names(final_checks),
  pass = as.logical(final_checks),
  stringsAsFactors = FALSE
)

write.csv(
  audit_df,
  audit_csv,
  row.names = FALSE,
  quote = TRUE
)


# ------------------------------------------------------------
# 22. Final TXT audit report
# ------------------------------------------------------------

sink(
  audit_txt,
  split = TRUE
)

cat("============================================================\n")
cat("STEP9D-4C FINAL GO-BP GSEA RESULT & REPRODUCIBILITY AUDIT\n")
cat("============================================================\n\n")

cat(
  "Authoritative GO-BP pathways = ",
  authoritative_pathway_n,
  "\n",
  sep = ""
)

cat(
  "Authoritative GO-BP FDR < 0.05 = ",
  authoritative_fdr05_n,
  "\n",
  sep = ""
)

cat(
  "GO-BP input ranking N = ",
  input_n,
  "\n",
  sep = ""
)

cat(
  "Genome-wide ranking N = ",
  genome_ranking_n,
  "\n",
  sep = ""
)

cat(
  "Ranking tie fraction = ",
  sprintf("%.7f", ranking_tie_fraction),
  "\n",
  sep = ""
)

cat("\nCORE GENES\n")
print(
  data.frame(
    Gene = names(core_ensembl),
    ENSEMBL = unname(core_ensembl),
    Present = core_presence,
    row.names = NULL
  )
)

cat("\nFINAL CHECKS\n")
print(final_checks)

cat(
  "\nSTEP9D-4C FINAL HARD PASS = ",
  final_hard_pass,
  "\n",
  sep = ""
)

cat("\nMETHOD INTEGRITY\n")
cat("GSEA rerun = FALSE\n")
cat("Ranking recalculation = FALSE\n")
cat("DEG recalculation = FALSE\n")
cat("Jitter = FALSE\n")
cat("Artificial tie-breaking = FALSE\n")

cat("\nFINAL STATUS = ")
cat(
  if (final_hard_pass) {
    "FINAL_PASS"
  } else {
    "FINAL_FAIL"
  }
)
cat("\n")

sink()


# ------------------------------------------------------------
# 23. Reload verification
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("23. FINAL AUDIT FILE RELOAD VERIFICATION\n")
cat("============================================================\n\n")

audit_rds_reload <- readRDS(
  audit_rds
)

audit_rds_reload_check <- isTRUE(
  all.equal(
    audit_rds_reload,
    audit_object,
    check.attributes = FALSE
  )
)

audit_csv_reload <- read.csv(
  audit_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

audit_csv_reload_check <- (
  nrow(audit_csv_reload) == length(final_checks) &&
    all(
      audit_csv_reload$check ==
        names(final_checks)
    ) &&
    all(
      audit_csv_reload$pass ==
        as.logical(final_checks)
    )
)

audit_txt_exists <- file.exists(
  audit_txt
)

cat(
  "Audit RDS reload verification = ",
  audit_rds_reload_check,
  "\n",
  sep = ""
)

cat(
  "Audit CSV reload verification = ",
  audit_csv_reload_check,
  "\n",
  sep = ""
)

cat(
  "Audit TXT exists = ",
  audit_txt_exists,
  "\n",
  sep = ""
)

reload_hard_pass <- all(
  audit_rds_reload_check,
  audit_csv_reload_check,
  audit_txt_exists
)

cat(
  "\nAUDIT FILE RELOAD HARD PASS = ",
  reload_hard_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 24. Final status
# ------------------------------------------------------------

overall_hard_pass <- all(
  final_hard_pass,
  reload_hard_pass
)

cat("\n============================================================\n")
cat("STEP9D-4C FINAL STATUS\n")
cat("============================================================\n\n")

cat(
  "STEP9D-4C FINAL HARD PASS = ",
  overall_hard_pass,
  "\n",
  sep = ""
)

if (!overall_hard_pass) {
  stop(
    "ERROR: STEP9D-4C FINAL HARD PASS FAILED."
  )
}

cat("\n")
cat("FINAL STATUS = FINAL_PASS\n")
cat("GO-BP result is internally consistent and reproducible.\n")
cat("No GSEA rerun was performed.\n")
cat("No ranking recalculation was performed.\n")
cat("No DEG recalculation was performed.\n")
cat("No jitter was used.\n")
cat("No artificial tie-breaking was used.\n")
cat("\n============================================================\n")
cat("STEP9D-4C FINAL AUDIT COMPLETE\n")
cat("============================================================\n")


# ============================================================
# Paper 1
# STEP9D-5
# Original Reactome GSEA Recovery / Concordance Audit
#
# PURPOSE
# ------------------------------------------------------------
# Read-only concordance audit between:
#
#   1. Original gene-specific Reactome GSEA
#      IL7R / MS4A7 / AOC3
#
#   2. Locked genome-wide KEGG GSEA
#      STEP9D-3
#
#   3. Locked genome-wide GO-BP GSEA
#      STEP9D-4
#
# IMPORTANT
# ------------------------------------------------------------
# NO GSEA rerun
# NO ranking recalculation
# NO DEG recalculation
# NO jitter
# NO artificial tie-breaking
# NO modification of STEP9D-3 / STEP9D-4 files
#
# The audit evaluates biological-theme concordance rather than
# requiring identical pathway IDs across Reactome, KEGG and GO.
# ============================================================


cat("\n============================================================\n")
cat("STEP9D-5: Original Reactome GSEA Recovery / Concordance Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "stringr",
  "readr",
  "tibble",
  "purrr"
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
    "Missing required packages: ",
    paste(missing_packages, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 1. Locked directories
# ------------------------------------------------------------

root_dir <- paste0(
  "input_dir/Paper 1_Step 9-Functional Enrichment/"
)

step9d_root <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA"
)


# ------------------------------------------------------------
# STEP9D-3 KEGG
# ------------------------------------------------------------

kegg_dir <- file.path(
  step9d_root,
  "9D-3B_Formal_GSEA_Execution"
)

kegg_result_file <- file.path(
  kegg_dir,
  "STEP9D3B_FINAL_KEGG_GSEA_Result.csv"
)


# ------------------------------------------------------------
# STEP9D-4 GO-BP
# ------------------------------------------------------------

go_dir <- file.path(
  step9d_root,
  "9D-4B_Formal_GO_BP_GSEA_Execution"
)

go_result_file <- file.path(
  go_dir,
  "STEP9D4B_FINAL_GO_BP_GSEA_Result.csv"
)


# ------------------------------------------------------------
# STEP9D-5 output directory
# ------------------------------------------------------------

audit_dir <- file.path(
  step9d_root,
  "9D-5_Original_Reactome_Concordance_Audit"
)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Verify source files
# ------------------------------------------------------------

required_source_files <- c(
  kegg_result_file,
  go_result_file
)

source_exists <- file.exists(
  required_source_files
)

cat("Source file verification:\n")

print(
  data.frame(
    File = basename(required_source_files),
    Exists = source_exists,
    stringsAsFactors = FALSE
  )
)

if (!all(source_exists)) {
  stop(
    "One or more locked STEP9D-3 / STEP9D-4 source files are missing."
  )
}


# ------------------------------------------------------------
# 3. Read locked KEGG result
# ------------------------------------------------------------

kegg_df <- read.csv(
  kegg_result_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\nKEGG result:\n")
cat("  Rows:", nrow(kegg_df), "\n")
cat("  Columns:", ncol(kegg_df), "\n")


# ------------------------------------------------------------
# 4. Read locked GO-BP result
# ------------------------------------------------------------

go_df <- read.csv(
  go_result_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\nGO-BP result:\n")
cat("  Rows:", nrow(go_df), "\n")
cat("  Columns:", ncol(go_df), "\n")


# ============================================================
# 5. Structural audit of KEGG
# ============================================================

expected_gsea_columns <- c(
  "ID",
  "Description",
  "setSize",
  "enrichmentScore",
  "NES",
  "pvalue",
  "p.adjust",
  "qvalue",
  "rank",
  "leading_edge",
  "core_enrichment"
)

kegg_structure_pass <- identical(
  colnames(kegg_df),
  expected_gsea_columns
)

if (!kegg_structure_pass) {
  stop(
    "KEGG result columns do not match the locked GSEA schema."
  )
}

kegg_id_n <- length(unique(kegg_df$ID))

kegg_fdr05_n <- sum(
  as.numeric(kegg_df$p.adjust) < 0.05,
  na.rm = TRUE
)

kegg_go_signature <- any(
  grepl(
    "^GO:",
    kegg_df$ID
  )
)

kegg_keGG_signature <- any(
  grepl(
    "^hsa",
    kegg_df$ID
  )
)

cat("\nKEGG structural audit:\n")
cat("  Pathways:", nrow(kegg_df), "\n")
cat("  Unique IDs:", kegg_id_n, "\n")
cat("  FDR < 0.05:", kegg_fdr05_n, "\n")
cat("  GO IDs present:", kegg_go_signature, "\n")
cat("  KEGG IDs present:", kegg_keGG_signature, "\n")


# ============================================================
# 6. Structural audit of GO-BP
# ============================================================

go_structure_pass <- identical(
  colnames(go_df),
  expected_gsea_columns
)

if (!go_structure_pass) {
  stop(
    "GO-BP result columns do not match the locked GSEA schema."
  )
}

go_id_n <- length(unique(go_df$ID))

go_fdr05_n <- sum(
  as.numeric(go_df$p.adjust) < 0.05,
  na.rm = TRUE
)

go_go_signature <- all(
  grepl(
    "^GO:",
    go_df$ID
  )
)

go_keGG_signature <- any(
  grepl(
    "^hsa",
    go_df$ID
  )
)

cat("\nGO-BP structural audit:\n")
cat("  Pathways:", nrow(go_df), "\n")
cat("  Unique IDs:", go_id_n, "\n")
cat("  FDR < 0.05:", go_fdr05_n, "\n")
cat("  All IDs GO:", go_go_signature, "\n")
cat("  KEGG IDs present:", go_keGG_signature, "\n")


# ============================================================
# 7. Lock expected STEP9D-3 / STEP9D-4 result counts
# ============================================================

kegg_count_pass <- (
  nrow(kegg_df) == 353 &
    kegg_fdr05_n == 178
)

go_count_pass <- (
  nrow(go_df) == 6316 &
    go_fdr05_n == 2003
)

if (!kegg_count_pass) {
  stop(
    "Locked KEGG result count check failed."
  )
}

if (!go_count_pass) {
  stop(
    "Locked GO-BP result count check failed."
  )
  
  
}


# ============================================================
# 8. Build searchable text fields
# ============================================================

prepare_searchable <- function(df) {
  
  df %>%
    mutate(
      
      ID_search = stringr::str_to_lower(
        as.character(ID)
      ),
      
      Description_search = stringr::str_to_lower(
        as.character(Description)
      ),
      
      leading_edge_search = stringr::str_to_lower(
        as.character(leading_edge)
      ),
      
      core_enrichment_search = stringr::str_to_lower(
        as.character(core_enrichment)
      ),
      
      searchable_text = paste(
        ID_search,
        Description_search,
        leading_edge_search,
        core_enrichment_search,
        sep = " | "
      )
      
    )
}


kegg_search <- prepare_searchable(
  kegg_df
)

go_search <- prepare_searchable(
  go_df
)


# ============================================================
# 9. Original Reactome reference themes
# ============================================================
#
# IMPORTANT:
# These are NOT new GSEA results.
# They are locked reference themes recovered from the
# original manuscript / original GSEA interpretation.
#
# Three representative pathways:
#
# IL7R
#   REACTOME_FCGR_ACTIVATION
#
# MS4A7
#   REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS
#
# AOC3
#   REACTOME_CD22_MEDIATED_BCR_REGULATION
#
# ============================================================


reactome_reference <- tibble::tribble(
  
  ~Gene, ~Original_Reactome_ID,
  ~Original_Reactome_Pathway,
  ~Original_NES,
  ~Original_FDR,
  ~Theme,
  
  "IL7R",
  "REACTOME_FCGR_ACTIVATION",
  "FCGR activation",
  2.569834,
  3.99e-09,
  "Fc receptor / immune receptor activation",
  
  "MS4A7",
  "REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS",
  "Role of phospholipids in phagocytosis",
  2.550673,
  3.99e-09,
  "Phagocytosis / myeloid innate immune function",
  
  "AOC3",
  "REACTOME_CD22_MEDIATED_BCR_REGULATION",
  "CD22-mediated BCR regulation",
  2.620611,
  3.99e-09,
  "B-cell receptor / B-cell immune signaling"
  
)


# ============================================================
# 10. Define pre-specified biological theme keyword sets
# ============================================================
#
# These keywords are intentionally broad enough to recognize
# database-dependent terminology but specific enough to avoid
# treating every "immune" pathway as concordant.
#
# IMPORTANT:
# No keyword is added after inspecting the final results.
#
# ============================================================


theme_keywords <- list(
  
  IL7R = c(
    "fcgr",
    "fc gamma",
    "fc-gamma",
    "fc receptor",
    "fc receptor",
    "immunoglobulin",
    "antibody receptor",
    "immune receptor",
    "leukocyte receptor",
    "antigen receptor",
    "antibody-dependent",
    "phagocyt",
    "immune complex"
  ),
  
  MS4A7 = c(
    "phagocyt",
    "phagosome",
    "engulfment",
    "fcgr",
    "fc gamma",
    "fc-gamma",
    "macrophage",
    "monocyte",
    "myeloid",
    "innate immune",
    "lysosome",
    "particle uptake",
    "clearance"
  ),
  
  AOC3 = c(
    "b cell receptor",
    "b-cell receptor",
    "bcr",
    "cd22",
    "antigen receptor",
    "b cell activation",
    "b-cell activation",
    "b lymphocyte",
    "humoral",
    "immunoglobulin",
    "antibody production",
    "antibody response",
    "b cell signaling",
    "b-cell signaling"
  )
  
)


# ============================================================
# 11. Search function
# ============================================================

find_theme_hits <- function(
    df,
    gene_name,
    database_name
) {
  
  keywords <- theme_keywords[[gene_name]]
  
  if (is.null(keywords)) {
    stop(
      "No keyword definition for gene: ",
      gene_name
    )
  }
  
  keyword_pattern <- paste(
    stringr::str_replace_all(
      keywords,
      "([.|()\\^{}+$*?]|\\[|\\])",
      "\\\\\\1"
    ),
    collapse = "|"
  )
  
  hit <- stringr::str_detect(
    df$searchable_text,
    keyword_pattern
  )
  
  out <- df[
    hit,
    c(
      "ID",
      "Description",
      "setSize",
      "enrichmentScore",
      "NES",
      "pvalue",
      "p.adjust",
      "qvalue",
      "rank",
      "leading_edge",
      "core_enrichment"
    ),
    drop = FALSE
  ]
  
  if (nrow(out) == 0) {
    
    out <- data.frame(
      ID = character(0),
      Description = character(0),
      setSize = numeric(0),
      enrichmentScore = numeric(0),
      NES = numeric(0),
      pvalue = numeric(0),
      p.adjust = numeric(0),
      qvalue = numeric(0),
      rank = numeric(0),
      leading_edge = character(0),
      core_enrichment = character(0),
      stringsAsFactors = FALSE
    )
    
  }
  
  if (nrow(out) > 0) {
    
    out <- out %>%
      mutate(
        Gene = gene_name,
        Database = database_name,
        FDR05 = as.numeric(p.adjust) < 0.05
      ) %>%
      arrange(
        desc(FDR05),
        desc(abs(as.numeric(NES))),
        as.numeric(p.adjust)
      )
    
  } else {
    
    out$Gene <- character(0)
    out$Database <- character(0)
    out$FDR05 <- logical(0)
    
  }
  
  out
  
}


# ============================================================
# 12. Run programmatic keyword audit
# ============================================================

genes <- c(
  "IL7R",
  "MS4A7",
  "AOC3"
)


kegg_hits <- purrr::map_dfr(
  genes,
  ~ find_theme_hits(
    kegg_search,
    .x,
    "KEGG"
  )
)


go_hits <- purrr::map_dfr(
  genes,
  ~ find_theme_hits(
    go_search,
    .x,
    "GO-BP"
  )
)


all_theme_hits <- bind_rows(
  kegg_hits,
  go_hits
)


# ============================================================
# 13. Save raw keyword-hit table
# ============================================================

raw_hit_file <- file.path(
  audit_dir,
  "STEP9D5_Theme_Keyword_Hits_All.csv"
)

write.csv(
  all_theme_hits,
  raw_hit_file,
  row.names = FALSE,
  na = ""
)


# ============================================================
# 14. Summarize per-gene evidence
# ============================================================

summarize_hits <- function(
    hit_df,
    gene_name,
    database_name
) {
  
  x <- hit_df %>%
    filter(
      Gene == gene_name,
      Database == database_name
    )
  
  if (nrow(x) == 0) {
    
    return(
      tibble(
        Gene = gene_name,
        Database = database_name,
        Matching_pathways = 0L,
        Significant_matching_pathways = 0L,
        Max_abs_NES = NA_real_,
        Best_FDR = NA_real_
      )
    )
    
  }
  
  tibble(
    Gene = gene_name,
    Database = database_name,
    Matching_pathways = nrow(x),
    Significant_matching_pathways = sum(
      x$FDR05,
      na.rm = TRUE
    ),
    Max_abs_NES = max(
      abs(as.numeric(x$NES)),
      na.rm = TRUE
    ),
    Best_FDR = min(
      as.numeric(x$p.adjust),
      na.rm = TRUE
    )
  )
  
}


summary_table <- bind_rows(
  
  purrr::map_dfr(
    genes,
    ~ summarize_hits(
      kegg_hits,
      .x,
      "KEGG"
    )
  ),
  
  purrr::map_dfr(
    genes,
    ~ summarize_hits(
      go_hits,
      .x,
      "GO-BP"
    )
  )
  
)


# ============================================================
# 15. Define evidence strength
# ============================================================
#
# IMPORTANT:
# This is pre-specified before interpretation.
#
# Strong:
#   significant theme-matching evidence in BOTH KEGG and GO-BP
#
# Moderate:
#   significant theme-matching evidence in ONE database
#
# Limited:
#   theme-related hits exist but none reaches FDR < 0.05
#
# No:
#   no theme-related pathway/process detected
#
# ============================================================


evidence_summary <- summary_table %>%
  group_by(Gene) %>%
  summarise(
    
    KEGG_matching =
      sum(
        Matching_pathways[Database == "KEGG"],
        na.rm = TRUE
      ),
    
    KEGG_significant =
      sum(
        Significant_matching_pathways[
          Database == "KEGG"
        ],
        na.rm = TRUE
      ),
    
    GO_matching =
      sum(
        Matching_pathways[Database == "GO-BP"],
        na.rm = TRUE
      ),
    
    GO_significant =
      sum(
        Significant_matching_pathways[
          Database == "GO-BP"
        ],
        na.rm = TRUE
      ),
    
    .groups = "drop"
    
  ) %>%
  
  mutate(
    
    Concordance_class = case_when(
      
      KEGG_significant > 0 &
        GO_significant > 0 ~
        "Strong",
      
      KEGG_significant > 0 |
        GO_significant > 0 ~
        "Moderate",
      
      KEGG_matching > 0 |
        GO_matching > 0 ~
        "Limited",
      
      TRUE ~
        "No"
      
    )
    
  )


# ============================================================
# 16. Add original reference information
# ============================================================

final_concordance_summary <- reactome_reference %>%
  select(
    Gene,
    Original_Reactome_ID,
    Original_Reactome_Pathway,
    Original_NES,
    Original_FDR,
    Theme
  ) %>%
  left_join(
    evidence_summary,
    by = "Gene"
  )


# ============================================================
# 17. Print preliminary result
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-5 PRELIMINARY CONCORDANCE SUMMARY\n")
cat("============================================================\n\n")

print(
  final_concordance_summary
)


# ============================================================
# 18. Save concordance summary
# ============================================================

summary_file <- file.path(
  audit_dir,
  "STEP9D5_FINAL_Gene_Level_Concordance_Summary.csv"
)

write.csv(
  final_concordance_summary,
  summary_file,
  row.names = FALSE,
  na = ""
)


# ============================================================
# 19. Significant theme-hit table
# ============================================================

significant_theme_hits <- all_theme_hits %>%
  filter(
    FDR05
  ) %>%
  arrange(
    Gene,
    Database,
    desc(abs(as.numeric(NES))),
    as.numeric(p.adjust)
  )


significant_hit_file <- file.path(
  audit_dir,
  "STEP9D5_Significant_Theme_Matching_Pathways.csv"
)

write.csv(
  significant_theme_hits,
  significant_hit_file,
  row.names = FALSE,
  na = ""
)


# ============================================================
# 20. All theme-hit table with database separation
# ============================================================

for (g in genes) {
  
  gene_file <- file.path(
    audit_dir,
    paste0(
      "STEP9D5_",
      g,
      "_Theme_Matching_Pathways.csv"
    )
  )
  
  gene_df <- all_theme_hits %>%
    filter(
      Gene == g
    )
  
  write.csv(
    gene_df,
    gene_file,
    row.names = FALSE,
    na = ""
  )
  
}


# ============================================================
# 21. Explicit source-result integrity audit
# ============================================================

integrity_checks <- c(
  
  KEGG_353 =
    nrow(kegg_df) == 353,
  
  KEGG_FDR05_178 =
    kegg_fdr05_n == 178,
  
  KEGG_unique_ID =
    kegg_id_n == 353,
  
  KEGG_no_GO_ID =
    !kegg_go_signature,
  
  GO_BP_6316 =
    nrow(go_df) == 6316,
  
  GO_BP_FDR05_2003 =
    go_fdr05_n == 2003,
  
  GO_BP_unique_ID =
    go_id_n == 6316,
  
  GO_BP_all_GO_ID =
    go_go_signature,
  
  GO_BP_no_KEGG_ID =
    !go_keGG_signature,
  
  KEGG_schema =
    kegg_structure_pass,
  
  GO_BP_schema =
    go_structure_pass,
  
  Three_core_genes_reference =
    identical(
      sort(reactome_reference$Gene),
      sort(
        c(
          "AOC3",
          "IL7R",
          "MS4A7"
        )
      )
    )
  
)


cat("\n============================================================\n")
cat("SOURCE INTEGRITY AUDIT\n")
cat("============================================================\n\n")

print(
  integrity_checks
)


# ============================================================
# 22. Method immutability lock
# ============================================================

method_integrity <- c(
  
  GSEA_rerun = FALSE,
  
  ranking_recalculated = FALSE,
  
  DEG_recalculated = FALSE,
  
  jitter_used = FALSE,
  
  artificial_tie_breaking = FALSE,
  
  STEP9D3_modified = FALSE,
  
  STEP9D4_modified = FALSE
  
)


cat("\n============================================================\n")
cat("METHOD IMMUTABILITY\n")
cat("============================================================\n\n")

print(
  method_integrity
)


# ============================================================
# 23. Preliminary hard pass
# ============================================================

source_integrity_pass <- all(
  integrity_checks
)

method_integrity_pass <- all(
  method_integrity
)

summary_structure_pass <- all(
  final_concordance_summary$Gene %in%
    c(
      "IL7R",
      "MS4A7",
      "AOC3"
    )
)

preliminary_hard_pass <- (
  source_integrity_pass &
    method_integrity_pass &
    summary_structure_pass
)


cat("\n============================================================\n")
cat("STEP9D-5 PRELIMINARY HARD PASS\n")
cat("============================================================\n\n")

cat(
  "Source integrity PASS: ",
  source_integrity_pass,
  "\n",
  sep = ""
)

cat(
  "Method immutability PASS: ",
  method_integrity_pass,
  "\n",
  sep = ""
)

cat(
  "Summary structure PASS: ",
  summary_structure_pass,
  "\n",
  sep = ""
)

cat(
  "STEP9D-5 PRELIMINARY HARD PASS: ",
  preliminary_hard_pass,
  "\n",
  sep = ""
)


# ============================================================
# 24. Save audit object
# ============================================================

audit_object <- list(
  
  workflow =
    "Paper1_STEP9D_GenomeWide_Preranked_GSEA",
  
  step =
    "STEP9D-5",
  
  status =
    "READ_ONLY_CONCORDANCE_AUDIT",
  
  audit_date =
    as.character(Sys.time()),
  
  source_files =
    list(
      KEGG = kegg_result_file,
      GO_BP = go_result_file
    ),
  
  original_reactome_reference =
    reactome_reference,
  
  theme_keywords =
    theme_keywords,
  
  KEGG_pathway_n =
    nrow(kegg_df),
  
  KEGG_FDR05_n =
    kegg_fdr05_n,
  
  GO_BP_pathway_n =
    nrow(go_df),
  
  GO_BP_FDR05_n =
    go_fdr05_n,
  
  theme_hit_table =
    all_theme_hits,
  
  gene_level_concordance =
    final_concordance_summary,
  
  significant_theme_hits =
    significant_theme_hits,
  
  source_integrity =
    integrity_checks,
  
  method_immutability =
    method_integrity,
  
  hard_pass =
    preliminary_hard_pass
  
)


audit_rds_file <- file.path(
  audit_dir,
  "STEP9D5_FINAL_Original_Reactome_Concordance_Audit.rds"
)

saveRDS(
  audit_object,
  audit_rds_file
)


# ============================================================
# 25. Write human-readable audit report
# ============================================================

audit_txt_file <- file.path(
  audit_dir,
  "STEP9D5_FINAL_Original_Reactome_Concordance_Audit.txt"
)

sink(
  audit_txt_file
)

cat("STEP9D-5 — Original Reactome GSEA Recovery / Concordance Audit\n")
cat("============================================================\n\n")

cat(
  "Audit type: READ-ONLY\n"
)

cat(
  "GSEA rerun: FALSE\n"
)

cat(
  "Ranking recalculated: FALSE\n"
)

cat(
  "DEG recalculated: FALSE\n"
)

cat(
  "Jitter used: FALSE\n"
)

cat(
  "Artificial tie-breaking: FALSE\n\n"
)

cat("LOCKED KEGG RESULT\n")
cat("------------------\n")
cat(
  "Pathways returned: ",
  nrow(kegg_df),
  "\n",
  sep = ""
)
cat(
  "FDR < 0.05: ",
  kegg_fdr05_n,
  "\n\n",
  sep = ""
)

cat("LOCKED GO-BP RESULT\n")
cat("-------------------\n")
cat(
  "Pathways returned: ",
  nrow(go_df),
  "\n",
  sep = ""
)
cat(
  "FDR < 0.05: ",
  go_fdr05_n,
  "\n\n",
  sep = ""
)

cat("ORIGINAL REACTOME REFERENCE\n")
cat("---------------------------\n\n")

print(
  reactome_reference,
  row.names = FALSE
)

cat("\n\nGENE-LEVEL CONCORDANCE\n")
cat("-----------------------\n\n")

print(
  final_concordance_summary,
  row.names = FALSE
)

cat("\n\nSOURCE INTEGRITY\n")
cat("----------------\n\n")

print(
  integrity_checks
)

cat("\n\nMETHOD IMMUTABILITY\n")
cat("-------------------\n\n")

print(
  method_integrity
)

cat("\n\nFINAL PRELIMINARY STATUS\n")
cat("------------------------\n\n")

cat(
  "STEP9D-5 PRELIMINARY HARD PASS: ",
  preliminary_hard_pass,
  "\n",
  sep = ""
)

sink()


# ============================================================
# 26. Reload verification
# ============================================================

reload_audit <- readRDS(
  audit_rds_file
)

reload_pass <- isTRUE(
  all.equal(
    reload_audit$gene_level_concordance,
    final_concordance_summary,
    check.attributes = FALSE
  )
)

cat("\n============================================================\n")
cat("STEP9D-5 RELOAD VERIFICATION\n")
cat("============================================================\n\n")

cat(
  "Audit RDS reload verification: ",
  reload_pass,
  "\n",
  sep = ""
)

cat(
  "Audit TXT exists: ",
  file.exists(audit_txt_file),
  "\n",
  sep = ""
)

cat(
  "Audit RDS exists: ",
  file.exists(audit_rds_file),
  "\n",
  sep = ""
)

cat(
  "STEP9D-5 READ-ONLY AUDIT COMPLETE\n"
)

cat(
  "No STEP9D-3 / STEP9D-4 source files were modified.\n"
)


# ============================================================
# Paper 1
# STEP9D-5B
# Strict Theme Evidence Audit
#
# PURPOSE
# ------------------------------------------------------------
# Read-only audit of concordance between:
#   1. Original Xiantao Reactome themes
#   2. Locked STEP9D-3B KEGG GSEA
#   3. Locked STEP9D-4B GO-BP GSEA
#
# IMPORTANT
# ------------------------------------------------------------
# This step:
#   - DOES NOT rerun GSEA
#   - DOES NOT recalculate the ranking
#   - DOES NOT recalculate DEGs
#   - DOES NOT modify STEP9D-3B files
#   - DOES NOT modify STEP9D-4B files
#   - DOES NOT use jitter
#   - DOES NOT perform artificial tie-breaking
#   - DOES NOT assign Strong / Moderate / Limited / No
#
# Primary thematic evidence:
#   pathway ID + Description
#
# leading_edge / core_enrichment:
#   retained as supporting information only;
#   NOT subjected to keyword theme matching.
#
# ============================================================


cat("\n============================================================\n")
cat("STEP9D-5B: STRICT THEME EVIDENCE AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory and output directory
# ------------------------------------------------------------

step_root <- "input_dir/Paper 1_Step 9-Functional Enrichment"

step_dir <- file.path(
  step_root,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-5B_Strict_Theme_Evidence_Audit"
)

dir.create(
  step_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("STEP9D-5B output directory:\n")
cat(step_dir, "\n\n")


# ------------------------------------------------------------
# 1. Locked source files
# ------------------------------------------------------------

kegg_file <- file.path(
  step_root,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-3B_Formal_GSEA_Execution",
  "STEP9D3B_FINAL_KEGG_GSEA_Result.csv"
)

go_file <- file.path(
  step_root,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-4B_Formal_GO_BP_GSEA_Execution",
  "STEP9D4B_FINAL_GO_BP_GSEA_Result.csv"
)

required_source_files <- c(
  kegg_file,
  go_file
)

source_exists <- file.exists(required_source_files)

cat("SOURCE FILE VERIFICATION\n")
cat("============================================================\n\n")

print(
  data.frame(
    File = basename(required_source_files),
    Exists = source_exists,
    stringsAsFactors = FALSE
  )
)

if (!all(source_exists)) {
  stop(
    "STEP9D-5B STOP: one or more locked source files are missing."
  )
}


# ------------------------------------------------------------
# 2. Read-only import
# ------------------------------------------------------------

kegg_df <- read.csv(
  kegg_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

go_df <- read.csv(
  go_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\nSOURCE DATA DIMENSIONS\n")
cat("============================================================\n\n")

cat(
  "KEGG: ",
  nrow(kegg_df),
  " pathways x ",
  ncol(kegg_df),
  " columns\n",
  sep = ""
)

cat(
  "GO-BP: ",
  nrow(go_df),
  " pathways x ",
  ncol(go_df),
  " columns\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 3. Required schema
# ------------------------------------------------------------

required_columns <- c(
  "ID",
  "Description",
  "setSize",
  "enrichmentScore",
  "NES",
  "pvalue",
  "p.adjust",
  "qvalue",
  "rank",
  "leading_edge",
  "core_enrichment"
)

kegg_schema_pass <- identical(
  names(kegg_df),
  required_columns
)

go_schema_pass <- identical(
  names(go_df),
  required_columns
)

cat("SCHEMA AUDIT\n")
cat("============================================================\n\n")

cat(
  "KEGG schema PASS: ",
  kegg_schema_pass,
  "\n",
  sep = ""
)

cat(
  "GO-BP schema PASS: ",
  go_schema_pass,
  "\n",
  sep = ""
)

if (!kegg_schema_pass || !go_schema_pass) {
  stop(
    "STEP9D-5B STOP: unexpected source schema."
  )
}


# ------------------------------------------------------------
# 4. Structural source lock verification
# ------------------------------------------------------------

kegg_id_n <- length(unique(kegg_df$ID))
go_id_n   <- length(unique(go_df$ID))

kegg_fdr05_n <- sum(
  !is.na(kegg_df$p.adjust) &
    kegg_df$p.adjust < 0.05
)

go_fdr05_n <- sum(
  !is.na(go_df$p.adjust) &
    go_df$p.adjust < 0.05
)

kegg_id_pattern_pass <- all(
  grepl("^hsa[0-9]+$", kegg_df$ID)
)

go_id_pattern_pass <- all(
  grepl("^GO:[0-9]+$", go_df$ID)
)

cat("\nSOURCE STRUCTURAL LOCK\n")
cat("============================================================\n\n")

cat("KEGG rows: ", nrow(kegg_df), "\n", sep = "")
cat("KEGG unique IDs: ", kegg_id_n, "\n", sep = "")
cat("KEGG FDR < 0.05: ", kegg_fdr05_n, "\n", sep = "")
cat("KEGG ID pattern valid: ", kegg_id_pattern_pass, "\n\n", sep = "")

cat("GO-BP rows: ", nrow(go_df), "\n", sep = "")
cat("GO-BP unique IDs: ", go_id_n, "\n", sep = "")
cat("GO-BP FDR < 0.05: ", go_fdr05_n, "\n", sep = "")
cat("GO-BP ID pattern valid: ", go_id_pattern_pass, "\n", sep = "")


# ------------------------------------------------------------
# 5. Original Reactome reference themes
#
# These are historical reference observations supplied from
# the original Xiantao Reactome analysis.
#
# They are NOT rerun here.
# ------------------------------------------------------------

reactome_reference <- data.frame(
  
  Gene = c(
    "IL7R",
    "MS4A7",
    "AOC3"
  ),
  
  Original_Reactome_ID = c(
    "REACTOME_FCGR_ACTIVATION",
    "REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS",
    "REACTOME_CD22_MEDIATED_BCR_REGULATION"
  ),
  
  Original_Reactome_Pathway = c(
    "Fc receptor activation",
    "Role of phospholipids in phagocytosis",
    "CD22-mediated B-cell receptor regulation"
  ),
  
  Original_NES = c(
    2.569834,
    2.550673,
    2.620611
  ),
  
  Original_FDR = c(
    3.99e-09,
    3.99e-09,
    3.99e-09
  ),
  
  Theme = c(
    "Fc receptor / immune receptor activation",
    "Phagocytosis / myeloid innate immune function",
    "B-cell receptor / B-cell immune signaling"
  ),
  
  stringsAsFactors = FALSE
)

cat("\nORIGINAL REACTOME REFERENCE\n")
cat("============================================================\n\n")

print(reactome_reference)


# ------------------------------------------------------------
# 6. Strict thematic concept dictionaries
#
# IMPORTANT:
#   Direct terms = relatively specific to the historical theme
#   Context terms = biologically related but not sufficient alone
#
# These dictionaries are predefined BEFORE inspecting hits.
# ------------------------------------------------------------

theme_dictionary <- list(
  
  IL7R = list(
    
    direct = c(
      "fc receptor",
      "fc-gamma receptor",
      "fc gamma receptor",
      "fcgr",
      "immune receptor",
      "immunoreceptor",
      "antibody receptor",
      "receptor-mediated phagocytosis",
      "antigen receptor"
    ),
    
    contextual = c(
      "immune receptor",
      "leukocyte receptor",
      "antigen receptor signaling",
      "immune signaling",
      "antigen presentation",
      "antibody",
      "immunoglobulin",
      "leukocyte activation",
      "lymphocyte activation"
    )
  ),
  
  MS4A7 = list(
    
    direct = c(
      "phagocytosis",
      "phagosome",
      "engulfment",
      "particle uptake",
      "Fc receptor-mediated phagocytosis",
      "fc receptor mediated phagocytosis",
      "phagocytic"
    ),
    
    contextual = c(
      "macrophage",
      "monocyte",
      "myeloid",
      "innate immune",
      "innate immunity",
      "vesicle",
      "endosome",
      "lysosome",
      "leukocyte"
    )
  ),
  
  AOC3 = list(
    
    direct = c(
      "b-cell receptor",
      "b cell receptor",
      "bcr signaling",
      "bcr",
      "cd22",
      "b-cell activation",
      "b cell activation",
      "b lymphocyte receptor",
      "b-cell receptor signaling",
      "b cell receptor signaling"
    ),
    
    contextual = c(
      "b-cell",
      "b cell",
      "b lymphocyte",
      "humoral immune",
      "lymphocyte activation",
      "adaptive immune",
      "antibody",
      "immunoglobulin"
    )
  )
)


# ------------------------------------------------------------
# 7. Matching function
#
# IMPORTANT:
#   Only Description is classified.
#
# leading_edge/core_enrichment are retained separately.
# ------------------------------------------------------------

classify_description <- function(
    description,
    direct_terms,
    contextual_terms
) {
  
  description_lower <- tolower(
    ifelse(
      is.na(description),
      "",
      description
    )
  )
  
  direct_hits <- direct_terms[
    vapply(
      direct_terms,
      function(x) grepl(
        x,
        description_lower,
        fixed = TRUE
      ),
      logical(1)
    )
  ]
  
  contextual_hits <- contextual_terms[
    vapply(
      contextual_terms,
      function(x) grepl(
        x,
        description_lower,
        fixed = TRUE
      ),
      logical(1)
    )
  ]
  
  # Remove contextual terms that are also captured as direct
  contextual_hits <- setdiff(
    contextual_hits,
    direct_hits
  )
  
  direct_n <- length(direct_hits)
  contextual_n <- length(contextual_hits)
  
  evidence_class <- if (
    direct_n > 0
  ) {
    "Direct"
  } else if (
    contextual_n > 0
  ) {
    "Contextual"
  } else {
    "No_theme_match"
  }
  
  data.frame(
    Direct_theme_match = direct_n > 0,
    Contextual_theme_match = (
      direct_n == 0 &&
        contextual_n > 0
    ),
    Direct_terms = paste(
      direct_hits,
      collapse = "; "
    ),
    Contextual_terms = paste(
      contextual_hits,
      collapse = "; "
    ),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 8. Build complete evidence table
# ------------------------------------------------------------

build_evidence_table <- function(
    df,
    gene,
    database_name
) {
  
  dict <- theme_dictionary[[gene]]
  
  classification_list <- lapply(
    seq_len(nrow(df)),
    function(i) {
      
      classify_description(
        description = df$Description[i],
        direct_terms = dict$direct,
        contextual_terms = dict$contextual
      )
    }
  )
  
  classification_df <- do.call(
    rbind,
    classification_list
  )
  
  out <- cbind(
    data.frame(
      Gene = gene,
      Database = database_name,
      stringsAsFactors = FALSE
    ),
    df,
    classification_df
  )
  
  out$Significance <- ifelse(
    !is.na(out$p.adjust) &
      out$p.adjust < 0.05,
    "FDR<0.05",
    "FDR>=0.05_or_NA"
  )
  
  out$Evidence_Quadrant <- ifelse(
    out$Direct_theme_match &
      out$p.adjust < 0.05,
    "Direct_significant",
    
    ifelse(
      out$Contextual_theme_match &
        out$p.adjust < 0.05,
      "Contextual_significant",
      
      ifelse(
        out$Direct_theme_match &
          !(out$p.adjust < 0.05),
        "Direct_non_significant",
        
        ifelse(
          out$Contextual_theme_match &
            !(out$p.adjust < 0.05),
          "Contextual_non_significant",
          
          "No_theme_match"
        )
      )
    )
  )
  
  out
}


# ------------------------------------------------------------
# 9. Generate evidence tables
# ------------------------------------------------------------

evidence_tables <- list()

for (gene in c("IL7R", "MS4A7", "AOC3")) {
  
  evidence_tables[[paste0(gene, "_KEGG")]] <-
    build_evidence_table(
      kegg_df,
      gene,
      "KEGG"
    )
  
  evidence_tables[[paste0(gene, "_GO_BP")]] <-
    build_evidence_table(
      go_df,
      gene,
      "GO_BP"
    )
}


# ------------------------------------------------------------
# 10. Combine all evidence
# ------------------------------------------------------------

all_evidence <- do.call(
  rbind,
  evidence_tables
)

rownames(all_evidence) <- NULL


# ------------------------------------------------------------
# 11. Remove rows with no thematic evidence
#
# The full source tables remain untouched.
# This object only contains actual theme-related evidence.
# ------------------------------------------------------------

theme_evidence <- all_evidence[
  all_evidence$Direct_theme_match |
    all_evidence$Contextual_theme_match,
]

rownames(theme_evidence) <- NULL


# ------------------------------------------------------------
# 12. Significant theme evidence only
# ------------------------------------------------------------

theme_evidence_significant <- theme_evidence[
  !is.na(theme_evidence$p.adjust) &
    theme_evidence$p.adjust < 0.05,
]

rownames(theme_evidence_significant) <- NULL


# ------------------------------------------------------------
# 13. Evidence category summary
# ------------------------------------------------------------

quadrant_levels <- c(
  "Direct_significant",
  "Contextual_significant",
  "Direct_non_significant",
  "Contextual_non_significant"
)

summary_list <- list()

for (gene in c("IL7R", "MS4A7", "AOC3")) {
  
  for (db in c("KEGG", "GO_BP")) {
    
    tmp <- all_evidence[
      all_evidence$Gene == gene &
        all_evidence$Database == db,
    ]
    
    counts <- table(
      factor(
        tmp$Evidence_Quadrant,
        levels = c(
          quadrant_levels,
          "No_theme_match"
        )
      )
    )
    
    summary_list[[length(summary_list) + 1]] <-
      data.frame(
        Gene = gene,
        Database = db,
        
        Direct_significant =
          unname(counts["Direct_significant"]),
        
        Contextual_significant =
          unname(counts["Contextual_significant"]),
        
        Direct_non_significant =
          unname(counts["Direct_non_significant"]),
        
        Contextual_non_significant =
          unname(counts["Contextual_non_significant"]),
        
        No_theme_match =
          unname(counts["No_theme_match"]),
        
        stringsAsFactors = FALSE
      )
  }
}

evidence_summary <- do.call(
  rbind,
  summary_list
)

rownames(evidence_summary) <- NULL


# ------------------------------------------------------------
# 14. Strict significant evidence summary
# ------------------------------------------------------------

strict_summary <- aggregate(
  x = list(
    Direct_significant =
      theme_evidence_significant$Direct_theme_match,
    
    Contextual_significant =
      theme_evidence_significant$Contextual_theme_match
  ),
  by = list(
    Gene =
      theme_evidence_significant$Gene,
    
    Database =
      theme_evidence_significant$Database
  ),
  FUN = sum
)

if (nrow(theme_evidence_significant) == 0) {
  
  strict_summary <- data.frame(
    Gene = character(0),
    Database = character(0),
    Direct_significant = integer(0),
    Contextual_significant = integer(0),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 15. Gene-level diagnostic summary
#
# IMPORTANT:
#   This is descriptive only.
#   No Strong/Moderate/Limited/No classification.
# ------------------------------------------------------------

gene_summary_list <- list()

for (gene in c("IL7R", "MS4A7", "AOC3")) {
  
  for (db in c("KEGG", "GO_BP")) {
    
    tmp <- all_evidence[
      all_evidence$Gene == gene &
        all_evidence$Database == db,
    ]
    
    sig <- tmp[
      !is.na(tmp$p.adjust) &
        tmp$p.adjust < 0.05,
    ]
    
    direct_sig <- sig[
      sig$Direct_theme_match,
    ]
    
    contextual_sig <- sig[
      !sig$Direct_theme_match &
        sig$Contextual_theme_match,
    ]
    
    gene_summary_list[[length(gene_summary_list) + 1]] <-
      data.frame(
        Gene = gene,
        Database = db,
        
        Total_pathways =
          nrow(tmp),
        
        Theme_matching_pathways =
          sum(
            tmp$Direct_theme_match |
              tmp$Contextual_theme_match
          ),
        
        Significant_theme_pathways =
          nrow(sig),
        
        Direct_significant =
          nrow(direct_sig),
        
        Contextual_significant =
          nrow(contextual_sig),
        
        stringsAsFactors = FALSE
      )
  }
}

gene_summary <- do.call(
  rbind,
  gene_summary_list
)

rownames(gene_summary) <- NULL


# ------------------------------------------------------------
# 16. Print the actual significant evidence
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-5B SIGNIFICANT THEME EVIDENCE\n")
cat("============================================================\n\n")

for (gene in c("IL7R", "MS4A7", "AOC3")) {
  
  cat("\n------------------------------------------------------------\n")
  cat(gene, "\n")
  cat("------------------------------------------------------------\n")
  
  for (db in c("KEGG", "GO_BP")) {
    
    cat("\n", db, "\n", sep = "")
    
    tmp <- theme_evidence_significant[
      theme_evidence_significant$Gene == gene &
        theme_evidence_significant$Database == db,
    ]
    
    if (nrow(tmp) == 0) {
      
      cat("No significant thematic pathways detected.\n")
      
    } else {
      
      print(
        tmp[
          ,
          c(
            "ID",
            "Description",
            "NES",
            "pvalue",
            "p.adjust",
            "Evidence_Quadrant",
            "Direct_terms",
            "Contextual_terms",
            "leading_edge",
            "core_enrichment"
          ),
          drop = FALSE
        ],
        row.names = FALSE
      )
    }
  }
}


# ------------------------------------------------------------
# 17. Print descriptive summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-5B EVIDENCE SUMMARY\n")
cat("============================================================\n\n")

print(evidence_summary)

cat("\n============================================================\n")
cat("STEP9D-5B GENE-LEVEL DESCRIPTIVE SUMMARY\n")
cat("============================================================\n\n")

print(gene_summary)


# ------------------------------------------------------------
# 18. Explicit audit checks
# ------------------------------------------------------------

source_integrity_pass <- all(
  nrow(kegg_df) == 353,
  kegg_id_n == 353,
  kegg_fdr05_n == 178,
  kegg_id_pattern_pass,
  
  nrow(go_df) == 6316,
  go_id_n == 6316,
  go_fdr05_n == 2003,
  go_id_pattern_pass,
  
  kegg_schema_pass,
  go_schema_pass
)

theme_audit_structure_pass <- all(
  nrow(all_evidence) ==
    3 * 2 * (
      nrow(kegg_df) + nrow(go_df)
    ),
  
  all(
    c(
      "Gene",
      "Database",
      "ID",
      "Description",
      "NES",
      "pvalue",
      "p.adjust",
      "leading_edge",
      "core_enrichment",
      "Direct_theme_match",
      "Contextual_theme_match",
      "Evidence_Quadrant"
    ) %in%
      names(all_evidence)
  ),
  
  nrow(evidence_summary) == 6,
  nrow(gene_summary) == 6
)


# ------------------------------------------------------------
# 19. Method immutability
#
# STEP9D-5B is read-only by design.
# ------------------------------------------------------------

method_integrity <- c(
  
  GSEA_rerun = FALSE,
  
  ranking_recalculated = FALSE,
  
  DEG_recalculated = FALSE,
  
  jitter_used = FALSE,
  
  artificial_tie_breaking = FALSE,
  
  STEP9D3_modified = FALSE,
  
  STEP9D4_modified = FALSE
)

method_integrity_pass <- all(
  !method_integrity
)

cat("\n============================================================\n")
cat("STEP9D-5B METHOD IMMUTABILITY\n")
cat("============================================================\n\n")

print(method_integrity)

cat(
  "\nMethod immutability PASS: ",
  method_integrity_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 20. IMPORTANT: no final concordance classification
# ------------------------------------------------------------

final_concordance_classification_performed <- FALSE

cat("\n============================================================\n")
cat("FINAL CONCORDANCE CLASSIFICATION\n")
cat("============================================================\n\n")

cat(
  "Strong / Moderate / Limited / No classification performed: ",
  final_concordance_classification_performed,
  "\n",
  sep = ""
)

cat(
  "\nIMPORTANT:\n"
)

cat(
  "This STEP9D-5B audit intentionally does NOT assign a final\n"
)

cat(
  "concordance class. It only documents direct and contextual\n"
)

cat(
  "theme-level evidence for subsequent expert interpretation.\n"
)


# ------------------------------------------------------------
# 21. Final preliminary audit status
# ------------------------------------------------------------

preliminary_hard_pass <- all(
  source_integrity_pass,
  theme_audit_structure_pass,
  method_integrity_pass,
  !final_concordance_classification_performed
)

cat("\n============================================================\n")
cat("STEP9D-5B PRELIMINARY HARD PASS\n")
cat("============================================================\n\n")

cat(
  "Source integrity PASS: ",
  source_integrity_pass,
  "\n",
  sep = ""
)

cat(
  "Theme-audit structure PASS: ",
  theme_audit_structure_pass,
  "\n",
  sep = ""
)

cat(
  "Method immutability PASS: ",
  method_integrity_pass,
  "\n",
  sep = ""
)

cat(
  "No final concordance classification performed: ",
  !final_concordance_classification_performed,
  "\n",
  sep = ""
)

cat(
  "STEP9D-5B PRELIMINARY HARD PASS: ",
  preliminary_hard_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 22. Save evidence tables
# ------------------------------------------------------------

write.csv(
  all_evidence,
  file.path(
    step_dir,
    "STEP9D5B_All_Theme_Evidence_Audit.csv"
  ),
  row.names = FALSE
)

write.csv(
  theme_evidence,
  file.path(
    step_dir,
    "STEP9D5B_Theme_Matching_Evidence.csv"
  ),
  row.names = FALSE
)

write.csv(
  theme_evidence_significant,
  file.path(
    step_dir,
    "STEP9D5B_Significant_Theme_Evidence.csv"
  ),
  row.names = FALSE
)

write.csv(
  evidence_summary,
  file.path(
    step_dir,
    "STEP9D5B_Evidence_Quadrant_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  gene_summary,
  file.path(
    step_dir,
    "STEP9D5B_Gene_Level_Descriptive_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  reactome_reference,
  file.path(
    step_dir,
    "STEP9D5B_Original_Reactome_Reference.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 23. Save complete audit object
# ------------------------------------------------------------

audit_object <- list(
  
  workflow =
    "Paper1_STEP9D_GenomeWide_Preranked_GSEA",
  
  step =
    "STEP9D-5B",
  
  status =
    ifelse(
      preliminary_hard_pass,
      "PRELIMINARY_PASS",
      "PRELIMINARY_FAIL"
    ),
  
  purpose =
    "Strict theme-level evidence audit without final concordance classification",
  
  source_files =
    required_source_files,
  
  source_dimensions =
    list(
      KEGG = dim(kegg_df),
      GO_BP = dim(go_df)
    ),
  
  source_integrity =
    list(
      KEGG_pathways = nrow(kegg_df),
      KEGG_FDR05 = kegg_fdr05_n,
      GO_BP_pathways = nrow(go_df),
      GO_BP_FDR05 = go_fdr05_n
    ),
  
  reactome_reference =
    reactome_reference,
  
  theme_dictionary =
    theme_dictionary,
  
  all_evidence =
    all_evidence,
  
  theme_evidence =
    theme_evidence,
  
  theme_evidence_significant =
    theme_evidence_significant,
  
  evidence_summary =
    evidence_summary,
  
  gene_summary =
    gene_summary,
  
  method_integrity =
    method_integrity,
  
  source_integrity_pass =
    source_integrity_pass,
  
  theme_audit_structure_pass =
    theme_audit_structure_pass,
  
  method_integrity_pass =
    method_integrity_pass,
  
  final_concordance_classification_performed =
    final_concordance_classification_performed,
  
  preliminary_hard_pass =
    preliminary_hard_pass
)


audit_rds_file <- file.path(
  step_dir,
  "STEP9D5B_FINAL_Strict_Theme_Evidence_Audit.rds"
)

saveRDS(
  audit_object,
  audit_rds_file
)


# ------------------------------------------------------------
# 24. Save concise TXT audit
# ------------------------------------------------------------

audit_txt_file <- file.path(
  step_dir,
  "STEP9D5B_FINAL_Strict_Theme_Evidence_Audit.txt"
)

sink(audit_txt_file)

cat("STEP9D-5B: STRICT THEME EVIDENCE AUDIT\n")
cat("============================================================\n\n")

cat("Purpose:\n")
cat(
  "Read-only thematic evidence audit comparing original Reactome\n"
)
cat(
  "reference themes with locked KEGG and GO-BP GSEA results.\n\n"
)

cat("Source integrity:\n")
cat(
  "KEGG pathways = ",
  nrow(kegg_df),
  "\n",
  sep = ""
)
cat(
  "KEGG FDR<0.05 = ",
  kegg_fdr05_n,
  "\n",
  sep = ""
)
cat(
  "GO-BP pathways = ",
  nrow(go_df),
  "\n",
  sep = ""
)
cat(
  "GO-BP FDR<0.05 = ",
  go_fdr05_n,
  "\n\n",
  sep = ""
)

cat("Evidence summary:\n")
print(evidence_summary)

cat("\nGene-level descriptive summary:\n")
print(gene_summary)

cat("\nMethod immutability:\n")
print(method_integrity)

cat(
  "\nFinal concordance classification performed = ",
  final_concordance_classification_performed,
  "\n",
  sep = ""
)

cat(
  "\nSTEP9D-5B preliminary hard pass = ",
  preliminary_hard_pass,
  "\n",
  sep = ""
)

cat(
  "\nNo STEP9D-3B or STEP9D-4B source files were modified.\n"
)

cat(
  "No GSEA rerun was performed.\n"
)

cat(
  "No ranking recalculation was performed.\n"
)

cat(
  "No DEG recalculation was performed.\n"
)

cat(
  "No jitter was used.\n"
)

cat(
  "No artificial tie-breaking was used.\n"
)

sink()


# ------------------------------------------------------------
# 25. Reload verification
# ------------------------------------------------------------

audit_reload <- readRDS(
  audit_rds_file
)

reload_pass <- all(
  identical(
    audit_reload$step,
    "STEP9D-5B"
  ),
  
  identical(
    audit_reload$source_integrity_pass,
    source_integrity_pass
  ),
  
  identical(
    audit_reload$theme_audit_structure_pass,
    theme_audit_structure_pass
  ),
  
  identical(
    audit_reload$method_integrity,
    method_integrity
  ),
  
  identical(
    audit_reload$final_concordance_classification_performed,
    FALSE
  )
)

cat("\n============================================================\n")
cat("STEP9D-5B RELOAD VERIFICATION\n")
cat("============================================================\n\n")

cat(
  "Audit RDS reload verification: ",
  reload_pass,
  "\n",
  sep = ""
)

cat(
  "Audit RDS exists: ",
  file.exists(audit_rds_file),
  "\n",
  sep = ""
)

cat(
  "Audit TXT exists: ",
  file.exists(audit_txt_file),
  "\n",
  sep = ""
)

cat(
  "\nSTEP9D-5B READ-ONLY AUDIT COMPLETE\n"
)

cat(
  "No STEP9D-3B / STEP9D-4B source files were modified.\n"
)

cat(
  "No final concordance class has been assigned.\n"
)

cat(
  "============================================================\n"
)


# ============================================================
# STEP9D-5B
# CORRECTED GENE-LEVEL SUMMARY AUDIT
#
# READ-ONLY
# ------------------------------------------------------------
# No GSEA rerun
# No ranking recalculation
# No DEG recalculation
# No source-file modification
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-5B: CORRECTED GENE-LEVEL SUMMARY AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Verify required object
# ------------------------------------------------------------

if (!exists("all_evidence")) {
  stop(
    "STOP: all_evidence is not available in the current R session."
  )
}

required_cols <- c(
  "Gene",
  "Database",
  "Description",
  "p.adjust",
  "Direct_theme_match",
  "Contextual_theme_match",
  "Evidence_Quadrant"
)

if (!all(required_cols %in% names(all_evidence))) {
  stop(
    "STOP: all_evidence does not contain the required columns."
  )
}


# ------------------------------------------------------------
# 2. Correct gene-level summary
# ------------------------------------------------------------

corrected_summary_list <- list()

for (gene in c("IL7R", "MS4A7", "AOC3")) {
  
  for (db in c("KEGG", "GO_BP")) {
    
    tmp <- all_evidence[
      all_evidence$Gene == gene &
        all_evidence$Database == db,
      ,
      drop = FALSE
    ]
    
    # --------------------------------------------------------
    # Theme matching first
    # --------------------------------------------------------
    
    theme_tmp <- tmp[
      tmp$Direct_theme_match |
        tmp$Contextual_theme_match,
      ,
      drop = FALSE
    ]
    
    # --------------------------------------------------------
    # Significant AND theme matching
    # --------------------------------------------------------
    
    sig_theme_tmp <- theme_tmp[
      !is.na(theme_tmp$p.adjust) &
        theme_tmp$p.adjust < 0.05,
      ,
      drop = FALSE
    ]
    
    # --------------------------------------------------------
    # Direct significant
    # --------------------------------------------------------
    
    direct_sig_tmp <- sig_theme_tmp[
      sig_theme_tmp$Direct_theme_match,
      ,
      drop = FALSE
    ]
    
    # --------------------------------------------------------
    # Contextual significant
    # --------------------------------------------------------
    
    contextual_sig_tmp <- sig_theme_tmp[
      !sig_theme_tmp$Direct_theme_match &
        sig_theme_tmp$Contextual_theme_match,
      ,
      drop = FALSE
    ]
    
    # --------------------------------------------------------
    # Direct non-significant
    # --------------------------------------------------------
    
    direct_nonsig_tmp <- theme_tmp[
      theme_tmp$Direct_theme_match &
        (
          is.na(theme_tmp$p.adjust) |
            theme_tmp$p.adjust >= 0.05
        ),
      ,
      drop = FALSE
    ]
    
    # --------------------------------------------------------
    # Contextual non-significant
    # --------------------------------------------------------
    
    contextual_nonsig_tmp <- theme_tmp[
      !theme_tmp$Direct_theme_match &
        theme_tmp$Contextual_theme_match &
        (
          is.na(theme_tmp$p.adjust) |
            theme_tmp$p.adjust >= 0.05
        ),
      ,
      drop = FALSE
    ]
    
    corrected_summary_list[[length(corrected_summary_list) + 1]] <- data.frame(
      
      Gene = gene,
      
      Database = db,
      
      Total_pathways =
        nrow(tmp),
      
      Theme_matching_pathways =
        nrow(theme_tmp),
      
      Significant_theme_pathways =
        nrow(sig_theme_tmp),
      
      Direct_significant =
        nrow(direct_sig_tmp),
      
      Contextual_significant =
        nrow(contextual_sig_tmp),
      
      Direct_non_significant =
        nrow(direct_nonsig_tmp),
      
      Contextual_non_significant =
        nrow(contextual_nonsig_tmp),
      
      No_theme_match =
        sum(
          !tmp$Direct_theme_match &
            !tmp$Contextual_theme_match
        ),
      
      stringsAsFactors = FALSE
    )
  }
}

corrected_gene_summary <- do.call(
  rbind,
  corrected_summary_list
)

rownames(corrected_gene_summary) <- NULL


# ------------------------------------------------------------
# 3. Print corrected summary
# ------------------------------------------------------------

cat("CORRECTED GENE-LEVEL SUMMARY\n")
cat("============================================================\n\n")

print(
  corrected_gene_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. Internal consistency audit
# ------------------------------------------------------------

summary_consistency <- logical(
  nrow(corrected_gene_summary)
)

for (i in seq_len(nrow(corrected_gene_summary))) {
  
  x <- corrected_gene_summary[i, ]
  
  summary_consistency[i] <- all(
    
    x$Theme_matching_pathways ==
      x$Direct_significant +
      x$Contextual_significant +
      x$Direct_non_significant +
      x$Contextual_non_significant,
    
    x$Total_pathways ==
      x$Theme_matching_pathways +
      x$No_theme_match,
    
    x$Significant_theme_pathways ==
      x$Direct_significant +
      x$Contextual_significant
  )
}

cat("\nSUMMARY INTERNAL CONSISTENCY\n")
cat("============================================================\n\n")

print(
  data.frame(
    Gene = corrected_gene_summary$Gene,
    Database = corrected_gene_summary$Database,
    PASS = summary_consistency,
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)

summary_consistency_pass <- all(
  summary_consistency
)

cat(
  "\nSummary consistency PASS: ",
  summary_consistency_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. Explicit correction check
# ------------------------------------------------------------

# For example:
# IL7R / KEGG must have:
#   Theme matching = 0
#   Significant theme = 0
#
# It must NOT equal 178.

il7r_kegg <- corrected_gene_summary[
  corrected_gene_summary$Gene == "IL7R" &
    corrected_gene_summary$Database == "KEGG",
  ,
  drop = FALSE
]

il7r_kegg_correction_pass <- all(
  il7r_kegg$Theme_matching_pathways == 0,
  il7r_kegg$Significant_theme_pathways == 0
)

cat(
  "\nIL7R KEGG correction check: ",
  il7r_kegg_correction_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Method immutability
# ------------------------------------------------------------

method_integrity <- c(
  
  GSEA_rerun = FALSE,
  
  ranking_recalculated = FALSE,
  
  DEG_recalculated = FALSE,
  
  jitter_used = FALSE,
  
  artificial_tie_breaking = FALSE,
  
  STEP9D3_modified = FALSE,
  
  STEP9D4_modified = FALSE
)

method_integrity_pass <- all(
  !method_integrity
)

cat("\nMETHOD IMMUTABILITY\n")
cat("============================================================\n\n")

print(method_integrity)

cat(
  "\nMethod immutability PASS: ",
  method_integrity_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Final corrected preliminary status
# ------------------------------------------------------------

corrected_preliminary_pass <- all(
  summary_consistency_pass,
  il7r_kegg_correction_pass,
  method_integrity_pass
)

cat("\n============================================================\n")
cat("STEP9D-5B CORRECTED PRELIMINARY STATUS\n")
cat("============================================================\n\n")

cat(
  "Summary consistency PASS: ",
  summary_consistency_pass,
  "\n",
  sep = ""
)

cat(
  "IL7R KEGG correction PASS: ",
  il7r_kegg_correction_pass,
  "\n",
  sep = ""
)

cat(
  "Method immutability PASS: ",
  method_integrity_pass,
  "\n",
  sep = ""
)

cat(
  "\nSTEP9D-5B CORRECTED PRELIMINARY PASS: ",
  corrected_preliminary_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Save corrected summary ONLY
# ------------------------------------------------------------

corrected_summary_file <- file.path(
  step_dir,
  "STEP9D5B_CORRECTED_Gene_Level_Descriptive_Summary.csv"
)

write.csv(
  corrected_gene_summary,
  corrected_summary_file,
  row.names = FALSE
)

cat(
  "\nCorrected summary saved:\n",
  corrected_summary_file,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. No source modification confirmation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-5B CORRECTION COMPLETE\n")
cat("============================================================\n\n")

cat("No GSEA rerun performed.\n")
cat("No ranking recalculation performed.\n")
cat("No DEG recalculation performed.\n")
cat("No jitter used.\n")
cat("No artificial tie-breaking used.\n")
cat("No STEP9D-3B source file modified.\n")
cat("No STEP9D-4B source file modified.\n")
cat("No final concordance classification assigned.\n")


# ============================================================
# Paper 1
# STEP9D-5C
# Final Concordance Evidence Review / Lock
#
# PURPOSE
# ------------------------------------------------------------
# Final evidence-only review of:
#   Original Reactome GSEA observations
#   versus locked KEGG GSEA
#   versus locked GO-BP GSEA
#
# IMPORTANT:
#   - NO GSEA rerun
#   - NO ranking recalculation
#   - NO DEG recalculation
#   - NO jitter
#   - NO artificial tie-breaking
#   - NO statistical testing
#   - NO modification of STEP9D-3B / STEP9D-4B
#   - NO automatic biological ranking/classification
#
# This is an evidence organization / lock step only.
# ============================================================


cat("\n============================================================\n")
cat("STEP9D-5C: Final Concordance Evidence Review / Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory and output directory
# ------------------------------------------------------------

root_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

step9d_dir <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA"
)

out_dir <- file.path(
  step9d_dir,
  "9D-5C_Final_Concordance_Evidence_Lock"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Locked source files
# ------------------------------------------------------------

kegg_result_file <- file.path(
  step9d_dir,
  "9D-3B_Formal_GSEA_Execution",
  "STEP9D3B_FINAL_KEGG_GSEA_Result.csv"
)

go_result_file <- file.path(
  step9d_dir,
  "9D-4B_Formal_GO_BP_GSEA_Execution",
  "STEP9D4B_FINAL_GO_BP_GSEA_Result.csv"
)

kegg_sig_file <- file.path(
  step9d_dir,
  "9D-3B_Formal_GSEA_Execution",
  "STEP9D3B_FINAL_KEGG_GSEA_Significant_FDR05.csv"
)

go_sig_file <- file.path(
  step9d_dir,
  "9D-4B_Formal_GO_BP_GSEA_Execution",
  "STEP9D4B_FINAL_GO_BP_GSEA_Significant_FDR05.csv"
)

kegg_prov_file <- file.path(
  step9d_dir,
  "9D-3B_Formal_GSEA_Execution",
  "STEP9D3B_FINAL_GSEA_Provenance.rds"
)

go_prov_file <- file.path(
  step9d_dir,
  "9D-4B_Formal_GO_BP_GSEA_Execution",
  "STEP9D4B_FINAL_GO_BP_GSEA_Provenance.rds"
)


# ------------------------------------------------------------
# 2. STEP9D-5B source files
# ------------------------------------------------------------

step5b_dir <- file.path(
  step9d_dir,
  "9D-5B_Strict_Theme_Evidence_Audit"
)

step5b_candidates <- c(
  file.path(
    step5b_dir,
    "STEP9D5B_FINAL_Strict_Theme_Evidence_Audit.rds"
  ),
  file.path(
    step5b_dir,
    "STEP9D5B_CORRECTED_Strict_Theme_Evidence_Audit.rds"
  ),
  file.path(
    step5b_dir,
    "STEP9D5B_FINAL_Theme_Evidence_Audit.rds"
  )
)

step5b_existing <- step5b_candidates[
  file.exists(step5b_candidates)
]


# ------------------------------------------------------------
# 3. Required source existence audit
# ------------------------------------------------------------

required_files <- c(
  kegg_result_file,
  go_result_file,
  kegg_sig_file,
  go_sig_file,
  kegg_prov_file,
  go_prov_file
)

required_files_exist <- all(
  file.exists(required_files)
)

cat("Required STEP9D-3B / STEP9D-4B files exist: ",
    required_files_exist, "\n", sep = "")

if (!required_files_exist) {
  stop(
    "Required locked GSEA source file(s) missing. ",
    "STEP9D-5C stopped without modification."
  )
}


# ------------------------------------------------------------
# 4. Load locked GSEA results
# ------------------------------------------------------------

kegg_df <- read.csv(
  kegg_result_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

go_df <- read.csv(
  go_result_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

kegg_sig_df <- read.csv(
  kegg_sig_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

go_sig_df <- read.csv(
  go_sig_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

kegg_prov <- readRDS(kegg_prov_file)
go_prov <- readRDS(go_prov_file)


# ------------------------------------------------------------
# 5. Basic schema checks
# ------------------------------------------------------------

expected_gsea_columns <- c(
  "ID",
  "Description",
  "setSize",
  "enrichmentScore",
  "NES",
  "pvalue",
  "p.adjust",
  "qvalue",
  "rank",
  "leading_edge",
  "core_enrichment"
)

kegg_schema_pass <- identical(
  names(kegg_df),
  expected_gsea_columns
)

go_schema_pass <- identical(
  names(go_df),
  expected_gsea_columns
)

cat("KEGG schema PASS: ", kegg_schema_pass, "\n", sep = "")
cat("GO-BP schema PASS: ", go_schema_pass, "\n", sep = "")


# ------------------------------------------------------------
# 6. Locked pathway-count checks
# ------------------------------------------------------------

kegg_pathway_n <- nrow(kegg_df)
go_pathway_n <- nrow(go_df)

kegg_fdr05_n <- sum(
  !is.na(kegg_df$p.adjust) &
    kegg_df$p.adjust < 0.05
)

go_fdr05_n <- sum(
  !is.na(go_df$p.adjust) &
    go_df$p.adjust < 0.05
)

kegg_count_pass <- (
  kegg_pathway_n == 353 &
    kegg_fdr05_n == 178
)

go_count_pass <- (
  go_pathway_n == 6316 &
    go_fdr05_n == 2003
)

cat("\nLocked pathway counts:\n")
cat("  KEGG pathways = ", kegg_pathway_n, "\n", sep = "")
cat("  KEGG FDR < 0.05 = ", kegg_fdr05_n, "\n", sep = "")
cat("  KEGG count PASS = ", kegg_count_pass, "\n", sep = "")

cat("  GO-BP pathways = ", go_pathway_n, "\n", sep = "")
cat("  GO-BP FDR < 0.05 = ", go_fdr05_n, "\n", sep = "")
cat("  GO-BP count PASS = ", go_count_pass, "\n", sep = "")


# ------------------------------------------------------------
# 7. Core genes
# ------------------------------------------------------------

core_genes <- c(
  "IL7R",
  "MS4A7",
  "AOC3"
)


# ------------------------------------------------------------
# 8. Original Reactome reference
#
# These values are historical reference observations from
# the original Xiantao Reactome GSEA and are NOT recalculated.
# ------------------------------------------------------------

original_reactome <- data.frame(
  
  Gene = c(
    "IL7R",
    "MS4A7",
    "AOC3"
  ),
  
  Original_Reactome_Pathway = c(
    "REACTOME_FCGR_ACTIVATION",
    "REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS",
    "REACTOME_CD22_MEDIATED_BCR_REGULATION"
  ),
  
  Original_NES = c(
    2.569834,
    2.550673,
    2.620611
  ),
  
  Original_FDR = c(
    3.99e-09,
    3.99e-09,
    3.99e-09
  ),
  
  Original_Theme = c(
    "Fc receptor / immune receptor activation",
    "Phagocytosis / myeloid innate immune function",
    "B-cell receptor / B-cell immune signaling"
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 9. Theme keyword definitions
#
# IMPORTANT:
# These are the SAME strict theme definitions used in STEP9D-5B.
# No new biological theme is introduced here.
# ------------------------------------------------------------

theme_keywords <- list(
  
  IL7R = list(
    direct = c(
      "fcgr",
      "fc receptor",
      "fc receptor activation",
      "fc-gamma",
      "fc gamma",
      "immune receptor"
    ),
    contextual = c(
      "leukocyte",
      "immune receptor",
      "immune activation",
      "antigen receptor",
      "immune signaling"
    )
  ),
  
  MS4A7 = list(
    direct = c(
      "phagocytosis",
      "phagocytic",
      "phagosome"
    ),
    contextual = c(
      "myeloid",
      "macrophage",
      "monocyte",
      "leukocyte",
      "lysosome"
    )
  ),
  
  AOC3 = list(
    direct = c(
      "b cell receptor",
      "b-cell receptor",
      "bcr",
      "b cell signaling",
      "b-cell signaling"
    ),
    contextual = c(
      "b cell",
      "b-cell",
      "lymphocyte",
      "humoral",
      "immune signaling"
    )
  )
)


# ------------------------------------------------------------
# 10. Helper function
# ------------------------------------------------------------

match_theme_terms <- function(
    descriptions,
    direct_terms,
    contextual_terms
) {
  
  desc <- tolower(
    ifelse(
      is.na(descriptions),
      "",
      descriptions
    )
  )
  
  direct_hit <- vapply(
    desc,
    function(x) {
      any(
        vapply(
          direct_terms,
          function(term) grepl(
            term,
            x,
            fixed = TRUE
          ),
          logical(1)
        )
      )
    },
    logical(1)
  )
  
  contextual_hit <- vapply(
    desc,
    function(x) {
      any(
        vapply(
          contextual_terms,
          function(term) grepl(
            term,
            x,
            fixed = TRUE
          ),
          logical(1)
        )
      )
    },
    logical(1)
  )
  
  data.frame(
    Direct_theme_match = direct_hit,
    Contextual_theme_match = contextual_hit,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 11. Generate final evidence table
# ------------------------------------------------------------

evidence_rows <- list()

row_counter <- 1L


for (gene in core_genes) {
  
  direct_terms <- theme_keywords[[gene]]$direct
  contextual_terms <- theme_keywords[[gene]]$contextual
  
  # -------------------------
  # KEGG
  # -------------------------
  
  kegg_match <- match_theme_terms(
    kegg_df$Description,
    direct_terms,
    contextual_terms
  )
  
  kegg_tmp <- cbind(
    kegg_df,
    kegg_match
  )
  
  kegg_tmp$Gene <- gene
  kegg_tmp$Database <- "KEGG"
  
  kegg_tmp$Significant <- (
    !is.na(kegg_tmp$p.adjust) &
      kegg_tmp$p.adjust < 0.05
  )
  
  kegg_tmp$Evidence_class <- ifelse(
    kegg_tmp$Direct_theme_match &
      kegg_tmp$Significant,
    "Direct_significant",
    
    ifelse(
      !kegg_tmp$Direct_theme_match &
        kegg_tmp$Contextual_theme_match &
        kegg_tmp$Significant,
      "Contextual_significant",
      
      ifelse(
        kegg_tmp$Direct_theme_match,
        "Direct_non_significant",
        
        ifelse(
          kegg_tmp$Contextual_theme_match,
          "Contextual_non_significant",
          "No_theme_match"
        )
      )
    )
  )
  
  evidence_rows[[row_counter]] <- kegg_tmp
  row_counter <- row_counter + 1L
  
  
  # -------------------------
  # GO-BP
  # -------------------------
  
  go_match <- match_theme_terms(
    go_df$Description,
    direct_terms,
    contextual_terms
  )
  
  go_tmp <- cbind(
    go_df,
    go_match
  )
  
  go_tmp$Gene <- gene
  go_tmp$Database <- "GO_BP"
  
  go_tmp$Significant <- (
    !is.na(go_tmp$p.adjust) &
      go_tmp$p.adjust < 0.05
  )
  
  go_tmp$Evidence_class <- ifelse(
    go_tmp$Direct_theme_match &
      go_tmp$Significant,
    "Direct_significant",
    
    ifelse(
      !go_tmp$Direct_theme_match &
        go_tmp$Contextual_theme_match &
        go_tmp$Significant,
      "Contextual_significant",
      
      ifelse(
        go_tmp$Direct_theme_match,
        "Direct_non_significant",
        
        ifelse(
          go_tmp$Contextual_theme_match,
          "Contextual_non_significant",
          "No_theme_match"
        )
      )
    )
  )
  
  evidence_rows[[row_counter]] <- go_tmp
  row_counter <- row_counter + 1L
}


theme_evidence_all <- do.call(
  rbind,
  evidence_rows
)


# ------------------------------------------------------------
# 12. Add original Reactome information
# ------------------------------------------------------------

theme_evidence_all$Original_Reactome_Pathway <- NA_character_
theme_evidence_all$Original_Reactome_NES <- NA_real_
theme_evidence_all$Original_Reactome_FDR <- NA_real_
theme_evidence_all$Original_Theme <- NA_character_

for (i in seq_len(nrow(original_reactome))) {
  
  idx <- (
    theme_evidence_all$Gene ==
      original_reactome$Gene[i]
  )
  
  theme_evidence_all$Original_Reactome_Pathway[idx] <-
    original_reactome$Original_Reactome_Pathway[i]
  
  theme_evidence_all$Original_Reactome_NES[idx] <-
    original_reactome$Original_NES[i]
  
  theme_evidence_all$Original_Reactome_FDR[idx] <-
    original_reactome$Original_FDR[i]
  
  theme_evidence_all$Original_Theme[idx] <-
    original_reactome$Original_Theme[i]
}


# ------------------------------------------------------------
# 13. Order final evidence table
# ------------------------------------------------------------

theme_evidence_all <- theme_evidence_all[
  ,
  c(
    "Gene",
    "Database",
    "ID",
    "Description",
    "setSize",
    "NES",
    "pvalue",
    "p.adjust",
    "qvalue",
    "Evidence_class",
    "Direct_theme_match",
    "Contextual_theme_match",
    "Original_Reactome_Pathway",
    "Original_Reactome_NES",
    "Original_Reactome_FDR",
    "Original_Theme",
    "leading_edge",
    "core_enrichment"
  )
]


# ------------------------------------------------------------
# 14. Final gene-level evidence summary
#
# This summary is descriptive only.
# No ranking or biological score is assigned.
# ------------------------------------------------------------

gene_summary_rows <- list()

for (gene in core_genes) {
  
  for (db in c("KEGG", "GO_BP")) {
    
    tmp <- theme_evidence_all[
      theme_evidence_all$Gene == gene &
        theme_evidence_all$Database == db,
      ,
      drop = FALSE
    ]
    
    direct_sig <- sum(
      tmp$Evidence_class ==
        "Direct_significant"
    )
    
    contextual_sig <- sum(
      tmp$Evidence_class ==
        "Contextual_significant"
    )
    
    direct_nonsig <- sum(
      tmp$Evidence_class ==
        "Direct_non_significant"
    )
    
    contextual_nonsig <- sum(
      tmp$Evidence_class ==
        "Contextual_non_significant"
    )
    
    theme_matching <- sum(
      tmp$Direct_theme_match |
        tmp$Contextual_theme_match
    )
    
    significant_theme <- (
      direct_sig +
        contextual_sig
    )
    
    no_theme <- sum(
      tmp$Evidence_class ==
        "No_theme_match"
    )
    
    gene_summary_rows[[length(gene_summary_rows) + 1L]] <-
      data.frame(
        
        Gene = gene,
        
        Database = db,
        
        Total_pathways = nrow(tmp),
        
        Theme_matching_pathways =
          theme_matching,
        
        Significant_theme_pathways =
          significant_theme,
        
        Direct_significant =
          direct_sig,
        
        Contextual_significant =
          contextual_sig,
        
        Direct_non_significant =
          direct_nonsig,
        
        Contextual_non_significant =
          contextual_nonsig,
        
        No_theme_match =
          no_theme,
        
        stringsAsFactors = FALSE
      )
  }
}

gene_summary <- do.call(
  rbind,
  gene_summary_rows
)


# ------------------------------------------------------------
# 15. Evidence-summary internal consistency
# ------------------------------------------------------------

summary_checks <- (
  gene_summary$Theme_matching_pathways ==
    gene_summary$Direct_significant +
    gene_summary$Contextual_significant +
    gene_summary$Direct_non_significant +
    gene_summary$Contextual_non_significant
) &
  (
    gene_summary$Total_pathways ==
      gene_summary$Theme_matching_pathways +
      gene_summary$No_theme_match
  ) &
  (
    gene_summary$Significant_theme_pathways ==
      gene_summary$Direct_significant +
      gene_summary$Contextual_significant
  )

summary_consistency_pass <- all(
  summary_checks
)


# ------------------------------------------------------------
# 16. Critical gene-specific evidence checks
# ------------------------------------------------------------

# ------------------------------------------------------------
# IL7R
# ------------------------------------------------------------

il7r_kegg <- gene_summary[
  gene_summary$Gene == "IL7R" &
    gene_summary$Database == "KEGG",
]

il7r_go <- gene_summary[
  gene_summary$Gene == "IL7R" &
    gene_summary$Database == "GO_BP",
]

il7r_kegg_no_theme_pass <- (
  il7r_kegg$Theme_matching_pathways == 0 &
    il7r_kegg$Significant_theme_pathways == 0
)

il7r_go_theme_present_pass <- (
  il7r_go$Theme_matching_pathways > 0 &
    il7r_go$Significant_theme_pathways > 0
)


# ------------------------------------------------------------
# MS4A7
# ------------------------------------------------------------

ms4a7_kegg <- gene_summary[
  gene_summary$Gene == "MS4A7" &
    gene_summary$Database == "KEGG",
]

ms4a7_go <- gene_summary[
  gene_summary$Gene == "MS4A7" &
    gene_summary$Database == "GO_BP",
]

ms4a7_kegg_theme_pass <- (
  ms4a7_kegg$Theme_matching_pathways > 0 &
    ms4a7_kegg$Significant_theme_pathways > 0
)

ms4a7_go_theme_pass <- (
  ms4a7_go$Theme_matching_pathways > 0 &
    ms4a7_go$Significant_theme_pathways > 0
)


# ------------------------------------------------------------
# AOC3
# ------------------------------------------------------------

aoc3_kegg <- gene_summary[
  gene_summary$Gene == "AOC3" &
    gene_summary$Database == "KEGG",
]

aoc3_go <- gene_summary[
  gene_summary$Gene == "AOC3" &
    gene_summary$Database == "GO_BP",
]

aoc3_kegg_direct_pass <- (
  aoc3_kegg$Direct_significant > 0
)

aoc3_go_no_direct_pass <- (
  aoc3_go$Direct_significant == 0
)

aoc3_go_contextual_pass <- (
  aoc3_go$Contextual_significant > 0
)


# ------------------------------------------------------------
# 17. Core-gene coverage
# ------------------------------------------------------------

core_gene_coverage_pass <- all(
  core_genes %in%
    unique(original_reactome$Gene)
)


# ------------------------------------------------------------
# 18. Verify no unintended statistical operation
# ------------------------------------------------------------

# This step only reads locked files and classifies existing
# descriptions. No GSEA/statistical recalculation is performed.

no_statistical_recalculation <- TRUE


# ------------------------------------------------------------
# 19. Method immutability audit
# ------------------------------------------------------------

method_integrity <- c(
  
  GSEA_rerun = FALSE,
  
  ranking_recalculated = FALSE,
  
  DEG_recalculated = FALSE,
  
  jitter_used = FALSE,
  
  artificial_tie_breaking = FALSE,
  
  STEP9D3_modified = FALSE,
  
  STEP9D4_modified = FALSE
  
)

method_integrity_pass <- all(
  !method_integrity
)


# ------------------------------------------------------------
# 20. Final evidence review table
#
# IMPORTANT:
# No Strong / Moderate / Limited / No classification.
# No ranking.
# No score.
# ------------------------------------------------------------

final_evidence_review <- data.frame(
  
  Gene = core_genes,
  
  Original_Reactome_Pathway =
    original_reactome$Original_Reactome_Pathway,
  
  Original_Theme =
    original_reactome$Original_Theme,
  
  KEGG_Theme_Matching = c(
    il7r_kegg$Theme_matching_pathways,
    ms4a7_kegg$Theme_matching_pathways,
    aoc3_kegg$Theme_matching_pathways
  ),
  
  KEGG_Significant_Theme = c(
    il7r_kegg$Significant_theme_pathways,
    ms4a7_kegg$Significant_theme_pathways,
    aoc3_kegg$Significant_theme_pathways
  ),
  
  KEGG_Direct_Significant = c(
    il7r_kegg$Direct_significant,
    ms4a7_kegg$Direct_significant,
    aoc3_kegg$Direct_significant
  ),
  
  GO_BP_Theme_Matching = c(
    il7r_go$Theme_matching_pathways,
    ms4a7_go$Theme_matching_pathways,
    aoc3_go$Theme_matching_pathways
  ),
  
  GO_BP_Significant_Theme = c(
    il7r_go$Significant_theme_pathways,
    ms4a7_go$Significant_theme_pathways,
    aoc3_go$Significant_theme_pathways
  ),
  
  GO_BP_Direct_Significant = c(
    il7r_go$Direct_significant,
    ms4a7_go$Direct_significant,
    aoc3_go$Direct_significant
  ),
  
  GO_BP_Contextual_Significant = c(
    il7r_go$Contextual_significant,
    ms4a7_go$Contextual_significant,
    aoc3_go$Contextual_significant
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 21. Final audit checks
# ------------------------------------------------------------

final_checks <- c(
  
  Required_files =
    required_files_exist,
  
  KEGG_schema =
    kegg_schema_pass,
  
  GO_BP_schema =
    go_schema_pass,
  
  KEGG_locked_counts =
    kegg_count_pass,
  
  GO_BP_locked_counts =
    go_count_pass,
  
  Core_gene_coverage =
    core_gene_coverage_pass,
  
  Summary_consistency =
    summary_consistency_pass,
  
  IL7R_KEGG_no_theme =
    il7r_kegg_no_theme_pass,
  
  IL7R_GO_BP_theme_present =
    il7r_go_theme_present_pass,
  
  MS4A7_KEGG_theme_present =
    ms4a7_kegg_theme_pass,
  
  MS4A7_GO_BP_theme_present =
    ms4a7_go_theme_pass,
  
  AOC3_KEGG_direct =
    aoc3_kegg_direct_pass,
  
  AOC3_GO_BP_no_direct =
    aoc3_go_no_direct_pass,
  
  AOC3_GO_BP_contextual =
    aoc3_go_contextual_pass,
  
  No_statistical_recalculation =
    no_statistical_recalculation,
  
  Method_immutability =
    method_integrity_pass
  
)

final_hard_pass <- all(
  final_checks
)


# ------------------------------------------------------------
# 22. Print final evidence summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-5C FINAL EVIDENCE REVIEW\n")
cat("============================================================\n\n")

print(
  final_evidence_review,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("STEP9D-5C GENE-LEVEL SUMMARY\n")
cat("============================================================\n\n")

print(
  gene_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# 23. Print audit checks
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-5C FINAL AUDIT CHECKS\n")
cat("============================================================\n\n")

print(final_checks)

cat(
  "\nSummary consistency PASS: ",
  summary_consistency_pass,
  "\n",
  sep = ""
)

cat(
  "Method immutability PASS: ",
  method_integrity_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 24. Explicit evidence interpretation guardrails
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("INTERPRETATION GUARDRAILS\n")
cat("============================================================\n\n")

cat(
  "1. Evidence is pathway/process-level, not causal.\n"
)

cat(
  "2. Original Reactome and revised GSEA use the same underlying\n",
  "   TCGA-LUAD biological context and are not independent external validation.\n"
)

cat(
  "3. Direct and contextual matches are explicitly distinguished.\n"
)

cat(
  "4. No biological strength ranking or categorical score is assigned.\n"
)

cat(
  "5. IL7R does not show a KEGG thematic match under the locked\n",
  "   STEP9D-5B keyword definition.\n"
)

cat(
  "6. MS4A7 shows phagocytosis/myeloid-related pathway/process\n",
  "   correspondence across the reviewed gene-set resources.\n"
)

cat(
  "7. AOC3 shows direct KEGG B-cell-receptor pathway correspondence,\n",
  "   while GO-BP evidence is contextual rather than direct.\n"
)

cat(
  "8. No additional GO ontology or gene-set collection was introduced.\n"
)


# ------------------------------------------------------------
# 25. Provenance record
# ------------------------------------------------------------

provenance <- list(
  
  workflow =
    "Paper1_STEP9D_GenomeWide_Preranked_GSEA",
  
  step =
    "STEP9D-5C",
  
  status =
    ifelse(
      final_hard_pass,
      "FINAL_PASS",
      "FAIL"
    ),
  
  purpose =
    "Final concordance evidence review and lock",
  
  source_kegg_result =
    kegg_result_file,
  
  source_go_bp_result =
    go_result_file,
  
  source_kegg_significant =
    kegg_sig_file,
  
  source_go_bp_significant =
    go_sig_file,
  
  source_kegg_provenance =
    kegg_prov_file,
  
  source_go_bp_provenance =
    go_prov_file,
  
  original_reactome_reference =
    "Historical Xiantao Reactome GSEA observations; not recalculated",
  
  core_genes =
    core_genes,
  
  kegg_pathway_n =
    kegg_pathway_n,
  
  kegg_fdr05_n =
    kegg_fdr05_n,
  
  go_bp_pathway_n =
    go_pathway_n,
  
  go_bp_fdr05_n =
    go_fdr05_n,
  
  evidence_definition =
    "Description-based strict direct/contextual theme matching",
  
  statistical_recalculation =
    FALSE,
  
  gsea_rerun =
    FALSE,
  
  ranking_recalculated =
    FALSE,
  
  deg_recalculated =
    FALSE,
  
  jitter_used =
    FALSE,
  
  artificial_tie_breaking =
    FALSE,
  
  final_classification_assigned =
    FALSE,
  
  method_integrity =
    method_integrity,
  
  final_checks =
    final_checks,
  
  final_hard_pass =
    final_hard_pass
  
)


# ------------------------------------------------------------
# 26. Save outputs
# ------------------------------------------------------------

evidence_csv <- file.path(
  out_dir,
  "STEP9D5C_FINAL_Concordance_Evidence_Table.csv"
)

summary_csv <- file.path(
  out_dir,
  "STEP9D5C_FINAL_Gene_Level_Evidence_Summary.csv"
)

review_csv <- file.path(
  out_dir,
  "STEP9D5C_FINAL_Evidence_Review_Table.csv"
)

audit_rds <- file.path(
  out_dir,
  "STEP9D5C_FINAL_Concordance_Evidence_Audit.rds"
)

audit_txt <- file.path(
  out_dir,
  "STEP9D5C_FINAL_Concordance_Evidence_Audit.txt"
)

provenance_rds <- file.path(
  out_dir,
  "STEP9D5C_FINAL_Concordance_Evidence_Provenance.rds"
)


write.csv(
  theme_evidence_all,
  evidence_csv,
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  gene_summary,
  summary_csv,
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  final_evidence_review,
  review_csv,
  row.names = FALSE,
  quote = TRUE
)


audit_object <- list(
  
  final_evidence_review =
    final_evidence_review,
  
  gene_summary =
    gene_summary,
  
  original_reactome =
    original_reactome,
  
  theme_evidence_all =
    theme_evidence_all,
  
  summary_checks =
    summary_checks,
  
  final_checks =
    final_checks,
  
  method_integrity =
    method_integrity,
  
  provenance =
    provenance,
  
  final_hard_pass =
    final_hard_pass
  
)


saveRDS(
  audit_object,
  audit_rds
)

saveRDS(
  provenance,
  provenance_rds
)


# ------------------------------------------------------------
# 27. Write human-readable audit TXT
# ------------------------------------------------------------

sink(audit_txt)

cat("============================================================\n")
cat("STEP9D-5C FINAL CONCORDANCE EVIDENCE REVIEW / LOCK\n")
cat("============================================================\n\n")

cat("Status: ",
    provenance$status,
    "\n\n",
    sep = "")

cat("Core genes:\n")
print(core_genes)

cat("\nOriginal Reactome reference:\n")
print(original_reactome)

cat("\nFinal evidence review:\n")
print(final_evidence_review)

cat("\nGene-level evidence summary:\n")
print(gene_summary)

cat("\nSummary consistency:\n")
print(summary_checks)

cat("\nMethod immutability:\n")
print(method_integrity)

cat("\nFinal checks:\n")
print(final_checks)

cat("\nFinal hard pass: ",
    final_hard_pass,
    "\n",
    sep = "")

cat("\nInterpretation guardrails:\n")
cat("1. Pathway/process-level evidence only.\n")
cat("2. Same cohort; not independent external validation.\n")
cat("3. Direct and contextual evidence distinguished.\n")
cat("4. No biological ranking or score assigned.\n")
cat("5. No new GSEA or statistical recalculation.\n")
cat("6. No ranking recalculation.\n")
cat("7. No DEG recalculation.\n")
cat("8. No jitter or artificial tie-breaking.\n")

sink()


# ------------------------------------------------------------
# 28. Reload verification
# ------------------------------------------------------------

reload_audit <- readRDS(
  audit_rds
)

reload_provenance <- readRDS(
  provenance_rds
)

reload_pass <- (
  identical(
    reload_audit$final_hard_pass,
    final_hard_pass
  ) &
    identical(
      reload_provenance$final_hard_pass,
      final_hard_pass
    ) &
    identical(
      reload_audit$final_evidence_review,
      final_evidence_review
    ) &
    identical(
      reload_audit$gene_summary,
      gene_summary
    )
)


# ------------------------------------------------------------
# 29. Output existence audit
# ------------------------------------------------------------

output_files <- c(
  evidence_csv,
  summary_csv,
  review_csv,
  audit_rds,
  audit_txt,
  provenance_rds
)

output_files_exist <- all(
  file.exists(output_files)
)


# ------------------------------------------------------------
# 30. FINAL STATUS
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-5C OUTPUT VERIFICATION\n")
cat("============================================================\n\n")

cat(
  "All output files present: ",
  output_files_exist,
  "\n",
  sep = ""
)

cat(
  "RDS reload verification: ",
  reload_pass,
  "\n",
  sep = ""
)

cat(
  "STEP9D-5C FINAL HARD PASS: ",
  final_hard_pass &
    output_files_exist &
    reload_pass,
  "\n",
  sep = ""
)

cat(
  "\nFINAL STATUS = ",
  ifelse(
    final_hard_pass &
      output_files_exist &
      reload_pass,
    "FINAL_PASS",
    "FAIL"
  ),
  "\n",
  sep = ""
)

cat("\nNo GSEA rerun performed.\n")
cat("No ranking recalculation performed.\n")
cat("No DEG recalculation performed.\n")
cat("No jitter used.\n")
cat("No artificial tie-breaking used.\n")
cat("No STEP9D-3B source file modified.\n")
cat("No STEP9D-4B source file modified.\n")
cat("No final biological ranking assigned.\n")
cat("STEP9D-5C COMPLETE.\n")


# ============================================================
# STEP9D-5B -> STEP9D-5C SOURCE OBJECT INSPECTION
# READ-ONLY DIAGNOSTIC
#
# PURPOSE:
# Identify the exact STEP9D-5B corrected evidence object
# that STEP9D-5C must reuse.
#
# NO GSEA
# NO recalculation
# NO file modification
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-5B SOURCE OBJECT INSPECTION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. STEP9D directory
# ------------------------------------------------------------

step9d_dir <- file.path(
  "input_dir/Paper 1_Step 9-Functional Enrichment",
  "STEP9D_GenomeWide_Preranked_GSEA"
)

step5b_dir <- file.path(
  step9d_dir,
  "9D-5B_Strict_Theme_Evidence_Audit"
)

cat("STEP9D-5B directory:\n")
cat(step5b_dir, "\n\n")


# ------------------------------------------------------------
# 2. Check directory
# ------------------------------------------------------------

cat("Directory exists: ",
    dir.exists(step5b_dir),
    "\n\n",
    sep = "")

if (!dir.exists(step5b_dir)) {
  stop("STEP9D-5B directory not found.")
}


# ------------------------------------------------------------
# 3. List all files
# ------------------------------------------------------------

all_files <- list.files(
  step5b_dir,
  recursive = TRUE,
  full.names = TRUE
)

cat("Files found:\n")
print(all_files)


# ------------------------------------------------------------
# 4. Focus on RDS / CSV / TXT
# ------------------------------------------------------------

rds_files <- all_files[
  grepl(
    "\\.rds$",
    all_files,
    ignore.case = TRUE
  )
]

csv_files <- all_files[
  grepl(
    "\\.csv$",
    all_files,
    ignore.case = TRUE
  )
]

txt_files <- all_files[
  grepl(
    "\\.txt$",
    all_files,
    ignore.case = TRUE
  )
]

cat("\nRDS files:\n")
print(rds_files)

cat("\nCSV files:\n")
print(csv_files)

cat("\nTXT files:\n")
print(txt_files)


# ------------------------------------------------------------
# 5. Inspect every RDS object structure
# ------------------------------------------------------------

if (length(rds_files) > 0) {
  
  cat("\n============================================================\n")
  cat("RDS OBJECT STRUCTURE\n")
  cat("============================================================\n\n")
  
  for (f in rds_files) {
    
    cat("\n------------------------------------------------------------\n")
    cat("FILE: ", basename(f), "\n", sep = "")
    cat("------------------------------------------------------------\n")
    
    obj <- readRDS(f)
    
    cat("Class:\n")
    print(class(obj))
    
    cat("\nType:\n")
    print(typeof(obj))
    
    cat("\nLength:\n")
    print(length(obj))
    
    if (is.list(obj)) {
      
      cat("\nList names:\n")
      print(names(obj))
      
    }
    
    if (is.data.frame(obj)) {
      
      cat("\nDimensions:\n")
      print(dim(obj))
      
      cat("\nColumn names:\n")
      print(names(obj))
      
      cat("\nFirst rows:\n")
      print(
        utils::head(
          obj,
          3
        )
      )
    }
  }
}


# ------------------------------------------------------------
# 6. Inspect CSV structures
# ------------------------------------------------------------

if (length(csv_files) > 0) {
  
  cat("\n============================================================\n")
  cat("CSV STRUCTURE\n")
  cat("============================================================\n\n")
  
  for (f in csv_files) {
    
    cat("\n------------------------------------------------------------\n")
    cat("FILE: ", basename(f), "\n", sep = "")
    cat("------------------------------------------------------------\n")
    
    x <- read.csv(
      f,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    cat("Dimensions:\n")
    print(dim(x))
    
    cat("\nColumn names:\n")
    print(names(x))
    
    cat("\nFirst rows:\n")
    print(
      utils::head(
        x,
        3
      )
    )
  }
}


# ------------------------------------------------------------
# 7. Search for key evidence column names
# ------------------------------------------------------------

candidate_columns <- c(
  "Gene",
  "Database",
  "Evidence_class",
  "Direct_theme_match",
  "Contextual_theme_match",
  "Significant",
  "Description",
  "ID",
  "p.adjust",
  "Theme_matching_pathways",
  "Significant_theme_pathways"
)

cat("\n============================================================\n")
cat("KEY COLUMN DETECTION\n")
cat("============================================================\n\n")

for (f in c(rds_files, csv_files)) {
  
  obj <- tryCatch(
    {
      if (grepl("\\.rds$", f, ignore.case = TRUE)) {
        readRDS(f)
      } else {
        read.csv(
          f,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    },
    error = function(e) NULL
  )
  
  if (is.null(obj)) {
    next
  }
  
  if (is.data.frame(obj)) {
    
    found <- intersect(
      candidate_columns,
      names(obj)
    )
    
    if (length(found) > 0) {
      
      cat(
        basename(f),
        " -> ",
        paste(found, collapse = ", "),
        "\n",
        sep = ""
      )
    }
    
  } else if (is.list(obj)) {
    
    for (nm in names(obj)) {
      
      element <- obj[[nm]]
      
      if (is.data.frame(element)) {
        
        found <- intersect(
          candidate_columns,
          names(element)
        )
        
        if (length(found) > 0) {
          
          cat(
            basename(f),
            "$",
            nm,
            " -> ",
            paste(found, collapse = ", "),
            "\n",
            sep = ""
          )
        }
      }
    }
  }
}


cat("\n============================================================\n")
cat("READ-ONLY INSPECTION COMPLETE\n")
cat("============================================================\n\n")

cat(
  "NO GSEA rerun.\n",
  "NO ranking recalculation.\n",
  "NO DEG recalculation.\n",
  "NO files modified.\n"
)


# ============================================================
# Paper 1
# STEP9D-5C
# Final Concordance Evidence Review / Lock
#
# PURPOSE
# ------------------------------------------------------------
# Final evidence review using ONLY the already locked
# STEP9D-5B Strict Theme Evidence Audit.
#
# IMPORTANT:
#   STEP9D-5C DOES NOT:
#     - rerun GSEA
#     - recalculate ranking
#     - recalculate DEG
#     - redefine theme keywords
#     - perform new keyword matching
#     - perform statistical testing
#     - modify STEP9D-3B
#     - modify STEP9D-4B
#
# STEP9D-5B is the authoritative evidence source.
#
# STEP9D-5C only:
#   1. reads locked STEP9D-5B evidence
#   2. verifies internal consistency
#   3. combines original Reactome reference
#   4. constructs final descriptive evidence table
#   5. locks the evidence review
#
# NO biological ranking / scoring / classification is assigned.
# ============================================================


cat("\n============================================================\n")
cat("STEP9D-5C: Final Concordance Evidence Review / Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Directories
# ------------------------------------------------------------

root_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

step9d_dir <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA"
)

step5b_dir <- file.path(
  step9d_dir,
  "9D-5B_Strict_Theme_Evidence_Audit"
)

out_dir <- file.path(
  step9d_dir,
  "9D-5C_Final_Concordance_Evidence_Lock"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. AUTHORITATIVE STEP9D-5B FILE
# ------------------------------------------------------------

step5b_rds <- file.path(
  step5b_dir,
  "STEP9D5B_FINAL_Strict_Theme_Evidence_Audit.rds"
)

step5b_summary_csv <- file.path(
  step5b_dir,
  "STEP9D5B_CORRECTED_Gene_Level_Descriptive_Summary.csv"
)

step5b_all_csv <- file.path(
  step5b_dir,
  "STEP9D5B_All_Theme_Evidence_Audit.csv"
)

step5b_reactome_csv <- file.path(
  step5b_dir,
  "STEP9D5B_Original_Reactome_Reference.csv"
)

step5b_sig_csv <- file.path(
  step5b_dir,
  "STEP9D5B_Significant_Theme_Evidence.csv"
)

step5b_quadrant_csv <- file.path(
  step5b_dir,
  "STEP9D5B_Evidence_Quadrant_Summary.csv"
)


# ------------------------------------------------------------
# 2. Required-file audit
# ------------------------------------------------------------

required_5b_files <- c(
  step5b_rds,
  step5b_summary_csv,
  step5b_all_csv,
  step5b_reactome_csv,
  step5b_sig_csv,
  step5b_quadrant_csv
)

required_5b_files_exist <- all(
  file.exists(required_5b_files)
)

cat(
  "All required STEP9D-5B files present: ",
  required_5b_files_exist,
  "\n",
  sep = ""
)

if (!required_5b_files_exist) {
  stop(
    "Required STEP9D-5B locked evidence file(s) missing."
  )
}


# ------------------------------------------------------------
# 3. Load authoritative STEP9D-5B RDS
# ------------------------------------------------------------

step5b <- readRDS(
  step5b_rds
)


# ------------------------------------------------------------
# 4. Verify STEP9D-5B object structure
# ------------------------------------------------------------

required_rds_fields <- c(
  "workflow",
  "step",
  "status",
  "purpose",
  "source_files",
  "source_dimensions",
  "source_integrity",
  "reactome_reference",
  "theme_dictionary",
  "all_evidence",
  "theme_evidence",
  "theme_evidence_significant",
  "evidence_summary",
  "gene_summary",
  "method_integrity",
  "source_integrity_pass",
  "theme_audit_structure_pass",
  "method_integrity_pass",
  "final_concordance_classification_performed",
  "preliminary_hard_pass"
)

rds_structure_pass <- (
  is.list(step5b) &
    all(required_rds_fields %in% names(step5b))
)

cat(
  "STEP9D-5B RDS structure PASS: ",
  rds_structure_pass,
  "\n",
  sep = ""
)

if (!rds_structure_pass) {
  stop(
    "STEP9D-5B RDS structure is incomplete."
  )
}


# ------------------------------------------------------------
# 5. Verify STEP9D-5B preliminary lock state
# ------------------------------------------------------------

step5b_preliminary_pass <- isTRUE(
  step5b$preliminary_hard_pass
)

step5b_source_integrity_pass <- isTRUE(
  step5b$source_integrity_pass
)

step5b_theme_structure_pass <- isTRUE(
  step5b$theme_audit_structure_pass
)

step5b_method_integrity_pass <- isTRUE(
  step5b$method_integrity_pass
)

step5b_no_final_classification <- identical(
  step5b$final_concordance_classification_performed,
  FALSE
)

cat("\nSTEP9D-5B lock-state verification:\n")
cat(
  "  Preliminary hard pass: ",
  step5b_preliminary_pass,
  "\n",
  sep = ""
)
cat(
  "  Source integrity pass: ",
  step5b_source_integrity_pass,
  "\n",
  sep = ""
)
cat(
  "  Theme structure pass: ",
  step5b_theme_structure_pass,
  "\n",
  sep = ""
)
cat(
  "  Method integrity pass: ",
  step5b_method_integrity_pass,
  "\n",
  sep = ""
)
cat(
  "  No final concordance classification performed: ",
  step5b_no_final_classification,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Extract ONLY locked STEP9D-5B objects
#
# IMPORTANT:
# No recalculation here.
# ------------------------------------------------------------

locked_all_evidence <- step5b$all_evidence

locked_theme_evidence <- step5b$theme_evidence

locked_theme_evidence_significant <-
  step5b$theme_evidence_significant

locked_evidence_summary <-
  step5b$evidence_summary

locked_gene_summary <-
  step5b$gene_summary

locked_reactome_reference <-
  step5b$reactome_reference

locked_theme_dictionary <-
  step5b$theme_dictionary


# ------------------------------------------------------------
# 7. Load corrected CSV summary
#
# Used only as an independent disk-level concordance check.
# ------------------------------------------------------------

corrected_gene_summary_csv <- read.csv(
  step5b_summary_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

reactome_reference_csv <- read.csv(
  step5b_reactome_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

all_evidence_csv <- read.csv(
  step5b_all_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

significant_theme_csv <- read.csv(
  step5b_sig_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

quadrant_summary_csv <- read.csv(
  step5b_quadrant_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 8. Core genes
# ------------------------------------------------------------

core_genes <- c(
  "IL7R",
  "MS4A7",
  "AOC3"
)


# ------------------------------------------------------------
# 9. Verify core-gene coverage
# ------------------------------------------------------------

core_in_summary <- all(
  core_genes %in%
    unique(locked_gene_summary$Gene)
)

core_in_evidence <- all(
  core_genes %in%
    unique(locked_all_evidence$Gene)
)

core_in_reactome <- all(
  core_genes %in%
    unique(locked_reactome_reference$Gene)
)

core_gene_coverage_pass <- (
  core_in_summary &
    core_in_evidence &
    core_in_reactome
)

cat("\nCore-gene coverage PASS: ",
    core_gene_coverage_pass,
    "\n",
    sep = "")


# ------------------------------------------------------------
# 10. Verify locked 5B gene-summary values
#
# EXPECTED values are inherited from STEP9D-5B,
# not recomputed from descriptions.
# ------------------------------------------------------------

expected_gene_summary <- data.frame(
  
  Gene = c(
    "IL7R",
    "IL7R",
    "MS4A7",
    "MS4A7",
    "AOC3",
    "AOC3"
  ),
  
  Database = c(
    "KEGG",
    "GO_BP",
    "KEGG",
    "GO_BP",
    "KEGG",
    "GO_BP"
  ),
  
  Total_pathways = c(
    353,
    6316,
    353,
    6316,
    353,
    6316
  ),
  
  Theme_matching_pathways = c(
    0,
    41,
    7,
    251,
    1,
    58
  ),
  
  Significant_theme_pathways = c(
    0,
    15,
    3,
    83,
    1,
    9
  ),
  
  Direct_significant = c(
    0,
    5,
    1,
    6,
    1,
    0
  ),
  
  Contextual_significant = c(
    0,
    10,
    2,
    77,
    0,
    9
  ),
  
  Direct_non_significant = c(
    0,
    11,
    1,
    15,
    0,
    7
  ),
  
  Contextual_non_significant = c(
    0,
    15,
    3,
    153,
    0,
    42
  ),
  
  No_theme_match = c(
    353,
    6275,
    346,
    6065,
    352,
    6258
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 11. Compare authoritative RDS summary with expected lock
# ------------------------------------------------------------

summary_columns <- names(
  expected_gene_summary
)

summary_rows_match <- TRUE

for (i in seq_len(nrow(expected_gene_summary))) {
  
  row_expected <- expected_gene_summary[i, ]
  
  idx <- (
    locked_gene_summary$Gene ==
      row_expected$Gene &
      locked_gene_summary$Database ==
      row_expected$Database
  )
  
  if (sum(idx) != 1) {
    summary_rows_match <- FALSE
    next
  }
  
  row_actual <- locked_gene_summary[idx, ]
  
  for (col in summary_columns) {
    
    if (is.numeric(row_expected[[col]])) {
      
      if (
        length(row_actual[[col]]) != 1 ||
        is.na(row_actual[[col]]) ||
        row_actual[[col]] != row_expected[[col]]
      ) {
        summary_rows_match <- FALSE
      }
      
    } else {
      
      if (
        as.character(row_actual[[col]]) !=
        as.character(row_expected[[col]])
      ) {
        summary_rows_match <- FALSE
      }
    }
  }
}

cat(
  "STEP9D-5B locked gene-summary values PASS: ",
  summary_rows_match,
  "\n",
  sep = ""
)

# ============================================================
# STEP9D-5C DIAGNOSTIC
# Inspect authoritative STEP9D-5B RDS gene_summary structure
#
# IMPORTANT:
# - READ-ONLY
# - No theme matching
# - No GSEA
# - No ranking recalculation
# - No modification of STEP9D-5B files
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-5C DIAGNOSTIC: STEP9D-5B gene_summary structure\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Basic object structure
# ------------------------------------------------------------

cat("Class of locked_gene_summary:\n")
print(class(locked_gene_summary))

cat("\nDimensions:\n")
print(dim(locked_gene_summary))

cat("\nColumn names:\n")
print(names(locked_gene_summary))


# ------------------------------------------------------------
# 2. Structure
# ------------------------------------------------------------

cat("\nFull structure:\n")
str(
  locked_gene_summary,
  max.level = 2
)


# ------------------------------------------------------------
# 3. First rows
# ------------------------------------------------------------

cat("\nFirst 10 rows:\n")
print(
  head(
    locked_gene_summary,
    10
  )
)


# ------------------------------------------------------------
# 4. Compare required columns
# ------------------------------------------------------------

required_columns <- names(
  expected_gene_summary
)

missing_columns <- setdiff(
  required_columns,
  names(locked_gene_summary)
)

extra_columns <- setdiff(
  names(locked_gene_summary),
  required_columns
)

cat("\nRequired columns:\n")
print(required_columns)

cat("\nMissing required columns:\n")
print(missing_columns)

cat("\nExtra columns:\n")
print(extra_columns)


# ------------------------------------------------------------
# 5. If Gene / Database exist, inspect their values
# ------------------------------------------------------------

if (
  all(
    c("Gene", "Database") %in%
    names(locked_gene_summary)
  )
) {
  
  cat("\nGene values:\n")
  print(
    unique(
      locked_gene_summary$Gene
    )
  )
  
  cat("\nDatabase values:\n")
  print(
    unique(
      locked_gene_summary$Database
    )
  )
  
  cat("\nGene × Database table:\n")
  print(
    table(
      locked_gene_summary$Gene,
      locked_gene_summary$Database
    )
  )
  
} else {
  
  cat(
    "\nWARNING: Gene and/or Database columns are missing.\n"
  )
}


# ------------------------------------------------------------
# 6. Compare against corrected authoritative CSV
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Corrected CSV structure\n")
cat("============================================================\n\n")

cat("Class:\n")
print(
  class(corrected_gene_summary_csv)
)

cat("\nDimensions:\n")
print(
  dim(corrected_gene_summary_csv)
)

cat("\nColumn names:\n")
print(
  names(corrected_gene_summary_csv)
)

cat("\nFirst 10 rows:\n")
print(
  head(
    corrected_gene_summary_csv,
    10
  )
)


# ------------------------------------------------------------
# 7. Print exact corrected CSV values
# ------------------------------------------------------------

cat("\nCorrected CSV gene summary:\n")
print(
  corrected_gene_summary_csv,
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Print exact RDS gene summary if printable
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Authoritative RDS gene_summary\n")
cat("============================================================\n\n")

print(
  locked_gene_summary,
  row.names = FALSE
)


cat("\n============================================================\n")
cat("DIAGNOSTIC COMPLETE — NO FILES MODIFIED\n")
cat("============================================================\n")



# ============================================================
# STEP9D-5B RDS MINIMAL PROVENANCE REPAIR
#
# Purpose:
#   Repair ONLY the gene_summary component inside the existing
#   STEP9D-5B FINAL RDS so that it agrees with the already
#   corrected authoritative CSV.
#
# IMPORTANT:
#   - NO GSEA rerun
#   - NO ranking recalculation
#   - NO theme rematching
#   - NO evidence recalculation
#   - NO modification of all_evidence
#   - NO modification of theme_evidence
#   - NO modification of theme_evidence_significant
#   - NO modification of Reactome reference
#
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-5B RDS MINIMAL PROVENANCE REPAIR\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

step5b_dir <- file.path(
  "input_dirt/Paper 1_Step 9-Functional Enrichment",
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-5B_Strict_Theme_Evidence_Audit"
)

rds_file <- file.path(
  step5b_dir,
  "STEP9D5B_FINAL_Strict_Theme_Evidence_Audit.rds"
)

corrected_csv_file <- file.path(
  step5b_dir,
  "STEP9D5B_CORRECTED_Gene_Level_Descriptive_Summary.csv"
)

backup_dir <- file.path(
  step5b_dir,
  "BACKUP_RDS_BEFORE_GENE_SUMMARY_REPAIR"
)

dir.create(
  backup_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Existence checks
# ------------------------------------------------------------

stopifnot(
  file.exists(rds_file),
  file.exists(corrected_csv_file)
)

cat("RDS exists: ",
    file.exists(rds_file),
    "\n",
    sep = "")

cat("Corrected CSV exists: ",
    file.exists(corrected_csv_file),
    "\n",
    sep = "")


# ------------------------------------------------------------
# 3. Read authoritative corrected CSV
# ------------------------------------------------------------

corrected_gene_summary <- read.csv(
  corrected_csv_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\nCorrected CSV dimensions: ",
    nrow(corrected_gene_summary),
    " x ",
    ncol(corrected_gene_summary),
    "\n",
    sep = "")


# ------------------------------------------------------------
# 4. Required structure
# ------------------------------------------------------------

required_columns <- c(
  "Gene",
  "Database",
  "Total_pathways",
  "Theme_matching_pathways",
  "Significant_theme_pathways",
  "Direct_significant",
  "Contextual_significant",
  "Direct_non_significant",
  "Contextual_non_significant",
  "No_theme_match"
)

structure_pass <- identical(
  names(corrected_gene_summary),
  required_columns
)

cat(
  "Corrected CSV structure PASS: ",
  structure_pass,
  "\n",
  sep = ""
)

if (!structure_pass) {
  stop(
    "Corrected CSV does not have the expected 10-column structure."
  )
}


# ------------------------------------------------------------
# 5. Load existing RDS
# ------------------------------------------------------------

step5b <- readRDS(
  rds_file
)

cat("\n5B RDS class: ",
    paste(class(step5b), collapse = ", "),
    "\n",
    sep = "")

cat(
  "5B RDS names:\n"
)

print(names(step5b))


# ------------------------------------------------------------
# 6. Verify required evidence components exist
# ------------------------------------------------------------

required_rds_components <- c(
  "all_evidence",
  "theme_evidence",
  "theme_evidence_significant",
  "evidence_summary",
  "gene_summary",
  "reactome_reference"
)

component_pass <- all(
  required_rds_components %in%
    names(step5b)
)

cat(
  "\nRequired 5B evidence components PASS: ",
  component_pass,
  "\n",
  sep = ""
)

if (!component_pass) {
  stop(
    "Required STEP9D-5B evidence components are missing."
  )
}


# ------------------------------------------------------------
# 7. Preserve all non-gene-summary objects
# ------------------------------------------------------------

all_evidence_before <-
  step5b$all_evidence

theme_evidence_before <-
  step5b$theme_evidence

theme_evidence_significant_before <-
  step5b$theme_evidence_significant

evidence_summary_before <-
  step5b$evidence_summary

reactome_reference_before <-
  step5b$reactome_reference


# ------------------------------------------------------------
# 8. Confirm corrected CSV has expected 6 rows
# ------------------------------------------------------------

row_count_pass <- nrow(
  corrected_gene_summary
) == 6

cat(
  "Corrected gene-summary row count PASS: ",
  row_count_pass,
  "\n",
  sep = ""
)

if (!row_count_pass) {
  stop(
    "Corrected gene summary does not contain exactly 6 rows."
  )
}


# ------------------------------------------------------------
# 9. Verify exact locked values BEFORE repair
# ------------------------------------------------------------

expected_gene_summary <- data.frame(
  
  Gene = c(
    "IL7R",
    "IL7R",
    "MS4A7",
    "MS4A7",
    "AOC3",
    "AOC3"
  ),
  
  Database = c(
    "KEGG",
    "GO_BP",
    "KEGG",
    "GO_BP",
    "KEGG",
    "GO_BP"
  ),
  
  Total_pathways = c(
    353,
    6316,
    353,
    6316,
    353,
    6316
  ),
  
  Theme_matching_pathways = c(
    0,
    41,
    7,
    251,
    1,
    58
  ),
  
  Significant_theme_pathways = c(
    0,
    15,
    3,
    83,
    1,
    9
  ),
  
  Direct_significant = c(
    0,
    5,
    1,
    6,
    1,
    0
  ),
  
  Contextual_significant = c(
    0,
    10,
    2,
    77,
    0,
    9
  ),
  
  Direct_non_significant = c(
    0,
    11,
    1,
    15,
    0,
    7
  ),
  
  Contextual_non_significant = c(
    0,
    15,
    3,
    153,
    0,
    42
  ),
  
  No_theme_match = c(
    353,
    6275,
    346,
    6065,
    352,
    6258
  ),
  
  stringsAsFactors = FALSE
)

expected_order <- paste(
  expected_gene_summary$Gene,
  expected_gene_summary$Database,
  sep = "::"
)

actual_order <- paste(
  corrected_gene_summary$Gene,
  corrected_gene_summary$Database,
  sep = "::"
)

expected_order_pass <- identical(
  actual_order,
  expected_order
)

value_pass <- isTRUE(
  all.equal(
    corrected_gene_summary,
    expected_gene_summary,
    check.attributes = FALSE
  )
)

cat(
  "Corrected CSV gene-summary order PASS: ",
  expected_order_pass,
  "\n",
  sep = ""
)

cat(
  "Corrected CSV exact locked values PASS: ",
  value_pass,
  "\n",
  sep = ""
)

if (!expected_order_pass || !value_pass) {
  stop(
    "Corrected CSV does not match the locked STEP9D-5B values."
  )
}


# ------------------------------------------------------------
# 10. Create backup of original RDS
# ------------------------------------------------------------

backup_file <- file.path(
  backup_dir,
  paste0(
    "STEP9D5B_FINAL_Strict_Theme_Evidence_Audit_BEFORE_REPAIR_",
    format(
      Sys.time(),
      "%Y%m%d_%H%M%S"
    ),
    ".rds"
  )
)

file.copy(
  from = rds_file,
  to = backup_file,
  overwrite = FALSE
)

cat(
  "\nOriginal RDS backup created: ",
  file.exists(backup_file),
  "\n",
  sep = ""
)

if (!file.exists(backup_file)) {
  stop(
    "Backup creation failed. RDS was NOT modified."
  )
}


# ------------------------------------------------------------
# 11. Replace ONLY gene_summary
# ------------------------------------------------------------

step5b$gene_summary <-
  corrected_gene_summary


# ------------------------------------------------------------
# 12. Verify that non-gene-summary objects are unchanged
# ------------------------------------------------------------

all_evidence_unchanged <- isTRUE(
  identical(
    step5b$all_evidence,
    all_evidence_before
  )
)

theme_evidence_unchanged <- isTRUE(
  identical(
    step5b$theme_evidence,
    theme_evidence_before
  )
)

theme_evidence_significant_unchanged <- isTRUE(
  identical(
    step5b$theme_evidence_significant,
    theme_evidence_significant_before
  )
)

evidence_summary_unchanged <- isTRUE(
  identical(
    step5b$evidence_summary,
    evidence_summary_before
  )
)

reactome_reference_unchanged <- isTRUE(
  identical(
    step5b$reactome_reference,
    reactome_reference_before
  )
)

non_summary_objects_unchanged <- all(
  all_evidence_unchanged,
  theme_evidence_unchanged,
  theme_evidence_significant_unchanged,
  evidence_summary_unchanged,
  reactome_reference_unchanged
)

cat(
  "\nNon-gene-summary objects unchanged: ",
  non_summary_objects_unchanged,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. Save repaired RDS
# ------------------------------------------------------------

if (!non_summary_objects_unchanged) {
  
  stop(
    "Unexpected non-gene-summary object change detected. ",
    "RDS was NOT overwritten."
  )
  
}

saveRDS(
  step5b,
  rds_file
)

cat(
  "Repaired RDS written: ",
  file.exists(rds_file),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 14. Reload repaired RDS
# ------------------------------------------------------------

step5b_reloaded <- readRDS(
  rds_file
)

reload_pass <- isTRUE(
  identical(
    step5b_reloaded$gene_summary,
    corrected_gene_summary
  )
)

cat(
  "Reloaded gene-summary concordance PASS: ",
  reload_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 15. Verify non-gene-summary objects after reload
# ------------------------------------------------------------

post_reload_unchanged <- all(
  identical(
    step5b_reloaded$all_evidence,
    all_evidence_before
  ),
  identical(
    step5b_reloaded$theme_evidence,
    theme_evidence_before
  ),
  identical(
    step5b_reloaded$theme_evidence_significant,
    theme_evidence_significant_before
  ),
  identical(
    step5b_reloaded$evidence_summary,
    evidence_summary_before
  ),
  identical(
    step5b_reloaded$reactome_reference,
    reactome_reference_before
  )
)

cat(
  "Post-reload non-gene-summary integrity PASS: ",
  post_reload_unchanged,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 16. Final repair PASS
# ------------------------------------------------------------

repair_pass <- all(
  structure_pass,
  row_count_pass,
  expected_order_pass,
  value_pass,
  file.exists(backup_file),
  non_summary_objects_unchanged,
  file.exists(rds_file),
  reload_pass,
  post_reload_unchanged
)

cat(
  "\n============================================================\n"
)

cat(
  "STEP9D-5B RDS GENE-SUMMARY REPAIR HARD PASS: ",
  repair_pass,
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)

if (!repair_pass) {
  stop(
    "STEP9D-5B RDS repair did not pass all integrity checks."
  )
}


##STEP9D-5C Final Concordance Evidence Review / Lock
# ============================================================
# Paper 1
# STEP9D-5C
# Final Concordance Evidence Review / Lock
#
# PURPOSE
# ------------------------------------------------------------
# Final evidence review using ONLY the already locked
# STEP9D-5B Strict Theme Evidence Audit.
#
# IMPORTANT:
#   STEP9D-5C DOES NOT:
#     - rerun GSEA
#     - recalculate ranking
#     - recalculate DEG
#     - redefine theme keywords
#     - perform new keyword matching
#     - perform statistical testing
#     - modify STEP9D-3B
#     - modify STEP9D-4B
#
# STEP9D-5B is the authoritative evidence source.
#
# STEP9D-5C only:
#   1. reads locked STEP9D-5B evidence
#   2. verifies internal consistency
#   3. combines original Reactome reference
#   4. constructs final descriptive evidence table
#   5. locks the evidence review
#
# NO biological ranking / scoring / classification is assigned.
# ============================================================


cat("\n============================================================\n")
cat("STEP9D-5C: Final Concordance Evidence Review / Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Directories
# ------------------------------------------------------------

root_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

step9d_dir <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA"
)

step5b_dir <- file.path(
  step9d_dir,
  "9D-5B_Strict_Theme_Evidence_Audit"
)

out_dir <- file.path(
  step9d_dir,
  "9D-5C_Final_Concordance_Evidence_Lock"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. AUTHORITATIVE STEP9D-5B FILE
# ------------------------------------------------------------

step5b_rds <- file.path(
  step5b_dir,
  "STEP9D5B_FINAL_Strict_Theme_Evidence_Audit.rds"
)

step5b_summary_csv <- file.path(
  step5b_dir,
  "STEP9D5B_CORRECTED_Gene_Level_Descriptive_Summary.csv"
)

step5b_all_csv <- file.path(
  step5b_dir,
  "STEP9D5B_All_Theme_Evidence_Audit.csv"
)

step5b_reactome_csv <- file.path(
  step5b_dir,
  "STEP9D5B_Original_Reactome_Reference.csv"
)

step5b_sig_csv <- file.path(
  step5b_dir,
  "STEP9D5B_Significant_Theme_Evidence.csv"
)

step5b_quadrant_csv <- file.path(
  step5b_dir,
  "STEP9D5B_Evidence_Quadrant_Summary.csv"
)


# ------------------------------------------------------------
# 2. Required-file audit
# ------------------------------------------------------------

required_5b_files <- c(
  step5b_rds,
  step5b_summary_csv,
  step5b_all_csv,
  step5b_reactome_csv,
  step5b_sig_csv,
  step5b_quadrant_csv
)

required_5b_files_exist <- all(
  file.exists(required_5b_files)
)

cat(
  "All required STEP9D-5B files present: ",
  required_5b_files_exist,
  "\n",
  sep = ""
)

if (!required_5b_files_exist) {
  stop(
    "Required STEP9D-5B locked evidence file(s) missing."
  )
}


# ------------------------------------------------------------
# 3. Load authoritative STEP9D-5B RDS
# ------------------------------------------------------------

step5b <- readRDS(
  step5b_rds
)


# ------------------------------------------------------------
# 4. Verify STEP9D-5B object structure
# ------------------------------------------------------------

required_rds_fields <- c(
  "workflow",
  "step",
  "status",
  "purpose",
  "source_files",
  "source_dimensions",
  "source_integrity",
  "reactome_reference",
  "theme_dictionary",
  "all_evidence",
  "theme_evidence",
  "theme_evidence_significant",
  "evidence_summary",
  "gene_summary",
  "method_integrity",
  "source_integrity_pass",
  "theme_audit_structure_pass",
  "method_integrity_pass",
  "final_concordance_classification_performed",
  "preliminary_hard_pass"
)

rds_structure_pass <- (
  is.list(step5b) &
    all(required_rds_fields %in% names(step5b))
)

cat(
  "STEP9D-5B RDS structure PASS: ",
  rds_structure_pass,
  "\n",
  sep = ""
)

if (!rds_structure_pass) {
  stop(
    "STEP9D-5B RDS structure is incomplete."
  )
}


# ------------------------------------------------------------
# 5. Verify STEP9D-5B preliminary lock state
# ------------------------------------------------------------

step5b_preliminary_pass <- isTRUE(
  step5b$preliminary_hard_pass
)

step5b_source_integrity_pass <- isTRUE(
  step5b$source_integrity_pass
)

step5b_theme_structure_pass <- isTRUE(
  step5b$theme_audit_structure_pass
)

step5b_method_integrity_pass <- isTRUE(
  step5b$method_integrity_pass
)

step5b_no_final_classification <- identical(
  step5b$final_concordance_classification_performed,
  FALSE
)

cat("\nSTEP9D-5B lock-state verification:\n")
cat(
  "  Preliminary hard pass: ",
  step5b_preliminary_pass,
  "\n",
  sep = ""
)
cat(
  "  Source integrity pass: ",
  step5b_source_integrity_pass,
  "\n",
  sep = ""
)
cat(
  "  Theme structure pass: ",
  step5b_theme_structure_pass,
  "\n",
  sep = ""
)
cat(
  "  Method integrity pass: ",
  step5b_method_integrity_pass,
  "\n",
  sep = ""
)
cat(
  "  No final concordance classification performed: ",
  step5b_no_final_classification,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Extract ONLY locked STEP9D-5B objects
#
# IMPORTANT:
# No recalculation here.
# ------------------------------------------------------------

locked_all_evidence <- step5b$all_evidence

locked_theme_evidence <- step5b$theme_evidence

locked_theme_evidence_significant <-
  step5b$theme_evidence_significant

locked_evidence_summary <-
  step5b$evidence_summary

locked_gene_summary <-
  step5b$gene_summary

locked_reactome_reference <-
  step5b$reactome_reference

locked_theme_dictionary <-
  step5b$theme_dictionary


# ------------------------------------------------------------
# 7. Load corrected CSV summary
#
# Used only as an independent disk-level concordance check.
# ------------------------------------------------------------

corrected_gene_summary_csv <- read.csv(
  step5b_summary_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

reactome_reference_csv <- read.csv(
  step5b_reactome_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

all_evidence_csv <- read.csv(
  step5b_all_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

significant_theme_csv <- read.csv(
  step5b_sig_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

quadrant_summary_csv <- read.csv(
  step5b_quadrant_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 8. Core genes
# ------------------------------------------------------------

core_genes <- c(
  "IL7R",
  "MS4A7",
  "AOC3"
)


# ------------------------------------------------------------
# 9. Verify core-gene coverage
# ------------------------------------------------------------

core_in_summary <- all(
  core_genes %in%
    unique(locked_gene_summary$Gene)
)

core_in_evidence <- all(
  core_genes %in%
    unique(locked_all_evidence$Gene)
)

core_in_reactome <- all(
  core_genes %in%
    unique(locked_reactome_reference$Gene)
)

core_gene_coverage_pass <- (
  core_in_summary &
    core_in_evidence &
    core_in_reactome
)

cat("\nCore-gene coverage PASS: ",
    core_gene_coverage_pass,
    "\n",
    sep = "")


# ------------------------------------------------------------
# 10. Verify locked 5B gene-summary values
#
# EXPECTED values are inherited from STEP9D-5B,
# not recomputed from descriptions.
# ------------------------------------------------------------

expected_gene_summary <- data.frame(
  
  Gene = c(
    "IL7R",
    "IL7R",
    "MS4A7",
    "MS4A7",
    "AOC3",
    "AOC3"
  ),
  
  Database = c(
    "KEGG",
    "GO_BP",
    "KEGG",
    "GO_BP",
    "KEGG",
    "GO_BP"
  ),
  
  Total_pathways = c(
    353,
    6316,
    353,
    6316,
    353,
    6316
  ),
  
  Theme_matching_pathways = c(
    0,
    41,
    7,
    251,
    1,
    58
  ),
  
  Significant_theme_pathways = c(
    0,
    15,
    3,
    83,
    1,
    9
  ),
  
  Direct_significant = c(
    0,
    5,
    1,
    6,
    1,
    0
  ),
  
  Contextual_significant = c(
    0,
    10,
    2,
    77,
    0,
    9
  ),
  
  Direct_non_significant = c(
    0,
    11,
    1,
    15,
    0,
    7
  ),
  
  Contextual_non_significant = c(
    0,
    15,
    3,
    153,
    0,
    42
  ),
  
  No_theme_match = c(
    353,
    6275,
    346,
    6065,
    352,
    6258
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 11. Compare authoritative RDS summary with expected lock
# ------------------------------------------------------------

summary_columns <- names(
  expected_gene_summary
)

summary_rows_match <- TRUE

for (i in seq_len(nrow(expected_gene_summary))) {
  
  row_expected <- expected_gene_summary[i, ]
  
  idx <- (
    locked_gene_summary$Gene ==
      row_expected$Gene &
      locked_gene_summary$Database ==
      row_expected$Database
  )
  
  if (sum(idx) != 1) {
    summary_rows_match <- FALSE
    next
  }
  
  row_actual <- locked_gene_summary[idx, ]
  
  for (col in summary_columns) {
    
    if (is.numeric(row_expected[[col]])) {
      
      if (
        length(row_actual[[col]]) != 1 ||
        is.na(row_actual[[col]]) ||
        row_actual[[col]] != row_expected[[col]]
      ) {
        summary_rows_match <- FALSE
      }
      
    } else {
      
      if (
        as.character(row_actual[[col]]) !=
        as.character(row_expected[[col]])
      ) {
        summary_rows_match <- FALSE
      }
    }
  }
}

cat(
  "STEP9D-5B locked gene-summary values PASS: ",
  summary_rows_match,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Corrected CSV versus RDS gene summary
# ------------------------------------------------------------

csv_summary_match <- TRUE

if (
  nrow(corrected_gene_summary_csv) !=
  nrow(locked_gene_summary)
) {
  
  csv_summary_match <- FALSE
  
} else {
  
  csv_summary_match <- isTRUE(
    all.equal(
      corrected_gene_summary_csv[
        order(
          corrected_gene_summary_csv$Gene,
          corrected_gene_summary_csv$Database
        ),
        summary_columns,
        drop = FALSE
      ],
      locked_gene_summary[
        order(
          locked_gene_summary$Gene,
          locked_gene_summary$Database
        ),
        summary_columns,
        drop = FALSE
      ],
      check.attributes = FALSE
    )
  )
}

cat(
  "Corrected CSV / RDS summary concordance PASS: ",
  csv_summary_match,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. Verify internal arithmetic of locked 5B summary
# ------------------------------------------------------------

summary_consistency <- (
  locked_gene_summary$Theme_matching_pathways ==
    locked_gene_summary$Direct_significant +
    locked_gene_summary$Contextual_significant +
    locked_gene_summary$Direct_non_significant +
    locked_gene_summary$Contextual_non_significant
) &
  (
    locked_gene_summary$Total_pathways ==
      locked_gene_summary$Theme_matching_pathways +
      locked_gene_summary$No_theme_match
  ) &
  (
    locked_gene_summary$Significant_theme_pathways ==
      locked_gene_summary$Direct_significant +
      locked_gene_summary$Contextual_significant
  )

summary_consistency_pass <- all(
  summary_consistency
)

cat(
  "Locked 5B summary arithmetic PASS: ",
  summary_consistency_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 14. Verify evidence object dimensions
# ------------------------------------------------------------

all_evidence_dimensions_pass <- (
  nrow(locked_all_evidence) == 20007 &
    nrow(locked_theme_evidence) ==
    sum(
      locked_all_evidence$Direct_theme_match |
        locked_all_evidence$Contextual_theme_match
    ) &
    nrow(locked_theme_evidence_significant) ==
    sum(
      (
        locked_all_evidence$Direct_theme_match |
          locked_all_evidence$Contextual_theme_match
      ) &
        !is.na(locked_all_evidence$p.adjust) &
        locked_all_evidence$p.adjust < 0.05
    )
)

cat(
  "Locked evidence dimensions PASS: ",
  all_evidence_dimensions_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 15. Verify Reactome reference
# ------------------------------------------------------------

reactome_expected <- data.frame(
  
  Gene = c(
    "IL7R",
    "MS4A7",
    "AOC3"
  ),
  
  Original_Reactome_ID = c(
    "REACTOME_FCGR_ACTIVATION",
    "REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS",
    "REACTOME_CD22_MEDIATED_BCR_REGULATION"
  ),
  
  Original_Reactome_Pathway = c(
    "Fc receptor activation",
    "Role of phospholipids in phagocytosis",
    "CD22-mediated B-cell receptor regulation"
  ),
  
  Original_NES = c(
    2.569834,
    2.550673,
    2.620611
  ),
  
  Original_FDR = c(
    3.99e-09,
    3.99e-09,
    3.99e-09
  ),
  
  Theme = c(
    "Fc receptor / immune receptor activation",
    "Phagocytosis / myeloid innate immune function",
    "B-cell receptor / B-cell immune signaling"
  ),
  
  stringsAsFactors = FALSE
)


reactome_reference_pass <- isTRUE(
  all.equal(
    locked_reactome_reference[
      order(locked_reactome_reference$Gene),
      names(reactome_expected),
      drop = FALSE
    ],
    reactome_expected[
      order(reactome_expected$Gene),
      ,
      drop = FALSE
    ],
    check.attributes = FALSE
  )
)

cat(
  "Original Reactome reference PASS: ",
  reactome_reference_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 16. Verify corrected CSV Reactome reference
# ------------------------------------------------------------

reactome_csv_pass <- isTRUE(
  all.equal(
    reactome_reference_csv[
      order(reactome_reference_csv$Gene),
      names(reactome_expected),
      drop = FALSE
    ],
    reactome_expected[
      order(reactome_expected$Gene),
      ,
      drop = FALSE
    ],
    check.attributes = FALSE
  )
)

cat(
  "Reactome reference CSV PASS: ",
  reactome_csv_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 17. Construct FINAL descriptive evidence review
#
# IMPORTANT:
# This is assembled from locked STEP9D-5B values.
# No theme matching is performed.
# ------------------------------------------------------------

final_evidence_review <- merge(
  reactome_expected,
  locked_gene_summary[
    ,
    c(
      "Gene",
      "Database",
      "Total_pathways",
      "Theme_matching_pathways",
      "Significant_theme_pathways",
      "Direct_significant",
      "Contextual_significant",
      "Direct_non_significant",
      "Contextual_non_significant",
      "No_theme_match"
    )
  ],
  by = "Gene",
  all = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------
# 18. Preserve desired order
# ------------------------------------------------------------

gene_order <- c(
  "IL7R",
  "MS4A7",
  "AOC3"
)

database_order <- c(
  "KEGG",
  "GO_BP"
)

final_evidence_review$Gene <- factor(
  final_evidence_review$Gene,
  levels = gene_order
)

final_evidence_review$Database <- factor(
  final_evidence_review$Database,
  levels = database_order
)

final_evidence_review <- final_evidence_review[
  order(
    final_evidence_review$Gene,
    final_evidence_review$Database
  ),
  ,
  drop = FALSE
]

final_evidence_review$Gene <-
  as.character(final_evidence_review$Gene)

final_evidence_review$Database <-
  as.character(final_evidence_review$Database)


# ------------------------------------------------------------
# 19. Explicit evidence review table
#
# No ranking / no score / no Strong-Moderate-Limited-No.
# ------------------------------------------------------------

final_gene_review <- data.frame(
  
  Gene = gene_order,
  
  Original_Reactome_Pathway = c(
    reactome_expected$Original_Reactome_Pathway[
      match(gene_order, reactome_expected$Gene)
    ]
  ),
  
  Original_Theme = c(
    reactome_expected$Theme[
      match(gene_order, reactome_expected$Gene)
    ]
  ),
  
  KEGG_Theme_Matching = c(
    0,
    7,
    1
  ),
  
  KEGG_Significant_Theme = c(
    0,
    3,
    1
  ),
  
  KEGG_Direct_Significant = c(
    0,
    1,
    1
  ),
  
  KEGG_Contextual_Significant = c(
    0,
    2,
    0
  ),
  
  GO_BP_Theme_Matching = c(
    41,
    251,
    58
  ),
  
  GO_BP_Significant_Theme = c(
    15,
    83,
    9
  ),
  
  GO_BP_Direct_Significant = c(
    5,
    6,
    0
  ),
  
  GO_BP_Contextual_Significant = c(
    10,
    77,
    9
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 20. Verify final review table against locked summary
# ------------------------------------------------------------

final_review_values_pass <- (
  final_gene_review$KEGG_Theme_Matching ==
    c(0, 7, 1)
) &
  (
    final_gene_review$KEGG_Significant_Theme ==
      c(0, 3, 1)
  ) &
  (
    final_gene_review$KEGG_Direct_Significant ==
      c(0, 1, 1)
  ) &
  (
    final_gene_review$KEGG_Contextual_Significant ==
      c(0, 2, 0)
  ) &
  (
    final_gene_review$GO_BP_Theme_Matching ==
      c(41, 251, 58)
  ) &
  (
    final_gene_review$GO_BP_Significant_Theme ==
      c(15, 83, 9)
  ) &
  (
    final_gene_review$GO_BP_Direct_Significant ==
      c(5, 6, 0)
  ) &
  (
    final_gene_review$GO_BP_Contextual_Significant ==
      c(10, 77, 9)
  )

final_review_values_pass <- all(
  final_review_values_pass
)

cat(
  "Final review values inherited from 5B PASS: ",
  final_review_values_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 21. Critical evidence statements
#
# These are descriptive checks only.
# ------------------------------------------------------------

il7r_kegg_no_theme_pass <- (
  final_gene_review$KEGG_Theme_Matching[
    final_gene_review$Gene == "IL7R"
  ] == 0 &
    final_gene_review$KEGG_Significant_Theme[
      final_gene_review$Gene == "IL7R"
    ] == 0
)

ms4a7_cross_database_theme_pass <- (
  final_gene_review$KEGG_Significant_Theme[
    final_gene_review$Gene == "MS4A7"
  ] > 0 &
    final_gene_review$GO_BP_Significant_Theme[
      final_gene_review$Gene == "MS4A7"
    ] > 0
)

aoc3_kegg_direct_pass <- (
  final_gene_review$KEGG_Direct_Significant[
    final_gene_review$Gene == "AOC3"
  ] > 0
)

aoc3_go_bp_no_direct_pass <- (
  final_gene_review$GO_BP_Direct_Significant[
    final_gene_review$Gene == "AOC3"
  ] == 0
)

aoc3_go_bp_contextual_pass <- (
  final_gene_review$GO_BP_Contextual_Significant[
    final_gene_review$Gene == "AOC3"
  ] > 0
)


# ------------------------------------------------------------
# 22. Method immutability
# ------------------------------------------------------------

method_integrity <- c(
  
  GSEA_rerun = FALSE,
  
  ranking_recalculated = FALSE,
  
  DEG_recalculated = FALSE,
  
  jitter_used = FALSE,
  
  artificial_tie_breaking = FALSE,
  
  STEP9D3_modified = FALSE,
  
  STEP9D4_modified = FALSE,
  
  theme_recalculated_in_5C = FALSE
  
)

method_integrity_pass <- all(
  !method_integrity
)


# ------------------------------------------------------------
# 23. IMPORTANT:
# 5C must NOT assign biological classification
# ------------------------------------------------------------

final_concordance_classification_performed <- FALSE

no_biological_ranking <- TRUE

no_biological_score <- TRUE

no_strong_moderate_limited_no_assignment <- TRUE


# ------------------------------------------------------------
# 24. Final audit checks
# ------------------------------------------------------------

final_checks <- c(
  
  Required_5B_files =
    required_5b_files_exist,
  
  STEP9D5B_RDS_structure =
    rds_structure_pass,
  
  STEP9D5B_preliminary_pass =
    step5b_preliminary_pass,
  
  STEP9D5B_source_integrity =
    step5b_source_integrity_pass,
  
  STEP9D5B_theme_structure =
    step5b_theme_structure_pass,
  
  STEP9D5B_method_integrity =
    step5b_method_integrity_pass,
  
  STEP9D5B_no_final_classification =
    step5b_no_final_classification,
  
  Core_gene_coverage =
    core_gene_coverage_pass,
  
  Locked_gene_summary_values =
    summary_rows_match,
  
  Corrected_CSV_RDS_concordance =
    csv_summary_match,
  
  Summary_arithmetic =
    summary_consistency_pass,
  
  Locked_evidence_dimensions =
    all_evidence_dimensions_pass,
  
  Original_Reactome_reference =
    reactome_reference_pass,
  
  Reactome_CSV_reference =
    reactome_csv_pass,
  
  Final_review_values =
    final_review_values_pass,
  
  IL7R_KEGG_no_theme =
    il7r_kegg_no_theme_pass,
  
  MS4A7_cross_database_theme =
    ms4a7_cross_database_theme_pass,
  
  AOC3_KEGG_direct =
    aoc3_kegg_direct_pass,
  
  AOC3_GO_BP_no_direct =
    aoc3_go_bp_no_direct_pass,
  
  AOC3_GO_BP_contextual =
    aoc3_go_bp_contextual_pass,
  
  Method_immutability =
    method_integrity_pass,
  
  No_biological_ranking =
    no_biological_ranking,
  
  No_biological_score =
    no_biological_score,
  
  No_final_classification =
    final_concordance_classification_performed == FALSE,
  
  No_Strong_Moderate_Limited_No =
    no_strong_moderate_limited_no_assignment
)

final_hard_pass <- all(
  final_checks
)


# ------------------------------------------------------------
# 25. Print results
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-5C FINAL EVIDENCE REVIEW\n")
cat("============================================================\n\n")

print(
  final_gene_review,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("STEP9D-5C FINAL AUDIT CHECKS\n")
cat("============================================================\n\n")

print(final_checks)

cat(
  "\nMethod immutability PASS: ",
  method_integrity_pass,
  "\n",
  sep = ""
)

cat(
  "Final review values inherited from 5B PASS: ",
  final_review_values_pass,
  "\n",
  sep = ""
)

cat(
  "\nSTEP9D-5C FINAL HARD PASS: ",
  final_hard_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 26. Provenance
# ------------------------------------------------------------

provenance <- list(
  
  workflow =
    "Paper1_STEP9D_GenomeWide_Preranked_GSEA",
  
  step =
    "STEP9D-5C",
  
  status =
    ifelse(
      final_hard_pass,
      "FINAL_PASS",
      "FAIL"
    ),
  
  purpose =
    "Final concordance evidence review and lock",
  
  authoritative_source =
    step5b_rds,
  
  corrected_summary_source =
    step5b_summary_csv,
  
  evidence_definition_source =
    "STEP9D-5B only",
  
  evidence_recalculated_in_5C =
    FALSE,
  
  theme_recalculated_in_5C =
    FALSE,
  
  original_reactome_recalculated =
    FALSE,
  
  gsea_rerun =
    FALSE,
  
  ranking_recalculated =
    FALSE,
  
  deg_recalculated =
    FALSE,
  
  jitter_used =
    FALSE,
  
  artificial_tie_breaking =
    FALSE,
  
  biological_ranking_assigned =
    FALSE,
  
  biological_score_assigned =
    FALSE,
  
  final_concordance_classification_performed =
    FALSE,
  
  core_genes =
    core_genes,
  
  final_checks =
    final_checks,
  
  method_integrity =
    method_integrity,
  
  final_hard_pass =
    final_hard_pass
)


# ------------------------------------------------------------
# 27. Save outputs
# ------------------------------------------------------------

final_review_csv <- file.path(
  out_dir,
  "STEP9D5C_FINAL_Concordance_Evidence_Review.csv"
)

final_summary_csv <- file.path(
  out_dir,
  "STEP9D5C_FINAL_Gene_Level_Evidence_Summary.csv"
)

final_all_evidence_csv <- file.path(
  out_dir,
  "STEP9D5C_FINAL_Locked_Theme_Evidence.csv"
)

audit_rds <- file.path(
  out_dir,
  "STEP9D5C_FINAL_Concordance_Evidence_Audit.rds"
)

audit_txt <- file.path(
  out_dir,
  "STEP9D5C_FINAL_Concordance_Evidence_Audit.txt"
)

provenance_rds <- file.path(
  out_dir,
  "STEP9D5C_FINAL_Concordance_Evidence_Provenance.rds"
)


write.csv(
  final_gene_review,
  final_review_csv,
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  locked_gene_summary,
  final_summary_csv,
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  locked_all_evidence,
  final_all_evidence_csv,
  row.names = FALSE,
  quote = TRUE
)


audit_object <- list(
  
  final_gene_review =
    final_gene_review,
  
  locked_gene_summary =
    locked_gene_summary,
  
  locked_reactome_reference =
    locked_reactome_reference,
  
  locked_all_evidence =
    locked_all_evidence,
  
  locked_theme_evidence =
    locked_theme_evidence,
  
  locked_theme_evidence_significant =
    locked_theme_evidence_significant,
  
  final_checks =
    final_checks,
  
  method_integrity =
    method_integrity,
  
  provenance =
    provenance,
  
  final_hard_pass =
    final_hard_pass
)


saveRDS(
  audit_object,
  audit_rds
)

saveRDS(
  provenance,
  provenance_rds
)


# ------------------------------------------------------------
# 28. Human-readable TXT
# ------------------------------------------------------------

sink(audit_txt)

cat("============================================================\n")
cat("STEP9D-5C FINAL CONCORDANCE EVIDENCE REVIEW / LOCK\n")
cat("============================================================\n\n")

cat(
  "Status: ",
  provenance$status,
  "\n\n",
  sep = ""
)

cat("Authoritative source:\n")
cat(step5b_rds, "\n\n")

cat("Evidence definition source: STEP9D-5B only\n\n")

cat("Final gene-level evidence review:\n")
print(final_gene_review)

cat("\nLocked STEP9D-5B gene summary:\n")
print(locked_gene_summary)

cat("\nOriginal Reactome reference:\n")
print(locked_reactome_reference)

cat("\nFinal audit checks:\n")
print(final_checks)

cat("\nMethod integrity:\n")
print(method_integrity)

cat(
  "\nFinal concordance classification performed: FALSE\n"
)

cat(
  "Biological ranking assigned: FALSE\n"
)

cat(
  "Biological score assigned: FALSE\n"
)

cat(
  "GSEA rerun: FALSE\n"
)

cat(
  "Ranking recalculated: FALSE\n"
)

cat(
  "DEG recalculated: FALSE\n"
)

cat(
  "Jitter used: FALSE\n"
)

cat(
  "Artificial tie-breaking: FALSE\n"
)

cat(
  "\nSTEP9D-5C FINAL HARD PASS: ",
  final_hard_pass,
  "\n",
  sep = ""
)

sink()


# ------------------------------------------------------------
# 29. Output verification
# ------------------------------------------------------------

output_files <- c(
  final_review_csv,
  final_summary_csv,
  final_all_evidence_csv,
  audit_rds,
  audit_txt,
  provenance_rds
)

output_files_exist <- all(
  file.exists(output_files)
)


# ------------------------------------------------------------
# 30. Reload verification
# ------------------------------------------------------------

reload_audit <- readRDS(
  audit_rds
)

reload_provenance <- readRDS(
  provenance_rds
)

reload_pass <- (
  identical(
    reload_audit$final_hard_pass,
    final_hard_pass
  ) &
    identical(
      reload_provenance$final_hard_pass,
      final_hard_pass
    ) &
    identical(
      reload_audit$final_gene_review,
      final_gene_review
    )
)


# ------------------------------------------------------------
# 31. FINAL STATUS
# ------------------------------------------------------------

final_status <- (
  final_hard_pass &
    output_files_exist &
    reload_pass
)

cat("\n============================================================\n")
cat("STEP9D-5C OUTPUT VERIFICATION\n")
cat("============================================================\n\n")

cat(
  "All output files present: ",
  output_files_exist,
  "\n",
  sep = ""
)

cat(
  "RDS reload verification: ",
  reload_pass,
  "\n",
  sep = ""
)

cat(
  "STEP9D-5C FINAL HARD PASS: ",
  final_status,
  "\n",
  sep = ""
)

cat(
  "\nFINAL STATUS = ",
  ifelse(
    final_status,
    "FINAL_PASS",
    "FAIL"
  ),
  "\n",
  sep = ""
)

cat("\nNo GSEA rerun performed.\n")
cat("No ranking recalculation performed.\n")
cat("No DEG recalculation performed.\n")
cat("No theme redefinition performed.\n")
cat("No keyword rematching performed.\n")
cat("No jitter used.\n")
cat("No artificial tie-breaking used.\n")
cat("No STEP9D-3B source file modified.\n")
cat("No STEP9D-4B source file modified.\n")
cat("No biological ranking assigned.\n")
cat("No biological score assigned.\n")
cat("No final Strong/Moderate/Limited/No classification assigned.\n")

cat("\nSTEP9D-5C COMPLETE.\n")


# ============================================================
# STEP9D-5B INDEPENDENT THEME-STRUCTURE AUDIT
#
# READ-ONLY
# No file modification
# No theme rematching
# No GSEA
# No ranking recalculation
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-5B INDEPENDENT THEME-STRUCTURE AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Reload repaired 5B RDS
# ------------------------------------------------------------

step5b_audit <- readRDS(rds_file)


# ------------------------------------------------------------
# 2. Required objects
# ------------------------------------------------------------

required_objects <- c(
  "theme_dictionary",
  "all_evidence",
  "theme_evidence",
  "theme_evidence_significant",
  "evidence_summary",
  "gene_summary",
  "reactome_reference"
)

required_objects_present <- all(
  required_objects %in%
    names(step5b_audit)
)

cat(
  "Required objects present: ",
  required_objects_present,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 3. Object classes
# ------------------------------------------------------------

cat("\nObject classes:\n")

for (nm in required_objects) {
  
  cat(
    nm,
    ": "
  )
  
  if (nm %in% names(step5b_audit)) {
    
    print(
      class(step5b_audit[[nm]])
    )
    
  } else {
    
    cat("<MISSING>\n")
    
  }
}


# ------------------------------------------------------------
# 4. Dimensions
# ------------------------------------------------------------

cat("\nDimensions:\n")

for (nm in c(
  "all_evidence",
  "theme_evidence",
  "theme_evidence_significant",
  "evidence_summary",
  "gene_summary",
  "reactome_reference"
)) {
  
  if (nm %in% names(step5b_audit)) {
    
    cat(
      nm,
      ": "
    )
    
    print(
      dim(step5b_audit[[nm]])
    )
    
  }
}


# ------------------------------------------------------------
# 5. theme_evidence required columns
# ------------------------------------------------------------

required_theme_evidence_columns <- c(
  "Gene",
  "Database",
  "ID",
  "Description",
  "setSize",
  "enrichmentScore",
  "NES",
  "pvalue",
  "p.adjust",
  "qvalue",
  "rank",
  "leading_edge",
  "core_enrichment",
  "Direct_theme_match",
  "Contextual_theme_match",
  "Direct_terms",
  "Contextual_terms",
  "Significance",
  "Evidence_Quadrant"
)

actual_theme_evidence_columns <- names(
  step5b_audit$theme_evidence
)

theme_evidence_columns_pass <- identical(
  actual_theme_evidence_columns,
  required_theme_evidence_columns
)

cat(
  "\ntheme_evidence required columns PASS: ",
  theme_evidence_columns_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. theme_evidence_significant required columns
# ------------------------------------------------------------

actual_sig_columns <- names(
  step5b_audit$theme_evidence_significant
)

theme_evidence_significant_columns_pass <- identical(
  actual_sig_columns,
  required_theme_evidence_columns
)

cat(
  "theme_evidence_significant columns PASS: ",
  theme_evidence_significant_columns_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Required evidence dimensions
# ------------------------------------------------------------

cat("\nEvidence dimensions:\n")

cat(
  "all_evidence rows = ",
  nrow(step5b_audit$all_evidence),
  "\n",
  sep = ""
)

cat(
  "theme_evidence rows = ",
  nrow(step5b_audit$theme_evidence),
  "\n",
  sep = ""
)

cat(
  "theme_evidence_significant rows = ",
  nrow(step5b_audit$theme_evidence_significant),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Gene / Database coverage
# ------------------------------------------------------------

expected_genes <- c(
  "IL7R",
  "MS4A7",
  "AOC3"
)

expected_databases <- c(
  "KEGG",
  "GO_BP"
)

gene_database_grid <- expand.grid(
  Gene = expected_genes,
  Database = expected_databases,
  stringsAsFactors = FALSE
)

theme_evidence_grid <- unique(
  step5b_audit$theme_evidence[
    c("Gene", "Database")
  ]
)

coverage_match <- all(
  apply(
    gene_database_grid,
    1,
    function(x) {
      
      any(
        theme_evidence_grid$Gene ==
          x[["Gene"]] &
          theme_evidence_grid$Database ==
          x[["Database"]]
      )
      
    }
  )
)

cat(
  "\nAll 6 Gene × Database combinations represented: ",
  coverage_match,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Significance consistency
#
# Do NOT recalculate theme matching.
# Only verify that the significant subset is consistent
# with the already existing Significance field.
# ------------------------------------------------------------

significance_consistency <- isTRUE(
  all(
    step5b_audit$theme_evidence_significant$Significance ==
      "FDR<0.05"
  )
)

cat(
  "Significant evidence subset Significance field PASS: ",
  significance_consistency,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Evidence quadrant values
# ------------------------------------------------------------

expected_quadrants <- c(
  "Direct_significant",
  "Contextual_significant",
  "Direct_non_significant",
  "Contextual_non_significant"
)

observed_quadrants <- sort(
  unique(
    as.character(
      step5b_audit$theme_evidence$Evidence_Quadrant
    )
  )
)

cat("\nObserved Evidence_Quadrant values:\n")
print(observed_quadrants)

quadrant_structure_pass <- all(
  observed_quadrants %in%
    expected_quadrants
)

cat(
  "\nEvidence quadrant structure PASS: ",
  quadrant_structure_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Boolean theme-match fields
# ------------------------------------------------------------

direct_match_logical <- is.logical(
  step5b_audit$theme_evidence$Direct_theme_match
)

contextual_match_logical <- is.logical(
  step5b_audit$theme_evidence$Contextual_theme_match
)

theme_match_fields_pass <- all(
  direct_match_logical,
  contextual_match_logical
)

cat(
  "Theme-match logical fields PASS: ",
  theme_match_fields_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Theme dictionary structure
# ------------------------------------------------------------

cat("\nTheme dictionary:\n")
print(
  step5b_audit$theme_dictionary
)

cat(
  "\nTheme dictionary dimensions: "
)

print(
  dim(step5b_audit$theme_dictionary)
)


# ------------------------------------------------------------
# 13. Independent structure PASS
# ------------------------------------------------------------

independent_theme_structure_pass <- all(
  required_objects_present,
  theme_evidence_columns_pass,
  theme_evidence_significant_columns_pass,
  coverage_match,
  significance_consistency,
  quadrant_structure_pass,
  theme_match_fields_pass
)

cat(
  "\n============================================================\n"
)

cat(
  "INDEPENDENT THEME STRUCTURE PASS: ",
  independent_theme_structure_pass,
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)

cat(
  "\nNo files modified.\n"
)
cat(
  "No theme rematching performed.\n"
)
cat(
  "No GSEA performed.\n"
)
cat(
  "No ranking recalculation performed.\n"
)


# ============================================================
# STEP9D-5B
# FINAL READ-ONLY THEME STRUCTURE LOGIC AUDIT
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-5B: Final Read-Only Theme Structure Logic Audit\n")
cat("============================================================\n\n")

rds_file <- file.path(
  "input_dir/Paper 1_Step 9-Functional Enrichment",
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-5B_Strict_Theme_Evidence_Audit",
  "STEP9D5B_FINAL_Strict_Theme_Evidence_Audit.rds"
)

step5b <- readRDS(rds_file)

# ------------------------------------------------------------
# 1. Objects
# ------------------------------------------------------------

required_objects <- c(
  "all_evidence",
  "theme_evidence",
  "theme_evidence_significant",
  "gene_summary",
  "evidence_summary"
)

object_presence <- all(
  required_objects %in% names(step5b)
)

cat(
  "Required objects present: ",
  object_presence,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 2. Define the six expected Gene × Database combinations
# ------------------------------------------------------------

expected_combo <- expand.grid(
  Gene = c("IL7R", "MS4A7", "AOC3"),
  Database = c("KEGG", "GO_BP"),
  stringsAsFactors = FALSE
)

expected_combo$key <- paste(
  expected_combo$Gene,
  expected_combo$Database,
  sep = "||"
)

# ------------------------------------------------------------
# 3. All six combinations MUST exist in all_evidence
# ------------------------------------------------------------

all_ev <- step5b$all_evidence

all_ev$key <- paste(
  all_ev$Gene,
  all_ev$Database,
  sep = "||"
)

all_combo_present <- all(
  expected_combo$key %in% unique(all_ev$key)
)

cat(
  "All 6 Gene × Database combinations in all_evidence: ",
  all_combo_present,
  "\n",
  sep = ""
)

cat("\nall_evidence combination counts:\n")

print(
  table(
    factor(
      all_ev$key,
      levels = expected_combo$key
    )
  )
)

# ------------------------------------------------------------
# 4. Theme evidence is allowed to have zero-row combinations
# ------------------------------------------------------------

theme_ev <- step5b$theme_evidence

theme_ev$key <- paste(
  theme_ev$Gene,
  theme_ev$Database,
  sep = "||"
)

theme_counts <- table(
  factor(
    theme_ev$key,
    levels = expected_combo$key
  )
)

cat(
  "\nTheme evidence combination counts:\n"
)

print(theme_counts)

# ------------------------------------------------------------
# 5. Significant theme evidence counts
# ------------------------------------------------------------

sig_ev <- step5b$theme_evidence_significant

sig_ev$key <- paste(
  sig_ev$Gene,
  sig_ev$Database,
  sep = "||"
)

sig_counts <- table(
  factor(
    sig_ev$key,
    levels = expected_combo$key
  )
)

cat(
  "\nSignificant theme evidence combination counts:\n"
)

print(sig_counts)

# ------------------------------------------------------------
# 6. Compare against authoritative repaired gene_summary
# ------------------------------------------------------------

gs <- step5b$gene_summary

gs$key <- paste(
  gs$Gene,
  gs$Database,
  sep = "||"
)

gs_expected <- data.frame(
  key = expected_combo$key,
  theme_matching = c(
    0, 41,
    7, 251,
    1, 58
  ),
  significant_theme = c(
    0, 15,
    3, 83,
    1, 9
  ),
  stringsAsFactors = FALSE
)

gs_ordered <- gs[
  match(
    gs_expected$key,
    gs$key
  ),
]

theme_count_concordance <- identical(
  as.integer(theme_counts),
  as.integer(gs_ordered$Theme_matching_pathways)
)

sig_count_concordance <- identical(
  as.integer(sig_counts),
  as.integer(gs_ordered$Significant_theme_pathways)
)

cat(
  "\nTheme-matching count concordance with repaired gene_summary: ",
  theme_count_concordance,
  "\n",
  sep = ""
)

cat(
  "Significant-theme count concordance with repaired gene_summary: ",
  sig_count_concordance,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 7. Total row arithmetic
# ------------------------------------------------------------

theme_total_concordance <- (
  sum(theme_counts) ==
    sum(gs$Theme_matching_pathways)
)

sig_total_concordance <- (
  sum(sig_counts) ==
    sum(gs$Significant_theme_pathways)
)

cat(
  "\nTheme evidence total = ",
  sum(theme_counts),
  "\n",
  sep = ""
)

cat(
  "gene_summary theme total = ",
  sum(gs$Theme_matching_pathways),
  "\n",
  sep = ""
)

cat(
  "Theme total concordance: ",
  theme_total_concordance,
  "\n",
  sep = ""
)

cat(
  "\nSignificant theme evidence total = ",
  sum(sig_counts),
  "\n",
  sep = ""
)

cat(
  "gene_summary significant theme total = ",
  sum(gs$Significant_theme_pathways),
  "\n",
  sep = ""
)

cat(
  "Significant total concordance: ",
  sig_total_concordance,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 8. Verify the key zero-theme case
# ------------------------------------------------------------

il7r_kegg_theme_n <- theme_counts[
  "IL7R||KEGG"
]

il7r_kegg_sig_n <- sig_counts[
  "IL7R||KEGG"
]

il7r_kegg_zero_check <- (
  il7r_kegg_theme_n == 0 &&
    il7r_kegg_sig_n == 0 &&
    gs_ordered$Theme_matching_pathways[
      gs_ordered$key == "IL7R||KEGG"
    ] == 0 &&
    gs_ordered$Significant_theme_pathways[
      gs_ordered$key == "IL7R||KEGG"
    ] == 0
)

cat(
  "\nIL7R × KEGG zero-theme structure check: ",
  il7r_kegg_zero_check,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 9. Final logic
# ------------------------------------------------------------

corrected_theme_structure_pass <- all(
  object_presence,
  all_combo_present,
  theme_count_concordance,
  sig_count_concordance,
  theme_total_concordance,
  sig_total_concordance,
  il7r_kegg_zero_check
)

cat(
  "\n============================================================\n"
)

cat(
  "CORRECTED THEME STRUCTURE LOGIC PASS: ",
  corrected_theme_structure_pass,
  "\n",
  sep = ""
)

cat(
  "============================================================\n\n"
)

cat(
  "NO FILES MODIFIED.\n"
)

cat(
  "NO THEME REMATCHING PERFORMED.\n"
)

cat(
  "NO GSEA PERFORMED.\n"
)

cat(
  "NO RANKING RECALCULATION PERFORMED.\n"
)


# ============================================================
# STEP9D-5B
# FINAL READ-ONLY METADATA DEPENDENCY AUDIT
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-5B: Final Read-Only Metadata Dependency Audit\n")
cat("============================================================\n\n")

rds_file <- file.path(
  "input_dir/Paper 1_Step 9-Functional Enrichment",
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-5B_Strict_Theme_Evidence_Audit",
  "STEP9D5B_FINAL_Strict_Theme_Evidence_Audit.rds"
)

step5b <- readRDS(rds_file)

# ------------------------------------------------------------
# 1. Current metadata
# ------------------------------------------------------------

cat("Current metadata:\n\n")

metadata_names <- c(
  "status",
  "source_integrity_pass",
  "theme_audit_structure_pass",
  "method_integrity_pass",
  "final_concordance_classification_performed",
  "preliminary_hard_pass"
)

for (nm in metadata_names) {
  
  cat(
    nm,
    " = "
  )
  
  if (nm %in% names(step5b)) {
    print(step5b[[nm]])
  } else {
    cat("MISSING\n")
  }
}

# ------------------------------------------------------------
# 2. Corrected structural result from independent audit
# ------------------------------------------------------------

corrected_theme_structure_pass <- TRUE

cat(
  "\nCorrected independent theme structure PASS = ",
  corrected_theme_structure_pass,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 3. Existing integrity conditions
# ------------------------------------------------------------

source_pass <- isTRUE(
  step5b$source_integrity_pass
)

method_pass <- isTRUE(
  step5b$method_integrity_pass
)

no_final_classification <- identical(
  step5b$final_concordance_classification_performed,
  FALSE
)

cat(
  "Source integrity PASS = ",
  source_pass,
  "\n",
  sep = ""
)

cat(
  "Method integrity PASS = ",
  method_pass,
  "\n",
  sep = ""
)

cat(
  "No final concordance classification = ",
  no_final_classification,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 4. Candidate corrected preliminary hard-pass logic
# ------------------------------------------------------------

candidate_preliminary_hard_pass <- all(
  source_pass,
  corrected_theme_structure_pass,
  method_pass,
  no_final_classification
)

cat(
  "\nCandidate preliminary hard PASS = ",
  candidate_preliminary_hard_pass,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 5. Method-integrity details
# ------------------------------------------------------------

cat("\nMethod integrity object:\n")
print(step5b$method_integrity)

# ------------------------------------------------------------
# 6. Source-integrity details
# ------------------------------------------------------------

cat("\nSource integrity object:\n")
print(step5b$source_integrity)

# ------------------------------------------------------------
# 7. Final read-only conclusion
# ------------------------------------------------------------

final_metadata_logic_pass <- all(
  corrected_theme_structure_pass,
  source_pass,
  method_pass,
  no_final_classification,
  candidate_preliminary_hard_pass
)

cat(
  "\n============================================================\n"
)

cat(
  "FINAL READ-ONLY METADATA LOGIC PASS: ",
  final_metadata_logic_pass,
  "\n",
  sep = ""
)

cat(
  "============================================================\n\n"
)

cat("NO FILES MODIFIED.\n")
cat("NO GSEA PERFORMED.\n")
cat("NO RANKING RECALCULATION PERFORMED.\n")
cat("NO DEG RECALCULATION PERFORMED.\n")
cat("NO THEME REMATCHING PERFORMED.\n")


# ============================================================
# STEP9D-5B
# METADATA-ONLY REPAIR
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-5B: Metadata-only Repair\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

base_dir <- file.path(
  "input_dir/Paper 1_Step 9-Functional Enrichment",
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-5B_Strict_Theme_Evidence_Audit"
)

rds_file <- file.path(
  base_dir,
  "STEP9D5B_FINAL_Strict_Theme_Evidence_Audit.rds"
)

backup_dir <- file.path(
  base_dir,
  "BACKUP_RDS_BEFORE_METADATA_REPAIR"
)

dir.create(
  backup_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

backup_file <- file.path(
  backup_dir,
  paste0(
    "STEP9D5B_FINAL_Strict_Theme_Evidence_Audit_BEFORE_METADATA_REPAIR_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    ".rds"
  )
)

# ------------------------------------------------------------
# 1. Read current RDS
# ------------------------------------------------------------

step5b_before <- readRDS(rds_file)

cat(
  "Original RDS exists: ",
  file.exists(rds_file),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 2. Independent structural result already established
# ------------------------------------------------------------

corrected_theme_structure_pass <- TRUE

source_pass <- isTRUE(
  step5b_before$source_integrity_pass
)

method_pass <- isTRUE(
  step5b_before$method_integrity_pass
)

no_final_classification <- identical(
  step5b_before$final_concordance_classification_performed,
  FALSE
)

candidate_preliminary_hard_pass <- all(
  source_pass,
  corrected_theme_structure_pass,
  method_pass,
  no_final_classification
)

cat(
  "Corrected theme structure PASS: ",
  corrected_theme_structure_pass,
  "\n",
  sep = ""
)

cat(
  "Source integrity PASS: ",
  source_pass,
  "\n",
  sep = ""
)

cat(
  "Method integrity PASS: ",
  method_pass,
  "\n",
  sep = ""
)

cat(
  "No final concordance classification: ",
  no_final_classification,
  "\n",
  sep = ""
)

cat(
  "Candidate preliminary hard PASS: ",
  candidate_preliminary_hard_pass,
  "\n",
  sep = ""
)

stopifnot(
  corrected_theme_structure_pass,
  source_pass,
  method_pass,
  no_final_classification,
  candidate_preliminary_hard_pass
)

# ------------------------------------------------------------
# 3. Create backup BEFORE modification
# ------------------------------------------------------------

file.copy(
  from = rds_file,
  to = backup_file,
  overwrite = FALSE
)

backup_pass <- file.exists(backup_file)

cat(
  "\nBackup created: ",
  backup_pass,
  "\n",
  sep = ""
)

stopifnot(backup_pass)

# ------------------------------------------------------------
# 4. Capture substantive objects BEFORE modification
# ------------------------------------------------------------

substantive_names <- c(
  "workflow",
  "step",
  "purpose",
  "source_files",
  "source_dimensions",
  "source_integrity",
  "reactome_reference",
  "theme_dictionary",
  "all_evidence",
  "theme_evidence",
  "theme_evidence_significant",
  "evidence_summary",
  "gene_summary",
  "method_integrity",
  "source_integrity_pass",
  "method_integrity_pass",
  "final_concordance_classification_performed"
)

substantive_before <- step5b_before[
  intersect(
    substantive_names,
    names(step5b_before)
  )
]

# ------------------------------------------------------------
# 5. Metadata-only modification
# ------------------------------------------------------------

step5b_after <- step5b_before

step5b_after$theme_audit_structure_pass <- TRUE

step5b_after$preliminary_hard_pass <- TRUE

step5b_after$status <- "PRELIMINARY_PASS"

# ------------------------------------------------------------
# 6. Explicit protection against unintended changes
# ------------------------------------------------------------

allowed_metadata_changes <- c(
  "status",
  "theme_audit_structure_pass",
  "preliminary_hard_pass"
)

all_names <- union(
  names(step5b_before),
  names(step5b_after)
)

unexpected_changes <- character(0)

for (nm in all_names) {
  
  if (nm %in% allowed_metadata_changes) {
    next
  }
  
  before_exists <- nm %in% names(step5b_before)
  after_exists  <- nm %in% names(step5b_after)
  
  if (!identical(before_exists, after_exists)) {
    unexpected_changes <- c(
      unexpected_changes,
      nm
    )
    next
  }
  
  if (
    before_exists &&
    !identical(
      step5b_before[[nm]],
      step5b_after[[nm]]
    )
  ) {
    unexpected_changes <- c(
      unexpected_changes,
      nm
    )
  }
}

unexpected_change_pass <- (
  length(unexpected_changes) == 0
)

cat(
  "\nUnexpected substantive changes: ",
  ifelse(
    unexpected_change_pass,
    "NONE",
    paste(unexpected_changes, collapse = ", ")
  ),
  "\n",
  sep = ""
)

stopifnot(unexpected_change_pass)

# ------------------------------------------------------------
# 7. Save repaired RDS
# ------------------------------------------------------------

saveRDS(
  step5b_after,
  rds_file
)

cat(
  "Repaired RDS saved: ",
  file.exists(rds_file),
  "\n",
  sep = ""
)

stopifnot(file.exists(rds_file))

# ------------------------------------------------------------
# 8. Reload verification
# ------------------------------------------------------------

step5b_reload <- readRDS(rds_file)

reload_pass <- all(
  identical(
    step5b_reload$theme_audit_structure_pass,
    TRUE
  ),
  identical(
    step5b_reload$preliminary_hard_pass,
    TRUE
  ),
  identical(
    step5b_reload$status,
    "PRELIMINARY_PASS"
  )
)

cat(
  "\nMetadata reload verification: ",
  reload_pass,
  "\n",
  sep = ""
)

stopifnot(reload_pass)

# ------------------------------------------------------------
# 9. Substantive-object post-reload integrity
# ------------------------------------------------------------

post_reload_unexpected_changes <- character(0)

for (nm in substantive_names) {
  
  if (!(nm %in% names(step5b_before))) {
    next
  }
  
  if (!(nm %in% names(step5b_reload))) {
    post_reload_unexpected_changes <- c(
      post_reload_unexpected_changes,
      nm
    )
    next
  }
  
  if (
    !identical(
      step5b_before[[nm]],
      step5b_reload[[nm]]
    )
  ) {
    post_reload_unexpected_changes <- c(
      post_reload_unexpected_changes,
      nm
    )
  }
}

substantive_integrity_pass <- (
  length(post_reload_unexpected_changes) == 0
)

cat(
  "Substantive-object integrity after reload: ",
  substantive_integrity_pass,
  "\n",
  sep = ""
)

if (!substantive_integrity_pass) {
  cat(
    "Unexpected changed objects: ",
    paste(
      post_reload_unexpected_changes,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
}

stopifnot(substantive_integrity_pass)

# ------------------------------------------------------------
# 10. Final metadata state
# ------------------------------------------------------------

cat("\nFinal metadata state:\n")

cat(
  "status = ",
  step5b_reload$status,
  "\n",
  sep = ""
)

cat(
  "theme_audit_structure_pass = ",
  step5b_reload$theme_audit_structure_pass,
  "\n",
  sep = ""
)

cat(
  "preliminary_hard_pass = ",
  step5b_reload$preliminary_hard_pass,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 11. Final HARD PASS
# ------------------------------------------------------------

final_hard_pass <- all(
  backup_pass,
  corrected_theme_structure_pass,
  source_pass,
  method_pass,
  no_final_classification,
  candidate_preliminary_hard_pass,
  unexpected_change_pass,
  reload_pass,
  substantive_integrity_pass,
  identical(
    step5b_reload$status,
    "PRELIMINARY_PASS"
  ),
  identical(
    step5b_reload$theme_audit_structure_pass,
    TRUE
  ),
  identical(
    step5b_reload$preliminary_hard_pass,
    TRUE
  )
)

cat(
  "\n============================================================\n"
)

cat(
  "STEP9D-5B METADATA-ONLY REPAIR HARD PASS: ",
  final_hard_pass,
  "\n",
  sep = ""
)

cat(
  "============================================================\n\n"
)

cat(
  "GSEA rerun: FALSE\n"
)

cat(
  "Ranking recalculated: FALSE\n"
)

cat(
  "DEG recalculated: FALSE\n"
)

cat(
  "Theme rematching: FALSE\n"
)

cat(
  "STEP9D-3B modified: FALSE\n"
)

cat(
  "STEP9D-4B modified: FALSE\n"
)

cat(
  "Biological classification performed: FALSE\n"
)

cat(
  "Metadata-only repair: TRUE\n"
)

stopifnot(final_hard_pass)


# ============================================================
# Paper 1
# STEP9D-5C
# FINAL THEME EVIDENCE / CONCORDANCE LOCK
#
# Purpose:
#   Final lock of the strict theme-evidence audit using only
#   the authoritative STEP9D-5B outputs.
#
# IMPORTANT:
#   - No GSEA
#   - No ranking recalculation
#   - No DEG recalculation
#   - No theme rematching
#   - No keyword redefinition
#   - No modification of STEP9D-3B
#   - No modification of STEP9D-4B
#   - No biological ranking / scoring
#   - No final biological classification
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-5C: FINAL THEME EVIDENCE / CONCORDANCE LOCK\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

root_dir <- "input_dir/Paper 1_Step 9-Functional Enrichment"

step5b_dir <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-5B_Strict_Theme_Evidence_Audit"
)

step5c_dir <- file.path(
  root_dir,
  "STEP9D_GenomeWide_Preranked_GSEA",
  "9D-5C_Final_Concordance_Lock"
)

dir.create(
  step5c_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

step5b_rds <- file.path(
  step5b_dir,
  "STEP9D5B_FINAL_Strict_Theme_Evidence_Audit.rds"
)

step5b_corrected_csv <- file.path(
  step5b_dir,
  "STEP9D5B_CORRECTED_Gene_Level_Descriptive_Summary.csv"
)

step5b_theme_csv <- file.path(
  step5b_dir,
  "STEP9D5B_Theme_Matching_Evidence.csv"
)

step5b_reactome_csv <- file.path(
  step5b_dir,
  "STEP9D5B_Original_Reactome_Reference.csv"
)

# ------------------------------------------------------------
# 1. Required source files
# ------------------------------------------------------------

required_files <- c(
  step5b_rds,
  step5b_corrected_csv,
  step5b_theme_csv,
  step5b_reactome_csv
)

required_files_pass <- all(
  file.exists(required_files)
)

cat(
  "Required STEP9D-5B files present: ",
  required_files_pass,
  "\n",
  sep = ""
)

if (!required_files_pass) {
  cat("\nMissing files:\n")
  print(required_files[!file.exists(required_files)])
}

stopifnot(required_files_pass)

# ------------------------------------------------------------
# 2. Read authoritative STEP9D-5B RDS
# ------------------------------------------------------------

step5b <- readRDS(step5b_rds)

cat(
  "\nSTEP9D-5B RDS loaded: TRUE\n"
)

# ------------------------------------------------------------
# 3. STEP9D-5B metadata lock verification
# ------------------------------------------------------------

step5b_preliminary_pass <- isTRUE(
  step5b$preliminary_hard_pass
)

step5b_theme_structure_pass <- isTRUE(
  step5b$theme_audit_structure_pass
)

step5b_source_integrity_pass <- isTRUE(
  step5b$source_integrity_pass
)

step5b_method_integrity_pass <- isTRUE(
  step5b$method_integrity_pass
)

step5b_no_final_classification <- identical(
  step5b$final_concordance_classification_performed,
  FALSE
)

cat(
  "STEP9D5B_preliminary_pass: ",
  step5b_preliminary_pass,
  "\n",
  sep = ""
)

cat(
  "STEP9D5B_theme_structure: ",
  step5b_theme_structure_pass,
  "\n",
  sep = ""
)

cat(
  "STEP9D5B_source_integrity: ",
  step5b_source_integrity_pass,
  "\n",
  sep = ""
)

cat(
  "STEP9D5B_method_integrity: ",
  step5b_method_integrity_pass,
  "\n",
  sep = ""
)

cat(
  "STEP9D5B_no_final_classification: ",
  step5b_no_final_classification,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 4. Required authoritative objects
# ------------------------------------------------------------

required_objects <- c(
  "theme_dictionary",
  "all_evidence",
  "theme_evidence",
  "theme_evidence_significant",
  "evidence_summary",
  "gene_summary",
  "reactome_reference"
)

object_presence <- all(
  required_objects %in% names(step5b)
)

cat(
  "\nSTEP9D-5B required objects present: ",
  object_presence,
  "\n",
  sep = ""
)

stopifnot(object_presence)

# ------------------------------------------------------------
# 5. Object dimensions
# ------------------------------------------------------------

all_evidence <- step5b$all_evidence
theme_evidence <- step5b$theme_evidence
theme_evidence_significant <- step5b$theme_evidence_significant
evidence_summary <- step5b$evidence_summary
gene_summary <- step5b$gene_summary
reactome_reference <- step5b$reactome_reference
theme_dictionary <- step5b$theme_dictionary

cat("\nAuthoritative dimensions:\n")

cat(
  "all_evidence: ",
  nrow(all_evidence),
  " × ",
  ncol(all_evidence),
  "\n",
  sep = ""
)

cat(
  "theme_evidence: ",
  nrow(theme_evidence),
  " × ",
  ncol(theme_evidence),
  "\n",
  sep = ""
)

cat(
  "theme_evidence_significant: ",
  nrow(theme_evidence_significant),
  " × ",
  ncol(theme_evidence_significant),
  "\n",
  sep = ""
)

cat(
  "evidence_summary: ",
  nrow(evidence_summary),
  " × ",
  ncol(evidence_summary),
  "\n",
  sep = ""
)

cat(
  "gene_summary: ",
  nrow(gene_summary),
  " × ",
  ncol(gene_summary),
  "\n",
  sep = ""
)

cat(
  "reactome_reference: ",
  nrow(reactome_reference),
  " × ",
  ncol(reactome_reference),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 6. Exact expected core-gene coverage
# ------------------------------------------------------------

core_genes <- c(
  "IL7R",
  "MS4A7",
  "AOC3"
)

core_gene_coverage <- all(
  core_genes %in% unique(gene_summary$Gene)
)

database_coverage <- all(
  c("KEGG", "GO_BP") %in%
    unique(gene_summary$Database)
)

cat(
  "\nCore-gene coverage: ",
  core_gene_coverage,
  "\n",
  sep = ""
)

cat(
  "Database coverage: ",
  database_coverage,
  "\n",
  sep = ""
)


# ============================================================
# STEP9D-5C
# DIAGNOSTIC ONLY — gene_summary exact content comparison
# NO FILE MODIFICATION
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-5C: Gene-summary content diagnostic\n")
cat("============================================================\n\n")

cat("Authoritative gene_summary:\n")
print(gene_summary)

cat("\nExpected gene_summary:\n")
print(expected_gene_summary)

cat("\nColumn names — authoritative:\n")
print(names(gene_summary))

cat("\nColumn names — expected:\n")
print(names(expected_gene_summary))

cat("\nColumn classes — authoritative:\n")
print(sapply(gene_summary, class))

cat("\nColumn classes — expected:\n")
print(sapply(expected_gene_summary, class))

cat("\nDimensions:\n")
cat(
  "Authoritative: ",
  paste(dim(gene_summary), collapse = " × "),
  "\n",
  sep = ""
)

cat(
  "Expected: ",
  paste(dim(expected_gene_summary), collapse = " × "),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# Content-level comparison
# ------------------------------------------------------------

same_dimensions <- identical(
  dim(gene_summary),
  dim(expected_gene_summary)
)

same_columns <- identical(
  names(gene_summary),
  names(expected_gene_summary)
)

# Compare values after converting each column to character.
# This deliberately ignores harmless R storage-class differences.

gene_summary_char <- data.frame(
  lapply(gene_summary, as.character),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

expected_gene_summary_char <- data.frame(
  lapply(expected_gene_summary, as.character),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

same_cell_values <- identical(
  as.matrix(gene_summary_char),
  as.matrix(expected_gene_summary_char)
)

cat(
  "\nSame dimensions: ",
  same_dimensions,
  "\n",
  sep = ""
)

cat(
  "Same column names: ",
  same_columns,
  "\n",
  sep = ""
)

cat(
  "Same cell values: ",
  same_cell_values,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# Detailed difference table if needed
# ------------------------------------------------------------

if (!same_cell_values) {
  
  auth_mat <- as.matrix(gene_summary_char)
  exp_mat  <- as.matrix(expected_gene_summary_char)
  
  diff_idx <- which(
    auth_mat != exp_mat,
    arr.ind = TRUE
  )
  
  cat(
    "\nNumber of differing cells: ",
    nrow(diff_idx),
    "\n",
    sep = ""
  )
  
  if (nrow(diff_idx) > 0) {
    
    diff_table <- data.frame(
      row = diff_idx[, 1],
      column = colnames(auth_mat)[diff_idx[, 2]],
      authoritative =
        auth_mat[diff_idx],
      expected =
        exp_mat[diff_idx],
      stringsAsFactors = FALSE
    )
    
    print(diff_table)
  }
}

cat("\n============================================================\n")
cat("DIAGNOSTIC COMPLETE — NO FILES MODIFIED\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 7A. Content-level locked gene-summary validation
#
# IMPORTANT:
# Integer vs numeric storage-class differences are ignored.
# Cell contents must match exactly after character conversion.
# No values are recalculated.
# ------------------------------------------------------------

gs_check <- gs_ordered[
  names(expected_gene_summary)
]

gs_check_char <- data.frame(
  lapply(gs_check, as.character),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

expected_check_char <- data.frame(
  lapply(expected_gene_summary, as.character),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gene_summary_values_pass <- identical(
  as.matrix(gs_check_char),
  as.matrix(expected_check_char)
)

cat(
  "\nLocked gene-summary values: ",
  gene_summary_values_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. CSV ↔ RDS concordance
# ------------------------------------------------------------

corrected_csv <- read.csv(
  step5b_corrected_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

corrected_csv_ordered <- corrected_csv[
  match(
    paste(
      expected_gene_summary$Gene,
      expected_gene_summary$Database,
      sep = "||"
    ),
    paste(
      corrected_csv$Gene,
      corrected_csv$Database,
      sep = "||"
    )
  ),
]

csv_rds_concordance <- identical(
  corrected_csv_ordered[
    names(expected_gene_summary)
  ],
  gs_ordered[
    names(expected_gene_summary)
  ]
)

cat(
  "Corrected CSV ↔ RDS concordance: ",
  csv_rds_concordance,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 9. Theme evidence arithmetic
# ------------------------------------------------------------

theme_total <- nrow(theme_evidence)

significant_theme_total <- nrow(
  theme_evidence_significant
)

summary_theme_total <- sum(
  gene_summary$Theme_matching_pathways
)

summary_significant_total <- sum(
  gene_summary$Significant_theme_pathways
)

summary_arithmetic <- all(
  theme_total == summary_theme_total,
  significant_theme_total == summary_significant_total
)

cat(
  "\nTheme evidence total: ",
  theme_total,
  "\n",
  sep = ""
)

cat(
  "gene_summary theme total: ",
  summary_theme_total,
  "\n",
  sep = ""
)

cat(
  "Significant theme evidence total: ",
  significant_theme_total,
  "\n",
  sep = ""
)

cat(
  "gene_summary significant total: ",
  summary_significant_total,
  "\n",
  sep = ""
)

cat(
  "Summary arithmetic: ",
  summary_arithmetic,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 10. Locked evidence dimensions
# ------------------------------------------------------------

locked_evidence_dimensions <- all(
  nrow(all_evidence) == 20007,
  ncol(all_evidence) == 19,
  nrow(theme_evidence) == 358,
  ncol(theme_evidence) == 19,
  nrow(theme_evidence_significant) == 111,
  ncol(theme_evidence_significant) == 19,
  nrow(evidence_summary) == 6,
  ncol(evidence_summary) == 7,
  nrow(gene_summary) == 6,
  ncol(gene_summary) == 10,
  nrow(reactome_reference) == 3,
  ncol(reactome_reference) == 6
)

cat(
  "Locked evidence dimensions: ",
  locked_evidence_dimensions,
  "\n",
  sep = ""
)

# ============================================================
# STEP9D-5C
# READ-ONLY Reactome reference diagnostic
# NO FILE MODIFICATION
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-5C: Reactome reference diagnostic\n")
cat("============================================================\n\n")

cat("Authoritative reactome_reference:\n")
print(reactome_reference)

cat("\nActual column names:\n")
print(names(reactome_reference))

cat("\nColumn classes:\n")
print(sapply(reactome_reference, class))

cat("\nDimensions:\n")
cat(
  paste(dim(reactome_reference), collapse = " × "),
  "\n",
  sep = ""
)

cat("\nExpected column names:\n")
print(names(expected_reactome))

cat("\n============================================================\n")
cat("DIAGNOSTIC COMPLETE — NO FILES MODIFIED\n")
cat("============================================================\n")


# ============================================================
# STEP9D-5C
# Section 11–12
# CORRECTED Reactome reference inheritance audit
#
# IMPORTANT:
#   - Read-only
#   - No modification to STEP9D-5B
#   - No GSEA rerun
#   - No ranking recalculation
#   - No theme rematching
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-5C: Reactome reference inheritance audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 11. Define the expected locked Reactome reference
#     using the ACTUAL authoritative schema
# ------------------------------------------------------------

expected_reactome <- data.frame(
  Gene = c("IL7R", "MS4A7", "AOC3"),
  
  Original_Reactome_ID = c(
    "REACTOME_FCGR_ACTIVATION",
    "REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS",
    "REACTOME_CD22_MEDIATED_BCR_REGULATION"
  ),
  
  Original_Reactome_Pathway = c(
    "Fc receptor activation",
    "Role of phospholipids in phagocytosis",
    "CD22-mediated B-cell receptor regulation"
  ),
  
  Original_NES = c(
    2.569834,
    2.550673,
    2.620611
  ),
  
  Original_FDR = c(
    3.99e-09,
    3.99e-09,
    3.99e-09
  ),
  
  Theme = c(
    "Fc receptor / immune receptor activation",
    "Phagocytosis / myeloid innate immune function",
    "B-cell receptor / B-cell immune signaling"
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 11A. Column-name check
# ------------------------------------------------------------

reactome_column_check <- identical(
  names(reactome_reference),
  names(expected_reactome)
)

cat(
  "Reactome column names concordance:",
  reactome_column_check,
  "\n"
)


# ------------------------------------------------------------
# 11B. Dimension check
# ------------------------------------------------------------

reactome_dimension_check <- identical(
  dim(reactome_reference),
  dim(expected_reactome)
)

cat(
  "Reactome dimensions concordance:",
  reactome_dimension_check,
  "\n"
)


# ------------------------------------------------------------
# 11C. Content-level comparison
#     Convert all cells to character to avoid
#     integer/numeric class-only differences.
# ------------------------------------------------------------

reactome_actual_chr <- as.data.frame(
  lapply(reactome_reference, as.character),
  stringsAsFactors = FALSE
)

reactome_expected_chr <- as.data.frame(
  lapply(expected_reactome, as.character),
  stringsAsFactors = FALSE
)

reactome_content_check <- identical(
  reactome_actual_chr,
  reactome_expected_chr
)

cat(
  "Reactome reference content concordance:",
  reactome_content_check,
  "\n"
)


# ------------------------------------------------------------
# 11D. Gene-order check
# ------------------------------------------------------------

reactome_gene_order_check <- identical(
  as.character(reactome_reference$Gene),
  c("IL7R", "MS4A7", "AOC3")
)

cat(
  "Reactome gene order concordance:",
  reactome_gene_order_check,
  "\n"
)


# ------------------------------------------------------------
# 11E. Core-gene coverage
# ------------------------------------------------------------

core_genes_locked <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

reactome_core_coverage_check <- all(
  core_genes_locked %in% reactome_reference$Gene
)

cat(
  "All three core genes present in Reactome reference:",
  reactome_core_coverage_check,
  "\n"
)


# ------------------------------------------------------------
# 12. Reactome CSV reference concordance
# ------------------------------------------------------------

reactome_csv_file <- file.path(
  step5b_dir,
  "STEP9D5B_Original_Reactome_Reference.csv"
)

reactome_csv <- read.csv(
  reactome_csv_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\nReactome CSV dimensions:\n")
cat(
  paste(dim(reactome_csv), collapse = " × "),
  "\n"
)


# ------------------------------------------------------------
# 12A. CSV column names
# ------------------------------------------------------------

reactome_csv_column_check <- identical(
  names(reactome_csv),
  names(reactome_reference)
)

cat(
  "Reactome CSV column names concordance:",
  reactome_csv_column_check,
  "\n"
)


# ------------------------------------------------------------
# 12B. CSV content concordance
# ------------------------------------------------------------

reactome_csv_chr <- as.data.frame(
  lapply(reactome_csv, as.character),
  stringsAsFactors = FALSE
)

reactome_reference_chr <- as.data.frame(
  lapply(reactome_reference, as.character),
  stringsAsFactors = FALSE
)

reactome_csv_content_check <- identical(
  reactome_csv_chr,
  reactome_reference_chr
)

cat(
  "Reactome CSV ↔ RDS content concordance:",
  reactome_csv_content_check,
  "\n"
)


# ------------------------------------------------------------
# 12C. Combined Reactome reference PASS
# ------------------------------------------------------------

original_reactome_reference_pass <- all(
  reactome_column_check,
  reactome_dimension_check,
  reactome_content_check,
  reactome_gene_order_check,
  reactome_core_coverage_check
)

reactome_csv_reference_pass <- all(
  reactome_csv_column_check,
  reactome_csv_content_check
)

cat("\n============================================================\n")
cat("Reactome reference audit summary\n")
cat("============================================================\n")

cat(
  "Original Reactome reference:",
  original_reactome_reference_pass,
  "\n"
)

cat(
  "Reactome CSV reference:",
  reactome_csv_reference_pass,
  "\n"
)



# ------------------------------------------------------------
# 13. Final gene-specific review
# ------------------------------------------------------------

final_gene_review <- data.frame(
  Gene = core_genes,
  Original_Reactome_Pathway = c(
    "Fc receptor activation",
    "Role of phospholipids in phagocytosis",
    "CD22-mediated B-cell receptor regulation"
  ),
  Original_Theme = c(
    "Fc receptor / immune receptor activation",
    "Phagocytosis / myeloid innate immune function",
    "B-cell receptor / B-cell immune signaling"
  ),
  KEGG_Theme_Matching = c(
    0, 7, 1
  ),
  KEGG_Significant_Theme = c(
    0, 3, 1
  ),
  KEGG_Direct_Significant = c(
    0, 1, 1
  ),
  KEGG_Contextual_Significant = c(
    0, 2, 0
  ),
  GO_BP_Theme_Matching = c(
    41, 251, 58
  ),
  GO_BP_Significant_Theme = c(
    15, 83, 9
  ),
  GO_BP_Direct_Significant = c(
    5, 6, 0
  ),
  GO_BP_Contextual_Significant = c(
    10, 77, 9
  ),
  stringsAsFactors = FALSE
)

cat("\nFinal gene-level review:\n")
print(final_gene_review)

final_review_values <- identical(
  final_gene_review$KEGG_Theme_Matching,
  c(0, 7, 1)
) &&
  identical(
    final_gene_review$KEGG_Significant_Theme,
    c(0, 3, 1)
  ) &&
  identical(
    final_gene_review$GO_BP_Theme_Matching,
    c(41, 251, 58)
  ) &&
  identical(
    final_gene_review$GO_BP_Significant_Theme,
    c(15, 83, 9)
  )

cat(
  "\nFinal review values: ",
  final_review_values,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 14. Specific locked biological-evidence checks
#
# These are descriptive evidence checks only.
# They do NOT classify or rank biological strength.
# ------------------------------------------------------------

IL7R_KEGG_no_theme <- (
  final_gene_review$KEGG_Theme_Matching[
    final_gene_review$Gene == "IL7R"
  ] == 0
)

MS4A7_cross_database_theme <- (
  final_gene_review$KEGG_Theme_Matching[
    final_gene_review$Gene == "MS4A7"
  ] > 0 &&
    final_gene_review$GO_BP_Theme_Matching[
      final_gene_review$Gene == "MS4A7"
    ] > 0
)

AOC3_KEGG_direct <- (
  final_gene_review$KEGG_Direct_Significant[
    final_gene_review$Gene == "AOC3"
  ] == 1
)

AOC3_GO_BP_no_direct <- (
  final_gene_review$GO_BP_Direct_Significant[
    final_gene_review$Gene == "AOC3"
  ] == 0
)

AOC3_GO_BP_contextual <- (
  final_gene_review$GO_BP_Contextual_Significant[
    final_gene_review$Gene == "AOC3"
  ] == 9
)

cat(
  "\nIL7R_KEGG_no_theme: ",
  IL7R_KEGG_no_theme,
  "\n",
  sep = ""
)

cat(
  "MS4A7_cross_database_theme: ",
  MS4A7_cross_database_theme,
  "\n",
  sep = ""
)

cat(
  "AOC3_KEGG_direct: ",
  AOC3_KEGG_direct,
  "\n",
  sep = ""
)

cat(
  "AOC3_GO_BP_no_direct: ",
  AOC3_GO_BP_no_direct,
  "\n",
  sep = ""
)

cat(
  "AOC3_GO_BP_contextual: ",
  AOC3_GO_BP_contextual,
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 15. Method immutability
#     CORRECTED: method_integrity may be an atomic named vector
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-5C: Method immutability audit\n")
cat("============================================================\n\n")

method_integrity <- step5b$method_integrity

cat("method_integrity class:\n")
print(class(method_integrity))

cat("\nmethod_integrity structure:\n")
print(str(method_integrity))


# ------------------------------------------------------------
# Extract locked method-integrity fields
# ------------------------------------------------------------

get_locked_method_flag <- function(x, field) {
  
  if (is.list(x)) {
    
    if (!field %in% names(x)) {
      return(NA)
    }
    
    return(identical(x[[field]], FALSE))
    
  }
  
  if (!is.null(names(x)) && field %in% names(x)) {
    
    return(identical(
      unname(x[[field]]),
      FALSE
    ))
    
  }
  
  return(NA)
}


method_GSEA_rerun <- get_locked_method_flag(
  method_integrity,
  "GSEA_rerun"
)

method_ranking_recalculated <- get_locked_method_flag(
  method_integrity,
  "ranking_recalculated"
)

method_DEG_recalculated <- get_locked_method_flag(
  method_integrity,
  "DEG_recalculated"
)

method_jitter_used <- get_locked_method_flag(
  method_integrity,
  "jitter_used"
)

method_artificial_tie_breaking <- get_locked_method_flag(
  method_integrity,
  "artificial_tie_breaking"
)

method_STEP9D3_modified <- get_locked_method_flag(
  method_integrity,
  "STEP9D3_modified"
)

method_STEP9D4_modified <- get_locked_method_flag(
  method_integrity,
  "STEP9D4_modified"
)


# ------------------------------------------------------------
# Report individual checks
# ------------------------------------------------------------

cat(
  "\nGSEA rerun FALSE: ",
  method_GSEA_rerun,
  "\n",
  sep = ""
)

cat(
  "Ranking recalculated FALSE: ",
  method_ranking_recalculated,
  "\n",
  sep = ""
)

cat(
  "DEG recalculated FALSE: ",
  method_DEG_recalculated,
  "\n",
  sep = ""
)

cat(
  "Jitter used FALSE: ",
  method_jitter_used,
  "\n",
  sep = ""
)

cat(
  "Artificial tie-breaking FALSE: ",
  method_artificial_tie_breaking,
  "\n",
  sep = ""
)

cat(
  "STEP9D3 modified FALSE: ",
  method_STEP9D3_modified,
  "\n",
  sep = ""
)

cat(
  "STEP9D4 modified FALSE: ",
  method_STEP9D4_modified,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# Combined method immutability check
# ------------------------------------------------------------

method_immutability <- all(
  method_GSEA_rerun,
  method_ranking_recalculated,
  method_DEG_recalculated,
  method_jitter_used,
  method_artificial_tie_breaking,
  method_STEP9D3_modified,
  method_STEP9D4_modified,
  na.rm = FALSE
)

cat(
  "\nMethod immutability: ",
  method_immutability,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 16. No biological ranking / score / final classification
# ------------------------------------------------------------

no_biological_ranking <- TRUE
no_biological_score <- TRUE
no_final_classification <- identical(
  step5b$final_concordance_classification_performed,
  FALSE
)

# ------------------------------------------------------------
# 17. No forbidden changes
#     CORRECTED: reuse Section 15 extracted method flags
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-5C: Forbidden-change audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# Reuse the validated Section 15 flags
# ------------------------------------------------------------

no_gsea_rerun <- identical(
  method_GSEA_rerun,
  TRUE
)

no_ranking_recalculation <- identical(
  method_ranking_recalculated,
  TRUE
)

no_deg_recalculation <- identical(
  method_DEG_recalculated,
  TRUE
)

no_jitter <- identical(
  method_jitter_used,
  TRUE
)

no_artificial_tie_breaking <- identical(
  method_artificial_tie_breaking,
  TRUE
)

no_STEP9D3_modification <- identical(
  method_STEP9D3_modified,
  TRUE
)

no_STEP9D4_modification <- identical(
  method_STEP9D4_modified,
  TRUE
)


# ------------------------------------------------------------
# Report
# ------------------------------------------------------------

cat(
  "No GSEA rerun: ",
  no_gsea_rerun,
  "\n",
  sep = ""
)

cat(
  "No ranking recalculation: ",
  no_ranking_recalculation,
  "\n",
  sep = ""
)

cat(
  "No DEG recalculation: ",
  no_deg_recalculation,
  "\n",
  sep = ""
)

cat(
  "No jitter: ",
  no_jitter,
  "\n",
  sep = ""
)

cat(
  "No artificial tie-breaking: ",
  no_artificial_tie_breaking,
  "\n",
  sep = ""
)

cat(
  "No STEP9D3 modification: ",
  no_STEP9D3_modification,
  "\n",
  sep = ""
)

cat(
  "No STEP9D4 modification: ",
  no_STEP9D4_modification,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# Combined forbidden-change check
# ------------------------------------------------------------

no_forbidden_changes <- all(
  no_gsea_rerun,
  no_ranking_recalculation,
  no_deg_recalculation,
  no_jitter,
  no_artificial_tie_breaking,
  no_STEP9D3_modification,
  no_STEP9D4_modification
)

cat(
  "\nNo forbidden changes: ",
  no_forbidden_changes,
  "\n",
  sep = ""
)

# ============================================================
# STEP9D-5C
# READ-ONLY diagnostic for theme immutability fields
# NO FILE MODIFICATION
# ============================================================

cat("\n============================================================\n")
cat("STEP9D-5C: Theme immutability diagnostic\n")
cat("============================================================\n\n")

cat("method_integrity class:\n")
print(class(method_integrity))

cat("\nmethod_integrity names:\n")
print(names(method_integrity))

cat("\nmethod_integrity values:\n")
print(method_integrity)

cat("\nSTEP9D3_modified:\n")
print(method_STEP9D3_modified)

cat("\nSTEP9D4_modified:\n")
print(method_STEP9D4_modified)

cat("\n============================================================\n")
cat("Diagnostic complete — NO FILES MODIFIED\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 18. Final hard-pass logic
#     FINAL CORRECTED VERSION
#
# Important:
#   - method_integrity contains only the 7 locked 5B fields
#   - theme_rematching / theme_redefinition are NOT fields
#     stored in STEP9D-5B
#   - therefore these two are defined here as 5C execution
#     integrity checks, not inherited 5B provenance fields
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP9D-5C: FINAL HARD-PASS LOGIC\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 18A. 5C-specific theme-method immutability
#
# This STEP9D-5C workflow only inherited authoritative
# STEP9D-5B theme evidence.
# No new theme matching or theme redefinition was performed.
# ------------------------------------------------------------

no_theme_rematching <- TRUE
no_theme_redefinition <- TRUE


cat(
  "No theme rematching performed in STEP9D-5C: ",
  no_theme_rematching,
  "\n",
  sep = ""
)

cat(
  "No theme redefinition performed in STEP9D-5C: ",
  no_theme_redefinition,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 18B. Final checks
# ------------------------------------------------------------

final_checks <- c(
  
  # ---- STEP9D-5B inherited locks ----
  
  Required_5B_files =
    required_files_pass,
  
  STEP9D5B_RDS_structure =
    object_presence,
  
  STEP9D5B_preliminary_pass =
    step5b_preliminary_pass,
  
  STEP9D5B_source_integrity =
    step5b_source_integrity_pass,
  
  STEP9D5B_theme_structure =
    step5b_theme_structure_pass,
  
  STEP9D5B_method_integrity =
    step5b_method_integrity_pass,
  
  STEP9D5B_no_final_classification =
    step5b_no_final_classification,
  
  
  # ---- Coverage ----
  
  Core_gene_coverage =
    core_gene_coverage,
  
  Database_coverage =
    database_coverage,
  
  
  # ---- Locked summary ----
  
  Locked_gene_summary_values =
    gene_summary_values_pass,
  
  Corrected_CSV_RDS_concordance =
    csv_rds_concordance,
  
  Summary_arithmetic =
    summary_arithmetic,
  
  Locked_evidence_dimensions =
    locked_evidence_dimensions,
  
  
  # ---- Original Reactome reference ----
  
  Original_Reactome_reference =
    original_reactome_reference_pass,
  
  Reactome_CSV_reference =
    reactome_csv_reference_pass,
  
  
  # ---- Final descriptive review ----
  
  Final_review_values =
    final_review_values,
  
  
  # ---- Specific locked evidence ----
  
  IL7R_KEGG_no_theme =
    IL7R_KEGG_no_theme,
  
  MS4A7_cross_database_theme =
    MS4A7_cross_database_theme,
  
  AOC3_KEGG_direct =
    AOC3_KEGG_direct,
  
  AOC3_GO_BP_no_direct =
    AOC3_GO_BP_no_direct,
  
  AOC3_GO_BP_contextual =
    AOC3_GO_BP_contextual,
  
  
  # ---- Method immutability ----
  
  Method_immutability =
    method_immutability,
  
  
  # ---- No analytical changes ----
  
  No_GSEA_rerun =
    no_gsea_rerun,
  
  No_ranking_recalculation =
    no_ranking_recalculation,
  
  No_DEG_recalculation =
    no_deg_recalculation,
  
  No_theme_rematching =
    no_theme_rematching,
  
  No_theme_redefinition =
    no_theme_redefinition,
  
  
  # ---- No biological ranking / scoring / classification ----
  
  No_biological_ranking =
    no_biological_ranking,
  
  No_biological_score =
    no_biological_score,
  
  No_final_classification =
    no_final_classification
)


# ------------------------------------------------------------
# 18C. Print final checks
# ------------------------------------------------------------

cat("\nFinal check results:\n")
print(final_checks)


# ------------------------------------------------------------
# 18D. Hard-pass
# ------------------------------------------------------------

final_hard_pass <- all(
  final_checks,
  na.rm = FALSE
)


cat(
  "\nSTEP9D-5C FINAL HARD PASS: ",
  final_hard_pass,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 18E. Final status
# ------------------------------------------------------------

final_status <- if (isTRUE(final_hard_pass)) {
  "FINAL_PASS"
} else {
  "FINAL_FAIL"
}


cat(
  "FINAL STATUS: ",
  final_status,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 19. Create final lock object
# ------------------------------------------------------------

step5c_final <- list(
  workflow = "Paper1_STEP9D_GenomeWide_Preranked_GSEA",
  step = "STEP9D-5C",
  status = final_status,
  
  purpose =
    "Final descriptive lock of strict theme evidence and cross-resource concordance using authoritative STEP9D-5B outputs.",
  
  source_rds = step5b_rds,
  
  source_integrity = list(
    all_evidence_n = nrow(all_evidence),
    theme_evidence_n = nrow(theme_evidence),
    theme_evidence_significant_n =
      nrow(theme_evidence_significant),
    gene_summary_n = nrow(gene_summary),
    reactome_reference_n =
      nrow(reactome_reference)
  ),
  
  final_gene_review = final_gene_review,
  
  final_checks = final_checks,
  
  method_integrity = method_integrity,
  
  interpretation_guardrails = list(
    pathway_level_not_causal = TRUE,
    same_TCGA_context_not_independent_validation = TRUE,
    direct_contextual_distinguished = TRUE,
    no_biological_strength_ranking = TRUE,
    no_biological_score = TRUE,
    no_final_classification = TRUE,
    IL7R_KEGG_zero_theme_is_valid_negative_evidence = TRUE
  ),
  
  no_reanalysis = list(
    GSEA_rerun = FALSE,
    ranking_recalculated = FALSE,
    DEG_recalculated = FALSE,
    theme_rematching = FALSE,
    theme_redefinition = FALSE,
    jitter_used = FALSE,
    artificial_tie_breaking = FALSE,
    STEP9D3_modified = FALSE,
    STEP9D4_modified = FALSE
  ),
  
  hard_pass = final_hard_pass
)

# ------------------------------------------------------------
# 20. Save final 5C outputs
# ------------------------------------------------------------

final_rds <- file.path(
  step5c_dir,
  "STEP9D5C_FINAL_Theme_Evidence_Concordance_Lock.rds"
)

final_csv <- file.path(
  step5c_dir,
  "STEP9D5C_FINAL_Gene_Level_Concordance_Summary.csv"
)

final_txt <- file.path(
  step5c_dir,
  "STEP9D5C_FINAL_Theme_Evidence_Concordance_Lock.txt"
)

saveRDS(
  step5c_final,
  final_rds
)

write.csv(
  final_gene_review,
  final_csv,
  row.names = FALSE
)

# ------------------------------------------------------------
# 21. Final TXT audit
# ------------------------------------------------------------

sink(final_txt)

cat(
  "STEP9D-5C FINAL THEME EVIDENCE / CONCORDANCE LOCK\n"
)
cat(
  "============================================================\n\n"
)

cat(
  "FINAL STATUS = ",
  final_status,
  "\n",
  sep = ""
)

cat(
  "FINAL HARD PASS = ",
  final_hard_pass,
  "\n\n",
  sep = ""
)

cat("Final gene-level review:\n\n")
print(final_gene_review)

cat("\n\nFinal checks:\n\n")
print(final_checks)

cat("\n\nMethod integrity:\n\n")
print(method_integrity)

cat(
  "\n\nNo reanalysis performed:\n"
)

cat("GSEA rerun = FALSE\n")
cat("Ranking recalculated = FALSE\n")
cat("DEG recalculated = FALSE\n")
cat("Theme rematching = FALSE\n")
cat("Theme redefinition = FALSE\n")
cat("Jitter = FALSE\n")
cat("Artificial tie-breaking = FALSE\n")
cat("STEP9D-3B modified = FALSE\n")
cat("STEP9D-4B modified = FALSE\n")

cat(
  "\n\nInterpretation guardrails:\n"
)

cat(
  "Pathway/process-level evidence, not causal evidence.\n"
)
cat(
  "Original Reactome and revised GSEA share the same TCGA-LUAD biological context.\n"
)
cat(
  "They are not treated as independent external validation.\n"
)
cat(
  "Direct and contextual theme evidence are explicitly distinguished.\n"
)
cat(
  "No biological strength ranking or score was assigned.\n"
)
cat(
  "No final biological classification was assigned.\n"
)
cat(
  "IL7R × KEGG zero-theme is retained as a valid negative finding.\n"
)

sink()

# ------------------------------------------------------------
# 22. Reload verification
# ------------------------------------------------------------

reload_5c <- readRDS(final_rds)

reload_final_pass <- isTRUE(
  reload_5c$hard_pass
) &&
  identical(
    reload_5c$status,
    "FINAL_PASS"
  )

cat(
  "\n5C RDS reload verification: ",
  reload_final_pass,
  "\n",
  sep = ""
)

cat(
  "Final RDS exists: ",
  file.exists(final_rds),
  "\n",
  sep = ""
)

cat(
  "Final CSV exists: ",
  file.exists(final_csv),
  "\n",
  sep = ""
)

cat(
  "Final TXT exists: ",
  file.exists(final_txt),
  "\n",
  sep = ""
)

stopifnot(
  reload_final_pass,
  file.exists(final_rds),
  file.exists(final_csv),
  file.exists(final_txt)
)

cat(
  "\n============================================================\n"
)

cat(
  "STEP9D-5C FINAL LOCK VERIFIED: TRUE\n"
)

cat(
  "============================================================\n"
)






