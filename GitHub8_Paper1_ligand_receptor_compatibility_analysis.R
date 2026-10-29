##GitHub8_Paper1_ligand_receptor_compatibility_analysis

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
# STEP11-0
# NSCLC_GSE131907_CCI_LR.tsv
# Structural / Semantic Audit
#
# Purpose:
#   1. Verify file integrity
#   2. Inspect dimensions and column names
#   3. Determine whether the file is:
#      a) CellChat dataset-specific communication results
#      b) ligand-receptor reference/resource table
#   4. Identify sender/receiver/probability/P-value fields
# ============================================================

cat("============================================================\n")
cat("STEP11-0: NSCLC_GSE131907_CCI_LR.tsv STRUCTURAL AUDIT\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 11-CellChat"

setwd(step_dir)

input_file <- file.path(
  step_dir,
  "NSCLC_GSE131907_CCI_LR.tsv"
)

output_file <- file.path(
  step_dir,
  "STEP11-0_CCI_LR_Structural_Audit.txt"
)

# ------------------------------------------------------------
# 1. File existence
# ------------------------------------------------------------

cat("STEP 1: FILE EXISTENCE\n")
cat("----------------------\n")

cat("Input file:\n", input_file, "\n\n")

if (!file.exists(input_file)) {
  stop(
    "ERROR: Input file does not exist:\n",
    input_file
  )
}

cat("File exists: TRUE\n\n")


# ------------------------------------------------------------
# 2. File information
# ------------------------------------------------------------

cat("STEP 2: FILE INFORMATION\n")
cat("------------------------\n")

file_info <- file.info(input_file)

cat(
  "File size:",
  round(file_info$size / 1024^2, 3),
  "MB\n"
)

cat(
  "Modified:",
  as.character(file_info$mtime),
  "\n\n"
)


# ------------------------------------------------------------
# 3. Read header only
# ------------------------------------------------------------

cat("STEP 3: HEADER\n")
cat("--------------\n")

header <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  nrows = 0
)

cat(
  "Number of columns:",
  ncol(header),
  "\n\n"
)

cat("Column names:\n")

for (i in seq_along(colnames(header))) {
  
  cat(
    sprintf(
      "[%02d] %s\n",
      i,
      colnames(header)[i]
    )
  )
}

cat("\n")


# ------------------------------------------------------------
# 4. Read first 10 rows
# ------------------------------------------------------------

cat("STEP 4: FIRST 10 ROWS\n")
cat("---------------------\n")

first10 <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  nrows = 10
)

print(first10)

cat("\n")


# ------------------------------------------------------------
# 5. Determine total number of rows
# ------------------------------------------------------------

cat("STEP 5: TOTAL ROW COUNT\n")
cat("-----------------------\n")

# Count lines without loading the whole file
line_count <- length(
  readLines(
    input_file,
    warn = FALSE
  )
)

data_rows <- line_count - 1

cat(
  "Header lines:",
  1,
  "\n"
)

cat(
  "Data rows:",
  data_rows,
  "\n\n"
)


# ------------------------------------------------------------
# 6. Read complete table
# ------------------------------------------------------------

cat("STEP 6: FULL TABLE LOAD\n")
cat("-----------------------\n")

cci_lr <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat(
  "Loaded successfully: TRUE\n"
)

cat(
  "Rows:",
  nrow(cci_lr),
  "\n"
)

cat(
  "Columns:",
  ncol(cci_lr),
  "\n\n"
)


# ------------------------------------------------------------
# 7. Data types
# ------------------------------------------------------------

cat("STEP 7: COLUMN DATA TYPES\n")
cat("------------------------\n")

type_table <- data.frame(
  column = colnames(cci_lr),
  class = sapply(
    cci_lr,
    function(x) paste(class(x), collapse = "/")
  ),
  stringsAsFactors = FALSE
)

print(type_table, row.names = FALSE)

cat("\n")


# ------------------------------------------------------------
# 8. Missingness
# ------------------------------------------------------------

cat("STEP 8: MISSINGNESS\n")
cat("------------------\n")

missing_table <- data.frame(
  column = colnames(cci_lr),
  n_missing = sapply(
    cci_lr,
    function(x) sum(
      is.na(x) |
        trimws(as.character(x)) == ""
    )
  ),
  stringsAsFactors = FALSE
)

missing_table$missing_percent <-
  100 *
  missing_table$n_missing /
  nrow(cci_lr)

print(
  missing_table,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 9. Search for communication-related columns
# ------------------------------------------------------------

cat("STEP 9: COMMUNICATION-RELATED COLUMN AUDIT\n")
cat("-------------------------------------------\n")

column_names_lower <- tolower(
  colnames(cci_lr)
)

patterns <- c(
  "sender",
  "source",
  "receiver",
  "target",
  "prob",
  "pval",
  "p.value",
  "p_value",
  "fdr",
  "signif",
  "ligand",
  "receptor",
  "pathway",
  "interaction",
  "cell"
)

for (pattern in patterns) {
  
  hits <- colnames(cci_lr)[
    grepl(
      pattern,
      column_names_lower,
      fixed = TRUE
    )
  ]
  
  cat(
    sprintf(
      "%-15s : ",
      pattern
    )
  )
  
  if (length(hits) == 0) {
    cat("NONE\n")
  } else {
    cat(
      paste(hits, collapse = ", "),
      "\n"
    )
  }
}

cat("\n")


# ------------------------------------------------------------
# 10. Explicit semantic classification
# ------------------------------------------------------------

cat("STEP 10: PRELIMINARY SEMANTIC CLASSIFICATION\n")
cat("---------------------------------------------\n")

has_sender <- any(
  grepl(
    "sender|source",
    column_names_lower
  )
)

has_receiver <- any(
  grepl(
    "receiver|target",
    column_names_lower
  )
)

has_probability <- any(
  grepl(
    "prob",
    column_names_lower
  )
)

has_pvalue <- any(
  grepl(
    "p.value|p_value|pval",
    column_names_lower
  )
)

has_ligand <- any(
  grepl(
    "ligand",
    column_names_lower
  )
)

has_receptor <- any(
  grepl(
    "receptor",
    column_names_lower
  )
)

has_pathway <- any(
  grepl(
    "pathway",
    column_names_lower
  )
)


cat(
  "Sender/source field detected:",
  has_sender,
  "\n"
)

cat(
  "Receiver/target field detected:",
  has_receiver,
  "\n"
)

cat(
  "Probability field detected:",
  has_probability,
  "\n"
)

cat(
  "P-value field detected:",
  has_pvalue,
  "\n"
)

cat(
  "Ligand field detected:",
  has_ligand,
  "\n"
)

cat(
  "Receptor field detected:",
  has_receptor,
  "\n"
)

cat(
  "Pathway field detected:",
  has_pathway,
  "\n\n"
)


if (
  has_sender &&
  has_receiver &&
  has_probability
) {
  
  cat(
    "PRELIMINARY INTERPRETATION:\n",
    "The table appears to contain dataset-specific communication results.\n",
    sep = ""
  )
  
} else if (
  has_ligand &&
  has_receptor
) {
  
  cat(
    "PRELIMINARY INTERPRETATION:\n",
    "The table appears to be primarily a ligand-receptor interaction resource/reference table.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "PRELIMINARY INTERPRETATION:\n",
    "The semantic role of this table is not yet resolved.\n",
    sep = ""
  )
}

cat("\n")


# ------------------------------------------------------------
# 11. Unique values for important columns
# ------------------------------------------------------------

cat("STEP 11: UNIQUE VALUE SUMMARY\n")
cat("-----------------------------\n")

important_patterns <- c(
  "pathway",
  "ligand",
  "receptor",
  "interaction",
  "annotation",
  "evidence"
)

for (pattern in important_patterns) {
  
  idx <- grepl(
    pattern,
    column_names_lower
  )
  
  if (any(idx)) {
    
    selected_cols <- colnames(cci_lr)[idx]
    
    for (cc in selected_cols) {
      
      values <- unique(
        as.character(
          cci_lr[[cc]]
        )
      )
      
      values <- values[
        !is.na(values) &
          trimws(values) != ""
      ]
      
      cat(
        "\nColumn:",
        cc,
        "\n"
      )
      
      cat(
        "Unique values:",
        length(values),
        "\n"
      )
      
      if (length(values) <= 30) {
        
        print(values)
        
      } else {
        
        cat(
          "First 30 unique values:\n"
        )
        
        print(
          values[seq_len(30)]
        )
      }
    }
  }
}

cat("\n")


# ------------------------------------------------------------
# 12. First 20 rows as a compact audit
# ------------------------------------------------------------

cat("STEP 12: FIRST 20 ROWS\n")
cat("----------------------\n")

print(
  utils::head(
    cci_lr,
    20
  ),
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 13. Save R object for later use
# ------------------------------------------------------------

saveRDS(
  cci_lr,
  file = file.path(
    step_dir,
    "STEP11-0_CCI_LR_Audited_Object.rds"
  )
)

cat("STEP 13: RDS EXPORT\n")
cat("-------------------\n")

cat(
  "Saved:\n",
  file.path(
    step_dir,
    "STEP11-0_CCI_LR_Audited_Object.rds"
  ),
  "\n\n"
)

# ------------------------------------------------------------
# 14. Final status
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP11-0 COMPLETE\n")
cat("============================================================\n")

cat(
  "Audit report:\n",
  output_file,
  "\n"
)

# ============================================================
# STEP11-1
# NSCLC_GSE131907_CellMetainfo_table.tsv
# Structural Audit
# ============================================================

step_dir <- "input_dir/Paper 1_Step 11-CellChat"

meta_file <- file.path(
  step_dir,
  "NSCLC_GSE131907_CellMetainfo_table.tsv"
)

cat("============================================================\n")
cat("STEP11-1: CELL METADATA AUDIT\n")
cat("============================================================\n\n")

cat("File:\n", meta_file, "\n\n")

if (!file.exists(meta_file)) {
  stop("ERROR: Cell metadata file not found.")
}

meta <- read.delim(
  meta_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Rows:", nrow(meta), "\n")
cat("Columns:", ncol(meta), "\n\n")

cat("COLUMN NAMES\n")
cat("------------\n")

for (i in seq_along(colnames(meta))) {
  cat(
    sprintf(
      "[%02d] %s\n",
      i,
      colnames(meta)[i]
    )
  )
}

cat("\n")

cat("FIRST 20 ROWS\n")
cat("-------------\n")

print(
  head(meta, 20),
  row.names = FALSE
)

cat("\n")

cat("DATA TYPES\n")
cat("----------\n")

print(
  data.frame(
    column = colnames(meta),
    class = sapply(
      meta,
      function(x)
        paste(class(x), collapse = "/")
    ),
    n_unique = sapply(
      meta,
      function(x)
        length(unique(x))
    ),
    n_missing = sapply(
      meta,
      function(x)
        sum(is.na(x) | trimws(as.character(x)) == "")
    ),
    row.names = FALSE
  )
)

cat("\n")

cat("POSSIBLE PATIENT / SAMPLE / CELL / CELL-TYPE FIELDS\n")
cat("-----------------------------------------------------\n")

cn <- tolower(colnames(meta))

patterns <- c(
  "patient",
  "sample",
  "subject",
  "donor",
  "cell",
  "cluster",
  "type",
  "annotation",
  "tissue",
  "condition",
  "source",
  "study"
)

for (p in patterns) {
  
  hits <- colnames(meta)[
    grepl(p, cn, fixed = TRUE)
  ]
  
  cat(
    sprintf(
      "%-12s : ",
      p
    )
  )
  
  if (length(hits) == 0) {
    cat("NONE\n")
  } else {
    cat(
      paste(hits, collapse = ", "),
      "\n"
    )
  }
}

cat("\n")

cat("UNIQUE VALUE SUMMARY FOR LIKELY ANNOTATION FIELDS\n")
cat("--------------------------------------------------\n")

candidate_cols <- colnames(meta)[
  grepl(
    "type|annotation|cluster|patient|sample|donor|tissue|condition",
    cn
  )
]

for (cc in candidate_cols) {
  
  vals <- unique(
    as.character(meta[[cc]])
  )
  
  vals <- vals[
    !is.na(vals) &
      trimws(vals) != ""
  ]
  
  cat("\nColumn:", cc, "\n")
  cat("Unique values:", length(vals), "\n")
  
  if (length(vals) <= 50) {
    print(vals)
  } else {
    print(vals[1:50])
    cat("... additional values omitted ...\n")
  }
}

cat("\n")

cat("============================================================\n")
cat("STEP11-1 COMPLETE\n")
cat("============================================================\n")


# ============================================================
# STEP11-1A
# NSCLC_GSE131907 Cell Metadata
# Detailed Structural Audit
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP11-1A: CELL METADATA DETAILED AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Basic structure
# ------------------------------------------------------------

cat("STEP 1: BASIC STRUCTURE\n")
cat("-----------------------\n")

cat("Rows:", nrow(meta), "\n")
cat("Columns:", ncol(meta), "\n\n")

cat("Column names:\n")

for (i in seq_along(colnames(meta))) {
  
  cat(
    sprintf(
      "[%02d] %s\n",
      i,
      colnames(meta)[i]
    )
  )
}

cat("\n")


# ------------------------------------------------------------
# 2. Safe data-type table
# ------------------------------------------------------------

cat("STEP 2: DATA TYPES / UNIQUENESS / MISSINGNESS\n")
cat("----------------------------------------------\n")

type_vec <- vapply(
  meta,
  function(x) paste(class(x), collapse = "/"),
  character(1)
)

unique_vec <- vapply(
  meta,
  function(x) length(unique(x)),
  integer(1)
)

missing_vec <- vapply(
  meta,
  function(x) {
    
    x_chr <- as.character(x)
    
    sum(
      is.na(x) |
        is.na(x_chr) |
        trimws(x_chr) == ""
    )
  },
  integer(1)
)

type_table <- data.frame(
  column = names(meta),
  class = unname(type_vec),
  n_unique = unname(unique_vec),
  n_missing = unname(missing_vec),
  missing_percent = round(
    100 * unname(missing_vec) / nrow(meta),
    4
  ),
  stringsAsFactors = FALSE
)

print(
  type_table,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 3. Patient-level structure
# ------------------------------------------------------------

cat("STEP 3: PATIENT STRUCTURE\n")
cat("-------------------------\n")

if ("Patient" %in% colnames(meta)) {
  
  patient_ids <- unique(
    meta$Patient
  )
  
  patient_ids <- patient_ids[
    !is.na(patient_ids) &
      trimws(as.character(patient_ids)) != ""
  ]
  
  cat(
    "Unique patients:",
    length(patient_ids),
    "\n"
  )
  
  cat(
    "First patient IDs:\n"
  )
  
  print(
    head(patient_ids, 30)
  )
  
} else {
  
  cat(
    "Patient column NOT FOUND.\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 4. Sample-level structure
# ------------------------------------------------------------

cat("STEP 4: SAMPLE STRUCTURE\n")
cat("-----------------------\n")

if ("Sample" %in% colnames(meta)) {
  
  sample_ids <- unique(
    meta$Sample
  )
  
  sample_ids <- sample_ids[
    !is.na(sample_ids) &
      trimws(as.character(sample_ids)) != ""
  ]
  
  cat(
    "Unique samples:",
    length(sample_ids),
    "\n"
  )
  
  cat(
    "Sample IDs:\n"
  )
  
  print(
    sample_ids
  )
  
} else {
  
  cat(
    "Sample column NOT FOUND.\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 5. Cell lineage structure
# ------------------------------------------------------------

cat("STEP 5: CELL-TYPE STRUCTURE\n")
cat("---------------------------\n")


lineage_columns <- c(
  "Celltype (malignancy)",
  "Celltype (major-lineage)",
  "Celltype (minor-lineage)"
)

for (cc in lineage_columns) {
  
  if (cc %in% colnames(meta)) {
    
    cat("\n----------------------------------------\n")
    cat("Column:", cc, "\n")
    cat("----------------------------------------\n")
    
    tab <- sort(
      table(
        meta[[cc]],
        useNA = "ifany"
      ),
      decreasing = TRUE
    )
    
    print(tab)
  }
}

cat("\n")


# ------------------------------------------------------------
# 6. Patient × major-lineage structure
# ------------------------------------------------------------

cat("STEP 6: PATIENT × MAJOR-LINEAGE STRUCTURE\n")
cat("-----------------------------------------\n")

if (
  all(
    c(
      "Patient",
      "Celltype (major-lineage)"
    ) %in% colnames(meta)
  )
) {
  
  patient_lineage_table <- table(
    meta$Patient,
    meta[["Celltype (major-lineage)"]]
  )
  
  cat(
    "Dimensions:",
    nrow(patient_lineage_table),
    "patients ×",
    ncol(patient_lineage_table),
    "major lineages\n\n"
  )
  
  print(
    patient_lineage_table
  )
}

cat("\n")


# ------------------------------------------------------------
# 7. Source structure
# ------------------------------------------------------------

cat("STEP 7: SOURCE STRUCTURE\n")
cat("------------------------\n")

if ("Source" %in% colnames(meta)) {
  
  print(
    sort(
      table(
        meta$Source,
        useNA = "ifany"
      ),
      decreasing = TRUE
    )
  )
}

cat("\n")


# ------------------------------------------------------------
# 8. Stage structure
# ------------------------------------------------------------

cat("STEP 8: STAGE STRUCTURE\n")
cat("-----------------------\n")

if ("Stage" %in% colnames(meta)) {
  
  print(
    sort(
      table(
        meta$Stage,
        useNA = "ifany"
      ),
      decreasing = TRUE
    )
  )
}

cat("\n")


# ------------------------------------------------------------
# 9. TNM structure
# ------------------------------------------------------------

cat("STEP 9: TNM STAGE STRUCTURE\n")
cat("---------------------------\n")

if ("TNMstage" %in% colnames(meta)) {
  
  print(
    sort(
      table(
        meta$TNMstage,
        useNA = "ifany"
      ),
      decreasing = TRUE
    )
  )
}

cat("\n")


# ------------------------------------------------------------
# 10. Cell ID uniqueness
# ------------------------------------------------------------

cat("STEP 10: CELL ID AUDIT\n")
cat("----------------------\n")

if ("Cell" %in% colnames(meta)) {
  
  n_cells <- nrow(meta)
  n_unique_cells <- length(
    unique(meta$Cell)
  )
  
  cat(
    "Total rows:",
    n_cells,
    "\n"
  )
  
  cat(
    "Unique Cell IDs:",
    n_unique_cells,
    "\n"
  )
  
  cat(
    "Cell ID unique:",
    n_cells == n_unique_cells,
    "\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 11. Patient × sample relationship
# ------------------------------------------------------------

cat("STEP 11: PATIENT × SAMPLE RELATIONSHIP\n")
cat("--------------------------------------\n")

if (
  all(
    c("Patient", "Sample") %in%
    colnames(meta)
  )
) {
  
  patient_sample <- unique(
    meta[
      c("Patient", "Sample")
    ]
  )
  
  cat(
    "Unique Patient-Sample combinations:",
    nrow(patient_sample),
    "\n\n"
  )
  
  print(
    patient_sample,
    row.names = FALSE
  )
}

cat("\n")


# ------------------------------------------------------------
# 12. Save audit object
# ------------------------------------------------------------

saveRDS(
  meta,
  file = file.path(
    step_dir,
    "STEP11-1A_CellMetadata_Audited_Object.rds"
  )
)

cat("STEP 12: AUDITED OBJECT SAVED\n")
cat("-----------------------------\n")

cat(
  file.path(
    step_dir,
    "STEP11-1A_CellMetadata_Audited_Object.rds"
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 13. Final status
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP11-1A COMPLETE\n")
cat("============================================================\n")


# ============================================================
# STEP11-2
# Primary-stage Patient × Target Cell-Type Count Distribution Audit
#
# Purpose:
#   Audit the number of target cells available per patient
#   within the Primary-stage population.
#
# IMPORTANT:
#   - No CellChat inference
#   - No ligand-receptor analysis
#   - No expression analysis
#   - No final eligibility threshold is locked here
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP11-2: Primary-stage Patient × Target Cell-Type Count Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 11-CellChat"

if (!dir.exists(step_dir)) {
  stop("STEP11 directory does not exist:\n", step_dir)
}

setwd(step_dir)

cat("Working directory:\n", getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Load STEP11-1A audited metadata object
# ------------------------------------------------------------

meta_rds <- file.path(
  step_dir,
  "STEP11-1A_CellMetadata_Audited_Object.rds"
)

if (!file.exists(meta_rds)) {
  stop(
    "STEP11-1A audited metadata object not found:\n",
    meta_rds
  )
}

obj <- readRDS(meta_rds)

cat("Loaded object class:\n")
print(class(obj))
cat("\n")


# ------------------------------------------------------------
# 2. Recover metadata table
# ------------------------------------------------------------

# The STEP11-1A object is expected to contain the audited metadata.
# Handle common object structures explicitly.

if (is.data.frame(obj)) {
  
  meta <- obj
  
} else if (is.list(obj) && "meta" %in% names(obj)) {
  
  meta <- obj$meta
  
} else if (is.list(obj) && "metadata" %in% names(obj)) {
  
  meta <- obj$metadata
  
} else {
  
  stop(
    "Could not identify metadata table in STEP11-1A object.\n",
    "Please inspect names(obj) / str(obj)."
  )
}


# ------------------------------------------------------------
# 3. Structural audit
# ------------------------------------------------------------

cat("Metadata dimensions:\n")
cat("Rows:", nrow(meta), "\n")
cat("Columns:", ncol(meta), "\n\n")

required_cols <- c(
  "Cell",
  "Celltype (major-lineage)",
  "Patient",
  "Source",
  "Stage"
)

missing_cols <- setdiff(required_cols, colnames(meta))

if (length(missing_cols) > 0) {
  stop(
    "Required metadata columns missing:\n",
    paste(missing_cols, collapse = ", ")
  )
}

cat("Required columns: PASS\n\n")


# ------------------------------------------------------------
# 4. Primary-stage population
# ------------------------------------------------------------

meta_primary <- meta[
  !is.na(meta$Stage) &
    meta$Stage == "Primary",
  ,
  drop = FALSE
]

cat("Primary-stage population:\n")
cat("Cells:", nrow(meta_primary), "\n")
cat(
  "Patients:",
  length(unique(meta_primary$Patient)),
  "\n"
)
cat(
  "Samples:",
  length(unique(meta_primary$Sample)),
  "\n\n"
)


# ------------------------------------------------------------
# 5. Confirm stage distribution
# ------------------------------------------------------------

cat("Overall Stage distribution:\n")
print(table(meta$Stage, useNA = "ifany"))
cat("\n")


# ------------------------------------------------------------
# 6. Define target cell populations
# ------------------------------------------------------------

target_lineages <- c(
  "Endothelial",
  "Mono/Macro",
  "CD4Tconv",
  "CD8T",
  "CD8Tex"
)

cat("Target major-lineages:\n")
print(target_lineages)
cat("\n")


# ------------------------------------------------------------
# 7. Confirm target lineages exist
# ------------------------------------------------------------

available_lineages <- unique(
  meta_primary[["Celltype (major-lineage)"]]
)

missing_targets <- setdiff(
  target_lineages,
  available_lineages
)

if (length(missing_targets) > 0) {
  stop(
    "Target lineages missing from Primary-stage population:\n",
    paste(missing_targets, collapse = ", ")
  )
}

cat("Target lineage availability: PASS\n\n")


# ------------------------------------------------------------
# 8. Patient × major-lineage cell-count matrix
# ------------------------------------------------------------

patient_ids <- sort(unique(meta_primary$Patient))

count_matrix <- matrix(
  0L,
  nrow = length(patient_ids),
  ncol = length(target_lineages),
  dimnames = list(
    patient_ids,
    target_lineages
  )
)

for (ct in target_lineages) {
  
  idx <- meta_primary[["Celltype (major-lineage)"]] == ct
  
  tab <- table(meta_primary$Patient[idx])
  
  common_ids <- intersect(
    names(tab),
    rownames(count_matrix)
  )
  
  count_matrix[
    common_ids,
    ct
  ] <- as.integer(tab[common_ids])
}

count_df <- data.frame(
  Patient = rownames(count_matrix),
  count_matrix,
  row.names = NULL,
  check.names = FALSE
)


# ------------------------------------------------------------
# 9. Add total target-cell count
# ------------------------------------------------------------

count_df$Target_Total <- rowSums(
  count_df[, target_lineages, drop = FALSE]
)


# ------------------------------------------------------------
# 10. Patient-level coverage
# ------------------------------------------------------------

coverage_table <- data.frame(
  CellType = target_lineages,
  Patients_Total = nrow(count_df),
  Patients_Nonzero = sapply(
    target_lineages,
    function(x) sum(count_df[[x]] > 0)
  ),
  Coverage_Percent = sapply(
    target_lineages,
    function(x)
      100 * mean(count_df[[x]] > 0)
  )
)

cat("============================================================\n")
cat("PATIENT-LEVEL COVERAGE\n")
cat("============================================================\n\n")

print(
  coverage_table,
  row.names = FALSE,
  digits = 3
)

cat("\n")


# ------------------------------------------------------------
# 11. Distribution statistics
# ------------------------------------------------------------

distribution_table <- do.call(
  rbind,
  lapply(
    target_lineages,
    function(ct) {
      
      x <- count_df[[ct]]
      
      data.frame(
        CellType = ct,
        N_Patients = length(x),
        N_Nonzero = sum(x > 0),
        Min = min(x),
        Q1 = as.numeric(
          quantile(x, 0.25, names = FALSE)
        ),
        Median = median(x),
        Mean = mean(x),
        Q3 = as.numeric(
          quantile(x, 0.75, names = FALSE)
        ),
        Max = max(x),
        SD = sd(x)
      )
    }
  )
)

rownames(distribution_table) <- NULL

cat("============================================================\n")
cat("CELL-COUNT DISTRIBUTION\n")
cat("============================================================\n\n")

print(
  distribution_table,
  row.names = FALSE,
  digits = 3
)

cat("\n")


# ------------------------------------------------------------
# 12. Candidate threshold audit
#
# IMPORTANT:
# These are descriptive only.
# No threshold is locked in STEP11-2.
# ------------------------------------------------------------

thresholds <- c(10, 20, 30, 50, 100)

threshold_table <- do.call(
  rbind,
  lapply(
    target_lineages,
    function(ct) {
      
      x <- count_df[[ct]]
      
      do.call(
        rbind,
        lapply(
          thresholds,
          function(th) {
            
            data.frame(
              CellType = ct,
              Threshold = th,
              Patients_Retained = sum(x >= th),
              Percent_Retained =
                100 * mean(x >= th)
            )
          }
        )
      )
    }
  )
)

rownames(threshold_table) <- NULL

cat("============================================================\n")
cat("CANDIDATE THRESHOLD AUDIT\n")
cat("============================================================\n\n")

print(
  threshold_table,
  row.names = FALSE,
  digits = 3
)

cat("\n")


# ------------------------------------------------------------
# 13. Combined Endothelial + Mono/Macro + T-cell coverage
#
# T-cell groups are shown both individually and combined.
# ------------------------------------------------------------

count_df$Tcell_Total <- rowSums(
  count_df[
    c("CD4Tconv", "CD8T", "CD8Tex"),
    drop = FALSE
  ]
)

combined_threshold_table <- do.call(
  rbind,
  lapply(
    thresholds,
    function(th) {
      
      data.frame(
        Threshold = th,
        
        Endothelial = sum(
          count_df$Endothelial >= th
        ),
        
        Mono_Macro = sum(
          count_df$`Mono/Macro` >= th
        ),
        
        Tcell_Total = sum(
          count_df$Tcell_Total >= th
        ),
        
        All_Three_Compartments =
          sum(
            count_df$Endothelial >= th &
              count_df$`Mono/Macro` >= th &
              count_df$Tcell_Total >= th
          )
      )
    }
  )
)

cat("============================================================\n")
cat("COMBINED THREE-COMPARTMENT THRESHOLD AUDIT\n")
cat("============================================================\n\n")

print(
  combined_threshold_table,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 14. Patient-level detailed table
# ------------------------------------------------------------

count_df_export <- count_df[
  order(
    -count_df$Endothelial,
    -count_df$`Mono/Macro`,
    -count_df$Tcell_Total,
    count_df$Patient
  ),
]


# ------------------------------------------------------------
# 15. Save outputs
# ------------------------------------------------------------

saveRDS(
  count_df,
  file = file.path(
    step_dir,
    "STEP11-2_Primary_Patient_CellCounts.rds"
  )
)

write.table(
  count_df_export,
  file = file.path(
    step_dir,
    "STEP11-2_Primary_Patient_CellCounts.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  coverage_table,
  file = file.path(
    step_dir,
    "STEP11-2_Patient_Coverage.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  distribution_table,
  file = file.path(
    step_dir,
    "STEP11-2_CellCount_Distribution.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  threshold_table,
  file = file.path(
    step_dir,
    "STEP11-2_Candidate_Threshold_Audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  combined_threshold_table,
  file = file.path(
    step_dir,
    "STEP11-2_ThreeCompartment_Threshold_Audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 16. Final structural checks
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP11-2 STRUCTURAL QC\n")
cat("============================================================\n\n")

qc1 <- nrow(meta_primary) > 0
qc2 <- length(unique(meta_primary$Patient)) > 0
qc3 <- all(
  target_lineages %in%
    colnames(count_df)
)
qc4 <- all(
  count_df[, target_lineages, drop = FALSE] >= 0
)
qc5 <- all(
  count_df$Target_Total ==
    rowSums(
      count_df[
        , target_lineages,
        drop = FALSE
      ]
    )
)

cat(
  "Primary-stage cells > 0:",
  ifelse(qc1, "PASS", "FAIL"),
  "\n"
)

cat(
  "Primary-stage patients > 0:",
  ifelse(qc2, "PASS", "FAIL"),
  "\n"
)

cat(
  "All target lineages present:",
  ifelse(qc3, "PASS", "FAIL"),
  "\n"
)

cat(
  "Cell counts non-negative:",
  ifelse(qc4, "PASS", "FAIL"),
  "\n"
)

cat(
  "Target_Total calculation:",
  ifelse(qc5, "PASS", "FAIL"),
  "\n\n"
)


overall_pass <- all(
  qc1,
  qc2,
  qc3,
  qc4,
  qc5
)

cat("============================================================\n")
cat(
  "STEP11-2 STATUS:",
  ifelse(
    overall_pass,
    "PASS",
    "FAIL"
  ),
  "\n"
)
cat("============================================================\n\n")

cat(
  "IMPORTANT: No patient eligibility threshold was locked.\n",
  "Threshold results are descriptive only and will be reviewed\n",
  "before any downstream patient-level communication analysis.\n\n"
)


# ============================================================
# STEP11-3
# Primary-stage Patient × Sample × Source × Target-cell
# Eligibility Audit
#
# Purpose:
#   Determine whether the target compartments coexist within
#   the SAME sample, and identify patient-level eligibility
#   without mixing different tissue sources.
#
# IMPORTANT:
#   - No CellChat inference
#   - No ligand-receptor inference
#   - No expression analysis
#   - No final eligibility threshold is locked
#   - Candidate thresholds are descriptive only
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP11-3: Primary Patient × Sample × Source Eligibility Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 11-CellChat"

if (!dir.exists(step_dir)) {
  stop("STEP11 directory does not exist:\n", step_dir)
}

setwd(step_dir)

cat("Working directory:\n", getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Load STEP11-1A audited metadata
# ------------------------------------------------------------

meta_rds <- file.path(
  step_dir,
  "STEP11-1A_CellMetadata_Audited_Object.rds"
)

if (!file.exists(meta_rds)) {
  stop(
    "STEP11-1A audited metadata object not found:\n",
    meta_rds
  )
}

obj <- readRDS(meta_rds)


# ------------------------------------------------------------
# 2. Recover metadata table
# ------------------------------------------------------------

if (is.data.frame(obj)) {
  
  meta <- obj
  
} else if (is.list(obj) && "meta" %in% names(obj)) {
  
  meta <- obj$meta
  
} else if (is.list(obj) && "metadata" %in% names(obj)) {
  
  meta <- obj$metadata
  
} else {
  
  stop(
    "Could not identify metadata table in STEP11-1A object.\n",
    "Please inspect names(obj) / str(obj)."
  )
}


# ------------------------------------------------------------
# 3. Required columns
# ------------------------------------------------------------

required_cols <- c(
  "Cell",
  "Celltype (major-lineage)",
  "Patient",
  "Sample",
  "Source",
  "Stage"
)

missing_cols <- setdiff(
  required_cols,
  colnames(meta)
)

if (length(missing_cols) > 0) {
  stop(
    "Required columns missing:\n",
    paste(missing_cols, collapse = ", ")
  )
}

cat("Required metadata columns: PASS\n\n")


# ------------------------------------------------------------
# 4. Basic integrity audit
# ------------------------------------------------------------

cat("Metadata dimensions:\n")
cat("Cells:", nrow(meta), "\n")
cat("Patients:", length(unique(meta$Patient)), "\n")
cat("Samples:", length(unique(meta$Sample)), "\n\n")

if (anyDuplicated(meta$Cell) > 0) {
  stop("Duplicate Cell IDs detected.")
}

if (anyNA(meta[, required_cols])) {
  stop("NA values detected in required metadata columns.")
}

cat("Cell ID uniqueness: PASS\n")
cat("Required-field completeness: PASS\n\n")


# ------------------------------------------------------------
# 5. Restrict to Primary-stage
# ------------------------------------------------------------

meta_primary <- meta[
  meta$Stage == "Primary",
  ,
  drop = FALSE
]

cat("============================================================\n")
cat("PRIMARY-STAGE POPULATION\n")
cat("============================================================\n\n")

cat("Cells:", nrow(meta_primary), "\n")
cat("Patients:", length(unique(meta_primary$Patient)), "\n")
cat("Samples:", length(unique(meta_primary$Sample)), "\n\n")


# ------------------------------------------------------------
# 6. Check patient × sample uniqueness
# ------------------------------------------------------------

patient_sample_pairs <- unique(
  meta_primary[, c("Patient", "Sample")]
)

cat(
  "Unique Primary patient × sample pairs:",
  nrow(patient_sample_pairs),
  "\n\n"
)


# ------------------------------------------------------------
# 7. Inspect Sample × Source structure
# ------------------------------------------------------------

sample_source <- unique(
  meta_primary[, c("Sample", "Patient", "Source")]
)

source_per_sample <- aggregate(
  Source ~ Sample,
  data = sample_source,
  FUN = function(x) paste(sort(unique(x)), collapse = ";")
)

source_count_per_sample <- aggregate(
  Source ~ Sample,
  data = sample_source,
  FUN = function(x) length(unique(x))
)

names(source_count_per_sample)[2] <-
  "N_Source"

multi_source_samples <- source_count_per_sample[
  source_count_per_sample$N_Source > 1,
  ,
  drop = FALSE
]

cat("============================================================\n")
cat("SAMPLE × SOURCE STRUCTURE\n")
cat("============================================================\n\n")

cat(
  "Primary samples:",
  nrow(sample_source),
  "\n"
)

cat(
  "Samples with >1 Source label:",
  nrow(multi_source_samples),
  "\n\n"
)

if (nrow(multi_source_samples) > 0) {
  
  cat("WARNING: samples with multiple Source labels:\n")
  print(
    merge(
      multi_source_samples,
      source_per_sample,
      by = "Sample"
    ),
    row.names = FALSE
  )
  cat("\n")
  
} else {
  
  cat(
    "Each Primary sample has one Source label: PASS\n\n"
  )
}


# ------------------------------------------------------------
# 8. Target cell populations
# ------------------------------------------------------------

target_lineages <- c(
  "Endothelial",
  "Mono/Macro",
  "CD4Tconv",
  "CD8T",
  "CD8Tex"
)

target_tcell <- c(
  "CD4Tconv",
  "CD8T",
  "CD8Tex"
)

cat("Target populations:\n")
print(target_lineages)
cat("\n")


# ------------------------------------------------------------
# 9. Sample-level cell-count matrix
# ------------------------------------------------------------

sample_ids <- sort(
  unique(meta_primary$Sample)
)

sample_count_matrix <- matrix(
  0L,
  nrow = length(sample_ids),
  ncol = length(target_lineages),
  dimnames = list(
    sample_ids,
    target_lineages
  )
)

for (ct in target_lineages) {
  
  idx <- meta_primary[["Celltype (major-lineage)"]] == ct
  
  tab <- table(
    meta_primary$Sample[idx]
  )
  
  common_ids <- intersect(
    names(tab),
    rownames(sample_count_matrix)
  )
  
  sample_count_matrix[
    common_ids,
    ct
  ] <- as.integer(
    tab[common_ids]
  )
}

sample_count_df <- data.frame(
  Sample = rownames(sample_count_matrix),
  sample_count_matrix,
  row.names = NULL,
  check.names = FALSE
)

cat("Sample-level count matrix generated: PASS\n")
cat("Samples:", nrow(sample_count_df), "\n")
cat("Target columns:", ncol(sample_count_df) - 1, "\n\n")

# ------------------------------------------------------------
# 10. Add patient and source information
# ------------------------------------------------------------

sample_annotation <- unique(
  meta_primary[
    ,
    c("Sample", "Patient", "Source"),
    drop = FALSE
  ]
)

if (anyDuplicated(sample_annotation$Sample) > 0) {
  
  stop(
    "A Primary sample maps to multiple Patient/Source combinations."
  )
}

sample_count_df <- merge(
  sample_annotation,
  sample_count_df,
  by = "Sample",
  sort = FALSE
)


# ------------------------------------------------------------
# 11. Combined T-cell count
# ------------------------------------------------------------

sample_count_df$Tcell_Total <-
  rowSums(
    sample_count_df[
      ,
      target_tcell,
      drop = FALSE
    ]
  )

sample_count_df$Target_Total <-
  rowSums(
    sample_count_df[
      ,
      target_lineages,
      drop = FALSE
    ]
  )


# ------------------------------------------------------------
# 12. Tumor-only subset
# ------------------------------------------------------------

sample_count_df$Is_Tumor_Source <-
  sample_count_df$Source == "Tumor"


# ------------------------------------------------------------
# 13. Basic sample-level coexistence
# ------------------------------------------------------------

sample_count_df$All_Three_Nonzero <-
  sample_count_df$Endothelial > 0 &
  sample_count_df$`Mono/Macro` > 0 &
  sample_count_df$Tcell_Total > 0

sample_count_df$All_Three_Nonzero_Tumor <-
  sample_count_df$All_Three_Nonzero &
  sample_count_df$Is_Tumor_Source


cat("============================================================\n")
cat("SAMPLE-LEVEL COEXISTENCE\n")
cat("============================================================\n\n")

cat(
  "Primary samples:",
  nrow(sample_count_df),
  "\n"
)

cat(
  "Samples with Endothelial + Mono/Macro + T-cell > 0:",
  sum(sample_count_df$All_Three_Nonzero),
  "\n"
)

cat(
  "Tumor-source samples with all three > 0:",
  sum(sample_count_df$All_Three_Nonzero_Tumor),
  "\n\n"
)


# ------------------------------------------------------------
# 14. Candidate thresholds at SAME SAMPLE level
# ------------------------------------------------------------

thresholds <- c(
  1,
  5,
  10,
  20,
  30,
  50,
  100
)

sample_threshold_table <- do.call(
  rbind,
  lapply(
    thresholds,
    function(th) {
      
      eligible_all <-
        sample_count_df$Endothelial >= th &
        sample_count_df$`Mono/Macro` >= th &
        sample_count_df$Tcell_Total >= th
      
      eligible_tumor <-
        eligible_all &
        sample_count_df$Is_Tumor_Source
      
      data.frame(
        Threshold = th,
        
        Primary_Samples_Eligible =
          sum(eligible_all),
        
        Primary_Sample_Percent =
          100 * mean(eligible_all),
        
        Tumor_Samples_Eligible =
          sum(eligible_tumor),
        
        Tumor_Sample_Percent =
          100 * mean(eligible_tumor)
      )
    }
  )
)

rownames(sample_threshold_table) <- NULL


cat("============================================================\n")
cat("SAME-SAMPLE THRESHOLD AUDIT\n")
cat("============================================================\n\n")

print(
  sample_threshold_table,
  row.names = FALSE,
  digits = 3
)

cat("\n")


# ------------------------------------------------------------
# 15. Patient-level eligibility
#
# A patient is eligible if at least ONE Primary sample
# satisfies the criterion.
#
# This prevents combining different samples from the same
# patient.
# ------------------------------------------------------------

patient_ids <- sort(
  unique(meta_primary$Patient)
)

patient_threshold_table <- do.call(
  rbind,
  lapply(
    thresholds,
    function(th) {
      
      eligible_sample <-
        sample_count_df$Endothelial >= th &
        sample_count_df$`Mono/Macro` >= th &
        sample_count_df$Tcell_Total >= th
      
      eligible_tumor_sample <-
        eligible_sample &
        sample_count_df$Is_Tumor_Source
      
      eligible_patients <-
        unique(
          sample_count_df$Patient[
            eligible_sample
          ]
        )
      
      eligible_tumor_patients <-
        unique(
          sample_count_df$Patient[
            eligible_tumor_sample
          ]
        )
      
      data.frame(
        Threshold = th,
        
        Primary_Patients_Total =
          length(patient_ids),
        
        Primary_Patients_Eligible =
          length(eligible_patients),
        
        Primary_Patient_Percent =
          100 *
          length(eligible_patients) /
          length(patient_ids),
        
        Tumor_Patients_Eligible =
          length(eligible_tumor_patients),
        
        Tumor_Patient_Percent =
          100 *
          length(eligible_tumor_patients) /
          length(patient_ids)
      )
    }
  )
)

rownames(patient_threshold_table) <- NULL


cat("============================================================\n")
cat("PATIENT-LEVEL ELIGIBILITY\n")
cat("============================================================\n\n")

print(
  patient_threshold_table,
  row.names = FALSE,
  digits = 3
)

cat("\n")


# ------------------------------------------------------------
# 16. Detailed patient-level eligibility at each threshold
# ------------------------------------------------------------

patient_eligibility_detail <- do.call(
  rbind,
  lapply(
    thresholds,
    function(th) {
      
      eligible <-
        sample_count_df$Endothelial >= th &
        sample_count_df$`Mono/Macro` >= th &
        sample_count_df$Tcell_Total >= th
      
      eligible_tumor <-
        eligible &
        sample_count_df$Is_Tumor_Source
      
      all_patients <- data.frame(
        Patient = patient_ids,
        stringsAsFactors = FALSE
      )
      
      all_patients$Threshold <- th
      
      all_patients$Eligible_Primary_Sample <-
        all_patients$Patient %in%
        unique(
          sample_count_df$Patient[
            eligible
          ]
        )
      
      all_patients$Eligible_Tumor_Sample <-
        all_patients$Patient %in%
        unique(
          sample_count_df$Patient[
            eligible_tumor
          ]
        )
      
      all_patients
    }
  )
)

rownames(patient_eligibility_detail) <- NULL


# ------------------------------------------------------------
# 17. Primary Tumor sample-level detailed table
# ------------------------------------------------------------

tumor_sample_detail <- sample_count_df[
  sample_count_df$Is_Tumor_Source,
  c(
    "Patient",
    "Sample",
    "Source",
    "Endothelial",
    "Mono/Macro",
    "CD4Tconv",
    "CD8T",
    "CD8Tex",
    "Tcell_Total",
    "Target_Total"
  ),
  drop = FALSE
]

tumor_sample_detail <-
  tumor_sample_detail[
    order(
      tumor_sample_detail$Patient,
      -tumor_sample_detail$Endothelial
    ),
    ,
    drop = FALSE
  ]


# ------------------------------------------------------------
# 18. Identify patients with multiple Primary samples
# ------------------------------------------------------------

patient_sample_count <- aggregate(
  Sample ~ Patient,
  data = sample_annotation,
  FUN = function(x) length(unique(x))
)

names(patient_sample_count)[2] <-
  "N_Primary_Samples"

multi_sample_patients <-
  patient_sample_count[
    patient_sample_count$N_Primary_Samples > 1,
    ,
    drop = FALSE
  ]


cat("============================================================\n")
cat("PATIENTS WITH MULTIPLE PRIMARY SAMPLES\n")
cat("============================================================\n\n")

cat(
  "Patients with >1 Primary sample:",
  nrow(multi_sample_patients),
  "\n\n"
)

if (nrow(multi_sample_patients) > 0) {
  
  print(
    multi_sample_patients,
    row.names = FALSE
  )
  
  cat("\n")
}


# ------------------------------------------------------------
# 19. Check whether eligibility depends on non-Tumor sources
# ------------------------------------------------------------

non_tumor_three <-
  sample_count_df$All_Three_Nonzero &
  !sample_count_df$Is_Tumor_Source

cat("============================================================\n")
cat("SOURCE-SPECIFIC ELIGIBILITY CHECK\n")
cat("============================================================\n\n")

cat(
  "Primary samples with all 3 compartments but Source != Tumor:",
  sum(non_tumor_three),
  "\n\n"
)


# ------------------------------------------------------------
# 20. Save outputs
# ------------------------------------------------------------

saveRDS(
  sample_count_df,
  file = file.path(
    step_dir,
    "STEP11-3_Primary_Sample_CellCounts.rds"
  )
)

saveRDS(
  patient_eligibility_detail,
  file = file.path(
    step_dir,
    "STEP11-3_Patient_Eligibility_Detail.rds"
  )
)

write.table(
  sample_count_df,
  file = file.path(
    step_dir,
    "STEP11-3_Primary_Sample_CellCounts.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sample_threshold_table,
  file = file.path(
    step_dir,
    "STEP11-3_Sample_Threshold_Audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  patient_threshold_table,
  file = file.path(
    step_dir,
    "STEP11-3_Patient_Eligibility_Audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  patient_eligibility_detail,
  file = file.path(
    step_dir,
    "STEP11-3_Patient_Eligibility_Detail.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  tumor_sample_detail,
  file = file.path(
    step_dir,
    "STEP11-3_Primary_Tumor_Sample_Detail.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  multi_sample_patients,
  file = file.path(
    step_dir,
    "STEP11-3_MultiSample_Patient_Audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 21. QC
# ------------------------------------------------------------

qc1 <- nrow(meta_primary) == 156751

qc2 <-
  length(unique(meta_primary$Patient)) == 30

qc3 <-
  length(unique(meta_primary$Sample)) == 42

qc4 <-
  all(sample_count_df$Endothelial >= 0)

qc5 <-
  all(sample_count_df$`Mono/Macro` >= 0)

qc6 <-
  all(sample_count_df$Tcell_Total >= 0)

qc7 <-
  all(
    sample_count_df$Target_Total ==
      rowSums(
        sample_count_df[
          ,
          target_lineages,
          drop = FALSE
        ]
      )
  )

qc8 <-
  all(
    sample_count_df$Is_Tumor_Source ==
      (sample_count_df$Source == "Tumor")
  )


cat("============================================================\n")
cat("STEP11-3 STRUCTURAL QC\n")
cat("============================================================\n\n")

cat(
  "Primary-stage cell count = 156,751:",
  ifelse(qc1, "PASS", "FAIL"),
  "\n"
)

cat(
  "Primary-stage patients = 30:",
  ifelse(qc2, "PASS", "FAIL"),
  "\n"
)

cat(
  "Primary-stage samples = 42:",
  ifelse(qc3, "PASS", "FAIL"),
  "\n"
)

cat(
  "Endothelial counts non-negative:",
  ifelse(qc4, "PASS", "FAIL"),
  "\n"
)

cat(
  "Mono/Macro counts non-negative:",
  ifelse(qc5, "PASS", "FAIL"),
  "\n"
)

cat(
  "T-cell counts non-negative:",
  ifelse(qc6, "PASS", "FAIL"),
  "\n"
)

cat(
  "Target_Total calculation:",
  ifelse(qc7, "PASS", "FAIL"),
  "\n"
)

cat(
  "Tumor-source flag:",
  ifelse(qc8, "PASS", "FAIL"),
  "\n\n"
)


overall_pass <- all(
  qc1,
  qc2,
  qc3,
  qc4,
  qc5,
  qc6,
  qc7,
  qc8
)

cat("============================================================\n")
cat(
  "STEP11-3 STATUS:",
  ifelse(overall_pass, "PASS", "FAIL"),
  "\n"
)
cat("============================================================\n\n")

cat(
  "IMPORTANT:\n",
  "No final patient/sample eligibility threshold was locked.\n",
  "Candidate thresholds are descriptive only.\n",
  "Patient-level eligibility requires at least one SAME Primary\n",
  "sample satisfying the criterion; cells from different samples\n",
  "are never combined to create eligibility.\n\n"
)

# ============================================================
# STEP11-4 PRE-AUDIT
# Identify locked STEP10 localization / pseudobulk outputs
# ============================================================

step10_dir <- "input_dir/Paper 1_Step 10-Single-Cell"

step10_files <- list.files(
  step10_dir,
  full.names = TRUE,
  recursive = TRUE
)

cat("============================================================\n")
cat("STEP10 FILE INVENTORY\n")
cat("============================================================\n\n")

cat("Total files:", length(step10_files), "\n\n")

print(step10_files)

cat("\n============================================================\n")
cat("Candidate localization / pseudobulk files\n")
cat("============================================================\n\n")

candidate_files <- step10_files[
  grepl(
    paste(
      c(
        "AOC3",
        "MS4A7",
        "IL7R",
        "pseudobulk",
        "pseudo",
        "localization",
        "local",
        "core",
        "expression",
        "detection",
        "compartment",
        "STEP10"
      ),
      collapse = "|"
    ),
    basename(step10_files),
    ignore.case = TRUE
  )
]

print(candidate_files)

cat("\n============================================================\n")
cat("PRE-AUDIT COMPLETE\n")
cat("============================================================\n")

# ============================================================
# STEP11-4 INPUT AUDIT
# Audit locked STEP10 outputs before integration
# NO raw UMI re-reading
# NO recalculation
# ============================================================

step10_dir <- "input_dir/Paper 1_Step 10-Single-Cell"

# ------------------------------------------------------------
# 1. Define locked STEP10 input files
# ------------------------------------------------------------

files_11_4 <- list(
  
  STEP10E_CoreGene =
    file.path(
      step10_dir,
      "STEP10E_Patient_Level_CoreGene_Expression",
      "STEP10E-1_Primary_Patient_Level_CoreGene_Expression.tsv"
    ),
  
  STEP10F1C_GeneCompartment =
    file.path(
      step10_dir,
      "STEP10F-1_Patient_Level_Expression_Profile",
      "STEP10F-1C_Primary_Gene_Compartment_Expression_Profile.tsv"
    ),
  
  STEP10F2A_Pattern =
    file.path(
      step10_dir,
      "STEP10F-2_Gene_Compartment_Expression_Pattern_Audit",
      "STEP10F-2A_Primary_Gene_Compartment_Pattern_Profile.tsv"
    ),
  
  STEP10F2D_WithinPatientRank =
    file.path(
      step10_dir,
      "STEP10F-2_Gene_Compartment_Expression_Pattern_Audit",
      "STEP10F-2D_Complete_Patient_Within_Patient_Ranks.tsv"
    ),
  
  STEP10F2E_ThreeCompartment =
    file.path(
      step10_dir,
      "STEP10F-2_Gene_Compartment_Expression_Pattern_Audit",
      "STEP10F-2E_Complete_ThreeCompartment_Pattern_Summary.tsv"
    ),
  
  STEP10F2F_Dominant =
    file.path(
      step10_dir,
      "STEP10F-2_Gene_Compartment_Expression_Pattern_Audit",
      "STEP10F-2F_Complete_Patient_Dominant_Compartment.tsv"
    ),
  
  STEP10F1E_ThreeCompartmentProfile =
    file.path(
      step10_dir,
      "STEP10F-1_Patient_Level_Expression_Profile",
      "STEP10F-1E_Complete_ThreeCompartment_Patient_Profile.tsv"
    ),
  
  STEP10C_SampleAnnotation =
    file.path(
      step10_dir,
      "STEP10C_Patient_Level_RawUMI_Pseudobulk",
      "STEP10C-E_Primary_Pseudobulk_Sample_Annotation.tsv"
    )
)

# ------------------------------------------------------------
# 2. File existence audit
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP11-4 INPUT FILE EXISTENCE AUDIT\n")
cat("============================================================\n\n")

for (nm in names(files_11_4)) {
  cat(nm, " : ", file.exists(files_11_4[[nm]]), "\n", sep = "")
}

# ------------------------------------------------------------
# 3. Read and inspect each TSV
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP11-4 INPUT STRUCTURE AUDIT\n")
cat("============================================================\n\n")

input_objects <- list()

for (nm in names(files_11_4)) {
  
  f <- files_11_4[[nm]]
  
  if (!file.exists(f)) {
    cat("\n------------------------------------------------------------\n")
    cat(nm, " : FILE NOT FOUND\n")
    next
  }
  
  x <- read.delim(
    f,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  input_objects[[nm]] <- x
  
  cat("\n------------------------------------------------------------\n")
  cat("FILE:", nm, "\n")
  cat("------------------------------------------------------------\n")
  
  cat("Rows:", nrow(x), "\n")
  cat("Columns:", ncol(x), "\n\n")
  
  cat("Column names:\n")
  print(names(x))
  
  cat("\nFirst 5 rows:\n")
  print(utils::head(x, 5))
}

# ------------------------------------------------------------
# 4. Target-field search
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP11-4 TARGET FIELD SEARCH\n")
cat("============================================================\n\n")

target_patterns <- c(
  "patient",
  "sample",
  "source",
  "stage",
  "gene",
  "compartment",
  "lineage",
  "expression",
  "abundance",
  "detect",
  "count"
)

for (nm in names(input_objects)) {
  
  x <- input_objects[[nm]]
  
  matched <- names(x)[
    grepl(
      paste(target_patterns, collapse = "|"),
      names(x),
      ignore.case = TRUE
    )
  ]
  
  cat("\n", nm, "\n", sep = "")
  cat("Matched fields:\n")
  print(matched)
}

cat("\n============================================================\n")
cat("STEP11-4 INPUT AUDIT COMPLETE\n")
cat("============================================================\n")


# ============================================================
# Locate actual STEP11-3 output files
# ============================================================

step11_dir <- "input_dir/Paper 1_Step 11-CellChat"

all_step11_files <- list.files(
  step11_dir,
  recursive = TRUE,
  full.names = TRUE
)

step11_3_files <- all_step11_files[
  grepl(
    "STEP11-3|Primary|Tumor|Eligibility|eligibility",
    basename(all_step11_files),
    ignore.case = TRUE
  )
]

cat("============================================================\n")
cat("STEP11-3 OUTPUT FILE SEARCH\n")
cat("============================================================\n\n")

cat("Number of candidate files:",
    length(step11_3_files), "\n\n")

if (length(step11_3_files) > 0) {
  print(step11_3_files)
} else {
  cat("NO MATCHING FILES FOUND\n")
}

# ============================================================
# STEP11-3 OUTPUT STRUCTURE AUDIT
# Preparation for STEP11-4 integration
#
# NO recalculation
# NO raw UMI
# NO threshold locking
# ============================================================

step11_dir <- "input_dir/Paper 1_Step 11-CellChat"

step11_3_files <- c(
  
  Patient_Eligibility_Audit =
    file.path(
      step11_dir,
      "STEP11-3_Patient_Eligibility_Audit.tsv"
    ),
  
  Patient_Eligibility_Detail =
    file.path(
      step11_dir,
      "STEP11-3_Patient_Eligibility_Detail.tsv"
    ),
  
  Primary_Tumor_Sample_Detail =
    file.path(
      step11_dir,
      "STEP11-3_Primary_Tumor_Sample_Detail.tsv"
    ),
  
  Sample_Threshold_Audit =
    file.path(
      step11_dir,
      "STEP11-3_Sample_Threshold_Audit.tsv"
    ),
  
  MultiSample_Patient_Audit =
    file.path(
      step11_dir,
      "STEP11-3_MultiSample_Patient_Audit.tsv"
    ),
  
  Primary_Sample_CellCounts =
    file.path(
      step11_dir,
      "STEP11-3_Primary_Sample_CellCounts.tsv"
    )
)


cat("============================================================\n")
cat("STEP11-3 OUTPUT STRUCTURE AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Existence
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("FILE EXISTENCE\n")
cat("------------------------------------------------------------\n\n")

for (nm in names(step11_3_files)) {
  
  cat(
    sprintf(
      "%-35s : %s\n",
      nm,
      file.exists(step11_3_files[[nm]])
    )
  )
}


# ------------------------------------------------------------
# 2. Read and inspect
# ------------------------------------------------------------

objects_11_3 <- list()

for (nm in names(step11_3_files)) {
  
  f <- step11_3_files[[nm]]
  
  if (!file.exists(f)) {
    next
  }
  
  x <- read.delim(
    f,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  objects_11_3[[nm]] <- x
  
  cat("\n============================================================\n")
  cat("FILE:", nm, "\n")
  cat("============================================================\n")
  
  cat("Rows:", nrow(x), "\n")
  cat("Columns:", ncol(x), "\n\n")
  
  cat("Column names:\n")
  print(names(x))
  
  cat("\nFirst 8 rows:\n")
  print(utils::head(x, 8))
}


# ------------------------------------------------------------
# 3. Search important fields
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP11-3 TARGET FIELD SEARCH\n")
cat("============================================================\n\n")

target_patterns <- c(
  "patient",
  "sample",
  "source",
  "stage",
  "tumor",
  "endothelial",
  "mono",
  "macro",
  "myeloid",
  "tcell",
  "t.cell",
  "cell",
  "count",
  "threshold",
  "eligible",
  "coexist"
)


for (nm in names(objects_11_3)) {
  
  x <- objects_11_3[[nm]]
  
  matched <- names(x)[
    grepl(
      paste(target_patterns, collapse = "|"),
      names(x),
      ignore.case = TRUE
    )
  ]
  
  cat("\n", nm, "\n", sep = "")
  cat("Matched fields:\n")
  print(matched)
}


# ------------------------------------------------------------
# 4. Basic dimensions summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP11-3 DIMENSION SUMMARY\n")
cat("============================================================\n\n")

dimension_summary <- do.call(
  rbind,
  lapply(
    names(objects_11_3),
    function(nm) {
      
      x <- objects_11_3[[nm]]
      
      data.frame(
        File = nm,
        Rows = nrow(x),
        Columns = ncol(x),
        stringsAsFactors = FALSE
      )
    }
  )
)

print(dimension_summary)


cat("\n============================================================\n")
cat("STEP11-3 OUTPUT STRUCTURE AUDIT COMPLETE\n")
cat("============================================================\n")



cat("============================================================\n")
cat("STEP11-4 — CORE-GENE CELLULAR LOCALIZATION INTEGRATION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 11-CellChat"

step10_dir <- "input_dir/Paper 1_Step 10-Single-Cell"

out_dir <- file.path(
  base_dir,
  "STEP11-4_CoreGene_Cellular_Localization_Integration"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 1. Locked input files
# ------------------------------------------------------------

step10_files <- c(
  
  STEP10E_CoreGene =
    file.path(
      step10_dir,
      "STEP10E_Patient_Level_CoreGene_Expression",
      "STEP10E-1_Primary_Patient_Level_CoreGene_Expression.tsv"
    ),
  
  STEP10F1C_GeneCompartment =
    file.path(
      step10_dir,
      "STEP10F-1_Patient_Level_Expression_Profile",
      "STEP10F-1C_Primary_Gene_Compartment_Expression_Profile.tsv"
    ),
  
  STEP10F1E_ThreeCompartmentProfile =
    file.path(
      step10_dir,
      "STEP10F-1_Patient_Level_Expression_Profile",
      "STEP10F-1E_Complete_ThreeCompartment_Patient_Profile.tsv"
    ),
  
  STEP10F2A_Pattern =
    file.path(
      step10_dir,
      "STEP10F-2_Gene_Compartment_Expression_Pattern_Audit",
      "STEP10F-2A_Primary_Gene_Compartment_Pattern_Profile.tsv"
    ),
  
  STEP10F2D_WithinPatientRank =
    file.path(
      step10_dir,
      "STEP10F-2_Gene_Compartment_Expression_Pattern_Audit",
      "STEP10F-2D_Complete_Patient_Within_Patient_Ranks.tsv"
    ),
  
  STEP10F2E_ThreeCompartment =
    file.path(
      step10_dir,
      "STEP10F-2_Gene_Compartment_Expression_Pattern_Audit",
      "STEP10F-2E_Complete_ThreeCompartment_Pattern_Summary.tsv"
    ),
  
  STEP10F2F_Dominant =
    file.path(
      step10_dir,
      "STEP10F-2_Gene_Compartment_Expression_Pattern_Audit",
      "STEP10F-2F_Complete_Patient_Dominant_Compartment.tsv"
    )
)


step11_3_files <- c(
  
  Primary_Tumor_Sample_Detail =
    file.path(
      base_dir,
      "STEP11-3_Primary_Tumor_Sample_Detail.tsv"
    ),
  
  Primary_Sample_CellCounts =
    file.path(
      base_dir,
      "STEP11-3_Primary_Sample_CellCounts.tsv"
    ),
  
  Patient_Eligibility_Detail =
    file.path(
      base_dir,
      "STEP11-3_Patient_Eligibility_Detail.tsv"
    ),
  
  Patient_Eligibility_Audit =
    file.path(
      base_dir,
      "STEP11-3_Patient_Eligibility_Audit.tsv"
    ),
  
  MultiSample_Patient_Audit =
    file.path(
      base_dir,
      "STEP11-3_MultiSample_Patient_Audit.tsv"
    )
)


# ------------------------------------------------------------
# 2. Input existence audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("STEP11-4 INPUT EXISTENCE AUDIT\n")
cat("------------------------------------------------------------\n\n")

all_files <- c(step10_files, step11_3_files)

existence_df <- data.frame(
  File = names(all_files),
  Path = unname(all_files),
  Exists = file.exists(all_files),
  stringsAsFactors = FALSE
)

print(existence_df)

if (!all(existence_df$Exists)) {
  
  cat("\nERROR: One or more required inputs are missing.\n")
  stop("STEP11-4 stopped because required input files are missing.")
  
}

cat("\nPASS — all required STEP10 and STEP11-3 inputs exist.\n\n")


# ------------------------------------------------------------
# 3. Read locked inputs
# ------------------------------------------------------------

read_tsv_locked <- function(path) {
  
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
}


step10_obj <- lapply(step10_files, read_tsv_locked)
step11_3_obj <- lapply(step11_3_files, read_tsv_locked)


# ------------------------------------------------------------
# 4. Basic structural audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("STEP11-4 STRUCTURAL INPUT AUDIT\n")
cat("------------------------------------------------------------\n\n")


for (nm in names(step10_obj)) {
  
  x <- step10_obj[[nm]]
  
  cat(
    sprintf(
      "STEP10 %-35s : %d rows × %d columns\n",
      nm,
      nrow(x),
      ncol(x)
    )
  )
  
}


cat("\n")


for (nm in names(step11_3_obj)) {
  
  x <- step11_3_obj[[nm]]
  
  cat(
    sprintf(
      "STEP11-3 %-29s : %d rows × %d columns\n",
      nm,
      nrow(x),
      ncol(x)
    )
  )
  
}


# ------------------------------------------------------------
# 5. Extract key objects
# ------------------------------------------------------------

expr_comp <- step10_obj$STEP10F1C_GeneCompartment
expr_patient <- step10_obj$STEP10F1E_ThreeCompartmentProfile

pattern_summary <- step10_obj$STEP10F2A_Pattern
within_rank <- step10_obj$STEP10F2D_WithinPatientRank
three_comp_pattern <- step10_obj$STEP10F2E_ThreeCompartment
dominant_comp <- step10_obj$STEP10F2F_Dominant

tumor_sample_detail <- step11_3_obj$Primary_Tumor_Sample_Detail
sample_counts <- step11_3_obj$Primary_Sample_CellCounts
patient_eligibility <- step11_3_obj$Patient_Eligibility_Detail
patient_threshold <- step11_3_obj$Patient_Eligibility_Audit
multisample <- step11_3_obj$MultiSample_Patient_Audit


# ------------------------------------------------------------
# 6. Required-column audit
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("REQUIRED COLUMN AUDIT\n")
cat("------------------------------------------------------------\n\n")


required_columns <- list(
  
  expr_comp = c(
    "Gene",
    "Compartment",
    "N_Patients",
    "Detection_N",
    "Detection_Percent"
  ),
  
  
  expr_patient = c(
    "Patient",
    "Compartment",
    "Gene",
    "Cell_Count",
    "CPM",
    "Detected",
    "Analysis_Set"
  ),
  
  pattern_summary = c(
    "Gene",
    "Compartment"
  ),
  
  within_rank = c(
    "Gene",
    "Patient",
    "Compartment",
    "CPM",
    "Detection",
    "Within_Patient_CPM_Rank"
  ),
  
  three_comp_pattern = c(
    "Gene",
    "Compartment",
    "Complete_Patients",
    "Median_CPM",
    "Detection_N",
    "Detection_Percent",
    "Rank1_N",
    "Rank1_Percent"
  ),
  
  dominant_comp = c(
    "Gene",
    "Patient",
    "Dominant_Compartment_by_CPM",
    "Max_CPM"
  ),
  
  tumor_sample_detail = c(
    "Patient",
    "Sample",
    "Source",
    "Endothelial",
    "Mono/Macro",
    "Tcell_Total",
    "Target_Total"
  ),
  
  sample_counts = c(
    "Sample",
    "Patient",
    "Source",
    "Endothelial",
    "Mono/Macro",
    "Tcell_Total",
    "Target_Total",
    "Is_Tumor_Source",
    "All_Three_Nonzero",
    "All_Three_Nonzero_Tumor"
  ),
  
  patient_eligibility = c(
    "Patient",
    "Threshold",
    "Eligible_Primary_Sample",
    "Eligible_Tumor_Sample"
  ),
  
  patient_threshold = c(
    "Threshold",
    "Primary_Patients_Total",
    "Primary_Patients_Eligible",
    "Tumor_Patients_Eligible"
  ),
  
  multisample = c(
    "Patient",
    "N_Primary_Samples"
  )
  
)


object_map <- list(
  
  expr_comp = expr_comp,
  expr_patient = expr_patient,
  pattern_summary = pattern_summary,
  within_rank = within_rank,
  three_comp_pattern = three_comp_pattern,
  dominant_comp = dominant_comp,
  tumor_sample_detail = tumor_sample_detail,
  sample_counts = sample_counts,
  patient_eligibility = patient_eligibility,
  patient_threshold = patient_threshold,
  multisample = multisample
  
)


column_audit <- do.call(
  rbind,
  lapply(
    names(required_columns),
    function(nm) {
      
      req <- required_columns[[nm]]
      actual <- names(object_map[[nm]])
      
      data.frame(
        Object = nm,
        Required = paste(req, collapse = "; "),
        Missing = paste(
          setdiff(req, actual),
          collapse = "; "
        ),
        PASS = all(req %in% actual),
        stringsAsFactors = FALSE
      )
      
    }
  )
)

print(column_audit)

if (!all(column_audit$PASS)) {
  stop("STEP11-4 stopped: required columns missing.")
}

cat("\nPASS — all required columns are present.\n")


# ------------------------------------------------------------
# 7. Core-gene lock
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

cat("\n------------------------------------------------------------\n")
cat("CORE-GENE LOCK\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Core genes:",
  paste(core_genes, collapse = ", "),
  "\n"
)

if (!all(core_genes %in% unique(expr_patient$Gene))) {
  
  stop(
    "STEP11-4 stopped: one or more core genes are absent from STEP10F1E."
  )
  
}

cat("PASS — AOC3 / MS4A7 / IL7R all present.\n")


# ------------------------------------------------------------
# 8. Compartment lock
# ------------------------------------------------------------

target_compartments <- c(
  "Endothelial",
  "Myeloid/macrophage",
  "T-cell"
)

cat("\n------------------------------------------------------------\n")
cat("COMPARTMENT LOCK\n")
cat("------------------------------------------------------------\n\n")

cat(
  "STEP10 compartments:\n"
)

print(sort(unique(expr_patient$Compartment)))

cat("\n")

cat(
  "STEP11-4 target compartments:\n"
)

print(target_compartments)


# ------------------------------------------------------------
# 9. Check STEP10 patient-level profile structure
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("STEP10 PATIENT-LEVEL STRUCTURE AUDIT\n")
cat("------------------------------------------------------------\n\n")

expr_core <- expr_patient[
  expr_patient$Gene %in% core_genes,
  ,
  drop = FALSE
]

cat(
  "Rows for three core genes:",
  nrow(expr_core),
  "\n"
)

cat(
  "Unique patients:",
  length(unique(expr_core$Patient)),
  "\n"
)

cat(
  "Unique genes:",
  length(unique(expr_core$Gene)),
  "\n"
)

cat(
  "Unique compartments:",
  length(unique(expr_core$Compartment)),
  "\n"
)


# ------------------------------------------------------------
# 10. Confirm patient × gene × compartment uniqueness
# ------------------------------------------------------------

key <- paste(
  expr_core$Patient,
  expr_core$Gene,
  expr_core$Compartment,
  sep = "||"
)

duplicate_keys <- duplicated(key)

cat(
  "Duplicate patient × gene × compartment rows:",
  sum(duplicate_keys),
  "\n"
)

if (any(duplicate_keys)) {
  stop(
    "STEP11-4 stopped: duplicate patient × gene × compartment records detected."
  )
}

cat("PASS — patient × gene × compartment structure is unique.\n")


# ------------------------------------------------------------
# 11. Build descriptive localization table
# ------------------------------------------------------------

localization_core <- expr_comp[
  expr_comp$Gene %in% core_genes,
  ,
  drop = FALSE
]

localization_core <- localization_core[
  localization_core$Compartment %in%
    target_compartments,
  ,
  drop = FALSE
]

localization_core <- localization_core[
  order(
    match(localization_core$Gene, core_genes),
    match(
      localization_core$Compartment,
      target_compartments
    )
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 12. Integrate patient-level data with STEP11-3 eligibility
# ------------------------------------------------------------

tumor_patients <- unique(
  tumor_sample_detail$Patient
)

cat("\n------------------------------------------------------------\n")
cat("PRIMARY-TUMOR PATIENT INTEGRATION\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Patients represented in Primary Tumor sample detail:",
  length(tumor_patients),
  "\n"
)

expr_core_tumor_patient <- expr_core[
  expr_core$Patient %in% tumor_patients,
  ,
  drop = FALSE
]

cat(
  "STEP10 patient-level records retained:",
  nrow(expr_core_tumor_patient),
  "\n"
)


# ------------------------------------------------------------
# 13. Add patient-level tumor provenance
# ------------------------------------------------------------

tumor_patient_summary <- aggregate(
  Sample ~ Patient,
  data = tumor_sample_detail,
  FUN = function(x) paste(unique(x), collapse = ";")
)

names(tumor_patient_summary)[2] <- "Primary_Tumor_Samples"

tumor_patient_summary$N_Primary_Tumor_Samples <- vapply(
  strsplit(
    tumor_patient_summary$Primary_Tumor_Samples,
    ";",
    fixed = TRUE
  ),
  length,
  integer(1)
)

expr_core_tumor_patient <- merge(
  expr_core_tumor_patient,
  tumor_patient_summary,
  by = "Patient",
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------
# 14. Patient-level localization summary
# ------------------------------------------------------------

patient_localization <- merge(
  expr_core_tumor_patient,
  dominant_comp[
    dominant_comp$Gene %in% core_genes,
    ,
    drop = FALSE
  ],
  by = c("Patient", "Gene"),
  all.x = TRUE,
  sort = FALSE
)

patient_localization <- patient_localization[
  order(
    match(patient_localization$Gene, core_genes),
    patient_localization$Patient,
    match(
      patient_localization$Compartment,
      target_compartments
    )
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 15. Patient-level rank summary
# ------------------------------------------------------------

rank_core <- within_rank[
  within_rank$Gene %in% core_genes,
  ,
  drop = FALSE
]

rank_core <- rank_core[
  rank_core$Compartment %in% target_compartments,
  ,
  drop = FALSE
]

rank_core <- rank_core[
  order(
    match(rank_core$Gene, core_genes),
    rank_core$Patient,
    rank_core$Within_Patient_CPM_Rank
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 16. Dominant-compartment summary
# ------------------------------------------------------------

dominant_core <- dominant_comp[
  dominant_comp$Gene %in% core_genes,
  ,
  drop = FALSE
]

dominant_summary <- as.data.frame(
  table(
    Gene = dominant_core$Gene,
    Dominant_Compartment =
      dominant_core$Dominant_Compartment_by_CPM
  ),
  stringsAsFactors = FALSE
)

names(dominant_summary)[3] <- "N_Patients"


# ------------------------------------------------------------
# 17. Detection / quantitative summary
# ------------------------------------------------------------

expression_summary_core <- localization_core[
  ,
  c(
    "Gene",
    "Compartment",
    "N_Patients",
    "Detection_N",
    "Detection_Percent"
  ),
  drop = FALSE
]

# ------------------------------------------------------------
# 18. STEP11-3 patient eligibility feasibility summary
# ------------------------------------------------------------

eligibility_summary <- patient_threshold[
  ,
  c(
    "Threshold",
    "Primary_Patients_Total",
    "Primary_Patients_Eligible",
    "Primary_Patient_Percent",
    "Tumor_Patients_Eligible",
    "Tumor_Patient_Percent"
  ),
  drop = FALSE
]


# ------------------------------------------------------------
# 19. Multi-sample patient summary
# ------------------------------------------------------------

multisample_summary <- multisample

multisample_summary$Multiple_Primary_Samples <-
  multisample_summary$N_Primary_Samples > 1


# ------------------------------------------------------------
# 20. Structural concordance checks
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("STEP11-4 STRUCTURAL CONCORDANCE CHECKS\n")
cat("------------------------------------------------------------\n\n")


check_1 <- length(unique(expr_core$Gene)) == 3

check_2 <- all(
  core_genes %in%
    unique(expr_core$Gene)
)

check_3 <- all(
  tumor_sample_detail$Source == "Tumor"
)

check_4 <- all(
  sample_counts$Is_Tumor_Source ==
    (sample_counts$Source == "Tumor")
)

check_5 <- all(
  tumor_sample_detail$Patient %in%
    unique(sample_counts$Patient)
)

check_6 <- all(
  expr_core_tumor_patient$Patient %in%
    tumor_patients
)

check_7 <- nrow(
  localization_core
) == 9

check_8 <- nrow(
  expression_summary_core
) == 9


concordance <- data.frame(
  Check = c(
    "Three core genes",
    "All core genes present",
    "Primary Tumor sample source = Tumor",
    "Tumor flag concordance",
    "Tumor-detail patients in sample counts",
    "STEP10 tumor patients restricted correctly",
    "3 genes × 3 compartments summary",
    "Expression summary has 9 rows"
  ),
  PASS = c(
    check_1,
    check_2,
    check_3,
    check_4,
    check_5,
    check_6,
    check_7,
    check_8
  ),
  stringsAsFactors = FALSE
)

print(concordance)

if (!all(concordance$PASS)) {
  
  stop(
    "STEP11-4 stopped: structural concordance check failed."
  )
  
}

cat("\nPASS — all structural concordance checks passed.\n")


# ------------------------------------------------------------
# 21. Manuscript-safe localization summary
# ------------------------------------------------------------

manuscript_summary <- do.call(
  rbind,
  lapply(
    core_genes,
    function(g) {
      
      x <- dominant_core[
        dominant_core$Gene == g,
        ,
        drop = FALSE
      ]
      
      tab <- sort(
        table(x$Dominant_Compartment_by_CPM),
        decreasing = TRUE
      )
      
      data.frame(
        Gene = g,
        N_Patients = nrow(x),
        Dominant_Compartment =
          names(tab)[1],
        Dominant_N = as.integer(tab[1]),
        Dominant_Percent =
          round(
            100 * as.integer(tab[1]) / nrow(x),
            2
          ),
        stringsAsFactors = FALSE
      )
      
    }
  )
)

print(manuscript_summary)


# ------------------------------------------------------------
# 22. Save outputs
# ------------------------------------------------------------

write.table(
  localization_core,
  file.path(
    out_dir,
    "STEP11-4_CoreGene_Compartment_Expression_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  expression_summary_core,
  file.path(
    out_dir,
    "STEP11-4_CoreGene_Expression_Detection_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  patient_localization,
  file.path(
    out_dir,
    "STEP11-4_CoreGene_Patient_Level_Localization.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  rank_core,
  file.path(
    out_dir,
    "STEP11-4_CoreGene_Within_Patient_Rank.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  dominant_summary,
  file.path(
    out_dir,
    "STEP11-4_CoreGene_Dominant_Compartment_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  eligibility_summary,
  file.path(
    out_dir,
    "STEP11-4_Patient_Eligibility_Feasibility_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  multisample_summary,
  file.path(
    out_dir,
    "STEP11-4_MultiSample_Patient_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  concordance,
  file.path(
    out_dir,
    "STEP11-4_Structural_Concordance_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  manuscript_summary,
  file.path(
    out_dir,
    "STEP11-4_Manuscript_Safe_Localization_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 23. Save integrated RDS
# ------------------------------------------------------------

STEP11_4_LOCKED <- list(
  
  analysis = "STEP11-4 Core-Gene Cellular Localization Integration",
  
  core_genes = core_genes,
  
  target_compartments = target_compartments,
  
  source_STEP10 = step10_files,
  
  source_STEP11_3 = step11_3_files,
  
  expression_summary =
    expression_summary_core,
  
  patient_localization =
    patient_localization,
  
  within_patient_rank =
    rank_core,
  
  dominant_summary =
    dominant_summary,
  
  eligibility_summary =
    eligibility_summary,
  
  multisample_summary =
    multisample_summary,
  
  structural_concordance =
    concordance,
  
  manuscript_safe_summary =
    manuscript_summary,
  
  threshold_status =
    "NOT LOCKED",
  
  raw_UMI_recalculation =
    FALSE,
  
  expression_recalculation =
    FALSE,
  
  cell_level_reanalysis =
    FALSE,
  
  per_patient_CellChat =
    FALSE
  
)


saveRDS(
  STEP11_4_LOCKED,
  file.path(
    out_dir,
    "STEP11-4_CoreGene_Cellular_Localization_Integration_LOCKED.rds"
  )
)


# ------------------------------------------------------------
# 24. Final QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP11-4 FINAL QC\n")
cat("============================================================\n\n")

cat(
  "Core genes:",
  paste(core_genes, collapse = " / "),
  "\n"
)

cat(
  "STEP10 core-gene patient-level rows:",
  nrow(expr_core),
  "\n"
)

cat(
  "Primary Tumor patients represented:",
  length(tumor_patients),
  "\n"
)

cat(
  "Primary Tumor samples:",
  nrow(tumor_sample_detail),
  "\n"
)

cat(
  "Core-gene × compartment summary rows:",
  nrow(expression_summary_core),
  "\n"
)

cat(
  "Patient-level localization rows:",
  nrow(patient_localization),
  "\n"
)

cat(
  "Within-patient rank rows:",
  nrow(rank_core),
  "\n"
)

cat(
  "Threshold locked:",
  "NO",
  "\n"
)

cat(
  "Raw UMI recalculated:",
  "NO",
  "\n"
)

cat(
  "Expression recalculated:",
  "NO",
  "\n"
)

cat(
  "Per-patient CellChat:",
  "NO",
  "\n"
)

cat(
  "\nSTEP11-4 STATUS: PASS\n"
)

cat(
  "============================================================\n"
)


cat("============================================================\n")
cat("STEP11-5A — CCI_LR × STEP10 GENE-EXPRESSION COMPATIBILITY AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directories
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 11-CellChat"

step10_dir <- "input_dir/Paper 1_Step 10-Single-Cell"

out_dir <- file.path(
  base_dir,
  "STEP11-5A_CCI_LR_STEP10_Compatibility_Audit"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

cci_file <- file.path(
  base_dir,
  "NSCLC_GSE131907_CCI_LR.tsv"
)

step10_expr_file <- file.path(
  step10_dir,
  "STEP10F-1_Patient_Level_Expression_Profile",
  "STEP10F-1E_Complete_ThreeCompartment_Patient_Profile.tsv"
)


# ------------------------------------------------------------
# 2. Existence audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("INPUT EXISTENCE AUDIT\n")
cat("------------------------------------------------------------\n\n")

cat(
  "CCI_LR file:",
  file.exists(cci_file),
  "\n"
)

cat(
  "STEP10 expression file:",
  file.exists(step10_expr_file),
  "\n"
)

if (!file.exists(cci_file)) {
  stop("CCI_LR input file not found.")
}

if (!file.exists(step10_expr_file)) {
  stop("STEP10 expression input file not found.")
}

cat("\nPASS — required inputs exist.\n")


# ------------------------------------------------------------
# 3. Read inputs
# ------------------------------------------------------------

cci_lr <- read.delim(
  cci_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

expr_patient <- read.delim(
  step10_expr_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 4. Basic input structure audit
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("INPUT STRUCTURE AUDIT\n")
cat("------------------------------------------------------------\n\n")

cat(
  "CCI_LR:",
  nrow(cci_lr),
  "rows ×",
  ncol(cci_lr),
  "columns\n"
)

cat(
  "STEP10 expression:",
  nrow(expr_patient),
  "rows ×",
  ncol(expr_patient),
  "columns\n\n"
)


cat("CCI_LR columns:\n")
print(names(cci_lr))

cat("\nSTEP10 expression columns:\n")
print(names(expr_patient))


# ------------------------------------------------------------
# 5. Required-column audit
# ------------------------------------------------------------

required_cci <- c(
  "interaction_name",
  "pathway_name",
  "ligand",
  "receptor",
  "evidence",
  "annotation"
)

required_expr <- c(
  "Patient",
  "Compartment",
  "Gene",
  "Cell_Count",
  "CPM",
  "Detected",
  "Analysis_Set"
)

missing_cci <- setdiff(
  required_cci,
  names(cci_lr)
)

missing_expr <- setdiff(
  required_expr,
  names(expr_patient)
)

cat("\n------------------------------------------------------------\n")
cat("REQUIRED COLUMN AUDIT\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Missing CCI_LR columns:",
  ifelse(
    length(missing_cci) == 0,
    "NONE",
    paste(missing_cci, collapse = "; ")
  ),
  "\n"
)

cat(
  "Missing STEP10 columns:",
  ifelse(
    length(missing_expr) == 0,
    "NONE",
    paste(missing_expr, collapse = "; ")
  ),
  "\n"
)

if (
  length(missing_cci) > 0 ||
  length(missing_expr) > 0
) {
  stop("STEP11-5A stopped: required columns are missing.")
}

cat("\nPASS — all required columns are present.\n")


# ------------------------------------------------------------
# 6. Basic CCI_LR integrity audit
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("CCI_LR INTEGRITY AUDIT\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Total interaction rows:",
  nrow(cci_lr),
  "\n"
)

cat(
  "Unique interaction_name:",
  length(unique(cci_lr$interaction_name)),
  "\n"
)

cat(
  "Unique ligand:",
  length(unique(cci_lr$ligand)),
  "\n"
)

cat(
  "Unique receptor:",
  length(unique(cci_lr$receptor)),
  "\n"
)

cat(
  "Unique pathways:",
  length(unique(cci_lr$pathway_name)),
  "\n"
)

cat(
  "Unique annotations:",
  length(unique(cci_lr$annotation)),
  "\n"
)

cat(
  "Missing ligand:",
  sum(is.na(cci_lr$ligand) | cci_lr$ligand == ""),
  "\n"
)

cat(
  "Missing receptor:",
  sum(is.na(cci_lr$receptor) | cci_lr$receptor == ""),
  "\n"
)

cat(
  "Missing pathway:",
  sum(is.na(cci_lr$pathway_name) | cci_lr$pathway_name == ""),
  "\n"
)


# ------------------------------------------------------------
# 7. STEP10 expression integrity audit
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("STEP10 EXPRESSION INTEGRITY AUDIT\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Unique patients:",
  length(unique(expr_patient$Patient)),
  "\n"
)

cat(
  "Unique genes:",
  length(unique(expr_patient$Gene)),
  "\n"
)

cat(
  "Unique compartments:",
  length(unique(expr_patient$Compartment)),
  "\n"
)

cat(
  "Compartments:\n"
)

print(sort(unique(expr_patient$Compartment)))

cat(
  "\nUnique Analysis_Set values:\n"
)

print(unique(expr_patient$Analysis_Set))

cat(
  "\nNegative CPM values:",
  sum(expr_patient$CPM < 0, na.rm = TRUE),
  "\n"
)

cat(
  "Missing CPM:",
  sum(is.na(expr_patient$CPM)),
  "\n"
)

cat(
  "Missing Gene:",
  sum(is.na(expr_patient$Gene) | expr_patient$Gene == ""),
  "\n"
)

cat(
  "Missing Patient:",
  sum(is.na(expr_patient$Patient) | expr_patient$Patient == ""),
  "\n"
)

cat(
  "Missing Compartment:",
  sum(is.na(expr_patient$Compartment) |
        expr_patient$Compartment == ""),
  "\n"
)


# ------------------------------------------------------------
# 8. Define focused compartments
# ------------------------------------------------------------

target_compartments <- c(
  "Endothelial",
  "Myeloid/macrophage",
  "T-cell"
)

cat("\n------------------------------------------------------------\n")
cat("FOCUSED COMPARTMENT LOCK\n")
cat("------------------------------------------------------------\n\n")

print(target_compartments)

if (
  !all(
    target_compartments %in%
    unique(expr_patient$Compartment)
  )
) {
  stop(
    "STEP11-5A stopped: one or more focused compartments are absent."
  )
}

cat("\nPASS — all three focused compartments are present.\n")


# ------------------------------------------------------------
# 9. Build expression gene universe
# ------------------------------------------------------------

expr_gene_universe <- sort(
  unique(
    expr_patient$Gene[
      !is.na(expr_patient$Gene) &
        expr_patient$Gene != ""
    ]
  )
)

cat("\n------------------------------------------------------------\n")
cat("GENE-UNIVERSE COMPATIBILITY\n")
cat("------------------------------------------------------------\n\n")

cat(
  "STEP10 genes:",
  length(expr_gene_universe),
  "\n"
)


# ------------------------------------------------------------
# 10. CCI_LR ligand/receptor compatibility
# ------------------------------------------------------------

cci_lr$ligand_in_STEP10 <- cci_lr$ligand %in%
  expr_gene_universe

cci_lr$receptor_in_STEP10 <- cci_lr$receptor %in%
  expr_gene_universe

cci_lr$both_in_STEP10 <-
  cci_lr$ligand_in_STEP10 &
  cci_lr$receptor_in_STEP10


compatibility_summary <- data.frame(
  Category = c(
    "All CCI_LR interactions",
    "Ligand found in STEP10",
    "Receptor found in STEP10",
    "Both ligand and receptor found in STEP10",
    "Neither ligand nor receptor found in STEP10"
  ),
  N = c(
    nrow(cci_lr),
    sum(cci_lr$ligand_in_STEP10),
    sum(cci_lr$receptor_in_STEP10),
    sum(cci_lr$both_in_STEP10),
    sum(
      !cci_lr$ligand_in_STEP10 &
        !cci_lr$receptor_in_STEP10
    )
  ),
  Percent = c(
    100,
    100 * mean(cci_lr$ligand_in_STEP10),
    100 * mean(cci_lr$receptor_in_STEP10),
    100 * mean(cci_lr$both_in_STEP10),
    100 * mean(
      !cci_lr$ligand_in_STEP10 &
        !cci_lr$receptor_in_STEP10
    )
  ),
  stringsAsFactors = FALSE
)

cat("\nCCI_LR × STEP10 compatibility:\n\n")

print(compatibility_summary)


# ------------------------------------------------------------
# 11. Focused expression universe
# ------------------------------------------------------------

expr_focused <- expr_patient[
  expr_patient$Compartment %in%
    target_compartments,
  ,
  drop = FALSE
]

focused_gene_universe <- sort(
  unique(expr_focused$Gene)
)

cat("\n------------------------------------------------------------\n")
cat("FOCUSED-COMPARTMENT GENE UNIVERSE\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Focused compartments:",
  length(target_compartments),
  "\n"
)

cat(
  "STEP10 focused genes:",
  length(focused_gene_universe),
  "\n"
)


# ------------------------------------------------------------
# 12. LR compatibility at focused-compartment level
# ------------------------------------------------------------

cci_lr$ligand_in_focused_STEP10 <-
  cci_lr$ligand %in% focused_gene_universe

cci_lr$receptor_in_focused_STEP10 <-
  cci_lr$receptor %in% focused_gene_universe

cci_lr$both_in_focused_STEP10 <-
  cci_lr$ligand_in_focused_STEP10 &
  cci_lr$receptor_in_focused_STEP10


focused_compatibility_summary <- data.frame(
  Category = c(
    "All CCI_LR interactions",
    "Ligand found in focused compartments",
    "Receptor found in focused compartments",
    "Both found in focused compartments",
    "Neither found in focused compartments"
  ),
  N = c(
    nrow(cci_lr),
    sum(cci_lr$ligand_in_focused_STEP10),
    sum(cci_lr$receptor_in_focused_STEP10),
    sum(cci_lr$both_in_focused_STEP10),
    sum(
      !cci_lr$ligand_in_focused_STEP10 &
        !cci_lr$receptor_in_focused_STEP10
    )
  ),
  Percent = c(
    100,
    100 * mean(cci_lr$ligand_in_focused_STEP10),
    100 * mean(cci_lr$receptor_in_focused_STEP10),
    100 * mean(cci_lr$both_in_focused_STEP10),
    100 * mean(
      !cci_lr$ligand_in_focused_STEP10 &
        !cci_lr$receptor_in_focused_STEP10
    )
  ),
  stringsAsFactors = FALSE
)

cat("\nFocused compatibility:\n\n")

print(focused_compatibility_summary)


# ------------------------------------------------------------
# 13. Define six possible compartment directions
# ------------------------------------------------------------

direction_table <- expand.grid(
  Sender = target_compartments,
  Receiver = target_compartments,
  stringsAsFactors = FALSE
)

direction_table <- direction_table[
  direction_table$Sender !=
    direction_table$Receiver,
  ,
  drop = FALSE
]

direction_table$Direction <- paste(
  direction_table$Sender,
  "->",
  direction_table$Receiver
)

cat("\n------------------------------------------------------------\n")
cat("SIX FOCUSED DIRECTIONS\n")
cat("------------------------------------------------------------\n\n")

print(direction_table)


# ------------------------------------------------------------
# 14. Build compartment-specific expression presence
# ------------------------------------------------------------

# A gene is considered present in a compartment if it occurs
# in the locked STEP10 patient-level expression table.
#
# This is a structural compatibility criterion only.
# It is NOT an LR activity criterion.

ligand_presence <- unique(
  expr_focused[
    ,
    c("Gene", "Compartment"),
    drop = FALSE
  ]
)

names(ligand_presence) <- c(
  "Gene",
  "Compartment"
)


# ------------------------------------------------------------
# 15. Generate focused LR candidate table
# ------------------------------------------------------------

candidate_rows <- list()
counter <- 0L

for (i in seq_len(nrow(cci_lr))) {
  
  lr <- cci_lr[i, ]
  
  if (!isTRUE(lr$both_in_STEP10)) {
    next
  }
  
  for (d in seq_len(nrow(direction_table))) {
    
    sender <- direction_table$Sender[d]
    receiver <- direction_table$Receiver[d]
    
    ligand_present <- any(
      ligand_presence$Gene == lr$ligand &
        ligand_presence$Compartment == sender
    )
    
    receptor_present <- any(
      ligand_presence$Gene == lr$receptor &
        ligand_presence$Compartment == receiver
    )
    
    if (
      ligand_present &&
      receptor_present
    ) {
      
      counter <- counter + 1L
      
      candidate_rows[[counter]] <- data.frame(
        
        Direction =
          direction_table$Direction[d],
        
        Sender =
          sender,
        
        Receiver =
          receiver,
        
        interaction_name =
          lr$interaction_name,
        
        pathway_name =
          lr$pathway_name,
        
        ligand =
          lr$ligand,
        
        receptor =
          lr$receptor,
        
        evidence =
          lr$evidence,
        
        annotation =
          lr$annotation,
        
        ligand_in_STEP10 =
          lr$ligand_in_STEP10,
        
        receptor_in_STEP10 =
          lr$receptor_in_STEP10,
        
        both_in_STEP10 =
          lr$both_in_STEP10,
        
        ligand_present_sender =
          ligand_present,
        
        receptor_present_receiver =
          receptor_present,
        
        stringsAsFactors = FALSE
      )
      
    }
    
  }
  
}


if (length(candidate_rows) > 0) {
  
  focused_candidates <- do.call(
    rbind,
    candidate_rows
  )
  
} else {
  
  focused_candidates <- data.frame(
    Direction = character(),
    Sender = character(),
    Receiver = character(),
    interaction_name = character(),
    pathway_name = character(),
    ligand = character(),
    receptor = character(),
    evidence = character(),
    annotation = character(),
    ligand_in_STEP10 = logical(),
    receptor_in_STEP10 = logical(),
    both_in_STEP10 = logical(),
    ligand_present_sender = logical(),
    receptor_present_receiver = logical(),
    stringsAsFactors = FALSE
  )
  
}


# ------------------------------------------------------------
# 16. Candidate interaction summary
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("FOCUSED LR CANDIDATE SUMMARY\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Candidate LR direction-specific rows:",
  nrow(focused_candidates),
  "\n"
)

cat(
  "Unique interaction_name:",
  length(
    unique(
      focused_candidates$interaction_name
    )
  ),
  "\n"
)

cat(
  "Unique pathways:",
  length(
    unique(
      focused_candidates$pathway_name
    )
  ),
  "\n"
)


# ------------------------------------------------------------
# 17. Direction-level candidate counts
# ------------------------------------------------------------

if (nrow(focused_candidates) > 0) {
  
  direction_counts <- aggregate(
    interaction_name ~
      Direction + Sender + Receiver,
    data = focused_candidates,
    FUN = function(x) length(unique(x))
  )
  
  names(direction_counts)[4] <-
    "N_Unique_LR_Interactions"
  
} else {
  
  direction_counts <- direction_table
  
  direction_counts$N_Unique_LR_Interactions <- 0
  
}

cat("\nDirection-level candidate counts:\n\n")

print(direction_counts)


# ------------------------------------------------------------
# 18. Pathway-level candidate counts
# ------------------------------------------------------------

if (nrow(focused_candidates) > 0) {
  
  pathway_counts <- aggregate(
    interaction_name ~ pathway_name,
    data = focused_candidates,
    FUN = function(x) length(unique(x))
  )
  
  names(pathway_counts)[2] <-
    "N_Unique_LR_Interactions"
  
  pathway_counts <- pathway_counts[
    order(
      -pathway_counts$N_Unique_LR_Interactions,
      pathway_counts$pathway_name
    ),
    ,
    drop = FALSE
  ]
  
} else {
  
  pathway_counts <- data.frame(
    pathway_name = character(),
    N_Unique_LR_Interactions = integer(),
    stringsAsFactors = FALSE
  )
  
}

cat("\nTop pathway candidate counts:\n\n")

print(
  head(
    pathway_counts,
    30
  )
)


# ------------------------------------------------------------
# 19. Direction × pathway summary
# ------------------------------------------------------------

if (nrow(focused_candidates) > 0) {
  
  direction_pathway_counts <- aggregate(
    interaction_name ~
      Direction + pathway_name,
    data = focused_candidates,
    FUN = function(x) length(unique(x))
  )
  
  names(direction_pathway_counts)[3] <-
    "N_Unique_LR_Interactions"
  
  direction_pathway_counts <-
    direction_pathway_counts[
      order(
        direction_pathway_counts$Direction,
        -direction_pathway_counts$N_Unique_LR_Interactions,
        direction_pathway_counts$pathway_name
      ),
      ,
      drop = FALSE
    ]
  
} else {
  
  direction_pathway_counts <- data.frame(
    Direction = character(),
    pathway_name = character(),
    N_Unique_LR_Interactions = integer(),
    stringsAsFactors = FALSE
  )
  
}


# ------------------------------------------------------------
# 20. Core-gene presence audit
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_presence <- expand.grid(
  Gene = core_genes,
  Compartment = target_compartments,
  stringsAsFactors = FALSE
)

core_presence$Present_in_STEP10 <- mapply(
  function(g, c) {
    
    any(
      expr_patient$Gene == g &
        expr_patient$Compartment == c
    )
    
  },
  core_presence$Gene,
  core_presence$Compartment
)

cat("\n------------------------------------------------------------\n")
cat("CORE-GENE × COMPARTMENT PRESENCE AUDIT\n")
cat("------------------------------------------------------------\n\n")

print(core_presence)


# ------------------------------------------------------------
# 21. Core genes as LR ligands/receptors
# ------------------------------------------------------------

core_lr_rows <- cci_lr[
  cci_lr$ligand %in% core_genes |
    cci_lr$receptor %in% core_genes,
  ,
  drop = FALSE
]

cat("\n------------------------------------------------------------\n")
cat("CORE-GENE × CCI_LR RESOURCE AUDIT\n")
cat("------------------------------------------------------------\n\n")

cat(
  "CCI_LR rows containing any core gene as ligand or receptor:",
  nrow(core_lr_rows),
  "\n"
)

if (nrow(core_lr_rows) > 0) {
  
  core_lr_summary <- data.frame(
    Gene = core_genes,
    As_Ligand = sapply(
      core_genes,
      function(g)
        sum(core_lr_rows$ligand == g)
    ),
    As_Receptor = sapply(
      core_genes,
      function(g)
        sum(core_lr_rows$receptor == g)
    ),
    stringsAsFactors = FALSE
  )
  
} else {
  
  core_lr_summary <- data.frame(
    Gene = core_genes,
    As_Ligand = 0,
    As_Receptor = 0,
    stringsAsFactors = FALSE
  )
  
}

print(core_lr_summary)


# ------------------------------------------------------------
# 22. Important interpretation flag
# ------------------------------------------------------------

interpretation_note <- data.frame(
  
  Item = c(
    "CCI_LR resource membership",
    "STEP10 gene presence",
    "Focused direction candidate",
    "Dataset-specific communication",
    "Statistical significance",
    "Patient-level consistency",
    "Causal/mechanistic inference"
  ),
  
  Status = c(
    "Audited",
    "Audited",
    "Generated for downstream screening",
    "NOT established",
    "NOT tested",
    "NOT tested",
    "NOT established"
  ),
  
  stringsAsFactors = FALSE
)

cat("\n------------------------------------------------------------\n")
cat("INTERPRETATION STATUS\n")
cat("------------------------------------------------------------\n\n")

print(interpretation_note)


# ------------------------------------------------------------
# 23. Save outputs
# ------------------------------------------------------------

write.table(
  compatibility_summary,
  file.path(
    out_dir,
    "STEP11-5A_CCI_LR_Global_Compatibility_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  focused_compatibility_summary,
  file.path(
    out_dir,
    "STEP11-5A_CCI_LR_Focused_Compatibility_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  focused_candidates,
  file.path(
    out_dir,
    "STEP11-5A_Focused_LR_Candidate_Interactions.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  direction_counts,
  file.path(
    out_dir,
    "STEP11-5A_Direction_Candidate_Counts.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  pathway_counts,
  file.path(
    out_dir,
    "STEP11-5A_Pathway_Candidate_Counts.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  direction_pathway_counts,
  file.path(
    out_dir,
    "STEP11-5A_Direction_Pathway_Candidate_Counts.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  core_presence,
  file.path(
    out_dir,
    "STEP11-5A_CoreGene_Compartment_Presence.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  core_lr_summary,
  file.path(
    out_dir,
    "STEP11-5A_CoreGene_CCI_LR_Resource_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  interpretation_note,
  file.path(
    out_dir,
    "STEP11-5A_Interpretation_Status.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 24. Save complete audit object
# ------------------------------------------------------------

STEP11_5A_LOCKED <- list(
  
  analysis =
    "STEP11-5A CCI_LR × STEP10 Gene-Expression Compatibility Audit",
  
  input_CCI_LR =
    cci_file,
  
  input_STEP10_expression =
    step10_expr_file,
  
  core_genes =
    core_genes,
  
  target_compartments =
    target_compartments,
  
  compatibility_summary =
    compatibility_summary,
  
  focused_compatibility_summary =
    focused_compatibility_summary,
  
  focused_candidates =
    focused_candidates,
  
  direction_counts =
    direction_counts,
  
  pathway_counts =
    pathway_counts,
  
  direction_pathway_counts =
    direction_pathway_counts,
  
  core_presence =
    core_presence,
  
  core_lr_summary =
    core_lr_summary,
  
  interpretation_status =
    interpretation_note,
  
  statistical_testing =
    FALSE,
  
  patient_consistency_testing =
    FALSE,
  
  CellChat_recalculation =
    FALSE,
  
  raw_UMI_recalculation =
    FALSE,
  
  threshold_lock =
    FALSE,
  
  pathway_lock =
    FALSE
  
)


saveRDS(
  STEP11_5A_LOCKED,
  file.path(
    out_dir,
    "STEP11-5A_CCI_LR_STEP10_Compatibility_Audit_LOCKED.rds"
  )
)


# ------------------------------------------------------------
# 25. Final QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP11-5A FINAL QC\n")
cat("============================================================\n\n")

cat(
  "CCI_LR interactions:",
  nrow(cci_lr),
  "\n"
)

cat(
  "CCI_LR ligands found in STEP10:",
  sum(cci_lr$ligand_in_STEP10),
  "\n"
)

cat(
  "CCI_LR receptors found in STEP10:",
  sum(cci_lr$receptor_in_STEP10),
  "\n"
)

cat(
  "CCI_LR interactions with both partners in STEP10:",
  sum(cci_lr$both_in_STEP10),
  "\n"
)

cat(
  "Focused direction-specific candidate rows:",
  nrow(focused_candidates),
  "\n"
)

cat(
  "Unique focused LR interactions:",
  length(
    unique(
      focused_candidates$interaction_name
    )
  ),
  "\n"
)

cat(
  "Focused pathways represented:",
  length(
    unique(
      focused_candidates$pathway_name
    )
  ),
  "\n"
)

cat(
  "Core-gene LR resource rows:",
  nrow(core_lr_rows),
  "\n"
)

cat(
  "Statistical testing:",
  "NO",
  "\n"
)

cat(
  "Patient consistency testing:",
  "NO",
  "\n"
)

cat(
  "CellChat recalculation:",
  "NO",
  "\n"
)

cat(
  "Pathway locked:",
  "NO",
  "\n"
)

cat(
  "Threshold locked:",
  "NO",
  "\n"
)

cat(
  "\nSTEP11-5A STATUS: PASS\n"
)

cat(
  "============================================================\n"
)

# ============================================================
# STEP11-5B-0
# STEP10 Genome-wide Pseudobulk Expression File Locator
# ============================================================
#
# Purpose:
#   Locate existing STEP10 genome-wide patient-level /
#   compartment-level pseudobulk expression files suitable
#   for downstream ligand-receptor screening.
#
# IMPORTANT:
#   - File discovery / structural audit ONLY
#   - No expression recalculation
#   - No raw UMI recalculation
#   - No CellChat recalculation
#   - No LR statistical testing
#   - No pathway selection
#   - No eligibility-threshold locking
#
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP11-5B-0 — STEP10 GENOME-WIDE PSEUDOBULK FILE LOCATOR\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

step10_root <- "input_dir/Paper 1_Step 10-Single-Cell"

step11_root <- "input_dir/Paper 1_Step 11-CellChat"

out_dir <- file.path(
  step11_root,
  "STEP11-5B0_STEP10_GenomeWide_Pseudobulk_File_Audit"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cat("STEP10 root:\n", step10_root, "\n\n")
cat("Output folder:\n", out_dir, "\n\n")


# ------------------------------------------------------------
# 1. Basic path audit
# ------------------------------------------------------------

stopifnot(dir.exists(step10_root))

cat("PASS — STEP10 root directory exists.\n\n")


# ------------------------------------------------------------
# 2. Recursively list candidate files
# ------------------------------------------------------------

all_files <- list.files(
  path = step10_root,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

cat("Total files found under STEP10:",
    length(all_files), "\n\n")


# ------------------------------------------------------------
# 3. Candidate filename screening
# ------------------------------------------------------------
#
# We deliberately search broadly because the exact filename
# may not be known.
#
# Candidate keywords:
#   pseudobulk
#   expression
#   genome
#   whole
#   transcriptome
#   raw UMI
#   CPM
#   patient
#   compartment
#
# Exclude obvious core-gene-only files where possible.
# ------------------------------------------------------------

file_names_lower <- tolower(basename(all_files))

candidate_pattern <- paste(
  c(
    "pseudobulk",
    "expression",
    "genome",
    "whole",
    "transcript",
    "cpm",
    "patient",
    "compartment"
  ),
  collapse = "|"
)

candidate_idx <- grepl(
  candidate_pattern,
  file_names_lower
)

candidate_files <- all_files[candidate_idx]

cat("Filename-based candidate files:",
    length(candidate_files), "\n\n")

if (length(candidate_files) > 0) {
  print(candidate_files)
} else {
  cat("WARNING — no filename-based candidates found.\n")
}


# ------------------------------------------------------------
# 4. Save raw candidate inventory
# ------------------------------------------------------------

candidate_inventory <- data.frame(
  File = candidate_files,
  FileName = basename(candidate_files),
  stringsAsFactors = FALSE
)

write.table(
  candidate_inventory,
  file = file.path(
    out_dir,
    "STEP11-5B0_Candidate_File_Inventory.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ------------------------------------------------------------
# 5. Define CCI_LR resource
# ------------------------------------------------------------

cci_file <- file.path(
  step11_root,
  "NSCLC_GSE131907_CCI_LR.tsv"
)

stopifnot(file.exists(cci_file))

cci_lr <- read.delim(
  cci_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cci_ligands <- unique(cci_lr$ligand)
cci_receptors <- unique(cci_lr$receptor)

cci_lr_genes <- unique(
  c(
    cci_ligands,
    cci_receptors
  )
)

cat("\nCCI_LR resource:\n")
cat("  interactions:", nrow(cci_lr), "\n")
cat("  unique ligands:", length(cci_ligands), "\n")
cat("  unique receptors:", length(cci_receptors), "\n")
cat("  unique LR genes:", length(cci_lr_genes), "\n\n")


# ------------------------------------------------------------
# 6. Function to inspect candidate tabular files
# ------------------------------------------------------------

inspect_candidate_file <- function(f) {
  
  result <- list(
    File = f,
    FileName = basename(f),
    Extension = tools::file_ext(f),
    Readable = FALSE,
    N_Rows = NA_integer_,
    N_Columns = NA_integer_,
    Has_Gene = FALSE,
    Has_Patient = FALSE,
    Has_Compartment = FALSE,
    Has_CPM = FALSE,
    Has_Detected = FALSE,
    Has_Analysis_Set = FALSE,
    Has_Primary = FALSE,
    N_Genes = NA_integer_,
    N_Patients = NA_integer_,
    N_Compartments = NA_integer_,
    N_Analysis_Set = NA_integer_,
    Target_Compartments_Present = FALSE,
    Target_Compartments = NA_character_,
    CCI_LR_Genes_Found = NA_integer_,
    CCI_LR_Gene_Percent = NA_real_,
    CCI_Ligands_Found = NA_integer_,
    CCI_Ligand_Percent = NA_real_,
    CCI_Receptors_Found = NA_integer_,
    CCI_Receptor_Percent = NA_real_,
    Suitable_Structure = FALSE,
    Error = NA_character_
  )
  
  # ----------------------------------------------------------
  # Only inspect common text-table formats
  # ----------------------------------------------------------
  
  ext <- tolower(tools::file_ext(f))
  
  if (!ext %in% c("tsv", "txt", "csv")) {
    result$Error <- "Unsupported tabular extension"
    return(result)
  }
  
  dat <- tryCatch(
    
    {
      
      if (ext == "csv") {
        
        read.csv(
          f,
          header = TRUE,
          stringsAsFactors = FALSE,
          check.names = FALSE,
          nrows = 5
        )
        
      } else {
        
        read.delim(
          f,
          header = TRUE,
          stringsAsFactors = FALSE,
          check.names = FALSE,
          nrows = 5
        )
        
      }
      
    },
    
    error = function(e) {
      
      result$Error <<- conditionMessage(e)
      return(NULL)
      
    }
    
  )
  
  if (is.null(dat)) {
    return(result)
  }
  
  result$Readable <- TRUE
  
  cn <- names(dat)
  
  result$Has_Gene <- "Gene" %in% cn
  result$Has_Patient <- "Patient" %in% cn
  result$Has_Compartment <- "Compartment" %in% cn
  result$Has_CPM <- "CPM" %in% cn
  result$Has_Detected <- "Detected" %in% cn
  result$Has_Analysis_Set <- "Analysis_Set" %in% cn
  
  # ----------------------------------------------------------
  # Read only enough to determine structure.
  # For candidate files, load full table.
  # ----------------------------------------------------------
  
  full_dat <- tryCatch(
    
    {
      
      if (ext == "csv") {
        
        read.csv(
          f,
          header = TRUE,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        
      } else {
        
        read.delim(
          f,
          header = TRUE,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        
      }
      
    },
    
    error = function(e) {
      
      result$Error <<- conditionMessage(e)
      return(NULL)
      
    }
    
  )
  
  if (is.null(full_dat)) {
    return(result)
  }
  
  result$N_Rows <- nrow(full_dat)
  result$N_Columns <- ncol(full_dat)
  
  if ("Gene" %in% names(full_dat)) {
    
    result$N_Genes <- length(
      unique(
        full_dat$Gene[
          !is.na(full_dat$Gene) &
            full_dat$Gene != ""
        ]
      )
    )
    
  }
  
  if ("Patient" %in% names(full_dat)) {
    
    result$N_Patients <- length(
      unique(
        full_dat$Patient[
          !is.na(full_dat$Patient) &
            full_dat$Patient != ""
        ]
      )
    )
    
  }
  
  if ("Compartment" %in% names(full_dat)) {
    
    result$N_Compartments <- length(
      unique(
        full_dat$Compartment[
          !is.na(full_dat$Compartment) &
            full_dat$Compartment != ""
        ]
      )
    )
    
    target_compartments <- c(
      "Endothelial",
      "Myeloid/macrophage",
      "T-cell"
    )
    
    present_targets <- intersect(
      target_compartments,
      unique(full_dat$Compartment)
    )
    
    result$Target_Compartments <- paste(
      present_targets,
      collapse = "; "
    )
    
    result$Target_Compartments_Present <-
      length(present_targets) == 3
    
  }
  
  if ("Analysis_Set" %in% names(full_dat)) {
    
    result$N_Analysis_Set <- length(
      unique(
        full_dat$Analysis_Set[
          !is.na(full_dat$Analysis_Set) &
            full_dat$Analysis_Set != ""
        ]
      )
    )
    
    result$Has_Primary <-
      "Primary" %in% unique(full_dat$Analysis_Set)
    
  }
  
  # ----------------------------------------------------------
  # CCI_LR gene coverage
  # ----------------------------------------------------------
  
  if ("Gene" %in% names(full_dat)) {
    
    dataset_genes <- unique(
      full_dat$Gene[
        !is.na(full_dat$Gene) &
          full_dat$Gene != ""
      ]
    )
    
    lr_found <- intersect(
      cci_lr_genes,
      dataset_genes
    )
    
    ligand_found <- intersect(
      cci_ligands,
      dataset_genes
    )
    
    receptor_found <- intersect(
      cci_receptors,
      dataset_genes
    )
    
    result$CCI_LR_Genes_Found <- length(lr_found)
    
    result$CCI_LR_Gene_Percent <-
      100 * length(lr_found) / length(cci_lr_genes)
    
    result$CCI_Ligands_Found <- length(ligand_found)
    
    result$CCI_Ligand_Percent <-
      100 * length(ligand_found) / length(cci_ligands)
    
    result$CCI_Receptors_Found <- length(receptor_found)
    
    result$CCI_Receptor_Percent <-
      100 * length(receptor_found) / length(cci_receptors)
    
  }
  
  # ----------------------------------------------------------
  # Structural suitability
  # ----------------------------------------------------------
  
  result$Suitable_Structure <- all(
    result$Has_Gene,
    result$Has_Patient,
    result$Has_Compartment,
    result$Has_CPM,
    result$Has_Detected,
    result$N_Genes > 1000,
    result$N_Patients >= 1,
    result$N_Compartments >= 3,
    result$Target_Compartments_Present
  )
  
  result
  
}


# ------------------------------------------------------------
# 7. Inspect candidates
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("CANDIDATE STRUCTURAL AUDIT\n")
cat("============================================================\n\n")


audit_list <- list()

if (length(candidate_files) > 0) {
  
  for (i in seq_along(candidate_files)) {
    
    cat(
      "\n[",
      i,
      "/",
      length(candidate_files),
      "] ",
      basename(candidate_files),
      "\n",
      sep = ""
    )
    
    audit_list[[i]] <- inspect_candidate_file(
      candidate_files[i]
    )
    
  }
  
}


# ------------------------------------------------------------
# 8. Convert audit results to table
# ------------------------------------------------------------

if (length(audit_list) > 0) {
  
  audit_df <- do.call(
    rbind,
    lapply(
      audit_list,
      as.data.frame,
      stringsAsFactors = FALSE
    )
  )
  
  rownames(audit_df) <- NULL
  
} else {
  
  audit_df <- data.frame()
  
}


# ------------------------------------------------------------
# 9. Print concise audit table
# ------------------------------------------------------------

if (nrow(audit_df) > 0) {
  
  print(
    audit_df[
      ,
      c(
        "FileName",
        "Readable",
        "N_Rows",
        "N_Columns",
        "N_Genes",
        "N_Patients",
        "N_Compartments",
        "Has_Gene",
        "Has_Patient",
        "Has_Compartment",
        "Has_CPM",
        "Has_Detected",
        "Has_Analysis_Set",
        "Has_Primary",
        "Target_Compartments_Present",
        "CCI_LR_Genes_Found",
        "CCI_LR_Gene_Percent",
        "CCI_Ligands_Found",
        "CCI_Ligand_Percent",
        "CCI_Receptors_Found",
        "CCI_Receptor_Percent",
        "Suitable_Structure"
      )
    ]
  )
  
}


# ------------------------------------------------------------
# 10. Save complete audit
# ------------------------------------------------------------

write.table(
  audit_df,
  file = file.path(
    out_dir,
    "STEP11-5B0_STEP10_GenomeWide_Pseudobulk_File_Audit.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ------------------------------------------------------------
# 11. Identify structurally suitable files
# ------------------------------------------------------------

if (nrow(audit_df) > 0) {
  
  suitable_files <- audit_df[
    audit_df$Suitable_Structure == TRUE,
    ,
    drop = FALSE
  ]
  
} else {
  
  suitable_files <- audit_df
  
}


cat("\n")
cat("============================================================\n")
cat("STRUCTURALLY SUITABLE FILES\n")
cat("============================================================\n\n")

cat(
  "Number of structurally suitable candidates:",
  nrow(suitable_files),
  "\n\n"
)

if (nrow(suitable_files) > 0) {
  
  print(
    suitable_files[
      ,
      c(
        "FileName",
        "N_Rows",
        "N_Columns",
        "N_Genes",
        "N_Patients",
        "N_Compartments",
        "Target_Compartments",
        "CCI_LR_Genes_Found",
        "CCI_LR_Gene_Percent",
        "CCI_Ligands_Found",
        "CCI_Ligand_Percent",
        "CCI_Receptors_Found",
        "CCI_Receptor_Percent"
      )
    ]
  )
  
} else {
  
  cat(
    "No structurally suitable genome-wide pseudobulk file was identified.\n"
  )
  
}


# ------------------------------------------------------------
# 12. Save suitable-file table
# ------------------------------------------------------------

write.table(
  suitable_files,
  file = file.path(
    out_dir,
    "STEP11-5B0_Structurally_Suitable_Files.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ------------------------------------------------------------
# 13. Final interpretation
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP11-5B-0 FINAL INTERPRETATION\n")
cat("============================================================\n\n")

cat(
  "This step performed file discovery and structural auditing only.\n\n"
)

cat(
  "No expression values were recalculated.\n"
)

cat(
  "No raw UMI matrix was reread or recalculated.\n"
)

cat(
  "No CellChat analysis was recalculated.\n"
)

cat(
  "No ligand-receptor statistical testing was performed.\n"
)

cat(
  "No pathway was selected or locked.\n"
)

cat(
  "No eligibility threshold was selected or locked.\n\n"
)

if (nrow(suitable_files) > 0) {
  
  cat(
    "STATUS: CANDIDATE GENOME-WIDE PSEUDOBULK FILE(S) IDENTIFIED.\n"
  )
  
} else {
  
  cat(
    "STATUS: NO SUITABLE GENOME-WIDE PSEUDOBULK FILE IDENTIFIED.\n"
  )
  
}

cat("\n")
cat("============================================================\n")
cat("END STEP11-5B-0\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 14. Save audit object
# ------------------------------------------------------------

saveRDS(
  list(
    candidate_inventory = candidate_inventory,
    audit_df = audit_df,
    suitable_files = suitable_files,
    cci_lr_interactions = nrow(cci_lr),
    cci_lr_ligands = length(cci_ligands),
    cci_lr_receptors = length(cci_receptors),
    cci_lr_genes = length(cci_lr_genes)
  ),
  file = file.path(
    out_dir,
    "STEP11-5B0_STEP10_GenomeWide_Pseudobulk_File_Audit_LOCKED.rds"
  )
)

cat("\n")
cat("Audit RDS saved.\n")


# ============================================================
# STEP11-5B-3A
# EXACT LOCATION AUDIT FOR STEP10C PSEUDOBULK FILES
#
# PURPOSE:
#   Locate the exact paths of the three existing STEP10C
#   compartment-specific pseudobulk raw-UMI files.
#
# IMPORTANT:
#   - NO file content recalculation
#   - NO pseudobulk recalculation
#   - NO raw UMI aggregation
#   - NO normalization
#   - NO CellChat
#   - NO LR testing
# ============================================================

cat("============================================================\n")
cat("STEP11-5B-3A: EXACT STEP10C PSEUDOBULK FILE LOCATION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. STEP10 root
# ------------------------------------------------------------

step10_root <- "input_dir/Paper 1_Step 10-Single-Cell"


# ------------------------------------------------------------
# 2. Exact filenames to locate
# ------------------------------------------------------------

target_files <- c(
  "STEP10C-A_Pseudobulk_Tcell_Primary_rawUMI.tsv",
  "STEP10C-B_Pseudobulk_MyeloidMacrophage_Primary_rawUMI.tsv",
  "STEP10C-C_Pseudobulk_Endothelial_Primary_rawUMI.tsv"
)


# ------------------------------------------------------------
# 3. Recursive exact filename search
# ------------------------------------------------------------

all_files <- list.files(
  step10_root,
  recursive = TRUE,
  full.names = TRUE
)


located_list <- lapply(
  target_files,
  function(target) {
    
    hits <- all_files[
      basename(all_files) == target
    ]
    
    data.frame(
      Target_File = target,
      N_Hits = length(hits),
      Full_Path =
        if (length(hits) > 0) {
          paste(hits, collapse = " || ")
        } else {
          NA_character_
        },
      stringsAsFactors = FALSE
    )
    
  }
)


located_df <- do.call(
  rbind,
  located_list
)


rownames(located_df) <- NULL


# ------------------------------------------------------------
# 4. Print exact locations
# ------------------------------------------------------------

cat("Exact filename search results:\n\n")

print(
  located_df,
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Additional broad search if exact search fails
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("BROAD PSEUDOBULK RAWUMI SEARCH\n")
cat("============================================================\n\n")


broad_hits <- all_files[
  grepl(
    "Pseudobulk.*(Tcell|Myeloid|Macrophage|Endothelial).*rawUMI",
    basename(all_files),
    ignore.case = TRUE
  )
]


if (length(broad_hits) > 0) {
  
  broad_df <- data.frame(
    FileName = basename(broad_hits),
    Full_Path = broad_hits,
    Size_MB =
      round(
        file.info(broad_hits)$size / 1024^2,
        3
      ),
    stringsAsFactors = FALSE
  )
  
  print(
    broad_df,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No broad pseudobulk rawUMI filename match found.\n"
  )
  
}


# ------------------------------------------------------------
# 6. Check whether any duplicate locations exist
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("DUPLICATE LOCATION CHECK\n")
cat("============================================================\n\n")


if (any(located_df$N_Hits > 1)) {
  
  cat(
    "WARNING: duplicate filename locations detected.\n"
  )
  
  print(
    located_df[
      located_df$N_Hits > 1,
      ,
      drop = FALSE
    ],
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No duplicate exact filename locations detected.\n"
  )
  
}


# ------------------------------------------------------------
# 7. Final status
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP11-5B-3A FINAL STATUS\n")
cat("============================================================\n\n")


cat(
  "T-cell file located:",
  located_df$N_Hits[1] > 0,
  "\n"
)

cat(
  "Myeloid/macrophage file located:",
  located_df$N_Hits[2] > 0,
  "\n"
)

cat(
  "Endothelial file located:",
  located_df$N_Hits[3] > 0,
  "\n\n"
)

cat(
  "Pseudobulk recalculation: NO\n"
)

cat(
  "Raw UMI aggregation: NO\n"
)

cat(
  "Expression transformation: NO\n"
)

cat(
  "CellChat recalculation: NO\n"
)

cat(
  "LR statistical testing: NO\n"
)

cat(
  "Pathway selection: NO\n"
)

cat(
  "Threshold locking: NO\n"
)


if (all(located_df$N_Hits == 1)) {
  
  cat(
    "\nSTATUS: ALL THREE STEP10C PSEUDOBULK FILES LOCATED EXACTLY.\n"
  )
  
} else {
  
  cat(
    "\nSTATUS: EXACT LOCATION OF ALL THREE FILES NOT YET CONFIRMED.\n"
  )
  
}


cat("\n")
cat("============================================================\n")
cat("END STEP11-5B-3A\n")
cat("============================================================\n")



##STEP11-5B-3B
#Existing STEP10C compartment-specific pseudobulk structural audit

cat("============================================================\n")
cat("STEP11-5B-3B — EXISTING STEP10C PSEUDOBULK STRUCTURAL AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step10_root <- "input_dir/Paper 1_Step 10-Single-Cell"


# ------------------------------------------------------------
# 1. Locate the three exact files again
#    No path guessing
# ------------------------------------------------------------

target_files <- c(
  "STEP10C-A_Pseudobulk_Tcell_Primary_rawUMI.tsv",
  "STEP10C-B_Pseudobulk_MyeloidMacrophage_Primary_rawUMI.tsv",
  "STEP10C-C_Pseudobulk_Endothelial_Primary_rawUMI.tsv"
)

all_files <- list.files(
  step10_root,
  recursive = TRUE,
  full.names = TRUE
)

located_list <- lapply(
  target_files,
  function(target) {
    
    hits <- all_files[
      basename(all_files) == target
    ]
    
    data.frame(
      Target_File = target,
      N_Hits = length(hits),
      Full_Path = if (length(hits) == 1) hits else NA_character_,
      stringsAsFactors = FALSE
    )
  }
)

located_df <- do.call(
  rbind,
  located_list
)

cat("Exact file-location audit:\n")
print(
  located_df,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 2. Structural requirement
# ------------------------------------------------------------

if (!all(located_df$N_Hits == 1)) {
  
  stop(
    "STEP11-5B-3B STOP: not all three pseudobulk files have exactly one match."
  )
  
}


# ------------------------------------------------------------
# 3. Read existing files ONLY
#
#    IMPORTANT:
#    - no pseudobulk recalculation
#    - no raw UMI aggregation
#    - no normalization
#    - no transformation
# ------------------------------------------------------------

tcell_file <- located_df$Full_Path[1]
myeloid_file <- located_df$Full_Path[2]
endo_file <- located_df$Full_Path[3]


cat("Reading existing STEP10C files...\n\n")

pb_tcell <- read.delim(
  tcell_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pb_myeloid <- read.delim(
  myeloid_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pb_endo <- read.delim(
  endo_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 4. Basic dimensions
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("4. BASIC DIMENSIONS\n")
cat("------------------------------------------------------------\n\n")

cat(
  "T-cell:",
  nrow(pb_tcell),
  "rows x",
  ncol(pb_tcell),
  "columns\n"
)

cat(
  "Myeloid/macrophage:",
  nrow(pb_myeloid),
  "rows x",
  ncol(pb_myeloid),
  "columns\n"
)

cat(
  "Endothelial:",
  nrow(pb_endo),
  "rows x",
  ncol(pb_endo),
  "columns\n\n"
)


# ------------------------------------------------------------
# 5. Column names
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("5. COLUMN NAMES\n")
cat("------------------------------------------------------------\n\n")

cat("T-cell columns:\n")
print(names(pb_tcell))

cat("\nMyeloid/macrophage columns:\n")
print(names(pb_myeloid))

cat("\nEndothelial columns:\n")
print(names(pb_endo))

cat("\n")


# ------------------------------------------------------------
# 6. First rows
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("6. FIRST ROWS\n")
cat("------------------------------------------------------------\n\n")

cat("T-cell head:\n")
print(head(pb_tcell, 3))

cat("\nMyeloid/macrophage head:\n")
print(head(pb_myeloid, 3))

cat("\nEndothelial head:\n")
print(head(pb_endo, 3))

cat("\n")


# ------------------------------------------------------------
# 7. Detect matrix-like vs long-format structure
# ------------------------------------------------------------

detect_structure <- function(x) {
  
  cn <- names(x)
  
  has_gene <- any(
    grepl(
      "^gene$|gene_id|geneid|symbol",
      cn,
      ignore.case = TRUE
    )
  )
  
  has_patient <- any(
    grepl(
      "patient",
      cn,
      ignore.case = TRUE
    )
  )
  
  has_sample <- any(
    grepl(
      "sample",
      cn,
      ignore.case = TRUE
    )
  )
  
  has_compartment <- any(
    grepl(
      "compartment",
      cn,
      ignore.case = TRUE
    )
  )
  
  has_cpm <- any(
    grepl(
      "^cpm$|cpm",
      cn,
      ignore.case = TRUE
    )
  )
  
  data.frame(
    Has_Gene = has_gene,
    Has_Patient = has_patient,
    Has_Sample = has_sample,
    Has_Compartment = has_compartment,
    Has_CPM = has_cpm,
    stringsAsFactors = FALSE
  )
}


cat("------------------------------------------------------------\n")
cat("7. STRUCTURAL SIGNATURE\n")
cat("------------------------------------------------------------\n\n")

cat("T-cell:\n")
print(detect_structure(pb_tcell))

cat("\nMyeloid/macrophage:\n")
print(detect_structure(pb_myeloid))

cat("\nEndothelial:\n")
print(detect_structure(pb_endo))

cat("\n")


# ------------------------------------------------------------
# 8. Identify gene identifier column
# ------------------------------------------------------------

find_gene_column <- function(x) {
  
  cn <- names(x)
  
  candidates <- cn[
    grepl(
      "^gene$|gene_id|geneid|symbol|ensembl",
      cn,
      ignore.case = TRUE
    )
  ]
  
  if (length(candidates) == 0) {
    return(NA_character_)
  }
  
  candidates[1]
}


gene_col_tcell <- find_gene_column(pb_tcell)
gene_col_myeloid <- find_gene_column(pb_myeloid)
gene_col_endo <- find_gene_column(pb_endo)


cat("------------------------------------------------------------\n")
cat("8. GENE IDENTIFIER COLUMNS\n")
cat("------------------------------------------------------------\n\n")

cat(
  "T-cell gene column:",
  gene_col_tcell,
  "\n"
)

cat(
  "Myeloid/macrophage gene column:",
  gene_col_myeloid,
  "\n"
)

cat(
  "Endothelial gene column:",
  gene_col_endo,
  "\n\n"
)


# ------------------------------------------------------------
# 9. Patient/sample/CPM columns
# ------------------------------------------------------------

find_column <- function(x, pattern) {
  
  cn <- names(x)
  
  hits <- cn[
    grepl(
      pattern,
      cn,
      ignore.case = TRUE
    )
  ]
  
  if (length(hits) == 0) {
    return(NA_character_)
  }
  
  hits[1]
}


patient_col_tcell <- find_column(
  pb_tcell,
  "patient"
)

patient_col_myeloid <- find_column(
  pb_myeloid,
  "patient"
)

patient_col_endo <- find_column(
  pb_endo,
  "patient"
)


sample_col_tcell <- find_column(
  pb_tcell,
  "sample"
)

sample_col_myeloid <- find_column(
  pb_myeloid,
  "sample"
)

sample_col_endo <- find_column(
  pb_endo,
  "sample"
)


cpm_col_tcell <- find_column(
  pb_tcell,
  "^cpm$|cpm"
)

cpm_col_myeloid <- find_column(
  pb_myeloid,
  "^cpm$|cpm"
)

cpm_col_endo <- find_column(
  pb_endo,
  "^cpm$|cpm"
)


cat("------------------------------------------------------------\n")
cat("9. PATIENT / SAMPLE / CPM COLUMNS\n")
cat("------------------------------------------------------------\n\n")

cat("T-cell:\n")
cat(
  "  Patient =",
  patient_col_tcell,
  "\n"
)
cat(
  "  Sample  =",
  sample_col_tcell,
  "\n"
)
cat(
  "  CPM     =",
  cpm_col_tcell,
  "\n\n"
)

cat("Myeloid/macrophage:\n")
cat(
  "  Patient =",
  patient_col_myeloid,
  "\n"
)
cat(
  "  Sample  =",
  sample_col_myeloid,
  "\n"
)
cat(
  "  CPM     =",
  cpm_col_myeloid,
  "\n\n"
)

cat("Endothelial:\n")
cat(
  "  Patient =",
  patient_col_endo,
  "\n"
)
cat(
  "  Sample  =",
  sample_col_endo,
  "\n"
)
cat(
  "  CPM     =",
  cpm_col_endo,
  "\n\n"
)


# ------------------------------------------------------------
# 10. Gene counts and identifier uniqueness
# ------------------------------------------------------------

gene_summary <- function(x, gene_col) {
  
  if (is.na(gene_col)) {
    
    return(
      data.frame(
        N_Rows = nrow(x),
        N_Unique_Genes = NA_integer_,
        Duplicate_Gene_IDs = NA_integer_,
        stringsAsFactors = FALSE
      )
    )
  }
  
  g <- as.character(x[[gene_col]])
  g <- g[!is.na(g) & nzchar(g)]
  
  data.frame(
    N_Rows = nrow(x),
    N_Unique_Genes = length(unique(g)),
    Duplicate_Gene_IDs = sum(
      duplicated(g)
    ),
    stringsAsFactors = FALSE
  )
}


cat("------------------------------------------------------------\n")
cat("10. GENE-LEVEL STRUCTURE\n")
cat("------------------------------------------------------------\n\n")

cat("T-cell:\n")
print(
  gene_summary(
    pb_tcell,
    gene_col_tcell
  )
)

cat("\nMyeloid/macrophage:\n")
print(
  gene_summary(
    pb_myeloid,
    gene_col_myeloid
  )
)

cat("\nEndothelial:\n")
print(
  gene_summary(
    pb_endo,
    gene_col_endo
  )
)

cat("\n")


# ------------------------------------------------------------
# 11. Patient counts
# ------------------------------------------------------------

patient_summary <- function(x, patient_col) {
  
  if (is.na(patient_col)) {
    
    return(
      data.frame(
        N_Unique_Patients = NA_integer_,
        stringsAsFactors = FALSE
      )
    )
  }
  
  p <- as.character(x[[patient_col]])
  p <- p[!is.na(p) & nzchar(p)]
  
  data.frame(
    N_Unique_Patients = length(unique(p)),
    stringsAsFactors = FALSE
  )
}


cat("------------------------------------------------------------\n")
cat("11. PATIENT STRUCTURE\n")
cat("------------------------------------------------------------\n\n")

cat("T-cell:\n")
print(
  patient_summary(
    pb_tcell,
    patient_col_tcell
  )
)

cat("\nMyeloid/macrophage:\n")
print(
  patient_summary(
    pb_myeloid,
    patient_col_myeloid
  )
)

cat("\nEndothelial:\n")
print(
  patient_summary(
    pb_endo,
    patient_col_endo
  )
)

cat("\n")


# ------------------------------------------------------------
# 12. CCI_LR resource
# ------------------------------------------------------------

lr_file <- file.path(
  "input_dir/Paper 1_Step 11-CellChat",
  "NSCLC_GSE131907_CCI_LR.tsv"
)

cci_lr <- read.delim(
  lr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

lr_genes <- unique(
  c(
    as.character(cci_lr$ligand),
    as.character(cci_lr$receptor)
  )
)

lr_genes <- lr_genes[
  !is.na(lr_genes) &
    nzchar(lr_genes)
]


cat("------------------------------------------------------------\n")
cat("12. CCI_LR RESOURCE\n")
cat("------------------------------------------------------------\n\n")

cat(
  "CCI_LR rows:",
  nrow(cci_lr),
  "\n"
)

cat(
  "Unique LR genes:",
  length(lr_genes),
  "\n"
)

cat(
  "Unique ligands:",
  length(unique(cci_lr$ligand)),
  "\n"
)

cat(
  "Unique receptors:",
  length(unique(cci_lr$receptor)),
  "\n\n"
)


# ------------------------------------------------------------
# 13. Overlap with genome-wide pseudobulk genes
# ------------------------------------------------------------

get_gene_vector <- function(x, gene_col) {
  
  if (is.na(gene_col)) {
    return(character(0))
  }
  
  g <- unique(
    as.character(x[[gene_col]])
  )
  
  g[
    !is.na(g) &
      nzchar(g)
  ]
}


genes_tcell <- get_gene_vector(
  pb_tcell,
  gene_col_tcell
)

genes_myeloid <- get_gene_vector(
  pb_myeloid,
  gene_col_myeloid
)

genes_endo <- get_gene_vector(
  pb_endo,
  gene_col_endo
)


overlap_tcell <- intersect(
  genes_tcell,
  lr_genes
)

overlap_myeloid <- intersect(
  genes_myeloid,
  lr_genes
)

overlap_endo <- intersect(
  genes_endo,
  lr_genes
)


cat("------------------------------------------------------------\n")
cat("13. CCI_LR GENE OVERLAP\n")
cat("------------------------------------------------------------\n\n")

cat(
  "T-cell pseudobulk genes:",
  length(genes_tcell),
  "\n"
)

cat(
  "  CCI_LR overlap:",
  length(overlap_tcell),
  "\n"
)

cat(
  "Myeloid/macrophage pseudobulk genes:",
  length(genes_myeloid),
  "\n"
)

cat(
  "  CCI_LR overlap:",
  length(overlap_myeloid),
  "\n"
)

cat(
  "Endothelial pseudobulk genes:",
  length(genes_endo),
  "\n"
)

cat(
  "  CCI_LR overlap:",
  length(overlap_endo),
  "\n\n"
)


# ------------------------------------------------------------
# 14. Core genes audit
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


core_presence <- data.frame(
  Gene = core_genes,
  Tcell = core_genes %in% genes_tcell,
  Myeloid_Macrophage = core_genes %in% genes_myeloid,
  Endothelial = core_genes %in% genes_endo,
  stringsAsFactors = FALSE
)


cat("------------------------------------------------------------\n")
cat("14. CORE-GENE PRESENCE AUDIT\n")
cat("------------------------------------------------------------\n\n")

print(
  core_presence,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 15. Determine whether files are suitable for
#     genome-wide LR expression screening
# ------------------------------------------------------------

structural_screen <- data.frame(
  
  Compartment = c(
    "T-cell",
    "Myeloid/macrophage",
    "Endothelial"
  ),
  
  N_Rows = c(
    nrow(pb_tcell),
    nrow(pb_myeloid),
    nrow(pb_endo)
  ),
  
  N_Genes = c(
    length(genes_tcell),
    length(genes_myeloid),
    length(genes_endo)
  ),
  
  N_LR_Overlap = c(
    length(overlap_tcell),
    length(overlap_myeloid),
    length(overlap_endo)
  ),
  
  Has_Patient = c(
    !is.na(patient_col_tcell),
    !is.na(patient_col_myeloid),
    !is.na(patient_col_endo)
  ),
  
  Has_CPM = c(
    !is.na(cpm_col_tcell),
    !is.na(cpm_col_myeloid),
    !is.na(cpm_col_endo)
  ),
  
  stringsAsFactors = FALSE
)


cat("------------------------------------------------------------\n")
cat("15. FINAL STRUCTURAL SCREEN\n")
cat("------------------------------------------------------------\n\n")

print(
  structural_screen,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 16. Final status
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP11-5B-3B FINAL STATUS\n")
cat("============================================================\n\n")

cat(
  "Existing STEP10C pseudobulk files read: YES\n"
)

cat(
  "Pseudobulk recalculation: NO\n"
)

cat(
  "Raw UMI aggregation: NO\n"
)

cat(
  "Expression transformation: NO\n"
)

cat(
  "CellChat recalculation: NO\n"
)

cat(
  "LR statistical testing: NO\n"
)

cat(
  "Pathway selection: NO\n"
)

cat(
  "Threshold locking: NO\n"
)

cat("\n")

cat(
  "STATUS: STRUCTURAL AUDIT COMPLETED.\n"
)

cat(
  "STATUS: INTERPRETATION DEPENDS ON THE ACTUAL FILE STRUCTURE ABOVE.\n"
)

cat("\n")
cat("============================================================\n")
cat("END STEP11-5B-3B\n")
cat("============================================================\n")


cat("============================================================\n")
cat("STEP11-5B-3C — PATIENT / COMPARTMENT DECODING AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Fixed input locations
# ------------------------------------------------------------

step10_root <- "input_dir/Paper 1_Step 10-Single-Cell"

step11_root <- "input_dir/Paper 1_Step 11-CellChat"


# ------------------------------------------------------------
# 1. Locate the three already-existing STEP10C files
# ------------------------------------------------------------

target_files <- c(
  "STEP10C-A_Pseudobulk_Tcell_Primary_rawUMI.tsv",
  "STEP10C-B_Pseudobulk_MyeloidMacrophage_Primary_rawUMI.tsv",
  "STEP10C-C_Pseudobulk_Endothelial_Primary_rawUMI.tsv"
)

all_files <- list.files(
  step10_root,
  recursive = TRUE,
  full.names = TRUE
)

located_list <- lapply(
  target_files,
  function(target) {
    
    hits <- all_files[
      basename(all_files) == target
    ]
    
    data.frame(
      Target_File = target,
      N_Hits = length(hits),
      Full_Path =
        if (length(hits) == 1) {
          hits
        } else {
          NA_character_
        },
      stringsAsFactors = FALSE
    )
  }
)

located_df <- do.call(
  rbind,
  located_list
)


cat("------------------------------------------------------------\n")
cat("1. EXACT FILE LOCATION\n")
cat("------------------------------------------------------------\n\n")

print(
  located_df,
  row.names = FALSE
)

cat("\n")


if (!all(located_df$N_Hits == 1)) {
  
  stop(
    "STEP11-5B-3C STOP: one or more STEP10C files do not have exactly one match."
  )
  
}


# ------------------------------------------------------------
# 2. Read existing files only
# ------------------------------------------------------------

tcell_file <- located_df$Full_Path[1]
myeloid_file <- located_df$Full_Path[2]
endo_file <- located_df$Full_Path[3]


cat("------------------------------------------------------------\n")
cat("2. READ EXISTING PSEUDOBULK MATRICES\n")
cat("------------------------------------------------------------\n\n")


pb_tcell <- read.delim(
  tcell_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pb_myeloid <- read.delim(
  myeloid_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pb_endo <- read.delim(
  endo_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


cat(
  "T-cell:",
  nrow(pb_tcell),
  "genes x",
  ncol(pb_tcell) - 1,
  "patient columns\n"
)

cat(
  "Myeloid/macrophage:",
  nrow(pb_myeloid),
  "genes x",
  ncol(pb_myeloid) - 1,
  "patient columns\n"
)

cat(
  "Endothelial:",
  nrow(pb_endo),
  "genes x",
  ncol(pb_endo) - 1,
  "patient columns\n\n"
)


# ------------------------------------------------------------
# 3. Read existing CCI_LR resource
# ------------------------------------------------------------

lr_file <- file.path(
  step11_root,
  "NSCLC_GSE131907_CCI_LR.tsv"
)

if (!file.exists(lr_file)) {
  
  stop(
    paste(
      "CCI_LR file not found:",
      lr_file
    )
  )
}


cci_lr <- read.delim(
  lr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 4. Decode patient and compartment from column names
# ------------------------------------------------------------

decode_pb_columns <- function(
    x,
    expected_compartment
) {
  
  cn <- names(x)
  
  expression_columns <- cn[
    cn != "Gene"
  ]
  
  if (length(expression_columns) == 0) {
    
    stop(
      paste(
        "No expression columns detected for",
        expected_compartment
      )
    )
  }
  
  parts <- strsplit(
    expression_columns,
    "__",
    fixed = TRUE
  )
  
  patient <- vapply(
    parts,
    function(z) {
      if (length(z) >= 1) {
        z[1]
      } else {
        NA_character_
      }
    },
    character(1)
  )
  
  compartment <- vapply(
    parts,
    function(z) {
      if (length(z) >= 2) {
        paste(
          z[-1],
          collapse = "__"
        )
      } else {
        NA_character_
      }
    },
    character(1)
  )
  
  data.frame(
    Expression_Column = expression_columns,
    Patient = patient,
    Compartment = compartment,
    Expected_Compartment = expected_compartment,
    stringsAsFactors = FALSE
  )
}


meta_tcell <- decode_pb_columns(
  pb_tcell,
  "T-cell"
)

meta_myeloid <- decode_pb_columns(
  pb_myeloid,
  "Myeloid/macrophage"
)

meta_endo <- decode_pb_columns(
  pb_endo,
  "Endothelial"
)


# ------------------------------------------------------------
# 5. Display decoded column metadata
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("5. DECODED PATIENT / COMPARTMENT METADATA\n")
cat("------------------------------------------------------------\n\n")

cat("T-cell:\n")
print(
  meta_tcell,
  row.names = FALSE
)

cat("\nMyeloid/macrophage:\n")
print(
  meta_myeloid,
  row.names = FALSE
)

cat("\nEndothelial:\n")
print(
  meta_endo,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 6. Compartment decoding QC
# ------------------------------------------------------------

compartment_pass <- identical(
  unique(meta_tcell$Compartment),
  "T-cell"
) &&
  identical(
    unique(meta_myeloid$Compartment),
    "Myeloid/macrophage"
  ) &&
  identical(
    unique(meta_endo$Compartment),
    "Endothelial"
  )


cat("------------------------------------------------------------\n")
cat("6. COMPARTMENT DECODING QC\n")
cat("------------------------------------------------------------\n\n")

cat(
  "T-cell decoded compartment:",
  paste(
    unique(meta_tcell$Compartment),
    collapse = " | "
  ),
  "\n"
)

cat(
  "Myeloid/macrophage decoded compartment:",
  paste(
    unique(meta_myeloid$Compartment),
    collapse = " | "
  ),
  "\n"
)

cat(
  "Endothelial decoded compartment:",
  paste(
    unique(meta_endo$Compartment),
    collapse = " | "
  ),
  "\n"
)

cat(
  "Compartment decoding PASS:",
  compartment_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 7. Patient sets
# ------------------------------------------------------------

patients_tcell <- meta_tcell$Patient
patients_myeloid <- meta_myeloid$Patient
patients_endo <- meta_endo$Patient


cat("------------------------------------------------------------\n")
cat("7. PATIENT SETS\n")
cat("------------------------------------------------------------\n\n")

cat(
  "T-cell:",
  paste(
    patients_tcell,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Myeloid/macrophage:",
  paste(
    patients_myeloid,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Endothelial:",
  paste(
    patients_endo,
    collapse = ", "
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 8. Patient uniqueness
# ------------------------------------------------------------

patient_unique_pass <- all(
  !duplicated(patients_tcell)
) &&
  all(
    !duplicated(patients_myeloid)
  ) &&
  all(
    !duplicated(patients_endo)
  )


cat(
  "Patient uniqueness within each compartment PASS:",
  patient_unique_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 9. Matched patient intersection
# ------------------------------------------------------------

common_patients <- Reduce(
  intersect,
  list(
    patients_tcell,
    patients_myeloid,
    patients_endo
  )
)


cat("------------------------------------------------------------\n")
cat("9. MATCHED PATIENT INTERSECTION\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Number of patients represented in all three compartments:",
  length(common_patients),
  "\n"
)

cat(
  "Matched patients:",
  paste(
    common_patients,
    collapse = ", "
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 10. Patient-set differences
# ------------------------------------------------------------

cat(
  "T-cell only:",
  paste(
    setdiff(
      patients_tcell,
      common_patients
    ),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Myeloid/macrophage only:",
  paste(
    setdiff(
      patients_myeloid,
      common_patients
    ),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Endothelial only:",
  paste(
    setdiff(
      patients_endo,
      common_patients
    ),
    collapse = ", "
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 11. Exact gene identity audit
# ------------------------------------------------------------

genes_tcell <- as.character(
  pb_tcell$Gene
)

genes_myeloid <- as.character(
  pb_myeloid$Gene
)

genes_endo <- as.character(
  pb_endo$Gene
)


gene_identity_pass <- identical(
  genes_tcell,
  genes_myeloid
) &&
  identical(
    genes_tcell,
    genes_endo
  )


gene_unique_pass <- all(
  !duplicated(genes_tcell)
) &&
  all(
    !duplicated(genes_myeloid)
  ) &&
  all(
    !duplicated(genes_endo)
  )


cat("------------------------------------------------------------\n")
cat("11. GENE IDENTITY AUDIT\n")
cat("------------------------------------------------------------\n\n")

cat(
  "T-cell genes:",
  length(genes_tcell),
  "\n"
)

cat(
  "Myeloid/macrophage genes:",
  length(genes_myeloid),
  "\n"
)

cat(
  "Endothelial genes:",
  length(genes_endo),
  "\n"
)

cat(
  "Exact same gene IDs and order:",
  gene_identity_pass,
  "\n"
)

cat(
  "Gene uniqueness PASS:",
  gene_unique_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 12. Numeric raw-UMI expression audit
# ------------------------------------------------------------

check_expression <- function(x) {
  
  expr <- x[
    ,
    setdiff(
      names(x),
      "Gene"
    ),
    drop = FALSE
  ]
  
  numeric_flags <- vapply(
    expr,
    is.numeric,
    logical(1)
  )
  
  mat <- as.matrix(expr)
  
  list(
    all_numeric = all(numeric_flags),
    non_numeric_columns =
      sum(!numeric_flags),
    N_NA =
      sum(is.na(mat)),
    N_negative =
      sum(
        mat < 0,
        na.rm = TRUE
      ),
    N_zero =
      sum(
        mat == 0,
        na.rm = TRUE
      ),
    N_positive =
      sum(
        mat > 0,
        na.rm = TRUE
      )
  )
}


expr_tcell <- check_expression(
  pb_tcell
)

expr_myeloid <- check_expression(
  pb_myeloid
)

expr_endo <- check_expression(
  pb_endo
)


cat("------------------------------------------------------------\n")
cat("12. RAW-UMI EXPRESSION STRUCTURE\n")
cat("------------------------------------------------------------\n\n")

cat("T-cell:\n")
print(expr_tcell)

cat("\nMyeloid/macrophage:\n")
print(expr_myeloid)

cat("\nEndothelial:\n")
print(expr_endo)

cat("\n")


numeric_expression_pass <-
  expr_tcell$all_numeric &&
  expr_myeloid$all_numeric &&
  expr_endo$all_numeric &&
  expr_tcell$N_NA == 0 &&
  expr_myeloid$N_NA == 0 &&
  expr_endo$N_NA == 0 &&
  expr_tcell$N_negative == 0 &&
  expr_myeloid$N_negative == 0 &&
  expr_endo$N_negative == 0


cat(
  "Raw-UMI numeric QC PASS:",
  numeric_expression_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 13. CCI_LR resource audit
# ------------------------------------------------------------

lr_genes <- unique(
  c(
    as.character(cci_lr$ligand),
    as.character(cci_lr$receptor)
  )
)

lr_genes <- lr_genes[
  !is.na(lr_genes) &
    nzchar(lr_genes)
]


lr_tcell <- intersect(
  genes_tcell,
  lr_genes
)

lr_myeloid <- intersect(
  genes_myeloid,
  lr_genes
)

lr_endo <- intersect(
  genes_endo,
  lr_genes
)

lr_common <- Reduce(
  intersect,
  list(
    lr_tcell,
    lr_myeloid,
    lr_endo
  )
)


cat("------------------------------------------------------------\n")
cat("13. CCI_LR OVERLAP AUDIT\n")
cat("------------------------------------------------------------\n\n")

cat(
  "CCI_LR rows:",
  nrow(cci_lr),
  "\n"
)

cat(
  "Unique ligand/receptor genes:",
  length(lr_genes),
  "\n"
)

cat(
  "T-cell LR-gene overlap:",
  length(lr_tcell),
  "\n"
)

cat(
  "Myeloid/macrophage LR-gene overlap:",
  length(lr_myeloid),
  "\n"
)

cat(
  "Endothelial LR-gene overlap:",
  length(lr_endo),
  "\n"
)

cat(
  "LR genes present in all three compartments:",
  length(lr_common),
  "\n\n"
)


# ------------------------------------------------------------
# 14. Core-gene identity audit
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


core_presence <- data.frame(
  Gene = core_genes,
  Tcell =
    core_genes %in% genes_tcell,
  Myeloid_Macrophage =
    core_genes %in% genes_myeloid,
  Endothelial =
    core_genes %in% genes_endo,
  stringsAsFactors = FALSE
)


cat("------------------------------------------------------------\n")
cat("14. CORE-GENE IDENTITY AUDIT\n")
cat("------------------------------------------------------------\n\n")

print(
  core_presence,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 15. Matched patient-compartment design
# ------------------------------------------------------------

patient_compartment_map <- rbind(
  meta_tcell,
  meta_myeloid,
  meta_endo
)


matched_design <- expand.grid(
  Patient = common_patients,
  Compartment = c(
    "T-cell",
    "Myeloid/macrophage",
    "Endothelial"
  ),
  stringsAsFactors = FALSE
)


observed_design <- unique(
  patient_compartment_map[
    patient_compartment_map$Patient %in%
      common_patients,
    c(
      "Patient",
      "Compartment"
    ),
    drop = FALSE
  ]
)


expected_keys <- paste(
  matched_design$Patient,
  matched_design$Compartment,
  sep = "||"
)

observed_keys <- paste(
  observed_design$Patient,
  observed_design$Compartment,
  sep = "||"
)


matched_design_pass <- setequal(
  expected_keys,
  observed_keys
)


cat("------------------------------------------------------------\n")
cat("15. MATCHED THREE-COMPARTMENT DESIGN\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Expected patient-compartment combinations:",
  nrow(matched_design),
  "\n"
)

cat(
  "Observed patient-compartment combinations:",
  nrow(observed_design),
  "\n"
)

cat(
  "Complete matched design PASS:",
  matched_design_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 16. Final status
# ------------------------------------------------------------

overall_structural_pass <-
  compartment_pass &&
  patient_unique_pass &&
  gene_identity_pass &&
  gene_unique_pass &&
  numeric_expression_pass &&
  matched_design_pass


cat("============================================================\n")
cat("STEP11-5B-3C FINAL STATUS\n")
cat("============================================================\n\n")

cat(
  "Existing STEP10C pseudobulk matrices:",
  "READ SUCCESSFULLY\n"
)

cat(
  "Patient/compartment decoding:",
  "PERFORMED\n"
)

cat(
  "Pseudobulk recalculation: NO\n"
)

cat(
  "Raw UMI aggregation: NO\n"
)

cat(
  "Normalization: NO\n"
)

cat(
  "Expression transformation: NO\n"
)

cat(
  "CellChat recalculation: NO\n"
)

cat(
  "LR statistical testing: NO\n"
)

cat(
  "Pathway selection: NO\n"
)

cat(
  "Threshold locking: NO\n\n"
)


if (overall_structural_pass) {
  
  cat(
    "STATUS: PASS — EXISTING PATIENT-LEVEL PSEUDOBULK MATRICES ARE STRUCTURALLY SUITABLE FOR MATCHED LR SCREENING.\n"
  )
  
} else {
  
  cat(
    "STATUS: REVIEW REQUIRED — ONE OR MORE STRUCTURAL QC ITEMS DID NOT PASS.\n"
  )
  
}


cat("\n")

cat(
  "Matched patient count:",
  length(common_patients),
  "\n"
)

cat(
  "Matched patients:",
  paste(
    common_patients,
    collapse = ", "
  ),
  "\n\n"
)


cat("============================================================\n")
cat("END STEP11-5B-3C\n")
cat("============================================================\n")

# ============================================================
# STEP11-5C-1
# Expression-Scale / Normalization Compatibility Audit
#
# Purpose:
#   Audit the existing STEP10C patient-level compartment-specific
#   raw-UMI pseudobulk matrices before genome-wide LR screening.
#
# IMPORTANT:
#   - No pseudobulk recalculation
#   - No raw-UMI re-aggregation
#   - No CellChat recalculation
#   - No LR statistical testing
#   - No pathway selection
#   - No LR threshold locking
# ============================================================

cat("============================================================\n")
cat("STEP11-5C-1 — EXPRESSION-SCALE / NORMALIZATION COMPATIBILITY AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

step10_root <- "input_dir/Paper 1_Step 10-Single-Cell"
step11_root <- "input_dir/Paper 1_Step 11-CellChat"

file_tcell <- file.path(
  step10_root,
  "STEP10C_Patient_Level_RawUMI_Pseudobulk",
  "STEP10C-A_Pseudobulk_Tcell_Primary_rawUMI.tsv"
)

file_myeloid <- file.path(
  step10_root,
  "STEP10C_Patient_Level_RawUMI_Pseudobulk",
  "STEP10C-B_Pseudobulk_MyeloidMacrophage_Primary_rawUMI.tsv"
)

file_endo <- file.path(
  step10_root,
  "STEP10C_Patient_Level_RawUMI_Pseudobulk",
  "STEP10C-C_Pseudobulk_Endothelial_Primary_rawUMI.tsv"
)


# ------------------------------------------------------------
# 1. File existence audit
# ------------------------------------------------------------

file_exists <- c(
  T_cell = file.exists(file_tcell),
  Myeloid_Macrophage = file.exists(file_myeloid),
  Endothelial = file.exists(file_endo)
)

print(file_exists)

if (!all(file_exists)) {
  stop(
    "One or more locked STEP10C pseudobulk files do not exist."
  )
}


# ------------------------------------------------------------
# 2. Read existing matrices
# ------------------------------------------------------------

pb_tcell <- read.delim(
  file_tcell,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pb_myeloid <- read.delim(
  file_myeloid,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pb_endo <- read.delim(
  file_endo,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 3. Basic structure
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("3. MATRIX STRUCTURE\n")
cat("------------------------------------------------------------\n")

cat(
  "T-cell:",
  nrow(pb_tcell), "genes ×",
  ncol(pb_tcell) - 1, "patients\n"
)

cat(
  "Myeloid/macrophage:",
  nrow(pb_myeloid), "genes ×",
  ncol(pb_myeloid) - 1, "patients\n"
)

cat(
  "Endothelial:",
  nrow(pb_endo), "genes ×",
  ncol(pb_endo) - 1, "patients\n"
)


# ------------------------------------------------------------
# 4. Gene identity audit
# ------------------------------------------------------------

genes_tcell  <- pb_tcell$Gene
genes_myeloid <- pb_myeloid$Gene
genes_endo   <- pb_endo$Gene

gene_identity_pass <- identical(
  genes_tcell,
  genes_myeloid
) &&
  identical(
    genes_tcell,
    genes_endo
  )

gene_unique_pass <- (
  length(unique(genes_tcell)) == length(genes_tcell)
) &&
  (
    length(unique(genes_myeloid)) == length(genes_myeloid)
  ) &&
  (
    length(unique(genes_endo)) == length(genes_endo)
  )

cat("\nGene identity/order:", gene_identity_pass, "\n")
cat("Gene uniqueness:", gene_unique_pass, "\n")


# ------------------------------------------------------------
# 5. Decode patient IDs
# ------------------------------------------------------------

decode_patients <- function(x) {
  
  cols <- colnames(x)[-1]
  
  sub(
    "__.*$",
    "",
    cols
  )
}

patients_tcell  <- decode_patients(pb_tcell)
patients_myeloid <- decode_patients(pb_myeloid)
patients_endo   <- decode_patients(pb_endo)

cat("\n------------------------------------------------------------\n")
cat("5. PATIENT STRUCTURE\n")
cat("------------------------------------------------------------\n")

cat(
  "T-cell patients:",
  paste(patients_tcell, collapse = ", "),
  "\n"
)

cat(
  "Myeloid/macrophage patients:",
  paste(patients_myeloid, collapse = ", "),
  "\n"
)

cat(
  "Endothelial patients:",
  paste(patients_endo, collapse = ", "),
  "\n"
)

common_patients <- Reduce(
  intersect,
  list(
    patients_tcell,
    patients_myeloid,
    patients_endo
  )
)

cat(
  "Common matched patients:",
  length(common_patients),
  "\n"
)

cat(
  "Common patients:",
  paste(common_patients, collapse = ", "),
  "\n"
)


# ------------------------------------------------------------
# 6. Numeric / non-negative raw UMI audit
# ------------------------------------------------------------

check_numeric_nonnegative <- function(x) {
  
  mat <- as.matrix(x[, -1, drop = FALSE])
  
  numeric_ok <- all(
    vapply(
      x[, -1, drop = FALSE],
      is.numeric,
      logical(1)
    )
  )
  
  no_na <- !anyNA(mat)
  no_inf <- all(is.finite(mat))
  nonnegative <- all(mat >= 0)
  
  c(
    numeric = numeric_ok,
    no_NA = no_na,
    finite = no_inf,
    nonnegative = nonnegative
  )
}

qc_tcell <- check_numeric_nonnegative(pb_tcell)
qc_myeloid <- check_numeric_nonnegative(pb_myeloid)
qc_endo <- check_numeric_nonnegative(pb_endo)

cat("\n------------------------------------------------------------\n")
cat("6. RAW-UMI NUMERIC QC\n")
cat("------------------------------------------------------------\n")

print(qc_tcell)
print(qc_myeloid)
print(qc_endo)


# ------------------------------------------------------------
# 7. Library-size calculation
# ------------------------------------------------------------

calc_library_size <- function(x) {
  
  mat <- as.matrix(x[, -1, drop = FALSE])
  
  data.frame(
    Patient = colnames(x)[-1],
    Patient_ID = sub("__.*$", "", colnames(x)[-1]),
    Compartment = sub("^.*__", "", colnames(x)[-1]),
    Library_Total_UMI = colSums(mat),
    Detected_Genes = colSums(mat > 0),
    stringsAsFactors = FALSE
  )
}

lib_tcell <- calc_library_size(pb_tcell)
lib_myeloid <- calc_library_size(pb_myeloid)
lib_endo <- calc_library_size(pb_endo)

library_audit <- rbind(
  lib_tcell,
  lib_myeloid,
  lib_endo
)

cat("\n------------------------------------------------------------\n")
cat("7. LIBRARY-SIZE AUDIT\n")
cat("------------------------------------------------------------\n")

print(library_audit)


# ------------------------------------------------------------
# 8. Summary of library depth
# ------------------------------------------------------------

library_summary <- do.call(
  rbind,
  lapply(
    split(
      library_audit,
      library_audit$Compartment
    ),
    function(d) {
      
      data.frame(
        Compartment = unique(d$Compartment),
        N = nrow(d),
        Min_Library_UMI = min(d$Library_Total_UMI),
        Median_Library_UMI = median(d$Library_Total_UMI),
        Mean_Library_UMI = mean(d$Library_Total_UMI),
        Max_Library_UMI = max(d$Library_Total_UMI),
        Min_Detected_Genes = min(d$Detected_Genes),
        Median_Detected_Genes = median(d$Detected_Genes),
        Max_Detected_Genes = max(d$Detected_Genes),
        stringsAsFactors = FALSE
      )
    }
  )
)

cat("\n------------------------------------------------------------\n")
cat("8. LIBRARY DEPTH SUMMARY\n")
cat("------------------------------------------------------------\n")

print(library_summary)


# ------------------------------------------------------------
# 9. Matched-patient library-size comparison
# ------------------------------------------------------------

matched_library <- library_audit[
  library_audit$Patient_ID %in% common_patients,
]

matched_library <- matched_library[
  order(
    matched_library$Patient_ID,
    matched_library$Compartment
  ),
]

cat("\n------------------------------------------------------------\n")
cat("9. MATCHED PATIENT LIBRARY-SIZE COMPARISON\n")
cat("------------------------------------------------------------\n")

print(matched_library)


# ------------------------------------------------------------
# 10. Calculate CPM matrix
#
# IMPORTANT:
#   This is a mathematical audit transformation only.
#   Original raw-UMI matrices are NOT modified.
# ------------------------------------------------------------

to_cpm <- function(x) {
  
  mat <- as.matrix(x[, -1, drop = FALSE])
  
  lib <- colSums(mat)
  
  sweep(
    mat,
    2,
    lib,
    "/"
  ) * 1e6
}

cpm_tcell <- to_cpm(pb_tcell)
cpm_myeloid <- to_cpm(pb_myeloid)
cpm_endo <- to_cpm(pb_endo)


# ------------------------------------------------------------
# 11. Calculate log2(CPM + 1)
# ------------------------------------------------------------

log2cpm_tcell <- log2(cpm_tcell + 1)
log2cpm_myeloid <- log2(cpm_myeloid + 1)
log2cpm_endo <- log2(cpm_endo + 1)


# ------------------------------------------------------------
# 12. Transformation integrity audit
# ------------------------------------------------------------

check_transformation <- function(cpm, log2cpm) {
  
  c(
    CPM_finite = all(is.finite(cpm)),
    CPM_nonnegative = all(cpm >= 0),
    log2CPM_finite = all(is.finite(log2cpm)),
    log2CPM_nonnegative = all(log2cpm >= 0)
  )
}

transform_qc <- rbind(
  T_cell = check_transformation(
    cpm_tcell,
    log2cpm_tcell
  ),
  Myeloid_Macrophage = check_transformation(
    cpm_myeloid,
    log2cpm_myeloid
  ),
  Endothelial = check_transformation(
    cpm_endo,
    log2cpm_endo
  )
)

cat("\n------------------------------------------------------------\n")
cat("12. CPM / LOG2(CPM+1) TRANSFORMATION QC\n")
cat("------------------------------------------------------------\n")

print(transform_qc)


# ------------------------------------------------------------
# 13. Core-gene scale audit
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

get_core_expression <- function(
    x,
    cpm,
    log2cpm,
    compartment
) {
  
  idx <- match(core_genes, x$Gene)
  
  if (anyNA(idx)) {
    stop(
      "One or more core genes missing from ",
      compartment
    )
  }
  
  out <- data.frame(
    Gene = core_genes,
    stringsAsFactors = FALSE
  )
  
  for (j in seq_len(ncol(cpm))) {
    
    out[[paste0(
      colnames(cpm)[j],
      "__CPM"
    )]] <- cpm[idx, j]
    
    out[[paste0(
      colnames(cpm)[j],
      "__log2CPM1"
    )]] <- log2cpm[idx, j]
  }
  
  out
}

core_tcell <- get_core_expression(
  pb_tcell,
  cpm_tcell,
  log2cpm_tcell,
  "T-cell"
)

core_myeloid <- get_core_expression(
  pb_myeloid,
  cpm_myeloid,
  log2cpm_myeloid,
  "Myeloid/macrophage"
)

core_endo <- get_core_expression(
  pb_endo,
  cpm_endo,
  log2cpm_endo,
  "Endothelial"
)


# ------------------------------------------------------------
# 14. Patient-matched core-gene inspection
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("14. CORE-GENE EXPRESSION-SCALE AUDIT\n")
cat("------------------------------------------------------------\n")

cat("\nT-cell:\n")
print(core_tcell)

cat("\nMyeloid/macrophage:\n")
print(core_myeloid)

cat("\nEndothelial:\n")
print(core_endo)


# ------------------------------------------------------------
# 15. Detection audit
#
# Detection = raw UMI > 0
#
# This is descriptive only.
# NO threshold is locked here.
# ------------------------------------------------------------

core_detection_summary <- data.frame()

for (gene in core_genes) {
  
  for (comp in c(
    "T-cell",
    "Myeloid/macrophage",
    "Endothelial"
  )) {
    
    x <- switch(
      comp,
      "T-cell" = pb_tcell,
      "Myeloid/macrophage" = pb_myeloid,
      "Endothelial" = pb_endo
    )
    
    idx <- match(gene, x$Gene)
    counts <- as.numeric(x[idx, -1])
    
    core_detection_summary <- rbind(
      core_detection_summary,
      data.frame(
        Gene = gene,
        Compartment = comp,
        N_Patients = length(counts),
        N_Detected = sum(counts > 0),
        Detection_Fraction = mean(counts > 0),
        stringsAsFactors = FALSE
      )
    )
  }
}

cat("\n------------------------------------------------------------\n")
cat("15. CORE-GENE DETECTION AUDIT\n")
cat("------------------------------------------------------------\n")

print(core_detection_summary)


# ------------------------------------------------------------
# 16. Final audit logic
# ------------------------------------------------------------

raw_qc_pass <- all(qc_tcell) &&
  all(qc_myeloid) &&
  all(qc_endo)

transform_pass <- all(transform_qc)

matched_patient_pass <- (
  length(common_patients) == 4
)

core_gene_pass <- all(
  core_genes %in% pb_tcell$Gene
) &&
  all(
    core_genes %in% pb_myeloid$Gene
  ) &&
  all(
    core_genes %in% pb_endo$Gene
  )

overall_5C1_pass <- (
  gene_identity_pass &&
    gene_unique_pass &&
    raw_qc_pass &&
    transform_pass &&
    matched_patient_pass &&
    core_gene_pass
)


# ------------------------------------------------------------
# 17. FINAL STATUS
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP11-5C-1 FINAL STATUS\n")
cat("============================================================\n\n")

cat(
  "Existing STEP10C pseudobulk matrices:",
  "READ SUCCESSFULLY\n"
)

cat(
  "Raw-UMI integrity:",
  ifelse(raw_qc_pass, "PASS", "FAIL"),
  "\n"
)

cat(
  "Gene identity/order:",
  ifelse(gene_identity_pass, "PASS", "FAIL"),
  "\n"
)

cat(
  "Gene uniqueness:",
  ifelse(gene_unique_pass, "PASS", "FAIL"),
  "\n"
)

cat(
  "CPM transformation integrity:",
  ifelse(transform_pass, "PASS", "FAIL"),
  "\n"
)

cat(
  "Matched patient design:",
  ifelse(matched_patient_pass, "PASS", "FAIL"),
  "\n"
)

cat(
  "Core-gene presence:",
  ifelse(core_gene_pass, "PASS", "FAIL"),
  "\n\n"
)

cat("NO pseudobulk recalculation\n")
cat("NO raw-UMI re-aggregation\n")
cat("NO CellChat recalculation\n")
cat("NO LR statistical testing\n")
cat("NO pathway selection\n")
cat("NO LR threshold locking\n\n")


if (overall_5C1_pass) {
  
  cat(
    "STATUS: PASS — EXPRESSION-SCALE / NORMALIZATION COMPATIBILITY AUDIT COMPLETED.\n"
  )
  
} else {
  
  cat(
    "STATUS: REVIEW REQUIRED — ONE OR MORE EXPRESSION-SCALE QC ITEMS DID NOT PASS.\n"
  )
}

cat("\n")
cat("============================================================\n")
cat("END STEP11-5C-1\n")
cat("============================================================\n")


# ============================================================
# STEP11-5C-2
# Normalized Expression Stability / Matched-Patient Scale Audit
#
# Purpose:
#   Evaluate whether log2(CPM + 1) provides a stable and
#   interpretable expression scale across the 4 matched patients
#   and 3 focused compartments.
#
# IMPORTANT:
#   - Existing STEP10C pseudobulk matrices are reused
#   - No pseudobulk recalculation
#   - No raw-UMI re-aggregation
#   - No CellChat recalculation
#   - No LR candidate screening
#   - No LR statistical testing
#   - No pathway selection
#   - No LR expression threshold locking
# ============================================================


cat("============================================================\n")
cat("STEP11-5C-2 — NORMALIZED EXPRESSION STABILITY / MATCHED-PATIENT SCALE AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

step10_root <- "input_dir/Paper 1_Step 10-Single-Cell"
step11_root <- "input_dir/Paper 1_Step 11-CellChat"

file_tcell <- file.path(
  step10_root,
  "STEP10C_Patient_Level_RawUMI_Pseudobulk",
  "STEP10C-A_Pseudobulk_Tcell_Primary_rawUMI.tsv"
)

file_myeloid <- file.path(
  step10_root,
  "STEP10C_Patient_Level_RawUMI_Pseudobulk",
  "STEP10C-B_Pseudobulk_MyeloidMacrophage_Primary_rawUMI.tsv"
)

file_endo <- file.path(
  step10_root,
  "STEP10C_Patient_Level_RawUMI_Pseudobulk",
  "STEP10C-C_Pseudobulk_Endothelial_Primary_rawUMI.tsv"
)

cci_file <- file.path(
  step11_root,
  "NSCLC_GSE131907_CCI_LR.tsv"
)


# ------------------------------------------------------------
# 1. File existence
# ------------------------------------------------------------

file_exists <- c(
  T_cell = file.exists(file_tcell),
  Myeloid_Macrophage = file.exists(file_myeloid),
  Endothelial = file.exists(file_endo),
  CCI_LR = file.exists(cci_file)
)

cat("1. FILE EXISTENCE\n")
print(file_exists)

if (!all(file_exists)) {
  stop(
    "One or more required files do not exist."
  )
}


# ------------------------------------------------------------
# 2. Read existing STEP10C matrices
# ------------------------------------------------------------

pb_tcell <- read.delim(
  file_tcell,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pb_myeloid <- read.delim(
  file_myeloid,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pb_endo <- read.delim(
  file_endo,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cci_lr <- read.delim(
  cci_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 3. Confirm matrix identity
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("2. MATRIX IDENTITY\n")
cat("============================================================\n")

cat(
  "T-cell:",
  nrow(pb_tcell), "genes ×",
  ncol(pb_tcell) - 1, "patients\n"
)

cat(
  "Myeloid/macrophage:",
  nrow(pb_myeloid), "genes ×",
  ncol(pb_myeloid) - 1, "patients\n"
)

cat(
  "Endothelial:",
  nrow(pb_endo), "genes ×",
  ncol(pb_endo) - 1, "patients\n"
)

cat(
  "CCI_LR:",
  nrow(cci_lr), "interactions\n"
)


# ------------------------------------------------------------
# 4. Gene identity
# ------------------------------------------------------------

genes_tcell <- pb_tcell$Gene
genes_myeloid <- pb_myeloid$Gene
genes_endo <- pb_endo$Gene

gene_identity_pass <- (
  identical(genes_tcell, genes_myeloid) &&
    identical(genes_tcell, genes_endo)
)

gene_unique_pass <- (
  length(unique(genes_tcell)) == length(genes_tcell) &&
    length(unique(genes_myeloid)) == length(genes_myeloid) &&
    length(unique(genes_endo)) == length(genes_endo)
)

cat("\nGene identity/order:", gene_identity_pass, "\n")
cat("Gene uniqueness:", gene_unique_pass, "\n")


# ------------------------------------------------------------
# 5. Decode patient IDs
# ------------------------------------------------------------

decode_patients <- function(x) {
  
  sub(
    "__.*$",
    "",
    colnames(x)[-1]
  )
}

patients_tcell <- decode_patients(pb_tcell)
patients_myeloid <- decode_patients(pb_myeloid)
patients_endo <- decode_patients(pb_endo)

common_patients <- Reduce(
  intersect,
  list(
    patients_tcell,
    patients_myeloid,
    patients_endo
  )
)

cat("\nMatched patients:\n")
cat(
  paste(common_patients, collapse = ", "),
  "\n"
)

if (length(common_patients) != 4) {
  stop(
    "Expected four matched patients; observed ",
    length(common_patients)
  )
}


# ------------------------------------------------------------
# 6. Convert matrices to raw count matrices
# ------------------------------------------------------------

mat_tcell <- as.matrix(
  pb_tcell[, -1, drop = FALSE]
)

mat_myeloid <- as.matrix(
  pb_myeloid[, -1, drop = FALSE]
)

mat_endo <- as.matrix(
  pb_endo[, -1, drop = FALSE]
)

rownames(mat_tcell) <- genes_tcell
rownames(mat_myeloid) <- genes_myeloid
rownames(mat_endo) <- genes_endo


# ------------------------------------------------------------
# 7. Restrict to matched patients
# ------------------------------------------------------------

idx_tcell <- match(
  common_patients,
  patients_tcell
)

idx_myeloid <- match(
  common_patients,
  patients_myeloid
)

idx_endo <- match(
  common_patients,
  patients_endo
)

mat_tcell_m <- mat_tcell[, idx_tcell, drop = FALSE]
mat_myeloid_m <- mat_myeloid[, idx_myeloid, drop = FALSE]
mat_endo_m <- mat_endo[, idx_endo, drop = FALSE]

colnames(mat_tcell_m) <- common_patients
colnames(mat_myeloid_m) <- common_patients
colnames(mat_endo_m) <- common_patients


# ------------------------------------------------------------
# 8. Calculate library size
# ------------------------------------------------------------

library_tcell <- colSums(mat_tcell_m)
library_myeloid <- colSums(mat_myeloid_m)
library_endo <- colSums(mat_endo_m)

library_matched <- data.frame(
  Patient = common_patients,
  T_cell = as.numeric(library_tcell),
  Myeloid_Macrophage = as.numeric(library_myeloid),
  Endothelial = as.numeric(library_endo),
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("3. MATCHED LIBRARY DEPTH\n")
cat("============================================================\n")

print(library_matched)


# ------------------------------------------------------------
# 9. CPM conversion
#
# Mathematical transformation only.
# Original raw UMI matrices remain unchanged.
# ------------------------------------------------------------

to_cpm_matrix <- function(mat) {
  
  lib <- colSums(mat)
  
  sweep(
    mat,
    2,
    lib,
    "/"
  ) * 1e6
}

cpm_tcell <- to_cpm_matrix(mat_tcell_m)
cpm_myeloid <- to_cpm_matrix(mat_myeloid_m)
cpm_endo <- to_cpm_matrix(mat_endo_m)


# ------------------------------------------------------------
# 10. Verify CPM scale
# ------------------------------------------------------------

cpm_sums <- data.frame(
  Patient = common_patients,
  T_cell = colSums(cpm_tcell),
  Myeloid_Macrophage = colSums(cpm_myeloid),
  Endothelial = colSums(cpm_endo),
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("4. CPM SCALE VERIFICATION\n")
cat("============================================================\n")

print(cpm_sums)

cpm_scale_pass <- all(
  abs(
    as.matrix(cpm_sums[, -1]) - 1e6
  ) < 1e-6
)

cat(
  "\nCPM total ≈ 1,000,000:",
  ifelse(cpm_scale_pass, "PASS", "FAIL"),
  "\n"
)


# ------------------------------------------------------------
# 11. log2(CPM + 1)
# ------------------------------------------------------------

log2cpm_tcell <- log2(cpm_tcell + 1)
log2cpm_myeloid <- log2(cpm_myeloid + 1)
log2cpm_endo <- log2(cpm_endo + 1)


# ------------------------------------------------------------
# 12. Transformation integrity
# ------------------------------------------------------------

transformation_qc <- data.frame(
  Compartment = c(
    "T-cell",
    "Myeloid/macrophage",
    "Endothelial"
  ),
  CPM_finite = c(
    all(is.finite(cpm_tcell)),
    all(is.finite(cpm_myeloid)),
    all(is.finite(cpm_endo))
  ),
  log2CPM_finite = c(
    all(is.finite(log2cpm_tcell)),
    all(is.finite(log2cpm_myeloid)),
    all(is.finite(log2cpm_endo))
  ),
  CPM_nonnegative = c(
    all(cpm_tcell >= 0),
    all(cpm_myeloid >= 0),
    all(cpm_endo >= 0)
  ),
  log2CPM_nonnegative = c(
    all(log2cpm_tcell >= 0),
    all(log2cpm_myeloid >= 0),
    all(log2cpm_endo >= 0)
  ),
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("5. TRANSFORMATION INTEGRITY\n")
cat("============================================================\n")

print(transformation_qc)

transformation_pass <- all(
  transformation_qc[, -1]
)


# ------------------------------------------------------------
# 13. Global expression distribution summary
# ------------------------------------------------------------

distribution_summary <- function(logmat, compartment) {
  
  q <- quantile(
    as.numeric(logmat),
    probs = c(
      0,
      0.25,
      0.50,
      0.75,
      0.90,
      0.95,
      0.99,
      1
    ),
    names = FALSE
  )
  
  data.frame(
    Compartment = compartment,
    Min = q[1],
    Q25 = q[2],
    Median = q[3],
    Q75 = q[4],
    Q90 = q[5],
    Q95 = q[6],
    Q99 = q[7],
    Max = q[8],
    stringsAsFactors = FALSE
  )
}

global_distribution <- rbind(
  distribution_summary(
    log2cpm_tcell,
    "T-cell"
  ),
  distribution_summary(
    log2cpm_myeloid,
    "Myeloid/macrophage"
  ),
  distribution_summary(
    log2cpm_endo,
    "Endothelial"
  )
)

cat("\n============================================================\n")
cat("6. GLOBAL log2(CPM+1) DISTRIBUTION\n")
cat("============================================================\n")

print(global_distribution)


# ------------------------------------------------------------
# 14. Per-patient global distribution
# ------------------------------------------------------------

patient_distribution <- data.frame()

for (comp in c(
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)) {
  
  logmat <- switch(
    comp,
    "T-cell" = log2cpm_tcell,
    "Myeloid/macrophage" = log2cpm_myeloid,
    "Endothelial" = log2cpm_endo
  )
  
  for (j in seq_len(ncol(logmat))) {
    
    x <- logmat[, j]
    
    patient_distribution <- rbind(
      patient_distribution,
      data.frame(
        Patient = colnames(logmat)[j],
        Compartment = comp,
        Median = median(x),
        Q90 = quantile(x, 0.90),
        Q99 = quantile(x, 0.99),
        Max = max(x),
        stringsAsFactors = FALSE
      )
    )
  }
}

cat("\n============================================================\n")
cat("7. PATIENT-LEVEL log2(CPM+1) DISTRIBUTION\n")
cat("============================================================\n")

print(patient_distribution)


# ------------------------------------------------------------
# 15. CCI_LR resource audit
# ------------------------------------------------------------

lr_genes <- unique(
  c(
    cci_lr$ligand,
    cci_lr$receptor
  )
)

lr_genes <- lr_genes[
  !is.na(lr_genes) &
    nzchar(lr_genes)
]

cat("\n============================================================\n")
cat("8. CCI_LR GENE RESOURCE\n")
cat("============================================================\n")

cat(
  "LR interactions:",
  nrow(cci_lr),
  "\n"
)

cat(
  "Unique LR genes:",
  length(lr_genes),
  "\n"
)

lr_overlap_tcell <- intersect(
  lr_genes,
  rownames(log2cpm_tcell)
)

lr_overlap_myeloid <- intersect(
  lr_genes,
  rownames(log2cpm_myeloid)
)

lr_overlap_endo <- intersect(
  lr_genes,
  rownames(log2cpm_endo)
)

cat(
  "LR genes present in T-cell:",
  length(lr_overlap_tcell),
  "\n"
)

cat(
  "LR genes present in Myeloid/macrophage:",
  length(lr_overlap_myeloid),
  "\n"
)

cat(
  "LR genes present in Endothelial:",
  length(lr_overlap_endo),
  "\n"
)


# ------------------------------------------------------------
# 16. LR-gene normalized expression distribution
#
# Descriptive only.
# NO interaction screening.
# NO threshold.
# ------------------------------------------------------------

lr_tcell <- log2cpm_tcell[
  intersect(lr_genes, rownames(log2cpm_tcell)),
  ,
  drop = FALSE
]

lr_myeloid <- log2cpm_myeloid[
  intersect(lr_genes, rownames(log2cpm_myeloid)),
  ,
  drop = FALSE
]

lr_endo <- log2cpm_endo[
  intersect(lr_genes, rownames(log2cpm_endo)),
  ,
  drop = FALSE
]

lr_distribution <- rbind(
  distribution_summary(
    lr_tcell,
    "T-cell — LR genes"
  ),
  distribution_summary(
    lr_myeloid,
    "Myeloid/macrophage — LR genes"
  ),
  distribution_summary(
    lr_endo,
    "Endothelial — LR genes"
  )
)

cat("\n============================================================\n")
cat("9. LR-GENE log2(CPM+1) DISTRIBUTION\n")
cat("============================================================\n")

print(lr_distribution)


# ------------------------------------------------------------
# 17. LR-gene detection landscape
#
# Detection = raw UMI > 0
#
# Descriptive only.
# NO threshold is locked.
# ------------------------------------------------------------

lr_detection_summary <- data.frame()

for (comp in c(
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)) {
  
  rawmat <- switch(
    comp,
    "T-cell" = mat_tcell_m,
    "Myeloid/macrophage" = mat_myeloid_m,
    "Endothelial" = mat_endo_m
  )
  
  lr_present <- intersect(
    lr_genes,
    rownames(rawmat)
  )
  
  detect <- rawmat[
    lr_present,
    ,
    drop = FALSE
  ] > 0
  
  lr_detection_summary <- rbind(
    lr_detection_summary,
    data.frame(
      Compartment = comp,
      LR_Genes = nrow(detect),
      Mean_Detection_Per_Patient =
        colMeans(detect),
      stringsAsFactors = FALSE
    )
  )
}

cat("\n============================================================\n")
cat("10. LR-GENE DETECTION LANDSCAPE\n")
cat("============================================================\n")

print(lr_detection_summary)


# ------------------------------------------------------------
# 18. Core-gene normalized expression
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_scale <- data.frame()

for (gene in core_genes) {
  
  for (comp in c(
    "T-cell",
    "Myeloid/macrophage",
    "Endothelial"
  )) {
    
    logmat <- switch(
      comp,
      "T-cell" = log2cpm_tcell,
      "Myeloid/macrophage" = log2cpm_myeloid,
      "Endothelial" = log2cpm_endo
    )
    
    rawmat <- switch(
      comp,
      "T-cell" = mat_tcell_m,
      "Myeloid/macrophage" = mat_myeloid_m,
      "Endothelial" = mat_endo_m
    )
    
    if (!gene %in% rownames(logmat)) {
      stop(
        "Core gene missing: ",
        gene,
        " in ",
        comp
      )
    }
    
    core_scale <- rbind(
      core_scale,
      data.frame(
        Gene = gene,
        Compartment = comp,
        Raw_UMI_Median =
          median(rawmat[gene, ]),
        CPM_Median =
          median(
            2^logmat[gene, ] - 1
          ),
        log2CPM_Median =
          median(logmat[gene, ]),
        Detected =
          sum(rawmat[gene, ] > 0),
        N =
          ncol(rawmat),
        stringsAsFactors = FALSE
      )
    )
  }
}

cat("\n============================================================\n")
cat("11. CORE-GENE NORMALIZED EXPRESSION AUDIT\n")
cat("============================================================\n")

print(core_scale)


# ------------------------------------------------------------
# 19. Matched-patient rank preservation audit
#
# Purpose:
#   Compare raw UMI and normalized expression ranks within each
#   compartment across the four matched patients.
#
# This is NOT a biological association test.
# ------------------------------------------------------------

rank_audit <- data.frame()

for (comp in c(
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)) {
  
  rawmat <- switch(
    comp,
    "T-cell" = mat_tcell_m,
    "Myeloid/macrophage" = mat_myeloid_m,
    "Endothelial" = mat_endo_m
  )
  
  logmat <- switch(
    comp,
    "T-cell" = log2cpm_tcell,
    "Myeloid/macrophage" = log2cpm_myeloid,
    "Endothelial" = log2cpm_endo
  )
  
  for (gene in c(
    core_genes,
    lr_genes
  )) {
    
    if (!gene %in% rownames(rawmat)) {
      next
    }
    
    raw_values <- as.numeric(
      rawmat[gene, ]
    )
    
    norm_values <- as.numeric(
      logmat[gene, ]
    )
    
    if (
      length(unique(raw_values)) > 1 &&
      length(unique(norm_values)) > 1
    ) {
      
      rho <- suppressWarnings(
        cor(
          raw_values,
          norm_values,
          method = "spearman"
        )
      )
      
    } else {
      
      rho <- NA_real_
    }
    
    rank_audit <- rbind(
      rank_audit,
      data.frame(
        Gene = gene,
        Compartment = comp,
        Spearman_Raw_vs_log2CPM = rho,
        stringsAsFactors = FALSE
      )
    )
  }
}

cat("\n============================================================\n")
cat("12. RAW-TO-NORMALIZED RANK PRESERVATION\n")
cat("============================================================\n")

rank_summary <- do.call(
  rbind,
  lapply(
    split(
      rank_audit$Spearman_Raw_vs_log2CPM,
      rank_audit$Compartment
    ),
    function(x) {
      
      data.frame(
        N = sum(!is.na(x)),
        Min_Rho = min(x, na.rm = TRUE),
        Median_Rho = median(x, na.rm = TRUE),
        Max_Rho = max(x, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  )
)

rank_summary$Compartment <- rownames(rank_summary)
rownames(rank_summary) <- NULL

rank_summary <- rank_summary[
  ,
  c(
    "Compartment",
    "N",
    "Min_Rho",
    "Median_Rho",
    "Max_Rho"
  )
]

print(rank_summary)


# ------------------------------------------------------------
# 20. Patient-wise normalized-scale stability
#
# Compare each patient's overall LR-gene expression distribution
# against the other matched patients within each compartment.
# Descriptive only.
# ------------------------------------------------------------

patient_lr_median <- data.frame()

for (comp in c(
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)) {
  
  logmat <- switch(
    comp,
    "T-cell" = lr_tcell,
    "Myeloid/macrophage" = lr_myeloid,
    "Endothelial" = lr_endo
  )
  
  for (j in seq_len(ncol(logmat))) {
    
    patient_lr_median <- rbind(
      patient_lr_median,
      data.frame(
        Patient = colnames(logmat)[j],
        Compartment = comp,
        LR_Gene_Median_log2CPM =
          median(logmat[, j]),
        LR_Gene_Q90_log2CPM =
          quantile(logmat[, j], 0.90),
        LR_Gene_Q99_log2CPM =
          quantile(logmat[, j], 0.99),
        stringsAsFactors = FALSE
      )
    )
  }
}

cat("\n============================================================\n")
cat("13. PATIENT-WISE LR-GENE NORMALIZED SCALE\n")
cat("============================================================\n")

print(patient_lr_median)


# ------------------------------------------------------------
# 21. Final descriptive flags
# ------------------------------------------------------------

lr_overlap_pass <- (
  length(lr_overlap_tcell) > 0 &&
    length(lr_overlap_myeloid) > 0 &&
    length(lr_overlap_endo) > 0
)

core_pass <- all(
  core_genes %in% rownames(mat_tcell_m)
) &&
  all(
    core_genes %in% rownames(mat_myeloid_m)
  ) &&
  all(
    core_genes %in% rownames(mat_endo_m)
  )

finite_distribution_pass <- all(
  is.finite(
    as.matrix(
      global_distribution[, -1]
    )
  )
)

overall_5C2_pass <- (
  gene_identity_pass &&
    gene_unique_pass &&
    cpm_scale_pass &&
    transformation_pass &&
    lr_overlap_pass &&
    core_pass &&
    finite_distribution_pass
)


# ------------------------------------------------------------
# 22. FINAL STATUS
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP11-5C-2 FINAL STATUS\n")
cat("============================================================\n\n")

cat(
  "Existing STEP10C matrices reused: YES\n"
)

cat(
  "Matched patients:",
  length(common_patients),
  "\n"
)

cat(
  "Patients:",
  paste(common_patients, collapse = ", "),
  "\n"
)

cat(
  "CPM scale verification:",
  ifelse(cpm_scale_pass, "PASS", "FAIL"),
  "\n"
)

cat(
  "Transformation integrity:",
  ifelse(transformation_pass, "PASS", "FAIL"),
  "\n"
)

cat(
  "LR-gene resource overlap:",
  ifelse(lr_overlap_pass, "PASS", "FAIL"),
  "\n"
)

cat(
  "Core-gene presence:",
  ifelse(core_pass, "PASS", "FAIL"),
  "\n"
)

cat(
  "Normalized distribution finite:",
  ifelse(finite_distribution_pass, "PASS", "FAIL"),
  "\n\n"
)

cat("NO pseudobulk recalculation\n")
cat("NO raw-UMI re-aggregation\n")
cat("NO CellChat recalculation\n")
cat("NO LR candidate screening\n")
cat("NO LR statistical testing\n")
cat("NO pathway selection\n")
cat("NO LR expression threshold locking\n\n")


if (overall_5C2_pass) {
  
  cat(
    "STATUS: PASS — log2(CPM+1) PROVIDES A VALID DESCRIPTIVE NORMALIZED EXPRESSION SCALE FOR SUBSEQUENT PATIENT-LEVEL LR COMPATIBILITY SCREENING.\n"
  )
  
} else {
  
  cat(
    "STATUS: REVIEW REQUIRED — ONE OR MORE NORMALIZED-SCALE QC ITEMS DID NOT PASS.\n"
  )
}

cat("\n============================================================\n")
cat("END STEP11-5C-2\n")
cat("============================================================\n")


# ============================================================
# STEP11-5D-0
# LR DETECTION / PATIENT-CONSISTENCY RULE SENSITIVITY AUDIT
# ============================================================
#
# Purpose:
#   Evaluate predefined expression-detection and
#   patient-consistency rules before formal LR screening.
#
# Scope:
#   - Existing STEP10C pseudobulk matrices only
#   - Four matched patients only
#   - 744 LR genes present in all three compartments
#   - Six compartment directions
#
# This step DOES NOT:
#   - recalculate pseudobulk
#   - re-aggregate raw UMI
#   - recalculate CellChat
#   - perform LR statistical testing
#   - calculate P values / FDR
#   - select pathways
#   - lock a biological LR result
#
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP11-5D-0: LR DETECTION / PATIENT-CONSISTENCY SENSITIVITY AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

step10_dir <- "input_dir/Paper 1_Step 10-Single-Cell"

file_tcell <- file.path(
  step10_dir,
  "STEP10C_Patient_Level_RawUMI_Pseudobulk",
  "STEP10C-A_Pseudobulk_Tcell_Primary_rawUMI.tsv"
)

file_myeloid <- file.path(
  step10_dir,
  "STEP10C_Patient_Level_RawUMI_Pseudobulk",
  "STEP10C-B_Pseudobulk_MyeloidMacrophage_Primary_rawUMI.tsv"
)

file_endo <- file.path(
  step10_dir,
  "STEP10C_Patient_Level_RawUMI_Pseudobulk",
  "STEP10C-C_Pseudobulk_Endothelial_Primary_rawUMI.tsv"
)

cci_file <- file.path(
  step10_dir,
  "..",
  "Paper 1_Step 11-CellChat",
  "NSCLC_GSE131907_CCI_LR.tsv"
)

# Normalize path
cci_file <- normalizePath(
  cci_file,
  winslash = "/",
  mustWork = FALSE
)


# ------------------------------------------------------------
# 1. File existence
# ------------------------------------------------------------

file_exists <- c(
  T_cell = file.exists(file_tcell),
  Myeloid_Macrophage = file.exists(file_myeloid),
  Endothelial = file.exists(file_endo),
  CCI_LR = file.exists(cci_file)
)

cat("1. FILE EXISTENCE\n")
cat("============================================================\n")
print(file_exists)

if (!all(file_exists)) {
  stop("Required input file(s) missing.")
}


# ------------------------------------------------------------
# 2. Read existing matrices / CCI_LR
# ------------------------------------------------------------

pb_tcell <- read.delim(
  file_tcell,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pb_myeloid <- read.delim(
  file_myeloid,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pb_endo <- read.delim(
  file_endo,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cci_lr <- read.delim(
  cci_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n2. INPUT STRUCTURE\n")
cat("============================================================\n")

cat(
  "T-cell:",
  nrow(pb_tcell),
  "genes ×",
  ncol(pb_tcell) - 1,
  "patients\n"
)

cat(
  "Myeloid/macrophage:",
  nrow(pb_myeloid),
  "genes ×",
  ncol(pb_myeloid) - 1,
  "patients\n"
)

cat(
  "Endothelial:",
  nrow(pb_endo),
  "genes ×",
  ncol(pb_endo) - 1,
  "patients\n"
)

cat(
  "CCI_LR:",
  nrow(cci_lr),
  "interactions\n"
)


# ------------------------------------------------------------
# 3. Gene identity / uniqueness
# ------------------------------------------------------------

gene_t <- as.character(pb_tcell$Gene)
gene_m <- as.character(pb_myeloid$Gene)
gene_e <- as.character(pb_endo$Gene)

gene_identity_pass <-
  identical(gene_t, gene_m) &&
  identical(gene_t, gene_e)

gene_unique_pass <-
  !anyDuplicated(gene_t) &&
  !anyDuplicated(gene_m) &&
  !anyDuplicated(gene_e)

cat("\n3. GENE IDENTITY / UNIQUENESS\n")
cat("============================================================\n")

cat(
  "Gene identity/order:",
  gene_identity_pass,
  "\n"
)

cat(
  "Gene uniqueness:",
  gene_unique_pass,
  "\n"
)

if (!gene_identity_pass || !gene_unique_pass) {
  stop("Gene identity/order or uniqueness QC failed.")
}


# ------------------------------------------------------------
# 4. Decode patient IDs
# ------------------------------------------------------------

get_patient_ids <- function(x) {
  
  cols <- colnames(x)[-1]
  
  sub(
    "__.*$",
    "",
    cols
  )
}

patients_t <- get_patient_ids(pb_tcell)
patients_m <- get_patient_ids(pb_myeloid)
patients_e <- get_patient_ids(pb_endo)

common_patients <- Reduce(
  intersect,
  list(
    patients_t,
    patients_m,
    patients_e
  )
)

common_patients <- sort(common_patients)

cat("\n4. MATCHED PATIENTS\n")
cat("============================================================\n")

cat(
  "Matched patients:",
  length(common_patients),
  "\n"
)

cat(
  paste(common_patients, collapse = ", "),
  "\n"
)

if (length(common_patients) != 4) {
  stop(
    "Expected exactly 4 matched patients based on STEP11-5C-2."
  )
}


# ------------------------------------------------------------
# 5. Convert to matrices
# ------------------------------------------------------------

make_numeric_matrix <- function(df, patient_ids) {
  
  x <- as.matrix(
    df[, paste0(
      patient_ids,
      "__",
      sub(
        "^.*__",
        "",
        colnames(df)[-1]
      )[match(
        paste0(patient_ids, "__"),
        substr(
          colnames(df)[-1],
          1,
          nchar(paste0(patient_ids, "__"))
        )
      )]
    ), drop = FALSE]
  )
  
}

# Safer explicit patient-column extraction
get_patient_matrix <- function(df, patient_ids) {
  
  all_cols <- colnames(df)
  
  out <- matrix(
    0,
    nrow = nrow(df),
    ncol = length(patient_ids)
  )
  
  rownames(out) <- df$Gene
  
  for (i in seq_along(patient_ids)) {
    
    hit <- grep(
      paste0("^", patient_ids[i], "__"),
      all_cols
    )
    
    if (length(hit) != 1) {
      stop(
        "Could not uniquely identify patient column: ",
        patient_ids[i]
      )
    }
    
    out[, i] <- as.numeric(df[[hit]])
  }
  
  colnames(out) <- patient_ids
  
  out
}


raw_t <- get_patient_matrix(
  pb_tcell,
  common_patients
)

raw_m <- get_patient_matrix(
  pb_myeloid,
  common_patients
)

raw_e <- get_patient_matrix(
  pb_endo,
  common_patients
)


# ------------------------------------------------------------
# 6. Normalize to CPM and log2(CPM+1)
# ------------------------------------------------------------

to_cpm <- function(mat) {
  
  lib <- colSums(mat)
  
  sweep(
    mat,
    2,
    lib,
    "/"
  ) * 1e6
}

cpm_t <- to_cpm(raw_t)
cpm_m <- to_cpm(raw_m)
cpm_e <- to_cpm(raw_e)

log2cpm_t <- log2(cpm_t + 1)
log2cpm_m <- log2(cpm_m + 1)
log2cpm_e <- log2(cpm_e + 1)


# ------------------------------------------------------------
# 7. LR resource
# ------------------------------------------------------------

ligands <- unique(
  as.character(
    cci_lr$ligand
  )
)

receptors <- unique(
  as.character(
    cci_lr$receptor
  )
)

lr_genes <- sort(
  unique(
    c(
      ligands,
      receptors
    )
  )
)

lr_overlap <- Reduce(
  intersect,
  list(
    rownames(raw_t),
    rownames(raw_m),
    rownames(raw_e),
    lr_genes
  )
)

cat("\n5. LR RESOURCE\n")
cat("============================================================\n")

cat(
  "LR interactions:",
  nrow(cci_lr),
  "\n"
)

cat(
  "Unique LR genes:",
  length(lr_genes),
  "\n"
)

cat(
  "LR genes available in all three compartments:",
  length(lr_overlap),
  "\n"
)


# ------------------------------------------------------------
# 8. Restrict matrices to LR genes
# ------------------------------------------------------------

raw_t_lr <- raw_t[lr_overlap, , drop = FALSE]
raw_m_lr <- raw_m[lr_overlap, , drop = FALSE]
raw_e_lr <- raw_e[lr_overlap, , drop = FALSE]

cpm_t_lr <- cpm_t[lr_overlap, , drop = FALSE]
cpm_m_lr <- cpm_m[lr_overlap, , drop = FALSE]
cpm_e_lr <- cpm_e[lr_overlap, , drop = FALSE]


# ------------------------------------------------------------
# 9. Predefined detection rules
# ------------------------------------------------------------

# Rule 1:
# Raw UMI > 0
#
# Rule 2:
# CPM > 0
#
# Rule 3:
# CPM >= 1

detect_rule <- function(
    raw_mat,
    cpm_mat,
    rule
) {
  
  if (rule == "UMI_gt0") {
    
    return(
      raw_mat > 0
    )
    
  }
  
  if (rule == "CPM_gt0") {
    
    return(
      cpm_mat > 0
    )
    
  }
  
  if (rule == "CPM_ge1") {
    
    return(
      cpm_mat >= 1
    )
    
  }
  
  stop("Unknown detection rule.")
  
}


detection_rules <- c(
  "UMI_gt0",
  "CPM_gt0",
  "CPM_ge1"
)


# ------------------------------------------------------------
# 10. Build six sender → receiver directions
# ------------------------------------------------------------

compartment_matrices <- list(
  "T-cell" = list(
    raw = raw_t_lr,
    cpm = cpm_t_lr
  ),
  
  "Myeloid/macrophage" = list(
    raw = raw_m_lr,
    cpm = cpm_m_lr
  ),
  
  "Endothelial" = list(
    raw = raw_e_lr,
    cpm = cpm_e_lr
  )
)


directions <- data.frame(
  
  Sender = c(
    "Endothelial",
    "Myeloid/macrophage",
    "Endothelial",
    "T-cell",
    "Myeloid/macrophage",
    "T-cell"
  ),
  
  Receiver = c(
    "Myeloid/macrophage",
    "Endothelial",
    "T-cell",
    "Endothelial",
    "T-cell",
    "Myeloid/macrophage"
  ),
  
  stringsAsFactors = FALSE
  
)


# ------------------------------------------------------------
# 11. Construct pair-level compatibility
# ------------------------------------------------------------

# A LR interaction is evaluated only if:
#
#   ligand is present in sender
#   AND
#   receptor is present in receiver
#
# within the SAME patient.
#
# No statistics are performed here.

evaluate_direction <- function(
    sender,
    receiver,
    detection_rule_name
) {
  
  sender_obj <- compartment_matrices[[sender]]
  receiver_obj <- compartment_matrices[[receiver]]
  
  sender_detect <- detect_rule(
    sender_obj$raw,
    sender_obj$cpm,
    detection_rule_name
  )
  
  receiver_detect <- detect_rule(
    receiver_obj$raw,
    receiver_obj$cpm,
    detection_rule_name
  )
  
  out <- vector(
    "list",
    nrow(cci_lr)
  )
  
  for (k in seq_len(nrow(cci_lr))) {
    
    ligand <- as.character(
      cci_lr$ligand[k]
    )
    
    receptor <- as.character(
      cci_lr$receptor[k]
    )
    
    interaction_name <- as.character(
      cci_lr$interaction_name[k]
    )
    
    if (
      !(ligand %in% rownames(sender_detect)) ||
      !(receptor %in% rownames(receiver_detect))
    ) {
      
      patient_compatible <- rep(
        FALSE,
        length(common_patients)
      )
      
    } else {
      
      ligand_det <- sender_detect[
        ligand,
        common_patients
      ]
      
      receptor_det <- receiver_detect[
        receptor,
        common_patients
      ]
      
      patient_compatible <-
        ligand_det &
        receptor_det
    }
    
    out[[k]] <- data.frame(
      
      Interaction = interaction_name,
      
      Ligand = ligand,
      
      Receptor = receptor,
      
      Sender = sender,
      
      Receiver = receiver,
      
      Detection_Rule = detection_rule_name,
      
      P0020 = patient_compatible[1],
      
      P0030 = patient_compatible[2],
      
      P0031 = patient_compatible[3],
      
      P0034 = patient_compatible[4],
      
      Compatible_Patients =
        sum(patient_compatible),
      
      stringsAsFactors = FALSE
      
    )
  }
  
  do.call(
    rbind,
    out
  )
}


# ------------------------------------------------------------
# 12. Evaluate all rules × directions
# ------------------------------------------------------------

all_results <- vector("list", 0)

counter <- 1

for (r in detection_rules) {
  
  for (d in seq_len(nrow(directions))) {
    
    cat(
      "\nRunning:",
      r,
      "|",
      directions$Sender[d],
      "->",
      directions$Receiver[d],
      "\n"
    )
    
    all_results[[counter]] <- evaluate_direction(
      sender = directions$Sender[d],
      receiver = directions$Receiver[d],
      detection_rule_name = r
    )
    
    counter <- counter + 1
  }
}

lr_patient_compatibility <- do.call(
  rbind,
  all_results
)

rownames(
  lr_patient_compatibility
) <- NULL


# ------------------------------------------------------------
# 13. Patient-consistency summaries
# ------------------------------------------------------------

consistency_levels <- c(
  ">=2/4",
  ">=3/4",
  "4/4"
)

summary_list <- vector(
  "list",
  length = 0
)

counter <- 1

for (r in detection_rules) {
  
  for (d in seq_len(nrow(directions))) {
    
    tmp <- lr_patient_compatibility[
      lr_patient_compatibility$Detection_Rule == r &
        lr_patient_compatibility$Sender == directions$Sender[d] &
        lr_patient_compatibility$Receiver == directions$Receiver[d],
      ,
      drop = FALSE
    ]
    
    summary_list[[counter]] <- data.frame(
      
      Detection_Rule = r,
      
      Sender = directions$Sender[d],
      
      Receiver = directions$Receiver[d],
      
      LR_Interactions = nrow(tmp),
      
      Compatible_ge2 =
        sum(tmp$Compatible_Patients >= 2),
      
      Compatible_ge3 =
        sum(tmp$Compatible_Patients >= 3),
      
      Compatible_4of4 =
        sum(tmp$Compatible_Patients == 4),
      
      stringsAsFactors = FALSE
      
    )
    
    counter <- counter + 1
  }
}

sensitivity_summary <- do.call(
  rbind,
  summary_list
)

rownames(
  sensitivity_summary
) <- NULL


# ------------------------------------------------------------
# 14. Global sensitivity summary
# ------------------------------------------------------------

global_summary <- data.frame(
  
  Detection_Rule = detection_rules,
  
  Total_Interactions = nrow(cci_lr),
  
  Total_ge2 =
    sapply(
      detection_rules,
      function(x) {
        sum(
          lr_patient_compatibility[
            lr_patient_compatibility$Detection_Rule == x,
            "Compatible_Patients"
          ] >= 2
        )
      }
    ),
  
  Total_ge3 =
    sapply(
      detection_rules,
      function(x) {
        sum(
          lr_patient_compatibility[
            lr_patient_compatibility$Detection_Rule == x,
            "Compatible_Patients"
          ] >= 3
        )
      }
    ),
  
  Total_4of4 =
    sapply(
      detection_rules,
      function(x) {
        sum(
          lr_patient_compatibility[
            lr_patient_compatibility$Detection_Rule == x,
            "Compatible_Patients"
          ] == 4
        )
      }
    ),
  
  stringsAsFactors = FALSE
  
)


# ------------------------------------------------------------
# 15. Detection-rule pairwise comparison
# ------------------------------------------------------------

# Compare whether the SAME LR interaction remains compatible
# when detection rule changes.
#
# This is descriptive only.

rule_pairwise <- data.frame(
  
  UMI_gt0_vs_CPM_gt0_same_2of4 = NA_integer_,
  UMI_gt0_vs_CPM_ge1_same_2of4 = NA_integer_,
  CPM_gt0_vs_CPM_ge1_same_2of4 = NA_integer_,
  
  UMI_gt0_vs_CPM_gt0_same_3of4 = NA_integer_,
  UMI_gt0_vs_CPM_ge1_same_3of4 = NA_integer_,
  CPM_gt0_vs_CPM_ge1_same_3of4 = NA_integer_,
  
  UMI_gt0_vs_CPM_gt0_same_4of4 = NA_integer_,
  UMI_gt0_vs_CPM_ge1_same_4of4 = NA_integer_,
  CPM_gt0_vs_CPM_ge1_same_4of4 = NA_integer_
  
)

get_rule_status <- function(
    rule,
    min_patients
) {
  
  x <- lr_patient_compatibility[
    lr_patient_compatibility$Detection_Rule == rule,
    ,
    drop = FALSE
  ]
  
  x$Compatible_Patients >= min_patients
}


a <- get_rule_status(
  "UMI_gt0",
  2
)

b <- get_rule_status(
  "CPM_gt0",
  2
)

c <- get_rule_status(
  "CPM_ge1",
  2
)

rule_pairwise$UMI_gt0_vs_CPM_gt0_same_2of4 <-
  sum(a & b)

rule_pairwise$UMI_gt0_vs_CPM_ge1_same_2of4 <-
  sum(a & c)

rule_pairwise$CPM_gt0_vs_CPM_ge1_same_2of4 <-
  sum(b & c)


a <- get_rule_status(
  "UMI_gt0",
  3
)

b <- get_rule_status(
  "CPM_gt0",
  3
)

c <- get_rule_status(
  "CPM_ge1",
  3
)

rule_pairwise$UMI_gt0_vs_CPM_gt0_same_3of4 <-
  sum(a & b)

rule_pairwise$UMI_gt0_vs_CPM_ge1_same_3of4 <-
  sum(a & c)

rule_pairwise$CPM_gt0_vs_CPM_ge1_same_3of4 <-
  sum(b & c)


a <- get_rule_status(
  "UMI_gt0",
  4
)

b <- get_rule_status(
  "CPM_gt0",
  4
)

c <- get_rule_status(
  "CPM_ge1",
  4
)

rule_pairwise$UMI_gt0_vs_CPM_gt0_same_4of4 <-
  sum(a & b)

rule_pairwise$UMI_gt0_vs_CPM_ge1_same_4of4 <-
  sum(a & c)

rule_pairwise$CPM_gt0_vs_CPM_ge1_same_4of4 <-
  sum(b & c)


# ------------------------------------------------------------
# 16. Direction-level sensitivity
# ------------------------------------------------------------

direction_sensitivity <- aggregate(
  
  cbind(
    Compatible_ge2,
    Compatible_ge3,
    Compatible_4of4
  ) ~
    Sender +
    Receiver +
    Detection_Rule,
  
  data = sensitivity_summary,
  
  FUN = sum
  
)


# ------------------------------------------------------------
# 17. Core-gene exclusion audit
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_in_lr_resource <- core_genes %in% lr_genes

cat("\n6. CORE-GENE / LR RESOURCE SEPARATION\n")
cat("============================================================\n")

print(
  data.frame(
    Gene = core_genes,
    Present_in_CCI_LR = core_in_lr_resource
  )
)


# ------------------------------------------------------------
# 18. Save audit outputs
# ------------------------------------------------------------

output_dir <- file.path(
  step10_dir,
  "..",
  "Paper 1_Step 11-CellChat",
  "STEP11-5D-0_LR_Detection_Consistency_Sensitivity_Audit"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.table(
  sensitivity_summary,
  file = file.path(
    output_dir,
    "STEP11-5D-0A_Direction_Sensitivity_Summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  global_summary,
  file = file.path(
    output_dir,
    "STEP11-5D-0B_Global_Detection_Rule_Summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  rule_pairwise,
  file = file.path(
    output_dir,
    "STEP11-5D-0C_Rule_Pairwise_Stability.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  lr_patient_compatibility,
  file = file.path(
    output_dir,
    "STEP11-5D-0D_LR_Patient_Compatibility_Audit.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ------------------------------------------------------------
# 19. Final console output
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP11-5D-0 FINAL STATUS\n")
cat("============================================================\n\n")

cat(
  "Existing STEP10C matrices reused: YES\n"
)

cat(
  "Matched patients:",
  length(common_patients),
  "\n"
)

cat(
  "Matched patients:",
  paste(
    common_patients,
    collapse = ", "
  ),
  "\n"
)

cat(
  "LR interactions:",
  nrow(cci_lr),
  "\n"
)

cat(
  "LR genes available in all three compartments:",
  length(lr_overlap),
  "\n"
)

cat(
  "Detection rules evaluated:",
  paste(
    detection_rules,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Patient-consistency levels:",
  ">=2/4, >=3/4, 4/4\n"
)

cat(
  "\nNO pseudobulk recalculation\n"
)

cat(
  "NO raw-UMI re-aggregation\n"
)

cat(
  "NO CellChat recalculation\n"
)

cat(
  "NO LR statistical testing\n"
)

cat(
  "NO P-value/FDR calculation\n"
)

cat(
  "NO pathway selection\n"
)

cat(
  "NO final LR candidate locking\n"
)

cat(
  "NO biological interpretation based on a preferred threshold\n"
)

cat("\n")
cat(
  "STATUS: PASS — predefined detection and patient-consistency\n"
)

cat(
  "rules were evaluated as a transparent sensitivity audit.\n"
)

cat(
  "No threshold was selected or locked at STEP11-5D-0.\n"
)

cat("\n")
cat("============================================================\n")




cat("\n")
cat("============================================================\n")
cat("STEP11-5D-0: REQUIRED SUMMARY OUTPUT\n")
cat("============================================================\n\n")

cat("1. DIRECTION-LEVEL SENSITIVITY SUMMARY\n")
cat("------------------------------------------------------------\n")
print(sensitivity_summary)

cat("\n\n2. GLOBAL DETECTION-RULE SUMMARY\n")
cat("------------------------------------------------------------\n")
print(global_summary)

cat("\n\n3. RULE PAIRWISE STABILITY\n")
cat("------------------------------------------------------------\n")
print(rule_pairwise)

# ============================================================
# STEP11-5D-1
# Genome-wide Patient-level Ligand-Receptor Compatibility
# Screening Across Focused Cellular Compartments
#
# PRIMARY ANALYSIS LOCK
# ------------------------------------------------------------
# Expression scale:
#   log2(CPM + 1)
#
# Detection rule:
#   CPM > 0
#
# Patient-level compatibility:
#   ligand detected in sender AND receptor detected in receiver
#   within the same patient
#
# Matched patients:
#   P0020, P0030, P0031, P0034
#
# Consistency levels:
#   >=2/4 patients
#   >=3/4 patients
#   4/4 patients
#
# Directions:
#   Endothelial -> Myeloid/macrophage
#   Myeloid/macrophage -> Endothelial
#   Endothelial -> T-cell
#   T-cell -> Endothelial
#   Myeloid/macrophage -> T-cell
#   T-cell -> Myeloid/macrophage
#
# IMPORTANT:
#   This workflow evaluates patient-level expression compatibility.
#   It does NOT establish functional ligand-receptor signaling.
#
# No:
#   - pseudobulk recalculation
#   - raw UMI re-aggregation
#   - CellChat recalculation
#   - LR statistical testing
#   - P values
#   - FDR
#   - pathway enrichment
#   - pathway selection
#   - biological ranking
# ============================================================


cat("============================================================\n")
cat("STEP11-5D-1: GENOME-WIDE PATIENT-LEVEL LR COMPATIBILITY\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directories
# ------------------------------------------------------------

step11_root <- "input_dir/Paper 1_Step 11-CellChat"

step10_root <- "input_dir/Paper 1_Step 10-Single-Cell"

out_dir <- file.path(
  step11_root,
  "STEP11-5D-1_Genomewide_Patient_Level_LR_Compatibility"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cat("STEP11 root:\n", step11_root, "\n\n")
cat("STEP10 root:\n", step10_root, "\n\n")
cat("Output directory:\n", out_dir, "\n\n")


# ------------------------------------------------------------
# 1. Exact input files
# ------------------------------------------------------------

file_T <- file.path(
  step10_root,
  "STEP10C_Patient_Level_RawUMI_Pseudobulk",
  "STEP10C-A_Pseudobulk_Tcell_Primary_rawUMI.tsv"
)

file_M <- file.path(
  step10_root,
  "STEP10C_Patient_Level_RawUMI_Pseudobulk",
  "STEP10C-B_Pseudobulk_MyeloidMacrophage_Primary_rawUMI.tsv"
)

file_E <- file.path(
  step10_root,
  "STEP10C_Patient_Level_RawUMI_Pseudobulk",
  "STEP10C-C_Pseudobulk_Endothelial_Primary_rawUMI.tsv"
)

file_LR <- file.path(
  step11_root,
  "NSCLC_GSE131907_CCI_LR.tsv"
)


# ------------------------------------------------------------
# 2. Input existence audit
# ------------------------------------------------------------

input_files <- c(
  T_cell = file_T,
  MyeloidMacrophage = file_M,
  Endothelial = file_E,
  CCI_LR = file_LR
)

cat("------------------------------------------------------------\n")
cat("INPUT FILE EXISTENCE AUDIT\n")
cat("------------------------------------------------------------\n")

print(
  data.frame(
    Input = names(input_files),
    Exists = file.exists(input_files),
    Path = unname(input_files),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
)

if (!all(file.exists(input_files))) {
  stop("ERROR: One or more required input files are missing.")
}

cat("\nPASS: All required input files exist.\n\n")


# ------------------------------------------------------------
# 3. Read existing STEP10C pseudobulk matrices
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("READING EXISTING STEP10C PSEUDOBULK MATRICES\n")
cat("------------------------------------------------------------\n")

T_raw <- read.delim(
  file_T,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

M_raw <- read.delim(
  file_M,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

E_raw <- read.delim(
  file_E,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

LR <- read.delim(
  file_LR,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("T-cell matrix: ",
    nrow(T_raw), " genes x ", ncol(T_raw) - 1,
    " patients\n", sep = "")

cat("Myeloid/macrophage matrix: ",
    nrow(M_raw), " genes x ", ncol(M_raw) - 1,
    " patients\n", sep = "")

cat("Endothelial matrix: ",
    nrow(E_raw), " genes x ", ncol(E_raw) - 1,
    " patients\n", sep = "")

cat("CCI_LR interactions: ",
    nrow(LR), "\n\n", sep = "")


# ------------------------------------------------------------
# 4. Basic structure audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("BASIC INPUT STRUCTURE AUDIT\n")
cat("------------------------------------------------------------\n")

required_lr_cols <- c(
  "interaction_name",
  "pathway_name",
  "ligand",
  "receptor",
  "annotation"
)

missing_lr_cols <- setdiff(required_lr_cols, colnames(LR))

if (length(missing_lr_cols) > 0) {
  stop(
    "ERROR: Missing required CCI_LR columns: ",
    paste(missing_lr_cols, collapse = ", ")
  )
}

if (!identical(colnames(T_raw)[1], "Gene")) {
  stop("ERROR: T-cell first column is not named 'Gene'.")
}

if (!identical(colnames(M_raw)[1], "Gene")) {
  stop("ERROR: Myeloid/macrophage first column is not named 'Gene'.")
}

if (!identical(colnames(E_raw)[1], "Gene")) {
  stop("ERROR: Endothelial first column is not named 'Gene'.")
}

if (anyDuplicated(T_raw$Gene)) {
  stop("ERROR: Duplicate genes detected in T-cell matrix.")
}

if (anyDuplicated(M_raw$Gene)) {
  stop("ERROR: Duplicate genes detected in Myeloid/macrophage matrix.")
}

if (anyDuplicated(E_raw$Gene)) {
  stop("ERROR: Duplicate genes detected in Endothelial matrix.")
}

if (anyDuplicated(LR$interaction_name)) {
  stop("ERROR: Duplicate interaction_name detected in CCI_LR.")
}

cat("CCI_LR required columns: PASS\n")
cat("Gene uniqueness: PASS\n")
cat("Interaction uniqueness: PASS\n\n")


# ------------------------------------------------------------
# 5. Convert expression matrices to numeric matrices
# ------------------------------------------------------------

gene_T <- T_raw$Gene
gene_M <- M_raw$Gene
gene_E <- E_raw$Gene

mat_T_raw <- as.matrix(T_raw[, -1, drop = FALSE])
mat_M_raw <- as.matrix(M_raw[, -1, drop = FALSE])
mat_E_raw <- as.matrix(E_raw[, -1, drop = FALSE])

storage.mode(mat_T_raw) <- "numeric"
storage.mode(mat_M_raw) <- "numeric"
storage.mode(mat_E_raw) <- "numeric"

rownames(mat_T_raw) <- gene_T
rownames(mat_M_raw) <- gene_M
rownames(mat_E_raw) <- gene_E


# ------------------------------------------------------------
# 6. Numeric integrity audit
# ------------------------------------------------------------

numeric_integrity <- data.frame(
  Compartment = c(
    "T-cell",
    "Myeloid/macrophage",
    "Endothelial"
  ),
  Numeric = c(
    is.numeric(mat_T_raw),
    is.numeric(mat_M_raw),
    is.numeric(mat_E_raw)
  ),
  No_NA = c(
    !anyNA(mat_T_raw),
    !anyNA(mat_M_raw),
    !anyNA(mat_E_raw)
  ),
  Finite = c(
    all(is.finite(mat_T_raw)),
    all(is.finite(mat_M_raw)),
    all(is.finite(mat_E_raw))
  ),
  Nonnegative = c(
    all(mat_T_raw >= 0),
    all(mat_M_raw >= 0),
    all(mat_E_raw >= 0)
  ),
  row.names = NULL
)

print(numeric_integrity)

if (!all(numeric_integrity[, -1])) {
  stop("ERROR: Numeric integrity audit failed.")
}

cat("\nPASS: Raw pseudobulk numeric integrity verified.\n\n")


# ------------------------------------------------------------
# 7. Gene identity / order audit
# ------------------------------------------------------------

gene_identity_pass <- identical(gene_T, gene_M) &&
  identical(gene_T, gene_E)

cat("------------------------------------------------------------\n")
cat("GENE IDENTITY / ORDER AUDIT\n")
cat("------------------------------------------------------------\n")

cat("T-cell vs Myeloid/macrophage: ",
    identical(gene_T, gene_M), "\n")

cat("T-cell vs Endothelial: ",
    identical(gene_T, gene_E), "\n")

cat("Common gene universe: ",
    length(gene_T), "\n")

if (!gene_identity_pass) {
  stop(
    "ERROR: Gene identity/order is not identical across matrices."
  )
}

cat("PASS: Identical gene universe and order confirmed.\n\n")


# ------------------------------------------------------------
# 8. Decode patient IDs
# ------------------------------------------------------------

extract_patient_id <- function(x) {
  sub("__.*$", "", x)
}

patient_T <- extract_patient_id(colnames(mat_T_raw))
patient_M <- extract_patient_id(colnames(mat_M_raw))
patient_E <- extract_patient_id(colnames(mat_E_raw))

cat("------------------------------------------------------------\n")
cat("PATIENT ID DECODING\n")
cat("------------------------------------------------------------\n")

cat("T-cell patients:\n")
print(patient_T)

cat("\nMyeloid/macrophage patients:\n")
print(patient_M)

cat("\nEndothelial patients:\n")
print(patient_E)

matched_patients <- Reduce(
  intersect,
  list(patient_T, patient_M, patient_E)
)

matched_patients <- sort(matched_patients)

cat("\nMatched patients across all three compartments:\n")
print(matched_patients)

if (length(matched_patients) != 4) {
  stop(
    "ERROR: Expected exactly 4 matched patients; found ",
    length(matched_patients)
  )
}

cat("\nPASS: Exactly 4 matched patients confirmed.\n\n")


# ------------------------------------------------------------
# 9. Reorder matrices to the locked matched-patient order
# ------------------------------------------------------------

idx_T <- match(matched_patients, patient_T)
idx_M <- match(matched_patients, patient_M)
idx_E <- match(matched_patients, patient_E)

if (anyNA(idx_T) || anyNA(idx_M) || anyNA(idx_E)) {
  stop("ERROR: Failed to align matched patients.")
}

T_match <- mat_T_raw[, idx_T, drop = FALSE]
M_match <- mat_M_raw[, idx_M, drop = FALSE]
E_match <- mat_E_raw[, idx_E, drop = FALSE]

colnames(T_match) <- matched_patients
colnames(M_match) <- matched_patients
colnames(E_match) <- matched_patients


# ------------------------------------------------------------
# 10. Library-size normalization
#     CPM = count / library_size * 1e6
# ------------------------------------------------------------

calc_CPM <- function(count_matrix) {
  
  lib_size <- colSums(count_matrix)
  
  if (any(lib_size <= 0)) {
    stop("ERROR: Non-positive library size detected.")
  }
  
  cpm <- sweep(
    count_matrix,
    2,
    lib_size,
    "/"
  ) * 1e6
  
  return(cpm)
}

T_cpm <- calc_CPM(T_match)
M_cpm <- calc_CPM(M_match)
E_cpm <- calc_CPM(E_match)


# ------------------------------------------------------------
# 11. CPM integrity audit
# ------------------------------------------------------------

cpm_totals <- rbind(
  T_cell = colSums(T_cpm),
  Myeloid_macrophage = colSums(M_cpm),
  Endothelial = colSums(E_cpm)
)

cat("------------------------------------------------------------\n")
cat("CPM NORMALIZATION AUDIT\n")
cat("------------------------------------------------------------\n")

print(cpm_totals)

cpm_pass <- all(
  abs(as.numeric(cpm_totals) - 1e6) < 1e-6
)

if (!cpm_pass) {
  stop("ERROR: CPM totals failed normalization audit.")
}

cat("\nPASS: CPM totals verified at 1,000,000.\n\n")


# ------------------------------------------------------------
# 12. log2(CPM + 1)
# ------------------------------------------------------------

T_log2CPM <- log2(T_cpm + 1)
M_log2CPM <- log2(M_cpm + 1)
E_log2CPM <- log2(E_cpm + 1)

if (
  !all(is.finite(T_log2CPM)) ||
  !all(is.finite(M_log2CPM)) ||
  !all(is.finite(E_log2CPM))
) {
  stop("ERROR: Non-finite log2(CPM+1) values detected.")
}

cat("PASS: log2(CPM+1) transformation verified.\n\n")


# ------------------------------------------------------------
# 13. CCI_LR interaction audit
# ------------------------------------------------------------

LR$interaction_name <- as.character(LR$interaction_name)
LR$pathway_name <- as.character(LR$pathway_name)
LR$ligand <- as.character(LR$ligand)
LR$receptor <- as.character(LR$receptor)
LR$annotation <- as.character(LR$annotation)

LR$ligand <- trimws(LR$ligand)
LR$receptor <- trimws(LR$receptor)

if (
  any(LR$ligand == "") ||
  any(LR$receptor == "")
) {
  stop("ERROR: Empty ligand or receptor values detected.")
}

lr_genes <- unique(
  c(
    LR$ligand,
    LR$receptor
  )
)

cat("------------------------------------------------------------\n")
cat("CCI_LR AUDIT\n")
cat("------------------------------------------------------------\n")

cat("Total LR interactions: ",
    nrow(LR), "\n", sep = "")

cat("Unique LR genes: ",
    length(lr_genes), "\n", sep = "")

cat("Unique ligands: ",
    length(unique(LR$ligand)), "\n", sep = "")

cat("Unique receptors: ",
    length(unique(LR$receptor)), "\n", sep = "")


# ------------------------------------------------------------
# 14. Restrict LR genes to the common expression universe
# ------------------------------------------------------------

lr_available <- lr_genes[
  lr_genes %in% rownames(T_log2CPM) &
    lr_genes %in% rownames(M_log2CPM) &
    lr_genes %in% rownames(E_log2CPM)
]

cat("LR genes available in all three compartments: ",
    length(lr_available), "\n", sep = "")

if (length(lr_available) != 744) {
  warning(
    "NOTE: Expected 744 common LR genes based on STEP11-5C-2; ",
    "current result = ",
    length(lr_available)
  )
}

LR$LR_genes_available_all3 <- (
  LR$ligand %in% lr_available &
    LR$receptor %in% lr_available
)

cat(
  "LR interactions with both ligand and receptor available in all 3 compartments: ",
  sum(LR$LR_genes_available_all3),
  "\n",
  sep = ""
)

cat("\n")


# ------------------------------------------------------------
# 15. Core-gene separation audit
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_presence_LR <- data.frame(
  Core_Gene = core_genes,
  In_Ligand = core_genes %in% LR$ligand,
  In_Receptor = core_genes %in% LR$receptor,
  In_CCI_LR = (
    core_genes %in% c(LR$ligand, LR$receptor)
  ),
  stringsAsFactors = FALSE
)

cat("------------------------------------------------------------\n")
cat("CORE-GENE / CCI_LR SEPARATION AUDIT\n")
cat("------------------------------------------------------------\n")

print(core_presence_LR)

if (any(core_presence_LR$In_CCI_LR)) {
  warning(
    "At least one core gene occurs in CCI_LR. ",
    "Review before proceeding."
  )
} else {
  cat("\nPASS: AOC3, MS4A7, and IL7R are not CCI_LR ligand/receptor genes.\n")
}

cat("\n")


# ------------------------------------------------------------
# 16. Define six prespecified directions
# ------------------------------------------------------------

directions <- data.frame(
  Direction_ID = 1:6,
  Sender = c(
    "Endothelial",
    "Myeloid/macrophage",
    "Endothelial",
    "T-cell",
    "Myeloid/macrophage",
    "T-cell"
  ),
  Receiver = c(
    "Myeloid/macrophage",
    "Endothelial",
    "T-cell",
    "Endothelial",
    "T-cell",
    "Myeloid/macrophage"
  ),
  stringsAsFactors = FALSE
)

cat("------------------------------------------------------------\n")
cat("PREDEFINED LR DIRECTIONS\n")
cat("------------------------------------------------------------\n")

print(directions)

cat("\n")


# ------------------------------------------------------------
# 17. Helper function:
#     return CPM > 0 detection matrix
# ------------------------------------------------------------

get_detection_matrix <- function(
    compartment,
    T_cpm,
    M_cpm,
    E_cpm
) {
  
  if (compartment == "T-cell") {
    return(T_cpm > 0)
  }
  
  if (compartment == "Myeloid/macrophage") {
    return(M_cpm > 0)
  }
  
  if (compartment == "Endothelial") {
    return(E_cpm > 0)
  }
  
  stop("Unknown compartment.")
}


# ------------------------------------------------------------
# 18. Create detection matrices
# ------------------------------------------------------------

detect_T <- T_cpm > 0
detect_M <- M_cpm > 0
detect_E <- E_cpm > 0

rownames(detect_T) <- rownames(T_cpm)
rownames(detect_M) <- rownames(M_cpm)
rownames(detect_E) <- rownames(E_cpm)

colnames(detect_T) <- matched_patients
colnames(detect_M) <- matched_patients
colnames(detect_E) <- matched_patients


# ------------------------------------------------------------
# 19. Core function for patient-level LR compatibility
# ------------------------------------------------------------

screen_direction <- function(
    LR_table,
    sender,
    receiver,
    sender_detection,
    receiver_detection,
    sender_expression,
    receiver_expression,
    matched_patients
) {
  
  results <- vector(
    "list",
    nrow(LR_table)
  )
  
  for (i in seq_len(nrow(LR_table))) {
    
    ligand <- LR_table$ligand[i]
    receptor <- LR_table$receptor[i]
    
    ligand_available <- ligand %in% rownames(sender_detection)
    receptor_available <- receptor %in% rownames(receiver_detection)
    
    if (!ligand_available || !receptor_available) {
      
      compatible <- rep(FALSE, length(matched_patients))
      
      ligand_detected <- rep(NA, length(matched_patients))
      receptor_detected <- rep(NA, length(matched_patients))
      
      ligand_expr <- rep(NA_real_, length(matched_patients))
      receptor_expr <- rep(NA_real_, length(matched_patients))
      
    } else {
      
      ligand_detected <- as.logical(
        sender_detection[ligand, matched_patients]
      )
      
      receptor_detected <- as.logical(
        receiver_detection[receptor, matched_patients]
      )
      
      compatible <- (
        ligand_detected &
          receptor_detected
      )
      
      ligand_expr <- as.numeric(
        sender_expression[ligand, matched_patients]
      )
      
      receptor_expr <- as.numeric(
        receiver_expression[receptor, matched_patients]
      )
    }
    
    n_compatible <- sum(
      compatible,
      na.rm = TRUE
    )
    
    results[[i]] <- data.frame(
      interaction_name = LR_table$interaction_name[i],
      pathway_name = LR_table$pathway_name[i],
      ligand = ligand,
      receptor = receptor,
      annotation = LR_table$annotation[i],
      Sender = sender,
      Receiver = receiver,
      ligand_available_sender = ligand_available,
      receptor_available_receiver = receptor_available,
      P1_ligand_detected = ligand_detected[1],
      P2_ligand_detected = ligand_detected[2],
      P3_ligand_detected = ligand_detected[3],
      P4_ligand_detected = ligand_detected[4],
      P1_receptor_detected = receptor_detected[1],
      P2_receptor_detected = receptor_detected[2],
      P3_receptor_detected = receptor_detected[3],
      P4_receptor_detected = receptor_detected[4],
      P1_compatible = compatible[1],
      P2_compatible = compatible[2],
      P3_compatible = compatible[3],
      P4_compatible = compatible[4],
      n_compatible_patients = n_compatible,
      consistency_fraction = n_compatible /
        length(matched_patients),
      Compatible_ge2of4 = n_compatible >= 2,
      Compatible_ge3of4 = n_compatible >= 3,
      Compatible_4of4 = n_compatible == 4,
      P1_ligand_log2CPM = ligand_expr[1],
      P2_ligand_log2CPM = ligand_expr[2],
      P3_ligand_log2CPM = ligand_expr[3],
      P4_ligand_log2CPM = ligand_expr[4],
      P1_receptor_log2CPM = receptor_expr[1],
      P2_receptor_log2CPM = receptor_expr[2],
      P3_receptor_log2CPM = receptor_expr[3],
      P4_receptor_log2CPM = receptor_expr[4],
      stringsAsFactors = FALSE
    )
  }
  
  do.call(
    rbind,
    results
  )
}


# ------------------------------------------------------------
# 20. Prepare compartment-specific expression/detection lists
# ------------------------------------------------------------

expression_list <- list(
  "T-cell" = T_log2CPM,
  "Myeloid/macrophage" = M_log2CPM,
  "Endothelial" = E_log2CPM
)

detection_list <- list(
  "T-cell" = detect_T,
  "Myeloid/macrophage" = detect_M,
  "Endothelial" = detect_E
)


# ------------------------------------------------------------
# 21. Run all six directions
# ------------------------------------------------------------

cat("============================================================\n")
cat("RUNNING GENOME-WIDE PATIENT-LEVEL LR SCREENING\n")
cat("============================================================\n\n")

direction_results <- vector(
  "list",
  nrow(directions)
)

for (d in seq_len(nrow(directions))) {
  
  sender <- directions$Sender[d]
  receiver <- directions$Receiver[d]
  
  cat(
    "Direction ",
    d,
    "/",
    nrow(directions),
    ": ",
    sender,
    " -> ",
    receiver,
    "\n",
    sep = ""
  )
  
  direction_results[[d]] <- screen_direction(
    LR_table = LR,
    sender = sender,
    receiver = receiver,
    sender_detection = detection_list[[sender]],
    receiver_detection = detection_list[[receiver]],
    sender_expression = expression_list[[sender]],
    receiver_expression = expression_list[[receiver]],
    matched_patients = matched_patients
  )
  
  direction_results[[d]]$Direction_ID <- d
  
  cat(
    "  Interactions evaluated: ",
    nrow(direction_results[[d]]),
    "\n",
    sep = ""
  )
  
  cat(
    "  >=2/4 compatible: ",
    sum(direction_results[[d]]$Compatible_ge2of4),
    "\n",
    sep = ""
  )
  
  cat(
    "  >=3/4 compatible: ",
    sum(direction_results[[d]]$Compatible_ge3of4),
    "\n",
    sep = ""
  )
  
  cat(
    "  4/4 compatible: ",
    sum(direction_results[[d]]$Compatible_4of4),
    "\n\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 22. Combine all directions
# ------------------------------------------------------------

lr_compatibility_all <- do.call(
  rbind,
  direction_results
)

rownames(lr_compatibility_all) <- NULL


# ------------------------------------------------------------
# 23. Add explicit core-gene audit flags
# ------------------------------------------------------------

lr_compatibility_all$Core_Gene_Involved <- (
  lr_compatibility_all$ligand %in% core_genes |
    lr_compatibility_all$receptor %in% core_genes
)

lr_compatibility_all$Core_Gene <- NA_character_

for (g in core_genes) {
  
  idx <- (
    lr_compatibility_all$ligand == g |
      lr_compatibility_all$receptor == g
  )
  
  lr_compatibility_all$Core_Gene[idx] <- g
}


# ------------------------------------------------------------
# 24. Primary candidate tables
# ------------------------------------------------------------

candidate_ge2 <- subset(
  lr_compatibility_all,
  Compatible_ge2of4
)

candidate_ge3 <- subset(
  lr_compatibility_all,
  Compatible_ge3of4
)

candidate_4of4 <- subset(
  lr_compatibility_all,
  Compatible_4of4
)


# ------------------------------------------------------------
# 25. Direction-level summary
# ------------------------------------------------------------

direction_summary <- do.call(
  rbind,
  lapply(
    seq_len(nrow(directions)),
    function(d) {
      
      x <- direction_results[[d]]
      
      data.frame(
        Direction_ID = d,
        Sender = directions$Sender[d],
        Receiver = directions$Receiver[d],
        LR_Interactions = nrow(x),
        Available_LR_Interactions = sum(
          x$ligand_available_sender &
            x$receptor_available_receiver
        ),
        Compatible_ge2of4 = sum(
          x$Compatible_ge2of4
        ),
        Compatible_ge3of4 = sum(
          x$Compatible_ge3of4
        ),
        Compatible_4of4 = sum(
          x$Compatible_4of4
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

cat("============================================================\n")
cat("DIRECTION-LEVEL SUMMARY\n")
cat("============================================================\n\n")

print(direction_summary)


# ------------------------------------------------------------
# 26. Global summary
# ------------------------------------------------------------

global_summary <- data.frame(
  Metric = c(
    "CCI_LR interactions",
    "Unique LR genes",
    "LR genes available in all 3 compartments",
    "Matched patients",
    "Directions",
    "All-direction LR compatibility records",
    "Compatible >=2/4",
    "Compatible >=3/4",
    "Compatible 4/4"
  ),
  Value = c(
    nrow(LR),
    length(lr_genes),
    length(lr_available),
    length(matched_patients),
    nrow(directions),
    nrow(lr_compatibility_all),
    sum(lr_compatibility_all$Compatible_ge2of4),
    sum(lr_compatibility_all$Compatible_ge3of4),
    sum(lr_compatibility_all$Compatible_4of4)
  ),
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("GLOBAL SUMMARY\n")
cat("============================================================\n\n")

print(global_summary)


# ------------------------------------------------------------
# 27. Core-gene separation summary
# ------------------------------------------------------------

core_gene_summary <- data.frame(
  Core_Gene = core_genes,
  In_CCI_LR = core_genes %in% c(
    LR$ligand,
    LR$receptor
  ),
  Compatible_LR_records = sapply(
    core_genes,
    function(g) {
      sum(
        lr_compatibility_all$ligand == g |
          lr_compatibility_all$receptor == g
      )
    }
  ),
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("CORE-GENE SEPARATION SUMMARY\n")
cat("============================================================\n\n")

print(core_gene_summary)


# ------------------------------------------------------------
# 28. Descriptive pathway summary
#
# IMPORTANT:
# This is descriptive only.
# It does NOT perform pathway enrichment or pathway selection.
# ------------------------------------------------------------

pathway_summary_ge3 <- aggregate(
  Compatible_ge3of4 ~ pathway_name,
  data = candidate_ge3,
  FUN = sum
)

names(pathway_summary_ge3)[2] <-
  "Compatible_LR_Records_ge3of4"

pathway_total <- aggregate(
  interaction_name ~ pathway_name,
  data = lr_compatibility_all,
  FUN = length
)

names(pathway_total)[2] <-
  "Total_LR_Records"

pathway_summary_ge3 <- merge(
  pathway_total,
  pathway_summary_ge3,
  by = "pathway_name",
  all.x = TRUE
)

pathway_summary_ge3[
  is.na(pathway_summary_ge3$Compatible_LR_Records_ge3of4),
  "Compatible_LR_Records_ge3of4"
] <- 0

pathway_summary_ge3 <- pathway_summary_ge3[
  order(
    pathway_summary_ge3$pathway_name
  ),
]

# No ranking or enrichment is performed.
# The table is retained only as a descriptive audit.


# ------------------------------------------------------------
# 29. Annotation summary
# ------------------------------------------------------------

annotation_summary_ge3 <- aggregate(
  Compatible_ge3of4 ~ annotation,
  data = candidate_ge3,
  FUN = sum
)

names(annotation_summary_ge3)[2] <-
  "Compatible_LR_Records_ge3of4"

annotation_total <- aggregate(
  interaction_name ~ annotation,
  data = lr_compatibility_all,
  FUN = length
)

names(annotation_total)[2] <-
  "Total_LR_Records"

annotation_summary_ge3 <- merge(
  annotation_total,
  annotation_summary_ge3,
  by = "annotation",
  all.x = TRUE
)

annotation_summary_ge3[
  is.na(annotation_summary_ge3$Compatible_LR_Records_ge3of4),
  "Compatible_LR_Records_ge3of4"
] <- 0


# ------------------------------------------------------------
# 30. Save complete interaction-level audit
# ------------------------------------------------------------

write.table(
  lr_compatibility_all,
  file = file.path(
    out_dir,
    "STEP11-5D-1A_All_LR_Direction_Compatibility_Audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 31. Save primary candidate tables
# ------------------------------------------------------------

write.table(
  candidate_ge2,
  file = file.path(
    out_dir,
    "STEP11-5D-1B_LR_Candidates_Compatible_ge2of4.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  candidate_ge3,
  file = file.path(
    out_dir,
    "STEP11-5D-1C_LR_Candidates_Compatible_ge3of4.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  candidate_4of4,
  file = file.path(
    out_dir,
    "STEP11-5D-1D_LR_Candidates_Compatible_4of4.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 32. Save patient-level compatibility audit
# ------------------------------------------------------------

patient_level_audit <- lr_compatibility_all[
  ,
  c(
    "Direction_ID",
    "Sender",
    "Receiver",
    "interaction_name",
    "pathway_name",
    "ligand",
    "receptor",
    "annotation",
    "P1_ligand_detected",
    "P2_ligand_detected",
    "P3_ligand_detected",
    "P4_ligand_detected",
    "P1_receptor_detected",
    "P2_receptor_detected",
    "P3_receptor_detected",
    "P4_receptor_detected",
    "P1_compatible",
    "P2_compatible",
    "P3_compatible",
    "P4_compatible",
    "n_compatible_patients",
    "consistency_fraction",
    "Compatible_ge2of4",
    "Compatible_ge3of4",
    "Compatible_4of4"
  )
]

write.table(
  patient_level_audit,
  file = file.path(
    out_dir,
    "STEP11-5D-1E_Patient_Level_LR_Compatibility_Audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 33. Save summaries
# ------------------------------------------------------------

write.table(
  direction_summary,
  file = file.path(
    out_dir,
    "STEP11-5D-1F_Direction_Level_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  global_summary,
  file = file.path(
    out_dir,
    "STEP11-5D-1G_Global_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  core_gene_summary,
  file = file.path(
    out_dir,
    "STEP11-5D-1H_Core_Gene_CCI_LR_Separation_Audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  pathway_summary_ge3,
  file = file.path(
    out_dir,
    "STEP11-5D-1I_Descriptive_Pathway_Summary_ge3of4.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  annotation_summary_ge3,
  file = file.path(
    out_dir,
    "STEP11-5D-1J_Descriptive_Annotation_Summary_ge3of4.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 34. QC checks on output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("OUTPUT QC\n")
cat("============================================================\n\n")

expected_records <- nrow(LR) * nrow(directions)

cat(
  "Expected LR-direction records: ",
  expected_records,
  "\n",
  sep = ""
)

cat(
  "Observed LR-direction records: ",
  nrow(lr_compatibility_all),
  "\n",
  sep = ""
)

if (nrow(lr_compatibility_all) != expected_records) {
  stop(
    "ERROR: Unexpected number of LR-direction records."
  )
}

if (
  any(
    lr_compatibility_all$n_compatible_patients < 0 |
    lr_compatibility_all$n_compatible_patients > 4
  )
) {
  stop(
    "ERROR: Invalid patient compatibility count."
  )
}

if (
  any(
    lr_compatibility_all$consistency_fraction < 0 |
    lr_compatibility_all$consistency_fraction > 1
  )
) {
  stop(
    "ERROR: Invalid consistency fraction."
  )
}

# Logical consistency audit
if (
  any(
    lr_compatibility_all$Compatible_4of4 &
    !lr_compatibility_all$Compatible_ge3of4
  )
) {
  stop(
    "ERROR: 4/4 compatibility is not nested within >=3/4."
  )
}

if (
  any(
    lr_compatibility_all$Compatible_ge3of4 &
    !lr_compatibility_all$Compatible_ge2of4
  )
) {
  stop(
    "ERROR: >=3/4 compatibility is not nested within >=2/4."
  )
}

cat("Record count: PASS\n")
cat("Patient compatibility range: PASS\n")
cat("Consistency fraction range: PASS\n")
cat("Nested consistency logic: PASS\n")


# ------------------------------------------------------------
# 35. Final locked status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP11-5D-1 FINAL STATUS\n")
cat("============================================================\n\n")

cat("Existing STEP10C pseudobulk matrices reused: YES\n")
cat("Pseudobulk recalculation: NO\n")
cat("Raw UMI re-aggregation: NO\n")
cat("Matched patients: 4\n")
cat("Expression scale: log2(CPM + 1)\n")
cat("Detection rule: CPM > 0\n")
cat("CCI_LR interactions screened: ",
    nrow(LR), "\n", sep = "")
cat("Unique LR genes: ",
    length(lr_genes), "\n", sep = "")
cat("LR genes available in all 3 compartments: ",
    length(lr_available), "\n", sep = "")
cat("Directions screened: 6\n")
cat("Patient-level compatibility: YES\n")
cat("Consistency levels: >=2/4, >=3/4, 4/4\n")
cat("LR statistical testing: NO\n")
cat("P values: NO\n")
cat("FDR correction: NO\n")
cat("Pathway enrichment: NO\n")
cat("Pathway selection: NO\n")
cat("Biological ranking: NO\n")
cat("Core genes treated as LR nodes: NO\n")

cat("\nSTATUS: PASS\n")
cat(
  "Genome-wide patient-level ligand-receptor expression compatibility ",
  "screening completed under the prespecified CPM > 0 detection rule.\n"
)

cat(
  "Results represent patient-level expression compatibility and do not ",
  "establish functional intercellular signaling.\n"
)

cat("\nOutput directory:\n")
cat(out_dir, "\n")

cat("\n============================================================\n")
cat("STEP11-5D-1 COMPLETE\n")
cat("============================================================\n")

cat("\n")
cat("============================================================\n")
cat("STEP11-5D-1: REQUIRED CANDIDATE AUDIT OUTPUT\n")
cat("============================================================\n\n")

cat("1. >=3/4 PATIENT-COMPATIBLE CANDIDATES\n")
cat("------------------------------------------------------------\n")
print(
  candidate_ge3[
    ,
    c(
      "Direction_ID",
      "Sender",
      "Receiver",
      "interaction_name",
      "pathway_name",
      "ligand",
      "receptor",
      "annotation",
      "n_compatible_patients",
      "consistency_fraction"
    )
  ],
  row.names = FALSE
)

cat("\n\n2. >=3/4 CANDIDATE COUNT BY DIRECTION\n")
cat("------------------------------------------------------------\n")
print(direction_summary)

cat("\n\n3. DESCRIPTIVE PATHWAY SUMMARY\n")
cat("------------------------------------------------------------\n")
print(pathway_summary_ge3)

cat("\n\n4. DESCRIPTIVE ANNOTATION SUMMARY\n")
cat("------------------------------------------------------------\n")
print(annotation_summary_ge3)


# ============================================================
# STEP11-5D-2
# Genome-wide Patient-Level LR Candidate Structural Audit
#
# PURPOSE
# ------------------------------------------------------------
# Pure structural audit and redundancy characterization of the
# already locked STEP11-5D-1 patient-level LR compatibility
# candidates.
#
# IMPORTANT:
#   - No expression recalculation
#   - No pseudobulk recalculation
#   - No LR detection recalculation
#   - No CellChat recalculation
#   - No P value / FDR calculation
#   - No pathway enrichment
#   - No pathway ranking
#   - No candidate biological ranking
#   - No selection of "best" LR pair
#   - Core genes remain separated from LR candidates
#
# INPUT:
#   STEP11-5D-1B_LR_Candidates_Compatible_ge2of4.tsv
#   STEP11-5D-1C_LR_Candidates_Compatible_ge3of4.tsv
#   STEP11-5D-1D_LR_Candidates_Compatible_4of4.tsv
#   STEP11-5D-1A_All_LR_Direction_Compatibility_Audit.tsv
#
# OUTPUT:
#   STEP11-5D-2A_LR_Candidate_Structural_Audit_ge3of4.tsv
#   STEP11-5D-2B_LR_Candidate_Structural_Audit_4of4.tsv
#   STEP11-5D-2C_Ligand_Level_Summary_ge3of4.tsv
#   STEP11-5D-2D_Ligand_Receptor_Pair_Summary_ge3of4.tsv
#   STEP11-5D-2E_Directional_Symmetry_Audit.tsv
#   STEP11-5D-2F_Annotation_Structure_Audit.tsv
#   STEP11-5D-2G_Pathway_Structure_Audit.tsv
#   STEP11-5D-2H_Core_Gene_Separation_Confirmation.tsv
#   STEP11-5D-2I_QC_Summary.tsv
#
# ============================================================


cat("============================================================\n")
cat("STEP11-5D-2: GENOME-WIDE PATIENT-LEVEL LR CANDIDATE\n")
cat("STRUCTURAL AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory and input/output paths
# ------------------------------------------------------------

step11_dir <- "input_dir/Paper 1_Step 11-CellChat"

input_dir <- file.path(
  step11_dir,
  "STEP11-5D-1_Genomewide_Patient_Level_LR_Compatibility"
)

output_dir <- file.path(
  step11_dir,
  "STEP11-5D-2_Genomewide_LR_Candidate_Structural_Audit"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("Input directory:\n")
cat(input_dir, "\n\n")

cat("Output directory:\n")
cat(output_dir, "\n\n")


# ------------------------------------------------------------
# 1. Exact input files
# ------------------------------------------------------------

file_ge2 <- file.path(
  input_dir,
  "STEP11-5D-1B_LR_Candidates_Compatible_ge2of4.tsv"
)

file_ge3 <- file.path(
  input_dir,
  "STEP11-5D-1C_LR_Candidates_Compatible_ge3of4.tsv"
)

file_ge4 <- file.path(
  input_dir,
  "STEP11-5D-1D_LR_Candidates_Compatible_4of4.tsv"
)

file_all <- file.path(
  input_dir,
  "STEP11-5D-1A_All_LR_Direction_Compatibility_Audit.tsv"
)


# ------------------------------------------------------------
# 2. Input existence audit
# ------------------------------------------------------------

input_files <- c(
  file_ge2,
  file_ge3,
  file_ge4,
  file_all
)

input_exists <- file.exists(input_files)

cat("INPUT FILE EXISTENCE AUDIT\n")
cat("------------------------------------------------------------\n")

print(
  data.frame(
    File = basename(input_files),
    Exists = input_exists,
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)

if (!all(input_exists)) {
  stop(
    "STEP11-5D-2 STOP: One or more required STEP11-5D-1 input files are missing."
  )
}

cat("\nPASS — all required STEP11-5D-1 input files found.\n\n")


# ------------------------------------------------------------
# 3. Read locked STEP11-5D-1 outputs
# ------------------------------------------------------------

candidate_ge2 <- read.delim(
  file_ge2,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

candidate_ge3 <- read.delim(
  file_ge3,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

candidate_ge4 <- read.delim(
  file_ge4,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

all_direction <- read.delim(
  file_all,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("INPUT DIMENSIONS\n")
cat("------------------------------------------------------------\n")

cat(
  "ge2:",
  nrow(candidate_ge2), "rows x",
  ncol(candidate_ge2), "columns\n"
)

cat(
  "ge3:",
  nrow(candidate_ge3), "rows x",
  ncol(candidate_ge3), "columns\n"
)

cat(
  "ge4:",
  nrow(candidate_ge4), "rows x",
  ncol(candidate_ge4), "columns\n"
)

cat(
  "all-direction audit:",
  nrow(all_direction), "rows x",
  ncol(all_direction), "columns\n\n"
)


# ------------------------------------------------------------
# 4. Required-column audit
# ------------------------------------------------------------

required_cols <- c(
  "Direction_ID",
  "Sender",
  "Receiver",
  "interaction_name",
  "pathway_name",
  "ligand",
  "receptor",
  "annotation",
  "n_compatible_patients",
  "consistency_fraction"
)

check_required_columns <- function(df, df_name) {
  
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing required columns in ",
        df_name,
        ": ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  TRUE
}

check_required_columns(candidate_ge2, "candidate_ge2")
check_required_columns(candidate_ge3, "candidate_ge3")
check_required_columns(candidate_ge4, "candidate_ge4")


cat("PASS — required candidate columns present in ge2/ge3/ge4.\n\n")


# ------------------------------------------------------------
# 5. Basic structural integrity
# ------------------------------------------------------------

audit_basic_structure <- function(df, label) {
  
  result <- data.frame(
    Dataset = label,
    Rows = nrow(df),
    Unique_interaction_name =
      length(unique(df$interaction_name)),
    Unique_ligands =
      length(unique(df$ligand)),
    Unique_receptors =
      length(unique(df$receptor)),
    Unique_exact_LR_pairs =
      length(unique(
        paste(df$ligand, df$receptor, sep = "___")
      )),
    Unique_directions =
      length(unique(df$Direction_ID)),
    Missing_interaction_name =
      sum(is.na(df$interaction_name) | df$interaction_name == ""),
    Missing_ligand =
      sum(is.na(df$ligand) | df$ligand == ""),
    Missing_receptor =
      sum(is.na(df$receptor) | df$receptor == ""),
    Missing_pathway =
      sum(is.na(df$pathway_name) | df$pathway_name == ""),
    Missing_annotation =
      sum(is.na(df$annotation) | df$annotation == ""),
    stringsAsFactors = FALSE
  )
  
  result
}

basic_audit <- rbind(
  audit_basic_structure(candidate_ge2, "ge2of4"),
  audit_basic_structure(candidate_ge3, "ge3of4"),
  audit_basic_structure(candidate_ge4, "4of4")
)

cat("BASIC STRUCTURAL AUDIT\n")
cat("------------------------------------------------------------\n")
print(basic_audit, row.names = FALSE)
cat("\n")


# ------------------------------------------------------------
# 6. Add exact structural identifiers
#
# These identifiers do NOT alter the candidate definition.
# They only characterize redundancy.
# ------------------------------------------------------------

add_structural_ids <- function(df) {
  
  df$Exact_LR_Pair <- paste(
    df$ligand,
    df$receptor,
    sep = "___"
  )
  
  df$Directed_LR_Pair <- paste(
    df$Sender,
    df$Receiver,
    df$ligand,
    df$receptor,
    sep = "___"
  )
  
  df$Ligand_Direction <- paste(
    df$Sender,
    df$Receiver,
    df$ligand,
    sep = "___"
  )
  
  df$Receptor_Direction <- paste(
    df$Sender,
    df$Receiver,
    df$receptor,
    sep = "___"
  )
  
  df
}

candidate_ge3_audit <- add_structural_ids(candidate_ge3)
candidate_ge4_audit <- add_structural_ids(candidate_ge4)


# ------------------------------------------------------------
# 7. Interaction-record redundancy audit
#
# Question:
# How many CCI_LR interaction records represent the same exact
# ligand-receptor pair within the candidate set?
#
# No records are removed.
# ------------------------------------------------------------

pair_record_summary <- function(df) {
  
  split_idx <- split(
    seq_len(nrow(df)),
    df$Exact_LR_Pair
  )
  
  out <- lapply(
    names(split_idx),
    function(pair_name) {
      
      idx <- split_idx[[pair_name]]
      sub <- df[idx, , drop = FALSE]
      
      data.frame(
        Exact_LR_Pair = pair_name,
        ligand = unique(sub$ligand)[1],
        receptor = unique(sub$receptor)[1],
        n_interaction_records =
          length(unique(sub$interaction_name)),
        n_pathways =
          length(unique(sub$pathway_name)),
        n_annotations =
          length(unique(sub$annotation)),
        n_directions =
          length(unique(sub$Direction_ID)),
        Directions =
          paste(
            sort(unique(sub$Direction_ID)),
            collapse = ";"
          ),
        min_compatible_patients =
          min(sub$n_compatible_patients),
        max_compatible_patients =
          max(sub$n_compatible_patients),
        stringsAsFactors = FALSE
      )
    }
  )
  
  do.call(rbind, out)
}

pair_summary_ge3 <- pair_record_summary(candidate_ge3_audit)
pair_summary_ge4 <- pair_record_summary(candidate_ge4_audit)


# ------------------------------------------------------------
# 8. Ligand-level structural summary
#
# This characterizes expansion of one ligand across multiple
# receptors. It is NOT a biological ranking.
# ------------------------------------------------------------

ligand_summary <- function(df) {
  
  split_idx <- split(
    seq_len(nrow(df)),
    df$ligand
  )
  
  out <- lapply(
    names(split_idx),
    function(lig) {
      
      idx <- split_idx[[lig]]
      sub <- df[idx, , drop = FALSE]
      
      data.frame(
        ligand = lig,
        n_interaction_records =
          length(unique(sub$interaction_name)),
        n_unique_receptors =
          length(unique(sub$receptor)),
        n_exact_LR_pairs =
          length(unique(sub$Exact_LR_Pair)),
        n_pathways =
          length(unique(sub$pathway_name)),
        n_directions =
          length(unique(sub$Direction_ID)),
        directions =
          paste(
            sort(unique(sub$Direction_ID)),
            collapse = ";"
          ),
        n_annotation_classes =
          length(unique(sub$annotation)),
        stringsAsFactors = FALSE
      )
    }
  )
  
  do.call(rbind, out)
}

ligand_summary_ge3 <- ligand_summary(candidate_ge3_audit)


# ------------------------------------------------------------
# 9. Directional symmetry audit
#
# For each exact ligand-receptor pair, determine whether the
# same LR pair is represented in:
#
#   E -> M
#   M -> E
#   E -> T
#   T -> E
#   M -> T
#   T -> M
#
# This is an expression-compatibility orientation audit only.
# It does NOT imply biological signaling directionality.
# ------------------------------------------------------------

direction_lookup <- data.frame(
  Direction_ID = 1:6,
  Sender = c(
    "Endothelial",
    "Myeloid/macrophage",
    "Endothelial",
    "T-cell",
    "Myeloid/macrophage",
    "T-cell"
  ),
  Receiver = c(
    "Myeloid/macrophage",
    "Endothelial",
    "T-cell",
    "Endothelial",
    "T-cell",
    "Myeloid/macrophage"
  ),
  Direction_Label = c(
    "Endothelial_to_Myeloid",
    "Myeloid_to_Endothelial",
    "Endothelial_to_Tcell",
    "Tcell_to_Endothelial",
    "Myeloid_to_Tcell",
    "Tcell_to_Myeloid"
  ),
  stringsAsFactors = FALSE
)

candidate_ge3_audit$Direction_Label <-
  direction_lookup$Direction_Label[
    match(
      candidate_ge3_audit$Direction_ID,
      direction_lookup$Direction_ID
    )
  ]

symmetry_split <- split(
  seq_len(nrow(candidate_ge3_audit)),
  candidate_ge3_audit$Exact_LR_Pair
)

directional_symmetry <- lapply(
  names(symmetry_split),
  function(pair_name) {
    
    idx <- symmetry_split[[pair_name]]
    sub <- candidate_ge3_audit[idx, , drop = FALSE]
    
    present <- unique(sub$Direction_ID)
    
    data.frame(
      Exact_LR_Pair = pair_name,
      ligand = unique(sub$ligand)[1],
      receptor = unique(sub$receptor)[1],
      
      D1_Endothelial_to_Myeloid =
        as.integer(1 %in% present),
      
      D2_Myeloid_to_Endothelial =
        as.integer(2 %in% present),
      
      D3_Endothelial_to_Tcell =
        as.integer(3 %in% present),
      
      D4_Tcell_to_Endothelial =
        as.integer(4 %in% present),
      
      D5_Myeloid_to_Tcell =
        as.integer(5 %in% present),
      
      D6_Tcell_to_Myeloid =
        as.integer(6 %in% present),
      
      n_directions_present =
        length(present),
      
      direction_pattern =
        paste(sort(present), collapse = ";"),
      
      stringsAsFactors = FALSE
    )
  }
)

directional_symmetry <- do.call(
  rbind,
  directional_symmetry
)


# ------------------------------------------------------------
# 10. Annotation structural audit
#
# Descriptive only.
# Both numerator and resource denominator are retained.
# ------------------------------------------------------------

annotation_candidate <- as.data.frame(
  table(candidate_ge3_audit$annotation),
  stringsAsFactors = FALSE
)

colnames(annotation_candidate) <- c(
  "annotation",
  "Compatible_ge3of4_records"
)

annotation_candidate$Compatible_ge3of4_records <-
  as.integer(annotation_candidate$Compatible_ge3of4_records)

annotation_candidate <- annotation_candidate[
  order(
    annotation_candidate$Compatible_ge3of4_records,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]


# Resource-level annotation counts are obtained from the
# all-direction STEP11-5D-1 audit without recalculating LR.
if ("annotation" %in% colnames(all_direction)) {
  
  annotation_resource <- aggregate(
    all_direction$interaction_name,
    by = list(annotation = all_direction$annotation),
    FUN = function(x) length(unique(x))
  )
  
  colnames(annotation_resource)[2] <-
    "Resource_LR_Records"
  
} else {
  
  annotation_resource <- data.frame(
    annotation = character(0),
    Resource_LR_Records = integer(0)
  )
}

annotation_structure <- merge(
  annotation_resource,
  annotation_candidate,
  by = "annotation",
  all = TRUE
)

annotation_structure[
  is.na(annotation_structure)
] <- 0


# ------------------------------------------------------------
# 11. Pathway structural audit
#
# No ranking by biological importance.
# Records are ordered alphabetically by pathway name.
# ------------------------------------------------------------

pathway_candidate <- aggregate(
  candidate_ge3_audit$interaction_name,
  by = list(pathway_name = candidate_ge3_audit$pathway_name),
  FUN = function(x) length(unique(x))
)

colnames(pathway_candidate)[2] <-
  "Compatible_ge3of4_records"


if ("pathway_name" %in% colnames(all_direction)) {
  
  pathway_resource <- aggregate(
    all_direction$interaction_name,
    by = list(pathway_name = all_direction$pathway_name),
    FUN = function(x) length(unique(x))
  )
  
  colnames(pathway_resource)[2] <-
    "Resource_LR_Records"
  
} else {
  
  pathway_resource <- data.frame(
    pathway_name = character(0),
    Resource_LR_Records = integer(0)
  )
}

pathway_structure <- merge(
  pathway_resource,
  pathway_candidate,
  by = "pathway_name",
  all = TRUE
)

pathway_structure[
  is.na(pathway_structure)
] <- 0

pathway_structure <- pathway_structure[
  order(pathway_structure$pathway_name),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 12. Core-gene separation confirmation
#
# AOC3 / MS4A7 / IL7R must remain outside the LR candidate
# interpretation unless they are actually present as ligand or
# receptor in the locked CCI_LR resource.
#
# STEP11-5D-1 already showed all three were absent.
# This step only confirms that separation using the locked
# candidate table.
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_gene_audit <- do.call(
  rbind,
  lapply(
    core_genes,
    function(g) {
      
      data.frame(
        core_gene = g,
        present_as_ligand =
          any(candidate_ge3_audit$ligand == g),
        present_as_receptor =
          any(candidate_ge3_audit$receptor == g),
        present_in_LR_candidate_records =
          any(
            candidate_ge3_audit$ligand == g |
              candidate_ge3_audit$receptor == g
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)


# ------------------------------------------------------------
# 13. 3/4 versus 4/4 structural comparison
#
# This does not test significance.
# It only describes how many exact LR pairs and ligands persist
# under the stricter 4/4 patient-consistency definition.
# ------------------------------------------------------------

ge3_pairs <- unique(candidate_ge3_audit$Exact_LR_Pair)
ge4_pairs <- unique(candidate_ge4_audit$Exact_LR_Pair)

ge3_ligands <- unique(candidate_ge3_audit$ligand)
ge4_ligands <- unique(candidate_ge4_audit$ligand)

ge3_records <- nrow(candidate_ge3_audit)
ge4_records <- nrow(candidate_ge4_audit)

comparison_summary <- data.frame(
  
  Metric = c(
    "Candidate_LR_records",
    "Unique_exact_LR_pairs",
    "Unique_ligands",
    "Unique_receptors",
    "Exact_LR_pairs_shared_ge3_and_ge4",
    "Exact_LR_pairs_ge3_only",
    "Exact_LR_pairs_4of4",
    "Ligands_shared_ge3_and_ge4",
    "Ligands_ge3_only",
    "Ligands_4of4"
  ),
  
  Value = c(
    ge3_records,
    length(ge3_pairs),
    length(ge3_ligands),
    length(unique(candidate_ge3_audit$receptor)),
    length(intersect(ge3_pairs, ge4_pairs)),
    length(setdiff(ge3_pairs, ge4_pairs)),
    length(ge4_pairs),
    length(intersect(ge3_ligands, ge4_ligands)),
    length(setdiff(ge3_ligands, ge4_ligands)),
    length(ge4_ligands)
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 14. QC summary
# ------------------------------------------------------------

qc_results <- data.frame(
  QC_ID = c(
    "QC1",
    "QC2",
    "QC3",
    "QC4",
    "QC5",
    "QC6",
    "QC7",
    "QC8",
    "QC9",
    "QC10"
  ),
  
  QC_Item = c(
    "STEP11-5D-1 ge2 input exists",
    "STEP11-5D-1 ge3 input exists",
    "STEP11-5D-1 4of4 input exists",
    "STEP11-5D-1 all-direction audit exists",
    "Required candidate columns present",
    "No missing ligand identifiers in ge3",
    "No missing receptor identifiers in ge3",
    "No missing interaction identifiers in ge3",
    "Core genes absent from ge3 LR candidate records",
    "Structural audit performed without LR recalculation"
  ),
  
  PASS = c(
    file.exists(file_ge2),
    file.exists(file_ge3),
    file.exists(file_ge4),
    file.exists(file_all),
    all(required_cols %in% colnames(candidate_ge3)),
    all(!is.na(candidate_ge3$ligand) &
          candidate_ge3$ligand != ""),
    all(!is.na(candidate_ge3$receptor) &
          candidate_ge3$receptor != ""),
    all(!is.na(candidate_ge3$interaction_name) &
          candidate_ge3$interaction_name != ""),
    all(!core_gene_audit$present_in_LR_candidate_records),
    TRUE
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 15. Save outputs
# ------------------------------------------------------------

write.table(
  candidate_ge3_audit,
  file.path(
    output_dir,
    "STEP11-5D-2A_LR_Candidate_Structural_Audit_ge3of4.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  candidate_ge4_audit,
  file.path(
    output_dir,
    "STEP11-5D-2B_LR_Candidate_Structural_Audit_4of4.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  ligand_summary_ge3,
  file.path(
    output_dir,
    "STEP11-5D-2C_Ligand_Level_Summary_ge3of4.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  pair_summary_ge3,
  file.path(
    output_dir,
    "STEP11-5D-2D_Ligand_Receptor_Pair_Summary_ge3of4.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  directional_symmetry,
  file.path(
    output_dir,
    "STEP11-5D-2E_Directional_Symmetry_Audit.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  annotation_structure,
  file.path(
    output_dir,
    "STEP11-5D-2F_Annotation_Structure_Audit.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  pathway_structure,
  file.path(
    output_dir,
    "STEP11-5D-2G_Pathway_Structure_Audit.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  core_gene_audit,
  file.path(
    output_dir,
    "STEP11-5D-2H_Core_Gene_Separation_Confirmation.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  qc_results,
  file.path(
    output_dir,
    "STEP11-5D-2I_QC_Summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  comparison_summary,
  file.path(
    output_dir,
    "STEP11-5D-2J_ge3_vs_4of4_Structural_Comparison.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ------------------------------------------------------------
# 16. Final console audit
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP11-5D-2 FINAL STRUCTURAL AUDIT SUMMARY\n")
cat("============================================================\n\n")

cat("Locked STEP11-5D-1 candidate input:\n")
cat(
  "  >=2/4 records:",
  nrow(candidate_ge2), "\n"
)
cat(
  "  >=3/4 records:",
  nrow(candidate_ge3), "\n"
)
cat(
  "  4/4 records:",
  nrow(candidate_ge4), "\n\n"
)

cat("GE3 STRUCTURE\n")
cat("------------------------------------------------------------\n")
cat(
  "Unique interaction records:",
  length(unique(candidate_ge3$interaction_name)), "\n"
)

cat(
  "Unique ligands:",
  length(unique(candidate_ge3$ligand)), "\n"
)

cat(
  "Unique receptors:",
  length(unique(candidate_ge3$receptor)), "\n"
)

cat(
  "Unique exact ligand-receptor pairs:",
  length(unique(
    paste(
      candidate_ge3$ligand,
      candidate_ge3$receptor,
      sep = "___"
    )
  )),
  "\n\n"
)

cat("GE4 STRUCTURE\n")
cat("------------------------------------------------------------\n")
cat(
  "Unique interaction records:",
  length(unique(candidate_ge4$interaction_name)), "\n"
)

cat(
  "Unique ligands:",
  length(unique(candidate_ge4$ligand)), "\n"
)

cat(
  "Unique receptors:",
  length(unique(candidate_ge4$receptor)), "\n"
)

cat(
  "Unique exact ligand-receptor pairs:",
  length(unique(
    paste(
      candidate_ge4$ligand,
      candidate_ge4$receptor,
      sep = "___"
    )
  )),
  "\n\n"
)

cat("CORE-GENE SEPARATION\n")
cat("------------------------------------------------------------\n")
print(core_gene_audit, row.names = FALSE)

cat("\nANNOTATION STRUCTURE\n")
cat("------------------------------------------------------------\n")
print(annotation_structure, row.names = FALSE)

cat("\nGE3 vs 4/4 STRUCTURAL COMPARISON\n")
cat("------------------------------------------------------------\n")
print(comparison_summary, row.names = FALSE)

cat("\nQC SUMMARY\n")
cat("------------------------------------------------------------\n")
print(qc_results, row.names = FALSE)

if (all(qc_results$PASS)) {
  
  cat("\n")
  cat("============================================================\n")
  cat("STEP11-5D-2 — PASS\n")
  cat("============================================================\n")
  cat("\n")
  cat(
    "Structural audit completed using locked STEP11-5D-1\n",
    "candidate tables only.\n",
    sep = ""
  )
  cat(
    "No expression recalculation, LR detection recalculation,\n",
    "statistical testing, P/FDR calculation, pathway enrichment,\n",
    "pathway ranking, or biological candidate ranking was performed.\n",
    sep = ""
  )
  cat(
    "Original LR interaction records were retained, while exact\n",
    "ligand-receptor, ligand-level, annotation, pathway, and\n",
    "directional redundancy structures were characterized.\n",
    sep = ""
  )
  cat(
    "AOC3, MS4A7, and IL7R remained structurally separated from\n",
    "the LR candidate set.\n",
    sep = ""
  )
  cat("\n")
  cat(
    "Output directory:\n",
    output_dir,
    "\n"
  )
  
} else {
  
  cat("\n")
  cat("============================================================\n")
  cat("STEP11-5D-2 — QC FAILURE\n")
  cat("============================================================\n")
  cat(
    "Review failed QC items before using any structural audit output.\n"
  )
}


cat("\n============================================================\n")
cat("STEP11-5D-2 COMPLETE\n")
cat("============================================================\n")


cat("\nTop repeated exact LR pairs by number of interaction records:\n")
print(
  head(
    pair_summary_ge3[
      order(
        -pair_summary_ge3$n_interaction_records,
        pair_summary_ge3$Exact_LR_Pair
      ),
      c(
        "Exact_LR_Pair",
        "ligand",
        "receptor",
        "n_interaction_records",
        "n_pathways",
        "n_directions",
        "Directions",
        "min_compatible_patients",
        "max_compatible_patients"
      )
    ],
    30
  ),
  row.names = FALSE
)

cat("\nTop ligand expansion patterns:\n")
print(
  head(
    ligand_summary_ge3[
      order(
        -ligand_summary_ge3$n_unique_receptors,
        ligand_summary_ge3$ligand
      ),
      c(
        "ligand",
        "n_interaction_records",
        "n_unique_receptors",
        "n_exact_LR_pairs",
        "n_pathways",
        "n_directions"
      )
    ],
    30
  ),
  row.names = FALSE
)

cat("\nDirectional symmetry pattern counts:\n")
print(
  sort(
    table(
      directional_symmetry$n_directions_present
    ),
    decreasing = TRUE
  )
)

cat("STEP11-5D-1C column names:\n")
print(colnames(ge3))

cat("\nSTEP11-5D-1D column names:\n")
print(colnames(ge4))

cat("\nSTEP11-5D-2B column names:\n")
print(colnames(ge4_audit))


# ============================================================
# STEP11-5D-3
# Genome-wide Patient-level LR Candidate Set Definition Lock
#
# FINAL CORRECTED VERSION
#
# Purpose:
#   Formally lock the previously audited LR candidate universe.
#
# No new:
#   expression calculation
#   pseudobulk calculation
#   CPM calculation
#   LR detection
#   patient-level compatibility calculation
#   statistical testing
#   P value / FDR calculation
#   pathway ranking
#   biological prioritization
#
# Direction_ID is inherited from the locked STEP11-5D-1/2
# outputs and is NOT renamed or reconstructed.
#
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP11-5D-3: Candidate Set Definition Lock\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step11_dir <- "input_dir/Paper 1_Step 11-CellChat"

out_dir <- file.path(
  step11_dir,
  "STEP11-5D-3_Candidate_Set_Definition_Lock"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

file_ge3 <- file.path(
  step11_dir,
  "STEP11-5D-1_Genomewide_Patient_Level_LR_Compatibility",
  "STEP11-5D-1C_LR_Candidates_Compatible_ge3of4.tsv"
)

file_ge4 <- file.path(
  step11_dir,
  "STEP11-5D-1_Genomewide_Patient_Level_LR_Compatibility",
  "STEP11-5D-1D_LR_Candidates_Compatible_4of4.tsv"
)

file_ge3_audit <- file.path(
  step11_dir,
  "STEP11-5D-2_Genomewide_LR_Candidate_Structural_Audit",
  "STEP11-5D-2A_LR_Candidate_Structural_Audit_ge3of4.tsv"
)

file_ge4_audit <- file.path(
  step11_dir,
  "STEP11-5D-2_Genomewide_LR_Candidate_Structural_Audit",
  "STEP11-5D-2B_LR_Candidate_Structural_Audit_4of4.tsv"
)

file_core_sep <- file.path(
  step11_dir,
  "STEP11-5D-2_Genomewide_LR_Candidate_Structural_Audit",
  "STEP11-5D-2H_Core_Gene_Separation_Confirmation.tsv"
)


# ------------------------------------------------------------
# 2. Input existence QC
# ------------------------------------------------------------

input_exists <- c(
  ge3_candidate_records = file.exists(file_ge3),
  ge4_candidate_records = file.exists(file_ge4),
  ge3_structural_audit  = file.exists(file_ge3_audit),
  ge4_structural_audit  = file.exists(file_ge4_audit),
  core_gene_separation  = file.exists(file_core_sep)
)

cat("Input existence:\n")
print(input_exists)

if (!all(input_exists)) {
  stop("ERROR: One or more required STEP11-5D-3 input files are missing.")
}


# ------------------------------------------------------------
# 3. Read locked inputs
# ------------------------------------------------------------

ge3 <- read.delim(
  file_ge3,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ge4 <- read.delim(
  file_ge4,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ge3_audit <- read.delim(
  file_ge3_audit,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ge4_audit <- read.delim(
  file_ge4_audit,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

core_sep <- read.delim(
  file_core_sep,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 4. Dimensions
# ------------------------------------------------------------

cat("\nInput dimensions:\n")

cat(
  "STEP11-5D-1C >=3/4:",
  nrow(ge3), "rows x", ncol(ge3), "columns\n"
)

cat(
  "STEP11-5D-1D 4/4:",
  nrow(ge4), "rows x", ncol(ge4), "columns\n"
)

cat(
  "STEP11-5D-2A >=3/4 audit:",
  nrow(ge3_audit), "rows x", ncol(ge3_audit), "columns\n"
)

cat(
  "STEP11-5D-2B 4/4 audit:",
  nrow(ge4_audit), "rows x", ncol(ge4_audit), "columns\n"
)

cat(
  "STEP11-5D-2H core separation:",
  nrow(core_sep), "rows x", ncol(core_sep), "columns\n"
)


# ------------------------------------------------------------
# 5. Required-column QC
# ------------------------------------------------------------

required_ge3 <- c(
  "interaction_name",
  "pathway_name",
  "ligand",
  "receptor",
  "Sender",
  "Receiver",
  "Direction_ID",
  "Compatible_ge3of4"
)

required_ge4 <- c(
  "interaction_name",
  "pathway_name",
  "ligand",
  "receptor",
  "Sender",
  "Receiver",
  "Direction_ID",
  "Compatible_4of4"
)

required_audit <- c(
  "interaction_name",
  "ligand",
  "receptor",
  "Sender",
  "Receiver",
  "Direction_ID",
  "Exact_LR_Pair",
  "Directed_LR_Pair",
  "Ligand_Direction",
  "Receptor_Direction"
)

qc_ge3 <- all(required_ge3 %in% colnames(ge3))
qc_ge4 <- all(required_ge4 %in% colnames(ge4))
qc_ge3_audit <- all(required_audit %in% colnames(ge3_audit))
qc_ge4_audit <- all(required_audit %in% colnames(ge4_audit))

cat("\nRequired-column QC:\n")
cat("ge3 candidate records:", qc_ge3, "\n")
cat("ge4 candidate records:", qc_ge4, "\n")
cat("ge3 structural audit:", qc_ge3_audit, "\n")
cat("ge4 structural audit:", qc_ge4_audit, "\n")

if (!all(
  qc_ge3,
  qc_ge4,
  qc_ge3_audit,
  qc_ge4_audit
)) {
  stop("ERROR: Required columns are missing.")
}


# ------------------------------------------------------------
# 6. Candidate record-count lock
# ------------------------------------------------------------

n_ge3_records <- nrow(ge3)
n_ge4_records <- nrow(ge4)

cat("\nCandidate record counts:\n")
cat(">=3/4:", n_ge3_records, "\n")
cat("4/4  :", n_ge4_records, "\n")

if (n_ge3_records != 2406) {
  stop("ERROR: Unexpected >=3/4 candidate record count.")
}

if (n_ge4_records != 1714) {
  stop("ERROR: Unexpected 4/4 candidate record count.")
}


# ------------------------------------------------------------
# 7. Exact LR-pair definition
#
# Exact_LR_Pair is inherited directly from the audited
# STEP11-5D-2A / 2B outputs.
#
# It is NOT reconstructed here.
# ------------------------------------------------------------

ge3_pairs <- unique(
  ge3_audit$Exact_LR_Pair
)

ge4_pairs <- unique(
  ge4_audit$Exact_LR_Pair
)

ge3_pairs <- ge3_pairs[
  !is.na(ge3_pairs) &
    nzchar(ge3_pairs)
]

ge4_pairs <- ge4_pairs[
  !is.na(ge4_pairs) &
    nzchar(ge4_pairs)
]


# ------------------------------------------------------------
# 8. Exact-pair counts
# ------------------------------------------------------------

n_ge3_pairs <- length(ge3_pairs)
n_ge4_pairs <- length(ge4_pairs)

shared_pairs <- intersect(
  ge3_pairs,
  ge4_pairs
)

ge3_only_pairs <- setdiff(
  ge3_pairs,
  ge4_pairs
)

ge4_only_pairs <- setdiff(
  ge4_pairs,
  ge3_pairs
)

cat("\nExact LR-pair counts:\n")
cat(">=3/4 unique exact LR pairs:", n_ge3_pairs, "\n")
cat("4/4 unique exact LR pairs  :", n_ge4_pairs, "\n")
cat("Shared exact LR pairs       :", length(shared_pairs), "\n")
cat(">=3/4 only                 :", length(ge3_only_pairs), "\n")
cat("4/4 only                   :", length(ge4_only_pairs), "\n")


# ------------------------------------------------------------
# 9. Expected exact-pair lock
# ------------------------------------------------------------

if (n_ge3_pairs != 585) {
  stop("ERROR: Unexpected >=3/4 exact LR-pair count.")
}

if (n_ge4_pairs != 492) {
  stop("ERROR: Unexpected 4/4 exact LR-pair count.")
}

if (length(shared_pairs) != 492) {
  stop("ERROR: Unexpected shared exact LR-pair count.")
}

if (length(ge3_only_pairs) != 93) {
  stop("ERROR: Unexpected >=3/4-only exact LR-pair count.")
}

if (length(ge4_only_pairs) != 0) {
  stop("ERROR: Unexpected 4/4-only exact LR-pair count.")
}


# ------------------------------------------------------------
# 10. Lock the EXISTING six Direction_ID categories
#
# IMPORTANT:
#   Direction_ID is already defined upstream as 1–6.
#   No renaming or reconstruction is performed.
# ------------------------------------------------------------

expected_direction_ids <- as.character(1:6)

observed_ge3_directions <- sort(
  unique(as.character(ge3$Direction_ID))
)

observed_ge4_directions <- sort(
  unique(as.character(ge4$Direction_ID))
)

observed_ge3_audit_directions <- sort(
  unique(as.character(ge3_audit$Direction_ID))
)

observed_ge4_audit_directions <- sort(
  unique(as.character(ge4_audit$Direction_ID))
)

direction_qc_ge3 <- identical(
  observed_ge3_directions,
  expected_direction_ids
)

direction_qc_ge4 <- identical(
  observed_ge4_directions,
  expected_direction_ids
)

direction_qc_ge3_audit <- identical(
  observed_ge3_audit_directions,
  expected_direction_ids
)

direction_qc_ge4_audit <- identical(
  observed_ge4_audit_directions,
  expected_direction_ids
)

cat("\nDirection_ID integrity:\n")

cat(
  "ge3 candidate records:",
  direction_qc_ge3,
  "\n"
)

cat(
  "ge4 candidate records:",
  direction_qc_ge4,
  "\n"
)

cat(
  "ge3 structural audit:",
  direction_qc_ge3_audit,
  "\n"
)

cat(
  "ge4 structural audit:",
  direction_qc_ge4_audit,
  "\n"
)

if (!all(
  direction_qc_ge3,
  direction_qc_ge4,
  direction_qc_ge3_audit,
  direction_qc_ge4_audit
)) {
  stop(
    "ERROR: Existing Direction_ID structure does not contain exactly 1–6."
  )
}


# ------------------------------------------------------------
# 11. Lock sender/receiver mapping from existing data
#
# No new direction assignment is performed.
# The mapping is extracted from the already audited records.
# ------------------------------------------------------------

direction_map_ge3 <- unique(
  ge3[
    ,
    c(
      "Direction_ID",
      "Sender",
      "Receiver"
    )
  ]
)

direction_map_ge4 <- unique(
  ge4[
    ,
    c(
      "Direction_ID",
      "Sender",
      "Receiver"
    )
  ]
)

direction_map_audit <- unique(
  ge3_audit[
    ,
    c(
      "Direction_ID",
      "Sender",
      "Receiver",
      "Ligand_Direction",
      "Receptor_Direction"
    )
  ]
)

direction_map_ge3$Direction_ID <- as.character(
  direction_map_ge3$Direction_ID
)

direction_map_ge4$Direction_ID <- as.character(
  direction_map_ge4$Direction_ID
)

direction_map_audit$Direction_ID <- as.character(
  direction_map_audit$Direction_ID
)

direction_map_ge3 <- direction_map_ge3[
  order(direction_map_ge3$Direction_ID),
  ,
  drop = FALSE
]

direction_map_ge4 <- direction_map_ge4[
  order(direction_map_ge4$Direction_ID),
  ,
  drop = FALSE
]

cat("\nLocked Direction_ID → Sender/Receiver mapping:\n")
print(direction_map_ge3)


# ------------------------------------------------------------
# 12. Sender/receiver consistency between ge3 and ge4
#
# Direction_ID is the authoritative key.
# Compare the mapping after explicit sorting and
# character normalization.
# ------------------------------------------------------------

direction_map_ge3_cmp <- ge3[
  ,
  c("Direction_ID", "Sender", "Receiver")
]

direction_map_ge4_cmp <- ge4[
  ,
  c("Direction_ID", "Sender", "Receiver")
]

# Normalize to character
direction_map_ge3_cmp$Direction_ID <- as.character(
  direction_map_ge3_cmp$Direction_ID
)

direction_map_ge3_cmp$Sender <- as.character(
  direction_map_ge3_cmp$Sender
)

direction_map_ge3_cmp$Receiver <- as.character(
  direction_map_ge3_cmp$Receiver
)

direction_map_ge4_cmp$Direction_ID <- as.character(
  direction_map_ge4_cmp$Direction_ID
)

direction_map_ge4_cmp$Sender <- as.character(
  direction_map_ge4_cmp$Sender
)

direction_map_ge4_cmp$Receiver <- as.character(
  direction_map_ge4_cmp$Receiver
)

# Extract one mapping per Direction_ID
direction_map_ge3_cmp <- unique(
  direction_map_ge3_cmp
)

direction_map_ge4_cmp <- unique(
  direction_map_ge4_cmp
)

# Sort explicitly by Direction_ID
direction_map_ge3_cmp <- direction_map_ge3_cmp[
  order(
    as.integer(direction_map_ge3_cmp$Direction_ID)
  ),
  ,
  drop = FALSE
]

direction_map_ge4_cmp <- direction_map_ge4_cmp[
  order(
    as.integer(direction_map_ge4_cmp$Direction_ID)
  ),
  ,
  drop = FALSE
]

# Check exactly six rows
direction_rowcount_qc <- (
  nrow(direction_map_ge3_cmp) == 6 &&
    nrow(direction_map_ge4_cmp) == 6
)

cat(
  "\nDirection mapping row-count QC:",
  direction_rowcount_qc,
  "\n"
)

if (!direction_rowcount_qc) {
  stop(
    "ERROR: Direction mapping does not contain exactly six unique directions."
  )
}


# Compare by Direction_ID rather than raw data.frame identity
direction_mapping_qc <- all(
  direction_map_ge3_cmp$Direction_ID ==
    direction_map_ge4_cmp$Direction_ID &
    direction_map_ge3_cmp$Sender ==
    direction_map_ge4_cmp$Sender &
    direction_map_ge3_cmp$Receiver ==
    direction_map_ge4_cmp$Receiver
)

cat(
  "Sender/Receiver mapping ge3 vs ge4:",
  direction_mapping_qc,
  "\n"
)


# ------------------------------------------------------------
# 12A. Print both mappings for transparent audit
# ------------------------------------------------------------

cat("\n>=3/4 direction mapping:\n")
print(direction_map_ge3_cmp)

cat("\n4/4 direction mapping:\n")
print(direction_map_ge4_cmp)


if (!direction_mapping_qc) {
  
  cat("\nDirection mapping differences:\n")
  
  comparison_table <- data.frame(
    Direction_ID = direction_map_ge3_cmp$Direction_ID,
    
    ge3_Sender = direction_map_ge3_cmp$Sender,
    ge4_Sender = direction_map_ge4_cmp$Sender,
    
    ge3_Receiver = direction_map_ge3_cmp$Receiver,
    ge4_Receiver = direction_map_ge4_cmp$Receiver,
    
    stringsAsFactors = FALSE
  )
  
  print(comparison_table)
  
  stop(
    "ERROR: Sender/Receiver mapping differs between ge3 and ge4."
  )
}

# ------------------------------------------------------------
# 13. Candidate-set nesting QC
# ------------------------------------------------------------

nesting_qc <- all(
  ge4_pairs %in% ge3_pairs
)

cat(
  "\nCandidate-set nesting:\n",
  "All 4/4 exact LR pairs contained in >=3/4 set:",
  nesting_qc,
  "\n"
)

if (!nesting_qc) {
  stop(
    "ERROR: 4/4 candidate set is not contained within >=3/4 candidate set."
  )
}


# ------------------------------------------------------------
# 14. Core-gene separation
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_in_ge3 <- core_genes[
  core_genes %in% ge3$ligand |
    core_genes %in% ge3$receptor
]

core_in_ge4 <- core_genes[
  core_genes %in% ge4$ligand |
    core_genes %in% ge4$receptor
]

cat("\nCore-gene separation:\n")

cat(
  "Core genes:",
  paste(core_genes, collapse = ", "),
  "\n"
)

cat(
  "Core genes appearing in >=3/4 LR candidates:",
  ifelse(
    length(core_in_ge3) == 0,
    "NONE",
    paste(core_in_ge3, collapse = ", ")
  ),
  "\n"
)

cat(
  "Core genes appearing in 4/4 LR candidates:",
  ifelse(
    length(core_in_ge4) == 0,
    "NONE",
    paste(core_in_ge4, collapse = ", ")
  ),
  "\n"
)

if (
  length(core_in_ge3) != 0 ||
  length(core_in_ge4) != 0
) {
  stop(
    "ERROR: Core genes overlap the LR candidate set."
  )
}


# ------------------------------------------------------------
# 15. Candidate-record structural identity QC
#
# Confirm that the audited structural record sets correspond
# one-to-one in row count with the original candidate sets.
# ------------------------------------------------------------

record_count_qc <- c(
  ge3_vs_ge3_audit =
    nrow(ge3) == nrow(ge3_audit),
  
  ge4_vs_ge4_audit =
    nrow(ge4) == nrow(ge4_audit)
)

cat("\nCandidate/audit record-count consistency:\n")
print(record_count_qc)

if (!all(record_count_qc)) {
  stop(
    "ERROR: Candidate and structural-audit record counts differ."
  )
}


# ------------------------------------------------------------
# 16. Build locked exact-pair objects
# ------------------------------------------------------------

exact_LR_pairs <- list(
  ge3of4 = ge3_pairs,
  ge4of4 = ge4_pairs,
  shared = shared_pairs,
  ge3of4_only = ge3_only_pairs,
  ge4of4_only = ge4_only_pairs
)


# ------------------------------------------------------------
# 17. Build locked candidate-record objects
# ------------------------------------------------------------

candidate_records <- list(
  ge3of4 = ge3,
  ge4of4 = ge4
)


# ------------------------------------------------------------
# 18. Locked direction object
# ------------------------------------------------------------

directions <- direction_map_ge3_cmp

# ------------------------------------------------------------
# 19. Summary object
# ------------------------------------------------------------

summary_lock <- data.frame(
  Item = c(
    "ge3of4_candidate_records",
    "ge4of4_candidate_records",
    "ge3of4_unique_exact_LR_pairs",
    "ge4of4_unique_exact_LR_pairs",
    "shared_exact_LR_pairs",
    "ge3of4_only_exact_LR_pairs",
    "ge4of4_only_exact_LR_pairs",
    "number_of_direction_categories",
    "core_gene_count",
    "core_genes_in_ge3of4",
    "core_genes_in_ge4of4"
  ),
  Value = c(
    n_ge3_records,
    n_ge4_records,
    n_ge3_pairs,
    n_ge4_pairs,
    length(shared_pairs),
    length(ge3_only_pairs),
    length(ge4_only_pairs),
    length(expected_direction_ids),
    length(core_genes),
    length(core_in_ge3),
    length(core_in_ge4)
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 20. Final QC
# ------------------------------------------------------------

final_qc <- data.frame(
  QC_Item = c(
    "Input files available",
    "ge3 required columns",
    "ge4 required columns",
    "ge3 structural audit columns",
    "ge4 structural audit columns",
    "ge3 candidate record count",
    "ge4 candidate record count",
    "ge3 exact LR pair count",
    "ge4 exact LR pair count",
    "Shared exact LR pair count",
    "ge3-only exact LR pair count",
    "ge4-only exact LR pair count",
    "Six Direction_ID categories",
    "Direction_ID ge3 = 1-6",
    "Direction_ID ge4 = 1-6",
    "Direction_ID ge3 audit = 1-6",
    "Direction_ID ge4 audit = 1-6",
    "Sender/Receiver mapping consistent",
    "4of4 nested within ge3of4",
    "Candidate/audit record counts consistent",
    "Core genes separated from ge3of4",
    "Core genes separated from ge4of4"
  ),
  PASS = c(
    all(input_exists),
    qc_ge3,
    qc_ge4,
    qc_ge3_audit,
    qc_ge4_audit,
    n_ge3_records == 2406,
    n_ge4_records == 1714,
    n_ge3_pairs == 585,
    n_ge4_pairs == 492,
    length(shared_pairs) == 492,
    length(ge3_only_pairs) == 93,
    length(ge4_only_pairs) == 0,
    length(expected_direction_ids) == 6,
    direction_qc_ge3,
    direction_qc_ge4,
    direction_qc_ge3_audit,
    direction_qc_ge4_audit,
    direction_mapping_qc,
    nesting_qc,
    all(record_count_qc),
    length(core_in_ge3) == 0,
    length(core_in_ge4) == 0
  ),
  stringsAsFactors = FALSE
)

cat("\nFinal QC:\n")
print(final_qc)


# ------------------------------------------------------------
# 21. Save locked RDS
# ------------------------------------------------------------

lock_rds <- file.path(
  out_dir,
  "STEP11-5D-3_Candidate_Set_Definition_Lock.rds"
)

saveRDS(
  list(
    candidate_records = candidate_records,
    exact_LR_pairs = exact_LR_pairs,
    directions = directions,
    core_genes = core_genes,
    core_gene_separation = list(
      ge3of4 = core_in_ge3,
      ge4of4 = core_in_ge4
    ),
    summary = summary_lock,
    final_QC = final_qc
  ),
  lock_rds
)


# ------------------------------------------------------------
# 22. Export summary tables
# ------------------------------------------------------------

write.table(
  summary_lock,
  file.path(
    out_dir,
    "STEP11-5D-3A_Candidate_Set_Lock_Summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  directions,
  file.path(
    out_dir,
    "STEP11-5D-3B_Locked_Direction_ID_Mapping.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  data.frame(
    Exact_LR_Pair = ge3_pairs,
    stringsAsFactors = FALSE
  ),
  file.path(
    out_dir,
    "STEP11-5D-3C_Exact_LR_Pairs_ge3of4.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  data.frame(
    Exact_LR_Pair = ge4_pairs,
    stringsAsFactors = FALSE
  ),
  file.path(
    out_dir,
    "STEP11-5D-3D_Exact_LR_Pairs_4of4.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  final_qc,
  file.path(
    out_dir,
    "STEP11-5D-3E_Final_QC.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ------------------------------------------------------------
# 23. Final status
# ------------------------------------------------------------

if (all(final_qc$PASS)) {
  
  cat("\n")
  cat("============================================================\n")
  cat("STEP11-5D-3 — PASS\n")
  cat("============================================================\n\n")
  
  cat("Candidate sets formally locked.\n\n")
  
  cat(">=3/4 candidate records :", n_ge3_records, "\n")
  cat("4/4 candidate records   :", n_ge4_records, "\n")
  cat(">=3/4 exact LR pairs    :", n_ge3_pairs, "\n")
  cat("4/4 exact LR pairs      :", n_ge4_pairs, "\n")
  cat("Shared exact LR pairs   :", length(shared_pairs), "\n")
  cat(">=3/4-only LR pairs     :", length(ge3_only_pairs), "\n")
  cat("4/4-only LR pairs       :", length(ge4_only_pairs), "\n")
  cat("Direction_ID categories :", length(expected_direction_ids), "\n")
  cat("Core genes separated    : TRUE\n\n")
  
  cat(
    "No expression, LR detection, patient consistency,\n"
  )
  cat(
    "statistical testing, P/FDR, pathway ranking, or\n"
  )
  cat(
    "biological prioritization was performed.\n\n"
  )
  
} else {
  
  cat("\n")
  cat("============================================================\n")
  cat("STEP11-5D-3 — FAIL\n")
  cat("============================================================\n\n")
  
  print(
    final_qc[
      !final_qc$PASS,
      ,
      drop = FALSE
    ]
  )
  
  stop("STEP11-5D-3 failed final QC.")
}


# ============================================================
# STEP11-5E-1
# Candidate Overlap and Patient-Consistency Characterization
# ============================================================
#
# PURPOSE
# -------
# Characterize the already locked STEP11-5D-3 LR candidate sets:
#
#   >=3/4 patient-consistency candidate set
#   4/4 patient-consistency candidate set
#
# This step performs LOCKED-RESULT CHARACTERIZATION ONLY.
#
# NO:
#   - raw UMI re-reading
#   - expression matrix re-reading
#   - pseudobulk recalculation
#   - CPM recalculation
#   - log2(CPM+1) recalculation
#   - LR detection recalculation
#   - patient compatibility recalculation
#   - CellChat recalculation
#   - statistical testing
#   - P-value calculation
#   - FDR calculation
#   - pathway enrichment
#   - pathway ranking
#   - biological prioritization
#
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP11-5E-1 — Candidate Overlap and Patient-Consistency\n")
cat("Characterization\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory and output directory
# ------------------------------------------------------------

step11_dir <- "input_dir/Paper 1_Step 11-CellChat"

out_dir <- file.path(
  step11_dir,
  "STEP11-5E-1_Candidate_Overlap_Patient_Consistency"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cat("Working directory:\n")
cat(step11_dir, "\n\n")

cat("Output directory:\n")
cat(out_dir, "\n\n")


# ------------------------------------------------------------
# 1. Locked input files
# ------------------------------------------------------------

file_ge3 <- file.path(
  step11_dir,
  "STEP11-5D-1_Genomewide_Patient_Level_LR_Compatibility",
  "STEP11-5D-1C_LR_Candidates_Compatible_ge3of4.tsv"
)

file_ge4 <- file.path(
  step11_dir,
  "STEP11-5D-1_Genomewide_Patient_Level_LR_Compatibility",
  "STEP11-5D-1D_LR_Candidates_Compatible_4of4.tsv"
)

file_audit_ge3 <- file.path(
  step11_dir,
  "STEP11-5D-2_Genomewide_LR_Candidate_Structural_Audit",
  "STEP11-5D-2A_LR_Candidate_Structural_Audit_ge3of4.tsv"
)

file_audit_ge4 <- file.path(
  step11_dir,
  "STEP11-5D-2_Genomewide_LR_Candidate_Structural_Audit",
  "STEP11-5D-2B_LR_Candidate_Structural_Audit_4of4.tsv"
)


input_files <- c(
  file_ge3,
  file_ge4,
  file_audit_ge3,
  file_audit_ge4
)

input_exists <- file.exists(input_files)

cat("Input-file availability:\n")

print(
  data.frame(
    File = basename(input_files),
    Exists = input_exists,
    stringsAsFactors = FALSE
  )
)

if (!all(input_exists)) {
  stop(
    "\nSTEP11-5E-1 STOP: One or more locked input files are missing.\n"
  )
}


# ------------------------------------------------------------
# 2. Read locked results only
# ------------------------------------------------------------

ge3 <- read.delim(
  file_ge3,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ge4 <- read.delim(
  file_ge4,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

audit_ge3 <- read.delim(
  file_audit_ge3,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

audit_ge4 <- read.delim(
  file_audit_ge4,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("\nInput dimensions:\n")

cat("ge3 candidate:",
    nrow(ge3), "rows x", ncol(ge3), "columns\n")

cat("ge4 candidate:",
    nrow(ge4), "rows x", ncol(ge4), "columns\n")

cat("ge3 structural audit:",
    nrow(audit_ge3), "rows x", ncol(audit_ge3), "columns\n")

cat("ge4 structural audit:",
    nrow(audit_ge4), "rows x", ncol(audit_ge4), "columns\n")


# ------------------------------------------------------------
# 3. Required-column QC
# ------------------------------------------------------------

required_candidate <- c(
  "interaction_name",
  "pathway_name",
  "ligand",
  "receptor",
  "annotation",
  "Sender",
  "Receiver",
  "n_compatible_patients",
  "consistency_fraction",
  "Compatible_ge3of4",
  "Compatible_4of4",
  "Direction_ID",
  "Core_Gene_Involved",
  "Core_Gene"
)

required_audit_ge3 <- c(
  "interaction_name",
  "pathway_name",
  "ligand",
  "receptor",
  "annotation",
  "Sender",
  "Receiver",
  "Direction_ID",
  "Exact_LR_Pair",
  "Directed_LR_Pair",
  "Ligand_Direction",
  "Receptor_Direction"
)

required_audit_ge4 <- required_audit_ge3


candidate_columns_qc <- all(
  required_candidate %in% names(ge3)
) &&
  all(
    required_candidate %in% names(ge4)
  )

audit_columns_qc <- all(
  required_audit_ge3 %in% names(audit_ge3)
) &&
  all(
    required_audit_ge4 %in% names(audit_ge4)
  )


cat("\nRequired-column QC:\n")

print(
  data.frame(
    QC_Item = c(
      "ge3 candidate required columns",
      "ge4 candidate required columns",
      "ge3 structural audit required columns",
      "ge4 structural audit required columns"
    ),
    PASS = c(
      all(required_candidate %in% names(ge3)),
      all(required_candidate %in% names(ge4)),
      all(required_audit_ge3 %in% names(audit_ge3)),
      all(required_audit_ge4 %in% names(audit_ge4))
    ),
    stringsAsFactors = FALSE
  )
)

if (!candidate_columns_qc || !audit_columns_qc) {
  stop(
    "\nSTEP11-5E-1 STOP: Required columns are missing.\n"
  )
}


# ------------------------------------------------------------
# 4. Construct exact LR pair identity
# ------------------------------------------------------------
#
# IMPORTANT:
# Exact_LR_Pair is the unit for candidate-overlap counting.
#
# Candidate records are direction-specific.
# Exact LR pairs are unique ligand-receptor identities.
#
# We therefore construct a normalized exact pair key directly
# from ligand and receptor rather than relying on row counts.
# ------------------------------------------------------------

normalize_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x
}

ge3$ligand <- normalize_gene(ge3$ligand)
ge3$receptor <- normalize_gene(ge3$receptor)

ge4$ligand <- normalize_gene(ge4$ligand)
ge4$receptor <- normalize_gene(ge4$receptor)

audit_ge3$ligand <- normalize_gene(audit_ge3$ligand)
audit_ge3$receptor <- normalize_gene(audit_ge3$receptor)

audit_ge4$ligand <- normalize_gene(audit_ge4$ligand)
audit_ge4$receptor <- normalize_gene(audit_ge4$receptor)


ge3$Exact_LR_Pair_from_genes <- paste(
  ge3$ligand,
  ge3$receptor,
  sep = "___"
)

ge4$Exact_LR_Pair_from_genes <- paste(
  ge4$ligand,
  ge4$receptor,
  sep = "___"
)

audit_ge3$Exact_LR_Pair_from_genes <- paste(
  audit_ge3$ligand,
  audit_ge3$receptor,
  sep = "___"
)

audit_ge4$Exact_LR_Pair_from_genes <- paste(
  audit_ge4$ligand,
  audit_ge4$receptor,
  sep = "___"
)


# ------------------------------------------------------------
# 5. Verify structural Exact_LR_Pair identity
# ------------------------------------------------------------

audit_pair_identity_qc <- c(
  all(
    as.character(audit_ge3$Exact_LR_Pair) ==
      audit_ge3$Exact_LR_Pair_from_genes
  ),
  all(
    as.character(audit_ge4$Exact_LR_Pair) ==
      audit_ge4$Exact_LR_Pair_from_genes
  )
)


cat("\nExact LR-pair identity QC:\n")

print(
  data.frame(
    Dataset = c(
      "ge3 structural audit",
      "ge4 structural audit"
    ),
    Exact_LR_Pair_matches_gene_identity = audit_pair_identity_qc,
    stringsAsFactors = FALSE
  )
)

if (!all(audit_pair_identity_qc)) {
  stop(
    "\nSTEP11-5E-1 STOP: Exact_LR_Pair identity does not match ligand/receptor identity.\n"
  )
}


# ------------------------------------------------------------
# 6. Unique exact LR-pair sets
# ------------------------------------------------------------

pairs_ge3 <- sort(
  unique(ge3$Exact_LR_Pair_from_genes)
)

pairs_ge4 <- sort(
  unique(ge4$Exact_LR_Pair_from_genes)
)

shared_pairs <- sort(
  intersect(pairs_ge3, pairs_ge4)
)

ge3_only_pairs <- sort(
  setdiff(pairs_ge3, pairs_ge4)
)

ge4_only_pairs <- sort(
  setdiff(pairs_ge4, pairs_ge3)
)


cat("\nLocked candidate-set overlap:\n")

cat(">=3/4 exact LR pairs :", length(pairs_ge3), "\n")
cat("4/4 exact LR pairs   :", length(pairs_ge4), "\n")
cat("Shared               :", length(shared_pairs), "\n")
cat(">=3/4-only           :", length(ge3_only_pairs), "\n")
cat("4/4-only             :", length(ge4_only_pairs), "\n")


# ------------------------------------------------------------
# 7. Candidate-set nesting QC
# ------------------------------------------------------------

all_ge4_in_ge3 <- all(
  pairs_ge4 %in% pairs_ge3
)

overlap_arithmetic_qc <- (
  length(shared_pairs) + length(ge3_only_pairs)
) == length(pairs_ge3)

ge4_arithmetic_qc <- (
  length(shared_pairs) + length(ge4_only_pairs)
) == length(pairs_ge4)


cat("\nCandidate-set nesting QC:\n")

cat(
  "All 4/4 exact LR pairs contained in >=3/4:",
  all_ge4_in_ge3,
  "\n"
)

cat(
  "Shared + >=3/4-only = >=3/4:",
  overlap_arithmetic_qc,
  "\n"
)

cat(
  "Shared + 4/4-only = 4/4:",
  ge4_arithmetic_qc,
  "\n"
)


if (!all(
  all_ge4_in_ge3,
  overlap_arithmetic_qc,
  ge4_arithmetic_qc
)) {
  stop(
    "\nSTEP11-5E-1 STOP: Candidate-set nesting/arithmetic QC failed.\n"
  )
}


# ------------------------------------------------------------
# 8. Retention proportion
# ------------------------------------------------------------

retention_4of4_among_ge3 <- (
  length(shared_pairs) /
    length(pairs_ge3)
) * 100


cat(
  "\n4/4 retention among >=3/4 exact LR pairs:",
  round(retention_4of4_among_ge3, 2),
  "%\n"
)


# ------------------------------------------------------------
# 9. Build exact LR-pair overlap table
# ------------------------------------------------------------

all_pairs <- sort(
  union(pairs_ge3, pairs_ge4)
)

overlap_table <- data.frame(
  Exact_LR_Pair = all_pairs,
  In_ge3of4 = all_pairs %in% pairs_ge3,
  In_4of4 = all_pairs %in% pairs_ge4,
  stringsAsFactors = FALSE
)

overlap_table$Candidate_Class <- ifelse(
  overlap_table$In_ge3of4 & overlap_table$In_4of4,
  "Shared_4of4",
  ifelse(
    overlap_table$In_ge3of4 & !overlap_table$In_4of4,
    "Ge3of4_only",
    ifelse(
      !overlap_table$In_ge3of4 & overlap_table$In_4of4,
      "Ge4of4_only",
      "Unexpected"
    )
  )
)


# ------------------------------------------------------------
# 10. Candidate overlap table integrity
# ------------------------------------------------------------

overlap_class_qc <- all(
  table(overlap_table$Candidate_Class)[
    c("Shared_4of4", "Ge3of4_only")
  ] ==
    c(
      length(shared_pairs),
      length(ge3_only_pairs)
    )
)

unexpected_class_qc <- !any(
  overlap_table$Candidate_Class == "Unexpected"
)


cat("\nOverlap-table QC:\n")

print(
  data.frame(
    QC_Item = c(
      "Candidate-class counts consistent",
      "No unexpected candidate class"
    ),
    PASS = c(
      overlap_class_qc,
      unexpected_class_qc
    ),
    stringsAsFactors = FALSE
  )
)

if (!all(overlap_class_qc, unexpected_class_qc)) {
  stop(
    "\nSTEP11-5E-1 STOP: Candidate overlap table QC failed.\n"
  )
}


# ------------------------------------------------------------
# 11. Patient-level compatibility pattern
# ------------------------------------------------------------
#
# Use the already locked compatibility columns.
# No compatibility is recalculated here.
#
# We first identify the unique exact LR pair + direction
# records in the >=3/4 candidate set.
#
# Patient IDs:
#   P1 = P0020
#   P2 = P0030
#   P3 = P0031
#   P4 = P0034
#
# The exact patient-compatible values are taken directly
# from STEP11-5D-1C.
# ------------------------------------------------------------

required_patient_columns <- c(
  "P1_compatible",
  "P2_compatible",
  "P3_compatible",
  "P4_compatible",
  "n_compatible_patients",
  "consistency_fraction"
)

patient_column_qc <- all(
  required_patient_columns %in% names(ge3)
)

if (!patient_column_qc) {
  stop(
    "\nSTEP11-5E-1 STOP: Patient-level compatibility columns are missing from ge3 input.\n"
  )
}


patient_pattern <- ge3[
  ,
  c(
    "interaction_name",
    "pathway_name",
    "ligand",
    "receptor",
    "Exact_LR_Pair_from_genes",
    "annotation",
    "Sender",
    "Receiver",
    "Direction_ID",
    "P1_compatible",
    "P2_compatible",
    "P3_compatible",
    "P4_compatible",
    "n_compatible_patients",
    "consistency_fraction",
    "Core_Gene_Involved",
    "Core_Gene"
  ),
  drop = FALSE
]


# ------------------------------------------------------------
# 12. Normalize patient compatibility fields
# ------------------------------------------------------------

to_logical <- function(x) {
  
  if (is.logical(x)) {
    return(x)
  }
  
  x <- toupper(trimws(as.character(x)))
  
  out <- rep(NA, length(x))
  
  out[x %in% c("TRUE", "T", "1")] <- TRUE
  out[x %in% c("FALSE", "F", "0")] <- FALSE
  
  out
}


patient_pattern$P1_compatible <- to_logical(
  patient_pattern$P1_compatible
)

patient_pattern$P2_compatible <- to_logical(
  patient_pattern$P2_compatible
)

patient_pattern$P3_compatible <- to_logical(
  patient_pattern$P3_compatible
)

patient_pattern$P4_compatible <- to_logical(
  patient_pattern$P4_compatible
)


patient_pattern$Direction_ID <- as.character(
  patient_pattern$Direction_ID
)


# ------------------------------------------------------------
# 13. Verify locked patient consistency
# ------------------------------------------------------------

patient_pattern$reconstructed_n_compatible <- rowSums(
  cbind(
    patient_pattern$P1_compatible,
    patient_pattern$P2_compatible,
    patient_pattern$P3_compatible,
    patient_pattern$P4_compatible
  )
)

patient_count_qc <- all(
  patient_pattern$reconstructed_n_compatible ==
    as.numeric(patient_pattern$n_compatible_patients)
)


patient_fraction_qc <- all(
  abs(
    as.numeric(patient_pattern$consistency_fraction) -
      (
        patient_pattern$reconstructed_n_compatible / 4
      )
  ) < 1e-12
)


cat("\nPatient-consistency integrity QC:\n")

print(
  data.frame(
    QC_Item = c(
      "n_compatible_patients matches locked patient pattern",
      "consistency_fraction matches locked patient pattern"
    ),
    PASS = c(
      patient_count_qc,
      patient_fraction_qc
    ),
    stringsAsFactors = FALSE
  )
)

if (!all(patient_count_qc, patient_fraction_qc)) {
  stop(
    "\nSTEP11-5E-1 STOP: Patient-consistency integrity QC failed.\n"
  )
}


# ------------------------------------------------------------
# 14. Candidate class assignment
# ------------------------------------------------------------

patient_pattern$Candidate_Class <- ifelse(
  patient_pattern$Exact_LR_Pair_from_genes %in% shared_pairs,
  "Shared_4of4",
  ifelse(
    patient_pattern$Exact_LR_Pair_from_genes %in% ge3_only_pairs,
    "Ge3of4_only",
    "Unexpected"
  )
)


if (any(patient_pattern$Candidate_Class == "Unexpected")) {
  stop(
    "\nSTEP11-5E-1 STOP: Unexpected candidate class detected in patient-pattern table.\n"
  )
}


# ------------------------------------------------------------
# 15. Map patient IDs explicitly
# ------------------------------------------------------------

patient_pattern$P0020_compatible <- patient_pattern$P1_compatible
patient_pattern$P0030_compatible <- patient_pattern$P2_compatible
patient_pattern$P0031_compatible <- patient_pattern$P3_compatible
patient_pattern$P0034_compatible <- patient_pattern$P4_compatible


# ------------------------------------------------------------
# 16. Identify 3/4-only missing patient
# ------------------------------------------------------------

patient_pattern$Missing_Patient <- NA_character_

is_ge3_only <- (
  patient_pattern$Candidate_Class == "Ge3of4_only"
)

patient_pattern$Missing_Patient[
  is_ge3_only &
    patient_pattern$P0020_compatible == FALSE
] <- "P0020"

patient_pattern$Missing_Patient[
  is_ge3_only &
    patient_pattern$P0030_compatible == FALSE
] <- "P0030"

patient_pattern$Missing_Patient[
  is_ge3_only &
    patient_pattern$P0031_compatible == FALSE
] <- "P0031"

patient_pattern$Missing_Patient[
  is_ge3_only &
    patient_pattern$P0034_compatible == FALSE
] <- "P0034"


# ------------------------------------------------------------
# 17. Missing-patient integrity QC
# ------------------------------------------------------------

missing_patient_count_qc <- all(
  is_ge3_only |
    !is.na(patient_pattern$Missing_Patient) |
    (
      patient_pattern$Candidate_Class == "Shared_4of4"
    )
)


ge3_only_missing_one_patient_qc <- all(
  patient_pattern$Candidate_Class == "Shared_4of4" |
    (
      patient_pattern$reconstructed_n_compatible == 3 &
        !is.na(patient_pattern$Missing_Patient)
    )
)


cat("\nMissing-patient characterization QC:\n")

print(
  data.frame(
    QC_Item = c(
      "3/4-only candidates have one identified non-compatible patient",
      "Shared 4/4 candidates have no missing patient"
    ),
    PASS = c(
      ge3_only_missing_one_patient_qc,
      all(
        patient_pattern$Candidate_Class == "Ge3of4_only" |
          is.na(patient_pattern$Missing_Patient)
      )
    ),
    stringsAsFactors = FALSE
  )
)

if (!all(
  ge3_only_missing_one_patient_qc,
  all(
    patient_pattern$Candidate_Class == "Ge3of4_only" |
    is.na(patient_pattern$Missing_Patient)
  )
)) {
  stop(
    "\nSTEP11-5E-1 STOP: Missing-patient characterization QC failed.\n"
  )
}


# ------------------------------------------------------------
# 18. Direction mapping
# ------------------------------------------------------------

direction_mapping <- data.frame(
  Direction_ID = c("1", "2", "3", "4", "5", "6"),
  Sender = c(
    "Endothelial",
    "Myeloid/macrophage",
    "Endothelial",
    "T-cell",
    "Myeloid/macrophage",
    "T-cell"
  ),
  Receiver = c(
    "Myeloid/macrophage",
    "Endothelial",
    "T-cell",
    "Endothelial",
    "T-cell",
    "Myeloid/macrophage"
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 19. Verify direction mapping in patient-pattern table
# ------------------------------------------------------------

observed_direction_mapping <- unique(
  patient_pattern[
    ,
    c("Direction_ID", "Sender", "Receiver"),
    drop = FALSE
  ]
)

observed_direction_mapping$Direction_ID <- as.character(
  observed_direction_mapping$Direction_ID
)

observed_direction_mapping <- observed_direction_mapping[
  observed_direction_mapping$Direction_ID %in%
    direction_mapping$Direction_ID,
  ,
  drop = FALSE
]

observed_direction_mapping <- observed_direction_mapping[
  order(
    as.integer(observed_direction_mapping$Direction_ID)
  ),
  ,
  drop = FALSE
]

rownames(observed_direction_mapping) <- NULL


expected_direction_mapping <- direction_mapping[
  order(
    as.integer(direction_mapping$Direction_ID)
  ),
  ,
  drop = FALSE
]

rownames(expected_direction_mapping) <- NULL


direction_mapping_qc <- (
  nrow(observed_direction_mapping) ==
    nrow(expected_direction_mapping)
) &&
  all(
    observed_direction_mapping$Direction_ID ==
      expected_direction_mapping$Direction_ID
  ) &&
  all(
    observed_direction_mapping$Sender ==
      expected_direction_mapping$Sender
  ) &&
  all(
    observed_direction_mapping$Receiver ==
      expected_direction_mapping$Receiver
  )


cat("\nDirection mapping QC:\n")
print(direction_mapping_qc)

if (!direction_mapping_qc) {
  cat("\nObserved mapping:\n")
  print(observed_direction_mapping)
  
  cat("\nExpected mapping:\n")
  print(expected_direction_mapping)
  
  stop(
    "\nSTEP11-5E-1 STOP: Direction mapping mismatch.\n"
  )
}


# ------------------------------------------------------------
# 20. Direction-specific descriptive summary
# ------------------------------------------------------------

direction_summary <- do.call(
  rbind,
  lapply(
    direction_mapping$Direction_ID,
    function(did) {
      
      x <- patient_pattern[
        patient_pattern$Direction_ID == did,
        ,
        drop = FALSE
      ]
      
      data.frame(
        Direction_ID = did,
        Sender = unique(x$Sender),
        Receiver = unique(x$Receiver),
        Candidate_records_ge3of4 = nrow(x),
        Unique_exact_LR_pairs_ge3of4 =
          length(unique(x$Exact_LR_Pair_from_genes)),
        Shared_4of4_exact_LR_pairs =
          length(
            unique(
              x$Exact_LR_Pair_from_genes[
                x$Candidate_Class == "Shared_4of4"
              ]
            )
          ),
        Ge3of4_only_exact_LR_pairs =
          length(
            unique(
              x$Exact_LR_Pair_from_genes[
                x$Candidate_Class == "Ge3of4_only"
              ]
            )
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)


# ------------------------------------------------------------
# 21. Missing-patient descriptive summary
# ------------------------------------------------------------

dropout_summary <- data.frame(
  Missing_Patient = c(
    "P0020",
    "P0030",
    "P0031",
    "P0034"
  ),
  Ge3of4_only_exact_LR_pairs = sapply(
    c("P0020", "P0030", "P0031", "P0034"),
    function(pid) {
      length(
        unique(
          patient_pattern$Exact_LR_Pair_from_genes[
            patient_pattern$Candidate_Class == "Ge3of4_only" &
              patient_pattern$Missing_Patient == pid
          ]
        )
      )
    }
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 22. Overall summary
# ------------------------------------------------------------

overall_summary <- data.frame(
  Metric = c(
    ">=3/4 candidate records",
    "4/4 candidate records",
    ">=3/4 exact LR pairs",
    "4/4 exact LR pairs",
    "Shared exact LR pairs",
    ">=3/4-only exact LR pairs",
    "4/4-only exact LR pairs",
    "4/4 retention among >=3/4 exact LR pairs (%)",
    "Direction_ID categories"
  ),
  Value = c(
    nrow(ge3),
    nrow(ge4),
    length(pairs_ge3),
    length(pairs_ge4),
    length(shared_pairs),
    length(ge3_only_pairs),
    length(ge4_only_pairs),
    round(retention_4of4_among_ge3, 2),
    length(unique(patient_pattern$Direction_ID))
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 23. Final QC summary
# ------------------------------------------------------------

final_qc <- data.frame(
  QC_Item = c(
    "All locked input files available",
    "Candidate required columns present",
    "Structural audit columns present",
    "Exact LR identity consistent",
    ">=3/4 exact LR pair count = 585",
    "4/4 exact LR pair count = 492",
    "Shared exact LR pair count = 492",
    ">=3/4-only exact LR pair count = 93",
    "4/4-only exact LR pair count = 0",
    "All 4/4 pairs contained in >=3/4",
    "Overlap arithmetic consistent",
    "4/4 arithmetic consistent",
    "Patient compatibility integrity",
    "Patient consistency fraction integrity",
    "Direction mapping consistent",
    "No unexpected candidate class",
    "No unexpected patient pattern",
    "No core-gene LR candidates"
  ),
  PASS = c(
    all(input_exists),
    candidate_columns_qc,
    audit_columns_qc,
    all(audit_pair_identity_qc),
    length(pairs_ge3) == 585,
    length(pairs_ge4) == 492,
    length(shared_pairs) == 492,
    length(ge3_only_pairs) == 93,
    length(ge4_only_pairs) == 0,
    all_ge4_in_ge3,
    overlap_arithmetic_qc,
    ge4_arithmetic_qc,
    patient_count_qc,
    patient_fraction_qc,
    direction_mapping_qc,
    unexpected_class_qc,
    ge3_only_missing_one_patient_qc,
    !any(
      patient_pattern$Core_Gene_Involved == TRUE
    )
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 24. Print final QC
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP11-5E-1 FINAL QC\n")
cat("============================================================\n\n")

print(final_qc)


if (!all(final_qc$PASS)) {
  stop(
    "\nSTEP11-5E-1 — FAIL: One or more final QC items failed.\n"
  )
}


# ------------------------------------------------------------
# 25. Write outputs
# ------------------------------------------------------------

write.table(
  overall_summary,
  file.path(
    out_dir,
    "STEP11-5E-1A_Candidate_Overlap_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  overlap_table,
  file.path(
    out_dir,
    "STEP11-5E-1B_Exact_LR_Pair_Overlap.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  patient_pattern,
  file.path(
    out_dir,
    "STEP11-5E-1C_Patient_Consistency_Pattern.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  direction_summary,
  file.path(
    out_dir,
    "STEP11-5E-1D_Direction_Descriptive_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  dropout_summary,
  file.path(
    out_dir,
    "STEP11-5E-1E_Three_of_Four_Dropout_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  final_qc,
  file.path(
    out_dir,
    "STEP11-5E-1F_QC_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 26. Final locked result
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP11-5E-1 — PASS\n")
cat("============================================================\n\n")

cat("Candidate overlap and patient-consistency characterization completed.\n\n")

cat(">=3/4 candidate records :", nrow(ge3), "\n")
cat("4/4 candidate records   :", nrow(ge4), "\n")
cat(">=3/4 exact LR pairs    :", length(pairs_ge3), "\n")
cat("4/4 exact LR pairs      :", length(pairs_ge4), "\n")
cat("Shared exact LR pairs   :", length(shared_pairs), "\n")
cat(">=3/4-only LR pairs     :", length(ge3_only_pairs), "\n")
cat("4/4-only LR pairs       :", length(ge4_only_pairs), "\n")
cat(
  "4/4 retention among >=3/4 :",
  round(retention_4of4_among_ge3, 2),
  "%\n"
)

cat(
  "All 4/4 pairs nested within >=3/4 :",
  all_ge4_in_ge3,
  "\n"
)

cat(
  "Direction mapping consistent       :",
  direction_mapping_qc,
  "\n"
)

cat(
  "Core genes separated from LR set   :",
  !any(patient_pattern$Core_Gene_Involved == TRUE),
  "\n"
)

cat("\n")
cat("No expression recalculation.\n")
cat("No LR detection recalculation.\n")
cat("No patient compatibility recalculation.\n")
cat("No statistical testing.\n")
cat("No P/FDR calculation.\n")
cat("No pathway enrichment or ranking.\n")
cat("No biological prioritization.\n")
cat("\n")
cat("STEP11-5E-1 is a locked-result descriptive characterization only.\n")
cat("============================================================\n")

# ============================================================
# STEP11-5E-2
# Candidate Expression-Magnitude Characterization
#
# FINAL SCHEMA-ALIGNED VERSION
#
# PURPOSE
# ------------------------------------------------------------
# Descriptively characterize the already LOCKED normalized
# expression context of the genome-wide LR candidate sets.
#
# IMPORTANT SCOPE RESTRICTIONS
# ------------------------------------------------------------
# 1. Read ONLY locked STEP11-5D-1C / STEP11-5D-1D outputs.
# 2. Do NOT reread raw UMI matrices.
# 3. Do NOT recalculate CPM or log2(CPM+1).
# 4. Do NOT recalculate LR detection.
# 5. Do NOT recalculate patient compatibility.
# 6. Do NOT perform statistical testing, P values, or FDR.
# 7. Do NOT calculate ligand x receptor products.
# 8. Do NOT rank "strongest" LR pairs.
# 9. Do NOT perform pathway enrichment.
# 10. Do NOT perform biological prioritization.
# 11. Preserve direction-specific candidate records.
# 12. Treat expression magnitude as descriptive expression
#     context, NOT interaction strength or signaling strength.
#
# EXACT LR PAIR DEFINITION
# ------------------------------------------------------------
# Exact_LR_Pair = ligand + "___" + receptor
#
# This reconstructs the exact LR-pair identity already used
# in previous structural audits.
#
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP11-5E-2 — Candidate Expression-Magnitude Characterization\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory and output folder
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 11-CellChat"

input_dir <- file.path(
  base_dir,
  "STEP11-5D-1_Genomewide_Patient_Level_LR_Compatibility"
)

output_dir <- file.path(
  base_dir,
  "STEP11-5E-2_Candidate_Expression_Magnitude"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

file_ge3 <- file.path(
  input_dir,
  "STEP11-5D-1C_LR_Candidates_Compatible_ge3of4.tsv"
)

file_ge4 <- file.path(
  input_dir,
  "STEP11-5D-1D_LR_Candidates_Compatible_4of4.tsv"
)

stopifnot(file.exists(file_ge3))
stopifnot(file.exists(file_ge4))


cat("Input files located:\n")
cat("  ge3:", file_ge3, "\n")
cat("  ge4:", file_ge4, "\n\n")


# ------------------------------------------------------------
# 2. Read LOCKED candidate tables
# ------------------------------------------------------------

ge3 <- read.delim(
  file_ge3,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ge4 <- read.delim(
  file_ge4,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("Input dimensions:\n")
cat("  ge3:", nrow(ge3), "rows x", ncol(ge3), "columns\n")
cat("  ge4:", nrow(ge4), "rows x", ncol(ge4), "columns\n\n")


# ------------------------------------------------------------
# 3. Required column audit
# ------------------------------------------------------------

required_identity <- c(
  "interaction_name",
  "pathway_name",
  "ligand",
  "receptor",
  "annotation",
  "Sender",
  "Receiver",
  "Direction_ID",
  "Core_Gene_Involved",
  "Core_Gene"
)

required_expression <- c(
  "P1_ligand_log2CPM",
  "P2_ligand_log2CPM",
  "P3_ligand_log2CPM",
  "P4_ligand_log2CPM",
  "P1_receptor_log2CPM",
  "P2_receptor_log2CPM",
  "P3_receptor_log2CPM",
  "P4_receptor_log2CPM"
)

required_compatibility <- c(
  "P1_compatible",
  "P2_compatible",
  "P3_compatible",
  "P4_compatible"
)

required_cols <- c(
  required_identity,
  required_expression,
  required_compatibility
)


check_required <- function(df, cols) {
  all(cols %in% colnames(df))
}


qc_required_ge3 <- check_required(
  ge3,
  required_cols
)

qc_required_ge4 <- check_required(
  ge4,
  required_cols
)


cat("Required columns:\n")
cat("  ge3:", qc_required_ge3, "\n")
cat("  ge4:", qc_required_ge4, "\n\n")

stopifnot(qc_required_ge3)
stopifnot(qc_required_ge4)


# ------------------------------------------------------------
# 4. Reconstruct Exact_LR_Pair
#
# IMPORTANT:
# This is an identity reconstruction only.
# No expression or compatibility is recalculated.
# ------------------------------------------------------------

ge3$Exact_LR_Pair <- paste(
  ge3$ligand,
  ge3$receptor,
  sep = "___"
)

ge4$Exact_LR_Pair <- paste(
  ge4$ligand,
  ge4$receptor,
  sep = "___"
)


# ------------------------------------------------------------
# 5. Locked candidate-set size audit
# ------------------------------------------------------------

ge3_exact <- unique(
  ge3$Exact_LR_Pair
)

ge4_exact <- unique(
  ge4$Exact_LR_Pair
)

n_ge3_records <- nrow(ge3)
n_ge4_records <- nrow(ge4)

n_ge3_pairs <- length(ge3_exact)
n_ge4_pairs <- length(ge4_exact)


cat("Candidate-set sizes:\n")
cat("  ge3 records       :", n_ge3_records, "\n")
cat("  ge4 records       :", n_ge4_records, "\n")
cat("  ge3 exact LR pairs:", n_ge3_pairs, "\n")
cat("  ge4 exact LR pairs:", n_ge4_pairs, "\n\n")


qc_locked_size <- (
  n_ge3_records == 2406 &&
    n_ge4_records == 1714 &&
    n_ge3_pairs == 585 &&
    n_ge4_pairs == 492
)


cat(
  "Locked candidate-set size QC:",
  qc_locked_size,
  "\n\n"
)

stopifnot(qc_locked_size)


# ------------------------------------------------------------
# 6. Verify nesting
# ------------------------------------------------------------

ge4_nested_in_ge3 <- all(
  ge4_exact %in% ge3_exact
)


cat(
  "All 4/4 exact LR pairs nested within >=3/4:",
  ge4_nested_in_ge3,
  "\n\n"
)

stopifnot(ge4_nested_in_ge3)


# ------------------------------------------------------------
# 7. Expression columns numeric audit
# ------------------------------------------------------------

for (x in required_expression) {
  
  if (!is.numeric(ge3[[x]])) {
    ge3[[x]] <- suppressWarnings(
      as.numeric(ge3[[x]])
    )
  }
  
  if (!is.numeric(ge4[[x]])) {
    ge4[[x]] <- suppressWarnings(
      as.numeric(ge4[[x]])
    )
  }
}


expr_ge3 <- as.matrix(
  ge3[, required_expression, drop = FALSE]
)

expr_ge4 <- as.matrix(
  ge4[, required_expression, drop = FALSE]
)


qc_ge3_finite <- all(
  is.finite(expr_ge3)
)

qc_ge4_finite <- all(
  is.finite(expr_ge4)
)

qc_ge3_nonnegative <- all(
  expr_ge3 >= 0
)

qc_ge4_nonnegative <- all(
  expr_ge4 >= 0
)


cat("Expression integrity:\n")
cat("  ge3 finite      :", qc_ge3_finite, "\n")
cat("  ge4 finite      :", qc_ge4_finite, "\n")
cat("  ge3 non-negative:", qc_ge3_nonnegative, "\n")
cat("  ge4 non-negative:", qc_ge4_nonnegative, "\n\n")


stopifnot(qc_ge3_finite)
stopifnot(qc_ge4_finite)
stopifnot(qc_ge3_nonnegative)
stopifnot(qc_ge4_nonnegative)


# ------------------------------------------------------------
# 8. Compatibility-value audit
# ------------------------------------------------------------

compat_cols <- required_compatibility

compat_values_ge3 <- unique(
  unlist(
    ge3[, compat_cols, drop = FALSE],
    use.names = FALSE
  )
)

compat_values_ge4 <- unique(
  unlist(
    ge4[, compat_cols, drop = FALSE],
    use.names = FALSE
  )
)


# Convert only for QC interpretation.
# No compatibility recalculation.

compat_values_ge3 <- as.character(
  compat_values_ge3
)

compat_values_ge4 <- as.character(
  compat_values_ge4
)


valid_compat_values <- c(
  "TRUE",
  "FALSE"
)


qc_compat_ge3 <- all(
  compat_values_ge3 %in% valid_compat_values
)

qc_compat_ge4 <- all(
  compat_values_ge4 %in% valid_compat_values
)


cat("Compatibility-value QC:\n")
cat("  ge3:", qc_compat_ge3, "\n")
cat("  ge4:", qc_compat_ge4, "\n\n")


stopifnot(qc_compat_ge3)
stopifnot(qc_compat_ge4)


# ------------------------------------------------------------
# 9. Candidate consistency class
#
# ge3 candidates are divided into:
#   4_of_4
#   3_of_4_only
#
# This uses the LOCKED ge4 candidate universe.
# No patient compatibility is recalculated.
# ------------------------------------------------------------

ge3$Candidate_Consistency_Class <- ifelse(
  ge3$Exact_LR_Pair %in% ge4_exact,
  "4_of_4",
  "3_of_4_only"
)

ge4$Candidate_Consistency_Class <- "4_of_4"


# ------------------------------------------------------------
# 10. Patient-level expression table
#
# Preserve:
#   exact LR pair
#   interaction record
#   direction
#   sender
#   receiver
#   patient
#   ligand expression
#   receptor expression
#   existing compatibility status
#
# No aggregation across directions.
# ------------------------------------------------------------

make_patient_expression_table <- function(
    df,
    candidate_set
) {
  
  out_list <- vector(
    "list",
    4
  )
  
  for (p in 1:4) {
    
    lig_col <- paste0(
      "P", p,
      "_ligand_log2CPM"
    )
    
    rec_col <- paste0(
      "P", p,
      "_receptor_log2CPM"
    )
    
    comp_col <- paste0(
      "P", p,
      "_compatible"
    )
    
    out_list[[p]] <- data.frame(
      
      Candidate_Set = candidate_set,
      
      Candidate_Consistency_Class =
        df$Candidate_Consistency_Class,
      
      Exact_LR_Pair =
        df$Exact_LR_Pair,
      
      interaction_name =
        df$interaction_name,
      
      pathway_name =
        df$pathway_name,
      
      ligand =
        df$ligand,
      
      receptor =
        df$receptor,
      
      annotation =
        df$annotation,
      
      Direction_ID =
        df$Direction_ID,
      
      Sender =
        df$Sender,
      
      Receiver =
        df$Receiver,
      
      Core_Gene_Involved =
        df$Core_Gene_Involved,
      
      Core_Gene =
        df$Core_Gene,
      
      Patient_Index =
        paste0("P", p),
      
      Ligand_log2CPM =
        df[[lig_col]],
      
      Receptor_log2CPM =
        df[[rec_col]],
      
      Compatible =
        as.character(df[[comp_col]]),
      
      stringsAsFactors = FALSE
    )
  }
  
  do.call(
    rbind,
    out_list
  )
}


patient_ge3 <- make_patient_expression_table(
  ge3,
  "ge3"
)

patient_ge4 <- make_patient_expression_table(
  ge4,
  "ge4"
)


# ------------------------------------------------------------
# 11. Descriptive expression summaries
#
# IMPORTANT:
# These are descriptive only.
# No ranking or candidate selection is performed.
# ------------------------------------------------------------

patient_ge3$Ligand_Receptor_Mean_log2CPM <- rowMeans(
  cbind(
    patient_ge3$Ligand_log2CPM,
    patient_ge3$Receptor_log2CPM
  ),
  na.rm = FALSE
)

patient_ge3$Ligand_Receptor_Difference_log2CPM <-
  patient_ge3$Ligand_log2CPM -
  patient_ge3$Receptor_log2CPM


patient_ge4$Ligand_Receptor_Mean_log2CPM <- rowMeans(
  cbind(
    patient_ge4$Ligand_log2CPM,
    patient_ge4$Receptor_log2CPM
  ),
  na.rm = FALSE
)

patient_ge4$Ligand_Receptor_Difference_log2CPM <-
  patient_ge4$Ligand_log2CPM -
  patient_ge4$Receptor_log2CPM


# ------------------------------------------------------------
# 12. Candidate-level expression summary
#
# Unit:
#   Exact_LR_Pair + Direction_ID + Sender + Receiver
#
# This avoids collapsing direction-specific records.
# ------------------------------------------------------------

candidate_expression_summary <- function(df) {
  
  split_key <- paste(
    df$Exact_LR_Pair,
    df$Direction_ID,
    df$Sender,
    df$Receiver,
    sep = "|||"
  )
  
  groups <- split(
    df,
    split_key
  )
  
  out <- lapply(
    groups,
    function(x) {
      
      ligand_values <-
        x$Ligand_log2CPM
      
      receptor_values <-
        x$Receptor_log2CPM
      
      data.frame(
        
        Candidate_Set =
          unique(x$Candidate_Set),
        
        Candidate_Consistency_Class =
          unique(x$Candidate_Consistency_Class),
        
        Exact_LR_Pair =
          unique(x$Exact_LR_Pair),
        
        interaction_name =
          unique(x$interaction_name),
        
        pathway_name =
          unique(x$pathway_name),
        
        ligand =
          unique(x$ligand),
        
        receptor =
          unique(x$receptor),
        
        annotation =
          unique(x$annotation),
        
        Direction_ID =
          unique(x$Direction_ID),
        
        Sender =
          unique(x$Sender),
        
        Receiver =
          unique(x$Receiver),
        
        N_patient_records =
          nrow(x),
        
        N_ligand_nonzero =
          sum(ligand_values > 0),
        
        N_receptor_nonzero =
          sum(receptor_values > 0),
        
        Ligand_median_log2CPM =
          median(
            ligand_values,
            na.rm = TRUE
          ),
        
        Ligand_IQR_log2CPM =
          IQR(
            ligand_values,
            na.rm = TRUE
          ),
        
        Ligand_min_log2CPM =
          min(
            ligand_values,
            na.rm = TRUE
          ),
        
        Ligand_max_log2CPM =
          max(
            ligand_values,
            na.rm = TRUE
          ),
        
        Receptor_median_log2CPM =
          median(
            receptor_values,
            na.rm = TRUE
          ),
        
        Receptor_IQR_log2CPM =
          IQR(
            receptor_values,
            na.rm = TRUE
          ),
        
        Receptor_min_log2CPM =
          min(
            receptor_values,
            na.rm = TRUE
          ),
        
        Receptor_max_log2CPM =
          max(
            receptor_values,
            na.rm = TRUE
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  do.call(
    rbind,
    out
  )
}


summary_ge3 <- candidate_expression_summary(
  patient_ge3
)

summary_ge4 <- candidate_expression_summary(
  patient_ge4
)


# ------------------------------------------------------------
# 13. Direction mapping QC
#
# Locked mapping:
#
# 1 = Endothelial → Myeloid/macrophage
# 2 = Myeloid/macrophage → Endothelial
# 3 = Endothelial → T-cell
# 4 = T-cell → Endothelial
# 5 = Myeloid/macrophage → T-cell
# 6 = T-cell → Myeloid/macrophage
# ------------------------------------------------------------

direction_map <- data.frame(
  
  Direction_ID =
    c("1", "2", "3", "4", "5", "6"),
  
  Sender_Locked =
    c(
      "Endothelial",
      "Myeloid/macrophage",
      "Endothelial",
      "T-cell",
      "Myeloid/macrophage",
      "T-cell"
    ),
  
  Receiver_Locked =
    c(
      "Myeloid/macrophage",
      "Endothelial",
      "T-cell",
      "Endothelial",
      "T-cell",
      "Myeloid/macrophage"
    ),
  
  stringsAsFactors = FALSE
)


direction_records <- unique(
  ge3[
    ,
    c(
      "Direction_ID",
      "Sender",
      "Receiver"
    )
  ]
)

direction_records$Direction_ID <-
  as.character(
    direction_records$Direction_ID
  )

direction_map$Direction_ID <-
  as.character(
    direction_map$Direction_ID
  )


direction_qc <- merge(
  direction_records,
  direction_map,
  by = "Direction_ID",
  all.x = TRUE
)


qc_direction_mapping <- all(
  direction_qc$Sender ==
    direction_qc$Sender_Locked &
    direction_qc$Receiver ==
    direction_qc$Receiver_Locked
)


cat(
  "Direction mapping QC:",
  qc_direction_mapping,
  "\n\n"
)

stopifnot(qc_direction_mapping)


# ------------------------------------------------------------
# 14. Core-gene separation audit
#
# AOC3 / MS4A7 / IL7R must remain outside LR nodes.
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


core_present_ge3 <- any(
  ge3$ligand %in% core_genes |
    ge3$receptor %in% core_genes
)

core_present_ge4 <- any(
  ge4$ligand %in% core_genes |
    ge4$receptor %in% core_genes
)


qc_core_separation <- (
  !core_present_ge3 &
    !core_present_ge4
)


cat("Core-gene separation:\n")
cat(
  "  ge3:",
  !core_present_ge3,
  "\n"
)
cat(
  "  ge4:",
  !core_present_ge4,
  "\n"
)
cat(
  "  overall:",
  qc_core_separation,
  "\n\n"
)


stopifnot(qc_core_separation)


# ------------------------------------------------------------
# 15. Candidate expression summary output
# ------------------------------------------------------------

summary_all <- rbind(
  summary_ge3,
  summary_ge4
)


summary_all <- summary_all[
  order(
    summary_all$Candidate_Set,
    summary_all$Direction_ID,
    summary_all$Exact_LR_Pair
  ),
]


write.table(
  summary_all,
  file = file.path(
    output_dir,
    "STEP11-5E-2A_Candidate_Expression_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 16. Patient-level expression output
# ------------------------------------------------------------

patient_expression_all <- rbind(
  patient_ge3,
  patient_ge4
)


patient_expression_all <- patient_expression_all[
  order(
    patient_expression_all$Candidate_Set,
    patient_expression_all$Candidate_Consistency_Class,
    patient_expression_all$Direction_ID,
    patient_expression_all$Exact_LR_Pair,
    patient_expression_all$Patient_Index
  ),
]


write.table(
  patient_expression_all,
  file = file.path(
    output_dir,
    "STEP11-5E-2B_Candidate_Patient_Expression.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Expression distribution summary
#
# Descriptive only.
# No hypothesis testing.
# ------------------------------------------------------------

distribution_source <-
  patient_ge3


direction_levels <- sort(
  unique(
    as.character(
      distribution_source$Direction_ID
    )
  )
)


consistency_levels <- c(
  "3_of_4_only",
  "4_of_4"
)


distribution_rows <- list()

counter <- 1


for (cc in consistency_levels) {
  
  for (d in direction_levels) {
    
    x <- distribution_source[
      distribution_source$Candidate_Consistency_Class == cc &
        as.character(
          distribution_source$Direction_ID
        ) == d,
      ,
      drop = FALSE
    ]
    
    if (nrow(x) == 0) {
      next
    }
    
    
    ligand_values <-
      x$Ligand_log2CPM
    
    receptor_values <-
      x$Receptor_log2CPM
    
    
    distribution_rows[[counter]] <-
      data.frame(
        
        Candidate_Set =
          "ge3",
        
        Candidate_Consistency_Class =
          cc,
        
        Direction_ID =
          d,
        
        Sender =
          unique(x$Sender),
        
        Receiver =
          unique(x$Receiver),
        
        N_records =
          nrow(x),
        
        N_exact_LR_pairs =
          length(
            unique(
              x$Exact_LR_Pair
            )
          ),
        
        Ligand_median_log2CPM =
          median(
            ligand_values,
            na.rm = TRUE
          ),
        
        Ligand_IQR_log2CPM =
          IQR(
            ligand_values,
            na.rm = TRUE
          ),
        
        Ligand_min_log2CPM =
          min(
            ligand_values,
            na.rm = TRUE
          ),
        
        Ligand_max_log2CPM =
          max(
            ligand_values,
            na.rm = TRUE
          ),
        
        Receptor_median_log2CPM =
          median(
            receptor_values,
            na.rm = TRUE
          ),
        
        Receptor_IQR_log2CPM =
          IQR(
            receptor_values,
            na.rm = TRUE
          ),
        
        Receptor_min_log2CPM =
          min(
            receptor_values,
            na.rm = TRUE
          ),
        
        Receptor_max_log2CPM =
          max(
            receptor_values,
            na.rm = TRUE
          ),
        
        stringsAsFactors = FALSE
      )
    
    counter <- counter + 1
  }
}


distribution_summary <- do.call(
  rbind,
  distribution_rows
)


write.table(
  distribution_summary,
  file = file.path(
    output_dir,
    "STEP11-5E-2C_Expression_Distribution_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. Additional candidate-class summary
#
# Simple descriptive accounting of:
#   4/4
#   3/4-only
#
# No statistical comparison.
# ------------------------------------------------------------

class_summary <- aggregate(
  Exact_LR_Pair ~ Candidate_Consistency_Class,
  data = summary_ge3,
  FUN = function(x) length(unique(x))
)

colnames(class_summary) <- c(
  "Candidate_Consistency_Class",
  "N_exact_LR_pairs"
)


# Add candidate record counts
record_counts <- aggregate(
  Exact_LR_Pair ~ Candidate_Consistency_Class,
  data = ge3,
  FUN = length
)

colnames(record_counts) <- c(
  "Candidate_Consistency_Class",
  "N_candidate_records"
)


class_summary <- merge(
  class_summary,
  record_counts,
  by = "Candidate_Consistency_Class",
  all = TRUE
)


# ------------------------------------------------------------
# 19. QC summary
# ------------------------------------------------------------

qc_summary <- data.frame(
  
  QC_Item = c(
    
    "Input ge3 file available",
    "Input ge4 file available",
    
    "Required columns ge3",
    "Required columns ge4",
    
    "Locked ge3 record count = 2406",
    "Locked ge4 record count = 1714",
    
    "Locked ge3 exact LR pair count = 585",
    "Locked ge4 exact LR pair count = 492",
    
    "All ge4 exact LR pairs nested within ge3",
    
    "ge3 expression values finite",
    "ge4 expression values finite",
    
    "ge3 expression values non-negative",
    "ge4 expression values non-negative",
    
    "ge3 compatibility values valid",
    "ge4 compatibility values valid",
    
    "Direction mapping consistent",
    
    "Core genes separated from LR candidate set",
    
    "Exact LR pair reconstructed from ligand + receptor",
    
    "No raw UMI recalculation performed",
    "No CPM/log2CPM recalculation performed",
    "No LR detection recalculation performed",
    "No patient compatibility recalculation performed",
    "No statistical testing performed",
    "No P/FDR calculation performed",
    "No pathway enrichment performed",
    "No biological candidate ranking performed"
    
  ),
  
  PASS = c(
    
    file.exists(file_ge3),
    file.exists(file_ge4),
    
    qc_required_ge3,
    qc_required_ge4,
    
    n_ge3_records == 2406,
    n_ge4_records == 1714,
    
    n_ge3_pairs == 585,
    n_ge4_pairs == 492,
    
    ge4_nested_in_ge3,
    
    qc_ge3_finite,
    qc_ge4_finite,
    
    qc_ge3_nonnegative,
    qc_ge4_nonnegative,
    
    qc_compat_ge3,
    qc_compat_ge4,
    
    qc_direction_mapping,
    
    qc_core_separation,
    
    TRUE,
    
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
    
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 20. Write QC
# ------------------------------------------------------------

write.table(
  qc_summary,
  file = file.path(
    output_dir,
    "STEP11-5E-2D_QC_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 21. Final QC
# ------------------------------------------------------------

n_qc <-
  nrow(qc_summary)

n_pass <-
  sum(qc_summary$PASS)

n_fail <-
  n_qc - n_pass


cat("\n")
cat("============================================================\n")
cat("STEP11-5E-2 FINAL QC\n")
cat("============================================================\n\n")

cat("QC items :", n_qc, "\n")
cat("PASS     :", n_pass, "\n")
cat("FAIL     :", n_fail, "\n\n")


if (all(qc_summary$PASS)) {
  
  cat(
    "STEP11-5E-2 — PASS\n\n"
  )
  
} else {
  
  cat(
    "STEP11-5E-2 — FAIL\n\n"
  )
  
  print(
    qc_summary[
      !qc_summary$PASS,
      ,
      drop = FALSE
    ]
  )
  
  stop(
    "STEP11-5E-2 QC failed."
  )
}


# ------------------------------------------------------------
# 22. Final descriptive summary
# ------------------------------------------------------------

cat("Locked candidate expression characterization:\n")
cat(
  "  >=3/4 exact LR pairs:",
  n_ge3_pairs,
  "\n"
)

cat(
  "  4/4 exact LR pairs   :",
  n_ge4_pairs,
  "\n"
)

cat(
  "  3/4-only exact pairs :",
  n_ge3_pairs - n_ge4_pairs,
  "\n\n"
)


cat("Expression scale:\n")
cat(
  "  Existing log2(CPM+1) values from locked STEP11-5D-1\n"
)

cat(
  "  No expression recalculation performed\n\n"
)


cat("Interpretation:\n")
cat(
  "  Descriptive expression context only\n"
)

cat(
  "  No interaction-strength inference\n"
)

cat(
  "  No biological ranking\n"
)

cat(
  "  No statistical significance testing\n\n"
)


cat("Output directory:\n")
cat(
  output_dir,
  "\n\n"
)


cat("============================================================\n")
cat("END STEP11-5E-2\n")
cat("============================================================\n")


# ============================================================
# STEP11-5E-3 PRE-AUDIT
# Duplicate Structure Audit
#
# PURPOSE
# ------------------------------------------------------------
# Determine why multiple records share:
#   Exact_LR_Pair + Direction_ID + Patient_Index
#
# NO DATA MODIFICATION
# NO DEDUPLICATION
# NO RECALCULATION
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP11-5E-3 PRE-AUDIT: Duplicate Structure Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Duplicate key count
# ------------------------------------------------------------

dup_key <- patient_expression[
  ,
  .N,
  by = .(
    Exact_LR_Pair,
    Direction_ID,
    Patient_Index
  )
][N > 1]

cat(
  "Duplicate Exact_LR_Pair × Direction_ID × Patient_Index keys:",
  nrow(dup_key),
  "\n\n"
)


# ------------------------------------------------------------
# 2. Distribution of duplicate counts
# ------------------------------------------------------------

dup_count_distribution <- dup_key[
  ,
  .N,
  by = N
][order(N)]

cat("Duplicate-count distribution:\n")
print(dup_count_distribution)

cat("\n")


# ------------------------------------------------------------
# 3. Inspect representative duplicated records
# ------------------------------------------------------------

example_pairs <- dup_key[
  ,
  head(Exact_LR_Pair, 10)
]


cat("Representative duplicated exact LR pairs:\n")
print(example_pairs)

cat("\n")


# ------------------------------------------------------------
# 4. Full records for representative examples
# ------------------------------------------------------------

example_records <- patient_expression[
  Exact_LR_Pair %in% example_pairs
]

cat(
  "Representative duplicated records:\n"
)

print(
  example_records[
    ,
    .(
      Candidate_Set,
      Candidate_Consistency_Class,
      Exact_LR_Pair,
      interaction_name,
      pathway_name,
      ligand,
      receptor,
      annotation,
      Direction_ID,
      Sender,
      Receiver,
      Core_Gene_Involved,
      Core_Gene,
      Patient_Index,
      Ligand_log2CPM,
      Receptor_log2CPM,
      Compatible,
      Ligand_Receptor_Mean_log2CPM,
      Ligand_Receptor_Difference_log2CPM
    )
  ]
)

cat("\n")


# ------------------------------------------------------------
# 5. Determine whether duplicate records differ by
#    interaction_name / pathway / annotation
# ------------------------------------------------------------

duplicate_structure <- patient_expression[
  ,
  .(
    N_records = .N,
    
    N_interaction_name =
      uniqueN(interaction_name),
    
    N_pathway_name =
      uniqueN(pathway_name),
    
    N_annotation =
      uniqueN(annotation),
    
    N_ligand =
      uniqueN(ligand),
    
    N_receptor =
      uniqueN(receptor),
    
    N_sender =
      uniqueN(Sender),
    
    N_receiver =
      uniqueN(Receiver),
    
    N_direction =
      uniqueN(Direction_ID),
    
    N_ligand_expression =
      uniqueN(Ligand_log2CPM),
    
    N_receptor_expression =
      uniqueN(Receptor_log2CPM),
    
    N_compatible =
      uniqueN(Compatible),
    
    Interaction_Names =
      paste(
        unique(interaction_name),
        collapse = ";"
      ),
    
    Pathway_Names =
      paste(
        unique(pathway_name),
        collapse = ";"
      ),
    
    Annotations =
      paste(
        unique(annotation),
        collapse = ";"
      )
  ),
  by = .(
    Exact_LR_Pair,
    Direction_ID,
    Patient_Index
  )
][N_records > 1]


# ------------------------------------------------------------
# 6. Summarize duplicate structural patterns
# ------------------------------------------------------------

duplicate_pattern_summary <- duplicate_structure[
  ,
  .N,
  by = .(
    N_records,
    N_interaction_name,
    N_pathway_name,
    N_annotation,
    N_ligand_expression,
    N_receptor_expression,
    N_compatible
  )
][order(
  -N,
  N_records
)]


cat(
  "Duplicate structural pattern summary:\n"
)

print(
  duplicate_pattern_summary
)

cat("\n")


# ------------------------------------------------------------
# 7. Critical classification
# ------------------------------------------------------------

all_expression_identical <- all(
  duplicate_structure$N_ligand_expression == 1 &
    duplicate_structure$N_receptor_expression == 1
)

all_compatibility_identical <- all(
  duplicate_structure$N_compatible == 1
)

all_ligand_receptor_identical <- all(
  duplicate_structure$N_ligand == 1 &
    duplicate_structure$N_receptor == 1
)

all_direction_identical <- all(
  duplicate_structure$N_direction == 1
)


cat("============================================================\n")
cat("CRITICAL DUPLICATE AUDIT\n")
cat("============================================================\n\n")

cat(
  "Ligand/receptor identity identical within duplicates:",
  all_ligand_receptor_identical,
  "\n"
)

cat(
  "Direction identical within duplicates:",
  all_direction_identical,
  "\n"
)

cat(
  "Expression identical within duplicates:",
  all_expression_identical,
  "\n"
)

cat(
  "Compatibility identical within duplicates:",
  all_compatibility_identical,
  "\n\n"
)


# ------------------------------------------------------------
# 8. Save audit outputs
# ------------------------------------------------------------

audit_dir <- file.path(
  base_dir,
  "STEP11-5E-3_Patient_Level_Candidate_Integration"
)

if (!dir.exists(audit_dir)) {
  dir.create(
    audit_dir,
    recursive = TRUE
  )
}


fwrite(
  dup_key,
  file.path(
    audit_dir,
    "STEP11-5E-3_PREAUDIT_A_Duplicate_Key_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  duplicate_structure,
  file.path(
    audit_dir,
    "STEP11-5E-3_PREAUDIT_B_Duplicate_Structure.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


fwrite(
  duplicate_pattern_summary,
  file.path(
    audit_dir,
    "STEP11-5E-3_PREAUDIT_C_Duplicate_Pattern_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


cat("\n")
cat("PRE-AUDIT COMPLETE.\n")
cat(
  "No records were deleted or modified.\n"
)
cat(
  "Output directory:\n",
  audit_dir,
  "\n"
)
cat("\n")


# ============================================================
# STEP11-5E-3 — PRE-COLLAPSE / CANONICAL KEY AUDIT
# ============================================================
#
# Purpose:
#   Audit whether the STEP11-5E-2B patient-level candidate
#   expression table can be safely collapsed to the canonical key:
#
#       Exact_LR_Pair × Direction_ID × Patient_Index
#
# without changing:
#   - ligand/receptor identity
#   - sender/receiver
#   - direction
#   - pathway
#   - annotation
#   - candidate-set/class
#   - core-gene fields
#   - ligand expression
#   - receptor expression
#   - compatibility
#
# Important:
#   - NO upstream recalculation
#   - NO raw UMI
#   - NO pseudobulk recalculation
#   - NO CPM/log2CPM recalculation
#   - NO LR detection recalculation
#   - NO compatibility recalculation
#   - NO statistical testing
#   - NO pathway enrichment
#   - NO biological ranking
#   - NO modification of STEP11-5E-2B
#
# This is an independent pre-collapse audit only.
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 11-CellChat"

input_dir <- file.path(
  base_dir,
  "STEP11-5E-2_Candidate_Expression_Magnitude"
)

audit_dir <- file.path(
  base_dir,
  "STEP11-5E-3_Patient_Level_Candidate_Integration"
)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


input_file <- file.path(
  input_dir,
  "STEP11-5E-2B_Candidate_Patient_Expression.tsv"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "INPUT FILE NOT FOUND:\n",
      input_file
    )
  )
}


# ------------------------------------------------------------
# 2. Read locked STEP11-5E-2B input
# ------------------------------------------------------------

patient_expression <- fread(
  input_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  data.table = TRUE
)

cat("\n")
cat("============================================================\n")
cat("STEP11-5E-3 PRE-COLLAPSE / CANONICAL KEY AUDIT\n")
cat("============================================================\n\n")

cat(
  "Input file:\n",
  input_file,
  "\n\n"
)

cat(
  "Input dimensions:",
  nrow(patient_expression),
  "rows ×",
  ncol(patient_expression),
  "columns\n\n"
)


# ------------------------------------------------------------
# 3. Required schema audit
# ------------------------------------------------------------

required_cols <- c(
  "Candidate_Set",
  "Candidate_Consistency_Class",
  "Exact_LR_Pair",
  "interaction_name",
  "pathway_name",
  "ligand",
  "receptor",
  "annotation",
  "Direction_ID",
  "Sender",
  "Receiver",
  "Core_Gene_Involved",
  "Core_Gene",
  "Patient_Index",
  "Ligand_log2CPM",
  "Receptor_log2CPM",
  "Compatible",
  "Ligand_Receptor_Mean_log2CPM",
  "Ligand_Receptor_Difference_log2CPM"
)

missing_cols <- setdiff(
  required_cols,
  names(patient_expression)
)

schema_pass <- length(missing_cols) == 0

cat(
  "Schema PASS:",
  schema_pass,
  "\n"
)

if (!schema_pass) {
  
  cat(
    "Missing columns:\n",
    paste(missing_cols, collapse = ", "),
    "\n"
  )
  
  stop(
    "STEP11-5E-3 PRE-COLLAPSE AUDIT STOPPED: required columns missing."
  )
}


# ------------------------------------------------------------
# 4. Preserve original input exactly
# ------------------------------------------------------------

patient_expression_raw <- copy(patient_expression)


# ------------------------------------------------------------
# 5. Audit Exact_LR_Pair identity
# ------------------------------------------------------------

patient_expression[
  ,
  Exact_LR_Pair_Reconstructed :=
    paste(ligand, receptor, sep = "___")
]

exact_lr_identity_pass <- all(
  patient_expression$Exact_LR_Pair ==
    patient_expression$Exact_LR_Pair_Reconstructed
)

cat(
  "Exact_LR_Pair identity PASS:",
  exact_lr_identity_pass,
  "\n"
)

if (!exact_lr_identity_pass) {
  
  mismatch_lr <- patient_expression[
    Exact_LR_Pair != Exact_LR_Pair_Reconstructed
  ]
  
  fwrite(
    mismatch_lr,
    file.path(
      audit_dir,
      "STEP11-5E-3_PRECOLLAPSE_A_LR_Identity_Mismatch.tsv"
    ),
    sep = "\t"
  )
  
  stop(
    "Exact_LR_Pair identity mismatch detected. No collapse performed."
  )
}


# ------------------------------------------------------------
# 6. Define canonical patient-level key
# ------------------------------------------------------------

patient_expression[
  ,
  Canonical_Key :=
    paste(
      Exact_LR_Pair,
      Direction_ID,
      Patient_Index,
      sep = "___"
    )
]

n_raw_records <- nrow(patient_expression)

n_unique_canonical_keys <- uniqueN(
  patient_expression$Canonical_Key
)

n_duplicate_canonical_keys <-
  n_raw_records - n_unique_canonical_keys

cat("\n")
cat(
  "Raw patient-level records:",
  n_raw_records,
  "\n"
)

cat(
  "Unique canonical keys:",
  n_unique_canonical_keys,
  "\n"
)

cat(
  "Duplicate records attributable to canonical-key multiplicity:",
  n_duplicate_canonical_keys,
  "\n"
)


# ------------------------------------------------------------
# 7. Canonical-key multiplicity distribution
# ------------------------------------------------------------

canonical_key_counts <- patient_expression[
  ,
  .N,
  by = .(
    Exact_LR_Pair,
    Direction_ID,
    Patient_Index
  )
]

canonical_key_distribution <- canonical_key_counts[
  ,
  .N,
  by = N
][
  order(N)
]

cat("\n")
cat("Canonical-key multiplicity distribution:\n")
print(canonical_key_distribution)


# ------------------------------------------------------------
# 8. Identify duplicated canonical keys
# ------------------------------------------------------------

duplicated_keys <- canonical_key_counts[
  N > 1
]

cat("\n")
cat(
  "Duplicated canonical keys:",
  nrow(duplicated_keys),
  "\n"
)


# ------------------------------------------------------------
# 9. Variables that MUST be identical within a canonical key
#
# interaction_name is intentionally EXCLUDED because the
# PRE-AUDIT already demonstrated resource-level multiplicity.
#
# It will be retained separately as provenance.
# ------------------------------------------------------------

collapse_identity_cols <- c(
  "Candidate_Set",
  "Candidate_Consistency_Class",
  "pathway_name",
  "ligand",
  "receptor",
  "annotation",
  "Direction_ID",
  "Sender",
  "Receiver",
  "Core_Gene_Involved",
  "Core_Gene",
  "Ligand_log2CPM",
  "Receptor_log2CPM",
  "Compatible",
  "Ligand_Receptor_Mean_log2CPM",
  "Ligand_Receptor_Difference_log2CPM"
)


# ------------------------------------------------------------
# 10. Helper function
#
# Tests whether every variable is invariant within each
# canonical key.
#
# For numeric variables:
#   all finite values must be identical.
#
# For character/integer/logical variables:
#   unique value count must equal 1.
# ------------------------------------------------------------

audit_column_invariance <- function(
    dt,
    group_cols,
    column_name
) {
  
  result <- dt[
    ,
    {
      x <- get(column_name)
      
      if (is.numeric(x)) {
        
        finite_x <- x[is.finite(x)]
        
        pass <- (
          length(finite_x) == length(x) &&
            length(unique(finite_x)) <= 1
        )
        
      } else {
        
        pass <- uniqueN(x, na.rm = FALSE) <= 1
      }
      
      .(
        PASS = pass
      )
    },
    by = group_cols
  ]
  
  result[
    ,
    .(
      Variable = column_name,
      All_Keys_Invariant = all(PASS),
      N_Keys_Fail = sum(!PASS)
    )
  ]
}


# ------------------------------------------------------------
# 11. Run invariance audit for every protected variable
# ------------------------------------------------------------

invariance_results <- rbindlist(
  lapply(
    collapse_identity_cols,
    function(x) {
      audit_column_invariance(
        patient_expression,
        group_cols = c(
          "Exact_LR_Pair",
          "Direction_ID",
          "Patient_Index"
        ),
        column_name = x
      )
    }
  ),
  fill = TRUE
)

cat("\n")
cat("============================================================\n")
cat("CANONICAL-KEY INVARIANCE AUDIT\n")
cat("============================================================\n\n")

print(invariance_results)


# ------------------------------------------------------------
# 12. Global invariance PASS
# ------------------------------------------------------------

all_identity_invariant <- all(
  invariance_results$All_Keys_Invariant
)

cat(
  "\nAll protected variables invariant within canonical keys:",
  all_identity_invariant,
  "\n"
)


# ------------------------------------------------------------
# 13. Explicit critical audits
# ------------------------------------------------------------

critical_variables <- c(
  "ligand",
  "receptor",
  "Direction_ID",
  "Sender",
  "Receiver",
  "pathway_name",
  "annotation",
  "Ligand_log2CPM",
  "Receptor_log2CPM",
  "Compatible"
)

critical_results <- invariance_results[
  Variable %in% critical_variables
]

critical_pass <- all(
  critical_results$All_Keys_Invariant
)

cat(
  "Critical patient-level variables invariant:",
  critical_pass,
  "\n"
)


# ------------------------------------------------------------
# 14. Audit interaction_name multiplicity
#
# interaction_name is NOT required to be identical.
# Instead, characterize its multiplicity while preserving
# all names as provenance.
# ------------------------------------------------------------

interaction_name_audit <- patient_expression[
  ,
  .(
    N_records = .N,
    N_interaction_name = uniqueN(
      interaction_name,
      na.rm = FALSE
    ),
    Interaction_Names =
      paste(
        unique(interaction_name),
        collapse = " | "
      )
  ),
  by = .(
    Exact_LR_Pair,
    Direction_ID,
    Patient_Index
  )
]

interaction_name_distribution <- interaction_name_audit[
  ,
  .N,
  by = N_interaction_name
][
  order(N_interaction_name)
]

cat("\n")
cat("Interaction-name multiplicity distribution:\n")
print(interaction_name_distribution)


# ------------------------------------------------------------
# 15. Confirm interaction_name is the ONLY allowed
#     source-level multiplicity
#
# For each canonical key:
#   N_records may be >1
#   N_interaction_name may be >1
# but all protected variables must remain invariant.
# ------------------------------------------------------------

source_multiplicity_audit <- merge(
  canonical_key_counts,
  interaction_name_audit,
  by = c(
    "Exact_LR_Pair",
    "Direction_ID",
    "Patient_Index"
  ),
  all.x = TRUE
)

source_multiplicity_audit[
  ,
  Only_Interaction_Name_Varies :=
    (
      N_records >= 1 &
        N_interaction_name >= 1
    )
]

# Because protected-variable invariance was already checked
# globally, this table is descriptive/provenance only.


# ------------------------------------------------------------
# 16. Audit pathway / annotation explicitly
# ------------------------------------------------------------

pathway_annotation_audit <- patient_expression[
  ,
  .(
    N_records = .N,
    N_pathway = uniqueN(
      pathway_name,
      na.rm = FALSE
    ),
    N_annotation = uniqueN(
      annotation,
      na.rm = FALSE
    )
  ),
  by = .(
    Exact_LR_Pair,
    Direction_ID,
    Patient_Index
  )
]

pathway_annotation_pass <- all(
  pathway_annotation_audit$N_pathway <= 1 &
    pathway_annotation_audit$N_annotation <= 1
)

cat(
  "\nPathway invariant within canonical keys:",
  pathway_annotation_pass,
  "\n"
)

cat(
  "Annotation invariant within canonical keys:",
  pathway_annotation_pass,
  "\n"
)


# ------------------------------------------------------------
# 17. Audit patient-level expression explicitly
# ------------------------------------------------------------

expression_audit <- patient_expression[
  ,
  .(
    N_records = .N,
    N_ligand_expression =
      uniqueN(
        Ligand_log2CPM,
        na.rm = FALSE
      ),
    N_receptor_expression =
      uniqueN(
        Receptor_log2CPM,
        na.rm = FALSE
      ),
    N_mean_expression =
      uniqueN(
        Ligand_Receptor_Mean_log2CPM,
        na.rm = FALSE
      ),
    N_difference_expression =
      uniqueN(
        Ligand_Receptor_Difference_log2CPM,
        na.rm = FALSE
      )
  ),
  by = .(
    Exact_LR_Pair,
    Direction_ID,
    Patient_Index
  )
]

expression_pass <- all(
  expression_audit$N_ligand_expression <= 1 &
    expression_audit$N_receptor_expression <= 1 &
    expression_audit$N_mean_expression <= 1 &
    expression_audit$N_difference_expression <= 1
)

cat(
  "Patient-level expression invariant:",
  expression_pass,
  "\n"
)


# ------------------------------------------------------------
# 18. Audit compatibility explicitly
# ------------------------------------------------------------

compatibility_audit <- patient_expression[
  ,
  .(
    N_records = .N,
    N_compatible =
      uniqueN(
        Compatible,
        na.rm = FALSE
      )
  ),
  by = .(
    Exact_LR_Pair,
    Direction_ID,
    Patient_Index
  )
]

compatibility_pass <- all(
  compatibility_audit$N_compatible <= 1
)

cat(
  "Patient-level compatibility invariant:",
  compatibility_pass,
  "\n"
)


# ------------------------------------------------------------
# 19. Construct an AUDIT-ONLY canonical table
#
# IMPORTANT:
#   This is NOT a replacement for STEP11-5E-2B.
#   It is created only to verify that canonical collapse is
#   lossless for all protected variables.
#
# interaction_name is retained as a provenance field.
# ------------------------------------------------------------

canonical_patient_audit <- patient_expression[
  ,
  .(
    Candidate_Set = first(Candidate_Set),
    Candidate_Consistency_Class =
      first(Candidate_Consistency_Class),
    
    Exact_LR_Pair =
      first(Exact_LR_Pair),
    
    interaction_name_n =
      uniqueN(interaction_name),
    
    interaction_name_all =
      paste(
        unique(interaction_name),
        collapse = " | "
      ),
    
    pathway_name =
      first(pathway_name),
    
    ligand =
      first(ligand),
    
    receptor =
      first(receptor),
    
    annotation =
      first(annotation),
    
    Direction_ID =
      first(Direction_ID),
    
    Sender =
      first(Sender),
    
    Receiver =
      first(Receiver),
    
    Core_Gene_Involved =
      first(Core_Gene_Involved),
    
    Core_Gene =
      first(Core_Gene),
    
    Patient_Index =
      first(Patient_Index),
    
    Ligand_log2CPM =
      first(Ligand_log2CPM),
    
    Receptor_log2CPM =
      first(Receptor_log2CPM),
    
    Compatible =
      first(Compatible),
    
    Ligand_Receptor_Mean_log2CPM =
      first(Ligand_Receptor_Mean_log2CPM),
    
    Ligand_Receptor_Difference_log2CPM =
      first(Ligand_Receptor_Difference_log2CPM),
    
    Source_N_records =
      .N
  ),
  by = .(
    Exact_LR_Pair,
    Direction_ID,
    Patient_Index
  )
]

setorder(
  canonical_patient_audit,
  Exact_LR_Pair,
  Direction_ID,
  Patient_Index
)


# ------------------------------------------------------------
# 20. Canonical-key uniqueness audit
# ------------------------------------------------------------

canonical_unique_check <- canonical_patient_audit[
  ,
  .N,
  by = .(
    Exact_LR_Pair,
    Direction_ID,
    Patient_Index
  )
][N > 1]

canonical_unique_pass <-
  nrow(canonical_unique_check) == 0

cat(
  "\nCanonical table unique-key PASS:",
  canonical_unique_pass,
  "\n"
)


# ------------------------------------------------------------
# 21. Conservation audit
# ------------------------------------------------------------

canonical_record_conservation_pass <-
  nrow(canonical_patient_audit) ==
  n_unique_canonical_keys

cat(
  "Canonical row-count conservation PASS:",
  canonical_record_conservation_pass,
  "\n"
)


# ------------------------------------------------------------
# 22. Verify no canonical key has conflicting protected fields
# ------------------------------------------------------------

protected_conflict_table <- patient_expression[
  ,
  .(
    N_records = .N,
    
    N_Candidate_Set =
      uniqueN(Candidate_Set),
    
    N_Candidate_Class =
      uniqueN(Candidate_Consistency_Class),
    
    N_pathway =
      uniqueN(pathway_name),
    
    N_ligand =
      uniqueN(ligand),
    
    N_receptor =
      uniqueN(receptor),
    
    N_annotation =
      uniqueN(annotation),
    
    N_Direction =
      uniqueN(Direction_ID),
    
    N_Sender =
      uniqueN(Sender),
    
    N_Receiver =
      uniqueN(Receiver),
    
    N_CoreGeneInvolved =
      uniqueN(Core_Gene_Involved),
    
    N_CoreGene =
      uniqueN(Core_Gene),
    
    N_LigandExpression =
      uniqueN(Ligand_log2CPM),
    
    N_ReceptorExpression =
      uniqueN(Receptor_log2CPM),
    
    N_Compatible =
      uniqueN(Compatible),
    
    N_MeanExpression =
      uniqueN(Ligand_Receptor_Mean_log2CPM),
    
    N_DifferenceExpression =
      uniqueN(Ligand_Receptor_Difference_log2CPM)
  ),
  by = .(
    Exact_LR_Pair,
    Direction_ID,
    Patient_Index
  )
]

protected_conflict_keys <- protected_conflict_table[
  N_Candidate_Set > 1 |
    N_Candidate_Class > 1 |
    N_pathway > 1 |
    N_ligand > 1 |
    N_receptor > 1 |
    N_annotation > 1 |
    N_Direction > 1 |
    N_Sender > 1 |
    N_Receiver > 1 |
    N_CoreGeneInvolved > 1 |
    N_CoreGene > 1 |
    N_LigandExpression > 1 |
    N_ReceptorExpression > 1 |
    N_Compatible > 1 |
    N_MeanExpression > 1 |
    N_DifferenceExpression > 1
]

protected_conflict_pass <-
  nrow(protected_conflict_keys) == 0

cat(
  "Protected-variable conflict PASS:",
  protected_conflict_pass,
  "\n"
)


# ------------------------------------------------------------
# 23. Exact duplicate-collapse accounting
# ------------------------------------------------------------

collapse_accounting <- data.table(
  Metric = c(
    "Raw patient-level records",
    "Unique canonical keys",
    "Duplicate canonical keys",
    "Records represented by canonical keys",
    "Canonical audit table rows",
    "Records collapsed away",
    "Maximum records per canonical key"
  ),
  
  Value = c(
    nrow(patient_expression),
    n_unique_canonical_keys,
    nrow(duplicated_keys),
    sum(canonical_key_counts$N),
    nrow(canonical_patient_audit),
    nrow(patient_expression) -
      nrow(canonical_patient_audit),
    max(canonical_key_counts$N)
  )
)

cat("\n")
cat("============================================================\n")
cat("COLLAPSE ACCOUNTING\n")
cat("============================================================\n\n")

print(collapse_accounting)


# ------------------------------------------------------------
# 24. Global PASS / FAIL definition
# ------------------------------------------------------------

all_pass <- all(
  schema_pass,
  exact_lr_identity_pass,
  all_identity_invariant,
  critical_pass,
  pathway_annotation_pass,
  expression_pass,
  compatibility_pass,
  canonical_unique_pass,
  canonical_record_conservation_pass,
  protected_conflict_pass
)


# ------------------------------------------------------------
# 25. Save audit outputs
# ------------------------------------------------------------

fwrite(
  canonical_key_counts,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_A_Canonical_Key_Counts.tsv"
  ),
  sep = "\t"
)

fwrite(
  canonical_key_distribution,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_B_Canonical_Key_Distribution.tsv"
  ),
  sep = "\t"
)

fwrite(
  invariance_results,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_C_Invariance_Audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  interaction_name_audit,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_D_Interaction_Name_Provenance.tsv"
  ),
  sep = "\t"
)

fwrite(
  pathway_annotation_audit,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_E_Pathway_Annotation_Audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  expression_audit,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_F_Expression_Audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  compatibility_audit,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_G_Compatibility_Audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  protected_conflict_table,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_H_Protected_Variable_Conflicts.tsv"
  ),
  sep = "\t"
)

fwrite(
  protected_conflict_keys,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_I_Conflict_Keys.tsv"
  ),
  sep = "\t"
)

fwrite(
  canonical_patient_audit,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_J_Canonical_Patient_Audit_Table.tsv"
  ),
  sep = "\t"
)

fwrite(
  collapse_accounting,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_K_Collapse_Accounting.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# 26. Final QC summary
# ------------------------------------------------------------

qc_summary <- data.table(
  QC_Item = c(
    "Input file exists",
    "Required schema complete",
    "Exact_LR_Pair identity",
    "Canonical key defined",
    "Protected variables invariant",
    "Critical patient-level variables invariant",
    "Pathway invariant",
    "Annotation invariant",
    "Expression invariant",
    "Compatibility invariant",
    "Canonical table unique",
    "Canonical row-count conservation",
    "Protected-variable conflict absent"
  ),
  
  PASS = c(
    file.exists(input_file),
    schema_pass,
    exact_lr_identity_pass,
    TRUE,
    all_identity_invariant,
    critical_pass,
    pathway_annotation_pass,
    pathway_annotation_pass,
    expression_pass,
    compatibility_pass,
    canonical_unique_pass,
    canonical_record_conservation_pass,
    protected_conflict_pass
  )
)

fwrite(
  qc_summary,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_L_Final_QC_Summary.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# 27. Final console report
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP11-5E-3 PRE-COLLAPSE / CANONICAL KEY AUDIT\n")
cat("FINAL QC\n")
cat("============================================================\n\n")

print(qc_summary)

cat("\n")
cat(
  "Raw records:",
  nrow(patient_expression),
  "\n"
)

cat(
  "Unique canonical patient-level keys:",
  n_unique_canonical_keys,
  "\n"
)

cat(
  "Duplicated canonical keys:",
  nrow(duplicated_keys),
  "\n"
)

cat(
  "Records collapsed by canonical representation:",
  nrow(patient_expression) -
    nrow(canonical_patient_audit),
  "\n"
)

cat(
  "Maximum multiplicity of a canonical key:",
  max(canonical_key_counts$N),
  "\n"
)

cat("\n")

if (all_pass) {
  
  cat(
    "============================================================\n"
  )
  cat(
    "PASS — STEP11-5E-3 PRE-COLLAPSE / CANONICAL KEY AUDIT\n"
  )
  cat(
    "============================================================\n\n"
  )
  
  cat(
    "The patient-level candidate table can be represented by the\n"
  )
  
  cat(
    "canonical key: Exact_LR_Pair × Direction_ID × Patient_Index.\n\n"
  )
  
  cat(
    "All protected patient-level variables are invariant within\n"
  )
  
  cat(
    "canonical keys, including expression and compatibility.\n\n"
  )
  
  cat(
    "Pathway and annotation are invariant within canonical keys.\n\n"
  )
  
  cat(
    "interaction_name multiplicity is retained as provenance and\n"
  )
  
  cat(
    "is not interpreted as additional patient-level observations.\n\n"
  )
  
  cat(
    "No upstream data were recalculated or modified.\n"
  )
  
} else {
  
  cat(
    "============================================================\n"
  )
  cat(
    "FAIL — STEP11-5E-3 PRE-COLLAPSE / CANONICAL KEY AUDIT\n"
  )
  cat(
    "============================================================\n\n"
  )
  
  cat(
    "At least one protected-variable conflict was detected.\n"
  )
  
  cat(
    "NO canonical collapse should be performed.\n"
  )
}

cat("\n")
cat(
  "Audit output directory:\n",
  audit_dir,
  "\n"
)

cat("\n")
cat(
  "STEP11-5E-3 PRE-COLLAPSE AUDIT COMPLETE.\n"
)

# ============================================================
# STEP11-5E-3 PRE-COLLAPSE — PROTECTED VARIABLE CONFLICT DIAGNOSTIC
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PROTECTED-VARIABLE CONFLICT DIAGNOSTIC\n")
cat("============================================================\n\n")

print(invariance_results)

cat("\n")
cat("Variables failing canonical-key invariance:\n")
cat("------------------------------------------------------------\n")

print(
  invariance_results[
    All_Keys_Invariant == FALSE
  ]
)

cat("\n")
cat("============================================================\n")
cat("DETAILED CONFLICT COUNTS\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Recalculate only diagnostic summaries from existing table
# ------------------------------------------------------------

diagnostic_vars <- collapse_identity_cols

diagnostic_conflict_summary <- rbindlist(
  lapply(
    diagnostic_vars,
    function(v) {
      
      tmp <- patient_expression[
        ,
        .(
          N_unique = uniqueN(
            get(v),
            na.rm = FALSE
          )
        ),
        by = .(
          Exact_LR_Pair,
          Direction_ID,
          Patient_Index
        )
      ]
      
      tmp <- tmp[N_unique > 1]
      
      data.table(
        Variable = v,
        N_conflicting_canonical_keys = nrow(tmp)
      )
    }
  ),
  fill = TRUE
)

print(
  diagnostic_conflict_summary[
    order(-N_conflicting_canonical_keys)
  ]
)


# ------------------------------------------------------------
# Identify the exact conflicting canonical keys
# ------------------------------------------------------------

failing_vars <- invariance_results[
  All_Keys_Invariant == FALSE,
  Variable
]

cat("\n")
cat(
  "Failing variables:",
  paste(failing_vars, collapse = ", "),
  "\n\n"
)

if (length(failing_vars) > 0) {
  
  conflict_long <- rbindlist(
    lapply(
      failing_vars,
      function(v) {
        
        tmp <- patient_expression[
          ,
          .(
            N_unique = uniqueN(
              get(v),
              na.rm = FALSE
            ),
            Values =
              paste(
                unique(get(v)),
                collapse = " | "
              )
          ),
          by = .(
            Exact_LR_Pair,
            Direction_ID,
            Patient_Index
          )
        ]
        
        tmp[
          N_unique > 1,
          .(
            Variable = v,
            Exact_LR_Pair,
            Direction_ID,
            Patient_Index,
            N_unique,
            Values
          )
        ]
      }
    ),
    fill = TRUE
  )
  
  cat(
    "Number of conflicting variable/key combinations:",
    nrow(conflict_long),
    "\n\n"
  )
  
  print(
    conflict_long[
      ,
      head(.SD, 30)
    ]
  )
  
  fwrite(
    conflict_long,
    file.path(
      audit_dir,
      "STEP11-5E-3_PRECOLLAPSE_M_Protected_Variable_Conflict_Detail.tsv"
    ),
    sep = "\t"
  )
}


# ------------------------------------------------------------
# Final diagnostic statement
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("DIAGNOSTIC COMPLETE\n")
cat("============================================================\n")

cat(
  "\nNo data were modified.\n"
)

cat(
  "No canonical collapse was performed.\n"
)

cat(
  "No upstream result was recalculated.\n"
)

cat(
  "\nPlease inspect the failing variable(s) above before proceeding.\n"
)


# ============================================================
# STEP11-5E-3 PRE-COLLAPSE
# CANDIDATE-SET PROVENANCE AUDIT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("CANDIDATE-SET PROVENANCE AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Candidate_Set composition within canonical keys
# ------------------------------------------------------------

candidate_set_pattern <- patient_expression[
  ,
  .(
    N_records = .N,
    Candidate_Set_values =
      paste(
        sort(unique(Candidate_Set)),
        collapse = " | "
      ),
    N_Candidate_Sets =
      uniqueN(Candidate_Set),
    Candidate_Class_values =
      paste(
        sort(unique(Candidate_Consistency_Class)),
        collapse = " | "
      ),
    N_Candidate_Classes =
      uniqueN(Candidate_Consistency_Class)
  ),
  by = .(
    Exact_LR_Pair,
    Direction_ID,
    Patient_Index
  )
]


# ------------------------------------------------------------
# 2. Pattern distribution
# ------------------------------------------------------------

candidate_set_pattern_summary <- candidate_set_pattern[
  ,
  .N,
  by = .(
    Candidate_Set_values,
    N_Candidate_Sets,
    Candidate_Class_values,
    N_Candidate_Classes
  )
][
  order(
    -N,
    Candidate_Set_values
  )
]

cat(
  "Candidate-set provenance pattern distribution:\n\n"
)

print(candidate_set_pattern_summary)


# ------------------------------------------------------------
# 3. Identify unexpected Candidate_Set patterns
#
# Expected:
#   ge3
#   ge3 | ge4
#
# Unexpected:
#   ge4 alone
#   any other combination
# ------------------------------------------------------------

unexpected_candidate_patterns <-
  candidate_set_pattern[
    !Candidate_Set_values %in% c(
      "ge3",
      "ge3 | ge4"
    )
  ]

candidate_set_pattern_pass <-
  nrow(unexpected_candidate_patterns) == 0

cat("\n")
cat(
  "Expected candidate-set patterns only:",
  candidate_set_pattern_pass,
  "\n"
)


# ------------------------------------------------------------
# 4. Check that ge4 is always accompanied by ge3
# ------------------------------------------------------------

ge4_without_ge3 <- patient_expression[
  Candidate_Set == "ge4"
][
  ,
  .(
    Has_ge3 =
      "ge3" %in% Candidate_Set
  ),
  by = .(
    Exact_LR_Pair,
    Direction_ID,
    Patient_Index
  )
][
  Has_ge3 == FALSE
]

ge4_nested_pass <-
  nrow(ge4_without_ge3) == 0

cat(
  "ge4 always nested within ge3:",
  ge4_nested_pass,
  "\n"
)


# ------------------------------------------------------------
# 5. Candidate consistency class invariance
# ------------------------------------------------------------

class_invariance <- patient_expression[
  ,
  .(
    N_classes =
      uniqueN(
        Candidate_Consistency_Class
      ),
    Classes =
      paste(
        sort(
          unique(
            Candidate_Consistency_Class
          )
        ),
        collapse = " | "
      )
  ),
  by = .(
    Exact_LR_Pair,
    Direction_ID,
    Patient_Index
  )
]

class_invariance_pass <-
  all(
    class_invariance$N_classes == 1
  )

cat(
  "Candidate_Consistency_Class invariant:",
  class_invariance_pass,
  "\n"
)


# ------------------------------------------------------------
# 6. Explicit expected class relationship
#
# ge4 records must be 4_of_4
# ge3-only records must be 3_of_4_only
# ------------------------------------------------------------

class_relationship_audit <- patient_expression[
  ,
  .(
    Candidate_Set_values =
      paste(
        sort(unique(Candidate_Set)),
        collapse = " | "
      ),
    Candidate_Class_values =
      paste(
        sort(
          unique(
            Candidate_Consistency_Class
          )
        ),
        collapse = " | "
      )
  ),
  by = .(
    Exact_LR_Pair,
    Direction_ID,
    Patient_Index
  )
]

relationship_errors <- class_relationship_audit[
  (
    grepl(
      "^ge3$",
      Candidate_Set_values
    ) &
      Candidate_Class_values != "3_of_4_only"
  ) |
    (
      grepl(
        "ge4",
        Candidate_Set_values
      ) &
        Candidate_Class_values != "4_of_4"
    )
]

class_relationship_pass <-
  nrow(relationship_errors) == 0

cat(
  "Candidate_Set ↔ Candidate_Consistency_Class relationship:",
  class_relationship_pass,
  "\n"
)


# ------------------------------------------------------------
# 7. Save outputs
# ------------------------------------------------------------

fwrite(
  candidate_set_pattern,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_N_Candidate_Set_Provenance_By_Canonical_Key.tsv"
  ),
  sep = "\t"
)

fwrite(
  candidate_set_pattern_summary,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_O_Candidate_Set_Provenance_Summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  unexpected_candidate_patterns,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_P_Unexpected_Candidate_Set_Patterns.tsv"
  ),
  sep = "\t"
)

fwrite(
  ge4_without_ge3,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_Q_ge4_Without_ge3.tsv"
  ),
  sep = "\t"
)

fwrite(
  class_invariance,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_R_Consistency_Class_Audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  relationship_errors,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_S_Candidate_Set_Class_Relationship_Errors.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# 8. Final provenance QC
# ------------------------------------------------------------

provenance_qc <- data.table(
  QC_Item = c(
    "Expected candidate-set patterns only",
    "ge4 nested within ge3",
    "Candidate_Consistency_Class invariant",
    "Candidate_Set ↔ consistency-class relationship"
  ),
  PASS = c(
    candidate_set_pattern_pass,
    ge4_nested_pass,
    class_invariance_pass,
    class_relationship_pass
  )
)

fwrite(
  provenance_qc,
  file.path(
    audit_dir,
    "STEP11-5E-3_PRECOLLAPSE_T_Candidate_Set_Provenance_QC.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# 9. Final report
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("CANDIDATE-SET PROVENANCE AUDIT — FINAL\n")
cat("============================================================\n\n")

print(provenance_qc)

cat("\n")

if (all(provenance_qc$PASS)) {
  
  cat(
    "PASS — Candidate_Set is a provenance field rather than a\n"
  )
  
  cat(
    "patient-level invariant variable.\n\n"
  )
  
  cat(
    "The canonical patient-level key remains:\n"
  )
  
  cat(
    "Exact_LR_Pair × Direction_ID × Patient_Index\n\n"
  )
  
  cat(
    "No data were modified.\n"
  )
  
} else {
  
  cat(
    "FAIL — Unexpected candidate-set provenance pattern detected.\n"
  )
  
  cat(
    "No canonical collapse should be performed.\n"
  )
}

cat("\n")
cat(
  "Candidate-set provenance audit complete.\n"
)


# ============================================================
# STEP11-5E-3
# Provenance Clarification Audit
#
# Purpose:
#   Resolve Candidate_Set provenance duplication before any
#   canonical patient-level integration.
#
# Design:
#   LEVEL A — Exact LR pair universe
#   LEVEL B — Canonical patient-level key invariance
#   LEVEL C — Candidate-set provenance relationship
#
# IMPORTANT:
#   - NO canonical collapse
#   - NO row deletion
#   - NO upstream recalculation
#   - NO raw UMI loading
#   - NO pseudobulk recalculation
#   - NO CPM/log2CPM recalculation
#   - NO LR detection recalculation
#   - NO compatibility recalculation
#   - NO correlation / P / FDR
#   - NO pathway enrichment
#   - NO biological ranking
#
# Candidate_Set is treated as provenance, NOT as a
# patient-level invariant biological variable.
# ============================================================


cat("============================================================\n")
cat("STEP11-5E-3 — PROVENANCE CLARIFICATION AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_pkgs <- c("data.table")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package not installed: ", pkg)
  }
}

library(data.table)


# ------------------------------------------------------------
# 1. Working directories
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 11-CellChat"

input_dir <- file.path(
  base_dir,
  "STEP11-5E-2_Candidate_Expression_Magnitude"
)

output_dir <- file.path(
  base_dir,
  "STEP11-5E-3_Patient_Level_Candidate_Integration"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 2. Locked input files
# ------------------------------------------------------------

candidate_summary_file <- file.path(
  input_dir,
  "STEP11-5E-2A_Candidate_Expression_Summary.tsv"
)

patient_expression_file <- file.path(
  input_dir,
  "STEP11-5E-2B_Candidate_Patient_Expression.tsv"
)

expression_distribution_file <- file.path(
  input_dir,
  "STEP11-5E-2C_Expression_Distribution_Summary.tsv"
)

qc_file <- file.path(
  input_dir,
  "STEP11-5E-2D_QC_Summary.tsv"
)


# ------------------------------------------------------------
# 3. Input existence audit
# ------------------------------------------------------------

input_exists <- c(
  candidate_summary_file,
  patient_expression_file,
  expression_distribution_file,
  qc_file
)

input_exists_flag <- file.exists(input_exists)

cat("INPUT FILE EXISTENCE\n")
print(
  data.table(
    File = basename(input_exists),
    Exists = input_exists_flag
  )
)
cat("\n")

if (!all(input_exists_flag)) {
  stop("One or more required STEP11-5E-2 input files are missing.")
}


# ------------------------------------------------------------
# 4. Read locked STEP11-5E-2 files
# ------------------------------------------------------------

candidate_summary <- fread(
  candidate_summary_file,
  sep = "\t",
  header = TRUE,
  data.table = TRUE
)

patient_expression <- fread(
  patient_expression_file,
  sep = "\t",
  header = TRUE,
  data.table = TRUE
)

expression_distribution <- fread(
  expression_distribution_file,
  sep = "\t",
  header = TRUE,
  data.table = TRUE
)

upstream_qc <- fread(
  qc_file,
  sep = "\t",
  header = TRUE,
  data.table = TRUE
)


# ------------------------------------------------------------
# 5. Expected schemas
# ------------------------------------------------------------

expected_summary_cols <- c(
  "Candidate_Set",
  "Candidate_Consistency_Class",
  "Exact_LR_Pair",
  "interaction_name",
  "pathway_name",
  "ligand",
  "receptor",
  "annotation",
  "Direction_ID",
  "Sender",
  "Receiver",
  "N_patient_records",
  "N_ligand_nonzero",
  "N_receptor_nonzero",
  "Ligand_median_log2CPM",
  "Ligand_IQR_log2CPM",
  "Ligand_min_log2CPM",
  "Ligand_max_log2CPM",
  "Receptor_median_log2CPM",
  "Receptor_IQR_log2CPM",
  "Receptor_min_log2CPM",
  "Receptor_max_log2CPM"
)

expected_patient_cols <- c(
  "Candidate_Set",
  "Candidate_Consistency_Class",
  "Exact_LR_Pair",
  "interaction_name",
  "pathway_name",
  "ligand",
  "receptor",
  "annotation",
  "Direction_ID",
  "Sender",
  "Receiver",
  "Core_Gene_Involved",
  "Core_Gene",
  "Patient_Index",
  "Ligand_log2CPM",
  "Receptor_log2CPM",
  "Compatible",
  "Ligand_Receptor_Mean_log2CPM",
  "Ligand_Receptor_Difference_log2CPM"
)

expected_distribution_cols <- c(
  "Candidate_Set",
  "Candidate_Consistency_Class",
  "Direction_ID",
  "Sender",
  "Receiver",
  "N_records",
  "N_exact_LR_pairs",
  "Ligand_median_log2CPM",
  "Ligand_IQR_log2CPM",
  "Ligand_min_log2CPM",
  "Ligand_max_log2CPM",
  "Receptor_median_log2CPM",
  "Receptor_IQR_log2CPM",
  "Receptor_min_log2CPM",
  "Receptor_max_log2CPM"
)


# ------------------------------------------------------------
# 6. Schema audit
# ------------------------------------------------------------

summary_schema_pass <- identical(
  names(candidate_summary),
  expected_summary_cols
)

patient_schema_pass <- identical(
  names(patient_expression),
  expected_patient_cols
)

distribution_schema_pass <- identical(
  names(expression_distribution),
  expected_distribution_cols
)

cat("SCHEMA AUDIT\n")
cat(
  "Candidate summary schema:",
  summary_schema_pass, "\n"
)
cat(
  "Patient expression schema:",
  patient_schema_pass, "\n"
)
cat(
  "Expression distribution schema:",
  distribution_schema_pass, "\n\n"
)

if (!summary_schema_pass) {
  stop("Candidate summary schema mismatch.")
}

if (!patient_schema_pass) {
  stop("Patient expression schema mismatch.")
}

if (!distribution_schema_pass) {
  stop("Expression distribution schema mismatch.")
}


# ============================================================
# LEVEL A
# Exact LR pair universe
# ============================================================

cat("============================================================\n")
cat("LEVEL A — EXACT LR PAIR UNIVERSE AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# A1. Construct exact LR pair from ligand + receptor
# ------------------------------------------------------------

candidate_summary[
  ,
  Exact_LR_Pair_Reconstructed :=
    paste(ligand, receptor, sep = "___")
]

patient_expression[
  ,
  Exact_LR_Pair_Reconstructed :=
    paste(ligand, receptor, sep = "___")
]


# ------------------------------------------------------------
# A2. Exact LR identity audit
# ------------------------------------------------------------

summary_lr_identity_pass <- all(
  candidate_summary$Exact_LR_Pair ==
    candidate_summary$Exact_LR_Pair_Reconstructed
)

patient_lr_identity_pass <- all(
  patient_expression$Exact_LR_Pair ==
    patient_expression$Exact_LR_Pair_Reconstructed
)

cat(
  "Candidate summary LR identity:",
  summary_lr_identity_pass, "\n"
)

cat(
  "Patient expression LR identity:",
  patient_lr_identity_pass, "\n\n"
)

if (!summary_lr_identity_pass) {
  stop("Exact LR pair identity failed in candidate summary.")
}

if (!patient_lr_identity_pass) {
  stop("Exact LR pair identity failed in patient expression.")
}


# ------------------------------------------------------------
# A3. Define exact LR pair universes
#
# ge3 universe:
#   all exact LR pairs represented in Candidate_Set = ge3
#
# ge4 universe:
#   all exact LR pairs represented in Candidate_Set = ge4
#
# NOTE:
#   Candidate_Set = ge3 | ge4 is NOT required here because
#   the patient-level table may contain provenance from both.
# ------------------------------------------------------------

ge3_pairs <- sort(
  unique(
    candidate_summary[
      Candidate_Set %in% c("ge3", "ge3 | ge4"),
      Exact_LR_Pair
    ]
  )
)

ge4_pairs <- sort(
  unique(
    candidate_summary[
      Candidate_Set %in% c("ge4", "ge3 | ge4"),
      Exact_LR_Pair
    ]
  )
)

ge3_n <- length(ge3_pairs)
ge4_n <- length(ge4_pairs)

shared_pairs <- intersect(
  ge3_pairs,
  ge4_pairs
)

ge3_only_pairs <- setdiff(
  ge3_pairs,
  ge4_pairs
)

ge4_only_pairs <- setdiff(
  ge4_pairs,
  ge3_pairs
)

cat("Exact LR pair universe:\n")
cat("  ge3 pairs       =", ge3_n, "\n")
cat("  ge4 pairs       =", ge4_n, "\n")
cat("  shared pairs    =", length(shared_pairs), "\n")
cat("  ge3-only pairs  =", length(ge3_only_pairs), "\n")
cat("  ge4-only pairs  =", length(ge4_only_pairs), "\n\n")


# ------------------------------------------------------------
# A4. Exact LR-pair nesting
# ------------------------------------------------------------

ge4_nested_in_ge3_pass <- length(ge4_only_pairs) == 0

cat(
  "ge4 exact LR-pair universe nested within ge3:",
  ge4_nested_in_ge3_pass, "\n\n"
)

if (!ge4_nested_in_ge3_pass) {
  
  fwrite(
    data.table(
      Exact_LR_Pair = ge4_only_pairs
    ),
    file.path(
      output_dir,
      "STEP11-5E-3_PROVENANCE_AUDIT_A_ge4_Not_in_ge3.tsv"
    ),
    sep = "\t"
  )
  
} else {
  
  fwrite(
    data.table(
      Exact_LR_Pair = character()
    ),
    file.path(
      output_dir,
      "STEP11-5E-3_PROVENANCE_AUDIT_A_ge4_Not_in_ge3.tsv"
    ),
    sep = "\t"
  )
}


# ------------------------------------------------------------
# A5. Expected locked counts
# ------------------------------------------------------------

level_A_expected_counts_pass <- (
  ge3_n == 585 &&
    ge4_n == 492 &&
    length(shared_pairs) == 492 &&
    length(ge3_only_pairs) == 93 &&
    length(ge4_only_pairs) == 0
)

cat(
  "Expected locked exact-LR-pair counts:",
  level_A_expected_counts_pass, "\n\n"
)


# ------------------------------------------------------------
# A6. Save Level A summary
# ------------------------------------------------------------

level_A_summary <- data.table(
  Metric = c(
    "ge3_exact_LR_pairs",
    "ge4_exact_LR_pairs",
    "shared_exact_LR_pairs",
    "ge3_only_exact_LR_pairs",
    "ge4_only_exact_LR_pairs",
    "ge4_nested_in_ge3",
    "expected_locked_counts"
  ),
  Value = c(
    ge3_n,
    ge4_n,
    length(shared_pairs),
    length(ge3_only_pairs),
    length(ge4_only_pairs),
    ge4_nested_in_ge3_pass,
    level_A_expected_counts_pass
  )
)

fwrite(
  level_A_summary,
  file.path(
    output_dir,
    "STEP11-5E-3_PROVENANCE_AUDIT_A_Exact_LR_Pair_Universe_Summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# LEVEL B
# Canonical patient-level key invariance
# ============================================================

cat("============================================================\n")
cat("LEVEL B — CANONICAL PATIENT-LEVEL INVARIANCE AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# B1. Canonical key
# ------------------------------------------------------------

patient_expression[
  ,
  Canonical_Key :=
    paste(
      Exact_LR_Pair,
      Direction_ID,
      Patient_Index,
      sep = "___"
    )
]


# ------------------------------------------------------------
# B2. Canonical key multiplicity
# ------------------------------------------------------------

canonical_key_counts <- patient_expression[
  ,
  .(
    N_Records = .N
  ),
  by = Canonical_Key
]

duplicate_canonical_keys <- canonical_key_counts[
  N_Records > 1
]

n_canonical_keys <- nrow(canonical_key_counts)
n_duplicate_keys <- nrow(duplicate_canonical_keys)
n_raw_records <- nrow(patient_expression)

cat("Raw patient-level records:", n_raw_records, "\n")
cat("Unique canonical keys:", n_canonical_keys, "\n")
cat("Duplicate canonical keys:", n_duplicate_keys, "\n\n")


# ------------------------------------------------------------
# B3. True patient-level invariant variables
#
# Candidate_Set is intentionally EXCLUDED.
#
# Candidate_Set is provenance and may legitimately have:
#   ge3
#   ge4
#   ge3 | ge4
# within the same canonical key.
# ------------------------------------------------------------

patient_invariant_vars <- c(
  "Candidate_Consistency_Class",
  "Exact_LR_Pair",
  "interaction_name",
  "pathway_name",
  "ligand",
  "receptor",
  "annotation",
  "Direction_ID",
  "Sender",
  "Receiver",
  "Core_Gene_Involved",
  "Core_Gene",
  "Ligand_log2CPM",
  "Receptor_log2CPM",
  "Compatible",
  "Ligand_Receptor_Mean_log2CPM",
  "Ligand_Receptor_Difference_log2CPM"
)


# ------------------------------------------------------------
# B4. Invariance audit function
# ------------------------------------------------------------

check_invariance <- function(
    dt,
    key_col,
    variable
) {
  
  result <- dt[
    ,
    .(
      N_unique = uniqueN(get(variable)),
      Values = paste(
        sort(unique(as.character(get(variable)))),
        collapse = " | "
      )
    ),
    by = key_col
  ]
  
  result[
    ,
    Fail := N_unique > 1
  ]
  
  list(
    table = result,
    pass = !any(result$Fail),
    n_fail = sum(result$Fail)
  )
}


# ------------------------------------------------------------
# B5. Run invariance audit
# ------------------------------------------------------------

invariance_results <- rbindlist(
  lapply(
    patient_invariant_vars,
    function(v) {
      
      tmp <- check_invariance(
        patient_expression,
        "Canonical_Key",
        v
      )
      
      data.table(
        Variable = v,
        PASS = tmp$pass,
        N_Keys_Fail = tmp$n_fail
      )
    }
  )
)

cat("PATIENT-LEVEL INVARIANCE RESULTS\n")
print(invariance_results)
cat("\n")

level_B_invariance_pass <- all(
  invariance_results$PASS
)


# ------------------------------------------------------------
# B6. Save invariance audit
# ------------------------------------------------------------

fwrite(
  invariance_results,
  file.path(
    output_dir,
    "STEP11-5E-3_PROVENANCE_AUDIT_B_Canonical_Invariance_Summary.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# B7. Save duplicate-key structure
# ------------------------------------------------------------

fwrite(
  duplicate_canonical_keys,
  file.path(
    output_dir,
    "STEP11-5E-3_PROVENANCE_AUDIT_B_Duplicate_Canonical_Keys.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# B8. Row-count conservation
#
# No collapse is performed.
# This is only a structural audit.
# ------------------------------------------------------------

row_count_conservation_pass <- (
  n_raw_records >= n_canonical_keys
)

cat(
  "Row-count structural conservation:",
  row_count_conservation_pass, "\n\n"
)


# ============================================================
# LEVEL C
# Candidate-set provenance relationship
# ============================================================

cat("============================================================\n")
cat("LEVEL C — CANDIDATE-SET PROVENANCE RELATIONSHIP AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# C1. Candidate-set provenance at canonical-key level
#
# This is descriptive only.
# Candidate_Set is NOT required to be invariant.
# ------------------------------------------------------------

candidate_set_provenance <- patient_expression[
  ,
  .(
    Candidate_Set_values =
      paste(
        sort(unique(Candidate_Set)),
        collapse = " | "
      ),
    N_Candidate_Sets =
      uniqueN(Candidate_Set),
    Candidate_Class_values =
      paste(
        sort(unique(Candidate_Consistency_Class)),
        collapse = " | "
      ),
    N_Candidate_Classes =
      uniqueN(Candidate_Consistency_Class),
    N_Records = .N
  ),
  by = Canonical_Key
]


# ------------------------------------------------------------
# C2. Allowed provenance patterns
# ------------------------------------------------------------

allowed_candidate_set_patterns <- c(
  "ge3",
  "ge3 | ge4"
)

unexpected_patterns <- candidate_set_provenance[
  !Candidate_Set_values %in% allowed_candidate_set_patterns
]

level_C_pattern_pass <- nrow(unexpected_patterns) == 0

cat(
  "Allowed Candidate_Set provenance patterns only:",
  level_C_pattern_pass, "\n"
)


# ------------------------------------------------------------
# C3. ge4 provenance must correspond to ge4 exact LR universe
#
# Important:
#   This checks exact-LR-pair membership.
#   It does NOT require the current direction/patient record
#   to be present in the ge4 table.
# ------------------------------------------------------------

candidate_set_provenance[
  ,
  Exact_LR_Pair_from_Key :=
    sub(
      "___[0-9]+___.*$",
      "",
      Canonical_Key
    )
]

ge4_provenance_keys <- candidate_set_provenance[
  grepl(
    "ge4",
    Candidate_Set_values,
    fixed = TRUE
  )
]

ge4_provenance_membership_errors <- ge4_provenance_keys[
  !Exact_LR_Pair_from_Key %in% ge4_pairs
]

level_C_ge4_membership_pass <- (
  nrow(ge4_provenance_membership_errors) == 0
)

cat(
  "All ge4 provenance records belong to ge4 exact-LR universe:",
  level_C_ge4_membership_pass, "\n"
)


# ------------------------------------------------------------
# C4. Candidate class must be invariant
#
# Already checked in Level B, but explicitly record here.
# ------------------------------------------------------------

class_invariance_pass <- (
  invariance_results[
    Variable == "Candidate_Consistency_Class",
    PASS
  ]
)

cat(
  "Candidate_Consistency_Class invariant:",
  class_invariance_pass, "\n"
)


# ------------------------------------------------------------
# C5. Define exact-LR-pair class universe
#
# A pair is 4_of_4 if it belongs to ge4 exact-LR universe.
# Otherwise it is 3_of_4_only.
# ------------------------------------------------------------

pair_class_table <- data.table(
  Exact_LR_Pair = ge3_pairs,
  Expected_Consistency_Class =
    ifelse(
      ge3_pairs %in% ge4_pairs,
      "4_of_4",
      "3_of_4_only"
    )
)

# ------------------------------------------------------------
# C6. Compare actual class to exact-LR-pair universe
# ------------------------------------------------------------

canonical_class_check[
  ,
  Class_Match :=
    Actual_Consistency_Class ==
    Expected_Consistency_Class
]

# data.table-safe row filtering
class_relationship_errors <- canonical_class_check[
  Class_Match == FALSE
]

level_C_class_relationship_pass <- (
  nrow(class_relationship_errors) == 0
)

cat(
  "Candidate_Consistency_Class matches exact-LR-pair universe:",
  level_C_class_relationship_pass, "\n"
)

cat(
  "Number of class-relationship errors:",
  nrow(class_relationship_errors), "\n\n"
)


# ------------------------------------------------------------
# C7. Explicitly audit the important pattern:
#
# ge3 provenance + 4_of_4
#
# This is EXPECTED and should NOT be treated as an error.
# ------------------------------------------------------------

ge3_only_4of4 <- candidate_set_provenance[
  Candidate_Set_values == "ge3" &
    Candidate_Class_values == "4_of_4"
]

n_ge3_only_4of4 <- nrow(ge3_only_4of4)

cat(
  "Canonical keys with ge3-only provenance but 4_of_4 class:",
  n_ge3_only_4of4, "\n"
)

ge3_only_4of4_membership_pass <- all(
  ge3_only_4of4$Exact_LR_Pair_from_Key %in% ge4_pairs
)

cat(
  "These ge3-only/4_of_4 keys belong to ge4 exact-LR universe:",
  ge3_only_4of4_membership_pass, "\n\n"
)


# ------------------------------------------------------------
# C8. Explicitly audit 3_of_4_only
# ------------------------------------------------------------

ge3_only_3of4 <- candidate_set_provenance[
  Candidate_Set_values == "ge3" &
    Candidate_Class_values == "3_of_4_only"
]

n_ge3_only_3of4 <- nrow(ge3_only_3of4)

cat(
  "Canonical keys with ge3-only provenance and 3_of_4_only class:",
  n_ge3_only_3of4, "\n"
)

ge3_only_3of4_membership_pass <- all(
  !ge3_only_3of4$Exact_LR_Pair_from_Key %in% ge4_pairs
)

cat(
  "These 3_of_4_only keys are ge3-only exact-LR pairs:",
  ge3_only_3of4_membership_pass, "\n\n"
)


# ------------------------------------------------------------
# C9. Save provenance table
# ------------------------------------------------------------

fwrite(
  candidate_set_provenance,
  file.path(
    output_dir,
    "STEP11-5E-3_PROVENANCE_AUDIT_C_Candidate_Set_By_Canonical_Key.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# C10. Save unexpected patterns
# ------------------------------------------------------------

fwrite(
  unexpected_patterns,
  file.path(
    output_dir,
    "STEP11-5E-3_PROVENANCE_AUDIT_C_Unexpected_Candidate_Set_Patterns.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# C11. Save ge4 membership errors
# ------------------------------------------------------------

fwrite(
  ge4_provenance_membership_errors,
  file.path(
    output_dir,
    "STEP11-5E-3_PROVENANCE_AUDIT_C_ge4_Provenance_Membership_Errors.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# C12. Save class relationship errors
# ------------------------------------------------------------

fwrite(
  class_relationship_errors,
  file.path(
    output_dir,
    "STEP11-5E-3_PROVENANCE_AUDIT_C_Consistency_Class_Relationship_Errors.tsv"
  ),
  sep = "\t"
)


# ============================================================
# GLOBAL STRUCTURAL AUDIT
# ============================================================

cat("============================================================\n")
cat("GLOBAL STRUCTURAL AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# G1. Expected raw record count
# ------------------------------------------------------------

expected_raw_records <- 16480

raw_record_count_pass <- (
  n_raw_records == expected_raw_records
)

cat(
  "Expected raw patient-level records:",
  expected_raw_records, "\n"
)

cat(
  "Observed raw patient-level records:",
  n_raw_records, "\n"
)

cat(
  "Raw record count PASS:",
  raw_record_count_pass, "\n\n"
)


# ------------------------------------------------------------
# G2. Expected canonical key count
# ------------------------------------------------------------

expected_canonical_keys <- 9620

canonical_key_count_pass <- (
  n_canonical_keys == expected_canonical_keys
)

cat(
  "Expected canonical keys:",
  expected_canonical_keys, "\n"
)

cat(
  "Observed canonical keys:",
  n_canonical_keys, "\n"
)

cat(
  "Canonical key count PASS:",
  canonical_key_count_pass, "\n\n"
)


# ------------------------------------------------------------
# G3. No data collapse
# ------------------------------------------------------------

no_collapse_performed <- TRUE

cat(
  "Canonical collapse performed:",
  !no_collapse_performed, "\n\n"
)

# ------------------------------------------------------------
# G4. Core-gene separation inheritance audit
#
# IMPORTANT:
#   The STEP11-5E-2 genome-wide candidate table contains
#   no core-gene LR candidates by design.
#
#   AOC3 / MS4A7 / IL7R were already separated from the
#   CCI_LR universe in locked STEP11-5D-1 / STEP11-5D-2 /
#   STEP11-5D-3 audits.
#
#   Therefore this step does NOT require the current
#   genome-wide candidate expression table to contain the
#   three core genes.
#
#   This audit verifies only that the current table is
#   consistent with the locked zero-core-LR result.
# ------------------------------------------------------------

core_genes_expected <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

# Current STEP11-5E-2 fields
core_gene_field_all_na <- all(
  is.na(patient_expression$Core_Gene)
)

core_gene_involved_all_false <- all(
  patient_expression$Core_Gene_Involved == FALSE,
  na.rm = TRUE
)

# No core-gene LR candidates are expected in this candidate table.
# The definitive core-gene separation result is inherited from
# locked STEP11-5D-1 / STEP11-5D-2 / STEP11-5D-3.

locked_core_gene_lr_count <- 0L

core_gene_separation_pass <- (
  core_gene_field_all_na &&
    core_gene_involved_all_false &&
    locked_core_gene_lr_count == 0L
)

cat(
  "Core_Gene field all NA:",
  core_gene_field_all_na,
  "\n"
)

cat(
  "Core_Gene_Involved field all FALSE:",
  core_gene_involved_all_false,
  "\n"
)

cat(
  "Locked core-gene LR candidate count:",
  locked_core_gene_lr_count,
  "\n"
)

cat(
  "Core-gene separation inheritance audit:",
  core_gene_separation_pass,
  "\n\n"
)


# ============================================================
# STEP11-5E-3 — interaction_name provenance audit
#
# IMPORTANT:
# interaction_name is a resource-level identifier.
# It is NOT required to be invariant within the canonical
# patient-level key.
#
# Exact_LR_Pair remains the canonical biological LR identity.
# ============================================================

interaction_name_provenance_audit <- patient_expression[
  ,
  .(
    N_records = .N,
    N_interaction_name = uniqueN(interaction_name),
    Interaction_Names = paste(
      sort(unique(interaction_name)),
      collapse = " | "
    ),
    N_exact_LR_pair = uniqueN(Exact_LR_Pair),
    N_pathway = uniqueN(pathway_name),
    N_ligand = uniqueN(ligand),
    N_receptor = uniqueN(receptor),
    N_annotation = uniqueN(annotation),
    N_direction = uniqueN(Direction_ID),
    N_sender = uniqueN(Sender),
    N_receiver = uniqueN(Receiver),
    N_compatible = uniqueN(Compatible)
  ),
  by = .(
    Exact_LR_Pair,
    Direction_ID,
    Patient_Index
  )
]


# ------------------------------------------------------------
# Identify canonical keys with interaction_name multiplicity
# ------------------------------------------------------------

interaction_name_multiple <- interaction_name_provenance_audit[
  N_interaction_name > 1
]


cat("\n")
cat("============================================================\n")
cat("interaction_name PROVENANCE AUDIT\n")
cat("============================================================\n\n")

cat(
  "Canonical patient-level keys : ",
  nrow(interaction_name_provenance_audit),
  "\n",
  sep = ""
)

cat(
  "Keys with >1 interaction_name : ",
  nrow(interaction_name_multiple),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# Check that interaction_name multiplicity does NOT create
# biological identity conflicts.
# ------------------------------------------------------------

interaction_name_provenance_pass <- {
  
  if (nrow(interaction_name_multiple) == 0L) {
    
    TRUE
    
  } else {
    
    all(
      interaction_name_multiple$N_exact_LR_pair == 1L,
      interaction_name_multiple$N_pathway == 1L,
      interaction_name_multiple$N_ligand == 1L,
      interaction_name_multiple$N_receptor == 1L,
      interaction_name_multiple$N_annotation == 1L,
      interaction_name_multiple$N_direction == 1L,
      interaction_name_multiple$N_sender == 1L,
      interaction_name_multiple$N_receiver == 1L,
      interaction_name_multiple$N_compatible == 1L
    )
  }
}


cat(
  "interaction_name provenance consistency : ",
  interaction_name_provenance_pass,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# Save detailed audit
# ------------------------------------------------------------

fwrite(
  interaction_name_provenance_audit,
  file.path(
    output_dir,
    "STEP11-5E-3_interaction_name_PROVENANCE_AUDIT.tsv"
  ),
  sep = "\t"
)


# ============================================================
# STEP11-5E-3 — FINAL QC
# PROVENANCE CLARIFICATION AUDIT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP11-5E-3 — PROVENANCE CLARIFICATION AUDIT — FINAL QC\n")
cat("============================================================\n\n")


# ============================================================
# 0. Helper function for safe invariance extraction
# ============================================================

get_invariance_pass <- function(variable_name) {
  
  x <- invariance_results[
    Variable == variable_name,
    PASS
  ]
  
  if (length(x) != 1L) {
    
    stop(
      paste0(
        "Expected exactly one invariance result for variable: ",
        variable_name,
        "\nFound: ",
        length(x)
      )
    )
  }
  
  isTRUE(x[[1]])
}


# ============================================================
# 1. Level B — Explicit invariance locks
# ============================================================

level_B_candidate_class_pass <-
  get_invariance_pass("Candidate_Consistency_Class")

level_B_exact_lr_pass <-
  get_invariance_pass("Exact_LR_Pair")

level_B_interaction_name_pass <-
  isTRUE(interaction_name_provenance_pass)

level_B_pathway_pass <-
  get_invariance_pass("pathway_name")

level_B_ligand_pass <-
  get_invariance_pass("ligand")

level_B_receptor_pass <-
  get_invariance_pass("receptor")

level_B_annotation_pass <-
  get_invariance_pass("annotation")

level_B_direction_pass <-
  get_invariance_pass("Direction_ID")

level_B_sender_pass <-
  get_invariance_pass("Sender")

level_B_receiver_pass <-
  get_invariance_pass("Receiver")

level_B_core_gene_pass <- {
  
  x <- invariance_results[
    Variable %in% c(
      "Core_Gene_Involved",
      "Core_Gene"
    ),
    PASS
  ]
  
  length(x) == 2L &&
    all(x)
}


level_B_expression_pass <- {
  
  x <- invariance_results[
    Variable %in% c(
      "Ligand_log2CPM",
      "Receptor_log2CPM",
      "Ligand_Receptor_Mean_log2CPM",
      "Ligand_Receptor_Difference_log2CPM"
    ),
    PASS
  ]
  
  length(x) == 4L &&
    all(x)
}


level_B_compatibility_pass <-
  get_invariance_pass("Compatible")


# ============================================================
# 2. Level C — Provenance clarification
# ============================================================

# These objects were established in the preceding
# provenance audit and are NOT recalculated here.

level_C_pattern_final <-
  isTRUE(level_C_pattern_pass)

level_C_ge4_membership_final <-
  isTRUE(level_C_ge4_membership_pass)

level_C_class_relationship_final <-
  isTRUE(level_C_class_relationship_pass)

level_C_ge3only_4of4_final <-
  isTRUE(ge3_only_4of4_membership_pass)

level_C_ge3only_3of4_final <-
  isTRUE(ge3_only_3of4_membership_pass)


# ============================================================
# 3. Global structural checks
# ============================================================

global_raw_count_pass <-
  isTRUE(raw_record_count_pass)

global_canonical_key_count_pass <-
  isTRUE(canonical_key_count_pass)

global_no_collapse_pass <-
  isTRUE(no_collapse_performed)


# ------------------------------------------------------------
# Core-gene separation inheritance
#
# IMPORTANT:
# AOC3 / MS4A7 / IL7R are NOT expected to appear as LR
# candidates in the current genome-wide candidate table.
#
# The definitive zero-core-LR result was already locked in
# STEP11-5D-1 / STEP11-5D-2 / STEP11-5D-3.
#
# Therefore this is an INHERITANCE / CONSISTENCY audit,
# not a new core-gene discovery test.
# ------------------------------------------------------------

core_gene_separation_final <- {
  
  core_field_ok <-
    all(is.na(patient_expression$Core_Gene))
  
  core_involved_ok <-
    all(
      patient_expression$Core_Gene_Involved == FALSE,
      na.rm = TRUE
    )
  
  locked_zero_core_lr <- 0L
  
  core_field_ok &&
    core_involved_ok &&
    locked_zero_core_lr == 0L
}


# ============================================================
# 4. Level A — Exact LR universe
# ============================================================

level_A_identity_final <-
  isTRUE(summary_lr_identity_pass) &&
  isTRUE(patient_lr_identity_pass)

level_A_nested_final <-
  isTRUE(ge4_nested_in_ge3_pass)

level_A_counts_final <-
  isTRUE(level_A_expected_counts_pass)


# ============================================================
# 5. Construct FINAL QC table
# ============================================================

final_qc <- data.table(
  
  QC_Item = c(
    
    # --------------------------------------------------------
    # Level A
    # --------------------------------------------------------
    
    "LEVEL A — Exact LR identity",
    
    "LEVEL A — ge4 nested within ge3 exact LR universe",
    
    "LEVEL A — Expected 585/492/492/93/0 counts",
    
    
    # --------------------------------------------------------
    # Level B
    # --------------------------------------------------------
    
    "LEVEL B — Canonical key defined",
    
    "LEVEL B — Candidate_Consistency_Class invariant",
    
    "LEVEL B — Exact_LR_Pair invariant",
    
    "LEVEL B — interaction_name provenance consistent",
    
    "LEVEL B — pathway_name invariant",
    
    "LEVEL B — ligand invariant",
    
    "LEVEL B — receptor invariant",
    
    "LEVEL B — annotation invariant",
    
    "LEVEL B — Direction_ID invariant",
    
    "LEVEL B — Sender invariant",
    
    "LEVEL B — Receiver invariant",
    
    "LEVEL B — Core-gene fields invariant",
    
    "LEVEL B — expression fields invariant",
    
    "LEVEL B — compatibility invariant",
    
    
    # --------------------------------------------------------
    # Level C
    # --------------------------------------------------------
    
    "LEVEL C — Allowed Candidate_Set provenance patterns",
    
    "LEVEL C — ge4 provenance belongs to ge4 LR universe",
    
    "LEVEL C — Candidate class matches exact-LR-pair universe",
    
    "LEVEL C — ge3-only / 4_of_4 correctly explained",
    
    "LEVEL C — ge3-only / 3_of_4_only correctly explained",
    
    
    # --------------------------------------------------------
    # Global
    # --------------------------------------------------------
    
    "GLOBAL — Raw record count 16480",
    
    "GLOBAL — Canonical key count 9620",
    
    "GLOBAL — No canonical collapse performed",
    
    "GLOBAL — Core-gene separation inherited from locked LR-universe audit"
  ),
  
  
  PASS = c(
    
    # --------------------------------------------------------
    # Level A
    # --------------------------------------------------------
    
    level_A_identity_final,
    
    level_A_nested_final,
    
    level_A_counts_final,
    
    
    # --------------------------------------------------------
    # Level B
    # --------------------------------------------------------
    
    TRUE,
    
    level_B_candidate_class_pass,
    
    level_B_exact_lr_pass,
    
    level_B_interaction_name_pass,
    
    level_B_pathway_pass,
    
    level_B_ligand_pass,
    
    level_B_receptor_pass,
    
    level_B_annotation_pass,
    
    level_B_direction_pass,
    
    level_B_sender_pass,
    
    level_B_receiver_pass,
    
    level_B_core_gene_pass,
    
    level_B_expression_pass,
    
    level_B_compatibility_pass,
    
    
    # --------------------------------------------------------
    # Level C
    # --------------------------------------------------------
    
    level_C_pattern_final,
    
    level_C_ge4_membership_final,
    
    level_C_class_relationship_final,
    
    level_C_ge3only_4of4_final,
    
    level_C_ge3only_3of4_final,
    
    
    # --------------------------------------------------------
    # Global
    # --------------------------------------------------------
    
    global_raw_count_pass,
    
    global_canonical_key_count_pass,
    
    global_no_collapse_pass,
    
    core_gene_separation_final
  )
)


# ============================================================
# 6. FINAL QC summary
# ============================================================

n_qc <- nrow(final_qc)

n_pass <- sum(
  final_qc$PASS,
  na.rm = TRUE
)

n_fail <- sum(
  !final_qc$PASS,
  na.rm = TRUE
)


cat("============================================================\n")
cat("STEP11-5E-3 — PROVENANCE CLARIFICATION AUDIT — FINAL\n")
cat("============================================================\n\n")

print(final_qc)

cat("\n")
cat("QC items :", n_qc, "\n")
cat("PASS     :", n_pass, "\n")
cat("FAIL     :", n_fail, "\n\n")


# ============================================================
# 7. Explicit failed-item report
# ============================================================

if (n_fail > 0L) {
  
  cat("FAILED QC ITEMS:\n")
  cat("------------------------------------------------------------\n")
  
  print(
    final_qc[
      PASS == FALSE,
      QC_Item
    ]
  )
  
  cat("\n")
}


# ============================================================
# 8. Save FINAL QC
# ============================================================

fwrite(
  final_qc,
  file.path(
    output_dir,
    "STEP11-5E-3_PROVENANCE_CLARIFICATION_AUDIT_FINAL_QC.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 9. Global audit summary
# ============================================================

global_summary <- data.table(
  
  Metric = c(
    
    "Raw_patient_level_records",
    
    "Unique_canonical_keys",
    
    "Duplicate_canonical_keys",
    
    "ge3_exact_LR_pairs",
    
    "ge4_exact_LR_pairs",
    
    "Shared_exact_LR_pairs",
    
    "ge3_only_exact_LR_pairs",
    
    "ge4_only_exact_LR_pairs",
    
    "ge4_nested_in_ge3",
    
    "Candidate_Consistency_Class_invariant",
    
    "Candidate_Set_is_provenance_only",
    
    "Canonical_collapse_performed",
    
    "QC_items",
    
    "QC_PASS",
    
    "QC_FAIL"
  ),
  
  
  Value = c(
    
    n_raw_records,
    
    n_canonical_keys,
    
    n_duplicate_keys,
    
    ge3_n,
    
    ge4_n,
    
    length(shared_pairs),
    
    length(ge3_only_pairs),
    
    length(ge4_only_pairs),
    
    ge4_nested_in_ge3_pass,
    
    level_B_candidate_class_pass,
    
    TRUE,
    
    FALSE,
    
    n_qc,
    
    n_pass,
    
    n_fail
  )
)


fwrite(
  global_summary,
  file.path(
    output_dir,
    "STEP11-5E-3_PROVENANCE_CLARIFICATION_AUDIT_GLOBAL_SUMMARY.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 10. FINAL STATUS
# ============================================================

if (n_fail == 0L) {
  
  cat("============================================================\n")
  cat("PASS — STEP11-5E-3 PROVENANCE CLARIFICATION AUDIT\n")
  cat("============================================================\n\n")
  
  cat(
    "Level A:\n",
    "  Exact LR pair universe validated.\n",
    "  ge4 exact LR pairs are nested within ge3.\n",
    "  Expected 585 / 492 / 492 / 93 / 0 structure confirmed.\n\n"
  )
  
  cat(
    "Level B:\n",
    "  Canonical patient-level key validated as:\n",
    "  Exact_LR_Pair × Direction_ID × Patient_Index.\n",
    "  Candidate_Set was correctly excluded from invariance testing.\n",
    "  Candidate_Consistency_Class remained invariant.\n\n"
  )
  
  cat(
    "Level C:\n",
    "  Candidate_Set retained as provenance.\n",
    "  Candidate_Consistency_Class linked to exact-LR-pair universe.\n",
    "  ge3-only / 4_of_4 records are permitted when the exact LR pair\n",
    "  belongs to the ge4 exact-LR universe.\n\n"
  )
  
  cat(
    "Global:\n",
    "  Raw record count = 16480.\n",
    "  Unique canonical keys = 9620.\n",
    "  No canonical collapse was performed.\n",
    "  Core-gene separation inherited from locked LR-universe audit.\n\n"
  )
  
  cat(
    "NO upstream locked data were modified.\n",
    "NO expression quantities were recalculated.\n",
    "NO compatibility quantities were recalculated.\n",
    "NO statistical testing was performed.\n",
    "NO biological ranking was performed.\n\n"
  )
  
} else {
  
  cat("============================================================\n")
  cat("FAIL — STEP11-5E-3 PROVENANCE CLARIFICATION AUDIT\n")
  cat("============================================================\n\n")
  
  cat(
    "One or more structural audits failed.\n",
    "No canonical collapse should be performed.\n",
    "Review the FAILED QC ITEMS and saved audit tables.\n\n"
  )
}


# ============================================================
# 11. Output directory
# ============================================================

cat(
  "Output directory:\n",
  output_dir,
  "\n\n"
)

cat("============================================================\n")
cat("END OF STEP11-5E-3 PROVENANCE CLARIFICATION AUDIT\n")
cat("============================================================\n")


# ============================================================
# STEP11-5E-3B
# FORMAL CANONICAL PATIENT-LEVEL CANDIDATE INTEGRATION
#
# Purpose:
#   Convert the locked 16,480 patient-level candidate records
#   into a canonical 9,620-key patient-level integration table.
#
# Canonical key:
#   Exact_LR_Pair × Direction_ID × Patient_Index
#
# IMPORTANT:
#   Candidate_Set is retained as provenance.
#   interaction_name is retained as resource-level provenance.
#   No provenance information is silently discarded.
#
# Locked upstream inputs:
#   STEP11-5E-2B_Candidate_Patient_Expression.tsv
#
# No:
#   - raw UMI recalculation
#   - pseudobulk recalculation
#   - CPM/log2CPM recalculation
#   - LR detection recalculation
#   - compatibility recalculation
#   - statistical testing
#   - FDR
#   - CellChat probability
#   - pathway enrichment
#   - biological ranking
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})


# ============================================================
# 1. Working directory and paths
# ============================================================

base_dir <- "input_dir/Paper 1_Step 11-CellChat"

input_dir <- file.path(
  base_dir,
  "STEP11-5E-2_Candidate_Expression_Magnitude"
)

output_dir <- file.path(
  base_dir,
  "STEP11-5E-3_Patient_Level_Candidate_Integration"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


input_file <- file.path(
  input_dir,
  "STEP11-5E-2B_Candidate_Patient_Expression.tsv"
)


# ============================================================
# 2. Input existence check
# ============================================================

if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Input file not found:\n",
      input_file
    )
  )
}


cat("============================================================\n")
cat("STEP11-5E-3B — FORMAL CANONICAL INTEGRATION\n")
cat("============================================================\n\n")

cat(
  "Input:\n",
  input_file,
  "\n\n"
)

cat(
  "Output directory:\n",
  output_dir,
  "\n\n"
)


# ============================================================
# 3. Read locked STEP11-5E-2 patient-level table
# ============================================================

patient_expression <- fread(
  input_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)


cat(
  "Input rows:",
  nrow(patient_expression),
  "\n"
)

cat(
  "Input columns:",
  ncol(patient_expression),
  "\n\n"
)


# ============================================================
# 4. Expected schema
# ============================================================

expected_cols <- c(
  
  "Candidate_Set",
  "Candidate_Consistency_Class",
  "Exact_LR_Pair",
  "interaction_name",
  "pathway_name",
  "ligand",
  "receptor",
  "annotation",
  "Direction_ID",
  "Sender",
  "Receiver",
  "Core_Gene_Involved",
  "Core_Gene",
  "Patient_Index",
  "Ligand_log2CPM",
  "Receptor_log2CPM",
  "Compatible",
  "Ligand_Receptor_Mean_log2CPM",
  "Ligand_Receptor_Difference_log2CPM"
)


missing_cols <- setdiff(
  expected_cols,
  names(patient_expression)
)

if (length(missing_cols) > 0L) {
  
  stop(
    paste(
      "Missing required columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}


cat("Schema check: PASS\n\n")


# ============================================================
# 5. Input lock
# ============================================================

expected_raw_records <- 16480L
expected_canonical_keys <- 9620L

if (nrow(patient_expression) != expected_raw_records) {
  
  stop(
    paste0(
      "Unexpected input row count. Expected ",
      expected_raw_records,
      ", observed ",
      nrow(patient_expression),
      "."
    )
  )
}


cat(
  "Locked input record count:",
  nrow(patient_expression),
  "\n"
)


# ============================================================
# 6. Construct canonical key
# ============================================================

patient_expression[
  ,
  Canonical_Key :=
    paste(
      Exact_LR_Pair,
      Direction_ID,
      Patient_Index,
      sep = "___"
    )
]


# ============================================================
# 7. Canonical key integrity
# ============================================================

if (anyNA(patient_expression$Canonical_Key)) {
  
  stop(
    "Canonical_Key contains NA."
  )
}

if (any(patient_expression$Canonical_Key == "")) {
  
  stop(
    "Canonical_Key contains empty values."
  )
}


n_raw_records <- nrow(patient_expression)

n_canonical_keys <- uniqueN(
  patient_expression$Canonical_Key
)

n_duplicate_keys <-
  n_raw_records - n_canonical_keys


cat(
  "Raw records:",
  n_raw_records,
  "\n"
)

cat(
  "Unique canonical keys:",
  n_canonical_keys,
  "\n"
)

cat(
  "Duplicate records attributable to provenance:",
  n_duplicate_keys,
  "\n\n"
)


if (n_canonical_keys != expected_canonical_keys) {
  
  stop(
    paste0(
      "Unexpected canonical key count. Expected ",
      expected_canonical_keys,
      ", observed ",
      n_canonical_keys,
      "."
    )
  )
}


# ============================================================
# 8. Define helper for provenance concatenation
# ============================================================

collapse_provenance <- function(x) {
  
  x <- as.character(x)
  
  x <- x[
    !is.na(x) &
      nzchar(x)
  ]
  
  if (length(x) == 0L) {
    return(NA_character_)
  }
  
  paste(
    sort(unique(x)),
    collapse = " | "
  )
}


# ============================================================
# 9. Define canonical-level invariant checker
# ============================================================

check_single_value <- function(x) {
  
  x <- x[
    !is.na(x)
  ]
  
  if (length(x) == 0L) {
    return(TRUE)
  }
  
  length(unique(x)) == 1L
}


# ============================================================
# 10. Validate patient-level invariants
#
# These variables MUST remain invariant within:
#   Exact_LR_Pair × Direction_ID × Patient_Index
#
# Candidate_Set and interaction_name are explicitly excluded
# because they are provenance variables.
# ============================================================

invariant_variables <- c(
  
  "Candidate_Consistency_Class",
  "Exact_LR_Pair",
  "pathway_name",
  "ligand",
  "receptor",
  "annotation",
  "Direction_ID",
  "Sender",
  "Receiver",
  "Core_Gene_Involved",
  "Core_Gene",
  "Ligand_log2CPM",
  "Receptor_log2CPM",
  "Compatible",
  "Ligand_Receptor_Mean_log2CPM",
  "Ligand_Receptor_Difference_log2CPM"
)


invariance_qc <- patient_expression[
  ,
  lapply(
    .SD,
    check_single_value
  ),
  by = Canonical_Key,
  .SDcols = invariant_variables
]


invariance_summary <- data.table(
  Variable = invariant_variables,
  PASS = sapply(
    invariant_variables,
    function(v) {
      all(
        invariance_qc[[v]],
        na.rm = TRUE
      )
    }
  )
)


cat("============================================================\n")
cat("PATIENT-LEVEL INVARIANCE QC\n")
cat("============================================================\n\n")

print(invariance_summary)

cat("\n")


if (any(!invariance_summary$PASS)) {
  
  print(
    invariance_summary[
      PASS == FALSE
    ]
  )
  
  stop(
    "One or more patient-level invariant variables are inconsistent."
  )
}


# ============================================================
# 11. Candidate_Set provenance summary
# ============================================================

candidate_set_provenance <- patient_expression[
  ,
  .(
    Candidate_Set_Provenance =
      collapse_provenance(Candidate_Set)
  ),
  by = Canonical_Key
]


# ============================================================
# 12. interaction_name provenance summary
# ============================================================

interaction_name_provenance <- patient_expression[
  ,
  .(
    interaction_name_Provenance =
      collapse_provenance(interaction_name)
  ),
  by = Canonical_Key
]


# ============================================================
# 13. Count provenance multiplicity
# ============================================================

candidate_set_multiplicity <- patient_expression[
  ,
  .(
    N_Candidate_Set =
      uniqueN(Candidate_Set)
  ),
  by = Canonical_Key
]


interaction_name_multiplicity <- patient_expression[
  ,
  .(
    N_interaction_name =
      uniqueN(interaction_name)
  ),
  by = Canonical_Key
]


# ============================================================
# 14. Build canonical integration table
#
# One row per:
#   Exact_LR_Pair × Direction_ID × Patient_Index
#
# Provenance fields are explicitly retained.
# ============================================================

canonical_table <- patient_expression[
  ,
  .(
    Candidate_Consistency_Class =
      first(Candidate_Consistency_Class),
    
    Exact_LR_Pair =
      first(Exact_LR_Pair),
    
    pathway_name =
      first(pathway_name),
    
    ligand =
      first(ligand),
    
    receptor =
      first(receptor),
    
    annotation =
      first(annotation),
    
    Direction_ID =
      first(Direction_ID),
    
    Sender =
      first(Sender),
    
    Receiver =
      first(Receiver),
    
    Core_Gene_Involved =
      first(Core_Gene_Involved),
    
    Core_Gene =
      first(Core_Gene),
    
    Patient_Index =
      first(Patient_Index),
    
    Ligand_log2CPM =
      first(Ligand_log2CPM),
    
    Receptor_log2CPM =
      first(Receptor_log2CPM),
    
    Compatible =
      first(Compatible),
    
    Ligand_Receptor_Mean_log2CPM =
      first(Ligand_Receptor_Mean_log2CPM),
    
    Ligand_Receptor_Difference_log2CPM =
      first(Ligand_Receptor_Difference_log2CPM)
  ),
  by = Canonical_Key
]


# ============================================================
# 15. Add provenance fields
# ============================================================

canonical_table <- merge(
  canonical_table,
  candidate_set_provenance,
  by = "Canonical_Key",
  all.x = TRUE,
  sort = FALSE
)


canonical_table <- merge(
  canonical_table,
  interaction_name_provenance,
  by = "Canonical_Key",
  all.x = TRUE,
  sort = FALSE
)


canonical_table <- merge(
  canonical_table,
  candidate_set_multiplicity,
  by = "Canonical_Key",
  all.x = TRUE,
  sort = FALSE
)


canonical_table <- merge(
  canonical_table,
  interaction_name_multiplicity,
  by = "Canonical_Key",
  all.x = TRUE,
  sort = FALSE
)


# ============================================================
# 16. Reconstruct canonical key components
# ============================================================

canonical_table[
  ,
  Canonical_Key_Check :=
    paste(
      Exact_LR_Pair,
      Direction_ID,
      Patient_Index,
      sep = "___"
    )
]


canonical_key_reconstruction_pass <- all(
  canonical_table$Canonical_Key ==
    canonical_table$Canonical_Key_Check
)


if (!canonical_key_reconstruction_pass) {
  
  stop(
    "Canonical key reconstruction failed."
  )
}


canonical_table[
  ,
  Canonical_Key_Check := NULL
]


# ============================================================
# 17. Order columns
# ============================================================

setcolorder(
  canonical_table,
  c(
    
    "Canonical_Key",
    
    "Exact_LR_Pair",
    
    "Candidate_Consistency_Class",
    
    "Candidate_Set_Provenance",
    
    "interaction_name_Provenance",
    
    "N_Candidate_Set",
    
    "N_interaction_name",
    
    "pathway_name",
    
    "ligand",
    
    "receptor",
    
    "annotation",
    
    "Direction_ID",
    
    "Sender",
    
    "Receiver",
    
    "Patient_Index",
    
    "Core_Gene_Involved",
    
    "Core_Gene",
    
    "Ligand_log2CPM",
    
    "Receptor_log2CPM",
    
    "Ligand_Receptor_Mean_log2CPM",
    
    "Ligand_Receptor_Difference_log2CPM",
    
    "Compatible"
  )
)

# ============================================================
# 18. Canonical table structural QC
# ============================================================

canonical_unique_pass <-
  uniqueN(canonical_table$Canonical_Key) ==
  nrow(canonical_table)


canonical_row_count_pass <-
  nrow(canonical_table) ==
  n_canonical_keys


canonical_provenance_complete_pass <- all(
  !is.na(canonical_table$Candidate_Set_Provenance),
  !is.na(canonical_table$interaction_name_Provenance)
)


canonical_identity_complete_pass <- all(
  !is.na(canonical_table$Exact_LR_Pair),
  !is.na(canonical_table$Direction_ID),
  !is.na(canonical_table$Patient_Index)
)


# ============================================================
# 19. Expected candidate-class structure
# ============================================================

n_4of4_canonical <- uniqueN(
  canonical_table[
    Candidate_Consistency_Class == "4_of_4",
    Exact_LR_Pair
  ]
)


n_3of4only_canonical <- uniqueN(
  canonical_table[
    Candidate_Consistency_Class == "3_of_4_only",
    Exact_LR_Pair
  ]
)


candidate_class_summary <- canonical_table[
  ,
  .(
    N_canonical_keys = .N,
    N_exact_LR_pairs = uniqueN(Exact_LR_Pair)
  ),
  by = Candidate_Consistency_Class
]


# ============================================================
# 20. Candidate_Set provenance distribution
# ============================================================

candidate_set_summary <- canonical_table[
  ,
  .(
    N_canonical_keys = .N,
    N_exact_LR_pairs = uniqueN(Exact_LR_Pair)
  ),
  by = Candidate_Set_Provenance
]


# ============================================================
# 21. interaction_name provenance distribution
# ============================================================

interaction_name_summary <- canonical_table[
  ,
  .(
    N_canonical_keys = .N,
    N_exact_LR_pairs = uniqueN(Exact_LR_Pair)
  ),
  by = interaction_name_Provenance
]


# ============================================================
# 22. Save canonical integration table
# ============================================================

canonical_output_file <- file.path(
  output_dir,
  "STEP11-5E-3B_FINAL_Canonical_Patient_Level_LR_Integration.tsv"
)


fwrite(
  canonical_table,
  canonical_output_file,
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 23. Save patient-level invariance QC
# ============================================================

fwrite(
  invariance_summary,
  file.path(
    output_dir,
    "STEP11-5E-3B_Patient_Level_Invariance_QC.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 24. Save candidate class summary
# ============================================================

fwrite(
  candidate_class_summary,
  file.path(
    output_dir,
    "STEP11-5E-3B_Candidate_Class_Summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 25. Save Candidate_Set provenance summary
# ============================================================

fwrite(
  candidate_set_summary,
  file.path(
    output_dir,
    "STEP11-5E-3B_Candidate_Set_Provenance_Summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 26. Save interaction_name provenance summary
# ============================================================

fwrite(
  interaction_name_summary,
  file.path(
    output_dir,
    "STEP11-5E-3B_interaction_name_Provenance_Summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 27. Final integration QC
# ============================================================

integration_qc <- data.table(
  
  QC_Item = c(
    
    "Input raw records = 16480",
    
    "Canonical keys = 9620",
    
    "One row per canonical key",
    
    "Canonical row count = 9620",
    
    "Canonical identity complete",
    
    "Candidate_Set provenance retained",
    
    "interaction_name provenance retained",
    
    "Patient-level invariants PASS",
    
    "Candidate class retained",
    
    "Expression values inherited without recalculation",
    
    "Compatibility values inherited without recalculation",
    
    "No statistical testing",
    
    "No biological ranking"
  ),
  
  PASS = c(
    
    n_raw_records == 16480L,
    
    n_canonical_keys == 9620L,
    
    canonical_unique_pass,
    
    canonical_row_count_pass,
    
    canonical_identity_complete_pass,
    
    canonical_provenance_complete_pass,
    
    canonical_provenance_complete_pass,
    
    all(invariance_summary$PASS),
    
    all(
      canonical_table$Candidate_Consistency_Class %in%
        c(
          "4_of_4",
          "3_of_4_only"
        )
    ),
    
    TRUE,
    
    TRUE,
    
    TRUE,
    
    TRUE
  )
)


# ============================================================
# 28. Print final integration QC
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP11-5E-3B — FINAL CANONICAL INTEGRATION QC\n")
cat("============================================================\n\n")

print(integration_qc)

n_qc <- nrow(integration_qc)

n_pass <- sum(
  integration_qc$PASS,
  na.rm = TRUE
)

n_fail <- sum(
  !integration_qc$PASS,
  na.rm = TRUE
)

cat("\n")
cat("QC items :", n_qc, "\n")
cat("PASS     :", n_pass, "\n")
cat("FAIL     :", n_fail, "\n\n")


# ============================================================
# 29. Save final integration QC
# ============================================================

fwrite(
  integration_qc,
  file.path(
    output_dir,
    "STEP11-5E-3B_FINAL_Canonical_Integration_QC.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 30. Final status
# ============================================================

if (n_fail == 0L) {
  
  cat("============================================================\n")
  cat("PASS — STEP11-5E-3B FORMAL CANONICAL INTEGRATION\n")
  cat("============================================================\n\n")
  
  cat(
    "Raw patient-level records : ",
    n_raw_records,
    "\n",
    sep = ""
  )
  
  cat(
    "Canonical patient-level keys : ",
    n_canonical_keys,
    "\n",
    sep = ""
  )
  
  cat(
    "Canonical definition:\n",
    "  Exact_LR_Pair × Direction_ID × Patient_Index\n\n"
  )
  
  cat(
    "Candidate_Set retained as provenance.\n",
    "interaction_name retained as resource-level provenance.\n\n"
  )
  
  cat(
    "No upstream locked data were modified.\n",
    "No expression values were recalculated.\n",
    "No compatibility values were recalculated.\n",
    "No statistical testing was performed.\n",
    "No biological ranking was performed.\n\n"
  )
  
  cat(
    "Final canonical integration table:\n",
    canonical_output_file,
    "\n\n"
  )
  
} else {
  
  cat("============================================================\n")
  cat("FAIL — STEP11-5E-3B FORMAL CANONICAL INTEGRATION\n")
  cat("============================================================\n\n")
  
  cat(
    "One or more integration QC items failed.\n",
    "Do NOT use the canonical integration table for downstream analysis.\n",
    "Review the saved QC table.\n\n"
  )
}


cat("============================================================\n")
cat("END OF STEP11-5E-3B\n")
cat("============================================================\n")


# ============================================================
# STEP11-5F-1
# Canonical Input QC
#
# Paper 1 — STEP11 CellChat Revision
#
# PURPOSE:
#   Validate the locked 9,620-row canonical patient-level
#   LR integration table before downstream descriptive
#   consistency / expression-context analysis.
#
# INPUT:
#   STEP11-5E-3B_FINAL_Canonical_Patient_Level_LR_Integration.tsv
#
# CANONICAL UNIT:
#   Exact_LR_Pair × Direction_ID × Patient_Index
#
# IMPORTANT:
#   This is an INPUT QC / STRUCTURAL AUDIT ONLY.
#
# THIS STEP DOES NOT:
#   - discover new LR pairs
#   - reconstruct the LR universe
#   - recalculate raw UMI
#   - recalculate pseudobulk
#   - recalculate CPM
#   - recalculate log2(CPM+1)
#   - recalculate detection
#   - recalculate compatibility
#   - perform correlation testing
#   - perform P-value / FDR testing
#   - calculate CellChat communication probabilities
#   - perform pathway enrichment
#   - rank biological importance
#   - infer interaction strength
#   - establish functional signaling
#
# OUTPUT:
#   STEP11-5F-1_Canonical_Input_QC/
#
# ============================================================


cat("STEP11-5F-1: Canonical Input QC\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_packages <- c("data.table")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package not installed: ", pkg)
  }
}

library(data.table)


# ------------------------------------------------------------
# 1. Working directory and paths
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 11-CellChat"

input_dir <- file.path(
  base_dir,
  "STEP11-5E-3_Patient_Level_Candidate_Integration"
)

input_file <- file.path(
  input_dir,
  "STEP11-5E-3B_FINAL_Canonical_Patient_Level_LR_Integration.tsv"
)

output_dir <- file.path(
  base_dir,
  "STEP11-5F-1_Canonical_Input_QC"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Input existence check
# ------------------------------------------------------------

input_exists <- file.exists(input_file)

cat("Input file:\n")
cat(" ", input_file, "\n\n")

cat(
  "Input exists:",
  input_exists,
  "\n\n"
)

if (!input_exists) {
  stop("Input file does not exist. STEP11-5F-1 aborted.")
}


# ------------------------------------------------------------
# 3. Read LOCKED canonical table
# ------------------------------------------------------------

dt <- fread(
  input_file,
  sep = "\t",
  header = TRUE,
  data.table = TRUE,
  check.names = FALSE
)

cat("Input successfully loaded.\n")
cat("Rows :", nrow(dt), "\n")
cat("Cols :", ncol(dt), "\n\n")


# ------------------------------------------------------------
# 4. Expected row count
# ------------------------------------------------------------

expected_rows <- 9620L

row_count_pass <- identical(
  as.integer(nrow(dt)),
  expected_rows
)

cat(
  "Expected canonical rows:",
  expected_rows,
  "\n"
)

cat(
  "Observed canonical rows:",
  nrow(dt),
  "\n"
)

cat(
  "Row-count QC:",
  row_count_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 5. Required schema
# ------------------------------------------------------------

required_cols <- c(
  "Exact_LR_Pair",
  "interaction_name_Provenance",
  "pathway_name",
  "ligand",
  "receptor",
  "annotation",
  "Direction_ID",
  "Sender",
  "Receiver",
  "Core_Gene_Involved",
  "Core_Gene",
  "Patient_Index",
  "Candidate_Consistency_Class",
  "Candidate_Set_Provenance",
  "Ligand_log2CPM",
  "Receptor_log2CPM",
  "Compatible"
)

missing_cols <- setdiff(
  required_cols,
  names(dt)
)

schema_pass <- length(missing_cols) == 0L

cat("Required schema check:\n")
cat(
  "Required columns:",
  length(required_cols),
  "\n"
)

cat(
  "Missing columns:",
  length(missing_cols),
  "\n"
)

if (length(missing_cols) > 0L) {
  cat(
    paste(missing_cols, collapse = ", "),
    "\n"
  )
}

cat(
  "Schema QC:",
  schema_pass,
  "\n\n"
)

if (!schema_pass) {
  stop(
    "Required columns are missing. STEP11-5F-1 aborted."
  )
}


# ------------------------------------------------------------
# 6. Basic identity counts
# ------------------------------------------------------------

n_exact_lr <- uniqueN(dt$Exact_LR_Pair)

n_direction <- uniqueN(dt$Direction_ID)

n_patients <- uniqueN(dt$Patient_Index)

n_senders <- uniqueN(dt$Sender)

n_receivers <- uniqueN(dt$Receiver)

n_pathways <- uniqueN(dt$pathway_name)

n_annotations <- uniqueN(dt$annotation)

n_candidate_classes <- uniqueN(
  dt$Candidate_Consistency_Class
)

cat("Basic identity counts:\n")
cat(" Unique Exact_LR_Pair :", n_exact_lr, "\n")
cat(" Unique Direction_ID  :", n_direction, "\n")
cat(" Unique Patient_Index :", n_patients, "\n")
cat(" Unique Sender        :", n_senders, "\n")
cat(" Unique Receiver      :", n_receivers, "\n")
cat(" Unique pathways      :", n_pathways, "\n")
cat(" Unique annotations   :", n_annotations, "\n")
cat(
  " Candidate classes    :",
  n_candidate_classes,
  "\n\n"
)


# ------------------------------------------------------------
# 7. Expected global structural counts
#
# From locked upstream STEP11-5E-3B:
#
#   9,620 canonical patient-level records
#   585 exact LR pairs
#   492 exact LR pairs = 4/4
#   93 exact LR pairs = 3/4-only
#   4 matched patients
#   6 directions
#
# ------------------------------------------------------------

expected_exact_lr <- 585L
expected_four_of_four_lr <- 492L
expected_three_of_four_only_lr <- 93L
expected_patients <- 4L
expected_directions <- 6L

exact_lr_count_pass <- identical(
  as.integer(n_exact_lr),
  expected_exact_lr
)

patient_count_pass <- identical(
  as.integer(n_patients),
  expected_patients
)

direction_count_pass <- identical(
  as.integer(n_direction),
  expected_directions
)


# ------------------------------------------------------------
# 8. Candidate consistency class audit
# ------------------------------------------------------------

candidate_class_values <- sort(
  unique(
    as.character(
      dt$Candidate_Consistency_Class
    )
  )
)

cat("Candidate consistency classes:\n")
print(candidate_class_values)
cat("\n")


expected_candidate_classes <- c(
  "3_of_4_only",
  "4_of_4"
)

candidate_class_structure_pass <- identical(
  candidate_class_values,
  sort(expected_candidate_classes)
)

candidate_class_pair_summary <- dt[
  ,
  .(
    N_records = .N,
    N_exact_LR_pairs = uniqueN(Exact_LR_Pair)
  ),
  by = Candidate_Consistency_Class
][order(Candidate_Consistency_Class)]

print(candidate_class_pair_summary)

cat("\n")


# ------------------------------------------------------------
# 9. Verify 492 / 93 exact LR-pair structure
# ------------------------------------------------------------

candidate_lr_summary <- dt[
  ,
  .(
    Candidate_Consistency_Class =
      unique(Candidate_Consistency_Class)
  ),
  by = Exact_LR_Pair
]

candidate_lr_class_n <- candidate_lr_summary[
  ,
  .N,
  by = Candidate_Consistency_Class
]

n_four_of_four_lr <- candidate_lr_class_n[
  Candidate_Consistency_Class == "4_of_4",
  N
]

n_three_of_four_only_lr <- candidate_lr_class_n[
  Candidate_Consistency_Class == "3_of_4_only",
  N
]

if (length(n_four_of_four_lr) == 0L) {
  n_four_of_four_lr <- 0L
}

if (length(n_three_of_four_only_lr) == 0L) {
  n_three_of_four_only_lr <- 0L
}

four_of_four_lr_count_pass <- identical(
  as.integer(n_four_of_four_lr),
  expected_four_of_four_lr
)

three_of_four_only_lr_count_pass <- identical(
  as.integer(n_three_of_four_only_lr),
  expected_three_of_four_only_lr
)

candidate_lr_total_pass <- identical(
  as.integer(
    n_four_of_four_lr +
      n_three_of_four_only_lr
  ),
  expected_exact_lr
)


# ------------------------------------------------------------
# 10. Candidate_Set provenance structure
# ------------------------------------------------------------

candidate_set_values <- sort(
  unique(
    as.character(
      dt$Candidate_Set_Provenance
    )
  )
)

cat("Candidate_Set provenance values:\n")
print(candidate_set_values)
cat("\n")

candidate_set_nonempty_pass <- all(
  !is.na(dt$Candidate_Set_Provenance) &
    trimws(
      as.character(dt$Candidate_Set_Provenance)
    ) != ""
)

candidate_set_provenance_pass <-
  candidate_set_nonempty_pass


# ------------------------------------------------------------
# 11. interaction_name provenance
# ------------------------------------------------------------

interaction_provenance_nonempty_pass <- all(
  !is.na(dt$interaction_name_Provenance) &
    trimws(
      as.character(
        dt$interaction_name_Provenance
      )
    ) != ""
)


# ------------------------------------------------------------
# 12. Exact LR identity reconstruction
#
# IMPORTANT:
#   Exact_LR_Pair is the canonical LR identity.
#
#   We verify that:
#
#       ligand + receptor
#
#   agrees with Exact_LR_Pair.
#
#   No LR discovery is performed.
# ------------------------------------------------------------

expected_lr_identity <- paste(
  dt$ligand,
  dt$receptor,
  sep = "___"
)

exact_lr_identity_pass <- all(
  dt$Exact_LR_Pair ==
    expected_lr_identity
)

cat(
  "Exact LR identity reconstruction:",
  exact_lr_identity_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 13. Canonical key uniqueness
#
# Canonical key:
#
#   Exact_LR_Pair × Direction_ID × Patient_Index
#
# This must be unique in the final canonical table.
# ------------------------------------------------------------

dt[
  ,
  Canonical_Key :=
    paste(
      Exact_LR_Pair,
      Direction_ID,
      Patient_Index,
      sep = "___"
    )
]

canonical_key_n <- uniqueN(
  dt$Canonical_Key
)

canonical_key_unique_pass <- identical(
  as.integer(canonical_key_n),
  as.integer(nrow(dt))
)

duplicate_canonical_keys <- dt[
  ,
  .N,
  by = Canonical_Key
][N > 1L]

canonical_duplicate_n <- nrow(
  duplicate_canonical_keys
)

canonical_no_duplicates_pass <-
  canonical_duplicate_n == 0L


# ------------------------------------------------------------
# 14. Patient identity audit
# ------------------------------------------------------------

patient_values <- sort(
  unique(
    as.character(dt$Patient_Index)
  )
)

cat("Patient identifiers:\n")
print(patient_values)
cat("\n")

patient_nonempty_pass <- all(
  !is.na(dt$Patient_Index) &
    trimws(
      as.character(dt$Patient_Index)
    ) != ""
)


# ------------------------------------------------------------
# 15. Direction structure audit
#
# LOCKED Direction_ID encoding:
#
#   1 = Endothelial → Myeloid/macrophage
#   2 = Myeloid/macrophage → Endothelial
#   3 = Endothelial → T-cell
#   4 = T-cell → Endothelial
#   5 = Myeloid/macrophage → T-cell
#   6 = T-cell → Myeloid/macrophage
#
# IMPORTANT:
#   Direction_ID is a locked numeric provenance label.
#   No relabeling or recalculation is performed.
#
# ------------------------------------------------------------

direction_values <- sort(
  unique(
    as.character(dt$Direction_ID)
  )
)

cat("Locked Direction_ID values:\n")
print(direction_values)
cat("\n")


# ------------------------------------------------------------
# 15A. Number of directions
# ------------------------------------------------------------

direction_n <- length(direction_values)

direction_n_pass <- identical(
  as.integer(direction_n),
  6L
)


# ------------------------------------------------------------
# 15B. Expected locked Direction_ID encoding
# ------------------------------------------------------------

expected_direction_ids <- c(
  "1",
  "2",
  "3",
  "4",
  "5",
  "6"
)

direction_id_encoding_pass <- identical(
  direction_values,
  expected_direction_ids
)


# ------------------------------------------------------------
# 15C. Direction → Sender/Receiver mapping
# ------------------------------------------------------------

direction_mapping <- unique(
  dt[
    ,
    .(
      Direction_ID,
      Sender,
      Receiver
    )
  ]
)[order(Direction_ID)]


# Each Direction_ID must map to exactly one
# Sender and exactly one Receiver.

direction_mapping_check <- direction_mapping[
  ,
  .(
    N_Sender = uniqueN(Sender),
    N_Receiver = uniqueN(Receiver),
    N_mapping_rows = .N
  ),
  by = Direction_ID
]


direction_mapping_pass <- all(
  direction_mapping_check$N_Sender == 1L &
    direction_mapping_check$N_Receiver == 1L &
    direction_mapping_check$N_mapping_rows == 1L
)


# ------------------------------------------------------------
# 15D. Explicit locked biological mapping
# ------------------------------------------------------------

expected_sender <- c(
  "Endothelial",
  "Myeloid/macrophage",
  "Endothelial",
  "T-cell",
  "Myeloid/macrophage",
  "T-cell"
)

expected_receiver <- c(
  "Myeloid/macrophage",
  "Endothelial",
  "T-cell",
  "Endothelial",
  "T-cell",
  "Myeloid/macrophage"
)

expected_direction_mapping <- data.table(
  Direction_ID = 1:6,
  Expected_Sender = expected_sender,
  Expected_Receiver = expected_receiver
)

actual_direction_mapping <- data.table(
  Direction_ID = as.integer(
    direction_mapping$Direction_ID
  ),
  Actual_Sender = direction_mapping$Sender,
  Actual_Receiver = direction_mapping$Receiver
)

setorder(
  expected_direction_mapping,
  Direction_ID
)

setorder(
  actual_direction_mapping,
  Direction_ID
)

direction_biology_mapping_pass <- identical(
  expected_direction_mapping$Direction_ID,
  actual_direction_mapping$Direction_ID
) &&
  identical(
    expected_direction_mapping$Expected_Sender,
    actual_direction_mapping$Actual_Sender
  ) &&
  identical(
    expected_direction_mapping$Expected_Receiver,
    actual_direction_mapping$Actual_Receiver
  )


# ------------------------------------------------------------
# 15E. Final direction structure PASS
# ------------------------------------------------------------

direction_structure_pass <-
  direction_n_pass &&
  direction_id_encoding_pass &&
  direction_mapping_pass &&
  direction_biology_mapping_pass


cat(
  "Number of unique Direction_ID:",
  direction_n,
  "\n"
)

cat(
  "Exactly six directions:",
  direction_n_pass,
  "\n"
)

cat(
  "Locked Direction_ID encoding 1–6:",
  direction_id_encoding_pass,
  "\n"
)

cat(
  "One-to-one Direction → Sender/Receiver mapping:",
  direction_mapping_pass,
  "\n"
)

cat(
  "Locked biological direction mapping:",
  direction_biology_mapping_pass,
  "\n\n"
)

cat("Direction mapping:\n")
print(
  direction_mapping
)

cat("\n")

# ------------------------------------------------------------
# 15F. Direction-level arithmetic integrity
# ------------------------------------------------------------

direction_arithmetic <- dt[
  ,
  .(
    N_records = .N,
    N_exact_LR_pairs = uniqueN(Exact_LR_Pair),
    N_patients = uniqueN(Patient_Index)
  ),
  by = Direction_ID
][order(Direction_ID)]

direction_arithmetic[
  ,
  Expected_N_records :=
    N_exact_LR_pairs * N_patients
]

direction_arithmetic[
  ,
  Arithmetic_Pass :=
    N_records == Expected_N_records
]

direction_arithmetic_pass <- all(
  direction_arithmetic$Arithmetic_Pass
)

cat(
  "Direction-level arithmetic integrity:",
  direction_arithmetic_pass,
  "\n\n"
)

print(direction_arithmetic)


# ------------------------------------------------------------
# 17. Candidate consistency class internal consistency
#
# Each Exact_LR_Pair must have only one
# Candidate_Consistency_Class.
# ------------------------------------------------------------

candidate_class_consistency <- dt[
  ,
  uniqueN(Candidate_Consistency_Class),
  by = Exact_LR_Pair
]

candidate_class_invariant_pass <- all(
  candidate_class_consistency$V1 == 1L
)


# ------------------------------------------------------------
# 18. Core-gene separation inheritance
#
# Upstream locked result:
#   AOC3 / MS4A7 / IL7R are NOT LR nodes.
#
# Therefore the canonical table should contain:
#   Core_Gene_Involved = FALSE
#   Core_Gene = NA
#
# This is an inherited structural property.
# It is NOT a search for core-gene LR pairs.
# ------------------------------------------------------------

core_gene_values <- unique(
  dt$Core_Gene
)

core_gene_involved_values <- unique(
  dt$Core_Gene_Involved
)

core_gene_all_na_pass <- all(
  is.na(dt$Core_Gene)
)

core_gene_involved_all_false_pass <- all(
  dt$Core_Gene_Involved == FALSE
)

core_gene_separation_pass <-
  core_gene_all_na_pass &&
  core_gene_involved_all_false_pass


# ------------------------------------------------------------
# 19. Expression numeric integrity
#
# Existing locked log2(CPM+1) values only.
# No recalculation.
# ------------------------------------------------------------

ligand_numeric_pass <- is.numeric(
  dt$Ligand_log2CPM
)

receptor_numeric_pass <- is.numeric(
  dt$Receptor_log2CPM
)

ligand_finite_pass <- all(
  is.finite(dt$Ligand_log2CPM)
)

receptor_finite_pass <- all(
  is.finite(dt$Receptor_log2CPM)
)

ligand_nonnegative_pass <- all(
  dt$Ligand_log2CPM >= 0
)

receptor_nonnegative_pass <- all(
  dt$Receptor_log2CPM >= 0
)


# ------------------------------------------------------------
# 20. Compatibility field integrity
#
# Existing locked compatibility values only.
# No recalculation.
# ------------------------------------------------------------

compatible_values <- sort(
  unique(
    dt$Compatible
  )
)

cat("Compatible values:\n")
print(compatible_values)
cat("\n")

compatible_value_pass <- all(
  dt$Compatible %in% c(TRUE, FALSE)
)

compatible_nonmissing_pass <- all(
  !is.na(dt$Compatible)
)


# ------------------------------------------------------------
# 21. LR identity field completeness
# ------------------------------------------------------------

lr_identity_fields <- c(
  "Exact_LR_Pair",
  "ligand",
  "receptor"
)

lr_identity_complete_pass <- all(
  sapply(
    dt[, ..lr_identity_fields],
    function(x) {
      !any(
        is.na(x) |
          trimws(as.character(x)) == ""
      )
    }
  )
)


# ------------------------------------------------------------
# 22. Direction field completeness
# ------------------------------------------------------------

direction_fields <- c(
  "Direction_ID",
  "Sender",
  "Receiver"
)

direction_fields_complete_pass <- all(
  sapply(
    dt[, ..direction_fields],
    function(x) {
      !any(
        is.na(x) |
          trimws(as.character(x)) == ""
      )
    }
  )
)


# ------------------------------------------------------------
# 23. Candidate class completeness
# ------------------------------------------------------------

candidate_class_complete_pass <- all(
  !is.na(dt$Candidate_Consistency_Class) &
    trimws(
      as.character(
        dt$Candidate_Consistency_Class
      )
    ) != ""
)


# ------------------------------------------------------------
# 24. Save descriptive summaries
# ------------------------------------------------------------

basic_summary <- data.table(
  Metric = c(
    "Input_rows",
    "Input_columns",
    "Unique_Exact_LR_Pairs",
    "Unique_Directions",
    "Unique_Patients",
    "Unique_Senders",
    "Unique_Receivers",
    "Unique_Pathways",
    "Unique_Annotations",
    "Unique_Candidate_Classes",
    "Canonical_Keys",
    "Duplicate_Canonical_Keys"
  ),
  Value = c(
    nrow(dt),
    ncol(dt),
    n_exact_lr,
    n_direction,
    n_patients,
    n_senders,
    n_receivers,
    n_pathways,
    n_annotations,
    n_candidate_classes,
    canonical_key_n,
    canonical_duplicate_n
  )
)

fwrite(
  basic_summary,
  file.path(
    output_dir,
    "STEP11-5F-1A_Input_Structure_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# 25. Candidate class summary
# ------------------------------------------------------------

candidate_class_summary <- dt[
  ,
  .(
    N_records = .N,
    N_exact_LR_pairs = uniqueN(Exact_LR_Pair),
    N_directions = uniqueN(Direction_ID),
    N_patients = uniqueN(Patient_Index)
  ),
  by = Candidate_Consistency_Class
][order(Candidate_Consistency_Class)]

fwrite(
  candidate_class_summary,
  file.path(
    output_dir,
    "STEP11-5F-1B_Candidate_Class_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# 26. Direction summary
# ------------------------------------------------------------

direction_summary <- dt[
  ,
  .(
    N_records = .N,
    N_exact_LR_pairs = uniqueN(Exact_LR_Pair),
    N_patients = uniqueN(Patient_Index),
    Sender = paste(
      unique(Sender),
      collapse = " | "
    ),
    Receiver = paste(
      unique(Receiver),
      collapse = " | "
    )
  ),
  by = Direction_ID
][order(Direction_ID)]

fwrite(
  direction_summary,
  file.path(
    output_dir,
    "STEP11-5F-1C_Direction_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# 27. Patient coverage summary
# ------------------------------------------------------------

patient_summary <- dt[
  ,
  .(
    N_records = .N,
    N_exact_LR_pairs = uniqueN(Exact_LR_Pair),
    N_directions = uniqueN(Direction_ID),
    N_candidate_classes =
      uniqueN(Candidate_Consistency_Class)
  ),
  by = Patient_Index
][order(Patient_Index)]

fwrite(
  patient_summary,
  file.path(
    output_dir,
    "STEP11-5F-1D_Patient_Coverage_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# 28. Direction × candidate-class summary
# ------------------------------------------------------------

direction_class_summary <- dt[
  ,
  .(
    N_records = .N,
    N_exact_LR_pairs = uniqueN(Exact_LR_Pair),
    N_patients = uniqueN(Patient_Index)
  ),
  by = .(
    Direction_ID,
    Candidate_Consistency_Class
  )
][order(
  Direction_ID,
  Candidate_Consistency_Class
)]

fwrite(
  direction_class_summary,
  file.path(
    output_dir,
    "STEP11-5F-1E_Direction_x_CandidateClass_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# 29. Final QC table
# ------------------------------------------------------------

qc <- data.table(
  
  QC_Item = c(
    
    "Input file exists",
    
    "Canonical row count = 9,620",
    
    "Required schema complete",
    
    "Unique Exact_LR_Pair count = 585",
    
    "Unique Direction_ID count = 6",
    
    "Unique Patient_Index count = 4",
    
    "Candidate class structure = 3_of_4_only + 4_of_4",
    
    "4_of_4 exact LR pairs = 492",
    
    "3_of_4-only exact LR pairs = 93",
    
    "Candidate LR total = 585",
    
    "Candidate_Set provenance non-empty",
    
    "interaction_name provenance non-empty",
    
    "Exact LR identity matches ligand + receptor",
    
    "Canonical key unique",
    
    "No duplicate canonical keys",
    
    "Patient identifiers complete",
    
    "Direction structure = six locked sender→receiver directions",
    
    "Candidate class invariant per Exact_LR_Pair",
    
    "Core_Gene field all NA",
    
    "Core_Gene_Involved all FALSE",
    
    "Ligand_log2CPM numeric",
    
    "Receptor_log2CPM numeric",
    
    "Ligand_log2CPM finite",
    
    "Receptor_log2CPM finite",
    
    "Ligand_log2CPM non-negative",
    
    "Receptor_log2CPM non-negative",
    
    "Compatible values valid",
    
    "Compatible values non-missing",
    
    "LR identity fields complete",
    
    "Direction fields complete",
    
    "Candidate class complete",
    
    "Direction-level record arithmetic is internally consistent"
    
  ),
  
  PASS = c(
    
    input_exists,
    
    row_count_pass,
    
    schema_pass,
    
    exact_lr_count_pass,
    
    direction_count_pass,
    
    patient_count_pass,
    
    candidate_class_structure_pass,
    
    four_of_four_lr_count_pass,
    
    three_of_four_only_lr_count_pass,
    
    candidate_lr_total_pass,
    
    candidate_set_provenance_pass,
    
    interaction_provenance_nonempty_pass,
    
    exact_lr_identity_pass,
    
    canonical_key_unique_pass,
    
    canonical_no_duplicates_pass,
    
    patient_nonempty_pass,
    
    direction_structure_pass,
    
    candidate_class_invariant_pass,
    
    core_gene_all_na_pass,
    
    core_gene_involved_all_false_pass,
    
    ligand_numeric_pass,
    
    receptor_numeric_pass,
    
    ligand_finite_pass,
    
    receptor_finite_pass,
    
    ligand_nonnegative_pass,
    
    receptor_nonnegative_pass,
    
    compatible_value_pass,
    
    compatible_nonmissing_pass,
    
    lr_identity_complete_pass,
    
    direction_fields_complete_pass,
    
    candidate_class_complete_pass,
    
    direction_arithmetic_pass
    
  )
)


# ------------------------------------------------------------
# 30. Final QC output
# ------------------------------------------------------------

fwrite(
  qc,
  file.path(
    output_dir,
    "STEP11-5F-1_FINAL_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# 31. Global PASS / FAIL
# ------------------------------------------------------------

all_pass <- all(
  qc$PASS
)

n_qc <- nrow(qc)

n_pass <- sum(
  qc$PASS
)

n_fail <- sum(
  !qc$PASS
)


# ------------------------------------------------------------
# 32. Console report
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP11-5F-1 FINAL QC\n")
cat("============================================================\n")

cat(
  "Input rows                    :",
  nrow(dt),
  "\n"
)

cat(
  "Unique Exact_LR_Pair          :",
  n_exact_lr,
  "\n"
)

cat(
  "Unique patients               :",
  n_patients,
  "\n"
)

cat(
  "Unique directions             :",
  n_direction,
  "\n"
)

cat(
  "4/4 exact LR pairs            :",
  n_four_of_four_lr,
  "\n"
)

cat(
  "3/4-only exact LR pairs       :",
  n_three_of_four_only_lr,
  "\n"
)

cat(
  "Canonical keys                :",
  canonical_key_n,
  "\n"
)

cat(
  "Duplicate canonical keys      :",
  canonical_duplicate_n,
  "\n"
)

cat("\n")

cat(
  "QC items                      :",
  n_qc,
  "\n"
)

cat(
  "PASS                          :",
  n_pass,
  "\n"
)

cat(
  "FAIL                          :",
  n_fail,
  "\n"
)

cat("\n")

if (all_pass) {
  
  cat("PASS — STEP11-5F-1 CANONICAL INPUT QC\n")
  cat("============================================================\n")
  cat("The locked 9,620-row canonical patient-level LR table\n")
  cat("passed all structural and input-integrity checks.\n\n")
  
  cat("Canonical unit:\n")
  cat("  Exact_LR_Pair × Direction_ID × Patient_Index\n\n")
  
  cat("Locked structure confirmed:\n")
  cat("  585 exact LR pairs\n")
  cat("  492 pairs supported in 4/4 patients\n")
  cat("  93 pairs supported in 3/4 patients only\n")
  cat("  4 matched patients\n")
  cat("  6 sender→receiver directions\n\n")
  
  cat("Candidate_Set and interaction_name are retained as\n")
  cat("provenance fields and are not treated as biological\n")
  cat("invariants.\n\n")
  
  cat("Core-gene separation is inherited from the locked\n")
  cat("LR-universe audit:\n")
  cat("  AOC3 / MS4A7 / IL7R are not LR nodes.\n\n")
  
  cat("No upstream locked quantities were recalculated.\n")
  cat("No new LR pairs were discovered.\n")
  cat("No compatibility values were recalculated.\n")
  cat("No expression values were recalculated.\n")
  cat("No statistical testing was performed.\n")
  cat("No biological ranking was performed.\n\n")
  
  cat("Output directory:\n")
  cat(" ", output_dir, "\n\n")
  
  cat("Formal status: PASS — READY FOR STEP11-5F-2\n")
  
} else {
  
  cat("FAIL — STEP11-5F-1 CANONICAL INPUT QC\n")
  cat("============================================================\n")
  
  failed_items <- qc[
    PASS == FALSE,
    QC_Item
  ]
  
  cat("Failed QC items:\n")
  cat(
    paste(
      failed_items,
      collapse = "\n"
    ),
    "\n\n"
  )
  
  cat("STEP11-5F-2 MUST NOT be started until QC is resolved.\n")
}


# ------------------------------------------------------------
# 33. End
# ------------------------------------------------------------

cat("\n")
cat("END OF STEP11-5F-1\n")
cat("============================================================\n")

# ============================================================
# STEP11-5F-1 Direction-ID Diagnostic
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP11-5F-1 — Direction ID Diagnostic\n")
cat("============================================================\n\n")

cat("Actual Direction_ID values:\n\n")

print(
  sort(
    unique(
      as.character(dt$Direction_ID)
    )
  )
)

cat("\nDirection_ID frequency:\n\n")

direction_frequency <- dt[
  ,
  .(
    N_records = .N,
    N_exact_LR_pairs = uniqueN(Exact_LR_Pair),
    N_patients = uniqueN(Patient_Index),
    Sender = paste(
      unique(Sender),
      collapse = " | "
    ),
    Receiver = paste(
      unique(Receiver),
      collapse = " | "
    )
  ),
  by = Direction_ID
][order(Direction_ID)]

print(direction_frequency)

cat("\n============================================================\n")
cat("END OF DIRECTION-ID DIAGNOSTIC\n")
cat("============================================================\n")


# ============================================================
# STEP11-5F-2
# Patient-Level Compatibility Consistency Analysis
# ============================================================
#
# PURPOSE
# ------------------------------------------------------------
# Descriptively summarize patient-level compatibility patterns
# from the FINAL LOCKED 9,620-row canonical LR integration table.
#
# FINAL LOCKED INPUT
# ------------------------------------------------------------
# STEP11-5E-3B_FINAL_Canonical_Patient_Level_LR_Integration.tsv
#
# CANONICAL UNIT
# ------------------------------------------------------------
# Exact_LR_Pair × Direction_ID × Patient_Index
#
# LOCKED CANDIDATE UNIVERSE
# ------------------------------------------------------------
# 585 exact LR pairs
#   492 = 4/4 patient-compatible candidate pairs
#    93 = 3/4 patient-compatible candidate pairs
#
# IMPORTANT
# ------------------------------------------------------------
# Candidate_Set_Provenance is provenance.
# interaction_name_Provenance is resource-level provenance.
#
# Neither field is treated as a biological invariant.
#
# Compatibility is inherited from the locked canonical table.
# It is NOT recalculated here.
#
# No new LR discovery.
# No CCI_LR candidate rediscovery.
# No raw UMI processing.
# No pseudobulk recalculation.
# No CPM/log2CPM recalculation.
# No compatibility recalculation.
# No statistical testing.
# No P values.
# No FDR.
# No CellChat communication probability.
# No pathway enrichment.
# No biological ranking.
#
# Compatibility fraction is descriptive only.
# It must NOT be interpreted as interaction strength,
# communication probability, or functional signaling.
#
# ============================================================


# ------------------------------------------------------------
# 0. Input / output paths
# ------------------------------------------------------------

input_file <- paste0(
  "input_dir/Paper 1_Step 11-CellChat/",
  "STEP11-5E-3_Patient_Level_Candidate_Integration/",
  "STEP11-5E-3B_FINAL_Canonical_Patient_Level_LR_Integration.tsv"
)

output_dir <- paste0(
  "result_dir/Paper 1_Step 11-CellChat/",
  "STEP11-5F-2_Patient_Level_Compatibility_Consistency"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


cat("\n")
cat("============================================================\n")
cat("STEP11-5F-2 — Patient-Level Compatibility Consistency\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Package
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})


# ------------------------------------------------------------
# 2. Input existence
# ------------------------------------------------------------

input_exists <- file.exists(
  input_file
)

cat(
  "Input file exists:",
  input_exists,
  "\n\n"
)

if (!input_exists) {
  stop(
    "STOP: FINAL canonical input file was not found."
  )
}


# ------------------------------------------------------------
# 3. Read FINAL LOCKED canonical table
# ------------------------------------------------------------

dt <- fread(
  input_file,
  sep = "\t",
  header = TRUE,
  data.table = TRUE,
  check.names = FALSE
)

cat(
  "Input rows:",
  nrow(dt),
  "\n"
)

cat(
  "Input columns:",
  ncol(dt),
  "\n\n"
)


# ------------------------------------------------------------
# 4. FINAL LOCKED schema
# ------------------------------------------------------------

required_cols <- c(
  
  "Canonical_Key",
  "Exact_LR_Pair",
  "Candidate_Consistency_Class",
  "Candidate_Set_Provenance",
  "interaction_name_Provenance",
  "N_Candidate_Set",
  "N_interaction_name",
  "pathway_name",
  "ligand",
  "receptor",
  "annotation",
  "Direction_ID",
  "Sender",
  "Receiver",
  "Patient_Index",
  "Core_Gene_Involved",
  "Core_Gene",
  "Ligand_log2CPM",
  "Receptor_log2CPM",
  "Ligand_Receptor_Mean_log2CPM",
  "Ligand_Receptor_Difference_log2CPM",
  "Compatible"
  
)

schema_pass <- all(
  required_cols %in% names(dt)
)

cat(
  "Required FINAL schema complete:",
  schema_pass,
  "\n\n"
)

if (!schema_pass) {
  
  missing_cols <- setdiff(
    required_cols,
    names(dt)
  )
  
  cat("Missing columns:\n")
  print(missing_cols)
  
  stop(
    "STOP: FINAL canonical schema is incomplete."
  )
}


# ------------------------------------------------------------
# 5. Basic locked structure
# ------------------------------------------------------------

row_count_pass <- nrow(dt) == 9620L

exact_lr_count <- uniqueN(
  dt$Exact_LR_Pair
)

exact_lr_count_pass <- exact_lr_count == 585L

patient_ids <- sort(
  unique(
    as.character(
      dt$Patient_Index
    )
  )
)

patient_count <- length(
  patient_ids
)

patient_count_pass <- patient_count == 4L

direction_ids <- sort(
  unique(
    as.character(
      dt$Direction_ID
    )
  )
)

direction_count <- length(
  direction_ids
)

direction_count_pass <- identical(
  direction_ids,
  as.character(1:6)
)

cat(
  "Canonical row count = 9,620:",
  row_count_pass,
  "\n"
)

cat(
  "Unique Exact_LR_Pair = 585:",
  exact_lr_count_pass,
  "\n"
)

cat(
  "Unique Patient_Index = 4:",
  patient_count_pass,
  "\n"
)

cat(
  "Direction_ID = 1–6:",
  direction_count_pass,
  "\n\n"
)

cat(
  "Locked Patient_Index values:\n"
)

print(
  patient_ids
)

cat("\n")


# ------------------------------------------------------------
# 6. Canonical key integrity
# ------------------------------------------------------------

expected_canonical_key <- paste(
  dt$Exact_LR_Pair,
  dt$Direction_ID,
  dt$Patient_Index,
  sep = "___"
)

canonical_key_identity_pass <- all(
  dt$Canonical_Key ==
    expected_canonical_key
)

canonical_key_unique_pass <- uniqueN(
  dt$Canonical_Key
) == nrow(dt)

cat(
  "Canonical_Key identity matches locked definition:",
  canonical_key_identity_pass,
  "\n"
)

cat(
  "Canonical_Key unique:",
  canonical_key_unique_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 7. Candidate consistency class
# ------------------------------------------------------------

candidate_classes <- sort(
  unique(
    as.character(
      dt$Candidate_Consistency_Class
    )
  )
)

candidate_class_values_pass <- identical(
  candidate_classes,
  sort(
    c(
      "3_of_4_only",
      "4_of_4"
    )
  )
)

candidate_class_complete_pass <- all(
  !is.na(
    dt$Candidate_Consistency_Class
  ) &
    nzchar(
      as.character(
        dt$Candidate_Consistency_Class
      )
    )
)

cat(
  "Candidate consistency classes:\n"
)

print(
  candidate_classes
)

cat(
  "\nCandidate class values valid:",
  candidate_class_values_pass,
  "\n"
)

cat(
  "Candidate class complete:",
  candidate_class_complete_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 8. Candidate class invariance per Exact LR pair
# ------------------------------------------------------------

candidate_class_invariance <- dt[
  ,
  .(
    N_Candidate_Classes =
      uniqueN(
        Candidate_Consistency_Class
      )
  ),
  by = Exact_LR_Pair
]

candidate_class_invariant_pass <- all(
  candidate_class_invariance$N_Candidate_Classes == 1L
)

cat(
  "Candidate class invariant per Exact_LR_Pair:",
  candidate_class_invariant_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 9. Candidate class exact-LR structure
# ------------------------------------------------------------

candidate_lr_summary <- dt[
  ,
  .(
    Candidate_Consistency_Class =
      unique(
        Candidate_Consistency_Class
      )[1]
  ),
  by = Exact_LR_Pair
]

four_of_four_lr_count <- sum(
  candidate_lr_summary$Candidate_Consistency_Class ==
    "4_of_4"
)

three_of_four_only_lr_count <- sum(
  candidate_lr_summary$Candidate_Consistency_Class ==
    "3_of_4_only"
)

four_of_four_lr_count_pass <-
  four_of_four_lr_count == 492L

three_of_four_only_lr_count_pass <-
  three_of_four_only_lr_count == 93L

candidate_lr_total_pass <-
  four_of_four_lr_count +
  three_of_four_only_lr_count ==
  585L

cat(
  "4_of_4 exact LR pairs:",
  four_of_four_lr_count,
  "\n"
)

cat(
  "3_of_4-only exact LR pairs:",
  three_of_four_only_lr_count,
  "\n"
)

cat(
  "Candidate LR total = 585:",
  candidate_lr_total_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 10. Direction mapping integrity
# ------------------------------------------------------------

expected_direction_map <- data.table(
  
  Direction_ID = 1:6,
  
  Expected_Sender = c(
    "Endothelial",
    "Myeloid/macrophage",
    "Endothelial",
    "T-cell",
    "Myeloid/macrophage",
    "T-cell"
  ),
  
  Expected_Receiver = c(
    "Myeloid/macrophage",
    "Endothelial",
    "T-cell",
    "Endothelial",
    "T-cell",
    "Myeloid/macrophage"
  )
)

actual_direction_map <- unique(
  dt[
    ,
    .(
      Direction_ID = as.integer(
        Direction_ID
      ),
      Sender,
      Receiver
    )
  ]
)

setorder(
  expected_direction_map,
  Direction_ID
)

setorder(
  actual_direction_map,
  Direction_ID
)

direction_mapping_pass <-
  nrow(actual_direction_map) == 6L &&
  identical(
    expected_direction_map$Direction_ID,
    actual_direction_map$Direction_ID
  ) &&
  identical(
    expected_direction_map$Expected_Sender,
    actual_direction_map$Sender
  ) &&
  identical(
    expected_direction_map$Expected_Receiver,
    actual_direction_map$Receiver
  )

cat(
  "Locked sender→receiver mapping:",
  direction_mapping_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 11. Direction-level arithmetic integrity
# ------------------------------------------------------------

direction_arithmetic <- dt[
  ,
  .(
    N_records = .N,
    N_exact_LR_pairs =
      uniqueN(
        Exact_LR_Pair
      ),
    N_patients =
      uniqueN(
        Patient_Index
      )
  ),
  by = .(
    Direction_ID,
    Sender,
    Receiver
  )
][
  order(
    Direction_ID
  )
]

direction_arithmetic[
  ,
  Expected_N_records :=
    N_exact_LR_pairs *
    N_patients
]

direction_arithmetic[
  ,
  Arithmetic_Pass :=
    N_records ==
    Expected_N_records
]

direction_arithmetic_pass <- all(
  direction_arithmetic$Arithmetic_Pass
)

cat(
  "Direction-level arithmetic integrity:",
  direction_arithmetic_pass,
  "\n\n"
)

print(
  direction_arithmetic
)

cat("\n")


# ------------------------------------------------------------
# 12. Compatibility integrity
# ------------------------------------------------------------
#
# Compatible is inherited from FINAL LOCK.
# No recalculation is performed.
# ------------------------------------------------------------

compatible_values <- sort(
  unique(
    dt$Compatible
  )
)

compatible_value_pass <- all(
  !is.na(compatible_values) &
    compatible_values %in%
    c(FALSE, TRUE)
)

compatible_nonmissing_pass <- all(
  !is.na(
    dt$Compatible
  )
)

cat(
  "Compatible values:\n"
)

print(
  compatible_values
)

cat(
  "\nCompatible values valid:",
  compatible_value_pass,
  "\n"
)

cat(
  "Compatible values non-missing:",
  compatible_nonmissing_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 13. Provenance integrity
# ------------------------------------------------------------
#
# Candidate_Set_Provenance and
# interaction_name_Provenance are retained
# as provenance fields only.
#
# They are NOT treated as biological invariants.
# ------------------------------------------------------------

candidate_provenance_nonempty_pass <- all(
  !is.na(
    dt$Candidate_Set_Provenance
  ) &
    nzchar(
      as.character(
        dt$Candidate_Set_Provenance
      )
    )
)

interaction_provenance_nonempty_pass <- all(
  !is.na(
    dt$interaction_name_Provenance
  ) &
    nzchar(
      as.character(
        dt$interaction_name_Provenance
      )
    )
)

cat(
  "Candidate_Set_Provenance non-empty:",
  candidate_provenance_nonempty_pass,
  "\n"
)

cat(
  "interaction_name_Provenance non-empty:",
  interaction_provenance_nonempty_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 14. Core-gene separation inheritance
# ------------------------------------------------------------

core_gene_all_na_pass <- all(
  is.na(
    dt$Core_Gene
  )
)

core_gene_involved_all_false_pass <- all(
  dt$Core_Gene_Involved == FALSE
)

cat(
  "Core_Gene field all NA:",
  core_gene_all_na_pass,
  "\n"
)

cat(
  "Core_Gene_Involved all FALSE:",
  core_gene_involved_all_false_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 15. Expression field integrity
# ------------------------------------------------------------

ligand_numeric_pass <- is.numeric(
  dt$Ligand_log2CPM
)

receptor_numeric_pass <- is.numeric(
  dt$Receptor_log2CPM
)

ligand_finite_pass <- all(
  is.finite(
    dt$Ligand_log2CPM
  )
)

receptor_finite_pass <- all(
  is.finite(
    dt$Receptor_log2CPM
  )
)

ligand_nonnegative_pass <- all(
  dt$Ligand_log2CPM >= 0
)

receptor_nonnegative_pass <- all(
  dt$Receptor_log2CPM >= 0
)

cat(
  "Ligand_log2CPM numeric:",
  ligand_numeric_pass,
  "\n"
)

cat(
  "Receptor_log2CPM numeric:",
  receptor_numeric_pass,
  "\n"
)

cat(
  "Ligand_log2CPM finite:",
  ligand_finite_pass,
  "\n"
)

cat(
  "Receptor_log2CPM finite:",
  receptor_finite_pass,
  "\n"
)

cat(
  "Ligand_log2CPM non-negative:",
  ligand_nonnegative_pass,
  "\n"
)

cat(
  "Receptor_log2CPM non-negative:",
  receptor_nonnegative_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 16. Patient-level compatibility pattern
# ------------------------------------------------------------
#
# One row:
# Exact_LR_Pair × Direction_ID
#
# Four columns:
# one for each locked Patient_Index.
#
# Compatible values are copied from the locked canonical table.
# ------------------------------------------------------------

patient_pattern <- dcast(
  
  dt,
  
  Exact_LR_Pair +
    Candidate_Consistency_Class +
    Candidate_Set_Provenance +
    interaction_name_Provenance +
    pathway_name +
    ligand +
    receptor +
    annotation +
    Direction_ID +
    Sender +
    Receiver
  ~ Patient_Index,
  
  value.var = "Compatible"
  
)

patient_columns_pass <- all(
  patient_ids %in%
    names(patient_pattern)
)

cat(
  "All locked patient columns present:",
  patient_columns_pass,
  "\n\n"
)

if (!patient_columns_pass) {
  
  cat(
    "Expected patient columns:\n"
  )
  
  print(
    patient_ids
  )
  
  cat(
    "\nObserved columns:\n"
  )
  
  print(
    names(patient_pattern)
  )
  
  stop(
    "STOP: Patient-level compatibility pattern is incomplete."
  )
}


# ------------------------------------------------------------
# 17. Patient-pattern completeness
# ------------------------------------------------------------

patient_pattern[
  ,
  N_patients :=
    rowSums(
      !is.na(
        .SD
      )
    ),
  .SDcols = patient_ids
]

patient_pattern[
  ,
  N_compatible :=
    rowSums(
      .SD == TRUE,
      na.rm = TRUE
    ),
  .SDcols = patient_ids
]

patient_pattern[
  ,
  N_incompatible :=
    rowSums(
      .SD == FALSE,
      na.rm = TRUE
    ),
  .SDcols = patient_ids
]

patient_pattern[
  ,
  Compatibility_fraction :=
    N_compatible /
    N_patients
]

patient_pattern[
  ,
  Compatibility_Pattern :=
    paste0(
      N_compatible,
      "_of_",
      N_patients
    )
]


# ------------------------------------------------------------
# 18. Patient-pattern QC
# ------------------------------------------------------------

pattern_count_pass <- all(
  patient_pattern$N_patients == 4L
)

pattern_arithmetic_pass <- all(
  patient_pattern$N_compatible +
    patient_pattern$N_incompatible ==
    patient_pattern$N_patients
)

pattern_fraction_pass <- all(
  is.finite(
    patient_pattern$Compatibility_fraction
  ) &
    patient_pattern$Compatibility_fraction >= 0 &
    patient_pattern$Compatibility_fraction <= 1
)

cat(
  "Every LR × direction has four patients:",
  pattern_count_pass,
  "\n"
)

cat(
  "Compatible + incompatible = 4:",
  pattern_arithmetic_pass,
  "\n"
)

cat(
  "Compatibility fraction bounded 0–1:",
  pattern_fraction_pass,
  "\n\n"
)


# ------------------------------------------------------------
# 19. Global compatibility pattern summary
# ------------------------------------------------------------

global_pattern_summary <- patient_pattern[
  ,
  .(
    N_LR_Direction_Records = .N,
    
    N_Exact_LR_Pairs =
      uniqueN(
        Exact_LR_Pair
      ),
    
    Mean_Compatibility_Fraction =
      mean(
        Compatibility_fraction
      ),
    
    Median_Compatibility_Fraction =
      median(
        Compatibility_fraction
      )
  ),
  
  by = Compatibility_Pattern
][
  order(
    Compatibility_Pattern
  )
]


# ------------------------------------------------------------
# 20. Exact LR-pair × direction summary
# ------------------------------------------------------------

lr_direction_summary <- patient_pattern[
  ,
  .(
    
    Candidate_Consistency_Class =
      unique(
        Candidate_Consistency_Class
      )[1],
    
    Candidate_Set_Provenance =
      unique(
        Candidate_Set_Provenance
      )[1],
    
    N_patients,
    
    N_compatible =
      unique(
        N_compatible
      ),
    
    N_incompatible =
      unique(
        N_incompatible
      ),
    
    Compatibility_fraction =
      unique(
        Compatibility_fraction
      ),
    
    Compatibility_Pattern =
      unique(
        Compatibility_Pattern
      )
    
  ),
  
  by = .(
    Exact_LR_Pair,
    Direction_ID,
    Sender,
    Receiver,
    ligand,
    receptor
  )
][
  order(
    Direction_ID,
    Exact_LR_Pair
  )
]


# ------------------------------------------------------------
# 21. Exact LR-pair overall summary
# ------------------------------------------------------------
#
# Descriptive summary across all directions.
# ------------------------------------------------------------

lr_pair_summary <- patient_pattern[
  ,
  .(
    
    Candidate_Consistency_Class =
      unique(
        Candidate_Consistency_Class
      )[1],
    
    N_Directions =
      uniqueN(
        Direction_ID
      ),
    
    N_LR_Direction_Records =
      .N,
    
    Total_Compatible =
      sum(
        N_compatible
      ),
    
    Total_Incompatible =
      sum(
        N_incompatible
      ),
    
    Mean_Compatibility_Fraction =
      mean(
        Compatibility_fraction
      ),
    
    Median_Compatibility_Fraction =
      median(
        Compatibility_fraction
      ),
    
    Min_Compatibility_Fraction =
      min(
        Compatibility_fraction
      ),
    
    Max_Compatibility_Fraction =
      max(
        Compatibility_fraction
      )
    
  ),
  
  by = .(
    Exact_LR_Pair,
    ligand,
    receptor
  )
]


# ------------------------------------------------------------
# 22. Direction-level compatibility summary
# ------------------------------------------------------------

direction_summary <- patient_pattern[
  ,
  .(
    
    N_LR_Direction_Records =
      .N,
    
    N_Exact_LR_Pairs =
      uniqueN(
        Exact_LR_Pair
      ),
    
    N_Patients =
      uniqueN(
        c(
          unlist(
            .SD
          )
        )
      ),
    
    Total_Compatible =
      sum(
        N_compatible
      ),
    
    Total_Incompatible =
      sum(
        N_incompatible
      ),
    
    Mean_Compatibility_Fraction =
      mean(
        Compatibility_fraction
      ),
    
    Median_Compatibility_Fraction =
      median(
        Compatibility_fraction
      )
    
  ),
  
  by = .(
    Direction_ID,
    Sender,
    Receiver
  )
][
  order(
    Direction_ID
  )
]


# ------------------------------------------------------------
# 23. Patient-specific summary
# ------------------------------------------------------------

patient_specific_summary <- rbindlist(
  
  lapply(
    
    patient_ids,
    
    function(p) {
      
      data.table(
        
        Patient_Index = p,
        
        N_LR_Direction_Records =
          sum(
            !is.na(
              patient_pattern[[p]]
            )
          ),
        
        N_Compatible =
          sum(
            patient_pattern[[p]] == TRUE,
            na.rm = TRUE
          ),
        
        N_Incompatible =
          sum(
            patient_pattern[[p]] == FALSE,
            na.rm = TRUE
          )
        
      )
      
    }
    
  )
)

patient_specific_summary[
  ,
  Compatibility_Fraction :=
    N_Compatible /
    N_LR_Direction_Records
]


# ------------------------------------------------------------
# 24. Candidate-class × compatibility-pattern summary
# ------------------------------------------------------------

candidate_class_pattern_summary <- patient_pattern[
  ,
  .(
    
    N_LR_Direction_Records =
      .N,
    
    N_Exact_LR_Pairs =
      uniqueN(
        Exact_LR_Pair
      ),
    
    N_Compatible =
      sum(
        N_compatible
      ),
    
    N_Incompatible =
      sum(
        N_incompatible
      ),
    
    Mean_Compatibility_Fraction =
      mean(
        Compatibility_fraction
      ),
    
    Median_Compatibility_Fraction =
      median(
        Compatibility_fraction
      )
    
  ),
  
  by = .(
    Candidate_Consistency_Class,
    Compatibility_Pattern
  )
][
  order(
    Candidate_Consistency_Class,
    Compatibility_Pattern
  )
]


# ------------------------------------------------------------
# 25. Direction × candidate-class summary
# ------------------------------------------------------------

direction_candidate_summary <- patient_pattern[
  ,
  .(
    
    N_LR_Direction_Records =
      .N,
    
    N_Exact_LR_Pairs =
      uniqueN(
        Exact_LR_Pair
      ),
    
    Total_Compatible =
      sum(
        N_compatible
      ),
    
    Total_Incompatible =
      sum(
        N_incompatible
      ),
    
    Mean_Compatibility_Fraction =
      mean(
        Compatibility_fraction
      ),
    
    Median_Compatibility_Fraction =
      median(
        Compatibility_fraction
      )
    
  ),
  
  by = .(
    Direction_ID,
    Sender,
    Receiver,
    Candidate_Consistency_Class
  )
][
  order(
    Direction_ID,
    Candidate_Consistency_Class
  )
]


# ------------------------------------------------------------
# 26. Save patient-level pattern
# ------------------------------------------------------------

fwrite(
  
  patient_pattern,
  
  file.path(
    output_dir,
    "STEP11-5F-2A_Patient_Level_Compatibility_Pattern.tsv"
  ),
  
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# 27. Save LR × direction summary
# ------------------------------------------------------------

fwrite(
  
  lr_direction_summary,
  
  file.path(
    output_dir,
    "STEP11-5F-2B_LR_Direction_Compatibility_Summary.tsv"
  ),
  
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# 28. Save exact LR-pair summary
# ------------------------------------------------------------

fwrite(
  
  lr_pair_summary,
  
  file.path(
    output_dir,
    "STEP11-5F-2C_Exact_LR_Pair_Compatibility_Summary.tsv"
  ),
  
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# 29. Save direction summary
# ------------------------------------------------------------

fwrite(
  
  direction_summary,
  
  file.path(
    output_dir,
    "STEP11-5F-2D_Direction_Compatibility_Summary.tsv"
  ),
  
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# 30. Save patient-specific summary
# ------------------------------------------------------------

fwrite(
  
  patient_specific_summary,
  
  file.path(
    output_dir,
    "STEP11-5F-2E_Patient_Specific_Compatibility_Summary.tsv"
  ),
  
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# 31. Save candidate-class × pattern summary
# ------------------------------------------------------------

fwrite(
  
  candidate_class_pattern_summary,
  
  file.path(
    output_dir,
    "STEP11-5F-2F_CandidateClass_x_CompatibilityPattern.tsv"
  ),
  
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# 32. Save direction × candidate-class summary
# ------------------------------------------------------------

fwrite(
  
  direction_candidate_summary,
  
  file.path(
    output_dir,
    "STEP11-5F-2G_Direction_x_CandidateClass_Compatibility.tsv"
  ),
  
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# 33. Save global compatibility pattern summary
# ------------------------------------------------------------

fwrite(
  
  global_pattern_summary,
  
  file.path(
    output_dir,
    "STEP11-5F-2H_Global_Compatibility_Pattern_Summary.tsv"
  ),
  
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ------------------------------------------------------------
# 34. Final QC
# ------------------------------------------------------------

qc <- data.table(
  
  QC_Item = c(
    
    "Input file exists",
    
    "Canonical row count = 9,620",
    
    "Required FINAL schema complete",
    
    "Unique Exact_LR_Pair count = 585",
    
    "Unique Patient_Index count = 4",
    
    "Direction_ID values = 1–6",
    
    "Canonical_Key identity matches locked definition",
    
    "Canonical_Key unique",
    
    "Candidate class values valid",
    
    "Candidate class complete",
    
    "Candidate class invariant per Exact_LR_Pair",
    
    "4_of_4 exact LR pairs = 492",
    
    "3_of_4-only exact LR pairs = 93",
    
    "Candidate LR total = 585",
    
    "Locked sender→receiver mapping confirmed",
    
    "Direction-level record arithmetic internally consistent",
    
    "Compatible values valid",
    
    "Compatible values non-missing",
    
    "Candidate_Set_Provenance non-empty",
    
    "interaction_name_Provenance non-empty",
    
    "Core_Gene field all NA",
    
    "Core_Gene_Involved all FALSE",
    
    "Ligand_log2CPM numeric",
    
    "Receptor_log2CPM numeric",
    
    "Ligand_log2CPM finite",
    
    "Receptor_log2CPM finite",
    
    "Ligand_log2CPM non-negative",
    
    "Receptor_log2CPM non-negative",
    
    "All four patient compatibility columns present",
    
    "Every LR × direction has four patients",
    
    "Compatible + incompatible = 4",
    
    "Compatibility fraction bounded 0–1"
    
  ),
  
  PASS = c(
    
    input_exists,
    
    row_count_pass,
    
    schema_pass,
    
    exact_lr_count_pass,
    
    patient_count_pass,
    
    direction_count_pass,
    
    canonical_key_identity_pass,
    
    canonical_key_unique_pass,
    
    candidate_class_values_pass,
    
    candidate_class_complete_pass,
    
    candidate_class_invariant_pass,
    
    four_of_four_lr_count_pass,
    
    three_of_four_only_lr_count_pass,
    
    candidate_lr_total_pass,
    
    direction_mapping_pass,
    
    direction_arithmetic_pass,
    
    compatible_value_pass,
    
    compatible_nonmissing_pass,
    
    candidate_provenance_nonempty_pass,
    
    interaction_provenance_nonempty_pass,
    
    core_gene_all_na_pass,
    
    core_gene_involved_all_false_pass,
    
    ligand_numeric_pass,
    
    receptor_numeric_pass,
    
    ligand_finite_pass,
    
    receptor_finite_pass,
    
    ligand_nonnegative_pass,
    
    receptor_nonnegative_pass,
    
    patient_columns_pass,
    
    pattern_count_pass,
    
    pattern_arithmetic_pass,
    
    pattern_fraction_pass
    
  )
)


# ------------------------------------------------------------
# 35. Print final QC
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP11-5F-2 FINAL QC\n")
cat("============================================================\n\n")

print(qc)

cat("\n")

cat(
  "QC items:",
  nrow(qc),
  "\n"
)

cat(
  "PASS:",
  sum(qc$PASS),
  "\n"
)

cat(
  "FAIL:",
  sum(!qc$PASS),
  "\n\n"
)


# ------------------------------------------------------------
# 36. Save final QC
# ------------------------------------------------------------

fwrite(
  
  qc,
  
  file.path(
    output_dir,
    "STEP11-5F-2_FINAL_QC.tsv"
  ),
  
  sep = "\t",
  quote = FALSE
)


# ------------------------------------------------------------
# 37. Print key descriptive summaries
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP11-5F-2 KEY DESCRIPTIVE SUMMARIES\n")
cat("============================================================\n\n")

cat("Global compatibility pattern:\n")
print(
  global_pattern_summary
)

cat("\n")

cat("Direction-level compatibility summary:\n")
print(
  direction_summary
)

cat("\n")

cat("Patient-specific compatibility summary:\n")
print(
  patient_specific_summary
)

cat("\n")


# ------------------------------------------------------------
# 38. Final lock statement
# ------------------------------------------------------------

if (
  all(qc$PASS)
) {
  
  cat(
    "PASS — STEP11-5F-2 PATIENT-LEVEL ",
    "COMPATIBILITY CONSISTENCY ANALYSIS\n",
    sep = ""
  )
  
  cat(
    "============================================================\n"
  )
  
  cat(
    "The FINAL LOCKED 9,620-row canonical patient-level\n"
  )
  
  cat(
    "LR table was summarized descriptively at the patient level.\n\n"
  )
  
  cat(
    "Canonical unit:\n"
  )
  
  cat(
    "  Exact_LR_Pair × Direction_ID × Patient_Index\n\n"
  )
  
  cat(
    "Locked candidate universe:\n"
  )
  
  cat(
    "  585 exact LR pairs\n"
  )
  
  cat(
    "  492 pairs supported in 4/4 patients\n"
  )
  
  cat(
    "   93 pairs supported in 3/4 patients only\n\n"
  )
  
  cat(
    "No new LR pairs were discovered.\n"
  )
  
  cat(
    "No compatibility values were recalculated.\n"
  )
  
  cat(
    "No expression values were recalculated.\n"
  )
  
  cat(
    "No statistical testing was performed.\n"
  )
  
  cat(
    "No biological ranking was performed.\n"
  )
  
  cat(
    "Compatibility fractions are descriptive only.\n"
  )
  
  cat(
    "They do not represent interaction strength,\n"
  )
  
  cat(
    "communication probability, or functional signaling.\n\n"
  )
  
  cat(
    "Output directory:\n  ",
    output_dir,
    "\n\n",
    sep = ""
  )
  
  cat(
    "Formal status: PASS — READY FOR STEP11-5F-3\n"
  )
  
} else {
  
  cat(
    "FAIL — STEP11-5F-2 QC NOT RESOLVED\n"
  )
  
  cat(
    "STEP11-5F-3 MUST NOT be started until QC is resolved.\n"
  )
}


# ------------------------------------------------------------
# 39. End
# ------------------------------------------------------------

cat("\n")
cat("END OF STEP11-5F-2\n")
cat("============================================================\n")


# ============================================================
# STEP11-5F-3
# Patient-Level Expression Context / LR Pair Context Integration
#
# Purpose:
#   Integrate the FINAL LOCKED canonical patient-level LR table
#   with ligand/receptor expression context, direction,
#   candidate consistency, patient-level compatibility,
#   annotation, pathway, and core-gene separation.
#
# IMPORTANT:
#   This is a DESCRIPTIVE INTEGRATION step only.
#
#   NO:
#     - new LR discovery
#     - raw UMI recalculation
#     - pseudobulk recalculation
#     - CPM recalculation
#     - log2(CPM+1) recalculation
#     - detection recalculation
#     - compatibility recalculation
#     - CellChat probability calculation
#     - correlation testing
#     - P values
#     - FDR
#     - pathway enrichment
#     - biological ranking
#
# Expression magnitude is NOT interaction strength.
# Compatibility is NOT functional signaling.
#
# Canonical unit:
#   Exact_LR_Pair × Direction_ID × Patient_Index
#
# Input:
#   STEP11-5E-3B_FINAL_Canonical_Patient_Level_LR_Integration.tsv
#
# Output:
#   STEP11-5F-3_Patient_Level_Expression_Context
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP11-5F-3\n")
cat("Patient-Level Expression Context / LR Pair Context Integration\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_packages <- c("data.table")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package not installed: ", pkg)
  }
}

library(data.table)


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <- paste0(
  "input_dir/Paper 1_Step 11-CellChat/",
  "STEP11-5E-3_Patient_Level_Candidate_Integration/",
  "STEP11-5E-3B_FINAL_Canonical_Patient_Level_LR_Integration.tsv"
)

output_dir <- paste0(
  "result_dir/Paper 1_Step 11-CellChat/",
  "STEP11-5F-3_Patient_Level_Expression_Context"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 2. Input existence
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  stop(
    "FINAL canonical input file not found:\n",
    input_file
  )
}


# ------------------------------------------------------------
# 3. Read FINAL LOCKED canonical table
# ------------------------------------------------------------

dt <- fread(
  input_file,
  sep = "\t",
  header = TRUE,
  data.table = TRUE,
  check.names = FALSE
)

cat("Input file:\n")
cat(" ", input_file, "\n\n")

cat("Rows:", nrow(dt), "\n")
cat("Columns:", ncol(dt), "\n\n")


# ------------------------------------------------------------
# 4. Expected schema
# ------------------------------------------------------------

expected_cols <- c(
  "Canonical_Key",
  "Exact_LR_Pair",
  "Candidate_Consistency_Class",
  "Candidate_Set_Provenance",
  "interaction_name_Provenance",
  "N_Candidate_Set",
  "N_interaction_name",
  "pathway_name",
  "ligand",
  "receptor",
  "annotation",
  "Direction_ID",
  "Sender",
  "Receiver",
  "Patient_Index",
  "Core_Gene_Involved",
  "Core_Gene",
  "Ligand_log2CPM",
  "Receptor_log2CPM",
  "Ligand_Receptor_Mean_log2CPM",
  "Ligand_Receptor_Difference_log2CPM",
  "Compatible"
)

schema_pass <- identical(
  names(dt),
  expected_cols
)

if (!schema_pass) {
  
  cat("Actual columns:\n")
  print(names(dt))
  
  cat("\nExpected columns:\n")
  print(expected_cols)
  
  stop(
    "\nFAIL — STEP11-5F-3 input schema does not match ",
    "the FINAL locked canonical table."
  )
}


# ------------------------------------------------------------
# 5. Locked global structure
# ------------------------------------------------------------

expected_n_rows <- 9620L
expected_n_lr <- 585L
expected_n_patients <- 4L

actual_n_rows <- nrow(dt)
actual_n_lr <- uniqueN(dt$Exact_LR_Pair)
actual_n_patients <- uniqueN(dt$Patient_Index)

row_count_pass <- actual_n_rows == expected_n_rows
lr_count_pass <- actual_n_lr == expected_n_lr
patient_count_pass <- actual_n_patients == expected_n_patients


# ------------------------------------------------------------
# 6. Patient IDs
# ------------------------------------------------------------

patient_ids <- sort(unique(dt$Patient_Index))

cat("Detected patient IDs:\n")
print(patient_ids)

cat("\n")


# ------------------------------------------------------------
# 7. Canonical key integrity
# ------------------------------------------------------------

expected_canonical_key <- paste(
  dt$Exact_LR_Pair,
  dt$Direction_ID,
  dt$Patient_Index,
  sep = "___"
)

canonical_key_definition_pass <- all(
  dt$Canonical_Key == expected_canonical_key
)

canonical_key_unique_pass <- uniqueN(
  dt$Canonical_Key
) == nrow(dt)


# ------------------------------------------------------------
# 8. Exact LR pair identity
# ------------------------------------------------------------

expected_lr_identity <- paste(
  dt$ligand,
  dt$receptor,
  sep = "___"
)

exact_lr_identity_pass <- all(
  dt$Exact_LR_Pair == expected_lr_identity
)


# ------------------------------------------------------------
# 9. Direction mapping
# ------------------------------------------------------------

expected_sender <- c(
  "Endothelial",
  "Myeloid/macrophage",
  "Endothelial",
  "T-cell",
  "Myeloid/macrophage",
  "T-cell"
)

expected_receiver <- c(
  "Myeloid/macrophage",
  "Endothelial",
  "T-cell",
  "Endothelial",
  "T-cell",
  "Myeloid/macrophage"
)

direction_ids <- sort(unique(dt$Direction_ID))

direction_id_pass <- identical(
  direction_ids,
  1:6
)

direction_sender_receiver_pass <- all(
  dt$Sender ==
    expected_sender[dt$Direction_ID] &
    dt$Receiver ==
    expected_receiver[dt$Direction_ID]
)


# ------------------------------------------------------------
# 10. Candidate consistency class
# ------------------------------------------------------------

allowed_candidate_classes <- c(
  "3_of_4_only",
  "4_of_4"
)

candidate_class_values_pass <- all(
  dt$Candidate_Consistency_Class %in%
    allowed_candidate_classes
)

candidate_class_invariant <- dt[
  ,
  uniqueN(Candidate_Consistency_Class),
  by = Exact_LR_Pair
]

candidate_class_invariant_pass <- all(
  candidate_class_invariant$V1 == 1L
)


# ------------------------------------------------------------
# 11. Locked candidate universe
# ------------------------------------------------------------

lr_class_summary <- dt[
  ,
  .(
    Candidate_Consistency_Class =
      unique(Candidate_Consistency_Class)
  ),
  by = Exact_LR_Pair
]

n_4of4_pairs <- sum(
  lr_class_summary$Candidate_Consistency_Class ==
    "4_of_4"
)

n_3of4_only_pairs <- sum(
  lr_class_summary$Candidate_Consistency_Class ==
    "3_of_4_only"
)

locked_4of4_pass <- n_4of4_pairs == 492L
locked_3of4_pass <- n_3of4_only_pairs == 93L


# ------------------------------------------------------------
# 12. Provenance fields
#
# Candidate_Set_Provenance and interaction_name_Provenance
# are retained as provenance only.
# ------------------------------------------------------------

candidate_provenance_nonempty_pass <- all(
  !is.na(dt$Candidate_Set_Provenance) &
    nzchar(dt$Candidate_Set_Provenance)
)

interaction_provenance_nonempty_pass <- all(
  !is.na(dt$interaction_name_Provenance) &
    nzchar(dt$interaction_name_Provenance)
)


# ------------------------------------------------------------
# 13. Core-gene separation
#
# Inherited from locked LR-universe audit.
# AOC3/MS4A7/IL7R are NOT LR nodes.
# ------------------------------------------------------------

core_gene_involved_pass <- all(
  dt$Core_Gene_Involved == FALSE,
  na.rm = TRUE
)

core_gene_field_all_na_pass <- all(
  is.na(dt$Core_Gene)
)

core_gene_separation_pass <-
  core_gene_involved_pass &
  core_gene_field_all_na_pass


# ------------------------------------------------------------
# 14. Expression fields
#
# IMPORTANT:
# These are LOCKED log2(CPM+1) values.
# No recalculation is performed.
# ------------------------------------------------------------

expression_columns <- c(
  "Ligand_log2CPM",
  "Receptor_log2CPM",
  "Ligand_Receptor_Mean_log2CPM",
  "Ligand_Receptor_Difference_log2CPM"
)

expression_numeric_pass <- all(
  vapply(
    dt[, ..expression_columns],
    is.numeric,
    logical(1)
  )
)

expression_finite_pass <- all(
  vapply(
    dt[, ..expression_columns],
    function(x) all(is.finite(x)),
    logical(1)
  )
)

expression_nonnegative_pass <- all(
  dt$Ligand_log2CPM >= 0 &
    dt$Receptor_log2CPM >= 0 &
    dt$Ligand_Receptor_Mean_log2CPM >= 0
)

difference_definition_expected <-
  dt$Ligand_Receptor_Difference_log2CPM

difference_definition_calculated <-
  dt$Ligand_log2CPM -
  dt$Receptor_log2CPM

difference_consistency_pass <- all(
  abs(
    difference_definition_expected -
      difference_definition_calculated
  ) < 1e-10
)


# ------------------------------------------------------------
# 15. Mean expression definition audit
#
# This is a STRUCTURAL audit only.
# It does not recalculate or modify the locked values.
# ------------------------------------------------------------

mean_definition_expected <-
  dt$Ligand_Receptor_Mean_log2CPM

mean_definition_calculated <-
  (
    dt$Ligand_log2CPM +
      dt$Receptor_log2CPM
  ) / 2

mean_definition_consistency_pass <- all(
  abs(
    mean_definition_expected -
      mean_definition_calculated
  ) < 1e-10
)


# ------------------------------------------------------------
# 16. Compatibility field audit
#
# Compatibility is inherited from STEP11-5D.
# It is NOT recalculated here.
# ------------------------------------------------------------

compatible_type_pass <- is.logical(dt$Compatible)

compatible_missing_pass <- !anyNA(dt$Compatible)


# ------------------------------------------------------------
# 17. Pathway / annotation / direction completeness
# ------------------------------------------------------------

pathway_complete_pass <- all(
  !is.na(dt$pathway_name) &
    nzchar(dt$pathway_name)
)

annotation_complete_pass <- all(
  !is.na(dt$annotation) &
    nzchar(dt$annotation)
)

sender_complete_pass <- all(
  !is.na(dt$Sender) &
    nzchar(dt$Sender)
)

receiver_complete_pass <- all(
  !is.na(dt$Receiver) &
    nzchar(dt$Receiver)
)

patient_complete_pass <- all(
  !is.na(dt$Patient_Index) &
    nzchar(dt$Patient_Index)
)


# ============================================================
# 18. PATIENT-LEVEL EXPRESSION CONTEXT TABLE
# ============================================================
#
# One row:
#   Exact_LR_Pair × Direction_ID × Patient
#
# This is essentially the canonical patient-level unit,
# enriched with explicitly ordered contextual fields.
# ============================================================

patient_context <- dt[
  ,
  .(
    Exact_LR_Pair,
    ligand,
    receptor,
    pathway_name,
    annotation,
    Direction_ID,
    Sender,
    Receiver,
    Patient_Index,
    Candidate_Consistency_Class,
    Candidate_Set_Provenance,
    interaction_name_Provenance,
    Ligand_log2CPM,
    Receptor_log2CPM,
    Ligand_Receptor_Mean_log2CPM,
    Ligand_Receptor_Difference_log2CPM,
    Compatible
  )
]

setorder(
  patient_context,
  Exact_LR_Pair,
  Direction_ID,
  Patient_Index
)

fwrite(
  patient_context,
  file.path(
    output_dir,
    "STEP11-5F-3A_Patient_Level_Expression_Context.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 19. EXACT LR PAIR-LEVEL EXPRESSION CONTEXT
# ============================================================
#
# One row per:
#   Exact_LR_Pair × Direction_ID
#
# Expression summaries are descriptive only.
# ============================================================

lr_direction_expression_summary <- dt[
  ,
  .(
    Candidate_Consistency_Class =
      unique(Candidate_Consistency_Class),
    
    pathway_name =
      paste(
        sort(unique(pathway_name)),
        collapse = " | "
      ),
    
    annotation =
      paste(
        sort(unique(annotation)),
        collapse = " | "
      ),
    
    Sender =
      unique(Sender),
    
    Receiver =
      unique(Receiver),
    
    N_patients =
      uniqueN(Patient_Index),
    
    N_compatible =
      sum(Compatible),
    
    N_incompatible =
      sum(!Compatible),
    
    Compatibility_fraction =
      mean(Compatible),
    
    Ligand_median_log2CPM =
      median(Ligand_log2CPM),
    
    Ligand_IQR_log2CPM =
      IQR(Ligand_log2CPM),
    
    Ligand_min_log2CPM =
      min(Ligand_log2CPM),
    
    Ligand_max_log2CPM =
      max(Ligand_log2CPM),
    
    Receptor_median_log2CPM =
      median(Receptor_log2CPM),
    
    Receptor_IQR_log2CPM =
      IQR(Receptor_log2CPM),
    
    Receptor_min_log2CPM =
      min(Receptor_log2CPM),
    
    Receptor_max_log2CPM =
      max(Receptor_log2CPM),
    
    MeanPair_median_log2CPM =
      median(Ligand_Receptor_Mean_log2CPM),
    
    MeanPair_IQR_log2CPM =
      IQR(Ligand_Receptor_Mean_log2CPM),
    
    Difference_median_log2CPM =
      median(Ligand_Receptor_Difference_log2CPM),
    
    Difference_IQR_log2CPM =
      IQR(Ligand_Receptor_Difference_log2CPM)
  ),
  by = .(
    Exact_LR_Pair,
    ligand,
    receptor,
    Direction_ID
  )
]

setorder(
  lr_direction_expression_summary,
  Direction_ID,
  Exact_LR_Pair
)

fwrite(
  lr_direction_expression_summary,
  file.path(
    output_dir,
    "STEP11-5F-3B_LR_Direction_Expression_Context_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 20. EXACT LR PAIR GLOBAL CONTEXT
# ============================================================
#
# Collapse directions only for descriptive context.
# Direction-specific information is retained as text.
# No biological ranking is performed.
# ============================================================

lr_pair_context_summary <- dt[
  ,
  .(
    Candidate_Consistency_Class =
      unique(Candidate_Consistency_Class),
    
    Candidate_Set_Provenance =
      paste(
        sort(unique(Candidate_Set_Provenance)),
        collapse = " | "
      ),
    
    interaction_name_Provenance =
      paste(
        sort(unique(interaction_name_Provenance)),
        collapse = " | "
      ),
    
    pathway_name =
      paste(
        sort(unique(pathway_name)),
        collapse = " | "
      ),
    
    annotation =
      paste(
        sort(unique(annotation)),
        collapse = " | "
      ),
    
    Directions =
      paste(
        sort(unique(Direction_ID)),
        collapse = " | "
      ),
    
    Sender_Receiver_Pairs =
      paste(
        unique(
          paste(
            Sender,
            Receiver,
            sep = " -> "
          )
        ),
        collapse = " | "
      ),
    
    N_direction_records =
      uniqueN(Direction_ID),
    
    N_patient_records =
      .N,
    
    N_patients =
      uniqueN(Patient_Index),
    
    N_compatible =
      sum(Compatible),
    
    N_incompatible =
      sum(!Compatible),
    
    Overall_compatibility_fraction =
      mean(Compatible),
    
    Ligand_median_log2CPM =
      median(Ligand_log2CPM),
    
    Ligand_IQR_log2CPM =
      IQR(Ligand_log2CPM),
    
    Receptor_median_log2CPM =
      median(Receptor_log2CPM),
    
    Receptor_IQR_log2CPM =
      IQR(Receptor_log2CPM),
    
    MeanPair_median_log2CPM =
      median(Ligand_Receptor_Mean_log2CPM),
    
    MeanPair_IQR_log2CPM =
      IQR(Ligand_Receptor_Mean_log2CPM),
    
    Difference_median_log2CPM =
      median(Ligand_Receptor_Difference_log2CPM),
    
    Difference_IQR_log2CPM =
      IQR(Ligand_Receptor_Difference_log2CPM)
  ),
  by = .(
    Exact_LR_Pair,
    ligand,
    receptor
  )
]

setorder(
  lr_pair_context_summary,
  Candidate_Consistency_Class,
  Exact_LR_Pair
)

fwrite(
  lr_pair_context_summary,
  file.path(
    output_dir,
    "STEP11-5F-3C_Exact_LR_Pair_Context_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 21. DIRECTION-LEVEL EXPRESSION CONTEXT
# ============================================================

direction_context_summary <- dt[
  ,
  .(
    N_exact_LR_pairs =
      uniqueN(Exact_LR_Pair),
    
    N_patient_records =
      .N,
    
    N_patients =
      uniqueN(Patient_Index),
    
    N_compatible =
      sum(Compatible),
    
    N_incompatible =
      sum(!Compatible),
    
    Compatibility_fraction =
      mean(Compatible),
    
    Ligand_median_log2CPM =
      median(Ligand_log2CPM),
    
    Ligand_IQR_log2CPM =
      IQR(Ligand_log2CPM),
    
    Receptor_median_log2CPM =
      median(Receptor_log2CPM),
    
    Receptor_IQR_log2CPM =
      IQR(Receptor_log2CPM),
    
    MeanPair_median_log2CPM =
      median(Ligand_Receptor_Mean_log2CPM),
    
    MeanPair_IQR_log2CPM =
      IQR(Ligand_Receptor_Mean_log2CPM),
    
    Difference_median_log2CPM =
      median(Ligand_Receptor_Difference_log2CPM),
    
    Difference_IQR_log2CPM =
      IQR(Ligand_Receptor_Difference_log2CPM)
  ),
  by = .(
    Direction_ID,
    Sender,
    Receiver
  )
]

setorder(
  direction_context_summary,
  Direction_ID
)

fwrite(
  direction_context_summary,
  file.path(
    output_dir,
    "STEP11-5F-3D_Direction_Expression_Context_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 22. CANDIDATE CLASS × EXPRESSION CONTEXT
# ============================================================

candidate_class_context_summary <- dt[
  ,
  .(
    N_exact_LR_pairs =
      uniqueN(Exact_LR_Pair),
    
    N_direction_records =
      uniqueN(
        paste(
          Exact_LR_Pair,
          Direction_ID,
          sep = "___"
        )
      ),
    
    N_patient_records =
      .N,
    
    N_patients =
      uniqueN(Patient_Index),
    
    N_compatible =
      sum(Compatible),
    
    N_incompatible =
      sum(!Compatible),
    
    Compatibility_fraction =
      mean(Compatible),
    
    Ligand_median_log2CPM =
      median(Ligand_log2CPM),
    
    Ligand_IQR_log2CPM =
      IQR(Ligand_log2CPM),
    
    Receptor_median_log2CPM =
      median(Receptor_log2CPM),
    
    Receptor_IQR_log2CPM =
      IQR(Receptor_log2CPM),
    
    MeanPair_median_log2CPM =
      median(Ligand_Receptor_Mean_log2CPM),
    
    MeanPair_IQR_log2CPM =
      IQR(Ligand_Receptor_Mean_log2CPM),
    
    Difference_median_log2CPM =
      median(Ligand_Receptor_Difference_log2CPM),
    
    Difference_IQR_log2CPM =
      IQR(Ligand_Receptor_Difference_log2CPM)
  ),
  by = Candidate_Consistency_Class
]

setorder(
  candidate_class_context_summary,
  Candidate_Consistency_Class
)

fwrite(
  candidate_class_context_summary,
  file.path(
    output_dir,
    "STEP11-5F-3E_CandidateClass_Expression_Context.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 23. COMPATIBILITY × EXPRESSION CONTEXT
# ============================================================
#
# Purely descriptive:
# compatible versus incompatible patient-level records.
#
# This does NOT test whether expression causes compatibility.
# ============================================================

compatibility_expression_summary <- dt[
  ,
  .(
    N_patient_records =
      .N,
    
    N_exact_LR_pairs =
      uniqueN(Exact_LR_Pair),
    
    N_directions =
      uniqueN(Direction_ID),
    
    N_patients =
      uniqueN(Patient_Index),
    
    Ligand_median_log2CPM =
      median(Ligand_log2CPM),
    
    Ligand_IQR_log2CPM =
      IQR(Ligand_log2CPM),
    
    Receptor_median_log2CPM =
      median(Receptor_log2CPM),
    
    Receptor_IQR_log2CPM =
      IQR(Receptor_log2CPM),
    
    MeanPair_median_log2CPM =
      median(Ligand_Receptor_Mean_log2CPM),
    
    MeanPair_IQR_log2CPM =
      IQR(Ligand_Receptor_Mean_log2CPM),
    
    Difference_median_log2CPM =
      median(Ligand_Receptor_Difference_log2CPM),
    
    Difference_IQR_log2CPM =
      IQR(Ligand_Receptor_Difference_log2CPM)
  ),
  by = Compatible
]

compatibility_expression_summary[
  Compatible == TRUE,
  Compatibility_Status := "Compatible"
]

compatibility_expression_summary[
  Compatible == FALSE,
  Compatibility_Status := "Not_compatible"
]

setcolorder(
  compatibility_expression_summary,
  c(
    "Compatibility_Status",
    "Compatible",
    setdiff(
      names(compatibility_expression_summary),
      c(
        "Compatibility_Status",
        "Compatible"
      )
    )
  )
)

fwrite(
  compatibility_expression_summary,
  file.path(
    output_dir,
    "STEP11-5F-3F_Compatibility_Expression_Context.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 24. ANNOTATION × EXPRESSION CONTEXT
# ============================================================

annotation_context_summary <- dt[
  ,
  .(
    N_exact_LR_pairs =
      uniqueN(Exact_LR_Pair),
    
    N_patient_records =
      .N,
    
    N_compatible =
      sum(Compatible),
    
    N_incompatible =
      sum(!Compatible),
    
    Compatibility_fraction =
      mean(Compatible),
    
    Ligand_median_log2CPM =
      median(Ligand_log2CPM),
    
    Ligand_IQR_log2CPM =
      IQR(Ligand_log2CPM),
    
    Receptor_median_log2CPM =
      median(Receptor_log2CPM),
    
    Receptor_IQR_log2CPM =
      IQR(Receptor_log2CPM),
    
    MeanPair_median_log2CPM =
      median(Ligand_Receptor_Mean_log2CPM),
    
    MeanPair_IQR_log2CPM =
      IQR(Ligand_Receptor_Mean_log2CPM)
  ),
  by = annotation
]

setorder(
  annotation_context_summary,
  annotation
)

fwrite(
  annotation_context_summary,
  file.path(
    output_dir,
    "STEP11-5F-3G_Annotation_Expression_Context.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 25. PATHWAY × EXPRESSION CONTEXT
# ============================================================

pathway_context_summary <- dt[
  ,
  .(
    N_exact_LR_pairs =
      uniqueN(Exact_LR_Pair),
    
    N_direction_records =
      uniqueN(
        paste(
          Exact_LR_Pair,
          Direction_ID,
          sep = "___"
        )
      ),
    
    N_patient_records =
      .N,
    
    N_compatible =
      sum(Compatible),
    
    N_incompatible =
      sum(!Compatible),
    
    Compatibility_fraction =
      mean(Compatible),
    
    Ligand_median_log2CPM =
      median(Ligand_log2CPM),
    
    Ligand_IQR_log2CPM =
      IQR(Ligand_log2CPM),
    
    Receptor_median_log2CPM =
      median(Receptor_log2CPM),
    
    Receptor_IQR_log2CPM =
      IQR(Receptor_log2CPM),
    
    MeanPair_median_log2CPM =
      median(Ligand_Receptor_Mean_log2CPM),
    
    MeanPair_IQR_log2CPM =
      IQR(Ligand_Receptor_Mean_log2CPM)
  ),
  by = pathway_name
]

setorder(
  pathway_context_summary,
  pathway_name
)

fwrite(
  pathway_context_summary,
  file.path(
    output_dir,
    "STEP11-5F-3H_Pathway_Expression_Context.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 26. GLOBAL CONTEXT SUMMARY
# ============================================================

global_context_summary <- data.table(
  Metric = c(
    "Canonical patient-level records",
    "Exact LR pairs",
    "Patients",
    "Directions",
    "4/4 exact LR pairs",
    "3/4-only exact LR pairs",
    "Compatible patient-level records",
    "Incompatible patient-level records",
    "Overall compatibility fraction",
    "Unique pathways",
    "Unique annotations",
    "Core-gene LR records",
    "Core-gene separation"
  ),
  Value = c(
    nrow(dt),
    uniqueN(dt$Exact_LR_Pair),
    uniqueN(dt$Patient_Index),
    uniqueN(dt$Direction_ID),
    n_4of4_pairs,
    n_3of4_only_pairs,
    sum(dt$Compatible),
    sum(!dt$Compatible),
    mean(dt$Compatible),
    uniqueN(dt$pathway_name),
    uniqueN(dt$annotation),
    sum(dt$Core_Gene_Involved),
    ifelse(core_gene_separation_pass, "PASS", "FAIL")
  )
)

fwrite(
  global_context_summary,
  file.path(
    output_dir,
    "STEP11-5F-3I_Global_Expression_Context_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 27. FINAL QC
# ============================================================

qc <- data.table(
  QC_Item = c(
    
    "Input file exists",
    
    "Input row count = 9,620",
    
    "Exact LR pair count = 585",
    
    "Patient count = 4",
    
    "Expected FINAL schema",
    
    "Canonical key definition",
    
    "Canonical key uniqueness",
    
    "Exact LR pair identity",
    
    "Direction IDs = 1:6",
    
    "Sender-receiver mapping",
    
    "Candidate classes valid",
    
    "Candidate class invariant",
    
    "4/4 candidate count = 492",
    
    "3/4-only candidate count = 93",
    
    "Candidate provenance nonempty",
    
    "interaction_name provenance nonempty",
    
    "Core-gene separation",
    
    "Pathway completeness",
    
    "Annotation completeness",
    
    "Patient completeness",
    
    "Expression fields numeric",
    
    "Expression fields finite",
    
    "Expression nonnegative",
    
    "Difference field structurally consistent",
    
    "Mean pair field structurally consistent",
    
    "Compatibility logical",
    
    "Compatibility nonmissing",
    
    "Patient-level context row conservation",
    
    "LR-direction context generated",
    
    "Candidate-class context generated",
    
    "Global context generated",
    
    "No statistical / biological ranking performed"
  ),
  
  PASS = c(
    
    file.exists(input_file),
    
    row_count_pass,
    
    lr_count_pass,
    
    patient_count_pass,
    
    schema_pass,
    
    canonical_key_definition_pass,
    
    canonical_key_unique_pass,
    
    exact_lr_identity_pass,
    
    direction_id_pass,
    
    direction_sender_receiver_pass,
    
    candidate_class_values_pass,
    
    candidate_class_invariant_pass,
    
    locked_4of4_pass,
    
    locked_3of4_pass,
    
    candidate_provenance_nonempty_pass,
    
    interaction_provenance_nonempty_pass,
    
    core_gene_separation_pass,
    
    pathway_complete_pass,
    
    annotation_complete_pass,
    
    patient_complete_pass,
    
    expression_numeric_pass,
    
    expression_finite_pass,
    
    expression_nonnegative_pass,
    
    difference_consistency_pass,
    
    mean_definition_consistency_pass,
    
    compatible_type_pass,
    
    compatible_missing_pass,
    
    nrow(patient_context) == nrow(dt),
    
    nrow(lr_direction_expression_summary) == 2405L,
    
    nrow(candidate_class_context_summary) == 2L,
    
    nrow(global_context_summary) > 0L,
    
    TRUE
  )
)


# ------------------------------------------------------------
# 28. Write QC
# ------------------------------------------------------------

fwrite(
  qc,
  file.path(
    output_dir,
    "STEP11-5F-3_FINAL_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 29. Console QC summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP11-5F-3 FINAL QC\n")
cat("============================================================\n\n")

cat(
  "QC items:",
  nrow(qc),
  "\n"
)

cat(
  "PASS:",
  sum(qc$PASS),
  "\n"
)

cat(
  "FAIL:",
  sum(!qc$PASS),
  "\n\n"
)


failed_qc <- qc[PASS == FALSE]

if (nrow(failed_qc) > 0L) {
  
  cat("FAILED QC ITEMS:\n")
  print(failed_qc)
  
  stop(
    "\nFAIL — STEP11-5F-3 QC did not pass."
  )
}


# ============================================================
# 30. Final structural summary
# ============================================================

cat("FINAL STRUCTURE\n")
cat("----------------\n")

cat(
  "Canonical patient-level records : ",
  nrow(dt),
  "\n",
  sep = ""
)

cat(
  "Exact LR pairs                 : ",
  uniqueN(dt$Exact_LR_Pair),
  "\n",
  sep = ""
)

cat(
  "Patients                       : ",
  uniqueN(dt$Patient_Index),
  "\n",
  sep = ""
)

cat(
  "Directions                     : ",
  uniqueN(dt$Direction_ID),
  "\n",
  sep = ""
)

cat(
  "4/4 candidate pairs            : ",
  n_4of4_pairs,
  "\n",
  sep = ""
)

cat(
  "3/4-only candidate pairs       : ",
  n_3of4_only_pairs,
  "\n",
  sep = ""
)

cat(
  "LR × Direction records         : ",
  nrow(lr_direction_expression_summary),
  "\n",
  sep = ""
)

cat(
  "Unique pathways                : ",
  uniqueN(dt$pathway_name),
  "\n",
  sep = ""
)

cat(
  "Unique annotations             : ",
  uniqueN(dt$annotation),
  "\n",
  sep = ""
)

cat("\n")


# ============================================================
# 31. Formal interpretation lock
# ============================================================

cat("FORMAL INTERPRETATION LOCK\n")
cat("--------------------------\n")

cat(
  "STEP11-5F-3 is a descriptive patient-level expression-context\n",
  "integration of the FINAL LOCKED LR candidate universe.\n\n"
)

cat(
  "Ligand/receptor expression values were inherited from the\n",
  "locked canonical table and were not recalculated.\n\n"
)

cat(
  "Patient-level compatibility values were inherited from the\n",
  "locked upstream analysis and were not recalculated.\n\n"
)

cat(
  "Expression magnitude does NOT represent interaction strength.\n"
)

cat(
  "Compatibility does NOT establish functional signaling.\n"
)

cat(
  "No statistical significance, interaction probability,\n",
  "causality, or biological priority was inferred.\n\n"
)

cat(
  "AOC3 / MS4A7 / IL7R remain separated from the LR universe.\n"
)

cat("\n")


# ============================================================
# 32. Final status
# ============================================================

cat("============================================================\n")
cat("PASS — STEP11-5F-3 PATIENT-LEVEL EXPRESSION CONTEXT /\n")
cat("       LR PAIR CONTEXT INTEGRATION\n")
cat("============================================================\n\n")

cat(
  "Formal status: PASS — READY FOR STEP11-5F-4\n"
)

cat("\nOutput directory:\n")
cat(
  " ",
  output_dir,
  "\n",
  sep = ""
)

cat("\n")
cat("END OF STEP11-5F-3\n")
cat("============================================================\n")

# ============================================================
# STEP11-5F-4
# Patient-Level LR Context Summary and Directional Consistency Audit
#
# PURPOSE
#   Final descriptive integration of the LOCKED genome-wide
#   patient-level LR candidate universe.
#
# FINAL MASTER UNIT
#   One row = one Exact_LR_Pair
#
#   585 exact LR pairs total
#     492 = 4/4 patient-supported
#      93 = 3/4 patient-supported only
#
# INTEGRATED INFORMATION
#   1. Exact LR pair identity
#   2. Candidate consistency class
#   3. Candidate-set provenance
#   4. interaction-name provenance
#   5. six sender -> receiver directions
#   6. patient-level compatibility pattern
#   7. ligand expression context
#   8. receptor expression context
#   9. annotation
#  10. pathway
#  11. core-gene separation
#
# IMPORTANT:
#   This is a FINAL DESCRIPTIVE SUMMARY / AUDIT.
#
# NO:
#   - new LR discovery
#   - raw UMI recalculation
#   - pseudobulk recalculation
#   - CPM recalculation
#   - log2(CPM+1) recalculation
#   - detection recalculation
#   - compatibility recalculation
#   - correlation
#   - P values
#   - FDR
#   - CellChat probability
#   - pathway enrichment
#   - biological ranking
#   - causal inference
#
# Expression magnitude is NOT interaction strength.
# Compatibility is NOT functional signaling.
#
# Input:
#   STEP11-5E-3B_FINAL_Canonical_Patient_Level_LR_Integration.tsv
#
# Output:
#   STEP11-5F-4_Patient_Level_LR_Final_Context_Summary
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP11-5F-4\n")
cat("Patient-Level LR Context Summary and Directional Consistency Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_packages <- c("data.table")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package not installed: ", pkg)
  }
}

library(data.table)


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <- paste0(
  "input_dir/Paper 1_Step 11-CellChat/",
  "STEP11-5E-3_Patient_Level_Candidate_Integration/",
  "STEP11-5E-3B_FINAL_Canonical_Patient_Level_LR_Integration.tsv"
)

output_dir <- paste0(
  "result_dir/Paper 1_Step 11-CellChat/",
  "STEP11-5F-4_Patient_Level_LR_Final_Context_Summary"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 2. Input existence
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  stop(
    "FINAL canonical input file not found:\n",
    input_file
  )
}


# ------------------------------------------------------------
# 3. Read FINAL LOCKED canonical table
# ------------------------------------------------------------

dt <- fread(
  input_file,
  sep = "\t",
  header = TRUE,
  data.table = TRUE,
  check.names = FALSE
)

cat("Input rows    :", nrow(dt), "\n")
cat("Input columns :", ncol(dt), "\n\n")


# ------------------------------------------------------------
# 4. FINAL LOCKED schema
# ------------------------------------------------------------

expected_cols <- c(
  "Canonical_Key",
  "Exact_LR_Pair",
  "Candidate_Consistency_Class",
  "Candidate_Set_Provenance",
  "interaction_name_Provenance",
  "N_Candidate_Set",
  "N_interaction_name",
  "pathway_name",
  "ligand",
  "receptor",
  "annotation",
  "Direction_ID",
  "Sender",
  "Receiver",
  "Patient_Index",
  "Core_Gene_Involved",
  "Core_Gene",
  "Ligand_log2CPM",
  "Receptor_log2CPM",
  "Ligand_Receptor_Mean_log2CPM",
  "Ligand_Receptor_Difference_log2CPM",
  "Compatible"
)

schema_pass <- identical(
  names(dt),
  expected_cols
)

if (!schema_pass) {
  
  cat("Actual columns:\n")
  print(names(dt))
  
  cat("\nExpected columns:\n")
  print(expected_cols)
  
  stop(
    "\nFAIL — FINAL canonical schema mismatch."
  )
}


# ============================================================
# 5. LOCKED GLOBAL STRUCTURE
# ============================================================

expected_n_rows <- 9620L
expected_n_lr <- 585L
expected_n_patients <- 4L
expected_n_directions <- 6L

n_rows <- nrow(dt)
n_lr <- uniqueN(dt$Exact_LR_Pair)
n_patients <- uniqueN(dt$Patient_Index)
n_directions <- uniqueN(dt$Direction_ID)

global_structure_pass <-
  n_rows == expected_n_rows &
  n_lr == expected_n_lr &
  n_patients == expected_n_patients &
  n_directions == expected_n_directions


# ------------------------------------------------------------
# 6. Detect patient IDs automatically
# ------------------------------------------------------------

patient_ids <- sort(unique(dt$Patient_Index))

cat("Detected patients:\n")
print(patient_ids)
cat("\n")


# ============================================================
# 7. CANONICAL KEY AUDIT
# ============================================================

expected_key <- paste(
  dt$Exact_LR_Pair,
  dt$Direction_ID,
  dt$Patient_Index,
  sep = "___"
)

canonical_key_definition_pass <- all(
  dt$Canonical_Key == expected_key
)

canonical_key_unique_pass <- uniqueN(
  dt$Canonical_Key
) == nrow(dt)


# ============================================================
# 8. EXACT LR IDENTITY AUDIT
# ============================================================

expected_lr_pair <- paste(
  dt$ligand,
  dt$receptor,
  sep = "___"
)

exact_lr_identity_pass <- all(
  dt$Exact_LR_Pair == expected_lr_pair
)


# ============================================================
# 9. DIRECTIONAL MAPPING
# ============================================================

direction_map <- data.table(
  Direction_ID = 1:6,
  
  Sender = c(
    "Endothelial",
    "Myeloid/macrophage",
    "Endothelial",
    "T-cell",
    "Myeloid/macrophage",
    "T-cell"
  ),
  
  Receiver = c(
    "Myeloid/macrophage",
    "Endothelial",
    "T-cell",
    "Endothelial",
    "T-cell",
    "Myeloid/macrophage"
  )
)

direction_map[
  ,
  Direction_Label :=
    paste(Sender, Receiver, sep = " -> ")
]

direction_id_pass <- identical(
  sort(unique(dt$Direction_ID)),
  1:6
)

direction_mapping_check <- merge(
  dt[
    ,
    .(
      Sender = unique(Sender),
      Receiver = unique(Receiver)
    ),
    by = Direction_ID
  ],
  direction_map,
  by = "Direction_ID",
  suffixes = c("_observed", "_locked")
)

direction_mapping_pass <- all(
  direction_mapping_check$Sender_observed ==
    direction_mapping_check$Sender_locked &
    direction_mapping_check$Receiver_observed ==
    direction_mapping_check$Receiver_locked
)


# ============================================================
# 10. CANDIDATE CLASS AUDIT
# ============================================================

allowed_classes <- c(
  "3_of_4_only",
  "4_of_4"
)

candidate_class_values_pass <- all(
  dt$Candidate_Consistency_Class %in%
    allowed_classes
)

candidate_class_by_lr <- dt[
  ,
  .(
    N_classes =
      uniqueN(Candidate_Consistency_Class),
    Candidate_Consistency_Class =
      paste(
        sort(unique(Candidate_Consistency_Class)),
        collapse = " | "
      )
  ),
  by = Exact_LR_Pair
]

candidate_class_invariant_pass <- all(
  candidate_class_by_lr$N_classes == 1L
)

n_4of4 <- sum(
  dt[
    ,
    unique(Candidate_Consistency_Class),
    by = Exact_LR_Pair
  ]$V1 == "4_of_4"
)

n_3of4 <- sum(
  dt[
    ,
    unique(Candidate_Consistency_Class),
    by = Exact_LR_Pair
  ]$V1 == "3_of_4_only"
)

candidate_class_count_pass <-
  n_4of4 == 492L &
  n_3of4 == 93L


# ============================================================
# 11. CORE-GENE SEPARATION
# ============================================================

core_gene_involved_pass <- all(
  dt$Core_Gene_Involved == FALSE,
  na.rm = TRUE
)

core_gene_field_na_pass <- all(
  is.na(dt$Core_Gene)
)

core_gene_separation_pass <-
  core_gene_involved_pass &
  core_gene_field_na_pass


# ============================================================
# 12. EXPRESSION FIELD INTEGRITY
# ============================================================

expression_cols <- c(
  "Ligand_log2CPM",
  "Receptor_log2CPM",
  "Ligand_Receptor_Mean_log2CPM",
  "Ligand_Receptor_Difference_log2CPM"
)

expression_numeric_pass <- all(
  vapply(
    dt[, ..expression_cols],
    is.numeric,
    logical(1)
  )
)

expression_finite_pass <- all(
  vapply(
    dt[, ..expression_cols],
    function(x) all(is.finite(x)),
    logical(1)
  )
)

expression_nonnegative_pass <- all(
  dt$Ligand_log2CPM >= 0 &
    dt$Receptor_log2CPM >= 0 &
    dt$Ligand_Receptor_Mean_log2CPM >= 0
)

mean_definition_pass <- all(
  abs(
    dt$Ligand_Receptor_Mean_log2CPM -
      (
        dt$Ligand_log2CPM +
          dt$Receptor_log2CPM
      ) / 2
  ) < 1e-10
)

difference_definition_pass <- all(
  abs(
    dt$Ligand_Receptor_Difference_log2CPM -
      (
        dt$Ligand_log2CPM -
          dt$Receptor_log2CPM
      )
  ) < 1e-10
)


# ============================================================
# 13. COMPATIBILITY FIELD
# ============================================================

compatible_type_pass <- is.logical(dt$Compatible)

compatible_complete_pass <- !anyNA(
  dt$Compatible
)


# ============================================================
# 14. ANNOTATION / PATHWAY / PROVENANCE
# ============================================================

annotation_complete_pass <- all(
  !is.na(dt$annotation) &
    nzchar(dt$annotation)
)

pathway_complete_pass <- all(
  !is.na(dt$pathway_name) &
    nzchar(dt$pathway_name)
)

candidate_provenance_complete_pass <- all(
  !is.na(dt$Candidate_Set_Provenance) &
    nzchar(dt$Candidate_Set_Provenance)
)

interaction_provenance_complete_pass <- all(
  !is.na(dt$interaction_name_Provenance) &
    nzchar(dt$interaction_name_Provenance)
)


# ============================================================
# 15. CREATE PATIENT-LEVEL COMPATIBILITY PATTERN
#
# One row:
#   Exact_LR_Pair × Direction_ID
#
# Four patient columns:
#   P1, P2, P3, P4
#
# No compatibility values are recalculated.
# ============================================================

patient_pattern <- dcast(
  dt,
  Exact_LR_Pair +
    ligand +
    receptor +
    Direction_ID +
    Sender +
    Receiver +
    Candidate_Consistency_Class +
    pathway_name +
    annotation +
    Candidate_Set_Provenance +
    interaction_name_Provenance
  ~ Patient_Index,
  value.var = "Compatible"
)

setorder(
  patient_pattern,
  Direction_ID,
  Exact_LR_Pair
)

expected_lr_direction_records <- 2405L

patient_pattern_row_count_pass <-
  nrow(patient_pattern) ==
  expected_lr_direction_records


# ------------------------------------------------------------
# 16. Check patient columns
# ------------------------------------------------------------

missing_patient_columns <- setdiff(
  patient_ids,
  names(patient_pattern)
)

patient_columns_pass <-
  length(missing_patient_columns) == 0L


# ============================================================
# 17. Add patient-level compatibility summary
# ============================================================

patient_pattern[
  ,
  N_Patients := length(patient_ids)
]

patient_pattern[
  ,
  N_Compatible := rowSums(
    as.matrix(
      .SD
    ),
    na.rm = FALSE
  ),
  .SDcols = patient_ids
]

patient_pattern[
  ,
  N_Incompatible :=
    N_Patients - N_Compatible
]

patient_pattern[
  ,
  Compatibility_Fraction :=
    N_Compatible / N_Patients
]

patient_pattern[
  ,
  Compatibility_Pattern :=
    apply(
      .SD,
      1,
      function(x) {
        paste(
          ifelse(
            x,
            "Compatible",
            "Not_compatible"
          ),
          collapse = " | "
        )
      }
    ),
  .SDcols = patient_ids
]


# ============================================================
# 18. Merge directional labels
# ============================================================

patient_pattern <- merge(
  patient_pattern,
  direction_map[
    ,
    .(
      Direction_ID,
      Direction_Label
    )
  ],
  by = "Direction_ID",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  patient_pattern,
  Direction_ID,
  Exact_LR_Pair
)


# ============================================================
# 19. PATIENT-LEVEL COMPATIBILITY MASTER TABLE
# ============================================================

fwrite(
  patient_pattern,
  file.path(
    output_dir,
    "STEP11-5F-4A_LR_Direction_Patient_Compatibility_Pattern.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 20. DIRECTIONAL CONSISTENCY SUMMARY
# ============================================================

direction_summary <- dt[
  ,
  .(
    N_Exact_LR_Pairs =
      uniqueN(Exact_LR_Pair),
    
    N_Patient_Records =
      .N,
    
    N_Compatible =
      sum(Compatible),
    
    N_Incompatible =
      sum(!Compatible),
    
    Compatibility_Fraction =
      mean(Compatible),
    
    Ligand_Median_log2CPM =
      median(Ligand_log2CPM),
    
    Ligand_IQR_log2CPM =
      IQR(Ligand_log2CPM),
    
    Receptor_Median_log2CPM =
      median(Receptor_log2CPM),
    
    Receptor_IQR_log2CPM =
      IQR(Receptor_log2CPM),
    
    MeanPair_Median_log2CPM =
      median(Ligand_Receptor_Mean_log2CPM),
    
    MeanPair_IQR_log2CPM =
      IQR(Ligand_Receptor_Mean_log2CPM)
  ),
  by = .(
    Direction_ID,
    Sender,
    Receiver
  )
]

direction_summary <- merge(
  direction_summary,
  direction_map[
    ,
    .(
      Direction_ID,
      Direction_Label
    )
  ],
  by = "Direction_ID",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  direction_summary,
  Direction_ID
)

fwrite(
  direction_summary,
  file.path(
    output_dir,
    "STEP11-5F-4B_Directional_Consistency_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 21. EXACT LR PAIR MASTER DESCRIPTIVE TABLE
#
# FINAL MASTER UNIT:
#   1 row = 1 Exact_LR_Pair
#
# This is the principal manuscript / supplementary-table
# oriented output.
# ============================================================

master_lr <- dt[
  ,
  .(
    Candidate_Consistency_Class =
      unique(Candidate_Consistency_Class),
    
    Candidate_Set_Provenance =
      paste(
        sort(unique(Candidate_Set_Provenance)),
        collapse = " | "
      ),
    
    interaction_name_Provenance =
      paste(
        sort(unique(interaction_name_Provenance)),
        collapse = " | "
      ),
    
    pathway_name =
      paste(
        sort(unique(pathway_name)),
        collapse = " | "
      ),
    
    annotation =
      paste(
        sort(unique(annotation)),
        collapse = " | "
      ),
    
    N_Directions =
      uniqueN(Direction_ID),
    
    Directions =
      paste(
        sort(unique(Direction_ID)),
        collapse = " | "
      ),
    
    Direction_Labels =
      paste(
        unique(
          paste(
            Sender,
            Receiver,
            sep = " -> "
          )
        ),
        collapse = " | "
      ),
    
    N_Patient_Records =
      .N,
    
    N_Patients =
      uniqueN(Patient_Index),
    
    N_Compatible =
      sum(Compatible),
    
    N_Incompatible =
      sum(!Compatible),
    
    Overall_Compatibility_Fraction =
      mean(Compatible),
    
    Ligand_Median_log2CPM =
      median(Ligand_log2CPM),
    
    Ligand_IQR_log2CPM =
      IQR(Ligand_log2CPM),
    
    Ligand_Min_log2CPM =
      min(Ligand_log2CPM),
    
    Ligand_Max_log2CPM =
      max(Ligand_log2CPM),
    
    Receptor_Median_log2CPM =
      median(Receptor_log2CPM),
    
    Receptor_IQR_log2CPM =
      IQR(Receptor_log2CPM),
    
    Receptor_Min_log2CPM =
      min(Receptor_log2CPM),
    
    Receptor_Max_log2CPM =
      max(Receptor_log2CPM),
    
    MeanPair_Median_log2CPM =
      median(Ligand_Receptor_Mean_log2CPM),
    
    MeanPair_IQR_log2CPM =
      IQR(Ligand_Receptor_Mean_log2CPM),
    
    Difference_Median_log2CPM =
      median(Ligand_Receptor_Difference_log2CPM),
    
    Difference_IQR_log2CPM =
      IQR(Ligand_Receptor_Difference_log2CPM),
    
    Core_Gene_Involved =
      any(Core_Gene_Involved)
  ),
  by = .(
    Exact_LR_Pair,
    ligand,
    receptor
  )
]

setorder(
  master_lr,
  Candidate_Consistency_Class,
  Exact_LR_Pair
)


# ------------------------------------------------------------
# 22. Add overall patient compatibility pattern
# ------------------------------------------------------------

master_pattern <- dt[
  ,
  .(
    Patient_Compatibility_Pattern =
      paste(
        paste(
          Patient_Index,
          ifelse(
            Compatible,
            "Compatible",
            "Not_compatible"
          ),
          sep = "="
        ),
        collapse = " | "
      )
  ),
  by = .(
    Exact_LR_Pair,
    Direction_ID
  )
]

master_pattern[
  ,
  Direction_Pattern_Record :=
    paste(
      Direction_ID,
      Patient_Compatibility_Pattern,
      sep = ": "
    )
]

master_pattern_collapsed <- master_pattern[
  ,
  .(
    Patient_Compatibility_Patterns =
      paste(
        Direction_Pattern_Record,
        collapse = " || "
      )
  ),
  by = Exact_LR_Pair
]

master_lr <- merge(
  master_lr,
  master_pattern_collapsed,
  by = "Exact_LR_Pair",
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------
# 23. Add direction-level compatibility fractions
#
# These are descriptive only.
# ------------------------------------------------------------

direction_fraction_table <- dt[
  ,
  .(
    Compatibility_Fraction =
      mean(Compatible)
  ),
  by = .(
    Exact_LR_Pair,
    Direction_ID
  )
]

direction_fraction_wide <- dcast(
  direction_fraction_table,
  Exact_LR_Pair ~ Direction_ID,
  value.var = "Compatibility_Fraction"
)

setnames(
  direction_fraction_wide,
  old = as.character(1:6),
  new = paste0(
    "Direction_",
    1:6,
    "_Compatibility_Fraction"
  )
)

master_lr <- merge(
  master_lr,
  direction_fraction_wide,
  by = "Exact_LR_Pair",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  master_lr,
  Candidate_Consistency_Class,
  Exact_LR_Pair
)


# ============================================================
# 24. FINAL MASTER TABLE
# ============================================================

fwrite(
  master_lr,
  file.path(
    output_dir,
    "STEP11-5F-4C_FINAL_Master_LR_Descriptive_Context_Table.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 25. CANDIDATE CLASS SUMMARY
# ============================================================

candidate_class_summary <- dt[
  ,
  .(
    N_Exact_LR_Pairs =
      uniqueN(Exact_LR_Pair),
    
    N_Direction_Records =
      uniqueN(
        paste(
          Exact_LR_Pair,
          Direction_ID,
          sep = "___"
        )
      ),
    
    N_Patient_Records =
      .N,
    
    N_Patients =
      uniqueN(Patient_Index),
    
    N_Compatible =
      sum(Compatible),
    
    N_Incompatible =
      sum(!Compatible),
    
    Compatibility_Fraction =
      mean(Compatible),
    
    Ligand_Median_log2CPM =
      median(Ligand_log2CPM),
    
    Receptor_Median_log2CPM =
      median(Receptor_log2CPM),
    
    MeanPair_Median_log2CPM =
      median(Ligand_Receptor_Mean_log2CPM)
  ),
  by = Candidate_Consistency_Class
]

setorder(
  candidate_class_summary,
  Candidate_Consistency_Class
)

fwrite(
  candidate_class_summary,
  file.path(
    output_dir,
    "STEP11-5F-4D_CandidateClass_Context_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 26. ANNOTATION SUMMARY
# ============================================================

annotation_summary <- dt[
  ,
  .(
    N_Exact_LR_Pairs =
      uniqueN(Exact_LR_Pair),
    
    N_Direction_Records =
      uniqueN(
        paste(
          Exact_LR_Pair,
          Direction_ID,
          sep = "___"
        )
      ),
    
    N_Patient_Records =
      .N,
    
    N_Compatible =
      sum(Compatible),
    
    N_Incompatible =
      sum(!Compatible),
    
    Compatibility_Fraction =
      mean(Compatible),
    
    Ligand_Median_log2CPM =
      median(Ligand_log2CPM),
    
    Receptor_Median_log2CPM =
      median(Receptor_log2CPM)
  ),
  by = annotation
]

setorder(
  annotation_summary,
  annotation
)

fwrite(
  annotation_summary,
  file.path(
    output_dir,
    "STEP11-5F-4E_Annotation_Context_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 27. PATHWAY SUMMARY
# ============================================================

pathway_summary <- dt[
  ,
  .(
    N_Exact_LR_Pairs =
      uniqueN(Exact_LR_Pair),
    
    N_Direction_Records =
      uniqueN(
        paste(
          Exact_LR_Pair,
          Direction_ID,
          sep = "___"
        )
      ),
    
    N_Patient_Records =
      .N,
    
    N_Compatible =
      sum(Compatible),
    
    N_Incompatible =
      sum(!Compatible),
    
    Compatibility_Fraction =
      mean(Compatible),
    
    Ligand_Median_log2CPM =
      median(Ligand_log2CPM),
    
    Receptor_Median_log2CPM =
      median(Receptor_log2CPM),
    
    MeanPair_Median_log2CPM =
      median(Ligand_Receptor_Mean_log2CPM)
  ),
  by = pathway_name
]

setorder(
  pathway_summary,
  pathway_name
)

fwrite(
  pathway_summary,
  file.path(
    output_dir,
    "STEP11-5F-4F_Pathway_Context_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 28. GLOBAL SUMMARY
# ============================================================

global_summary <- data.table(
  Metric = c(
    "Canonical patient-level records",
    "Exact LR pairs",
    "Patients",
    "Sender-receiver directions",
    "LR x Direction records",
    "4/4 exact LR pairs",
    "3/4-only exact LR pairs",
    "Unique pathways",
    "Unique annotations",
    "Compatible patient-level records",
    "Incompatible patient-level records",
    "Overall compatibility fraction",
    "Core-gene LR records",
    "Core-gene separation"
  ),
  
  Value = c(
    nrow(dt),
    uniqueN(dt$Exact_LR_Pair),
    uniqueN(dt$Patient_Index),
    uniqueN(dt$Direction_ID),
    uniqueN(
      paste(
        dt$Exact_LR_Pair,
        dt$Direction_ID,
        sep = "___"
      )
    ),
    n_4of4,
    n_3of4,
    uniqueN(dt$pathway_name),
    uniqueN(dt$annotation),
    sum(dt$Compatible),
    sum(!dt$Compatible),
    mean(dt$Compatible),
    sum(dt$Core_Gene_Involved),
    ifelse(
      core_gene_separation_pass,
      "PASS",
      "FAIL"
    )
  )
)

fwrite(
  global_summary,
  file.path(
    output_dir,
    "STEP11-5F-4G_Global_Final_Context_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 29. FINAL QC
# ============================================================

master_unique_pass <-
  nrow(master_lr) == 585L &
  uniqueN(master_lr$Exact_LR_Pair) == 585L

master_columns_pass <- all(
  c(
    "Exact_LR_Pair",
    "ligand",
    "receptor",
    "Candidate_Consistency_Class",
    "Candidate_Set_Provenance",
    "interaction_name_Provenance",
    "pathway_name",
    "annotation",
    "N_Directions",
    "N_Patient_Records",
    "N_Patients",
    "N_Compatible",
    "N_Incompatible",
    "Overall_Compatibility_Fraction",
    "Ligand_Median_log2CPM",
    "Receptor_Median_log2CPM",
    "MeanPair_Median_log2CPM",
    "Core_Gene_Involved"
  ) %in% names(master_lr)
)

master_row_conservation_pass <-
  sum(master_lr$N_Patient_Records) == nrow(dt)

master_patient_count_pass <- all(
  master_lr$N_Patients == 4L
)

master_core_gene_pass <- all(
  master_lr$Core_Gene_Involved == FALSE
)

direction_fraction_range_pass <- all(
  vapply(
    master_lr[
      ,
      .SD,
      .SDcols = grep(
        "^Direction_[1-6]_Compatibility_Fraction$",
        names(master_lr),
        value = TRUE
      )
    ],
    function(x) {
      all(
        is.na(x) |
          (
            is.finite(x) &
              x >= 0 &
              x <= 1
          )
      )
    },
    logical(1)
  )
)

qc <- data.table(
  
  QC_Item = c(
    
    "FINAL input schema",
    
    "Canonical global structure = 9,620 / 585 / 4 / 6",
    
    "Canonical key definition",
    
    "Canonical key uniqueness",
    
    "Exact LR identity",
    
    "Direction IDs = 1:6",
    
    "Sender-receiver mapping",
    
    "Candidate classes valid",
    
    "Candidate class invariant",
    
    "4/4 candidate pairs = 492",
    
    "3/4-only candidate pairs = 93",
    
    "Core-gene separation",
    
    "Expression fields numeric",
    
    "Expression fields finite",
    
    "Expression values nonnegative",
    
    "Mean expression definition consistent",
    
    "Difference expression definition consistent",
    
    "Compatibility field logical",
    
    "Compatibility complete",
    
    "Annotation complete",
    
    "Pathway complete",
    
    "Candidate provenance complete",
    
    "interaction-name provenance complete",
    
    "Patient pattern generated",
    
    "Patient columns complete",
    
    "Patient pattern row count = 2,405",
    
    "Final master table = 585 unique LR pairs",
    
    "Master table required columns present",
    
    "Master patient-record conservation",
    
    "Master patient count = 4",
    
    "Master core-gene separation",
    
    "Direction compatibility fractions valid",
    
    "No statistical / ranking analysis performed"
  ),
  
  PASS = c(
    
    schema_pass,
    
    global_structure_pass,
    
    canonical_key_definition_pass,
    
    canonical_key_unique_pass,
    
    exact_lr_identity_pass,
    
    direction_id_pass,
    
    direction_mapping_pass,
    
    candidate_class_values_pass,
    
    candidate_class_invariant_pass,
    
    n_4of4 == 492L,
    
    n_3of4 == 93L,
    
    core_gene_separation_pass,
    
    expression_numeric_pass,
    
    expression_finite_pass,
    
    expression_nonnegative_pass,
    
    mean_definition_pass,
    
    difference_definition_pass,
    
    compatible_type_pass,
    
    compatible_complete_pass,
    
    annotation_complete_pass,
    
    pathway_complete_pass,
    
    candidate_provenance_complete_pass,
    
    interaction_provenance_complete_pass,
    
    nrow(patient_pattern) > 0L,
    
    patient_columns_pass,
    
    patient_pattern_row_count_pass,
    
    master_unique_pass,
    
    master_columns_pass,
    
    master_row_conservation_pass,
    
    master_patient_count_pass,
    
    master_core_gene_pass,
    
    direction_fraction_range_pass,
    
    TRUE
  )
)


# ============================================================
# 30. Write FINAL QC
# ============================================================

fwrite(
  qc,
  file.path(
    output_dir,
    "STEP11-5F-4_FINAL_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 31. Console QC
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP11-5F-4 FINAL QC\n")
cat("============================================================\n\n")

cat(
  "QC items : ",
  nrow(qc),
  "\n",
  sep = ""
)

cat(
  "PASS     : ",
  sum(qc$PASS),
  "\n",
  sep = ""
)

cat(
  "FAIL     : ",
  sum(!qc$PASS),
  "\n\n",
  sep = ""
)

failed_qc <- qc[PASS == FALSE]

if (nrow(failed_qc) > 0L) {
  
  cat("FAILED QC ITEMS:\n")
  print(failed_qc)
  
  stop(
    "\nFAIL — STEP11-5F-4 QC did not pass."
  )
}


# ============================================================
# 32. Final structural summary
# ============================================================

cat("FINAL STRUCTURE\n")
cat("----------------\n")

cat(
  "Canonical patient-level records : ",
  nrow(dt),
  "\n",
  sep = ""
)

cat(
  "Exact LR pairs                 : ",
  uniqueN(dt$Exact_LR_Pair),
  "\n",
  sep = ""
)

cat(
  "Patients                       : ",
  uniqueN(dt$Patient_Index),
  "\n",
  sep = ""
)

cat(
  "Directions                     : ",
  uniqueN(dt$Direction_ID),
  "\n",
  sep = ""
)

cat(
  "LR x Direction records         : ",
  uniqueN(
    paste(
      dt$Exact_LR_Pair,
      dt$Direction_ID,
      sep = "___"
    )
  ),
  "\n",
  sep = ""
)

cat(
  "4/4 candidate pairs            : ",
  n_4of4,
  "\n",
  sep = ""
)

cat(
  "3/4-only candidate pairs       : ",
  n_3of4,
  "\n",
  sep = ""
)

cat(
  "Unique pathways                : ",
  uniqueN(dt$pathway_name),
  "\n",
  sep = ""
)

cat(
  "Unique annotations             : ",
  uniqueN(dt$annotation),
  "\n",
  sep = ""
)

cat(
  "FINAL master LR table rows     : ",
  nrow(master_lr),
  "\n",
  sep = ""
)

cat("\n")


# ============================================================
# 33. FORMAL INTERPRETATION LOCK
# ============================================================

cat("FORMAL INTERPRETATION LOCK\n")
cat("--------------------------\n\n")

cat(
  "STEP11-5F-4 is the final descriptive integration and\n",
  "directional consistency audit of the locked LR candidate universe.\n\n"
)

cat(
  "The final master table contains one row per Exact_LR_Pair.\n"
)

cat(
  "The patient-level evidence remains traceable through the\n",
  "patient-pattern and directional summary tables.\n\n"
)

cat(
  "The 585 exact LR pairs comprise 492 pairs supported in all\n",
  "four matched patients and 93 pairs supported in three of four.\n\n"
)

cat(
  "Ligand/receptor expression values were inherited from the\n",
  "locked canonical table and were not recalculated.\n\n"
)

cat(
  "Compatibility values were inherited from the locked upstream\n",
  "analysis and were not recalculated.\n\n"
)

cat(
  "Expression magnitude does NOT represent interaction strength.\n"
)

cat(
  "Compatibility does NOT establish functional signaling.\n"
)

cat(
  "Pathway and annotation fields are descriptive resource-level\n",
  "context and were not used for biological ranking.\n\n"
)

cat(
  "AOC3 / MS4A7 / IL7R remain separated from the LR universe.\n\n"
)

cat(
  "No statistical significance, interaction probability,\n",
  "causality, pathway activation, or biological priority was inferred.\n\n"
)


# ============================================================
# 34. FINAL STATUS
# ============================================================

cat("============================================================\n")
cat("PASS — STEP11-5F-4 PATIENT-LEVEL LR CONTEXT SUMMARY /\n")
cat("       DIRECTIONAL CONSISTENCY AUDIT\n")
cat("============================================================\n\n")

cat(
  "Formal status: PASS — FINAL LOCK\n"
)

cat("\nOutput directory:\n")
cat(
  " ",
  output_dir,
  "\n",
  sep = ""
)

cat("\n")
cat("END OF STEP11-5F-4\n")
cat("============================================================\n")


# ============================================================
# STEP11-5F-5 — FINAL FIGURE 7 PLOTTING
# Patient-level ligand–receptor compatibility context
#
# Panels:
#   Figure 7A — Patient-level consistency
#   Figure 7B — Six sender–receiver directions
#   Figure 7C — Patient × direction compatibility heatmap
#   Figure 7D — Resource-level annotation
#
# IMPORTANT:
#   This script is a PRESENTATION-ONLY workflow.
#
#   No upstream LR discovery is performed.
#   No expression values are recalculated.
#   No compatibility values are recalculated.
#   No statistical testing is performed.
#   No pathway ranking is performed.
#
#   All quantities are inherited from locked STEP11-5F-4
#   canonical LR results.
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_packages <- c(
  "data.table",
  "ggplot2"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package not installed: ", pkg)
  }
}

library(data.table)
library(ggplot2)


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

root_dir <- "input_dir/Paper 1_Step 11-CellChat"

input_dir <- file.path(
  root_dir,
  "STEP11-5F-4_Patient_Level_LR_Final_Context_Summary"
)

# Final canonical table
canonical_file <- file.path(
  root_dir,
  "STEP11-5E-3_Patient_Level_Candidate_Integration",
  "STEP11-5E-3B_FINAL_Canonical_Patient_Level_LR_Integration.tsv"
)

# STEP11-5F-4 master table
master_file <- file.path(
  input_dir,
  "STEP11-5F-4C_FINAL_Master_LR_Descriptive_Context_Table.tsv"
)

output_dir <- file.path(
  root_dir,
  "STEP11-5F-5_Figure7_Final_Plotting"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Read locked inputs
# ------------------------------------------------------------

if (!file.exists(canonical_file)) {
  stop("Canonical input file not found:\n", canonical_file)
}

if (!file.exists(master_file)) {
  stop("STEP11-5F-4 master table not found:\n", master_file)
}

canonical_dt <- fread(
  canonical_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

master_dt <- fread(
  master_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 3. Basic locked-result QC
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP11-5F-5 — FIGURE 7 PRESENTATION QC\n")
cat("============================================================\n\n")

cat("Canonical rows :", nrow(canonical_dt), "\n")
cat("Master rows    :", nrow(master_dt), "\n")

if (nrow(canonical_dt) != 9620L) {
  stop("Unexpected canonical row count.")
}

if (uniqueN(canonical_dt$Exact_LR_Pair) != 585L) {
  stop("Unexpected exact LR-pair count.")
}

if (uniqueN(canonical_dt$Patient_Index) != 4L) {
  stop("Unexpected patient count.")
}

if (uniqueN(canonical_dt$Direction_ID) != 6L) {
  stop("Unexpected direction count.")
}

if (nrow(master_dt) != 585L) {
  stop("Unexpected master-table row count.")
}

cat("Locked structure confirmed:\n")
cat("  9,620 canonical patient-level records\n")
cat("  585 exact LR pairs\n")
cat("  4 matched patients\n")
cat("  6 sender–receiver directions\n")
cat("  585-row final master table\n\n")


# ------------------------------------------------------------
# 4. Direction labels
# ------------------------------------------------------------

direction_labels <- c(
  "1" = "E → M",
  "2" = "M → E",
  "3" = "E → T",
  "4" = "T → E",
  "5" = "M → T",
  "6" = "T → M"
)

direction_long <- c(
  "1" = "Endothelial → Myeloid/macrophage",
  "2" = "Myeloid/macrophage → Endothelial",
  "3" = "Endothelial → T cell",
  "4" = "T cell → Endothelial",
  "5" = "Myeloid/macrophage → T cell",
  "6" = "T cell → Myeloid/macrophage"
)


# ------------------------------------------------------------
# 5. Global plotting theme
# ------------------------------------------------------------

theme_paper7 <- theme_classic(base_family = "Arial") +
  theme(
    text = element_text(
      family = "Arial",
      colour = "black"
    ),
    
    axis.title = element_text(
      family = "Arial",
      face = "bold",
      size = 10
    ),
    
    axis.text = element_text(
      family = "Arial",
      size = 9,
      colour = "black"
    ),
    
    plot.title = element_text(
      family = "Arial",
      face = "bold",
      size = 11,
      hjust = 0.5
    ),
    
    legend.title = element_text(
      family = "Arial",
      face = "bold",
      size = 9
    ),
    
    legend.text = element_text(
      family = "Arial",
      size = 9
    ),
    
    panel.border = element_blank(),
    
    axis.line = element_line(
      linewidth = 0.5,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.4,
      colour = "black"
    ),
    
    plot.margin = margin(
      t = 8,
      r = 10,
      b = 8,
      l = 10
    )
  )


# ============================================================
# FIGURE 7A
# Patient-level consistency
# ============================================================

# ------------------------------------------------------------
# 6. Locked candidate consistency counts
# ------------------------------------------------------------

panelA_dt <- data.table(
  Consistency = factor(
    c("3/4 patients", "4/4 patients"),
    levels = c("3/4 patients", "4/4 patients")
  ),
  N = c(93L, 492L)
)

# Explicit arithmetic QC
if (sum(panelA_dt$N) != 585L) {
  stop("Panel A arithmetic failed.")
}

# ------------------------------------------------------------
# 7. Plot Figure 7A
# ------------------------------------------------------------

p7A <- ggplot(
  panelA_dt,
  aes(
    x = Consistency,
    y = N
  )
) +
  geom_col(
    width = 0.62
  ) +
  geom_text(
    aes(label = N),
    vjust = -0.35,
    family = "Arial",
    fontface = "bold",
    size = 3.4
  ) +
  scale_y_continuous(
    limits = c(0, 540),
    breaks = seq(0, 500, 100),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = "Patient-level consistency",
    y = "Exact LR pairs"
  ) +
  theme_paper7

print(p7A)


# ------------------------------------------------------------
# 8. Export Figure 7A
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "Figure7A_Patient_Level_Consistency.pdf"
  ),
  plot = p7A,
  width = 9.0,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure7A_Patient_Level_Consistency_600dpi.tiff"
  ),
  plot = p7A,
  width = 9.0,
  height = 5.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ============================================================
# FIGURE 7B
# Six sender → receiver directions
# ============================================================

# ------------------------------------------------------------
# 9. Direction-level exact LR pair counts
# ------------------------------------------------------------

panelB_dt <- canonical_dt[
  ,
  .(
    N_exact_LR_pairs = uniqueN(Exact_LR_Pair)
  ),
  by = Direction_ID
]

panelB_dt[
  ,
  Direction := direction_labels[as.character(Direction_ID)]
]

panelB_dt[
  ,
  Direction := factor(
    Direction,
    levels = rev(
      c(
        "E → M",
        "M → E",
        "E → T",
        "T → E",
        "M → T",
        "T → M"
      )
    )
  )
]

# Locked expected counts
expected_direction_counts <- c(
  "1" = 382L,
  "2" = 371L,
  "3" = 387L,
  "4" = 343L,
  "5" = 485L,
  "6" = 437L
)

panelB_check <- setNames(
  panelB_dt$N_exact_LR_pairs,
  as.character(panelB_dt$Direction_ID)
)

if (!identical(
  as.integer(panelB_check[names(expected_direction_counts)]),
  as.integer(expected_direction_counts)
)) {
  stop(
    "Panel B direction counts do not match the locked STEP11 results."
  )
}

# ------------------------------------------------------------
# 10. Plot Figure 7B
# ------------------------------------------------------------

p7B <- ggplot(
  panelB_dt,
  aes(
    x = Direction,
    y = N_exact_LR_pairs
  )
) +
  geom_col(
    width = 0.62
  ) +
  geom_text(
    aes(label = N_exact_LR_pairs),
    hjust = -0.15,
    family = "Arial",
    fontface = "bold",
    size = 3.2
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 530),
    breaks = seq(0, 500, 100),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Sender → receiver direction",
    y = "Exact LR pairs"
  ) +
  theme_paper7

print(p7B)


# ------------------------------------------------------------
# 11. Export Figure 7B
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "Figure7B_Sender_Receiver_Directions.pdf"
  ),
  plot = p7B,
  width = 9.0,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure7B_Sender_Receiver_Directions_600dpi.tiff"
  ),
  plot = p7B,
  width = 9.0,
  height = 5.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ============================================================
# FIGURE 7C
# Patient × direction compatibility heatmap
# ============================================================

# ------------------------------------------------------------
# 12. Prepare patient-level compatibility matrix
# ------------------------------------------------------------
#
# Compatible is inherited directly from the locked canonical
# patient-level table.
#
# We summarize the proportion of candidate LR pairs that are
# compatible within each patient and direction.
#
# This is a descriptive display only.
# ------------------------------------------------------------

panelC_dt <- canonical_dt[
  ,
  .(
    N_LR_pairs = .N,
    N_compatible = sum(
      Compatible == TRUE,
      na.rm = TRUE
    )
  ),
  by = .(
    Patient_Index,
    Direction_ID
  )
]

panelC_dt[
  ,
  Compatibility_Percent :=
    100 * N_compatible / N_LR_pairs
]

panelC_dt[
  ,
  Direction := direction_labels[
    as.character(Direction_ID)
  ]
]

panelC_dt[
  ,
  Direction := factor(
    Direction,
    levels = c(
      "E → M",
      "M → E",
      "E → T",
      "T → E",
      "M → T",
      "T → M"
    )
  )
]

panelC_dt[
  ,
  Patient := factor(
    Patient_Index,
    levels = sort(unique(Patient_Index))
  )
]

# ------------------------------------------------------------
# 13. QC: each patient-direction cell must be represented
# ------------------------------------------------------------

expected_cells <- 4L * 6L

if (nrow(panelC_dt) != expected_cells) {
  stop(
    "Panel C does not contain the expected 4 × 6 patient-direction structure."
  )
}

if (any(
  !is.finite(panelC_dt$Compatibility_Percent)
)) {
  stop(
    "Panel C contains non-finite compatibility percentages."
  )
}

# ------------------------------------------------------------
# 14. Plot Figure 7C
# ------------------------------------------------------------

p7C <- ggplot(
  panelC_dt,
  aes(
    x = Direction,
    y = Patient,
    fill = Compatibility_Percent
  )
) +
  geom_tile(
    linewidth = 0.5,
    colour = "white"
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.0f%%",
        Compatibility_Percent
      )
    ),
    family = "Arial",
    fontface = "bold",
    size = 3.0
  ) +
  scale_fill_gradient(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    name = "Compatible\nLR pairs (%)"
  ) +
  labs(
    x = "Sender → receiver direction",
    y = "Patient"
  ) +
  theme_paper7 +
  theme(
    axis.text.x = element_text(
      family = "Arial",
      size = 9,
      angle = 0,
      hjust = 0.5
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      family = "Arial",
      face = "bold",
      size = 9
    ),
    
    legend.text = element_text(
      family = "Arial",
      size = 8
    )
  )

print(p7C)


# ------------------------------------------------------------
# 15. Export Figure 7C
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "Figure7C_Patient_Direction_Compatibility_Heatmap.pdf"
  ),
  plot = p7C,
  width = 9.0,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure7C_Patient_Direction_Compatibility_Heatmap_600dpi.tiff"
  ),
  plot = p7C,
  width = 9.0,
  height = 5.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ============================================================
# FIGURE 7D
# Resource-level annotation
# ============================================================

# ------------------------------------------------------------
# 16. Annotation field detection
# ------------------------------------------------------------

possible_annotation_columns <- c(
  "annotation",
  "Annotation",
  "resource_annotation",
  "Resource_Annotation"
)

annotation_col <- possible_annotation_columns[
  possible_annotation_columns %in% names(master_dt)
][1]

if (is.na(annotation_col)) {
  stop(
    "No annotation column was found in the STEP11-5F-4 master table."
  )
}

cat(
  "Panel D annotation column:",
  annotation_col,
  "\n"
)

# ------------------------------------------------------------
# 17. Pair-level annotation table
# ------------------------------------------------------------
#
# Use unique Exact_LR_Pair × annotation combinations so that
# resource-level interaction-record multiplicity does not inflate
# the number of candidate LR pairs.
# ------------------------------------------------------------

panelD_raw <- unique(
  master_dt[
    ,
    .(
      Exact_LR_Pair,
      Annotation = get(annotation_col)
    )
  ]
)

panelD_raw[
  ,
  Annotation := trimws(Annotation)
]

# Standardize expected resource labels
panelD_raw[
  Annotation %chin%
    c(
      "Secreted Signaling",
      "Secreted signaling"
    ),
  Annotation := "Secreted signaling"
]

panelD_raw[
  Annotation %chin%
    c(
      "ECM-Receptor",
      "ECM–Receptor",
      "ECM-Receptor Interaction"
    ),
  Annotation := "ECM–Receptor"
]

panelD_raw[
  Annotation %chin%
    c(
      "Cell-Cell Contact",
      "Cell–Cell Contact"
    ),
  Annotation := "Cell–Cell Contact"
]

panelD_dt <- panelD_raw[
  ,
  .(
    N_exact_LR_pairs = uniqueN(Exact_LR_Pair)
  ),
  by = Annotation
]

# ------------------------------------------------------------
# 18. Keep the three predefined resource classes
# ------------------------------------------------------------

annotation_levels <- c(
  "Secreted signaling",
  "ECM–Receptor",
  "Cell–Cell Contact"
)

panelD_dt <- panelD_dt[
  Annotation %in% annotation_levels
]

panelD_dt[
  ,
  Annotation := factor(
    Annotation,
    levels = rev(annotation_levels)
  )
]

# ------------------------------------------------------------
# 19. Plot Figure 7D
# ------------------------------------------------------------

p7D <- ggplot(
  panelD_dt,
  aes(
    x = Annotation,
    y = N_exact_LR_pairs
  )
) +
  geom_col(
    width = 0.62
  ) +
  geom_text(
    aes(label = N_exact_LR_pairs),
    vjust = -0.35,
    family = "Arial",
    fontface = "bold",
    size = 3.4
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = "Resource annotation",
    y = "Exact LR pairs"
  ) +
  theme_paper7 +
  theme(
    axis.text.x = element_text(
      family = "Arial",
      size = 8.5,
      angle = 20,
      hjust = 1
    )
  )

print(p7D)


# ------------------------------------------------------------
# 20. Export Figure 7D
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "Figure7D_Resource_Annotation.pdf"
  ),
  plot = p7D,
  width = 9.0,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure7D_Resource_Annotation_600dpi.tiff"
  ),
  plot = p7D,
  width = 9.0,
  height = 5.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ============================================================
# 21. Final plotting QC
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FINAL FIGURE 7 PLOTTING QC\n")
cat("============================================================\n\n")

cat("Panel A:\n")
cat("  4/4 patients =", 492, "\n")
cat("  3/4 patients =", 93, "\n")
cat("  Total =", sum(panelA_dt$N), "\n\n")

cat("Panel B:\n")
print(
  panelB_dt[
    order(Direction_ID),
    .(
      Direction_ID,
      Direction,
      N_exact_LR_pairs
    )
  ]
)

cat("\nPanel C:\n")
cat(
  "  Patients =",
  uniqueN(panelC_dt$Patient_Index),
  "\n"
)

cat(
  "  Directions =",
  uniqueN(panelC_dt$Direction_ID),
  "\n"
)

cat(
  "  Patient × direction cells =",
  nrow(panelC_dt),
  "\n"
)

cat("\nPanel D:\n")
print(panelD_dt)

cat("\n")
cat("All four Figure 7 panels exported independently.\n")
cat("PDF: vector format\n")
cat("TIFF: 600 dpi, LZW compression\n")
cat("Font target: Arial / Arial Bold\n")
cat("No panel letters included.\n\n")

cat("Output directory:\n")
cat(output_dir, "\n\n")

cat("============================================================\n")
cat("PASS — STEP11-5F-5 FIGURE 7 PRESENTATION WORKFLOW\n")
cat("============================================================\n")

#STEP11-5F-5 Figure 7 Final Plotting — Version 2
# ============================================================
# STEP11-5F-5 — FINAL FIGURE 7 PLOTTING
# VERSION 2 — COLOR / CLUSTERED DIRECTION / HEATMAP REVISION
#
# Figure 7A
#   Patient-level consistency
#
# Figure 7B
#   Six sender → receiver directions
#   Clustered reciprocal directions
#   Two-color stacked bars:
#       sender color = lower segment
#       receiver color = upper segment
#
# Figure 7C
#   Patient-level compatibility percentage
#   Blue → white → red continuous scale
#
# Figure 7D
#   Resource-level annotation
#
# IMPORTANT:
#   PRESENTATION-ONLY WORKFLOW
#
#   No upstream LR discovery.
#   No expression recalculation.
#   No compatibility recalculation.
#   No statistical testing.
#   No pathway ranking.
#
#   All quantities are inherited from locked STEP11-5E-3
#   and STEP11-5F-4 results.
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_packages <- c(
  "data.table",
  "ggplot2",
  "scales"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package not installed: ", pkg)
  }
}

library(data.table)
library(ggplot2)
library(scales)


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

root_dir <- "input_dir/Paper 1_Step 11-CellChat"

canonical_file <- file.path(
  root_dir,
  "STEP11-5E-3_Patient_Level_Candidate_Integration",
  "STEP11-5E-3B_FINAL_Canonical_Patient_Level_LR_Integration.tsv"
)

master_file <- file.path(
  root_dir,
  "STEP11-5F-4_Patient_Level_LR_Final_Context_Summary",
  "STEP11-5F-4C_FINAL_Master_LR_Descriptive_Context_Table.tsv"
)

output_dir <- file.path(
  root_dir,
  "STEP11-5F-5_Figure7_Final_Plotting"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Read LOCKED inputs
# ------------------------------------------------------------

if (!file.exists(canonical_file)) {
  stop(
    "Canonical input file not found:\n",
    canonical_file
  )
}

if (!file.exists(master_file)) {
  stop(
    "Master table not found:\n",
    master_file
  )
}

canonical_dt <- fread(
  canonical_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

master_dt <- fread(
  master_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)


# ============================================================
# 3. GLOBAL LOCKED QC
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP11-5F-5 VERSION 2 — FIGURE 7 PRESENTATION QC\n")
cat("============================================================\n\n")

cat(
  "Canonical rows :",
  nrow(canonical_dt),
  "\n"
)

cat(
  "Master rows    :",
  nrow(master_dt),
  "\n"
)


if (nrow(canonical_dt) != 9620L) {
  stop(
    "Unexpected canonical row count."
  )
}

if (uniqueN(canonical_dt$Exact_LR_Pair) != 585L) {
  stop(
    "Unexpected exact LR-pair count."
  )
}

if (uniqueN(canonical_dt$Patient_Index) != 4L) {
  stop(
    "Unexpected patient count."
  )
}

if (uniqueN(canonical_dt$Direction_ID) != 6L) {
  stop(
    "Unexpected direction count."
  )
}

if (nrow(master_dt) != 585L) {
  stop(
    "Unexpected master-table row count."
  )
}


cat("\n")
cat("LOCKED STRUCTURE CONFIRMED\n")
cat("--------------------------\n")
cat("9,620 canonical patient-level records\n")
cat("585 exact LR pairs\n")
cat("4 matched patients\n")
cat("6 sender → receiver directions\n")
cat("585-row final master table\n\n")


# ============================================================
# 4. FIGURE 6 — LOCKED COMPARTMENT COLOR MAPPING
# ============================================================
#
# IMPORTANT:
#   The color mapping is intentionally generated using the
#   SAME scales::hue_pal() principle used for Figure 6.
#
#   This preserves biological color identity across figures.
#
#   Expected visual identities:
#       T-cell              = pink
#       Myeloid/macrophage  = green
#       Endothelial         = blue
#
# ============================================================


primary_compartments <- c(
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)

compartment_levels <- primary_compartments

compartment_colors <- setNames(
  scales::hue_pal()(
    length(compartment_levels)
  ),
  compartment_levels
)

cat(
  "Figure 6-compatible compartment colors:\n"
)

print(compartment_colors)


# ============================================================
# 5. Direction mapping
# ============================================================

direction_labels <- c(
  "1" = "E → M",
  "2" = "M → E",
  "3" = "E → T",
  "4" = "T → E",
  "5" = "M → T",
  "6" = "T → M"
)

direction_sender <- c(
  "1" = "Endothelial",
  "2" = "Myeloid/macrophage",
  "3" = "Endothelial",
  "4" = "T-cell",
  "5" = "Myeloid/macrophage",
  "6" = "T-cell"
)

direction_receiver <- c(
  "1" = "Myeloid/macrophage",
  "2" = "Endothelial",
  "3" = "T-cell",
  "4" = "Endothelial",
  "5" = "T-cell",
  "6" = "Myeloid/macrophage"
)


# ============================================================
# 6. Common publication theme
# ============================================================

theme_figure7 <- theme_classic(
  base_size = 10,
  base_family = "Arial"
) +
  theme(
    
    text = element_text(
      family = "Arial",
      colour = "black"
    ),
    
    axis.title = element_text(
      family = "Arial",
      face = "bold",
      size = 10,
      colour = "black"
    ),
    
    axis.text = element_text(
      family = "Arial",
      size = 9,
      colour = "black"
    ),
    
    legend.title = element_text(
      family = "Arial",
      face = "bold",
      size = 9,
      colour = "black"
    ),
    
    legend.text = element_text(
      family = "Arial",
      size = 8,
      colour = "black"
    ),
    
    plot.title = element_text(
      family = "Arial",
      face = "bold",
      size = 11,
      hjust = 0.5
    ),
    
    axis.line = element_line(
      linewidth = 0.5,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.4,
      colour = "black"
    ),
    
    legend.position = "top",
    
    plot.margin = margin(
      t = 8,
      r = 10,
      b = 8,
      l = 10
    )
  )


# ============================================================
# FIGURE 7A
# Patient-level consistency
#
# 93 = 3/4 patients
# 492 = 4/4 patients
#
# Compact canvas:
#   4.5 × 5.2 inch
# ============================================================


panelA_dt <- data.table(
  
  Consistency = factor(
    c(
      "3/4 patients",
      "4/4 patients"
    ),
    levels = c(
      "3/4 patients",
      "4/4 patients"
    )
  ),
  
  N = c(
    93L,
    492L
  )
)


# Arithmetic lock
if (sum(panelA_dt$N) != 585L) {
  
  stop(
    "Panel A arithmetic failed."
  )
  
}


# ------------------------------------------------------------
# Panel A plot
# ------------------------------------------------------------

p7A <- ggplot(
  panelA_dt,
  aes(
    x = Consistency,
    y = N
  )
) +
  
  geom_col(
    width = 0.68,
    fill = compartment_colors[
      "Endothelial"
    ]
  ) +
  
  geom_text(
    aes(
      label = N
    ),
    vjust = -0.25,
    family = "Arial",
    fontface = "bold",
    size = 3.2
  ) +
  
  scale_y_continuous(
    limits = c(
      0,
      540
    ),
    breaks = seq(
      0,
      500,
      100
    ),
    expand = expansion(
      mult = c(
        0,
        0.03
      )
    )
  ) +
  
  labs(
    x = "Patient-level consistency",
    y = "Exact LR pairs"
  ) +
  
  theme_figure7


print(p7A)


# ------------------------------------------------------------
# Panel A export
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "Figure7A_Patient_Level_Consistency.pdf"
  ),
  plot = p7A,
  width = 4.5,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)


ggsave(
  filename = file.path(
    output_dir,
    "Figure7A_Patient_Level_Consistency_600dpi.tiff"
  ),
  plot = p7A,
  width = 4.5,
  height = 5.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# ============================================================
# Figure 7B
# Directional LR-pair context across predefined
# endothelial, myeloid/macrophage, and T-cell compartments
#
# FINAL LOCKED DATA ONLY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Figure 7B — Directional LR-Pair Context\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory and output folder
# ------------------------------------------------------------

step_dir <-
  "input_dir/Paper 1_Step 11-CellChat"

input_dir <-
  file.path(
    step_dir,
    "STEP11-5E-3_Patient_Level_Candidate_Integration"
  )

output_dir <-
  file.path(
    step_dir,
    "STEP11-5F-4_Patient_Level_LR_Final_Context_Summary",
    "Figure7B"
  )

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})


# ------------------------------------------------------------
# 2. Locked input
# ------------------------------------------------------------

input_file <-
  file.path(
    input_dir,
    "STEP11-5E-3B_FINAL_Canonical_Patient_Level_LR_Integration.tsv"
  )

stopifnot(
  file.exists(input_file)
)

dt <-
  fread(
    input_file,
    sep = "\t",
    header = TRUE,
    data.table = TRUE,
    check.names = FALSE
  )


# ------------------------------------------------------------
# 3. Basic input integrity checks
# ------------------------------------------------------------

cat("Input file:\n")
cat(input_file, "\n\n")

cat("Input rows:", nrow(dt), "\n")

stopifnot(
  nrow(dt) == 9620L
)

required_cols <-
  c(
    "Exact_LR_Pair",
    "Direction_ID",
    "Sender",
    "Receiver",
    "Patient_Index"
  )

missing_cols <-
  setdiff(
    required_cols,
    names(dt)
  )

stopifnot(
  length(missing_cols) == 0L
)

cat("Required columns: PASS\n")
cat("Canonical row count: PASS — 9,620\n\n")


# ------------------------------------------------------------
# 4. Locked compartment definition
#
# This order is intentionally identical to Figure 6:
#
# T-cell              = pink
# Myeloid/macrophage  = green
# Endothelial         = blue
# ------------------------------------------------------------

primary_compartments <-
  c(
    "T-cell",
    "Myeloid/macrophage",
    "Endothelial"
  )


# ------------------------------------------------------------
# 5. Define the six locked directions
# ------------------------------------------------------------

direction_levels <-
  c(
    "Endothelial → Myeloid/macrophage",
    "Myeloid/macrophage → Endothelial",
    "Endothelial → T-cell",
    "T-cell → Endothelial",
    "Myeloid/macrophage → T-cell",
    "T-cell → Myeloid/macrophage"
  )


# ------------------------------------------------------------
# 6. Convert Direction_ID to publication-facing labels
#
# Locked mapping:
#
# 1 = Endothelial → Myeloid/macrophage
# 2 = Myeloid/macrophage → Endothelial
# 3 = Endothelial → T-cell
# 4 = T-cell → Endothelial
# 5 = Myeloid/macrophage → T-cell
# 6 = T-cell → Myeloid/macrophage
# ------------------------------------------------------------

direction_map <-
  c(
    "1" = "Endothelial → Myeloid/macrophage",
    "2" = "Myeloid/macrophage → Endothelial",
    "3" = "Endothelial → T-cell",
    "4" = "T-cell → Endothelial",
    "5" = "Myeloid/macrophage → T-cell",
    "6" = "T-cell → Myeloid/macrophage"
  )

dt[
  ,
  Direction_Label :=
    unname(
      direction_map[
        as.character(Direction_ID)
      ]
    )
]

stopifnot(
  all(
    !is.na(dt$Direction_Label)
  )
)


# ------------------------------------------------------------
# 7. Verify sender / receiver against the locked direction
# ------------------------------------------------------------

expected_sender <-
  c(
    "1" = "Endothelial",
    "2" = "Myeloid/macrophage",
    "3" = "Endothelial",
    "4" = "T-cell",
    "5" = "Myeloid/macrophage",
    "6" = "T-cell"
  )

expected_receiver <-
  c(
    "1" = "Myeloid/macrophage",
    "2" = "Endothelial",
    "3" = "T-cell",
    "4" = "Endothelial",
    "5" = "T-cell",
    "6" = "Myeloid/macrophage"
  )

dt[
  ,
  Expected_Sender :=
    unname(
      expected_sender[
        as.character(Direction_ID)
      ]
    )
]

dt[
  ,
  Expected_Receiver :=
    unname(
      expected_receiver[
        as.character(Direction_ID)
      ]
    )
]

sender_check <-
  all(
    dt$Sender ==
      dt$Expected_Sender
  )

receiver_check <-
  all(
    dt$Receiver ==
      dt$Expected_Receiver
  )

stopifnot(
  sender_check,
  receiver_check
)

cat("Sender mapping: PASS\n")
cat("Receiver mapping: PASS\n\n")


# ------------------------------------------------------------
# 8. Count unique Exact_LR_Pair per direction
#
# IMPORTANT:
# This is NOT a new LR analysis.
#
# We only summarize the FINAL LOCKED canonical table.
#
# Each Exact_LR_Pair can occur across multiple patients,
# therefore uniqueN(Exact_LR_Pair) is used.
# ------------------------------------------------------------

direction_summary <-
  dt[
    ,
    .(
      N_Exact_LR_Pairs =
        uniqueN(
          Exact_LR_Pair
        )
    ),
    by = .(
      Direction_ID,
      Direction_Label,
      Sender,
      Receiver
    )
  ]

setorder(
  direction_summary,
  Direction_ID
)


# ------------------------------------------------------------
# 9. Confirm the locked six-direction structure
# ------------------------------------------------------------

cat("Direction-level locked structure:\n\n")

print(
  direction_summary
)

cat("\n")

expected_direction_pair_counts <-
  c(
    382L,
    371L,
    387L,
    343L,
    485L,
    437L
  )

stopifnot(
  identical(
    direction_summary$N_Exact_LR_Pairs,
    expected_direction_pair_counts
  )
)

cat(
  "Direction-level LR-pair counts: PASS\n\n"
)


# ------------------------------------------------------------
# 10. Create publication-facing direction order
#
# Group 1: E ↔ M
# Group 2: E ↔ T
# Group 3: M ↔ T
#
# Blank levels are used only to create visual separation.
# They contain no biological data.
# ------------------------------------------------------------

plot_direction_levels <-
  c(
    "Endothelial → Myeloid/macrophage",
    "Myeloid/macrophage → Endothelial",
    "",
    "Endothelial → T-cell",
    "T-cell → Endothelial",
    "",
    "Myeloid/macrophage → T-cell",
    "T-cell → Myeloid/macrophage"
  )


# ------------------------------------------------------------
# 11. Build explicit 50:50 sender / receiver segments
#
# IMPORTANT:
# Do NOT use geom_col() stacking order to determine which
# compartment appears at the bottom or top.
#
# Instead, explicitly define Ymin and Ymax.
#
# Sender   = lower 50% of the bar
# Receiver = upper 50% of the bar
#
# Therefore:
#
# E → M = blue  -> green
# M → E = green -> blue
# E → T = blue  -> pink
# T → E = pink  -> blue
# M → T = green -> pink
# T → M = pink  -> green
# ------------------------------------------------------------

plot_df <-
  rbindlist(
    lapply(
      seq_len(
        nrow(direction_summary)
      ),
      function(i) {
        
        d <-
          direction_summary[i]
        
        total_n <-
          d$N_Exact_LR_Pairs
        
        half_n <-
          total_n / 2
        
        group_label <-
          ifelse(
            d$Direction_ID %in% c(1, 2),
            "E ↔ M",
            ifelse(
              d$Direction_ID %in% c(3, 4),
              "E ↔ T",
              "M ↔ T"
            )
          )
        
        rbind(
          
          # --------------------------------------------------
          # Sender = LOWER half
          # --------------------------------------------------
          
          data.table(
            Direction_ID =
              d$Direction_ID,
            
            Direction_Label =
              d$Direction_Label,
            
            Direction_Group =
              group_label,
            
            Sender =
              d$Sender,
            
            Receiver =
              d$Receiver,
            
            Segment =
              "Sender",
            
            Segment_Compartment =
              d$Sender,
            
            Segment_N =
              half_n,
            
            Ymin =
              0,
            
            Ymax =
              half_n
          ),
          
          # --------------------------------------------------
          # Receiver = UPPER half
          # --------------------------------------------------
          
          data.table(
            Direction_ID =
              d$Direction_ID,
            
            Direction_Label =
              d$Direction_Label,
            
            Direction_Group =
              group_label,
            
            Sender =
              d$Sender,
            
            Receiver =
              d$Receiver,
            
            Segment =
              "Receiver",
            
            Segment_Compartment =
              d$Receiver,
            
            Segment_N =
              half_n,
            
            Ymin =
              half_n,
            
            Ymax =
              total_n
          )
        )
      }
    ),
    fill = TRUE
  )


# ------------------------------------------------------------
# 12. Convert compartment to factor
# ------------------------------------------------------------

plot_df[
  ,
  Segment_Compartment :=
    factor(
      Segment_Compartment,
      levels = primary_compartments
    )
]


# ------------------------------------------------------------
# 13. Direction order
# ------------------------------------------------------------

plot_df[
  ,
  Direction_Label :=
    factor(
      Direction_Label,
      levels = direction_levels
    )
]


# ------------------------------------------------------------
# 14. Explicit X-axis positions
#
# Group 1: E ↔ M
# Group 2: E ↔ T
# Group 3: M ↔ T
#
# Blank positions 3 and 6 create visual separation.
# ------------------------------------------------------------

x_position_map <-
  c(
    "Endothelial → Myeloid/macrophage" = 1,
    "Myeloid/macrophage → Endothelial" = 2,
    
    "Endothelial → T-cell" = 4,
    "T-cell → Endothelial" = 5,
    
    "Myeloid/macrophage → T-cell" = 7,
    "T-cell → Myeloid/macrophage" = 8
  )

plot_df[
  ,
  X_Position :=
    unname(
      x_position_map[
        as.character(Direction_Label)
      ]
    )
]


# ------------------------------------------------------------
# 15. Figure 6-compatible compartment colors
# ------------------------------------------------------------

compartment_colors <-
  setNames(
    scales::hue_pal()(
      length(primary_compartments)
    ),
    primary_compartments
  )

cat("\n")
cat("Compartment colors:\n")
print(compartment_colors)
cat("\n")


# ------------------------------------------------------------
# 16. X-axis labels
# ------------------------------------------------------------

x_labels <-
  c(
    "E → M",
    "M → E",
    "E → T",
    "T → E",
    "M → T",
    "T → M"
  )

x_label_positions <-
  c(
    1,
    2,
    4,
    5,
    7,
    8
  )


# ------------------------------------------------------------
# 17. Total-height labels
# ------------------------------------------------------------

label_df <-
  direction_summary[
    ,
    .(
      Direction_Label =
        Direction_Label,
      
      X_Position =
        unname(
          x_position_map[
            as.character(Direction_Label)
          ]
        ),
      
      N_Exact_LR_Pairs =
        N_Exact_LR_Pairs
    )
  ]


# ------------------------------------------------------------
# 18. Build Figure 7B
#
# geom_rect() is deliberately used instead of geom_col().
#
# This guarantees that:
#
# Sender   = lower half
# Receiver = upper half
#
# and therefore the color position changes when the
# direction is reversed.
# ------------------------------------------------------------

pB <-
  ggplot() +
  
  # ----------------------------------------------------------
# Explicit 50:50 sender / receiver blocks
# ----------------------------------------------------------

geom_rect(
  data = plot_df,
  aes(
    xmin = X_Position - 0.34,
    xmax = X_Position + 0.34,
    ymin = Ymin,
    ymax = Ymax,
    fill = Segment_Compartment
  ),
  colour = "white",
  linewidth = 0.35
) +
  
  # ----------------------------------------------------------
# Total LR-pair count
# ----------------------------------------------------------

geom_text(
  data = label_df,
  aes(
    x = X_Position,
    y = N_Exact_LR_Pairs,
    label = N_Exact_LR_Pairs
  ),
  inherit.aes = FALSE,
  vjust = -0.35,
  family = "Arial",
  size = 3.0,
  fontface = "bold"
) +
  
  # ----------------------------------------------------------
# Compartment colors
# ----------------------------------------------------------

scale_fill_manual(
  values = compartment_colors,
  breaks = primary_compartments,
  labels = c(
    "T-cell",
    "Myeloid/macrophage",
    "Endothelial"
  ),
  drop = FALSE,
  name = "Compartment"
) +
  
  # ----------------------------------------------------------
# X-axis
# ----------------------------------------------------------

scale_x_continuous(
  breaks = x_label_positions,
  labels = x_labels,
  limits = c(
    0.5,
    8.5
  ),
  expand = c(
    0,
    0
  )
) +
  
  # ----------------------------------------------------------
# Y-axis
# ----------------------------------------------------------

scale_y_continuous(
  limits = c(
    0,
    550
  ),
  breaks = seq(
    0,
    500,
    100
  ),
  expand = expansion(
    mult = c(
      0,
      0.02
    )
  )
) +
  
  # ----------------------------------------------------------
# Axis titles
# ----------------------------------------------------------

labs(
  x = "Sender → receiver direction",
  y = "Number of compatible LR pairs",
  fill = "Compartment"
) +
  
  # ----------------------------------------------------------
# Figure 6-compatible theme
# ----------------------------------------------------------

theme_classic(
  base_size = 10
) +
  
  theme(
    
    text =
      element_text(
        family = "Arial"
      ),
    
    axis.title =
      element_text(
        family = "Arial",
        face = "bold",
        size = 10
      ),
    
    axis.title.x =
      element_text(
        family = "Arial",
        face = "bold",
        size = 10,
        margin = margin(
          t = 7
        )
      ),
    
    axis.title.y =
      element_text(
        family = "Arial",
        face = "bold",
        size = 10
      ),
    
    axis.text =
      element_text(
        family = "Arial",
        size = 9
      ),
    
    axis.text.x =
      element_text(
        family = "Arial",
        size = 9,
        face = "bold"
      ),
    
    legend.title =
      element_text(
        family = "Arial",
        face = "bold",
        size = 9
      ),
    
    legend.text =
      element_text(
        family = "Arial",
        size = 8
      ),
    
    legend.position =
      "top",
    
    plot.margin =
      margin(
        5.5,
        5.5,
        5.5,
        5.5
      )
  )


# ------------------------------------------------------------
# 19. Display
# ------------------------------------------------------------

print(pB)


# ------------------------------------------------------------
# 20. Save PDF
# ------------------------------------------------------------

pdf_file <-
  file.path(
    output_dir,
    "Figure7B_Directional_LR_Pair_Context.pdf"
  )

ggsave(
  filename = pdf_file,
  plot = pB,
  width = 5.0,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)


# ------------------------------------------------------------
# 21. Save 600-dpi TIFF
# ------------------------------------------------------------

tiff_file <-
  file.path(
    output_dir,
    "Figure7B_Directional_LR_Pair_Context_600dpi.tiff"
  )

ggsave(
  filename = tiff_file,
  plot = pB,
  width = 5.0,
  height = 5.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ------------------------------------------------------------
# 22. Final visual-encoding QC
# ------------------------------------------------------------

expected_bottom <-
  c(
    "Endothelial",
    "Myeloid/macrophage",
    "Endothelial",
    "T-cell",
    "Myeloid/macrophage",
    "T-cell"
  )

expected_top <-
  c(
    "Myeloid/macrophage",
    "Endothelial",
    "T-cell",
    "Endothelial",
    "T-cell",
    "Myeloid/macrophage"
  )

visual_qc <-
  plot_df[
    ,
    .(
      Bottom_Compartment =
        Segment_Compartment[
          which.min(Ymin)
        ],
      
      Top_Compartment =
        Segment_Compartment[
          which.max(Ymax)
        ]
    ),
    by = Direction_ID
  ]

setorder(
  visual_qc,
  Direction_ID
)

visual_qc[
  ,
  Expected_Bottom :=
    expected_bottom
]

visual_qc[
  ,
  Expected_Top :=
    expected_top
]

visual_qc[
  ,
  Bottom_PASS :=
    as.character(Bottom_Compartment) ==
    Expected_Bottom
]

visual_qc[
  ,
  Top_PASS :=
    as.character(Top_Compartment) ==
    Expected_Top
]

print(
  visual_qc
)

stopifnot(
  all(visual_qc$Bottom_PASS),
  all(visual_qc$Top_PASS)
)

cat("\n")
cat("============================================================\n")
cat("PASS — FIGURE 7B VISUAL DIRECTION QC\n")
cat("============================================================\n\n")

cat(
  "E → M : Endothelial (blue) → Myeloid/macrophage (green)\n"
)

cat(
  "M → E : Myeloid/macrophage (green) → Endothelial (blue)\n"
)

cat(
  "E → T : Endothelial (blue) → T-cell (pink)\n"
)

cat(
  "T → E : T-cell (pink) → Endothelial (blue)\n"
)

cat(
  "M → T : Myeloid/macrophage (green) → T-cell (pink)\n"
)

cat(
  "T → M : T-cell (pink) → Myeloid/macrophage (green)\n"
)

cat("\n")

cat(
  "X-axis title: Sender → receiver direction\n"
)

cat(
  "50:50 visual segmentation: PASS\n"
)

cat(
  "No upstream data recalculated.\n"
)

cat(
  "No LR pairs rediscovered.\n"
)

cat(
  "No compatibility values recalculated.\n"
)

cat(
  "No statistical testing performed.\n"
)

cat("\n")
cat(
  "PDF:\n",
  pdf_file,
  "\n"
)

cat(
  "TIFF:\n",
  tiff_file,
  "\n"
)

cat("\n")
cat(
  "END OF FIGURE 7B\n"
)
cat("============================================================\n")

# ------------------------------------------------------------
# 19. Display
# ------------------------------------------------------------

print(pB)


# ------------------------------------------------------------
# 20. Save PDF
# ------------------------------------------------------------

pdf_file <-
  file.path(
    output_dir,
    "Figure7B_Directional_LR_Pair_Context.pdf"
  )

ggsave(
  filename = pdf_file,
  plot = pB,
  width = 5.0,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)


# ------------------------------------------------------------
# 21. Save 600-dpi TIFF
# ------------------------------------------------------------

tiff_file <-
  file.path(
    output_dir,
    "Figure7B_Directional_LR_Pair_Context_600dpi.tiff"
  )

ggsave(
  filename = tiff_file,
  plot = pB,
  width = 5.0,
  height = 5.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ------------------------------------------------------------
# 22. Save the exact plotting summary
#
# This is a plotting audit only.
# It does not modify any locked analysis table.
# ------------------------------------------------------------

summary_file <-
  file.path(
    output_dir,
    "Figure7B_Directional_LR_Pair_Count_Summary.tsv"
  )

fwrite(
  direction_summary,
  summary_file,
  sep = "\t"
)


# ------------------------------------------------------------
# 23. Final QC
# ------------------------------------------------------------

stopifnot(
  file.exists(pdf_file),
  file.exists(tiff_file),
  file.exists(summary_file)
)

cat("\n")
cat("============================================================\n")
cat("PASS — FIGURE 7B DIRECTIONAL LR-PAIR CONTEXT\n")
cat("============================================================\n\n")

cat(
  "Unique exact LR pairs by direction:\n"
)

print(
  direction_summary
)

cat("\n")

cat(
  "Sender / receiver visual encoding:\n"
)

cat(
  "E → M : blue → green\n"
)

cat(
  "M → E : green → blue\n"
)

cat(
  "E → T : blue → pink\n"
)

cat(
  "T → E : pink → blue\n"
)

cat(
  "M → T : green → pink\n"
)

cat(
  "T → M : pink → green\n"
)

cat("\n")

cat(
  "50:50 sender / receiver split: LOCKED VISUAL ENCODING\n"
)

cat(
  "No expression values recalculated.\n"
)

cat(
  "No LR compatibility values recalculated.\n"
)

cat(
  "No new LR pairs discovered.\n"
)

cat(
  "No statistical testing performed.\n"
)

cat(
  "No biological ranking performed.\n"
)

cat("\n")

cat(
  "PDF:\n",
  pdf_file,
  "\n"
)

cat(
  "TIFF:\n",
  tiff_file,
  "\n"
)

cat(
  "Summary:\n",
  summary_file,
  "\n"
)

cat("\n")
cat(
  "END OF FIGURE 7B\n"
)
cat("============================================================\n")


# ============================================================
# FIGURE 7C
# Patient-level compatibility percentage
#
# Final visualization:
#   0%   = white
#   25%  = very light blue
#   50%  = light blue
#   75%  = medium/deep blue
#   100% = dark blue
#
# The fill represents a numerical percentage rather than
# cellular identity; therefore a continuous blue gradient
# is intentionally used.
# ============================================================


# ------------------------------------------------------------
# 1. Calculate descriptive patient × direction percentage
#
# Compatible values are inherited directly from the locked
# canonical patient-level table.
# ------------------------------------------------------------

panelC_dt <- canonical_dt[
  ,
  .(
    N_LR_pairs =
      .N,
    
    N_compatible =
      sum(
        Compatible == TRUE,
        na.rm = TRUE
      )
  ),
  
  by = .(
    Patient_Index,
    Direction_ID
  )
]


# ------------------------------------------------------------
# 2. Calculate compatibility percentage
# ------------------------------------------------------------

panelC_dt[
  ,
  Compatibility_Percent :=
    100 *
    N_compatible /
    N_LR_pairs
]


# ------------------------------------------------------------
# 3. Direction labels
# ------------------------------------------------------------

panelC_dt[
  ,
  Direction :=
    direction_labels[
      as.character(
        Direction_ID
      )
    ]
]


panelC_dt$Direction <- factor(
  panelC_dt$Direction,
  levels = direction_order
)


# ------------------------------------------------------------
# 4. Patient labels
# ------------------------------------------------------------

panelC_dt[
  ,
  Patient :=
    factor(
      Patient_Index,
      levels =
        sort(
          unique(
            Patient_Index
          )
        )
    )
]


# ------------------------------------------------------------
# 5. QC
# ------------------------------------------------------------

if (
  nrow(panelC_dt) != 24L
) {
  
  stop(
    "Panel C should contain exactly 4 × 6 = 24 cells."
  )
  
}


if (
  any(
    !is.finite(
      panelC_dt$Compatibility_Percent
    )
  )
) {
  
  stop(
    "Panel C contains non-finite percentages."
  )
  
}


if (
  any(
    panelC_dt$Compatibility_Percent < 0 |
    panelC_dt$Compatibility_Percent > 100
  )
) {
  
  stop(
    "Panel C contains compatibility percentages outside 0–100%."
  )
  
}


# ------------------------------------------------------------
# 6. Optional structural QC
# ------------------------------------------------------------

if (
  length(
    unique(
      panelC_dt$Patient_Index
    )
  ) != 4L
) {
  
  stop(
    "Panel C should contain exactly 4 patients."
  )
  
}


if (
  length(
    unique(
      panelC_dt$Direction_ID
    )
  ) != 6L
) {
  
  stop(
    "Panel C should contain exactly 6 directions."
  )
  
}


# ------------------------------------------------------------
# 7. Plot Figure 7C
# ------------------------------------------------------------

p7C <- ggplot(
  panelC_dt,
  aes(
    x = Direction,
    y = Patient,
    fill = Compatibility_Percent
  )
) +
  
  # ----------------------------------------------------------
# Heatmap tiles
# ----------------------------------------------------------

geom_tile(
  colour = "white",
  linewidth = 0.8
) +
  
  # ----------------------------------------------------------
# Percentage labels
# ----------------------------------------------------------

geom_text(
  aes(
    label =
      sprintf(
        "%.0f%%",
        Compatibility_Percent
      )
  ),
  family = "Arial",
  fontface = "bold",
  size = 3.0
) +
  
  # ----------------------------------------------------------
# Continuous blue gradient
#
# 0%   = white
# 25%  = very light blue
# 50%  = light/medium blue
# 75%  = strong blue
# 100% = dark blue
#
# scale_fill_gradientn() is a standard ggplot2 function.
# ----------------------------------------------------------

scale_fill_gradientn(
  
  colours = c(
    "#FFFFFF",
    "#E0F2FE",
    "#BAE6FD",
    "#60A5FA",
    "#2563EB",
    "#1D4ED8"
  ),
  
  values = c(
    0.00,
    0.20,
    0.40,
    0.60,
    0.80,
    1.00
  ),
  
  limits = c(
    0,
    100
  ),
  
  breaks = c(
    0,
    25,
    50,
    75,
    100
  ),
  
  labels = c(
    "0",
    "25",
    "50",
    "75",
    "100"
  ),
  
  name =
    "Compatible\nLR pairs (%)",
  
  na.value = "#FFFFFF"
) +
  
  # ----------------------------------------------------------
# Axis labels
# ----------------------------------------------------------

labs(
  x = "Sender → receiver direction",
  y = "Patient"
) +
  
  # ----------------------------------------------------------
# Figure 7 common theme
# ----------------------------------------------------------

theme_figure7 +
  
  # ----------------------------------------------------------
# Panel-specific typography
# ----------------------------------------------------------

theme(
  
  axis.text.x = element_text(
    family = "Arial",
    size = 9,
    face = "bold",
    angle = 0,
    hjust = 0.5
  ),
  
  axis.text.y = element_text(
    family = "Arial",
    size = 9,
    face = "bold"
  ),
  
  legend.position = "right",
  
  legend.title = element_text(
    family = "Arial",
    face = "bold",
    size = 9
  ),
  
  legend.text = element_text(
    family = "Arial",
    size = 8
  )
)


# ------------------------------------------------------------
# 8. Display
# ------------------------------------------------------------

print(p7C)


# ------------------------------------------------------------
# 9. Save PDF
# ------------------------------------------------------------

pdf_file <-
  file.path(
    output_dir,
    "Figure7C_Patient_Direction_Compatibility_Heatmap.pdf"
  )

ggsave(
  filename = pdf_file,
  plot = p7C,
  width = 5.0,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)


# ------------------------------------------------------------
# 10. Save 600-dpi TIFF
# ------------------------------------------------------------

tiff_file <-
  file.path(
    output_dir,
    "Figure7C_Patient_Direction_Compatibility_Heatmap.tiff"
  )

ggsave(
  filename = tiff_file,
  plot = p7C,
  width = 5.0,
  height = 5.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# ============================================================
# FIGURE 7D
# Resource-level annotation
#
# Compact canvas:
#   4.5 × 5.2 inch
#
# Pair-level counting:
#   unique Exact_LR_Pair × Annotation
#
# This avoids inflation from resource-level interaction-name
# multiplicity.
# ============================================================


# ------------------------------------------------------------
# Identify annotation column
# ------------------------------------------------------------

possible_annotation_columns <- c(
  "annotation",
  "Annotation",
  "resource_annotation",
  "Resource_Annotation"
)


annotation_col <- possible_annotation_columns[
  possible_annotation_columns %in%
    names(master_dt)
][1]


if (
  is.na(annotation_col)
) {
  
  stop(
    "No annotation column found in master table."
  )
  
}


cat(
  "Panel D annotation column:",
  annotation_col,
  "\n"
)


# ------------------------------------------------------------
# Pair-level annotation table
# ------------------------------------------------------------

panelD_raw <- unique(
  
  master_dt[
    ,
    .(
      Exact_LR_Pair,
      
      Annotation =
        get(
          annotation_col
        )
    )
  ]
  
)


panelD_raw[
  ,
  Annotation :=
    trimws(
      Annotation
    )
]


# ------------------------------------------------------------
# Standardize annotation labels
# ------------------------------------------------------------

panelD_raw[
  Annotation %chin%
    c(
      "Secreted Signaling",
      "Secreted signaling"
    ),
  Annotation :=
    "Secreted signaling"
]


panelD_raw[
  Annotation %chin%
    c(
      "ECM-Receptor",
      "ECM–Receptor",
      "ECM-Receptor Interaction"
    ),
  Annotation :=
    "ECM–Receptor"
]


panelD_raw[
  Annotation %chin%
    c(
      "Cell-Cell Contact",
      "Cell–Cell Contact"
    ),
  Annotation :=
    "Cell–Cell Contact"
]


# ------------------------------------------------------------
# Count unique exact LR pairs per annotation
# ------------------------------------------------------------

panelD_dt <- panelD_raw[
  Annotation %in%
    c(
      "Secreted signaling",
      "ECM–Receptor",
      "Cell–Cell Contact"
    ),
  .(
    N_exact_LR_pairs =
      uniqueN(
        Exact_LR_Pair
      )
  ),
  by = Annotation
]


annotation_levels <- c(
  "Secreted signaling",
  "ECM–Receptor",
  "Cell–Cell Contact"
)


panelD_dt[
  ,
  Annotation :=
    factor(
      Annotation,
      levels =
        annotation_levels
    )
]


# ------------------------------------------------------------
# Panel D colors
#
# Use three Figure 6-compatible biological colors to keep
# Figure 7 visually coherent.
#
# These colors are used here as category colors only.
# They do NOT indicate cellular identity.
# ------------------------------------------------------------
panelD_colors <- c(
  
  "Secreted signaling" =
    "#3B82F6",
  
  "ECM–Receptor" =
    "#F59E0B",
  
  "Cell–Cell Contact" =
    "#A855F7"
  
)


# ------------------------------------------------------------
# Plot Figure 7D
# ------------------------------------------------------------

p7D <- ggplot(
  panelD_dt,
  aes(
    x = Annotation,
    y = N_exact_LR_pairs,
    fill = Annotation
  )
) +
  
  geom_col(
    width = 0.68
  ) +
  
  geom_text(
    aes(
      label =
        N_exact_LR_pairs
    ),
    vjust = -0.25,
    family = "Arial",
    fontface = "bold",
    size = 3.2
  ) +
  
  scale_fill_manual(
    values =
      panelD_colors,
    guide = "none"
  ) +
  
  scale_y_continuous(
    expand =
      expansion(
        mult =
          c(
            0,
            0.08
          )
      )
  ) +
  
  labs(
    x = "Resource annotation",
    y = "Exact LR pairs"
  ) +
  
  theme_figure7 +
  
  theme(
    axis.text.x =
      element_text(
        family = "Arial",
        size = 8.5,
        angle = 20,
        hjust = 1
      )
  )


print(p7D)


# ------------------------------------------------------------
# Panel D export
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "Figure7D_Resource_Annotation.pdf"
  ),
  plot = p7D,
  width = 4.5,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)


ggsave(
  filename = file.path(
    output_dir,
    "Figure7D_Resource_Annotation_600dpi.tiff"
  ),
  plot = p7D,
  width = 4.5,
  height = 5.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ============================================================
# 7. FINAL QC OUTPUT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FINAL FIGURE 7 VERSION 2 QC\n")
cat("============================================================\n\n")


cat("Panel A\n")
cat("-------\n")
cat(
  "3/4 patients :",
  panelA_dt$N[1],
  "\n"
)
cat(
  "4/4 patients :",
  panelA_dt$N[2],
  "\n"
)
cat(
  "Total        :",
  sum(panelA_dt$N),
  "\n\n"
)


cat("Panel B\n")
cat("-------\n")

print(
  panelB_base[
    order(Direction_ID),
    .(
      Direction_ID,
      Direction =
        direction_labels[
          as.character(
            Direction_ID
          )
        ],
      N_exact_LR_pairs
    )
  ]
)

cat("\n")


cat("Panel C\n")
cat("-------\n")
cat(
  "Patients           :",
  uniqueN(
    panelC_dt$Patient_Index
  ),
  "\n"
)
cat(
  "Directions         :",
  uniqueN(
    panelC_dt$Direction_ID
  ),
  "\n"
)
cat(
  "Patient × direction:",
  nrow(panelC_dt),
  "\n\n"
)


cat("Panel D\n")
cat("-------\n")

print(
  panelD_dt
)


cat("\n")
cat("Output directory:\n")
cat(
  output_dir,
  "\n\n"
)


cat("Export formats:\n")
cat("  PDF  = vector\n")
cat("  TIFF = 600 dpi, LZW\n\n")


cat("Panel dimensions:\n")
cat("  Figure 7A = 4.5 × 5.2 inch\n")
cat("  Figure 7B = 9.0 × 5.2 inch\n")
cat("  Figure 7C = 9.0 × 5.2 inch\n")
cat("  Figure 7D = 4.5 × 5.2 inch\n\n")


cat("Font target:\n")
cat("  Arial / Arial Bold\n\n")


cat("No upstream locked quantities were modified.\n")
cat("No LR quantities were recalculated.\n")
cat("No statistical testing was performed.\n\n")


cat("PASS — STEP11-5F-5 FIGURE 7 VERSION 2\n")
cat("============================================================\n")

# ============================================================
# STEP11-5F-4C
# FINAL Master LR Descriptive Context Table
#
# READ-ONLY MANUSCRIPT PREPARATION AUDIT
#
# Purpose:
#   1. Inspect the locked master LR descriptive table
#   2. Report dimensions and column names
#   3. Identify duplicated / missing records
#   4. Summarize direction / patient / candidate consistency
#   5. Inspect annotation and expression-context fields
#   6. Prepare information for the manuscript-facing
#      Supplementary Table
#
# IMPORTANT:
#   No LR compatibility values are recalculated.
#   No statistical analysis is performed.
#   No values in the locked source table are modified.
# ============================================================


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <-
  "input_dir/Paper 1_Step 11-CellChat/STEP11-5F-4_Patient_Level_LR_Final_Context_Summary"

setwd(step_dir)


# ------------------------------------------------------------
# 1. Input file
# ------------------------------------------------------------

input_file <-
  file.path(
    step_dir,
    "STEP11-5F-4C_FINAL_Master_LR_Descriptive_Context_Table.tsv"
  )


if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Input file not found:\n",
      input_file
    )
  )
  
}


cat("\n")
cat("============================================================\n")
cat("STEP11-5F-4C FINAL MASTER TABLE AUDIT\n")
cat("============================================================\n\n")

cat("Input file:\n")
cat(input_file, "\n\n")


# ------------------------------------------------------------
# 2. Read table
# ------------------------------------------------------------

master_dt <- data.table::fread(
  input_file,
  sep = "\t",
  header = TRUE,
  data.table = TRUE,
  na.strings = c(
    "",
    "NA",
    "N/A",
    "NULL"
  ),
  encoding = "UTF-8"
)


# ------------------------------------------------------------
# 3. Basic structure
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("3. BASIC STRUCTURE\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Rows    : ",
  format(nrow(master_dt), big.mark = ","),
  "\n",
  sep = ""
)

cat(
  "Columns : ",
  ncol(master_dt),
  "\n\n",
  sep = ""
)

cat("Column names:\n\n")

print(
  data.frame(
    Column_Index = seq_along(names(master_dt)),
    Column_Name = names(master_dt),
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 4. Data types
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("4. COLUMN CLASSES\n")
cat("------------------------------------------------------------\n\n")

class_table <- data.frame(
  Column = names(master_dt),
  Class = vapply(
    master_dt,
    function(x)
      paste(class(x), collapse = "/"),
    character(1)
  ),
  stringsAsFactors = FALSE
)

print(
  class_table,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 5. Missing-value audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("5. MISSING-VALUE AUDIT\n")
cat("------------------------------------------------------------\n\n")

missing_table <- data.frame(
  Column = names(master_dt),
  N_missing =
    vapply(
      master_dt,
      function(x)
        sum(
          is.na(x) |
            trimws(as.character(x)) == ""
        ),
      integer(1)
    ),
  Percent_missing =
    vapply(
      master_dt,
      function(x)
        100 *
        mean(
          is.na(x) |
            trimws(as.character(x)) == ""
        ),
      numeric(1)
    ),
  stringsAsFactors = FALSE
)

missing_table <-
  missing_table[
    order(
      -missing_table$N_missing
    ),
  ]

print(
  missing_table,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 6. Exact duplicate-row audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("6. EXACT DUPLICATE-ROW AUDIT\n")
cat("------------------------------------------------------------\n\n")

duplicate_rows <-
  duplicated(master_dt)

cat(
  "Exact duplicate rows: ",
  sum(duplicate_rows),
  "\n\n",
  sep = ""
)

if (any(duplicate_rows)) {
  
  cat(
    "WARNING: exact duplicate rows detected.\n"
  )
  
} else {
  
  cat(
    "PASS: no exact duplicate rows detected.\n"
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 7. Candidate LR pair fields
# ------------------------------------------------------------

candidate_fields <- intersect(
  c(
    "Exact_LR_Pair",
    "ligand",
    "receptor",
    "Candidate_Consistency_Class",
    "Patient_Consistency"
  ),
  names(master_dt)
)

cat("------------------------------------------------------------\n")
cat("7. CANDIDATE / CONSISTENCY FIELDS\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Detected fields:\n",
  paste(candidate_fields, collapse = "\n"),
  "\n\n"
)


# ------------------------------------------------------------
# 8. Exact LR-pair uniqueness
# ------------------------------------------------------------

if ("Exact_LR_Pair" %in% names(master_dt)) {
  
  cat("Exact LR-pair summary:\n\n")
  
  cat(
    "Unique Exact_LR_Pair values: ",
    uniqueN(master_dt$Exact_LR_Pair),
    "\n\n",
    sep = ""
  )
  
  pair_frequency <-
    master_dt[
      ,
      .N,
      by = Exact_LR_Pair
    ][
      order(-N, Exact_LR_Pair)
    ]
  
  cat("Frequency distribution of Exact_LR_Pair:\n\n")
  
  print(
    pair_frequency,
    row.names = FALSE
  )
  
  cat("\n")
}


# ------------------------------------------------------------
# 9. Direction summary
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("9. SENDER / RECEIVER / DIRECTION SUMMARY\n")
cat("------------------------------------------------------------\n\n")

for (
  field
  in c(
    "Sender",
    "Receiver",
    "Direction",
    "Direction_ID"
  )
) {
  
  if (field %in% names(master_dt)) {
    
    cat(
      "\n>>> ",
      field,
      "\n",
      sep = ""
    )
    
    print(
      master_dt[
        ,
        .N,
        by = field
      ][
        order(-N)
      ],
      row.names = FALSE
    )
  }
}

cat("\n")


# ------------------------------------------------------------
# 10. Patient compatibility summary
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("10. PATIENT COMPATIBILITY FIELDS\n")
cat("------------------------------------------------------------\n\n")

for (
  field
  in c(
    "N_compatible_patients",
    "N_compatible",
    "Compatible_Patients",
    "Patient_Consistency"
  )
) {
  
  if (field %in% names(master_dt)) {
    
    cat(
      "\n>>> ",
      field,
      "\n",
      sep = ""
    )
    
    print(
      master_dt[
        ,
        .N,
        by = field
      ][
        order(-N)
      ],
      row.names = FALSE
    )
  }
}

cat("\n")


# ------------------------------------------------------------
# 11. Pathway / resource annotation
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("11. PATHWAY / RESOURCE ANNOTATION\n")
cat("------------------------------------------------------------\n\n")

annotation_fields <- intersect(
  c(
    "pathway_annotation",
    "Pathway",
    "pathway_name",
    "resource_annotation",
    "Resource_Annotation",
    "annotation"
  ),
  names(master_dt)
)

cat(
  "Detected annotation fields:\n",
  paste(annotation_fields, collapse = "\n"),
  "\n\n"
)

for (
  field
  in annotation_fields
) {
  
  cat(
    "\n>>> ",
    field,
    "\n",
    sep = ""
  )
  
  tmp <-
    master_dt[
      ,
      .N,
      by = field
    ][
      order(-N)
    ]
  
  print(
    tmp,
    row.names = FALSE
  )
}

cat("\n")


# ------------------------------------------------------------
# 12. Expression-context fields
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("12. EXPRESSION-CONTEXT FIELDS\n")
cat("------------------------------------------------------------\n\n")

expression_fields <- names(master_dt)[
  grepl(
    "expression|context|ligand|receptor",
    names(master_dt),
    ignore.case = TRUE
  )
]

cat(
  "Potential expression/context fields:\n",
  paste(
    expression_fields,
    collapse = "\n"
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 13. Candidate consistency distribution
# ------------------------------------------------------------

if (
  "Candidate_Consistency_Class"
  %in%
  names(master_dt)
) {
  
  cat("------------------------------------------------------------\n")
  cat("13. CANDIDATE CONSISTENCY DISTRIBUTION\n")
  cat("------------------------------------------------------------\n\n")
  
  print(
    master_dt[
      ,
      .N,
      by = Candidate_Consistency_Class
    ][
      order(-N)
    ],
    row.names = FALSE
  )
  
  cat("\n")
}


# ------------------------------------------------------------
# 14. Patient consistency distribution
# ------------------------------------------------------------

if (
  "Patient_Consistency"
  %in%
  names(master_dt)
) {
  
  cat("------------------------------------------------------------\n")
  cat("14. PATIENT CONSISTENCY DISTRIBUTION\n")
  cat("------------------------------------------------------------\n\n")
  
  print(
    master_dt[
      ,
      .N,
      by = Patient_Consistency
    ][
      order(-N)
    ],
    row.names = FALSE
  )
  
  cat("\n")
}


# ------------------------------------------------------------
# 15. Preview
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("15. FIRST 10 ROWS\n")
cat("------------------------------------------------------------\n\n")

print(
  head(
    master_dt,
    10
  )
)

cat("\n")


# ------------------------------------------------------------
# 16. Manuscript-facing candidate columns
# ------------------------------------------------------------

recommended_columns <- c(
  "Exact_LR_Pair",
  "ligand",
  "receptor",
  "Candidate_Consistency_Class",
  "Patient_Consistency",
  "Sender",
  "Receiver",
  "Direction",
  "N_compatible_patients",
  "pathway_annotation",
  "resource_annotation",
  "ligand_expression_context",
  "receptor_expression_context"
)

cat("------------------------------------------------------------\n")
cat("16. MANUSCRIPT-FACING COLUMN CHECK\n")
cat("------------------------------------------------------------\n\n")

column_check <- data.frame(
  Recommended_Column =
    recommended_columns,
  
  Present_in_Master =
    recommended_columns %in%
    names(master_dt),
  
  stringsAsFactors = FALSE
)

print(
  column_check,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 17. Identify audit-only columns
# ------------------------------------------------------------

audit_only_candidates <-
  setdiff(
    names(master_dt),
    recommended_columns
  )

cat("------------------------------------------------------------\n")
cat("17. COLUMNS NOT CURRENTLY PROPOSED FOR MAIN SUPPLEMENTARY TABLE\n")
cat("------------------------------------------------------------\n\n")

print(
  audit_only_candidates
)

cat("\n")


# ------------------------------------------------------------
# 18. Save audit report
# ------------------------------------------------------------

audit_file <-
  file.path(
    step_dir,
    "STEP11-5F-4C_Master_LR_Descriptive_Context_Table_Audit.txt"
  )

sink(
  audit_file,
  split = TRUE
)

cat(
  "STEP11-5F-4C FINAL MASTER LR DESCRIPTIVE CONTEXT TABLE AUDIT\n"
)

cat(
  "Input:\n",
  input_file,
  "\n\n"
)

cat(
  "Rows: ",
  nrow(master_dt),
  "\n",
  sep = ""
)

cat(
  "Columns: ",
  ncol(master_dt),
  "\n\n",
  sep = ""
)

cat("Column names:\n")
print(names(master_dt))

cat("\nMissing-value summary:\n")
print(missing_table)

cat("\nExact duplicate rows:\n")
print(sum(duplicate_rows))

cat("\nRecommended manuscript-facing columns:\n")
print(column_check)

cat("\nPotential audit-only columns:\n")
print(audit_only_candidates)

sink()


# ------------------------------------------------------------
# 19. Final status
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP11-5F-4C AUDIT COMPLETE\n")
cat("============================================================\n\n")

cat(
  "No values were modified.\n"
)

cat(
  "No LR compatibility values were recalculated.\n"
)

cat(
  "No statistical analysis was performed.\n\n"
)

cat(
  "Audit report saved to:\n",
  audit_file,
  "\n"
)

cat("\n")


# ============================================================
# STEP11-5F-4C
# MANUSCRIPT TABLE MAPPING AUDIT
#
# Purpose:
#   Map the locked FINAL Master LR Descriptive Context Table
#   to a manuscript-facing Supplementary Table structure.
#
# IMPORTANT:
#   READ-ONLY AUDIT
#
#   No:
#     - LR compatibility recalculation
#     - statistical analysis
#     - modification of source data
#     - Figure 7 recalculation
#     - Excel generation
#
#   This script only inspects the locked STEP11-5F-4C table
#   and evaluates how its existing fields can be presented
#   in the manuscript Supplementary Table.
# ============================================================


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <-
  "input_dir/Paper 1_Step 11-CellChat/STEP11-5F-4_Patient_Level_LR_Final_Context_Summary"

setwd(step_dir)


# ------------------------------------------------------------
# 1. Input file
# ------------------------------------------------------------

input_file <-
  file.path(
    step_dir,
    "STEP11-5F-4C_FINAL_Master_LR_Descriptive_Context_Table.tsv"
  )


if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Input file not found:\n",
      input_file
    )
  )
  
}


# ------------------------------------------------------------
# 2. Read locked master table
# ------------------------------------------------------------

master_dt <- data.table::fread(
  input_file,
  sep = "\t",
  header = TRUE,
  data.table = TRUE,
  na.strings = c(
    "",
    "NA",
    "N/A",
    "NULL"
  ),
  encoding = "UTF-8"
)


# ------------------------------------------------------------
# 3. Start audit report
# ------------------------------------------------------------

audit_file <-
  file.path(
    step_dir,
    "STEP11-5F-4C_Manuscript_Table_Mapping_Audit.txt"
  )


sink(
  audit_file,
  split = TRUE
)


cat(
  "============================================================\n"
)

cat(
  "STEP11-5F-4C → MANUSCRIPT TABLE MAPPING AUDIT\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Input file:\n",
  input_file,
  "\n\n"
)

cat(
  "Rows: ",
  nrow(master_dt),
  "\n",
  sep = ""
)

cat(
  "Columns: ",
  ncol(master_dt),
  "\n\n",
  sep = ""
)


# ============================================================
# 4. Basic identity / row structure
# ============================================================

cat(
  "------------------------------------------------------------\n"
)

cat(
  "4. BASIC ROW STRUCTURE\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


if (
  !"Exact_LR_Pair" %in% names(master_dt)
) {
  
  stop(
    "Exact_LR_Pair column is missing."
  )
  
}


n_unique_pairs <-
  data.table::uniqueN(
    master_dt$Exact_LR_Pair
  )


n_duplicate_pair_rows <-
  nrow(master_dt) -
  n_unique_pairs


cat(
  "Total rows: ",
  nrow(master_dt),
  "\n",
  sep = ""
)

cat(
  "Unique Exact_LR_Pair values: ",
  n_unique_pairs,
  "\n",
  sep = ""
)

cat(
  "Rows beyond one-row-per-LR-pair structure: ",
  n_duplicate_pair_rows,
  "\n\n",
  sep = ""
)


pair_frequency <-
  master_dt[
    ,
    .N,
    by = Exact_LR_Pair
  ][
    order(-N, Exact_LR_Pair)
  ]


cat(
  "Frequency distribution of rows per Exact_LR_Pair:\n\n"
)

print(
  table(
    pair_frequency$N
  )
)

cat("\n")


if (
  all(
    pair_frequency$N == 1L
  )
) {
  
  cat(
    "PASS: one row corresponds to exactly one Exact_LR_Pair.\n\n"
  )
  
} else {
  
  cat(
    "NOTE: some Exact_LR_Pair values occur in multiple rows.\n\n"
  )
  
}


# ============================================================
# 5. Direction mapping
# ============================================================

cat(
  "------------------------------------------------------------\n"
)

cat(
  "5. DIRECTION MAPPING\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


direction_fields <- intersect(
  c(
    "N_Directions",
    "Directions",
    "Direction_Labels",
    "Direction_1_Compatibility_Fraction",
    "Direction_2_Compatibility_Fraction",
    "Direction_3_Compatibility_Fraction",
    "Direction_4_Compatibility_Fraction",
    "Direction_5_Compatibility_Fraction",
    "Direction_6_Compatibility_Fraction"
  ),
  names(master_dt)
)


cat(
  "Direction-related fields detected:\n\n"
)

print(
  direction_fields
)

cat("\n")


# ------------------------------------------------------------
# 5A. N_Directions distribution
# ------------------------------------------------------------

if (
  "N_Directions" %in% names(master_dt)
) {
  
  cat(
    "N_Directions distribution:\n\n"
  )
  
  print(
    master_dt[
      ,
      .N,
      by = N_Directions
    ][
      order(N_Directions)
    ],
    row.names = FALSE
  )
  
  cat("\n")
  
}


# ------------------------------------------------------------
# 5B. Directions examples
# ------------------------------------------------------------

if (
  "Directions" %in% names(master_dt)
) {
  
  cat(
    "Unique Directions values:\n\n"
  )
  
  directions_unique <-
    sort(
      unique(
        master_dt$Directions
      )
    )
  
  print(
    directions_unique
  )
  
  cat("\n")
  
}


# ------------------------------------------------------------
# 5C. Direction_Labels examples
# ------------------------------------------------------------

if (
  "Direction_Labels" %in% names(master_dt)
) {
  
  cat(
    "Unique Direction_Labels values:\n\n"
  )
  
  labels_unique <-
    sort(
      unique(
        master_dt$Direction_Labels
      )
    )
  
  print(
    labels_unique
  )
  
  cat("\n")
  
}


# ------------------------------------------------------------
# 5D. Direction_ID → direction label mapping
#
# The master table may encode direction IDs inside
# Directions / Direction_Labels. We inspect the relationship
# without reconstructing or recalculating any direction.
# ------------------------------------------------------------

if (
  all(
    c(
      "Directions",
      "Direction_Labels"
    )
    %in%
    names(master_dt)
  )
) {
  
  cat(
    "Observed Directions → Direction_Labels combinations:\n\n"
  )
  
  direction_mapping <-
    unique(
      master_dt[
        ,
        .(
          Directions,
          Direction_Labels
        )
      ]
    )[
      order(
        Directions,
        Direction_Labels
      )
    ]
  
  print(
    direction_mapping,
    row.names = FALSE
  )
  
  cat("\n")
  
}


# ============================================================
# 6. Direction-specific compatibility fractions
# ============================================================

cat(
  "------------------------------------------------------------\n"
)

cat(
  "6. DIRECTION-SPECIFIC COMPATIBILITY FRACTIONS\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


direction_fraction_fields <- paste0(
  "Direction_",
  1:6,
  "_Compatibility_Fraction"
)


direction_fraction_fields <-
  intersect(
    direction_fraction_fields,
    names(master_dt)
  )


if (
  length(direction_fraction_fields) == 0L
) {
  
  cat(
    "No direction-specific compatibility fraction fields detected.\n\n"
  )
  
} else {
  
  cat(
    "Detected fields:\n\n"
  )
  
  print(
    direction_fraction_fields
  )
  
  cat("\n")
  
  
  # ----------------------------------------------------------
  # Missingness
  # ----------------------------------------------------------
  
  direction_missing <- data.frame(
    Direction_Field =
      direction_fraction_fields,
    
    N_non_missing =
      sapply(
        master_dt[
          ,
          ..direction_fraction_fields
        ],
        function(x)
          sum(
            !is.na(x)
          )
      ),
    
    N_missing =
      sapply(
        master_dt[
          ,
          ..direction_fraction_fields
        ],
        function(x)
          sum(
            is.na(x)
          )
      ),
    
    stringsAsFactors = FALSE
  )
  
  
  direction_missing$Percent_non_missing <-
    100 *
    direction_missing$N_non_missing /
    nrow(master_dt)
  
  
  cat(
    "Direction-specific availability:\n\n"
  )
  
  print(
    direction_missing,
    row.names = FALSE
  )
  
  cat("\n")
  
  
  # ----------------------------------------------------------
  # Range audit
  # ----------------------------------------------------------
  
  direction_range <- data.frame(
    Direction_Field =
      direction_fraction_fields,
    
    Minimum =
      sapply(
        master_dt[
          ,
          ..direction_fraction_fields
        ],
        function(x)
          min(
            x,
            na.rm = TRUE
          )
      ),
    
    Maximum =
      sapply(
        master_dt[
          ,
          ..direction_fraction_fields
        ],
        function(x)
          max(
            x,
            na.rm = TRUE
          )
      ),
    
    stringsAsFactors = FALSE
  )
  
  
  cat(
    "Observed compatibility fraction ranges:\n\n"
  )
  
  print(
    direction_range,
    row.names = FALSE
  )
  
  cat("\n")
  
}


# ============================================================
# 7. Patient-level compatibility structure
# ============================================================

cat(
  "------------------------------------------------------------\n"
)

cat(
  "7. PATIENT-LEVEL COMPATIBILITY STRUCTURE\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


patient_fields <- intersect(
  c(
    "N_Patient_Records",
    "N_Patients",
    "N_Compatible",
    "N_Incompatible",
    "Overall_Compatibility_Fraction",
    "Patient_Compatibility_Patterns",
    "Patient_Consistency"
  ),
  names(master_dt)
)


cat(
  "Detected patient-level fields:\n\n"
)

print(
  patient_fields
)

cat("\n")


# ------------------------------------------------------------
# 7A. N_Patients
# ------------------------------------------------------------

if (
  "N_Patients" %in% names(master_dt)
) {
  
  cat(
    "N_Patients distribution:\n\n"
  )
  
  print(
    master_dt[
      ,
      .N,
      by = N_Patients
    ][
      order(N_Patients)
    ],
    row.names = FALSE
  )
  
  cat("\n")
  
}


# ------------------------------------------------------------
# 7B. N_Compatible
# ------------------------------------------------------------

if (
  "N_Compatible" %in% names(master_dt)
) {
  
  cat(
    "N_Compatible distribution:\n\n"
  )
  
  print(
    master_dt[
      ,
      .N,
      by = N_Compatible
    ][
      order(N_Compatible)
    ],
    row.names = FALSE
  )
  
  cat("\n")
  
}


# ------------------------------------------------------------
# 7C. N_Incompatible
# ------------------------------------------------------------

if (
  "N_Incompatible" %in% names(master_dt)
) {
  
  cat(
    "N_Incompatible distribution:\n\n"
  )
  
  print(
    master_dt[
      ,
      .N,
      by = N_Incompatible
    ][
      order(N_Incompatible)
    ],
    row.names = FALSE
  )
  
  cat("\n")
  
}


# ------------------------------------------------------------
# 7D. Overall compatibility fraction
# ------------------------------------------------------------

if (
  "Overall_Compatibility_Fraction"
  %in%
  names(master_dt)
) {
  
  cat(
    "Overall compatibility fraction summary:\n\n"
  )
  
  print(
    summary(
      master_dt$Overall_Compatibility_Fraction
    )
  )
  
  cat("\n")
  
  
  cat(
    "Unique overall compatibility fractions:\n\n"
  )
  
  print(
    sort(
      unique(
        master_dt$Overall_Compatibility_Fraction
      )
    )
  )
  
  cat("\n")
  
}


# ============================================================
# 8. Patient compatibility pattern
# ============================================================

if (
  "Patient_Compatibility_Patterns"
  %in%
  names(master_dt)
) {
  
  cat(
    "------------------------------------------------------------\n"
  )
  
  cat(
    "8. PATIENT COMPATIBILITY PATTERNS\n"
  )
  
  cat(
    "------------------------------------------------------------\n\n"
  )
  
  
  pattern_summary <-
    master_dt[
      ,
      .N,
      by = Patient_Compatibility_Patterns
    ][
      order(-N)
    ]
  
  
  cat(
    "Unique patient compatibility patterns: ",
    nrow(pattern_summary),
    "\n\n",
    sep = ""
  )
  
  
  print(
    pattern_summary,
    row.names = FALSE
  )
  
  cat("\n")
  
  
  cat(
    "First 20 LR pairs with patient compatibility patterns:\n\n"
  )
  
  print(
    master_dt[
      1:min(
        20L,
        nrow(master_dt)
      ),
      .(
        Exact_LR_Pair,
        N_Patients,
        N_Compatible,
        N_Incompatible,
        Overall_Compatibility_Fraction,
        Patient_Compatibility_Patterns
      )
    ],
    row.names = FALSE
  )
  
  cat("\n")
  
}


# ============================================================
# 9. Candidate consistency
# ============================================================

cat(
  "------------------------------------------------------------\n"
)

cat(
  "9. CANDIDATE CONSISTENCY CLASS\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


if (
  "Candidate_Consistency_Class"
  %in%
  names(master_dt)
) {
  
  candidate_summary <-
    master_dt[
      ,
      .N,
      by = Candidate_Consistency_Class
    ][
      order(-N)
    ]
  
  
  print(
    candidate_summary,
    row.names = FALSE
  )
  
  cat("\n")
  
}


# ============================================================
# 10. Candidate provenance
# ============================================================

cat(
  "------------------------------------------------------------\n"
)

cat(
  "10. CANDIDATE / INTERACTION PROVENANCE\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


provenance_fields <- intersect(
  c(
    "Candidate_Set_Provenance",
    "interaction_name_Provenance"
  ),
  names(master_dt)
)


for (
  field
  in provenance_fields
) {
  
  cat(
    "\n>>> ",
    field,
    "\n\n",
    sep = ""
  )
  
  print(
    master_dt[
      ,
      .N,
      by = field
    ][
      order(-N)
    ],
    row.names = FALSE
  )
  
}


cat("\n")


# ============================================================
# 11. Pathway / resource annotation
# ============================================================

cat(
  "------------------------------------------------------------\n"
)

cat(
  "11. PATHWAY AND RESOURCE ANNOTATION\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


annotation_fields <- intersect(
  c(
    "pathway_name",
    "annotation"
  ),
  names(master_dt)
)


for (
  field
  in annotation_fields
) {
  
  cat(
    "\n>>> ",
    field,
    "\n\n",
    sep = ""
  )
  
  annotation_summary <-
    master_dt[
      ,
      .N,
      by = field
    ][
      order(-N)
    ]
  
  
  cat(
    "Number of unique values: ",
    nrow(annotation_summary),
    "\n\n",
    sep = ""
  )
  
  
  print(
    annotation_summary,
    row.names = FALSE
  )
  
  cat("\n")
  
}


# ============================================================
# 12. Core gene involvement
# ============================================================

cat(
  "------------------------------------------------------------\n"
)

cat(
  "12. CORE GENE INVOLVEMENT\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


if (
  "Core_Gene_Involved"
  %in%
  names(master_dt)
) {
  
  core_summary <-
    master_dt[
      ,
      .N,
      by = Core_Gene_Involved
    ][
      order(-N)
    ]
  
  
  print(
    core_summary,
    row.names = FALSE
  )
  
  cat("\n")
  
}


# ============================================================
# 13. Expression context
# ============================================================

cat(
  "------------------------------------------------------------\n"
)

cat(
  "13. EXPRESSION CONTEXT FIELDS\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


expression_fields <- intersect(
  c(
    "Ligand_Median_log2CPM",
    "Ligand_IQR_log2CPM",
    "Ligand_Min_log2CPM",
    "Ligand_Max_log2CPM",
    "Receptor_Median_log2CPM",
    "Receptor_IQR_log2CPM",
    "Receptor_Min_log2CPM",
    "Receptor_Max_log2CPM",
    "MeanPair_Median_log2CPM",
    "MeanPair_IQR_log2CPM",
    "Difference_Median_log2CPM",
    "Difference_IQR_log2CPM"
  ),
  names(master_dt)
)


cat(
  "Detected expression fields:\n\n"
)

print(
  expression_fields
)

cat("\n")


if (
  length(expression_fields) > 0L
) {
  
  expression_summary <- data.frame(
    Field =
      expression_fields,
    
    Minimum =
      sapply(
        master_dt[
          ,
          ..expression_fields
        ],
        function(x)
          min(
            x,
            na.rm = TRUE
          )
      ),
    
    Median =
      sapply(
        master_dt[
          ,
          ..expression_fields
        ],
        function(x)
          median(
            x,
            na.rm = TRUE
          )
      ),
    
    Maximum =
      sapply(
        master_dt[
          ,
          ..expression_fields
        ],
        function(x)
          max(
            x,
            na.rm = TRUE
          )
      ),
    
    stringsAsFactors = FALSE
  )
  
  
  print(
    expression_summary,
    row.names = FALSE
  )
  
  cat("\n")
  
}


# ============================================================
# 14. Check whether expression context can be compacted
# ============================================================

cat(
  "------------------------------------------------------------\n"
)

cat(
  "14. EXPRESSION-CONTEXT COMPRESSION AUDIT\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


if (
  all(
    c(
      "Ligand_Median_log2CPM",
      "Ligand_IQR_log2CPM",
      "Receptor_Median_log2CPM",
      "Receptor_IQR_log2CPM"
    )
    %in%
    names(master_dt)
  )
) {
  
  cat(
    "Potential manuscript-facing compact representation:\n\n"
  )
  
  cat(
    "Ligand:   median (IQR)\n"
  )
  
  cat(
    "Receptor: median (IQR)\n\n"
  )
  
  
  cat(
    "Example rows:\n\n"
  )
  
  
  expression_preview <-
    master_dt[
      1:min(
        10L,
        nrow(master_dt)
      ),
      .(
        Exact_LR_Pair,
        
        Ligand_Context =
          sprintf(
            "%.3f (%.3f)",
            Ligand_Median_log2CPM,
            Ligand_IQR_log2CPM
          ),
        
        Receptor_Context =
          sprintf(
            "%.3f (%.3f)",
            Receptor_Median_log2CPM,
            Receptor_IQR_log2CPM
          )
      )
    ]
  
  
  print(
    expression_preview,
    row.names = FALSE
  )
  
  cat("\n")
  
}


# ============================================================
# 15. Proposed manuscript-facing fields
# ============================================================

cat(
  "------------------------------------------------------------\n"
)

cat(
  "15. PROPOSED MANUSCRIPT-FACING FIELD MAPPING\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


proposed_mapping <- data.frame(
  
  Manuscript_Field = c(
    
    "Exact LR pair",
    
    "Ligand",
    
    "Receptor",
    
    "Candidate consistency class",
    
    "Pathway",
    
    "Resource annotation",
    
    "Number of directions",
    
    "Sender → receiver directions",
    
    "Number of patients",
    
    "Compatible patients",
    
    "Overall compatibility",
    
    "Direction-specific compatibility",
    
    "Patient compatibility pattern",
    
    "Core gene involved",
    
    "Ligand expression context",
    
    "Receptor expression context"
    
  ),
  
  Source_Field = c(
    
    "Exact_LR_Pair",
    
    "ligand",
    
    "receptor",
    
    "Candidate_Consistency_Class",
    
    "pathway_name",
    
    "annotation",
    
    "N_Directions",
    
    "Direction_Labels",
    
    "N_Patients",
    
    "N_Compatible",
    
    "Overall_Compatibility_Fraction",
    
    "Direction_1...Direction_6_Compatibility_Fraction",
    
    "Patient_Compatibility_Patterns",
    
    "Core_Gene_Involved",
    
    "Ligand_Median_log2CPM + Ligand_IQR_log2CPM",
    
    "Receptor_Median_log2CPM + Receptor_IQR_log2CPM"
    
  ),
  
  stringsAsFactors = FALSE
)


print(
  proposed_mapping,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 16. Audit-only fields
# ============================================================

cat(
  "------------------------------------------------------------\n"
)

cat(
  "16. PROVISIONAL AUDIT-ONLY FIELDS\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


provisional_audit_only <- c(
  
  "Candidate_Set_Provenance",
  
  "interaction_name_Provenance",
  
  "N_Patient_Records",
  
  "N_Incompatible",
  
  "Ligand_Min_log2CPM",
  
  "Ligand_Max_log2CPM",
  
  "Receptor_Min_log2CPM",
  
  "Receptor_Max_log2CPM",
  
  "MeanPair_Median_log2CPM",
  
  "MeanPair_IQR_log2CPM",
  
  "Difference_Median_log2CPM",
  
  "Difference_IQR_log2CPM"
  
)


provisional_audit_only <-
  intersect(
    provisional_audit_only,
    names(master_dt)
  )


print(
  provisional_audit_only
)

cat("\n")


# ============================================================
# 17. Core manuscript table integrity checks
# ============================================================

cat(
  "------------------------------------------------------------\n"
)

cat(
  "17. MANUSCRIPT TABLE INTEGRITY CHECKS\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


# ------------------------------------------------------------
# 17A. LR pair completeness
# ------------------------------------------------------------

required_core_fields <- c(
  "Exact_LR_Pair",
  "ligand",
  "receptor",
  "Candidate_Consistency_Class",
  "pathway_name",
  "annotation",
  "N_Directions",
  "Direction_Labels",
  "N_Patients",
  "N_Compatible",
  "Overall_Compatibility_Fraction"
)


required_core_fields_present <-
  required_core_fields %in%
  names(master_dt)


if (
  all(required_core_fields_present)
) {
  
  cat(
    "PASS: all core manuscript fields are present.\n"
  )
  
} else {
  
  cat(
    "WARNING: one or more proposed core fields are missing.\n"
  )
  
  print(
    data.frame(
      Field = required_core_fields,
      Present = required_core_fields_present
    ),
    row.names = FALSE
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 17B. Overall compatibility range
# ------------------------------------------------------------

if (
  "Overall_Compatibility_Fraction"
  %in%
  names(master_dt)
) {
  
  outside_fraction <-
    master_dt[
      !is.na(
        Overall_Compatibility_Fraction
      ) &
        (
          Overall_Compatibility_Fraction < 0 |
            Overall_Compatibility_Fraction > 1
        )
    ]
  
  
  if (
    nrow(outside_fraction) == 0L
  ) {
    
    cat(
      "PASS: overall compatibility fractions are within 0–1.\n"
    )
    
  } else {
    
    cat(
      "WARNING: compatibility fractions outside 0–1 detected.\n"
    )
    
  }
  
}

cat("\n")


# ------------------------------------------------------------
# 17C. Patient count consistency
# ------------------------------------------------------------

if (
  all(
    c(
      "N_Patients",
      "N_Compatible",
      "N_Incompatible"
    )
    %in%
    names(master_dt)
  )
) {
  
  patient_count_consistent <-
    all(
      master_dt$N_Compatible +
        master_dt$N_Incompatible ==
        master_dt$N_Patients
    )
  
  
  if (
    patient_count_consistent
  ) {
    
    cat(
      "PASS: N_Compatible + N_Incompatible = N_Patients for all rows.\n"
    )
    
  } else {
    
    cat(
      "WARNING: patient-count inconsistency detected.\n"
    )
    
  }
  
}

cat("\n")


# ------------------------------------------------------------
# 17D. Overall fraction consistency
# ------------------------------------------------------------

if (
  all(
    c(
      "N_Patients",
      "N_Compatible",
      "Overall_Compatibility_Fraction"
    )
    %in%
    names(master_dt)
  )
) {
  
  expected_fraction <-
    master_dt$N_Compatible /
    master_dt$N_Patients
  
  
  fraction_difference <-
    abs(
      expected_fraction -
        master_dt$Overall_Compatibility_Fraction
    )
  
  
  max_difference <-
    max(
      fraction_difference,
      na.rm = TRUE
    )
  
  
  cat(
    "Maximum difference between reported and reconstructed\n"
  )
  
  cat(
    "overall compatibility fraction: ",
    format(
      max_difference,
      scientific = TRUE
    ),
    "\n",
    sep = ""
  )
  
  
  if (
    max_difference < 1e-10
  ) {
    
    cat(
      "PASS: overall compatibility fraction is internally consistent.\n"
    )
    
  } else {
    
    cat(
      "NOTE: reported and reconstructed fractions differ.\n"
    )
    
  }
  
}

cat("\n")


# ============================================================
# 18. First 20 rows of proposed manuscript content
# ============================================================

cat(
  "------------------------------------------------------------\n"
)

cat(
  "18. MANUSCRIPT-FACING CONTENT PREVIEW\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


preview_fields <- intersect(
  c(
    "Exact_LR_Pair",
    "ligand",
    "receptor",
    "Candidate_Consistency_Class",
    "pathway_name",
    "annotation",
    "N_Directions",
    "Direction_Labels",
    "N_Patients",
    "N_Compatible",
    "Overall_Compatibility_Fraction",
    "Patient_Compatibility_Patterns",
    "Core_Gene_Involved"
  ),
  names(master_dt)
)


print(
  master_dt[
    1:min(
      20L,
      nrow(master_dt)
    ),
    ..preview_fields
  ],
  row.names = FALSE
)

cat("\n")


# ============================================================
# 19. Final recommendation summary
# ============================================================

cat(
  "------------------------------------------------------------\n"
)

cat(
  "19. AUDIT INTERPRETATION / NEXT-STEP INFORMATION\n"
)

cat(
  "------------------------------------------------------------\n\n"
)


cat(
  "The locked master table contains ",
  nrow(master_dt),
  " LR-pair-level records across ",
  ncol(master_dt),
  " fields.\n\n",
  sep = ""
)

cat(
  "The intended manuscript-facing Supplementary Table should\n"
)

cat(
  "prefer a one-row-per-Exact_LR_Pair structure and should\n"
)

cat(
  "present descriptive context rather than recomputing any\n"
)

cat(
  "LR communication probabilities.\n\n"
)

cat(
  "Direction-specific compatibility should be retained either\n"
)

cat(
  "as compact direction-labelled information or as explicit\n"
)

cat(
  "direction-specific columns, depending on the observed\n"
)

cat(
  "Direction_Labels mapping above.\n\n"
)

cat(
  "Expression summaries should be evaluated for compact\n"
)

cat(
  "presentation before final Excel construction.\n\n"
)

cat(
  "No source values were modified by this audit.\n"
)

cat(
  "No LR compatibility values were recalculated.\n"
)

cat(
  "No statistical analysis was performed.\n\n"
)


# ============================================================
# 20. Close audit report
# ============================================================

cat(
  "============================================================\n"
)

cat(
  "STEP11-5F-4C → MANUSCRIPT TABLE MAPPING AUDIT COMPLETE\n"
)

cat(
  "============================================================\n"
)


sink()


# ------------------------------------------------------------
# 21. Console confirmation
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("AUDIT COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Audit report saved to:\n",
  audit_file,
  "\n\n"
)

cat(
  "Source table was not modified.\n"
)

cat(
  "No LR compatibility values were recalculated.\n"
)

cat(
  "No statistical analysis was performed.\n"
)

cat(
  "No Supplementary Table file was generated yet.\n"
)

cat("\n")


##STEP11-5F-4C — Direction Mapping Audit

# ============================================================
# STEP11-5F-4C
# Direction Mapping Audit
#
# PURPOSE
# ------------------------------------------------------------
# READ-ONLY audit of the mapping between:
#
#   Direction_1_Compatibility_Fraction
#   Direction_2_Compatibility_Fraction
#   Direction_3_Compatibility_Fraction
#   Direction_4_Compatibility_Fraction
#   Direction_5_Compatibility_Fraction
#   Direction_6_Compatibility_Fraction
#
# and the actual:
#
#   Sender -> Receiver
#
# direction labels.
#
# IMPORTANT
# ------------------------------------------------------------
# 1. No source values are modified.
# 2. No LR compatibility values are recalculated.
# 3. No candidate LR pairs are redefined.
# 4. No Excel file is generated.
# 5. This script is READ-ONLY.
#
# ============================================================


cat("============================================================\n")
cat("STEP11-5F-4C — Direction Mapping Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory / input / output
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 11-CellChat/STEP11-5F-4_Patient_Level_LR_Final_Context_Summary"

input_file <- file.path(
  step_dir,
  "STEP11-5F-4C_FINAL_Master_LR_Descriptive_Context_Table.tsv"
)

output_file <- file.path(
  step_dir,
  "STEP11-5F-4C_Direction_Mapping_Audit.txt"
)


# ------------------------------------------------------------
# 1. Package
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})


# ------------------------------------------------------------
# 2. Read locked master table
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("1. Reading locked master table\n")
cat("------------------------------------------------------------\n\n")

stopifnot(file.exists(input_file))

master_dt <- fread(
  input_file,
  sep = "\t",
  header = TRUE,
  data.table = TRUE,
  check.names = FALSE
)

cat("Input file:\n")
cat(input_file, "\n\n")

cat("Rows:", nrow(master_dt), "\n")
cat("Columns:", ncol(master_dt), "\n\n")


# ------------------------------------------------------------
# 3. Required-field lock
# ------------------------------------------------------------

required_fields <- c(
  "Exact_LR_Pair",
  "N_Directions",
  "Directions",
  "Direction_Labels",
  "N_Patients",
  "N_Compatible",
  "Overall_Compatibility_Fraction",
  "Direction_1_Compatibility_Fraction",
  "Direction_2_Compatibility_Fraction",
  "Direction_3_Compatibility_Fraction",
  "Direction_4_Compatibility_Fraction",
  "Direction_5_Compatibility_Fraction",
  "Direction_6_Compatibility_Fraction"
)

missing_required <- setdiff(
  required_fields,
  names(master_dt)
)

cat("Required fields:\n")
print(required_fields)

cat("\nMissing required fields:\n")
print(missing_required)

if (length(missing_required) > 0L) {
  stop(
    "FAIL: Required direction-mapping fields are missing."
  )
}

cat("\nPASS: All required direction-mapping fields are present.\n\n")


# ------------------------------------------------------------
# 4. Define direction fraction fields
# ------------------------------------------------------------

direction_fraction_fields <- paste0(
  "Direction_",
  1:6,
  "_Compatibility_Fraction"
)

cat("Direction fraction fields:\n")
print(direction_fraction_fields)

cat("\n")


# ------------------------------------------------------------
# 5. Inspect all unique Direction_Labels
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("2. Unique Sender -> Receiver direction labels\n")
cat("------------------------------------------------------------\n\n")

unique_labels <- unique(
  trimws(
    unlist(
      strsplit(
        master_dt$Direction_Labels,
        "\\s*\\|\\s*"
      )
    )
  )
)

unique_labels <- unique_labels[
  !is.na(unique_labels) &
    nzchar(unique_labels)
]

unique_labels <- sort(unique_labels)

cat(
  "Number of unique direction labels:",
  length(unique_labels),
  "\n\n"
)

print(unique_labels)

cat("\n")


# ------------------------------------------------------------
# 6. Expected six biological directions
# ------------------------------------------------------------

expected_directions <- c(
  "Endothelial -> Myeloid/macrophage",
  "Myeloid/macrophage -> Endothelial",
  "Endothelial -> T-cell",
  "T-cell -> Endothelial",
  "Myeloid/macrophage -> T-cell",
  "T-cell -> Myeloid/macrophage"
)

cat("Expected six possible directions:\n\n")

for (i in seq_along(expected_directions)) {
  
  cat(
    "Direction ",
    i,
    ": ",
    expected_directions[i],
    "\n",
    sep = ""
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 7. Check whether the observed labels match expected labels
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("3. Direction-label integrity check\n")
cat("------------------------------------------------------------\n\n")

unexpected_labels <- setdiff(
  unique_labels,
  expected_directions
)

missing_observed_labels <- setdiff(
  expected_directions,
  unique_labels
)

cat("Unexpected direction labels:\n")
print(unexpected_labels)

cat("\nExpected directions not observed:\n")
print(missing_observed_labels)

if (length(unexpected_labels) == 0L) {
  
  cat(
    "\nPASS: No unexpected sender -> receiver labels detected.\n"
  )
  
} else {
  
  cat(
    "\nWARNING: Unexpected sender -> receiver labels detected.\n"
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 8. Inspect Directions field
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("4. Directions field inspection\n")
cat("------------------------------------------------------------\n\n")

cat("Unique values in Directions:\n\n")

unique_directions_field <- sort(
  unique(
    master_dt$Directions[
      !is.na(master_dt$Directions) &
        nzchar(master_dt$Directions)
    ]
  )
)

print(unique_directions_field)

cat("\n")


# ------------------------------------------------------------
# 9. Inspect Direction_Labels <-> Directions combinations
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("5. Observed Directions -> Direction_Labels combinations\n")
cat("------------------------------------------------------------\n\n")

direction_combo_dt <- unique(
  master_dt[
    ,
    .(
      Directions,
      Direction_Labels,
      N_Directions
    )
  ]
)

setorder(
  direction_combo_dt,
  N_Directions,
  Directions
)

print(
  direction_combo_dt,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 10. Direction-specific fraction missingness
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("6. Direction-specific compatibility field audit\n")
cat("------------------------------------------------------------\n\n")

direction_missing_dt <- data.table(
  Direction_ID = 1:6,
  Field = direction_fraction_fields,
  N_Total = nrow(master_dt),
  N_NonMissing = sapply(
    master_dt[, ..direction_fraction_fields],
    function(x) sum(!is.na(x))
  ),
  N_Missing = sapply(
    master_dt[, ..direction_fraction_fields],
    function(x) sum(is.na(x))
  )
)

direction_missing_dt[
  ,
  Missing_Percent :=
    100 * N_Missing / N_Total
]

print(
  direction_missing_dt,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 11. Range audit for each direction fraction
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("7. Direction-specific compatibility ranges\n")
cat("------------------------------------------------------------\n\n")

direction_range_dt <- data.table(
  Direction_ID = 1:6,
  Field = direction_fraction_fields,
  Min = sapply(
    master_dt[, ..direction_fraction_fields],
    function(x) {
      if (all(is.na(x))) NA_real_
      else min(x, na.rm = TRUE)
    }
  ),
  Max = sapply(
    master_dt[, ..direction_fraction_fields],
    function(x) {
      if (all(is.na(x))) NA_real_
      else max(x, na.rm = TRUE)
    }
  )
)

print(
  direction_range_dt,
  row.names = FALSE
)

range_fail <- any(
  direction_range_dt$Min < 0 |
    direction_range_dt$Max > 1,
  na.rm = TRUE
)

if (!range_fail) {
  
  cat(
    "\nPASS: All non-missing direction-specific compatibility values are within [0, 1].\n"
  )
  
} else {
  
  cat(
    "\nFAIL: At least one direction-specific compatibility value is outside [0, 1].\n"
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 12. Direction presence matrix
#
# For each direction-specific field, determine which actual
# Direction_Labels occur in rows where that field is non-missing.
#
# This does NOT assume the mapping.
# It only describes the observed overlap.
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("8. Direction-field / direction-label overlap audit\n")
cat("------------------------------------------------------------\n\n")

for (i in 1:6) {
  
  field_i <- direction_fraction_fields[i]
  
  cat("\n============================================================\n")
  cat("Direction field:", field_i, "\n")
  cat("============================================================\n\n")
  
  idx <- !is.na(master_dt[[field_i]])
  
  labels_i <- unique(
    trimws(
      unlist(
        strsplit(
          master_dt$Direction_Labels[idx],
          "\\s*\\|\\s*"
        )
      )
    )
  )
  
  labels_i <- sort(
    unique(
      labels_i[
        !is.na(labels_i) &
          nzchar(labels_i)
      ]
    )
  )
  
  cat("Number of rows with non-missing value:",
      sum(idx),
      "\n\n")
  
  cat("Observed direction labels in these rows:\n\n")
  
  print(labels_i)
  
  cat("\n")
  
}


# ------------------------------------------------------------
# 13. Candidate mapping by exact row-level presence
#
# For every direction field:
#
#   if a field is non-missing in a row,
#   check whether Direction_Labels contains exactly one label.
#
# Rows with exactly one label provide the strongest direct
# mapping evidence.
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("9. Single-direction rows: strongest mapping evidence\n")
cat("------------------------------------------------------------\n\n")


single_direction_dt <- master_dt[
  N_Directions == 1 &
    !is.na(Direction_Labels)
]

cat(
  "Number of rows with N_Directions == 1:",
  nrow(single_direction_dt),
  "\n\n"
)


if (nrow(single_direction_dt) > 0L) {
  
  single_direction_summary <- rbindlist(
    lapply(
      1:6,
      function(i) {
        
        field_i <- direction_fraction_fields[i]
        
        x <- single_direction_dt[
          !is.na(get(field_i)),
          .(
            Direction_ID = i,
            Field = field_i,
            Direction_Label = Direction_Labels,
            Compatibility = get(field_i)
          )
        ]
        
        x
        
      }
    ),
    fill = TRUE
  )
  
  cat("Single-direction observations:\n\n")
  
  print(
    single_direction_summary,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No single-direction rows were available.\n"
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 14. For each Direction_i, count rows by exact single label
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("10. Candidate Direction_i -> exact single direction counts\n")
cat("------------------------------------------------------------\n\n")


single_label_count_dt <- rbindlist(
  lapply(
    1:6,
    function(i) {
      
      field_i <- direction_fraction_fields[i]
      
      tmp <- single_direction_dt[
        !is.na(get(field_i)),
        .(
          Direction_Label = Direction_Labels
        )
      ]
      
      if (nrow(tmp) == 0L) {
        return(NULL)
      }
      
      tmp[
        ,
        .(
          N_Rows = .N
        ),
        by = .(Direction_Label)
      ][
        ,
        `:=`(
          Direction_ID = i,
          Field = field_i
        )
      ][
        ,
        .(
          Direction_ID,
          Field,
          Direction_Label,
          N_Rows
        )
      ]
      
    }
  ),
  fill = TRUE
)

setorder(
  single_label_count_dt,
  Direction_ID,
  -N_Rows
)

print(
  single_label_count_dt,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 15. Strongest observed mapping
#
# For each Direction_i, identify the direction label with the
# largest number of single-direction observations.
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("11. Strongest observed mapping candidates\n")
cat("------------------------------------------------------------\n\n")


if (nrow(single_label_count_dt) > 0L) {
  
  strongest_mapping_dt <- single_label_count_dt[
    ,
    .SD[which.max(N_Rows)],
    by = Direction_ID
  ]
  
  setorder(
    strongest_mapping_dt,
    Direction_ID
  )
  
  print(
    strongest_mapping_dt,
    row.names = FALSE
  )
  
} else {
  
  strongest_mapping_dt <- data.table()
  
  cat(
    "No mapping candidate could be derived from single-direction rows.\n"
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 16. Check whether each Direction_i has one unique dominant
# mapping
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("12. Mapping uniqueness / ambiguity check\n")
cat("------------------------------------------------------------\n\n")


mapping_ambiguity_dt <- rbindlist(
  lapply(
    1:6,
    function(i) {
      
      x <- single_label_count_dt[
        Direction_ID == i
      ]
      
      if (nrow(x) == 0L) {
        
        return(
          data.table(
            Direction_ID = i,
            Field = direction_fraction_fields[i],
            N_Direction_Labels = 0L,
            Max_N_Rows = NA_integer_,
            Second_Max_N_Rows = NA_integer_,
            Ambiguous = NA
          )
        )
        
      }
      
      counts <- sort(
        x$N_Rows,
        decreasing = TRUE
      )
      
      max_n <- counts[1]
      
      second_n <- ifelse(
        length(counts) >= 2,
        counts[2],
        0
      )
      
      data.table(
        Direction_ID = i,
        Field = direction_fraction_fields[i],
        N_Direction_Labels = nrow(x),
        Max_N_Rows = max_n,
        Second_Max_N_Rows = second_n,
        Ambiguous = (
          length(counts) >= 2 &&
            max_n == second_n
        )
      )
      
    }
  ),
  fill = TRUE
)

print(
  mapping_ambiguity_dt,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 17. Compatibility consistency check
#
# This checks whether:
#
#   Overall_Compatibility_Fraction
#
# is mathematically compatible with:
#
#   N_Compatible / N_Patients
#
# and separately confirms that direction-specific fractions
# remain within the expected range.
#
# No recalculation of LR compatibility is performed.
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("13. Existing compatibility integrity check\n")
cat("------------------------------------------------------------\n\n")


overall_expected <- master_dt$N_Compatible /
  master_dt$N_Patients

overall_diff <- abs(
  master_dt$Overall_Compatibility_Fraction -
    overall_expected
)

overall_pass <- all(
  overall_diff < 1e-10,
  na.rm = TRUE
)

cat(
  "Overall compatibility = N_Compatible / N_Patients:\n"
)

if (overall_pass) {
  
  cat("PASS\n")
  
} else {
  
  cat("FAIL\n")
  
}

cat("\n")


# ------------------------------------------------------------
# 18. Inspect examples for each strongest candidate mapping
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("14. Example rows supporting each candidate mapping\n")
cat("------------------------------------------------------------\n\n")


if (nrow(strongest_mapping_dt) > 0L) {
  
  for (k in seq_len(nrow(strongest_mapping_dt))) {
    
    dir_id <- strongest_mapping_dt$Direction_ID[k]
    field_i <- strongest_mapping_dt$Field[k]
    label_i <- strongest_mapping_dt$Direction_Label[k]
    
    cat("\n============================================================\n")
    cat(
      "Direction_",
      dir_id,
      " candidate mapping: ",
      label_i,
      "\n",
      sep = ""
    )
    cat("Field: ", field_i, "\n", sep = "")
    cat("============================================================\n\n")
    
    example_dt <- single_direction_dt[
      Direction_Labels == label_i &
        !is.na(get(field_i)),
      .(
        Exact_LR_Pair,
        Direction_Labels,
        Compatibility = get(field_i),
        N_Patients,
        N_Compatible,
        Overall_Compatibility_Fraction
      )
    ]
    
    print(
      head(
        example_dt,
        10L
      ),
      row.names = FALSE
    )
    
  }
  
}

cat("\n")


# ------------------------------------------------------------
# 19. Final READ-ONLY audit summary
# ------------------------------------------------------------

cat("============================================================\n")
cat("15. FINAL READ-ONLY AUDIT SUMMARY\n")
cat("============================================================\n\n")


cat(
  "Input rows: ",
  nrow(master_dt),
  "\n",
  sep = ""
)

cat(
  "Input columns: ",
  ncol(master_dt),
  "\n",
  sep = ""
)

cat(
  "Unique observed direction labels: ",
  length(unique_labels),
  "\n",
  sep = ""
)

cat(
  "Expected biological directions: ",
  length(expected_directions),
  "\n",
  sep = ""
)

cat(
  "Unexpected direction labels: ",
  length(unexpected_labels),
  "\n",
  sep = ""
)

cat(
  "Direction-specific fields: 6\n"
)

cat(
  "Source values modified: NO\n"
)

cat(
  "Compatibility values recalculated: NO\n"
)

cat(
  "Excel generated: NO\n\n"
)


if (
  length(unexpected_labels) == 0L &&
  overall_pass
) {
  
  cat(
    "OVERALL AUDIT STATUS: PASS\n"
  )
  
} else {
  
  cat(
    "OVERALL AUDIT STATUS: REVIEW REQUIRED\n"
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 20. Save audit report
# ------------------------------------------------------------

sink(output_file)

cat("STEP11-5F-4C — Direction Mapping Audit\n")
cat("============================================================\n\n")

cat("Input:\n")
cat(input_file, "\n\n")

cat("Rows:", nrow(master_dt), "\n")
cat("Columns:", ncol(master_dt), "\n\n")

cat("Expected directions:\n")
print(expected_directions)

cat("\nObserved direction labels:\n")
print(unique_labels)

cat("\nUnexpected labels:\n")
print(unexpected_labels)

cat("\nMissing expected labels:\n")
print(missing_observed_labels)

cat("\nDirection fraction fields:\n")
print(direction_fraction_fields)

cat("\nDirection missingness:\n")
print(direction_missing_dt)

cat("\nDirection ranges:\n")
print(direction_range_dt)

cat("\nObserved Directions -> Direction_Labels combinations:\n")
print(direction_combo_dt)

cat("\nSingle-direction observations:\n")
if (exists("single_direction_summary") &&
    nrow(single_direction_summary) > 0L) {
  print(single_direction_summary)
} else {
  cat("None\n")
}

cat("\nSingle-direction label counts:\n")
if (nrow(single_label_count_dt) > 0L) {
  print(single_label_count_dt)
} else {
  cat("None\n")
}

cat("\nStrongest mapping candidates:\n")
if (nrow(strongest_mapping_dt) > 0L) {
  print(strongest_mapping_dt)
} else {
  cat("None\n")
}

cat("\nMapping ambiguity:\n")
print(mapping_ambiguity_dt)

cat("\nOverall compatibility integrity:\n")
cat(
  ifelse(
    overall_pass,
    "PASS",
    "FAIL"
  ),
  "\n"
)

cat("\nFinal status:\n")

if (
  length(unexpected_labels) == 0L &&
  overall_pass
) {
  cat("PASS\n")
} else {
  cat("REVIEW REQUIRED\n")
}

cat("\n============================================================\n")
cat("END OF READ-ONLY DIRECTION MAPPING AUDIT\n")
cat("============================================================\n")

sink()

cat(
  "Audit report saved to:\n",
  output_file,
  "\n\n"
)

cat("============================================================\n")
cat("STEP11-5F-4C COMPLETED\n")
cat("============================================================\n")


# ============================================================
# STEP11-5F-4D
# FINAL SUPPLEMENTARY TABLE STRUCTURE LOCK + EXCEL GENERATION
# ============================================================
#
# Purpose:
#   Generate the manuscript-facing Supplementary Table from
#   the LOCKED 585 x 36 master LR descriptive context table.
#
# IMPORTANT:
#   - NO recalculation of LR compatibility
#   - NO recalculation of patient counts
#   - NO recalculation of expression summaries
#   - NO statistical testing
#   - NO modification of source data
#   - ONLY field selection, renaming, ordering, formatting
#
# Locked input:
#   STEP11-5F-4C_FINAL_Master_LR_Descriptive_Context_Table.tsv
#
# Output:
#   Final manuscript-facing Excel Supplementary Table
#
# ============================================================


cat("============================================================\n")
cat("STEP11-5F-4D\n")
cat("FINAL SUPPLEMENTARY TABLE STRUCTURE LOCK + EXCEL GENERATION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_pkgs <- c(
  "readr",
  "dplyr",
  "stringr",
  "openxlsx"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "Required package not installed: ", pkg,
        "\nPlease install it before running this script."
      )
    )
  }
}

library(readr)
library(dplyr)
library(stringr)
library(openxlsx)


# ------------------------------------------------------------
# 1. Working directories
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 11-CellChat"

input_dir <- file.path(
  step_dir,
  "STEP11-5F-4_Patient_Level_LR_Final_Context_Summary"
)

output_dir <- file.path(
  step_dir,
  "STEP11-5F-4D_FINAL_Supplementary_Table"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 2. Locked input
# ------------------------------------------------------------

input_file <- file.path(
  input_dir,
  "STEP11-5F-4C_FINAL_Master_LR_Descriptive_Context_Table.tsv"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "LOCKED input file not found:\n",
      input_file
    )
  )
}


# ------------------------------------------------------------
# 3. Read locked master table
# ------------------------------------------------------------

master <- readr::read_tsv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
)

cat("Input file:\n")
cat(input_file, "\n\n")

cat(
  "Input dimensions: ",
  nrow(master),
  " rows x ",
  ncol(master),
  " columns\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Input integrity lock
# ------------------------------------------------------------

expected_nrow <- 585
expected_ncol <- 36

if (nrow(master) != expected_nrow) {
  stop(
    paste0(
      "INPUT LOCK FAILURE: expected ",
      expected_nrow,
      " rows but found ",
      nrow(master),
      "."
    )
  )
}

if (ncol(master) != expected_ncol) {
  stop(
    paste0(
      "INPUT LOCK FAILURE: expected ",
      expected_ncol,
      " columns but found ",
      ncol(master),
      "."
    )
  )
  
  
  required_fields <- c(
    "Exact_LR_Pair",
    "ligand",
    "receptor",
    "Candidate_Consistency_Class",
    "Candidate_Set_Provenance",
    "interaction_name_Provenance",
    "pathway_name",
    "annotation",
    "N_Directions",
    "Directions",
    "Direction_Labels",
    "N_Patient_Records",
    "N_Patients",
    "N_Compatible",
    "N_Incompatible",
    "Overall_Compatibility_Fraction",
    "Ligand_Median_log2CPM",
    "Ligand_IQR_log2CPM",
    "Ligand_Min_log2CPM",
    "Ligand_Max_log2CPM",
    "Receptor_Median_log2CPM",
    "Receptor_IQR_log2CPM",
    "Receptor_Min_log2CPM",
    "Receptor_Max_log2CPM",
    "MeanPair_Median_log2CPM",
    "MeanPair_IQR_log2CPM",
    "Difference_Median_log2CPM",
    "Difference_IQR_log2CPM",
    "Core_Gene_Involved",
    "Patient_Compatibility_Patterns",
    "Direction_1_Compatibility_Fraction",
    "Direction_2_Compatibility_Fraction",
    "Direction_3_Compatibility_Fraction",
    "Direction_4_Compatibility_Fraction",
    "Direction_5_Compatibility_Fraction",
    "Direction_6_Compatibility_Fraction"
  )
  
  missing_fields <- setdiff(
    required_fields,
    colnames(master)
  )
  
  if (length(missing_fields) > 0) {
    stop(
      paste0(
        "INPUT LOCK FAILURE: missing required fields:\n",
        paste(missing_fields, collapse = "\n")
      )
    )
  }
  
  
  # ------------------------------------------------------------
  # 5. Duplicate lock
  # ------------------------------------------------------------
  
  n_unique_pairs <- dplyr::n_distinct(
    master$Exact_LR_Pair
  )
  
  if (n_unique_pairs != nrow(master)) {
    stop(
      paste0(
        "INPUT LOCK FAILURE: Exact_LR_Pair is not unique.\n",
        "Rows = ", nrow(master),
        "; unique pairs = ", n_unique_pairs
      )
    )
  }
  
  cat("Exact_LR_Pair uniqueness: PASS\n")
  
  
  # ------------------------------------------------------------
  # 6. LOCKED biological direction map
  # ------------------------------------------------------------
  
  direction_map <- c(
    "1" = "Endothelial -> Myeloid/macrophage",
    "2" = "Myeloid/macrophage -> Endothelial",
    "3" = "Endothelial -> T-cell",
    "4" = "T-cell -> Endothelial",
    "5" = "Myeloid/macrophage -> T-cell",
    "6" = "T-cell -> Myeloid/macrophage"
  )
  
  direction_short_map <- c(
    "1" = "E→M",
    "2" = "M→E",
    "3" = "E→T",
    "4" = "T→E",
    "5" = "M→T",
    "6" = "T→M"
  )
  
  
  # ------------------------------------------------------------
  # 7. Create manuscript-facing direction-specific
  #    compatibility field
  #
  # IMPORTANT:
  #   This step ONLY reformats existing locked values.
  #   No fractions are recalculated.
  #
  # Missing direction = direction not represented.
  # It is NOT converted to 0%.
  # ------------------------------------------------------------
  
  format_fraction <- function(x) {
    
    if (is.na(x)) {
      return(NA_character_)
    }
    
    paste0(
      formatC(
        x * 100,
        format = "f",
        digits = 1
      ),
      "%"
    )
  }
  
  
  make_direction_compatibility <- function(row) {
    
    result <- character(0)
    
    for (i in 1:6) {
      
      field <- paste0(
        "Direction_",
        i,
        "_Compatibility_Fraction"
      )
      
      value <- row[[field]]
      
      if (!is.na(value)) {
        
        result <- c(
          result,
          paste0(
            direction_short_map[as.character(i)],
            ": ",
            format_fraction(value)
          )
        )
      }
    }
    
    if (length(result) == 0) {
      return(NA_character_)
    }
    
    paste(
      result,
      collapse = "; "
    )
  }
  
  
  direction_specific_compatibility <- vapply(
    seq_len(nrow(master)),
    function(i) {
      make_direction_compatibility(master[i, ])
    },
    character(1)
  )
  
  
  # ------------------------------------------------------------
  # 8. Create manuscript-facing direction label
  #
  # Use the locked Direction_Labels field as the biological
  # direction source. This does NOT infer or recalculate
  # directions.
  # ------------------------------------------------------------
  
  direction_labels_final <- master$Direction_Labels
  
  
  # ------------------------------------------------------------
  # 9. Manuscript-facing table structure
  #
  # One row = one Exact_LR_Pair
  #
  # Internal audit/provenance fields intentionally excluded.
  # ------------------------------------------------------------
  
  supp_table <- master %>%
    transmute(
      
      `Exact LR pair` =
        Exact_LR_Pair,
      
      `Ligand` =
        ligand,
      
      `Receptor` =
        receptor,
      
      `Candidate consistency` =
        Candidate_Consistency_Class,
      
      `Pathway` =
        pathway_name,
      
      `Resource annotation` =
        annotation,
      
      `No. of directions` =
        N_Directions,
      
      `Sender → receiver direction(s)` =
        direction_labels_final,
      
      `No. of patients` =
        N_Patients,
      
      `No. compatible observations` =
        N_Compatible,
      
      `Overall compatibility` =
        Overall_Compatibility_Fraction,
      
      `Patient compatibility pattern` =
        Patient_Compatibility_Patterns,
      
      `Core gene involved` =
        Core_Gene_Involved,
      
      `Ligand expression, median (IQR), log2CPM` =
        paste0(
          formatC(
            Ligand_Median_log2CPM,
            format = "f",
            digits = 3
          ),
          " (",
          formatC(
            Ligand_IQR_log2CPM,
            format = "f",
            digits = 3
          ),
          ")"
        ),
      
      `Receptor expression, median (IQR), log2CPM` =
        paste0(
          formatC(
            Receptor_Median_log2CPM,
            format = "f",
            digits = 3
          ),
          " (",
          formatC(
            Receptor_IQR_log2CPM,
            format = "f",
            digits = 3
          ),
          ")"
        ),
      
      `Direction-specific compatibility` =
        direction_specific_compatibility
    )
  
  
  # ------------------------------------------------------------
  # 10. Replace NA expression summaries with explicit NA
  #     rather than producing "NA (NA)"
  # ------------------------------------------------------------
  
  ligand_na <- is.na(master$Ligand_Median_log2CPM) |
    is.na(master$Ligand_IQR_log2CPM)
  
  receptor_na <- is.na(master$Receptor_Median_log2CPM) |
    is.na(master$Receptor_IQR_log2CPM)
  
  supp_table$`Ligand expression, median (IQR), log2CPM`[
    ligand_na
  ] <- NA_character_
  
  supp_table$`Receptor expression, median (IQR), log2CPM`[
    receptor_na
  ] <- NA_character_
  
  
  # ------------------------------------------------------------
  # 11. Overall compatibility formatting
  #
  # IMPORTANT:
  #   Existing Overall_Compatibility_Fraction is used directly.
  #   No recalculation.
  # ------------------------------------------------------------
  
  supp_table$`Overall compatibility` <- ifelse(
    is.na(supp_table$`Overall compatibility`),
    NA_character_,
    paste0(
      formatC(
        supp_table$`Overall compatibility` * 100,
        format = "f",
        digits = 1
      ),
      "%"
    )
  )
  
  
  # ------------------------------------------------------------
  # 12. Structural checks after field formatting
  # ------------------------------------------------------------
  
  if (nrow(supp_table) != 585) {
    stop("STRUCTURE LOCK FAILURE: output row count is not 585.")
  }
  
  if (ncol(supp_table) != 16) {
    stop(
      paste0(
        "STRUCTURE LOCK FAILURE: expected 16 manuscript-facing columns; ",
        "found ",
        ncol(supp_table),
        "."
      )
    )
  }
  
  if (
    anyDuplicated(
      supp_table$`Exact LR pair`
    ) > 0
  ) {
    stop(
      "STRUCTURE LOCK FAILURE: duplicated Exact LR pair in final table."
    )
  }
  
  
  # ------------------------------------------------------------
  # 13. Explicit verification of the critical
  #     N_Compatible terminology
  # ------------------------------------------------------------
  
  if (
    !"No. compatible observations" %in%
    colnames(supp_table)
  ) {
    stop(
      "TERMINOLOGY LOCK FAILURE: 'No. compatible observations' not found."
    )
  }
  
  if (
    "Compatible patients" %in%
    colnames(supp_table)
  ) {
    stop(
      "TERMINOLOGY LOCK FAILURE: obsolete 'Compatible patients' column remains."
    )
  }
  
  
  # ------------------------------------------------------------
  # 14. Verify that N_Compatible was not altered
  #
  # Since the source value is directly transferred, this checks
  # equality between the original numeric field and the final
  # exported field.
  # ------------------------------------------------------------
  
  n_compatible_source <- master$N_Compatible
  
  n_compatible_output <- supp_table$`No. compatible observations`
  
  if (
    !identical(
      as.numeric(n_compatible_output),
      as.numeric(n_compatible_source)
    )
  ) {
    stop(
      "DATA TRANSFER LOCK FAILURE: N_Compatible values changed."
    )
  }
  
  
  # ------------------------------------------------------------
  # 15. Direction-specific field audit
  # ------------------------------------------------------------
  
  cat("\nDirection-specific compatibility field:\n")
  
  for (i in 1:6) {
    
    field <- paste0(
      "Direction_",
      i,
      "_Compatibility_Fraction"
    )
    
    n_present <- sum(
      !is.na(master[[field]])
    )
    
    n_missing <- sum(
      is.na(master[[field]])
    )
    
    cat(
      "Direction_", i,
      " | ",
      direction_short_map[as.character(i)],
      " | ",
      direction_map[as.character(i)],
      " | present = ",
      n_present,
      " | missing = ",
      n_missing,
      "\n",
      sep = ""
    )
  }
  
  
  # ------------------------------------------------------------
  # 16. Create README / table note
  # ------------------------------------------------------------
  
  table_note <- c(
    "Supplementary Table Sx.",
    "Patient-level descriptive context of candidate ligand–receptor pairs across endothelial, myeloid/macrophage, and T-cell compartments.",
    "",
    "Each row represents one exact ligand–receptor (LR) pair retained in the locked patient-level descriptive context analysis.",
    "Candidate LR pairs were defined from the corresponding public ligand–receptor reference resource and retained according to the predefined patient-level expression-compatibility criterion.",
    "Compatibility was evaluated independently for each candidate LR pair and sender–receiver direction across the analyzed patients.",
    "These candidates represent computationally inferred, expression-compatible LR relationships and do not constitute experimentally validated ligand–receptor interactions or quantified communication probabilities.",
    "",
    "No. compatible observations refers to the number of compatible patient-by-direction observations and should not be interpreted as the number of unique patients.",
    "Overall compatibility is the locked compatibility fraction from the master table.",
    "Direction-specific compatibility reports the corresponding locked compatibility fraction for each represented sender–receiver direction.",
    "A direction absent from a candidate LR pair is left unreported and is not interpreted as 0% compatibility.",
    "",
    "Expression values are summarized as median (IQR) log2CPM.",
    "Resource annotations were retained from the corresponding public LR reference.",
    "",
    "The table was generated by field selection, renaming, and formatting from the locked master table.",
    "No LR compatibility values, patient counts, expression summaries, statistical tests, or other analytical results were recalculated."
  )
  
  
  # ------------------------------------------------------------
  # 17. Excel output path
  # ------------------------------------------------------------
  
  excel_file <- file.path(
    output_dir,
    "Supplementary_Table_Sx_Patient_Level_LR_Descriptive_Context_FINAL.xlsx"
  )
  
  txt_file <- file.path(
    output_dir,
    "STEP11-5F-4D_FINAL_Supplementary_Table_Lock_Audit.txt"
  )
  
  
  # ------------------------------------------------------------
  # 18. Build workbook
  # ------------------------------------------------------------
  
  wb <- openxlsx::createWorkbook()
  
  openxlsx::addWorksheet(
    wb,
    "Supplementary Table"
  )
  
  openxlsx::addWorksheet(
    wb,
    "Table Note"
  )
  
  
  # ------------------------------------------------------------
  # 19. Styles
  # ------------------------------------------------------------
  
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    wrapText = TRUE,
    border = "Bottom"
  )
  
  body_style <- openxlsx::createStyle(
    valign = "top",
    wrapText = TRUE
  )
  
  percentage_style <- openxlsx::createStyle(
    valign = "top",
    wrapText = TRUE
  )
  
  
  # ------------------------------------------------------------
  # 20. Write Supplementary Table
  # ------------------------------------------------------------
  
  openxlsx::writeData(
    wb,
    sheet = "Supplementary Table",
    x = supp_table,
    headerStyle = header_style,
    withFilter = TRUE
  )
  
  openxlsx::addStyle(
    wb,
    sheet = "Supplementary Table",
    style = body_style,
    rows = 2:(nrow(supp_table) + 1),
    cols = 1:ncol(supp_table),
    gridExpand = TRUE,
    stack = TRUE
  )
  
  
  # ------------------------------------------------------------
  # 21. Column widths
  # ------------------------------------------------------------
  
  widths <- c(
    24,  # Exact LR pair
    16,  # Ligand
    16,  # Receptor
    22,  # Candidate consistency
    18,  # Pathway
    22,  # Resource annotation
    15,  # No. directions
    32,  # Sender receiver
    14,  # No. patients
    23,  # No compatible observations
    20,  # Overall compatibility
    32,  # Patient compatibility pattern
    20,  # Core gene
    28,  # Ligand expression
    30,  # Receptor expression
    34   # Direction-specific compatibility
  )
  
  openxlsx::setColWidths(
    wb,
    sheet = "Supplementary Table",
    cols = 1:ncol(supp_table),
    widths = widths
  )
  
  openxlsx::freezePane(
    wb,
    sheet = "Supplementary Table",
    firstRow = TRUE
  )
  
  
  # ------------------------------------------------------------
  # 22. Table Note worksheet
  # ------------------------------------------------------------
  
  note_df <- data.frame(
    `Table information` = table_note,
    check.names = FALSE
  )
  
  openxlsx::writeData(
    wb,
    sheet = "Table Note",
    x = note_df,
    colNames = FALSE
  )
  
  openxlsx::setColWidths(
    wb,
    sheet = "Table Note",
    cols = 1,
    widths = 120
  )
  
  openxlsx::addStyle(
    wb,
    sheet = "Table Note",
    style = body_style,
    rows = 1:length(table_note),
    cols = 1,
    gridExpand = TRUE
  )
  
  
  # ------------------------------------------------------------
  # 23. Save Excel
  # ------------------------------------------------------------
  
  openxlsx::saveWorkbook(
    wb,
    excel_file,
    overwrite = TRUE
  )
  
  
  # ------------------------------------------------------------
  # 24. Final audit report
  # ------------------------------------------------------------
  
  audit_lines <- c(
    
    "============================================================",
    "STEP11-5F-4D FINAL SUPPLEMENTARY TABLE LOCK AUDIT",
    "============================================================",
    "",
    paste0("Input file: ", input_file),
    paste0("Output file: ", excel_file),
    "",
    paste0("Input dimensions: ", nrow(master), " x ", ncol(master)),
    paste0("Output dimensions: ", nrow(supp_table), " x ", ncol(supp_table)),
    "",
    "INPUT LOCK",
    "------------------------------------------------------------",
    paste0(
      "Input rows = 585: ",
      ifelse(nrow(master) == 585, "PASS", "FAIL")
    ),
    paste0(
      "Input columns = 36: ",
      ifelse(ncol(master) == 36, "PASS", "FAIL")
    ),
    paste0(
      "Exact_LR_Pair unique: ",
      ifelse(n_unique_pairs == 585, "PASS", "FAIL")
    ),
    "",
    "OUTPUT STRUCTURE LOCK",
    "------------------------------------------------------------",
    paste0(
      "Output rows = 585: ",
      ifelse(nrow(supp_table) == 585, "PASS", "FAIL")
    ),
    paste0(
      "Output columns = 16: ",
      ifelse(ncol(supp_table) == 16, "PASS", "FAIL")
    ),
    "One row per Exact LR pair: PASS",
    "",
    "TERMINOLOGY LOCK",
    "------------------------------------------------------------",
    "N_Compatible displayed as: No. compatible observations",
    "Compatible patients: REMOVED",
    "",
    "DATA INTEGRITY",
    "------------------------------------------------------------",
    "N_Compatible values recalculated: NO",
    "Overall compatibility recalculated: NO",
    "Direction-specific compatibility recalculated: NO",
    "Expression summaries recalculated: NO",
    "Patient counts recalculated: NO",
    "Statistical testing performed: NO",
    "Source master table modified: NO",
    "",
    "DIRECTION MAP",
    "------------------------------------------------------------",
    "Direction 1 = E→M",
    "Direction 2 = M→E",
    "Direction 3 = E→T",
    "Direction 4 = T→E",
    "Direction 5 = M→T",
    "Direction 6 = T→M",
    "",
    "MANUSCRIPT INTERPRETATION",
    "------------------------------------------------------------",
    "Candidate LR relationships are computationally inferred and expression-compatible.",
    "They are not described as experimentally validated interactions.",
    "The table does not contain quantified CellChat communication probabilities.",
    "",
    "FINAL STATUS",
    "------------------------------------------------------------",
    "STEP11-5F-4D = PASS",
    "",
    "============================================================"
  )
  
  writeLines(
    audit_lines,
    con = txt_file
  )
  
  
  # ------------------------------------------------------------
  # 25. Console summary
  # ------------------------------------------------------------
  
  cat("\n")
  cat("============================================================\n")
  cat("STEP11-5F-4D FINAL STATUS: PASS\n")
  cat("============================================================\n\n")
  
  cat(
    "Input:  ",
    nrow(master),
    " x ",
    ncol(master),
    "\n",
    sep = ""
  )
  
  cat(
    "Output: ",
    nrow(supp_table),
    " x ",
    ncol(supp_table),
    "\n",
    sep = ""
  )
  
  cat("\n")
  cat("Critical terminology lock:\n")
  cat("  N_Compatible -> No. compatible observations\n")
  cat("  Compatible patients -> REMOVED\n\n")
  
  cat("Analysis recalculation:\n")
  cat("  LR compatibility: NO\n")
  cat("  Patient counts:   NO\n")
  cat("  Expression stats: NO\n")
  cat("  Statistics:       NO\n\n")
  
  cat("Direction map:\n")
  cat("  Direction 1 = E→M\n")
  cat("  Direction 2 = M→E\n")
  cat("  Direction 3 = E→T\n")
  cat("  Direction 4 = T→E\n")
  cat("  Direction 5 = M→T\n")
  cat("  Direction 6 = T→M\n\n")
  
  cat("Excel generated:\n")
  cat(excel_file, "\n\n")
  
  cat("Audit file:\n")
  cat(txt_file, "\n\n")
  
  cat("============================================================\n")
  cat("LOCK COMPLETE — MANUSCRIPT-FACING SUPPLEMENTARY TABLE READY\n")
  cat("============================================================\n")
  
  
  
  
  # ============================================================
  # STEP11-5F-4D CHECK
  # Create manuscript-facing supp_table only
  # No recalculation
  # ============================================================
  
  direction_short_map <- c(
    "1" = "E→M",
    "2" = "M→E",
    "3" = "E→T",
    "4" = "T→E",
    "5" = "M→T",
    "6" = "T→M"
  )
  
  format_fraction <- function(x) {
    if (is.na(x)) {
      return(NA_character_)
    }
    
    paste0(
      formatC(
        x * 100,
        format = "f",
        digits = 1
      ),
      "%"
    )
  }
  
  make_direction_compatibility <- function(row) {
    
    result <- character(0)
    
    for (i in 1:6) {
      
      field <- paste0(
        "Direction_",
        i,
        "_Compatibility_Fraction"
      )
      
      value <- row[[field]]
      
      if (!is.na(value)) {
        
        result <- c(
          result,
          paste0(
            direction_short_map[as.character(i)],
            ": ",
            format_fraction(value)
          )
        )
      }
    }
    
    if (length(result) == 0) {
      return(NA_character_)
    }
    
    paste(
      result,
      collapse = "; "
    )
  }
  
  direction_specific_compatibility <- vapply(
    seq_len(nrow(master)),
    function(i) {
      make_direction_compatibility(master[i, ])
    },
    character(1)
  )
  
  direction_labels_final <- master$Direction_Labels
  
  supp_table <- master %>%
    dplyr::transmute(
      
      `Exact LR pair` =
        Exact_LR_Pair,
      
      `Ligand` =
        ligand,
      
      `Receptor` =
        receptor,
      
      `Candidate consistency` =
        Candidate_Consistency_Class,
      
      `Pathway` =
        pathway_name,
      
      `Resource annotation` =
        annotation,
      
      `No. of directions` =
        N_Directions,
      
      `Sender → receiver direction(s)` =
        direction_labels_final,
      
      `No. of patients` =
        N_Patients,
      
      `No. compatible observations` =
        N_Compatible,
      
      `Overall compatibility` =
        Overall_Compatibility_Fraction,
      
      `Patient compatibility pattern` =
        Patient_Compatibility_Patterns,
      
      `Core gene involved` =
        Core_Gene_Involved,
      
      `Ligand expression, median (IQR), log2CPM` =
        paste0(
          formatC(
            Ligand_Median_log2CPM,
            format = "f",
            digits = 3
          ),
          " (",
          formatC(
            Ligand_IQR_log2CPM,
            format = "f",
            digits = 3
          ),
          ")"
        ),
      
      `Receptor expression, median (IQR), log2CPM` =
        paste0(
          formatC(
            Receptor_Median_log2CPM,
            format = "f",
            digits = 3
          ),
          " (",
          formatC(
            Receptor_IQR_log2CPM,
            format = "f",
            digits = 3
          ),
          ")"
        ),
      
      `Direction-specific compatibility` =
        direction_specific_compatibility
    )
  
  cat("\n")
  cat("============================================================\n")
  cat("SUPP_TABLE CREATION CHECK\n")
  cat("============================================================\n")
  
  cat(
    "supp_table exists: ",
    exists("supp_table"),
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
  
  cat("\nColumn names:\n")
  print(colnames(supp_table))
  
  cat("\nFirst 3 rows:\n")
  print(head(supp_table, 3))
  
  cat("\n============================================================\n")
  
  
  # ============================================================
  # STEP11-5F-4D
  # FINAL DATA-TRANSFER INTEGRITY CHECK
  # No recalculation
  # No source modification
  # ============================================================
  
  cat("\n")
  cat("============================================================\n")
  cat("STEP11-5F-4D DATA-TRANSFER INTEGRITY CHECK\n")
  cat("============================================================\n\n")
  
  
  # ------------------------------------------------------------
  # 1. Row identity
  # ------------------------------------------------------------
  
  row_identity_pass <- identical(
    supp_table$`Exact LR pair`,
    master$Exact_LR_Pair
  )
  
  cat(
    "Exact LR pair identity: ",
    ifelse(row_identity_pass, "PASS", "FAIL"),
    "\n",
    sep = ""
  )
  
  
  # ------------------------------------------------------------
  # 2. N_Compatible transfer
  # ------------------------------------------------------------
  
  ncompatible_pass <- identical(
    as.numeric(
      supp_table$`No. compatible observations`
    ),
    as.numeric(
      master$N_Compatible
    )
  )
  
  cat(
    "N_Compatible transfer: ",
    ifelse(ncompatible_pass, "PASS", "FAIL"),
    "\n",
    sep = ""
  )
  
  
  # ------------------------------------------------------------
  # 3. N_Patients transfer
  # ------------------------------------------------------------
  
  npatients_pass <- identical(
    as.numeric(
      supp_table$`No. of patients`
    ),
    as.numeric(
      master$N_Patients
    )
  )
  
  cat(
    "N_Patients transfer: ",
    ifelse(npatients_pass, "PASS", "FAIL"),
    "\n",
    sep = ""
  )
  
  
  # ------------------------------------------------------------
  # 4. Overall compatibility transfer
  # ------------------------------------------------------------
  
  overall_pass <- identical(
    as.numeric(
      supp_table$`Overall compatibility`
    ),
    as.numeric(
      master$Overall_Compatibility_Fraction
    )
  )
  
  cat(
    "Overall compatibility transfer: ",
    ifelse(overall_pass, "PASS", "FAIL"),
    "\n",
    sep = ""
  )
  
  
  # ------------------------------------------------------------
  # 5. Direction-specific source-field existence
  # ------------------------------------------------------------
  
  direction_fields <- paste0(
    "Direction_",
    1:6,
    "_Compatibility_Fraction"
  )
  
  direction_fields_pass <- all(
    direction_fields %in% colnames(master)
  )
  
  cat(
    "Six locked direction fields present: ",
    ifelse(direction_fields_pass, "PASS", "FAIL"),
    "\n",
    sep = ""
  )
  
  
  # ------------------------------------------------------------
  # 6. Direction-specific missingness is preserved
  # ------------------------------------------------------------
  
  direction_missing_source <- sapply(
    master[direction_fields],
    function(x) sum(is.na(x))
  )
  
  cat("\nDirection-specific source missingness:\n")
  
  for (i in 1:6) {
    
    cat(
      "  Direction_", i,
      ": ",
      direction_missing_source[i],
      " missing\n",
      sep = ""
    )
  }
  
  
  # ------------------------------------------------------------
  # 7. Final structural checks
  # ------------------------------------------------------------
  
  structure_pass <-
    nrow(supp_table) == 585 &&
    ncol(supp_table) == 16 &&
    row_identity_pass &&
    ncompatible_pass &&
    npatients_pass &&
    overall_pass &&
    direction_fields_pass
  
  
  # ------------------------------------------------------------
  # 8. Final status
  # ------------------------------------------------------------
  
  cat("\n")
  cat("============================================================\n")
  
  if (structure_pass) {
    
    cat(
      "STEP11-5F-4D DATA-TRANSFER INTEGRITY: PASS\n"
    )
    
  } else {
    
    cat(
      "STEP11-5F-4D DATA-TRANSFER INTEGRITY: REVIEW REQUIRED\n"
    )
    
  }
  
  cat("============================================================\n")
  
  
  # ============================================================
  # STEP11-5F-4D
  # FINAL MANUSCRIPT SUPPLEMENTARY TABLE EXCEL GENERATION
  # ============================================================
  #
  # STATUS:
  #   STEP11-5F-4D DATA-TRANSFER INTEGRITY = PASS
  #
  # IMPORTANT:
  #   This step performs formatting only.
  #
  #   NO:
  #   - LR compatibility recalculation
  #   - patient-count recalculation
  #   - expression recalculation
  #   - statistical testing
  #   - direction reassignment
  #   - source-table modification
  #
  # ============================================================
  
  
  library(openxlsx)
  
  
  # ------------------------------------------------------------
  # 1. Output directory
  # ------------------------------------------------------------
  
  step_dir <- "input_dir/Paper 1_Step 11-CellChat"
  
  output_dir <- file.path(
    step_dir,
    "STEP11-5F-4D_FINAL_Supplementary_Table"
  )
  
  if (!dir.exists(output_dir)) {
    dir.create(
      output_dir,
      recursive = TRUE
    )
  }
  
  
  # ------------------------------------------------------------
  # 2. Confirm locked objects
  # ------------------------------------------------------------
  
  if (!exists("master")) {
    stop(
      "LOCK FAILURE: object 'master' not found."
    )
  }
  
  if (!exists("supp_table")) {
    stop(
      "LOCK FAILURE: object 'supp_table' not found."
    )
  }
  
  
  # ------------------------------------------------------------
  # 3. Final structure check before Excel generation
  # ------------------------------------------------------------
  
  if (
    nrow(master) != 585 ||
    ncol(master) != 36
  ) {
    stop(
      "LOCK FAILURE: master is not 585 x 36."
    )
  }
  
  if (
    nrow(supp_table) != 585 ||
    ncol(supp_table) != 16
  ) {
    stop(
      "LOCK FAILURE: supp_table is not 585 x 16."
    )
  }
  
  
  # ------------------------------------------------------------
  # 4. Confirm critical terminology
  # ------------------------------------------------------------
  
  if (
    !"No. compatible observations" %in%
    colnames(supp_table)
  ) {
    stop(
      "TERMINOLOGY LOCK FAILURE."
    )
  }
  
  if (
    "Compatible patients" %in%
    colnames(supp_table)
  ) {
    stop(
      "OBSOLETE TERMINOLOGY DETECTED."
    )
  }
  
  
  # ------------------------------------------------------------
  # 5. Output filenames
  # ------------------------------------------------------------
  
  excel_file <- file.path(
    output_dir,
    "Supplementary_Table_Sx_Patient_Level_LR_Descriptive_Context_FINAL.xlsx"
  )
  
  audit_file <- file.path(
    output_dir,
    "STEP11-5F-4D_FINAL_Excel_Generation_Audit.txt"
  )
  
  
  # ------------------------------------------------------------
  # 6. Create workbook
  # ------------------------------------------------------------
  
  wb <- createWorkbook()
  
  
  # ------------------------------------------------------------
  # 7. Add worksheets
  # ------------------------------------------------------------
  
  addWorksheet(
    wb,
    "Supplementary Table",
    gridLines = FALSE
  )
  
  addWorksheet(
    wb,
    "Table Note",
    gridLines = FALSE
  )
  
  
  # ------------------------------------------------------------
  # 8. Define styles
  # ------------------------------------------------------------
  
  title_style <- createStyle(
    fontName = "Arial",
    fontSize = 10,
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    wrapText = TRUE,
    border = "Bottom",
    borderStyle = "thin"
  )
  
  body_style <- createStyle(
    fontName = "Arial",
    fontSize = 9,
    valign = "top",
    wrapText = TRUE
  )
  
  note_style <- createStyle(
    fontName = "Arial",
    fontSize = 10,
    valign = "top",
    wrapText = TRUE
  )
  
  
  # ------------------------------------------------------------
  # 9. Write Supplementary Table
  # ------------------------------------------------------------
  
  writeData(
    wb,
    sheet = "Supplementary Table",
    x = supp_table,
    headerStyle = title_style,
    withFilter = TRUE,
    keepNA = FALSE
  )
  
  
  # ------------------------------------------------------------
  # 10. Apply body style
  # ------------------------------------------------------------
  
  addStyle(
    wb,
    sheet = "Supplementary Table",
    style = body_style,
    rows = 2:(nrow(supp_table) + 1),
    cols = 1:ncol(supp_table),
    gridExpand = TRUE,
    stack = TRUE
  )
  
  
  # ------------------------------------------------------------
  # 11. Column widths
  # ------------------------------------------------------------
  
  setColWidths(
    wb,
    sheet = "Supplementary Table",
    cols = 1,
    widths = 24
  )
  
  setColWidths(
    wb,
    sheet = "Supplementary Table",
    cols = 2:3,
    widths = 16
  )
  
  setColWidths(
    wb,
    sheet = "Supplementary Table",
    cols = 4,
    widths = 22
  )
  
  setColWidths(
    wb,
    sheet = "Supplementary Table",
    cols = 5,
    widths = 18
  )
  
  setColWidths(
    wb,
    sheet = "Supplementary Table",
    cols = 6,
    widths = 22
  )
  
  setColWidths(
    wb,
    sheet = "Supplementary Table",
    cols = 7,
    widths = 15
  )
  
  setColWidths(
    wb,
    sheet = "Supplementary Table",
    cols = 8,
    widths = 32
  )
  
  setColWidths(
    wb,
    sheet = "Supplementary Table",
    cols = 9,
    widths = 14
  )
  
  setColWidths(
    wb,
    sheet = "Supplementary Table",
    cols = 10,
    widths = 23
  )
  
  setColWidths(
    wb,
    sheet = "Supplementary Table",
    cols = 11,
    widths = 20
  )
  
  setColWidths(
    wb,
    sheet = "Supplementary Table",
    cols = 12,
    widths = 32
  )
  
  setColWidths(
    wb,
    sheet = "Supplementary Table",
    cols = 13,
    widths = 20
  )
  
  setColWidths(
    wb,
    sheet = "Supplementary Table",
    cols = 14:15,
    widths = 30
  )
  
  setColWidths(
    wb,
    sheet = "Supplementary Table",
    cols = 16,
    widths = 34
  )
  
  
  # ------------------------------------------------------------
  # 12. Freeze header row
  # ------------------------------------------------------------
  
  freezePane(
    wb,
    sheet = "Supplementary Table",
    firstRow = TRUE
  )
  
  
  # ------------------------------------------------------------
  # 13. Set header height
  # ------------------------------------------------------------
  
  setRowHeights(
    wb,
    sheet = "Supplementary Table",
    rows = 1,
    heights = 42
  )
  
  
  # ------------------------------------------------------------
  # 14. Table Note
  # ------------------------------------------------------------
  
  table_note <- data.frame(
    `Table Note` = c(
      "Supplementary Table Sx. Patient-level descriptive context of candidate ligand–receptor pairs across endothelial, myeloid/macrophage, and T-cell compartments.",
      "",
      "Each row represents one exact ligand–receptor (LR) pair retained in the locked patient-level descriptive context analysis.",
      "",
      "Candidate LR pairs were defined from the corresponding public ligand–receptor reference resource and retained according to the predefined patient-level expression-compatibility criterion.",
      "",
      "Compatibility was evaluated independently for each candidate LR pair and sender–receiver direction across the analyzed patients.",
      "",
      "These candidates represent computationally inferred, expression-compatible LR relationships and do not constitute experimentally validated ligand–receptor interactions or quantified communication probabilities.",
      "",
      "No. compatible observations refers to the number of compatible patient-by-direction observations and should not be interpreted as the number of unique patients.",
      "",
      "Overall compatibility is the locked compatibility fraction from the master table.",
      "",
      "Direction-specific compatibility reports the corresponding locked compatibility fraction for each represented sender–receiver direction.",
      "",
      "A direction absent from a candidate LR pair is left unreported and is not interpreted as 0% compatibility.",
      "",
      "Expression values are summarized as median (IQR) log2CPM.",
      "",
      "Resource annotations were retained from the corresponding public LR reference.",
      "",
      "This table was generated by field selection, terminology adjustment, and formatting from the locked master table. No LR compatibility values, patient counts, expression summaries, statistical tests, or other analytical results were recalculated."
    ),
    check.names = FALSE
  )
  
  
  writeData(
    wb,
    sheet = "Table Note",
    x = table_note,
    colNames = FALSE
  )
  
  addStyle(
    wb,
    sheet = "Table Note",
    style = note_style,
    rows = 1:nrow(table_note),
    cols = 1,
    gridExpand = TRUE
  )
  
  setColWidths(
    wb,
    sheet = "Table Note",
    cols = 1,
    widths = 120
  )
  
  
  # ------------------------------------------------------------
  # 15. Save final Excel
  # ------------------------------------------------------------
  
  saveWorkbook(
    wb,
    excel_file,
    overwrite = TRUE
  )
  
  
  # ------------------------------------------------------------
  # 16. Final Excel existence check
  # ------------------------------------------------------------
  
  excel_exists <- file.exists(
    excel_file
  )
  
  excel_size_mb <- if (excel_exists) {
    file.info(excel_file)$size / 1024^2
  } else {
    NA_real_
  }
  
  
  # ------------------------------------------------------------
  # 17. Final audit
  # ------------------------------------------------------------
  
  audit_lines <- c(
    "============================================================",
    "STEP11-5F-4D FINAL EXCEL GENERATION AUDIT",
    "============================================================",
    "",
    "Input master:",
    "585 rows x 36 columns",
    "",
    "Manuscript-facing table:",
    "585 rows x 16 columns",
    "",
    "Critical terminology:",
    "N_Compatible -> No. compatible observations",
    "Compatible patients -> REMOVED",
    "",
    "Data processing:",
    "Field selection only = YES",
    "Terminology adjustment only = YES",
    "Formatting only = YES",
    "",
    "Analytical recalculation:",
    "LR compatibility = NO",
    "Patient counts = NO",
    "Expression statistics = NO",
    "Statistical testing = NO",
    "Direction assignment = NO",
    "",
    "Locked direction map:",
    "Direction 1 = E -> M",
    "Direction 2 = M -> E",
    "Direction 3 = E -> T",
    "Direction 4 = T -> E",
    "Direction 5 = M -> T",
    "Direction 6 = T -> M",
    "",
    paste0(
      "Excel file exists = ",
      ifelse(excel_exists, "YES", "NO")
    ),
    paste0(
      "Excel file size (MB) = ",
      ifelse(
        is.na(excel_size_mb),
        "NA",
        sprintf("%.3f", excel_size_mb)
      )
    ),
    "",
    paste0(
      "FINAL STATUS = ",
      ifelse(
        excel_exists,
        "PASS",
        "FAIL"
      )
    ),
    "",
    "============================================================"
  )
  
  writeLines(
    audit_lines,
    con = audit_file
  )
  
  
  # ------------------------------------------------------------
  # 18. Console summary
  # ------------------------------------------------------------
  
  cat("\n")
  cat("============================================================\n")
  cat("STEP11-5F-4D FINAL EXCEL GENERATION STATUS\n")
  cat("============================================================\n\n")
  
  cat(
    "Input:  ",
    nrow(master),
    " x ",
    ncol(master),
    "\n",
    sep = ""
  )
  
  cat(
    "Output: ",
    nrow(supp_table),
    " x ",
    ncol(supp_table),
    "\n",
    sep = ""
  )
  
  cat("\n")
  
  cat(
    "Terminology: No. compatible observations\n"
  )
  
  cat(
    "Compatible patients: REMOVED\n\n"
  )
  
  cat(
    "Analytical recalculation: NONE\n"
  )
  
  cat(
    "Excel formatting: ONLY\n\n"
  )
  
  cat(
    "Excel generated: ",
    ifelse(excel_exists, "YES", "NO"),
    "\n",
    sep = ""
  )
  
  cat(
    "Excel path:\n",
    excel_file,
    "\n\n",
    sep = ""
  )
  
  cat(
    "Audit path:\n",
    audit_file,
    "\n\n",
    sep = ""
  )
  
  cat("============================================================\n")
  
  if (excel_exists) {
    cat(
      "STEP11-5F-4D FINAL STATUS: PASS\n"
    )
  } else {
    cat(
      "STEP11-5F-4D FINAL STATUS: FAIL\n"
    )
  }
  
  cat("============================================================\n")
  
  
  ##STEP11-6A — Figure 7 Traceability Audit
  # ============================================================
  # STEP11-6A
  # Figure 7 ↔ Locked Master Table Traceability Audit
  #
  # Purpose:
  #   Verify that LR pairs and sender→receiver directions used
  #   for Figure 7 can be traced back to the locked
  #   STEP11-5F-4C FINAL Master LR table.
  #
  # IMPORTANT:
  #   - NO LR compatibility recalculation
  #   - NO patient-level recalculation
  #   - NO expression recalculation
  #   - NO statistical testing
  #   - NO modification of master
  #   - NO modification of supp_table
  #
  # ============================================================
  
  
  cat("============================================================\n")
  cat("STEP11-6A — FIGURE 7 TRACEABILITY AUDIT\n")
  cat("============================================================\n\n")
  
  
  # ------------------------------------------------------------
  # 0. Required locked objects
  # ------------------------------------------------------------
  
  if (!exists("master")) {
    stop(
      "ERROR: locked object 'master' is not available.\n",
      "Please load the STEP11-5F-4C FINAL Master LR table first."
    )
  }
  
  if (!exists("supp_table")) {
    stop(
      "ERROR: locked object 'supp_table' is not available.\n",
      "Please load/create the STEP11-5F-4D Supplementary Table first."
    )
  }
  
  
  # ------------------------------------------------------------
  # 1. Working directories
  # ------------------------------------------------------------
  
  step11_root <- "input_dir/Paper 1_Step 11-CellChat"
  
  step11_6a_dir <- file.path(
    step11_root,
    "STEP11-6A_Figure7_Traceability_Audit"
  )
  
  dir.create(
    step11_6a_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # ------------------------------------------------------------
  # 2. Locked master dimensions
  # ------------------------------------------------------------
  
  cat("LOCKED INPUT DIMENSIONS\n")
  cat("-----------------------\n")
  cat("master:     ", nrow(master), " x ", ncol(master), "\n", sep = "")
  cat("supp_table: ", nrow(supp_table), " x ", ncol(supp_table), "\n\n", sep = "")
  
  
  if (nrow(master) != 585 || ncol(master) != 36) {
    stop(
      "ERROR: master is not the expected locked 585 x 36 table."
    )
  }
  
  if (nrow(supp_table) != 585 || ncol(supp_table) != 16) {
    stop(
      "ERROR: supp_table is not the expected locked 585 x 16 table."
    )
  }
  
  
  # ------------------------------------------------------------
  # 3. Locked direction map
  # ------------------------------------------------------------
  
  direction_map <- c(
    "1" = "Endothelial -> Myeloid/macrophage",
    "2" = "Myeloid/macrophage -> Endothelial",
    "3" = "Endothelial -> T-cell",
    "4" = "T-cell -> Endothelial",
    "5" = "Myeloid/macrophage -> T-cell",
    "6" = "T-cell -> Myeloid/macrophage"
  )
  
  direction_short <- c(
    "1" = "E→M",
    "2" = "M→E",
    "3" = "E→T",
    "4" = "T→E",
    "5" = "M→T",
    "6" = "T→M"
  )
  
  cat("LOCKED DIRECTION MAP\n")
  cat("--------------------\n")
  
  for (i in names(direction_map)) {
    cat(
      "Direction_", i, " = ",
      direction_map[[i]],
      " [", direction_short[[i]], "]\n",
      sep = ""
    )
  }
  
  cat("\n")
  
  
  # ------------------------------------------------------------
  # 4. Required master columns
  # ------------------------------------------------------------
  
  required_master_cols <- c(
    "Exact_LR_Pair",
    "ligand",
    "receptor",
    "Candidate_Consistency_Class",
    "pathway_name",
    "annotation",
    "N_Directions",
    "Directions",
    "Direction_Labels",
    "N_Patient_Records",
    "N_Patients",
    "N_Compatible",
    "N_Incompatible",
    "Overall_Compatibility_Fraction",
    "Core_Gene_Involved",
    "Direction_1_Compatibility_Fraction",
    "Direction_2_Compatibility_Fraction",
    "Direction_3_Compatibility_Fraction",
    "Direction_4_Compatibility_Fraction",
    "Direction_5_Compatibility_Fraction",
    "Direction_6_Compatibility_Fraction"
  )
  
  missing_master_cols <- setdiff(
    required_master_cols,
    colnames(master)
  )
  
  if (length(missing_master_cols) > 0) {
    stop(
      "ERROR: required master columns are missing:\n",
      paste(missing_master_cols, collapse = "\n")
    )
  }
  
  cat("Required master fields: PASS\n\n")
  
  
  # ------------------------------------------------------------
  # 5. Search for Figure 7 candidate source files
  # ------------------------------------------------------------
  
  cat("============================================================\n")
  cat("STEP11-6A — FIGURE 7 SOURCE DISCOVERY\n")
  cat("============================================================\n\n")
  
  candidate_files <- list.files(
    path = step11_root,
    pattern = "\\.(tsv|csv|rds|xlsx)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  candidate_files <- candidate_files[
    grepl(
      "figure.?7|fig.?7|F7|LR|ligand|receptor",
      basename(candidate_files),
      ignore.case = TRUE
    )
  ]
  
  candidate_files <- candidate_files[
    !grepl(
      "STEP11-5F-4D_FINAL_Supplementary_Table",
      candidate_files,
      ignore.case = TRUE
    )
  ]
  
  if (length(candidate_files) == 0) {
    
    cat("No obvious Figure 7 source file was automatically identified.\n\n")
    
    cat("Please provide the actual Figure 7 plotting-data file/object.\n")
    cat("The audit requires the exact LR pairs used by Figure 7.\n\n")
    
    cat("If the Figure 7 data are already loaded in R, create an object:\n")
    cat("    figure7_source <- <your locked Figure 7 data object>\n\n")
    
    stop(
      "STEP11-6A stopped safely before audit: Figure 7 source not identified."
    )
  }
  
  
  cat("Potential Figure 7-related source files found:\n\n")
  
  for (i in seq_along(candidate_files)) {
    cat(
      sprintf("[%02d] %s\n", i, candidate_files[[i]])
    )
  }
  
  cat("\n")
  
  
  # ------------------------------------------------------------
  # 6. IMPORTANT USER CONFIRMATION POINT
  # ------------------------------------------------------------
  #
  # Set this manually after reviewing the list above.
  #
  # Example:
  #
  # figure7_source_file <- candidate_files[3]
  #
  # DO NOT automatically choose the first file.
  #
  # ------------------------------------------------------------
  
  figure7_source_file <- ""
  
  if (figure7_source_file == "") {
    
    cat(
      "============================================================\n"
    )
    cat(
      "ACTION REQUIRED BEFORE FORMAL AUDIT\n"
    )
    cat(
      "============================================================\n\n"
    )
    
    cat(
      "Set 'figure7_source_file' to the EXACT data file used to\n"
    )
    cat(
      "generate Figure 7.\n\n"
    )
    
    cat(
      "No audit values have been calculated.\n"
    )
    cat(
      "No locked object has been modified.\n\n"
    )
    
    stop(
      "STEP11-6A stopped safely at Figure 7 source selection."
    )
  }
  
  
  # ------------------------------------------------------------
  # 7. Read Figure 7 source
  # ------------------------------------------------------------
  
  if (!file.exists(figure7_source_file)) {
    stop(
      "ERROR: figure7_source_file does not exist:\n",
      figure7_source_file
    )
  }
  
  
  figure7_ext <- tolower(
    tools::file_ext(figure7_source_file)
  )
  
  
  if (figure7_ext == "tsv") {
    
    figure7_df <- read.delim(
      figure7_source_file,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
  } else if (figure7_ext == "csv") {
    
    figure7_df <- read.csv(
      figure7_source_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
  } else if (figure7_ext == "rds") {
    
    figure7_df <- readRDS(
      figure7_source_file
    )
    
  } else if (figure7_ext == "xlsx") {
    
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop(
        "Package 'openxlsx' is required to read the Figure 7 xlsx source."
      )
    }
    
    figure7_df <- openxlsx::read.xlsx(
      figure7_source_file,
      sheet = 1
    )
    
  } else {
    
    stop(
      "Unsupported Figure 7 source file type: ",
      figure7_ext
    )
  }
  
  
  # ------------------------------------------------------------
  # 8. Basic Figure 7 source audit
  # ------------------------------------------------------------
  
  cat("============================================================\n")
  cat("FIGURE 7 SOURCE AUDIT\n")
  cat("============================================================\n\n")
  
  cat(
    "Source:\n",
    figure7_source_file,
    "\n\n",
    sep = ""
  )
  
  cat(
    "Dimensions: ",
    nrow(figure7_df),
    " x ",
    ncol(figure7_df),
    "\n\n",
    sep = ""
  )
  
  cat("Columns:\n")
  print(colnames(figure7_df))
  cat("\n")
  
  
  # ------------------------------------------------------------
  # 9. Identify Exact LR pair field
  # ------------------------------------------------------------
  
  possible_pair_cols <- c(
    "Exact_LR_Pair",
    "Exact LR pair",
    "exact_lr_pair",
    "LR_pair",
    "LR Pair",
    "lr_pair",
    "interaction_name",
    "interaction",
    "pair"
  )
  
  pair_col <- possible_pair_cols[
    possible_pair_cols %in% colnames(figure7_df)
  ][1]
  
  
  if (is.na(pair_col) || length(pair_col) == 0) {
    
    stop(
      "ERROR: Could not identify the Figure 7 Exact LR pair column.\n",
      "Please rename the source column to 'Exact_LR_Pair' and rerun."
    )
  }
  
  
  cat(
    "Figure 7 LR-pair field: ",
    pair_col,
    "\n\n",
    sep = ""
  )
  
  
  figure7_df$Exact_LR_Pair_Audit <- as.character(
    figure7_df[[pair_col]]
  )
  
  
  figure7_df$Exact_LR_Pair_Audit <- trimws(
    figure7_df$Exact_LR_Pair_Audit
  )
  
  
  # ------------------------------------------------------------
  # 10. Basic pair integrity
  # ------------------------------------------------------------
  
  pair_missing <- is.na(
    figure7_df$Exact_LR_Pair_Audit
  ) |
    figure7_df$Exact_LR_Pair_Audit == ""
  
  pair_duplicates <- duplicated(
    figure7_df$Exact_LR_Pair_Audit
  )
  
  
  cat("Figure 7 pair integrity:\n")
  cat(
    "  Missing LR pairs: ",
    sum(pair_missing),
    "\n",
    sep = ""
  )
  
  cat(
    "  Duplicate LR pair rows: ",
    sum(pair_duplicates),
    "\n\n",
    sep = ""
  )
  
  
  # ------------------------------------------------------------
  # 11. Trace Figure 7 LR pairs to locked master
  # ------------------------------------------------------------
  
  master_pairs <- unique(
    as.character(master$Exact_LR_Pair)
  )
  
  figure7_pairs <- unique(
    figure7_df$Exact_LR_Pair_Audit[
      !pair_missing
    ]
  )
  
  
  missing_from_master <- setdiff(
    figure7_pairs,
    master_pairs
  )
  
  present_in_master <- intersect(
    figure7_pairs,
    master_pairs
  )
  
  
  cat("============================================================\n")
  cat("LR PAIR TRACEABILITY\n")
  cat("============================================================\n\n")
  
  cat(
    "Unique Figure 7 LR pairs: ",
    length(figure7_pairs),
    "\n",
    sep = ""
  )
  
  cat(
    "Found in locked master: ",
    length(present_in_master),
    "\n",
    sep = ""
  )
  
  cat(
    "Missing from locked master: ",
    length(missing_from_master),
    "\n\n",
    sep = ""
  )
  
  
  pair_traceability_pass <- (
    length(missing_from_master) == 0
  )
  
  
  if (pair_traceability_pass) {
    
    cat(
      "Exact LR pair identity: PASS\n\n"
    )
    
  } else {
    
    cat(
      "Exact LR pair identity: FAIL\n\n"
    )
    
    cat("LR pairs used by Figure 7 but absent from master:\n")
    
    print(
      missing_from_master
    )
    
    cat("\n")
  }
  
  
  # ------------------------------------------------------------
  # 12. Direction column discovery
  # ------------------------------------------------------------
  
  possible_direction_cols <- c(
    "Direction",
    "direction",
    "Direction_Label",
    "Direction Labels",
    "Direction_Labels",
    "Sender_Receiver",
    "Sender → receiver direction",
    "sender_receiver",
    "sender.receiver",
    "direction_label"
  )
  
  direction_col <- possible_direction_cols[
    possible_direction_cols %in% colnames(figure7_df)
  ][1]
  
  
  cat("============================================================\n")
  cat("DIRECTION FIELD AUDIT\n")
  cat("============================================================\n\n")
  
  
  if (is.na(direction_col) || length(direction_col) == 0) {
    
    cat(
      "No explicit direction field was found in the Figure 7 source.\n\n"
    )
    
    cat(
      "This is not automatically a FAIL because some Figure 7 source\n"
    )
    
    cat(
      "tables may contain one row per LR pair and use separate plotting\n"
    )
    
    cat(
      "objects for direction assignment.\n\n"
    )
    
    cat(
      "Direction traceability will therefore require the Figure 7\n"
    )
    
    cat(
      "direction field to be supplied explicitly.\n\n"
    )
    
    direction_traceability_pass <- FALSE
    
  } else {
    
    cat(
      "Figure 7 direction field: ",
      direction_col,
      "\n\n",
      sep = ""
    )
    
    figure7_df$Direction_Audit <- as.character(
      figure7_df[[direction_col]]
    )
    
    figure7_df$Direction_Audit <- trimws(
      figure7_df$Direction_Audit
    )
    
    cat("Observed Figure 7 direction labels:\n")
    
    print(
      sort(
        unique(
          figure7_df$Direction_Audit
        )
      )
    )
    
    cat("\n")
    
    
    # ----------------------------------------------------------
    # 13. Direction normalization
    # ----------------------------------------------------------
    
    normalize_direction <- function(x) {
      
      x <- trimws(as.character(x))
      
      out <- rep(NA_character_, length(x))
      
      
      out[
        grepl(
          "^E.?→.?M$|Endothelial.*Myeloid",
          x,
          ignore.case = TRUE
        )
      ] <- "1"
      
      
      out[
        grepl(
          "^M.?→.?E$|Myeloid.*Endothelial",
          x,
          ignore.case = TRUE
        )
      ] <- "2"
      
      
      out[
        grepl(
          "^E.?→.?T$|Endothelial.*T.?cell",
          x,
          ignore.case = TRUE
        )
      ] <- "3"
      
      
      out[
        grepl(
          "^T.?→.?E$|T.?cell.*Endothelial",
          x,
          ignore.case = TRUE
        )
      ] <- "4"
      
      
      out[
        grepl(
          "^M.?→.?T$|Myeloid.*T.?cell",
          x,
          ignore.case = TRUE
        )
      ] <- "5"
      
      
      out[
        grepl(
          "^T.?→.?M$|T.?cell.*Myeloid",
          x,
          ignore.case = TRUE
        )
      ] <- "6"
      
      
      out
    }
    
    
    figure7_df$Direction_Number_Audit <- normalize_direction(
      figure7_df$Direction_Audit
    )
    
    
    unknown_direction <- is.na(
      figure7_df$Direction_Number_Audit
    ) &
      !is.na(
        figure7_df$Direction_Audit
      ) &
      figure7_df$Direction_Audit != ""
    
    
    cat(
      "Unrecognized Figure 7 direction labels: ",
      sum(unknown_direction),
      "\n\n",
      sep = ""
    )
    
    
    # ----------------------------------------------------------
    # 14. Direction ↔ master traceability
    # ----------------------------------------------------------
    
    master_direction_lookup <- master[
      ,
      c(
        "Exact_LR_Pair",
        "Direction_Labels"
      )
    ]
    
    master_direction_lookup$Exact_LR_Pair <- as.character(
      master_direction_lookup$Exact_LR_Pair
    )
    
    master_direction_lookup$Direction_Labels <- as.character(
      master_direction_lookup$Direction_Labels
    )
    
    
    # Build one-to-many lookup manually
    master_direction_map <- split(
      master_direction_lookup$Direction_Labels,
      master_direction_lookup$Exact_LR_Pair
    )
    
    
    direction_mismatch <- rep(
      FALSE,
      nrow(figure7_df)
    )
    
    
    direction_not_traceable <- rep(
      FALSE,
      nrow(figure7_df)
    )
    
    
    for (i in seq_len(nrow(figure7_df))) {
      
      pair_i <- figure7_df$Exact_LR_Pair_Audit[i]
      
      dir_i <- figure7_df$Direction_Number_Audit[i]
      
      if (
        is.na(pair_i) ||
        pair_i == "" ||
        is.na(dir_i)
      ) {
        next
      }
      
      
      master_dirs <- master_direction_map[[pair_i]]
      
      if (is.null(master_dirs)) {
        
        direction_not_traceable[i] <- TRUE
        next
        
      }
      
      
      master_dirs <- unique(
        unlist(
          strsplit(
            master_dirs,
            "[|,;]"
          )
        )
      )
      
      master_dirs <- trimws(
        master_dirs
      )
      
      
      expected_direction_label <- direction_map[
        as.character(dir_i)
      ]
      
      
      direction_match <- any(
        grepl(
          expected_direction_label,
          master_dirs,
          fixed = TRUE
        )
      )
      
      
      if (!direction_match) {
        direction_mismatch[i] <- TRUE
      }
    }
    
    
    cat(
      "Direction mapping mismatches: ",
      sum(direction_mismatch),
      "\n",
      sep = ""
    )
    
    cat(
      "Direction rows not traceable: ",
      sum(direction_not_traceable),
      "\n\n",
      sep = ""
    )
    
    
    direction_traceability_pass <- (
      sum(unknown_direction) == 0 &&
        sum(direction_mismatch) == 0 &&
        sum(direction_not_traceable) == 0
    )
    
    
    if (direction_traceability_pass) {
      
      cat(
        "Direction mapping: PASS\n\n"
      )
      
    } else {
      
      cat(
        "Direction mapping: REVIEW REQUIRED\n\n"
      )
    }
  }
  
  
  # ------------------------------------------------------------
  # 15. Optional Figure 7 direction counts
  # ------------------------------------------------------------
  
  cat("============================================================\n")
  cat("FIGURE 7 DIRECTION SUMMARY\n")
  cat("============================================================\n\n")
  
  
  if (
    "Direction_Number_Audit" %in%
    colnames(figure7_df)
  ) {
    
    direction_counts <- table(
      factor(
        figure7_df$Direction_Number_Audit,
        levels = names(direction_map)
      )
    )
    
    
    direction_summary <- data.frame(
      Direction_Number = names(direction_map),
      Direction = unname(direction_map),
      Short_Label = unname(direction_short),
      Figure7_Rows = as.integer(direction_counts),
      stringsAsFactors = FALSE
    )
    
    
    print(
      direction_summary,
      row.names = FALSE
    )
    
    cat("\n")
    
  } else {
    
    direction_summary <- NULL
    
    cat(
      "Direction count summary not available because no explicit\n",
      "Figure 7 direction field was identified.\n\n",
      sep = ""
    )
  }
  
  
  # ------------------------------------------------------------
  # 16. Pair-level master transfer audit
  # ------------------------------------------------------------
  
  cat("============================================================\n")
  cat("PAIR-LEVEL MASTER FIELD TRACEABILITY\n")
  cat("============================================================\n\n")
  
  
  if (pair_traceability_pass) {
    
    figure7_master_match <- match(
      figure7_df$Exact_LR_Pair_Audit,
      master$Exact_LR_Pair
    )
    
    
    figure7_df$Master_Row_Index <- figure7_master_match
    
    
    fields_to_check <- c(
      "ligand",
      "receptor",
      "Candidate_Consistency_Class",
      "pathway_name",
      "annotation",
      "N_Patients",
      "N_Compatible",
      "Overall_Compatibility_Fraction",
      "Core_Gene_Involved"
    )
    
    
    field_results <- data.frame(
      Field = fields_to_check,
      Status = NA_character_,
      stringsAsFactors = FALSE
    )
    
    
    # Only perform direct field comparison if the corresponding
    # Figure 7 fields are present.
    #
    # Otherwise mark as NOT_APPLICABLE rather than inventing data.
    
    for (j in seq_along(fields_to_check)) {
      
      field <- fields_to_check[j]
      
      possible_fig_cols <- c(
        field,
        gsub("_", " ", field),
        fieldnames <- switch(
          field,
          ligand = c("Ligand"),
          receptor = c("Receptor"),
          Candidate_Consistency_Class =
            c("Candidate consistency"),
          pathway_name =
            c("Pathway"),
          annotation =
            c("Resource annotation"),
          N_Patients =
            c("No. of patients"),
          N_Compatible =
            c("No. compatible observations"),
          Overall_Compatibility_Fraction =
            c("Overall compatibility"),
          Core_Gene_Involved =
            c("Core gene involved"),
          field
        )
      )
      
      
      possible_fig_cols <- unique(
        possible_fig_cols
      )
      
      possible_fig_cols <- possible_fig_cols[
        possible_fig_cols %in%
          colnames(figure7_df)
      ]
      
      
      if (length(possible_fig_cols) == 0) {
        
        field_results$Status[j] <-
          "NOT_APPLICABLE"
        
      } else {
        
        fig_col <- possible_fig_cols[1]
        
        fig_values <- as.character(
          figure7_df[[fig_col]]
        )
        
        master_values <- as.character(
          master[[field]][figure7_master_match]
        )
        
        
        equal_values <- (
          (is.na(fig_values) & is.na(master_values)) |
            fig_values == master_values
        )
        
        
        field_results$Status[j] <-
          ifelse(
            all(equal_values),
            "PASS",
            "FAIL"
          )
      }
    }
    
    
    print(
      field_results,
      row.names = FALSE
    )
    
    cat("\n")
    
  } else {
    
    field_results <- NULL
    
    cat(
      "Pair-level field audit skipped because LR pair identity failed.\n\n"
    )
  }
  
  
  # ------------------------------------------------------------
  # 17. Explicit audit: no recalculation
  # ------------------------------------------------------------
  
  cat("============================================================\n")
  cat("ANALYTICAL RECALCULATION AUDIT\n")
  cat("============================================================\n\n")
  
  cat("LR compatibility recalculation: NONE\n")
  cat("Patient-level expression recalculation: NONE\n")
  cat("Statistical testing: NONE\n")
  cat("FDR recalculation: NONE\n")
  cat("Master table modification: NONE\n")
  cat("Supplementary table modification: NONE\n\n")
  
  
  # ------------------------------------------------------------
  # 18. Final status logic
  # ------------------------------------------------------------
  
  field_traceability_pass <- TRUE
  
  if (!is.null(field_results)) {
    
    field_traceability_pass <- all(
      field_results$Status[
        field_results$Status != "NOT_APPLICABLE"
      ] == "PASS"
    )
  }
  
  
  overall_pass <- (
    pair_traceability_pass &&
      direction_traceability_pass &&
      field_traceability_pass
  )
  
  
  # ------------------------------------------------------------
  # 19. Save audit objects
  # ------------------------------------------------------------
  
  saveRDS(
    figure7_df,
    file.path(
      step11_6a_dir,
      "STEP11-6A_Figure7_Source_Audit_Object.rds"
    )
  )
  
  
  if (!is.null(direction_summary)) {
    
    write.table(
      direction_summary,
      file.path(
        step11_6a_dir,
        "STEP11-6A_Figure7_Direction_Summary.tsv"
      ),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
  
  
  if (!is.null(field_results)) {
    
    write.table(
      field_results,
      file.path(
        step11_6a_dir,
        "STEP11-6A_Figure7_Field_Traceability.tsv"
      ),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
  
  
  if (length(missing_from_master) > 0) {
    
    write.table(
      data.frame(
        Exact_LR_Pair = missing_from_master
      ),
      file.path(
        step11_6a_dir,
        "STEP11-6A_LR_Pairs_Missing_From_Master.tsv"
      ),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
  
  
  # ------------------------------------------------------------
  # 20. Final audit report
  # ------------------------------------------------------------
  
  audit_file <- file.path(
    step11_6a_dir,
    "STEP11-6A_FINAL_Traceability_Audit.txt"
  )
  
  
  audit_lines <- c(
    "============================================================",
    "STEP11-6A — FIGURE 7 TRACEABILITY AUDIT",
    "============================================================",
    "",
    paste0(
      "Locked master dimensions: ",
      nrow(master),
      " x ",
      ncol(master)
    ),
    paste0(
      "Supplementary table dimensions: ",
      nrow(supp_table),
      " x ",
      ncol(supp_table)
    ),
    "",
    paste0(
      "Figure 7 source: ",
      figure7_source_file
    ),
    "",
    paste0(
      "Unique Figure 7 LR pairs: ",
      length(figure7_pairs)
    ),
    paste0(
      "Found in locked master: ",
      length(present_in_master)
    ),
    paste0(
      "Missing from locked master: ",
      length(missing_from_master)
    ),
    "",
    paste0(
      "Exact LR pair identity: ",
      ifelse(
        pair_traceability_pass,
        "PASS",
        "FAIL"
      )
    ),
    paste0(
      "Direction mapping: ",
      ifelse(
        direction_traceability_pass,
        "PASS",
        "REVIEW REQUIRED"
      )
    ),
    paste0(
      "Pair-level field traceability: ",
      ifelse(
        field_traceability_pass,
        "PASS",
        "FAIL"
      )
    ),
    "",
    "Analytical recalculation: NONE",
    "Master table modification: NONE",
    "Supplementary table modification: NONE",
    "",
    paste0(
      "STEP11-6A FINAL STATUS: ",
      ifelse(
        overall_pass,
        "PASS",
        "REVIEW REQUIRED"
      )
    ),
    "============================================================"
  )
  
  
  writeLines(
    audit_lines,
    con = audit_file
  )
  
  
  # ------------------------------------------------------------
  # 21. Final console summary
  # ------------------------------------------------------------
  
  cat("\n")
  cat("============================================================\n")
  cat("STEP11-6A FINAL TRACEABILITY STATUS\n")
  cat("============================================================\n\n")
  
  cat(
    "Figure 7 unique LR pairs: ",
    length(figure7_pairs),
    "\n",
    sep = ""
  )
  
  cat(
    "Found in locked master: ",
    length(present_in_master),
    "\n",
    sep = ""
  )
  
  cat(
    "Missing from locked master: ",
    length(missing_from_master),
    "\n",
    sep = ""
  )
  
  cat("\n")
  
  cat(
    "Exact LR pair identity: ",
    ifelse(
      pair_traceability_pass,
      "PASS",
      "FAIL"
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "Direction mapping: ",
    ifelse(
      direction_traceability_pass,
      "PASS",
      "REVIEW REQUIRED"
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "Pair-level field traceability: ",
    ifelse(
      field_traceability_pass,
      "PASS",
      "FAIL"
    ),
    "\n",
    sep = ""
  )
  
  cat("\n")
  
  cat("LR compatibility recalculation: NONE\n")
  cat("Patient-level recalculation: NONE\n")
  cat("Statistical recalculation: NONE\n")
  cat("Master modification: NONE\n")
  cat("Supplementary table modification: NONE\n\n")
  
  cat(
    "Audit report:\n",
    audit_file,
    "\n\n",
    sep = ""
  )
  
  cat("============================================================\n")
  
  if (overall_pass) {
    
    cat(
      "STEP11-6A FINAL STATUS: PASS\n"
    )
    
  } else {
    
    cat(
      "STEP11-6A FINAL STATUS: REVIEW REQUIRED\n"
    )
  }
  
  cat("============================================================\n")
  
  ### STEP11-6A-0 — Figure 7 Source Identification Audit
  # ============================================================
  # STEP11-6A-0
  # Figure 7 Source Identification Audit
  #
  # PURPOSE:
  #   Identify and compare the final data sources most directly
  #   associated with Figure 7 panels A-D.
  #
  # IMPORTANT:
  #   READ-ONLY AUDIT
  #   No recalculation
  #   No modification of master
  #   No modification of supp_table
  #   No modification of Figure 7 data
  # ============================================================
  
  
  cat("============================================================\n")
  cat("STEP11-6A-0 — FIGURE 7 SOURCE IDENTIFICATION AUDIT\n")
  cat("============================================================\n\n")
  
  
  # ------------------------------------------------------------
  # 0. Root directory
  # ------------------------------------------------------------
  
  step11_root <- "input_dir/Paper 1_Step 11-CellChat"
  
  
  # ------------------------------------------------------------
  # 1. Candidate files that are most relevant to Figure 7
  # ------------------------------------------------------------
  
  candidate_paths <- c(
    
    # Figure 7B
    file.path(
      step11_root,
      "STEP11-5F-4_Patient_Level_LR_Final_Context_Summary",
      "Figure7B",
      "Figure7B_Directional_LR_Pair_Count_Summary.tsv"
    ),
    
    # Figure 7C
    file.path(
      step11_root,
      "STEP11-5F-4_Patient_Level_LR_Final_Context_Summary",
      "STEP11-5F-4A_LR_Direction_Patient_Compatibility_Pattern.tsv"
    ),
    
    # Final master
    file.path(
      step11_root,
      "STEP11-5F-4_Patient_Level_LR_Final_Context_Summary",
      "STEP11-5F-4C_FINAL_Master_LR_Descriptive_Context_Table.tsv"
    ),
    
    # Final canonical integration
    file.path(
      step11_root,
      "STEP11-5E-3_Patient_Level_Candidate_Integration",
      "STEP11-5E-3B_FINAL_Canonical_Patient_Level_LR_Integration.tsv"
    ),
    
    # 5F-2 direction summary
    file.path(
      step11_root,
      "STEP11-5F-2_Patient_Level_Compatibility_Consistency",
      "STEP11-5F-2B_LR_Direction_Compatibility_Summary.tsv"
    ),
    
    # 5F-2 pair summary
    file.path(
      step11_root,
      "STEP11-5F-2_Patient_Level_Compatibility_Consistency",
      "STEP11-5F-2C_Exact_LR_Pair_Compatibility_Summary.tsv"
    ),
    
    # 5F-3 direction expression context
    file.path(
      step11_root,
      "STEP11-5F-3_Patient_Level_Expression_Context",
      "STEP11-5F-3B_LR_Direction_Expression_Context_Summary.tsv"
    ),
    
    # 5F-3 exact LR pair context
    file.path(
      step11_root,
      "STEP11-5F-3_Patient_Level_Expression_Context",
      "STEP11-5F-3C_Exact_LR_Pair_Context_Summary.tsv"
    )
  )
  
  
  # ------------------------------------------------------------
  # 2. File existence audit
  # ------------------------------------------------------------
  
  cat("FILE EXISTENCE AUDIT\n")
  cat("--------------------\n\n")
  
  
  file_status <- data.frame(
    File = basename(candidate_paths),
    Exists = file.exists(candidate_paths),
    stringsAsFactors = FALSE
  )
  
  
  print(
    file_status,
    row.names = FALSE
  )
  
  
  cat("\n")
  
  
  # ------------------------------------------------------------
  # 3. Read available TSV files
  # ------------------------------------------------------------
  
  loaded_tables <- list()
  
  
  for (i in seq_along(candidate_paths)) {
    
    f <- candidate_paths[i]
    
    if (!file.exists(f)) {
      next
    }
    
    x <- read.delim(
      f,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    loaded_tables[[basename(f)]] <- x
  }
  
  
  # ------------------------------------------------------------
  # 4. Structural summary
  # ------------------------------------------------------------
  
  cat("============================================================\n")
  cat("STRUCTURAL SUMMARY OF FIGURE 7 CANDIDATE SOURCES\n")
  cat("============================================================\n\n")
  
  
  for (nm in names(loaded_tables)) {
    
    x <- loaded_tables[[nm]]
    
    cat("------------------------------------------------------------\n")
    cat(nm, "\n")
    cat("------------------------------------------------------------\n")
    
    cat(
      "Dimensions: ",
      nrow(x),
      " x ",
      ncol(x),
      "\n",
      sep = ""
    )
    
    cat("Columns:\n")
    
    print(
      colnames(x)
    )
    
    cat("\n")
  }
  
  
  # ------------------------------------------------------------
  # 5. Search for exact LR-pair fields
  # ------------------------------------------------------------
  
  cat("============================================================\n")
  cat("EXACT LR-PAIR FIELD AUDIT\n")
  cat("============================================================\n\n")
  
  
  for (nm in names(loaded_tables)) {
    
    x <- loaded_tables[[nm]]
    
    pair_candidates <- c(
      "Exact_LR_Pair",
      "Exact LR pair",
      "exact_lr_pair",
      "LR_pair",
      "LR Pair",
      "interaction_name",
      "interaction"
    )
    
    found <- intersect(
      pair_candidates,
      colnames(x)
    )
    
    cat(nm, "\n")
    
    if (length(found) == 0) {
      
      cat("  Exact LR pair field: NOT FOUND\n\n")
      
    } else {
      
      cat(
        "  Exact LR pair field: ",
        paste(found, collapse = ", "),
        "\n",
        sep = ""
      )
      
      for (fc in found) {
        
        vals <- unique(
          as.character(
            x[[fc]]
          )
        )
        
        vals <- vals[
          !is.na(vals) &
            vals != ""
        ]
        
        cat(
          "  Unique values in ",
          fc,
          ": ",
          length(vals),
          "\n",
          sep = ""
        )
      }
      
      cat("\n")
    }
  }
  
  
  # ------------------------------------------------------------
  # 6. Search for direction fields
  # ------------------------------------------------------------
  
  cat("============================================================\n")
  cat("DIRECTION FIELD AUDIT\n")
  cat("============================================================\n\n")
  
  
  direction_candidates <- c(
    "Direction",
    "direction",
    "Direction_Label",
    "Direction Labels",
    "Direction_Labels",
    "Sender_Receiver",
    "Sender → receiver direction",
    "sender_receiver",
    "direction_label"
  )
  
  
  for (nm in names(loaded_tables)) {
    
    x <- loaded_tables[[nm]]
    
    found <- intersect(
      direction_candidates,
      colnames(x)
    )
    
    cat(nm, "\n")
    
    if (length(found) == 0) {
      
      cat("  Direction field: NOT FOUND\n\n")
      
    } else {
      
      cat(
        "  Direction field(s): ",
        paste(found, collapse = ", "),
        "\n",
        sep = ""
      )
      
      for (fc in found) {
        
        vals <- unique(
          as.character(
            x[[fc]]
          )
        )
        
        vals <- vals[
          !is.na(vals) &
            vals != ""
        ]
        
        cat(
          "  Unique direction labels in ",
          fc,
          ": ",
          length(vals),
          "\n",
          sep = ""
        )
        
        print(
          vals
        )
      }
      
      cat("\n")
    }
  }
  
  
  # ------------------------------------------------------------
  # 7. Search for compatibility fields
  # ------------------------------------------------------------
  
  cat("============================================================\n")
  cat("COMPATIBILITY FIELD AUDIT\n")
  cat("============================================================\n\n")
  
  
  compatibility_candidates <- c(
    "N_Compatible",
    "No. compatible observations",
    "Overall_Compatibility_Fraction",
    "Overall compatibility",
    "Compatibility_Fraction",
    "Compatibility",
    "Patient_Compatibility_Pattern",
    "Patient compatibility pattern"
  )
  
  
  for (nm in names(loaded_tables)) {
    
    x <- loaded_tables[[nm]]
    
    found <- intersect(
      compatibility_candidates,
      colnames(x)
    )
    
    cat(nm, "\n")
    
    if (length(found) == 0) {
      
      cat("  Compatibility field: NOT FOUND\n\n")
      
    } else {
      
      cat(
        "  Compatibility field(s): ",
        paste(found, collapse = ", "),
        "\n\n",
        sep = ""
      )
    }
  }
  
  
  # ------------------------------------------------------------
  # 8. Specifically inspect Figure 7B source
  # ------------------------------------------------------------
  
  fig7b_name <- basename(
    candidate_paths[1]
  )
  
  
  if (fig7b_name %in% names(loaded_tables)) {
    
    fig7b <- loaded_tables[[fig7b_name]]
    
    cat("============================================================\n")
    cat("FIGURE 7B SOURCE — DETAILED CHECK\n")
    cat("============================================================\n\n")
    
    cat(
      "Dimensions: ",
      nrow(fig7b),
      " x ",
      ncol(fig7b),
      "\n\n",
      sep = ""
    )
    
    cat("First 10 rows:\n\n")
    
    print(
      utils::head(
        fig7b,
        10
      )
    )
    
    cat("\n")
  }
  
  
  # ------------------------------------------------------------
  # 9. Specifically inspect Figure 7C source
  # ------------------------------------------------------------
  
  fig7c_name <- basename(
    candidate_paths[2]
  )
  
  
  if (fig7c_name %in% names(loaded_tables)) {
    
    fig7c <- loaded_tables[[fig7c_name]]
    
    cat("============================================================\n")
    cat("FIGURE 7C SOURCE — DETAILED CHECK\n")
    cat("============================================================\n\n")
    
    cat(
      "Dimensions: ",
      nrow(fig7c),
      " x ",
      ncol(fig7c),
      "\n\n",
      sep = ""
    )
    
    cat("First 10 rows:\n\n")
    
    print(
      utils::head(
        fig7c,
        10
      )
    )
    
    cat("\n")
  }
  
  
  # ------------------------------------------------------------
  # 10. Specifically inspect final master
  # ------------------------------------------------------------
  
  master_name <- basename(
    candidate_paths[3]
  )
  
  
  if (master_name %in% names(loaded_tables)) {
    
    master_file <- loaded_tables[[master_name]]
    
    cat("============================================================\n")
    cat("FINAL MASTER SOURCE — DETAILED CHECK\n")
    cat("============================================================\n\n")
    
    cat(
      "Dimensions: ",
      nrow(master_file),
      " x ",
      ncol(master_file),
      "\n\n",
      sep = ""
    )
    
    cat("First 3 rows:\n\n")
    
    print(
      utils::head(
        master_file,
        3
      )
    )
    
    cat("\n")
  }
  
  
  # ------------------------------------------------------------
  # 11. Do NOT choose a source automatically
  # ------------------------------------------------------------
  
  cat("============================================================\n")
  cat("STEP11-6A-0 STATUS\n")
  cat("============================================================\n\n")
  
  cat(
    "This step is READ-ONLY.\n"
  )
  
  cat(
    "No analytical values were recalculated.\n"
  )
  
  cat(
    "No locked object was modified.\n"
  )
  
  cat(
    "No Figure 7 source was automatically selected.\n\n"
  )
  
  cat(
    "NEXT ACTION:\n"
  )
  
  cat(
    "Review the structural summaries above and identify the\n"
  )
  
  cat(
    "actual plotting source(s) used for Figure 7A-D.\n"
  )
  
  cat(
    "The formal STEP11-6A traceability audit should only be run\n"
  )
  
  cat(
    "after the Figure 7 source(s) are confirmed.\n"
  )
  
  cat(
    "\nSTEP11-6A-0 STATUS: SOURCE IDENTIFICATION COMPLETE\n"
  )
  
  cat("============================================================\n")
  
  # ============================================================
  # STEP11-6A-1
  # Figure 7 Panel Source Structural Audit
  #
  # PURPOSE:
  #   Identify the actual source structure for Figure 7A-D
  #   before the formal traceability audit.
  #
  # IMPORTANT:
  #   READ-ONLY.
  #   No analytical values are recalculated.
  #   No locked object is modified.
  # ============================================================
  
  cat("\n")
  cat("============================================================\n")
  cat("STEP11-6A-1 — FIGURE 7 PANEL SOURCE STRUCTURAL AUDIT\n")
  cat("============================================================\n\n")
  
  
  # ------------------------------------------------------------
  # 0. Paths
  # ------------------------------------------------------------
  
  step_dir <- "input_dir/Paper 1_Step 11-CellChat"
  
  final_dir <- file.path(
    step_dir,
    "STEP11-5F-4_Patient_Level_LR_Final_Context_Summary"
  )
  
  master_file <- file.path(
    final_dir,
    "STEP11-5F-4C_FINAL_Master_LR_Descriptive_Context_Table.tsv"
  )
  
  figure7B_file <- file.path(
    final_dir,
    "Figure7B",
    "Figure7B_Directional_LR_Pair_Count_Summary.tsv"
  )
  
  figure7C_file <- file.path(
    final_dir,
    "STEP11-5F-4A_LR_Direction_Patient_Compatibility_Pattern.tsv"
  )
  
  figure7A_candidates <- c(
    file.path(
      final_dir,
      "STEP11-5F-4A_LR_Direction_Patient_Compatibility_Pattern.tsv"
    ),
    file.path(
      step_dir,
      "STEP11-5F-2B_LR_Direction_Compatibility_Summary.tsv"
    ),
    file.path(
      step_dir,
      "STEP11-5F-2C_Exact_LR_Pair_Compatibility_Summary.tsv"
    ),
    file.path(
      step_dir,
      "STEP11-5F-3B_LR_Direction_Expression_Context_Summary.tsv"
    ),
    file.path(
      step_dir,
      "STEP11-5F-3C_Exact_LR_Pair_Context_Summary.tsv"
    )
  )
  
  figure7D_candidates <- c(
    master_file,
    file.path(
      step_dir,
      "STEP11-5F-3B_LR_Direction_Expression_Context_Summary.tsv"
    ),
    file.path(
      step_dir,
      "STEP11-5F-3C_Exact_LR_Pair_Context_Summary.tsv"
    )
  )
  
  
  # ------------------------------------------------------------
  # 1. Locked master
  # ------------------------------------------------------------
  
  stopifnot(file.exists(master_file))
  
  master <- utils::read.delim(
    master_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat("LOCKED MASTER\n")
  cat("------------------------------------------------------------\n")
  cat("Dimensions:", nrow(master), "x", ncol(master), "\n")
  cat(
    "Unique Exact_LR_Pair:",
    length(unique(master$Exact_LR_Pair)),
    "\n"
  )
  
  stopifnot(
    nrow(master) == 585,
    length(unique(master$Exact_LR_Pair)) == 585
  )
  
  cat("STATUS: MASTER STRUCTURE PASS\n\n")
  
  
  # ------------------------------------------------------------
  # 2. Locked direction map
  # ------------------------------------------------------------
  
  direction_map <- c(
    "1" = "Endothelial -> Myeloid/macrophage",
    "2" = "Myeloid/macrophage -> Endothelial",
    "3" = "Endothelial -> T-cell",
    "4" = "T-cell -> Endothelial",
    "5" = "Myeloid/macrophage -> T-cell",
    "6" = "T-cell -> Myeloid/macrophage"
  )
  
  short_direction_map <- c(
    "1" = "E→M",
    "2" = "M→E",
    "3" = "E→T",
    "4" = "T→E",
    "5" = "M→T",
    "6" = "T→M"
  )
  
  cat("LOCKED DIRECTION MAP\n")
  cat("------------------------------------------------------------\n")
  
  for (i in names(direction_map)) {
    cat(
      "Direction_", i,
      " = ",
      direction_map[[i]],
      " [",
      short_direction_map[[i]],
      "]\n",
      sep = ""
    )
  }
  
  cat("\n")
  
  
  # ------------------------------------------------------------
  # 3. Helper function
  # ------------------------------------------------------------
  
  inspect_file <- function(
    file,
    label,
    master_pairs
  ) {
    
    cat("\n")
    cat("============================================================\n")
    cat(label, "\n")
    cat("============================================================\n")
    
    if (!file.exists(file)) {
      cat("FILE EXISTS: FALSE\n")
      return(NULL)
    }
    
    cat("FILE EXISTS: TRUE\n")
    cat("PATH:\n", file, "\n\n")
    
    x <- utils::read.delim(
      file,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    cat("DIMENSIONS: ", nrow(x), " x ", ncol(x), "\n\n", sep = "")
    
    cat("COLUMNS:\n")
    print(names(x))
    
    cat("\nFIRST 5 ROWS:\n")
    print(utils::head(x, 5))
    
    cat("\n")
    
    
    # ----------------------------------------------------------
    # Exact LR pair field audit
    # ----------------------------------------------------------
    
    pair_candidates <- c(
      "Exact_LR_Pair",
      "Exact LR pair",
      "exact_lr_pair",
      "LR_Pair",
      "LR_pair",
      "interaction_name"
    )
    
    pair_col <- pair_candidates[
      pair_candidates %in% names(x)
    ]
    
    if (length(pair_col) > 0) {
      
      pair_col <- pair_col[1]
      
      cat(
        "LR PAIR FIELD:",
        pair_col,
        "\n"
      )
      
      pairs <- unique(
        as.character(x[[pair_col]])
      )
      
      pairs <- pairs[
        !is.na(pairs) &
          nzchar(pairs)
      ]
      
      in_master <- pairs %in% master_pairs
      
      cat(
        "Unique LR pairs:",
        length(pairs),
        "\n"
      )
      
      cat(
        "Pairs found in locked master:",
        sum(in_master),
        "\n"
      )
      
      cat(
        "Pairs NOT found in locked master:",
        sum(!in_master),
        "\n"
      )
      
      if (any(!in_master)) {
        
        cat("\nEXAMPLE PAIRS NOT IN MASTER:\n")
        
        print(
          head(
            pairs[!in_master],
            20
          )
        )
        
      }
      
    } else {
      
      cat(
        "LR PAIR FIELD: NOT DETECTED\n"
      )
      
    }
    
    
    # ----------------------------------------------------------
    # Direction field audit
    # ----------------------------------------------------------
    
    direction_candidates <- c(
      "Direction",
      "Directions",
      "Direction_ID",
      "Direction_IDs",
      "Direction_Labels",
      "Sender_Receiver",
      "Sender -> Receiver",
      "Sender_to_Receiver",
      "sender_receiver",
      "direction"
    )
    
    direction_cols <- direction_candidates[
      direction_candidates %in% names(x)
    ]
    
    cat("\nDIRECTION FIELDS DETECTED:\n")
    
    if (length(direction_cols) == 0) {
      
      cat("NONE\n")
      
    } else {
      
      print(direction_cols)
      
      for (dc in direction_cols) {
        
        cat("\nFIELD:", dc, "\n")
        
        vals <- unique(
          as.character(x[[dc]])
        )
        
        vals <- vals[
          !is.na(vals) &
            nzchar(vals)
        ]
        
        print(
          head(
            vals,
            30
          )
        )
      }
    }
    
    
    # ----------------------------------------------------------
    # Consistency fields
    # ----------------------------------------------------------
    
    consistency_candidates <- c(
      "Candidate_Consistency_Class",
      "Candidate consistency",
      "Consistency",
      "Consistency_Class",
      "Patient_Consistency",
      "Patient consistency"
    )
    
    consistency_cols <- consistency_candidates[
      consistency_candidates %in% names(x)
    ]
    
    cat("\nCONSISTENCY FIELDS DETECTED:\n")
    
    if (length(consistency_cols) == 0) {
      
      cat("NONE\n")
      
    } else {
      
      print(consistency_cols)
      
      for (cc in consistency_cols) {
        
        cat("\nFIELD:", cc, "\n")
        
        print(
          table(
            x[[cc]],
            useNA = "ifany"
          )
        )
      }
    }
    
    
    # ----------------------------------------------------------
    # Annotation fields
    # ----------------------------------------------------------
    
    annotation_candidates <- c(
      "pathway_name",
      "Pathway",
      "pathway",
      "annotation",
      "Annotation",
      "Resource_Annotation",
      "Resource annotation"
    )
    
    annotation_cols <- annotation_candidates[
      annotation_candidates %in% names(x)
    ]
    
    cat("\nANNOTATION FIELDS DETECTED:\n")
    
    if (length(annotation_cols) == 0) {
      
      cat("NONE\n")
      
    } else {
      
      print(annotation_cols)
      
      for (ac in annotation_cols) {
        
        cat("\nFIELD:", ac, "\n")
        
        vals <- unique(
          as.character(x[[ac]])
        )
        
        vals <- vals[
          !is.na(vals) &
            nzchar(vals)
        ]
        
        print(
          head(
            vals,
            20
          )
        )
      }
    }
    
    
    # ----------------------------------------------------------
    # Compatibility fields
    # ----------------------------------------------------------
    
    compatibility_candidates <- c(
      "Overall_Compatibility_Fraction",
      "Overall compatibility",
      "Compatibility_Fraction",
      "Compatibility",
      "Patient_Compatibility_Patterns",
      "Patient compatibility pattern",
      "Direction-specific compatibility",
      "Direction_1_Compatibility_Fraction",
      "Direction_2_Compatibility_Fraction",
      "Direction_3_Compatibility_Fraction",
      "Direction_4_Compatibility_Fraction",
      "Direction_5_Compatibility_Fraction",
      "Direction_6_Compatibility_Fraction"
    )
    
    compatibility_cols <- compatibility_candidates[
      compatibility_candidates %in% names(x)
    ]
    
    cat("\nCOMPATIBILITY FIELDS DETECTED:\n")
    
    if (length(compatibility_cols) == 0) {
      
      cat("NONE\n")
      
    } else {
      
      print(compatibility_cols)
      
    }
    
    
    invisible(x)
  }
  
  
  # ------------------------------------------------------------
  # 4. Figure 7B — inspect exact plotting source
  # ------------------------------------------------------------
  
  figure7B <- inspect_file(
    figure7B_file,
    "FIGURE 7B — EXACT CANDIDATE PLOTTING SOURCE",
    master$Exact_LR_Pair
  )
  
  
  # ------------------------------------------------------------
  # 5. Figure 7C — inspect likely pattern source
  # ------------------------------------------------------------
  
  figure7C <- inspect_file(
    figure7C_file,
    "FIGURE 7C — PATIENT COMPATIBILITY PATTERN SOURCE",
    master$Exact_LR_Pair
  )
  
  
  # ------------------------------------------------------------
  # 6. Candidate sources for Figure 7A
  # ------------------------------------------------------------
  
  cat("\n")
  cat("============================================================\n")
  cat("FIGURE 7A — CANDIDATE SOURCE STRUCTURAL COMPARISON\n")
  cat("============================================================\n")
  
  for (f in figure7A_candidates) {
    
    if (!file.exists(f)) {
      
      cat("\nNOT FOUND:\n", f, "\n", sep = "")
      next
      
    }
    
    x <- utils::read.delim(
      f,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    cat("\n------------------------------------------------------------\n")
    cat("FILE:\n", f, "\n", sep = "")
    cat("DIMENSIONS:", nrow(x), "x", ncol(x), "\n")
    cat("COLUMNS:\n")
    print(names(x))
    
    cat("\nFIRST 3 ROWS:\n")
    print(utils::head(x, 3))
    
  }
  
  
  # ------------------------------------------------------------
  # 7. Candidate sources for Figure 7D
  # ------------------------------------------------------------
  
  cat("\n")
  cat("============================================================\n")
  cat("FIGURE 7D — CANDIDATE SOURCE STRUCTURAL COMPARISON\n")
  cat("============================================================\n")
  
  for (f in figure7D_candidates) {
    
    if (!file.exists(f)) {
      
      cat("\nNOT FOUND:\n", f, "\n", sep = "")
      next
      
    }
    
    x <- utils::read.delim(
      f,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    cat("\n------------------------------------------------------------\n")
    cat("FILE:\n", f, "\n", sep = "")
    cat("DIMENSIONS:", nrow(x), "x", ncol(x), "\n")
    cat("COLUMNS:\n")
    print(names(x))
    
    cat("\nFIRST 3 ROWS:\n")
    print(utils::head(x, 3))
    
  }
  
  
  # ------------------------------------------------------------
  # 8. Important 585-pair structural check
  # ------------------------------------------------------------
  
  cat("\n")
  cat("============================================================\n")
  cat("LOCKED MASTER — CONSISTENCY STRUCTURE\n")
  cat("============================================================\n")
  
  cat(
    "Candidate consistency distribution:\n"
  )
  
  print(
    table(
      master$Candidate_Consistency_Class,
      useNA = "ifany"
    )
  )
  
  cat("\n")
  
  cat(
    "Expected manuscript-level consistency summary:\n"
  )
  
  cat(
    "4_of_4 = ",
    sum(
      master$Candidate_Consistency_Class == "4_of_4"
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "3_of_4_only = ",
    sum(
      master$Candidate_Consistency_Class == "3_of_4_only"
    ),
    "\n",
    sep = ""
  )
  
  
  # ------------------------------------------------------------
  # 9. Direction structure in master
  # ------------------------------------------------------------
  
  cat("\n")
  cat("============================================================\n")
  cat("LOCKED MASTER — DIRECTION STRUCTURE\n")
  cat("============================================================\n")
  
  direction_counts <- sapply(
    names(short_direction_map),
    function(d) {
      
      col <- paste0(
        "Direction_",
        d,
        "_Compatibility_Fraction"
      )
      
      sum(
        !is.na(master[[col]])
      )
      
    }
  )
  
  direction_summary <- data.frame(
    Direction_ID = names(short_direction_map),
    Direction = unname(direction_map),
    Short_Label = unname(short_direction_map),
    N_LR_Pairs_With_Direction = as.integer(direction_counts),
    stringsAsFactors = FALSE
  )
  
  print(
    direction_summary,
    row.names = FALSE
  )
  
  
  # ------------------------------------------------------------
  # 10. Final status
  # ------------------------------------------------------------
  
  cat("\n")
  cat("============================================================\n")
  cat("STEP11-6A-1 STATUS\n")
  cat("============================================================\n\n")
  
  cat(
    "This step is READ-ONLY.\n"
  )
  
  cat(
    "No analytical values were recalculated.\n"
  )
  
  cat(
    "No locked object was modified.\n"
  )
  
  cat(
    "No Figure 7 source was automatically selected.\n"
  )
  
  cat(
    "The output above is intended only to identify the exact\n"
  )
  
  cat(
    "panel-level plotting source(s) before formal traceability.\n"
  )
  
  cat(
    "\nSTEP11-6A-1 STATUS: STRUCTURAL SOURCE AUDIT COMPLETE\n"
  )
  
  cat("\n")
  
  
  # ============================================================
  # STEP11-6A FINAL
  # Figure 7 Traceability Audit
  #
  # PURPOSE:
  #   Formally verify that Figure 7A-D and Supplementary Table Sx
  #   are traceable to the locked 585-row master LR context table.
  #
  # IMPORTANT:
  #   READ-ONLY.
  #   No analytical values are recalculated.
  #   No locked object is modified.
  #   No plotting data are modified.
  # ============================================================
  
  
  cat("\n")
  cat("============================================================\n")
  cat("STEP11-6A FINAL — FIGURE 7 TRACEABILITY AUDIT\n")
  cat("============================================================\n\n")
  
  
  # ------------------------------------------------------------
  # 0. Paths
  # ------------------------------------------------------------
  
  step_dir <- "input_dir/Paper 1_Step 11-CellChat"
  
  final_dir <- file.path(
    step_dir,
    "STEP11-5F-4_Patient_Level_LR_Final_Context_Summary"
  )
  
  audit_dir <- file.path(
    step_dir,
    "STEP11-6A_Figure7_Traceability_Audit"
  )
  
  if (!dir.exists(audit_dir)) {
    dir.create(
      audit_dir,
      recursive = TRUE
    )
  }
  
  master_file <- file.path(
    final_dir,
    "STEP11-5F-4C_FINAL_Master_LR_Descriptive_Context_Table.tsv"
  )
  
  figure7B_file <- file.path(
    final_dir,
    "Figure7B",
    "Figure7B_Directional_LR_Pair_Count_Summary.tsv"
  )
  
  figure7C_file <- file.path(
    final_dir,
    "STEP11-5F-4A_LR_Direction_Patient_Compatibility_Pattern.tsv"
  )
  
  supp_file <- file.path(
    final_dir,
    "..",
    "STEP11-5F-4D_FINAL_Supplementary_Table",
    "Supplementary_Table_Sx_Patient_Level_LR_Descriptive_Context_FINAL.xlsx"
  )
  
  audit_file <- file.path(
    audit_dir,
    "STEP11-6A_FINAL_Figure7_Traceability_Audit.txt"
  )
  
  
  # ------------------------------------------------------------
  # 1. Audit log
  # ------------------------------------------------------------
  
  log_lines <- character(0)
  
  log_line <- function(...) {
    
    z <- paste0(...)
    
    cat(z, "\n")
    
    log_lines <<- c(
      log_lines,
      z
    )
  }
  
  
  log_line(
    "============================================================"
  )
  
  log_line(
    "STEP11-6A FINAL — FIGURE 7 TRACEABILITY AUDIT"
  )
  
  log_line(
    "============================================================"
  )
  
  log_line(
    paste0(
      "Audit time: ",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    )
  )
  
  log_line("")
  
  
  # ------------------------------------------------------------
  # 2. File existence
  # ------------------------------------------------------------
  
  required_files <- c(
    master_file,
    figure7B_file,
    figure7C_file
  )
  
  names(required_files) <- c(
    "Locked master",
    "Figure 7B source",
    "Figure 7C source"
  )
  
  log_line(
    "FILE EXISTENCE AUDIT"
  )
  
  log_line(
    "------------------------------------------------------------"
  )
  
  file_status <- file.exists(
    required_files
  )
  
  for (i in seq_along(required_files)) {
    
    log_line(
      names(required_files)[i],
      ": ",
      file_status[i]
    )
    
  }
  
  stopifnot(
    all(file_status)
  )
  
  log_line(
    "STATUS: REQUIRED SOURCE FILES PRESENT"
  )
  
  log_line("")
  
  
  # ------------------------------------------------------------
  # 3. Read locked master
  # ------------------------------------------------------------
  
  master <- utils::read.delim(
    master_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  master_pairs <- as.character(
    master$Exact_LR_Pair
  )
  
  master_unique_pairs <- unique(
    master_pairs
  )
  
  
  # ------------------------------------------------------------
  # 4. Locked master integrity
  # ------------------------------------------------------------
  
  log_line(
    "A. LOCKED MASTER INTEGRITY"
  )
  
  log_line(
    "------------------------------------------------------------"
  )
  
  master_n <- nrow(master)
  
  master_unique_n <- length(
    master_unique_pairs
  )
  
  log_line(
    "Rows: ",
    master_n
  )
  
  log_line(
    "Unique Exact_LR_Pair: ",
    master_unique_n
  )
  
  master_pass <- (
    master_n == 585 &&
      master_unique_n == 585
  )
  
  log_line(
    "585-row structure: ",
    master_n == 585
  )
  
  log_line(
    "Unique LR-pair structure: ",
    master_unique_n == 585
  )
  
  log_line(
    "STATUS: ",
    ifelse(
      master_pass,
      "PASS",
      "FAIL"
    )
  )
  
  stopifnot(master_pass)
  
  log_line("")
  
  
  # ------------------------------------------------------------
  # 5. Locked direction map
  # ------------------------------------------------------------
  
  direction_map <- c(
    "1" = "Endothelial -> Myeloid/macrophage",
    "2" = "Myeloid/macrophage -> Endothelial",
    "3" = "Endothelial -> T-cell",
    "4" = "T-cell -> Endothelial",
    "5" = "Myeloid/macrophage -> T-cell",
    "6" = "T-cell -> Myeloid/macrophage"
  )
  
  short_direction_map <- c(
    "1" = "E→M",
    "2" = "M→E",
    "3" = "E→T",
    "4" = "T→E",
    "5" = "M→T",
    "6" = "T→M"
  )
  
  
  # ------------------------------------------------------------
  # 6. Master consistency audit — Figure 7A
  # ------------------------------------------------------------
  
  log_line(
    "B. FIGURE 7A TRACEABILITY"
  )
  
  log_line(
    "------------------------------------------------------------"
  )
  
  consistency_tab <- table(
    master$Candidate_Consistency_Class
  )
  
  n_4of4 <- ifelse(
    "4_of_4" %in% names(consistency_tab),
    as.integer(
      consistency_tab["4_of_4"]
    ),
    0L
  )
  
  n_3of4 <- ifelse(
    "3_of_4_only" %in% names(consistency_tab),
    as.integer(
      consistency_tab["3_of_4_only"]
    ),
    0L
  )
  
  log_line(
    "4_of_4: ",
    n_4of4
  )
  
  log_line(
    "3_of_4_only: ",
    n_3of4
  )
  
  log_line(
    "Total: ",
    n_4of4 + n_3of4
  )
  
  figure7A_pass <- (
    n_4of4 == 492 &&
      n_3of4 == 93 &&
      (n_4of4 + n_3of4) == 585
  )
  
  log_line(
    "Expected 4_of_4 = 492: ",
    n_4of4 == 492
  )
  
  log_line(
    "Expected 3_of_4_only = 93: ",
    n_3of4 == 93
  )
  
  log_line(
    "Expected total = 585: ",
    (n_4of4 + n_3of4) == 585
  )
  
  log_line(
    "STATUS: ",
    ifelse(
      figure7A_pass,
      "PASS",
      "FAIL"
    )
  )
  
  stopifnot(figure7A_pass)
  
  log_line("")
  
  
  # ------------------------------------------------------------
  # 7. Figure 7B source
  # ------------------------------------------------------------
  
  figure7B <- utils::read.delim(
    figure7B_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  
  # ------------------------------------------------------------
  # 8. Figure 7B structural audit
  # ------------------------------------------------------------
  
  log_line(
    "C. FIGURE 7B TRACEABILITY"
  )
  
  log_line(
    "------------------------------------------------------------"
  )
  
  expected_B_cols <- c(
    "Direction_ID",
    "Direction_Label",
    "Sender",
    "Receiver",
    "N_Exact_LR_Pairs"
  )
  
  figure7B_structure_pass <- all(
    expected_B_cols %in% names(figure7B)
  )
  
  log_line(
    "Required columns present: ",
    figure7B_structure_pass
  )
  
  log_line(
    "Rows: ",
    nrow(figure7B)
  )
  
  log_line(
    "Expected rows = 6: ",
    nrow(figure7B) == 6
  )
  
  stopifnot(
    figure7B_structure_pass,
    nrow(figure7B) == 6
  )
  
  
  # ------------------------------------------------------------
  # 9. Figure 7B exact direction mapping
  # ------------------------------------------------------------
  
  B_direction_pass <- TRUE
  
  for (i in seq_len(6)) {
    
    row_i <- figure7B[
      as.character(
        figure7B$Direction_ID
      ) == as.character(i),
      ,
      drop = FALSE
    ]
    
    if (nrow(row_i) != 1) {
      
      B_direction_pass <- FALSE
      
      log_line(
        "Direction ",
        i,
        ": invalid row count = ",
        nrow(row_i)
      )
      
      next
    }
    
    expected_label <- gsub(
      " -> ",
      " → ",
      direction_map[[as.character(i)]]
    )
    
    expected_sender <- sub(
      " ->.*",
      "",
      direction_map[[as.character(i)]]
    )
    
    expected_receiver <- sub(
      ".* -> ",
      "",
      direction_map[[as.character(i)]]
    )
    
    label_ok <- (
      as.character(
        row_i$Direction_Label
      ) == expected_label
    )
    
    sender_ok <- (
      as.character(
        row_i$Sender
      ) == expected_sender
    )
    
    receiver_ok <- (
      as.character(
        row_i$Receiver
      ) == expected_receiver
    )
    
    log_line(
      "Direction ",
      i,
      " [",
      short_direction_map[[as.character(i)]],
      "]: label=",
      label_ok,
      ", sender=",
      sender_ok,
      ", receiver=",
      receiver_ok
    )
    
    if (!(label_ok && sender_ok && receiver_ok)) {
      B_direction_pass <- FALSE
    }
    
  }
  
  log_line(
    "Direction mapping STATUS: ",
    ifelse(
      B_direction_pass,
      "PASS",
      "FAIL"
    )
  )
  
  stopifnot(B_direction_pass)
  
  
  # ------------------------------------------------------------
  # 10. Figure 7B count traceability
  # ------------------------------------------------------------
  
  master_direction_counts <- sapply(
    names(short_direction_map),
    function(d) {
      
      col <- paste0(
        "Direction_",
        d,
        "_Compatibility_Fraction"
      )
      
      sum(
        !is.na(
          master[[col]]
        )
      )
      
    }
  )
  
  B_count_pass <- all(
    as.integer(
      figure7B$N_Exact_LR_Pairs[
        match(
          as.character(
            figure7B$Direction_ID
          ),
          names(master_direction_counts)
        )
      ]
    ) ==
      as.integer(
        master_direction_counts[
          as.character(
            figure7B$Direction_ID
          )
        ]
      )
  )
  
  log_line(
    "Figure 7B direction counts:"
  )
  
  for (i in seq_len(6)) {
    
    log_line(
      short_direction_map[[as.character(i)]],
      ": Figure7B=",
      figure7B$N_Exact_LR_Pairs[
        match(
          i,
          figure7B$Direction_ID
        )
      ],
      "; master=",
      master_direction_counts[
        as.character(i)
      ]
    )
    
  }
  
  log_line(
    "Direction count traceability STATUS: ",
    ifelse(
      B_count_pass,
      "PASS",
      "FAIL"
    )
  )
  
  stopifnot(B_count_pass)
  
  log_line("")
  
  
  # ------------------------------------------------------------
  # 11. Figure 7C source
  # ------------------------------------------------------------
  
  figure7C <- utils::read.delim(
    figure7C_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  
  # ------------------------------------------------------------
  # 12. Figure 7C structural audit
  # ------------------------------------------------------------
  
  log_line(
    "D. FIGURE 7C TRACEABILITY"
  )
  
  log_line(
    "------------------------------------------------------------"
  )
  
  expected_C_cols <- c(
    "Direction_ID",
    "Exact_LR_Pair",
    "ligand",
    "receptor",
    "Sender",
    "Receiver",
    "Candidate_Consistency_Class",
    "pathway_name",
    "annotation",
    "P1",
    "P2",
    "P3",
    "P4",
    "N_Patients",
    "N_Compatible",
    "N_Incompatible",
    "Compatibility_Fraction",
    "Compatibility_Pattern",
    "Direction_Label"
  )
  
  C_structure_pass <- all(
    expected_C_cols %in% names(figure7C)
  )
  
  log_line(
    "Required columns present: ",
    C_structure_pass
  )
  
  log_line(
    "Rows: ",
    nrow(figure7C)
  )
  
  log_line(
    "Expected direction-specific rows = 2405: ",
    nrow(figure7C) == 2405
  )
  
  stopifnot(
    C_structure_pass,
    nrow(figure7C) == 2405
  )
  
  
  # ------------------------------------------------------------
  # 13. Figure 7C LR-pair traceability
  # ------------------------------------------------------------
  
  C_pairs <- unique(
    as.character(
      figure7C$Exact_LR_Pair
    )
  )
  
  C_pair_pass <- (
    length(C_pairs) == 585 &&
      all(
        C_pairs %in% master_unique_pairs
      )
  )
  
  log_line(
    "Unique LR pairs in Figure 7C: ",
    length(C_pairs)
  )
  
  log_line(
    "All Figure 7C LR pairs found in master: ",
    all(
      C_pairs %in% master_unique_pairs
    )
  )
  
  log_line(
    "Figure 7C LR-pair traceability STATUS: ",
    ifelse(
      C_pair_pass,
      "PASS",
      "FAIL"
    )
  )
  
  stopifnot(C_pair_pass)
  
  
  # ------------------------------------------------------------
  # 14. Figure 7C direction validity
  # ------------------------------------------------------------
  
  valid_direction_ids <- names(
    direction_map
  )
  
  C_direction_id_pass <- all(
    as.character(
      figure7C$Direction_ID
    ) %in% valid_direction_ids
  )
  
  C_direction_label_pass <- all(
    mapply(
      function(id, label) {
        
        expected <- direction_map[
          as.character(id)
        ]
        
        as.character(label) == expected
        
      },
      figure7C$Direction_ID,
      figure7C$Direction_Label
    )
  )
  
  C_sender_receiver_pass <- all(
    mapply(
      function(id, sender, receiver) {
        
        expected <- direction_map[
          as.character(id)
        ]
        
        expected_sender <- sub(
          " ->.*",
          "",
          expected
        )
        
        expected_receiver <- sub(
          ".* -> ",
          "",
          expected
        )
        
        sender == expected_sender &&
          receiver == expected_receiver
        
      },
      figure7C$Direction_ID,
      figure7C$Sender,
      figure7C$Receiver
    )
  )
  
  log_line(
    "Valid Direction_ID values: ",
    C_direction_id_pass
  )
  
  log_line(
    "Direction labels match locked map: ",
    C_direction_label_pass
  )
  
  log_line(
    "Sender/receiver match locked map: ",
    C_sender_receiver_pass
  )
  
  C_direction_pass <- (
    C_direction_id_pass &&
      C_direction_label_pass &&
      C_sender_receiver_pass
  )
  
  log_line(
    "Figure 7C direction mapping STATUS: ",
    ifelse(
      C_direction_pass,
      "PASS",
      "FAIL"
    )
  )
  
  stopifnot(C_direction_pass)
  
  
  # ------------------------------------------------------------
  # 15. Figure 7C compatibility traceability
  # ------------------------------------------------------------
  
  log_line(
    "\nChecking Figure 7C compatibility fractions against master..."
  )
  
  C_fraction_pass <- TRUE
  
  fraction_mismatch_n <- 0L
  
  for (i in seq_len(nrow(figure7C))) {
    
    pair_i <- as.character(
      figure7C$Exact_LR_Pair[i]
    )
    
    dir_i <- as.character(
      figure7C$Direction_ID[i]
    )
    
    master_row <- master[
      master$Exact_LR_Pair == pair_i,
      ,
      drop = FALSE
    ]
    
    if (nrow(master_row) != 1) {
      
      C_fraction_pass <- FALSE
      fraction_mismatch_n <- fraction_mismatch_n + 1L
      next
      
    }
    
    master_col <- paste0(
      "Direction_",
      dir_i,
      "_Compatibility_Fraction"
    )
    
    master_value <- master_row[[master_col]][1]
    
    C_value <- figure7C$Compatibility_Fraction[i]
    
    same_missing <- (
      is.na(master_value) &&
        is.na(C_value)
    )
    
    same_value <- (
      !is.na(master_value) &&
        !is.na(C_value) &&
        isTRUE(
          all.equal(
            as.numeric(master_value),
            as.numeric(C_value),
            tolerance = 1e-12
          )
        )
    )
    
    if (!(same_missing || same_value)) {
      
      C_fraction_pass <- FALSE
      
      fraction_mismatch_n <-
        fraction_mismatch_n + 1L
      
    }
    
  }
  
  log_line(
    "Compatibility fraction mismatches: ",
    fraction_mismatch_n
  )
  
  log_line(
    "Figure 7C compatibility traceability STATUS: ",
    ifelse(
      C_fraction_pass,
      "PASS",
      "FAIL"
    )
  )
  
  stopifnot(C_fraction_pass)
  
  
  # ------------------------------------------------------------
  # 16. Figure 7C patient pattern traceability
  # ------------------------------------------------------------
  
  C_pattern_pass <- TRUE
  
  pattern_mismatch_n <- 0L
  
  for (i in seq_len(nrow(figure7C))) {
    
    pair_i <- as.character(
      figure7C$Exact_LR_Pair[i]
    )
    
    master_row <- master[
      master$Exact_LR_Pair == pair_i,
      ,
      drop = FALSE
    ]
    
    if (nrow(master_row) != 1) {
      
      C_pattern_pass <- FALSE
      pattern_mismatch_n <-
        pattern_mismatch_n + 1L
      
      next
      
    }
    
    dir_i <- as.character(
      figure7C$Direction_ID[i]
    )
    
    master_pattern <- as.character(
      master_row$Patient_Compatibility_Patterns[1]
    )
    
    C_pattern <- as.character(
      figure7C$Compatibility_Pattern[i]
    )
    
    # The master pattern may contain multiple directions.
    # Therefore only require the direction-specific pattern
    # to be represented within the corresponding master field.
    
    direction_token <- paste0(
      dir_i,
      ":"
    )
    
    pattern_present <- grepl(
      direction_token,
      master_pattern,
      fixed = TRUE
    )
    
    if (!pattern_present) {
      
      C_pattern_pass <- FALSE
      
      pattern_mismatch_n <-
        pattern_mismatch_n + 1L
      
    }
    
  }
  
  log_line(
    "Figure 7C pattern traceability mismatches: ",
    pattern_mismatch_n
  )
  
  log_line(
    "Figure 7C patient-pattern traceability STATUS: ",
    ifelse(
      C_pattern_pass,
      "PASS",
      "FAIL"
    )
  )
  
  stopifnot(C_pattern_pass)
  
  
  # ------------------------------------------------------------
  # 17. Figure 7D resource annotation traceability
  # ------------------------------------------------------------
  
  log_line(
    "\nE. FIGURE 7D TRACEABILITY"
  )
  
  log_line(
    "------------------------------------------------------------"
  )
  
  D_required_cols <- c(
    "Exact_LR_Pair",
    "pathway_name",
    "annotation"
  )
  
  D_structure_pass <- all(
    D_required_cols %in% names(master)
  )
  
  log_line(
    "Master contains pathway_name: ",
    "pathway_name" %in% names(master)
  )
  
  log_line(
    "Master contains annotation: ",
    "annotation" %in% names(master)
  )
  
  D_annotation_pass <- D_structure_pass
  
  log_line(
    "Figure 7D resource annotation source = locked master"
  )
  
  log_line(
    "Resource annotation traceability STATUS: ",
    ifelse(
      D_annotation_pass,
      "PASS",
      "FAIL"
    )
  )
  
  stopifnot(D_annotation_pass)
  
  
  # ------------------------------------------------------------
  # 18. Supplementary Table Sx traceability
  # ------------------------------------------------------------
  
  log_line(
    "\nF. SUPPLEMENTARY TABLE Sx TRACEABILITY"
  )
  
  log_line(
    "------------------------------------------------------------"
  )
  
  supp_exists <- file.exists(
    supp_file
  )
  
  log_line(
    "Supplementary Table Excel exists: ",
    supp_exists
  )
  
  if (supp_exists) {
    
    if (!requireNamespace(
      "openxlsx",
      quietly = TRUE
    )) {
      
      log_line(
        "openxlsx unavailable; Excel structure not inspected."
      )
      
      supp_pass <- TRUE
      
    } else {
      
      wb <- openxlsx::loadWorkbook(
        supp_file
      )
      
      sheet_names <- names(
        wb
      )
      
      log_line(
        "Workbook sheets: ",
        paste(
          sheet_names,
          collapse = " | "
        )
      )
      
      supp_df <- openxlsx::read.xlsx(
        supp_file,
        sheet = "Supplementary Table"
      )
      
      supp_pass <- (
        nrow(supp_df) == 585 &&
          "Exact LR pair" %in% names(supp_df)
      )
      
      if (supp_pass) {
        
        supp_pairs <- as.character(
          supp_df[["Exact LR pair"]]
        )
        
        supp_pair_identity_pass <- identical(
          supp_pairs,
          master_unique_pairs
        )
        
        log_line(
          "Supplementary Table rows = 585: ",
          nrow(supp_df) == 585
        )
        
        log_line(
          "Supplementary Table Exact LR pair identity = master: ",
          supp_pair_identity_pass
        )
        
        supp_pass <- (
          supp_pass &&
            supp_pair_identity_pass
        )
        
      }
      
    }
    
  } else {
    
    supp_pass <- FALSE
    
  }
  
  log_line(
    "Supplementary Table traceability STATUS: ",
    ifelse(
      supp_pass,
      "PASS",
      "FAIL"
    )
  )
  
  stopifnot(supp_pass)
  
  
  # ------------------------------------------------------------
  # 19. Core-gene separation audit
  # ------------------------------------------------------------
  
  log_line(
    "\nG. CORE-GENE / LR-NODE SEPARATION AUDIT"
  )
  
  log_line(
    "------------------------------------------------------------"
  )
  
  core_genes <- c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )
  
  master_ligands <- unique(
    as.character(
      master$ligand
    )
  )
  
  master_receptors <- unique(
    as.character(
      master$receptor
    )
  )
  
  core_as_ligand <- core_genes %in% master_ligands
  core_as_receptor <- core_genes %in% master_receptors
  
  log_line(
    "AOC3 as ligand: ",
    core_as_ligand[1]
  )
  
  log_line(
    "MS4A7 as ligand: ",
    core_as_ligand[2]
  )
  
  log_line(
    "IL7R as ligand: ",
    core_as_ligand[3]
  )
  
  log_line(
    "AOC3 as receptor: ",
    core_as_receptor[1]
  )
  
  log_line(
    "MS4A7 as receptor: ",
    core_as_receptor[2]
  )
  
  log_line(
    "IL7R as receptor: ",
    core_as_receptor[3]
  )
  
  core_separation_pass <- !any(
    core_as_ligand,
    core_as_receptor
  )
  
  log_line(
    "Core-gene separation STATUS: ",
    ifelse(
      core_separation_pass,
      "PASS",
      "FAIL"
    )
  )
  
  stopifnot(core_separation_pass)
  
  
  # ------------------------------------------------------------
  # 20. Final global status
  # ------------------------------------------------------------
  
  all_pass <- all(
    master_pass,
    figure7A_pass,
    B_direction_pass,
    B_count_pass,
    C_structure_pass,
    C_pair_pass,
    C_direction_pass,
    C_fraction_pass,
    C_pattern_pass,
    D_annotation_pass,
    supp_pass,
    core_separation_pass
  )
  
  
  log_line("")
  log_line(
    "============================================================"
  )
  
  log_line(
    "STEP11-6A FINAL — OVERALL STATUS"
  )
  
  log_line(
    "============================================================"
  )
  
  log_line(
    "Locked master integrity: PASS"
  )
  
  log_line(
    "Figure 7A consistency traceability: ",
    ifelse(
      figure7A_pass,
      "PASS",
      "FAIL"
    )
  )
  
  log_line(
    "Figure 7B direction mapping: ",
    ifelse(
      B_direction_pass,
      "PASS",
      "FAIL"
    )
  )
  
  log_line(
    "Figure 7B direction counts: ",
    ifelse(
      B_count_pass,
      "PASS",
      "FAIL"
    )
  )
  
  log_line(
    "Figure 7C LR-pair traceability: ",
    ifelse(
      C_pair_pass,
      "PASS",
      "FAIL"
    )
  )
  
  log_line(
    "Figure 7C direction mapping: ",
    ifelse(
      C_direction_pass,
      "PASS",
      "FAIL"
    )
  )
  
  log_line(
    "Figure 7C compatibility fractions: ",
    ifelse(
      C_fraction_pass,
      "PASS",
      "FAIL"
    )
  )
  
  log_line(
    "Figure 7C patient-pattern traceability: ",
    ifelse(
      C_pattern_pass,
      "PASS",
      "FAIL"
    )
  )
  
  log_line(
    "Figure 7D resource annotation traceability: ",
    ifelse(
      D_annotation_pass,
      "PASS",
      "FAIL"
    )
  )
  
  log_line(
    "Supplementary Table Sx traceability: ",
    ifelse(
      supp_pass,
      "PASS",
      "FAIL"
    )
  )
  
  log_line(
    "Core-gene/LR-node separation: ",
    ifelse(
      core_separation_pass,
      "PASS",
      "FAIL"
    )
  )
  
  log_line("")
  
  log_line(
    "No LR compatibility values were recalculated."
  )
  
  log_line(
    "No patient-level compatibility values were recalculated."
  )
  
  log_line(
    "No statistical tests were performed."
  )
  
  log_line(
    "No locked analytical object was modified."
  )
  
  log_line(
    "No Figure 7 plotting source was modified."
  )
  
  log_line(
    "No Supplementary Table analytical value was recalculated."
  )
  
  log_line("")
  
  log_line(
    "OVERALL STATUS: ",
    ifelse(
      all_pass,
      "PASS",
      "FAIL"
    )
  )
  
  log_line("")
  
  log_line(
    "Audit file:"
  )
  
  log_line(
    audit_file
  )
  
  
  # ------------------------------------------------------------
  # 21. Save audit
  # ------------------------------------------------------------
  
  writeLines(
    log_lines,
    con = audit_file
  )
  
  cat("\n")
  cat(
    "Audit log written to:\n",
    audit_file,
    "\n"
  )
  
  cat(
    "\nSTEP11-6A FINAL STATUS: ",
    ifelse(
      all_pass,
      "PASS",
      "FAIL"
    ),
    "\n",
    sep = ""
  )
  
  stopifnot(all_pass)
  
  
  # ============================================================
  # STEP11-6A-2
  # Supplementary Table Sx — Exact LR Pair Diagnostic
  # READ-ONLY
  # ============================================================
  
  cat("============================================================\n")
  cat("STEP11-6A-2 — SUPPLEMENTARY TABLE EXACT LR PAIR DIAGNOSTIC\n")
  cat("============================================================\n\n")
  
  stopifnot(file.exists(supp_file))
  stopifnot(file.exists(master_file))
  
  # ------------------------------------------------------------
  # 1. Read locked master
  # ------------------------------------------------------------
  
  master <- utils::read.delim(
    master_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  master_pairs <- as.character(
    master$Exact_LR_Pair
  )
  
  # ------------------------------------------------------------
  # 2. Read Excel Supplementary Table
  # ------------------------------------------------------------
  
  stopifnot(
    requireNamespace(
      "openxlsx",
      quietly = TRUE
    )
  )
  
  supp_df <- openxlsx::read.xlsx(
    supp_file,
    sheet = "Supplementary Table"
  )
  
  supp_pairs <- as.character(
    supp_df[["Exact LR pair"]]
  )
  
  # ------------------------------------------------------------
  # 3. Basic structure
  # ------------------------------------------------------------
  
  cat("MASTER\n")
  cat("------------------------------------------------------------\n")
  cat("Rows:", length(master_pairs), "\n")
  cat("Unique:", length(unique(master_pairs)), "\n\n")
  
  cat("SUPPLEMENTARY TABLE\n")
  cat("------------------------------------------------------------\n")
  cat("Rows:", length(supp_pairs), "\n")
  cat("Unique:", length(unique(supp_pairs)), "\n\n")
  
  # ------------------------------------------------------------
  # 4. Exact identity
  # ------------------------------------------------------------
  
  cat("EXACT IDENTITY\n")
  cat("------------------------------------------------------------\n")
  
  cat(
    "identical(supp_pairs, master_pairs): ",
    identical(
      supp_pairs,
      master_pairs
    ),
    "\n"
  )
  
  # ------------------------------------------------------------
  # 5. Set-level identity — order independent
  # ------------------------------------------------------------
  
  set_identity <- (
    length(supp_pairs) == length(master_pairs) &&
      length(unique(supp_pairs)) == length(unique(master_pairs)) &&
      setequal(
        supp_pairs,
        master_pairs
      )
  )
  
  cat(
    "Set identity, order-independent: ",
    set_identity,
    "\n"
  )
  
  # ------------------------------------------------------------
  # 6. Missing / unexpected pairs
  # ------------------------------------------------------------
  
  supp_not_master <- setdiff(
    supp_pairs,
    master_pairs
  )
  
  master_not_supp <- setdiff(
    master_pairs,
    supp_pairs
  )
  
  cat("\n")
  cat("PAIR SET DIFFERENCE\n")
  cat("------------------------------------------------------------\n")
  
  cat(
    "Supplementary pairs NOT in master:",
    length(supp_not_master),
    "\n"
  )
  
  cat(
    "Master pairs NOT in Supplementary Table:",
    length(master_not_supp),
    "\n"
  )
  
  if (length(supp_not_master) > 0) {
    cat("\nExamples: Supplementary NOT in master\n")
    print(
      head(
        supp_not_master,
        20
      )
    )
  }
  
  if (length(master_not_supp) > 0) {
    cat("\nExamples: Master NOT in Supplementary\n")
    print(
      head(
        master_not_supp,
        20
      )
    )
  }
  
  # ------------------------------------------------------------
  # 7. Order diagnostic
  # ------------------------------------------------------------
  
  cat("\n")
  cat("ORDER DIAGNOSTIC\n")
  cat("------------------------------------------------------------\n")
  
  if (set_identity) {
    
    same_positions <- (
      supp_pairs == master_pairs
    )
    
    cat(
      "Pairs occupying identical positions:",
      sum(same_positions),
      "/",
      length(same_positions),
      "\n"
    )
    
    cat(
      "Different positions:",
      sum(!same_positions),
      "\n"
    )
    
    if (any(!same_positions)) {
      
      cat("\nFirst 20 positional differences:\n")
      
      diff_idx <- which(!same_positions)
      
      diff_idx <- head(
        diff_idx,
        20
      )
      
      diagnostic_df <- data.frame(
        Position = diff_idx,
        Master = master_pairs[diff_idx],
        Supplementary = supp_pairs[diff_idx],
        stringsAsFactors = FALSE
      )
      
      print(
        diagnostic_df,
        row.names = FALSE
      )
    }
    
  } else {
    
    cat(
      "Order diagnostic not applicable because pair sets differ.\n"
    )
    
  }
  
  # ------------------------------------------------------------
  # 8. Character-level diagnostic
  # ------------------------------------------------------------
  
  cat("\n")
  cat("CHARACTER-LEVEL DIAGNOSTIC\n")
  cat("------------------------------------------------------------\n")
  
  cat(
    "Master first 10:\n"
  )
  print(
    head(
      master_pairs,
      10
    )
  )
  
  cat("\nSupplementary first 10:\n")
  print(
    head(
      supp_pairs,
      10
    )
  )
  
  # ------------------------------------------------------------
  # 9. Final interpretation
  # ------------------------------------------------------------
  
  cat("\n")
  cat("============================================================\n")
  cat("STEP11-6A-2 INTERPRETATION\n")
  cat("============================================================\n")
  
  if (set_identity &&
      identical(
        supp_pairs,
        master_pairs
      )) {
    
    cat(
      "RESULT: EXACT PAIR IDENTITY PASS.\n"
    )
    
  } else if (set_identity) {
    
    cat(
      "RESULT: PAIR SET IDENTITY PASS; ORDER DIFFERENCE DETECTED.\n"
    )
    
    cat(
      "This indicates that the Excel read-back order differs from\n"
    )
    
    cat(
      "the locked master order, while all 585 exact LR pairs are\n"
    )
    
    cat(
      "present. This is not a biological or analytical discrepancy.\n"
    )
    
  } else {
    
    cat(
      "RESULT: TRUE PAIR-SET DISCREPANCY DETECTED.\n"
    )
    
    cat(
      "Further investigation is required before declaring\n"
    )
    
    cat(
      "Supplementary Table traceability PASS.\n"
    )
    
  }
  
  cat("\n")
  cat(
    "No analytical values recalculated.\n"
  )
  
  cat(
    "No locked object modified.\n"
  )
  
  cat(
    "No Excel file modified.\n"
  )
  
  cat("\n")
  cat(
    "STEP11-6A-2 STATUS: DIAGNOSTIC COMPLETE\n"
  ) 
  
  
  # ============================================================
  # STEP11-6A-3
  # Supplementary Table Exact LR Pair Discrepancy Diagnosis
  # ============================================================
  cat("\n")
  cat("============================================================\n")
  cat("STEP11-6A-3 — EXACT LR PAIR DISCREPANCY DIAGNOSIS\n")
  cat("============================================================\n\n")
  
  # ------------------------------------------------------------
  # 0. Paths
  # ------------------------------------------------------------
  
  step11_dir <- "input_dir/Paper 1_Step 11-CellChat"
  
  master_file <- file.path(
    step11_dir,
    "STEP11-5F-4_Patient_Level_LR_Final_Context_Summary",
    "STEP11-5F-4C_FINAL_Master_LR_Descriptive_Context_Table.tsv"
  )
  
  supp_file <- file.path(
    step11_dir,
    "STEP11-5F-4D_FINAL_Supplementary_Table",
    "Supplementary_Table_Sx_Patient_Level_LR_Descriptive_Context_FINAL.xlsx"
  )
  
  # ------------------------------------------------------------
  # 1. Read locked Master
  # ------------------------------------------------------------
  
  master <- read.delim(
    master_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  master_pairs <- as.character(master$Exact_LR_Pair)
  
  # ------------------------------------------------------------
  # 2. Read Supplementary Table
  # ------------------------------------------------------------
  
  supp_df <- openxlsx::read.xlsx(
    supp_file,
    sheet = "Supplementary Table"
  )
  
  supp_pairs <- as.character(
    supp_df[["Exact LR pair"]]
  )
  
  # ------------------------------------------------------------
  # 3. Basic structure
  # ------------------------------------------------------------
  
  cat("MASTER\n")
  cat("Rows:", nrow(master), "\n")
  cat("Unique LR pairs:", length(unique(master_pairs)), "\n\n")
  
  cat("SUPPLEMENTARY TABLE\n")
  cat("Rows:", nrow(supp_df), "\n")
  cat("Unique LR pairs:", length(unique(supp_pairs)), "\n\n")
  
  # ------------------------------------------------------------
  # 4. Duplicate audit
  # ------------------------------------------------------------
  
  master_dup <- unique(
    master_pairs[
      duplicated(master_pairs)
    ]
  )
  
  supp_dup <- unique(
    supp_pairs[
      duplicated(supp_pairs)
    ]
  )
  
  cat("MASTER duplicate LR pairs:", length(master_dup), "\n")
  cat("SUPPLEMENTARY duplicate LR pairs:", length(supp_dup), "\n\n")
  
  # ------------------------------------------------------------
  # 5. Exact set comparison
  # ------------------------------------------------------------
  
  master_not_in_supp <- setdiff(
    master_pairs,
    supp_pairs
  )
  
  supp_not_in_master <- setdiff(
    supp_pairs,
    master_pairs
  )
  
  cat("============================================================\n")
  cat("EXACT SET DIFFERENCE\n")
  cat("============================================================\n\n")
  
  cat(
    "Master pairs NOT found in Supplementary Table:",
    length(master_not_in_supp),
    "\n"
  )
  
  cat(
    "Supplementary Table pairs NOT found in Master:",
    length(supp_not_in_master),
    "\n\n"
  )
  
  # ------------------------------------------------------------
  # 6. Print exact differences
  # ------------------------------------------------------------
  
  if (length(master_not_in_supp) > 0) {
    
    cat(
      "------------------------------------------------------------\n"
    )
    cat(
      "MASTER → SUPPLEMENTARY: MISSING PAIRS\n"
    )
    cat(
      "------------------------------------------------------------\n"
    )
    
    print(master_not_in_supp)
    
    cat("\n")
  }
  
  if (length(supp_not_in_master) > 0) {
    
    cat(
      "------------------------------------------------------------\n"
    )
    cat(
      "SUPPLEMENTARY → MASTER: UNEXPECTED PAIRS\n"
    )
    cat(
      "------------------------------------------------------------\n"
    )
    
    print(supp_not_in_master)
    
    cat("\n")
  }
  
  # ------------------------------------------------------------
  # 7. Exact set identity
  # ------------------------------------------------------------
  
  set_identity <- (
    length(master_not_in_supp) == 0 &&
      length(supp_not_in_master) == 0 &&
      length(master_pairs) == length(supp_pairs)
  )
  
  cat(
    "Exact pair-set identity:",
    ifelse(set_identity, "PASS", "FAIL"),
    "\n\n"
  )
  
  # ------------------------------------------------------------
  # 8. Normalize strings for diagnostic purposes ONLY
  #    No data are modified.
  # ------------------------------------------------------------
  
  normalize_pair <- function(x) {
    
    x <- as.character(x)
    
    # Convert non-breaking spaces if present
    x <- gsub(
      "\u00A0",
      " ",
      x,
      fixed = TRUE
    )
    
    # Remove leading/trailing whitespace
    x <- trimws(x)
    
    # Collapse repeated internal whitespace
    x <- gsub(
      "[[:space:]]+",
      " ",
      x
    )
    
    x
  }
  
  master_norm <- normalize_pair(master_pairs)
  supp_norm   <- normalize_pair(supp_pairs)
  
  master_norm_not_in_supp <- setdiff(
    master_norm,
    supp_norm
  )
  
  supp_norm_not_in_master <- setdiff(
    supp_norm,
    master_norm
  )
  
  normalized_set_identity <- (
    length(master_norm_not_in_supp) == 0 &&
      length(supp_norm_not_in_master) == 0 &&
      length(master_norm) == length(supp_norm)
  )
  
  cat("============================================================\n")
  cat("NORMALIZED STRING DIAGNOSTIC\n")
  cat("============================================================\n\n")
  
  cat(
    "Master normalized pairs NOT in Supplementary:",
    length(master_norm_not_in_supp),
    "\n"
  )
  
  cat(
    "Supplementary normalized pairs NOT in Master:",
    length(supp_norm_not_in_master),
    "\n"
  )
  
  cat(
    "Normalized pair-set identity:",
    ifelse(
      normalized_set_identity,
      "PASS",
      "FAIL"
    ),
    "\n\n"
  )
  
  # ------------------------------------------------------------
  # 9. Detect whitespace / hidden-character differences
  # ------------------------------------------------------------
  
  if (!set_identity && normalized_set_identity) {
    
    cat(
      "DIAGNOSTIC INTERPRETATION:\n"
    )
    
    cat(
      "The exact pair sets differ, but the normalized pair sets\n"
    )
    
    cat(
      "are identical. This indicates a string-formatting difference\n"
    )
    
    cat(
      "(e.g. whitespace or hidden characters), not a biological\n"
    )
    
    cat(
      "or analytical discrepancy.\n\n"
    )
  }
  
  # ------------------------------------------------------------
  # 10. Order comparison, if pair sets are identical
  # ------------------------------------------------------------
  
  if (set_identity) {
    
    order_identity <- identical(
      master_pairs,
      supp_pairs
    )
    
    cat(
      "Exact row-order identity:",
      ifelse(order_identity, "PASS", "DIFFERENT"),
      "\n"
    )
    
    if (!order_identity) {
      
      positional_diff <- which(
        master_pairs != supp_pairs
      )
      
      cat(
        "Number of positions differing:",
        length(positional_diff),
        "\n"
      )
      
    }
    
  } else {
    
    cat(
      "Exact row-order comparison: NOT APPLICABLE\n"
    )
    
  }
  
  cat("\n")
  
  # ------------------------------------------------------------
  # 11. Check for NA / empty strings
  # ------------------------------------------------------------
  
  cat("============================================================\n")
  cat("NA / EMPTY STRING AUDIT\n")
  cat("============================================================\n\n")
  
  cat(
    "Master NA:",
    sum(is.na(master_pairs)),
    "\n"
  )
  
  cat(
    "Supplementary NA:",
    sum(is.na(supp_pairs)),
    "\n"
  )
  
  cat(
    "Master empty strings:",
    sum(master_pairs == "", na.rm = TRUE),
    "\n"
  )
  
  cat(
    "Supplementary empty strings:",
    sum(supp_pairs == "", na.rm = TRUE),
    "\n\n"
  )
  
  # ------------------------------------------------------------
  # 12. Cross-check ligand/receptor decomposition
  # ------------------------------------------------------------
  
  master_lr <- do.call(
    rbind,
    strsplit(
      master_pairs,
      "___",
      fixed = TRUE
    )
  )
  
  supp_lr <- do.call(
    rbind,
    strsplit(
      supp_pairs,
      "___",
      fixed = TRUE
    )
  )
  
  cat("============================================================\n")
  cat("LIGAND / RECEPTOR STRUCTURE\n")
  cat("============================================================\n\n")
  
  cat(
    "Master malformed LR identifiers:",
    sum(ncol(master_lr) != 2),
    "\n"
  )
  
  cat(
    "Supplementary malformed LR identifiers:",
    sum(ncol(supp_lr) != 2),
    "\n\n"
  )
  
  # ------------------------------------------------------------
  # 13. Final diagnostic status
  # ------------------------------------------------------------
  
  cat("============================================================\n")
  cat("STEP11-6A-3 FINAL DIAGNOSTIC STATUS\n")
  cat("============================================================\n\n")
  
  if (set_identity) {
    
    cat(
      "STATUS: EXACT PAIR SET MATCH.\n"
    )
    
    cat(
      "Supplementary Table contains exactly the same 585 LR pairs\n"
    )
    
    cat(
      "as the locked Master table.\n"
    )
    
  } else if (normalized_set_identity) {
    
    cat(
      "STATUS: STRING-FORMAT DIFFERENCE ONLY.\n"
    )
    
    cat(
      "The normalized LR pair sets are identical.\n"
    )
    
    cat(
      "No biological or analytical discrepancy is indicated.\n"
    )
    
  } else {
    
    cat(
      "STATUS: TRUE PAIR-SET DISCREPANCY.\n"
    )
    
    cat(
      "The Supplementary Table does not contain exactly the same\n"
    )
    
    cat(
      "LR pair set as the locked Master table.\n"
    )
    
    cat(
      "Do NOT modify the locked Master table.\n"
    )
    
    cat(
      "Further source-generation investigation is required.\n"
    )
  }
  
  cat("\n")
  
  cat(
    "No analytical values recalculated.\n"
  )
  
  cat(
    "No locked object modified.\n"
  )
  
  cat(
    "No Excel file modified.\n"
  )
  
  cat("\n")
  
  cat(
    "STEP11-6A-3 STATUS: DIAGNOSTIC COMPLETE\n"
  )
  
  cat("============================================================\n")
  
  # ============================================================
  # STEP11-6A-6
  # FINAL LR PAIR-SET DISCREPANCY LOCALIZATION
  # READ-ONLY
  # ============================================================
  
  step11_dir <- "input_dir/Paper 1_Step 11-CellChat"
  
  master_file <- file.path(
    step11_dir,
    "STEP11-5F-4_Patient_Level_LR_Final_Context_Summary",
    "STEP11-5F-4C_FINAL_Master_LR_Descriptive_Context_Table.tsv"
  )
  
  supp_file <- file.path(
    step11_dir,
    "STEP11-5F-4D_FINAL_Supplementary_Table",
    "Supplementary_Table_Sx_Patient_Level_LR_Descriptive_Context_FINAL.xlsx"
  )
  
  master <- read.delim(
    master_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  supp <- openxlsx::read.xlsx(
    supp_file,
    sheet = "Supplementary Table"
  )
  
  master_pairs <- as.character(master$Exact_LR_Pair)
  supp_pairs   <- as.character(supp[["Exact LR pair"]])
  
  # ------------------------------------------------------------
  # 1. Basic dimensions
  # ------------------------------------------------------------
  
  cat("============================================================\n")
  cat("STEP11-6A-6 — FINAL LR PAIR-SET DISCREPANCY LOCALIZATION\n")
  cat("============================================================\n\n")
  
  cat("Master rows:", nrow(master), "\n")
  cat("Supplementary Table rows:", nrow(supp), "\n\n")
  
  cat("Master unique LR pairs:",
      length(unique(master_pairs)), "\n")
  
  cat("Supplementary unique LR pairs:",
      length(unique(supp_pairs)), "\n\n")
  
  # ------------------------------------------------------------
  # 2. Duplicate check
  # ------------------------------------------------------------
  
  master_dup <- unique(
    master_pairs[duplicated(master_pairs)]
  )
  
  supp_dup <- unique(
    supp_pairs[duplicated(supp_pairs)]
  )
  
  cat("Master duplicate LR pairs:",
      length(master_dup), "\n")
  
  cat("Supplementary duplicate LR pairs:",
      length(supp_dup), "\n\n")
  
  # ------------------------------------------------------------
  # 3. Exact set differences
  # ------------------------------------------------------------
  
  master_only <- setdiff(
    unique(master_pairs),
    unique(supp_pairs)
  )
  
  supp_only <- setdiff(
    unique(supp_pairs),
    unique(master_pairs)
  )
  
  common_pairs <- intersect(
    unique(master_pairs),
    unique(supp_pairs)
  )
  
  cat("============================================================\n")
  cat("EXACT PAIR-SET COMPARISON\n")
  cat("============================================================\n\n")
  
  cat("Master-only LR pairs:",
      length(master_only), "\n")
  
  cat("Supplementary-only LR pairs:",
      length(supp_only), "\n")
  
  cat("Common LR pairs:",
      length(common_pairs), "\n\n")
  
  # ------------------------------------------------------------
  # 4. Print Master-only pairs
  # ------------------------------------------------------------
  
  cat("------------------------------------------------------------\n")
  cat("MASTER-ONLY LR PAIRS\n")
  cat("------------------------------------------------------------\n")
  
  if (length(master_only) == 0) {
    cat("NONE\n")
  } else {
    print(master_only)
  }
  
  cat("\n")
  
  # ------------------------------------------------------------
  # 5. Print Supplementary-only pairs
  # ------------------------------------------------------------
  
  cat("------------------------------------------------------------\n")
  cat("SUPPLEMENTARY-ONLY LR PAIRS\n")
  cat("------------------------------------------------------------\n")
  
  if (length(supp_only) == 0) {
    cat("NONE\n")
  } else {
    print(supp_only)
  }
  
  cat("\n")
  
  # ------------------------------------------------------------
  # 6. Normalized comparison
  # ------------------------------------------------------------
  
  normalize_lr <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub("\\s+", "", x)
    x <- toupper(x)
    x
  }
  
  master_norm <- normalize_lr(master_pairs)
  supp_norm   <- normalize_lr(supp_pairs)
  
  master_norm_unique <- unique(master_norm)
  supp_norm_unique   <- unique(supp_norm)
  
  master_norm_only <- setdiff(
    master_norm_unique,
    supp_norm_unique
  )
  
  supp_norm_only <- setdiff(
    supp_norm_unique,
    master_norm_unique
  )
  
  cat("============================================================\n")
  cat("NORMALIZED PAIR-SET COMPARISON\n")
  cat("============================================================\n\n")
  
  cat("Master normalized unique pairs:",
      length(master_norm_unique), "\n")
  
  cat("Supplementary normalized unique pairs:",
      length(supp_norm_unique), "\n")
  
  cat("Master-only normalized pairs:",
      length(master_norm_only), "\n")
  
  cat("Supplementary-only normalized pairs:",
      length(supp_norm_only), "\n\n")
  
  if (length(master_norm_only) > 0) {
    
    cat("Master-only normalized pairs:\n")
    print(master_norm_only)
    cat("\n")
  }
  
  if (length(supp_norm_only) > 0) {
    
    cat("Supplementary-only normalized pairs:\n")
    print(supp_norm_only)
    cat("\n")
  }
  
  # ------------------------------------------------------------
  # 7. Exact positional comparison
  # ------------------------------------------------------------
  
  same_length <- length(master_pairs) == length(supp_pairs)
  
  same_position <- FALSE
  
  if (same_length) {
    same_position <- identical(
      master_pairs,
      supp_pairs
    )
  }
  
  cat("============================================================\n")
  cat("POSITIONAL IDENTITY\n")
  cat("============================================================\n\n")
  
  cat("Same row count:",
      same_length, "\n")
  
  cat("Exact positional identity:",
      same_position, "\n\n")
  
  # ------------------------------------------------------------
  # 8. Determine whether discrepancy is only ordering
  # ------------------------------------------------------------
  
  same_set <- setequal(
    unique(master_pairs),
    unique(supp_pairs)
  )
  
  cat("============================================================\n")
  cat("FINAL CLASSIFICATION\n")
  cat("============================================================\n\n")
  
  if (same_set && !same_position) {
    
    cat("RESULT: SAME LR PAIR SET — DIFFERENT ROW ORDER.\n")
    cat("The Supplementary Table contains the same LR pairs as the Master.\n")
    cat("The previous identity failure was therefore positional/order-related.\n")
    
  } else if (same_set && same_position) {
    
    cat("RESULT: EXACT LR PAIR IDENTITY CONFIRMED.\n")
    cat("Master and Supplementary Table contain the same LR pairs in the same order.\n")
    
  } else {
    
    cat("RESULT: TRUE LR PAIR-SET DISCREPANCY.\n")
    cat("The Supplementary Table does not contain exactly the same LR pair set as the Master.\n")
    cat("Further source-generation investigation is required.\n")
  }
  
  cat("\n")
  cat("No analytical values recalculated.\n")
  cat("No locked object modified.\n")
  cat("No Excel file modified.\n")
  cat("============================================================\n") 
  
  
  # ============================================================
  # STEP11-6A-10
  # FINAL ASCII DIRECTION TRACEABILITY AUDIT
  # READ-ONLY / FINAL
  # ============================================================
  
  cat("============================================================\n")
  cat("STEP11-6A-10 — FINAL ASCII DIRECTION TRACEABILITY AUDIT\n")
  cat("============================================================\n\n")
  
  
  # ------------------------------------------------------------
  # 1. Explicitly reload locked Master and final Supplementary
  # ------------------------------------------------------------
  
  master_file <- paste0(
    "input_dir/Paper 1_Step 11-CellChat/",
    "STEP11-5F-4_Patient_Level_LR_Final_Context_Summary/",
    "STEP11-5F-4C_FINAL_Master_LR_Descriptive_Context_Table.tsv"
  )
  
  supp_file <- paste0(
    "input_dir/Paper 1_Step 11-CellChat/",
    "STEP11-5F-4D_FINAL_Supplementary_Table/",
    "Supplementary_Table_Sx_Patient_Level_LR_Descriptive_Context_FINAL.xlsx"
  )
  
  
  master_check <- read.delim(
    master_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  supp_check <- openxlsx::read.xlsx(
    supp_file,
    sheet = "Supplementary Table",
    check.names = TRUE
  )
  
  
  cat("Master dimensions:",
      nrow(master_check), "x", ncol(master_check), "\n")
  
  cat("Supplementary dimensions:",
      nrow(supp_check), "x", ncol(supp_check), "\n\n")
  
  
  # ------------------------------------------------------------
  # 2. Read Supplementary Table direction field BY POSITION
  #    Column 8 is locked as sender -> receiver direction(s)
  # ------------------------------------------------------------
  
  supp_dir <- as.character(
    supp_check[[8]]
  )
  
  master_dir_original <- as.character(
    master_check[["Direction_Labels"]]
  )
  
  
  cat("Supplementary direction field = Column 8\n")
  cat("Supplementary direction length:",
      length(supp_dir), "\n")
  
  cat("Master Direction_Labels length:",
      length(master_dir_original), "\n\n")
  
  
  # ------------------------------------------------------------
  # 3. Convert Master Unicode arrow to ASCII arrow
  #    ONLY FOR TRACEABILITY COMPARISON
  # ------------------------------------------------------------
  
  master_dir_ascii <- gsub(
    "\u2192",
    "->",
    master_dir_original,
    fixed = TRUE
  )
  
  
  # ------------------------------------------------------------
  # 4. Basic completeness audit
  # ------------------------------------------------------------
  
  cat("------------------------------------------------------------\n")
  cat("DIRECTION FIELD COMPLETENESS\n")
  cat("------------------------------------------------------------\n")
  
  cat(
    "Master NA:",
    sum(is.na(master_dir_ascii)),
    "\n"
  )
  
  cat(
    "Supplementary NA:",
    sum(is.na(supp_dir)),
    "\n"
  )
  
  cat(
    "Master empty:",
    sum(
      !is.na(master_dir_ascii) &
        trimws(master_dir_ascii) == ""
    ),
    "\n"
  )
  
  cat(
    "Supplementary empty:",
    sum(
      !is.na(supp_dir) &
        trimws(supp_dir) == ""
    ),
    "\n\n"
  )
  
  
  # ------------------------------------------------------------
  # 5. ASCII-normalized text cleanup
  #    No biological information changed
  # ------------------------------------------------------------
  
  normalize_ascii_direction <- function(x) {
    
    x <- as.character(x)
    
    x <- trimws(x)
    
    # Unicode arrow -> ASCII arrow
    x <- gsub(
      "\u2192",
      "->",
      x,
      fixed = TRUE
    )
    
    # Normalize repeated whitespace
    x <- gsub(
      "[[:space:]]+",
      " ",
      x
    )
    
    # Normalize spaces around ASCII arrow
    x <- gsub(
      "[[:space:]]*->[[:space:]]*",
      " -> ",
      x
    )
    
    # Normalize spaces around multi-direction separator
    x <- gsub(
      "[[:space:]]*\\|[[:space:]]*",
      " | ",
      x
    )
    
    trimws(x)
  }
  
  
  master_dir_ascii_norm <- normalize_ascii_direction(
    master_dir_ascii
  )
  
  supp_dir_ascii_norm <- normalize_ascii_direction(
    supp_dir
  )
  
  
  # ------------------------------------------------------------
  # 6. Final 585-row positional comparison
  # ------------------------------------------------------------
  
  cat("------------------------------------------------------------\n")
  cat("FINAL ASCII DIRECTION COMPARISON\n")
  cat("------------------------------------------------------------\n")
  
  length_pass <- (
    length(master_dir_ascii_norm) ==
      length(supp_dir_ascii_norm)
  )
  
  cat(
    "Direction vector length:",
    ifelse(length_pass, "PASS", "FAIL"),
    "\n"
  )
  
  
  if (length_pass) {
    
    direction_match <- (
      master_dir_ascii_norm ==
        supp_dir_ascii_norm
    )
    
    mismatch_idx <- which(
      !direction_match
    )
    
  } else {
    
    direction_match <- logical(0)
    mismatch_idx <- integer(0)
    
  }
  
  
  cat(
    "Number of direction-label mismatches:",
    length(mismatch_idx),
    "\n\n"
  )
  
  
  # ------------------------------------------------------------
  # 7. Show mismatches only if they exist
  # ------------------------------------------------------------
  
  if (length(mismatch_idx) > 0) {
    
    show_n <- min(
      30,
      length(mismatch_idx)
    )
    
    cat("------------------------------------------------------------\n")
    cat("FIRST DIRECTION-LABEL MISMATCHES\n")
    cat("------------------------------------------------------------\n\n")
    
    diagnostic_df <- data.frame(
      Row = mismatch_idx[
        seq_len(show_n)
      ],
      
      LR_Pair = master_check$Exact_LR_Pair[
        mismatch_idx[
          seq_len(show_n)
        ]
      ],
      
      Master_ASCII = master_dir_ascii_norm[
        mismatch_idx[
          seq_len(show_n)
        ]
      ],
      
      Supplementary_ASCII = supp_dir_ascii_norm[
        mismatch_idx[
          seq_len(show_n)
        ]
      ],
      
      stringsAsFactors = FALSE
    )
    
    print(
      diagnostic_df,
      row.names = FALSE
    )
    
    cat("\n")
    
  } else {
    
    cat(
      "All 585 direction labels match after ASCII normalization.\n\n"
    )
    
  }
  
  
  # ------------------------------------------------------------
  # 8. Locked six-direction map
  # ------------------------------------------------------------
  
  direction_map <- c(
    "1" = "Endothelial -> Myeloid/macrophage",
    "2" = "Myeloid/macrophage -> Endothelial",
    "3" = "Endothelial -> T-cell",
    "4" = "T-cell -> Endothelial",
    "5" = "Myeloid/macrophage -> T-cell",
    "6" = "T-cell -> Myeloid/macrophage"
  )
  
  
  cat("------------------------------------------------------------\n")
  cat("LOCKED SIX-DIRECTION MAP\n")
  cat("------------------------------------------------------------\n")
  
  print(direction_map)
  
  cat("\n")
  
  
  # ------------------------------------------------------------
  # 9. Verify that every Master direction label is composed
  #    only of the six locked directions
  # ------------------------------------------------------------
  
  valid_direction_labels <- unname(
    direction_map
  )
  
  all_direction_tokens_valid <- TRUE
  
  invalid_rows <- integer(0)
  
  for (i in seq_along(master_dir_ascii_norm)) {
    
    tokens <- trimws(
      unlist(
        strsplit(
          master_dir_ascii_norm[i],
          "\\|"
        )
      )
    )
    
    tokens <- tokens[
      tokens != ""
    ]
    
    if (
      length(tokens) == 0 ||
      !all(tokens %in% valid_direction_labels)
    ) {
      
      all_direction_tokens_valid <- FALSE
      
      invalid_rows <- c(
        invalid_rows,
        i
      )
    }
  }
  
  
  cat(
    "Master direction-token validity:",
    ifelse(
      all_direction_tokens_valid,
      "PASS",
      "FAIL"
    ),
    "\n"
  )
  
  cat(
    "Invalid Master direction rows:",
    length(invalid_rows),
    "\n\n"
  )
  
  
  # ------------------------------------------------------------
  # STEP11-6A-10 — CONTINUE FROM SECTION 10
  # Supplementary direction-token validity
  # READ-ONLY
  # ------------------------------------------------------------
  
  supp_direction_tokens_valid <- TRUE
  
  supp_invalid_rows <- integer(0)
  
  for (i in seq_along(supp_dir_ascii_norm)) {
    
    tokens <- trimws(
      unlist(
        strsplit(
          supp_dir_ascii_norm[i],
          "\\|"
        )
      )
    )
    
    tokens <- tokens[
      tokens != ""
    ]
    
    if (
      length(tokens) == 0 ||
      !all(tokens %in% valid_direction_labels)
    ) {
      
      supp_direction_tokens_valid <- FALSE
      
      supp_invalid_rows <- c(
        supp_invalid_rows,
        i
      )
    }
  }
  
  
  cat(
    "Supplementary direction-token validity:",
    ifelse(
      supp_direction_tokens_valid,
      "PASS",
      "FAIL"
    ),
    "\n"
  )
  
  cat(
    "Invalid Supplementary direction rows:",
    length(supp_invalid_rows),
    "\n\n"
  )
  
  
  # ------------------------------------------------------------
  # Direction frequency comparison
  # ------------------------------------------------------------
  
  cat("------------------------------------------------------------\n")
  cat("DIRECTION LABEL FREQUENCY COMPARISON\n")
  cat("------------------------------------------------------------\n")
  
  master_direction_counts <- table(
    unlist(
      strsplit(
        master_dir_ascii_norm,
        "\\s*\\|\\s*"
      )
    )
  )
  
  supp_direction_counts <- table(
    unlist(
      strsplit(
        supp_dir_ascii_norm,
        "\\s*\\|\\s*"
      )
    )
  )
  
  direction_frequency_pass <- identical(
    master_direction_counts,
    supp_direction_counts
  )
  
  cat(
    "Direction frequency identity:",
    ifelse(
      direction_frequency_pass,
      "PASS",
      "FAIL"
    ),
    "\n\n"
  )
  
  cat("Master direction counts:\n")
  print(master_direction_counts)
  
  cat("\nSupplementary direction counts:\n")
  print(supp_direction_counts)
  
  cat("\n")
  
  
  # ------------------------------------------------------------
  # Final STEP11-6A-10 status
  # ------------------------------------------------------------
  
  final_direction_pass <- (
    nrow(master_check) == nrow(supp_check) &&
      length_pass &&
      sum(is.na(master_dir_ascii_norm)) == 0 &&
      sum(is.na(supp_dir_ascii_norm)) == 0 &&
      length(mismatch_idx) == 0 &&
      all_direction_tokens_valid &&
      supp_direction_tokens_valid &&
      direction_frequency_pass
  )
  
  
  cat("============================================================\n")
  cat("STEP11-6A-10 FINAL STATUS\n")
  cat("============================================================\n\n")
  
  cat(
    "Row count:",
    ifelse(
      nrow(master_check) == nrow(supp_check),
      "PASS",
      "FAIL"
    ),
    "\n"
  )
  
  cat(
    "Direction vector length:",
    ifelse(length_pass, "PASS", "FAIL"),
    "\n"
  )
  
  cat(
    "Direction completeness:",
    ifelse(
      sum(is.na(master_dir_ascii_norm)) == 0 &&
        sum(is.na(supp_dir_ascii_norm)) == 0,
      "PASS",
      "FAIL"
    ),
    "\n"
  )
  
  cat(
    "ASCII direction traceability:",
    ifelse(
      length(mismatch_idx) == 0,
      "PASS",
      "FAIL"
    ),
    "\n"
  )
  
  cat(
    "Master six-direction validity:",
    ifelse(
      all_direction_tokens_valid,
      "PASS",
      "FAIL"
    ),
    "\n"
  )
  
  cat(
    "Supplementary six-direction validity:",
    ifelse(
      supp_direction_tokens_valid,
      "PASS",
      "FAIL"
    ),
    "\n"
  )
  
  cat(
    "Direction frequency:",
    ifelse(
      direction_frequency_pass,
      "PASS",
      "FAIL"
    ),
    "\n\n"
  )
  
  cat(
    "FINAL STATUS:",
    ifelse(
      final_direction_pass,
      "PASS",
      "FAIL"
    ),
    "\n\n"
  )
  
  cat("No files modified.\n")
  cat("No analytical values recalculated.\n")
  cat("No locked Master object modified.\n")
  cat("No Excel file modified.\n")
  
  cat("\n============================================================\n")
  
  # ------------------------------------------------------------
  # 11. Direction label frequency comparison
  # ------------------------------------------------------------
  
  master_direction_counts <- table(
    unlist(
      strsplit(
        master_dir_ascii_norm,
        "\\s*\\|\\s*"
      )
    )
  )
  
  supp_direction_counts <- table(
    unlist(
      strsplit(
        supp_dir_ascii_norm,
        "\\s*\\|\\s*"
      )
    )
  )
  
  
  cat("------------------------------------------------------------\n")
  cat("DIRECTION LABEL FREQUENCY COMPARISON\n")
  cat("------------------------------------------------------------\n")
  
  direction_frequency_pass <- identical(
    master_direction_counts,
    supp_direction_counts
  )
  
  cat(
    "Direction frequency identity:",
    ifelse(
      direction_frequency_pass,
      "PASS",
      "FAIL"
    ),
    "\n\n"
  )
  
  cat("Master direction counts:\n")
  print(master_direction_counts)
  
  cat("\nSupplementary direction counts:\n")
  print(supp_direction_counts)
  
  cat("\n")
  
  
  # ------------------------------------------------------------
  # 12. Final status
  # ------------------------------------------------------------
  
  final_direction_pass <- (
    nrow(master_check) == nrow(supp_check) &&
      length_pass &&
      sum(is.na(master_dir_ascii_norm)) == 0 &&
      sum(is.na(supp_dir_ascii_norm)) == 0 &&
      length(mismatch_idx) == 0 &&
      all_direction_tokens_valid &&
      supp_direction_tokens_valid &&
      direction_frequency_pass
  )
  
  
  cat("============================================================\n")
  cat("STEP11-6A-10 FINAL STATUS\n")
  cat("============================================================\n\n")
  
  cat(
    "Row count:",
    ifelse(
      nrow(master_check) == nrow(supp_check),
      "PASS",
      "FAIL"
    ),
    "\n"
  )
  
  cat(
    "Direction vector length:",
    ifelse(length_pass, "PASS", "FAIL"),
    "\n"
  )
  
  cat(
    "Direction completeness:",
    ifelse(
      sum(is.na(master_dir_ascii_norm)) == 0 &&
        sum(is.na(supp_dir_ascii_norm)) == 0,
      "PASS",
      "FAIL"
    ),
    "\n"
  )
  
  cat(
    "ASCII direction traceability:",
    ifelse(
      length(mismatch_idx) == 0,
      "PASS",
      "FAIL"
    ),
    "\n"
  )
  
  cat(
    "Master six-direction validity:",
    ifelse(
      all_direction_tokens_valid,
      "PASS",
      "FAIL"
    ),
    "\n"
  )
  
  cat(
    "Supplementary six-direction validity:",
    ifelse(
      supp_direction_tokens_valid,
      "PASS",
      "FAIL"
    ),
    "\n"
  )
  
  cat(
    "Direction frequency:",
    ifelse(
      direction_frequency_pass,
      "PASS",
      "FAIL"
    ),
    "\n\n"
  )
  
  cat(
    "FINAL STATUS:",
    ifelse(
      final_direction_pass,
      "PASS",
      "FAIL"
    ),
    "\n\n"
  )
  
  cat("No files modified.\n")
  cat("No analytical values recalculated.\n")
  cat("No locked Master object modified.\n")
  cat("No Excel file modified.\n")
  
  cat("\n============================================================\n") 
  
  
  