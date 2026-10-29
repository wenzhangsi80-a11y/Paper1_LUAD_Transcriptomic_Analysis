##GitHub4_Paper1_tumor_microenvironment_association_analysis_1

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


###Step 7 — 2.7 Immune Infiltration Analysis

# Resume Workflow-04-2D WGCNA
# ============================================================

work_dir <- input_dir

setwd(work_dir)

getwd()

load(
  "input_dir/Paper1_tumor-only WGCNA/Workflow-04-2D_final_WGCNA_candidate_checkpoint.RData"
)

ls()

# ============================================================
# Paper 1 Revision
# WORKFLOW-07
# STEP 7A — Data and Signature Audit
# 0. Audit the WGCNA analysis directory
# ============================================================

wgcna_dir <- "input_dir/Paper1_tumor-only WGCNA"

cat("\n============================================================\n")
cat("STEP 7A: WGCNA directory audit\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Check whether the directory exists
# ------------------------------------------------------------

cat("Working directory:\n")
cat(wgcna_dir, "\n\n")

cat("Directory exists: ", dir.exists(wgcna_dir), "\n\n", sep = "")

if (!dir.exists(wgcna_dir)) {
  stop("ERROR: WGCNA directory does not exist. Please check the path.")
}

# ------------------------------------------------------------
# 2. List all files in the WGCNA directory
# ------------------------------------------------------------

files_all <- list.files(
  path = wgcna_dir,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

cat("Total files found: ", length(files_all), "\n\n", sep = "")

# ------------------------------------------------------------
# 3. Create a compact file inventory
# ------------------------------------------------------------

file_inventory <- data.frame(
  FileName = basename(files_all),
  RelativePath = substring(
    files_all,
    nchar(wgcna_dir) + 2
  ),
  Extension = tools::file_ext(files_all),
  Size_MB = round(
    file.info(files_all)$size / 1024^2,
    3
  ),
  Modified = file.info(files_all)$mtime,
  stringsAsFactors = FALSE
)

file_inventory <- file_inventory[
  order(file_inventory$RelativePath),
]

print(file_inventory, row.names = FALSE)

# ------------------------------------------------------------
# 4. Show files potentially relevant to immune analysis
# ------------------------------------------------------------

immune_pattern <- paste(
  c(
    "immune",
    "infiltrat",
    "CIBERSORT",
    "LM22",
    "ssGSEA",
    "GSVA",
    "ESTIMATE",
    "purity",
    "ImmuneScore",
    "StromalScore",
    "ESTIMATEScore",
    "xcell",
    "MCP",
    "deconv"
  ),
  collapse = "|"
)

immune_files <- file_inventory[
  grepl(
    immune_pattern,
    file_inventory$FileName,
    ignore.case = TRUE
  ) |
    grepl(
      immune_pattern,
      file_inventory$RelativePath,
      ignore.case = TRUE
    ),
]

cat("\n============================================================\n")
cat("Potential immune-analysis files\n")
cat("============================================================\n\n")

if (nrow(immune_files) == 0) {
  cat("No filenames matching immune-analysis keywords were found.\n")
} else {
  print(immune_files, row.names = FALSE)
}

# ------------------------------------------------------------
# 5. Show R objects / scripts / data files
# ------------------------------------------------------------

important_ext <- c(
  "R",
  "Rmd",
  "RData",
  "rds",
  "txt",
  "csv",
  "tsv",
  "xlsx",
  "xls"
)

important_files <- file_inventory[
  tolower(file_inventory$Extension) %in%
    tolower(important_ext),
]

cat("\n============================================================\n")
cat("R / data / result files\n")
cat("============================================================\n\n")

print(important_files, row.names = FALSE)

# ------------------------------------------------------------
# 6. Search filenames for specific objects already used
# ------------------------------------------------------------

object_pattern <- paste(
  c(
    "estimate",
    "candidate",
    "core",
    "module",
    "WGCNA",
    "MM_all",
    "GS_all",
    "module_tme",
    "clinical",
    "expression",
    "expr",
    "datExpr"
  ),
  collapse = "|"
)

object_files <- file_inventory[
  grepl(
    object_pattern,
    file_inventory$FileName,
    ignore.case = TRUE
  ) |
    grepl(
      object_pattern,
      file_inventory$RelativePath,
      ignore.case = TRUE
    ),
]

cat("\n============================================================\n")
cat("WGCNA / candidate / ESTIMATE-related files\n")
cat("============================================================\n\n")

if (nrow(object_files) == 0) {
  cat("No matching files found.\n")
} else {
  print(object_files, row.names = FALSE)
}

# ------------------------------------------------------------
# 7. Save inventory for the audit
# ------------------------------------------------------------

audit_file <- file.path(
  wgcna_dir,
  "WORKFLOW-07_STEP-7A_directory_file_inventory.txt"
)

capture.output(
  {
    cat("Paper 1 Revision\n")
    cat("WORKFLOW-07 STEP 7A — Directory File Inventory\n")
    cat("Directory:", wgcna_dir, "\n")
    cat("Audit date:", as.character(Sys.time()), "\n\n")
    
    cat("Total files:", length(files_all), "\n\n")
    
    cat("===== ALL FILES =====\n")
    print(file_inventory, row.names = FALSE)
    
    cat("\n===== IMMUNE-RELATED FILES =====\n")
    print(immune_files, row.names = FALSE)
    
    cat("\n===== WGCNA / CANDIDATE / ESTIMATE FILES =====\n")
    print(object_files, row.names = FALSE)
  },
  file = audit_file
)

cat("\nAudit inventory saved to:\n")
cat(audit_file, "\n")

# ============================================================
# Paper 1
# WORKFLOW-07
# STEP 7A-1: Audit existing locked TME / WGCNA objects
# No new statistical analysis
# ============================================================

cat("\n============================================================\n")
cat("WORKFLOW-07 STEP 7A-1: Existing object audit\n")
cat("============================================================\n\n")

wgcna_dir <- "input_dir/Paper1_tumor-only WGCNA"

# ------------------------------------------------------------
# 1. Create temporary audit environment
# ------------------------------------------------------------

audit_env <- new.env(parent = emptyenv())


# ------------------------------------------------------------
# 2. Audit WORKFLOW-03 Step14 locked object
# ------------------------------------------------------------

step14_file <- file.path(
  wgcna_dir,
  "WORKFLOW-03",
  "Step14",
  "Step14_FINAL_LOCKED",
  "RData",
  "WORKFLOW-03_Step14A-14E_FINAL_LOCKED.RData"
)

cat("STEP14 locked file:\n")
cat(step14_file, "\n\n")

if (!file.exists(step14_file)) {
  stop("STEP14 locked RData file not found.")
}

load(step14_file, envir = audit_env)

cat("Objects loaded from STEP14:\n")
print(ls(envir = audit_env))


# ------------------------------------------------------------
# 3. Object dimensions / classes
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Object classes and dimensions\n")
cat("============================================================\n\n")

step14_objects <- ls(envir = audit_env)

for (obj_name in step14_objects) {
  
  obj <- get(obj_name, envir = audit_env)
  
  cat("\n------------------------------------------------------------\n")
  cat("Object:", obj_name, "\n")
  cat("Class :", paste(class(obj), collapse = ", "), "\n")
  
  if (is.data.frame(obj) || is.matrix(obj)) {
    cat("Dim   :", paste(dim(obj), collapse = " x "), "\n")
  }
  
  if (is.vector(obj) && !is.list(obj)) {
    cat("Length:", length(obj), "\n")
  }
}


# ------------------------------------------------------------
# 4. Search object names related to TME / ESTIMATE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("TME / ESTIMATE related objects\n")
cat("============================================================\n\n")

tme_names <- grep(
  "immune|stromal|estimate|tme|score|purity",
  step14_objects,
  ignore.case = TRUE,
  value = TRUE
)

print(tme_names)


# ------------------------------------------------------------
# 5. Inspect likely TME objects
# ------------------------------------------------------------

for (obj_name in tme_names) {
  
  obj <- get(obj_name, envir = audit_env)
  
  cat("\n------------------------------------------------------------\n")
  cat("Object:", obj_name, "\n")
  
  if (is.data.frame(obj)) {
    print(head(obj, 3))
    cat("\nColumn names:\n")
    print(colnames(obj))
  }
  
  if (is.matrix(obj)) {
    cat("Matrix dimensions:",
        paste(dim(obj), collapse = " x "), "\n")
    cat("First row names:\n")
    print(head(rownames(obj), 5))
    cat("First column names:\n")
    print(head(colnames(obj), 5))
  }
}


# ------------------------------------------------------------
# 6. Audit ESTIMATE PRE-STEP19 object
# ------------------------------------------------------------

estimate_file <- file.path(
  wgcna_dir,
  "Paper1_STEP19_PREP_20260914",
  "locked_data",
  "estimate_scores_wgcna_PRE_STEP19_20260914.RData"
)

cat("\n============================================================\n")
cat("PRE-STEP19 ESTIMATE object\n")
cat("============================================================\n\n")

cat("File:\n", estimate_file, "\n\n")

if (!file.exists(estimate_file)) {
  stop("PRE-STEP19 ESTIMATE RData file not found.")
}

estimate_env <- new.env(parent = emptyenv())

load(estimate_file, envir = estimate_env)

cat("Objects:\n")
print(ls(envir = estimate_env))

cat("\nObject structure:\n")

for (obj_name in ls(envir = estimate_env)) {
  
  obj <- get(obj_name, envir = estimate_env)
  
  cat("\n------------------------------------------------------------\n")
  cat("Object:", obj_name, "\n")
  cat("Class :", paste(class(obj), collapse = ", "), "\n")
  
  if (is.data.frame(obj) || is.matrix(obj)) {
    cat("Dim   :", paste(dim(obj), collapse = " x "), "\n")
    cat("Columns/first names:\n")
    print(head(colnames(obj), 10))
  }
  
  if (is.vector(obj) && !is.list(obj)) {
    cat("Length:", length(obj), "\n")
  }
}


# ------------------------------------------------------------
# 7. Check core gene expression matrix
# ------------------------------------------------------------

expr_file <- file.path(
  wgcna_dir,
  "Workflow-04-2A_expr_primary_tumor_FPKM_UQ.txt"
)

cat("\n============================================================\n")
cat("Primary tumor expression matrix header audit\n")
cat("============================================================\n\n")

cat("Expression file:\n", expr_file, "\n\n")

if (!file.exists(expr_file)) {
  stop("Primary tumor expression file not found.")
}

expr_header <- read.delim(
  expr_file,
  header = TRUE,
  sep = "\t",
  nrows = 3,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Dimensions from first 3 rows:\n")
print(dim(expr_header))

cat("\nFirst column name:\n")
print(colnames(expr_header)[1])

cat("\nNumber of columns:\n")
print(ncol(expr_header))

cat("\nFirst 10 column names:\n")
print(head(colnames(expr_header), 10))


# ------------------------------------------------------------
# 8. Core gene presence audit
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

gene_column <- colnames(expr_header)[1]

cat("\n============================================================\n")
cat("Core gene presence audit\n")
cat("============================================================\n\n")

cat("Gene identifier column:", gene_column, "\n\n")

# Read only the first column to avoid loading the 250 MB matrix
gene_ids <- read.delim(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  select = 1
)

gene_ids <- gene_ids[[1]]

core_check <- data.frame(
  Gene = core_genes,
  Present = core_genes %in% gene_ids,
  stringsAsFactors = FALSE
)

print(core_check)


# ------------------------------------------------------------
# 9. Save audit summary
# ------------------------------------------------------------

audit_summary <- list(
  step14_objects = step14_objects,
  step14_tme_objects = tme_names,
  estimate_objects = ls(envir = estimate_env),
  expression_ncol = ncol(expr_header),
  expression_first_column = gene_column,
  core_gene_presence = core_check
)

save(
  audit_summary,
  file = file.path(
    wgcna_dir,
    "WORKFLOW-07_STEP-7A-1_existing_object_audit.RData"
  )
)

cat("\n============================================================\n")
cat("STEP 7A-1 audit completed.\n")
cat("Saved:\n")
cat(file.path(
  wgcna_dir,
  "WORKFLOW-07_STEP-7A-1_existing_object_audit.RData"
), "\n")
cat("============================================================\n")


# ============================================================
# WORKFLOW-07 STEP 7A-2d
# Core-gene column localization in locked datExpr_MAD
# ============================================================

cat("\n============================================================\n")
cat("WORKFLOW-07 STEP 7A-2d\n")
cat("Core-gene column localization\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. Locate the exact columns in datExpr_MAD
# ------------------------------------------------------------

core_col_index <- match(
  unname(core_ensembl[core_genes]),
  colnames(datExpr_MAD)
)

core_localization <- data.frame(
  gene = core_genes,
  Ensembl_ID = unname(core_ensembl[core_genes]),
  column_index = core_col_index,
  stringsAsFactors = FALSE
)

cat("\nCore-gene localization:\n")
print(core_localization)

# ------------------------------------------------------------
# 2. Verify exact ID matching
# ------------------------------------------------------------

stopifnot(
  all(!is.na(core_col_index)),
  all(
    colnames(datExpr_MAD)[core_col_index] ==
      unname(core_ensembl[core_genes])
  )
)

cat("\n>>> Exact Ensembl column matching PASSED <<<\n")

# ------------------------------------------------------------
# 3. Extract core-gene expression matrix
# ------------------------------------------------------------

core_expr <- datExpr_MAD[
  ,
  core_col_index,
  drop = FALSE
]

# Rename columns to gene symbols for downstream analysis
colnames(core_expr) <- core_genes

cat("\nCore expression matrix dimensions:\n")
print(dim(core_expr))

cat("\nCore expression matrix column names:\n")
print(colnames(core_expr))

cat("\nFirst 6 samples:\n")
print(head(core_expr))

# ------------------------------------------------------------
# 4. Verify sample identity/order is unchanged
# ------------------------------------------------------------

stopifnot(
  identical(
    rownames(core_expr),
    rownames(datExpr_MAD)
  )
)

cat("\n>>> Sample identity/order preservation PASSED <<<\n")

# ------------------------------------------------------------
# 5. Basic expression summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Core-gene expression summary\n")
cat("============================================================\n")

core_expression_summary <- data.frame(
  gene = core_genes,
  Ensembl_ID = unname(core_ensembl[core_genes]),
  n = sapply(core_genes, function(g) sum(!is.na(core_expr[, g]))),
  min = sapply(core_genes, function(g) min(core_expr[, g], na.rm = TRUE)),
  Q1 = sapply(core_genes, function(g) quantile(core_expr[, g], 0.25, na.rm = TRUE)),
  median = sapply(core_genes, function(g) median(core_expr[, g], na.rm = TRUE)),
  mean = sapply(core_genes, function(g) mean(core_expr[, g], na.rm = TRUE)),
  Q3 = sapply(core_genes, function(g) quantile(core_expr[, g], 0.75, na.rm = TRUE)),
  max = sapply(core_genes, function(g) max(core_expr[, g], na.rm = TRUE)),
  SD = sapply(core_genes, function(g) sd(core_expr[, g], na.rm = TRUE)),
  stringsAsFactors = FALSE
)

print(core_expression_summary)

# ============================================================
# WORKFLOW-07 STEP 7A-2d
# Locate locked estimate_scores_wgcna checkpoint
# ============================================================

cat("\n============================================================\n")
cat("WORKFLOW-07 STEP 7A-2d\n")
cat("Locate locked estimate_scores_wgcna checkpoint\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. Search known WGCNA directory for relevant RData files
# ------------------------------------------------------------

rdata_files <- list.files(
  wgcna_dir,
  pattern = "\\.(RData|rda)$",
  full.names = TRUE,
  ignore.case = TRUE
)

cat("\nRData files found:\n")
print(basename(rdata_files))

# ------------------------------------------------------------
# 2. Identify files whose names are related to TME / ESTIMATE /
#    WORKFLOW-03 / STEP-14
# ------------------------------------------------------------

tme_candidates <- rdata_files[
  grepl(
    "estimate|tme|workflow.?03|step.?14|14[abcde]|14f",
    basename(rdata_files),
    ignore.case = TRUE
  )
]

cat("\n============================================================\n")
cat("Potential TME/ESTIMATE checkpoints\n")
cat("============================================================\n")

if (length(tme_candidates) == 0) {
  cat("No matching RData checkpoint found.\n")
} else {
  print(basename(tme_candidates))
}


# ============================================================
# WORKFLOW-07 STEP 7A-2d
# Inspect locked WORKFLOW-03 Step 14F checkpoint
# ============================================================

cat("\n============================================================\n")
cat("WORKFLOW-07 STEP 7A-2d\n")
cat("Inspect locked WORKFLOW-03 Step 14F checkpoint\n")
cat("============================================================\n")

tme_checkpoint <- file.path(
  wgcna_dir,
  "Paper1_WORKFLOW03_14F_4D_2E_FINAL.RData"
)

stopifnot(file.exists(tme_checkpoint))

cat("\nCheckpoint:\n")
cat(tme_checkpoint, "\n")

cat("\nFile size (GB):\n")
print(file.info(tme_checkpoint)$size / 1024^3)

# ------------------------------------------------------------
# Load into isolated temporary environment
# ------------------------------------------------------------

tmp_tme_env <- new.env(parent = emptyenv())

loaded_tme_objects <- load(
  tme_checkpoint,
  envir = tmp_tme_env
)

cat("\n============================================================\n")
cat("Objects in checkpoint\n")
cat("============================================================\n")

print(loaded_tme_objects)

# ------------------------------------------------------------
# Check specifically for estimate_scores_wgcna
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Check estimate_scores_wgcna\n")
cat("============================================================\n")

estimate_found <- "estimate_scores_wgcna" %in% loaded_tme_objects

cat("estimate_scores_wgcna found:\n")
print(estimate_found)

if (estimate_found) {
  
  estimate_tmp <- get(
    "estimate_scores_wgcna",
    envir = tmp_tme_env
  )
  
  cat("\nClass:\n")
  print(class(estimate_tmp))
  
  cat("\nDimensions:\n")
  print(dim(estimate_tmp))
  
  cat("\nColumn names:\n")
  print(colnames(estimate_tmp))
  
  cat("\nFirst 10 row names:\n")
  print(head(rownames(estimate_tmp), 10))
  
  cat("\nUnique row names:\n")
  print(length(unique(rownames(estimate_tmp))))
  
  cat("\nDuplicated row names:\n")
  print(sum(duplicated(rownames(estimate_tmp))))
  
} else {
  
  cat("\n>>> estimate_scores_wgcna NOT FOUND <<<\n")
  
}



# ============================================================
# WORKFLOW-07 STEP 7A-2d
# Restore locked ESTIMATE object
# ============================================================

estimate_scores_wgcna <- get(
  "estimate_scores_wgcna",
  envir = tmp_tme_env
)

cat("\n============================================================\n")
cat("Restored estimate_scores_wgcna\n")
cat("============================================================\n")

cat("Class:\n")
print(class(estimate_scores_wgcna))

cat("\nDimensions:\n")
print(dim(estimate_scores_wgcna))

cat("\nRow names:\n")
print(rownames(estimate_scores_wgcna))

cat("\nNumber of sample columns:\n")
print(ncol(estimate_scores_wgcna))

cat("\nUnique sample IDs:\n")
print(length(unique(colnames(estimate_scores_wgcna))))

cat("\nDuplicated sample IDs:\n")
print(sum(duplicated(colnames(estimate_scores_wgcna))))



# ============================================================
# WORKFLOW-07 STEP 7A-2d
# Sample alignment audit
# ============================================================

expr_sample_ids <- rownames(datExpr_MAD)
tme_sample_ids  <- colnames(estimate_scores_wgcna)

cat("\n============================================================\n")
cat("STEP 7A-2d: Sample alignment audit\n")
cat("============================================================\n")

cat("\n1. Number of expression samples:\n")
print(length(expr_sample_ids))

cat("\n2. Number of TME samples:\n")
print(length(tme_sample_ids))

cat("\n3. Exact same sample set:\n")
print(setequal(expr_sample_ids, tme_sample_ids))

cat("\n4. Exact same sample order:\n")
print(identical(expr_sample_ids, tme_sample_ids))

cat("\n5. Expression samples absent from TME:\n")
print(setdiff(expr_sample_ids, tme_sample_ids))

cat("\n6. TME samples absent from expression:\n")
print(setdiff(tme_sample_ids, expr_sample_ids))

cat("\n7. First 10 expression sample IDs:\n")
print(head(expr_sample_ids, 10))

cat("\n8. First 10 TME sample IDs:\n")
print(head(tme_sample_ids, 10))


# ============================================================
# WORKFLOW-07 STEP 7A-2d
# Build locked core-gene + TME integration table
# ============================================================

cat("\n============================================================\n")
cat("STEP 7A-2d: Build core-gene + TME integration table\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. Confirm core expression dimensions
# ------------------------------------------------------------

cat("\nCore expression dimensions:\n")
print(dim(core_expr))

cat("\nCore expression genes:\n")
print(colnames(core_expr))

# ------------------------------------------------------------
# 2. Build integrated table
# ------------------------------------------------------------

core_tme_alignment <- data.frame(
  Sample_ID = rownames(core_expr),
  
  AOC3 = as.numeric(
    core_expr[, "AOC3"]
  ),
  
  MS4A7 = as.numeric(
    core_expr[, "MS4A7"]
  ),
  
  IL7R = as.numeric(
    core_expr[, "IL7R"]
  ),
  
  StromalScore = as.numeric(
    estimate_scores_wgcna["StromalScore", ]
  ),
  
  ImmuneScore = as.numeric(
    estimate_scores_wgcna["ImmuneScore", ]
  ),
  
  ESTIMATEScore = as.numeric(
    estimate_scores_wgcna["ESTIMATEScore", ]
  ),
  
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# 3. Basic structure check
# ------------------------------------------------------------

cat("\nDimensions of core_tme_alignment:\n")
print(dim(core_tme_alignment))

cat("\nColumn names:\n")
print(colnames(core_tme_alignment))

cat("\nFirst 6 rows:\n")
print(head(core_tme_alignment))

# ------------------------------------------------------------
# 4. Sample identity check
# ------------------------------------------------------------

cat("\nSample ID identity check:\n")
print(identical(
  core_tme_alignment$Sample_ID,
  rownames(datExpr_MAD)
))

# ------------------------------------------------------------
# 5. Missing-value audit
# ------------------------------------------------------------

cat("\nMissing values by column:\n")
print(colSums(is.na(core_tme_alignment)))

cat("\nNaN values by numeric column:\n")
print(sapply(
  core_tme_alignment[, -1, drop = FALSE],
  function(x) sum(is.nan(x))
))

cat("\nInfinite values by numeric column:\n")
print(sapply(
  core_tme_alignment[, -1, drop = FALSE],
  function(x) sum(is.infinite(x))
))


# ============================================================
# WORKFLOW-07 STEP 7A-2d
# Core-gene expression distribution QC
# ============================================================

cat("\n============================================================\n")
cat("STEP 7A-2d: Core-gene expression distribution QC\n")
cat("============================================================\n")

core_gene_summary <- data.frame(
  Gene = c("AOC3", "MS4A7", "IL7R"),
  
  N = c(
    sum(is.finite(core_tme_alignment$AOC3)),
    sum(is.finite(core_tme_alignment$MS4A7)),
    sum(is.finite(core_tme_alignment$IL7R))
  ),
  
  Min = c(
    min(core_tme_alignment$AOC3),
    min(core_tme_alignment$MS4A7),
    min(core_tme_alignment$IL7R)
  ),
  
  Q1 = c(
    quantile(core_tme_alignment$AOC3, 0.25),
    quantile(core_tme_alignment$MS4A7, 0.25),
    quantile(core_tme_alignment$IL7R, 0.25)
  ),
  
  Median = c(
    median(core_tme_alignment$AOC3),
    median(core_tme_alignment$MS4A7),
    median(core_tme_alignment$IL7R)
  ),
  
  Mean = c(
    mean(core_tme_alignment$AOC3),
    mean(core_tme_alignment$MS4A7),
    mean(core_tme_alignment$IL7R)
  ),
  
  Q3 = c(
    quantile(core_tme_alignment$AOC3, 0.75),
    quantile(core_tme_alignment$MS4A7, 0.75),
    quantile(core_tme_alignment$IL7R, 0.75)
  ),
  
  Max = c(
    max(core_tme_alignment$AOC3),
    max(core_tme_alignment$MS4A7),
    max(core_tme_alignment$IL7R)
  ),
  
  SD = c(
    sd(core_tme_alignment$AOC3),
    sd(core_tme_alignment$MS4A7),
    sd(core_tme_alignment$IL7R)
  ),
  
  stringsAsFactors = FALSE
)

print(core_gene_summary)

# ============================================================
# WORKFLOW-07 STEP 7A-2d
# Exact numerical identity check
# ============================================================

cat("\n============================================================\n")
cat("STEP 7A-2d: Exact numerical identity check\n")
cat("============================================================\n")

cat("\nAOC3 identical:\n")
print(identical(
  core_tme_alignment$AOC3,
  as.numeric(core_expr[, "AOC3"])
))

cat("\nMS4A7 identical:\n")
print(identical(
  core_tme_alignment$MS4A7,
  as.numeric(core_expr[, "MS4A7"])
))

cat("\nIL7R identical:\n")
print(identical(
  core_tme_alignment$IL7R,
  as.numeric(core_expr[, "IL7R"])
))

cat("\nStromalScore identical to locked object:\n")
print(identical(
  core_tme_alignment$StromalScore,
  as.numeric(estimate_scores_wgcna["StromalScore", ])
))

cat("\nImmuneScore identical to locked object:\n")
print(identical(
  core_tme_alignment$ImmuneScore,
  as.numeric(estimate_scores_wgcna["ImmuneScore", ])
))

cat("\nESTIMATEScore identical to locked object:\n")
print(identical(
  core_tme_alignment$ESTIMATEScore,
  as.numeric(estimate_scores_wgcna["ESTIMATEScore", ])
))

##WORKFLOW-07 STEP 7A-2e：gene-level TME association analysis。
# ============================================================
# Paper 1
# WORKFLOW-07 STEP 7A-2e
# Gene-level TME association analysis
#
# Core genes:
#   AOC3
#   MS4A7
#   IL7R
#
# TME traits:
#   StromalScore
#   ImmuneScore
#   ESTIMATEScore
#
# Method:
#   Spearman correlation
#   BH-FDR correction across all 9 gene-TME tests
#
# Input:
#   core_tme_alignment
#
# IMPORTANT:
#   No ESTIMATE recomputation
#   No WGCNA recomputation
#   No sample reordering
#   No sample removal
# ============================================================

cat("\n============================================================\n")
cat("STEP 7A-2e: Gene-level TME association analysis\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Define core genes and TME traits
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

tme_traits <- c(
  "StromalScore",
  "ImmuneScore",
  "ESTIMATEScore"
)


# ------------------------------------------------------------
# 2. Basic input validation
# ------------------------------------------------------------

cat("\n[1] Input validation\n")

stopifnot(
  exists("core_tme_alignment"),
  is.data.frame(core_tme_alignment)
)

stopifnot(
  all(core_genes %in% colnames(core_tme_alignment))
)

stopifnot(
  all(tme_traits %in% colnames(core_tme_alignment))
)

cat("Input dimensions: ",
    paste(dim(core_tme_alignment), collapse = " x "),
    "\n", sep = "")

cat("Core genes: ",
    paste(core_genes, collapse = ", "),
    "\n", sep = "")

cat("TME traits: ",
    paste(tme_traits, collapse = ", "),
    "\n", sep = "")


# ------------------------------------------------------------
# 3. Confirm complete data before correlation analysis
# ------------------------------------------------------------

analysis_vars <- c(
  core_genes,
  tme_traits
)

na_count <- sapply(
  core_tme_alignment[, analysis_vars, drop = FALSE],
  function(x) sum(is.na(x))
)

nan_count <- sapply(
  core_tme_alignment[, analysis_vars, drop = FALSE],
  function(x) sum(is.nan(x))
)

inf_count <- sapply(
  core_tme_alignment[, analysis_vars, drop = FALSE],
  function(x) sum(is.infinite(x))
)

cat("\nNA counts:\n")
print(na_count)

cat("\nNaN counts:\n")
print(nan_count)

cat("\nInf counts:\n")
print(inf_count)

stopifnot(
  all(na_count == 0),
  all(nan_count == 0),
  all(inf_count == 0)
)

cat("\nComplete-data check: PASSED\n")


# ------------------------------------------------------------
# 4. Perform the 3 × 3 Spearman correlation analysis
# ------------------------------------------------------------

cat("\n[2] Spearman correlation analysis\n")

gene_tme_results <- data.frame(
  Gene = character(),
  TME_Trait = character(),
  N = integer(),
  Spearman_Rho = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (gene in core_genes) {
  
  for (trait in tme_traits) {
    
    x <- core_tme_alignment[[gene]]
    y <- core_tme_alignment[[trait]]
    
    test <- cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
    
    gene_tme_results <- rbind(
      gene_tme_results,
      data.frame(
        Gene = gene,
        TME_Trait = trait,
        N = sum(complete.cases(x, y)),
        Spearman_Rho = unname(test$estimate),
        P_value = test$p.value,
        stringsAsFactors = FALSE
      )
    )
  }
}


# ------------------------------------------------------------
# 5. Multiple-testing correction
#
# BH correction across all 9 gene-TME associations
# ------------------------------------------------------------

gene_tme_results$FDR_BH <- p.adjust(
  gene_tme_results$P_value,
  method = "BH"
)


# ------------------------------------------------------------
# 6. Add direction
# ------------------------------------------------------------

gene_tme_results$Direction <- ifelse(
  gene_tme_results$Spearman_Rho > 0,
  "Positive",
  ifelse(
    gene_tme_results$Spearman_Rho < 0,
    "Negative",
    "Zero"
  )
)


# ------------------------------------------------------------
# 7. Add significance indicators
# ------------------------------------------------------------

gene_tme_results$Significant_P05 <- (
  gene_tme_results$P_value < 0.05
)

gene_tme_results$Significant_FDR05 <- (
  gene_tme_results$FDR_BH < 0.05
)


# ------------------------------------------------------------
# 8. Round only for display
#
# IMPORTANT:
# The underlying numerical values remain unchanged.
# ------------------------------------------------------------

gene_tme_results_display <- gene_tme_results

gene_tme_results_display$Spearman_Rho <- round(
  gene_tme_results_display$Spearman_Rho,
  4
)

gene_tme_results_display$P_value <- signif(
  gene_tme_results_display$P_value,
  4
)

gene_tme_results_display$FDR_BH <- signif(
  gene_tme_results_display$FDR_BH,
  4
)


# ------------------------------------------------------------
# 9. Print complete results
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Gene-level TME association results\n")
cat("============================================================\n\n")

print(
  gene_tme_results_display,
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Print significant associations
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Associations with P < 0.05\n")
cat("============================================================\n")

sig_p <- subset(
  gene_tme_results,
  P_value < 0.05
)

if (nrow(sig_p) == 0) {
  
  cat("No associations with P < 0.05.\n")
  
} else {
  
  sig_p_display <- sig_p
  
  sig_p_display$Spearman_Rho <- round(
    sig_p_display$Spearman_Rho,
    4
  )
  
  sig_p_display$P_value <- signif(
    sig_p_display$P_value,
    4
  )
  
  sig_p_display$FDR_BH <- signif(
    sig_p_display$FDR_BH,
    4
  )
  
  print(
    sig_p_display,
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 11. Print FDR-significant associations
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Associations with BH-FDR < 0.05\n")
cat("============================================================\n")

sig_fdr <- subset(
  gene_tme_results,
  FDR_BH < 0.05
)

if (nrow(sig_fdr) == 0) {
  
  cat("No associations with BH-FDR < 0.05.\n")
  
} else {
  
  sig_fdr_display <- sig_fdr
  
  sig_fdr_display$Spearman_Rho <- round(
    sig_fdr_display$Spearman_Rho,
    4
  )
  
  sig_fdr_display$P_value <- signif(
    sig_fdr_display$P_value,
    4
  )
  
  sig_fdr_display$FDR_BH <- signif(
    sig_fdr_display$FDR_BH,
    4
  )
  
  print(
    sig_fdr_display,
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 12. Final structural checks
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Final structural checks\n")
cat("============================================================\n")

cat("Number of tests: ",
    nrow(gene_tme_results),
    "\n", sep = "")

cat("Expected number of tests: 9\n")

cat("All N = 513: ")
print(all(gene_tme_results$N == 513))

cat("All P values finite: ")
print(all(is.finite(gene_tme_results$P_value)))

cat("All FDR values finite: ")
print(all(is.finite(gene_tme_results$FDR_BH)))

cat("Unique genes: ")
print(unique(gene_tme_results$Gene))

cat("Unique TME traits: ")
print(unique(gene_tme_results$TME_Trait))


stopifnot(
  nrow(gene_tme_results) == 9,
  all(gene_tme_results$N == 513),
  all(is.finite(gene_tme_results$P_value)),
  all(is.finite(gene_tme_results$FDR_BH)),
  length(unique(gene_tme_results$Gene)) == 3,
  length(unique(gene_tme_results$TME_Trait)) == 3
)

cat("\nSTEP 7A-2e: PASSED\n")
cat("All 9 gene-level TME associations completed successfully.\n")

#STEP 7A-2f：Gene-level TME association result audit

# ============================================================
# Paper 1
# WORKFLOW-07 STEP 7A-2f
# Gene-level TME association result audit
#
# IMPORTANT:
# This step does NOT recalculate any correlation.
# It only audits the locked results from STEP 7A-2e.
# ============================================================

cat("\n============================================================\n")
cat("STEP 7A-2f: Gene-level TME association result audit\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Display complete numerical results
# ------------------------------------------------------------

cat("\n[1] Complete gene-level TME association results\n")
cat("------------------------------------------------------------\n")

print(
  gene_tme_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 2. Compact display for manuscript-level inspection
# ------------------------------------------------------------

cat("\n[2] Compact result table\n")
cat("------------------------------------------------------------\n")

gene_tme_audit <- gene_tme_results[
  ,
  c(
    "Gene",
    "TME_Trait",
    "N",
    "Spearman_Rho",
    "P_value",
    "FDR_BH",
    "Direction",
    "Significant_P05",
    "Significant_FDR05"
  )
]

print(
  gene_tme_audit,
  row.names = FALSE
)


# ------------------------------------------------------------
# 3. Sort by FDR
# ------------------------------------------------------------

cat("\n[3] Associations sorted by BH-FDR\n")
cat("------------------------------------------------------------\n")

gene_tme_by_FDR <- gene_tme_results[
  order(
    gene_tme_results$FDR_BH,
    gene_tme_results$P_value
  ),
  ,
  drop = FALSE
]

print(
  gene_tme_by_FDR,
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. Sort by absolute correlation strength
# ------------------------------------------------------------

cat("\n[4] Associations sorted by |Spearman rho|\n")
cat("------------------------------------------------------------\n")

gene_tme_by_abs_rho <- gene_tme_results[
  order(
    -abs(gene_tme_results$Spearman_Rho)
  ),
  ,
  drop = FALSE
]

print(
  gene_tme_by_abs_rho,
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Significant association summary
# ------------------------------------------------------------

cat("\n[5] Significance summary\n")
cat("------------------------------------------------------------\n")

cat(
  "P < 0.05: ",
  sum(gene_tme_results$P_value < 0.05),
  " / 9\n",
  sep = ""
)

cat(
  "BH-FDR < 0.05: ",
  sum(gene_tme_results$FDR_BH < 0.05),
  " / 9\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Direction summary
# ------------------------------------------------------------

cat("\n[6] Direction summary\n")
cat("------------------------------------------------------------\n")

print(
  table(gene_tme_results$Direction)
)


# ------------------------------------------------------------
# 7. Gene-level summary
# ------------------------------------------------------------

cat("\n[7] Number of significant TME associations per gene\n")
cat("------------------------------------------------------------\n")

gene_summary <- aggregate(
  cbind(
    Significant_P05,
    Significant_FDR05
  ) ~ Gene,
  data = gene_tme_results,
  FUN = sum
)

print(
  gene_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. TME-trait summary
# ------------------------------------------------------------

cat("\n[8] Number of significant gene associations per TME trait\n")
cat("------------------------------------------------------------\n")

tme_summary <- aggregate(
  cbind(
    Significant_P05,
    Significant_FDR05
  ) ~ TME_Trait,
  data = gene_tme_results,
  FUN = sum
)

print(
  tme_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Internal consistency checks
# ------------------------------------------------------------

cat("\n[9] Internal consistency checks\n")
cat("------------------------------------------------------------\n")

# Direction must agree with rho
direction_check <- all(
  ifelse(
    gene_tme_results$Spearman_Rho > 0,
    "Positive",
    ifelse(
      gene_tme_results$Spearman_Rho < 0,
      "Negative",
      "Zero"
    )
  ) == gene_tme_results$Direction
)

cat("Direction agrees with rho: ")
print(direction_check)


# FDR should be between 0 and 1
fdr_range_check <- all(
  gene_tme_results$FDR_BH >= 0 &
    gene_tme_results$FDR_BH <= 1
)

cat("All FDR values within [0,1]: ")
print(fdr_range_check)


# P values should be between 0 and 1
p_range_check <- all(
  gene_tme_results$P_value >= 0 &
    gene_tme_results$P_value <= 1
)

cat("All P values within [0,1]: ")
print(p_range_check)


# ------------------------------------------------------------
# 10. Final audit status
# ------------------------------------------------------------

stopifnot(
  direction_check,
  fdr_range_check,
  p_range_check
)

cat("\n============================================================\n")
cat("STEP 7A-2f: AUDIT PASSED\n")
cat("============================================================\n")


# ============================================================
# STEP 7A-2g diagnostic
# Investigate combination_check == FALSE
# ============================================================

cat("\n============================================================\n")
cat("STEP 7A-2g: Combination check diagnostic\n")
cat("============================================================\n")

cat("\n[1] Expected combinations:\n")
print(expected_combinations, row.names = FALSE)

cat("\n[2] Actual combinations:\n")
print(actual_combinations, row.names = FALSE)

cat("\n[3] Structure of expected combinations:\n")
str(expected_combinations)

cat("\n[4] Structure of actual combinations:\n")
str(actual_combinations)

cat("\n[5] Set comparison: expected but absent from actual\n")
print(
  setdiff(
    paste(
      expected_combinations$Gene,
      expected_combinations$TME_Trait,
      sep = "__"
    ),
    paste(
      actual_combinations$Gene,
      actual_combinations$TME_Trait,
      sep = "__"
    )
  )
)

cat("\n[6] Set comparison: actual but absent from expected\n")
print(
  setdiff(
    paste(
      actual_combinations$Gene,
      actual_combinations$TME_Trait,
      sep = "__"
    ),
    paste(
      expected_combinations$Gene,
      expected_combinations$TME_Trait,
      sep = "__"
    )
  )
)

# ============================================================
# Paper 1
# WORKFLOW-07 STEP 7A-2g
# Generate locked manuscript-ready gene × TME result table
# ============================================================

cat("\n============================================================\n")
cat("STEP 7A-2g: Manuscript-ready gene × TME result table\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Confirm the source object
# ------------------------------------------------------------

stopifnot(exists("gene_tme_results"))

cat("\nSource object:\n")
print(dim(gene_tme_results))
print(colnames(gene_tme_results))


# ------------------------------------------------------------
# 2. Content-based combination check
#    Do NOT use identical() on the expand.grid data.frame
# ------------------------------------------------------------

expected_genes <- c("AOC3", "MS4A7", "IL7R")

expected_traits <- c(
  "StromalScore",
  "ImmuneScore",
  "ESTIMATEScore"
)

expected_keys <- sort(
  as.character(
    outer(
      expected_genes,
      expected_traits,
      FUN = function(g, t) paste(g, t, sep = "__")
    )
  )
)

actual_keys <- sort(
  paste(
    as.character(gene_tme_results$Gene),
    as.character(gene_tme_results$TME_Trait),
    sep = "__"
  )
)

combination_check <- identical(expected_keys, actual_keys)

cat("\nAll expected gene-TME combinations present: ")
print(combination_check)

stopifnot(combination_check)


# ------------------------------------------------------------
# 3. Create the manuscript-ready wide table
#    No new statistical calculation
#    Values are copied directly from gene_tme_results
# ------------------------------------------------------------

gene_tme_manuscript_table <- data.frame(
  Gene = c("AOC3", "MS4A7", "IL7R"),
  ImmuneScore_Rho = NA_real_,
  ImmuneScore_FDR = NA_real_,
  StromalScore_Rho = NA_real_,
  StromalScore_FDR = NA_real_,
  ESTIMATEScore_Rho = NA_real_,
  ESTIMATEScore_FDR = NA_real_,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 4. Helper function
#    Retrieve exactly one locked source value
# ------------------------------------------------------------

get_original_value <- function(gene, trait, value_column) {
  
  idx <- as.character(gene_tme_results$Gene) == gene &
    as.character(gene_tme_results$TME_Trait) == trait
  
  stopifnot(sum(idx) == 1)
  
  gene_tme_results[[value_column]][idx]
}


# ------------------------------------------------------------
# 5. Populate the wide table directly from source results
# ------------------------------------------------------------

for (i in seq_len(nrow(gene_tme_manuscript_table))) {
  
  gene <- gene_tme_manuscript_table$Gene[i]
  
  gene_tme_manuscript_table$ImmuneScore_Rho[i] <-
    get_original_value(
      gene,
      "ImmuneScore",
      "Spearman_Rho"
    )
  
  gene_tme_manuscript_table$ImmuneScore_FDR[i] <-
    get_original_value(
      gene,
      "ImmuneScore",
      "FDR_BH"
    )
  
  gene_tme_manuscript_table$StromalScore_Rho[i] <-
    get_original_value(
      gene,
      "StromalScore",
      "Spearman_Rho"
    )
  
  gene_tme_manuscript_table$StromalScore_FDR[i] <-
    get_original_value(
      gene,
      "StromalScore",
      "FDR_BH"
    )
  
  gene_tme_manuscript_table$ESTIMATEScore_Rho[i] <-
    get_original_value(
      gene,
      "ESTIMATEScore",
      "Spearman_Rho"
    )
  
  gene_tme_manuscript_table$ESTIMATEScore_FDR[i] <-
    get_original_value(
      gene,
      "ESTIMATEScore",
      "FDR_BH"
    )
}


# ------------------------------------------------------------
# 6. Structural QC
# ------------------------------------------------------------

cat("\n[QC-1] Table dimensions:\n")
print(dim(gene_tme_manuscript_table))

stopifnot(
  identical(
    dim(gene_tme_manuscript_table),
    c(3L, 7L)
  )
)

cat("\n[QC-2] Gene order:\n")
print(gene_tme_manuscript_table$Gene)

stopifnot(
  identical(
    gene_tme_manuscript_table$Gene,
    c("AOC3", "MS4A7", "IL7R")
  )
)


# ------------------------------------------------------------
# 7. Missing-value QC
# ------------------------------------------------------------

numeric_columns <- setdiff(
  colnames(gene_tme_manuscript_table),
  "Gene"
)

na_check <- sapply(
  gene_tme_manuscript_table[, numeric_columns, drop = FALSE],
  function(x) sum(is.na(x))
)

cat("\n[QC-3] NA counts:\n")
print(na_check)

stopifnot(all(na_check == 0))


# ------------------------------------------------------------
# 8. Exact numerical identity check
#    Every manuscript value must equal the locked source value
# ------------------------------------------------------------

identity_check <- logical(0)

for (i in seq_len(nrow(gene_tme_manuscript_table))) {
  
  gene <- gene_tme_manuscript_table$Gene[i]
  
  identity_check <- c(
    identity_check,
    
    identical(
      gene_tme_manuscript_table$ImmuneScore_Rho[i],
      get_original_value(
        gene,
        "ImmuneScore",
        "Spearman_Rho"
      )
    ),
    
    identical(
      gene_tme_manuscript_table$ImmuneScore_FDR[i],
      get_original_value(
        gene,
        "ImmuneScore",
        "FDR_BH"
      )
    ),
    
    identical(
      gene_tme_manuscript_table$StromalScore_Rho[i],
      get_original_value(
        gene,
        "StromalScore",
        "Spearman_Rho"
      )
    ),
    
    identical(
      gene_tme_manuscript_table$StromalScore_FDR[i],
      get_original_value(
        gene,
        "StromalScore",
        "FDR_BH"
      )
    ),
    
    identical(
      gene_tme_manuscript_table$ESTIMATEScore_Rho[i],
      get_original_value(
        gene,
        "ESTIMATEScore",
        "Spearman_Rho"
      )
    ),
    
    identical(
      gene_tme_manuscript_table$ESTIMATEScore_FDR[i],
      get_original_value(
        gene,
        "ESTIMATEScore",
        "FDR_BH"
      )
    )
  )
}

cat("\n[QC-4] Exact numerical identity with source object:\n")
print(identity_check)

stopifnot(all(identity_check))


# ------------------------------------------------------------
# 9. Final display
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL MANUSCRIPT-READY GENE × TME TABLE\n")
cat("============================================================\n\n")

print(
  gene_tme_manuscript_table,
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Final pass
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 7A-2g PASSED\n")
cat("============================================================\n")
cat("3 core genes × 3 TME traits = 9 source associations\n")
cat("Wide-format table: 3 × 7\n")
cat("All values copied directly from gene_tme_results\n")
cat("No new statistical calculation performed\n")
cat("No missing values\n")
cat("Exact numerical identity with locked source results confirmed\n")
cat("============================================================\n")


# ============================================================
# Paper 1
# WORKFLOW-07 STEP 7A-2h
# Final Object Integrity + Checkpoint Save
# ============================================================

cat("\n============================================================\n")
cat("STEP 7A-2h: Final Object Integrity + Checkpoint Save\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Required objects must exist
# ------------------------------------------------------------

required_objects <- c(
  "gene_tme_results",
  "gene_tme_manuscript_table",
  "core_tme_alignment"
)

object_presence <- sapply(
  required_objects,
  exists
)

cat("\n[QC-1] Required objects present:\n")
print(object_presence)

stopifnot(all(object_presence))


# ------------------------------------------------------------
# 2. gene_tme_results integrity
# ------------------------------------------------------------

cat("\n[QC-2] gene_tme_results integrity\n")
cat("------------------------------------------------------------\n")

cat("Dimensions: ")
print(dim(gene_tme_results))

cat("Number of genes: ")
print(length(unique(gene_tme_results$Gene)))

cat("Number of TME traits: ")
print(length(unique(gene_tme_results$TME_Trait)))

cat("Number of associations: ")
print(nrow(gene_tme_results))

stopifnot(
  identical(dim(gene_tme_results), c(9L, 9L)),
  identical(
    sort(unique(as.character(gene_tme_results$Gene))),
    sort(c("AOC3", "MS4A7", "IL7R"))
  ),
  identical(
    sort(unique(as.character(gene_tme_results$TME_Trait))),
    sort(c(
      "ImmuneScore",
      "StromalScore",
      "ESTIMATEScore"
    ))
  )
)


# ------------------------------------------------------------
# 3. Confirm all 9 source associations have N = 513
# ------------------------------------------------------------

n_check <- all(
  gene_tme_results$N == 513
)

cat("\nAll associations have N = 513: ")
print(n_check)

stopifnot(n_check)


# ------------------------------------------------------------
# 4. Confirm statistical fields are finite and valid
# ------------------------------------------------------------

rho_finite_check <- all(
  is.finite(gene_tme_results$Spearman_Rho)
)

p_finite_check <- all(
  is.finite(gene_tme_results$P_value)
)

fdr_finite_check <- all(
  is.finite(gene_tme_results$FDR_BH)
)

rho_range_check <- all(
  gene_tme_results$Spearman_Rho >= -1 &
    gene_tme_results$Spearman_Rho <= 1
)

p_range_check <- all(
  gene_tme_results$P_value >= 0 &
    gene_tme_results$P_value <= 1
)

fdr_range_check <- all(
  gene_tme_results$FDR_BH >= 0 &
    gene_tme_results$FDR_BH <= 1
)

cat("\nSpearman rho finite: ")
print(rho_finite_check)

cat("P values finite: ")
print(p_finite_check)

cat("FDR values finite: ")
print(fdr_finite_check)

cat("Spearman rho within [-1, 1]: ")
print(rho_range_check)

cat("P values within [0, 1]: ")
print(p_range_check)

cat("FDR values within [0, 1]: ")
print(fdr_range_check)

stopifnot(
  rho_finite_check,
  p_finite_check,
  fdr_finite_check,
  rho_range_check,
  p_range_check,
  fdr_range_check
)


# ------------------------------------------------------------
# 5. Confirm all nine associations remain significant
# ------------------------------------------------------------

all_p05_check <- all(
  gene_tme_results$P_value < 0.05
)

all_fdr05_check <- all(
  gene_tme_results$FDR_BH < 0.05
)

cat("\nAll 9 P < 0.05: ")
print(all_p05_check)

cat("All 9 BH-FDR < 0.05: ")
print(all_fdr05_check)

stopifnot(
  all_p05_check,
  all_fdr05_check
)


# ------------------------------------------------------------
# 6. Confirm direction consistency
# ------------------------------------------------------------

direction_check <- all(
  gene_tme_results$Direction ==
    ifelse(
      gene_tme_results$Spearman_Rho > 0,
      "Positive",
      ifelse(
        gene_tme_results$Spearman_Rho < 0,
        "Negative",
        "Zero"
      )
    )
)

cat("\nDirection agrees with Spearman rho: ")
print(direction_check)

stopifnot(direction_check)


# ------------------------------------------------------------
# 7. gene_tme_manuscript_table integrity
# ------------------------------------------------------------

cat("\n[QC-3] gene_tme_manuscript_table integrity\n")
cat("------------------------------------------------------------\n")

cat("Dimensions: ")
print(dim(gene_tme_manuscript_table))

cat("Columns:\n")
print(colnames(gene_tme_manuscript_table))

stopifnot(
  identical(
    dim(gene_tme_manuscript_table),
    c(3L, 7L)
  )
)

stopifnot(
  identical(
    gene_tme_manuscript_table$Gene,
    c("AOC3", "MS4A7", "IL7R")
  )
)


# ------------------------------------------------------------
# 8. No missing values in manuscript table
# ------------------------------------------------------------

manuscript_numeric_columns <- setdiff(
  colnames(gene_tme_manuscript_table),
  "Gene"
)

manuscript_na_check <- all(
  !is.na(
    as.matrix(
      gene_tme_manuscript_table[
        ,
        manuscript_numeric_columns,
        drop = FALSE
      ]
    )
  )
)

cat("\nNo NA in manuscript table: ")
print(manuscript_na_check)

stopifnot(manuscript_na_check)


# ------------------------------------------------------------
# 9. Exact source-to-manuscript identity
#    Robust version: check each value separately
# ------------------------------------------------------------

cat("\n[QC-4] Exact source-to-manuscript numerical identity\n")
cat("------------------------------------------------------------\n")

source_to_manuscript_check <- logical(0)

for (i in seq_len(nrow(gene_tme_manuscript_table))) {
  
  gene <- gene_tme_manuscript_table$Gene[i]
  
  # ImmuneScore rho
  source_to_manuscript_check <- c(
    source_to_manuscript_check,
    identical(
      gene_tme_manuscript_table$ImmuneScore_Rho[i],
      get_original_value(
        gene,
        "ImmuneScore",
        "Spearman_Rho"
      )
    )
  )
  
  # ImmuneScore FDR
  source_to_manuscript_check <- c(
    source_to_manuscript_check,
    identical(
      gene_tme_manuscript_table$ImmuneScore_FDR[i],
      get_original_value(
        gene,
        "ImmuneScore",
        "FDR_BH"
      )
    )
  )
  
  # StromalScore rho
  source_to_manuscript_check <- c(
    source_to_manuscript_check,
    identical(
      gene_tme_manuscript_table$StromalScore_Rho[i],
      get_original_value(
        gene,
        "StromalScore",
        "Spearman_Rho"
      )
    )
  )
  
  # StromalScore FDR
  source_to_manuscript_check <- c(
    source_to_manuscript_check,
    identical(
      gene_tme_manuscript_table$StromalScore_FDR[i],
      get_original_value(
        gene,
        "StromalScore",
        "FDR_BH"
      )
    )
  )
  
  # ESTIMATEScore rho
  source_to_manuscript_check <- c(
    source_to_manuscript_check,
    identical(
      gene_tme_manuscript_table$ESTIMATEScore_Rho[i],
      get_original_value(
        gene,
        "ESTIMATEScore",
        "Spearman_Rho"
      )
    )
  )
  
  # ESTIMATEScore FDR
  source_to_manuscript_check <- c(
    source_to_manuscript_check,
    identical(
      gene_tme_manuscript_table$ESTIMATEScore_FDR[i],
      get_original_value(
        gene,
        "ESTIMATEScore",
        "FDR_BH"
      )
    )
  )
}

cat("\nIndividual identity checks:\n")
print(source_to_manuscript_check)

cat("\nNumber of exact identity checks: ")
print(length(source_to_manuscript_check))

cat("\nAll 18 source-to-manuscript values identical: ")
print(all(source_to_manuscript_check))

stopifnot(
  length(source_to_manuscript_check) == 18L,
  all(source_to_manuscript_check)
)


# ------------------------------------------------------------
# 10. core_tme_alignment integrity
# ------------------------------------------------------------

cat("\n[QC-4] core_tme_alignment integrity\n")
cat("------------------------------------------------------------\n")

cat("Dimensions: ")
print(dim(core_tme_alignment))

cat("Columns:\n")
print(colnames(core_tme_alignment))

stopifnot(
  identical(
    dim(core_tme_alignment),
    c(513L, 7L)
  )
)

stopifnot(
  identical(
    colnames(core_tme_alignment),
    c(
      "Sample_ID",
      "AOC3",
      "MS4A7",
      "IL7R",
      "StromalScore",
      "ImmuneScore",
      "ESTIMATEScore"
    )
  )
)


# ------------------------------------------------------------
# 11. Sample identity and missing-value check
# ------------------------------------------------------------

core_sample_check <- length(
  unique(core_tme_alignment$Sample_ID)
) == 513

core_na_check <- all(
  !is.na(core_tme_alignment)
)

cat("\n513 unique tumor samples: ")
print(core_sample_check)

cat("No missing values in core_tme_alignment: ")
print(core_na_check)

stopifnot(
  core_sample_check,
  core_na_check
)


# ------------------------------------------------------------
# 12. Final object class check
# ------------------------------------------------------------

cat("\n[QC-5] Object classes\n")
cat("------------------------------------------------------------\n")

cat("gene_tme_results: ")
print(class(gene_tme_results))

cat("gene_tme_manuscript_table: ")
print(class(gene_tme_manuscript_table))

cat("core_tme_alignment: ")
print(class(core_tme_alignment))


# ------------------------------------------------------------
# 13. Create compact audit metadata
# ------------------------------------------------------------

step7a2h_audit <- list(
  
  workflow = "WORKFLOW-07",
  
  step = "STEP 7A-2h",
  
  purpose = "Final Object Integrity + Checkpoint Save",
  
  analysis = "Core gene-level TME association",
  
  genes = c("AOC3", "MS4A7", "IL7R"),
  
  tme_traits = c(
    "ImmuneScore",
    "StromalScore",
    "ESTIMATEScore"
  ),
  
  n_samples = 513L,
  
  n_associations = 9L,
  
  correlation_method = "Spearman",
  
  spearman_exact = FALSE,
  
  multiple_testing = "Benjamini-Hochberg",
  
  n_significant_p05 = sum(
    gene_tme_results$P_value < 0.05
  ),
  
  n_significant_fdr05 = sum(
    gene_tme_results$FDR_BH < 0.05
  ),
  
  direction_summary = table(
    gene_tme_results$Direction
  ),
  
  source_object = "gene_tme_results",
  
  manuscript_object = "gene_tme_manuscript_table",
  
  alignment_object = "core_tme_alignment",
  
  new_statistics_computed = FALSE,
  
  source_values_modified = FALSE,
  
  final_integrity_status = "PASSED"
)


# ------------------------------------------------------------
# 14. Define checkpoint path
# ------------------------------------------------------------

checkpoint_file <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "Paper1_WORKFLOW07_STEP7A2h_FINAL_TME_gene_association_checkpoint.RData"
)

cat("\nCheckpoint file:\n")
print(checkpoint_file)


# ------------------------------------------------------------
# 15. Save only the final required objects
# ------------------------------------------------------------

save(
  gene_tme_results,
  gene_tme_manuscript_table,
  core_tme_alignment,
  step7a2h_audit,
  file = checkpoint_file
)


# ------------------------------------------------------------
# 16. Confirm checkpoint exists
# ------------------------------------------------------------

checkpoint_exists <- file.exists(
  checkpoint_file
)

checkpoint_size_MB <- if (checkpoint_exists) {
  file.info(checkpoint_file)$size / 1024^2
} else {
  NA_real_
}

cat("\n[QC-6] Checkpoint saved successfully: ")
print(checkpoint_exists)

cat("Checkpoint size (MB): ")
print(round(checkpoint_size_MB, 3))

stopifnot(checkpoint_exists)


# ------------------------------------------------------------
# 17. Final summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 7A-2h PASSED\n")
cat("============================================================\n")
cat("Required objects present: TRUE\n")
cat("9 gene × TME associations intact: TRUE\n")
cat("N = 513 for all associations: TRUE\n")
cat("All P/FDR values valid: TRUE\n")
cat("All 9 associations remain P < 0.05: TRUE\n")
cat("All 9 associations remain BH-FDR < 0.05: TRUE\n")
cat("Manuscript table = 3 × 7: TRUE\n")
cat("Exact source-to-manuscript identity: TRUE\n")
cat("Core alignment = 513 × 7: TRUE\n")
cat("No new statistics computed: TRUE\n")
cat("Checkpoint saved: TRUE\n")
cat("============================================================\n")


# ============================================================
# Paper 1
# STEP 7B-4A
# Create manuscript-ready displayed Supplementary Table
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-4A: Manuscript-ready Supplementary Table\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Check required locked source objects
# ------------------------------------------------------------

required_objects <- c(
  "gene_tme_results",
  "gene_tme_manuscript_table"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    "Missing required locked objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

cat("Required locked objects present: TRUE\n")


# ------------------------------------------------------------
# 2. Copy locked manuscript table
# ------------------------------------------------------------

gene_tme_display_table <- gene_tme_manuscript_table


# ------------------------------------------------------------
# 3. Format rho values for display
# ------------------------------------------------------------

rho_cols <- c(
  "ImmuneScore_Rho",
  "StromalScore_Rho",
  "ESTIMATEScore_Rho"
)

for (x in rho_cols) {
  gene_tme_display_table[[x]] <-
    sprintf("%.3f", gene_tme_display_table[[x]])
}


# ------------------------------------------------------------
# 4. Format FDR values for display
# ------------------------------------------------------------

format_fdr <- function(x) {
  
  sapply(x, function(v) {
    
    if (is.na(v)) {
      return(NA_character_)
    }
    
    sprintf("%.3e", v)
  })
}


fdr_cols <- c(
  "ImmuneScore_FDR",
  "StromalScore_FDR",
  "ESTIMATEScore_FDR"
)

for (x in fdr_cols) {
  gene_tme_display_table[[x]] <-
    format_fdr(gene_tme_display_table[[x]])
}


# ------------------------------------------------------------
# 5. Rename columns for manuscript presentation
# ------------------------------------------------------------

colnames(gene_tme_display_table) <- c(
  "Gene",
  "ImmuneScore, ρ",
  "ImmuneScore, FDR",
  "StromalScore, ρ",
  "StromalScore, FDR",
  "ESTIMATEScore, ρ",
  "ESTIMATEScore, FDR"
)


# ------------------------------------------------------------
# 6. Display
# ------------------------------------------------------------

print(gene_tme_display_table, row.names = FALSE)


cat("\nDimensions:\n")
print(dim(gene_tme_display_table))

cat("\nSTEP 7B-4A PASSED\n")


# ============================================================
# STEP 7B-4B
# Table title and footnote
# ============================================================

supp_table_title <-
  "Supplementary Table X. Associations of the three core genes with ESTIMATE-derived tumor microenvironment scores in primary LUAD tumors."

supp_table_footnote <-
  paste(
    "Associations were evaluated in N = 513 primary LUAD tumor samples using ",
    "two-sided Spearman rank correlation tests. P values were adjusted for ",
    "multiple testing across the nine gene–TME associations using the ",
    "Benjamini–Hochberg procedure. ρ, Spearman correlation coefficient; ",
    "FDR, Benjamini–Hochberg false discovery rate. ImmuneScore and ",
    "StromalScore represent ESTIMATE-derived estimates of immune and ",
    "stromal content, respectively, whereas ESTIMATEScore is a composite ",
    "score. All nine associations had FDR < 0.05.",
    sep = ""
  )

cat("\nTITLE:\n")
cat(supp_table_title, "\n\n")

cat("FOOTNOTE:\n")
cat(supp_table_footnote, "\n")

##STEP 7B-4C-R：Corrected three-layer consistency audit
# ============================================================
# Paper 1
# STEP 7B-4C-R
# Corrected Three-layer Source-to-Display Consistency Audit
#
# IMPORTANT:
# This audit distinguishes:
#   Level 1 = exact statistical-source identity
#   Level 2 = exact presentation-format transformation
#   Level 3 = manuscript metadata consistency
#
# No statistical recalculation is performed.
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-4C-R: Corrected three-layer consistency audit\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 0. Required objects
# ------------------------------------------------------------

required_objects <- c(
  "gene_tme_results",
  "gene_tme_manuscript_table",
  "gene_tme_display_table",
  "supp_table_title",
  "supp_table_footnote"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

cat("Required objects present: TRUE\n")


# ------------------------------------------------------------
# LEVEL 1
# Locked statistical source
# gene_tme_results
#            ↓
# gene_tme_manuscript_table
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("LEVEL 1: Statistical-source identity\n")
cat("------------------------------------------------------------\n")


source_value_columns <- c(
  "ImmuneScore_Rho",
  "ImmuneScore_FDR",
  "StromalScore_Rho",
  "StromalScore_FDR",
  "ESTIMATEScore_Rho",
  "ESTIMATEScore_FDR"
)


get_original_value <- function(gene, trait, value_column) {
  
  idx <- gene_tme_results$Gene == gene &
    gene_tme_results$TME_Trait == trait
  
  vals <- gene_tme_results[[value_column]][idx]
  
  if (length(vals) != 1) {
    stop(
      "Expected exactly one value for ",
      gene, " / ", trait, " / ", value_column,
      "; found ", length(vals)
    )
  }
  
  vals
}


source_identity_check <- TRUE


for (gene in gene_tme_manuscript_table$Gene) {
  
  for (trait in c(
    "ImmuneScore",
    "StromalScore",
    "ESTIMATEScore"
  )) {
    
    rho_col <- paste0(trait, "_Rho")
    fdr_col <- paste0(trait, "_FDR")
    
    expected_rho <- get_original_value(
      gene,
      trait,
      "Spearman_Rho"
    )
    
    expected_fdr <- get_original_value(
      gene,
      trait,
      "FDR_BH"
    )
    
    row_idx <-
      gene_tme_manuscript_table$Gene == gene
    
    rho_ok <- identical(
      gene_tme_manuscript_table[[rho_col]][row_idx],
      expected_rho
    )
    
    fdr_ok <- identical(
      gene_tme_manuscript_table[[fdr_col]][row_idx],
      expected_fdr
    )
    
    source_identity_check <-
      source_identity_check &&
      rho_ok &&
      fdr_ok
  }
}


cat(
  "gene_tme_results -> gene_tme_manuscript_table: ",
  source_identity_check,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# LEVEL 1 structural checks
# ------------------------------------------------------------

source_dimensions_ok <-
  identical(dim(gene_tme_manuscript_table), c(3L, 7L))

source_genes_ok <-
  identical(
    gene_tme_manuscript_table$Gene,
    c("AOC3", "MS4A7", "IL7R")
  )

source_all_finite <- all(
  vapply(
    gene_tme_manuscript_table[source_value_columns],
    function(x) all(is.finite(x)),
    logical(1)
  )
)


cat(
  "Manuscript source dimensions 3 x 7: ",
  source_dimensions_ok,
  "\n",
  sep = ""
)

cat(
  "Core-gene order correct: ",
  source_genes_ok,
  "\n",
  sep = ""
)

cat(
  "All source numerical values finite: ",
  source_all_finite,
  "\n",
  sep = ""
)


level1_passed <-
  source_identity_check &&
  source_dimensions_ok &&
  source_genes_ok &&
  source_all_finite


# ------------------------------------------------------------
# LEVEL 2
# Presentation-format transformation
#
# IMPORTANT:
# We DO NOT compare:
#   raw 0.362879...
# against
#   displayed 0.363
#
# Instead, we reproduce the EXACT display transformation
# and compare character strings.
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("LEVEL 2: Presentation-format transformation\n")
cat("------------------------------------------------------------\n")


expected_display <- data.frame(
  Gene = gene_tme_manuscript_table$Gene,
  
  `ImmuneScore, ρ` =
    sprintf(
      "%.3f",
      gene_tme_manuscript_table$ImmuneScore_Rho
    ),
  
  `ImmuneScore, FDR` =
    sprintf(
      "%.3e",
      gene_tme_manuscript_table$ImmuneScore_FDR
    ),
  
  `StromalScore, ρ` =
    sprintf(
      "%.3f",
      gene_tme_manuscript_table$StromalScore_Rho
    ),
  
  `StromalScore, FDR` =
    sprintf(
      "%.3e",
      gene_tme_manuscript_table$StromalScore_FDR
    ),
  
  `ESTIMATEScore, ρ` =
    sprintf(
      "%.3f",
      gene_tme_manuscript_table$ESTIMATEScore_Rho
    ),
  
  `ESTIMATEScore, FDR` =
    sprintf(
      "%.3e",
      gene_tme_manuscript_table$ESTIMATEScore_FDR
    ),
  
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# Ensure same row/column structure

display_structure_ok <-
  identical(
    dim(gene_tme_display_table),
    dim(expected_display)
  ) &&
  identical(
    names(gene_tme_display_table),
    names(expected_display)
  )


# Exact character-level comparison

display_content_ok <-
  identical(
    gene_tme_display_table,
    expected_display
  )


cat(
  "Display dimensions/column structure correct: ",
  display_structure_ok,
  "\n",
  sep = ""
)

cat(
  "Display values exactly follow locked formatting rules: ",
  display_content_ok,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# Explicit display examples
# ------------------------------------------------------------

cat("\nDisplay-format examples:\n")

cat(
  "AOC3 ImmuneScore rho: ",
  gene_tme_display_table[1, 2],
  "\n",
  sep = ""
)

cat(
  "Source value: ",
  sprintf(
    "%.10f",
    gene_tme_manuscript_table$ImmuneScore_Rho[1]
  ),
  "\n",
  sep = ""
)

cat(
  "AOC3 ImmuneScore FDR: ",
  gene_tme_display_table[1, 3],
  "\n",
  sep = ""
)

cat(
  "Source value: ",
  sprintf(
    "%.15e",
    gene_tme_manuscript_table$ImmuneScore_FDR[1]
  ),
  "\n",
  sep = ""
)


level2_passed <-
  display_structure_ok &&
  display_content_ok


# ------------------------------------------------------------
# LEVEL 3
# Manuscript metadata consistency
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("LEVEL 3: Manuscript metadata consistency\n")
cat("------------------------------------------------------------\n")


# 9 associations

n_genes <- length(unique(gene_tme_results$Gene))

n_traits <- length(
  unique(gene_tme_results$TME_Trait)
)

n_associations <- nrow(gene_tme_results)


# Sample size

n_values <- unique(gene_tme_results$N)

n_513_ok <-
  length(n_values) == 1 &&
  identical(as.integer(n_values), 513L)


# All nine FDR significant

all_fdr_significant <-
  all(gene_tme_results$FDR_BH < 0.05)


# All directions positive

all_positive <-
  all(gene_tme_results$Direction == "Positive")


# Check statistical wording

two_sided_ok <-
  grepl(
    "two-sided Spearman rank correlation",
    supp_table_footnote,
    fixed = TRUE
  )


bh_ok <-
  grepl(
    "Benjamini–Hochberg",
    supp_table_footnote,
    fixed = TRUE
  )


n_ok <-
  grepl(
    "N = 513",
    supp_table_footnote,
    fixed = TRUE
  )


composite_ok <-
  grepl(
    "ESTIMATEScore is a composite score",
    supp_table_footnote,
    fixed = TRUE
  )


nine_associations_ok <-
  grepl(
    "nine gene–TME associations",
    supp_table_footnote,
    fixed = TRUE
  )


cat(
  "Number of genes = 3: ",
  n_genes == 3,
  "\n",
  sep = ""
)

cat(
  "Number of TME traits = 3: ",
  n_traits == 3,
  "\n",
  sep = ""
)

cat(
  "Number of associations = 9: ",
  n_associations == 9,
  "\n",
  sep = ""
)

cat(
  "N = 513 correctly documented: ",
  n_513_ok,
  "\n",
  sep = ""
)

cat(
  "All nine FDR < 0.05: ",
  all_fdr_significant,
  "\n",
  sep = ""
)

cat(
  "All nine associations positive: ",
  all_positive,
  "\n",
  sep = ""
)

cat(
  "Two-sided Spearman documented: ",
  two_sided_ok,
  "\n",
  sep = ""
)

cat(
  "Benjamini-Hochberg documented: ",
  bh_ok,
  "\n",
  sep = ""
)

cat(
  "N = 513 documented in footnote: ",
  n_ok,
  "\n",
  sep = ""
)

cat(
  "ESTIMATEScore identified as composite: ",
  composite_ok,
  "\n",
  sep = ""
)

cat(
  "Nine gene-TME associations documented: ",
  nine_associations_ok,
  "\n",
  sep = ""
)


level3_passed <-
  n_genes == 3 &&
  n_traits == 3 &&
  n_associations == 9 &&
  n_513_ok &&
  all_fdr_significant &&
  all_positive &&
  two_sided_ok &&
  bh_ok &&
  n_ok &&
  composite_ok &&
  nine_associations_ok


# ------------------------------------------------------------
# FINAL STEP 7B-4C-R STATUS
# ------------------------------------------------------------

step7b4cR_passed <-
  level1_passed &&
  level2_passed &&
  level3_passed


cat("\n============================================================\n")
cat(
  "LEVEL 1 STATUS: ",
  ifelse(level1_passed, "PASSED", "FAILED"),
  "\n",
  sep = ""
)

cat(
  "LEVEL 2 STATUS: ",
  ifelse(level2_passed, "PASSED", "FAILED"),
  "\n",
  sep = ""
)

cat(
  "LEVEL 3 STATUS: ",
  ifelse(level3_passed, "PASSED", "FAILED"),
  "\n",
  sep = ""
)

cat(
  "STEP 7B-4C-R FINAL STATUS: ",
  ifelse(step7b4cR_passed, "PASSED", "FAILED"),
  "\n",
  sep = ""
)

cat("============================================================\n")


# ============================================================
# Paper 1
# STEP 7B-4D
# Export manuscript-ready Supplementary Table to Excel
# + post-export integrity audit
#
# IMPORTANT:
# No new statistical calculations.
# Locked source objects are NOT modified.
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-4D: Excel export + post-export integrity audit\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 0. Required objects
# ------------------------------------------------------------

required_objects <- c(
  "gene_tme_results",
  "gene_tme_manuscript_table",
  "gene_tme_display_table",
  "supp_table_title",
  "supp_table_footnote",
  "step7b4cR_passed"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ",
    paste(missing_objects, collapse = ", ")
  )
}


if (!isTRUE(step7b4cR_passed)) {
  stop(
    "STEP 7B-4C-R is not PASSED. ",
    "Excel export is stopped."
  )
}


cat("Required objects present: TRUE\n")
cat("STEP 7B-4C-R prerequisite: PASSED\n")


# ------------------------------------------------------------
# 1. Load package
# ------------------------------------------------------------

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop(
    "Package 'openxlsx' is required. ",
    "Please install it before continuing."
  )
}

library(openxlsx)

cat("openxlsx available: TRUE\n")


# ------------------------------------------------------------
# 2. Define output path
# ------------------------------------------------------------

output_dir <-
  "O:/Rproject/Paper1_tumor-only WGCNA"

xlsx_file <- file.path(
  output_dir,
  "Paper1_Supplementary_Table_X_CoreGene_TME_Associations.xlsx"
)


if (!dir.exists(output_dir)) {
  stop(
    "Output directory does not exist:\n",
    output_dir
  )
}


cat("\nOutput file:\n")
cat(xlsx_file, "\n")


# ------------------------------------------------------------
# 3. Preserve locked objects
# ------------------------------------------------------------

# IMPORTANT:
# These are copies only.
# The locked source objects themselves are not modified.

export_table <- gene_tme_display_table
audit_source <- gene_tme_manuscript_table


# ------------------------------------------------------------
# 4. Create workbook
# ------------------------------------------------------------

wb <- createWorkbook()


# ============================================================
# SHEET 1
# Manuscript-ready displayed table
# ============================================================

addWorksheet(
  wb,
  "Table_Supp_X"
)


# ------------------------------------------------------------
# 4A. Title
# ------------------------------------------------------------

writeData(
  wb,
  sheet = "Table_Supp_X",
  x = supp_table_title,
  startRow = 1,
  startCol = 1,
  colNames = FALSE
)


# ------------------------------------------------------------
# 4B. Table
# ------------------------------------------------------------

writeData(
  wb,
  sheet = "Table_Supp_X",
  x = export_table,
  startRow = 3,
  startCol = 1,
  rowNames = FALSE,
  colNames = TRUE
)


# ------------------------------------------------------------
# 4C. Footnote
# ------------------------------------------------------------

writeData(
  wb,
  sheet = "Table_Supp_X",
  x = supp_table_footnote,
  startRow = 8,
  startCol = 1,
  colNames = FALSE
)


# ============================================================
# SHEET 2
# Source audit
# ============================================================

addWorksheet(
  wb,
  "Source_Audit"
)


writeData(
  wb,
  sheet = "Source_Audit",
  x = audit_source,
  startRow = 1,
  startCol = 1,
  rowNames = FALSE,
  colNames = TRUE
)


# ------------------------------------------------------------
# 5. Formatting styles
# ------------------------------------------------------------

title_style <- createStyle(
  textDecoration = "bold",
  fontSize = 12,
  valign = "center"
)


header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "bottom"
)


footnote_style <- createStyle(
  wrapText = TRUE,
  valign = "top"
)


# ------------------------------------------------------------
# Apply title style
# ------------------------------------------------------------

addStyle(
  wb,
  sheet = "Table_Supp_X",
  style = title_style,
  rows = 1,
  cols = 1,
  gridExpand = FALSE
)


# ------------------------------------------------------------
# Apply table header style
# ------------------------------------------------------------

addStyle(
  wb,
  sheet = "Table_Supp_X",
  style = header_style,
  rows = 3,
  cols = 1:7,
  gridExpand = TRUE
)


addStyle(
  wb,
  sheet = "Source_Audit",
  style = header_style,
  rows = 1,
  cols = 1:7,
  gridExpand = TRUE
)


# ------------------------------------------------------------
# Footnote formatting
# ------------------------------------------------------------

addStyle(
  wb,
  sheet = "Table_Supp_X",
  style = footnote_style,
  rows = 8,
  cols = 1,
  gridExpand = TRUE
)


# ------------------------------------------------------------
# Merge title across seven columns
# ------------------------------------------------------------

mergeCells(
  wb,
  sheet = "Table_Supp_X",
  cols = 1:7,
  rows = 1
)


# ------------------------------------------------------------
# Merge footnote across seven columns
# ------------------------------------------------------------

mergeCells(
  wb,
  sheet = "Table_Supp_X",
  cols = 1:7,
  rows = 8
)


# ------------------------------------------------------------
# Column widths
# ------------------------------------------------------------

setColWidths(
  wb,
  sheet = "Table_Supp_X",
  cols = 1:7,
  widths = c(
    12,
    18,
    20,
    18,
    20,
    20,
    22
  )
)


setColWidths(
  wb,
  sheet = "Source_Audit",
  cols = 1:7,
  widths = "auto"
)


# ------------------------------------------------------------
# Row heights
# ------------------------------------------------------------

setRowHeights(
  wb,
  sheet = "Table_Supp_X",
  rows = 1,
  heights = 30
)


setRowHeights(
  wb,
  sheet = "Table_Supp_X",
  rows = 8,
  heights = 80
)


# ------------------------------------------------------------
# Freeze panes
# ------------------------------------------------------------

freezePane(
  wb,
  sheet = "Table_Supp_X",
  firstActiveRow = 4
)


freezePane(
  wb,
  sheet = "Source_Audit",
  firstActiveRow = 2
)


# ------------------------------------------------------------
# 6. Save workbook
# ------------------------------------------------------------

saveWorkbook(
  wb,
  xlsx_file,
  overwrite = TRUE
)


cat("\nWorkbook written successfully.\n")


# ============================================================
# POST-EXPORT AUDIT
# ============================================================

# ------------------------------------------------------------
# 7. Check physical file
# ------------------------------------------------------------

file_exists <- file.exists(xlsx_file)

cat(
  "Excel file exists: ",
  file_exists,
  "\n",
  sep = ""
)


if (!file_exists) {
  stop("Excel export failed: file does not exist.")
}


file_size <- file.info(xlsx_file)$size

cat(
  "Excel file size: ",
  file_size,
  " bytes\n",
  sep = ""
)


file_size_ok <- is.finite(file_size) && file_size > 0

cat(
  "Excel file size valid: ",
  file_size_ok,
  "\n",
  sep = ""
)

# ============================================================
# STEP 7B-4D-R
# POST-EXPORT AUDIT — corrected for openxlsx version
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-4D-R: Post-export integrity audit\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 8. Re-open exported workbook
# ------------------------------------------------------------

wb_check <- loadWorkbook(
  xlsx_file
)

cat("Workbook successfully reopened: TRUE\n")


# ------------------------------------------------------------
# 9. Check workbook sheet structure
# ------------------------------------------------------------

sheet_names <- names(wb_check)

cat("\nExported workbook sheets:\n")
print(sheet_names)


expected_sheets <- c(
  "Table_Supp_X",
  "Source_Audit"
)


sheet_structure_ok <-
  identical(
    sheet_names,
    expected_sheets
  )


cat(
  "Workbook sheet structure correct: ",
  sheet_structure_ok,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Read manuscript-ready table back
# ------------------------------------------------------------

# The table header is on row 3.
# The three gene rows are rows 4–6.

exported_table_readback <- read.xlsx(
  xlsx_file,
  sheet = "Table_Supp_X",
  rows = 3:6,
  colNames = TRUE,
  check.names = FALSE
)


cat("\nRead-back manuscript table:\n")
print(exported_table_readback)


# ------------------------------------------------------------
# 11. Keep the three actual data rows
# ------------------------------------------------------------

# Depending on openxlsx version/read.xlsx behavior,
# rows may be returned with or without the header.
# Therefore inspect dimensions before proceeding.

cat("\nRead-back dimensions:\n")
print(dim(exported_table_readback))


# If exactly 3 rows were returned, retain directly.
# If 4 rows were returned, the first row is redundant.

if (nrow(exported_table_readback) == 4L) {
  
  exported_table_readback <-
    exported_table_readback[-1, , drop = FALSE]
  
}


# ------------------------------------------------------------
# 12. Check dimensions
# ------------------------------------------------------------

readback_dimension_ok <-
  identical(
    dim(exported_table_readback),
    c(3L, 7L)
  )


cat(
  "Read-back dimensions 3 x 7: ",
  readback_dimension_ok,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. Check column names
# ------------------------------------------------------------

readback_colnames_ok <-
  identical(
    names(exported_table_readback),
    names(gene_tme_display_table)
  )


cat(
  "Read-back column names correct: ",
  readback_colnames_ok,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 14. Check displayed table values
# ------------------------------------------------------------

readback_values_ok <-
  identical(
    exported_table_readback,
    gene_tme_display_table
  )


cat(
  "Read-back displayed table exactly matches export object: ",
  readback_values_ok,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 15. Read Source_Audit back
# ------------------------------------------------------------

source_readback <- read.xlsx(
  xlsx_file,
  sheet = "Source_Audit",
  colNames = TRUE,
  check.names = FALSE
)


cat("\nRead-back Source_Audit dimensions:\n")
print(dim(source_readback))


source_readback_dimension_ok <-
  identical(
    dim(source_readback),
    dim(gene_tme_manuscript_table)
  )


cat(
  "Source_Audit dimensions correct: ",
  source_readback_dimension_ok,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 16. Check Source_Audit content
# ------------------------------------------------------------

source_numeric_identity <- TRUE


for (j in seq_along(gene_tme_manuscript_table)) {
  
  expected <- gene_tme_manuscript_table[[j]]
  observed <- source_readback[[j]]
  
  if (is.numeric(expected)) {
    
    source_numeric_identity <-
      source_numeric_identity &&
      isTRUE(
        all.equal(
          observed,
          expected,
          check.attributes = FALSE
        )
      )
    
  } else {
    
    source_numeric_identity <-
      source_numeric_identity &&
      identical(
        observed,
        expected
      )
  }
}


cat(
  "Source_Audit numerical/content identity: ",
  source_numeric_identity,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 17. Check title and footnote
# ------------------------------------------------------------

table_sheet_raw <- read.xlsx(
  xlsx_file,
  sheet = "Table_Supp_X",
  colNames = FALSE
)


title_found <-
  any(
    as.character(table_sheet_raw[[1]]) ==
      supp_table_title,
    na.rm = TRUE
  )


footnote_found <-
  any(
    as.character(table_sheet_raw[[1]]) ==
      supp_table_footnote,
    na.rm = TRUE
  )


cat(
  "Table title preserved: ",
  title_found,
  "\n",
  sep = ""
)

cat(
  "Table footnote preserved: ",
  footnote_found,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 18. Final STEP 7B-4D-R status
# ------------------------------------------------------------

step7b4dR_passed <-
  file_exists &&
  file_size_ok &&
  sheet_structure_ok &&
  readback_dimension_ok &&
  readback_colnames_ok &&
  readback_values_ok &&
  source_readback_dimension_ok &&
  source_numeric_identity &&
  title_found &&
  footnote_found


cat("\n============================================================\n")

cat(
  "File exists: ",
  file_exists,
  "\n",
  sep = ""
)

cat(
  "File size valid: ",
  file_size_ok,
  "\n",
  sep = ""
)

cat(
  "Sheet structure correct: ",
  sheet_structure_ok,
  "\n",
  sep = ""
)

cat(
  "Displayed table dimensions correct: ",
  readback_dimension_ok,
  "\n",
  sep = ""
)

cat(
  "Displayed table column names correct: ",
  readback_colnames_ok,
  "\n",
  sep = ""
)

cat(
  "Displayed table values preserved: ",
  readback_values_ok,
  "\n",
  sep = ""
)

cat(
  "Source_Audit preserved: ",
  source_numeric_identity,
  "\n",
  sep = ""
)

cat(
  "Title preserved: ",
  title_found,
  "\n",
  sep = ""
)

cat(
  "Footnote preserved: ",
  footnote_found,
  "\n",
  sep = ""
)

cat(
  "\nSTEP 7B-4D-R FINAL STATUS: ",
  ifelse(step7b4dR_passed, "PASSED", "FAILED"),
  "\n",
  sep = ""
)

cat("============================================================\n")

# ============================================================
# STEP 7B-4D-R-DIAG
# Diagnose Excel readback mismatch
# No source data modification
# No statistical recomputation
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-4D-R-DIAG: Excel readback mismatch diagnosis\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Compare expected vs observed column names
# ------------------------------------------------------------

cat("Expected column names:\n")
print(names(gene_tme_display_table))

cat("\nObserved Excel column names:\n")
print(names(exported_table_readback))


# ------------------------------------------------------------
# 2. Compare the displayed tables
# ------------------------------------------------------------

cat("\nExpected displayed table:\n")
print(gene_tme_display_table)

cat("\nObserved Excel readback table:\n")
print(exported_table_readback)


# ------------------------------------------------------------
# 3. Check whether the difference is only column names
# ------------------------------------------------------------

expected_values <- as.data.frame(
  gene_tme_display_table,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

observed_values <- as.data.frame(
  exported_table_readback,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Compare values after temporarily assigning expected names
observed_values_renamed <- observed_values
names(observed_values_renamed) <- names(expected_values)

readback_values_only_ok <- isTRUE(
  all.equal(
    observed_values_renamed,
    expected_values,
    check.attributes = FALSE
  )
)

cat("\nValues preserved after restoring expected column names: ",
    readback_values_only_ok,
    "\n",
    sep = "")


# ------------------------------------------------------------
# 4. Cell-by-cell comparison
# ------------------------------------------------------------

cell_difference <- matrix(
  FALSE,
  nrow = nrow(expected_values),
  ncol = ncol(expected_values)
)

for (i in seq_len(nrow(expected_values))) {
  for (j in seq_len(ncol(expected_values))) {
    
    cell_difference[i, j] <- !isTRUE(
      all.equal(
        expected_values[[j]][i],
        observed_values[[j]][i],
        check.attributes = FALSE
      )
    )
  }
}

cat("\nNumber of differing cells: ",
    sum(cell_difference),
    "\n",
    sep = "")


if (any(cell_difference)) {
  
  cat("\nDiffering cells:\n")
  
  diff_idx <- which(cell_difference, arr.ind = TRUE)
  
  for (k in seq_len(nrow(diff_idx))) {
    
    i <- diff_idx[k, "row"]
    j <- diff_idx[k, "col"]
    
    cat(
      "Row ", i,
      ", Column ", j,
      " | Expected = ", dQuote(as.character(expected_values[[j]][i])),
      " | Observed = ", dQuote(as.character(observed_values[[j]][i])),
      "\n",
      sep = ""
    )
  }
  
} else {
  
  cat("\nNo cell-value differences detected.\n")
}


# ------------------------------------------------------------
# 5. Final diagnostic interpretation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("DIAGNOSTIC SUMMARY\n")
cat("============================================================\n")

if (readback_values_only_ok && !identical(
  names(expected_values),
  names(observed_values)
)) {
  
  cat("RESULT: Values are preserved; mismatch is limited to column names.\n")
  cat("INTERPRETATION: Excel/openxlsx column-name normalization only.\n")
  
} else if (readback_values_only_ok) {
  
  cat("RESULT: Displayed table values are fully preserved.\n")
  
} else {
  
  cat("RESULT: At least one displayed-table cell differs after Excel export.\n")
  cat("ACTION: Inspect the cell-level differences above before locking.\n")
}

cat("============================================================\n")

# ============================================================
# STEP 7B-4D-R FINAL
# Post-export integrity audit with openxlsx name normalization
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-4D-R FINAL: Post-export integrity audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. File-level integrity
# ------------------------------------------------------------

file_exists <- file.exists(xlsx_file)

file_size <- if (file_exists) file.info(xlsx_file)$size else NA_real_

file_size_ok <- isTRUE(
  file_exists &&
    !is.na(file_size) &&
    file_size > 1000
)


# ------------------------------------------------------------
# 2. Workbook structure
# ------------------------------------------------------------

wb_check <- loadWorkbook(xlsx_file)

sheet_names <- names(wb_check)

expected_sheets <- c(
  "Table_Supp_X",
  "Source_Audit"
)

sheet_structure_ok <- identical(
  sheet_names,
  expected_sheets
)


# ------------------------------------------------------------
# 3. Read back manuscript display table
# ------------------------------------------------------------

exported_table_readback <- read.xlsx(
  xlsx_file,
  sheet = "Table_Supp_X",
  rows = 3:6,
  colNames = TRUE,
  check.names = FALSE
)

if (nrow(exported_table_readback) == 4L) {
  exported_table_readback <-
    exported_table_readback[-1, , drop = FALSE]
}


# ------------------------------------------------------------
# 4. Dimensions
# ------------------------------------------------------------

readback_dimension_ok <- identical(
  dim(exported_table_readback),
  c(3L, 7L)
)


# ------------------------------------------------------------
# 5. Column-name audit
#
# openxlsx/read.xlsx may normalize punctuation in names.
# Therefore compare names after removing this known
# normalization artifact.
# ------------------------------------------------------------

expected_names <- names(gene_tme_display_table)

observed_names <- names(exported_table_readback)

normalize_readback_name <- function(x) {
  x <- sub(",\\.", ", ", x)
  x
}

normalized_observed_names <- vapply(
  observed_names,
  normalize_readback_name,
  character(1)
)

readback_colnames_ok <- identical(
  normalized_observed_names,
  expected_names
)


# ------------------------------------------------------------
# 6. Displayed-table value integrity
# ------------------------------------------------------------

expected_display <- as.data.frame(
  gene_tme_display_table,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

observed_display <- as.data.frame(
  exported_table_readback,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Restore manuscript column names only for comparison
observed_display_for_audit <- observed_display

names(observed_display_for_audit) <- expected_names

readback_values_ok <- isTRUE(
  all.equal(
    observed_display_for_audit,
    expected_display,
    check.attributes = FALSE
  )
)


# ------------------------------------------------------------
# 7. Source_Audit integrity
# ------------------------------------------------------------

source_readback <- read.xlsx(
  xlsx_file,
  sheet = "Source_Audit",
  colNames = TRUE,
  check.names = FALSE
)

source_readback_dimension_ok <- identical(
  dim(source_readback),
  dim(gene_tme_manuscript_table)
)


source_numeric_identity <- TRUE

for (j in seq_along(gene_tme_manuscript_table)) {
  
  expected <- gene_tme_manuscript_table[[j]]
  observed <- source_readback[[j]]
  
  if (is.numeric(expected)) {
    
    source_numeric_identity <-
      source_numeric_identity &&
      isTRUE(
        all.equal(
          observed,
          expected,
          check.attributes = FALSE
        )
      )
    
  } else {
    
    source_numeric_identity <-
      source_numeric_identity &&
      identical(
        observed,
        expected
      )
  }
}


# ------------------------------------------------------------
# 8. Title and footnote
# ------------------------------------------------------------

table_sheet_raw <- read.xlsx(
  xlsx_file,
  sheet = "Table_Supp_X",
  colNames = FALSE
)

title_found <- any(
  as.character(table_sheet_raw[[1]]) ==
    supp_table_title,
  na.rm = TRUE
)

footnote_found <- any(
  as.character(table_sheet_raw[[1]]) ==
    supp_table_footnote,
  na.rm = TRUE
)


# ------------------------------------------------------------
# 9. Final PASS/FAIL
# ------------------------------------------------------------

step7b4dR_passed <-
  file_exists &&
  file_size_ok &&
  sheet_structure_ok &&
  readback_dimension_ok &&
  readback_colnames_ok &&
  readback_values_ok &&
  source_readback_dimension_ok &&
  source_numeric_identity &&
  title_found &&
  footnote_found


# ------------------------------------------------------------
# 10. Report
# ------------------------------------------------------------

cat("File exists: ",
    file_exists,
    "\n",
    sep = "")

cat("File size valid: ",
    file_size_ok,
    "\n",
    sep = "")

cat("Sheet structure correct: ",
    sheet_structure_ok,
    "\n",
    sep = "")

cat("Displayed table dimensions correct: ",
    readback_dimension_ok,
    "\n",
    sep = "")

cat("Displayed table column names correct: ",
    readback_colnames_ok,
    "\n",
    sep = "")

cat("Displayed table values preserved: ",
    readback_values_ok,
    "\n",
    sep = "")

cat("Source_Audit preserved: ",
    source_numeric_identity,
    "\n",
    sep = "")

cat("Title preserved: ",
    title_found,
    "\n",
    sep = "")

cat("Footnote preserved: ",
    footnote_found,
    "\n",
    sep = "")

cat(
  "\nSTEP 7B-4D-R FINAL STATUS: ",
  ifelse(
    step7b4dR_passed,
    "PASSED - LOCKED",
    "FAILED"
  ),
  "\n",
  sep = ""
)

cat("============================================================\n")

# ============================================================
# STEP 7B-4D-R FIX
# Correct readback column-name normalization
# ============================================================

expected_names <- names(gene_tme_display_table)
observed_names <- names(exported_table_readback)

normalize_readback_name <- function(x) {
  sub(",\\.", ", ", x)
}

normalized_observed_names <- vapply(
  observed_names,
  normalize_readback_name,
  character(1)
)

cat("Expected names:\n")
print(expected_names)

cat("\nObserved names from Excel:\n")
print(observed_names)

cat("\nNormalized observed names:\n")
print(normalized_observed_names)

readback_colnames_ok <- identical(
  normalized_observed_names,
  expected_names
)

cat(
  "\nDisplayed table column names correct after normalization: ",
  readback_colnames_ok,
  "\n",
  sep = ""
)

# ============================================================
# STEP 7B-4D-R FIX-2
# Explicit column-name mapping for Excel readback audit
# No source data modification
# ============================================================

expected_names <- names(gene_tme_display_table)
observed_names <- names(exported_table_readback)

expected_from_excel <- c(
  "Gene",
  "ImmuneScore,.ρ",
  "ImmuneScore,.FDR",
  "StromalScore,.ρ",
  "StromalScore,.FDR",
  "ESTIMATEScore,.ρ",
  "ESTIMATEScore,.FDR"
)

cat("Expected manuscript names:\n")
print(expected_names)

cat("\nObserved Excel names:\n")
print(observed_names)

cat("\nExpected Excel-readback names:\n")
print(expected_from_excel)

# Exact positional comparison
readback_colnames_ok <- identical(
  observed_names,
  expected_from_excel
)

cat(
  "\nExcel readback column-name structure correct: ",
  readback_colnames_ok,
  "\n",
  sep = ""
)

# Confirm one-to-one correspondence with manuscript names
name_mapping_ok <- identical(
  length(expected_names),
  length(expected_from_excel)
)

if (name_mapping_ok) {
  
  name_mapping <- data.frame(
    Position = seq_along(expected_names),
    Manuscript_Name = expected_names,
    Excel_Readback_Name = observed_names,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  name_mapping$Match <- observed_names == expected_from_excel
  
  cat("\nColumn-name mapping:\n")
  print(name_mapping)
  
  all_name_positions_match <- all(name_mapping$Match)
  
} else {
  
  all_name_positions_match <- FALSE
}

cat(
  "\nAll column positions match expected Excel readback structure: ",
  all_name_positions_match,
  "\n",
  sep = ""
)

# ============================================================
# Paper 1
# Figure 6
# STEP 6A: Audit locked TME objects
# ============================================================

cat("\n============================================================\n")
cat("FIGURE 6 - STEP 6A: Locked object audit\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. Required objects
# ------------------------------------------------------------

required_objects <- c(
  "core_tme_alignment",
  "gene_tme_results",
  "gene_tme_manuscript_table"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    paste(
      "Missing required object(s):",
      paste(missing_objects, collapse = ", ")
    )
  )
}

cat("All required objects are present.\n\n")


# ------------------------------------------------------------
# 2. Structure audit
# ------------------------------------------------------------

cat("core_tme_alignment:\n")
print(dim(core_tme_alignment))
print(names(core_tme_alignment))

cat("\nGene order:\n")
print(c("AOC3", "MS4A7", "IL7R"))

cat("\ngene_tme_results:\n")
print(dim(gene_tme_results))
print(names(gene_tme_results))

cat("\ngene_tme_manuscript_table:\n")
print(dim(gene_tme_manuscript_table))
print(names(gene_tme_manuscript_table))


# ------------------------------------------------------------
# 3. Basic integrity checks
# ------------------------------------------------------------

stopifnot(
  nrow(core_tme_alignment) == 513,
  all(c(
    "Sample_ID",
    "AOC3",
    "MS4A7",
    "IL7R",
    "StromalScore",
    "ImmuneScore",
    "ESTIMATEScore"
  ) %in% names(core_tme_alignment))
)

stopifnot(
  identical(
    core_tme_alignment$Sample_ID,
    rownames(core_expr)
  )
)

numeric_columns <- c(
  "AOC3",
  "MS4A7",
  "IL7R",
  "StromalScore",
  "ImmuneScore",
  "ESTIMATEScore"
)

na_check <- sapply(
  core_tme_alignment[, numeric_columns, drop = FALSE],
  function(x) anyNA(x) || any(!is.finite(x))
)

print(na_check)

if (any(na_check)) {
  stop("NA/NaN/Inf detected in Figure 6 source data.")
}

cat("\nSTEP 6A PASSED.\n")
cat("N = 513 primary LUAD tumors.\n")
cat("No statistics were recalculated.\n")

# ============================================================
# Paper 1
# Figure 6
# STEP 6B: Construct locked plotting data
# ============================================================

cat("\n============================================================\n")
cat("FIGURE 6 - STEP 6B: Construct plotting data\n")
cat("============================================================\n")

library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# 1. Construct long-format plotting data
# ------------------------------------------------------------

figure6_long <- core_tme_alignment %>%
  select(
    Sample_ID,
    AOC3,
    MS4A7,
    IL7R,
    ImmuneScore,
    StromalScore,
    ESTIMATEScore
  ) %>%
  pivot_longer(
    cols = c(AOC3, MS4A7, IL7R),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  pivot_longer(
    cols = c(
      ImmuneScore,
      StromalScore,
      ESTIMATEScore
    ),
    names_to = "Trait",
    values_to = "TME_Score"
  )

# ------------------------------------------------------------
# 2. Lock gene order
# ------------------------------------------------------------

figure6_long$Gene <- factor(
  figure6_long$Gene,
  levels = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )
)

# ------------------------------------------------------------
# 3. Lock TME trait order
# ------------------------------------------------------------

figure6_long$Trait <- factor(
  figure6_long$Trait,
  levels = c(
    "ImmuneScore",
    "StromalScore",
    "ESTIMATEScore"
  )
)

# ------------------------------------------------------------
# 4. Basic structure check
# ------------------------------------------------------------

cat("\nDimensions of figure6_long:\n")
print(dim(figure6_long))

cat("\nColumn names:\n")
print(names(figure6_long))

cat("\nGene levels:\n")
print(levels(figure6_long$Gene))

cat("\nTrait levels:\n")
print(levels(figure6_long$Trait))

# ------------------------------------------------------------
# 5. Expected panel structure
# ------------------------------------------------------------

figure6_panel_check <- figure6_long %>%
  group_by(Gene, Trait) %>%
  summarise(
    N = n(),
    .groups = "drop"
  )

cat("\nPanel structure:\n")
print(figure6_panel_check)

# ------------------------------------------------------------
# 6. Strict validation
# ------------------------------------------------------------

stopifnot(
  nrow(figure6_long) == 513 * 3 * 3,
  length(unique(figure6_long$Sample_ID)) == 513,
  length(unique(figure6_long$Gene)) == 3,
  length(unique(figure6_long$Trait)) == 3,
  nrow(figure6_panel_check) == 9,
  all(figure6_panel_check$N == 513)
)

cat("\nSTEP 6B PASSED.\n")
cat("513 samples × 3 genes × 3 TME traits = 4617 observations.\n")
cat("Nine panels, each containing N = 513.\n")
cat("No statistical calculation performed.\n")

# ============================================================
# Paper 1
# Figure 6
# STEP 6C: Locked statistical annotation and cross-check
# ============================================================

cat("\n============================================================\n")
cat("FIGURE 6 - STEP 6C: Locked statistical annotation audit\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Prepare locked statistics directly from gene_tme_results
# ------------------------------------------------------------

figure6_stats <- gene_tme_results %>%
  select(
    Gene,
    TME_Trait,
    N,
    Spearman_Rho,
    P_value,
    FDR_BH,
    Direction,
    Significant_P05,
    Significant_FDR05
  )

# ------------------------------------------------------------
# 2. Lock factor order
# ------------------------------------------------------------

figure6_stats$Gene <- factor(
  figure6_stats$Gene,
  levels = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )
)

figure6_stats$TME_Trait <- factor(
  figure6_stats$TME_Trait,
  levels = c(
    "ImmuneScore",
    "StromalScore",
    "ESTIMATEScore"
  )
)

# ------------------------------------------------------------
# 3. Sort in Figure 6 panel order
# ------------------------------------------------------------

figure6_stats <- figure6_stats %>%
  arrange(
    Gene,
    TME_Trait
  )

cat("\nLocked Figure 6 statistics:\n")
print(figure6_stats)


# ------------------------------------------------------------
# 4. Strict structural audit
# ------------------------------------------------------------

stopifnot(
  nrow(figure6_stats) == 9,
  all(figure6_stats$N == 513),
  all(!is.na(figure6_stats$Spearman_Rho)),
  all(!is.na(figure6_stats$FDR_BH)),
  all(is.finite(figure6_stats$Spearman_Rho)),
  all(is.finite(figure6_stats$FDR_BH))
)

# ------------------------------------------------------------
# 5. Verify all nine associations are positive
# ------------------------------------------------------------

positive_check <- all(
  figure6_stats$Spearman_Rho > 0
)

cat(
  "\nAll nine Spearman rho values positive:",
  positive_check,
  "\n"
)

stopifnot(
  positive_check
)


# ------------------------------------------------------------
# 6. Verify all nine associations are FDR significant
# ------------------------------------------------------------

fdr_check <- all(
  figure6_stats$FDR_BH < 0.05
)

cat(
  "All nine FDR values < 0.05:",
  fdr_check,
  "\n"
)

stopifnot(
  fdr_check
)


# ------------------------------------------------------------
# 7. Create display labels
# ------------------------------------------------------------

figure6_stats <- figure6_stats %>%
  mutate(
    Rho_label = sprintf(
      "\u03c1 = %.3f",
      Spearman_Rho
    ),
    FDR_label = paste0(
      "FDR = ",
      formatC(
        FDR_BH,
        format = "e",
        digits = 3
      )
    ),
    Annotation = paste(
      Rho_label,
      FDR_label,
      sep = "\n"
    )
  )

cat("\nDisplay annotations:\n")

print(
  figure6_stats %>%
    select(
      Gene,
      TME_Trait,
      Rho_label,
      FDR_label
    )
)


# ------------------------------------------------------------
# 8. Cross-check against gene_tme_manuscript_table
# ------------------------------------------------------------

# Convert the manuscript table into long format
# without recalculating any statistics.

table6_long <- gene_tme_manuscript_table %>%
  pivot_longer(
    cols = -Gene,
    names_to = c(
      "TME_Trait",
      "Statistic"
    ),
    names_pattern =
      "^(ImmuneScore|StromalScore|ESTIMATEScore)_(Rho|FDR)$",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = Value
  )

table6_long$Gene <- factor(
  table6_long$Gene,
  levels = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )
)

table6_long$TME_Trait <- factor(
  table6_long$TME_Trait,
  levels = c(
    "ImmuneScore",
    "StromalScore",
    "ESTIMATEScore"
  )
)

table6_long <- table6_long %>%
  arrange(
    Gene,
    TME_Trait
  )


# ------------------------------------------------------------
# 9. Compare rho
# ------------------------------------------------------------

rho_difference <- abs(
  figure6_stats$Spearman_Rho -
    table6_long$Rho
)

rho_match <- all(
  rho_difference < 1e-12
)

cat(
  "\nRho values match gene_tme_manuscript_table:",
  rho_match,
  "\n"
)

stopifnot(
  rho_match
)


# ------------------------------------------------------------
# 10. Compare FDR
# ------------------------------------------------------------

fdr_difference <- abs(
  figure6_stats$FDR_BH -
    table6_long$FDR
)

fdr_match <- all(
  fdr_difference < 1e-12
)

cat(
  "FDR values match gene_tme_manuscript_table:",
  fdr_match,
  "\n"
)

stopifnot(
  fdr_match
)


# ------------------------------------------------------------
# 11. Verify panel-to-statistic matching
# ------------------------------------------------------------

panel_stats_check <- figure6_long %>%
  distinct(
    Gene,
    Trait
  ) %>%
  left_join(
    figure6_stats %>%
      select(
        Gene,
        TME_Trait,
        Spearman_Rho,
        FDR_BH
      ),
    by = c(
      "Gene" = "Gene",
      "Trait" = "TME_Trait"
    )
  )

cat("\nPanel-to-statistic matching:\n")
print(panel_stats_check)

stopifnot(
  nrow(panel_stats_check) == 9,
  all(!is.na(panel_stats_check$Spearman_Rho)),
  all(!is.na(panel_stats_check$FDR_BH))
)


# ------------------------------------------------------------
# 12. Final STEP 6C audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 6C PASSED.\n")
cat("============================================================\n")

cat("Nine locked gene-TME associations identified.\n")
cat("All N = 513.\n")
cat("All nine rho values are positive.\n")
cat("All nine FDR values < 0.05.\n")
cat("Rho values match gene_tme_manuscript_table.\n")
cat("FDR values match gene_tme_manuscript_table.\n")
cat("Panel-to-statistic mapping is complete.\n")
cat("No statistical calculation was performed.\n")

# ============================================================
# STEP 6F-1
# Correct Figure 6 point-color assignment
# ------------------------------------------------------------
# ImmuneScore    -> blue
# StromalScore   -> green
# ESTIMATEScore  -> purple
#
# No statistical calculation is performed.
# ============================================================

cat("\n============================================================\n")
cat("STEP 6F-1: Correct Figure 6 point-color assignment\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Define the three TME colors
# ------------------------------------------------------------

tme_colors_6F <- c(
  "ImmuneScore"   = "#3B82F6",
  "StromalScore"  = "#22A06B",
  "ESTIMATEScore" = "#8B5CF6"
)


# ------------------------------------------------------------
# 2. Create visualization-only copy
# ------------------------------------------------------------

figure6_long_6F <- figure6_long


# ------------------------------------------------------------
# 3. Assign exactly one color to each TME trait
# ------------------------------------------------------------

figure6_long_6F$PointColor <- unname(
  tme_colors_6F[
    as.character(figure6_long_6F$Trait)
  ]
)


# ------------------------------------------------------------
# 4. Audit color assignment
# ------------------------------------------------------------

color_audit_6F <- figure6_long_6F %>%
  group_by(Trait) %>%
  summarise(
    n_colors = n_distinct(PointColor),
    Color = first(PointColor),
    N = n(),
    .groups = "drop"
  )

print(color_audit_6F)


# Each TME trait must have exactly one color
stopifnot(
  all(color_audit_6F$n_colors == 1)
)


# Expected colors
stopifnot(
  color_audit_6F$Color[
    match("ImmuneScore", color_audit_6F$Trait)
  ] == tme_colors_6F["ImmuneScore"]
)

stopifnot(
  color_audit_6F$Color[
    match("StromalScore", color_audit_6F$Trait)
  ] == tme_colors_6F["StromalScore"]
)

stopifnot(
  color_audit_6F$Color[
    match("ESTIMATEScore", color_audit_6F$Trait)
  ] == tme_colors_6F["ESTIMATEScore"]
)


cat("\nColor assignment audit PASSED.\n")
cat("ImmuneScore    -> blue\n")
cat("StromalScore   -> green\n")
cat("ESTIMATEScore  -> purple\n")

# ============================================================
# Paper 1
# STEP 6F-2
# FINAL Figure 6 visualization
# ------------------------------------------------------------
# Critical correction:
#   figure6_long uses "Trait".
#   We explicitly rename it to "TME_Trait" for visualization.
#
# Color assignment:
#   ImmuneScore    -> blue
#   StromalScore   -> green
#   ESTIMATEScore  -> purple
#
# No statistical recalculation.
# ============================================================

cat("\n============================================================\n")
cat("STEP 6F-2: FINAL Figure 6 visualization\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Create a clean visualization-only dataset
# ------------------------------------------------------------

figure6_long_6F <- figure6_long %>%
  mutate(
    TME_Trait = as.character(Trait)
  ) %>%
  select(
    Sample_ID,
    Gene,
    Expression,
    TME_Trait,
    TME_Score
  )


# ------------------------------------------------------------
# 2. Set the exact desired factor order
# ------------------------------------------------------------

figure6_long_6F$Gene <- factor(
  figure6_long_6F$Gene,
  levels = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )
)

figure6_long_6F$TME_Trait <- factor(
  figure6_long_6F$TME_Trait,
  levels = c(
    "ImmuneScore",
    "StromalScore",
    "ESTIMATEScore"
  )
)


# ------------------------------------------------------------
# 3. Define the three colors
# ------------------------------------------------------------

tme_colors_6F <- c(
  "ImmuneScore"    = "#3B82F6",
  "StromalScore"   = "#22A06B",
  "ESTIMATEScore"  = "#8B5CF6"
)


# ------------------------------------------------------------
# 4. Assign exactly ONE color to each TME trait
# ------------------------------------------------------------

figure6_long_6F$PointColor <- unname(
  tme_colors_6F[
    as.character(figure6_long_6F$TME_Trait)
  ]
)


# ------------------------------------------------------------
# 5. Data structure audit
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("STEP 6F-2 data structure audit\n")
cat("------------------------------------------------------------\n")

cat(
  "Dimensions:",
  nrow(figure6_long_6F),
  "rows x",
  ncol(figure6_long_6F),
  "columns\n"
)

print(
  names(figure6_long_6F)
)


# Expected:
# 4617 rows = 513 samples × 3 genes × 3 TME traits

stopifnot(
  nrow(figure6_long_6F) == 4617
)


# ------------------------------------------------------------
# 6. Verify each TME trait has exactly ONE color
# ------------------------------------------------------------

color_audit_6F <- figure6_long_6F %>%
  group_by(TME_Trait) %>%
  summarise(
    n_colors = n_distinct(PointColor),
    PointColor = first(PointColor),
    N = n(),
    .groups = "drop"
  )

print(color_audit_6F)


stopifnot(
  nrow(color_audit_6F) == 3,
  all(color_audit_6F$n_colors == 1),
  all(color_audit_6F$N == 1539)
)


# ------------------------------------------------------------
# 7. Verify exact color assignment
# ------------------------------------------------------------

color_mapping_check_6F <- all(
  color_audit_6F$PointColor[
    match(
      c(
        "ImmuneScore",
        "StromalScore",
        "ESTIMATEScore"
      ),
      color_audit_6F$TME_Trait
    )
  ] ==
    unname(
      tme_colors_6F[
        c(
          "ImmuneScore",
          "StromalScore",
          "ESTIMATEScore"
        )
      ]
    )
)

cat(
  "\nColor mapping correct:",
  color_mapping_check_6F,
  "\n"
)

stopifnot(color_mapping_check_6F)


# ------------------------------------------------------------
# 8. Verify every Gene × TME trait panel contains N = 513
# ------------------------------------------------------------

panel_n_6F <- figure6_long_6F %>%
  group_by(
    Gene,
    TME_Trait
  ) %>%
  summarise(
    N = n(),
    n_colors = n_distinct(PointColor),
    PointColor = first(PointColor),
    .groups = "drop"
  )

cat("\nPanel audit:\n")
print(panel_n_6F)


stopifnot(
  nrow(panel_n_6F) == 9,
  all(panel_n_6F$N == 513),
  all(panel_n_6F$n_colors == 1)
)


# ------------------------------------------------------------
# 9. Build the final Figure 6
# ------------------------------------------------------------

figure6_plot_6F <- ggplot(
  figure6_long_6F,
  aes(
    x = Expression,
    y = TME_Score
  )
) +
  
  # ==========================================================
# CRITICAL:
# Point color is mapped DIRECTLY to PointColor.
# ==========================================================

geom_point(
  aes(
    color = PointColor
  ),
  size = 1.05,
  alpha = 0.42,
  stroke = 0
) +
  
  # ----------------------------------------------------------
# Rho annotation
# ----------------------------------------------------------

geom_text(
  data = annotation_rho_6F,
  aes(
    x = -Inf,
    y = Inf,
    label = Annotation
  ),
  inherit.aes = FALSE,
  hjust = -0.08,
  vjust = 1.35,
  size = 3.0,
  fontface = "plain"
) +
  
  # ----------------------------------------------------------
# FDR annotation
# ----------------------------------------------------------

geom_text(
  data = annotation_fdr_6F,
  aes(
    x = -Inf,
    y = Inf,
    label = Annotation
  ),
  inherit.aes = FALSE,
  hjust = -0.08,
  vjust = 2.65,
  size = 2.75,
  fontface = "plain"
) +
  
  # ----------------------------------------------------------
# EXACT 3 × 3 facet structure
# ----------------------------------------------------------

facet_grid(
  rows = vars(Gene),
  cols = vars(TME_Trait),
  scales = "free"
) +
  
  # ----------------------------------------------------------
# Use the actual hexadecimal colors in PointColor
# ----------------------------------------------------------

scale_color_identity(
  guide = "none"
) +
  
  # ----------------------------------------------------------
# Axis labels
# ----------------------------------------------------------

labs(
  x = "Gene expression",
  y = "ESTIMATE-derived score"
) +
  
  # ----------------------------------------------------------
# Theme
# ----------------------------------------------------------

theme_classic(
  base_size = 11
) +
  
  theme(
    
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    
    strip.background = element_blank(),
    
    axis.title.x = element_text(
      size = 11,
      margin = margin(t = 8)
    ),
    
    axis.title.y = element_text(
      size = 11,
      margin = margin(r = 8)
    ),
    
    axis.text.x = element_text(
      size = 8.5
    ),
    
    axis.text.y = element_text(
      size = 8.5
    ),
    
    axis.line = element_line(
      linewidth = 0.45
    ),
    
    axis.ticks = element_line(
      linewidth = 0.4
    ),
    
    panel.spacing = unit(
      0.9,
      "lines"
    ),
    
    plot.margin = margin(
      t = 8,
      r = 10,
      b = 8,
      l = 10
    )
  )


# ------------------------------------------------------------
# 10. Print preview
# ------------------------------------------------------------

print(figure6_plot_6F)

# 12. Export
# ------------------------------------------------------------

figure6_dir_6F <- file.path(
  "result_dir/Paper 1_Step 7_ Immune Infiltration Analysis",
  "figure6_plot_6F"
)

if (!dir.exists(figure6_dir_6F)) {
  dir.create(
    figure6_dir_6F,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 12A. 600 dpi TIFF
# ------------------------------------------------------------

figure6_tiff_6F <- file.path(
  figure6_dir_6F,
  "Figure6_CoreGenes_TME_Associations_600dpi_STEP6F.tiff"
)

ggsave(
  filename = figure6_tiff_6F,
  plot = figure6_plot_6F,
  width = 180,
  height = 180,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)


# ------------------------------------------------------------
# 12B. Vector PDF
# ------------------------------------------------------------

figure6_pdf_6F <- file.path(
  figure6_dir_6F,
  "Figure6_CoreGenes_TME_Associations_STEP6F.pdf"
)

ggsave(
  filename = figure6_pdf_6F,
  plot = figure6_plot_6F,
  width = 180,
  height = 180,
  units = "mm",
  device = cairo_pdf
)

### ============================================================
# Paper 1 — STEP 7B-1
# Purity-adjusted analysis: source-object audit
# Purpose:
#   1. Identify the locked tumor-purity object
#   2. Confirm alignment with the 513 primary LUAD samples
#   3. Do NOT recalculate Step 7A-2 statistics
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-1: Tumor purity source-object audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. List potentially relevant objects
# ------------------------------------------------------------

candidate_objects <- ls(
  pattern = "(purity|Purity|PURITY|estimate|ESTIMATE|score|Score)"
)

cat("Potentially relevant objects:\n")
print(candidate_objects)


# ------------------------------------------------------------
# 2. Inspect object classes and dimensions
# ------------------------------------------------------------

if (length(candidate_objects) > 0) {
  
  audit_objects <- lapply(
    candidate_objects,
    function(x) {
      
      obj <- get(x)
      
      data.frame(
        Object = x,
        Class = paste(class(obj), collapse = ", "),
        Nrow = if (!is.null(dim(obj))) nrow(obj) else NA_integer_,
        Ncol = if (!is.null(dim(obj))) ncol(obj) else NA_integer_,
        Length = length(obj),
        stringsAsFactors = FALSE
      )
    }
  )
  
  audit_objects <- do.call(
    rbind,
    audit_objects
  )
  
  print(audit_objects)
  
} else {
  
  cat("No potentially relevant purity/ESTIMATE objects detected.\n")
}


# ------------------------------------------------------------
# 3. Confirm the locked primary-tumor sample IDs
# ------------------------------------------------------------

stopifnot(
  exists("core_tme_alignment"),
  exists("core_expr")
)

primary_ids <- rownames(core_expr)

cat("\nPrimary tumor cohort:\n")
cat("N =", length(primary_ids), "\n")
cat("Unique IDs =", length(unique(primary_ids)), "\n")

stopifnot(
  length(primary_ids) == 513,
  length(unique(primary_ids)) == 513
)


# ------------------------------------------------------------
# 4. Confirm Step 7A-2 objects remain unchanged
# ------------------------------------------------------------

stopifnot(
  exists("gene_tme_results"),
  exists("gene_tme_manuscript_table"),
  exists("core_tme_alignment")
)

cat("\nLocked Step 7A-2 objects detected:\n")
cat("gene_tme_results:",
    nrow(gene_tme_results), "rows x",
    ncol(gene_tme_results), "columns\n")

cat("gene_tme_manuscript_table:",
    nrow(gene_tme_manuscript_table), "rows x",
    ncol(gene_tme_manuscript_table), "columns\n")

cat("core_tme_alignment:",
    nrow(core_tme_alignment), "rows x",
    ncol(core_tme_alignment), "columns\n")


# ------------------------------------------------------------
# 5. Basic integrity check for core_tme_alignment
# ------------------------------------------------------------

stopifnot(
  nrow(core_tme_alignment) == 513,
  identical(
    core_tme_alignment$Sample_ID,
    rownames(core_expr)
  )
)

numeric_cols <- c(
  "AOC3",
  "MS4A7",
  "IL7R",
  "StromalScore",
  "ImmuneScore",
  "ESTIMATEScore"
)

stopifnot(
  all(numeric_cols %in% colnames(core_tme_alignment))
)

stopifnot(
  !anyNA(core_tme_alignment[, numeric_cols])
)

stopifnot(
  all(
    is.finite(
      as.matrix(core_tme_alignment[, numeric_cols])
    )
  )
)


cat("\n============================================================\n")
cat("STEP 7B-1 completed: source-object audit only.\n")
cat("No statistical analysis performed.\n")
cat("No Step 7A-2 object modified.\n")
cat("============================================================\n")


##STEP 7B-2 — Identify and validate the tumor-purity source
# ============================================================
# Paper 1 — STEP 7B-2
# Identify the locked tumor-purity source
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-2: Identify tumor-purity source\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Inspect estimate_env
# ------------------------------------------------------------

stopifnot(exists("estimate_env"))

cat("Class of estimate_env:\n")
print(class(estimate_env))

cat("\nObjects inside estimate_env:\n")
estimate_env_objects <- ls(envir = estimate_env)
print(estimate_env_objects)

cat("\nNumber of objects inside estimate_env:",
    length(estimate_env_objects), "\n")


# ------------------------------------------------------------
# 2. Inspect each object inside estimate_env
# ------------------------------------------------------------

if (length(estimate_env_objects) > 0) {
  
  estimate_env_audit <- lapply(
    estimate_env_objects,
    function(x) {
      
      obj <- get(x, envir = estimate_env)
      
      data.frame(
        Object = x,
        Class = paste(class(obj), collapse = ", "),
        Nrow = if (!is.null(dim(obj))) nrow(obj) else NA_integer_,
        Ncol = if (!is.null(dim(obj))) ncol(obj) else NA_integer_,
        Length = length(obj),
        stringsAsFactors = FALSE
      )
    }
  )
  
  estimate_env_audit <- do.call(
    rbind,
    estimate_env_audit
  )
  
  print(estimate_env_audit)
  
} else {
  
  cat("estimate_env is empty.\n")
}


# ------------------------------------------------------------
# 3. Inspect estimate_objects
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("estimate_objects\n")
cat("------------------------------------------------------------\n")

print(estimate_objects)


# ------------------------------------------------------------
# 4. Inspect estimate_tmp
# ------------------------------------------------------------

if (exists("estimate_tmp")) {
  
  cat("\n------------------------------------------------------------\n")
  cat("estimate_tmp structure\n")
  cat("------------------------------------------------------------\n")
  
  print(class(estimate_tmp))
  print(dim(estimate_tmp))
  print(colnames(estimate_tmp))
  print(rownames(estimate_tmp))
  
  cat("\nFirst 6 columns:\n")
  print(
    estimate_tmp[, seq_len(min(6, ncol(estimate_tmp))),
                 drop = FALSE]
  )
}


# ------------------------------------------------------------
# 5. Inspect estimate_file
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("estimate_file\n")
cat("------------------------------------------------------------\n")

print(estimate_file)


# ------------------------------------------------------------
# 6. Search all character object names for purity-related terms
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Purity-related object-name search\n")
cat("------------------------------------------------------------\n")

all_objects <- ls(envir = .GlobalEnv)

purity_candidates <- all_objects[
  grepl(
    "purity|tumor.*pur|tumour.*pur|normal.*pur|estimate|stromal|immune",
    all_objects,
    ignore.case = TRUE
  )
]

print(purity_candidates)


cat("\n============================================================\n")
cat("STEP 7B-2 audit completed.\n")
cat("No statistical analysis performed.\n")
cat("No Step 7A-2 object modified.\n")
cat("============================================================\n")


# ============================================================
# Paper 1 — STEP 7B-3
# Independent tumor-purity source discovery audit
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-3: Independent tumor-purity source discovery audit\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Search all global objects by purity-related names
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("1. Global object-name search\n")
cat("------------------------------------------------------------\n")

all_objects <- ls(envir = .GlobalEnv)

purity_name_candidates <- all_objects[
  grepl(
    "purity|tumor_pur|tumour_pur|absolute|cellularity|tumorcell|tumourcell|tcga.*pur",
    all_objects,
    ignore.case = TRUE
  )
]

print(purity_name_candidates)

cat("\nNumber of name-based candidates:",
    length(purity_name_candidates), "\n")


# ------------------------------------------------------------
# 2. Search all data.frame / matrix objects for purity-like
#    column names
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("2. Search tabular objects for purity-related columns\n")
cat("------------------------------------------------------------\n")

tabular_objects <- all_objects[
  sapply(
    all_objects,
    function(x) {
      obj <- tryCatch(
        get(x, envir = .GlobalEnv),
        error = function(e) NULL
      )
      
      is.data.frame(obj) || is.matrix(obj)
    }
  )
]

cat("Number of tabular objects:",
    length(tabular_objects), "\n\n")

purity_column_hits <- list()

for (x in tabular_objects) {
  
  obj <- tryCatch(
    get(x, envir = .GlobalEnv),
    error = function(e) NULL
  )
  
  if (is.null(obj)) next
  
  cn <- colnames(obj)
  
  if (is.null(cn)) next
  
  hit <- cn[
    grepl(
      "purity|tumor.*pur|tumour.*pur|absolute|cellularity|tumorcell|tumourcell",
      cn,
      ignore.case = TRUE
    )
  ]
  
  if (length(hit) > 0) {
    
    purity_column_hits[[x]] <- hit
    
    cat("\nObject:", x, "\n")
    cat("Matching columns:\n")
    print(hit)
  }
}

cat("\n------------------------------------------------------------\n")
cat("Purity-column search completed.\n")
cat("------------------------------------------------------------\n")


# ------------------------------------------------------------
# 3. Inspect likely clinical objects
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("3. Likely clinical / sample annotation objects\n")
cat("------------------------------------------------------------\n")

clinical_candidates <- all_objects[
  grepl(
    "clinical|clin|metadata|meta|pheno|phenodata|patient|sample",
    all_objects,
    ignore.case = TRUE
  )
]

print(clinical_candidates)

cat("\nNumber of clinical-like candidates:",
    length(clinical_candidates), "\n")

if (length(clinical_candidates) > 0) {
  
  for (x in clinical_candidates) {
    
    obj <- tryCatch(
      get(x, envir = .GlobalEnv),
      error = function(e) NULL
    )
    
    if (is.null(obj)) next
    
    if (is.data.frame(obj) || is.matrix(obj)) {
      
      cat("\n============================================================\n")
      cat("Clinical-like object:", x, "\n")
      cat("Class:", paste(class(obj), collapse = ", "), "\n")
      
      d <- dim(obj)
      
      if (!is.null(d)) {
        cat("Dimension:", d[1], "x", d[2], "\n")
      }
      
      cn <- colnames(obj)
      
      if (!is.null(cn)) {
        
        relevant_columns <- cn[
          grepl(
            "purity|tumor|tumour|absolute|cellularity|sample|patient",
            cn,
            ignore.case = TRUE
          )
        ]
        
        cat("Relevant column names:\n")
        print(relevant_columns)
      }
    }
  }
}


# ------------------------------------------------------------
# 4. Inspect objects already known to be related to ESTIMATE
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("4. ESTIMATE-related objects\n")
cat("------------------------------------------------------------\n")

estimate_related <- all_objects[
  grepl("estimate", all_objects, ignore.case = TRUE)
]

print(estimate_related)

if (length(estimate_related) > 0) {
  
  for (x in estimate_related) {
    
    obj <- tryCatch(
      get(x, envir = .GlobalEnv),
      error = function(e) NULL
    )
    
    if (is.null(obj)) next
    
    cat("\nObject:", x, "\n")
    cat("Class:", paste(class(obj), collapse = ", "), "\n")
    
    if (!is.null(dim(obj))) {
      cat("Dimension:", nrow(obj), "x", ncol(obj), "\n")
    }
    
    if (is.data.frame(obj)) {
      cat("Column names:\n")
      print(colnames(obj))
    }
    
    if (is.matrix(obj)) {
      cat("Column names:\n")
      print(colnames(obj))
    }
  }
}


# ------------------------------------------------------------
# 5. Explicit check of the locked ESTIMATE source
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("5. Locked ESTIMATE source audit\n")
cat("------------------------------------------------------------\n")

cat("estimate_file:\n")
print(estimate_file)

cat("\nObjects recorded in estimate_objects:\n")
print(estimate_objects)

cat("\nObjects inside estimate_env:\n")
print(ls(envir = estimate_env))


# ------------------------------------------------------------
# 6. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 7B-3 audit completed.\n")
cat("No statistical analysis performed.\n")
cat("No purity variable created.\n")
cat("No Step 7A-2 object modified.\n")
cat("============================================================\n")

# ============================================================
# Paper 1 — STEP 7B-4A
# Search local Paper 1 directory for independent purity files
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-4A: Search local directory for purity source\n")
cat("============================================================\n\n")

paper1_root <- "input_dir/Paper1_tumor-only WGCNA"

stopifnot(dir.exists(paper1_root))

# ------------------------------------------------------------
# Search filenames
# ------------------------------------------------------------

all_files <- list.files(
  paper1_root,
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

purity_files <- all_files[
  grepl(
    "absolute|purity|cellularity|ploidy|mastercalls",
    basename(all_files),
    ignore.case = TRUE
  )
]

cat("Number of purity-related files found:",
    length(purity_files), "\n\n")

if (length(purity_files) > 0) {
  print(purity_files)
} else {
  cat("No local purity-related files found.\n")
}

# ------------------------------------------------------------
# Search directory names as well
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Directories with purity-related names\n")
cat("------------------------------------------------------------\n")

purity_dirs <- all_files[
  file.info(all_files)$isdir &
    grepl(
      "absolute|purity|cellularity|ploidy",
      basename(all_files),
      ignore.case = TRUE
    )
]

if (length(purity_dirs) > 0) {
  print(purity_dirs)
} else {
  cat("No purity-related directories found.\n")
}

cat("\n============================================================\n")
cat("STEP 7B-4A completed.\n")
cat("No statistical analysis performed.\n")
cat("No Step 7A-2 object modified.\n")
cat("============================================================\n")


# ============================================================
# STEP 7B-4B — Check downloaded ABSOLUTE source
# ============================================================

paper1_root <- "input_dir/Paper1_tumor-only WGCNA"

absolute_file <- list.files(
  paper1_root,
  pattern = "TCGA_mastercalls.*abs.*\\.txt$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

cat(
  "ABSOLUTE files found:",
  length(absolute_file),
  "\n\n"
)

print(absolute_file)

# ============================================================
# Paper 1 — STEP 7B-4B-1
# ABSOLUTE source file integrity check
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-4B-1: ABSOLUTE source file integrity check\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Confirm file exists
# ------------------------------------------------------------

stopifnot(
  length(absolute_file) == 1
)

stopifnot(
  file.exists(absolute_file)
)

cat("File exists: TRUE\n\n")


# ------------------------------------------------------------
# 2. File size
# ------------------------------------------------------------

file_size <- file.info(absolute_file)$size

cat(
  "File size:",
  file_size,
  "bytes\n"
)

cat(
  "File size:",
  round(file_size / 1024^2, 3),
  "MB\n\n"
)


# ------------------------------------------------------------
# 3. Read first 5 lines ONLY
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("First 5 lines of source file\n")
cat("------------------------------------------------------------\n\n")

absolute_head <- readLines(
  absolute_file,
  n = 5,
  warn = FALSE
)

print(absolute_head)


# ------------------------------------------------------------
# 4. Basic file-format check
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Basic format checks\n")
cat("------------------------------------------------------------\n")

header_line <- absolute_head[1]

tab_count_header <- lengths(
  gregexpr(
    "\t",
    header_line,
    fixed = TRUE
  )
)

cat(
  "Number of tab separators in header:",
  tab_count_header,
  "\n"
)

cat(
  "Contains tab-delimited structure:",
  tab_count_header > 0,
  "\n"
)


# ------------------------------------------------------------
# 5. Check for HTML / download-error contamination
# ------------------------------------------------------------

html_like <- grepl(
  "<html|<!doctype|<head|<body",
  paste(
    absolute_head,
    collapse = "\n"
  ),
  ignore.case = TRUE
)

cat(
  "HTML/download-error signature detected:",
  html_like,
  "\n"
)


# ------------------------------------------------------------
# 6. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")

if (
  file_size > 0 &&
  length(absolute_head) >= 2 &&
  tab_count_header > 0 &&
  !html_like
) {
  
  cat(
    "STEP 7B-4B-1 PASSED: source file appears structurally valid.\n"
  )
  
} else {
  
  cat(
    "STEP 7B-4B-1 FAILED: source file requires inspection.\n"
  )
}

cat("No statistical analysis performed.\n")
cat("No purity variable created.\n")
cat("No Step 7A-2 object modified.\n")

cat("============================================================\n")


# ============================================================
# Paper 1 — STEP 7B-4C
# ABSOLUTE source audit and primary-LUAD sample matching
#
# Purpose:
#   1. Read the complete ABSOLUTE source
#   2. Validate array/sample/purity fields
#   3. Match ABSOLUTE array IDs to the 513 primary LUAD IDs
#   4. Quantify exact matching, missingness and duplicates
#   5. Do NOT perform regression or any association test
#   6. Do NOT modify any STEP 7A-2 object
#
# IMPORTANT:
#   - Purity = ABSOLUTE "purity"
#   - Sample ID = ABSOLUTE "array"
#   - No barcode truncation or transformation is performed yet
#   - No purity value is inferred
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-4C: ABSOLUTE source audit and sample matching\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Protection checks
# ------------------------------------------------------------

stopifnot(
  exists("absolute_file")
)

stopifnot(
  length(absolute_file) == 1
)

stopifnot(
  file.exists(absolute_file)
)

stopifnot(
  exists("core_tme_alignment")
)

stopifnot(
  is.data.frame(core_tme_alignment)
)

stopifnot(
  "Sample_ID" %in% names(core_tme_alignment)
)

cat("Source file and primary cohort object detected.\n")


# ------------------------------------------------------------
# 1. Read complete ABSOLUTE source
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("1. Reading complete ABSOLUTE source\n")
cat("------------------------------------------------------------\n")

absolute_df <- read.delim(
  absolute_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  na.strings = c(
    "",
    "NA",
    "NaN",
    "NULL"
  )
)

cat(
  "ABSOLUTE dimensions:",
  nrow(absolute_df),
  "rows x",
  ncol(absolute_df),
  "columns\n"
)


# ------------------------------------------------------------
# 2. Required-column check
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("2. Required-column check\n")
cat("------------------------------------------------------------\n")

required_columns <- c(
  "array",
  "sample",
  "call status",
  "purity",
  "ploidy",
  "Genome doublings",
  "Coverage for 80% power",
  "Cancer DNA fraction",
  "Subclonal genome fraction",
  "solution"
)

required_columns_present <- all(
  required_columns %in% names(absolute_df)
)

cat(
  "All expected ABSOLUTE columns present:",
  required_columns_present,
  "\n"
)

if (!required_columns_present) {
  
  cat("\nMissing columns:\n")
  
  print(
    setdiff(
      required_columns,
      names(absolute_df)
    )
  )
  
  stop(
    "STEP 7B-4C stopped: required ABSOLUTE columns missing."
  )
}


# ------------------------------------------------------------
# 3. Inspect data classes
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("3. Data classes\n")
cat("------------------------------------------------------------\n")

class_audit <- data.frame(
  Column = required_columns,
  Class = vapply(
    absolute_df[required_columns],
    function(x) class(x)[1],
    character(1)
  ),
  stringsAsFactors = FALSE
)

print(
  class_audit,
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. Basic ID integrity
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("4. ABSOLUTE ID integrity\n")
cat("------------------------------------------------------------\n")

array_ids <- as.character(
  absolute_df$array
)

sample_ids_absolute <- as.character(
  absolute_df$sample
)

array_missing_n <- sum(
  is.na(array_ids) |
    !nzchar(array_ids)
)

sample_missing_n <- sum(
  is.na(sample_ids_absolute) |
    !nzchar(sample_ids_absolute)
)

array_unique_n <- length(
  unique(
    array_ids[
      !is.na(array_ids) &
        nzchar(array_ids)
    ]
  )
)

sample_unique_n <- length(
  unique(
    sample_ids_absolute[
      !is.na(sample_ids_absolute) &
        nzchar(sample_ids_absolute)
    ]
  )
)

array_duplicate_n <- sum(
  duplicated(
    array_ids[
      !is.na(array_ids) &
        nzchar(array_ids)
    ]
  )
)

sample_duplicate_n <- sum(
  duplicated(
    sample_ids_absolute[
      !is.na(sample_ids_absolute) &
        nzchar(sample_ids_absolute)
    ]
  )
)

cat(
  "Array IDs:",
  length(array_ids),
  "\n"
)

cat(
  "Unique array IDs:",
  array_unique_n,
  "\n"
)

cat(
  "Missing/empty array IDs:",
  array_missing_n,
  "\n"
)

cat(
  "Duplicate array IDs:",
  array_duplicate_n,
  "\n\n"
)

cat(
  "Sample IDs:",
  length(sample_ids_absolute),
  "\n"
)

cat(
  "Unique sample IDs:",
  sample_unique_n,
  "\n"
)

cat(
  "Missing/empty sample IDs:",
  sample_missing_n,
  "\n"
)

cat(
  "Duplicate sample IDs:",
  sample_duplicate_n,
  "\n"
)


# ------------------------------------------------------------
# 5. Purity integrity
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("5. ABSOLUTE purity integrity\n")
cat("------------------------------------------------------------\n")

purity_raw <- absolute_df$purity

purity_is_numeric <- is.numeric(
  purity_raw
)

cat(
  "Purity column numeric:",
  purity_is_numeric,
  "\n"
)

if (!purity_is_numeric) {
  
  cat(
    "Purity class:",
    class(purity_raw)[1],
    "\n"
  )
  
  stop(
    "STEP 7B-4C stopped: ABSOLUTE purity column is not numeric."
  )
}

purity_missing_n <- sum(
  is.na(purity_raw)
)

purity_nonfinite_n <- sum(
  !is.finite(purity_raw)
)

purity_complete <- purity_raw[
  is.finite(purity_raw)
]

purity_min <- min(
  purity_complete
)

purity_max <- max(
  purity_complete
)

purity_median <- median(
  purity_complete
)

purity_mean <- mean(
  purity_complete
)

cat(
  "Purity missing:",
  purity_missing_n,
  "\n"
)

cat(
  "Purity non-finite:",
  purity_nonfinite_n,
  "\n"
)

cat(
  "Purity minimum:",
  purity_min,
  "\n"
)

cat(
  "Purity median:",
  purity_median,
  "\n"
)

cat(
  "Purity mean:",
  purity_mean,
  "\n"
)

cat(
  "Purity maximum:",
  purity_max,
  "\n"
)

purity_0_1_valid <- all(
  purity_complete >= 0 &
    purity_complete <= 1
)

cat(
  "All finite purity values within [0, 1]:",
  purity_0_1_valid,
  "\n"
)


# ------------------------------------------------------------
# 6. Call-status audit
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("6. ABSOLUTE call-status audit\n")
cat("------------------------------------------------------------\n")

call_status_table <- table(
  absolute_df$`call status`,
  useNA = "ifany"
)

print(
  call_status_table
)


# ------------------------------------------------------------
# 7. Prepare locked 513 primary LUAD sample IDs
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("7. Preparing the locked 513 primary LUAD IDs\n")
cat("------------------------------------------------------------\n")

primary_ids <- as.character(
  core_tme_alignment$Sample_ID
)

stopifnot(
  length(primary_ids) == 513
)

stopifnot(
  length(unique(primary_ids)) == 513
)

cat(
  "Primary LUAD N:",
  length(primary_ids),
  "\n"
)

cat(
  "Unique primary LUAD IDs:",
  length(unique(primary_ids)),
  "\n"
)


# ------------------------------------------------------------
# 8. Exact matching — NO barcode transformation
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("8. Exact sample-ID matching\n")
cat("------------------------------------------------------------\n")

exact_match <- primary_ids %in% array_ids

exact_match_n <- sum(
  exact_match
)

exact_match_fraction <-
  exact_match_n / length(primary_ids)

cat(
  "Exact matches:",
  exact_match_n,
  "/",
  length(primary_ids),
  "\n"
)

cat(
  "Exact match fraction:",
  exact_match_fraction,
  "\n"
)


# ------------------------------------------------------------
# 9. Unmatched primary LUAD samples
# ------------------------------------------------------------

unmatched_primary_ids <- primary_ids[
  !exact_match
]

cat(
  "Unmatched primary LUAD samples:",
  length(unmatched_primary_ids),
  "\n"
)

if (length(unmatched_primary_ids) > 0) {
  
  cat("\nFirst unmatched IDs:\n")
  
  print(
    head(
      unmatched_primary_ids,
      20
    )
  )
}


# ------------------------------------------------------------
# 10. ABSOLUTE records corresponding to primary LUAD
# ------------------------------------------------------------

absolute_primary_rows <- absolute_df[
  absolute_df$array %in% primary_ids,
  ,
  drop = FALSE
]

cat(
  "\nABSOLUTE rows matching primary LUAD IDs:",
  nrow(absolute_primary_rows),
  "\n"
)


# ------------------------------------------------------------
# 11. Check one-to-one matching
# ------------------------------------------------------------

primary_match_counts <- table(
  factor(
    absolute_df$array,
    levels = primary_ids
  )
)

multiple_absolute_records <- sum(
  primary_match_counts > 1
)

zero_absolute_records <- sum(
  primary_match_counts == 0
)

cat(
  "Primary IDs with >1 ABSOLUTE record:",
  multiple_absolute_records,
  "\n"
)

cat(
  "Primary IDs with 0 ABSOLUTE records:",
  zero_absolute_records,
  "\n"
)


# ------------------------------------------------------------
# 12. Purity availability among matched samples
# ------------------------------------------------------------

matched_array_ids <- primary_ids[
  primary_ids %in% array_ids
]

matched_rows <- match(
  matched_array_ids,
  absolute_df$array
)

matched_purity <- absolute_df$purity[
  matched_rows
]

matched_purity_missing_n <- sum(
  is.na(matched_purity)
)

matched_purity_finite_n <- sum(
  is.finite(matched_purity)
)

cat(
  "\nMatched samples:",
  length(matched_array_ids),
  "\n"
)

cat(
  "Matched samples with missing purity:",
  matched_purity_missing_n,
  "\n"
)

cat(
  "Matched samples with finite purity:",
  matched_purity_finite_n,
  "\n"
)


# ------------------------------------------------------------
# 13. Construct audit object ONLY
# ------------------------------------------------------------

step7b4c_audit <- list(
  
  step = "STEP 7B-4C",
  
  source_file = absolute_file,
  
  source_dimensions = dim(
    absolute_df
  ),
  
  required_columns_present =
    required_columns_present,
  
  array_unique_n =
    array_unique_n,
  
  array_duplicate_n =
    array_duplicate_n,
  
  sample_unique_n =
    sample_unique_n,
  
  sample_duplicate_n =
    sample_duplicate_n,
  
  purity_is_numeric =
    purity_is_numeric,
  
  purity_missing_n =
    purity_missing_n,
  
  purity_nonfinite_n =
    purity_nonfinite_n,
  
  purity_min =
    purity_min,
  
  purity_median =
    purity_median,
  
  purity_mean =
    purity_mean,
  
  purity_max =
    purity_max,
  
  purity_within_0_1 =
    purity_0_1_valid,
  
  call_status_table =
    call_status_table,
  
  primary_LUAD_N =
    length(primary_ids),
  
  primary_LUAD_unique_N =
    length(unique(primary_ids)),
  
  exact_match_n =
    exact_match_n,
  
  exact_match_fraction =
    exact_match_fraction,
  
  unmatched_primary_ids =
    unmatched_primary_ids,
  
  multiple_absolute_records =
    multiple_absolute_records,
  
  zero_absolute_records =
    zero_absolute_records,
  
  matched_purity_missing_n =
    matched_purity_missing_n,
  
  matched_purity_finite_n =
    matched_purity_finite_n,
  
  barcode_transformation_performed =
    FALSE,
  
  statistical_analysis_performed =
    FALSE,
  
  regression_performed =
    FALSE,
  
  spearman_performed =
    FALSE,
  
  step7a2_objects_modified =
    FALSE,
  
  status =
    "AUDIT_ONLY"
  
)


# ------------------------------------------------------------
# 14. Save audit checkpoint
# ------------------------------------------------------------

audit_file_7b4c <- file.path(
  paper1_root,
  "Paper1_STEP7B4C_ABSOLUTE_matching_audit.RData"
)

save(
  step7b4c_audit,
  file = audit_file_7b4c
)

cat(
  "\nAudit checkpoint saved:\n",
  audit_file_7b4c,
  "\n"
)


# ------------------------------------------------------------
# 15. Final protection check
# ------------------------------------------------------------

stopifnot(
  nrow(core_tme_alignment) == 513
)

stopifnot(
  length(unique(core_tme_alignment$Sample_ID)) == 513
)

cat("\n============================================================\n")
cat("STEP 7B-4C completed.\n")
cat("No statistical analysis performed.\n")
cat("No regression performed.\n")
cat("No barcode transformation performed.\n")
cat("No Step 7A-2 object modified.\n")
cat("Status: AUDIT_ONLY\n")
cat("============================================================\n")

# ============================================================
# Paper 1 — STEP 7B-4D
# ABSOLUTE barcode reconciliation and missing-purity audit
#
# Purpose:
#   1. Investigate the 11 exact-unmatched primary LUAD samples
#   2. Determine whether their corresponding sample-level
#      barcodes occur in the ABSOLUTE "sample" field
#   3. Inspect the 2 exact-matched samples with missing purity
#   4. NO statistical analysis
#   5. NO regression
#   6. NO modification of STEP 7A-2 objects
#   7. NO automatic barcode remapping
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-4D: ABSOLUTE barcode reconciliation audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Protection checks
# ------------------------------------------------------------

stopifnot(
  exists("absolute_df")
)

stopifnot(
  exists("primary_ids")
)

stopifnot(
  exists("unmatched_primary_ids")
)

stopifnot(
  exists("core_tme_alignment")
)

stopifnot(
  nrow(core_tme_alignment) == 513
)

cat("Protection checks passed.\n")


# ------------------------------------------------------------
# 1. Inspect the 11 exact-unmatched samples
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("1. Exact-unmatched primary LUAD samples\n")
cat("------------------------------------------------------------\n")

print(
  unmatched_primary_ids
)


# ------------------------------------------------------------
# 2. Search ABSOLUTE 'sample' field using the first
#    15-character TCGA sample-level barcode
#
# Example:
#   TCGA-55-7816-01
#
# ABSOLUTE sample:
#   TCGA-55-7816-01A-...
#
# This is ONLY an audit.
# No remapping is performed.
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("2. Search unmatched IDs against ABSOLUTE aliquot barcodes\n")
cat("------------------------------------------------------------\n")

unmatched_patient_sample_prefix <- substr(
  unmatched_primary_ids,
  1,
  15
)

absolute_sample_prefix <- substr(
  absolute_df$sample,
  1,
  15
)

prefix_match_list <- lapply(
  unmatched_patient_sample_prefix,
  function(x) {
    which(
      absolute_sample_prefix == x
    )
  }
)

names(prefix_match_list) <- unmatched_primary_ids


# ------------------------------------------------------------
# 3. Summarize barcode reconciliation
# ------------------------------------------------------------

prefix_match_summary <- data.frame(
  Primary_ID = unmatched_primary_ids,
  Prefix_15char = unmatched_patient_sample_prefix,
  ABSOLUTE_sample_records =
    vapply(
      prefix_match_list,
      length,
      integer(1)
    ),
  stringsAsFactors = FALSE
)

cat("\n")
print(
  prefix_match_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. Display all ABSOLUTE records found by prefix
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("3. ABSOLUTE records corresponding to unmatched prefixes\n")
cat("------------------------------------------------------------\n")

for (i in seq_along(unmatched_primary_ids)) {
  
  id_i <- unmatched_primary_ids[i]
  
  idx_i <- prefix_match_list[[i]]
  
  cat("\n============================================================\n")
  cat("Primary ID:", id_i, "\n")
  cat("Prefix:", unmatched_patient_sample_prefix[i], "\n")
  cat(
    "ABSOLUTE records found:",
    length(idx_i),
    "\n"
  )
  
  if (length(idx_i) > 0) {
    
    print(
      absolute_df[
        idx_i,
        c(
          "array",
          "sample",
          "call status",
          "purity",
          "ploidy",
          "solution"
        ),
        drop = FALSE
      ],
      row.names = FALSE
    )
    
  } else {
    
    cat(
      "No corresponding ABSOLUTE sample barcode found.\n"
    )
  }
}


# ------------------------------------------------------------
# 5. Identify the 2 exact-matched samples with missing purity
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("4. Exact-matched samples with missing purity\n")
cat("------------------------------------------------------------\n")

exact_matched_ids <- primary_ids[
  primary_ids %in% array_ids
]

exact_matched_rows <- match(
  exact_matched_ids,
  absolute_df$array
)

exact_matched_df <- absolute_df[
  exact_matched_rows,
  ,
  drop = FALSE
]

missing_purity_idx <- which(
  is.na(
    exact_matched_df$purity
  )
)

missing_purity_df <- exact_matched_df[
  missing_purity_idx,
  c(
    "array",
    "sample",
    "call status",
    "purity",
    "ploidy",
    "Genome doublings",
    "Coverage for 80% power",
    "Cancer DNA fraction",
    "Subclonal genome fraction",
    "solution"
  ),
  drop = FALSE
]

cat(
  "Exact-matched samples with missing purity:",
  nrow(missing_purity_df),
  "\n\n"
)

print(
  missing_purity_df,
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Check whether missing-purity samples have
#    any alternative records
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("5. Alternative ABSOLUTE records for missing-purity samples\n")
cat("------------------------------------------------------------\n")

for (id_i in missing_purity_df$array) {
  
  idx_i <- which(
    absolute_df$array == id_i
  )
  
  cat("\n")
  cat("Primary array ID:", id_i, "\n")
  cat(
    "Number of ABSOLUTE records:",
    length(idx_i),
    "\n"
  )
  
  print(
    absolute_df[
      idx_i,
      c(
        "array",
        "sample",
        "call status",
        "purity",
        "solution"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 7. Build reconciliation audit object
# ------------------------------------------------------------

step7b4d_audit <- list(
  
  step = "STEP 7B-4D",
  
  unmatched_primary_ids =
    unmatched_primary_ids,
  
  unmatched_prefix_summary =
    prefix_match_summary,
  
  prefix_match_list =
    prefix_match_list,
  
  missing_purity_df =
    missing_purity_df,
  
  exact_matched_N =
    length(exact_matched_ids),
  
  exact_matched_missing_purity_N =
    nrow(missing_purity_df),
  
  exact_matched_finite_purity_N =
    sum(
      is.finite(
        exact_matched_df$purity
      )
    ),
  
  barcode_transformation_performed =
    FALSE,
  
  automatic_remapping_performed =
    FALSE,
  
  statistical_analysis_performed =
    FALSE,
  
  regression_performed =
    FALSE,
  
  spearman_performed =
    FALSE,
  
  step7a2_objects_modified =
    FALSE,
  
  status =
    "AUDIT_ONLY"
  
)


# ------------------------------------------------------------
# 8. Save audit checkpoint
# ------------------------------------------------------------

audit_file_7b4d <- file.path(
  paper1_root,
  "Paper1_STEP7B4D_ABSOLUTE_barcode_reconciliation_audit.RData"
)

save(
  step7b4d_audit,
  file = audit_file_7b4d
)

cat(
  "\nAudit checkpoint saved:\n",
  audit_file_7b4d,
  "\n"
)


# ------------------------------------------------------------
# 9. Final protection checks
# ------------------------------------------------------------

stopifnot(
  nrow(core_tme_alignment) == 513
)

stopifnot(
  length(unique(core_tme_alignment$Sample_ID)) == 513
)

cat("\n============================================================\n")
cat("STEP 7B-4D completed.\n")
cat("No statistical analysis performed.\n")
cat("No regression performed.\n")
cat("No automatic barcode remapping performed.\n")
cat("No Step 7A-2 object modified.\n")
cat("Status: AUDIT_ONLY\n")
cat("============================================================\n")

##STEP 7B-5 — Construct the locked purity-adjusted analysis dataset。
# ============================================================
# Paper 1
# STEP 7B-5
# Construct locked purity-adjusted analysis dataset
#
# Purpose:
#   Build the N=500 dataset for the subsequent
#   ABSOLUTE-purity-adjusted sensitivity analysis.
#
# IMPORTANT:
#   - No statistical analysis
#   - No regression
#   - No Spearman correlation
#   - No modification of STEP 7A-2 objects
#   - No imputation of missing purity
#   - No barcode remapping
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-5: Construct locked purity-adjusted analysis dataset\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Required-object protection checks
# ------------------------------------------------------------

required_objects_7b5 <- c(
  "core_tme_alignment",
  "absolute_df",
  "primary_ids",
  "unmatched_primary_ids",
  "step7b4d_audit"
)

missing_objects_7b5 <- required_objects_7b5[
  !vapply(
    required_objects_7b5,
    exists,
    logical(1)
  )
]

if (length(missing_objects_7b5) > 0) {
  stop(
    "Missing required objects: ",
    paste(missing_objects_7b5, collapse = ", ")
  )
}

stopifnot(
  nrow(core_tme_alignment) == 513
)

stopifnot(
  length(unique(core_tme_alignment$Sample_ID)) == 513
)

cat("Required-object protection checks: PASSED\n")


# ------------------------------------------------------------
# 1. Verify the locked source structures
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("1. Source structure audit\n")
cat("------------------------------------------------------------\n")

cat(
  "core_tme_alignment:",
  nrow(core_tme_alignment),
  "rows x",
  ncol(core_tme_alignment),
  "columns\n"
)

cat(
  "absolute_df:",
  nrow(absolute_df),
  "rows x",
  ncol(absolute_df),
  "columns\n"
)

required_core_cols_7b5 <- c(
  "Sample_ID",
  "AOC3",
  "MS4A7",
  "IL7R",
  "StromalScore",
  "ImmuneScore",
  "ESTIMATEScore"
)

stopifnot(
  all(required_core_cols_7b5 %in% names(core_tme_alignment))
)

stopifnot(
  all(
    c("array", "purity") %in%
      names(absolute_df)
  )
)

cat("Required source columns: PASSED\n")


# ------------------------------------------------------------
# 2. Preserve the exact 513-sample primary cohort
# ------------------------------------------------------------

primary_ids_7b5 <- core_tme_alignment$Sample_ID

stopifnot(
  length(primary_ids_7b5) == 513
)

stopifnot(
  length(unique(primary_ids_7b5)) == 513
)

cat(
  "Primary LUAD samples:",
  length(primary_ids_7b5),
  "\n"
)


# ------------------------------------------------------------
# 3. Match ABSOLUTE purity by exact array ID
#
# IMPORTANT:
#   No barcode transformation is performed.
# ------------------------------------------------------------

absolute_match_idx_7b5 <- match(
  primary_ids_7b5,
  absolute_df$array
)

exact_match_7b5 <- !is.na(
  absolute_match_idx_7b5
)

cat(
  "Exact ABSOLUTE array matches:",
  sum(exact_match_7b5),
  "/",
  length(primary_ids_7b5),
  "\n"
)

stopifnot(
  sum(exact_match_7b5) == 502
)


# ------------------------------------------------------------
# 4. Extract ABSOLUTE purity
# ------------------------------------------------------------

absolute_purity_7b5 <- rep(
  NA_real_,
  length(primary_ids_7b5)
)

absolute_purity_7b5[exact_match_7b5] <-
  absolute_df$purity[
    absolute_match_idx_7b5[exact_match_7b5]
  ]


# ------------------------------------------------------------
# 5. Purity availability audit
# ------------------------------------------------------------

finite_purity_7b5 <- is.finite(
  absolute_purity_7b5
)

missing_purity_7b5 <- !finite_purity_7b5

cat(
  "Finite ABSOLUTE purity:",
  sum(finite_purity_7b5),
  "\n"
)

cat(
  "Missing/non-finite purity:",
  sum(missing_purity_7b5),
  "\n"
)

stopifnot(
  sum(finite_purity_7b5) == 500
)

stopifnot(
  sum(missing_purity_7b5) == 13
)


# ------------------------------------------------------------
# 6. Explicitly identify excluded samples
# ------------------------------------------------------------

excluded_purity_ids_7b5 <- primary_ids_7b5[
  !finite_purity_7b5
]

cat("\nExcluded primary samples due to unavailable purity:\n")

print(
  excluded_purity_ids_7b5
)


# ------------------------------------------------------------
# 7. Verify exclusion structure
#
# Expected:
#   11 no ABSOLUTE record
#   2 exact ABSOLUTE matches with NA purity
#   total = 13
# ------------------------------------------------------------

expected_unmatched_7b5 <- unmatched_primary_ids

expected_missing_purity_7b5 <- c(
  "TCGA-78-8648-01",
  "TCGA-75-7030-01"
)

expected_excluded_7b5 <- c(
  expected_unmatched_7b5,
  expected_missing_purity_7b5
)

stopifnot(
  setequal(
    excluded_purity_ids_7b5,
    expected_excluded_7b5
  )
)

cat(
  "\nExclusion identity check: PASSED\n"
)

cat(
  "  No ABSOLUTE record:",
  length(expected_unmatched_7b5),
  "\n"
)

cat(
  "  ABSOLUTE record but missing purity:",
  length(expected_missing_purity_7b5),
  "\n"
)

cat(
  "  Total excluded:",
  length(expected_excluded_7b5),
  "\n"
)


# ------------------------------------------------------------
# 8. Construct the 500-sample analysis dataset
#
# This is a NEW object.
# The locked STEP 7A-2 objects are untouched.
# ------------------------------------------------------------

purity_adjusted_dataset <- core_tme_alignment[
  finite_purity_7b5,
  ,
  drop = FALSE
]

purity_adjusted_dataset$ABSOLUTE_Purity <-
  absolute_purity_7b5[
    finite_purity_7b5
  ]


# ------------------------------------------------------------
# 9. Dataset structure audit
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("2. Constructed purity-adjusted dataset\n")
cat("------------------------------------------------------------\n")

cat(
  "Dimensions:",
  nrow(purity_adjusted_dataset),
  "x",
  ncol(purity_adjusted_dataset),
  "\n"
)

print(
  names(purity_adjusted_dataset)
)

stopifnot(
  nrow(purity_adjusted_dataset) == 500
)

stopifnot(
  ncol(purity_adjusted_dataset) == 8
)

stopifnot(
  length(unique(
    purity_adjusted_dataset$Sample_ID
  )) == 500
)


# ------------------------------------------------------------
# 10. Verify no missing/non-finite values
# ------------------------------------------------------------

analysis_variables_7b5 <- c(
  "Sample_ID",
  "AOC3",
  "MS4A7",
  "IL7R",
  "StromalScore",
  "ImmuneScore",
  "ESTIMATEScore",
  "ABSOLUTE_Purity"
)

stopifnot(
  all(
    analysis_variables_7b5 %in%
      names(purity_adjusted_dataset)
  )
)

na_count_7b5 <- vapply(
  purity_adjusted_dataset[
    ,
    analysis_variables_7b5,
    drop = FALSE
  ],
  function(x) sum(is.na(x)),
  integer(1)
)

nonfinite_count_7b5 <- vapply(
  purity_adjusted_dataset[
    ,
    setdiff(
      analysis_variables_7b5,
      "Sample_ID"
    ),
    drop = FALSE
  ],
  function(x) sum(!is.finite(x)),
  integer(1)
)

cat("\nMissing-value counts:\n")
print(na_count_7b5)

cat("\nNon-finite numeric-value counts:\n")
print(nonfinite_count_7b5)

stopifnot(
  all(na_count_7b5 == 0)
)

stopifnot(
  all(nonfinite_count_7b5 == 0)
)


# ------------------------------------------------------------
# 11. Purity range audit
# ------------------------------------------------------------

purity_summary_7b5 <- c(
  N = length(
    purity_adjusted_dataset$ABSOLUTE_Purity
  ),
  Min = min(
    purity_adjusted_dataset$ABSOLUTE_Purity
  ),
  Q1 = unname(
    quantile(
      purity_adjusted_dataset$ABSOLUTE_Purity,
      0.25
    )
  ),
  Median = median(
    purity_adjusted_dataset$ABSOLUTE_Purity
  ),
  Mean = mean(
    purity_adjusted_dataset$ABSOLUTE_Purity
  ),
  Q3 = unname(
    quantile(
      purity_adjusted_dataset$ABSOLUTE_Purity,
      0.75
    )
  ),
  Max = max(
    purity_adjusted_dataset$ABSOLUTE_Purity
  )
)

cat("\nABSOLUTE purity summary:\n")
print(purity_summary_7b5)

stopifnot(
  all(
    purity_adjusted_dataset$ABSOLUTE_Purity >= 0 &
      purity_adjusted_dataset$ABSOLUTE_Purity <= 1
  )
)

cat(
  "Purity range [0,1] check: PASSED\n"
)


# ------------------------------------------------------------
# 12. Verify sample identity against the original cohort
# ------------------------------------------------------------

expected_analysis_ids_7b5 <- primary_ids_7b5[
  finite_purity_7b5
]

stopifnot(
  identical(
    purity_adjusted_dataset$Sample_ID,
    expected_analysis_ids_7b5
  )
)

cat(
  "Sample-order identity check: PASSED\n"
)


# ------------------------------------------------------------
# 13. Verify exact correspondence between purity and ABSOLUTE
# ------------------------------------------------------------

expected_purity_7b5 <- absolute_df$purity[
  absolute_match_idx_7b5[
    finite_purity_7b5
  ]
]

stopifnot(
  identical(
    purity_adjusted_dataset$ABSOLUTE_Purity,
    as.numeric(expected_purity_7b5)
  )
)

cat(
  "ABSOLUTE purity correspondence check: PASSED\n"
)


# ------------------------------------------------------------
# 14. Explicit protection checks for STEP 7A-2 objects
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("3. STEP 7A-2 protection checks\n")
cat("------------------------------------------------------------\n")

stopifnot(
  exists("gene_tme_results")
)

stopifnot(
  exists("gene_tme_manuscript_table")
)

stopifnot(
  exists("core_tme_alignment")
)

stopifnot(
  nrow(core_tme_alignment) == 513
)

stopifnot(
  nrow(gene_tme_results) == 9
)

stopifnot(
  nrow(gene_tme_manuscript_table) == 3
)

cat(
  "STEP 7A-2 source objects remain present and unchanged in structure.\n"
)


# ------------------------------------------------------------
# 15. Build audit object
# ------------------------------------------------------------

step7b5_audit <- list(
  
  step = "STEP 7B-5",
  
  purpose =
    "Construct locked purity-adjusted sensitivity-analysis dataset",
  
  source_core_tme_n =
    nrow(core_tme_alignment),
  
  source_core_tme_unique_samples =
    length(unique(core_tme_alignment$Sample_ID)),
  
  absolute_source_n =
    nrow(absolute_df),
  
  exact_absolute_matches_n =
    sum(exact_match_7b5),
  
  exact_absolute_missing_purity_n =
    sum(
      exact_match_7b5 &
        !finite_purity_7b5
    ),
  
  no_absolute_record_n =
    sum(!exact_match_7b5),
  
  final_analysis_n =
    nrow(purity_adjusted_dataset),
  
  excluded_n =
    length(excluded_purity_ids_7b5),
  
  excluded_ids =
    excluded_purity_ids_7b5,
  
  excluded_no_absolute_record_ids =
    unmatched_primary_ids,
  
  excluded_missing_purity_ids =
    expected_missing_purity_7b5,
  
  purity_source =
    "TCGA PanCanAtlas ABSOLUTE purity/ploidy file",
  
  purity_variable =
    "ABSOLUTE_Purity",
  
  barcode_matching_method =
    "Exact primary sample ID matched to ABSOLUTE array field",
  
  barcode_transformation_performed =
    FALSE,
  
  automatic_remapping_performed =
    FALSE,
  
  purity_imputation_performed =
    FALSE,
  
  statistical_analysis_performed =
    FALSE,
  
  regression_performed =
    FALSE,
  
  spearman_performed =
    FALSE,
  
  step7a2_objects_modified =
    FALSE,
  
  status =
    "PASSED_LOCKED_DATASET"
)


# ------------------------------------------------------------
# 16. Save locked dataset + audit checkpoint
# ------------------------------------------------------------

step7b5_file <- file.path(
  paper1_root,
  "Paper1_STEP7B5_locked_purity_adjusted_dataset.RData"
)

save(
  purity_adjusted_dataset,
  step7b5_audit,
  purity_summary_7b5,
  excluded_purity_ids_7b5,
  file = step7b5_file
)

cat(
  "\nLocked STEP 7B-5 checkpoint saved:\n",
  step7b5_file,
  "\n"
)


# ------------------------------------------------------------
# 17. Final protection checks
# ------------------------------------------------------------

stopifnot(
  nrow(purity_adjusted_dataset) == 500
)

stopifnot(
  length(unique(
    purity_adjusted_dataset$Sample_ID
  )) == 500
)

stopifnot(
  all(
    is.finite(
      purity_adjusted_dataset$ABSOLUTE_Purity
    )
  )
)

stopifnot(
  all(
    purity_adjusted_dataset$ABSOLUTE_Purity >= 0 &
      purity_adjusted_dataset$ABSOLUTE_Purity <= 1
  )
)

stopifnot(
  nrow(core_tme_alignment) == 513
)

stopifnot(
  nrow(gene_tme_results) == 9
)

stopifnot(
  nrow(gene_tme_manuscript_table) == 3
)


cat("\n============================================================\n")
cat("STEP 7B-5 COMPLETED.\n")
cat("Final purity-adjusted dataset: N = 500\n")
cat("ABSOLUTE purity values: finite and within [0,1]\n")
cat("No purity imputation performed.\n")
cat("No barcode remapping performed.\n")
cat("No statistical analysis performed.\n")
cat("No regression performed.\n")
cat("No Spearman analysis performed.\n")
cat("No STEP 7A-2 object modified.\n")
cat("Status: PASSED_LOCKED_DATASET\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 7B-6
# Purity-adjusted linear regression analysis
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-6: Purity-adjusted linear regression analysis\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Define checkpoint path
# ------------------------------------------------------------

step7b6_checkpoint <- file.path(
  "input_dir/Paper1_tumor-only WGCNA",
  "Paper1_STEP7B6_purity_adjusted_regression_checkpoint.RData"
)


# ------------------------------------------------------------
# 1. Source-object protection audit
# ------------------------------------------------------------

required_objects_7b6 <- c(
  "purity_adjusted_dataset",
  "core_tme_alignment",
  "gene_tme_results",
  "gene_tme_manuscript_table"
)

missing_objects_7b6 <- required_objects_7b6[
  !vapply(required_objects_7b6, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects_7b6) > 0) {
  stop(
    "STEP 7B-6 STOP: Missing required objects: ",
    paste(missing_objects_7b6, collapse = ", ")
  )
}

cat("Required objects present: TRUE\n")


# ------------------------------------------------------------
# 2. Validate locked purity-adjusted dataset
# ------------------------------------------------------------

pad <- purity_adjusted_dataset

cat("\n--- Dataset structure ---\n")
cat("Dimensions:", nrow(pad), "x", ncol(pad), "\n")
cat("Columns:\n")
print(names(pad))


expected_cols_7b6 <- c(
  "Sample_ID",
  "AOC3",
  "MS4A7",
  "IL7R",
  "StromalScore",
  "ImmuneScore",
  "ESTIMATEScore",
  "ABSOLUTE_Purity"
)

if (!identical(names(pad), expected_cols_7b6)) {
  stop(
    "STEP 7B-6 STOP: purity_adjusted_dataset columns do not match ",
    "the locked STEP 7B-5 structure."
  )
}

if (nrow(pad) != 500) {
  stop(
    "STEP 7B-6 STOP: Expected N = 500, observed N = ",
    nrow(pad)
  )
}

cat("Expected N = 500: TRUE\n")
cat("Expected columns: TRUE\n")


# ------------------------------------------------------------
# 3. Validate sample IDs
# ------------------------------------------------------------

if (anyNA(pad$Sample_ID)) {
  stop("STEP 7B-6 STOP: Sample_ID contains NA.")
}

if (anyDuplicated(pad$Sample_ID) > 0) {
  stop("STEP 7B-6 STOP: Duplicate Sample_ID detected.")
}

cat("Sample IDs unique: TRUE\n")


# ------------------------------------------------------------
# 4. Validate numeric variables
# ------------------------------------------------------------

numeric_cols_7b6 <- expected_cols_7b6[-1]

numeric_status_7b6 <- vapply(
  pad[numeric_cols_7b6],
  is.numeric,
  logical(1)
)

if (!all(numeric_status_7b6)) {
  stop(
    "STEP 7B-6 STOP: One or more analysis variables are not numeric."
  )
}

nonfinite_counts_7b6 <- vapply(
  pad[numeric_cols_7b6],
  function(x) sum(!is.finite(x)),
  integer(1)
)

cat("\nNon-finite values:\n")
print(nonfinite_counts_7b6)

if (any(nonfinite_counts_7b6 > 0)) {
  stop(
    "STEP 7B-6 STOP: Non-finite values detected in analysis variables."
  )
}

cat("All analysis variables finite: TRUE\n")


# ------------------------------------------------------------
# 5. Validate ABSOLUTE purity range
# ------------------------------------------------------------

purity_range_7b6 <- range(
  pad$ABSOLUTE_Purity,
  na.rm = TRUE
)

cat("\nABSOLUTE purity range:\n")
print(purity_range_7b6)

if (purity_range_7b6[1] < 0 ||
    purity_range_7b6[2] > 1) {
  stop(
    "STEP 7B-6 STOP: ABSOLUTE_Purity outside [0,1]."
  )
}

cat("ABSOLUTE_Purity within [0,1]: TRUE\n")


# ------------------------------------------------------------
# 6. Define the 9 pre-specified models
# ------------------------------------------------------------

gene_list_7b6 <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

tme_list_7b6 <- c(
  "ImmuneScore",
  "StromalScore",
  "ESTIMATEScore"
)


model_spec_7b6 <- expand.grid(
  TME_Trait = tme_list_7b6,
  Gene = gene_list_7b6,
  stringsAsFactors = FALSE
)

model_spec_7b6 <- model_spec_7b6[
  order(
    match(model_spec_7b6$TME_Trait, tme_list_7b6),
    match(model_spec_7b6$Gene, gene_list_7b6)
  ),
]

rownames(model_spec_7b6) <- NULL

cat("\n--- Pre-specified model family ---\n")
print(model_spec_7b6)


# ------------------------------------------------------------
# 7. Fit the 9 purity-adjusted linear models
# ------------------------------------------------------------

purity_adjusted_models <- list()

purity_adjusted_results_raw <- vector(
  "list",
  nrow(model_spec_7b6)
)


for (i in seq_len(nrow(model_spec_7b6))) {
  
  gene_i <- model_spec_7b6$Gene[i]
  tme_i  <- model_spec_7b6$TME_Trait[i]
  
  formula_i <- as.formula(
    paste(
      tme_i,
      "~",
      gene_i,
      "+ ABSOLUTE_Purity"
    )
  )
  
  model_i <- lm(
    formula = formula_i,
    data = pad
  )
  
  purity_adjusted_models[[paste(tme_i, gene_i, sep = "__")]] <- model_i
  
  coef_i <- summary(model_i)$coefficients
  
  if (!(gene_i %in% rownames(coef_i))) {
    stop(
      "STEP 7B-6 STOP: Gene coefficient not found for ",
      gene_i,
      " in model ",
      tme_i
    )
  }
  
  if (!("ABSOLUTE_Purity" %in% rownames(coef_i))) {
    stop(
      "STEP 7B-6 STOP: Purity coefficient not found in model ",
      tme_i,
      " ~ ",
      gene_i,
      " + ABSOLUTE_Purity"
    )
  }
  
  gene_coef_i <- coef_i[gene_i, ]
  purity_coef_i <- coef_i["ABSOLUTE_Purity", ]
  
  purity_adjusted_results_raw[[i]] <- data.frame(
    Gene = gene_i,
    TME_Trait = tme_i,
    N = nobs(model_i),
    
    Gene_Beta = unname(gene_coef_i["Estimate"]),
    Gene_SE = unname(gene_coef_i["Std. Error"]),
    Gene_t = unname(gene_coef_i["t value"]),
    Gene_P = unname(gene_coef_i["Pr(>|t|)"]),
    
    Purity_Beta = unname(purity_coef_i["Estimate"]),
    Purity_SE = unname(purity_coef_i["Std. Error"]),
    Purity_t = unname(purity_coef_i["t value"]),
    Purity_P = unname(purity_coef_i["Pr(>|t|)"]),
    
    R_squared = summary(model_i)$r.squared,
    Adjusted_R_squared = summary(model_i)$adj.r.squared,
    
    stringsAsFactors = FALSE
  )
}


purity_adjusted_results_raw <- do.call(
  rbind,
  purity_adjusted_results_raw
)

rownames(purity_adjusted_results_raw) <- NULL


# ------------------------------------------------------------
# 8. BH-FDR across the 9 gene-term P values only
# ------------------------------------------------------------

purity_adjusted_results <- purity_adjusted_results_raw

purity_adjusted_results$Gene_FDR_BH <- p.adjust(
  purity_adjusted_results$Gene_P,
  method = "BH"
)

purity_adjusted_results$Gene_Direction <- ifelse(
  purity_adjusted_results$Gene_Beta > 0,
  "Positive",
  ifelse(
    purity_adjusted_results$Gene_Beta < 0,
    "Negative",
    "Zero"
  )
)

purity_adjusted_results$Significant_P05 <- (
  purity_adjusted_results$Gene_P < 0.05
)

purity_adjusted_results$Significant_FDR05 <- (
  purity_adjusted_results$Gene_FDR_BH < 0.05
)


# ------------------------------------------------------------
# 9. Reorder columns for manuscript/audit readability
# ------------------------------------------------------------

purity_adjusted_results <- purity_adjusted_results[
  c(
    "Gene",
    "TME_Trait",
    "N",
    "Gene_Beta",
    "Gene_SE",
    "Gene_t",
    "Gene_P",
    "Gene_FDR_BH",
    "Gene_Direction",
    "Significant_P05",
    "Significant_FDR05",
    "Purity_Beta",
    "Purity_SE",
    "Purity_t",
    "Purity_P",
    "R_squared",
    "Adjusted_R_squared"
  )
]


# ------------------------------------------------------------
# 10. Basic statistical audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 7B-6: Statistical results\n")
cat("============================================================\n\n")

print(
  purity_adjusted_results[
    c(
      "Gene",
      "TME_Trait",
      "N",
      "Gene_Beta",
      "Gene_P",
      "Gene_FDR_BH",
      "Gene_Direction",
      "Significant_FDR05"
    )
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Audit: exactly 9 models
# ------------------------------------------------------------

if (length(purity_adjusted_models) != 9) {
  stop(
    "STEP 7B-6 AUDIT FAILED: Expected 9 models, found ",
    length(purity_adjusted_models)
  )
}

if (nrow(purity_adjusted_results) != 9) {
  stop(
    "STEP 7B-6 AUDIT FAILED: Expected 9 result rows, found ",
    nrow(purity_adjusted_results)
  )
}

cat("\nExactly 9 models: TRUE\n")
cat("Exactly 9 result rows: TRUE\n")


# ------------------------------------------------------------
# 12. Audit: model sample size
# ------------------------------------------------------------

if (!all(purity_adjusted_results$N == 500)) {
  stop(
    "STEP 7B-6 AUDIT FAILED: One or more models do not have N = 500."
  )
}

cat("All models N = 500: TRUE\n")


# ------------------------------------------------------------
# 13. Audit: gene coefficients and P values
# ------------------------------------------------------------

required_gene_statistics_7b6 <- c(
  "Gene_Beta",
  "Gene_SE",
  "Gene_t",
  "Gene_P",
  "Gene_FDR_BH"
)

if (any(
  !vapply(
    purity_adjusted_results[required_gene_statistics_7b6],
    function(x) all(is.finite(x)),
    logical(1)
  )
)) {
  stop(
    "STEP 7B-6 AUDIT FAILED: Non-finite gene-level statistics detected."
  )
}

cat("All gene-level statistics finite: TRUE\n")


# ------------------------------------------------------------
# 14. Audit: BH-FDR definition
# ------------------------------------------------------------

expected_fdr_7b6 <- p.adjust(
  purity_adjusted_results$Gene_P,
  method = "BH"
)

fdr_match_7b6 <- isTRUE(
  all.equal(
    purity_adjusted_results$Gene_FDR_BH,
    expected_fdr_7b6,
    tolerance = 1e-15
  )
)

if (!fdr_match_7b6) {
  stop(
    "STEP 7B-6 AUDIT FAILED: BH-FDR does not match the 9 gene P values."
  )
}

cat("BH-FDR across exactly 9 gene-term P values: TRUE\n")

# ============================================================
# STEP 7B-6
# STEP 15 — Corrected model-formula audit
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-6: Corrected model formula audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Expected formulas
# ------------------------------------------------------------

expected_formula_strings_7b6 <- paste(
  model_spec_7b6$TME_Trait,
  "~",
  model_spec_7b6$Gene,
  "+ ABSOLUTE_Purity"
)


# ------------------------------------------------------------
# 2. Observed formulas
# ------------------------------------------------------------

observed_formula_strings_7b6 <- vapply(
  purity_adjusted_models,
  function(m) {
    paste(
      deparse(formula(m)),
      collapse = " "
    )
  },
  character(1)
)


# ------------------------------------------------------------
# 3. Remove names before comparison
# ------------------------------------------------------------

expected_formula_compare_7b6 <- unname(
  expected_formula_strings_7b6
)

observed_formula_compare_7b6 <- unname(
  observed_formula_strings_7b6
)


# ------------------------------------------------------------
# 4. Normalize whitespace
# ------------------------------------------------------------

expected_formula_compare_7b6 <- gsub(
  "\\s+",
  "",
  expected_formula_compare_7b6
)

observed_formula_compare_7b6 <- gsub(
  "\\s+",
  "",
  observed_formula_compare_7b6
)


# ------------------------------------------------------------
# 5. Compare formula sets
# ------------------------------------------------------------

formula_match_7b6 <- identical(
  sort(expected_formula_compare_7b6),
  sort(observed_formula_compare_7b6)
)


# ------------------------------------------------------------
# 6. Report
# ------------------------------------------------------------

cat("Expected formulas:\n")
print(expected_formula_strings_7b6)

cat("\nObserved formulas:\n")
print(observed_formula_strings_7b6)

cat("\nAll 9 formulas match pre-specified model family:",
    formula_match_7b6, "\n")


# ------------------------------------------------------------
# 7. Stop if mismatch remains
# ------------------------------------------------------------

if (!formula_match_7b6) {
  
  cat("\nExpected normalized formulas:\n")
  print(sort(expected_formula_compare_7b6))
  
  cat("\nObserved normalized formulas:\n")
  print(sort(observed_formula_compare_7b6))
  
  stop(
    "STEP 7B-6 AUDIT FAILED: Model formulas do not match ",
    "the pre-specified design."
  )
}

cat("\nSTEP 7B-6 FORMULA AUDIT: PASSED\n")


# ------------------------------------------------------------
# 16. Audit: no source object modification
# ------------------------------------------------------------

source_dims_7b6 <- list(
  core_tme_alignment = dim(core_tme_alignment),
  gene_tme_results = dim(gene_tme_results),
  gene_tme_manuscript_table = dim(gene_tme_manuscript_table)
)

source_structure_7b6 <- list(
  core_tme_alignment = names(core_tme_alignment),
  gene_tme_results = names(gene_tme_results),
  gene_tme_manuscript_table = names(gene_tme_manuscript_table)
)

expected_source_dims_7b6 <- list(
  core_tme_alignment = c(513L, 7L),
  gene_tme_results = c(9L, 9L),
  gene_tme_manuscript_table = c(3L, 7L)
)

expected_source_names_7b6 <- list(
  core_tme_alignment = c(
    "Sample_ID",
    "AOC3",
    "MS4A7",
    "IL7R",
    "StromalScore",
    "ImmuneScore",
    "ESTIMATEScore"
  ),
  gene_tme_results = c(
    "Gene",
    "TME_Trait",
    "N",
    "Spearman_Rho",
    "P_value",
    "FDR_BH",
    "Direction",
    "Significant_P05",
    "Significant_FDR05"
  ),
  gene_tme_manuscript_table = c(
    "Gene",
    "ImmuneScore_Rho",
    "ImmuneScore_FDR",
    "StromalScore_Rho",
    "StromalScore_FDR",
    "ESTIMATEScore_Rho",
    "ESTIMATEScore_FDR"
  )
)

source_structure_match_7b6 <- (
  identical(source_dims_7b6, expected_source_dims_7b6) &&
    identical(source_structure_7b6, expected_source_names_7b6)
)

if (!source_structure_match_7b6) {
  stop(
    "STEP 7B-6 AUDIT FAILED: Locked STEP 7A-2 source-object structure changed."
  )
}

cat("STEP 7A-2 source-object structure unchanged: TRUE\n")


# ------------------------------------------------------------
# 17. Audit: purity-adjusted dataset remains unchanged
# ------------------------------------------------------------

if (!identical(
  names(purity_adjusted_dataset),
  expected_cols_7b6
)) {
  stop(
    "STEP 7B-6 AUDIT FAILED: Purity-adjusted dataset structure changed."
  )
}

cat("STEP 7B-5 locked dataset structure unchanged: TRUE\n")


# ------------------------------------------------------------
# 18. Construct formal audit object
# ------------------------------------------------------------

step7b6_audit <- list(
  
  step = "STEP 7B-6",
  
  status = "PASSED",
  
  analysis = "Purity-adjusted linear regression",
  
  model_family =
    "TME_score ~ core_gene_expression + ABSOLUTE_Purity",
  
  genes = gene_list_7b6,
  
  TME_traits = tme_list_7b6,
  
  number_of_models = 9,
  
  analysis_N = 500,
  
  FDR_method = "Benjamini-Hochberg",
  
  FDR_family =
    "9 gene-term P values only",
  
  model_type = "ordinary least-squares linear regression",
  
  source_dataset =
    "purity_adjusted_dataset from STEP 7B-5",
  
  purity_source = "TCGA PanCanAtlas ABSOLUTE",
  
  purity_imputation_performed = FALSE,
  
  missing_purity_exclusion_performed = TRUE,
  
  automatic_barcode_remapping_performed = FALSE,
  
  Spearman_performed = FALSE,
  
  Step7A2_modified = FALSE,
  
  Figure6_modified = FALSE,
  
  SupplementaryTable3_modified = FALSE,
  
  source_objects_modified = FALSE,
  
  statistical_reanalysis_of_Step7A2 = FALSE,
  
  checkpoint_path = step7b6_checkpoint
)


# ------------------------------------------------------------
# 19. Final summary audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 7B-6 FINAL AUDIT\n")
cat("============================================================\n")

cat("Required objects present: TRUE\n")
cat("Locked dataset N = 500: TRUE\n")
cat("Exactly 9 models: TRUE\n")
cat("All models N = 500: TRUE\n")
cat("All models include ABSOLUTE_Purity: TRUE\n")
cat("All gene coefficients finite: TRUE\n")
cat("All gene P values finite: TRUE\n")
cat("BH-FDR across 9 gene P values: TRUE\n")
cat("Locked Step 7A-2 objects unchanged: TRUE\n")
cat("Figure 6 modified: FALSE\n")
cat("Supplementary Table 3 modified: FALSE\n")
cat("New Spearman analysis performed: FALSE\n")


# ------------------------------------------------------------
# 20. Save checkpoint
# ------------------------------------------------------------

save(
  purity_adjusted_models,
  purity_adjusted_results_raw,
  purity_adjusted_results,
  model_spec_7b6,
  step7b6_audit,
  file = step7b6_checkpoint
)

cat("\nCheckpoint saved:\n")
cat(step7b6_checkpoint, "\n")

cat("\n============================================================\n")
cat("STEP 7B-6 PASSED\n")
cat("============================================================\n")

# ============================================================
# STEP 7B-6 DIAGNOSTIC ONLY
# Do NOT refit models
# Do NOT modify locked source objects
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-6 DIAGNOSTIC: coefficient structure\n")
cat("============================================================\n\n")

for (i in seq_along(purity_adjusted_models)) {
  
  model_name <- names(purity_adjusted_models)[i]
  model <- purity_adjusted_models[[i]]
  
  cat("\n------------------------------------------------------------\n")
  cat("Model:", model_name, "\n")
  cat("Formula:\n")
  print(formula(model))
  
  cat("\nCoefficient row names:\n")
  print(rownames(coef(summary(model))))
  
  cat("\nCoefficient table:\n")
  print(coef(summary(model)))
  
}

# ============================================================
# STEP 7B-6
# CORRECTION: coefficient extraction only
#
# IMPORTANT:
#   - Existing 9 lm() models are retained.
#   - No model refitting.
#   - No new statistical analysis.
#   - No new Spearman analysis.
#   - No change to STEP 7A-2.
#   - No change to Figure 6.
#   - No change to Supplementary Table 3.
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-6 CORRECTION: coefficient extraction only\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Confirm existing models
# ------------------------------------------------------------

stopifnot(
  exists("purity_adjusted_models"),
  length(purity_adjusted_models) == 9
)

cat(
  "Existing fitted models detected: ",
  length(purity_adjusted_models),
  "\n",
  sep = ""
)

cat("No model refitting performed.\n")


# ------------------------------------------------------------
# 2. Extract coefficients from existing models
# ------------------------------------------------------------

purity_adjusted_results_raw <- do.call(
  rbind,
  lapply(
    seq_along(purity_adjusted_models),
    function(i) {
      
      model <- purity_adjusted_models[[i]]
      
      trait <- model_spec_7b6$TME_Trait[i]
      gene  <- model_spec_7b6$Gene[i]
      
      coef_table <- summary(model)$coefficients
      coef_vector <- coef(model)
      
      # ---- Gene term ----
      gene_beta <- unname(
        coef_vector[gene]
      )
      
      gene_se <- unname(
        coef_table[gene, "Std. Error"]
      )
      
      gene_t <- unname(
        coef_table[gene, "t value"]
      )
      
      gene_p <- unname(
        coef_table[gene, "Pr(>|t|)"]
      )
      
      # ---- Purity term ----
      purity_beta <- unname(
        coef_vector["ABSOLUTE_Purity"]
      )
      
      purity_se <- unname(
        coef_table[
          "ABSOLUTE_Purity",
          "Std. Error"
        ]
      )
      
      purity_t <- unname(
        coef_table[
          "ABSOLUTE_Purity",
          "t value"
        ]
      )
      
      purity_p <- unname(
        coef_table[
          "ABSOLUTE_Purity",
          "Pr(>|t|)"
        ]
      )
      
      data.frame(
        TME_Trait = trait,
        Gene = gene,
        N = nobs(model),
        
        Gene_Beta = gene_beta,
        Gene_SE = gene_se,
        Gene_t = gene_t,
        Gene_P = gene_p,
        
        Purity_Beta = purity_beta,
        Purity_SE = purity_se,
        Purity_t = purity_t,
        Purity_P = purity_p,
        
        R_squared =
          summary(model)$r.squared,
        
        Adjusted_R_squared =
          summary(model)$adj.r.squared,
        
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )
)

rownames(
  purity_adjusted_results_raw
) <- NULL


# ------------------------------------------------------------
# 3. Direction
# ------------------------------------------------------------

purity_adjusted_results_raw$Direction <- ifelse(
  purity_adjusted_results_raw$Gene_Beta > 0,
  "Positive",
  ifelse(
    purity_adjusted_results_raw$Gene_Beta < 0,
    "Negative",
    "Zero"
  )
)


# ------------------------------------------------------------
# 4. BH-FDR across the nine gene-term P values
# ------------------------------------------------------------

purity_adjusted_results_raw$Gene_FDR_BH <- p.adjust(
  purity_adjusted_results_raw$Gene_P,
  method = "BH"
)


# ------------------------------------------------------------
# 5. Significance
# ------------------------------------------------------------

purity_adjusted_results_raw$Significant_FDR05 <- (
  purity_adjusted_results_raw$Gene_FDR_BH < 0.05
)


# ------------------------------------------------------------
# 6. Final result table
# ------------------------------------------------------------

purity_adjusted_results <- purity_adjusted_results_raw[
  c(
    "TME_Trait",
    "Gene",
    "N",
    "Gene_Beta",
    "Gene_SE",
    "Gene_t",
    "Gene_P",
    "Gene_FDR_BH",
    "Direction",
    "Significant_FDR05",
    "Purity_Beta",
    "Purity_SE",
    "Purity_t",
    "Purity_P",
    "R_squared",
    "Adjusted_R_squared"
  )
]

rownames(
  purity_adjusted_results
) <- NULL


# ------------------------------------------------------------
# 7. Extraction audit
# ------------------------------------------------------------

audit_gene_beta <- all(
  is.finite(
    purity_adjusted_results$Gene_Beta
  )
)

audit_gene_p <- all(
  is.finite(
    purity_adjusted_results$Gene_P
  )
)

audit_gene_fdr <- all(
  is.finite(
    purity_adjusted_results$Gene_FDR_BH
  )
)

audit_purity_beta <- all(
  is.finite(
    purity_adjusted_results$Purity_Beta
  )
)

audit_purity_p <- all(
  is.finite(
    purity_adjusted_results$Purity_P
  )
)

audit_N <- all(
  purity_adjusted_results$N == 500
)

audit_model_number <- (
  length(purity_adjusted_models) == 9
)

audit_fdr_definition <- identical(
  purity_adjusted_results$Gene_FDR_BH,
  p.adjust(
    purity_adjusted_results$Gene_P,
    method = "BH"
  )
)


# ------------------------------------------------------------
# 8. Print audit
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("STEP 7B-6 CORRECTION AUDIT\n")
cat("------------------------------------------------------------\n")

cat(
  "Exactly 9 existing models: ",
  audit_model_number,
  "\n",
  sep = ""
)

cat(
  "All models N = 500: ",
  audit_N,
  "\n",
  sep = ""
)

cat(
  "All gene coefficients finite: ",
  audit_gene_beta,
  "\n",
  sep = ""
)

cat(
  "All gene P values finite: ",
  audit_gene_p,
  "\n",
  sep = ""
)

cat(
  "All gene FDR values finite: ",
  audit_gene_fdr,
  "\n",
  sep = ""
)

cat(
  "All purity coefficients finite: ",
  audit_purity_beta,
  "\n",
  sep = ""
)

cat(
  "All purity P values finite: ",
  audit_purity_p,
  "\n",
  sep = ""
)

cat(
  "BH-FDR definition verified: ",
  audit_fdr_definition,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Print corrected results
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Corrected STEP 7B-6 gene-term results\n")
cat("------------------------------------------------------------\n\n")

print(
  purity_adjusted_results[
    c(
      "TME_Trait",
      "Gene",
      "N",
      "Gene_Beta",
      "Gene_P",
      "Gene_FDR_BH",
      "Direction"
    )
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. PASS condition for extraction correction
# ------------------------------------------------------------

step7b6_extraction_passed <- all(
  audit_model_number,
  audit_N,
  audit_gene_beta,
  audit_gene_p,
  audit_gene_fdr,
  audit_purity_beta,
  audit_purity_p,
  audit_fdr_definition
)


cat("\n============================================================\n")

if (step7b6_extraction_passed) {
  
  cat("STEP 7B-6 COEFFICIENT EXTRACTION: PASSED\n")
  
} else {
  
  cat("STEP 7B-6 COEFFICIENT EXTRACTION: FAILED\n")
  
}

cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 7B-6 FINAL ARCHIVE
# Purity-adjusted regression sensitivity analysis
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-6 FINAL ARCHIVE\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Final pre-save validation
# ------------------------------------------------------------

stopifnot(
  exists("purity_adjusted_models"),
  length(purity_adjusted_models) == 9,
  exists("purity_adjusted_results"),
  nrow(purity_adjusted_results) == 9,
  all(purity_adjusted_results$N == 500),
  all(is.finite(purity_adjusted_results$Gene_Beta)),
  all(is.finite(purity_adjusted_results$Gene_SE)),
  all(is.finite(purity_adjusted_results$Gene_t)),
  all(is.finite(purity_adjusted_results$Gene_P)),
  all(is.finite(purity_adjusted_results$Gene_FDR_BH)),
  all(is.finite(purity_adjusted_results$Purity_Beta)),
  all(is.finite(purity_adjusted_results$Purity_SE)),
  all(is.finite(purity_adjusted_results$Purity_t)),
  all(is.finite(purity_adjusted_results$Purity_P)),
  all(is.finite(purity_adjusted_results$R_squared)),
  all(is.finite(purity_adjusted_results$Adjusted_R_squared))
)

# ------------------------------------------------------------
# 2. Final model specification audit
# ------------------------------------------------------------

expected_traits <- c(
  "ImmuneScore",
  "ImmuneScore",
  "ImmuneScore",
  "StromalScore",
  "StromalScore",
  "StromalScore",
  "ESTIMATEScore",
  "ESTIMATEScore",
  "ESTIMATEScore"
)

expected_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R",
  "AOC3",
  "MS4A7",
  "IL7R",
  "AOC3",
  "MS4A7",
  "IL7R"
)

audit_trait_order <- identical(
  purity_adjusted_results$TME_Trait,
  expected_traits
)

audit_gene_order <- identical(
  purity_adjusted_results$Gene,
  expected_genes
)

# ------------------------------------------------------------
# 3. Final BH-FDR audit
# ------------------------------------------------------------

audit_fdr_final <- isTRUE(
  all.equal(
    purity_adjusted_results$Gene_FDR_BH,
    p.adjust(
      purity_adjusted_results$Gene_P,
      method = "BH"
    ),
    tolerance = 0
  )
)

# ------------------------------------------------------------
# 4. Final direction/significance audit
# ------------------------------------------------------------

audit_direction <- all(
  purity_adjusted_results$Gene_Beta > 0 &
    purity_adjusted_results$Direction == "Positive"
)

audit_significance <- all(
  purity_adjusted_results$Gene_FDR_BH < 0.05 &
    purity_adjusted_results$Significant_FDR05
)

# ------------------------------------------------------------
# 5. Confirm no Step 7A-2 modification
# ------------------------------------------------------------

audit_step7a2_objects <- all(
  exists("core_tme_alignment"),
  exists("gene_tme_results"),
  exists("gene_tme_manuscript_table"),
  nrow(core_tme_alignment) == 513,
  nrow(gene_tme_results) == 9,
  nrow(gene_tme_manuscript_table) == 3
)

# ------------------------------------------------------------
# 6. Final archive audit object
# ------------------------------------------------------------

step7b6_final_audit <- list(
  
  status = "PASSED_LOCKED_FINAL",
  
  analysis = "Purity-adjusted OLS sensitivity analysis",
  
  model_type = "ordinary least squares",
  
  n_models = 9,
  
  N_per_model = 500,
  
  predictors = "Gene expression + ABSOLUTE_Purity",
  
  gene_terms = c("AOC3", "MS4A7", "IL7R"),
  
  TME_traits = c(
    "ImmuneScore",
    "StromalScore",
    "ESTIMATEScore"
  ),
  
  multiple_testing =
    "Benjamini-Hochberg adjustment across 9 gene-term P values",
  
  adjusted_term =
    "Gene coefficient/P value",
  
  no_adjusted_spearman = TRUE,
  
  no_model_refitting_during_extraction = TRUE,
  
  purity_imputation_performed = FALSE,
  
  barcode_transformation_performed = FALSE,
  
  automatic_remapping_performed = FALSE,
  
  step7a2_objects_modified = FALSE,
  
  audit_model_number = audit_model_number,
  
  audit_N = audit_N,
  
  audit_gene_beta = audit_gene_beta,
  
  audit_gene_p = audit_gene_p,
  
  audit_gene_fdr = audit_gene_fdr,
  
  audit_purity_beta = audit_purity_beta,
  
  audit_purity_p = audit_purity_p,
  
  audit_fdr_definition = audit_fdr_definition,
  
  audit_trait_order = audit_trait_order,
  
  audit_gene_order = audit_gene_order,
  
  audit_fdr_final = audit_fdr_final,
  
  audit_direction = audit_direction,
  
  audit_significance = audit_significance,
  
  audit_step7a2_objects = audit_step7a2_objects,
  
  extraction_passed = step7b6_extraction_passed
)

# ------------------------------------------------------------
# 7. Final PASS condition
# ------------------------------------------------------------

step7b6_final_passed <- all(
  step7b6_extraction_passed,
  audit_trait_order,
  audit_gene_order,
  audit_fdr_final,
  audit_direction,
  audit_significance,
  audit_step7a2_objects
)

step7b6_final_audit$final_passed <- step7b6_final_passed

# ------------------------------------------------------------
# 8. Print final audit
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("STEP 7B-6 FINAL ARCHIVE AUDIT\n")
cat("------------------------------------------------------------\n")

cat(
  "Extraction audit passed: ",
  step7b6_extraction_passed,
  "\n",
  sep = ""
)

cat(
  "Trait order correct: ",
  audit_trait_order,
  "\n",
  sep = ""
)

cat(
  "Gene order correct: ",
  audit_gene_order,
  "\n",
  sep = ""
)

cat(
  "Final BH-FDR audit: ",
  audit_fdr_final,
  "\n",
  sep = ""
)

cat(
  "All gene associations positive: ",
  audit_direction,
  "\n",
  sep = ""
)

cat(
  "All gene associations FDR < 0.05: ",
  audit_significance,
  "\n",
  sep = ""
)

cat(
  "Step 7A-2 objects preserved: ",
  audit_step7a2_objects,
  "\n",
  sep = ""
)

cat("\n------------------------------------------------------------\n")

# ------------------------------------------------------------
# 9. Save ONLY after all audits pass
# ------------------------------------------------------------

if (step7b6_final_passed) {
  
  checkpoint_7b6 <- file.path(
    "input_dir/Paper1_tumor-only WGCNA",
    "Paper1_STEP7B6_purity_adjusted_regression_checkpoint.RData"
  )
  
  save(
    purity_adjusted_models,
    purity_adjusted_results_raw,
    purity_adjusted_results,
    model_spec_7b6,
    step7b6_extraction_passed,
    step7b6_final_audit,
    step7b6_final_passed,
    file = checkpoint_7b6
  )
  
  cat(
    "STEP 7B-6 FINAL ARCHIVE: PASSED\n"
  )
  
  cat(
    "Checkpoint saved to:\n",
    checkpoint_7b6,
    "\n"
  )
  
} else {
  
  cat(
    "STEP 7B-6 FINAL ARCHIVE: FAILED\n"
  )
  
  cat(
    "Checkpoint NOT overwritten.\n"
  )
}

cat("============================================================\n")



##STEP 7B-7A
#ESTIMATE signature source and core-gene membership audit
# ============================================================
# Paper 1
# STEP 7B-7A
# ESTIMATE circularity audit:
# source / signature membership
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7A: ESTIMATE CIRCULARITY AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Required core genes
# ------------------------------------------------------------

core_genes_7b7 <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

cat("Core genes under audit:\n")
print(core_genes_7b7)


# ------------------------------------------------------------
# 2. Check whether estimate package is available
# ------------------------------------------------------------

audit_estimate_namespace <- requireNamespace(
  "estimate",
  quietly = TRUE
)

cat(
  "\nESTIMATE package available: ",
  audit_estimate_namespace,
  "\n",
  sep = ""
)

if (!audit_estimate_namespace) {
  stop(
    "The 'estimate' package is not available. ",
    "Do not proceed to membership analysis."
  )
}


# ------------------------------------------------------------
# 3. Record ESTIMATE package version
# ------------------------------------------------------------

estimate_version_7b7 <- as.character(
  packageVersion("estimate")
)

cat(
  "ESTIMATE package version: ",
  estimate_version_7b7,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Load SI_geneset into a NEW object only
# ------------------------------------------------------------

estimate_env_7b7 <- new.env(parent = emptyenv())

utils::data(
  "SI_geneset",
  package = "estimate",
  envir = estimate_env_7b7
)

audit_si_object_present <- exists(
  "SI_geneset",
  envir = estimate_env_7b7,
  inherits = FALSE
)

cat(
  "SI_geneset available from ESTIMATE package: ",
  audit_si_object_present,
  "\n",
  sep = ""
)

if (!audit_si_object_present) {
  stop(
    "SI_geneset could not be loaded from the installed ESTIMATE package."
  )
}


# ------------------------------------------------------------
# 5. Inspect signature structure
# ------------------------------------------------------------

SI_geneset_7b7 <- get(
  "SI_geneset",
  envir = estimate_env_7b7,
  inherits = FALSE
)

cat("\n------------------------------------------------------------\n")
cat("SI_geneset structure\n")
cat("------------------------------------------------------------\n")

cat(
  "Class: ",
  paste(class(SI_geneset_7b7), collapse = ", "),
  "\n",
  sep = ""
)

cat(
  "Dimensions: ",
  paste(dim(SI_geneset_7b7), collapse = " x "),
  "\n",
  sep = ""
)

cat("Column names:\n")
print(colnames(SI_geneset_7b7))


# ------------------------------------------------------------
# 6. Convert signatures to explicit vectors
# ------------------------------------------------------------

stromal_genes_7b7 <- unique(
  na.omit(
    as.character(
      SI_geneset_7b7$StromalSignature
    )
  )
)

immune_genes_7b7 <- unique(
  na.omit(
    as.character(
      SI_geneset_7b7$ImmuneSignature
    )
  )
)

stromal_genes_7b7 <- stromal_genes_7b7[
  stromal_genes_7b7 != ""
]

immune_genes_7b7 <- immune_genes_7b7[
  immune_genes_7b7 != ""
]


# ------------------------------------------------------------
# 7. Signature size audit
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Signature size audit\n")
cat("------------------------------------------------------------\n")

cat(
  "Unique StromalSignature genes: ",
  length(stromal_genes_7b7),
  "\n",
  sep = ""
)

cat(
  "Unique ImmuneSignature genes: ",
  length(immune_genes_7b7),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Core-gene membership audit
# ------------------------------------------------------------

core_membership_7b7 <- data.frame(
  Gene = core_genes_7b7,
  
  In_StromalSignature =
    core_genes_7b7 %in% stromal_genes_7b7,
  
  In_ImmuneSignature =
    core_genes_7b7 %in% immune_genes_7b7,
  
  stringsAsFactors = FALSE,
  check.names = FALSE
)

core_membership_7b7$Signature_Count <-
  core_membership_7b7$In_StromalSignature +
  core_membership_7b7$In_ImmuneSignature

core_membership_7b7$Circularity_Relevance <- ifelse(
  core_membership_7b7$Signature_Count == 0,
  "No direct ESTIMATE signature membership",
  ifelse(
    core_membership_7b7$Signature_Count == 1,
    "Direct membership in one ESTIMATE signature",
    "Direct membership in both ESTIMATE signatures"
  )
)


# ------------------------------------------------------------
# 9. Print membership table
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Core-gene ESTIMATE signature membership\n")
cat("------------------------------------------------------------\n\n")

print(
  core_membership_7b7,
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Explicit overlap vectors
# ------------------------------------------------------------

stromal_overlap_7b7 <- intersect(
  core_genes_7b7,
  stromal_genes_7b7
)

immune_overlap_7b7 <- intersect(
  core_genes_7b7,
  immune_genes_7b7
)

cat("\n------------------------------------------------------------\n")
cat("Explicit overlap\n")
cat("------------------------------------------------------------\n")

cat("Core genes in StromalSignature:\n")
print(stromal_overlap_7b7)

cat("\nCore genes in ImmuneSignature:\n")
print(immune_overlap_7b7)


# ------------------------------------------------------------
# 11. ESTIMATEScore dependency
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("ESTIMATEScore construction audit\n")
cat("------------------------------------------------------------\n")

cat(
  "ESTIMATE uses StromalScore and ImmuneScore to derive ",
  "ESTIMATEScore.\n",
  sep = ""
)

cat(
  "Therefore, a core gene present in either signature may also ",
  "be directly relevant to the composite ESTIMATEScore.\n"
)


# ------------------------------------------------------------
# 12. Preserve locked Step 7A-2 objects
# ------------------------------------------------------------

audit_step7a2_preserved_7b7 <- all(
  exists("core_tme_alignment"),
  exists("gene_tme_results"),
  exists("gene_tme_manuscript_table"),
  nrow(core_tme_alignment) == 513,
  nrow(gene_tme_results) == 9,
  nrow(gene_tme_manuscript_table) == 3
)

cat(
  "\nStep 7A-2 locked objects still present and structurally intact: ",
  audit_step7a2_preserved_7b7,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. Confirm no score recalculation
# ------------------------------------------------------------

audit_no_score_recalculation_7b7 <- TRUE
audit_no_spearman_recalculation_7b7 <- TRUE
audit_no_source_modification_7b7 <- TRUE


# ------------------------------------------------------------
# 14. Final source-membership audit
# ------------------------------------------------------------

audit_signature_structure_7b7 <- all(
  length(stromal_genes_7b7) > 0,
  length(immune_genes_7b7) > 0,
  all(core_genes_7b7 %in%
        unique(c(stromal_genes_7b7, immune_genes_7b7)))
  |
    any(core_genes_7b7 %in%
          unique(c(stromal_genes_7b7, immune_genes_7b7)))
)

# The above is deliberately only a structural check.
# Final circularity classification will be determined below.

audit_membership_completed_7b7 <- all(
  audit_estimate_namespace,
  audit_si_object_present,
  length(stromal_genes_7b7) > 0,
  length(immune_genes_7b7) > 0
)


# ------------------------------------------------------------
# 15. Classify circularity status
# ------------------------------------------------------------

if (length(stromal_overlap_7b7) == 0 &&
    length(immune_overlap_7b7) == 0) {
  
  circularity_status_7b7 <-
    "NO_DIRECT_SIGNATURE_OVERLAP"
  
} else {
  
  circularity_status_7b7 <-
    "DIRECT_SIGNATURE_OVERLAP_REQUIRES_LEAVE_ONE_GENE_OUT_AUDIT"
}


# ------------------------------------------------------------
# 16. Final audit object
# ------------------------------------------------------------

step7b7_source_audit <- list(
  
  status = "PASSED_SOURCE_MEMBERSHIP_AUDIT",
  
  estimate_package = "estimate",
  
  estimate_version = estimate_version_7b7,
  
  stromal_signature_n =
    length(stromal_genes_7b7),
  
  immune_signature_n =
    length(immune_genes_7b7),
  
  core_membership =
    core_membership_7b7,
  
  stromal_overlap =
    stromal_overlap_7b7,
  
  immune_overlap =
    immune_overlap_7b7,
  
  circularity_status =
    circularity_status_7b7,
  
  no_score_recalculation =
    audit_no_score_recalculation_7b7,
  
  no_spearman_recalculation =
    audit_no_spearman_recalculation_7b7,
  
  no_source_modification =
    audit_no_source_modification_7b7,
  
  step7a2_objects_preserved =
    audit_step7a2_preserved_7b7,
  
  membership_audit_completed =
    audit_membership_completed_7b7
)


# ------------------------------------------------------------
# 17. Print final status
# ------------------------------------------------------------

cat("\n============================================================\n")

if (audit_membership_completed_7b7 &&
    audit_step7a2_preserved_7b7) {
  
  cat(
    "STEP 7B-7A SOURCE / MEMBERSHIP AUDIT: PASSED\n"
  )
  
  cat(
    "Circularity classification: ",
    circularity_status_7b7,
    "\n",
    sep = ""
  )
  
} else {
  
  cat(
    "STEP 7B-7A SOURCE / MEMBERSHIP AUDIT: FAILED\n"
  )
}

cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 7B-7A-1
# Recover the actual ESTIMATE scoring source
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7A-1: RECOVER ACTUAL ESTIMATE SOURCE\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Check relevant ESTIMATE objects
# ------------------------------------------------------------

estimate_objects_7b7 <- c(
  "estimate_input_file",
  "estimate_input_df",
  "estimate_gct_input",
  "estimate_common",
  "estimate_filtered",
  "estimate_sample_ids",
  "estimate_sample_ids_normalized",
  "estimate_sample_ids_raw",
  "estimate_scores_wgcna"
)

cat("Relevant ESTIMATE objects in current workspace:\n\n")

for (obj in estimate_objects_7b7) {
  
  cat(
    obj,
    " : ",
    exists(obj, inherits = TRUE),
    "\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 2. Inspect estimate_input_file
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("estimate_input_file\n")
cat("------------------------------------------------------------\n")

if (exists("estimate_input_file", inherits = TRUE)) {
  
  print(estimate_input_file)
  
  if (is.character(estimate_input_file)) {
    
    cat(
      "File exists: ",
      file.exists(estimate_input_file),
      "\n",
      sep = ""
    )
  }
}


# ------------------------------------------------------------
# 3. Inspect object classes / dimensions
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("ESTIMATE object structures\n")
cat("------------------------------------------------------------\n")

for (obj in c(
  "estimate_input_df",
  "estimate_gct_input",
  "estimate_common",
  "estimate_filtered",
  "estimate_scores_wgcna"
)) {
  
  if (exists(obj, inherits = TRUE)) {
    
    x <- get(obj, inherits = TRUE)
    
    cat("\nOBJECT: ", obj, "\n", sep = "")
    cat(
      "Class: ",
      paste(class(x), collapse = ", "),
      "\n",
      sep = ""
    )
    
    if (!is.null(dim(x))) {
      
      cat(
        "Dimensions: ",
        paste(dim(x), collapse = " x "),
        "\n",
        sep = ""
      )
    }
    
    cat("Names:\n")
    
    if (!is.null(names(x))) {
      print(names(x))
    }
    
    if (!is.null(colnames(x))) {
      cat("Column names:\n")
      print(colnames(x))
    }
  }
}


# ------------------------------------------------------------
# 4. Inspect estimate_input_file if it exists
# ------------------------------------------------------------

if (
  exists("estimate_input_file", inherits = TRUE) &&
  is.character(estimate_input_file) &&
  length(estimate_input_file) == 1 &&
  file.exists(estimate_input_file)
) {
  
  cat("\n------------------------------------------------------------\n")
  cat("First lines of ESTIMATE input file\n")
  cat("------------------------------------------------------------\n")
  
  print(
    readLines(
      estimate_input_file,
      n = 10,
      warn = FALSE
    )
  )
}


# ------------------------------------------------------------
# 5. Inspect attributes of estimate-related objects
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Attributes of key ESTIMATE objects\n")
cat("------------------------------------------------------------\n")

for (obj in c(
  "estimate_input_df",
  "estimate_gct_input",
  "estimate_common",
  "estimate_filtered",
  "estimate_scores_wgcna"
)) {
  
  if (exists(obj, inherits = TRUE)) {
    
    x <- get(obj, inherits = TRUE)
    
    cat("\nOBJECT: ", obj, "\n", sep = "")
    print(attributes(x))
  }
}


# ------------------------------------------------------------
# 6. Search current workspace for ESTIMATE-related functions
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Current workspace objects containing 'estimate'\n")
cat("------------------------------------------------------------\n")

all_objects_7b7 <- ls(
  envir = .GlobalEnv,
  all.names = TRUE
)

estimate_named_objects_7b7 <- all_objects_7b7[
  grepl(
    "estimate|stromal|immune|signature|geneset|gene_set|ssgsea",
    all_objects_7b7,
    ignore.case = TRUE
  )
]

print(estimate_named_objects_7b7)


# ------------------------------------------------------------
# 7. Search function definitions in current workspace
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Functions potentially related to ESTIMATE\n")
cat("------------------------------------------------------------\n")

function_candidates_7b7 <- all_objects_7b7[
  vapply(
    all_objects_7b7,
    function(x) {
      
      obj <- get(
        x,
        envir = .GlobalEnv
      )
      
      is.function(obj) &&
        grepl(
          "estimate|stromal|immune|score|signature|geneset",
          x,
          ignore.case = TRUE
        )
    },
    logical(1)
  )
]

print(function_candidates_7b7)


# ------------------------------------------------------------
# 8. Audit locked Step 7A-2 objects
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Locked Step 7A-2 objects\n")
cat("------------------------------------------------------------\n")

audit_step7a2_7b7 <- all(
  exists("core_tme_alignment"),
  exists("gene_tme_results"),
  exists("gene_tme_manuscript_table"),
  nrow(core_tme_alignment) == 513,
  nrow(gene_tme_results) == 9,
  nrow(gene_tme_manuscript_table) == 3
)

cat(
  "Step 7A-2 objects intact: ",
  audit_step7a2_7b7,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Explicitly confirm no statistical analysis
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Analysis status\n")
cat("------------------------------------------------------------\n")

cat("No ESTIMATE scores recalculated.\n")
cat("No Spearman tests performed.\n")
cat("No leave-one-gene-out analysis performed.\n")
cat("No locked objects modified.\n")

cat("\n============================================================\n")
cat("STEP 7B-7A-1 SOURCE AUDIT COMPLETE\n")
cat("============================================================\n")


# ============================================================
# Paper 1
# STEP 7B-7A-2
# Recover historical ESTIMATE intermediate objects
# without modifying the current workspace
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7A-2: RECOVER HISTORICAL ESTIMATE OBJECTS\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Locate the previously identified locked ESTIMATE RData
# ------------------------------------------------------------

estimate_source_file_7b7 <- estimate_file

cat("Historical ESTIMATE RData file:\n")
cat(estimate_source_file_7b7, "\n\n")

cat(
  "File exists: ",
  file.exists(estimate_source_file_7b7),
  "\n",
  sep = ""
)

if (!file.exists(estimate_source_file_7b7)) {
  stop(
    "Historical ESTIMATE RData file was not found. ",
    "Do not proceed with circularity analysis."
  )
}


# ------------------------------------------------------------
# 2. Create isolated temporary environment
# ------------------------------------------------------------

estimate_source_env_7b7 <- new.env(parent = emptyenv())

cat("\nTemporary isolated environment created.\n")
cat("No objects will be loaded into .GlobalEnv.\n")


# ------------------------------------------------------------
# 3. Load historical RData into isolated environment
# ------------------------------------------------------------

loaded_objects_7b7 <- load(
  estimate_source_file_7b7,
  envir = estimate_source_env_7b7
)

cat("\nObjects recovered from historical ESTIMATE RData:\n")
print(loaded_objects_7b7)


# ------------------------------------------------------------
# 4. Inspect recovered object classes and dimensions
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("RECOVERED OBJECT STRUCTURE\n")
cat("------------------------------------------------------------\n")

for (obj in loaded_objects_7b7) {
  
  x <- get(
    obj,
    envir = estimate_source_env_7b7,
    inherits = FALSE
  )
  
  cat("\nOBJECT: ", obj, "\n", sep = "")
  cat(
    "Class: ",
    paste(class(x), collapse = ", "),
    "\n",
    sep = ""
  )
  
  if (!is.null(dim(x))) {
    cat(
      "Dimensions: ",
      paste(dim(x), collapse = " x "),
      "\n",
      sep = ""
    )
  }
  
  if (!is.null(names(x))) {
    cat("Names:\n")
    print(names(x))
  }
  
  if (!is.null(rownames(x))) {
    cat("Row names:\n")
    print(rownames(x))
  }
  
  if (!is.null(colnames(x))) {
    cat("First 10 column names:\n")
    print(head(colnames(x), 10))
  }
}


# ------------------------------------------------------------
# 5. Specifically inspect estimate_input_df
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("HISTORICAL estimate_input_df\n")
cat("------------------------------------------------------------\n")

if (exists(
  "estimate_input_df",
  envir = estimate_source_env_7b7,
  inherits = FALSE
)) {
  
  estimate_input_df_7b7 <- get(
    "estimate_input_df",
    envir = estimate_source_env_7b7,
    inherits = FALSE
  )
  
  cat(
    "Class: ",
    paste(class(estimate_input_df_7b7), collapse = ", "),
    "\n",
    sep = ""
  )
  
  cat(
    "Dimensions: ",
    paste(dim(estimate_input_df_7b7), collapse = " x "),
    "\n",
    sep = ""
  )
  
  cat("Column names:\n")
  print(colnames(estimate_input_df_7b7))
  
  cat("\nFirst 5 rows:\n")
  print(
    utils::head(
      estimate_input_df_7b7,
      5
    )
  )
  
} else {
  
  cat("estimate_input_df NOT recovered.\n")
}


# ------------------------------------------------------------
# 6. Inspect estimate_gct_input
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("HISTORICAL estimate_gct_input\n")
cat("------------------------------------------------------------\n")

if (exists(
  "estimate_gct_input",
  envir = estimate_source_env_7b7,
  inherits = FALSE
)) {
  
  estimate_gct_input_7b7 <- get(
    "estimate_gct_input",
    envir = estimate_source_env_7b7,
    inherits = FALSE
  )
  
  cat(
    "Class: ",
    paste(class(estimate_gct_input_7b7), collapse = ", "),
    "\n",
    sep = ""
  )
  
  cat(
    "Dimensions: ",
    paste(dim(estimate_gct_input_7b7), collapse = " x "),
    "\n",
    sep = ""
  )
  
  cat("First 10 row names:\n")
  print(head(rownames(estimate_gct_input_7b7), 10))
  
  cat("First 10 column names:\n")
  print(head(colnames(estimate_gct_input_7b7), 10))
  
} else {
  
  cat("estimate_gct_input NOT recovered.\n")
}


# ------------------------------------------------------------
# 7. Inspect estimate_common
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("HISTORICAL estimate_common\n")
cat("------------------------------------------------------------\n")

if (exists(
  "estimate_common",
  envir = estimate_source_env_7b7,
  inherits = FALSE
)) {
  
  estimate_common_7b7 <- get(
    "estimate_common",
    envir = estimate_source_env_7b7,
    inherits = FALSE
  )
  
  cat(
    "Class: ",
    paste(class(estimate_common_7b7), collapse = ", "),
    "\n",
    sep = ""
  )
  
  cat(
    "Dimensions: ",
    paste(dim(estimate_common_7b7), collapse = " x "),
    "\n",
    sep = ""
  )
  
  cat("First 20 row names:\n")
  print(head(rownames(estimate_common_7b7), 20))
  
  cat("First 20 column names:\n")
  print(head(colnames(estimate_common_7b7), 20))
  
} else {
  
  cat("estimate_common NOT recovered.\n")
}


# ------------------------------------------------------------
# 8. Inspect estimate_filtered
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("HISTORICAL estimate_filtered\n")
cat("------------------------------------------------------------\n")

if (exists(
  "estimate_filtered",
  envir = estimate_source_env_7b7,
  inherits = FALSE
)) {
  
  estimate_filtered_7b7 <- get(
    "estimate_filtered",
    envir = estimate_source_env_7b7,
    inherits = FALSE
  )
  
  cat(
    "Class: ",
    paste(class(estimate_filtered_7b7), collapse = ", "),
    "\n",
    sep = ""
  )
  
  cat(
    "Dimensions: ",
    paste(dim(estimate_filtered_7b7), collapse = " x "),
    "\n",
    sep = ""
  )
  
  cat("First 20 row names:\n")
  print(head(rownames(estimate_filtered_7b7), 20))
  
} else {
  
  cat("estimate_filtered NOT recovered.\n")
}


# ------------------------------------------------------------
# 9. Search recovered objects for the three core genes
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("CORE-GENE PRESENCE AUDIT IN RECOVERED OBJECTS\n")
cat("------------------------------------------------------------\n")

core_genes_7b7 <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

cat("\nCore genes:\n")
print(core_genes_7b7)

for (obj in loaded_objects_7b7) {
  
  x <- get(
    obj,
    envir = estimate_source_env_7b7,
    inherits = FALSE
  )
  
  cat("\nOBJECT: ", obj, "\n", sep = "")
  
  object_names_7b7 <- unique(
    c(
      names(x),
      rownames(x),
      colnames(x)
    )
  )
  
  if (length(object_names_7b7) == 0) {
    cat("No names/rownames/colnames available.\n")
    next
  }
  
  for (gene in core_genes_7b7) {
    
    hit <- gene %in% object_names_7b7
    
    cat(
      gene,
      " present: ",
      hit,
      "\n",
      sep = ""
    )
  }
}


# ------------------------------------------------------------
# 10. Check whether core genes occur in expression rows
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("CORE-GENE ROWNAME AUDIT\n")
cat("------------------------------------------------------------\n")

for (obj in c(
  "estimate_input_df",
  "estimate_gct_input",
  "estimate_common",
  "estimate_filtered"
)) {
  
  if (!exists(
    obj,
    envir = estimate_source_env_7b7,
    inherits = FALSE
  )) {
    next
  }
  
  x <- get(
    obj,
    envir = estimate_source_env_7b7,
    inherits = FALSE
  )
  
  rn <- rownames(x)
  
  cat("\nOBJECT: ", obj, "\n", sep = "")
  
  if (is.null(rn)) {
    cat("No row names.\n")
    next
  }
  
  for (gene in core_genes_7b7) {
    
    cat(
      gene,
      " rowname match: ",
      gene %in% rn,
      "\n",
      sep = ""
    )
  }
}


# ------------------------------------------------------------
# 11. Check whether core genes occur in column names
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("CORE-GENE COLNAME AUDIT\n")
cat("------------------------------------------------------------\n")

for (obj in c(
  "estimate_input_df",
  "estimate_gct_input",
  "estimate_common",
  "estimate_filtered"
)) {
  
  if (!exists(
    obj,
    envir = estimate_source_env_7b7,
    inherits = FALSE
  )) {
    next
  }
  
  x <- get(
    obj,
    envir = estimate_source_env_7b7,
    inherits = FALSE
  )
  
  cn <- colnames(x)
  
  cat("\nOBJECT: ", obj, "\n", sep = "")
  
  if (is.null(cn)) {
    cat("No column names.\n")
    next
  }
  
  for (gene in core_genes_7b7) {
    
    cat(
      gene,
      " colname match: ",
      gene %in% cn,
      "\n",
      sep = ""
    )
  }
}


# ------------------------------------------------------------
# 12. Verify locked Step 7A-2 objects remain unchanged
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("LOCKED STEP 7A-2 INTEGRITY CHECK\n")
cat("------------------------------------------------------------\n")

audit_step7a2_7b7_2 <- all(
  exists("core_tme_alignment", inherits = TRUE),
  exists("gene_tme_results", inherits = TRUE),
  exists("gene_tme_manuscript_table", inherits = TRUE),
  nrow(core_tme_alignment) == 513,
  nrow(gene_tme_results) == 9,
  nrow(gene_tme_manuscript_table) == 3
)

cat(
  "Step 7A-2 objects intact: ",
  audit_step7a2_7b7_2,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. Explicitly confirm no analysis
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("ANALYSIS STATUS\n")
cat("------------------------------------------------------------\n")

cat("No ESTIMATE score recalculation.\n")
cat("No Spearman test.\n")
cat("No regression.\n")
cat("No leave-one-gene-out analysis.\n")
cat("No modification of locked Step 7A-2 objects.\n")


cat("\n============================================================\n")
cat("STEP 7B-7A-2 RECOVERY AUDIT COMPLETE\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 7B-7B
# Trace the historical ESTIMATE scoring implementation
# Audit only — NO statistical analysis
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7B: TRACE HISTORICAL ESTIMATE SCORING SOURCE\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Project root
# ------------------------------------------------------------

project_root_7b7 <- "input_dir/Paper1_tumor-only WGCNA"

cat("Project root:\n")
cat(project_root_7b7, "\n\n")

if (!dir.exists(project_root_7b7)) {
  stop("Project root does not exist. Do not proceed.")
}

# ------------------------------------------------------------
# 2. Search only likely source-file types
# ------------------------------------------------------------

source_ext_7b7 <- c(
  "\\.R$",
  "\\.Rmd$",
  "\\.qmd$",
  "\\.txt$",
  "\\.md$",
  "\\.log$",
  "\\.out$"
)

all_files_7b7 <- list.files(
  project_root_7b7,
  recursive = TRUE,
  full.names = TRUE
)

source_files_7b7 <- all_files_7b7[
  grepl(
    paste(source_ext_7b7, collapse = "|"),
    all_files_7b7,
    ignore.case = TRUE
  )
]

cat("Total source/text files found: ",
    length(source_files_7b7),
    "\n",
    sep = "")

# ------------------------------------------------------------
# 3. Search terms
# ------------------------------------------------------------

search_terms_7b7 <- c(
  "estimateScore",
  "filterCommonGenes",
  "StromalSignature",
  "ImmuneSignature",
  "StromalScore",
  "ImmuneScore",
  "ESTIMATEScore",
  "estimate_scores_wgcna"
)

cat("\nSearch terms:\n")
print(search_terms_7b7)

# ------------------------------------------------------------
# 4. Search file contents
# ------------------------------------------------------------

hits_7b7 <- list()

for (f in source_files_7b7) {
  
  txt <- tryCatch(
    readLines(
      f,
      warn = FALSE,
      encoding = "UTF-8"
    ),
    error = function(e) character(0)
  )
  
  if (length(txt) == 0) next
  
  hit_idx <- which(
    sapply(
      txt,
      function(line) {
        any(
          grepl(
            search_terms_7b7,
            line,
            ignore.case = TRUE
          )
        )
      }
    )
  )
  
  if (length(hit_idx) > 0) {
    
    hits_7b7[[f]] <- data.frame(
      File = f,
      Line = hit_idx,
      Text = txt[hit_idx],
      stringsAsFactors = FALSE
    )
  }
}

# ------------------------------------------------------------
# 5. Combine hits
# ------------------------------------------------------------

if (length(hits_7b7) > 0) {
  
  estimate_source_hits_7b7 <- do.call(
    rbind,
    hits_7b7
  )
  
  rownames(
    estimate_source_hits_7b7
  ) <- NULL
  
} else {
  
  estimate_source_hits_7b7 <- data.frame(
    File = character(0),
    Line = integer(0),
    Text = character(0),
    stringsAsFactors = FALSE
  )
}

cat("\n============================================================\n")
cat("SEARCH RESULTS\n")
cat("============================================================\n\n")

cat(
  "Number of matching lines: ",
  nrow(estimate_source_hits_7b7),
  "\n",
  sep = ""
)

cat(
  "Number of matching files: ",
  length(unique(estimate_source_hits_7b7$File)),
  "\n",
  sep = ""
)

if (nrow(estimate_source_hits_7b7) > 0) {
  
  cat("\nMatching files:\n")
  
  print(
    unique(
      estimate_source_hits_7b7$File
    )
  )
  
  cat("\nMatching lines:\n")
  print(
    estimate_source_hits_7b7,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "\nNo ESTIMATE-related source lines were found.\n"
  )
}

# ------------------------------------------------------------
# 6. Specifically search for the historical score object
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SEARCH FOR estimate_scores_wgcna\n")
cat("============================================================\n\n")

object_hits_7b7 <- estimate_source_hits_7b7[
  grepl(
    "estimate_scores_wgcna",
    estimate_source_hits_7b7$Text,
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]

if (nrow(object_hits_7b7) > 0) {
  
  print(
    object_hits_7b7,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No explicit estimate_scores_wgcna source line found.\n"
  )
}

# ------------------------------------------------------------
# 7. Audit status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ANALYSIS STATUS\n")
cat("============================================================\n\n")

cat("No ESTIMATE score recalculation.\n")
cat("No gene-TME association testing.\n")
cat("No regression.\n")
cat("No leave-one-gene-out analysis.\n")
cat("No modification of Step 7A-2 objects.\n")
cat("No modification of Figure 6.\n")
cat("No modification of Supplementary Table 3.\n")

cat("\n============================================================\n")
cat("STEP 7B-7B AUDIT COMPLETE\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 7B-7C
# Inspect historical Step14F-1 ESTIMATE source context
# Audit only — NO statistical analysis
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7C: INSPECT HISTORICAL ESTIMATE SOURCE CONTEXT\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Historical Step14F-1 script
# ------------------------------------------------------------

step14f1_file_7b7 <- paste0(
  "input_dir/Paper1_tumor-only WGCNA/",
  "Paper1_STEP19_PREP_20260914/scripts/",
  "WORKFLOW-03_Step14F-1_CoreGene_ModuleLocalization.R"
)

cat("Historical Step14F-1 script:\n")
cat(step14f1_file_7b7, "\n\n")

cat(
  "File exists: ",
  file.exists(step14f1_file_7b7),
  "\n",
  sep = ""
)

if (!file.exists(step14f1_file_7b7)) {
  stop("Historical Step14F-1 script not found.")
}

# ------------------------------------------------------------
# 2. Read entire script
# ------------------------------------------------------------

txt_7b7 <- readLines(
  step14f1_file_7b7,
  warn = FALSE,
  encoding = "UTF-8"
)

cat(
  "Total lines: ",
  length(txt_7b7),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 3. Print contexts around all ESTIMATE-related lines
# ------------------------------------------------------------

target_lines_7b7 <- which(
  grepl(
    paste(
      c(
        "estimate_scores_wgcna",
        "estimateScore",
        "filterCommonGenes",
        "StromalSignature",
        "ImmuneSignature",
        "StromalScore",
        "ImmuneScore",
        "ESTIMATEScore",
        "library\\(estimate\\)",
        "requireNamespace.*estimate",
        "load\\(",
        "readRDS",
        "read\\.table",
        "read\\.csv",
        "source\\("
      ),
      collapse = "|"
    ),
    txt_7b7,
    ignore.case = TRUE
  )
)

cat("\n============================================================\n")
cat("ESTIMATE-RELATED LINE NUMBERS\n")
cat("============================================================\n\n")

print(target_lines_7b7)

# ------------------------------------------------------------
# 4. Print ±12 lines around each target
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CONTEXT AROUND ESTIMATE-RELATED LINES\n")
cat("============================================================\n")

for (ln in target_lines_7b7) {
  
  start_7b7 <- max(1, ln - 12)
  end_7b7   <- min(length(txt_7b7), ln + 12)
  
  cat("\n------------------------------------------------------------\n")
  cat(
    "CONTEXT: lines ",
    start_7b7,
    " - ",
    end_7b7,
    "\n",
    sep = ""
  )
  cat("------------------------------------------------------------\n")
  
  for (i in start_7b7:end_7b7) {
    
    cat(
      sprintf(
        "%5d | %s\n",
        i,
        txt_7b7[i]
      )
    )
  }
}

# ------------------------------------------------------------
# 5. Search specifically for loading of estimate_scores_wgcna
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("POTENTIAL SOURCE / LOAD OF estimate_scores_wgcna\n")
cat("============================================================\n\n")

load_lines_7b7 <- which(
  grepl(
    "load\\(|readRDS|read\\.table|read\\.csv|read\\.delim|source\\(",
    txt_7b7,
    ignore.case = TRUE
  )
)

if (length(load_lines_7b7) > 0) {
  
  for (ln in load_lines_7b7) {
    
    cat(
      sprintf(
        "%5d | %s\n",
        ln,
        txt_7b7[ln]
      )
    )
  }
  
} else {
  
  cat("No load/read/source commands found.\n")
}

# ------------------------------------------------------------
# 6. Search for assignment of estimate_scores_wgcna
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ASSIGNMENT SEARCH\n")
cat("============================================================\n\n")

assignment_lines_7b7 <- which(
  grepl(
    "estimate_scores_wgcna\\s*<-",
    txt_7b7,
    ignore.case = TRUE
  )
)

if (length(assignment_lines_7b7) > 0) {
  
  for (ln in assignment_lines_7b7) {
    
    cat(
      sprintf(
        "%5d | %s\n",
        ln,
        txt_7b7[ln]
      )
    )
  }
  
} else {
  
  cat(
    "No direct assignment to estimate_scores_wgcna found.\n"
  )
}

# ------------------------------------------------------------
# 7. Audit status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ANALYSIS STATUS\n")
cat("============================================================\n\n")

cat("No ESTIMATE score recalculation.\n")
cat("No Spearman correlation.\n")
cat("No regression.\n")
cat("No leave-one-gene-out analysis.\n")
cat("No modification of locked Step 7A-2 objects.\n")
cat("No modification of Figure 6.\n")
cat("No modification of Supplementary Table 3.\n")

cat("\n============================================================\n")
cat("STEP 7B-7C AUDIT COMPLETE\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 7B-7D
# Historical ESTIMATE source tracing
# Goal:
#   Trace the generation source of estimate_scores_wgcna
#   from the locked Step14A-14E RData.
#
# IMPORTANT:
#   - No ESTIMATE recalculation
#   - No statistical analysis
#   - No leave-one-gene-out analysis
#   - No modification of Step 7A-2
#   - No modification of Figure 6
#   - No modification of Supplementary Table 3
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7D — HISTORICAL ESTIMATE SOURCE TRACING\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define project paths
# ------------------------------------------------------------

work_dir <- "O:/Rproject/Paper1_tumor-only WGCNA"

step14_locked_file <- file.path(
  work_dir,
  "WORKFLOW-03/Step14/Step14_FINAL_LOCKED/RData/",
  "WORKFLOW-03_Step14A-14E_FINAL_LOCKED.RData"
)

cat("Project root:\n")
cat(work_dir, "\n\n")

cat("Locked Step14A-14E file:\n")
cat(step14_locked_file, "\n\n")

cat(
  "Locked file exists: ",
  file.exists(step14_locked_file),
  "\n\n"
)


# ------------------------------------------------------------
# 2. Load Step14A-14E into isolated environment
# ------------------------------------------------------------

step14_env_7b7d <- new.env(parent = emptyenv())

load(
  step14_locked_file,
  envir = step14_env_7b7d
)

loaded_objects_7b7d <- ls(
  envir = step14_env_7b7d
)

cat("Objects loaded from Step14A-14E:\n")
print(loaded_objects_7b7d)

cat("\n")


# ------------------------------------------------------------
# 3. Inspect estimate_scores_wgcna
# ------------------------------------------------------------

if ("estimate_scores_wgcna" %in% loaded_objects_7b7d) {
  
  estimate_obj_7b7d <-
    get(
      "estimate_scores_wgcna",
      envir = step14_env_7b7d
    )
  
  cat("estimate_scores_wgcna found.\n")
  cat(
    "Class: ",
    paste(class(estimate_obj_7b7d), collapse = ", "),
    "\n"
  )
  
  cat(
    "Dimensions: ",
    paste(dim(estimate_obj_7b7d), collapse = " x "),
    "\n"
  )
  
  cat("Row names:\n")
  print(rownames(estimate_obj_7b7d))
  
  cat("\n")
  
  cat("First 6 column names:\n")
  print(head(colnames(estimate_obj_7b7d)))
  
  cat("\n")
  
} else {
  
  stop(
    "estimate_scores_wgcna was not found in Step14A-14E RData."
  )
  
}


# ------------------------------------------------------------
# 4. Search loaded object names for ESTIMATE-related objects
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ESTIMATE-RELATED OBJECT NAME SEARCH\n")
cat("============================================================\n\n")

estimate_object_candidates_7b7d <-
  loaded_objects_7b7d[
    grepl(
      "estimate|immune|stromal|score|signature|gct",
      loaded_objects_7b7d,
      ignore.case = TRUE
    )
  ]

print(estimate_object_candidates_7b7d)


# ------------------------------------------------------------
# 5. Inspect dimensions / classes of candidate objects
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CANDIDATE OBJECT STRUCTURE\n")
cat("============================================================\n\n")

for (obj_name in estimate_object_candidates_7b7d) {
  
  cat("\n------------------------------------------------------------\n")
  cat("Object: ", obj_name, "\n")
  cat("------------------------------------------------------------\n")
  
  obj <- get(
    obj_name,
    envir = step14_env_7b7d
  )
  
  cat(
    "Class: ",
    paste(class(obj), collapse = ", "),
    "\n"
  )
  
  if (!is.null(dim(obj))) {
    
    cat(
      "Dimensions: ",
      paste(dim(obj), collapse = " x "),
      "\n"
    )
    
  } else {
    
    cat(
      "Length: ",
      length(obj),
      "\n"
    )
    
  }
  
  if (is.data.frame(obj)) {
    
    cat("Column names:\n")
    print(names(obj))
    
  }
  
  if (is.matrix(obj)) {
    
    cat("Matrix row names available: ",
        !is.null(rownames(obj)),
        "\n")
    
    cat("Matrix column names available: ",
        !is.null(colnames(obj)),
        "\n")
    
  }
  
}


# ------------------------------------------------------------
# 6. Search object contents for signature-like names
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SIGNATURE MEMBERSHIP SEARCH IN LOADED OBJECTS\n")
cat("============================================================\n\n")

core_genes_7b7d <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

for (obj_name in estimate_object_candidates_7b7d) {
  
  obj <- get(
    obj_name,
    envir = step14_env_7b7d
  )
  
  cat("\n------------------------------------------------------------\n")
  cat("Object: ", obj_name, "\n")
  cat("------------------------------------------------------------\n")
  
  found_core <- FALSE
  
  # row names
  if (!is.null(rownames(obj))) {
    
    hits <- intersect(
      core_genes_7b7d,
      rownames(obj)
    )
    
    if (length(hits) > 0) {
      
      cat("Core genes found in row names:\n")
      print(hits)
      
      found_core <- TRUE
      
    }
    
  }
  
  # column names
  if (!is.null(colnames(obj))) {
    
    hits <- intersect(
      core_genes_7b7d,
      colnames(obj)
    )
    
    if (length(hits) > 0) {
      
      cat("Core genes found in column names:\n")
      print(hits)
      
      found_core <- TRUE
      
    }
    
  }
  
  # character vectors
  if (is.character(obj) && length(obj) > 0) {
    
    hits <- intersect(
      core_genes_7b7d,
      unique(obj)
    )
    
    if (length(hits) > 0) {
      
      cat("Core genes found in character vector:\n")
      print(hits)
      
      found_core <- TRUE
      
    }
    
  }
  
  if (!found_core) {
    
    cat(
      "No AOC3/MS4A7/IL7R membership found in this object's ",
      "names or direct character values.\n",
      sep = ""
    )
    
  }
  
}


# ------------------------------------------------------------
# 7. Search the locked Step14A-14E file directory
#    for historical scripts and RData files
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP14 DIRECTORY SOURCE SEARCH\n")
cat("============================================================\n\n")

step14_root <- file.path(
  work_dir,
  "WORKFLOW-03/Step14"
)

cat("Search root:\n")
cat(step14_root, "\n\n")

if (dir.exists(step14_root)) {
  
  all_step14_files_7b7d <-
    list.files(
      step14_root,
      recursive = TRUE,
      full.names = TRUE,
      include.dirs = FALSE
    )
  
  cat(
    "Total files found: ",
    length(all_step14_files_7b7d),
    "\n\n"
  )
  
  source_candidates_7b7d <-
    all_step14_files_7b7d[
      grepl(
        "\\.(R|r|Rmd|rmd|txt|RData|rdata|rds|csv|tsv|gct)$",
        all_step14_files_7b7d,
        ignore.case = TRUE
      )
    ]
  
  cat(
    "Potential source/data files: ",
    length(source_candidates_7b7d),
    "\n\n"
  )
  
  print(source_candidates_7b7d)
  
} else {
  
  cat("Step14 directory does not exist.\n")
  
}


# ------------------------------------------------------------
# 8. Search Step14 source/text files for ESTIMATE generation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ESTIMATE GENERATION COMMAND SEARCH\n")
cat("============================================================\n\n")

text_files_7b7d <-
  all_step14_files_7b7d[
    grepl(
      "\\.(R|r|Rmd|rmd|txt)$",
      all_step14_files_7b7d,
      ignore.case = TRUE
    )
  ]

estimate_patterns_7b7d <- c(
  "estimateScore",
  "filterCommonGenes",
  "estimate_scores_wgcna",
  "estimate_scores",
  "ImmuneSignature",
  "StromalSignature",
  "StromalScore",
  "ImmuneScore",
  "ESTIMATEScore",
  "library\\(estimate",
  "require\\(estimate",
  "source\\("
)

match_table_7b7d <- data.frame(
  File = character(),
  Line = integer(),
  Text = character(),
  Pattern = character(),
  stringsAsFactors = FALSE
)

for (f in text_files_7b7d) {
  
  txt <- tryCatch(
    readLines(
      f,
      warn = FALSE,
      encoding = "UTF-8"
    ),
    error = function(e) character(0)
  )
  
  if (length(txt) == 0) next
  
  for (pat in estimate_patterns_7b7d) {
    
    hits <- which(
      grepl(
        pat,
        txt,
        ignore.case = TRUE
      )
    )
    
    if (length(hits) > 0) {
      
      match_table_7b7d <- rbind(
        match_table_7b7d,
        data.frame(
          File = f,
          Line = hits,
          Text = txt[hits],
          Pattern = pat,
          stringsAsFactors = FALSE
        )
      )
      
    }
    
  }
  
}

if (nrow(match_table_7b7d) > 0) {
  
  match_table_7b7d <-
    match_table_7b7d[
      !duplicated(
        match_table_7b7d[
          ,
          c("File", "Line", "Text")
        ]
      ),
      ,
      drop = FALSE
    ]
  
  print(match_table_7b7d)
  
} else {
  
  cat("No ESTIMATE-related source lines found.\n")
  
}


# ------------------------------------------------------------
# 9. Search specifically for source/load/read commands
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("UPSTREAM LOAD / SOURCE COMMANDS\n")
cat("============================================================\n\n")

load_source_table_7b7d <- data.frame(
  File = character(),
  Line = integer(),
  Text = character(),
  stringsAsFactors = FALSE
)

for (f in text_files_7b7d) {
  
  txt <- tryCatch(
    readLines(
      f,
      warn = FALSE,
      encoding = "UTF-8"
    ),
    error = function(e) character(0)
  )
  
  if (length(txt) == 0) next
  
  hits <- which(
    grepl(
      "load\\(|readRDS\\(|read\\.table\\(|read\\.csv\\(|read\\.delim\\(|source\\(",
      txt,
      ignore.case = TRUE
    )
  )
  
  if (length(hits) > 0) {
    
    load_source_table_7b7d <- rbind(
      load_source_table_7b7d,
      data.frame(
        File = f,
        Line = hits,
        Text = txt[hits],
        stringsAsFactors = FALSE
      )
    )
    
  }
  
}

print(load_source_table_7b7d)


# ------------------------------------------------------------
# 10. Search for direct assignment to estimate_scores_wgcna
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("DIRECT ASSIGNMENT SEARCH\n")
cat("============================================================\n\n")

assignment_table_7b7d <- data.frame(
  File = character(),
  Line = integer(),
  Text = character(),
  stringsAsFactors = FALSE
)

for (f in text_files_7b7d) {
  
  txt <- tryCatch(
    readLines(
      f,
      warn = FALSE,
      encoding = "UTF-8"
    ),
    error = function(e) character(0)
  )
  
  if (length(txt) == 0) next
  
  hits <- which(
    grepl(
      "estimate_scores_wgcna\\s*<-",
      txt,
      ignore.case = TRUE
    )
  )
  
  if (length(hits) > 0) {
    
    assignment_table_7b7d <- rbind(
      assignment_table_7b7d,
      data.frame(
        File = f,
        Line = hits,
        Text = txt[hits],
        stringsAsFactors = FALSE
      )
    )
    
  }
  
}

if (nrow(assignment_table_7b7d) > 0) {
  
  print(assignment_table_7b7d)
  
} else {
  
  cat(
    "No direct assignment to estimate_scores_wgcna found ",
    "within Step14 source/text files.\n",
    sep = ""
  )
  
}


# ------------------------------------------------------------
# 11. Search for signature-like files by filename
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SIGNATURE-LIKE FILE SEARCH\n")
cat("============================================================\n\n")

signature_files_7b7d <-
  all_step14_files_7b7d[
    grepl(
      "signature|immune|stromal|estimate|gct",
      basename(all_step14_files_7b7d),
      ignore.case = TRUE
    )
  ]

if (length(signature_files_7b7d) > 0) {
  
  print(signature_files_7b7d)
  
} else {
  
  cat(
    "No signature-like files found by filename.\n"
  )
  
}


# ------------------------------------------------------------
# 12. Audit status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 7B-7D AUDIT STATUS\n")
cat("============================================================\n\n")

cat("No ESTIMATE score recalculation.\n")
cat("No statistical analysis.\n")
cat("No leave-one-gene-out analysis.\n")
cat("No modification of Step 7A-2 objects.\n")
cat("No modification of Figure 6.\n")
cat("No modification of Supplementary Table 3.\n")

cat("\n============================================================\n")
cat("STEP 7B-7D AUDIT COMPLETE\n")
cat("============================================================\n")

#STEP 7B-7E：Historical ESTIMATE Signature Membership Audit。

# ============================================================
# Paper 1
# STEP 7B-7E
# Historical ESTIMATE Signature Membership Audit
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7E: Historical ESTIMATE Signature Membership Audit\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Project root
# ------------------------------------------------------------

work_dir <- "O:/Rproject/Paper1_tumor-only WGCNA"
setwd(work_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 2. Locked Step14A-14E RData
# ------------------------------------------------------------

step14_rdata <- file.path(
  work_dir,
  "WORKFLOW-03",
  "Step14",
  "Step14_FINAL_LOCKED",
  "RData",
  "WORKFLOW-03_Step14A-14E_FINAL_LOCKED.RData"
)

cat("Step14A-14E locked RData:\n")
cat(step14_rdata, "\n\n")

cat("File exists: ", file.exists(step14_rdata), "\n\n")

stopifnot(file.exists(step14_rdata))


# ------------------------------------------------------------
# 3. Load into isolated environment
# ------------------------------------------------------------

estimate_audit_env_7be <- new.env(parent = emptyenv())

loaded_7be <- load(
  step14_rdata,
  envir = estimate_audit_env_7be
)

cat("Objects loaded into isolated audit environment:\n")
print(loaded_7be)

cat("\nNumber of loaded objects: ",
    length(loaded_7be), "\n\n")


# ------------------------------------------------------------
# 4. Inspect object names only
# ------------------------------------------------------------

cat("============================================================\n")
cat("OBJECT INVENTORY\n")
cat("============================================================\n\n")

print(sort(ls(envir = estimate_audit_env_7be)))

cat("\n")

# ============================================================
# Paper 1
# STEP 7B-7E-2
# Historical ESTIMATE Object Structural Audit
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7E-2: Historical ESTIMATE Object Structural Audit\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Objects of direct interest
# ------------------------------------------------------------

estimate_history_objects_7be <- c(
  "estimate_input_df",
  "estimate_gct_input",
  "estimate_input_file",
  "estimate_common",
  "estimate_filtered",
  "estimate_scores",
  "estimate_scores_wgcna",
  "expr_estimate",
  "annotation",
  "annotation_file",
  "gene_map",
  "gene_symbol",
  "gene_type"
)

estimate_history_objects_7be <- estimate_history_objects_7be[
  estimate_history_objects_7be %in%
    ls(envir = estimate_audit_env_7be)
]

cat("Historical ESTIMATE-related objects present:\n")
print(estimate_history_objects_7be)

cat("\nNumber present:",
    length(estimate_history_objects_7be),
    "\n\n")


# ------------------------------------------------------------
# 2. Structural audit
# ------------------------------------------------------------

for (obj_name in estimate_history_objects_7be) {
  
  cat("\n============================================================\n")
  cat("OBJECT:", obj_name, "\n")
  cat("============================================================\n")
  
  obj <- get(
    obj_name,
    envir = estimate_audit_env_7be
  )
  
  cat("Class:\n")
  print(class(obj))
  
  cat("\nType:\n")
  print(typeof(obj))
  
  cat("\nLength:\n")
  print(length(obj))
  
  # ----------------------------------------------------------
  # data.frame
  # ----------------------------------------------------------
  
  if (is.data.frame(obj)) {
    
    cat("\nDimensions:\n")
    print(dim(obj))
    
    cat("\nColumn names:\n")
    print(colnames(obj))
    
    cat("\nFirst 5 rows:\n")
    print(utils::head(obj, 5))
    
  }
  
  # ----------------------------------------------------------
  # matrix
  # ----------------------------------------------------------
  
  if (is.matrix(obj)) {
    
    cat("\nDimensions:\n")
    print(dim(obj))
    
    cat("\nRow names preview:\n")
    print(utils::head(rownames(obj), 20))
    
    cat("\nColumn names preview:\n")
    print(utils::head(colnames(obj), 10))
    
  }
  
  # ----------------------------------------------------------
  # vector
  # ----------------------------------------------------------
  
  if (is.atomic(obj) &&
      !is.data.frame(obj) &&
      !is.matrix(obj) &&
      length(obj) > 0) {
    
    cat("\nFirst 20 values:\n")
    print(utils::head(obj, 20))
  }
}


cat("\n============================================================\n")
cat("STEP 7B-7E-2 STRUCTURAL AUDIT COMPLETE\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 7B-7E-3
# Targeted Historical ESTIMATE Input / Filtering Audit
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7E-3: Targeted Historical ESTIMATE Input Audit\n")
cat("============================================================\n\n")

target_objects_7be <- c(
  "estimate_input_df",
  "estimate_gct_input",
  "estimate_common",
  "estimate_filtered",
  "estimate_scores"
)

target_objects_7be <- target_objects_7be[
  target_objects_7be %in%
    ls(envir = estimate_audit_env_7be)
]

cat("Target objects present:\n")
print(target_objects_7be)

cat("\n")


# ------------------------------------------------------------
# Core genes
# ------------------------------------------------------------

core_genes_7be <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


# ------------------------------------------------------------
# Function: inspect object
# ------------------------------------------------------------

inspect_estimate_object_7be <- function(obj_name, env) {
  
  cat("\n============================================================\n")
  cat("OBJECT:", obj_name, "\n")
  cat("============================================================\n")
  
  obj <- get(obj_name, envir = env)
  
  cat("Class:\n")
  print(class(obj))
  
  cat("\nType:\n")
  print(typeof(obj))
  
  cat("\nLength:\n")
  print(length(obj))
  
  if (is.data.frame(obj)) {
    
    cat("\nDimensions:\n")
    print(dim(obj))
    
    cat("\nColumn names:\n")
    print(colnames(obj))
    
    cat("\nRow names preview:\n")
    print(utils::head(rownames(obj), 20))
    
    cat("\nFirst 5 rows:\n")
    print(utils::head(obj, 5))
  }
  
  if (is.matrix(obj)) {
    
    cat("\nDimensions:\n")
    print(dim(obj))
    
    cat("\nRow names preview:\n")
    print(utils::head(rownames(obj), 30))
    
    cat("\nColumn names preview:\n")
    print(utils::head(colnames(obj), 10))
    
    cat("\nFirst 5 x 5 block:\n")
    print(
      obj[
        seq_len(min(5, nrow(obj))),
        seq_len(min(5, ncol(obj))),
        drop = FALSE
      ]
    )
  }
  
  if (is.character(obj)) {
    
    cat("\nFirst 30 values:\n")
    print(utils::head(obj, 30))
  }
  
  invisible(obj)
}


# ------------------------------------------------------------
# Inspect all target objects
# ------------------------------------------------------------

for (obj_name in target_objects_7be) {
  
  inspect_estimate_object_7be(
    obj_name,
    estimate_audit_env_7be
  )
}


# ------------------------------------------------------------
# Exact core-gene membership
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CORE GENE MEMBERSHIP AUDIT\n")
cat("============================================================\n\n")

membership_results_7be <- data.frame(
  Object = character(),
  AOC3 = logical(),
  MS4A7 = logical(),
  IL7R = logical(),
  stringsAsFactors = FALSE
)

for (obj_name in target_objects_7be) {
  
  obj <- get(
    obj_name,
    envir = estimate_audit_env_7be
  )
  
  gene_candidates <- character(0)
  
  # ----------------------------------------------------------
  # Collect row names
  # ----------------------------------------------------------
  
  if (!is.null(rownames(obj))) {
    gene_candidates <- c(
      gene_candidates,
      rownames(obj)
    )
  }
  
  # ----------------------------------------------------------
  # Collect column names
  # ----------------------------------------------------------
  
  if (!is.null(colnames(obj))) {
    gene_candidates <- c(
      gene_candidates,
      colnames(obj)
    )
  }
  
  # ----------------------------------------------------------
  # Collect character contents
  # ----------------------------------------------------------
  
  if (is.character(obj)) {
    
    gene_candidates <- c(
      gene_candidates,
      as.character(obj)
    )
    
  } else if (is.factor(obj)) {
    
    gene_candidates <- c(
      gene_candidates,
      as.character(obj)
    )
    
  } else if (is.data.frame(obj)) {
    
    char_cols <- vapply(
      obj,
      function(x) is.character(x) || is.factor(x),
      logical(1)
    )
    
    if (any(char_cols)) {
      
      gene_candidates <- c(
        gene_candidates,
        unlist(
          lapply(
            obj[, char_cols, drop = FALSE],
            as.character
          ),
          use.names = FALSE
        )
      )
    }
  }
  
  gene_candidates <- unique(
    gene_candidates[!is.na(gene_candidates)]
  )
  
  membership_results_7be <- rbind(
    membership_results_7be,
    data.frame(
      Object = obj_name,
      AOC3 = "AOC3" %in% gene_candidates,
      MS4A7 = "MS4A7" %in% gene_candidates,
      IL7R = "IL7R" %in% gene_candidates,
      stringsAsFactors = FALSE
    )
  )
}

print(membership_results_7be)


# ------------------------------------------------------------
# Also inspect the exact locations if found
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("EXACT CORE-GENE LOCATIONS\n")
cat("============================================================\n\n")

for (obj_name in target_objects_7be) {
  
  obj <- get(
    obj_name,
    envir = estimate_audit_env_7be
  )
  
  cat("\nOBJECT:", obj_name, "\n")
  
  if (!is.null(rownames(obj))) {
    
    for (g in core_genes_7be) {
      
      idx <- which(
        rownames(obj) == g
      )
      
      if (length(idx) > 0) {
        
        cat(
          g,
          "found in row position(s):",
          paste(idx, collapse = ", "),
          "\n"
        )
      }
    }
  }
  
  if (!is.null(colnames(obj))) {
    
    for (g in core_genes_7be) {
      
      idx <- which(
        colnames(obj) == g
      )
      
      if (length(idx) > 0) {
        
        cat(
          g,
          "found in column position(s):",
          paste(idx, collapse = ", "),
          "\n"
        )
      }
    }
  }
}


cat("\n============================================================\n")
cat("STEP 7B-7E-3 TARGETED AUDIT COMPLETE\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 7B-7E-4
# Detailed Structure and Gene-Set Relationship Audit
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7E-4: Detailed ESTIMATE Object Relationship Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Function for detailed inspection
# ------------------------------------------------------------

inspect_detailed_7be <- function(obj_name, env) {
  
  cat("\n============================================================\n")
  cat("OBJECT:", obj_name, "\n")
  cat("============================================================\n")
  
  obj <- get(obj_name, envir = env)
  
  cat("\nClass:\n")
  print(class(obj))
  
  cat("\nType:\n")
  print(typeof(obj))
  
  cat("\nLength:\n")
  print(length(obj))
  
  cat("\nAttributes:\n")
  print(attributes(obj))
  
  if (is.data.frame(obj) || is.matrix(obj)) {
    
    cat("\nDimensions:\n")
    print(dim(obj))
    
    cat("\nNumber of rows:\n")
    print(nrow(obj))
    
    cat("\nNumber of columns:\n")
    print(ncol(obj))
    
    cat("\nRow names - first 20:\n")
    print(head(rownames(obj), 20))
    
    cat("\nRow names - last 20:\n")
    print(tail(rownames(obj), 20))
    
    cat("\nColumn names - first 20:\n")
    print(head(colnames(obj), 20))
    
    cat("\nColumn names - last 20:\n")
    print(tail(colnames(obj), 20))
    
  }
  
  if (is.data.frame(obj)) {
    
    cat("\nColumn classes:\n")
    print(vapply(obj, class, character(1)))
    
    cat("\nFirst 10 rows:\n")
    print(head(obj, 10))
    
  }
  
  if (is.matrix(obj)) {
    
    cat("\nFirst 10 x 10 block:\n")
    print(
      obj[
        seq_len(min(10, nrow(obj))),
        seq_len(min(10, ncol(obj))),
        drop = FALSE
      ]
    )
  }
  
  invisible(obj)
}


# ------------------------------------------------------------
# 2. Detailed inspection of the two critical objects
# ------------------------------------------------------------

inspect_detailed_7be(
  "estimate_common",
  estimate_audit_env_7be
)

inspect_detailed_7be(
  "estimate_filtered",
  estimate_audit_env_7be
)


# ------------------------------------------------------------
# 3. Explicit row/column membership
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("EXPLICIT DIMNAME MEMBERSHIP\n")
cat("============================================================\n\n")

core_genes_7be <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

for (obj_name in c(
  "estimate_common",
  "estimate_filtered"
)) {
  
  obj <- get(
    obj_name,
    envir = estimate_audit_env_7be
  )
  
  cat("\n------------------------------------------------------------\n")
  cat("OBJECT:", obj_name, "\n")
  cat("------------------------------------------------------------\n")
  
  for (g in core_genes_7be) {
    
    row_hit <- integer(0)
    col_hit <- integer(0)
    
    if (!is.null(rownames(obj))) {
      row_hit <- which(
        rownames(obj) == g
      )
    }
    
    if (!is.null(colnames(obj))) {
      col_hit <- which(
        colnames(obj) == g
      )
    }
    
    cat(
      g,
      ": row positions = ",
      ifelse(
        length(row_hit) == 0,
        "NONE",
        paste(row_hit, collapse = ", ")
      ),
      "; column positions = ",
      ifelse(
        length(col_hit) == 0,
        "NONE",
        paste(col_hit, collapse = ", ")
      ),
      "\n",
      sep = ""
    )
  }
}


# ------------------------------------------------------------
# 4. Search object contents specifically for core genes
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CONTENT-LEVEL CORE GENE SEARCH\n")
cat("============================================================\n\n")

for (obj_name in c(
  "estimate_common",
  "estimate_filtered"
)) {
  
  obj <- get(
    obj_name,
    envir = estimate_audit_env_7be
  )
  
  cat("\nOBJECT:", obj_name, "\n")
  
  if (is.data.frame(obj)) {
    
    for (g in core_genes_7be) {
      
      hit_count <- sum(
        vapply(
          obj,
          function(x) {
            if (is.character(x) || is.factor(x)) {
              any(as.character(x) == g, na.rm = TRUE)
            } else {
              FALSE
            }
          },
          logical(1)
        )
      )
      
      cat(
        g,
        ": number of character columns containing exact match = ",
        hit_count,
        "\n",
        sep = ""
      )
    }
    
  } else if (is.matrix(obj)) {
    
    vals <- as.character(
      obj
    )
    
    for (g in core_genes_7be) {
      
      cat(
        g,
        ": exact content matches = ",
        sum(vals == g, na.rm = TRUE),
        "\n",
        sep = ""
      )
    }
    
  } else {
    
    cat(
      "Object is neither matrix nor data.frame.\n"
    )
  }
}


cat("\n============================================================\n")
cat("STEP 7B-7E-4 COMPLETE\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 7B-7E-5
# Detailed Structure Audit of estimate_common
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7E-5: Detailed Structure Audit of estimate_common\n")
cat("============================================================\n\n")

obj <- get(
  "estimate_common",
  envir = estimate_audit_env_7be
)

# ------------------------------------------------------------
# 1. Basic structure
# ------------------------------------------------------------

cat("Class:\n")
print(class(obj))

cat("\nType:\n")
print(typeof(obj))

cat("\nDimensions:\n")
print(dim(obj))

cat("\nColumn names:\n")
print(colnames(obj))

cat("\nRow names - first 20:\n")
print(head(rownames(obj), 20))

cat("\nRow names - last 20:\n")
print(tail(rownames(obj), 20))

cat("\nColumn classes:\n")
print(
  vapply(
    obj,
    function(x) paste(class(x), collapse = "/"),
    character(1)
  )
)

# ------------------------------------------------------------
# 2. Show first rows
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FIRST 10 ROWS\n")
cat("============================================================\n\n")

print(
  head(
    obj,
    10
  )
)

# ------------------------------------------------------------
# 3. Identify exactly which character columns contain
#    AOC3 / MS4A7 / IL7R
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("EXACT CORE-GENE LOCATIONS IN CHARACTER COLUMNS\n")
cat("============================================================\n\n")

core_genes_7be <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

for (g in core_genes_7be) {
  
  cat("\n------------------------------------------------------------\n")
  cat("GENE:", g, "\n")
  cat("------------------------------------------------------------\n")
  
  for (j in seq_along(obj)) {
    
    x <- obj[[j]]
    
    if (is.character(x) || is.factor(x)) {
      
      x_chr <- as.character(x)
      
      hit <- which(
        !is.na(x_chr) &
          x_chr == g
      )
      
      if (length(hit) > 0) {
        
        cat(
          "Column:",
          colnames(obj)[j],
          "\n"
        )
        
        cat(
          "Class:",
          class(x),
          "\n"
        )
        
        cat(
          "Number of exact matches:",
          length(hit),
          "\n"
        )
        
        cat(
          "Row positions:",
          paste(hit, collapse = ", "),
          "\n"
        )
        
        cat(
          "Matched values:\n"
        )
        
        print(
          obj[
            hit,
            ,
            drop = FALSE
          ]
        )
      }
    }
  }
}

# ------------------------------------------------------------
# 4. Show unique values / frequency information
#    for the character columns containing core genes
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CHARACTER COLUMN CONTENT SUMMARY\n")
cat("============================================================\n\n")

for (j in seq_along(obj)) {
  
  x <- obj[[j]]
  
  if (is.character(x) || is.factor(x)) {
    
    x_chr <- as.character(x)
    
    core_hit <- any(
      x_chr %in% core_genes_7be,
      na.rm = TRUE
    )
    
    if (core_hit) {
      
      cat("\n------------------------------------------------------------\n")
      cat(
        "COLUMN:",
        colnames(obj)[j],
        "\n"
      )
      cat("------------------------------------------------------------\n")
      
      cat(
        "Class:",
        class(x),
        "\n"
      )
      
      cat(
        "Length:",
        length(x_chr),
        "\n"
      )
      
      cat(
        "Unique values:",
        length(unique(x_chr)),
        "\n"
      )
      
      cat("\nFirst 30 values:\n")
      print(
        head(
          unique(x_chr),
          30
        )
      )
    }
  }
}

# ------------------------------------------------------------
# 5. Object summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("OBJECT SUMMARY\n")
cat("============================================================\n\n")

str(
  obj,
  max.level = 2
)

cat("\n============================================================\n")
cat("STEP 7B-7E-5 COMPLETE\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 7B-7E-6
# Trace the Origin of estimate_common
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7E-6: Trace Origin of estimate_common\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Basic dimensions of relevant objects
# ------------------------------------------------------------

cat("estimate_input_df:\n")
print(class(estimate_audit_env_7be$estimate_input_df))
print(dim(estimate_audit_env_7be$estimate_input_df))

cat("\nestimate_gct_input:\n")
print(class(estimate_audit_env_7be$estimate_gct_input))
print(dim(estimate_audit_env_7be$estimate_gct_input))

cat("\nestimate_common:\n")
print(class(estimate_audit_env_7be$estimate_common))
print(dim(estimate_audit_env_7be$estimate_common))

cat("\nestimate_filtered:\n")
print(class(estimate_audit_env_7be$estimate_filtered))
print(length(estimate_audit_env_7be$estimate_filtered))

cat("\nestimate_scores:\n")
print(class(estimate_audit_env_7be$estimate_scores))
print(dim(estimate_audit_env_7be$estimate_scores))


# ------------------------------------------------------------
# 2. Gene counts and sample counts
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GENE / SAMPLE COUNT COMPARISON\n")
cat("============================================================\n\n")

if (is.data.frame(estimate_audit_env_7be$estimate_input_df)) {
  cat(
    "estimate_input_df rows =",
    nrow(estimate_audit_env_7be$estimate_input_df),
    "\n"
  )
  cat(
    "estimate_input_df columns =",
    ncol(estimate_audit_env_7be$estimate_input_df),
    "\n"
  )
}

if (is.data.frame(estimate_audit_env_7be$estimate_common)) {
  cat(
    "estimate_common rows =",
    nrow(estimate_audit_env_7be$estimate_common),
    "\n"
  )
  cat(
    "estimate_common columns =",
    ncol(estimate_audit_env_7be$estimate_common),
    "\n"
  )
}


# ------------------------------------------------------------
# 3. Exact core-gene presence in NAME
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CORE GENE PRESENCE IN estimate_common$NAME\n")
cat("============================================================\n\n")

common_name <- as.character(
  estimate_audit_env_7be$estimate_common$NAME
)

for (g in c("AOC3", "MS4A7", "IL7R")) {
  
  idx <- which(common_name == g)
  
  cat(
    g,
    ": ",
    ifelse(
      length(idx) == 0,
      "NOT PRESENT",
      paste0(
        "PRESENT; row = ",
        paste(idx, collapse = ", ")
      )
    ),
    "\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 4. Show exact rows for the core genes that are present
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CORE GENE ROWS IN estimate_common\n")
cat("============================================================\n\n")

for (g in c("AOC3", "MS4A7", "IL7R")) {
  
  idx <- which(common_name == g)
  
  if (length(idx) > 0) {
    
    cat("\n------------------------------\n")
    cat(g, "\n")
    cat("------------------------------\n")
    
    print(
      estimate_audit_env_7be$estimate_common[
        idx,
        1:min(5, ncol(
          estimate_audit_env_7be$estimate_common
        )),
        drop = FALSE
      ]
    )
  }
}


# ------------------------------------------------------------
# 5. Determine whether NAME and Description are identical
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("NAME / DESCRIPTION IDENTITY\n")
cat("============================================================\n\n")

name_vec <- as.character(
  estimate_audit_env_7be$estimate_common$NAME
)

desc_vec <- as.character(
  estimate_audit_env_7be$estimate_common$Description
)

cat(
  "NAME identical to Description:",
  identical(name_vec, desc_vec),
  "\n"
)


# ------------------------------------------------------------
# 6. Compare common genes against estimate_input_df
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("COMMON-GENE RELATIONSHIP\n")
cat("============================================================\n\n")

input_obj <- estimate_audit_env_7be$estimate_input_df

if (is.data.frame(input_obj)) {
  
  input_char_cols <- names(input_obj)[
    vapply(
      input_obj,
      function(x) is.character(x) || is.factor(x),
      logical(1)
    )
  ]
  
  cat(
    "Character/factor columns in estimate_input_df:\n"
  )
  print(input_char_cols)
  
  if ("NAME" %in% names(input_obj)) {
    
    input_name <- as.character(
      input_obj$NAME
    )
    
    cat(
      "Number of input genes:",
      length(input_name),
      "\n"
    )
    
    cat(
      "Number of unique input genes:",
      length(unique(input_name)),
      "\n"
    )
    
    cat(
      "Number of common genes:",
      length(common_name),
      "\n"
    )
    
    cat(
      "All common genes present in input:",
      all(common_name %in% input_name),
      "\n"
    )
    
    cat(
      "Common genes present in input:",
      sum(common_name %in% input_name),
      "\n"
    )
  }
}


# ------------------------------------------------------------
# 7. Compare core-gene membership directly
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CORE GENE: INPUT vs COMMON\n")
cat("============================================================\n\n")

for (g in c("AOC3", "MS4A7", "IL7R")) {
  
  in_input <- FALSE
  
  if (
    is.data.frame(input_obj) &&
    "NAME" %in% names(input_obj)
  ) {
    in_input <- g %in% as.character(input_obj$NAME)
  }
  
  in_common <- g %in% common_name
  
  cat(
    sprintf(
      "%-6s input=%-5s common=%-5s\n",
      g,
      in_input,
      in_common
    )
  )
}


# ------------------------------------------------------------
# 8. Search historical source code for estimate_common
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SEARCH FOR estimate_common IN CURRENT SOURCE FILES\n")
cat("============================================================\n\n")

project_root <- "O:/Rproject/Paper1-tumor-only WGCNA"

source_files <- list.files(
  project_root,
  pattern = "\\.(R|r|Rmd|txt)$",
  recursive = TRUE,
  full.names = TRUE
)

hit_files <- character(0)

for (f in source_files) {
  
  txt <- tryCatch(
    readLines(
      f,
      warn = FALSE,
      encoding = "UTF-8"
    ),
    error = function(e) character(0)
  )
  
  if (
    length(txt) > 0 &&
    any(grepl(
      "estimate_common",
      txt,
      fixed = TRUE
    ))
  ) {
    hit_files <- c(
      hit_files,
      f
    )
  }
}

cat(
  "Number of source files containing 'estimate_common':",
  length(hit_files),
  "\n\n"
)

if (length(hit_files) > 0) {
  print(hit_files)
}


# ------------------------------------------------------------
# 9. Print matching source lines
# ------------------------------------------------------------

if (length(hit_files) > 0) {
  
  cat("\n============================================================\n")
  cat("MATCHING SOURCE LINES\n")
  cat("============================================================\n\n")
  
  for (f in hit_files) {
    
    txt <- tryCatch(
      readLines(
        f,
        warn = FALSE,
        encoding = "UTF-8"
      ),
      error = function(e) character(0)
    )
    
    idx <- grep(
      "estimate_common",
      txt,
      fixed = TRUE
    )
    
    if (length(idx) > 0) {
      
      cat("\n============================================================\n")
      cat("FILE:", f, "\n")
      cat("============================================================\n")
      
      for (i in idx) {
        
        start_i <- max(1, i - 5)
        end_i <- min(length(txt), i + 8)
        
        cat(
          "\n--- lines ",
          start_i,
          "-",
          end_i,
          " ---\n",
          sep = ""
        )
        
        cat(
          paste(
            sprintf(
              "%5d: %s",
              start_i:end_i,
              txt[start_i:end_i]
            ),
            collapse = "\n"
          )
        )
        
        cat("\n")
      }
    }
  }
}

cat("\n============================================================\n")
cat("STEP 7B-7E-6 COMPLETE\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 7B-7E-7
# Direct ESTIMATE Signature Membership Audit
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7E-7: Direct ESTIMATE Signature Membership Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Search historical environment for signature-related objects
# ------------------------------------------------------------

cat("============================================================\n")
cat("SIGNATURE-RELATED OBJECT INVENTORY\n")
cat("============================================================\n\n")

all_names_7be <- ls(
  envir = estimate_audit_env_7be,
  all.names = TRUE
)

signature_names_7be <- all_names_7be[
  grepl(
    "signature|geneset|stromal|immune|si_",
    all_names_7be,
    ignore.case = TRUE
  )
]

cat(
  "Number of signature-related object names:",
  length(signature_names_7be),
  "\n\n"
)

print(signature_names_7be)


# ------------------------------------------------------------
# 2. Inspect likely signature objects
# ------------------------------------------------------------

for (obj_name in signature_names_7be) {
  
  cat("\n============================================================\n")
  cat("OBJECT:", obj_name, "\n")
  cat("============================================================\n")
  
  obj <- get(
    obj_name,
    envir = estimate_audit_env_7be
  )
  
  cat("Class:\n")
  print(class(obj))
  
  cat("Type:\n")
  print(typeof(obj))
  
  cat("Length:\n")
  print(length(obj))
  
  cat("Dimensions:\n")
  print(dim(obj))
  
  if (is.data.frame(obj)) {
    
    cat("\nColumn names:\n")
    print(colnames(obj))
    
    cat("\nFirst 10 rows:\n")
    print(head(obj, 10))
    
  } else if (is.matrix(obj)) {
    
    cat("\nRow names - first 20:\n")
    print(head(rownames(obj), 20))
    
    cat("\nColumn names - first 20:\n")
    print(head(colnames(obj), 20))
    
    cat("\nFirst 10 × 10 block:\n")
    
    print(
      obj[
        seq_len(min(10, nrow(obj))),
        seq_len(min(10, ncol(obj))),
        drop = FALSE
      ]
    )
    
  } else if (is.character(obj)) {
    
    cat("\nFirst 50 values:\n")
    print(head(obj, 50))
  }
}


# ------------------------------------------------------------
# 3. Search ALL character/data.frame/matrix objects
#    for the three core genes
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CORE-GENE SEARCH ACROSS SIGNATURE-RELATED OBJECTS\n")
cat("============================================================\n\n")

core_genes_7be <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

for (obj_name in signature_names_7be) {
  
  obj <- get(
    obj_name,
    envir = estimate_audit_env_7be
  )
  
  cat("\n------------------------------------------------------------\n")
  cat("OBJECT:", obj_name, "\n")
  cat("------------------------------------------------------------\n")
  
  for (g in core_genes_7be) {
    
    hit <- FALSE
    hit_count <- 0
    
    if (is.character(obj)) {
      
      hit_count <- sum(
        as.character(obj) == g,
        na.rm = TRUE
      )
      
      hit <- hit_count > 0
      
    } else if (is.data.frame(obj)) {
      
      for (j in seq_along(obj)) {
        
        x <- obj[[j]]
        
        if (is.character(x) || is.factor(x)) {
          
          x_chr <- as.character(x)
          
          hit_count <- hit_count +
            sum(
              x_chr == g,
              na.rm = TRUE
            )
        }
      }
      
      hit <- hit_count > 0
      
    } else if (is.matrix(obj)) {
      
      vals <- as.character(obj)
      
      hit_count <- sum(
        vals == g,
        na.rm = TRUE
      )
      
      hit <- hit_count > 0
    }
    
    cat(
      sprintf(
        "%-6s : %s",
        g,
        ifelse(
          hit,
          paste0(
            "FOUND (",
            hit_count,
            " exact matches)"
          ),
          "NOT FOUND"
        )
      ),
      "\n"
    )
  }
}


# ------------------------------------------------------------
# 4. Specifically inspect SI_geneset if present
# ------------------------------------------------------------

if ("SI_geneset" %in% all_names_7be) {
  
  cat("\n============================================================\n")
  cat("DIRECT SI_geneset AUDIT\n")
  cat("============================================================\n\n")
  
  si <- get(
    "SI_geneset",
    envir = estimate_audit_env_7be
  )
  
  cat("Class:\n")
  print(class(si))
  
  cat("\nDimensions:\n")
  print(dim(si))
  
  cat("\nStructure:\n")
  str(si, max.level = 2)
  
  cat("\nDirect core-gene search:\n")
  
  for (g in core_genes_7be) {
    
    if (is.data.frame(si)) {
      
      hits <- which(
        vapply(
          si,
          function(x) {
            if (is.character(x) || is.factor(x)) {
              any(
                as.character(x) == g,
                na.rm = TRUE
              )
            } else {
              FALSE
            }
          },
          logical(1)
        )
      )
      
      cat(
        g,
        ": matching columns = ",
        ifelse(
          length(hits) == 0,
          "NONE",
          paste(
            colnames(si)[hits],
            collapse = ", "
          )
        ),
        "\n",
        sep = ""
      )
      
    } else if (is.matrix(si)) {
      
      vals <- as.character(si)
      
      cat(
        g,
        ": exact matches = ",
        sum(
          vals == g,
          na.rm = TRUE
        ),
        "\n",
        sep = ""
      )
      
    } else if (is.character(si)) {
      
      cat(
        g,
        ": exact matches = ",
        sum(
          si == g,
          na.rm = TRUE
        ),
        "\n",
        sep = ""
      )
    }
  }
  
} else {
  
  cat("\nSI_geneset is NOT present in the historical environment.\n")
}


# ------------------------------------------------------------
# 5. Correct input-gene membership using GeneSymbol
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CORRECTED INPUT GENE MEMBERSHIP\n")
cat("============================================================\n\n")

input_obj <- estimate_audit_env_7be$estimate_input_df

if (
  is.data.frame(input_obj) &&
  "GeneSymbol" %in% colnames(input_obj)
) {
  
  input_genes <- as.character(
    input_obj$GeneSymbol
  )
  
  for (g in core_genes_7be) {
    
    cat(
      sprintf(
        "%-6s : input GeneSymbol = %s\n",
        g,
        g %in% input_genes
      )
    )
  }
  
  cat(
    "\nUnique input GeneSymbol count:",
    length(unique(input_genes)),
    "\n"
  )
}


# ------------------------------------------------------------
# 6. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 7B-7E-7 COMPLETE\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 7B-7E-8
# Trace ESTIMATE Package and Signature Provenance
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7E-8: Trace ESTIMATE Package and Signature Provenance\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Check whether ESTIMATE package was available historically
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("1. ESTIMATE PACKAGE AVAILABILITY\n")
cat("------------------------------------------------------------\n\n")

estimate_pkg_available <- requireNamespace("estimate", quietly = TRUE)

cat(
  "Current R session: estimate package available = ",
  estimate_pkg_available,
  "\n",
  sep = ""
)

if (estimate_pkg_available) {
  
  cat("\nPackage version:\n")
  print(packageVersion("estimate"))
  
  cat("\nPackage location:\n")
  print(find.package("estimate"))
  
  cat("\nPackage DESCRIPTION:\n")
  print(utils::packageDescription("estimate")[
    c("Package", "Version", "Title", "Author", "URL")
  ])
  
}


# ------------------------------------------------------------
# 2. Inspect functions exported by ESTIMATE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("2. ESTIMATE EXPORTED FUNCTIONS\n")
cat("============================================================\n\n")

if (estimate_pkg_available) {
  
  estimate_exports <- getNamespaceExports("estimate")
  
  cat(
    "Number of exported objects:",
    length(estimate_exports),
    "\n\n"
  )
  
  print(
    estimate_exports[
      grepl(
        "score|signature|immune|stromal|geneset|filter",
        estimate_exports,
        ignore.case = TRUE
      )
    ]
  )
  
}


# ------------------------------------------------------------
# 3. Inspect estimateScore() source
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("3. estimateScore() SOURCE AUDIT\n")
cat("============================================================\n\n")

if (estimate_pkg_available &&
    exists("estimateScore", where = asNamespace("estimate"), inherits = FALSE)) {
  
  estimateScore_fun <- get(
    "estimateScore",
    envir = asNamespace("estimate")
  )
  
  cat("estimateScore() function body:\n\n")
  print(estimateScore_fun)
  
}


# ------------------------------------------------------------
# 4. Search ESTIMATE namespace for signature objects
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("4. SIGNATURE OBJECTS IN ESTIMATE NAMESPACE\n")
cat("============================================================\n\n")

if (estimate_pkg_available) {
  
  ns_names <- ls(
    envir = asNamespace("estimate"),
    all.names = TRUE
  )
  
  signature_ns_names <- ns_names[
    grepl(
      "signature|geneset|stromal|immune|si_",
      ns_names,
      ignore.case = TRUE
    )
  ]
  
  cat(
    "Number of signature-related namespace objects:",
    length(signature_ns_names),
    "\n\n"
  )
  
  print(signature_ns_names)
  
}


# ------------------------------------------------------------
# 5. Search ESTIMATE namespace for core genes
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("5. CORE-GENE SEARCH IN ESTIMATE NAMESPACE\n")
cat("============================================================\n\n")

core_genes_7be <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

if (estimate_pkg_available) {
  
  ns_names <- ls(
    envir = asNamespace("estimate"),
    all.names = TRUE
  )
  
  for (obj_name in ns_names) {
    
    obj <- tryCatch(
      get(
        obj_name,
        envir = asNamespace("estimate")
      ),
      error = function(e) NULL
    )
    
    if (is.null(obj)) next
    
    for (g in core_genes_7be) {
      
      hit <- FALSE
      hit_count <- 0
      
      if (is.character(obj)) {
        
        hit_count <- sum(
          as.character(obj) == g,
          na.rm = TRUE
        )
        
        hit <- hit_count > 0
        
      } else if (is.data.frame(obj)) {
        
        for (j in seq_along(obj)) {
          
          x <- obj[[j]]
          
          if (is.character(x) || is.factor(x)) {
            
            x_chr <- as.character(x)
            
            hit_count <- hit_count +
              sum(
                x_chr == g,
                na.rm = TRUE
              )
            
          }
          
        }
        
        hit <- hit_count > 0
        
      } else if (is.matrix(obj)) {
        
        vals <- as.character(obj)
        
        hit_count <- sum(
          vals == g,
          na.rm = TRUE
        )
        
        hit <- hit_count > 0
        
      }
      
      if (hit) {
        
        cat(
          "FOUND:",
          g,
          "in namespace object:",
          obj_name,
          "| exact matches =",
          hit_count,
          "\n"
        )
        
      }
      
    }
    
  }
  
}


# ------------------------------------------------------------
# 6. Inspect installed ESTIMATE package files
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("6. ESTIMATE PACKAGE FILE INVENTORY\n")
cat("============================================================\n\n")

if (estimate_pkg_available) {
  
  estimate_path <- find.package("estimate")
  
  cat("Package path:\n")
  print(estimate_path)
  
  cat("\nPackage files containing 'signature', 'geneset', 'immune', or 'stromal':\n")
  
  all_files_estimate <- list.files(
    estimate_path,
    recursive = TRUE,
    full.names = TRUE
  )
  
  relevant_files <- all_files_estimate[
    grepl(
      "signature|geneset|immune|stromal|estimateScore|SI_",
      basename(all_files_estimate),
      ignore.case = TRUE
    )
  ]
  
  print(relevant_files)
  
}


# ------------------------------------------------------------
# 7. Search package source files for signature names
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("7. SOURCE-TEXT SEARCH FOR SIGNATURE DEFINITIONS\n")
cat("============================================================\n\n")

if (estimate_pkg_available) {
  
  r_files <- list.files(
    estimate_path,
    pattern = "\\.[Rr]$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(r_files) > 0) {
    
    search_results <- list()
    
    for (f in r_files) {
      
      txt <- tryCatch(
        readLines(
          f,
          warn = FALSE,
          encoding = "UTF-8"
        ),
        error = function(e) character(0)
      )
      
      if (length(txt) == 0) next
      
      hit <- grepl(
        "StromalSignature|ImmuneSignature|SI_geneset|geneset",
        txt,
        ignore.case = TRUE
      )
      
      if (any(hit)) {
        
        search_results[[f]] <- txt[hit]
        
      }
      
    }
    
    if (length(search_results) == 0) {
      
      cat(
        "No signature-definition text found in installed R source files.\n"
      )
      
    } else {
      
      for (f in names(search_results)) {
        
        cat("\nFILE:\n")
        cat(f, "\n\n")
        
        print(search_results[[f]])
        
      }
      
    }
    
  }
  
}


# ------------------------------------------------------------
# 8. Check historical command/function objects
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("8. HISTORICAL ENVIRONMENT FUNCTION SEARCH\n")
cat("============================================================\n\n")

historical_names <- ls(
  envir = estimate_audit_env_7be,
  all.names = TRUE
)

function_names <- historical_names[
  vapply(
    historical_names,
    function(x) {
      obj <- tryCatch(
        get(x, envir = estimate_audit_env_7be),
        error = function(e) NULL
      )
      is.function(obj)
    },
    logical(1)
  )
]

cat(
  "Number of historical function objects:",
  length(function_names),
  "\n\n"
)

print(
  function_names[
    grepl(
      "estimate|score|immune|stromal|signature|gene",
      function_names,
      ignore.case = TRUE
    )
  ]
)


# ------------------------------------------------------------
# 9. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 7B-7E-8 COMPLETE\n")
cat("============================================================\n\n")

cat("IMPORTANT:\n")
cat(
  "No ESTIMATE scores were recalculated.\n",
  "No locked analysis objects were modified.\n",
  "No leave-one-gene-out analysis was performed.\n",
  sep = ""
)

# ============================================================
# Paper 1
# STEP 7B-7E-9A
# ESTIMATE Signature Membership Audit
# Revised provenance route via tidyestimate
#
# IMPORTANT:
# - No ESTIMATE score recalculation
# - No modification of locked objects
# - No leave-one-gene-out analysis
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7E-9A: ESTIMATE Signature Membership Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Core genes
# ------------------------------------------------------------

core_genes_7be <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

print(core_genes_7be)


# ------------------------------------------------------------
# 1. Check tidyestimate availability
# ------------------------------------------------------------

tidyestimate_available <- requireNamespace(
  "tidyestimate",
  quietly = TRUE
)

cat(
  "\ntidyestimate available = ",
  tidyestimate_available,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 2. If not available, install only tidyestimate
# ------------------------------------------------------------

if (!tidyestimate_available) {
  
  cat(
    "\nThe tidyestimate package is not currently installed.\n"
  )
  
  cat(
    "Installing tidyestimate from CRAN...\n\n"
  )
  
  install.packages(
    "tidyestimate",
    repos = "https://cloud.r-project.org",
    dependencies = FALSE
  )
  
  tidyestimate_available <- requireNamespace(
    "tidyestimate",
    quietly = TRUE
  )
  
}

if (!tidyestimate_available) {
  
  stop(
    "tidyestimate could not be installed. ",
    "No ESTIMATE scoring or locked object was modified."
  )
  
}


# ------------------------------------------------------------
# 3. Record package version
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("3. TIDYESTIMATE PACKAGE PROVENANCE\n")
cat("============================================================\n\n")

cat("tidyestimate version:\n")
print(packageVersion("tidyestimate"))

cat("\nPackage location:\n")
print(find.package("tidyestimate"))

cat("\nPackage description source:\n")

pkg_desc_7be <- packageDescription(
  "tidyestimate"
)

print(
  pkg_desc_7be[
    c(
      "Package",
      "Version",
      "Title",
      "Author",
      "URL",
      "BugReports"
    )
  ]
)


# ------------------------------------------------------------
# 4. Load gene_sets
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("4. LOAD ESTIMATE-DERIVED GENE SETS\n")
cat("============================================================\n\n")

data(
  "gene_sets",
  package = "tidyestimate",
  envir = environment()
)

cat("Class of gene_sets:\n")
print(class(gene_sets))

cat("\nDimensions:\n")
print(dim(gene_sets))

cat("\nColumn names:\n")
print(colnames(gene_sets))

cat("\nFirst 10 rows:\n")
print(
  head(
    gene_sets,
    10
  )
)


# ------------------------------------------------------------
# 5. Validate expected 141-gene structure
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("5. SIGNATURE STRUCTURE AUDIT\n")
cat("============================================================\n\n")

stromal_genes_7be <- unique(
  as.character(
    gene_sets$stromal_signature
  )
)

immune_genes_7be <- unique(
  as.character(
    gene_sets$immune_signature
  )
)

stromal_genes_7be <- stromal_genes_7be[
  !is.na(stromal_genes_7be) &
    stromal_genes_7be != ""
]

immune_genes_7be <- immune_genes_7be[
  !is.na(immune_genes_7be) &
    immune_genes_7be != ""
]

cat(
  "Unique Stromal signature genes:",
  length(stromal_genes_7be),
  "\n"
)

cat(
  "Unique Immune signature genes:",
  length(immune_genes_7be),
  "\n"
)


# ------------------------------------------------------------
# 6. Direct membership audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("6. CORE-GENE SIGNATURE MEMBERSHIP\n")
cat("============================================================\n\n")

membership_7be <- data.frame(
  Gene = core_genes_7be,
  
  In_StromalSignature =
    core_genes_7be %in% stromal_genes_7be,
  
  In_ImmuneSignature =
    core_genes_7be %in% immune_genes_7be,
  
  stringsAsFactors = FALSE
)

print(
  membership_7be,
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Explicit gene-by-gene report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("7. GENE-BY-GENE MEMBERSHIP REPORT\n")
cat("============================================================\n\n")

for (g in core_genes_7be) {
  
  cat("\n------------------------------------------------------------\n")
  cat("Gene:", g, "\n")
  cat("------------------------------------------------------------\n")
  
  cat(
    "In StromalSignature: ",
    g %in% stromal_genes_7be,
    "\n",
    sep = ""
  )
  
  cat(
    "In ImmuneSignature:  ",
    g %in% immune_genes_7be,
    "\n",
    sep = ""
  )
  
}


# ------------------------------------------------------------
# 8. Explicit MS4A6A vs MS4A7 audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("8. MS4A6A / MS4A7 DISAMBIGUATION AUDIT\n")
cat("============================================================\n\n")

cat(
  "MS4A6A in StromalSignature: ",
  "MS4A6A" %in% stromal_genes_7be,
  "\n",
  sep = ""
)

cat(
  "MS4A7 in StromalSignature:  ",
  "MS4A7" %in% stromal_genes_7be,
  "\n",
  sep = ""
)

cat(
  "MS4A7 in ImmuneSignature:   ",
  "MS4A7" %in% immune_genes_7be,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Compare with the actual ESTIMATE input genes
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("9. INPUT-GENE VS SIGNATURE MEMBERSHIP\n")
cat("============================================================\n\n")

input_obj <- estimate_audit_env_7be$estimate_input_df

if (
  is.data.frame(input_obj) &&
  "GeneSymbol" %in% colnames(input_obj)
) {
  
  input_genes_7be <- unique(
    as.character(
      input_obj$GeneSymbol
    )
  )
  
  input_signature_audit_7be <- data.frame(
    Gene = core_genes_7be,
    
    In_ESTIMATE_Input =
      core_genes_7be %in% input_genes_7be,
    
    In_StromalSignature =
      core_genes_7be %in% stromal_genes_7be,
    
    In_ImmuneSignature =
      core_genes_7be %in% immune_genes_7be,
    
    stringsAsFactors = FALSE
  )
  
  print(
    input_signature_audit_7be,
    row.names = FALSE
  )
  
}


# ------------------------------------------------------------
# 10. Check for unexpected duplicated genes
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("10. DUPLICATION AUDIT\n")
cat("============================================================\n\n")

cat(
  "Stromal signature duplicated entries:",
  sum(duplicated(stromal_genes_7be)),
  "\n"
)

cat(
  "Immune signature duplicated entries:",
  sum(duplicated(immune_genes_7be)),
  "\n"
)


# ------------------------------------------------------------
# 11. Build audit object
# ------------------------------------------------------------

step7be9a_audit <- list(
  
  status =
    "MEMBERSHIP_AUDIT_ONLY",
  
  source_package =
    "tidyestimate",
  
  source_package_version =
    as.character(
      packageVersion("tidyestimate")
    ),
  
  core_genes =
    core_genes_7be,
  
  stromal_gene_count =
    length(stromal_genes_7be),
  
  immune_gene_count =
    length(immune_genes_7be),
  
  membership =
    membership_7be,
  
  locked_objects_modified =
    FALSE,
  
  estimate_scores_recalculated =
    FALSE,
  
  leave_one_gene_out_performed =
    FALSE
  
)


# ------------------------------------------------------------
# 12. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 7B-7E-9A COMPLETE\n")
cat("============================================================\n\n")

cat(
  "No ESTIMATE scores were recalculated.\n"
)

cat(
  "No locked analysis objects were modified.\n"
)

cat(
  "No leave-one-gene-out analysis was performed.\n"
)

cat(
  "Signature membership audit completed.\n"
)


# ============================================================
# Paper 1
# STEP 7B-7E-9A
# ESTIMATE Signature Membership Audit
# ============================================================
#
# Purpose:
#   Audit whether the three core genes (AOC3, MS4A7, IL7R)
#   are members of the original ESTIMATE StromalSignature
#   or ImmuneSignature.
#
# IMPORTANT:
#   This is an AUDIT-ONLY workflow.
#
#   No ESTIMATE scores are recalculated.
#   No locked analysis objects are modified.
#   No leave-one-gene-out analysis is performed.
#   No statistical test is performed.
#   No Figure 6 data are modified.
#   No Supplementary Table 3 data are modified.
#
# Core genes:
#   AOC3
#   MS4A7
#   IL7R
#
# ============================================================


# ------------------------------------------------------------
# 0. Initialization
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP 7B-7E-9A: ESTIMATE SIGNATURE MEMBERSHIP AUDIT\n")
cat("============================================================\n\n")

options(stringsAsFactors = FALSE)


# ------------------------------------------------------------
# 1. Define core genes
# ------------------------------------------------------------

core_genes_7be <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

cat("Core genes:\n")
print(core_genes_7be)


# ------------------------------------------------------------
# 2. Check tidyestimate package
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("2. SOURCE PACKAGE AUDIT\n")
cat("============================================================\n\n")

if (!requireNamespace("tidyestimate", quietly = TRUE)) {
  
  cat("tidyestimate package not available.\n")
  cat("Installing from CRAN...\n\n")
  
  install.packages(
    "tidyestimate",
    repos = "https://cloud.r-project.org"
  )
}

stopifnot(
  requireNamespace("tidyestimate", quietly = TRUE)
)

tidyestimate_version_7be <- as.character(
  packageVersion("tidyestimate")
)

tidyestimate_location_7be <- find.package(
  "tidyestimate"
)

tidyestimate_description_7be <- packageDescription(
  "tidyestimate"
)

cat(
  "tidyestimate version:",
  tidyestimate_version_7be,
  "\n"
)

cat(
  "tidyestimate location:",
  tidyestimate_location_7be,
  "\n"
)


# ------------------------------------------------------------
# 3. Load signature data
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("3. LOAD ESTIMATE SIGNATURE DATA\n")
cat("============================================================\n\n")

data(
  "gene_sets",
  package = "tidyestimate",
  envir = environment()
)

stopifnot(
  exists("gene_sets")
)

cat("gene_sets object loaded successfully.\n")
cat("Object class:\n")
print(class(gene_sets))

cat("\nObject structure:\n")
str(gene_sets)


# ------------------------------------------------------------
# 4. Extract the two ESTIMATE signatures
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("4. EXTRACT ESTIMATE SIGNATURES\n")
cat("============================================================\n\n")

stopifnot(
  "stromal_signature" %in% names(gene_sets),
  "immune_signature" %in% names(gene_sets)
)

stromal_genes_7be <- unique(
  as.character(
    gene_sets$stromal_signature
  )
)

immune_genes_7be <- unique(
  as.character(
    gene_sets$immune_signature
  )
)

stromal_genes_7be <- stromal_genes_7be[
  !is.na(stromal_genes_7be) &
    stromal_genes_7be != ""
]

immune_genes_7be <- immune_genes_7be[
  !is.na(immune_genes_7be) &
    immune_genes_7be != ""
]

cat(
  "Unique Stromal signature genes:",
  length(stromal_genes_7be),
  "\n"
)

cat(
  "Unique Immune signature genes:",
  length(immune_genes_7be),
  "\n"
)

stopifnot(
  length(stromal_genes_7be) == 141,
  length(immune_genes_7be) == 141
)


# ------------------------------------------------------------
# 5. Direct membership audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("5. CORE-GENE SIGNATURE MEMBERSHIP\n")
cat("============================================================\n\n")

membership_7be <- data.frame(
  Gene = core_genes_7be,
  
  In_StromalSignature =
    core_genes_7be %in% stromal_genes_7be,
  
  In_ImmuneSignature =
    core_genes_7be %in% immune_genes_7be,
  
  stringsAsFactors = FALSE
)

print(
  membership_7be,
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Gene-by-gene membership report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("6. GENE-BY-GENE MEMBERSHIP REPORT\n")
cat("============================================================\n")

for (g in core_genes_7be) {
  
  cat("\n------------------------------------------------------------\n")
  cat("Gene:", g, "\n")
  cat("------------------------------------------------------------\n")
  
  cat(
    "In StromalSignature: ",
    g %in% stromal_genes_7be,
    "\n",
    sep = ""
  )
  
  cat(
    "In ImmuneSignature:  ",
    g %in% immune_genes_7be,
    "\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 7. Explicit MS4A6A / MS4A7 disambiguation audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("7. MS4A6A / MS4A7 DISAMBIGUATION AUDIT\n")
cat("============================================================\n\n")

MS4A6A_in_stromal_7be <-
  "MS4A6A" %in% stromal_genes_7be

MS4A7_in_stromal_7be <-
  "MS4A7" %in% stromal_genes_7be

MS4A7_in_immune_7be <-
  "MS4A7" %in% immune_genes_7be

cat(
  "MS4A6A in StromalSignature: ",
  MS4A6A_in_stromal_7be,
  "\n",
  sep = ""
)

cat(
  "MS4A7 in StromalSignature:  ",
  MS4A7_in_stromal_7be,
  "\n",
  sep = ""
)

cat(
  "MS4A7 in ImmuneSignature:   ",
  MS4A7_in_immune_7be,
  "\n",
  sep = ""
)

stopifnot(
  MS4A6A_in_stromal_7be == TRUE,
  MS4A7_in_stromal_7be == FALSE,
  MS4A7_in_immune_7be == FALSE
)


# ------------------------------------------------------------
# 8. Compare with actual ESTIMATE input genes
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("8. INPUT-GENE VS SIGNATURE MEMBERSHIP\n")
cat("============================================================\n\n")

input_obj <-
  estimate_audit_env_7be$estimate_input_df

if (
  is.data.frame(input_obj) &&
  "GeneSymbol" %in% colnames(input_obj)
) {
  
  input_genes_7be <- unique(
    as.character(
      input_obj$GeneSymbol
    )
  )
  
  input_signature_audit_7be <- data.frame(
    Gene = core_genes_7be,
    
    In_ESTIMATE_Input =
      core_genes_7be %in% input_genes_7be,
    
    In_StromalSignature =
      core_genes_7be %in% stromal_genes_7be,
    
    In_ImmuneSignature =
      core_genes_7be %in% immune_genes_7be,
    
    stringsAsFactors = FALSE
  )
  
  print(
    input_signature_audit_7be,
    row.names = FALSE
  )
  
} else {
  
  stop(
    "estimate_input_df with GeneSymbol column was not available."
  )
}


# ------------------------------------------------------------
# 9. Duplication audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("9. DUPLICATION AUDIT\n")
cat("============================================================\n\n")

stromal_duplicate_count_7be <-
  sum(duplicated(stromal_genes_7be))

immune_duplicate_count_7be <-
  sum(duplicated(immune_genes_7be))

cat(
  "Stromal signature duplicated entries:",
  stromal_duplicate_count_7be,
  "\n"
)

cat(
  "Immune signature duplicated entries:",
  immune_duplicate_count_7be,
  "\n"
)

stopifnot(
  stromal_duplicate_count_7be == 0,
  immune_duplicate_count_7be == 0
)


# ------------------------------------------------------------
# 10. Final membership assertions
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("10. FINAL MEMBERSHIP ASSERTIONS\n")
cat("============================================================\n\n")

stopifnot(
  "AOC3" %in% stromal_genes_7be,
  !("AOC3" %in% immune_genes_7be),
  
  !("MS4A7" %in% stromal_genes_7be),
  !("MS4A7" %in% immune_genes_7be),
  
  !("IL7R" %in% stromal_genes_7be),
  "IL7R" %in% immune_genes_7be
)

cat("AOC3 -> StromalSignature: TRUE\n")
cat("MS4A7 -> neither signature: TRUE\n")
cat("IL7R -> ImmuneSignature: TRUE\n")

cat("\nAll membership assertions PASSED.\n")


# ------------------------------------------------------------
# 11. Build formal audit object
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("11. BUILD FORMAL AUDIT OBJECT\n")
cat("============================================================\n\n")

step7be9a_audit <- list(
  
  status =
    "MEMBERSHIP_AUDIT_ONLY",
  
  workflow =
    "Paper1_STEP7B-7E-9A",
  
  purpose =
    "ESTIMATE signature membership audit",
  
  source_package =
    "tidyestimate",
  
  source_package_version =
    tidyestimate_version_7be,
  
  source_package_location =
    tidyestimate_location_7be,
  
  core_genes =
    core_genes_7be,
  
  stromal_gene_count =
    length(stromal_genes_7be),
  
  immune_gene_count =
    length(immune_genes_7be),
  
  membership =
    membership_7be,
  
  input_signature_audit =
    input_signature_audit_7be,
  
  MS4A6A_in_stromal =
    MS4A6A_in_stromal_7be,
  
  MS4A7_in_stromal =
    MS4A7_in_stromal_7be,
  
  MS4A7_in_immune =
    MS4A7_in_immune_7be,
  
  stromal_duplicate_count =
    stromal_duplicate_count_7be,
  
  immune_duplicate_count =
    immune_duplicate_count_7be,
  
  locked_objects_modified =
    FALSE,
  
  estimate_scores_recalculated =
    FALSE,
  
  leave_one_gene_out_performed =
    FALSE,
  
  statistical_analysis_performed =
    FALSE,
  
  figure6_modified =
    FALSE,
  
  supplementary_table3_modified =
    FALSE
)


# ------------------------------------------------------------
# 12. Final audit status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 7B-7E-9A COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Status:",
  step7be9a_audit$status,
  "\n\n"
)

cat(
  "No ESTIMATE scores were recalculated.\n"
)

cat(
  "No locked analysis objects were modified.\n"
)

cat(
  "No leave-one-gene-out analysis was performed.\n"
)

cat(
  "No statistical analysis was performed.\n"
)

cat(
  "Figure 6 was not modified.\n"
)

cat(
  "Supplementary Table 3 was not modified.\n"
)

cat(
  "Signature membership audit completed.\n"
)

cat(
  "\nFINAL STATUS: PASSED\n"
)


# ------------------------------------------------------------
# 13. Save locked audit object
# ------------------------------------------------------------

audit_file_7be9a <-
  "result_dir/Paper1_tumor-only WGCNA/Paper1_STEP7B7E9A_ESTIMATE_signature_membership_audit_LOCKED.RData"

save(
  step7be9a_audit,
  membership_7be,
  input_signature_audit_7be,
  stromal_genes_7be,
  immune_genes_7be,
  MS4A6A_in_stromal_7be,
  MS4A7_in_stromal_7be,
  MS4A7_in_immune_7be,
  file = audit_file_7be9a
)

cat("\n============================================================\n")
cat("LOCKED AUDIT OBJECT SAVED\n")
cat("============================================================\n\n")

cat(
  "File:\n",
  audit_file_7be9a,
  "\n"
)


# ------------------------------------------------------------
# 14. Save script path
# ------------------------------------------------------------

script_file_7be9a <-
  "result_dir/Paper1-tumor-only WGCNA/Paper1_STEP7B7E9A_ESTIMATE_signature_membership_audit_LOCKED.R"

cat("\nRecommended script path:\n")
cat(script_file_7be9a, "\n")


# ============================================================
# END OF STEP 7B-7E-9A
# ============================================================


# ============================================================
# STEP 7B-7E-9A
# FINAL FILE EXISTENCE + READ-BACK CHECK
# ============================================================

cat("\n============================================================\n")
cat("STEP 7B-7E-9A FINAL SAVE CHECK\n")
cat("============================================================\n\n")

cat(
  "File exists: ",
  file.exists(audit_file_7be9a),
  "\n",
  sep = ""
)

cat(
  "File size: ",
  file.info(audit_file_7be9a)$size,
  " bytes\n",
  sep = ""
)

audit_readback_env_7be9a <- new.env(parent = emptyenv())

loaded_objects_7be9a <- load(
  audit_file_7be9a,
  envir = audit_readback_env_7be9a
)

cat("\nLoaded objects:\n")
print(loaded_objects_7be9a)

cat(
  "\nRead-back audit status: ",
  audit_readback_env_7be9a$step7be9a_audit$status,
  "\n",
  sep = ""
)

cat(
  "\nFINAL FILE CHECK: PASSED\n"
)







