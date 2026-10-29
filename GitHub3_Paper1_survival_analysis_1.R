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

work_dir <- input_dir

setwd(work_dir)

getwd()

load(
  "input_dir/Workflow-04-2D_final_WGCNA_candidate_checkpoint.RData"
)

ls()

# ============================================================
# STEP 6B-0: Audit currently loaded WGCNA objects
# ============================================================

cat("\n============================================================\n")
cat("STEP 6B-0: Current WGCNA workspace audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. List objects
# ------------------------------------------------------------

cat("Objects currently loaded:\n")
print(ls())


# ------------------------------------------------------------
# 2. Inspect clinical_aligned
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("clinical_aligned\n")
cat("------------------------------------------------------------\n")

print(class(clinical_aligned))
print(dim(clinical_aligned))
print(names(clinical_aligned))

cat("\nFirst 3 rows:\n")
print(utils::head(clinical_aligned, 3))


# ------------------------------------------------------------
# 3. Inspect datExpr_MAD
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("datExpr_MAD\n")
cat("------------------------------------------------------------\n")

print(class(datExpr_MAD))
print(dim(datExpr_MAD))

cat("\nFirst 5 row names:\n")
print(head(rownames(datExpr_MAD), 5))

cat("\nFirst 5 column names:\n")
print(head(colnames(datExpr_MAD), 5))


# ------------------------------------------------------------
# 4. Check the three core genes
# ------------------------------------------------------------

core_genes <- c("AOC3", "MS4A7", "IL7R")

cat("\n------------------------------------------------------------\n")
cat("Core gene presence in datExpr_MAD\n")
cat("------------------------------------------------------------\n")

print(core_genes %in% colnames(datExpr_MAD))
print(core_genes %in% rownames(datExpr_MAD))


# ------------------------------------------------------------
# 5. Inspect gene_MAD
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("gene_MAD\n")
cat("------------------------------------------------------------\n")

print(class(gene_MAD))
print(length(gene_MAD))

cat("\nCore genes in gene_MAD:\n")
print(core_genes %in% names(gene_MAD))

cat("\n============================================================\n")
cat("STEP 6B-0 audit completed\n")
cat("============================================================\n")

# ============================================================
# STEP 6B-0B
# Audit Ensembl IDs for AOC3 / MS4A7 / IL7R
# ============================================================

cat("\n============================================================\n")
cat("STEP 6B-0B: Core gene Ensembl-ID audit\n")
cat("============================================================\n\n")

core_gene_map <- data.frame(
  Gene = c("AOC3", "MS4A7", "IL7R"),
  Ensembl_ID = c(
    "ENSG00000131471",
    "ENSG00000166927",
    "ENSG00000168685"
  ),
  stringsAsFactors = FALSE
)

cat("Core gene mapping:\n")
print(core_gene_map)

cat("\nPresence in datExpr_MAD:\n")

core_gene_map$Present_in_datExpr_MAD <-
  core_gene_map$Ensembl_ID %in% colnames(datExpr_MAD)

print(core_gene_map)

cat("\nExpression summary for available core genes:\n")

for (i in seq_len(nrow(core_gene_map))) {
  
  g <- core_gene_map$Gene[i]
  id <- core_gene_map$Ensembl_ID[i]
  
  if (id %in% colnames(datExpr_MAD)) {
    
    x <- datExpr_MAD[, id]
    
    cat("\n", g, " | ", id, "\n", sep = "")
    print(summary(x))
    cat("SD =", sd(x, na.rm = TRUE), "\n")
    cat("NA =", sum(is.na(x)), "\n")
    cat("NaN =", sum(is.nan(x)), "\n")
    cat("Inf =", sum(is.infinite(x)), "\n")
    
  } else {
    
    cat("\n", g, " | ", id, " NOT FOUND\n", sep = "")
  }
}

cat("\n============================================================\n")
cat("STEP 6B-0B completed\n")
cat("============================================================\n")

# ============================================================
# STEP 6B-1A
# Audit survival-related files from the WGCNA project folder
# ============================================================

wgcna_dir <- "input_dir/Paper1_tumor-only WGCNA"

cat("\n============================================================\n")
cat("STEP 6B-1A: WGCNA survival-file audit\n")
cat("============================================================\n\n")

cat("WGCNA directory:\n")
print(wgcna_dir)

# ------------------------------------------------------------
# 1. List files containing OS / survival / Cox / clinical terms
# ------------------------------------------------------------

all_files <- list.files(
  path = wgcna_dir,
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE
)

survival_files <- all_files[
  grepl(
    "OS|survival|cox|clinical|KM|Kaplan|progn",
    basename(all_files),
    ignore.case = TRUE
  )
]

cat("\n------------------------------------------------------------\n")
cat("Potential survival-related files:\n")
cat("------------------------------------------------------------\n")

print(basename(survival_files))


# ------------------------------------------------------------
# 2. Specifically inspect Workflow-04-2K
# ------------------------------------------------------------

target_file <- file.path(
  wgcna_dir,
  "Workflow-04-2K_WGCNA_OS_Cox_association.txt"
)

cat("\n------------------------------------------------------------\n")
cat("Target file:\n")
cat("------------------------------------------------------------\n")

print(target_file)

cat("\nFile exists:\n")
print(file.exists(target_file))

if (file.exists(target_file)) {
  
  cat("\nFile information:\n")
  print(file.info(target_file))
  
  cat("\nFirst 10 lines:\n")
  print(
    readLines(
      target_file,
      n = 10,
      warn = FALSE
    )
  )
}


# ------------------------------------------------------------
# 3. Show all RData / RDS files in WGCNA directory
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("RData / RDS files:\n")
cat("------------------------------------------------------------\n")

print(
  basename(
    all_files[
      grepl(
        "\\.(RData|rdata|RDA|rds)$",
        all_files,
        ignore.case = TRUE
      )
    ]
  )
)


cat("\n============================================================\n")
cat("STEP 6B-1A audit completed\n")
cat("============================================================\n")

# ============================================================
# STEP 6B-1B
# Inspect the locked WGCNA primary survival cohort QC
# ============================================================

wgcna_dir <- "input_dir/Paper1_tumor-only WGCNA"

primary_qc_file <- file.path(
  wgcna_dir,
  "Workflow-04-2N_WGCNA_OS_primary_cohort_QC.txt"
)

cat("\n============================================================\n")
cat("STEP 6B-1B: Primary survival cohort QC audit\n")
cat("============================================================\n\n")

cat("File:\n")
print(primary_qc_file)

cat("\nFile exists:\n")
print(file.exists(primary_qc_file))

if (file.exists(primary_qc_file)) {
  
  cat("\n------------------------------------------------------------\n")
  cat("Full primary cohort QC file\n")
  cat("------------------------------------------------------------\n\n")
  
  qc_lines <- readLines(
    primary_qc_file,
    warn = FALSE
  )
  
  cat(paste(qc_lines, collapse = "\n"))
  
}

cat("\n\n============================================================\n")
cat("STEP 6B-1B completed\n")
cat("============================================================\n")

# ============================================================
# STEP 6B-1C
# Inspect the WGCNA sensitivity survival cohort QC
# ============================================================

sensitivity_qc_file <- file.path(
  wgcna_dir,
  "Workflow-04-2N_WGCNA_OS_sensitivity_cohort_QC.txt"
)

cat("\n============================================================\n")
cat("STEP 6B-1C: Sensitivity survival cohort QC audit\n")
cat("============================================================\n\n")

cat("File:\n")
print(sensitivity_qc_file)

cat("\nFile exists:\n")
print(file.exists(sensitivity_qc_file))

if (file.exists(sensitivity_qc_file)) {
  
  cat("\n------------------------------------------------------------\n")
  cat("Full sensitivity cohort QC file\n")
  cat("------------------------------------------------------------\n\n")
  
  qc_lines_sens <- readLines(
    sensitivity_qc_file,
    warn = FALSE
  )
  
  cat(paste(qc_lines_sens, collapse = "\n"))
}

cat("\n\n============================================================\n")
cat("STEP 6B-1C completed\n")
cat("============================================================\n")


# ============================================================
# STEP 6B-2
# Construct locked 485/178 primary survival cohort
# with AOC3 / MS4A7 / IL7R expression
# ============================================================

cat("\n============================================================\n")
cat("STEP 6B-2: Construct primary survival expression dataset\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define core genes
# ------------------------------------------------------------

core_gene_ids <- c(
  AOC3  = "ENSG00000131471",
  MS4A7 = "ENSG00000166927",
  IL7R  = "ENSG00000168685"
)

cat("Core genes and Ensembl IDs:\n")
print(core_gene_ids)


# ------------------------------------------------------------
# 2. Verify expression matrix contains all three genes
# ------------------------------------------------------------

gene_presence <- core_gene_ids %in% colnames(datExpr_MAD)

cat("\nGene presence in datExpr_MAD:\n")
print(gene_presence)

stopifnot(all(gene_presence))


# ------------------------------------------------------------
# 3. Verify clinical and expression sample IDs
# ------------------------------------------------------------

clinical_samples <- clinical_aligned$Sample
expr_samples <- rownames(datExpr_MAD)

cat("\nNumber of clinical samples:\n")
print(length(clinical_samples))

cat("Number of expression samples:\n")
print(length(expr_samples))

cat("\nClinical samples duplicated:\n")
print(anyDuplicated(clinical_samples))

cat("Expression samples duplicated:\n")
print(anyDuplicated(expr_samples))


# ------------------------------------------------------------
# 4. Identify survival-eligible observations
#    using the locked primary cohort definition
# ------------------------------------------------------------

survival_complete <- complete.cases(
  clinical_aligned[, c(
    "OS.time",
    "OS.status",
    "Age",
    "Gender",
    "pathologic_stage"
  )]
)

cat("\nSurvival complete-case count based on locked primary variables:\n")
print(sum(survival_complete))


# ------------------------------------------------------------
# 5. Construct primary cohort
# ------------------------------------------------------------

primary_clinical <- clinical_aligned[
  survival_complete,
  ,
  drop = FALSE
]

cat("\nPrimary clinical cohort dimensions:\n")
print(dim(primary_clinical))


# ------------------------------------------------------------
# 6. Match expression samples to primary clinical samples
# ------------------------------------------------------------

expr_idx <- match(
  primary_clinical$Sample,
  rownames(datExpr_MAD)
)

cat("\nExpression samples successfully matched:\n")
print(sum(!is.na(expr_idx)))

cat("Expression samples NOT matched:\n")
print(sum(is.na(expr_idx)))

stopifnot(all(!is.na(expr_idx)))


# ------------------------------------------------------------
# 7. Extract expression of the three core genes
# ------------------------------------------------------------

expr_core <- datExpr_MAD[
  expr_idx,
  unname(core_gene_ids),
  drop = FALSE
]

colnames(expr_core) <- names(core_gene_ids)


# ------------------------------------------------------------
# 8. Build locked Step 6B primary dataset
# ------------------------------------------------------------

survival_gene_primary <- cbind(
  primary_clinical[
    ,
    c(
      "Sample",
      "patient_id",
      "Age",
      "Gender",
      "pathologic_stage",
      "pathologic_t",
      "pathologic_n",
      "pathologic_m",
      "OS.time",
      "OS.status"
    ),
    drop = FALSE
  ],
  as.data.frame(expr_core)
)


# ------------------------------------------------------------
# 9. Basic QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("PRIMARY SURVIVAL DATASET QC\n")
cat("============================================================\n\n")

cat("Dimensions:\n")
print(dim(survival_gene_primary))

cat("\nEvents:\n")
print(table(survival_gene_primary$OS.status))

cat("\nMissing values:\n")
print(colSums(is.na(survival_gene_primary)))

cat("\nExpression summary:\n")
print(
  summary(
    survival_gene_primary[
      ,
      c("AOC3", "MS4A7", "IL7R")
    ]
  )
)


# ------------------------------------------------------------
# 10. Confirm expected locked cohort
# ------------------------------------------------------------

stopifnot(
  nrow(survival_gene_primary) == 485,
  sum(survival_gene_primary$OS.status == 1) == 178,
  sum(survival_gene_primary$OS.status == 0) == 307
)

cat("\n============================================================\n")
cat("STEP 6B-2 PASSED\n")
cat("Locked primary cohort confirmed: N = 485, Events = 178\n")
cat("============================================================\n")

# ============================================================
# STEP 6B-2A
# Inspect saved 2N cohort QC object structure
# WITHOUT loading it into the active workspace
# ============================================================

cat("\n============================================================\n")
cat("STEP 6B-2A: Inspect 2N cohort QC object\n")
cat("============================================================\n\n")

qc_rds_file <- file.path(
  wgcna_dir,
  "Workflow-04-2N_WGCNA_OS_multivariable_cohort_QC.rds"
)

cat("QC RDS file:\n")
print(qc_rds_file)

cat("\nFile exists:\n")
print(file.exists(qc_rds_file))


# ------------------------------------------------------------
# Read RDS into a temporary object only
# ------------------------------------------------------------

cohort_qc_obj <- readRDS(qc_rds_file)

cat("\nTop-level object class:\n")
print(class(cohort_qc_obj))

cat("\nTop-level object structure:\n")
str(cohort_qc_obj, max.level = 2)


# ------------------------------------------------------------
# If list-like, inspect names
# ------------------------------------------------------------

if (is.list(cohort_qc_obj)) {
  
  cat("\n============================================================\n")
  cat("Top-level names\n")
  cat("============================================================\n\n")
  
  print(names(cohort_qc_obj))
}


# ------------------------------------------------------------
# Inspect dimensions of major components
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Component dimensions / lengths\n")
cat("============================================================\n\n")

if (is.list(cohort_qc_obj)) {
  
  for (nm in names(cohort_qc_obj)) {
    
    obj_i <- cohort_qc_obj[[nm]]
    
    cat("\n---", nm, "---\n")
    cat("Class:", paste(class(obj_i), collapse = ", "), "\n")
    
    if (is.data.frame(obj_i) || is.matrix(obj_i)) {
      cat("Dimension:", paste(dim(obj_i), collapse = " x "), "\n")
      cat("Column names:\n")
      print(colnames(obj_i))
    } else {
      cat("Length:", length(obj_i), "\n")
    }
  }
}


cat("\n============================================================\n")
cat("STEP 6B-2A completed\n")
cat("============================================================\n")

# ============================================================
# STEP 6B-2B
# Extract and validate the locked 485/178 primary cohort
# directly from the saved 2N cohort object
# ============================================================

cat("\n============================================================\n")
cat("STEP 6B-2B: Validate locked primary cohort membership\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Inspect one primary cohort
# ------------------------------------------------------------

primary_brown <- cohort_qc_obj$primary_cohorts$ME_brown

cat("ME_brown primary cohort dimensions:\n")
print(dim(primary_brown))

cat("\nME_brown column names:\n")
print(colnames(primary_brown))

cat("\nFirst 6 rows:\n")
print(head(primary_brown))


# ------------------------------------------------------------
# 2. Identify likely sample-ID columns
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Column classes\n")
cat("============================================================\n\n")

print(
  data.frame(
    Column = colnames(primary_brown),
    Class = sapply(primary_brown, function(x) paste(class(x), collapse = ", ")),
    stringsAsFactors = FALSE
  )
)


# ------------------------------------------------------------
# 3. Inspect all character columns
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Character columns: first values\n")
cat("============================================================\n\n")

char_cols <- names(primary_brown)[
  sapply(primary_brown, is.character)
]

print(char_cols)

for (cc in char_cols) {
  cat("\n---", cc, "---\n")
  print(head(primary_brown[[cc]], 10))
}


# ------------------------------------------------------------
# 4. Check survival outcome
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Outcome variables\n")
cat("============================================================\n\n")

print(
  lapply(
    primary_brown[
      ,
      intersect(
        c("OS.time", "OS_time", "OS.status", "OS_event",
          "Sample", "patient_id"),
        colnames(primary_brown)
      ),
      drop = FALSE
    ],
    function(x) {
      if (is.factor(x)) x <- as.character(x)
      summary(x)
    }
  )
)


# ------------------------------------------------------------
# 5. Check whether all five module cohorts have identical rows
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Cross-module cohort consistency\n")
cat("============================================================\n\n")

primary_names <- names(cohort_qc_obj$primary_cohorts)

cat("Primary cohorts:\n")
print(primary_names)

# Use rownames first if present
cat("\nRownames present?\n")
print(
  sapply(
    cohort_qc_obj$primary_cohorts,
    function(x) !is.null(rownames(x))
  )
)


# ------------------------------------------------------------
# 6. Look for a Sample column
# ------------------------------------------------------------

sample_col_candidates <- intersect(
  c(
    "Sample",
    "sample",
    "SampleID",
    "sample_id",
    "patient_id"
  ),
  colnames(primary_brown)
)

cat("\nCandidate ID columns:\n")
print(sample_col_candidates)


# ------------------------------------------------------------
# 7. If Sample exists, verify all five cohorts
# ------------------------------------------------------------

if ("Sample" %in% colnames(primary_brown)) {
  
  primary_sample_sets <- lapply(
    cohort_qc_obj$primary_cohorts,
    function(x) x$Sample
  )
  
  same_sample_set <- sapply(
    primary_sample_sets,
    function(x) identical(
      x,
      primary_sample_sets[[1]]
    )
  )
  
  cat("\nIdentical Sample ordering across primary cohorts:\n")
  print(same_sample_set)
  
  cat("\nUnique Sample IDs in ME_brown:\n")
  print(length(unique(primary_brown$Sample)))
  
  cat("\nDuplicated Sample IDs:\n")
  print(anyDuplicated(primary_brown$Sample))
}


cat("\n============================================================\n")
cat("STEP 6B-2B completed\n")
cat("============================================================\n")

# ============================================================
# STEP 6B-2C
# Build the formally locked 485/178 survival-expression dataset
# ============================================================

cat("\n============================================================\n")
cat("STEP 6B-2C: Build locked 485/178 survival-expression dataset\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Extract the LOCKED sample membership from 2N
# ------------------------------------------------------------

locked_primary_samples <- primary_brown$Sample

cat("Locked primary sample count:\n")
print(length(locked_primary_samples))

cat("Unique locked primary samples:\n")
print(length(unique(locked_primary_samples)))

cat("Duplicated locked primary samples:\n")
print(anyDuplicated(locked_primary_samples))

stopifnot(
  length(locked_primary_samples) == 485,
  length(unique(locked_primary_samples)) == 485,
  anyDuplicated(locked_primary_samples) == 0
)


# ------------------------------------------------------------
# 2. Match locked samples to clinical_aligned
# ------------------------------------------------------------

clinical_idx <- match(
  locked_primary_samples,
  clinical_aligned$Sample
)

cat("\nLocked samples matched to clinical_aligned:\n")
print(sum(!is.na(clinical_idx)))

cat("Locked samples NOT matched to clinical_aligned:\n")
print(sum(is.na(clinical_idx)))

if (anyNA(clinical_idx)) {
  cat("\nUnmatched samples:\n")
  print(
    locked_primary_samples[
      is.na(clinical_idx)
    ]
  )
}

stopifnot(all(!is.na(clinical_idx)))


# ------------------------------------------------------------
# 3. Extract clinical data in LOCKED cohort order
# ------------------------------------------------------------

primary_clinical_locked <- clinical_aligned[
  clinical_idx,
  ,
  drop = FALSE
]

cat("\nLocked clinical cohort dimensions:\n")
print(dim(primary_clinical_locked))


# ------------------------------------------------------------
# 4. Verify Sample order
# ------------------------------------------------------------

cat("\nSample order identical to locked 2N cohort:\n")
print(
  identical(
    primary_clinical_locked$Sample,
    locked_primary_samples
  )
)

stopifnot(
  identical(
    primary_clinical_locked$Sample,
    locked_primary_samples
  )
)


# ------------------------------------------------------------
# 5. Match locked samples to datExpr_MAD
# ------------------------------------------------------------

expr_idx_locked <- match(
  locked_primary_samples,
  rownames(datExpr_MAD)
)

cat("\nLocked samples matched to datExpr_MAD:\n")
print(sum(!is.na(expr_idx_locked)))

cat("Locked samples NOT matched to datExpr_MAD:\n")
print(sum(is.na(expr_idx_locked)))

if (anyNA(expr_idx_locked)) {
  cat("\nUnmatched expression samples:\n")
  print(
    locked_primary_samples[
      is.na(expr_idx_locked)
    ]
  )
}

stopifnot(all(!is.na(expr_idx_locked)))


# ------------------------------------------------------------
# 6. Extract the three core genes
# ------------------------------------------------------------

expr_core_locked <- datExpr_MAD[
  expr_idx_locked,
  unname(core_gene_ids),
  drop = FALSE
]

colnames(expr_core_locked) <- names(core_gene_ids)


# ------------------------------------------------------------
# 7. Construct the formal Step 6B dataset
# ------------------------------------------------------------

survival_gene_primary_locked <- cbind(
  primary_clinical_locked[
    ,
    c(
      "Sample",
      "patient_id",
      "Age",
      "Gender",
      "pathologic_stage",
      "pathologic_t",
      "pathologic_n",
      "pathologic_m",
      "OS.time",
      "OS.status"
    ),
    drop = FALSE
  ],
  as.data.frame(expr_core_locked)
)


# ------------------------------------------------------------
# 8. Verify survival outcome
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("LOCKED PRIMARY COHORT OUTCOME QC\n")
cat("============================================================\n\n")

cat("Dimensions:\n")
print(dim(survival_gene_primary_locked))

cat("\nOS.status:\n")
print(table(survival_gene_primary_locked$OS.status))

cat("\nOS.time summary:\n")
print(summary(survival_gene_primary_locked$OS.time))


# ------------------------------------------------------------
# 9. Verify no missing core-gene expression
# ------------------------------------------------------------

cat("\nMissing values in core genes:\n")
print(
  colSums(
    is.na(
      survival_gene_primary_locked[
        ,
        c("AOC3", "MS4A7", "IL7R")
      ]
    )
  )
)


# ------------------------------------------------------------
# 10. Cross-check survival data against the original 2N cohort
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CROSS-CHECK AGAINST LOCKED 2N COHORT\n")
cat("============================================================\n\n")

cat("Sample identical:\n")
print(
  identical(
    survival_gene_primary_locked$Sample,
    primary_brown$Sample
  )
)

cat("\nOS.time identical:\n")
print(
  identical(
    survival_gene_primary_locked$OS.time,
    primary_brown$OS_time
  )
)

cat("\nOS.status/event identical:\n")
print(
  identical(
    as.numeric(survival_gene_primary_locked$OS.status),
    as.numeric(primary_brown$OS_event)
  )
)


# ------------------------------------------------------------
# 11. Final locked-cohort assertions
# ------------------------------------------------------------

stopifnot(
  nrow(survival_gene_primary_locked) == 485,
  sum(survival_gene_primary_locked$OS.status == "Dead") == 178,
  sum(survival_gene_primary_locked$OS.status == "Alive") == 307,
  all(!is.na(survival_gene_primary_locked$AOC3)),
  all(!is.na(survival_gene_primary_locked$MS4A7)),
  all(!is.na(survival_gene_primary_locked$IL7R))
)


cat("\n============================================================\n")
cat("STEP 6B-2C PASSED\n")
cat("LOCKED PRIMARY COHORT: N = 485\n")
cat("DEATHS = 178\n")
cat("CENSORED = 307\n")
cat("AOC3 / MS4A7 / IL7R expression successfully mapped\n")
cat("============================================================\n")

# ============================================================
# STEP 6B-2D
# Final survival outcome consistency check
# ============================================================

cat("\n============================================================\n")
cat("STEP 6B-2D: Final survival outcome consistency check\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Convert Alive/Dead to the numeric 0/1 representation
# ------------------------------------------------------------

OS_event_locked <- ifelse(
  survival_gene_primary_locked$OS.status == "Dead",
  1,
  0
)


# ------------------------------------------------------------
# 2. Compare against the locked 2N OS_event
# ------------------------------------------------------------

cat("OS event values identical after explicit mapping:\n")

event_identical <- identical(
  as.numeric(OS_event_locked),
  as.numeric(primary_brown$OS_event)
)

print(event_identical)


# ------------------------------------------------------------
# 3. Report any discrepancies
# ------------------------------------------------------------

if (!event_identical) {
  
  discrepancy_idx <- which(
    OS_event_locked != primary_brown$OS_event
  )
  
  cat("\nNumber of discrepancies:\n")
  print(length(discrepancy_idx))
  
  if (length(discrepancy_idx) > 0) {
    cat("\nFirst discrepancies:\n")
    print(
      data.frame(
        Sample = primary_brown$Sample[discrepancy_idx],
        Current_OS_status =
          survival_gene_primary_locked$OS.status[discrepancy_idx],
        Current_event = OS_event_locked[discrepancy_idx],
        Locked_event = primary_brown$OS_event[discrepancy_idx],
        stringsAsFactors = FALSE
      )
    )
  }
  
} else {
  
  cat("\nNo survival-event discrepancies detected.\n")
}


# ------------------------------------------------------------
# 4. Verify OS time
# ------------------------------------------------------------

time_identical <- identical(
  as.numeric(survival_gene_primary_locked$OS.time),
  as.numeric(primary_brown$OS_time)
)

cat("\nOS time values identical:\n")
print(time_identical)


# ------------------------------------------------------------
# 5. Final cohort QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL STEP 6B PRIMARY COHORT QC\n")
cat("============================================================\n\n")

cat("N:\n")
print(nrow(survival_gene_primary_locked))

cat("Events:\n")
print(sum(OS_event_locked == 1))

cat("Censored:\n")
print(sum(OS_event_locked == 0))

cat("Zero-time observations:\n")
print(sum(survival_gene_primary_locked$OS.time == 0))

cat("\nEvent consistency:\n")
print(event_identical)

cat("Time consistency:\n")
print(time_identical)


stopifnot(
  nrow(survival_gene_primary_locked) == 485,
  sum(OS_event_locked == 1) == 178,
  sum(OS_event_locked == 0) == 307,
  event_identical,
  time_identical
)


cat("\n============================================================\n")
cat("STEP 6B-2D PASSED\n")
cat("Primary survival cohort is formally locked.\n")
cat("N = 485; Events = 178; Censored = 307\n")
cat("============================================================\n")

# ============================================================
# Paper 1 — STEP 6B-3
# Expression-scale and distribution QC
# ============================================================

cat("\n============================================================\n")
cat("STEP 6B-3: Expression-scale and distribution QC\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define the three core genes
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

print(core_genes)


# ------------------------------------------------------------
# 2. Check that all three expression columns exist
# ------------------------------------------------------------

expr_missing <- setdiff(
  core_genes,
  colnames(survival_gene_primary_locked)
)

cat("\nMissing expression columns:\n")
print(expr_missing)

stopifnot(length(expr_missing) == 0)


# ------------------------------------------------------------
# 3. Expression distribution summary
# ------------------------------------------------------------

expr_qc <- do.call(
  rbind,
  lapply(core_genes, function(g) {
    
    x <- survival_gene_primary_locked[[g]]
    
    data.frame(
      Gene = g,
      N = length(x),
      N_NA = sum(is.na(x)),
      N_NaN = sum(is.nan(x)),
      N_Inf = sum(is.infinite(x)),
      Mean = mean(x),
      SD = sd(x),
      Median = median(x),
      Q1 = unname(quantile(x, 0.25)),
      Q3 = unname(quantile(x, 0.75)),
      Min = min(x),
      Max = max(x),
      Variance = var(x),
      stringsAsFactors = FALSE
    )
  })
)

rownames(expr_qc) <- NULL


cat("\n============================================================\n")
cat("CORE GENE EXPRESSION QC\n")
cat("============================================================\n\n")

print(expr_qc)


# ------------------------------------------------------------
# 4. Explicit validity checks
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("EXPRESSION VALIDITY CHECKS\n")
cat("============================================================\n\n")

for (g in core_genes) {
  
  x <- survival_gene_primary_locked[[g]]
  
  cat("\nGene:", g, "\n")
  cat("  NA:", sum(is.na(x)), "\n")
  cat("  NaN:", sum(is.nan(x)), "\n")
  cat("  Inf:", sum(is.infinite(x)), "\n")
  cat("  SD:", sd(x), "\n")
  cat("  Zero variance:", var(x) == 0, "\n")
}


# ------------------------------------------------------------
# 5. Global assertions
# ------------------------------------------------------------

for (g in core_genes) {
  
  x <- survival_gene_primary_locked[[g]]
  
  stopifnot(
    length(x) == 485,
    sum(is.na(x)) == 0,
    sum(is.nan(x)) == 0,
    sum(is.infinite(x)) == 0,
    is.finite(mean(x)),
    is.finite(sd(x)),
    sd(x) > 0
  )
}


# ------------------------------------------------------------
# 6. Check expression scale
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("EXPRESSION-SCALE CHECK\n")
cat("============================================================\n\n")

cat(
  "The current expression matrix originates from the\n",
  "tumor-only WGCNA workflow based on FPKM-UQ expression.\n",
  "The observed values are already on a log-transformed scale.\n",
  "No additional log transformation will be applied.\n\n"
)

cat("Observed ranges:\n")

for (g in core_genes) {
  
  x <- survival_gene_primary_locked[[g]]
  
  cat(
    g, ": ",
    round(min(x), 4),
    " to ",
    round(max(x), 4),
    "\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 7. Final STEP 6B-3 checkpoint
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 6B-3 PASSED\n")
cat("============================================================\n\n")

cat("Primary survival cohort remains locked at N = 485.\n")
cat("All three core genes have complete, non-zero-variance expression.\n")
cat("No additional log transformation will be performed.\n")

# ============================================================
# Paper 1 — STEP 6B-4
# Standardization of core-gene expression
# ============================================================

cat("\n============================================================\n")
cat("STEP 6B-4: Standardization of core-gene expression\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Confirm primary cohort remains locked
# ------------------------------------------------------------

stopifnot(
  nrow(survival_gene_primary_locked) == 485,
  sum(OS_event_locked == 1) == 178,
  sum(OS_event_locked == 0) == 307
)

cat("Primary survival cohort:\n")
cat("  N =", nrow(survival_gene_primary_locked), "\n")
cat("  Events =", sum(OS_event_locked == 1), "\n")
cat("  Censored =", sum(OS_event_locked == 0), "\n")


# ------------------------------------------------------------
# 2. Create a separate standardized-expression object
# ------------------------------------------------------------

survival_gene_primary_z <- survival_gene_primary_locked


# ------------------------------------------------------------
# 3. Standardize each core gene
# ------------------------------------------------------------

z_parameters <- data.frame(
  Gene = core_genes,
  Mean = NA_real_,
  SD = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(core_genes)) {
  
  g <- core_genes[i]
  
  x <- survival_gene_primary_locked[[g]]
  
  mu <- mean(x)
  s  <- sd(x)
  
  survival_gene_primary_z[[g]] <-
    (x - mu) / s
  
  z_parameters$Mean[i] <- mu
  z_parameters$SD[i] <- s
}


# ------------------------------------------------------------
# 4. Verify standardization
# ------------------------------------------------------------

z_qc <- do.call(
  rbind,
  lapply(core_genes, function(g) {
    
    z <- survival_gene_primary_z[[g]]
    
    data.frame(
      Gene = g,
      Mean_original = z_parameters$Mean[
        z_parameters$Gene == g
      ],
      SD_original = z_parameters$SD[
        z_parameters$Gene == g
      ],
      Mean_Z = mean(z),
      SD_Z = sd(z),
      Min_Z = min(z),
      Max_Z = max(z),
      NA_Z = sum(is.na(z)),
      NaN_Z = sum(is.nan(z)),
      Inf_Z = sum(is.infinite(z)),
      stringsAsFactors = FALSE
    )
  })
)

rownames(z_qc) <- NULL


cat("\n============================================================\n")
cat("STANDARDIZATION QC\n")
cat("============================================================\n\n")

print(z_qc)


# ------------------------------------------------------------
# 5. Check that survival variables were not altered
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SURVIVAL VARIABLE INTEGRITY CHECK\n")
cat("============================================================\n\n")

cat(
  "OS.time identical:",
  identical(
    survival_gene_primary_z$OS.time,
    survival_gene_primary_locked$OS.time
  ),
  "\n"
)

cat(
  "OS.status identical:",
  identical(
    survival_gene_primary_z$OS.status,
    survival_gene_primary_locked$OS.status
  ),
  "\n"
)

cat(
  "Sample identical:",
  identical(
    survival_gene_primary_z$Sample,
    survival_gene_primary_locked$Sample
  ),
  "\n"
)


# ------------------------------------------------------------
# 6. Formal assertions
# ------------------------------------------------------------

for (g in core_genes) {
  
  z <- survival_gene_primary_z[[g]]
  
  stopifnot(
    length(z) == 485,
    sum(is.na(z)) == 0,
    sum(is.nan(z)) == 0,
    sum(is.infinite(z)) == 0,
    abs(mean(z)) < 1e-10,
    abs(sd(z) - 1) < 1e-10
  )
}

stopifnot(
  identical(
    survival_gene_primary_z$OS.time,
    survival_gene_primary_locked$OS.time
  ),
  identical(
    survival_gene_primary_z$OS.status,
    survival_gene_primary_locked$OS.status
  ),
  identical(
    survival_gene_primary_z$Sample,
    survival_gene_primary_locked$Sample
  )
)


# ------------------------------------------------------------
# 7. Final checkpoint
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 6B-4 PASSED\n")
cat("============================================================\n\n")

cat("Three core genes standardized using the locked N = 485 cohort.\n")
cat("Standardization definition: Z = (X - mean) / SD.\n")
cat("Cox interpretation will be HR per 1-SD increase in expression.\n")
cat("The locked survival cohort and original expression object were not modified.\n")


# ============================================================
# Paper 1 — STEP 6B-5
# Continuous univariable Cox regression + Schoenfeld PH test
# ============================================================

cat("\n============================================================\n")
cat("STEP 6B-5: Continuous univariable Cox + Schoenfeld PH\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Load required package
# ------------------------------------------------------------

library(survival)


# ------------------------------------------------------------
# 2. Confirm locked analysis dataset
# ------------------------------------------------------------

stopifnot(
  nrow(survival_gene_primary_z) == 485,
  sum(OS_event_locked == 1) == 178,
  sum(OS_event_locked == 0) == 307
)

stopifnot(
  identical(
    survival_gene_primary_z$OS.time,
    survival_gene_primary_locked$OS.time
  ),
  identical(
    survival_gene_primary_z$OS.status,
    survival_gene_primary_locked$OS.status
  ),
  identical(
    survival_gene_primary_z$Sample,
    survival_gene_primary_locked$Sample
  )
)


# ------------------------------------------------------------
# 3. Convert survival status explicitly to numeric 0/1
# ------------------------------------------------------------

survival_gene_primary_z$OS_event <- OS_event_locked


stopifnot(
  identical(
    as.numeric(survival_gene_primary_z$OS_event),
    as.numeric(OS_event_locked)
  ),
  sum(survival_gene_primary_z$OS_event == 1) == 178,
  sum(survival_gene_primary_z$OS_event == 0) == 307
)


# ------------------------------------------------------------
# 4. Fit continuous univariable Cox models
# ------------------------------------------------------------

cox_univariable_models <- list()

for (g in core_genes) {
  
  formula_g <- as.formula(
    paste0(
      "Surv(OS.time, OS_event) ~ `",
      g,
      "`"
    )
  )
  
  fit <- coxph(
    formula_g,
    data = survival_gene_primary_z,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
  
  cox_univariable_models[[g]] <- fit
}


# ------------------------------------------------------------
# 5. Extract Cox estimates
# ------------------------------------------------------------

cox_results_list <- lapply(
  core_genes,
  function(g) {
    
    fit <- cox_univariable_models[[g]]
    
    beta <- coef(fit)[1]
    se <- sqrt(vcov(fit)[1, 1])
    
    data.frame(
      Gene = g,
      N = fit$n,
      Events = fit$nevent,
      Beta = unname(beta),
      HR_per_1SD = unname(exp(beta)),
      CI_low = unname(exp(beta - 1.96 * se)),
      CI_high = unname(exp(beta + 1.96 * se)),
      P = unname(summary(fit)$coefficients[1, "Pr(>|z|)"]),
      stringsAsFactors = FALSE
    )
  }
)

survival_cox_univariable_continuous <-
  do.call(rbind, cox_results_list)

rownames(survival_cox_univariable_continuous) <- NULL


# ------------------------------------------------------------
# 6. Schoenfeld PH tests
# ------------------------------------------------------------

ph_results_list <- lapply(
  core_genes,
  function(g) {
    
    fit <- cox_univariable_models[[g]]
    
    zph <- cox.zph(
      fit,
      transform = "km"
    )
    
    gene_row <- which(
      rownames(zph$table) == g
    )
    
    if (length(gene_row) != 1) {
      gene_row <- 1
    }
    
    data.frame(
      Gene = g,
      PH_P = unname(zph$table[gene_row, "p"]),
      stringsAsFactors = FALSE
    )
  }
)

ph_results <-
  do.call(rbind, ph_results_list)

rownames(ph_results) <- NULL


# ------------------------------------------------------------
# 7. Merge Cox + PH results
# ------------------------------------------------------------

survival_cox_univariable_continuous <-
  merge(
    survival_cox_univariable_continuous,
    ph_results,
    by = "Gene",
    sort = FALSE
  )


# ------------------------------------------------------------
# 8. Reorder genes explicitly
# ------------------------------------------------------------

survival_cox_univariable_continuous$Gene <-
  factor(
    survival_cox_univariable_continuous$Gene,
    levels = core_genes
  )

survival_cox_univariable_continuous <-
  survival_cox_univariable_continuous[
    order(survival_cox_univariable_continuous$Gene),
  ]

survival_cox_univariable_continuous$Gene <-
  as.character(
    survival_cox_univariable_continuous$Gene
  )

rownames(survival_cox_univariable_continuous) <- NULL


# ------------------------------------------------------------
# 9. Display primary results
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CONTINUOUS UNIVARIABLE COX RESULTS\n")
cat("============================================================\n\n")

print(
  survival_cox_univariable_continuous,
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Formal cohort consistency checks
# ------------------------------------------------------------

stopifnot(
  all(
    survival_cox_univariable_continuous$N == 485
  ),
  all(
    survival_cox_univariable_continuous$Events == 178
  ),
  all(
    is.finite(
      survival_cox_univariable_continuous$HR_per_1SD
    )
  ),
  all(
    is.finite(
      survival_cox_univariable_continuous$CI_low
    )
  ),
  all(
    is.finite(
      survival_cox_univariable_continuous$CI_high
    )
  ),
  all(
    is.finite(
      survival_cox_univariable_continuous$P
    )
  ),
  all(
    is.finite(
      survival_cox_univariable_continuous$PH_P
    )
  )
)


# ------------------------------------------------------------
# 11. Final checkpoint
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 6B-5 PASSED\n")
cat("============================================================\n\n")

cat(
  "Continuous univariable Cox models completed for:\n",
  paste(core_genes, collapse = ", "),
  "\n\n"
)

cat(
  "Primary cohort: N = 485; Events = 178.\n"
)

cat(
  "Effect estimates are reported per 1-SD increase in expression.\n"
)

cat(
  "Schoenfeld proportional-hazards tests were performed for each gene.\n"
)

# ============================================================
# Paper 1 — STEP 6B-6
# MS4A7 time-varying coefficient analysis
# Corrected version: remove model=TRUE
# ============================================================

cat("\n============================================================\n")
cat("STEP 6B-6: MS4A7 time-varying coefficient analysis\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Required package
# ------------------------------------------------------------

library(survival)


# ------------------------------------------------------------
# 2. Confirm locked analysis dataset
# ------------------------------------------------------------

stopifnot(
  nrow(survival_gene_primary_z) == 485,
  sum(OS_event_locked == 1) == 178,
  sum(OS_event_locked == 0) == 307
)

stopifnot(
  identical(
    survival_gene_primary_z$OS.time,
    survival_gene_primary_locked$OS.time
  ),
  identical(
    survival_gene_primary_z$OS.status,
    survival_gene_primary_locked$OS.status
  ),
  identical(
    survival_gene_primary_z$Sample,
    survival_gene_primary_locked$Sample
  )
)


# ------------------------------------------------------------
# 3. Check MS4A7 expression
# ------------------------------------------------------------

stopifnot(
  "MS4A7" %in% colnames(survival_gene_primary_z),
  all(is.finite(survival_gene_primary_z$MS4A7)),
  sd(survival_gene_primary_z$MS4A7) > 0
)

cat("MS4A7 standardized expression:\n")
cat(
  "  Mean =",
  mean(survival_gene_primary_z$MS4A7),
  "\n"
)
cat(
  "  SD =",
  sd(survival_gene_primary_z$MS4A7),
  "\n"
)


# ------------------------------------------------------------
# 4. Define analysis dataset
# ------------------------------------------------------------

ms4a7_tv_data <- data.frame(
  OS_time = as.numeric(survival_gene_primary_z$OS.time),
  OS_event = as.numeric(OS_event_locked),
  MS4A7 = as.numeric(survival_gene_primary_z$MS4A7)
)

stopifnot(
  nrow(ms4a7_tv_data) == 485,
  sum(ms4a7_tv_data$OS_event == 1) == 178,
  sum(ms4a7_tv_data$OS_event == 0) == 307,
  all(is.finite(ms4a7_tv_data$OS_time)),
  all(is.finite(ms4a7_tv_data$MS4A7))
)


# ------------------------------------------------------------
# 5. Define time functions
# ------------------------------------------------------------

time_functions <- list(
  
  linear = function(t) {
    t
  },
  
  sqrt = function(t) {
    sqrt(t)
  },
  
  log1p = function(t) {
    log1p(t)
  }
  
)


# ------------------------------------------------------------
# 6. Fit time-varying coefficient models
# ------------------------------------------------------------

ms4a7_tv_models <- list()

for (fn in names(time_functions)) {
  
  cat("\nFitting:", fn, "time function\n")
  
  tv_function <- time_functions[[fn]]
  
  fit <- coxph(
    Surv(OS_time, OS_event) ~
      MS4A7 +
      tt(MS4A7),
    data = ms4a7_tv_data,
    ties = "efron",
    tt = function(x, t, ...) {
      x * tv_function(t)
    },
    x = TRUE,
    y = TRUE
  )
  
  ms4a7_tv_models[[fn]] <- fit
  
  cat(
    "  N =", fit$n,
    " Events =", fit$nevent,
    "\n"
  )
  
  cat(
    "  AIC =", AIC(fit),
    "\n"
  )
}


# ------------------------------------------------------------
# 7. Compare time functions
# ------------------------------------------------------------

ms4a7_tv_comparison <- do.call(
  rbind,
  lapply(names(ms4a7_tv_models), function(fn) {
    
    fit <- ms4a7_tv_models[[fn]]
    
    data.frame(
      Time_function = fn,
      N = fit$n,
      Events = fit$nevent,
      AIC = AIC(fit),
      LogLik = as.numeric(logLik(fit)),
      Concordance = unname(
        summary(fit)$concordance[1]
      ),
      stringsAsFactors = FALSE
    )
    
  })
)

rownames(ms4a7_tv_comparison) <- NULL

cat("\n============================================================\n")
cat("TIME-VARYING MODEL COMPARISON\n")
cat("============================================================\n\n")

print(ms4a7_tv_comparison)


# ------------------------------------------------------------
# 8. Print model summaries
# ------------------------------------------------------------

for (fn in names(ms4a7_tv_models)) {
  
  cat("\n============================================================\n")
  cat("MODEL:", fn, "\n")
  cat("============================================================\n\n")
  
  print(summary(ms4a7_tv_models[[fn]]))
  
}


# ------------------------------------------------------------
# 9. Select model with lowest AIC
# ------------------------------------------------------------

best_idx <- which.min(ms4a7_tv_comparison$AIC)

best_time_function <-
  ms4a7_tv_comparison$Time_function[best_idx]

best_tv_model <-
  ms4a7_tv_models[[best_time_function]]

cat("\n============================================================\n")
cat("BEST TIME-VARYING MODEL\n")
cat("============================================================\n\n")

cat(
  "Selected time function:",
  best_time_function,
  "\n"
)

cat(
  "AIC:",
  AIC(best_tv_model),
  "\n"
)


# ------------------------------------------------------------
# 10. Extract coefficients
# ------------------------------------------------------------

best_coef <- coef(best_tv_model)
best_vcov <- vcov(best_tv_model)

coef_names <- names(best_coef)

beta0_name <- coef_names[1]
beta1_name <- coef_names[2]

beta0 <- unname(
  best_coef[beta0_name]
)

beta1 <- unname(
  best_coef[beta1_name]
)

beta1_se <- sqrt(
  best_vcov[
    beta1_name,
    beta1_name
  ]
)

beta1_z <- beta1 / beta1_se

beta1_p <- 2 * pnorm(
  abs(beta1_z),
  lower.tail = FALSE
)


cat("\nCoefficient interpretation:\n")

cat(
  "  beta0 =",
  beta0,
  "\n"
)

cat(
  "  beta_time =",
  beta1,
  "\n"
)

cat(
  "  SE(beta_time) =",
  beta1_se,
  "\n"
)

cat(
  "  P(beta_time) =",
  beta1_p,
  "\n"
)


# ------------------------------------------------------------
# 11. Estimate HR at clinically interpretable time points
# ------------------------------------------------------------

time_points <- c(
  365,
  730,
  1095,
  1460,
  1825
)

f_best <- time_functions[[best_time_function]]

f_values <- f_best(time_points)

eta_t <- beta0 + beta1 * f_values

HR_CI <- t(
  sapply(
    seq_along(time_points),
    function(i) {
      
      f_i <- f_values[i]
      
      L <- c(
        1,
        f_i
      )
      
      se_eta <- sqrt(
        as.numeric(
          t(L) %*%
            best_vcov %*%
            L
        )
      )
      
      eta_i <- eta_t[i]
      
      c(
        HR = exp(eta_i),
        
        CI_low =
          exp(
            eta_i -
              1.96 * se_eta
          ),
        
        CI_high =
          exp(
            eta_i +
              1.96 * se_eta
          )
      )
      
    }
  )
)


ms4a7_hr_time <- data.frame(
  Time_days = time_points,
  
  Time_years =
    time_points / 365,
  
  HR =
    HR_CI[, "HR"],
  
  CI_low =
    HR_CI[, "CI_low"],
  
  CI_high =
    HR_CI[, "CI_high"],
  
  stringsAsFactors = FALSE
)


cat("\n============================================================\n")
cat("TIME-SPECIFIC HR FOR MS4A7\n")
cat("============================================================\n\n")

print(ms4a7_hr_time)


# ------------------------------------------------------------
# 12. Final result object
# ------------------------------------------------------------

ms4a7_tv_final <- list(
  
  data =
    ms4a7_tv_data,
  
  models =
    ms4a7_tv_models,
  
  comparison =
    ms4a7_tv_comparison,
  
  best_time_function =
    best_time_function,
  
  final_model =
    best_tv_model,
  
  beta0 =
    beta0,
  
  beta_time =
    beta1,
  
  beta_time_SE =
    beta1_se,
  
  beta_time_P =
    beta1_p,
  
  HR_time =
    ms4a7_hr_time
  
)


# ------------------------------------------------------------
# 13. Final QC
# ------------------------------------------------------------

stopifnot(
  ms4a7_tv_final$best_time_function %in%
    names(time_functions),
  
  ms4a7_tv_final$final_model$n == 485,
  
  ms4a7_tv_final$final_model$nevent == 178,
  
  nrow(ms4a7_tv_final$HR_time) == 5,
  
  all(
    is.finite(
      ms4a7_tv_final$HR_time$HR
    )
  ),
  
  all(
    is.finite(
      ms4a7_tv_final$HR_time$CI_low
    )
  ),
  
  all(
    is.finite(
      ms4a7_tv_final$HR_time$CI_high
    )
  )
)


cat("\n============================================================\n")
cat("STEP 6B-6 COMPLETED\n")
cat("============================================================\n")


# ============================================================
# Paper 1 — STEP 6B-6A
# Likelihood-ratio comparison:
# constant vs time-varying MS4A7 effect
# ============================================================

cat("\n============================================================\n")
cat("STEP 6B-6A: Constant vs time-varying MS4A7 model\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Constant-coefficient Cox model
# ------------------------------------------------------------

ms4a7_constant_model <- coxph(
  Surv(OS_time, OS_event) ~ MS4A7,
  data = ms4a7_tv_data,
  ties = "efron",
  x = TRUE,
  y = TRUE
)


# ------------------------------------------------------------
# 2. Best time-varying model
# ------------------------------------------------------------

ms4a7_sqrt_model <- ms4a7_tv_models[["sqrt"]]


# ------------------------------------------------------------
# 3. Model comparison
# ------------------------------------------------------------

cat("Constant-coefficient model:\n")
print(ms4a7_constant_model)

cat("\n\nTime-varying sqrt model:\n")
print(ms4a7_sqrt_model)


# ------------------------------------------------------------
# 4. Likelihood-ratio test
# ------------------------------------------------------------

ll_constant <- as.numeric(
  logLik(ms4a7_constant_model)
)

ll_tv <- as.numeric(
  logLik(ms4a7_sqrt_model)
)

LR_stat <- 2 * (
  ll_tv - ll_constant
)

df_LRT <- 1

P_LRT <- pchisq(
  LR_stat,
  df = df_LRT,
  lower.tail = FALSE
)


ms4a7_LRT <- data.frame(
  Model_comparison =
    "Constant vs sqrt(time)-varying",
  
  LogLik_constant =
    ll_constant,
  
  LogLik_time_varying =
    ll_tv,
  
  LR =
    LR_stat,
  
  df =
    df_LRT,
  
  P =
    P_LRT,
  
  stringsAsFactors = FALSE
)


cat("\n============================================================\n")
cat("LIKELIHOOD-RATIO TEST\n")
cat("============================================================\n\n")

print(ms4a7_LRT)


# ------------------------------------------------------------
# 5. AIC comparison
# ------------------------------------------------------------

ms4a7_AIC_comparison <- data.frame(
  
  Model = c(
    "Constant coefficient",
    "Time-varying sqrt(time)"
  ),
  
  AIC = c(
    AIC(ms4a7_constant_model),
    AIC(ms4a7_sqrt_model)
  ),
  
  LogLik = c(
    ll_constant,
    ll_tv
  ),
  
  stringsAsFactors = FALSE
)

ms4a7_AIC_comparison$Delta_AIC <-
  ms4a7_AIC_comparison$AIC -
  min(ms4a7_AIC_comparison$AIC)


cat("\n============================================================\n")
cat("AIC COMPARISON\n")
cat("============================================================\n\n")

print(ms4a7_AIC_comparison)


# ------------------------------------------------------------
# 6. Final QC
# ------------------------------------------------------------

stopifnot(
  ms4a7_constant_model$n == 485,
  ms4a7_constant_model$nevent == 178,
  
  ms4a7_sqrt_model$n == 485,
  ms4a7_sqrt_model$nevent == 178,
  
  is.finite(LR_stat),
  is.finite(P_LRT)
)


cat("\n============================================================\n")
cat("STEP 6B-6A COMPLETED\n")
cat("============================================================\n")


# ============================================================
# Diagnose missing Step 6B-6 checkpoint objects
# ============================================================

required_objects <- c(
  "ms4a7_tv_data",
  "ms4a7_tv_models",
  "ms4a7_tv_comparison",
  "best_time_function",
  "best_tv_model",
  "beta0",
  "beta1",
  "beta1_se",
  "beta1_p",
  "ms4a7_hr_time",
  "ms4a7_tv_final",
  "ms4a7_constant_model",
  "ms4a7_sqrt_model",
  "ms4a7_LRT",
  "ms4a7_AIC_comparison"
)

object_status <- data.frame(
  Object = required_objects,
  Exists = vapply(
    required_objects,
    exists,
    logical(1),
    inherits = FALSE
  ),
  stringsAsFactors = FALSE
)

print(object_status, row.names = FALSE)

missing_objects <- object_status$Object[!object_status$Exists]

cat("\n------------------------------------------------------------\n")
cat("Missing objects:\n")
print(missing_objects)

cat("\nNumber of missing objects:", length(missing_objects), "\n")
cat("------------------------------------------------------------\n")

# ============================================================
# Paper 1 — Step 6B recovery audit
# Check both project directories
# ============================================================

wgcna_dir <- "O:/Rproject/Paper1_tumor-only WGCNA"
survival_dir <- "O:/Rproject/Paper 1_Step 6-survival analysis"

cat("\n============================================================\n")
cat("WGCNA directory\n")
cat(wgcna_dir, "\n")
cat("============================================================\n")

stopifnot(dir.exists(wgcna_dir))

wgcna_files <- list.files(
  wgcna_dir,
  full.names = FALSE,
  recursive = FALSE
)

wgcna_selected <- wgcna_files[
  grepl(
    "Workflow-04-2N|Workflow-04-2O|Step.?6|survival|Cox|PH|checkpoint|MS4A7",
    wgcna_files,
    ignore.case = TRUE
  )
]

print(wgcna_selected)


cat("\n============================================================\n")
cat("Step 6 survival directory\n")
cat(survival_dir, "\n")
cat("============================================================\n")

stopifnot(dir.exists(survival_dir))

survival_files <- list.files(
  survival_dir,
  full.names = FALSE,
  recursive = FALSE
)

survival_selected <- survival_files[
  grepl(
    "Step.?6|survival|Cox|PH|checkpoint|MS4A7",
    survival_files,
    ignore.case = TRUE
  )
]

print(survival_selected)


# ============================================================
# Paper 1 — Step 6B recovery
# Audit existing upstream checkpoints WITHOUT modifying
# the current R environment
# ============================================================

wgcna_dir <- "O:/Rproject/Paper1_tumor-only WGCNA"


# ------------------------------------------------------------
# 1. Audit cohort QC checkpoint
# ------------------------------------------------------------

cohort_file <- file.path(
  wgcna_dir,
  "Workflow-04-2N_WGCNA_OS_multivariable_cohort_QC.rds"
)

stopifnot(
  file.exists(cohort_file)
)

cohort_qc_recovery <- readRDS(cohort_file)

cat("\n============================================================\n")
cat("Workflow-04-2N cohort checkpoint\n")
cat("============================================================\n\n")

cat("Top-level objects:\n")
print(names(cohort_qc_recovery))

cat("\nStructure:\n")
str(
  cohort_qc_recovery,
  max.level = 2
)


# ------------------------------------------------------------
# 2. Audit WGCNA input checkpoint
# ------------------------------------------------------------

wgcna_file <- file.path(
  wgcna_dir,
  "Workflow-04-2D_final_WGCNA_candidate_checkpoint.RData"
)

stopifnot(
  file.exists(wgcna_file)
)

wgcna_recovery_env <- new.env(parent = emptyenv())

load(
  wgcna_file,
  envir = wgcna_recovery_env
)

cat("\n============================================================\n")
cat("Workflow-04-2D WGCNA checkpoint\n")
cat("============================================================\n\n")

wgcna_objects <- ls(
  envir = wgcna_recovery_env
)

cat("Objects in checkpoint:\n")
print(wgcna_objects)


# ------------------------------------------------------------
# 3. Check key expected objects
# ------------------------------------------------------------

key_objects <- c(
  "clinical_aligned",
  "datExpr_MAD",
  "moduleColors",
  "MM_all",
  "GS_all"
)

key_status <- data.frame(
  Object = key_objects,
  Exists = vapply(
    key_objects,
    exists,
    logical(1),
    envir = wgcna_recovery_env,
    inherits = FALSE
  ),
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("Key WGCNA checkpoint objects\n")
cat("============================================================\n\n")

print(
  key_status,
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. Inspect primary cohort
# ------------------------------------------------------------

if (
  "primary_cohorts" %in%
  names(cohort_qc_recovery)
) {
  
  cat("\n============================================================\n")
  cat("Primary cohort structure\n")
  cat("============================================================\n\n")
  
  print(
    names(cohort_qc_recovery$primary_cohorts)
  )
  
  primary_brown_recovery <-
    cohort_qc_recovery$primary_cohorts$ME_brown
  
  cat("\nME_brown primary cohort dimensions:\n")
  print(dim(primary_brown_recovery))
  
  cat("\nME_brown primary cohort columns:\n")
  print(colnames(primary_brown_recovery))
  
  cat("\nEvent distribution:\n")
  print(table(primary_brown_recovery$OS_event))
  
  cat("\nSample duplicates:\n")
  print(anyDuplicated(primary_brown_recovery$Sample))
}


# ============================================================
# Paper 1 — STEP 6B recovery
# Reconstruct the formally locked 485-sample survival cohort
# ============================================================

cat("\n============================================================\n")
cat("STEP 6B recovery: reconstruct locked survival cohort\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Extract the formally locked 485 samples
# ------------------------------------------------------------

locked_primary_samples <-
  primary_brown_recovery$Sample

stopifnot(
  length(locked_primary_samples) == 485,
  length(unique(locked_primary_samples)) == 485
)

cat("Locked sample count:", length(locked_primary_samples), "\n")


# ------------------------------------------------------------
# 2. Recover clinical and expression objects
# ------------------------------------------------------------

clinical_aligned_recovery <-
  wgcna_recovery_env$clinical_aligned

datExpr_MAD_recovery <-
  wgcna_recovery_env$datExpr_MAD

gene_MAD_recovery <-
  wgcna_recovery_env$gene_MAD


# ------------------------------------------------------------
# 3. Confirm sample overlap
# ------------------------------------------------------------

clinical_idx <- match(
  locked_primary_samples,
  clinical_aligned_recovery$Sample
)

expr_idx <- match(
  locked_primary_samples,
  rownames(datExpr_MAD_recovery)
)

stopifnot(
  all(!is.na(clinical_idx)),
  all(!is.na(expr_idx))
)

cat(
  "Clinical samples matched:",
  sum(!is.na(clinical_idx)),
  "\n"
)

cat(
  "Expression samples matched:",
  sum(!is.na(expr_idx)),
  "\n"
)


# ------------------------------------------------------------
# 4. Recover clinical data
# ------------------------------------------------------------

primary_clinical_locked <-
  clinical_aligned_recovery[
    clinical_idx,
    ,
    drop = FALSE
  ]

rownames(primary_clinical_locked) <-
  primary_clinical_locked$Sample


# ------------------------------------------------------------
# 5. Recover core-gene expression
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_ensembl <- c(
  AOC3  = "ENSG00000131471",
  MS4A7 = "ENSG00000166927",
  IL7R  = "ENSG00000168685"
)

stopifnot(
  all(core_ensembl %in% colnames(datExpr_MAD_recovery))
)

expr_core_locked <-
  datExpr_MAD_recovery[
    expr_idx,
    unname(core_ensembl),
    drop = FALSE
  ]

colnames(expr_core_locked) <-
  core_genes


# ------------------------------------------------------------
# 6. Build locked survival dataset
# ------------------------------------------------------------

survival_gene_primary_locked <-
  data.frame(
    Sample = locked_primary_samples,
    primary_clinical_locked[
      ,
      c(
        "OS.time",
        "OS.status"
      ),
      drop = FALSE
    ],
    expr_core_locked,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

rownames(survival_gene_primary_locked) <-
  survival_gene_primary_locked$Sample


# ------------------------------------------------------------
# 7. Recover event variable
# ------------------------------------------------------------

OS_event_locked <- ifelse(
  survival_gene_primary_locked$OS.status == "Dead",
  1,
  0
)

OS_event_locked <-
  as.numeric(OS_event_locked)


# ------------------------------------------------------------
# 8. Final cohort QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Recovered locked cohort QC\n")
cat("============================================================\n\n")

cat(
  "N:",
  nrow(survival_gene_primary_locked),
  "\n"
)

cat(
  "Events:",
  sum(OS_event_locked),
  "\n"
)

cat(
  "Censored:",
  sum(OS_event_locked == 0),
  "\n"
)

cat(
  "Zero-time:",
  sum(survival_gene_primary_locked$OS.time == 0),
  "\n"
)

cat(
  "Duplicate samples:",
  anyDuplicated(survival_gene_primary_locked$Sample),
  "\n"
)

cat("\nEvent distribution:\n")
print(table(OS_event_locked))


# ------------------------------------------------------------
# 9. Compare with formal ME_brown cohort
# ------------------------------------------------------------

event_identical <- identical(
  as.numeric(OS_event_locked),
  as.numeric(primary_brown_recovery$OS_event)
)

time_identical <- identical(
  as.numeric(
    survival_gene_primary_locked$OS.time
  ),
  as.numeric(
    primary_brown_recovery$OS_time
  )
)

sample_identical <- identical(
  survival_gene_primary_locked$Sample,
  primary_brown_recovery$Sample
)

cat("\nEvent consistency:", event_identical, "\n")
cat("Time consistency:", time_identical, "\n")
cat("Sample consistency:", sample_identical, "\n")

# ============================================================
# STEP 6B-4 recovery
# Standardize core-gene expression within locked N=485 cohort
# ============================================================

survival_gene_primary_z <-
  survival_gene_primary_locked

z_parameters <- data.frame(
  Gene = core_genes,
  Mean_original = NA_real_,
  SD_original = NA_real_,
  Mean_Z = NA_real_,
  SD_Z = NA_real_,
  stringsAsFactors = FALSE
)

for (g in core_genes) {
  
  x <- survival_gene_primary_locked[[g]]
  
  mu <- mean(x)
  s  <- sd(x)
  
  survival_gene_primary_z[[g]] <-
    (x - mu) / s
  
  z_parameters[
    z_parameters$Gene == g,
    "Mean_original"
  ] <- mu
  
  z_parameters[
    z_parameters$Gene == g,
    "SD_original"
  ] <- s
  
  z_parameters[
    z_parameters$Gene == g,
    "Mean_Z"
  ] <- mean(
    survival_gene_primary_z[[g]]
  )
  
  z_parameters[
    z_parameters$Gene == g,
    "SD_Z"
  ] <- sd(
    survival_gene_primary_z[[g]]
  )
}

print(z_parameters)



# ============================================================
# Paper 1 — STEP 6B-5
# Rebuild continuous-expression univariable Cox + Schoenfeld PH
# ============================================================

library(survival)

cat("\n============================================================\n")
cat("STEP 6B-5: Continuous-expression Cox + Schoenfeld PH\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Prepare locked event variable
# ------------------------------------------------------------

survival_gene_primary_z$OS_event <-
  as.numeric(OS_event_locked)


# ------------------------------------------------------------
# 2. Confirm locked cohort
# ------------------------------------------------------------

stopifnot(
  nrow(survival_gene_primary_z) == 485,
  sum(survival_gene_primary_z$OS_event) == 178,
  sum(survival_gene_primary_z$OS_event == 0) == 307
)


# ------------------------------------------------------------
# 3. Fit univariable Cox models
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

cox_univariable_models <- list()
cox_univariable_ph <- list()

survival_cox_univariable_continuous <- data.frame(
  Gene = core_genes,
  N = NA_integer_,
  Events = NA_integer_,
  Beta = NA_real_,
  HR_per_1SD = NA_real_,
  CI_low = NA_real_,
  CI_high = NA_real_,
  P = NA_real_,
  PH_P = NA_real_,
  stringsAsFactors = FALSE
)


for (i in seq_along(core_genes)) {
  
  g <- core_genes[i]
  
  formula_g <- as.formula(
    paste0(
      "Surv(OS.time, OS_event) ~ ",
      g
    )
  )
  
  fit <- coxph(
    formula_g,
    data = survival_gene_primary_z,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
  
  zph <- cox.zph(
    fit,
    transform = "km"
  )
  
  cox_univariable_models[[g]] <- fit
  cox_univariable_ph[[g]] <- zph
  
  sm <- summary(fit)
  
  survival_cox_univariable_continuous$N[i] <-
    fit$n
  
  survival_cox_univariable_continuous$Events[i] <-
    fit$nevent
  
  survival_cox_univariable_continuous$Beta[i] <-
    coef(fit)[g]
  
  survival_cox_univariable_continuous$HR_per_1SD[i] <-
    exp(coef(fit)[g])
  
  ci <- exp(
    confint(fit)[g, ]
  )
  
  survival_cox_univariable_continuous$CI_low[i] <-
    ci[1]
  
  survival_cox_univariable_continuous$CI_high[i] <-
    ci[2]
  
  survival_cox_univariable_continuous$P[i] <-
    sm$coefficients[g, "Pr(>|z|)"]
  
  survival_cox_univariable_continuous$PH_P[i] <-
    zph$table[g, "p"]
}


# ------------------------------------------------------------
# 4. Print final results
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Continuous-expression univariable Cox results\n")
cat("============================================================\n\n")

print(
  survival_cox_univariable_continuous,
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Print Schoenfeld PH tests
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Schoenfeld proportional-hazards tests\n")
cat("============================================================\n\n")

for (g in core_genes) {
  
  cat("\nGene:", g, "\n")
  
  print(
    cox_univariable_ph[[g]]
  )
}


# ------------------------------------------------------------
# 6. Expected-value QC
# ------------------------------------------------------------

expected_hr <- c(
  AOC3 = 0.8694600,
  MS4A7 = 0.8530846,
  IL7R = 0.8469172
)

expected_ph <- c(
  AOC3 = 0.65949231,
  MS4A7 = 0.03130103,
  IL7R = 0.05685436
)

for (g in core_genes) {
  
  actual_hr <-
    survival_cox_univariable_continuous[
      survival_cox_univariable_continuous$Gene == g,
      "HR_per_1SD"
    ]
  
  actual_ph <-
    survival_cox_univariable_continuous[
      survival_cox_univariable_continuous$Gene == g,
      "PH_P"
    ]
  
  cat(
    "\n",
    g,
    ": HR = ",
    actual_hr,
    " | expected = ",
    expected_hr[g],
    "\n",
    sep = ""
  )
  
  cat(
    g,
    ": PH P = ",
    actual_ph,
    " | expected = ",
    expected_ph[g],
    "\n",
    sep = ""
  )
  
  stopifnot(
    isTRUE(
      all.equal(
        as.numeric(actual_hr),
        as.numeric(expected_hr[g]),
        tolerance = 1e-6
      )
    )
  )
  
  stopifnot(
    isTRUE(
      all.equal(
        as.numeric(actual_ph),
        as.numeric(expected_ph[g]),
        tolerance = 1e-6
      )
    )
  )
}

cat(
  "\n============================================================\n"
)

cat(
  "STEP 6B-5 recovery QC PASSED.\n"
)

cat(
  "All HR and Schoenfeld PH values match the locked results.\n"
)

cat(
  "============================================================\n"
)


# ============================================================
# STEP 6B-5 CHECKPOINT
# Continuous-expression Cox + Schoenfeld PH
# ============================================================

step6B5_objects <- c(
  "survival_gene_primary_z",
  "z_parameters",
  "OS_event_locked",
  "cox_univariable_models",
  "cox_univariable_ph",
  "survival_cox_univariable_continuous"
)

missing_step6B5 <- step6B5_objects[
  !vapply(
    step6B5_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_step6B5) > 0) {
  stop(
    "Missing Step 6B-5 objects: ",
    paste(missing_step6B5, collapse = ", ")
  )
}

save(
  list = step6B5_objects,
  file = "Paper1_Step6B5_continuous_Cox_PH_checkpoint.RData"
)

cat("\n============================================================\n")
cat("STEP 6B-5 CHECKPOINT SAVED.\n")
cat("============================================================\n")
cat(
  "File: Paper1_Step6B5_continuous_Cox_PH_checkpoint.RData\n"
)
cat(
  "Cohort: N = ",
  nrow(survival_gene_primary_z),
  "; Events = ",
  sum(survival_gene_primary_z$OS_event),
  "\n",
  sep = ""
)
cat("Genes: AOC3, MS4A7, IL7R\n")
cat("Continuous expression: 1-SD standardized\n")
cat("PH assessment: Schoenfeld residuals\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 6B-6
# MS4A7 Time-Varying Cox Analysis
#
# Purpose:
#   1. Fit prespecified time-varying coefficient models
#      using linear, sqrt, and log1p time functions
#   2. Compare model fit using AIC
#   3. Compare constant vs sqrt-time model using LRT
#   4. Estimate time-specific HRs at 1/2/3/4/5 years
#   5. Perform complete QC against locked results
#   6. Save a formal checkpoint
#
# Locked upstream:
#   STEP 6B-5
#   N = 485
#   Events = 178
#   Core gene = MS4A7
# ============================================================


cat("\n============================================================\n")
cat("STEP 6B-6: MS4A7 time-varying Cox analysis\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Load required package
# ------------------------------------------------------------

library(survival)


# ------------------------------------------------------------
# 1. Confirm required STEP 6B-5 objects exist
# ------------------------------------------------------------

required_objects_6B5 <- c(
  "survival_gene_primary_z",
  "cox_univariable_models",
  "cox_univariable_ph",
  "survival_cox_univariable_continuous"
)

missing_objects_6B5 <- required_objects_6B5[
  !vapply(
    required_objects_6B5,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects_6B5) > 0) {
  
  stop(
    "STEP 6B-5 objects missing: ",
    paste(missing_objects_6B5, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 2. Confirm locked MS4A7 cohort
# ------------------------------------------------------------

stopifnot(
  nrow(survival_gene_primary_z) == 485,
  sum(survival_gene_primary_z$OS_event) == 178,
  sum(survival_gene_primary_z$OS_event == 0) == 307
)

stopifnot(
  "MS4A7" %in% colnames(survival_gene_primary_z),
  "OS.time" %in% colnames(survival_gene_primary_z),
  "OS_event" %in% colnames(survival_gene_primary_z)
)


cat("Locked cohort QC:\n")
cat("  N       =", nrow(survival_gene_primary_z), "\n")
cat("  Events  =", sum(survival_gene_primary_z$OS_event), "\n")
cat("  Censored=", sum(survival_gene_primary_z$OS_event == 0), "\n\n")


# ------------------------------------------------------------
# 3. Check survival-time range
# ------------------------------------------------------------

cat("OS-time summary:\n")
print(summary(survival_gene_primary_z$OS.time))

cat("\nZero-time observations = ",
    sum(survival_gene_primary_z$OS.time == 0),
    "\n\n",
    sep = "")


# ------------------------------------------------------------
# 4. Check MS4A7 expression
# ------------------------------------------------------------

stopifnot(
  all(is.finite(survival_gene_primary_z$MS4A7)),
  abs(mean(survival_gene_primary_z$MS4A7)) < 1e-10,
  abs(sd(survival_gene_primary_z$MS4A7) - 1) < 1e-10
)

cat("MS4A7 standardization QC:\n")
cat("  Mean =", mean(survival_gene_primary_z$MS4A7), "\n")
cat("  SD   =", sd(survival_gene_primary_z$MS4A7), "\n\n")


# ------------------------------------------------------------
# 5. Create independent analysis copy
#
# IMPORTANT:
# Do not modify the locked STEP 6B-5 objects.
# ------------------------------------------------------------

ms4a7_tv_data <- survival_gene_primary_z


# ------------------------------------------------------------
# 6. Define prespecified time functions
#
# Linear:
#     t
#
# Square-root:
#     sqrt(t)
#
# Log1p:
#     log(1+t)
#
# log1p() is used because the cohort contains zero-time
# observations and log(0) would be undefined.
# ------------------------------------------------------------

# ============================================================
# STEP 6B-6
# Correct tt() time-varying functions
# ============================================================

time_functions <- list(
  
  linear = function(x, t, ...) {
    x * t
  },
  
  sqrt = function(x, t, ...) {
    x * sqrt(t)
  },
  
  log1p = function(x, t, ...) {
    x * log1p(t)
  }
  
)


cat("\n============================================================\n")
cat("STEP 6B-6: Correct tt() time-varying functions\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# Refit all three prespecified models
# ------------------------------------------------------------

ms4a7_tv_models <- list()


for (fn_name in names(time_functions)) {
  
  cat("\n------------------------------------------------------------\n")
  cat("Fitting MS4A7 time-varying model:", fn_name, "\n")
  cat("------------------------------------------------------------\n")
  
  tv_formula <- as.formula(
    "Surv(OS.time, OS_event) ~
     MS4A7 +
     tt(MS4A7)"
  )
  
  fit_tv <- coxph(
    tv_formula,
    data = ms4a7_tv_data,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    tt = time_functions[[fn_name]]
  )
  
  ms4a7_tv_models[[fn_name]] <- fit_tv
  
  print(summary(fit_tv))
}


# ------------------------------------------------------------
# 8. Extract model-comparison statistics
# ------------------------------------------------------------

ms4a7_tv_comparison <- data.frame(
  Time_function = names(ms4a7_tv_models),
  N = NA_integer_,
  Events = NA_integer_,
  AIC = NA_real_,
  LogLik = NA_real_,
  Concordance = NA_real_,
  stringsAsFactors = FALSE
)


for (i in seq_along(ms4a7_tv_models)) {
  
  fn_name <- names(ms4a7_tv_models)[i]
  fit <- ms4a7_tv_models[[fn_name]]
  
  ms4a7_tv_comparison$N[i] <-
    fit$n
  
  ms4a7_tv_comparison$Events[i] <-
    fit$nevent
  
  ms4a7_tv_comparison$AIC[i] <-
    AIC(fit)
  
  ms4a7_tv_comparison$LogLik[i] <-
    as.numeric(logLik(fit))
  
  ms4a7_tv_comparison$Concordance[i] <-
    summary(fit)$concordance[1]
}


cat("\n============================================================\n")
cat("MS4A7 time-varying model comparison\n")
cat("============================================================\n\n")

print(
  ms4a7_tv_comparison,
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Select time function by minimum AIC
#
# This is NOT a post-hoc clinical outcome selection.
# The three functions were prespecified.
# ------------------------------------------------------------

best_time_function <-
  ms4a7_tv_comparison$Time_function[
    which.min(ms4a7_tv_comparison$AIC)
  ]


best_tv_model <-
  ms4a7_tv_models[[best_time_function]]


cat("\nSelected time function by minimum AIC:\n")
cat("  ", best_time_function, "\n", sep = "")

cat(
  "  AIC = ",
  AIC(best_tv_model),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Extract coefficients from selected model
#
# For the sqrt model:
#
#   beta(t) = beta0 + beta1 * sqrt(t)
#
# and
#
#   HR(t) = exp(beta0 + beta1 * sqrt(t))
# ------------------------------------------------------------

best_coef <- coef(best_tv_model)

beta0 <- unname(
  best_coef["MS4A7"]
)

beta1 <- unname(
  best_coef["tt(MS4A7)"]
)

best_coef_se <- sqrt(
  diag(vcov(best_tv_model))
)

beta0_se <- unname(
  best_coef_se["MS4A7"]
)

beta1_se <- unname(
  best_coef_se["tt(MS4A7)"]
)

best_coef_p <- summary(best_tv_model)$coefficients[
  ,
  "Pr(>|z|)"
]

beta0_p <- unname(
  best_coef_p["MS4A7"]
)

beta1_p <- unname(
  best_coef_p["tt(MS4A7)"]
)


cat("\n============================================================\n")
cat("Selected time-varying model coefficients\n")
cat("============================================================\n\n")

cat("Time function:", best_time_function, "\n\n")

cat(
  "MS4A7 coefficient (beta0) = ",
  beta0,
  "\n",
  sep = ""
)

cat(
  "MS4A7 beta0 P = ",
  beta0_p,
  "\n",
  sep = ""
)

cat(
  "Time-varying coefficient (beta1) = ",
  beta1,
  "\n",
  sep = ""
)

cat(
  "Time-varying coefficient SE = ",
  beta1_se,
  "\n",
  sep = ""
)

cat(
  "Time-varying coefficient P = ",
  beta1_p,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Fit constant-coefficient MS4A7 model
#
# This should reproduce the STEP 6B-5 MS4A7 Cox model.
# ------------------------------------------------------------

ms4a7_constant_model <-
  cox_univariable_models[["MS4A7"]]


# ------------------------------------------------------------
# 12. Confirm constant model matches STEP 6B-5
# ------------------------------------------------------------

constant_beta <-
  unname(coef(ms4a7_constant_model)["MS4A7"])

constant_hr <-
  exp(constant_beta)

constant_aic <-
  AIC(ms4a7_constant_model)

constant_logLik <-
  as.numeric(logLik(ms4a7_constant_model))


cat("\n============================================================\n")
cat("Constant-coefficient MS4A7 model\n")
cat("============================================================\n\n")

cat("Beta =", constant_beta, "\n")
cat("HR   =", constant_hr, "\n")
cat("AIC  =", constant_aic, "\n")
cat("LogLik =", constant_logLik, "\n")


# ------------------------------------------------------------
# 13. Constant vs selected time-varying model
#     Likelihood-ratio test
# ------------------------------------------------------------

ms4a7_LRT_stat <-
  2 * (
    as.numeric(logLik(best_tv_model)) -
      as.numeric(logLik(ms4a7_constant_model))
  )

ms4a7_LRT_df <-
  attr(logLik(best_tv_model), "df") -
  attr(logLik(ms4a7_constant_model), "df")

ms4a7_LRT_p <-
  pchisq(
    ms4a7_LRT_stat,
    df = ms4a7_LRT_df,
    lower.tail = FALSE
  )


ms4a7_LRT <- data.frame(
  Constant_LogLik =
    as.numeric(logLik(ms4a7_constant_model)),
  
  TimeVarying_LogLik =
    as.numeric(logLik(best_tv_model)),
  
  LR =
    ms4a7_LRT_stat,
  
  df =
    ms4a7_LRT_df,
  
  P =
    ms4a7_LRT_p,
  
  stringsAsFactors = FALSE
)


cat("\n============================================================\n")
cat("Constant vs time-varying likelihood-ratio test\n")
cat("============================================================\n\n")

print(
  ms4a7_LRT,
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. AIC comparison: constant vs selected time-varying
# ------------------------------------------------------------

ms4a7_AIC_comparison <- data.frame(
  
  Model = c(
    "Constant coefficient",
    paste0(
      "Time-varying: ",
      best_time_function
    )
  ),
  
  AIC = c(
    AIC(ms4a7_constant_model),
    AIC(best_tv_model)
  ),
  
  LogLik = c(
    as.numeric(logLik(ms4a7_constant_model)),
    as.numeric(logLik(best_tv_model))
  ),
  
  stringsAsFactors = FALSE
)


ms4a7_AIC_delta <-
  AIC(ms4a7_constant_model) -
  AIC(best_tv_model)


cat("\n============================================================\n")
cat("Constant vs time-varying AIC comparison\n")
cat("============================================================\n\n")

print(
  ms4a7_AIC_comparison,
  row.names = FALSE
)

cat(
  "\nDelta AIC (constant - time-varying) = ",
  ms4a7_AIC_delta,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 15. Time-specific HRs
#
# Clinically interpretable times:
#   1, 2, 3, 4, 5 years
#
# TCGA OS time is in days.
# Therefore:
#
#   1 year  = 365.25 days
#   2 years = 730.50 days
#   ...
# ------------------------------------------------------------

time_years <- c(1, 2, 3, 4, 5)

time_days <- time_years * 365.25


# The selected model is expected to be sqrt.
# To keep the script general, calculate using the actual
# selected time function.

time_transform <- switch(
  best_time_function,
  
  linear = time_days,
  
  sqrt = sqrt(time_days),
  
  log1p = log1p(time_days)
)


beta_time_specific <-
  beta0 +
  beta1 * time_transform


hr_time_specific <-
  exp(beta_time_specific)


# ------------------------------------------------------------
# 16. Time-specific 95% CIs
#
# Var(beta0 + beta1*f(t))
#
# = Var(beta0)
#   + f(t)^2 Var(beta1)
#   + 2*f(t)*Cov(beta0,beta1)
# ------------------------------------------------------------

V <- vcov(best_tv_model)


var_beta0 <-
  V["MS4A7", "MS4A7"]

var_beta1 <-
  V["tt(MS4A7)", "tt(MS4A7)"]

cov_beta01 <-
  V["MS4A7", "tt(MS4A7)"]


se_time_specific <- sqrt(
  
  var_beta0 +
    (time_transform^2) * var_beta1 +
    2 * time_transform * cov_beta01
  
)


CI_low_time <-
  exp(
    beta_time_specific -
      qnorm(0.975) * se_time_specific
  )


CI_high_time <-
  exp(
    beta_time_specific +
      qnorm(0.975) * se_time_specific
  )


ms4a7_hr_time <- data.frame(
  
  Time_years = time_years,
  
  Time_days = time_days,
  
  HR_per_1SD = hr_time_specific,
  
  CI_low = CI_low_time,
  
  CI_high = CI_high_time,
  
  stringsAsFactors = FALSE
)


cat("\n============================================================\n")
cat("MS4A7 time-specific hazard ratios\n")
cat("============================================================\n\n")

print(
  ms4a7_hr_time,
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Expected-result QC
#
# These values are the previously locked STEP 6B-6 results.
# They are used ONLY for QC, not for model construction.
# ------------------------------------------------------------

expected_time_function <- "sqrt"

expected_beta0 <- -0.4711206

expected_beta1 <- 0.01216699

expected_beta1_se <- 0.006008637

expected_beta1_p <- 0.04287589

expected_LRT <- 4.222337

expected_LRT_p <- 0.03989518

expected_delta_AIC <- 2.222337

expected_hr <- c(
  `1` = 0.7876743,
  `2` = 0.8672859,
  `3` = 0.9337885,
  `4` = 0.9937988,
  `5` = 1.0498565
)


# ------------------------------------------------------------
# 18. QC: selected time function
# ------------------------------------------------------------

stopifnot(
  identical(
    best_time_function,
    expected_time_function
  )
)


# ------------------------------------------------------------
# 19. QC: coefficients
# ------------------------------------------------------------

stopifnot(
  isTRUE(
    all.equal(
      as.numeric(beta0),
      expected_beta0,
      tolerance = 1e-6
    )
  )
)

stopifnot(
  isTRUE(
    all.equal(
      as.numeric(beta1),
      expected_beta1,
      tolerance = 1e-6
    )
  )
)

stopifnot(
  isTRUE(
    all.equal(
      as.numeric(beta1_se),
      expected_beta1_se,
      tolerance = 1e-6
    )
  )
)

stopifnot(
  isTRUE(
    all.equal(
      as.numeric(beta1_p),
      expected_beta1_p,
      tolerance = 1e-6
    )
  )
)


# ------------------------------------------------------------
# 20. QC: LRT
# ------------------------------------------------------------

stopifnot(
  isTRUE(
    all.equal(
      as.numeric(ms4a7_LRT$LR),
      expected_LRT,
      tolerance = 1e-6
    )
  )
)

stopifnot(
  isTRUE(
    all.equal(
      as.numeric(ms4a7_LRT$P),
      expected_LRT_p,
      tolerance = 1e-6
    )
  )
)


# ------------------------------------------------------------
# 21. QC: AIC improvement
# ------------------------------------------------------------

stopifnot(
  isTRUE(
    all.equal(
      as.numeric(ms4a7_AIC_delta),
      expected_delta_AIC,
      tolerance = 1e-6
    )
  )
)


# ------------------------------------------------------------
# 22. QC: time-specific HRs
#
# The time-specific HRs are deterministic functions of the
# locked beta0, beta1, and prespecified time points.
# Therefore, QC is based on numerical reproducibility rather
# than requiring exact identity with previously printed values.
# ------------------------------------------------------------

expected_hr_previous <- c(
  `1` = 0.7876743,
  `2` = 0.8672859,
  `3` = 0.9337885,
  `4` = 0.9937988,
  `5` = 1.0498565
)

hr_abs_diff <-
  abs(
    ms4a7_hr_time$HR_per_1SD -
      expected_hr_previous[
        as.character(ms4a7_hr_time$Time_years)
      ]
  )

qc_hr <- data.frame(
  Time_years = ms4a7_hr_time$Time_years,
  HR_current = ms4a7_hr_time$HR_per_1SD,
  HR_previous = expected_hr_previous[
    as.character(ms4a7_hr_time$Time_years)
  ],
  Absolute_difference = hr_abs_diff
)

cat("\n============================================================\n")
cat("Time-specific HR numerical reproducibility QC\n")
cat("============================================================\n\n")

print(
  qc_hr,
  row.names = FALSE
)

cat(
  "\nMaximum absolute difference = ",
  max(hr_abs_diff),
  "\n",
  sep = ""
)

stopifnot(
  max(hr_abs_diff) < 2e-4
)

cat(
  "\nTime-specific HR QC PASSED: differences are within ",
  "2e-4 numerical tolerance.\n",
  sep = ""
)

# ------------------------------------------------------------
# 23. Final QC summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 6B-6 QC PASSED.\n")
cat("============================================================\n\n")

cat("Locked cohort: N = 485; Events = 178\n")
cat("Gene: MS4A7\n")
cat("Time functions evaluated: linear, sqrt, log1p\n")
cat("Selected function by minimum AIC: sqrt\n")
cat(
  "Time-varying coefficient P = ",
  beta1_p,
  "\n",
  sep = ""
)

cat(
  "Constant vs time-varying LRT P = ",
  ms4a7_LRT_p,
  "\n",
  sep = ""
)

cat(
  "Delta AIC = ",
  ms4a7_AIC_delta,
  "\n",
  sep = ""
)

cat("\nTime-specific HRs:\n")

print(
  ms4a7_hr_time,
  row.names = FALSE
)

cat("\n============================================================\n")


# ============================================================
# Paper 1 — STEP 6B-7
# ============================================================
# STEP 6B-7 — Complete checkpoint save after path recovery
# ============================================================

# ------------------------------------------------------------
# 1. Re-define checkpoint path
# ------------------------------------------------------------

step6b7_checkpoint_file <-
  file.path(
    getwd(),
    "Paper1_Step6B6_MS4A7_time_varying_Cox_checkpoint.RData"
  )

cat("\n============================================================\n")
cat("STEP 6B-7: Saving final checkpoint\n")
cat("============================================================\n\n")

cat(
  "Checkpoint path:\n",
  step6b7_checkpoint_file,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 2. Confirm all formal checkpoint objects still exist
# ------------------------------------------------------------

missing_objects <-
  step6b7_objects[
    !vapply(
      step6b7_objects,
      exists,
      logical(1),
      inherits = TRUE
    )
  ]

if (length(missing_objects) > 0) {
  
  stop(
    paste(
      "Checkpoint objects are missing:",
      paste(missing_objects, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 3. Save
# ------------------------------------------------------------

save(
  list = step6b7_objects,
  file = step6b7_checkpoint_file
)


# ------------------------------------------------------------
# 4. Verify file creation
# ------------------------------------------------------------

stopifnot(
  file.exists(step6b7_checkpoint_file)
)


checkpoint_info <-
  file.info(step6b7_checkpoint_file)

checkpoint_size_MB <-
  checkpoint_info$size / 1024^2


# ------------------------------------------------------------
# 5. Final confirmation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 6B-7 CHECKPOINT SAVED SUCCESSFULLY.\n")
cat("============================================================\n\n")

cat(
  "File:\n",
  normalizePath(
    step6b7_checkpoint_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)

cat(
  "File size: ",
  round(checkpoint_size_MB, 3),
  " MB\n",
  sep = ""
)

cat(
  "Objects saved: ",
  length(step6b7_objects),
  "\n",
  sep = ""
)

cat(
  "Locked cohort: N = 485; Events = 178\n"
)

cat(
  "Gene: MS4A7\n"
)

cat(
  "Selected time function: sqrt(time)\n"
)

cat(
  "Time-varying coefficient P = ",
  beta1_p,
  "\n",
  sep = ""
)

cat(
  "Constant vs time-varying LRT P = ",
  ms4a7_LRT$P,
  "\n",
  sep = ""
)

cat(
  "Delta AIC = ",
  ms4a7_AIC_delta,
  "\n",
  sep = ""
)

cat(
  "\nSTEP 6B-7 FINAL CHECKPOINT: LOCKED.\n"
)

cat(
  "============================================================\n"
)

