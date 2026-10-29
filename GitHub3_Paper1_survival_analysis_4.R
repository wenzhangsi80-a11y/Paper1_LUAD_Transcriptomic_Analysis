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


# ============================================================
# STEP 6E — KM CURVE HARMONIZATION
# Phase 1: Load and inspect the locked survival cohort
# ============================================================


base_dir <- "input_dir/Paper 1_Step 6-survival analysis"

report_dir <- file.path(
  base_dir,
  "STEP_6D1_survival_reporting"
)

out_dir <- file.path(
  base_dir,
  "STEP_6E_KM_curve_harmonization"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

checkpoint <- file.path(
  report_dir,
  "Paper1_Step6D1_Survival_Reporting_checkpoint.RData"
)

if (!file.exists(checkpoint)) {
  stop("Checkpoint not found: ", checkpoint)
}

loaded_names <- load(checkpoint)

cat("\n========== LOADED OBJECTS ==========\n")
print(loaded_names)

# Inspect the expected downstream cohort object
if (!exists("survival_gene_primary_z_downstream")) {
  stop(
    "survival_gene_primary_z_downstream was not found. ",
    "Inspect the loaded objects and identify the saved 485-case cohort."
  )
}

dat <- as.data.frame(survival_gene_primary_z_downstream)

cat("\n========== DATA DIMENSIONS ==========\n")
print(dim(dat))

cat("\n========== COLUMN NAMES ==========\n")
print(names(dat))

cat("\n========== FIRST SIX ROWS ==========\n")
print(head(dat))

cat("\n========== DATA STRUCTURE ==========\n")
str(dat)

cat("\n========== MISSING VALUES ==========\n")
print(colSums(is.na(dat)))

write.csv(
  data.frame(
    column = names(dat),
    class = vapply(dat, function(x) paste(class(x), collapse = "/"),
                   character(1)),
    missing_n = colSums(is.na(dat))
  ),
  file.path(out_dir, "STEP6E_input_column_audit.csv"),
  row.names = FALSE
)

saveRDS(
  dat,
  file.path(out_dir, "STEP6E_input_cohort_checkpoint.rds")
)

cat("\nPhase 1 completed. Verify the column mapping before plotting.\n")



# ============================================================
# STEP 6E-1 — Inspect nested STEP 6D1 checkpoint
# ============================================================

cat("\n========== TOP-LEVEL OBJECTS ==========\n")
print(ls())

cat("\n========== CHECKPOINT CLASS ==========\n")
print(class(step6d1_checkpoint))

cat("\n========== CHECKPOINT STRUCTURE ==========\n")
str(step6d1_checkpoint, max.level = 2, give.attr = FALSE)

# Recursively inventory nested data frames and matrices
inventory_nested <- function(x, path = "root", depth = 0,
                             max_depth = 6) {
  
  out <- list()
  
  if (depth > max_depth) return(out)
  
  if (is.data.frame(x) || is.matrix(x)) {
    out[[length(out) + 1L]] <- data.frame(
      path = path,
      class = paste(class(x), collapse = "/"),
      nrow = nrow(x),
      ncol = ncol(x),
      columns = paste(colnames(x), collapse = " | "),
      stringsAsFactors = FALSE
    )
    return(out)
  }
  
  if (is.list(x)) {
    nm <- names(x)
    
    for (i in seq_along(x)) {
      child_name <- if (!is.null(nm) &&
                        nzchar(nm[i])) nm[i] else paste0("[[", i, "]]")
      
      child_path <- paste0(path, "$", child_name)
      
      child_result <- inventory_nested(
        x[[i]],
        path = child_path,
        depth = depth + 1L,
        max_depth = max_depth
      )
      
      out <- c(out, child_result)
    }
  }
  
  out
}

nested_inventory <- do.call(
  rbind,
  inventory_nested(step6d1_checkpoint)
)

cat("\n========== NESTED DATA-FRAME / MATRIX INVENTORY ==========\n")

if (is.null(nested_inventory) || nrow(nested_inventory) == 0) {
  cat("No nested data frames or matrices found within search depth.\n")
} else {
  print(nested_inventory, row.names = FALSE)
  
  write.csv(
    nested_inventory,
    file.path(
      "O:/Rproject/Paper 1_Step 6-survival analysis",
      "STEP6E_checkpoint_nested_inventory.csv"
    ),
    row.names = FALSE
  )
}

# Check whether the locked MS4A7 cohort is available
cat("\n========== MODEL A COHORT ==========\n")

if (exists("modelA_cc")) {
  cat("Dimensions:", paste(dim(modelA_cc), collapse = " x "), "\n")
  cat("Deaths:", sum(modelA_cc$OS_event == 1, na.rm = TRUE), "\n")
  cat("Censored:", sum(modelA_cc$OS_event == 0, na.rm = TRUE), "\n")
  cat("Unique IDs:", length(unique(modelA_cc$Sample)), "\n")
  print(head(modelA_cc))
} else {
  cat("modelA_cc is not currently loaded as a top-level object.\n")
}


# ============================================================
# STEP 6E-2 — Locate the locked survival cohort and expression data
# Read-only diagnostic; does not recalculate survival models
# ============================================================

cat("\n========== A. INSPECT THE OTHER CHECKPOINT OBJECT ==========\n")

if (exists("checkpoint")) {
  cat("Class:\n")
  print(class(checkpoint))
  
  cat("\nNames:\n")
  print(names(checkpoint))
  
  cat("\nStructure (maximum depth 2):\n")
  str(checkpoint, max.level = 2, give.attr = FALSE)
} else {
  cat("Object 'checkpoint' is absent.\n")
}

cat("\n========== B. SEARCH SURVIVAL ANALYSIS FILES ==========\n")

search_dirs <- unique(c(
  "O:/Rproject/Paper 1_Step 6-survival analysis",
  "O:/Rproject/Paper 1_Step 6-survival analysis/STEP_6D1_survival_reporting"
))

for (d in search_dirs) {
  cat("\nDirectory:", d, "\n")
  
  if (dir.exists(d)) {
    ff <- list.files(
      d,
      recursive = TRUE,
      full.names = TRUE,
      all.files = FALSE
    )
    
    candidate <- ff[
      grepl(
        "\\.(rds|rda|rdata|csv|tsv|xlsx|r|rmd|qs)$",
        ff,
        ignore.case = TRUE
      )
    ]
    
    cat("Candidate file count:", length(candidate), "\n")
    
    if (length(candidate) > 0) {
      print(data.frame(
        file = candidate,
        size_MB = round(
          file.info(candidate)$size / 1024^2, 3
        ),
        modified = as.character(
          file.info(candidate)$mtime
        ),
        row.names = NULL
      ))
    }
  } else {
    cat("Directory not found.\n")
  }
}

cat("\n========== C. SEARCH FILE NAMES FOR COHORT / EXPRESSION TERMS ==========\n")

root_dir <- "O:/Rproject/Paper 1_Step 6-survival analysis"

if (dir.exists(root_dir)) {
  all_files <- list.files(
    root_dir,
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE
  )
  
  relevant_files <- all_files[
    grepl(
      "survival_gene_primary_z|modelA_cc|primary_z|survival.*cohort|cohort.*survival|expression.*survival|survival.*expression|step6c5|step6d",
      basename(all_files),
      ignore.case = TRUE
    )
  ]
  
  print(relevant_files)
}

cat("\n========== D. REPORT CHECKPOINT REPORTING TABLES ==========\n")

for (nm in c(
  "univariable_results",
  "modelA_results",
  "modelB_results",
  "continuous_cox_table",
  "figure4D_data"
)) {
  if (!is.null(step6d1_checkpoint[[nm]])) {
    cat("\n---", nm, "---\n")
    print(step6d1_checkpoint[[nm]])
  }
}

cat("\n========== STEP 6E-2 COMPLETE ==========\n")



# ============================================================
# STEP 6E-3 — Trace KM inputs and inspect smaller checkpoints
# Read-only diagnostic
# ============================================================

root <- "input_dir/Paper 1_Step 6-survival analysis"

script_file <- file.path(
  root,
  "Paper 1 – STEP_6E_KM_curve_harmonization.R"
)

cat("\n========== A. INSPECT STEP 6E SCRIPT ==========\n")

if (file.exists(script_file)) {
  
  script_lines <- readLines(
    script_file,
    warn = FALSE,
    encoding = "UTF-8"
  )
  
  # Show lines likely to identify input data, checkpoint objects,
  # cohort filtering, expression scaling, or KM construction.
  pattern <- paste(
    "load\\(",
    "readRDS\\(",
    "read\\.csv\\(",
    "read\\.table\\(",
    "read_excel\\(",
    "read\\.xlsx\\(",
    "saveRDS\\(",
    "save\\(",
    "survfit\\(",
    "Surv\\(",
    "survival_gene_primary_z",
    "modelA_cc",
    "OS_time",
    "OS_event",
    "AOC3",
    "MS4A7",
    "IL7R",
    "median",
    "checkpoint",
    sep = "|"
  )
  
  hits <- grep(
    pattern,
    script_lines,
    ignore.case = TRUE
  )
  
  cat("Script lines:", length(script_lines), "\n")
  cat("Matching lines:", length(hits), "\n\n")
  
  if (length(hits) > 0) {
    # Print matching lines with line numbers; cap output
    # to avoid flooding the console.
    hits_to_show <- head(hits, 160)
    
    cat(paste(
      sprintf("%4d: %s", hits_to_show, script_lines[hits_to_show]),
      collapse = "\n"
    ), "\n")
  }
  
} else {
  cat("STEP 6E script not found at the expected path.\n")
}

cat("\n========== B. INSPECT SMALLER CHECKPOINTS ==========\n")

checkpoint_files <- c(
  step6B5 = file.path(
    root,
    "Paper1_Step6B5_continuous_Cox_PH_checkpoint.RData"
  ),
  step6B6 = file.path(
    root,
    "Paper1_Step6B6_MS4A7_time_varying_Cox_checkpoint.RData"
  ),
  step6C5 = file.path(
    root,
    "Paper1_Step6C5_adjusted_MS4A7_time_varying_Cox_checkpoint.RData"
  ),
  step6D1_continuous = file.path(
    root,
    "STEP_6D1_survival_reporting",
    "Paper1_Step6D1_continuous_Cox_reporting_checkpoint.RData"
  )
)

inventory_objects <- function(x, prefix = "") {
  
  result <- list()
  
  for (nm in names(x)) {
    
    obj <- x[[nm]]
    object_path <- if (nzchar(prefix)) {
      paste0(prefix, "$", nm)
    } else {
      nm
    }
    
    if (is.data.frame(obj) || is.matrix(obj)) {
      
      result[[length(result) + 1L]] <- data.frame(
        object = object_path,
        class = paste(class(obj), collapse = "/"),
        nrow = nrow(obj),
        ncol = ncol(obj),
        columns = paste(colnames(obj), collapse = " | "),
        stringsAsFactors = FALSE
      )
      
    } else if (is.list(obj) && length(obj) > 0) {
      
      nested <- inventory_objects(obj, object_path)
      result <- c(result, nested)
    }
  }
  
  result
}

for (label in names(checkpoint_files)) {
  
  f <- checkpoint_files[[label]]
  
  cat("\n---", label, "---\n")
  cat("File:", f, "\n")
  
  if (!file.exists(f)) {
    cat("File not found.\n")
    next
  }
  
  # Load into a separate environment, not the global workspace.
  e <- new.env(parent = emptyenv())
  
  loaded <- tryCatch(
    load(f, envir = e),
    error = function(err) {
      cat("LOAD ERROR:", conditionMessage(err), "\n")
      character(0)
    }
  )
  
  if (length(loaded) == 0) next
  
  cat("Top-level objects:\n")
  print(loaded)
  
  inv <- inventory_objects(e)
  
  if (length(inv) > 0) {
    inv_df <- do.call(rbind, inv)
    print(inv_df, row.names = FALSE)
  } else {
    cat("No nested data frames or matrices found.\n")
  }
  
  # Show useful object names without printing large expression matrices.
  candidate_names <- loaded[
    grepl(
      "survival_gene_primary_z|modelA_cc|primary_z|cohort|clinical|expression|survival",
      loaded,
      ignore.case = TRUE
    )
  ]
  
  if (length(candidate_names) > 0) {
    cat("Potentially relevant objects:\n")
    print(candidate_names)
  }
  
  rm(e)
}

cat("\n========== STEP 6E-3 COMPLETE ==========\n")



# ============================================================
# STEP 6E-4 — Validate the locked three-gene survival dataset
# Read-only; no re-standardization and no model refitting
# ============================================================

f <- file.path(
  "input_dir/Paper 1_Step 6-survival analysis",
  "Paper1_Step6B5_continuous_Cox_PH_checkpoint.RData"
)

stopifnot(file.exists(f))

e <- new.env(parent = emptyenv())
loaded <- load(f, envir = e)

cat("\n========== A. SOURCE OBJECTS ==========\n")
print(loaded)

required <- c("survival_gene_primary_z", "z_parameters")

cat("\nRequired objects available:\n")
print(setNames(required %in% loaded, required))

stopifnot(all(required %in% loaded))

dat <- e$survival_gene_primary_z
zp  <- e$z_parameters

cat("\n========== B. DATA STRUCTURE ==========\n")
print(dim(dat))
print(names(dat))
print(head(dat, 3))

cat("\n========== C. SAMPLE / EVENT AUDIT ==========\n")

cat("Rows:", nrow(dat), "\n")
cat("Unique sample IDs:", length(unique(dat$Sample)), "\n")
cat("Duplicated sample IDs:", sum(duplicated(dat$Sample)), "\n")
cat("Missing sample IDs:", sum(is.na(dat$Sample)), "\n")
cat("Missing OS.time:", sum(is.na(dat$OS.time)), "\n")
cat("Missing OS.status:", sum(is.na(dat$OS.status)), "\n")
cat("Missing OS_event:", sum(is.na(dat$OS_event)), "\n")

cat("\nOS_event counts:\n")
print(table(dat$OS_event, useNA = "ifany"))

cat("\nOS.status counts:\n")
print(table(dat$OS.status, useNA = "ifany"))

cat("\nEvent/status cross-tabulation:\n")
print(table(
  OS_status = dat$OS.status,
  OS_event = dat$OS_event,
  useNA = "ifany"
))

cat("\nNon-positive survival times:\n")
print(sum(dat$OS.time <= 0, na.rm = TRUE))

cat("\n========== D. THREE-GENE EXPRESSION AUDIT ==========\n")

genes <- c("AOC3", "MS4A7", "IL7R")

expr_summary <- do.call(rbind, lapply(genes, function(g) {
  x <- dat[[g]]
  
  data.frame(
    Gene = g,
    N_nonmissing = sum(!is.na(x)),
    Mean = mean(x, na.rm = TRUE),
    SD = sd(x, na.rm = TRUE),
    Median = median(x, na.rm = TRUE),
    Min = min(x, na.rm = TRUE),
    Max = max(x, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

print(expr_summary, row.names = FALSE)

cat("\nSaved z-score parameters:\n")
print(zp)

cat("\n========== E. SAMPLE ID / ROW-NAME CHECK ==========\n")
cat("Has row names identical to Sample column:",
    identical(rownames(dat), dat$Sample), "\n")
cat("First 10 sample IDs:\n")
print(head(dat$Sample, 10))

cat("\n========== F. EXPORT AUDIT SUMMARY ==========\n")

out_dir <- file.path(
  "result_dir/Paper 1_Step 6-survival analysis",
  "STEP_6E_KM_harmonization_audit"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(
  expr_summary,
  file.path(out_dir, "STEP6E4_expression_summary.csv"),
  row.names = FALSE
)

write.csv(
  zp,
  file.path(out_dir, "STEP6E4_saved_z_parameters.csv"),
  row.names = FALSE
)

cat("Audit summaries saved to:\n", out_dir, "\n")
cat("\nSTEP 6E-4 COMPLETE\n")


# ============================================================
# STEP 6E-5B — Zero-time cohort impact audit
# Uses the already-locked STEP 6B5 survival cohort
# Read-only: does not alter or overwrite the source data
# ============================================================

source_file <- file.path(
  "input_dir/Paper 1_Step 6-survival analysis",
  "Paper1_Step6B5_continuous_Cox_PH_checkpoint.RData"
)

e <- new.env(parent = emptyenv())
load(source_file, envir = e)

dat <- e$survival_gene_primary_z

stopifnot(
  nrow(dat) == 485,
  length(unique(dat$Sample)) == 485,
  sum(dat$OS_event == 1) == 178,
  sum(dat$OS_event == 0) == 307
)

zero_dat <- dat[dat$OS.time == 0, , drop = FALSE]
positive_dat <- dat[dat$OS.time > 0, , drop = FALSE]

cat("\n========== A. LOCKED COHORT ==========\n")
cat("N =", nrow(dat), "\n")
cat("Events =", sum(dat$OS_event == 1), "\n")
cat("Censored =", sum(dat$OS_event == 0), "\n")

cat("\n========== B. ZERO-TIME OBSERVATIONS ==========\n")
print(
  zero_dat[, c("Sample", "OS.time", "OS.status", "OS_event")],
  row.names = FALSE
)

cat("\nZero-time N =", nrow(zero_dat), "\n")
cat("Zero-time events =", sum(zero_dat$OS_event == 1), "\n")
cat("Zero-time censored =", sum(zero_dat$OS_event == 0), "\n")

cat("\n========== C. POSITIVE-TIME SUBSET (AUDIT ONLY) ==========\n")
cat("N =", nrow(positive_dat), "\n")
cat("Events =", sum(positive_dat$OS_event == 1), "\n")
cat("Censored =", sum(positive_dat$OS_event == 0), "\n")

cat("\n========== D. ZERO-TIME SAMPLE IDS ==========\n")
print(zero_dat$Sample)

cat("\nSTEP 6E-5B COMPLETE — source data unchanged.\n")



script_file <- file.path(
  "input_dir/Paper 1_Step 6-survival analysis",
  "Paper 1 – STEP_6E_KM_curve_harmonization.R"
)

script_lines <- readLines(
  script_file,
  warn = FALSE,
  encoding = "UTF-8"
)

cat(
  script_lines[seq_len(min(85, length(script_lines)))],
  sep = "\n"
)


# ============================================================
# STEP 6E — KM CURVE HARMONIZATION
# Phase 1: Load and audit the locked 485-case survival cohort
#
# Authoritative source:
#   Paper1_Step6B5_continuous_Cox_PH_checkpoint.RData
#   Object: survival_gene_primary_z
#
# This phase does NOT:
#   - modify OS.time
#   - remove zero-time observations
#   - re-standardize gene expression
#   - calculate or plot KM curves
# ============================================================

# ------------------------------------------------------------
# 0. Paths and output directory
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 6-survival analysis"

source_checkpoint <- file.path(
  base_dir,
  "Paper1_Step6B5_continuous_Cox_PH_checkpoint.RData"
)

out_dir <- file.path(
  base_dir,
  "STEP_6E_KM_curve_harmonization"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(source_checkpoint)) {
  stop("Authoritative checkpoint not found: ", source_checkpoint)
}

cat("\n============================================================\n")
cat("STEP 6E — PHASE 1: LOCKED COHORT AUDIT\n")
cat("============================================================\n")

cat("\nSource checkpoint:\n", source_checkpoint, "\n")
cat("Output directory:\n", out_dir, "\n")

# ------------------------------------------------------------
# 1. Load checkpoint into an isolated environment
# ------------------------------------------------------------

source_env <- new.env(parent = emptyenv())

loaded_names <- load(
  source_checkpoint,
  envir = source_env
)

cat("\n========== A. LOADED OBJECTS ==========\n")
print(loaded_names)

required_object <- "survival_gene_primary_z"

if (!required_object %in% loaded_names) {
  stop(
    "Required object not found in checkpoint: ",
    required_object
  )
}

dat <- as.data.frame(
  source_env$survival_gene_primary_z
)

# ------------------------------------------------------------
# 2. Required columns and data types
# ------------------------------------------------------------

required_columns <- c(
  "Sample",
  "OS.time",
  "OS.status",
  "AOC3",
  "MS4A7",
  "IL7R",
  "OS_event"
)

missing_columns <- setdiff(
  required_columns,
  names(dat)
)

if (length(missing_columns) > 0) {
  stop(
    "Required columns missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (!is.numeric(dat$OS.time)) {
  stop("OS.time is not numeric. Do not coerce silently.")
}

if (!is.numeric(dat$OS_event)) {
  stop("OS_event is not numeric. Inspect source coding first.")
}

gene_columns <- c("AOC3", "MS4A7", "IL7R")

if (!all(vapply(dat[gene_columns], is.numeric, logical(1)))) {
  stop("At least one gene-expression column is not numeric.")
}

cat("\n========== B. DATA DIMENSIONS ==========\n")
print(dim(dat))

cat("\nColumn names:\n")
print(names(dat))

# ------------------------------------------------------------
# 3. Locked cohort integrity checks
# ------------------------------------------------------------

if (nrow(dat) != 485L) {
  stop("Expected 485 cases; observed ", nrow(dat))
}

if (anyNA(dat$Sample) || any(dat$Sample == "")) {
  stop("Missing or empty sample IDs detected.")
}

if (anyDuplicated(dat$Sample) > 0L) {
  stop("Duplicate sample IDs detected.")
}

if (anyNA(dat[, c(
  "OS.time", "OS.status", "OS_event",
  "AOC3", "MS4A7", "IL7R"
)])) {
  stop(
    "Missing values detected in required survival/gene fields. ",
    "Do not silently drop cases."
  )
}

if (any(dat$OS.time < 0)) {
  stop("Negative survival times detected. Inspect source data.")
}

if (!all(dat$OS_event %in% c(0, 1))) {
  stop("OS_event contains values other than 0 and 1.")
}

if (!all(dat$OS.status %in% c("Alive", "Dead"))) {
  stop("Unexpected OS.status labels detected.")
}

if (!all(
  (dat$OS.status == "Dead") == (dat$OS_event == 1)
)) {
  stop("OS.status and OS_event are inconsistent.")
}

if (sum(dat$OS_event == 1) != 178L ||
    sum(dat$OS_event == 0) != 307L) {
  stop("Event/censor counts do not match the locked cohort.")
}

cat("\n========== C. LOCKED COHORT SUMMARY ==========\n")
cat("N =", nrow(dat), "\n")
cat("Unique sample IDs =", length(unique(dat$Sample)), "\n")
cat("Events =", sum(dat$OS_event == 1), "\n")
cat("Censored =", sum(dat$OS_event == 0), "\n")
cat("Missing values in required fields = 0\n")
cat("Negative survival times =", sum(dat$OS.time < 0), "\n")

cat("\nOS.status counts:\n")
print(table(dat$OS.status))

cat("\nOS_event counts:\n")
print(table(dat$OS_event))

# ------------------------------------------------------------
# 4. Zero-time audit — retain all original observations
# ------------------------------------------------------------

zero_dat <- dat[
  dat$OS.time == 0,
  c(
    "Sample", "OS.time", "OS.status", "OS_event",
    "AOC3", "MS4A7", "IL7R"
  ),
  drop = FALSE
]

cat("\n========== D. ZERO-TIME AUDIT ==========\n")
cat("Zero-time observations =", nrow(zero_dat), "\n")
cat("Zero-time events =", sum(zero_dat$OS_event == 1), "\n")
cat("Zero-time censored =", sum(zero_dat$OS_event == 0), "\n")

print(zero_dat, row.names = FALSE)

if (nrow(zero_dat) != 4L ||
    sum(zero_dat$OS_event == 1) != 1L ||
    sum(zero_dat$OS_event == 0) != 3L) {
  stop(
    "Zero-time counts differ from the previously locked audit. ",
    "Review provenance before proceeding."
  )
}

# ------------------------------------------------------------
# 5. Verify that gene values are already standardized
# ------------------------------------------------------------

cat("\n========== E. GENE EXPRESSION AUDIT ==========\n")

gene_summary <- do.call(
  rbind,
  lapply(gene_columns, function(g) {
    x <- dat[[g]]
    
    data.frame(
      Gene = g,
      N = sum(!is.na(x)),
      Mean = mean(x),
      SD = sd(x),
      Median = median(x),
      Min = min(x),
      Max = max(x),
      stringsAsFactors = FALSE
    )
  })
)

print(gene_summary, row.names = FALSE)

# These checks verify the expected Z-score scale.
# They do not transform or re-standardize the data.
if (any(abs(gene_summary$Mean) > 1e-6) ||
    any(abs(gene_summary$SD - 1) > 1e-6)) {
  stop(
    "Gene values do not match the expected standardized scale. ",
    "Do not re-standardize automatically."
  )
}

# ------------------------------------------------------------
# 6. Export audit tables and unchanged input cohort
# ------------------------------------------------------------

write.csv(
  gene_summary,
  file.path(out_dir, "STEP6E_phase1_gene_expression_audit.csv"),
  row.names = FALSE
)

write.csv(
  zero_dat,
  file.path(out_dir, "STEP6E_phase1_zero_time_audit.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    Sample = dat$Sample,
    OS.time = dat$OS.time,
    OS.status = dat$OS.status,
    OS_event = dat$OS_event,
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "STEP6E_phase1_survival_fields_audit.csv"),
  row.names = FALSE
)

# Save the unchanged source cohort for downstream use.
saveRDS(
  dat,
  file.path(out_dir, "STEP6E_input_cohort_checkpoint.rds")
)

# Record source provenance and the handling rule.
provenance_lines <- c(
  "STEP 6E Phase 1 — input provenance",
  paste("Source checkpoint:", source_checkpoint),
  paste("Source object:", required_object),
  paste("Rows:", nrow(dat)),
  paste("Events:", sum(dat$OS_event == 1)),
  paste("Censored:", sum(dat$OS_event == 0)),
  paste("Zero-time observations:", nrow(zero_dat)),
  "OS.time was not modified.",
  "No observations were removed.",
  "Gene expression was not re-standardized.",
  "KM curves have not yet been calculated."
)

writeLines(
  provenance_lines,
  file.path(out_dir, "STEP6E_phase1_provenance.txt")
)

cat("\n========== F. OUTPUT FILES ==========\n")
print(list.files(out_dir, full.names = TRUE))

cat("\nSTEP 6E PHASE 1 COMPLETE.\n")
cat("Input cohort retained: N = 485; events = 178; censored = 307.\n")
cat("Zero-time observations retained: N = 4.\n")
cat("No KM curves generated in this phase.\n")




# ============================================================
# STEP 6E — Phase 2A
# Audit median-split rules before KM estimation
# ============================================================

base_dir <- "input_dir/Paper 1_Step 6-survival analysis"

out_dir <- file.path(
  base_dir,
  "STEP_6E_KM_curve_harmonization"
)

input_rds <- file.path(
  out_dir,
  "STEP6E_input_cohort_checkpoint.rds"
)

stopifnot(file.exists(input_rds))

dat <- readRDS(input_rds)

genes <- c("AOC3", "MS4A7", "IL7R")

stopifnot(
  nrow(dat) == 485,
  sum(dat$OS_event == 1) == 178,
  sum(dat$OS_event == 0) == 307,
  all(dat$OS.time >= 0)
)

cat("\n========== A. R AND PACKAGE VERSIONS ==========\n")
cat("R:", R.version.string, "\n")

if (requireNamespace("survival", quietly = TRUE)) {
  cat(
    "survival:",
    as.character(packageVersion("survival")),
    "\n"
  )
} else {
  stop("Package 'survival' is not installed.")
}

if (requireNamespace("survminer", quietly = TRUE)) {
  cat(
    "survminer:",
    as.character(packageVersion("survminer")),
    "\n"
  )
} else {
  cat("survminer is not installed; this audit does not require it.\n")
}

cat("\n========== B. MEDIAN-SPLIT AUDIT ==========\n")

split_summary <- do.call(
  rbind,
  lapply(genes, function(g) {
    
    x <- dat[[g]]
    med <- median(x)
    
    # Prespecified deterministic rule:
    # values <= median -> Low; values > median -> High
    group <- ifelse(x > med, "High", "Low")
    
    data.frame(
      Gene = g,
      Median = med,
      N_Low = sum(group == "Low"),
      N_High = sum(group == "High"),
      N_Equal_Median = sum(x == med),
      stringsAsFactors = FALSE
    )
  })
)

print(split_summary, row.names = FALSE)

cat("\n========== C. ZERO-TIME OBSERVATIONS BY GROUP ==========\n")

zero_time_split <- do.call(
  rbind,
  lapply(genes, function(g) {
    
    med <- median(dat[[g]])
    group <- ifelse(dat[[g]] > med, "High", "Low")
    
    z <- dat$OS.time == 0
    
    data.frame(
      Gene = g,
      Group = c("Low", "High"),
      N_zero_time = c(
        sum(z & group == "Low"),
        sum(z & group == "High")
      ),
      N_zero_time_events = c(
        sum(z & group == "Low" & dat$OS_event == 1),
        sum(z & group == "High" & dat$OS_event == 1)
      ),
      stringsAsFactors = FALSE
    )
  })
)

print(zero_time_split, row.names = FALSE)

cat("\n========== D. CHECK KM ESTIMATION SUPPORT ==========\n")

# The survival package is used for the later KM step.
# Do not fit the curves in this audit stage.
cat("survival::Surv available:",
    exists("Surv", envir = asNamespace("survival")),
    "\n")
cat("survival::survfit available:",
    exists("survfit", envir = asNamespace("survival")),
    "\n")
cat("Zero-time observations retained:",
    sum(dat$OS.time == 0),
    "\n")

write.csv(
  split_summary,
  file.path(out_dir, "STEP6E_phase2A_median_split_audit.csv"),
  row.names = FALSE
)

write.csv(
  zero_time_split,
  file.path(out_dir, "STEP6E_phase2A_zero_time_by_group.csv"),
  row.names = FALSE
)

cat("\nSTEP 6E PHASE 2A COMPLETE.\n")
cat("No KM models fitted and no survival P values calculated.\n")



# ============================================================
# STEP 6E — PHASE 2B
# Fit and export harmonized Kaplan-Meier curves
#
# Locked input:
#   STEP6E_input_cohort_checkpoint.rds
#   N = 485; events = 178; censored = 307
#
# Rules:
#   - median split within the locked 485-case cohort
#   - expression <= median: Low
#   - expression > median: High
#   - retain all four OS.time = 0 observations
#   - do not re-standardize gene expression
# ============================================================

# ------------------------------------------------------------
# 0. Paths and packages
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 6-survival analysis"

out_dir <- file.path(
  base_dir,
  "STEP_6E_KM_curve_harmonization"
)

input_rds <- file.path(
  out_dir,
  "STEP6E_input_cohort_checkpoint.rds"
)

plot_dir <- file.path(out_dir, "KM_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_rds)) {
  stop("Phase 1 input RDS not found: ", input_rds)
}

required_packages <- c("survival", "survminer", "ggplot2")

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
    "Required packages not installed: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(survival)
library(survminer)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load and validate the locked cohort
# ------------------------------------------------------------

dat <- readRDS(input_rds)

genes <- c("AOC3", "MS4A7", "IL7R")

stopifnot(
  nrow(dat) == 485L,
  length(unique(dat$Sample)) == 485L,
  sum(dat$OS_event == 1) == 178L,
  sum(dat$OS_event == 0) == 307L,
  sum(dat$OS.time == 0) == 4L,
  all(dat$OS.time >= 0),
  all(dat$OS_event %in% c(0, 1)),
  all(genes %in% names(dat))
)

cat("\n========== A. SOFTWARE VERSIONS ==========\n")
cat("R:", R.version.string, "\n")
cat("survival:", as.character(packageVersion("survival")), "\n")
cat("survminer:", as.character(packageVersion("survminer")), "\n")
cat("ggplot2:", as.character(packageVersion("ggplot2")), "\n")

cat("\n========== B. LOCKED COHORT ==========\n")
cat("N =", nrow(dat), "\n")
cat("Events =", sum(dat$OS_event == 1), "\n")
cat("Censored =", sum(dat$OS_event == 0), "\n")
cat("Zero-time observations =", sum(dat$OS.time == 0), "\n")

# ------------------------------------------------------------
# 2. Fixed median split
# ------------------------------------------------------------

split_summary <- list()
group_membership <- list()
km_fits <- list()
km_plots <- list()
logrank_results <- list()

for (g in genes) {
  
  cat("\n============================================================\n")
  cat("GENE:", g, "\n")
  cat("============================================================\n")
  
  med <- median(dat[[g]])
  
  gene_dat <- dat
  gene_dat$Gene <- g
  gene_dat$Expression <- gene_dat[[g]]
  
  # Fixed rule: equal to median belongs to Low.
  gene_dat$Group <- factor(
    ifelse(gene_dat$Expression > med, "High", "Low"),
    levels = c("Low", "High")
  )
  
  # Confirm no cases have been dropped.
  stopifnot(
    nrow(gene_dat) == 485L,
    sum(gene_dat$Group == "Low") == 243L,
    sum(gene_dat$Group == "High") == 242L
  )
  
  cat("Median threshold:", format(med, digits = 10), "\n")
  cat("Group counts:\n")
  print(table(gene_dat$Group))
  
  cat("Events by group:\n")
  print(table(
    Group = gene_dat$Group,
    Event = gene_dat$OS_event
  ))
  
  cat("Zero-time observations by group:\n")
  print(table(
    Group = gene_dat$Group[gene_dat$OS.time == 0]
  ))
  
  # Save exact group assignment for reproducibility.
  group_membership[[g]] <- data.frame(
    Sample = gene_dat$Sample,
    Gene = g,
    Expression = gene_dat$Expression,
    Median_threshold = med,
    Group = as.character(gene_dat$Group),
    OS.time = gene_dat$OS.time,
    OS.status = gene_dat$OS.status,
    OS_event = gene_dat$OS_event,
    stringsAsFactors = FALSE
  )
  
  # ----------------------------------------------------------
  # 3. Kaplan-Meier fit
  # ----------------------------------------------------------
  fit <- survival::survfit(
    survival::Surv(OS.time, OS_event) ~ Group,
    data = gene_dat,
    conf.type = "log"
  )
  
  lr <- survival::survdiff(
    survival::Surv(OS.time, OS_event) ~ Group,
    data = gene_dat,
    rho = 0
  )
  
  # ----------------------------------------------------------
  # 4. Log-rank test
  # ----------------------------------------------------------
  
  lr <- survival::survdiff(
    km_formula,
    data = gene_dat,
    rho = 0
  )
  
  lr_df <- length(lr$n) - 1L
  
  lr_p <- stats::pchisq(
    lr$chisq,
    df = lr_df,
    lower.tail = FALSE
  )
  
  logrank_results[[g]] <- data.frame(
    Gene = g,
    Median_threshold = med,
    N_total = nrow(gene_dat),
    N_Low = sum(gene_dat$Group == "Low"),
    N_High = sum(gene_dat$Group == "High"),
    Events_total = sum(gene_dat$OS_event == 1),
    Events_Low = sum(
      gene_dat$OS_event == 1 & gene_dat$Group == "Low"
    ),
    Events_High = sum(
      gene_dat$OS_event == 1 & gene_dat$Group == "High"
    ),
    Zero_time_total = sum(gene_dat$OS.time == 0),
    Zero_time_events = sum(
      gene_dat$OS.time == 0 & gene_dat$OS_event == 1
    ),
    Logrank_chisq = unname(lr$chisq),
    Logrank_df = lr_df,
    Logrank_P = lr_p,
    stringsAsFactors = FALSE
  )
  
  cat("\nLog-rank test:\n")
  print(logrank_results[[g]], row.names = FALSE)
  
  # ----------------------------------------------------------
  # 5. Plot with confidence intervals and risk table
  # ----------------------------------------------------------
  
  p <- survminer::ggsurvplot(
    fit,
    data = gene_dat,
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.title = "Number at risk",
    conf.int = TRUE,
    pval = FALSE,
    censor = TRUE,
    censor.shape = 124,
    censor.size = 3,
    xlab = "Overall survival time (days)",
    ylab = "Overall survival probability",
    title = paste0(g, " expression and overall survival"),
    legend.title = "Expression",
    legend.labs = c("Low", "High"),
    palette = c("#3C5488", "#E64B35"),
    ggtheme = ggplot2::theme_classic(
      base_size = 12,
      base_family = "Arial"
    ),
    tables.theme = ggplot2::theme_classic(
      base_size = 10,
      base_family = "Arial"
    )
  )
  
  p$plot <- p$plot +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      legend.title = ggplot2::element_text(face = "bold")
    )
  
  p$table <- p$table +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black")
    )
  
  km_plots[[g]] <- p
  
  # ----------------------------------------------------------
  # 6. Export combined curve + risk table
  # ----------------------------------------------------------
  
  tiff_file <- file.path(
    plot_dir,
    paste0("STEP6E_KM_", g, "_485case.tiff")
  )
  
  grDevices::tiff(
    filename = tiff_file,
    width = 7.5,
    height = 7.0,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  
  print(p)
  grDevices::dev.off()
  
  pdf_file <- file.path(
    plot_dir,
    paste0("STEP6E_KM_", g, "_485case.pdf")
  )
  
  grDevices::cairo_pdf(
    filename = pdf_file,
    width = 7.5,
    height = 7.0,
    family = "Arial"
  )
  
  print(p)
  grDevices::dev.off()
  
  cat("\nSaved:\n")
  cat(tiff_file, "\n")
  cat(pdf_file, "\n")
}

# Refit AOC3 using an explicit formula
fit_test <- survival::survfit(
  survival::Surv(OS.time, OS_event) ~ Group,
  data = gene_dat,
  conf.type = "log"
)

print(fit_test)
print(fit_test$call)

# Test plotting separately
p_test <- survminer::ggsurvplot(
  fit_test,
  data = gene_dat,
  risk.table = TRUE,
  conf.int = TRUE,
  pval = FALSE,
  censor = TRUE,
  xlab = "Overall survival time (days)",
  ylab = "Overall survival probability",
  legend.title = "Expression",
  legend.labs = c("Low", "High")
)

print(p_test)



# ------------------------------------------------------------
# 7. Export tables and model objects
# ------------------------------------------------------------

split_table <- do.call(rbind, group_membership)
logrank_table <- do.call(rbind, logrank_results)

write.csv(
  split_table,
  file.path(out_dir, "STEP6E_phase2B_locked_group_membership.csv"),
  row.names = FALSE
)

write.csv(
  logrank_table,
  file.path(out_dir, "STEP6E_phase2B_logrank_results.csv"),
  row.names = FALSE
)

saveRDS(
  km_fits,
  file.path(out_dir, "STEP6E_phase2B_KM_fit_objects.rds")
)

saveRDS(
  km_plots,
  file.path(out_dir, "STEP6E_phase2B_KM_plot_objects.rds")
)

writeLines(
  c(
    "STEP 6E Phase 2B — KM analysis record",
    paste("Input:", input_rds),
    paste("R:", R.version.string),
    paste("survival:", as.character(packageVersion("survival"))),
    paste("survminer:", as.character(packageVersion("survminer"))),
    paste("ggplot2:", as.character(packageVersion("ggplot2"))),
    "Cohort: 485 cases; 178 events; 307 censored.",
    "Median split: expression <= median = Low; > median = High.",
    "Zero-time observations retained without modification.",
    "KM: Kaplan-Meier estimate; confidence intervals displayed.",
    "Group comparison: two-sided log-rank test.",
    "Continuous Cox results were not recalculated or changed."
  ),
  file.path(out_dir, "STEP6E_phase2B_analysis_record.txt")
)

cat("\n========== PHASE 2B COMPLETE ==========\n")
cat("KM plots saved to:\n", plot_dir, "\n")

cat("\nLog-rank summary:\n")
print(logrank_table, row.names = FALSE)

cat("\nPlease inspect all three TIFF/PDF plots before final figure assembly.\n")


# ============================================================
# STEP 6E Phase 2B — Output completeness audit
# ============================================================

cat("\nRows in logrank_table:\n")
print(nrow(logrank_table))

cat("\nGenes in logrank_table:\n")
print(logrank_table$Gene)

cat("\nRows in group_membership:\n")
print(nrow(split_table))

cat("\nGenes in group_membership:\n")
print(unique(split_table$Gene))

cat("\nNumber of fitted KM objects:\n")
print(length(km_fits))
print(names(km_fits))

cat("\nNumber of plot objects:\n")
print(length(km_plots))
print(names(km_plots))

cat("\nFiles currently in KM plot directory:\n")
print(
  list.files(
    plot_dir,
    pattern = "\\.(tif|tiff|pdf)$",
    ignore.case = TRUE,
    full.names = TRUE
  )
)



# ============================================================
# STEP 6E Phase 2B — Pre-rerun environment check
# ============================================================

required_objects <- c(
  "input_rds",
  "out_dir",
  "plot_dir",
  "gene_list",
  "km_fits",
  "km_plots",
  "group_membership",
  "logrank_results"
)

object_status <- data.frame(
  object = required_objects,
  exists = vapply(
    required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
)

print(object_status)

cat("\nGene list:\n")
if (exists("gene_list")) print(gene_list)

cat("\nOutput directory:\n")
if (exists("out_dir")) print(out_dir)

cat("\nPlot directory:\n")
if (exists("plot_dir")) print(plot_dir)


# ============================================================
# STEP 6E — PHASE 2B + PHASE 2C
# Complete KM rerun + grouped Cox HR annotation
#
# Locked cohort: N = 485
# Fixed split: Expression <= median = Low; > median = High
#
# Preserve:
#   - original cohort and median-split rule
#   - KM confidence intervals
#   - number-at-risk table
#   - original curve colours
#   - original time unit: days
#
# Modify only:
#   - add log-rank P and grouped Cox HR (95% CI)
#   - annotation at lower-left, without border
#   - KM axis font: Arial, 10 pt, bold
# ============================================================


# ------------------------------------------------------------
# 0. Packages and input validation
# ------------------------------------------------------------

library(survival)
library(survminer)
library(ggplot2)
library(grid)

if (!exists("dat")) {
  stop(
    "Object 'dat' not found. Load STEP6E_input_cohort_checkpoint.rds first."
  )
}

gene_list <- c("AOC3", "MS4A7", "IL7R")

required_columns <- c(
  "Sample", "OS.time", "OS_event",
  "AOC3", "MS4A7", "IL7R"
)

missing_columns <- setdiff(required_columns, names(dat))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (nrow(dat) != 485L) {
  stop("Expected locked cohort N = 485; observed N = ", nrow(dat))
}

if (anyDuplicated(dat$Sample)) {
  stop("Duplicate Sample IDs detected. Stop before plotting.")
}

if (anyNA(dat[, required_columns])) {
  stop("Missing values detected in required locked-cohort fields.")
}

if (!all(dat$OS_event %in% c(0, 1))) {
  stop("OS_event must contain only 0 and 1.")
}

if (any(dat$OS.time < 0)) {
  stop("Negative OS.time detected.")
}


# ------------------------------------------------------------
# 1. Output directory
# ------------------------------------------------------------

if (exists("out_dir")) {
  base_out_dir <- out_dir
} else {
  base_out_dir <- file.path(
    "O:/Rproject",
    "Paper 1_Step 6-survival analysis",
    "STEP_6E_KM_curve_harmonization"
  )
}

final_out_dir <- file.path(
  base_out_dir,
  "KM_plots_Phase2B_C_minimal_revision"
)

dir.create(
  final_out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Initialise result objects
# ------------------------------------------------------------

km_fits <- list()
km_plot_objects <- list()
group_membership <- list()
logrank_results <- list()
grouped_cox_results <- list()
file_audit_list <- list()


# ------------------------------------------------------------
# 3. Run the same prespecified median-split KM analysis
# ------------------------------------------------------------

for (g in gene_list) {
  cat("\n--------------------------------------------\n")
  cat("Processing:", g, "\n")
  cat("--------------------------------------------\n")
  
  gene_dat <- data.frame(
    Sample = dat$Sample,
    OS.time = dat$OS.time,
    OS_event = dat$OS_event,
    Expression = dat[[g]],
    stringsAsFactors = FALSE
  )
  
  if (!is.numeric(gene_dat$Expression)) {
    stop(g, ": expression column is not numeric.")
  }
  
  if (anyNA(gene_dat$Expression)) {
    stop(g, ": missing expression values.")
  }
  
  # Prespecified median split:
  # expression <= median -> Low
  # expression > median  -> High
  
  med <- median(gene_dat$Expression)
  
  gene_dat$Group <- factor(
    ifelse(gene_dat$Expression <= med, "Low", "High"),
    levels = c("Low", "High")
  )
  
  group_counts <- table(gene_dat$Group)
  
  if (!identical(as.integer(group_counts), c(243L, 242L))) {
    stop(
      g, ": unexpected median-split group counts: ",
      paste(
        names(group_counts),
        as.integer(group_counts),
        collapse = "; "
      )
    )
  }
  
  event_counts <- tapply(
    gene_dat$OS_event,
    gene_dat$Group,
    sum
  )
  
  zero_counts <- tapply(
    gene_dat$OS.time == 0,
    gene_dat$Group,
    sum
  )
  
  cat("Median threshold:", format(med, digits = 12), "\n")
  cat(
    "Group sizes:",
    paste(names(group_counts), group_counts, collapse = "; "),
    "\n"
  )
  cat(
    "Events:",
    paste(names(event_counts), event_counts, collapse = "; "),
    "\n"
  )
  cat(
    "Zero-time cases:",
    paste(names(zero_counts), zero_counts, collapse = "; "),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 3A. Save locked group membership
  # ----------------------------------------------------------
  
  group_membership[[g]] <- data.frame(
    Gene = g,
    Sample = gene_dat$Sample,
    Expression = gene_dat$Expression,
    Median_threshold = med,
    Group = as.character(gene_dat$Group),
    OS.time = gene_dat$OS.time,
    OS_event = gene_dat$OS_event,
    stringsAsFactors = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 3B. Kaplan-Meier fit
  # ----------------------------------------------------------
  
  fit <- survival::survfit(
    survival::Surv(OS.time, OS_event) ~ Group,
    data = gene_dat,
    conf.type = "log"
  )
  
  km_fits[[g]] <- fit
  
  
  # ----------------------------------------------------------
  # 3C. Log-rank test
  # ----------------------------------------------------------
  
  lr <- survival::survdiff(
    survival::Surv(OS.time, OS_event) ~ Group,
    data = gene_dat,
    rho = 0
  )
  
  p_lr <- stats::pchisq(
    lr$chisq,
    df = length(lr$n) - 1,
    lower.tail = FALSE
  )
  
  logrank_results[[g]] <- data.frame(
    Gene = g,
    Median_threshold = med,
    N_total = nrow(gene_dat),
    N_Low = unname(group_counts["Low"]),
    N_High = unname(group_counts["High"]),
    Events_total = sum(gene_dat$OS_event),
    Events_Low = unname(event_counts["Low"]),
    Events_High = unname(event_counts["High"]),
    Zero_time_total = sum(gene_dat$OS.time == 0),
    Zero_time_events = sum(
      gene_dat$OS.time == 0 & gene_dat$OS_event == 1
    ),
    Logrank_chisq = unname(lr$chisq),
    Logrank_df = length(lr$n) - 1,
    Logrank_P = p_lr,
    stringsAsFactors = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 3D. Grouped Cox HR: High versus Low
  #
  # IMPORTANT:
  # This is a median-group HR, not the continuous-expression HR.
  # ----------------------------------------------------------
  
  grouped_cox_fit <- survival::coxph(
    survival::Surv(OS.time, OS_event) ~ Group,
    data = gene_dat,
    ties = "efron"
  )
  
  grouped_cox_summary <- summary(grouped_cox_fit)
  grouped_cox_ci <- exp(stats::confint(grouped_cox_fit))
  
  hr <- unname(
    exp(stats::coef(grouped_cox_fit)["GroupHigh"])
  )
  
  ci_low <- unname(
    grouped_cox_ci["GroupHigh", 1]
  )
  
  ci_high <- unname(
    grouped_cox_ci["GroupHigh", 2]
  )
  
  grouped_cox_p <- unname(
    grouped_cox_summary$coefficients[
      "GroupHigh", "Pr(>|z|)"
    ]
  )
  
  grouped_cox_results[[g]] <- data.frame(
    Gene = g,
    Median_threshold = med,
    N = nrow(gene_dat),
    Events_Low = unname(event_counts["Low"]),
    Events_High = unname(event_counts["High"]),
    HR_High_vs_Low = hr,
    CI95_lower = ci_low,
    CI95_upper = ci_high,
    Grouped_Cox_P = grouped_cox_p,
    stringsAsFactors = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 3E. Statistical annotation
  # ----------------------------------------------------------
  
  lr_label <- if (p_lr < 0.001) {
    "Log-rank P < 0.001"
  } else {
    sprintf("Log-rank P = %.3f", p_lr)
  }
  
  hr_label <- sprintf(
    "HR (High vs Low) = %.2f (95%% CI %.2f-%.2f)",
    hr, ci_low, ci_high
  )
  
  stat_label <- paste(
    lr_label,
    hr_label,
    sep = "\n"
  )
  
  
  # ----------------------------------------------------------
  # 3F. Construct the original-style KM plot
  #
  # Keep the original risk table, confidence intervals,
  # palette, title and labels.
  # ----------------------------------------------------------
  
  p <- survminer::ggsurvplot(
    fit,
    data = gene_dat,
    risk.table = TRUE,
    conf.int = TRUE,
    pval = FALSE,
    censor = TRUE,
    palette = c("#3C5488", "#E64B35"),
    legend.title = "Expression",
    legend.labs = c("Low", "High"),
    xlab = "Overall survival time (days)",
    ylab = "Overall survival probability",
    title = paste0(g, " expression and overall survival"),
    risk.table.title = "Number at risk",
    risk.table.height = 0.25,
    ggtheme = ggplot2::theme_classic(
      base_family = "Arial"
    ),
    break.time.by = 2000
  )
  
  
  # ----------------------------------------------------------
  # 3G. Add plain statistical text at lower-left
  #
  # No border, no background.
  # Left aligned near X = 0, Y = 0.
  # ----------------------------------------------------------
  
  p$plot <- p$plot +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0.035,
      label = stat_label,
      hjust = 0,
      vjust = 0,
      family = "Arial",
      fontface = "plain",
      size = 3
    )
  
  
  # ----------------------------------------------------------
  # 3H. Set KM axis fonts directly
  #
  # Direct theme-component assignment avoids the theme-merging
  # error previously encountered with ggplot2 4.0.3.
  # ----------------------------------------------------------
  
  p$plot$theme$axis.title.x <- ggplot2::element_text(
    family = "Arial",
    face = "bold",
    size = 10
  )
  
  p$plot$theme$axis.title.y <- ggplot2::element_text(
    family = "Arial",
    face = "bold",
    size = 10,
    angle = 90
  )
  
  p$plot$theme$axis.text.x <- ggplot2::element_text(
    family = "Arial",
    face = "bold",
    size = 10
  )
  
  p$plot$theme$axis.text.y <- ggplot2::element_text(
    family = "Arial",
    face = "bold",
    size = 10,
    angle = 90
  )
  
  # Keep the risk table; apply matching typography to its axes.
  if (!is.null(p$table)) {
    
    p$table$theme$axis.text.x <- ggplot2::element_text(
      family = "Arial",
      face = "bold",
      size = 10
    )
    
    p$table$theme$axis.text.y <- ggplot2::element_text(
      family = "Arial",
      face = "bold",
      size = 10,
      angle = 90
    )
  }
  
  km_plot_objects[[g]] <- p
  
  
  # ----------------------------------------------------------
  # 3I. Combine KM plot and risk table into one-page object
  # ----------------------------------------------------------
  
  combined_plot <- survminer::arrange_ggsurvplots(
    list(p),
    print = FALSE,
    ncol = 1,
    nrow = 1,
    surv.plot.height = 0.75,
    risk.table.height = 0.25
  )
  
  
  # ----------------------------------------------------------
  # 3J. Export TIFF and PDF
  # ----------------------------------------------------------
  
  tiff_file <- file.path(
    final_out_dir,
    paste0("STEP6E_KM_", g, "_locked485_final.tiff")
  )
  
  pdf_file <- file.path(
    final_out_dir,
    paste0("STEP6E_KM_", g, "_locked485_final.pdf")
  )
  
  ggplot2::ggsave(
    filename = tiff_file,
    plot = combined_plot,
    device = "tiff",
    width = 7,
    height = 6,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
  
  ggplot2::ggsave(
    filename = pdf_file,
    plot = combined_plot,
    device = grDevices::cairo_pdf,
    width = 7,
    height = 6,
    units = "in",
    bg = "white"
  )
  
  
  # ----------------------------------------------------------
  # 3K. File audit
  # ----------------------------------------------------------
  
  file_audit_list[[g]] <- data.frame(
    Gene = g,
    TIFF = tiff_file,
    TIFF_exists = file.exists(tiff_file),
    TIFF_bytes = if (file.exists(tiff_file)) {
      file.info(tiff_file)$size
    } else {
      NA_real_
    },
    PDF = pdf_file,
    PDF_exists = file.exists(pdf_file),
    PDF_bytes = if (file.exists(pdf_file)) {
      file.info(pdf_file)$size
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
  
  cat("Log-rank P:", signif(p_lr, 6), "\n")
  cat(
    "Grouped Cox HR:",
    round(hr, 3),
    "(95% CI",
    round(ci_low, 3), "-",
    round(ci_high, 3), ")\n"
  )
  
}

# ------------------------------------------------------------
# 4. Save results and audit tables
# ------------------------------------------------------------

group_membership_table <- do.call(
  rbind,
  group_membership
)

logrank_table <- do.call(
  rbind,
  logrank_results
)

grouped_cox_table <- do.call(
  rbind,
  grouped_cox_results
)

file_audit <- do.call(
  rbind,
  file_audit_list
)

rownames(group_membership_table) <- NULL
rownames(logrank_table) <- NULL
rownames(grouped_cox_table) <- NULL
rownames(file_audit) <- NULL

write.csv(
  group_membership_table,
  file.path(
    final_out_dir,
    "STEP6E_Phase2D_locked_group_membership.csv"
  ),
  row.names = FALSE
)

write.csv(
  logrank_table,
  file.path(
    final_out_dir,
    "STEP6E_Phase2D_logrank_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  grouped_cox_table,
  file.path(
    final_out_dir,
    "STEP6E_Phase2D_grouped_Cox_HR_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  file_audit,
  file.path(
    final_out_dir,
    "STEP6E_Phase2D_plot_file_audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Final console output
# ------------------------------------------------------------

cat("\n============================================\n")
cat("FINAL LOG-RANK RESULTS\n")
cat("============================================\n")
print(logrank_table)

cat("\n============================================\n")
cat("FINAL GROUPED COX HR RESULTS\n")
cat("============================================\n")
print(grouped_cox_table)

cat("\n============================================\n")
cat("OUTPUT FILE AUDIT\n")
cat("============================================\n")
print(file_audit)

cat("\nFinal output directory:\n")
cat(normalizePath(final_out_dir), "\n")









