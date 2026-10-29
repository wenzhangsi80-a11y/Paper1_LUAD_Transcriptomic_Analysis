
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
# STEP 18A PRE-CHECK
# Identify existing STEP 14 objects
# ============================================================

cat("\n============================================================\n")
cat("STEP 18A PRE-CHECK: Existing objects\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. List objects currently in the Global Environment
# ------------------------------------------------------------

all_objects <- ls(envir = .GlobalEnv)

cat("Number of objects in Global Environment:",
    length(all_objects), "\n\n")

print(all_objects)


# ------------------------------------------------------------
# 2. Show objects whose names are potentially related to
#    candidate genes / TME / GS / correlation / WGCNA
# ------------------------------------------------------------

pattern_objects <- grep(
  pattern = paste(
    c(
      "candidate",
      "cand",
      "TME",
      "tme",
      "GS",
      "gene",
      "cor",
      "spearman",
      "FDR",
      "WGCNA",
      "module"
    ),
    collapse = "|"
  ),
  x = all_objects,
  value = TRUE
)

cat("\n============================================================\n")
cat("Potentially relevant objects\n")
cat("============================================================\n\n")

print(pattern_objects)


# ------------------------------------------------------------
# 3. Show object classes and dimensions
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Object class / dimensions\n")
cat("============================================================\n\n")

for (obj_name in pattern_objects) {
  
  obj <- get(
    obj_name,
    envir = .GlobalEnv
  )
  
  cat("\n------------------------------------------------------------\n")
  cat("OBJECT:", obj_name, "\n")
  cat("CLASS :", paste(class(obj), collapse = ", "), "\n")
  
  if (is.data.frame(obj) || is.matrix(obj)) {
    
    cat(
      "DIM   :",
      paste(dim(obj), collapse = " x "),
      "\n"
    )
    
    cat("COLS  :\n")
    print(colnames(obj))
    
  } else if (is.vector(obj) || is.factor(obj)) {
    
    cat(
      "LENGTH:",
      length(obj),
      "\n"
    )
    
  }
  
}

# ============================================================
# Paper 1
# STEP 18A
# 30-candidate gene-level TME association audit
#
# FINAL LOCKED VERSION
#
# INPUT:
#   candidate_30
#   candidate_master
#
# PURPOSE:
#   Audit the already-computed gene-level TME association
#   results for the 30 candidate genes.
#
# IMPORTANT:
#   - NO new Spearman correlation
#   - NO recalculation of P values
#   - NO recalculation of FDR
#   - NO new candidate selection
#   - NO modification of STEP 14 results
#
# ============================================================


cat("\n============================================================\n")
cat("STEP 18A: 30-candidate gene-level TME association audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Check required objects
# ------------------------------------------------------------

cat("STEP 18A.1: Checking required objects...\n\n")

required_objects <- c(
  "candidate_30",
  "candidate_master"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    envir = .GlobalEnv
  )
]

if (length(missing_objects) > 0) {
  
  stop(
    paste0(
      "ERROR: Required object(s) missing:\n",
      paste(
        missing_objects,
        collapse = ", "
      ),
      "\n"
    )
  )
  
}

cat("PASS: Required objects are present.\n\n")


# ------------------------------------------------------------
# 2. Verify candidate_30
# ------------------------------------------------------------

cat("STEP 18A.2: Checking candidate_30...\n\n")

cat(
  "Number of candidate genes:",
  length(candidate_30),
  "\n"
)

cat(
  "Unique candidate genes:",
  length(unique(candidate_30)),
  "\n\n"
)

if (length(candidate_30) != 30) {
  
  stop(
    "ERROR: candidate_30 does not contain exactly 30 genes."
  )
  
}

if (length(unique(candidate_30)) != 30) {
  
  stop(
    "ERROR: candidate_30 contains duplicated genes."
  )
  
}

cat("PASS: candidate_30 contains 30 unique genes.\n\n")


# ------------------------------------------------------------
# 3. Verify candidate_master structure
# ------------------------------------------------------------

cat("STEP 18A.3: Checking candidate_master structure...\n\n")

cat(
  "Dimensions:",
  paste(
    dim(candidate_master),
    collapse = " x "
  ),
  "\n\n"
)

required_columns <- c(
  "Candidate_Order",
  "Symbol",
  "ENSEMBL",
  "Core3",
  "TumorOnly_Module",
  "TumorOnly_ME",
  "MM",
  "GS_Immune",
  "P_Immune",
  "FDR_Immune",
  "GS_Stromal",
  "P_Stromal",
  "FDR_Stromal",
  "GS_ESTIMATE",
  "P_ESTIMATE",
  "FDR_ESTIMATE"
)

missing_columns <- setdiff(
  required_columns,
  colnames(candidate_master)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste0(
      "ERROR: candidate_master is missing columns:\n",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
  
}

cat("PASS: candidate_master contains all expected columns.\n\n")


# ------------------------------------------------------------
# 4. Verify 30-gene identity
# ------------------------------------------------------------

cat("STEP 18A.4: Verifying candidate identity...\n\n")

candidate_master_genes <- as.character(
  candidate_master$Symbol
)

missing_from_master <- setdiff(
  candidate_30,
  candidate_master_genes
)

extra_in_master <- setdiff(
  candidate_master_genes,
  candidate_30
)

if (length(missing_from_master) > 0) {
  
  cat(
    "WARNING: Candidate genes missing from candidate_master:\n"
  )
  
  print(missing_from_master)
  
} else {
  
  cat(
    "PASS: All candidate_30 genes are present in candidate_master.\n"
  )
  
}

if (length(extra_in_master) > 0) {
  
  cat(
    "\nWARNING: candidate_master contains genes not present in candidate_30:\n"
  )
  
  print(extra_in_master)
  
}

cat("\n")


# ------------------------------------------------------------
# 5. Create STEP 18A audit table
#
# IMPORTANT:
#   All GS/P/FDR values are COPIED from candidate_master.
#   Nothing is recalculated.
# ------------------------------------------------------------

cat("STEP 18A.5: Creating audit table...\n\n")

audit_30 <- candidate_master[
  ,
  required_columns,
  drop = FALSE
]


# Preserve candidate_30 order
audit_30$Candidate_Order <- match(
  audit_30$Symbol,
  candidate_30
)

audit_30 <- audit_30[
  order(audit_30$Candidate_Order),
  ,
  drop = FALSE
]

rownames(audit_30) <- NULL


# ------------------------------------------------------------
# 6. Confirm locked GS/P/FDR values are present
# ------------------------------------------------------------

cat("STEP 18A.6: Checking GS/P/FDR completeness...\n\n")

association_columns <- c(
  "GS_Immune",
  "P_Immune",
  "FDR_Immune",
  "GS_Stromal",
  "P_Stromal",
  "FDR_Stromal",
  "GS_ESTIMATE",
  "P_ESTIMATE",
  "FDR_ESTIMATE"
)

missing_count <- sapply(
  audit_30[
    ,
    association_columns,
    drop = FALSE
  ],
  function(x) sum(is.na(x))
)

print(missing_count)

cat("\n")


# ------------------------------------------------------------
# 7. Descriptive direction classification
#
# IMPORTANT:
#   This does NOT calculate anything new.
#   It only labels the existing GS values.
# ------------------------------------------------------------

cat("STEP 18A.7: Classifying existing GS directions...\n\n")

audit_30$Direction_Immune <- ifelse(
  audit_30$GS_Immune > 0,
  "Positive",
  ifelse(
    audit_30$GS_Immune < 0,
    "Negative",
    "Zero"
  )
)

audit_30$Direction_Stromal <- ifelse(
  audit_30$GS_Stromal > 0,
  "Positive",
  ifelse(
    audit_30$GS_Stromal < 0,
    "Negative",
    "Zero"
  )
)

audit_30$Direction_ESTIMATE <- ifelse(
  audit_30$GS_ESTIMATE > 0,
  "Positive",
  ifelse(
    audit_30$GS_ESTIMATE < 0,
    "Negative",
    "Zero"
  )
)


# ------------------------------------------------------------
# 8. Existing FDR significance flags
#
# Descriptive classification only.
# No FDR is recalculated.
# ------------------------------------------------------------

FDR_threshold <- 0.05

audit_30$FDR_Immune_status <- ifelse(
  !is.na(audit_30$FDR_Immune) &
    audit_30$FDR_Immune < FDR_threshold,
  "FDR < 0.05",
  "FDR >= 0.05"
)

audit_30$FDR_Stromal_status <- ifelse(
  !is.na(audit_30$FDR_Stromal) &
    audit_30$FDR_Stromal < FDR_threshold,
  "FDR < 0.05",
  "FDR >= 0.05"
)

audit_30$FDR_ESTIMATE_status <- ifelse(
  !is.na(audit_30$FDR_ESTIMATE) &
    audit_30$FDR_ESTIMATE < FDR_threshold,
  "FDR < 0.05",
  "FDR >= 0.05"
)


# ------------------------------------------------------------
# 9. Calculate descriptive mean absolute GS
#
# This is NOT a new statistical test.
# It is only a summary of existing GS values.
# ------------------------------------------------------------

audit_30$Mean_abs_GS <- rowMeans(
  abs(
    audit_30[
      ,
      c(
        "GS_Immune",
        "GS_Stromal",
        "GS_ESTIMATE"
      )
    ]
  ),
  na.rm = TRUE
)


# ------------------------------------------------------------
# 10. Final audit table
# ------------------------------------------------------------

audit_30_final <- audit_30[
  ,
  c(
    "Candidate_Order",
    "Symbol",
    "ENSEMBL",
    "Core3",
    "TumorOnly_Module",
    "TumorOnly_ME",
    "MM",
    "GS_Immune",
    "P_Immune",
    "FDR_Immune",
    "GS_Stromal",
    "P_Stromal",
    "FDR_Stromal",
    "GS_ESTIMATE",
    "P_ESTIMATE",
    "FDR_ESTIMATE",
    "Direction_Immune",
    "Direction_Stromal",
    "Direction_ESTIMATE",
    "FDR_Immune_status",
    "FDR_Stromal_status",
    "FDR_ESTIMATE_status",
    "Mean_abs_GS"
  )
]

rownames(audit_30_final) <- NULL


# ------------------------------------------------------------
# 11. Display complete audit table
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 18A: COMPLETE 30-CANDIDATE AUDIT TABLE\n")
cat("============================================================\n\n")

print(
  audit_30_final,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 12. Core 3 audit
#
# Locked Paper 1 core genes:
#   AOC3
#   MS4A7
#   IL7R
# ------------------------------------------------------------

cat("STEP 18A.8: Locked core-gene audit\n\n")

core_genes_locked <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_audit <- audit_30_final[
  audit_30_final$Symbol %in%
    core_genes_locked,
  ,
  drop = FALSE
]

core_audit <- core_audit[
  match(
    core_genes_locked,
    core_audit$Symbol
  ),
  ,
  drop = FALSE
]

print(
  core_audit,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 13. Count FDR-significant candidate associations
#
# This only summarizes existing FDR values.
# ------------------------------------------------------------

cat("STEP 18A.9: FDR significance summary\n\n")

cat(
  "Immune FDR < 0.05:",
  sum(
    audit_30$FDR_Immune < FDR_threshold,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Stromal FDR < 0.05:",
  sum(
    audit_30$FDR_Stromal < FDR_threshold,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "ESTIMATE FDR < 0.05:",
  sum(
    audit_30$FDR_ESTIMATE < FDR_threshold,
    na.rm = TRUE
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 14. Show candidates with any significant TME association
# ------------------------------------------------------------

any_FDR_significant <- (
  audit_30$FDR_Immune < FDR_threshold |
    audit_30$FDR_Stromal < FDR_threshold |
    audit_30$FDR_ESTIMATE < FDR_threshold
)

audit_any_FDR <- audit_30_final[
  any_FDR_significant,
  ,
  drop = FALSE
]

cat("STEP 18A.10: Candidates with any FDR < 0.05\n\n")

if (nrow(audit_any_FDR) == 0) {
  
  cat(
    "No candidate gene has FDR < 0.05 for any of the three TME traits.\n"
  )
  
} else {
  
  print(
    audit_any_FDR,
    row.names = FALSE
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 15. Rank candidates by mean absolute GS
#
# Descriptive ranking only.
# ------------------------------------------------------------

cat("STEP 18A.11: Descriptive ranking by mean |GS|\n\n")

audit_ranked <- audit_30_final[
  order(
    -audit_30_final$Mean_abs_GS,
    audit_30_final$Candidate_Order
  ),
  ,
  drop = FALSE
]

print(
  audit_ranked[
    ,
    c(
      "Candidate_Order",
      "Symbol",
      "Core3",
      "GS_Immune",
      "GS_Stromal",
      "GS_ESTIMATE",
      "Mean_abs_GS"
    )
  ],
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 16. Final integrity checks
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP 18A.12: FINAL INTEGRITY CHECK\n")
cat("============================================================\n\n")

cat(
  "candidate_30 length:",
  length(candidate_30),
  "\n"
)

cat(
  "candidate_master rows:",
  nrow(candidate_master),
  "\n"
)

cat(
  "audit_30_final rows:",
  nrow(audit_30_final),
  "\n"
)

cat(
  "Core genes detected:",
  sum(
    audit_30_final$Symbol %in%
      core_genes_locked
  ),
  "/ 3\n"
)

cat("\n")

cat("NO correlation recalculation performed.\n")
cat("NO P-value recalculation performed.\n")
cat("NO FDR recalculation performed.\n")
cat("NO new candidate selection performed.\n")

cat("\n")


# ------------------------------------------------------------
# 17. Export
# ------------------------------------------------------------

cat("STEP 18A.13: Exporting audit results...\n\n")

write.csv(
  audit_30_final,
  file = "STEP18A_30_candidate_gene_TME_audit.csv",
  row.names = FALSE
)

write.csv(
  core_audit,
  file = "STEP18A_core_gene_TME_audit.csv",
  row.names = FALSE
)

write.csv(
  audit_ranked,
  file = "STEP18A_30_candidate_gene_TME_audit_ranked.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. Save R objects
# ------------------------------------------------------------

saveRDS(
  audit_30_final,
  file = "STEP18A_30_candidate_gene_TME_audit.rds"
)

saveRDS(
  core_audit,
  file = "STEP18A_core_gene_TME_audit.rds"
)

saveRDS(
  audit_ranked,
  file = "STEP18A_30_candidate_gene_TME_audit_ranked.rds"
)


cat("\n============================================================\n")
cat("STEP 18A COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat(
  "STEP 18A is AUDIT-ONLY.\n",
  "All GS/P/FDR values were read directly from candidate_master.\n",
  "No gene-level correlation was recalculated.\n",
  "No P value was recalculated.\n",
  "No FDR was recalculated.\n",
  "No candidate gene was newly selected.\n",
  sep = ""
)

cat("\n============================================================\n")


##STEP 18B — 30-candidate WGCNA module/MM mapping integrity audit
# ============================================================
# Paper 1
# STEP 18B
# 30-candidate WGCNA module/MM mapping integrity audit
#
# PURPOSE:
#   Audit the integrity of existing WGCNA module/MM mappings
#   for the locked 30 candidate genes.
#
# STRICTLY AUDIT-ONLY:
#   - NO WGCNA recalculation
#   - NO module eigengene recalculation
#   - NO MM recalculation
#   - NO GS recalculation
#   - NO P-value recalculation
#   - NO FDR recalculation
#   - NO candidate reselection
#   - NO modification of Step 14 results
#
# MAIN QUESTION:
#   Why are TumorOnly_Module / TumorOnly_ME / MM all NA
#   in candidate_master for the 30 candidates?
#
# ============================================================

cat("\n============================================================\n")
cat("STEP 18B: 30-candidate WGCNA module/MM mapping integrity audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Define the locked candidate set
# ------------------------------------------------------------

core_genes_expected <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

stopifnot(exists("candidate_30"))
stopifnot(exists("candidate_master"))

candidate_30_locked <- unique(as.character(candidate_30))

cat("Number of candidate genes:", length(candidate_30_locked), "\n")
cat("Unique candidate genes:", length(unique(candidate_30_locked)), "\n")

if (length(candidate_30_locked) != 30) {
  stop("ERROR: candidate_30 does not contain exactly 30 unique genes.")
}

cat("\nLocked core genes:\n")
print(core_genes_expected)


# ------------------------------------------------------------
# 1. Required object audit
# ------------------------------------------------------------

required_objects <- c(
  "candidate_30",
  "candidate_master",
  "candidate_master_16D",
  "candidate_map_raw",
  "moduleColors",
  "module_genes",
  "module_names_colors",
  "module_names_MM",
  "core_step17",
  "core_position_final",
  "core_final_integrated_save",
  "core_visualization_table"
)

object_status <- data.frame(
  Object = required_objects,
  Exists = vapply(
    required_objects,
    exists,
    logical(1)
  ),
  stringsAsFactors = FALSE
)

cat("\n------------------------------------------------------------\n")
cat("1. Required object existence audit\n")
cat("------------------------------------------------------------\n")

print(object_status, row.names = FALSE)

if (!all(object_status$Exists)) {
  missing_objects <- object_status$Object[!object_status$Exists]
  
  stop(
    paste(
      "ERROR: Required object(s) missing:",
      paste(missing_objects, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 2. Inspect candidate_master WGCNA-related columns
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("2. candidate_master WGCNA-related columns\n")
cat("------------------------------------------------------------\n")

candidate_master_wgcna_cols <- intersect(
  c(
    "Symbol",
    "TumorOnly_Module",
    "TumorOnly_ME",
    "MM",
    "Original_Module"
  ),
  colnames(candidate_master)
)

print(candidate_master[, candidate_master_wgcna_cols, drop = FALSE])


# ------------------------------------------------------------
# 3. Quantify NA pattern in candidate_master
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("3. candidate_master NA pattern\n")
cat("------------------------------------------------------------\n")

wgcna_mapping_cols <- c(
  "TumorOnly_Module",
  "TumorOnly_ME",
  "MM"
)

na_summary <- data.frame(
  Variable = wgcna_mapping_cols,
  N_total = nrow(candidate_master),
  N_nonNA = sapply(
    candidate_master[wgcna_mapping_cols],
    function(x) sum(!is.na(x))
  ),
  N_NA = sapply(
    candidate_master[wgcna_mapping_cols],
    function(x) sum(is.na(x))
  ),
  stringsAsFactors = FALSE
)

print(na_summary, row.names = FALSE)


# ------------------------------------------------------------
# 4. Check candidate_map_raw
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("4. candidate_map_raw audit\n")
cat("------------------------------------------------------------\n")

print(
  head(
    candidate_map_raw,
    10
  )
)

required_candidate_map_cols <- c(
  "SYMBOL",
  "ENSEMBL",
  "In_WGCNA"
)

if (!all(required_candidate_map_cols %in% colnames(candidate_map_raw))) {
  
  stop(
    paste(
      "ERROR: candidate_map_raw is missing required columns:",
      paste(
        setdiff(
          required_candidate_map_cols,
          colnames(candidate_map_raw)
        ),
        collapse = ", "
      )
    )
  )
}

candidate_map_audit <- candidate_map_raw[
  match(
    candidate_30_locked,
    candidate_map_raw$SYMBOL
  ),
  ,
  drop = FALSE
]

candidate_map_audit$Candidate <- candidate_30_locked

candidate_map_audit <- candidate_map_audit[
  ,
  c(
    "Candidate",
    "SYMBOL",
    "ENSEMBL",
    "In_WGCNA"
  ),
  drop = FALSE
]

cat("\nCandidate → candidate_map_raw mapping:\n")
print(
  candidate_map_audit,
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Check whether all 30 candidates are represented in
#    moduleColors
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("5. moduleColors membership audit\n")
cat("------------------------------------------------------------\n")

cat("Length(moduleColors):", length(moduleColors), "\n")
cat("Number of unique moduleColors:",
    length(unique(moduleColors)),
    "\n")

cat("\nmoduleColors object structure:\n")
print(
  head(
    moduleColors,
    10
  )
)


# IMPORTANT:
# moduleColors itself is only a vector of module assignments.
# We need to determine whether its names identify genes.

cat("\nDoes moduleColors have gene names?\n")
cat(
  "names(moduleColors) is NULL: ",
  is.null(names(moduleColors)),
  "\n",
  sep = ""
)

if (!is.null(names(moduleColors))) {
  
  candidate_in_moduleColors <- candidate_30_locked %in%
    names(moduleColors)
  
  moduleColors_membership <- data.frame(
    Symbol = candidate_30_locked,
    In_moduleColors = candidate_in_moduleColors,
    Module = unname(
      moduleColors[
        match(
          candidate_30_locked,
          names(moduleColors)
        )
      ]
    ),
    stringsAsFactors = FALSE
  )
  
} else {
  
  moduleColors_membership <- data.frame(
    Symbol = candidate_30_locked,
    In_moduleColors = NA,
    Module = NA_character_,
    stringsAsFactors = FALSE
  )
  
  cat(
    "\nNOTE: moduleColors has no names.\n",
    "Direct gene-level lookup is therefore not possible from this object alone.\n",
    sep = ""
  )
}

print(
  moduleColors_membership,
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Audit module_genes
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("6. module_genes audit\n")
cat("------------------------------------------------------------\n")

cat("Length(module_genes):", length(module_genes), "\n")

cat("\nFirst entries:\n")
print(
  head(
    module_genes,
    20
  )
)

candidate_in_module_genes <- candidate_30_locked %in%
  module_genes

module_genes_audit <- data.frame(
  Symbol = candidate_30_locked,
  In_module_genes = candidate_in_module_genes,
  stringsAsFactors = FALSE
)

print(
  module_genes_audit,
  row.names = FALSE
)

cat(
  "\nNumber of candidates in module_genes:",
  sum(candidate_in_module_genes),
  "/",
  length(candidate_30_locked),
  "\n"
)


# ------------------------------------------------------------
# 7. Audit module_names_colors
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("7. module_names_colors audit\n")
cat("------------------------------------------------------------\n")

cat("Length(module_names_colors):",
    length(module_names_colors),
    "\n")

print(module_names_colors)


# ------------------------------------------------------------
# 8. Audit module_names_MM
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("8. module_names_MM audit\n")
cat("------------------------------------------------------------\n")

cat("Length(module_names_MM):",
    length(module_names_MM),
    "\n")

print(module_names_MM)


# ------------------------------------------------------------
# 9. Compare module_names_colors and module_names_MM
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("9. Module color ↔ MM mapping consistency\n")
cat("------------------------------------------------------------\n")

module_name_mapping_audit <- data.frame(
  Module = module_names_colors,
  MM_value = module_names_MM[
    match(
      module_names_colors,
      names(module_names_MM)
    )
  ],
  stringsAsFactors = FALSE
)

print(
  module_name_mapping_audit,
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Inspect core_step17
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("10. Locked core_step17 audit\n")
cat("------------------------------------------------------------\n")

core_step17_cols <- intersect(
  c(
    "Symbol",
    "TumorOnly_Module",
    "Module_Size",
    "MM",
    "MM_Rank",
    "MM_Percentile"
  ),
  colnames(core_step17)
)

print(
  core_step17[, core_step17_cols, drop = FALSE],
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Inspect core_position_final
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("11. Locked core_position_final audit\n")
cat("------------------------------------------------------------\n")

core_position_cols <- intersect(
  c(
    "Symbol",
    "TumorOnly_Module",
    "Module_Size",
    "MM",
    "MM_Rank",
    "MM_Percentile"
  ),
  colnames(core_position_final)
)

print(
  core_position_final[, core_position_cols, drop = FALSE],
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Inspect core_final_integrated_save
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("12. Locked core_final_integrated_save audit\n")
cat("------------------------------------------------------------\n")

core_integrated_cols <- intersect(
  c(
    "Symbol",
    "TumorOnly_Module",
    "TumorOnly_ME",
    "MM"
  ),
  colnames(core_final_integrated_save)
)

print(
  core_final_integrated_save[
    ,
    core_integrated_cols,
    drop = FALSE
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Inspect core_visualization_table
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("13. Locked core_visualization_table audit\n")
cat("------------------------------------------------------------\n")

core_visual_cols <- intersect(
  c(
    "Symbol",
    "TumorOnly_Module",
    "TumorOnly_ME",
    "MM"
  ),
  colnames(core_visualization_table)
)

print(
  core_visualization_table[
    ,
    core_visual_cols,
    drop = FALSE
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Direct comparison:
#     candidate_master vs locked core tables
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("14. Direct core-gene mapping comparison\n")
cat("------------------------------------------------------------\n")

core_compare <- data.frame(
  Symbol = core_genes_expected,
  stringsAsFactors = FALSE
)

core_compare$candidate_master_Module <-
  candidate_master$TumorOnly_Module[
    match(
      core_compare$Symbol,
      candidate_master$Symbol
    )
  ]

core_compare$candidate_master_ME <-
  candidate_master$TumorOnly_ME[
    match(
      core_compare$Symbol,
      candidate_master$Symbol
    )
  ]

core_compare$candidate_master_MM <-
  candidate_master$MM[
    match(
      core_compare$Symbol,
      candidate_master$Symbol
    )
  ]

core_compare$core_step17_Module <-
  core_step17$TumorOnly_Module[
    match(
      core_compare$Symbol,
      core_step17$Symbol
    )
  ]

core_compare$core_step17_MM <-
  core_step17$MM[
    match(
      core_compare$Symbol,
      core_step17$Symbol
    )
  ]

core_compare$core_position_Module <-
  core_position_final$TumorOnly_Module[
    match(
      core_compare$Symbol,
      core_position_final$Symbol
    )
  ]

core_compare$core_position_MM <-
  core_position_final$MM[
    match(
      core_compare$Symbol,
      core_position_final$Symbol
    )
  ]

print(
  core_compare,
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. Determine whether the discrepancy is an integration issue
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("15. Mapping discrepancy diagnosis\n")
cat("------------------------------------------------------------\n")

core_candidate_master_all_na <-
  all(
    is.na(core_compare$candidate_master_Module) &
      is.na(core_compare$candidate_master_ME) &
      is.na(core_compare$candidate_master_MM)
  )

core_locked_tables_have_mapping <-
  all(
    !is.na(core_compare$core_step17_Module) &
      !is.na(core_compare$core_step17_MM)
  )

if (
  core_candidate_master_all_na &&
  core_locked_tables_have_mapping
) {
  
  cat(
    "\nDIAGNOSIS: candidate_master integration/mapping problem is strongly indicated.\n",
    sep = ""
  )
  
  cat(
    "The locked core tables contain valid module/MM values,\n",
    "whereas candidate_master contains NA for the same core genes.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "\nDIAGNOSIS: discrepancy requires further object-level inspection.\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 16. Build a consolidated 30-candidate mapping audit
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("16. Consolidated 30-candidate mapping audit\n")
cat("------------------------------------------------------------\n")

mapping_audit_30 <- data.frame(
  Symbol = candidate_30_locked,
  stringsAsFactors = FALSE
)

mapping_audit_30$candidate_master_Module <-
  candidate_master$TumorOnly_Module[
    match(
      mapping_audit_30$Symbol,
      candidate_master$Symbol
    )
  ]

mapping_audit_30$candidate_master_ME <-
  candidate_master$TumorOnly_ME[
    match(
      mapping_audit_30$Symbol,
      candidate_master$Symbol
    )
  ]

mapping_audit_30$candidate_master_MM <-
  candidate_master$MM[
    match(
      mapping_audit_30$Symbol,
      candidate_master$Symbol
    )
  ]

mapping_audit_30$candidate_map_In_WGCNA <-
  candidate_map_raw$In_WGCNA[
    match(
      mapping_audit_30$Symbol,
      candidate_map_raw$SYMBOL
    )
  ]

if (!is.null(names(moduleColors))) {
  
  mapping_audit_30$moduleColors_Module <-
    unname(
      moduleColors[
        match(
          mapping_audit_30$Symbol,
          names(moduleColors)
        )
      ]
    )
  
  mapping_audit_30$moduleColors_Mapped <-
    !is.na(
      mapping_audit_30$moduleColors_Module
    )
  
} else {
  
  mapping_audit_30$moduleColors_Module <-
    NA_character_
  
  mapping_audit_30$moduleColors_Mapped <-
    NA
}

mapping_audit_30$In_module_genes <-
  mapping_audit_30$Symbol %in% module_genes


print(
  mapping_audit_30,
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Check consistency of candidate_master vs candidate_map_raw
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("17. candidate_master ↔ candidate_map_raw consistency\n")
cat("------------------------------------------------------------\n")

candidate_master_symbols <- as.character(
  candidate_master$Symbol
)

candidate_map_symbols <- as.character(
  candidate_map_raw$SYMBOL
)

cat(
  "Candidates missing from candidate_master:",
  sum(
    !candidate_30_locked %in%
      candidate_master_symbols
  ),
  "\n"
)

cat(
  "Candidates missing from candidate_map_raw:",
  sum(
    !candidate_30_locked %in%
      candidate_map_symbols
  ),
  "\n"
)

cat(
  "candidate_master contains duplicate candidate symbols:",
  anyDuplicated(candidate_master_symbols),
  "\n"
)

cat(
  "candidate_map_raw contains duplicate candidate symbols:",
  anyDuplicated(candidate_map_symbols),
  "\n"
)


# ------------------------------------------------------------
# 18. Check exact candidate order consistency
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("18. Candidate order consistency\n")
cat("------------------------------------------------------------\n")

candidate_order_check <- data.frame(
  Position = seq_along(candidate_30_locked),
  candidate_30 = candidate_30_locked,
  candidate_master = candidate_master$Symbol[
    match(
      candidate_30_locked,
      candidate_master$Symbol
    )
  ],
  stringsAsFactors = FALSE
)

candidate_order_check$Same <-
  candidate_order_check$candidate_30 ==
  candidate_order_check$candidate_master

print(
  candidate_order_check,
  row.names = FALSE
)

cat(
  "\nAll candidate symbols correctly mapped:",
  all(candidate_order_check$Same),
  "\n"
)


# ------------------------------------------------------------
# 19. Core-gene identity check
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("19. Locked core-gene identity check\n")
cat("------------------------------------------------------------\n")

core_identity <- data.frame(
  Core_Gene = core_genes_expected,
  In_candidate_30 =
    core_genes_expected %in% candidate_30_locked,
  In_candidate_master =
    core_genes_expected %in% candidate_master$Symbol,
  In_candidate_map_raw =
    core_genes_expected %in% candidate_map_raw$SYMBOL,
  stringsAsFactors = FALSE
)

print(
  core_identity,
  row.names = FALSE
)


# ------------------------------------------------------------
# 20. Final audit summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 18B FINAL AUDIT SUMMARY\n")
cat("============================================================\n")

n_candidates <- length(candidate_30_locked)

n_master_module <- sum(
  !is.na(candidate_master$TumorOnly_Module)
)

n_master_me <- sum(
  !is.na(candidate_master$TumorOnly_ME)
)

n_master_mm <- sum(
  !is.na(candidate_master$MM)
)

n_candidate_map <- sum(
  candidate_30_locked %in%
    candidate_map_raw$SYMBOL
)

n_module_genes <- sum(
  candidate_30_locked %in%
    module_genes
)

cat(
  "30 candidates:",
  n_candidates,
  "\n"
)

cat(
  "candidate_master with module:",
  n_master_module,
  "/",
  n_candidates,
  "\n"
)

cat(
  "candidate_master with ME:",
  n_master_me,
  "/",
  n_candidates,
  "\n"
)

cat(
  "candidate_master with MM:",
  n_master_mm,
  "/",
  n_candidates,
  "\n"
)

cat(
  "candidate_map_raw represented:",
  n_candidate_map,
  "/",
  n_candidates,
  "\n"
)

cat(
  "module_genes represented:",
  n_module_genes,
  "/",
  n_candidates,
  "\n"
)


# ------------------------------------------------------------
# 21. PASS / REVIEW logic
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("21. STEP 18B audit status\n")
cat("------------------------------------------------------------\n")

if (
  core_candidate_master_all_na &&
  core_locked_tables_have_mapping
) {
  
  cat(
    "STATUS: REVIEW REQUIRED — likely candidate_master mapping/integration issue.\n"
  )
  
} else {
  
  cat(
    "STATUS: COMPLETE — no definitive integration discrepancy identified at this stage.\n"
  )
}


cat("\n============================================================\n")
cat("STEP 18B completed.\n")
cat("NO WGCNA / GS / P / FDR recalculation was performed.\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 18C
# ENSEMBL-based WGCNA module membership integrity audit
#
# PURPOSE:
#   Audit whether the 30 locked candidates can be mapped through
#   SYMBOL -> ENSEMBL -> moduleColors and whether the original
#   WGCNA module assignments of AOC3 / MS4A7 / IL7R can be
#   independently recovered.
#
# IMPORTANT:
#   - AUDIT ONLY
#   - NO WGCNA recalculation
#   - NO MM recalculation
#   - NO GS / P / FDR recalculation
#   - NO candidate reselection
#   - NO modification of candidate_master
#   - NO modification of locked core results
# ============================================================

cat("\n============================================================\n")
cat("STEP 18C: ENSEMBL-based WGCNA module membership audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Required objects
# ------------------------------------------------------------

required_objects <- c(
  "candidate_30",
  "candidate_master",
  "candidate_map_raw",
  "moduleColors",
  "module_genes",
  "core_step17"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    "Missing required object(s): ",
    paste(missing_objects, collapse = ", ")
  )
}

cat("All required objects are present.\n\n")


# ------------------------------------------------------------
# 1. Lock candidate order
# ------------------------------------------------------------

candidate_30_locked <- as.character(candidate_30)

cat("Number of locked candidates: ",
    length(candidate_30_locked), "\n", sep = "")

cat("Candidate order preserved: ",
    identical(
      candidate_30_locked,
      as.character(candidate_master$Symbol)
    ),
    "\n\n",
    sep = "")


# ------------------------------------------------------------
# 2. Basic integrity of moduleColors
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("2. moduleColors integrity\n")
cat("------------------------------------------------------------\n")

cat("Length(moduleColors): ",
    length(moduleColors), "\n", sep = "")

cat("Number of unique module colors: ",
    length(unique(unname(moduleColors))),
    "\n", sep = "")

cat("names(moduleColors) is NULL: ",
    is.null(names(moduleColors)),
    "\n", sep = "")

cat("Number of moduleColors names: ",
    length(names(moduleColors)),
    "\n\n", sep = "")

if (is.null(names(moduleColors))) {
  stop(
    "moduleColors has no names. ENSEMBL-based mapping cannot be audited."
  )
}

cat("First 10 moduleColors names:\n")
print(head(names(moduleColors), 10))

cat("\nFirst 10 moduleColors values:\n")
print(head(unname(moduleColors), 10))

cat("\n")


# ------------------------------------------------------------
# 3. Inspect candidate_map_raw structure
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("3. candidate_map_raw structure\n")
cat("------------------------------------------------------------\n")

cat("Dimensions: ",
    paste(dim(candidate_map_raw), collapse = " x "),
    "\n", sep = "")

cat("Columns:\n")
print(names(candidate_map_raw))

cat("\n")


# ------------------------------------------------------------
# 4. Restrict mapping audit to the 30 locked candidates
# ------------------------------------------------------------

candidate_map_30 <- candidate_map_raw[
  candidate_map_raw$SYMBOL %in% candidate_30_locked,
  ,
  drop = FALSE
]

cat("Rows in candidate_map_raw corresponding to locked candidates: ",
    nrow(candidate_map_30),
    "\n\n",
    sep = "")

print(candidate_map_30)


# ------------------------------------------------------------
# 5. Audit duplicate SYMBOL mappings
#
# IMPORTANT:
# Do NOT silently use match() before understanding duplicates.
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("5. Duplicate SYMBOL mapping audit\n")
cat("------------------------------------------------------------\n")

candidate_map_30$SYMBOL <- as.character(candidate_map_30$SYMBOL)
candidate_map_30$ENSEMBL <- as.character(candidate_map_30$ENSEMBL)

dup_symbols <- unique(
  candidate_map_30$SYMBOL[
    duplicated(candidate_map_30$SYMBOL) |
      duplicated(candidate_map_30$SYMBOL, fromLast = TRUE)
  ]
)

if (length(dup_symbols) == 0) {
  
  cat("No duplicated candidate SYMBOL mappings detected.\n\n")
  
} else {
  
  cat("Duplicated candidate SYMBOL(s):\n")
  print(dup_symbols)
  
  cat("\nDetailed duplicate mapping rows:\n")
  
  duplicate_mapping_audit <- candidate_map_30[
    candidate_map_30$SYMBOL %in% dup_symbols,
    c("SYMBOL", "ENSEMBL", "In_WGCNA"),
    drop = FALSE
  ]
  
  print(duplicate_mapping_audit)
  
  cat("\n")
}


# ------------------------------------------------------------
# 6. Identify candidates with missing SYMBOL -> ENSEMBL mapping
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("6. Missing SYMBOL -> ENSEMBL mapping audit\n")
cat("------------------------------------------------------------\n")

mapped_symbols <- unique(
  candidate_map_30$SYMBOL[
    !is.na(candidate_map_30$ENSEMBL) &
      candidate_map_30$ENSEMBL != ""
  ]
)

missing_symbol_mapping <- setdiff(
  candidate_30_locked,
  mapped_symbols
)

cat("Candidates with available ENSEMBL mapping: ",
    length(mapped_symbols),
    "/",
    length(candidate_30_locked),
    "\n", sep = "")

cat("Candidates without available ENSEMBL mapping: ",
    length(missing_symbol_mapping),
    "\n", sep = "")

if (length(missing_symbol_mapping) > 0) {
  print(missing_symbol_mapping)
}

cat("\n")


# ------------------------------------------------------------
# 7. Build candidate -> ENSEMBL mapping
#
# For unique SYMBOLs:
#   use the single available ENSEMBL.
#
# For duplicated SYMBOLs:
#   keep all mappings for audit rather than silently collapsing.
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("7. Candidate SYMBOL -> ENSEMBL mapping\n")
cat("------------------------------------------------------------\n")

mapping_list <- lapply(
  candidate_30_locked,
  function(sym) {
    
    x <- candidate_map_30[
      candidate_map_30$SYMBOL == sym,
      "ENSEMBL",
      drop = TRUE
    ]
    
    x <- unique(
      x[
        !is.na(x) &
          x != ""
      ]
    )
    
    x
  }
)

names(mapping_list) <- candidate_30_locked

mapping_summary <- data.frame(
  Symbol = candidate_30_locked,
  ENSEMBL = vapply(
    mapping_list,
    function(x) {
      if (length(x) == 0) {
        NA_character_
      } else if (length(x) == 1) {
        x
      } else {
        paste(x, collapse = ";")
      }
    },
    character(1)
  ),
  ENSEMBL_n = vapply(
    mapping_list,
    length,
    integer(1)
  ),
  stringsAsFactors = FALSE
)

print(mapping_summary)

cat("\n")


# ------------------------------------------------------------
# 8. Check ENSEMBL IDs against names(moduleColors)
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("8. ENSEMBL -> moduleColors membership audit\n")
cat("------------------------------------------------------------\n")

moduleColors_ids <- names(moduleColors)

candidate_module_audit <- data.frame(
  Symbol = candidate_30_locked,
  ENSEMBL = NA_character_,
  ENSEMBL_in_moduleColors = NA,
  WGCNA_Module = NA_character_,
  stringsAsFactors = FALSE
)


for (i in seq_along(candidate_30_locked)) {
  
  sym <- candidate_30_locked[i]
  ids <- mapping_list[[sym]]
  
  if (length(ids) == 0) {
    next
  }
  
  # IDs represented in moduleColors
  ids_present <- ids[ids %in% moduleColors_ids]
  
  if (length(ids_present) == 0) {
    
    candidate_module_audit$ENSEMBL[i] <-
      paste(ids, collapse = ";")
    
    candidate_module_audit$ENSEMBL_in_moduleColors[i] <-
      FALSE
    
  } else {
    
    candidate_module_audit$ENSEMBL[i] <-
      paste(ids_present, collapse = ";")
    
    candidate_module_audit$ENSEMBL_in_moduleColors[i] <-
      TRUE
    
    modules_found <- unique(
      unname(
        moduleColors[ids_present]
      )
    )
    
    candidate_module_audit$WGCNA_Module[i] <-
      paste(modules_found, collapse = ";")
  }
}

print(candidate_module_audit)

cat("\n")


# ------------------------------------------------------------
# 9. Summary of ENSEMBL representation
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("9. ENSEMBL representation summary\n")
cat("------------------------------------------------------------\n")

cat(
  "Candidates with ENSEMBL mapping: ",
  sum(!is.na(candidate_module_audit$ENSEMBL)),
  "/",
  nrow(candidate_module_audit),
  "\n",
  sep = ""
)

cat(
  "Candidates represented in moduleColors: ",
  sum(
    candidate_module_audit$ENSEMBL_in_moduleColors %in% TRUE,
    na.rm = TRUE
  ),
  "/",
  nrow(candidate_module_audit),
  "\n",
  sep = ""
)

cat(
  "Candidates not represented in moduleColors after mapping: ",
  sum(
    candidate_module_audit$ENSEMBL_in_moduleColors %in% FALSE,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat("\n")


# ------------------------------------------------------------
# 10. Module size audit directly from moduleColors
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("10. Module size audit from original moduleColors\n")
cat("------------------------------------------------------------\n")

module_size_from_moduleColors <- sort(
  table(unname(moduleColors)),
  decreasing = TRUE
)

module_size_df <- data.frame(
  Module = names(module_size_from_moduleColors),
  Module_Size = as.integer(module_size_from_moduleColors),
  row.names = NULL,
  stringsAsFactors = FALSE
)

cat("Blue module size:\n")
print(
  module_size_df[
    module_size_df$Module == "blue",
    ,
    drop = FALSE
  ]
)

cat("\nMagenta module size:\n")
print(
  module_size_df[
    module_size_df$Module == "magenta",
    ,
    drop = FALSE
  ]
)

cat("\n")


# ------------------------------------------------------------
# 11. Explicit core-gene recovery audit
#
# This is the key test:
#   AOC3  -> magenta
#   MS4A7 -> blue
#   IL7R  -> blue
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("11. Core-gene module recovery\n")
cat("------------------------------------------------------------\n")

core_genes_locked <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_recovery <- candidate_module_audit[
  match(core_genes_locked, candidate_module_audit$Symbol),
  ,
  drop = FALSE
]

print(core_recovery)

cat("\n")


# ------------------------------------------------------------
# 12. Compare recovered core modules with locked core_step17
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("12. Comparison with locked core_step17\n")
cat("------------------------------------------------------------")

cat("\n")

cat("core_step17 columns:\n")
print(names(core_step17))

cat("\nLocked core_step17 records:\n")
print(
  core_step17[
    core_step17$Symbol %in% core_genes_locked,
    ,
    drop = FALSE
  ]
)

cat("\n")


# Detect the module column in core_step17
possible_module_columns <- c(
  "TumorOnly_Module",
  "Module",
  "module",
  "WGCNA_Module"
)

core_module_column <- possible_module_columns[
  possible_module_columns %in% names(core_step17)
]

if (length(core_module_column) >= 1) {
  
  core_module_column <- core_module_column[1]
  
  locked_core_modules <- core_step17[
    match(core_genes_locked, core_step17$Symbol),
    c("Symbol", core_module_column),
    drop = FALSE
  ]
  
  names(locked_core_modules)[2] <- "Locked_Module"
  
  core_module_comparison <- merge(
    core_recovery[c("Symbol", "WGCNA_Module"),drop = FALSE],
    locked_core_modules,
    by = "Symbol",
    all = TRUE,
    sort = FALSE
  )
  
  core_module_comparison$Module_Match <-
    core_module_comparison$WGCNA_Module ==
    core_module_comparison$Locked_Module
  
  cat("Core module recovery comparison:\n")
  print(core_module_comparison)
  
} else {
  
  cat(
    "No recognizable module column found in core_step17.\n",
    "Module comparison will therefore use the printed locked records only.\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 13. Reconstruct blue module directly from moduleColors
#
# Purpose:
#   Test whether module_genes (1798 ENSEMBL IDs) is actually
#   the blue-module gene list.
#
# No WGCNA calculation is performed.
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("13. module_genes identity audit\n")
cat("------------------------------------------------------------\n")

blue_ids_from_moduleColors <- names(
  moduleColors[
    unname(moduleColors) == "blue"
  ]
)

module_genes_ids <- as.character(module_genes)

cat(
  "Blue module size from moduleColors: ",
  length(blue_ids_from_moduleColors),
  "\n",
  sep = ""
)

cat(
  "Length(module_genes): ",
  length(module_genes_ids),
  "\n",
  sep = ""
)

cat(
  "module_genes all contained in blue module: ",
  all(module_genes_ids %in% blue_ids_from_moduleColors),
  "\n",
  sep = ""
)

cat(
  "All blue-module IDs contained in module_genes: ",
  all(blue_ids_from_moduleColors %in% module_genes_ids),
  "\n",
  sep = ""
)

cat(
  "Exact same blue-module ENSEMBL set: ",
  setequal(
    module_genes_ids,
    blue_ids_from_moduleColors
  ),
  "\n",
  sep = ""
)

cat("\n")


# ------------------------------------------------------------
# 14. Blue and magenta module membership of core genes
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("14. Explicit core-gene membership checks\n")
cat("------------------------------------------------------------\n")

core_expected <- data.frame(
  Symbol = core_genes_locked,
  Expected_Module = c(
    "magenta",
    "blue",
    "blue"
  ),
  stringsAsFactors = FALSE
)

core_expected$Recovered_Module <-
  candidate_module_audit$WGCNA_Module[
    match(
      core_expected$Symbol,
      candidate_module_audit$Symbol
    )
  ]

core_expected$Match <-
  core_expected$Expected_Module ==
  core_expected$Recovered_Module

print(core_expected)

cat("\n")


# ------------------------------------------------------------
# 15. AOC3 / MS4A7 / IL7R ENSEMBL IDs and moduleColors values
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("15. Direct original-object lookup for core genes\n")
cat("------------------------------------------------------------\n")

core_direct_lookup <- do.call(
  rbind,
  lapply(
    core_genes_locked,
    function(sym) {
      
      ids <- mapping_list[[sym]]
      
      if (length(ids) == 0) {
        
        return(
          data.frame(
            Symbol = sym,
            ENSEMBL = NA_character_,
            Module = NA_character_,
            stringsAsFactors = FALSE
          )
        )
        
      }
      
      ids_present <- ids[
        ids %in% names(moduleColors)
      ]
      
      if (length(ids_present) == 0) {
        
        return(
          data.frame(
            Symbol = sym,
            ENSEMBL = paste(ids, collapse = ";"),
            Module = NA_character_,
            stringsAsFactors = FALSE
          )
        )
        
      }
      
      data.frame(
        Symbol = sym,
        ENSEMBL = ids_present,
        Module = unname(
          moduleColors[ids_present]
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

print(core_direct_lookup)

cat("\n")


# ------------------------------------------------------------
# 16. Candidate-level module distribution
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("16. Candidate module distribution\n")
cat("------------------------------------------------------------\n")

candidate_module_distribution <- table(
  candidate_module_audit$WGCNA_Module,
  useNA = "ifany"
)

print(candidate_module_distribution)

cat("\n")


# ------------------------------------------------------------
# 17. Final audit status
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP 18C FINAL AUDIT SUMMARY\n")
cat("============================================================\n")

n_candidates <- length(candidate_30_locked)

n_symbol_mapped <- sum(
  !is.na(candidate_module_audit$ENSEMBL)
)

n_in_moduleColors <- sum(
  candidate_module_audit$ENSEMBL_in_moduleColors %in% TRUE,
  na.rm = TRUE
)

core_module_ok <- all(
  core_expected$Match,
  na.rm = FALSE
)

blue_module_identity_ok <- setequal(
  module_genes_ids,
  blue_ids_from_moduleColors
)

cat("Locked candidates: ",
    n_candidates,
    "\n", sep = "")

cat("Candidates with SYMBOL -> ENSEMBL mapping: ",
    n_symbol_mapped,
    "/",
    n_candidates,
    "\n", sep = "")

cat("Candidates represented in moduleColors: ",
    n_in_moduleColors,
    "/",
    n_candidates,
    "\n", sep = "")

cat("Core genes recovered as AOC3=magenta, MS4A7=blue, IL7R=blue: ",
    core_module_ok,
    "\n", sep = "")

cat("module_genes exactly matches blue module from moduleColors: ",
    blue_module_identity_ok,
    "\n", sep = "")

cat("\n")

if (
  core_module_ok &&
  blue_module_identity_ok
) {
  
  cat(
    "STATUS: PASS\n",
    "The original module assignments of the three core genes can be\n",
    "independently recovered through SYMBOL -> ENSEMBL -> moduleColors.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "STATUS: REVIEW REQUIRED\n",
    "At least one core module recovery or module_genes identity check\n",
    "did not pass. Inspect the detailed audit tables above.\n",
    sep = ""
  )
}

cat("\n============================================================\n")
cat("STEP 18C completed. No WGCNA statistics were recalculated.\n")
cat("============================================================\n")

# ============================================================
# Paper 1
# STEP 18D
# Pure WGCNA / TME integration backfill and repair audit
#
# PURPOSE:
#   Safely integrate already-locked WGCNA module membership,
#   core-gene MM information, and TME association statistics
#   into a NEW 30-candidate audit table.
#
# IMPORTANT:
#   - PURE INTEGRATION / BACKFILL ONLY
#   - NO WGCNA recalculation
#   - NO module eigengene recalculation
#   - NO MM recalculation
#   - NO GS recalculation
#   - NO P-value recalculation
#   - NO FDR recalculation
#   - NO candidate reselection
#   - DO NOT modify candidate_master
#   - DO NOT overwrite locked objects
#
# DATA SOURCES:
#   1. candidate_master
#      -> locked candidate order
#      -> TME GS / P / FDR
#
#   2. candidate_module_audit
#      -> STEP 18C verified SYMBOL -> ENSEMBL -> moduleColors
#         module membership
#
#   3. core_step17
#      -> locked core-gene module size / MM / MM rank / percentile
#
# OUTPUT:
#   candidate_master_18D
#   candidate_master_18D_audit
# ============================================================


cat("\n============================================================\n")
cat("STEP 18D: Pure WGCNA / TME integration backfill\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Required objects
# ------------------------------------------------------------

required_objects <- c(
  "candidate_30",
  "candidate_master",
  "candidate_module_audit",
  "core_step17"
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
    "Missing required object(s): ",
    paste(missing_objects, collapse = ", ")
  )
  
}

cat("All required STEP 18D objects are present.\n\n")


# ------------------------------------------------------------
# 1. Preserve locked candidate order
# ------------------------------------------------------------

candidate_30_locked <- as.character(candidate_30)

cat("Locked candidate count: ",
    length(candidate_30_locked),
    "\n",
    sep = "")

candidate_master_order_ok <- identical(
  candidate_30_locked,
  as.character(candidate_master$Symbol)
)

cat(
  "candidate_master order matches candidate_30: ",
  candidate_master_order_ok,
  "\n\n",
  sep = ""
)

if (!candidate_master_order_ok) {
  
  stop(
    "Candidate order mismatch between candidate_30 and candidate_master."
  )
  
}


# ------------------------------------------------------------
# 2. Confirm STEP 18C module audit is internally consistent
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("2. STEP 18C module audit integrity check\n")
cat("------------------------------------------------------------\n")

required_module_columns <- c(
  "Symbol",
  "ENSEMBL",
  "ENSEMBL_in_moduleColors",
  "WGCNA_Module"
)

missing_module_columns <- setdiff(
  required_module_columns,
  names(candidate_module_audit)
)

if (length(missing_module_columns) > 0) {
  
  stop(
    "candidate_module_audit is missing required column(s): ",
    paste(missing_module_columns, collapse = ", ")
  )
  
}

candidate_module_audit_30 <- candidate_module_audit[
  match(
    candidate_30_locked,
    candidate_module_audit$Symbol
  ),
  required_module_columns,
  drop = FALSE
]

cat(
  "Rows recovered from candidate_module_audit: ",
  nrow(candidate_module_audit_30),
  "\n",
  sep = ""
)

cat(
  "All 30 candidate symbols recovered: ",
  identical(
    as.character(candidate_module_audit_30$Symbol),
    candidate_30_locked
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 3. Confirm core module assignments against locked core_step17
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("3. Core module consistency check\n")
cat("------------------------------------------------------------\n")

core_genes_locked <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_step17_core <- core_step17[
  core_step17$Symbol %in% core_genes_locked,
  ,
  drop = FALSE
]

core_module_compare <- data.frame(
  Symbol = core_genes_locked,
  Module_18C = candidate_module_audit_30$WGCNA_Module[
    match(
      core_genes_locked,
      candidate_module_audit_30$Symbol
    )
  ],
  Module_STEP17 = core_step17$TumorOnly_Module[
    match(
      core_genes_locked,
      core_step17$Symbol
    )
  ],
  stringsAsFactors = FALSE
)

core_module_compare$Module_Match <-
  core_module_compare$Module_18C ==
  core_module_compare$Module_STEP17

print(core_module_compare)

if (!all(core_module_compare$Module_Match)) {
  
  stop(
    "Core module assignment mismatch detected. ",
    "STEP 18D stopped before any integration."
  )
  
}

cat("\nCore module assignments are fully consistent.\n\n")


# ------------------------------------------------------------
# 4. Build WGCNA module backfill table
#
# IMPORTANT:
#   Module assignment comes ONLY from STEP 18C.
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("4. Construct module backfill table\n")
cat("------------------------------------------------------------\n")

module_backfill <- data.frame(
  Symbol = candidate_module_audit_30$Symbol,
  ENSEMBL = candidate_module_audit_30$ENSEMBL,
  TumorOnly_Module =
    candidate_module_audit_30$WGCNA_Module,
  Module_Source =
    rep(
      "STEP18C_moduleColors_audit",
      nrow(candidate_module_audit_30)
    ),
  stringsAsFactors = FALSE
)

module_backfill$Module_Mapping_Status <- ifelse(
  is.na(module_backfill$TumorOnly_Module),
  "Mapping_unavailable",
  "Recovered_from_moduleColors"
)

print(module_backfill)

cat("\n")


# ------------------------------------------------------------
# 5. Build locked core MM backfill table
#
# IMPORTANT:
#   MM is NOT recalculated.
#   Only existing locked STEP 17 values are imported.
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("5. Construct locked core MM backfill table\n")
cat("------------------------------------------------------------\n")

required_core_columns <- c(
  "Symbol",
  "TumorOnly_Module",
  "Module_Size",
  "MM",
  "MM_Rank",
  "MM_Percentile"
)

missing_core_columns <- setdiff(
  required_core_columns,
  names(core_step17)
)

if (length(missing_core_columns) > 0) {
  
  stop(
    "core_step17 is missing required column(s): ",
    paste(missing_core_columns, collapse = ", ")
  )
  
}

core_mm_backfill <- core_step17[
  core_step17$Symbol %in% core_genes_locked,
  required_core_columns,
  drop = FALSE
]

core_mm_backfill <- core_mm_backfill[
  match(
    core_genes_locked,
    core_mm_backfill$Symbol
  ),
  ,
  drop = FALSE
]

core_mm_backfill$MM_Source <-
  "Locked_STEP17"

print(core_mm_backfill)

cat("\n")


# ------------------------------------------------------------
# 6. Create complete 30-candidate integration table
#
# Start from candidate_master so all existing TME fields are
# preserved exactly.
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("6. Build candidate_master_18D\n")
cat("------------------------------------------------------------\n")

candidate_master_18D <- candidate_master

# ------------------------------------------------------------
# Add / replace module information in NEW object only
# ------------------------------------------------------------

candidate_master_18D$TumorOnly_Module <-
  module_backfill$TumorOnly_Module[
    match(
      candidate_master_18D$Symbol,
      module_backfill$Symbol
    )
  ]

candidate_master_18D$Module_Source <-
  module_backfill$Module_Source[
    match(
      candidate_master_18D$Symbol,
      module_backfill$Symbol
    )
  ]

candidate_master_18D$Module_Mapping_Status <-
  module_backfill$Module_Mapping_Status[
    match(
      candidate_master_18D$Symbol,
      module_backfill$Symbol
    )
  ]


# ------------------------------------------------------------
# Add ENSEMBL ID from STEP 18C
# ------------------------------------------------------------

candidate_master_18D$ENSEMBL_18D <-
  module_backfill$ENSEMBL[
    match(
      candidate_master_18D$Symbol,
      module_backfill$Symbol
    )
  ]


# ------------------------------------------------------------
# Add locked core MM information
#
# Non-core candidates remain NA.
# ------------------------------------------------------------

candidate_master_18D$MM_18D <- NA_real_
candidate_master_18D$MM_Rank_18D <- NA_integer_
candidate_master_18D$MM_Percentile_18D <- NA_real_
candidate_master_18D$MM_Source <- NA_character_

core_match <- match(
  candidate_master_18D$Symbol,
  core_mm_backfill$Symbol
)

has_locked_mm <- !is.na(core_match)

candidate_master_18D$MM_18D[has_locked_mm] <-
  core_mm_backfill$MM[core_match[has_locked_mm]]

candidate_master_18D$MM_Rank_18D[has_locked_mm] <-
  core_mm_backfill$MM_Rank[core_match[has_locked_mm]]

candidate_master_18D$MM_Percentile_18D[has_locked_mm] <-
  core_mm_backfill$MM_Percentile[core_match[has_locked_mm]]

candidate_master_18D$MM_Source[has_locked_mm] <-
  "Locked_STEP17"


# ------------------------------------------------------------
# 7. Add explicit integration status
# ------------------------------------------------------------

candidate_master_18D$Integration_Status <- "Integrated"

candidate_master_18D$Integration_Status[
  is.na(candidate_master_18D$TumorOnly_Module)
] <- "Module_mapping_unavailable"

candidate_master_18D$Integration_Status[
  candidate_master_18D$Symbol %in% core_genes_locked &
    is.na(candidate_master_18D$MM_18D)
] <- "Core_MM_missing"

cat(
  "candidate_master_18D dimensions: ",
  paste(dim(candidate_master_18D), collapse = " x "),
  "\n\n",
  sep = ""
)

print(candidate_master_18D)


# ------------------------------------------------------------
# 8. TME statistics preservation audit
#
# Verify that no TME association values changed.
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("8. TME statistics preservation audit\n")
cat("------------------------------------------------------------\n")

tme_columns <- c(
  "GS_Immune",
  "P_Immune",
  "FDR_Immune",
  "GS_Stromal",
  "P_Stromal",
  "FDR_Stromal",
  "GS_ESTIMATE",
  "P_ESTIMATE",
  "FDR_ESTIMATE"
)

missing_tme_columns <- setdiff(
  tme_columns,
  names(candidate_master)
)

if (length(missing_tme_columns) > 0) {
  
  stop(
    "candidate_master is missing TME column(s): ",
    paste(missing_tme_columns, collapse = ", ")
  )
  
}

tme_preservation_results <- data.frame(
  Column = tme_columns,
  Identical = vapply(
    tme_columns,
    function(x) {
      identical(
        candidate_master[[x]],
        candidate_master_18D[[x]]
      )
    },
    logical(1)
  ),
  stringsAsFactors = FALSE
)

print(tme_preservation_results)

if (!all(tme_preservation_results$Identical)) {
  
  stop(
    "One or more locked TME association columns changed."
  )
  
}

cat("\nAll locked TME association values preserved exactly.\n\n")


# ------------------------------------------------------------
# 9. Verify candidate order and count
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("9. Candidate order / count audit\n")
cat("------------------------------------------------------------\n")

candidate_order_final <- identical(
  as.character(candidate_master_18D$Symbol),
  candidate_30_locked
)

candidate_count_final <-
  nrow(candidate_master_18D) == 30

cat(
  "Exactly 30 candidates: ",
  candidate_count_final,
  "\n",
  sep = ""
)

cat(
  "Candidate order preserved: ",
  candidate_order_final,
  "\n",
  sep = ""
)

if (!candidate_count_final || !candidate_order_final) {
  
  stop(
    "Candidate count or order changed during integration."
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 10. Module recovery summary
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("10. Module recovery summary\n")
cat("------------------------------------------------------------\n")

module_recovery_summary <- data.frame(
  Metric = c(
    "Total candidates",
    "Module recovered",
    "Module unavailable",
    "Blue",
    "Brown",
    "Darkgrey",
    "Magenta"
  ),
  N = c(
    nrow(candidate_master_18D),
    sum(!is.na(candidate_master_18D$TumorOnly_Module)),
    sum(is.na(candidate_master_18D$TumorOnly_Module)),
    sum(candidate_master_18D$TumorOnly_Module == "blue",
        na.rm = TRUE),
    sum(candidate_master_18D$TumorOnly_Module == "brown",
        na.rm = TRUE),
    sum(candidate_master_18D$TumorOnly_Module == "darkgrey",
        na.rm = TRUE),
    sum(candidate_master_18D$TumorOnly_Module == "magenta",
        na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

print(module_recovery_summary)

cat("\n")


# ------------------------------------------------------------
# 11. Core-gene final integration audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("11. Core-gene final integration audit\n")
cat("------------------------------------------------------------\n")

core_final <- candidate_master_18D[
  match(
    core_genes_locked,
    candidate_master_18D$Symbol
  ),
  c(
    "Symbol",
    "ENSEMBL_18D",
    "TumorOnly_Module",
    "Module_Mapping_Status",
    "MM_18D",
    "MM_Rank_18D",
    "MM_Percentile_18D",
    "MM_Source",
    "GS_Immune",
    "FDR_Immune",
    "GS_Stromal",
    "FDR_Stromal",
    "GS_ESTIMATE",
    "FDR_ESTIMATE"
  ),
  drop = FALSE
]

print(core_final)

cat("\n")


# ------------------------------------------------------------
# 12. Verify core module and MM values against locked STEP 17
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("12. Core MM preservation audit\n")
cat("------------------------------------------------------------\n")

core_mm_compare <- data.frame(
  Symbol = core_genes_locked,
  
  Module_18D =
    core_final$TumorOnly_Module[
      match(
        core_genes_locked,
        core_final$Symbol
      )
    ],
  
  Module_STEP17 =
    core_step17$TumorOnly_Module[
      match(
        core_genes_locked,
        core_step17$Symbol
      )
    ],
  
  MM_18D =
    core_final$MM_18D[
      match(
        core_genes_locked,
        core_final$Symbol
      )
    ],
  
  MM_STEP17 =
    core_step17$MM[
      match(
        core_genes_locked,
        core_step17$Symbol
      )
    ],
  
  MM_Rank_18D =
    core_final$MM_Rank_18D[
      match(
        core_genes_locked,
        core_final$Symbol
      )
    ],
  
  MM_Rank_STEP17 =
    core_step17$MM_Rank[
      match(
        core_genes_locked,
        core_step17$Symbol
      )
    ],
  
  MM_Percentile_18D =
    core_final$MM_Percentile_18D[
      match(
        core_genes_locked,
        core_final$Symbol
      )
    ],
  
  MM_Percentile_STEP17 =
    core_step17$MM_Percentile[
      match(
        core_genes_locked,
        core_step17$Symbol
      )
    ],
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# Equality checks
# ------------------------------------------------------------

core_mm_compare$Module_Match <-
  core_mm_compare$Module_18D ==
  core_mm_compare$Module_STEP17

core_mm_compare$MM_Match <-
  abs(
    core_mm_compare$MM_18D -
      core_mm_compare$MM_STEP17
  ) < 1e-12

core_mm_compare$MM_Rank_Match <-
  core_mm_compare$MM_Rank_18D ==
  core_mm_compare$MM_Rank_STEP17

core_mm_compare$Percentile_Match <-
  abs(
    core_mm_compare$MM_Percentile_18D -
      core_mm_compare$MM_Percentile_STEP17
  ) < 1e-12


print(core_mm_compare)


# ------------------------------------------------------------
# Final audit
# ------------------------------------------------------------

if (!all(
  core_mm_compare$Module_Match,
  core_mm_compare$MM_Match,
  core_mm_compare$MM_Rank_Match,
  core_mm_compare$Percentile_Match
)) {
  
  stop(
    "Core STEP 17 values were not preserved."
  )
  
}

cat(
  "\nPASS: Core module, MM, MM rank, and MM percentile ",
  "match locked STEP 17 values.\n",
  sep = ""
)



# ------------------------------------------------------------
# 13. Explicit FAM189A2 audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("13. FAM189A2 mapping status audit\n")
cat("------------------------------------------------------------\n")

fam189a2_row <- candidate_master_18D[
  candidate_master_18D$Symbol == "FAM189A2",
  ,
  drop = FALSE
]

print(
  fam189a2_row[
    ,
    c(
      "Symbol",
      "ENSEMBL_18D",
      "TumorOnly_Module",
      "Module_Mapping_Status",
      "MM_18D",
      "MM_Source",
      "GS_Immune",
      "GS_Stromal",
      "GS_ESTIMATE"
    ),
    drop = FALSE
  ]
)

cat("\n")


# ------------------------------------------------------------
# 14. Confirm original candidate_master was NOT modified
# ------------------------------------------------------------
#
# We cannot prove historical object mutation without a saved
# pre-STEP-18D copy. However, we can confirm that STEP 18D
# did not assign into candidate_master in this script.
#
# The following checks are structural and informational only.
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("14. Original candidate_master structural audit\n")
cat("------------------------------------------------------------\n")

cat(
  "candidate_master dimensions: ",
  paste(dim(candidate_master), collapse = " x "),
  "\n",
  sep = ""
)

cat(
  "candidate_master still contains original WGCNA fields:\n"
)

print(
  intersect(
    c(
      "TumorOnly_Module",
      "TumorOnly_ME",
      "MM"
    ),
    names(candidate_master)
  )
)

cat("\n")


# ------------------------------------------------------------
# 15. Save a separate audit copy
#
# This does NOT overwrite candidate_master.
# ------------------------------------------------------------

candidate_master_18D_audit <- candidate_master_18D

cat("------------------------------------------------------------\n")
cat("15. STEP 18D output objects\n")
cat("------------------------------------------------------------\n")

cat(
  "Created object: candidate_master_18D\n"
)

cat(
  "Created object: candidate_master_18D_audit\n\n"
)

cat(
  "Original candidate_master has NOT been overwritten.\n\n"
)


# ------------------------------------------------------------
# 16. Final STEP 18D summary
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP 18D FINAL SUMMARY\n")
cat("============================================================\n")

cat(
  "30 candidates integrated: ",
  nrow(candidate_master_18D) == 30,
  "\n",
  sep = ""
)

cat(
  "Candidate order preserved: ",
  candidate_order_final,
  "\n",
  sep = ""
)

cat(
  "29 candidates have recovered WGCNA module membership: ",
  sum(!is.na(candidate_master_18D$TumorOnly_Module)) == 29,
  "\n",
  sep = ""
)

cat(
  "FAM189A2 remains unresolved due to missing mapping: ",
  is.na(
    candidate_master_18D$TumorOnly_Module[
      candidate_master_18D$Symbol == "FAM189A2"
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "All locked TME statistics preserved exactly: ",
  all(tme_preservation_results$Identical),
  "\n",
  sep = ""
)

cat(
  "All locked core module/MM values preserved exactly: ",
  all(
    core_mm_compare$Module_Match,
    core_mm_compare$MM_Match,
    core_mm_compare$MM_Rank_Match,
    core_mm_compare$Percentile_Match
  ),
  "\n",
  sep = ""
)

cat(
  "AOC3 = magenta: ",
  candidate_master_18D$TumorOnly_Module[
    candidate_master_18D$Symbol == "AOC3"
  ] == "magenta",
  "\n",
  sep = ""
)

cat(
  "MS4A7 = blue: ",
  candidate_master_18D$TumorOnly_Module[
    candidate_master_18D$Symbol == "MS4A7"
  ] == "blue",
  "\n",
  sep = ""
)

cat(
  "IL7R = blue: ",
  candidate_master_18D$TumorOnly_Module[
    candidate_master_18D$Symbol == "IL7R"
  ] == "blue",
  "\n",
  sep = ""
)

cat("\n")

if (
  nrow(candidate_master_18D) == 30 &&
  candidate_order_final &&
  sum(!is.na(candidate_master_18D$TumorOnly_Module)) == 29 &&
  all(tme_preservation_results$Identical) &&
  all(
    core_mm_compare$Module_Match,
    core_mm_compare$MM_Match,
    core_mm_compare$MM_Rank_Match,
    core_mm_compare$Percentile_Match
  ) &&
  candidate_master_18D$TumorOnly_Module[
    candidate_master_18D$Symbol == "AOC3"
  ] == "magenta" &&
  candidate_master_18D$TumorOnly_Module[
    candidate_master_18D$Symbol == "MS4A7"
  ] == "blue" &&
  candidate_master_18D$TumorOnly_Module[
    candidate_master_18D$Symbol == "IL7R"
  ] == "blue"
) {
  
  cat(
    "STATUS: PASS\n",
    "STEP 18D successfully created a new integrated 30-candidate\n",
    "WGCNA/TME audit table using only previously locked or audited\n",
    "values. No WGCNA or statistical recalculation was performed.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "STATUS: REVIEW REQUIRED\n",
    "One or more integration checks failed. No original locked\n",
    "analysis object should be replaced.\n",
    sep = ""
  )
}

cat("\n============================================================\n")
cat("STEP 18D completed.\n")
cat("============================================================\n")


# ============================================================
# Paper 1 Revision
# STEP 18F
# Pre-STEP19 Full Workspace Checkpoint
# Date: 
#
# Purpose:
#   Preserve the complete R workspace before starting STEP 19
#   and Supplementary Figure 4 plotting.
#
# IMPORTANT:
#   This step does NOT recalculate, overwrite, or modify
#   analytical results.
# ============================================================

cat("\n============================================================\n")
cat("STEP 18F: PRE-STEP19 FULL WORKSPACE CHECKPOINT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Create checkpoint directories
# ------------------------------------------------------------

checkpoint_root <- file.path(
  getwd(),
  "Paper1_STEP19_PREP_20260914"
)

dir.create(
  checkpoint_root,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(checkpoint_root, "R_workspace"),
  showWarnings = FALSE
)

dir.create(
  file.path(checkpoint_root, "scripts"),
  showWarnings = FALSE
)

dir.create(
  file.path(checkpoint_root, "locked_data"),
  showWarnings = FALSE
)

dir.create(
  file.path(checkpoint_root, "audit"),
  showWarnings = FALSE
)

dir.create(
  file.path(checkpoint_root, "README"),
  showWarnings = FALSE
)


cat("Checkpoint root:\n")
cat(checkpoint_root, "\n\n")


# ------------------------------------------------------------
# 2. Record current working directory
# ------------------------------------------------------------

cat("Current working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 3. Record R session information
# ------------------------------------------------------------

session_file <- file.path(
  checkpoint_root,
  "README",
  "R_sessionInfo.txt"
)

capture.output(
  sessionInfo(),
  file = session_file
)

cat("Saved sessionInfo():\n")
cat(session_file, "\n\n")


# ------------------------------------------------------------
# 4. Save COMPLETE current workspace
# ------------------------------------------------------------

workspace_file <- file.path(
  checkpoint_root,
  "R_workspace",
  "Paper1_PRE_STEP19_WORKSPACE_20260914.RData"
)

save.image(
  file = workspace_file
)

cat("Complete workspace saved:\n")
cat(workspace_file, "\n\n")


# ------------------------------------------------------------
# 5. Save current object inventory
# ------------------------------------------------------------

object_inventory <- data.frame(
  Object = ls(envir = .GlobalEnv),
  Class = sapply(
    ls(envir = .GlobalEnv),
    function(x) {
      paste(class(get(x, envir = .GlobalEnv)), collapse = ";")
    }
  ),
  Size_MB = sapply(
    ls(envir = .GlobalEnv),
    function(x) {
      round(
        as.numeric(
          object.size(get(x, envir = .GlobalEnv))
        ) / 1024^2,
        3
      )
    }
  ),
  stringsAsFactors = FALSE
)

object_inventory <- object_inventory[
  order(-object_inventory$Size_MB),
]

write.csv(
  object_inventory,
  file = file.path(
    checkpoint_root,
    "audit",
    "STEP18F_object_inventory_20260914.csv"
  ),
  row.names = FALSE
)

cat("Object inventory saved.\n\n")


# ------------------------------------------------------------
# 6. Record important locked objects
# ------------------------------------------------------------

locked_objects <- c(
  "module_tme_results",
  "GS_all",
  "MM_all",
  "moduleColors",
  "ME_color_formal_raw",
  "datExpr_MAD",
  "estimate_scores_wgcna",
  "module_genes",
  "core_visualization_table",
  "core_position_final",
  "core_descriptive_profile",
  "candidate_master_18D",
  "module_tme_lookup_save",
  "brown_survival_reference_save",
  "final_module_integration_save"
)

locked_status <- data.frame(
  Object = locked_objects,
  Exists = sapply(
    locked_objects,
    exists,
    envir = .GlobalEnv
  ),
  stringsAsFactors = FALSE
)

write.csv(
  locked_status,
  file = file.path(
    checkpoint_root,
    "audit",
    "STEP18F_locked_object_status_20260914.csv"
  ),
  row.names = FALSE
)

print(locked_status)

cat("\n")


# ------------------------------------------------------------
# 7. Basic dimensions of key locked objects
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("KEY OBJECT DIMENSIONS\n")
cat("------------------------------------------------------------\n")

if (exists("module_tme_results")) {
  cat(
    "module_tme_results:",
    paste(dim(module_tme_results), collapse = " x "),
    "\n"
  )
}

if (exists("GS_all")) {
  cat(
    "GS_all:",
    paste(dim(GS_all), collapse = " x "),
    "\n"
  )
}

if (exists("MM_all")) {
  cat(
    "MM_all:",
    paste(dim(MM_all), collapse = " x "),
    "\n"
  )
}

if (exists("moduleColors")) {
  cat(
    "moduleColors length:",
    length(moduleColors),
    "\n"
  )
}

if (exists("module_genes")) {
  cat(
    "module_genes length:",
    length(module_genes),
    "\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 8. Save key locked objects individually
# ------------------------------------------------------------

save_object_if_exists <- function(
    object_name,
    destination
) {
  
  if (exists(object_name, envir = .GlobalEnv)) {
    
    save(
      list = object_name,
      file = file.path(
        destination,
        paste0(
          object_name,
          "_PRE_STEP19_20260914.RData"
        )
      ),
      envir = .GlobalEnv
    )
    
    cat(
      "Saved:",
      object_name,
      "\n"
    )
    
  } else {
    
    cat(
      "NOT FOUND:",
      object_name,
      "\n"
    )
  }
}


for (obj in locked_objects) {
  
  save_object_if_exists(
    obj,
    file.path(
      checkpoint_root,
      "locked_data"
    )
  )
}


cat("\n============================================================\n")
cat("STEP 18F STATUS: CHECKPOINT SAVED\n")
cat("============================================================\n")

# ------------------------------------------------------------
# Copy today's analysis scripts
# ------------------------------------------------------------

script_files <- c(
  "Paper1_WORKFLOW03_14A.R",
  "Paper1_WORKFLOW03_14B.R",
  "Paper1_WORKFLOW03_14C.R",
  "Paper1_WORKFLOW03_14D.R",
  "Paper1_WORKFLOW03_14E.R",
  "Paper1_STEP18A.R",
  "Paper1_STEP18B.R",
  "Paper1_STEP18C.R",
  "Paper1_STEP18D.R"
)

for (f in script_files) {
  
  if (file.exists(f)) {
    
    file.copy(
      from = f,
      to = file.path(
        checkpoint_root,
        "scripts",
        basename(f)
      ),
      overwrite = TRUE
    )
    
    cat("Copied:", f, "\n")
    
  } else {
    
    cat("NOT FOUND:", f, "\n")
  }
}


# ============================================================
# STEP 18F-2
# Locked result summary
# ============================================================

summary_file <- file.path(
  checkpoint_root,
  "README",
  "STEP18F_locked_summary_20260914.txt"
)

sink(summary_file)

cat("============================================================\n")
cat("Paper 1 Revision\n")
cat("Pre-STEP19 Locked State Summary\n")
cat("Date: 2026-09-14\n")
cat("============================================================\n\n")

cat("Original Paper 1 core genes:\n")
cat("AOC3, MS4A7, IL7R\n\n")

cat("Tumor-only WGCNA:\n")
cat("Samples: 513 TCGA-LUAD tumors\n")
cat("Modules: 50\n")
cat("TME traits: ImmuneScore, StromalScore, ESTIMATEScore\n\n")

cat("Key modules:\n")
cat("ME2 / blue: ImmuneScore rho = 0.969\n")
cat("ME2 / blue: StromalScore rho = 0.708\n")
cat("ME2 / blue: ESTIMATEScore rho = 0.918\n")
cat("ME9 / magenta: ImmuneScore rho = 0.421\n")
cat("ME9 / magenta: StromalScore rho = 0.645\n")
cat("ME9 / magenta: ESTIMATEScore rho = 0.570\n\n")

cat("Core gene module positioning:\n")
cat("AOC3: ME9 / magenta, module size 534, MM rank 7, percentile 98.87430\n")
cat("IL7R: ME2 / blue, module size 1798, MM rank 120, percentile 93.37785\n")
cat("MS4A7: ME2 / blue, module size 1798, MM rank 267, percentile 85.19755\n\n")

cat("Workflow status:\n")
cat("STEP 18D candidate integration: PASS / LOCKED\n")
cat("Original 30-candidate WGCNA provenance: historical Sangerbox workflow\n")
cat("Tumor-only WGCNA: revision sensitivity/refinement analysis\n")
cat("STEP 19: NOT YET STARTED\n")
cat("Supplementary Figure 4 plotting: NOT YET STARTED\n\n")

cat("Important methodological rule:\n")
cat("Do not use tumor-only WGCNA to redefine the historical 30 candidates.\n")
cat("Do not recalculate locked GS/MM/TME statistics during plotting.\n")
cat("Supplementary Figure 4 is a visualization of locked results.\n")

sink()

cat(
  "Saved locked summary:\n",
  summary_file,
  "\n"
)

# ============================================================
# STEP 18F-3
# Checkpoint verification
# ============================================================

cat("\n============================================================\n")
cat("STEP 18F-3: CHECKPOINT VERIFICATION\n")
cat("============================================================\n\n")

required_files <- c(
  file.path(
    checkpoint_root,
    "R_workspace",
    "Paper1_PRE_STEP19_WORKSPACE_20260914.RData"
  ),
  file.path(
    checkpoint_root,
    "audit",
    "STEP18F_object_inventory_20260914.csv"
  ),
  file.path(
    checkpoint_root,
    "audit",
    "STEP18F_locked_object_status_20260914.csv"
  ),
  file.path(
    checkpoint_root,
    "README",
    "R_sessionInfo.txt"
  ),
  file.path(
    checkpoint_root,
    "README",
    "STEP18F_locked_summary_20260914.txt"
  )
)

verification <- data.frame(
  File = required_files,
  Exists = file.exists(required_files),
  Size_MB = round(
    file.info(required_files)$size / 1024^2,
    3
  ),
  stringsAsFactors = FALSE
)

print(verification)

if (all(verification$Exists)) {
  
  cat("\nPASS: All primary checkpoint files exist.\n")
  
} else {
  
  cat("\nWARNING: One or more checkpoint files are missing.\n")
}

cat("\n============================================================\n")










