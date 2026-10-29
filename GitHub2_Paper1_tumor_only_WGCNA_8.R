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

# Resume Workflow-04-2D WGCNA
# ============================================================

work_dir <- input_dir

setwd(work_dir)

getwd()

# ============================================================
# Paper 1
# Original 30 candidates × NEW tumor-only WGCNA integration
#
# Purpose:
#   1. Map original 30 candidates to tumor-only WGCNA modules
#   2. Extract MM / MM percentile
#   3. Extract GS for ImmuneScore / StromalScore / ESTIMATEScore
#   4. Extract P / FDR where available
#   5. Integrate original STRING/PPI degree
#   6. Generate a manuscript-ready integrated table
#
# IMPORTANT:
#   This script does NOT re-run WGCNA.
#   It uses the LOCKED objects generated in WORKFLOW03.
# ============================================================

# ------------------------------------------------------------
# 0. Load LOCKED tumor-only WGCNA results
# ------------------------------------------------------------

rdata_file <- "input_dir/Paper1_WORKFLOW03_14F_4D_2E_FINAL.RData"

load(rdata_file)

ls()
# ------------------------------------------------------------
# 1. Check that required objects exist
# ------------------------------------------------------------

required_objects <- c(
  "MM_all",
  "GS_all",
  "module_tme_results",
  "moduleColors",
  "core_position_final",
  "core_descriptive_profile",
  "core_visualization_table"
)

object_check <- data.frame(
  Object = required_objects,
  Exists = sapply(required_objects, exists),
  stringsAsFactors = FALSE
)

print(object_check)


# ------------------------------------------------------------
# 2. Inspect object structures
# ------------------------------------------------------------

cat("\n================ MM_all =================\n")
print(class(MM_all))
print(dim(MM_all))
print(head(MM_all))
print(colnames(MM_all))


cat("\n================ GS_all =================\n")
print(class(GS_all))
print(dim(GS_all))
print(head(GS_all))
print(colnames(GS_all))


cat("\n================ module_tme_results =================\n")
print(class(module_tme_results))
print(dim(module_tme_results))
print(head(module_tme_results))
print(colnames(module_tme_results))


cat("\n================ moduleColors =================\n")
print(class(moduleColors))
print(length(moduleColors))
print(head(moduleColors))


cat("\n================ core_position_final =================\n")
print(class(core_position_final))
print(dim(core_position_final))
print(head(core_position_final))
print(colnames(core_position_final))

### ============================================================
# Paper 1 - Tumor-only WGCNA
# Step 15A: Inspect remaining locked objects
# ============================================================

cat("\n============================================================\n")
cat("STEP 15A: Inspect remaining locked objects\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. core_descriptive_profile
# ------------------------------------------------------------

cat("\n================ core_descriptive_profile =================\n")

print(class(core_descriptive_profile))
print(dim(core_descriptive_profile))
print(colnames(core_descriptive_profile))
print(core_descriptive_profile)


# ------------------------------------------------------------
# 2. core_visualization_table
# ------------------------------------------------------------

cat("\n================ core_visualization_table =================\n")

print(class(core_visualization_table))
print(dim(core_visualization_table))
print(colnames(core_visualization_table))
print(head(core_visualization_table, 10))


# ------------------------------------------------------------
# 3. Check rowname consistency
# ------------------------------------------------------------

cat("\n================ Gene ID consistency =================\n")

cat("MM_all rownames:", nrow(MM_all), "\n")
cat("GS_all rownames:", nrow(GS_all), "\n")
cat("moduleColors length:", length(moduleColors), "\n")

cat(
  "MM_all vs GS_all rownames identical: ",
  identical(rownames(MM_all), rownames(GS_all)),
  "\n"
)

cat(
  "MM_all rownames vs GS_all$Gene identical: ",
  identical(rownames(MM_all), GS_all$Gene),
  "\n"
)


# ------------------------------------------------------------
# 4. Check module names
# ------------------------------------------------------------

cat("\n================ Module consistency =================\n")

module_names_MM <- sub("^MM", "", colnames(MM_all))

module_names_colors <- sort(unique(moduleColors))

cat("Number of MM modules:", length(module_names_MM), "\n")
cat("Number of moduleColors:", length(module_names_colors), "\n")

cat("\nModules in MM_all:\n")
print(module_names_MM)

cat("\nModules in moduleColors:\n")
print(module_names_colors)

cat(
  "\nMM modules and moduleColors modules identical: ",
  setequal(module_names_MM, module_names_colors),
  "\n"
)


# ============================================================
# Paper 1
# STEP 16A
# Original 30 candidate genes × tumor-only WGCNA
#
# Purpose:
#   Map the original 30 candidate genes onto the locked
#   tumor-only WGCNA framework.
#
# IMPORTANT:
#   This is an integration/audit analysis.
#   It does NOT redefine the three study genes.
#   It does NOT calculate a composite score.
# ============================================================


cat("\n============================================================\n")
cat("STEP 16A: Original 30 candidates × tumor-only WGCNA\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Original 30 candidate genes
# ------------------------------------------------------------

candidate_30 <- c(
  
  # Original BLUE module candidates
  "VIPR1",
  "CA4",
  "PLA2G4F",
  "SPARCL1",
  "GIMAP7",
  "GIMAP8",
  "ANOS1",
  "FAM189A2",
  "CD34",
  "CGNL1",
  "FHL1",
  "GPIHBP1",
  "FAM107A",
  "DPYSL2",
  "CLEC3B",
  "IL7R",
  "VEGFD",
  "GIMAP6",
  "MFAP4",
  "AGER",
  "TMEM204",
  "SFTPC",
  "CAT",
  "AOC3",
  
  # Original GREEN module candidates
  "MS4A7",
  "CD300LF",
  "GLIPR2",
  "CYP27A1",
  "DOK2",
  "CD52"
)


core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


stopifnot(length(candidate_30) == 30)
stopifnot(length(unique(candidate_30)) == 30)


# ------------------------------------------------------------
# 2. Check annotation packages
# ------------------------------------------------------------

if (!requireNamespace(
  "org.Hs.eg.db",
  quietly = TRUE
)) {
  stop(
    "Package 'org.Hs.eg.db' is required."
  )
}

if (!requireNamespace(
  "AnnotationDbi",
  quietly = TRUE
)) {
  stop(
    "Package 'AnnotationDbi' is required."
  )
}


library(AnnotationDbi)
library(org.Hs.eg.db)


# ------------------------------------------------------------
# 3. Symbol -> Ensembl mapping
# ------------------------------------------------------------

candidate_map_raw <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = candidate_30,
  keytype = "SYMBOL",
  columns = c(
    "SYMBOL",
    "ENSEMBL"
  )
)


candidate_map_raw <- candidate_map_raw[
  !is.na(candidate_map_raw$SYMBOL) &
    !is.na(candidate_map_raw$ENSEMBL),
  c(
    "SYMBOL",
    "ENSEMBL"
  )
]


candidate_map_raw <- unique(
  candidate_map_raw
)


# ------------------------------------------------------------
# 4. Mapping completeness
# ------------------------------------------------------------

mapped_symbols <- unique(
  candidate_map_raw$SYMBOL
)

unmapped_symbols <- setdiff(
  candidate_30,
  mapped_symbols
)


cat("\n================ Mapping QC =================\n")

cat(
  "Original candidates: ",
  length(candidate_30),
  "\n",
  sep = ""
)

cat(
  "Mapped symbols:      ",
  length(mapped_symbols),
  "\n",
  sep = ""
)

cat(
  "Unmapped symbols:    ",
  length(unmapped_symbols),
  "\n",
  sep = ""
)


if (length(unmapped_symbols) > 0) {
  
  cat("\nUnmapped symbols:\n")
  print(unmapped_symbols)
  
}


# ------------------------------------------------------------
# 5. One-to-many mapping audit
# ------------------------------------------------------------

mapping_count <- aggregate(
  ENSEMBL ~ SYMBOL,
  data = candidate_map_raw,
  FUN = function(x) {
    length(unique(x))
  }
)


multi_mapping <- mapping_count[
  mapping_count$ENSEMBL > 1,
]


cat(
  "\nSymbols with >1 Ensembl ID: ",
  nrow(multi_mapping),
  "\n",
  sep = ""
)


if (nrow(multi_mapping) > 0) {
  
  print(multi_mapping)
  
}


# ------------------------------------------------------------
# 6. Check whether Ensembl IDs are present in WGCNA
# ------------------------------------------------------------

candidate_map_raw$In_WGCNA <-
  candidate_map_raw$ENSEMBL %in%
  rownames(MM_all)


absent_wgcna <- unique(
  candidate_map_raw$ENSEMBL[
    !candidate_map_raw$In_WGCNA
  ]
)


cat(
  "\nEnsembl IDs absent from WGCNA: ",
  length(absent_wgcna),
  "\n",
  sep = ""
)


if (length(absent_wgcna) > 0) {
  
  print(absent_wgcna)
  
}


# ------------------------------------------------------------
# 7. Build candidate master table
# ------------------------------------------------------------

candidate_master <- data.frame(
  Candidate_Order = seq_along(candidate_30),
  Symbol = candidate_30,
  stringsAsFactors = FALSE
)


candidate_master$ENSEMBL <- candidate_map_raw$ENSEMBL[
  match(
    candidate_master$Symbol,
    candidate_map_raw$SYMBOL
  )
]


# ------------------------------------------------------------
# 8. Core-gene flag
# ------------------------------------------------------------

candidate_master$Core3 <-
  candidate_master$Symbol %in%
  core_genes


# ------------------------------------------------------------
# 9. Tumor-only WGCNA module
# ------------------------------------------------------------

candidate_master$TumorOnly_Module <-
  unname(
    moduleColors[
      match(
        candidate_master$ENSEMBL,
        names(moduleColors)
      )
    ]
  )


# ------------------------------------------------------------
# 10. Correct ME label
# ------------------------------------------------------------

candidate_master$TumorOnly_ME <- NA_character_


for (i in seq_len(nrow(candidate_master))) {
  
  module_i <-
    candidate_master$TumorOnly_Module[i]
  
  if (!is.na(module_i)) {
    
    mm_col <-
      paste0(
        "MM",
        module_i
      )
    
    if (mm_col %in% colnames(MM_all)) {
      
      candidate_master$TumorOnly_ME[i] <-
        sub(
          "^MM",
          "ME",
          mm_col
        )
      
    }
    
  }
  
}


# ------------------------------------------------------------
# 11. Module membership
# ------------------------------------------------------------

candidate_master$MM <- NA_real_


for (i in seq_len(nrow(candidate_master))) {
  
  gene_i <-
    candidate_master$ENSEMBL[i]
  
  module_i <-
    candidate_master$TumorOnly_Module[i]
  
  if (
    !is.na(gene_i) &&
    !is.na(module_i)
  ) {
    
    mm_col <-
      paste0(
        "MM",
        module_i
      )
    
    if (mm_col %in% colnames(MM_all)) {
      
      candidate_master$MM[i] <-
        MM_all[
          gene_i,
          mm_col
        ]
      
    }
    
  }
  
}


# ------------------------------------------------------------
# 12. Gene significance
# ------------------------------------------------------------

gs_index <- match(
  candidate_master$ENSEMBL,
  rownames(GS_all)
)


candidate_master$GS_Immune <-
  GS_all$GS_Immune[gs_index]

candidate_master$P_Immune <-
  GS_all$P_Immune[gs_index]

candidate_master$FDR_Immune <-
  GS_all$FDR_Immune[gs_index]


candidate_master$GS_Stromal <-
  GS_all$GS_Stromal[gs_index]

candidate_master$P_Stromal <-
  GS_all$P_Stromal[gs_index]

candidate_master$FDR_Stromal <-
  GS_all$FDR_Stromal[gs_index]


candidate_master$GS_ESTIMATE <-
  GS_all$GS_ESTIMATE[gs_index]

candidate_master$P_ESTIMATE <-
  GS_all$P_ESTIMATE[gs_index]

candidate_master$FDR_ESTIMATE <-
  GS_all$FDR_ESTIMATE[gs_index]


# ------------------------------------------------------------
# 13. Final column order
# ------------------------------------------------------------

candidate_master <- candidate_master[
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
    "FDR_ESTIMATE"
  )
]


# ------------------------------------------------------------
# 14. Final 30-gene completeness audit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL 30-GENE COMPLETENESS AUDIT\n")
cat("============================================================\n")


cat(
  "Rows in candidate_master: ",
  nrow(candidate_master),
  "\n",
  sep = ""
)

cat(
  "Unique symbols:           ",
  length(unique(candidate_master$Symbol)),
  "\n",
  sep = ""
)

cat(
  "Ensembl IDs present:      ",
  sum(!is.na(candidate_master$ENSEMBL)),
  "\n",
  sep = ""
)

cat(
  "Tumor-only modules found: ",
  sum(!is.na(candidate_master$TumorOnly_Module)),
  "\n",
  sep = ""
)

cat(
  "MM values found:          ",
  sum(!is.na(candidate_master$MM)),
  "\n",
  sep = ""
)

cat(
  "GS values found:          ",
  sum(!is.na(candidate_master$GS_Immune)),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 15. Print complete table
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ORIGINAL 30 CANDIDATES × TUMOR-ONLY WGCNA\n")
cat("============================================================\n")

print(
  candidate_master,
  row.names = FALSE
)


# ============================================================
# Paper 1
# STEP 16B
# QC: verify moduleColors ↔ MM_all gene identity
# ============================================================

cat("\n================ STEP 16B: moduleColors QC ================\n")

cat("\n1. Basic dimensions:\n")
cat("length(moduleColors) =", length(moduleColors), "\n")
cat("nrow(MM_all)         =", nrow(MM_all), "\n")

cat("\n2. Names of moduleColors:\n")
cat("is.null(names(moduleColors)) =", is.null(names(moduleColors)), "\n")
cat("length(names(moduleColors))  =", length(names(moduleColors)), "\n")

cat("\n3. First 10 moduleColors names:\n")
print(head(names(moduleColors), 10))

cat("\n4. First 10 MM_all rownames:\n")
print(head(rownames(MM_all), 10))

cat("\n5. Identity check:\n")

if (!is.null(names(moduleColors))) {
  
  cat(
    "identical(names(moduleColors), rownames(MM_all)) = ",
    identical(names(moduleColors), rownames(MM_all)),
    "\n",
    sep = ""
  )
  
  cat(
    "sum(names(moduleColors) == rownames(MM_all)) = ",
    sum(names(moduleColors) == rownames(MM_all)),
    "\n",
    sep = ""
  )
  
} else {
  
  cat("moduleColors currently has NO names.\n")
  
}

cat("\n6. Candidate gene matching test:\n")

test_symbols <- c("AOC3", "MS4A7", "IL7R")

test_ens <- candidate_master$ENSEMBL[
  match(test_symbols, candidate_master$Symbol)
]

print(
  data.frame(
    Symbol = test_symbols,
    ENSEMBL = test_ens,
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)

cat("\n7. Match against rownames(MM_all):\n")

print(
  data.frame(
    Symbol = test_symbols,
    ENSEMBL = test_ens,
    In_MM_all = test_ens %in% rownames(MM_all),
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)

cat("\n============================================================\n")



# ============================================================
# Paper 1
# STEP 16C
# Restore moduleColors gene identifiers
# and perform strict identity QC
# ============================================================

cat("\n================ STEP 16C: Restore moduleColors names ================\n")


# ------------------------------------------------------------
# 1. Basic structural QC
# ------------------------------------------------------------

stopifnot(
  length(moduleColors) == nrow(MM_all)
)

stopifnot(
  length(moduleColors) == nrow(GS_all)
)

cat("\n1. Dimension checks:\n")
cat("length(moduleColors) =", length(moduleColors), "\n")
cat("nrow(MM_all)         =", nrow(MM_all), "\n")
cat("nrow(GS_all)         =", nrow(GS_all), "\n")


# ------------------------------------------------------------
# 2. Restore gene identifiers
# ------------------------------------------------------------

if (is.null(names(moduleColors))) {
  
  names(moduleColors) <- rownames(MM_all)
  
  cat("\n2. moduleColors names were NULL.\n")
  cat("Names restored from rownames(MM_all).\n")
  
} else {
  
  cat("\n2. moduleColors already had names.\n")
  
}


# ------------------------------------------------------------
# 3. Strict identity checks
# ------------------------------------------------------------

cat("\n3. Identity checks:\n")

cat(
  "is.null(names(moduleColors)) = ",
  is.null(names(moduleColors)),
  "\n",
  sep = ""
)

cat(
  "length(names(moduleColors)) = ",
  length(names(moduleColors)),
  "\n",
  sep = ""
)

cat(
  "identical(names(moduleColors), rownames(MM_all)) = ",
  identical(names(moduleColors), rownames(MM_all)),
  "\n",
  sep = ""
)

cat(
  "identical(names(moduleColors), rownames(GS_all)) = ",
  identical(names(moduleColors), rownames(GS_all)),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Check duplicated / missing gene identifiers
# ------------------------------------------------------------

cat("\n4. Gene identifier integrity:\n")

cat(
  "Duplicated moduleColors names = ",
  sum(duplicated(names(moduleColors))),
  "\n",
  sep = ""
)

cat(
  "Missing moduleColors names = ",
  sum(is.na(names(moduleColors)) | names(moduleColors) == ""),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. Check module set consistency
# ------------------------------------------------------------

module_set_colors <- sort(unique(moduleColors))

module_set_MM <- sort(
  sub("^MM", "", colnames(MM_all))
)

cat("\n5. Module set consistency:\n")

cat(
  "Number of modules in moduleColors = ",
  length(module_set_colors),
  "\n",
  sep = ""
)

cat(
  "Number of MM module columns       = ",
  length(module_set_MM),
  "\n",
  sep = ""
)

cat(
  "identical(module sets) = ",
  identical(module_set_colors, module_set_MM),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Candidate-gene module mapping test
# ------------------------------------------------------------

test_symbols <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

test_ens <- candidate_master$ENSEMBL[
  match(test_symbols, candidate_master$Symbol)
]

test_module <- moduleColors[
  match(test_ens, names(moduleColors))
]

test_table <- data.frame(
  Symbol = test_symbols,
  ENSEMBL = test_ens,
  TumorOnly_Module = unname(test_module),
  stringsAsFactors = FALSE
)

cat("\n6. Core gene module mapping:\n")

print(
  test_table,
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Expected locked module check
# ------------------------------------------------------------

expected_modules <- c(
  AOC3 = "magenta",
  MS4A7 = "blue",
  IL7R = "blue"
)

observed_modules <- setNames(
  test_table$TumorOnly_Module,
  test_table$Symbol
)

cat("\n7. Expected locked module check:\n")

print(
  data.frame(
    Symbol = names(expected_modules),
    Expected_Module = unname(expected_modules),
    Observed_Module = unname(
      observed_modules[names(expected_modules)]
    ),
    Match = unname(
      observed_modules[names(expected_modules)] ==
        expected_modules
    ),
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Final hard-stop QC
# ------------------------------------------------------------

stopifnot(
  identical(names(moduleColors), rownames(MM_all))
)

stopifnot(
  identical(names(moduleColors), rownames(GS_all))
)

stopifnot(
  identical(module_set_colors, module_set_MM)
)

stopifnot(
  all(
    observed_modules[names(expected_modules)] ==
      expected_modules
  )
)

cat("\n================ STEP 16C QC PASSED ================\n")
cat("moduleColors is now correctly indexed by Ensembl ID.\n")
cat("AOC3  -> magenta / ME9\n")
cat("MS4A7 -> blue    / ME2\n")
cat("IL7R  -> blue    / ME2\n")
cat("=====================================================\n")


# ============================================================
# Paper 1
# STEP 16D
# Rebuild the 30-candidate tumor-only WGCNA master table
# after restoring moduleColors gene identifiers
# ============================================================

cat("\n================ STEP 16D: Rebuild candidate master ================\n")


# ------------------------------------------------------------
# 1. Fixed candidate set
# ------------------------------------------------------------

candidate_30 <- c(
  "VIPR1","CA4","PLA2G4F","SPARCL1","GIMAP7","GIMAP8","ANOS1",
  "FAM189A2","CD34","CGNL1","FHL1","GPIHBP1","FAM107A","DPYSL2",
  "CLEC3B","IL7R","VEGFD","GIMAP6","MFAP4","AGER","TMEM204",
  "SFTPC","CAT","AOC3","MS4A7","CD300LF","GLIPR2","CYP27A1",
  "DOK2","CD52"
)

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

stopifnot(length(candidate_30) == 30)
stopifnot(length(unique(candidate_30)) == 30)


# ------------------------------------------------------------
# 2. Original WGCNA module assignment
# ------------------------------------------------------------

original_module <- c(
  rep("blue", 24),
  rep("green", 6)
)

stopifnot(length(original_module) == 30)


# ------------------------------------------------------------
# 3. Build base candidate table
# ------------------------------------------------------------

candidate_master_16D <- data.frame(
  Candidate_Order = seq_along(candidate_30),
  Symbol = candidate_30,
  Original_Module = original_module,
  Core3 = candidate_30 %in% core_genes,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 4. Add Ensembl mapping
# ------------------------------------------------------------

candidate_master_16D$ENSEMBL <-
  candidate_master$ENSEMBL[
    match(
      candidate_master_16D$Symbol,
      candidate_master$Symbol
    )
  ]


# ------------------------------------------------------------
# 5. Map tumor-only WGCNA module
# ------------------------------------------------------------

candidate_master_16D$TumorOnly_Module <-
  unname(
    moduleColors[
      match(
        candidate_master_16D$ENSEMBL,
        names(moduleColors)
      )
    ]
  )


# ------------------------------------------------------------
# 6. Derive correct ME label
# ------------------------------------------------------------

candidate_master_16D$TumorOnly_ME <- NA_character_

for (i in seq_len(nrow(candidate_master_16D))) {
  
  module_i <- candidate_master_16D$TumorOnly_Module[i]
  
  if (!is.na(module_i)) {
    
    mm_col <- paste0("MM", module_i)
    
    if (mm_col %in% colnames(MM_all)) {
      
      candidate_master_16D$TumorOnly_ME[i] <-
        sub("^MM", "ME", mm_col)
    }
  }
}


# ------------------------------------------------------------
# 7. Extract module membership (MM)
# ------------------------------------------------------------

candidate_master_16D$MM <- NA_real_

for (i in seq_len(nrow(candidate_master_16D))) {
  
  gene_i <- candidate_master_16D$ENSEMBL[i]
  module_i <- candidate_master_16D$TumorOnly_Module[i]
  
  if (
    !is.na(gene_i) &&
    !is.na(module_i)
  ) {
    
    mm_col <- paste0("MM", module_i)
    
    if (mm_col %in% colnames(MM_all)) {
      
      candidate_master_16D$MM[i] <-
        MM_all[gene_i, mm_col]
    }
  }
}


# ------------------------------------------------------------
# 8. Extract gene significance for all three TME traits
# ------------------------------------------------------------

gs_index <- match(
  candidate_master_16D$ENSEMBL,
  rownames(GS_all)
)

candidate_master_16D$GS_Immune <-
  GS_all$GS_Immune[gs_index]

candidate_master_16D$P_Immune <-
  GS_all$P_Immune[gs_index]

candidate_master_16D$FDR_Immune <-
  GS_all$FDR_Immune[gs_index]

candidate_master_16D$GS_Stromal <-
  GS_all$GS_Stromal[gs_index]

candidate_master_16D$P_Stromal <-
  GS_all$P_Stromal[gs_index]

candidate_master_16D$FDR_Stromal <-
  GS_all$FDR_Stromal[gs_index]

candidate_master_16D$GS_ESTIMATE <-
  GS_all$GS_ESTIMATE[gs_index]

candidate_master_16D$P_ESTIMATE <-
  GS_all$P_ESTIMATE[gs_index]

candidate_master_16D$FDR_ESTIMATE <-
  GS_all$FDR_ESTIMATE[gs_index]


# ------------------------------------------------------------
# 9. Reorder columns
# ------------------------------------------------------------

candidate_master_16D <- candidate_master_16D[
  ,
  c(
    "Candidate_Order",
    "Symbol",
    "ENSEMBL",
    "Core3",
    "Original_Module",
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
]


# ------------------------------------------------------------
# 10. Print complete table
# ------------------------------------------------------------

cat("\n10. Rebuilt 30-candidate master table:\n\n")

print(
  candidate_master_16D,
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Basic QC
# ------------------------------------------------------------

cat("\n11. QC summary:\n")

cat(
  "Number of candidates = ",
  nrow(candidate_master_16D),
  "\n",
  sep = ""
)

cat(
  "Unique symbols = ",
  length(unique(candidate_master_16D$Symbol)),
  "\n",
  sep = ""
)

cat(
  "Ensembl mapped = ",
  sum(!is.na(candidate_master_16D$ENSEMBL)),
  "\n",
  sep = ""
)

cat(
  "Tumor-only module mapped = ",
  sum(!is.na(candidate_master_16D$TumorOnly_Module)),
  "\n",
  sep = ""
)

cat(
  "Tumor-only ME mapped = ",
  sum(!is.na(candidate_master_16D$TumorOnly_ME)),
  "\n",
  sep = ""
)

cat(
  "MM available = ",
  sum(!is.na(candidate_master_16D$MM)),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Core gene QC
# ------------------------------------------------------------

core_check_16D <- candidate_master_16D[
  candidate_master_16D$Core3,
  c(
    "Symbol",
    "ENSEMBL",
    "Original_Module",
    "TumorOnly_Module",
    "TumorOnly_ME",
    "MM",
    "GS_Immune",
    "GS_Stromal",
    "GS_ESTIMATE"
  )
]

cat("\n12. Core three genes:\n\n")

print(
  core_check_16D,
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Hard-stop checks
# ------------------------------------------------------------

stopifnot(
  nrow(candidate_master_16D) == 30
)

stopifnot(
  length(unique(candidate_master_16D$Symbol)) == 30
)

stopifnot(
  candidate_master_16D$TumorOnly_Module[
    match("AOC3", candidate_master_16D$Symbol)
  ] == "magenta"
)

stopifnot(
  candidate_master_16D$TumorOnly_Module[
    match("MS4A7", candidate_master_16D$Symbol)
  ] == "blue"
)

stopifnot(
  candidate_master_16D$TumorOnly_Module[
    match("IL7R", candidate_master_16D$Symbol)
  ] == "blue"
)

cat("\n================ STEP 16D COMPLETED ================\n")
cat("30-candidate tumor-only WGCNA master table rebuilt.\n")
cat("=====================================================\n")


# ============================================================
# Paper 1
# STEP 16E
# Validate / establish formal tumor-only module labels
# ============================================================

cat("\n================ STEP 16E: Module label QC ================\n")


# ------------------------------------------------------------
# 1. Inspect locked core_position_final
# ------------------------------------------------------------

cat("\n1. Locked core_position_final:\n")

print(
  core_position_final[
    ,
    c(
      "Symbol",
      "Module",
      "ME",
      "Module_Size",
      "MM",
      "MM_percentile",
      "MM_rank"
    )
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 2. Build validated color → ME mapping for the three
#    prioritized study genes
# ------------------------------------------------------------

core_me_mapping <- unique(
  core_position_final[
    ,
    c("Module", "ME")
  ]
)

cat("\n2. Validated module color → ME mapping:\n")

print(
  core_me_mapping,
  row.names = FALSE
)


# ------------------------------------------------------------
# 3. Check expected mappings
# ------------------------------------------------------------

expected_core_mapping <- data.frame(
  Module = c("magenta", "blue"),
  ME = c("ME9", "ME2"),
  stringsAsFactors = FALSE
)

cat("\n3. Expected core mappings:\n")

print(
  expected_core_mapping,
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. Hard-stop validation
# ------------------------------------------------------------

stopifnot(
  core_position_final$ME[
    match("AOC3", core_position_final$Symbol)
  ] == "ME9"
)

stopifnot(
  core_position_final$ME[
    match("MS4A7", core_position_final$Symbol)
  ] == "ME2"
)

stopifnot(
  core_position_final$ME[
    match("IL7R", core_position_final$Symbol)
  ] == "ME2"
)

cat("\n================ STEP 16E QC PASSED ================\n")
cat("AOC3  -> magenta / ME9\n")
cat("MS4A7 -> blue    / ME2\n")
cat("IL7R  -> blue    / ME2\n")
cat("=====================================================\n")


# ============================================================
# Paper 1
# STEP 17
# Module-specific MM rank and percentile for the 30 candidates
# ============================================================

cat("\n============================================================\n")
cat("STEP 17: Module-specific MM rank / percentile\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Basic QC
# ------------------------------------------------------------

cat("\n1. Basic QC\n")

stopifnot(
  identical(rownames(MM_all), rownames(GS_all))
)

stopifnot(
  identical(rownames(MM_all), names(moduleColors))
)

cat("MM_all rows:", nrow(MM_all), "\n")
cat("Number of modules:", ncol(MM_all), "\n")


# ------------------------------------------------------------
# 2. Extract candidate-level module and MM
# ------------------------------------------------------------

step17_df <- candidate_master

step17_df$TumorOnly_Module <- NA_character_
step17_df$MM <- NA_real_


for (i in seq_len(nrow(step17_df))) {
  
  gene_id <- step17_df$ENSEMBL[i]
  
  if (is.na(gene_id)) {
    next
  }
  
  # Tumor-only module
  module_i <- moduleColors[gene_id]
  
  if (length(module_i) == 1 && !is.na(module_i)) {
    
    step17_df$TumorOnly_Module[i] <- unname(module_i)
    
    # Corresponding MM column
    mm_col <- paste0("MM", module_i)
    
    if (mm_col %in% colnames(MM_all)) {
      
      step17_df$MM[i] <- MM_all[
        gene_id,
        mm_col
      ]
    }
  }
}


# ------------------------------------------------------------
# 3. Calculate module-specific MM rank
# ------------------------------------------------------------

step17_df$MM_Rank <- NA_integer_
step17_df$MM_Percentile <- NA_real_
step17_df$Module_Size <- NA_integer_


for (i in seq_len(nrow(step17_df))) {
  
  gene_id <- step17_df$ENSEMBL[i]
  module_i <- step17_df$TumorOnly_Module[i]
  
  if (
    is.na(gene_id) ||
    is.na(module_i) ||
    is.na(step17_df$MM[i])
  ) {
    next
  }
  
  mm_col <- paste0("MM", module_i)
  
  if (!(mm_col %in% colnames(MM_all))) {
    next
  }
  
  
  # All MM values in this gene's own module
  module_genes <- names(moduleColors)[
    moduleColors == module_i
  ]
  
  module_mm <- MM_all[
    module_genes,
    mm_col
  ]
  
  module_mm <- module_mm[
    !is.na(module_mm)
  ]
  
  
  # Module size
  n_module <- length(module_mm)
  
  step17_df$Module_Size[i] <- n_module
  
  
  # Rank genes by raw MM, descending
  # Higher MM = higher rank
  rank_i <- rank(
    -module_mm,
    ties.method = "min"
  )[
    match(gene_id, names(module_mm))
  ]
  
  step17_df$MM_Rank[i] <- as.integer(rank_i)
  
  
  # Percentile definition used in the locked core analysis:
  #
  # percentile = 100 * (N - rank + 1) / N
  #
  step17_df$MM_Percentile[i] <-
    100 * (
      n_module - rank_i + 1
    ) / n_module
}


# ------------------------------------------------------------
# 4. Print all 30 candidates
# ------------------------------------------------------------

cat("\n2. Complete 30-candidate MM positioning:\n")

print(
  step17_df[
    ,
    c(
      "Symbol",
      "TumorOnly_Module",
      "Module_Size",
      "MM",
      "MM_Rank",
      "MM_Percentile"
    )
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Validate against locked core results
# ------------------------------------------------------------

cat("\n3. Validation against locked core results:\n")

core_check <- merge(
  step17_df[
    step17_df$Symbol %in% c(
      "AOC3",
      "MS4A7",
      "IL7R"
    ),
    c(
      "Symbol",
      "TumorOnly_Module",
      "Module_Size",
      "MM",
      "MM_Rank",
      "MM_Percentile"
    )
  ],
  core_position_final[
    ,
    c(
      "Symbol",
      "Module",
      "Module_Size",
      "MM",
      "MM_rank",
      "MM_percentile"
    )
  ],
  by = "Symbol",
  suffixes = c("_new", "_locked")
)

print(
  core_check,
  row.names = FALSE
)


# ============================================================
# STEP 17A
# Inspect candidate_master / step17_df column structure
# ============================================================

cat("\n============================================================\n")
cat("STEP 17A: Inspect candidate table columns\n")
cat("============================================================\n")

cat("\n1. candidate_master columns:\n")
print(colnames(candidate_master))

cat("\n2. step17_df columns:\n")
print(colnames(step17_df))

cat("\n3. Dimensions:\n")
cat("candidate_master:", nrow(candidate_master), "x", ncol(candidate_master), "\n")
cat("step17_df:", nrow(step17_df), "x", ncol(step17_df), "\n")

cat("\n4. Current step17_df preview:\n")
print(
  step17_df,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("STEP 17A COMPLETE\n")
cat("============================================================\n")

colnames(step17_df)

# ============================================================
# Paper 1
# STEP 17B
# Fix module-specific MM rank / percentile calculation
# ============================================================

cat("\n============================================================\n")
cat("STEP 17B: Recalculate MM rank / percentile\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Reset rank-related columns
# ------------------------------------------------------------

step17_df$MM_Rank <- NA_integer_
step17_df$MM_Percentile <- NA_real_
step17_df$Module_Size <- NA_integer_


# ------------------------------------------------------------
# 2. Calculate rank within each gene's own tumor-only module
# ------------------------------------------------------------

for (i in seq_len(nrow(step17_df))) {
  
  gene_id <- step17_df$ENSEMBL[i]
  module_i <- step17_df$TumorOnly_Module[i]
  
  # Skip genes without valid Ensembl/module/MM
  if (
    is.na(gene_id) ||
    is.na(module_i) ||
    is.na(step17_df$MM[i])
  ) {
    next
  }
  
  
  # Corresponding MM column
  mm_col <- paste0("MM", module_i)
  
  if (!(mm_col %in% colnames(MM_all))) {
    next
  }
  
  
  # ----------------------------------------------------------
  # Genes belonging to this tumor-only module
  # ----------------------------------------------------------
  
  module_genes <- rownames(MM_all)[
    moduleColors == module_i
  ]
  
  
  # ----------------------------------------------------------
  # Extract MM values
  # ----------------------------------------------------------
  
  module_mm <- MM_all[
    module_genes,
    mm_col,
    drop = TRUE
  ]
  
  
  # Explicitly restore gene IDs as names
  names(module_mm) <- module_genes
  
  
  # Remove missing MM values
  keep <- !is.na(module_mm)
  
  module_mm <- module_mm[keep]
  
  
  # ----------------------------------------------------------
  # Module size
  # ----------------------------------------------------------
  
  n_module <- length(module_mm)
  
  step17_df$Module_Size[i] <- n_module
  
  
  # ----------------------------------------------------------
  # Module-specific MM rank
  #
  # Higher raw MM = better rank
  #
  # rank 1 = highest MM
  # ----------------------------------------------------------
  
  rank_all <- rank(
    -module_mm,
    ties.method = "min"
  )
  
  
  # Locate current gene by Ensembl ID
  rank_i <- rank_all[
    match(gene_id, names(rank_all))
  ]
  
  
  # ----------------------------------------------------------
  # Save rank
  # ----------------------------------------------------------
  
  step17_df$MM_Rank[i] <- as.integer(rank_i)
  
  
  # ----------------------------------------------------------
  # MM percentile
  #
  # percentile =
  # 100 * (N - rank + 1) / N
  # ----------------------------------------------------------
  
  step17_df$MM_Percentile[i] <-
    100 * (
      n_module - rank_i + 1
    ) / n_module
}


# ------------------------------------------------------------
# 3. Print complete result
# ------------------------------------------------------------

cat("\n1. Complete 30-candidate MM positioning:\n")

print(
  step17_df[
    ,
    c(
      "Candidate_Order",
      "Symbol",
      "TumorOnly_Module",
      "Module_Size",
      "MM",
      "MM_Rank",
      "MM_Percentile"
    )
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. Check number of successful rank calculations
# ------------------------------------------------------------

cat("\n2. Rank calculation completeness:\n")

cat(
  "Valid MM:",
  sum(!is.na(step17_df$MM)),
  "/",
  nrow(step17_df),
  "\n"
)

cat(
  "Valid MM rank:",
  sum(!is.na(step17_df$MM_Rank)),
  "/",
  nrow(step17_df),
  "\n"
)

cat(
  "Valid MM percentile:",
  sum(!is.na(step17_df$MM_Percentile)),
  "/",
  nrow(step17_df),
  "\n"
)


# ------------------------------------------------------------
# 5. Core-gene validation
# ------------------------------------------------------------

cat("\n3. Core-gene validation:\n")

core_step17 <- step17_df[
  step17_df$Symbol %in%
    c("AOC3", "MS4A7", "IL7R"),
  c(
    "Symbol",
    "TumorOnly_Module",
    "Module_Size",
    "MM",
    "MM_Rank",
    "MM_Percentile"
  )
]

print(
  core_step17,
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Hard-stop validation against locked results
# ------------------------------------------------------------

aoc3_new <- step17_df[
  step17_df$Symbol == "AOC3",
]

ms4a7_new <- step17_df[
  step17_df$Symbol == "MS4A7",
]

il7r_new <- step17_df[
  step17_df$Symbol == "IL7R",
]


# AOC3
stopifnot(
  aoc3_new$TumorOnly_Module == "magenta",
  aoc3_new$Module_Size == 534,
  aoc3_new$MM_Rank == 7
)

stopifnot(
  abs(
    aoc3_new$MM_Percentile -
      98.87430
  ) < 1e-5
)


# MS4A7
stopifnot(
  ms4a7_new$TumorOnly_Module == "blue",
  ms4a7_new$Module_Size == 1798,
  ms4a7_new$MM_Rank == 267
)

stopifnot(
  abs(
    ms4a7_new$MM_Percentile -
      85.19755
  ) < 1e-5
)


# IL7R
stopifnot(
  il7r_new$TumorOnly_Module == "blue",
  il7r_new$Module_Size == 1798,
  il7r_new$MM_Rank == 120
)

stopifnot(
  abs(
    il7r_new$MM_Percentile -
      93.37785
  ) < 1e-5
)


# ------------------------------------------------------------
# 7. Final QC
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "STEP 17B QC PASSED\n"
)

cat(
  "AOC3  : magenta / 534 / rank 7 / 98.87430%\n"
)

cat(
  "MS4A7 : blue    / 1798 / rank 267 / 85.19755%\n"
)

cat(
  "IL7R  : blue    / 1798 / rank 120 / 93.37785%\n"
)

cat(
  "============================================================\n"
)


# ============================================================
# Paper 1
# STEP 17C
# Correct MM percentile to reproduce locked core results
# ============================================================

cat("\n============================================================\n")
cat("STEP 17C: Correct MM percentile calculation\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 1. Recalculate percentile using the locked definition
#
# percentile =
# 100 * (N - rank) / (N - 1)
#
# This gives:
# AOC3  = 98.87430%
# IL7R  = 93.37785%
# MS4A7 = 85.19755%
# ------------------------------------------------------------

step17_df$MM_Percentile <- NA_real_


valid_rank <- !is.na(step17_df$MM_Rank) &
  !is.na(step17_df$Module_Size) &
  step17_df$Module_Size > 1


step17_df$MM_Percentile[valid_rank] <-
  100 * (
    step17_df$Module_Size[valid_rank] -
      step17_df$MM_Rank[valid_rank]
  ) /
  (
    step17_df$Module_Size[valid_rank] - 1
  )


# ------------------------------------------------------------
# 2. Print all 30 candidates
# ------------------------------------------------------------

cat("\n1. Complete 30-candidate MM positioning:\n")

print(
  step17_df[
    ,
    c(
      "Candidate_Order",
      "Symbol",
      "TumorOnly_Module",
      "Module_Size",
      "MM",
      "MM_Rank",
      "MM_Percentile"
    )
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 3. Check completeness
# ------------------------------------------------------------

cat("\n2. Completeness:\n")

cat(
  "Valid MM:",
  sum(!is.na(step17_df$MM)),
  "/",
  nrow(step17_df),
  "\n"
)

cat(
  "Valid MM rank:",
  sum(!is.na(step17_df$MM_Rank)),
  "/",
  nrow(step17_df),
  "\n"
)

cat(
  "Valid MM percentile:",
  sum(!is.na(step17_df$MM_Percentile)),
  "/",
  nrow(step17_df),
  "\n"
)


# ------------------------------------------------------------
# 4. Core-gene validation
# ------------------------------------------------------------

cat("\n3. Core-gene validation:\n")

core_step17 <- step17_df[
  step17_df$Symbol %in%
    c("AOC3", "MS4A7", "IL7R"),
  c(
    "Symbol",
    "TumorOnly_Module",
    "Module_Size",
    "MM",
    "MM_Rank",
    "MM_Percentile"
  )
]

print(
  core_step17,
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Extract core genes
# ------------------------------------------------------------

aoc3_new <- step17_df[
  step17_df$Symbol == "AOC3",
]

ms4a7_new <- step17_df[
  step17_df$Symbol == "MS4A7",
]

il7r_new <- step17_df[
  step17_df$Symbol == "IL7R",
]


# ------------------------------------------------------------
# 6. Hard-stop validation
# ------------------------------------------------------------

# AOC3
stopifnot(
  aoc3_new$TumorOnly_Module == "magenta",
  aoc3_new$Module_Size == 534,
  aoc3_new$MM_Rank == 7
)

stopifnot(
  abs(
    aoc3_new$MM_Percentile -
      98.87430
  ) < 1e-5
)


# MS4A7
stopifnot(
  ms4a7_new$TumorOnly_Module == "blue",
  ms4a7_new$Module_Size == 1798,
  ms4a7_new$MM_Rank == 267
)

stopifnot(
  abs(
    ms4a7_new$MM_Percentile -
      85.19755
  ) < 1e-5
)


# IL7R
stopifnot(
  il7r_new$TumorOnly_Module == "blue",
  il7r_new$Module_Size == 1798,
  il7r_new$MM_Rank == 120
)

stopifnot(
  abs(
    il7r_new$MM_Percentile -
      93.37785
  ) < 1e-5
)


# ------------------------------------------------------------
# 7. Final QC
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "STEP 17C QC PASSED\n"
)

cat(
  "MM rank and percentile exactly reproduce the locked core results.\n"
)

cat(
  "AOC3  : rank 7 / 534 / 98.87430%\n"
)

cat(
  "MS4A7 : rank 267 / 1798 / 85.19755%\n"
)

cat(
  "IL7R  : rank 120 / 1798 / 93.37785%\n"
)

cat(
  "============================================================\n"
)

### ============================================================
# Paper 1
# STEP 18A
# 30-candidate gene-level TME association audit
#
# PURPOSE:
#   Audit the already-computed gene-level GS/P/FDR results
#   for the 30 candidate genes.
#
# IMPORTANT:
#   - DO NOT recompute Spearman correlation
#   - DO NOT change GS
#   - DO NOT change P values
#   - DO NOT change FDR
#   - DO NOT perform any new statistical test
#   - Only subset, merge, classify, rank, and summarize
#
# Locked baseline:
#   STEP 14 gene-level TME association results
#
# ============================================================

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
    core_recovery[
      c("Symbol", "WGCNA_Module"),
      drop = FALSE
    ],
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
# Paper 1
# STEP 19A
# Candidate-selection logic consistency audit
# 30 -> 24 -> 6 -> 3
#
# PURPOSE:
#   Audit the historical candidate-selection chain without
#   recalculating WGCNA, GS, MM, P-values, FDR, or candidate
#   selection.
#
# IMPORTANT:
#   - candidate_master is NOT modified
#   - candidate_master_18D is NOT modified
#   - No WGCNA recalculation
#   - No statistical recalculation
#   - No new candidate filtering
# ============================================================

cat("\n============================================================\n")
cat("STEP 19A: Candidate-selection logic consistency audit\n")
cat("30 -> 24 -> 6 -> 3\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Required objects
# ------------------------------------------------------------

required_objects <- c(
  "candidate_30",
  "candidate_master",
  "candidate_master_18D",
  "core_genes",
  "core_step17",
  "candidate_module_audit"
)

cat("1. Checking required objects...\n\n")

object_status <- sapply(
  required_objects,
  exists,
  inherits = FALSE
)

print(
  data.frame(
    Object = required_objects,
    Exists = object_status,
    row.names = NULL
  )
)

if (!all(object_status)) {
  stop(
    "\nSTOP: One or more required objects are missing.\n",
    "No audit was performed."
  )
}


# ------------------------------------------------------------
# 2. Lock the 30-candidate starting set
# ------------------------------------------------------------

candidate30 <- as.character(candidate_30)

cat("\n------------------------------------------------------------\n")
cat("2. Starting candidate set\n")
cat("------------------------------------------------------------\n\n")

cat("Number of candidates:", length(candidate30), "\n")
cat("Unique symbols:", length(unique(candidate30)), "\n")
cat("Duplicated symbols:\n")

print(
  candidate30[
    duplicated(candidate30) |
      duplicated(candidate30, fromLast = TRUE)
  ]
)

cat("\nCandidate order:\n")
print(candidate30)


# ============================================================
# Paper 1
# STEP 19-0
# Workspace recovery audit
#
# PURPOSE:
#   Check whether the locked STEP 18A-18D objects are
#   available in the current R session or can be recovered
#   from saved R objects.
#
# IMPORTANT:
#   No WGCNA recalculation
#   No statistical recalculation
#   No candidate reselection
# ============================================================

cat("\n============================================================\n")
cat("STEP 19-0: Workspace recovery audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Current working directory
# ------------------------------------------------------------

cat("1. Current working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 2. Check current R session
# ------------------------------------------------------------

required_objects <- c(
  "candidate_30",
  "candidate_master",
  "candidate_master_18D",
  "core_genes",
  "core_step17",
  "candidate_module_audit"
)

current_status <- sapply(
  required_objects,
  exists,
  inherits = FALSE
)

cat("2. Current R-session object status:\n\n")

print(
  data.frame(
    Object = required_objects,
    Exists = current_status,
    row.names = NULL
  )
)


# ------------------------------------------------------------
# 3. List RData files in working directory
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("3. RData files in current directory\n")
cat("------------------------------------------------------------\n\n")

rdata_files <- list.files(
  path = ".",
  pattern = "\\.RData$|\\.rdata$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(rdata_files) == 0) {
  
  cat("No .RData files found in the current directory.\n")
  
} else {
  
  print(rdata_files)
}


# ------------------------------------------------------------
# 4. List RDS files in working directory
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("4. RDS files in current directory\n")
cat("------------------------------------------------------------\n\n")

rds_files <- list.files(
  path = ".",
  pattern = "\\.rds$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(rds_files) == 0) {
  
  cat("No .rds files found in the current directory.\n")
  
} else {
  
  print(rds_files)
}


# ------------------------------------------------------------
# 5. Search recursively for likely STEP 18 files
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("5. Search for likely Paper 1 / STEP 18 files\n")
cat("------------------------------------------------------------\n\n")

all_files <- list.files(
  path = ".",
  pattern = "\\.(RData|rdata|rds|R|r)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(all_files) == 0) {
  
  cat("No R-related files found.\n")
  
} else {
  
  likely_files <- all_files[
    grepl(
      paste(
        c(
          "Paper.?1",
          "STEP.?18",
          "18A",
          "18B",
          "18C",
          "18D",
          "candidate",
          "WGCNA",
          "core"
        ),
        collapse = "|"
      ),
      basename(all_files),
      ignore.case = TRUE
    )
  ]
  
  if (length(likely_files) == 0) {
    
    cat(
      "No likely STEP 18 / candidate / WGCNA files ",
      "were identified by filename.\n",
      sep = ""
    )
    
  } else {
    
    print(likely_files)
  }
}


# ------------------------------------------------------------
# 6. Check .RData objects without loading them
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("6. Inspect available .RData files\n")
cat("------------------------------------------------------------\n\n")

if (length(rdata_files) > 0) {
  
  for (f in rdata_files) {
    
    cat("\nFile:\n")
    cat(normalizePath(f), "\n")
    
    objs <- load(
      f,
      envir = new.env(parent = emptyenv())
    )
    
    cat("Objects contained in this file:\n")
    print(objs)
  }
}


# ------------------------------------------------------------
# 7. Final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 19-0 FINAL STATUS\n")
cat("============================================================\n\n")

if (all(current_status)) {
  
  cat(
    "STATUS: READY\n",
    "All required STEP 19 objects are already available.\n",
    "STEP 19A can be run directly.\n",
    sep = ""
  )
  
} else if (length(rdata_files) > 0 || length(rds_files) > 0) {
  
  cat(
    "STATUS: RECOVERY REQUIRED\n",
    "The current R session does not contain the required objects,\n",
    "but saved RData/RDS files are available for inspection.\n",
    "Do NOT recalculate WGCNA or candidate selection.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "STATUS: SAVED-OBJECT SEARCH REQUIRED\n",
    "The required objects are absent from the current session,\n",
    "and no directly visible RData/RDS files were found in the\n",
    "current working directory.\n",
    "Do NOT recalculate WGCNA or candidate selection yet.\n",
    sep = ""
  )
}

### ============================================================
# Paper 1
# STEP 19A
# ORIGINAL 30-CANDIDATE ORIGIN AUDIT
#
# Purpose:
#   Audit the historical provenance of the original
#   30 integrated candidate genes.
#
# Historical workflow to verify:
#
#   TCGA tumor + normal
#          ↓
#        DEGs
#          +
#   Original WGCNA clinically relevant modules
#          ↓
#      intersection
#          ↓
#   30 initial integrated candidates
#          ↓
#     STRING / PPI
#          ↓
#   Cytoscape / network characterization
#          ↓
#   biological prioritization
#          ↓
#   AOC3 / MS4A7 / IL7R
#
# IMPORTANT:
#   READ-ONLY AUDIT
#   NO WGCNA RECALCULATION
#   NO DEG RECALCULATION
#   NO PPI RECALCULATION
#   NO GLOBAL ENVIRONMENT MODIFICATION
# ============================================================


cat("\n============================================================\n")
cat("STEP 19A: ORIGINAL 30-CANDIDATE ORIGIN AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Current workspace inventory
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("1. CURRENT WORKSPACE INVENTORY\n")
cat("------------------------------------------------------------\n")

current_objects <- ls(envir = .GlobalEnv)

cat("Objects currently in .GlobalEnv:", length(current_objects), "\n\n")

cat(
  paste(
    current_objects,
    collapse = ", "
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 2. Define audit directory
# ------------------------------------------------------------

audit_dir <- getwd()

cat("Audit directory:\n")
cat(audit_dir, "\n\n")


# ------------------------------------------------------------
# 3. File inventory
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("3. FILE INVENTORY\n")
cat("------------------------------------------------------------\n")

all_files <- list.files(
  path = audit_dir,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

cat("Total files found:", length(all_files), "\n\n")


# ------------------------------------------------------------
# 4. Candidate-related file search
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("4. CANDIDATE-RELATED FILE SEARCH\n")
cat("------------------------------------------------------------\n")

candidate_pattern <- paste(
  c(
    "candidate",
    "30",
    "DEG",
    "WGCNA",
    "module",
    "STRING",
    "string",
    "PPI",
    "Cytoscape",
    "hub"
  ),
  collapse = "|"
)

candidate_files <- all_files[
  grepl(
    candidate_pattern,
    basename(all_files),
    ignore.case = TRUE
  )
]

cat(
  "Candidate/history-related files:",
  length(candidate_files),
  "\n\n"
)

if (length(candidate_files) > 0) {
  print(
    data.frame(
      File = candidate_files,
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )
}

cat("\n")


# ------------------------------------------------------------
# 5. Search R scripts for historical workflow keywords
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("5. HISTORICAL WORKFLOW SCRIPT SEARCH\n")
cat("------------------------------------------------------------\n")

r_files <- all_files[
  grepl(
    "\\.R$|\\.r$",
    all_files,
    ignore.case = TRUE
  )
]

cat("R scripts found:", length(r_files), "\n\n")


workflow_keywords <- c(
  "30",
  "candidate",
  "intersection",
  "intersect",
  "DEG",
  "differential",
  "blue",
  "green",
  "module",
  "WGCNA",
  "STRING",
  "PPI",
  "Cytoscape",
  "hub",
  "AOC3",
  "MS4A7",
  "IL7R"
)


# ------------------------------------------------------------
# 5A. Read-only script keyword audit
# ------------------------------------------------------------

script_audit_list <- list()

if (length(r_files) > 0) {
  
  for (f in r_files) {
    
    txt <- tryCatch(
      readLines(
        f,
        warn = FALSE,
        encoding = "UTF-8"
      ),
      error = function(e) {
        character(0)
      }
    )
    
    if (length(txt) == 0) next
    
    hit <- rep(FALSE, length(txt))
    
    for (kw in workflow_keywords) {
      
      hit <- hit |
        grepl(
          kw,
          txt,
          ignore.case = TRUE,
          fixed = TRUE
        )
    }
    
    idx <- which(hit)
    
    if (length(idx) > 0) {
      
      script_audit_list[[f]] <- data.frame(
        File = f,
        Line = idx,
        Code = txt[idx],
        stringsAsFactors = FALSE
      )
    }
  }
}


if (length(script_audit_list) > 0) {
  
  script_audit <- do.call(
    rbind,
    script_audit_list
  )
  
  rownames(script_audit) <- NULL
  
  cat(
    "Keyword-containing script lines:",
    nrow(script_audit),
    "\n\n"
  )
  
  print(
    head(
      script_audit,
      300
    ),
    row.names = FALSE
  )
  
} else {
  
  script_audit <- data.frame()
  
  cat(
    "No keyword-containing script lines found.\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 6. Specifically inspect the 30-candidate integration script
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("6. SPECIFIC 30-CANDIDATE SCRIPT AUDIT\n")
cat("------------------------------------------------------------\n")

specific_candidates_script <- r_files[
  grepl(
    "30Candidates|30_Candidates|30.*candidate|candidate.*30",
    basename(r_files),
    ignore.case = TRUE
  )
]

cat(
  "Potential 30-candidate scripts:",
  length(specific_candidates_script),
  "\n\n"
)

if (length(specific_candidates_script) > 0) {
  
  print(
    specific_candidates_script
  )
  
  cat("\n")
  
  for (f in specific_candidates_script) {
    
    cat("\n========================================================\n")
    cat("SCRIPT:\n")
    cat(f, "\n")
    cat("========================================================\n")
    
    txt <- tryCatch(
      readLines(
        f,
        warn = FALSE,
        encoding = "UTF-8"
      ),
      error = function(e) {
        character(0)
      }
    )
    
    if (length(txt) == 0) {
      cat("Unable to read script.\n")
      next
    }
    
    cat(
      paste(
        seq_along(txt),
        txt,
        sep = ": ",
        collapse = "\n"
      )
    )
    
    cat("\n")
  }
  
} else {
  
  cat(
    "No specifically named 30-candidate script found.\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 7. Locate historical WGCNA / candidate RData files
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("7. HISTORICAL RDATA FILE AUDIT\n")
cat("------------------------------------------------------------\n")

rdata_files <- all_files[
  grepl(
    "\\.RData$|\\.rdata$|\\.rda$",
    all_files,
    ignore.case = TRUE
  )
]

cat(
  "RData/RDA files:",
  length(rdata_files),
  "\n\n"
)

if (length(rdata_files) > 0) {
  
  print(
    data.frame(
      File = rdata_files,
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )
}

cat("\n")


# ------------------------------------------------------------
# 8. Inspect RData contents WITHOUT loading into GlobalEnv
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("8. RDATA OBJECT INVENTORY — ISOLATED ENVIRONMENTS\n")
cat("------------------------------------------------------------\n")

rdata_inventory <- list()

if (length(rdata_files) > 0) {
  
  for (f in rdata_files) {
    
    e <- new.env(
      parent = emptyenv()
    )
    
    loaded <- tryCatch(
      load(
        f,
        envir = e
      ),
      error = function(err) {
        NULL
      }
    )
    
    if (is.null(loaded)) next
    
    objs <- ls(
      envir = e,
      all.names = TRUE
    )
    
    if (length(objs) == 0) next
    
    rdata_inventory[[f]] <- data.frame(
      File = f,
      Object = objs,
      stringsAsFactors = FALSE
    )
  }
}


if (length(rdata_inventory) > 0) {
  
  rdata_inventory_df <- do.call(
    rbind,
    rdata_inventory
  )
  
  rownames(rdata_inventory_df) <- NULL
  
  cat(
    "Total RData objects identified:",
    nrow(rdata_inventory_df),
    "\n\n"
  )
  
  print(
    rdata_inventory_df,
    row.names = FALSE
  )
  
} else {
  
  rdata_inventory_df <- data.frame()
  
  cat(
    "No RData objects could be inventoried.\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 9. Search RData object names for historical candidates
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("9. RDATA OBJECT NAME SEARCH\n")
cat("------------------------------------------------------------\n")

if (nrow(rdata_inventory_df) > 0) {
  
  relevant_objects <- rdata_inventory_df[
    grepl(
      paste(
        c(
          "candidate",
          "DEG",
          "gene",
          "module",
          "color",
          "blue",
          "green",
          "WGCNA",
          "STRING",
          "PPI",
          "hub",
          "intersection",
          "intersect",
          "AOC3",
          "MS4A7",
          "IL7R"
        ),
        collapse = "|"
      ),
      rdata_inventory_df$Object,
      ignore.case = TRUE
    ),
  ]
  
  cat(
    "Relevant historical objects:",
    nrow(relevant_objects),
    "\n\n"
  )
  
  if (nrow(relevant_objects) > 0) {
    print(
      relevant_objects,
      row.names = FALSE
    )
  }
  
} else {
  
  relevant_objects <- data.frame()
}

cat("\n")


# ------------------------------------------------------------
# 10. Inspect known historical checkpoint files
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("10. KNOWN HISTORICAL CHECKPOINTS\n")
cat("------------------------------------------------------------\n")

known_patterns <- c(
  "Workflow-04-2D_final_WGCNA_candidate_checkpoint.RData",
  "Workflow-04-2D_current_workspace.RData",
  "Paper1_WORKFLOW03_14F_4D_2E_FINAL.RData",
  "Paper1_ME_brown_checkpoint_Step10G.RData",
  "Step13G_FINAL_checkpoint.RData",
  "Workflow-03F-9_provenance_audit.RData"
)


known_files <- rdata_files[
  basename(rdata_files) %in% known_patterns
]

cat(
  "Known checkpoint files located:",
  length(known_files),
  "\n\n"
)

if (length(known_files) > 0) {
  
  for (f in known_files) {
    
    cat("\n========================================================\n")
    cat("CHECKPOINT:\n")
    cat(f, "\n")
    cat("========================================================\n")
    
    e <- new.env(
      parent = emptyenv()
    )
    
    loaded <- tryCatch(
      load(
        f,
        envir = e
      ),
      error = function(err) {
        NULL
      }
    )
    
    if (!is.null(loaded)) {
      
      objs <- ls(
        envir = e,
        all.names = TRUE
      )
      
      cat(
        "Objects:",
        length(objs),
        "\n\n"
      )
      
      print(objs)
      
    } else {
      
      cat("Unable to load checkpoint in isolated environment.\n")
    }
  }
}

cat("\n")


# ------------------------------------------------------------
# 11. RDS file inventory
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("11. RDS FILE INVENTORY\n")
cat("------------------------------------------------------------\n")

rds_files <- all_files[
  grepl(
    "\\.rds$",
    all_files,
    ignore.case = TRUE
  )
]

cat(
  "RDS files:",
  length(rds_files),
  "\n\n"
)

if (length(rds_files) > 0) {
  
  print(
    data.frame(
      File = rds_files,
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )
}

cat("\n")


# ------------------------------------------------------------
# 12. Inspect candidate-related RDS objects
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("12. CANDIDATE-RELATED RDS AUDIT\n")
cat("------------------------------------------------------------\n")

candidate_rds <- rds_files[
  grepl(
    "candidate|30|DEG|WGCNA|core|module",
    basename(rds_files),
    ignore.case = TRUE
  )
]

cat(
  "Candidate-related RDS files:",
  length(candidate_rds),
  "\n\n"
)

if (length(candidate_rds) > 0) {
  
  for (f in candidate_rds) {
    
    cat("\n========================================================\n")
    cat("RDS:\n")
    cat(f, "\n")
    cat("========================================================\n")
    
    obj <- tryCatch(
      readRDS(f),
      error = function(e) {
        NULL
      }
    )
    
    if (is.null(obj)) {
      cat("Unable to read RDS.\n")
      next
    }
    
    cat(
      "Class:",
      paste(class(obj), collapse = ", "),
      "\n"
    )
    
    cat(
      "Length:",
      length(obj),
      "\n"
    )
    
    if (is.data.frame(obj)) {
      
      cat(
        "Dimensions:",
        paste(dim(obj), collapse = " × "),
        "\n"
      )
      
      cat(
        "Columns:\n"
      )
      
      print(
        names(obj)
      )
      
    } else if (is.matrix(obj)) {
      
      cat(
        "Dimensions:",
        paste(dim(obj), collapse = " × "),
        "\n"
      )
      
    } else {
      
      cat(
        "Object structure:\n"
      )
      
      print(
        str(
          obj,
          max.level = 1
        )
      )
    }
  }
}

cat("\n")


# ------------------------------------------------------------
# 13. Search all TSV/CSV files for STRING/PPI/Cytoscape evidence
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("13. STRING / PPI / CYTOSCAPE FILE SEARCH\n")
cat("------------------------------------------------------------\n")

network_files <- all_files[
  grepl(
    "string|ppi|cyto|cytoscape|interaction|degree|network|hub",
    basename(all_files),
    ignore.case = TRUE
  )
]

cat(
  "Potential network-analysis files:",
  length(network_files),
  "\n\n"
)

if (length(network_files) > 0) {
  
  print(
    data.frame(
      File = network_files,
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )
}

cat("\n")


# ------------------------------------------------------------
# 14. Specifically inspect STRING node-degree / interaction files
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("14. STRING NODE / INTERACTION FILE AUDIT\n")
cat("------------------------------------------------------------\n")

string_files <- network_files[
  grepl(
    "string_node|string_interaction|interaction_short|degree",
    basename(network_files),
    ignore.case = TRUE
  )
]

cat(
  "Potential STRING files:",
  length(string_files),
  "\n\n"
)

if (length(string_files) > 0) {
  
  for (f in string_files) {
    
    cat("\n========================================================\n")
    cat("NETWORK FILE:\n")
    cat(f, "\n")
    cat("========================================================\n")
    
    dat <- tryCatch(
      read.delim(
        f,
        header = TRUE,
        sep = "\t",
        check.names = FALSE,
        stringsAsFactors = FALSE,
        quote = ""
      ),
      error = function(e) {
        NULL
      }
    )
    
    if (is.null(dat)) {
      
      cat(
        "Could not read as tab-delimited file.\n"
      )
      
      next
    }
    
    cat(
      "Dimensions:",
      paste(dim(dat), collapse = " × "),
      "\n\n"
    )
    
    cat(
      "Columns:\n"
    )
    
    print(
      names(dat)
    )
    
    cat("\nFirst 10 rows:\n")
    
    print(
      head(
        dat,
        10
      ),
      row.names = FALSE
    )
  }
}

cat("\n")


# ------------------------------------------------------------
# 15. Search all text-like files for exact core genes
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("15. CORE-GENE HISTORICAL TRACE\n")
cat("------------------------------------------------------------\n")

core_gene_pattern <- "AOC3|MS4A7|IL7R"

text_files <- all_files[
  grepl(
    "\\.txt$|\\.tsv$|\\.csv$|\\.R$|\\.r$|\\.md$|\\.log$",
    all_files,
    ignore.case = TRUE
  )
]

core_gene_hits <- list()

if (length(text_files) > 0) {
  
  for (f in text_files) {
    
    txt <- tryCatch(
      readLines(
        f,
        warn = FALSE,
        encoding = "UTF-8"
      ),
      error = function(e) {
        character(0)
      }
    )
    
    if (length(txt) == 0) next
    
    idx <- grep(
      core_gene_pattern,
      txt,
      ignore.case = FALSE
    )
    
    if (length(idx) > 0) {
      
      core_gene_hits[[f]] <- data.frame(
        File = f,
        Line = idx,
        Text = txt[idx],
        stringsAsFactors = FALSE
      )
    }
  }
}


if (length(core_gene_hits) > 0) {
  
  core_gene_hits_df <- do.call(
    rbind,
    core_gene_hits
  )
  
  rownames(core_gene_hits_df) <- NULL
  
  cat(
    "Core-gene-containing text lines:",
    nrow(core_gene_hits_df),
    "\n\n"
  )
  
  print(
    head(
      core_gene_hits_df,
      300
    ),
    row.names = FALSE
  )
  
} else {
  
  core_gene_hits_df <- data.frame()
  
  cat(
    "No core-gene text hits found.\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 16. Search text files for the original 30-candidate workflow
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("16. ORIGINAL 30-CANDIDATE WORKFLOW TRACE\n")
cat("------------------------------------------------------------\n")

historical_terms <- c(
  "30 candidate",
  "30 candidates",
  "candidate genes",
  "DEG",
  "intersect",
  "intersection",
  "blue",
  "green",
  "module",
  "STRING",
  "Cytoscape",
  "PPI"
)

historical_hits <- list()

if (length(text_files) > 0) {
  
  for (f in text_files) {
    
    txt <- tryCatch(
      readLines(
        f,
        warn = FALSE,
        encoding = "UTF-8"
      ),
      error = function(e) {
        character(0)
      }
    )
    
    if (length(txt) == 0) next
    
    hit <- rep(FALSE, length(txt))
    
    for (term in historical_terms) {
      
      hit <- hit |
        grepl(
          term,
          txt,
          ignore.case = TRUE,
          fixed = TRUE
        )
    }
    
    idx <- which(hit)
    
    if (length(idx) > 0) {
      
      historical_hits[[f]] <- data.frame(
        File = f,
        Line = idx,
        Text = txt[idx],
        stringsAsFactors = FALSE
      )
    }
  }
}


if (length(historical_hits) > 0) {
  
  historical_hits_df <- do.call(
    rbind,
    historical_hits
  )
  
  rownames(historical_hits_df) <- NULL
  
  cat(
    "Historical workflow hits:",
    nrow(historical_hits_df),
    "\n\n"
  )
  
  print(
    head(
      historical_hits_df,
      500
    ),
    row.names = FALSE
  )
  
} else {
  
  historical_hits_df <- data.frame()
  
  cat(
    "No historical workflow text hits found.\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 17. Explicit audit questions
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP 19A — AUDIT QUESTIONS\n")
cat("============================================================\n\n")

cat(
  "Q1. Can the original 30-gene candidate set be located?\n"
)

cat(
  "Q2. Can the original tumor+normal WGCNA modules be identified?\n"
)

cat(
  "Q3. Can the historical blue + green module assignments be verified?\n"
)

cat(
  "Q4. Can the DEG set be identified?\n"
)

cat(
  "Q5. Can DEG ∩ module genes be independently documented as the source of 30 candidates?\n"
)

cat(
  "Q6. Can the 30 genes be traced into STRING/PPI analysis?\n"
)

cat(
  "Q7. Can Cytoscape/network characterization be traced after STRING?\n"
)

cat(
  "Q8. Can AOC3/MS4A7/IL7R be traced as the final prioritized study genes?\n\n"
)


# ------------------------------------------------------------
# 18. IMPORTANT STATUS MESSAGE
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP 19A STATUS\n")
cat("============================================================\n\n")

cat(
  "This step is an AUDIT ONLY.\n"
)

cat(
  "No WGCNA, DEG, MM, GS, PPI, or candidate-selection\n"
)

cat(
  "recalculation has been performed.\n\n"
)

cat(
  "Do NOT interpret current tumor-only blue/magenta modules\n"
)

cat(
  "as the historical original blue/green modules.\n\n"
)

cat(
  "Historical provenance must be established from archived\n"
)

cat(
  "scripts, checkpoints, candidate files, and network files.\n\n"
)

cat(
  "STEP 19A COMPLETE — WAIT FOR OUTPUT REVIEW.\n"
)

cat("============================================================\n")


###STEP 18F — Pre-STEP19 Full Workspace Checkpoint
# ============================================================
# Paper 1 Revision
# STEP 18F
# Pre-STEP19 Full Workspace Checkpoint
# Date: 2026-09-14
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

list.files(
  getwd(),
  pattern = "\\.R$",
  full.names = FALSE
)


# ============================================================
# Paper 1 Revision
# STEP 18G
# Archive all existing R scripts before STEP 19
# Date: 
#
# Purpose:
#   Preserve the exact R scripts currently present in the
#   working directory before starting STEP 19 / SF4 plotting.
#
# IMPORTANT:
#   This step ONLY copies existing script files.
#   It does NOT execute or modify any analysis.
# ============================================================

cat("\n============================================================\n")
cat("STEP 18G: PRE-STEP19 SCRIPT ARCHIVE\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Checkpoint directory
# ------------------------------------------------------------

checkpoint_root <- file.path(
  getwd(),
  "Paper1_STEP19_PREP_20260914"
)

script_archive_dir <- file.path(
  checkpoint_root,
  "scripts"
)

dir.create(
  script_archive_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Exact list of existing R scripts
# ------------------------------------------------------------

script_files <- c(
  "Paper1_ME_brown_checkpoint_Step10G.R",
  "Paper1_ME_brown_Gender_checkpoint_Final.R",
  "paper1_Step13G_FINAL_checkpoint.R",
  "Paper1_tumor-only WGCNA.R",
  "Paper1_WORKFLOW04_30Candidates_TumorOnly_Integration.R",
  "WORKFLOW-03_Step14F-1_CoreGene_ModuleLocalization.R",
  "Workflow-03F-9_provenance_audit.R",
  "Workflow-04-1F-1_clinical_cleaning_checkpoint.R",
  "Workflow-04-2D-5 Soft-threshold checkpoint.R"
)


# ------------------------------------------------------------
# 3. Verify that every expected script exists
# ------------------------------------------------------------

script_status <- data.frame(
  Script = script_files,
  Exists = file.exists(
    file.path(
      getwd(),
      script_files
    )
  ),
  stringsAsFactors = FALSE
)

cat("Script verification:\n")
print(script_status)

cat("\n")


# ------------------------------------------------------------
# 4. Stop if any expected script is missing
# ------------------------------------------------------------

if (!all(script_status$Exists)) {
  
  missing_scripts <- script_status$Script[
    !script_status$Exists
  ]
  
  stop(
    paste(
      "ARCHIVE STOPPED.\n",
      "The following expected scripts were not found:\n",
      paste(missing_scripts, collapse = "\n")
    )
  )
}


# ------------------------------------------------------------
# 5. Copy scripts WITHOUT modification
# ------------------------------------------------------------

copy_result <- file.copy(
  from = file.path(
    getwd(),
    script_files
  ),
  to = file.path(
    script_archive_dir,
    script_files
  ),
  overwrite = TRUE
)

archive_status <- data.frame(
  Script = script_files,
  Copied = copy_result,
  stringsAsFactors = FALSE
)

cat("Archive result:\n")
print(archive_status)

cat("\n")


# ------------------------------------------------------------
# 6. Save archive manifest
# ------------------------------------------------------------

script_manifest <- data.frame(
  Script = script_files,
  Original_Path = file.path(
    getwd(),
    script_files
  ),
  Archived_Path = file.path(
    script_archive_dir,
    script_files
  ),
  File_Size_KB = round(
    file.info(
      file.path(
        getwd(),
        script_files
      )
    )$size / 1024,
    2
  ),
  stringsAsFactors = FALSE
)

write.csv(
  script_manifest,
  file = file.path(
    checkpoint_root,
    "audit",
    "STEP18G_script_archive_manifest_20260914.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Final verification
# ------------------------------------------------------------

final_check <- file.exists(
  file.path(
    script_archive_dir,
    script_files
  )
)

cat("Final archive verification:\n")
print(
  data.frame(
    Script = script_files,
    Archived = final_check,
    stringsAsFactors = FALSE
  )
)

cat("\n")


if (all(final_check)) {
  
  cat("============================================================\n")
  cat("STEP 18G STATUS: PASS\n")
  cat("All 9 R scripts were archived successfully.\n")
  cat("============================================================\n")
  
} else {
  
  stop(
    "STEP 18G FAILED: one or more scripts were not archived."
  )
}


list.files(
  file.path(
    getwd(),
    "Paper1_STEP19_PREP_20260914"
  ),
  recursive = TRUE,
  full.names = FALSE
)

list.files(
  getwd(),
  pattern = "\\.(RData|rds)$",
  full.names = FALSE
)



# ============================================================
# Paper 1 Revision
# STEP 18H
# Archive all existing RData / RDS files
# Date: 2026-09-14
#
# Purpose:
#   Preserve all existing analysis data/checkpoint files before
#   starting STEP 19 and Supplementary Figure 4 plotting.
#
# IMPORTANT:
#   This step ONLY copies files.
#   No analysis is recalculated.
#   No R objects are modified.
# ============================================================

cat("\n============================================================\n")
cat("STEP 18H: PRE-STEP19 RDATA / RDS ARCHIVE\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define checkpoint directories
# ------------------------------------------------------------

checkpoint_root <- file.path(
  getwd(),
  "Paper1_STEP19_PREP_20260914"
)

locked_data_dir <- file.path(
  checkpoint_root,
  "locked_data"
)

audit_dir <- file.path(
  checkpoint_root,
  "audit"
)

dir.create(
  locked_data_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Identify ALL existing RData / RDS files
# ------------------------------------------------------------

data_files <- list.files(
  getwd(),
  pattern = "\\.(RData|rds)$",
  full.names = FALSE
)

cat("Number of RData/RDS files detected:\n")
cat(length(data_files), "\n\n")

cat("Files to archive:\n")
print(data_files)

cat("\n")


# ------------------------------------------------------------
# 3. Safety check
# ------------------------------------------------------------

if (length(data_files) == 0) {
  
  stop(
    "STEP 18H STOPPED: No RData/RDS files were detected."
  )
  
}


# ------------------------------------------------------------
# 4. Copy all files WITHOUT modification
# ------------------------------------------------------------

copy_result <- file.copy(
  from = file.path(
    getwd(),
    data_files
  ),
  to = file.path(
    locked_data_dir,
    data_files
  ),
  overwrite = TRUE
)


# ------------------------------------------------------------
# 5. Build archive status table
# ------------------------------------------------------------

archive_status <- data.frame(
  File = data_files,
  Copied = copy_result,
  Original_Size_MB = round(
    file.info(
      file.path(
        getwd(),
        data_files
      )
    )$size / 1024^2,
    3
  ),
  Archived_Size_MB = round(
    file.info(
      file.path(
        locked_data_dir,
        data_files
      )
    )$size / 1024^2,
    3
  ),
  stringsAsFactors = FALSE
)


cat("Archive status:\n")
print(archive_status)

cat("\n")


# ------------------------------------------------------------
# 6. Verify copied files
# ------------------------------------------------------------

archive_status$Exists_After_Copy <- file.exists(
  file.path(
    locked_data_dir,
    data_files
  )
)


# ------------------------------------------------------------
# 7. Verify file sizes
# ------------------------------------------------------------

archive_status$Size_Match <-
  archive_status$Original_Size_MB ==
  archive_status$Archived_Size_MB


# ------------------------------------------------------------
# 8. Save manifest
# ------------------------------------------------------------

manifest_file <- file.path(
  audit_dir,
  "STEP18H_RData_RDS_archive_manifest_20260914.csv"
)

write.csv(
  archive_status,
  file = manifest_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Final PASS / FAIL
# ------------------------------------------------------------

all_exist <- all(
  archive_status$Exists_After_Copy
)

all_size_match <- all(
  archive_status$Size_Match
)


cat("\n============================================================\n")
cat("STEP 18H FINAL VERIFICATION\n")
cat("============================================================\n\n")

cat(
  "All files exist after copy:",
  all_exist,
  "\n"
)

cat(
  "All file sizes match:",
  all_size_match,
  "\n"
)

cat(
  "Number of files archived:",
  nrow(archive_status),
  "\n\n"
)


if (all_exist && all_size_match) {
  
  cat("============================================================\n")
  cat("STEP 18H STATUS: PASS\n")
  cat("All RData/RDS files were archived successfully.\n")
  cat("============================================================\n")
  
} else {
  
  stop(
    "STEP 18H FAILED: archive verification did not pass."
  )
  
}


