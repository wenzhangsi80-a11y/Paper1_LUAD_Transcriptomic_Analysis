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

load("input_dir/WORKFLOW-03_Step14A-14E_FINAL_LOCKED.RData")


ls()

#then check：

dim(datExpr_MAD)

dim(MEs)

dim(expr_estimate)

dim(estimate_scores_wgcna)

dim(ME_color_formal_raw)

dim(module_tme_results)

#Re-examine the three most critical results：

# 1. 150 tests
nrow(module_tme_results)

# 2. FDR significant
sum(module_tme_results$FDR < 0.05)

# 3. Brown = ME3
cor(
  ME_color_formal_raw[, "MEbrown"],
  MEs[, "ME3"],
  method = "pearson"
)

#Final check of Brown–TME：

subset(
  module_tme_results,
  Module == "ME3"
)

############################################################
# WORKFLOW-03
# Step 14F-1: Core-gene localization in tumor-only WGCNA
#
# Core genes:
#   AOC3
#   MS4A7
#   IL7R
#
# Purpose:
#   1. Map core gene symbols to Ensembl IDs
#   2. Confirm their presence in datExpr_MAD
#   3. Determine WGCNA module membership
#   4. Map module color -> internal ME number
#   5. Perform basic expression / module QC
#
# IMPORTANT:
#   This step only performs module localization.
#   MM / GS / TME integration will be performed AFTER
#   the localization results are confirmed.
############################################################


############################################################
# 0. BASIC SETTINGS
############################################################

options(stringsAsFactors = FALSE)

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


############################################################
# 1. CHECK REQUIRED OBJECTS
############################################################

required_objects <- c(
  "datExpr_MAD",
  "moduleColors",
  "MEs",
  "annotation"
)

missing_objects <- required_objects[
  !required_objects %in% ls(envir = .GlobalEnv)
]

if (length(missing_objects) > 0) {
  
  stop(
    paste0(
      "The following required objects are missing:\n",
      paste(missing_objects, collapse = ", "),
      "\n\nPlease load the Step 14A-14E locked workspace first."
    )
  )
}

cat("\n============================================\n")
cat("STEP 14F-1: CORE-GENE MODULE LOCALIZATION\n")
cat("============================================\n\n")


############################################################
# 2. BASIC WGCNA DATA QC
############################################################

cat("2. WGCNA data dimensions:\n")

print(dim(datExpr_MAD))

cat("\n")

cat("Number of module assignments:\n")
print(length(moduleColors))

cat("\n")

cat("Number of module eigengenes:\n")
print(ncol(MEs))

cat("\n")

if (length(moduleColors) != ncol(datExpr_MAD)) {
  
  stop(
    "moduleColors length does not equal the number of genes in datExpr_MAD."
  )
}


############################################################
# 3. INSPECT ANNOTATION OBJECT
############################################################

cat("\n============================================\n")
cat("3. ANNOTATION QC\n")
cat("============================================\n\n")

cat("Annotation dimensions:\n")
print(dim(annotation))

cat("\nAnnotation columns:\n")
print(colnames(annotation))


############################################################
# 4. STANDARDIZE ANNOTATION COLUMN NAMES
############################################################

# Expected annotation:
#
#   Ensembl_ID
#   gene_name
#   gene_type
#   seqname
#   start
#   end
#   strand
#   Annotation_status
#
# The following checks prevent accidental use of
# the wrong columns.

required_annotation_columns <- c(
  "Ensembl_ID",
  "gene_name"
)

missing_annotation_columns <- required_annotation_columns[
  !required_annotation_columns %in% colnames(annotation)
]

if (length(missing_annotation_columns) > 0) {
  
  stop(
    paste0(
      "Missing annotation columns: ",
      paste(
        missing_annotation_columns,
        collapse = ", "
      )
    )
  )
  
  
}


############################################################
# 5. CREATE GENE -> ENSEMBL LOOKUP FOR CORE GENES
############################################################

core_annotation <- annotation[
  annotation$gene_name %in% core_genes,
  c(
    "Ensembl_ID",
    "gene_name"
  ),
  drop = FALSE
]

core_annotation <- core_annotation[
  order(
    match(
      core_annotation$gene_name,
      core_genes
    )
  ),
  ,
  drop = FALSE
]

cat("\n============================================\n")
cat("5. CORE GENE ANNOTATION\n")
cat("============================================\n\n")

print(core_annotation)


############################################################
# 6. CHECK FOR MULTIPLE ENSEMBL IDs PER CORE GENE
############################################################

core_id_count <- aggregate(
  Ensembl_ID ~ gene_name,
  data = core_annotation,
  FUN = length
)

colnames(core_id_count)[2] <- "n_Ensembl_ID"

cat("\nNumber of Ensembl IDs per core gene:\n")
print(core_id_count)


############################################################
# 7. IDENTIFY CORE GENE IDs PRESENT IN WGCNA
############################################################

wgcna_ensembl_ids <- colnames(datExpr_MAD)

core_annotation$In_WGCNA <- (
  core_annotation$Ensembl_ID %in%
    wgcna_ensembl_ids
)

cat("\n============================================\n")
cat("7. CORE GENE PRESENCE IN WGCNA\n")
cat("============================================\n\n")

print(core_annotation)


############################################################
# 8. STOP IF ANY CORE GENE IS NOT FOUND
############################################################

missing_core_genes <- core_genes[
  !core_genes %in%
    core_annotation$gene_name[
      core_annotation$In_WGCNA
    ]
]

if (length(missing_core_genes) > 0) {
  
  warning(
    paste0(
      "The following core genes were not found in ",
      "the tumor-only WGCNA matrix: ",
      paste(
        missing_core_genes,
        collapse = ", "
      )
    )
  )
}


############################################################
# 9. KEEP ONLY CORE IDs PRESENT IN WGCNA
############################################################

core_in_wgcna <- core_annotation[
  core_annotation$In_WGCNA,
  ,
  drop = FALSE
]

if (nrow(core_in_wgcna) == 0) {
  
  stop(
    "None of the three core genes was found in datExpr_MAD."
  )
}


############################################################
# 10. ADD WGCNA MODULE COLOR
############################################################

# moduleColors corresponds exactly to the columns
# of datExpr_MAD.

core_in_wgcna$ModuleColor <- moduleColors[
  match(
    core_in_wgcna$Ensembl_ID,
    wgcna_ensembl_ids
  )
]


############################################################
# 11. BASIC EXPRESSION SUMMARY
############################################################

core_expression_summary <- data.frame(
  
  Gene = core_in_wgcna$gene_name,
  
  Ensembl_ID = core_in_wgcna$Ensembl_ID,
  
  ModuleColor = core_in_wgcna$ModuleColor,
  
  MeanExpression = NA_real_,
  
  SDExpression = NA_real_,
  
  MedianExpression = NA_real_,
  
  MinExpression = NA_real_,
  
  MaxExpression = NA_real_,
  
  stringsAsFactors = FALSE
  
)

for (i in seq_len(nrow(core_expression_summary))) {
  
  gene_id <- core_expression_summary$Ensembl_ID[i]
  
  x <- datExpr_MAD[
    ,
    gene_id
  ]
  
  core_expression_summary$MeanExpression[i] <-
    mean(x, na.rm = TRUE)
  
  core_expression_summary$SDExpression[i] <-
    sd(x, na.rm = TRUE)
  
  core_expression_summary$MedianExpression[i] <-
    median(x, na.rm = TRUE)
  
  core_expression_summary$MinExpression[i] <-
    min(x, na.rm = TRUE)
  
  core_expression_summary$MaxExpression[i] <-
    max(x, na.rm = TRUE)
}


############################################################
# 12. MAP MODULE COLOR -> INTERNAL ME NUMBER
############################################################

# The formal mapping was validated in Step 14E.
#
# We reconstruct it directly from the formal
# moduleEigengenes object rather than using the
# previously discarded rowMeans(scale()) approach.

ME_color_names <- colnames(ME_color_formal_raw)

ME_internal_names <- colnames(MEs)

module_mapping <- data.frame(
  ModuleColor = sub(
    "^ME",
    "",
    ME_color_names
  ),
  Internal_ME = NA_character_,
  Pearson_r = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(module_mapping))) {
  
  color_ME <- ME_color_names[i]
  
  correlations <- sapply(
    ME_internal_names,
    function(me_name) {
      
      cor(
        ME_color_formal_raw[, color_ME],
        MEs[, me_name],
        method = "pearson",
        use = "pairwise.complete.obs"
      )
      
    }
  )
  
  best_index <- which.max(
    abs(correlations)
  )
  
  module_mapping$Internal_ME[i] <-
    ME_internal_names[best_index]
  
  module_mapping$Pearson_r[i] <-
    correlations[best_index]
}


############################################################
# 13. CHECK FORMAL MAPPING
############################################################

cat("\n============================================\n")
cat("13. FORMAL MODULE MAPPING QC\n")
cat("============================================\n\n")

print(module_mapping)

cat("\nNumber of mappings with |r| < 0.999999:\n")

print(
  sum(
    abs(module_mapping$Pearson_r) <
      0.999999
  )
)

cat("\nNumber of unique internal ME assignments:\n")

print(
  length(
    unique(
      module_mapping$Internal_ME
    )
  )
)


############################################################
# 14. ADD INTERNAL ME NUMBER TO CORE GENES
############################################################

core_expression_summary$Internal_ME <-
  module_mapping$Internal_ME[
    match(
      core_expression_summary$ModuleColor,
      module_mapping$ModuleColor
    )
  ]

core_expression_summary$ModuleMapping_r <-
  module_mapping$Pearson_r[
    match(
      core_expression_summary$ModuleColor,
      module_mapping$ModuleColor
    )
  ]


############################################################
# 15. FINAL CORE-GENE MODULE LOCALIZATION TABLE
############################################################

core_module_localization <- core_expression_summary[
  ,
  c(
    "Gene",
    "Ensembl_ID",
    "ModuleColor",
    "Internal_ME",
    "ModuleMapping_r",
    "MeanExpression",
    "SDExpression",
    "MedianExpression",
    "MinExpression",
    "MaxExpression"
  )
]

core_module_localization <- core_module_localization[
  order(
    match(
      core_module_localization$Gene,
      core_genes
    )
  ),
  ,
  drop = FALSE
]


############################################################
# 16. PRINT FINAL RESULT
############################################################

cat("\n============================================\n")
cat("16. FINAL CORE-GENE MODULE LOCALIZATION\n")
cat("============================================\n\n")

print(
  core_module_localization,
  row.names = FALSE
)


############################################################
# 17. CHECK WHETHER CORE GENES SHARE A MODULE
############################################################

core_module_table <- table(
  core_module_localization$ModuleColor
)

cat("\n============================================\n")
cat("17. CORE-GENE MODULE DISTRIBUTION\n")
cat("============================================\n\n")

print(core_module_table)

shared_modules <- names(
  core_module_table[
    core_module_table > 1
  ]
)

if (length(shared_modules) == 0) {
  
  cat(
    "\nNo WGCNA module contains more than one of ",
    "the three core genes.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "\nModules containing multiple core genes:\n"
  )
  
  print(shared_modules)
}


############################################################
# 18. CHECK WHETHER ALL THREE ARE IN SAME MODULE
############################################################

same_module <- length(
  unique(
    core_module_localization$ModuleColor
  )
) == 1

cat("\nDo AOC3, MS4A7 and IL7R share one module?\n")
print(same_module)


############################################################
# 19. SAVE RESULTS TO OBJECTS
############################################################

Step14F1_core_annotation <- core_annotation

Step14F1_core_module_localization <-
  core_module_localization

Step14F1_module_mapping <-
  module_mapping

Step14F1_core_module_table <-
  core_module_table

Step14F1_same_module <-
  same_module


############################################################
# 20. SAVE TO CONSOLE-FRIENDLY SUMMARY
############################################################

cat("\n============================================\n")
cat("STEP 14F-1 SUMMARY\n")
cat("============================================\n\n")

cat(
  "Core genes analyzed: ",
  paste(core_genes, collapse = ", "),
  "\n",
  sep = ""
)

cat(
  "Genes successfully located in WGCNA: ",
  nrow(core_module_localization),
  "/",
  length(core_genes),
  "\n",
  sep = ""
)

cat("\nModule localization:\n")

print(
  core_module_localization[
    ,
    c(
      "Gene",
      "Ensembl_ID",
      "ModuleColor",
      "Internal_ME",
      "ModuleMapping_r"
    )
  ],
  row.names = FALSE
)

cat("\n============================================\n")
cat("IMPORTANT:\n")
cat("Do NOT proceed to MM / GS integration yet.\n")
cat("First inspect and confirm the output above.\n")
cat("============================================\n")


module_tme_results


############################################################
# WORKFLOW-03
# Step 14F-2:
# Core-module TME and survival integration
#
# Locked core-gene localization:
#   AOC3  -> magenta / ME9
#   MS4A7 -> blue    / ME2
#   IL7R  -> blue    / ME2
#
# Survival-associated reference module:
#   brown / ME3
#
# Purpose:
#   1. Extract locked TME associations for ME2, ME9 and ME3
#   2. Summarize module-level TME associations
#   3. Compare core-gene modules with the brown survival module
#   4. Integrate previously established survival information
#
# IMPORTANT:
#   - Do NOT recompute the 150 module-TME tests.
#   - Do NOT perform MM or GS in this step.
#   - Do NOT infer gene regulation or mechanism.
############################################################


############################################################
# 0. BASIC SETTINGS
############################################################

options(stringsAsFactors = FALSE)

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_module_map <- data.frame(
  Gene = core_genes,
  ModuleColor = c(
    "magenta",
    "blue",
    "blue"
  ),
  Internal_ME = c(
    "ME9",
    "ME2",
    "ME2"
  ),
  stringsAsFactors = FALSE
)

reference_module <- data.frame(
  ModuleColor = "brown",
  Internal_ME = "ME3",
  Role = "survival-associated reference module",
  stringsAsFactors = FALSE
)


############################################################
# 1. CHECK REQUIRED OBJECTS
############################################################

required_objects <- c(
  "module_tme_results",
  "core_module_localization"
)

missing_objects <- required_objects[
  !required_objects %in% ls(envir = .GlobalEnv)
]

if (length(missing_objects) > 0) {
  
  stop(
    paste0(
      "The following required objects are missing:\n",
      paste(
        missing_objects,
        collapse = ", "
      ),
      "\n\nPlease load the Step 14A-14E locked workspace ",
      "and the Step 14F-1 results first."
    )
  )
}


############################################################
# 2. BASIC QC OF LOCKED TME RESULTS
############################################################

cat("\n============================================\n")
cat("STEP 14F-2\n")
cat("CORE-MODULE TME AND SURVIVAL INTEGRATION\n")
cat("============================================\n\n")

cat("module_tme_results dimensions:\n")
print(dim(module_tme_results))

cat("\nExpected: 150 x 6\n")


if (nrow(module_tme_results) != 150) {
  
  stop(
    "module_tme_results does not contain exactly 150 tests. ",
    "Do not proceed without checking the locked Step 14D result."
  )
}


############################################################
# 3. DEFINE THE THREE MODULES OF INTEREST
############################################################

modules_of_interest <- c(
  "ME2",
  "ME9",
  "ME3"
)

module_color_lookup <- data.frame(
  Internal_ME = c(
    "ME2",
    "ME9",
    "ME3"
  ),
  ModuleColor = c(
    "blue",
    "magenta",
    "brown"
  ),
  ModuleRole = c(
    "MS4A7/IL7R core-gene module",
    "AOC3 core-gene module",
    "survival-associated reference module"
  ),
  stringsAsFactors = FALSE
)


############################################################
# 4. EXTRACT TME RESULTS FOR THE THREE MODULES
############################################################

module_tme_selected <- module_tme_results[
  module_tme_results$Module %in%
    modules_of_interest,
  ,
  drop = FALSE
]


############################################################
# 5. ADD COLOR AND BIOLOGICAL ROLE
############################################################

module_tme_selected$ModuleColor <-
  module_color_lookup$ModuleColor[
    match(
      module_tme_selected$Module,
      module_color_lookup$Internal_ME
    )
  ]

module_tme_selected$ModuleRole <-
  module_color_lookup$ModuleRole[
    match(
      module_tme_selected$Module,
      module_color_lookup$Internal_ME
    )
  ]


############################################################
# 6. ORDER RESULTS
############################################################

module_order <- c(
  "ME2",
  "ME9",
  "ME3"
)

trait_order <- c(
  "ImmuneScore",
  "StromalScore",
  "ESTIMATEScore"
)

module_tme_selected$Module <-
  factor(
    module_tme_selected$Module,
    levels = module_order
  )

module_tme_selected$Trait <-
  factor(
    module_tme_selected$Trait,
    levels = trait_order
  )

module_tme_selected <- module_tme_selected[
  order(
    module_tme_selected$Module,
    module_tme_selected$Trait
  ),
  ,
  drop = FALSE
]

module_tme_selected$Module <-
  as.character(module_tme_selected$Module)

module_tme_selected$Trait <-
  as.character(module_tme_selected$Trait)


############################################################
# 7. PRINT THE THREE-MODULE TME RESULTS
############################################################

cat("\n============================================\n")
cat("7. THREE-MODULE TME ASSOCIATIONS\n")
cat("============================================\n\n")

print(
  module_tme_selected,
  row.names = FALSE
)


############################################################
# 8. VERIFY THAT ALL 9 EXPECTED RESULTS ARE PRESENT
############################################################

expected_n <- 3 * 3

if (nrow(module_tme_selected) != expected_n) {
  
  stop(
    paste0(
      "Expected 9 module-TME results, but found ",
      nrow(module_tme_selected),
      "."
    )
  )
}


############################################################
# 9. CREATE A COMPACT MANUSCRIPT-READY TABLE
############################################################

module_tme_compact <- module_tme_selected[
  ,
  c(
    "ModuleColor",
    "Module",
    "ModuleRole",
    "Trait",
    "N",
    "Spearman_rho",
    "P_value",
    "FDR"
  )
]

module_tme_compact$Spearman_rho <-
  round(
    module_tme_compact$Spearman_rho,
    4
  )

module_tme_compact$P_value <-
  signif(
    module_tme_compact$P_value,
    4
  )

module_tme_compact$FDR <-
  signif(
    module_tme_compact$FDR,
    4
  )


############################################################
# 10. CREATE WIDE TME SUMMARY
############################################################

tme_wide <- reshape(
  module_tme_selected[
    ,
    c(
      "Module",
      "ModuleColor",
      "Trait",
      "Spearman_rho",
      "FDR"
    )
  ],
  idvar = c(
    "Module",
    "ModuleColor"
  ),
  timevar = "Trait",
  direction = "wide"
)

rownames(tme_wide) <- NULL


############################################################
# 11. PRINT WIDE SUMMARY
############################################################

cat("\n============================================\n")
cat("11. MODULE-LEVEL TME SUMMARY\n")
cat("============================================\n\n")

print(
  tme_wide,
  row.names = FALSE
)


############################################################
# 12. CALCULATE MEAN ABSOLUTE TME ASSOCIATION
#
# This is descriptive only.
# It is NOT a new statistical test.
############################################################

module_tme_strength <- aggregate(
  abs(Spearman_rho) ~ Module + ModuleColor,
  data = module_tme_selected,
  FUN = mean
)

colnames(module_tme_strength)[3] <-
  "MeanAbsRho"

module_tme_strength <- module_tme_strength[
  order(
    -module_tme_strength$MeanAbsRho
  ),
  ,
  drop = FALSE
]


############################################################
# 13. ADD MODULE ROLE
############################################################

module_tme_strength$ModuleRole <-
  module_color_lookup$ModuleRole[
    match(
      module_tme_strength$Module,
      module_color_lookup$Internal_ME
    )
  ]

module_tme_strength <- module_tme_strength[
  ,
  c(
    "Module",
    "ModuleColor",
    "ModuleRole",
    "MeanAbsRho"
  )
]


############################################################
# 14. PRINT TME STRENGTH
############################################################

cat("\n============================================\n")
cat("14. DESCRIPTIVE TME ASSOCIATION STRENGTH\n")
cat("============================================\n\n")

print(
  module_tme_strength,
  row.names = FALSE
)


############################################################
# 15. EXTRACT CORE-GENE MODULES ONLY
############################################################

core_modules <- unique(
  core_module_map$Internal_ME
)

core_module_tme <- module_tme_selected[
  module_tme_selected$Module %in%
    core_modules,
  ,
  drop = FALSE
]


############################################################
# 16. ADD CORE-GENE LABEL
############################################################

core_module_tme$CoreGenes <- sapply(
  core_module_tme$Module,
  function(me) {
    
    paste(
      core_module_map$Gene[
        core_module_map$Internal_ME == me
      ],
      collapse = ", "
    )
    
  }
)


############################################################
# 17. PRINT CORE-GENE MODULE TME RESULTS
############################################################

cat("\n============================================\n")
cat("17. CORE-GENE MODULE TME ASSOCIATIONS\n")
cat("============================================\n\n")

print(
  core_module_tme[
    ,
    c(
      "CoreGenes",
      "ModuleColor",
      "Module",
      "Trait",
      "Spearman_rho",
      "P_value",
      "FDR"
    )
  ],
  row.names = FALSE
)


############################################################
# 18. BROWN MODULE RESULTS
############################################################

brown_module_tme <- module_tme_selected[
  module_tme_selected$Module == "ME3",
  ,
  drop = FALSE
]

cat("\n============================================\n")
cat("18. BROWN / ME3 TME ASSOCIATIONS\n")
cat("============================================\n\n")

print(
  brown_module_tme[
    ,
    c(
      "ModuleColor",
      "Module",
      "Trait",
      "Spearman_rho",
      "P_value",
      "FDR"
    )
  ],
  row.names = FALSE
)


############################################################
# 19. COMPARE CORE MODULES WITH BROWN
############################################################

# Calculate difference in rho for each TME trait:
#
# ME2 - ME3
# ME9 - ME3
#
# This is descriptive only.
# It does NOT constitute a formal statistical comparison
# between dependent correlations.

compare_tme <- merge(
  module_tme_selected[
    module_tme_selected$Module %in%
      c("ME2", "ME9"),
    c(
      "Module",
      "ModuleColor",
      "Trait",
      "Spearman_rho"
    )
  ],
  brown_module_tme[
    ,
    c(
      "Trait",
      "Spearman_rho"
    )
  ],
  by = "Trait",
  suffixes = c(
    "_CoreModule",
    "_Brown"
  )
)

compare_tme$DeltaRho <-
  compare_tme$Spearman_rho_CoreModule -
  compare_tme$Spearman_rho_Brown


############################################################
# 20. PRINT DESCRIPTIVE COMPARISON
############################################################

cat("\n============================================\n")
cat("20. CORE MODULE vs BROWN MODULE\n")
cat("============================================\n\n")

print(
  compare_tme,
  row.names = FALSE
)


############################################################
# 21. SURVIVAL REFERENCE INFORMATION
############################################################

# Brown / ME3 is the already validated survival-associated
# module from the preceding tumor-only Cox workflow.
#
# No survival model is refitted here.

survival_reference <- data.frame(
  ModuleColor = "brown",
  Internal_ME = "ME3",
  SurvivalRole =
    "Previously validated survival-associated module",
  RefitPerformed = FALSE,
  stringsAsFactors = FALSE
)

cat("\n============================================\n")
cat("21. SURVIVAL REFERENCE MODULE\n")
cat("============================================\n\n")

print(
  survival_reference,
  row.names = FALSE
)


############################################################
# 22. CORE-GENE -> MODULE -> TME INTEGRATION TABLE
############################################################

core_module_integration <- core_module_map

core_module_integration$MeanAbsTME_rho <-
  module_tme_strength$MeanAbsRho[
    match(
      core_module_integration$Internal_ME,
      module_tme_strength$Module
    )
  ]

core_module_integration$TME_FDR05_count <-
  sapply(
    core_module_integration$Internal_ME,
    function(me) {
      
      sum(
        module_tme_results$Module == me &
          module_tme_results$FDR < 0.05
      )
      
    }
  )


############################################################
# 23. PRINT CORE-GENE INTEGRATION SUMMARY
############################################################

cat("\n============================================\n")
cat("23. CORE-GENE -> MODULE -> TME SUMMARY\n")
cat("============================================\n\n")

print(
  core_module_integration,
  row.names = FALSE
)


############################################################
# 24. SAVE ANALYSIS OBJECTS
############################################################

Step14F2_module_tme_selected <-
  module_tme_selected

Step14F2_module_tme_compact <-
  module_tme_compact

Step14F2_module_tme_wide <-
  tme_wide

Step14F2_module_tme_strength <-
  module_tme_strength

Step14F2_core_module_tme <-
  core_module_tme

Step14F2_brown_module_tme <-
  brown_module_tme

Step14F2_compare_core_vs_brown <-
  compare_tme

Step14F2_survival_reference <-
  survival_reference

Step14F2_core_module_integration <-
  core_module_integration


############################################################
# 25. FINAL CONSOLE SUMMARY
############################################################

cat("\n============================================\n")
cat("STEP 14F-2 FINAL SUMMARY\n")
cat("============================================\n\n")

cat(
  "AOC3      -> magenta / ME9\n"
)

cat(
  "MS4A7     -> blue    / ME2\n"
)

cat(
  "IL7R      -> blue    / ME2\n"
)

cat(
  "Reference -> brown   / ME3\n"
)

cat("\nNo TME correlations were recomputed.\n")
cat("No survival model was refitted.\n")
cat("No MM/GS analysis was performed.\n")

cat("\n============================================\n")
cat("STEP 14F-2 COMPLETE\n")
cat("============================================\n")

# ============================================================
# Step 14F-3
# Module-level survival integration
# ME2 / ME9 / ME3
# ============================================================

cat("\n================ Step 14F-3 ================\n")
cat("Module-level survival integration\n")
cat("==============================================\n\n")


# ------------------------------------------------------------
# 1. Check required objects
# ------------------------------------------------------------

required_objects <- c(
  "MEs",
  "clinical_aligned",
  "cox_design_complete",
  "moduleColors",
  "module_tme_results"
)

missing_objects <- required_objects[
  !sapply(required_objects, exists)
]

if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

cat("All required objects are available.\n\n")

# ------------------------------------------------------------
# 2. Check module eigengenes
# ------------------------------------------------------------

target_modules <- c(
  "ME2",
  "ME9",
  "ME3"
)

target_modules %in% colnames(MEs)

cat("\nTarget module eigengenes:\n")
print(
  colnames(MEs)[
    colnames(MEs) %in% target_modules
  ]
)

cat("\nDimensions of MEs:\n")
print(dim(MEs))

cat("\nNumber of NA values:\n")
print(
  sapply(
    MEs[, target_modules, drop = FALSE],
    function(x) sum(is.na(x))
  )
)

# ------------------------------------------------------------
# 3. Sample alignment
# ------------------------------------------------------------

cat("\n================ Sample alignment ================\n")

cat(
  "MEs samples:",
  nrow(MEs),
  "\n"
)

cat(
  "Clinical samples:",
  nrow(clinical_aligned),
  "\n"
)

cat(
  "cox_design_complete samples:",
  nrow(cox_design_complete),
  "\n"
)

cat(
  "MEs rownames identical to clinical_aligned:",
  identical(
    rownames(MEs),
    rownames(clinical_aligned)
  ),
  "\n"
)

cat(
  "MEs rownames identical to cox_design_complete:",
  identical(
    rownames(MEs),
    rownames(cox_design_complete)
  ),
  "\n"
)

# ============================================================
# Step 14F-3A
# Sample-ID alignment QC
# ============================================================

cat("\n================ Step 14F-3A ================\n")
cat("Sample-ID alignment QC\n")
cat("==============================================\n\n")


# ------------------------------------------------------------
# 1. Basic ID inspection
# ------------------------------------------------------------

cat("First 10 MEs sample IDs:\n")
print(
  head(rownames(MEs), 10)
)

cat("\nFirst 10 clinical_aligned sample IDs:\n")
print(
  head(rownames(clinical_aligned), 10)
)

cat("\nFirst 10 cox_design_complete sample IDs:\n")
print(
  head(rownames(cox_design_complete), 10)
)

# ------------------------------------------------------------
# 2. Number of unique IDs
# ------------------------------------------------------------

cat("\n================ Unique ID QC ================\n")

cat(
  "MEs unique IDs:",
  length(unique(rownames(MEs))),
  "\n"
)

cat(
  "clinical_aligned unique IDs:",
  length(unique(rownames(clinical_aligned))),
  "\n"
)

cat(
  "cox_design_complete unique IDs:",
  length(unique(rownames(cox_design_complete))),
  "\n"
)

cat("\nDuplicated IDs:\n")

cat(
  "MEs:",
  sum(duplicated(rownames(MEs))),
  "\n"
)

cat(
  "clinical_aligned:",
  sum(duplicated(rownames(clinical_aligned))),
  "\n"
)

cat(
  "cox_design_complete:",
  sum(duplicated(rownames(cox_design_complete))),
  "\n"
)

# ------------------------------------------------------------
# 3. Set comparison
# ------------------------------------------------------------

ids_MEs <- rownames(MEs)

ids_clinical <- rownames(clinical_aligned)

ids_cox <- rownames(cox_design_complete)


cat("\n================ Set comparison ================\n")


cat(
  "MEs ∩ clinical_aligned:",
  length(
    intersect(
      ids_MEs,
      ids_clinical
    )
  ),
  "\n"
)

cat(
  "MEs ∩ cox_design_complete:",
  length(
    intersect(
      ids_MEs,
      ids_cox
    )
  ),
  "\n"
)

cat(
  "clinical_aligned ∩ cox_design_complete:",
  length(
    intersect(
      ids_clinical,
      ids_cox
    )
  ),
  "\n"
)

#4. Identify exactly which 10 samples caused the shift from 513 to 503.

# ------------------------------------------------------------
# 4. IDs present in clinical_aligned but absent from cox
# ------------------------------------------------------------

clinical_not_in_cox <- setdiff(
  ids_clinical,
  ids_cox
)

cat(
  "\nSamples in clinical_aligned but NOT in cox_design_complete:",
  length(clinical_not_in_cox),
  "\n"
)

print(
  clinical_not_in_cox
)

#Then, check in the reverse direction:
# ------------------------------------------------------------
# 5. IDs present in cox but absent from clinical
# ------------------------------------------------------------

cox_not_in_clinical <- setdiff(
  ids_cox,
  ids_clinical
)

cat(
  "\nSamples in cox_design_complete but NOT in clinical_aligned:",
  length(cox_not_in_clinical),
  "\n"
)

print(
  cox_not_in_clinical
)
#Check again：
# ------------------------------------------------------------
# 6. IDs present in MEs but absent from clinical
# ------------------------------------------------------------

MEs_not_in_clinical <- setdiff(
  ids_MEs,
  ids_clinical
)

cat(
  "\nSamples in MEs but NOT in clinical_aligned:",
  length(MEs_not_in_clinical),
  "\n"
)

print(
  MEs_not_in_clinical
)


# ------------------------------------------------------------
# 7. IDs present in clinical but absent from MEs
# ------------------------------------------------------------

clinical_not_in_MEs <- setdiff(
  ids_clinical,
  ids_MEs
)

cat(
  "\nSamples in clinical_aligned but NOT in MEs:",
  length(clinical_not_in_MEs),
  "\n"
)

print(
  clinical_not_in_MEs
)


# ------------------------------------------------------------
# 8. Is cox_design_complete a subset of clinical_aligned?
# ------------------------------------------------------------

cat(
  "\ncox_design_complete is subset of clinical_aligned:",
  all(
    ids_cox %in% ids_clinical
  ),
  "\n"
)

cat(
  "cox_design_complete is subset of MEs:",
  all(
    ids_cox %in% ids_MEs
  ),
  "\n"
)

# ------------------------------------------------------------
# 9. Same IDs but different order?
# ------------------------------------------------------------

cat(
  "\nSame MEs and clinical IDs:",
  setequal(
    ids_MEs,
    ids_clinical
  ),
  "\n"
)

cat(
  "Same MEs and cox IDs:",
  setequal(
    ids_MEs,
    ids_cox
  ),
  "\n"
)

cat(
  "Same clinical and cox IDs:",
  setequal(
    ids_clinical,
    ids_cox
  ),
  "\n"
)

# ============================================================
# Step 14F-3B
# Recover and validate sample IDs
# ============================================================

cat("\n================ Step 14F-3B ================\n")
cat("Recover and validate sample IDs\n")
cat("==============================================\n\n")

cat("clinical_aligned dimensions:\n")
print(dim(clinical_aligned))

cat("\nclinical_aligned column names:\n")
print(colnames(clinical_aligned))

# ------------------------------------------------------------
# 2. Inspect first rows of clinical_aligned
# ------------------------------------------------------------

cat("\nFirst 5 rows of clinical_aligned:\n")

print(
  clinical_aligned[
    1:min(5, nrow(clinical_aligned)),
    ,
    drop = FALSE
  ]
)

cat("\ncox_design_complete column names:\n")
print(colnames(cox_design_complete))

cat("\nFirst 5 rows of cox_design_complete:\n")

print(
  cox_design_complete[
    1:min(5, nrow(cox_design_complete)),
    ,
    drop = FALSE
  ]
)

# ------------------------------------------------------------
# 3. Identify possible sample-ID columns
# ------------------------------------------------------------

candidate_id_cols <- colnames(clinical_aligned)[
  grepl(
    "sample|patient|barcode|submitter|id",
    colnames(clinical_aligned),
    ignore.case = TRUE
  )
]

cat(
  "\nPossible ID columns in clinical_aligned:\n"
)

print(candidate_id_cols)


candidate_id_cols_cox <- colnames(cox_design_complete)[
  grepl(
    "sample|patient|barcode|submitter|id",
    colnames(cox_design_complete),
    ignore.case = TRUE
  )
]

cat(
  "\nPossible ID columns in cox_design_complete:\n"
)

print(candidate_id_cols_cox)


# ------------------------------------------------------------
# 4. Check attributes and structure
# ------------------------------------------------------------

cat("\nStructure of clinical_aligned:\n")
str(
  clinical_aligned,
  max.level = 1
)

cat("\nAttributes of clinical_aligned:\n")
print(attributes(clinical_aligned))


cat("\nStructure of cox_design_complete:\n")

str(
  cox_design_complete,
  max.level = 1
)

cat("\nAttributes of cox_design_complete:\n")

print(
  attributes(cox_design_complete)
)

# ============================================================
# Step 14F-3C
# Explicitly validate TCGA Sample IDs against WGCNA MEs
# ============================================================

cat("\n================ Step 14F-3C ================\n")
cat("Explicit TCGA Sample ID validation\n")
cat("==============================================\n\n")


# ------------------------------------------------------------
# 1. Extract explicit sample IDs
# ------------------------------------------------------------

ME_sample_ids <- rownames(MEs)

clinical_sample_ids <- clinical_aligned$Sample

cox_sample_ids <- cox_design_complete$Sample


cat("Number of MEs sample IDs:\n")
print(length(ME_sample_ids))

cat("\nNumber of clinical Sample IDs:\n")
print(length(clinical_sample_ids))

cat("\nNumber of Cox Sample IDs:\n")
print(length(cox_sample_ids))


# ------------------------------------------------------------
# 2. Check uniqueness
# ------------------------------------------------------------

cat("\nUnique MEs IDs:\n")
print(length(unique(ME_sample_ids)))

cat("\nUnique clinical Sample IDs:\n")
print(length(unique(clinical_sample_ids)))

cat("\nUnique Cox Sample IDs:\n")
print(length(unique(cox_sample_ids)))


# ------------------------------------------------------------
# 3. Check whether clinical Sample IDs exactly match MEs
# ------------------------------------------------------------

cat("\nMEs IDs == clinical Sample IDs as sets:\n")
print(setequal(ME_sample_ids, clinical_sample_ids))

cat("\nMEs IDs == clinical Sample IDs in exact order:\n")
print(identical(ME_sample_ids, clinical_sample_ids))


# ------------------------------------------------------------
# 4. If not identical, identify differences
# ------------------------------------------------------------

cat("\nMEs IDs not found in clinical Sample IDs:\n")
print(
  setdiff(
    ME_sample_ids,
    clinical_sample_ids
  )
)

cat("\nClinical Sample IDs not found in MEs IDs:\n")
print(
  setdiff(
    clinical_sample_ids,
    ME_sample_ids
  )
)


# ------------------------------------------------------------
# 5. Check exact order position-by-position
# ------------------------------------------------------------

same_position <- ME_sample_ids == clinical_sample_ids

cat("\nNumber of identical positions:\n")
print(sum(same_position))

cat("\nNumber of different positions:\n")
print(sum(!same_position))


# ------------------------------------------------------------
# 6. Check first and last 10 samples
# ------------------------------------------------------------

cat("\nFirst 10 MEs IDs:\n")
print(ME_sample_ids[1:10])

cat("\nFirst 10 clinical Sample IDs:\n")
print(clinical_sample_ids[1:10])

cat("\nLast 10 MEs IDs:\n")
print(tail(ME_sample_ids, 10))

cat("\nLast 10 clinical Sample IDs:\n")
print(tail(clinical_sample_ids, 10))


# ------------------------------------------------------------
# 7. Verify cox_design_complete against clinical Sample IDs
# ------------------------------------------------------------

cat("\nCox IDs are subset of clinical IDs:\n")
print(
  all(
    cox_sample_ids %in%
      clinical_sample_ids
  )
)

cat("\nCox IDs are subset of MEs IDs:\n")
print(
  all(
    cox_sample_ids %in%
      ME_sample_ids
  )
)


# ------------------------------------------------------------
# 8. Check Cox sample order against the corresponding
#    original clinical rows
# ------------------------------------------------------------

cox_original_position <- match(
  cox_sample_ids,
  clinical_sample_ids
)

cat("\nCox original clinical positions:\n")
print(
  head(
    cox_original_position,
    20
  )
)

cat("\nAre Cox samples in the same order as their original\nclinical_aligned positions?\n")

print(
  identical(
    cox_original_position,
    sort(cox_original_position)
  )
)


# ------------------------------------------------------------
# 9. Verify Cox samples correspond to the same MEs rows
# ------------------------------------------------------------

ME_position <- match(
  cox_sample_ids,
  ME_sample_ids
)

cat("\nFirst 20 Cox -> MEs positions:\n")
print(
  head(
    ME_position,
    20
  )
)

cat("\nAre all Cox sample IDs uniquely matched to MEs?\n")
print(
  all(!is.na(ME_position)) &&
    length(unique(ME_position)) == length(ME_position)
)


# ------------------------------------------------------------
# 10. Final diagnostic summary
# ------------------------------------------------------------

cat("\n================ FINAL DIAGNOSTIC ================\n")

cat(
  "MEs vs clinical set equality: ",
  setequal(ME_sample_ids, clinical_sample_ids),
  "\n",
  sep = ""
)

cat(
  "MEs vs clinical exact order: ",
  identical(ME_sample_ids, clinical_sample_ids),
  "\n",
  sep = ""
)

cat(
  "Cox samples all found in clinical: ",
  all(cox_sample_ids %in% clinical_sample_ids),
  "\n",
  sep = ""
)

cat(
  "Cox samples all found in MEs: ",
  all(cox_sample_ids %in% ME_sample_ids),
  "\n",
  sep = ""
)

cat(
  "Cox order preserved relative to clinical: ",
  identical(
    cox_original_position,
    sort(cox_original_position)
  ),
  "\n",
  sep = ""
)

cat(
  "Cox samples uniquely matched to MEs: ",
  all(!is.na(ME_position)) &&
    length(unique(ME_position)) == length(ME_position),
  "\n",
  sep = ""
)

# ============================================================
# Step 14F-3D
# Build validated survival dataset for ME2 / ME9 / ME3
# ============================================================

cat("\n================ Step 14F-3D ================\n")
cat("Build validated survival dataset\n")
cat("==============================================\n\n")


# ------------------------------------------------------------
# 1. Define validated sample IDs
# ------------------------------------------------------------

ME_sample_ids <- rownames(MEs)

clinical_sample_ids <- clinical_aligned$Sample

cox_sample_ids <- cox_design_complete$Sample


# ------------------------------------------------------------
# 2. Confirm the three key module columns
# ------------------------------------------------------------

required_ME_cols <- c(
  "ME2",
  "ME9",
  "ME3"
)

cat("Required module columns:\n")
print(required_ME_cols)

cat("\nAvailable MEs columns:\n")
print(required_ME_cols %in% colnames(MEs))


# ------------------------------------------------------------
# 3. Create explicit module data using TCGA Sample IDs
# ------------------------------------------------------------

module_survival_df <- data.frame(
  Sample = ME_sample_ids,
  ME2 = MEs[, "ME2"],
  ME9 = MEs[, "ME9"],
  ME3 = MEs[, "ME3"],
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 4. Create clinical survival data using explicit Sample ID
# ------------------------------------------------------------

clinical_survival_df <- data.frame(
  Sample = clinical_aligned$Sample,
  Age = clinical_aligned$Age,
  Gender = clinical_aligned$Gender,
  Stage = clinical_aligned$pathologic_stage,
  OS_time = clinical_aligned$OS.time,
  OS_status = clinical_aligned$OS.status,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 5. Verify exact ID matching before merge
# ------------------------------------------------------------

cat("\nModule survival IDs == clinical survival IDs:\n")
print(
  identical(
    module_survival_df$Sample,
    clinical_survival_df$Sample
  )
)

cat("\nModule survival ID set equality:\n")
print(
  setequal(
    module_survival_df$Sample,
    clinical_survival_df$Sample
  )
)


# ------------------------------------------------------------
# 6. Merge by explicit TCGA Sample ID
# ------------------------------------------------------------

surv_module_all <- merge(
  module_survival_df,
  clinical_survival_df,
  by = "Sample",
  sort = FALSE
)


# ------------------------------------------------------------
# 7. Restore MEs / clinical original order
# ------------------------------------------------------------

surv_module_all <- surv_module_all[
  match(
    ME_sample_ids,
    surv_module_all$Sample
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 8. Verify final order
# ------------------------------------------------------------

cat("\nFinal survival dataset dimensions:\n")
print(dim(surv_module_all))

cat("\nFinal Sample IDs exactly match MEs order:\n")
print(
  identical(
    surv_module_all$Sample,
    ME_sample_ids
  )
)


# ------------------------------------------------------------
# 9. Convert survival variables
# ------------------------------------------------------------

surv_module_all$OS_event <- ifelse(
  surv_module_all$OS_status == "Dead",
  1,
  0
)


# ------------------------------------------------------------
# 10. Convert categorical variables to factors
# ------------------------------------------------------------

surv_module_all$Gender <- factor(
  surv_module_all$Gender,
  levels = c("female", "male")
)

surv_module_all$Stage <- factor(
  surv_module_all$Stage,
  levels = c(
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IV"
  )
)


# ------------------------------------------------------------
# 11. Check complete cases for the planned Cox model
# ------------------------------------------------------------

analysis_vars <- c(
  "Sample",
  "Age",
  "Gender",
  "Stage",
  "OS_time",
  "OS_event",
  "ME2",
  "ME9",
  "ME3"
)

surv_module_complete <- surv_module_all[
  complete.cases(
    surv_module_all[, analysis_vars]
  ),
  analysis_vars,
  drop = FALSE
]


# ------------------------------------------------------------
# 12. Set TCGA Sample ID as rownames
# ------------------------------------------------------------

rownames(surv_module_complete) <-
  surv_module_complete$Sample


# ------------------------------------------------------------
# 13. Basic QC
# ------------------------------------------------------------

cat("\n================ COMPLETE-CASE QC ================\n")

cat("\nComplete-case dimensions:\n")
print(dim(surv_module_complete))

cat("\nNumber of OS events:\n")
print(table(surv_module_complete$OS_event))

cat("\nGender distribution:\n")
print(table(surv_module_complete$Gender))

cat("\nStage distribution:\n")
print(table(surv_module_complete$Stage))


# ------------------------------------------------------------
# 14. Missingness summary
# ------------------------------------------------------------

cat("\nMissing values in all 513 samples:\n")

missing_summary <- sapply(
  surv_module_all[, analysis_vars],
  function(x) sum(is.na(x))
)

print(missing_summary)


# ------------------------------------------------------------
# 15. Compare complete-case IDs with existing
#     cox_design_complete
# ------------------------------------------------------------

cat("\nComplete-case IDs == cox_design_complete IDs:\n")

print(
  setequal(
    rownames(surv_module_complete),
    cox_design_complete$Sample
  )
)

cat("\nComplete-case IDs in identical order:\n")

print(
  identical(
    rownames(surv_module_complete),
    cox_design_complete$Sample
  )
)


# ------------------------------------------------------------
# 16. Verify ME3 equals the previously used ME_brown
# ------------------------------------------------------------

brown_check <- cor(
  surv_module_complete$ME3,
  cox_design_complete$ME_brown,
  method = "pearson"
)

cat("\nPearson correlation: ME3 vs existing ME_brown:\n")
print(brown_check)


# ------------------------------------------------------------
# 17. Final summary
# ------------------------------------------------------------

cat("\n================ FINAL SUMMARY ================\n")

cat(
  "All 513 module/clinical samples aligned: ",
  identical(
    surv_module_all$Sample,
    ME_sample_ids
  ),
  "\n",
  sep = ""
)

cat(
  "Complete-case N: ",
  nrow(surv_module_complete),
  "\n",
  sep = ""
)

cat(
  "Complete-case events: ",
  sum(surv_module_complete$OS_event),
  "\n",
  sep = ""
)

cat(
  "Complete-case non-events: ",
  sum(surv_module_complete$OS_event == 0),
  "\n",
  sep = ""
)

cat(
  "Complete-case IDs identical to previous Cox cohort: ",
  setequal(
    rownames(surv_module_complete),
    cox_design_complete$Sample
  ),
  "\n",
  sep = ""
)

cat(
  "ME3 identical to previous brown module: ",
  isTRUE(
    abs(brown_check - 1) < 1e-12
  ),
  "\n",
  sep = ""
)


# ============================================================
# Step 14F-3E-0
# Diagnose discrepancy between the 503-case and 480-case cohorts
# ============================================================

cat("\n================ Step 14F-3E-0 ================\n")
cat("Diagnose survival cohort discrepancy\n")
cat("==============================================\n\n")


# ------------------------------------------------------------
# 1. Extract the two cohorts
# ------------------------------------------------------------

new_ids <- surv_module_complete$Sample
old_ids <- cox_design_complete$Sample


cat("New cohort N:\n")
print(length(new_ids))

cat("\nPrevious Cox cohort N:\n")
print(length(old_ids))


# ------------------------------------------------------------
# 2. Identify cohort differences
# ------------------------------------------------------------

new_not_old <- setdiff(
  new_ids,
  old_ids
)

old_not_new <- setdiff(
  old_ids,
  new_ids
)


cat("\nSamples in new cohort but NOT previous Cox cohort:\n")
print(length(new_not_old))

print(new_not_old)


cat("\nSamples in previous Cox cohort but NOT new cohort:\n")
print(length(old_not_new))

print(old_not_new)


# ------------------------------------------------------------
# 3. Check whether old cohort is a subset of the
#    513 WGCNA samples
# ------------------------------------------------------------

cat("\nPrevious Cox cohort is subset of MEs:\n")
print(
  all(old_ids %in% rownames(MEs))
)


# ------------------------------------------------------------
# 4. Check missingness in clinical variables specifically
# ------------------------------------------------------------

clinical_check <- data.frame(
  Sample = clinical_aligned$Sample,
  Age = clinical_aligned$Age,
  Gender = clinical_aligned$Gender,
  Stage = clinical_aligned$pathologic_stage,
  OS_time = clinical_aligned$OS.time,
  OS_status = clinical_aligned$OS.status,
  stringsAsFactors = FALSE
)


cat("\nMissingness in clinical variables:\n")

print(
  sapply(
    clinical_check[, -1],
    function(x) sum(is.na(x))
  )
)


# ------------------------------------------------------------
# 5. Show samples with ANY missing clinical variable
# ------------------------------------------------------------

clinical_missing_any <- !complete.cases(
  clinical_check[, -1]
)


cat("\nNumber of samples with any missing clinical variable:\n")
print(sum(clinical_missing_any))


cat("\nSamples with any missing clinical variable:\n")

print(
  clinical_check[
    clinical_missing_any,
    ,
    drop = FALSE
  ]
)


# ------------------------------------------------------------
# 6. Compare clinical variables between the two cohorts
# ------------------------------------------------------------

old_clinical <- clinical_check[
  match(
    old_ids,
    clinical_check$Sample
  ),
  ,
  drop = FALSE
]


cat("\nMissingness among previous 503-case Cox cohort:\n")

print(
  sapply(
    old_clinical[, -1],
    function(x) sum(is.na(x))
  )
)


# ------------------------------------------------------------
# 7. Compare Stage definitions directly
# ------------------------------------------------------------

stage_comparison <- data.frame(
  Sample = old_ids,
  Stage_old = cox_design_complete$Stage_trait,
  Stage_new = clinical_check$Stage[
    match(
      old_ids,
      clinical_check$Sample
    )
  ],
  stringsAsFactors = FALSE
)


cat("\nStage comparison in previous Cox cohort:\n")

print(
  table(
    stage_comparison$Stage_old,
    stage_comparison$Stage_new,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 8. Compare OS time
# ------------------------------------------------------------

os_comparison <- data.frame(
  Sample = old_ids,
  OS_old = cox_design_complete$OS_time,
  OS_new = clinical_check$OS_time[
    match(
      old_ids,
      clinical_check$Sample
    )
  ],
  stringsAsFactors = FALSE
)


cat("\nOS time differences in previous Cox cohort:\n")

print(
  sum(
    os_comparison$OS_old != os_comparison$OS_new,
    na.rm = TRUE
  )
)


# ------------------------------------------------------------
# 9. Compare OS event
# ------------------------------------------------------------

os_event_comparison <- data.frame(
  Sample = old_ids,
  Event_old = cox_design_complete$OS_event,
  Status_new = clinical_check$OS_status[
    match(
      old_ids,
      clinical_check$Sample
    )
  ],
  Event_new = ifelse(
    clinical_check$OS_status[
      match(
        old_ids,
        clinical_check$Sample
      )
    ] == "Dead",
    1,
    0
  ),
  stringsAsFactors = FALSE
)


cat("\nOS event differences in previous Cox cohort:\n")

print(
  sum(
    os_event_comparison$Event_old !=
      os_event_comparison$Event_new,
    na.rm = TRUE
  )
)


# ------------------------------------------------------------
# 10. Compare Age
# ------------------------------------------------------------

age_comparison <- data.frame(
  Sample = old_ids,
  Age_old = cox_design_complete$Age_trait,
  Age_new = clinical_check$Age[
    match(
      old_ids,
      clinical_check$Sample
    )
  ],
  stringsAsFactors = FALSE
)


cat("\nAge differences in previous Cox cohort:\n")

print(
  sum(
    age_comparison$Age_old != age_comparison$Age_new,
    na.rm = TRUE
  )
)


# ------------------------------------------------------------
# 11. Final diagnostic
# ------------------------------------------------------------

cat("\n================ FINAL DIAGNOSTIC ================\n")

cat(
  "New cohort N: ",
  length(new_ids),
  "\n",
  sep = ""
)

cat(
  "Previous Cox cohort N: ",
  length(old_ids),
  "\n",
  sep = ""
)

cat(
  "Old cohort missing in new cohort: ",
  length(old_not_new),
  "\n",
  sep = ""
)

cat(
  "New cohort missing from old cohort: ",
  length(new_not_old),
  "\n",
  sep = ""
)

# ============================================================
# Step 14F-3E-1
# Trace the origin of Stage_trait
# DO NOT modify any existing object
# ============================================================

cat("\n================ Step 14F-3E-1 ================\n")
cat("Trace Stage_trait against original clinical T/N/M\n")
cat("=================================================\n\n")

# ------------------------------------------------------------
# 1. Build clinical reference using Sample ID
# ------------------------------------------------------------

clinical_stage_trace <- data.frame(
  Sample = clinical_aligned$Sample,
  Stage_raw = clinical_aligned$pathologic_stage,
  T_raw = clinical_aligned$pathologic_t,
  N_raw = clinical_aligned$pathologic_n,
  M_raw = clinical_aligned$pathologic_m,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Match to previous 503-case cohort
# ------------------------------------------------------------

stage_trace <- data.frame(
  Sample = cox_design_complete$Sample,
  
  Stage_old = cox_design_complete$Stage_trait,
  
  Stage_raw = clinical_stage_trace$Stage_raw[
    match(
      cox_design_complete$Sample,
      clinical_stage_trace$Sample
    )
  ],
  
  T_raw = clinical_stage_trace$T_raw[
    match(
      cox_design_complete$Sample,
      clinical_stage_trace$Sample
    )
  ],
  
  N_raw = clinical_stage_trace$N_raw[
    match(
      cox_design_complete$Sample,
      clinical_stage_trace$Sample
    )
  ],
  
  M_raw = clinical_stage_trace$M_raw[
    match(
      cox_design_complete$Sample,
      clinical_stage_trace$Sample
    )
  ],
  
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 3. Overall distributions
# ------------------------------------------------------------

cat("Stage_old distribution:\n")
print(table(stage_trace$Stage_old, useNA = "ifany"))

cat("\nStage_raw distribution:\n")
print(table(stage_trace$Stage_raw, useNA = "ifany"))

cat("\nT_raw distribution:\n")
print(table(stage_trace$T_raw, useNA = "ifany"))

cat("\nN_raw distribution:\n")
print(table(stage_trace$N_raw, useNA = "ifany"))

cat("\nM_raw distribution:\n")
print(table(stage_trace$M_raw, useNA = "ifany"))

# ------------------------------------------------------------
# 4. Focus specifically on old Stage IV
# ------------------------------------------------------------

cat("\n================ OLD STAGE IV ================\n")

old_stage_IV <- stage_trace[
  stage_trace$Stage_old == "Stage IV",
  ,
  drop = FALSE
]

cat("Number of old Stage IV samples:\n")
print(nrow(old_stage_IV))

cat("\nClinical variables of old Stage IV samples:\n")
print(old_stage_IV)

# ------------------------------------------------------------
# 5. Cross-tabulation of old Stage vs raw T/N/M
# ------------------------------------------------------------

cat("\n================ STAGE IV T/N/M =================\n")

cat("\nStage_old vs T_raw:\n")
print(
  table(
    stage_trace$Stage_old,
    stage_trace$T_raw,
    useNA = "ifany"
  )
)

cat("\nStage_old vs N_raw:\n")
print(
  table(
    stage_trace$Stage_old,
    stage_trace$N_raw,
    useNA = "ifany"
  )
)

cat("\nStage_old vs M_raw:\n")
print(
  table(
    stage_trace$Stage_old,
    stage_trace$M_raw,
    useNA = "ifany"
  )
)

# ------------------------------------------------------------
# 6. Samples where old Stage and raw Stage differ
# ------------------------------------------------------------

stage_mismatch <- stage_trace[
  !is.na(stage_trace$Stage_old) &
    !is.na(stage_trace$Stage_raw) &
    stage_trace$Stage_old != stage_trace$Stage_raw,
  ,
  drop = FALSE
]

cat("\n================ STAGE MISMATCH ================\n")

cat("Number of Stage_old vs Stage_raw mismatches:\n")
print(nrow(stage_mismatch))

cat("\nMismatch table:\n")
print(
  table(
    stage_mismatch$Stage_old,
    stage_mismatch$Stage_raw,
    useNA = "ifany"
  )
)

cat("\nMismatch samples:\n")
print(stage_mismatch)

# ------------------------------------------------------------
# 7. Final diagnostic
# ------------------------------------------------------------

cat("\n================ FINAL DIAGNOSTIC ================\n")

cat(
  "Old Stage IV N = ",
  sum(stage_trace$Stage_old == "Stage IV", na.rm = TRUE),
  "\n",
  sep = ""
)

cat(
  "Old Stage IV with raw Stage I = ",
  sum(
    stage_trace$Stage_old == "Stage IV" &
      stage_trace$Stage_raw == "Stage I",
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Old Stage IV with M1 = ",
  sum(
    stage_trace$Stage_old == "Stage IV" &
      stage_trace$M_raw == "M1",
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

# ============================================================
# Step 14F-3E-2
# Inspect how Stage_trait may have been constructed
# DO NOT modify any existing object
# ============================================================

cat("\n================ Step 14F-3E-2 ================\n")
cat("Inspect Stage_trait construction and related objects\n")
cat("====================================================\n\n")

# ------------------------------------------------------------
# 1. Inspect cox_design_complete structure
# ------------------------------------------------------------

cat("cox_design_complete structure:\n")
str(cox_design_complete)

# ------------------------------------------------------------
# 2. Inspect Stage_trait levels and attributes
# ------------------------------------------------------------

cat("\nStage_trait class:\n")
print(class(cox_design_complete$Stage_trait))

cat("\nStage_trait levels:\n")
print(levels(cox_design_complete$Stage_trait))

cat("\nStage_trait attributes:\n")
print(attributes(cox_design_complete$Stage_trait))

# ------------------------------------------------------------
# 3. Inspect all objects in the current environment
#    whose names contain Stage / stage / clinical / cox
# ------------------------------------------------------------

all_objects <- ls(envir = .GlobalEnv)

stage_related_objects <- all_objects[
  grepl(
    "stage|Stage|clinical|Clinical|cox|Cox",
    all_objects
  )
]

cat("\nStage/clinical/Cox-related objects in Global Environment:\n")
print(stage_related_objects)

# ------------------------------------------------------------
# 4. Inspect dimensions/classes of potentially relevant objects
# ------------------------------------------------------------

cat("\nPotentially relevant objects:\n")

for (obj_name in stage_related_objects) {
  
  obj <- get(obj_name, envir = .GlobalEnv)
  
  cat(
    "\n--- ",
    obj_name,
    " ---\n",
    sep = ""
  )
  
  cat("Class: ")
  print(class(obj))
  
  if (is.data.frame(obj) || is.matrix(obj)) {
    cat("Dimensions: ")
    print(dim(obj))
    
    if (is.data.frame(obj)) {
      cat("Columns:\n")
      print(names(obj))
    }
  }
  
  if (is.factor(obj)) {
    cat("Levels:\n")
    print(levels(obj))
  }
}

# ------------------------------------------------------------
# 5. Search current function/source history if available
# ------------------------------------------------------------

cat("\nSearch history for Stage_trait / pathologic_stage:\n")

if (exists(".Last.value")) {
  cat("`.Last.value` exists.\n")
}

cat("\nUse history() manually if needed:\n")
cat("history(max.show = 500)\n")

cat("\n====================================================\n")
cat("END Step 14F-3E-2\n")

# ============================================================
# Step 14F-3E-3
# Validate clinical_cox as the original survival clinical source
# DO NOT modify any existing object
# ============================================================

cat("\n================ Step 14F-3E-3 ================\n")
cat("Validate clinical_cox as primary survival clinical source\n")
cat("=========================================================\n\n")

# ------------------------------------------------------------
# 1. Dimensions and columns
# ------------------------------------------------------------

cat("clinical_cox dimensions:\n")
print(dim(clinical_cox))

cat("\nclinical_cox columns:\n")
print(names(clinical_cox))

# ------------------------------------------------------------
# 2. Check Sample IDs
# ------------------------------------------------------------

cat("\nclinical_cox Sample ID summary:\n")
cat("N = ", length(clinical_cox$Sample), "\n", sep = "")
cat(
  "Unique N = ",
  length(unique(clinical_cox$Sample)),
  "\n",
  sep = ""
)

cat(
  "All clinical_cox samples in MEs = ",
  all(clinical_cox$Sample %in% rownames(MEs)),
  "\n",
  sep = ""
)

cat(
  "Same Sample IDs/order as clinical_aligned = ",
  identical(
    clinical_cox$Sample,
    clinical_aligned$Sample
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 3. Compare clinical_cox with cox_design_complete
# ------------------------------------------------------------

old_ids <- cox_design_complete$Sample

clinical_cox_old <- clinical_cox[
  match(
    old_ids,
    clinical_cox$Sample
  ),
  ,
  drop = FALSE
]

cat("\n================ Compare clinical_cox vs cox_design_complete ================\n")

cat("\nSample ID differences:\n")
print(
  sum(
    clinical_cox_old$Sample !=
      cox_design_complete$Sample,
    na.rm = TRUE
  )
)

cat("\nAge_trait differences:\n")
print(
  sum(
    clinical_cox_old$Age_trait !=
      cox_design_complete$Age_trait,
    na.rm = TRUE
  )
)

cat("\nGender_trait differences:\n")
print(
  sum(
    clinical_cox_old$Gender_trait !=
      cox_design_complete$Gender_trait,
    na.rm = TRUE
  )
)

cat("\nStage_trait differences:\n")
print(
  sum(
    clinical_cox_old$Stage_trait !=
      cox_design_complete$Stage_trait,
    na.rm = TRUE
  )
)

cat("\nT_trait differences:\n")
print(
  sum(
    clinical_cox_old$T_trait !=
      cox_design_complete$T_trait,
    na.rm = TRUE
  )
)

cat("\nN_trait differences:\n")
print(
  sum(
    clinical_cox_old$N_trait !=
      cox_design_complete$N_trait,
    na.rm = TRUE
  )
)

cat("\nM_trait differences:\n")
print(
  sum(
    clinical_cox_old$M_trait !=
      cox_design_complete$M_trait,
    na.rm = TRUE
  )
)

cat("\nOS_time differences:\n")
print(
  sum(
    clinical_cox_old$OS_time !=
      cox_design_complete$OS_time,
    na.rm = TRUE
  )
)

cat("\nOS_event differences:\n")
print(
  sum(
    clinical_cox_old$OS_event !=
      cox_design_complete$OS_event,
    na.rm = TRUE
  )
)

# ------------------------------------------------------------
# 4. Missingness in clinical_cox
# ------------------------------------------------------------

cat("\n================ Missingness in clinical_cox ================\n")

print(
  sapply(
    clinical_cox[, -1],
    function(x) sum(is.na(x))
  )
)

# ------------------------------------------------------------
# 5. Distribution of all harmonized clinical variables
# ------------------------------------------------------------

cat("\n================ Harmonized clinical distributions ================\n")

cat("\nAge_trait:\n")
print(summary(clinical_cox$Age_trait))

cat("\nGender_trait:\n")
print(table(
  clinical_cox$Gender_trait,
  useNA = "ifany"
))

cat("\nStage_trait:\n")
print(table(
  clinical_cox$Stage_trait,
  useNA = "ifany"
))

cat("\nT_trait:\n")
print(table(
  clinical_cox$T_trait,
  useNA = "ifany"
))

cat("\nN_trait:\n")
print(table(
  clinical_cox$N_trait,
  useNA = "ifany"
))

cat("\nM_trait:\n")
print(table(
  clinical_cox$M_trait,
  useNA = "ifany"
))

cat("\nOS_event:\n")
print(table(
  clinical_cox$OS_event,
  useNA = "ifany"
))

# ------------------------------------------------------------
# 6. Compare clinical_cox with clinical_WGCNA_traits
# ------------------------------------------------------------

cat("\n================ clinical_cox vs clinical_WGCNA_traits ================\n")

clinical_wgcna_match <- clinical_WGCNA_traits[
  match(
    clinical_cox$Sample,
    clinical_WGCNA_traits$Sample
  ),
  ,
  drop = FALSE
]

cat("\nSample ID/order identical:\n")
print(
  identical(
    clinical_cox$Sample,
    clinical_wgcna_match$Sample
  )
)

cat("\nAge_trait differences:\n")
print(
  sum(
    clinical_cox$Age_trait !=
      clinical_wgcna_match$Age_trait,
    na.rm = TRUE
  )
)

cat("\nGender_trait differences:\n")
print(
  sum(
    clinical_cox$Gender_trait !=
      clinical_wgcna_match$Gender_trait,
    na.rm = TRUE
  )
)

cat("\nStage_trait differences:\n")
print(
  sum(
    clinical_cox$Stage_trait !=
      clinical_wgcna_match$Stage_trait,
    na.rm = TRUE
  )
)

cat("\nT_trait differences:\n")
print(
  sum(
    clinical_cox$T_trait !=
      clinical_wgcna_match$T_trait,
    na.rm = TRUE
  )
)

cat("\nN_trait differences:\n")
print(
  sum(
    clinical_cox$N_trait !=
      clinical_wgcna_match$N_trait,
    na.rm = TRUE
  )
)

cat("\nM_trait differences:\n")
print(
  sum(
    clinical_cox$M_trait !=
      clinical_wgcna_match$M_trait,
    na.rm = TRUE
  )
)

cat("\nOS_time differences:\n")
print(
  sum(
    clinical_cox$OS_time !=
      clinical_wgcna_match$OS_time,
    na.rm = TRUE
  )
)

cat("\nOS_event differences:\n")
print(
  sum(
    clinical_cox$OS_event !=
      clinical_wgcna_match$OS_event,
    na.rm = TRUE
  )
)

# ------------------------------------------------------------
# 7. Identify the exact 503-case selection rule
# ------------------------------------------------------------

cat("\n================ 503-case selection ================\n")

cat("\nNumber of complete cases using all clinical_cox variables:\n")

cox_required_vars <- c(
  "Age_trait",
  "Gender_trait",
  "Stage_trait",
  "T_trait",
  "N_trait",
  "M_trait",
  "OS_time",
  "OS_event"
)

complete_clinical_cox <- complete.cases(
  clinical_cox[, c("Sample", cox_required_vars)]
)

print(sum(complete_clinical_cox))

cat("\nNumber of complete cases using only variables in the original conventional Cox design:\n")

cox_primary_vars <- c(
  "Age_trait",
  "Gender_trait",
  "Stage_trait",
  "OS_time",
  "OS_event"
)

complete_primary_clinical <- complete.cases(
  clinical_cox[, c("Sample", cox_primary_vars)]
)

print(sum(complete_primary_clinical))

# ------------------------------------------------------------
# 8. Compare exact 503-case IDs
# ------------------------------------------------------------

clinical_cox_503_ids <- clinical_cox$Sample[
  complete_primary_clinical
]

cat("\nExact 503-case match with cox_design_complete:\n")

cat(
  "Same set = ",
  setequal(
    clinical_cox_503_ids,
    cox_design_complete$Sample
  ),
  "\n",
  sep = ""
)

cat(
  "Same order = ",
  identical(
    clinical_cox_503_ids,
    cox_design_complete$Sample
  ),
  "\n",
  sep = ""
)

cat(
  "clinical_cox complete-case N = ",
  length(clinical_cox_503_ids),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 9. Compare exact 485-case final Cox cohort
# ------------------------------------------------------------

cat("\n================ 485-case final Cox cohort ================\n")

if (exists("stage_analysis_data")) {
  
  stage_485_ids <- stage_analysis_data$Sample
  
  cat(
    "stage_analysis_data N = ",
    length(stage_485_ids),
    "\n",
    sep = ""
  )
  
  cat(
    "503-case clinical cohort contains all 485 = ",
    all(
      stage_485_ids %in%
        cox_design_complete$Sample
    ),
    "\n",
    sep = ""
  )
  
} else {
  
  cat("stage_analysis_data not found.\n")
  
}

# ------------------------------------------------------------
# 10. Final diagnostic summary
# ------------------------------------------------------------

cat("\n================ FINAL DIAGNOSTIC ================\n")

cat(
  "clinical_cox N = ",
  nrow(clinical_cox),
  "\n",
  sep = ""
)

cat(
  "cox_design_complete N = ",
  nrow(cox_design_complete),
  "\n",
  sep = ""
)

cat(
  "clinical_cox primary-variable complete cases = ",
  sum(complete_primary_clinical),
  "\n",
  sep = ""
)

cat(
  "clinical_cox == cox_design_complete clinical variables: ",
  all(
    clinical_cox_old$Age_trait ==
      cox_design_complete$Age_trait
  ) &&
    all(
      clinical_cox_old$Gender_trait ==
        cox_design_complete$Gender_trait
    ) &&
    all(
      clinical_cox_old$Stage_trait ==
        cox_design_complete$Stage_trait
    ),
  "\n",
  sep = ""
)


# ============================================================
# Step 14F-3E-4
# Identify the exact 18 samples included in cox_design_complete
# but excluded from clinical_cox conventional complete cases
#
# DO NOT modify any existing object
# ============================================================

cat("\n================ Step 14F-3E-4 ================\n")
cat("Identify exact 503 vs 485 cohort discrepancy\n")
cat("================================================\n\n")

# ------------------------------------------------------------
# 1. Define the two cohorts
# ------------------------------------------------------------

ids_503 <- cox_design_complete$Sample

ids_485 <- clinical_cox$Sample[
  complete_primary_clinical
]

cat("503 cohort N = ", length(ids_503), "\n", sep = "")
cat("485 cohort N = ", length(ids_485), "\n", sep = "")

# ------------------------------------------------------------
# 2. Exact set difference
# ------------------------------------------------------------

extra_503 <- setdiff(
  ids_503,
  ids_485
)

missing_503 <- setdiff(
  ids_485,
  ids_503
)

cat("\nSamples in 503 but NOT in 485:\n")
print(extra_503)

cat(
  "\nN in 503 but NOT in 485 = ",
  length(extra_503),
  "\n",
  sep = ""
)

cat("\nSamples in 485 but NOT in 503:\n")
print(missing_503)

cat(
  "\nN in 485 but NOT in 503 = ",
  length(missing_503),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 3. Inspect the 18 extra samples directly from clinical_cox
# ------------------------------------------------------------

extra_503_clinical <- clinical_cox[
  match(
    extra_503,
    clinical_cox$Sample
  ),
  ,
  drop = FALSE
]

cat("\n================ Clinical data of extra 503 samples ================\n")

print(extra_503_clinical)

# ------------------------------------------------------------
# 4. Missingness pattern of the 18 samples
# ------------------------------------------------------------

cat("\n================ Missingness pattern ================\n")

primary_vars <- c(
  "Age_trait",
  "Gender_trait",
  "Stage_trait",
  "OS_time",
  "OS_event"
)

extra_missing <- sapply(
  extra_503_clinical[, primary_vars, drop = FALSE],
  function(x) is.na(x)
)

print(extra_missing)

cat("\nNumber of missing values per variable among extra 503 samples:\n")

print(
  colSums(extra_missing)
)

# ------------------------------------------------------------
# 5. Missingness pattern by sample
# ------------------------------------------------------------

extra_missing_df <- data.frame(
  Sample = extra_503_clinical$Sample,
  Age_missing = is.na(extra_503_clinical$Age_trait),
  Gender_missing = is.na(extra_503_clinical$Gender_trait),
  Stage_missing = is.na(extra_503_clinical$Stage_trait),
  OS_time_missing = is.na(extra_503_clinical$OS_time),
  OS_event_missing = is.na(extra_503_clinical$OS_event)
)

cat("\nMissingness by sample:\n")

print(extra_missing_df)

# ------------------------------------------------------------
# 6. Compare rownames of cox_design_complete
# ------------------------------------------------------------

cat("\n================ cox_design_complete rownames ================\n")

cat(
  "First 20 rownames:\n"
)

print(
  head(
    rownames(cox_design_complete),
    20
  )
)

cat(
  "\nAre rownames numeric? ",
  all(
    grepl(
      "^[0-9]+$",
      rownames(cox_design_complete)
    )
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 7. Compare Sample IDs with original clinical positions
# ------------------------------------------------------------

clinical_cox_position <- match(
  ids_503,
  clinical_cox$Sample
)

cat("\n================ Original clinical positions ================\n")

print(
  head(
    clinical_cox_position,
    30
  )
)

cat(
  "\nAre 503 samples in original clinical order? ",
  identical(
    clinical_cox_position,
    sort(clinical_cox_position)
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 8. Compare 503 cohort against all five primary variables
# ------------------------------------------------------------

cat("\n================ 503 cohort missingness ================\n")

cox_503_from_clinical <- clinical_cox[
  match(
    ids_503,
    clinical_cox$Sample
  ),
  ,
  drop = FALSE
]

print(
  sapply(
    cox_503_from_clinical[, primary_vars, drop = FALSE],
    function(x) sum(is.na(x))
  )
)

# ------------------------------------------------------------
# 9. Compare exact 503 cohort to complete-case definitions
# ------------------------------------------------------------

cat("\n================ Complete-case comparison ================\n")

cc_primary <- complete.cases(
  clinical_cox[, primary_vars]
)

cc_age_stage_os <- complete.cases(
  clinical_cox[
    ,
    c(
      "Age_trait",
      "Stage_trait",
      "OS_time",
      "OS_event"
    )
  ]
)

cc_age_os <- complete.cases(
  clinical_cox[
    ,
    c(
      "Age_trait",
      "OS_time",
      "OS_event"
    )
  ]
)

cat(
  "Age + Gender + Stage + OS complete = ",
  sum(cc_primary),
  "\n",
  sep = ""
)

cat(
  "Age + Stage + OS complete = ",
  sum(cc_age_stage_os),
  "\n",
  sep = ""
)

cat(
  "Age + OS complete = ",
  sum(cc_age_os),
  "\n",
  sep = ""
)

cat(
  "503 cohort == Age + Stage + OS complete? ",
  setequal(
    ids_503,
    clinical_cox$Sample[cc_age_stage_os]
  ),
  "\n",
  sep = ""
)

cat(
  "503 cohort == Age + OS complete? ",
  setequal(
    ids_503,
    clinical_cox$Sample[cc_age_os]
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 10. Final summary
# ------------------------------------------------------------

cat("\n================ FINAL SUMMARY ================\n")

cat(
  "cox_design_complete N = ",
  length(ids_503),
  "\n",
  sep = ""
)

cat(
  "clinical_cox conventional complete-case N = ",
  length(ids_485),
  "\n",
  sep = ""
)

cat(
  "Extra samples in 503 cohort = ",
  length(extra_503),
  "\n",
  sep = ""
)

cat(
  "Samples in 485 but not 503 = ",
  length(missing_503),
  "\n",
  sep = ""
)


# ============================================================
# Step 14F-3E-5
# Final identity validation of the 485-sample survival cohort
#
# DO NOT modify any existing object
# ============================================================

cat("\n================ Step 14F-3E-5 ================\n")
cat("Final validation of the 485-sample survival cohort\n")
cat("=================================================\n\n")

# ------------------------------------------------------------
# 1. Clinical complete-case 485
# ------------------------------------------------------------

primary_vars <- c(
  "Age_trait",
  "Gender_trait",
  "Stage_trait",
  "OS_time",
  "OS_event"
)

clinical_485_ids <- clinical_cox$Sample[
  complete.cases(
    clinical_cox[, primary_vars, drop = FALSE]
  )
]

cat("clinical_cox complete-case N = ",
    length(clinical_485_ids),
    "\n",
    sep = "")

# ------------------------------------------------------------
# 2. Samples actually complete for the Brown Cox model
# ------------------------------------------------------------

brown_required <- c(
  "Sample",
  "ME_brown",
  "Age_trait",
  "Gender_trait",
  "Stage_trait",
  "OS_time",
  "OS_event"
)

brown_model_df <- cox_design_complete[
  complete.cases(
    cox_design_complete[, brown_required, drop = FALSE]
  ),
  ,
  drop = FALSE
]

brown_485_ids <- brown_model_df$Sample

cat(
  "Brown Cox complete-case N = ",
  length(brown_485_ids),
  "\n",
  sep = ""
)

cat(
  "Brown Cox events = ",
  sum(brown_model_df$OS_event),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 3. Compare clinical 485 vs Brown Cox 485
# ------------------------------------------------------------

cat("\n================ clinical 485 vs Brown Cox 485 ================\n")

cat(
  "Same set = ",
  setequal(
    clinical_485_ids,
    brown_485_ids
  ),
  "\n",
  sep = ""
)

cat(
  "Same order = ",
  identical(
    clinical_485_ids,
    brown_485_ids
  ),
  "\n",
  sep = ""
)

cat(
  "Clinical 485 but not Brown 485 = ",
  length(
    setdiff(
      clinical_485_ids,
      brown_485_ids
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Brown 485 but not Clinical 485 = ",
  length(
    setdiff(
      brown_485_ids,
      clinical_485_ids
    )
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 4. Compare stage_analysis_data with the same 485 cohort
# ------------------------------------------------------------

cat("\n================ stage_analysis_data comparison ================\n")

stage_485_ids <- stage_analysis_data$Sample

cat(
  "stage_analysis_data N = ",
  length(stage_485_ids),
  "\n",
  sep = ""
)

cat(
  "stage_analysis_data == clinical 485 (same set) = ",
  setequal(
    stage_485_ids,
    clinical_485_ids
  ),
  "\n",
  sep = ""
)

cat(
  "stage_analysis_data == Brown 485 (same set) = ",
  setequal(
    stage_485_ids,
    brown_485_ids
  ),
  "\n",
  sep = ""
)

cat(
  "Clinical 485 but not stage_analysis_data = ",
  length(
    setdiff(
      clinical_485_ids,
      stage_485_ids
    )
  ),
  "\n",
  sep = ""
)

cat(
  "stage_analysis_data but not clinical 485 = ",
  length(
    setdiff(
      stage_485_ids,
      clinical_485_ids
    )
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 5. Exact identity of Brown Cox cohort
# ------------------------------------------------------------

cat("\n================ Brown cohort identity ================\n")

cat(
  "Brown 485 is subset of cox_design_complete 503 = ",
  all(
    brown_485_ids %in%
      cox_design_complete$Sample
  ),
  "\n",
  sep = ""
)

cat(
  "Clinical 485 is subset of cox_design_complete 503 = ",
  all(
    clinical_485_ids %in%
      cox_design_complete$Sample
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 6. Check event identity
# ------------------------------------------------------------

clinical_485_check <- clinical_cox[
  match(
    brown_485_ids,
    clinical_cox$Sample
  ),
  ,
  drop = FALSE
]

cat("\n================ Event/time identity ================\n")

cat(
  "OS_time identical = ",
  identical(
    clinical_485_check$OS_time,
    brown_model_df$OS_time
  ),
  "\n",
  sep = ""
)

cat(
  "OS_event identical = ",
  identical(
    clinical_485_check$OS_event,
    brown_model_df$OS_event
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 7. Final diagnostic
# ------------------------------------------------------------

cat("\n================ FINAL DIAGNOSTIC ================\n")

if (
  setequal(clinical_485_ids, brown_485_ids) &&
  all(clinical_485_ids %in% cox_design_complete$Sample)
) {
  
  cat(
    "FINAL CONCLUSION: clinical_cox 485 and Brown Cox 485 ",
    "represent the same sample set.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "FINAL CONCLUSION: 485-sample cohort identity requires further investigation.\n"
  )
}


# ============================================================
# Step 14F-3F-1
# Build and validate the exact 485-sample survival dataset
# for ME2 / ME9 integration
#
# DO NOT modify any existing object
# ============================================================

cat("\n================ Step 14F-3F-1 ================\n")
cat("Build exact 485-sample survival dataset for ME2/ME9\n")
cat("====================================================\n\n")

# ------------------------------------------------------------
# 1. Locked 485-sample cohort from Brown Cox
# ------------------------------------------------------------

brown_485_ids_locked <- brown_485_ids

cat(
  "Locked Brown survival cohort N = ",
  length(brown_485_ids_locked),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 2. Extract clinical variables from clinical_cox
# ------------------------------------------------------------

clinical_ME_surv <- clinical_cox[
  match(
    brown_485_ids_locked,
    clinical_cox$Sample
  ),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 3. Extract WGCNA module eigengenes
# ------------------------------------------------------------

ME_ME_surv <- MEs[
  match(
    brown_485_ids_locked,
    rownames(MEs)
  ),
  c(
    "ME2",
    "ME9",
    "ME3"
  ),
  drop = FALSE
]

# ------------------------------------------------------------
# 4. Construct final 485-sample dataset
# ------------------------------------------------------------

surv_ME2_ME9_485 <- data.frame(
  Sample = brown_485_ids_locked,
  ME2 = ME_ME_surv$ME2,
  ME9 = ME_ME_surv$ME9,
  ME3 = ME_ME_surv$ME3,
  Age_trait = clinical_ME_surv$Age_trait,
  Gender_trait = clinical_ME_surv$Gender_trait,
  Stage_trait = clinical_ME_surv$Stage_trait,
  OS_time = clinical_ME_surv$OS_time,
  OS_event = clinical_ME_surv$OS_event,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 5. Basic dimensions
# ------------------------------------------------------------

cat("\n================ Dataset dimensions ================\n")

print(
  dim(
    surv_ME2_ME9_485
  )
)

print(
  names(
    surv_ME2_ME9_485
  )
)

# ------------------------------------------------------------
# 6. Sample ID validation
# ------------------------------------------------------------

cat("\n================ Sample ID validation ================\n")

cat(
  "N unique Sample IDs = ",
  length(
    unique(
      surv_ME2_ME9_485$Sample
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Sample IDs identical to Brown 485 = ",
  identical(
    surv_ME2_ME9_485$Sample,
    brown_485_ids_locked
  ),
  "\n",
  sep = ""
)

cat(
  "All Samples in MEs = ",
  all(
    surv_ME2_ME9_485$Sample %in%
      rownames(MEs)
  ),
  "\n",
  sep = ""
)

cat(
  "All Samples in clinical_cox = ",
  all(
    surv_ME2_ME9_485$Sample %in%
      clinical_cox$Sample
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 7. Missingness
# ------------------------------------------------------------

cat("\n================ Missingness ================\n")

print(
  sapply(
    surv_ME2_ME9_485,
    function(x) sum(is.na(x))
  )
)

# ------------------------------------------------------------
# 8. Event distribution
# ------------------------------------------------------------

cat("\n================ Survival distribution ================\n")

cat(
  "N = ",
  nrow(surv_ME2_ME9_485),
  "\n",
  sep = ""
)

cat(
  "Events = ",
  sum(
    surv_ME2_ME9_485$OS_event
  ),
  "\n",
  sep = ""
)

cat(
  "Non-events = ",
  sum(
    surv_ME2_ME9_485$OS_event == 0
  ),
  "\n",
  sep = ""
)

cat("\nOS_event:\n")

print(
  table(
    surv_ME2_ME9_485$OS_event,
    useNA = "ifany"
  )
)

# ------------------------------------------------------------
# 9. Clinical distributions
# ------------------------------------------------------------

cat("\n================ Clinical variables ================\n")

cat("\nGender_trait:\n")

print(
  table(
    surv_ME2_ME9_485$Gender_trait,
    useNA = "ifany"
  )
)

cat("\nStage_trait:\n")

print(
  table(
    surv_ME2_ME9_485$Stage_trait,
    useNA = "ifany"
  )
)

cat("\nAge_trait:\n")

print(
  summary(
    surv_ME2_ME9_485$Age_trait
  )
)

# ------------------------------------------------------------
# 10. Module summaries
# ------------------------------------------------------------

cat("\n================ Module summaries ================\n")

cat("\nME2:\n")
print(
  summary(
    surv_ME2_ME9_485$ME2
  )
)

cat("\nME9:\n")
print(
  summary(
    surv_ME2_ME9_485$ME9
  )
)

cat("\nME3:\n")
print(
  summary(
    surv_ME2_ME9_485$ME3
  )
)

# ------------------------------------------------------------
# 11. Verify ME3 exactly matches Brown model data
# ------------------------------------------------------------

cat("\n================ ME3 identity check ================\n")

brown_ME3_check <- brown_model_df[
  match(
    surv_ME2_ME9_485$Sample,
    brown_model_df$Sample
  ),
  ,
  drop = FALSE
]

cat(
  "ME3 identical to Brown model ME_brown = ",
  identical(
    surv_ME2_ME9_485$ME3,
    brown_ME3_check$ME_brown
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 12. Final diagnostic
# ------------------------------------------------------------

cat("\n================ FINAL DIAGNOSTIC ================\n")

if (
  nrow(surv_ME2_ME9_485) == 485 &&
  sum(surv_ME2_ME9_485$OS_event) == 178 &&
  identical(
    surv_ME2_ME9_485$Sample,
    brown_485_ids_locked
  ) &&
  all(
    surv_ME2_ME9_485$Sample %in%
    rownames(MEs)
  ) &&
  all(
    surv_ME2_ME9_485$Sample %in%
    clinical_cox$Sample
  )
) {
  
  cat(
    "FINAL CONCLUSION: ME2/ME9 survival dataset is correctly aligned ",
    "to the locked 485-sample Brown Cox cohort.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "FINAL CONCLUSION: alignment requires further investigation.\n"
  )
}


# ============================================================
# Step 14F-3F-2
# Conventional Cox models for ME2 and ME9
#
# Locked cohort:
#   485 samples
#   178 events
#
# Model:
#   Module + Age + Gender + Stage
#
# DO NOT modify any existing object
# ============================================================

cat("\n================ Step 14F-3F-2 ================\n")
cat("Conventional Cox models for ME2 and ME9\n")
cat("Locked cohort: 485 samples / 178 events\n")
cat("================================================\n\n")


# ------------------------------------------------------------
# 1. Confirm locked cohort before modeling
# ------------------------------------------------------------

stopifnot(
  nrow(surv_ME2_ME9_485) == 485,
  sum(surv_ME2_ME9_485$OS_event) == 178,
  all(!is.na(surv_ME2_ME9_485$OS_time)),
  all(!is.na(surv_ME2_ME9_485$OS_event)),
  all(!is.na(surv_ME2_ME9_485$ME2)),
  all(!is.na(surv_ME2_ME9_485$ME9)),
  all(!is.na(surv_ME2_ME9_485$ME3)),
  all(!is.na(surv_ME2_ME9_485$Age_trait)),
  all(!is.na(surv_ME2_ME9_485$Gender_trait)),
  all(!is.na(surv_ME2_ME9_485$Stage_trait))
)


# ------------------------------------------------------------
# 2. Standardize ME2 / ME9 by 1 SD
#
# This improves HR interpretability:
# HR = change in hazard per 1-SD increase in module eigengene
#
# Existing ME2 / ME9 are NOT modified.
# ------------------------------------------------------------

ME2_sd_locked <- sd(
  surv_ME2_ME9_485$ME2
)

ME9_sd_locked <- sd(
  surv_ME2_ME9_485$ME9
)

surv_ME2_ME9_485$ME2_SD <- (
  surv_ME2_ME9_485$ME2 /
    ME2_sd_locked
)

surv_ME2_ME9_485$ME9_SD <- (
  surv_ME2_ME9_485$ME9 /
    ME9_sd_locked
)


# ------------------------------------------------------------
# 3. Check SD scaling
# ------------------------------------------------------------

cat("\n================ SD scaling ================\n")

cat(
  "ME2 SD = ",
  ME2_sd_locked,
  "\n",
  sep = ""
)

cat(
  "ME9 SD = ",
  ME9_sd_locked,
  "\n",
  sep = ""
)

cat(
  "SD(ME2_SD) = ",
  sd(surv_ME2_ME9_485$ME2_SD),
  "\n",
  sep = ""
)

cat(
  "SD(ME9_SD) = ",
  sd(surv_ME2_ME9_485$ME9_SD),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. ME2 conventional Cox
# ------------------------------------------------------------

ME2_cox_model <- survival::coxph(
  survival::Surv(
    OS_time,
    OS_event
  ) ~
    ME2_SD +
    Age_trait +
    Gender_trait +
    Stage_trait,
  data = surv_ME2_ME9_485,
  ties = "efron"
)


# ------------------------------------------------------------
# 5. ME9 conventional Cox
# ------------------------------------------------------------

ME9_cox_model <- survival::coxph(
  survival::Surv(
    OS_time,
    OS_event
  ) ~
    ME9_SD +
    Age_trait +
    Gender_trait +
    Stage_trait,
  data = surv_ME2_ME9_485,
  ties = "efron"
)


# ------------------------------------------------------------
# 6. Model summaries
# ------------------------------------------------------------

cat("\n================ ME2 Cox model ================\n")

print(
  summary(
    ME2_cox_model
  )
)


cat("\n================ ME9 Cox model ================\n")

print(
  summary(
    ME9_cox_model
  )
)


# ------------------------------------------------------------
# 7. Extract compact HR / CI / P table
# ------------------------------------------------------------

extract_cox_results <- function(
    model,
    module_name
) {
  
  sm <- summary(model)
  
  coef_df <- as.data.frame(
    sm$coefficients
  )
  
  ci_df <- as.data.frame(
    sm$conf.int
  )
  
  out <- data.frame(
    Module = module_name,
    Variable = rownames(coef_df),
    Beta = coef_df[, "coef"],
    HR = ci_df[, "exp(coef)"],
    Lower95CI = ci_df[, "lower .95"],
    Upper95CI = ci_df[, "upper .95"],
    SE = coef_df[, "se(coef)"],
    P = coef_df[, "Pr(>|z|)"],
    row.names = NULL,
    check.names = FALSE
  )
  
  out
}


ME2_ME9_cox_results <- rbind(
  extract_cox_results(
    ME2_cox_model,
    "ME2"
  ),
  extract_cox_results(
    ME9_cox_model,
    "ME9"
  )
)


# ------------------------------------------------------------
# 8. Print compact results
# ------------------------------------------------------------

cat("\n================ Compact Cox results ================\n")

print(
  ME2_ME9_cox_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Extract model-level statistics
# ------------------------------------------------------------

cat("\n================ Model-level statistics ================\n")

cat("\nME2:\n")

cat(
  "N = ",
  nobs(ME2_cox_model),
  "\n",
  sep = ""
)

cat(
  "Events = ",
  ME2_cox_model$nevent,
  "\n",
  sep = ""
)

cat(
  "C-index = ",
  summary(ME2_cox_model)$concordance[1],
  "\n",
  sep = ""
)

cat(
  "LRT P = ",
  summary(ME2_cox_model)$logtest["pvalue"],
  "\n",
  sep = ""
)

cat(
  "AIC = ",
  AIC(ME2_cox_model),
  "\n",
  sep = ""
)


cat("\nME9:\n")

cat(
  "N = ",
  nobs(ME9_cox_model),
  "\n",
  sep = ""
)

cat(
  "Events = ",
  ME9_cox_model$nevent,
  "\n",
  sep = ""
)

cat(
  "C-index = ",
  summary(ME9_cox_model)$concordance[1],
  "\n",
  sep = ""
)

cat(
  "LRT P = ",
  summary(ME9_cox_model)$logtest["pvalue"],
  "\n",
  sep = ""
)

cat(
  "AIC = ",
  AIC(ME9_cox_model),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Primary module results only
# ------------------------------------------------------------

cat("\n================ Primary module results ================\n")

print(
  ME2_ME9_cox_results[
    ME2_ME9_cox_results$Variable %in%
      c("ME2_SD", "ME9_SD"),
    ,
    drop = FALSE
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Final diagnostic
# ------------------------------------------------------------

cat("\n================ FINAL DIAGNOSTIC ================\n")

if (
  nobs(ME2_cox_model) == 485 &&
  nobs(ME9_cox_model) == 485 &&
  ME2_cox_model$nevent == 178 &&
  ME9_cox_model$nevent == 178
) {
  
  cat(
    "FINAL CONCLUSION: ME2 and ME9 conventional Cox models ",
    "were fitted on the identical locked 485-sample / 178-event ",
    "survival cohort.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "FINAL CONCLUSION: model cohort requires further investigation.\n"
  )
}


# ============================================================
# Step 14F-3F-3
# Proportional hazards assessment for ME2 and ME9
#
# Existing Cox models:
#   ME2_cox_model
#   ME9_cox_model
#
# DO NOT refit or modify existing models
# ============================================================

cat("\n================ Step 14F-3F-3 ================\n")
cat("PH assessment for ME2 and ME9\n")
cat("================================================\n\n")


# ------------------------------------------------------------
# 1. PH test for ME2
# ------------------------------------------------------------

ME2_ph_test <- survival::cox.zph(
  ME2_cox_model,
  transform = "km"
)

cat("\n================ ME2 PH test ================\n")

print(
  ME2_ph_test
)


# ------------------------------------------------------------
# 2. PH test for ME9
# ------------------------------------------------------------

ME9_ph_test <- survival::cox.zph(
  ME9_cox_model,
  transform = "km"
)

cat("\n================ ME9 PH test ================\n")

print(
  ME9_ph_test
)


# ------------------------------------------------------------
# 3. Extract compact PH-test tables
# ------------------------------------------------------------

ME2_ph_table <- as.data.frame(
  ME2_ph_test$table
)

ME2_ph_table$Variable <- rownames(
  ME2_ph_table
)

rownames(
  ME2_ph_table
) <- NULL


ME9_ph_table <- as.data.frame(
  ME9_ph_test$table
)

ME9_ph_table$Variable <- rownames(
  ME9_ph_table
)

rownames(
  ME9_ph_table
) <- NULL


# ------------------------------------------------------------
# 4. Print compact tables
# ------------------------------------------------------------

cat("\n================ ME2 PH compact table ================\n")

print(
  ME2_ph_table,
  row.names = FALSE
)


cat("\n================ ME9 PH compact table ================\n")

print(
  ME9_ph_table,
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Extract module-specific PH P values
# ------------------------------------------------------------

ME2_module_ph_p <- ME2_ph_test$table[
  "ME2_SD",
  "p"
]

ME9_module_ph_p <- ME9_ph_test$table[
  "ME9_SD",
  "p"
]


ME2_global_ph_p <- ME2_ph_test$table[
  "GLOBAL",
  "p"
]

ME9_global_ph_p <- ME9_ph_test$table[
  "GLOBAL",
  "p"
]


# ------------------------------------------------------------
# 6. Print primary PH results
# ------------------------------------------------------------

cat("\n================ Primary PH results ================\n")

cat(
  "ME2 module PH P = ",
  ME2_module_ph_p,
  "\n",
  sep = ""
)

cat(
  "ME2 global PH P = ",
  ME2_global_ph_p,
  "\n",
  sep = ""
)

cat(
  "ME9 module PH P = ",
  ME9_module_ph_p,
  "\n",
  sep = ""
)

cat(
  "ME9 global PH P = ",
  ME9_global_ph_p,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Simple diagnostic
# ------------------------------------------------------------

cat("\n================ PH diagnostic ================\n")

cat(
  "ME2 module PH assumption: ",
  ifelse(
    ME2_module_ph_p < 0.05,
    "potential violation",
    "no evidence of violation"
  ),
  "\n",
  sep = ""
)

cat(
  "ME9 module PH assumption: ",
  ifelse(
    ME9_module_ph_p < 0.05,
    "potential violation",
    "no evidence of violation"
  ),
  "\n",
  sep = ""
)

cat(
  "ME2 global PH assumption: ",
  ifelse(
    ME2_global_ph_p < 0.05,
    "potential violation",
    "no evidence of violation"
  ),
  "\n",
  sep = ""
)

cat(
  "ME9 global PH assumption: ",
  ifelse(
    ME9_global_ph_p < 0.05,
    "potential violation",
    "no evidence of violation"
  ),
  "\n",
  sep = ""
)


# ============================================================
# Step 14F-3F-4
# Sensitivity analysis for global PH violation
#
# Purpose:
#   Allow Gender effect to vary with time because:
#     ME2 model: Gender PH P = 0.021
#     ME9 model: Gender PH P = 0.029
#
# Primary module effects remain time-constant because:
#     ME2 PH P = 0.151
#     ME9 PH P = 0.307
#
# DO NOT modify existing conventional Cox models.
# ============================================================
cat("\n================ Step 14F-3F-4A ================\n")
cat("Prepare numeric Gender variable for time-varying sensitivity analysis\n")
cat("================================================\n\n")

# ------------------------------------------------------------
# 1. Check original Gender coding
# ------------------------------------------------------------

cat("Original Gender_trait:\n")
print(table(surv_ME2_ME9_485$Gender_trait, useNA = "ifany"))

cat("\nGender_trait class:\n")
print(class(surv_ME2_ME9_485$Gender_trait))


# ------------------------------------------------------------
# 2. Create numeric Gender_Male variable
#    Female = 0
#    Male   = 1
# ------------------------------------------------------------

surv_ME2_ME9_485$Gender_Male <- ifelse(
  surv_ME2_ME9_485$Gender_trait == "Male",
  1,
  0
)

cat("\nNumeric Gender_Male coding:\n")
print(table(
  surv_ME2_ME9_485$Gender_Male,
  useNA = "ifany"
))

cat("\nCross-tabulation:\n")
print(table(
  surv_ME2_ME9_485$Gender_trait,
  surv_ME2_ME9_485$Gender_Male,
  useNA = "ifany"
))


# ------------------------------------------------------------
# 3. Basic validation
# ------------------------------------------------------------

stopifnot(
  nrow(surv_ME2_ME9_485) == 485,
  sum(is.na(surv_ME2_ME9_485$Gender_Male)) == 0,
  all(
    surv_ME2_ME9_485$Gender_Male[
      surv_ME2_ME9_485$Gender_trait == "Female"
    ] == 0
  ),
  all(
    surv_ME2_ME9_485$Gender_Male[
      surv_ME2_ME9_485$Gender_trait == "Male"
    ] == 1
  )
)

cat("\nValidation passed.\n")
cat("N =", nrow(surv_ME2_ME9_485), "\n")
cat("Events =", sum(surv_ME2_ME9_485$OS_event), "\n")
cat("Gender_Male: Female=0, Male=1\n")


cat("\n================ Step 14F-3F-4B ================\n")
cat("Gender time-varying sensitivity models\n")
cat("================================================\n\n")


# ------------------------------------------------------------
# 1. ME2 + time-varying Gender
# ------------------------------------------------------------

ME2_gender_td_model <- survival::coxph(
  survival::Surv(
    OS_time,
    OS_event
  ) ~
    ME2_SD +
    Age_trait +
    Gender_Male +
    tt(Gender_Male) +
    Stage_trait,
  data = surv_ME2_ME9_485,
  ties = "efron",
  tt = function(x, t, ...) {
    x * (sqrt(t) / sqrt(1000))
  }
)


# ------------------------------------------------------------
# 2. ME9 + time-varying Gender
# ------------------------------------------------------------

ME9_gender_td_model <- survival::coxph(
  survival::Surv(
    OS_time,
    OS_event
  ) ~
    ME9_SD +
    Age_trait +
    Gender_Male +
    tt(Gender_Male) +
    Stage_trait,
  data = surv_ME2_ME9_485,
  ties = "efron",
  tt = function(x, t, ...) {
    x * (sqrt(t) / sqrt(1000))
  }
)


# ------------------------------------------------------------
# 3. Print model summaries
# ------------------------------------------------------------

cat("\n================ ME2 Gender-TD model ================\n")
print(summary(ME2_gender_td_model))


cat("\n================ ME9 Gender-TD model ================\n")
print(summary(ME9_gender_td_model))



cat("\n================ Step 14F-3F-4C ================\n")
cat("Conventional vs Gender-TD sensitivity comparison\n")
cat("================================================\n\n")


# ------------------------------------------------------------
# Helper function
# ------------------------------------------------------------

extract_cox_module <- function(model, term) {
  
  sm <- summary(model)
  
  coef_table <- sm$coefficients
  ci_table   <- sm$conf.int
  
  beta <- coef_table[term, "coef"]
  pval <- coef_table[term, "Pr(>|z|)"]
  
  hr <- ci_table[term, "exp(coef)"]
  lower <- ci_table[term, "lower .95"]
  upper <- ci_table[term, "upper .95"]
  
  data.frame(
    N = nrow(model.frame(model)),
    Events = model$nevent,
    Beta = beta,
    HR = hr,
    CI_lower = lower,
    CI_upper = upper,
    P = pval,
    C_index = unname(sm$concordance[1]),
    C_index_SE = unname(sm$concordance[2]),
    AIC = AIC(model),
    LRT_P = unname(sm$logtest["pvalue"]),
    row.names = NULL
  )
}


# ------------------------------------------------------------
# Extract four models
# ------------------------------------------------------------

ME2_conventional_summary <- extract_cox_module(
  ME2_cox_model,
  "ME2_SD"
)

ME2_gender_td_summary <- extract_cox_module(
  ME2_gender_td_model,
  "ME2_SD"
)

ME9_conventional_summary <- extract_cox_module(
  ME9_cox_model,
  "ME9_SD"
)

ME9_gender_td_summary <- extract_cox_module(
  ME9_gender_td_model,
  "ME9_SD"
)


# ------------------------------------------------------------
# Combine
# ------------------------------------------------------------

ME2_ME9_sensitivity_comparison <- rbind(
  cbind(
    Model = "ME2 conventional Cox",
    ME2_conventional_summary
  ),
  cbind(
    Model = "ME2 Gender-TD sensitivity",
    ME2_gender_td_summary
  ),
  cbind(
    Model = "ME9 conventional Cox",
    ME9_conventional_summary
  ),
  cbind(
    Model = "ME9 Gender-TD sensitivity",
    ME9_gender_td_summary
  )
)


# ------------------------------------------------------------
# Print
# ------------------------------------------------------------

print(ME2_ME9_sensitivity_comparison)


cat("\n================ Step 14F-4A ================\n")
cat("Prepare formal MM / GS analysis objects\n")
cat("================================================\n\n")


# ------------------------------------------------------------
# 1. Check required objects
# ------------------------------------------------------------

required_objects <- c(
  "datExpr_MAD",
  "moduleColors",
  "ME_color_formal_raw",
  "annotation",
  "expr_estimate",
  "module_tme_results",
  "surv_ME2_ME9_485"
)

cat("Checking required objects:\n")

for (obj in required_objects) {
  cat(
    obj,
    ":",
    exists(obj),
    "\n"
  )
}

stopifnot(
  exists("datExpr_MAD"),
  exists("moduleColors"),
  exists("ME_color_formal_raw"),
  exists("annotation")
)


# ------------------------------------------------------------
# 2. Validate WGCNA expression / module alignment
# ------------------------------------------------------------

cat("\nWGCNA dimensions:\n")
cat(
  "datExpr_MAD:",
  nrow(datExpr_MAD),
  "samples x",
  ncol(datExpr_MAD),
  "genes\n"
)

cat(
  "moduleColors:",
  length(moduleColors),
  "genes\n"
)

cat(
  "Formal MEs:",
  nrow(ME_color_formal_raw),
  "samples x",
  ncol(ME_color_formal_raw),
  "modules\n"
)


stopifnot(
  length(moduleColors) == ncol(datExpr_MAD),
  nrow(ME_color_formal_raw) == nrow(datExpr_MAD),
  identical(
    rownames(ME_color_formal_raw),
    rownames(datExpr_MAD)
  )
)


# ------------------------------------------------------------
# 3. Check formal module eigengenes
# ------------------------------------------------------------

required_MEs <- c(
  "MEblue",
  "MEmagenta",
  "MEbrown"
)

cat("\nRequired formal module eigengenes:\n")

print(
  required_MEs %in% colnames(ME_color_formal_raw)
)

stopifnot(
  all(required_MEs %in% colnames(ME_color_formal_raw))
)


# ------------------------------------------------------------
# 4. Check core genes
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

cat("\nCore genes:\n")
print(core_genes)


# ------------------------------------------------------------
# 5. Confirm gene annotation
# ------------------------------------------------------------

core_annotation <- annotation[
  annotation$gene_name %in% core_genes,
  c(
    "Ensembl_ID",
    "gene_name",
    "gene_type",
    "Annotation_status"
  ),
  drop = FALSE
]

core_annotation <- core_annotation[
  match(
    core_genes,
    core_annotation$gene_name
  ),
  ,
  drop = FALSE
]

cat("\nCore gene annotation:\n")
print(core_annotation)


stopifnot(
  all(
    core_genes %in%
      core_annotation$gene_name
  )
)


# ------------------------------------------------------------
# 6. Locate core genes in datExpr_MAD
# ------------------------------------------------------------

core_gene_ensembl <- core_annotation$Ensembl_ID

core_gene_location <- data.frame(
  Gene = core_genes,
  Ensembl_ID = core_gene_ensembl,
  ModuleColor = moduleColors[
    match(
      core_gene_ensembl,
      colnames(datExpr_MAD)
    )
  ],
  stringsAsFactors = FALSE
)

cat("\nCore gene module locations:\n")
print(core_gene_location)


# ------------------------------------------------------------
# 7. Validate expected module assignments
# ------------------------------------------------------------

expected_modules <- c(
  AOC3 = "magenta",
  MS4A7 = "blue",
  IL7R = "blue"
)

actual_modules <- setNames(
  core_gene_location$ModuleColor,
  core_gene_location$Gene
)

print(actual_modules)

stopifnot(
  all(
    actual_modules[names(expected_modules)] ==
      expected_modules
  )
)


cat("\nStep 14F-4A validation PASSED.\n")
cat("AOC3  -> magenta / ME9\n")
cat("MS4A7 -> blue / ME2\n")
cat("IL7R  -> blue / ME2\n")



cat("\n================ Step 14F-4B ================\n")
cat("Calculate formal Module Membership (MM / kME)\n")
cat("================================================\n\n")


# ------------------------------------------------------------
# 1. Calculate genome-wide MM / kME
# ------------------------------------------------------------

MM_all <- WGCNA::signedKME(
  datExpr_MAD,
  ME_color_formal_raw,
  outputColumnName = "MM"
)


# ------------------------------------------------------------
# 2. Basic dimensions
# ------------------------------------------------------------

cat("MM_all dimensions:\n")
print(dim(MM_all))

cat("\nExpected:\n")
cat(
  "Genes:",
  ncol(datExpr_MAD),
  "\n"
)

cat(
  "Modules:",
  ncol(ME_color_formal_raw),
  "\n"
)


stopifnot(
  nrow(MM_all) == ncol(datExpr_MAD),
  ncol(MM_all) == ncol(ME_color_formal_raw)
)


# ------------------------------------------------------------
# 3. Align MM rows with WGCNA gene order
# ------------------------------------------------------------

rownames(MM_all) <- colnames(datExpr_MAD)

stopifnot(
  identical(
    rownames(MM_all),
    colnames(datExpr_MAD)
  )
)


# ------------------------------------------------------------
# 4. Inspect actual MM column names
# ------------------------------------------------------------

cat("\nMM column names:\n")
print(colnames(MM_all))


# ------------------------------------------------------------
# 5. Check for expected module MM columns
# ------------------------------------------------------------

expected_MM_columns <- c(
  "MMblue",
  "MEmagenta",
  "MMbrown"
)

cat("\nChecking expected MM columns:\n")
print(
  expected_MM_columns %in%
    colnames(MM_all)
)


# ------------------------------------------------------------
# 6. Identify MM columns corresponding to formal MEs
# ------------------------------------------------------------

MM_module_map <- data.frame(
  ME = colnames(ME_color_formal_raw),
  MM_column = paste0(
    "MM",
    sub(
      "^ME",
      "",
      colnames(ME_color_formal_raw)
    )
  ),
  stringsAsFactors = FALSE
)

cat("\nFormal ME -> MM column mapping:\n")
print(MM_module_map)


# ------------------------------------------------------------
# 7. Confirm required columns exist
# ------------------------------------------------------------

required_MM_columns <- c(
  "MMblue",
  "MMmagenta",
  "MMbrown"
)

stopifnot(
  all(
    required_MM_columns %in%
      colnames(MM_all)
  )
)


# ------------------------------------------------------------
# 8. Extract MM for the three core genes
# ------------------------------------------------------------

MM_core <- data.frame(
  Gene = core_genes,
  
  Ensembl_ID = core_gene_ensembl,
  
  Module = c(
    "magenta",
    "blue",
    "blue"
  ),
  
  Target_ME = c(
    "MEmagenta",
    "MEblue",
    "MEblue"
  ),
  
  MM_target = c(
    MM_all[
      core_gene_ensembl[1],
      "MMmagenta"
    ],
    
    MM_all[
      core_gene_ensembl[2],
      "MMblue"
    ],
    
    MM_all[
      core_gene_ensembl[3],
      "MMblue"
    ]
  ),
  
  MM_brown = MM_all[
    core_gene_ensembl,
    "MMbrown"
  ],
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 9. Add absolute MM
# ------------------------------------------------------------

MM_core$Abs_MM_target <- abs(
  MM_core$MM_target
)


# ------------------------------------------------------------
# 10. Print core-gene MM results
# ------------------------------------------------------------

cat("\nCore-gene Module Membership:\n")
print(MM_core)


# ------------------------------------------------------------
# 11. Basic sanity checks
# ------------------------------------------------------------

stopifnot(
  all(
    is.finite(
      MM_core$MM_target
    )
  ),
  all(
    is.finite(
      MM_core$MM_brown
    )
  )
)


cat("\nStep 14F-4B completed successfully.\n")


# ============================================================
# Step 14F-4C
# Gene Significance (GS) for tumor microenvironment traits
# ============================================================

cat("\n================ R Session Diagnosis ================\n\n")

cat("Current working directory:\n")
print(getwd())

cat("\nNumber of objects in current environment:\n")
print(length(ls()))

cat("\nFirst objects in current environment:\n")
print(head(ls(), 30))

cat("\nCheck key Workflow-03 objects:\n")

check_objects <- c(
  "datExpr_MAD",
  "expr_estimate",
  "moduleColors",
  "ME_color_formal_raw",
  "MM_all",
  "MM_core",
  "module_tme_results",
  "surv_ME2_ME9_485",
  "core_genes",
  "core_gene_ensembl"
)

object_status <- data.frame(
  Object = check_objects,
  Exists = vapply(
    check_objects,
    exists,
    logical(1),
    inherits = FALSE
  ),
  stringsAsFactors = FALSE
)

print(object_status)

cat("\n=====================================================\n")


cat("\n============================================================\n")
cat("Workflow-03 object recovery diagnosis\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Search for WGCNA expression-related objects
# ------------------------------------------------------------

cat("### 1. Expression / WGCNA candidate objects ###\n\n")

all_object_names <- ls()

wgcna_candidates <- all_object_names[
  grepl(
    "datexpr|expr|expression|wgcna|mes|module",
    tolower(all_object_names)
  )
]

print(wgcna_candidates)

# ------------------------------------------------------------
# 2. Search for TME / ESTIMATE-related objects
# ------------------------------------------------------------

cat("\n### 2. TME / ESTIMATE candidate objects ###\n\n")

tme_candidates <- all_object_names[
  grepl(
    "estimate|immune|stromal|tme|score",
    tolower(all_object_names)
  )
]

print(tme_candidates)


# ------------------------------------------------------------
# 3. Search for MM / kME-related objects
# ------------------------------------------------------------

cat("\n### 3. MM / kME candidate objects ###\n\n")

mm_candidates <- all_object_names[
  grepl(
    "mm|kme|membership",
    tolower(all_object_names)
  )
]

print(mm_candidates)


# ------------------------------------------------------------
# 4. Search for module-color-related objects
# ------------------------------------------------------------

cat("\n### 4. Module-color candidate objects ###\n\n")

color_candidates <- all_object_names[
  grepl(
    "color|colour",
    tolower(all_object_names)
  )
]

print(color_candidates)


# ------------------------------------------------------------
# 5. Search for core-gene objects
# ------------------------------------------------------------

cat("\n### 5. Core-gene candidate objects ###\n\n")

core_candidates <- all_object_names[
  grepl(
    "core|aoc3|ms4a7|il7r",
    tolower(all_object_names)
  )
]

print(core_candidates)


# ------------------------------------------------------------
# 6. Search for survival / cohort objects
# ------------------------------------------------------------

cat("\n### 6. Survival / cohort candidate objects ###\n\n")

survival_candidates <- all_object_names[
  grepl(
    "surv|cox|clinical|brown|485|me2|me9",
    tolower(all_object_names)
  )
]

print(survival_candidates)

# ------------------------------------------------------------
# 7. Find objects with key dimensions
# ------------------------------------------------------------

cat("\n### 7. Objects with key WGCNA dimensions ###\n\n")

dimension_candidates <- list()

for (obj_name in all_object_names) {
  
  obj <- tryCatch(
    get(obj_name),
    error = function(e) NULL
  )
  
  if (is.null(obj)) {
    next
  }
  
  d <- dim(obj)
  
  if (is.null(d) || length(d) != 2) {
    next
  }
  
  if (
    any(
      d %in% c(
        31609,
        513,
        503,
        485,
        50
      )
    )
  ) {
    
    dimension_candidates[[obj_name]] <- d
    
  }
}

if (length(dimension_candidates) > 0) {
  
  for (nm in names(dimension_candidates)) {
    
    cat(
      nm,
      " : ",
      paste(
        dimension_candidates[[nm]],
        collapse = " x "
      ),
      "\n",
      sep = ""
    )
    
  }
  
} else {
  
  cat("No obvious dimension match found.\n")
  
}


# ============================================================
# Step 14F-4C — Pre-check
# Validate expression, TME scores, MM, and core-gene objects
# ============================================================

cat("============================================================\n")
cat("Step 14F-4C — Object validation before GS calculation\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Expression matrix
# ------------------------------------------------------------

cat("### 1. datExpr_MAD ###\n")

cat("Dimensions: ",
    paste(dim(datExpr_MAD), collapse = " x "),
    "\n")

cat("Samples: ",
    length(rownames(datExpr_MAD)),
    "\n")

cat("Genes: ",
    length(colnames(datExpr_MAD)),
    "\n")

cat("Unique sample IDs: ",
    length(unique(rownames(datExpr_MAD))),
    "\n")

cat("Unique gene IDs: ",
    length(unique(colnames(datExpr_MAD))),
    "\n")

cat("Any NA: ",
    anyNA(datExpr_MAD),
    "\n\n")


# ------------------------------------------------------------
# 2. TME score matrix
# ------------------------------------------------------------

cat("### 2. estimate_scores_wgcna ###\n")

cat("Dimensions: ",
    paste(dim(estimate_scores_wgcna), collapse = " x "),
    "\n")

cat("Row names:\n")
print(rownames(estimate_scores_wgcna))

cat("\nFirst 6 column names:\n")
print(head(colnames(estimate_scores_wgcna)))

cat("\nLast 6 column names:\n")
print(tail(colnames(estimate_scores_wgcna)))

cat("\nAny NA: ",
    anyNA(estimate_scores_wgcna),
    "\n\n")


# ------------------------------------------------------------
# 3. Sample ID alignment
# ------------------------------------------------------------

expr_sample_ids <- rownames(datExpr_MAD)
tme_sample_ids  <- colnames(estimate_scores_wgcna)

cat("### 3. Sample ID alignment ###\n")

cat(
  "Expression samples = TME samples: ",
  length(expr_sample_ids) == length(tme_sample_ids),
  "\n"
)

cat(
  "Same set: ",
  setequal(expr_sample_ids, tme_sample_ids),
  "\n"
)

cat(
  "Identical order: ",
  identical(expr_sample_ids, tme_sample_ids),
  "\n"
)

cat(
  "Different positions: ",
  sum(expr_sample_ids != tme_sample_ids),
  "\n\n"
)


# ------------------------------------------------------------
# 4. TME trait names
# ------------------------------------------------------------

cat("### 4. TME trait names ###\n")

required_traits <- c(
  "ImmuneScore",
  "StromalScore",
  "ESTIMATEScore"
)

cat(
  "Required traits present: ",
  all(required_traits %in% rownames(estimate_scores_wgcna)),
  "\n"
)

print(
  setdiff(
    required_traits,
    rownames(estimate_scores_wgcna)
  )
)

cat("\n")


# ------------------------------------------------------------
# 5. MM_all
# ------------------------------------------------------------

cat("### 5. MM_all ###\n")

cat("Dimensions: ",
    paste(dim(MM_all), collapse = " x "),
    "\n")

cat(
  "Same genes as datExpr_MAD: ",
  identical(
    rownames(MM_all),
    colnames(datExpr_MAD)
  ),
  "\n"
)

cat(
  "All finite: ",
  all(is.finite(as.matrix(MM_all))),
  "\n\n"
)


# ------------------------------------------------------------
# 6. Core genes
# ------------------------------------------------------------

cat("### 6. Core genes ###\n")

print(core_genes)

cat("\nCore gene Ensembl IDs:\n")
print(core_gene_ensembl)

cat("\n")

cat(
  "Core genes present in datExpr_MAD: ",
  all(core_gene_ensembl %in% colnames(datExpr_MAD)),
  "\n"
)

cat(
  "Core genes present in MM_all: ",
  all(core_gene_ensembl %in% rownames(MM_all)),
  "\n"
)

cat("\n============================================================\n")
cat("Pre-check completed.\n")
cat("============================================================\n")


# ============================================================
# Step 14F-4C
# Gene Significance (GS) for tumor microenvironment traits
# ============================================================

cat("============================================================\n")
cat("Step 14F-4C — Gene Significance analysis\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Prepare TME traits
# ------------------------------------------------------------

tme_scores <- data.frame(
  ImmuneScore = as.numeric(
    estimate_scores_wgcna["ImmuneScore", ]
  ),
  StromalScore = as.numeric(
    estimate_scores_wgcna["StromalScore", ]
  ),
  ESTIMATEScore = as.numeric(
    estimate_scores_wgcna["ESTIMATEScore", ]
  ),
  row.names = colnames(estimate_scores_wgcna),
  check.names = FALSE
)


# ------------------------------------------------------------
# 2. Final sample-order validation
# ------------------------------------------------------------

stopifnot(
  identical(
    rownames(datExpr_MAD),
    rownames(tme_scores)
  )
)

cat("Sample alignment confirmed: TRUE\n")
cat("N samples: ", nrow(tme_scores), "\n")
cat("N WGCNA genes: ", ncol(datExpr_MAD), "\n\n")


# ------------------------------------------------------------
# 3. Function for genome-wide Spearman GS
# ------------------------------------------------------------

calculate_gs <- function(expr_matrix, trait_vector) {
  
  n_genes <- ncol(expr_matrix)
  
  rho <- numeric(n_genes)
  p_value <- numeric(n_genes)
  
  for (i in seq_len(n_genes)) {
    
    test <- suppressWarnings(
      cor.test(
        expr_matrix[, i],
        trait_vector,
        method = "spearman",
        exact = FALSE
      )
    )
    
    rho[i] <- unname(test$estimate)
    p_value[i] <- test$p.value
  }
  
  fdr <- p.adjust(
    p_value,
    method = "BH"
  )
  
  data.frame(
    Gene = colnames(expr_matrix),
    rho = rho,
    P = p_value,
    FDR = fdr,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 4. Calculate GS for ImmuneScore
# ------------------------------------------------------------

cat("Calculating GS for ImmuneScore...\n")

GS_Immune <- calculate_gs(
  datExpr_MAD,
  tme_scores$ImmuneScore
)

cat("Completed.\n\n")


# ------------------------------------------------------------
# 5. Calculate GS for StromalScore
# ------------------------------------------------------------

cat("Calculating GS for StromalScore...\n")

GS_Stromal <- calculate_gs(
  datExpr_MAD,
  tme_scores$StromalScore
)

cat("Completed.\n\n")


# ------------------------------------------------------------
# 6. Calculate GS for ESTIMATEScore
# ------------------------------------------------------------

cat("Calculating GS for ESTIMATEScore...\n")

GS_ESTIMATE <- calculate_gs(
  datExpr_MAD,
  tme_scores$ESTIMATEScore
)

cat("Completed.\n\n")


# ------------------------------------------------------------
# 7. Rename columns for integration
# ------------------------------------------------------------

GS_Immune <- GS_Immune[
  ,
  c("Gene", "rho", "P", "FDR")
]

colnames(GS_Immune) <- c(
  "Gene",
  "GS_Immune",
  "P_Immune",
  "FDR_Immune"
)


GS_Stromal <- GS_Stromal[
  ,
  c("Gene", "rho", "P", "FDR")
]

colnames(GS_Stromal) <- c(
  "Gene",
  "GS_Stromal",
  "P_Stromal",
  "FDR_Stromal"
)


GS_ESTIMATE <- GS_ESTIMATE[
  ,
  c("Gene", "rho", "P", "FDR")
]

colnames(GS_ESTIMATE) <- c(
  "Gene",
  "GS_ESTIMATE",
  "P_ESTIMATE",
  "FDR_ESTIMATE"
)


# ------------------------------------------------------------
# 8. Verify gene-order consistency
# ------------------------------------------------------------

stopifnot(
  identical(
    GS_Immune$Gene,
    GS_Stromal$Gene
  )
)

stopifnot(
  identical(
    GS_Immune$Gene,
    GS_ESTIMATE$Gene
  )
)

stopifnot(
  identical(
    GS_Immune$Gene,
    colnames(datExpr_MAD)
  )
)

cat("Gene-order consistency confirmed: TRUE\n\n")


# ------------------------------------------------------------
# 9. Integrate three GS results
# ------------------------------------------------------------

GS_all <- data.frame(
  Gene = GS_Immune$Gene,
  GS_Immune = GS_Immune$GS_Immune,
  P_Immune = GS_Immune$P_Immune,
  FDR_Immune = GS_Immune$FDR_Immune,
  
  GS_Stromal = GS_Stromal$GS_Stromal,
  P_Stromal = GS_Stromal$P_Stromal,
  FDR_Stromal = GS_Stromal$FDR_Stromal,
  
  GS_ESTIMATE = GS_ESTIMATE$GS_ESTIMATE,
  P_ESTIMATE = GS_ESTIMATE$P_ESTIMATE,
  FDR_ESTIMATE = GS_ESTIMATE$FDR_ESTIMATE,
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 10. Add absolute GS values
# ------------------------------------------------------------

GS_all$Abs_GS_Immune <-
  abs(GS_all$GS_Immune)

GS_all$Abs_GS_Stromal <-
  abs(GS_all$GS_Stromal)

GS_all$Abs_GS_ESTIMATE <-
  abs(GS_all$GS_ESTIMATE)


# ------------------------------------------------------------
# 11. Extract core genes
# ------------------------------------------------------------

GS_core <- GS_all[
  GS_all$Gene %in% core_gene_ensembl,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 12. Add gene symbols
# ------------------------------------------------------------

core_symbol_map <- data.frame(
  Gene = core_gene_ensembl,
  Symbol = core_genes,
  stringsAsFactors = FALSE
)

GS_core <- merge(
  core_symbol_map,
  GS_core,
  by = "Gene",
  all.x = TRUE,
  sort = FALSE
)


# Restore the desired core-gene order
GS_core$Symbol <- factor(
  GS_core$Symbol,
  levels = core_genes
)

GS_core <- GS_core[
  order(GS_core$Symbol),
  ,
  drop = FALSE
]

GS_core$Symbol <- as.character(GS_core$Symbol)


# ------------------------------------------------------------
# 13. Compact core-gene table
# ------------------------------------------------------------

GS_core_compact <- GS_core[
  ,
  c(
    "Symbol",
    "Gene",
    "GS_Immune",
    "FDR_Immune",
    "GS_Stromal",
    "FDR_Stromal",
    "GS_ESTIMATE",
    "FDR_ESTIMATE"
  )
]


# ------------------------------------------------------------
# 14. Basic validation
# ------------------------------------------------------------

cat("============================================================\n")
cat("GS calculation completed.\n")
cat("============================================================\n\n")

cat("GS_all dimensions: ",
    paste(dim(GS_all), collapse = " x "),
    "\n")

cat(
  "Expected dimensions: 31609 x 13\n"
)

cat(
  "Core genes recovered: ",
  nrow(GS_core_compact),
  "\n"
)

cat("\nCore-gene GS results:\n")
print(GS_core_compact)

cat("\n")


# ------------------------------------------------------------
# 15. Genome-wide FDR summary
# ------------------------------------------------------------

cat("### Genome-wide FDR summary ###\n")

cat(
  "ImmuneScore FDR < 0.05: ",
  sum(GS_all$FDR_Immune < 0.05),
  "\n"
)

cat(
  "StromalScore FDR < 0.05: ",
  sum(GS_all$FDR_Stromal < 0.05),
  "\n"
)

cat(
  "ESTIMATEScore FDR < 0.05: ",
  sum(GS_all$FDR_ESTIMATE < 0.05),
  "\n"
)

cat("\n============================================================\n")
cat("Step 14F-4C finished.\n")
cat("============================================================\n")


# ============================================================
# Step 14F-4D-1 — CORRECTED
# Module-internal MM ranking for core genes
# ============================================================

cat("============================================================\n")
cat("Step 14F-4D-1 — Corrected module-internal MM ranking\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Validate core WGCNA objects
# ------------------------------------------------------------

stopifnot(
  length(moduleColors) == ncol(datExpr_MAD)
)

stopifnot(
  identical(
    names(moduleColors),
    colnames(datExpr_MAD)
  ) ||
    identical(
      as.character(moduleColors),
      as.character(moduleColors)
    )
)

stopifnot(
  identical(
    rownames(MM_all),
    colnames(datExpr_MAD)
  )
)

stopifnot(
  all(c("MMblue", "MMmagenta") %in% colnames(MM_all))
)

cat("WGCNA object validation: PASSED\n\n")


# ------------------------------------------------------------
# 2. Build gene-to-module table
# ------------------------------------------------------------

gene_module_df <- data.frame(
  Gene = colnames(datExpr_MAD),
  Module = as.character(moduleColors),
  stringsAsFactors = FALSE
)

cat("Total WGCNA genes: ",
    nrow(gene_module_df),
    "\n")

cat(
  "Unique modules: ",
  length(unique(gene_module_df$Module)),
  "\n\n"
)


# ------------------------------------------------------------
# 3. Confirm target module sizes
# ------------------------------------------------------------

blue_genes <- gene_module_df$Gene[
  gene_module_df$Module == "blue"
]

magenta_genes <- gene_module_df$Gene[
  gene_module_df$Module == "magenta"
]

cat("Blue module size: ",
    length(blue_genes),
    "\n")

cat("Magenta module size: ",
    length(magenta_genes),
    "\n\n")


# ------------------------------------------------------------
# 4. Function for TRUE module-internal MM ranking
# ------------------------------------------------------------

rank_module_MM_corrected <- function(
    MM_matrix,
    module_genes,
    module_column,
    module_name) {
  
  # Check that all module genes exist in MM matrix
  stopifnot(
    all(module_genes %in% rownames(MM_matrix))
  )
  
  df <- data.frame(
    Gene = module_genes,
    MM = MM_matrix[
      module_genes,
      module_column
    ],
    stringsAsFactors = FALSE
  )
  
  # Keep finite values
  df <- df[
    is.finite(df$MM),
    ,
    drop = FALSE
  ]
  
  # Absolute MM is used ONLY for ranking
  df$Abs_MM <- abs(df$MM)
  
  # Rank within this module only
  df <- df[
    order(
      -df$Abs_MM,
      -df$MM
    ),
    ,
    drop = FALSE
  ]
  
  df$MM_rank <- seq_len(nrow(df))
  
  # Percentile within the module
  if (nrow(df) > 1) {
    
    df$MM_percentile <- (
      1 -
        (df$MM_rank - 1) /
        (nrow(df) - 1)
    ) * 100
    
  } else {
    
    df$MM_percentile <- 100
    
  }
  
  df$Module <- module_name
  
  df <- df[
    ,
    c(
      "Gene",
      "Module",
      "MM",
      "Abs_MM",
      "MM_rank",
      "MM_percentile"
    )
  ]
  
  return(df)
}


# ------------------------------------------------------------
# 5. Blue / ME2
# ------------------------------------------------------------

MM_ME2_rank_corrected <- rank_module_MM_corrected(
  MM_matrix = MM_all,
  module_genes = blue_genes,
  module_column = "MMblue",
  module_name = "blue / ME2"
)


# ------------------------------------------------------------
# 6. Magenta / ME9
# ------------------------------------------------------------

MM_ME9_rank_corrected <- rank_module_MM_corrected(
  MM_matrix = MM_all,
  module_genes = magenta_genes,
  module_column = "MMmagenta",
  module_name = "magenta / ME9"
)


# ------------------------------------------------------------
# 7. Verify that only true module members were ranked
# ------------------------------------------------------------

stopifnot(
  all(
    MM_ME2_rank_corrected$Gene %in% blue_genes
  )
)

stopifnot(
  all(
    MM_ME9_rank_corrected$Gene %in% magenta_genes
  )
)

cat("Module membership filtering: PASSED\n\n")


# ------------------------------------------------------------
# 8. Extract core genes from correct module rankings
# ------------------------------------------------------------

core_symbol_map <- data.frame(
  Gene = core_gene_ensembl,
  Symbol = core_genes,
  stringsAsFactors = FALSE
)


MM_core_ME2 <- merge(
  core_symbol_map,
  MM_ME2_rank_corrected,
  by = "Gene",
  all.x = FALSE,
  sort = FALSE
)


MM_core_ME9 <- merge(
  core_symbol_map,
  MM_ME9_rank_corrected,
  by = "Gene",
  all.x = FALSE,
  sort = FALSE
)


# ------------------------------------------------------------
# 9. Combine core-gene rankings
# ------------------------------------------------------------

MM_core_rank_corrected <- rbind(
  MM_core_ME9,
  MM_core_ME2
)


# ------------------------------------------------------------
# 10. Keep only biologically correct module assignments
# ------------------------------------------------------------

MM_core_rank_corrected <- MM_core_rank_corrected[
  (
    MM_core_rank_corrected$Symbol == "AOC3" &
      MM_core_rank_corrected$Module == "magenta / ME9"
  ) |
    (
      MM_core_rank_corrected$Symbol %in% c("MS4A7", "IL7R") &
        MM_core_rank_corrected$Module == "blue / ME2"
    ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 11. Order core genes
# ------------------------------------------------------------

MM_core_rank_corrected$Symbol <- factor(
  MM_core_rank_corrected$Symbol,
  levels = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )
)

MM_core_rank_corrected <- MM_core_rank_corrected[
  order(MM_core_rank_corrected$Symbol),
  ,
  drop = FALSE
]

MM_core_rank_corrected$Symbol <- as.character(
  MM_core_rank_corrected$Symbol
)


# ------------------------------------------------------------
# 12. Print module sizes
# ------------------------------------------------------------

cat("============================================================\n")
cat("TRUE module sizes\n")
cat("============================================================\n\n")

cat(
  "Blue / ME2: ",
  nrow(MM_ME2_rank_corrected),
  " genes\n"
)

cat(
  "Magenta / ME9: ",
  nrow(MM_ME9_rank_corrected),
  " genes\n\n"
)


# ------------------------------------------------------------
# 13. Print core-gene module-internal ranking
# ------------------------------------------------------------

cat("============================================================\n")
cat("Core-gene TRUE module-internal MM ranking\n")
cat("============================================================\n\n")

print(
  MM_core_rank_corrected[
    ,
    c(
      "Symbol",
      "Gene",
      "Module",
      "MM",
      "Abs_MM",
      "MM_rank",
      "MM_percentile"
    )
  ]
)


# ------------------------------------------------------------
# 14. Top 20 genes — blue / ME2
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Blue / ME2 — Top 20 TRUE module members by |MM|\n")
cat("============================================================\n\n")

print(
  head(
    MM_ME2_rank_corrected,
    20
  )
)


# ------------------------------------------------------------
# 15. Top 20 genes — magenta / ME9
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Magenta / ME9 — Top 20 TRUE module members by |MM|\n")
cat("============================================================\n\n")

print(
  head(
    MM_ME9_rank_corrected,
    20
  )
)



# ============================================================
# Diagnostic — GS_all vs datExpr_MAD gene ID alignment
# ============================================================

cat("============================================================\n")
cat("Diagnostic: GS_all vs datExpr_MAD alignment\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Dimensions
# ------------------------------------------------------------

cat("datExpr_MAD dimensions:\n")
print(dim(datExpr_MAD))

cat("\nGS_all dimensions:\n")
print(dim(GS_all))


# ------------------------------------------------------------
# 2. Basic ID information
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("ID counts\n")
cat("------------------------------------------------------------\n\n")

cat(
  "datExpr_MAD gene IDs: ",
  ncol(datExpr_MAD),
  "\n"
)

cat(
  "GS_all gene IDs: ",
  nrow(GS_all),
  "\n"
)

cat(
  "Unique datExpr_MAD IDs: ",
  length(unique(colnames(datExpr_MAD))),
  "\n"
)

cat(
  "Unique GS_all IDs: ",
  length(unique(rownames(GS_all))),
  "\n"
)


# ------------------------------------------------------------
# 3. Exact identity
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Exact identity\n")
cat("------------------------------------------------------------\n\n")

cat(
  "identical IDs/order: ",
  identical(
    rownames(GS_all),
    colnames(datExpr_MAD)
  ),
  "\n"
)


# ------------------------------------------------------------
# 4. Set equality
# ------------------------------------------------------------

cat(
  "setequal IDs: ",
  setequal(
    rownames(GS_all),
    colnames(datExpr_MAD)
  ),
  "\n"
)


# ------------------------------------------------------------
# 5. IDs in datExpr but absent from GS_all
# ------------------------------------------------------------

missing_from_GS <- setdiff(
  colnames(datExpr_MAD),
  rownames(GS_all)
)

cat(
  "\nGenes in datExpr_MAD but absent from GS_all: ",
  length(missing_from_GS),
  "\n"
)

if (length(missing_from_GS) > 0) {
  
  cat("\nFirst 20 missing IDs:\n")
  print(
    head(
      missing_from_GS,
      20
    )
  )
}


# ------------------------------------------------------------
# 6. IDs in GS_all but absent from datExpr
# ------------------------------------------------------------

extra_in_GS <- setdiff(
  rownames(GS_all),
  colnames(datExpr_MAD)
)

cat(
  "\nGenes in GS_all but absent from datExpr_MAD: ",
  length(extra_in_GS),
  "\n"
)

if (length(extra_in_GS) > 0) {
  
  cat("\nFirst 20 extra GS IDs:\n")
  print(
    head(
      extra_in_GS,
      20
    )
  )
}


# ------------------------------------------------------------
# 7. If sets are equal, check ordering
# ------------------------------------------------------------

if (
  setequal(
    rownames(GS_all),
    colnames(datExpr_MAD)
  )
) {
  
  same_position <- (
    rownames(GS_all) ==
      colnames(datExpr_MAD)
  )
  
  cat(
    "\nSame-position IDs: ",
    sum(same_position),
    "/",
    length(same_position),
    "\n"
  )
  
  cat(
    "Different-position IDs: ",
    sum(!same_position),
    "\n"
  )
  
  if (any(!same_position)) {
    
    cat("\nFirst 20 position mismatches:\n")
    
    mismatch_idx <- which(
      !same_position
    )
    
    print(
      data.frame(
        Position = head(
          mismatch_idx,
          20
        ),
        datExpr_ID = colnames(
          datExpr_MAD
        )[
          head(
            mismatch_idx,
            20
          )
        ],
        GS_all_ID = rownames(
          GS_all
        )[
          head(
            mismatch_idx,
            20
          )
        ],
        stringsAsFactors = FALSE
      )
    )
  }
}


# ------------------------------------------------------------
# 8. Check core genes specifically
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Core-gene ID presence\n")
cat("------------------------------------------------------------\n\n")

core_check <- data.frame(
  Gene = core_gene_ensembl,
  In_datExpr = core_gene_ensembl %in%
    colnames(datExpr_MAD),
  In_GS_all = core_gene_ensembl %in%
    rownames(GS_all),
  stringsAsFactors = FALSE
)

print(
  core_check,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("Diagnostic completed\n")
cat("============================================================\n")


#### ============================================================
# Diagnostic — inspect GS_all column structure
# ============================================================

cat("============================================================\n")
cat("Inspect GS_all column structure\n")
cat("============================================================\n\n")

cat("Column names of GS_all:\n")
print(
  colnames(GS_all)
)

cat("\nFirst 6 rows of GS_all:\n")
print(
  head(
    GS_all
  )
)

cat("\nStructure of GS_all:\n")
str(GS_all)

cat("\n============================================================\n")
cat("Diagnostic completed\n")
cat("============================================================\n")


# ============================================================
# Step 14F-4D-2A
# Restore and validate GS_all gene identifiers
# IMPORTANT: do NOT recalculate GS
# ============================================================

cat("============================================================\n")
cat("Step 14F-4D-2A — Restore and validate GS_all IDs\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Basic validation
# ------------------------------------------------------------

stopifnot(
  "Gene" %in% colnames(GS_all)
)

stopifnot(
  nrow(GS_all) == ncol(datExpr_MAD)
)

stopifnot(
  !anyNA(GS_all$Gene)
)

stopifnot(
  !anyDuplicated(GS_all$Gene)
)


# ------------------------------------------------------------
# 2. Restore rownames from the original Gene column
# ------------------------------------------------------------

GS_all$Gene <- as.character(GS_all$Gene)

rownames(GS_all) <- GS_all$Gene


# ------------------------------------------------------------
# 3. Validate exact gene-ID matching
# ------------------------------------------------------------

cat("Gene ID validation:\n\n")

cat(
  "Number of genes in datExpr_MAD:",
  ncol(datExpr_MAD),
  "\n"
)

cat(
  "Number of genes in GS_all:",
  nrow(GS_all),
  "\n"
)

cat(
  "Unique GS_all Gene IDs:",
  length(unique(GS_all$Gene)),
  "\n"
)

cat(
  "Exact rowname/column-name identity:",
  identical(
    rownames(GS_all),
    colnames(datExpr_MAD)
  ),
  "\n"
)

cat(
  "Same gene set:",
  setequal(
    rownames(GS_all),
    colnames(datExpr_MAD)
  ),
  "\n"
)


# ------------------------------------------------------------
# 4. Check order explicitly
# ------------------------------------------------------------

same_order <- all(
  rownames(GS_all) == colnames(datExpr_MAD)
)

cat(
  "Same gene order:",
  same_order,
  "\n"
)


# ------------------------------------------------------------
# 5. Check missing / extra genes
# ------------------------------------------------------------

missing_from_GS <- setdiff(
  colnames(datExpr_MAD),
  rownames(GS_all)
)

extra_in_GS <- setdiff(
  rownames(GS_all),
  colnames(datExpr_MAD)
)

cat(
  "Genes missing from GS_all:",
  length(missing_from_GS),
  "\n"
)

cat(
  "Extra genes in GS_all:",
  length(extra_in_GS),
  "\n"
)


# ------------------------------------------------------------
# 6. Core-gene validation
# ------------------------------------------------------------

core_genes <- c(
  "ENSG00000131471",  # AOC3
  "ENSG00000166927",  # MS4A7
  "ENSG00000168685"   # IL7R
)

core_symbols <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_check <- data.frame(
  Symbol = core_symbols,
  Ensembl_ID = core_genes,
  In_datExpr = core_genes %in% colnames(datExpr_MAD),
  In_GS_all = core_genes %in% rownames(GS_all),
  stringsAsFactors = FALSE
)

cat("\nCore-gene ID validation:\n")
print(core_check)


# ------------------------------------------------------------
# 7. Validate core-gene GS values
# ------------------------------------------------------------

GS_core_validation <- GS_all[
  core_genes,
  c(
    "Gene",
    "GS_Immune",
    "FDR_Immune",
    "GS_Stromal",
    "FDR_Stromal",
    "GS_ESTIMATE",
    "FDR_ESTIMATE"
  ),
  drop = FALSE
]

GS_core_validation$Symbol <- core_symbols[
  match(
    rownames(GS_core_validation),
    core_genes
  )
]

GS_core_validation <- GS_core_validation[
  ,
  c(
    "Symbol",
    "Gene",
    "GS_Immune",
    "FDR_Immune",
    "GS_Stromal",
    "FDR_Stromal",
    "GS_ESTIMATE",
    "FDR_ESTIMATE"
  )
]

cat("\nCore-gene GS validation:\n")
print(
  GS_core_validation,
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Final structural validation
# ------------------------------------------------------------

stopifnot(
  identical(
    rownames(GS_all),
    colnames(datExpr_MAD)
  )
)

stopifnot(
  identical(
    rownames(GS_all),
    names(moduleColors)
  ) ||
    identical(
      rownames(GS_all),
      colnames(datExpr_MAD)
    )
)

stopifnot(
  all(core_genes %in% rownames(GS_all))
)


cat("\n============================================================\n")
cat("GS_all ID restoration and validation completed\n")
cat("GS values were NOT recalculated\n")
cat("============================================================\n")


# ============================================================
# Step 14F-4D-2B
# Construct module-specific MM–GS datasets
# ============================================================

cat("============================================================\n")
cat("Step 14F-4D-2B — Construct MM–GS datasets\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Build integrated MM–GS table
# ------------------------------------------------------------

MM_GS_all <- data.frame(
  Gene = colnames(datExpr_MAD),
  Module = as.character(moduleColors),
  
  MM_blue = MM_all[
    colnames(datExpr_MAD),
    "MMblue"
  ],
  
  MM_magenta = MM_all[
    colnames(datExpr_MAD),
    "MMmagenta"
  ],
  
  GS_Immune = GS_all[
    colnames(datExpr_MAD),
    "GS_Immune"
  ],
  
  GS_Stromal = GS_all[
    colnames(datExpr_MAD),
    "GS_Stromal"
  ],
  
  GS_ESTIMATE = GS_all[
    colnames(datExpr_MAD),
    "GS_ESTIMATE"
  ],
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 2. Restore gene IDs as rownames
# ------------------------------------------------------------

rownames(MM_GS_all) <- MM_GS_all$Gene


# ------------------------------------------------------------
# 3. Structural validation
# ------------------------------------------------------------

stopifnot(
  nrow(MM_GS_all) == ncol(datExpr_MAD)
)

stopifnot(
  identical(
    rownames(MM_GS_all),
    colnames(datExpr_MAD)
  )
)

stopifnot(
  !anyNA(MM_GS_all)
)

stopifnot(
  !anyDuplicated(MM_GS_all$Gene)
)


# ------------------------------------------------------------
# 4. Extract actual module members
# ------------------------------------------------------------

ME2_df <- MM_GS_all[
  MM_GS_all$Module == "blue",
  ,
  drop = FALSE
]

ME9_df <- MM_GS_all[
  MM_GS_all$Module == "magenta",
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 5. Validate module sizes
# ------------------------------------------------------------

cat("Module sizes:\n\n")

cat(
  "Blue / ME2:",
  nrow(ME2_df),
  "genes\n"
)

cat(
  "Magenta / ME9:",
  nrow(ME9_df),
  "genes\n"
)


# ------------------------------------------------------------
# 6. Validate MM values
# ------------------------------------------------------------

cat("\nMM value validation:\n\n")

cat(
  "Blue MM finite:",
  all(is.finite(ME2_df$MM_blue)),
  "\n"
)

cat(
  "Magenta MM finite:",
  all(is.finite(ME9_df$MM_magenta)),
  "\n"
)


# ------------------------------------------------------------
# 7. Validate core-gene module membership
# ------------------------------------------------------------

core_module_check <- MM_GS_all[
  core_genes,
  c(
    "Gene",
    "Module",
    "MM_blue",
    "MM_magenta",
    "GS_Immune",
    "GS_Stromal",
    "GS_ESTIMATE"
  ),
  drop = FALSE
]

core_module_check$Symbol <- core_symbols[
  match(
    rownames(core_module_check),
    core_genes
  )
]

core_module_check <- core_module_check[
  ,
  c(
    "Symbol",
    "Gene",
    "Module",
    "MM_blue",
    "MM_magenta",
    "GS_Immune",
    "GS_Stromal",
    "GS_ESTIMATE"
  )
]

cat("\nCore-gene MM–GS values:\n")

print(
  core_module_check,
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Final confirmation
# ------------------------------------------------------------

stopifnot(
  nrow(ME2_df) == 1798
)

stopifnot(
  nrow(ME9_df) == 534
)

stopifnot(
  ME2_df["ENSG00000166927", "Module"] == "blue"
)

stopifnot(
  ME2_df["ENSG00000168685", "Module"] == "blue"
)

stopifnot(
  ME9_df["ENSG00000131471", "Module"] == "magenta"
)


cat("\n============================================================\n")
cat("MM–GS datasets successfully constructed\n")
cat("Blue / ME2 = 1,798 genes\n")
cat("Magenta / ME9 = 534 genes\n")
cat("No MM or GS values were recalculated\n")
cat("============================================================\n")

### ============================================================
# Step 14F-4D-2C
# Module-wide MM–GS Spearman correlations
# ============================================================

cat("============================================================\n")
cat("Step 14F-4D-2C — Module-wide MM–GS correlations\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Define correlation function
# ------------------------------------------------------------

calc_mm_gs_cor <- function(
    df,
    mm_col,
    gs_col,
    module_name,
    trait_name
) {
  
  x <- df[[mm_col]]
  y <- df[[gs_col]]
  
  test <- suppressWarnings(
    cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )
  
  data.frame(
    Module = module_name,
    Trait = trait_name,
    N = length(x),
    MM_variable = mm_col,
    GS_variable = gs_col,
    Spearman_rho = unname(test$estimate),
    P_value = test$p.value,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 2. Six predefined MM–GS tests
# ------------------------------------------------------------

mm_gs_results_raw <- rbind(
  
  # ----------------------------------------------------------
  # Blue / ME2
  # ----------------------------------------------------------
  
  calc_mm_gs_cor(
    df = ME2_df,
    mm_col = "MM_blue",
    gs_col = "GS_Immune",
    module_name = "blue / ME2",
    trait_name = "ImmuneScore"
  ),
  
  calc_mm_gs_cor(
    df = ME2_df,
    mm_col = "MM_blue",
    gs_col = "GS_Stromal",
    module_name = "blue / ME2",
    trait_name = "StromalScore"
  ),
  
  calc_mm_gs_cor(
    df = ME2_df,
    mm_col = "MM_blue",
    gs_col = "GS_ESTIMATE",
    module_name = "blue / ME2",
    trait_name = "ESTIMATEScore"
  ),
  
  # ----------------------------------------------------------
  # Magenta / ME9
  # ----------------------------------------------------------
  
  calc_mm_gs_cor(
    df = ME9_df,
    mm_col = "MM_magenta",
    gs_col = "GS_Immune",
    module_name = "magenta / ME9",
    trait_name = "ImmuneScore"
  ),
  
  calc_mm_gs_cor(
    df = ME9_df,
    mm_col = "MM_magenta",
    gs_col = "GS_Stromal",
    module_name = "magenta / ME9",
    trait_name = "StromalScore"
  ),
  
  calc_mm_gs_cor(
    df = ME9_df,
    mm_col = "MM_magenta",
    gs_col = "GS_ESTIMATE",
    module_name = "magenta / ME9",
    trait_name = "ESTIMATEScore"
  )
  
)


# ------------------------------------------------------------
# 3. Unified BH/FDR correction across all six tests
# ------------------------------------------------------------

mm_gs_results <- mm_gs_results_raw

mm_gs_results$FDR <- p.adjust(
  mm_gs_results$P_value,
  method = "BH"
)


# ------------------------------------------------------------
# 4. Add absolute correlation and direction
# ------------------------------------------------------------

mm_gs_results$Abs_rho <- abs(
  mm_gs_results$Spearman_rho
)

mm_gs_results$Direction <- ifelse(
  mm_gs_results$Spearman_rho > 0,
  "Positive",
  "Negative"
)

mm_gs_results$FDR_significant <- (
  mm_gs_results$FDR < 0.05
)


# ------------------------------------------------------------
# 5. Arrange columns
# ------------------------------------------------------------

mm_gs_results <- mm_gs_results[
  ,
  c(
    "Module",
    "Trait",
    "N",
    "MM_variable",
    "GS_variable",
    "Spearman_rho",
    "Abs_rho",
    "P_value",
    "FDR",
    "Direction",
    "FDR_significant"
  )
]


# ------------------------------------------------------------
# 6. Print complete results
# ------------------------------------------------------------

cat("Complete MM–GS correlation results:\n\n")

print(
  mm_gs_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Compact numerical summary
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Numerical summary\n")
cat("------------------------------------------------------------\n\n")

for (i in seq_len(nrow(mm_gs_results))) {
  
  cat(
    mm_gs_results$Module[i],
    "|",
    mm_gs_results$Trait[i],
    "\n"
  )
  
  cat(
    "  rho =",
    sprintf(
      "%.6f",
      mm_gs_results$Spearman_rho[i]
    ),
    "\n"
  )
  
  cat(
    "  P =",
    format.pval(
      mm_gs_results$P_value[i],
      digits = 6,
      eps = 1e-300
    ),
    "\n"
  )
  
  cat(
    "  FDR =",
    format.pval(
      mm_gs_results$FDR[i],
      digits = 6,
      eps = 1e-300
    ),
    "\n"
  )
  
  cat(
    "  Direction =",
    mm_gs_results$Direction[i],
    "\n\n"
  )
}


# ------------------------------------------------------------
# 8. Final validation
# ------------------------------------------------------------

stopifnot(
  nrow(mm_gs_results) == 6
)

stopifnot(
  all(
    mm_gs_results$N ==
      c(
        nrow(ME2_df),
        nrow(ME2_df),
        nrow(ME2_df),
        nrow(ME9_df),
        nrow(ME9_df),
        nrow(ME9_df)
      )
  )
)

stopifnot(
  all(is.finite(mm_gs_results$Spearman_rho))
)

stopifnot(
  all(is.finite(mm_gs_results$P_value))
)

stopifnot(
  all(is.finite(mm_gs_results$FDR))
)

stopifnot(
  all(
    mm_gs_results$FDR ==
      p.adjust(
        mm_gs_results$P_value,
        method = "BH"
      )
  )
)


cat("============================================================\n")
cat("Step 14F-4D-2C completed\n")
cat("Six MM–GS tests with unified BH/FDR correction\n")
cat("============================================================\n")


# ============================================================
# Step 14F-4D-2D
# Core-gene integrated positioning
#
# Purpose:
# Integrate:
#   1. module localization
#   2. module size
#   3. module membership (MM)
#   4. MM rank / percentile
#   5. gene significance (GS)
#   6. corresponding module-level TME correlations
#   7. corresponding module-wide MM–GS correlations
#
# IMPORTANT:
# No new GS/MM correlations are calculated here.
# All statistics are retrieved from LOCKED objects/results.
# ============================================================

cat("============================================================\n")
cat("Step 14F-4D-2D — Core-gene integrated positioning\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Core-gene definition
# ------------------------------------------------------------

core_info <- data.frame(
  Symbol = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  
  Gene = c(
    "ENSG00000131471",
    "ENSG00000166927",
    "ENSG00000168685"
  ),
  
  Module = c(
    "magenta",
    "blue",
    "blue"
  ),
  
  ME = c(
    "ME9",
    "ME2",
    "ME2"
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 2. Validate core-gene IDs
# ------------------------------------------------------------

stopifnot(
  all(
    core_info$Gene %in%
      colnames(datExpr_MAD)
  )
)

stopifnot(
  all(
    core_info$Gene %in%
      rownames(MM_all)
  )
)

stopifnot(
  all(
    core_info$Gene %in%
      rownames(GS_all)
  )
)


# ------------------------------------------------------------
# 3. Define module-level TME correlations
#    from LOCKED Step 14F-2 results
# ------------------------------------------------------------

module_tme_lookup <- data.frame(
  
  Module = c(
    "blue",
    "blue",
    "blue",
    "magenta",
    "magenta",
    "magenta"
  ),
  
  ME = c(
    "ME2",
    "ME2",
    "ME2",
    "ME9",
    "ME9",
    "ME9"
  ),
  
  Trait = c(
    "ImmuneScore",
    "StromalScore",
    "ESTIMATEScore",
    "ImmuneScore",
    "StromalScore",
    "ESTIMATEScore"
  ),
  
  Module_TME_rho = c(
    0.9689631,
    0.7081174,
    0.9179035,
    0.4209524,
    0.6451482,
    0.5695171
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 4. Define module-wide MM–GS correlations
#    from LOCKED Step 14F-4D-2C
# ------------------------------------------------------------

mm_gs_lookup <- data.frame(
  
  Module = c(
    "blue",
    "blue",
    "blue",
    "magenta",
    "magenta",
    "magenta"
  ),
  
  ME = c(
    "ME2",
    "ME2",
    "ME2",
    "ME9",
    "ME9",
    "ME9"
  ),
  
  Trait = c(
    "ImmuneScore",
    "StromalScore",
    "ESTIMATEScore",
    "ImmuneScore",
    "StromalScore",
    "ESTIMATEScore"
  ),
  
  MM_GS_rho = c(
    0.9433047,
    0.7068876,
    0.9081796,
    0.4023558,
    0.5221497,
    0.5088140
  ),
  
  MM_GS_FDR = c(
    0,
    4.394097e-272,
    0,
    3.373681e-22,
    1.660438e-38,
    1.984393e-36
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 5. Calculate module size
# ------------------------------------------------------------

module_size_lookup <- data.frame(
  
  Module = c(
    "blue",
    "magenta"
  ),
  
  ME = c(
    "ME2",
    "ME9"
  ),
  
  Module_Size = c(
    sum(moduleColors == "blue"),
    sum(moduleColors == "magenta")
  ),
  
  stringsAsFactors = FALSE
)


# Validate known sizes
stopifnot(
  module_size_lookup[
    module_size_lookup$Module == "blue",
    "Module_Size"
  ] == 1798
)

stopifnot(
  module_size_lookup[
    module_size_lookup$Module == "magenta",
    "Module_Size"
  ] == 534
)


# ------------------------------------------------------------
# 6. Calculate core-gene MM
# ------------------------------------------------------------

core_position <- core_info

core_position$MM <- NA_real_

core_position$MM[core_position$Module == "blue"] <-
  MM_all[
    core_position$Gene[core_position$Module == "blue"],
    "MMblue"
  ]

core_position$MM[core_position$Module == "magenta"] <-
  MM_all[
    core_position$Gene[core_position$Module == "magenta"],
    "MMmagenta"
  ]


# ------------------------------------------------------------
# 7. Calculate module-internal MM rank and percentile
#
# Rank is based on absolute MM, descending.
# This reproduces the previously LOCKED ranking framework.
# ------------------------------------------------------------

core_position$MM_rank <- NA_integer_

core_position$MM_percentile <- NA_real_


# -----------------------------
# Blue / ME2
# -----------------------------

blue_members <- names(
  moduleColors[
    moduleColors == "blue"
  ]
)

blue_MM <- MM_all[
  blue_members,
  "MMblue"
]

blue_abs_MM <- abs(
  blue_MM
)

blue_order <- order(
  blue_abs_MM,
  decreasing = TRUE
)

blue_rank <- integer(
  length(blue_members)
)

blue_rank[blue_order] <-
  seq_along(blue_members)

names(blue_rank) <-
  blue_members


# -----------------------------
# Magenta / ME9
# -----------------------------

magenta_members <- names(
  moduleColors[
    moduleColors == "magenta"
  ]
)

magenta_MM <- MM_all[
  magenta_members,
  "MMmagenta"
]

magenta_abs_MM <- abs(
  magenta_MM
)

magenta_order <- order(
  magenta_abs_MM,
  decreasing = TRUE
)

magenta_rank <- integer(
  length(magenta_members)
)

magenta_rank[magenta_order] <-
  seq_along(magenta_members)

names(magenta_rank) <-
  magenta_members


# ------------------------------------------------------------
# 8. Assign ranks to core genes
# ------------------------------------------------------------

for (i in seq_len(nrow(core_position))) {
  
  gene_i <- core_position$Gene[i]
  
  if (core_position$Module[i] == "blue") {
    
    core_position$MM_rank[i] <-
      blue_rank[gene_i]
    
    core_position$MM_percentile[i] <-
      100 *
      (
        1 -
          (
            core_position$MM_rank[i] - 1
          ) /
          (
            module_size_lookup[
              module_size_lookup$Module == "blue",
              "Module_Size"
            ] - 1
          )
      )
    
  }
  
  if (core_position$Module[i] == "magenta") {
    
    core_position$MM_rank[i] <-
      magenta_rank[gene_i]
    
    core_position$MM_percentile[i] <-
      100 *
      (
        1 -
          (
            core_position$MM_rank[i] - 1
          ) /
          (
            module_size_lookup[
              module_size_lookup$Module == "magenta",
              "Module_Size"
            ] - 1
          )
      )
  }
}


# ------------------------------------------------------------
# 9. Retrieve gene-level GS
# ------------------------------------------------------------

core_position$GS_Immune <-
  GS_all[
    core_position$Gene,
    "GS_Immune"
  ]

core_position$GS_Stromal <-
  GS_all[
    core_position$Gene,
    "GS_Stromal"
  ]

core_position$GS_ESTIMATE <-
  GS_all[
    core_position$Gene,
    "GS_ESTIMATE"
  ]


# ------------------------------------------------------------
# 10. Add module size
# ------------------------------------------------------------

core_position$Module_Size <-
  module_size_lookup$Module_Size[
    match(
      core_position$Module,
      module_size_lookup$Module
    )
  ]


# ------------------------------------------------------------
# 11. Create final long-form module-level TME table
# ------------------------------------------------------------

core_module_tme_long <- merge(
  core_position[
    ,
    c(
      "Symbol",
      "Gene",
      "Module",
      "ME",
      "Module_Size",
      "MM",
      "MM_rank",
      "MM_percentile",
      "GS_Immune",
      "GS_Stromal",
      "GS_ESTIMATE"
    )
  ],
  module_tme_lookup,
  by = c(
    "Module",
    "ME"
  ),
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------
# 12. Add corresponding MM–GS module-level rho
# ------------------------------------------------------------

core_module_tme_long <- merge(
  core_module_tme_long,
  mm_gs_lookup[
    ,
    c(
      "Module",
      "ME",
      "Trait",
      "MM_GS_rho",
      "MM_GS_FDR"
    )
  ],
  by = c(
    "Module",
    "ME",
    "Trait"
  ),
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------
# 13. Add gene-specific GS corresponding to each trait
# ------------------------------------------------------------

core_module_tme_long$Gene_GS <- NA_real_


core_module_tme_long$Gene_GS[
  core_module_tme_long$Trait ==
    "ImmuneScore"
] <-
  core_module_tme_long$GS_Immune[
    core_module_tme_long$Trait ==
      "ImmuneScore"
  ]

core_module_tme_long$Gene_GS[
  core_module_tme_long$Trait ==
    "StromalScore"
] <-
  core_module_tme_long$GS_Stromal[
    core_module_tme_long$Trait ==
      "StromalScore"
  ]

core_module_tme_long$Gene_GS[
  core_module_tme_long$Trait ==
    "ESTIMATEScore"
] <-
  core_module_tme_long$GS_ESTIMATE[
    core_module_tme_long$Trait ==
      "ESTIMATEScore"
  ]


# ------------------------------------------------------------
# 14. Order rows
# ------------------------------------------------------------

core_module_tme_long$Gene_Order <-
  match(
    core_module_tme_long$Symbol,
    c(
      "AOC3",
      "MS4A7",
      "IL7R"
    )
  )

core_module_tme_long$Trait_Order <-
  match(
    core_module_tme_long$Trait,
    c(
      "ImmuneScore",
      "StromalScore",
      "ESTIMATEScore"
    )
  )

core_module_tme_long <-
  core_module_tme_long[
    order(
      core_module_tme_long$Gene_Order,
      core_module_tme_long$Trait_Order
    ),
    ,
    drop = FALSE
  ]

core_module_tme_long$Gene_Order <- NULL
core_module_tme_long$Trait_Order <- NULL


# ------------------------------------------------------------
# 15. Create requested compact final table
# ------------------------------------------------------------

core_position_final <- core_position[
  ,
  c(
    "Symbol",
    "Gene",
    "Module",
    "ME",
    "Module_Size",
    "MM",
    "MM_percentile",
    "MM_rank",
    "GS_Immune",
    "GS_Stromal",
    "GS_ESTIMATE"
  )
]


# ------------------------------------------------------------
# 16. Print core-gene positioning table
# ------------------------------------------------------------

cat("============================================================\n")
cat("Core-gene integrated positioning\n")
cat("============================================================\n\n")

print(
  core_position_final,
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Print module-level TME / MM–GS integration
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Core-gene × module-level TME / MM–GS integration\n")
cat("============================================================\n\n")

print(
  core_module_tme_long[
    ,
    c(
      "Symbol",
      "Module",
      "ME",
      "Module_Size",
      "MM",
      "MM_percentile",
      "MM_rank",
      "Trait",
      "Gene_GS",
      "Module_TME_rho",
      "MM_GS_rho",
      "MM_GS_FDR"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 14F-4D-2D FIX
# Recalculate within-module MM rank and percentile
# using actual module members
# ============================================================

# ------------------------------------------------------------
# 1. Basic validation
# ------------------------------------------------------------

stopifnot(
  nrow(MM_GS_all) == ncol(datExpr_MAD)
)

stopifnot(
  identical(
    rownames(MM_GS_all),
    colnames(datExpr_MAD)
  )
)

stopifnot(
  all(core_info$Gene %in% rownames(MM_GS_all))
)

# ------------------------------------------------------------
# 2. Function for within-module MM ranking
# ------------------------------------------------------------

get_module_mm_rank <- function(
    module_name,
    me_column,
    module_colors,
    MM_matrix
) {
  
  # Actual genes belonging to this module
  module_genes <- colnames(datExpr_MAD)[
    module_colors == module_name
  ]
  
  # Extract MM for actual module members only
  mm_values <- MM_matrix[
    module_genes,
    me_column
  ]
  
  # Rank by absolute MM, descending
  rank_order <- order(
    abs(mm_values),
    decreasing = TRUE,
    na.last = NA
  )
  
  ranked_genes <- module_genes[rank_order]
  
  ranked_mm <- mm_values[rank_order]
  
  module_size <- length(ranked_genes)
  
  rank_df <- data.frame(
    Gene = ranked_genes,
    MM = as.numeric(ranked_mm),
    MM_abs = abs(as.numeric(ranked_mm)),
    MM_rank = seq_along(ranked_genes),
    MM_percentile =
      100 * (
        1 -
          (seq_along(ranked_genes) - 1) /
          (module_size - 1)
      ),
    stringsAsFactors = FALSE
  )
  
  rank_df$Module <- module_name
  rank_df$Module_Size <- module_size
  
  rank_df
}

# ------------------------------------------------------------
# 3. Calculate TRUE module-internal rankings
# ------------------------------------------------------------

blue_rank_df <- get_module_mm_rank(
  module_name = "blue",
  me_column = "MMblue",
  module_colors = moduleColors,
  MM_matrix = MM_all
)

magenta_rank_df <- get_module_mm_rank(
  module_name = "magenta",
  me_column = "MMmagenta",
  module_colors = moduleColors,
  MM_matrix = MM_all
)

# ------------------------------------------------------------
# 4. Check module sizes
# ------------------------------------------------------------

cat("Blue module size   :", nrow(blue_rank_df), "\n")
cat("Magenta module size:", nrow(magenta_rank_df), "\n")

# ------------------------------------------------------------
# 5. Extract core-gene positions
# ------------------------------------------------------------

core_mm_rank <- rbind(
  blue_rank_df[
    blue_rank_df$Gene %in% core_info$Gene[
      core_info$Module == "blue"
    ],
    ,
    drop = FALSE
  ],
  
  magenta_rank_df[
    magenta_rank_df$Gene %in% core_info$Gene[
      core_info$Module == "magenta"
    ],
    ,
    drop = FALSE
  ]
)

# ------------------------------------------------------------
# 6. Add gene symbols
# ------------------------------------------------------------

core_mm_rank$Symbol <- core_info$Symbol[
  match(
    core_mm_rank$Gene,
    core_info$Gene
  )
]

# ------------------------------------------------------------
# 7. Reorder
# ------------------------------------------------------------

core_mm_rank <- core_mm_rank[
  match(
    core_info$Symbol,
    core_mm_rank$Symbol
  ),
  ,
  drop = FALSE
]

rownames(core_mm_rank) <- NULL

# ------------------------------------------------------------
# 8. Print core-gene MM positions
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("TRUE within-module MM ranking of core genes\n")
cat("============================================================\n\n")

print(
  core_mm_rank[
    ,
    c(
      "Symbol",
      "Gene",
      "Module",
      "Module_Size",
      "MM",
      "MM_rank",
      "MM_percentile"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 14F-4D-2D
# Rebuild final core-gene integrated positioning table
# ============================================================

# ------------------------------------------------------------
# 1. Merge MM rank information with locked MM/GS information
# ------------------------------------------------------------

core_position_final <- core_mm_rank[
  ,
  c(
    "Symbol",
    "Gene",
    "Module",
    "Module_Size",
    "MM",
    "MM_percentile",
    "MM_rank"
  )
]

# Add ME label
core_position_final$ME <- ifelse(
  core_position_final$Module == "blue",
  "ME2",
  ifelse(
    core_position_final$Module == "magenta",
    "ME9",
    NA_character_
  )
)

# ------------------------------------------------------------
# 2. Add locked GS values
# ------------------------------------------------------------

gs_match <- match(
  core_position_final$Gene,
  rownames(GS_all)
)

core_position_final$GS_Immune <-
  GS_all$GS_Immune[gs_match]

core_position_final$GS_Stromal <-
  GS_all$GS_Stromal[gs_match]

core_position_final$GS_ESTIMATE <-
  GS_all$GS_ESTIMATE[gs_match]

# ------------------------------------------------------------
# 3. Reorder columns
# ------------------------------------------------------------

core_position_final <- core_position_final[
  ,
  c(
    "Symbol",
    "Gene",
    "Module",
    "ME",
    "Module_Size",
    "MM",
    "MM_percentile",
    "MM_rank",
    "GS_Immune",
    "GS_Stromal",
    "GS_ESTIMATE"
  )
]

# ------------------------------------------------------------
# 4. Final validation
# ------------------------------------------------------------

stopifnot(
  all(!is.na(core_position_final$MM_rank))
)

stopifnot(
  all(!is.na(core_position_final$MM_percentile))
)

stopifnot(
  core_position_final$MM_rank[
    core_position_final$Symbol == "AOC3"
  ] == 7
)

stopifnot(
  core_position_final$MM_rank[
    core_position_final$Symbol == "MS4A7"
  ] == 267
)

stopifnot(
  core_position_final$MM_rank[
    core_position_final$Symbol == "IL7R"
  ] == 120
)

# ------------------------------------------------------------
# 5. Print final table
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL core-gene integrated positioning\n")
cat("============================================================\n\n")

print(
  core_position_final,
  row.names = FALSE
)


# ============================================================
# 14F-4D-2E-2 FINAL FIX
# Direct module-level mapping using match()
# No merge() is used
# ============================================================

# ------------------------------------------------------------
# 1. Start from the validated core-position table
# ------------------------------------------------------------

core_final_integrated <- core_position_final[
  ,
  c(
    "Symbol",
    "Gene",
    "Module",
    "ME",
    "Module_Size",
    "MM",
    "MM_percentile",
    "MM_rank",
    "GS_Immune",
    "GS_Stromal",
    "GS_ESTIMATE"
  ),
  drop = FALSE
]

# ------------------------------------------------------------
# 2. Module-level TME correlations
#    LOCKED values from 14F-4D-2C
# ------------------------------------------------------------

module_tme_lookup_correct <- data.frame(
  
  Module = c(
    "blue",
    "magenta"
  ),
  
  Module_TME_Immune_rho = c(
    0.9689631,
    0.4209524
  ),
  
  Module_TME_Stromal_rho = c(
    0.7081174,
    0.6451482
  ),
  
  Module_TME_ESTIMATE_rho = c(
    0.9179035,
    0.5695171
  ),
  
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 3. Module-level MM-GS correlations
#    LOCKED values from 14F-4D-2C
# ------------------------------------------------------------

mmgs_lookup_correct <- data.frame(
  
  Module = c(
    "blue",
    "magenta"
  ),
  
  MM_GS_Immune_rho = c(
    0.9433047,
    0.4023558
  ),
  
  MM_GS_Stromal_rho = c(
    0.7068876,
    0.5221497
  ),
  
  MM_GS_ESTIMATE_rho = c(
    0.9081796,
    0.5088140
  ),
  
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 4. Match module-level TME values by Module
# ------------------------------------------------------------

idx_tme <- match(
  core_final_integrated$Module,
  module_tme_lookup_correct$Module
)

stopifnot(
  !any(is.na(idx_tme))
)

core_final_integrated$Module_TME_Immune_rho <-
  module_tme_lookup_correct$Module_TME_Immune_rho[
    idx_tme
  ]

core_final_integrated$Module_TME_Stromal_rho <-
  module_tme_lookup_correct$Module_TME_Stromal_rho[
    idx_tme
  ]

core_final_integrated$Module_TME_ESTIMATE_rho <-
  module_tme_lookup_correct$Module_TME_ESTIMATE_rho[
    idx_tme
  ]

# ------------------------------------------------------------
# 5. Match module-level MM-GS values by Module
# ------------------------------------------------------------

idx_mmgs <- match(
  core_final_integrated$Module,
  mmgs_lookup_correct$Module
)

stopifnot(
  !any(is.na(idx_mmgs))
)

core_final_integrated$MM_GS_Immune_rho <-
  mmgs_lookup_correct$MM_GS_Immune_rho[
    idx_mmgs
  ]

core_final_integrated$MM_GS_Stromal_rho <-
  mmgs_lookup_correct$MM_GS_Stromal_rho[
    idx_mmgs
  ]

core_final_integrated$MM_GS_ESTIMATE_rho <-
  mmgs_lookup_correct$MM_GS_ESTIMATE_rho[
    idx_mmgs
  ]

# ------------------------------------------------------------
# 6. Restore exact core-gene order
# ------------------------------------------------------------

core_order <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_final_integrated <- core_final_integrated[
  match(
    core_order,
    core_final_integrated$Symbol
  ),
  ,
  drop = FALSE
]

rownames(core_final_integrated) <- NULL

# ------------------------------------------------------------
# 7. Add descriptive TME profiles
# ------------------------------------------------------------

core_final_integrated$Mean_abs_GS <- rowMeans(
  abs(
    core_final_integrated[
      ,
      c(
        "GS_Immune",
        "GS_Stromal",
        "GS_ESTIMATE"
      )
    ]
  )
)

core_final_integrated$Mean_abs_Module_TME_rho <- rowMeans(
  abs(
    core_final_integrated[
      ,
      c(
        "Module_TME_Immune_rho",
        "Module_TME_Stromal_rho",
        "Module_TME_ESTIMATE_rho"
      )
    ]
  )
)

# ------------------------------------------------------------
# 8. Final column order
# ------------------------------------------------------------

core_final_integrated <- core_final_integrated[
  ,
  c(
    "Symbol",
    "Gene",
    "Module",
    "ME",
    "Module_Size",
    "MM",
    "MM_percentile",
    "MM_rank",
    "GS_Immune",
    "GS_Stromal",
    "GS_ESTIMATE",
    "Module_TME_Immune_rho",
    "Module_TME_Stromal_rho",
    "Module_TME_ESTIMATE_rho",
    "MM_GS_Immune_rho",
    "MM_GS_Stromal_rho",
    "MM_GS_ESTIMATE_rho",
    "Mean_abs_GS",
    "Mean_abs_Module_TME_rho"
),
drop = FALSE
]

# ------------------------------------------------------------
# 9. Print final table
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL CORE-GENE INTEGRATED TABLE\n")
cat("============================================================\n\n")

print(
  core_final_integrated,
  row.names = FALSE
)

# ============================================================
# 14F-4D-2E-2 FINAL VALIDATION
# ============================================================

# Core genes
stopifnot(
  identical(
    core_final_integrated$Symbol,
    c("AOC3", "MS4A7", "IL7R")
  )
)

# Module localization
stopifnot(
  core_final_integrated$Module[
    core_final_integrated$Symbol == "AOC3"
  ] == "magenta"
)

stopifnot(
  all(
    core_final_integrated$Module[
      core_final_integrated$Symbol %in%
        c("MS4A7", "IL7R")
    ] == "blue"
  )
)

# MM ranks
stopifnot(
  core_final_integrated$MM_rank[
    core_final_integrated$Symbol == "AOC3"
  ] == 7
)

stopifnot(
  core_final_integrated$MM_rank[
    core_final_integrated$Symbol == "MS4A7"
  ] == 267
)

stopifnot(
  core_final_integrated$MM_rank[
    core_final_integrated$Symbol == "IL7R"
  ] == 120
)

# Module-level TME mapping
stopifnot(
  core_final_integrated$Module_TME_Immune_rho[
    core_final_integrated$Symbol == "AOC3"
  ] == 0.4209524
)

stopifnot(
  core_final_integrated$Module_TME_Immune_rho[
    core_final_integrated$Symbol == "MS4A7"
  ] == 0.9689631
)

stopifnot(
  core_final_integrated$Module_TME_Immune_rho[
    core_final_integrated$Symbol == "IL7R"
  ] == 0.9689631
)

# Blue module must be identical for MS4A7 and IL7R
stopifnot(
  identical(
    core_final_integrated$Module_TME_Immune_rho[
      core_final_integrated$Symbol == "MS4A7"
    ],
    core_final_integrated$Module_TME_Immune_rho[
      core_final_integrated$Symbol == "IL7R"
    ]
  )
)

# Expected mean module-level TME associations
stopifnot(
  abs(
    core_final_integrated$Mean_abs_Module_TME_rho[
      core_final_integrated$Symbol == "AOC3"
    ] - 0.5452059
  ) < 1e-6
)

stopifnot(
  all(
    abs(
      core_final_integrated$Mean_abs_Module_TME_rho[
        core_final_integrated$Symbol %in%
          c("MS4A7", "IL7R")
      ] - 0.8649947
    ) < 1e-6
  )
)

cat("\n============================================================\n")
cat("14F-4D-2E-2 FINAL VALIDATION: PASSED\n")
cat("============================================================\n")


# ============================================================

# 14F-4D-2E-2 FINAL FIXED VERSION
# Core-gene integrated positioning table
# ============================================================

# IMPORTANT:
# Rebuild from LOCKED core_position_final.
# Do NOT use the previously merged core_final_integrated.

# ------------------------------------------------------------
# 1. Start from LOCKED core-gene positioning table
# ------------------------------------------------------------

core_final_integrated <- core_position_final[
  ,
  c(
    "Symbol",
    "Gene",
    "Module",
    "ME",
    "Module_Size",
    "MM",
    "MM_percentile",
    "MM_rank",
    "GS_Immune",
    "GS_Stromal",
    "GS_ESTIMATE"
  ),
  drop = FALSE
]

# ------------------------------------------------------------
# 2. Match module-level TME correlations
# ------------------------------------------------------------

idx_tme <- match(
  core_final_integrated$Module,
  module_tme_lookup_correct$Module
)

core_final_integrated$Module_TME_Immune_rho <-
  module_tme_lookup_correct$Module_TME_Immune_rho[
    idx_tme
  ]

core_final_integrated$Module_TME_Stromal_rho <-
  module_tme_lookup_correct$Module_TME_Stromal_rho[
    idx_tme
  ]

core_final_integrated$Module_TME_ESTIMATE_rho <-
  module_tme_lookup_correct$Module_TME_ESTIMATE_rho[
    idx_tme
  ]

# ------------------------------------------------------------
# 3. Match module-level MM-GS correlations
# ------------------------------------------------------------

idx_mmgs <- match(
  core_final_integrated$Module,
  mmgs_lookup_correct$Module
)

core_final_integrated$MM_GS_Immune_rho <-
  mmgs_lookup_correct$MM_GS_Immune_rho[
    idx_mmgs
  ]

core_final_integrated$MM_GS_Stromal_rho <-
  mmgs_lookup_correct$MM_GS_Stromal_rho[
    idx_mmgs
  ]

core_final_integrated$MM_GS_ESTIMATE_rho <-
  mmgs_lookup_correct$MM_GS_ESTIMATE_rho[
    idx_mmgs
  ]

# ------------------------------------------------------------
# 4. Restore explicit core-gene order
# ------------------------------------------------------------

core_order <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_final_integrated <- core_final_integrated[
  match(
    core_order,
    core_final_integrated$Symbol
  ),
  ,
  drop = FALSE
]

rownames(core_final_integrated) <- NULL

# ------------------------------------------------------------
# 5. Check columns
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Final columns\n")
cat("============================================================\n\n")

print(
  names(core_final_integrated)
)

# ------------------------------------------------------------
# 6. Print final integrated table
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL CORE-GENE INTEGRATED TABLE\n")
cat("============================================================\n\n")

print(
  core_final_integrated,
  row.names = FALSE
)


# ============================================================
# 14F-4D-2E FINAL CORE TABLE VALIDATION
# ============================================================

# ------------------------------------------------------------
# Basic structure
# ------------------------------------------------------------

stopifnot(
  nrow(core_final_integrated) == 3
)

stopifnot(
  identical(
    core_final_integrated$Symbol,
    c(
      "AOC3",
      "MS4A7",
      "IL7R"
    )
  )
)

# ------------------------------------------------------------
# Module localization
# ------------------------------------------------------------

stopifnot(
  core_final_integrated$Module[
    core_final_integrated$Symbol == "AOC3"
  ] == "magenta"
)

stopifnot(
  core_final_integrated$Module[
    core_final_integrated$Symbol == "MS4A7"
  ] == "blue"
)

stopifnot(
  core_final_integrated$Module[
    core_final_integrated$Symbol == "IL7R"
  ] == "blue"
)

# ------------------------------------------------------------
# Module size
# ------------------------------------------------------------

stopifnot(
  core_final_integrated$Module_Size[
    core_final_integrated$Symbol == "AOC3"
  ] == 534
)

stopifnot(
  all(
    core_final_integrated$Module_Size[
      core_final_integrated$Symbol %in%
        c("MS4A7", "IL7R")
    ] == 1798
  )
)

# ------------------------------------------------------------
# MM rank
# ------------------------------------------------------------

stopifnot(
  core_final_integrated$MM_rank[
    core_final_integrated$Symbol == "AOC3"
  ] == 7
)

stopifnot(
  core_final_integrated$MM_rank[
    core_final_integrated$Symbol == "MS4A7"
  ] == 267
)

stopifnot(
  core_final_integrated$MM_rank[
    core_final_integrated$Symbol == "IL7R"
  ] == 120
)

# ------------------------------------------------------------
# MM percentile
# ------------------------------------------------------------

stopifnot(
  abs(
    core_final_integrated$MM_percentile[
      core_final_integrated$Symbol == "AOC3"
    ] - 98.87430
  ) < 1e-5
)

stopifnot(
  abs(
    core_final_integrated$MM_percentile[
      core_final_integrated$Symbol == "MS4A7"
    ] - 85.19755
  ) < 1e-5
)

stopifnot(
  abs(
    core_final_integrated$MM_percentile[
      core_final_integrated$Symbol == "IL7R"
    ] - 93.37785
  ) < 1e-5
)

# ------------------------------------------------------------
# Module-level TME validation
# ------------------------------------------------------------

# AOC3 / magenta
stopifnot(
  all(
    abs(
      core_final_integrated[
        core_final_integrated$Symbol == "AOC3",
        c(
          "Module_TME_Immune_rho",
          "Module_TME_Stromal_rho",
          "Module_TME_ESTIMATE_rho"
        )
      ] -
        c(
          0.4209524,
          0.6451482,
          0.5695171
        )
    ) < 1e-7
  )
)

# MS4A7 / blue
stopifnot(
  all(
    abs(
      core_final_integrated[
        core_final_integrated$Symbol == "MS4A7",
        c(
          "Module_TME_Immune_rho",
          "Module_TME_Stromal_rho",
          "Module_TME_ESTIMATE_rho"
        )
      ] -
        c(
          0.9689631,
          0.7081174,
          0.9179035
        )
    ) < 1e-7
  )
)

# IL7R / blue
stopifnot(
  all(
    abs(
      core_final_integrated[
        core_final_integrated$Symbol == "IL7R",
        c(
          "Module_TME_Immune_rho",
          "Module_TME_Stromal_rho",
          "Module_TME_ESTIMATE_rho"
        )
      ] -
        c(
          0.9689631,
          0.7081174,
          0.9179035
        )
    ) < 1e-7
  )
)

# ------------------------------------------------------------
# MM-GS module-level validation
# ------------------------------------------------------------

stopifnot(
  all(
    abs(
      core_final_integrated[
        core_final_integrated$Symbol == "AOC3",
        c(
          "MM_GS_Immune_rho",
          "MM_GS_Stromal_rho",
          "MM_GS_ESTIMATE_rho"
        )
      ] -
        c(
          0.4023558,
          0.5221497,
          0.5088140
        )
    ) < 1e-7
  )
)

stopifnot(
  all(
    abs(
      core_final_integrated[
        core_final_integrated$Symbol == "MS4A7",
        c(
          "MM_GS_Immune_rho",
          "MM_GS_Stromal_rho",
          "MM_GS_ESTIMATE_rho"
        )
      ] -
        c(
          0.9433047,
          0.7068876,
          0.9081796
        )
    ) < 1e-7
  )
)

stopifnot(
  all(
    abs(
      core_final_integrated[
        core_final_integrated$Symbol == "IL7R",
        c(
          "MM_GS_Immune_rho",
          "MM_GS_Stromal_rho",
          "MM_GS_ESTIMATE_rho"
        )
      ] -
        c(
          0.9433047,
          0.7068876,
          0.9081796
        )
    ) < 1e-7
  )
)

cat("\n============================================================\n")
cat("14F-4D-2E CORE TABLE VALIDATION: PASSED\n")
cat("============================================================\n")


# ============================================================
# 14F-4D-2E-3
# Descriptive TME relevance profile
# ============================================================

core_final_integrated$Mean_abs_GS <- rowMeans(
  abs(
    core_final_integrated[
      ,
      c(
        "GS_Immune",
        "GS_Stromal",
        "GS_ESTIMATE"
      )
    ]
  )
)

core_final_integrated$Mean_abs_Module_TME_rho <- rowMeans(
  abs(
    core_final_integrated[
      ,
      c(
        "Module_TME_Immune_rho",
        "Module_TME_Stromal_rho",
        "Module_TME_ESTIMATE_rho"
      )
    ]
  )
)

cat("\n============================================================\n")
cat("DESCRIPTIVE CORE-GENE TME PROFILE\n")
cat("============================================================\n\n")

print(
  core_final_integrated[
    ,
    c(
      "Symbol",
      "Module",
      "MM_percentile",
      "MM_rank",
      "GS_Immune",
      "GS_Stromal",
      "GS_ESTIMATE",
      "Mean_abs_GS",
      "Mean_abs_Module_TME_rho"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 14F-4D-2E-4
# Survival-associated Brown / ME3 reference
# ============================================================

brown_survival_reference <- data.frame(
  
  Module = "brown",
  ME = "ME3",
  
  Module_Size = sum(
    moduleColors == "brown"
  ),
  
  Survival_N = 485,
  Events = 178,
  
  Cox_HR_1SD = 1.804223,
  Cox_P = 0.0023,
  
  TimeVarying_Term_HR = 0.674764,
  TimeVarying_Term_P = 0.0720,
  
  C_index = 0.6978267,
  C_index_SE = 0.0223648,
  
  AIC = 1822.820,
  
  LRT_P = 5.767664e-13,
  
  ImmuneScore_rho = -0.1621890,
  StromalScore_rho = -0.1924747,
  ESTIMATEScore_rho = -0.1863114,
  
  Survival_role =
    "Survival-associated reference module",
  
  TME_role =
    "Weak inverse TME association",
  
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("BROWN / ME3 SURVIVAL REFERENCE\n")
cat("============================================================\n\n")

print(
  brown_survival_reference,
  row.names = FALSE
)


# ============================================================
# 14F-4D-2E-5
# Final manuscript-level analytical summary
# ============================================================

manuscript_level_summary <- data.frame(
  
  Feature = c(
    "AOC3",
    "MS4A7",
    "IL7R",
    "Blue / ME2",
    "Magenta / ME9",
    "Brown / ME3"
  ),
  
  Module_or_role = c(
    "Magenta module",
    "Blue module",
    "Blue module",
    "Core TME-associated module",
    "Core TME-associated module",
    "Survival-associated reference module"
  ),
  
  Key_finding = c(
    "Highest within-module MM percentile; relatively stronger stromal GS",
    "High MM and strong TME associations",
    "High MM and strongest gene-level TME associations",
    "Strongest module-level TME association and MM–GS concordance",
    "Positive TME association, particularly stromal",
    "Survival-associated with weak inverse TME associations"
  ),
  
  Interpretation = c(
    "Strong module positioning with stromal relevance",
    "TME-associated blue-module component",
    "Strongest gene-level TME relevance",
    "Immune/TME-associated transcriptional module",
    "Stromal-associated transcriptional module",
    "Survival reference rather than core-gene module"
  ),
  
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("MANUSCRIPT-LEVEL ANALYTICAL SUMMARY\n")
cat("============================================================\n\n")

print(
  manuscript_level_summary,
  row.names = FALSE
)


# ============================================================
# FINAL VALIDATION OF 14F-4D-2E
# ============================================================

# Core genes
stopifnot(
  identical(
    core_final_integrated$Symbol,
    c("AOC3", "MS4A7", "IL7R")
  )
)

# Module localization
stopifnot(
  core_final_integrated$Module[
    core_final_integrated$Symbol == "AOC3"
  ] == "magenta"
)

stopifnot(
  all(
    core_final_integrated$Module[
      core_final_integrated$Symbol %in% c(
        "MS4A7",
        "IL7R"
      )
    ] == "blue"
  )
)

# Locked MM ranks
stopifnot(
  core_final_integrated$MM_rank[
    core_final_integrated$Symbol == "AOC3"
  ] == 7
)

stopifnot(
  core_final_integrated$MM_rank[
    core_final_integrated$Symbol == "MS4A7"
  ] == 267
)

stopifnot(
  core_final_integrated$MM_rank[
    core_final_integrated$Symbol == "IL7R"
  ] == 120
)

# Locked MM percentiles
stopifnot(
  abs(
    core_final_integrated$MM_percentile[
      core_final_integrated$Symbol == "AOC3"
    ] - 98.87430
  ) < 1e-5
)

stopifnot(
  abs(
    core_final_integrated$MM_percentile[
      core_final_integrated$Symbol == "MS4A7"
    ] - 85.19755
  ) < 1e-5
)

stopifnot(
  abs(
    core_final_integrated$MM_percentile[
      core_final_integrated$Symbol == "IL7R"
    ] - 93.37785
  ) < 1e-5
)

# Brown survival reference
stopifnot(
  brown_survival_reference$Survival_N == 485
)

stopifnot(
  brown_survival_reference$Events == 178
)

cat("\n============================================================\n")
cat("14F-4D-2E FINAL VALIDATION: PASSED\n")
cat("============================================================\n")


# ============================================================
# 14F-4D-2E-SAVE
# Final integrated objects for downstream visualization
# ============================================================

# ------------------------------------------------------------
# 0. Basic validation
# ------------------------------------------------------------

stopifnot(
  all(
    core_position_final$Symbol %in%
      c("AOC3", "MS4A7", "IL7R")
  )
)

stopifnot(
  identical(
    core_position_final$Symbol,
    c("AOC3", "MS4A7", "IL7R")
  )
)

# ============================================================
# 1. Correct module-level lookup
# ============================================================

module_tme_lookup_save <- data.frame(
  
  Module = c(
    "blue",
    "magenta"
  ),
  
  ME = c(
    "ME2",
    "ME9"
  ),
  
  Module_TME_Immune_rho = c(
    0.9689631,
    0.4209524
  ),
  
  Module_TME_Stromal_rho = c(
    0.7081174,
    0.6451482
  ),
  
  Module_TME_ESTIMATE_rho = c(
    0.9179035,
    0.5695171
  ),
  
  MM_GS_Immune_rho = c(
    0.9433047,
    0.4023558
  ),
  
  MM_GS_Stromal_rho = c(
    0.7068876,
    0.5221497
  ),
  
  MM_GS_ESTIMATE_rho = c(
    0.9081796,
    0.5088140
  ),
  
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Check lookup
# ------------------------------------------------------------

print(module_tme_lookup_save)


# ============================================================
# 2. Build final core-gene integrated table
# ============================================================

core_final_integrated_save <- core_position_final

# Match each core gene's module to the lookup table
module_idx <- match(
  core_final_integrated_save$Module,
  module_tme_lookup_save$Module
)

# ------------------------------------------------------------
# Add module-level TME correlations
# ------------------------------------------------------------

core_final_integrated_save$Module_TME_Immune_rho <-
  module_tme_lookup_save$Module_TME_Immune_rho[module_idx]

core_final_integrated_save$Module_TME_Stromal_rho <-
  module_tme_lookup_save$Module_TME_Stromal_rho[module_idx]

core_final_integrated_save$Module_TME_ESTIMATE_rho <-
  module_tme_lookup_save$Module_TME_ESTIMATE_rho[module_idx]

# ------------------------------------------------------------
# Add module-level MM-GS correlations
# ------------------------------------------------------------

core_final_integrated_save$MM_GS_Immune_rho <-
  module_tme_lookup_save$MM_GS_Immune_rho[module_idx]

core_final_integrated_save$MM_GS_Stromal_rho <-
  module_tme_lookup_save$MM_GS_Stromal_rho[module_idx]

core_final_integrated_save$MM_GS_ESTIMATE_rho <-
  module_tme_lookup_save$MM_GS_ESTIMATE_rho[module_idx]

# ------------------------------------------------------------
# Reorder columns
# ------------------------------------------------------------

core_final_integrated_save <- core_final_integrated_save[
  ,
  c(
    "Symbol",
    "Gene",
    "Module",
    "ME",
    "Module_Size",
    "MM",
    "MM_percentile",
    "MM_rank",
    "GS_Immune",
    "GS_Stromal",
    "GS_ESTIMATE",
    "Module_TME_Immune_rho",
    "Module_TME_Stromal_rho",
    "Module_TME_ESTIMATE_rho",
    "MM_GS_Immune_rho",
    "MM_GS_Stromal_rho",
    "MM_GS_ESTIMATE_rho"
  )
]

# ------------------------------------------------------------
# Restore exact gene order
# ------------------------------------------------------------

core_order <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_final_integrated_save <- core_final_integrated_save[
  match(
    core_order,
    core_final_integrated_save$Symbol
  ),
  ,
  drop = FALSE
]

rownames(core_final_integrated_save) <- NULL


# ============================================================
# 3. Validate final integrated core-gene table
# ============================================================

# ------------------------------------------------------------
# Core gene order
# ------------------------------------------------------------

stopifnot(
  identical(
    core_final_integrated_save$Symbol,
    c("AOC3", "MS4A7", "IL7R")
  )
)

# ------------------------------------------------------------
# Module localization
# ------------------------------------------------------------

stopifnot(
  core_final_integrated_save$Module[
    core_final_integrated_save$Symbol == "AOC3"
  ] == "magenta"
)

stopifnot(
  core_final_integrated_save$Module[
    core_final_integrated_save$Symbol == "MS4A7"
  ] == "blue"
)

stopifnot(
  core_final_integrated_save$Module[
    core_final_integrated_save$Symbol == "IL7R"
  ] == "blue"
)

# ------------------------------------------------------------
# MM rank
# ------------------------------------------------------------

stopifnot(
  core_final_integrated_save$MM_rank[
    core_final_integrated_save$Symbol == "AOC3"
  ] == 7
)

stopifnot(
  core_final_integrated_save$MM_rank[
    core_final_integrated_save$Symbol == "MS4A7"
  ] == 267
)

stopifnot(
  core_final_integrated_save$MM_rank[
    core_final_integrated_save$Symbol == "IL7R"
  ] == 120
)

# ------------------------------------------------------------
# Module-level TME validation
# ------------------------------------------------------------

stopifnot(
  core_final_integrated_save$Module_TME_Immune_rho[
    core_final_integrated_save$Symbol == "AOC3"
  ] == 0.4209524
)

stopifnot(
  core_final_integrated_save$Module_TME_Immune_rho[
    core_final_integrated_save$Symbol == "MS4A7"
  ] == 0.9689631
)

stopifnot(
  core_final_integrated_save$Module_TME_Immune_rho[
    core_final_integrated_save$Symbol == "IL7R"
  ] == 0.9689631
)

# ------------------------------------------------------------
# MM-GS validation
# ------------------------------------------------------------

stopifnot(
  core_final_integrated_save$MM_GS_Immune_rho[
    core_final_integrated_save$Symbol == "AOC3"
  ] == 0.4023558
)

stopifnot(
  core_final_integrated_save$MM_GS_Immune_rho[
    core_final_integrated_save$Symbol == "MS4A7"
  ] == 0.9433047
)

stopifnot(
  core_final_integrated_save$MM_GS_Immune_rho[
    core_final_integrated_save$Symbol == "IL7R"
  ] == 0.9433047
)

cat("\n============================================================\n")
cat("FINAL CORE-GENE INTEGRATED TABLE: VALIDATION PASSED\n")
cat("============================================================\n")


# ============================================================
# 4. Compact core-gene table for visualization
# ============================================================

core_visualization_table <- core_final_integrated_save[
  ,
  c(
    "Symbol",
    "Module",
    "ME",
    "Module_Size",
    "MM",
    "MM_percentile",
    "MM_rank",
    "GS_Immune",
    "GS_Stromal",
    "GS_ESTIMATE",
    "Module_TME_Immune_rho",
    "Module_TME_Stromal_rho",
    "Module_TME_ESTIMATE_rho",
    "MM_GS_Immune_rho",
    "MM_GS_Stromal_rho",
    "MM_GS_ESTIMATE_rho"
  )
]

print(
  core_visualization_table,
  row.names = FALSE
)


# ============================================================
# 5. Brown / ME3 survival reference table
# ============================================================

brown_survival_reference_save <- data.frame(
  
  Module = "brown",
  ME = "ME3",
  
  Module_Size = sum(
    moduleColors == "brown"
  ),
  
  Survival_N = 485,
  Events = 178,
  
  Cox_HR_1SD = 1.804223,
  Cox_P = 0.0023,
  
  TimeVarying_Term_HR = 0.674764,
  TimeVarying_Term_P = 0.0720,
  
  C_index = 0.6978267,
  C_index_SE = 0.0223648,
  
  AIC = 1822.820,
  
  LRT_P = 5.767664e-13,
  
  ImmuneScore_rho = -0.1621890,
  StromalScore_rho = -0.1924747,
  ESTIMATEScore_rho = -0.1863114,
  
  Survival_role =
    "Survival-associated reference module",
  
  TME_role =
    "Weak inverse TME association",
  
  stringsAsFactors = FALSE
)

print(
  brown_survival_reference_save,
  row.names = FALSE
)


# ============================================================
# 6. Final module-level integration table
# ============================================================

final_module_integration_save <- data.frame(
  
  Module = c(
    "blue",
    "magenta",
    "brown"
  ),
  
  ME = c(
    "ME2",
    "ME9",
    "ME3"
  ),
  
  Module_Size = c(
    sum(moduleColors == "blue"),
    sum(moduleColors == "magenta"),
    sum(moduleColors == "brown")
  ),
  
  Core_Genes = c(
    "MS4A7; IL7R",
    "AOC3",
    NA_character_
  ),
  
  ImmuneScore_rho = c(
    0.9689631,
    0.4209524,
    -0.1621890
  ),
  
  StromalScore_rho = c(
    0.7081174,
    0.6451482,
    -0.1924747
  ),
  
  ESTIMATEScore_rho = c(
    0.9179035,
    0.5695171,
    -0.1863114
  ),
  
  MM_GS_Immune_rho = c(
    0.9433047,
    0.4023558,
    NA_real_
  ),
  
  MM_GS_Stromal_rho = c(
    0.7068876,
    0.5221497,
    NA_real_
  ),
  
  MM_GS_ESTIMATE_rho = c(
    0.9081796,
    0.5088140,
    NA_real_
  ),
  
  Survival_reference = c(
    "No",
    "No",
    "Yes"
  ),
  
  TME_pattern = c(
    "Strong positive TME association",
    "Positive TME association, relatively stromal",
    "Weak inverse TME association"
  ),
  
  MM_GS_pattern = c(
    "Very strong positive concordance",
    "Moderate positive concordance",
    NA_character_
  ),
  
  stringsAsFactors = FALSE
)

print(
  final_module_integration_save,
  row.names = FALSE
)


# ============================================================
# 7. Descriptive TME profile
# ============================================================

core_visualization_table$Mean_abs_GS <- rowMeans(
  abs(
    core_visualization_table[
      ,
      c(
        "GS_Immune",
        "GS_Stromal",
        "GS_ESTIMATE"
      )
    ]
  )
)

core_visualization_table$Mean_abs_Module_TME_rho <- rowMeans(
  abs(
    core_visualization_table[
      ,
      c(
        "Module_TME_Immune_rho",
        "Module_TME_Stromal_rho",
        "Module_TME_ESTIMATE_rho"
      )
    ]
  )
)

core_descriptive_profile <- core_visualization_table[
  ,
  c(
    "Symbol",
    "Module",
    "MM_percentile",
    "MM_rank",
    "GS_Immune",
    "GS_Stromal",
    "GS_ESTIMATE",
    "Mean_abs_GS",
    "Mean_abs_Module_TME_rho"
  )
]

print(
  core_descriptive_profile,
  row.names = FALSE
)


# ============================================================
# 8. Save all 14F-4D-2E objects
# ============================================================

save(
  
  # Core-gene tables
  core_position_final,
  core_final_integrated_save,
  core_visualization_table,
  core_descriptive_profile,
  
  # Module-level tables
  final_module_integration_save,
  module_tme_lookup_save,
  
  # Survival reference
  brown_survival_reference_save,
  
  # Original locked supporting objects
  MM_all,
  GS_all,
  moduleColors,
  ME_color_formal_raw,
  
  # WGCNA / TME core objects
  datExpr_MAD,
  estimate_scores_wgcna,
  module_tme_results,
  
  file = "Paper1_WORKFLOW03_14F_4D_2E_FINAL.RData"
)

cat("\n============================================================\n")
cat("14F-4D-2E FINAL OBJECTS SAVED\n")
cat("File: Paper1_WORKFLOW03_14F_4D_2E_FINAL.RData\n")
cat("============================================================\n")


# ============================================================
# 9. Save CSV files for visualization / manuscript preparation
# ============================================================

write.csv(
  core_final_integrated_save,
  "Paper1_14F_4D_2E_core_final_integrated.csv",
  row.names = FALSE
)

write.csv(
  core_visualization_table,
  "Paper1_14F_4D_2E_core_visualization.csv",
  row.names = FALSE
)

write.csv(
  core_descriptive_profile,
  "Paper1_14F_4D_2E_core_descriptive_profile.csv",
  row.names = FALSE
)

write.csv(
  final_module_integration_save,
  "Paper1_14F_4D_2E_module_integration.csv",
  row.names = FALSE
)

write.csv(
  brown_survival_reference_save,
  "Paper1_14F_4D_2E_brown_survival_reference.csv",
  row.names = FALSE
)

cat("\nCSV files saved successfully.\n")


# ============================================================
# Paper 1
# WORKFLOW-03
# 14F-4D-2E FINAL RELOAD
# ============================================================

rm(list = ls())

load(
  "Paper1_WORKFLOW03_14F_4D_2E_FINAL.RData"
)

cat("\n============================================================\n")
cat("Paper 1 WORKFLOW-03 / 14F-4D-2E FINAL OBJECTS LOADED\n")
cat("============================================================\n\n")

cat(
  "Core genes:",
  paste(core_visualization_table$Symbol, collapse = ", "),
  "\n"
)

cat(
  "Modules:",
  paste(
    unique(core_visualization_table$Module),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Blue module size:",
  sum(moduleColors == "blue"),
  "\n"
)

cat(
  "Magenta module size:",
  sum(moduleColors == "magenta"),
  "\n"
)

cat(
  "Brown module size:",
  sum(moduleColors == "brown"),
  "\n"
)

cat("\nReady for downstream visualization.\n")


# ============================================================
# 10. Final save verification
# ============================================================

cat("\n================ FINAL CHECK ================\n\n")

print(
  core_visualization_table,
  row.names = FALSE
)

cat("\n--- Brown survival reference ---\n\n")

print(
  brown_survival_reference_save,
  row.names = FALSE
)

cat("\n--- Module integration ---\n\n")

print(
  final_module_integration_save,
  row.names = FALSE
)

cat("\n=============================================\n")
cat("14F-4D-2E READY FOR TOMORROW'S FIGURES\n")
cat("=============================================\n")



