
# 00_setup.R

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Please install the 'here' package.")
}

project_root <- here::here()

input_dir <- file.path(project_root, "data")
result_dir <- file.path(project_root, "results")
figure_dir <- file.path(project_root, "figures")

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)


setwd("input_dir/GitHub2_Paper1_tumor-only WGCNA")
getwd()

# ============================================================
# Step 14A
# Check existing tumor-only WGCNA objects
# ============================================================

# Check available objects
ls()

# Check dimensions of the WGCNA expression matrix
dim(datExpr)

# Check module colors
table(moduleColors)

# Check module eigengenes
dim(MEs)

# Inspect first few columns
head(MEs)

# Check sample IDs
head(rownames(MEs))

# Check expression matrix orientation
dim(datExpr)
head(rownames(datExpr))
head(colnames(datExpr))

# ============================================================
# Step 14A-1
# Identify the actual tumor-only WGCNA expression object
# ============================================================

dim(datExpr_MAD)

dim(MEs)

length(moduleColors)

dim(gene_MAD)

dim(net)

##then：
head(rownames(datExpr_MAD))
head(colnames(datExpr_MAD))

head(rownames(MEs))
head(colnames(MEs))

head(moduleColors)

class(datExpr_MAD)

typeof(datExpr_MAD)

summary(datExpr_MAD[, 1])

range(datExpr_MAD[, 1], na.rm = TRUE)


anyNA(datExpr_MAD)

sum(is.na(datExpr_MAD))

##Step 14A-2: locate and verify the complete tumor-only expression matrix 
##Workflow-04-2C_expr_primary_tumor_log2_FPKM_UQ.txt.

file.exists(
  "input_dir/Workflow-04-2C_expr_primary_tumor_log2_FPKM_UQ.txt"
)

readLines(
  "input_dir/Workflow-04-2C_expr_primary_tumor_log2_FPKM_UQ.txt",
  n = 5
)

file.info(
  "input_dir/Workflow-04-2C_expr_primary_tumor_log2_FPKM_UQ.txt"
)

# ============================================================
# Step 14A-3
# Verify sample IDs in the complete tumor-only expression file
# ============================================================

expr_file <- "input_dir/Workflow-04-2C_expr_primary_tumor_log2_FPKM_UQ.txt"

header <- strsplit(
  readLines(expr_file, n = 1),
  "\t",
  fixed = TRUE
)[[1]]

# Number of samples
length(header) - 1

# First 10 sample IDs
head(header[-1], 10)

# Last 10 sample IDs
tail(header[-1], 10)

# Compare with WGCNA sample IDs
identical(
  header[-1],
  rownames(datExpr_MAD)
)

# Number of overlapping samples
length(
  intersect(
    header[-1],
    rownames(datExpr_MAD)
  )
)

# Inspect first few genes and expression values
expr_preview <- read.delim(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  nrows = 5
)

dim(expr_preview)

head(
  expr_preview[, 1:min(6, ncol(expr_preview))]
)

# Check gene ID column
head(expr_preview[[1]])

# Check expression values
summary(
  as.numeric(
    unlist(expr_preview[, -1])
  )
)

# ============================================================
# Step 14A-4
# Inspect GENCODE v22 annotation
# ============================================================

annotation_file <- "input_dir/Workflow-04-2C-1A_GENCODE_v22_complete_annotation.txt"

# ------------------------------------------------------------
# 1. Check file existence
# ------------------------------------------------------------

file.exists(annotation_file)

# ------------------------------------------------------------
# 2. Inspect header
# ------------------------------------------------------------

annotation_header <- strsplit(
  readLines(annotation_file, n = 1),
  "\t",
  fixed = TRUE
)[[1]]

annotation_header

length(annotation_header)

# ------------------------------------------------------------
# 3. Read first 10 rows only
# ------------------------------------------------------------

annotation_preview <- read.delim(
  annotation_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  nrows = 10
)

dim(annotation_preview)

head(annotation_preview)

# ------------------------------------------------------------
# 4. Show column names
# ------------------------------------------------------------

colnames(annotation_preview)

# ------------------------------------------------------------
# 5. Show first few values of every column
# ------------------------------------------------------------

lapply(
  annotation_preview,
  head
)

# ============================================================
# Step 14A-5
# QC of Ensembl ID -> Gene Symbol mapping
# ============================================================

# ------------------------------------------------------------
# 1. Read complete GENCODE v22 annotation
# ------------------------------------------------------------

annotation <- read.delim(
  annotation_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

dim(annotation)

# ------------------------------------------------------------
# 2. Basic structure
# ------------------------------------------------------------

str(annotation)

# ------------------------------------------------------------
# 3. Check duplicated Ensembl IDs
# ------------------------------------------------------------

sum(duplicated(annotation$Ensembl_ID))

# ------------------------------------------------------------
# 4. Check missing / empty Ensembl IDs
# ------------------------------------------------------------

sum(
  is.na(annotation$Ensembl_ID) |
    annotation$Ensembl_ID == ""
)

# ------------------------------------------------------------
# 5. Check missing / empty gene symbols
# ------------------------------------------------------------

sum(
  is.na(annotation$gene_name) |
    annotation$gene_name == ""
)

# ------------------------------------------------------------
# 6. Number of unique Ensembl IDs
# ------------------------------------------------------------

length(unique(annotation$Ensembl_ID))

# ------------------------------------------------------------
# 7. Number of unique gene symbols
# ------------------------------------------------------------

length(unique(annotation$gene_name))

# ------------------------------------------------------------
# 8. Number of duplicated gene symbols
# ------------------------------------------------------------

sum(duplicated(annotation$gene_name))

# ------------------------------------------------------------
# 9. Most frequent duplicated gene symbols
# ------------------------------------------------------------

gene_name_freq <- sort(
  table(annotation$gene_name),
  decreasing = TRUE
)

head(
  gene_name_freq[gene_name_freq > 1],
  20
)

# ------------------------------------------------------------
# 10. Gene type distribution
# ------------------------------------------------------------

sort(
  table(annotation$gene_type),
  decreasing = TRUE
)

# ============================================================
# Step 14A-5B
# Compare expression-matrix Ensembl IDs with GENCODE v22
# ============================================================

expr_gene_ids <- read.delim(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  colClasses = c("character", rep("NULL", 513))
)

# ------------------------------------------------------------
# 1. Dimension
# ------------------------------------------------------------

dim(expr_gene_ids)

# ------------------------------------------------------------
# 2. Number of expression genes
# ------------------------------------------------------------

length(expr_gene_ids$Tag)

# ------------------------------------------------------------
# 3. Duplicate Ensembl IDs in expression matrix
# ------------------------------------------------------------

sum(
  duplicated(expr_gene_ids$Tag)
)

# ------------------------------------------------------------
# 4. Expression genes found in GENCODE annotation
# ------------------------------------------------------------

matched_ids <- intersect(
  expr_gene_ids$Tag,
  annotation$Ensembl_ID
)

length(matched_ids)

# ------------------------------------------------------------
# 5. Expression genes NOT found in annotation
# ------------------------------------------------------------

unmatched_ids <- setdiff(
  expr_gene_ids$Tag,
  annotation$Ensembl_ID
)

length(unmatched_ids)

# ------------------------------------------------------------
# 6. Mapping percentage
# ------------------------------------------------------------

100 *
  length(matched_ids) /
  length(expr_gene_ids$Tag)


# ============================================================
# Step 14A-6
# Inspect installed ESTIMATE package and functions
# ============================================================

# ------------------------------------------------------------
# 1. Check whether estimate package is installed
# ------------------------------------------------------------

"estimate" %in% rownames(installed.packages())

# ------------------------------------------------------------
# 2. Package version
# ------------------------------------------------------------

if ("estimate" %in% rownames(installed.packages())) {
  
  packageVersion("estimate")
  
}

# ------------------------------------------------------------
# 3. Load package if available
# ------------------------------------------------------------

if ("estimate" %in% rownames(installed.packages())) {
  
  library(estimate)
  
}

# ------------------------------------------------------------
# 4. Check available ESTIMATE functions
# ------------------------------------------------------------

if ("estimate" %in% rownames(installed.packages())) {
  
  ls("package:estimate")
  
}

# ------------------------------------------------------------
# 5. Inspect function arguments
# ------------------------------------------------------------

if ("estimate" %in% rownames(installed.packages())) {
  
  args(filterCommonGenes)
  args(estimateScore)
  
}

# ------------------------------------------------------------
# 6. Show function source
# ------------------------------------------------------------

if ("estimate" %in% rownames(installed.packages())) {
  
  print(filterCommonGenes)
  
  print(estimateScore)
  
}

# ============================================================
# Step 14B-1
# Prepare Ensembl -> Gene Symbol mapping
# ============================================================

# ------------------------------------------------------------
# 1. Create mapping table
# ------------------------------------------------------------

gene_map <- annotation[
  ,
  c(
    "Ensembl_ID",
    "gene_name",
    "gene_type"
  )
]

# ------------------------------------------------------------
# 2. Basic QC
# ------------------------------------------------------------

dim(gene_map)

sum(duplicated(gene_map$Ensembl_ID))

sum(
  is.na(gene_map$gene_name) |
    gene_map$gene_name == ""
)

# ------------------------------------------------------------
# 3. Number of unique gene symbols
# ------------------------------------------------------------

length(unique(gene_map$gene_name))

# ------------------------------------------------------------
# 4. Number of duplicated gene symbols
# ------------------------------------------------------------

sum(duplicated(gene_map$gene_name))

# ------------------------------------------------------------
# 5. Save mapping table
# ------------------------------------------------------------

write.table(
  gene_map,
  file = "Workflow-04-2C-1A_Ensembl_to_GeneSymbol_GENCODE_v22.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Step 14B-2
# Read complete tumor-only expression matrix
# ============================================================

expr_full <- read.delim(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

##Then check immediately：
# ------------------------------------------------------------
# 1. Dimension
# ------------------------------------------------------------

dim(expr_full)

# ------------------------------------------------------------
# 2. First column
# ------------------------------------------------------------

colnames(expr_full)[1]

# ------------------------------------------------------------
# 3. Sample identity
# ------------------------------------------------------------

identical(
  colnames(expr_full)[-1],
  rownames(datExpr_MAD)
)

# ------------------------------------------------------------
# 4. Gene identity
# ------------------------------------------------------------

identical(
  expr_full$Tag,
  annotation$Ensembl_ID
)

# ------------------------------------------------------------
# 5. Missing values
# ------------------------------------------------------------

sum(is.na(expr_full))

# ------------------------------------------------------------
# 6. Numeric expression check
# ------------------------------------------------------------

is.numeric(expr_full[[2]])

# ============================================================
# Step 14B-3
# Map Ensembl IDs to Gene Symbols
# ============================================================

# ------------------------------------------------------------
# 1. Extract expression matrix
# ------------------------------------------------------------

expr_mat <- as.matrix(
  expr_full[, -1]
)

rownames(expr_mat) <- expr_full$Tag

# ------------------------------------------------------------
# 2. Ensure numeric matrix
# ------------------------------------------------------------

storage.mode(expr_mat) <- "numeric"

# ------------------------------------------------------------
# 3. Match annotation to expression matrix
# ------------------------------------------------------------

map_idx <- match(
  rownames(expr_mat),
  annotation$Ensembl_ID
)

stopifnot(
  all(!is.na(map_idx))
)

gene_symbol <- annotation$gene_name[map_idx]
gene_type <- annotation$gene_type[map_idx]

# ------------------------------------------------------------
# 4. Check mapping
# ------------------------------------------------------------

sum(is.na(gene_symbol))

sum(gene_symbol == "")

length(gene_symbol)

length(unique(gene_symbol))

# ============================================================
# Step 14B-4
# Resolve duplicated Gene Symbols
# ============================================================

# ------------------------------------------------------------
# 1. Median expression for each Ensembl gene
# ------------------------------------------------------------

gene_median <- apply(
  expr_mat,
  1,
  median,
  na.rm = TRUE
)

# ------------------------------------------------------------
# 2. Mapping table for duplicate resolution
# ------------------------------------------------------------

mapping_for_selection <- data.frame(
  Ensembl_ID = rownames(expr_mat),
  GeneSymbol = gene_symbol,
  GeneType = gene_type,
  MedianExpression = gene_median,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 3. Order by GeneSymbol and descending median expression
# ------------------------------------------------------------

mapping_for_selection <- mapping_for_selection[
  order(
    mapping_for_selection$GeneSymbol,
    -mapping_for_selection$MedianExpression
  ),
]

# ------------------------------------------------------------
# 4. Keep the highest-median-expression Ensembl entry
#    for each Gene Symbol
# ------------------------------------------------------------

keep_idx <- !duplicated(
  mapping_for_selection$GeneSymbol
)

selected_mapping <- mapping_for_selection[
  keep_idx,
]

# ------------------------------------------------------------
# 5. QC
# ------------------------------------------------------------

dim(selected_mapping)

sum(
  duplicated(selected_mapping$GeneSymbol)
)

length(
  unique(selected_mapping$GeneSymbol)
)

# ============================================================
# Step 14B-5
# Create unique Gene Symbol expression matrix
# ============================================================

selected_ensembl <- selected_mapping$Ensembl_ID

expr_estimate <- expr_mat[
  selected_ensembl,
  ,
  drop = FALSE
]

# Replace Ensembl rownames with Gene Symbols

rownames(expr_estimate) <- selected_mapping$GeneSymbol

# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

dim(expr_estimate)

sum(
  duplicated(rownames(expr_estimate))
)

sum(is.na(expr_estimate))

# ------------------------------------------------------------
# Sample identity
# ------------------------------------------------------------

identical(
  colnames(expr_estimate),
  rownames(datExpr_MAD)
)

# ------------------------------------------------------------
# Preview
# ------------------------------------------------------------

expr_estimate[
  1:10,
  1:6
]

# ============================================================
# Step 14B-6
# Save ESTIMATE input matrix
# ============================================================

estimate_input_file <-
  "Workflow-04-2C-1A_ESTIMATE_input_GeneSymbol_log2_FPKM_UQ.txt"

write.table(
  cbind(
    GeneSymbol = rownames(expr_estimate),
    as.data.frame(expr_estimate)
  ),
  file = estimate_input_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

file.exists(estimate_input_file)

file.info(
  estimate_input_file
)$size

# ============================================================
# Step 14B-7
# Final QC before ESTIMATE scoring
# ============================================================

# Gene number
nrow(expr_estimate)

# Sample number
ncol(expr_estimate)

# Unique Gene Symbols
length(unique(rownames(expr_estimate)))

# Unique samples
length(unique(colnames(expr_estimate)))

# Missing values
sum(is.na(expr_estimate))

# Duplicate genes
sum(duplicated(rownames(expr_estimate)))

# Sample identity
identical(
  colnames(expr_estimate),
  rownames(datExpr_MAD)
)

# Expression range
summary(
  as.numeric(expr_estimate)
)


# ============================================================
# Step 14C-1
# Prepare ESTIMATE-compatible GCT input
# ============================================================

estimate_gct_input <-
  "Workflow-04-2C-1A_ESTIMATE_input_GCT.txt"

estimate_input_df <- data.frame(
  GeneSymbol = rownames(expr_estimate),
  Description = rownames(expr_estimate),
  expr_estimate,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.table(
  estimate_input_df,
  file = estimate_gct_input,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Check file
# ------------------------------------------------------------

file.exists(estimate_gct_input)

file.info(
  estimate_gct_input
)$size

##Then check the first few lines:
# ------------------------------------------------------------
# Preview
# ------------------------------------------------------------

read.delim(
  estimate_gct_input,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  nrows = 5
)[, 1:6]


# ============================================================
# Step 14C-2A
# Diagnose ESTIMATE package / namespace
# ============================================================

# ------------------------------------------------------------
# 1. Check package installation
# ------------------------------------------------------------

"estimate" %in% rownames(installed.packages())

# ------------------------------------------------------------
# 2. Show installed version
# ------------------------------------------------------------

packageVersion("estimate")

# ------------------------------------------------------------
# 3. Check whether package is currently attached
# ------------------------------------------------------------

"package:estimate" %in% search()

# ------------------------------------------------------------
# 4. Load estimate package explicitly
# ------------------------------------------------------------

library(estimate)

# ------------------------------------------------------------
# 5. Check whether functions exist in namespace
# ------------------------------------------------------------

exists(
  "filterCommonGenes",
  where = asNamespace("estimate"),
  inherits = FALSE
)

exists(
  "estimateScore",
  where = asNamespace("estimate"),
  inherits = FALSE
)

# ------------------------------------------------------------
# 6. Check exported functions
# ------------------------------------------------------------

ls("package:estimate")

# ============================================================
# Step 14C-2B
# Locate ESTIMATE functions
# ============================================================

find("filterCommonGenes")

find("estimateScore")

getAnywhere(filterCommonGenes)

getAnywhere(estimateScore)

exists(
  "filterCommonGenes",
  where = asNamespace("estimate"),
  inherits = FALSE
)

find("filterCommonGenes")


# ============================================================
# Step 14C-2C
# Check R / Bioconductor environment
# ============================================================

R.version.string

BiocManager::version()

"estimate" %in% rownames(installed.packages())

if ("estimate" %in% rownames(installed.packages())) {
  packageVersion("estimate")
} else {
  cat("ESTIMATE package is NOT installed.\n")
}

.libPaths()

sessionInfo()


# ============================================================
# Step 14C-2D
# Install ESTIMATE
# ============================================================

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "Moonerss/ESTIMATE"
)

# ============================================================
# Step 14C-2E
# Verify ESTIMATE installation
# ============================================================

library(ESTIMATE)

packageVersion("ESTIMATE")

find("filterCommonGenes")
find("estimateScore")
find("plotPurity")

exists(
  "filterCommonGenes",
  where = asNamespace("ESTIMATE"),
  inherits = FALSE
)

exists(
  "estimateScore",
  where = asNamespace("ESTIMATE"),
  inherits = FALSE
)

exists(
  "plotPurity",
  where = asNamespace("ESTIMATE"),
  inherits = FALSE
)

# ============================================================
# Step 14C-2F
# Inspect ESTIMATE 0.1.0 function interfaces
# ============================================================

packageVersion("ESTIMATE")

args(filterCommonGenes)

args(estimateScore)

args(plotPurity)

getAnywhere(filterCommonGenes)

getAnywhere(estimateScore)

# ============================================================
# Step 14C-3
# ESTIMATE common-gene filtering
# ============================================================

library(ESTIMATE)

# ------------------------------------------------------------
# 3.1 Confirm input matrix
# ------------------------------------------------------------

stopifnot(is.matrix(expr_estimate))
stopifnot(is.numeric(expr_estimate))

dim(expr_estimate)

sum(is.na(expr_estimate))

sum(duplicated(rownames(expr_estimate)))

sum(duplicated(colnames(expr_estimate)))

# ------------------------------------------------------------
# 3.2 Filter to ESTIMATE common genes
# ------------------------------------------------------------

estimate_common <-
  filterCommonGenes(
    expr = expr_estimate,
    id = "GeneSymbol"
  )

# ------------------------------------------------------------
# 3.3 Basic QC
# ------------------------------------------------------------

dim(estimate_common)

head(estimate_common[, 1:6])

sum(is.na(estimate_common))

sum(duplicated(estimate_common$NAME))

length(unique(estimate_common$NAME))

# ============================================================
# Step 14C-4A
# Calculate ESTIMATE TME scores
# ============================================================

estimate_scores <-
  estimateScore(
    merged_expr = estimate_common,
    platform = "illumina"
  )

# ------------------------------------------------------------
# Basic structure
# ------------------------------------------------------------

dim(estimate_scores)

rownames(estimate_scores)

colnames(estimate_scores)[1:10]

colnames(estimate_scores)[501:513]

# ------------------------------------------------------------
# Check missing values
# ------------------------------------------------------------

sum(is.na(estimate_scores))

# ------------------------------------------------------------
# Inspect first samples
# ------------------------------------------------------------

estimate_scores[, 1:6]

# ============================================================
# Step 14C-4B
# ESTIMATE score distribution
# ============================================================

summary(
  estimate_scores["StromalScore", ]
)

summary(
  estimate_scores["ImmuneScore", ]
)

summary(
  estimate_scores["ESTIMATEScore", ]
)

# ------------------------------------------------------------
# Score variability
# ------------------------------------------------------------

sd(
  estimate_scores["StromalScore", ]
)

sd(
  estimate_scores["ImmuneScore", ]
)

sd(
  estimate_scores["ESTIMATEScore", ]
)

# ============================================================
# Step 14C-4C
# Verify sample alignment
# ============================================================

estimate_sample_ids <-
  colnames(estimate_scores)

wgcna_sample_ids <-
  rownames(MEs)

length(estimate_sample_ids)

length(wgcna_sample_ids)

length(intersect(
  estimate_sample_ids,
  wgcna_sample_ids
))

identical(
  estimate_sample_ids,
  wgcna_sample_ids
)

head(estimate_sample_ids)

head(wgcna_sample_ids)

tail(estimate_sample_ids)

tail(wgcna_sample_ids)


# ============================================================
# Step 14C-4D
# Normalize ESTIMATE sample IDs
# ============================================================

# ------------------------------------------------------------
# 1. Preserve original ESTIMATE IDs
# ------------------------------------------------------------

estimate_sample_ids_raw <-
  colnames(estimate_scores)

# ------------------------------------------------------------
# 2. Convert "." back to "-" for TCGA IDs
# ------------------------------------------------------------

estimate_sample_ids_normalized <-
  gsub(
    "\\.",
    "-",
    estimate_sample_ids_raw
  )

# ------------------------------------------------------------
# 3. Check normalized IDs
# ------------------------------------------------------------

head(estimate_sample_ids_normalized)

tail(estimate_sample_ids_normalized)

# ------------------------------------------------------------
# 4. Check duplicates
# ------------------------------------------------------------

sum(
  duplicated(estimate_sample_ids_normalized)
)

# ------------------------------------------------------------
# 5. Compare with WGCNA samples
# ------------------------------------------------------------

length(
  intersect(
    estimate_sample_ids_normalized,
    wgcna_sample_ids
  )
)

length(
  setdiff(
    estimate_sample_ids_normalized,
    wgcna_sample_ids
  )
)

length(
  setdiff(
    wgcna_sample_ids,
    estimate_sample_ids_normalized
  )
)

# ============================================================
# Step 14C-4E
# Reorder ESTIMATE scores to WGCNA sample order
# ============================================================

# ------------------------------------------------------------
# Assign normalized TCGA sample IDs
# ------------------------------------------------------------

colnames(estimate_scores) <-
  estimate_sample_ids_normalized

# ------------------------------------------------------------
# Reorder columns according to WGCNA sample order
# ------------------------------------------------------------

estimate_scores_wgcna <-
  estimate_scores[
    ,
    wgcna_sample_ids,
    drop = FALSE
  ]

# ------------------------------------------------------------
# Final alignment QC
# ------------------------------------------------------------

dim(estimate_scores_wgcna)

identical(
  colnames(estimate_scores_wgcna),
  rownames(MEs)
)

sum(
  is.na(estimate_scores_wgcna)
)

sum(
  duplicated(colnames(estimate_scores_wgcna))
)

# ------------------------------------------------------------
# Preview
# ------------------------------------------------------------

head(colnames(estimate_scores_wgcna))

head(rownames(MEs))

tail(colnames(estimate_scores_wgcna))

tail(rownames(MEs))

##Step 14D：Tumor-only WGCNA module–TME association
# ============================================================
# Step 14D-1
# Confirm WGCNA module and TME structures
# ============================================================

# ------------------------------------------------------------
# WGCNA module eigengenes
# ------------------------------------------------------------

dim(MEs)

head(colnames(MEs))

tail(colnames(MEs))

# ------------------------------------------------------------
# ESTIMATE TME scores
# ------------------------------------------------------------

dim(estimate_scores_wgcna)

rownames(estimate_scores_wgcna)

# ------------------------------------------------------------
# Sample alignment
# ------------------------------------------------------------

identical(
  rownames(MEs),
  colnames(estimate_scores_wgcna)
)

# ------------------------------------------------------------
# Missingness
# ------------------------------------------------------------

sum(is.na(MEs))

sum(is.na(estimate_scores_wgcna))

# ============================================================
# Step 14D-2
# Tumor-only WGCNA module–TME association
# Spearman correlation + BH/FDR
# ============================================================

# ------------------------------------------------------------
# 1. Prepare objects
# ------------------------------------------------------------

module_names <-
  colnames(MEs)

tme_traits <-
  rownames(estimate_scores_wgcna)

# ------------------------------------------------------------
# 2. Initialize result container
# ------------------------------------------------------------

module_tme_results <-
  data.frame(
    Module = character(),
    Trait = character(),
    N = integer(),
    Spearman_rho = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )

# ------------------------------------------------------------
# 3. Spearman correlation for each module × TME trait
# ------------------------------------------------------------

for (module in module_names) {
  
  for (trait in tme_traits) {
    
    x <- MEs[, module]
    y <- as.numeric(
      estimate_scores_wgcna[trait, ]
    )
    
    complete_idx <-
      complete.cases(x, y)
    
    x_complete <-
      x[complete_idx]
    
    y_complete <-
      y[complete_idx]
    
    test <-
      cor.test(
        x_complete,
        y_complete,
        method = "spearman",
        exact = FALSE
      )
    
    module_tme_results <-
      rbind(
        module_tme_results,
        data.frame(
          Module = module,
          Trait = trait,
          N = length(x_complete),
          Spearman_rho = unname(
            test$estimate
          ),
          P_value = test$p.value,
          stringsAsFactors = FALSE
        )
      )
  }
}

# ------------------------------------------------------------
# 4. BH/FDR correction across all 150 tests
# ------------------------------------------------------------

module_tme_results$FDR <-
  p.adjust(
    module_tme_results$P_value,
    method = "BH"
  )

# ------------------------------------------------------------
# 5. Sort by absolute correlation
# ------------------------------------------------------------

module_tme_results <-
  module_tme_results[
    order(
      -abs(module_tme_results$Spearman_rho)
    ),
  ]

rownames(module_tme_results) <- NULL

# ------------------------------------------------------------
# 6. Basic QC
# ------------------------------------------------------------

dim(module_tme_results)

sum(is.na(module_tme_results$Spearman_rho))

sum(is.na(module_tme_results$P_value))

sum(is.na(module_tme_results$FDR))

head(
  module_tme_results,
  20
)

# ============================================================
# Step 14D-3
# Significant module–TME associations
# ============================================================

module_tme_sig <-
  subset(
    module_tme_results,
    FDR < 0.05
  )

dim(module_tme_sig)

module_tme_sig[
  order(
    module_tme_sig$FDR
  ),
]

# ------------------------------------------------------------
# Significant associations by TME trait
# ------------------------------------------------------------

table(
  module_tme_sig$Trait
)

# ------------------------------------------------------------
# Top associations for each TME trait
# ------------------------------------------------------------

lapply(
  split(
    module_tme_results,
    module_tme_results$Trait
  ),
  function(df) {
    df[
      order(df$FDR),
      c(
        "Module",
        "N",
        "Spearman_rho",
        "P_value",
        "FDR"
      )
    ][1:min(10, nrow(df)), ]
  }
)

# ============================================================
# Step 14E-1
# Convert module–TME results to a 50 × 3 correlation matrix
# ============================================================

# ------------------------------------------------------------
# 1. Correlation matrix
# ------------------------------------------------------------

module_tme_rho_matrix <-
  reshape(
    module_tme_results[
      c(
        "Module",
        "Trait",
        "Spearman_rho"
      )
    ],
    idvar = "Module",
    timevar = "Trait",
    direction = "wide"
  )

# ------------------------------------------------------------
# 2. Clean column names
# ------------------------------------------------------------

colnames(module_tme_rho_matrix) <-
  sub(
    "^Spearman_rho\\.",
    "",
    colnames(module_tme_rho_matrix)
  )

# ------------------------------------------------------------
# 3. Restore Module as row names
# ------------------------------------------------------------

rownames(module_tme_rho_matrix) <-
  module_tme_rho_matrix$Module

module_tme_rho_matrix$Module <-
  NULL

# ------------------------------------------------------------
# 4. Arrange traits in a fixed order
# ------------------------------------------------------------

module_tme_rho_matrix <-
  module_tme_rho_matrix[
    ,
    c(
      "ImmuneScore",
      "StromalScore",
      "ESTIMATEScore"
    ),
    drop = FALSE
  ]

# ------------------------------------------------------------
# 5. QC
# ------------------------------------------------------------

dim(module_tme_rho_matrix)

head(
  module_tme_rho_matrix,
  10
)

module_tme_rho_matrix[
  "ME_brown",
  ,
  drop = FALSE
]

sum(
  is.na(module_tme_rho_matrix)
)

# ============================================================
# Step 14E-2
# Rank modules by their overall TME association
# ============================================================

module_tme_strength <-
  data.frame(
    Module = rownames(module_tme_rho_matrix),
    ImmuneScore =
      module_tme_rho_matrix[, "ImmuneScore"],
    StromalScore =
      module_tme_rho_matrix[, "StromalScore"],
    ESTIMATEScore =
      module_tme_rho_matrix[, "ESTIMATEScore"],
    stringsAsFactors = FALSE
  )

module_tme_strength$Mean_abs_rho <-
  rowMeans(
    abs(
      module_tme_strength[
        ,
        c(
          "ImmuneScore",
          "StromalScore",
          "ESTIMATEScore"
        )
      ]
    )
  )

module_tme_strength <-
  module_tme_strength[
    order(
      -module_tme_strength$Mean_abs_rho
    ),
  ]

rownames(module_tme_strength) <- NULL

head(
  module_tme_strength,
  15
)

# ============================================================
# Step 14E-3
# Identify the actual WGCNA module name corresponding to brown
# ============================================================

# 1. Inspect all ME column names
colnames(MEs)

# 2. Search for "brown" in module names
grep(
  "brown",
  colnames(MEs),
  ignore.case = TRUE,
  value = TRUE
)

# 3. Inspect moduleColors
sort(
  unique(moduleColors)
)

# 4. Confirm the number of genes in each module
table(moduleColors)["brown"]

# 5. Show all ME columns containing brown, if any
brown_ME_candidates <-
  grep(
    "brown",
    colnames(MEs),
    ignore.case = TRUE,
    value = TRUE
  )

brown_ME_candidates

# ============================================================
# Step 14E-4
# Map WGCNA module colors to the corresponding ME columns
# ============================================================

# ------------------------------------------------------------
# 1. Check sample alignment
# ------------------------------------------------------------

identical(
  rownames(datExpr_MAD),
  rownames(MEs)
)

# ------------------------------------------------------------
# 2. Get all module colors
# ------------------------------------------------------------

module_color_names <-
  sort(
    unique(moduleColors)
  )

length(module_color_names)

# ------------------------------------------------------------
# 3. Calculate eigengene for each module color
# ------------------------------------------------------------

ME_color_list <- list()

for (color in module_color_names) {
  
  gene_idx <-
    which(moduleColors == color)
  
  expr_module <-
    datExpr_MAD[
      ,
      gene_idx,
      drop = FALSE
    ]
  
  ME_color_list[[color]] <-
    rowMeans(
      scale(expr_module),
      na.rm = TRUE
    )
}

ME_color_df <-
  as.data.frame(
    ME_color_list,
    check.names = FALSE
  )

rownames(ME_color_df) <-
  rownames(datExpr_MAD)

# ------------------------------------------------------------
# 4. QC
# ------------------------------------------------------------

dim(ME_color_df)

sum(is.na(ME_color_df))

head(
  colnames(ME_color_df)
)

# ------------------------------------------------------------
# 5. Correlate color-derived MEs with existing MEs
# ------------------------------------------------------------

ME_mapping_results <-
  data.frame(
    ModuleColor = character(),
    ME_column = character(),
    Correlation = numeric(),
    stringsAsFactors = FALSE
  )

for (color in module_color_names) {
  
  for (ME_name in colnames(MEs)) {
    
    r <-
      cor(
        ME_color_df[, color],
        MEs[, ME_name],
        method = "pearson",
        use = "complete.obs"
      )
    
    ME_mapping_results <-
      rbind(
        ME_mapping_results,
        data.frame(
          ModuleColor = color,
          ME_column = ME_name,
          Correlation = r,
          stringsAsFactors = FALSE
        )
      )
  }
}

# ------------------------------------------------------------
# 6. Find the best ME match for each module color
# ------------------------------------------------------------

ME_mapping_best <-
  do.call(
    rbind,
    lapply(
      split(
        ME_mapping_results,
        ME_mapping_results$ModuleColor
      ),
      function(df) {
        
        df[
          which.max(
            abs(df$Correlation)
          ),
          ,
          drop = FALSE
        ]
      }
    )
  )

rownames(ME_mapping_best) <- NULL

# ------------------------------------------------------------
# 7. Sort by module color
# ------------------------------------------------------------

ME_mapping_best <-
  ME_mapping_best[
    order(
      ME_mapping_best$ModuleColor
    ),
  ]

# ------------------------------------------------------------
# 8. Inspect the mapping
# ------------------------------------------------------------

ME_mapping_best

# ============================================================
# Step 14E-4A
# Inspect original WGCNA eigengene object
# ============================================================

colnames(net$MEs)

dim(net$MEs)

head(
  colnames(net$MEs),
  20
)

grep(
  "brown",
  colnames(net$MEs),
  ignore.case = TRUE,
  value = TRUE
)

grep(
  "blue",
  colnames(net$MEs),
  ignore.case = TRUE,
  value = TRUE
)

# ============================================================
# Step 14E-4B
# Reconstruct formal WGCNA module eigengenes by module color
# ============================================================
# ============================================================
# Step 14E-4B-1
# Check WGCNA availability
# ============================================================

find("moduleEigengenes")

packageVersion("WGCNA")

sessionInfo()

# ------------------------------------------------------------
# 1. Calculate formal module eigengenes
# ------------------------------------------------------------

ME_color_formal_raw <-
  moduleEigengenes(
    datExpr_MAD,
    colors = moduleColors
  )$eigengenes


# ------------------------------------------------------------
# 2. Inspect the resulting eigengene names
# ------------------------------------------------------------

dim(ME_color_formal_raw)

colnames(ME_color_formal_raw)

grep(
  "brown",
  colnames(ME_color_formal_raw),
  ignore.case = TRUE,
  value = TRUE
)

grep(
  "blue",
  colnames(ME_color_formal_raw),
  ignore.case = TRUE,
  value = TRUE
)


# ------------------------------------------------------------
# 3. Sample alignment
# ------------------------------------------------------------

identical(
  rownames(ME_color_formal_raw),
  rownames(MEs)
)


# ------------------------------------------------------------
# 4. Check missing values
# ------------------------------------------------------------

sum(
  is.na(ME_color_formal_raw)
)

##Once it has run successfully, please proceed with the following QC.
dim(ME_color_formal_raw)

colnames(ME_color_formal_raw)

grep(
  "brown",
  colnames(ME_color_formal_raw),
  ignore.case = TRUE,
  value = TRUE
)

grep(
  "blue",
  colnames(ME_color_formal_raw),
  ignore.case = TRUE,
  value = TRUE
)

identical(
  rownames(ME_color_formal_raw),
  rownames(MEs)
)

sum(
  is.na(ME_color_formal_raw)
)


# ============================================================
# Step 14E-5A
# Trace the brown-module variable used in the Cox workflow
# ============================================================

# ------------------------------------------------------------
# 1. Inspect cox_design_complete
# ------------------------------------------------------------

dim(cox_design_complete)

colnames(cox_design_complete)

grep(
  "brown",
  colnames(cox_design_complete),
  ignore.case = TRUE,
  value = TRUE
)


# ------------------------------------------------------------
# 2. Inspect stage_complete_data
# ------------------------------------------------------------

dim(stage_complete_data)

colnames(stage_complete_data)

grep(
  "brown",
  colnames(stage_complete_data),
  ignore.case = TRUE,
  value = TRUE
)


# ------------------------------------------------------------
# 3. Inspect final model data variables
# ------------------------------------------------------------

all.vars(
  formula(final_td_model)
)

# ============================================================
# Step 14E-5B
# Match Cox ME_brown with the formal WGCNA MEbrown
# ============================================================

# ------------------------------------------------------------
# 1. Build a sample-indexed formal brown eigengene
# ------------------------------------------------------------

MEbrown_formal_df <-
  data.frame(
    Sample = rownames(ME_color_formal_raw),
    MEbrown_formal = MEbrown_formal,
    stringsAsFactors = FALSE
  )


# ------------------------------------------------------------
# 2. Match formal MEbrown to cox_design_complete
# ------------------------------------------------------------

brown_validation <-
  merge(
    cox_design_complete[
      ,
      c("Sample", "ME_brown"),
      drop = FALSE
    ],
    MEbrown_formal_df,
    by = "Sample",
    all = FALSE,
    sort = FALSE
  )


# ------------------------------------------------------------
# 3. Basic QC
# ------------------------------------------------------------

dim(brown_validation)

sum(
  is.na(brown_validation$ME_brown)
)

sum(
  is.na(brown_validation$MEbrown_formal)
)

length(
  unique(brown_validation$Sample)
)


# ------------------------------------------------------------
# 4. Check sample identity
# ------------------------------------------------------------

setequal(
  brown_validation$Sample,
  cox_design_complete$Sample
)


setequal(
  brown_validation$Sample,
  rownames(ME_color_formal_raw)[
    rownames(ME_color_formal_raw) %in%
      cox_design_complete$Sample
  ]
)


# ------------------------------------------------------------
# 5. Correlation
# ------------------------------------------------------------

cor_brown <-
  cor(
    brown_validation$ME_brown,
    brown_validation$MEbrown_formal,
    method = "pearson",
    use = "complete.obs"
  )

cor_brown


# ------------------------------------------------------------
# 6. Linear relationship
# ------------------------------------------------------------

summary(
  lm(
    ME_brown ~ MEbrown_formal,
    data = brown_validation
  )
)


# ------------------------------------------------------------
# 7. Compare means and SDs
# ------------------------------------------------------------

rbind(
  Cox_ME_brown =
    c(
      Mean = mean(
        brown_validation$ME_brown
      ),
      SD = sd(
        brown_validation$ME_brown
      )
    ),
  
  Formal_MEbrown =
    c(
      Mean = mean(
        brown_validation$MEbrown_formal
      ),
      SD = sd(
        brown_validation$MEbrown_formal
      )
    )
)

# ============================================================
# Step 14E-6
# Direct verification of brown module–TME associations
# ============================================================

# ------------------------------------------------------------
# 1. Prepare brown module eigengene
# ------------------------------------------------------------

MEbrown_formal <-
  ME_color_formal_raw[, "MEbrown"]


# ------------------------------------------------------------
# 2. Check sample alignment
# ------------------------------------------------------------

identical(
  names(MEbrown_formal),
  colnames(estimate_scores_wgcna)
)

identical(
  names(MEbrown_formal),
  rownames(MEs)
)


# ------------------------------------------------------------
# 3. Calculate Spearman correlations
# ------------------------------------------------------------

brown_tme_results <-
  data.frame(
    Module = "brown",
    Trait = rownames(estimate_scores_wgcna),
    Spearman_rho = NA_real_,
    P_value = NA_real_,
    N = NA_integer_,
    stringsAsFactors = FALSE
  )


for (i in seq_len(nrow(estimate_scores_wgcna))) {
  
  trait_values <-
    as.numeric(
      estimate_scores_wgcna[i, ]
    )
  
  test <-
    cor.test(
      MEbrown_formal,
      trait_values,
      method = "spearman",
      exact = FALSE
    )
  
  brown_tme_results$Spearman_rho[i] <-
    unname(test$estimate)
  
  brown_tme_results$P_value[i] <-
    test$p.value
  
  brown_tme_results$N[i] <-
    sum(
      complete.cases(
        MEbrown_formal,
        trait_values
      )
    )
}


# ------------------------------------------------------------
# 4. Multiple-testing correction
# ------------------------------------------------------------
#
# Important:
# The primary 150-test analysis already used BH correction
# across all 50 modules × 3 TME traits.
#
# Therefore, we do NOT treat these three brown-module tests
# as a new independent 3-test discovery set.
#
# We will later retrieve the corresponding FDR values
# from the locked 150-test result.
# ------------------------------------------------------------


# ------------------------------------------------------------
# 5. Inspect results
# ------------------------------------------------------------

brown_tme_results


# ============================================================
# Step 14E-7
# Formal sample alignment and brown-module label mapping
# ============================================================

# ------------------------------------------------------------
# 1. Recover formal sample IDs of MEbrown
# ------------------------------------------------------------

brown_sample_ids <-
  rownames(ME_color_formal_raw)

estimate_sample_ids <-
  colnames(estimate_scores_wgcna)

MEs_sample_ids <-
  rownames(MEs)


# ------------------------------------------------------------
# 2. Check dimensions and sample IDs
# ------------------------------------------------------------

cat("MEbrown sample number:",
    length(brown_sample_ids), "\n")

cat("ESTIMATE sample number:",
    length(estimate_sample_ids), "\n")

cat("MEs sample number:",
    length(MEs_sample_ids), "\n")


cat(
  "Brown vs ESTIMATE setequal:",
  setequal(
    brown_sample_ids,
    estimate_sample_ids
  ),
  "\n"
)

cat(
  "Brown vs MEs setequal:",
  setequal(
    brown_sample_ids,
    MEs_sample_ids
  ),
  "\n"
)


# ------------------------------------------------------------
# 3. Check exact order
# ------------------------------------------------------------

cat(
  "Brown vs ESTIMATE identical:",
  identical(
    brown_sample_ids,
    estimate_sample_ids
  ),
  "\n"
)

cat(
  "Brown vs MEs identical:",
  identical(
    brown_sample_ids,
    MEs_sample_ids
  ),
  "\n"
)


# ------------------------------------------------------------
# 4. Check whether reordering is required
# ------------------------------------------------------------

brown_order_in_estimate <-
  match(
    brown_sample_ids,
    estimate_sample_ids
  )

brown_order_in_MEs <-
  match(
    brown_sample_ids,
    MEs_sample_ids
  )

cat(
  "All Brown IDs found in ESTIMATE:",
  all(!is.na(brown_order_in_estimate)),
  "\n"
)

cat(
  "All Brown IDs found in MEs:",
  all(!is.na(brown_order_in_MEs)),
  "\n"
)


# ------------------------------------------------------------
# 5. Formal Brown eigengene vs existing MEs
# ------------------------------------------------------------

brown_vs_MEs_cor <-
  sapply(
    colnames(MEs),
    function(x) {
      
      cor(
        MEbrown_formal,
        MEs[, x],
        method = "pearson"
      )
      
    }
  )

brown_vs_MEs_cor_sorted <-
  sort(
    abs(brown_vs_MEs_cor),
    decreasing = TRUE
  )

brown_vs_MEs_cor_sorted[1:10]


# ------------------------------------------------------------
# 6. Identify the exact MExx corresponding to brown
# ------------------------------------------------------------

brown_match_MExx <-
  names(
    brown_vs_MEs_cor_sorted[
      brown_vs_MEs_cor_sorted > 0.999999
    ]
  )

brown_match_MExx


# ------------------------------------------------------------
# 7. Display the exact correlation
# ------------------------------------------------------------

brown_vs_MEs_cor[
  brown_match_MExx
]

# ============================================================
# Step 14E-8
# Extract the locked 150-test FDR for the brown module
# ============================================================

# ------------------------------------------------------------
# 1. Inspect the structure of the locked 150-test result
# ------------------------------------------------------------

dim(module_tme_results)

colnames(module_tme_results)

head(module_tme_results)


# ------------------------------------------------------------
# 2. Identify the ME3 rows
# ------------------------------------------------------------

brown_locked_results <-
  module_tme_results[
    module_tme_results$Module == "ME3",
    ,
    drop = FALSE
  ]

brown_locked_results


# ------------------------------------------------------------
# 3. Compare with direct formal-brown results
# ------------------------------------------------------------

brown_tme_results


# ------------------------------------------------------------
# 4. Check that rho and P values agree
# ------------------------------------------------------------

comparison_brown <-
  merge(
    brown_locked_results,
    brown_tme_results,
    by = "Trait",
    suffixes = c("_locked150", "_formal_brown"),
    all = TRUE,
    sort = FALSE
  )

comparison_brown


# ------------------------------------------------------------
# 5. Quantify numerical differences
# ------------------------------------------------------------

comparison_brown$rho_difference <-
  comparison_brown$Spearman_rho_locked150 -
  comparison_brown$Spearman_rho_formal_brown

comparison_brown$P_difference <-
  comparison_brown$P_value_locked150 -
  comparison_brown$P_value_formal_brown

comparison_brown

# ============================================================
# Step 14E-9
# Final module–TME association summary
# ============================================================

# ------------------------------------------------------------
# 1. Add absolute correlation strength
# ------------------------------------------------------------

module_tme_results$Abs_rho <-
  abs(module_tme_results$Spearman_rho)


# ------------------------------------------------------------
# 2. Statistical significance
# ------------------------------------------------------------

module_tme_results$FDR_significant <-
  module_tme_results$FDR < 0.05


# ------------------------------------------------------------
# 3. Direction
# ------------------------------------------------------------

module_tme_results$Direction <-
  ifelse(
    module_tme_results$Spearman_rho > 0,
    "Positive",
    "Negative"
  )


# ------------------------------------------------------------
# 4. Significant associations by trait
# ------------------------------------------------------------

tme_significant_summary <-
  aggregate(
    FDR_significant ~ Trait,
    data = module_tme_results,
    FUN = sum
  )

tme_significant_summary


# ------------------------------------------------------------
# 5. Strongest associations for each TME trait
# ------------------------------------------------------------

tme_top10_by_trait <-
  do.call(
    rbind,
    lapply(
      unique(module_tme_results$Trait),
      function(trait_name) {
        
        tmp <-
          module_tme_results[
            module_tme_results$Trait == trait_name,
            ,
            drop = FALSE
          ]
        
        tmp <-
          tmp[
            order(
              -tmp$Abs_rho
            ),
            ,
            drop = FALSE
          ]
        
        head(tmp, 10)
        
      }
    )
  )

tme_top10_by_trait


# ------------------------------------------------------------
# 6. Modules with strongest overall TME associations
# ------------------------------------------------------------

module_tme_overall <-
  aggregate(
    Abs_rho ~ Module,
    data = module_tme_results,
    FUN = mean
  )

module_tme_overall <-
  module_tme_overall[
    order(
      -module_tme_overall$Abs_rho
    ),
    ,
    drop = FALSE
  ]

head(
  module_tme_overall,
  20
)


# ------------------------------------------------------------
# 7. Identify modules significant for all three TME traits
# ------------------------------------------------------------

module_all3 <-
  unique(
    module_tme_results$Module[
      ave(
        module_tme_results$FDR < 0.05,
        module_tme_results$Module,
        FUN = sum
      ) == 3
    ]
  )

module_all3

# ============================================================
# Step 14E-10
# Formal mapping between MExx labels and WGCNA module colors
# ============================================================

# ------------------------------------------------------------
# 1. Calculate Pearson correlation between every formal
#    color-based eigengene and every existing MExx eigengene
# ------------------------------------------------------------

ME_color_names <-
  colnames(ME_color_formal_raw)

ME_internal_names <-
  colnames(MEs)


ME_mapping_cor <-
  matrix(
    NA_real_,
    nrow = length(ME_color_names),
    ncol = length(ME_internal_names),
    dimnames = list(
      ME_color_names,
      ME_internal_names
    )
  )


for (i in seq_along(ME_color_names)) {
  
  for (j in seq_along(ME_internal_names)) {
    
    ME_mapping_cor[i, j] <-
      cor(
        ME_color_formal_raw[, i],
        MEs[, j],
        method = "pearson"
      )
  }
}


# ------------------------------------------------------------
# 2. Find the best internal MExx for each formal color
# ------------------------------------------------------------

ME_color_to_internal <-
  data.frame(
    Color_ME = ME_color_names,
    Internal_ME = NA_character_,
    Pearson_r = NA_real_,
    stringsAsFactors = FALSE
  )


for (i in seq_len(nrow(ME_color_to_internal))) {
  
  cor_values <-
    ME_mapping_cor[i, ]
  
  best_index <-
    which.max(
      abs(cor_values)
    )
  
  ME_color_to_internal$Internal_ME[i] <-
    names(cor_values)[best_index]
  
  ME_color_to_internal$Pearson_r[i] <-
    cor_values[best_index]
}


# ------------------------------------------------------------
# 3. Inspect mapping quality
# ------------------------------------------------------------

ME_color_to_internal[
  order(
    ME_color_to_internal$Color_ME
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 4. Check whether all mappings are essentially exact
# ------------------------------------------------------------

summary(
  ME_color_to_internal$Pearson_r
)

sum(
  abs(
    ME_color_to_internal$Pearson_r
  ) < 0.999999
)


# ------------------------------------------------------------
# 5. Check uniqueness of internal MExx labels
# ------------------------------------------------------------

length(
  unique(
    ME_color_to_internal$Internal_ME
  )
)

nrow(
  ME_color_to_internal
)


# ------------------------------------------------------------
# 6. Display the complete mapping
# ------------------------------------------------------------

ME_color_to_internal


############################################################
# WORKFLOW-03
# Step 14A–14E FINAL LOCKED SAVE
#
# Paper 1: TCGA-LUAD tumor-only WGCNA revision
#
# Purpose:
#   Save the completed and validated Step 14A–14E analysis
#   as a reproducible checkpoint.
#
# Locked analysis:
#   Step 14A  Expression source + annotation QC
#   Step 14B  Ensembl -> gene-symbol expression matrix
#   Step 14C  ESTIMATE scoring
#   Step 14D  Module-TME correlations
#   Step 14E  Module color / ME mapping + Brown validation
#
# IMPORTANT:
#   Do NOT modify or recompute the locked results unless
#   a methodological correction is explicitly required.
############################################################


############################################################
# 0. BASIC SETTINGS
############################################################

options(stringsAsFactors = FALSE)

analysis_name <- "WORKFLOW-03_Step14A-14E_FINAL_LOCKED"

save_date <- format(Sys.time(), "%Y-%m-%d_%H%M%S")


############################################################
# 1. CREATE FINAL SAVE DIRECTORIES
############################################################

# Modify this path if your WORKFLOW-03 root directory is different

step14_root <- file.path(
  getwd(),
  "WORKFLOW-03",
  "Step14"
)

final_root <- file.path(
  step14_root,
  "Step14_FINAL_LOCKED"
)

dir.create(final_root, recursive = TRUE, showWarnings = FALSE)

dir.create(
  file.path(final_root, "RData"),
  showWarnings = FALSE
)

dir.create(
  file.path(final_root, "Tables"),
  showWarnings = FALSE
)

dir.create(
  file.path(final_root, "QC"),
  showWarnings = FALSE
)

dir.create(
  file.path(final_root, "Scripts"),
  showWarnings = FALSE
)

dir.create(
  file.path(final_root, "SessionInfo"),
  showWarnings = FALSE
)


############################################################
# 2. SAVE COMPLETE R WORKSPACE
############################################################

# Save all current objects.
#
# This is the most important recovery file.

save(
  list = ls(envir = .GlobalEnv),
  file = file.path(
    final_root,
    "RData",
    "WORKFLOW-03_Step14A-14E_FINAL_LOCKED.RData"
  )
)


############################################################
# 3. SAVE KEY ANALYSIS OBJECTS AS A SMALLER CHECKPOINT
############################################################

# These are the core objects required for Step 14 continuation.

key_objects <- c(
  "datExpr_MAD",
  "MEs",
  "moduleColors",
  "moduleSizes",
  "net",
  "clinical_aligned",
  "expr_estimate",
  "estimate_common",
  "estimate_scores",
  "estimate_scores_wgcna",
  "ME_color_formal_raw",
  "MEbrown_formal",
  "module_tme_results"
)

# Keep only objects that actually exist

key_objects <- key_objects[
  key_objects %in% ls(envir = .GlobalEnv)
]

save(
  list = key_objects,
  file = file.path(
    final_root,
    "RData",
    "WORKFLOW-03_Step14_KEY_OBJECTS.RData"
  )
)


############################################################
# 4. SAVE MODULE-TME RESULTS
############################################################

if (exists("module_tme_results")) {
  
  write.csv(
    module_tme_results,
    file = file.path(
      final_root,
      "Tables",
      "Step14D_module_TME_all_150_tests.csv"
    ),
    row.names = FALSE
  )
}


############################################################
# 5. SAVE SIGNIFICANT MODULE-TME ASSOCIATIONS
############################################################

if (exists("module_tme_results")) {
  
  module_tme_sig <- subset(
    module_tme_results,
    FDR < 0.05
  )
  
  write.csv(
    module_tme_sig,
    file = file.path(
      final_root,
      "Tables",
      "Step14D_module_TME_FDR05.csv"
    ),
    row.names = FALSE
  )
}


############################################################
# 6. SAVE BROWN MODULE VALIDATION
############################################################

if (exists("brown_validation")) {
  
  write.csv(
    brown_validation,
    file = file.path(
      final_root,
      "Tables",
      "Step14E_Brown_module_validation.csv"
    ),
    row.names = FALSE
  )
}


############################################################
# 7. SAVE FORMAL MODULE COLOR -> INTERNAL ME MAPPING
############################################################

# Reconstruct the formal mapping if the mapping object exists.
#
# Expected object:
#   module_color_to_ME_mapping
#
# If your actual mapping object has a different name,
# replace it here with that object.

mapping_candidates <- c(
  "module_color_to_ME_mapping",
  "color_ME_mapping",
  "module_mapping"
)

mapping_name <- mapping_candidates[
  mapping_candidates %in% ls(envir = .GlobalEnv)
]

if (length(mapping_name) > 0) {
  
  mapping_object <- get(mapping_name[1])
  
  write.csv(
    mapping_object,
    file = file.path(
      final_root,
      "Tables",
      "Step14E_module_color_to_MExx_mapping.csv"
    ),
    row.names = FALSE
  )
}


############################################################
# 8. CREATE A DEFINITIVE COLOR -> MExx TABLE
############################################################

# The following mapping was formally validated:
# Pearson r = 1 for all 50 modules.

formal_mapping <- data.frame(
  ModuleColor = c(
    "bisque4",
    "black",
    "blue",
    "brown",
    "brown4",
    "cyan",
    "darkgreen",
    "darkgrey",
    "darkmagenta",
    "darkolivegreen",
    "darkorange",
    "darkorange2",
    "darkred",
    "darkslateblue",
    "darkturquoise",
    "floralwhite",
    "green",
    "greenyellow",
    "grey",
    "grey60",
    "ivory",
    "lightcyan",
    "lightcyan1",
    "lightgreen",
    "lightsteelblue1",
    "lightyellow",
    "magenta",
    "mediumpurple3",
    "midnightblue",
    "orange",
    "orangered4",
    "paleturquoise",
    "pink",
    "plum1",
    "plum2",
    "purple",
    "red",
    "royalblue",
    "saddlebrown",
    "salmon",
    "sienna3",
    "skyblue",
    "skyblue3",
    "steelblue",
    "tan",
    "turquoise",
    "violet",
    "white",
    "yellow",
    "yellowgreen"
  ),
  Internal_ME = c(
    "ME47",
    "ME7",
    "ME2",
    "ME3",
    "ME46",
    "ME14",
    "ME22",
    "ME24",
    "ME34",
    "ME33",
    "ME26",
    "ME45",
    "ME21",
    "ME48",
    "ME23",
    "ME44",
    "ME5",
    "ME11",
    "ME0",
    "ME17",
    "ME43",
    "ME16",
    "ME42",
    "ME18",
    "ME41",
    "ME19",
    "ME9",
    "ME40",
    "ME15",
    "ME25",
    "ME39",
    "ME31",
    "ME8",
    "ME38",
    "ME49",
    "ME10",
    "ME6",
    "ME20",
    "ME29",
    "ME13",
    "ME35",
    "ME28",
    "ME37",
    "ME30",
    "ME12",
    "ME1",
    "ME32",
    "ME27",
    "ME4",
    "ME36"
  ),
  Pearson_r = 1
)

write.csv(
  formal_mapping,
  file = file.path(
    final_root,
    "Tables",
    "Step14E_DEFINITIVE_module_color_to_MExx_mapping.csv"
  ),
  row.names = FALSE
)


############################################################
# 9. SAVE BROWN-TME RESULTS
############################################################

if (exists("module_tme_results")) {
  
  brown_tme <- subset(
    module_tme_results,
    Module == "ME3"
  )
  
  # Add manuscript-friendly module name
  
  brown_tme$ModuleColor <- "brown"
  
  brown_tme <- brown_tme[
    c(
      "ModuleColor",
      "Module",
      "Trait",
      "N",
      "Spearman_rho",
      "P_value",
      "FDR"
    )
  ]
  
  write.csv(
    brown_tme,
    file = file.path(
      final_root,
      "Tables",
      "Step14E_Brown_ME3_TME_locked.csv"
    ),
    row.names = FALSE
  )
}


############################################################
# 10. SAVE SUMMARY OF SIGNIFICANT MODULES
############################################################

if (exists("module_tme_results")) {
  
  sig_summary <- aggregate(
    FDR < 0.05 ~ Trait,
    data = module_tme_results,
    FUN = sum
  )
  
  names(sig_summary) <- c(
    "Trait",
    "Number_FDR05"
  )
  
  write.csv(
    sig_summary,
    file = file.path(
      final_root,
      "Tables",
      "Step14D_TME_significance_summary.csv"
    ),
    row.names = FALSE
  )
}


############################################################
# 11. SAVE MODULES SIGNIFICANT FOR ALL THREE TME TRAITS
############################################################

if (exists("module_tme_results")) {
  
  sig_all3 <- subset(
    module_tme_results,
    FDR < 0.05
  )
  
  all3_modules <- unique(
    sig_all3$Module[
      ave(
        sig_all3$Module,
        sig_all3$Module,
        FUN = length
      ) == 3
    ]
  )
  
  all3_table <- data.frame(
    Module = sort(all3_modules)
  )
  
  write.csv(
    all3_table,
    file = file.path(
      final_root,
      "Tables",
      "Step14D_modules_significant_all_three_TME_traits.csv"
    ),
    row.names = FALSE
  )
}


############################################################
# 12. SAVE SAMPLE ALIGNMENT QC
############################################################

qc_alignment <- data.frame(
  Check = c(
    "WGCNA samples",
    "ESTIMATE samples",
    "ME formal samples",
    "WGCNA vs ESTIMATE identical",
    "WGCNA vs ME formal identical",
    "ESTIMATE NA",
    "Module-TME result rows"
  ),
  Value = c(
    nrow(MEs),
    ncol(estimate_scores_wgcna),
    nrow(ME_color_formal_raw),
    identical(
      rownames(MEs),
      colnames(estimate_scores_wgcna)
    ),
    identical(
      rownames(MEs),
      rownames(ME_color_formal_raw)
    ),
    sum(is.na(estimate_scores_wgcna)),
    ifelse(
      exists("module_tme_results"),
      nrow(module_tme_results),
      NA
    )
  )
)

write.csv(
  qc_alignment,
  file = file.path(
    final_root,
    "QC",
    "Step14_sample_alignment_QC.csv"
  ),
  row.names = FALSE
)


############################################################
# 13. SAVE BASIC DATASET DIMENSIONS
############################################################

dimension_summary <- data.frame(
  Object = c(
    "datExpr_MAD",
    "MEs",
    "expr_estimate",
    "estimate_common",
    "estimate_scores_wgcna",
    "ME_color_formal_raw",
    "module_tme_results"
  ),
  Rows = c(
    nrow(datExpr_MAD),
    nrow(MEs),
    nrow(expr_estimate),
    nrow(estimate_common),
    nrow(estimate_scores_wgcna),
    nrow(ME_color_formal_raw),
    ifelse(
      exists("module_tme_results"),
      nrow(module_tme_results),
      NA
    )
  ),
  Columns = c(
    ncol(datExpr_MAD),
    ncol(MEs),
    ncol(expr_estimate),
    ncol(estimate_common),
    ncol(estimate_scores_wgcna),
    ncol(ME_color_formal_raw),
    ifelse(
      exists("module_tme_results"),
      ncol(module_tme_results),
      NA
    )
  )
)

write.csv(
  dimension_summary,
  file = file.path(
    final_root,
    "QC",
    "Step14_object_dimensions.csv"
  ),
  row.names = FALSE
)


############################################################
# 14. SAVE LOCKED KEY NUMERICAL RESULTS
############################################################

locked_results <- data.frame(
  
  Item = c(
    "Tumor samples",
    "WGCNA modules",
    "TME traits",
    "Total module-TME tests",
    "FDR < 0.05 associations",
    "ESTIMATEScore FDR < 0.05",
    "ImmuneScore FDR < 0.05",
    "StromalScore FDR < 0.05",
    "Brown internal module",
    "Brown vs ME3 Pearson r",
    "Brown-StromalScore rho",
    "Brown-StromalScore FDR",
    "Brown-ImmuneScore rho",
    "Brown-ImmuneScore FDR",
    "Brown-ESTIMATEScore rho",
    "Brown-ESTIMATEScore FDR"
  ),
  
  Value = c(
    513,
    50,
    3,
    150,
    sum(module_tme_results$FDR < 0.05),
    sum(
      module_tme_results$Trait == "ESTIMATEScore" &
        module_tme_results$FDR < 0.05
    ),
    sum(
      module_tme_results$Trait == "ImmuneScore" &
        module_tme_results$FDR < 0.05
    ),
    sum(
      module_tme_results$Trait == "StromalScore" &
        module_tme_results$FDR < 0.05
    ),
    "ME3",
    1,
    -0.1924747,
    2.616581e-05,
    -0.1621890,
    4.020444e-04,
    -0.1863114,
    4.783196e-05
  )
)

write.csv(
  locked_results,
  file = file.path(
    final_root,
    "QC",
    "Step14A-14E_LOCKED_KEY_RESULTS.csv"
  ),
  row.names = FALSE
)


############################################################
# 15. SAVE SESSION INFORMATION
############################################################

session_file <- file.path(
  final_root,
  "SessionInfo",
  paste0(
    "sessionInfo_",
    save_date,
    ".txt"
  )
)

sink(session_file)

cat("============================================\n")
cat("WORKFLOW-03 Step 14A-14E FINAL LOCKED\n")
cat("============================================\n\n")

cat("Save time:\n")
print(Sys.time())

cat("\n\nR version and package information:\n")
print(sessionInfo())

sink()


############################################################
# 16. COPY THIS SCRIPT INTO THE FINAL DIRECTORY
############################################################

# IMPORTANT:
# Before running this section, save this script itself as:
#
# WORKFLOW-03_Step14_SAVE_FINAL_LOCKED.R
#
# Then modify the source path below if necessary.

# Example:
#
# file.copy(
#   from = "WORKFLOW-03_Step14_SAVE_FINAL_LOCKED.R",
#   to = file.path(
#     final_root,
#     "Scripts",
#     "WORKFLOW-03_Step14_SAVE_FINAL_LOCKED.R"
#   ),
#   overwrite = TRUE
# )


############################################################
# 17. FINAL CHECK
############################################################

cat("\n")
cat("============================================\n")
cat("STEP 14A-14E FINAL LOCKED SAVE COMPLETED\n")
cat("============================================\n\n")

cat("Final directory:\n")
cat(normalizePath(final_root), "\n\n")

cat("Saved objects:\n")
print(key_objects)

cat("\n")
cat("Module-TME tests:", nrow(module_tme_results), "\n")
cat(
  "FDR < 0.05:",
  sum(module_tme_results$FDR < 0.05),
  "\n"
)

cat("\n")
cat("Brown module = ME3\n")
cat("Brown vs ME3 Pearson r = 1\n")

cat("\n")
cat("Brown-TME locked results:\n")
print(
  subset(
    module_tme_results,
    Module == "ME3"
  )
)

cat("\n")
cat("============================================\n")
cat("DO NOT RECOMPUTE STEP 14A-14E UNLESS NEEDED\n")
cat("FOR AN EXPLICIT METHODOLOGICAL CORRECTION.\n")
cat("============================================\n")













