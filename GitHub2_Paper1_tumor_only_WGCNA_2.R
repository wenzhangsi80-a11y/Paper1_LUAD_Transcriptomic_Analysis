
##Step 1 setup.R

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

work_dir <- "input_dir/GitHub2_Paper1_tumor-only WGCNA"

setwd(work_dir)

getwd()

##Step 2：Confirm that the two checkpoint files exist.
candidate_checkpoint <- file.path(
  input_dir,
  "Workflow-04-2D_final_WGCNA_candidate_checkpoint.RData"
)

done_file <- file.path(
  input_dir,
  "Workflow-04-2D_WGCNA_input_checkpoint_DONE.txt"
)

file.exists(candidate_checkpoint)
file.exists(done_file)

##Step 3：Read the DONE file
cat(
  readLines(done_file),
  sep = "\n"
)

##Step 4：Loading .RData
load(candidate_checkpoint)

##Step 5：Perform a quick recovery verification
cat(
  "\nDimensions of datExpr_MAD:",
  nrow(datExpr_MAD),
  "samples ×",
  ncol(datExpr_MAD),
  "genes\n"
)

cat(
  "Dimensions of clinical_aligned:",
  nrow(clinical_aligned),
  "rows ×",
  ncol(clinical_aligned),
  "columns\n"
)

cat(
  "\nNA:",
  sum(is.na(datExpr_MAD)),
  "\n"
)

cat(
  "NaN:",
  sum(is.nan(datExpr_MAD)),
  "\n"
)

cat(
  "Inf:",
  sum(is.infinite(datExpr_MAD)),
  "\n"
)

cat(
  "Duplicated samples:",
  sum(duplicated(rownames(datExpr_MAD))),
  "\n"
)

cat(
  "Duplicated genes:",
  sum(duplicated(colnames(datExpr_MAD))),
  "\n"
)

cat(
  "Sample matching:",
  identical(
    rownames(datExpr_MAD),
    clinical_aligned$Sample
  ),
  "\n"
)

cat(
  "MAD vector length:",
  length(gene_MAD),
  "\n"
)

cat(
  "MAD-positive genes:",
  sum(gene_MAD > 0),
  "\n"
)

############################################################
# Step 6 — WGCNA soft-threshold selection
# Primary analysis: tumor-only TCGA-LUAD
############################################################

library(WGCNA)

options(stringsAsFactors = FALSE)

# Allow WGCNA to use multiple cores if available
allowWGCNAThreads(nThreads = 2)
enableWGCNAThreads()
##In a Windows environment, multi-process or multi-threading scenarios may encounter issues during computation：serialize(data, node$con)
disableWGCNAThreads()
WGCNAnThreads()

#-----------------------------------------------------------
# 6.1 Confirm input
#-----------------------------------------------------------

cat("\n========================================\n")
cat("WGCNA Soft-threshold Selection\n")
cat("========================================\n")

cat(
  "Expression matrix:",
  nrow(datExpr_MAD),
  "samples ×",
  ncol(datExpr_MAD),
  "genes\n"
)

cat(
  "Network type: unsigned\n"
)

cat(
  "Correlation: Pearson\n"
)

cat(
  "Missing values:",
  sum(is.na(datExpr_MAD)),
  "\n"
)

#-----------------------------------------------------------
# 6.2 Candidate soft-thresholding powers
#-----------------------------------------------------------

powers <- c(
  1:10,
  seq(12, 30, by = 2)
)

cat(
  "\nCandidate powers:",
  paste(powers, collapse = ", "),
  "\n"
)

#-----------------------------------------------------------
# 6.3 pickSoftThreshold
#-----------------------------------------------------------

sft <- pickSoftThreshold(
  datExpr_MAD,
  dataIsExpr = TRUE,
  corFnc = "cor",
  corOptions = list(use = "p"),
  networkType = "unsigned",
  powerVector = powers,
  RsquaredCut = 0.85,
  moreNetworkConcepts = TRUE,
  verbose = 5
)

#-----------------------------------------------------------
# 6.4 Save result
#-----------------------------------------------------------

save(
  sft,
  powers,
  file = "result_dir/Workflow-04-3A_WGCNA_soft_threshold.RData"
)

cat(
  "\nSoft-threshold analysis completed.\n"
)

cat(
  "Result saved as:\n",
  "Workflow-04-3A_WGCNA_soft_threshold.RData\n"
)


#===========================================================
# Step 6. Soft-threshold power selection
#===========================================================

selected_power <- 7

cat("\n========================================\n")
cat("Selected soft-threshold power\n")
cat("========================================\n")

cat("Selected power (β):", selected_power, "\n")
cat(
  "Scale-free topology fit (R²):",
  sft$fitIndices[
    sft$fitIndices$Power == selected_power,
    "SFT.R.sq"
  ],
  "\n"
)

cat(
  "Mean connectivity:",
  sft$fitIndices[
    sft$fitIndices$Power == selected_power,
    "mean.k."
  ],
  "\n"
)

cat(
  "Median connectivity:",
  sft$fitIndices[
    sft$fitIndices$Power == selected_power,
    "median.k."
  ],
  "\n"
)

cat(
  "Maximum connectivity:",
  sft$fitIndices[
    sft$fitIndices$Power == selected_power,
    "max.k."
  ],
  "\n"
)
##6. Save Step 6 checkpoint
saveRDS(
  sft,
  file = "result_dir/Workflow-04-2E_WGCNA_soft_threshold.rds"
)

write.table(
  sft$fitIndices,
  file = "result_dir/Workflow-04-2E_WGCNA_soft_threshold_fitIndices.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
##Perform a simple integrity check

file.exists("Workflow-04-2E_WGCNA_soft_threshold.rds")

sft_check <- readRDS(
  "Workflow-04-2E_WGCNA_soft_threshold.rds"
)

dim(sft_check$fitIndices)

sft_check$fitIndices
sft_check$powerEstimate

##Let’s refine the Step 6 checkpoint a bit.

WGCNA_Step6 <- list(
  sft = sft,
  candidate_powers = powers,
  selected_power = 7,
  networkType = "unsigned",
  correlation = "Pearson",
  RsquaredCut = 0.85,
  n_samples = nrow(datExpr_MAD),
  n_genes = ncol(datExpr_MAD)
)

saveRDS(
  WGCNA_Step6,
  file = "result_dir/Workflow-04-2E_WGCNA_soft_threshold_checkpoint.rds"
)

###7. It is recommended to also generate a formal soft-thresholding plot.

fit <- sft$fitIndices

plot(
  fit$Power,
  -sign(fit$slope) * fit$SFT.R.sq,
  type = "n",
  xlab = "Soft-threshold power (β)",
  ylab = "Scale-free topology fit (signed R²)",
  main = "Soft-threshold selection"
)

text(
  fit$Power,
  -sign(fit$slope) * fit$SFT.R.sq,
  labels = fit$Power,
  cex = 0.8
)

abline(h = 0.85, lty = 2)

abline(
  v = selected_power,
  lty = 2
)

##然后 connectivity：
plot(
  fit$Power,
  fit$mean.k.,
  type = "n",
  xlab = "Soft-threshold power (β)",
  ylab = "Mean connectivity",
  main = "Mean connectivity across powers"
)

text(
  fit$Power,
  fit$mean.k.,
  labels = fit$Power,
  cex = 0.8
)

abline(
  v = selected_power,
  lty = 2
)

#===========================================================
# Step 7. Tumor-only WGCNA module construction
#===========================================================

library(WGCNA)

options(stringsAsFactors = FALSE)

cat("\n========================================\n")
cat("Step 7: Tumor-only WGCNA\n")
cat("Module construction using blockwiseModules()\n")
cat("========================================\n\n")

#-----------------------------------------------------------
# 7.1 Confirm input
#-----------------------------------------------------------
#7.1 Initialization and Data Confirmation

cat(
  "Expression matrix:",
  nrow(datExpr_MAD),
  "samples ×",
  ncol(datExpr_MAD),
  "genes\n"
)

cat(
  "Missing values:",
  sum(is.na(datExpr_MAD)),
  "\n"
)

cat(
  "Duplicated samples:",
  sum(duplicated(rownames(datExpr_MAD))),
  "\n"
)

cat(
  "Duplicated genes:",
  sum(duplicated(colnames(datExpr_MAD))),
  "\n"
)

#-----------------------------------------------------------
# 7.2 WGCNA parameters
#-----------------------------------------------------------

softPower <- 7

networkType <- "unsigned"
TOMType <- "unsigned"

minModuleSize <- 30
deepSplit <- 3
mergeCutHeight <- 0.25

cat("\nWGCNA parameters:\n")
cat("Soft-threshold power (β):", softPower, "\n")
cat("Network type:", networkType, "\n")
cat("TOM type:", TOMType, "\n")
cat("Minimum module size:", minModuleSize, "\n")
cat("Deep split:", deepSplit, "\n")
cat("Module merge cut height:", mergeCutHeight, "\n")

#-----------------------------------------------------------
# 7.3 Keep WGCNA single-threaded
#-----------------------------------------------------------

disableWGCNAThreads()

cat(
  "\nWGCNA threads:",
  WGCNAnThreads(),
  "\n"
)

#-----------------------------------------------------------
# 7.4 Confirm sample alignment
#-----------------------------------------------------------

stopifnot(
  identical(
    rownames(datExpr_MAD),
    clinical_aligned$Sample
  )
)

cat(
  "Sample order check: PASSED\n"
)


#-----------------------------------------------------------
# 7.5 Construct co-expression modules
#-----------------------------------------------------------

cat("\n========================================\n")
cat("Running blockwiseModules()\n")
cat("This may take some time...\n")
cat("========================================\n\n")

net <- blockwiseModules(
  datExpr_MAD,
  
  # Network construction
  power = softPower,
  networkType = networkType,
  TOMType = TOMType,
  
  # Module detection
  minModuleSize = minModuleSize,
  deepSplit = deepSplit,
  
  # Module merging
  mergeCutHeight = mergeCutHeight,
  
  # Module eigengenes
  reassignThreshold = 0,
  pamRespectsDendro = FALSE,
  
  # TOM calculation
  saveTOMs = FALSE,
  loadTOMs = FALSE,
  
  # Computational settings
  numericLabels = TRUE,
  verbose = 3
)

#-----------------------------------------------------------
# 7.6 Extract module colors
#-----------------------------------------------------------

moduleColors <- labels2colors(net$colors)

cat("\n========================================\n")
cat("Module summary\n")
cat("========================================\n\n")

cat(
  "Number of detected modules:",
  length(unique(moduleColors)),
  "\n"
)

cat(
  "Number of genes assigned to modules:",
  sum(moduleColors != "grey"),
  "\n"
)

cat(
  "Number of grey/unassigned genes:",
  sum(moduleColors == "grey"),
  "\n"
)

#-----------------------------------------------------------
# 7.7 Module size distribution
#-----------------------------------------------------------

moduleSizes <- sort(
  table(moduleColors),
  decreasing = TRUE
)

moduleSizes

#-----------------------------------------------------------
# 7.8 Module size plot
#-----------------------------------------------------------

moduleSize_df <- data.frame(
  Module = names(moduleSizes),
  GeneCount = as.numeric(moduleSizes)
)

moduleSize_df$Module <- factor(
  moduleSize_df$Module,
  levels = names(moduleSizes)
)

library(ggplot2)

p_module_size <- ggplot(
  moduleSize_df,
  aes(
    x = Module,
    y = GeneCount
  )
) +
  geom_col() +
  theme_classic() +
  labs(
    x = "WGCNA module",
    y = "Number of genes"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_module_size


#-----------------------------------------------------------
# 7.9 Extract module eigengenes
#-----------------------------------------------------------

MEs0 <- net$MEs

MEs <- orderMEs(MEs0)

cat(
  "\nModule eigengene matrix:",
  nrow(MEs),
  "samples ×",
  ncol(MEs),
  "modules\n"
)

#-----------------------------------------------------------
# 7.10 Check module eigengenes
#-----------------------------------------------------------

cat(
  "Missing values in module eigengenes:",
  sum(is.na(MEs)),
  "\n"
)

#-----------------------------------------------------------
# 7.11 Save Step 7 checkpoint
#-----------------------------------------------------------

WGCNA_Step7 <- list(
  
  net = net,
  
  moduleColors = moduleColors,
  
  moduleSizes = moduleSizes,
  
  MEs = MEs,
  
  softPower = softPower,
  
  networkType = networkType,
  
  TOMType = TOMType,
  
  minModuleSize = minModuleSize,
  
  deepSplit = deepSplit,
  
  mergeCutHeight = mergeCutHeight,
  
  nSamples = nrow(datExpr_MAD),
  
  nGenes = ncol(datExpr_MAD)
)

saveRDS(
  WGCNA_Step7,
  file = "Workflow-04-2F_WGCNA_tumor_only_modules.rds"
)

cat(
  "\nStep 7 checkpoint saved successfully.\n"
)

#-----------------------------------------------------------
# 7.12 Save gene-module assignments
#-----------------------------------------------------------

gene_module_assignment <- data.frame(
  Gene = colnames(datExpr_MAD),
  Module = moduleColors,
  stringsAsFactors = FALSE
)

write.table(
  gene_module_assignment,
  file = "result_dir/Workflow-04-2F_gene_module_assignment.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


subset(
  gene_module_assignment,
  Gene %in% c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )
)
# 1. Check whether the three genes are present in the current WGCNA input matrix.
c(
  AOC3 = "AOC3" %in% colnames(datExpr_MAD),
  MS4A7 = "MS4A7" %in% colnames(datExpr_MAD),
  IL7R = "IL7R" %in% colnames(datExpr_MAD)
)

# 2. Search for column names that contain these three gene names.
grep(
  "AOC3|MS4A7|IL7R",
  colnames(datExpr_MAD),
  value = TRUE,
  ignore.case = TRUE
)
# 3. View the assignment data frame structure
head(gene_module_assignment)

str(gene_module_assignment)

nrow(gene_module_assignment)

# Check whether the Ensembl IDs in the current expression matrix include version numbers.
head(colnames(datExpr_MAD), 20)

# Check for the presence of an Ensembl ID + version number format.
grep(
  "\\.",
  colnames(datExpr_MAD)[1:100],
  value = TRUE
)

# View the column names of the original/candidate expression matrix object
head(colnames(datExpr_MAD))

# View objects in the current workspace
ls()

# Step 7A-1：Check the structure of gene_MAD
str(gene_MAD)

# View the first 20 gene_MADs
head(gene_MAD, 20)

# Check whether the target gene is contained within gene_MAD.
gene_MAD[
  grepl(
    "AOC3|MS4A7|IL7R",
    gene_MAD,
    ignore.case = TRUE
  )
]
#Step 7A-2: Check the structure of `candidate_checkpoint`
str(candidate_checkpoint)

# View its name
names(candidate_checkpoint)

# Step 7A-3：Check for the target gene within the Ensembl IDs in gene_MAD
target_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

grep(
  paste(target_genes, collapse = "|"),
  names(gene_MAD),
  value = TRUE,
  ignore.case = TRUE
)


# Step 7A-4：Check which objects are saved in the original checkpoint
tmp_env <- new.env()

load(
  candidate_checkpoint,
  envir = tmp_env
)

ls(tmp_env)

# Examine the basic structure of the objects within a checkpoint
lapply(
  ls(tmp_env),
  function(x) {
    obj <- get(x, envir = tmp_env)
    c(
      object = x,
      class = class(obj)[1],
      length = length(obj)
    )
  }
)

# Step 7A-5：Locate the original expression file
list.files(
  recursive = TRUE,
  pattern = "input_dir/Merge_RNA_seq_FPKM_UQ",
  full.names = TRUE
)

# Step 7A-6：Check only the header and the first few rows of the raw expression file.
raw_file <- "input_dir/Merge_RNA_seq_FPKM_UQ.txt"

raw_header <- readLines(
  raw_file,
  n = 3
)

cat(
  paste(raw_header, collapse = "\n")
)

# Step 7A-7：Check whether the human gene annotation package is available
requireNamespace(
  "org.Hs.eg.db",
  quietly = TRUE
)

# Step 7A-8：Load human gene annotations
library(
  AnnotationDbi
)

library(
  org.Hs.eg.db
)
# Step 7A-9：Check the Ensembl IDs for the three target genes.
AnnotationDbi::select(
  org.Hs.eg.db,
  keys = c(
    "ENSG00000131471",
    "ENSG00000166927",
    "ENSG00000168685"
  ),
  keytype = "ENSEMBL",
  columns = c(
    "ENSEMBL",
    "SYMBOL"
  )
)

# Step 7A-10：WGCNA module to which AOC3 / MS4A7 / IL7R are assigned

target_ensembl <- c(
  "ENSG00000131471",  # AOC3
  "ENSG00000166927",  # MS4A7
  "ENSG00000168685"   # IL7R
)

target_symbols <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

target_annotation <- data.frame(
  Gene = target_symbols,
  Ensembl = target_ensembl,
  stringsAsFactors = FALSE
)

target_annotation$Module <- moduleColors[
  match(
    target_annotation$Ensembl,
    colnames(datExpr_MAD)
  )
]

target_annotation

# Step 7A-11：Confirm whether the three genes are included in the current WGCNA input matrix.

data.frame(
  Gene = target_symbols,
  Ensembl = target_ensembl,
  In_WGCNA = target_ensembl %in% colnames(datExpr_MAD),
  stringsAsFactors = FALSE
)

# Step 7A-12：Create a comprehensive Ensembl ID–Gene Symbol–Module mapping table

all_ensembl <- colnames(datExpr_MAD)

gene_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = all_ensembl,
  keytype = "ENSEMBL",
  columns = c(
    "ENSEMBL",
    "SYMBOL"
  )
)

head(gene_annotation)

# Check the number of mappings
nrow(gene_annotation)

# Check for duplicate Ensembl IDs
sum(duplicated(gene_annotation$ENSEMBL))

# Check the number of entries lacking a Gene Symbol
sum(
  is.na(gene_annotation$SYMBOL) |
    gene_annotation$SYMBOL == ""
)

# Step 7A-13：Calculate the Module Membership for AOC3 / MS4A7 / IL7R

# Ensure that the column names of the module eigengenes match `moduleColors`.
MEs_ordered <- orderMEs(MEs)

#The module containing the three core genes
target_modules <- target_annotation$Module

# Calculate the correlation between all genes and each module eigengene.
MM_all <- cor(
  datExpr_MAD,
  MEs_ordered,
  use = "p",
  method = "pearson"
)

# Extract three core genes
target_MM <- data.frame(
  Gene = target_symbols,
  Ensembl = target_ensembl,
  Module = target_modules,
  stringsAsFactors = FALSE
)

target_MM$MM <- mapply(
  function(gene, module) {
    
    module_ME <- paste0("ME", module)
    
    MM_all[
      gene,
      module_ME
    ]
    
  },
  target_MM$Ensembl,
  target_MM$Module
)

target_MM


# Step 7A-13a：Check the column names of the Module Eigengene.

colnames(MEs_ordered)

# Examine the modules corresponding to the three core genes.
target_annotation[, c("Gene", "Ensembl", "Module")]

# Check the dimensions and column names of MM_all
dim(MM_all)
colnames(MM_all)

# Step 7A-13b：Establish a mapping between module colors and numeric labels.

module_label_color <- data.frame(
  Label = sort(unique(net$colors)),
  Color = labels2colors(sort(unique(net$colors))),
  stringsAsFactors = FALSE
)

module_label_color

# Step 7A-13c：Convert the module color into its corresponding numeric label.

target_MM <- data.frame(
  Gene = target_symbols,
  Ensembl = target_ensembl,
  Module = target_modules,
  stringsAsFactors = FALSE
)

target_MM$Label <- module_label_color$Label[
  match(target_MM$Module, module_label_color$Color)
]

target_MM$ME <- paste0("ME", target_MM$Label)

target_MM

# Step 7A-13d：Extract the Module Membership for AOC3 / MS4A7 / IL7R

target_MM$MM <- mapply(
  function(gene, ME) {
    MM_all[gene, ME]
  },
  target_MM$Ensembl,
  target_MM$ME
)

target_MM


# Step 8-1：Check the structure of clinical_aligned

dim(clinical_aligned)

colnames(clinical_aligned)

str(clinical_aligned)

summary(clinical_aligned)

# Step 8-2：Check the unique values for each clinical trait.

lapply(
  clinical_aligned,
  function(x) {
    if (is.factor(x)) {
      levels(x)
    } else {
      unique(x)
    }
  }
)

# Step 8-pre：Save current WGCNA + MM checkpoint

WGCNA_Step8_pre <- list(
  net = net,
  moduleColors = moduleColors,
  module_label_color = module_label_color,
  MEs = MEs_ordered,
  MM_all = MM_all,
  target_MM = target_MM,
  clinical_aligned = clinical_aligned,
  datExpr_MAD = datExpr_MAD,
  softPower = softPower,
  networkType = networkType,
  TOMType = TOMType,
  minModuleSize = minModuleSize,
  deepSplit = deepSplit,
  mergeCutHeight = mergeCutHeight
)

saveRDS(
  WGCNA_Step8_pre,
  file = "result_dir/Workflow-04-2G_WGCNA_tumor_only_MM_checkpoint.rds"
)


colnames(clinical_aligned)

table(clinical_aligned$Gender, useNA = "ifany")

table(clinical_aligned$pathologic_stage_2, useNA = "ifany")

table(clinical_aligned$pathologic_t, useNA = "ifany")

table(clinical_aligned$pathologic_n, useNA = "ifany")

table(clinical_aligned$pathologic_m, useNA = "ifany")

table(clinical_aligned$OS.status, useNA = "ifany")

summary(clinical_aligned$Age)

summary(clinical_aligned$OS.time)


# ============================================================
# Step 8A. Construct tumor-only WGCNA clinical trait matrix
# ============================================================

clinical_WGCNA_traits <- clinical_aligned

# ------------------------------------------------------------
# 1. Age
# ------------------------------------------------------------

clinical_WGCNA_traits$Age_trait <-
  clinical_WGCNA_traits$Age


# ------------------------------------------------------------
# 2. Gender
# female = 0
# male   = 1
# ------------------------------------------------------------

clinical_WGCNA_traits$Gender_trait <-
  ifelse(
    clinical_WGCNA_traits$Gender == "male",
    1,
    ifelse(
      clinical_WGCNA_traits$Gender == "female",
      0,
      NA
    )
  )


# ------------------------------------------------------------
# 3. Pathological stage
# Combine Stage I/IA/IB, Stage II/IIA/IIB,
# Stage IIIA/IIIB, and Stage IV
# ------------------------------------------------------------

clinical_WGCNA_traits$Stage_trait <- NA_character_

clinical_WGCNA_traits$Stage_trait[
  clinical_WGCNA_traits$pathologic_stage_2 %in%
    c("Stage I", "Stage IA", "Stage IB")
] <- "Stage I"

clinical_WGCNA_traits$Stage_trait[
  clinical_WGCNA_traits$pathologic_stage_2 %in%
    c("Stage II", "Stage IIA", "Stage IIB")
] <- "Stage II"

clinical_WGCNA_traits$Stage_trait[
  clinical_WGCNA_traits$pathologic_stage_2 %in%
    c("Stage IIIA", "Stage IIIB")
] <- "Stage III"

clinical_WGCNA_traits$Stage_trait[
  clinical_WGCNA_traits$pathologic_stage_2 == "Stage IV"
] <- "Stage IV"

clinical_WGCNA_traits$Stage_trait <-
  factor(
    clinical_WGCNA_traits$Stage_trait,
    levels = c(
      "Stage I",
      "Stage II",
      "Stage III",
      "Stage IV"
    )
  )


# ------------------------------------------------------------
# 4. T stage
# Combine T1/T1a/T1b and T2/T2a/T2b
# Exclude TX
# ------------------------------------------------------------

clinical_WGCNA_traits$T_trait <- NA_character_

clinical_WGCNA_traits$T_trait[
  clinical_WGCNA_traits$pathologic_t %in%
    c("T1", "T1a", "T1b")
] <- "T1"

clinical_WGCNA_traits$T_trait[
  clinical_WGCNA_traits$pathologic_t %in%
    c("T2", "T2a", "T2b")
] <- "T2"

clinical_WGCNA_traits$T_trait[
  clinical_WGCNA_traits$pathologic_t == "T3"
] <- "T3"

clinical_WGCNA_traits$T_trait[
  clinical_WGCNA_traits$pathologic_t == "T4"
] <- "T4"

clinical_WGCNA_traits$T_trait <-
  factor(
    clinical_WGCNA_traits$T_trait,
    levels = c("T1", "T2", "T3", "T4")
  )


# ------------------------------------------------------------
# 5. N stage
# Keep N0/N1/N2/N3
# Exclude NX
# ------------------------------------------------------------

clinical_WGCNA_traits$N_trait <- NA_character_

clinical_WGCNA_traits$N_trait[
  clinical_WGCNA_traits$pathologic_n %in%
    c("N0", "N1", "N2", "N3")
] <- clinical_WGCNA_traits$pathologic_n[
  clinical_WGCNA_traits$pathologic_n %in%
    c("N0", "N1", "N2", "N3")
]

clinical_WGCNA_traits$N_trait <-
  factor(
    clinical_WGCNA_traits$N_trait,
    levels = c("N0", "N1", "N2", "N3")
  )


# ------------------------------------------------------------
# 6. M stage
# Combine M1/M1a/M1b into M1
# Exclude MX
# ------------------------------------------------------------

clinical_WGCNA_traits$M_trait <- NA_character_

clinical_WGCNA_traits$M_trait[
  clinical_WGCNA_traits$pathologic_m == "M0"
] <- "M0"

clinical_WGCNA_traits$M_trait[
  clinical_WGCNA_traits$pathologic_m %in%
    c("M1", "M1a", "M1b")
] <- "M1"

clinical_WGCNA_traits$M_trait <-
  factor(
    clinical_WGCNA_traits$M_trait,
    levels = c("M0", "M1")
  )


# ------------------------------------------------------------
# 7. Survival variables
# ------------------------------------------------------------

clinical_WGCNA_traits$OS_time <-
  clinical_WGCNA_traits$OS.time

clinical_WGCNA_traits$OS_event <-
  ifelse(
    clinical_WGCNA_traits$OS.status == "Dead",
    1,
    ifelse(
      clinical_WGCNA_traits$OS.status == "Alive",
      0,
      NA
    )
  )

##After running it, check only these few results first.
table(
  clinical_WGCNA_traits$Gender_trait,
  useNA = "ifany"
)

table(
  clinical_WGCNA_traits$Stage_trait,
  useNA = "ifany"
)

table(
  clinical_WGCNA_traits$T_trait,
  useNA = "ifany"
)

table(
  clinical_WGCNA_traits$N_trait,
  useNA = "ifany"
)

table(
  clinical_WGCNA_traits$M_trait,
  useNA = "ifany"
)

table(
  clinical_WGCNA_traits$OS_event,
  useNA = "ifany"
)

summary(clinical_WGCNA_traits$Age_trait)

summary(clinical_WGCNA_traits$OS_time)

# ============================================================
# Step 8B-1. Prepare biological module eigengenes
# ============================================================

# Recalculate ordered module eigengenes
MEs_ordered <- orderMEs(MEs)

# Build numeric-label to module-color mapping
module_label_color <- data.frame(
  Label = sort(unique(net$colors)),
  Color = labels2colors(sort(unique(net$colors))),
  stringsAsFactors = FALSE
)

# Match ME columns to module colors
ME_labels <- as.integer(
  sub("^ME", "", colnames(MEs_ordered))
)

ME_colors <- module_label_color$Color[
  match(ME_labels, module_label_color$Label)
]

# Remove grey module
keep_ME <- ME_colors != "grey"

MEs_biological <- MEs_ordered[, keep_ME, drop = FALSE]

# Rename columns using module colors
colnames(MEs_biological) <- paste0(
  "ME_",
  ME_colors[keep_ME]
)

# Check dimensions
dim(MEs_biological)

# Check module names
colnames(MEs_biological)

##3. Perform another very important check on the sample sequence.
identical(
  rownames(MEs_biological),
  clinical_WGCNA_traits$Sample
)

dim(MEs_biological)

colnames(MEs_biological)

identical(
  rownames(MEs_biological),
  clinical_WGCNA_traits$Sample
)

# ============================================================
# Step 8B-2A. Module eigengene association with Age
# ============================================================

age_results <- data.frame(
  Module = colnames(MEs_biological),
  Correlation = NA_real_,
  P_value = NA_real_,
  N = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_len(ncol(MEs_biological))) {
  
  module_values <- MEs_biological[, i]
  age_values <- clinical_WGCNA_traits$Age_trait
  
  complete_idx <- complete.cases(
    module_values,
    age_values
  )
  
  age_results$N[i] <- sum(complete_idx)
  
  if (sum(complete_idx) >= 3) {
    
    test <- cor.test(
      module_values[complete_idx],
      age_values[complete_idx],
      method = "pearson"
    )
    
    age_results$Correlation[i] <- unname(
      test$estimate
    )
    
    age_results$P_value[i] <- test$p.value
  }
}

# Multiple-testing correction
age_results$FDR <- p.adjust(
  age_results$P_value,
  method = "BH"
)

# Sort by P value
age_results <- age_results[
  order(age_results$P_value),
]

rownames(age_results) <- NULL

age_results

# ============================================================
# Step 8B-2B. Module eigengene association with Gender
# ============================================================

gender_results <- data.frame(
  Module = colnames(MEs_biological),
  P_value = NA_real_,
  Female_median = NA_real_,
  Male_median = NA_real_,
  N_female = NA_integer_,
  N_male = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_len(ncol(MEs_biological))) {
  
  module_values <- MEs_biological[, i]
  gender_values <- clinical_WGCNA_traits$Gender_trait
  
  complete_idx <- complete.cases(
    module_values,
    gender_values
  )
  
  female_values <- module_values[
    complete_idx &
      gender_values == 0
  ]
  
  male_values <- module_values[
    complete_idx &
      gender_values == 1
  ]
  
  gender_results$N_female[i] <-
    length(female_values)
  
  gender_results$N_male[i] <-
    length(male_values)
  
  if (
    length(female_values) >= 2 &&
    length(male_values) >= 2
  ) {
    
    test <- wilcox.test(
      female_values,
      male_values,
      alternative = "two.sided"
    )
    
    gender_results$P_value[i] <-
      test$p.value
    
    gender_results$Female_median[i] <-
      median(female_values)
    
    gender_results$Male_median[i] <-
      median(male_values)
  }
}

# Multiple-testing correction
gender_results$FDR <- p.adjust(
  gender_results$P_value,
  method = "BH"
)

# Sort by P value
gender_results <- gender_results[
  order(gender_results$P_value),
]

rownames(gender_results) <- NULL

gender_results

# ============================================================
# Step 8B-3A. Module eigengene association with pathological stage
# Tumor-only WGCNA
# ============================================================

# ------------------------------------------------------------
# 1. Initialize result table
# ------------------------------------------------------------

stage_results <- data.frame(
  Module = colnames(MEs_biological),
  P_value = NA_real_,
  N = NA_integer_,
  N_Stage_I = NA_integer_,
  N_Stage_II = NA_integer_,
  N_Stage_III = NA_integer_,
  N_Stage_IV = NA_integer_,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 2. Kruskal–Wallis test for each module
# ------------------------------------------------------------

for (i in seq_len(ncol(MEs_biological))) {
  
  module_values <- MEs_biological[, i]
  stage_values <- clinical_WGCNA_traits$Stage_trait
  
  complete_idx <- complete.cases(
    module_values,
    stage_values
  )
  
  module_complete <- module_values[complete_idx]
  stage_complete <- droplevels(
    stage_values[complete_idx]
  )
  
  # Total number of complete samples
  stage_results$N[i] <- length(module_complete)
  
  # Number of samples in each stage
  stage_results$N_Stage_I[i] <-
    sum(stage_complete == "Stage I")
  
  stage_results$N_Stage_II[i] <-
    sum(stage_complete == "Stage II")
  
  stage_results$N_Stage_III[i] <-
    sum(stage_complete == "Stage III")
  
  stage_results$N_Stage_IV[i] <-
    sum(stage_complete == "Stage IV")
  
  # Kruskal–Wallis test
  if (
    length(unique(stage_complete)) >= 2
  ) {
    
    test <- kruskal.test(
      module_complete ~ stage_complete
    )
    
    stage_results$P_value[i] <-
      test$p.value
  }
}


# ------------------------------------------------------------
# 3. Multiple-testing correction
# ------------------------------------------------------------

stage_results$FDR <- p.adjust(
  stage_results$P_value,
  method = "BH"
)


# ------------------------------------------------------------
# 4. Sort by nominal P value
# ------------------------------------------------------------

stage_results <- stage_results[
  order(stage_results$P_value),
]

rownames(stage_results) <- NULL


# ------------------------------------------------------------
# 5. Display complete results
# ------------------------------------------------------------

stage_results

##6. It is recommended to formally save the current stage results.

write.table(
  stage_results,
  file = "result_dir/Workflow-04-2G_WGCNA_stage_association.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  stage_results,
  file = "result_dir/Workflow-04-2G_WGCNA_stage_association.rds"
)

#It is also recommended to maintain a dedicated table of significant modules：
stage_sig <- stage_results[
  stage_results$FDR < 0.05,
  ,
  drop = FALSE
]

write.table(
  stage_sig,
  file = "result_dir/Workflow-04-2G_WGCNA_stage_significant_modules.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

stage_sig


# ============================================================
# Step 8B-3B. Module eigengene association with pathological T stage
# Tumor-only WGCNA
#
# T stage:
#   T1
#   T2
#   T3
#   T4
#
# TX and NA are excluded
#
# Statistical test:
#   Kruskal-Wallis test
#
# Multiple-testing correction:
#   Benjamini-Hochberg (BH)
# ============================================================


# ------------------------------------------------------------
# 1. Check T-stage distribution
# ------------------------------------------------------------

table(
  clinical_WGCNA_traits$T_trait,
  useNA = "ifany"
)


# ------------------------------------------------------------
# 2. Confirm that T_trait contains the expected categories
# ------------------------------------------------------------

levels(
  clinical_WGCNA_traits$T_trait
)

unique(
  clinical_WGCNA_traits$T_trait
)


# ------------------------------------------------------------
# 3. Initialize results table
# ------------------------------------------------------------

T_results <- data.frame(
  Module = colnames(MEs_biological),
  P_value = NA_real_,
  N = NA_integer_,
  N_T1 = NA_integer_,
  N_T2 = NA_integer_,
  N_T3 = NA_integer_,
  N_T4 = NA_integer_,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 4. Kruskal-Wallis test for each biological module
# ------------------------------------------------------------

for (i in seq_len(ncol(MEs_biological))) {
  
  module_values <- MEs_biological[, i]
  
  T_values <- clinical_WGCNA_traits$T_trait
  
  
  # Complete cases only
  complete_idx <- complete.cases(
    module_values,
    T_values
  )
  
  
  module_complete <- module_values[complete_idx]
  
  T_complete <- droplevels(
    T_values[complete_idx]
  )
  
  
  # Sample size
  T_results$N[i] <- length(module_complete)
  
  
  # Sample size of each T category
  T_results$N_T1[i] <-
    sum(T_complete == "T1")
  
  T_results$N_T2[i] <-
    sum(T_complete == "T2")
  
  T_results$N_T3[i] <-
    sum(T_complete == "T3")
  
  T_results$N_T4[i] <-
    sum(T_complete == "T4")
  
  
  # Kruskal-Wallis test
  if (
    length(unique(T_complete)) >= 2
  ) {
    
    test <- kruskal.test(
      module_complete ~ T_complete
    )
    
    T_results$P_value[i] <-
      test$p.value
  }
}


# ------------------------------------------------------------
# 5. Benjamini-Hochberg correction
# ------------------------------------------------------------

T_results$FDR <- p.adjust(
  T_results$P_value,
  method = "BH"
)


# ------------------------------------------------------------
# 6. Sort by nominal P value
# ------------------------------------------------------------

T_results <- T_results[
  order(T_results$P_value),
]

rownames(T_results) <- NULL


# ------------------------------------------------------------
# 7. Display complete results
# ------------------------------------------------------------

T_results

# ============================================================
# Step 8B-3B QC
# ============================================================

# Number of modules
nrow(T_results)


# Sample size used for each module
unique(T_results$N)


# Distribution of T stages used in the analysis
T_results[
  ,
  c(
    "Module",
    "N",
    "N_T1",
    "N_T2",
    "N_T3",
    "N_T4"
  )
]


# Original T-stage distribution
table(
  clinical_WGCNA_traits$T_trait,
  useNA = "ifany"
)

unique(T_results$N)

# ============================================================
# Check ME_blue
# ============================================================

T_results[
  T_results$Module == "ME_blue",
  ,
  drop = FALSE
]

# ============================================================
# Check ME_magenta
# ============================================================

T_results[
  T_results$Module == "ME_magenta",
  ,
  drop = FALSE
]

# ============================================================
# Significant modules after BH correction
# ============================================================

T_sig <- T_results[
  T_results$FDR < 0.05,
  ,
  drop = FALSE
]

T_sig

nrow(T_sig)

T_sig[
  ,
  c(
    "Module",
    "P_value",
    "FDR",
    "N",
    "N_T1",
    "N_T2",
    "N_T3",
    "N_T4"
  )
]


# ============================================================
# 14.Save T-stage association results
# ============================================================

write.table(
  T_results,
  file = "result_dir/Workflow-04-2H_WGCNA_T_stage_association.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


saveRDS(
  T_results,
  file = "result_dir/Workflow-04-2H_WGCNA_T_stage_association.rds"
)


# Save significant modules separately

write.table(
  T_sig,
  file = "result_dir/Workflow-04-2H_WGCNA_T_stage_significant_modules.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Save Step 8B-3B results
# ============================================================

write.table(
  T_results,
  file = "result_dir/Workflow-04-2H_WGCNA_T_stage_association.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  T_results,
  file = "result_dir/Workflow-04-2H_WGCNA_T_stage_association.rds"
)

write.table(
  T_sig,
  file = "result_dir/Workflow-04-2H_WGCNA_T_stage_significant_modules.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Step 8B-3C. Module eigengene association with pathological N stage
# Tumor-only WGCNA
#
# N stage:
#   N0
#   N1
#   N2
#   N3
#
# NX and NA are excluded
#
# Statistical test:
#   Kruskal-Wallis test
#
# Multiple-testing correction:
#   Benjamini-Hochberg (BH)
#
# Important:
#   N3 contains only 2 samples.
#   No pairwise post-hoc analysis is performed at this stage.
# ============================================================


# ------------------------------------------------------------
# 1. Check N-stage distribution
# ------------------------------------------------------------

table(
  clinical_WGCNA_traits$N_trait,
  useNA = "ifany"
)


# ------------------------------------------------------------
# 2. Confirm N-stage categories
# ------------------------------------------------------------

levels(
  clinical_WGCNA_traits$N_trait
)

unique(
  clinical_WGCNA_traits$N_trait
)


# ------------------------------------------------------------
# 3. Initialize results table
# ------------------------------------------------------------

N_results <- data.frame(
  Module = colnames(MEs_biological),
  P_value = NA_real_,
  N = NA_integer_,
  N_N0 = NA_integer_,
  N_N1 = NA_integer_,
  N_N2 = NA_integer_,
  N_N3 = NA_integer_,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 4. Kruskal-Wallis test for each biological module
# ------------------------------------------------------------

for (i in seq_len(ncol(MEs_biological))) {
  
  module_values <- MEs_biological[, i]
  
  N_values <- clinical_WGCNA_traits$N_trait
  
  
  # Complete cases only
  complete_idx <- complete.cases(
    module_values,
    N_values
  )
  
  
  module_complete <- module_values[complete_idx]
  
  N_complete <- droplevels(
    N_values[complete_idx]
  )
  
  
  # Sample size used for this module
  N_results$N[i] <- length(module_complete)
  
  
  # Sample size of each N category
  N_results$N_N0[i] <-
    sum(N_complete == "N0")
  
  N_results$N_N1[i] <-
    sum(N_complete == "N1")
  
  N_results$N_N2[i] <-
    sum(N_complete == "N2")
  
  N_results$N_N3[i] <-
    sum(N_complete == "N3")
  
  
  # Kruskal-Wallis test
  if (
    length(unique(N_complete)) >= 2
  ) {
    
    test <- kruskal.test(
      module_complete ~ N_complete
    )
    
    N_results$P_value[i] <-
      test$p.value
  }
}


# ------------------------------------------------------------
# 5. Benjamini-Hochberg correction
# ------------------------------------------------------------

N_results$FDR <- p.adjust(
  N_results$P_value,
  method = "BH"
)


# ------------------------------------------------------------
# 6. Sort by nominal P value
# ------------------------------------------------------------

N_results <- N_results[
  order(N_results$P_value),
]

rownames(N_results) <- NULL


# ------------------------------------------------------------
# 7. Display complete results
# ------------------------------------------------------------

N_results

# ============================================================
# Step 8B-3C QC
# ============================================================

# Number of modules
nrow(N_results)


# Sample size used for each module
unique(N_results$N)


# Distribution of N stages used in the analysis
N_results[
  ,
  c(
    "Module",
    "N",
    "N_N0",
    "N_N1",
    "N_N2",
    "N_N3"
  )
]


# Original N-stage distribution
table(
  clinical_WGCNA_traits$N_trait,
  useNA = "ifany"
)

unique(N_results$N)

# ============================================================
# Check ME_blue
# ============================================================

N_results[
  N_results$Module == "ME_blue",
  ,
  drop = FALSE
]

# ============================================================
# Check ME_magenta
# ============================================================

N_results[
  N_results$Module == "ME_magenta",
  ,
  drop = FALSE
]

# ============================================================
# Significant modules after BH correction
# ============================================================

N_sig <- N_results[
  N_results$FDR < 0.05,
  ,
  drop = FALSE
]

N_sig

nrow(N_sig)

N_sig[
  ,
  c(
    "Module",
    "P_value",
    "FDR",
    "N",
    "N_N0",
    "N_N1",
    "N_N2",
    "N_N3"
  )
]

table(
  clinical_WGCNA_traits$N_trait,
  useNA = "ifany"
)

# ============================================================
# Save N-stage association results
# ============================================================

write.table(
  N_results,
  file = "result_dir/Workflow-04-2I_WGCNA_N_stage_association.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


saveRDS(
  N_results,
  file = "result_dir/Workflow-04-2I_WGCNA_N_stage_association.rds"
)


# Save significant modules separately

write.table(
  N_sig,
  file = "result_dir/Workflow-04-2I_WGCNA_N_stage_significant_modules.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Save Step 8B-3C results
# ============================================================

write.table(
  N_results,
  file = "result_dir/Workflow-04-2I_WGCNA_N_stage_association.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  N_results,
  file = "result_dir/Workflow-04-2I_WGCNA_N_stage_association.rds"
)

write.table(
  N_sig,
  file = "result_dir/Workflow-04-2I_WGCNA_N_stage_significant_modules.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# Step 8B-3D
# Module eigengene association with pathological M stage
#
# M0 vs M1
# MX / NA excluded
# Two-sided Wilcoxon rank-sum test
# BH correction
############################################################

# ----------------------------------------------------------
# 1. Check M-stage distribution
# ----------------------------------------------------------

table(
  clinical_WGCNA_traits$M_trait,
  useNA = "ifany"
)

# ----------------------------------------------------------
# 2. Initialize results table
# ----------------------------------------------------------

M_results <- data.frame(
  Module = colnames(MEs_biological),
  P_value = NA_real_,
  N = NA_integer_,
  N_M0 = NA_integer_,
  N_M1 = NA_integer_,
  stringsAsFactors = FALSE
)


# ----------------------------------------------------------
# 3. Two-sided Wilcoxon rank-sum test
# ----------------------------------------------------------

for (i in seq_len(ncol(MEs_biological))) {
  
  module_values <- MEs_biological[, i]
  M_values <- clinical_WGCNA_traits$M_trait
  
  # Complete cases only
  complete_idx <- complete.cases(
    module_values,
    M_values
  )
  
  module_complete <- module_values[complete_idx]
  M_complete <- droplevels(
    M_values[complete_idx]
  )
  
  # Record sample size
  M_results$N[i] <- length(module_complete)
  
  # Record group sizes
  M_results$N_M0[i] <- sum(
    M_complete == "M0"
  )
  
  M_results$N_M1[i] <- sum(
    M_complete == "M1"
  )
  
  # Wilcoxon rank-sum test
  if (length(unique(M_complete)) >= 2) {
    
    test <- wilcox.test(
      module_complete ~ M_complete,
      alternative = "two.sided",
      exact = FALSE
    )
    
    M_results$P_value[i] <- test$p.value
  }
}


# ----------------------------------------------------------
# 4. BH correction
# ----------------------------------------------------------

M_results$FDR <- p.adjust(
  M_results$P_value,
  method = "BH"
)


# ----------------------------------------------------------
# 5. Sort by nominal P value
# ----------------------------------------------------------

M_results <- M_results[
  order(M_results$P_value),
]

rownames(M_results) <- NULL


############################################################
# 8B-3D.3QC
############################################################

# Number of modules
nrow(M_results)

# Unique complete-case sample sizes
unique(M_results$N)

# M0 sample size
unique(M_results$N_M0)

# M1 sample size
unique(M_results$N_M1)

##8B-3D.4 Check whether M0/M1 are entirely correct.

M_complete_idx <- complete.cases(
  clinical_WGCNA_traits$M_trait
)

table(
  clinical_WGCNA_traits$M_trait[M_complete_idx],
  useNA = "ifany"
)

sum(M_complete_idx)

#8B-3D.5 Check core modules: blue and magenta
M_results[
  M_results$Module %in% c(
    "ME_blue",
    "ME_magenta"
  ),
]

#8B-3D.6 Extract modules with an FDR < 0.05
M_sig <- M_results[
  M_results$FDR < 0.05,
  ,
  drop = FALSE
]

M_sig

nrow(M_sig)

#8B-3D.7 Check whether the core module has entered the set of significant modules.
M_results[
  M_results$Module %in% c(
    "ME_blue",
    "ME_magenta"
  ),
  c(
    "Module",
    "P_value",
    "FDR",
    "N",
    "N_M0",
    "N_M1"
  )
]


############################################################
# 8B-3D.8 Save complete M-stage association results
############################################################

write.table(
  M_results,
  file = "result_dir/Workflow-04-2J_WGCNA_M_stage_association.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


saveRDS(
  M_results,
  file = "result_dir/Workflow-04-2J_WGCNA_M_stage_association.rds"
)

#8B-3D.9 Save significant modules
write.table(
  M_sig,
  file = "result_dir/Workflow-04-2J_WGCNA_M_stage_significant_modules.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 8B-3D.10 Final QC summary
############################################################

cat(
  "Number of biological modules:",
  nrow(M_results),
  "\n"
)

cat(
  "Complete-case sample size:",
  unique(M_results$N),
  "\n"
)

cat(
  "M0 sample size:",
  unique(M_results$N_M0),
  "\n"
)

cat(
  "M1 sample size:",
  unique(M_results$N_M1),
  "\n"
)

cat(
  "Number of significant modules (FDR < 0.05):",
  nrow(M_sig),
  "\n"
)

cat(
  "Blue module:\n"
)

print(
  M_results[
    M_results$Module == "ME_blue",
    ,
    drop = FALSE
  ]
)

cat(
  "Magenta module:\n"
)

print(
  M_results[
    M_results$Module == "ME_magenta",
    ,
    drop = FALSE
  ]
)


############################################################
# Step 8C
# Module eigengene association with overall survival
#
# Univariable Cox proportional hazards regression
# 49 biological modules
# BH correction
############################################################

# ----------------------------------------------------------
# 1. Check survival variables
# ----------------------------------------------------------

table(
  clinical_WGCNA_traits$OS_event,
  useNA = "ifany"
)

summary(
  clinical_WGCNA_traits$OS_time
)

# ----------------------------------------------------------
# 2. Basic survival data QC
# ----------------------------------------------------------

sum(
  is.na(clinical_WGCNA_traits$OS_time)
)

sum(
  is.na(clinical_WGCNA_traits$OS_event)
)

sum(
  clinical_WGCNA_traits$OS_time <= 0,
  na.rm = TRUE
)

table(
  clinical_WGCNA_traits$OS_event,
  useNA = "ifany"
)

############################################################
# Step 8C — Survival time QC
# Identify OS.time <= 0 samples
############################################################

OS_time_problem_idx <- which(
  !is.na(clinical_WGCNA_traits$OS_time) &
    clinical_WGCNA_traits$OS_time <= 0
)

length(OS_time_problem_idx)

clinical_WGCNA_traits[
  OS_time_problem_idx,
  c(
    "Sample",
    "patient_id",
    "OS_time",
    "OS_event"
  ),
  drop = FALSE
]

clinical_aligned[
  clinical_aligned$Sample %in%
    clinical_WGCNA_traits$Sample[OS_time_problem_idx],
  c(
    "Sample",
    "patient_id",
    "OS.time",
    "OS.status",
    "Age",
    "Gender",
    "pathologic_stage_2",
    "pathologic_t",
    "pathologic_n",
    "pathologic_m"
  ),
  drop = FALSE
]

############################################################
# Step 8C — Cox model technical QC
# Test with ME_blue before running all 49 modules
############################################################

blue_values <- MEs_biological[, "ME_blue"]

complete_idx_blue <- complete.cases(
  blue_values,
  clinical_WGCNA_traits$OS_time,
  clinical_WGCNA_traits$OS_event
)

blue_complete <- blue_values[complete_idx_blue]
OS_time_blue <- clinical_WGCNA_traits$OS_time[complete_idx_blue]
OS_event_blue <- clinical_WGCNA_traits$OS_event[complete_idx_blue]

# QC
length(blue_complete)

sum(OS_event_blue == 1)

sum(OS_time_blue <= 0)

# Test Cox model
blue_cox <- survival::coxph(
  survival::Surv(
    OS_time_blue,
    OS_event_blue
  ) ~ blue_complete
)

summary(blue_cox)

############################################################
# Step 8C
# Module eigengene association with overall survival
#
# Univariable Cox proportional hazards regression
# 49 biological modules
# Module eigengene treated as a continuous variable
# BH correction
############################################################

library(survival)

# ----------------------------------------------------------
# 1. Initialize results table
# ----------------------------------------------------------

Cox_results <- data.frame(
  Module = colnames(MEs_biological),
  N = NA_integer_,
  Events = NA_integer_,
  HR = NA_real_,
  CI_lower = NA_real_,
  CI_upper = NA_real_,
  P_value = NA_real_,
  stringsAsFactors = FALSE
)


# ----------------------------------------------------------
# 2. Run univariable Cox regression for each module
# ----------------------------------------------------------

for (i in seq_len(ncol(MEs_biological))) {
  
  module_name <- colnames(MEs_biological)[i]
  
  module_values <- MEs_biological[, i]
  
  OS_time <- clinical_WGCNA_traits$OS_time
  OS_event <- clinical_WGCNA_traits$OS_event
  
  # Complete cases
  complete_idx <- complete.cases(
    module_values,
    OS_time,
    OS_event
  )
  
  module_complete <- module_values[complete_idx]
  OS_time_complete <- OS_time[complete_idx]
  OS_event_complete <- OS_event[complete_idx]
  
  # Sample size
  Cox_results$N[i] <- length(
    module_complete
  )
  
  # Number of events
  Cox_results$Events[i] <- sum(
    OS_event_complete == 1
  )
  
  # Cox proportional hazards model
  if (
    length(module_complete) > 0 &&
    length(unique(OS_event_complete)) >= 2
  ) {
    
    cox_model <- survival::coxph(
      survival::Surv(
        OS_time_complete,
        OS_event_complete
      ) ~ module_complete
    )
    
    cox_summary <- summary(
      cox_model
    )
    
    # Hazard ratio
    Cox_results$HR[i] <- exp(
      coef(cox_model)[1]
    )
    
    # 95% confidence interval
    Cox_results$CI_lower[i] <- exp(
      confint(cox_model)[1, 1]
    )
    
    Cox_results$CI_upper[i] <- exp(
      confint(cox_model)[1, 2]
    )
    
    # P value
    Cox_results$P_value[i] <-
      cox_summary$coefficients[
        1,
        "Pr(>|z|)"
      ]
  }
}


# ----------------------------------------------------------
# 3. BH correction
# ----------------------------------------------------------

Cox_results$FDR <- p.adjust(
  Cox_results$P_value,
  method = "BH"
)


# ----------------------------------------------------------
# 4. Sort by nominal P value
# ----------------------------------------------------------

Cox_results <- Cox_results[
  order(Cox_results$P_value),
]

rownames(Cox_results) <- NULL


############################################################
##8C.2 QC
############################################################

# Number of biological modules
nrow(Cox_results)

# Sample size
unique(Cox_results$N)

# Number of events
unique(Cox_results$Events)

# Number of valid Cox models
sum(
  !is.na(Cox_results$P_value)
)

# Number of FDR-significant modules
sum(
  Cox_results$FDR < 0.05,
  na.rm = TRUE
)

#8C.3 View all 49 Cox results
Cox_results

#8C.4 Check the two core modules individually
Cox_results[
  Cox_results$Module %in% c(
    "ME_blue",
    "ME_magenta"
  ),
  c(
    "Module",
    "N",
    "Events",
    "HR",
    "CI_lower",
    "CI_upper",
    "P_value",
    "FDR"
  ),
]

#8C.5 Extract modules with an FDR < 0.05
Cox_sig <- Cox_results[
  Cox_results$FDR < 0.05,
  ,
  drop = FALSE
]

Cox_sig

nrow(Cox_sig)


############################################################
#8C.6 Save complete Cox results
############################################################

write.table(
  Cox_results,
  file = "result_dir/Workflow-04-2K_WGCNA_OS_Cox_association.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  Cox_results,
  file = "result_dir/Workflow-04-2K_WGCNA_OS_Cox_association.rds"
)

#8C.7 Save significant modules
write.table(
  Cox_sig,
  file = "result_dir/Workflow-04-2K_WGCNA_OS_Cox_significant_modules.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 8C.8 Final QC summary
############################################################

cat(
  "Number of biological modules:",
  nrow(Cox_results),
  "\n"
)

cat(
  "Unique sample size(s):",
  paste(
    unique(Cox_results$N),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Unique event number(s):",
  paste(
    unique(Cox_results$Events),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Number of valid Cox models:",
  sum(!is.na(Cox_results$P_value)),
  "\n"
)

cat(
  "Number of significant modules (FDR < 0.05):",
  nrow(Cox_sig),
  "\n"
)

cat("\nBlue module:\n")

print(
  Cox_results[
    Cox_results$Module == "ME_blue",
    ,
    drop = FALSE
  ]
)

cat("\nMagenta module:\n")

print(
  Cox_results[
    Cox_results$Module == "ME_magenta",
    ,
    drop = FALSE
  ]
)

##Step 8C-1 — Proportional Hazards Assumption
############################################################
# Step 8C-1A
# Proportional hazards assumption check
# Five FDR-significant WGCNA modules
############################################################

library(survival)

#-----------------------------------------------------------
# 1. Define FDR-significant modules from Step 8C
#-----------------------------------------------------------

ph_modules <- c(
  "ME_lightgreen",
  "ME_grey60",
  "ME_brown",
  "ME_darkgreen",
  "ME_darkgrey"
)

#-----------------------------------------------------------
# 2. Prepare survival data
#-----------------------------------------------------------

surv_df <- data.frame(
  Sample   = rownames(MEs_biological),
  OS_time  = clinical_WGCNA_traits$OS_time,
  OS_event = clinical_WGCNA_traits$OS_event,
  MEs_biological[, ph_modules, drop = FALSE],
  check.names = FALSE
)

#-----------------------------------------------------------
# 3. Complete-case analysis
#    Same survival cohort as Step 8C
#-----------------------------------------------------------

surv_df_complete <- surv_df[
  complete.cases(surv_df[, c(
    "OS_time",
    "OS_event",
    ph_modules
  )]),
  ,
  drop = FALSE
]

#-----------------------------------------------------------
# 4. QC
#-----------------------------------------------------------

cat("Number of complete survival cases:",
    nrow(surv_df_complete), "\n")

cat("Number of events:",
    sum(surv_df_complete$OS_event == 1), "\n")

cat("Number of censored cases:",
    sum(surv_df_complete$OS_event == 0), "\n")

cat("Number of OS.time = 0 samples:",
    sum(surv_df_complete$OS_time == 0), "\n")

cat("Number of modules tested:",
    length(ph_modules), "\n")


############################################################
# Step 8C-1B
# Cox PH assumption test using cox.zph()
############################################################

ph_results <- list()

for (module in ph_modules) {
  
  #---------------------------------------------------------
  # Univariable Cox model
  #---------------------------------------------------------
  
  formula_module <- as.formula(
    paste0(
      "Surv(OS_time, OS_event) ~ ",
      module
    )
  )
  
  cox_model <- coxph(
    formula_module,
    data = surv_df_complete
  )
  
  #---------------------------------------------------------
  # Proportional hazards assumption
  #---------------------------------------------------------
  
  zph_test <- cox.zph(cox_model)
  
  #---------------------------------------------------------
  # Extract PH test results
  #---------------------------------------------------------
  
  zph_table <- as.data.frame(
    zph_test$table
  )
  
  zph_table$Term <- rownames(zph_table)
  rownames(zph_table) <- NULL
  
  #---------------------------------------------------------
  # Store results
  #---------------------------------------------------------
  
  ph_results[[module]] <- list(
    module = module,
    cox_model = cox_model,
    zph = zph_test,
    zph_table = zph_table
  )
}


############################################################
# Step 8C-1C
# Extract PH assumption results
############################################################

ph_summary <- do.call(
  rbind,
  lapply(ph_modules, function(module) {
    
    zph_table <- ph_results[[module]]$zph_table
    
    module_row <- zph_table[
      zph_table$Term == module,
      ,
      drop = FALSE
    ]
    
    global_row <- zph_table[
      zph_table$Term == "GLOBAL",
      ,
      drop = FALSE
    ]
    
    data.frame(
      Module = module,
      PH_chisq = module_row$chisq,
      PH_df = module_row$df,
      PH_P = module_row$p,
      Global_chisq = global_row$chisq,
      Global_df = global_row$df,
      Global_P = global_row$p,
      stringsAsFactors = FALSE
    )
  })
)

print(ph_summary)


############################################################
# Step 8C-1D
# Interpretation of PH assumption
############################################################

ph_summary$PH_interpretation <- ifelse(
  ph_summary$Global_P >= 0.05,
  "No evidence of PH violation",
  "Potential PH violation"
)

print(ph_summary)


############################################################
# Step 8C-1E
# Save PH assumption results
############################################################

write.table(
  ph_summary,
  file = "result_dir/Workflow-04-2L_WGCNA_OS_PH_assumption.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  ph_results,
  file = "result_dir/Workflow-04-2L_WGCNA_OS_PH_assumption.rds"
)

cat(
  "PH assumption analysis completed.\n",
  "Results saved as:\n",
  "  Workflow-04-2L_WGCNA_OS_PH_assumption.txt\n",
  "  Workflow-04-2L_WGCNA_OS_PH_assumption.rds\n"
)


############################################################
# Step 8C-1F
# Visual inspection of PH assumption
############################################################

plot(
  ph_results[["ME_lightgreen"]]$zph,
  var = 1
)

abline(
  h = 0,
  lty = 2
)

#Then check the other 4 individually.：
plot(
  ph_results[["ME_grey60"]]$zph,
  var = 1
)

plot(
  ph_results[["ME_brown"]]$zph,
  var = 1
)

plot(
  ph_results[["ME_darkgreen"]]$zph,
  var = 1
)

plot(
  ph_results[["ME_darkgrey"]]$zph,
  var = 1
)


############################################################
#8C-1G. Final QC
############################################################

cat("\n========== Step 8C-1 Final QC ==========\n")

cat("Modules tested:", nrow(ph_summary), "\n")
cat("Complete cases:", nrow(surv_df_complete), "\n")
cat("Events:",
    sum(surv_df_complete$OS_event == 1), "\n")
cat("OS.time = 0:",
    sum(surv_df_complete$OS_time == 0), "\n")

cat("\nPH assumption results:\n")
print(
  ph_summary[, c(
    "Module",
    "PH_chisq",
    "PH_df",
    "PH_P",
    "PH_interpretation"
  )],
  row.names = FALSE
)

cat("\n========================================\n")

############################################################
# Step 8C-1G
# Final PH assumption classification
############################################################

ph_summary$PH_class <- ifelse(
  ph_summary$PH_P < 0.05,
  "Potential PH violation",
  "No evidence of PH violation"
)

ph_summary$Trend <- c(
  "Upward",
  "Approximately horizontal",
  "Downward",
  "Upward",
  "Upward"
)

ph_summary

write.table(
  ph_summary,
  file = "result_dir/Workflow-04-2L_WGCNA_OS_PH_assumption_final.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# Step 8C-2A
# Clinical covariate design for multivariable Cox analysis
# No final multivariable Cox model is fitted at this step
############################################################

library(survival)

#-----------------------------------------------------------
# 1. Clinical covariates
#-----------------------------------------------------------

clinical_cox <- clinical_WGCNA_traits[, c(
  "Sample",
  "Age_trait",
  "Gender_trait",
  "Stage_trait",
  "T_trait",
  "N_trait",
  "M_trait",
  "OS_time",
  "OS_event"
)]

#-----------------------------------------------------------
# 2. Convert categorical clinical variables to factors
#-----------------------------------------------------------

clinical_cox$Gender_trait <- factor(
  clinical_cox$Gender_trait,
  levels = c(0, 1),
  labels = c("Female", "Male")
)

clinical_cox$Stage_trait <- factor(
  clinical_cox$Stage_trait,
  levels = c("Stage I", "Stage II", "Stage III", "Stage IV"),
  ordered = FALSE
)

clinical_cox$T_trait <- factor(
  clinical_cox$T_trait,
  levels = c("T1", "T2", "T3", "T4"),
  ordered = FALSE
)

clinical_cox$N_trait <- factor(
  clinical_cox$N_trait,
  levels = c("N0", "N1", "N2", "N3"),
  ordered = FALSE
)

clinical_cox$M_trait <- factor(
  clinical_cox$M_trait,
  levels = c("M0", "M1"),
  ordered = FALSE
)

#-----------------------------------------------------------
# 3. Merge with five FDR-significant modules
#-----------------------------------------------------------

ph_modules <- c(
  "ME_lightgreen",
  "ME_grey60",
  "ME_brown",
  "ME_darkgreen",
  "ME_darkgrey"
)

cox_design <- merge(
  clinical_cox,
  data.frame(
    Sample = rownames(MEs_biological),
    MEs_biological[, ph_modules, drop = FALSE],
    check.names = FALSE
  ),
  by = "Sample",
  sort = FALSE
)

#-----------------------------------------------------------
# 4. Restrict to complete survival cases
#-----------------------------------------------------------

cox_design_complete <- cox_design[
  complete.cases(
    cox_design[, c(
      "OS_time",
      "OS_event",
      ph_modules
    )]
  ),
  ,
  drop = FALSE
]

#-----------------------------------------------------------
# 5. QC
#-----------------------------------------------------------

cat("\n========== Step 8C-2A QC ==========\n")

cat("Samples:", nrow(cox_design_complete), "\n")
cat("Events:",
    sum(cox_design_complete$OS_event == 1),
    "\n")

cat("OS.time = 0:",
    sum(cox_design_complete$OS_time == 0),
    "\n")

cat("\nClinical variable missingness:\n")

print(
  sapply(
    cox_design_complete[, c(
      "Age_trait",
      "Gender_trait",
      "Stage_trait",
      "T_trait",
      "N_trait",
      "M_trait"
    )],
    function(x) sum(is.na(x))
  )
)

cat("\n====================================\n")


############################################################
# Step 8C-2B
# Clinical covariate overlap
############################################################

cat("\n========== Clinical Covariate Overlap ==========\n")

cat("\nStage vs T:\n")
print(
  table(
    cox_design_complete$Stage_trait,
    cox_design_complete$T_trait,
    useNA = "ifany"
  )
)

cat("\nStage vs N:\n")
print(
  table(
    cox_design_complete$Stage_trait,
    cox_design_complete$N_trait,
    useNA = "ifany"
  )
)

cat("\nStage vs M:\n")
print(
  table(
    cox_design_complete$Stage_trait,
    cox_design_complete$M_trait,
    useNA = "ifany"
  )
)

cat("\nT vs N:\n")
print(
  table(
    cox_design_complete$T_trait,
    cox_design_complete$N_trait,
    useNA = "ifany"
  )
)

cat("\nT vs M:\n")
print(
  table(
    cox_design_complete$T_trait,
    cox_design_complete$M_trait,
    useNA = "ifany"
  )
)

cat("\nN vs M:\n")
print(
  table(
    cox_design_complete$N_trait,
    cox_design_complete$M_trait,
    useNA = "ifany"
  )
)

cat("\n===============================================\n")


############################################################
# Step 8C-2C
# Age and Gender QC
############################################################

cat("\n========== Age and Gender QC ==========\n")

cat("\nAge summary:\n")
print(
  summary(cox_design_complete$Age_trait)
)

cat("\nGender distribution:\n")
print(
  table(cox_design_complete$Gender_trait)
)

cat("\nAge by Gender:\n")
print(
  tapply(
    cox_design_complete$Age_trait,
    cox_design_complete$Gender_trait,
    summary
  )
)

cat("\n========================================\n")


##Step 8C-2D — Identify two candidate clinical adjustment frameworks
#Model A：overall pathological stage model
#Model B：anatomic component model

#Step 8C-2E — Modules in different pH states should be handled differently.
#Conventional Cox candidate
conventional_modules <- c(
  "ME_grey60",
  "ME_brown"
)

#Time-dependent / sensitivity candidate
time_dependent_modules <- c(
  "ME_lightgreen",
  "ME_darkgreen",
  "ME_darkgrey"
)

############################################################
# Step 8C-2E
# Module classification according to PH assumption
############################################################

conventional_modules <- c(
  "ME_grey60",
  "ME_brown"
)

time_dependent_modules <- c(
  "ME_lightgreen",
  "ME_darkgreen",
  "ME_darkgrey"
)

module_design <- data.frame(
  Module = c(
    conventional_modules,
    time_dependent_modules
  ),
  PH_status = c(
    rep("No evidence of PH violation",
        length(conventional_modules)),
    rep("Potential PH violation",
        length(time_dependent_modules))
  ),
  Recommended_model = c(
    rep("Conventional Cox",
        length(conventional_modules)),
    rep("Time-dependent / sensitivity analysis",
        length(time_dependent_modules))
  ),
  stringsAsFactors = FALSE
)

print(module_design)


############################################################
# Step 8C-2F
# Save model design only
############################################################

model_design <- list(
  
  survival_cohort_n = nrow(cox_design_complete),
  
  events = sum(
    cox_design_complete$OS_event == 1
  ),
  
  zero_time_samples = sum(
    cox_design_complete$OS_time == 0
  ),
  
  primary_clinical_model = "Age + Gender + Stage",
  
  sensitivity_clinical_model = "Age + Gender + T + N + M",
  
  conventional_modules = conventional_modules,
  
  time_dependent_modules = time_dependent_modules,
  
  PH_results = ph_summary,
  
  note = paste(
    "Final multivariable Cox models were not fitted at this step."
  )
)

saveRDS(
  model_design,
  file = "result_dir/Workflow-04-2M_WGCNA_OS_multivariable_model_design.rds"
)

write.table(
  module_design,
  file = "result_dir/Workflow-04-2M_WGCNA_OS_multivariable_model_design.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nStep 8C-2 model design saved.\n"
)

############################################################
# Step 8C-2G
# Define complete-case cohorts for each clinical model
# No Cox model is fitted at this step
############################################################

#-----------------------------------------------------------
# 1. Define clinical variables
#-----------------------------------------------------------

primary_covariates <- c(
  "Age_trait",
  "Gender_trait",
  "Stage_trait"
)

sensitivity_covariates <- c(
  "Age_trait",
  "Gender_trait",
  "T_trait",
  "N_trait",
  "M_trait"
)

#-----------------------------------------------------------
# 2. Function to define complete-case cohort
#-----------------------------------------------------------

get_complete_cohort <- function(
    data,
    module,
    covariates
) {
  
  required_vars <- c(
    "OS_time",
    "OS_event",
    module,
    covariates
  )
  
  keep <- complete.cases(
    data[, required_vars, drop = FALSE]
  )
  
  data[
    keep,
    ,
    drop = FALSE
  ]
}

#-----------------------------------------------------------
# 3. Primary model cohorts
#-----------------------------------------------------------

primary_cohorts <- lapply(
  c(
    conventional_modules,
    time_dependent_modules
  ),
  function(module) {
    
    get_complete_cohort(
      data = cox_design,
      module = module,
      covariates = primary_covariates
    )
  }
)

names(primary_cohorts) <- c(
  conventional_modules,
  time_dependent_modules
)

#-----------------------------------------------------------
# 4. Sensitivity model cohorts
#-----------------------------------------------------------

sensitivity_cohorts <- lapply(
  c(
    conventional_modules,
    time_dependent_modules
  ),
  function(module) {
    
    get_complete_cohort(
      data = cox_design,
      module = module,
      covariates = sensitivity_covariates
    )
  }
)

names(sensitivity_cohorts) <- c(
  conventional_modules,
  time_dependent_modules
)


############################################################
# Step 8C-2G QC
############################################################

module_order <- c(
  conventional_modules,
  time_dependent_modules
)

cohort_summary <- do.call(
  rbind,
  lapply(
    module_order,
    function(module) {
      
      primary_data <- primary_cohorts[[module]]
      sensitivity_data <- sensitivity_cohorts[[module]]
      
      data.frame(
        Module = module,
        
        Primary_N = nrow(primary_data),
        Primary_Events = sum(
          primary_data$OS_event == 1
        ),
        Primary_ZeroTime = sum(
          primary_data$OS_time == 0
        ),
        
        Sensitivity_N = nrow(sensitivity_data),
        Sensitivity_Events = sum(
          sensitivity_data$OS_event == 1
        ),
        Sensitivity_ZeroTime = sum(
          sensitivity_data$OS_time == 0
        ),
        
        Sensitivity_M0 = sum(
          sensitivity_data$M_trait == "M0"
        ),
        Sensitivity_M1 = sum(
          sensitivity_data$M_trait == "M1"
        ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)

print(
  cohort_summary,
  row.names = FALSE
)

##############################################################
# Step 8C-2H
# Final QC of complete-case cohorts before multivariable Cox
############################################################

cat("\n==============================\n")
cat("Step 8C-2H: Final cohort QC\n")
cat("==============================\n\n")

module_order <- c(
  conventional_modules,
  time_dependent_modules
)

############################################################
# 1. Primary cohort QC
############################################################

primary_qc <- do.call(
  rbind,
  lapply(
    module_order,
    function(module) {
      
      dat <- primary_cohorts[[module]]
      
      data.frame(
        Module = module,
        N = nrow(dat),
        Events = sum(dat$OS_event == 1),
        Censored = sum(dat$OS_event == 0),
        ZeroTime = sum(dat$OS_time == 0),
        Age_N = sum(!is.na(dat$Age_trait)),
        Gender_N = sum(!is.na(dat$Gender_trait)),
        Stage_N = sum(!is.na(dat$Stage_trait)),
        stringsAsFactors = FALSE
      )
    }
  )
)

print(
  primary_qc,
  row.names = FALSE
)


############################################################
# 2. Primary cohort: Stage distribution
############################################################

cat("\n--- Primary cohort: Stage distribution ---\n")

table(
  primary_cohorts[[module_order[1]]]$Stage_trait,
  useNA = "ifany"
)


############################################################
# 3. Primary cohort: Gender distribution
############################################################

cat("\n--- Primary cohort: Gender distribution ---\n")

table(
  primary_cohorts[[module_order[1]]]$Gender_trait,
  useNA = "ifany"
)


############################################################
# 4. Sensitivity cohort QC
############################################################

sensitivity_qc <- do.call(
  rbind,
  lapply(
    module_order,
    function(module) {
      
      dat <- sensitivity_cohorts[[module]]
      
      data.frame(
        Module = module,
        N = nrow(dat),
        Events = sum(dat$OS_event == 1),
        Censored = sum(dat$OS_event == 0),
        ZeroTime = sum(dat$OS_time == 0),
        Age_N = sum(!is.na(dat$Age_trait)),
        Gender_N = sum(!is.na(dat$Gender_trait)),
        T_N = sum(!is.na(dat$T_trait)),
        N_N = sum(!is.na(dat$N_trait)),
        M_N = sum(!is.na(dat$M_trait)),
        stringsAsFactors = FALSE
      )
    }
  )
)

cat("\n--- Sensitivity cohort QC ---\n")

print(
  sensitivity_qc,
  row.names = FALSE
)


############################################################
# 5. Sensitivity cohort: Stage distribution
############################################################

cat("\n--- Sensitivity cohort: Stage distribution ---\n")

table(
  sensitivity_cohorts[[module_order[1]]]$Stage_trait,
  useNA = "ifany"
)


############################################################
# 6. Sensitivity cohort: T distribution
############################################################

cat("\n--- Sensitivity cohort: T distribution ---\n")

table(
  sensitivity_cohorts[[module_order[1]]]$T_trait,
  useNA = "ifany"
)


############################################################
# 7. Sensitivity cohort: N distribution
############################################################

cat("\n--- Sensitivity cohort: N distribution ---\n")

table(
  sensitivity_cohorts[[module_order[1]]]$N_trait,
  useNA = "ifany"
)


############################################################
# 8. Sensitivity cohort: M distribution
############################################################

cat("\n--- Sensitivity cohort: M distribution ---\n")

table(
  sensitivity_cohorts[[module_order[1]]]$M_trait,
  useNA = "ifany"
)


############################################################
# 9. Factor levels
############################################################

cat("\n--- Factor levels ---\n")

cat("Gender levels:\n")
print(levels(cox_design$Gender_trait))

cat("\nStage levels:\n")
print(levels(cox_design$Stage_trait))

cat("\nT levels:\n")
print(levels(cox_design$T_trait))

cat("\nN levels:\n")
print(levels(cox_design$N_trait))

cat("\nM levels:\n")
print(levels(cox_design$M_trait))


############################################################
# 10. Save final cohort QC
############################################################

saveRDS(
  list(
    primary_qc = primary_qc,
    sensitivity_qc = sensitivity_qc,
    primary_cohorts = primary_cohorts,
    sensitivity_cohorts = sensitivity_cohorts
  ),
  file = "result_dir/Workflow-04-2N_WGCNA_OS_multivariable_cohort_QC.rds"
)

write.table(
  primary_qc,
  file = "result_dir/Workflow-04-2N_WGCNA_OS_primary_cohort_QC.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sensitivity_qc,
  file = "result_dir/Workflow-04-2N_WGCNA_OS_sensitivity_cohort_QC.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nStep 8C-2H completed.\n")


###Step 8C-2I：Primary conventional multivariable Cox

#Model：
#ME_grey60 + Age + Gender + Stage
#and
#ME_brown + Age + Gender + Stage
############################################################
# Step 8C-2I
# Primary multivariable Cox models
# Conventional PH models
#
# Modules:
#   ME_grey60
#   ME_brown
#
# Covariates:
#   Age + Gender + Stage
############################################################

library(survival)

cat("\n========================================\n")
cat("Step 8C-2I: Primary multivariable Cox\n")
cat("========================================\n\n")


############################################################
# 1. Define conventional modules
############################################################

conventional_modules <- c(
  "ME_grey60",
  "ME_brown"
)


############################################################
# 2. Primary model function
############################################################

fit_primary_cox <- function(module, data) {
  
  formula <- as.formula(
    paste(
      "Surv(OS_time, OS_event) ~",
      module,
      "+ Age_trait",
      "+ Gender_trait",
      "+ Stage_trait"
    )
  )
  
  model <- coxph(
    formula,
    data = data,
    ties = "efron"
  )
  
  return(model)
}


############################################################
# 3. Fit models
############################################################

primary_models <- lapply(
  conventional_modules,
  function(module) {
    
    dat <- primary_cohorts[[module]]
    
    fit_primary_cox(
      module = module,
      data = dat
    )
  }
)

names(primary_models) <- conventional_modules


############################################################
# 4. Print model summaries
############################################################

for (module in conventional_modules) {
  
  cat("\n----------------------------------------\n")
  cat("Module:", module, "\n")
  cat("----------------------------------------\n")
  
  print(
    summary(primary_models[[module]])
  )
}


############################################################
# 5. Extract HR / 95% CI / P
############################################################

extract_cox_results <- function(model, module) {
  
  s <- summary(model)
  
  coef_table <- as.data.frame(
    s$coefficients
  )
  
  conf_table <- as.data.frame(
    s$conf.int
  )
  
  result <- data.frame(
    Module = module,
    Variable = rownames(coef_table),
    Beta = coef_table[, "coef"],
    HR = conf_table[, "exp(coef)"],
    Lower95CI = conf_table[, "lower .95"],
    Upper95CI = conf_table[, "upper .95"],
    P = coef_table[, "Pr(>|z|)"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  return(result)
}


primary_results <- do.call(
  rbind,
  lapply(
    conventional_modules,
    function(module) {
      
      extract_cox_results(
        model = primary_models[[module]],
        module = module
      )
    }
  )
)


############################################################
# 6. Print extracted results
############################################################

cat("\n========================================\n")
cat("Primary multivariable Cox results\n")
cat("========================================\n\n")

print(
  primary_results,
  row.names = FALSE
)


############################################################
# 7. Model-level QC
############################################################

primary_model_qc <- do.call(
  rbind,
  lapply(
    conventional_modules,
    function(module) {
      
      model <- primary_models[[module]]
      s <- summary(model)
      
      data.frame(
        Module = module,
        N = model$n,
        Events = model$nevent,
        Concordance = unname(
          s$concordance[1]
        ),
        Concordance_SE = unname(
          s$concordance[2]
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

print(
  primary_model_qc,
  row.names = FALSE
)



############################################################
# 8. Save objects
############################################################

saveRDS(
  list(
    models = primary_models,
    results = primary_results,
    model_qc = primary_model_qc
  ),
  file = "result_dir/Workflow-04-2O_WGCNA_OS_primary_conventional_Cox.rds"
)

write.table(
  primary_results,
  file = "result_dir/Workflow-04-2O_WGCNA_OS_primary_conventional_Cox.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  primary_model_qc,
  file = "result_dir/Workflow-04-2O_WGCNA_OS_primary_conventional_Cox_QC.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nStep 8C-2I completed.\n")



############################################################
# Step 8C-2J
# PH assumption for primary multivariable Cox models
############################################################

library(survival)

cat("\n========================================\n")
cat("Step 8C-2J: Multivariable PH assumption\n")
cat("========================================\n\n")


############################################################
# 1. cox.zph for ME_grey60 model
############################################################

zph_grey60 <- cox.zph(
  primary_models[["ME_grey60"]],
  transform = "km"
)

cat("\n--- ME_grey60 multivariable model ---\n")
print(zph_grey60)


############################################################
# 2. cox.zph for ME_brown model
############################################################

zph_brown <- cox.zph(
  primary_models[["ME_brown"]],
  transform = "km"
)

cat("\n--- ME_brown multivariable model ---\n")
print(zph_brown)


############################################################
# 3. Extract PH test results
############################################################

extract_ph_results <- function(zph_object, module) {
  
  tab <- as.data.frame(
    zph_object$table
  )
  
  tab$Variable <- rownames(tab)
  
  result <- data.frame(
    Module = module,
    Variable = tab$Variable,
    ChiSquare = tab[, "chisq"],
    P = tab[, "p"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  return(result)
}


ph_results <- rbind(
  extract_ph_results(
    zph_grey60,
    "ME_grey60"
  ),
  extract_ph_results(
    zph_brown,
    "ME_brown"
  )
)

cat("\n========================================\n")
cat("Multivariable PH test results\n")
cat("========================================\n\n")

print(
  ph_results,
  row.names = FALSE
)


############################################################
# 4. Save PH results
############################################################

saveRDS(
  list(
    ME_grey60 = zph_grey60,
    ME_brown = zph_brown,
    results = ph_results
  ),
  file = "result_dir/Workflow-04-2P_WGCNA_OS_multivariable_PH.rds"
)

write.table(
  ph_results,
  file = "result_dir/Workflow-04-2P_WGCNA_OS_multivariable_PH.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nStep 8C-2J completed.\n")


##############################################################
# Step 8C-2K
# Schoenfeld residual diagnostics
############################################################

par(
  mfrow = c(2, 2)
)

plot(
  zph_grey60,
  var = "ME_grey60",
  main = "ME_grey60"
)

plot(
  zph_grey60,
  var = "Gender_trait",
  main = "Gender — ME_grey60 model"
)

plot(
  zph_brown,
  var = "ME_brown",
  main = "ME_brown"
)

plot(
  zph_brown,
  var = "Stage_trait",
  main = "Stage — ME_brown model"
)

par(
  mfrow = c(1, 1)
)

############################################################
# Detailed PH diagnostic summaries
############################################################

cat("\n========================================\n")
cat("ME_grey60 PH diagnostic\n")
cat("========================================\n")

print(zph_grey60)

cat("\n========================================\n")
cat("ME_brown PH diagnostic\n")
cat("========================================\n")

print(zph_brown)

cat("\n========================================\n")
cat("ME_brown residual structure\n")
cat("========================================\n")

print(
  zph_brown$y[, "ME_brown"]
  |> summary()
)

# ============================================================
# Paper 1 LUAD WGCNA
# End-of-day checkpoint
# Date: 
# ============================================================

checkpoint_dir <- result_dir

if (!dir.exists(checkpoint_dir)) {
  dir.create(checkpoint_dir, recursive = TRUE)
}

# ------------------------------------------------------------
# 1. Save complete R workspace
# ------------------------------------------------------------

save.image(
  file = file.path(
    checkpoint_dir,
    "Workspace_20260911_WGCNA_Step8C.RData"
  )
)


# ------------------------------------------------------------
# 2. Save key Paper 1 WGCNA analysis objects
# ------------------------------------------------------------

save(
  clinical_aligned,
  datExpr_MAD,
  gene_MAD,
  net,
  moduleColors,
  moduleSizes,
  MEs,
  MEs_biological,
  clinical_WGCNA_traits,
  clinical_cox,
  cox_design_complete,
  primary_model_qc,
  primary_models,
  primary_results,
  zph_grey60,
  zph_brown,
  file = file.path(
    checkpoint_dir,
    "Workflow-04-2N_WGCNA_PH_multivariable_checkpoint_20260911.RData"
  )
)

# ------------------------------------------------------------
# 3. Save R session information
# ------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    checkpoint_dir,
    "SessionInfo_20260911.txt"
  )
)

# ------------------------------------------------------------
# 4. Save the list of current workspace objects
# ------------------------------------------------------------

writeLines(
  capture.output(ls()),
  con = file.path(
    checkpoint_dir,
    "Workspace_objects_20260911.txt"
  )
)

# ------------------------------------------------------------
# 5. Verify saved files
# ------------------------------------------------------------

list.files(
  checkpoint_dir,
  full.names = TRUE
)

file.exists(
  "result_dir/Workspace_20260911_WGCNA_Step8C.RData"
)

file.exists(
  "result_dir/Workflow-04-2N_WGCNA_PH_multivariable_checkpoint_20260911.RData"
)

file.exists(
  "result_dir/SessionInfo_20260911.txt"
)





